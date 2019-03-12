
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c revision log :
c
c 09.01.1994	ggu	(from scratch)
c 19.01.1994	ggu	$$flux - flux conserving property
c 20.01.1994	ggu	$$iclin - iclin not used to compute volume
c 20.01.1994	ggu	$$lumpc - evaluate conz. nodewise
c 03.02.1994	ggu	$$itot0 - exception for itot=0 or 3
c 04.02.1994	ggu	$$fact3 - factor 3 missing in transport
c 04.02.1994	ggu	$$azpar - azpar used to compute transport
c 04.02.1994	ggu	$$condry - comute conz also in dry areas
c 07.02.1994	ggu	$$istot - istot for fractional time step
c 01.06.1994	ggu	restructured for 3-d model
c 18.07.1994	ggu	$$htop - use htop instead of htopo for mass cons.
c 09.04.1996	ggu	$$rvadj adjust rv in certain areas
c 20.05.2015    erp     transformed for OMP
c 30.09.2015    ggu     routine cleaned, no reals in conz3d
c 20.11.2015    ggu&erp chunk size introduced, omp finalized
c 20.10.2016    ccf     pass rtauv for differential nudging
c 11.05.2018    ggu     compute only unique nodes (needed for zeta layers)
c 11.10.2018    ggu     code adjusted for sediment deposition (negative loads)
c 01.02.2019    ggu     bug fix for conz==0 with negative loading
c 14.02.2019    ggu     bug fix for conz<0 with negative loading
c
c**************************************************************

        subroutine conz3d_omp(cn1,co1
     +			,ddt
     +                  ,rkpar,difhv,difv
     +			,difmol,cbound
     +		 	,itvd,itvdv,gradxv,gradyv
     +			,cobs,robs,rtauv
     +			,wsink,wsinkv
     +			,rload,load
     +			,azpar,adpar,aapar
     +			,istot,isact,nlvddi
     +                  ,nlev)
     
c computes concentration
c
c cn     new concentration
c co     old concentration              !not used !FIXME
c caux   aux vector
c clow	 lower diagonal of vertical system
c chig	 upper diagonal of vertical system
c ddt    time step
c rkpar  horizontal turbulent diffusivity
c difhv  horizontal turbulent diffusivity (variable between elements)
c difv   vertical turbulent diffusivity
c difmol vertical molecular diffusivity
c cbound boundary condition (mass flux) [kg/s] -> now concentration [kg/m**3]
c itvd	 type of horizontal transport algorithm used
c itvdv	 type of vertical transport algorithm used
c gradxv,gradyv  gradient vectors for TVD algorithm
c cobs	 observations for nudging
c robs	 use observations for nuding (real)
c rtauv	 variable relaxation coefficient (real)
c wsink	 factor for settling velocity
c wsinkv variable settling velocity [m/s]
c rload	 factor for loading
c load   load (source or sink) [kg/s]
c azpar  time weighting parameter
c adpar  time weighting parameter for vertical diffusion (ad)
c aapar  time weighting parameter for vertical advection (aa)
c istot	 total inter time steps
c isact	 actual inter time step
c nlvddi	 dimension in z direction
c nlv	 actual needed levels
c
c solution of purely diffusional part :
c
c dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
c
c C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )
c
c for n-dimensions and
c
c C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )
c
c for 1 dimension
c
c the solution is normalized, i.e.  int(C(x,t)dx) = 1 over the whole area
c
c DPGGU -> introduced double precision to stabilize solution

	use mod_bound_geom
	use mod_geom
	use mod_depth
	use mod_diff_aux
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
!$	use omp_lib
	use mod_subset
	use shympi

	implicit none

	integer, intent(in) :: nlvddi,nlev,itvd,itvdv,istot,isact
	real, intent(in) :: difmol,robs,wsink,rload,ddt,rkpar
	real, intent(in) :: azpar,adpar,aapar
	real,dimension(nlvddi,nkn),intent(inout) :: cn1
	real,dimension(nlvddi,nkn),intent(in) :: co1,cbound
	real,dimension(nlvddi,nel),intent(in) :: difhv
	real,dimension(nlvddi,nkn),intent(in) :: gradxv,gradyv
	real,dimension(nlvddi,nkn),intent(in) :: cobs,rtauv
	real,dimension(nlvddi,nkn),intent(inout) :: load		!LLL
	real,dimension(0:nlvddi,nkn),intent(in) :: difv,wsinkv
        !double precision,dimension(nlvddi,nkn),intent(out) :: cn
        
	logical :: btvdv
	integer :: ie,k,ilevel,ibase,ii,l,n,i,j,x,ies,iend,kl,kend,ntot
	integer :: myid,numthreads,j_init,j_end,knod,k_end,jel
	integer,allocatable,dimension(:) :: subset_l
	real :: time1,time2
	double precision :: dtime1,dtime2
	integer :: nchunk,nthreads,nelems,nnodes
	double precision :: dt
	double precision :: az,ad,aa,azt,adt,aat,an,ant
	double precision :: rstot,rso,rsn,rsot,rsnt
	double precision :: timer,timer1,chunk,rest
	
! 	double precision,dimension(nlvddi,nkn) :: cn
! 	double precision,dimension(nlvddi,nkn) :: co        
!         double precision,dimension(nlvddi,nkn) :: cdiag
! 	double precision,dimension(nlvddi,nkn) :: clow
! 	double precision,dimension(nlvddi,nkn) :: chigh
	
	double precision,dimension(:,:),allocatable :: cn        
	double precision,dimension(:,:),allocatable :: co        
        double precision,dimension(:,:),allocatable :: cdiag
	double precision,dimension(:,:),allocatable :: clow
	double precision,dimension(:,:),allocatable :: chigh

        if(nlv.ne.nlev) stop 'error stop conz3d_omp: nlv/=nlev'
	
c----------------------------------------------------------------
c initialize variables and parameters
c----------------------------------------------------------------

!	call cpu_time(time1)
!!$	dtime1 = omp_get_wtime()
	
	ALLOCATE(cn(nlvddi,nkn))
	ALLOCATE(co(nlvddi,nkn))
	ALLOCATE(cdiag(nlvddi,nkn))
	ALLOCATE(clow(nlvddi,nkn))
	ALLOCATE(chigh(nlvddi,nkn))
	
	az = azpar
	ad = adpar
	aa = aapar
	an = 0.			!implicit parameter nudging
	
	azt=1.-az
	adt=1.-ad
	aat=1.-aa
	ant=1.-an

	rstot = istot			!ERIC - what a brown paper bag bug
	rso=(isact-1)/rstot
	rsn=(isact)/rstot
	rsot=1.-rso
	rsnt=1.-rsn

	dt=ddt/rstot
	
	btvdv = itvdv .gt. 0
	if( btvdv .and. aapar .ne. 0. ) then
	  write(6,*) 'aapar = ',aapar,'  itvdv = ',itvdv
	  write(6,*) 'Cannot use implicit vertical advection'
	  write(6,*) 'together with vertical TVD scheme.'
	  write(6,*) 'Please set either aapar = 0 (explicit) or'
	  write(6,*) 'itvdv = 0 (no vertical TVD) in the STR file.'
	  stop 'error stop conz3d: vertical tvd scheme'
	end if

        cn=0.
	co=cn1
        cdiag=0.
        clow=0.
        chigh=0.

        nchunk = 1
	nthreads = 1
!$	nthreads = omp_get_num_threads()
 
      do i=1,subset_num 	! loop over indipendent subset
       
!$     nchunk = subset_el(i) / ( nthreads * 10 )
       nchunk = max(nchunk,1)

!$OMP TASKWAIT 
!!!$OMP TASKGROUP 
       do jel=1,subset_el(i),nchunk

!$OMP TASK FIRSTPRIVATE(jel,i) DEFAULT(NONE)
!$OMP& PRIVATE(j,ie)
!$OMP& SHARED(nlvddi,nlev,itvd,itvdv,istot,isact,aa,nchunk)
!$OMP& SHARED(difmol,robs,wsink,rload,ddt,rkpar,az,ad)
!$OMP& SHARED(an,ant)
!$OMP& SHARED(azt,adt,aat,rso,rsn,rsot,rsnt,dt,nkn)
!$OMP& SHARED(cn,co,cdiag,clow,chigh,subset_el,cn1,co1) 
!$OMP& SHARED(subset_num,indipendent_subset) 
!$OMP& SHARED(difhv,cbound,gradxv,gradyv,cobs,rtauv,load,difv,wsinkv)

       do j=jel,jel+nchunk-1 	! loop over elements in subset
		if(j .le. subset_el(i)) then
	        ie = indipendent_subset(j,i)
	        !print *,i,ie
                call conz3d_element(ie,cdiag,clow,chigh,cn,cn1
     +			,dt
     +                  ,rkpar,difhv,difv
     +			,difmol,cbound
     +		 	,itvd,itvdv,gradxv,gradyv
     +			,cobs,robs,rtauv
     +			,wsink,wsinkv
     +			,rload,load
     +			,az,ad,aa,azt,adt,aat,an,ant
     +			,rso,rsn,rsot,rsnt
     +			,nlvddi,nlev)
		end if
	end do ! end loop over el in subset
!$OMP END TASK
      end do

!!!$OMP END TASKGROUP       
!$OMP TASKWAIT       

       end do ! end loop over subset
       
       if( shympi_partition_on_elements() ) then
         !call shympi_comment('shympi_elem: exchange scalar')
         call shympi_exchange_and_sum_3d_nodes(cn)
         call shympi_exchange_and_sum_3d_nodes(cdiag)
         call shympi_exchange_and_sum_3d_nodes(clow)
         call shympi_exchange_and_sum_3d_nodes(chigh)
       end if

       ntot = nkn
       if( shympi_partition_on_nodes() ) ntot = nkn_unique
!$     nchunk = ntot / ( nthreads * 10 )
       nchunk = max(nchunk,1)

!$OMP TASKWAIT
!!!$OMP TASKGROUP
       do knod=1,ntot,nchunk
!$OMP TASK FIRSTPRIVATE(knod) PRIVATE(k) DEFAULT(NONE)
!$OMP& SHARED(cn,cdiag,clow,chigh,cn1,cbound,load,nchunk,
!$OMP&           rload,ad,aa,dt,nlvddi,ntot)
	 do k=knod,knod+nchunk-1
	 if(k .le. ntot) then
	   call conz3d_nodes(k,cn,cdiag(:,k),clow(:,k),chigh(:,k),
     +                          cn1,cbound,load,rload,
     +                          ad,aa,dt,nlvddi)
         endif
         enddo
!$OMP END TASK 	      
	end do

!!!$OMP END TASKGROUP
!$OMP TASKWAIT       

	!cn1 = 0.
	!cn1 = cn
	cn1 = real(cn)
	
	DEALLOCATE(cn)
	DEALLOCATE(co)
	DEALLOCATE(cdiag)
	DEALLOCATE(clow)
	DEALLOCATE(chigh)
	
!	call cpu_time(time2)
!!$	dtime2 = omp_get_wtime()
!	write(6,*) time2-time1,dtime2-dtime1

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*****************************************************************

       subroutine conz3d_element(ie
     +			,cdiag,clow,chigh,cn,cn1
     +			,dt
     +                  ,rkpar,difhv,difv
     +			,difmol,cbound
     +		 	,itvd,itvdv,gradxv,gradyv
     +			,cobs,robs,rtauv
     +			,wsink,wsinkv
     +			,rload,load
     +			,az,ad,aa,azt,adt,aat,an,ant
     +			,rso,rsn,rsot,rsnt
     +			,nlvddi,nlev)
     
        use mod_bound_geom
	use mod_geom
	use mod_depth
	use mod_diff_aux
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
	use mod_layer_thickness
      
      implicit none
      
      integer,intent(in) :: ie,nlvddi,nlev,itvd,itvdv
      real,intent(in) :: difmol,robs,wsink,rload,rkpar
      real,dimension(nlvddi,nkn),intent(in) :: cn1,cbound
      real,dimension(nlvddi,nel),intent(in) :: difhv
      real,dimension(nlvddi,nkn),intent(in) :: gradxv,gradyv
      real,dimension(nlvddi,nkn),intent(in) :: cobs,rtauv,load
      real,intent(in),dimension(0:nlvddi,nkn) :: wsinkv,difv
      double precision,intent(in) :: dt
      double precision,intent(in) :: az,ad,aa,azt,adt,aat,an,ant
      double precision,intent(in) :: rso,rsn,rsot,rsnt
      double precision,dimension(nlvddi,nkn),intent(inout) :: cdiag
      double precision,dimension(nlvddi,nkn),intent(inout) :: clow
      double precision,dimension(nlvddi,nkn),intent(inout) :: chigh
      double precision,dimension(nlvddi,nkn),intent(inout) :: cn
        
        logical :: btvdv,btvd,bgradup
	integer :: k,ii,l,iii,ll,ibase,lstart,ilevel,itot,isum
	integer :: n,i,iext
	integer, dimension(3) :: kn
        double precision :: cexpl,cbm,ccm,waux,loading,wws,us,vs
        double precision :: aj,rk3,aj4,aj12
        double precision :: hmed,hmbot,hmtop,hmotop,hmobot
        double precision :: hmntop,hmnbot,rvptop,rvpbot,w,aux
        double precision :: flux_tot,flux_tot1,flux_top,flux_bot
        double precision :: rstot,hn,ho,cdummy,alow,adiag,ahigh
        double precision :: rkmin,rkmax,cconz
      double precision,dimension(3) :: fw,fd,fl,fnudge
      double precision,dimension(3) :: b,c,f,wdiff
      double precision,dimension(0:nlvddi+1) :: hdv,haver,presentl
      double precision,dimension(0:nlvddi+1,3) :: hnew,htnew,rtau,cob
      double precision,dimension(0:nlvddi+1,3) :: hold,htold,vflux,wl
      double precision,dimension(0:nlvddi+1,3) :: cl
      double precision,dimension(0:nlvddi+1,3) :: finu
      double precision,dimension(nlvddi,3) :: clc,clm,clp,cle
	
	if(nlv.ne.nlev) stop 'error stop conz3d_element: nlv/=nlev'

! ----------------------------------------------------------------
!  initialize variables and parameters
! ----------------------------------------------------------------

	btvd = itvd .gt. 0
	bgradup = itvd .eq. 2	!use upwind gradient for tvd scheme
	btvdv = itvdv .gt. 0

! ----------------------------------------------------------------
! global arrays for accumulation of implicit terms
! ----------------------------------------------------------------

! 	 ALLOCATE(fw(3),fd(3),fl(3),fnudge(3),wdiff(3))
! 	 ALLOCATE(b(3),c(3),f(3))
! 	 ALLOCATE(hdv(0:nlvddi+1),haver(0:nlvddi+1))
! 	 ALLOCATE(presentl(0:nlvddi+1))
! 	 ALLOCATE(hnew(0:nlvddi+1,3),htnew(0:nlvddi+1,3))
! 	 ALLOCATE(rtau(0:nlvddi+1,3),cob(0:nlvddi+1,3))
! 	 ALLOCATE(hold(0:nlvddi+1,3),htold(0:nlvddi+1,3))
! 	 ALLOCATE(vflux(0:nlvddi+1,3),wl(0:nlvddi+1,3))
! 	 ALLOCATE(cl(0:nlvddi+1,3))
! 	 ALLOCATE(clc(nlvddi,3),clm(nlvddi,3))
! 	 ALLOCATE(clp(nlvddi,3),cle(nlvddi,3))
	 
          hdv = 0.		!layer thickness
          haver = 0.
	  presentl = 0.		!1. if layer is present
	  hnew = 0.		!as hreal but with zeta_new
	  hold = 0.		!as hreal but with zeta_old
	  cl = 0.		!concentration in layer
	  wl = 0.		!vertical velocity
	  vflux = 0.		!vertical flux
	
!	these are the local arrays for accumulation of implicit terms
!	(maybe we do not need them, but just to be sure...)
!	after accumulation we copy them onto the global arrays

	    cle = 0.
	    clc = 0.
	    clm = 0.
	    clp = 0.
      
	do ii=1,3
          k=nen3v(ii,ie)
	  kn(ii)=k
	  b(ii)=ev(ii+3,ie)
	  c(ii)=ev(ii+6,ie)
	end do

	aj=ev(10,ie)    !area of triangle / 12
	aj4=4.*aj
	aj12=12.*aj
        ilevel=ilhv(ie)

! 	----------------------------------------------------------------
! 	set up vectors for use in assembling contributions
! 	----------------------------------------------------------------

        do l=1,ilevel
	  hdv(l) = hdeov(l,ie)		!use old time step -> FIXME
          !haver(l) = 0.5 * ( hdeov(l,ie) + hdenv(l,ie) )
          haver(l) = rso*hdenv(l,ie) + rsot*hdeov(l,ie)
	  presentl(l) = 1.
	  do ii=1,3
	    k=kn(ii)
	    hn = hdknv(l,k)		! there are never more layers in ie
	    ho = hdkov(l,k)		! ... than in k
            htold(l,ii) = ho
            htnew(l,ii) = hn
	    hold(l,ii) = rso * hn + rsot * ho
	    hnew(l,ii) = rsn * hn + rsnt * ho
	    cl(l,ii) = cn1(l,k)
	    cob(l,ii) = cobs(l,k)	!observations
	    rtau(l,ii) = rtauv(l,k)	!observations
	    wl(l,ii) = wlnv(l,k) - wsink * wsinkv(l,k)
	  end do
	end do

	do l=ilevel+1,nlv
	  presentl(l) = 0.
	end do

! 	----------------------------------------------------------------
! 	set vertical velocities in surface and bottom layer
! 	----------------------------------------------------------------
! 
! 	we do not set wl(0,ii) because otherwise we loose concentration
! 	through surface
! 
! 	we set wl(ilevel,ii) to 0 because we are on the bottom
! 	and there should be no contribution from this element
! 	to the vertical velocity

	do ii=1,3
	  wl(ilevel,ii) = 0.
	end do

! 	----------------------------------------------------------------
! 	compute vertical fluxes (w/o vertical TVD scheme)
! 	----------------------------------------------------------------

	wws = 0.	!sinking already in wl
	call vertical_flux_ie(btvdv,ie,ilevel,dt,wws,cl,wl,hold,vflux)

! ----------------------------------------------------------------
!  loop over levels
! ----------------------------------------------------------------

        do l=1,ilevel

        us=az*utlnv(l,ie)+azt*utlov(l,ie)             !$$azpar
        vs=az*vtlnv(l,ie)+azt*vtlov(l,ie)

        rk3 = 3. * rkpar * difhv(l,ie)

	cbm=0.
	ccm=0.
	itot=0
	isum=0
	do ii=1,3
	  k=kn(ii)
	  f(ii)=us*b(ii)+vs*c(ii)	!$$azpar
	  if(f(ii).lt.0.) then	!flux out of node
	    itot=itot+1
	    isum=isum+ii
	  end if
	  cbm=cbm+b(ii)*cl(l,ii)
	  ccm=ccm+c(ii)*cl(l,ii)

! 	  ----------------------------------------------------------------
! 	  initialization to be sure we are in a clean state
! 	  ----------------------------------------------------------------

	  fw(ii) = 0.
	  !cle(l,ii) = 0.	!ERIC
	  !clc(l,ii) = 0.
	  !clm(l,ii) = 0.
	  !clp(l,ii) = 0.

! 	  ----------------------------------------------------------------
! 	  contributions from horizontal diffusion
! 	  ----------------------------------------------------------------

          waux = 0.
          do iii=1,3
            waux = waux + wdifhv(iii,ii,ie) * cl(l,iii)
          end do
          wdiff(ii) = waux

! 	  ----------------------------------------------------------------
! 	  contributions from vertical diffusion
! 	  ----------------------------------------------------------------
! 
! 	  in fd(ii) is explicit contribution
! 	  the sign is for the term on the left side, therefore
! 	  fd(ii) must be subtracted from the right side
! 
! 	  maybe we should use real layer thickness, or even the
! 	  time dependent layer thickness

	  rvptop = difv(l-1,k) + difmol
	  rvpbot = difv(l,k) + difmol
	  !hmtop = 2. * rvptop * presentl(l-1) / (hdv(l-1)+hdv(l))
	  !hmbot = 2. * rvpbot * presentl(l+1) / (hdv(l)+hdv(l+1))
	  hmotop =2.*rvptop*presentl(l-1)/(hold(l-1,ii)+hold(l,ii))
	  hmobot =2.*rvpbot*presentl(l+1)/(hold(l,ii)+hold(l+1,ii))
	  hmntop =2.*rvptop*presentl(l-1)/(hnew(l-1,ii)+hnew(l,ii))
	  hmnbot =2.*rvpbot*presentl(l+1)/(hnew(l,ii)+hnew(l+1,ii))

	  fd(ii) = adt * ( 
     +			(cl(l,ii)-cl(l+1,ii))*hmobot -
     +			(cl(l-1,ii)-cl(l,ii))*hmotop
     +			  )

	  clc(l,ii) = clc(l,ii) + ad * ( hmntop + hmnbot )
	  clm(l,ii) = clm(l,ii) - ad * ( hmntop )
	  clp(l,ii) = clp(l,ii) - ad * ( hmnbot )

! 	  ----------------------------------------------------------------
! 	  contributions from vertical advection
! 	  ----------------------------------------------------------------
! 
! 	  in fw(ii) is explicit contribution
! 	  the sign is for the term on the left side, therefore
! 	  fw(ii) must be subtracted from the right side
! 
! 	  if we are in last layer, w(l,ii) is zero
! 	  if we are in first layer, w(l-1,ii) is zero (see above)

	  w = wl(l-1,ii)		!top of layer
	  if( l .eq. 1 ) w = 0.		!surface -> no transport (WZERO)
	  if( w .ge. 0. ) then
	    fw(ii) = aat*w*cl(l,ii)
	    flux_top = w*cl(l,ii)
	    clc(l,ii) = clc(l,ii) + aa*w
	  else
	    fw(ii) = aat*w*cl(l-1,ii)
	    flux_top = w*cl(l-1,ii)
	    clm(l,ii) = clm(l,ii) + aa*w
	  end if

	  w = wl(l,ii)			!bottom of layer
	  if( l .eq. ilevel ) w = 0.	!bottom -> handle flux elsewhere (WZERO)
	  if( w .gt. 0. ) then
	    fw(ii) = fw(ii) - aat*w*cl(l+1,ii)
	    flux_bot = w*cl(l+1,ii)
	    clp(l,ii) = clp(l,ii) - aa*w
	  else
	    fw(ii) = fw(ii) - aat*w*cl(l,ii)
	    flux_bot = w*cl(l,ii)
	    clc(l,ii) = clc(l,ii) - aa*w
	  end if

	  flux_tot1 = aat * ( flux_top - flux_bot )
	  flux_tot = aat * ( vflux(l-1,ii) - vflux(l,ii) )

	  fw(ii) = flux_tot
	end do

! 	----------------------------------------------------------------
! 	contributions from horizontal advection (only explicit)
! 	----------------------------------------------------------------
! 
! 	f(ii) > 0 ==> flux into node ii
! 	itot=1 -> flux out of one node
! 		compute flux with concentration of this node
! 	itot=2 -> flux into one node
! 		for flux use conz. of the other two nodes and
! 		minus the sum of these nodes for the flux of this node

	if(itot.eq.1) then	!$$flux
	  fl(1)=f(1)*cl(l,isum)
	  fl(2)=f(2)*cl(l,isum)
	  fl(3)=f(3)*cl(l,isum)
	else if(itot.eq.2) then
	  isum=6-isum
	  fl(1)=f(1)*cl(l,1)
	  fl(2)=f(2)*cl(l,2)
	  fl(3)=f(3)*cl(l,3)
	  fl(isum) = 0.
	  fl(isum) = -(fl(1)+fl(2)+fl(3))
	  isum=6-isum		!reset to original value
	else			!exception	$$itot0
	  fl(1)=0.
	  fl(2)=0.
	  fl(3)=0.
	end if

! 	----------------------------------------------------------------
! 	horizontal TVD scheme start
! 	----------------------------------------------------------------

        if( btvd ) then
	  iext = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( is_external_boundary(k) ) iext = iext + 1
	  end do

          if( iext .eq. 0 ) then
	    call tvd_fluxes(ie,l,itot,isum,dt,cl,cn1,gradxv,gradyv,f,fl)
	  end if
	end if

! 	----------------------------------------------------------------
! 	horizontal TVD scheme finish
! 	----------------------------------------------------------------

! 	----------------------------------------------------------------
! 	contributions from nudging
! 	----------------------------------------------------------------

	do ii=1,3
	  fnudge(ii) = robs * rtau(l,ii) * ( cob(l,ii) - ant * cl(l,ii) )
	  finu(l,ii) = an * robs * rtau(l,ii)	!implicit contribution
	end do

! 	----------------------------------------------------------------
! 	sum explicit contributions
! 	----------------------------------------------------------------

	do ii=1,3
	  k=kn(ii)
          hmed = haver(l)                    !new ggu   !HACK
	  cexpl = aj4 * ( hold(l,ii)*cl(l,ii)
     +				+ dt *  ( 
     +					    hold(l,ii)*fnudge(ii)
     +					  + 3.*fl(ii) 
     +					  - fw(ii)
     +					  - rk3*hmed*wdiff(ii)
     +					  - fd(ii)
     +					)
     +			               )
	  
	  !clm(1,ii) = 0.		!ERIC
	  !clp(ilevel,ii) = 0.
	  ! next check to be deleted
	  if( clm(1,ii) /= 0. .or. clp(ilevel,ii) /= 0. ) then
	    write(6,*) ie,ii,ilevel
	    write(6,*) clm(1,ii),clp(ilevel,ii)
	    stop 'error stop: assumption violated'
	  end if
	  
	  alow  = aj4 * dt * clm(l,ii)
	  ahigh = aj4 * dt * clp(l,ii)
	  adiag = aj4 * dt * clc(l,ii) 
     +			+ aj4 * (1.+dt*finu(l,ii)) * hnew(l,ii)
	  cn(l,k)    = cn(l,k)    + cexpl
	  clow(l,k)  = clow(l,k)  + alow
	  chigh(l,k) = chigh(l,k) + ahigh   
          cdiag(l,k) = cdiag(l,k) + adiag
	end do

	end do		! loop over l
	
! ----------------------------------------------------------------
!  end of loop over l
! ----------------------------------------------------------------

! 	deallocate(fw,fd,fl,fnudge)
! 	deallocate(b,c,f,wdiff)
! 	deallocate(hdv,haver,presentl)
! 	deallocate(hnew,htnew,rtau,cob)
! 	deallocate(hold,htold,vflux,wl,cl)
! 	deallocate(clc,clm,clp,cle)
! 	
! ----------------------------------------------------------------
!  end of routine
! ----------------------------------------------------------------

      end subroutine conz3d_element

! *****************************************************************
      
       subroutine conz3d_nodes(k,cn,cdiag,clow,chigh,cn1,cbound,
     +                         load,rload,ad,aa,dt,nlvddi)

      	use mod_bound_geom
	use mod_geom
	use mod_depth
	use mod_diff_aux
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi
	
	implicit none
	
	integer,intent(in) :: k,nlvddi
	real,intent(in) :: rload
	real,dimension(nlvddi,nkn),intent(in) :: cn1,cbound
	real,dimension(nlvddi,nkn),intent(inout) :: load 		!LLL
	double precision, intent(in) :: dt
	double precision, intent(in) :: ad,aa

	double precision,dimension(nlvddi,nkn),intent(inout) :: cn
	double precision,dimension(nlvddi),intent(inout) :: cdiag
	double precision,dimension(nlvddi),intent(inout) :: clow
	double precision,dimension(nlvddi),intent(inout) :: chigh

	integer :: l,ilevel,lstart,i,ii,ie,n,ibase
	double precision :: mflux,qflux,cconz
	double precision :: loading,aux,cload

	double precision, parameter :: d_tiny = tiny(1.d+0)
	double precision, parameter :: r_tiny = tiny(1.)
      
! ----------------------------------------------------------------
!  handle boundary (flux) conditions
! ----------------------------------------------------------------

      	  ilevel = ilhkv(k)

	  do l=1,ilevel

            !mflux = cbound(l,k)		!mass flux has been passed
	    cconz = cbound(l,k)			!concentration has been passed
	    qflux = mfluxv(l,k)
	    if( qflux .lt. 0. .and. is_boundary(k) ) cconz = cn1(l,k)
	    mflux = qflux * cconz

            cn(l,k) = cn(l,k) + dt * mflux	!explicit treatment

	    loading = rload*load(l,k)
            if( loading == 0 ) then			!no loading
              !nothing
            else if ( loading < 0.d0 ) then		!excess deposition
	      cload = 0.
	      if( cn(l,k) > 0. ) then
                cload = - dt * loading
                cload = cn(l,k) * ( 1. - exp(-cload/cn(l,k)) )
	      else if( cn(l,k) < 0. ) then
	        cn(l,k) = 0.
	      end if
              if( cload > cn(l,k) ) goto 98
              loading = -cload / dt
              if( rload > 0. ) load(l,k) = loading / rload
              cn(l,k) = cn(l,k) + dt*loading
            else					!erosion
              cn(l,k) = cn(l,k) + dt*loading
            end if
 
	  end do

! ----------------------------------------------------------------
!  compute concentration for each node (solve system)
! ----------------------------------------------------------------

	if((aa .eq. 0. .and. ad .eq. 0.).or.(nlv .eq. 1)) then

	  if( nlv .gt. 1 ) then
	    write(6,*) 'conz: computing explicitly ',nlv
	  end if

	  ilevel = ilhkv(k)
	  do l=1,ilevel
	    if(cdiag(l).ne.0.) then
	      cn(l,k)=cn(l,k)/cdiag(l)
	    end if
	  end do

	else

	  ilevel = ilhkv(k)
	  aux=1./cdiag(1)
	  chigh(1)=chigh(1)*aux
	  cn(1,k)=cn(1,k)*aux
	  do l=2,ilevel
	    if( cdiag(l) == 0. ) goto 99
	    aux=1./(cdiag(l)-clow(l)*chigh(l-1))
	    chigh(l)=chigh(l)*aux
	    cn(l,k)=(cn(l,k)-clow(l)*cn(l-1,k))*aux
	  end do
	  lstart = ilevel-1
	  do l=lstart,1,-1	!$$LEV0 bug 14.08.1998 -> ran to 0
	    cn(l,k)=cn(l,k)-cn(l+1,k)*chigh(l)
	  end do
	end if
	
! ----------------------------------------------------------------
!  end of routine
! ----------------------------------------------------------------

	return
   98	continue
	write(6,*) 'error computing loading: ',l,k
	write(6,*) 'loading,cload,cn(l,k): ',loading,cload,cn(l,k)
	stop 'error stop conz3d_nodes: internal error (1)'
   99	continue
	write(6,*) k,l,ilevel
	write(6,*) nkn_inner,nkn_local,nkn
	write(6,*) cdiag(l),clow(l),chigh(l-1)
	stop 'error stop conz3d_nodes (omp): diag == 0'
      end subroutine conz3d_nodes

c*****************************************************************

