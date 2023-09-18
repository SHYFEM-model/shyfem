!
! revision log :
!
! 20.05.2015    erp     transformed for OMP
! 30.09.2015    ggu     routine cleaned, no reals in conz3d
! 20.11.2015    ggu&erp chunk size introduced, omp finalized
!
!**************************************************************
!------------------------------------------------------------------------------------------------
        module concentration_omp
!------------------------------------------------------------------------------------------------
        contains
!------------------------------------------------------------------------------------------------

        subroutine conz3d_omp(cn1,co1,ddt,rkpar,difhv,difv,difmol,cbound,itvd,itvdv,gradxv,gradyv       &
     &			,cobs,robs,wsink,wsinkv,rload,load,azpar,adpar,aapar,istot,isact,nlvddi,nlev)
     
! computes concentration
!
! cn     new concentration
! co     old concentration              !not used !FIXME
! caux   aux vector
! clow	 lower diagonal of vertical system
! chig	 upper diagonal of vertical system
! ddt    time step
! rkpar  horizontal turbulent diffusivity
! difhv  horizontal turbulent diffusivity (variable between elements)
! difv   vertical turbulent diffusivity
! difmol vertical molecular diffusivity
! cbound boundary condition (mass flux) [kg/s] -> now concentration [kg/m**3]
! itvd	 type of horizontal transport algorithm used
! itvdv	 type of vertical transport algorithm used
! gradxv,gradyv  gradient vectors for TVD algorithm
! cobs	 observations for nudging
! robs	 use observations for nuding (double precision)
! wsink	 factor for settling velocity
! wsinkv variable settling velocity [m/s]
! rload	 factor for loading
! load   load (source or sink) [kg/s]
! azpar  time weighting parameter
! adpar  time weighting parameter for vertical diffusion (ad)
! aapar  time weighting parameter for vertical advection (aa)
! istot	 total inter time steps
! isact	 actual inter time step
! nlvddi	 dimension in z direction
! nlv	 actual needed levels
!
! written 09.01.94 by ggu  (from scratch)
! revised 19.01.94 by ggu  $$flux - flux conserving property
! revised 20.01.94 by ggu  $$iclin - iclin not used to compute volume
! revised 20.01.94 by ggu  $$lumpc - evaluate conz. nodewise
! revised 03.02.94 by ggu  $$itot0 - exception for itot=0 or 3
! revised 04.02.94 by ggu  $$fact3 - factor 3 missing in transport
! revised 04.02.94 by ggu  $$azpar - azpar used to compute transport
! revised 04.02.94 by ggu  $$condry - comute conz also in dry areas
! revised 07.02.94 by ggu  $$istot - istot for fractional time step
! revised 01.06.94 by ggu  restructured for 3-d model
! revised 18.07.94 by ggu  $$htop - use htop instead of htopo for mass cons.
! revised 09.04.96 by ggu  $$rvadj adjust rv in certain areas
!
! solution of purely diffusional part :
!
! dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
!
! C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )
!
! for n-dimensions and
!
! C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )
!
! for 1 dimension
!
! the solution is normalized, i.e.  int(C(x,t)dx) = 1 over the whole area
!
! DPGGU -> introduced double precision to stabilize solution

	use bnd_geom
	use geom
	use depth
	use diff_aux
	use bnd_dynamic
	use area
	use ts
	use hydro_vel
	use hydro_admin
	use evgeom
	use levels
	use basin
!$	use omp_lib
	use subset
	use shympi
        use mpi_io_admin
	use layer_thickness
        use timing

	implicit none

	integer, intent(in) :: nlvddi,nlev,itvd,itvdv,istot,isact
	double precision, intent(in) :: difmol,robs,wsink,rload,ddt,rkpar
	double precision, intent(in) :: azpar,adpar,aapar
	double precision,dimension(nlvddi,nkn),intent(in) :: co1,cbound
	double precision,dimension(nlvddi,nkn),intent(in) :: load
        !double precision,dimension(nlvddi,nkn),intent(out) :: cn
        
	logical :: btvdv
	integer :: ie,k,ilevel,ibase,ii,l,n,i,j,x,ies,iend,kl,kend
	integer :: myid,numthreads,j_init,j_end,knod,k_end,jel
	integer,allocatable,dimension(:) :: subset_l
	double precision :: time1,time2
	double precision :: dtime1,dtime2
	integer :: nchunk,nthreads,nelems,nnodes
	double precision :: dt
	double precision :: az,ad,aa,azt,adt,aat
	double precision :: rstot,rso,rsn,rsot,rsnt
	double precision :: timer,timer1,chunk,rest
	
! 	double precision,dimension(nlvddi,nkn) :: cn
! 	double precision,dimension(nlvddi,nkn) :: co        
!         double precision,dimension(nlvddi,nkn) :: cdiag
! 	double precision,dimension(nlvddi,nkn) :: clow
! 	double precision,dimension(nlvddi,nkn) :: chigh
	
	double precision,dimension(:,:),allocatable :: cn        
        double precision,dimension(:,:),allocatable :: cdiag
	double precision,dimension(:,:),allocatable :: clow
	double precision,dimension(:,:),allocatable :: chigh
#ifdef DEBUGON
        integer e
	double precision,dimension(nlvddi,nkn_local),intent(in) :: cobs
	double precision,dimension(nlvddi,nkn_local),intent(in) :: gradxv,gradyv
	double precision,dimension(0:nlvddi,nkn_local),intent(in) :: difv,wsinkv
	double precision,dimension(nlvddi,nel_local),intent(in) :: difhv
	double precision,dimension(nlvddi,nkn_local),intent(inout) :: cn1
#else
	double precision,dimension(nlvddi,nkn),intent(in) :: cobs
	double precision,dimension(nlvddi,nkn),intent(in) :: gradxv,gradyv
	double precision,dimension(0:nlvddi,nkn),intent(in) :: difv,wsinkv
	double precision,dimension(nlvddi,nel),intent(in) :: difhv
	double precision,dimension(nlvddi,nkn),intent(inout) :: cn1
#endif

        if(nlv.ne.nlev) stop 'error stop conzstab: level'
	
!----------------------------------------------------------------
! initialize variables and parameters
!----------------------------------------------------------------

!	call cpu_time(time1)
!!$	dtime1 = omp_get_wtime()
	

#ifdef DEBUGON
        call shympi_exchange_halo_3d0_nodes(wlnv)
        call shympi_exchange_halo_3d0_nodes(difv)
        call shympi_exchange_halo_3d0_nodes(wsinkv)
        call shympi_exchange_halo_2d_elems(ilhv)
        call shympi_exchange_halo_3d_elems(hdeov)
        call shympi_exchange_halo_3d_elems(hdenv)
        call shympi_exchange_halo_3d_nodes(hdknv)
        call shympi_exchange_halo_3d_nodes(hdkov)
        call shympi_exchange_halo_3d_elems(utlnv)
        call shympi_exchange_halo_3d_elems(vtlnv)
        call shympi_exchange_halo_3d_elems(utlov)
        call shympi_exchange_halo_3d_elems(vtlov)
        call shympi_exchange_halo_3d_elems(difhv)
        call shympi_exchange_halo_4d_elems(3,3,wdifhv)
!        call shympi_exchange_halo_3d_nodes(mfluxv)
        call shympi_exchange_halo_3d_nodes(rtauv)
        call shympi_exchange_halo_3d_nodes(cobs)
        call shympi_exchange_halo_3d_nodes(cn1)
        call shympi_exchange_halo_3d_nodes(gradxv)
        call shympi_exchange_halo_3d_nodes(gradyv)
        call shympi_exchange_halo_2d_nodes(iopbnd)
        call shympi_exchange_halo_3d_elems(ulnv)
        call shympi_exchange_halo_3d_elems(vlnv)
	ALLOCATE(cn(nlvddi,nkn_local))
	ALLOCATE(cdiag(nlvddi,nkn_local))
	ALLOCATE(clow(nlvddi,nkn_local))
	ALLOCATE(chigh(nlvddi,nkn_local))
#else
	ALLOCATE(cn(nlvddi,nkn))
	ALLOCATE(cdiag(nlvddi,nkn))
	ALLOCATE(clow(nlvddi,nkn))
	ALLOCATE(chigh(nlvddi,nkn))
#endif

	az = azpar
	ad = adpar
	aa = aapar
	
	azt=1.-az
	adt=1.-ad
	aat=1.-aa

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
        cdiag=0.
        clow=0.
        chigh=0.

        nchunk = 1
	nthreads = 1
!$	nthreads = omp_get_num_threads()

#ifdef DEBUGON
         do e=1,nel_local
            if(bmpi) then
               ie=domain%elems%mapID(e)
            else
               ie=e
            end if
#else
      do i=1,subset_num 	! loop over indipendent subset
       
!$     nchunk = subset_el(i) / ( nthreads * 10 )
       nchunk = max(nchunk,1)

!$OMP TASKGROUP 
       do jel=1,subset_el(i),nchunk

!$OMP TASK FIRSTPRIVATE(jel,i) DEFAULT(NONE)                       &
!$OMP   & PRIVATE(j,ie)                                            &
!$OMP   & SHARED(nlvddi,nlev,itvd,itvdv,istot,isact,aa,nchunk)     &
!$OMP   & SHARED(difmol,robs,wsink,rload,ddt,rkpar,az,ad)          &
!$OMP   & SHARED(azt,adt,aat,rso,rsn,rsot,rsnt,dt,nkn)             &
!$OMP   & SHARED(cn,co,cdiag,clow,chigh,subset_el,cn1,co1)         &
!$OMP   & SHARED(subset_num,indipendent_subset)                    &
!$OMP   & SHARED(difhv,cbound,gradxv,gradyv,cobs,load,difv,wsinkv) 

       do j=jel,jel+nchunk-1 	! loop over elements in subset
		if(j .le. subset_el(i)) then
	        ie = indipendent_subset(j,i)
	        !print *,i,ie
#endif
                call conz3d_element(ie,cdiag,clow,chigh,cn,cn1,dt,rkpar,difhv,difv,difmol       &
     &			,itvd,itvdv,gradxv,gradyv,cobs,robs,wsink,wsinkv                 &             
     &                  ,az,ad,aa,azt,adt,aat,rso,rsn,rsot,rsnt,nlvddi,nlev)
#ifndef DEBUGON
		end if
	end do ! end loop over el in subset
!$OMP END TASK
      end do

!$OMP END TASKGROUP       
#endif
       end do ! end loop over subset

#ifndef DEBUGON
       if( shympi_partition_on_elements() ) then
          if(ln_timing) time1 = shympi_wtime()
         call shympi_exchange_and_sum_3d_nodes(cn)
         call shympi_exchange_and_sum_3d_nodes(cdiag)
         call shympi_exchange_and_sum_3d_nodes(clow)
         call shympi_exchange_and_sum_3d_nodes(chigh)
         if(ln_timing) then
            comm_scalar_time = comm_scalar_time + shympi_wtime() - time1
            !comm_ts_time = comm_ts_time + shympi_wtime() - time1
         end if
       end if
#endif




!$     nchunk = nkn / ( nthreads * 10 )
       nchunk = max(nchunk,1)

!$OMP TASKGROUP
       do knod=1,nkn,nchunk
!$OMP TASK FIRSTPRIVATE(knod) PRIVATE(k) DEFAULT(NONE)          &
!$OMP   & SHARED(cn,cdiag,clow,chigh,cn1,cbound,load,nchunk,    &
!$OMP   &           rload,ad,aa,dt,nlvddi,nkn)
	 do k=knod,knod+nchunk-1
	 if(k .le. nkn) then
	   call conz3d_nodes(k,cn,cdiag(:,k),clow(:,k),chigh(:,k),cn1,cbound,load,rload,ad,aa,dt,nlvddi)
         endif
         enddo
!$OMP END TASK 	      
	end do

!$OMP END TASKGROUP

	!cn1 = 0.
	!cn1 = cn
	!cn1 = double precision(cn)
	cn1 = cn
	
	DEALLOCATE(cn)
	DEALLOCATE(cdiag)
	DEALLOCATE(clow)
	DEALLOCATE(chigh)
	
!	call cpu_time(time2)
!!$	dtime2 = omp_get_wtime()
!	write(6,*) time2-time1,dtime2-dtime1

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!*****************************************************************

       subroutine conz3d_element(ie,cdiag,clow,chigh,cn,cn1,dt,rkpar,difhv,difv,difmol   &
     &			,itvd,itvdv,gradxv,gradyv,cobs,robs,wsink,wsinkv,az,ad,aa,azt,adt,aat   &
     &			,rso,rsn,rsot,rsnt,nlvddi,nlev)
     
        use bnd_geom
        use geom
        use depth
        use diff_aux
        use bnd_dynamic
        use area
        use ts
        use hydro_vel
        use hydro_admin
        use evgeom
        use levels
        use basin
        use layer_thickness
        use tvd_admin
        use mpi_common_struct
      
        implicit none
      
        integer,intent(in) :: ie,nlvddi,nlev,itvd,itvdv
        double precision,intent(in) :: difmol,robs,wsink,rkpar
        double precision,intent(in) :: dt
        double precision,intent(in) :: az,ad,aa,azt,adt,aat
        double precision,intent(in) :: rso,rsn,rsot,rsnt
        
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
        double precision,dimension(nlvddi,3) :: clc,clm,clp,cle
#ifdef DEBUGON
        double precision,dimension(nlvddi,nkn_local),intent(in) :: cn1
        double precision,dimension(nlvddi,nkn_local),intent(in) :: cobs
        double precision,dimension(nlvddi,nkn_local),intent(in) :: gradxv,gradyv
        double precision,dimension(nlvddi,nel_local),intent(in) :: difhv
        double precision,dimension(nlvddi,nkn_local),intent(inout) :: cdiag
        double precision,dimension(nlvddi,nkn_local),intent(inout) :: clow
        double precision,dimension(nlvddi,nkn_local),intent(inout) :: chigh
        double precision,dimension(nlvddi,nkn_local),intent(inout) :: cn
        double precision,intent(in),dimension(0:nlvddi,nkn_local) :: wsinkv,difv
#else
        double precision,dimension(nlvddi,nkn),intent(in) :: cn1
        double precision,dimension(nlvddi,nkn),intent(in) :: cobs
        double precision,dimension(nlvddi,nkn),intent(in) :: gradxv,gradyv
        double precision,dimension(nlvddi,nel),intent(in) :: difhv
        double precision,dimension(nlvddi,nkn),intent(inout) :: cdiag
        double precision,dimension(nlvddi,nkn),intent(inout) :: clow
        double precision,dimension(nlvddi,nkn),intent(inout) :: chigh
        double precision,dimension(nlvddi,nkn),intent(inout) :: cn
        double precision,intent(in),dimension(0:nlvddi,nkn) :: wsinkv,difv
#endif
        
        if(nlv.ne.nlev) stop 'error stop conzstab: level'

! ----------------------------------------------------------------
!  initialize variables and parameters
! ----------------------------------------------------------------

        btvd    = itvd .gt. 0
        bgradup = itvd .eq. 2   !use upwind gradient for tvd scheme
        btvdv   = itvdv .gt. 0

        wws = 0.

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
         
        hdv = 0.        !layer thickness
        haver = 0.
        presentl = 0.   !1. if layer is present
        hnew = 0.       !as hreal but with zeta_new
        hold = 0.       !as hreal but with zeta_old
        cl = 0.         !concentration in layer
        wl = 0.         !vertical velocity
        vflux = 0.      !vertical flux
        
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
          hdv(l) = hdeov(l,ie)          !use old time step -> FIXME
          !haver(l) = 0.5 * ( hdeov(l,ie) + hdenv(l,ie) )
          haver(l) = rso*hdenv(l,ie) + rsot*hdeov(l,ie)
          presentl(l) = 1.
          do ii=1,3
            k=kn(ii)
            hn = hdknv(l,k)             ! there are never more layers in ie
            ho = hdkov(l,k)             ! ... than in k
            htold(l,ii) = ho
            htnew(l,ii) = hn
            hold(l,ii) = rso * hn + rsot * ho
            hnew(l,ii) = rsn * hn + rsnt * ho
            cl(l,ii) = cn1(l,k)
            cob(l,ii) = cobs(l,k)       !observations
            rtau(l,ii) = rtauv(l,k)     !observations
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
          f(ii)=us*b(ii)+vs*c(ii)       !$$azpar
          if(f(ii).lt.0.) then          !flux out of node
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
! 	  maybe we should use double precision layer thickness, or even the
! 	  time dependent layer thickness

          rvptop = difv(l-1,k) + difmol
          rvpbot = difv(l,k) + difmol
          !hmtop = 2. * rvptop * presentl(l-1) / (hdv(l-1)+hdv(l))
          !hmbot = 2. * rvpbot * presentl(l+1) / (hdv(l)+hdv(l+1))
          hmotop =2.*rvptop*presentl(l-1)/(hold(l-1,ii)+hold(l,ii))
          hmobot =2.*rvpbot*presentl(l+1)/(hold(l,ii)+hold(l+1,ii))
          hmntop =2.*rvptop*presentl(l-1)/(hnew(l-1,ii)+hnew(l,ii))
          hmnbot =2.*rvpbot*presentl(l+1)/(hnew(l,ii)+hnew(l+1,ii))

          fd(ii) = adt * ((cl(l,ii)-cl(l+1,ii))*hmobot - (cl(l-1,ii)-cl(l,ii))*hmotop)

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

          w = wl(l-1,ii) - wws		!top of layer
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

          w = wl(l,ii) - wws            !bottom of layer
          if( l .eq. ilevel ) w = 0.    !bottom -> handle flux elsewhere (WZERO)
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
          isum=6-isum           !reset to original value
        else                    !exception	$$itot0
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
          fnudge(ii) = robs * rtau(l,ii) * ( cob(l,ii) - cl(l,ii) )
        end do

! 	----------------------------------------------------------------
! 	sum explicit contributions
! 	----------------------------------------------------------------

        do ii=1,3
          k=kn(ii)
          hmed = haver(l)                    !new ggu   !HACK
          cexpl = aj4 * ( hold(l,ii)*cl(l,ii)+ dt * ( hold(l,ii)*fnudge(ii)     &
     &                  + 3.*fl(ii) - fw(ii) - rk3*hmed*wdiff(ii) - fd(ii)))
          
          !clm(1,ii) = 0.		!ERIC
          !clp(ilevel,ii) = 0.
          ! next check to be deleted
          if( (clm(1,ii) .ne. 0.) .or. (clp(ilevel,ii) .ne. 0.) ) then
            write(6,*) ie,ii,ilevel
            write(6,*) clm(1,ii),clp(ilevel,ii)
            stop 'error stop: assumption violated'
          end if
          
          alow  = aj4 * dt * clm(l,ii)
          ahigh = aj4 * dt * clp(l,ii)
          adiag = aj4 * dt * clc(l,ii) + aj4 * hnew(l,ii)
          cn(l,k)    = cn(l,k)    + cexpl
          clow(l,k)  = clow(l,k)  + alow
          chigh(l,k) = chigh(l,k) + ahigh   
          cdiag(l,k) = cdiag(l,k) + adiag
           
        end do          ! loop over 3 nodes

        end do          ! loop over l
        
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
      
       subroutine conz3d_nodes(k,cn,cdiag,clow,chigh,cn1,cbound,load,rload,ad,aa,dt,nlvddi)

      	use bnd_geom
	use geom
	use depth
	use diff_aux
	use bnd_dynamic
	use area
	use ts
	use hydro_vel
	use hydro_admin
	use evgeom
	use levels
	use basin
	
	implicit none
	
	integer,intent(in) :: k,nlvddi
	double precision,intent(in) :: rload
	double precision,dimension(nlvddi,nkn),intent(in) :: cn1,cbound,load
	double precision, intent(in) :: dt
	double precision, intent(in) :: ad,aa

	double precision,dimension(nlvddi,nkn),intent(inout) :: cn
	double precision,dimension(nlvddi),intent(inout) :: cdiag
	double precision,dimension(nlvddi),intent(inout) :: clow
	double precision,dimension(nlvddi),intent(inout) :: chigh

	integer :: l,ilevel,lstart,i,ii,ie,n,ibase
	double precision :: mflux,qflux,cconz
	double precision :: loading,aux

	double precision, parameter :: d_tiny = tiny(1.d+0)
	double precision, parameter :: r_tiny = tiny(1.)
      
! ----------------------------------------------------------------
!  handle boundary (flux) conditions
! ----------------------------------------------------------------

      	  ilevel = ilhkv(k)

	  do l=1,ilevel
            !mflux = cbound(l,k)		!mass flux has been passed
	    cconz = cbound(l,k)		!concentration has been passed
	    qflux = mfluxv(l,k)
	    if( qflux .lt. 0. .and. is_boundary(k) ) cconz = cn1(l,k)
	    mflux = qflux * cconz

            cn(l,k) = cn(l,k) + dt * mflux	!explicit treatment

	    loading = rload*load(l,k)
            if( loading .eq. 0. ) then
	      !nothing
	    else if( loading .gt. 0. ) then    		!treat explicit
              cn(l,k) = cn(l,k) + dt * loading
            else !if( loading .lt. 0. ) then		!treat quasi implicit
	      if( cn1(l,k) .gt. 0. ) then
                cdiag(l) = cdiag(l) - dt * loading/cn1(l,k)
	      end if
            end if
	  end do

! ----------------------------------------------------------------
!  compute concentration for each node (solve system)
! ----------------------------------------------------------------

	if((aa .eq. 0. .and. ad .eq. 0.).or.(nlv .eq. 1)) then

	if( nlv .gt. 1 ) then
	  write(6,*) 'conz: computing explicitly ',nlv
	end if

	!do k=1,nkn
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

      end subroutine conz3d_nodes

!*****************************************************************

!------------------------------------------------------------------------------------------------
        end module concentration_omp
!------------------------------------------------------------------------------------------------
