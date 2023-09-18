!
! $Id: newstab.f,v 1.3 2010-03-08 17:46:45 georg Exp $
!! subroutine conzstab(cn1,co1
!     +                  ,ddt
!     +                  ,robs,wsink,wsinkv
!     +                  ,rkpar,difhv,difv
!     +                  ,difmol,azpar
!     +                  ,adpar,aapar
!     +                  ,sindex
!     +                  ,istot,isact
!     +                  ,nlvddi,nlv)
!
! routines for stability computations
!
! revision log :
!
! 11.10.2002    ggu     con3sh removed, conzstab better commented
! 02.12.2004    ggu     return also sindex in conzstab
! 11.11.2008    ggu     conzstab cleaned
! 06.12.2008    ggu     in conzstab changed wprv, new routine write_elem_info()
! 19.02.2010    ggu     new file to contain stability computations
! 26.02.2010    ggu     internal_stability restructured (on element)
! 26.02.2010	ggu	new call to momentum_viscous_stability()
! 26.02.2010	ggu	new momentum_advective_stability()
! 08.03.2010    ggu     run only down to avail layers in stability (bug fix)
! 26.01.2011    ggu     robs for nudging implemented
! 16.02.2011    ggu     pass robs to subroutines, write stb-ind to nos file
! 20.05.2011    ggu     allow for elimination of elems due to high rstab
! 01.06.2011    ggu     wsink for stability integrated
! 12.07.2011    ggu     new routine output_stability()
! 14.07.2011    ggu     new routine output_stability_node()
! 21.06.2012    ggu&ccf variable vertical sinking velocity integrated
! 08.04.2014    ggu	use rlin to determine advective stability
! 20.05.2015    ggu	always compute stability, call to conzstab changed
! 20.10.2015    ggu	in output_stability() bug that icall was not adjouned
!
!*****************************************************************
!*****************************************************************
!*****************************************************************
!
! typical usage :
!
! call reset_stability at beginning of time loop
! call make_stability before advection of similar variables
!
! notes :
!
!	scal3sh							newcon
!		info_stability					newstab
!		make_stability					newstab
!			compute_stability			newstab
!				conzstab			newcon
!
!	set_timestep						subtime
!		hydro_stability					newstab
!			internal_stability			newstab
!				momentum_advective_stability	newexpl
!				momentum_viscous_stability	newexpl
!
!*****************************************************************
!--------------------------------------------------------------------------------
        module stability
!--------------------------------------------------------------------------------
        contains
!--------------------------------------------------------------------------------

	subroutine compute_stability(robs,wsink,wsinkv,rkpar,azpar,rindex,saux)

! computes stability index

	!use conz_common
	use diffusion
	use levels, only : nlvdi,nlv
	use basin
        use para
#ifdef DEBUGON
        use mpi_common_struct
#endif

	implicit none

        include 'param.h'

	double precision robs
	double precision wsink
#ifdef DEBUGON
	double precision wsinkv(0:nlvdi,nkn_local)
#else
	double precision wsinkv(0:nlvdi,nkn)
#endif
        double precision rkpar
        double precision azpar
        double precision rindex
        double precision saux(nlvdi,nkn)

        double precision adpar,aapar
        double precision difmol
	double precision ddt
        integer isact,istot

!----------------------------------------------------------------
! set parameters
!----------------------------------------------------------------

	adpar=getpar('adpar')
	aapar=getpar('aapar')

        isact = 1
        istot = 1
        difmol = 0.d0
	ddt = 1.d0		!always for 1 sec
	rindex = 0.d0

!----------------------------------------------------------------
! call conzstab
!----------------------------------------------------------------

        !call conzstab(cnv,saux
!     +          ,ddt,robs,wsink,wsinkv,rkpar,difhv,difv
        call conzstab(ddt,robs,wsink,wsinkv,rkpar,difhv,difv,difmol,azpar       &
     &                  ,adpar,aapar,rindex,istot,isact,nlvdi,nlv)

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!*****************************************************************

        subroutine make_stability(dt,robs,wsink,wsinkv,rkpar,rindex,istot,saux)

! gets stability index (if necessary computes it)

	use levels, only : nlvdi,nlv
	use basin
        use para
#ifdef DEBUGON
        use mpi_common_struct
#endif

        implicit none

	include 'param.h'

	double precision dt
	double precision robs
	double precision wsink
#ifdef DEBUGON
	double precision wsinkv(0:nlvdi,nkn_local)
#else
	double precision wsinkv(0:nlvdi,nkn)
#endif
        double precision rkpar
        double precision rindex
        integer istot
	double precision saux(nlvdi,nkn)

	double precision azpar
	logical exist_stability

!----------------------------------------------------------------
! see if already computed (only for wsink == 0)
!----------------------------------------------------------------

	!if( wsink .eq. 0. .and. exist_stability(rkpar,rindex) ) goto 1
	  
!----------------------------------------------------------------
! compute stability index
!----------------------------------------------------------------

	call getaz(azpar)
	call compute_stability(robs,wsink,wsinkv,rkpar,azpar,rindex,saux)

!----------------------------------------------------------------
! insert stability index (only for wsink == 0)
!----------------------------------------------------------------

	!if( wsink .eq. 0. ) call insert_stability(rkpar,rindex)

!----------------------------------------------------------------
! scale to double precision time step dt
!----------------------------------------------------------------

    1	continue
	rindex = dt * rindex
	istot = 1 + rindex

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

        end

!*****************************************************************

        subroutine info_stability(dt,robs,wsink,wsinkv,rkpar,rindex,istot,saux)

! gets stability index (if necessary computes it)

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
        use conz_util
        use para
#ifdef DEBUGON
        use mpi_common_struct
#endif

        implicit none

	include 'param.h'

        double precision dt
	double precision robs
	double precision wsink
#ifdef DEBUGON
	double precision wsinkv(0:nlvdi,nkn_local)
#else
	double precision wsinkv(0:nlvdi,nkn)
#endif
        double precision rkpar
        double precision rindex
        integer istot
	double precision saux(nlvdi,nkn)

	integer ia,iustab
	integer l,k
	double precision aindex
	double precision azpar
	logical exist_stability

!----------------------------------------------------------------
! compute stability index
!----------------------------------------------------------------

	call getaz(azpar)
	call compute_stability(robs,wsink,wsinkv,rkpar,azpar,rindex,saux)
	rindex = dt * rindex
	istot = 1 + rindex

!----------------------------------------------------------------
! write to terminal
!----------------------------------------------------------------

	ia = 1
	aindex = saux(1,1)

	do k=1,nkn
	  do l=1,nlv
	    if( saux(l,k) .gt. aindex ) then
	      aindex = saux(l,k)
	      ia = k
	    end if
	  end do
	end do

!ggu protect
	write(6,*) 'info_stability'
	write(6,*) rkpar,azpar,rindex,istot
	write(6,*) ia,aindex,dt*aindex
	iustab = 0
	call conwrite(iustab,'.sta',1,778,nlvdi,saux)
!ggu protect

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine reset_stability

! rests stability index

        implicit none

	include 'stab.h'

	nentry = 0

        end

!**********************************************************************

	subroutine insert_stability(rkpar,rindex)

! inserts stability index

	implicit none

	double precision rkpar,rindex

	include 'stab.h'

	integer i

	do i=1,nentry
	  if( rkind(1,i) .eq. rkpar ) return
	end do
	if( i .gt. ndim_stab ) goto 99

	!return		!uncomment for debug -> never insert

!ggu protect
	nentry = i
	rkind(1,i) = rkpar
	rkind(2,i) = rindex
!ggu protect

	return
   99	continue
	write(6,*) 'nentry,ndim: ',i,ndim_stab
	stop 'error stop insert_stability: ndim'
	end

!**********************************************************************

	function exist_stability(rkpar,rindex)

! tests if stability index has already been computed and returns it

	implicit none

	logical exist_stability
	double precision rkpar,rindex

	include 'stab.h'

	integer i

	do i=1,nentry
	  if( rkind(1,i) .eq. rkpar ) goto 1
	end do
    1	continue

	exist_stability = i .le. nentry

	if( exist_stability ) then
	  rindex = rkind(2,i)
	else
	  rindex = 0
	end if

	end

!**********************************************************************
!**********************************************************************
!**********************************************************************

        subroutine hydro_stability(dt,rindex)

! computes stability index for hydro timestep - no error output

        implicit none

        double precision dt
        double precision rindex

	call internal_stability(0,dt,rindex)

	end

!**********************************************************************

        subroutine error_stability(dt,rindex)

! computes stability index for hydro timestep - with error output
!
! after this call the program should abort

        implicit none

        double precision dt
        double precision rindex

	call internal_stability(1,dt,rindex)

	end

!**********************************************************************

        subroutine eliminate_stability(rmax)

! eliminates elements with stability index higher than rmax

        implicit none

        double precision rmax

	integer mode
        double precision rindex,dt

	mode = 2
	dt = 0.
	rindex = rmax
	call internal_stability(mode,dt,rindex)

	end

!**********************************************************************

        subroutine internal_stability(mode,dt,rindex)

! computes stability index for hydro timestep (internal)
!
! mode = 0		normal call, compute stability
! mode = 1		error call, compute stability and write error message
! mode = 2		eliminate elements with r>rindex

	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi
        use para
        use wetdry
        use openmp_admin
        use check

        implicit none

	include 'param.h'

	integer mode		!0: normal call  1:error output
        double precision dt			!time step to be used
        double precision rindex		!stability index (return)

	integer ie,l,lmax,iweg
        double precision rkpar,azpar,ahpar,rlin
	double precision dindex,aindex,tindex,sindex
	double precision rmax


	!double precision, allocatable :: sauxe1(:,:)
	!double precision, allocatable :: sauxe2(:,:)
        double precision, allocatable :: sauxe1(:,:)
        double precision, allocatable :: sauxe2(:,:)

	logical is_i_nan

        rkpar = 0.d0
	azpar = 1.d0
	ahpar = getpar('ahpar')
	rlin = getpar('rlin')

	allocate(sauxe1(nlvdi,nel),sauxe2(nlvdi,nel))
	sauxe1 = 0.
	sauxe2 = 0.

	rmax = 1.d+30
	if( mode .eq. 2 ) rmax = rindex
	if( mode .eq. 2 ) then
		write(6,*) 'eliminating rmax: ',rmax
	end if

	call momentum_advective_stability(rlin,aindex,sauxe1)
	call momentum_viscous_stability(ahpar,dindex,sauxe2)

	if( .not. openmp_in_parallel() ) then
          call output_stability(dt,sauxe1,sauxe2)	!in case write to file
	end if

	aindex = aindex*dt
	dindex = dindex*dt

	tindex = 0.d0
	do ie=1,nel
	  lmax = ilhv(ie)
	  iweg = 0
	  do l=1,lmax
	    sindex = sauxe1(l,ie)+sauxe2(l,ie)
	    if( sindex .ge. rmax ) iweg = 1
	    tindex = max(tindex,sindex)
	  end do
	  if( iweg .gt. 0 ) then
	    write(6,*) 'eliminating element for stability: ',ie
	    write(569,*) 'eliminating element for stability: ',ie
	    call check_set_unit(570)
	    call check_elem(ie)
	    call set_element_dry(ie)
	  end if
	end do

	tindex = shympi_max(tindex)

	tindex = tindex*dt

	if( mode .eq. 1 ) then		!error output
	  write(6,*) 'internal_stability: '
	  write(6,*) aindex,dindex,tindex
	  call output_errout_stability(dt,sauxe1,sauxe2)
	end if

	rindex = tindex

	deallocate(sauxe1,sauxe2)

	!write(6,*) 'rindex = ',rindex,aindex,dindex

        end

!*****************************************************************

	subroutine momentum_viscous_stability(ahpar,rindex,dstab)

! computes stability for viscosity
!
! stability is computed for dt == 1

	use geom
	use internal
	use diffusion
	use evgeom
	use levels
	use basin
        use shympi

	implicit none

        include 'param.h'

	double precision ahpar
	double precision rindex
        !double precision dstab(nlvdi,nel)
        double precision dstab(nlvdi,nel)

	integer ie,ii,iei,l,lmax,k
	double precision u,v,ui,vi
	double precision anu,ax,ay
	!double precision area,areai
        double precision area,areai
	double precision dt
	double precision afact,r
        double precision a,ai,amax

	rindex = 0.
	if( ahpar .le. 0 ) return

	amax = 0.

	do ie=1,nel

	  lmax = ilhv(ie)
	  area = 12. * ev(10,ie)

	  do l=1,lmax

	    a = 0.
	    r = 0.
	    do ii=1,3

              iei = auxv_iei(ii,ie)
	      k = nen3v(ii,ie)
              if( iei .le. 0 ) iei = ie

              areai = 12. * ev(10,iei)

	      r = r + rdistv(k)
	      anu = ahpar * difhv(l,ie)
              ai = 2. * anu / ( area + areai )
              a = a + ai
	    end do
	    r = r/3.

	    a = a * r
	    amax = max(amax,a)
	    dstab(l,ie) = a

          end do

	end do

	!rindex = amax
	rindex = shympi_max(amax)

	end

!******************************************************************

	subroutine momentum_advective_stability(rlin,rindex,astab)

! computes courant number of advective terms in momentum equation
!
! stability is computed for dt == 1

	use internal
	use geom_dynamic
	use layer_thickness
	use hydro_admin
	use evgeom
	use levels
	use basin
        use shympi

	implicit none

        include 'param.h'

	double precision rlin		   !factor for advective terms - normally 1
	double precision rindex		   !stability index (return)
	!double precision astab(nlvdi,nel)      !stability matrix (return)
        double precision astab(nlvdi,nel)      !stability matrix (return)

	integer ie,l,ii,k,lmax,iweg
	!double precision cc,cmax
        double precision cc,cmax
	double precision ut,vt
	double precision h
        double precision area,vol
	double precision r
        double precision b,c,f,ftot

	cmax = 0.
	!call compute_stability_stats(-1,cc)

	do ie=1,nel
	  area = 12. * ev(10,ie)
	  lmax = ilhv(ie)
	  iweg = iwegv(ie)
	  do l=1,lmax

            h = hdenv(l,ie)
	    vol = area * h

  	    ut = utlnv(l,ie)
  	    vt = vtlnv(l,ie)

	    ftot = 0.
	    r = 0.
            do ii=1,3
                k = nen3v(ii,ie)
                b = ev(3+ii,ie)
                c = ev(6+ii,ie)
		r = r + rdistv(k)
                f = ut * b + vt * c
                if( f .lt. 0. ) ftot = ftot - f
            end do
	    r = r/3.

	    cc = rlin*r*area*ftot/vol
	    if( iweg .gt. 0 ) cc = 0.	! dry element
	    astab(l,ie) = cc
	    cmax = max(cmax,cc)
	    !call compute_stability_stats(0,cc)

	  end do
	end do

	!rindex = cmax
        rindex = shympi_max(cmax)
	!call compute_stability_stats(1,cc)

	end

!******************************************************************

        subroutine output_errout_stability(dt,sauxe1,sauxe2)

! outputs stability index for hydro timestep (internal) (error handling)

	use levels
	use basin, only : nkn,nel,ngr,mbw
        use conz_util
        use check
        use transforms

        implicit none

	include 'param.h'

        double precision dt
	double precision sauxe1(nlvdi,nel)
	double precision sauxe2(nlvdi,nel)
	double precision sauxn(nlvdi,nkn)

	logical bnos
	integer ie,l,lmax
	integer ia,id,it
	double precision aindex,dindex,tindex

! set ifnos in order to have output to nos file

	integer icall,iustab,ifnos
	save icall,iustab,ifnos
	data icall,iustab,ifnos /0,0,0/

	icall = icall + 1

	ia = 0
	id = 0
	it = 0
	aindex = 0.
	dindex = 0.
	tindex = 0.

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    if( sauxe1(l,ie) .gt. aindex ) then
	      aindex = sauxe1(l,ie)
	      ia = ie
	    end if
	    if( sauxe2(l,ie) .gt. dindex ) then
	      dindex = sauxe2(l,ie)
	      id = ie
	    end if
	    if( sauxe1(l,ie)+sauxe2(l,ie) .gt. tindex ) then
	      tindex = sauxe1(l,ie)+sauxe2(l,ie)
	      it = ie
	    end if
	  end do
	end do

	write(6,*) 'errout_stability: int-node stab-index stab-index*dt'
	write(6,*) 'advective: ',ia,aindex,aindex*dt
	write(6,*) 'diffusive: ',id,dindex,dindex*dt
	write(6,*) 'total:     ',it,tindex,tindex*dt

	call check_set_unit(6)
	if( ia .ne. 0 ) then
	  call check_elem(ia)
	end if
	if( id .ne. 0 .and. id .ne. ia ) then
	  call check_elem(id)
	end if
	if( it .ne. 0 .and. it .ne. id .and. it .ne. ia ) then
	  call check_elem(it)
	end if

	if( ifnos .gt. 0 .and. mod(icall,ifnos) .eq. 0 ) then
	  call e2n3d_minmax(+1,nlvdi,sauxe1,sauxn)
	  call conwrite(iustab,'.sta',1,778,nlvdi,sauxn)
	  call e2n3d_minmax(+1,nlvdi,sauxe2,sauxn)
	  call conwrite(iustab,'.sta',1,778,nlvdi,sauxn)
	  do ie=1,nel
	    do l=1,nlv
	      sauxe1(l,ie) = sauxe1(l,ie) + sauxe2(l,ie)
	    end do
	  end do
	  call e2n3d_minmax(+1,nlvdi,sauxe1,sauxn)
	  call conwrite(iustab,'.sta',1,778,nlvdi,sauxn)
	end if

	end

!*****************************************************************

        subroutine output_stability_node(dt,cwrite)

! outputs stability index for hydro timestep (internal)

	use levels
	use basin
        use para
        use output
        use conz_util

        implicit none

	include 'param.h'

        double precision dt
	double precision cwrite(nlvdi,nkn)

	double precision, save, allocatable :: smax(:)

	logical bnos
	integer ie,ii,k,l,lmax
	integer ia,id
	double precision sindex,smin

	include 'femtime.h'

	integer icall,iustab,ia_out(4)
	save icall,iustab,ia_out
	data icall,iustab /0,0/

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  call init_output('itmsti','idtsti',ia_out)
	  call increase_output(ia_out)
	  if( .not. has_output(ia_out) ) icall = -1
	  if( icall .lt. 0 ) return
	  call open_scalar_file(ia_out,1,1,'.stb')
	  allocate(smax(nkn))
	  smax = 0.
	end if

	icall = icall + 1

	if( .not. is_over_output(ia_out) ) return 

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    sindex = cwrite(l,k)
	    if( sindex .gt. smax(k) ) smax(k) = sindex
	  end do
	end do

	if( next_output(ia_out) ) then
	  do k=1,nkn				!convert to time step
	    if( smax(k) > 0 ) then
	      smax(k) = 1./smax(k)
	      if( smax(k) > dt_orig ) smax(k) = dt_orig
	    else
	      smax(k) = dt_orig
	    end if 
	  end do
	  call write_scalar_file(ia_out,778,1,smax)
	  smax = 0.
	end if

	end

!*****************************************************************

        subroutine output_stability(dt,sauxe1,sauxe2)

! outputs stability index for hydro timestep (internal)

	use levels
	use basin
        use para
        use output
        use conz_util

        implicit none

	include 'param.h'

        double precision dt
	double precision sauxe1(nlvdi,nel)	!advective stability index
	double precision sauxe2(nlvdi,nel)	!diffusive stability index

	double precision, save, allocatable :: smax(:)

	logical bnos
	integer ie,ii,k,l,lmax
	integer ia,id
	double precision sindex,smin
	double precision sx,sn

	include 'femtime.h'

	integer icall,ia_out(4)
	save icall,ia_out
	data icall /0/

!	idtsti = 3600
!	idtsti = 0
!	itmsti = -1

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  if( icall .lt. 0 ) return
	  call init_output('itmsti','idtsti',ia_out)
	  call increase_output(ia_out)
	  if( .not. has_output(ia_out) ) icall = -1
	  if( icall .lt. 0 ) return
	  call open_scalar_file(ia_out,1,1,'sti')
	  allocate(smax(nkn))
	  smax = 0.
	end if

	icall = icall + 1

	!write(111,*) 'sti: ',it,t_act,ia_out

	if( .not. is_over_output(ia_out) ) return 

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    sindex = sauxe1(l,ie)+sauxe2(l,ie)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( sindex .gt. smax(k) ) smax(k) = sindex
	    end do
	  end do
	end do

	if( next_output(ia_out) ) then
	  do k=1,nkn				!convert to time step
	    if( smax(k) > 0 ) then
	      smax(k) = 1./smax(k)
	      if( smax(k) > dt_orig ) smax(k) = dt_orig
	    else
	      smax(k) = dt_orig
	    end if 
	  end do
	  call write_scalar_file(ia_out,779,1,smax)
	  smax = 0.
	end if

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine parallel_test

! tests parallel implementation

	use levels, only : nlvdi,nlv
	use basin

	implicit none

	include 'param.h'

	double precision dt,rkpar,azpar,rindex
	double precision robs,wsink
	double precision wsinkv(0:nlvdi,nkn)
	double precision saux(nlvdi,nkn)

	azpar = 0.
	rkpar = 0.
	robs = 0.
	wsink = 0.

	write(6,*) 'parallel test...'
	call compute_stability(robs,wsink,wsinkv,rkpar,azpar,rindex,saux)
	write(6,*) 'parallel is ok.'

	end

!*****************************************************************

!        subroutine conzstab(cn1,co1
!     +			,ddt
        subroutine conzstab(ddt,robs,wsink,wsinkv,rkpar,difhv,difv,difmol,azpar &
     &			,adpar,aapar,sindex,istot,isact,nlvddi,nlev)
!
! checks stability
!
! cn     new concentration
! co     old concentration
! caux   aux vector
! clow	 lower diagonal of vertical system
! chig	 upper diagonal of vertical system
! ddt    time step
! robs	 factor for nudging
! rkpar  horizontal turbulent diffusivity
! difhv  horizontal turbulent diffusivity (variable between elements)
! difv   vertical turbulent diffusivity
! difmol vertical molecular diffusivity
! azpar  time weighting parameter
! adpar  time weighting parameter for vertical diffusion (ad)
! aapar  time weighting parameter for vertical advection (aa)
! sindex stability index
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
	use depth
	use layer_thickness
	use diff_aux
	use bnd_dynamic
	use area
	use ts
	use hydro_vel
	use hydro_admin
	use evgeom
	use levels
	use basin
	use shympi
        use para
        use bndo_admin
        use openmp_admin
        use time_util
        use timing
#ifdef DEBUGON
        use mpi_io_admin
#endif

	implicit none
!
! arguments
	integer nlvddi,nlev
        !double precision cn1(nlvddi,nkn),co1(nlvddi,nkn)		!DPGGU
	double precision difmol
        double precision ddt,rkpar,azpar,adpar,aapar			!$$azpar
	double precision robs,wsink
	integer istot,isact
#ifdef DEBUGON
        integer e
        double precision difhv(nlvddi,nel_local)
	double precision wsinkv(0:nlvddi,nkn_local)
        double precision difv(0:nlvddi,nkn_local)
#else
        double precision difhv(nlvddi,nel)
	double precision wsinkv(0:nlvdi,nkn)
        double precision difv(0:nlvddi,nkn)
#endif
! common
	include 'femtime.h'
	include 'mkonst.h'

! local
	logical bdebug,bdebug1,debug
	integer k,ie,ii,l,iii
	integer lstart
	integer ilevel
	integer itot,isum	!$$flux
	logical berror
	integer kn(3)
        double precision sindex,rstol
	double precision us,vs
	double precision az,azt
	double precision aa,aat,ad,adt
	double precision aj,rk3,rv,aj4
	double precision hmed,hmbot,hmtop
	double precision rvptop,rvpbot
	double precision dt,w,aux
	double precision aux1,aux2,aux3,aux4,aux5
	double precision b(3),c(3),f(3)
	double precision rso,rsn,rsot,rsnt,rstot
	double precision hn,ho
        double precision wdiff(3)

	double precision difabs,difrel,volold,volnew,flxin,flxtot,diff
	double precision stabind,stabadv,stabdiff,stabvert,stabpoint
	double precision voltot

!------------------------------------------------------------
! big arrays
!------------------------------------------------------------
	!double precision cn(nlvddi,nkn)		!DPGGU	!FIXME
	!double precision co(nlvddi,nkn)
	!double precision cdiag(nlvddi,nkn)
	!double precision clow(nlvddi,nkn)
	!double precision chigh(nlvddi,nkn)
        !double precision cwrite(nlvddi,nkn)
        !double precision saux(nlvddi,nkn)
	double precision, allocatable :: cn(:,:)
	double precision, allocatable :: co(:,:)
	double precision, allocatable :: cdiag(:,:)
	double precision, allocatable :: clow(:,:)
	double precision, allocatable :: chigh(:,:)
        double precision, allocatable :: cwrite(:,:)
        double precision, allocatable :: saux(:,:)
!------------------------------------------------------------
! end of big arrays
!------------------------------------------------------------

	double precision cexpl
	double precision fw(3),fd(3)
	double precision fl(3)
! local (new)
	double precision clc(nlvddi,3), clm(nlvddi,3), clp(nlvddi,3)
	!double precision cl(0:nlvddi+1,3)
	double precision wl(0:nlvddi+1,3)
!
	double precision hdv(0:nlvddi+1)
	double precision haver(0:nlvddi+1)
	double precision hnew(0:nlvddi+1,3)
	double precision hold(0:nlvddi+1,3)
	double precision present(0:nlvddi+1)

        integer kstab
	double precision dtorig,time1

        !integer iustab
        !save iustab
        !data iustab /0/
! functions

	!write(6,*) 'conzstab called...'

        if(nlv.ne.nlev) stop 'error stop conzstab: level'

!-----------------------------------------------------------------
! allocation
!-----------------------------------------------------------------

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
        call shympi_exchange_halo_3d_nodes(mfluxv)
        call shympi_exchange_halo_3d_nodes(rtauv)

	allocate(cn(nlvddi,nkn_local),co(nlvddi,nkn_local),cdiag(nlvddi,nkn_local))
	allocate(clow(nlvddi,nkn_local),chigh(nlvddi,nkn_local))
	allocate(cwrite(nlvddi,nkn_local),saux(nlvddi,nkn_local))
#else
	allocate(cn(nlvddi,nkn),co(nlvddi,nkn),cdiag(nlvddi,nkn))
	allocate(clow(nlvddi,nkn),chigh(nlvddi,nkn))
	allocate(cwrite(nlvddi,nkn),saux(nlvddi,nkn))
#endif

!-----------------------------------------------------------------
! initialization
!-----------------------------------------------------------------

        bdebug1 = .true.
        bdebug1 = .false.
        debug = .true.
        debug = .false.
	bdebug=.false.
	berror=.false.

        if( bdebug1 ) then
                write(6,*) 'debug parameters in conz3d'
		write(6,*) ddt,rkpar,difmol,azpar,adpar,aapar
                write(6,*) istot,isact,nlvddi,nlv
                write(6,*) nkn,nel
        end if

	az=azpar		!$$azpar
	azt=1.-az
	ad=adpar
	adt=1.-ad
	aa=aapar
	aat=1.-aa

	if( aa .ne. 0. .and. nlv .gt. 1 ) then
	  write(6,*) 'aapar = ',aapar
	  write(6,*) 'Cannot use implicit vertical advection.'
	  write(6,*) 'This might be resolved in a future version.'
	  write(6,*) 'Please set aapar = 0 in the STR file.'
	  stop 'error stop conzstab: implicit vertical advection'
	end if

!	-----------------------------------------------------------------
!	 fractional time step
!	-----------------------------------------------------------------

	rstot=istot		!$$istot
	rso=(isact-1)/rstot
	rsn=(isact)/rstot
	rsot=1.-rso
	rsnt=1.-rsn

	dt=ddt/rstot

!	-----------------------------------------------------------------
!	 initialize global arrays for accumulation of implicit terms
!	-----------------------------------------------------------------

#ifdef DEBUGON
	do k=1,nkn_local
#else
	do k=1,nkn
#endif
          do l=1,nlv
	    !co(l,k)=cn1(l,k)	!DPGGU	!not used for stability
            cn(l,k)=0.          !Malta
            co(l,k)=0.
	    if( mfluxv(l,k) .gt. 0. ) co(l,k) = mfluxv(l,k)	!point sources
            cdiag(l,k)=0.
            clow(l,k)=0.
            chigh(l,k)=0.
            cwrite(l,k)=0.
          end do
	end do

!	-----------------------------------------------------------------
!	these are aux arrays (bigger than needed) to avoid checking for
!	what layer we are in -> we never get out of bounds
!	-----------------------------------------------------------------

        do l=0,nlv+1
	  hdv(l) = 0.		!layer thickness
          haver(l) = 0.
	  present(l) = 0.	!1. if layer is present
	  do ii=1,3
	    hnew(l,ii) = 0.	!as hreal but with zeta_new
	    hold(l,ii) = 0.	!as hreal but with zeta_old
	    !cl(l,ii) = 0.	!concentration in layer
	    wl(l,ii) = 0.	!vertical velocity
	  end do
	end do

!	-----------------------------------------------------------------
!	these are the local arrays for accumulation of implicit terms
!	(maybe we do not need them, but just to be sure...)
!	after accumulation we copy them on the global arrays
!	-----------------------------------------------------------------

        do l=1,nlv
	  do ii=1,3
	    clc(l,ii) = 0.
	    clm(l,ii) = 0.
	    clp(l,ii) = 0.
	  end do
	end do

!	-----------------------------------------------------------------
!	vertical velocities
!	-----------------------------------------------------------------

!-----------------------------------------------------------------
! loop over elements
!-----------------------------------------------------------------

#ifdef DEBUGON
        do e=1,nel_local
            if(bmpi) then
               ie=domain%elems%mapID(e)
            else
               ie=e
            end if
#else
        do ie=1,nel
#endif
	do ii=1,3
          k=nen3v(ii,ie)
	  kn(ii)=k
	  b(ii)=ev(ii+3,ie)
	  c(ii)=ev(ii+6,ie)
	end do

	aj=ev(10,ie)    !area of triangle / 12
	aj4=4.*aj
        ilevel=ilhv(ie)

! set up vectors for use in assembling contributions

        do l=1,ilevel
	  hdv(l) = hdeov(l,ie)		!use old time step -> FIXME
          haver(l) = 0.5 * ( hdeov(l,ie) + hdenv(l,ie) )
	  present(l) = 1.
	  do ii=1,3
	    k=kn(ii)
	    hn = hdknv(l,k)		! there are never more layers in ie
	    ho = hdkov(l,k)		! ... than in k
	    hold(l,ii) = rso * hn + rsot * ho
	    hnew(l,ii) = rsn * hn + rsnt * ho
	    !cl(l,ii) = co(l,k)
	    wl(l,ii) = wlnv(l,k) - wsink * wsinkv(l,k)
	  end do
	end do

	do l=ilevel+1,nlv
	  present(l) = 0.
	end do

!	set vertical velocities in surface and bottom layer
!
!	we do not set wl(0,ii) because otherwise we loose concentration
!	through surface
!
!	we set wl(ilevel,ii) to 0 because we are on the bottom
!	and there should be no contribution from this element
!	to the vertical velocity

	do ii=1,3
	  wl(ilevel,ii) = 0.
	end do

!-----------------------------------------------------------------
! loop over levels
!-----------------------------------------------------------------

        do l=1,ilevel

        us=az*utlnv(l,ie)+azt*utlov(l,ie)             !$$azpar
        vs=az*vtlnv(l,ie)+azt*vtlov(l,ie)

        rk3 = 3. * rkpar * difhv(l,ie)

	itot=0
	isum=0
	do ii=1,3
	  k=kn(ii)
	  f(ii)=us*b(ii)+vs*c(ii)	!$$azpar
	  if(f(ii).lt.0.) then	!flux out of node
	    itot=itot+1
	    isum=isum+ii
	  end if

! new weights for diffusion

          wdiff(ii) = wdifhv(ii,ii,ie)

!	  initialization to be sure we are in a clean state

	  fw(ii) = 0.
	  clc(l,ii) = 0.
	  clm(l,ii) = 0.
	  clp(l,ii) = 0.

!	  contributions from vertical advection
!
!	  in fw(ii) is explicit contribution
!	  the sign is for the term on the left side, therefore
!	  fw(ii) must be subtracted from the right side
!
!	  if we are in last layer, w(l,ii) is zero
!	  if we are in first layer, w(l-1,ii) is zero (see above)

	  w = wl(l-1,ii)                !top of layer
	  if( w .gt. 0. ) then          !out
	    fw(ii) = aat*w
	    clc(l,ii) = clc(l,ii) + aa*w
	  else
	    fw(ii) = 0.
	    clm(l,ii) = clm(l,ii) + aa*w
	  end if

	  w = wl(l,ii)                  !bottom of layer
	  if( w .gt. 0. ) then
	    clp(l,ii) = clp(l,ii) - aa*w
	  else
	    fw(ii) = fw(ii) - aat*w
	    clc(l,ii) = clc(l,ii) - aa*w
	  end if

!	  contributions from vertical diffusion
!
!	  in fd(ii) is explicit contribution
!	  the sign is for the term on the left side, therefore
!	  fd(ii) must be subtracted from the right side
!
!	  maybe we should use double precision layer thickness, or even the
!	  time dependent layer thickness

	  rvptop = difv(l-1,k) + difmol
	  rvpbot = difv(l,k) + difmol
	  hmtop = 2. * rvptop * present(l-1) / (hdv(l-1)+hdv(l))
	  hmbot = 2. * rvpbot * present(l+1) / (hdv(l)+hdv(l+1))

          fd(ii) = adt * ( hmtop + hmbot )

	  clc(l,ii) = clc(l,ii) + ad * ( hmtop + hmbot )
	  clm(l,ii) = clm(l,ii) - ad * ( hmtop )
	  clp(l,ii) = clp(l,ii) - ad * ( hmbot )
	end do

! sum explicit contributions

	do ii=1,3
	  k=kn(ii)
          hmed = hold(l,ii)                      !new ggu   !HACK
          cexpl = dt * aj4 * rk3 * hmed * wdiff(ii)	!bug fix 12.2.2010
	  clow(l,k) = clow(l,k) + cexpl
          cexpl = dt * aj4 * 3. * f(ii)
          if( cexpl .lt. 0. ) then             !flux out of node
	    !chigh(l,k) = chigh(l,k) - cexpl
          end if
          if( cexpl .gt. 0. ) then             !flux into node
	    chigh(l,k) = chigh(l,k) + cexpl
          end if
          cn(l,k) = cn(l,k) + dt * aj4 * ( fw(ii) + fd(ii) )
          co(l,k) = co(l,k) + dt * aj4 * hmed * robs * rtauv(l,k) !nudging
	end do

	end do		! loop over l

! set up implicit contributions
!
! cdiag is diagonal of tri-diagonal system
! chigh is high (right) part of tri-diagonal system
! clow is low (left) part of tri-diagonal system

	do ii=1,3
	  clm(1,ii) = 0.
	  clp(ilevel,ii) = 0.
	end do

        do l=1,ilevel
	  do ii=1,3
	    k=kn(ii)
	    !clow(l,k)  = clow(l,k)  + aj4 * dt * clm(l,ii)
	    !chigh(l,k) = chigh(l,k) + aj4 * dt * clp(l,ii)
	    !cdiag(l,k) = cdiag(l,k) + aj4 * dt * clc(l,ii)
	    !clow(l,k)  = clow(l,k)  + aj4 * hold(l,ii)
            hmed = min(hold(l,ii),hnew(l,ii))
	    cdiag(l,k) = cdiag(l,k) + aj4 * hmed
	  end do
	end do

	end do		! loop over ie



#ifndef DEBUGON
        if(shympi_partition_on_elements()) then
          !call shympi_comment('shympi_elem: exchange and compute stabind')
          if(ln_timing) time1 = shympi_wtime()
          call shympi_exchange_and_sum_3D_nodes(cn)
          call shympi_exchange_and_sum_3D_nodes(co)
          call shympi_exchange_and_sum_3D_nodes(cdiag)
          call shympi_exchange_and_sum_3D_nodes(clow)
          call shympi_exchange_and_sum_3D_nodes(chigh)
          if(ln_timing) comm_scalar_time = comm_scalar_time + shympi_wtime() - time1
        end if
#endif

!-----------------------------------------------------------------
! compute stability
!
! cdiag		volume of cell
! chigh		flux due to horizontal advection
! clow		flux due to horizontal diffusion
! cn		flux due to vertical advection and diffusion (explicit)
! co		flux due to point sources and nudging
!-----------------------------------------------------------------

        stabind = 0.		!total max stability index
        stabadv = 0.		!advective max stability index
        stabdiff = 0.		!diffusive max stability index
        stabvert = 0.		!vertical max stability index
        stabpoint = 0.		!point source max stability index
        kstab = 0		!node with highest stabind

	do k=1,nkn
	  !bdebug1 = k .eq. 1402
	  !bdebug1 = k .eq. 1405
	  bdebug1 = k .eq. -1
	  ilevel = ilhkv(k)
          if( .not. is_zeta_bound(k) ) then	!FIXME
          !if( is_inner(k) ) then	!FIXME
          !if( .true. ) then	!FIXME
	   do l=1,ilevel
            voltot = cdiag(l,k)
            flxtot = chigh(l,k) + clow(l,k) + cn(l,k) + co(l,k)
	    if( bdebug1 ) write(99,*) k,l,voltot,flxtot
            if( voltot .gt. 0. ) then
                  aux1 = flxtot / voltot
                  if( aux1 .gt. stabind ) kstab = k
                  stabind = max(stabind,aux1)
		  cwrite(l,k) = aux1		!save for write
                  aux2 = chigh(l,k) / voltot
                  stabadv = max(stabadv,aux2)
		  saux(l,k) = aux2		!for adv. stab.
                  aux3 = clow(l,k) / voltot
                  stabdiff = max(stabdiff,aux3)
                  aux4 = cn(l,k) / voltot
                  stabvert = max(stabvert,aux4)
                  aux5 = co(l,k) / voltot
                  stabpoint = max(stabpoint,aux5)
	          if( bdebug1 ) write(99,*) aux1,aux2,aux3,aux4,aux5

!		  aux=flxtot / voltot
!		  if( 300*aux/dt .gt. 1000 ) then
!			  write(6,*) is_boundary(k)
!			  write(6,*) is_external_boundary(k)
!			  write(6,*) is_internal_boundary(k)
!			  write(6,*) is_inner(k)
!			  write(6,*) stabind,stabadv,stabdiff,stabvert
!			  call check_set_unit(6)
!			  call check_node(k)
!		  end if

            else
		  cwrite(l,k) = 0
		  saux(l,k) = 0.
            end if
	   end do
          else
	   do l=1,ilevel
		  cwrite(l,k) = 0
		  saux(l,k) = 0.
           end do
          end if
	end do

        if(ln_timing) time1 = shympi_wtime()
        stabind = shympi_max(stabind)
        if(ln_timing) comm_scalar_time = comm_scalar_time + shympi_wtime() - time1


!        write(6,*) 'stab check: ',nkn,nlv
!        call check2Dr(nlvddi,nlv,nkn,cwrite,0.,0.,"NaN check","cstab")

!-----------------------------------------------------------------
! in stabind is stability index (advection and diffusion)
! in cdiag is the local value of the stability index
! in cwrite is the value of the stability index for each node
!
! istot  is saved and returned from subroutine (number of iterations)
! sindex is saved and returned from subroutine (stability index)
!-----------------------------------------------------------------

	rstol = getpar('rstol')
        istot = 1 + stabind / rstol
        sindex = stabind

	call get_orig_timestep(dtorig)
	if( .not. openmp_in_parallel() ) then
	  call output_stability_node(dtorig,cwrite)
	end if

!        if( .false. ) then
!        !if( idt .le. 3 ) then
!	  write(6,*) 'kstab = ',kstab,'  stabind = ',stabind
!          call conwrite(iustab,'.stb',1,777,nlvddi,cwrite)
!        end if

!        call stb_histo(it,nlvddi,nkn,ilhkv,cwrite)

!-----------------------------------------------------------------
! allocation
!-----------------------------------------------------------------

	deallocate(cn,co,cdiag)
	deallocate(clow,chigh)
	deallocate(cwrite,saux)

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
	end


!*****************************************************************

!--------------------------------------------------------------------------------
        end module stability
!--------------------------------------------------------------------------------
