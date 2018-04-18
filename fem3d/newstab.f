c
c $Id: newstab.f,v 1.3 2010-03-08 17:46:45 georg Exp $
c
c routines for stability computations
c
c revision log :
c
c 19.02.2010    ggu     new file to contain stability computations
c 26.02.2010    ggu     internal_stability restructured (on element)
c 08.03.2010    ggu     run only down to avail layers in stability (bug fix)
c 26.01.2011    ggu     robs for nudging implemented
c 16.02.2011    ggu     pass robs to subroutines, write stb-ind to nos file
c 20.05.2011    ggu     allow for elimination of elems due to high rstab
c 01.06.2011    ggu     wsink for stability integrated
c 12.07.2011    ggu     new routine output_stability()
c 14.07.2011    ggu     new routine output_stability_node()
c 21.06.2012    ggu&ccf variable vertical sinking velocity integrated
c 08.04.2014    ggu	use rlin to determine advective stability
c 20.05.2015    ggu	always compute stability, call to conzstab changed
c 20.10.2015    ggu	in output_stability() bug that icall was not adjouned
c 20.10.2016    ccf     pass rtauv for differential nudging
c 13.04.2018    ggu	re-structured, included gravity wave stability
c 16.04.2018    ggu	use also ilin to compute stability
c
c*****************************************************************
c*****************************************************************
c*****************************************************************
c
c typical usage :
c
c call scalar_stability before advection of similar variables
c
c notes :
c
c	scal3sh							newcon
c		scalar_info_stability				newstab
c		scalar_stability				newstab
c			scalar_compute_stability		newstab
c				conzstab			newcon
c
c	set_timestep						subtime
c		hydro_stability					newstab
c			hydro_internal_stability		newstab
c				momentum_advective_stability	newexpl
c				momentum_viscous_stability	newexpl
c				gravity_wave_stability		newstab
c
c*****************************************************************

!==================================================================
        module stab
!==================================================================

	implicit none

	integer, parameter :: ndim_stab = 50
	integer, save :: nentry = 0
	real, save :: rkind(2,ndim_stab) = 0.

!==================================================================
        end module stab
!==================================================================

	subroutine scalar_compute_stability(robs,rtauv,wsink,wsinkv
     +					,rkpar,azpar,rindex,saux)

c computes stability index

	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin

	implicit none

	real robs
	real rtauv(nlvdi,nkn)
	real wsink
	real wsinkv(0:nlvdi,nkn)
        real rkpar
        real azpar
        real rindex
        real saux(nlvdi,nkn)

        real adpar,aapar
        real difmol
	real ddt
        integer isact,istot

        real getpar

c----------------------------------------------------------------
c set parameters
c----------------------------------------------------------------

	adpar=getpar('adpar')
	aapar=getpar('aapar')

        isact = 1
        istot = 1
        difmol = 0.
	ddt = 1.		!always for 1 sec
	rindex = 0.

c----------------------------------------------------------------
c call conzstab
c----------------------------------------------------------------

        call conzstab(
     +          ddt,robs,rtauv,wsink,wsinkv,rkpar,difhv,difv
     +		,difmol,azpar,adpar,aapar
     +          ,rindex,istot,isact,nlvdi,nlv)

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*****************************************************************

        subroutine scalar_basic_stability(dt,rkpar,rindex)

c computes scalar stability without nudging and sinking

	use levels, only : nlvdi,nlv
	use basin

        implicit none

	real dt
        real rkpar
        real rindex

        integer istot
	real azpar
	real robs
	real wsink
	real, allocatable :: wsinkv(:,:)
	real, allocatable :: rtauv(:,:)
	real, allocatable :: saux(:,:)

	allocate(wsinkv(0:nlvdi,nkn))
	allocate(rtauv(nlvdi,nkn))
	allocate(saux(nlvdi,nkn))

	wsinkv = 0.
	rtauv = 0.
	saux = 0.
	robs = 0.
	wsink = 0.

	call getaz(azpar)
	call scalar_compute_stability(robs,rtauv,wsink,wsinkv,rkpar,azpar,
     +					rindex,saux)

	end

c*****************************************************************

        subroutine scalar_stability(dt,robs,rtauv,wsink,wsinkv
     +					,rkpar,rindex,istot,saux)

c gets stability index (if necessary computes it)

	use levels, only : nlvdi,nlv
	use basin

        implicit none

	real dt
	real robs
	real rtauv(nlvdi,nkn)
	real wsink
	real wsinkv(0:nlvdi,nkn)
        real rkpar
        real rindex
        integer istot
	real saux(nlvdi,nkn)

	real azpar

c----------------------------------------------------------------
c compute stability index
c----------------------------------------------------------------

	call getaz(azpar)
	call scalar_compute_stability(robs,rtauv,wsink,wsinkv,rkpar,azpar,
     +					rindex,saux)

c----------------------------------------------------------------
c scale to real time step dt
c----------------------------------------------------------------

	rindex = dt * rindex
	istot = 1 + rindex

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

        end

c*****************************************************************

        subroutine scalar_info_stability(dt,robs,rtauv,wsink,wsinkv
     +					,rkpar,rindex,istot,saux)

c gets stability index (if necessary computes it)

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        real dt
	real robs
	real rtauv(nlvdi,nkn)
	real wsink
	real wsinkv(0:nlvdi,nkn)
        real rkpar
        real rindex
        integer istot
	real saux(nlvdi,nkn)

	integer ia,iustab
	integer l,k
	real aindex
	real azpar

c----------------------------------------------------------------
c compute stability index
c----------------------------------------------------------------

	call getaz(azpar)
	call scalar_compute_stability(robs,rtauv,wsink,wsinkv,rkpar,azpar,
     +					rindex,saux)
	rindex = dt * rindex
	istot = 1 + rindex

c----------------------------------------------------------------
c write to terminal
c----------------------------------------------------------------

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

cggu protect
	write(6,*) 'scalar_info_stability:'
	write(6,*) rkpar,azpar,rindex,istot
	write(6,*) ia,aindex,dt*aindex
	iustab = 0
	call conwrite(iustab,'.sta',1,778,nlvdi,saux)
cggu protect

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine hydro_stability(dt,rindex)

c computes stability index for hydro timestep - no error output

        implicit none

        real dt
        real rindex

	call hydro_internal_stability(0,dt,rindex)

	end

c**********************************************************************

        subroutine error_stability(dt,rindex)

c computes stability index for hydro timestep - with error output
c
c after this call the program should abort

        implicit none

        real dt
        real rindex

	call hydro_internal_stability(1,dt,rindex)

	end

c**********************************************************************

        subroutine eliminate_stability(rmax)

c eliminates elements with stability index higher than rmax

        implicit none

        real rmax

	integer mode
        real rindex,dt

	mode = 2
	dt = 0.
	rindex = rmax
	call hydro_internal_stability(mode,dt,rindex)

	end

c**********************************************************************

        subroutine hydro_internal_stability(mode,dt,rindex)

c computes stability index for hydro timestep (internal)
c
c mode = 0		normal call, compute stability
c mode = 1		error call, compute stability and write error message
c mode = 2		eliminate elements with r>rindex

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

	integer mode		!0: normal call  1:error output
        real dt			!time step to be used
        real rindex		!stability index (return)

	integer ie,l,lmax,iweg,ilin,ibarcl
        real rkpar,azpar,ahpar,rlin
	real dindex,aindex,tindex,sindex,gindex
	real rmax

	logical openmp_in_parallel

	real, allocatable :: sauxe1(:,:)
	real, allocatable :: sauxe2(:,:)
	real, allocatable :: sauxe3(:)
	real, allocatable :: sauxe(:,:)

	real getpar
	logical is_i_nan

        rkpar = 0.
	azpar = 1.
	ahpar = getpar('ahpar')
	ibarcl = nint(getpar('ibarcl'))
	rlin = getpar('rlin')
	ilin = nint(getpar('ilin'))
	if( ibarcl == 1 ) ilin = 0	!for baroclinic use advection stability
	rlin = rlin * (1-ilin)

	allocate(sauxe1(nlvdi,nel),sauxe2(nlvdi,nel))
	allocate(sauxe(nlvdi,nel),sauxe3(nel))
	sauxe1 = 0.
	sauxe2 = 0.

	rmax = 1.e+30
	if( mode .eq. 2 ) rmax = rindex
	if( mode .eq. 2 ) then
		write(6,*) 'eliminating rmax: ',rmax
	end if

	call momentum_advective_stability(rlin,aindex,sauxe1)
	call momentum_viscous_stability(ahpar,dindex,sauxe2)
	call gravity_wave_stability(gindex,sauxe3)

	do ie=1,nel
	  sauxe(:,ie) = sauxe1(:,ie) + sauxe2(:,ie) + sauxe3(ie)
	end do

	if( .not. openmp_in_parallel() ) then
          call output_stability(dt,sauxe)	!in case write to file
	end if

	tindex = 0.
	do ie=1,nel
	  lmax = ilhv(ie)
	  iweg = 0
	  do l=1,lmax
	    sindex = sauxe(l,ie)
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

	tindex = tindex*dt
	aindex = aindex*dt
	dindex = dindex*dt
	gindex = gindex*dt

	if( mode .eq. 1 ) then		!error output
	  write(6,*) 'hydro_internal_stability: '
	  write(6,*) aindex,dindex,gindex,tindex
	  call output_errout_stability(dt,sauxe)
	end if

	rindex = tindex

	deallocate(sauxe1,sauxe2,sauxe3,sauxe)

	!write(6,*) 'rindex = ',rindex,aindex,dindex

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine output_errout_stability(dt,sauxe)

c outputs stability index for hydro timestep (internal) (error handling)

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        real dt
	real sauxe(nlvdi,nel)

	logical bnos
	integer ie,l,lmax
	integer ia
	real tindex
	real sauxn(nlvdi,nkn)

c set ifnos in order to have output to nos file

	integer icall,iustab,ifnos
	save icall,iustab,ifnos
	data icall,iustab,ifnos /0,0,0/

	icall = icall + 1

	ia = 0
	tindex = 0.

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    if( sauxe(l,ie) .gt. tindex ) then
	      tindex = sauxe(l,ie)
	      ia = ie
	    end if
	  end do
	end do

	write(6,*) 'errout_stability: int-node stab-index stab-index*dt'
	write(6,*) 'total:     ',ia,tindex,tindex*dt

	if( ifnos .gt. 0 .and. mod(icall,ifnos) .eq. 0 ) then
	  call e2n3d_minmax(+1,nlvdi,sauxe,sauxn)
	  call conwrite(iustab,'.sta',1,778,nlvdi,sauxn)
	end if

	end

c*****************************************************************

        subroutine output_stability_node(dt,cwrite)

c outputs stability index for hydro timestep (internal)

	use levels
	use basin

        implicit none

        real dt
	real cwrite(nlvdi,nkn)

	real, save, allocatable :: smax(:)

	logical bnos
	integer ie,ii,k,l,lmax
	integer ia,id
	real sindex,smin
	logical has_output,next_output,is_over_output

	include 'femtime.h'

	integer icall,iustab,ia_out(4)
	save icall,iustab,ia_out
	data icall,iustab /0,0/

	real getpar

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

c*****************************************************************

        subroutine output_stability(dt,sauxe)

c outputs stability index for hydro timestep (internal)

	use levels
	use basin

        implicit none

        real dt
	real sauxe(nlvdi,nel)	!stability index for element

	real, save, allocatable :: smax(:)

	logical bnos
	integer ie,ii,k,l,lmax
	integer ia,id
	real sindex,smin
	real sx,sn
	logical next_output,has_output,is_over_output

	include 'femtime.h'

	integer icall,ia_out(4)
	save icall,ia_out
	data icall /0/

	real getpar

c	idtsti = 3600
c	idtsti = 0
c	itmsti = -1

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
	    sindex = sauxe(l,ie)
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

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine gravity_wave_stability(gindex,garray)

	use basin
	use mod_hydro

	implicit none

	real gindex
	real garray(nel)

	include 'pkonst.h'

	integer ie,ii,ii1,k1,k2
	real distmin,d,dx,dy
	real am,az
	real hz,ri

	integer, save :: icall = 0
	real,save,allocatable :: dist(:)
        real, parameter :: high = 1.e+30

	real getpar

	gindex = 0.
	garray = 0.

	if( icall < 0 ) return

	if( icall == 0 ) then
	  az = getpar('azpar')
	  am = getpar('ampar')
	  if( az >= 0.5 .and. am >= 0.5 ) then
	    if( az == am ) then	!unconditionally stable
	      icall = -1
	    else
	      goto 99
	    end if
	  else if( az == 0. .and. am == 1. ) then
	    !ok
	  else if( az == 1. .and. am == 0. ) then
	    !ok
	  else
	    goto 99
	  end if
	  if( icall < 0 ) return
	  icall = 1

	  allocate(dist(nel))
          do ie=1,nel
           distmin = high
           do ii=1,3
            ii1 = mod(ii,3) + 1
            k1 = nen3v(ii,ie)
            k2 = nen3v(ii1,ie)
            call compute_distance(xgv(k1),ygv(k1),xgv(k2),ygv(k2),dx,dy)
            d = dx**2 + dy**2
            distmin = min(distmin,d)
	   end do
	   dist(ie) = sqrt(distmin)
          end do
	end if

	gindex = 0.
	do ie=1,nel
	  hz = maxval( hm3v(:,ie) + zenv(:,ie) )
	  hz = maxval( hm3v(:,ie) )
	  ri = sqrt(grav*hz) / dist(ie)
	  gindex = max(gindex,ri)
	  garray(ie) = ri
	end do

	!write(6,*) nel,gindex,1./gindex
	!stop

	return
   99	continue
	write(6,*) 'azpar,ampar: ',az,am
	write(6,*) 'this combination of parameters is not allowed'
	write(6,*) 'for explicit runs please use'
	write(6,*) '  either az=1 and am=0 or az=0 and am=1'
	write(6,*) 'for (semi-)implicit runs please use'
	write(6,*) '  az==am and az>=0.5'
	stop 'error stop gravity_wave_stability: az and am'
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine parallel_test

c tests parallel implementation

	use levels, only : nlvdi,nlv
	use basin

	implicit none

	real dt,rkpar,azpar,rindex
	real robs,wsink
	real rtauv(nlvdi,nkn)
	real wsinkv(0:nlvdi,nkn)
	real saux(nlvdi,nkn)

	azpar = 0.
	rkpar = 0.
	robs = 0.
	wsink = 0.
	rtauv = 0.

	write(6,*) 'parallel test...'
	call scalar_compute_stability(robs,rtauv,wsink,wsinkv,rkpar,azpar,
     +					rindex,saux)
	write(6,*) 'parallel is ok.'

	end

c*****************************************************************

