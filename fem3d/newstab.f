
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010-2012,2014-2019  Georg Umgiesser
!    Copyright (C) 2012,2016  Christian Ferrarin
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

c routines for stability computations
c
c revision log :
c
c 19.02.2010	ggu	new file to contain stability computations
c 26.02.2010	ggu	internal_stability restructured (on element)
c 08.03.2010	ggu	run only down to avail layers in stability (bug fix)
c 23.03.2010	ggu	changed v6.1.1
c 26.01.2011	ggu	robs for nudging implemented
c 16.02.2011	ggu	pass robs to subroutines, write stb-ind to nos file
c 20.05.2011	ggu	allow for elimination of elems due to high rstab
c 31.05.2011	ggu	changed VERS_6_1_23
c 01.06.2011	ggu	wsink for stability integrated
c 12.07.2011	ggu	new routine output_stability()
c 14.07.2011	ggu	new routine output_stability_node()
c 21.06.2012	ggu&ccf	variable vertical sinking velocity integrated
c 07.03.2014	ggu	changed VERS_6_1_72
c 08.04.2014	ggu	use rlin to determine advective stability
c 05.05.2014	ggu	changed VERS_6_1_74
c 26.11.2014	ggu	changed VERS_7_0_7
c 23.12.2014	ggu	changed VERS_7_0_11
c 09.01.2015	ggu	changed VERS_7_0_12
c 19.01.2015	ggu	changed VERS_7_1_3
c 20.05.2015	ggu	always compute stability, call to conzstab changed
c 05.06.2015	ggu	changed VERS_7_1_12
c 13.07.2015	ggu	changed VERS_7_1_51
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 23.09.2015	ggu	changed VERS_7_2_4
c 20.10.2015	ggu	in output_stability() bug that icall was not adjouned
c 23.10.2015	ggu	changed VERS_7_3_9
c 09.11.2015	ggu	changed VERS_7_3_13
c 25.05.2016	ggu	changed VERS_7_5_10
c 20.10.2016	ccf	pass rtauv for differential nudging
c 12.01.2017	ggu	changed VERS_7_5_21
c 05.12.2017	ggu	changed VERS_7_5_39
c 13.04.2018	ggu	re-structured, included gravity wave stability
c 16.04.2018	ggu	use also ilin to compute stability
c 16.02.2019	ggu	changed VERS_7_5_60
c 06.11.2019	ggu	femtime eliminated
c 30.03.2021	ggu	better error output
c 20.03.2022	ggu	upgraded to da_out
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
c
c it is assumed that the program will exit after this call

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
	double precision dtime
	real saux(nlvdi,nkn)

	integer ia,id
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
	id = 0
	call get_act_dtime(dtime)
	call shy_write_scalar(id,'sta',dtime,1,778,nlvdi,saux)
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
	  write(6,*) 'aindex,dindex,gindex,tindex: '
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
	integer iemax,iee
	real tindex
	real sauxn(nlvdi,nkn)

c set ifnos in order to have output to nos file

	integer icall,iustab,ifnos
	save icall,iustab,ifnos
	data icall,iustab,ifnos /0,0,0/

	integer ieext

	icall = icall + 1

	iemax = 0
	tindex = 0.

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    if( sauxe(l,ie) .gt. tindex ) then
	      tindex = sauxe(l,ie)
	      iemax = ie
	    end if
	  end do
	end do

	ie = iemax
	iee = ieext(ie)
	write(6,*) 'errout_stability: '
	write(6,*) '  ie-intern   ie-extern' //
     +			'   stab-index    dt*stab-index'
	write(6,*) ie,iee,tindex,tindex*dt

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
	integer ia,id,idc,nvar
	real sindex,smin
	real dtorig
	double precision dtime
	logical next_output_d,has_output_d,is_over_output_d

	integer, save :: icall = 0
	double precision, save :: da_out(4)

	real getpar

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  call init_output_d('itmsti','idtsti',da_out)
	  call increase_output_d(da_out)
	  if( .not. has_output_d(da_out) ) icall = -1
	  if( icall .lt. 0 ) return
	  nvar = 1
          call shyfem_init_scalar_file('stb',nvar,.true.,id)
          da_out(4) = id
	  allocate(smax(nkn))
	  smax = 0.
	end if

	icall = icall + 1

	if( .not. is_over_output_d(da_out) ) return 

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    sindex = cwrite(l,k)
	    if( sindex .gt. smax(k) ) smax(k) = sindex
	  end do
	end do

	if( next_output_d(da_out) ) then
	  call get_orig_timestep(dtorig)
	  do k=1,nkn				!convert to time step
	    if( smax(k) > 0 ) then
	      smax(k) = 1./smax(k)
	      if( smax(k) > dtorig ) smax(k) = dtorig
	    else
	      smax(k) = dtorig
	    end if 
	  end do
	  call get_act_dtime(dtime)
	  id = nint(da_out(4))
          idc = 778
          call shy_write_scalar_record2d(id,dtime,idc,smax)
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
	integer ia,id,idc,nvar
	real sindex,smin
	real sx,sn
	real dtorig
	double precision dtime
	logical next_output_d,has_output_d,is_over_output_d

	integer, save :: icall = 0
	double precision, save :: da_out(4)

	real getpar

c	idtsti = 3600
c	idtsti = 0
c	itmsti = -1

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  call init_output_d('itmsti','idtsti',da_out)
	  call increase_output_d(da_out)
	  if( .not. has_output_d(da_out) ) icall = -1
	  if( icall .lt. 0 ) return
	  nvar = 1
          call shyfem_init_scalar_file('sti',nvar,.true.,id)
          da_out(4) = id
	  allocate(smax(nkn))
	  smax = 0.
	end if

	icall = icall + 1

	!write(111,*) 'sti: ',it,t_act,ia_out

	if( .not. is_over_output_d(da_out) ) return 

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

	if( next_output_d(da_out) ) then
	  call get_orig_timestep(dtorig)
	  do k=1,nkn				!convert to time step
	    if( smax(k) > 0 ) then
	      smax(k) = 1./smax(k)
	      if( smax(k) > dtorig ) smax(k) = dtorig
	    else
	      smax(k) = dtorig
	    end if 
	  end do
	  call get_act_dtime(dtime)
	  id = nint(da_out(4))
          idc = 779
          call shy_write_scalar_record2d(id,dtime,idc,smax)
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

