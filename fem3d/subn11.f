
!--------------------------------------------------------------------------
!
!    Copyright (C) 1992-1994,1996-1999,2001,2003-2015  Georg Umgiesser
!    Copyright (C) 2017-2019  Georg Umgiesser
!    Copyright (C) 2008,2014  Christian Ferrarin
!    Copyright (C) 2010  Debora Bellafiore
!    Copyright (C) 2018-2019  Marco Bajo
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

c boundary condition routines
c
c contents :
c
c subroutine sp111(mode)		set up boundary/initial cond.
c
c subroutine z_tilt			tilting of boundary surface (fixed)
c subroutine c_tilt			tilting of boundary surface (Coriolis)
c
c subroutine init_tilt(ibc)       	finds tilting node in node list
c subroutine init_flux(ibc)		initializes flux boundary
c
c subroutine set_mass_flux		sets up (water) mass flux array mfluxv
c subroutine adjust_mass_flux		adjusts mass flux for dry nodes
c subroutine make_scal_flux(what,r3v,scal,sflux,sconz,ssurf)	sets scal flux
c subroutine flux_debug(what,mfluxv,sflux,sconz)
c subroutine check_scal_flux(what,scal,sconz)			checks scal flux
c 
c subroutine init_scal_bc(r3v)		initializes array for scalar BC
c subroutine mult_scal_bc(r3v,value)	multiplies array for scalar BC by value
c 
c subroutine dist_3d(nlvddi,r3v,kn,nbdim,values)
c subroutine dist_horizontal(nlvddi,r3v,n,value)
c subroutine aver_horizontal(nlvddi,r3v,n,value)
c
c subroutine print_scal_bc(r3v)		prints non-flag entries of scalar BC
c subroutine get_bflux(k,flux)		returns boundary flux of node k
c 
c subroutine level_flux(dtime,levflx,kn,rw) compute discharge from water level
c subroutine z_smooth(z)		smooths z values
c subroutine flow_out_piece_new(z,rout)
c
c function get_discharge(ibc)		returns discharge through boundary ibc
c
c revision log :
c
c 05.08.1992	ggu	$$ibtyp3 - implementation of ibtyp=3
c 31.08.1992	ggu	$$rqv1   - rqv initialized in sp159b, not here
c 05.09.1992	ggu	$$close1 - rqv initialized here
c 27.10.1993	ggu	$$roger  - ibtyp=70 (nwe-shelf)
c 05.11.1993	ggu	$$roger  - call to roger has been commented
c 11.01.1994	ggu	$$restart - restart is reading all variables
c 20.01.1994	ggu	$$zeov - initialize zeov
c 20.01.1994	ggu	$$conz - impl. of conz with bnd(12,.)
c 31.01.1994	ggu	$$conzbc - impl. of bc for conz with rcv
c 24.03.1994	ggu	$$surel - impl. of distributed source/sink
c 02.04.1996	ggu	$$exxqq - interpolation for more files ok
c 25.06.1997	ggu	complete restructured -> new subroutines
c 18.09.1997	ggu	$$FLUX3 - special treatment of type 3 boundary
c 23.09.1997	ggu	concentration boundary as file implemented
c 03.12.1997	ggu	$$TS - temp/salt implemented
c 20.05.1998	ggu	rain from file implemented
c 22.05.1998	ggu	local variable t introduced
c 22.05.1998	ggu	corrected bug for rain (surel called twice)
c 28.05.1998	ggu	new ruv, rvv (momentum input) (incomplete)
c 20.06.1998	ggu	more on momentum input (mvalue)
c 29.06.1998	ggu	!$$momin1 - bug fix for momentum input
c 13.07.1998	ggu	!initialize ruv, rvv by node
c 14.07.1998	ggu	finally momentum input finished -> new exxqq
c 20.08.1998	ggu	iextpo finally eliminated
c 20.08.1998	ggu	2d dependent routines copied to subini.f
c 21.08.1998	ggu	xv eliminated
c 24.08.1998	ggu	BC for concentration is bnd(20,..)
c 24.08.1998	ggu	BC for maximum input level is bnd(12,..) -> levmax
c 05.11.1998	ggu	slightly restructured restart
c 18.02.1999	ggu	allow for boundary type 0
c 20.04.1999	ggu	converted to stress instead of wind (tauxnv...)
c 01.12.1999	ggu	handle negative boundary type
c 07.05.2001	ggu	introduced variable zfact
c 24.10.2001	ggu	ignore boundary with type 0
c 10.08.2003	ggu	better commented, new routines meteo_init, meteo_force
c 14.08.2003	ggu	initialization of z and uv made explicit
c 10.03.2004	ggu	RQVDT - value in rqv is now discharge [m**3/s]
c 03.09.2004	ggu	restart taken out to ht
c 11.10.2004	ggu	new ibtyp=31, multiply with zfact also for sin (ZFACT)
c 11.03.2005	ggu	new boundary routines b3dvalue, c3dvalue (3D bounds)
c 17.02.2006	ggu	new routine get_bflux()
c 23.03.2006	ggu	changed time step to real
c 18.09.2007	ggu	new set_mass_flux, bug fix in dist_3d
c 25.09.2007	ggu	routines deleted: [mbc]value, testbc
c 02.10.2007	ggu	bug fix in make_scal_flux: surface flux only in layer 1
c 08.10.2007	ggu	deleted commented lines
c 17.03.2008	ggu	name of some variables changed, new scalar values
c 07.04.2008	ggu	deleted c2dvalue, set_scal_bc
c 07.04.2008	ggu	differential input introduced (-5555)
c 09.04.2008	ggu	only level boundary array, re-arranged
c 10.04.2008	ggu&ccf	new call to init_z0b()
c 17.04.2008	ggu	lmin introduced in set_mass_flux (negative levmax)
c 17.04.2008	ggu	evaporation introduced, rain adjusted
c 18.04.2008	ggu	rain simplified, bugfix
c 22.04.2008	ggu	in make_scal_flux do not alter mfluxv (parallel code)
c 03.06.2008	ggu	levmin introduced
c 24.03.2009	ggu	bug fix for rain; rain0d; new rain2distributed()
c 27.03.2009	ggu	new routine adjust_mass_flux() for dry nodes
c 31.03.2009	ggu	in make_scal_flux() bug fix for negative fluxes
c 02.04.2009	ggu	use get_bnd_(i)par() for special routines
c 03.04.2009	ggu	set intpol depending on ibtyp if intpol==0
c 20.04.2009	ggu	new routine z_tilt (tilting around a fixed value)
c 27.04.2009	ggu	can use ktilt with ztilt
c 19.01.2010	ggu	in make_scal_flux() return also sconz at bound
c 26.02.2010	ggu	bug fix in make_scal_flux() - sconz was out of l-loop
c 01.03.2010	dbf	tramp introduced
c 10.03.2010	ggu	new routine check_scal_flux()
c 22.03.2010	ggu	in make_scal_flux() forgotten to initialize sconz
c 14.04.2010	ggu	account for negative levmin/max
c 16.02.2011	ggu	copied meteo routines to submet.f
c 23.02.2011	ggu	new parameters tramp and levflx implemented
c 01.03.2011	ggu	implement smoothing for levflx
c 14.04.2011	ggu	changed VERS_6_1_22
c 14.05.2011	ggu	new routine get_discharge()
c 31.05.2011	ggu	changed VERS_6_1_23
c 07.06.2011	ggu	changed VERS_6_1_25
c 30.03.2012	ggu	changed VERS_6_1_51
c 13.06.2013	ggu	changed VERS_6_1_65
c 29.11.2013	ggu	prepared for ibtyp == 0
c 05.12.2013	ggu	changed VERS_6_1_70
c 05.05.2014	ggu	changed VERS_6_1_74
c 27.06.2014	ggu	changed VERS_6_1_78
c 07.07.2014	ggu	changed VERS_6_1_79
c 10.07.2014	ggu	only new file format allowed
c 18.07.2014	ggu	changed VERS_7_0_1
c 30.10.2014	ggu	in c_tilt() compute distance also for lat/lon
c 31.10.2014	ccf	new call to init_z0 instead than init_z0b
c 05.11.2014	ggu	changed VERS_7_0_5
c 07.11.2014	ggu	bug fix for distance computation in z_tilt, c_tilt
c 05.12.2014	ggu	changed VERS_7_0_8
c 19.12.2014	ggu	changed VERS_7_0_10
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 10.02.2015	ggu	new call to iff_read_and_interpolate()
c 26.02.2015	ggu	changed VERS_7_1_5
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 30.07.2015	ggu	changed VERS_7_1_83
c 05.11.2015	ggu	changed VERS_7_3_12
c 16.12.2015	ggu	changed VERS_7_3_16
c 18.12.2015	ggu	changed VERS_7_3_17
c 12.01.2017	ggu	changed VERS_7_5_21
c 05.12.2017	ggu	changed VERS_7_5_39
c 03.04.2018	ggu	changed VERS_7_5_43
c 04.04.2018	ggu	initialization for zeta revised for mpi
c 19.04.2018	ggu	changed VERS_7_5_45
c 01.06.2018	mbj	added tramp for sea level boundary condition
c 06.07.2018	ggu	changed VERS_7_5_48
c 11.10.2018	ggu	zconst setting adjourned
c 18.12.2018	ggu	changed VERS_7_5_52
c 16.02.2019	ggu	changed VERS_7_5_60
c 06.03.2019	mbj	added tramp = -2, (average, Bajo et al., 2019)
c 13.03.2019	ggu	changed VERS_7_5_61
c 04.07.2019	ggu	setweg & setznv introduced before make_new_depth
c 05.06.2020	ggu	bug fix: set_area was not called (MPI_SET_AREA)
c 30.03.2021	ggu	in set_mass_flux() use new time step data
c 31.03.2021	ggu	some debug code (iudbg)
c
c***************************************************************

	subroutine sp111(mode)

c set up boundary and initial conditions
c
c mode		1 : first call, initialize b.c.
c		2 : read in b.c.

	use mod_bound_geom
	use mod_bnd_aux
	use mod_bound_dynamic
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
	use intp_fem_file
	use shympi

	implicit none

	integer mode

	include 'mkonst.h'

	real rwv2(nkn)
	integer, save, allocatable :: ids(:)		!id values for BC

	integer nodes(nkn)
	real vconst(nkn)

	character*12 auxname
	logical bimpose
	integer kranf,krend,k,kn
	integer ibc,ibtyp
        integer nk,i,kk,kindex,iv,iw
        integer nsize,il
        integer iunrad,ktilt
	integer ip,l,lmax,ivar
	integer levflx
	integer nbc
	integer id,intpol,nvar,ierr
	integer iudbg
	double precision dtime0,dtime,ddtime
	real rw,zconst,aux
	real dt
	real conz,temp,salt
	real conzdf,tempdf,saltdf
c	real dz
	real rmu,rmv
	real getpar,rwint
	real conz3,temp3,salt3
	real tramp,alpha
	character*80 zfile
	logical, save ::  bdebug = .false.
	real, parameter :: zflag = -999.

	integer iround
        integer nkbnds,kbnds,itybnd,nbnds
	integer ipext,kbndind

	real, allocatable, save :: dzbav(:)

	call get_timestep(dt)

	if( mode .eq. 1 ) goto 1
	if( mode .eq. 2 ) goto 2

	write(6,*) 'mode = ',mode
	stop 'error stop: internal error sp111: mode'

c---------------------------------------------------------------
c first call %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c---------------------------------------------------------------

    1	continue

c	-----------------------------------------------------
c       initialize boundary ids
c	-----------------------------------------------------

	nbc = nbnds()
	allocate(ids(nbc))
	allocate(dzbav(nbc))	!This is very small - mbj
	ids = 0
	dzbav = zflag
	
c	-----------------------------------------------------
c       initialize meteo
c	-----------------------------------------------------

	call meteo_init

c	-----------------------------------------------------
c       initialize tilted and flux boundaries
c	-----------------------------------------------------

	nbc = nbnds()

	do ibc=1,nbc
	  ibtyp=itybnd(ibc)
          nk = nkbnds(ibc)

          bimpose = ibtyp .ge. 1 .and. ibtyp .le. 2

          if( bimpose .and. nk .le. 1 ) goto 95 !$$ibtyp3

	  if( bimpose ) then    		!$$FLUX3 - not for ibtyp=3
	    call init_tilt(ibc)           	!find nodes for tilting
	    call init_flux(ibc)   		!set up rlv,rhv,rrv,ierv
	  end if
	end do

c	-----------------------------------------------------
c       initialization of aux array
c	-----------------------------------------------------

c	-----------------------------------------------------
c       initialization of fem_intp
c	-----------------------------------------------------

	call get_first_dtime(dtime0)

	nvar = 1
	vconst = 0.
	ids = 0

	do ibc=1,nbc
          nk = nkbnds(ibc)
	  if( nk .le. 0 ) cycle
	  do i=1,nk
            nodes(i) = kbnds(ibc,i)
	  end do
	  call get_boundary_file(ibc,'zeta',zfile)
	  call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	  call get_bnd_ipar(ibc,'intpol',intpol)
	  if( intpol .le. 0 ) then
	    intpol = 2
	    if( ibtyp .eq. 1 ) intpol = 4
	  end if
	  write(6,'(a,3i5)') 'opening boundary file: (ibc,ibtyp,intpol) '
     +				,ibc,ibtyp,intpol
          call iff_init(dtime0,zfile,nvar,nk,0,intpol
     +                          ,nodes,vconst,id)
	  if( ibtyp .eq. 0 ) then
	    write(auxname,'(a8,1x,i3)') 'not used',ibtyp
	  else if( ibtyp .lt. 0 ) then
	    write(auxname,'(a8,1x,i3)') 'closed  ',ibtyp
	  else if( ibtyp .eq. 1 ) then
	    write(auxname,'(a8,1x,i3)') 'zlevel  ',ibtyp
	  else if( ibtyp .eq. 2 .or. ibtyp .eq. 3 ) then
	    write(auxname,'(a8,1x,i3)') 'disch   ',ibtyp
	  else
	    write(auxname,'(a8,1x,i3)') 'bound   ',ibtyp
	  end if
	  call iff_set_description(id,ibc,auxname)
	  ids(ibc) = id
	  write(6,'(a,2i5,4a)') ' boundary file opened: '
     +			,ibc,id,' ',auxname,' ',trim(zfile)
	end do

	!call iff_print_info(ids(1))

c	-----------------------------------------------------
c       determine constant for z initialization
c	-----------------------------------------------------

	zconst=getpar('zconst')	!constant initial z value
	dtime = dtime0
	ivar = 1
	lmax = 1

	do ibc=1,nbc
	  ibtyp=itybnd(ibc)
	  id = ids(ibc)
	  if( id .le. 0 ) cycle
	  if( ibtyp == 1 ) exit
	end do

	ibc = shympi_min(ibc)	!choose boundary with minimum index
	if( ibc > nbc ) ibc = 0

	if( zconst /= zflag ) then
	  !use this value
	else if( ibc > 0 ) then
	  ibtyp=itybnd(ibc)
	  nk = nkbnds(ibc)
	  id = ids(ibc)
	  call iff_read_and_interpolate(id,dtime)
	  call iff_time_interpolate(id,dtime,ivar,nk,lmax,rwv2)
	  call adjust_bound(id,ibc,dtime,nk,rwv2)
	  zconst = rwv2(1)
	else
	  zconst = 0.
	end if

	call putpar('zconst',zconst)

c	-----------------------------------------------------
c       initialize variables or restart
c	...the variables that have to be set are zenv, utlnv, vtlnv
c	-----------------------------------------------------

	call init_z(zconst)	!initializes znv and zenv
	call setweg(0,iw)	!sets minimum depth
	call setznv		!adjusts znv

	call set_area		!initializes area	!bugfix MPI_SET_AREA
	call make_new_depth	!initializes layer thickness

	call init_uvt		!initializes utlnv, vtlnv
	call init_z0		!initializes surface and bottom roughness

	!call uvint
	call init_uv		!computes velocities
	call copy_uvz		!copies to old time level
	call copy_depth		!copies layer thickness

c	-----------------------------------------------------
c       finish
c	-----------------------------------------------------

c next only for radiation condition
c
c        write(78,*) 46728645,1
c        write(78,*) 50,0
c        write(78,*) 0,50,-0.01,0.01
c        write(79,*) 46728645,1
c        write(79,*) 50,0
c        write(79,*) 0,50,-0.01,0.01

	return

c---------------------------------------------------------------
c normal call for boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%
c---------------------------------------------------------------

    2	continue

	call get_act_dtime(dtime)

c	-----------------------------------------------------
c	initialize node vectors with boundary conditions
c	-----------------------------------------------------

	do k=1,nkn
          rzv(k)=flag
          rqv(k)=0.	!$$rqv1 !$$close1	[m**3/s]
          rqpsv(k)=0.	!fluxes - point sources [m**3/s]
          rqdsv(k)=0.	!fluxes - distributed sources through surface [m**3/s]
          ruv(k)=0.	!momentum input
          rvv(k)=0.
	end do

c	-----------------------------------------------------
c	loop over boundaries
c	-----------------------------------------------------

        !call bndo_radiat(it,rzv)
	nbc = nbnds()

	ivar = 1
	lmax = 1

	iudbg = 789
	iudbg = 0
	if(iudbg>0) write(iudbg,*) dtime,nbc

	do ibc=1,nbc

          call get_bnd_ipar(ibc,'ibtyp',ibtyp)
          call get_bnd_ipar(ibc,'levflx',levflx)
          call get_bnd_par(ibc,'tramp',tramp)
	  id = ids(ibc)

	  if( ibtyp .eq. 0 ) cycle
	  if( id .le. 0 ) cycle

          nk = nkbnds(ibc)   !total number of nodes of this boundary

	  rmu = 0.
	  rmv = 0.

	  call iff_read_and_interpolate(id,dtime)
	  call iff_time_interpolate(id,dtime,ivar,nk,lmax,rwv2)
	  call adjust_bound(id,ibc,dtime,nk,rwv2)

	  if( abs(ibtyp) == 1 ) call setbnds(ibc,rwv2(1))	!for closure

          alpha = 1.
          if( tramp .gt. 0. ) then
            call get_passed_dtime(ddtime)
            alpha = ddtime/tramp
            if( alpha .gt. 1. ) alpha = 1.
	  else if( nint(tramp) == -1 ) then	!subtract single values
            alpha = -1.
	  else if( nint(tramp) == -2 ) then	!subtract average over ibc bound
	    alpha = -2.
	    call z_ramp_avinit(nbc,ibc,nk,rwv2,dzbav(ibc))
	  end if

	  do i=1,nk

             kn = kbnds(ibc,i)
	     rw = rwv2(i)
	     !write(6,*) ibc,i,kn,ipext(kn),rw
	     if( kn <= 0 ) cycle

	     if(ibtyp.eq.1) then		!z boundary
               if (nint(tramp) /= 0) call z_ramp(kn,alpha,dzbav(ibc),rw)
               rzv(kn) = rw
	     else if(ibtyp.eq.2) then		!q boundary
	       call level_flux(dtime,levflx,kn,rw)	!zeta to flux
	       kindex = kbndind(ibc,i)
               rqpsv(kn)=alpha*rw*rrv(kindex)	!BUGFIX 21-08-2002, RQVDT
             else if(ibtyp.eq.3) then		!$$ibtyp3 - volume inlet
               rqpsv(kn) = alpha*rw
             else if(ibtyp.eq.4) then		!momentum input
	       ruv(kn) = rmu
	       rvv(kn) = rmv
             else if(ibtyp.eq.31) then		!zero gradient for z
               ! already done...
	       call get_bnd_ipar(ibc,'ktilt',ktilt)
               rzv(ktilt)=rw
               if( ibc .eq. 1 ) then
                 iunrad = 78
               else
                 iunrad = 79
               end if
               if( i .eq. 1 ) write(iunrad,*) dtime,nk,ktilt,rw
               write(iunrad,*) rzv(kn)
             else if(ibtyp.eq.32) then		!for malta 
c	       nothing
             else if(ibtyp.eq.0) then		!switched off
c	       nothing
             else if(ibtyp.lt.0) then		!closed
c	       nothing
	     else
c               kranf,krend not available...
c	       call zspeci(ibtyp,kranf,krend,rw)	!for radiation...
	       write(6,*) 'boundary = ',ibc,'   type = ',ibtyp
	       stop 'error stop sp111: Unknown boundary type'
	     end if

	     if(iudbg>0) write(iudbg,*) ibc,i,id,ibtyp,rw
	  end do

	end do

	!write(99,*) '======================================'
	!write(99,*) ' it = ',it
	!write(99,*) '======================================'
	!call iff_print_info(12,0,.true.)

c	-----------------------------------------------------
c	tilting
c	-----------------------------------------------------

	call z_tilt
	call c_tilt

c	-----------------------------------------------------
c	meteo forcing					!$$surel
c	-----------------------------------------------------

	call meteo_force

c	-----------------------------------------------------
c	set mass flux -> fills mfluxv and integrates to rqv
c	-----------------------------------------------------

	call set_mass_flux

c	-----------------------------------------------------
c	testing
c	-----------------------------------------------------

c	call tsbnds

c -----------------------------------------------------------
c end of routine
c -----------------------------------------------------------


	return
   95	continue
	write(6,*) 'One node boundary not allowed'
	write(6,*) 'Boundary :',ibc
	write(6,*) 'type     :',ibtyp
	stop 'error stop : sp111'
	end

c**************************************************************

	subroutine z_ramp(kn,alpha,bav,rw)

! implements ramping of water levels for Ensemble Kalman Filter

	use basin, only : nkn
	use mod_hydro, only : znv

	implicit none

	integer, intent(in) :: kn
	real, intent(in) :: alpha, bav
	real, intent(inout) :: rw

	real, parameter :: flag = -999.
	integer, save :: icall = 0
	real, save, allocatable ::  dzb(:)

	if( icall == 0 ) then
	  allocate(dzb(nkn))
	  dzb = flag
	  icall = 1
	end if

        if ( alpha >= 0 ) then		! from the z-restart to the z-file
          rw = znv(kn) + alpha * (rw - znv(kn))
        else if ( nint(alpha) == -1 ) then	
	  ! subtract the initial bias from the z-file records
	  if ( dzb(kn) == flag ) dzb(kn) = znv(kn) - rw
          rw = rw + dzb(kn)
	else if ( nint(alpha) == -2 ) then
	  ! subtract the initial bias between the averages from z-file records
	  if( nint(bav) == flag ) 
     +      stop 'error stop z_ramp: error in average boundary value.'
	  rw = rw + bav
        end if

	end

c**************************************************************

	subroutine z_ramp_avinit(nbc,ib,nkb,rwv2,bav)

	use basin, only : nkn
	use mod_hydro, only : znv

	implicit none

	integer, intent(in) :: nbc,ib,nkb
	real, intent(in) :: rwv2(nkn)
	real, intent(out) :: bav

	integer i,kn,kbnds
	real rw
	integer, save :: icall = 0

	if( icall /= 0 ) return

	bav = 0.
	do i = 1,nkb
	   kn = kbnds(ib,i)
	   rw = rwv2(i)
	   bav = bav + znv(kn) - rw
	end do
	bav = bav/nkb

	if (ib == nbc) icall = 1

	end subroutine z_ramp_avinit


c**************************************************************

	subroutine z_tilt

c artificial tilting of boundary surface - uses ktilt and ztilt
c
c the first boundary node is set to -ztilt, and the last to +ztilt
c the total water level difference is therefore 2*ztilt
c if ktilt is not given then the other nodes are linearily interpolated
c	between these two values
c if ktilt is given then this node will be set to z=0 and the other
c	nodes are linearly interpolated between start-ktilt and ktilt-end

	use mod_bound_geom
	use mod_bound_dynamic
	use basin

	implicit none

	integer ibc,ibtyp,ktilt
	integer nbc
	integer k,kranf,krend,kn1,kn2
	real dx,dy,ztilt,z
	double precision rltot,rltot1,rltot2,rl

	integer itybnd,nbnds

	nbc = nbnds()

	do ibc=1,nbc
          ibtyp = itybnd(ibc)
	  if( ibtyp <= 0 ) cycle
          call kanfend(ibc,kranf,krend)
	  call get_bnd_par(ibc,'ztilt',ztilt)
	  call get_bnd_ipar(ibc,'ktilt',ktilt)
	  if( ztilt .ne. 0. .and. ibtyp .eq. 1 ) then
	    rltot = 0.
	    rltot1 = 0.
	    do k=kranf+1,krend
		kn2=irv(k)
		kn1=irv(k-1)
		call compute_distance(xgv(kn1),ygv(kn1)
     +				,xgv(kn2),ygv(kn2),dx,dy)
	        rltot = rltot + sqrt(dx*dx+dy*dy)
		if( k .eq. ktilt ) rltot1 = rltot	!BUG 3.12.2013
	    end do
	    rltot2 = rltot - rltot1

c in rltot the whole length of boundary is stored
c in rltot1 and rltot2 are first and second part of length of boundary
c rltot1 from start to ktilt, and rltot2 from ktilt to end
c if no ktilt is given rltot1/rltot2 are not used

	    rl = 0.
            rzv(irv(kranf)) = rzv(irv(kranf)) - ztilt
	    do k=kranf+1,krend
		kn2=irv(k)
		kn1=irv(k-1)
		call compute_distance(xgv(kn1),ygv(kn1)
     +				,xgv(kn2),ygv(kn2),dx,dy)
	        rl = rl + sqrt(dx*dx+dy*dy)
		if( ktilt .le. 0 ) then			!no fixed node
	          z = - ztilt + (rl/rltot) * 2. * ztilt
		else
		  if( rl .le. rltot1 ) then
	            z = - ztilt + (rl/rltot1) * ztilt
		  else
	            z = ((rl-rltot1)/rltot2) * ztilt
		  end if
		end if
                rzv(kn2) = rzv(kn2) + z
	    end do
	  end if
	end do

	end

c**************************************************************

	subroutine c_tilt

c tilting of boundary surface due to Coriolis acceleration - needs ktilt
c
c if ztilt is given then z_tilt() is used to tilt water level
c if ktilt is not given nothing is tilted

	use mod_bound_geom
	use mod_bound_dynamic
	use mod_hydro_print
	use mod_hydro
	use basin

	implicit none

	include 'pkonst.h'
	include 'mkonst.h'

	integer ibc,ibtyp,kranf,krend,ktilt,k,kn1,kn2
	integer nbc
	real roinv,f,ginv,dx,dy,taux,tauy
	real taux1,taux2,tauy1,tauy2,wx,wy
	real u,v,z,h,p1,p2,b,hh,ztilt
	real getpar
	integer iround,itybnd,nbnds

	roinv=1./rowass
	f=fcor
	ginv=1./grav

	nbc = nbnds()

	do ibc=1,nbc
         ibtyp = itybnd(ibc)
	 if( ibtyp <= 0 ) cycle
         call kanfend(ibc,kranf,krend)
	 call get_bnd_ipar(ibc,'ktilt',ktilt)
	 call get_bnd_par(ibc,'ztilt',ztilt)

	 if( ztilt .ne. 0 ) then	!is handled in z_tilt
		!nothing
	 else if(ktilt.gt.0.and.ibtyp.eq.1) then
	   do k=ktilt+1,krend,1
		kn2=irv(k)
		kn1=irv(k-1)
		call compute_distance(xgv(kn1),ygv(kn1)
     +				,xgv(kn2),ygv(kn2),dx,dy)
		call get_meteo_forcing(kn1,wx,wy,taux1,tauy1,p1)
		call get_meteo_forcing(kn2,wx,wy,taux2,tauy2,p2)
		taux = 0.5*(taux1+taux2)
		tauy = 0.5*(tauy1+tauy2)
		u=(up0v(kn1)+up0v(kn2))*.5
		v=(vp0v(kn1)+vp0v(kn2))*.5
		z=(znv(kn1)+znv(kn2))*.5
		h=(rhv(k-1)+rhv(k))*.5
		b=1./(h+z)
		hh=-roinv*(p2-p1)+b*(taux*dx+tauy*dy)+f*(v*dx-u*dy)
		rzv(kn2)=rzv(kn1)+ginv*hh
	   end do

	   do k=ktilt-1,kranf,-1
		kn2=irv(k+1)
		kn1=irv(k)
		call compute_distance(xgv(kn1),ygv(kn1)
     +				,xgv(kn2),ygv(kn2),dx,dy)
		call get_meteo_forcing(kn1,wx,wy,taux1,tauy1,p1)
		call get_meteo_forcing(kn2,wx,wy,taux2,tauy2,p2)
		taux = 0.5*(taux1+taux2)
		tauy = 0.5*(tauy1+tauy2)
		u=(up0v(kn1)+up0v(kn2))*.5
		v=(vp0v(kn1)+vp0v(kn2))*.5
		z=(znv(kn1)+znv(kn2))*.5
		h=(rhv(k+1)+rhv(k))*.5
		b=1./(h+z)
		hh=-roinv*(p2-p1)+b*(taux*dx+tauy*dy)+f*(v*dx-u*dy)
		rzv(kn1)=rzv(kn2)-ginv*hh
	   end do
	 end if
	end do

	end

c**************************************************************

	subroutine init_tilt(ibc)

c finds tilting node in boundary node list

	implicit none

	integer ibc

	logical berr
	integer kranf,krend
	integer ktilt,i
        integer kb
	integer ibtyp

	integer ipext,iround,kbnd,itybnd

        ibtyp = itybnd(ibc)
	if( ibtyp <= 0 ) return

	call get_bnd_ipar(ibc,'ktilt',ktilt)
	if(ktilt.le.0) return

	call kanfend(ibc,kranf,krend)

	berr = .true.
	do i=kranf,krend
           kb = kbnd(i)
	   if( kb .eq. ktilt ) then
		call set_bnd_ipar(ibc,'ktilt',i)
		berr = .false.
	   end if
	end do

	if( berr ) then
	  write(6,*) 'Node number for tilting not in boundary node list'
	  write(6,*) 'ktilt :',ipext(ktilt)
	  stop 'error stop init_tilt: no node'
	end if

	end

c******************************************************************

	subroutine init_flux(ibc)

c initializes flux boundary

	use mod_bound_geom
	use basin

	implicit none

        integer ibc

	integer kranf,krend
	integer ie,i,k1,k2,kk1,kk2,ii1,ii2
	integer ibtyp
	real fm,dx,dy,rl,h1,h2,fl
	integer ipext,itybnd

        ibtyp = itybnd(ibc)
	if( ibtyp <= 0 ) return

        call kmanfend(ibc,kranf,krend)

	if( krend-kranf .le. 0 ) return

	ii1 = 0
	ii2 = 0
	kk1 = 0
	kk2 = 0

	fm=0
	do i=kranf,krend	!these are already set to 0 at initialization
	   rrv(i)=0.
	   rhv(i)=0.
	end do

	do i=kranf,krend-1

	  k1=irv(i)
	  k2=irv(i+1)

	  do ie=1,nel
	   do ii1=1,3
	      ii2=mod(ii1,3)+1
	      kk1=nen3v(ii1,ie)
	      kk2=nen3v(ii2,ie)
	      if(k1.eq.kk1.and.k2.eq.kk2) goto 33
	   end do
	  end do

   33	  continue

	  if( k1 .ne. kk1 .or. k2 .ne. kk2 ) then
	    write(6,*) 'Cannot locate boundary nodes in element index'
	    write(6,*) 'node 1,node 2 :',ipext(k1),ipext(k2)
	    write(6,*) '(Are you sure that boundary nodes are given'
	    write(6,*) '   in anti-clockwise sense ?)'
	    stop 'error stop init_flux: no node'
	  end if

	  dx=xgv(k1)-xgv(k2)
	  dy=ygv(k1)-ygv(k2)
	  rl=sqrt(dx*dx+dy*dy)
	  h1=hm3v(ii1,ie)
	  h2=hm3v(ii2,ie)
	  fl=rl*(h1+h2)/2.
	  fm=fm+fl

	  rrv(i)=rrv(i)+rl*(2.*h1+h2)/6.
	  rrv(i+1)=rrv(i+1)+rl*(h1+2.*h2)/6.

	  rlv(i)=rl

	  rhv(i)=rhv(i)+h1
	  rhv(i+1)=rhv(i+1)+h2

	  ierv(1,i)=ie
	  ierv(2,i)=mod(ii2,3)+1

	end do

	do i=kranf,krend
	   rrv(i)=rrv(i)/fm
	end do
	do i=kranf+1,krend-1
	   rhv(i)=rhv(i)/2.
	end do

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine adjust_bound(id,ibc,dtime,nk,rw)

	use intp_fem_file

	implicit none

	integer id
	integer ibc
	double precision dtime
	integer nk
	real rw(nk)

	integer i
	real rit,rw0,zfact

	call get_bnd_par(ibc,'zfact',zfact)

	if( iff_has_file(id) ) then
	  do i=1,nk
	    rw(i) = rw(i) * zfact
	  end do
	else
	  call get_oscil(ibc,dtime,rw0)
	  do i=1,nk
	    rw(i) = rw0 * zfact
	  end do
	end if

	end

c*******************************************************************

	subroutine set_mass_flux

c sets up (water) mass flux array mfluxv (3d) and rqv (vertically integrated)

	use mod_bound_dynamic
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	logical debug,bdebug
	integer i,k,l,lmin,lmax,nk,ibc,mode,iu
	integer ibtyp,levmax,levmin
	integer nbc
	double precision dtime
	real flux,vol,voltot,fluxtot,fluxnode
	real vols(nkn)

	integer nkbnds,kbnds,nbnds
	real volnode		!function to compute volume of node

c------------------------------------------------------------------
c initialize arrays and parameter
c------------------------------------------------------------------

	mode = -1		!old time step
	mode = +1		!old time step
	debug = .true.
	debug = .false.

	iu = 0
	call get_act_dtime(dtime)
	debug = ( iu > 0 .and. dtime >= 86400. )
	if( debug ) write(iu,*) 'computing mass flux at time: ',dtime

	mfluxv = 0.

c------------------------------------------------------------------
c loop over boundaries for point sources -> use volumes as weight
c------------------------------------------------------------------

	nbc = nbnds()

	do ibc=1,nbc

          nk = nkbnds(ibc)
	  call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	  call get_bnd_ipar(ibc,'levmin',levmin)
	  call get_bnd_ipar(ibc,'levmax',levmax)

	  if( ibtyp .lt. 2 .or. ibtyp .gt. 3 ) nk = 0		!skip

	  if(debug) then
	    write(iu,*) 'computing mass flux: ',ibc,ibtyp,nk,levmax
	  end if

	  do i=1,nk
            k = kbnds(ibc,i)
	    if( k <= 0 ) cycle
	    bdebug = ( debug .and. k == 3239 ) 
	    lmax = ilhkv(k)
	    if( levmax .gt. 0 ) lmax = min(lmax,levmax)
	    if( levmax .lt. 0 ) lmax = min(lmax,lmax+1+levmax)
	    lmin = 1
	    if( levmin .gt. 0 ) lmin = max(lmin,levmin)
	    if( levmin .lt. 0 ) lmin = max(lmin,lmax+1+levmin)

	    if( lmin .gt. lmax ) goto 98

	    voltot = 0.
	    do l=lmin,lmax
	      vol = volnode(l,k,mode)
	      vols(l) = vol
	      voltot = voltot + vol
	    end do

	    if( voltot .le. 0. ) goto 99

	    flux = rqpsv(k)
	    if(bdebug) write(iu,*) '   ',k,lmin,lmax,flux,voltot
	    do l=lmin,lmax
	      mfluxv(l,k) = flux * vols(l) / voltot
	    end do
	  end do

	end do

c------------------------------------------------------------------
c add distributed sources
c------------------------------------------------------------------

	do k=1,nkn
	  mfluxv(1,k) = mfluxv(1,k) + rqdsv(k)	!rain, evaporation
	  !lmax = ilhkv(k)
	  !mfluxv(lmax,k) = gwf		!here distributed ground water flow
	end do

c------------------------------------------------------------------
c compute total flux for check and integrate flux into rqv
c------------------------------------------------------------------

	fluxtot = 0.
	do k=1,nkn
	  fluxnode = 0.
	  lmax = ilhkv(k)
	  do l=1,lmax
	    flux = mfluxv(l,k)
	    !if( debug .and. flux .gt. 0. ) write(6,*) '  flux: ',k,l,flux
	    fluxnode = fluxnode + flux
	  end do
	  rqv(k) = fluxnode
	  fluxtot = fluxtot + fluxnode
	end do

	if( debug ) write(iu,*) '  total flux: ',fluxtot

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	return
   98	continue
	write(6,*) 'lmin > lmax'
   99	continue
	write(6,*) 'ibc = ',ibc
	write(6,*) 'i = ',i
	write(6,*) 'k = ',k
	write(6,*) 'ilhkv(k) = ',ilhkv(k)
	write(6,*) 'levmin = ',levmin
	write(6,*) 'levmax = ',levmax
	write(6,*) 'lmin = ',lmin
	write(6,*) 'lmax = ',lmax
	stop 'error stop set_mass_flux: voltot = 0'
	end

c**********************************************************************

	subroutine adjust_mass_flux

c adjusts mass flux for dry nodes

	use mod_geom_dynamic
	use mod_bound_dynamic
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,l

        logical iskout
        iskout(k) = inodv(k).eq.-2

	do k=1,nkn
	  if( iskout(k) ) then
	    do l=1,nlv
	      mfluxv(l,k) = 0.
	    end do
	    rqv(k) = 0.
	    rqdsv(k) = 0.
	  end if
	end do

	end

c**********************************************************************

	subroutine make_scal_flux(what,r3v,scal,sflux,sconz,ssurf)

c computes scalar flux from fluxes and concentrations

	use mod_bound_dynamic
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) what
	real r3v(nlvdi,nkn)	!concentration for boundary condition
	real scal(nlvdi,nkn)	!concentration of scalar
	real sflux(nlvdi,nkn)	!mass flux for each finite volume (return)
	real sconz(nlvdi,nkn)	!concentration for each finite volume (return)
	real ssurf		!value of scalar for surface flux

	include 'mkonst.h'


	integer k,l,lmax,ks
	real flux,conz
	real surf_flux
	real getpar

	ks = 2827
	ks = 2757
	ks = 2831
	ks = -1
	!ks = nint(getpar('kref'))	!not working - here global, but local

	do k=1,nkn
	  lmax = ilhkv(k)
	  surf_flux = rqdsv(k)
	  do l=1,lmax
	    sflux(l,k) = 0.
	    sconz(l,k) = 0.
	    flux = mfluxv(l,k)
	    if( l .eq. 1 ) flux = flux - surf_flux	!without surface flux
	    conz = r3v(l,k)
	    if( flux .ne. 0. .or. conz .ne. flag ) then
	      if( flux .ne. 0. .and. conz .eq. flag ) goto 99
	      if( conz .le. -990. ) conz = scal(l,k)	!ambient value
	      if( conz .le. -5555. ) conz = scal(l,ks) - 10000. - conz !diff
	      if( flux .lt. 0. ) conz = scal(l,k)	!ambient value (BUG)
	      sflux(l,k) = flux * conz
	      sconz(l,k) = conz			!bug fix - was out of loop
	    end if
	  end do
	  conz = ssurf
	  if( ssurf .le. -990 ) conz = scal(1,k)
	  if( ssurf .le. -5555 ) conz = scal(1,k) - 10000. - ssurf !diff
	  sflux(1,k) = sflux(1,k) + surf_flux * conz
	  ! next should be sconz(1,k) = conz if surf_flux is eliminated
	  !sconz(1,k) = 0.
	  if( mfluxv(1,k) .ne. 0 ) then
	    sconz(1,k) = sflux(1,k) / mfluxv(1,k)
	  end if
	end do

	!if( what .eq. 'salt' ) then
	!  call flux_debug(what,mfluxv,sflux,sconz)
	!end if

	return
   99	continue
	write(6,*) 'what: ',what
	write(6,*) 'k,l: ',k,l
	write(6,*) 'flux,conz: ',flux,conz
	stop 'error stop make_scal_flux: boundary condition mismatch'
	end

c**********************************************************************

	subroutine flux_debug(what,mfluxv,sflux,sconz)

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) what
	real mfluxv(nlvdi,nkn)	!mass flux
	real sflux(nlvdi,nkn)	!scalar flux
	real sconz(nlvdi,nkn)	!concentration for each finite volume

	integer k,l,lmax
	integer ifemop
	real qtot,stot
	double precision dtime

	integer, save :: iunit = 0

	call get_act_dtime(dtime)

	if( iunit .eq. 0 ) then
	  iunit = ifemop('.ggg','formatted','unknown')
	end if

	do k=1,nkn
	  lmax = ilhkv(k)
	  stot = 0.
	  qtot = 0.
	  do l=1,lmax
	    qtot = qtot + mfluxv(l,k)
	    stot = stot + mfluxv(l,k) * sconz(l,k)
	  end do
	  if( qtot .ne. 0 ) then
	    write(iunit,*) dtime,k,qtot,stot
	  end if
	end do

	end

c**********************************************************************

	subroutine check_scal_flux(what,scal,sconz)

c checks scalar flux

	use mod_bound_dynamic
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) what
	real scal(nlvdi,nkn)	!concentration of scalar
	real sconz(nlvdi,nkn)	!concentration for each finite volume

	include 'mkonst.h'

	integer k,l,lmax,ks
	real cconz,qflux,mflux
	double precision dtime

	call get_act_dtime(dtime)

	write(46,*) 'check_scal_flux ',what,dtime

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cconz = sconz(l,k)         !concentration has been passed
            qflux = mfluxv(l,k)
            if( qflux .lt. 0. ) cconz = scal(l,k)
            mflux = qflux * cconz
	    if( qflux .ne. 0 ) then
	      write(146,1000) k,l,mflux,qflux,cconz,scal(l,k)
	    end if
          end do
        end do

	return
 1000	format(2i10,4f10.4)
	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine init_scal_bc(r3v)

c initializes array for scalar boundary condition

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real r3v(nlvdi,nkn)

	include 'mkonst.h'

	integer k,l

	do k=1,nkn
	  do l=1,nlv
	    r3v(l,k) = flag
	  end do
	end do

	end

c*******************************************************************

	subroutine mult_scal_bc(r3v,value)

c multiplies array for scalar boundary condition with value

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real r3v(nlvdi,nkn)
	real value

	include 'mkonst.h'

	integer k,l

	do k=1,nkn
	  do l=1,nlv
	    if( r3v(l,k) .ne. flag ) r3v(l,k) = r3v(l,k) * value
	  end do
	end do

	end


c*******************************************************************

	subroutine dist_3d(nlvddi,r3v,kn,nbdim,values)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nlvddi
	real r3v(nlvddi,nkn)
	integer kn
	integer nbdim
	real values(nkn)

	integer l,lmax

	if( kn <= 0 ) return

	if( nbdim .eq. 0 ) then
	  lmax = 1
	else
	  lmax = min(nbdim,nlvddi)
	end if

	do l=1,lmax
	  r3v(l,kn) = values(l)
	end do
	  
	do l=lmax+1,nlvddi
	  !r3v(l,kn) = values(nbdim)	!BUGFIX
	  r3v(l,kn) = r3v(lmax,kn)
	end do
	  
	end

c**********************************************************************

	subroutine dist_horizontal(nlvddi,r3v,n,value)

	implicit none

	integer nlvddi
	real r3v(nlvddi,n)
	integer n
	real value

	integer k

	do k=1,n
	  r3v(1,k) = value
	end do
	  
	end

c**********************************************************************

        subroutine aver_horizontal(nlvddi,r3v,n,value)

        implicit none

        integer nlvddi
        real r3v(nlvddi,n)
        integer n
        real value

        integer k

	value = 0.
        do k=1,n
          value = value + r3v(1,k)
        end do
	value = value / n

	end

c**********************************************************************

	subroutine print_scal_bc(r3v)

c prints non-flag entries of array for scalar boundary condition

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real r3v(nlvdi,nkn)

	include 'mkonst.h'

	integer k,l
	real value

	do k=1,nkn
	  do l=1,nlv
	    value = r3v(l,k)
	    if( value .ne. flag ) then
		write(6,*) 'print_scal_bc: ',k,l,value
	    end if
	  end do
	end do

	end

c**********************************************************************

	subroutine get_bflux(k,flux)

c returns boundary flux of node k

	use mod_bound_dynamic

	implicit none

	integer k	!node
	real flux	!flux in node k (return)

	flux = rqv(k)

	end

c**********************************************************************

	subroutine level_flux(dtime,levflx,kn,rw)

c compute discharge from water level

	use mod_hydro

	implicit none

	double precision dtime	!time
	integer levflx		!type of function
	integer kn		!node number
	real rw			!discharge computed

	real z,a,b,c,z0

c use fit 1 for compatibility to old simulations
c use fit 3 for best results without using piecewise fit

	if( levflx .eq. 0 ) then
	  !nothing changed -> return same rw
	else if( levflx .eq. 1 ) then
	  z = znv(kn)
	  call z_smooth(z)

	  !a = 1648		!no outliers
	  !c = -2684		!no outliers

c-------------- old fit ----------------------------	1
	  z0 = 5.		!reference level
	  a = 1808
	  c = -2875
	  rw = a*log(z0+z) + c
c---------------------------------------------------

c-------------- fit based on mass balance -------------  2
	  !a = 211.521555200211
	  !b = 0.855968510497031     
	  !rw = a * z**b
c---------------------------------------------------

c-------------- new fit based on mass balance ----------  3
	  !a = 223.529239458563		!new calibration
	  !b = 0.862816507081288
	  !rw = a * z**b
c---------------------------------------------------

c-------------- best fit based on mass balance -----------  4
          !call flow_out_piece_new(z,rw)	!this is best (if it works)
c---------------------------------------------------

	  rw = -rw
	  write(134,*) dtime,z,rw
	else
	  write(6,*) 'levflx = ',levflx
	  stop 'error stop level_flux: levflx'
	end if

	end

c**********************************************************************

	subroutine z_smooth(z)

c smooths z values

	implicit none

	integer ndim
	parameter (ndim=86400/100)

	real z

	integer i

	integer ipz
	real za(0:ndim)
	save ipz,za
	data ipz / 0 /

	if( ipz .le. 0 ) then	!initialize
	  do i=0,ndim
	    za(i) = z
	  end do
	  ipz = 1
	end if

	za(0) = za(0) + (z-za(ipz))/ndim
	za(ipz) = z
	ipz = mod(ipz,ndim) + 1

	z = za(0)

	end

c**********************************************************************

        subroutine flow_out_piece_new(z,rout)

        implicit none

        real z,rout
        real rw,a,b

        if( z .le. 1. ) then
          a = 100.594731852550
          b = 78.7626400203164
        else if( z .le. 2. ) then
          a = -123.192333413729
          b = 285.797516554523
        else if( z .le. 3. ) then
          a = 21.2143155529982
          b = 207.240814424064
        else
          a = 468.988395301776
          b = 102.075378828782
        end if

        rw = a + z*b
        !rout = -rw
        rout = rw

        end

c**********************************************************************
 
	function get_discharge(ibc)

c returns discharge through boundary ibc for points sources
c for z-boundaries 0 is returned

	use mod_bound_dynamic

	implicit none

	real get_discharge
	integer ibc

	integer itype,nk,i,k
	real acc

	integer itybnd,nkbnds,kbnds

	get_discharge = 0.

	itype = itybnd(ibc)
	if( itype .le. 1 .or. itype .gt. 3 ) return

	acc = 0.
        nk = nkbnds(ibc)
        do i=1,nk
          k = kbnds(ibc,i)
	  acc = acc + rqpsv(k)
        end do

	get_discharge = acc

	end

c**********************************************************************

