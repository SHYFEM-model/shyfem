
!--------------------------------------------------------------------------
!
!    Copyright (C) 1988,1991-1993,1995-2001,2003-2020  Georg Umgiesser
!    Copyright (C) 2006,2011,2019  Christian Ferrarin
!    Copyright (C) 2007,2013  Debora Bellafiore
!    Copyright (C) 2008  Marco Bajo
!    Copyright (C) 2015  Erik Pascolo
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

c assembling linear system routine
c
c contents :
c
c subroutine hydro			administrates one time step
c subroutine hydro_zeta(vqv)		assemble matrix
c subroutine hydro_transports		computes transports (temporary)
c subroutine hydro_transports_final	computes transports (final)
c subroutine hydro_vertical(dzeta)	computes vertical velocities
c subroutine correct_zeta(dzeta)	corrects zeta values
c
c notes :
c
c look for "!ccc" to see important changes
c
c ASYM			passage to asymmetrix matrix
c ASYM_OPSPLT		new version without operator splitting (can leave)
c ASYM_OPSPLT_CH	new version -> change here
c
c nknddi	dimension for total number of nodes
c nelddi	dimension for total number of elements
c nrbddi	dimension for total number of boundary condition nodes
c nbcddi	dimension for total number of open boundaries
c mbwddi	dimension for bandwidth
c ngrddi	dimension for grade of nodes (number of elements attached
c		...to one node)
c narddi	dimension for total number of area codes
c nexddi	dimension for total number of extra nodes for output
c
c nkn,nel	total number of nodes/elements
c nrz,nrq	total number of nodes with water level/flux boundary conditions
c nrb		total number of nodes with boundary conditions
c nbc		total number of open boundaries
c ngr		maximum grade of nodes
c mbw		bandwidth of system matrix
c flag		flag value (to recognize boundary conditions)
c grav,dcor	gravitational accel./medium latitude of basin
c rowass,roluft	density of water/air
c idt,nits	time step/total iterations to go
c niter,it	actual number of iterations/actual time
c nlvdi,nlv	dimension for levels, number of used levels
c
c nen3v(..,ie)	element index - node numbers (3) of element ie
c
c ipv(k)	external node number of node k
c ipev(ie)	external element number of element ie
c
c rqv(k),rzv(k)	flux/water level boundary conditions for node k
c		...if(rzv(k).eq.flag) --> no b.c. is specified for node k
c bnd(..,ib)	specification of boundary conditions for open boundary ib
c ev(..,ie)	geometric parameters for element ie
c		...1-3 = area, 4-6 = b, 7-9 = c, 10 = Aomega, 11-13 = angle
c uov(ie)	depth integrated transport (old time level)
c vov(ie)	...
c unv(ie)	depth integrated transport (new time level)
c vnv(ie)	...
c zov(k)	water level of old/new time level
c znv(k)	...
c ulov(l,ie)	velocity of old time level in x direction of layer l and elem ie
c vlov(l,ie)	velocity of old time level in y direction of layer l and elem ie
c wlov(l,i)	velocity of old time level in z direction of layer l and node i
c ulnv(l,ie)	velocity of new time level in x direction of layer l and elem ie
c vlnv(l,ie)	velocity of new time level in y direction of layer l and elem ie
c wlnv(l,i)	velocity of new time level in z direction of layer l and node i
c utlov(l,ie)	transport of old time level in x direction of layer l and el. ie
c vtlov(l,ie)	transport of old time level in y direction of layer l and el. ie
c utlnv(l,ie)	transport of new time level in x direction of layer l and el. ie
c vtlnv(l,ie)	transport of new time level in y direction of layer l and el. ie
c uprv(l,k)	velocity (averaged) in x direction of layer l and node k
c vprv(l,k)	velocity (averaged) in y direction of layer l and node k
c wprv(l,k)	velocity (averaged) in z direction of layer l and node k
c up0v(k)	total velocity (averaged) in x direction of node k
c vp0v(k)	total velocity (averaged) in y direction of node k
c tauxnv(k)	normalized stress in x-direction at node k
c tauynv(k)	normalized stress in y-direction at node k
c rhov(l,k)	density for level l and node k
c fcorv(ie)	coriolis parameter for elem ie
c
c visv(l,k)	vertical turbulent viscosity for layer l and node k (alv)
c difv(l,k)	vertical turbulent diffusivity for layer l and node k (slv)
c
c xgv(k),ygv(k)	coordinates of node k
c
c hldv(l)	thickness of layer l
c hlv(l)	absolute depth of bottom of layer l
c ilhv(ie)	number of levels for element ie
c ilhkv(k)	number of levels for node k
c hlhv(ie)	thickness of last layer for element ie
c hev(ie)	total depth at element ie (no water level)
c hkv(k)	total depth at node k (no water level)
c hm3v(3,ie)	depth of three nodes in element ie
c
c v1v,v2v...	auxiliary vectors
c rmat		band matrix for one vertical system (as above)
c rvec		constant vector for vertical system
c
c $$h1new	use new water level to compute velocities
c		(only for printing velocities important, but
c		algorithm has to be checked if consistent)
c $$rtmax	use maximal friction coefficient of rdt (=1./dt)
c
c revision log :
c
c 27.07.1988	ggu	(from sp159f)
c 18.02.1991	ggu	(from scratch)
c 04.06.1991	ggu	(c=(1) : friction term has been corrected) (hydro_zeta)
c 27.08.1991	ggu	(from scratch) (hydro_vertical)
c 01.10.1992	ggu	(staggered FE - completely restructured) (hydro_zeta)
c 01.07.1993	ggu	$$UVBARO - u/vov introduced for	iteration on rad cond
c 03.11.1993	ggu	$$cmplerr - compiler warnings hydro
c 05.11.1993	ggu	$$fric - normal friction
c 05.11.1993	ggu	$$crador - crador call commented
c 05.11.1993	ggu	subroutine crador in file newcra.f
c 05.11.1993	ggu	$$VBARO-ERR - unv(ie)=vov(ie)
c 28.08.1995	ggu	$$BAROC_AREA - do baroc only for iarv(ie) = 0
c 30.08.1995	ggu	$$AUST - austausch coefficient introduced
c 01.09.1995	ggu	$$AWEIGH - area weighting of austausch coefficient
c 06.03.1996	ggu	$$BAROC_AREA0 - introduced baroc0
c 06.03.1996	ggu	$$VERT_AUST_ADJUST - adjustment of vert. aust. coef.
c 06.06.1996	ggu	$$BCHAO - modifications for vel. profile (temp.)
c 10.06.1996	ggu	$$UVPADV - modifications for advective term
c 23.07.1997	ggu	(from scratch) (hydro_transports_final)
c 14.08.1998	ggu	set w = 0 at open boundary nodes
c 20.08.1998	ggu	some documentation for sp256w
c 08.04.1999	ggu	equilibrium tide introduced (zeqv)
c 20.04.1999	ggu	converted to stress instead of wind (tauxnv...)
c 24.06.1999	ggu	call to rescur commented (use 2D call resid)
c 07.03.2000	ggu	eliminated VERT_AUST_ADJUST
c 07.03.2000	ggu	arrays for vertical turbulent diffusion coefficient
c 12.01.2001	ggu	solve for znv and not z difference (ZNEW) (hydro_zeta)
c 09.11.2001	ggu	BCHAO commented out, no compiler directives
c 14.11.2001	ggu	all compiler directives eliminated
c 10.08.2003	ggu	deleted commented parts (radiation, etc..)
c 10.08.2003	ggu	some code in new subroutines: make_prvel, copy_uvz
c 13.08.2003	ggu	delete useless vars (nrand, epsit), no iteration
c 13.08.2003	ggu	call setnod, set_link_info
c 10.03.2004	ggu	RQVDT - value in rqv is now discharge [m**3/s]
c 10.01.2005	ggu	BAROC_DIST - scale baroclinic term with distance
c 25.01.2005	ggu	BUGADV - bugfix for advective terms
c 24.02.2005	ggu	new reynolds stresses -> green
c 15.03.2005	ggu	austv,aust eliminated, austau() in diff_h_set()
c 29.04.2005	ggu	semi-lagrangian advection (backadv)
c 29.04.2005	ggu	no baroclinic terms close to step (ilevmin)
c 04.11.2005	ggu	use itlin to decide about semi-lagrangian advection
c 23.03.2006	ggu	changed time step to real
c 31.05.2006	ggu	new friction type ireib=7
c 18.10.2006	ccf	radx,rady for radiation stress introduced
c 28.11.2006	ggu	in u/vadv is now difference and not absolute transport
c 02.04.2007	ggu	new algorithm (look for ASYM)
c 08.06.2007	ggu&dbf	restructured for new explicit terms
c 28.09.2007	ggu	deleted chao, semi-lagrange to newexp.f, no indov
c 24.06.2008	ggu	bpresv deleted
c 10.10.2008	ggu&mbj	prepared for pardiso -> modify system (gguexclude)
c 10.12.2008	ggu	use rfricv for bottom friction -> other can be deleted
c 18.12.2008	ggu	more debug info for error in vertical system
c 13.01.2009	ggu	Pardiso lib integrated
c 04.03.2009	ggu	matrix amat deleted from file -> only locally used
c 27.03.2009	ggu	call new routine adjust_mass_flux() for dry nodes
c 06.04.2009	ggu	deleted routine write_elem_vel_info()
c 07.05.2009	ggu	new routines for scalar interpolation (not finished)
c 10.03.2010	ggu	bug fix in sp256w() for ibtyp=2
c 11.03.2010	ggu	new routine check_volume() to check for negative vol
c 23.03.2010	ggu	changed v6.1.1
c 12.04.2010	ggu	ad hoc routine for Yaron
c 22.04.2010	ggu	changed VERS_6_1_5
c 08.10.2010	ggu	changed VERS_6_1_13
c 15.12.2010	ggu	changed VERS_6_1_14
c 16.12.2010	ggu	in sp256w() account for changing volume (sigma)
c 19.02.2011	ccf	3D radiation stress
c 01.03.2011	ggu	changed VERS_6_1_20
c 31.05.2011	ggu	changed VERS_6_1_23
c 18.10.2011	ggu	changed VERS_6_1_33
c 04.11.2011	ggu	deleted computation of firction term (in subn35.f)
c 29.03.2012	ggu	cleaned up, sp256v prepared for OpenMP
c 10.05.2013	dbf&ggu	new routines for non-hydro
c 13.06.2013	ggu	changed VERS_6_1_65
c 12.09.2013	ggu	changed VERS_6_1_67
c 25.10.2013	ggu	changed VERS_6_1_68
c 29.10.2013	ggu	nudging implemented
c 12.11.2013	ggu	changed VERS_6_1_69
c 29.11.2013	ggu	zeta correction
c 05.12.2013	ggu	changed VERS_6_1_70
c 25.03.2014	ggu	new offline
c 10.04.2014	ggu	cleaning up of a lot of stuff
c 05.05.2014	ggu	changed VERS_6_1_74
c 18.06.2014	ggu	changed VERS_6_1_77
c 05.11.2014	ggu	changed VERS_7_0_5
c 05.12.2014	ggu	changed VERS_7_0_8
c 12.12.2014	ggu	changed VERS_7_0_9
c 19.12.2014	ggu	changed VERS_7_0_10
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.05.2015	ggu	changed VERS_7_1_10
c 06.05.2015	ggu	cleaning up of sp256f
c 20.05.2015	ggu&erp	sp256v parallelized
c 05.06.2015	ggu	changed VERS_7_1_12
c 10.07.2015	ggu	changed VERS_7_1_50
c 13.07.2015	ggu	changed VERS_7_1_51
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 24.07.2015	ggu	changed VERS_7_1_82
c 17.09.2015	ggu	sp256w renamed to hydro_vertical
c 17.09.2015	ggu	sp259f renamed to hydro
c 18.09.2015	ggu	sp256 renamed to hydro_transports, file cleaned
c 10.10.2015	ggu	changed VERS_7_3_2
c 22.10.2015	ggu	changed VERS_7_3_7
c 05.11.2015	ggu	changed VERS_7_3_12
c 20.11.2015	ggu&erp	chunk size introduced, omp finalized
c 16.12.2015	ggu	changed VERS_7_3_16
c 10.03.2016	ggu	in sp256v_intern() b/cpres in double precision
c 11.03.2016	ggu	most variables passed in double precision
c 14.06.2016	ggu	changed VERS_7_5_14
c 17.06.2016	ggu	changed VERS_7_5_15
c 31.10.2016	ggu	parallel part modified
c 12.01.2017	ggu	changed VERS_7_5_21
c 24.03.2017	ggu	new treatment for afix in sp256v_intern()
c 31.03.2017	ggu	changed VERS_7_5_24
c 04.11.2017	ggu	changed VERS_7_5_34
c 05.12.2017	ggu	changed VERS_7_5_39
c 24.01.2018	ggu	changed VERS_7_5_41
c 03.04.2018	ggu	changed VERS_7_5_43
c 19.04.2018	ggu	changed VERS_7_5_45
c 26.04.2018	ggu	changed VERS_7_5_46
c 11.05.2018	ggu	changed VERS_7_5_47
c 18.12.2018	ggu	changed VERS_7_5_52
c 19.12.2018	ggu	error in experimental code with penta solver
c 18.01.2019	ggu	finished implementing and testing penta solver
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 16.04.2019	ggu	introduced rcomputev for excluding elements (rcomp)
c 21.05.2019	ggu	changed VERS_7_5_62
c 02.07.2019	ggu	switched completely to penta solver
c 04.07.2019	ccf	for offline also compute horizontal diffusion params
c 16.07.2019	ggu	rmsdiff was not set to 0 (bug)
c 26.03.2020	ggu	adjust viscosity in case of closure (rcomp)
c 26.05.2020	ggu	new variable ruseterm to shut off selected terms
c 04.06.2020	ggu	debug_new3di() for selected debug
c 30.03.2021	ggu	copy to old into shyfem routine
c 30.04.2021    clr&ggu adapted for petsc solver
c 31.05.2021    ggu	eleminated bug for closing in hydro_transports_final()
c
c******************************************************************

	subroutine hydro

c administrates one hydrodynamic time step for system to solve

	use mod_depth
	use mod_bound_dynamic
	use mod_area
	use mod_hydro_baro
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv
	!use basin, only : nkn,nel,ngr,mbw
	use basin
	use shympi
        use mod_zeta_system, only : solver_type

	implicit none

	logical boff,bdebout
	logical bzcorr
	integer i,l,k,ie,ier,ii
	integer nrand
	integer iw,iwa,iloop
	integer nmat
	integer kspecial
	integer iwhat
	real azpar,ampar
	real dzeta(nkn)
	double precision dtime

	integer iround
	real getpar
	logical bnohyd

        real epseps
        parameter (epseps = 1.e-6)

        integer iwvel !DWNH
	kspecial = 0
	bdebout = .false.

	if( nint(getpar('ihydro')) <= 0 ) return	!only for debug

c-----------------------------------------------------------------
c set parameter for hydro or non hydro 
c-----------------------------------------------------------------

	call nonhydro_get_flag(bnohyd)
        iwvel = nint(getpar('iwvel')) !DWNH
	call get_act_dtime(dtime)

	azpar = getpar('azpar')
	ampar = getpar('ampar')
	if( azpar == 0. .and. ampar == 1. ) then
	  call system_set_explicit
	else if( azpar == 1. .and. ampar == 0. ) then
	  call system_set_explicit
	else if( shympi_is_parallel() ) then
	  !if( shympi_is_master() ) then
	  !  write(6,*) 'system is not solved explicitly'
	  !  write(6,*) 'cannot solve semi-implicitly with MPI'
	  !  write(6,*) 'azpar,ampar: ',azpar,ampar
	  !  write(6,*) 'please use azpar=1 and ampar=0'
	  !end if
	  !call shympi_stop('no semi-implicit solution')
	end if

c-----------------------------------------------------------------
c dry areas
c-----------------------------------------------------------------

	iw=0
	call sp136(iw)

c-----------------------------------------------------------------
c copy variables to old time level
c-----------------------------------------------------------------

	!call copy_uvz		!copies uvz to old time level
	!call nonhydro_copy	!copies non hydrostatic pressure terms
	!call copy_depth	!copies layer depth to old

	call set_diffusivity	!horizontal viscosity/diffusivity (needs uvprv)

c-----------------------------------------------------------------
c offline
c-----------------------------------------------------------------

	call is_offline(1,boff)
	!if( boff ) write(6,*) 'hydro reading from offline...'
	if( boff ) return

c-----------------------------------------------------------------
c solve for hydrodynamic variables
c-----------------------------------------------------------------

	iloop = 0

	do 				!loop over changing domain

	  iloop = iloop + 1

	  call hydro_transports		!compute intermediate transports

	  call setnod			!set info on dry nodes
	  call set_link_info		!information on areas, islands, etc..
	  call adjust_mass_flux		!cope with dry nodes

	  call system_init		!initializes matrix
	  call hydro_zeta(rqv)		!assemble system matrix for z
	  call system_solve(nkn,znv)	!solves system matrix for z
	  call system_get(nkn,znv)	!copies solution to new z

	  if(trim(solver_type) /= 'PETSc')then
            call shympi_exchange_2d_node(znv)
          endif

	  call setweg(1,iw)		!controll intertidal flats
	  !write(6,*) 'hydro: iw = ',iw,iloop,my_id
	  iw = shympi_sum(iw)
	  !if( iw > 0 .and. shympi_is_parallel() ) goto 99
	  if( iw == 0 ) exit

	end do

	call hydro_transports_final	!final transports (also barotropic)

c-----------------------------------------------------------------
c end of soulution for hydrodynamic variables
c-----------------------------------------------------------------

        call setzev			!copy znv to zenv
        call setuvd			!set velocities in dry areas
	call baro2l 			!sets transports in dry areas

	call make_new_depth
	call check_volume		!checks for negative volume 
        call arper

c-----------------------------------------------------------------
c vertical velocities and non-hydrostatic step
c-----------------------------------------------------------------

	if (bnohyd) then
	  call sp256wnh
	  call nonhydro_adjust
	end if

	call hydro_vertical(dzeta)		!compute vertical velocities

	!if (bnohyd) call nh_handle_output(dtime)
	if (bnohyd .or. (iwvel .eq. 1)) then
	  call nh_handle_output(dtime)!DWNH
	end if

c-----------------------------------------------------------------
c correction for zeta
c-----------------------------------------------------------------

	bzcorr = .true.
	bzcorr = .false.
	if( bzcorr ) then
	  call correct_zeta(dzeta)
          call setzev     !znv -> zenv
	  call make_new_depth
	  call hydro_vertical(dzeta)		!$$VERVEL
	end if

c-----------------------------------------------------------------
c some checks
c-----------------------------------------------------------------

	call vol_mass(1)		!computes and writes total volume
	if( bdebout ) call debug_output(dtime)
	call mass_conserve		!check mass balance

c-----------------------------------------------------------------
c compute velocities on elements and nodes
c-----------------------------------------------------------------

	call compute_velocities

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
   99	continue
	write(6,*) 'no wetting and drying for mpi yet...'
	stop 'error stop hydro: not yet ready'
	end

c******************************************************************

	subroutine hydro_zeta(vqv)

c assembles linear system matrix
c
c vqv		flux boundary condition vector
c
c semi-implicit scheme for 3d model

	use mod_nudging
	use mod_internal
	use mod_geom_dynamic
	use mod_depth
	use mod_bound_dynamic
	use mod_hydro_baro
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi
        use mod_zeta_system, only : kn,hia,hik,solver_type
         
	implicit none

	real vqv(nkn)

	double precision drittl
	parameter (drittl=1./3.)

	include 'mkonst.h'
	include 'pkonst.h'

        integer afix            !chao dbf
	logical bcolin
	logical bdebug
	integer ie,i,j,j1,j2,n,m,kk,l,k,ie_mpi
	integer ngl
	integer ilevel
	integer ju,jv
        integer nel_loop

	real azpar,ampar
	real dt

	!real az,am,af
	!real zm
	!real ht
	!real amatr(3,3)
	!real delta,h11,hh999
	!real z(3)
	!real andg,zndg(3)
	!real b(3),c(3)
	!real acu
	!real uold,vold
	!real ut,vt,uhat,vhat
	!real dbb,dbc,dcb,dcc,abn,acn

	double precision aj,rw,ddt
	double precision amatr(3,3)
	double precision delta,h11,hh999
	double precision z(3)
	double precision andg,zndg(3)
	double precision zm
	double precision az,am,af
	double precision b(3),c(3)
	double precision acu
	double precision uold,vold
	double precision ut,vt,uhat,vhat
	double precision dbb,dbc,dcb,dcc,abn,acn

c	data amatr / 2.,1.,1.,1.,2.,1.,1.,1.,2. /	!original
	data amatr / 4.,0.,0.,0.,4.,0.,0.,0.,4. /	!lumped

        integer locsps,loclp,iround
	real getpar
	!logical iskbnd,iskout,iseout
	logical iskbnd,iseout
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
        !iskout(k) = inodv(k).eq.-2
        iseout(ie) = iwegv(ie).ne.0

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	bcolin=iround(getpar('iclin')).ne.0

	call getazam(azpar,ampar)
	az=azpar
	am=ampar
	af=getpar('afpar')
	call get_timestep(dt)
	ddt = dt

	ngl=nkn

c-------------------------------------------------------------
c loop over elements
c-------------------------------------------------------------
        if(trim(solver_type)=='PETSc')then
          nel_loop=nel_unique
        else
          nel_loop=nel
        endif

	do ie_mpi=1,nel_loop
          if(trim(solver_type)=='PETSc')then
            ie = ie_mpi
          else
            ie = ip_sort_elem(ie_mpi)
          endif
	!write(6,*) ie_mpi,ie,ipev(ie),nel

c	------------------------------------------------------
c	compute level gradient
c	------------------------------------------------------

	zm=0.
	do i=1,3
		kk=nen3v(i,ie)
		kn(i)=kk
		b(i)=ev(i+3,ie)
		c(i)=ev(i+6,ie)
		!z(i)=zov(kk)
		z(i)=zeov(i,ie)		!ZEONV
		zndg(i) = andgzv(kk)	!nudging
		zm=zm+z(i)
	end do

	zm=zm*drittl

	!if(bcolin) then
	!	ht=hev(ie)
	!else
	!	ht=hev(ie)+zm
	!end if

	ilevel=ilhv(ie)
	aj=ev(10,ie)
        afix=1-iuvfix(ie)      !chao dbf

        delta=ddt*ddt*az*am*grav*afix         !ASYM_OPSPLT        !chao dbf

c	------------------------------------------------------
c	compute contribution from H^x and H^y
c	------------------------------------------------------

	dbb = 0.
	dbc = 0.
	dcb = 0.
	dcc = 0.
	do l=1,ilevel			!ASYM_OPSPLT
	  jv=l+l
	  ju=jv-1
	  dbb = dbb + ddxv(ju,ie)
	  dbc = dbc + ddyv(ju,ie)
	  dcb = dcb + ddxv(jv,ie)
	  dcc = dcc + ddyv(jv,ie)
	end do

c	------------------------------------------------------
c	compute barotropic transport
c	------------------------------------------------------

	uold = 0.
	vold = 0.
	uhat = 0.
	vhat = 0.

	do l=1,ilevel
	  uold = uold + utlov(l,ie)
	  vold = vold + vtlov(l,ie)
	  uhat = uhat + utlnv(l,ie)
	  vhat = vhat + vtlnv(l,ie)
	end do

	ut = az * uhat + (1.-az) * uold
	vt = az * vhat + (1.-az) * vold

c	------------------------------------------------------
c	set element matrix and RHS
c	------------------------------------------------------

	do n=1,3
	  do m=1,3
	    abn = b(n) * ( b(m) * dbb + c(m) * dbc )
	    acn = c(n) * ( b(m) * dcb + c(m) * dcc )
	    h11 = delta*( abn + acn )			!ASYM_OPSPLT_CH
	    hia(n,m) = aj * (amatr(n,m) + 12.*h11)
	  end do
	  acu = hia(n,1)*z(1) + hia(n,2)*z(2) + hia(n,3)*z(3)
	  andg = 4.*aj*ddt*zndg(n)
	  !hia(n,n) = hia(n,n) + 4 * ddt * aj / tau
	  hik(n) = acu + andg + 12.*aj*ddt*( ut*b(n) + vt*c(n) )	!ZNEW
	end do

c	------------------------------------------------------
c	level boundary conditions
c	------------------------------------------------------

	do i=1,3
	  if(rzv(kn(i)).ne.flag) then
		!rw=rzv(kn(i))-zov(kn(i))	!this for dz
		rw=rzv(kn(i))			!FIXME !ZNEW (this for znew)
		j1=mod(i,3)+1
		j2=mod(i+1,3)+1
		hik(j1)=hik(j1)-rw*hia(j1,i)
		hik(j2)=hik(j2)-rw*hia(j2,i)
		hia(i,j1)=0.
		hia(i,j2)=0.
		hia(j1,i)=0.
		hia(j2,i)=0.
		!hia(i,i)=12.*aj
		hik(i)=rw*hia(i,i)
	  end if
	  !call handle_ship_boundary(it,i,k,hia,hik)
	end do

c	------------------------------------------------------
c	excluded areas
c	------------------------------------------------------

          if( iseout(ie) ) then	!ZEONV
            hh999=aj*12.
            do n=1,3
              do m=1,3
                hia(n,m)=hh999*(b(n)*b(m)+c(n)*c(m))
              end do
              hik(n)=0.
            end do

            do n=1,3
              if( iskbnd(kn(n)) ) then	!not internal and not out of system
                do m=1,3
                  hia(n,m)=0.
                  !hia(m,n)=0.		!gguexclude - comment
                end do
                hik(n)=0.
              end if
            end do
          end if

c	------------------------------------------------------
c	in hia(i,j),hik(i),i,j=1,3 is system
c	------------------------------------------------------
	  !call system_assemble(ie,nkn,mbw,kn,hia,hik)
	  call system_assemble(ie)

	  !call debug_new3di('zeta',0,ie,hia,hik)

	end do

c-------------------------------------------------------------
c end of loop over elements
c-------------------------------------------------------------

c-------------------------------------------------------------
c Add additional flux boundary condition values to the rhs vector
c-------------------------------------------------------------

	call system_add_rhs(dt,nkn,vqv)

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c******************************************************************

	subroutine hydro_transports

	use basin, only : nkn,nel,ngr,mbw
!$	use omp_lib

	implicit none

	include 'pkonst.h'

	integer ie
	integer ies,iend
	integer ith
	integer count0,dcount,chunk,nt
	integer ibaroc
	integer ilin,itlin
	integer num_threads,myid,el_do,rest_do,init_do,end_do
	integer nchunk,nthreads
	logical bcolin,baroc
	real az,am,af,at,av,azpar,ampar
	real rlin,radv
	real vismol,rrho0
	real dt

	double precision rmsdif,rmsmax
	double precision tempo
	double precision openmp_get_wtime
	!integer openmp_get_num_threads,openmp_get_thread_num
	real getpar

c-------------------------------------------------------------
c initialize
c-------------------------------------------------------------

	ibaroc = nint(getpar('ibarcl'))		! baroclinic contributions
        vismol  = getpar('vismol')		! molecular viscosity
	bcolin = nint(getpar('iclin')).ne.0	! linearized conti
	itlin = nint(getpar('itlin'))		! advection scheme
	ilin = nint(getpar('ilin'))		! non-linear terms?
	rlin = getpar('rlin')			! non-linear strength?

	baroc = ibaroc .eq. 1 .or. ibaroc .eq. 2

	call getazam(azpar,ampar)
	az=azpar			! weighting in continuity
	am=ampar			! weighting in momentum
	af=getpar('afpar')		! weighting of coriolis term
	at=getpar('atpar')		! weighting of vertical viscosity
	av=getpar('avpar')		! weighting of advective terms

	radv = 0.
	if( ilin .eq. 0 .and. itlin .eq. 0 ) then	!need non-lin terms
	  radv = rlin * av	!strength * implicit factor
	end if

	call get_timestep(dt)
    
	rrho0=1./rowass
	if( .not. baroc ) rrho0 = 0.

c-------------------------------------------------------------
c computation of explicit part (sets arrays fxv(l,ie),fyv(l,ie)
c-------------------------------------------------------------

	call bottom_friction	!set bottom friction
        call set_explicit       !new HYDRO dbf
	!call set_yaron

c-------------------------------------------------------------
c loop over elements
c-------------------------------------------------------------

!$OMP PARALLEL 
!$OMP SINGLE

        call omp_compute_chunk(nel,nchunk)

	do ie=1,nel,nchunk

!$OMP TASK FIRSTPRIVATE(ie,bcolin,baroc,az,am,af,at,radv
!$OMP& 	   ,vismol,rrho0,dt) PRIVATE(ies,iend,rmsdif,rmsmax)
!$OMP&     SHARED(nel,nchunk)	 DEFAULT(NONE)
	 
	  rmsmax = 0.
 	  iend = ie+nchunk-1
 	  if(iend .gt. nel) iend = nel

 	  do ies=ie,iend
	    call hydro_intern(ies,bcolin,baroc,az,am,af,at,radv
     +			,vismol,rrho0,dt,rmsdif)
	    rmsmax = max(rmsmax,rmsdif)
	  end do

	  if( rmsmax > 1.D-10 ) then
	    write(6,*) 'rmsmax: ',rmsmax
	    stop 'error stop hydro_transports: rms too high'
	  end if

!$OMP END TASK

	end do

!$OMP END SINGLE
!$OMP TASKWAIT	
!$OMP END PARALLEL      

c-------------------------------------------------------------
c end of loop over elements
c-------------------------------------------------------------

ccc	call get_clock_count_diff(count0,dcount)
ccc	write(6,*) 'count: ',dcount

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c******************************************************************

	subroutine hydro_intern(ie,bcolin,baroc,az,am,af,at,radv
     +			,vismol,rrho0,dt,rmsdif)

c assembles vertical system matrix
c
c semi-implicit scheme for 3d model

	use tide
	use mod_meteo
	use mod_waves
	use mod_fluidmud
	use mod_internal
	use mod_depth
	use mod_layer_thickness
	use mod_roughness
	use mod_diff_visc_fric
	use mod_hydro_baro
	use mod_hydro_print
	use mod_hydro
	use evgeom
	use levels
	use basin
	!use ieee_exceptions

	implicit none

	integer ie
	logical bcolin,baroc
	real az,am,af,at
	real radv			!non-linear implicit contribution
	real vismol,rrho0
	real dt
	double precision rmsdif

c parameters
	double precision drittl
	parameter (drittl=1./3.)
c common
	include 'mkonst.h'
	include 'pkonst.h'
c local

	logical bbaroc,barea0                  !$$BAROC_AREA0
	logical bnewpenta

        integer afix             !chao dbf
	logical bfirst,blast
	logical debug,bdebug
        logical bdebggu
	integer kn(3)
	integer kk,ii,l,ju,jv,ierr
	integer ngl,mbb
	integer ilevel,ier,ilevmin
	integer lp,lm
	integer k1,k2,k3,k
        integer imin,imax
	real hlh
	real xbcl,ybcl
        real xexpl,yexpl
	real ulm,vlm
	real gamma,gammat
        real hhi,hhim,hhip,uui,uuim,uuip,vvi,vvim,vvip
	real bb,bbt,cc,cct,aa,aat,aux
	real rfric
	real aust
	real fact                       !$$BCHAO - not used
        real rhp,rhm,aus
	real hzg,gcz
        real xmin,xmax
        real xadv,yadv,fm,uc,vc,f,um,vm,up,vp
	real rraux,cdf,dtafix
	real ss
	logical b2d
	logical, parameter :: debug_mpi = .false.

	double precision b(3),c(3)
	double precision bpres,cpres,presx,presy
	double precision zz,zm,zmm
	double precision bz,cz
	double precision taux,tauy,rdist,rcomp,ruseterm
	double precision gravx,gravy,wavex,wavey
	double precision vis
	double precision uuadv,uvadv,vuadv,vvadv

c-----------------------------------------
	real hact(0:nlvdi+1)
	real rhact(0:nlvdi+1)
	real alev(0:nlvdi)
c-----------------------------------------
	double precision rmat(10*nlvdi)
	double precision smat(-2:2,2*nlvdi)
	double precision s2dmat(-1:1,2)		!for 2D
	double precision rvec(6*nlvdi)		!ASYM (3 systems to solve)
	double precision rvecp(6*nlvdi)		!ASYM (3 systems to solve)
	double precision solv(6*nlvdi)		!ASYM (3 systems to solve)
	double precision ppx,ppy
c-----------------------------------------
c function
	integer locssp,iround

        real epseps
        parameter (epseps = 1.e-6)

c-------------------------------------------------------------
c initialization and baroclinic terms
c-------------------------------------------------------------

	bnewpenta = .true.
	bdebug=.false.
	debug=.false.
        barea0 = .false.     ! baroclinic only with ia = 0 (HACK - do not use)

        bbaroc = baroc
	if( barea0 ) then               !$$BAROC_AREA $$BAROC_AREA0
	  if( iarv(ie) .ne. 0 ) bbaroc = .false.
        end if

c-------------------------------------------------------------
c dimensions of vertical system
c-------------------------------------------------------------

	ilevel=ilhv(ie)
	!ilevmin=ilmv(ie)
	ngl=2*ilevel
	mbb=2
	if(ngl.eq.2) mbb=1
	b2d = (ngl == 2)

c-------------------------------------------------------------
c compute barotropic terms (wind, atmospheric pressure, water level
c-------------------------------------------------------------

        rdist = rdistv(ie)		!use terms (distance from OB)
	rcomp = rcomputev(ie)		!use terms (custom elements)
        ruseterm = min(rcomp,rdist)	!use terms (both)

	bz=0.
	cz=0.
	bpres=0.
	cpres=0.
	zm=0.
	zmm=0.
	taux=0.
	tauy=0.
	do ii=1,3
	  kk=nen3v(ii,ie)
	  kn(ii)=kk
	  b(ii)=ev(ii+3,ie)
	  c(ii)=ev(ii+6,ie)

	  zz = zeov(ii,ie) - zeqv(kk)	!tide

          zm = zm + zz
	  zmm = zmm + zeov(ii,ie)		!ZEONV

	  bz=bz+zz*b(ii)
	  cz=cz+zz*c(ii)
	  bpres=bpres+ppv(kk)*b(ii)
	  cpres=cpres+ppv(kk)*c(ii)
	  taux=taux+tauxnv(kk)
	  tauy=tauy+tauynv(kk)
	end do

	zm=zm*drittl
	zmm=zmm*drittl
	taux=rcomp*taux*drittl
	tauy=rcomp*tauy*drittl

c-------------------------------------------------------------
c coriolis parameter
c-------------------------------------------------------------

	gammat=fcorv(ie)*ruseterm
        gamma=af*dt*gammat

c-------------------------------------------------------------
c reset vertical system 
c
c may be not the whole matrix every time
c ...size of matrix : ngl*(2*mbw+1) with mbw=2
c-------------------------------------------------------------

	do ii=1,ngl*5
	  rmat(ii)=0.
	end do
	if( b2d ) then
	  s2dmat = 0.
	else
	  smat(:,1:ngl) = 0.
	end if
	rvec = 0.
	solv = 0.

c-------------------------------------------------------------
c compute layer thicknes and store in hact and rhact
c-------------------------------------------------------------

	hact(0) = 0.
	do l=1,ilevel
	  hact(l) = hdeov(l,ie)
	end do
	hact(ilevel+1) = 0.
	hact(nlvdi+1) = 0.

	if( bcolin ) then
	  hact(1) = hact(1) - zmm		!FIXME
	end if

	do l=0,ilevel+1
	  if( hact(l) .le. 0. ) then
	    rhact(l) = 0.
	  else
	    rhact(l) = 1. / hact(l)
	  end if
	end do

c-------------------------------------------------------------
c compute element averaged turbulent viscosity
c-------------------------------------------------------------

	k1 = nen3v(1,ie)
	k2 = nen3v(2,ie)
	k3 = nen3v(3,ie)
	do l=0,ilevel
	    vis = vismol
	    vis = vis + (visv(l,k1)+visv(l,k2)+visv(l,k3))/3.
	    vis = vis + (vts(l,k1)+vts(l,k2)+vts(l,k3))/3.
	    if( rcomp /= 1. ) vis = vis_max * (1.-rcomp) + vis * rcomp
	    alev(l) = vis
	end do

c-------------------------------------------------------------
c start of vertical loop
c
c first set depth values
c
c hhi/hhip/hhim		thickness of i/i+1/i-1 layer
c uui/uuip/uuim		transport in x in i/i+1/i-1 layer
c vvi/vvip/vvim		transport in y in i/i+1/i-1 layer
c
c in case of a layer that does not exist (i-1 of first layer) give any
c ...value because the corrisponding a/b/c will be 0
c-------------------------------------------------------------

	do l=1,ilevel

	bfirst = l .eq. 1
	blast  = l .eq. ilevel
	
	lp = min(l+1,ilevel)
	lm = max(l-1,1)

	uui = utlov(l,ie)
	uuip = utlov(lp,ie)
	uuim = utlov(lm,ie)

	vvi = vtlov(l,ie)
	vvip = vtlov(lp,ie)
	vvim = vtlov(lm,ie)
        
	hhi = hact(l)
	hhip = hact(l+1)
	hhim = hact(l-1)

c	------------------------------------------------------
c	set up contributions of vertical viscosity
c	------------------------------------------------------

	rhp = 0.
	rhm = 0.
	if( hhi .gt. 0. ) then		!may be redundant
	  if( hhip .gt. 0. ) then	!lower interface
	    rhp = 2.0 * alev(l) / ( hhi + hhip )
	  end if
	  if( hhim .gt. 0. ) then	!upper interface
	    rhm = 2.0 * alev(l-1) / ( hhi + hhim )
	  end if
	end if
        
c	aus = afact * alev(l)
c	aux = dt * at * aus

	aus = 1.
	aux = dt * at

	aa  = aux * rhact(l) * ( rhm + rhp )
	aat = aus * rhact(l) * ( rhm + rhp )
	bb  = aux * rhact(l+1) * rhp
	bbt = aus * rhact(l+1) * rhp
	cc  = aux * rhact(l-1) * rhm
	cct = aus * rhact(l-1) * rhm

c	------------------------------------------------------
c	boundary conditions for stress on surface and bottom
c	------------------------------------------------------

	ppx = 0.
	ppy = 0.
	if( bfirst ) then
	  ppx = ppx - taux
	  ppy = ppy - tauy
	end if
	if( blast ) then
	  rfric = rfricv(ie)
	  if( rcomp /= 1. ) rfric = rfric_max * (1.-rcomp) + rfric * rcomp
	  aa  = aa + dt * rfric
	  aat = aat + rfric
	end if

c	------------------------------------------------------
c	implicit advective contribution
c	------------------------------------------------------

	uuadv = 0.
	uvadv = 0.
	vuadv = 0.
	vvadv = 0.

	aux = dt * radv * ruseterm	!implicit contribution

	if( aux .gt. 0. ) then		!implict treatment of non-linear terms

	uc = uui/hhi
	vc = vvi/hhi

	do ii=1,3
          k = nen3v(ii,ie)
          up = momentxv(l,k) / hhi
          vp = momentyv(l,k) / hhi
          f = uui * b(ii) + vvi * c(ii)
          if( f .lt. 0. ) then    !flux out of node => into element
	    uuadv = uuadv + aux*b(ii)*( up - uc )
	    uvadv = uvadv + aux*c(ii)*( up - uc )
	    vuadv = vuadv + aux*b(ii)*( vp - vc )
	    vvadv = vvadv + aux*c(ii)*( vp - vc )
            !xadv = xadv + f * ( up - uc )
            !yadv = yadv + f * ( vp - vc )
          end if
	end do

	end if

c	------------------------------------------------------
c	explicit contribution (non-linear, baroclinic, diffusion)
c	------------------------------------------------------
        
        xexpl = fxv(l,ie)			!explicit terms
        yexpl = fyv(l,ie)

	wavex = ruseterm * wavefx(l,ie)		!wave radiation stresss
	wavey = ruseterm * wavefy(l,ie)

	presx = rcomp * bpres			!atmospheric pressure
	presy = rcomp * cpres

	gravx = rcomp * grav*hhi*bz		!barotropic pressure
	gravy = rcomp * grav*hhi*cz

c	------------------------------------------------------
c	ppx/ppy is contribution on the left side of equation
c	ppx corresponds to -F^x_l in the documentation
c	ppy corresponds to -F^y_l in the documentation
c	------------------------------------------------------

	ppx = ppx + aat*uui - bbt*uuip - cct*uuim - gammat*vvi 
     +			+ gravx + (hhi/rowass)*presx + xexpl 
     +  		+ wavex
	ppy = ppy + aat*vvi - bbt*vvip - cct*vvim + gammat*uui 
     +			+ gravy + (hhi/rowass)*presy + yexpl 
     +  		+ wavey

c	------------------------------------------------------
c	set up matrix A
c	------------------------------------------------------

	jv=l+l
	ju=jv-1

	rmat(locssp(ju,ju,ngl,mbb)) = 1. + aa + uuadv
	rmat(locssp(jv,jv,ngl,mbb)) = 1. + aa + vvadv
	rmat(locssp(jv,ju,ngl,mbb)) =  gamma  + vuadv
	rmat(locssp(ju,jv,ngl,mbb)) = -gamma  + uvadv

	if( b2d ) then
	  s2dmat(0,ju) = 1. + aa + uuadv
	  s2dmat(0,jv) = 1. + aa + vvadv
	  s2dmat(-1,jv) =  gamma  + vuadv
	  s2dmat(+1,ju) = -gamma  + uvadv
	  !s2dmat(-1,jv) = -gamma  + vuadv
	  !s2dmat(+1,ju) =  gamma  + uvadv
	else
	  smat(0,ju) = 1. + aa + uuadv
	  smat(0,jv) = 1. + aa + vvadv
	  smat(-1,jv) =  gamma  + vuadv
	  smat(+1,ju) = -gamma  + uvadv
	end if

	if(.not.blast) then
		rmat(locssp(ju,ju+2,ngl,mbb)) = -bb
		rmat(locssp(jv,jv+2,ngl,mbb)) = -bb
		smat(+2,ju) = -bb
		smat(+2,jv) = -bb
        end if
	if(.not.bfirst) then
		rmat(locssp(ju,ju-2,ngl,mbb)) = -cc
		rmat(locssp(jv,jv-2,ngl,mbb)) = -cc
		smat(-2,ju) = -cc
		smat(-2,jv) = -cc
        end if

c	------------------------------------------------------
c	set up right hand side -F^x and -F^y 
c	------------------------------------------------------

	rvec(ju) = ppx
	rvec(jv) = ppy

c	------------------------------------------------------
c	set up H^x and H^y
c	------------------------------------------------------

	rvec(ngl+ju) = hhi		!ASYM_OPSPLT
	rvec(ngl+jv) = 0.d+0
	rvec(2*ngl+ju) = 0.d+0
	rvec(2*ngl+jv) = hhi

	end do

c-------------------------------------------------------------
c end of vertical loop
c-------------------------------------------------------------

c-------------------------------------------------------------
c solution of vertical system (we solve 3 systems in one call)
c-------------------------------------------------------------

c	----------------------------
c	new solution with penta
c	----------------------------

	rvecp = rvec
	if( b2d ) then
	  call tria_multi(ngl,3,s2dmat,rvecp,solv)
	else
	  call penta_fact(ngl,smat)
	  call penta_solve(ngl,smat,rvecp,solv)
	  call penta_solve(ngl,smat,rvecp(ngl+1),solv(ngl+1))
	  call penta_solve(ngl,smat,rvecp(2*ngl+1),solv(2*ngl+1))
	end if

c	----------------------------
c	old solution with dgelb
c	----------------------------

        !call gelb(rvec,rmat,ngl,1,mbb,mbb,epseps,ier)
        !call dgelb(rvec,rmat,ngl,1,mbb,mbb,epseps,ier)
        !call dgelb(rvec,rmat,ngl,3,mbb,mbb,epseps,ier)		!ASYM_OPSPLT

c	----------------------------
c	check difference
c	----------------------------

	rmsdif = 0.
	!rmsdif = sum((rvec-solv)**2)
	!rmsdif = sqrt(rmsdif/ngl)
	!if( rmsdif > 0.001 ) then
	!  write(6,*) 'rmsdif: ',ie,b2d,ngl,ilevel,rmsdif
	!end if

c	----------------------------
c	copy new solution
c	----------------------------

	rvec = solv

c	----------------------------
c	check for error
c	----------------------------

	!if(ier.ne.0) then
	!  call vel_matrix_error(ier,ie,ilevel,rvec,rmat,hact,alev)
	!  stop 'error stop hydro_intern: inverting vertical matrix'
	!end if

c-------------------------------------------------------------
c compute u^hat (negative sign because ppx/ppy was -F^x/-F^y)
c-------------------------------------------------------------

        afix=1-iuvfix(ie)       !chao dbf

	if( afix /= 0. ) then
	  dtafix = dt * afix
	  !ierr = ieee_handler( 'set', 'exception', SIGFPE_IGNORE )
	  do l=1,ilevel
	    utlnv(l,ie) = utlov(l,ie) - dtafix * rvec(2*l-1)
	    vtlnv(l,ie) = vtlov(l,ie) - dtafix * rvec(2*l)
	  end do
	  !ierr = ieee_handler( 'set', 'exception', SIGFPE_ABORT )
	end if

c-------------------------------------------------------------
c save contribution A^{-1} H^x and A^{-1} H^y
c-------------------------------------------------------------

	do l=1,ngl						!ASYM_OPSPLT
	  ddxv(l,ie) = rvec(ngl+l)
	  ddyv(l,ie) = rvec(2*ngl+l)
	end do

c-------------------------------------------------------------
c special information
c-------------------------------------------------------------

	if( ie .eq. 1 .and. barea0 .and. baroc ) then  !$$BAROC_AREA0
	  write(6,*) 'hydro_intern: BAROC_AREA0 active '
	end if

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c******************************************************************

	subroutine vel_matrix_error(ier,ie,lmax,rvec,rmat,hact,alev)

	implicit none

	integer ier,ie,lmax
	double precision rmat(10*lmax)
	double precision rvec(6*lmax)
	real hact(0:lmax+1)
	real alev(0:lmax)

	integer ii,l,kk,ngl,mbb
	real auxaux(-2:2)

	integer locssp

	ngl=2*lmax
	mbb=2
	if(ngl.eq.2) mbb=1

	write(6,*) 'Error in inverting matrix (vertical system)'
	write(6,*) 'ier : ',ier
	write(6,*) 'ie,lmax,ngl,mbb: ',ie,lmax,ngl,mbb
	write(6,*) 'rvec: ',(rvec(l),l=1,ngl)
	write(6,*) 'matrix: '

	do ii=1,ngl
	  do l=-2,2
	    kk=locssp(ii,ii+l,ngl,mbb)
	    if(kk.eq.0) then
		auxaux(l)=0.
	    else
		auxaux(l)=rmat(kk)
	    end if
	  end do
	  write(6,*) auxaux
	end do

	write(6,*) 'hact: ',(hact(l),l=0,lmax)
	write(6,*) 'alev: ',(alev(l),l=0,lmax)

	call check_set_unit(6)
	call check_elem(ie)
	call check_nodes_in_elem(ie)

	end

c******************************************************************

	subroutine hydro_transports_final

c post processing of time step

	use mod_internal
	use mod_depth
	use mod_hydro_baro
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi

	implicit none

	include 'mkonst.h'
	include 'pkonst.h'

	logical bcolin,bdebug
	integer ie,ii,l,kk,ie_mpi
	integer ilevel
	integer ju,jv
        integer afix            !chao dbf
	real dt,azpar,ampar
	double precision az,am,beta
	double precision bz,cz,um,vm
	double precision dz
	double precision du,dv
	double precision rfix,rcomp
c function
	integer iround
	real getpar

c-------------------------------------------------------------
c initialize
c-------------------------------------------------------------

	bcolin=iround(getpar('iclin')).ne.0	! linearized conti
	bdebug = .false.

	call get_timestep(dt)
	call getazam(azpar,ampar)
	az=azpar
	am=ampar

	beta = dt * grav * am 

c-------------------------------------------------------------
c start loop on elements
c-------------------------------------------------------------

	do ie_mpi=1,nel

	ie = ie_mpi
	ie = ip_sort_elem(ie_mpi)

	ilevel=ilhv(ie)

	rcomp = rcomputev(ie)		!use terms in element
        afix = 1 - iuvfix(ie)       	!chao dbf
	!rfix = afix * rcomp		!bug for closing
	rfix = afix

c	------------------------------------------------------
c	compute barotropic pressure term
c	------------------------------------------------------

	bz=0.
	cz=0.
	do ii=1,3
	  kk=nen3v(ii,ie)
	  dz = znv(kk) - zeov(ii,ie)
	  bz = bz + dz * ev(ii+3,ie)
	  cz = cz + dz * ev(ii+6,ie)
	end do

c	------------------------------------------------------
c	new transports from u/v hat variable
c	------------------------------------------------------

	do l=1,ilevel

	  jv=l+l
	  ju=jv-1

	  du = beta * ( ddxv(ju,ie)*bz + ddyv(ju,ie)*cz )	!ASYM_OPSPLT_CH
	  dv = beta * ( ddxv(jv,ie)*bz + ddyv(jv,ie)*cz )	!ASYM_OPSPLT_CH

	  utlnv(l,ie) = utlnv(l,ie) - du*rfix   		!chao dbf
	  vtlnv(l,ie) = vtlnv(l,ie) - dv*rfix   		!chao dbf

	end do

c	------------------------------------------------------
c	barotropic transports
c	------------------------------------------------------

	um = 0.
	vm = 0.
	do l=1,ilevel
	  um = um + utlnv(l,ie)
	  vm = vm + vtlnv(l,ie)
	end do
	unv(ie) = um
	vnv(ie) = vm

	end do

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c******************************************************************

	subroutine hydro_vertical(dzeta)

c computes vertical velocities
c
c velocities are computed on S/T points (top and bottom of layer)
c bottom velocity of the whole column is assumed to be 0
c -> maybe change this
c
c computes volume differences and from these computes vertical
c velocities at every time step so that the computed velocities
c satisfy the continuity equation for every single time step
c
c wlnv is computed horizontally at a node and vertically
c it is at the center of the layer -> there are nlv velocities
c computed
c
c b,c are 1/m, (phi is dimensionless)
c aj is m**2
c utlnv... is m**2/s
c dvol is in m**3/s
c vv is m**2 (area)
c
c wlnv is first used to accumulate volume difference -> dvol
c at the end it receives the vertical velocity
c
c wlnv (dvol)   aux array for volume difference
c vv            aux array for area

	use mod_bound_geom
	use mod_geom_dynamic
	use mod_bound_dynamic
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi

	implicit none

c arguments
	real dzeta(nkn)
c local
	logical debug
	integer k,ie,ii,kk,l,lmax
	integer ilevel
        integer ibc,ibtyp
	real aj,wbot,wdiv,ff,atop,abot,wfold
	real b,c
	real am,az,azt,azpar,ampar
	real ffn,ffo
	real volo,voln,dt,dvdt,q
	real dzmax,dz
	real, allocatable :: vf(:,:)
	real, allocatable :: va(:,:)
c statement functions

	logical is_zeta_bound
	real volnode

	!logical isein
        !isein(ie) = iwegv(ie).eq.0

c 2d -> nothing to be done

	dzeta = 0.
	wlnv = 0.
	if( nlvdi == 1 ) return

c initialize

	call getazam(azpar,ampar)
	az=azpar
	am=ampar
	azt = 1. - az
	call get_timestep(dt)

	allocate(vf(nlvdi,nkn),va(nlvdi,nkn))
	vf = 0.
	va = 0.

c compute difference of velocities for each layer
c
c f(ii) > 0 ==> flux into node ii
c aj * ff -> [m**3/s]     ( ff -> [m/s]   aj -> [m**2]    b,c -> [1/m] )

	do ie=1,nel
	 !if( isein(ie) ) then		!FIXME
	  aj=4.*ev(10,ie)		!area of triangle / 3
	  ilevel = ilhv(ie)
	  do l=1,ilevel
	    do ii=1,3
		kk=nen3v(ii,ie)
		b = ev(ii+3,ie)
		c = ev(ii+6,ie)
		ffn = utlnv(l,ie)*b + vtlnv(l,ie)*c
		ffo = utlov(l,ie)*b + vtlov(l,ie)*c
		ff = ffn * az + ffo * azt
		!ff = ffn
		vf(l,kk) = vf(l,kk) + 3. * aj * ff
		va(l,kk) = va(l,kk) + aj
	    end do
	  end do
	 !end if
	end do

	if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange vf,va')
          call shympi_exchange_and_sum_3d_nodes(vf)
          call shympi_exchange_and_sum_3d_nodes(va)
	end if

c from vel difference get absolute velocity (w_bottom = 0)
c	-> wlnv(nlv,k) is already in place !
c	-> wlnv(nlv,k) = 0 + wlnv(nlv,k)
c w of bottom of last layer must be 0 ! -> shift everything up
c wlnv(nlv,k) is always 0
c
c dividing wlnv [m**3/s] by area [vv] gives vertical velocity
c
c in va(l,k) is the area of the upper interface: a(l) = a_i(l-1)
c =>  w(l-1) = flux(l-1) / a_i(l-1)  =>  w(l-1) = flux(l-1) / a(l)

	dzmax = 0.

	do k=1,nkn
	  lmax = ilhkv(k)
	  wlnv(lmax,k) = 0.
	  debug = k .eq. 0
	  abot = 0.
	  do l=lmax,1,-1
	    atop = va(l,k)
            voln = volnode(l,k,+1)
            volo = volnode(l,k,-1)
	    dvdt = (voln-volo)/dt
	    q = mfluxv(l,k)
	    wdiv = vf(l,k) + q
	    !wfold = azt * (atop*wlov(l-1,k)-abot*wlov(l,k))
	    !wlnv(l-1,k) = wlnv(l,k) + (wdiv-dvdt+wfold)/az
	    wlnv(l-1,k) = wlnv(l,k) + wdiv - dvdt
	    abot = atop
	    if( debug ) write(6,*) k,l,wdiv,wlnv(l,k),wlnv(l-1,k)
	  end do
	  dz = dt * wlnv(0,k) / va(1,k)
	  dzmax = max(dzmax,abs(dz))
	  wlnv(0,k) = 0.	! ensure no flux across surface - is very small
	  dzeta(k) = dz
	end do

	!write(6,*) 'hydro_vertical: dzmax = ',dzmax

	do k=1,nkn
	  lmax = ilhkv(k)
	  debug = k .eq. 0
	  do l=2,lmax
	    atop = va(l,k)
	    if( atop .gt. 0. ) then
	      wlnv(l-1,k) = wlnv(l-1,k) / atop
	      if( debug ) write(6,*) k,l,atop,wlnv(l-1,k)
	    end if
	  end do
	end do

c set w to zero at open boundary nodes (new 14.08.1998)
c
c FIXME	-> only for ibtyp = 1,2 !!!!

	do k=1,nkn
          !if( is_external_boundary(k) ) then	!bug fix 10.03.2010
          if( is_zeta_bound(k) ) then
	    wlnv(:,k) = 0.
	    dzeta(k) = 0.
          end if
	end do

	deallocate(vf,va)

	if( shympi_partition_on_nodes() ) then
	  !call shympi_comment('exchanging wlnv')
          call shympi_exchange_3d0_node(wlnv)
	  !call shympi_barrier
	end if

	end

c*******************************************************************

	subroutine correct_zeta(dzeta)

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real dzeta(nkn)		!zeta correction

	znv = znv + dzeta

	end

c*******************************************************************

	subroutine debug_new3di(text,k,ie,hia,hik)

	use shympi
	use basin
	use mod_hydro
	use mod_hydro_baro
	use mod_system

	implicit none

	character*(*) text
	integer k,ie
	integer iunit
	real hia(3,3)
	real hik(3)
	type(smatrix), pointer :: mm

	logical bggu
	integer i,kk,ke,ike
	integer kn(3)
	double precision dtime

	integer ipext,ieext

	return
	call get_act_dtime(dtime)
	!if( dtime < 1038. ) return
	if( dtime < 1034. ) return
	if( dtime > 1042. ) return

	bggu = .false.
	if( ie > 0 ) then
	  do i=1,3
	    kk=nen3v(i,ie)
	    kn(i) = kk
	    ke = ipext(kk)
	    if( ke == 1934 ) then
	      bggu = .true.
	      ike = i
	    end if
	  end do
	  ke = 0
	else
	  !ke = ipext(k)
	  !if( ke == 1934 ) bggu = .true.
	end if

	if( .not. bggu ) return

	iunit = 166 + my_id

	mm => l_matrix

	write(iunit,*) '-----------------------'
	write(iunit,*) trim(text),dtime
	write(iunit,*) my_id,ke,k,ie,ieext(ie)

	if( ie > 0 ) then
	  write(iunit,*) nen3v(:,ie)
	  write(iunit,*) (ipext(nen3v(i,ie)),i=1,3)
	  write(iunit,*) zeov(:,ie)
	  write(iunit,*) zenv(:,ie)
	  write(iunit,*) unv(ie),vnv(ie)
	  !write(iunit,*) (hia(i,i),i=1,3)
	  !write(iunit,*) hik
	  write(iunit,*) hia(ike,ike)
	  write(iunit,*) hik(ike)
	  kk = kn(ike)
	  write(iunit,*) kk,ipext(kk),mm%raux2d(kk),mm%rvec2d(kk)
	end if

	if( k > 0 ) then
	  write(iunit,*) zov(k),znv(k)
	end if

	write(iunit,*) '-----------------------'
	  
	end

c*******************************************************************

