
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003-2004,2006,2008,2010-2017,2019  Georg Umgiesser
!    Copyright (C) 2007-2008,2015  Christian Ferrarin
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

c gotm module
c
c contents :
c
c revision log :
c
c 10.08.2003	ggu	call gotm_init
c 05.10.2004	ggu	new administration routine turb_closure, Munk-And.
c 23.03.2006	ggu	changed time step to real
c 20.11.2006	ggu	new version of keps, for gotm changed common blocks
c 20.10.2007	ccf	new version of gotm (4.0.0)
c 10.04.2008	ggu	integrated in main branch
c 18.09.2008	ccf	bug fix for m2 in setm2n2
c 02.12.2008	ggu	bug in gotm_init: no limiting values for initialization
c 18.12.2008	ggu	bug in GOTM module and setm2n2() corrected
c 23.03.2010	ggu	changed v6.1.1
c 16.02.2011	ggu	write n2max to info file, profiles in special node
c 30.03.2012	ggu	changed VERS_6_1_51
c 29.03.2013	ggu	avoid call to areaele -> ev(10,ie)
c 13.06.2013	ggu	changed VERS_6_1_65
c 25.03.2014	ggu	new offline
c 18.06.2014	ggu	changed VERS_6_1_77
c 27.06.2014	ggu	changed VERS_6_1_78
c 05.11.2014	ggu	changed VERS_7_0_5
c 05.12.2014	ggu	changed VERS_7_0_8
c 19.12.2014	ggu	changed VERS_7_0_10
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 10.07.2015	ggu	changed VERS_7_1_50
c 13.07.2015	ggu	changed VERS_7_1_51
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 24.07.2015	ggu	changed VERS_7_1_82
c 30.07.2015	ggu	changed VERS_7_1_83
c 03.12.2015	ccf	levdbg introduced for checka
c 16.12.2015	ggu	changed VERS_7_3_16
c 18.12.2015	ggu	changed VERS_7_3_17
c 25.05.2016	ggu	changed VERS_7_5_10
c 05.12.2017	ggu	changed VERS_7_5_39
c 03.02.2019	ggu	in gotm_shell check for 0layer and z0s/bmin (GGUZ0)
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 01.04.2021	ggu	new routine handle_gotm_init()
c 01.04.2021	ggu	debug code - look for iudbg
c 29.03.2022	ggu	do not need gotm support if running 2d simulation
c 07.04.2022	ggu	ie_mpi for bnstress, exchange visv,difv, debug code
c 09.04.2022	ggu	lots of debug code added
c 02.05.2023    ggu     fix mpi bug for nlv==1
c 09.05.2023    lrp     introduce top layer index variable
c 06.06.2023    ggu     minor change writing n2max
c
c**************************************************************

	subroutine turb_closure

c administers turbulence closure

	implicit none

	real getpar
	logical boff

	integer iturb
	save iturb
	data iturb / 0 /

	if( iturb .lt. 0 ) return

	call is_offline(4,boff)
	if( boff ) return

	if( iturb .eq. 0 ) then
	  iturb = nint(getpar('iturb'))
	  if( iturb .le. 0 ) iturb = -1
	  if( iturb .lt. 0 ) return
	  write(*,*) 'starting turbulence model: iturb = ',iturb
	end if

	if( iturb .eq. 1 ) then		!Gotm
	  call gotm_shell
	else if( iturb .eq. 2 ) then	!keps
	  call keps_shell
	else if( iturb .eq. 3 ) then	!Munk Anderson
	  call munk_anderson_shell
	else
	  write(6,*) 'Value iturb not possible: ',iturb
	  stop 'error stop turb_closure: iturb'
	end if

	end

c**************************************************************

	subroutine munk_anderson_shell

c computes turbulent quantities with Munk - Anderson model

	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'femtime.h'
	include 'pkonst.h'

c---------------------------------------------------------------
c aux arrays superposed onto other aux arrays
c---------------------------------------------------------------

	real shearf2(nlvdi,nkn)
	real buoyf2(nlvdi,nkn)
	real richard(nlvdi,nkn)
	real h(nlvdi)


	integer k,l
	integer nlev,flev
	integer mode
	real ri,vis,dif
	real diftur,vistur
	real a,b,alpha,beta

	real getpar

	integer icall
	save icall
	data icall / 0 /

c------------------------------------------------------
c initialization
c------------------------------------------------------

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  write(*,*) 'starting Munk Anderson turbulence model'
	end if

	icall = icall + 1

c------------------------------------------------------
c set up buoyancy frequency and shear frequency
c------------------------------------------------------

	call setm2n2(nlvdi,buoyf2,shearf2)

c------------------------------------------------------
c set up parameters
c------------------------------------------------------

	mode = +1
	a = 10.
	b = 3.33
	alpha = -1./2.
	beta = -3./2.
	vistur = getpar('vistur')
	diftur = getpar('diftur')

c------------------------------------------------------
c compute richardson number for each water column
c------------------------------------------------------

	do k=1,nkn

	    nlev = nlvdi
	    call dep3dnod(k,mode,flev,nlev,h)

	    do l=flev,nlev-1
	      ri = buoyf2(l,k) / shearf2(l,k)
	      vis = vistur*(1.+a*ri)**alpha
	      dif = diftur*(1.+b*ri)**beta

	      visv(l,k) = vis
	      difv(l,k) = dif
	      richard(l,k) = ri
	    end do

            if( nlev .eq. flev ) then
	      visv(flev-1,k) = vistur
	      difv(flev-1,k) = diftur
	      visv(flev,k) = vistur
	      difv(flev,k) = diftur
	      richard(flev,k) = 0.
            else
	      visv(flev-1,k) = visv(flev,k)
	      difv(flev-1,k) = difv(flev,k)
	      visv(nlev,k) = visv(nlev-1,k)
	      difv(nlev,k) = difv(nlev-1,k)
	      richard(nlev,k) = richard(nlev-1,k)
            end if

	end do

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	end

c**************************************************************

	subroutine gotm_shell

c computes turbulent quantities with GOTM model

	use mod_meteo
	use mod_gotm_aux
	use mod_ts
	use mod_roughness
	use mod_diff_visc_fric
	use mod_hydro_print
	use levels, only : nlvdi,nlv
	use basin
	use shympi

	implicit none

	include 'femtime.h'
	include 'pkonst.h'

	double precision dt
	double precision u_taus,u_taub

	double precision hh(0:nlvdi)
	double precision nn(0:nlvdi), ss(0:nlvdi)

	double precision num(0:nlvdi), nuh(0:nlvdi)
	double precision ken(0:nlvdi), dis(0:nlvdi)
	double precision len(0:nlvdi)

	double precision num_old(0:nlvdi), nuh_old(0:nlvdi)
	double precision ken_old(0:nlvdi), dis_old(0:nlvdi)
	double precision len_old(0:nlvdi)
 
c---------------------------------------------------------------
c aux arrays superposed onto other aux arrays
c---------------------------------------------------------------

	real taub(nkn)

	integer iunit,iudbg,kdebug,iuout,iwhat
	integer ioutfreq,ks
	integer k,l
	integer laux
	integer nlev,flev,numOfLev
	real g
	real czdef
	save czdef

	real h(nlvdi)
	double precision depth		!total depth [m]
	double precision z0s,z0b	!surface/bottom roughness length [m]
	double precision rlmax,dz0
	integer nltot,iudeb
	logical bwrite

	double precision, parameter :: dz0min = 1.1	!min value for dz0=d/z0
	real, parameter :: charnock_val=1400.	!emp. Charnock constant

	real ubot,vbot,rr

	real dtreal
	real getpar
	integer ipext

	logical bwave,has_waves,bgotm,bdeb
	save bwave

	character*80 fn	
	integer, save :: icall = 0

	integer, save	:: levdbg

c------------------------------------------------------
c documentation
c------------------------------------------------------

c total stress :  tau_x = a |u| u   tau_y = a |u| v
c tau_tot^2 = a^2 |u|^2 ( u^2 + v^2 ) a^2 |u|^4
c tau_tot = sqrt(tau_x^2+tau_y^2) = a |u|^2
c
c u_taus=sqrt(sqrt(tx*tx+ty*ty))	(friction velocity)

c------------------------------------------------------
c initialization
c------------------------------------------------------

	if( icall .lt. 0 ) return

	kdebug = 8
	kdebug = -1
	iudeb = 0
	iwhat = 0
	if( kdebug > 0 ) iudeb = 760 + my_id

	if( icall .eq. 0 ) then

	  if( nlv_global <= 1 ) icall = -1	!avoid nlv==1 mpi bug
	  if( icall < 0 ) return

	  call has_gotm(bgotm)
	  if( .not. bgotm ) then
	    write(6,*) 'the model has been compiled without GOTM support'
	    write(6,*) 'please eneable GOTM=true in the Rules.make file'
	    write(6,*) 'or set iturb to another value'
	    stop 'error stop gotm_shell: no GOTM support'
	  end if

	  czdef = getpar('czdef')
	  bwave = has_waves()
          levdbg = nint(getpar('levdbg'))

c         --------------------------------------------------------
c         Initializes gotm arrays 
c         --------------------------------------------------------

	  write(*,*) 'starting initializing GOTM turbulence model'

	  call handle_gotm_init

c         --------------------------------------------------------
c         Internal GOTM turbulence model initialization
c         --------------------------------------------------------

          call getfnm('gotmpa',fn)
	  iunit = 10
	  call init_gotm_turb(iunit,fn,nlvdi)

	  write(*,*) 'finished initializing GOTM turbulence model'

	  icall = 1
	  call shympi_barrier
	end if

	call get_timestep(dtreal)
	dt = dtreal
	g = grav

c------------------------------------------------------
c set up bottom stress on nodes
c------------------------------------------------------

	call bnstress(czdef,taub)

	!call shympi_comment('exchanging taub')
	call shympi_exchange_2d_node(taub)

c------------------------------------------------------
c set up buoyancy frequency and shear frequency
c------------------------------------------------------

	shearf2 = 0.
 	buoyf2 = 0.
	call setm2n2(nlvdi,buoyf2,shearf2)

	!write(6,*) 'm2n2: ',maxval(buoyf2),maxval(shearf2)

c------------------------------------------------------
c call gotm for each water column
c------------------------------------------------------

	rlmax = 0.
	nltot = 0
	bdeb = .false.

	do k=1,nkn

	    bdeb = ( kdebug > 0 .and. ipext(k) == kdebug )
	    call save_gotm_set_debug(bdeb)

	    nlev = nlvdi
	    call dep3dnod(k,+1,flev,nlev,h)		!here nlev,flev is passed back
            numOfLev = nlev-flev+1

	    if( count( h(flev:nlev) <= 0. ) > 0 ) goto 97
            if( nlev .eq. flev ) goto 1

	    if( numOfLev < 1 .or. numOfLev > 10000 ) goto 96

c           ------------------------------------------------------
c           update boyancy and shear-frequency vectors
c           ------------------------------------------------------

	    do l=flev,nlev-1
	      laux = nlev - l
	      nn(laux) = buoyf2(l,k)
	      ss(laux) = shearf2(l,k)
	    end do
	    nn(0) = 0.
	    nn(numOfLev) = 0.
	    ss(0) = ss(1)
	    ss(numOfLev) = ss(numOfLev-1)

c           ------------------------------------------------------
c           compute layer thickness and total depth
c           ------------------------------------------------------

	    depth = 0.
	    do l=flev,nlev
	      hh(nlev-l+1) = h(l)
	      depth = depth + h(l)
	    end do
	    hh(0) = 0.

c           ------------------------------------------------------
c           compute surface friction velocity (m/s)
c           ------------------------------------------------------

	    u_taus = sqrt( sqrt( tauxnv(k)**2 + tauynv(k)**2 ) )

	    if ( bwave ) then
	      z0s = z0sk(k)
	    else
	      z0s = charnock_val*u_taus**2/g
	    end if
	    z0s = max(z0s,z0smin)			!GGUZ0

c           ------------------------------------------------------
c           compute bottom friction velocity (m/s)
c           ------------------------------------------------------

	    z0b = z0bk(k)
	    z0b = max(z0bmin,z0b)			!GGUZ0
	    u_taub = sqrt( taub(k) )

	    ubot = uprv(nlev,k)
  	    vbot = vprv(nlev,k)
	    dz0 = hh(1)/z0b
	    if( dz0 < dz0min ) dz0 = dz0min		!GGUZ0
	    rr = 0.4/( log( 0.5*(1.+dz0) ) )
	    !rr = 0.4/(log((z0b+hh(1)/2)/z0b))
            u_taub = rr*sqrt( ubot*ubot + vbot*vbot )

c           ------------------------------------------------------
c           update 1-dimensional vectors
c           ------------------------------------------------------

	    num(:) = numv_gotm(:,k)
	    nuh(:) = nuhv_gotm(:,k)
	    ken(:) = tken_gotm(:,k)
	    dis(:) = eps_gotm(:,k)
	    len(:) = rls_gotm(:,k)
	    num_old(:) = numv_gotm(:,k)
	    nuh_old(:) = nuhv_gotm(:,k)
	    ken_old(:) = tken_gotm(:,k)
	    dis_old(:) = eps_gotm(:,k)
	    len_old(:) = rls_gotm(:,k)

	    !call save_gotm_init

c           ------------------------------------------------------
c           call GOTM turbulence routine
c           ------------------------------------------------------

	    !write(6,*) 'before do_gotm_turb: nlev = ',k,nlev,nlvdi
	    if( bdeb .and. iudeb > 0 ) then
	      iwhat = iwhat + 1
	      iuout = iudeb
              call gotm_debug_write_f(iuout,iwhat,k,nlev
     +			,dt,depth
     +			,u_taus,u_taub,z0s,z0b
     +			,hh,nn,ss,num,nuh,ken,dis,len)
	      iuout = iudeb + 100
              call gotm_debug_write_u(iuout,iwhat,k,nlev
     +			,dt,depth
     +			,u_taus,u_taub,z0s,z0b
     +			,hh,nn,ss,num,nuh,ken,dis,len)
	    end if

 	    call do_gotm_turb   (
     &				  numOfLev,dt,depth
     &				 ,u_taus,u_taub
     &				 ,z0s,z0b,hh
     &	                         ,nn,ss
     &				 ,num,nuh,ken,dis,len
     &				 )

	    if( bdeb .and. iudeb > 0 ) then
	      iwhat = iwhat + 1
	      iuout = iudeb
              call gotm_debug_write_f(iuout,iwhat,k,nlev
     +			,dt,depth
     +			,u_taus,u_taub,z0s,z0b
     +			,hh,nn,ss,num,nuh,ken,dis,len)
	      iuout = iudeb + 100
              call gotm_debug_write_u(iuout,iwhat,k,nlev
     +			,dt,depth
     +			,u_taus,u_taub,z0s,z0b
     +			,hh,nn,ss,num,nuh,ken,dis,len)
	    end if

c           ------------------------------------------------------
c           copy back to node vectors
c           ------------------------------------------------------

	    numv_gotm(:,k) = num(:)
	    nuhv_gotm(:,k) = nuh(:)
	    tken_gotm(:,k) = ken(:)
	    eps_gotm(:,k)  = dis(:)
	    rls_gotm(:,k)  = len(:)

	    bwrite = .false.
	    do l=0,numOfLev
	      rlmax = max(rlmax,len(l))	!ggu
	      if( len(l) .gt. 100. ) then
		nltot = nltot + 1
		bwrite = .true.
	      end if

	    end do

	    iudbg = 489
	    iudbg = 0
	    bwrite = .false.
	    bwrite = ( iudbg > 0 .and. k == 1659 .and. it > 14400 )

	    if( bwrite ) then

	      !write(45,*) it,k
	      !call save_gotm_write

	      write(iudbg,*) '==========================='
	      write(iudbg,*) 'time___: ',it
	      write(iudbg,*) 'node___: ',k
	      write(iudbg,*) 'nlev___: ',nlev
	      write(iudbg,*) 'dt_____: ',dt
	      write(iudbg,*) 'depth__: ',depth
	      write(iudbg,*) 'utaus_b: ',u_taus,u_taub
	      write(iudbg,*) 'z0s_z0b: ',z0s,z0b
	      write(iudbg,*) 'hh_____: ',(hh     (l),l=0,nlev)
	      write(iudbg,*) 'nn_____: ',(nn     (l),l=0,nlev)
	      write(iudbg,*) 'ss_____: ',(ss     (l),l=0,nlev)
	      write(iudbg,*) 'num_old: ',(num_old(l),l=0,nlev)
	      write(iudbg,*) 'nuh_old: ',(nuh_old(l),l=0,nlev)
	      write(iudbg,*) 'ken_old: ',(ken_old(l),l=0,nlev)
	      write(iudbg,*) 'dis_old: ',(dis_old(l),l=0,nlev)
	      write(iudbg,*) 'len_old: ',(len_old(l),l=0,nlev)
	      write(iudbg,*) 'num____: ',(num    (l),l=0,nlev)
	      write(iudbg,*) 'nuh____: ',(nuh    (l),l=0,nlev)
	      write(iudbg,*) 'ken____: ',(ken    (l),l=0,nlev)
	      write(iudbg,*) 'dis____: ',(dis    (l),l=0,nlev)
	      write(iudbg,*) 'len____: ',(len    (l),l=0,nlev)
	      write(iudbg,*) '==========================='
	    end if

c           ------------------------------------------------------
c           update viscosity and diffusivity
c           ------------------------------------------------------

	    do l=flev-1,nlev
	      laux = nlev - l
	      visv(l,k) = num(laux)
	      difv(l,k) = nuh(laux)
	    end do

    1     continue
	end do

	call shympi_exchange_3d0_node(visv)
	call shympi_exchange_3d0_node(difv)

	iudbg = 654
	iudbg = 0
	if( iudbg > 0 ) then
	  k=2201
	  nlev = nlvdi
	  call dep3dnod(k,+1,flev,nlev,h)
	  write(iudbg,*) 'shell: ',it,k,nlev,numv_gotm(:,k)
	  !write(iudbg,*) 'depth: ',k,nlev,h(1:nlev)
	end if

	if( levdbg >= 2 ) call checka(nlvdi,shearf2,buoyf2,taub)
 
	ks = 0			!internal node number
	ioutfreq = 3600		!output frequency
	if( ks .gt. 0 .and. mod(it,ioutfreq) .eq. 0 ) then
	  write(188,*) it,nlev,(visv(l,ks),l=flev,nlev)
	  write(189,*) it,nlev,(difv(l,ks),l=flev,nlev)
	end if

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	return
   96	continue
	write(6,*) 'wrong number of layers: ',numOfLev
	stop 'error stop gotm_shell: wrong nlev'
   97	continue
	write(6,*) 'layers without depth...'
	write(6,*) k,nlev
	write(6,*) h(flev:nlev)
	stop 'error stop gotm_shell: no layer'
	end

c**************************************************************

	subroutine gotm_debug(iu,iwhat,k,nlev,dt,depth
     +			,u_taus,u_taub,z0s,z0b
     +			,hh,nn,ss,num,nuh,ken,dis,len)

	implicit none

	integer iu,iwhat,k,nlev
	double precision dt,depth
	double precision u_taus,u_taub
	double precision z0s,z0b
	double precision :: hh(0:nlev)
	double precision :: nn(0:nlev), ss(0:nlev)

	double precision :: num(0:nlev), nuh(0:nlev)
	double precision :: ken(0:nlev), dis(0:nlev)
	double precision :: len(0:nlev)

	integer iuu

	iuu = iu + 100

	      write(iu,*) k,iwhat,nlev
	      write(iu,*) dt,depth
	      write(iu,*) u_taus,u_taub
	      write(iu,*) z0s,z0b
	      write(iu,*) hh(0:nlev)
	      write(iu,*) nn(0:nlev)
	      write(iu,*) ss(0:nlev)
	      write(iu,*) num(0:nlev)
	      write(iu,*) nuh(0:nlev)
	      write(iu,*) ken(0:nlev)
	      write(iu,*) dis(0:nlev)
	      write(iu,*) len(0:nlev)
	      flush(iu)

	      write(iuu) k,iwhat,nlev
	      write(iuu) dt,depth
	      write(iuu) u_taus,u_taub
	      write(iuu) z0s,z0b
	      write(iuu) hh(0:nlev)
	      write(iuu) nn(0:nlev)
	      write(iuu) ss(0:nlev)
	      write(iuu) num(0:nlev)
	      write(iuu) nuh(0:nlev)
	      write(iuu) ken(0:nlev)
	      write(iuu) dis(0:nlev)
	      write(iuu) len(0:nlev)
	      flush(iuu)

	end

c**************************************************************

	subroutine gotm_internal_init

c initializes gotm arrays

	use mod_gotm_aux

	implicit none

        double precision, parameter    :: num_min  = 1.e-6
        double precision, parameter    :: nuh_min  = 1.e-6
        double precision, parameter    :: tken_min = 1.e-10
        double precision, parameter    :: eps_min  = 1.e-12
        double precision, parameter    :: rls_min  = 1.e-10

        numv_gotm = num_min
        nuhv_gotm = nuh_min
        tken_gotm = tken_min
        eps_gotm  = eps_min
        rls_gotm  = rls_min

	end

c**************************************************************

	subroutine handle_gotm_init

c handles initialization of gotm

	use basin
	use levels
	use mod_gotm_aux
	use mod_diff_visc_fric

	implicit none

	logical bdebug
	integer k,l,lmax,lmin,laux
	integer, save :: iudbg = 0
	integer, save :: icall = 0

	real h(10)
	logical rst_use_restart

	if( icall > 0 ) return
	icall = 1

	bdebug = ( iudbg > 0 )

	write(6,*) 'handle_gotm_init: ',rst_use_restart(8)

	if( bdebug ) then
	  write(iudbg,*) 'handle_gotm_init: ',rst_use_restart(8)
	  write(iudbg,*) 'init: ',numv_gotm(:,1)
	end if

	call mod_gotm_aux_init(nkn,nlvdi)	!probably useless

        if( .not. rst_use_restart(8) ) then
	  call gotm_internal_init
        end if

	visv = 0.
	difv = 0.

	do k=1,nkn
	  lmax = ilhkv(k)
	  lmin = jlhkv(k)
	  do l=lmin-1,lmax
	    laux = lmax - l
	    visv(l,k) = numv_gotm(laux,k)
	    difv(l,k) = nuhv_gotm(laux,k)
	  end do
	end do

	if( bdebug ) then
	!do k=1,nkn,200
	!write(iudbg,*) 'init: ',k,ilhkv(k),numv_gotm(:,k)
	!end do
	write(iudbg,*) '----------------'
	k=2201
	lmax = 10
	call dep3dnod(k,+1,lmin,lmax,h)
	write(iudbg,*) 'init: ',k,ilhkv(k),numv_gotm(:,k)
	write(iudbg,*) '----------------'
	end if

	end

c**************************************************************

	subroutine gotm_get(k,nlev,num,nuh,tk,ep,rl)

c returns internal parameters from turbulence closure

	use mod_gotm_aux

	implicit none

	integer k		!node
	integer nlev		!number of levels to return
	real num(0:nlev)	!viscosity
	real nuh(0:nlev)	!diffusivity
	real tk (0:nlev)	!kinetic energy
	real ep (0:nlev)	!dissipation
	real rl (0:nlev)	!length scale

	integer l,laux

        do l=0,nlev
	  laux = nlev - l
          num(l) = numv_gotm(laux,k)
          nuh(l) = nuhv_gotm(laux,k)
          tk(l)  = tken_gotm(laux,k)
          ep(l)  = eps_gotm (laux,k)
          rl(l)  = rls_gotm (laux,k)
        end do

	end

c**************************************************************

      subroutine Yevol(Nmx,dt,h,nuh,difmol,Qsour,Yvar,Yd,
     &                 Bcup,Yup,Bcdw,Ydw,Taur)

      implicit none

	integer Nmx		!number of levels
	real dt			!time step [s]
	real h(Nmx)		!level thickness [m]
	real nuh(0:Nmx)		!turbulent diffusivity [m**2/s]
	real difmol		!molecular diffusivity [m**2/s]
	real Qsour(Nmx)		!internal source term
	real Yvar(Nmx)		!variable to diffuse
	real Yd(Nmx)		!relaxation value for variable
	integer Bcup		!type of boundary condition (upper)
	real Yup		!value for boundary condition (upper)
	integer Bcdw		!type of boundary condition (lower)
	real Ydw		!value for boundary condition (lower)
	real Taur		!time scale to use with Yd [s]

c if Taur <= 0	-> do not use Yd
c Bcup/Bcdw :   1 = Neuman condition   2 = Dirichlet condition
c Yup,Ydw   :   value to use for boundary condition
c
c example: Bcup = 1  and Yup = 0  --> no flux condition across surface

!     This subroutine computes the evolution of any 
!     state variable "Y" after Mixing/Advection/Source
!     - with Neuman    Boundary Condition: Bc=1
!     - with Dirichlet Boundary Condition: Bc=2
!     Convention: fluxex are taken positive upward

	integer MaxNmx
	parameter(MaxNmx = 500)

	integer i
	double precision au(0:MaxNmx),bu(0:MaxNmx)
	double precision cu(0:MaxNmx),du(0:MaxNmx)
	double precision avh(0:MaxNmx),Y(0:MaxNmx)
	double precision cnpar,a,c

	if( Nmx .gt. MaxNmx ) stop 'error stop Tridiagonal: MaxNmx'

      cnpar = 1.      		!Crank Nicholson parameter (implicit)

!
! copy to auxiliary values
!

	 Y(0) = 0.
	 avh(0) = 0.
	 do i=1,Nmx
	   Y(i) = Yvar(i)
	   avh(i) = nuh(i) + difmol
	 end do

! Water column 
! => [a b c].X=[d] where X=[2, ...,i,i+1,...Nmx-1]

         do i=2,Nmx-1
          c    =2*dt*avh(i)  /(h(i)+h(i+1))/h(i) 
          a    =2*dt*avh(i-1)/(h(i)+h(i-1))/h(i)
          cu(i)=-cnpar*c                                 !i+1,n+1
          au(i)=-cnpar*a                                 !i-1,n+1
          bu(i)=1-au(i)-cu(i)                            !i  ,n+1
          du(i)=Y(i)+dt*Qsour(i)                         !i  ,n
     &          +(1-cnpar)*(a*Y(i-1)-(a+c)*Y(i)+c*Y(i+1))
         end do

! Surface 
! => [a b /].X=[d] where X=[Nmx-1,Nmx,/]

         if (Bcup.eq.1) then                       !BC Neuman
          a      =2*dt*avh(Nmx-1)/(h(Nmx)+h(Nmx-1))/h(Nmx)
          au(Nmx)=-cnpar*a
          bu(Nmx)=1-au(Nmx)
          du(Nmx)=Y(Nmx)
     &            +dt*(Qsour(Nmx)-Yup/h(Nmx))   
     &            +(1-cnpar)*a*(Y(Nmx-1)-Y(Nmx))
         else if (Bcup.eq.2) then                    !BC Dirichlet
          au(Nmx)=0.
          bu(Nmx)=1.
          du(Nmx)=Yup
         end if

! Bottom  
! => [/ b c].X=[d] where X=[/,1,2]

         if (Bcdw.eq.1) then                        !BC Neuman              
          c    =2*dt*avh(1)/(h(1)+h(2))/h(1)
          cu(1)=-cnpar*c 
          bu(1)=1-cu(1)
          du(1)=Y(1)             
     &          +dt*(Qsour(1)+Ydw/h(1))
     &          +(1-cnpar)*c*(Y(2)-Y(1))
         else if (Bcdw.eq.2) then                    !BC Dirichlet
          cu(1)=0.
          bu(1)=1.
          du(1)=Ydw
         end if
!
! Implicit internal restoring
!
         if ( Taur .gt. 0. .and. Taur .lt. 1.E10 ) then
          do i=1,Nmx
           bu(i)=bu(i)+dt/Taur
           du(i)=du(i)+dt/Taur*Yd(i)
          end do
         end if
!
! Implicit vertical mixing
!

         call Tridiagonal(Nmx,1,Nmx,au,bu,cu,du,Y)

	 do i=1,Nmx
	   Yvar(i) = Y(i)
	 end do

      return
      end

c**************************************************************

      subroutine Tridiagonal(Nmx,fi,lt,au,bu,cu,du,value)
      implicit none

	integer MaxNmx
	parameter(MaxNmx=500)
!
! !USES:
!
! !INPUT PARAMETERS:
      integer Nmx
      integer fi,lt
      double precision au(0:Nmx),bu(0:Nmx),cu(0:Nmx),du(0:Nmx)
!
! !OUTPUT PARAMETERS:
      double precision value(0:Nmx)
!
! !DESCRIPTION:
!     Used trough out the program.
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  25Nov98   The GOTM developers  Initial code.
!
! !LOCAL VARIABLES:
      double precision ru(0:MaxNmx),qu(0:MaxNmx)
      integer i
!
!EOP
!BOC
	if( Nmx .gt. MaxNmx ) stop 'error stop Tridiagonal: MaxNmx'

      ru(lt)=au(lt)/bu(lt)
      qu(lt)=du(lt)/bu(lt)

      do i=lt-1,fi+1,-1
         ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
         qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
      end do

      qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

      value(fi)=qu(fi)
      do i=fi+1,lt
         value(i)=qu(i)-ru(i)*value(i-1)
      end do

      return
      end
!EOC

c**************************************************************

	subroutine setm2n2(nldim,buoyf2,shearf2)

c sets buoyancy and shear frequency
c
c this is done directly for each node
c in rhov is already rho^prime = rho - rho0 (deviation)
!  Discretisation of vertical shear squared according to Burchard (2002)
!  in order to guarantee conservation of kinetic energy when transformed
!  from mean kinetic energy to turbulent kinetic energy.
c
c bug fix in computation of shearf2 -> abs() statements to avoid negative vals

	use mod_ts
	use mod_hydro_print
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

	integer nldim
	real buoyf2(nldim,nkn)
	real shearf2(nldim,nkn)

	include 'pkonst.h'
	include 'femtime.h'

	integer k,l,nlev,flev
	real aux,dh,du,dv,m2,dbuoy
	real h(nldim)
	real cnpar			!numerical "implicitness" parameter
	real n2max,n2
	real nfreq,nperiod

	integer iuinfo
	save iuinfo
	data iuinfo / 0 /
 
        aux = -grav / rowass
	cnpar = 1
	n2max = 0.
 
        do k=1,nkn
	  nlev = nldim
          call dep3dnod(k,+1,flev,nlev,h)
          do l=flev,nlev-1
            dh = 0.5 * ( h(l) + h(l+1) )
            dbuoy = aux * ( rhov(l,k) - rhov(l+1,k) )
            n2 = dbuoy / dh
	    n2max = max(n2max,n2)
            buoyf2(l,k) = n2

            du = 0.5*(
     &       (cnpar*abs((uprv(l+1,k)-uprv(l,k))
     &		*(uprv(l+1,k)-upro(l,k)))+
     &       (1.-cnpar)*abs((upro(l+1,k)-upro(l,k))
     &		*(upro(l+1,k)-uprv(l,k))))
     &       /dh/h(l)
     &      +(cnpar*abs((uprv(l+1,k)-uprv(l,k))
     &		*(upro(l+1,k)-uprv(l,k)))+
     &       (1.-cnpar)*abs((upro(l+1,k)-upro(l,k))
     &		*(uprv(l+1,k)-upro(l,k))))
     &       /dh/h(l+1)
     &                )

            dv = 0.5*(
     &       (cnpar*abs((vprv(l+1,k)-vprv(l,k))
     &		*(vprv(l+1,k)-vpro(l,k)))+
     &       (1.-cnpar)*abs((vpro(l+1,k)-vpro(l,k))
     &		*(vpro(l+1,k)-vprv(l,k))))
     &       /dh/h(l)
     &      +(cnpar*abs((vprv(l+1,k)-vprv(l,k))
     &		*(vpro(l+1,k)-vprv(l,k)))+
     &       (1.-cnpar)*abs((vpro(l+1,k)-vpro(l,k))
     &		*(vprv(l+1,k)-vpro(l,k))))
     &       /dh/h(l+1)
     &                )

            !m2 = du**2 + dv**2
            m2 = du + dv
            shearf2(l,k) = m2
          end do
        end do

	n2max = shympi_max(n2max)

	nfreq = sqrt(n2max)
	nperiod = 0.
	if( nfreq .gt. 0. ) nperiod = 1. / nfreq
	if( iuinfo .le. 0 ) call getinfo(iuinfo)
	if(shympi_is_master()) then
	  write(iuinfo,*) 'n2max: ',n2max,nfreq,nperiod
	end if

	end

c**************************************************************

	subroutine bnstress(czdef,taub)

c computes bottom stress at nodes
c
c this is evaluated for every element and then averaged for each node
c taub (stress at bottom) is accumulated and weighted by area
 
	use mod_hydro_vel
	use evgeom
	use levels
	use basin
	use shympi
	use mod_area

	implicit none

	real czdef
	real taub(nkn)

	integer k,ie,ii,n,nlev,lmin,ie_mpi
	real aj,taubot

c	---------------------------------------------------
c	initialize arrays
c	---------------------------------------------------

        taub = 0.
 
c	---------------------------------------------------
c	accumulate
c	---------------------------------------------------

        do ie_mpi=1,nel

	  ie = ip_sort_elem(ie_mpi)
 
          !call elebase(ie,n,ibase)
	  n = 3
          aj = 4. * ev(10,ie)
	  nlev = ilhv(ie)

          taubot = czdef * ( ulnv(nlev,ie)**2 + vlnv(nlev,ie)**2 )
          do ii=1,n
            k = nen3v(ii,ie)
            taub(k) = taub(k) + taubot * aj
          end do

	end do

!       shympi_elem: exchange taub
!       for area as weight use surface value areakv(lmin,k)
        !call shympi_comment('shympi_elem: exchange taub')
        call shympi_exchange_and_sum_2d_nodes(taub)

c	---------------------------------------------------
c	compute bottom stress
c	---------------------------------------------------

	do k=1,nkn
	  lmin = jlhkv(k)
	  if( areakv(lmin,k) <= 0. ) stop 'error stop bnstress: (2)'
	  taub(k) = taub(k) / areakv(lmin,k)
	end do

	end

c**************************************************************

	subroutine checka(nldim,shearf2,buoyf2,taub)

c checks arrays for nan or other strange values

	use mod_meteo
	use mod_diff_visc_fric
	use mod_hydro_print
	use mod_hydro_vel
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nldim
	real buoyf2(nldim,nkn)
	real shearf2(nldim,nkn)
	real taub(nkn)

	call nantest(nkn*nldim,shearf2,'shearf2')
	call nantest(nkn*nldim,buoyf2,'buoyf2')
	call nantest(nkn,taub,'taub')

	call nantest(nkn*(nldim+1),visv,'visv')
	call nantest(nkn*(nldim+1),difv,'difv')

	call nantest(nkn*nldim,uprv,'uprv')
	call nantest(nkn*nldim,vprv,'vprv')
	call nantest(nel*nldim,ulnv,'ulnv')
	call nantest(nel*nldim,vlnv,'vlnv')

	call nantest(nkn,tauxnv,'tauxnv')
	call nantest(nkn,tauynv,'tauynv')
 
	end

c**************************************************************

	subroutine checkb(text)

c checks arrays for strange values

	use mod_meteo
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) text

	integer one,three
	real zero,valmax
	real amin,amax

	one = 1
	three = 3
	zero = 0.

	valmax = 5.
	call valtest(one,nkn,-valmax,valmax,znv,text,'znv')
	call valtest(three,nel,-valmax,valmax,zenv,text,'zenv')

	valmax = 10000.
	call valtest(nlvdi+1,nkn,zero,valmax,visv,text,'visv')
	call valtest(nlvdi+1,nkn,zero,valmax,difv,text,'difv')

	valmax = 100.
	call valtest(nlvdi,nkn,-valmax,valmax,rhov,text,'rhov')

	valmax = 3.
	call valtest(nlvdi,nel,-valmax,valmax,ulnv,text,'ulnv')
	call valtest(nlvdi,nel,-valmax,valmax,vlnv,text,'vlnv')

	valmax = 100.
	call valtest(nlvdi,nel,-valmax,valmax,utlnv,text,'utlnv')
	call valtest(nlvdi,nel,-valmax,valmax,vtlnv,text,'vtlnv')

	end

c**************************************************************

	subroutine keps_shell

	use mod_keps
	use mod_layer_thickness
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,lmax,lmin,numOfLev,l
	real rho0,rhoair
	real cd,cb
	real wx,wy
	real ubot,vbot
	real dtreal,dt
	real taus,taub

	integer icall
	save icall
	data icall /0/

	if( icall .eq. 0 ) then
	  call keps_init
	  icall = 1
	end if

	call get_timestep(dtreal)
	dt = dtreal
	rho0 = 1024.
	rhoair = 1.225
	cd = 2.5e-3
	cb = 2.5e-3

	do k=1,nkn

	  lmax = ilhkv(k)
	  lmin = jlhkv(k)
	  numOfLev = lmax-lmin+1
	  ubot = uprv(lmax,k)
	  vbot = vprv(lmax,k)
	  call get_wind(k,wx,wy)
	  taus = rhoair * cd * (wx*wx+wy*wy)
	  taub = rho0 * cb * (ubot*ubot+vbot*vbot)

          call keps(numOfLev,dt,rho0
     +		,taus,taub
     +		,hdknv(lmin,k),uprv(lmin,k),vprv(lmin,k)
     +		,rhov(1,k),visv(0,k),difv(0,k),tken(0,k),eps(0,k))

	end do

	end

c**************************************************************

        subroutine keps_init

c initializes arrays for keps routine

	use mod_keps
	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

        implicit none

	real kmin,epsmin,lenmin,avumol,avtmol,avsmol
c       parameter(kmin=1.e-10,epsmin=1.e-12,lenmin=0.01)
        parameter(kmin=3.e-6,epsmin=5.e-10,lenmin=0.01)
        parameter(avumol=1.3e-6,avtmol=1.4e-7,avsmol=1.1e-9)

        integer k,l

	do k=1,nkn
          do l=0,nlvdi
            tken(l,k) = kmin
            eps(l,k) = epsmin
            rls(l,k) = lenmin
            visv(l,k) = avumol
            difv(l,k) = avsmol
	  end do
        end do

        end

c**************************************************************
c**************************************************************
c debug routines
c**************************************************************
c**************************************************************

	subroutine write_node_vel_info(iunit,it,k,ndim,depth
     +			,u_taus,u_taub,z0s,z0b
     +			,num,nuh,ken,dis,len,nn,ss,hh)

	implicit none

	integer iunit,it,k,ndim
	double precision depth,u_taus,u_taub,z0s,z0b

	double precision hh(0:ndim)
	double precision nn(0:ndim), ss(0:ndim)

	double precision num(0:ndim), nuh(0:ndim)
	double precision ken(0:ndim), dis(0:ndim)
	double precision len(0:ndim)

	integer l

	write(iunit,*) '----------------------------------'
	write(iunit,*) 'time,node,nlev: ',it,k,ndim
	write(iunit,*) 'depth: ',depth
	write(iunit,*) 'u_taus,u_taub: ',u_taus,u_taub
	write(iunit,*) 'z0s,z0b: ',z0s,z0b
	write(iunit,*) 'num:',(num(l),l=0,ndim)
	write(iunit,*) 'nuh:',(nuh(l),l=0,ndim)
	write(iunit,*) 'ken:',(ken(l),l=0,ndim)
	write(iunit,*) 'dis:',(dis(l),l=0,ndim)
	write(iunit,*) 'len:',(len(l),l=0,ndim)
	write(iunit,*) 'nn :',(nn(l),l=0,ndim)
	write(iunit,*) 'ss :',(ss(l),l=0,ndim)
	write(iunit,*) 'hh :',(hh(l),l=0,ndim)

	end

c**************************************************************

	subroutine write_node_bin_info(iunit,it,k,ndim,depth
     +			,u_taus,u_taub,z0s,z0b
     +			,num,nuh,ken,dis,len,nn,ss,hh)

	implicit none

	integer iunit,it,k,ndim
	double precision depth,u_taus,u_taub,z0s,z0b

	double precision hh(0:ndim)
	double precision nn(0:ndim), ss(0:ndim)

	double precision num(0:ndim), nuh(0:ndim)
	double precision ken(0:ndim), dis(0:ndim)
	double precision len(0:ndim)

	integer l

	write(iunit) it,k,ndim
	write(iunit) depth
	write(iunit) u_taus,u_taub
	write(iunit) z0s,z0b
	write(iunit) (num(l),l=0,ndim)
	write(iunit) (nuh(l),l=0,ndim)
	write(iunit) (ken(l),l=0,ndim)
	write(iunit) (dis(l),l=0,ndim)
	write(iunit) (len(l),l=0,ndim)
	write(iunit) (nn(l),l=0,ndim)
	write(iunit) (ss(l),l=0,ndim)
	write(iunit) (hh(l),l=0,ndim)

	end

c**************************************************************

        subroutine gotm_debug_write_f(iu,iwhat,k,nlev
     +			,dt,depth
     +                  ,u_taus,u_taub,z0s,z0b
     +                  ,hh,nn,ss,num,nuh,ken,dis,len)

        implicit none

        integer iu,iwhat,k,nlev
        double precision dt,depth
        double precision u_taus,u_taub
        double precision z0s,z0b
        double precision :: hh(0:nlev)
        double precision :: nn(0:nlev), ss(0:nlev)

        double precision :: num(0:nlev), nuh(0:nlev)
        double precision :: ken(0:nlev), dis(0:nlev)
        double precision :: len(0:nlev)

	integer ndim,l

	ndim = nlev

	write(iu,*) iwhat,k,nlev
	write(iu,*) dt,depth
	write(iu,*) u_taus,u_taub
	write(iu,*) z0s,z0b
	write(iu,*) (hh(l),l=0,ndim)
	write(iu,*) (nn(l),l=0,ndim)
	write(iu,*) (ss(l),l=0,ndim)
	write(iu,*) (num(l),l=0,ndim)
	write(iu,*) (nuh(l),l=0,ndim)
	write(iu,*) (ken(l),l=0,ndim)
	write(iu,*) (dis(l),l=0,ndim)
	write(iu,*) (len(l),l=0,ndim)
	flush(iu)

	end

c**************************************************************

        subroutine gotm_debug_write_u(iu,iwhat,k,nlev
     +			,dt,depth
     +                  ,u_taus,u_taub,z0s,z0b
     +                  ,hh,nn,ss,num,nuh,ken,dis,len)

        implicit none

        integer iu,iwhat,k,nlev
        double precision dt,depth
        double precision u_taus,u_taub
        double precision z0s,z0b
        double precision :: hh(0:nlev)
        double precision :: nn(0:nlev), ss(0:nlev)

        double precision :: num(0:nlev), nuh(0:nlev)
        double precision :: ken(0:nlev), dis(0:nlev)
        double precision :: len(0:nlev)

	integer ndim,l

	ndim = nlev

	write(iu) iwhat,k,nlev
	write(iu) dt,depth
	write(iu) u_taus,u_taub
	write(iu) z0s,z0b
	write(iu) (hh(l),l=0,ndim)
	write(iu) (nn(l),l=0,ndim)
	write(iu) (ss(l),l=0,ndim)
	write(iu) (num(l),l=0,ndim)
	write(iu) (nuh(l),l=0,ndim)
	write(iu) (ken(l),l=0,ndim)
	write(iu) (dis(l),l=0,ndim)
	write(iu) (len(l),l=0,ndim)
	flush(iu)

	end

c**************************************************************

        subroutine gotm_debug_read_u(iu,iwhat,k,nlev
     +			,dt,depth
     +                  ,u_taus,u_taub,z0s,z0b
     +                  ,hh,nn,ss,num,nuh,ken,dis,len,ierr)

        implicit none

        integer iu,iwhat,k,nlev,ierr
        double precision dt,depth
        double precision u_taus,u_taub
        double precision z0s,z0b
        double precision :: hh(0:nlev)
        double precision :: nn(0:nlev), ss(0:nlev)

        double precision :: num(0:nlev), nuh(0:nlev)
        double precision :: ken(0:nlev), dis(0:nlev)
        double precision :: len(0:nlev)

	integer ndim

	ndim = nlev
	ierr = 0

	read(iu,end=1) iwhat,k,nlev
	if( nlev > ndim ) goto 99
	read(iu) dt,depth
	read(iu) u_taus,u_taub
	read(iu) z0s,z0b
	read(iu) hh(0:nlev)
	read(iu) nn(0:nlev)
	read(iu) ss(0:nlev)
	read(iu) num(0:nlev)
	read(iu) nuh(0:nlev)
	read(iu) ken(0:nlev)
	read(iu) dis(0:nlev)
	read(iu) len(0:nlev)

    1	continue
	ierr = -1
	return
   99	continue
	write(6,*) ndim,nlev
	stop 'error stop gotm_run_debug: dimensions'
	end

c**************************************************************

