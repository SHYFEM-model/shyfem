
!--------------------------------------------------------------------------
!
!    Copyright (C) 1988,1990,1992,1998-1999,2001,2003-2004  Georg Umgiesser
!    Copyright (C) 2008,2019-2020  Georg Umgiesser
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

! closing routines
!
! contents :
!
! subroutine sp136(ic)		opens and closes sections & inlets
!
! subroutine inclos		initializes closing sections
! subroutine rdclos(isc)	reads closing sections
! subroutine ckclos		post-processes closing sections
! subroutine prclos		prints info on closing sections
! subroutine tsclos		tests closing sections
!
! revision log :
!
! 26.07.1988	ggu	(introduction of ic)
! 15.12.1988	ggu	(bfluss,bspec)
! 20.12.1988	ggu	(iop=6,7,8 using level out of lagoon)
! 14.03.1990	ggu	(completely restructured)
! 31.03.1990	ggu	(test : change chezy values ^^^^)
! 27.08.1992	ggu	$$0 - not used
! 27.08.1992	ggu	$$1 - for new algorithm (const. form func.)
! 31.08.1992	ggu	$$impli - implicit time step
! 24.09.1992	ggu	$$2 - special technital (fluxes) -> deleted
! 24.09.1992	ggu	$$3 - chezy closing
! 29.09.1992	ggu	$$4 - remove writing of vectors (NAN)
! 25.03.1998	ggu	integrated changes from technital version
! 27.03.1998	ggu	utility routines for reading etc...
! 27.03.1998	ggu	dead code deleted, xv(1) -> xv(3,1)
! 27.03.1998	ggu	/bnd/ substituted by utility routine
! 29.04.1998	ggu	uses module for semi-implicit time-step
! 22.10.1999	ggu	volag copied to this file (only used here)
! 05.12.2001	ggu	fixed compiler error with -Wall -pedantic
! 09.12.2003	ggu	fix for icl=10 (FIX)
! 10.03.2004	ggu	RQVDT - value in rqv is now discharge [m**3/s]
! 11.10.2008	ggu	bug in call to nrdnxt (real instead of double p.)
! 14.03.2019	ggu	re-written for new framework
! 12.04.2019	ggu	first finished draft
! 17.04.2019	ggu	tested on old venlag62
! 21.05.2019	ggu	changed VERS_7_5_62
! 05.03.2020	ggu	output streamlined, set rfmax here
! 26.03.2020	ggu	set vmax here
! 09.04.2020	ggu	increase non computing elements with ndist
! 06.12.2020	ggu	deleted close.h
! 18.05.2023	ccf	include reading closure from file
!
!************************************************************************

!==================================================================
        module close
!==================================================================

        implicit none

        type, private :: entry

          integer :: isc
          integer, allocatable :: kboc(:)
          double precision, allocatable :: itb(:)
          integer :: kref
          integer :: kdir
          integer :: kout
          integer :: kin
          integer :: ibndz
          integer :: ibnd
          double precision :: dsoft
          double precision :: dwait
          real :: zdate
          real :: vdate
          real :: zdiff
          character*80 :: cfile         !closure file name
          integer :: idfile             !closure file id

          integer :: ioper		!operation mode (see below is_*)
          integer :: iact		!where are we in itb
          integer :: imode
          real :: scal
          real :: geyer
          double precision :: dstart

          real, allocatable :: hdep(:)
          integer, allocatable :: ieboc(:)

          real, allocatable :: distfact(:)

        end type entry

	integer, save :: iclose = 0	!closing enabled?

	integer, save :: nclose = 0	!total number of closing sects
	integer, save :: nb13 = 0	!unit for output

	integer, save :: ndist = 0	!distance for not computing terms

        integer, save, private :: ndim = 0
        type(entry), save, allocatable :: pentry(:)

	integer, parameter :: is_closing = -1
	integer, parameter :: is_closed_waiting = -2
	integer, parameter :: is_closed = -3
	integer, parameter :: is_opening = +1
	integer, parameter :: is_open_waiting = +2
	integer, parameter :: is_open = +3

	integer, parameter :: is_operating = 1
	integer, parameter :: is_waiting = 2
	integer, parameter :: is_ready = 3

!==================================================================
        contains
!==================================================================

        subroutine close_init_alloc

        type(entry), allocatable :: paux(:)

        if( ndim == 0 ) then
          ndim = 10
          allocate(pentry(ndim))
          return
        else
          ndim = ndim*2
          allocate(paux(ndim))
          paux(1:ndim/2) = pentry(1:ndim/2)
          call move_alloc(paux,pentry)
        end if

        end subroutine close_init_alloc

!******************************************************************

        subroutine close_init_new_id(id)

        integer id

        nclose = nclose + 1
        if( nclose > ndim ) then
          call close_init_alloc
        end if
        id = nclose

        call close_init_id(id)

        end subroutine close_init_new_id

!******************************************************************

        subroutine close_init_id(id)

        integer id

        if( id > ndim ) then
          stop 'error stop close_init_id: ndim'
        end if

        pentry(id)%isc = 0
        pentry(id)%kref = 0
        pentry(id)%kout = 0
        pentry(id)%kin = 0
        pentry(id)%kdir = 0
        pentry(id)%ibndz = 0
        pentry(id)%ibnd = 0
        pentry(id)%dsoft = -1.
        pentry(id)%dwait = -1.
        pentry(id)%zdate = 0.
        pentry(id)%vdate = 0.
        pentry(id)%zdiff = 0.
        pentry(id)%cfile = ''
        pentry(id)%idfile = 0

        pentry(id)%ioper = 0
        pentry(id)%iact = 0
        pentry(id)%imode = 0
        pentry(id)%scal = 0.
        pentry(id)%dstart = 0.

        end subroutine close_init_id

!==================================================================
        end module close
!==================================================================

	subroutine close_init

	use arrays
	use close
	use basin

	implicit none

	logical belem
	integer idim,i,k,ie,id
	integer nkboc,ibnd,nitb,ibndz,nieboc
	integer ibtyp
	integer nvert,nbc
	integer kvert(10)

	integer nkbnds,nbnds,itybnd
	real, save :: rfmax = 1./300.		!time scale 5 minutes
	real, save :: vmax = 10.		!viscosity
	real hkv(nkn)				!aux array
	integer index(nkn)
	integer color(nkn)
	real dfact(nkn)

        logical bcfile
        double precision dtime
        character*80 file       !closure file name
        integer idfile          !closure file unit

	real getpar

	logical, save :: binit = .false.

	if( binit ) return			!already initialized
	binit = .true.

	iclose = nint(getpar('iclose'))
        if(iclose.le.0) return          !no closing enabled

	nbc = nbnds()
	call makehkv_minmax(hkv,0)		!get average depths

	call set_fric_max(rfmax)  !set maximum friction (inverse tiemscale)
	call set_vis_max(vmax)    !set maximum viscosity

	do id=1,nclose

        nkboc = size(pentry(id)%kboc)
	ibnd = pentry(id)%ibnd
        nitb = size(pentry(id)%itb)

!	-----------------------------------------------
!	setup nodes in closing section
!	-----------------------------------------------

	if( ibnd > 0 ) then	!use nodes from boundary
	  nkboc = nkbnds(ibnd)
	  call init_array(pentry(id)%kboc,nkboc)
	  call irbnds(ibnd,nkboc,idim,pentry(id)%kboc)
	  if( idim /= nkboc ) then
	    stop 'error stop close_init: internal error (1)'
	  end if
	end if

!	-----------------------------------------------
!	setup other nodes
!	-----------------------------------------------

	if( pentry(id)%kref <= 0 ) pentry(id)%kref = pentry(id)%kboc(1)
	if( pentry(id)%kdir <= 0 ) pentry(id)%kdir = pentry(id)%kref
	if( pentry(id)%kout <= 0 ) pentry(id)%kout = pentry(id)%kref
	if( pentry(id)%kin  <= 0 ) pentry(id)%kin  = pentry(id)%kref

	ibtyp = 0
	ibndz = pentry(id)%ibndz
	if(ibndz.gt.0) then
	  ibtyp=itybnd(ibndz)
	  if(iabs(ibtyp).eq.1) pentry(id)%kout = -ibndz
	end if

!	-----------------------------------------------
!	find depth on nodes and flag nodes in closing section
!	-----------------------------------------------

	call init_array(pentry(id)%hdep,nkboc)
	index = 0

	do i=1,nkboc
	  k = pentry(id)%kboc(i)
	  pentry(id)%hdep(i) = hkv(k)
	  index(k) = 1
	end do

!	-----------------------------------------------
!	find elements containing closing section nodes
!	-----------------------------------------------

	call init_array(pentry(id)%ieboc)
	nieboc = 0

	do ie=1,nel
	  call nindex(ie,nvert,kvert)
	  if( nvert > 3 ) stop 'error stop close_init: internal error (2)'
	  belem = .false.
	  do i=1,nvert
	    k = kvert(i)
	    if( index(k) > 0 ) belem = .true.
	  end do
	  if( belem ) then
	    call append_to_array_i(pentry(id)%ieboc,nieboc,ie)
	  end if
	end do

	call trim_array(pentry(id)%ieboc,nieboc)

!	-----------------------------------------------
!	compute distance array around closing nodes
!	-----------------------------------------------

	call init_array(pentry(id)%distfact,nel)
	color = 0.

	do i=1,nkboc
	  k = pentry(id)%kboc(i)
	  color(k) = 1
	end do

	call flood_fill_progessive(color)
	dfact = color
	call adjust_distfact(ndist,dfact,pentry(id)%distfact)

!	-----------------------------------------------
!	set info if nodes are open boundary
!	-----------------------------------------------

	ibndz = pentry(id)%ibndz

	if(ibnd.gt.nbc) goto 82
	if(ibndz.gt.nbc) goto 82
	if(ibndz.le.0.and.ibnd.gt.0) pentry(id)%ibndz = ibnd

!	-----------------------------------------------
!	initialize changing parameters
!	-----------------------------------------------

	pentry(id)%iact = 1
	pentry(id)%ioper = is_open	!always start open
	pentry(id)%imode = 0		!nothing before first event

!	-----------------------------------------------
!	initialize other parameters
!	-----------------------------------------------

        if( pentry(id)%dsoft < 0 ) pentry(id)%dsoft = 0.
        if( pentry(id)%dwait < 0 ) pentry(id)%dwait = pentry(id)%dsoft

!	-----------------------------------------------
!	initialize file
!	-----------------------------------------------

        call get_act_dtime(dtime)

        file = pentry(id)%cfile
        bcfile = ( file /= ' ' )
        if (bcfile) then
           call iff_ts_init(dtime,file,2,1,idfile)
           pentry(id)%idfile = idfile
        end if

	call close_info(id)

	end do

	write(6,*) 'closing sections initialized: ',nclose

	call check_dist		!checks and plots

!	-----------------------------------------------
!	end of routine
!	-----------------------------------------------

	return
   82	continue
	write(6,*) 'ibnd,ibndz,nbc: ',ibnd,ibndz,nbc
	stop 'error stop close_init: parameters'
	end

!******************************************************************

	subroutine close_info(id)

	use close

	implicit none

	integer id

	integer ii
	integer kref,kdir,kout,kin,nkboc,nieboc
        character*80 file
	integer ipext,ieext,ideffi

        if( nb13 == 0 ) then        !open file
          nb13=ideffi('datdir','runnam','.cls','form','new')
          if( nb13 .le. 0 ) then
            stop 'error stop close_info : Cannot open CLS file'
          end if
        end if

	kref = pentry(id)%kref
	kdir = pentry(id)%kdir
	kout = pentry(id)%kout
	kin  = pentry(id)%kin
	nkboc = size(pentry(id)%kboc)
	nieboc = size(pentry(id)%ieboc)
	file = pentry(id)%cfile

	write(nb13,*)
	write(nb13,*) 'first call for closing section nr. : ',id
	write(nb13,*)
	write(nb13,*) 'kout,kin  : ',ipext(kout),ipext(kin)
	write(nb13,*) 'kref,kdir : ',ipext(kref),ipext(kdir)
	write(nb13,*) 'kboc : ',nkboc
	write(nb13,*) (ipext(pentry(id)%kboc(ii)),ii=1,nkboc)
	write(nb13,*) 'hboc : ',nkboc
	write(nb13,*) (pentry(id)%hdep(ii),ii=1,nkboc)
	write(nb13,*) 'ieboc : ',nieboc
	write(nb13,*) (ieext(pentry(id)%ieboc(ii)),ii=1,nieboc)
	write(nb13,*) 'cfile : ',trim(file)
	write(nb13,*)

	flush(nb13)

	end

!******************************************************************

	subroutine close_handle(ic)

!---------------------------------------------------------------------------
!
! opens and closes sections & inlets (section $close#)
!
! iclose	1 : closing by diminishing depth
!		2 : closing by diminishing flux
!		3 : closing by diminishing flux that depends on
!			water volumn in basin
!		4 : partial closing by changing chezy
!		...2 and 3 may be used only with ibnd != 0
!
! kboc		vector containing node numbers that define
!		...opening/closing section
! ibnd		number of open boundary to use for opening/closing
!		...(kboc and ibnd are mutually esclusive)
!
! kref		reference node for water level (used for mode)
! kin,kout	inner and outer node for section (used for mode)
! kdir		node defining direction for icl=2,3,10: direction = kdir-kin
!		...default for nodes kref=kin=kout=kdir and
!		...kref is middle node in section
!
! ibndz		number of boundary that is used to establish the value
!		...of zout instead of taking it from node kout
!		...(when a section is closed at an open boundary
!		...the value of zout at node kout is similar to zin
!		...and not to the value of z outside of the section)
!
! itb		vector containing: it1,imode1,it2,imode2...
!		...where it. is the time from when imode. is valid
! imode		= 100 * icl + iop
!
! icl		0:no closing  1:forced closing
!		2:vin=0 and change to positive dir. (empty basin)
!		3:vin=0 and change to negative dir. (full basin)
!		4:vin>vdate   5:vin<vdate
!		6:zout>zdate  7:zout<zdate
!		8:zin>zdate   9:zin<zdate
!		10:zref>zdate and positive velocity direction
!		icl>20:immediate (e.g. 25 <=> icl=5 + immediate)
!
! iop		0:no opening      1:forced opening
!		2:zin>zout        3:zin<zout	(same as 4,5 with zdiff=0)
!		4:zin-zout>zdiff  5:zin-zout<zdiff
!		6:zout>zdate      7:zout<zdate
!		iop>20:immediate (e.g. 25 <=> icl=5 + immediate)
!
! dsoft		time [s] over which opening/closing
!		...has to be distributed
! dwait		time [s] after end of opening/closing operations
!		...where no other operation can be performed
!
! zdate,zdiff	water level variables used in mode
! vdate		velocity variable used in mode
!
! cfile 	Name of the file for the closure. In file is provided, the
!		reference nodes and levels are ignored. Format is:
!  		   datetime flag (with flag=0 open; flag=1 close)
!
!---------------------------------------------------------------------------

	use close
	use shympi
	use mod_internal

	implicit none

	integer ic		!0 if no change in configuration   (out)

	include 'mkonst.h'

	logical bnewmode,bfirst
	logical bclos,bopen,bimm

	integer id,ioper
        integer nsc
        integer iact,imode

	integer icltot,ioptot,icl,iop
	real geyer
	character*20 aline

	double precision dtime,dstart,dend,dsoft,dwait
	double precision, parameter :: dflag = -999.

	integer idfile
        logical, save :: bcfile
        real fclose

	if(iclose.le.0) return		!no closing enabled

	if( shympi_is_parallel() ) then
	  stop 'error stop sp136: cannot run in mpi mode'
	end if

!	----------------------------------------------
!	initialize
!	----------------------------------------------

	ic=0				!&&&&   do not compute new uv-matrix
	nsc = nclose

	rcomputev = 1			!set all elements to compute

	call get_act_dtime(dtime)
	call get_act_timeline(aline)

	write(nb13,*) '********** ',aline,' **********'
	write(nb13,*) '***********************************'
	write(nb13,*)

!	----------------------------------------------
!	start loop on closing sections
!	----------------------------------------------

	do id=1,nsc

! get parameters

        ioper  = pentry(id)%ioper
        iact   = pentry(id)%iact
        imode  = pentry(id)%imode
        dsoft  = pentry(id)%dsoft
        dwait  = pentry(id)%dwait
        idfile = pentry(id)%idfile
        bcfile = ( idfile /= 0 )

!	----------------------------------------------
!	new open & close mode
!	----------------------------------------------

	call print_closing_info(id,'std')
	bnewmode=.false.
        call get_new_mode(id,dtime,iact,imode,bnewmode)
	!if( bnewmode ) call print_closing_info(id,'***')

        if ( bcfile ) call iff_ts_intp1(idfile,dtime,fclose)

	write(nb13,*)
	write(nb13,'(1x,a,4i5)') 'id,iact,imode,ioper :'
     +				,id,iact,imode,ioper
	write(nb13,*)

!	----------------------------------------------
!	decide about closing & opening
!	----------------------------------------------

	icltot=imode/100
	ioptot=imode-100*icltot
	icl=icltot
	if(icl.gt.20) icl=icl-20
	iop=ioptot
	if(iop.gt.20) iop=iop-20

	bopen=.false.
	bclos=.false.
	bimm=.false.
	if(icltot.gt.20) bimm=.true.
	if(ioptot.gt.20) bimm=.true.
	bimm = bimm .or. (dsoft.le.0.)

	bfirst = .false.

	ioper = pentry(id)%ioper

	if( abs(ioper) == is_operating ) then	!in operation phase
		if( ioper == is_closing ) bclos=.true.
		if( ioper == is_opening ) bopen=.true.
	else if( abs(ioper) == is_waiting ) then	!in waiting phase
		!nothing
		if( ioper == is_closed_waiting ) bclos=.true.
		if( ioper == is_open_waiting ) bopen=.true.
	else if( ioper == is_open ) then	!inlet is open
                if (bcfile) then
                   bclos = (fclose == 1)
                else
                   call closing(id,icl,bclos)
                end if
	else if( ioper == is_closed ) then	!inlet is closed
               if (bcfile) then
                  bopen = (fclose == 0)
                else
                  call opening(id,iop,bopen)
                end if
	end if

!	----------------------------------------------
!	closing inlets
!	----------------------------------------------

	if(bclos) then		!close inlet
	  if( ioper == is_open ) then		!first step of closure
	    bfirst = .true.
	    ioper = is_closing
	    pentry(id)%dstart = dtime
	  end if
	  dstart = pentry(id)%dstart
	  dend = dstart + pentry(id)%dsoft
	  if( bimm ) dend = dstart
	  dwait = dend + pentry(id)%dwait

	  if( dstart == dend ) then
	    geyer = 0.
	  else
	    geyer = (dend-dtime)/(dend-dstart)
	  end if

	  if( geyer <= 0. ) ioper = is_closed_waiting
	  if( dtime >= dwait ) ioper = is_closed
	  call set_geyer(id,geyer)
	  !write(6,*) 'ioper,geyer : ',ioper,geyer
	  !write(67,*) 'ioper,geyer : ',ioper,geyer
	  write(nb13,*) 'ioper,geyer : ',ioper,geyer

	  if( bfirst ) then
       		call get_new_mode(id,dflag,iact,imode,bnewmode)
		write(6,*) 'inlet ',id,' closed at ',aline
		write(nb13 ,*) 'inlet ',id,' closed at ',aline
	  end if
	end if

!	----------------------------------------------
!	opening inlets
!	----------------------------------------------

	if(bopen) then		!open inlet
	  if( ioper == is_closed ) then		!first step of opening
	    bfirst = .true.
	    ioper = is_opening
	    pentry(id)%dstart = dtime
	  end if
	  dstart = pentry(id)%dstart
	  dend = dstart + pentry(id)%dsoft
	  if( bimm ) dend = dstart
	  dwait = dend + pentry(id)%dwait

	  if( dstart == dend ) then
	    geyer = 0.
	  else
	    geyer = (dend-dtime)/(dend-dstart)
	  end if

	  if( geyer <= 0. ) ioper = is_open_waiting
	  if( dtime >= dwait ) ioper = is_open
	  geyer = 1. - geyer
	  call set_geyer(id,geyer)
	  !write(6,*) 'ioper,geyer : ',ioper,geyer
	  !write(67,*) 'ioper,geyer : ',ioper,geyer
	  write(nb13,*) 'ioper,geyer : ',ioper,geyer

	  if( bfirst ) then
       		call get_new_mode(id,dflag,iact,imode,bnewmode)
		write(6,*) 'inlet ',id,' opened at ',aline
		write(nb13 ,*) 'inlet ',id,' opened at ',aline
	  end if
	end if

!	----------------------------------------------
!	handle no operations
!	----------------------------------------------

	if( bclos .and. bopen ) then
	  write(6,*) 'both closing and opening: ',bclos,bopen
	  stop 'error stop sp136: internal error (1)'
	else if( .not. bclos .and. .not. bopen ) then	!out of operation
	  if( ioper == is_closed ) then			!inlet is closed
	    call set_geyer(id,0.)			!keep closed
	  end if
	end if

!	----------------------------------------------
!	remember variable parameters
!	----------------------------------------------

	call compute_aver_geyer(id,pentry(id)%geyer)

        pentry(id)%ioper = ioper
        pentry(id)%iact = iact
        pentry(id)%imode = imode

	end do

	write(nb13,*) '***********************************'

!	----------------------------------------------
!	end of routine
!	----------------------------------------------

	end

!********************************************************************
!********************************************************************
!********************************************************************

	subroutine closing(id,icl,bclose)

! decides if to close inlet

	use close
	use mod_hydro
	use mod_hydro_print
	use basin

	implicit none

	integer id
	integer icl
	logical bclose

	integer kref,kout,kin,kdir
	real zdiff,zdate
	real zin,zout,zref
	real u,v,uvref2,dx,dy,scal,scalo
	real vdate

	real zvbnds

        kref = pentry(id)%kref
        kdir = pentry(id)%kdir
        kout = pentry(id)%kout
        kin  = pentry(id)%kin

	zdiff = pentry(id)%zdiff
	zdate = pentry(id)%zdate
	vdate = pentry(id)%vdate
	scalo = pentry(id)%scal

	zin=znv(kin)
        zref=znv(kref)
	if( kout > 0 ) then
	  zout=znv(kout)
	else
	  zout=zvbnds(-kout)
	end if

	u=up0v(kin)
	v=vp0v(kin)
	uvref2=u*u+v*v
	dx=xgv(kdir)-xgv(kin)
	dy=ygv(kdir)-ygv(kin)
	scal=u*dx+v*dy

	if(icl.eq.0) then
		!nothing
	else if(icl.eq.1) then
		bclose=.true.
	else if(icl.eq.2) then
		if(scal*scalo.lt.0.and.scal.gt.0) bclose=.true.
	else if(icl.eq.3) then
		if(scal*scalo.lt.0.and.scal.lt.0) bclose=.true.
	else if(icl.eq.4) then
		if(uvref2.gt.vdate*vdate) bclose=.true.
	else if(icl.eq.5) then
		if(uvref2.lt.vdate*vdate) bclose=.true.
	else if(icl.eq.6) then
		if(zout.gt.zdate) bclose=.true.
	else if(icl.eq.7) then
		if(zout.lt.zdate) bclose=.true.
	else if(icl.eq.8) then
		if(zin.gt.zdate) bclose=.true.
	else if(icl.eq.9) then
		if(zin.lt.zdate) bclose=.true.
	else if(icl.eq.10) then         !FIX 9.12.2003
		if(zref.gt.zdate.and.scal.gt.0.) bclose=.true.
        else
		write(6,*) 'icl = ',icl
		stop 'error stop closing: no such code for icl'
        end if

	pentry(id)%scal = scal

	end

!********************************************************************

	subroutine opening(id,iop,bopen)

! decides if to open inlet

	use close
	use mod_hydro

	implicit none

	integer id
	integer iop
	logical bopen

	integer kout,kin
	real zdiff,zdate
	real zin,zout,zref

	real zvbnds

        kout = pentry(id)%kout
        kin  = pentry(id)%kin

	zdiff = pentry(id)%zdiff
	zdate = pentry(id)%zdate

	zin=znv(kin)
	if( kout > 0 ) then
	  zout=znv(kout)
	else
	  zout=zvbnds(-kout)
	end if

        if(iop.eq.0) then
                !nothing
        else if(iop.eq.1) then
                bopen=.true.
        else if(iop.eq.2) then
                if(zin.gt.zout) bopen=.true.
        else if(iop.eq.3) then
                if(zin.lt.zout) bopen=.true.
        else if(iop.eq.4) then
                if(zin-zout.gt.zdiff) bopen=.true.
        else if(iop.eq.5) then
                if(zin-zout.lt.zdiff) bopen=.true.
        else if(iop.eq.6) then
                if(zout.gt.zdate) bopen=.true.
        else if(iop.eq.7) then
                if(zout.lt.zdate) bopen=.true.
        else
		write(6,*) 'iop = ',iop
		stop 'error stop opening: no such code for iop'
        end if

	end

!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************

	subroutine inclos

! initializes closing sections

	use close

	implicit none

	call close_init_alloc

	end

!********************************************************************

	subroutine rdclos(isc)

! reads closing sections

	use close
	use arrays

	implicit none

	integer isc		!number of actual section to read (in)

	character*6 name
	character*80 text
	real value
	double precision dvalue
	integer nsc
	integer id
	integer iweich
	integer ival
	integer ikboc,iitb

	integer nrdnxt

        call close_init_new_id(id)
	pentry(id)%isc = isc

	ikboc = 0
	iitb = 0

	call init_array(pentry(id)%kboc)
	call init_array(pentry(id)%itb)

        do

            iweich = nrdnxt(name,dvalue,text)
	    value = dvalue
            call to_lower(name)

            if( iweich .eq. 2 ) then
              if( name .ne. 'kboc' .and. name .ne. 'itb' ) goto 93
	    end if
            if( iweich .eq. 3 ) then
              if( name .ne. 'cfile' ) goto 93
	    end if

	    if( iweich == 0 ) exit
	    if( iweich > 3 ) goto 98
	    if( iweich < 0 ) goto 98

              if( name .eq. 'kboc' ) then
	        call append_to_array(pentry(id)%kboc,ikboc,nint(value))
              else if( name .eq. 'itb' ) then
	        call append_to_array(pentry(id)%itb,iitb,dvalue)
	      else if(name.eq.'kref') then
	        pentry(id)%kref=nint(value)
	      else if(name.eq.'kdir') then
	        pentry(id)%kdir=nint(value)
	      else if(name.eq.'kout') then
	        pentry(id)%kout=nint(value)
	      else if(name.eq.'kin') then
	        pentry(id)%kin=nint(value)
	      else if(name.eq.'ibndz') then
	        pentry(id)%ibndz=nint(value)
	      else if(name.eq.'ibnd') then
	        pentry(id)%ibnd=nint(value)
	      else if(name.eq.'dsoft') then
	        pentry(id)%dsoft=dvalue
	      else if(name.eq.'dwait') then
	        pentry(id)%dwait=dvalue
	      else if(name.eq.'zdate') then
	        pentry(id)%zdate=dvalue
	      else if(name.eq.'vdate') then
	        pentry(id)%vdate=dvalue
	      else if(name.eq.'zdiff') then
	        pentry(id)%zdiff=dvalue
	      else if(name.eq.'cfile') then
	        pentry(id)%cfile=text
	      else
		goto 96
	      end if

        end do

	call trim_array(pentry(id)%kboc,ikboc)
	call trim_array(pentry(id)%itb,iitb)

	return
   93   continue
        write(6,*) 'Variable not allowed in this context : ',name
        stop 'error stop : rdclos'
   96   continue
        write(6,*) 'Not recognized variable name : ',name
        stop 'error stop : rdclos'
   98   continue
	write(6,*) 'iweich = ',iweich
        write(6,*) 'Read error in closing section : ',isc
        stop 'error stop : rdclos'
	end

!********************************************************************

	subroutine ckclos

! post-processes closing sections

	use close

	implicit none

	logical bstop
	integer nsc,ivdim,ivful
	integer nbc
	integer i,j,k,id
	integer nkboc,ibnd,nitb
	integer knode,kint
	integer jkboc
	integer ibndz

	integer ipint
	integer nbnds,nkbnds

	bstop = .false.

	nsc = nclose
	nbc = nbnds()

        do j=1,nsc

	   id = j

           nkboc = size(pentry(id)%kboc)
	   ibnd = pentry(id)%ibnd
           nitb = size(pentry(id)%itb)

!	   nodes of boundary

           if(ibnd.gt.0.and.nkboc.gt.0) then
             write(6,'(a,i2,a)') ' section CLOSE ',j,' :'
             write(6,*) '   ibnd and kboc are mutually exclusive'
             bstop=.true.
           else if(ibnd.gt.0) then
             if(ibnd.gt.nbc) then
                write(6,'(a,i2,a)') ' section CLOSE ',j,' :'
                write(6,*) '   ibnd must have value of open boundary'
                write(6,*) '   ibnd : ',ibnd
                bstop=.true.
             end if
           else if(nkboc.gt.0.and.ibnd.eq.0) then
             do i=1,nkboc
                knode=pentry(id)%kboc(i)
                kint=ipint(knode)
                if( knode <= 0 .or. kint <= 0 ) then
                   write(6,'(a,i2,a)') ' section CLOSE ',j,' :'
                   write(6,*) '   node not found ',knode
                   bstop=.true.
                end if
                pentry(id)%kboc(i) = kint
             end do
           else
             write(6,'(a,i2,a)') ' section CLOSE ',j,' :'
             write(6,*) '   ibnd = ',ibnd,'   nkboc = ',nkboc
             write(6,*) '   No data read for kboc'
             bstop=.true.
           end if

!	   various other checks

	   ibndz = pentry(id)%ibndz
           if(ibndz.gt.nbc) then
                write(6,'(a,i2,a)') ' section CLOSE ',j,' :'
                write(6,*) '   ibndz must have value of open boundary'
                write(6,*) '   ibndz : ',ibndz
                bstop=.true.
           end if

           if(nitb.le.0) then
             write(6,'(a,i2,a)') ' section CLOSE ',j,' :'
             write(6,*) '   No data read for itb'
             bstop=.true.
           end if

           if(nitb.ne.(nitb/2)*2) then
             write(6,'(a,i2,a)') ' section CLOSE ',j,' :'
             write(6,*) '   Odd number of data read for itb'
             write(6,*) '   nitb = ',nitb
             bstop=.true.
           end if

	   call convk(pentry(id)%kref,bstop)
	   call convk(pentry(id)%kdir,bstop)
	   call convk(pentry(id)%kout,bstop)
	   call convk(pentry(id)%kin,bstop)
        end do

        if(bstop) stop 'error stop : ckclos'

	end

!********************************************************************

	subroutine prclos

! prints info on closing sections

	use close

	implicit none

	integer j,i,id,nsc

	nsc = nclose

	if( nsc .le. 0 ) return

	write(6,*)
        write(6,*) '====== info on closing sections ========='

        write(6,*) ' number of close sections:' ,nsc

        do j=1,nsc

	id = j

        write(6,*) ' close section number: ', pentry(id)%isc
        write(6,*) '   kref:  ', pentry(id)%kref
        write(6,*) '   kout:  ', pentry(id)%kout
        write(6,*) '   kin:   ', pentry(id)%kin
        write(6,*) '   ibndz: ', pentry(id)%ibndz
        write(6,*) '   ibnd:  ', pentry(id)%ibnd
        write(6,*) '   dsoft: ', pentry(id)%dsoft
        write(6,*) '   dwait: ', pentry(id)%dwait
        write(6,*) '   zdate: ', pentry(id)%zdate
        write(6,*) '   vdate: ', pentry(id)%vdate
        write(6,*) '   zdiff: ', pentry(id)%zdiff
        write(6,*) '   cfile: ', trim(pentry(id)%cfile)

        end do

	end

!********************************************************************

	subroutine tsclos

! tests closing sections

	implicit none

        write(6,*) '/close/'

	call prclos

	end

!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************

	subroutine set_zdate(iduse,zdate)

	use close

	implicit none

	integer iduse
	real zdate

        integer id,idstart,idend

	if( iduse <= 0 ) then
	  idstart = 1
	  idend = nclose
	else
	  idstart = iduse
	  idend = iduse
	end if

	do id=idstart,idend
	  pentry(id)%zdate = zdate
	end do

	end

!*****************************************************************

	subroutine get_zdate(iduse,zdate)

	use close

	implicit none

        integer iduse
	real zdate

	integer id

	id = iduse
	if( id <= 0 ) id = 1

	zdate = pentry(id)%zdate

	end

!*****************************************************************

        subroutine get_new_mode(id,dtime,iact,imode,bchange)

	use close

        implicit none

        integer id,iact,imode
	double precision dtime
	logical bchange

        double precision, parameter :: dflag = -999.
	integer nitb
	double precision dnext

	if( bchange ) return		!already changed mode in this time step
	if( iact <= 0 ) return

	if( dtime == dflag ) then
	    dnext = pentry(id)%itb(iact)
	    if( dnext /= dflag ) return
	    imode = nint(pentry(id)%itb(iact+1))
	    call get_new_iact(id,iact)
	    bchange = .true.
	else
	  do
	    dnext = pentry(id)%itb(iact)
	    if( dnext == dflag ) exit
	    if( dtime < dnext ) exit
	    imode = nint(pentry(id)%itb(iact+1))
	    call get_new_iact(id,iact)
	    bchange = .true.
	    if( iact <= 0 ) exit
	  end do
	end if

        end

!*****************************************************************

	subroutine get_new_iact(id,iact)

	use close

	implicit none

	integer id
	integer iact

        integer, parameter :: icycle = 4
	integer nitb

	if( iact <= 0 ) return

	nitb = size(pentry(id)%itb)

	iact = iact + 2

	if( iact .gt. nitb ) then	!no more data
	  if( icycle .gt. 0 ) then	!use old data
	    iact = iact - icycle
	  else				!no more closings
	    iact = 0
	  end if
	end if

	end

!*****************************************************************

	subroutine compute_aver_geyer(id,geyer)

! geyer is 0 for closed, 1 for open

	use basin
	use close
	use mod_internal

	implicit none

	integer id
	real geyer

	integer ie,n
	real f,dist

	n = 0
	f = 0.

	do ie=1,nel
	  dist = pentry(id)%distfact(ie)
	  if( dist <= 0. ) cycle
	  f = f + rcomputev(ie)
	  n = n + 1
	end do

	geyer = f / n

	end

!*****************************************************************

	subroutine set_geyer(id,geyer)

! geyer is 0 for closed, 1 for open

	use basin
	use close
	use mod_internal

	implicit none

	integer id
	real geyer

	integer nieboc,i,ie
	real g,f,dist

	g = geyer
	g = min(g,1.)
	g = max(g,0.)

	do ie=1,nel
	  dist = pentry(id)%distfact(ie)
	  if( dist <= 0. ) cycle
	  f = 1. + dist*(g-1.)
	  rcomputev(ie) = f
	end do

	!nieboc = size(pentry(id)%ieboc)
	!pentry(id)%geyer = g

	!do i=1,nieboc
	!  ie = pentry(id)%ieboc(i)
	!  rcomputev(ie) = g
	!end do

	end

!*****************************************************************

	subroutine convk(k,bstop)

	implicit none

	integer k
	logical bstop

	integer kint
	integer ipint

	bstop = .false.
	if( k <= 0 ) return

	kint = ipint(k)
	if( kint <= 0 ) then
	  write(6,*) 'cannot find node: ',k
	  bstop = .true.
	end if
	k = kint

	end

!*****************************************************************

	subroutine print_closing_info(id,text)

	use close
	use mod_hydro

	integer id
	character*(*) text

	integer iact,imode,ioper
	integer kref,kin,kout
	integer izin,izout,izref,izdate,iflux
	real zin,zout,zref,zdate,geyer,flux
	double precision dtime
	character*80 string1,string2
	character*20 aline
	integer, save :: icall = 0

	call get_act_dtime(dtime)

	iact = pentry(id)%iact
	imode = pentry(id)%imode
	ioper = pentry(id)%ioper

	kref = pentry(id)%kref
	kin = pentry(id)%kin
	kout = pentry(id)%kout

	zdate = pentry(id)%zdate
	geyer = pentry(id)%geyer
	scal  = pentry(id)%scal
	if( scal > 0. ) scal = 1.
	if( scal < 0. ) scal = -1.

	call get_barotropic_flux(id,flux)
	iflux = nint(flux)

	zin=znv(kin)
        zref=znv(kref)
	if( kout > 0 ) then
	  zout=znv(kout)
	else
	  zout=zvbnds(-kout)
	end if

	izin = nint(100.*zin)
	izout = nint(100.*zout)
	izref = nint(100.*zref)
	izdate = nint(100.*zdate)

!	if( id == 1 ) icall = icall + 1

!	string1 = '     time    id  iact imode ioper'
!     +			//' zdate  zref   zin  zout geyer  scal'

!	string2 = '      date and time   id iact mode oper'
!     +			//' zctr zref  zin zout geyer  scal   flux'

!	it = nint(dtime)
!	call get_timeline(dtime,aline)
!	if( mod(icall-1,50) == 0 ) write(70+id,*) trim(string1)
!	write(70+id,1000) it,id,iact,imode,ioper
!     +			,izdate,izref,izin,izout,geyer,scal

!	if( id == 1 ) write(66,*) trim(string1)
!	write(66,1000) it,id,iact,imode,ioper
!     +			,izdate,izref,izin,izout,geyer,scal
! 1000	format(i10,8i6,2f6.2)

!	if( id == 1 ) write(68,*) trim(string2)
!	write(68,1100) aline,id,iact,imode,ioper
!     +			,izdate,izref,izin,izout,geyer,scal,iflux
! 1100	format(a20,8i5,2f6.2,i7)

	end

!*****************************************************************

	subroutine adjust_distfact(ndist,fact,efact)

! sets up array efact - will be between 1 (barriers) and 0 (far field)

	use basin

	implicit none

	integer ndist
	real fact(nkn)
	real efact(nel)

	integer k,ie,ii
	real fdist,fdist2,f,fmax,fdiv

	fdist = ndist + 1
	!fdist2 = fdist/2
	fdist2 = 0
	fdiv = fdist - fdist2

	do k=1,nkn
	  if( fact(k) == 1 ) then		!starting point - closing 
	    fact(k) = 1
	  else if( fact(k) < fdist2 ) then
	    fact(k) = 1
	  else if( fact(k) > fdist ) then
	    fact(k) = 0
	  else
	    f = fact(k) - fdist2
	    fact(k) = 1. - f/fdiv
	  end if
	end do

	where( fact > 1 ) fact = 1
	where( fact < 0 ) fact = 0

	do ie=1,nel
	  fmax = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    f = fact(k)
	    fmax = max(fmax,f)
	  end do
	  efact(ie) = fmax
	end do

	end

!*****************************************************************

	subroutine check_dist

	use basin
	use close

	implicit none

	integer id
	real dist(nel)
	integer check(nel)

	check = 0
	do id=1,nclose
	  dist = dist + pentry(id)%distfact
	  where( pentry(id)%distfact > 0 ) check = check + 1
	end do

	if( maxval(check) > 1 ) then
	  write(6,*) 'barriers are too close for chosen value of ndist'
	  stop 'error stop check_dist: ndist too big'
	end if

	call basin_to_grd
	call grd_flag_depth
	call grd_set_element_depth(dist)
	call grd_write('dist_close.grd')

	end

!*****************************************************************

