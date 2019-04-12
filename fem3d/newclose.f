
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

c closing routines
c
c contents :
c
c subroutine sp136(ic)		opens and closes sections & inlets
c
c subroutine inclos		initializes closing sections
c subroutine rdclos(isc)	reads closing sections
c subroutine ckclos		post-processes closing sections
c subroutine prclos		prints info on closing sections
c subroutine tsclos		tests closing sections
c
c revision log :
c
c 26.07.1988	ggu	(introduction of ic)
c 15.12.1988	ggu	(bfluss,bspec)
c 20.12.1988	ggu	(iop=6,7,8 using level out of lagoon)
c 14.03.1990	ggu	(completely restructured)
c 31.03.1990	ggu	(test : change chezy values ^^^^)
c 27.08.1992	ggu	$$0 - not used
c 27.08.1992	ggu	$$1 - for new algorithm (const. form func.)
c 31.08.1992	ggu	$$impli - implicit time step
c 24.09.1992	ggu	$$2 - special technital (fluxes) -> deleted
c 24.09.1992	ggu	$$3 - chezy closing
c 29.09.1992	ggu	$$4 - remove writing of vectors (NAN)
c 25.03.1998	ggu	integrated changes from technital version
c 27.03.1998	ggu	utility routines for reading etc...
c 27.03.1998	ggu	dead code deleted, xv(1) -> xv(3,1)
c 27.03.1998	ggu	/bnd/ substituted by utility routine
c 29.04.1998    ggu     uses module for semi-implicit time-step
c 22.10.1999    ggu     volag copied to this file (only used here)
c 05.12.2001    ggu     fixed compiler error with -Wall -pedantic
c 09.12.2003    ggu     fix for icl=10 (FIX)
c 10.03.2004    ggu     RQVDT - value in rqv is now discharge [m**3/s]
c 11.10.2008	ggu	bug in call to nrdnxt (real instead of double p.)
c 14.03.2019	ggu	re-written for new framework
c
c************************************************************************

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

          integer :: ioper
          integer :: iact
          integer :: imode
          real :: scal
          double precision :: dstart

          real, allocatable :: hdep(:)
          integer, allocatable :: ieboc(:)

        end type entry

	integer, save :: nclose = 0	!total number of closing sects
	integer, save :: nb13 = 0	!unit for output

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
        pentry(id)%dsoft = 0.
        pentry(id)%dwait = 0.
        pentry(id)%zdate = 0.
        pentry(id)%vdate = 0.
        pentry(id)%zdiff = 0.

        pentry(id)%ioper = 0
        pentry(id)%iact = 0
        pentry(id)%imode = 0
        pentry(id)%scal = 0.
        pentry(id)%dstart = 0.

        end subroutine close_init_id

!******************************************************************

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
	real hkv(nkn)		!aux array
	integer index(nkn)

	logical, save :: binit = .false.

	if( binit ) return		!already initialized
	binit = .true.

	nbc = nbnds()
	call makehkv_minmax(hkv,0)	!get average depths

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

	end do

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
	write(nb13,*)

	end

!******************************************************************

	subroutine sp136(ic)

c opens and closes sections & inlets
c
c iclose	1 : closing by diminishing depth
c		2 : closing by diminishing flux
c		3 : closing by diminishing flux that depends on
c			water volumn in basin
c		4 : partial closing by changing chezy
c		...2 and 3 may be used only with ibnd != 0
c
c kboc		vector containing node numbers that define
c		...opening/closing section
c ibnd		number of open boundary to use for opening/closing
c		...(kboc and ibnd are mutually esclusive)
c
c kref		reference node in section (used for mode)
c kin,kout	inner and outer node for section (used for mode)
c kdir		node defining direction for icl=2,3 : direction = kdir-kref
c		...default for nodes kref=kin=kout=kdir and
c		...kref is middle node in section
c
c ibndz		number of boundary that is used to establish the value
c		...of zout instead of taking it from node kout
c		...(when a section is closed at an open boundary
c		...the value of zout at node kout is similar to zin
c		...and not to the value of z outside of the section)
c
c itb		vector containing: it1,imode1,it2,imode2...
c		...where it. is the time from when imode. is valid
c imode		= 100 * icl + iop
c
c icl		0:no closing  1:forced closing
c		2:vref=0 and change to positive dir. (empty basin)
c		3:vref=0 and change to negative dir. (full basin)
c		4:vref>vdate  5:vref<vdate
c		6:zout>zdate  7:zout<zdate
c		8:zin>zdate  9:zin<zdate
c		icl>20:immediate (e.g. 25 <=> icl=5 + immediate)
c
c iop		0:no opening  1:forced opening
c		2:zin>zout  3:zin<zout	(same as 4,5 with zdiff=0.)
c		4:zin-zout>zdiff  5:zin-zout<zdiff
c		6:zout>zdate  7:zout<zdate
c		iop>20:immediate (e.g. 25 <=> icl=5 + immediate)
c
c isoft		number of timesteps for which opening/closing
c		...has to be distributed
c mnstp		number of timesteps for which no opening/closing
c		...can be performed after last action
c
c zdate,zdiff	water level variables used in mode
c vdate		velocity variable used in mode

	use close
	use shympi

	implicit none

	integer ic		!0 if no change in configuration   (out)

	include 'mkonst.h'

	logical bnewmode,bfirst
	logical bclos,bopen,bimm

	integer id,ioper
        integer nsc
        integer iact,imode

	integer iclose
	integer icltot,ioptot,icl,iop
	real geyer
	character*20 aline

	double precision dtime,dstart,dend,dsoft,dwait
	double precision, parameter :: dflag = -999.

	real getpar

	iclose=nint(getpar('iclose'))
	if(iclose.le.0) return		!no closing enabled

	if( shympi_is_parallel() ) then
	  stop 'error stop sp136: cannot run in mpi mode'
	end if

	call close_init			!internally called only once

	ic=0				!&&&&   do not compute new uv-matrix
	nsc = nclose

	call get_act_dtime(dtime)
	call get_act_timeline(aline)

	write(nb13,*) '********** ',aline,' **********'
	write(nb13,*) '***********************************'
	write(nb13,*)

	do id=1,nsc

c get parameters

        ioper  = pentry(id)%ioper
        iact  = pentry(id)%iact
        imode  = pentry(id)%imode
        dsoft  = pentry(id)%dsoft
        dwait  = pentry(id)%dwait

c	new open & close mode

	bnewmode=.false.
        call get_new_mode(id,dtime,iact,imode,bnewmode)

!	----------------------------------------------
	write(nb13,*)
	write(nb13,'(1x,a,4i5)') 'id,iact,imode,ioper :'
     +				,id,iact,imode,ioper
	write(nb13,*)
!	----------------------------------------------

c	decide closing & opening

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
	else if( ioper == is_open ) then	!inlet is open
		call closing(id,icl,bclos)
	else if( ioper == is_closed ) then	!inlet is closed
		call opening(id,iop,bopen)
	end if

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
	  write(6,*) 'ioper,geyer : ',ioper,geyer
	  write(nb13,*) 'ioper,geyer : ',ioper,geyer

	  if( bfirst ) then
       		call get_new_mode(id,dflag,iact,imode,bnewmode)
		write(6,*) 'inlet ',id,' closed at ',aline
		write(nb13 ,*) 'inlet ',id,' closed at ',aline
	  end if
	end if

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
	  write(6,*) 'ioper,geyer : ',ioper,geyer
	  write(nb13,*) 'ioper,geyer : ',ioper,geyer

	  if( bfirst ) then
       		call get_new_mode(id,dflag,iact,imode,bnewmode)
		write(6,*) 'inlet ',id,' opened at ',aline
		write(nb13 ,*) 'inlet ',id,' opened at ',aline
	  end if
	end if

        pentry(id)%ioper = ioper
        pentry(id)%iact = iact
        pentry(id)%imode = imode

	end do

	write(nb13,*) '***********************************'

	return
	end

c********************************************************************
c********************************************************************
c********************************************************************

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

c********************************************************************

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

c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************

	subroutine inclos

c initializes closing sections

	use close

	implicit none

	call close_init_alloc

	end

c********************************************************************

	subroutine rdclos(isc)

c reads closing sections

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

	    if( iweich == 0 ) exit
	    if( iweich > 2 ) goto 98
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

c********************************************************************

	subroutine ckclos

c post-processes closing sections

	use close

	implicit none

	include 'close.h'

	logical bstop
	integer nsc,ivdim,ivful
	integer nbc
	integer i,j,k,id
	integer nkboc,ibnd,nitb
	integer knode,kint
	integer jkboc
	integer ibndz

	integer ipint
c	integer iround
	integer nbnds,nkbnds

	bstop = .false.

	nsc = nclose
	nbc = nbnds()

        do j=1,nsc

	   id = j

           nkboc = size(pentry(id)%kboc)
	   ibnd = pentry(id)%ibnd
           nitb = size(pentry(id)%itb)

c	   nodes of boundary

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

c	   various other checks

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

c********************************************************************

	subroutine prclos

c prints info on closing sections

	use close

	implicit none

	integer j,i,id,nsc

	nsc = nclose

	if( nsc .le. 0 ) return

	write(6,*)
        write(6,*) 'closing sections :'

        write(6,*) '...header :' ,nsc

        do j=1,nsc

	id = j

        write(6,*) pentry(id)%isc
        write(6,*) pentry(id)%kref
        write(6,*) pentry(id)%kout
        write(6,*) pentry(id)%kin
        write(6,*) pentry(id)%ibndz
        write(6,*) pentry(id)%ibnd
        write(6,*) pentry(id)%dsoft
        write(6,*) pentry(id)%dwait
        write(6,*) pentry(id)%zdate
        write(6,*) pentry(id)%vdate
        write(6,*) pentry(id)%zdiff

        end do

	end

c********************************************************************

	subroutine tsclos

c tests closing sections

	implicit none

        write(6,*) '/close/'

	call prclos

	end

c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************

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

c*****************************************************************

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

c*****************************************************************

        subroutine get_new_mode(id,dtime,iact,imode,bchange)

	use close

        implicit none

        integer id,iact,imode
	double precision dtime
	logical bchange

        integer, parameter :: icycle = 4
        double precision, parameter :: dflag = -999.
	integer nitb
	double precision dnext

	if( bchange ) return		!already changed mode in this time step

	nitb = size(pentry(id)%itb)

	do
		if( iact <= 0 ) exit
		dnext = pentry(id)%itb(iact)
		if( dtime == dflag ) then
		  if( dnext /= dflag ) exit
		else
		  if( dtime < dnext ) exit
		end if

		bchange = .true.
		imode = nint(pentry(id)%itb(iact+1))
		iact=iact+2

		if( iact .gt. nitb ) then	!no more data
		  if( icycle .gt. 0 ) then	!use old data
		    iact = iact - icycle
		  else				!no more closings
		    iact = 0
		  end if
		end if
	end do

        end

c*****************************************************************

	subroutine set_geyer(id,geyer)

	use close
	use mod_internal

	implicit none

	integer id
	real geyer

	integer nieboc,i,ie
	real g

	g = geyer
	g = min(g,1.)
	g = max(g,0.)

	nieboc = size(pentry(id)%ieboc)

	do i=1,nieboc
	  ie = pentry(id)%ieboc(i)
	  rcomputev(ie) = g
	end do

	end

c*****************************************************************

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

c*****************************************************************

