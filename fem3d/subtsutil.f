
!--------------------------------------------------------------------------
!
!    Copyright (C) 2014-2015,2017-2019  Georg Umgiesser
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

! utility routines for time series handling
!
! revision log :
!
! 08.12.2020	ggu	started from subtsfile.f
! 12.12.2020	ggu	allow setting absolute time
! 24.11.2021	ggu	documentation added
!
!*************************************************************

! notes :
!
! no change in time step is allowed in files
! atime0 is set by default to 0 (dtime is absolute time)
!
!        subroutine ts_util_set_zero_time(atime0)
!
!        subroutine ts_util_open_file(file,intpol,dtime0,id,ierr)
!        subroutine ts_util_read_records(id,dtime0,ierr)
!        subroutine ts_util_close_file(id)
!
!        subroutine ts_util_get_max(id,ivar,dtime0,period,rmax)
!        subroutine ts_util_get_val(id,ivar,nval,dtimes,vals)
!        subroutine ts_util_interpolate(id,ivar,dtime0,rintp)

! subiso8601.f subtsfile.f subscn.f subfil.f subdts.f

!=============================================================
	module ts_util
!=============================================================

        implicit none

        type, private :: entry

          integer :: iunit
          integer :: nvar
          integer :: ntime
          integer :: istart
          double precision :: ddt
          integer :: period
          double precision :: atime0
          double precision, allocatable :: dtime(:)
          real, allocatable :: data(:,:)

        end type entry

        integer, save, private :: nutil = 0
        integer, save, private :: ndim = 0
        double precision, save :: ts_atime0 = 0.
        type(entry), save, allocatable :: pentry(:)
 
	!logical, parameter :: btsdebug = .true.
	logical, parameter :: btsdebug = .false.
	logical, parameter :: btsdebmin = .true. .or. btsdebug
	!logical, parameter :: btsdebmin = .false. .or. btsdebug

!=============================================================
	contains
!=============================================================

        subroutine ts_util_init_alloc

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

        end subroutine ts_util_init_alloc

!*************************************************************

        subroutine ts_util_init_new_id(id)

        integer id

        nutil = nutil + 1
        if( nutil > ndim ) then
          call ts_util_init_alloc
        end if
        id = nutil

        call ts_util_init_id(id)

        end subroutine ts_util_init_new_id

!*************************************************************

        subroutine ts_util_init_id(id)

        integer id

        if( id > ndim ) then
          stop 'error stop ts_util_init_id: ndim'
        end if

        pentry(id)%iunit = 0
        pentry(id)%nvar = 0
        pentry(id)%period = 0
        pentry(id)%ddt = 0
        pentry(id)%ntime = 0
        pentry(id)%istart = 0
        pentry(id)%atime0 = 0

        end subroutine ts_util_init_id

!=============================================================
	end module ts_util
!=============================================================

!*************************************************************

        subroutine ts_util_set_zero_time(atime0)

! sets zero time to which dtime in time series is referred to

	use ts_util

	implicit none

	double precision atime0

	ts_atime0 = atime0

	end

!*************************************************************
!*************************************************************
!*************************************************************

        subroutine ts_util_open_file(file,intpol,dtime0,id,ierr)

! opens time series file for read
!
! intpol is either positive or negative
!	intpol > 0	only 2 or 4 for linear and cubic interpolation
!	intpol < 0	period=-intpol, in secs the period to be read in
! negative intpol is used to determin max in this period

	use ts_util

	implicit none

	character*(*) file		!file name
	integer intpol			!interpolation (if>0) or period (if<0)
	double precision dtime0		!time to start with
	integer id			!file id (return)
	integer ierr			!error code (return)

	integer iunit
	integer nvar
	integer ntime
	integer period
	integer istart,i
	integer datetime(2)
	double precision ddt,dtime,dtime1
	character*80 varline
	character*20 aline

	real, allocatable :: f(:),f1(:)

!----------------------------------------
! open file
!----------------------------------------

	id = 0
	ierr = 95
	call ts_open_file(file,nvar,datetime,varline,iunit)
	if( iunit <= 0 ) return
	if( nvar <= 0 ) return

!----------------------------------------
! find time step of file data and allocate arrays
!----------------------------------------

	allocate(f(nvar),f1(nvar))
	call ts_read_next_record(iunit,nvar,dtime,f,datetime,ierr)
	dtime = dtime - ts_atime0
	if( ierr /= 0 ) return
	call ts_peek_next_record(iunit,nvar,dtime1,f1,datetime,ierr)
	dtime1 = dtime1 - ts_atime0
	if( ierr /= 0 ) return

	if( btsdebmin ) then
	  write(6,*) 'open file:'
	  call dts_format_abs_time(dtime,aline)
	  write(6,*) dtime,aline
	  call dts_format_abs_time(dtime1,aline)
	  write(6,*) dtime1,aline
	end if


	if( dtime > dtime0 ) then
	  write(6,*) 'error: required time not in file'
	  call dts_format_abs_time(dtime,aline)
	  write(6,*) 'dtime;  ',dtime,aline
	  call dts_format_abs_time(dtime0,aline)
	  write(6,*) 'dtime0: ',dtime0,aline
	  ierr = 77
	  return
	end if

	istart = 1
	ddt = dtime1 - dtime
	if( intpol > 0 ) then			!number of points to keep
	  ntime = intpol
	  period = 0
	  if( intpol /= 2 .and. intpol /= 4 ) then
	    write(6,*) 'intpol = ',intpol
	    stop 'error stop ts_util_open_file: intpol must be 2 or 4'
	  end if
	else if( intpol < 0 ) then		!period of time given
	  period = -intpol
	  ntime = 2 + nint(period/ddt)
	else
	  write(6,*) 'error in intpol. Cannot be 0'
	  ierr = 79
	  return
	end if

        call ts_util_init_new_id(id)

        pentry(id)%iunit = iunit
        pentry(id)%nvar = nvar
        pentry(id)%period = period
        pentry(id)%ddt = ddt
        pentry(id)%ntime = ntime
        pentry(id)%istart = 1
	allocate(pentry(id)%dtime(ntime))
	allocate(pentry(id)%data(nvar,ntime))

	if( btsdebmin ) then
	  write(6,*) 'parameters:'
	  write(6,*) nvar,ntime
	  write(6,*) period,nint(ddt)
	end if

!----------------------------------------
! re-open file
!----------------------------------------

	close(iunit)
	call ts_open_file(file,nvar,datetime,varline,iunit)

!----------------------------------------
! fill data arrays
!----------------------------------------

	if( btsdebug ) then
	  write(6,*) 'filling arrays:'
	end if

	do i=1,ntime
	  call ts_read_next_record(iunit,nvar,dtime,f,datetime,ierr)
	  if( ierr /= 0 ) return
	  dtime = dtime - ts_atime0
	  pentry(id)%dtime(i) = dtime
	  pentry(id)%data(:,i) = f(1:nvar)
	  if( btsdebug ) then
	    call dts_format_abs_time(dtime,aline)
	    write(6,*) i,dtime,aline
	  end if
	end do

!----------------------------------------
! populate records
!----------------------------------------

	call ts_util_read_records(id,dtime0,ierr)

!----------------------------------------
! end of routine
!----------------------------------------

	end

!*************************************************************

	subroutine ts_util_read_records(id,dtime0,ierr)

! reads new values for time dtime0

	use ts_util

	implicit none

	integer id			!file id
	double precision dtime0		!read for this new time value
	integer ierr			!error code (return)

	logical bcubic
	integer istart,ilast,i,i1,i2
	integer iunit,nvar,ntime
	integer period
	integer datetime(2)
	real f(pentry(id)%nvar)
	double precision ddt,dtime,dtime1
	character*20 aline

!----------------------------------------
! set parameters
!----------------------------------------

	iunit = pentry(id)%iunit
	nvar = pentry(id)%nvar
	ntime = pentry(id)%ntime
	period = pentry(id)%period
	istart = pentry(id)%istart
	bcubic = ( period == 0 .and. ntime == 4 )
        ddt = pentry(id)%ddt

!----------------------------------------
! iterate to desired time
!----------------------------------------

	do
	  i = 2
	  if( bcubic ) i = 3
	  if( pentry(id)%dtime(i) >= dtime0 ) exit
	  call ts_read_next_record(iunit,nvar,dtime,f,datetime,ierr)
	  if( ierr /= 0 ) return
	  dtime = dtime - ts_atime0
	  do i=1,ntime-1
	    pentry(id)%dtime(i) = pentry(id)%dtime(i+1)
	    pentry(id)%data(:,i) = pentry(id)%data(:,i+1)
	  end do
	  pentry(id)%dtime(ntime) = dtime
	  pentry(id)%data(:,ntime) = f(1:nvar)
	  if( btsdebug ) then
	    call dts_format_abs_time(dtime,aline)
	    write(6,*) 'reading: ',dtime,aline
	  end if
	end do

!----------------------------------------
! sanity check
!----------------------------------------

	if( dtime0 < pentry(id)%dtime(1) ) goto 99
	if( dtime0 + period > pentry(id)%dtime(ntime) ) goto 99

	dtime = pentry(id)%dtime(1)
	do i=2,ntime
	  dtime1 = pentry(id)%dtime(i)
	  if( dtime1-dtime /= ddt ) then
	    write(6,*) 'change in time step not allowed'
	    write(6,*) ddt,dtime1-dtime
	    call dts_format_abs_time(dtime,aline)
	    write(6,*) 'at position: ',dtime,aline
	    goto 99
	  end if
	  dtime = dtime1
	end do

	if( bcubic ) then
	  if( dtime0 < pentry(id)%dtime(2) ) goto 98
	  if( dtime0 > pentry(id)%dtime(3) ) goto 98
	end if

	if (btsdebug ) then
	  write(6,*) 'sanity check:'
	  do i=1,ntime
	    dtime = pentry(id)%dtime(i)
	    call dts_format_abs_time(dtime,aline)
	    write(6,*) i,dtime,aline
	  end do
	end if

!----------------------------------------
! set final error code and handle errors
!----------------------------------------

	ierr = 0

	return
   99	continue
	dtime = dtime0
	call dts_format_abs_time(dtime,aline)
	write(6,*) 'looking for time: ',dtime,aline
	if( period /= 0 ) then
	  dtime = dtime0 + period
	  call dts_format_abs_time(dtime,aline)
	  write(6,*) 'maximum time needed: ',dtime,aline
	end if
	write(6,*) 'available times: '
	do i=1,ntime
	  dtime = pentry(id)%dtime(i)
	  call dts_format_abs_time(dtime,aline)
	  write(6,*) i,dtime,aline
	end do
	ierr = 55
	return
   98	continue
	write(6,*) 'cubic interpolation: out of time window'
	dtime = pentry(id)%dtime(2)
	call dts_format_abs_time(dtime,aline)
	write(6,*) 'minimum time: ',dtime,aline
	dtime = pentry(id)%dtime(3)
	call dts_format_abs_time(dtime,aline)
	write(6,*) 'maximum time: ',dtime,aline
	ierr = 57
	return

!----------------------------------------
! end of routine
!----------------------------------------

	end

!*************************************************************

	subroutine ts_util_close_file(id)

! closes time series file

	use ts_util

	implicit none

	integer id			!file id

	pentry(id)%iunit = 0

	end

!*************************************************************
!*************************************************************
!*************************************************************

	subroutine ts_util_get_max(id,ivar,dtime0,period,rmax)

! gets maximum of variable ivar

	use ts_util

	implicit none

	integer id			!file id
	integer ivar			!number of variable to look at
	double precision dtime0		!start of time to look for maximum
	integer period			!period to look for (0=>all data) [s]
	real rmax			!maximum found (return)

	integer i,ntime,nvar
	double precision dtime,dtime1
	real, parameter :: high = 1.e+30

	nvar = pentry(id)%nvar
	ntime = pentry(id)%ntime
	dtime1 = dtime0 + period
	if( period == 0 ) dtime1 = dtime0 + high
	rmax = -high
	!write(6,*) ntime,period,ivar

	if( ivar>nvar .or. ivar<1 ) then
	  write(6,*) 'ivar=',ivar,'  nvar=',nvar
	  stop 'error stop ts_util_get_max: ivar'
	else if( period < 0 ) then
	  write(6,*) 'period=',period
	  stop 'error stop ts_util_get_max: period'
	end if

	do i=1,ntime
	  dtime = pentry(id)%dtime(i)
	  if( dtime < dtime0 ) cycle
	  if( dtime > dtime1 ) exit
	  !write(6,*) i,dtime,dtime0,dtime1,rmax
	  rmax = max(rmax,pentry(id)%data(ivar,i))
	end do

	end

!*************************************************************

	subroutine ts_util_get_val(id,ivar,nval,dtimes,vals)

! returns time and values of variable ivar

	use ts_util

	implicit none

	integer id			!file id
	integer ivar			!number of variable to look at
	integer nval			!dimension on entry, filling on exit
	double precision dtimes(nval)	!time values (return)
	real vals(nval)			!values of variable ivar (return)

	integer i,ntime,nvar

	ntime = pentry(id)%ntime
	nvar = pentry(id)%nvar
	nval = min(nval,ntime)

	if( ivar>nvar .or. ivar<1 ) then
	  write(6,*) 'ivar=',ivar,'  nvar=',nvar
	  stop 'error stop ts_util_get_val: ivar'
	else if( nval < 1 ) then
	  write(6,*) 'nval=',nval
	  stop 'error stop ts_util_get_val: nval'
	end if

	do i=1,nval
	  dtimes(i) = pentry(id)%dtime(i)
	  vals(i) = pentry(id)%data(ivar,i)
	end do

	end

!*************************************************************

        subroutine ts_util_interpolate(id,ivar,dtime0,rintp)

	use ts_util

	implicit none

	integer id			!file id
	integer ivar			!number of variable to look at
	double precision dtime0		!time of which to interpolate
	real rintp			!interpolated value (return)

	integer nvar,nintp,period,i
	double precision dtimes(pentry(id)%ntime)
	double precision vals(pentry(id)%ntime)

	double precision neville_intern

	nvar = pentry(id)%nvar
	nintp = pentry(id)%ntime
	period = pentry(id)%period

	if( ivar>nvar .or. ivar<1 ) then
	  write(6,*) 'ivar=',ivar,'  nvar=',nvar
	  stop 'error stop ts_util_interpolate: ivar'
	else if( period /= 0 ) then
	  write(6,*) 'period=',period
	  write(6,*) 'file has been opened to look for maximum'
	  stop 'error stop ts_util_interpolate: period'
	end if

	do i=1,nintp
	  dtimes(i) = pentry(id)%dtime(i)
	  vals(i) = pentry(id)%data(ivar,i)
	end do

        rintp = neville_intern(nintp,dtimes,vals,dtime0)

	end

!*************************************************************

        function neville_intern(nintp,xa,ya,x)

! use Neville algorithm for Lagrangian interpolation (double precision)
!
! internal routine
!
! use nintp=4 for cubic interpolation
! use nintp=2 for linear interpolation

        implicit none

        double precision neville_intern         !interpolated value (return)
        integer nintp                           !grade of interpolation
        double precision xa(0:nintp-1)          !x points to use
        double precision ya(0:nintp-1)          !y points to use
        double precision x                      !x value where to interpolate

        integer i,k,n
        double precision xl,xh
        double precision p(0:nintp-1)

        n = nintp - 1
        p = ya

        do k=1,n
          do i=n,k,-1
            xl = xa(i-k)
            xh = xa(i)
            p(i) = ( (x-xl)*p(i) - (x-xh)*p(i-1) ) / (xh-xl)
          end do
        end do

        neville_intern = p(n)

        end

!*************************************************************
!*************************************************************
!*************************************************************

	subroutine ts_util_test

	use iso8601

	implicit none

	integer intpol,id,ierr,ih,imax,ivar,nval,i
	integer date,time
	character*80 file
	character*20 aline
	real rmax
	double precision dtime,dtime0,atime0

	integer, parameter :: ndim = 100
	double precision dtimes(ndim)
	real vals(ndim)

	real r
	integer icm
	icm(r) = nint(100.*r)

	call dtsyear(2020)

	file = 'Punta_Salute.txt'
	file = 'test.txt'
	intpol = 4			!cubic interpolation
	intpol = -4*3600		!4 hours
	intpol = -2*3600		!2 hours

	ivar = 1

	call dts_to_abs_time(20200101,0,atime0)	!zero time is 2020-01-01
	call dts_to_abs_time(20200701,0,dtime0)	!start time is 2020-07-01
	dtime0 = dtime0 - atime0

        call ts_util_set_zero_time(atime0)
        call ts_util_open_file(file,intpol,dtime0,id,ierr)
	if( ierr /= 0 ) stop 'error stop opening file'

	dtime = dtime0
	call dts_format_abs_time(dtime,aline)
	write(6,*) 'first time to read: ',dtime,aline

	imax = 10
	imax = 365*24

	do ih=1,imax
	  dtime = dtime0 + ih*3600
	  call dts_format_abs_time(dtime,aline)
	  !write(6,*) 'next time to read: ',dtime,aline
	  call ts_util_read_records(id,dtime,ierr)
	  if( ierr < 0 ) exit
	  if( ierr > 0 ) stop 'error stop reading file'
	  call ts_util_get_max(id,ivar,dtime,0,rmax)
	  nval = ndim
	  call ts_util_get_val(id,ivar,nval,dtimes,vals)
	  !write(6,*) aline,nval,rmax
	  !write(6,1000) icm(rmax),(icm(vals(i)),i=1,nval)
	  !write(6,2000) aline(1:14),icm(rmax),(icm(vals(i)),i=1,nval)
	end do

	call dts_format_abs_time(dtime,aline)
	write(6,*) 'last time read: ',dtime,aline

	call ts_util_close_file(id)

	stop
 1000	format(i4,2x,100i4)
 2000	format(a,i4,2x,100i4)
	end

!*************************************************************
!	program ts_util_main
!	call ts_util_test
!	end
!*************************************************************

