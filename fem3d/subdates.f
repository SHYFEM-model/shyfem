
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016,2019  Georg Umgiesser
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

! custom dates for reset or other usage
!
! revision log :
!
! 09.09.2016	ggu	changed VERS_7_5_17
! 13.02.2019	ggu	isolated from more subroutines
! 16.02.2019	ggu	changed VERS_7_5_60
!
!**************************************************************

!==============================================================
	module custom_dates
!==============================================================

	implicit none

	logical, save, private :: bdebug = .true.
	integer, save, private :: idate = 0
	integer, save, private :: ndate = 0
	double precision, save, private, allocatable :: restime(:)

!==============================================================
	contains
!==============================================================

	subroutine custom_dates_init(atime,file)

	double precision atime
	character*(*) file

	integer ndim

	!write(6,*) 'custom: ',it,file

	if( idate /= 0 ) return		!called more than once

	idate = -1
	if( file == ' ' ) return			!no file given

	call get_custom_dates(file,0,ndim,restime)
	allocate(restime(ndim))
	call get_custom_dates(file,ndim,ndate,restime)

	write(6,*) 'custom dates used: ',ndate,'  ',trim(file)

	idate = 0

	do
	  idate = idate + 1
	  if( idate > ndate ) exit
	  if( atime < restime(idate) ) exit
	end do

	end subroutine custom_dates_init

!**************************************************************

	subroutine custom_dates_over(atime,bover)

	double precision atime
	logical bover

	character*80 file

!---------------------------------------------------------------
! initialize - convert date to relative time
!---------------------------------------------------------------

	bover = .false.

	if( idate == -1 ) return
	if( idate > ndate ) return

!---------------------------------------------------------------
! see if we have to reset
!---------------------------------------------------------------

	if( atime < restime(idate) ) return

!---------------------------------------------------------------
! ok, reset needed - advance to next reset time
!---------------------------------------------------------------

	do
	  idate = idate + 1
	  if( idate > ndate ) exit
	  if( atime < restime(idate) ) exit
	end do

	bover = .true.

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	end subroutine custom_dates_over

!**************************************************************

	subroutine get_custom_dates(file,ndim,n,atimes)

! gets dates from file and converts them to absolute time

	use iso8601

	character*(*) file
	integer ndim		!ndim==0 => check how many dates are given
	integer n		!on return total number of dates given
	double precision atimes(ndim)	!on return absolute dates given

	integer ianz,ios,nline,i,ierr
	integer date,time
	double precision d(2)
	double precision atime,atime_old
	character*80 line
	logical bcheck

	integer iscand

	n = 0
	nline = 0
	bcheck = ( ndim == 0 )		!only check, no dates returned

	open(1,file=file,status='old',form='formatted',iostat=ios)

	if( ios /= 0 ) then
	  write(6,*) 'cannot open custom reset file: ',trim(file)
	  stop 'error stop get_custom_dates: opening file'
	else if( .not. bcheck ) then
	  write(6,*) 'reading custom reset file: ',trim(file)
	end if

	do
	  read(1,'(a)',iostat=ios) line
	  nline = nline + 1
	  if( ios /= 0 ) exit
	  ianz = iscand(line,d,2)
	  if( ianz == 0 ) then
	    cycle
	  else if( ianz == 1 ) then
	    date = nint(d(1))
	    time = 0
	  else if( ianz == 2 ) then
	    date = nint(d(1))
	    time = nint(d(2))
	  else
            call string2date(line,date,time,ierr)
            if( ierr /= 0 ) then
              write(6,*) 'parse error in date string:'
              write(6,*) 'line: ',trim(line)
              write(6,*) 'file: ',trim(file)
	      stop 'error stop get_custom_dates: parse error'
            end if
	  end if

	  n = n + 1
	  if( bcheck ) cycle
	  if( n > ndim ) then
	    write(6,*) 'n,ndim: ',n,ndim
	    stop 'error stop get_custom_dates: dimension error ndim'
	  end if

	  call dts_to_abs_time(date,time,atime)
	  atimes(n) = atime		!insert absolute time
	end do

	if( ios > 0 ) then
	  write(6,*) 'read error...'
	  write(6,*) 'file: ',trim(file)
	  write(6,*) 'line number: ',nline
	  stop 'error stop get_custom_dates: read error'
	end if

	close(1)

	if( bcheck ) return	!nothing in atimes to be debuged

	if( bdebug ) then
	  write(6,*) 'custom reset times: ',n
	  atime_old = atimes(1) - 1
	  do i=1,n
	    atime = atimes(i)
	    call dts_format_abs_time(atime,line)
	    write(6,*) trim(line)
	    if( atime <= atime_old ) then
	      write(6,*) 'times in custom reset must be ascending...'
	      stop 'error stop get_custom_dates: wrong order'
	    end if
	  end do
	end if

	end subroutine get_custom_dates

!==============================================================
	end module custom_dates
!==============================================================

