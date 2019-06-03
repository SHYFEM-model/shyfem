
!--------------------------------------------------------------------------
!
!    Copyright (C) 2019  Georg Umgiesser
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

! shyelab_extract.f: utility for extracting records
!
! revision log :
!
! 24.04.2019	ggu	written from scratch
! 21.05.2019	ggu	changed VERS_7_5_62
!
!***************************************************************

!===============================================================
	module shy_extract
!===============================================================

! parses string and decides which records to extract
!
! string must be in one of two formats:
!
!	i1,i2,i3,etc.		explicit records to be extracted are given
!	is..if..ie		start, frequency and end of records are given
!				if ie is not given, extracts to end

	implicit none

	logical, save, private :: debug = .true.

	integer, save, private :: mode = 0
	integer, save, private :: nextract = 0
	integer, save, private :: is = 0
	integer, save, private :: if = 0
	integer, save, private :: ie = 0
	integer, save, allocatable, private :: irecs(:)

!===============================================================
	contains
!===============================================================

	subroutine initialize_extract(sextract)

	implicit none

	character*(*) sextract

	integer ianz
	double precision daux(1)
	double precision, allocatable :: dval(:)

	integer i,irec
	integer iscand

	if( sextract == ' ' ) return		! no extract information
	
	ianz = iscand(sextract,daux,0)		! count numbers in string

	if( ianz < 0 ) then			! range given
	  mode = 2
	  call parse_extract_range(sextract,is,if,ie)
	else if( ianz == 0 ) then
	  goto 99
	else
	  mode = 1
	  nextract = ianz
	  allocate(dval(ianz))
	  ianz = iscand(sextract,dval,ianz)	! get numbers
	  nextract = nint(maxval(dval))
	  if( nextract <= 0 ) goto 99
	  allocate(irecs(nextract))
	  irecs = 0
	  do i=1,ianz
	    irec = nint(dval(i))
	    if( irec > nextract ) goto 99
	    if( irec < 1 ) goto 99
	    irecs(irec) = 1
	  end do
	end if

	if( debug ) then
	  write(6,*) 'debug information for extract:'
	  write(6,*) 'sextract: ',trim(sextract)
	  write(6,*) mode,nextract,is,if,ie
	  if( allocated(irecs) ) write(6,*) irecs
	end if

	return
  99	continue
	write(6,*) trim(sextract)
	stop 'error stop initialize_extract: error in record list'
	end subroutine initialize_extract

!***************************************************************

	function must_extract(nread)

	implicit none

	logical must_extract
	integer nread

	must_extract = .false.

	if( mode == 0 ) then
	  must_extract = .true.
	else if( mode == 1 ) then
	  if( nread > nextract ) return
	  if( irecs(nread) == 0 ) return
	  must_extract = .true.
	else if( mode == 2 ) then
	  if( ie > 0 .and. nread > ie ) return
	  if( mod(nread-is,if) /= 0 ) return
	  must_extract = .true.
	else
	  write(6,*) 'mode = ',mode
	  stop 'error stop must_extract: no such mode'
	end if

	end function must_extract

!***************************************************************

	subroutine parse_extract_range(string,is,if,ie)

	implicit none

	character*(*) string
	integer is			!start record
	integer if			!frequency of record
	integer ie			!end record

	integer ianz,i1,i2
	double precision d(1)
	integer iscand

	i1 = index(string,'..')	
	if( i1 == 0 ) goto 99
	i2 = index(string,'..',.true.)	
	if( i1 == i2 ) i2 = 0

	ianz = iscand(string(1:i1-1),d,1)
	if( ianz <= 0 ) goto 99
	is = nint(d(1))
	if( i2 > 0 ) then
	  ianz = iscand(string(i2+2:),d,1)
	  if( ianz <= 0 ) goto 99
	  ie = nint(d(1))
	  ianz = iscand(string(i1+2:i2-1),d,1)
	  if( ianz <= 0 ) goto 99
	  if = nint(d(1))
	else
	  ianz = iscand(string(i1+2:),d,1)
	  if( ianz <= 0 ) goto 99
	  if = nint(d(1))
	  ie = 0
	end if

	if( is <= 0 ) goto 98
	if( if <= 0 ) goto 98
	if( ie /= 0 .and. ie < is ) goto 98

	return
  98	continue
	write(6,*) trim(string)
	write(6,*) is,if,ie
	stop 'error stop parse_extract_range: error in parsed parameters'
  99	continue
	write(6,*) trim(string)
	stop 'error stop parse_extract_range: error in parsing'
	end subroutine parse_extract_range

!===============================================================
	end module shy_extract
!===============================================================

