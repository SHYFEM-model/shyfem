
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

c*******************************************************/

	subroutine qfont00( font )

	character*(*) font
	character*257 nline
	common /nline/ nline

	call cline( font )
	!call qfont0( nline )

	end

c*******************************************************/

	subroutine qtext00( x , y , s )

	real x,y
	character*(*) s
	character*257 nline
	common /nline/ nline

	call cline( s )
	!call qtext0( x , y , nline )

	end

c*******************************************************/

	subroutine qtsize00( s , w , h )   ! NEW 

	real w,h
	character*(*) s
	character*257 nline
	common /nline/ nline

	call cline( s )
	!call qtsiz0( nline , w , h )   ! NEW

	end

c*******************************************************/

	subroutine qcomm( s )

	character*(*) s
	character*257 nline
	common /nline/ nline

	call cline( s )
	call qcomm0( nline )

	end

c*******************************************************/

	subroutine cline( line )

c converts line to c style (string zero terminated)

	implicit none

	character*(*) line
	character*1 c
	integer n,i
	character*257 nline
	common /nline/ nline

	n = len(line)
	do i=n,1,-1
	  c = line(i:i)
	  if( c .ne. '	' .and. c .ne. ' ' ) goto 1
	end do
    1	continue

	if( i .gt. 0 ) nline = line(1:i)
	if( i .lt. 257 ) i = i + 1
	nline(i:i) = char(0)

	end

c*******************************************************/

	subroutine string_length( line , length )

c determines length of fortran string

	implicit none

	character*(*) line
	integer length

	integer i,n
	character*1 c

	n = len(line)

	do i=n,1,-1
	  c = line(i:i)
	  if( c .ne. '	' .and. c .ne. ' ' ) goto 1
	end do
    1	continue

	length = i

	write(6,*) 'string_length: ',n,length,i

	end

c*******************************************************/

