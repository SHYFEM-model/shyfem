
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

c***********************************************************

	subroutine memnod(x,y)

c saves deleted node

	implicit none

	include 'param.h'

	real x,y

	integer next
	common /nkone/ next
	real xegv(nkndim), yegv(nkndim)
	common /xegv/xegv, /yegv/yegv

	integer icall	!FIXME
	save icall
	data icall / 0 /

	if( icall .eq. 0 ) then	!FIXME
	  icall = 1
	  next = 0
	end if

	next = next + 1
	if( next .gt. nkndim ) then
	  stop 'error stop memnod: dimension nkndim'
	end if

	xegv(next) = x
	yegv(next) = y

	end

c***********************************************************

	subroutine primem(nbase)

c prints deleted nodes

	implicit none

	include 'param.h'

	integer nbase

	integer next
	common /nkone/ next
	real xegv(nkndim), yegv(nkndim)
	common /xegv/xegv, /yegv/yegv

	integer k,n,it

	next = 0	!do not print nodes

	do k=1,next
	    n = nbase + k
	    it = 65
	    write(1,'(i1,2i8,2f14.4)') 1,n,it,xegv(k),yegv(k)
	end do

	end

c***********************************************************

c	block data
c		?? initialize next = 0		!FIXME
c	end

c***********************************************************

