
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2019  Georg Umgiesser
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

c*********************************************************************

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60

c*********************************************************************

	subroutine dlist_init(n,ix,bcirc)

	implicit none

	integer n
	integer ix(2,1)
	logical bcirc

	integer i,is,ie

	if( n .eq. 0 ) return

	do i=2,n-1
	  ix(1,i) = i-1
	  ix(2,i) = i+1
	end do

	if( bcirc ) then
	  is = n
	  ie = 1
	else
	  is = 0
	  ie = 0
	end if

	if( n .eq. 1 ) then
	  ix(1,1) = is
	  ix(2,1) = ie
	else
	  ix(1,1) = is
	  ix(2,1) = 2
	  ix(1,n) = n-1
	  ix(2,n) = ie
	end if

	end

c*********************************************************************

	subroutine dlist_delete(n,ix,i)

	implicit none

	integer n
	integer ix(2,1)
	integer i

	integer ia,ib

	if( n .eq. 0 ) return

	if( n .gt. 1 ) then
	  ib = ix(1,i)
	  ia = ix(2,i)
	  if( ib .gt. 0 ) ix(2,ib) = ia
	  if( ia .gt. 0 ) ix(1,ia) = ib
	end if

	ix(1,i) = 0
	ix(2,i) = 0

	end

c*********************************************************************

