
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
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

! revision log :
!
! 11.05.2018	ggu	changed VERS_7_5_47
! 14.02.2019	ggu	changed VERS_7_5_56

	program array_assign

! tests array assignments

	implicit none

	interface
	subroutine write_array(a1,a2)
        integer a1(:),a2(:)
	end
	end interface


	integer i
	integer, allocatable :: a(:)
	integer, allocatable :: b(:)
	integer, allocatable :: c(:)

	allocate(a(5),b(7),c(9))

	do i=1,7
	  b(i) = i
	end do

	c = 0
	call write_array(a,b)
	call write_array(c,b)

	end

	subroutine write_array(a1,a2)

	implicit none

	integer a1(:),a2(:)

	a1 = a2
	write(6,*) a1

	end

