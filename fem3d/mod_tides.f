
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

! revision log :
!
! 23.01.2015	ggu	changed VERS_7_1_4
! 10.07.2015	ggu	changed VERS_7_1_50
! 28.04.2016	ggu	changed VERS_7_5_9
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

	module mod_tides

	implicit none

	integer, private, save :: nkn_tides = 0

	real, allocatable, save :: xgeov(:)
	real, allocatable, save :: ygeov(:)
	real, allocatable, save :: xcartv(:)
	real, allocatable, save :: ycartv(:)
	real, allocatable, save :: zeqv(:)

	contains

!************************************************************

	subroutine mod_tides_init(nkn)

	integer nkn

	if( nkn == nkn_tides ) return

	if( nkn_tides > 0 ) then
	  deallocate(xgeov)
	  deallocate(ygeov)
	  deallocate(xcartv)
	  deallocate(ycartv)
	  deallocate(zeqv)
	end if

	nkn_tides = nkn

	if( nkn == 0 ) return

	allocate(xgeov(nkn))
	allocate(ygeov(nkn))
	allocate(xcartv(nkn))
	allocate(ycartv(nkn))
	allocate(zeqv(nkn))

	end subroutine mod_tides_init

!************************************************************

	end module mod_tides

