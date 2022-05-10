
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

! routines for coupling WW3 model
!
! revision log :
!
! 04.07.2019    ggu     written from scratch
! 24.03.2022    ggu     newly started with Aron
! 06.05.2022    ggu     module renamed to subww3

!===========================================================
        module subww3
!===========================================================

	implicit none

	logical, save :: bww3 = .false.

!===========================================================
        end module subww3
!===========================================================

!***********************************************************

	subroutine ww3_init
	end

!***********************************************************

	subroutine ww3_loop
	end

!***********************************************************

	subroutine ww3_finalize
	end

!***********************************************************

