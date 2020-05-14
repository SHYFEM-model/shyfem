
!--------------------------------------------------------------------------
!
!    Copyright (C) 2019-2020  Georg Umgiesser
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

! routines for BFM module - dummy routine
!
! revision log :
!
! 02.10.2019	ggu	dummy routine introduced
! 22.04.2020	ggu	error if called

!******************************************************************

	subroutine bfm_init_internal

	write(6,*) 'bfm routines have not been linked into shyfem'
	write(6,*) 'please correct Rules.make file'
	stop 'error stop bfm_init_internal: no bfm routines'

	end

	subroutine bfm_reactor_internal
	end

!******************************************************************

