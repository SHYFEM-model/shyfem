
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

!===========================================================
        module ww3_shyfem
!===========================================================

	implicit none

	logical, save :: bww3 = .false.

!===========================================================
        end module ww3_shyfem
!===========================================================

!***********************************************************

	subroutine ww3_init

	use ww3_shyfem
	use basin

	implicit none

	integer ier
	integer iwave

	real getpar

	iwave = nint(getpar('iwave'))
	bww3 = ( iwave == 11 )
	bww3 = .true.

	if( .not. bww3 ) return

! here setup of WW3 model

	end

!***********************************************************

	subroutine ww3_loop

	use ww3_shyfem
	use basin
	use mod_hydro
	use mod_meteo

	implicit none

	integer n,nsave

	if( .not. bww3 ) return

! here in time loop - exchange arrays

	call ww3_exchange_wind(nkn,wxv,wyv)

	n = 1
	nsave = n
	call ww3_exchange_info(n)
	write(6,*) 'ww3 before and after n  = ',nsave,n

	end

!***********************************************************

	subroutine ww3_finalize

	use ww3_shyfem
	use basin

	implicit none

	if( .not. bww3 ) return

! here finalize ww3

	end

!***********************************************************

