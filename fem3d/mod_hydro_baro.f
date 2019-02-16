
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

	module mod_hydro_baro

	implicit none

	integer, private, save :: nel_hydro_baro = 0

	real, allocatable, save :: uov(:), vov(:)
	real, allocatable, save :: unv(:), vnv(:)

	contains

!************************************************************

        subroutine mod_hydro_baro_init(nel)

        integer nel

        if( nel == nel_hydro_baro ) return

        if( nel_hydro_baro > 0 ) then
          deallocate(uov)
          deallocate(vov)

          deallocate(unv)
          deallocate(vnv)
        end if

        nel_hydro_baro = nel

        if( nel == 0 ) return

        allocate(uov(nel))
        allocate(vov(nel))

        allocate(unv(nel))
        allocate(vnv(nel))

        end subroutine mod_hydro_baro_init

!************************************************************

        end module mod_hydro_baro
