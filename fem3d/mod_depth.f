
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

!====================================================
	module mod_depth
!====================================================

	implicit none

        integer, private, save :: nkn_depth = 0
        integer, private, save :: nel_depth = 0
        
        real, allocatable, save :: hkv(:)
        real, allocatable, save :: hev(:)

        real, allocatable, save :: hkv_min(:)
        real, allocatable, save :: hkv_max(:)

!====================================================
        contains
!====================================================

	subroutine mod_depth_init(nkn,nel)

	integer nkn
        integer nel

        if( nkn == nkn_depth .and. nel == nel_depth ) return

        if( nkn > 0 .or. nel > 0 ) then
          if( nkn == 0 .or. nel == 0 ) then
            write(6,*) 'nkn,nel: ',nkn,nel
            stop 'error stop mod_depth_init: incompatible params'
          end if
        end if

	if( nkn_depth > 0 ) then
          deallocate(hev)
          deallocate(hkv)
          deallocate(hkv_min)
          deallocate(hkv_max)
        end if

        nkn_depth = nkn 
        nel_depth = nel
        
        if( nkn == 0 ) return
        
        allocate(hkv(nkn))
        allocate(hev(nel))
        allocate(hkv_min(nkn))
        allocate(hkv_max(nkn))

        end subroutine mod_depth_init 

!====================================================
        end module mod_depth
!====================================================

