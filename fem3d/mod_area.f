
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

        module mod_area

        implicit none

        integer, private, save :: nkn_area = 0
        integer, private, save :: nlv_area = 0

        real, allocatable, save :: areakv(:,:)

	contains


!************************************************************

        subroutine mod_area_init(nkn,nlv)

        integer nkn
        integer nlv

        if( nkn == nkn_area .and. nlv == nlv_area) return   !!!!!giusto

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
            stop 'error stop mod_area_init: incompatible parameters'
          end if
        end if

        if( nkn_area > 0 ) then
          deallocate(areakv)
        end if

        nkn_area = nkn
        nlv_area = nlv

        if( nkn == 0 ) return

        allocate(areakv(nlv,nkn))

        end subroutine mod_area_init

!************************************************************

        end module mod_area

