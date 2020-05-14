
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2016,2019  Georg Umgiesser
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
! 10.07.2015	ggu	changed VERS_7_1_50
! 28.04.2016	ggu	changed VERS_7_5_9
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

        module mod_sinking

        implicit none

        integer, private, save  :: nkn_sink = 0
        integer, private, save  :: nlv_sink = 0

        real, allocatable, save :: wsinkv(:,:)       ! sinking velocity [m/s]

        contains

!************************************************************

        subroutine mod_sinking_init(nkn,nlv)

        integer  :: nkn
        integer  :: nlv

        if( nkn == nkn_sink .and. nlv == nlv_sink ) return

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
            stop 'error stop mod_sinking_init: incompatible parameters'
          end if
        end if

        if( nkn_sink > 0 ) then
	  deallocate(wsinkv)
	end if

        nkn_sink = nkn
        nlv_sink = nlv

        if( nkn == 0 ) return

        allocate(wsinkv(0:nlv,nkn))

        end subroutine mod_sinking_init

!************************************************************

        end module mod_sinking


