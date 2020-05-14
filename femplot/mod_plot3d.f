
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2016,2018-2019  Georg Umgiesser
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
! 17.07.2015	ggu	changed VERS_7_1_52
! 28.04.2016	ggu	changed VERS_7_5_9
! 18.12.2018	ggu	changed VERS_7_5_52
! 21.05.2019	ggu	changed VERS_7_5_62

!==================================================================
	module mod_plot3d
!==================================================================

	implicit none

	integer, private, save  :: nkn_plot3d = 0
	integer, private, save  :: nel_plot3d = 0
	integer, private, save  :: nlv_plot3d = 0
	integer, private, save  :: np_plot3d = 0

	real, allocatable, save :: fvlv(:,:)
	real, allocatable, save :: wauxv(:,:)
	real, allocatable, save :: het3v(:,:)
	real, allocatable, save :: p3(:,:)

!==================================================================
	contains
!==================================================================

	subroutine mod_plot3d_init(nkn,nel,nlv,np)

	integer nkn,nel,nlv,np

        if( nkn == nkn_plot3d .and. nel == nel_plot3d .and.
     +		nlv == nlv_plot3d .and. np == np_plot3d ) return

        if( nel > 0 .or. nkn > 0 .or. nlv > 0 .or. np > 0 ) then
          if( nel == 0 .or. nkn == 0 .or. nlv == 0 .or. np == 0 ) then
            write(6,*) 'nel,nkn,nlv,np: ',nel,nkn,nlv,np
            stop 'error stop mod_plot3d_init: incompatible parameters'
          end if
        end if

        if( nkn_plot3d > 0 ) then
          deallocate(fvlv)
          deallocate(wauxv)
          deallocate(het3v)
          deallocate(p3)
        end if

        nkn_plot3d = nkn
        nel_plot3d = nel
        nlv_plot3d = nlv
        np_plot3d = np

	if( nkn == 0 ) return

        allocate(fvlv(nlv,nkn))
        allocate(wauxv(0:nlv,nkn))
        allocate(het3v(nlv,nel))
        allocate(p3(nlv,np))

	end subroutine mod_plot3d_init

!==================================================================
	end module mod_plot3d
!==================================================================

