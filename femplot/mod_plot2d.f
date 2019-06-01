
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2016,2018  Georg Umgiesser
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
	module mod_plot2d
!==================================================================

	implicit none

	integer, private, save  :: nkn_plot2d = 0
	integer, private, save  :: nel_plot2d = 0
	integer, private, save  :: np_plot2d = 0

	real, allocatable, save :: arfvlv(:)
	real, allocatable, save :: hetv(:)
	real, allocatable, save :: parray(:)

	logical, allocatable, save :: bwater(:)
	logical, allocatable, save :: bkwater(:)

!==================================================================
	contains
!==================================================================

	subroutine mod_plot2d_init(nkn,nel,np)

	integer nkn,nel,np

        if( nkn == nkn_plot2d .and. nel == nel_plot2d 
     +		.and. np == np_plot2d ) return

        if( nel > 0 .or. nkn > 0 .and. np > 0 ) then
          if( nel == 0 .or. nkn == 0 .and. np == 0 ) then
            write(6,*) 'nel,nkn,np: ',nel,nkn,np
            stop 'error stop mod_plot2d_init: incompatible parameters'
          end if
        end if

        if( nkn_plot2d > 0 ) then
          deallocate(arfvlv)
          deallocate(hetv)
          deallocate(parray)
          deallocate(bwater)
          deallocate(bkwater)
        end if

        nel_plot2d = nel
        nkn_plot2d = nkn
        np_plot2d = np

	if( nkn == 0 ) return

        allocate(arfvlv(nkn))
        allocate(hetv(nel))
        allocate(parray(np))
        allocate(bwater(nel))
        allocate(bkwater(nkn))

	end subroutine mod_plot2d_init

!==================================================================
	end module mod_plot2d
!==================================================================

