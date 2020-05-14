
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2017,2019  Georg Umgiesser
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
! 05.12.2017	ggu	changed VERS_7_5_39
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

	module mod_bnd_aux

	implicit none

	integer, private, save :: nkn_bnd_aux = 0
        integer, private, save :: nel_bnd_aux = 0
        
        real, allocatable, save :: rvv(:)	!momentum input (2D)
        real, allocatable, save :: ruv(:)
        real, allocatable, save :: crad(:)	!radiation condition

	contains

*******************************************************************

	subroutine mod_bnd_aux_init(nkn,nel)

	integer nkn
        integer nel

        if( nkn == nkn_bnd_aux .and. nel == nel_bnd_aux ) return

        if( nkn > 0 .or. nel > 0 ) then
          if( nkn == 0 .or. nel == 0 ) then
            write(6,*) 'nkn,nel: ',nkn,nel
            stop 'error stop mod_bnd_aux_init: incompatible params'
          end if
        end if

	if( nkn_bnd_aux > 0 ) then
          deallocate(ruv)
          deallocate(rvv)
          deallocate(crad)
        end if

        nel_bnd_aux = nel
        nkn_bnd_aux = nkn        
        
        if( nkn == 0 .or. nel == 0 ) return
        
        allocate(ruv(nkn))
        allocate(rvv(nkn))
        allocate(crad(nel))
        
	ruv = 0.
	rvv = 0.
	crad = 0.

        end subroutine mod_bnd_aux_init 

!*****************************************************

        end module mod_bnd_aux

