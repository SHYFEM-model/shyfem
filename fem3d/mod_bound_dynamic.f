
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2019  Georg Umgiesser
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
! 29.09.2015	ggu	changed VERS_7_2_5
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 11.10.2022	ggu	initialize values

!**************************************************************************

!==================================================================
        module mod_bound_dynamic
!==================================================================

        implicit none

        integer, private, save :: nkn_bound_dynamic = 0
        integer, private, save :: nlv_bound_dynamic = 0

        real, parameter, private :: flag = -9988765.0

        real, allocatable, save :: rzv(:)
        real, allocatable, save :: rqv(:)
        real, allocatable, save :: rqpsv(:)
        real, allocatable, save :: rqdsv(:)
        real, allocatable, save :: mfluxv(:,:)

!==================================================================
        contains
!==================================================================

        subroutine mod_bound_dynamic_init(nkn,nlv)

        integer nkn
        integer nlv

        if( nkn == nkn_bound_dynamic .and. 
     +		nlv == nlv_bound_dynamic ) return

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
	    stop 'error stop mod_bound_dynamic_init: incompatible params'
          end if
        end if

        if( nkn_bound_dynamic > 0 ) then
          deallocate(rzv)
          deallocate(rqv)
          deallocate(rqpsv)
          deallocate(rqdsv)
          deallocate(mfluxv)
        end if

        nkn_bound_dynamic = nkn
        nlv_bound_dynamic = nlv

        if( nkn == 0 ) return

        allocate(rzv(nkn))
        allocate(rqv(nkn))
        allocate(rqpsv(nkn))
        allocate(rqdsv(nkn))
        allocate(mfluxv(nlv,nkn))

	rzv = 0.
	rqv = 0.
	rqpsv = 0.
	rqdsv = 0.
	mfluxv = 0.

        end subroutine mod_bound_dynamic_init

!******************************************************************

        function is_zeta_boundary(k)

        logical is_zeta_boundary
        integer k

        is_zeta_boundary = rzv(k) /= flag

        end function is_zeta_boundary

!==================================================================
        end module mod_bound_dynamic
!==================================================================


