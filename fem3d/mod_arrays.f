
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

! revision log :
!
! 14.03.2019	ggu	written from scratch
! 12.04.2019	ggu	added double precision

!==================================================================
        module arrays
!==================================================================

        implicit none

        INTERFACE init_array
        MODULE PROCEDURE         init_array_i,init_array_in
     +                          ,init_array_r,init_array_rn
     +                          ,init_array_d,init_array_dn
        END INTERFACE

        INTERFACE extend_array
        MODULE PROCEDURE         extend_array_i
     +                          ,extend_array_r
     +                          ,extend_array_d
        END INTERFACE

        INTERFACE append_to_array
        MODULE PROCEDURE         append_to_array_i
     +                          ,append_to_array_r
     +                          ,append_to_array_d
        END INTERFACE

        INTERFACE trim_array
        MODULE PROCEDURE         trim_array_i
     +                          ,trim_array_r
     +                          ,trim_array_d
        END INTERFACE

!==================================================================
        contains
!==================================================================

        subroutine init_array_i(array)

	integer, allocatable :: array(:)

	if( allocated(array) ) deallocate(array)
	allocate(array(10))

        end subroutine init_array_i

!*********************

        subroutine init_array_r(array)

	real, allocatable :: array(:)

	if( allocated(array) ) deallocate(array)
	allocate(array(10))

        end subroutine init_array_r

!*********************

        subroutine init_array_d(array)

	double precision, allocatable :: array(:)

	if( allocated(array) ) deallocate(array)
	allocate(array(10))

        end subroutine init_array_d

!*********************

        subroutine init_array_in(array,ndim)

	integer, allocatable :: array(:)
	integer ndim

	if( allocated(array) ) deallocate(array)
	allocate(array(ndim))

        end subroutine init_array_in

!*********************

        subroutine init_array_rn(array,ndim)

	real, allocatable :: array(:)
	integer ndim

	if( allocated(array) ) deallocate(array)
	allocate(array(ndim))

        end subroutine init_array_rn

!*********************

        subroutine init_array_dn(array,ndim)

	double precision, allocatable :: array(:)
	integer ndim

	if( allocated(array) ) deallocate(array)
	allocate(array(ndim))

        end subroutine init_array_dn

!******************************************************************

        subroutine extend_array_i(array)

	integer, allocatable :: array(:)

	integer :: ndim
	integer, allocatable :: aux(:)

	ndim = 0
	if( allocated(array) ) ndim = size(array)

        if( ndim == 0 ) then
          ndim = 10
          allocate(array(ndim))
          return
        else
          ndim = ndim*2
          allocate(aux(ndim))
          aux(1:ndim/2) = array(1:ndim/2)
          call move_alloc(aux,array)
        end if

        end subroutine extend_array_i

!*********************

        subroutine extend_array_r(array)

	real, allocatable :: array(:)

	integer :: ndim
	real, allocatable :: aux(:)

	ndim = 0
	if( allocated(array) ) ndim = size(array)

        if( ndim == 0 ) then
          ndim = 10
          allocate(array(ndim))
          return
        else
          ndim = ndim*2
          allocate(aux(ndim))
          aux(1:ndim/2) = array(1:ndim/2)
          call move_alloc(aux,array)
        end if

        end subroutine extend_array_r

!*********************

        subroutine extend_array_d(array)

	double precision, allocatable :: array(:)

	integer :: ndim
	double precision, allocatable :: aux(:)

	ndim = 0
	if( allocated(array) ) ndim = size(array)

        if( ndim == 0 ) then
          ndim = 10
          allocate(array(ndim))
          return
        else
          ndim = ndim*2
          allocate(aux(ndim))
          aux(1:ndim/2) = array(1:ndim/2)
          call move_alloc(aux,array)
        end if

        end subroutine extend_array_d

!******************************************************************

        subroutine append_to_array_i(array,ifill,value)

	integer, allocatable :: array(:)
	integer ifill
	integer value

	ifill = ifill + 1
	if( ifill > size(array) ) call extend_array_i(array)
	array(ifill) = value

        end subroutine append_to_array_i

!*********************

        subroutine append_to_array_r(array,ifill,value)

	real, allocatable :: array(:)
	integer ifill
	real value

	ifill = ifill + 1
	if( ifill > size(array) ) call extend_array_r(array)
	array(ifill) = value

        end subroutine append_to_array_r

!*********************

        subroutine append_to_array_d(array,ifill,value)

	double precision, allocatable :: array(:)
	integer ifill
	double precision value

	ifill = ifill + 1
	if( ifill > size(array) ) call extend_array_d(array)
	array(ifill) = value

        end subroutine append_to_array_d

!******************************************************************

        subroutine trim_array_i(array,ifill)

	integer, allocatable :: array(:)
	integer ifill

	integer :: ndim
	integer, allocatable :: aux(:)

        ndim = size(array)
	if( ifill > ndim ) stop 'error stop trim_array: internal error'
	if( ifill == ndim ) return

        allocate(aux(ifill))
        aux(1:ifill) = array(1:ifill)
        call move_alloc(aux,array)

        end subroutine trim_array_i

!*********************

        subroutine trim_array_r(array,ifill)

	real, allocatable :: array(:)
	integer ifill

	integer :: ndim
	real, allocatable :: aux(:)

        ndim = size(array)
	if( ifill > ndim ) stop 'error stop trim_array: internal error'
	if( ifill == ndim ) return

        allocate(aux(ifill))
        aux(1:ifill) = array(1:ifill)
        call move_alloc(aux,array)

        end subroutine trim_array_r

!*********************

        subroutine trim_array_d(array,ifill)

	double precision, allocatable :: array(:)
	integer ifill

	integer :: ndim
	double precision, allocatable :: aux(:)

        ndim = size(array)
	if( ifill > ndim ) stop 'error stop trim_array: internal error'
	if( ifill == ndim ) return

        allocate(aux(ifill))
        aux(1:ifill) = array(1:ifill)
        call move_alloc(aux,array)

        end subroutine trim_array_d

!******************************************************************

!==================================================================
        end module arrays
!==================================================================

