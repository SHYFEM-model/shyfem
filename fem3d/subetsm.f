
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

c ETS file administration routines (deals with section $EXTTS)
c
c revision log :
c
c 24.01.2014	ggu	copied from subexta.f
c 21.05.2015	ggu	changed VERS_7_1_11
c 01.04.2016	ggu	changed VERS_7_5_7
c 16.02.2019	ggu	changed VERS_7_5_60
c
c******************************************************************
c******************************************************************
c******************************************************************

!==================================================================
        module ets
!==================================================================

        implicit none

        integer, save :: nets = 0
        integer, save, allocatable :: nkets(:)		!node numbers
        character*80, save, allocatable :: chets(:)	!description

        real, save, allocatable :: xets(:)		!x-coordinates
        real, save, allocatable :: yets(:)		!y-coordinates

        integer, save, allocatable :: ilets(:)		!layers
        real, save, allocatable :: hets(:)		!depth

        real, save, allocatable :: outets(:,:)		!aux array

        integer, save, allocatable :: il4ets(:)		!layers
        real, save, allocatable :: out4ets(:,:)		!aux array

!==================================================================
        contains
!==================================================================

	subroutine ets_init_module(n)

	integer n

	nets = n

        allocate(nkets(n))
        allocate(chets(n))

        allocate(xets(n))
        allocate(yets(n))
        allocate(ilets(n))
        allocate(il4ets(n))
        allocate(hets(n))

	end subroutine ets_init_module

!==================================================================
        end module ets
!==================================================================

