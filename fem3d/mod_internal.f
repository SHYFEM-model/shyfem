
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013,2015,2019-2020  Georg Umgiesser
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
! 01.01.2013	ggu	written from scratch
! 10.07.2015	ggu	changed VERS_7_1_50
! 18.09.2015	ggu	changed VERS_7_2_3
! 16.02.2019	ggu	changed VERS_7_5_60
! 11.04.2019	ggu	added rcomputev
! 21.05.2019	ggu	changed VERS_7_5_62
! 06.03.2020	ggu	custom routine set_fric_max()
! 26.03.2020	ggu	new variable vis_max and routine set_vis_max()
! 26.05.2020	ggu	rdistv is now defined on elements

!==========================================================================
	module mod_internal
!==========================================================================

	implicit none

        integer, private, save :: nkn_internal = 0
        integer, private, save :: nel_internal = 0
        integer, private, save :: nlv_internal = 0
        
        real, allocatable, save :: rcomputev(:)	!compute all terms
        real, allocatable, save :: rdistv(:)	!compute expl terms by dist
        real, allocatable, save :: fcorv(:)
        real, allocatable, save :: fxv(:,:)
        real, allocatable, save :: fyv(:,:)
        real, allocatable, save :: momentxv(:,:)
        real, allocatable, save :: momentyv(:,:)
        real, allocatable, save :: iuvfix(:)
        double precision, allocatable, save :: ddxv(:,:)
        double precision, allocatable, save :: ddyv(:,:)

        real, save :: rfric_max = 1./1800.     !half hour time scale
        real, save :: vis_max = 10.            !possibly highest value

!==========================================================================
        contains
!==========================================================================

	subroutine mod_internal_init(nkn,nel,nlv)

	integer nkn
        integer nel
        integer nlv

        if( nkn == nkn_internal .and. nel == nel_internal
     +  .and. nlv == nlv_internal) return

        if( nkn > 0 .or. nel > 0 .or. nlv > 0) then
          if( nkn == 0 .or. nel == 0 .or. nlv == 0) then
            write(6,*) 'nkn,nel,nlv: ',nkn,nel,nlv
            stop 'error stop mod_internal_init: incompatible params'
          end if
        end if

	if( nkn_internal > 0 ) then
          deallocate(rcomputev)
          deallocate(rdistv)
          deallocate(fcorv)
          deallocate(fxv)
          deallocate(fyv)
          deallocate(momentxv)
          deallocate(momentyv)
          deallocate(iuvfix)
          deallocate(ddxv)
          deallocate(ddyv)
        end if

        nel_internal = nel
        nkn_internal = nkn 
        nlv_internal = nlv       
        
        if( nkn == 0 ) return
        
        allocate (rcomputev(nel))
        allocate (rdistv(nel))
        allocate (fcorv(nel))
        allocate (fxv(nlv,nel))
        allocate (fyv(nlv,nel))
        allocate (momentxv(nlv,nkn))
        allocate (momentyv(nlv,nkn))
        allocate (iuvfix(nel))
        allocate (ddxv(2*nlv,nel))
        allocate (ddyv(2*nlv,nel))
        
	rcomputev = 1.
	rdistv = 1.

        end subroutine mod_internal_init 

!==========================================================================
        end module mod_internal
!==========================================================================

	subroutine set_fric_max(fmax)

	use mod_internal

	implicit none

	real fmax

	rfric_max = fmax

	end

!**************************************************************************

	subroutine set_vis_max(fmax)

	use mod_internal

	implicit none

	real fmax

	vis_max = fmax

	end

!**************************************************************************

