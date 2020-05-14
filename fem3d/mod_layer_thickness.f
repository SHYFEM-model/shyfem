
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
! 24.07.2015	ggu	changed VERS_7_1_82
! 09.09.2016	ggu	changed VERS_7_5_17
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

	module mod_layer_thickness

	implicit none

        integer, private, save :: nkn_layer_thickness = 0
        integer, private, save :: nel_layer_thickness = 0
        integer, private, save :: nlv_layer_thickness = 0
        
        real, allocatable, save :: hdknv(:,:)
        real, allocatable, save :: hdkov(:,:)
        real, allocatable, save :: hdenv(:,:)
        real, allocatable, save :: hdeov(:,:)

        contains

*******************************************************************

	subroutine mod_layer_thickness_init(nkn,nel,nlv)

	integer nkn
        integer nel
        integer nlv

        if( nkn == nkn_layer_thickness .and. nel == nel_layer_thickness
     +  .and. nlv == nlv_layer_thickness ) return

        if( nkn > 0 .or. nel > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nel == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nel,nlv: ',nkn,nel,nlv
            stop 'error stop mod_layer_thickness_init: '//
     +				'incompatible params'
          end if
        end if

	if( nkn_layer_thickness > 0 ) then
          deallocate(hdknv)
          deallocate(hdkov)
          deallocate(hdenv)
          deallocate(hdeov)
        end if

        nkn_layer_thickness = nkn 
        nel_layer_thickness = nel
        nlv_layer_thickness = nlv       
        
        if( nkn == 0 ) return
        
        allocate(hdknv(nlv,nkn))
        allocate(hdkov(nlv,nkn))
        allocate(hdenv(nlv,nel))
        allocate(hdeov(nlv,nel))

	hdknv = 0.
	hdkov = 0.
	hdenv = 0.
	hdeov = 0.

        end subroutine mod_layer_thickness_init 

!*****************************************************

	function mod_layer_thickness_is_initialized()

	logical mod_layer_thickness_is_initialized

	mod_layer_thickness_is_initialized = nkn_layer_thickness > 0

	end function mod_layer_thickness_is_initialized

!*****************************************************

        end module mod_layer_thickness

