
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

!==================================================================
        module mod_geom_dynamic
!==================================================================

        implicit none

	integer, private, save :: nkn_geom_dynamic = 0
        integer, private, save :: nel_geom_dynamic = 0
	
	integer, allocatable, save :: iwegv(:)
	integer, allocatable, save :: iwetv(:)
	integer, allocatable, save :: inodv(:)
	integer, allocatable, save :: inode_static(:)

!==================================================================
	contains
!==================================================================

	subroutine mod_geom_dynamic_init(nkn,nel)

	integer nkn
	integer nel

	if( nkn == nkn_geom_dynamic .and. 
     +  nel == nel_geom_dynamic ) return 

	if( nkn > 0 .or. nel > 0 ) then
          if( nkn == 0 .or. nel == 0 ) then
            write(6,*) 'nkn,nel: ',nkn,nel
            stop 'error stop mod_geom_dynamic_init: 
     +            incompatible parameters'
          end if
        end if

	if( nel_geom_dynamic > 0 ) then
	  deallocate(iwegv)
	  deallocate(iwetv)
          deallocate(inodv)
          deallocate(inode_static)
        end if
	
	nel_geom_dynamic = nel
	nkn_geom_dynamic = nkn	 
	 
	if( nkn == 0 ) return
    	
	allocate(iwegv(nel))
	allocate(iwetv(nel))
	allocate(inodv(nkn))
	allocate(inode_static(nkn))
	
	iwegv = 0
	iwetv = 0
	inodv = 0
	inode_static = 0

       	end subroutine mod_geom_dynamic_init	

!==================================================================
	end module mod_geom_dynamic
!==================================================================


	

	


