
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
! 24.07.2015	ggu	changed VERS_7_1_82
! 16.12.2015	ggu	changed VERS_7_3_16
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 13.04.2022	ggu	new array iboundv (indicator of boundary node)
! 26.04.2022	ggu	dxv, dyv eliminated

!**************************************************************************

!==================================================================
	module mod_geom
!==================================================================

	implicit none

	integer, private, save  :: nkn_geom = 0
	integer, private, save  :: nel_geom = 0
	integer, private, save  :: ngr_geom = 0
	integer, private, save  :: nlk_geom = 0

	integer, save :: maxlnk = 0

!	ilinkv(1) == 0
!	ilinkv(k+1) == max entries in lenkv and linkv
!	ibase = ilinkv(k)
!	n = ilinkv(k+1)
!	n - ibase -> nodes around k
!	if( lenkv(n) == 0 ) boundary node

	integer, allocatable, save :: ilinkv(:)  !pointer into arrays below
	integer, allocatable, save :: lenkv(:)   !element numbers
	integer, allocatable, save :: lenkiiv(:) !vertex number of node in elem
	integer, allocatable, save :: linkv(:)   !node numbers

	integer, allocatable, save :: iboundv(:)
	integer, allocatable, save :: ieltv(:,:)
	integer, allocatable, save :: kantv(:,:)

!==================================================================
	contains
!==================================================================

	subroutine mod_geom_init(nkn,nel,ngr)

	integer nkn,nel,ngr

	integer nlk

        if( ngr == ngr_geom .and. nel == nel_geom .and.
     +      nkn == nkn_geom ) return

        if( ngr > 0 .or. nel > 0 .or. nkn > 0 ) then
          if( ngr == 0 .or. nel == 0 .or. nkn == 0 ) then
            write(6,*) 'ngr,nel,nkn: ',ngr,nel,nkn
            stop 'error stop mod_geom_init: incompatible parameters'
          end if
        end if

        if( nkn_geom > 0 ) then
          deallocate(ilinkv)
          deallocate(lenkv)
          deallocate(lenkiiv)
          deallocate(linkv)
          deallocate(iboundv)
          deallocate(ieltv)
          deallocate(kantv)
        end if

	nlk = 3*nel + 2*nkn
	maxlnk = ngr

        ngr_geom = ngr
        nel_geom = nel
        nkn_geom = nkn
        nlk_geom = nlk

	if( nkn == 0 ) return

        allocate(ilinkv(nkn+1))
        allocate(lenkv(nlk))
        allocate(lenkiiv(nlk))
        allocate(linkv(nlk))
        allocate(iboundv(nkn))
        allocate(ieltv(3,nel))
        allocate(kantv(2,nkn))

	end subroutine mod_geom_init

!------------------------------------------------------------------

	pure function is_boundary_node(k)

	logical is_boundary_node
	integer, intent(in) :: k

	is_boundary_node = ( iboundv(k) /= 0 )

	end function is_boundary_node

!------------------------------------------------------------------

	pure function n_nodes_around(k)

	integer n_nodes_around
	integer, intent(in) :: k

	integer n

	n = ilinkv(k+1) - ilinkv(k)
	n_nodes_around = n

	end function n_nodes_around

!------------------------------------------------------------------

	pure function n_elems_around(k)

	integer n_elems_around
	integer, intent(in) :: k

	integer n

	n = ilinkv(k+1) - ilinkv(k)
	if( is_boundary_node(k) ) n = n - 1
	n_elems_around = n

	end function n_elems_around

!==================================================================
	end module mod_geom
!==================================================================

