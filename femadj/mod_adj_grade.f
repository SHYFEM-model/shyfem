
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015  Georg Umgiesser
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

!--------------------------------------------------------------------------
!
! revision log :
!
! 30.07.2015    ggu     written
!
!--------------------------------------------------------------------------

!===================================================================
	module mod_adj_grade
!===================================================================

	integer, save :: ngrdi = 0

	integer, save, allocatable :: ngrade(:)
	integer, save, allocatable :: nbound(:)
	integer, save, allocatable :: ngri(:,:)

!===================================================================
	contains
!===================================================================

	subroutine mod_adj_grade_init(nkn,ngr)

	integer nkn,ngr

	ngrdi = ngr

	allocate(ngrade(nkn))
	allocate(nbound(nkn))
	allocate(ngri(2*ngr,nkn))

	end subroutine mod_adj_grade_init

!===================================================================
	end module mod_adj_grade
!===================================================================

! how to mark static nodes (not moveable nodes)
!
!	iastatic gives the node type that has to be considered static
!	if you do not want to use static nodes, set it to -1
!	do not use the value 0
!	if set to a positive values it marks all nodes with this
!	node type as static (no exchange, no smooth)
!	boundary nodes are considered static in any case
!
!	nbstatic is the flag used in nbound to indicate static nodes
!	nbstatic cannot be 0 (internal) or 1 (boundary)
!	you can leave the default value below

!===================================================================
	module mod_adj_static
!===================================================================

        integer, parameter :: iastatic = -1
        integer, parameter :: nbstatic = 4

!===================================================================
	end module mod_adj_static
!===================================================================

