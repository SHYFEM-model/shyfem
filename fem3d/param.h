
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

!------------------------------------------------------------------------
! param.h - parameter file for SHYFEM
!------------------------------------------------------------------------

	integer nkndim			!maximum number of nodes
	integer neldim			!maximum number of elements
	integer nlvdim			!maximum number of vertical levels

	integer mbwdim			!maximum bandwidth
	integer ngrdim			!maximum grade of nodes

	integer nbcdim			!maximum number of open boundaries
	integer nrbdim			!maximum number of open boundary nodes
	integer nb3dim			!maximum storage for boundary info

	integer nardim			!maximum area code
	integer nexdim			!maximum number of extra points
	integer nfxdim			!maximum number of flux points
	integer ncsdim			!maximum concentration variables

	integer nbdydim			!maximum particles for lagrange model

	parameter ( nkndim = 1 )
	parameter ( neldim = 1 )
	parameter ( nlvdim = 1 )

	parameter ( mbwdim = 1 )
	parameter ( ngrdim = 1 )

	parameter ( nbcdim = 1 )
	parameter ( nrbdim = 1 )
	parameter ( nb3dim = 1 )

	parameter ( nardim = 1 )
	parameter ( nexdim = 1 )
	parameter ( nfxdim = 1 )
	parameter ( ncsdim = 1 )

	parameter ( nbdydim = 100000 )

!------------------------------------------------------------------------
! do not change anything beyond this line
!------------------------------------------------------------------------

	integer nlkdim			!dimension for side index
        parameter (nlkdim=3*neldim+2*nkndim)

!------------------------------------------------------------------------
! end of parameter file
!------------------------------------------------------------------------

