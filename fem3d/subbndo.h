
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

c----------------------------------------------------------------------
c data structures for open boundary conditions
c----------------------------------------------------------------------

	integer kbcdim
	integer kopdim
	parameter ( kbcdim = nrbdim )	!total number of open boundary nodes
	parameter ( kopdim = ngrdim )	!maximum number of nodes close to OB

	integer nbndo			!total number of OB nodes
	integer ndebug			!unit number for debug messages

	integer nopnod(kbcdim)		!number of internal nodes close to OB
	integer ibcnod(kbcdim)		!number of boundary
	integer kbcnod(kbcdim)		!number of boundary node
	integer itynod(kbcdim)		!type of boundary

	integer nopnodes(kopdim,kbcdim)	!nodes close to OB

	real xynorm(2,kbcdim)		!normal direction for OB node
	real wopnodes(kopdim,kbcdim)	!weights of nodes close to OB

	common /ibndoc/ nbndo,ndebug,nopnod,ibcnod,kbcnod,itynod,nopnodes
	common /rbndoc/ xynorm,wopnodes

	save /ibndoc/, /rbndoc/

c----------------------------------------------------------------------
c next array is global and can be used to check for open boundary nodes
c----------------------------------------------------------------------
c
c       integer iopbnd(nkndim)          !if >0 pointer into array irv
c                                       !if <0 internal boundary (= -ibc)
c
c----------------------------------------------------------------------
c end of file
c----------------------------------------------------------------------

c possible tests:
c
c iopbnd(k) .ne. 0 		boundary node (external or internal)
c iopbnd(k) .gt. 0 		external boundary node (ibtyp = 1,2)
c iopbnd(k) .lt. 0 		internal boundary node (ibtyp = 3)
 
