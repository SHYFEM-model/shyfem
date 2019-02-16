
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
c statement functions to test for nodes with boundary condition
c----------------------------------------------------------------------
c
c iopbnd(k) = 0		no open BC
c iopbnd(k) > 0		external open BC (ibtyp=1,2)
c iopbnd(k) < 0		internal open BC (ibtyp=3)
c
c----------------------------------------------------------------------

	integer iopbnd(nkndim)
	common /iopbnd/iopbnd
	save /iopbnd/

	integer k_n
	logical is_boundary, is_external_boundary, is_internal_boundary
	logical is_inner

	is_boundary(k_n) = iopbnd(k_n) .ne. 0
	is_external_boundary(k_n) = iopbnd(k_n) .gt. 0
	is_internal_boundary(k_n) = iopbnd(k_n) .lt. 0

	is_inner(k_n) = iopbnd(k_n) .eq. 0

c----------------------------------------------------------------------
c end of file
c----------------------------------------------------------------------

