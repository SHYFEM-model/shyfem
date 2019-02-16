
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

	double precision numv_gotm(0:nlvdim,nkndim)
	double precision nuhv_gotm(0:nlvdim,nkndim)
	double precision tken_gotm(0:nlvdim,nkndim)
	double precision eps_gotm(0:nlvdim,nkndim)
	double precision rls_gotm(0:nlvdim,nkndim)

        common /numv_gotm/ numv_gotm
        common /nuhv_gotm/ nuhv_gotm
        common /tken_gotm/ tken_gotm
        common /eps_gotm/ eps_gotm
        common /rls_gotm/ rls_gotm

	save /numv_gotm/,/nuhv_gotm/
	save /tken_gotm/,/eps_gotm/,/rls_gotm/

        real shearf2(nlvdim,nkndim)
        common /shearf2/shearf2
        real buoyf2(nlvdim,nkndim)
        common /buoyf2/buoyf2
	save /shearf2/,/buoyf2/

