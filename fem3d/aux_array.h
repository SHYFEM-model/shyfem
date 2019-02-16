
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

        real v1v(nkndim)
        common /v1v/v1v
        real v2v(nkndim)
        common /v2v/v2v
        real v3v(nkndim)
        common /v3v/v3v
        real ve1v(neldim)
        common /ve1v/ve1v
        real saux1(nlvdim,nkndim)
        common /saux1/saux1
        real saux2(nlvdim,nkndim)
        common /saux2/saux2
        real saux3(nlvdim,nkndim)
        common /saux3/saux3
        real saux4(nlvdim,nkndim)
        common /saux4/saux4
        real sauxe1(nlvdim,neldim)
        common /sauxe1/sauxe1
        real sauxe2(nlvdim,neldim)
        common /sauxe2/sauxe2

	save /v1v/,/v2v/,/v3v/,/ve1v/
	save /saux1/,/saux2/,/saux3/,/saux4/
	save /sauxe1/,/sauxe2/

