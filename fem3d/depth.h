
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
! 19.01.2015	ggu	changed VERS_7_1_3
! 16.02.2019	ggu	changed VERS_7_5_60

        real hkv(nkndim)
        common /hkv/hkv
        real hev(neldim)
        common /hev/hev

        real hdknv(nlvdim,nkndim)
        common /hdknv/hdknv
        real hdkov(nlvdim,nkndim)
        common /hdkov/hdkov

        real hdenv(nlvdim,neldim)
        common /hdenv/hdenv
        real hdeov(nlvdim,neldim)
        common /hdeov/hdeov

	save /hkv/,/hev/,/hdknv/,/hdkov/,/hdenv/,/hdeov/

	real hkv_min(nkndim), hkv_max(nkndim)
	common /hkv_min/hkv_min, /hkv_max/hkv_max
	save /hkv_min/,/hkv_max/

