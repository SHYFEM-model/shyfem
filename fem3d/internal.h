
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

        real rdistv(nkndim)
        common /rdistv/rdistv

        real fcorv(neldim)
        common /fcorv/fcorv
        real fxv(nlvdim,neldim)          !new HYDRO debora
        common /fxv/fxv
        real fyv(nlvdim,neldim)
        common /fyv/fyv

	save /rdistv/,/fcorv/,/fxv/,/fyv/

        integer iuvfix(neldim)
        common /iuvfix/iuvfix
	save /iuvfix/

        double precision ddxv(2*nlvdim,neldim)  !ASYM
        double precision ddyv(2*nlvdim,neldim)  !ASYM
        common /ddxv/ddxv, /ddyv/ddyv
        save /ddxv/, /ddyv/

