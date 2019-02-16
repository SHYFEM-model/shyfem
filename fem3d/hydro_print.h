
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

        real uprv(nlvdim,nkndim)
        common /uprv/uprv
        real vprv(nlvdim,nkndim)
        common /vprv/vprv
        real upro(nlvdim,nkndim)
        common /upro/upro
        real vpro(nlvdim,nkndim)
        common /vpro/vpro
        real wprv(0:nlvdim,nkndim)
        common /wprv/wprv

	save /uprv/,/vprv/,/upro/,/vpro/,/wprv/

        real up0v(nkndim)
        common /up0v/up0v
        real vp0v(nkndim)
        common /vp0v/vp0v

	save /up0v/,/vp0v/

        real xv(3,nkndim)
        common /xv/xv

	save /xv/

