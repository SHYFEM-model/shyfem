
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
! 09.05.2023    lrp     introduce top layer index variable

        integer ilhkv(nkndim)
        common /ilhkv/ilhkv

        integer ilhv(neldim)
        common /ilhv/ilhv

        integer jlhkv(nkndim)
        common /jlhkv/jlhkv

        integer jlhv(neldim)
        common /jlhv/jlhv

        real hlv(nlvdim), hldv(nlvdim)
        common /hlv/hlv, /hldv/hldv

        integer ilmv(neldim)
        common /ilmv/ilmv

        integer ilmkv(nkndim)
        common /ilmkv/ilmkv

	save /ilhkv/,/ilhv/,/hlv/,/hldv/,/ilmv/,/ilmkv/
        save /jlhkv/,/jlhv/
