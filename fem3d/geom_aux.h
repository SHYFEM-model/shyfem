
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

!        integer ilinkv(nkndim+1)
!        common /ilinkv/ilinkv

!        integer lenkv(nlkdim)
!        common /lenkv/lenkv
!        integer linkv(nlkdim)
!        common /linkv/linkv

! revision log :
!
! 19.01.2015	ggu	changed VERS_7_1_2
! 16.02.2019	ggu	changed VERS_7_5_60

        integer ieltv(3,neldim)
        common /ieltv/ieltv

        integer kantv(2,nkndim)
        common /kantv/kantv

        real dxv(nkndim), dyv(nkndim)
        common /dxv/dxv, /dyv/dyv

!        save /ilinkv/,/lenkv/,/linkv/
        save /ieltv/,/kantv/,/dxv/,/dyv/

