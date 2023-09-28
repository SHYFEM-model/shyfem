
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2012,2014,2019  Georg Umgiesser
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
! 23.03.2010	ggu	changed v6.1.1
! 16.12.2010	ggu	changed VERS_6_1_15
! 24.01.2012	ggu	changed VERS_6_1_41
! 18.06.2014	ggu	changed VERS_6_1_77
! 05.11.2014	ggu	changed VERS_7_0_5
! 16.02.2019	ggu	changed VERS_7_5_60

        integer ielfix(0:3,neldim)
        common /ielfix/ielfix
        save /ielfix/

        integer tnudgev(neldim)
        common /tnudgev/tnudgev
        save /tnudgev/

        real ubound(nlvdim,nkndim)
        real vbound(nlvdim,nkndim)
        common /ubound/ubound
        common /vbound/vbound
	save /ubound/,/vbound/
