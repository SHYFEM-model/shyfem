
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2018  Georg Umgiesser
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
! 26.02.2015	ggu	changed VERS_7_1_5
! 18.12.2018	ggu	changed VERS_7_5_52

	logical bdebug_out
	parameter (bdebug_out = .false.)

	integer nunit

	integer iformat,iwave
        common /supout1/ iformat,iwave

	integer nunit_wave,nunit_ous,nunit_nos,nunit_fvl
     +			,nunit_eos,nunit_fem
        common /supout2/ nunit_wave,nunit_ous,nunit_nos,nunit_fvl
     +                  ,nunit_eos,nunit_fem

	real regp(7)
	common /supout3/regp

	save /supout1/,/supout2/,/supout3/

