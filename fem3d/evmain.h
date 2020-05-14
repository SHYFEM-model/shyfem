
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010-2011,2015,2019  Georg Umgiesser
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

c------------------------------------------------
c finite element information for each element
c------------------------------------------------
c
c please change this info also in main

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 17.02.2011	ggu	changed VERS_6_1_18
! 05.06.2015	ggu	changed VERS_7_1_12
! 16.02.2019	ggu	changed VERS_7_5_60

c------------------------------------------------


	integer evdim
	parameter ( evdim = 19 )

	!real ev(evdim,neldim)
	double precision ev(evdim,neldim)
	common /ev/ev

	save /ev/

c------------------------------------------------

