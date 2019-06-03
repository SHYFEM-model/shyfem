
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

c dummy subroutines
c
c revision log :
c
c 17.06.1998	ggu	nrdskp introduced (needed to skip section)
c 23.03.2010	ggu	changed v6.1.1
c 05.12.2014	ggu	changed VERS_7_0_8
c 16.02.2019	ggu	changed VERS_7_5_60
c
c***********************************************************

	subroutine colrd
	implicit none
	call nrdskp
	end

c***********************************************************

	subroutine legrd
	implicit none
	call nrdskp
	end

c***********************************************************
 
