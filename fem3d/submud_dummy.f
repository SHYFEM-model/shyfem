
!--------------------------------------------------------------------------
!
!    Copyright (C) 2012-2019  Georg Umgiesser
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
! 26.06.2012	ggu	changed VERS_6_1_55
! 25.10.2013	ggu	changed VERS_6_1_68
! 05.12.2013	ggu	changed VERS_6_1_70
! 30.05.2014	ggu	changed VERS_6_1_76
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.06.2015	ggu	changed VERS_7_1_12
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 12.10.2015	ggu	changed VERS_7_3_3
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

c*******************************************************************
c
c if the mud module is not needed we do the minimum needed
c
c we have to make sure that vts is initialized -> submud.f
c
c questions: 
c	do we really need an extra array vts?
c	can we not use the original array visv for this?
c	what about diffusivity?
c
c*******************************************************************

	subroutine readmud

	use mod_fluidmud
	use basin
	use levels

	implicit none

	call nrdskp

	end

c*******************************************************************

	subroutine submud

	implicit none

	end

c*******************************************************************

	subroutine submud_init

	use mod_fluidmud
	use basin
	use levels

	implicit none

	call mod_fluidmud_dummy_init(nkn,nlv)
	vts = 0.

	end

c*******************************************************************

	subroutine set_mud_roughness(k,l,alpha)

	implicit none

	integer k,l
	real alpha

	alpha = 1.

	end

c*******************************************************************

	subroutine set_rhomud(k,l,rhop)

	implicit none

	integer k,l
	real rhop

	end

c*******************************************************************

