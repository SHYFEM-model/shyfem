
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017,2019  Georg Umgiesser
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

! shyelab_average.f: utility for averaging
!
! revision log :
!
! 07.10.2017	ggu	started
! 16.02.2019	ggu	changed VERS_7_5_60
!
!***********************************************************

	subroutine average_vertical_node(lmax,hlv,z,htot,values,aver)

! averages vertically a profile of values

	implicit none

	integer lmax
	real hlv(lmax)
	real z,htot
	real values(lmax)
	real aver		!return

	integer l
	integer nlvaux,nsigma
	real hsigma
	real h
	real hd(lmax)
	double precision vaccum,haccum

	aver = 0.
	if( lmax == 1 ) aver = values(1)
	if( lmax <= 1 ) return

        call get_sigma_info(nlvaux,nsigma,hsigma)
        call get_layer_thickness(lmax,nsigma,hsigma,z,htot,hlv,hd)

	vaccum = 0.
	haccum = 0.
	do l=1,lmax
	  h = hd(l)
	  vaccum = vaccum + values(l) * h
	  haccum = haccum + h
	end do

	aver = vaccum / haccum

	end

!***********************************************************

