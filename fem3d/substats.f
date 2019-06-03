
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

! revision log :
!
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**********************************************************************

	subroutine least_square_regression(n,x,y,b0,b1)

! least square regression - returns b0, b1 for y = b0 + b1 * x

	implicit none

	integer n
	real x(n),y(n)
	real b0,b1

	integer i
	real xm,ym
	double precision sxy,sxx,dx,dy

	b0 = 0.
	b1 = 0.

	if( n == 0 ) return

	xm = sum(x)/n
	ym = sum(y)/n

	sxy = 0.
	sxx = 0.

	do i=1,n
	  dx = (x(i)-xm)
	  dy = (y(i)-ym)
	  sxy = sxy + dx*dy
	  sxx = sxx + dx*dx
	end do

	b1 = sxy / sxx
	b0 = ym - b1 * xm

	end

!**********************************************************************

