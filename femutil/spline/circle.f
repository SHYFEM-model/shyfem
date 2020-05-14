
!--------------------------------------------------------------------------
!
!    Copyright (C) 2019  Georg Umgiesser
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

c creates circle with random fluctuation

! revision log :
!
! 13.03.2019	ggu	changed VERS_7_5_61

	implicit none

	integer ndim
	parameter(ndim = 360)
	real x(ndim)
	real y(ndim)

	real pi,rad,alpha
	real fluct,radius
	real r,rand
	integer i,iseed

	real grand

	radius = 1.
	fluct = 0.2

	pi = 4. * atan(1.)
	rad = pi / 180.
	alpha = rad * 360./ndim
	iseed = 0

	do i=1,ndim
	  call random_number(rand)
	  r = radius + fluct * radius * ( rand - 0.5 )
	  x(i) = r * cos( i * alpha )
	  y(i) = r * sin( i * alpha )
	  write(6,1) 1,i,0,x(i),y(i)
	end do

	write(6,3) 3,1,0,ndim+1
	write(6,*) (i,i=1,ndim),1

	stop
    1	format(i1,2i8,2f12.4)
    3	format(i1,3i8)
	end

