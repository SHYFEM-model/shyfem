
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

! revision log :
!
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62

!--------------------------------------------------------------------------

c    Red, Green, Blue in [0,1]
c    Hue, Saturation, Intensity in [0,1]
c
c    x = Red-0.5*(Green+Blue)
c    y = 0.866*(Green-Blue)
c
c    Hue = arctan2(x,y)/(2*PI) 
c    Saturation = (x^2+y^2)^0.5
c    Intensity = (Red+Green+Blue)/3
c

c*****************************************************************

	subroutine rgb2hsi(r,g,b,h,s,i)

c computes HSI from RGB

	implicit none

	real r,g,b
	real h,s,i

	real a,pi
	parameter( a = 0.8660254 , pi = 3.1415927 )

	real x,y

	x = r - 0.5 * ( g + b )
	y =       a * ( g - b )

	if( x .eq. 0. .and. y .eq. 0. ) then
	  h = 0.
	else
	  h = atan2(x,y) / ( 2. * pi )
	end if

c	h in [-0.5,0.5] -> bring to [0,1]

	if( h .lt. 0. ) h = h + 1.

	s = sqrt( x*x + y*y )
	i = ( r + g + b ) / 3.

	end

c*****************************************************************

	subroutine hsi2rgb(h,s,i,r,g,b)

c computes RGB from HSI

	implicit none

	real h,s,i
	real r,g,b

	real a,pi
	parameter( a = 0.8660254 , pi = 3.1415927 )

	real x,y

	x = s * sin(h)
	y = s * cos(h)

	if( x .eq. 0. .and. y .eq. 0. ) then
	  h = 0.
	else
	  h = atan2(x,y) / ( 2. * pi )
	end if

c	h in [-0.5,0.5] -> bring to [0,1]

	if( h .lt. 0. ) h = h + 1.

	s = sqrt( x*x + y*y )
	i = ( r + g + b ) / 3.

	end

c*****************************************************************



