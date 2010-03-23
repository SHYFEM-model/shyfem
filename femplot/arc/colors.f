c
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



