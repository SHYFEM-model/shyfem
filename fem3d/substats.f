
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

