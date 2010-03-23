
c creates circle with random fluctuation

	implicit none

	integer ndim
	parameter(ndim = 360)
	real x(ndim)
	real y(ndim)

	real pi,rad,alpha
	real fluct,radius
	real r
	integer i,iseed

	real grand

	radius = 1.
	fluct = 0.2

	pi = 4. * atan(1.)
	rad = pi / 180.
	alpha = rad * 360./ndim
	iseed = 0

	do i=1,ndim
	  r = radius + fluct * radius * ( grand(iseed) - 0.5 )
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

