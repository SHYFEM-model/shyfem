
	program gauss

	implicit none

	real pi
	real mu,sigma
	real a,b
	real val,valtot,x
	real eps,fact
	integer i

	mu = 0.
	sigma = 3.
	pi = 4.*atan(1.)
	eps = 1.e-7

	write(6,*) 'mu,pi,eps : ',mu,pi,eps
	write(6,*) 'Enter sigma : '
	read(5,*) sigma

	a = 1. / ( sqrt(2.*pi) * sigma )
	b = 1. / ( 2 * sigma * sigma )

	fact = 1.	!used to include x=0 only once
	valtot = 0.
	val = 1.
	i = 0

	do while( val .ge. eps )
	  x = i
	  val = a * exp( -(x-mu)**2 * b )
	  valtot = valtot + fact * val
	  write(6,*) x,val,valtot
	  fact = 2.
	  i = i + 1
	end do

	end
