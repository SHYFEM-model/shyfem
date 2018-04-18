
	program nudge

	implicit none

	integer niter
	real c,co,cn,cobs
	real tau,dt,alfa
	real eps,high

	eps = 1.e-4
	high = 1.e+10

	niter = 0
	c = 100
	cobs = 50

	tau = 100
	dt = 750
	alfa = dt/tau

	do
	  niter = niter + 1
	  if( abs(c-cobs) < eps ) exit
	  if( abs(c) > high ) exit
	  !if( abs(cn-co) < eps ) exit
	  !cn = alfa*cobs + (1.-alfa)*c		!explicit
	  cn = (c + alfa*cobs)/(1.+alfa)	!implicit
	  write(6,*) niter,c,cn,cobs
	  c = cn
	end do

	write(6,*) 'program finished'

	end
