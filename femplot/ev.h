
c------------------------------------------------
c finite element information for each element
c------------------------------------------------
c
c please change this info also in main

	integer evdim
	parameter ( evdim = 19 )

	!real ev(evdim,1)
	double precision ev(evdim,1)
	common /ev/ev

	save /ev/

c------------------------------------------------

