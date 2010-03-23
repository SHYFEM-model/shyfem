
	f = 1.

	do i=1,10
	  f = f * 1.5
	  rfl = roundm(f,-1)
	  rfh = roundm(f,+1)
	  write(6,*) i,f,rfl,rfh
	  rf = divdist(f,7,0)
	end do

	end


c***************************************************
c***************************************************
c***************************************************

	real function roundm(r,mode)

c rounds r to the closest given value
c
c r		value to round
c mode
c		 0: do not round, return r
c		 1: round to higher value
c		-1: round to lower value

	implicit none

	integer ndim
	parameter (ndim = 4)
	real r
	integer mode

	double precision raux,fact,sign,rr
	integer i

	real rmaster(ndim),eps
	save rmaster,eps
	data eps /1.e-5/
	data rmaster /1.,2.,5.,10./
c	data rmaster /1.,2.,4.,5.,8.,10./

	roundm = r
	if( mode .eq. 0 ) return

	fact = 1.0d+0
	if( r .lt. 0. ) then
	  sign = -1.
	  rr = -r
	else
	  sign = 1.
	  rr = r
	end if

	do while( rr*fact .gt. 10. )
	  fact = 0.1d+0 * fact
	end do
	do while( rr*fact .lt. 1. )
	  fact = 10.d+0 * fact
	end do

	rr = rr * fact		!rr is between 1. and 10.

	if( mode .gt. 0 ) then
	  do i=ndim,1,-1
	    if( rr .le. rmaster(i) + eps ) raux = rmaster(i)
	  end do
	else if( mode .lt. 0 ) then
	  do i=1,ndim
	    if( rr .ge. rmaster(i) - eps ) raux = rmaster(i)
	  end do
	end if

	roundm = sign * raux / fact

	end

c****************************************************

	real function rround(r,rmaster,mode)

c rounds r to next rmaster value
c
c r		value to round
c rmaster	value to which r is rounded (must be positive)
c mode
c		 0: do not round, return r
c		 1: round to higher value
c		-1: round to lower value
c
c negative values are respected

	implicit none

	real r,rmaster
	integer mode

	integer iaux
	real raux

	iaux = r/rmaster
	raux = iaux * rmaster

	if( mode .gt. 0 ) then
	  if( raux .lt. r ) raux = raux + rmaster
	else if( mode .lt. 0 ) then
	  if( raux .gt. r ) raux = raux - rmaster
	else
	  raux = r
	end if

	rround = raux

	return
	end

c****************************************************

        real function rdist(xmin,ymin,xmax,ymax)

c computes gridspacing for frame (4-7 grid lines)
c
c xmin,ymin     coordinates of lower left point
c xmax,ymax     coordinates of upper rigth point

	implicit none

	real xmin,ymin,xmax,ymax

	real xdist,ydist,dist,fdist
	integer istell,lines

        xdist=xmax-xmin
        ydist=ymax-ymin

        if( xdist .gt. ydist ) then
		dist = xdist
	else
		dist = ydist
	end if

        dist=nint(dist)

        istell=log10(dist)
        fdist=10**istell
        lines=dist/fdist

        if(lines.le.3) fdist=fdist*0.5
        if(lines.ge.8) fdist=fdist*2.

        rdist=fdist

        return
	end

c**************************************************************

	real function divdist(x,n,mode)

c divides x into n similar pieces (with rounding)

	implicit none

	real x		!length to be divided
	integer n	!number of pieces to create
	integer mode	!0: get closest to n -1:get less  +1: get more 

	integer lines,lines_high,lines_low
	real dist_high,dist_low
	real dist,fdist

	real roundm

        dist=nint(x)

	dist_high = roundm(dist/n,+1)
	dist_low  = roundm(dist/n,-1)

	lines_high = dist / dist_low
	lines_low  = dist / dist_high

	if( mode .lt. 0 ) then
	  fdist = dist_high
	else if( mode .gt. 0 ) then
	  fdist = dist_low
	else
	  if( lines_high-n .lt. n-lines_low ) then
	    fdist = dist_low
	  else
	    fdist = dist_high
	  end if
	end if

	write(6,*) '-------------------'
	write(6,*) x,dist,n,mode
	write(6,*) lines_low,lines_high
	write(6,*) dist_high,dist_low
	write(6,*) fdist,int(dist/fdist)
	write(6,*) '-------------------'

        divdist = fdist

	end

c**************************************************************

