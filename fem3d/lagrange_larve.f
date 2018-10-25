c
c $Id: lagrange_larve.f,v 1.2 2008-11-03 10:42:26 georg Exp $
c
c simulates continuous release over open boundaries
c
c revision log :
c
c 12.12.2007    ggu	written from scratch
c
c*******************************************************************

	subroutine lgr_larvae

c manages larvae

	use mod_lagrange

	implicit none

	include 'param.h'

	integer i,ie,ii
	real rlinit,x,y,z,rl
        double precision xx,yy
        double precision xi(3)

	real getpar

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. 0 ) then
	  rlinit = 0.
	  do i=1,nbdy
	    lgr_ar(i)%custom(1) = rlinit
	  end do
	end if

	icall = icall + 1

	do i=1,nbdy
          do ii=1,3
            xi(ii) = lgr_ar(i)%actual%xi(ii)
          end do
	  ie = lgr_ar(i)%actual%ie
          z  = lgr_ar(i)%actual%z	!rel. depth   0=surface  1=bottom
	  rl = lgr_ar(i)%custom(1)
          call xi2xy(abs(ie),xx,yy,xi)
          x   = xx
          y   = yy
	  
	  if( ie .gt. 0 ) then
	    call treat_larva(x,y,z,ie,rl)
	  end if

	  lgr_ar(i)%custom(1) = rl
          lgr_ar(i)%actual%z = z
	end do

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*******************************************************************

	subroutine treat_larva(x,y,z,ie,rl)

	implicit none

	real x,y,z
	integer ie
	real rl,r,perc

	real ggrand

	perc = 0.1
	perc = 0.
	r = ggrand(0)

	if( r .lt. perc ) z = 1.0

	end

c*******************************************************************

