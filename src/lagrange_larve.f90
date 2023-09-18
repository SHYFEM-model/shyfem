!
! $Id: lagrange_larve.f,v 1.2 2008-11-03 10:42:26 georg Exp $
!
! simulates continuous release over open boundaries
!
! revision log :
!
! 12.12.2007    ggu	written from scratch
!
!*******************************************************************
!-------------------------------------------------------------------------
        module lagrange_larve
!-------------------------------------------------------------------------
        contains
!-------------------------------------------------------------------------

	subroutine lgr_larvae(it)

! manages larvae

	use lagrange_data

	implicit none

	integer it

	include 'param.h'

	integer i,ie
	double precision rlinit,x,y,z,rl

	double precision getpar

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. 0 ) then
	  rlinit = 0.
	  do i=1,nbdy
	    lgr_ar(i)%c = rlinit
	  end do
	end if

	icall = icall + 1

	do i=1,nbdy
          x  = lgr_ar(i)%x
          y  = lgr_ar(i)%y
	  ie = lgr_ar(i)%ie
          z  = lgr_ar(i)%z		!rel. depth   0=surface  1=bottom
	  rl = lgr_ar(i)%c
	  
	  if( ie .gt. 0 ) then
	    call treat_larva(x,y,z,ie,rl)
	  end if

	  lgr_ar(i)%c = rl
          lgr_ar(i)%z = z
	end do

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!*******************************************************************

	subroutine treat_larva(x,y,z,ie,rl)

        use random_gen

	implicit none

	double precision x,y,z
	integer ie
	double precision rl,r,perc

	perc = 0.1
	perc = 0.
	r = ggrand(0)

	if( r .lt. perc ) z = 1.0

	end

!*******************************************************************

!-------------------------------------------------------------------------
        end module lagrange_larve
!-------------------------------------------------------------------------
