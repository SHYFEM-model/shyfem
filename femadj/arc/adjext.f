c
c $Id: adjext.f,v 1.3 2007-03-20 13:19:42 georg Exp $
c
c***********************************************************

	subroutine memnod(x,y)

c saves deleted node

	implicit none

	include 'param.h'

	real x,y

	integer next
	common /nkone/ next
	real xegv(nkndim), yegv(nkndim)
	common /xegv/xegv, /yegv/yegv

	integer icall	!FIXME
	save icall
	data icall / 0 /

	if( icall .eq. 0 ) then	!FIXME
	  icall = 1
	  next = 0
	end if

	next = next + 1
	if( next .gt. nkndim ) then
	  stop 'error stop memnod: dimension nkndim'
	end if

	xegv(next) = x
	yegv(next) = y

	end

c***********************************************************

	subroutine primem(nbase)

c prints deleted nodes

	implicit none

	include 'param.h'

	integer nbase

	integer next
	common /nkone/ next
	real xegv(nkndim), yegv(nkndim)
	common /xegv/xegv, /yegv/yegv

	integer k,n,it

	next = 0	!do not print nodes

	do k=1,next
	    n = nbase + k
	    it = 65
	    write(1,'(i1,2i8,2f14.4)') 1,n,it,xegv(k),yegv(k)
	end do

	end

c***********************************************************

c	block data
c		?? initialize next = 0		!FIXME
c	end

c***********************************************************

