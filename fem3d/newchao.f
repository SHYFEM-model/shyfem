c
c $Id: newchao.f,v 1.1 2001/11/16 07:35:43 georg Exp $
c
c explicitly sets the velocities (transports) to prescribed value
c
c revision log :
c
c 06.06.1996    ggu     written (from sp159f)
c 09.11.2001    ggu     compiler directives removed
c
c******************************************************************

	subroutine chao

c written on 06.06.96 by ggu   (from sp159f)

	implicit none

c parameter
	include 'param.h'
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv
	include 'femtime.h'
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real unv(1),vnv(1)
	real utlnv(nlvdim,1),vtlnv(nlvdim,1)
	integer ipev(1)
	common /unv/unv, /vnv/vnv
	common /utlnv/utlnv, /vtlnv/vtlnv
	common /ipev/ipev

	real umax,dz,fact
	integer ie,l,last,iex
	integer ichange

	ichange = 0
	umax = 0.1
	dz = 3.

	do ie=1,nel

	  iex = ipev(ie)
	  last = mod(iex,100)

	  if( last .le. 2 ) then
	    ichange = ichange + 1
	    do l=1,nlv
		fact = - 0.5
		if( l .eq. 1 ) fact = 1.
		if( l .eq. 2 ) fact = 0.5
		utlnv(l,ie) = umax * dz * fact
		vtlnv(l,ie) = 0.
	    end do
	    unv(ie) = 0.
	    vnv(ie) = 0.
	  end if

	end do

	write(6,*) 'chao: ',umax,nlv,ichange

	return
	end

