c
c $Id: sigmautil.f,v 1.3 2009-09-14 08:20:58 georg Exp $
c
c sigma utilities for output files
c
c revision log :
c
c 07.11.2011	ggu	layer thickness for hybrid coordinates
c 14.11.2011	ggu	new sigma routines copied to this file
c 02.12.2011	ggu	bug fix in init_sigma_info() for nlv == 1
c
c******************************************************************

	subroutine get_sigma_info(nlv,nsigma,hsigma)

	implicit none

	integer nlv
	integer nsigma
	real hsigma

	integer nlv_com,nsigma_com
	common /nsigma_com/nlv_com,nsigma_com
	real hsigma_com
	common /hsigma_com/hsigma_com
	save /nsigma_com/,/hsigma_com/

	nlv    = nlv_com
	nsigma = nsigma_com
	hsigma = hsigma_com

	end

c******************************************************************

	subroutine set_sigma_info(nlv,nsigma,hsigma)

	implicit none

	integer nlv
	integer nsigma
	real hsigma

	integer nlv_com,nsigma_com
	common /nsigma_com/nlv_com,nsigma_com
	real hsigma_com
	common /hsigma_com/hsigma_com
	save /nsigma_com/,/hsigma_com/

	nlv_com    = nlv
	nsigma_com = nsigma
	hsigma_com = hsigma

	end

c******************************************************************

	subroutine init_sigma_info(nlv,hlv)

	implicit none

	integer nlv
	real hlv(1)

	integer l
	integer nsigma
	real hsigma

c---------------------------------------------------------
c scan depth structure
c---------------------------------------------------------

	hsigma = 10000.
        l = 2                           !HACK for nlv == 1
        if( nlv .eq. 1 ) goto 1

	do l=2,nlv
	  if( hlv(l) .gt. hlv(l-1) ) goto 1
	end do

c---------------------------------------------------------
c only sigma layers found
c---------------------------------------------------------

	if( hlv(nlv) .ne. -1 ) then
          write(6,*) 'nlv,hlv(nlv): ',nlv,hlv(nlv)
	  stop 'error stop init_sigma_info: internal error (1)'
	end if
	nsigma = nlv
	return

c---------------------------------------------------------
c zeta or hybrid levels found
c
c this algorithm cannot handle hybrid levels with only 2 sigma layers
c---------------------------------------------------------

    1	continue
	if( l .eq. 2 ) then	!only zeta levels
	  nsigma = 0
	else			!hybrid levels
	  nsigma = l
	  hsigma = hlv(l)
	end if

c---------------------------------------------------------
c store information in common block
c---------------------------------------------------------

	call set_sigma_info(nlv,nsigma,hsigma)

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c******************************************************************

	subroutine get_layer_thickness(ie,lmax,bzeta,hl)

c returns layer thickness - works also for lmax higher than actual layers
c
c works also for lmax higher than actual layers
c in this case the last values for hl are 0

	implicit none

	integer ie
	integer lmax
	logical bzeta
	real hl(1)

	real hlv(1)
	common /hlv/hlv
	real hev(1)
	common /hev/hev
	real zenv(3,1)
	common /zenv/zenv

	logical bdebug
	integer ii,l
	integer nlv,nsigma
	real hsigma
	real zmed
	real htot,hsig,htop,hbot

	!bdebug = ie .eq. 100
	bdebug = .false.

c---------------------------------------------------------
c get sigma info
c---------------------------------------------------------

	call get_sigma_info(nlv,nsigma,hsigma)

	if( bdebug ) write(6,*) nlv,nsigma,hsigma,ie,lmax,hev(ie)

c---------------------------------------------------------
c compute contribution of zeta
c---------------------------------------------------------

	zmed = 0.
	if( bzeta ) then
          zmed = 0.
          do ii=1,3
            zmed = zmed + zenv(ii,ie)
          end do
          zmed = zmed / 3.
	end if

c---------------------------------------------------------
c compute level structure of sigma levels
c---------------------------------------------------------

	htot = hev(ie)
	hsig = min(htot,hsigma) + zmed

	hbot = 0.
	do l=1,nsigma
	  htop = hbot
	  hbot = hlv(l)
	  if( l .eq. nsigma ) hbot = -1.
	  hl(l) = -hsig * (hbot-htop)
	end do

	if( bdebug ) write(6,*) l,hsig,lmax,nsigma

c---------------------------------------------------------
c compute level structure of zeta and/or hybrid levels
c---------------------------------------------------------

	if( lmax .gt. nsigma ) then		!also zeta coordinates
	  if( lmax .eq. 1 ) then		!just one layer
	    hl(1) = htot + zmed
	  else
	    hbot = hsigma
	    if( nsigma .eq. 0 ) hbot = -zmed
	    do l=nsigma+1,lmax
	      htop = hbot
	      hbot = hlv(l)
	      if( bdebug ) write(6,*) l,htop,hbot,htot
	      if( hbot .gt. htot ) hbot = htot	!last layer
	      hl(l) = hbot - htop
	    end do
	  end if
	end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c******************************************************************

