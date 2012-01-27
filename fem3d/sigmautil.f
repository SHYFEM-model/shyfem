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
c 16.12.2011	ggu	check for non-initialized data structure (blockdata)
c 19.12.2011	ggu	bug fix in init_sigma_info(): call set_sigma_info()
c 27.01.2012	deb&ggu	changes to get_layer_thickness()
c 27.01.2012	deb&ggu	new routine compute_sigma_info()
c
c******************************************************************

	blockdata sigma_out

	implicit none

	integer nlv_com,nsigma_com
	common /nsigma_com/nlv_com,nsigma_com
	real hsigma_com
	common /hsigma_com/hsigma_com
	save /nsigma_com/,/hsigma_com/

	data nlv_com,nsigma_com /-1,-1/
	data hsigma_com /10000./

	end

c******************************************************************

	subroutine check_sigma_initialized

	implicit none

	integer nlv_com,nsigma_com
	common /nsigma_com/nlv_com,nsigma_com
	real hsigma_com
	common /hsigma_com/hsigma_com
	save /nsigma_com/,/hsigma_com/

	if( nlv_com .le. 0 ) then
	  write(6,*) 'nlv_com: ',nlv_com
	  stop 'error stop check_sigma_initialized: not initialized'
	end if

	end

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

	call check_sigma_initialized

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

	integer nsigma
	real hsigma

	call compute_sigma_info(nlv,hlv,nsigma,hsigma)
	call set_sigma_info(nlv,nsigma,hsigma)

	end

c******************************************************************

	subroutine compute_sigma_info(nlv,hlv,nsigma,hsigma)

	implicit none

	integer nlv
	real hlv(1)
	integer nsigma		!number of sigma layers (return)
	real hsigma		!depth of hybrid level (return)

	integer l

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
c end of routine
c---------------------------------------------------------

	end

c******************************************************************

	subroutine get_layer_thickness(ie,lmax,bzeta,nsigma,hsigma,hl)

c returns layer thickness - works also for lmax higher than actual layers
c
c works also for lmax higher than actual layers
c in this case the last values for hl are 0

	implicit none

	integer ie
	integer lmax
	logical bzeta
	integer nsigma
	real hsigma
	real hl(1)

	real hlv(1)
	common /hlv/hlv
	real hev(1)
	common /hev/hev
	real zenv(3,1)
	common /zenv/zenv

	logical bdebug
	integer ii,l
	real zmed
	real htot,hsig,htop,hbot

	bdebug = .false.

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
	    if( bdebug ) write(6,*) nsigma,lmax,zmed,hbot
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

