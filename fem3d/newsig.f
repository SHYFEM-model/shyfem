c
c routines for sigma levels
c
c revision log :
c
c 16.12.2010    ggu     program partially finished
c 19.09.2011    ggu     new routine set_bsigma()
c 04.11.2011    ggu     new routines for hybrid levels
c
c notes :
c
c important files where sigma levels are explicitly needed:
c
c	newini.f		set up of structure
c	newbcl.f		for computation of rho
c	newexpl.f		for baroclinic term
c	subele.f		set new layer thickness
c
c	lagrange_flux.f		limit zeta layers to surface layer
c
c********************************************************************

	subroutine set_sigma_hkhe

c this is called in newini.f

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v
        real hm3v(3,1)
        common /hm3v/hm3v
        real hkv(1)
        common /hkv/hkv
        real hev(1)
        common /hev/hev
        real v1v(1)
        common /v1v/v1v

	logical berror
	integer k,ie,ii
	integer nsigma
	integer ihmin,ihmax
	real hsigma,ha
	real flag,h,hm

c-------------------------------------------------------
c initialize with flag
c-------------------------------------------------------

	call get_sigma(nsigma,hsigma)

	berror = .false.
	flag = -999.

	do k=1,nkn
	  hkv(k) = flag
	  v1v(k) = 0.
	end do

c-------------------------------------------------------
c check hm3v
c-------------------------------------------------------

	do ie=1,nel
	  ihmin = 0
	  ihmax = 0
	  do ii=1,3
	    h = hm3v(ii,ie)
	    if( h .lt. hsigma ) then
	      ihmin = ihmin + 1 
	    else if( h .gt. hsigma ) then
	      ihmax = ihmax + 1 
	    end if
	  end do
	  if( ihmin .gt. 0 .and. ihmax .gt. 0 ) then
	    write(6,*) ie,(hm3v(ii,ie),ii=1,3)
	    berror = .true.
	  end if
	end do

	if( berror ) stop 'error stop set_sigma_hkhe: hsigma crossing'

c-------------------------------------------------------
c set hkv and hev
c-------------------------------------------------------

	do ie=1,nel
	  hm = 0.
	  do ii=1,3
	    h = hm3v(ii,ie)
	    k = nen3v(ii,ie)
	    if( hkv(k) .eq. flag ) then
	      hkv(k) = h
	      v1v(k) = 1.
	    else
	      ha = hkv(k)/v1v(k)
	      if( ha .ge. hsigma .and. h .ge. hsigma ) then
		hkv(k) = hkv(k) + h
		v1v(k) = v1v(k) + 1.
	      else if( ha .ne. h ) then
	        write(6,*) 'different depth: ',k,ha,h
	        berror = .true.
	      end if
	    end if
	    hm = hm + h
	  end do
	  hev(ie) = hm / 3.
	end do

	if( berror ) stop 'error stop set_sigma_hkhe: depth values'

c-------------------------------------------------------
c check hkv
c-------------------------------------------------------

	do k=1,nkn
	  if( hkv(k) .eq. flag ) then
	    write(6,*) 'flag in depth: ',k
	    berror = .true.
	  else
	    hkv(k) = hkv(k)/v1v(k)
	  end if
	end do

	if( berror ) stop 'error stop set_sigma_hkhe: flag in depth'

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c********************************************************************

	subroutine get_bsigma(bsigma)

c returns bsigma which is true if sigma layers are used

	implicit none

	logical bsigma

        real hldv(1)
        common /hldv/hldv

	bsigma = hldv(1) .lt. 0.

	end

c********************************************************************

	subroutine get_sigma(nsigma,hsigma)

	implicit none

	integer nsigma
	real hsigma

	real getpar

	nsigma = nint(getpar('nsigma'))
	hsigma = getpar('hsigma')

	end

c********************************************************************

	subroutine set_sigma(nsigma,hsigma)

	implicit none

	integer nsigma
	real hsigma

	real getpar

	call putpar('nsigma',float(nsigma))
	call putpar('hsigma',hsigma)

	end 

c********************************************************************

	subroutine make_sigma_levels(nsigma)

	implicit none

	integer nsigma

        real hlv(1)
        common /hlv/hlv

	integer l
	real hl

	if( nsigma .le. 0 ) stop 'error stop make_sigma_levels'

        hl = -1. / nsigma
        do l=1,nsigma
          hlv(l) = l * hl
        end do

	end

c********************************************************************

	subroutine make_zeta_levels(lmin,hmin,dzreg)

	implicit none

	integer lmin
	real hmin,dzreg

        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real hlv(1)
        common /hlv/hlv

	integer l
	real hbot

	if( dzreg .le. 0. ) stop 'error stop make_zeta_levels'

        hbot = hmin
	if( lmin .gt. 0 ) hlv(lmin) = hbot

        do l=lmin+1,nlvdi
          hbot = hbot + dzreg
          hlv(l) = hbot
        end do

	end

c********************************************************************

	subroutine set_hybrid_depth(lmax,zeta,htot,hlfem)

c sets depth structure and passes it back in hdep

	implicit none

	integer lmax
	real zeta
	real htot
	real hlfem(1)

	integer nsigma,l,i
	real hsigma,hsig

        real hlv(1)
        common /hlv/hlv

        call get_sigma(nsigma,hsigma)

        hsig = min(htot,hsigma) + zeta

	do l=1,nsigma-1
          hlfem(l) = -zeta - hsig * hlv(l)
	end do

	hlfem(nsigma) = -zeta + hsig

        do l=nsigma+1,lmax
          hlfem(l) = hlv(l)
        end do

        hlfem(lmax) = htot

c check ... may be deleted

	do l=2,lmax
	  if( hlfem(l) - hlfem(l-1) .le. 0. ) then
	    write(6,*) (hlfem(i),i=1,lmax)
	    stop 'error stop set_hybrid_depth: hlfem'
	  end if
	end do

	end

c********************************************************************






