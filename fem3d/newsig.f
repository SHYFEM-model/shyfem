c
c routines for sigma levels
c
c revision log :
c
c 16.12.2010    ggu     program partially finished
c 19.09.2011    ggu     new routine set_bsigma()
c 04.11.2011    ggu     new routines for hybrid levels
c 10.11.2011    ggu     adjust depth for hybrid levels
c 11.11.2011    ggu     error check in set_hkv_and_hev()
c 11.11.2011    ggu     in check_hsigma_crossing set zeta levels to const depth
c 18.11.2011    ggu     restructured hybrid - adjustment to bashsigma
c 12.12.2011    ggu     eliminated (stupid) compiler bug (getpar)
c 27.01.2012    deb&ggu adapted for hybrid levels
c 23.02.2012    ccf	bug fix in set_hybrid_depth (no call to get_sigma)
c
c notes :
c
c important files where sigma levels are explicitly needed:
c
c	newini.f		set up of structure
c	subele.f		set new layer thickness
c
c	newbcl.f		for computation of rho
c	newexpl.f		for baroclinic term
c
c	lagrange_flux.f		limit zeta layers to surface layer
c
c********************************************************************

	subroutine set_sigma_hkhe

c this is called in newini.f

	implicit none

	logical berror,badjust,bzetaconti
	integer k,ie,ii
	integer nsigma,iadjust,idisc
	integer ihmin,ihmax
	integer idm,idp
	real hsigma,ha
	real h,hm
	real f(3)

        real hkv(1)
        common /hkv/hkv
        real hev(1)
        common /hev/hev

c-------------------------------------------------------
c initialize
c-------------------------------------------------------

	call get_sigma(nsigma,hsigma)

c-------------------------------------------------------
c check continuous and hsigma crossing
c-------------------------------------------------------

	call set_hkv_and_hev(hsigma,idm,idp)

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c********************************************************************

	subroutine set_hkv_and_hev(hsigma,idm,idp)

c sets hkv and hev

	implicit none

	real hsigma
	integer idm	!total number of non-continuous nodes < hsigma
	integer idp	!total number of non-continuous nodes > hsigma

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
	integer inc,ihmin,ihmax
	real hm,h

c-------------------------------------------------------
c initialize
c-------------------------------------------------------

	do k=1,nkn
	  hkv(k) = 0.
	  v1v(k) = 0.
	end do

	berror = .false.
	inc = 0

c-------------------------------------------------------
c set hkv and hev
c-------------------------------------------------------

	do ie=1,nel
	  hm = 0.
	  do ii=1,3
	    h = hm3v(ii,ie)
	    k = nen3v(ii,ie)
	    if( v1v(k) .eq. 0. ) then
	      hkv(k) = h
	      v1v(k) = 1.
	    else
	      if( h .ne. hkv(k) ) then
		write(6,*) 'depth of node not unique: ',ie,k,h,hkv(k)
	        inc = inc + 1
	      end if
	    end if
	    hm = hm + h
	  end do
	  hev(ie) = hm / 3.
	end do

	if( inc .gt. 0 ) then
	  write(6,*) 'number of occurences found: ',inc
	  stop 'error stop set_hkv_and_hev: depth not unique'
	end if

c-------------------------------------------------------
c check hkv
c-------------------------------------------------------

	do k=1,nkn
	  if( v1v(k) .le. 0. ) then
	    write(6,*) 'no depth in node: ',k
	    berror = .true.
	  end if
	end do

	if( berror ) stop 'error stop set_hkv_and_hev: no depth in node'

c-------------------------------------------------------
c check hsigma crossing
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
	    write(6,*) 'hsigma crossing: ',ie,(hm3v(ii,ie),ii=1,3)
	    berror = .true.
	  end if
	end do

	if( berror ) then
	  write(6,*) 'elements with hsigma crossing depths'
	  stop 'error stop set_hkv_and_hev: hsigma crossing'
	end if

c-------------------------------------------------------
c flatten elements of zeta coordinates
c-------------------------------------------------------

	do ie=1,nel
	  hm = hev(ie)
	  if( hm .gt. hsigma ) then
	    do ii=1,3
	      hm3v(ii,ie) = hm
	    end do
	  end if
	end do

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine get_bsigma(bsigma)

c returns bsigma which is true if sigma layers are used

	implicit none

	logical bsigma

	real getpar

	bsigma = nint(getpar('nsigma')) .gt. 0

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

	subroutine set_hybrid_depth(lmax,zeta,htot
     +					,hlv,nsigma,hsigma,hlfem)

c sets depth structure and passes it back in hlfem

	implicit none

	integer lmax
	real zeta
	real htot
	real hlv(1)
	integer nsigma
	real hsigma
	real hlfem(1)

	logical bsigma
	integer l,i
	real hsig

	bsigma = nsigma .gt. 0

	if( nsigma .gt. 0 ) then
          hsig = min(htot,hsigma) + zeta

	  do l=1,nsigma-1
            hlfem(l) = -zeta - hsig * hlv(l)
	  end do

	  hlfem(nsigma) = -zeta + hsig
	end if

        do l=nsigma+1,lmax
          hlfem(l) = hlv(l)
        end do

	if( nsigma .lt. lmax ) hlfem(lmax) = htot	!zeta or hybrid

c check ... may be deleted

	do l=2,lmax
	  if( hlfem(l) - hlfem(l-1) .le. 0. ) then
	    write(6,*) (hlfem(i),i=1,lmax)
	    stop 'error stop set_hybrid_depth: hlfem'
	  end if
	end do

	end

c********************************************************************






