c
c routines for sigma levels
c
c revision log :
c
c 16.12.2010    ggu     program partially finished
c 19.09.2011    ggu     new routine set_bsigma()
c 04.11.2011    ggu     new routines for hybrid levels
c 10.11.2011    ggu     adjust depth for hybrid levels
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

	badjust = .false.
	badjust = .true.	!adjust depth to hsigma value
	bzetaconti = .true.	!make also zeta values continuous

	call get_sigma(nsigma,hsigma)

	call set_hkv_and_hev(hsigma,idm,idp)

	if( idm .gt. 0 .or. idp .gt. 0 ) then
c	  this might be relaxed later
	  write(6,*) 'basin has non continuous depth values: ',idm,idp
	  stop 'error stop set_sigma_hkhe: depth values not on nodes'
	end if

c-------------------------------------------------------
c eliminate hsigma crossing
c-------------------------------------------------------

	iadjust = 0
	call check_hsigma_crossing(hsigma,iadjust)

	write(6,*) 'must adjust elements for hybrid levels: ',iadjust
	if( .not. badjust ) then
	  stop 'error stop set_sigma_hkhe: hsigma crossing'
	end if

	call check_hsigma_crossing(hsigma,iadjust)

	write(6,*) 'must adjust elements for hybrid levels: ',iadjust

c-------------------------------------------------------
c make continuous depth
c-------------------------------------------------------

	write(6,*) 'finalizing depths...'

	call set_hkv_and_hev(hsigma,idm,idp)
	write(6,*) 'non-continuous nodes: ',idm,idp

	call make_nodes_continuous(hsigma,bzetaconti)
	
	call set_hkv_and_hev(hsigma,idm,idp)
	write(6,*) 'non-continuous nodes: ',idm,idp
	!call check_hsigma_crossing(hsigma,iadjust)
	write(6,*) 'must adjust elements for hybrid levels: ',iadjust

c-------------------------------------------------------
c write out grid
c-------------------------------------------------------

	open(1,file='hsigma.grd',status='unknown',form='formatted')
	call wrgrd(1,hkv,hev,0)
	close(1)

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c********************************************************************

        subroutine wrgrd(iunit,hkv,hev,ike)

c writes grd file from bas
c
c hev or hkv must be set

        implicit none

        integer iunit
        integer ike             !ike==1 -> depth on elements
        real hkv(1)
        real hev(1)

        include 'param.h'
        include 'basin.h'

        integer k,ie,ii

        do k=1,nkn
          if( ike .eq. 1 ) then
            write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k)
          else
            write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k),hkv(k)
          end if
        end do

        write(iunit,*)

        do ie=1,nel
          if( ike .eq. 1 ) then
            write(iunit,1100) 2,ipev(ie),iarv(ie)
     +          ,3,(ipv(nen3v(ii,ie)),ii=1,3),hev(ie)
          else
            write(iunit,1100) 2,ipev(ie),iarv(ie)
     +          ,3,(ipv(nen3v(ii,ie)),ii=1,3)
          end if
        end do

        return
 1000   format(i1,2i10,3e16.8)
 1100   format(i1,2i10,i4,3i10,e16.8)
        end

c********************************************************************

	subroutine adjust_for_hybrid(hsigma,h,f)

c adjusts depth values for hybrid levels

	implicit none

	real hsigma
	real h(3),f(3)

	logical baver
	integer ii,iii,ip
	real hm
	real dh1,dh2

	baver = .false.

	hm = 0.
	do ii=1,3
	  hm = hm + h(ii)
	  f(ii) = 0.
	end do
	hm = hm / 3.

	if( baver ) then		!decide based on average

	if( hm .gt. hsigma ) then
	  do ii=1,3
	    if( h(ii) .lt. hsigma ) then
	      h(ii) = hsigma
	      f(ii) = 1.
	    end if
	  end do
	else
	  do ii=1,3
	    if( h(ii) .gt. hsigma ) then
	      h(ii) = hsigma
	      f(ii) = 1.
	    end if
	  end do
	end if

	else				!decide based on min/max

	do ii=1,3
	  iii = mod(ii,3) + 1
	  dh1 = h(ii)-hsigma
	  dh2 = h(iii)-hsigma
	  if( dh1*dh2 .lt. 0. ) then
	    if( abs(dh1) .lt. abs(dh2) ) then
	      ip = ii
	    else
	      ip = iii
	    end if
	    h(ip) = hsigma
	    f(ip) = 1.
	  end if
	end do

	end if

	end

c********************************************************************

	subroutine check_hsigma_crossing(hsigma,iadjust)

c checks and adjusts hsigma crossing

	implicit none

	real hsigma
	integer iadjust

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v
        real hm3v(3,1)
        common /hm3v/hm3v
        real v1v(1)
        common /v1v/v1v

	logical berror,bdebug
	integer k,ie,ii
	integer ihmin,ihmax
	real h
	real f(3)

	bdebug = .true.
	bdebug = .false.
	iadjust = 0

	do k=1,nkn
	  v1v(k) = 0.
	end do

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
	    iadjust = iadjust + 1
	    if(bdebug) write(6,*) 'before ',ie,(hm3v(ii,ie),ii=1,3)
	    call adjust_for_hybrid(hsigma,hm3v(1,ie),f)
	    if(bdebug) write(6,*) 'after  ',ie,(hm3v(ii,ie),ii=1,3)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      v1v(k) = v1v(k) + f(ii)
	    end do
	  end if
	end do

	do ie=1,nel
	  do ii=1,3
	    h = hm3v(ii,ie)
	    k = nen3v(ii,ie)
	    if( v1v(k) .ne. 0. ) then
	      hm3v(ii,ie) = hsigma
	    end if
	  end do
	end do

	end

c********************************************************************

	subroutine make_nodes_continuous(hsigma,bzetaconti)

c makes nodes continuous

	implicit none

	real hsigma
	logical bzetaconti

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
        real v2v(1)
        common /v2v/v2v

	logical berror,bdebug,bkdebug
	integer k,ie,ii
	real hm,h

c-------------------------------------------------------
c initialize
c-------------------------------------------------------

	berror = .false.

	do k=1,nkn
	  hkv(k) = 0.
	  v1v(k) = 0.
	  v2v(k) = 0.
	end do

c-------------------------------------------------------
c set hkv and hev
c-------------------------------------------------------

	do ie=1,nel
	  hm = 0.
	  bdebug = ie .eq. 4385
	  bdebug = ie .eq. -1
	  do ii=1,3
	    h = hm3v(ii,ie)
	    hm = hm + h
	    k = nen3v(ii,ie)
	    hkv(k) = hkv(k) + h
	    v1v(k) = v1v(k) + 1.
	    if( h .le. hsigma ) then
	      v2v(k) = v2v(k) - 1.
	    else
	      v2v(k) = v2v(k) + 1.
	    end if
	  end do
	if( bdebug ) write(6,*) ie,(hm3v(ii,ie),ii=1,3)
	  hev(ie) = hm / 3.
	end do

c-------------------------------------------------------
c check hkv
c-------------------------------------------------------

	do k=1,nkn
	  if( v1v(k) .le. 0. ) then
	    write(6,*) 'no depth in node: ',k
	    berror = .true.
	  else
	    hkv(k) = hkv(k)/v1v(k)
	  end if
	end do

	if( berror ) stop 'error stop set_hkv_and_hev: no depth in node'

c-------------------------------------------------------
c make continuous
c-------------------------------------------------------

	do ie=1,nel
	  hm = hev(ie)
	  bdebug = ie .eq. 4385 .or. ie .eq. 4384
	  bdebug = ie .eq. -1
	  if( bdebug ) write(6,*) ie,(hm3v(ii,ie),ii=1,3)
	  do ii=1,3
	    h = hm3v(ii,ie)
	    k = nen3v(ii,ie)
	    bkdebug = k .eq. 2542
	    bkdebug = k .eq. -1
	    if( bdebug ) write(6,*) ii,k,h,v1v(k),v2v(k)
	    if( v1v(k)-abs(v2v(k)) .gt. 0.5 ) then  !discontinuous over hsigma
	      if( hm .le. hsigma ) then
	        hm3v(ii,ie) = hsigma
	      else if( bzetaconti ) then	!make all continuous
	        hm3v(ii,ie) = hsigma
	      end if
	    else if( v2v(k) .lt. 0. ) then	!in sigma levels
	      hm3v(ii,ie) = hkv(k)
	    else if( bzetaconti ) then		!in zeta levels
	      hm3v(ii,ie) = hkv(k)
	    end if
	    if( bkdebug ) write(6,*) '== ',k,ie,hm3v(ii,ie)
	  end do
	  if( bdebug ) write(6,*) ie,(hm3v(ii,ie),ii=1,3)
	end do

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
	real hm,h

c-------------------------------------------------------
c initialize
c-------------------------------------------------------

	berror = .false.

	do k=1,nkn
	  hkv(k) = 0.
	  v1v(k) = 0.
	end do

c-------------------------------------------------------
c set hkv and hev
c-------------------------------------------------------

	do ie=1,nel
	  hm = 0.
	  do ii=1,3
	    h = hm3v(ii,ie)
	    k = nen3v(ii,ie)
	    hkv(k) = hkv(k) + h
	    v1v(k) = v1v(k) + 1.
	    hm = hm + h
	  end do
	  hev(ie) = hm / 3.
	end do

c-------------------------------------------------------
c check hkv
c-------------------------------------------------------

	do k=1,nkn
	  if( v1v(k) .le. 0. ) then
	    write(6,*) 'no depth in node: ',k
	    berror = .true.
	  else
	    hkv(k) = hkv(k)/v1v(k)
	  end if
	end do

	if( berror ) stop 'error stop set_hkv_and_hev: no depth in node'

c-------------------------------------------------------
c count non-continous nodes
c-------------------------------------------------------

	idm = 0
	idp = 0

	do ie=1,nel
	  hm = hev(ie)
	  do ii=1,3
	    h = hm3v(ii,ie)
	    k = nen3v(ii,ie)
	    if( abs(hkv(k)-h) .gt. 1.e-3 ) then
	      if( hm .le. hsigma ) then
	        idm = idm + 1
	      else
	        idp = idp + 1
	      end if
	    end if
	  end do
	end do

	write(6,*) 'non-continuous depth: ',idm,idp

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






