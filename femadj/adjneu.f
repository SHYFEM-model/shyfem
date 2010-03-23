c
c $Id: adjneu.f,v 1.5 2009-04-07 09:23:33 georg Exp $
c
c description :
c
c node and element utility routines
c
c contents :
c
c subroutine delele(ie,nkn,nel,ngrdim,ngrade,ngri)
c               deletes element
c subroutine subnod(k,knew,nkn,nel,ngrdim,ngrade,ngri)
c               substitutes node in indices (k subst by knew)
c subroutine delnod(k,nkn,nel,ngrdim,ngrade,ngri)
c               deletes node
c function ifindel(k1,k2,k3)
c               finds element given three nodes
c subroutine setele(ie,k1,k2,k3,nen3v)
c               inserts grade in index (in index of node k after node k1)
c
c***********************************************************

	subroutine delele(ie,nkn,nel,ngrdim,ngrade,ngri)

c deletes element

	implicit none

	integer ie,ngrdim
	integer nkn,nel
	integer ngrade(1)
	integer ngri(2*ngrdim,1)

	integer ipv(1)
	integer ipev(1)
	integer iarv(1)
	integer nen3v(3,1)
	real xgv(1), ygv(1)
	real hev(1), hkv(1)

c	common /nkon/ nkn,nel
	common /ipv/ipv
	common /ipev/ipev
	common /iarv/iarv
	common /nen3v/nen3v
	common /xgv/xgv, /ygv/ygv
	common /hev/hev, /hkv/hkv

	integer ieold,ienew,ii

	ienew = ie
	ieold = nel
	nel = nel - 1

	ipev(ienew) = ipev(ieold)
	iarv(ienew) = iarv(ieold)
	hev(ienew) = hev(ieold)

	do ii=1,3
	  nen3v(ii,ienew) = nen3v(ii,ieold)
	end do

	end

c***********************************************************

	subroutine newele(nel)

c new element

	implicit none

	integer nel

	integer ipv(1)
	integer ipev(1)
	integer iarv(1)
	integer nen3v(3,1)
	real xgv(1), ygv(1)
	real hev(1), hkv(1)

	common /ipv/ipv
	common /ipev/ipev
	common /iarv/iarv
	common /nen3v/nen3v
	common /xgv/xgv, /ygv/ygv
	common /hev/hev, /hkv/hkv

	integer ie
	integer iemax

	iemax = 0
	do ie=1,nel
	  if( ipev(ie) .gt. iemax ) iemax = ipev(ie)
	end do

	nel = nel + 1

	ipev(nel) = iemax
	iarv(nel) = 0
	hev(nel) = 0.

	end

c***********************************************************

	subroutine subnod(k,knew,nkn,nel,ngrdim,ngrade,ngri)

c substitutes node in indices (k subst by knew)

	implicit none

	integer k,ngrdim
	integer nkn,nel
	integer ngrade(1)
	integer ngri(2*ngrdim,1)

	integer ipv(1)
	integer ipev(1)
	integer iarv(1)
	integer nen3v(3,1)
	real xgv(1), ygv(1)
	real hev(1), hkv(1)

	common /ipv/ipv
	common /ipev/ipev
	common /iarv/iarv
	common /nen3v/nen3v
	common /xgv/xgv, /ygv/ygv
	common /hev/hev, /hkv/hkv

	integer knew,ie,ii,n,kk,i

	do ie=1,nel
	 do ii=1,3
	  if( nen3v(ii,ie) .eq. k ) nen3v(ii,ie) = knew
	 end do
	end do

	do kk=1,nkn
	 n = ngrade(kk)
	 do i=1,n
	  if( ngri(i,kk) .eq. k ) ngri(i,kk) = knew
	 end do
	end do

	end

c***********************************************************

	subroutine delnod(k,nkn,nel,ngrdim,ngrade,ngri)

c deletes node

	implicit none

	integer k,ngrdim
	integer nkn,nel
	integer ngrade(1)
	integer ngri(2*ngrdim,1)

	integer ipv(1)
	integer ipev(1)
	integer iarv(1)
	integer nen3v(3,1)
	real xgv(1), ygv(1)
	real hev(1), hkv(1)

	integer nbound(1)
	common /nbound/nbound

	common /ipv/ipv
	common /ipev/ipev
	common /iarv/iarv
	common /nen3v/nen3v
	common /xgv/xgv, /ygv/ygv
	common /hev/hev, /hkv/hkv

	logical bdebug
	integer knew,kold,ie,ii,n,kk,i
	integer kspecial

	kspecial = 1138
	bdebug = .false.
	if( k .eq. kspecial ) bdebug = .true.

	knew = k
	kold = nkn
	nkn = nkn - 1

	!call memnod(xgv(knew),ygv(knew))	!save deleted node for later

	ipv(knew) = ipv(kold)
	xgv(knew) = xgv(kold)
	ygv(knew) = ygv(kold)
	hkv(knew) = hkv(kold)

	if( bdebug ) then
	  write(6,*) 'delnod: ',knew,kold,nkn
	end if

	do ie=1,nel
	 do ii=1,3
	  if( nen3v(ii,ie) .eq. kold ) nen3v(ii,ie) = knew
	 end do
	end do

	ngrade(knew) = ngrade(kold)
	nbound(knew) = nbound(kold)
	n = ngrade(kold)
	do i = 1,n
	  ngri(i,knew) = ngri(i,kold)
	end do

	do kk=1,nkn
	 n = ngrade(kk)
	 do i=1,n
	  if( ngri(i,kk) .eq. kold ) ngri(i,kk) = knew
	 end do
	end do

	end

c***********************************************************

	subroutine newnod(nkn,ngrade,nbound)

c new node

	implicit none

	integer nkn
	integer ngrade(1)
	integer nbound(1)

	integer ipv(1)
	integer ipev(1)
	integer iarv(1)
	integer nen3v(3,1)
	real xgv(1), ygv(1)
	real hev(1), hkv(1)

	common /ipv/ipv
	common /ipev/ipev
	common /iarv/iarv
	common /nen3v/nen3v
	common /xgv/xgv, /ygv/ygv
	common /hev/hev, /hkv/hkv

	integer k,kmax

	kmax = 0
	do k=1,nkn
	  if( ipv(k) .gt. kmax ) kmax = ipv(k)
	end do

	nkn = nkn + 1

	ipv(nkn) = kmax
	xgv(nkn) = 0.
	ygv(nkn) = 0.
	hkv(nkn) = 0.

	ngrade(nkn) = 0
	nbound(nkn) = 0

	end

c***********************************************************

	function ifindel(k1,k2,k3)

c finds element given three nodes

	implicit none

	integer ifindel
	integer k1,k2,k3
	integer ie,ii,iii,kn2,kn3

	integer nkn,nel
	integer nen3v(3,1)
	common /nkon/ nkn,nel
	common /nen3v/nen3v

	ifindel = 0

	do ie=1,nel
	  do ii=1,3
	    if( nen3v(ii,ie) .eq. k1 ) then
		iii = mod(ii,3) +1
		kn2 = nen3v(iii,ie)
		iii = mod(iii,3) +1
		kn3 = nen3v(iii,ie)
		if( k2 .eq. kn2 .and. k3 .eq. kn3 ) then
		  ifindel = ie
		  return
		end if
	    end if
	  end do
	end do

	end

c**************************************************************

	subroutine setele(ie,k1,k2,k3,nen3v)

c inserts grade in index (in index of node k after node k1)

	implicit none

	integer ie,k1,k2,k3
	integer nen3v(3,1)

	nen3v(1,ie) = k1
	nen3v(2,ie) = k2
	nen3v(3,ie) = k3

	end

c**************************************************************

	subroutine subval(n,iarray,kold,knew)

c substitutes in array iarray kold with knew

	implicit none

	integer n
	integer iarray(1)
	integer kold,knew
	integer i

	do i=1,n
	  if( iarray(i) .eq. kold ) iarray(i) = knew
	end do

	end

c**************************************************************

