
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c
c revision log :
c
c 01.01.2003    ggu     written
c
c description :
c
c node and element utility routines
c
c contents :
c
c subroutine delele(ie)		deletes element
c subroutine newele(nnew)	new element
c subroutine subnod(k,knew)     substitutes node in indices (k subst by knew)
c subroutine delnod(k)		deletes node
c subroutine newnod(nnew)	new node
c
c function ifindel(k1,k2,k3)
c               finds element given three nodes
c subroutine setele(ie,k1,k2,k3,nen3v)
c               inserts grade in index (in index of node k after node k1)
c subroutine subval(n,iarray,kold,knew)
c		substitutes in array iarray kold with knew
c
c***********************************************************

	subroutine delele(ie)

c deletes element

	use mod_adj_grade
	use mod_depth
	use basin

	implicit none

	integer ie

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

	subroutine newele(nnew)

c new element

	use mod_adj_grade
	use mod_depth
	use basin

	implicit none

	integer nnew

	integer ie,ii
	integer iemax

	iemax = 0
	do ie=1,nel
	  if( ipev(ie) .gt. iemax ) iemax = ipev(ie)
	end do

	nel = nel + 1
	if( nel .gt. neldi ) stop 'error stop newele: neldi'

	ipev(nel) = iemax + 1
	iarv(nel) = 0
	hev(nel) = 0.
	do ii=1,3
	  hm3v(ii,nel) = 0.
	end do

	nnew = nel

	end

c***********************************************************
c***********************************************************
c***********************************************************

	subroutine subnod(k,knew)

c substitutes node in indices (k subst by knew)

	use mod_adj_grade
	use mod_depth
	use basin

	implicit none

	integer k,knew

	integer ie,ii,n,kk,i

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

	subroutine delnod(k)

c deletes node

	use mod_adj_grade
	use mod_depth
	use basin

	implicit none

	integer k

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

	subroutine newnod(nnew)

c new node

	use mod_adj_grade
	use mod_depth
	use basin

	implicit none

	integer nnew
	integer k,kmax

	kmax = 0
	do k=1,nkn
	  if( ipv(k) .gt. kmax ) kmax = ipv(k)
	end do

	nkn = nkn + 1
	if( nkn .gt. nkndi ) stop 'error stop newnod: nkndi'

	ipv(nkn) = kmax + 1
	iarnv(nkn) = 0
	xgv(nkn) = 0.
	ygv(nkn) = 0.
	hkv(nkn) = 0.

	ngrade(nkn) = 0
	nbound(nkn) = 0

	nnew = nkn

	end

c***********************************************************
c***********************************************************
c***********************************************************

	function ifindel(k1,k2,k3)

c finds element given three nodes

	use basin

	implicit none

	integer ifindel
	integer k1,k2,k3
	integer ie,ii,iii,kn2,kn3

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

c inserts nodes in index

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

