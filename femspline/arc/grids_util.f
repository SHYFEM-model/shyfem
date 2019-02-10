
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c**********************************************************

	subroutine extrli(l,nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv,hkv
     +		,xt,yt,ht,nl,nt)

	implicit none

	integer l
	integer nli
	integer iplv(1)
	integer ialrv(1)
	integer ipntlv(0:1)
	integer inodlv(1)
	real xgv(1)
	real ygv(1)
	real hkv(1)
	integer nl
	real xt(1)
	real yt(1)
	real ht(1)
	integer nt

	integer nvert,i,ibase
	integer node,ier,k

	integer retriv

	if( l .gt. nli ) stop 'error stop extrli : l > nli'

	nvert = ipntlv(l) - ipntlv(l-1)
	ibase = ipntlv(l-1)

	if( nvert .gt. nl ) goto 98

	do i=1,nvert
	    node = inodlv(ibase+i)
	    ier = retriv(node,k)
	    if( ier .lt. 0 ) goto 99
	    xt(i) = xgv(k)
	    yt(i) = ygv(k)
	    ht(i) = hkv(k)
	end do

	nl = nvert
	nt = ialrv(l)

	write(6,*) 'extracting line ',l,iplv(l),nvert,nt

	return
   98	continue
	write(6,*) 'nvert = ',nvert
	write(6,*) 'nl = ',nl
	stop 'error stop extrli : dimension nl'
   99	continue
	stop 'error stop extrli: hash routines'
	end

c**********************************************************

	subroutine mkperiod(xt,yt,nl,bperiod)

c decides if line is periodic or not

	implicit none

	real xt(1)
	real yt(1)
	integer nl
	logical bperiod

	logical is_periodic

	bperiod = is_periodic(xt,yt,nl)
	if( bperiod ) nl = nl - 1

	end

c**********************************************************

	function is_periodic(xt,yt,nl)

c checks if line is periodic or not

	implicit none

	real xt(1)
	real yt(1)
	integer nl
	logical is_periodic

	if( xt(1) .eq. xt(nl) .and. yt(1) .eq. yt(nl) ) then
	  is_periodic = .true.
	else
	  is_periodic = .false.
	end if

	end

c**********************************************************

	subroutine prlisxy(l,nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv)

	implicit none

	integer l
	integer nli
	integer iplv(1)
	integer ialrv(1)
	integer ipntlv(0:1)
	integer inodlv(1)
	real xgv(1)
	real ygv(1)

	integer nvert,i,ibase
	integer node,ier,k

	integer retriv

	  nvert = ipntlv(l) - ipntlv(l-1)
	  ibase = ipntlv(l-1)
	  write(6,*) 'line : ',l,iplv(l),ialrv(l),nvert
	  do i=1,nvert
	    node = inodlv(ibase+i)
	    ier = retriv(node,k)
	    if( ier .lt. 0 ) goto 99
	    write(6,*) i,node,xgv(k),ygv(k)
	  end do

	return
   99	continue
	stop 'error stop prlixy: hash routines'
	end

c*********************************************************************

	subroutine prlixy(nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv)

	implicit none

	integer nli
	integer iplv(1)
	integer ialrv(1)
	integer ipntlv(0:1)
	integer inodlv(1)
	real xgv(1)
	real ygv(1)

	integer l

	do l=1,nli
	  call prlisxy(l,nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv)
	end do

	end

c**********************************************************

	subroutine prline(nli,iplv,ialrv,ipntlv,inodlv)

	implicit none

	integer nli
	integer iplv(1)
	integer ialrv(1)
	integer ipntlv(0:1)
	integer inodlv(1)

	integer l,nvert,i,ibase

	do l=1,nli
	  nvert = ipntlv(l) - ipntlv(l-1)
	  ibase = ipntlv(l-1)
	  write(6,*) 'line : ',l,iplv(l),ialrv(l),nvert
	  write(6,*) (inodlv(ibase+i),i=1,nvert)
	end do

	end

c**********************************************************

	subroutine insnod(nkn,ipv)

c inserts nodes into hash table

	implicit none

	integer nkn
	integer ipv(1)

	integer k,ier
	integer insert

	call init

	do k=1,nkn
	  ier = insert(ipv(k),k)
	  if( ier .lt. 0 ) goto 99
	end do

	return
   99	continue
	stop 'error stop insnod: hash routines'
	end

c********************************************************

	subroutine wrline(iunit,nline,nnode,nl,xt,yt,ht,nt,bperiod)

c write line

	integer iunit
	integer nline
	integer nnode
	integer nl
	real xt(1)
	real yt(1)
	real ht(1)
	integer nt
	logical bperiod

	integer i
	integer ntot

	if( nl .le. 0 ) return

	do i=1,nl
	  write(iunit,'(i1,2i10,3e14.6)') 1,nnode+i,nt,xt(i),yt(i),ht(i)
	end do
	  
c	nline = nline + 1

	ntot = nl
	if( bperiod ) ntot = ntot + 1

	write(iunit,'(i1,3i10)') 3,nline,nt,ntot
	write(iunit,*) (nnode+i,i=1,nl)
	if( bperiod ) then
	  write(iunit,*) nnode+1
	end if

	nnode = nnode + nl

	end

c********************************************************

