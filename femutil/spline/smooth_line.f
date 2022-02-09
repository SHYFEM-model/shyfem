
!--------------------------------------------------------------------------
!
!    Copyright (C) 2002,2004-2005,2011-2012  Georg Umgiesser
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

!**********************************************************

	subroutine extrli(l,nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv,hkv
     +		,xt,yt,ht,kt,nl,nt)

! extracts line from list

	implicit none

	integer l		!actual line
	integer nli		!total number of lines
	integer iplv(1)		!external line number
	integer ialrv(1)	!line type
	integer ipntlv(0:1)
	integer inodlv(1)
	real xgv(1)
	real ygv(1)
	real hkv(1)
	integer nl		!on entry dim, on return number of nodes in line
	real xt(1)
	real yt(1)
	real ht(1)
	integer kt(1)		!external node number
	integer nt		!type of line

	integer nvert,i,ibase
	integer node,ier,k,nint,next

	integer retriv
	integer ipext

	if( l .gt. nli ) stop 'error stop extrli : l > nli'

	nvert = ipntlv(l) - ipntlv(l-1)
	ibase = ipntlv(l-1)

	if( nvert .gt. nl ) goto 98

	write(6,*) 'extracting line ',l,iplv(l),nvert,nt

	do i=1,nvert
	    nint = inodlv(ibase+i)
	    next = ipext(nint)
	    ier = retriv(next,k)
	    if( ier .le. 0 ) goto 99
	    xt(i) = xgv(k)
	    yt(i) = ygv(k)
	    ht(i) = hkv(nint)
	    kt(i) = next
	end do

	nl = nvert
	nt = ialrv(l)

	return
   98	continue
	write(6,*) '*** Dimension error for nl'
	write(6,*) 'Number of vertices in line: nvert = ',nvert
	write(6,*) 'Dimension for vertices:        nl = ',nl
	write(6,*) 'Please adjust dimension of ndim'
	stop 'error stop extrli : dimension nl'
   99	continue
	write(6,*) 'error retrieving node: ',node,ier
	stop 'error stop extrli: hash routines'
	end

!**********************************************************

	subroutine mkperiod(xt,yt,nl,bperiod)

! decides if line is periodic or not

	implicit none

	real xt(1)
	real yt(1)
	integer nl
	logical bperiod

	bperiod = .false.

	if( xt(1) .eq. xt(nl) .and. yt(1) .eq. yt(nl) ) then
	  bperiod = .true.
	  nl = nl - 1
	end if

	end

!**********************************************************

	subroutine prlixy(nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv)

! print info on lines

	implicit none

	integer nli
	integer iplv(1)
	integer ialrv(1)
	integer ipntlv(0:1)
	integer inodlv(1)
	real xgv(1)
	real ygv(1)

	integer l,nvert,i,ibase
	integer node,ier,k

	integer retriv

	do l=1,nli
	  nvert = ipntlv(l) - ipntlv(l-1)
	  ibase = ipntlv(l-1)
	  write(6,*) 'line : ',l,iplv(l),ialrv(l),nvert
	  do i=1,nvert
	    node = inodlv(ibase+i)
	    ier = retriv(node,k)
	    if( ier .lt. 0 ) goto 99
	    write(6,*) i,node,xgv(k),ygv(k)
	  end do
	end do

	return
   99	continue
	stop 'error stop prlixy: hash routines'
	end

!**********************************************************

	subroutine prline(nli,iplv,ialrv,ipntlv,inodlv)

! print info on lines

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

!**********************************************************

	subroutine insnod(nkn,ipv)

! inserts nodes into hash table

	implicit none

	integer nkn
	integer ipv(1)

	integer k,ier,kk
	integer insert,retriv

	call init

	do k=1,nkn
	  ier = insert(ipv(k),k)
	  if( ier .le. 0 ) goto 99
	  ier = retriv(ipv(k),kk)
	  if( ier .le. 0 ) goto 99
	  if( k /= kk ) then
	    write(6,*) 'k/=kk ',k,kk
	    goto 99
	  end if
	end do

	return
   99	continue
	stop 'error stop insnod: hash routines'
	end

!********************************************************

	subroutine wrline(iunit,nline,nnode,nl,xt,yt,ht,nt,bperiod)

! write line

	integer iunit
	integer nline
	integer nnode
	integer nl
	real xt(nl)
	real yt(nl)
	real ht(nl)
	integer nt
	logical bperiod

	integer i
	integer ntot

	if( nl .le. 0 ) return

	do i=1,nl
	  write(iunit,'(i1,2i10,3e14.6)') 1,nnode+i,nt,xt(i),yt(i),ht(i)
	end do
	  
!	nline = nline + 1

	ntot = nl
	if( bperiod ) ntot = ntot + 1

	write(iunit,'(i1,3i10)') 3,nline,nt,ntot
	write(iunit,'(7i10)') (nnode+i,i=1,nl)
	if( bperiod ) write(iunit,'(i10)') nnode+1

	nnode = nnode + nl

	end

!********************************************************

