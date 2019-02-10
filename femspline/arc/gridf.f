
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

	program gridf

	implicit none

	integer nkndim,neldim,nlidim,nlndim,ndim
	parameter( nkndim = 50000 )
	parameter( neldim = 100000 )
	parameter( nlidim = 10000 )
	parameter( nlndim = 100000 )
	parameter( ndim = 100000 )

	integer ner
	integer nco,nkn,nknh,nel,nelh,nli
	logical bstop
	character*80 file

	integer ipv(nkndim)
	integer ipev(neldim)
	integer iarv(neldim)
	integer nen3v(3,neldim)

	real xgv(nkndim)
	real ygv(nkndim)
	real hkv(nkndim)
	real hev(neldim)

	integer iplv(nlidim)
	integer ialrv(nlidim)
	integer ipntlv(0:nlidim)
	integer inodlv(nlndim)

	integer nl,nt
	real xt(ndim)
	real yt(ndim)

	integer l
	integer nline
	integer nnode
	real sigma
	real reduct

	write(6,*) 'program reads from stdin'
	write(6,*) 'Please use as: gridf < file.grd'

c-----------------------------------
c sigma		smoothing parameter (size of gaussian kernel)
c reduct	reduction of points in line (2-> half points)
c-----------------------------------
	sigma = 2.0
	sigma = 1.5
	reduct = 4.0
c-----------------------------------

	nline = 0
	nnode = 0

	ner = 6
	bstop = .false.
	file = ' '

        call rdgrd(file,ner,bstop,nco,nkn,nknh,nel,nelh,nli
     +                  ,nkndim,neldim,nlidim
     +                  ,ipv,ipev,iarv,nen3v,xgv,ygv,hev,hkv
     +                  ,nlndim,iplv,ialrv,ipntlv,inodlv)

	write(6,*) 'comments : ',nco
	write(6,*) 'nodes    : ',nkn,nknh
	write(6,*) 'elements : ',nel,nelh
	write(6,*) 'lines    : ',nli

	call insnod(nkn,ipv)
c	call prline(nli,iplv,ialrv,ipntlv,inodlv)
c	call prlixy(nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv)

	open(99,file='smooth.grd',status='unknown',form='formatted')
	open(98,file='reduce.grd',status='unknown',form='formatted')

	do l=1,nli
	  nl = ndim
	  call extrli(l,nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv,xt,yt,nl,nt)
	  call smooth(sigma,xt,yt,nl)
	  call wrline(99,nline,nnode,nl,xt,yt,nt)
	  call reduce(reduct,xt,yt,nl)
	  call wrline(98,nline,nnode,nl,xt,yt,nt)
	end do

	close(99)
	close(98)

	write(6,*) 'files smooth.grd and reduce.grd written'

	end

c**********************************************************

	subroutine extrli(l,nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv
     +		,xt,yt,nl,nt)

	implicit none

	integer l
	integer nli
	integer iplv(1)
	integer ialrv(1)
	integer ipntlv(0:1)
	integer inodlv(1)
	real xgv(1)
	real ygv(1)
	integer nl
	real xt(1)
	real yt(1)
	integer nt

	integer nvert,i,ibase
	integer node,ier,k

	integer retriv

	if( l .gt. nli ) stop 'error stop extrli : l > nli'

	nvert = ipntlv(l) - ipntlv(l-1)
	ibase = ipntlv(l-1)

	if( nvert .gt. nl ) stop 'error stop extrli : dimension nl'

	write(6,*) 'extracting line ',l,nvert,nt

	do i=1,nvert
	    node = inodlv(ibase+i)
	    ier = retriv(node,k)
	    xt(i) = xgv(k)
	    yt(i) = ygv(k)
	end do

	nl = nvert
	nt = ialrv(l)

	end

c**********************************************************

	subroutine prlixy(nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv)

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
	    write(6,*) i,node,xgv(k),ygv(k)
	  end do
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
	end do

	end

c********************************************************

	subroutine smooth(sigma,xt,yt,nl)

c smoothing

	real sigma
	real xt(1)
	real yt(1)
	integer nl

	integer ndim
	parameter(ndim=10000)

	integer ngk
	real gk(-ndim:ndim)
	real raux(-ndim:2*ndim)

c handle closed curve

	if( xt(1) .eq. xt(nl) .and. yt(1) .eq. yt(nl) ) nl = nl - 1

	if( nl .gt. ndim ) stop 'error stop smooth: ndim'

c create smoothing kernel

	ngk = ndim
	call gkernel(sigma,ndim,ngk,gk)
	write(6,*) 'Gaussian kernel : ',sigma,ngk

c smooth x and y

	call gsmooth(ndim,nl,xt,raux,ngk,gk)
	call gsmooth(ndim,nl,yt,raux,ngk,gk)

	end

c********************************************************

	subroutine gkernel(sigma,ndim,ngk,gk)

c gaussian smoothing

	real sigma
	integer ndim
	integer ngk
	real gk(-ndim:ndim)

	integer i
	real pi,eps,a,b,x,val

        pi = 4.*atan(1.)
        eps = 1.e-7

        a = 1. / ( sqrt(2.*pi) * sigma )
        b = 1. / ( 2 * sigma * sigma )

	do i=0,ndim
          x = i
          val = a * exp( -x**2 * b )
	  gk(i) = val
	  gk(-i) = val
	  if( val .lt. eps ) goto 1
        end do
    1	continue
	ngk = i-1

c	do i=-ngk,ngk
c	  write(6,*) i,gk(i)
c	end do

	end

c********************************************************

	subroutine gsmooth(ndim,nl,rt,raux,ngk,gk)

c gaussian smoothing

	integer ndim
	integer nl
	real rt(1)
	real raux(-ndim:2*ndim)
	integer ngk
	real gk(-ndim:ndim)

	integer i,j
	real val

c set up circular array

	do i=1,nl
	  val = rt(i)
	  raux(i) = val
	  raux(i+nl) = val
	  raux(i-nl) = val
	end do

c convolution

	do i=1,nl
	  val = 0.
	  do j=-ngk,ngk
	    val = val + raux(i+j) * gk(j)
	  end do
	  rt(i) = val
	end do

	end

c********************************************************

	subroutine wrline(iunit,nline,nnode,nl,xt,yt,nt)

c write line

	integer iunit
	integer nline
	integer nnode
	integer nl
	real xt(1)
	real yt(1)
	integer nt

	integer i

	do i=1,nl
	  write(iunit,'(i1,2i10,2e14.6)') 1,nnode+i,nt,xt(i),yt(i)
	end do
	  
	nline = nline + 1

	write(iunit,'(i1,3i10)') 3,nline,nt,nl+1
	write(iunit,*) (nnode+i,i=1,nl),nnode+1

	nnode = nnode + nl

	end

c********************************************************

	  subroutine reduce(reduct,xt,yt,nl)

c reduces points in line

	implicit none

	real reduct
	real xt(1), yt(1)
	integer nl

	integer i,nnew
	real rr,rtot

c very simplicistic approach

	rr = reduct
	rtot = 1.

	if( nl .lt. 3 * rr ) then
	  rr = nl / 3.
	end if

	nnew = 0
	do i=1,nl
	  if( i .ge. rtot ) then
	    nnew = nnew + 1
	    xt(nnew) = xt(i)
	    yt(nnew) = yt(i)
	    rtot = rtot + rr
	  end if
	end do

	if( nnew .lt. 3 ) stop 'error stop reduce: line too short...'

	nl = nnew

	end

