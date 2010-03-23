c
c $Id: adj575.f,v 1.5 2007-03-20 13:19:42 georg Exp $
c
c description :
c
c 5-7-5 grade routines
c
c contents :
c
c subroutine elim57(nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
c			eliminates special 575 connections (shell)
c subroutine elim575(k,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
c			eliminates 5-7-5 connections
c
c***********************************************************

	subroutine elim57(nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)

c eliminates low grades

	implicit none

        integer nkn,nel,ngrdim
        integer ngrade(1)
        integer nbound(1)
        integer ngri(2*ngrdim,1)
        integer nen3v(3,1)

        integer k,n

        write(6,*) 'eliminating grades for grade 5-7-5... '

        do k=1,nkn
          n = ngrade(k)
          if( n .eq. 7 .and. nbound(k) .eq. 0 ) then
            call elim575(k,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
	    call chkgrd
          end if
        end do

	end

c***********************************************************

	subroutine elim575(k,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)

c eliminates 5-7-5 connections

	implicit none

	integer k
        integer nkn,nel,ngrdim
        integer ngrade(1)
        integer nbound(1)
        integer ngri(2*ngrdim,1)
        integer nen3v(3,1)

	logical bdebug
        integer n,i,nc,ii
	integer ie,kk,iii
	integer ip1,ip2
	integer ip
	integer ng,idp
	integer nga(30)
	integer ngr(30)
	integer nba(30)
	integer iau(30)
	real x,y,xm,ym

        real xgv(1), ygv(1)
        common /xgv/xgv, /ygv/ygv

	if( k .gt. nkn ) return

	bdebug = .true.
	bdebug = .false.

	if( bdebug ) write(6,*) 'new node: ',k

c make list

        n = ngrade(k)
	do i=1,n
	  nga(i) = ngri(i,k)
	end do

	do i=1,n
	  ngr(i) = ngrade(nga(i))
	  nba(i) = 0
	  if( nbound(nga(i)) .ne. 0 ) then
	    ngr(i) = 6
	    nba(i) = 1
	  end if
	end do

	if( bdebug ) then
	  do i=1,n
c	    write(6,*) nga(i-1),nga(i),nga(i+1)
	  end do
	end if

c check if exchange is possible

	nc = 0
	do i=1,n
	  ng = ngr(i)
	  if( ng .eq. 5 ) nc = nc + 1
	end do

	if( nc .ne. 2 ) return	!if not exactly 2 cannot proceed

c find out distance of 5 grades

	nc = 0
	do i=1,n
	  ng = ngr(i)
	  if( ng .eq. 5 ) then
	    nc = nc + 1
	    if( nc .eq. 1 ) then
		ip1 = i
	    else
		ip2 = i
	    end if
	  end if
	end do

	idp = ip2 - ip1
	if( idp .le. 2 .or. idp .ge. 5 ) return

	write(6,*) k,nc,ip1,ip2

	if( bdebug ) then
	  write(6,*) nga(ip1),nga(ip2)
	  write(6,'(7i10)') (nga(i),i=1,7)
	  write(6,'(7i10)') (ngr(i),i=1,7)
	  write(6,'(7i10)') (nba(i),i=1,7)
	end if

c	call plosno(k)

c reorder node list
c node 1 is a 5-grade, and node 5 is a 5-grade

	ip = ip1
	if( idp .eq. 3 ) ip = ip2
	call nshift(ip,n,nga,iau)
	call nshift(ip,n,ngr,iau)
	call nshift(ip,n,nba,iau)

	if( bdebug ) then
	  write(6,'(7i10)') (nga(i),i=1,7)
	  write(6,'(7i10)') (ngr(i),i=1,7)
	  write(6,'(7i10)') (nba(i),i=1,7)
	end if

c new node 

	call newnod(nkn,ngrade,nbound)

c substitute new node for old one in node index

	do ie=1,nel
	  do ii=1,3
	    kk = nen3v(ii,ie)
	    if( kk .eq. k ) then
	      do iii=1,3
		if( nen3v(iii,ie) .eq. nga(2) ) nen3v(ii,ie) = nkn
		if( nen3v(iii,ie) .eq. nga(4) ) nen3v(ii,ie) = nkn
	      end do
	    end if
	  end do
	end do

c new elements

	call newele(nel)
	call setele(nel,nga(1),nkn,k,nen3v)
	call newele(nel)
	call setele(nel,nga(5),k,nkn,nen3v)

c adjust grade index of old node (5 grade)

	call delgr(k,nga(2),ngrdim,ngrade,ngri)
	call delgr(k,nga(3),ngrdim,ngrade,ngri)
	call delgr(k,nga(4),ngrdim,ngrade,ngri)
	call insgr(k,nga(1),nkn,ngrdim,ngrade,ngri)

c adjust grade index of new node (6 grade)

	do i=1,5
	  ngri(i,nkn) = nga(i)
	end do
	ngri(6,nkn) = k
	ngrade(nkn) = 6

c adjust grade index of 5-5 nodes

	call insgrb(nga(1),k,nkn,ngrdim,ngrade,ngri)
	call insgr(nga(5),k,nkn,ngrdim,ngrade,ngri)

c substitute new node in grade index of nodes close to new node

	call exchgr(nga(2),k,nkn,ngrdim,ngrade,ngri)
	call exchgr(nga(3),k,nkn,ngrdim,ngrade,ngri)
	call exchgr(nga(4),k,nkn,ngrdim,ngrade,ngri)

	if( bdebug ) then
	  call prgr(k,ngrdim,ngrade,ngri)
	  call prgr(nkn,ngrdim,ngrade,ngri)
	  call prgr(nga(1),ngrdim,ngrade,ngri)
	  call prgr(nga(5),ngrdim,ngrade,ngri)
	end if

c adjust coordinates

	xm = 0.5 * ( xgv(nga(6)) + xgv(nga(7)) )
	ym = 0.5 * ( ygv(nga(6)) + ygv(nga(7)) )
	x = xgv(nga(3))
	y = ygv(nga(3))

	xgv(k) = xm + (1./3.) * ( x - xm )
	ygv(k) = ym + (1./3.) * ( y - ym )
	xgv(nkn) = xm + (2./3.) * ( x - xm )
	ygv(nkn) = ym + (2./3.) * ( y - ym )

c	call plosel2(nel-1,nel,nkn,nel,ngrdim,nen3v,ngrade,ngri)

	end

c*******************************************************
