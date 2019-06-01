
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2010  Georg Umgiesser
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

c description :
c
c various utility routines for adj
c
c contents :
c
c subroutine smooth_grid(npass,omega,nkn,nel,nen3v,nbound,xgv,ygv,dx,dy,ic)
c               smoothing of internal nodes
c subroutine chkgrd
c               check of grid
c
c revision log :
c
c 01.01.2003	ggu	written
c 19.05.2003	ggu	some more utility routines
c 10.03.2010	ggu	area computation changed in checkarea (bug in 64 bit)
c 23.03.2010	ggu	changed v6.1.1
c 25.10.2013	ggu	changed VERS_6_1_68
c 19.01.2015	ggu	changed VERS_7_1_2
c 10.07.2015	ggu	changed VERS_7_1_50
c 24.07.2015	ggu	changed VERS_7_1_82
c 30.07.2015	ggu	changed VERS_7_1_83
c 12.10.2015	ggu	changed VERS_7_3_3
c 19.10.2015	ggu	changed VERS_7_3_6
c 18.12.2018	ggu	changed VERS_7_5_52
c 21.05.2019	ggu	changed VERS_7_5_62
c
c**************************************************************

	subroutine smooth_grid(npass,omega)

c smoothing of internal nodes

	use mod_adj_grade
	use basin

	implicit none

	integer npass
	real omega

	integer n,k,ie,ii
	integer ic(nkn)
	real dx(nkn)
	real dy(nkn)
	real xm,ym

	write(6,*) 'smoothing grid '

	call checkarea('before smoothing')

	do n=1,npass

	  !if( mod(n,10) .eq. 0 ) write(6,*) 'smoothing... pass ',n
	  if( mod(n,1) .eq. 0 ) write(6,*) 'smoothing... pass ',n

c initialize

	  do k=1,nkn
	    dx(k) = 0.
	    dy(k) = 0.
	    ic(k) = 0
	  end do

c loop over elements -> accumulate

	  do ie=1,nel
	    xm = 0.
	    ym = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      xm = xm + xgv(k)
	      ym = ym + ygv(k)
	    end do
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( nbound(k) .eq. 0 ) then
		ic(k) = ic(k) + 2
		dx(k) = dx(k) + 3. * xgv(k) - xm
		dy(k) = dy(k) + 3. * ygv(k) - ym
	      end if
	    end do
	  end do

c change coordinate

	  do k=1,nkn
	    if( ic(k) .gt. 0 ) then
		xgv(k) = xgv(k) - omega * dx(k) / ic(k)
		ygv(k) = ygv(k) - omega * dy(k) / ic(k)
	    end if
	  end do

c check area

	  call checkarea('during smoothing')

	end do

	write(6,*) 'smoothing finished - passes ',npass
	call checkarea('after smoothing')

	end

c**************************************************************

	subroutine mkstatic(nkn,ianv,nbound)

c marks nodes as static (not moveable)

	use mod_adj_static

	implicit none

	integer nkn
	integer ianv(1)
        integer nbound(1)

        integer k

        do k=1,nkn
          if( ianv(k) .eq. iastatic .and. nbound(k) .eq. 0 ) then
	    write(6,*) '*** marking node ',k,' (internal) as static'
            nbound(k) = nbstatic
          end if
        end do

        end

c**************************************************************

	subroutine chkgrd(text)

c check of grid
c
c checks consistency of node and grade index

	use mod_adj_grade
	use mod_adj_static
	use basin

	implicit none

	character*(*) text

	integer iaux(nkn)

	logical bstop,bverb
	integer n,k,ie,ii,kk,i,ke
	integer i1,i2,k1,k2
        integer nb

	integer nextgr

	bverb = .true.
	bverb = .false.

	bstop = .false.

c--------------------------------------------------
c check element index for unknown nodes
c--------------------------------------------------

	if( text .ne. ' ' ) write(6,*) 'chkgrd: '//trim(text)

        if( bverb ) write(6,*) 'checking for unknown nodes ...'

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( k .le. 0 .or. k .gt. nkn ) then
		bstop =.true.
		write(6,*) 'chkgrd (1): ',ie,ii,k
	    end if
	  end do
	end do

	if( bstop ) then
		call wr0grd
		write(6,*) 'nkn,nel: ',nkn,nel
		stop 'error stop chkgrd (1)'
	end if

c--------------------------------------------------
c check grades for strange grades
c--------------------------------------------------

        if( bverb ) write(6,*) 'checking for strange grades ...'

	do k=1,nkn
	  n = ngrade(k)
	  if( n .le. 0 ) then
		bstop =.true.
		write(6,*) 'chkgrd (2): ',k,n
	  end if
	end do

	if( bstop ) then
		call wr0grd
		write(6,*) 'nkn,nel: ',nkn,nel
		stop 'error stop chkgrd (2)'
	end if

c--------------------------------------------------
c check grade index for unknown nodes
c--------------------------------------------------

        if( bverb ) write(6,*) 'checking grade index for nodes ...'

	do k=1,nkn
	  n = ngrade(k)
	  do i=1,n
	    kk = ngri(i,k)
	    if( kk .le. 0 .or. kk .gt. nkn ) then
		bstop =.true.
		write(6,*) 'chkgrd (3): ',k,n,i,kk
	    end if
	  end do
	end do

	if( bstop ) then
		call wr0grd
		write(6,*) 'nkn,nel: ',nkn,nel
		stop 'error stop chkgrd (3)'
	end if

c--------------------------------------------------
c check consistency of grade index
c--------------------------------------------------

        if( bverb ) write(6,*) 'consistency check  ...'

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    i1 = mod(ii,3) + 1
	    k1 = nen3v(i1,ie)
	    i2 = mod(i1,3) + 1
	    k2 = nen3v(i2,ie)

	    kk = nextgr(k,k1,ngrdi,ngrade,ngri)
	    if( kk .ne. k2 ) then
		bstop =.true.
		write(6,*) 'chkgrd (4): ',ie,k,k1,k2,kk
	    end if
	  end do
	end do

	if( bstop ) then
		call wr0grd
		write(6,*) 'nkn,nel: ',nkn,nel
		stop 'error stop chkgrd (4)'
	end if

c--------------------------------------------------
c is grade still ok?
c--------------------------------------------------

        if( bverb ) write(6,*) 'checking for final grades ...'

	do k=1,nkn
	  iaux(k) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    iaux(k) = iaux(k) + 1
	  end do
	end do

	do k=1,nkn
          nb = nbound(k)
	  if( nb .ne. 0 .and. nb .ne. nbstatic ) then	!boundary node
	    if( iaux(k) + 1 .ne. ngrade(k) ) then
		write(6,*) 'chkgrd (5a): ',k,nbound(k),iaux(k),ngrade(k)
		ke = ipv(k)
		write(6,*) '   (problem in boundary node ',ke,' )'
		bstop =.true.
	    end if
	  else				!internal node
	    if( iaux(k) .ne. ngrade(k) ) then
		write(6,*) 'chkgrd (5b): ',k,nbound(k),iaux(k),ngrade(k)
		ke = ipv(k)
		write(6,*) '   (problem in internal node ',ke,' )'
		bstop =.true.
	    end if
	  end if
	end do

	if( bstop ) then
		call wr0grd
		write(6,*) 'nkn,nel: ',nkn,nel
		stop 'error stop chkgrd (5)'
	end if

c--------------------------------------------------
c checking area
c--------------------------------------------------

	call checkarea(text)

c--------------------------------------------------
c end of routine
c--------------------------------------------------

        if( bverb ) write(6,*) 'check ok ...'

	end

c*******************************************************

	subroutine nshift(ip,n,nga,iau)

c shifts array - index ip will be first index
c iau is auxiliary array

	implicit none

	integer ip,n
	integer nga(n),iau(n)

	integer i,ipa

	ipa = ip - 1 
	do i=1,n
	  ipa = mod(ipa,n) + 1
	  iau(i) = nga(ipa)
	end do

	do i=1,n
	  nga(i) = iau(i)
	end do

	end

c*******************************************************

	subroutine checkarea(text)

c check if area is positive

	use basin

        implicit none

	character*(*) text

	logical bstop
	integer ie,ii,i1,i2,k1,k2
        integer ieext
	real x1,y1,x2,y2,x3,y3
	real aj				!is twice the area
	character*80 string

	real areat

	bstop = .false.
	string = text

	do ie=1,nel
	  x1 = xgv(nen3v(1,ie))
	  x2 = xgv(nen3v(2,ie))
	  x3 = xgv(nen3v(3,ie))
	  y1 = ygv(nen3v(1,ie))
	  y2 = ygv(nen3v(2,ie))
	  y3 = ygv(nen3v(3,ie))
	  aj = areat(x1,y1,x2,y2,x3,y3)
	  if( aj .le. 0. ) then
            call eint2ext(ie,ieext)
	    write(6,*) 'element ',ie,' (extern ',ieext
     +                          ,') has negative area...'
	    call elem_info(ie)
	    bstop = .true.
	  end if
        end do

	if( bstop ) then
	    write(6,*) 'error while checking: ',trim(string)
	    call wr0grd
	    write(6,*) 'nkn,nel: ',nkn,nel
	    stop 'error stop checkarea'
	end if

	end

c*******************************************************

	subroutine elem_info(ie)

	use mod_adj_grade
	use basin

        implicit none

	integer ie,ieext,ii,k

        call eint2ext(ie,ieext)

	write(6,*) 'info on element ie = ',ie,' extern = ',ieext

	do ii=1,3
	  k = nen3v(ii,ie)
	  call node_info(k)
	end do

	end

c*******************************************************

	subroutine node_info(k)

c writes info on node

	use mod_adj_grade
	use basin

        implicit none

	integer k,kext

	integer i

	if( k .le. 0 ) return

        call nint2ext(k,kext)

	write(6,*) 'info on node k = ',k,' extern = ',kext
	write(6,*) ngrade(k),nbound(k),xgv(k),ygv(k)
	write(6,*) nkn,nel,ngrdi
	write(6,*) (ngri(i,k),i=1,ngrade(k))

	end

c****************************************************************

	function rangle(k1,k2,k3)

c gives angle between nodes
c
c rangle	angle [degrees] ,     rangle < 180 => right turn
c k1,k2,k3	node numbers

	use basin

	implicit none

	real rangle
	integer k1,k2,k3

	real x1,y1,x2,y2,x3,y3
	real angle

	x1 = xgv(k1)
	y1 = ygv(k1)
	x2 = xgv(k2)
	y2 = ygv(k2)
	x3 = xgv(k3)
	y3 = ygv(k3)

	rangle = angle(x1,y1,x2,y2,x3,y3)

	end

c************************************************************

        subroutine nint2ext(kint,kext)

	use mod_adj_grade
	use basin

        implicit none

        integer kint,kext

        kext = ipv(kint)

        end

c************************************************************

        subroutine eint2ext(ieint,ieext)

	use mod_adj_grade
	use basin

        implicit none

        integer ieint,ieext

        ieext = ipev(ieint)

        end

c************************************************************

