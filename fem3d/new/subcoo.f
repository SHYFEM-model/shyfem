
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2019  Georg Umgiesser
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

c******************************************************************

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60

c******************************************************************

	subroutine coo_init(nel,nkn,mbw,nen3v,ndim,nnot0,ijp,ip,jp)

c construct pointers for coo matrix

	implicit none

	integer nel,nkn,mbw
	integer nen3v(3,1)
	integer ndim			!dimension of arrays ip, jp
	integer nnot0			!number of non 0 elements
	integer ijp(-mbw:mbw,nkn)	!index pointer into matrix
	integer ip(ndim)		!row index of non zero element
	integer jp(ndim)		!col index of non zero element

	integer k,m,n,ie,ii,ii1,k1,k2

	do k=1,nkn
	  do m=-mbw,mbw
	    ijp(m,k) = 0
	  end do
	end do

	n = 0

	do ie=1,nel
	  do ii=1,3
	    k1 = nen3v(ii,ie)
	    ii1 = mod(ii,3) + 1
	    k2 = nen3v(ii1,ie)
	    call coo_init_insert(k1,k2,nkn,mbw,ijp,n)	!out of diagonal
	    call coo_init_insert(k1,k1,nkn,mbw,ijp,n)	!in diagonal
	  end do
	end do

	if( n .gt. ndim) goto 99

	nnot0 = n

	do k=1,nkn
	  do m=-mbw,mbw
	    n = ijp(m,k)
	    jp(n) = k
	    ip(n) = k + m
	  end do
	end do

	call coo_check(nkn,mbw,ijp,ip,jp)

	return
   99	continue
	write(6,*) 'nnot0,ndim: ',nnot0,ndim
	stop 'error stop coo_init: ndim'
	end

c******************************************************************

	subroutine coo_find(i,j,mbw,ijp,n)

c finds position of non zero element in arrays

	implicit none

	integer i,j			!row and col
	integer mbw
	integer ijp(-mbw:mbw,1)
	integer n			!position

	n = ijp(i-j,j)

	end

c******************************************************************

	subroutine coo_check(nkn,mbw,ijp,ip,jp)

c checks sanity of ijp, ip, jp arrays

	implicit none

	integer nkn,mbw
	integer ijp(-mbw:mbw,1)
	integer ip(1)
	integer jp(1)

	integer i,j,k,m,n

	do k=1,nkn
	  do m=-mbw,mbw
	    n = ijp(m,k)
	    if( n .gt. 0 ) then
	      j = k
	      i = j + m
	      if( ip(n) .ne. i ) goto 99
	      if( jp(n) .ne. j ) goto 99
	    end if
	  end do
	end do

	return
   99	continue
	write(6,*) n,k,m
	write(6,*) i,ip(n)
	write(6,*) j,jp(n)
	stop 'error stop coo_check: internal error'
	end

c******************************************************************

	subroutine coo_init_insert(k1,k2,nkn,mbw,ip,n)

c internal routine for insertion of non 0 elements

	implicit none

	integer k1,k2,nkn,mbw,n
	integer ip(-mbw:mbw,nkn)

	integer idk,idk1,idk2
	integer ifact
	integer ip1,ip2

	idk = k1 - k2

	ifact = 0
	if( idk .lt. 0 ) ifact = -1
	if( idk .gt. 0 ) ifact = +1

	if( ifact .ne. 0 ) then
	  idk1 = -ifact * idk
	  idk2 = +ifact * idk
	  ip1 = ip(idk1,k2)
	  ip2 = ip(idk2,k1)
	  if( ip1 .ne. 0 .and. ip2 .ne. 0 ) then
	    return
	  else if( ip1 .eq. 0 .and. ip2 .eq. 0 ) then
	    ip(idk1,k2) = n + 1
	    ip(idk2,k1) = n + 2
	    n = n + 2
	  else
	    write(6,*) 'internal error: ',ip1,ip2
	    stop 'error stop coo_insert: internal error (1)'
	  end if
	else
	  ip1 = ip(0,k1)
	  if( ip1 .eq. 0 ) then
	    n = n + 1
	    ip(0,k1) = n
	  end if
	end if

	end
	  
c******************************************************************

