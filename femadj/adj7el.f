
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

c description :
c
c 7 grade routines
c
c contents :
c
c subroutine elimh(nmax,nkn,nel,ngrddi,ngrade,nbound,ngri,nen3v)
c               eliminates high grades
c subroutine elim77(k,nkn,nel,ngrddi,ngrade,nbound,ngri,nen3v)
c               eliminates high grades
c
c revision log :
c
c 01.01.2003	ggu	written
c 19.05.2003	ggu	in elim77 check angle of diagonals to exchange (area<0)
c 23.03.2010	ggu	changed v6.1.1
c 19.01.2015	ggu	changed VERS_7_1_2
c 24.07.2015	ggu	changed VERS_7_1_82
c 30.07.2015	ggu	changed VERS_7_1_83
c 12.10.2015	ggu	changed VERS_7_3_3
c 18.12.2018	ggu	changed VERS_7_5_52
c 13.03.2019	ggu	changed VERS_7_5_61
c 21.05.2019	ggu	changed VERS_7_5_62
c
c***********************************************************

	subroutine elimhigh(nmax)

c eliminates high grades

	use mod_adj_grade
	use basin

	implicit none

	integer nmax

        integer k,n

	if( nmax < 7 ) return

        write(6,*) 'eliminating grades for grades higher ... ',nmax

        do k=1,nkn
          n = ngrade(k)
          if( n .ge. nmax .and. nbound(k) .eq. 0 ) then
            call elim77(k)
          end if
        end do

	end

c***********************************************************

	subroutine elim77(k)

c eliminates high grades by switching diagonal

	use mod_adj_grade
	use basin

	implicit none

	integer k

	logical bdebug
        integer n,i,nc,nmax,nb,ii
	integer ie1,ie2
	integer np,nt,nn
	integer nval,ip
	integer ngav(0:ngrdi+1)
	integer ngrv(0:ngrdi+1)
	integer nbav(0:ngrdi+1)

	real a1,a2

	integer ifindel
	real rangle

	if( k .gt. nkn ) return

	bdebug = .false.
c	if( k .eq. 543 ) bdebug = .true.
c	if( k .eq. 2643 ) bdebug = .true.
c	if( k .eq. 554 ) bdebug = .true.

c make circular list

        n = ngrade(k)
	ngav(0) = ngri(n,k)
	do i=1,n
	  ngav(i) = ngri(i,k)
	end do
	ngav(n+1) = ngri(1,k)

	do i=0,n+1
	  ngrv(i) = ngrade(ngav(i))
	  nbav(i) = 0
	  if( nbound(ngav(i)) .ne. 0 ) then
	    ngrv(i) = 6
	    nbav(i) = 1
	  end if
	end do

	if( bdebug ) then
	  do i=1,n
	    write(6,*) ngav(i-1),ngav(i),ngav(i+1)
	  end do
	  call plosno(k)
	end if

c check if exchange is possible

	nc = 0
	nmax = 0
	ip = 0
	do i=1,n
	  np = ngrv(i-1)
	  nt = ngrv(i)
	  nn = ngrv(i+1)

	  nb = nbav(i-1) + nbav(i+1)

	  if( bdebug) write(6,*) '   ',n,nt,np,nn,nb

	  nval = n+nt - (np+nn)
	  if( nt .le. 5 ) nval = 0
	  if( np .ge. 7 .or. nn .ge. 7 ) nval = 0
	  a1 = rangle(ngav(i+1),k,ngav(i-1))
	  a2 = rangle(ngav(i-1),ngav(i),ngav(i+1))
	  if( nval .gt. nmax ) then
	    if( a1 .le. 180. .or. a2 .le. 180. ) then
	      write(6,*) '****** not convex ',k,ngav(i),a1,a2
	    else
	      nc = 1
	      ip = i
	      nmax = nval
	    end if
	  else if( nval .eq. nmax ) then
	    nc = nc + 1
	  end if
	end do

	if( nmax .ge. 3 ) write(6,*) k,n,nmax,nc,ip

c we decide to take the first choice

	if( nmax .ge. 3 ) then

	ie1 = ifindel(k,ngav(ip),ngav(ip+1))
	ie2 = ifindel(k,ngav(ip-1),ngav(ip))

	if( ie1 .eq. 0 .or. ie2 .eq. 0 ) then
	  stop 'error stop elim77: internal error (2)'
	end if

	if( bdebug ) then
	  write(6,*) ie1,k,ngav(ip),ngav(ip+1)
	  write(6,*) (nen3v(ii,ie1),ii=1,3)
	  write(6,*) ie2,k,ngav(ip-1),ngav(ip)
	  write(6,*) (nen3v(ii,ie2),ii=1,3)
	end if

	call setele(ie1,k,ngav(ip-1),ngav(ip+1),nen3v)
	call setele(ie2,ngav(ip),ngav(ip+1),ngav(ip-1),nen3v)

	if( bdebug ) then
	  call prgr(k,ngrdi,ngrade,ngri)
	  call prgr(ngav(ip),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip-1),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip+1),ngrdi,ngrade,ngri)
	end if

	call insgr(ngav(ip-1),ngav(ip),ngav(ip+1),ngrdi,ngrade,ngri)
	call insgr(ngav(ip+1),k,ngav(ip-1),ngrdi,ngrade,ngri)
	call delgr(k,ngav(ip),ngrdi,ngrade,ngri)
	call delgr(ngav(ip),k,ngrdi,ngrade,ngri)

	if( bdebug ) then
	  call prgr(k,ngrdi,ngrade,ngri)
	  call prgr(ngav(ip),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip-1),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip+1),ngrdi,ngrade,ngri)
	  call plosno(k)
	end if

	end if

	end

c************************************************

