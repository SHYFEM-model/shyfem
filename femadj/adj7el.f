c
c $Id: adj7el.f,v 1.6 2007-03-20 13:19:42 georg Exp $
c
c description :
c
c 7 grade routines
c
c contents :
c
c subroutine elimh(nmax,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
c               eliminates high grades
c subroutine elim77(k,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
c               eliminates high grades
c
c revision log :
c
c 19.05.2003	ggu	in elim77 check angle of diagonals to exchange
c				...avoids negative area
c
c***********************************************************

	subroutine elimh(nmax,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)

c eliminates high grades

	implicit none

	integer nmax
        integer nkn,nel,ngrdim
        integer ngrade(1)
        integer nbound(1)
        integer ngri(2*ngrdim,1)
        integer nen3v(3,1)

        integer k,n

        write(6,*) 'eliminating grades for grades higher ... ',nmax

        do k=1,nkn
          n = ngrade(k)
          if( n .ge. nmax .and. nbound(k) .eq. 0 ) then
            call elim77(k,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
          end if
        end do

	end

c***********************************************************

	subroutine elim77(k,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)

c eliminates high grades

	implicit none

	integer k
        integer nkn,nel,ngrdim
        integer ngrade(1)
        integer nbound(1)
        integer ngri(2*ngrdim,1)
        integer nen3v(3,1)

	logical bdebug
        integer n,i,nc,nmax,nb,ii
	integer ie1,ie2
	integer np,nt,nn
	integer nval,ip
	integer nga(0:30)
	integer ngr(0:30)
	integer nba(0:30)

	real a1,a2

	integer ifindel
	real rangle

	if( k .gt. nkn ) return

	bdebug = .false.
c	if( k .eq. 543 ) bdebug = .true.
	if( k .eq. 2643 ) bdebug = .true.
	if( k .eq. 554 ) bdebug = .true.

c make circular list

        n = ngrade(k)
	nga(0) = ngri(n,k)
	do i=1,n
	  nga(i) = ngri(i,k)
	end do
	nga(n+1) = ngri(1,k)

	do i=0,n+1
	  ngr(i) = ngrade(nga(i))
	  nba(i) = 0
	  if( nbound(nga(i)) .ne. 0 ) then
	    ngr(i) = 6
	    nba(i) = 1
	  end if
	end do

	if( bdebug ) then
	  do i=1,n
	    write(6,*) nga(i-1),nga(i),nga(i+1)
	  end do
	  call plosno(k)
	end if

c check if exchange is possible

	nc = 0
	nmax = 0
	ip = 0
	do i=1,n
	  np = ngr(i-1)
	  nt = ngr(i)
	  nn = ngr(i+1)

	  nb = nba(i-1) + nba(i+1)

	  if( bdebug) write(6,*) '   ',n,nt,np,nn,nb

	  nval = n+nt - (np+nn)
	  if( nt .le. 5 ) nval = 0
	  if( np .ge. 7 .or. nn .ge. 7 ) nval = 0
	  a1 = rangle(nga(i+1),k,nga(i-1))
	  a2 = rangle(nga(i-1),nga(i),nga(i+1))
	  if( nval .gt. nmax ) then
	    if( a1 .le. 180. .or. a2 .le. 180. ) then
	      write(6,*) '************* not convex ',k,nga(i),a1,a2
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

	ie1 = ifindel(k,nga(ip),nga(ip+1))
	ie2 = ifindel(k,nga(ip-1),nga(ip))

	if( ie1 .eq. 0 .or. ie2 .eq. 0 ) then
	  stop 'error stop elim77: internal error (2)'
	end if

	if( bdebug ) then
	  write(6,*) ie1,k,nga(ip),nga(ip+1)
	  write(6,*) (nen3v(ii,ie1),ii=1,3)
	  write(6,*) ie2,k,nga(ip-1),nga(ip)
	  write(6,*) (nen3v(ii,ie2),ii=1,3)
	end if

	call setele(ie1,k,nga(ip-1),nga(ip+1),nen3v)
	call setele(ie2,nga(ip),nga(ip+1),nga(ip-1),nen3v)

	if( bdebug ) then
	  call prgr(k,ngrdim,ngrade,ngri)
	  call prgr(nga(ip),ngrdim,ngrade,ngri)
	  call prgr(nga(ip-1),ngrdim,ngrade,ngri)
	  call prgr(nga(ip+1),ngrdim,ngrade,ngri)
	end if

	call insgr(nga(ip-1),nga(ip),nga(ip+1),ngrdim,ngrade,ngri)
	call insgr(nga(ip+1),k,nga(ip-1),ngrdim,ngrade,ngri)
	call delgr(k,nga(ip),ngrdim,ngrade,ngri)
	call delgr(nga(ip),k,ngrdim,ngrade,ngri)

	if( bdebug ) then
	  call prgr(k,ngrdim,ngrade,ngri)
	  call prgr(nga(ip),ngrdim,ngrade,ngri)
	  call prgr(nga(ip-1),ngrdim,ngrade,ngri)
	  call prgr(nga(ip+1),ngrdim,ngrade,ngri)
	  call plosno(k)
	end if

	end if

	end

c************************************************

