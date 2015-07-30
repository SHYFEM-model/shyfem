c
c $Id: adj5el.f,v 1.5 2009-01-26 13:27:24 georg Exp $
c
c description :
c
c 5 grade routines
c
c contents :
c
c subroutine elim5(nkn,nel,ngrddi,ngrade,nbound,ngri,nen3v)
c			eliminates low grades
c subroutine elim55(k,nkn,nel,ngrddi,ngrade,nbound,ngri,nen3v)
c			eliminates 5-5 connections
c
c***********************************************************

	subroutine elim5

c eliminates low grades

	use mod_adj_grade
	use basin

	implicit none

        include 'param.h'

        integer k,n

        write(6,*) 'eliminating grades for grade 5... '

        do k=1,nkn
          n = ngrade(k)
          if( n .eq. 5 .and. nbound(k) .eq. 0 ) then
            call elim55(k)
c	    call node_info(1290)
	    call chkgrd(' ')
          end if
        end do

	end

c***********************************************************

	subroutine elim55(k)

c eliminates 5-5 connections

	use mod_adj_grade
	use basin

	implicit none

        include 'param.h'

	integer k

	logical bdebug
        integer n,i,nc,nmax,ii
	integer ie1,ie2
	integer ip1,ip2
	integer np,nt,nn
	integer nval,ip
	integer ngav(0:ngrdi+1)
	integer ngrv(0:ngrdi+1)
	integer nbav(0:ngrdi+1)

	integer ifindel

	if( k .gt. nkn ) return

	bdebug = .true.
	bdebug = .false.
	if( k .eq. 1138 ) bdebug = .true.

	if( bdebug ) then
	  write(6,*) '==============================================='
	  write(6,*) 'debug of new node: ',k
	end if

c make circular list
c
c ngav 	node numbers around k
c ngrv	grades of node numbers around k
c nbav	boundary  flag for nodes around k

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
	    ngrv(i) = 6	!FIXME
	    nbav(i) = 1
	  end if
	end do

c check if exchange is possible

	nc = 0		!how many of this nmax value
	nmax = 0	!maximum sum of grades -> must be at least 3
	ip = 0		!pointer to node in list that has been chosen
	do i=1,n
	  np = ngrv(i-1)
	  nt = ngrv(i)
	  nn = ngrv(i+1)

	  nval = np + nn - n - nt

	  if( nval .gt. nmax ) then
	    nc = 1
	    ip = i
	    nmax = nval
	  else if( nval .eq. nmax ) then
	    nc = nc + 1
	  end if
	end do

	if( nmax .lt. 3 ) return

	write(6,*) k,n,nmax,nc,ip

c nc gives number of occurences of this value of nmax ...
c ip is the pointer to the node to be exchanged
c
c we decide to take the first choice
c
c k is eliminated, ngav(ip) is retained (to account for boundary node)

	if( bdebug ) then
	    write(6,*) 'exchanging with node ... ',ngav(ip)
	    write(6,'(7i10)') (ngav(i),i=0,n+1)
	    write(6,'(7i10)') (ngrv(i),i=0,n+1)
	    write(6,'(7i10)') (nbav(i),i=0,n+1)
	    call plosno(k)
	    call plosno(ngav(ip))
	end if

c find elements that have to be deleted

	ie1 = ifindel(k,ngav(ip),ngav(ip+1))
	ie2 = ifindel(k,ngav(ip-1),ngav(ip))

	if( ie1 .eq. 0 .or. ie2 .eq. 0 ) then
	  stop 'error stop elim55: internal error (2)'
	end if

	if( bdebug ) then
	  write(6,*) 'elements to be deleted... ',ie1,ie2
	  write(6,*) ie1,k,ngav(ip),ngav(ip+1)
	  write(6,*) (nen3v(ii,ie1),ii=1,3)
	  write(6,*) ie2,k,ngav(ip-1),ngav(ip)
	  write(6,*) (nen3v(ii,ie2),ii=1,3)
	  call plosel2(ie1,ie2)
	end if

c delete elements

	if( ie1 .gt. ie2 ) then		!to avoid bug
	  call delele(ie1)
	  call delele(ie2)
	else
	  call delele(ie2)
	  call delele(ie1)
	end if

	if( bdebug ) then
	  write(6,*) 'grade index befor manipulation:'
	  call prgr(k,ngrdi,ngrade,ngri)
	  call prgr(ngav(ip),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip-1),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip+1),ngrdi,ngrade,ngri)
	end if

c new coordinates for node

	if( nbav(ip) .le. 0 ) then	!no boundary node
	  xgv(k) = 0.5 * ( xgv(k) + xgv(ngav(ip)) )
	  ygv(k) = 0.5 * ( ygv(k) + ygv(ngav(ip)) )
	end if

c substitute all occurrences of k with ngav(ip)

	call subnod(k,ngav(ip))

	if( bdebug ) then
	  write(6,*) 'after substitution...'
	  call prgr(ngav(ip),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip-1),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip+1),ngrdi,ngrade,ngri)
	end if

c adjourn grade (delete) for nodes ip-1, ip+1

	call delgr(ngav(ip-1),ngav(ip),ngrdi,ngrade,ngri)
	call delgr(ngav(ip+1),ngav(ip),ngrdi,ngrade,ngri)

	if( bdebug ) then
	  write(6,*) 'after deleting ip-1,ip+1...'
	  call prgr(ngav(ip),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip-1),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip+1),ngrdi,ngrade,ngri)
	end if

c adjourn grade index for ip and delete node k finally

	call delgr(ngav(ip),ngav(ip),ngrdi,ngrade,ngri)
	call delnod(k)
	call subval(n+2,ngav(0),nkn+1,k)	!if nkn is in ngav

	if( bdebug ) then
	  write(6,*) 'after deleting ip...'
	  call prgr(ngav(ip),ngrdi,ngrade,ngri)
	end if

	ip1 = mod(ip+2,n)
	ip2 = mod(ip+3,n)
	call insgr(ngav(ip),ngav(ip+1),ngav(ip1),ngrdi,ngrade,ngri)
	call insgr(ngav(ip),ngav(ip1),ngav(ip2),ngrdi,ngrade,ngri)

	if( bdebug ) then
	  write(6,*) 'grade index after manipulation:'
	  call prgr(ngav(ip),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip-1),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip+1),ngrdi,ngrade,ngri)
	end if

	if( bdebug ) then
	  call plosno(ngav(ip))
	  call plosno(ngav(ip-1))
	  call plosno(ngav(ip+1))
	end if

	if( bdebug ) then
	  write(6,*) 'end of debug of node ',k
	  write(6,*) '==============================================='
	end if

	end

