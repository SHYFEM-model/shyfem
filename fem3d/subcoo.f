
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009-2010,2012,2014-2015,2018-2019  Georg Umgiesser
!    Copyright (C) 2015  Debora Bellafiore
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

c generic coo and csr handling routines
c
c revision log :
c
c 31.01.2009	ggu	cleaned, only generic coo routines here
c 23.03.2010	ggu	changed v6.1.1
c 29.03.2012	ggu	in loccoo avoid call to coo_find (speed)
c 05.11.2012	ggu	changed VERS_6_1_60
c 12.12.2014	ggu	changed VERS_7_0_9
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 18.09.2015	ggu	bug fix in coo_init - do not set ip/jp for n==0
c 04.12.2015	ggu	new approach for constructing matrix - also 3d
c 15.12.2015	ggu&dbf	finished and validated 3d approach
c 23.04.2018	ggu	new matrix type
c 21.12.2018	ggu	changed VERS_7_5_53
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c
!******************************************************************

	subroutine coo_init_new(matrix)

! construct pointers for coo matrix format
!
! 3d system will have size: nkn * (nlv+2)

	use mod_system

	implicit none

        type(smatrix), target :: matrix

	logical bnohydro
	logical bcheck
	integer k,n,i,ie,ii,m,l,iii
	integer kn(3),ki,kj,j
	integer kic,kjc,kijc,iiis,iiie
	integer ipp,ipp0
	integer nn
	integer, allocatable :: nodes(:)
	integer, allocatable :: nlist(:,:)
	integer nlv,nkn,nel,ngr
	integer n2zero,n3zero,n3max

	integer loccoo3d

        type(smatrix), pointer :: mm

        mm => matrix

	bcheck = .not. mm%bglobal	!check only if not global matrix

!------------------------------------------------------------------
! contruct grades
!------------------------------------------------------------------

	nkn = mm%nkn_system
	nel = mm%nel_system
	nlv = mm%nlv_system
	ngr = mm%ngr_system

	mm%ntg(0) = 0
	mm%nt3g(0) = 0

	allocate(nodes(ngr+1))
	allocate(nlist(0:2*ngr+1,nkn))
	call construct_node_list(nkn,nel,ngr,mm%nen3v_system,nlist)

	do k=1,nkn
	  n = nlist(0,k)
	  nodes(1:n) = nlist(1:n,k)
	  mm%ng(k) = n
	  mm%ntg(k) = mm%ntg(k-1) + n
	  nn = (n-1)*nlv + 2 + 3*nlv		!3d entries in row k
	  mm%n3g(k) = nn
	  mm%nt3g(k) = mm%nt3g(k-1) + nn
	  mm%iorder(1:n,k) = nodes(1:n)
	  do i=1,n
	    if( nodes(i) == k ) exit
	  end do
	  if( i > n ) goto 99
	  mm%diag(k) = i
	end do

!------------------------------------------------------------------
! check grades - can be deleted later
!------------------------------------------------------------------

	if( bcheck ) then			!can be deleted later
	 write(6,*) 'checking node list...'
	 do k=1,nkn
	  call get_nodes_around(k,ngr,n,nodes)
	  n = n + 1
	  nodes(n) = k				!add diagonal node
	  call sort_int(n,nodes)
	  if( n /= nlist(0,k) 
     +			.or. any( nodes(1:n) /= nlist(1:n,k) ) ) then
	    write(6,*) nkn,nel,ngr
	    write(6,*) k,n,nlist(0,k)
	    write(6,*) nodes(1:n)
	    write(6,*) nlist(1:n,k)
	    stop 'error stop coo_init_new: internal error (2)'
	  end if
	 end do
	end if

!------------------------------------------------------------------
! contruct 2d pointer
!------------------------------------------------------------------

	n2zero = mm%ntg(nkn)	!set global 2d value
	mm%n2zero = n2zero
	mm%ijp_ie = 0
	mm%i2coo = 0
	mm%j2coo = 0

	do ie=1,nel
	  do ii=1,3
	    kn(ii) = mm%nen3v_system(ii,ie)
	  end do
	  do i=1,3
	    ki = kn(i)				!row
	    ipp0 = mm%ntg(ki-1)			!last entry in row ki-1
	    n = mm%ng(ki)				!entries in row ki
	    nodes(1:n) = mm%iorder(1:n,ki)		!nodes of row ki
	    do j=1,3
	      kj = kn(j)			!col
	      do m=1,n
	        if( nodes(m) == kj ) exit	!find kj in nodes
	      end do
	      if( m > n ) goto 98
	      ipp = ipp0 + m
	      if( ipp > n2zero ) goto 97
	      mm%ijp_ie(i,j,ie) = ipp
	      mm%i2coo(ipp) = ki			!row
	      mm%j2coo(ipp) = kj			!col
	    end do
	  end do
	end do

	if( .not. bsys3d ) return

!------------------------------------------------------------------
! contruct 3d pointer
!------------------------------------------------------------------

	n3max = mm%n3max
	n3zero = mm%ntg(nkn)*nlv + nkn*(2+3*nlv)	!set global 3d value
	if( n3zero /= mm%nt3g(nkn) ) goto 91
	if( n3zero > n3max ) goto 91
	mm%n3zero = n3zero

	mm%i3coo = -99
	mm%j3coo = -99
	mm%back3coo = -88

	do ie=1,nel
	  do ii=1,3
	    kn(ii) = mm%nen3v_system(ii,ie)
	  end do
	  do i=1,3
	    ki = kn(i)				!row
	    do j=1,3
	      kj = kn(j)			!col
	      do l=1,nlv
		ipp = loccoo3d(i,j,kn,l,ie)
	        if( ipp > n3max ) goto 96
	        if( ipp < 0 ) goto 96
	        kic = (ki-1)*(nlv+2) + l + 1
	        kjc = (kj-1)*(nlv+2) + l + 1
		iiis = 0
		iiie = 0
		if( i == j ) then
		  iiis = -1
		  iiie = +1
		end if
		do iii=iiis,iiie
	          mm%i3coo(ipp+iii) = kic
	          mm%j3coo(ipp+iii) = kjc + iii
		  mm%back3coo(1,ipp+iii) = ki
		  mm%back3coo(2,ipp+iii) = kj
		  mm%back3coo(3,ipp+iii) = l
		  mm%back3coo(4,ipp+iii) = iii
		end do
	      end do
	    end do
	  end do
	end do

	do k=1,nkn
	  ipp = mm%nt3g(k-1) + 1
	        kijc = (k-1)*(nlv+2) + 1
	          mm%i3coo(ipp) = kijc
	          mm%j3coo(ipp) = kijc
		  mm%back3coo(1,ipp) = k
		  mm%back3coo(2,ipp) = k
		  mm%back3coo(3,ipp) = 0
		  mm%back3coo(4,ipp) = 0
	  ipp = mm%nt3g(k)
	        kijc = (k-1)*(nlv+2) + nlv + 2
	          mm%i3coo(ipp) = kijc
	          mm%j3coo(ipp) = kijc
		  mm%back3coo(1,ipp) = k
		  mm%back3coo(2,ipp) = k
		  mm%back3coo(3,ipp) = nlv + 1
		  mm%back3coo(4,ipp) = 0
	end do

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	return
   91	continue
	write(6,*) n3zero,mm%nt3g(nkn),n3max
	stop 'error stop coo_init_new: internal error (5)'
   96	continue
	write(6,*) ipp,n3zero
	stop 'error stop coo_init_new: internal error (4)'
   97	continue
	write(6,*) ipp,n2zero,i,j,ki,kj,n,m,nodes(1:n)
	stop 'error stop coo_init_new: internal error (3)'
   98	continue
	write(6,*) n,kj,nodes(1:n)
	stop 'error stop coo_init_new: internal error (2)'
   99	continue
	write(6,*) n,k,nodes(1:n)
	stop 'error stop coo_init_new: internal error (1)'
	end

!******************************************************************

	subroutine coo_adjust_3d(matrix)

! adjusts zero values in coo matrix (used for 3d matrix)

	use mod_system

	implicit none

        type(smatrix), target :: matrix

	integer ie,ii,k,i,j,l,ipp
	integer nlv,nkn,nel
	integer kn(3)

	integer loccoo3d

        type(smatrix), pointer :: mm

        mm => matrix

!------------------------------------------------------------------
! adjusts zero values in case not all layers are in system
!------------------------------------------------------------------

	nlv = mm%nlv_system
	nkn = mm%nkn_system
	nel = mm%nel_system

	do ie=1,nel
	  do ii=1,3
	    kn(ii) = mm%nen3v_system(ii,ie)
	  end do
	  do i=1,3
	      do l=1,nlv
		ipp = loccoo3d(i,i,kn,l,ie)
	        if( mm%c3coo(ipp) == 0. ) then
		  mm%c3coo(ipp) = 1.
		end if
	      end do
	  end do
	end do

!------------------------------------------------------------------
! adjusts zero values for layer 0 and nlv+1
!------------------------------------------------------------------

	do k=1,nkn
	  ipp = mm%nt3g(k-1) + 1
	  mm%c3coo(ipp) = 1.
	  ipp = mm%nt3g(k)
	  mm%c3coo(ipp) = 1.
	end do

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	end

!******************************************************************

	subroutine construct_node_list(nkn,nel,ngr,nen3v,nlist)

! makes node list around nodes with central node inserted - list is ordered

	implicit none

	integer nkn,nel,ngr
	integer nen3v(3,nel)
	integer nlist(0:2*ngr+1,nkn)	!index 0 is total number of nodes

	integer ie,ii,i2,k,k1,k2,n,i,ip
	integer nodes(2*ngr+1)

	nlist = 0

!----------------------------------------------------
! insert nodes in list - inner nodes are inserted twice
!----------------------------------------------------

	do ie=1,nel
	  do ii=1,3
	    i2 = mod(ii,3) + 1
	    k1 = nen3v(ii,ie)
	    k2 = nen3v(i2,ie)
	    n = nlist(0,k1) + 1
	    nlist(n,k1) = k2
	    nlist(0,k1) = n
	    n = nlist(0,k2) + 1
	    nlist(n,k2) = k1
	    nlist(0,k2) = n
	  end do
	end do

!----------------------------------------------------
! add proper node to list
!----------------------------------------------------

	do k=1,nkn
	  n = nlist(0,k) + 1
	  if( n > 2*ngr+1 ) goto 99
	  nlist(n,k) = k
	  nlist(0,k) = n
	end do

!----------------------------------------------------
! sort and make unique
!----------------------------------------------------

	do k=1,nkn
	  n = nlist(0,k)
	  nodes(1:n) = nlist(1:n,k)
	  call sort_int(n,nodes)
	  ip = 1
	  do i=2,n
	    if( nodes(i) == nodes(ip) ) cycle
	    ip = ip + 1
	    if( ip /= i ) nodes(ip) = nodes(i)
	  end do
	  n = ip
	  if( n > ngr+1 ) goto 98
	  nlist(0,k) = n
	  nlist(1:n,k) = nodes(1:n)
	  nlist(n+1:,k) = 0
	end do

!----------------------------------------------------
! end of routine
!----------------------------------------------------

	return
   98	continue
	write(6,*) 'internal error: n,ngr: ',n,ngr
	stop 'error stop construct_node_list: n > ngr+1'
   99	continue
	write(6,*) 'internal error: n,ngr: ',n,ngr
	stop 'error stop construct_node_list: n > 2*ngr+1'
	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine coo_init(nnel,nnkn,mmbw,nnen3v,ndim,nnot0,iijp,ip,jp)

! construct pointers for coo matrix format
!
! old routine - not used anymore

	implicit none

	integer nnel,nnkn,mmbw		!size of elements, nodes, bandwidth
	integer nnen3v(3,1)		!element index
	integer ndim			!dimension of arrays ip, jp (non zero)
	integer nnot0			!number of non 0 elements (return)
	integer iijp(-mmbw:mmbw,nnkn)	!index pointer into matrix (return)
	integer ip(ndim)		!row index of non zero element (return)
	integer jp(ndim)		!col index of non zero element (return)

	integer k,m,n,ie,ii,ii1,k1,k2

	n = 0
	iijp = 0

	do ie=1,nnel
	  do ii=1,3
	    k1 = nnen3v(ii,ie)
	    ii1 = mod(ii,3) + 1
	    k2 = nnen3v(ii1,ie)
	    call coo_init_insert(k1,k2,nnkn,mmbw,iijp,n)	!out of diagonal
	    call coo_init_insert(k1,k1,nnkn,mmbw,iijp,n)	!in diagonal
	  end do
	end do

	if( n .gt. ndim) goto 99

	nnot0 = n

	do k=1,nnkn
	  do m=-mmbw,mmbw
	    n = iijp(m,k)
	    if( n .ne. 0 ) then
	      jp(n) = k
	      ip(n) = k + m
	    end if
	  end do
	end do

	!call coo_debug(nnkn,mmbw,iijp)
	call coo_check(nnkn,mmbw,iijp,ip,jp)

	return
   99	continue
	write(6,*) 'nnot0,ndim: ',nnot0,ndim
	stop 'error stop coo_init: ndim'
	end

!******************************************************************

	subroutine coo_find(i,j,mmbw,iijp,n)

! finds position of non zero element in arrays
!
! old routine - not used anymore

	implicit none

	integer i,j			!row and col
	integer mmbw
	integer iijp(-mmbw:mmbw,1)
	integer n			!position

	n = iijp(i-j,j)

	end

!******************************************************************

	subroutine coo_debug(nnkn,mmbw,iijp)

! checks sanity of iijp, ip, jp arrays
!
! old routine - not used anymore

	implicit none

	integer nnkn,mmbw
	integer iijp(-mmbw:mmbw,1)

	integer k,m,n,nn,nmax

	n = 0
	nn = 0
	nmax = 0

	do k=1,nnkn
	  do m=-mmbw,mmbw
	    n = iijp(m,k)
	    nmax = max(n,nmax)
	    if( n .gt. 0 ) then
		    nn = nn + 1
	    	    !write(6,*) 'ggg: ',m,k,n
		    if( iijp(-m,k+m) .eq. 0 ) then !check symmetric element
			    write(6,*) m,k,nn,n,iijp(-m,k+m)
			    stop 'error stop coo_debug: internal error (7)'
		    end if
	    end if
	  end do
	end do

	if( nn .ne. nmax ) then
		write(6,*) 'coo_debug: ',nmax,nn
		stop 'error stop coo_debug: nn /= nmax'
	end if

	write(6,*) 'coo_debug: ',nmax,nn

	end

!******************************************************************

	subroutine coo_check(nnkn,mmbw,iijp,ip,jp)

! checks sanity of iijp, ip, jp arrays
!
! not used anymore, only called by not used routine

	implicit none

	integer nnkn,mmbw
	integer iijp(-mmbw:mmbw,1)
	integer ip(1)
	integer jp(1)

	integer i,j,k,m,n
	integer n0
	logical debug

	debug = .true.
	n0 = 0
	n = 0
	j = 0
	i = 0

	do k=1,nnkn
	  do m=-mmbw,mmbw
	    n = iijp(m,k)
	    if( n .gt. 0 ) then
	      n0 = n0 + 1
	      j = k
	      i = j + m
	      if( ip(n) .ne. i ) goto 99
	      if( jp(n) .ne. j ) goto 99
	    end if
	  end do
	end do

	if( debug ) then
	  write(6,*) 'coo routine setup:'
	  write(6,*) '  non zeros   = ',n0
	  write(6,*) '  real band   = ',nnkn + mmbw * (2*nnkn-mmbw-1)
	  write(6,*) '  full band   = ',(1+2*mmbw) * nnkn
	  write(6,*) '  full matrix = ',nnkn*nnkn
	end if

	return
   99	continue
	write(6,*) n,k,m
	write(6,*) i,ip(n)
	write(6,*) j,jp(n)
	stop 'error stop coo_check: internal error'
	end

!******************************************************************

	subroutine coo_init_insert(k1,k2,nnkn,mmbw,ip,n)

! internal routine for insertion of non 0 elements
!
! not used anymore, only called by not used routine

	implicit none

	integer k1,k2,n,nnkn,mmbw
	integer ip(-mmbw:mmbw,nnkn)

	integer idk    
	integer ip1,ip2

	idk = k1 - k2

	if( idk .ne. 0 ) then		!out of diagonal
	  ip1 = ip(idk,k2)
	  ip2 = ip(-idk,k2+idk)
	  !write(6,*) k1,k2,idk,n
	  if( ip1 .ne. 0 .and. ip2 .ne. 0 ) then	!already inserted
	    return
	  else if( ip1 .eq. 0 .and. ip2 .eq. 0 ) then	!must insert
	    ip(idk,k2) = n + 1
	    ip(-idk,k2+idk) = n + 2
	    n = n + 2
	  else						!not possible
	    write(6,*) k1,k2,idk,n
	    write(6,*) 'internal error: ',ip1,ip2
	    stop 'error stop coo_insert: internal error (1)'
	  end if
	else				!on diagonal
	  !write(6,*) k1,k2,idk,n
	  ip1 = ip(0,k1)
	  if( ip1 .eq. 0 ) then
	    n = n + 1
	    ip(0,k1) = n
	  end if
	end if

	end
	  
!******************************************************************

	function loccoo(i,j,nnkn,mmbw)

! localize for COO routines
!
! not used anymore, only called by not used routine

	implicit none

	integer loccoo			!position
	integer i,j			!row and col
	integer nnkn			!size of system
	integer mmbw			!bandwidth

	integer ip,loccoo1

	stop 'error stop loccoo: do not call'

	!loccoo = ijp(i-j,j)
	!ip = j*(2*mmbw)+i-mmbw

	!ip = j*(2*mmbw+1)+i-j-mmbw
	!loccoo=ijp(ip)
	loccoo = 0

	!call coo_find(i,j,mmbw,ijp,loccoo1)

	!if( loccoo .ne. loccoo1) then
	!  write(6,*) 'loccoo...',i,j,mmbw,loccoo,loccoo1
	!  stop 'error stop loccoo'
	!end if

	end

!******************************************************************
!******************************************************************
!******************************************************************

	function loccoo3d(i,j,kn,l,ie)

! localize for COO routines (3d version)

	use mod_system

	implicit none

	integer loccoo3d		!position
	integer i,j,l,ie		!row and col
	integer kn(3)

	integer ki,idiag,irel,icorr
	integer ipp1,ipp2,ipp3

        type(smatrix), pointer :: mm

        mm => l_matrix

	ki = kn(i)
	ipp1 = mm%nt3g(ki-1)				!before row i
	ipp2 = 1 + (l-1) * (mm%ng(ki)+2)		!in row i before level l
	idiag = mm%diag(ki) + 1
	irel = mm%ijp_ie(i,j,ie) - mm%ijp_ie(i,i,ie)	!diff kj - diag
	icorr = 0
	if( irel /= 0 ) icorr = irel/abs(irel)		!icorr is -1,0,+1
	ipp3 = idiag + irel + icorr

	loccoo3d = ipp1 + ipp2 + ipp3

	end

!******************************************************************

