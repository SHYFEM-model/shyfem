
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
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

! revision log :
!
! 07.12.2017	ggu	changed VERS_7_5_40
! 24.01.2018	ggu	changed VERS_7_5_41
! 22.02.2018	ggu	changed VERS_7_5_42
! 19.04.2018	ggu	changed VERS_7_5_45
! 26.04.2018	ggu	changed VERS_7_5_46
! 11.05.2018	ggu	changed VERS_7_5_47
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 05.04.2022	ggu	new routine check_global_indices()
! 12.04.2022	ggu	possibility to partition online now
! 10.03.2023	ggu	new call to ghost_debug()
! 20.03.2023	ggu	ghost node calls transferred to submpi_ghost.f
! 20.03.2023	ggu	call to write_grd_domain() at end of shympi_setup()
! 12.04.2023	ggu	protect from nkn == nel

! notes :
!
! shympi_setup
!
!	handle_partition
!	make_domain
!
!      	make_index(my_id,nkn,n_lk,nodes,nindex)
!	make_index(my_id,nel,n_le,elems,eindex)
!      	shympi_alloc_id(n_lk,n_le)
!      	shympi_alloc_sort(n_lk,n_le)
!      	adjust_indices(n_lk,n_le,nodes,elems,nindex,eindex)
!
!	transfer_domain
!
!	ghost_handle
!       	ghost_make         !here also call to shympi_alloc_ghost()
!       	ghost_check
!       	ghost_debug
!       	ghost_write
!
!       	ghost_exchange
!
!       shympi_univocal_nodes
!
!*****************************************************************

	subroutine shympi_setup

! this sets up the single domains
! it should be called right after shympi_init
! it does nothing if there is only one thread

	use basin
	use shympi

	implicit none

	integer k,ie,i,nc
	integer icust
	integer n,nn
	integer n_lk,n_le
	integer nkn_tot,nel_tot
	integer nodes(nkn)
	integer elems(nel)
	integer nindex(nkn)
	integer eindex(nel)
	integer area_node(nkn)
	integer ivals(n_threads)

	if( .not. bmpi ) return

	if( shympi_partition_on_elements() ) then
	  stop 'error stop shympi_setup: cannot yet partition on elems'
	end if

	if( shympi_is_master() ) then
	  write(6,*) 'setting up mpi with number of threads: ',n_threads
	end if

!      -----------------------------------------------------
!      check if nkn == nel (global)
!      -----------------------------------------------------

       if( nkn /= nkn_global ) stop 'error stop nkn /= nkn_global'
       if( nel /= nel_global ) stop 'error stop nel /= nel_global'
       if( nkn == nel ) then
         write(6,*) 'my_id,nkn,nel: ',my_id,nkn,nel
         stop 'error stop nkn == nel (global)'
       end if

!	-----------------------------------------------------
!	-----------------------------------------------------
!	do partitioning
!	-----------------------------------------------------

	call handle_partition(area_node)

!	=====================================================================
!	the next call is custom call
!	sets array area_node(), with values from 0 to n_threads-1
!	was used only for testing ... not used anymore
!	=====================================================================

!	call make_custom_domain_area(area_node)

!	=====================================================================
!	from here on everything is general
!	=====================================================================

	nc = maxval(area_node)
	if( nc+1 /= n_threads ) then
	  if( shympi_is_master() ) then
	    write(6,*) 'number of threads = ',n_threads
	    write(6,*) 'number of domains = ',nc+1
	  end if
	  call shympi_stop('error stop: thread/domain mismatch')
	end if

!	-----------------------------------------------------
!	set up domain
!	-----------------------------------------------------

	call make_domain(my_id,area_node,nodes,elems,nc)

	call make_index(my_id,nkn,n_lk,nodes,nindex)
	call make_index(my_id,nel,n_le,elems,eindex)
	call shympi_alloc_id(n_lk,n_le)
	call shympi_alloc_sort(n_lk,n_le)
	call adjust_indices(n_lk,n_le,nodes,elems,nindex,eindex)

!	-----------------------------------------------------
!	debug output
!	-----------------------------------------------------

	if( my_unit > 0 ) then
	  write(my_unit,*) 'nodes in domain: ',nkn_local
	  write(my_unit,'(10i7)') (nindex(i),i=1,nkn_local)
	  write(my_unit,*) 'elems in domain: ',nel_local
	  write(my_unit,'(10i7)') (eindex(i),i=1,nel_local)
	end if

!	-----------------------------------------------------
!	transfers global domain to local domain
!	-----------------------------------------------------

	call transfer_domain(nkn_local,nel_local,nindex,eindex)

!	-----------------------------------------------------
!	debug output
!	-----------------------------------------------------

	if( my_unit > 0 ) then
	  write(my_unit,*) 'my_id: ',my_id
	  write(my_unit,*) 'nkn,nel :',nkn,nel
	  write(my_unit,*) 'nkn_inner,nel_inner :',nkn_inner,nel_inner
	  write(my_unit,*) 'nkn_unique,nel_unique :',nkn_unique,nel_unique
	  do k=1,nkn
	    write(my_unit,*) k,ipv(k),id_node(k)
	  end do
	  do ie=1,nel
	    write(my_unit,*) ie,ipev(ie),id_elem(:,ie)
	  end do
	end if

!	-----------------------------------------------------
!	set up ghost nodes
!	-----------------------------------------------------

	call ghost_handle

!	-----------------------------------------------------
!	other stuff
!	-----------------------------------------------------

        call shympi_univocal_nodes

!	-----------------------------------------------------
!	gather info on domains
!	-----------------------------------------------------

	call shympi_gather(nkn_local,ivals)
	nkn_domains = ivals
	call shympi_gather(nel_local,ivals)
	nel_domains = ivals

	call shympi_gather(nkn_unique,ivals)
	nkn_domains_u = ivals
	call shympi_gather(nel_unique,ivals)
	nel_domains_u = ivals

	nk_max = maxval(nkn_domains)
	ne_max = maxval(nel_domains)
	nn_max = max(nk_max,ne_max)

	do i=1,n_threads
	  nkn_cum_domains(i) = nkn_cum_domains(i-1) + nkn_domains(i)
	  nel_cum_domains(i) = nel_cum_domains(i-1) + nel_domains(i)
	end do

!	-----------------------------------------------------
!	write to terminal and final check
!	-----------------------------------------------------

	call shympi_syncronize

	nkn_tot = shympi_sum(nkn_unique)
	nel_tot = shympi_sum(nel_unique)

	write(6,*) 'info on nkn/nel: ',my_id
	write(6,*) 'nkn: ',nkn_global,nkn_local,nkn_unique,nkn_inner
	write(6,*) 'nel: ',nel_global,nel_local,nel_unique,nel_inner
	write(6,*) 'tot: ',nkn_tot,nel_tot

	if( nkn_global /= nkn_tot .or. nel_global /= nel_tot ) then
	  write(6,*) nkn_global,nkn_tot
	  write(6,*) nel_global,nel_tot
	  stop 'error stop shympi_setup: internal error (5)'
	end if

	call shympi_syncronize

!      -----------------------------------------------------
!      check if nkn == nel (local)
!      -----------------------------------------------------

       if( nkn /= nkn_local ) stop 'error stop nkn /= nkn_local'
       if( nel /= nel_local ) stop 'error stop nel /= nel_local'
       if( nkn == nel ) then
         write(6,*) 'my_id,nkn,nel: ',my_id,nkn,nel
         !stop 'error stop nkn == nel (local)'
       end if

!	-----------------------------------------------------
!	write domain*.grd files
!	-----------------------------------------------------

	if( bmpi_debug ) call write_grd_domain

!	-----------------------------------------------------
!	end of routine
!	-----------------------------------------------------

	!call shympi_stop('forced stop in shympi_setup')

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine make_domain(my_id,n_area,nodes,elems,nc)

! sets nodes, elems, and nc
!
! nodes is my_id for proper nodes, >= 0 for ghost nodes and -1 else
! elems is my_id for proper elems, -2 for border elems and -1 else
!
! nkn and nel are still global values

	use basin

	implicit none

	integer my_id
	integer n_area(nkn)
	integer nodes(nkn)	!nodes information (return)
	integer elems(nel)	!elems information (return)
	integer nc		!maximum number of colors (return)

	integer ie,k,ii
	integer n,n_my,n_ghost
	integer, allocatable :: ncs(:)

	nodes = -1
	elems = -1

!	-----------------------------------------------------
!	set up ncs - number of nodes in areas
!	-----------------------------------------------------

	nc = maxval(n_area)
	allocate(ncs(0:nc))
	ncs = 0

	if( my_id < 0 .or. my_id > nc ) then
	  write(6,*) 'my_id,nc : ',my_id,nc
	  stop 'error stop make_domain: my_id out of range'
	end if

	do k=1,nkn
	  n = n_area(k)
	  if( n < 0 .or. n > nc ) then
	    write(6,*) 'n_area = ',n,' in node ',k
	    stop 'error stop make_domain: out of range'
	  end if
	  ncs(n) = ncs(n) + 1
	end do

	do n=0,nc
	  if( ncs(n) == 0 ) then
	    write(6,*) 'n_area = ',n,' has 0 nodes'
	    stop 'error stop make_domain:  no nodes'
	  end if
	end do

!	-----------------------------------------------------
!	flag elements - my_id if internal, -2 if border, else -1
!	-----------------------------------------------------

	do ie=1,nel
	  n = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( n_area(k) == my_id ) n = n + 1
	  end do
	  if( n == 3 ) then		!internal elem
	    elems(ie) = my_id
	  else if( n > 0 ) then		!border elem
	    elems(ie) = -2
	    do ii=1,3
	      k = nen3v(ii,ie)
	      nodes(k) = -2		!flag temporarily
	    end do
	  else
	    !nothing, leave -1
	  end if
	end do

!	-----------------------------------------------------
!	flag nodes - my_id if internal, area if border, else -1
!	-----------------------------------------------------

	do k=1,nkn
	  if( n_area(k) == my_id ) then
	    nodes(k) = n_area(k)
	  else if( nodes(k) == -2 ) then
	    nodes(k) = n_area(k)
	  else
	    !nothing, leave -1
	  end if
	end do

!	-----------------------------------------------------
!	count proper and ghost nodes
!	-----------------------------------------------------

	n_my = 0
	n_ghost = 0
	do k=1,nkn
	  if( nodes(k) == my_id ) then		!proper node
	    n_my = n_my + 1
	  else if( nodes(k) >= 0 ) then		!ghost node
	    n_ghost = n_ghost + 1
	  else if( nodes(k) == -1 ) then	!any other node
	    !nothing
	  else
	    write(6,*) nodes(k)
	    stop 'error stop make_domain:  impossible value for nodes'
	  end if
	end do

!	-----------------------------------------------------
!	final checks
!	-----------------------------------------------------

	n = ncs(my_id)
	if( n /= n_my ) then
	  write(6,*) n,n_my,n_ghost,n_my+n_ghost
	  stop 'error stop make_domain:  internal error (2)'
	end if

	!write(6,*) 'domain = ',my_id,n_my,n_ghost

	deallocate(ncs)

!	-----------------------------------------------------
!	end of routine
!	-----------------------------------------------------

	end

!*****************************************************************

	subroutine adjust_indices(n_lk,n_le
     +				,nodes,elems,nindex,eindex)

! computes nkn_local/unique/inner and nel_local/unique/inner
! also rearranges eindex to keep track of this
! nen3v is still global

!--------------------------------------------------------------
!
! id_node: 
!
!    either my_id for inner node or other for ghost node
!
! id_elem:
!
!    (1,my_id,-1,-1) for inner
!    (2,my_id,id,-1) for border to color id (2 nodes my_id, 1 node id)
!    (3,my_id,id1,id2) for border to two colors id1 and id2 (tripple point)
!
! if main domain for element is not my_id, then here are the possibilities:
!
!    (2,id,my_id,-1) for two nodes id, one node my_id
!    (3,id1,my_id,id2) for triple point, main color is id1
!    (3,id1,id2,my_id) for triple point, main color is id1
!
! impossible constellation:
!
!    (2,my_id,id,id) -> must be (2,my_id,id,-1)
!    (3,my_id,id,id) -> must be (3,my_id,id1,id2)
!
! in id_elem(1) is always the main domain of the element
!
!--------------------------------------------------------------

	use basin
	use shympi

	implicit none

	integer n_lk,n_le
	integer nodes(nkn)	!nkn is still global
	integer elems(nel)
	integer nindex(n_lk)
	integer eindex(n_le)

	integer ie,k,ii,n,i
	integer it,is
	integer is1,is2,k1,k2,kmin
	integer iu,id
	integer n1,n2
	integer iunique(n_le)
	integer idiff(n_le)
	integer id_aux(0:3,nel)		!total elements
	logical bthis,bexchange,bwrite

	bexchange = .false.
	bexchange = .true.

!	-----------------------------------------------------
!	deal with node index
!	-----------------------------------------------------

	nkn_local = n_lk
	id_node = -1

	do i=1,nkn_local
	  k = nindex(i)
	  n = nodes(k)
	  if( n /= my_id ) exit
	  id_node(i) = n
	end do

	nkn_inner = i-1
	nkn_unique = nkn_inner
	
	do i=nkn_inner+1,nkn_local
	  k = nindex(i)
	  n = nodes(k)
	  id_node(i) = n
	  if( n == my_id ) then
	    stop 'error stop adjust_indices: internal error (7)'
	  end if
	end do

!	-----------------------------------------------------
!	deal with elem index
!	-----------------------------------------------------

	nel_local = n_le
	id_aux = -1

	do i=1,nel_local
	  ie = eindex(i)
	  n = elems(ie)		!is -2 if border element
	  if( n /= my_id ) exit
	  id_aux(0,ie) = 1
	  id_aux(1,ie) = n
	end do

	nel_inner = i-1

	do i=nel_inner+1,nel_local
	  ie = eindex(i)
	  n = elems(ie)
	  if( n == my_id ) then
	    stop 'error stop adjust_indices: internal error (8)'
	  end if
	end do

!	-----------------------------------------------------
!	compute unique elem index
!	-----------------------------------------------------

	iu = 0
	id = 0
	!in = 0
	do i=nel_inner+1,nel_local
	  ie = eindex(i)
	  it = 0
	  is = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( nodes(k) == my_id ) then
	      is = is + ii
	      it = it + 1
	    end if
	  end do
	  if( it == 2 ) then			!two nodes with my_id
	    id_aux(0,ie) = 2
	    id_aux(1,ie) = my_id
	    is = 6 - is
	    k = nen3v(is,ie)
	    n = nodes(k)
	    id_aux(2,ie) = n
	  else if( it == 1 ) then		!one node only with my_id
	    is1 = mod(is,3) + 1
	    is2 = mod(is1,3) + 1
	    k1 = nen3v(is1,ie)
	    k2 = nen3v(is2,ie)
	    n1 = nodes(k1)
	    n2 = nodes(k2)
	    if( n1 /= n2 ) then			!all three nodes are different
	      id_aux(0,ie) = 3
	      kmin = minval(nen3v(:,ie))
	      k = nen3v(is,ie)
	      if( kmin == k ) then
	        id_aux(1,ie) = my_id
	        id_aux(2,ie) = n1
	        id_aux(3,ie) = n2
	      else if( kmin == k1 ) then
	        id_aux(1,ie) = n1
	        id_aux(2,ie) = my_id
	        id_aux(3,ie) = n2
	      else
	        id_aux(1,ie) = n2
	        id_aux(2,ie) = my_id
	        id_aux(3,ie) = n1
	      end if
	    else				!two nodes of other color
	      n = nodes(k1)
	      id_aux(0,ie) = 2
	      id_aux(1,ie) = n
	      id_aux(2,ie) = my_id
	    end if
	  else
	    stop 'error stop adjust_indices: internal error (1)'
	  end if
	  if( id_aux(1,ie) == my_id ) then
	    iu = iu + 1
	    iunique(iu) = ie
	  else
	    id = id + 1
	    idiff(id) = ie
	  end if
	end do

	if( nel_local-nel_inner /= iu+id ) then
	  stop 'error stop adjust_indices: internal error (2)'
	end if

	nel_unique = nel_inner + iu

	if( nel_unique + id /= nel_local ) then
	  stop 'error stop adjust_indices: internal error (3)'
	end if

	!write(6,*) 'eindex: ',my_id,eindex(1:nel_inner)
	!write(6,*) 'eindex before: ',my_id,eindex(nel_inner+1:nel_local)

	if( bexchange ) then
	  do i=1,nel_inner
	    ie = eindex(i)
	    id_elem(:,i) = id_aux(:,ie)
	  end do
	  do i=1,iu
	    ie = iunique(i)
	    eindex(nel_inner+i) = ie
	    id_elem(:,nel_inner+i) = id_aux(:,ie)
	  end do
	  do i=1,id
	    ie = idiff(i)
	    eindex(nel_unique+i) = ie
	    id_elem(:,nel_unique+i) = id_aux(:,ie)
	  end do
	end if

!	-----------------------------------------------------
!	check computed element index
!	-----------------------------------------------------

	do i=1,nel_local
	  n = id_elem(0,i)
	  if( .not. ( any(id_elem(1:3,i) == my_id ) ) ) goto 99
	  if( n == 1 ) then
	    if( id_elem(1,i) /= my_id ) goto 99
	    if( any(id_elem(2:3,i) /= -1 ) ) goto 99
	  else if( n == 2 ) then
	    if( any(id_elem(1:2,i) == -1 ) ) goto 99
	    if( id_elem(3,i) /= -1 ) goto 99
	  else if( n == 3 ) then
	    if( any(id_elem(1:3,i) == -1 ) ) goto 99
	  else
	    goto 99
	  end if
	end do

	iu = 400+my_id
	iu = 0
	if( iu > 0 ) then
	write(iu,*) '========================================'
	write(iu,*) 'adjust_indices'
	write(iu,*) '========================================'
	write(iu,*) 'item index for color',my_id
	write(iu,*) 'border nodes',nkn_inner,nkn_unique,nkn_local
	write(iu,*) 'inner nodes'
	bwrite = .false.
	do k=1,nkn_local
	  if( k > nkn_inner ) bwrite = .true.
	  if( k <= nkn_inner .and. id_node(k) /= my_id ) goto 98
	  if( bwrite ) write(iu,*) k,id_node(k)
	  if( k == nkn_inner ) write(iu,*) 'unique nodes'
	  if( k == nkn_unique ) write(iu,*) 'outer nodes'
	end do
	write(iu,*) 'border elements',nel_inner,nel_unique,nel_local
	write(iu,*) 'inner elements'
	bwrite = .false.
	do ie=1,nel_local
	  if( ie > nel_inner ) bwrite = .true.
	  if( ie <= nel_inner ) then
	    if( id_elem(1,ie) /= my_id .or. 
     +			any(id_elem(2:3,ie)/=-1) ) goto 98
	  end if
	  if( bwrite ) write(iu,*) ie,id_elem(:,ie)
	  if( ie == nel_inner ) write(iu,*) 'unique elements'
	  if( ie == nel_unique ) write(iu,*) 'outer elements'
	end do
	flush(iu)
	end if

	!write(my_unit,*) 'debug id_elem: ',my_id
	!do i=1,nel_local
	!  write(my_unit,*) i,eindex(i),id_elem(:,i)
	!end do

	call shympi_syncronize
	!write(6,*) 'finished running adjust_indices'
	!flush(6)

!	-----------------------------------------------------
!	end of routine
!	-----------------------------------------------------

	return
   98	continue
	stop 'error 400'
   99	continue
	write(6,*) 'error in element index: '
	write(6,*) my_id,nel_inner,nel_unique,nel_local,nel
	write(6,*) i,id_elem(:,i)
	stop 'error stop adjust_indices: internal error (7)'
	end

!*****************************************************************

	subroutine make_index(my_id,n_g,n_l,items,index)

! returns index of items, sorted first by proper and then ghost items

	implicit none

	integer my_id
	integer n_g		!global number of items
	integer n_l		!local number of items (return)
	integer items(n_g)	!color of items
	integer index(n_g)	!index of items (return)

	integer i,k,n

	i = 0

!	-----------------------------------------------------
!	first lists inner (proper) items
!	-----------------------------------------------------

	do k=1,n_g
	  n = items(k)
	  if( n == my_id ) then
	    i = i + 1
	    index(i) = k
	  end if
	end do

!	-----------------------------------------------------
!	now lists ghost items
!	-----------------------------------------------------

	do k=1,n_g
	  n = items(k)
	  if( n /= my_id .and. n /= -1 ) then
	    i = i + 1
	    index(i) = k
	  end if
	end do

!	-----------------------------------------------------
!	end of routine
!	-----------------------------------------------------

	n_l = i

	end

!*****************************************************************

	subroutine transfer_domain(n_lk,n_le,nindex,eindex)

! transfers global domain to local domain
!
! nindex/eindex contain global internal node/elem numbers -> save

	use basin
	use shympi

	implicit none

	integer n_lk,n_le,nmax
	integer nindex(n_lk)
	integer eindex(n_le)

	integer ie,k,i,ii,kk
	integer inverse(nkn)		!aux array for check

        integer, allocatable :: nen3v_aux(:,:)
        integer, allocatable :: ipev_aux(:)
        integer, allocatable :: ipv_aux(:)
        integer, allocatable :: iarv_aux(:)
        integer, allocatable :: iarnv_aux(:)

        real, allocatable :: xgv_aux(:)
        real, allocatable :: ygv_aux(:)
        real, allocatable :: hm3v_aux(:,:)

!	----------------------------------
!	allocate aux arrays
!	----------------------------------

	allocate(nen3v_aux(3,n_le))
	allocate(ipev_aux(n_le))
	allocate(ipv_aux(n_lk))
	allocate(iarv_aux(n_le))
	allocate(iarnv_aux(n_lk))

	allocate(xgv_aux(n_lk))
	allocate(ygv_aux(n_lk))
	allocate(hm3v_aux(3,n_le))

!	----------------------------------
!	set up inverse information
!	----------------------------------

	inverse = 0
	do i=1,n_lk
	  k = nindex(i)
	  inverse(k) = i
	end do

!	----------------------------------
!	create and set up auxiliary arrays
!	----------------------------------

	do i=1,n_le
	  ie = eindex(i)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    kk = inverse(k)
	    if( kk <= 0 ) then
	      write(6,*) ie,k,kk
	      stop 'error stop transfer_domain: internal error (3)'
	    end if
	    nen3v_aux(ii,i) = kk
	  end do
	  hm3v_aux(:,i) = hm3v(:,ie)
	  ipev_aux(i) = ipev(ie)
	  iarv_aux(i) = iarv(ie)
	end do

	do i=1,n_lk
	  k = nindex(i)
	  ipv_aux(i) = ipv(k)
	  iarnv_aux(i) = iarnv(k)
	  xgv_aux(i) = xgv(k)
	  ygv_aux(i) = ygv(k)
	end do

!	----------------------------------
!	create new basin and copy
!	----------------------------------

	call basin_init(n_lk,n_le)

	nen3v = nen3v_aux
	ipev = ipev_aux
	ipv = ipv_aux
	iarv = iarv_aux
	iarnv = iarnv_aux
	xgv = xgv_aux
	ygv = ygv_aux
	hm3v = hm3v_aux

!	----------------------------------
!	deallocate aux arrays
!	----------------------------------

	deallocate(nen3v_aux)
	deallocate(ipev_aux)
	deallocate(ipv_aux)
	deallocate(iarv_aux)
	deallocate(iarnv_aux)

	deallocate(xgv_aux)
	deallocate(ygv_aux)
	deallocate(hm3v_aux)

!	----------------------------------
!	save global internal node/elem numbers
!	----------------------------------

	deallocate(ip_int_node,ip_int_elem)
	allocate(ip_int_node(n_lk),ip_int_elem(n_le))
	ip_int_node = nindex
	ip_int_elem = eindex

!	----------------------------------
!	make pointer from local to global arrays
!	----------------------------------

	call make_intern2global_index(n_lk,n_le)
	
!	----------------------------------
!	sort index to nodes and elems
!	----------------------------------

	call mpi_sort_index(n_lk,n_le)

!	----------------------------------
!	final check
!	----------------------------------

	call check_global_indices(n_lk,n_le)

!	----------------------------------
!	end routine
!	----------------------------------

	write(6,*) 'finished transfer_domain'

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine make_intern2global_index(n_lk,n_le)

	use basin
	use shympi

	implicit none

	integer n_lk,n_le,nmax

	integer ie,k,i,ii,kk
	integer ia,iext,ip,n,iiext
	integer iwhat

        integer, allocatable :: index(:)
        integer, allocatable :: ipaux(:)
        integer, allocatable :: ipextern(:,:)

	integer locate

!	----------------------------------
!	compute max values which are still not available
!	----------------------------------

	!write(6,*) 'n_lk,n_le',n_lk,n_le
	!write(6,*) 'nkn_domains: ',size(nkn_domains)

        call shympi_gather(n_lk,nkn_domains)
        call shympi_gather(n_le,nel_domains)

	nk_max = maxval(nkn_domains)
	ne_max = maxval(nel_domains)
	nn_max = max(nk_max,ne_max)

!	----------------------------------
!	re-allocate arrays
!	----------------------------------

	deallocate(ip_int_nodes,ip_int_elems)
	allocate(ip_int_nodes(nk_max,n_threads))
	allocate(ip_int_elems(ne_max,n_threads))
	ip_int_nodes = 0
	ip_int_elems = 0

!	----------------------------------
!	create node index
!	----------------------------------

	iwhat = 1

	allocate(index(nkn_global))
	allocate(ipaux(nk_max))
        allocate(ipextern(nk_max,n_threads))
	ipaux = 0
	ipextern = 0
	index = 0

        call isort(nkn_global,ip_ext_node,index)

	ipaux(1:n_lk) = ipv(1:n_lk)
        call shympi_gather(ipaux,ipextern)

	do ia=1,n_threads
	  n=nkn_domains(ia)
	  do i=1,n
	    iext = ipextern(i,ia)
	    ip = locate(nkn_global,ip_ext_node,index,iext)
	    if( ip <= 0 ) goto 99
	    ip_int_nodes(i,ia) = ip
	  end do
	end do

	write(6,*) 'size of ipextern: ',size(ipextern,1),n_lk,my_id
	write(6,*) nk_max,ne_max,nn_max
	write(6,*) size(ipextern,1),size(ipextern,2)

	deallocate(ipaux)
        deallocate(ipextern)
        deallocate(index)

	!call shympi_syncronize
	!stop

!	----------------------------------
!	create elem index
!	----------------------------------

	iwhat = 2

	allocate(index(nel_global))
	allocate(ipaux(ne_max))
        allocate(ipextern(ne_max,n_threads))
	ipaux = 0
	ipextern = 0
	index = 0

        call isort(nel_global,ip_ext_elem,index)

	ipaux(1:n_le) = ipev(1:n_le)
        call shympi_gather(ipaux,ipextern)

	do ia=1,n_threads
	  n=nel_domains(ia)
	  do i=1,n
	    iext = ipextern(i,ia)
	    ip = locate(nel_global,ip_ext_elem,index,iext)
	    if( ip <= 0 ) goto 99
	    ip_int_elems(i,ia) = ip
	  end do
	end do

	deallocate(ipaux)
        deallocate(ipextern)
        deallocate(index)
	
!	----------------------------------
!	local check
!	----------------------------------

	ia = my_id + 1

	do i=1,n_lk
	  iext = ipv(i)
	  ip = ip_int_nodes(i,ia)
	  if( ip <= 0 ) goto 98
	  iiext = ip_ext_node(ip)
	  if( iext /= iiext ) goto 98
	end do

	do i=1,n_le
	  iext = ipev(i)
	  ip = ip_int_elems(i,ia)
	  if( ip <= 0 ) goto 98
	  iiext = ip_ext_elem(ip)
	  if( iext /= iiext ) goto 98
	end do

!	----------------------------------
!	end routine
!	----------------------------------

	!stop 'successfull completion of make_inttern2global_index'

	return
   98	continue
	write(6,*) 'cannot verify item: ',ia,n_le,i,iext,iiext
	stop 'error stop make_inttern2global_index: wrong item'
   99	continue
	write(6,*) 'cannot find item: ',iwhat
	write(6,*) ia,n,i,iext
	write(6,*) 'nkn_domains: ',nkn_domains
	write(6,*) 'nel_domains: ',nel_domains
	write(6,*) n_lk,n_le
	stop 'error stop make_inttern2global_index: ip == 0'
	end

!*****************************************************************

	subroutine check_global_indices(n_lk,n_le)

	use shympi

	implicit none

	integer n_lk,n_le

	logical bstop
	integer ia,id,k,ie
	integer iu,idiff

	iu = 444 + my_id
	idiff = 0
	bstop = .false.

	ia = my_id + 1

	do k=1,n_lk
	  if( ip_int_node(k) /= ip_int_nodes(k,ia) ) then
	    write(iu,*) 'node differences: ',k,ip_int_node(k)
     +				,ip_int_nodes(k,ia)
	    idiff = idiff + 1
	    bstop = .true.
	  end if
	end do

	do ie=1,n_le
	  if( ip_int_elem(ie) /= ip_int_elems(ie,ia) ) then
	    write(iu,*) 'elem differences: ',ie,ip_int_elem(ie)
     +				,ip_int_elems(ie,ia)
	    idiff = idiff + 1
	    bstop = .true.
	  end if
	end do
	
	if( bstop ) then
	  write(6,*) n_lk,n_le,my_id,idiff
	  write(6,*) 'more info in files 444+'
	  stop 'error stop check_global_indices: inconsistency'
	end if

	end

!*****************************************************************

	subroutine handle_partition(area_node)

	use basin
	use shympi

	implicit none

	integer area_node(nkn)

	integer nparts
	integer nnp,nep
	integer nmin,nmax
	integer area_elem(nel)
	integer ierr1,ierr2

	ierr1 = 0
	ierr2 = 0

	call basin_get_partition(nkn,nel,nnp,nep,area_node,area_elem)

	nnp = nnp + 1
	nparts = n_threads
 
	if( .not. bmpi ) then
	  stop 'error stop handle_partition: internal error (1)'
	end if

	if( nnp == 1 ) then
	  if( shympi_is_master() ) then
	    write(6,*) 'no partitiones contained in basin...'
	    write(6,*) 'we will do partitioning for domains: ',nparts
	  end if
	  call do_partition(nkn,nel,nen3v,nparts,area_node,area_elem)
	  call check_partition(area_node,area_elem,ierr1,ierr2)
	  area_node = area_node - 1	!gives back 1-nparts
	else if( nnp == nparts ) then
	  if( shympi_is_master() ) then
	    write(6,*) 'partitiones contained in basin: ',nnp
	    write(6,*) 'using these partitiones...'
	  end if
	else
	  if( shympi_is_master() ) then
	    write(6,*) 'partitiones contained in basin: ',nnp
	    write(6,*) 'partitiones required: ',nparts
	    stop 'error stop handle_partition: domains not compatible'
	  end if
	end if

	if( shympi_is_master() ) then
	  nmin = minval(area_node)
	  nmax = maxval(area_node)

	  !write(6,*) nkn,nel
	  !write(6,*) nnp,nep
	  !write(6,*) nmin,nmax
	  write(6,*) 'domains: ',nmin,nmax

	  call info_partition(nparts,area_node)
	end if

	if( ierr1 /= 0 .or. ierr2 /= 0 ) then
	  write(6,*) 'error in partitioning: ',ierr1,ierr2
	  stop 'error stop handle_partition: partitioning error'
	end if

	end

!*****************************************************************

