
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

!*****************************************************************

	subroutine shympi_setup

	use basin
	use shympi

	implicit none

	integer k,ie,i,nc
	integer icust
	integer n,nn
	integer n_lk,n_le
	integer nnp,nep
	integer nkn_tot,nel_tot
	integer nodes(nkn)
	integer elems(nel)
	integer nindex(nkn)
	integer eindex(nel)
	integer area_node(nkn)
	integer area_elem(nel)
	integer vals(n_threads)

	if( .not. bmpi ) return

	if( shympi_partition_on_elements() ) then
	  stop 'error stop shympi_setup: cannot yet partition on elems'
	end if

	if( shympi_is_master() ) then
	  write(6,*) 'setting up mpi with number of threads: ',n_threads
	end if

	call basin_get_partition(nkn,nel,nnp,nep,area_node,area_elem)

!	=====================================================================
!	the next call is custom call
!	sets array area_node(), with values from 0 to n_threads-1
!	=====================================================================

	call make_custom_domain_area(area_node)

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

	call make_domain(my_id,area_node,nodes,elems,nc)

	call make_index(my_id,nkn,n_lk,nodes,nindex)
	call make_index(my_id,nel,n_le,elems,eindex)
	call shympi_alloc_id(n_lk,n_le)
	call shympi_alloc_sort(n_lk,n_le)
	call adjust_indices(n_lk,n_le,nodes,elems,nindex,eindex)


	if( my_unit > 0 ) then
	  write(my_unit,*) 'nodes in domain: ',nkn_local
	  write(my_unit,'(10i7)') (nindex(i),i=1,nkn_local)
	  write(my_unit,*) 'elems in domain: ',nel_local
	  write(my_unit,'(10i7)') (eindex(i),i=1,nel_local)
	end if

	call transfer_domain(nkn_local,nel_local,nindex,eindex)
	!call make_domain_final(area_node,nindex,eindex)

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

	call ghost_make		!here also call to shympi_alloc_ghost()
	call ghost_check
	call ghost_write

	write(6,*) 'mpi my_unit: ',my_unit
	write(6,'(a,9i7)') ' mpi domain: ',my_id,n_ghost_areas
	write(6,'(a,5i7)') 'nkn: '
     +			,nkn_global,nkn_local,nkn_unique
     +			,nkn_inner,nkn_local-nkn_inner
	write(6,'(a,5i7)') 'nel: '
     +			,nel_global,nel_local,nel_unique
     +			,nel_inner,nel_local-nel_inner

	call shympi_syncronize

	call shympi_alloc_buffer(n_ghost_max)	!should probably be 1
	call ghost_exchange

        call shympi_univocal_nodes

!	-----------------------------------------------------
!	exchange info on domains
!	-----------------------------------------------------

	call shympi_gather(nkn_local,vals)
	nkn_domains = vals
	call shympi_gather(nel_local,vals)
	nel_domains = vals

	nk_max = maxval(nkn_domains)
	ne_max = maxval(nel_domains)
	nn_max = max(nk_max,ne_max)

	do i=1,n_threads
	  nkn_cum_domains(i) = nkn_cum_domains(i-1) + nkn_domains(i)
	  nel_cum_domains(i) = nel_cum_domains(i-1) + nel_domains(i)
	end do

	!write(6,*) 'domain gather: ',n_threads,my_id,nkn,nel
	!write(6,*) nkn_domains
	!write(6,*) nkn_cum_domains
	!write(6,*) nel_domains
	!write(6,*) nel_cum_domains
	!call shympi_finalize
	!stop

!	-----------------------------------------------------
!	write to terminal
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

	!call shympi_stop('forced stop in shympi_setup')

!	-----------------------------------------------------
!	end of routine
!	-----------------------------------------------------

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine make_domain(my_id,n_area,nodes,elems,nc)

! sets nodes, elems, and nc
!
! nodes is my_id for proper nodes, >= 0 for ghost nodes and -1 else
! elems is my_id for proper elems, -2 for border elems and -1 else

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
!	flag elements - my_id if internal, -2 if border
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
	  end if
	end do

!	-----------------------------------------------------
!	flag nodes - my_id if internal, area if border
!	-----------------------------------------------------

	do k=1,nkn
	  if( n_area(k) == my_id ) then
	    nodes(k) = my_id
	  else if( nodes(k) == -2 ) then
	    nodes(k) = n_area(k)
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

	subroutine make_domain_final(area_node,nindex,eindex)

! sets id_node, id_elem
!
! only for check - not used anymore
!
! id_node: either my_id for inner node or other for ghost node
! id_elem:
!    (-1,-1) for inner
!    (id,-1) or (id,id) for border to color id (1 or two nodes)
!    (id1,id2) for border to two colors id1 and id2

	use basin
	use shympi

	implicit none

	integer area_node(nkn)
	integer nindex(nkn_local)
	integer eindex(nel_local)

	integer ie,k,ii,n,i
	integer iaux

	!id_node = -1
	!id_elem = -1

!	-----------------------------------------------------
!	sets id_node
!	-----------------------------------------------------

	do i=1,nkn_local
	  k = nindex(i)
	  iaux = area_node(k)
	  if( iaux /= id_node(i) ) then
	    write(6,*) 'difference...: ',i,k,iaux,id_node(i)
	    stop 'error stop make_domain_final'
	  end if
	end do

!	-----------------------------------------------------
!	sets id_elem
!	-----------------------------------------------------

	do ie=1,nel
	  n = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( id_node(k) == my_id ) n = n + 1
	  end do
	  if( n == 3 ) then		!internal elem
	    !nothing - leave at -1
	    !id_elem(1,ie) = my_id
	    if( id_elem(0,ie) /= my_id ) goto 99
	    if( any(id_elem(1:2,ie)/=-1) ) goto 99
	  else if( n > 0 ) then		!border elem
	    n = 0
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( id_node(k) /= my_id ) then
	        n = n + 1
		if( id_elem(n,ie) /= id_node(k) ) goto 99
		id_elem(n,ie) = id_node(k)
	      end if
	    end do
	    !if( id_elem(1,ie) == id_elem(2,ie) ) id_elem(2,ie) = -1
	  else				!error
	    write(6,*) 'writing error message...'
	    write(6,*) n,ie,my_id
	    do ii=1,3
	      k = nen3v(ii,ie)
	      write(6,*) ii,k,id_node(k)
	    end do
	    stop 'error stop make_domain_final: internal error (1)'
	  end if
	end do

!	-----------------------------------------------------
!	end of routine
!	-----------------------------------------------------

	return
   99	continue
	write(6,*) my_id,ie,n,k,id_node(k)
	write(6,*) id_elem(:,ie)
	stop 'error stop make_domain_final: elems...'
	end

!*****************************************************************

	subroutine adjust_indices(n_lk,n_le
     +				,nodes,elems,nindex,eindex)

! computes nkn_local/unique/inner and nel_local/unique/inner
! also rearranges eindex to keep track of this

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
	integer id_aux(0:2,nel)		!total elements
	logical bthis,bexchange

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
	  n = elems(ie)
	  if( n /= my_id ) exit
	  id_aux(0,ie) = n
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
	    id_aux(0,ie) = my_id
	    is = 6 - is
	    k = nen3v(is,ie)
	    n = nodes(k)
	    id_aux(1,ie) = n
	  else if( it == 1 ) then		!one node only with my_id
	    is1 = mod(is,3) + 1
	    is2 = mod(is1,3) + 1
	    k1 = nen3v(is1,ie)
	    k2 = nen3v(is2,ie)
	    n1 = nodes(k1)
	    n2 = nodes(k2)
	    if( n1 /= n2 ) then			!all three nodes are different
	      kmin = minval(nen3v(:,ie))
	      k = nen3v(is,ie)
	      if( kmin == k ) then
	        id_aux(0,ie) = my_id
	      else if( kmin == k1 ) then
	        id_aux(0,ie) = n1
	      else
	        id_aux(0,ie) = n2
	      end if
	      id_aux(1,ie) = n1
	      id_aux(2,ie) = n2
	    else				!two nodes of other color
	      n = nodes(k1)
	      id_aux(0,ie) = n
	      id_aux(1:2,ie) = n
	    end if
	  else
	    stop 'error stop adjust_indices: internal error (1)'
	  end if
	  if( id_aux(0,ie) == my_id ) then
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

	write(my_unit,*) 'debug id_elem: ',my_id
	do i=1,nel_local
	  write(my_unit,*) i,eindex(i),id_elem(:,i)
	end do
	!write(6,*) 'eindex after: ',my_id,eindex(nel_inner+1:nel_local)

	if( nel_unique + id /= nel_local ) then
	  stop 'error stop adjust_indices: internal error (3)'
	end if

!	-----------------------------------------------------
!	end of routine
!	-----------------------------------------------------

	end

!*****************************************************************

	subroutine make_index(my_id,n_g,n_l,items,index)

! returns index of items, sorted first by proper and then ghost items

	implicit none

	integer my_id
	integer n_g		!global number of items
	integer n_l		!local number of items (return)
	integer items(n_g)	!color of items
	integer index(n_g)	!index of local nodes (return)

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

