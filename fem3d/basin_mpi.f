
!*****************************************************************

	subroutine shympi_setup

	use basin
	use shympi

	implicit none

	integer k,ie,i,nc
	integer n_lk,n_le
	integer nodes(nkn)
	integer elems(nel)
	integer nindex(nkn)
	integer eindex(nel)

	call make_domain_area

	call make_domain(my_id,node_area,nodes,elems,nc)

	if( nc+1 /= n_threads ) then
	  write(6,*) 'number of threads = ',n_threads
	  write(6,*) 'number of domains = ',nc+1
	  stop 'error stop: thread/domain mismatch'
	end if

	call make_index(my_id,nkn,n_lk,nodes,nindex)
	call make_index(my_id,nel,n_le,elems,eindex)

	!write(6,*) 'nodes in domain: ',n_lk
	!write(6,*) 'elems in domain: ',n_le

	nkn_local = n_lk
	nel_local = n_le

	write(my_unit,*) 'nodes in domain: ',nkn_local
	write(my_unit,'(10i7)') (nindex(i),i=1,nkn_local)
	write(my_unit,*) 'elems in domain: ',nel_local
	write(my_unit,'(10i7)') (eindex(i),i=1,nel_local)

	call shympi_alloc
	call transfer_domain(nkn_local,nel_local,nindex,eindex)
	id_elem = -1
	call make_domain_final

	write(my_unit,*) 'my_id: ',my_id
	write(my_unit,*) 'nkn,nel :',nkn,nel
	write(my_unit,*) 'nkn_inner,nel_inner :',nkn_inner,nel_inner
	do k=1,nkn
	  write(my_unit,*) k,ipv(k),id_node(k)
	end do
	do ie=1,nel
	  write(my_unit,*) ie,ipev(ie),id_elem(:,ie)
	end do

	call ghost_make
	call ghost_check
	call ghost_write

	write(6,'(a,9i7)') 'domain: ',my_id,n_ghost_areas,nkn_global
     +			,nkn_local,nkn_inner,nkn_local-nkn_inner
     +			,nel_local,nel_inner,nel_local-nel_inner

	call shympi_syncronize

	call shympi_alloc_buffer(n_ghost_max)
	call ghost_exchange

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine make_domain_area

	use basin
	use shympi

	implicit none

	integer ie,ii
	real r,h

	node_area = 0

	if( nkn /= 225 .or. nel /= 384 ) then
	  if( n_threads > 1 ) then
	    write(6,*) 'nkn,nel: ',nkn,nel
	    write(6,*) 'expecting: ',225,384
	    stop 'error stop make_domain_area: wrong basin'
	  end if
	end if

	if( n_threads == 1 ) then
	  return
	else if( n_threads == 2 ) then
	  call make_domain_area_2
	else if( n_threads == 3 ) then
	  call make_domain_area_3
	else if( n_threads == 4 ) then
	  call make_domain_area_4
	else
	  write(6,*) 'n_threads = ',n_threads
	  stop 'error stop make_domain_area: cannot handle'
	end if

	do ie=1,nel
	  do ii=1,3
	    call random_number(r)
	    h = 10.* r
	    h = 2.* r
	    hm3v(ii,ie) = h
	    !write(6,*) ie,ii,r,h
	  end do
	end do

	end

!*****************************************************************

	subroutine make_domain_area_2

	use basin
	use shympi

	implicit none

	integer k

	do k=1,nkn
	  if( ygv(k) > 3100.  ) then
	    node_area(k) = 1
	  end if
	end do

	end

!*****************************************************************

	subroutine make_domain_area_3

	use basin
	use shympi

	implicit none

	integer k

	do k=1,nkn
	  if( ygv(k) > 4100.  ) then
	    node_area(k) = 2
	  else if( ygv(k) > 2100.  ) then
	    node_area(k) = 1
	  end if
	end do

	end

!*****************************************************************

	subroutine make_domain_area_4

	use basin
	use shympi

	implicit none

	integer k

	do k=1,nkn
	  if( ygv(k) > 3100.  ) then
	    node_area(k) = 2
	  end if
	  if( xgv(k) > 100.  ) then
	    node_area(k) = node_area(k) + 1
	  end if
	end do

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine make_domain(my_id,node_area,nodes,elems,nc)

	use basin

	implicit none

	integer my_id
	integer node_area(nkn)
	integer nodes(nkn)
	integer elems(nel)

	integer ie,k,ii,nc
	integer n,n_my,n_ghost
	integer, allocatable :: ncs(:)

	nodes = -1
	elems = -1

	nc = maxval(node_area)
	allocate(ncs(0:nc))
	ncs = 0

	if( my_id < 0 .or. my_id > nc ) then
	  write(6,*) 'my_id,nc : ',my_id,nc
	  stop 'error stop make_domain: my_id out of range'
	end if

	do k=1,nkn
	  n = node_area(k)
	  if( n < 0 .or. n > nc ) then
	    write(6,*) 'node_area = ',n,' in node ',k
	    stop 'error stop make_domain: out of range'
	  end if
	  ncs(n) = ncs(n) + 1
	end do

	do n=0,nc
	  if( ncs(n) == 0 ) then
	    write(6,*) 'node_area = ',n,' has 0 nodes'
	    stop 'error stop make_domain:  no nodes'
	  end if
	end do

	do ie=1,nel
	  n = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( node_area(k) == my_id ) n = n + 1
	  end do
	  if( n == 3 ) then		!internal elem
	    elems(ie) = my_id
	  else if( n > 0 ) then		!border elem
	    elems(ie) = -2
	    do ii=1,3
	      k = nen3v(ii,ie)
	      nodes(k) = -2
	    end do
	  end if
	end do

	do k=1,nkn
	  if( node_area(k) == my_id ) then
	    nodes(k) = my_id
	  else if( nodes(k) == -2 ) then
	    nodes(k) = node_area(k)
	  end if
	end do

	n_my = 0
	n_ghost = 0
	do k=1,nkn
	  if( nodes(k) == my_id ) then
	    n_my = n_my + 1
	  else if( nodes(k) >= 0 ) then
	    n_ghost = n_ghost + 1
	  else if( nodes(k) == -1 ) then
	    !nothing
	  else
	    write(6,*) nodes(k)
	    stop 'error stop make_domain:  impossible value for nodes'
	  end if
	end do

	n = ncs(my_id)
	if( n /= n_my ) then
	  write(6,*) n,n_my,n_ghost,n_my+n_ghost
	  stop 'error stop make_domain:  internal error (2)'
	end if

	!write(6,*) 'domain = ',my_id,n_my,n_ghost

	deallocate(ncs)

	end

!*****************************************************************

	subroutine make_domain_final

	use basin
	use shympi

	implicit none

	integer ie,k,ii,n

	id_elem = -1

	do ie=1,nel
	  n = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( id_node(k) == my_id ) n = n + 1
	  end do
	  if( n == 3 ) then		!internal elem
	    id_elem(1,ie) = my_id
	  else if( n > 0 ) then		!border elem
	    n = 0
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( id_node(k) /= my_id ) then
	        n = n + 1
		id_elem(n,ie) = id_node(k)
	      end if
	    end do
	    if( id_elem(1,ie) == id_elem(2,ie) ) id_elem(2,ie) = -1
	  else
	    write(6,*) 'writing error message...'
	    write(6,*) n,ie,my_id
	    do ii=1,3
	      k = nen3v(ii,ie)
	      write(6,*) ii,k,id_node(k)
	    end do
	    stop 'error stop make_domain_final: internal error (1)'
	  end if
	end do

	do k=1,nkn
	  is_inner_node(k) = (id_node(k) == my_id)
	  if( is_inner_node(k) ) nkn_inner = k
	end do

	do ie=1,nel
	  is_inner_elem(ie) = (id_elem(1,ie) == my_id)
	  if( is_inner_elem(ie) ) nel_inner = ie
	end do

	end

!*****************************************************************

	subroutine make_index(my_id,n_g,n_l,nodes,index)

	implicit none

	integer my_id
	integer n_g,n_l
	integer nodes(n_g)
	integer index(n_g)

	integer i,k,n

	i = 0
	do k=1,n_g
	  n = nodes(k)
	  if( n == my_id ) then
	    i = i + 1
	    index(i) = k
	  end if
	end do

	do k=1,n_g
	  n = nodes(k)
	  if( n /= my_id .and. n /= -1 ) then
	    i = i + 1
	    index(i) = k
	  end if
	end do

	n_l = i
	if( i /= n_l ) then
	  write(6,*) i,n_l
	  stop 'error stop transfer_internal:  internal error (1)'
	end if

	end

!*****************************************************************

	subroutine transfer_domain(n_lk,n_le,nindex,eindex)

! transfers global domain to local domain

	use basin
	use shympi

	implicit none

	integer n_lk,n_le
	integer nindex(n_lk)
	integer eindex(n_le)

	integer ie,k,i,ii,kk
	integer inverse(nkn)

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
!	transfer information
!	----------------------------------

	inverse = 0
	do i=1,n_lk
	  k = nindex(i)
	  inverse(k) = i
	end do

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
	  id_node(i) = node_area(k)
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
!	end routine
!	----------------------------------

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ghost_make

! makes list of ghost nodes

	use basin
	use shympi

	implicit none

	integer nc,k,id,i,n,ncsmax,ia,ic,ie,ii,id1,id2
	integer iaux(nkn)
	integer, allocatable :: ncs(:)
	integer, allocatable :: ga(:)

!	--------------------------------------------------
!	find total number of neibor ghost areas
!	--------------------------------------------------

	nc = n_threads - 1
	allocate(ncs(0:nc))
	ncs = 0

	do k=1,nkn
	  id = id_node(k)
	  if( id /= my_id ) ncs(id) = ncs(id) + 1
	end do

	n_ghost_areas = 0
	do n=0,nc
	  if( ncs(n) > 0 ) n_ghost_areas = n_ghost_areas + 1
	end do

!	--------------------------------------------------
!	collect info on neibor ghost areas
!	--------------------------------------------------

	allocate(ga(n_ghost_areas))	!temporary
	ga = 0

	i = 0
	ncsmax = 0
	do n=0,nc
	  if( ncs(n) > 0 ) then
	    i = i + 1
	    ncsmax = max(ncsmax,ncs(n))
	    ga(i) = n	!what id
	  end if
	end do

	n_ghost_nodes_max = ncsmax	!outer ghost nodes

!	--------------------------------------------------
!	find maximum of inner ghost nodes
!	--------------------------------------------------

	ncsmax = 0
	do ia=1,n_ghost_areas
	  ic = ga(ia)
	  iaux = 0
	  do ie=1,nel
	    if( is_inner_elem(ie) ) cycle
	    if( id_elem(1,ie) /= ic .and. id_elem(2,ie) /= ic ) cycle
	    do ii=1,3
	      k = nen3v(ii,ie)
	      iaux(k) = iaux(k) + 1
	    end do
	  end do
	  nc = 0
	  do k=1,nkn
	    if( is_inner_node(k) ) cycle
	    if( iaux(k) == 0 ) cycle
	    nc = nc + 1
	  end do
	  ncsmax = max(ncsmax,nc)
	end do

	n_ghost_nodes_max = max(n_ghost_nodes_max,ncsmax)

!	--------------------------------------------------
!	find maximum of ghost elements
!	--------------------------------------------------

	ncs = 0
	do ie=1,nel
	  do i=1,2
	    id = id_elem(i,ie)
	    if( id /= my_id .and. id /= -1 ) then
	      ncs(id) = ncs(id) + 1
	    end if
	  end do
	end do

	ncsmax = 0
	do ia=1,n_ghost_areas
	  ic = ga(ia)
	  ncsmax = max(ncsmax,ncs(ic))
	end do
	n_ghost_elems_max = ncsmax

!	--------------------------------------------------
!	allocate ghost arrays
!	--------------------------------------------------

	n_ghost_max = max(n_ghost_nodes_max,n_ghost_elems_max)
	call shympi_alloc_ghost(n_ghost_max)
	ghost_areas(1,:) = ga(:)
	deallocate(ga)

!	--------------------------------------------------
!	set up list of outer ghost nodes
!	--------------------------------------------------

	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  nc = 0
	  do k=1,nkn
	    id = id_node(k)
	    if( id /= ic ) cycle
	    nc = nc + 1
	    if( nc > ncsmax ) then
	      write(6,*) ia,id,nc,ncsmax
	      stop 'error stop ghost_make: internal error (1)'
	    end if
	    ghost_nodes_out(nc,ia) = k
	  end do
	  ghost_areas(2,ia) = nc
	end do

!	--------------------------------------------------
!	set up list of inner ghost nodes
!	--------------------------------------------------

	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  iaux = 0
	  do ie=1,nel
	    if( is_inner_elem(ie) ) cycle
	    if( id_elem(1,ie) /= ic .and. id_elem(2,ie) /= ic ) cycle
	    do ii=1,3
	      k = nen3v(ii,ie)
	      iaux(k) = iaux(k) + 1
	    end do
	  end do
	  nc = 0
	  do k=1,nkn
	    if( .not. is_inner_node(k) ) cycle
	    if( iaux(k) == 0 ) cycle
	    nc = nc + 1
	    if( nc > ncsmax ) then
	      write(6,*) ia,ic,nc,ncsmax
	      stop 'error stop ghost_make: internal error (2)'
	    end if
	    ghost_nodes_in(nc,ia) = k
	  end do
	  ghost_areas(3,ia) = nc
	end do

!	--------------------------------------------------
!	set up list of ghost elements
!	--------------------------------------------------

	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  nc = 0
	  do ie=1,nel
	    id1 = id_elem(1,ie)
	    id2 = id_elem(2,ie)
	    if( id1 /= ic .and. id2 /= ic ) cycle
	    nc = nc + 1
	    if( nc > ncsmax ) then
	      write(6,*) ia,id,nc,ncsmax
	      stop 'error stop ghost_make: internal error (3)'
	    end if
	    ghost_elems(nc,ia) = ie
	  end do
	  ghost_areas(4,ia) = nc
	end do

!	--------------------------------------------------
!	end of routine
!	--------------------------------------------------

	deallocate(ncs)

	end

!*****************************************************************

	subroutine ghost_check

	use shympi
	use basin

	implicit none

	integer ia,i,ic,k,nc

	do ia=1,n_ghost_areas

	  ic = ghost_areas(1,ia)
	  if( my_id == ic ) then
	    write(6,*) ia,ic,my_id
	    stop 'error stop ghost_check: internal error (1)'
	  end if

	  nc = ghost_areas(2,ia)
	  do i=1,nc
	    k = ghost_nodes_out(i,ia)
	    if( id_node(k) /= ic ) then
	      write(6,*) ia,i,k,id_node(k),ic
	      stop 'error stop ghost_check: internal error (2)'
	    end if
	  end do

	  nc = ghost_areas(3,ia)
	  do i=1,nc
	    k = ghost_nodes_in(i,ia)
	    if( id_node(k) /= my_id ) then
	      write(6,*) ia,i,k,id_node(k),my_id
	      stop 'error stop ghost_check: internal error (3)'
	    end if
	  end do

	end do

	end

!*****************************************************************

	subroutine ghost_write

	use shympi
	use basin

	implicit none

	integer ia,ic,nc,i,k,ie

	write(my_unit,*) 'n_ghost_areas = ',n_ghost_areas
	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  nc = ghost_areas(2,ia)
	  write(my_unit,*) 'outer: ',ic,nc
	  do i=1,nc
	    k = ghost_nodes_out(i,ia)
	    write(my_unit,*) k,ipv(k),xgv(k),ygv(k)
	  end do
	  nc = ghost_areas(3,ia)
	  write(my_unit,*) 'inner: ',ic,nc
	  do i=1,nc
	    k = ghost_nodes_in(i,ia)
	    write(my_unit,*) k,ipv(k),xgv(k),ygv(k)
	  end do
	  nc = ghost_areas(4,ia)
	  write(my_unit,*) 'elems: ',ic,nc
	  do i=1,nc
	    ie = ghost_elems(i,ia)
	    write(my_unit,*) ie,ipev(ie)
	  end do
	end do

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine print_ghost_nodes_r(val,text)

	use shympi
	use basin

	implicit none

	real val(nkn)
	character*(*) text

	integer ia,ic,nc,k,i

	write(my_unit,*) 'printing ghost nodes: ' // text
	write(my_unit,*) 'n_ghost_areas = ',n_ghost_areas,my_id

	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  nc = ghost_areas(2,ia)
	  write(my_unit,*) 'outer: ',ic,nc
	  do i=1,nc
	    k = ghost_nodes_out(i,ia)
	    write(my_unit,*) k,ipv(k),val(k)
	  end do
	  nc = ghost_areas(3,ia)
	  write(my_unit,*) 'inner: ',ic,nc
	  do i=1,nc
	    k = ghost_nodes_in(i,ia)
	    write(my_unit,*) k,ipv(k),val(k)
	  end do
	end do

	end

!*****************************************************************

	subroutine print_ghost_nodes_i(val,text)

	use shympi
	use basin

	implicit none

	integer val(nkn)
	character*(*) text

	integer ia,ic,nc,k,i

	write(my_unit,*) 'printing ghost nodes: ' // text
	write(my_unit,*) 'n_ghost_areas = ',n_ghost_areas,my_id

	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  nc = ghost_areas(2,ia)
	  write(my_unit,*) 'outer: ',ic,nc
	  do i=1,nc
	    k = ghost_nodes_out(i,ia)
	    write(my_unit,*) k,ipv(k),val(k)
	  end do
	  nc = ghost_areas(3,ia)
	  write(my_unit,*) 'inner: ',ic,nc
	  do i=1,nc
	    k = ghost_nodes_in(i,ia)
	    write(my_unit,*) k,ipv(k),val(k)
	  end do
	end do

	end

!*****************************************************************

	subroutine print_ghost_elems_i(val,text)

	use shympi
	use basin

	implicit none

	integer val(nel)
	character*(*) text

	integer ia,ic,nc,k,i,ie

	write(my_unit,*) 'printing ghost elems: ' // text
	write(my_unit,*) 'n_ghost_areas = ',n_ghost_areas,my_id

	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  nc = ghost_areas(4,ia)
	  write(my_unit,*) 'elems: ',ic,nc
	  do i=1,nc
	    ie = ghost_elems(i,ia)
	    write(my_unit,*) ie,ipev(ie),val(ie)
	  end do
	end do

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ghost_exchange

	use shympi
	use basin

	implicit none

	integer i
	integer num_elems(nel)
	integer num_nodes(nkn)

	num_elems = ipev
	num_nodes = ipv

	!call shympi_exchange_2d_elem_i(ipev)
	call shympi_exchange_2d_node_i(ipv)

	call shympi_check_2d_elem_i(ipev,'ghost ipev')
	call shympi_check_2d_node_i(ipv,'ghost ipv')

	call shympi_check_array_i(nel,num_elems,ipev,'ghost ipev')
	call shympi_check_array_i(nkn,num_nodes,ipv,'ghost ipv')

	write(6,*) 'finished exchange...'

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

