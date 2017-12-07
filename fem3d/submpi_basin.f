
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

	if( .not. bmpi ) return

	if( shympi_partition_on_elements() ) then
	  stop 'error stop shympi_setup: cannot yet partition on elems'
	end if

	write(6,*) 'setting up mpi with number of threads: ',n_threads

!	=====================================================================
!	the next call sets array node_area(), with values from 0 to n_threads-1
!	=====================================================================

	call make_custom_domain_area

!	=====================================================================
!	from here on everything is general
!	=====================================================================

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

	if( my_unit > 0 ) then
	  write(my_unit,*) 'nodes in domain: ',nkn_local
	  write(my_unit,'(10i7)') (nindex(i),i=1,nkn_local)
	  write(my_unit,*) 'elems in domain: ',nel_local
	  write(my_unit,'(10i7)') (eindex(i),i=1,nel_local)
	end if

	call transfer_domain(nkn_local,nel_local,nindex,eindex)
	call shympi_alloc
	call set_id_node(nkn_local,nindex)
	id_elem = -1
	call make_domain_final

	if( my_unit > 0 ) then
	  write(my_unit,*) 'my_id: ',my_id
	  write(my_unit,*) 'nkn,nel :',nkn,nel
	  write(my_unit,*) 'nkn_inner,nel_inner :',nkn_inner,nel_inner
	  do k=1,nkn
	    write(my_unit,*) k,ipv(k),id_node(k)
	  end do
	  do ie=1,nel
	    write(my_unit,*) ie,ipev(ie),id_elem(:,ie)
	  end do
	end if

	call ghost_make
	call ghost_check
	call ghost_write

	write(6,*) 'mpi my_unit: ',my_unit
	write(6,'(a,9i7)') ' mpi domain: '
     +			,my_id,n_ghost_areas,nkn_global
     +			,nkn_local,nkn_inner,nkn_local-nkn_inner
     +			,nel_local,nel_inner,nel_local-nel_inner

	call shympi_syncronize

	call shympi_alloc_buffer(n_ghost_max)
	call ghost_exchange

        call shympi_univocal_nodes

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

	subroutine set_id_node(n_lk,nindex)

	use basin
	use shympi

	implicit none

	integer n_lk
	integer nindex(n_lk)

	integer i,k

	do i=1,n_lk
	  k = nindex(i)
	  id_node(i) = node_area(k)
	end do

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
	  !id_node(i) = node_area(k)
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

