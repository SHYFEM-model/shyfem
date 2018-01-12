
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ghost_make

! makes list of ghost nodes - basin has been already transfered

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
	    ghost_elems_in(nc,ia) = ie
	    ghost_elems_out(nc,ia) = ie
	  end do
	  ghost_areas(4,ia) = nc
	  ghost_areas(5,ia) = nc
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
	  write(my_unit,*) 'nodes outer: ',ic,nc
	  do i=1,nc
	    k = ghost_nodes_out(i,ia)
	    write(my_unit,*) k,ipv(k),xgv(k),ygv(k)
	  end do
	  nc = ghost_areas(3,ia)
	  write(my_unit,*) 'nodes inner: ',ic,nc
	  do i=1,nc
	    k = ghost_nodes_in(i,ia)
	    write(my_unit,*) k,ipv(k),xgv(k),ygv(k)
	  end do
	  nc = ghost_areas(4,ia)
	  write(my_unit,*) 'elems outer: ',ic,nc
	  do i=1,nc
	    ie = ghost_elems_out(i,ia)
	    write(my_unit,*) ie,ipev(ie)
	  end do
	  nc = ghost_areas(5,ia)
	  write(my_unit,*) 'elems inner: ',ic,nc
	  do i=1,nc
	    ie = ghost_elems_in(i,ia)
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
	  write(my_unit,*) 'elems outer: ',ic,nc
	  do i=1,nc
	    ie = ghost_elems_out(i,ia)
	    write(my_unit,*) ie,ipev(ie),val(ie)
	  end do
	  nc = ghost_areas(5,ia)
	  write(my_unit,*) 'elems inner: ',ic,nc
	  do i=1,nc
	    ie = ghost_elems_in(i,ia)
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

	call shympi_exchange_2d_node(ipv)
	call shympi_exchange_2d_elem(ipev)

	call shympi_check_2d_node(ipv,'ghost ipv')
	call shympi_check_2d_elem(ipev,'ghost ipev')

	call shympi_check_array(nkn,num_nodes,ipv,'ghost ipv')
	call shympi_check_array(nel,num_elems,ipev,'ghost ipev')

	write(6,*) 'finished exchange...'

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

