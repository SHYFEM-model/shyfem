
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
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 10.04.2022	ggu	adjourned call to shympi_check_array()
! 10.03.2023    ggu     new routine ghost_debug()
! 18.03.2023	ggu	resolved problem exchanging elements (ghost nodes)
! 20.03.2023	ggu	new subroutine ghost_handle()
! 24.03.2023	ggu	bug fix... ic not defined
! 27.03.2023	ggu	bug fix... iloop == 4 eliminated
! 27.03.2023	ggu	bug fix... no ieaux
! 19.04.2023	ggu	in ghost_exchange adapt call to shympi_check_array()

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine ghost_handle

! handles ghost nodes (general routine

	use shympi

	implicit none

        call ghost_make         !here also call to shympi_alloc_ghost()
        !call ghost_debug
        call ghost_write
        call ghost_check

!       -----------------------------------------------------
!       debug output for ghost nodes
!       -----------------------------------------------------

        write(6,*) 'mpi my_unit: ',my_unit
        write(6,'(a,9i7)') ' mpi domain: ',my_id,n_ghost_areas
        write(6,'(a,5i7)') 'nkn: '
     +                  ,nkn_global,nkn_local,nkn_unique
     +                  ,nkn_inner,nkn_local-nkn_inner
        write(6,'(a,5i7)') 'nel: '
     +                  ,nel_global,nel_local,nel_unique
     +                  ,nel_inner,nel_local-nel_inner

        call shympi_syncronize

        call shympi_alloc_buffer(n_ghost_max)		!really not needed

        call ghost_exchange

	end

!*****************************************************************

	subroutine ghost_make

! makes list of ghost nodes - basin has been already transfered

	use basin
	use shympi

	implicit none

	integer k,id,i,n,ncsmax,ia,ic,ie,ii,iu,iu5,iu6
	integer nc,nc_in,nc_out
	integer iea,ies,iloop,id0
	integer iext,kext
	integer iaux(nkn)
	integer, allocatable :: ncs(:)
	integer, allocatable :: ga(:)

	integer ipext,ieext

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

	write(my_unit,*) 'debug ghost: ',my_id
	write(my_unit,*) n_ghost_areas
	write(my_unit,*) ncs

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

	write(my_unit,*) ga
	write(my_unit,*) 'outer_max: ',ncsmax

	n_ghost_nodes_max = ncsmax	!outer ghost nodes

!	--------------------------------------------------
!	find maximum of inner ghost nodes
!	--------------------------------------------------

	ncsmax = 0
	do ia=1,n_ghost_areas
	  ic = ga(ia)			!color of ghost area
	  iaux = 0
	  do ie=1,nel
	    if( shympi_is_inner_elem(ie) ) cycle
	    if( id_elem(2,ie) /= ic .and. id_elem(3,ie) /= ic ) cycle
	    !if( all(id_elem(1:3,ie) /= ic ) ) cycle	!GGUghost
	    do ii=1,3
	      k = nen3v(ii,ie)
	      iaux(k) = iaux(k) + 1
	    end do
	  end do
	  nc = 0
	  do k=1,nkn
	    if( shympi_is_inner_node(k) ) cycle
	    if( iaux(k) == 0 ) cycle
	    nc = nc + 1
	  end do
	  ncsmax = max(ncsmax,nc)
	end do

	write(my_unit,*) 'inner_max: ',ncsmax

	n_ghost_nodes_max = max(n_ghost_nodes_max,ncsmax)

!	--------------------------------------------------
!	find maximum of ghost elements
!	--------------------------------------------------

	ncs = 0
	do ie=1,nel
	  n = id_elem(0,ie)
	  do i=1,n
	    id = id_elem(i,ie)
	    if( id == my_id ) cycle
	    if( id == -1 ) cycle
	    ncs(id) = ncs(id) + 1
	  end do
	end do

	ncsmax = 0
	do ia=1,n_ghost_areas
	  ic = ga(ia)
	  ncsmax = max(ncsmax,ncs(ic))
	end do
	n_ghost_elems_max = ncsmax
	write(my_unit,*) 'maximum elem: ',ncsmax

!	--------------------------------------------------
!	allocate ghost arrays
!	--------------------------------------------------

	n_ghost_max = max(n_ghost_nodes_max,n_ghost_elems_max)
	write(my_unit,*) 'n_ghost_max: ',n_ghost_max
	call shympi_alloc_ghost(n_ghost_max)
	ghost_areas(1,:) = ga(:)
	deallocate(ga)
	ncsmax = n_ghost_max

	iu5 = 500 + my_id
	!write(iu5,*) 'ncsmax: ',ncsmax,nkn,nel

!	--------------------------------------------------
!	set up list of outer ghost nodes
!	--------------------------------------------------

	iloop = 1

	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  nc = 0
	  do k=1,nkn
	    id = id_node(k)
	    if( id /= ic ) cycle
	    if( id == my_id ) cycle
	    nc = nc + 1
	    if( nc > ncsmax ) goto 99
	    ghost_nodes_out(nc,ia) = k
	  end do
	  ghost_areas(2,ia) = nc
	  if( nc > ncsmax ) goto 99	!just to be sure...
	end do

!	--------------------------------------------------
!	set up list of inner ghost nodes
!	--------------------------------------------------

	iloop = 2

	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  iaux = 0
	  do ie=1,nel
	    if( shympi_is_inner_elem(ie) ) cycle
	    !if( id_elem(2,ie) /= ic .and. id_elem(3,ie) /= ic ) cycle
	    if( any(id_elem(1:3,ie) == ic ) ) then
	    do ii=1,3
	      k = nen3v(ii,ie)
	      iaux(k) = iaux(k) + 1
	    end do
	    end if
	  end do
	  nc = 0
	  do k=1,nkn
	    if( .not. shympi_is_inner_node(k) ) cycle
	    if( iaux(k) == 0 ) cycle
	    nc = nc + 1
	    if( nc > ncsmax ) goto 99
	    ghost_nodes_in(nc,ia) = k
	  end do
	  ghost_areas(3,ia) = nc
	end do

!	--------------------------------------------------
!	set up list of ghost elements
!	--------------------------------------------------

	iloop = 3

	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  nc_in = 0
	  nc_out = 0
	  do ie=1,nel
	    if( shympi_is_inner_elem(ie) ) cycle
	    id0 = id_elem(1,ie)
	    if( .not. any(id_elem(1:3,ie)==ic) ) cycle
	    if( id0 == my_id ) then
	      nc_in = nc_in + 1
	      if( nc_in > ncsmax ) goto 99
	      ghost_elems_in(nc_in,ia) = ie
	    else
	      if( id0 /= ic ) cycle	!no inner element for this outer
	      nc_out = nc_out + 1
	      if( nc_out > ncsmax ) goto 99
	      ghost_elems_out(nc_out,ia) = ie
	    end if
	    !write(iu5,'(7i8)') nc_in,nc_out,ie,id_elem(:,ie)
	  end do
	  ghost_areas(4,ia) = nc_out
	  ghost_areas(5,ia) = nc_in
	end do

!	--------------------------------------------------
!	write debug information
!	--------------------------------------------------

	if( bmpi_debug ) then

	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  nc = ghost_areas(4,ia)
	  if( nc > ncsmax ) goto 99
	  write(my_unit,*) 'node test: ',my_id,nkn
	  do k=1,nkn
	    write(my_unit,*) k,ipv(k),id_node(k)
	  end do
	  write(my_unit,*) 'elem test: ',my_id,nel
	  do ie=1,nel
	    write(my_unit,*) ie,ipev(ie),id_elem(:,ie)
	  end do
	  write(my_unit,*) 'ghost test: ',my_id,ia,ic
	  nc = ghost_areas(5,ia)
	  write(my_unit,*) '   inner elems...',nc
	  do i=1,nc
	    ie = ghost_elems_in(i,ia)
	    if( ie < 1 ) goto 97
	    if( ie > nel ) goto 97
	    write(my_unit,*) i,ie,ipev(ie),id_elem(1,ie)
	  end do
	  nc = ghost_areas(4,ia)
	  write(my_unit,*) '   outer elems...',nc
	  do i=1,nc
	    ie = ghost_elems_out(i,ia)
	    if( ie < 1 ) goto 97
	    if( ie > nel ) goto 97
	    write(my_unit,*) i,ie,ipev(ie),id_elem(1,ie)
	  end do
	end do

	iu = 300 + my_id
	write(iu,'(a,6i10)') 'looking for tripple points',my_id
	do ie=1,nel
	  if( id_elem(0,ie) == 3 ) then
	    iext = ieext(ie)
	    write(iu,'(a,6i10)') 'tripple: ',ie,iext,id_elem(:,ie)
	  end if
	end do
	flush(iu)

	iu6 = 600 + my_id
	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  nc = ghost_areas(2,ia)
	  write(iu6,*) '----------------------'
	  write(iu6,*) 'color',my_id,ic
	  write(iu6,*) '----------------------'
	  write(iu6,*) 'outer nodes',my_id,ic,nc
	  do i=1,nc
	    k = ghost_nodes_out(i,ia)
	    if( id_node(k) == my_id ) stop 'error stop 600 2'
	    kext = ipext(k)
	    write(iu6,*) i,k,kext
	  end do
	  nc = ghost_areas(3,ia)
	  write(iu6,*) 'inner nodes',my_id,ic,nc
	  do i=1,nc
	    k = ghost_nodes_in(i,ia)
	    if( id_node(k) /= my_id ) stop 'error stop 600 3'
	    kext = ipext(k)
	    write(iu6,*) i,k,kext
	  end do
	  nc = ghost_areas(4,ia)
	  write(iu6,*) 'outer elems',my_id,ic,nc
	  do i=1,nc
	    ie = ghost_elems_out(i,ia)
	    if( id_elem(1,ie) == my_id ) stop 'error stop 600 4'
	    iext = ieext(ie)
	    write(iu6,*) i,ie,iext
	  end do
	  nc = ghost_areas(5,ia)
	  write(iu6,*) 'inner elems',my_id,ic,nc
	  do i=1,nc
	    ie = ghost_elems_in(i,ia)
	    if( id_elem(1,ie) /= my_id ) stop 'error stop 600 5'
	    iext = ieext(ie)
	    write(iu6,*) i,ie,iext
	  end do
	end do
	flush(iu6)

	end if

!	--------------------------------------------------
!	end of routine
!	--------------------------------------------------

	deallocate(ncs)

	write(my_unit,*) 'finished setting up ghost: ',my_id
	flush(my_unit)

	return
   97	continue
	write(6,*) 'no elements in list',ia,nc,ie
	write(6,*) ghost_elems_in(1:nc,ia)
	write(6,*) ghost_elems_out(1:nc,ia)
	flush(iu)
	flush(iu5)
	stop 'error stop ghost_make: internal error (8)'
   99	continue
	write(6,*) 'iloop = ',iloop
	write(6,*) 'my_id = ',my_id
	write(6,*) ia,id,nc,ncsmax
	flush(iu)
	flush(iu5)
	stop 'error stop ghost_make: internal error (7)'
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

	write(my_unit,*) 'writing ghost node info for my_id: ',my_id

	write(my_unit,*) 'n_ghost_areas = ',n_ghost_areas
	write(my_unit,*) 'ghost_areas: ',ghost_areas

	do ia=1,n_ghost_areas
	  write(my_unit,*) 'handling ghost_area: ',ia
	  write(my_unit,*) 'ghost_area: ',ghost_areas(:,ia)
	  ic = ghost_areas(1,ia)
	  nc = ghost_areas(2,ia)
	  write(my_unit,*) 'ghost_area color: ',ic
	  write(my_unit,*) 'nkn inner,unique,local: '
     +     				,nkn_inner,nkn_unique,nkn_local
	  write(my_unit,*) 'nodes outer: ',ic,nc
	  flush(my_unit)
	  do i=1,nc
	    write(my_unit,*) i,ia
	  flush(my_unit)
	    k = ghost_nodes_out(i,ia)
	    write(my_unit,*) k
	  flush(my_unit)
	    write(my_unit,*) k,ipv(k),xgv(k),ygv(k)
	  flush(my_unit)
	  end do
	  nc = ghost_areas(3,ia)
	  write(my_unit,*) 'nodes inner: ',ic,nc
	  flush(my_unit)
	  do i=1,nc
	    k = ghost_nodes_in(i,ia)
	    write(my_unit,*) k,ipv(k),xgv(k),ygv(k)
	  flush(my_unit)
	  end do
	  nc = ghost_areas(4,ia)
	  write(my_unit,*) 'nel inner,unique,local: '
     +     				,nel_inner,nel_unique,nel_local
	  flush(my_unit)
	  write(my_unit,*) 'elems outer: ',ic,nc
	  flush(my_unit)
	  do i=1,nc
	    ie = ghost_elems_out(i,ia)
	    write(my_unit,*) ie,ipev(ie)
	  flush(my_unit)
	  end do
	  nc = ghost_areas(5,ia)
	  write(my_unit,*) 'elems inner: ',ic,nc
	  flush(my_unit)
	  do i=1,nc
	    ie = ghost_elems_in(i,ia)
	    write(my_unit,*) ie,ipev(ie)
	  flush(my_unit)
	  end do
	end do

	write(my_unit,*) 'finished writing ghost node info'
	  flush(my_unit)

	end

!*****************************************************************

	subroutine ghost_debug

	use shympi
	use basin

	implicit none

	integer ia,ic,nc,i,k,ie,iu
	integer kext,iext
	integer ipext,ieext

	iu = 200 + my_id

	write(iu,*) '=================================='
	write(iu,*) 'writing ghost_debug: ',my_id
	write(iu,*) '=================================='

	write(iu,*) 'n_ghost_areas = ',n_ghost_areas
	write(iu,*) 'ghost_areas: ',ghost_areas
	write(iu,*) 'ghost_nodes_out: ',ghost_nodes_out
	write(iu,*) 'ghost_nodes_in: ',ghost_nodes_in
	write(iu,*) 'ghost_elems_out: ',ghost_elems_out
	write(iu,*) 'ghost_elems_in: ',ghost_elems_in

	do ia=1,n_ghost_areas
	  ic = ghost_areas(1,ia)
	  write(iu,*) '-----------------------------'
	  write(iu,*) 'ghost color:',ia,ic
	  write(iu,*) '-----------------------------'
	  nc = ghost_areas(2,ia)
	  write(iu,*) 'outer nodes',nc
	  do i=1,nc
	    k = ghost_nodes_out(i,ia)
	    kext = ipext(k)
	    write(iu,*) i,k,kext
	  end do
	  nc = ghost_areas(3,ia)
	  write(iu,*) 'inner nodes',nc
	  do i=1,nc
	    k = ghost_nodes_in(i,ia)
	    kext = ipext(k)
	    write(iu,*) i,k,kext
	  end do
	  nc = ghost_areas(4,ia)
	  write(iu,*) 'outer elems',nc
	  do i=1,nc
	    ie = ghost_elems_out(i,ia)
	    iext = ieext(ie)
	    write(iu,*) i,ie,iext
	  end do
	  nc = ghost_areas(5,ia)
	  write(iu,*) 'inner elems',nc
	  do i=1,nc
	    ie = ghost_elems_in(i,ia)
	    iext = ieext(ie)
	    write(iu,*) i,ie,iext
	  end do
	end do

	write(iu,*) 'finished writing ghost_debug'
	flush(iu)

	end

!*****************************************************************

        subroutine ghost_debug_1

        use shympi
        use basin

        implicit none

        integer ia,ic,nc,i,k,ie,iu

        iu = 200 + my_id

        write(iu,*) 'writing ghost debug: ',my_id

        write(iu,*) 'n_ghost_areas = ',n_ghost_areas
        write(iu,*) 'ghost_areas: ',ghost_areas
        write(iu,*) 'ghost_nodes_out: ',ghost_nodes_out
        write(iu,*) 'ghost_nodes_in: ',ghost_nodes_in
        write(iu,*) 'ghost_elems_out: ',ghost_elems_out
        write(iu,*) 'ghost_elems_in: ',ghost_elems_in

        write(iu,*) 'finished writing ghost debug'
        flush(iu)

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

	logical, parameter :: be = .true.
	logical, parameter :: bn = .false.
	integer i
	integer, allocatable :: num_elems(:)
	integer, allocatable :: num_nodes(:)
	integer, allocatable :: num_e(:)
	integer, allocatable :: num_n(:)
	
	allocate(num_elems(nel))
	allocate(num_nodes(nkn))
	num_elems = ipev
	num_nodes = ipv

	write(6,*) 'start exchange ghost',my_id
	flush(6)
	call shympi_syncronize

	call shympi_exchange_2d_node(ipv)
	write(6,*) 'exchange ghost nodes',my_id,nkn
	flush(6)
	call shympi_syncronize
	call shympi_exchange_2d_elem(ipev)
	write(6,*) 'exchange ghost elems',my_id,nel
	flush(6)
	call shympi_syncronize

	call shympi_check_2d_node(ipv,'ghost ipv')
	write(6,*) 'check ghost nodes',my_id
	flush(6)
	call shympi_syncronize
	call shympi_check_2d_elem(ipev,'ghost ipev')
	write(6,*) 'check ghost elems',my_id
	flush(6)
	call shympi_syncronize

	allocate(num_e(nel))
	allocate(num_n(nkn))
	num_e = ipev
	num_n = ipv

	call shympi_syncronize
	write(6,*) 'extra checks',my_id
	flush(6)

	i = count( num_nodes /= num_n )
	write(6,*) 'different nodes: ',my_id,i
	i = count( num_elems /= num_e )
	write(6,*) 'different elems: ',my_id,i
	flush(6)
	call shympi_syncronize

	call shympi_check_array(bn,1,nkn,nkn,num_nodes,ipv,'ghost ipv')
	call shympi_check_array(be,1,nel,nel,num_elems,ipev,'ghost ipev')
	call shympi_check_array(bn,1,nkn,nkn,num_n,ipv,'ghost ipv')
	call shympi_check_array(be,1,nel,nel,num_e,ipev,'ghost ipev')

	call shympi_syncronize
	write(6,*) 'ghost_exchange: finished exchange...',my_id
	flush(6)
	!stop

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

