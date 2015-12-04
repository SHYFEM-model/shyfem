!
! mpi routines
!
! contents :
!
! revision log :
!
! 24.11.2015    ggu     project started
!
!******************************************************************

!==================================================================
        module shympi
!==================================================================

	implicit none

	public

	logical, save :: bmpi = .false.

	integer,save :: n_threads = 1
	integer,save :: my_id = 0
	integer,save :: my_unit = 0

	integer,save :: nkn_global = 0		!total basin
	integer,save :: nel_global = 0
	integer,save :: nkn_local = 0		!this domain
	integer,save :: nel_local = 0
	integer,save :: nkn_inner = 0		!only proper, no ghost
	integer,save :: nel_inner = 0

	integer,save :: n_ghost_areas = 0
	integer,save :: n_ghost_nodes_max = 0
	integer,save :: n_ghost_elems_max = 0
	integer,save :: n_ghost_max = 0
	integer,save :: n_buffer = 0

	integer,save,allocatable :: ghost_areas(:,:)
	integer,save,allocatable :: ghost_nodes_in(:,:)
	integer,save,allocatable :: ghost_nodes_out(:,:)
	integer,save,allocatable :: ghost_elems(:,:)

	integer,save,allocatable :: i_buffer_in(:,:)
	integer,save,allocatable :: i_buffer_out(:,:)
	real,save,allocatable    :: r_buffer_in(:,:)
	real,save,allocatable    :: r_buffer_out(:,:)
	
	integer,save,allocatable :: node_area(:)	!global
	integer,save,allocatable :: request(:)		!for exchange
	integer,save,allocatable :: status(:,:)		!for exchange
	integer,save,allocatable :: ival(:)

	logical,save,allocatable :: is_inner_node(:)
	logical,save,allocatable :: is_inner_elem(:)
	integer,save,allocatable :: id_node(:)
	integer,save,allocatable :: id_elem(:,:)

        INTERFACE shympi_exchange_3d_node
        MODULE PROCEDURE  shympi_exchange_3d_node_r
     +                   ,shympi_exchange_3d_node_d
     +                   ,shympi_exchange_3d_node_i
        END INTERFACE

        INTERFACE shympi_exchange_2d_node
        MODULE PROCEDURE  shympi_exchange_2d_node_r
     +                   ,shympi_exchange_2d_node_d
     +                   ,shympi_exchange_2d_node_i
        END INTERFACE

        INTERFACE shympi_exchange_2d_elem
        MODULE PROCEDURE  shympi_exchange_2d_elem_r
     +                   ,shympi_exchange_2d_elem_d
     +                   ,shympi_exchange_2d_elem_i
        END INTERFACE

        INTERFACE shympi_check_2d_node
        MODULE PROCEDURE  shympi_check_2d_node_r
     +                   ,shympi_check_2d_node_d
     +                   ,shympi_check_2d_node_i
        END INTERFACE

        INTERFACE shympi_check_2d_elem
        MODULE PROCEDURE  shympi_check_2d_elem_r
     +                   ,shympi_check_2d_elem_d
     +                   ,shympi_check_2d_elem_i
        END INTERFACE

!==================================================================
        contains
!==================================================================

	subroutine shympi_init

	use basin

	integer ierr,size
	character*10 cunit
	character*80 file

	call shympi_init_internal(my_id,n_threads)

	nkn_global = nkn
	nel_global = nel
	nkn_local = nkn
	nel_local = nel
	nkn_inner = nkn
	nel_inner = nel

	bmpi = n_threads > 1

	call shympi_get_status_size_internal(size)

	allocate(node_area(nkn_global))
	allocate(ival(n_threads))
	allocate(request(2*n_threads))
	allocate(status(size,2*n_threads))

	if( bmpi ) then
	  write(cunit,'(i10)') my_id
	  cunit = adjustl(cunit)
	  file = 'mpi_debug_' // trim(cunit) // '.txt'
	  open(newunit=my_unit,file=file,status='unknown')
	  write(my_unit,*) 'shympi initialized: ',my_id,n_threads
	end if

	write(6,*) 'shympi initialized: ',my_id,n_threads
	flush(6)

	call shympi_barrier_internal
	call shympi_syncronize_initial
	call shympi_syncronize_internal

	end subroutine shympi_init

!******************************************************************

	subroutine shympi_alloc

	use basin

	allocate(is_inner_node(nkn))
	allocate(is_inner_elem(nel))
	allocate(id_node(nkn))
	allocate(id_elem(2,nel))

	is_inner_node = .true.
	is_inner_elem = .true.
	id_node = my_id
	id_elem = my_id

	end subroutine shympi_alloc

!******************************************************************

	subroutine shympi_alloc_ghost(n)

	use basin

	integer n

	allocate(ghost_areas(4,n_ghost_areas))
        allocate(ghost_nodes_out(n,n_ghost_areas))
        allocate(ghost_nodes_in(n,n_ghost_areas))
        allocate(ghost_elems(n,n_ghost_areas))

	ghost_areas = 0
        ghost_nodes_out = 0
        ghost_nodes_in = 0
        ghost_elems = 0

	end subroutine shympi_alloc_ghost

!******************************************************************

        subroutine shympi_alloc_buffer(n)

        integer n

        if( n_buffer >= n ) return

        if( n_buffer > 0 ) then
          deallocate(i_buffer_in)
          deallocate(i_buffer_out)
          deallocate(r_buffer_in)
          deallocate(r_buffer_out)
        end if

        n_buffer = n

        allocate(i_buffer_in(n_buffer,n_ghost_areas))
        allocate(i_buffer_out(n_buffer,n_ghost_areas))

        allocate(r_buffer_in(n_buffer,n_ghost_areas))
        allocate(r_buffer_out(n_buffer,n_ghost_areas))

        end subroutine shympi_alloc_buffer

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_barrier

	call shympi_barrier_internal

	end subroutine shympi_barrier

!******************************************************************

	subroutine shympi_stop(text)

	character*(*) text

	call shympi_finalize_internal
	write(6,*) text
	stop 'shympi_stop'

	end subroutine shympi_stop

!******************************************************************

	subroutine shympi_finalize

	call shympi_finalize_internal

	end subroutine shympi_finalize

!******************************************************************

	subroutine shympi_syncronize

	call shympi_syncronize_internal

	end subroutine shympi_syncronize

!******************************************************************

	subroutine shympi_abort

	call shympi_abort_internal

	end subroutine shympi_abort

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_3d_node_i(val)

	use basin
	use levels

	integer val(nlvdi,nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_i(belem,nlvdi,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_3d_node_i

!*******************************

	subroutine shympi_exchange_3d_node_r(val)

	use basin
	use levels

	real val(nlvdi,nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_r(belem,nlvdi,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_3d_node_r

!*******************************

	subroutine shympi_exchange_3d_node_d(val)

	use basin
	use levels

	double precision val(nlvdi,nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_d(belem,nlvdi,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_3d_node_d

!******************************************************************

	subroutine shympi_exchange_2d_node_i(val)

	use basin
	use levels

	integer val(nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_i(belem,1,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_2d_node_i

!*******************************

	subroutine shympi_exchange_2d_node_r(val)

	use basin
	use levels

	real val(nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_r(belem,1,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_2d_node_r

!*******************************

	subroutine shympi_exchange_2d_node_d(val)

	use basin
	use levels

	double precision val(nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_d(belem,1,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_2d_node_d

!******************************************************************

	subroutine shympi_exchange_2d_elem_i(val)

	use basin
	use levels

	integer val(nel)
	logical, parameter :: belem = .true.

	call shympi_exchange_internal_i(belem,1,nel,ilhv
     +			,ghost_elems,ghost_elems,val)

	end subroutine shympi_exchange_2d_elem_i

!*******************************

	subroutine shympi_exchange_2d_elem_r(val)

	use basin
	use levels

	real val(nel)
	logical, parameter :: belem = .true.

	call shympi_exchange_internal_r(belem,1,nel,ilhv
     +			,ghost_elems,ghost_elems,val)

	end subroutine shympi_exchange_2d_elem_r

!*******************************

	subroutine shympi_exchange_2d_elem_d(val)

	use basin
	use levels

	double precision val(nel)
	logical, parameter :: belem = .true.

	call shympi_exchange_internal_d(belem,1,nel,ilhv
     +			,ghost_elems,ghost_elems,val)

	end subroutine shympi_exchange_2d_elem_d

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine count_buffer(nlvddi,n,nc,il,nodes,nb)

        integer nlvddi,n,nc
        integer il(n)
        integer nodes(nc)
        integer nb

        integer i,k,l,lmax

        if( nlvddi == 1 ) then
          nb = nc
        else
          nb = 0
          do i=1,nc
            k = nodes(i)
            lmax = il(k)
            nb = nb + lmax
          end do
        end if

        end subroutine count_buffer

!******************************************************************

        subroutine to_buffer_i(nlvddi,n,nc,il,nodes,val,nb,buffer)

        integer nlvddi,n,nc
        integer il(n)
        integer nodes(nc)
        integer val(nlvddi,n)
        integer nb
        integer buffer(:)

        integer i,k,l,lmax

        if( nlvddi == 1 ) then
          do i=1,nc
            k = nodes(i)
            buffer(i) = val(1,k)
          end do
          nb = nc
        else
          nb = 0
          do i=1,nc
            k = nodes(i)
            lmax = il(k)
            do l=1,lmax
              nb = nb + 1
              buffer(nb) = val(l,k)
            end do
          end do
        end if

        end subroutine to_buffer_i

!******************************************************************

        subroutine from_buffer_i(nlvddi,n,nc,il,nodes,val,nb,buffer)

        integer nlvddi,n,nc
        integer il(n)
        integer nodes(nc)
        integer val(nlvddi,n)
        integer nb
        integer buffer(:)

        integer i,k,l,lmax

        if( nlvddi == 1 ) then
          do i=1,nc
            k = nodes(i)
            val(1,k) = buffer(i)
          end do
          nb = nc
        else
          nb = 0
          do i=1,nc
            k = nodes(i)
            lmax = il(k)
            do l=1,lmax
              nb = nb + 1
              val(l,k) = buffer(nb)
            end do
          end do
        end if

        end subroutine from_buffer_i

!******************************************************************

        subroutine to_buffer_r(nlvddi,n,nc,il,nodes,val,nb,buffer)

        integer nlvddi,n,nc
        integer il(n)
        integer nodes(nc)
        real val(nlvddi,n)
        integer nb
        real buffer(:)

        integer i,k,l,lmax

        if( nlvddi == 1 ) then
          do i=1,nc
            k = nodes(i)
            buffer(i) = val(1,k)
          end do
          nb = nc
        else
          nb = 0
          do i=1,nc
            k = nodes(i)
            lmax = il(k)
            do l=1,lmax
              nb = nb + 1
              buffer(nb) = val(l,k)
            end do
          end do
        end if

        end subroutine to_buffer_r

!******************************************************************

        subroutine from_buffer_r(nlvddi,n,nc,il,nodes,val,nb,buffer)

        integer nlvddi,n,nc
        integer il(n)
        integer nodes(nc)
        real val(nlvddi,n)
        integer nb
        real buffer(:)

        integer i,k,l,lmax

        if( nlvddi == 1 ) then
          do i=1,nc
            k = nodes(i)
            val(1,k) = buffer(i)
          end do
          nb = nc
        else
          nb = 0
          do i=1,nc
            k = nodes(i)
            lmax = il(k)
            do l=1,lmax
              nb = nb + 1
              val(l,k) = buffer(nb)
            end do
          end do
        end if

        end subroutine from_buffer_r

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_check_2d_node_i(val,text)

	use basin

	integer val(nkn)
	character*(*) text

	integer aux(nkn)

	aux = val
	call shympi_exchange_2d_node_i(aux)
	call shympi_check_array_i(nkn,val,aux,text)

	end subroutine shympi_check_2d_node_i

!*******************************

	subroutine shympi_check_2d_node_r(val,text)

	use basin

	real val(nkn)
	character*(*) text

	real aux(nkn)

	aux = val
	call shympi_exchange_2d_node_r(aux)
	call shympi_check_array_r(nkn,val,aux,text)

	end subroutine shympi_check_2d_node_r

!*******************************

	subroutine shympi_check_2d_node_d(val,text)

	use basin

	double precision val(nkn)
	character*(*) text

	double precision aux(nkn)

	aux = val
	call shympi_exchange_2d_node_d(aux)
	call shympi_check_array_d(nkn,val,aux,text)

	end subroutine shympi_check_2d_node_d

!******************************************************************

	subroutine shympi_check_2d_elem_i(val,text)

	use basin

	integer val(nel)
	character*(*) text

	integer aux(nel)

	aux = val
	call shympi_exchange_2d_elem_i(aux)
	call shympi_check_array_i(nel,val,aux,text)

	end subroutine shympi_check_2d_elem_i

!*******************************

	subroutine shympi_check_2d_elem_r(val,text)

	use basin

	real val(nel)
	character*(*) text

	real aux(nel)

	aux = val
	call shympi_exchange_2d_elem_r(aux)
	call shympi_check_array_r(nel,val,aux,text)

	end subroutine shympi_check_2d_elem_r

!*******************************

	subroutine shympi_check_2d_elem_d(val,text)

	use basin

	double precision val(nel)
	character*(*) text

	double precision aux(nel)

	aux = val
	call shympi_exchange_2d_elem_d(aux)
	call shympi_check_array_d(nel,val,aux,text)

	end subroutine shympi_check_2d_elem_d

!******************************************************************

	subroutine shympi_check_array_i(n,a1,a2,text)

	integer n
	integer a1(n),a2(n)
	character*(*) text

	integer i

        if( .not. all( a1 == a2 ) ) then
          write(6,*) 'arrays are different: ' // text
          write(6,*) 'process: ',my_id
	  do i=1,n
	    if( a1(i) /= a2(i) ) then
	      write(6,*) i,a1(i),a2(i)
	    end if
	  end do
	  call shympi_finalize
          stop 'error stop shympi_check_array_i'
        end if

	end subroutine shympi_check_array_i

!*******************************

	subroutine shympi_check_array_r(n,a1,a2,text)

	integer n
	real a1(n),a2(n)
	character*(*) text

	integer i

        if( .not. all( a1 == a2 ) ) then
          write(6,*) 'arrays are different: ' // text
          write(6,*) 'process: ',my_id
	  do i=1,n
	    if( a1(i) /= a2(i) ) then
	      write(6,*) i,a1(i),a2(i)
	    end if
	  end do
	  call shympi_finalize
          stop 'error stop shympi_check_array_i'
        end if

	end subroutine shympi_check_array_r

!*******************************

	subroutine shympi_check_array_d(n,a1,a2,text)

	integer n
	double precision a1(n),a2(n)
	character*(*) text

	integer i

        if( .not. all( a1 == a2 ) ) then
          write(6,*) 'arrays are different: ' // text
          write(6,*) 'process: ',my_id
	  do i=1,n
	    if( a1(i) /= a2(i) ) then
	      write(6,*) i,a1(i),a2(i)
	    end if
	  end do
	  call shympi_finalize
          stop 'error stop shympi_check_array_i'
        end if

	end subroutine shympi_check_array_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_gather_i(val)

	integer val

	call shympi_gather_i_internal(val)

	end subroutine shympi_gather_i

!*******************************

	subroutine shympi_bcast_i(val)

	integer val

	call shympi_bcast_i_internal(val)

	end subroutine shympi_bcast_i

!*******************************

	subroutine shympi_reduce_r(what,vals,val)

	character*(*) what
	real vals(:)
	real val

	if( what == 'min' ) then
	  val = MINVAL(vals)
	  call shympi_reduce_r_internal(what,val)
	else if( what == 'max' ) then
	  val = MAXVAL(vals)
	  call shympi_reduce_r_internal(what,val)
	else
	  write(6,*) 'what = ',what
	  stop 'error stop shympi_reduce_r: not ready'
	end if

	end subroutine shympi_reduce_r

!******************************************************************

	function shympi_output()

	logical shympi_output

	shympi_output = my_id == 0

	end function shympi_output

!******************************************************************

	function shympi_is_master()

	logical shympi_is_master

	shympi_is_master = my_id == 0

	end function shympi_is_master

!==================================================================
        end module shympi
!==================================================================

