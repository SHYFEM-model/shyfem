!
! mpi routines
!
! contents :
!
! revision log :
!
! 24.11.2015    ggu     project started
! 04.04.2018    ggu     new routine shympi_exit
! 04.04.2018    ggu     bug fix on scalar reduction (argument was changed)
! 10.04.2018    ggu     code to exchange arrays
!
!******************************************************************

!==================================================================
        module shympi
!==================================================================

	implicit none

	public

	logical, save :: bmpi = .false.
	logical, save :: bmpi_debug = .false.

	integer,save :: n_threads = 1
	integer,save :: my_id = 0
	integer,save :: my_unit = 0

	integer,save :: nkn_global = 0		!total basin
	integer,save :: nel_global = 0
	integer,save :: nkn_local = 0		!this domain
	integer,save :: nel_local = 0
	integer,save :: nkn_unique = 0		!this domain unique
	integer,save :: nel_unique = 0
	integer,save :: nkn_inner = 0		!only proper, no ghost
	integer,save :: nel_inner = 0

	integer,save :: n_ghost_areas = 0
	integer,save :: n_ghost_nodes_max = 0
	integer,save :: n_ghost_elems_max = 0
	integer,save :: n_ghost_max = 0
	integer,save :: n_buffer = 0

	integer,save,allocatable :: nkn_domains(:)
	integer,save,allocatable :: nel_domains(:)
	integer,save,allocatable :: nkn_cum_domains(:)
	integer,save,allocatable :: nel_cum_domains(:)

	integer,save,allocatable :: ghost_areas(:,:)
	integer,save,allocatable :: ghost_nodes_in(:,:)
	integer,save,allocatable :: ghost_nodes_out(:,:)
	integer,save,allocatable :: ghost_elems_in(:,:)
	integer,save,allocatable :: ghost_elems_out(:,:)

	integer,save,allocatable :: i_buffer_in(:,:)
	integer,save,allocatable :: i_buffer_out(:,:)
	real,save,allocatable    :: r_buffer_in(:,:)
	real,save,allocatable    :: r_buffer_out(:,:)
	
	integer,save,allocatable :: request(:)		!for exchange
	integer,save,allocatable :: status(:,:)		!for exchange

	integer,save,allocatable :: id_node(:)
	integer,save,allocatable :: id_elem(:,:)

        type communication_info
          integer, public :: numberID
          integer, public :: totalID
          !contains the global ID of the elements belonging to the process
          integer, public, dimension(:), allocatable :: globalID,localID
          !contains the rank ID of the process to which belong the global ID
          !elements with which to communicate
          integer, public, dimension(:,:), allocatable :: rankID
          integer, public, dimension(:,:), allocatable :: neighborID
        end type communication_info

	type (communication_info), SAVE :: univocal_nodes
        integer, allocatable, save, dimension(:) :: allPartAssign

!---------------------

        INTERFACE shympi_exchange_3d_node
        	MODULE PROCEDURE  
     +			  shympi_exchange_3d_node_r
     +                   ,shympi_exchange_3d_node_d
     +                   ,shympi_exchange_3d_node_i
        END INTERFACE

        INTERFACE shympi_exchange_3d0_node
        	MODULE PROCEDURE  
     +			  shympi_exchange_3d0_node_r
!     +                   ,shympi_exchange_3d0_node_d
!     +                   ,shympi_exchange_3d0_node_i
        END INTERFACE

        INTERFACE shympi_exchange_2d_node
        	MODULE PROCEDURE  
     +			  shympi_exchange_2d_node_r
     +                   ,shympi_exchange_2d_node_d
     +                   ,shympi_exchange_2d_node_i
        END INTERFACE

        INTERFACE shympi_exchange_3d_elem
        	MODULE PROCEDURE  
     +			  shympi_exchange_3d_elem_r
!     +                   ,shympi_exchange_3d_elem_d
!     +                   ,shympi_exchange_3d_elem_i
        END INTERFACE

        INTERFACE shympi_exchange_2d_elem
        	MODULE PROCEDURE  
     +			  shympi_exchange_2d_elem_r
     +                   ,shympi_exchange_2d_elem_d
     +                   ,shympi_exchange_2d_elem_i
        END INTERFACE

!---------------------

        INTERFACE shympi_check_elem
        	MODULE PROCEDURE  
     +			  shympi_check_2d_elem_r
     +                   ,shympi_check_2d_elem_d
     +                   ,shympi_check_2d_elem_i
     +			 ,shympi_check_3d_elem_r
!     +                   ,shympi_check_3d_elem_d
!     +                   ,shympi_check_3d_elem_i
        END INTERFACE

        INTERFACE shympi_check_node
        	MODULE PROCEDURE  
     +			  shympi_check_2d_node_r
     +                   ,shympi_check_2d_node_d
     +                   ,shympi_check_2d_node_i
     +			 ,shympi_check_3d_node_r
!     +                   ,shympi_check_3d_node_d
!     +                   ,shympi_check_3d_node_i
        END INTERFACE

        INTERFACE shympi_check_2d_node
        	MODULE PROCEDURE  
     +			  shympi_check_2d_node_r
     +                   ,shympi_check_2d_node_d
     +                   ,shympi_check_2d_node_i
        END INTERFACE

        INTERFACE shympi_check_2d_elem
        	MODULE PROCEDURE  
     +			  shympi_check_2d_elem_r
     +                   ,shympi_check_2d_elem_d
     +                   ,shympi_check_2d_elem_i
        END INTERFACE

        INTERFACE shympi_check_3d_node
        	MODULE PROCEDURE  
     +			  shympi_check_3d_node_r
!     +                   ,shympi_check_3d_node_d
!     +                   ,shympi_check_3d_node_i
        END INTERFACE

        INTERFACE shympi_check_3d0_node
        	MODULE PROCEDURE  
     +			  shympi_check_3d0_node_r
!     +                   ,shympi_check_3d0_node_d
!     +                   ,shympi_check_3d0_node_i
        END INTERFACE

        INTERFACE shympi_check_3d_elem
        	MODULE PROCEDURE  
     +			  shympi_check_3d_elem_r
!     +                   ,shympi_check_3d_elem_d
!     +                   ,shympi_check_3d_elem_i
        END INTERFACE

        INTERFACE shympi_check_array
                MODULE PROCEDURE
     +                    shympi_check_array_i
     +                   ,shympi_check_array_r
     +                   ,shympi_check_array_d
        END INTERFACE

!---------------------

        INTERFACE shympi_gather
                MODULE PROCEDURE
     +                    shympi_gather_scalar_i
     +                   ,shympi_gather_array_r
     +                   ,shympi_gather_array_d
        END INTERFACE

        INTERFACE shympi_gather_and_sum
                MODULE PROCEDURE
     +                    shympi_gather_and_sum_i
     +                   ,shympi_gather_and_sum_r
     +                   ,shympi_gather_and_sum_d
        END INTERFACE

        INTERFACE shympi_bcast
                MODULE PROCEDURE
     +                    shympi_bcast_scalar_i
     +                   ,shympi_bcast_array_r
        END INTERFACE

        INTERFACE shympi_collect_node_value
                MODULE PROCEDURE
     +                     shympi_collect_node_value_2d_i
     +                    ,shympi_collect_node_value_2d_r
     +                    ,shympi_collect_node_value_3d_r
!     +                    ,shympi_collect_node_value_2d_i
        END INTERFACE

        INTERFACE shympi_reduce
                MODULE PROCEDURE
     +                    shympi_reduce_r
!     +                   ,shympi_reduce_i
        END INTERFACE

!---------------------

        INTERFACE shympi_min
        	MODULE PROCEDURE  
     +			   shympi_min_r
     +			  ,shympi_min_i
     +			  ,shympi_min_d
     +			  ,shympi_min_0_r
     +			  ,shympi_min_0_i
     +			  ,shympi_min_0_d
        END INTERFACE

        INTERFACE shympi_max
        	MODULE PROCEDURE  
     +			   shympi_max_r
     +			  ,shympi_max_i
!     +			  ,shympi_max_d
     +			  ,shympi_max_0_r
     +			  ,shympi_max_0_i
!     +			  ,shympi_max_0_d
        END INTERFACE

        INTERFACE shympi_sum
        	MODULE PROCEDURE  
     +			   shympi_sum_r
     +			  ,shympi_sum_i
     +			  ,shympi_sum_d
     +			  ,shympi_sum_0_r
     +			  ,shympi_sum_0_i
     +			  ,shympi_sum_0_d
        END INTERFACE

!---------------------

        INTERFACE shympi_exchange_array
        	MODULE PROCEDURE  
     +			   shympi_exchange_array_2d_r
     +			  ,shympi_exchange_array_2d_i
     +			  ,shympi_exchange_array_3d_r
     +			  ,shympi_exchange_array_3d_i
        END INTERFACE

        INTERFACE shympi_get_array
        	MODULE PROCEDURE  
     +			   shympi_get_array_2d_r
     +			  ,shympi_get_array_2d_i
!     +			  ,shympi_get_array_3d_r
!     +			  ,shympi_get_array_3d_i
        END INTERFACE

        INTERFACE shympi_getvals
        	MODULE PROCEDURE  
     +			   shympi_getvals_2d_node_r
     +			  ,shympi_getvals_2d_node_i
     +			  ,shympi_getvals_3d_node_r
     +			  ,shympi_getvals_3d_node_i
        END INTERFACE

!---------------------

        INTERFACE shympi_exchange_and_sum_3d_nodes
        	MODULE PROCEDURE  
     +			   shympi_exchange_and_sum_3d_nodes_r
     +			  ,shympi_exchange_and_sum_3d_nodes_d
        END INTERFACE

        INTERFACE shympi_exchange_and_sum_2d_nodes
        	MODULE PROCEDURE  
     +			   shympi_exchange_and_sum_2d_nodes_r
     +			  ,shympi_exchange_and_sum_2d_nodes_d
        END INTERFACE

        INTERFACE shympi_exchange_2d_nodes_min
        	MODULE PROCEDURE  
     +			   shympi_exchange_2d_nodes_min_i
     +			  ,shympi_exchange_2d_nodes_min_r
        END INTERFACE

        INTERFACE shympi_exchange_2d_nodes_max
        	MODULE PROCEDURE  
     +			   shympi_exchange_2d_nodes_max_i
     +			  ,shympi_exchange_2d_nodes_max_r
        END INTERFACE

!==================================================================
        contains
!==================================================================

	subroutine shympi_init(b_want_mpi)

	use basin

	logical b_want_mpi

	integer ierr,size
	character*10 cunit
	character*80 file

	!-----------------------------------------------------
	! initializing
	!-----------------------------------------------------

	nkn_global = nkn
	nel_global = nel
	nkn_local = nkn
	nel_local = nel
	nkn_inner = nkn
	nel_inner = nel
	nkn_unique = nkn
	nel_unique = nel

	write(6,*) 'calling shympi_init_internal'
	call shympi_init_internal(my_id,n_threads)
	!call check_part_basin('nodes')

	bmpi = n_threads > 1

	if( b_want_mpi .and. .not. bmpi ) then
	  write(6,*) 'program wants mpi but only one thread available'
	  stop 'error stop shympi_init'
	end if

	if( .not. b_want_mpi .and. bmpi ) then
	  write(6,*) 'program does not want mpi '//
     +			'but program running in mpi mode'
	  stop 'error stop shympi_init'
	end if

	!-----------------------------------------------------
	! allocate important arrays
	!-----------------------------------------------------

	call shympi_get_status_size_internal(size)

	allocate(request(2*n_threads))
	allocate(status(size,2*n_threads))
        allocate(nkn_domains(n_threads))
        allocate(nel_domains(n_threads))
        allocate(nkn_cum_domains(0:n_threads))
        allocate(nel_cum_domains(0:n_threads))
	nkn_domains(1) = nkn
	nel_domains(1) = nel
	nkn_cum_domains(0) = 0
	nkn_cum_domains(1) = nkn
	nel_cum_domains(0) = 0
	nel_cum_domains(1) = nel

	!-----------------------------------------------------
	! next is needed if program is not running in mpi mode
	!-----------------------------------------------------

	if( .not. bmpi ) call shympi_alloc_id(nkn,nel)

	!-----------------------------------------------------
	! output to terminal
	!-----------------------------------------------------

	if( bmpi ) then
	  write(cunit,'(i10)') my_id
	  cunit = adjustl(cunit)
	  file = 'mpi_debug_' // trim(cunit) // '.txt'
	  call shympi_get_new_unit(my_unit)
	  open(unit=my_unit,file=file,status='unknown')
	  write(my_unit,*) 'shympi initialized: ',my_id,n_threads,my_unit
	  !write(6,*) '######### my_unit ',my_unit,my_id
	else
	  write(6,*) 'shympi initialized: ',my_id,n_threads
	  write(6,*) 'shympi is not running in mpi mode'
	end if

	flush(6)

	!-----------------------------------------------------
	! finish up
	!-----------------------------------------------------

	call shympi_barrier_internal
	call shympi_syncronize_initial
	call shympi_syncronize_internal

	!-----------------------------------------------------
	! end of routine
	!-----------------------------------------------------

	end subroutine shympi_init

!******************************************************************

	subroutine shympi_alloc_id(nk,ne)

	integer nk,ne

	write(6,*) 'shympi_alloc_id: ',nk,ne

	allocate(id_node(nk))
	allocate(id_elem(0:2,ne))

	id_node = my_id
	id_elem = my_id

	end subroutine shympi_alloc_id

!******************************************************************

	subroutine shympi_alloc_ghost(n)

	use basin

	integer n

	allocate(ghost_areas(5,n_ghost_areas))
        allocate(ghost_nodes_out(n,n_ghost_areas))
        allocate(ghost_nodes_in(n,n_ghost_areas))
        allocate(ghost_elems_out(n,n_ghost_areas))
        allocate(ghost_elems_in(n,n_ghost_areas))

	ghost_areas = 0
        ghost_nodes_out = 0
        ghost_nodes_in = 0
        ghost_elems_out = 0
        ghost_elems_in = 0

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

	subroutine shympi_get_new_unit(iunit)

	integer iunit

	integer iu,iumin,iumax,ios
        logical opened

	iumin = 20
	iumax = 1000

	do iu=iumin,iumax
          inquire (unit=iu, opened=opened, iostat=ios)
          if (ios.ne.0) cycle
          if (.not.opened) exit
	end do

	if( iu > iumax ) then
	  iu = 0
	  stop 'error stop shympi_get_new_unit: no new unit'
	end if

	iunit = iu

	end subroutine shympi_get_new_unit

!******************************************************************
!******************************************************************
!******************************************************************

	function shympi_partition_on_elements()

	logical shympi_partition_on_elements

	shympi_partition_on_elements = .false.

	end function shympi_partition_on_elements

!******************************************************************

        function shympi_partition_on_nodes()

        logical shympi_partition_on_nodes

        shympi_partition_on_nodes = .true.

        end function shympi_partition_on_nodes

!******************************************************************
!******************************************************************
!******************************************************************

        function shympi_is_inner_node(k)

        logical shympi_is_inner_node
        integer k

        shympi_is_inner_node = ( k <= nkn_inner )

        end function shympi_is_inner_node

!******************************************************************

        function shympi_is_inner_elem(ie)

        logical shympi_is_inner_elem
        integer ie

        shympi_is_inner_elem = ( ie <= nel_inner )

        end function shympi_is_inner_elem

!******************************************************************

        function shympi_is_unique_node(k)

        logical shympi_is_unique_node
        integer k

        shympi_is_unique_node = ( k <= nkn_unique )

        end function shympi_is_unique_node

!******************************************************************

        function shympi_is_unique_elem(ie)

        logical shympi_is_unique_elem
        integer ie

        shympi_is_unique_elem = ( ie <= nel_unique )

        end function shympi_is_unique_elem

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_barrier

	call shympi_barrier_internal

	end subroutine shympi_barrier

!******************************************************************

	subroutine shympi_syncronize

	call shympi_syncronize_internal

	end subroutine shympi_syncronize

!******************************************************************

	subroutine shympi_finalize

	call shympi_barrier_internal
	call shympi_finalize_internal
	stop

	end subroutine shympi_finalize

!******************************************************************

	subroutine shympi_exit(ierr)

	integer ierr

	call shympi_barrier_internal
	call shympi_finalize_internal
	call exit(ierr)
	stop

	end subroutine shympi_exit

!******************************************************************

	subroutine shympi_stop(text)

	character*(*) text

	if( shympi_is_master() ) then
	  write(6,*) 'error stop shympi_stop: ',trim(text)
	end if
	call shympi_barrier_internal
	call shympi_finalize_internal
	stop

	end subroutine shympi_stop

!******************************************************************

	subroutine shympi_abort

	call shympi_abort_internal(33)
	stop

	end subroutine shympi_abort

!******************************************************************

	function shympi_wtime()

	double precision shympi_wtime
	double precision shympi_wtime_internal

	shympi_wtime = shympi_wtime_internal()

	end function shympi_wtime

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_3d_node_i(val)

	use basin
	use levels

	integer val(nlvdi,nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_i(belem,1,nlvdi,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_3d_node_i

!*******************************

	subroutine shympi_exchange_3d_node_r(val)

	use basin
	use levels

	real val(nlvdi,nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_r(belem,1,nlvdi,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_3d_node_r

!*******************************

	subroutine shympi_exchange_3d_node_d(val)

	use basin
	use levels

	double precision val(nlvdi,nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_d(belem,1,nlvdi,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_3d_node_d

!******************************************************************

	subroutine shympi_exchange_3d0_node_r(val)

	use basin
	use levels

	real val(0:nlvdi,nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_r(belem,0,nlvdi,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_3d0_node_r

!******************************************************************

	subroutine shympi_exchange_3d_elem_r(val)

	use basin
	use levels

	real val(nlvdi,nel)
	logical, parameter :: belem = .true.

	call shympi_exchange_internal_r(belem,1,nlvdi,nel,ilhv
     +			,ghost_elems_in,ghost_elems_out,val)

	end subroutine shympi_exchange_3d_elem_r

!******************************************************************

	subroutine shympi_exchange_2d_node_i(val)

	use basin
	use levels

	integer val(nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_i(belem,1,1,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_2d_node_i

!*******************************

	subroutine shympi_exchange_2d_node_r(val)

	use basin
	use levels

	real val(nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_r(belem,1,1,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_2d_node_r

!*******************************

	subroutine shympi_exchange_2d_node_d(val)

	use basin
	use levels

	double precision val(nkn)
	logical, parameter :: belem = .false.

	call shympi_exchange_internal_d(belem,1,1,nkn,ilhkv
     +			,ghost_nodes_in,ghost_nodes_out,val)

	end subroutine shympi_exchange_2d_node_d

!******************************************************************

	subroutine shympi_exchange_2d_elem_i(val)

	use basin
	use levels

	integer val(nel)
	logical, parameter :: belem = .true.

	call shympi_exchange_internal_i(belem,1,1,nel,ilhv
     +			,ghost_elems_in,ghost_elems_out,val)

	end subroutine shympi_exchange_2d_elem_i

!*******************************

	subroutine shympi_exchange_2d_elem_r(val)

	use basin
	use levels

	real val(nel)
	logical, parameter :: belem = .true.

	call shympi_exchange_internal_r(belem,1,1,nel,ilhv
     +			,ghost_elems_in,ghost_elems_out,val)

	end subroutine shympi_exchange_2d_elem_r

!*******************************

	subroutine shympi_exchange_2d_elem_d(val)

	use basin
	use levels

	double precision val(nel)
	logical, parameter :: belem = .true.

	call shympi_exchange_internal_d(belem,1,1,nel,ilhv
     +			,ghost_elems_in,ghost_elems_out,val)

	end subroutine shympi_exchange_2d_elem_d

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

	subroutine shympi_check_3d_node_r(val,text)

	use basin
	use levels

	real val(nlvdi,nkn)
	character*(*) text

	real aux(nlvdi,nkn)

	aux = val
	call shympi_exchange_3d_node_r(aux)
	call shympi_check_array_r(nlvdi*nkn,val,aux,text)

	end subroutine shympi_check_3d_node_r

!******************************************************************

	subroutine shympi_check_3d0_node_r(val,text)

	use basin
	use levels

	real val(0:nlvdi,nkn)
	character*(*) text

	real aux(0:nlvdi,nkn)

	aux = val
	call shympi_exchange_3d0_node_r(aux)
	call shympi_check_array_r((nlvdi+1)*nkn,val,aux,text)

	end subroutine shympi_check_3d0_node_r

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

	subroutine shympi_check_3d_elem_r(val,text)

	use basin
	use levels

	real val(nlvdi,nel)
	character*(*) text

	real aux(nlvdi,nel)

	aux = val
	call shympi_exchange_3d_elem_r(aux)
	call shympi_check_array_r(nlvdi*nel,val,aux,text)

	end subroutine shympi_check_3d_elem_r

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_check_array_i(n,a1,a2,text)

	integer n
	integer a1(n),a2(n)
	character*(*) text

	integer i

        if( .not. all( a1 == a2 ) ) then
          write(6,*) 'arrays are different: ' // text
          write(6,*) 'process id: ',my_id
          write(6,*) 'total array size: ',n
	  do i=1,n
	    if( a1(i) /= a2(i) ) then
	      write(6,*) my_id,i,a1(i),a2(i)
	    end if
	  end do
	  call shympi_abort
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
          write(6,*) 'process id: ',my_id
          write(6,*) 'total array size: ',n
	  do i=1,n
	    if( a1(i) /= a2(i) ) then
	      write(6,*) my_id,i,a1(i),a2(i)
	    end if
	  end do
	  call shympi_abort
          stop 'error stop shympi_check_array_r'
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
          write(6,*) 'process id: ',my_id
          write(6,*) 'total array size: ',n
	  do i=1,n
	    if( a1(i) /= a2(i) ) then
	      write(6,*) my_id,i,a1(i),a2(i)
	    end if
	  end do
	  call shympi_abort
          stop 'error stop shympi_check_array_d'
        end if

	end subroutine shympi_check_array_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_gather_scalar_i(val,vals)

	integer val
	integer vals(n_threads)

	integer n

	n = 1
	call shympi_allgather_i_internal(n,val,vals)

	end subroutine shympi_gather_scalar_i

!*******************************

	subroutine shympi_gather_array_r(val,vals)

	real val(:)
	real vals(size(val),n_threads)

	integer n

	n = size(val)
	call shympi_allgather_r_internal(n,val,vals)

	end subroutine shympi_gather_array_r

!*******************************

	subroutine shympi_gather_array_d(val,vals)

	double precision val(:)
	double precision vals(size(val),n_threads)

	integer n

	n = size(val)
	call shympi_allgather_d_internal(n,val,vals)

	end subroutine shympi_gather_array_d

!*******************************

	subroutine shympi_gather_and_sum_i(val)

	integer val(:)

	integer n
	integer vals(size(val),n_threads)

	n = size(val)
	call shympi_allgather_i_internal(n,val,vals)
	val(:) = SUM(vals,dim=2)

	end subroutine shympi_gather_and_sum_i

!*******************************

	subroutine shympi_gather_and_sum_r(val)

	real val(:)

	integer n
	real vals(size(val),n_threads)

	n = size(val)
	call shympi_allgather_r_internal(n,val,vals)
	val(:) = SUM(vals,dim=2)

	end subroutine shympi_gather_and_sum_r

!*******************************

	subroutine shympi_gather_and_sum_d(val)

	double precision val(:)

	integer n
	double precision vals(size(val),n_threads)

	n = size(val)
	call shympi_allgather_d_internal(n,val,vals)
	val(:) = SUM(vals,dim=2)

	end subroutine shympi_gather_and_sum_d

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_bcast_scalar_i(val)

	integer val

	integer n

	n = 1
	call shympi_bcast_i_internal(n,val)

	end subroutine shympi_bcast_scalar_i

!*******************************

	subroutine shympi_bcast_array_r(val)

	integer val(:)

	integer n

	n = size(val)
	call shympi_bcast_r_internal(n,val)

	end subroutine shympi_bcast_array_r

!*******************************

	subroutine shympi_collect_node_value_2d_r(k,vals,val)

	integer k
	real vals(:)
	real val

	val = 0
	if( k > 0 .and. k <= nkn_unique ) val = vals(k)
	if( bmpi ) val = shympi_sum(val)

	end subroutine shympi_collect_node_value_2d_r

!*******************************

	subroutine shympi_collect_node_value_2d_i(k,vals,val)

	integer k
	integer vals(:)
	integer val

	val = 0
	if( k > 0 .and. k <= nkn_unique ) val = vals(k)
	if( bmpi ) val = shympi_sum(val)

	end subroutine shympi_collect_node_value_2d_i

!*******************************

	subroutine shympi_collect_node_value_3d_r(k,vals,val)

	integer k
	real vals(:,:)
	real val(:)

	real vaux(size(val),n_threads)

	val = 0
	if( k > 0 .and. k <= nkn_unique ) val(:) = vals(:,k)
	if( bmpi ) then
	  call shympi_gather(val,vaux)
	  val = SUM(vaux,dim=2)
	end if

	end subroutine shympi_collect_node_value_3d_r

!*******************************

!	subroutine shympi_collect_node_value_3d_i(k,vals,val)
!
!	integer k
!	integer vals(:)
!	integer val
!
!	val = 0
!	if( k > 0 .and. k <= nkn_unique ) val = vals(k)
!	if( bmpi ) val = shympi_sum(val)
!
!	end subroutine shympi_collect_node_value_3d_i

!*******************************

	subroutine shympi_find_node(ke,ki,id)

	use basin

	integer ke,ki,id

	integer k,ic,i
	integer vals(n_threads)

	do k=1,nkn_unique
	  if( ipv(k) == ke ) exit
	end do
	if( k > nkn_unique ) k = 0

	call shympi_allgather_i_internal(1,k,vals)

	ic = count( vals /= 0 )
	if( ic /= 1 ) then
	  write(6,*) 'node found in more than one domain: '
	  write(6,*) '==========================='
	  write(6,*) n_threads,my_id
	  write(6,*) vals
	  write(6,*) '==========================='
	  call shympi_finalize
	  stop 'error stop shympi_find_node: more than one domain'
	end if

	do i=1,n_threads
	  if( vals(i) /= 0 ) exit
	end do

	ki = vals(i)
	id = i

	end subroutine shympi_find_node

!*******************************

	subroutine shympi_reduce_r(what,vals,val)

	character*(*) what
	real vals(:)
	real val

	if( what == 'min' ) then
	  val = MINVAL(vals)
	  if( bmpi ) call shympi_reduce_r_internal(what,val)
	else if( what == 'max' ) then
	  val = MAXVAL(vals)
	  if( bmpi ) call shympi_reduce_r_internal(what,val)
	else if( what == 'sum' ) then
	  val = SUM(vals)
	  if( bmpi ) call shympi_reduce_r_internal(what,val)
	else
	  write(6,*) 'what = ',what
	  stop 'error stop shympi_reduce_r: not ready'
	end if

	end subroutine shympi_reduce_r

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_get_array_2d_r(n,vals,val_out)

	use basin
	use levels

	integer n
	real vals(:)
	real val_out(n)

!	call shympi_getvals_internal_r(kind,1,nkn
!     +                                    ,vals,val)
	stop 'error stop shympi_get_array_2d_r: not ready'

	end subroutine shympi_get_array_2d_r

!*******************************

	subroutine shympi_get_array_2d_i(n,vals,val_out)

	use basin
	use levels

	integer n
	integer vals(:)
	integer val_out(n)

	call shympi_get_array_internal_i(1,n
     +                                    ,vals,val_out)

	end subroutine shympi_get_array_2d_i

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_array_3d_r(vals,val_out)

	real vals(:,:)
	real val_out(:,:)

	integer n,ni1,no1

	ni1 = size(vals,1)
	no1 = size(val_out,1)
	n = size(val_out,2)

	if( ni1 /= no1 ) stop 'error stop exchange: first dimension'

	call shympi_exchange_array_internal_r(ni1,n
     +                                    ,vals,val_out)

	end subroutine shympi_exchange_array_3d_r

!*******************************

	subroutine shympi_exchange_array_3d_i(vals,val_out)

	integer vals(:,:)
	integer val_out(:,:)

	integer n,ni1,no1

	ni1 = size(vals,1)
	no1 = size(val_out,1)
	n = size(val_out,2)

	if( ni1 /= no1 ) stop 'error stop exchange: first dimension'

	call shympi_exchange_array_internal_i(ni1,n
     +                                    ,vals,val_out)

	end subroutine shympi_exchange_array_3d_i

!*******************************

	subroutine shympi_exchange_array_2d_r(vals,val_out)

	real vals(:)
	real val_out(:)

	integer n

	n = size(val_out)

	call shympi_exchange_array_internal_r(1,n
     +                                    ,vals,val_out)

	end subroutine shympi_exchange_array_2d_r

!*******************************

	subroutine shympi_exchange_array_2d_i(vals,val_out)

	integer vals(:)
	integer val_out(:)

	integer n

	n = size(val_out)

	call shympi_exchange_array_internal_i(1,n
     +                                    ,vals,val_out)

	end subroutine shympi_exchange_array_2d_i

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_getvals_2d_node_r(kind,vals,val)

	use basin
	use levels

	integer kind(2)
	real vals(nkn)
	real val

	call shympi_getvals_internal_r(kind,1,nkn
     +                                    ,vals,val)

	end subroutine shympi_getvals_2d_node_r

!*******************************

	subroutine shympi_getvals_2d_node_i(kind,vals,val)

	use basin
	use levels

	integer kind(2)
	integer vals(nkn)
	integer val

	call shympi_getvals_internal_i(kind,1,nkn
     +                                    ,vals,val)

	end subroutine shympi_getvals_2d_node_i

!*******************************

	subroutine shympi_getvals_3d_node_r(kind,vals,val)

	use basin
	use levels

	integer kind(2)
	real vals(nlvdi,nkn)
	real val(nlvdi)

	call shympi_getvals_internal_r(kind,nlvdi,nkn
     +                                    ,vals,val)

	end subroutine shympi_getvals_3d_node_r

!*******************************

	subroutine shympi_getvals_3d_node_i(kind,vals,val)

	use basin
	use levels

	integer kind(2)
	integer vals(nlvdi,nkn)
	integer val(nlvdi)

	call shympi_getvals_internal_i(kind,nlvdi,nkn
     +                                    ,vals,val)

	end subroutine shympi_getvals_3d_node_i

!******************************************************************
!******************************************************************
!******************************************************************

	function shympi_min_d(vals)

	double precision shympi_min_d
	double precision vals(:)
	double precision val

	val = MINVAL(vals)
	if( bmpi ) call shympi_reduce_d_internal('min',val)

	shympi_min_d = val

	end function shympi_min_d

!******************************************************************

	function shympi_min_r(vals)

	real shympi_min_r
	real vals(:)
	real val

	val = MINVAL(vals)
	if( bmpi ) call shympi_reduce_r_internal('min',val)

	shympi_min_r = val

	end function shympi_min_r

!******************************************************************

	function shympi_min_i(vals)

	integer shympi_min_i
	integer vals(:)
	integer val

	val = MINVAL(vals)
	if( bmpi ) call shympi_reduce_i_internal('min',val)

	shympi_min_i = val

	end function shympi_min_i

!******************************************************************

	function shympi_min_0_d(val0)

	double precision shympi_min_0_d
	double precision val0
	double precision val

	val = val0
	if( bmpi ) call shympi_reduce_d_internal('min',val)

	shympi_min_0_d = val

	end function shympi_min_0_d

!******************************************************************

	function shympi_min_0_r(val0)

	real shympi_min_0_r
	real val0
	real val

	val = val0
	if( bmpi ) call shympi_reduce_r_internal('min',val)

	shympi_min_0_r = val

	end function shympi_min_0_r

!******************************************************************

	function shympi_min_0_i(val0)

	integer shympi_min_0_i
	integer val0
	integer val

	val = val0
	if( bmpi ) call shympi_reduce_i_internal('min',val)

	shympi_min_0_i = val

	end function shympi_min_0_i

!******************************************************************

	function shympi_max_r(vals)

	real shympi_max_r
	real vals(:)
	real val

	val = MAXVAL(vals)
	if( bmpi ) call shympi_reduce_r_internal('max',val)

	shympi_max_r = val

	end function shympi_max_r

!******************************************************************

	function shympi_max_i(vals)

	integer shympi_max_i
	integer vals(:)
	integer val

	val = MAXVAL(vals)
	if( bmpi ) call shympi_reduce_i_internal('max',val)

	shympi_max_i = val

	end function shympi_max_i

!******************************************************************

	function shympi_max_0_i(val0)

! routine for val that is scalar

	integer shympi_max_0_i
	integer val0
	integer val

	val = val0
	if( bmpi ) call shympi_reduce_i_internal('max',val)

	shympi_max_0_i = val

	end function shympi_max_0_i

!******************************************************************

	function shympi_max_0_r(val0)

! routine for val that is scalar

	real shympi_max_0_r
	real val0
	real val

	val = val0
	if( bmpi ) call shympi_reduce_r_internal('max',val)

	shympi_max_0_r = val

	end function shympi_max_0_r

!******************************************************************

	function shympi_sum_d(vals)

	double precision shympi_sum_d
	double precision vals(:)
	double precision val

	val = SUM(vals)
	if( bmpi ) call shympi_reduce_d_internal('sum',val)

	shympi_sum_d = val

	end function shympi_sum_d

!******************************************************************

	function shympi_sum_r(vals)

	real shympi_sum_r
	real vals(:)
	real val

	val = SUM(vals)
	if( bmpi ) call shympi_reduce_r_internal('sum',val)

	shympi_sum_r = val

	end function shympi_sum_r

!******************************************************************

	function shympi_sum_i(vals)

	integer shympi_sum_i
	integer vals(:)
	integer val

	val = SUM(vals)
	if( bmpi ) call shympi_reduce_i_internal('sum',val)

	shympi_sum_i = val

	end function shympi_sum_i

!******************************************************************

	function shympi_sum_0_d(val0)

	double precision shympi_sum_0_d
	double precision val0
	double precision val

	val = val0
	if( bmpi ) call shympi_reduce_d_internal('sum',val)

	shympi_sum_0_d = val

	end function shympi_sum_0_d

!******************************************************************

	function shympi_sum_0_r(val0)

	real shympi_sum_0_r
	real val0
	real val

	val = val0
	if( bmpi ) call shympi_reduce_r_internal('sum',val)

	shympi_sum_0_r = val

	end function shympi_sum_0_r

!******************************************************************

	function shympi_sum_0_i(val0)

	integer shympi_sum_0_i
	integer val0
	integer val

	val = val0
	if( bmpi ) call shympi_reduce_i_internal('sum',val)

	shympi_sum_0_i = val

	end function shympi_sum_0_i

!******************************************************************
!******************************************************************
!******************************************************************
! next are routines for partition on nodes - empty here
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_and_sum_3d_nodes_r(a)
	real a(:,:)
	end subroutine shympi_exchange_and_sum_3d_nodes_r

	subroutine shympi_exchange_and_sum_3d_nodes_d(a)
	double precision a(:,:)
	end subroutine shympi_exchange_and_sum_3d_nodes_d

	subroutine shympi_exchange_and_sum_2d_nodes_r(a)
	real a(:)
	end subroutine shympi_exchange_and_sum_2d_nodes_r

	subroutine shympi_exchange_and_sum_2d_nodes_d(a)
	double precision a(:)
	end subroutine shympi_exchange_and_sum_2d_nodes_d

	subroutine shympi_exchange_2d_nodes_min_i(a)
	integer a(:)
	end subroutine shympi_exchange_2d_nodes_min_i

	subroutine shympi_exchange_2d_nodes_max_i(a)
	integer a(:)
	end subroutine shympi_exchange_2d_nodes_max_i

	subroutine shympi_exchange_2d_nodes_min_r(a)
	real a(:)
	end subroutine shympi_exchange_2d_nodes_min_r

	subroutine shympi_exchange_2d_nodes_max_r(a)
	real a(:)
	end subroutine shympi_exchange_2d_nodes_max_r

!******************************************************************
!******************************************************************
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

!******************************************************************

	subroutine shympi_comment(text)

	character*(*) text

	!if( bmpi .and. bmpi_debug .and. my_id == 0 ) then
	if( bmpi_debug .and. my_id == 0 ) then
	  write(6,*) 'shympi_comment: ' // trim(text)
	  write(299,*) 'shympi_comment: ' // trim(text)
	end if

	end subroutine shympi_comment

!******************************************************************

	function shympi_is_parallel()

	logical shympi_is_parallel

	shympi_is_parallel = bmpi

	end function shympi_is_parallel

!******************************************************************

	function shympi_can_parallel()

	logical shympi_can_parallel

	shympi_can_parallel = .true.

	end function shympi_can_parallel

!******************************************************************

	subroutine shympi_univocal_nodes

	use basin, only : nkn

	implicit none

	integer ierr

	integer i

	univocal_nodes%numberID = nkn

	allocate(univocal_nodes%localID(nkn))

	do i=1,nkn
	  univocal_nodes%localID(i) = i
	end do

	end subroutine shympi_univocal_nodes

!******************************************************************

        subroutine check_part_basin(what)

        use basin

        implicit none

        character*(5) what
        integer pnkn,pnel,pn_threads,ierr
        integer control
        character*(14) filename
        character*(11) pwhat
        integer i

        call shympi_get_filename(filename,what)

        write(6,*) filename,what

        open(unit=108, file=filename, form="formatted"
     +   , iostat=control, status="old", action="read")
        
        if(control .ne. 0) then
          if(my_id .eq. 0) then
            write(6,*)'error stop: partitioning file not found'
          end if
          call shympi_barrier
        stop
        end if

        read(unit=108, fmt="(i12,i12,i12,A12)") pnkn,pnel,pn_threads
     +                  ,pwhat

        if(pnkn .ne. nkndi .or. pnel .ne. neldi .or. pn_threads
     &          .ne. n_threads .or. pwhat .ne. what) then
         if(my_id .eq. 0) then
          write(6,*)'basin file does not match'
          write(6,*)'partitioning file:nkn,nel,n_threads,partitioning'
     +          ,pnkn,pnel,pn_threads,pwhat
          write(6,*)'basin in str file:nkn,nel,n_threads,partitioning'
     +          ,nkndi,neldi,n_threads,what
         end if
         call shympi_barrier
         stop
        end if

        if(what .eq. 'nodes') then
          allocate(allPartAssign(nkndi))
          read(unit=108,fmt="(i12,i12,i12,i12,i12,i12)")
     +          (allPartAssign(i),i=1,nkndi)
        else 
          write(6,*)'error partitioning file on nodes'
          stop
        end if

        close(108)

        return

        end subroutine check_part_basin        

!******************************************************************

        subroutine shympi_get_filename(filename,what)

          implicit none

          character*(5) what
          character*(14) filename,format_string

          if(n_threads .gt. 999) then
            format_string = "(A5,A5,I4)"
            write(filename,format_string)'part_',what,n_threads
          else if(n_threads .gt. 99) then
            format_string = "(A5,A5,I3)"
            write(filename,format_string)'part_',what,n_threads
          else if(n_threads .gt. 9) then
            format_string = "(A5,A5,I2)"
            write(filename,format_string)'part_',what,n_threads
          else
            format_string = "(A5,A5,I1)"
            write(filename,format_string)'part_',what,n_threads
          end if

          return

        end subroutine shympi_get_filename


!==================================================================
        end module shympi
!==================================================================

