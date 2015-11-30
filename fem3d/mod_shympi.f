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

	integer,save :: nkn_global = 0
	integer,save :: nel_global = 0
	integer,save :: nkn_local = 0
	integer,save :: nel_local = 0
	integer,save :: nkn_inner = 0
	integer,save :: nel_inner = 0

	integer,save :: n_ghost_areas = 0
	integer,save :: n_ghost_nodes_max = 0
	integer,save :: n_ghost_elems_max = 0
	integer,save :: n_ghost_max = 0
	integer,save,allocatable :: ghost_areas(:,:)
	integer,save,allocatable :: ghost_nodes_in(:,:)
	integer,save,allocatable :: ghost_nodes_out(:,:)
	integer,save,allocatable :: ghost_elems(:,:)

	integer,save,allocatable :: i_buffer_in(:,:)
	integer,save,allocatable :: i_buffer_out(:,:)
	
	integer,save,allocatable :: node_area(:)	!global
	integer,save,allocatable :: request(:)		!for exchange
	integer,save,allocatable :: status(:,:)		!for exchange

	logical,save,allocatable :: is_inner_node(:)
	logical,save,allocatable :: is_inner_elem(:)
	integer,save,allocatable :: id_node(:)
	integer,save,allocatable :: id_elem(:,:)

	integer,save,allocatable :: ival(:)

        INTERFACE shympi_exchange_node
        MODULE PROCEDURE  shympi_exchange_node_r
     +                   ,shympi_exchange_node_d
     +                   ,shympi_exchange_node_i
        END INTERFACE

        INTERFACE shympi_exchange_2d_node
        MODULE PROCEDURE  shympi_exchange_2d_node_r
     +                   ,shympi_exchange_2d_node_d
     +                   ,shympi_exchange_2d_node_i
        END INTERFACE

        INTERFACE shympi_check_2d_node
        MODULE PROCEDURE  shympi_check_2d_node_r
     +                   ,shympi_check_2d_node_d
     +                   ,shympi_check_2d_node_i
        END INTERFACE

!==================================================================
        contains
!==================================================================

	subroutine shympi_init

	use basin

	integer ierr
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

	allocate(node_area(nkn_global))
	allocate(ival(n_threads))
	!allocate(request(2*n_threads))
	!allocate(status(MPI_STATUS_SIZE,2*n_threads))

	if( bmpi ) then
	  write(cunit,'(i10)') my_id
	  cunit = adjustl(cunit)
	  file = 'mpi_debug_' // trim(cunit) // '.txt'
	  open(newunit=my_unit,file=file,status='unknown')
	  write(my_unit,*) 'shympi initialized: ',my_id,n_threads
	end if

	write(6,*) 'shympi initialized: ',my_id,n_threads

	call shympi_barrier

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

	subroutine shympi_barrier

	integer ierr

	call shympi_barrier_internal

	end subroutine shympi_barrier

!******************************************************************

	subroutine shympi_finalize

	integer ierr

	call shympi_finalize_internal

	end subroutine shympi_finalize

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_node_r(val)

	use basin
	use levels

	real val(nlvdi,nkn)

	call shympi_exchange_internal_r(nlvdi,nkn,ilhkv,val)

	end subroutine shympi_exchange_node_r

!*******************************

	subroutine shympi_exchange_node_d(val)

	use basin
	use levels

	double precision val(nlvdi,nkn)

	call shympi_exchange_internal_d(nlvdi,nkn,ilhkv,val)

	end subroutine shympi_exchange_node_d

!*******************************

	subroutine shympi_exchange_node_i(val)

	use basin
	use levels

	integer val(nlvdi,nkn)

	call shympi_exchange_internal_i(nlvdi,nkn,ilhkv,val)

	end subroutine shympi_exchange_node_i

!******************************************************************

	subroutine shympi_exchange_2d_node_r(val)

	use basin
	use levels

	real val(nkn)

	call shympi_exchange_internal_r(1,nkn,ilhkv,val)

	end subroutine shympi_exchange_2d_node_r

!*******************************

	subroutine shympi_exchange_2d_node_d(val)

	use basin
	use levels

	double precision val(nkn)

	call shympi_exchange_internal_d(1,nkn,ilhkv,val)

	end subroutine shympi_exchange_2d_node_d

!*******************************

	subroutine shympi_exchange_2d_node_i(val)

	use basin
	use levels

	integer val(nkn)

	call shympi_exchange_internal_i(1,nkn,ilhkv,val)

	end subroutine shympi_exchange_2d_node_i

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_check_2d_node_r(val)

	use basin

	real val(nkn)

        stop 'error stop shympi_check_2d_node_r: not ready'

	end subroutine shympi_check_2d_node_r

!*******************************

	subroutine shympi_check_2d_node_d(val)

	use basin

	double precision val(nkn)

        stop 'error stop shympi_check_2d_node_d: not ready'

	end subroutine shympi_check_2d_node_d

!*******************************

	subroutine shympi_check_2d_node_i(val,text)

	use basin

	integer val(nkn)
	character*(*) text

	integer iaux(nkn)

	iaux = ival
	call shympi_exchange_2d_node_i(iaux)

        if( .not. all( val == iaux ) ) then
          write(6,*) 'arrays are different: ' // text
          stop 'error stop shympi_check_2d_node_i'
        end if

	end subroutine shympi_check_2d_node_i

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_gather_i(val)

	integer val

!	MPI_GATHER	-> gather into ival
	ival(1) = val

	end subroutine shympi_gather_i

!*******************************

	subroutine shympi_bcast_i(val)

	integer val

!	MPI_BCAST	-> broadcast value to all tasks

	end subroutine shympi_bcast_i

!*******************************

	subroutine shympi_reduce_r(what,vals,val)

	character*(*) what
	real vals(:)
	real val

!	MPI_REDUCE	-> reduce over all tasks

	if( what == 'min' ) then
	  val = MINVAL(vals)
	else if( what == 'max' ) then
	  val = MAXVAL(vals)
	else
	  write(6,*) 'what = ',what
	  stop 'error stop shympi_reduce_r: not ready'
	end if

	end subroutine shympi_reduce_r

!******************************************************************

	function shympi_output()

	logical shympi_output

	shympi_output = .true.

	end function shympi_output

!==================================================================
        end module shympi
!==================================================================

