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

	include "mpif.h"

	integer,save :: n_threads = 0
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

        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, n_threads, ierr )

	nkn_global = nkn
	nel_global = nel

	allocate(ival(n_threads))
	allocate(request(2*n_threads))
	allocate(status(MPI_STATUS_SIZE,2*n_threads))
	allocate(node_area(nkn_global))

	!allocate(is_inner_node(nkn))
	!allocate(is_inner_elem(nel))
	!allocate(id_node(nkn))
	!allocate(id_elem(2,nel))

	!is_inner_node = .true.
	!is_inner_elem = .true.
	!id_node = my_id
	!id_elem = my_id

	write(cunit,'(i10)') my_id
	cunit = adjustl(cunit)
	file = 'mpi_debug_' // trim(cunit) // '.txt'

	open(newunit=my_unit,file=file,status='unknown')

	write(6,*) 'shympi initialized: ',my_id,n_threads
	write(my_unit,*) 'shympi initialized: ',my_id,n_threads

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	end subroutine shympi_init

!******************************************************************

	subroutine shympi_alloc

	use basin

	allocate(is_inner_node(nkn))
	allocate(is_inner_elem(nel))
	allocate(id_node(nkn))
	allocate(id_elem(2,nel))

	end subroutine shympi_alloc

!******************************************************************

	subroutine shympi_finalize

	integer ierr

	call MPI_FINALIZE(ierr)

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

	subroutine shympi_exchange_internal_i(nlvddi,n,il,val)

	implicit none

	integer nlvddi,n
	integer il(n)
	integer val(nlvddi,n)

	integer tag,ir,ia,id
	integer i,k,nc,ierr
	integer nb

        tag=1234
	ir = 0

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(2,ia)
	  call count_buffer(nlvddi,n,nc,il,ghost_nodes_out,nb)
          call MPI_Irecv(i_buffer_out,nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(3,ia)
	  call to_buffer_i(nlvddi,n,nc,il
     +		,ghost_nodes_in(:,ia),val,nb,i_buffer_in(:,ia))
          call MPI_Isend(i_buffer_in,nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

        call MPI_WaitAll(ir,request,status,ierr)

	do ia=1,n_ghost_areas
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(2,ia)
	  call from_buffer_i(nlvddi,n,nc,il
     +		,ghost_nodes_out(:,ia),val,nb,i_buffer_out(:,ia))
	end do

	end subroutine shympi_exchange_internal_i

!******************************************************************

	subroutine shympi_exchange_internal_r(nlvddi,n,il,val)
	
	integer nlvddi,n
	integer il(n)
	real val(nlvddi,n)

        stop 'error stop shympi_exchange_internal_r: not ready'

	end subroutine shympi_exchange_internal_r

!******************************************************************

	subroutine shympi_exchange_internal_d(nlvddi,n,il,val)
	
	integer nlvddi,n
	integer il(n)
	double precision val(nlvddi,n)

        stop 'error stop shympi_exchange_internal_d: not ready'

	end subroutine shympi_exchange_internal_d

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

	subroutine to_buffer_i(nlvddi,n,nc,il,nodes,val,nb,i_buffer)

	integer nlvddi,n,nc
	integer il(n)
	integer nodes(nc)
	integer val(nlvddi,n)
	integer nb
	integer i_buffer(:)

	integer i,k,l,lmax

	if( nlvddi == 1 ) then
	  do i=1,nc
	    k = nodes(i)
	    i_buffer(i) = val(1,k)
	  end do
	  nb = nc
	else
	  nb = 0
	  do i=1,nc
	    k = nodes(i)
	    lmax = il(k)
	    do l=1,lmax
	      nb = nb + 1
	      i_buffer(nb) = val(l,k)
	    end do
	  end do
	end if

	end subroutine to_buffer_i

!******************************************************************

	subroutine from_buffer_i(nlvddi,n,nc,il,nodes,val,nb,i_buffer)

	integer nlvddi,n,nc
	integer il(n)
	integer nodes(nc)
	integer val(nlvddi,n)
	integer nb
	integer i_buffer(:)

	integer i,k,l,lmax

	if( nlvddi == 1 ) then
	  do i=1,nc
	    k = nodes(i)
	    val(1,k) = i_buffer(i)
	  end do
	  nb = nc
	else
	  nb = 0
	  do i=1,nc
	    k = nodes(i)
	    lmax = il(k)
	    do l=1,lmax
	      nb = nb + 1
	      val(l,k) = i_buffer(nb)
	    end do
	  end do
	end if

	end subroutine from_buffer_i

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

