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

        module shympi_aux

        use mpi

        implicit none

        !include 'mpif.h'

        end module shympi_aux

!********************************

	subroutine shympi_init_internal(my_id,n_threads)

	use shympi_aux

	implicit none

	integer my_id,n_threads

	integer ierr

        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, n_threads, ierr )

	end subroutine shympi_init_internal

!******************************************************************

	subroutine shympi_barrier_internal

	use shympi_aux

	implicit none

	integer ierr

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	end subroutine shympi_barrier_internal

!******************************************************************

        subroutine shympi_abort_internal

	use shympi_aux

        implicit none

        integer ierr,ierr_code

        ierr_code = 33
	call MPI_ABORT(MPI_COMM_WORLD,ierr_code,ierr)

        end subroutine shympi_abort_internal

!******************************************************************

        subroutine shympi_finalize_internal

	use shympi_aux

        implicit none

        integer ierr

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)
	call MPI_FINALIZE(ierr)

        end subroutine shympi_finalize_internal

!******************************************************************

        subroutine shympi_get_status_size_internal(size)

	use shympi_aux

        implicit none

        integer size

	size = MPI_STATUS_SIZE

        end subroutine shympi_get_status_size_internal

!******************************************************************

	subroutine shympi_syncronize_internal

	use shympi_aux

	implicit none

	integer ierr

	flush(6)
	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	end subroutine shympi_syncronize_internal

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_syncronize_initial

	use shympi_aux

        implicit none

	integer my_id,nt
	integer root,ierr,i
	integer count
	integer local
	integer, allocatable :: buf(:)

        call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, nt, ierr )

	allocate(buf(nt))

	root = my_id
	count = 1
	root = 0

	if( my_id == root ) then
	  do i=1,nt
	    buf(i) = i
	  end do
	end if

	call MPI_SCATTER (buf,count,MPI_INT
     +			,local,count,MPI_INT
     +			,root,MPI_COMM_WORLD,ierr)

	local = local * 2
	root = 0

	call MPI_GATHER (local,count,MPI_INT
     +			,buf,count,MPI_INT
     +			,root,MPI_COMM_WORLD,ierr)

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	if( my_id == root ) then
	  !write(6,*) 'mpi sync: ',nt,root,(buf(i),i=1,nt)
	  !write(6,*) 'mpi sync: ',nt,root
	end if

	deallocate(buf)

	end subroutine shympi_syncronize_initial

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_internal_i(belem,n0,nlvddi,n,il
     +						,g_in,g_out,val)

	use shympi_aux

	use shympi

	implicit none

	logical belem
	integer n0,nlvddi,n
	integer il(n)
	integer g_in(n_ghost_max,n_ghost_areas)
	integer g_out(n_ghost_max,n_ghost_areas)
	integer val(n0:nlvddi,n)

	integer tag,ir,ia,id
	integer i,k,nc,ierr
	integer nb
	integer iout,iin

        tag=1234
	ir = 0

	iout = 2
	iin = 3
	if( belem ) then
	  iout = 4
	  iin = 4
	end if

	nb = (nlvddi-n0+1) * n_ghost_max
	call shympi_alloc_buffer(nb)

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iout,ia)
	  call count_buffer(n0,nlvddi,n,nc,il,g_out(:,ia),nb)
	  !write(6,*) 'ex1: ',my_id,ia,id,nc,n,nb
          call MPI_Irecv(i_buffer_out(:,ia),nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iin,ia)
	  call to_buffer_i(n0,nlvddi,n,nc,il
     +		,g_in(:,ia),val,nb,i_buffer_in(:,ia))
          call MPI_Isend(i_buffer_in(:,ia),nb,MPI_INTEGER,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

        call MPI_WaitAll(ir,request,status,ierr)

	do ia=1,n_ghost_areas
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iout,ia)
	  call from_buffer_i(n0,nlvddi,n,nc,il
     +		,g_out(:,ia),val,nb,i_buffer_out(:,ia))
	end do

	end subroutine shympi_exchange_internal_i

!******************************************************************

	subroutine shympi_exchange_internal_r(belem,n0,nlvddi,n,il
     +						,g_in,g_out,val)

	use shympi_aux

	use shympi

	implicit none

	logical belem
	integer n0,nlvddi,n
	integer il(n)
	integer g_in(n_ghost_max,n_ghost_areas)
	integer g_out(n_ghost_max,n_ghost_areas)
	real val(n0:nlvddi,n)

	integer tag,ir,ia,id
	integer i,k,nc,ierr
	integer nb
	integer iout,iin

        tag=1234
	ir = 0

	iout = 2
	iin = 3
	if( belem ) then
	  iout = 4
	  iin = 4
	end if

	nb = (nlvddi-n0+1) * n_ghost_max
	call shympi_alloc_buffer(nb)

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iout,ia)
	  call count_buffer(n0,nlvddi,n,nc,il,g_out(:,ia),nb)
	  !write(6,*) 'ex1: ',my_id,ia,id,nc,n,nb
          call MPI_Irecv(r_buffer_out(:,ia),nb,MPI_REAL,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

	do ia=1,n_ghost_areas
	  ir = ir + 1
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iin,ia)
	  call to_buffer_r(n0,nlvddi,n,nc,il
     +		,g_in(:,ia),val,nb,r_buffer_in(:,ia))
          call MPI_Isend(r_buffer_in(:,ia),nb,MPI_REAL,id
     +	          ,tag,MPI_COMM_WORLD,request(ir),ierr)
	end do

        call MPI_WaitAll(ir,request,status,ierr)

	do ia=1,n_ghost_areas
	  id = ghost_areas(1,ia)
	  nc = ghost_areas(iout,ia)
	  call from_buffer_r(n0,nlvddi,n,nc,il
     +		,g_out(:,ia),val,nb,r_buffer_out(:,ia))
	end do

	end subroutine shympi_exchange_internal_r

!******************************************************************

	subroutine shympi_exchange_internal_d(belem,n0,nlvddi,n,il
     +						,g_in,g_out,val)

	use shympi_aux
	use shympi
	
	logical belem
	integer n0,nlvddi,n
	integer il(n)
	integer g_in(n_ghost_max,n_ghost_areas)
	integer g_out(n_ghost_max,n_ghost_areas)
	double precision val(n0:nlvddi,n)

        stop 'error stop shympi_exchange_internal_d: not ready'

	end subroutine shympi_exchange_internal_d

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_gather_i_internal(val)

	use shympi_aux

	use shympi

	implicit none

        integer val

        integer ierr

        call MPI_GATHER (val,1,MPI_INT
     +                  ,ival,1,MPI_INT
     +                  ,0,MPI_COMM_WORLD,ierr)

        end subroutine shympi_gather_i_internal

!*******************************

        subroutine shympi_bcast_i_internal(val)

	use shympi_aux

	implicit none

        integer val

        integer ierr

        call MPI_BCAST(val,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

        end subroutine shympi_bcast_i_internal

!*******************************

	subroutine shympi_reduce_r_internal(what,val)

	use shympi_aux

	implicit none

	character*(*) what
	real val

        integer ierr
	real valout

        if( what == 'min' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_REAL,MPI_MIN
     +				,MPI_COMM_WORLD,ierr)
	  val = valout
        else if( what == 'max' ) then
	  call MPI_ALLREDUCE(val,valout,1,MPI_REAL,MPI_MAX
     +				,MPI_COMM_WORLD,ierr)
	  val = valout
        else
          write(6,*) 'what = ',what
          stop 'error stop shympi_reduce_r_internal: not ready'
        end if

	end subroutine shympi_reduce_r_internal

!******************************************************************
!******************************************************************
!******************************************************************

