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

	subroutine shympi_init_internal(my_id,n_threads)

	use shympi

	implicit none

	include "mpif.h"

	integer ierr

        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, n_threads, ierr )

	end subroutine shympi_init_internal

!******************************************************************

	subroutine shympi_barrier_internal

	implicit none

	include "mpif.h"

	integer ierr

	call MPI_BARRIER( MPI_COMM_WORLD, ierr)

	end subroutine shympi_barrier_internal

!******************************************************************

        subroutine shympi_finalize_internal

        implicit none

	include "mpif.h"

        integer ierr

	call MPI_FINALIZE(ierr)

        end subroutine shympi_finalize_internal

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_internal_i(nlvddi,n,il,val)

	use shympi

	implicit none

	include "mpif.h"

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


