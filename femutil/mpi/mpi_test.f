
	program mpi_test

	implicit none

	integer myid,numprocs,ierr

!        use fmpi
        include "mpif.h"

        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

        write (*,*) "Hello from ",myid,' / ',numprocs

        call MPI_FINALIZE(ierr)

	end
