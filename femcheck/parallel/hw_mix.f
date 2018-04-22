
! Fortran MPI example  

	program hello_world

	use mpi
        use omp_lib

	implicit none

	!include 'mpif.h'

	integer :: rank, size, ierror, tag, status(MPI_STATUS_SIZE)
        integer :: iam
        integer :: nthreads = 3

        call omp_set_num_threads(nthreads)
   
	call MPI_INIT(ierror)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

!$omp parallel private(iam)
        iam=omp_get_thread_num()
!$omp critical
        !write(*,*) 'Hello World (OMP): ', iam,nthreads
	write(6,'(a,a,i4,a,i4)') 'Hello World (MPI/OMP): '
     +			,'   mpi_rank = ',rank
     +			,'   omp_thread = ',iam
!$omp end critical
!$omp end parallel

	call MPI_FINALIZE(ierror)

	end





