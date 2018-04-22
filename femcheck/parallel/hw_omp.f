
	Program Hello_from_Threads

	use omp_lib

	implicit none

	integer :: iam
	integer :: nthreads = 3

	call omp_set_num_threads(nthreads) 

!$omp parallel private(iam)
	iam=omp_get_thread_num()

!$omp critical
	write(*,*) 'Hello World (OMP): ', iam,nthreads
!$omp end critical

!$omp end parallel

	end program Hello_from_Threads

