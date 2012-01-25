
c***************************************************************

	subroutine openmp_get_max_threads(n)

	implicit none

	integer n
	integer OMP_GET_MAX_THREADS

	n = OMP_GET_MAX_THREADS()

	end

c***************************************************************

	subroutine openmp_set_num_threads(n)

	implicit none

	integer n

	call OMP_SET_NUM_THREADS(n)

	end

c***************************************************************

	subroutine openmp_get_num_threads(n)

	implicit none

	integer n
	integer OMP_GET_NUM_THREADS

	n = OMP_GET_NUM_THREADS()

	end

c***************************************************************

	subroutine openmp_get_thread_num(it)

	implicit none

	integer it
	integer OMP_GET_THREAD_NUM

	it = OMP_GET_THREAD_NUM()

	end

c***************************************************************

