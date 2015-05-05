
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

        function openmp_is_parallel()

        implicit none

        logical openmp_is_parallel

	openmp_is_parallel = .true.

        end

c***************************************************************

        function openmp_in_parallel()

c true if in parallel region

        implicit none

        logical openmp_in_parallel
	logical OMP_IN_PARALLEL

	openmp_in_parallel = OMP_IN_PARALLEL()

        end

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine get_clock_count(count)

	implicit none

	integer count
	integer count_rate,count_max

	call system_clock(count,count_rate,count_max)

	end

c***************************************************************

	subroutine get_clock_count_diff(count0,dcount)

	implicit none

	integer count0,dcount
	integer count,count_rate,count_max

	call system_clock(count,count_rate,count_max)

	dcount = count - count0

	end

c***************************************************************

        function openmp_get_wtime()

c gets time

        implicit none

        double precision openmp_get_wtime
	double precision OMP_GET_WTIME

	openmp_get_wtime = OMP_GET_WTIME()
	
	end

c***************************************************************

