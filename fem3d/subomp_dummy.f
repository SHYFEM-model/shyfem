
c***************************************************************

	subroutine openmp_get_max_threads(n)

	implicit none

	integer n

	n = 1

	end

c***************************************************************

	subroutine openmp_set_num_threads(n)

	implicit none

	integer n

	end

c***************************************************************

	subroutine openmp_get_num_threads(n)

	implicit none

	integer n

	n = 1

	end

c***************************************************************

	subroutine openmp_get_thread_num(it)

	implicit none

	integer it

	it = 0

	end

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

