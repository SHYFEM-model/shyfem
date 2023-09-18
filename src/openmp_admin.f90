!---------------------------------------------------------------
        module openmp_admin
!---------------------------------------------------------------
        contains
!---------------------------------------------------------------
        
!***************************************************************

	subroutine openmp_get_max_threads(n)

	implicit none

	integer n

	n = 1

	end

!***************************************************************

	subroutine openmp_set_num_threads(n)

	implicit none

	integer n

	end

!***************************************************************

	subroutine openmp_get_num_threads(n)

	implicit none

	integer n

	n = 1

	end

!***************************************************************

	subroutine openmp_get_thread_num(it)

	implicit none

	integer it

	it = 0

	end

!***************************************************************

	function openmp_is_parallel()

	implicit none

	logical openmp_is_parallel

	openmp_is_parallel = .false.

	end

!***************************************************************

        function openmp_in_parallel()

! true if in parallel region

        implicit none

        logical openmp_in_parallel

        openmp_in_parallel = .false.

        end

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine get_clock_count(count)

        implicit none

        integer count
        integer count_rate,count_max

        call system_clock(count,count_rate,count_max)

        end

!***************************************************************

        subroutine get_clock_count_diff(count0,dcount)

        implicit none

        integer count0,dcount
        integer count,count_rate,count_max

        call system_clock(count,count_rate,count_max)

        dcount = count - count0

        end

!***************************************************************

        function openmp_get_wtime()

! gets time

        implicit none

        double precision openmp_get_wtime

        openmp_get_wtime = 0.

        end

!***************************************************************

	subroutine setup_omp_parallel

        use para

	implicit none

	integer n,nomp

	write(6,*) 'start of setup of parallel OMP threads'

	if( openmp_is_parallel() ) then
	  write(6,*) 'the program can run in OMP parallel mode'
	else
	  write(6,*) 'the program can run only in serial mode'
	end if
	  
	call openmp_get_max_threads(n)
	write(6,*) 'maximum available threads: ',n

	nomp = nint(getpar('nomp'))
	if( nomp .gt. 0 ) then
	  nomp = min(nomp,n)
	  call openmp_set_num_threads(nomp)
	else
	  nomp = n
	end if
	call putpar('nomp',dble(nomp))

	write(6,*) 'available threads: ',nomp

	write(6,*) 'end of setup of parallel OMP threads'

	end

!**********************************************************************
!---------------------------------------------------------------
        end module openmp_admin
!---------------------------------------------------------------
