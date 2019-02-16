
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

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

        function openmp_is_master()

        implicit none
    
        logical openmp_is_master

        openmp_is_master = .true.

        end

c***************************************************************

	function openmp_is_parallel()

	implicit none

	logical openmp_is_parallel

	openmp_is_parallel = .false.

	end

c***************************************************************

        function openmp_in_parallel()

c true if in parallel region

        implicit none

        logical openmp_in_parallel

        openmp_in_parallel = .false.

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

        openmp_get_wtime = 0.

        end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine omp_compute_chunk(imax,nchunk)

        implicit none

        integer imax,nchunk

        nchunk = 1

        end

c***************************************************************

        subroutine omp_compute_minmax(nchunk,imax,i,iend)

        implicit none

        integer nchunk,imax,i,iend

        iend = i

        end

c***************************************************************

