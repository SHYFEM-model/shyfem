
!--------------------------------------------------------------------------
!
!    Copyright (C) 2012,2015-2017,2019  Georg Umgiesser
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

! revision log :
!
! 25.01.2012	ggu	changed VERS_6_1_42
! 30.03.2012	ggu	changed VERS_6_1_51
! 29.08.2012	ggu	changed VERS_6_1_56
! 26.02.2015	ggu	changed VERS_7_1_5
! 05.05.2015	ggu	changed VERS_7_1_10
! 09.09.2016	ggu	changed VERS_7_5_17
! 12.01.2017	ggu	changed VERS_7_5_21
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

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

	subroutine openmp_parallel_code(text)

	implicit none

	character*(*) text

	text = 'serial'

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

