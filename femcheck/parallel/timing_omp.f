
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
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
! 26.04.2018	ggu	changed VERS_7_5_46
! 14.02.2019	ggu	changed VERS_7_5_56

	Program Timing

	call do_work(1)
	call do_work(3)
	call do_work(5)

	end program

!--------------------------------------------------------------------------

	subroutine do_work(nthreads)

	use omp_lib

	implicit none

	integer :: nthreads

	integer :: iam
	integer :: i
	integer, parameter :: nloop = 1000000
	double precision :: timer
	double precision :: val
	double precision :: array(nloop)

	call omp_set_num_threads(nthreads) 

!$omp parallel private(iam,timer,i,val)
	iam=omp_get_thread_num()
!$omp critical
	write(*,*) 'Hello World (OMP): ', iam,nthreads
!$omp end critical
!$omp barrier

!$      timer = omp_get_wtime() 

!$omp do
	do i=1,nloop
	  val = i
	  array(i) = log(val)
	end do
!$omp end do

!$      timer = omp_get_wtime() - timer
!$      print *,"TIME TO SOLUTION (OMP)  = ",iam,timer
!$omp barrier

!$omp end parallel

	end subroutine

!--------------------------------------------------------------------------

