
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

! Fortran MPI example  

! revision log :
!
! 26.04.2018	ggu	changed VERS_7_5_46
! 14.02.2019	ggu	changed VERS_7_5_56

	program hello_world

	use mpi

	implicit none

	!include 'mpif.h'

	integer rank, size, ierror, tag, status(MPI_STATUS_SIZE)
	double precision timer
   
	!time = MPI_Wtime();

	call MPI_INIT(ierror)

	call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
	print*, 'Hello World (MPI): ',rank,size
	call MPI_BARRIER(MPI_COMM_WORLD, ierror)

	timer = MPI_Wtime();
	call do_work(rank,size)
	timer = MPI_Wtime() - timer;
	print*, 'timing: ',rank,timer
	call MPI_BARRIER(MPI_COMM_WORLD, ierror)

	call MPI_FINALIZE(ierror)

	end

!--------------------------------------------------------------------------

	subroutine do_work(rank,size)

	implicit none

	integer :: rank,size

        integer :: i
        integer, parameter :: nloop = 1000000
        double precision :: val
        double precision :: array(nloop)

        do i=1+rank,nloop,size
          val = i
          array(i) = log(val)
        end do

	end subroutine

!--------------------------------------------------------------------------

