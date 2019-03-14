
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2018  Georg Umgiesser
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
! 05.12.2017    ggu     changed VERS_7_5_39
! 19.04.2018    ggu     changed VERS_7_5_45
! 18.12.2018    ggu     changed VERS_7_5_52

!--------------------------------------------------------------------------

	subroutine mpi_init(ierr)
	integer ierr
	ierr = 0
	end

	subroutine mpi_comm_rank(comm,my_id,ierr)
	integer comm,my_id,ierr
	my_id = 0
	ierr = 0
	end

	subroutine mpi_comm_size(comm,n_threads,ierr)
	integer comm,n_threads,ierr
	n_threads = 1
	ierr = 0
	end

	subroutine mpi_barrier(comm,ierr)
	integer comm,ierr
	ierr = 0
	end

	subroutine mpi_finalize
	end

	subroutine mpi_irecv
	end

	subroutine mpi_isend
	end

	subroutine mpi_waitall
	end

	subroutine mpi_abort
	end

	subroutine mpi_scatter
	end

	subroutine mpi_gather
	end

	subroutine mpi_allgather
	end

	subroutine mpi_bcast
	end

	subroutine mpi_allreduce(what,val)
	end

	subroutine mpi_wtime
	end
