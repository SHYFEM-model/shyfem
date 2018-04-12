
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
