!
! mpi routines
!
! contents :
!
! revision log :
!
! 24.11.2015    ggu     project started
!
!******************************************************************

	subroutine shympi_init_internal(my_id,n_threads)

	implicit none

	integer ierr

	integer my_id,n_threads

	my_id = 0
	n_threads = 1

	end subroutine shympi_init_internal

!******************************************************************

	subroutine shympi_barrier_internal

	implicit none

	integer ierr

	end subroutine shympi_barrier_internal

!******************************************************************

	subroutine shympi_finalize_internal

	implicit none

	integer ierr

	end subroutine shympi_finalize_internal

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine shympi_exchange_internal_i(nlvddi,n,il,val)

	implicit none

	integer nlvddi,n
	integer il(n)
	integer val(nlvddi,n)

	end subroutine shympi_exchange_internal_i

!******************************************************************

	subroutine shympi_exchange_internal_r(nlvddi,n,il,val)
	
	integer nlvddi,n
	integer il(n)
	real val(nlvddi,n)

        stop 'error stop shympi_exchange_internal_r: not ready'

	end subroutine shympi_exchange_internal_r

!******************************************************************

	subroutine shympi_exchange_internal_d(nlvddi,n,il,val)
	
	integer nlvddi,n
	integer il(n)
	double precision val(nlvddi,n)

        stop 'error stop shympi_exchange_internal_d: not ready'

	end subroutine shympi_exchange_internal_d

!******************************************************************
!******************************************************************
!******************************************************************

