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

	subroutine shympi_abort_internal

	implicit none

	integer ierr

	end subroutine shympi_abort_internal

!******************************************************************

	subroutine shympi_finalize_internal

	implicit none

	integer ierr

	end subroutine shympi_finalize_internal

!******************************************************************

        subroutine shympi_get_status_size_internal(size)

        implicit none

        integer size

        size = 1

        end subroutine shympi_get_status_size_internal

!******************************************************************

	subroutine shympi_syncronize_internal

	implicit none

	integer ierr

	end subroutine shympi_syncronize_internal

!******************************************************************

	subroutine shympi_syncronize_initial

	implicit none

	integer ierr

	end subroutine shympi_syncronize_initial

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
	
	implicit none

	integer nlvddi,n
	integer il(n)
	real val(nlvddi,n)

	end subroutine shympi_exchange_internal_r

!******************************************************************

	subroutine shympi_exchange_internal_d(nlvddi,n,il,val)
	
	implicit none

	integer nlvddi,n
	integer il(n)
	double precision val(nlvddi,n)

	end subroutine shympi_exchange_internal_d

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_gather_i_internal(val)

        use shympi

        implicit none

        integer val

	ival(1) = val

        end subroutine shympi_gather_i_internal

!******************************************************************

        subroutine shympi_bcast_i_internal(val)

        use shympi

        implicit none

        integer val

        end subroutine shympi_bcast_i_internal

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine shympi_reduce_r_internal(what,val)
        
        implicit none

        character*(*) what
        real val

        end subroutine shympi_reduce_r_internal

!******************************************************************

