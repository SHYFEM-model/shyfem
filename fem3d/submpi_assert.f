!
! checks if all assumptions on variables are true
!
!***************************************************************

	subroutine mpi_assert_all

	call mpi_assert_coriolis

	end

!***************************************************************

	subroutine mpi_assert_coriolis

	use shympi

	implicit none

	integer isphe

	call get_coords_ev(isphe)

	call shympi_gather_i(isphe)

	if( shympi_is_master() ) then
	  if( any(ival/=isphe) ) then
	    write(6,*) 'error in isphe: ',isphe,ival
	    stop 'error stop mpi_assert_coriolis: isphe'
	  end if
	end if

	end

!***************************************************************

