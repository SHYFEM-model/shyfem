!
! handles timing information
!
!***************************************************************

!===============================================================
	module shympi_time
!===============================================================

	implicit none

	integer, parameter :: ndim = 10

	double precision, save :: time(ndim) = 0.

!===============================================================
	end module shympi_time
!===============================================================

	subroutine shympi_time_reset(id)

	use shympi_time

	implicit none

	integer id

	if( id < 1 .or. id > ndim ) call shympi_time_error(id)

	time(id) = 0.

	end

!***************************************************************

	subroutine shympi_time_accum(id,dt)

	use shympi_time

	implicit none

	integer id
	double precision dt

	if( id < 1 .or. id > ndim ) call shympi_time_error(id)

	time(id) = time(id) + dt

	end

!***************************************************************

	subroutine shympi_time_get(id,at)

	use shympi_time

	implicit none

	integer id
	double precision at

	if( id < 1 .or. id > ndim ) call shympi_time_error(id)

	at = time(id)

	end

!***************************************************************

	subroutine shympi_time_error(id)

	use shympi_time

	implicit none

	integer id

	write(6,*) 'id = ',id,'  ndim = ',ndim
	stop 'error stop shympi_time_error: id out of range'

	end

!***************************************************************

