
	subroutine nonhydro_solve_matrix

! gaussian solver
	
	implicit none
	
	write(6,*) '*** cannot use Gauss solver for nonhydro'
	stop 'error stop nonhydro_solve_matrix: no gauss'

	end
