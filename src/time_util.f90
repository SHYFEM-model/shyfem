!----------------------------------------------------------------------------------
        module time_util
!----------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
        contains
!----------------------------------------------------------------------------------

!**********************************************************************

        subroutine get_timestep(dt)

! returns double precision time step (in real seconds)

        implicit none

	include 'femtime.h'

	double precision dt		!time step (return)

	!dt = idt/float(itunit)
	dt = dt_act

	end

!**********************************************************************

        subroutine get_act_time(itact)

! returns actual time

        implicit none

	include 'femtime.h'

	integer itact

	itact = it

	end

!**********************************************************************

        subroutine get_orig_timestep(dt)

! returns original double precision time step (in real seconds)

        implicit none

	double precision dt		!time step (return)

	include 'femtime.h'

	dt = idtorig/float(itunit)

	end

!**********************************************************************

!----------------------------------------------------------------------------------
        end module time_util
!----------------------------------------------------------------------------------
