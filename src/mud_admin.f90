
!*******************************************************************
!
! if the mud module is not needed we do the minimum needed
!
! we have to make sure that vts is initialized -> submud.f
!
! questions: 
!	do we really need an extra array vts?
!	can we not use the original array visv for this?
!	what about diffusivity?
!
!*******************************************************************
!-------------------------------------------------------------------
        module mud_admin
!-------------------------------------------------------------------
        contains
!-------------------------------------------------------------------

	subroutine readmud

	use fluidmud
	use basin
	use levels
        use nls

	implicit none

	call nrdskp

	end

!*******************************************************************

	subroutine submud(it,dt)

        implicit none
        integer it
        double precision dt

	end

!*******************************************************************

	subroutine submud_init

	use fluidmud
	use basin
	use levels

	implicit none

	call mod_fluidmud_dummy_init(nkn,nlv)
	vts = 0.

	end

!*******************************************************************

	subroutine set_mud_roughness(k,l,alpha)

	implicit none

	integer k,l
	double precision alpha

	alpha = 1.

	end

!*******************************************************************

	subroutine set_rhomud(k,l,rhop)

	implicit none

	integer k,l
	double precision rhop

	end

!*******************************************************************

!-------------------------------------------------------------------
        end module mud_admin
!-------------------------------------------------------------------
