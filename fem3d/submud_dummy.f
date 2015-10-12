
c*******************************************************************
c
c if the mud module is not needed we do the minimum needed
c
c we have to make sure that vts is initialized -> submud.f
c
c questions: 
c	do we really need an extra array vts?
c	can we not use the original array visv for this?
c	what about diffusivity?
c
c*******************************************************************

	subroutine readmud

	use mod_fluidmud
	use basin
	use levels

	implicit none

	call nrdskp

	end

c*******************************************************************

	subroutine submud(it,dt)
	end

c*******************************************************************

	subroutine submud_init

	use mod_fluidmud
	use basin
	use levels

	implicit none

	call mod_fluidmud_dummy_init(nkn,nlv)
	vts = 0.

	end

c*******************************************************************

	subroutine set_mud_roughness(k,l,alpha)

	implicit none

	integer k,l
	real alpha

	alpha = 1.

	end

c*******************************************************************

	subroutine set_rhomud(k,l,rhop)

	implicit none

	integer k,l
	real rhop

	end

c*******************************************************************

