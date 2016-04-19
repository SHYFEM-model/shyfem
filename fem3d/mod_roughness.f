	module mod_roughness

	implicit none

	integer, private, save :: nkn_roughness = 0

	real, allocatable, save :: z0bk(:)	!bottom roughenss on nodes
	real, allocatable, save :: z0sk(:)	!surface roughenss on nodes

	contains

!************************************************************

	subroutine mod_roughness_init(nkn)

	integer nkn

	if( nkn == nkn_roughness ) return

	if( nkn_roughness > 0 ) then
	  deallocate(z0bk)
	  deallocate(z0sk)
	end if

	nkn_roughness = nkn

	if( nkn == 0 ) return

	allocate(z0bk(nkn))
	allocate(z0sk(nkn))

	end subroutine mod_roughness_init

!************************************************************

	end module mod_roughness
