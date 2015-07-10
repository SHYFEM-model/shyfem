	module mod_roughness

	implicit none

        !real z0bk(nkndim)                   !bottom roughenss on nodes
        !common /z0bk/z0bk
        !real z0sk(nkndim)                   !surface roughenss on nodes
        !common /z0sk/z0sk
        !save /z0bk/,/z0sk/

	integer, private, save :: nkn_roughness = 0

	real, allocatable, save :: z0bk(:)
	real, allocatable, save :: z0sk(:)

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
