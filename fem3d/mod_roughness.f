
!============================================================
	module mod_roughness
!============================================================

	implicit none

	integer, private, save :: nkn_roughness = 0

	real, allocatable, save :: z0bk(:)	!bottom roughness on nodes
	real, allocatable, save :: z0sk(:)	!surface roughness on nodes

	real, parameter :: z0bmin = 1.e-4
	real, parameter :: z0smin = 0.02

	real, parameter :: z0bini = 0.03*0.03
	real, parameter :: z0sini = 0.02

!============================================================
	contains
!============================================================

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

	z0bk = z0bini
	z0sk = z0sini

	end subroutine mod_roughness_init

!============================================================
	end module mod_roughness
!============================================================

