
	module mod_tides

	implicit none

	integer, private, save :: nkn_tides = 0

	real, allocatable, save :: xgeov(:)
	real, allocatable, save :: ygeov(:)
	real, allocatable, save :: xcartv(:)
	real, allocatable, save :: ycartv(:)
	real, allocatable, save :: zeqv(:)

	contains

!************************************************************

	subroutine mod_tides_init(nkn)

	integer nkn

	if( nkn == nkn_tides ) return

	if( nkn_tides > 0 ) then
	  deallocate(xgeov)
	  deallocate(ygeov)
	  deallocate(xcartv)
	  deallocate(ycartv)
	  deallocate(zeqv)
	end if

	nkn_tides = nkn

	if( nkn == 0 ) return

	allocate(xgeov(nkn))
	allocate(ygeov(nkn))
	allocate(xcartv(nkn))
	allocate(ycartv(nkn))
	allocate(zeqv(nkn))

	end subroutine mod_tides_init

!************************************************************

	end module mod_tides

