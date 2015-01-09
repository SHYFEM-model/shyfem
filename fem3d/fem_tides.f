
	module fem_tides

	implicit none

        !real xgeov(nkndim), ygeov(nkndim)
        !common /xgeov/xgeov, /ygeov/ygeov
        !real zeqv(nkndim)
        !common /zeqv/zeqv
	!save /xgeov/, /ygeov/, /zeqv/

	integer, private, save :: nkn_tides = 0

	real, allocatable, save :: xgeov(:)
	real, allocatable, save :: ygeov(:)
	real, allocatable, save :: zeqv(:)

	contains

!************************************************************

	subroutine fem_tides_init(nkn)

	integer nkn

	if( nkn == nkn_tides ) return

	if( nkn_tides > 0 ) then
	  deallocate(xgeov)
	  deallocate(ygeov)
	  deallocate(zeqv)
	end if

	if( nkn == 0 ) return

	nkn_tides = nkn

	allocate(xgeov(nkn))
	allocate(ygeov(nkn))
	allocate(zeqv(nkn))

	end subroutine fem_tides_init

!************************************************************

	end module fem_tides

