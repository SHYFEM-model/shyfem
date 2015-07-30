
	module mod_adj_grade

	integer, save :: ngrdi = 0

	integer, save, allocatable :: ngrade(:)
	integer, save, allocatable :: nbound(:)
	integer, save, allocatable :: ngri(:,:)

	contains

	subroutine mod_adj_grade_init(nkn,ngr)

	integer nkn,ngr

	ngrdi = ngr

	allocate(ngrade(nkn))
	allocate(nbound(nkn))
	allocate(ngri(2*ngr,nkn))

	end subroutine mod_adj_grade_init

	end module mod_adj_grade
