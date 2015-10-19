
!==================================================================
        module levels
!==================================================================

        implicit none

        integer, save :: nlv = 0
        integer, save :: nlvdi = 0

        integer, save, private :: nlv_alloc = 0

        integer, save, allocatable :: ilhv(:)
        integer, save, allocatable :: ilhkv(:)
        integer, save, allocatable :: ilmv(:)
        integer, save, allocatable :: ilmkv(:)

        real, save, allocatable :: hlv(:)
        real, save, allocatable :: hldv(:)

!==================================================================
        contains
!==================================================================

	subroutine levels_init(nkn,nel,nl)

	integer nkn,nel,nl

	if( nl == 0 ) stop 'error stop levels_init: nl == 0'

	nlvdi = nl
	nlv_alloc = nl

	allocate(ilhv(nel))
	allocate(ilmv(nel))
	allocate(ilhkv(nkn))
	allocate(ilmkv(nkn))

	allocate(hlv(nl))
	allocate(hldv(nl))
	
	hlv = 0.
	hldv = 0.

	!write(6,*) 'levels allocated: ',nkn,nel,nl

	end subroutine levels_init

!******************************************************************

	subroutine levels_hlv_init(nl)

! allocates only hlv

	integer nl

	if( nlv_alloc == nl ) return

        if( nlv_alloc > 0 ) then
          deallocate(hlv)
        end if

	nlvdi = nl
        nlv_alloc = nl

        if( nl == 0 ) return

        allocate(hlv(nl))

	end subroutine levels_hlv_init

!******************************************************************

	subroutine levels_reinit(nl)

! re-allocates arrays depending on nl

	integer nl

	integer n
	real, allocatable :: hlv_aux(:)
	real, allocatable :: hldv_aux(:)

	if( nlv_alloc == nl ) return

	n = min(nl,nlv_alloc)
	write(6,*) 'levels_reinit: ',nl,nlv_alloc,n

	if( nlv_alloc > 0 ) then
	  if( nl > 0 ) then
	    allocate(hlv_aux(nl))
	    allocate(hldv_aux(nl))
	    hlv_aux = 0.
	    hldv_aux = 0.
	    hlv_aux(1:n) = hlv(1:n)
	    hldv_aux(1:n) = hldv(1:n)
	  end if
	  deallocate(hlv)
	  deallocate(hldv)
	end if

	write(6,*) 'levels_reinit: ',nl,nlv_alloc,n
	nlvdi = nl
	nlv_alloc = nl

	if( nl > 0 ) then
	  allocate(hlv(nl))
	  allocate(hldv(nl))
	  hlv=0.
	  hldv=0.
	  if( n > 0 ) then
	    hlv(1:n) = hlv_aux(1:n)
	    hldv(1:n) = hldv_aux(1:n)
	    deallocate(hlv_aux)
	    deallocate(hldv_aux)
	  end if
	end if
	
	end subroutine levels_reinit

!******************************************************************

	subroutine levels_get_dimension(nl)

	integer nl

	nl = nlvdi

	end subroutine levels_get_dimension

!==================================================================
        end module levels
!==================================================================

