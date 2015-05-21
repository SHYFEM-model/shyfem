
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

        integer, save, private :: nlv_aux = 0
        real, save, allocatable :: hlv_aux(:)

!==================================================================
        contains
!==================================================================

	subroutine levels_aux_init(nl)

	integer nl

	if( nlv_aux == nl ) return

	if( nlv_aux > 0 ) then
	  deallocate(hlv_aux)
	end if

	nlv_aux = nl

	if( nl > 0 ) then
	  allocate(hlv_aux(nl))
	end if

	end subroutine levels_aux_init

!******************************************************************

	subroutine levels_init(nkn,nel,nl)

	integer nkn,nel,nl

	nlv = nl
	nlvdi = nl
	nlv_alloc = nl

	allocate(ilhv(nel))
	allocate(ilmv(nel))
	allocate(ilhkv(nkn))
	allocate(ilmkv(nkn))

	allocate(hlv(nl))
	allocate(hldv(nl))
	
	end subroutine levels_init

!******************************************************************

	subroutine levels_reinit(nl)

! re-allocates arrays depending on nl

	integer nl

	if( nlv_alloc == nl ) return

	if( nlv_alloc > 0 ) then
	  deallocate(hlv)
	  deallocate(hldv)
	end if

	nlv = nl
	nlvdi = nl
	nlv_alloc = nl

	if( nl > 0 ) then
	  allocate(hlv(nl))
	  allocate(hldv(nl))
	end if
	
	end subroutine levels_reinit

!******************************************************************

	subroutine copy_hlv(nl,hl)

	integer nl
	real hl(nl)

	integer l

	if( nl < nlv_aux ) then
	  write(6,*) 'nl,nlv: ',nl,nlv_aux
	  stop 'error stop copy_hlv: nlv'
	end if

	nl = nlv_aux
	do l=1,nlv_aux
	  hl(l) = hlv_aux(l)
	end do

	end subroutine copy_hlv

!==================================================================
        end module levels
!==================================================================

	subroutine read_hlv

	use levels
	use nls

	implicit none

	integer n

	n = nls_read_vector()
	call levels_aux_init(n)
	call nls_copy_real_vect(n,hlv_aux)

	end subroutine read_hlv

!******************************************************************

	subroutine copy_hlv1(nl,hl)

	use levels

	implicit none

	integer nl
	real hl(nl)

	integer l

	call copy_hlv(nl,hl)

	write(6,*) 'hlv copied: ',nl
	write(6,*) (hl(l),l=1,nl)

	end

!******************************************************************

