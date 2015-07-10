
!==================================================================
        module mod_bnd
!==================================================================

        implicit none

        integer, private, save :: nbc_bnd = 0
	integer, parameter :: nbvdim = 25
	!integer, parameter :: nbndim = 25

	integer, save :: nbc = 0
	integer, save :: nrb = 0

        real, allocatable, save :: bnd(:,:)

!==================================================================
	contains
!==================================================================

        subroutine mod_bnd_init(nb)

        integer nb

        if( nb == nbc_bnd ) return

        if( nbc_bnd > 0 ) then
          deallocate(bnd)
        end if

        nbc_bnd = nb

        if( nb == 0 ) return

        allocate(bnd(nbvdim,nb))

        end subroutine mod_bnd_init

!***************************************************************

	subroutine mod_bnd_adjust(nb)

        integer nb

        integer ndim
        real, allocatable :: bnd_aux(:,:)

        ndim = nbc_bnd

        if( ndim == 0 ) then
          ndim = max(10,2*nb)
          allocate(bnd(nbvdim,ndim))
        else if( nb > ndim ) then
          ndim = ndim*2
          allocate(bnd_aux(nbvdim,ndim))
          bnd_aux(:,1:ndim/2) = bnd(:,1:ndim/2)
          call move_alloc(bnd_aux,bnd)
        end if

	nbc_bnd = ndim

	end subroutine mod_bnd_adjust

!***************************************************************

        subroutine mod_bnd_reinit(nb)

        integer nb

	real, allocatable :: bnd_aux(:,:)

	if( nb > nbc_bnd ) then
	  write(6,*) 'nb,nbc_bnd: ',nb,nbc_bnd
	  stop 'error stop mod_bnd_reinit: nb > nbc_bnd'
	end if

	allocate(bnd_aux(nbvdim,nb))
	bnd_aux(:,1:nb) = bnd(:,1:nb)

        call mod_bnd_init(nb)

	bnd = 0
	bnd(:,1:nb) = bnd_aux(:,1:nb)
	deallocate(bnd_aux)

        end subroutine mod_bnd_reinit

!==================================================================
        end module mod_bnd
!==================================================================

