
!==================================================================
        module bnd_geom
!==================================================================

!----------------------------------------------------------------------
! iopbnd(k) = 0         no open BC
! iopbnd(k) > 0         external open BC (ibtyp=1,2)
! iopbnd(k) < 0         internal open BC (ibtyp=3)
!----------------------------------------------------------------------

        implicit none

        integer, private, save :: nkn_bnd_geom = 0
        integer, private, save :: nrb_bnd_geom = 0

        integer, allocatable, save :: irv(:)
        integer, allocatable, save :: ierv(:,:)
        double precision, allocatable, save :: rhv(:)
        double precision, allocatable, save :: rlv(:)
        double precision, allocatable, save :: rrv(:)
        integer, allocatable, save :: iopbnd(:)

!==================================================================
	contains
!==================================================================

        subroutine mod_bnd_geom_init(nkn,nrb,nkex)

        integer nkn
        integer nrb
	integer,optional :: nkex
	integer naux

	if( nkn == nkn_bnd_geom .and. nrb == nrb_bnd_geom ) return

        if( nkn_bnd_geom > 0 ) then
          deallocate(irv)
          deallocate(ierv)
          deallocate(rhv)
          deallocate(rlv)
          deallocate(rrv)
          deallocate(iopbnd)
        end if

        nkn_bnd_geom = nkn
        nrb_bnd_geom = nrb

        if( nkn == 0 ) return

	naux = max(1,nrb)

        allocate(irv(naux))
        allocate(ierv(2,naux))
        allocate(rhv(naux))
        allocate(rlv(naux))
        allocate(rrv(naux))
	if(present(nkex)) then
          allocate(iopbnd(nkex))
	else
          allocate(iopbnd(nkn))
	end if

	irv = 0
	ierv = 0
	rhv = 0.
	rlv = 0.
	rrv = 0.
	iopbnd = 0

        end subroutine mod_bnd_geom_init

!***************************************************************

        subroutine mod_irv_init(nrb)

	integer nrb

	integer ndim
        integer, allocatable :: irv_aux(:)

	ndim = nrb_bnd_geom

        if( ndim == 0 ) then
          ndim = 10
          allocate(irv(ndim))
        else if( nrb > ndim ) then
          ndim = ndim*2
          allocate(irv_aux(ndim))
          irv_aux(1:ndim/2) = irv(1:ndim/2)
          call move_alloc(irv_aux,irv)
        end if

	nrb_bnd_geom = ndim
	!write(6,*) 'mod_irv_init: ',nrb,ndim

        end subroutine mod_irv_init

!***************************************************************

        subroutine mod_bnd_geom_reinit(nkn,nrb,nkex)

        integer nkn
        integer nrb
	integer,optional :: nkex

	integer, allocatable :: irv_aux(:)

	if( nrb > nrb_bnd_geom ) then
	  write(6,*) 'nrb,nrb_bnd_geom: ',nrb,nrb_bnd_geom
	  stop 'error stop mod_bnd_geom_reinit: nrb > nrb_bnd_geom'
	end if

	allocate(irv_aux(nrb))
	irv_aux(1:nrb) = irv(1:nrb)

	write(6,*) 'mod_bnd_geom_reinit: ',nkn,nrb,nrb_bnd_geom

	nrb_bnd_geom = 0
	deallocate(irv)

	if(present(nkex)) then
          call mod_bnd_geom_init(nkn,nrb,nkex)
	else
          call mod_bnd_geom_init(nkn,nrb)
	end if

	irv = 0
	irv(1:nrb) = irv_aux(1:nrb)
	deallocate(irv_aux)

	!write(6,*) 'mod_bnd_geom_reinit: ',nrb
	!write(6,'(8i9)') irv

        end subroutine mod_bnd_geom_reinit

!***************************************************************

        subroutine mod_bnd_geom_info

	integer iu,i
	integer nkn,nrb

	iu = 88
	nkn = nkn_bnd_geom
	nrb = nrb_bnd_geom

	write(iu,*) 'mod_bnd_geom_info: ',nkn,nrb
        write(iu,*) 'irv: ',nrb,(irv(i),i=1,nrb)
        write(iu,*) 'ierv: ',(ierv(1,i),i=1,nrb)
        write(iu,*) 'ierv: ',(ierv(2,i),i=1,nrb)
        write(iu,*) 'rhv: ',(rhv(i),i=1,nrb)
        write(iu,*) 'rlv: ',(rlv(i),i=1,nrb)
        write(iu,*) 'rrv: ',(rrv(i),i=1,nrb)
        write(iu,*) 'iopbnd: '
        do i=1,nkn
          if( iopbnd(i) .ne. 0 ) write(iu,*) i,iopbnd(i)
        end do
	write(iu,*) 'mod_bnd_geom_info end'

        end subroutine mod_bnd_geom_info

!***************************************************************

	function is_boundary(k)

	logical is_boundary
	integer k

	is_boundary = iopbnd(k) .ne. 0

	end function is_boundary

!***************************************************************

	function is_external_boundary(k)

	logical is_external_boundary
	integer k

	is_external_boundary = iopbnd(k) .gt. 0

	end function is_external_boundary

!***************************************************************

	function is_internal_boundary(k)

	logical is_internal_boundary
	integer k

	is_internal_boundary = iopbnd(k) .lt. 0

	end function is_internal_boundary

!***************************************************************

	function is_inner(k)

	logical is_inner
	integer k

	is_inner = iopbnd(k) .eq. 0

	end function is_inner

!==================================================================

        subroutine mod_irv_initialize

	implicit none

        call mod_irv_init(0)

	end subroutine mod_irv_initialize

!***************************************************************

!==================================================================
        end module bnd_geom
!==================================================================
