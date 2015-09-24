
!==================================================================
	module mod_bndo
!==================================================================

	implicit none

	integer, private, save  :: nrb_bndo = 0
	integer, private, save  :: ngr_bndo = 0

	integer, save :: nbndo = 0	!total number of OB nodes
	integer, save :: ndebug = 0	!unit number for debug messages

	integer, save :: kbcdim = 0	!to be deleted later
	integer, save :: kopdim = 0	!to be deleted later

	integer, allocatable, save :: nopnod(:)	!number of nodes close to OB
	integer, allocatable, save :: ibcnod(:)	!number of boundary
	integer, allocatable, save :: kbcnod(:)	!number of boundary node
	integer, allocatable, save :: itynod(:)	!type of boundary
	integer, allocatable, save :: nopnodes(:,:)	!nodes close to OB
	real, allocatable, save :: xynorm(:,:)		!normal direction
	real, allocatable, save :: wopnodes(:,:)	!weights

!==================================================================
	contains
!==================================================================

	subroutine mod_bndo_init(ngr,nrb)

	integer ngr,nrb

	integer nlk,naux

        if( ngr == ngr_bndo .and. nrb == nrb_bndo ) return

        if( ngr_bndo > 0 ) then
          deallocate(nopnod)
          deallocate(ibcnod)
          deallocate(kbcnod)
          deallocate(itynod)
          deallocate(nopnodes)
          deallocate(xynorm)
          deallocate(wopnodes)
        end if

        ngr_bndo = ngr
        nrb_bndo = nrb

	kbcdim = nrb		!to be deleted later
	kopdim = ngr		!to be deleted later

	if( ngr == 0 ) return

	naux = max(1,nrb)

        allocate(nopnod(naux))
        allocate(ibcnod(naux))
        allocate(kbcnod(naux))
        allocate(itynod(naux))
        allocate(nopnodes(ngr,naux))
        allocate(xynorm(2,naux))
        allocate(wopnodes(ngr,naux))

	end subroutine mod_bndo_init

	subroutine mod_bndo_info

	write(6,*) 'mod_bndo_info ================='
	write(6,*) 'ngr_bndo: ',ngr_bndo
	write(6,*) 'nrb_bndo: ',nrb_bndo
	write(6,*) 'nbndo: ',nbndo
	write(6,*) 'mod_bndo_info end ================='

	end subroutine mod_bndo_info

!==================================================================
	end module mod_bndo
!==================================================================

