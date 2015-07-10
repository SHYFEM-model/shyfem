
!==================================================================
        module mod_tvd
!==================================================================

        implicit none

        integer, private, save :: nel_tvd = 0

	integer, save :: itvd_type = 0

        real, allocatable, save :: tvdupx(:,:,:)
        real, allocatable, save :: tvdupy(:,:,:)
        integer, allocatable, save :: ietvdup(:,:,:)

!==================================================================
	contains
!==================================================================

        subroutine mod_tvd_init(nel)

        integer nel

        if( nel == nel_tvd ) return

        if( nel_tvd > 0 ) then
          deallocate(tvdupx)
          deallocate(tvdupy)
          deallocate(ietvdup)
        end if

        nel_tvd = nel

        if( nel == 0 ) return

        allocate(tvdupx(3,3,nel))
        allocate(tvdupy(3,3,nel))
        allocate(ietvdup(3,3,nel))

        end subroutine mod_tvd_init

!==================================================================
        end module mod_tvd
!==================================================================

