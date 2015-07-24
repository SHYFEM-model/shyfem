
	module mod_depth

	implicit none

        integer, private, save :: nkn_depth = 0
        integer, private, save :: nel_depth = 0
        
        real, allocatable, save :: hkv(:)
        real, allocatable, save :: hev(:)

        real, allocatable, save :: hkv_min(:)
        real, allocatable, save :: hkv_max(:)

        contains

*******************************************************************

	subroutine mod_depth_init(nkn,nel)

	integer nkn
        integer nel

        if( nkn == nkn_depth .and. nel == nel_depth ) return

        if( nkn > 0 .or. nel > 0 ) then
          if( nkn == 0 .or. nel == 0 ) then
            write(6,*) 'nkn,nel: ',nkn,nel
            stop 'error stop mod_depth_init: incompatible params'
          end if
        end if

	if( nkn_depth > 0 ) then
          deallocate(hev)
          deallocate(hkv)
          deallocate(hkv_min)
          deallocate(hkv_max)
        end if

        nkn_depth = nkn 
        nel_depth = nel
        
        if( nkn == 0 ) return
        
        allocate(hkv(nkn))
        allocate(hev(nel))
        allocate(hkv_min(nkn))
        allocate(hkv_max(nkn))

        end subroutine mod_depth_init 

!*****************************************************

        end module mod_depth

