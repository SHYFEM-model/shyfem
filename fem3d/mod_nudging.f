
	module mod_nudging

	implicit none

        !real andgzv(nkndim)             !contribution to z-computation
        !common /andgzv/andgzv
        !save /andgzv/

	integer, private, save :: nkn_nudging = 0
        
        real, allocatable, save :: andgzv(:)

	contains

*******************************************************************

	subroutine mod_nudging_init(nkn)

	integer nkn

        if( nkn == nkn_nudging ) return

	if( nkn_nudging > 0 ) then
          deallocate(andgzv)
        end if

        nkn_nudging = nkn        
        
        if( nkn == 0 ) return
        
        allocate(andgzv(nkn))
        
        end subroutine mod_nudging_init 

!*****************************************************

        end module mod_nudging

