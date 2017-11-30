
	module mod_bnd_aux

	implicit none

	integer, private, save :: nkn_bnd_aux = 0
        integer, private, save :: nel_bnd_aux = 0
        
        real, allocatable, save :: rvv(:)	!momentum input (2D)
        real, allocatable, save :: ruv(:)
        real, allocatable, save :: crad(:)	!radiation condition

	contains

*******************************************************************

	subroutine mod_bnd_aux_init(nkn,nel)

	integer nkn
        integer nel

        if( nkn == nkn_bnd_aux .and. nel == nel_bnd_aux ) return

        if( nkn > 0 .or. nel > 0 ) then
          if( nkn == 0 .or. nel == 0 ) then
            write(6,*) 'nkn,nel: ',nkn,nel
            stop 'error stop mod_bnd_aux_init: incompatible params'
          end if
        end if

	if( nkn_bnd_aux > 0 ) then
          deallocate(ruv)
          deallocate(rvv)
          deallocate(crad)
        end if

        nel_bnd_aux = nel
        nkn_bnd_aux = nkn        
        
        if( nkn == 0 .or. nel == 0 ) return
        
        allocate(ruv(nkn))
        allocate(rvv(nkn))
        allocate(crad(nel))
        
	ruv = 0.
	rvv = 0.
	crad = 0.

        end subroutine mod_bnd_aux_init 

!*****************************************************

        end module mod_bnd_aux

