

        module mod_diff_aux

        implicit none

        integer, private, save :: nel_diff_aux = 0

        real, allocatable, save :: wdifhv(:,:,:)  !weights for horizontal diff

        contains


!************************************************************

        subroutine mod_diff_aux_init(nel)

        integer nel

        if( nel == nel_diff_aux ) return   

        if( nel_diff_aux > 0 ) then        
          deallocate(wdifhv)
        end if

        nel_diff_aux = nel

        if( nel == 0 ) return              

        allocate(wdifhv(3,3,nel))

        end subroutine mod_diff_aux_init

!************************************************************

        end module mod_diff_aux


