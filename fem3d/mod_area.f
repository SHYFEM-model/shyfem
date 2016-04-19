
        module mod_area

        implicit none

        integer, private, save :: nkn_area = 0
        integer, private, save :: nlv_area = 0

        real, allocatable, save :: areakv(:,:)

	contains


!************************************************************

        subroutine mod_area_init(nkn,nlv)

        integer nkn
        integer nlv

        if( nkn == nkn_area .and. nlv == nlv_area) return   !!!!!giusto

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
            stop 'error stop mod_area_init: incompatible parameters'
          end if
        end if

        if( nkn_area > 0 ) then
          deallocate(areakv)
        end if

        nkn_area = nkn
        nlv_area = nlv

        if( nkn == 0 ) return

        allocate(areakv(nlv,nkn))

        end subroutine mod_area_init

!************************************************************

        end module mod_area

