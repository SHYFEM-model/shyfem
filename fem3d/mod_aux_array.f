
        module mod_aux_array

        implicit none

        integer, private, save :: nkn_aux_array = 0
        integer, private, save :: nlv_aux_array = 0
        integer, private, save :: nel_aux_array = 0

        real, allocatable, save :: v1v(:)
        real, allocatable, save :: v2v(:)
        real, allocatable, save :: v3v(:)
        real, allocatable, save :: ve1v(:)
        real, allocatable, save :: saux1(:,:)
        real, allocatable, save :: sauxe1(:,:)
        real, allocatable, save :: sauxe2(:,:)

	contains


!************************************************************

        subroutine mod_aux_array_init(nkn,nel,nlv)

        integer nkn
        integer nel
        integer nlv

        if( nkn == nkn_aux_array .and. nlv == nlv_aux_array 
     +		.and. nel == nel_aux_array) return 

        if( nkn > 0 .or. nlv > 0 .or. nel > 0 ) then
          if( nkn == 0 .or. nlv == 0 .or. nel == 0 ) then
            write(6,*) 'nkn,nel,nlv: ',nkn,nel,nlv
            stop 'error stop mod_aux_array_init: incompatible params'
          end if
        end if

        if( nkn_aux_array > 0 ) then
          deallocate(v1v)
          deallocate(v2v)
          deallocate(v3v)
          deallocate(ve1v)
          deallocate(saux1)
          deallocate(sauxe1)
          deallocate(sauxe2)
        end if

        nkn_aux_array = nkn
        nlv_aux_array = nlv
	nel_aux_array = nel

        if( nkn == 0 ) return

          allocate(v1v(nkn))
          allocate(v2v(nkn))
          allocate(v3v(nkn))
	  allocate(ve1v(nel))

          allocate(saux1(nlv,nkn))
	
          allocate(sauxe1(nlv,nel))
          allocate(sauxe2(nlv,nel))

        end subroutine mod_aux_array_init

!************************************************************

        end module mod_aux_array
