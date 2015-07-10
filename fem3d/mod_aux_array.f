
        module mod_aux_array

        implicit none

        !real v1v(nkndim)
        !common /v1v/v1v
        !real v2v(nkndim)
        !common /v2v/v2v
        !real v3v(nkndim)
        !common /v3v/v3v
        !real ve1v(neldim)
        !common /ve1v/ve1v
        !real saux1(nlvdim,nkndim)
        !common /saux1/saux1
        !real saux2(nlvdim,nkndim)
        !common /saux2/saux2
        !real saux3(nlvdim,nkndim)
        !common /saux3/saux3
        !real saux4(nlvdim,nkndim)
        !common /saux4/saux4
        !real sauxe1(nlvdim,neldim)
        !common /sauxe1/sauxe1
        !real sauxe2(nlvdim,neldim)
        !common /sauxe2/sauxe2
	!save /v1v/,/v2v/,/v3v/,/ve1v/
	!save /saux1/,/saux2/,/saux3/,/saux4/
	!save /sauxe1/,/sauxe2/

        integer, private, save :: nkn_aux_array = 0
        integer, private, save :: nlv_aux_array = 0
        integer, private, save :: nel_aux_array = 0

        real, allocatable, save :: v1v(:)
        real, allocatable, save :: v2v(:)
        real, allocatable, save :: v3v(:)
        real, allocatable, save :: ve1v(:)
        real, allocatable, save :: saux1(:,:)
        real, allocatable, save :: saux2(:,:)
        real, allocatable, save :: saux3(:,:)
        real, allocatable, save :: saux4(:,:)
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
          deallocate(saux2)
          deallocate(saux3)
          deallocate(saux4)
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
          allocate(saux2(nlv,nkn))
          allocate(saux3(nlv,nkn))
          allocate(saux4(nlv,nkn))	
	
          allocate(sauxe1(nlv,nel))
          allocate(sauxe2(nlv,nel))
  

        end subroutine mod_aux_array_init

!************************************************************

        end module mod_aux_array
