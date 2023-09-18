
!==================================================================
        module bnd_dynamic
!==================================================================

        implicit none

        integer, private, save :: nkn_bnd_dynamic = 0
        integer, private, save :: nlv_bnd_dynamic = 0

        double precision, parameter, private :: flag = -9988765.0

        double precision, allocatable, save :: rzv(:)
        double precision, allocatable, save :: rqv(:)
        double precision, allocatable, save :: rqpsv(:)
        double precision, allocatable, save :: rqdsv(:)
        double precision, allocatable, save :: mfluxv(:,:)

!==================================================================
        contains
!==================================================================

        subroutine mod_bnd_dynamic_init(nkn,nlv)

        integer nkn
        integer nlv

        if( nkn == nkn_bnd_dynamic .and. nlv == nlv_bnd_dynamic ) return

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
	    stop 'error stop mod_bnd_dynamic_init: incompatible params'
          end if
        end if

        if( nkn_bnd_dynamic > 0 ) then
          deallocate(rzv)
          deallocate(rqv)
          deallocate(rqpsv)
          deallocate(rqdsv)
          deallocate(mfluxv)
        end if

        nkn_bnd_dynamic = nkn
        nlv_bnd_dynamic = nlv

        if( nkn == 0 ) return

          allocate(rzv(nkn))
          allocate(rqv(nkn))
          allocate(rqpsv(nkn))
          allocate(rqdsv(nkn))
          allocate(mfluxv(nlv,nkn))

          mfluxv = 0.

        end subroutine mod_bnd_dynamic_init

!******************************************************************

        function is_zeta_boundary(k)

        logical is_zeta_boundary
        integer k

        is_zeta_boundary = rzv(k) /= flag

        end function is_zeta_boundary

!==================================================================
        end module bnd_dynamic
!==================================================================


