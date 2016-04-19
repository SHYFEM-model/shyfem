

        module mod_gotm_aux

        implicit none

        integer, private, save :: nkn_gotm_aux = 0
        integer, private, save :: nlv_gotm_aux = 0

        double precision, allocatable, save :: numv_gotm(:,:)
        double precision, allocatable, save :: nuhv_gotm(:,:)
        double precision, allocatable, save :: tken_gotm(:,:)
        double precision, allocatable, save :: eps_gotm(:,:)
        double precision, allocatable, save :: rls_gotm(:,:)

        real, allocatable, save :: shearf2(:,:)
        real, allocatable, save :: buoyf2(:,:)

        contains


!************************************************************

        subroutine mod_gotm_aux_init(nkn,nlv)

        integer nlv
        integer nkn


	if( nkn == nkn_gotm_aux .and. nlv == nlv_gotm_aux ) return   

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
            stop 'error stop mod_gotm_aux_init: incompatible parameters'
          end if
        end if


        if( nkn_gotm_aux > 0 ) then   
          deallocate(numv_gotm)
          deallocate(nuhv_gotm)
          deallocate(tken_gotm)
          deallocate(eps_gotm)
          deallocate(rls_gotm)
          deallocate(shearf2)
          deallocate(buoyf2)
        end if

        nkn_gotm_aux = nkn
	nlv_gotm_aux = nlv

        if( nkn == 0 ) return              

          allocate(numv_gotm(0:nlv,nkn))
          allocate(nuhv_gotm(0:nlv,nkn))
          allocate(tken_gotm(0:nlv,nkn))
          allocate(eps_gotm(0:nlv,nkn))
          allocate(rls_gotm(0:nlv,nkn))
          allocate(shearf2(nlv,nkn))
          allocate(buoyf2(nlv,nkn))

        end subroutine mod_gotm_aux_init

!************************************************************

        end module mod_gotm_aux


