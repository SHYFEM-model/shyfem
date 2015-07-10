
        module mod_bound_dynamic

        implicit none

        !real rzv(nkndim), rqv(nkndim)
        !common /rzv/rzv, /rqv/rqv

        !real rqpsv(nkndim), rqdsv(nkndim)
        !common /rqpsv/rqpsv, /rqdsv/rqdsv
        	!real evapv(nkndim)	!!!QUESTI ERANO PRECOMMENTATI
        	!common /evapv/evapv	!!!VANNO MESSI O NO?
        !real mfluxv(nlvdim,nkndim)
        !common /mfluxv/mfluxv
	!save /rzv/,/rqv/
	!save /rqpsv/,/rqdsv/,/mfluxv/

        integer, private, save :: nkn_bound_dynamic = 0
        integer, private, save :: nlv_bound_dynamic = 0

        real, allocatable, save :: rzv(:)
        real, allocatable, save :: rqv(:)
        real, allocatable, save :: rqpsv(:)
        real, allocatable, save :: rqdsv(:)
        real, allocatable, save :: mfluxv(:,:)

        contains

!************************************************************

        subroutine mod_bound_dynamic_init(nkn,nlv)

        integer nkn
        integer nlv

        if( nkn == nkn_bound_dynamic .and. 
     +		nlv == nlv_bound_dynamic ) return

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
	    stop 'error stop mod_bound_dynamic_init: incompatible params'
          end if
        end if

        if( nkn_bound_dynamic > 0 ) then
          deallocate(rzv)
          deallocate(rqv)
          deallocate(rqpsv)
          deallocate(rqdsv)
          deallocate(mfluxv)
        end if

        nkn_bound_dynamic = nkn
        nlv_bound_dynamic = nlv

        if( nkn == 0 ) return

          allocate(rzv(nkn))
          allocate(rqv(nkn))
          allocate(rqpsv(nkn))
          allocate(rqdsv(nkn))
          allocate(mfluxv(nlv,nkn))

        end subroutine mod_bound_dynamic_init

!************************************************************

        end module mod_bound_dynamic


