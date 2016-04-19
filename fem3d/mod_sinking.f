
        module mod_sinking

        implicit none

        integer, private, save  :: nkn_sink = 0
        integer, private, save  :: nlv_sink = 0

        real, allocatable, save :: wsinkv(:,:)       ! sinking velocity [m/s]

        contains

!************************************************************

        subroutine mod_sinking_init(nkn,nlv)

        integer  :: nkn
        integer  :: nlv

        if( nkn == nkn_sink .and. nlv == nlv_sink ) return

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
            stop 'error stop mod_sinking_init: incompatible parameters'
          end if
        end if

        if( nkn_sink > 0 ) then
	  deallocate(wsinkv)
	end if

        nkn_sink = nkn
        nlv_sink = nlv

        if( nkn == 0 ) return

        allocate(wsinkv(0:nlv,nkn))

        end subroutine mod_sinking_init

!************************************************************

        end module mod_sinking


