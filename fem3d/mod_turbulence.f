
        module mod_turbulence

        implicit none

        !real tken(0:nlvdim,nkndim)      !turbulent kinetic energ
        !common /tken/tken
        !real eps(0:nlvdim,nkndim)        !dissipation rate
        !common /eps/eps
        !real rls(0:nlvdim,nkndim)        !length scale
        !common /rls/rls
        !save /tken/,/eps/,/rls/

        integer, private, save  :: nkn_turb = 0
        integer, private, save  :: nlv_turb = 0

        real, allocatable, save :: tken(:,:)       ! turbulent kinetic energy 
        real, allocatable, save :: eps(:,:)        ! dissipation rate
        real, allocatable, save :: rls(:,:)        ! length scale

        contains

!************************************************************

        subroutine mod_turbulence_init(nkn,nlv)

        integer  :: nkn
        integer  :: nlv

        if( nkn == nkn_turb .and. nlv == nlv_turb ) return

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
            stop 'error stop mod_turbulence_init: incompatible parms'
          end if
        end if

        if( nkn_turb > 0 ) then
          deallocate(tken)
          deallocate(eps)
          deallocate(rls)
        end if

        nkn_turb = nkn
        nlv_turb = nlv

        if( nkn == 0 ) return

        allocate(tken(0:nlv,nkn))
        allocate(rls(0:nlv,nkn))
        allocate(eps(0:nlv,nkn))

        end subroutine mod_turbulence_init

!************************************************************

        end module mod_turbulence


