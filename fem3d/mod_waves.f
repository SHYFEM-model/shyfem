
!==================================================================
        module mod_waves
!==================================================================

        implicit none

        integer, private, save  :: nkn_waves = 0
        integer, private, save  :: nlv_waves = 0
        integer, private, save  :: nel_waves = 0

        integer, save  :: iwave  = 0	! call parameter to the wave model
        integer, save  :: iwwm   = 0	! type of shyfem-wwm coupling
        integer, save  :: idcoup = 0	! shyfem-wwm coupling time step [s]

        real, allocatable, save :: waveh(:)  !significant wave height [m]
        real, allocatable, save :: wavep(:)  !wave mean period [s]
        real, allocatable, save :: wavepp(:) !wave peak period [s]
        real, allocatable, save :: waved(:)  !mean wave direction [deg]
        real, allocatable, save :: waveov(:) !wave bottom orbital velocity [m/s]

        real, allocatable, save :: wavefx(:,:)	! wave forcing term in x
        real, allocatable, save :: wavefy(:,:)	! wave forcing term in y

!==================================================================
        contains
!==================================================================

        subroutine mod_waves_init(nkn,nel,nlv)

        integer  :: nkn
        integer  :: nel
        integer  :: nlv

        if( nlv == nlv_waves .and. nel == nel_waves .and.
     +      nkn == nkn_waves ) return

        if( nlv > 0 .or. nel > 0 .or. nkn > 0 ) then
          if( nlv == 0 .or. nel == 0 .or. nkn == 0 ) then
            write(6,*) 'nlv,nel,nkn: ',nlv,nel,nkn
            stop 'error stop mod_waves_init: incompatible parameters'
          end if
        end if

        if( nkn_waves > 0 ) then
          deallocate(waveh)
          deallocate(wavep)
          deallocate(wavepp)
          deallocate(waved)
          deallocate(waveov)
          deallocate(wavefx)
          deallocate(wavefy)
        end if

        nlv_waves = nlv
        nel_waves = nel
        nkn_waves = nkn

        if( nkn == 0 ) return

        allocate(waveh(nkn))
        allocate(wavep(nkn))
        allocate(wavepp(nkn))
        allocate(waved(nkn))
        allocate(waveov(nkn))
        allocate(wavefx(nlv,nel))
        allocate(wavefy(nlv,nel))

	waveh = 0.
	wavep = 0.
	wavepp = 0.
	waved = 0.
	waveov = 0.
	wavefx = 0.
	wavefy = 0.

        end subroutine mod_waves_init

!==================================================================
        end module mod_waves
!==================================================================

