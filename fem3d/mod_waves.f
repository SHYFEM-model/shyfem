
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

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

        double precision, save  :: da_wav(4) = 0

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

