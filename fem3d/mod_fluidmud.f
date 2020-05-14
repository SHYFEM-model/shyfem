
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2016,2019  Georg Umgiesser
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

! revision log :
!
! 10.07.2015	ggu	changed VERS_7_1_50
! 12.10.2015	ggu	changed VERS_7_3_3
! 28.04.2016	ggu	changed VERS_7_5_9
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

        module mod_fluidmud

        implicit none

        integer, private, save  :: nkn_fmud = 0
        integer, private, save  :: nlv_fmud = 0

        double precision, private, save     :: rhosed  ! Mud primary particle density (kg/m3)
        double precision, private, save     :: dm0     ! ccf viene definita anche in submud.f ???
        double precision, private, save     :: nf      ! ccf scalare o array???

        real, allocatable, save :: z0bkmud(:)          ! bottom roughenss on nodes for mud
        real, allocatable, save :: mudc(:,:)           ! Fluid mud concentrationarray (kg/m3)
        double precision, allocatable, save :: rhomud(:,:) ! Mud floc part. density (kg/m3)
        real, allocatable, save :: visv_yield(:,:)     ! viscosity (mud)                     
        real, allocatable, save :: diff_yield(:,:)     ! diffusivity (mud)                     
        real, allocatable, save :: lambda(:,:)         ! Structural parameter                  
        real, allocatable, save :: vts(:,:)            ! Rheological Viscosity [m2/s]          
        real, allocatable, save :: dmf_mud(:,:)        ! Floc size array.                     
        real, allocatable, save :: wprvs(:,:)          ! Water density (kg/m3)              

        contains

!************************************************************

        subroutine mod_fluidmud_init(nkn,nlv)

        integer  :: nkn
        integer  :: nlv

        if( nkn == nkn_fmud .and. nlv == nlv_fmud ) return

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
            stop 'error stop mod_fluidmud_init: incompatible parms'
          end if
        end if

        if( nkn_fmud > 0 ) then
          deallocate(z0bkmud)
          deallocate(mudc)
          deallocate(rhomud)
          deallocate(visv_yield)
          deallocate(diff_yield)
          deallocate(lambda)
          deallocate(vts)
          deallocate(dmf_mud)
          deallocate(wprvs)
        end if

        nkn_fmud = nkn
        nlv_fmud = nlv

        if( nkn == 0 ) return

        allocate(z0bkmud(nkn))
        allocate(mudc(nlv,nkn))
        allocate(rhomud(nlv,nkn))
        allocate(visv_yield(0:nlv,nkn))
        allocate(diff_yield(0:nlv,nkn))
        allocate(lambda(nlv,nkn))
        allocate(vts(0:nlv,nkn))
        allocate(dmf_mud(nlv,nkn))
        allocate(wprvs(0:nlv,nkn))

        end subroutine mod_fluidmud_init

!************************************************************

        subroutine mod_fluidmud_dummy_init(nkn,nlv)

	integer nkn,nlv

        allocate(vts(0:nlv,nkn))
	vts = 0.

        end subroutine mod_fluidmud_dummy_init

!************************************************************

        end module mod_fluidmud

