
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
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 16.12.2015	ggu	changed VERS_7_3_16
! 28.04.2016	ggu	changed VERS_7_5_9
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 09.10.2022	ggu	set wlnv, wlov to zero

!**************************************************************************

	module mod_hydro_vel

	implicit none

	integer, private, save :: nkn_hydro_vel = 0
	integer, private, save :: nel_hydro_vel = 0
	integer, private, save :: nlv_hydro_vel = 0

	real, allocatable, save :: ulov(:,:)
	real, allocatable, save :: ulnv(:,:)
	real, allocatable, save :: vlov(:,:)
	real, allocatable, save :: vlnv(:,:)
	real, allocatable, save :: wlov(:,:)
	real, allocatable, save :: wlnv(:,:)

	contains

!************************************************************

        subroutine mod_hydro_vel_init(nkn,nel,nlv)

        integer nkn, nel, nlv

        if( nkn == nkn_hydro_vel .and. nel == nel_hydro_vel .and.
     +      nlv == nlv_hydro_vel ) return

        if( nel > 0 .or. nkn > 0 .or. nlv > 0 ) then
          if( nel == 0 .or. nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nel,nkn,nlv: ',nel,nkn,nlv
            stop 'error stop mod_hydro_vel_init: incompatible params'
          end if
        end if

        if( nkn_hydro_vel > 0 ) then
          deallocate(ulov)
          deallocate(ulnv)
          deallocate(vlov)
          deallocate(vlnv)
          deallocate(wlov)
          deallocate(wlnv)
        end if

        nkn_hydro_vel = nkn
        nel_hydro_vel = nel
        nlv_hydro_vel = nlv

        if( nkn == 0 ) return

        allocate(ulov(nlv,nel))
        allocate(ulnv(nlv,nel))
        allocate(vlov(nlv,nel))
        allocate(vlnv(nlv,nel))
        allocate(wlov(0:nlv,nkn))
        allocate(wlnv(0:nlv,nkn))

	!ulnv = -999.
	!vlnv = -999.
	!ulov = -999.
	!vlov = -999.
	ulnv = 0.
	vlnv = 0.
	ulov = 0.
	vlov = 0.
	wlnv = 0.
	wlov = 0.

        end subroutine mod_hydro_vel_init

!************************************************************

        end module mod_hydro_vel
