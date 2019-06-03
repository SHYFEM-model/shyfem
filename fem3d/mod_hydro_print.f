
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

! module for hydrodynamic print values
!
! revision log :
!
! 10.07.2015	ggu	changed VERS_7_1_50
! 28.04.2016	ggu	changed VERS_7_5_9
! 17.06.2016	ggu	wprv is now from 1:nlv
! 16.02.2019	ggu	changed VERS_7_5_60
!
!============================================================
	module mod_hydro_print
!============================================================

	implicit none

        integer, private, save :: nkn_hydro_print = 0
        integer, private, save :: nlv_hydro_print = 0

        real, allocatable, save :: uprv(:,:)
        real, allocatable, save :: vprv(:,:)
        real, allocatable, save :: upro(:,:)
        real, allocatable, save :: vpro(:,:)
        real, allocatable, save :: wprv(:,:)
        real, allocatable, save :: up0v(:), vp0v(:)
        real, allocatable, save :: xv(:,:)

!============================================================
	contains
!============================================================

        subroutine mod_hydro_print_init(nkn,nlv)

        integer nkn, nlv

	if( nkn == nkn_hydro_print .and. nlv == nlv_hydro_print ) return

        if( nlv > 0 .or. nkn > 0 ) then
          if( nlv == 0 .or. nkn == 0 ) then
            write(6,*) 'nlv,nkn: ',nlv,nkn
            stop 'error stop mod_hydro_print_init: incompatible params'
          end if
        end if

	if( nkn_hydro_print > 0 ) then
          deallocate(uprv)
          deallocate(vprv)
          deallocate(upro)
          deallocate(vpro)
          deallocate(wprv)

          deallocate(up0v)
          deallocate(vp0v)

          deallocate(xv)
        end if

        nkn_hydro_print = nkn
        nlv_hydro_print = nlv

        if( nkn == 0 ) return

        allocate(uprv(nlv,nkn))
        allocate(vprv(nlv,nkn))
        allocate(upro(nlv,nkn))
        allocate(vpro(nlv,nkn))
        allocate(wprv(nlv,nkn))

	allocate(up0v(nkn))
	allocate(vp0v(nkn))

	allocate(xv(3,nkn))

        end subroutine mod_hydro_print_init

!============================================================
        end module mod_hydro_print
!============================================================

