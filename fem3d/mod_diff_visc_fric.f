
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

	module mod_diff_visc_fric

	implicit none

	integer, private, save :: nkn_diff_visc_fric = 0
	integer, private, save :: nel_diff_visc_fric = 0
	integer, private, save :: nlv_diff_visc_fric = 0

	real, allocatable, save :: rfricv(:)	!friction term
	real, allocatable, save :: rcdv(:)	!bottom drag coefficient
	real, allocatable, save :: bnstressv(:)	!normalized bottom stress
	real, allocatable, save :: czv(:)	!friction parameter (as given)

	real, allocatable, save :: difhv(:,:)
	real, allocatable, save :: visv(:,:)
	real, allocatable, save :: difv(:,:)

	contains

!************************************************************

        subroutine mod_diff_visc_fric_init(nkn,nel,nlv)

        integer nkn, nel, nlv

        if( nkn == nkn_diff_visc_fric .and. nel == nel_diff_visc_fric
     +      .and. nlv == nlv_diff_visc_fric ) return

        if( nel > 0 .or. nkn > 0 .or. nlv > 0 ) then
          if( nel == 0 .or. nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nel,nkn,nlv: ',nel,nkn,nlv
	    stop 'error stop mod_diff_visc_fric_init: incompatible params'
          end if
        end if

        if( nkn_diff_visc_fric > 0 ) then
          deallocate(rfricv)
          deallocate(rcdv)
          deallocate(bnstressv)
          deallocate(czv)
          deallocate(difhv)
          deallocate(visv)
          deallocate(difv)
        end if

        nkn_diff_visc_fric = nkn
        nel_diff_visc_fric = nel
        nlv_diff_visc_fric = nlv

        if( nkn == 0 ) return

        allocate(rfricv(nel))
        allocate(rcdv(nel))
        allocate(bnstressv(nel))
        allocate(czv(nel))
        allocate(difhv(nlv,nel))
        allocate(visv(0:nlv,nkn))
        allocate(difv(0:nlv,nkn))

        end subroutine mod_diff_visc_fric_init

!************************************************************

        end module mod_diff_visc_fric
