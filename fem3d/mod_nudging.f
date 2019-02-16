
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
	module mod_nudging
!==================================================================

	implicit none

	integer, private, save :: nkn_nudging = 0
	integer, private, save :: nel_nudging = 0
	integer, private, save :: nlv_nudging = 0
        
	real, save :: anpar = 0.		!implicit parameter for nudging
	real, save :: taudefvel = 0.		!default tau for velocities

        real, allocatable, save :: andgzv(:)	!contribution to zeta
        real, allocatable, save :: tauvel(:,:)	!weighting for vel
        real, allocatable, save :: uobs(:,:)	!observations for x vel
        real, allocatable, save :: vobs(:,:)	!observations for y vel

!==================================================================
	contains
!==================================================================

	subroutine mod_nudging_init(nkn,nel,nlv)

	integer nkn,nel,nlv

        if( nkn == nkn_nudging .and. nel == nel_nudging .and.
     +      nlv == nlv_nudging ) return

        if( nel > 0 .or. nkn > 0 .or. nlv > 0 ) then
          if( nel == 0 .or. nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nel,nkn,nlv: ',nel,nkn,nlv
            stop 'error stop mod_nudging_init: incompatible parameters'
          end if
        end if

        if( nkn == nkn_nudging ) return

	if( nkn_nudging > 0 ) then
          deallocate(andgzv)
          deallocate(tauvel)
          deallocate(uobs)
          deallocate(vobs)
        end if

        nkn_nudging = nkn
        nel_nudging = nel
        nlv_nudging = nlv
        
        if( nkn == 0 ) return
        
        allocate(andgzv(nkn))
        allocate(tauvel(nlv,nel))
        allocate(uobs(nlv,nel))
        allocate(vobs(nlv,nel))

	andgzv = 0.
	tauvel = 0.
	uobs = 0.
	vobs = 0.
        
        end subroutine mod_nudging_init 

!==================================================================
        end module mod_nudging
!==================================================================

