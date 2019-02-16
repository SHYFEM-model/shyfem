
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

	module mod_internal

	implicit none

        integer, private, save :: nkn_internal = 0
        integer, private, save :: nel_internal = 0
        integer, private, save :: nlv_internal = 0
        
        real, allocatable, save :: rdistv(:)
        real, allocatable, save :: fcorv(:)
        real, allocatable, save :: fxv(:,:)
        real, allocatable, save :: fyv(:,:)
        real, allocatable, save :: momentxv(:,:)
        real, allocatable, save :: momentyv(:,:)
        real, allocatable, save :: iuvfix(:)
        double precision, allocatable, save :: ddxv(:,:)
        double precision, allocatable, save :: ddyv(:,:)

        contains

*******************************************************************

	subroutine mod_internal_init(nkn,nel,nlv)

	integer nkn
        integer nel
        integer nlv

        if( nkn == nkn_internal .and. nel == nel_internal
     +  .and. nlv == nlv_internal) return

        if( nkn > 0 .or. nel > 0 .or. nlv > 0) then
          if( nkn == 0 .or. nel == 0 .or. nlv == 0) then
            write(6,*) 'nkn,nel,nlv: ',nkn,nel,nlv
            stop 'error stop mod_internal_init: incompatible params'
          end if
        end if

	if( nkn_internal > 0 ) then
         deallocate(rdistv)
         deallocate(fcorv)
         deallocate(fxv)
         deallocate(fyv)
         deallocate(momentxv)
         deallocate(momentyv)
         deallocate(iuvfix)
         deallocate(ddxv)
         deallocate(ddyv)
        end if

        nel_internal = nel
        nkn_internal = nkn 
        nlv_internal = nlv       
        
        if( nkn == 0 ) return
        
         allocate (rdistv(nkn))
         allocate (fcorv(nel))
         allocate (fxv(nlv,nel))
         allocate (fyv(nlv,nel))
         allocate (momentxv(nlv,nkn))
         allocate (momentyv(nlv,nkn))
         allocate (iuvfix(nel))
         allocate (ddxv(2*nlv,nel))
         allocate (ddyv(2*nlv,nel))
        
        end subroutine mod_internal_init 

!*****************************************************

        end module mod_internal

