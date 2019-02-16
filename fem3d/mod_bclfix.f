
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

        module mod_bclfix

        implicit none

        integer, private, save :: nel_bclfix = 0
        integer, private, save :: nkn_bclfix = 0
        integer, private, save :: nlv_bclfix = 0

        real, allocatable, save :: ielfix(:,:)	! total number of nodes for ele
        real, allocatable, save :: tnudgev(:)	! nudging coefficient array
        real, allocatable, save :: ubound(:,:)	! current velocity in x for boundary [m/s]
        real, allocatable, save :: vbound(:,:)	! current velocity in y for boundary [m/s]

        contains

!************************************************************

        subroutine mod_bclfix_init(nkn,nel,nlv)

        integer  :: nkn
        integer  :: nel
        integer  :: nlv

        if( nlv == nlv_bclfix .and. nel == nel_bclfix .and. 
     +      nkn == nkn_bclfix ) return

        if( nlv > 0 .or. nel > 0 .or. nkn > 0 ) then
          if( nlv == 0 .or. nel == 0 .or. nkn == 0 ) then
            write(6,*) 'nlv,nel,nkn: ',nlv,nel,nkn
            stop 'error stop mod_bclfix_init: incompatible parameters'
          end if
        end if

        if( nkn_bclfix > 0 ) then 
          deallocate(ielfix)
          deallocate(tnudgev)
          deallocate(ubound)
          deallocate(vbound)
        end if

        nlv_bclfix = nlv
        nel_bclfix = nel
        nkn_bclfix = nkn

        if( nkn == 0 ) return

        allocate(ielfix(0:3,nel))
        allocate(tnudgev(nel))
        allocate(ubound(nlv,nkn))
        allocate(vbound(nlv,nkn))

        end subroutine mod_bclfix_init

!************************************************************

        end module mod_bclfix

