
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2017,2019  Georg Umgiesser
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
! 18.09.2015	ggu	changed VERS_7_2_3
! 16.12.2015	ggu	changed VERS_7_3_16
! 12.01.2017	ggu	changed VERS_7_5_21
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

        module mod_ts

        implicit none

	integer, private, save :: nkn_ts = 0
	integer, private, save :: nlv_ts = 0

        real, allocatable, save :: rhov(:,:)
        real, allocatable, save :: saltv(:,:)
        real, allocatable, save :: tempv(:,:)

        real, allocatable, save :: sobsv(:,:)
        real, allocatable, save :: stauv(:,:)
        real, allocatable, save :: tobsv(:,:)
        real, allocatable, save :: ttauv(:,:)

        real, allocatable, save :: bpresv(:,:)
        real, allocatable, save :: bpresxv(:,:)
        real, allocatable, save :: bpresyv(:,:)

        contains

!************************************************************

        subroutine mod_ts_init(nkn,nlv)

        integer nkn
        integer nlv

        if( nkn == nkn_ts .and. nlv == nlv_ts) return

	if( nkn > 0 .or. nlv > 0 ) then
	  if( nkn == 0 .or. nlv == 0 ) then
	    write(6,*) 'nkn,nlv: ',nkn,nlv
	    stop 'error stop mod_ts_init: incompatible parameters'
	  end if
	end if

        if( nkn_ts > 0 ) then
          deallocate(rhov)
          deallocate(saltv)
          deallocate(tempv)
          deallocate(sobsv)
          deallocate(tobsv)
          deallocate(stauv)
          deallocate(ttauv)
          deallocate(bpresv)
          deallocate(bpresxv)
          deallocate(bpresyv)
        end if

        nkn_ts = nkn
	nlv_ts = nlv

        if( nkn == 0 ) return

          allocate(rhov(nlv,nkn))
          allocate(saltv(nlv,nkn))
          allocate(tempv(nlv,nkn))
          allocate(sobsv(nlv,nkn))
          allocate(tobsv(nlv,nkn))
          allocate(stauv(nlv,nkn))
          allocate(ttauv(nlv,nkn))
          allocate(bpresv(nlv,nkn))
          allocate(bpresxv(nlv,nkn))
          allocate(bpresyv(nlv,nkn))

	rhov = 0.
	saltv = 0.
	tempv = 0.
	sobsv = 0.
	tobsv = 0.
	stauv = 0.
	ttauv = 0.
	bpresv = 0.
	bpresxv = 0.
	bpresyv = 0.

        end subroutine mod_ts_init

!************************************************************

        end module mod_ts

