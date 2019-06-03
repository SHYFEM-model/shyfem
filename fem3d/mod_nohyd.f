
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

! revision log :
!
! 10.07.2015	ggu	changed VERS_7_1_50
! 14.06.2016	ggu	changed VERS_7_5_14
! 17.06.2016	ggu	changed VERS_7_5_15
! 18.12.2018	ggu	changed VERS_7_5_52
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

!===================================================================
	module mod_nohyd
!===================================================================

	implicit none

	logical, save :: bnohydro = .false.

        integer, private, save :: nkn_nohyd = 0
        integer, private, save :: nlv_nohyd = 0
        
        real, allocatable, save :: qpnv(:,:)
        real, allocatable, save :: qpov(:,:)
        real, allocatable, save :: qdistv(:) !DWNH

!===================================================================
        contains
!===================================================================

!*******************************************************************

	subroutine mod_nohyd_init(nkn,nlv)

	integer nkn
        integer nlv

        if( nkn == nkn_nohyd .and. nlv == nlv_nohyd ) return

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
            stop 'error stop mod_nohyd_init: incompatible params'
          end if
        end if

	if( nkn_nohyd > 0 ) then
          deallocate(qpnv)
	  deallocate(qpov)
          deallocate(qdistv) !DWNH
        end if

        nkn_nohyd = nkn 
        nlv_nohyd = nlv       
        
        if( nkn == 0 ) return
        
        allocate (qpnv(nlv,nkn))
        allocate (qpov(nlv,nkn))
        allocate (qdistv(nkn)) !DWNH

	qpnv = 0.
	qpov = 0.
        
        end subroutine mod_nohyd_init 

!*******************************************************************

!===================================================================
        end module mod_nohyd
!===================================================================

