
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
! 28.04.2016	ggu	changed VERS_7_5_9
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 01.04.2021	ggu	create restart routines here
! 13.10.2022	ggu	restart also with mpi
! 28.04.2023    ggu     update function calls for belem

!**************************************************************************

!==================================================================
        module mod_gotm_aux
!==================================================================

        implicit none

        integer, private, save :: nkn_gotm_aux = 0
        integer, private, save :: nlv_gotm_aux = 0

        double precision, allocatable, save :: numv_gotm(:,:)
        double precision, allocatable, save :: nuhv_gotm(:,:)
        double precision, allocatable, save :: tken_gotm(:,:)
        double precision, allocatable, save :: eps_gotm(:,:)
        double precision, allocatable, save :: rls_gotm(:,:)

        real, allocatable, save :: shearf2(:,:)
        real, allocatable, save :: buoyf2(:,:)

!==================================================================
        contains
!==================================================================

        subroutine mod_gotm_aux_init(nkn,nlv)

        integer nlv
        integer nkn

	if( nkn == nkn_gotm_aux .and. nlv == nlv_gotm_aux ) return   

        if( nkn > 0 .or. nlv > 0 ) then
          if( nkn == 0 .or. nlv == 0 ) then
            write(6,*) 'nkn,nlv: ',nkn,nlv
            stop 'error stop mod_gotm_aux_init: incompatible parameters'
          end if
        end if

        if( nkn_gotm_aux > 0 ) then   
          deallocate(numv_gotm)
          deallocate(nuhv_gotm)
          deallocate(tken_gotm)
          deallocate(eps_gotm)
          deallocate(rls_gotm)
          deallocate(shearf2)
          deallocate(buoyf2)
        end if

        nkn_gotm_aux = nkn
	nlv_gotm_aux = nlv

        if( nkn == 0 ) return              

          allocate(numv_gotm(0:nlv,nkn))
          allocate(nuhv_gotm(0:nlv,nkn))
          allocate(tken_gotm(0:nlv,nkn))
          allocate(eps_gotm(0:nlv,nkn))
          allocate(rls_gotm(0:nlv,nkn))
          allocate(shearf2(nlv,nkn))
          allocate(buoyf2(nlv,nkn))

        end subroutine mod_gotm_aux_init

!==================================================================
        end module mod_gotm_aux
!==================================================================

!**************************************************************

        subroutine write_restart_gotm(iunit)
        use basin
        use levels
        use mod_gotm_aux
	use mod_restart
        implicit none
        integer iunit
	logical, parameter :: belem = .false.
        !integer l,k
        !write(iunit) ((numv_gotm(l,k),l=0,nlv),k=1,nkn)
        !write(iunit) ((nuhv_gotm(l,k),l=0,nlv),k=1,nkn)
        !write(iunit) ((tken_gotm(l,k),l=0,nlv),k=1,nkn)
        !write(iunit) ((eps_gotm(l,k),l=0,nlv),k=1,nkn)
        !write(iunit) ((rls_gotm(l,k),l=0,nlv),k=1,nkn)
	call restart_write_value(iunit,belem,numv_gotm)
	call restart_write_value(iunit,belem,nuhv_gotm)
	call restart_write_value(iunit,belem,tken_gotm)
	call restart_write_value(iunit,belem,eps_gotm)
	call restart_write_value(iunit,belem,rls_gotm)
        end

        subroutine skip_restart_gotm(iunit)
        implicit none
        integer iunit
        integer i
        do i=1,5
          read(iunit)
        end do
        end

        subroutine read_restart_gotm(iunit)
        use basin
        use levels
        use mod_gotm_aux
	use mod_restart
        implicit none
        integer iunit
	logical, parameter :: belem = .false.
        !integer l,k
	call restart_read_value(iunit,belem,numv_gotm)
	call restart_read_value(iunit,belem,nuhv_gotm)
	call restart_read_value(iunit,belem,tken_gotm)
	call restart_read_value(iunit,belem,eps_gotm)
	call restart_read_value(iunit,belem,rls_gotm)
        !read(iunit) ((numv_gotm(l,k),l=0,nlv),k=1,nkn)
        !read(iunit) ((nuhv_gotm(l,k),l=0,nlv),k=1,nkn)
        !read(iunit) ((tken_gotm(l,k),l=0,nlv),k=1,nkn)
        !read(iunit) ((eps_gotm(l,k),l=0,nlv),k=1,nkn)
        !read(iunit) ((rls_gotm(l,k),l=0,nlv),k=1,nkn)
        end

!**************************************************************

