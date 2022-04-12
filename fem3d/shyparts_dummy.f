
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
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
! 19.05.2020	ccf	started from scratch
! 12.04.2022	ggu	adapted
!
!****************************************************************

        subroutine do_partition(nkn,nel,nen3v,nparts,npart,epart)

! shyparts dummy routine

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer nparts
        integer npart(nkn)
        integer epart(nel)

        write(6,*)' For automatic partitioning of a grid install' 
        write(6,*)' one of the following libraries and set the'
        write(6,*)' parameters PARTS and PARTSDIR in the'
        write(6,*)' Rules.make configuration file' 
        write(6,*)'   - METIS'
        write(6,*)' Then recompile: "make fem"'

	stop 'error stop do_partition: no metis available'

	end

!*******************************************************************

	subroutine check_partition(npart,epart,ierr1,ierr2)

	use basin

	implicit none

        integer npart(nkn)
        integer epart(nel)
	integer ierr1,ierr2

	end

!*******************************************************************

