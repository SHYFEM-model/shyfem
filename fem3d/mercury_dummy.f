
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016-2020  Georg Umgiesser
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

! mercury dummy routines
!
! contents :
!
! revision log :
!
! 15.05.2016	ggu	started mercury from bio3d
! 25.05.2016	ggu	changed VERS_7_5_10
! 07.06.2016	ggu	changed VERS_7_5_12
! 10.06.2016	ggu	changed VERS_7_5_13
! 17.06.2016	ggu	changed VERS_7_5_15
! 13.04.2017	ggu	changed VERS_7_5_25
! 09.05.2017	ggu	changed VERS_7_5_26
! 13.06.2017	ggu	changed VERS_7_5_29
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 09.01.2020	ggu	dummy routine written
! 09.03.2020	ggu	dummy restart routines
!
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!
! notes :
!
! State variables used: (mercury) -> Donata, please adjourn
!
! State variables used: (mer
! Hg0           81      1
! Hg2           82      2
! Hg3           83      3
!
!********************************************************************

        subroutine mercury_module

! general interface to mercury module

        implicit none

        integer, save :: icall = 0

	integer imerc
	real getpar

        if( icall .le. -1 ) return

        if( icall .eq. 0 ) then
          imerc = nint(getpar('imerc'))
          if( imerc .le. 0 ) icall = -1
          if( icall .le. -1 ) return
          stop 'error stop mercury_init: mercury module is not linked'
        end if

        end

!********************************************************************

	subroutine mercury_init

! initializes mercury routines

	implicit none

	end

!*************************************************************

	subroutine write_restart_mercury(iunit)
	implicit none
	integer iunit
	end

	subroutine skip_restart_mercury(iunit)
	implicit none
	integer iunit
	end

	subroutine read_restart_mercury(iunit)
	implicit none
	integer iunit
	end

!*************************************************************

