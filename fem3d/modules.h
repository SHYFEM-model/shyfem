
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2019  Georg Umgiesser
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

c------------------------------------------------------------
c parameter definition for module framework
c------------------------------------------------------------

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60

c------------------------------------------------------------


	integer M_INIT		!initialization stand-alone
	integer M_READ		!read in section
	integer M_CHECK		!check data read in
	integer M_SETUP		!set up data structures -> depend on externals
	integer M_PRINT		!print out details for log file
	integer M_TEST		!test output for debug
	integer M_BEFOR		!do at beginning of time loop
	integer M_AFTER		!do at end of time loop

	parameter(
     +			  M_INIT	= 1
     +			, M_READ	= 2
     +			, M_CHECK	= 3
     +			, M_SETUP	= 4
     +			, M_PRINT	= 5
     +			, M_TEST	= 6
     +			, M_BEFOR	= 7
     +			, M_AFTER	= 8
     +		 )

