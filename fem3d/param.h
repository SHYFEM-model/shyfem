
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010-2015,2019  Georg Umgiesser
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

!------------------------------------------------------------------------
! param.h - parameter file for SHYFEM
!------------------------------------------------------------------------

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 26.03.2010	ggu	changed v6.1.2
! 14.04.2010	ggu	changed v6.1.4
! 26.04.2010	ggu	changed VERS_6_1_6
! 22.07.2010	ggu	changed VERS_6_1_9
! 29.09.2010	ggu	changed VERS_6_1_12
! 20.12.2010	ggu	changed VERS_6_1_16
! 27.01.2011	ggu	changed VERS_6_1_17
! 17.02.2011	ggu	changed VERS_6_1_18
! 18.02.2011	ggu	changed VERS_6_1_19
! 01.03.2011	ggu	changed VERS_6_1_20
! 23.03.2011	ggu	changed VERS_6_1_21
! 14.04.2011	ggu	changed VERS_6_1_22
! 31.05.2011	ggu	changed VERS_6_1_23
! 31.05.2011	ggu	changed VERS_6_1_24
! 07.06.2011	ggu	changed VERS_6_1_25
! 08.06.2011	ggu	changed VERS_6_1_26
! 14.07.2011	ggu	changed VERS_6_1_27
! 15.07.2011	ggu	changed VERS_6_1_28
! 19.08.2011	ggu	changed VERS_6_1_29
! 18.10.2011	ggu	changed VERS_6_1_33
! 24.01.2012	ggu	changed VERS_6_1_41
! 09.03.2012	ggu	changed VERS_6_1_47
! 16.03.2012	ggu	changed VERS_6_1_48
! 19.03.2012	ggu	changed VERS_6_1_49
! 21.03.2012	ggu	changed VERS_6_1_50
! 01.06.2012	ggu	changed VERS_6_1_53
! 26.06.2012	ggu	changed VERS_6_1_55
! 29.08.2012	ggu	changed VERS_6_1_56
! 12.09.2012	ggu	changed VERS_6_1_57
! 08.10.2012	ggu	changed VERS_6_1_58
! 05.11.2012	ggu	changed VERS_6_1_60
! 19.11.2012	ggu	changed VERS_6_1_61
! 25.01.2013	ggu	changed VERS_6_1_62
! 03.05.2013	ggu	changed VERS_6_1_63
! 10.05.2013	ggu	changed VERS_6_1_64
! 13.06.2013	ggu	changed VERS_6_1_65
! 25.10.2013	ggu	changed VERS_6_1_68
! 12.11.2013	ggu	changed VERS_6_1_69
! 05.12.2013	ggu	changed VERS_6_1_70
! 27.03.2014	ggu	changed VERS_6_1_73
! 30.05.2014	ggu	changed VERS_6_1_76
! 18.06.2014	ggu	changed VERS_6_1_77
! 07.07.2014	ggu	changed VERS_6_1_79
! 07.11.2014	ggu	changed VERS_7_0_6
! 26.11.2014	ggu	changed VERS_7_0_7
! 01.04.2015	ggu	changed VERS_7_1_7
! 23.04.2015	ggu	changed VERS_7_1_8
! 30.04.2015	ggu	changed VERS_7_1_9
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_53
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 31.07.2015	ggu	changed VERS_7_1_84
! 16.02.2019	ggu	changed VERS_7_5_60

!------------------------------------------------------------------------


	integer nkndim			!maximum number of nodes
	integer neldim			!maximum number of elements
	integer nlvdim			!maximum number of vertical levels

	integer mbwdim			!maximum bandwidth
	integer ngrdim			!maximum grade of nodes

	integer nbcdim			!maximum number of open boundaries
	integer nrbdim			!maximum number of open boundary nodes
	integer nb3dim			!maximum storage for boundary info

	integer nardim			!maximum area code
	integer nexdim			!maximum number of extra points
	integer nfxdim			!maximum number of flux points
	integer ncsdim			!maximum concentration variables

	integer nbdydim			!maximum particles for lagrange model

	parameter ( nkndim = 1 )
	parameter ( neldim = 1 )
	parameter ( nlvdim = 1 )

	parameter ( mbwdim = 1 )
	parameter ( ngrdim = 1 )

	parameter ( nbcdim = 1 )
	parameter ( nrbdim = 1 )
	parameter ( nb3dim = 1 )

	parameter ( nardim = 1 )
	parameter ( nexdim = 1 )
	parameter ( nfxdim = 1 )
	parameter ( ncsdim = 1 )

	parameter ( nbdydim = 100000 )

!------------------------------------------------------------------------
! do not change anything beyond this line
!------------------------------------------------------------------------

	integer nlkdim			!dimension for side index
        parameter (nlkdim=3*neldim+2*nkndim)

!------------------------------------------------------------------------
! end of parameter file
!------------------------------------------------------------------------

