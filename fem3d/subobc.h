
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013,2019  Georg Umgiesser
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
! 13.06.2013	ggu	changed VERS_6_1_65
! 16.02.2019	ggu	changed VERS_7_5_60

	integer ip_iunit		!unit number, -1 for not open
	integer ip_nintp		!interpolation, 2: linear, 4: cubic
	integer ip_nvar			!number of variables stored
	integer ip_nsize		!total size for one variable
	integer ip_ndata		!total size for all variables
	integer ip_ndim			!dimension of array
	integer ip_nextra		!extra header information
	integer ip_ires			!pointer to where results are stored
	integer ip_nspace		!space needed for all info in array
	integer ip_np			!number of horizontal points
	integer ip_lmax			!max number of levels
	integer ip_iformat		!is formatted?

c np = 0, nsize = 0	-> time series
c np > 0, nsize > 0	-> fem file format

	parameter( ip_iunit   =  1 )
	parameter( ip_nintp   =  2 )
	parameter( ip_nvar    =  3 )
	parameter( ip_nsize   =  4 )
	parameter( ip_ndata   =  5 )
	parameter( ip_ndim    =  6 )
	parameter( ip_nextra  =  7 )
	parameter( ip_ires    =  8 )
	parameter( ip_nspace  =  9 )
	parameter( ip_np      = 10 )
	parameter( ip_lmax    = 11 )
	parameter( ip_iformat = 12 )

        integer nextra
        parameter( nextra = 13 )

        real rguard
        parameter( rguard = 1.234543e+20 )

