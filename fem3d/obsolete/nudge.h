
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
! 12.11.2013	ggu	changed VERS_6_1_69
! 12.11.2013	ggu	changed VERS_6_1_69a
! 16.02.2019	ggu	changed VERS_7_5_60

	integer ndgdim			!max size of variables
	parameter (ndgdim = 50)

	integer ndgdatdim		!max size for BC file
	parameter (ndgdatdim = 10*ndgdim)

	real andg_data(ndgdatdim)	!data of observations
	integer ndg_nodelist(ndgdim)	!nodes of obs
	integer ndg_use(ndgdim)		!use observations
	character*40 ndg_names(ndgdim)	!name of stations

	real andg_dist(nkndim)		!distance to be used in algorithm
	real andg_weight(nkndim)	!weight to be used in algorithm
	real andg_obs(nkndim)		!observations to be used in algorithm

	integer ndg_nodes(nkndim)	!nodes of influence
	integer ndg_area(nkndim)	!area of influence

	integer nvar
	real tramp

	common /i_nudge/ nvar,ndg_nodes,ndg_area,ndg_nodelist,ndg_use
	save /i_nudge/

	common /r_nudge/ tramp,andg_data,andg_dist,andg_weight,andg_obs
	save /r_nudge/

	common /c_nudge/ ndg_names
	save /c_nudge/

