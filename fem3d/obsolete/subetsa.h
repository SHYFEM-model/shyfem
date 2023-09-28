
!--------------------------------------------------------------------------
!
!    Copyright (C) 2014,2019  Georg Umgiesser
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
! 28.01.2014	ggu	changed VERS_6_1_71
! 16.02.2019	ggu	changed VERS_7_5_60

	integer nets			!total number of extra nodes
	integer nkets(nexdim)		!node numbers
	real xets(nexdim)		!x-coordinates of nodes
	real yets(nexdim)		!y-coordinates of nodes
	character*80 chets(nexdim)	!description of nodes

	integer ilets(nexdim)		!layers for nodes (needed for header)
	real hets(nexdim)		!depth for nodes (needed for header)

	common /ets_int/ nets,nkets,ilets
	common /ets_real/ xets,yets,hets
	common /ets_char/ chets

	save /ets_int/,/ets_real/,/ets_char/

