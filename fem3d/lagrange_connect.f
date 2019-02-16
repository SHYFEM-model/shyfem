
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

c computes connectivity matrix
c
c revision log :
c
c 23.01.2012    ggu	new routine for release in point, connectivity
c
c*******************************************************************

	subroutine lagr_connect_continuous_points(brelease)

c continuous release from points

	implicit none

	logical brelease

	end

c*******************************************************************

	subroutine lagr_connect_count(ibdy,ie,ieorig,time)

	implicit none

	integer ibdy,ie,ieorig
	real time

	end

c******************************************************************

	subroutine lagr_connect_get_station(ie,ip_station,r_station)

	implicit none

	integer ie,ip_station
	real r_station

	ip_station = 0
	r_station = 0.5

	end

c******************************************************************

        subroutine lagr_connect_bitmap_init(ip)

        implicit none

        integer ip

        end

c******************************************************************

        subroutine lagr_connect_bitmap_copy(ifrom,ito)

        implicit none

        integer ifrom,ito

        end

c******************************************************************

