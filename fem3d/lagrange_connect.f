c
c $Id: lagrange_cont.f,v 1.5 2009-09-14 08:20:57 georg Exp $
c
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

