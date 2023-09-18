!
! $Id: lagrange_cont.f,v 1.5 2009-09-14 08:20:57 georg Exp $
!
! computes connectivity matrix
!
! revision log :
!
! 23.01.2012    ggu	new routine for release in point, connectivity
!
!*******************************************************************
!------------------------------------------------------------------------
        module lagrange_connect
!------------------------------------------------------------------------
        contains
!------------------------------------------------------------------------

	subroutine lagr_connect_continuous_points(brelease)

! continuous release from points

	implicit none

	logical brelease

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine lagr_connect_count(ibdy,ie,ieorig,time,ic)

	implicit none

	integer ibdy,ie,ieorig
	double precision time
	integer ic

	end

!*******************************************************************

!******************************************************************

	subroutine lagr_connect_get_station(ie,ip_station,r_station)

	implicit none

	integer ie,ip_station
	double precision r_station

	  ip_station = 0
	  r_station = 0.5

	end

!******************************************************************

!------------------------------------------------------------------------
        end module lagrange_connect
!------------------------------------------------------------------------
