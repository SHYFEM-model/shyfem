c
c $Id: lagrange_util.f,v 1.1 2009-02-13 17:22:44 georg Exp $
c
c general utilities for lagrangian model
c
c revision log :
c
c 05.02.2009    ggu     copied from other files
c 28.03.2014    ggu     compute total length of open boundary
c
c*******************************************************************

        function dist_node(k1,k2)

c returns distance between two nodes

	use basin

        implicit none

        real dist_node
        integer k1,k2

	include 'param.h'

        real x1,y1,x2,y2,dx,dy

        x1 = xgv(k1)
        y1 = ygv(k1)
        x2 = xgv(k2)
        y2 = ygv(k2)

        dx = x1 - x2
        dy = y1 - y2

        dist_node = sqrt( dx*dx + dy*dy )

        end

c*******************************************************************

	subroutine dist_total(ibc,totdist)

c returns total length of open boundary

	implicit none

	integer ibc
	real totdist

	integer nk,i,k1,k2
	real dxy

	integer nkbnds,kbnds
	real dist_node

	nk = nkbnds(ibc)

        totdist = 0.
        do i=2,nk
          k1 = kbnds(ibc,i-1)
          k2 = kbnds(ibc,i)
          dxy = dist_node(k1,k2)
          totdist = totdist + dxy
        end do

	end

c*******************************************************************

	subroutine basin_center(xm,ym)

c returns center of gravity of total basin

	use basin

	implicit none

	real xm,ym

	include 'param.h'

	call xy_center(nkn,xgv,ygv,xm,ym)

	end
	
c*******************************************************************

	subroutine compute_total_area(area)

	use evgeom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real area



	integer ie
	double precision a,tot_area

	tot_area = 0.

	do ie=1,nel
	  a = ev(10,ie)
	  tot_area = tot_area + a
	end do

	area = 12. * tot_area

	end

c*******************************************************************

