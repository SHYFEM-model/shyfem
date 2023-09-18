!
! $Id: lagrange_util.f,v 1.1 2009-02-13 17:22:44 georg Exp $
!
! general utilities for lagrangian model
!
! revision log :
!
! 05.02.2009    ggu     copied from other files
! 28.03.2014    ggu     compute total length of open boundary
!
!*******************************************************************
!---------------------------------------------------------------------
        module  lagrange_util
!---------------------------------------------------------------------
        contains
!---------------------------------------------------------------------

        function dist_node(k1,k2)

! returns distance between two nodes

	use basin

        implicit none

        double precision dist_node
        integer k1,k2

	include 'param.h'

        double precision x1,y1,x2,y2,dx,dy

        x1 = xgv(k1)
        y1 = ygv(k1)
        x2 = xgv(k2)
        y2 = ygv(k2)

        dx = x1 - x2
        dy = y1 - y2

        dist_node = sqrt( dx*dx + dy*dy )

        end

!*******************************************************************

	subroutine dist_total(ibc,totdist)

! returns total length of open boundary

        use bnd_admin

	implicit none

	integer ibc
	double precision totdist

	integer nk,i,k1,k2
	double precision dxy

	nk = nkbnds(ibc)

        totdist = 0.
        do i=2,nk
          k1 = kbnds(ibc,i-1)
          k2 = kbnds(ibc,i)
          dxy = dist_node(k1,k2)
          totdist = totdist + dxy
        end do

	end

!*******************************************************************

	subroutine xy_minmax(n,x,y,xmin,xmax,ymin,ymax)

! returns min/max coordinates of polygon

	implicit none

	integer n
	double precision x(1),y(1)
	double precision xmin,xmax
	double precision ymin,ymax

	integer i

	xmin = x(1)
	xmax = x(1)
	ymin = y(1)
	ymax = y(1)

	do i=2,n
	  xmin = min(xmin,x(i))
	  xmax = max(xmax,x(i))
	  ymin = min(ymin,y(i))
	  ymax = max(ymax,y(i))
	end do

	end

!*******************************************************************

	subroutine basin_center(xm,ym)

! returns center of gravity of total basin

	use basin

	implicit none

	double precision xm,ym

	include 'param.h'

	call xy_center(nkn,xgv,ygv,xm,ym)

	end
	
!*******************************************************************

	subroutine xy_center(n,x,y,xm,ym)

! returns center of gravity of polygon

	implicit none

	integer n
	double precision x(1),y(1)
	double precision xm,ym

	integer i

	xm = 0.
	ym = 0.

	do i=1,n
	  xm = xm + x(i)
	  ym = ym + y(i)
	end do

	xm = xm / n
	ym = ym / n

	end

!*******************************************************************

	subroutine xy_and_minmax(ie,x,y,xmin,xmax,ymin,ymax)

! returns x,y and min/max coordinates of vertices of element ie

	use basin

	implicit none

	integer ie
	double precision x(3),y(3)
	double precision xmin,xmax
	double precision ymin,ymax

	include 'param.h'

	integer ii,k

	do ii=1,3
          k=nen3v(ii,ie)
          x(ii)=xgv(k)
          y(ii)=ygv(k)
        end do

        xmin=min(x(1),x(2),x(3))
        xmax=max(x(1),x(2),x(3))
        ymin=min(y(1),y(2),y(3))
        ymax=max(y(1),y(2),y(3))

	end

!*******************************************************************

	subroutine compute_total_area(area)

	use evgeom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	double precision area



	integer ie
	double precision a,tot_area

	tot_area = 0.

	do ie=1,nel
	  a = ev(10,ie)
	  tot_area = tot_area + a
	end do

	area = 12. * tot_area

	end

!*******************************************************************

!---------------------------------------------------------------------
        end module lagrange_util
!---------------------------------------------------------------------
