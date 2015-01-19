c
c $Id: basinf.f,v 1.19 2010-03-08 17:46:45 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	completely restructured
c 04.06.1999	ggu	new statistics are computed
c 28.11.2005	ggu	new call to makehkv
c 31.05.2007	ggu	added area and volume frequency curve
c 24.08.2007	ggu	added new routine write_grd_from_bas
c 06.04.2009    ggu     read param.h
c 12.06.2009    ggu     areatr in double precision - new algorithm
c 01.03.2010    ggu     new routine basqual() to compute grid quality
c
c****************************************************************

        program wincrea

c creates special wind file

	implicit none

	include 'param.h'
	include 'basin.h'

	real haux(nkndim)
	include 'depth.h'
	include 'evmain.h'

	logical bnode,belem
	integer it,it0,idt,it_spinup,nt
	integer i
	real dn,ds,speedn,speeds
	real speed,dir
	real wx(nkndim), wy(nkndim)

	integer iapini

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        if( iapini(1,nkndim,neldim,0) .le. 0 ) stop

	call bas_info

	call set_ev

c-----------------------------------------------------------------
c interpolation in time -> please change parameters here
c-----------------------------------------------------------------

	it_spinup = 2000000	!time for spin up
	it0 = 0			!time where wind changes
	idt = 3600		!time step for rotating wind
	nt = 6			!number of time steps to rotate

	dn = 45.		!direction where wind comes from
	ds = -45.		!n: north  s: south
	speedn = 15.		!speed of wind
	speeds = 12.		!n: north  s: south

	it = it0 - it_spinup
	call make_wind(it,speeds,ds,speeds,ds,wx,wy)
	it = it0
	call make_wind(it,speeds,ds,speeds,ds,wx,wy)

	do i=1,nt
	  it = it + idt
	  speed = ( i*speedn + (nt-i)*speeds ) / nt
	  dir = ( i*dn + (nt-i)*ds ) / nt
	  call make_wind(it,speed,dir,speeds,ds,wx,wy)
	end do

	it = it + it_spinup
	call make_wind(it,speedn,dn,speeds,ds,wx,wy)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

	subroutine make_wind(it,speedn,dn,speeds,ds,wx,wy)

	implicit none

	include 'param.h'

	integer it
	real speeds,ds
	real speedn,dn
	real wx(1), wy(1)


	include 'basin.h'

	logical left,lefton

	real rn1x,rn1y
	real rn2x,rn2y
	real rs1x,rs1y
	real rs2x,rs2y

	integer k,i
	real x,y,qx,qy
	real d1,d2,dd
	real speed,dir

c------------------------------------------------------
c define north and south line
c------------------------------------------------------

c 935 0 -279348.000000 343894.000000 148.862305
c 2549 0 -143277.000000 432092.000000 10.000000
c 415 0 -187531.000000 252456.000000 375.171387
c 689 0 -67511.601562 346310.000000 142.882690

c northern line

	rn1x = -279348.000000
	rn1y = 343894.000000
	rn2x = -143277.000000
	rn2y = 432092.000000

c southern line

	rs1x = -187531.000000
	rs1y = 252456.000000
	rs2x = -67511.601562
	rs2y = 346310.000000

c------------------------------------------------------
c do not change anything beyond this point
c------------------------------------------------------

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)

	  if( lefton(rn1x,rn1y,rn2x,rn2y,x,y) ) then	! (x,y) is to the north
	    call convert_wind(speedn,dn,wx(k),wy(k))
	  else if( .not. left(rs1x,rs1y,rs2x,rs2y,x,y) ) then	! south
	    call convert_wind(speeds,ds,wx(k),wy(k))
	  else						! intermediate band
	    call dist_point_to_line(x,y,rn1x,rn1y,rn2x,rn2y,qx,qy,d1)
	    call dist_point_to_line(x,y,rs1x,rs1y,rs2x,rs2y,qx,qy,d2)
	    dd = d1 + d2
	    speed = (d1*speeds + d2*speedn)/dd
	    dir = (d1*ds + d2*dn)/dd
	    call convert_wind(speed,dir,wx(k),wy(k))
	  end if
	end do

	write(1) it,nkn
	write(1) (wx(i),wy(i),i=1,nkn)

	end

c*******************************************************************

	subroutine convert_wind(speed,dir,wx,wy)

	implicit none

	real speed,dir
	real wx,wy

	real pi,rad
	parameter ( pi=3.14159 , rad = pi/180. )

	real d

	d = rad * (dir+180.)
	wx = speed * cos(d)
	wy = speed * sin(d)

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine dist_point_to_line0(px,py,ax,ay,bx,by,qx,qy,d)

c computes distance from point (px,py) to line given by (ax,ay) - (bx,by)
c returns point (qx,qy) on line closest to p and its distance d

	implicit none

	real px,py		! single point P
	real ax,ay		! point A defining line
	real bx,by		! point B defining line
	real qx,qy		! closest point Q on line to P (return)
	real d			! distance between Q and P (return)

	real rx,ry,sx,sy
	real lambda

	rx = bx - ax
	ry = by - ay
	sx = px - ax
	sy = py - ay

	lambda = (rx*sx + ry*sy) / (rx*rx + ry*ry)

	qx = ax + lambda * rx
	qy = ay + lambda * ry

	d = sqrt( (qx-px)**2 + (qy-py)**2 )

	end

c*******************************************************************

