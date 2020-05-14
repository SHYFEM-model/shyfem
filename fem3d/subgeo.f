
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2001,2010-2011,2014-2015,2014-2015  Georg Umgiesser
!    Copyright (C) 2017-2019  Georg Umgiesser
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

c geometrical routines
c
c contents :
c
c function areat(x1,y1,x2,y2,x3,y3)
c function areapoly(n,x,y)
c
c subroutine centert(x1,y1,x2,y2,x3,y3,xc,yc)
c subroutine centerpoly(n,x,y,xc,yc)
c
c function left(x1,y1,x2,y2,x3,y3)
c function lefton(x1,y1,x2,y2,x3,y3)
c function collinear(x1,y1,x2,y2,x3,y3)
c
c function between(x1,y1,x2,y2,x3,y3)
c function betweens(x1,y1,x2,y2,x3,y3)
c function angle(x1,y1,x2,y2,x3,y3)
c
c subroutine dist_point_to_line(px,py,ax,ay,bx,by,qx,qy,d)
c
c function segsegint(x1,y1,x2,y2,x3,y3,x4,y4,xi,yi)
c function parallelint(x1,y1,x2,y2,x3,y3,x4,y4,xi,yi)
c
c function isconvex(n,x,y)
c function inconvex(n,x,y,x0,y0)
c function polyinpoly(n1,x1,y1,n2,x2,y2)
c
c function inpoly(n,x,y,x0,y0)		!shell for inpoly0 or inpoly1
c function inpoly1(n,x,y,x0,y0)		!fast algorithm
c function inpoly0(n,x,y,x0,y0)		!using winding number
c function winding(n,x,y,x0,y0)
c
c subroutine c2p(u,v,s,d)		converts cartesian to polar coordinates
c
c revision log :
c
c 18.03.1999	ggu	written from scratch (based on Joseph O-Rourke, 1998)
c 28.03.1999	ggu	bug fix : area2 -> areat (area was 2x)
c 28.03.1999	ggu	new routine polyinpoly
c 29.03.1999	ggu	segsegint is working with double precision internally
c 19.11.1999	ggu	new routines inpoly and winding
c 06.02.2001	ggu	new routine angle
c 22.03.2010	ggu	new routine inpoly() has been implemented
c 05.07.2010	ggu	new routine dist_point_to_line(), bug fix in inpoly1()
c 22.07.2010	ggu	changed VERS_6_1_9
c 14.07.2011	ggu	changed VERS_6_1_27
c 21.02.2014	ggu	new routine centert and centerpoly
c 07.03.2014	ggu	changed VERS_6_1_72
c 14.09.2015	ggu	changed VERS_7_2_2
c 07.10.2017	ggu	new routine c2p_ocean()
c 31.08.2018	ggu	changed VERS_7_5_49
c 25.10.2018	ggu	changed VERS_7_5_51
c 18.12.2018	ggu	changed VERS_7_5_52
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c
c**************************************************************
c**************************************************************
c**************************************************************

	function areat(x1,y1,x2,y2,x3,y3)

c computes area of triangle

	implicit none

	real areat
	real x1,y1,x2,y2,x3,y3

	areat = 0.5 * ( (x2-x1) * (y3-y1) - (x3-x1) * (y2-y1) )

	end

c**************************************************************

	function areapoly(n,x,y)

c computes area of polygon

	implicit none

	real areapoly
	integer n
	real x(n),y(n)

	integer i
	real sum
	real x1,y1,xo,yo,xn,yn

	real areat

	sum = 0.

	x1 = x(1)
	y1 = y(1)
	xn = x(2)
	yn = y(2)

	do i=3,n

	  xo = xn
	  yo = yn
	  xn = x(i)
	  yn = y(i)

	  sum = sum + areat(x1,y1,xo,yo,xn,yn)

	end do

	areapoly = sum

	end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine centert(x1,y1,x2,y2,x3,y3,xc,yc)

c computes center of triangle

	implicit none

	real x1,y1,x2,y2,x3,y3
	real xc,yc

	xc = (x1+x2+x3) / 3
	yc = (y1+y2+y3) / 3

	end

c**************************************************************

	subroutine centerpoly(n,x,y,xc,yc)

c computes center of polygon (center of gravity)

	implicit none

	integer n
	real x(n),y(n)
	real xc,yc

	integer i
	real x1,y1,xo,yo,xn,yn
	real area
	double precision areasum,xd,yd

        real areat

	xc = 0.
	yc = 0.

	if( n .le. 0 ) return

        areasum = 0.
	xd = 0.
	yd = 0.

        x1 = x(1)
        y1 = y(1)
        xn = x(2)
        yn = y(2)

        do i=3,n

          xo = xn
          yo = yn
          xn = x(i)
          yn = y(i)

          area = areat(x1,y1,xo,yo,xn,yn)

	  areasum = areasum + area
	  xd = xd + area*(x1+xo+xn)
	  yd = yd + area*(y1+yo+yn)

        end do

	xc = xd / (3.*areasum)
	yc = yd / (3.*areasum)

	end

c**************************************************************
c**************************************************************
c**************************************************************

	function left(x1,y1,x2,y2,x3,y3)

c left turn ?

	implicit none

	logical left
	real x1,y1,x2,y2,x3,y3

	real areat

	left = areat(x1,y1,x2,y2,x3,y3) .gt. 0.

	end

c**************************************************************

	function lefton(x1,y1,x2,y2,x3,y3)

c left turn or straight ?

	implicit none

	logical lefton
	real x1,y1,x2,y2,x3,y3

	real areat

	lefton = areat(x1,y1,x2,y2,x3,y3) .ge. 0.

	end

c**************************************************************

	function collinear(x1,y1,x2,y2,x3,y3)

c left turn ?

	implicit none

	logical collinear
	real x1,y1,x2,y2,x3,y3

	real areat

	collinear = areat(x1,y1,x2,y2,x3,y3) .eq. 0.

	end

c**************************************************************
c**************************************************************
c**************************************************************

	function between(x1,y1,x2,y2,x3,y3)

c is 3 between 1 and 2 ?

	implicit none

	logical between
	real x1,y1,x2,y2,x3,y3

	logical b1,b2
	logical collinear

	if( .not. collinear(x1,y1,x2,y2,x3,y3) ) then	!not on one line
	  between = .false.
	else if( x1 .ne. x2 ) then			!check on x coordinate
	  b1 = x1 .le. x3 .and. x3 .le. x2
	  b2 = x1 .ge. x3 .and. x3 .ge. x2
	  between = b1 .or. b2
	else						!check on y coordinate
	  b1 = y1 .le. y3 .and. y3 .le. y2
	  b2 = y1 .ge. y3 .and. y3 .ge. y2
	  between = b1 .or. b2
	end if

	end

c**************************************************************

	function betweens(x1,y1,x2,y2,x3,y3)

c is 3 between 1 and 2 ? (we already know that points are collinear)

	implicit none

	logical betweens
	real x1,y1,x2,y2,x3,y3

	logical b1,b2

	if( x1 .ne. x2 ) then				!check on x coordinate
	  b1 = x1 .le. x3 .and. x3 .le. x2
	  b2 = x1 .ge. x3 .and. x3 .ge. x2
	  betweens = b1 .or. b2
	else						!check on y coordinate
	  b1 = y1 .le. y3 .and. y3 .le. y2
	  b2 = y1 .ge. y3 .and. y3 .ge. y2
	  betweens = b1 .or. b2
	end if

	end

c**************************************************************

	function angle(x1,y1,x2,y2,x3,y3)

c angle between 1-2-3 in degrees 
c
c a < 180 => right turn
c a = 180 => straight
c a > 180 => left turn

	implicit none

	real angle
	real x1,y1,x2,y2,x3,y3

	double precision null,full,pi,rad
	parameter( null = 0. , full = 360. )
	parameter( pi = 3.14159265 , rad = 180. / pi )

	double precision norm,aux,ang
	double precision dx1,dy1,dx2,dy2

        dx1=x1-x2
        dy1=y1-y2
        dx2=x3-x2
        dy2=y3-y2
 
        norm = sqrt( (dx1*dx1+dy1*dy1) * (dx2*dx2+dy2*dy2) )
        aux = ( dx1*dx2 + dy1*dy2 ) / norm
 
        ang = acos( aux ) * rad
        if( dx1*dy2 - dy1*dx2 .lt. null ) ang = full - ang

	angle = ang

	end

c**************************************************************
c**************************************************************
c**************************************************************

        subroutine dist_point_to_line(px,py,ax,ay,bx,by,qx,qy,d)

c computes distance from point P to line given by A - B
c returns point Q on line closest to P and its distance d

        implicit none

        real px,py              ! single point P
        real ax,ay              ! point A defining line
        real bx,by              ! point B defining line
        real qx,qy              ! closest point Q on line to P (return)
        real d                  ! distance between Q and P (return)

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

c**************************************************************
c**************************************************************
c**************************************************************

	function segsegint(x1r,y1r,x2r,y2r,x3r,y3r,x4r,y4r,xi,yi)

c intersects two lines ( 1-2 and 3-4 ) and returns intersection in (xi,yi)
c
c return code :
c			0	no intersection
c			1	proper intersection
c			2	segments collinearly overlap
c			3	endpoint of one segment is on other
c
c test for .ne. 0 to find out on any intersection

	implicit none

	integer segsegint
	real x1r,y1r,x2r,y2r,x3r,y3r,x4r,y4r,xi,yi

	real*8 x1,y1,x2,y2,x3,y3,x4,y4
	real*8 denom,nums,numt
	real*8 s,t
	logical bsi,bti,bso,bto
	integer parallelint

	logical bdebug
	bdebug = .false.

c convert to double precision

	x1 = x1r
	x2 = x2r
	x3 = x3r
	x4 = x4r
	y1 = y1r
	y2 = y2r
	y3 = y3r
	y4 = y4r

c compute denominator

	denom = (x1-x2)*(y4-y3) + (x4-x3)*(y2-y1)

c if denominator segments are parallel -> handle seperately

	if( denom .eq. 0. ) then
	  segsegint = parallelint(x1r,y1r,x2r,y2r,x3r,y3r,x4r,y4r,xi,yi)
	  return
	end if

c compute intersection

	nums = x1*(y4-y3) + x3*(y1-y4) + x4*(y3-y1)
	s = nums / denom

	if( nums .eq. 0.0 .or. nums .eq. denom ) segsegint = 3
	if( s .eq. 0.0 .or. s .eq. 1.0 ) segsegint = 3

	numt = - ( x1*(y3-y2) + x2*(y1-y3) + x3*(y2-y1) )
	t = numt / denom

	if( numt .eq. 0.0 .or. numt .eq. denom ) segsegint = 3
	if( t .eq. 0.0 .or. t .eq. 1.0 ) segsegint = 3

c see if segments intersect

	bsi = 0.0 .lt. s .and. s .lt. 1.0
	bti = 0.0 .lt. t .and. t .lt. 1.0

	bso = 0.0 .gt. s .or. s .gt. 1.0
	bto = 0.0 .gt. t .or. t .gt. 1.0

	if( bsi .and. bti ) then
	  segsegint = 1
	else if( bso .or. bto ) then
	  segsegint = 0
	else
	  if( segsegint .ne. 3 ) then
		stop 'error stop segsegint: internal error'
	  end if
	end if

c compute point of intersection

	xi = x1 + s * ( x2 - x1 )
	yi = y1 + s * ( y2 - y1 )

	if( bdebug ) then
	  write(6,*) '---------- segsegint -------------'
	  write(6,*) x1,y1,x2,y2
	  write(6,*) x3,y3,x4,y4
	  write(6,*) denom,nums,numt,s,t
	  write(6,*) bsi,bti,bso,bto
	  write(6,*) segsegint,xi,yi
	  write(6,*) '----------------------------------'
	end if

	end

c**************************************************************

	function parallelint(x1,y1,x2,y2,x3,y3,x4,y4,xi,yi)

c intersects two parallel lines ( 1-2 and 3-4 ) and returns intersection
c
c return code :
c			0	no intersection
c			2	segments collinearly overlap

	implicit none

	integer parallelint
	real x1,y1,x2,y2,x3,y3,x4,y4,xi,yi

	logical collinear,betweens

	if( .not. collinear(x1,y1,x2,y2,x3,y3) ) then
	  parallelint = 0
	else if( betweens(x1,y1,x2,y2,x3,y3) ) then
	  xi = x3
	  yi = y3
	  parallelint = 2
	else if( betweens(x1,y1,x2,y2,x4,y4) ) then
	  xi = x4
	  yi = y4
	  parallelint = 2
	else if( betweens(x3,y3,x4,y4,x1,y1) ) then
	  xi = x1
	  yi = y1
	  parallelint = 2
	else if( betweens(x3,y3,x4,y4,x2,y2) ) then
	  xi = x2
	  yi = y2
	  parallelint = 2
	else
	  parallelint = 0
	end if

	end

c**************************************************************
c**************************************************************
c**************************************************************

	function isconvex(n,x,y)

c is polygon convex ?

	implicit none

	logical isconvex
	integer n
	real x(n),y(n)

	integer i
	real xl,yl,xm,ym,xn,yn
	logical lefton

	isconvex = .false.

	xm = x(n-1)
	ym = y(n-1)
	xn = x(n)
	yn = y(n)

	do i=1,n
	  xl = xm
	  yl = ym
	  xm = xn
	  ym = yn
	  xn = x(i)
	  yn = y(i)

	  if( .not. lefton(xl,yl,xm,ym,xn,yn) ) then
c	    write(6,*) 'isconvex: ',i-1,xm,ym
	    return
	  end if
	end do

	isconvex = .true.

	end

c**************************************************************

	function inconvex(n,x,y,x0,y0)

c checks if point (x0,y0) is in convex polygon (border is inside)

	implicit none

        logical inconvex !true if (x0,y0) is inside closed line (return)
	integer n        !total number of points in line
	real x(n),y(n)   !coordinates of line
	real x0,y0       !coordinates of point to check

	integer i
	real xm,ym,xn,yn
	logical lefton

	inconvex = .false.

	xn = x(n)
	yn = y(n)

	do i=1,n
	  xm = xn
	  ym = yn
	  xn = x(i)
	  yn = y(i)

	  if( .not. lefton(xm,ym,xn,yn,x0,y0) ) return
	end do

	inconvex = .true.

	end

c**************************************************************

	function polyinpoly(n1,x1,y1,n2,x2,y2)

c polygon 1 in convex polygon 2 ? (border is inside)

	implicit none

	logical polyinpoly
	integer n1,n2
	real x1(n1),y1(n1)
	real x2(n2),y2(n2)

	integer i
	logical inconvex

	polyinpoly = .false.

	do i=1,n1
	  if( .not. inconvex(n2,x2,y2,x1(i),y1(i)) ) return
	end do

	polyinpoly = .true.

	end

c**************************************************************
c**************************************************************
c**************************************************************

	function inpoly(n,x,y,x0,y0)

c shell for in-polygon check

	implicit none

        logical inpoly   !true if (x0,y0) is inside closed line (return)
	integer n        !total number of points in line
	real x(n),y(n)   !coordinates of line
	real x0,y0       !coordinates of point to check

	logical inpoly0	!old routine -> use classical winding number
	logical inpoly1	!new routine -> uses new algorithm from Dan Sunday

	if( n .gt. 0 ) then
	  inpoly = inpoly1(n,x,y,x0,y0)		!use new by default
	else if( n .lt. 0 ) then
	  inpoly = inpoly0(-n,x,y,x0,y0)
	else
	  inpoly = .false.
	end if

	end

c**************************************************************

	function inpoly1(n,x,y,x0,y0)

c checks if point (x0,y0) is in general polygon (border ?)
c
c fast algorithm from Dan Sunday
c
c http://softsurfer.com/Archive/algorithm_0103/algorithm_0103.htm

	implicit none

        logical inpoly1  !true if (x0,y0) is inside closed line (return)
	integer n        !total number of points in line
	real x(n),y(n)   !coordinates of line
	real x0,y0       !coordinates of point to check

	integer i,iw
	real xo,yo,xn,yn
	real area

	real areat

	iw = 0
	xo=x(n)
	yo=y(n)

	do i=1,n
	  xn=x(i)
	  yn=y(i)
	  area = areat(xo,yo,xn,yn,x0,y0)

	  if( yo .le. y0 ) then
	    if( yn .gt. y0 ) then
	      if( area .gt. 0. ) iw = iw + 1
	    end if
	  else
	    if( yn .le. y0 ) then
	      if( area .lt. 0. ) iw = iw - 1
	    end if
	  end if

	  xo = xn
	  yo = yn
	end do

	inpoly1 = iw .ne. 0

	end

c**************************************************************

	function inpoly0(n,x,y,x0,y0)

c checks if point (x0,y0) is in general polygon (border is inside)
c
c use classical winding number

	implicit none

        logical inpoly0  !true if (x0,y0) is inside closed line (return)
	integer n        !total number of points in line
	real x(n),y(n)   !coordinates of line
	real x0,y0       !coordinates of point to check

	integer i
	real xm,ym,xn,yn
	real xmin,xmax,ymin,ymax
	logical bdebug

	logical isconvex,inconvex
	integer winding

	bdebug = .false.

	if(bdebug) write(6,*) 'debug is set in inpoly...'

c-------------------------------------------------------
c initialize inpoly
c-------------------------------------------------------

	inpoly0 = .false.

	if( n .le. 0 ) return

c-------------------------------------------------------
c maybe polygon is convex ?
c-------------------------------------------------------

	if( isconvex(n,x,y) ) then
	  inpoly0 = inconvex(n,x,y,x0,y0)
	  return
	end if

c-------------------------------------------------------
c no it is not -> get min/max and make first decisions
c-------------------------------------------------------

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

	if( bdebug ) then
	  write(6,*) n,x0,y0
	  write(6,*) xmin,xmax,ymin,ymax
	end if

	if( x0 .lt. xmin .or. x0 .gt. xmax ) return
	if( y0 .lt. ymin .or. y0 .gt. ymax ) return

c-------------------------------------------------------
c not yet decided -> proceed with heavy algorithm
c-------------------------------------------------------

	if( winding(n,x,y,x0,y0) .eq. 0 ) return

c-------------------------------------------------------
c ok, it is inside
c-------------------------------------------------------

	inpoly0 = .true.

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	end

c**************************************************************

	function winding(n,x,y,x0,y0)

c compute winding number for x0/y0 (integer)

	implicit none

	integer winding
	integer n
	real x(n),y(n)
	real x0,y0

	integer i
	logical bdebug
	real xm,ym,xn,yn
	double precision dxn,dyn,dxm,dym
	double precision scal,cros,modu,arg
	double precision totangle,angle
	real pi2

	bdebug = .false.

	pi2 = 2. * 4. * atan(1.)
	totangle = 0.

	xn=x(n)
	yn=y(n)

	do i=1,n
	  xm=x(i)
	  ym=y(i)

	  dxn = xn-x0
	  dyn = yn-y0
	  dxm = xm-x0
	  dym = ym-y0

	  scal = dxn * dxm + dyn * dym
	  cros = dxn * dym - dxm * dyn
	  modu = sqrt( (dxn**2 + dyn**2) * (dxm**2 + dym**2) )
	  arg = scal / modu

	  if( arg .gt.  1. ) arg =  1.
	  if( arg .lt. -1. ) arg = -1.

	  angle = acos( arg )	! arg in [-1,+1], acos in [0,pi]

	  if( cros .lt. 0. ) angle = -angle

	  totangle = totangle + angle

	  xn = xm
	  yn = ym
	end do

	winding = nint( totangle / pi2 )

	if( bdebug ) then
	  write(6,*) 'winding... ',totangle,winding
	end if

	end

c**************************************************************
c**************************************************************
c**************************************************************

	subroutine c2p_ocean(u,v,s,d)

c converts cartesian to polar coordinates
c oceanographic convention: 0 -> current to north

	real u,v,s,d

	call c2p(u,v,s,d)

        d = d + 180.
        if( d > 360. ) d = d - 360.

	end

c**************************************************************

	subroutine c2p(u,v,s,d)

c converts cartesian to polar coordinates
c
c meteorological convention: 0 -> wind from north
c if oceanographic convention is needed: d=d+180; if(d>360) d=d-360

	implicit none

	real u,v
	real s,d

	real rad,a

c----------------------------------------------------------------
c Pi
c----------------------------------------------------------------

	rad = 45. / atan (1.)

c----------------------------------------------------------------
c compute speed
c----------------------------------------------------------------

	s = sqrt( u**2 + v**2 )

c----------------------------------------------------------------
c compute mathematical angle
c----------------------------------------------------------------

	if( u .eq. 0. ) then
	  if( v .gt. 0. ) then
	    a = 90.
	  else
	    a = -90.
	  end if
	else
	  a = rad * atan(v/u)
	end if

	if( u .lt. 0. ) a = a + 180.

c----------------------------------------------------------------
c convert to wind directions
c----------------------------------------------------------------

	d = 270 - a

	if( d .gt. 360. ) d = d - 360.
	if( d .lt. 0.   ) d = d + 360.

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c**************************************************************
c**************************************************************
c**************************************************************

