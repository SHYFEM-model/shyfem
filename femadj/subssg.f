c
c $Id: subssg.f,v 1.2 2007-03-20 13:19:42 georg Exp $
c
c geometrical routines
c
c contents :
c
c function vproj(u1,v1,u2,v2)				projection of vector
c function icut(x1,y1,x2,y2,x3,y3,x4,y4,xout,yout)	intersects two lines
c function intria(x,y,xp,yp)				point in triangle or not
c subroutine mirr(x1,y1,x2,y2,xp,yp)			reflection of point
c subroutine intrfk(x,y,xp,yp,fk)			interpolates point (fk)
c function rintrz(x,y,xp,yp,zp)				interpolates point (zp)
c function icutri(x,y,x1,y1,x2,y2,ds)			intersects line & tri.
c
c***************************************************************
c
	function vproj(u1,v1,u2,v2)
c
c projects a vector (u1,v1) onto another vector (u2,v2)
c
c (u1,v1)	vector to project
c (u2,v2)	vector on whom to project (u1,v1)
c vproj		return value, giving the length of the projected
c		...vector in terms of the second one. The
c		...projected vector can be obtained multiplying
c		...the second vector by vproj :
c		...   (up,vp) = vproj * (u2,v2)
c
c formulas :	u is to be projected onto v
c		... |uproj| = |u| * cos(alfa)    and
c		...  uproj  = |uproj| * v / |v|    with
c		... cos(alfa) = (u * v) / (|u| * |v|)
c
	real vproj

	uv2 = u2*u2 + v2*v2
c
	if(uv2.lt.1.e-6) then
		write(6,*) 'Cannot project on nill vector'
		write(6,*) 'Length of second vector : ',uv2
		vproj = 0.	!$$ALPHA
		return
	end if
c
	vproj = (u1*u2 + v1*v2) / uv2
c
	return
	end
c
c*************************************************************
c
	function icut(x1,y1,x2,y2,x3,y3,x4,y4,xout,yout)
c
c intersects two lines and gives back intersection point
c
c x1,y1		starting point of first line
c x2,y2		ending point of first line
c x3,y3		starting point of second line
c x4,y4		ending point od second line
c xout,yout	point of intersection (if any)
c icut		return code :
c		-1 : lines parallel (but not equal)
c		-2 : lines equal
c		0  : lines intersect
c		1  : intersection inside 1 but outside 2
c		2  : intersection inside 2 but outside 1
c		3  : intersection outside of both lines
c		...(not-negative return code returns intersection
c		...point in xout,yout)
c
c formulas :	linear system
c				a*r1 + b*r2 = e
c				c*r1 + d*r2 = f
c
	parameter (eps1=1.e-6,eps2=1.e-6)
c	parameter (rmin=0.,rmax=1.)
	parameter (rmin=0.-eps1,rmax=1.+eps1)
c
	a=x2-x1
	b=x3-x4
	c=y2-y1
	d=y3-y4
	e=x3-x1
	f=y3-y1
c
	det = a*d-c*b
c
	if(abs(det).lt.eps2) then	!lines parallel
		if(abs(a*f-c*e).lt.eps2) then	!lines equivalent
			icut=-2
		else
			icut=-1
		end if
		return
	else				!compute intersection
		r1=(e*d-f*b)/det
		r2=(a*f-c*e)/det
c
		if(r1.ge.rmin.and.r1.le.rmax) then
			if(r2.ge.rmin.and.r2.le.rmax) then
				icut=0
			else
				icut=1
			end if
			xout=x1+r1*(x2-x1)
			yout=y1+r1*(y2-y1)
		else
			if(r2.ge.rmin.and.r2.le.rmax) then
				icut=2
				xout=x3+r2*(x4-x3)
				yout=y3+r2*(y4-y3)
			else
				icut=3
				if(abs(r1).lt.abs(r2)) then
					xout=x1+r1*(x2-x1)
					yout=y1+r1*(y2-y1)
				else
					xout=x3+r2*(x4-x3)
					yout=y3+r2*(y4-y3)
				end if
			end if
		end if
	end if
c
	return
	end
c
c***********************************************************
c
	function intria(x,y,xp,yp)
c
c point in triangle or not
c
c x,y		array of coordinates of vertices of triangle
c xp,yp		coordinates of point
c intria	1: point is in triangle  0: point outside (return value)
c
	dimension x(3),y(3)
	data eps /1.e-7/
c
	xs=0.
	ys=0.
	do i=1,3
	   xs=xs+x(i)
	   ys=ys+y(i)
	end do
	xs=xs/3.
	ys=ys/3.
c
	intria=0
c
	do k1=1,3
	   k2=mod(k1,3)+1
	   x12=x(k1)-x(k2)
	   y12=y(k1)-y(k2)
	   det=(xs-xp)*y12-(ys-yp)*x12
	   if(abs(det).ge.eps) then
		detlam=(x(k1)-xp)*y12-(y(k1)-yp)*x12
		rlamb=detlam/det
		if(rlamb.gt.0..and.rlamb.lt.1.) return	!outside
	   end if
	end do
c
	intria=1	!inside
c
	return
	end
c
c****************************************************************
c
	function mirr(x1,y1,x2,y2,xp,yp)
c
c reflection of one point on a line
c
c x1,y1,x2,y2	points defining line
c xp,yp		on entry : point to be reflected
c 		on return : reflected point
c mirr		0 : ok     -1 : error
c
	u=x2-x1
	v=y2-y1
c
	uv = u*u + v*v
c
	if(uv.lt.1.e-6) then
		write(6,*) 'Cannot reflect on nill vector'
		write(6,*) 'Length of second vector : ',uv
		mirr=-1
		return
	end if
c
	r=(u*(xp-x1)+v*(yp-y1))/uv
c
	xp=2.*(r*u+x1)-xp
	yp=2.*(r*v+y1)-yp
c
	mirr=0
c
	return
	end
c
c*********************************************************************
c
	subroutine intrfk(x,y,xp,yp,fk)
c
c interpolates point in triangle and returns fk
c
c x,y		vertices of triangle
c xp,yp		coordinates of point
c fk		form functions (return value)
c
c if c(x(i),y(i)),i=1,3 are the values at vertices,
c ... then c(xp,yp) can be obtained by :
c
c		c(xp,yp) = c(1)*fk(1) + c(2)*fk(2) + c(3)*fk(3)
c
	dimension x(3),y(3),fk(3)
c
	f=0.
	do i=1,3
	   ii=mod(i,3)+1
	   iii=mod(ii,3)+1
	   a=x(ii)*y(iii)-x(iii)*y(ii)
	   b=y(ii)-y(iii)
	   c=x(iii)-x(ii)
	   fk(i) = a + xp*b + yp*c
	   f=f+a
	end do
c
	do i=1,3
	   fk(i) = fk(i)/f
	end do
c
	return
	end
c
c*********************************************************************
c
	function rintrz(x,y,xp,yp,zp)
c
c interpolates point in triangle and returns value
c
c x,y		vertices of triangle
c xp,yp		coordinates of point
c zp		values at vertices
c rintrz	interpolated value at (xp,yp)
c
	dimension x(3),y(3),zp(3)
c
	f=0.
	z=0.
	do i=1,3
	   ii=mod(i,3)+1
	   iii=mod(ii,3)+1
	   a=x(ii)*y(iii)-x(iii)*y(ii)
	   b=y(ii)-y(iii)
	   c=x(iii)-x(ii)
	   fk = a + xp*b + yp*c
	   z = z + zp(i)*fk
	   f=f+a
	end do
c
	rintrz = z/f
c
	return
	end
c
c********************************************************************
c
	function icutri(x,y,x1,y1,x2,y2,ds)
c
c intersects line with sides of triangle
c
c x,y		vertices of triangle
c x1,y1		coordinates of starting point of line
c x2,y2		coordinates of ending point of line
c ds		distance of intersection with vertex
c		...ds(i) = [0...1] intersection with line i
c		...ds(i) = -1 no intersection with line i
c		...(line i is opposit to vertex i)
c		...if ds(i) = [0...1] and xp,yp is intersection point
c		...then ds(i) = sqrt ( (xp-x(i+1))**2 + (yp-y(i+1))**2 )
c icutri	number of sides intersect with line
c
	dimension x(3),y(3),ds(3)
c
	icutri=0
c
	do i=1,3
	   ii=mod(i,3)+1
	   iii=mod(ii,3)+1
	   x3=x(ii)
	   y3=y(ii)
	   x4=x(iii)
	   y4=y(iii)
	   isw=icut(x1,y1,x2,y2,x3,y3,x4,y4,xp,yp)
	   if(isw.eq.0) then
		rpunt = (xp-x3)*(xp-x3) + (yp-y3)*(yp-y3)
		rside = (x4-x3)*(x4-x3) + (y4-y3)*(y4-y3)
		ds(i) = sqrt( rpunt/rside )
		icutri = icutri +1
	   else
		ds(i)=-1.
	   end if
	end do
c
	return
	end
