c
c $Id: subreg.f,v 1.17 2009-09-14 08:20:58 georg Exp $
c
c routines for interpolation onto regular grid
c
c contents :
c
c subroutine setgeo(x0,y0,dx,dy,flag)
c		sets grid values for interpolation
c subroutine getgeo(x0,y0,dx,dy,flag)
c		gets grid values for interpolation
c subroutine getgeoflag(flag)
c		gets flag value for interpolation
c
c subroutine av2am(av,am,ip,jp)
c		interpolation of av onto a regular net
c subroutine av2amk(bwater,av,am,ip,jp)
c		interpolation of av onto a regular net
c function intri(x,y,xp,yp)
c		point in triangle or not
c function intrid(x,y,xp,yp)
c		point in triangle or not (double precision version)
c subroutine am2av(am,av,ip,jp)
c		interpolation of am onto finite element mesh
c function am2val(am,ip,jp,xx,yy)
c		interpolation of am onto finite element mesh
c subroutine ave2am(av,am,ip,jp)
c		interpolation of av (elementwise) onto a regular net
c
c subroutine mkmask(bwater,zv,href,hzoff)
c		makes mask in element
c
c subroutine mimareg(am,ip,jp,amin,amax)
c		computes min/max of regular matrix (without flag values)
c subroutine a2char(am,ac,ip,jp)
c		creates 1 char representation of matrix
c subroutine prchar(ac,ip,jp)
c		prints 1 char representation of matrix
c
c subroutine femintp(ie,z,xp,yp,zp)
c               interpolation in element (with ev)
c subroutine elemintp(x,y,z,xp,yp,zp)
c               interpolation in element (no ev)
c
c subroutine find_elem_from_old(ieold,xp,yp,ielem)
c		finds element for point (xp,yp) starting from ieold
c subroutine find_element(xp,yp,ielem)
c		finds element for point (xp,yp)
c function in_element(ie,xp,yp)
c		checks if point (xp,yp) is in element ie
c subroutine get_xy_elem(ie,x,y)
c		returns x,y of vertices of element ie
c
c revision log :
c
c 18.11.1998	ggu	routine commented
c 18.11.1998	ggu	routine setgeo introduced
c 19.11.1998	ggu	routines a2char, prchar added
c 19.10.1999	ggu	routine mkmask added from subutl
c 25.11.2004	ggu	new routines femintp and elemintp for interpolation
c 14.03.2005	ggu	new routines for interpolation in element
c 11.03.2009	ggu	new helper routine getgeoflag()
c 12.06.2009	ggu	passing to double precision, intrid, bug bug_f_64bit
c
c notes :
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points
c
c******************************************************

	subroutine setgeo(x0,y0,dx,dy,flag)

c sets grid values for interpolation
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points

	implicit none

	real x0,y0	!coordinates of lower,left point
	real dx,dy	!grid spacing in x/y direction
	real flag	!flag for land points

        real pxareg,pyareg,pxdreg,pydreg,pzlreg
        common /ppp20/ pxareg,pyareg,pxdreg,pydreg,pzlreg
	save /ppp20/

	pxareg = x0
	pyareg = y0
	pxdreg = dx
	pydreg = dy
	pzlreg = flag

	end

c******************************************************

	subroutine getgeo(x0,y0,dx,dy,flag)

c gets grid values for interpolation
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points

	implicit none

	real x0,y0	!coordinates of lower,left point
	real dx,dy	!grid spacing in x/y direction
	real flag	!flag for land points

        real pxareg,pyareg,pxdreg,pydreg,pzlreg
        common /ppp20/ pxareg,pyareg,pxdreg,pydreg,pzlreg
	save /ppp20/

	x0   = pxareg
	y0   = pyareg
	dx   = pxdreg
	dy   = pydreg
	flag = pzlreg

	end

c******************************************************

	subroutine getgeoflag(flag)

c gets flag value for interpolation
c
c flag                value for land points

	implicit none

	real flag	!flag for land points

        real pxareg,pyareg,pxdreg,pydreg,pzlreg
        common /ppp20/ pxareg,pyareg,pxdreg,pydreg,pzlreg
	save /ppp20/

	flag = pzlreg

	end

c******************************************************

	subroutine av2am(av,am,ip,jp)

c interpolation of av onto a regular net (nodal values)
c
c av                    array to be interpolated
c am                    matrices of interpolated values (u,v,z,h)
c ip,jp                 dimension of matrices
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points

	implicit none

c arguments
	integer ip,jp
	real av(1)
	real am(ip,jp)
c parameter
	double precision eps
	parameter ( eps = 1.d-14 )
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	real xgv(1), ygv(1)
	real pxareg,pyareg,pxdreg,pydreg,pzlreg
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /xgv/xgv, /ygv/ygv
	common /nen3v/nen3v
	common /ppp20/ pxareg,pyareg,pxdreg,pydreg,pzlreg
c local
	integer i,j,ii,iii,ie,k,kn,iin
	integer imin,imax,jmin,jmax
	double precision x(3),y(3),z(3),a(3),b(3),c(3)
	double precision zh,fh,f,xp,yp
	double precision xmin,xmax,ymin,ymax
c function
	integer intrid

	do j=1,jp
	    do i=1,ip
		am(i,j)=pzlreg
	    end do
	end do

	do ie=1,nel
	    do i=1,3
		kn=nen3v(i,ie)
		x(i)=xgv(kn)
		y(i)=ygv(kn)
		z(i)=av(kn)
	    end do

	    !f=0.
	    do i=1,3
		ii=mod(i,3)+1
		iii=mod(ii,3)+1
		a(i)=x(ii)*y(iii)-x(iii)*y(ii)
		b(i)=y(ii)-y(iii)
		c(i)=x(iii)-x(ii)
		!f=f+a(i)
	    end do
	    f = c(3)*b(2) - c(2)*b(3)		!bug_f_64bit
	    if( f .le. eps ) goto 99

	    xmax=max(x(1),x(2),x(3))
	    xmin=min(x(1),x(2),x(3))
	    ymin=min(y(1),y(2),y(3))
	    ymax=max(y(1),y(2),y(3))

	    imin=(xmin-pxareg)/pxdreg+1.99
	    imax=(xmax-pxareg)/pxdreg+1.01
	    jmin=(ymin-pyareg)/pydreg+1.99
	    jmax=(ymax-pyareg)/pydreg+1.01

	    if(imin.lt.1) imin=1
	    if(imax.gt.ip)imax=ip
	    if(jmin.lt.1) jmin=1
	    if(jmax.gt.jp)jmax=jp

	    do i=imin,imax
		do j=jmin,jmax
		    xp=(i-1)*pxdreg+pxareg
		    yp=(j-1)*pydreg+pyareg

		    iin=intrid(x,y,xp,yp)

		    if(iin.ne.0) then
			zh=0.
			do k=1,3
			   fh=(a(k)+xp*b(k)+yp*c(k))/f
			   zh=zh+z(k)*fh
			end do
			am(i,j)=zh
		    end if
		end do
	    end do
	end do

	return
   99	continue
	write(6,*) ie,f
	write(6,*) x
	write(6,*) y
	write(6,*) a
	write(6,*) b
	write(6,*) c
	stop 'error stop av2am: area of element'
	end

c******************************************************

	subroutine av2amk(bwater,av,am,ip,jp)

c interpolation of av onto a regular net (nodal values) with mask
c
c bwater		mask for water points
c av                    array to be interpolated
c am                    matrices of interpolated values (u,v,z,h)
c ip,jp                 dimension of matrices
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points

	implicit none

c arguments
	integer ip,jp
	real av(1)
	real am(ip,jp)
	logical bwater(1)
c parameter
	double precision eps
	parameter ( eps = 1.d-14 )
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	real xgv(1), ygv(1)
	real pxareg,pyareg,pxdreg,pydreg,pzlreg
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /xgv/xgv, /ygv/ygv
	common /nen3v/nen3v
	common /ppp20/ pxareg,pyareg,pxdreg,pydreg,pzlreg
c local
	integer i,j,ii,iii,ie,k,kn,iin
	integer imin,imax,jmin,jmax
	double precision x(3),y(3),z(3),a(3),b(3),c(3)
	double precision zh,fh,f,xp,yp
	double precision xmin,xmax,ymin,ymax
c function
	integer intrid

	do j=1,jp
	    do i=1,ip
		am(i,j)=pzlreg
	    end do
	end do

	do ie=1,nel
	  if( bwater(ie) ) then	!wet
	    do i=1,3
		kn=nen3v(i,ie)
		x(i)=xgv(kn)
		y(i)=ygv(kn)
		z(i)=av(kn)
	    end do

	    !f=0.
	    do i=1,3
		ii=mod(i,3)+1
		iii=mod(ii,3)+1
		a(i)=x(ii)*y(iii)-x(iii)*y(ii)
		b(i)=y(ii)-y(iii)
		c(i)=x(iii)-x(ii)
		!f=f+a(i)
	    end do
	    f = c(3)*b(2) - c(2)*b(3)		!bug_f_64bit
	    if( f .le. eps ) goto 99

	    xmax=max(x(1),x(2),x(3))
	    xmin=min(x(1),x(2),x(3))
	    ymin=min(y(1),y(2),y(3))
	    ymax=max(y(1),y(2),y(3))

	    imin=(xmin-pxareg)/pxdreg+1.99
	    imax=(xmax-pxareg)/pxdreg+1.01
	    jmin=(ymin-pyareg)/pydreg+1.99
	    jmax=(ymax-pyareg)/pydreg+1.01

	    if(imin.lt.1) imin=1
	    if(imax.gt.ip)imax=ip
	    if(jmin.lt.1) jmin=1
	    if(jmax.gt.jp)jmax=jp

	    do i=imin,imax
		do j=jmin,jmax
		    xp=(i-1)*pxdreg+pxareg
		    yp=(j-1)*pydreg+pyareg

		    iin=intrid(x,y,xp,yp)

		    if(iin.ne.0) then
			zh=0.
			do k=1,3
			   fh=(a(k)+xp*b(k)+yp*c(k))/f
			   zh=zh+z(k)*fh
			end do
			am(i,j)=zh
		    end if
		end do
	    end do
	  end if
	end do

	return
   99	continue
	write(6,*) ie,f
	write(6,*) x
	write(6,*) y
	write(6,*) a
	write(6,*) b
	write(6,*) c
	stop 'error stop av2amk: area of element'
	end

c************************************************

	function intri(x,y,xp,yp)

c point in triangle or not
c
c x,y		array of coordinates of vertices of triangle
c xp,yp		coordinates of point
c intri		1: point is in triangle  0: point outside (return value)

	implicit none

c arguments
	integer intri
	real x(3),y(3),xp,yp
c local
	integer k1,k2
	double precision x21,y21,xn,yn
	double precision scal,eps
c save
	save eps
	data eps /1.e-13/

	intri=0

	do k1=1,3
	   k2=mod(k1,3)+1
	   x21=x(k2)-x(k1)
	   y21=y(k2)-y(k1)
	   yn = x21
	   xn = -y21
	   scal=(xp-x(k1))*xn+(yp-y(k1))*yn
	   if(scal.lt.0.) return
	end do

	intri=1	!inside

	end

c************************************************

	function intrid(x,y,xp,yp)

c point in triangle or not (double precision version)
c
c x,y		array of coordinates of vertices of triangle
c xp,yp		coordinates of point
c intri		1: point is in triangle  0: point outside (return value)

	implicit none

c arguments
	integer intrid
	double precision x(3),y(3),xp,yp
c local
	integer k1,k2
	double precision x21,y21,xn,yn
	double precision scal,eps
c save
	save eps
	data eps /1.e-13/

	intrid=0

	do k1=1,3
	   k2=mod(k1,3)+1
	   x21=x(k2)-x(k1)
	   y21=y(k2)-y(k1)
	   yn = x21
	   xn = -y21
	   scal=(xp-x(k1))*xn+(yp-y(k1))*yn
	   if(scal.lt.0.) return
	end do

	intrid=1	!inside

	end

c****************************************************************
c
	function intri0(x,y,xp,yp)
c
c point in triangle or not
c
c x,y		array of coordinates of vertices of triangle
c xp,yp		coordinates of point
c intri		1: point is in triangle  0: point outside (return value)
c
	implicit none
c
c arguments
	integer intri0
	real x(3),y(3),xp,yp
c local
	integer i,k1,k2
	double precision xs,ys,x12,y12
	double precision det,detlam,rlamb
	double precision eps
c save
	save eps
	data eps /1.e-13/
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
	intri0=0
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
	intri0=1	!inside
c
	return
	end
c
c****************************************************************
c
	subroutine am2av(am,av,ip,jp)
c
c interpolation of am onto finite element mesh
c
c am                    matrix of value
c av                    finite element array on return
c ip,jp                 dimension of matrix
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points
c
	implicit none
c
c arguments
	integer ip,jp
	real av(1)
	real am(ip,jp)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real xgv(1), ygv(1)
	real pxareg,pyareg,pxdreg,pydreg,pzlreg
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /xgv/xgv, /ygv/ygv
	common /ppp20/ pxareg,pyareg,pxdreg,pydreg,pzlreg
c local
	integer k
	integer imin,jmin
	real xx,yy,z1,z2,z3,z4,x1,y1,t,u
 
	do k=1,nkn
	    xx=xgv(k)
	    yy=ygv(k)
 
	    av(k)=pzlreg
 
	    imin=(xx-pxareg)/pxdreg
	    jmin=(yy-pyareg)/pydreg
	    imin=imin+1
	    jmin=jmin+1

	    if( imin.lt.1 .or. jmin.lt.1 ) goto 1
	    if( imin+1.gt.ip .or. jmin+1.gt.jp ) goto 1

	    z1=am(imin,jmin)
	    z2=am(imin+1,jmin)
	    z3=am(imin+1,jmin+1)
	    z4=am(imin,jmin+1)

	    if( z1.eq.pzlreg .or. z2.eq.pzlreg ) goto 1
	    if( z3.eq.pzlreg .or. z4.eq.pzlreg ) goto 1

	    x1=pxareg+(imin-1)*pxdreg
	    y1=pyareg+(jmin-1)*pydreg
	    t=(xx-x1)/pxdreg
	    u=(yy-y1)/pydreg

	if( u.gt.1. .or. u.lt.0. ) write(6,*) 'error am2av',u,t,k
	if( t.gt.1. .or. t.lt.0. ) write(6,*) 'error am2av',u,t,k

	    av(k)=(1-t)*(1-u)*z1+t*(1-u)*z2+t*u*z3+(1-t)*u*z4

    1	    continue
	end do
 
	return
	end
c
c****************************************************************
c
	function am2val(am,ip,jp,xx,yy)
c
c interpolation of am onto finite element mesh
c
c am                    matrix of value
c ip,jp                 dimension of matrix
c xx,yy			coordinates for desired value
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points
c
	implicit none
c
c arguments
	real am2val
	integer ip,jp
	real am(ip,jp)
	real xx,yy
c common
	real pxareg,pyareg,pxdreg,pydreg,pzlreg
	common /ppp20/ pxareg,pyareg,pxdreg,pydreg,pzlreg
c local
	integer k
	integer imin,jmin
	real z1,z2,z3,z4,x1,y1,t,u
 
	    am2val=pzlreg
 
	    imin=(xx-pxareg)/pxdreg
	    jmin=(yy-pyareg)/pydreg
	    imin=imin+1
	    jmin=jmin+1

	    if( imin.lt.1 .or. jmin.lt.1 ) goto 1
	    if( imin+1.gt.ip .or. jmin+1.gt.jp ) goto 1

	    z1=am(imin,jmin)
	    z2=am(imin+1,jmin)
	    z3=am(imin+1,jmin+1)
	    z4=am(imin,jmin+1)

	    if( z1.eq.pzlreg .or. z2.eq.pzlreg ) goto 1
	    if( z3.eq.pzlreg .or. z4.eq.pzlreg ) goto 1

	    x1=pxareg+(imin-1)*pxdreg
	    y1=pyareg+(jmin-1)*pydreg
	    t=(xx-x1)/pxdreg
	    u=(yy-y1)/pydreg

	if( u.gt.1. .or. u.lt.0. ) write(6,*) 'error am2val',u,t,k
	if( t.gt.1. .or. t.lt.0. ) write(6,*) 'error am2val',u,t,k

	    am2val=(1-t)*(1-u)*z1+t*(1-u)*z2+t*u*z3+(1-t)*u*z4

    1	    continue
 
	return
	end
c
c******************************************************
c
	subroutine ave2am(av,am,ip,jp)
c
c interpolation of av onto a regular net (element values)
c
c av                    array to be interpolated
c am                    matrices of interpolated values 
c ip,jp                 dimension of matrices
c
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points
c
	implicit none
c
c arguments
	integer ip,jp
	real av(1)
	real am(ip,jp)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	real xgv(1), ygv(1)
	real pxareg,pyareg,pxdreg,pydreg,pzlreg
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /xgv/xgv, /ygv/ygv
	common /nen3v/nen3v
	common /ppp20/ pxareg,pyareg,pxdreg,pydreg,pzlreg
c local
c	integer i,j,ii,iii,ie,k,kn,iin
	integer i,j,ie,kn,iin
	integer imin,imax,jmin,jmax
	real x(3),y(3)
	real zh,xp,yp
	real xmin,xmax,ymin,ymax
c function
	integer intri
c
	do j=1,jp
	    do i=1,ip
		am(i,j)=pzlreg
	    end do
	end do
c
	do ie=1,nel
	    zh=av(ie)
	    do i=1,3
		kn=nen3v(i,ie)
		x(i)=xgv(kn)
		y(i)=ygv(kn)
	    end do
c
	    xmax=max(x(1),x(2),x(3))
	    xmin=min(x(1),x(2),x(3))
	    ymin=min(y(1),y(2),y(3))
	    ymax=max(y(1),y(2),y(3))
c
	    imin=(xmin-pxareg)/pxdreg+1.99
	    imax=(xmax-pxareg)/pxdreg+1.01
	    jmin=(ymin-pyareg)/pydreg+1.99
	    jmax=(ymax-pyareg)/pydreg+1.01
c
	    if(imin.lt.1) imin=1
	    if(imax.gt.ip)imax=ip
	    if(jmin.lt.1) jmin=1
	    if(jmax.gt.jp)jmax=jp
c
	    do i=imin,imax
		do j=jmin,jmax
		    xp=(i-1)*pxdreg+pxareg
		    yp=(j-1)*pydreg+pyareg
c
		    iin=intri(x,y,xp,yp)
c
		    if(iin.ne.0) then
			am(i,j)=zh
		    end if
		end do
	    end do
	end do
c
	return
	end

c******************************************************

	subroutine drymask(bwater,zv,href,hzoff)

c makes mask for dry and wet areas
c
c bwater is elementwise mask:	true = water point

	implicit none

c arguments
	logical bwater(1)
	real zv(1)
	real href,hzoff
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	real hm3v(3,1)
	real zenv(3,1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nen3v/nen3v, /hm3v/hm3v
	common /zenv/zenv
c local
	integer idry,idry1,idryd
	integer itot,itot1
	integer ie,ii

        do ie=1,nel

          itot=0
          do ii=1,3
            if( hm3v(ii,ie)+zenv(ii,ie)-href .gt. hzoff ) itot=itot+1    !wet
          end do

          itot1=0
          do ii=1,3
            if(zv(nen3v(ii,ie)).eq.zenv(ii,ie)) itot1=itot1+1
          end do

          if(itot.ne.3.or.itot1.ne.3)  bwater(ie) = .false.

        end do

	end

c******************************************************

	subroutine initmask(bwater)

c initializes mask for water points
c
c bwater is elementwise mask:	true = water point

	implicit none

c arguments
	logical bwater(1)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
c local
	integer ie

	do ie=1,nel
	  bwater(ie) = .true.
	end do

	end

c******************************************************

	subroutine levelmask(bwater,ilhv,level)

c makes mask for water points (level)
c
c bwater is elementwise mask:	true = water point

	implicit none

c arguments
	logical bwater(1)
	integer ilhv(1)
	integer level
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
c local
	integer ie

	do ie=1,nel
	  if( level .gt. ilhv(ie) ) bwater(ie) = .false.
	end do

	end

c******************************************************

	subroutine nodemask(bwater,bkwater)

c makes node mask from element mask
c
c bwater is elementwise mask:	true = water point

	implicit none

c arguments
	logical bwater(1)
	logical bkwater(1)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
c local
	integer ie,ii,k

	do k=1,nkn
	  bkwater(k) = .false.
	end do

	do ie=1,nel
	  if( bwater(ie) ) then
	    do ii=1,3
	      k = nen3v(ii,ie)
	      bkwater(k) = .true.
	    end do
	  end if
	end do

	end

c******************************************************

	subroutine mimareg(am,ip,jp,amin,amax)

c computes min/max of regular matrix (without flag values)

	implicit none

	integer ip,jp
	real am(ip,jp)
	real amin,amax

        real pxareg,pyareg,pxdreg,pydreg,pzlreg
        common /ppp20/ pxareg,pyareg,pxdreg,pydreg,pzlreg

	integer i,j
	real a,high

	high = 1.e+30

	amin =  high
	amax = -high

	do j=1,jp
	  do i=1,ip
	    a = am(i,j)
	    if( a .ne. pzlreg ) then
	      if( a .gt. amax ) amax = a
	      if( a .lt. amin ) amin = a
	    end if
	  end do
	end do

	if( amin .eq.  high ) amin = pzlreg
	if( amax .eq. -high ) amax = pzlreg

	end

c******************************************************

	subroutine a2char(am,ac,ip,jp)

c creates 1 char representation of matrix

	implicit none

	integer ip,jp		!dimension of matrix
	real am(ip,jp)		!matrix containing data
	character*1 ac(ip,jp)	!matrix containing chars on return

        real pxareg,pyareg,pxdreg,pydreg,pzlreg
        common /ppp20/ pxareg,pyareg,pxdreg,pydreg,pzlreg

	integer i,j

	do j=1,jp
	  do i=1,ip
	    if( am(i,j) .eq. pzlreg ) then	!land
		ac(i,j) = '.'
	    else
		ac(i,j) = '*'			!data
	    end if
	  end do
	end do

	end

c******************************************************

	subroutine prchar(ac,ip,jp)

c prints 1 char representation of matrix

	implicit none

	integer ip,jp		!dimension of matrix
	character*1 ac(ip,jp)	!matrix containing chars

	character*256 line
	integer ipmax
	integer i,j

	ipmax = min(256,ip)

	do j=jp,1,-1
	  line = ' '
	  do i=1,ipmax
	    line(i:i) = ac(i,j)
	  end do
	  write(6,'(a)') line(1:ipmax)
	end do

	end

c******************************************************

        subroutine femintp(ie,z,xp,yp,zp)

c interpolation in element (with ev)
c
c interpolates in element ie from nodal values z to point xp,yp
c result is in zp
c
c needs array ev

        integer ie      !element
        real z(3)       !values on nodes
        real xp,yp      !coordinates of point
        real zp         !interpolated value (return)

        integer ii
        double precision zh,a,b,c,w

	include 'ev.h'

        zh=0.
        do ii=1,3
          a = ev(ii,ie)
          b = ev(3+ii,ie)
          c = ev(6+ii,ie)
          w = a + b*xp + c*yp
          zh = zh + z(ii) * w
        end do

        zp = zh

        end

c******************************************************

        subroutine elemintp(x,y,z,xp,yp,zp)

c interpolation in element (no ev)
c
c interpolates in element given by x,y nodal values z to point xp,yp
c result is in zp
c
c needs no other vectors but needs x,y of nodes

        real x(3),y(3)  !coordinates of nodes
        real z(3)       !values on nodes
        real xp,yp      !coordinates of point
        real zp         !interpolated value (return)

	double precision eps
	parameter ( eps = 1.d-14 )

        integer i,ii,iii
        double precision zh,f,fh
        double precision a(3),b(3),c(3)

        f = 0.
        do i=1,3
           ii=mod(i,3)+1
           iii=mod(ii,3)+1
           a(i)=x(ii)*y(iii)-x(iii)*y(ii)
           b(i)=y(ii)-y(iii)
           c(i)=x(iii)-x(ii)
           f=f+a(i)
        end do
	f = c(3)*b(2) - c(2)*b(3)		!bug_f_64bit
	if( f .le. eps ) goto 99

        zh=0.
        do i=1,3
           fh=(a(i)+xp*b(i)+yp*c(i))/f
           zh=zh+z(i)*fh
        end do

        zp = zh

	return
   99	continue
	write(6,*) 0,f
	write(6,*) x
	write(6,*) y
	write(6,*) a
	write(6,*) b
	write(6,*) c
	stop 'error stop elemintp: area of element'
        end

c******************************************************

	subroutine find_elem_from_old(ieold,xp,yp,ielem)

c finds element for point (xp,yp) starting from ieold

	implicit none

	integer ieold
	real xp,yp
	integer ielem	!element number on return

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	logical in_element
	integer iem,iep

c-------------------------------------------------------------
c check if old element is given -> if not test all elements
c-------------------------------------------------------------

	if( ieold .le. 0 .or. ieold .gt. nel ) then
	  call find_element(xp,yp,ielem)
	  return
	end if

c-------------------------------------------------------------
c check if in old element
c-------------------------------------------------------------

	if( in_element(ieold,xp,yp) ) then
	  ielem = ieold
	  return
	end if

c-------------------------------------------------------------
c start from old element going upwards and downwards
c-------------------------------------------------------------

	iem = ieold-1
	iep = ieold+1

	do while( iem .ne. ieold .and. iep .ne. ieold )
	  if( in_element(iem,xp,yp) ) then
	    ielem = iem
	    return
	  end if
	  iem = iem - 1
	  if( iem .lt. 1 ) iem = nel

	  if( in_element(iep,xp,yp) ) then
	    ielem = iep
	    return
	  end if
	  iep = iep + 1
	  if( iep .gt. nel ) iep = 1
	end do

c-------------------------------------------------------------
c no element found
c-------------------------------------------------------------

	ielem = 0

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c******************************************************

	subroutine find_element(xp,yp,ielem)

c finds element for point (xp,yp)

	implicit none

	real xp,yp
	integer ielem	!element number on return

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ie
	logical in_element

	do ie=1,nel
	  if( in_element(ie,xp,yp) ) then
		  ielem = ie
		  return
	  end if
	end do

	ielem = 0

	end

c******************************************************

	function in_element(ie,xp,yp)

c checks if point (xp,yp) is in element ie

	implicit none

	logical in_element
	integer ie
	real xp,yp

	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	integer nen3v(3,1)
	common /nen3v/nen3v

	integer ii,k,in
	real xmin,ymin,xmax,ymax
	real x(3),y(3)

	integer intri

	in_element = .false.

	do ii=1,3
	    k = nen3v(ii,ie)
	    x(ii) = xgv(k)
	    y(ii) = ygv(k)
	end do

	xmin = min(x(1),x(2),x(3))
	ymin = min(y(1),y(2),y(3))
	xmax = max(x(1),x(2),x(3))
	ymax = max(y(1),y(2),y(3))

	if( xp .ge. xmin .and. xp .le. xmax ) then
	  if( yp .ge. ymin .and. yp .le. ymax ) then
		in = intri(x,y,xp,yp)
		if( in .gt. 0 ) in_element = .true.
	  end if
	end if

	end

c******************************************************

	subroutine get_xy_elem(ie,x,y)

c returns x,y of vertices of element ie

	implicit none

	integer ie
	real x(3), y(3)

	integer nen3v(3,1)
	common /nen3v/nen3v
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv

	integer ii,k

	do ii=1,3
	  k = nen3v(ii,ie)
	  x(ii) = xgv(k)
	  y(ii) = ygv(k)
	end do

	end

c******************************************************

