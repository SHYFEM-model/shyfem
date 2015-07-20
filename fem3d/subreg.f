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
c 26.01.2011	ggu&mb	handling extrapolation in am2av()
c 27.01.2011	ggu&ccf	bug fix in find_elem_from_old() BUG_27.01.2011
c 31.03.2011	ggu	new routine elemmask()
c 24.11.2011	ggu	new routine find_close_elem()
c 20.06.2012	ggu	new routine get_scal_elem()
c 07.10.2012	ggu	new routine av2fm()
c 10.10.2012	ggu	new routine fm2am2d() and fm2am3d()
c 26.10.2012	ggu	bug fix: do not access not existing storage
c 30.05.2014	ggu	in av2amk() do not interpolate for flag values
c 07.07.2014	ggu	new routine intp_reg()
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

	include 'reg.h'

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

	include 'reg.h'

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

	include 'reg.h'

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

	use basin

	implicit none

c arguments
	integer ip,jp
	real av(1)
	real am(ip,jp)
c parameter
	double precision eps
	parameter ( eps = 1.d-14 )
c common
	include 'param.h'
	include 'reg.h'
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

	use basin

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
	include 'param.h'
	include 'reg.h'
c local
	integer i,j,ii,iii,ie,k,kn,iin
	integer imin,imax,jmin,jmax
	integer iflag
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
	    iflag = 0
	    do i=1,3
		kn=nen3v(i,ie)
		x(i)=xgv(kn)
		y(i)=ygv(kn)
		z(i)=av(kn)
	        if( z(i) > pzlreg ) iflag = iflag + 1
	    end do
	    if( iflag .ne. 3 ) cycle

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
c************************************************
c************************************************

	subroutine av2fm(fm,ip,jp)

c computation of interpolation matrix of regular net (nodal values)

	implicit none

	real fm(4,ip,jp)	!values for interpolation (fm(4,i,j) = ie)
	integer ip,jp		!dimension of matrices

	logical bw		!use wet mask?
	logical bwater(1)	!wet mask for each element

	bw = .false.
	bwater(1) = .true.

	call av2fmk(bw,bwater,fm,ip,jp)

	end

c************************************************

	subroutine av2fmk(bw,bwater,fm,ip,jp)

c computation of interpolation matrix of regular net (nodal values) with mask
c
c the interpolation can be carried out as
c
c	do j=1,jp
c	  do i=1,ip
c	    ie = nint(fm(4,i,j))
c	    if( ie .gt. 0 ) then
c	      a = 0.
c	      do ii=1,3
c	        k = nen3v(ii,ie)
c		a = a + val(k) * fm(ii,i,j)
c	      end do
c	    else
c	      a = flag
c	    end if
c	    am(i,j) = a
c	  end do
c	end do
	        
	use basin

	implicit none

c arguments
	logical bw		!use wet mask?
	logical bwater(1)	!wet mask for each element
	real fm(4,ip,jp)	!values for interpolation (fm(4,i,j) = ie)
	integer ip,jp		!dimension of matrices
c parameter
	double precision eps
	parameter ( eps = 1.d-14 )
c common
	include 'param.h'
c pxareg,pyareg         coordinates of lower left point of matrix
c pxdreg,pydreg         grid size of matrix
c pzlreg                value of z for land points
	include 'reg.h'
c local
	integer i,j,ii,iii,ie,k,kn,iin
	integer imin,imax,jmin,jmax
	logical bok
	double precision x(3),y(3),z(3),a(3),b(3),c(3)
	double precision zh,fh,f,xp,yp
	double precision xmin,xmax,ymin,ymax
c function
	integer intrid

	do j=1,jp
	    do i=1,ip
		fm(4,i,j) = 0.
	    end do
	end do

	do ie=1,nel
	  bok = .true.
	  if( bw ) bok = bwater(ie)
	  if( bok ) then			!wet
	    do i=1,3
		kn=nen3v(i,ie)
		x(i)=xgv(kn)
		y(i)=ygv(kn)
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
			do ii=1,3
			   fh=(a(ii)+xp*b(ii)+yp*c(ii))/f
			   fm(ii,i,j) = fh
			end do
			fm(4,i,j) = ie
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
	stop 'error stop av2fm: area of element'
	end

c************************************************

        subroutine fm2am2d(femval,nx,ny,fm,am)

c interpolation 2d of fem values to regular grid using fm matrix

        implicit none

        real femval(1)			!values of fem array
        integer nx,ny			!dimension of regular matrix
        real fm(4,nx,ny)		!interpolation matrix
        real am(nx,ny)			!interpolated values (return)

	integer nlvdi,nlv,ilhv(1)

	nlvdi = 1
	nlv = 1
	ilhv(1) = 1

        call fm2am3d(nlvdi,ilhv,femval,nlv,nx,ny,fm,am)

	end

c************************************************

        subroutine fm2am3d(nlvdi,ilhv,femval,nlv,nx,ny,fm,am)

c interpolation 3d of fem values to regular grid using fm matrix

	use basin

        implicit none

	integer nlvdi			!vertical dimension of fem array
        integer ilhv(1)			!vertical discretization (element!!)
        real femval(nlvdi,1)		!values of fem array
        integer nlv,nx,ny		!dimension of regular matrix
        real fm(4,nx,ny)		!interpolation matrix
        real am(nlv,nx,ny)		!interpolated values (return)

	include 'param.h'

	include 'reg.h'

        integer i,j,l,lmax,ie,ii,k
        real a
        real flag

	flag = pzlreg

        do j=1,ny
          do i=1,nx
            ie = nint(fm(4,i,j))
            lmax = 0
            if( ie .gt. 0 ) then
              lmax = 1
              if( nlvdi .gt. 1 ) lmax = ilhv(ie)
            end if 
            do l=1,lmax
              a = 0.
              do ii=1,3
                k = nen3v(ii,ie)
                a = a + femval(l,k) * fm(ii,i,j)
              end do
              am(l,i,j) = a
            end do
            do l=lmax+1,nlv
              am(l,i,j) = flag
            end do
          end do
        end do

        end

c************************************************
c************************************************
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

	subroutine intp_reg(nx,ny,x0,y0,dx,dy,flag,regval,femval,ierr)

c interpolation of regular array onto fem grid
c
c ierr:
c		0	no errors
c		< 0	interpolation out of domain (extrapolation)
c		> 0	values of flag used in interpolation

	use basin

	implicit none

	integer nx,ny
	real x0,y0,dx,dy
	real flag
	real regval(nx,ny)
	real femval(*)		!interpolated values on fem grid (return)
	integer ierr		!error code (return)

	include 'param.h'

	logical bextra
	integer k
	integer imin,jmin
	integer iflag,iout
	real xx,yy,z1,z2,z3,z4,x1,y1,t,u
	real xn,yn
 
	iflag = 0	!used flag for interpolation
	iout = 0	!used outside point for interpolation

	xn = x0 + (nx-1)*dx
	yn = y0 + (ny-1)*dy

	do k=1,nkn
	    xx = xgv(k)
	    yy = ygv(k)
 
	    femval(k) = flag
 
	    if( xx .le. x0 ) then
	      imin = 1
	    else if( xx .ge. xn ) then
	      imin = nx-1
	    else
	      imin=1+(xx-x0)/dx
	    end if
	    if( yy .le. y0 ) then
	      jmin = 1
	    else if( yy .ge. yn ) then
	      jmin = ny-1
	    else
	      jmin=1+(yy-y0)/dy
	    end if

	    if( imin.lt.1 .or. jmin.lt.1 ) goto 99
	    if( imin+1.gt.nx .or. jmin+1.gt.ny ) goto 99

	    z1 = regval(imin,jmin)
	    z2 = regval(imin+1,jmin)
	    z3 = regval(imin+1,jmin+1)
	    z4 = regval(imin,jmin+1)

	    if( z1.eq.flag .or. z2.eq.flag ) iflag = iflag + 1
	    if( z3.eq.flag .or. z4.eq.flag ) iflag = iflag + 1

	    x1 = x0+(imin-1)*dx
	    y1 = y0+(jmin-1)*dy
	    t = (xx-x1)/dx
	    u = (yy-y1)/dy

	    if( u.gt.1. .or. u.lt.0. ) iout = iout + 1
	    if( t.gt.1. .or. t.lt.0. ) iout = iout + 1

	    femval(k)=(1-t)*(1-u)*z1+t*(1-u)*z2+t*u*z3+(1-t)*u*z4
	end do
 
	ierr = 0
	if( iout .gt. 0 ) ierr = -iout
	if( iflag .gt. 0 ) ierr = iflag

	return
   99	continue
	write(6,*) imin,jmin,nx,ny
	stop 'error stop intp_reg: internal error (1)'
	end

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
	use basin

	implicit none
c
c arguments
	integer ip,jp
	real av(1)
	real am(ip,jp)
c common
	include 'param.h'
	include 'reg.h'
c local
	logical bextra
	integer k
	integer imin,jmin
	real xx,yy,z1,z2,z3,z4,x1,y1,t,u
 
	bextra = .false.   !extrapolate if out of regular domain (else error)

	do k=1,nkn
	    xx=xgv(k)
	    yy=ygv(k)
 
	    av(k)=pzlreg
 
	    imin=(xx-pxareg)/pxdreg
	    jmin=(yy-pyareg)/pydreg
	    imin=imin+1
	    jmin=jmin+1

	    !if( imin.lt.1 .or. jmin.lt.1 ) goto 1
	    !if( imin+1.gt.ip .or. jmin+1.gt.jp ) goto 1

	    if( bextra ) then
	      if( imin.lt.1 ) imin = 1
	      if( jmin.lt.1 ) jmin = 1
	      if( imin+1.gt.ip ) imin = ip - 1
	      if( jmin+1.gt.jp ) jmin = jp - 1
	    else
	      if( imin.lt.1 .or. jmin.lt.1 ) goto 1
	      if( imin+1.gt.ip .or. jmin+1.gt.jp ) goto 1
	    end if

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
c parameter
	real eps,zero,one
	parameter( eps = 1.e-4 )
	!parameter( eps = 0. )		!to use this pass in double precision
	parameter( zero = 0. - eps , one = 1. + eps )
c common
	include 'reg.h'
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

	    if( u.gt.one .or. u.lt.zero .or.
     +			 t.gt.one .or. t.lt.zero ) then
		write(6,*) 'error am2val',u,t,xx,yy
		write(6,*) '  ',imin,jmin,ip,jp
	    end if

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
	use basin

	implicit none
c
c arguments
	integer ip,jp
	real av(1)
	real am(ip,jp)
c common
	include 'param.h'
	include 'reg.h'
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
c******************************************************
c******************************************************

	subroutine set_dry_mask(bwater,zv,href,hzoff)

c makes mask for dry and wet areas - zenv must be available
c
c bwater is elementwise mask:	true = water point

	use mod_hydro
	use basin

	implicit none

c arguments
	logical bwater(1)
	real zv(1)
	real href,hzoff
c common
	include 'param.h'
c local
	integer itot,itot1
	integer ie,ii

        do ie=1,nel

          itot=0
          do ii=1,3
            if( hm3v(ii,ie)+zenv(ii,ie)-href .gt. hzoff ) then
		itot=itot+1    !wet
	    end if
          end do

          itot1=0
          do ii=1,3
            if(zv(nen3v(ii,ie)).eq.zenv(ii,ie)) itot1=itot1+1
          end do

          if(itot.ne.3.or.itot1.ne.3)  bwater(ie) = .false.

        end do

	end

c******************************************************

	subroutine init_dry_mask(bwater)

c initializes mask for water points
c
c bwater is elementwise mask:	true = water point

	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	logical bwater(1)
c common
c local
	integer ie

	do ie=1,nel
	  bwater(ie) = .true.
	end do

	end

c******************************************************

	subroutine set_level_mask(bwater,ilhv,level)

c makes mask for water points (level)
c
c bwater is elementwise mask:	true = water point

	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	logical bwater(1)
	integer ilhv(1)
	integer level
c common
c local
	integer ie,nedry

	nedry = 0

	do ie=1,nel
	  if( level .gt. ilhv(ie) ) then
	    bwater(ie) = .false.
	    nedry = nedry + 1
	  end if
	end do

	write(6,*) 'levelmask: no such level =',nedry,nel,level

	end

c******************************************************

	subroutine make_dry_node_mask(bwater,bkwater)

c makes node mask from element mask
c
c bwater is elementwise mask:	true = water point

	use basin

	implicit none

c arguments
	logical bwater(1)
	logical bkwater(1)
c common
	include 'param.h'
c local
	integer ie,ii,k
	integer nndry,nedry

	nndry = 0
	nedry = 0

	do k=1,nkn
	  bkwater(k) = .false.
	end do

	do ie=1,nel
	  if( bwater(ie) ) then
	    do ii=1,3
	      k = nen3v(ii,ie)
	      bkwater(k) = .true.
	    end do
	  else
	    nedry = nedry + 1
	  end if
	end do

	do k=1,nkn
	  if( .not. bkwater(k) ) nndry = nndry + 1
	end do

	write(6,*) 'make_dry_node_mask: dry elements =',nedry,nel
	write(6,*) 'make_dry_node_mask: dry nodes    =',nndry,nkn

	end

c******************************************************

	subroutine make_dry_elem_mask(bwater,bkwater)

c makes elem mask from node mask
c
c bwater is elementwise mask:	true = water point

	use basin

	implicit none

c arguments
	logical bwater(1)
	logical bkwater(1)
c common
	include 'param.h'
c local
	integer ie,ii,k
	integer nedry

	nedry = 0

	do ie=1,nel
	  bwater(ie) = .true.
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( .not. bkwater(k) ) then
	      bwater(ie) = .false.
	    end if
	  end do
	  if( .not. bwater(ie) ) nedry = nedry + 1
	end do

	write(6,*) 'make_dry_elem_mask: dry elements =',nedry,nel

	end

c******************************************************
c******************************************************
c******************************************************

	subroutine mimareg(am,ip,jp,amin,amax)

c computes min/max of regular matrix (without flag values)

	implicit none

	integer ip,jp
	real am(ip,jp)
	real amin,amax

	include 'reg.h'

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

	include 'reg.h'

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

	use evgeom

        integer ie      !element
        real z(3)       !values on nodes
        real xp,yp      !coordinates of point
        real zp         !interpolated value (return)

        integer ii
        double precision zh,a,b,c,w


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

	subroutine find_close_elem(ieold,xp,yp,ielem)

c finds element for point (xp,yp) starting from ieold
c
c uses data structure ev and ieltv

	use mod_geom
	use evgeom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	integer ieold
	real xp,yp
	integer ielem	!element number on return



	logical binit,bdebug
	integer ie,ii,iside,lmax,loop
	real xi,ximin
	double precision a(3),b(3),c(3)

	logical in_element

	bdebug = .false.
	lmax = 10

	call is_init_ev(binit)

c-------------------------------------------------------------
c check if old element is given -> if not test all elements
c-------------------------------------------------------------

	if( bdebug ) write(6,*) 'ggu_xi (1) ',ieold
	if( ieold .le. 0 .or. ieold .gt. nel ) then
	  call find_element(xp,yp,ielem)
	  return
	end if

	if( bdebug ) write(6,*) 'ggu_xi (2) ',ieold
	if( .not. binit ) then
	  call find_elem_from_old(ieold,xp,yp,ielem)
	  return
	end if

c-------------------------------------------------------------
c start from old element
c-------------------------------------------------------------

	if( bdebug ) write(6,*) 'ggu_xi (3) ',ieold

	loop = 0
	ie = ieold
	do while( ie .gt. 0 )
	  iside = 0
	  ximin = 1.e+30
	  call xi_abc(ie,a,b,c)
	  do ii=1,3
	    xi = a(ii) + b(ii)*xp + c(ii)*yp
	    if( bdebug ) write(6,*) 'ggu_xiii ',ie,ii,xi
	    if( xi .lt. ximin ) then
	      ximin = xi
	      iside = ii
	    end if
	  end do
	  if( ximin .ge. 0. ) then
	    ielem = ie
	    return
	  end if
	  if( iside .le. 0 .or. iside .gt. 3 ) then
	    if( bdebug ) write(6,*) '******** ',iside,ie,ximin,xi
	    ie = 0
	  else
	    ie = ieltv(iside,ie)
	    if( bdebug ) write(6,*) 'ggu_xiii iterate',ie,iside,ximin
	  end if
	  loop = loop + 1
	  if( loop .gt. lmax ) ie = 0
	end do

	ielem = 0

	end

c******************************************************

	subroutine find_elem_from_old(ieold,xp,yp,ielem)

c finds element for point (xp,yp) starting from ieold

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ieold
	real xp,yp
	integer ielem	!element number on return


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
	if( iem .lt. 1 ) iem = nel		!BUG_27.01.2011
	iep = ieold+1
	if( iep .gt. nel ) iep = 1		!BUG_27.01.2011

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

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real xp,yp
	integer ielem	!element number on return


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

	use basin

	implicit none

	logical in_element
	integer ie
	real xp,yp

	include 'param.h'

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

	use basin

	implicit none

	integer ie
	real x(3), y(3)

	include 'param.h'

	integer ii,k

	do ii=1,3
	  k = nen3v(ii,ie)
	  x(ii) = xgv(k)
	  y(ii) = ygv(k)
	end do

	end

c******************************************************

	subroutine get_scal_elem(ie,sv,s)

c returns s at vertices of element ie

	use basin

	implicit none

	integer ie
	real sv(1)
	real s(3)

	include 'param.h'

	integer ii,k

	do ii=1,3
	  k = nen3v(ii,ie)
	  s(ii) = sv(k)
	end do

	end

c******************************************************

