c
c $Id: subsse.f,v 1.25 2009-04-03 16:38:23 georg Exp $
c
c interpolation routines
c
c contents :
c
c subroutine intp_lsqr(nn,ndiff,xx,xee)	least square interpolation
c subroutine matinv0(a,ip,hv,n,nd)	inversion of square matrix
c function intp_neville(nintp,xa,ya,x)	Neville algorithm for Lagrangian intp
c function intp_lagr(n,x,y,xe)		lagrangian interpolation
c function exxpp(nintp,nmax,x,y,xe,is)	interpolation from array
c function rbilin(z,xdelta,ydelta,flag)	bilinear interpolation in a square
c function rintq(...)			interpolation in a square
c
c revision log :
c
c 29.06.1998	ggu	error check in exxqq9
c 14.07.1998	ggu	new routine exxqq
c 21.08.1998	ggu	xv renamed to xx
c 31.05.2001	ggu	revised exxp (style) and exxpp (error fix and grade)
c 12.02.2003	ggu	return constant if t is out of range in exxqq
c 11.03.2005	ggu	new interpolation also for 3D time series
c 16.05.2005	ggu	more debug and new error messages
c 01.02.2006	ggu	new routine exfpres() for pointer to results
c 17.02.2006	ggu	set iunit to -1 if no file given (exffil)
c 29.02.2008	ggu	name change: exxp -> intp_lagr, extrp -> intp_lsqr
c 29.02.2008	ggu	deleted exxqq, exxqq9
c 17.03.2008	ggu	completely restructured
c 18.04.2008	ggu	new exffild, adjusted exffil, bugfix in exfini, exfintp
c 03.09.2008	ggu	bug fix in intp_ts, bug fix in exfini (aux array)
c 08.10.2008	ggu	introduced matinv0 for independence
c 08.11.2008	ggu	better error handling
c 02.04.2009	ggu	if less data given lower interpolation (REDINT)
c 03.04.2009	ggu	new routine intp_neville() (stable lagrange interpol.)
c 18.05.2014	ggu	new routines d_intp_neville(), rd_intp_neville()
c
c**************************************************
c
	subroutine intp_lsqr(nn,ndiff,xx,xee)
c
c interpolation by least square method
c
c nn		number of points used (4=cubic)
c ndiff		number of grades the interpolating polynome
c		...is lowered (normally 1)
c xx		on entry  : x-values used for interpolation
c		on return : c-values, that have to be used
c		...in computing interpolated y-value -->
c		...ye=c(i)*y(i) ; i=1,n
c xee		x-value for which y-value has to be computed
c
	parameter(iexpmx=5,idim=iexpmx+1)
	real x(idim,idim),a(idim,idim)
	real b(idim,idim)
	real xx(idim),xe(idim)
	integer ip(idim)
c
	if(nn.le.ndiff) goto 99
	if(ndiff.le.0) goto 98
	if(nn.gt.iexpmx+1) goto 97
c
	n=nn-ndiff
c
	do ii=1,nn
	x(ii,1)=1.
	do i=2,n
	x(ii,i)=x(ii,i-1)*xx(ii)
	end do
	end do
c
	do i=1,n
	do j=1,n
	h=0.
	do ii=1,nn
	h=h+x(ii,j)*x(ii,i)
	end do
	a(i,j)=h
	end do
	end do
c
	call matinv0(a,ip,xe,n,idim)	!xx,xe aux vectors
c
	do i=1,n
	do ii=1,nn
	h=0.
	do j=1,n
	h=h+a(i,j)*x(ii,j)
	end do
	b(i,ii)=h
	end do
	end do
c
	xe(1)=1.
	do i=2,n
	xe(i)=xe(i-1)*xee
	end do
c
	do ii=1,nn
	h=0.
	do i=1,n
	h=h+b(i,ii)*xe(i)
	end do
	xx(ii)=h
	end do
c
	return
c
   97	continue
	write(6,*) 'dimension error'
	write(6,*) 'iexpmx,nn :',iexpmx,nn
	stop 'error stop : intp_lsqr'
c
   98	continue
	write(6,*) 'ndiff must be greater than 0'
	write(6,*) 'ndiff =',ndiff
	stop 'error stop : intp_lsqr'
c
   99	continue
	write(6,*) 'nn must be greater than ndiff'
	write(6,*) 'nn,ndiff =',nn,ndiff
	stop 'error stop : intp_lsqr'
c
	end

c**********************************************

        subroutine matinv0(a,ip,hv,n,nd)

c inversion of square matrix (no band matrix) (original in subssm.f)
c
c a             square matrix
c ip            aux vector for pivot search
c hv            aux vector for resolution of matrix
c n             actual dimension of matrix
c nd            formal dimension of matrix

        parameter (eps=1.e-30)
        real a(nd,nd),hv(nd)
        integer ip(nd)

        do j=1,n
          ip(j)=j
        end do

        do j=1,n

c look for pivot

        amax=abs(a(j,j))
        ir=j
        do i=j+1,n
          if(abs(a(i,j)).gt.amax) then
            amax=abs(a(i,j))
            ir=i
          end if
        end do
        if(amax.lt.eps) goto 99

c change row

        if(ir.gt.j) then
          do k=1,n
            hr=a(j,k)
            a(j,k)=a(ir,k)
            a(ir,k)=hr
          end do
          ihi=ip(j)
          ip(j)=ip(ir)
          ip(ir)=ihi
        end if

c transformation

        hr=1./a(j,j)
        do i=1,n
          a(i,j)=hr*a(i,j)
        end do
        a(j,j)=hr
        do k=1,n
          if(k.ne.j) then
            do i=1,n
              if(i.ne.j) a(i,k)=a(i,k)-a(i,j)*a(j,k)
            end do
            a(j,k)=-hr*a(j,k)
          end if
        end do

        end do

c change column

        do i=1,n
          do k=1,n
            ipk=ip(k)
            hv(ipk)=a(i,k)
          end do
          do k=1,n
            a(i,k)=hv(k)
          end do
        end do

        return

   99   continue
        write(6,*) 'matrix singular. cannot invert matrix'
        write(6,*) 'i,j,amax(i,j) :',ir,j,amax
        stop 'error stop : matinv0'

        end

c**********************************************
c**********************************************
c**********************************************

        function intp_neville(nintp,xa,ya,x)

c use Neville algorithm for Lagrangian interpolation
c
c use nintp=4 for cubic interpolation
c use nintp=2 for linear interpolation

        implicit none

        real intp_neville			!interpolated value (return)
        integer nintp				!grade of interpolation
        real xa(0:nintp-1), ya(0:nintp-1)	!points to use
        real x					!where to interpolate

        integer i,k,n
        double precision xl,xh
        double precision p(0:nintp-1)

        n = nintp - 1
	p = ya

        do k=1,n
          do i=n,k,-1
            xl = xa(i-k)
            xh = xa(i)
            p(i) = ( (x-xl)*p(i) - (x-xh)*p(i-1) ) / (xh-xl)
          end do
        end do

        intp_neville = p(n)

        end

c**********************************************

        function d_intp_neville(nintp,xa,ya,x)

c use Neville algorithm for Lagrangian interpolation (double precision)
c
c use nintp=4 for cubic interpolation
c use nintp=2 for linear interpolation

        implicit none

        double precision d_intp_neville		!interpolated value (return)
        integer nintp				!grade of interpolation
        double precision xa(0:nintp-1)		!x points to use
        double precision ya(0:nintp-1)		!y points to use
        double precision x			!x value where to interpolate

        integer i,k,n
        double precision xl,xh
        double precision p(0:nintp-1)

        n = nintp - 1
	p = ya

        do k=1,n
          do i=n,k,-1
            xl = xa(i-k)
            xh = xa(i)
            p(i) = ( (x-xl)*p(i) - (x-xh)*p(i-1) ) / (xh-xl)
          end do
        end do

        d_intp_neville = p(n)

        end

c**********************************************

        function rd_intp_neville(nintp,xa,ya,x)

c use Neville algorithm for Lagrangian interpolation (double precision)
c
c use nintp=4 for cubic interpolation
c use nintp=2 for linear interpolation

        implicit none

        double precision rd_intp_neville	!interpolated value (return)
        integer nintp				!grade of interpolation
        double precision xa(0:nintp-1)		!x points to use
        real ya(0:nintp-1)			!y points to use
        double precision x			!x value where to interpolate

        integer i,k,n
        double precision xl,xh
        double precision p(0:nintp-1)

        n = nintp - 1
	p = ya

        do k=1,n
          do i=n,k,-1
            xl = xa(i-k)
            xh = xa(i)
            p(i) = ( (x-xl)*p(i) - (x-xh)*p(i-1) ) / (xh-xl)
          end do
        end do

        rd_intp_neville = p(n)

        end

c**********************************************
c**********************************************
c**********************************************

	function intp_lagr(nintp,x,y,xe)

c interpolation routine (Lagrangian interpolation)
c
c this algorithm is unstable ... please use intp_neville()
c
c intp_lagr	interpolated y-value for xe
c nintp		number of points used for interpolation (4=cubic)
c x,y		x/y-values (vector) ( x(i) != x(j) for i != j )
c xe		x-value for which y-value has to computed

	implicit none

	real intp_lagr
	integer nintp
	real x(nintp),y(nintp)
	real xe

	integer i,ii
	double precision f,g

	f=0.
	do i=1,nintp
	  g=1.
	  do ii=1,nintp
	    if(i.ne.ii) then
		g=g*(xe-x(ii))/(x(i)-x(ii))
	    end if
	  end do
	  f=f+y(i)*g
	end do

	intp_lagr=f

	end

c*************************************************************

	function exxpp(nintp,nmax,x,y,xe,iact)

c interpolation from array
c
c from given values x,y a value ye corresponding
c to xe is interpolated. a cubic interpolation is used
c
c the program is looking for the closest x-value
c only in foreward direction. for this reason
c xe-values must be passed in an increasing sequence
c
c the program does not check, if the value of xe
c is in the bounds of x(1) - x(nmax)
c eventually an extrapolated value is returned in exxpp

	implicit none

	real exxpp	!extrapolated values
	integer nintp	!number of values to use (4 for cubic, 2 for linear)
	integer nmax	!length of data arrays
        real x(1),y(1)	!data arrays
	real xe		!x-value for which y-value has to be interpolated
	integer iact	!element closest to xe (of last call on entry)
			!must be 0 for initialization

	logical bdebug
	integer nanf,nend,i
	integer min,max
	real xlow,xhigh
	real ye

	real intp_neville

	bdebug = .true.
	bdebug = .false.

c----------------------------------------------------------
c start searching from first element in x
c----------------------------------------------------------

        if( iact .le. 0 ) iact=1

c----------------------------------------------------------
c find to xe closest x-value in vector x
c----------------------------------------------------------

	do while( iact .lt. nmax )
          xlow  = abs(x(iact)-xe)
          xhigh = abs(x(iact+1)-xe)
          if( xhigh .ge. xlow ) goto 1
          iact = iact + 1
	end do
    1   continue

c----------------------------------------------------------
c x(iact) is closest value to xe ...now get closest points around xe
c----------------------------------------------------------

	if( mod(nintp,2) .eq. 0 ) then	!even
		max = nintp / 2
		min = max - 1
	else
		max = nintp / 2
		min = max
	end if

        if( x(iact) .gt. xe ) then
                nanf=iact-max
                nend=iact+min
        else
                nanf=iact-min
                nend=iact+max
        end if

c----------------------------------------------------------
c handling for the beginning or the end of array x
c----------------------------------------------------------

        if( nanf .lt. 1 ) then
                nanf=1
                nend=nintp
        else if(nend.gt.nmax) then
                nanf=nmax-nintp+1
                nend=nmax
        end if

c----------------------------------------------------------
c interpolation
c----------------------------------------------------------

        ye=intp_neville(nintp,x(nanf),y(nanf),xe)

c----------------------------------------------------------
c debug
c----------------------------------------------------------

	if( bdebug ) then
	  write(6,*) '-------------------------------- debug exxpp'
	  write(6,*) iact,nanf,nend,nintp,nmax
	  write(6,*) (x(i),i=nanf,nend)
	  write(6,*) (y(i),i=nanf,nend)
	  write(6,*) xe,ye
	  write(6,*) '--------------------------------'
	end if

c----------------------------------------------------------
c in ye is interpolated value
c----------------------------------------------------------

        exxpp=ye
	
        end

c***************************************************************
c***************************************************************
c***************************************************************

	function rbilin(z,xdelta,ydelta,flag)

c bilinear interpolation in a square
c
c z 		vector containing the values at the nodes
c xdelta	relative x-coordinate (0...1) of interpolation point in square
c ydelta	   "     y-    "         "            "          "        "
c		...(e.g. : xdelta=0.5=ydelta is centre of square)
c flag		value at node that has not to be used for interpolation
c		...these values are extrapolated from the other nodes
c		...in case all 4 nodes equal to flag ==> rbilin = flag
c rbilin	interpolated value at return
c
c numeration of square
c			  (0,1)            (1,1)
c				+--------+
c				| 3    4 |
c				|        |
c				| 1    2 |
c				+--------+
c			  (0,0)            (1,0)
c
c diagonal sum = 5
c
c formula for interpolation :
c
c	z = a + b*x + c*y + d*x*y
c
c	a = z(0,0)
c	b = z(1,0) - z(0,0)
c	c = z(0,1) - z(0,0)
c	d = z(1,1) + z(0,0) - z(1,0) - z(0,1)

	implicit none

	real rbilin
	real z(4)
	real xdelta,ydelta
	real flag

	integer nout,nin,i,idiag,iih
	integer iout(4),iin(4)
	real zh(4),zhh,ztot
	real a,b,c,d

c----------------------------------------------------------------
c get inner nodes (with value) and outer nodes (without value)
c----------------------------------------------------------------

	nout=0
	nin=0
	ztot=0.
	do i=1,4
	  zhh = z(i)
	  zh(i) = zhh
	  if(zhh.eq.flag) then
		nout=nout+1
		iout(nout)=i
	  else
		nin=nin+1
		iin(nin)=i
		ztot=ztot+zhh
	  end if
	end do

c----------------------------------------------------------------
c extrapolate good (inner) to bad (outer) nodes
c----------------------------------------------------------------

	if(nout.eq.0) then	!no outer points
c		nothing
	else if(nin.eq.0) then	!no inner point
		rbilin=flag
		return
	else if(nin.eq.1) then	!only 1 inner point
		do i=1,4
		  zh(i)=ztot
		end do
	else if(nout.eq.1) then	!only 1 outer point
		iih=iout(1)
		idiag=5-iih
		zh(iih)=ztot-2.*zh(idiag)	!extrapolation from inner
						!...triangel to ext. point
	else			!2 inner/outer points
		if(iin(1)+iin(2).eq.5) then	!diagonal out of area
			zhh=ztot*0.5
			zh(iout(1))=zhh
			zh(iout(2))=zhh
		else				!side out of area
			iih=5-iin(2)		!to find point to be extrapol.
			zh(iih)=zh(iin(1))	!...get second inner point
			iih=5-iin(1)		!...and go to the diagonal
			zh(iih)=zh(iin(2))	!...
		end if
	end if

c----------------------------------------------------------------
c interpolation in square
c----------------------------------------------------------------

	a=zh(1)
	b=zh(2)-a
	c=zh(3)-a
	d=zh(4)-a-b-c

	rbilin = a + b*xdelta + c*ydelta + d*xdelta*ydelta

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c***************************************************************

	function rintq(z,idim,jdim,ipos,jpos,xdelta,ydelta,flag,ier)

c interpolation in a square
c
c z 		matrix containing the values at the nodes
c (idim,jdim)	dimension of z
c ipos,jpos	position of local node (0,0) in z
c xdelta	relative x-coordinate (0...1) of interpolation point in square
c ydelta	   "     y-    "         "            "          "        "
c		...(e.g. : xdelta=0.5=ydelta is centre of square)
c flag		value at node that has not to be used for interpolation
c		...these values are extrapolated from the other nodes
c		...in case all 4 nodes equal to flag ==> rintq = flag
c ier		error status (return value)
c		... 0 : ok
c		... 1 : coordinates have been adjusted because out of grid
c rintq		interpolated value at return
c
c uses rbilin to do bilinear interpolation in square
c
c numeration of square
c			  (0,1)            (1,1)
c				+--------+
c				| 3    4 |
c				|        |
c				| 1    2 |
c				+--------+
c			  (0,0)            (1,0)
c
c-----------------------------------------------------------------------

	implicit none

	real rintq
	integer idim,jdim
	integer ipos,jpos
	integer ier
	real xdelta,ydelta
	real flag
	real z(idim,jdim)

	integer i0,j0,i
	real xd,yd
	real zh(4)
	real rbilin

	integer iv(4),jv(4)
	data iv,jv /0,1,0,1,0,0,1,1/	!translates matrix into vector coord.

	i0=ipos
	j0=jpos
	xd=xdelta
	yd=ydelta

c----------------------------------------------------------------
c adjust coordinates out of grid
c----------------------------------------------------------------

	ier=0
	if(i0.lt.1) then
		ier=1
		i0=1
		xd=0.
	else if(i0.ge.idim) then
		ier=1
		i0=idim-1
		xd=1.
	end if
	if(j0.lt.1) then
		ier=1
		j0=1
		yd=0.
	else if(j0.ge.jdim) then
		ier=1
		j0=jdim-1
		yd=1.
	end if

c----------------------------------------------------------------
c copy values to array
c----------------------------------------------------------------

	do i=1,4
	  zh(i)=z(i0+iv(i),j0+jv(i))
	end do

c----------------------------------------------------------------
c bilinear interpolation
c----------------------------------------------------------------

	rintq = rbilin(zh,xd,yd,flag)

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c***************************************************************

