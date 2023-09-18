!
! $Id: subsse.f,v 1.25 2009-04-03 16:38:23 georg Exp $
!
! interpolation routines
!
! contents :
!
! subroutine intp_lsqr(nn,ndiff,xx,xee)	least square interpolation
! subroutine matinv0(a,ip,hv,n,nd)	inversion of square matrix
! function intp_neville(nintp,xa,ya,x)	Neville algorithm for Lagrangian intp
! function intp_lagr(n,x,y,xe)		lagrangian interpolation
! function exxpp(nintp,nmax,x,y,xe,is)	interpolation from array
! function rbilin(z,xdelta,ydelta,flag)	bilinear interpolation in a square
! function rintq(...)			interpolation in a square
!
! revision log :
!
! 29.06.1998	ggu	error check in exxqq9
! 14.07.1998	ggu	new routine exxqq
! 21.08.1998	ggu	xv renamed to xx
! 31.05.2001	ggu	revised exxp (style) and exxpp (error fix and grade)
! 12.02.2003	ggu	return constant if t is out of range in exxqq
! 11.03.2005	ggu	new interpolation also for 3D time series
! 16.05.2005	ggu	more debug and new error messages
! 01.02.2006	ggu	new routine exfpres() for pointer to results
! 17.02.2006	ggu	set iunit to -1 if no file given (exffil)
! 29.02.2008	ggu	name change: exxp -> intp_lagr, extrp -> intp_lsqr
! 29.02.2008	ggu	deleted exxqq, exxqq9
! 17.03.2008	ggu	completely restructured
! 18.04.2008	ggu	new exffild, adjusted exffil, bugfix in exfini, exfintp
! 03.09.2008	ggu	bug fix in intp_ts, bug fix in exfini (aux array)
! 08.10.2008	ggu	introduced matinv0 for independence
! 08.11.2008	ggu	better error handling
! 02.04.2009	ggu	if less data given lower interpolation (REDINT)
! 03.04.2009	ggu	new routine intp_neville() (stable lagrange interpol.)
! 18.05.2014	ggu	new routines d_intp_neville(), rd_intp_neville()
!
!**************************************************

!---------------------------------------------------------
        module intp
!---------------------------------------------------------
        contains
!---------------------------------------------------------

	subroutine intp_lsqr(nn,ndiff,xx,xee)
!
! interpolation by least square method
!
! nn		number of points used (4=cubic)
! ndiff		number of grades the interpolating polynome
!		...is lowered (normally 1)
! xx		on entry  : x-values used for interpolation
!		on return : c-values, that have to be used
!		...in computing interpolated y-value -->
!		...ye=c(i)*y(i) ; i=1,n
! xee		x-value for which y-value has to be computed
!
	parameter(iexpmx=5,idim=iexpmx+1)
	double precision x(idim,idim),a(idim,idim)
	double precision b(idim,idim)
	double precision xx(idim),xe(idim)
	integer ip(idim)
!
	if(nn.le.ndiff) goto 99
	if(ndiff.le.0) goto 98
	if(nn.gt.iexpmx+1) goto 97
!
	n=nn-ndiff
!
	do ii=1,nn
	x(ii,1)=1.
	do i=2,n
	x(ii,i)=x(ii,i-1)*xx(ii)
	end do
	end do
!
	do i=1,n
	do j=1,n
	h=0.
	do ii=1,nn
	h=h+x(ii,j)*x(ii,i)
	end do
	a(i,j)=h
	end do
	end do
!
	call matinv0(a,ip,xe,n,idim)	!xx,xe aux vectors
!
	do i=1,n
	do ii=1,nn
	h=0.
	do j=1,n
	h=h+a(i,j)*x(ii,j)
	end do
	b(i,ii)=h
	end do
	end do
!
	xe(1)=1.
	do i=2,n
	xe(i)=xe(i-1)*xee
	end do
!
	do ii=1,nn
	h=0.
	do i=1,n
	h=h+b(i,ii)*xe(i)
	end do
	xx(ii)=h
	end do
!
	return
!
   97	continue
	write(6,*) 'dimension error'
	write(6,*) 'iexpmx,nn :',iexpmx,nn
	stop 'error stop : intp_lsqr'
!
   98	continue
	write(6,*) 'ndiff must be greater than 0'
	write(6,*) 'ndiff =',ndiff
	stop 'error stop : intp_lsqr'
!
   99	continue
	write(6,*) 'nn must be greater than ndiff'
	write(6,*) 'nn,ndiff =',nn,ndiff
	stop 'error stop : intp_lsqr'
!
	end

!**********************************************

        subroutine matinv0(a,ip,hv,n,nd)

! inversion of square matrix (no band matrix) (original in subssm.f)
!
! a             square matrix
! ip            aux vector for pivot search
! hv            aux vector for resolution of matrix
! n             actual dimension of matrix
! nd            formal dimension of matrix

        parameter (eps=1.e-30)
        double precision a(nd,nd),hv(nd)
        integer ip(nd)

        do j=1,n
          ip(j)=j
        end do

        do j=1,n

! look for pivot

        amax=abs(a(j,j))
        ir=j
        do i=j+1,n
          if(abs(a(i,j)).gt.amax) then
            amax=abs(a(i,j))
            ir=i
          end if
        end do
        if(amax.lt.eps) goto 99

! change row

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

! transformation

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

! change column

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

!**********************************************
!**********************************************
!**********************************************

        function intp_neville(nintp,xa,ya,x)

! use Neville algorithm for Lagrangian interpolation
!
! use nintp=4 for cubic interpolation
! use nintp=2 for linear interpolation

        implicit none

        double precision intp_neville			!interpolated value (return)
        integer nintp				!grade of interpolation
        double precision xa(0:nintp-1), ya(0:nintp-1)	!points to use
        double precision x					!where to interpolate

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

!**********************************************

        function d_intp_neville(nintp,xa,ya,x)

! use Neville algorithm for Lagrangian interpolation (double precision)
!
! use nintp=4 for cubic interpolation
! use nintp=2 for linear interpolation

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

!**********************************************

        function rd_intp_neville(nintp,xa,ya,x)

! use Neville algorithm for Lagrangian interpolation (double precision)
!
! use nintp=4 for cubic interpolation
! use nintp=2 for linear interpolation

        implicit none

        double precision rd_intp_neville        !interpolated value (return)
        integer nintp	                        !grade of interpolation
        double precision xa(0:nintp-1)          !x points to use
!        double precision ya(0:nintp-1)                      !y points to use
        double precision ya(0:nintp-1)          !y points to use (ivb)
        double precision x                      !x value where to interpolate

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

!**********************************************
!**********************************************
!**********************************************

	function intp_lagr(nintp,x,y,xe)

! interpolation routine (Lagrangian interpolation)
!
! this algorithm is unstable ... please use intp_neville()
!
! intp_lagr	interpolated y-value for xe
! nintp		number of points used for interpolation (4=cubic)
! x,y		x/y-values (vector) ( x(i) != x(j) for i != j )
! xe		x-value for which y-value has to computed

	implicit none

	double precision intp_lagr
	integer nintp
	double precision x(nintp),y(nintp)
	double precision xe

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

!*************************************************************

	function exxpp(nintp,nmax,x,y,xe,iact)

! interpolation from array
!
! from given values x,y a value ye corresponding
! to xe is interpolated. a cubic interpolation is used
!
! the program is looking for the closest x-value
! only in foreward direction. for this reason
! xe-values must be passed in an increasing sequence
!
! the program does not check, if the value of xe
! is in the bounds of x(1) - x(nmax)
! eventually an extrapolated value is returned in exxpp

	implicit none

	double precision exxpp	!extrapolated values
	integer nintp	!number of values to use (4 for cubic, 2 for linear)
	integer nmax	!length of data arrays
        double precision x(1),y(1)	!data arrays
	double precision xe		!x-value for which y-value has to be interpolated
	integer iact	!element closest to xe (of last call on entry)
			!must be 0 for initialization

	logical bdebug
	integer nanf,nend,i
	integer min,max
	double precision xlow,xhigh
	double precision ye

	bdebug = .true.
	bdebug = .false.

!----------------------------------------------------------
! start searching from first element in x
!----------------------------------------------------------

        if( iact .le. 0 ) iact=1

!----------------------------------------------------------
! find to xe closest x-value in vector x
!----------------------------------------------------------

	do while( iact .lt. nmax )
          xlow  = abs(x(iact)-xe)
          xhigh = abs(x(iact+1)-xe)
          if( xhigh .ge. xlow ) goto 1
          iact = iact + 1
	end do
    1   continue

!----------------------------------------------------------
! x(iact) is closest value to xe ...now get closest points around xe
!----------------------------------------------------------

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

!----------------------------------------------------------
! handling for the beginning or the end of array x
!----------------------------------------------------------

        if( nanf .lt. 1 ) then
                nanf=1
                nend=nintp
        else if(nend.gt.nmax) then
                nanf=nmax-nintp+1
                nend=nmax
        end if

!----------------------------------------------------------
! interpolation
!----------------------------------------------------------

        ye=intp_neville(nintp,x(nanf),y(nanf),xe)

!----------------------------------------------------------
! debug
!----------------------------------------------------------

	if( bdebug ) then
	  write(6,*) '-------------------------------- debug exxpp'
	  write(6,*) iact,nanf,nend,nintp,nmax
	  write(6,*) (x(i),i=nanf,nend)
	  write(6,*) (y(i),i=nanf,nend)
	  write(6,*) xe,ye
	  write(6,*) '--------------------------------'
	end if

!----------------------------------------------------------
! in ye is interpolated value
!----------------------------------------------------------

        exxpp=ye
	
        end

!***************************************************************
!***************************************************************
!***************************************************************

	function rbilin(z,xdelta,ydelta,flag)

! bilinear interpolation in a square
!
! z 		vector containing the values at the nodes
! xdelta	relative x-coordinate (0...1) of interpolation point in square
! ydelta	   "     y-    "         "            "          "        "
!		...(e.g. : xdelta=0.5=ydelta is centre of square)
! flag		value at node that has not to be used for interpolation
!		...these values are extrapolated from the other nodes
!		...in case all 4 nodes equal to flag ==> rbilin = flag
! rbilin	interpolated value at return
!
! numeration of square
!			  (0,1)            (1,1)
!				+--------+
!				| 3    4 |
!				|        |
!				| 1    2 |
!				+--------+
!			  (0,0)            (1,0)
!
! diagonal sum = 5
!
! formula for interpolation :
!
!	z = a + b*x + c*y + d*x*y
!
!	a = z(0,0)
!	b = z(1,0) - z(0,0)
!	c = z(0,1) - z(0,0)
!	d = z(1,1) + z(0,0) - z(1,0) - z(0,1)

	implicit none

	double precision rbilin
	double precision z(4)
	double precision xdelta,ydelta
	double precision flag

	integer nout,nin,i,idiag,iih
	integer iout(4),iin(4)
	double precision zh(4),zhh,ztot
	double precision a,b,c,d

!----------------------------------------------------------------
! get inner nodes (with value) and outer nodes (without value)
!----------------------------------------------------------------

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

!----------------------------------------------------------------
! extrapolate good (inner) to bad (outer) nodes
!----------------------------------------------------------------

	if(nout.eq.0) then	!no outer points
!		nothing
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

!----------------------------------------------------------------
! interpolation in square
!----------------------------------------------------------------

	a=zh(1)
	b=zh(2)-a
	c=zh(3)-a
	d=zh(4)-a-b-c

	rbilin = a + b*xdelta + c*ydelta + d*xdelta*ydelta

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!***************************************************************

	function rintq(z,idim,jdim,ipos,jpos,xdelta,ydelta,flag,ier)

! interpolation in a square
!
! z 		matrix containing the values at the nodes
! (idim,jdim)	dimension of z
! ipos,jpos	position of local node (0,0) in z
! xdelta	relative x-coordinate (0...1) of interpolation point in square
! ydelta	   "     y-    "         "            "          "        "
!		...(e.g. : xdelta=0.5=ydelta is centre of square)
! flag		value at node that has not to be used for interpolation
!		...these values are extrapolated from the other nodes
!		...in case all 4 nodes equal to flag ==> rintq = flag
! ier		error status (return value)
!		... 0 : ok
!		... 1 : coordinates have been adjusted because out of grid
! rintq		interpolated value at return
!
! uses rbilin to do bilinear interpolation in square
!
! numeration of square
!			  (0,1)            (1,1)
!				+--------+
!				| 3    4 |
!				|        |
!				| 1    2 |
!				+--------+
!			  (0,0)            (1,0)
!
!-----------------------------------------------------------------------

	implicit none

	double precision rintq
	integer idim,jdim
	integer ipos,jpos
	integer ier
	double precision xdelta,ydelta
	double precision flag
	double precision z(idim,jdim)

	integer i0,j0,i
	double precision xd,yd
	double precision zh(4)

	integer iv(4),jv(4)
	data iv,jv /0,1,0,1,0,0,1,1/	!translates matrix into vector coord.

	i0=ipos
	j0=jpos
	xd=xdelta
	yd=ydelta

!----------------------------------------------------------------
! adjust coordinates out of grid
!----------------------------------------------------------------

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

!----------------------------------------------------------------
! copy values to array
!----------------------------------------------------------------

	do i=1,4
	  zh(i)=z(i0+iv(i),j0+jv(i))
	end do

!----------------------------------------------------------------
! bilinear interpolation
!----------------------------------------------------------------

	rintq = rbilin(zh,xd,yd,flag)

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!***************************************************************

!---------------------------------------------------------
        end module intp
!---------------------------------------------------------
