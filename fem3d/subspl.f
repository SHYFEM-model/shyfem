c
c $Id: subspl.f,v 1.2 2007-06-08 11:22:00 georg Exp $
c
c spline routines
c
c contents :
c
c subroutine spline(n,x,y,y2,u)		prepares cubic spline
c subroutine splint(n,xa,ya,y2a,x,y)	evaluates spline -> gives back y(x)
c subroutine spltst			test spline
c
c revision log :
c
c 06.03.1999    ggu     routines written from scratch (Numerical receipes)
c
c***************************************************************

      subroutine spline(n,x,y,y2,u)

c prepares cubic spline

      implicit none

      integer n		!total number of points in arrays
      real x(1)		!x-values in ascending order
      real y(1)		!y-values
      real y2(1)	!values of second derivative on return
      real u(1)		!auxiliary array

      integer i
      real sig,p

c     we treat only natural boundary conditions

      y2(1)=0.
      y2(n)=0.
      u(1)=0.

c     decomposition of tridiagonal loop

      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do

c     back substitution

      do i=n-1,1,-1
        y2(i)=y2(i)*y2(i+1)+u(i)
      end do

      end

c***************************************************************

      subroutine splint(n,xa,ya,y2a,x,y)

c evaluates spline -> gives back y(x)

      implicit none

      integer n		!total number of points
      real xa(1)	!x-values in ascending order
      real ya(1)	!y-values
      real y2a(1)	!values of second derivative
      real x		!x-value where spline has to be evaluated
      real y		!y-value to given x -> y(x)

      integer klo,khi,k
      real a,b,h

c     search for x by bisection

      klo=1
      khi=n

      do while (khi-klo.gt.1)
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      end do

      h=xa(khi)-xa(klo)
      if (h.le.0.) stop 'error stop splint: bad xa array'

      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

      end

c***************************************************************

	subroutine spltst

c test spline

	implicit none

	integer ndim
	parameter (ndim = 8)

	real x(ndim), y(ndim)
	real y2(ndim), aux(ndim)

	integer i,n
	real xa,ya

	data x /0.,1.,2.,4.,6.,7.,9.,12./
	data y /0.,1.,3.,4.,2.,0.,2.,5./

	n = ndim

	do i=1,n
	  write(70,'(2f12.4)') x(i),y(i)
	end do
	
	call spline(n,x,y,y2,aux)

	do i=0,12*4
	  xa = i * 0.25
	  if( xa .gt. 12. ) xa = 12.
	  
	  call splint(n,x,y,y2,xa,ya)

	  write(71,'(2f12.4)') xa,ya
	end do

	end
	  
c***************************************************************

c	program spltest
c	call spltst
c	end

c***************************************************************

