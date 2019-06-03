
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c spline routines
c
c contents :
c
c subroutine spline(n,x,y,y2)		prepares cubic spline
c subroutine splint(n,xa,ya,y2a,x,y)	evaluates spline -> gives back y(x)
c subroutine spltst			test spline
c
c revision log :
c
c 06.03.1999	ggu	routines written from scratch (Numerical receipes)
c 23.03.2010	ggu	changed v6.1.1
c 18.12.2018	ggu	changed VERS_7_5_52
c 16.02.2019	ggu	changed VERS_7_5_60
c
c***************************************************************

      subroutine spline(n,x,y,y2)

c prepares cubic spline

      implicit none

      integer n		!total number of points in arrays
      real x(n)		!x-values in ascending order
      real y(n)		!y-values
      real y2(n)	!values of second derivative on return

      integer i
      real sig,p
      real u(n)		!auxiliary array

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
      real xa(n)	!x-values in ascending order
      real ya(n)	!y-values
      real y2a(n)	!values of second derivative
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
c***************************************************************
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
	
	call spline(n,x,y,y2)

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

