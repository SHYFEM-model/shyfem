
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

c routines for bilinear interpolation in quadrilateral
c
c********************************************************************
c
c numeration of square
c
c                         (0,1)            (1,1)
c                               +--------+
c                               | 4    3 |
c                               |        |
c                               | 1    2 |
c                               +--------+
c                         (0,0)            (1,0)
c
c please note difference in numbering from old routines
c numbering must be in anti-clockwise sense
c
c********************************************************************

	subroutine bilin_intp(xq,yq,vals,x,y,val)

c performs bilinear interpolation in quadrilateral
c coordinates can be arbitrary, but must be convex
c numbering of nodes counterclockwise

	real xq(4)	!x coordinates of quadrilateral
	real yq(4)	!y coordinates of quadrilateral
	real vals(4)	!values on corners
	real x,y	!x/y coordinates of interpolation point (inside quad)
	real val	!computed value

	real u,v

	call bilin_inverse(xq,yq,u,v,x,y)	!find u,v
	call bilin_intp_square(vals,u,v,val)	!interpolate

	end

c********************************************************************

	subroutine bilin_intp_square(vals,u,v,val)

c performs bilinear interpolation in square [0-1]
c numbering of nodes counterclockwise

	real vals(4)	!values on corners
	real u,v	!relative coordinates in square [0-1]
	real val	!computed value

	double precision a0,a1,a2,a3

	a0 = vals(1)
	a1 = vals(2) - vals(1)
	a2 = vals(4) - vals(1)
	a3 = vals(1) + vals(3) - vals(2) - vals(4)

	val = a0 + a1*u + a2*v + a3*u*v

	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine bilin_direct(xq,yq,u,v,x,y)

c direct transformation from unit square to arbitrary quadrilateral

	implicit none

	real xq(4)	!x coordinates of quadrilateral
	real yq(4)	!y coordinates of quadrilateral
	real u,v	!relative coordinates in square [0-1]
	real x,y	!computed x/y coordinates in quadrilateral

	double precision a0,a1,a2,a3
	double precision b0,b1,b2,b3

	a0 = xq(1)
	b0 = yq(1)
	a1 = xq(2) - xq(1)
	b1 = yq(2) - yq(1)
	a2 = xq(4) - xq(1)
	b2 = yq(4) - yq(1)
	a3 = xq(1) + xq(3) - xq(2) - xq(4)
	b3 = yq(1) + yq(3) - yq(2) - yq(4)

	x = a0 + a1*u + a2*v + a3*u*v
	y = b0 + b1*u + b2*v + b3*u*v

	end

c********************************************************************

	subroutine bilin_inverse(xq,yq,u,v,x,y)

c inverse transformation from unit square to arbitrary quadrilateral
c computes relative coordinates in square

	implicit none

	real xq(4)	!x coordinates of quadrilateral
	real yq(4)	!y coordinates of quadrilateral
	real u,v	!computed relative coordinates in square [0-1]
	real x,y	!x/y coordinates in quadrilateral

	double precision a0,a1,a2,a3
	double precision b0,b1,b2,b3
	double precision aa,cc1,cc,bb1,bb,dd
	double precision aa_1,aa_2
	double precision uu,vv
	double precision uu1,vv1,uu2,vv2
	double precision xx1,yy1,xx2,yy2,dd1,dd2

	double precision eps

	eps = 1.e-6
	eps = 5.e-6

	a0 = xq(1)
	b0 = yq(1)
	a1 = xq(2) - xq(1)
	b1 = yq(2) - yq(1)
	a2 = xq(4) - xq(1)
	b2 = yq(4) - yq(1)
	a3 = xq(1) + xq(3) - xq(2) - xq(4)
	b3 = yq(1) + yq(3) - yq(2) - yq(4)

	aa_1 = b2*a3 - b3*a2
	aa_2 = b1*a3 - b3*a1

	aa = aa_1
	cc1 = b0*a1 - b1*a0
	cc = cc1 + b1*x - a1*y
	bb1 = b0*a3 - b3*a0 + b2*a1 - b1*a2
	bb = bb1 + b3*x - a3*y
	dd = sqrt(bb*bb - 4*aa*cc)

	vv1 = (-bb + dd) / (2*aa)
	uu1 = (x - a0 - a2*vv1)/(a1+a3*vv1)

	aa = aa_2
	cc1 = b0*a2 - b2*a0
	cc = cc1 + b2*x - a2*y
	bb1 = b0*a3 - b3*a0 + b1*a2 - b2*a1
	bb = bb1 + b3*x - a3*y
	dd = sqrt(bb*bb - 4*aa*cc)

	uu2 = (-bb - dd) / (2*aa)
	vv2 = (x - a0 - a1*uu2)/(a2+a3*uu2)

c	if( uu1 .lt. 0. .or. uu1 .gt. 1. .or.
c     +			vv1 .lt. 0. .or. vv1 .gt. 1. ) then
c	  uu = uu2
c	  vv = vv2
c	else if( uu2 .lt. 0. .or. uu2 .gt. 1. .or.
c     +			vv2 .lt. 0. .or. vv2 .gt. 1. ) then
c	  uu = uu1
c	  vv = vv1
c	else		!check x/y
c	  xx1 = a0 + a1*uu1 + a2*vv1 + a3*uu1*vv1
c	  yy1 = b0 + b1*uu1 + b2*vv1 + b3*uu1*vv1
c	  xx2 = a0 + a1*uu2 + a2*vv2 + a3*uu2*vv2
c	  yy2 = b0 + b1*uu2 + b2*vv2 + b3*uu2*vv2
c	  dd1 = max(abs(xx1-x),abs(yy1-y))
c	  dd2 = max(abs(xx2-x),abs(yy2-y))
c	  write(6,*) 'diff inverse: ',dd1,dd2
c	  if( dd1 .gt. dd2 ) then
c	    uu = uu2
c	    vv = vv2
c	  else
c	    uu = uu1
c	    vv = vv1
c	  end if
c	end if

	if( abs(aa_1) .gt. abs(aa_2) ) then
	  uu = uu2
	  vv = vv2
	else
	  uu = uu1
	  vv = vv1
	end if

        uu = uu1
        vv = vv1

	if( uu .lt. 0. .and. uu .gt. -eps ) uu = 0.
	if( uu .gt. 1. .and. uu .lt. 1+eps ) uu = 1.
	if( vv .lt. 0. .and. vv .gt. -eps ) vv = 0.
	if( vv .gt. 1. .and. vv .lt. 1+eps ) vv = 1.

	u = uu
	v = vv

	end

c********************************************************************
c********************************************************************
c********************************************************************
c
c from here on testing routines
c
c********************************************************************
c********************************************************************
c********************************************************************

	subroutine bilin_test_one(xq,yq,u,v,diff)

	implicit none

	real xq(4),yq(4)
	real u,v
	real eps_warn,eps_error

	real x,y,uu,vv,diff

	eps_warn = 1.e-5
	eps_warn = 5.e-5
	eps_error = 1.e-3

	call bilin_direct(xq,yq,u,v,x,y)
	call bilin_inverse(xq,yq,uu,vv,x,y)

	!write(6,*) u,v,uu,vv,x,y

	if( u .lt. 0. .or. u .gt. 1. ) goto 99
	if( v .lt. 0. .or. v .gt. 1. ) goto 99
	if( uu .lt. 0. .or. uu .gt. 1. ) goto 99
	if( vv .lt. 0. .or. vv .gt. 1. ) goto 99

	diff = max(abs(uu-u),abs(vv-v))
	!write(6,*) u,v,uu,vv,diff
	if( diff .gt. eps_warn ) then
		call bilin_plot_grid(66,x,y,xq,yq)
		write(6,*) '        diff: ',diff
	end if
	if( diff .gt. eps_error ) goto 98

	return
   98	continue
	write(6,*) 'diff: ',diff
   99	continue
	write(6,*) 'error...'
	write(6,*) xq
	write(6,*) yq
	write(6,*) u,v
	write(6,*) uu,vv
	call bilin_plot_grid(77,x,y,xq,yq)
	stop
	end

c********************************************************************

	subroutine bilin_test

	implicit none

	integer n,i,idum,j0,j,ntot,nexe,ndiff
	real xq(4),yq(4)
	real u,v
	real fact
	real eps_warn,diff,perc,diff_max

	integer xsign(0:3)
	integer ysign(0:3)
	data xsign /-1,+1,+1,-1/
	data ysign /-1,-1,+1,+1/

	real bilin_rand
	logical bilin_check_convex

	j0 = 0
	idum = 99
	fact = 1.
	eps_warn = 5.e-5
	diff_max = 0.
	ntot = 10000000
	nexe = 0
	ndiff = 0

	do n=1,ntot
	  j0 = mod(j0+1,4)
	  do i=0,3
	    j = mod(j0+i,4)
	    j = i
	    !xq(1+i) = xsign(j)*fact*bilin_rand(idum)
	    !yq(1+i) = ysign(j)*fact*bilin_rand(idum)
	    xq(1+i) = fact*bilin_rand(idum)
	    yq(1+i) = fact*bilin_rand(idum)
	  end do
	  if(  bilin_check_convex(xq,yq) ) then
	    u = bilin_rand(idum)
	    v = bilin_rand(idum)
	    call bilin_test_one(xq,yq,u,v,diff)
	    nexe = nexe + 1
	    if( diff .gt. eps_warn ) then
	      ndiff = ndiff + 1
	      diff_max = max(diff_max,diff)
	    end if
	  end if
	  perc = (100.*n)/ntot
	  if( mod(n,ntot/100) .eq. 0 ) write(6,*) n,nexe,ntot,perc,' %'
	end do

	write(6,*) 'test finished: ',ndiff,nexe,ntot,(100.*nexe)/ntot,' %'
	write(6,*) 'differences: ',ndiff,diff_max

	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine bilin_plot_grid(iunit,xx,yy,xq,yq)

	implicit none

	integer iunit
	real xx,yy
	real xq(4), yq(4)

	integer ix,iy,n
	real x,y,u,v

	integer nn,nl
	save nn,nl
	data nn,nl /0,0/
	real x0
	save x0
	data x0 /0./

	x0 = x0 + 2.

	write(iunit,*) 1,nn+1,3,x0+xq(1),yq(1)
	write(iunit,*) 1,nn+2,3,x0+xq(2),yq(2)
	write(iunit,*) 1,nn+3,3,x0+xq(3),yq(3)
	write(iunit,*) 1,nn+4,3,x0+xq(4),yq(4)
	write(iunit,*) 1,nn+5,3,x0+xx,yy
	write(iunit,*) 3,nl+1,3,5,nn+1,nn+2,nn+3,nn+4,nn+1

	nn = nn + 5
	nl = nl + 1

	n = nn
	do ix=0,10
	  do iy=0,10
	    u = ix/10.
	    v = iy/10.
	    call bilin_direct(xq,yq,u,v,x,y)
	    n = n + 1
	    write(iunit+1,*) 1,n,4,x0+x,y
	  end do
	end do
	nn = n

	end

c********************************************************************

	function bilin_check_convex(xq,yq)

	implicit none

	logical bilin_check_convex
	real xq(4), yq(4)

	integer i,n
	real xl,yl,xm,ym,xn,yn

	logical bilin_lefton,bconv

	n = 4
	bconv = .true.

        xm = xq(n-1)
        ym = yq(n-1)
        xn = xq(n)
        yn = yq(n)

        do i=1,n
          xl = xm
          yl = ym
          xm = xn
          ym = yn
          xn = xq(i)
          yn = yq(i)

          if( .not. bilin_lefton(xl,yl,xm,ym,xn,yn) ) bconv = .false.
        end do

	bilin_check_convex = bconv

	end

c********************************************************************

        function bilin_lefton(x1,y1,x2,y2,x3,y3)

c left turn or straight ?

        implicit none

        logical bilin_lefton
        real x1,y1,x2,y2,x3,y3

        real bilin_areat

        bilin_lefton = bilin_areat(x1,y1,x2,y2,x3,y3) .ge. 0.

	end

c********************************************************************

        function bilin_areat(x1,y1,x2,y2,x3,y3)

c computes area of triangle

        implicit none

        real bilin_areat
        real x1,y1,x2,y2,x3,y3

        bilin_areat = 0.5 * ( (x2-x1) * (y3-y1) - (x3-x1) * (y2-y1) )

        end

c********************************************************************

      function bilin_rand(idum)
      parameter (m=714025,ia=1366,ic=150889,rm=1.4005112e-6)
      dimension ir(97)
      save iy,ir
      data iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        idum=mod(ic-idum,m)
        do 11 j=1,97
          idum=mod(ia*idum+ic,m)
          ir(j)=idum
11      continue
        idum=mod(ia*idum+ic,m)
        iy=idum
      endif
      j=1+(97*iy)/m
      if(j.gt.97.or.j.lt.1) stop 'error stop ran2: internal error'
      iy=ir(j)
      bilin_rand=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      return
      end

c********************************************************************

	program bilin_test_main
	call bilin_test
	end

c********************************************************************

