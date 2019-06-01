
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2014,2018  Georg Umgiesser
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

c routines to transform between color spaces

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 05.12.2014	ggu	changed VERS_7_0_8
! 18.12.2018	ggu	changed VERS_7_5_52
! 21.05.2019	ggu	changed VERS_7_5_62

c****************************************************************

	subroutine rgb2cmy(r,g,b,c,m,y)

c rgb -> cmy

	implicit none

	real r,g,b	![0-1]
	real c,m,y	![0-1]

	c = 1. - r
	m = 1. - g
	y = 1. - b

	end

c****************************************************************

	subroutine cmy2rgb(c,m,y,r,g,b)

c cmy -> rgb

	implicit none

	real c,m,y	![0-1]
	real r,g,b	![0-1]

	r = 1. - c
	g = 1. - m
	b = 1. - y

	end

c****************************************************************

	subroutine rgb2hsv(r,g,b,h,s,v)

c rgb -> hsv		note: h is not [0-360] but [0-1]

	implicit none

	real r,g,b	![0-1]
	real h,s,v	![0-1]

	real hmax			! 360 or 1
c	parameter( hmax = 360. )	! h [0-360]
	parameter( hmax = 1. )		! h [0-1]
	real hconv
	parameter( hconv = hmax / 6. )
	real undef
	parameter( undef = 0. )

	real maxv,minv
	real diff
	real rdist,gdist,bdist

	maxv = max(r,g,b)
	minv = min(r,g,b)
	diff = maxv - minv

	v = maxv
	s = 0.
	if( maxv .gt. 0. ) s = diff/maxv

	if( s .eq. 0. ) then
	  h = undef
	else
	  rdist = (maxv-r)/diff
	  gdist = (maxv-g)/diff
	  bdist = (maxv-b)/diff
	  if( r .eq. maxv ) then
	    h = bdist - gdist
	  else if( g .eq. maxv ) then
	    h = 2. + rdist - bdist
	  else if( b .eq. maxv ) then
	    h = 4. + gdist - rdist
	  else
	    stop 'error stop rgb2hsv: internal error (1)'
	  end if
	  h = h * hconv
	  if( h .lt. 0. ) h = h + hmax
	end if

	end

c****************************************************************

	subroutine hsv2rgb(h,s,v,r,g,b)

c hsv -> rgb		note: h is not [0-360] but [0-1]

	implicit none

	real hconv,hdist,hmax
c	parameter( hmax = 360. )
	parameter( hmax = 1. )
	parameter( hdist = hmax/6. )
	parameter( hconv = 1./hdist )

	real h,s,v	![0-1]
	real r,g,b	![0-1]

	integer i
	real p,q

	i = hconv * h
	i = mod(i,6)
	p = v * (h-i*hdist) / hdist	!rising
	q = v - p			!falling

	if( i .eq. 0 ) then
	    r=v
	    g=p
	    b=0
	else if( i .eq. 1 ) then
	    r=q
	    g=v
	    b=0
	else if( i .eq. 2 ) then
	    r=0
	    g=v
	    b=p
	else if( i .eq. 3 ) then
	    r=0
	    g=q
	    b=v
	else if( i .eq. 4 ) then
	    r=p
	    g=0
	    b=v
	else if( i .eq. 5 ) then
	    r=v
	    g=0
	    b=q
	else
	    stop 'error stop hsv2rgb: internal error (1)'
	end if

	r = v + (r-v) * s
	g = v + (g-v) * s
	b = v + (b-v) * s

	end

c****************************************************************

      function randomtb(idum)

	implicit none
	real randomtb
	integer idum
	integer m,ia,ic
	integer ir,iy,iff,j
	real rm

      parameter (m=714025,ia=1366,ic=150889,rm=1.4005112e-6)
      dimension ir(97)
      save iy,ir
      data iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        idum=mod(ic-idum,m)
        do j=1,97
          idum=mod(ia*idum+ic,m)
          ir(j)=idum
	end do
        idum=mod(ia*idum+ic,m)
        iy=idum
      endif
      j=1+(97*iy)/m
      if(j.gt.97.or.j.lt.1) stop 'error stop ran2: internal error'
      iy=ir(j)
      randomtb=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      return
      end

c****************************************************************

	subroutine test_ct

c testst conversion routine

	implicit none

	real hmax
c	parameter( hmax = 360. )
	parameter( hmax = 1. )

	integer i,j,k
	integer itot,iseed
	logical berror
	real r,g,b
	real rn,gn,bn
	real h,s,v
	real hn,sn,vn
	real eps

	real randomtb

	eps = 0.0001
	itot = 10000
	iseed = 124357
	berror = .false.

	write(6,*) '--------------'
c-----------------------------------------------

	do i=1,itot
	      
	  h = randomtb(iseed)
	  s = randomtb(iseed)
	  v = randomtb(iseed)

	  h = h * hmax
	  call hsv2rgb(h,s,v,r,g,b)
	  call rgb2hsv(r,g,b,hn,sn,vn)

	  write(6,1000) h,hn,s,sn,v,vn,r,g,b
	  call equaltb(h,hn,eps,berror)
	  call equaltb(s,sn,eps,berror)
	  call equaltb(v,vn,eps,berror)
	  if( berror ) goto 99

	end do
c-----------------------------------------------
	write(6,*) '--------------'
c-----------------------------------------------
	      
	do i=1,itot

	  r = randomtb(iseed)
	  g = randomtb(iseed)
	  b = randomtb(iseed)

	  call rgb2hsv(r,g,b,h,s,v)
	  call hsv2rgb(h,s,v,rn,gn,bn)

	  write(6,1000) r,rn,g,gn,b,bn,h,s,v
	  call equaltb(r,rn,eps,berror)
	  call equaltb(g,gn,eps,berror)
	  call equaltb(b,bn,eps,berror)
	  if( berror ) goto 99

	end do
c-----------------------------------------------
	write(6,*) '--------------'

	return
   99	continue
	write(6,*) r,g,b
	write(6,*) h,s,v
	stop 'error stop...'
 1000	format(9f8.3)
	end

c****************************************************************
	
	subroutine equaltb(a,b,eps,berror)

c tests for nearly equality

	implicit none

	real a,b
	real eps
	logical berror

	if( abs(a-b) .gt. eps ) then
	  write(6,*) '*** ',a,b,abs(a-b),eps
	  berror = .true.
	end if

	end

c****************************************************************

	subroutine test_ii

	implicit none

	real h,s,v
	real r,g,b
	logical h2r,r2h

	h2r = .false.
	h2r = .true.
	r2h = .true.
	r2h = .false.

	do while(h2r)
	  write(6,*) 'Enter h/s/v: '
	  read(5,*) h,s,v
	  call hsv2rgb(h,s,v,r,g,b)
	  write(6,*) 'r/g/b: ',r,g,b
	end do

	do while(r2h)
	  write(6,*) 'Enter r/g/b: '
	  read(5,*) r,g,b
	  call rgb2hsv(r,g,b,h,s,v)
	  write(6,*) 'h/s/v: ',h,s,v
	end do

	end

c****************************************************************
c
c	program test
c	call test_ct	!tests color transform
c	!call test_ii	!tests color transform interactively
c	end
c
c****************************************************************

