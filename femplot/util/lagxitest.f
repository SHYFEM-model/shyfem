
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
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

c*****************************************************************
c
c test internal routines
c
c*****************************************************************

! revision log :
!
! 18.12.2018	ggu	changed VERS_7_5_52
! 14.02.2019	ggu	changed VERS_7_5_56
! 21.05.2019	ggu	changed VERS_7_5_62

c*****************************************************************


	program xi_test_main

	use evgeom
	use basin

	integer iapini

        if(iapini(1,0,0,0).eq.0) then
                stop 'error stop : iapini'
        end if
	call set_ev
	write(6,*) nkn,nel

	call lagxi_path_plot_test
	call lagxi_path_test
	call lagxi_dist_test
	call lagxi_conv_test		!conversion routines

	end

c*****************************************************************

	subroutine lagxi_dist_test

c tests the internal coordinates and distances

	use evgeom
	use basin

	implicit none

	logical bdebug
	integer ie,i,ii,n,ifreq
	double precision xia(3),xib(3)
	double precision xa,ya,xb,yb
	double precision l,lxy,r,eps

	n = 100000
	ifreq = n / 10
	ifreq = max(ifreq,1)
	eps = 1.e-5
	bdebug = .false.

	write(6,*) 'running lagxi_dist_test: ',n

	do i=1,n
	  call random_number(r)
	  ie = r*nel + 1
	  ie = min(ie,nel)
	  call make_side_point(xia)
	  call make_side_point(xib)
	  call xi_dist(ie,xia,xib,l)
	  call xi2xy(ie,xa,ya,xia)
	  call xi2xy(ie,xb,yb,xib)
	  lxy = sqrt((xa-xb)**2+(ya-yb)**2)
	  if( mod(i,ifreq) == 0 ) write(6,*) i,ie,l,lxy
	  if( abs(l-lxy) > eps ) goto 99
	end do

	return
   99	continue
	  write(6,*) 'error in length: ',l,lxy
	  write(6,*) ie,nel
	  write(6,*) xia
	  write(6,*) xib
	  write(6,*) (ev(10+ii,ie),ii=1,3)
	  write(6,*) (ev(16+ii,ie),ii=1,3)
	stop 'error stop: distances'
	end

c*****************************************************************

	subroutine make_side_point(xi)

	implicit none

	double precision xi(3)

	integer i
	double precision r,eps

	eps = 0.05

	call random_number(r)
	i = 1 + 3.*r
	i = min(i,3)
	
	xi(i) = 0.

	i = mod(i,3) + 1
	call random_number(r)
	if( r < eps ) r = 0.
	xi(i) = r

	i = mod(i,3) + 1
	xi(i) = 1. - r

	end

c*****************************************************************

	subroutine lagxi_path_test

c tests the internal coordinates and paths in one triangle

	use evgeom
	use basin

	implicit none

	integer n,ie,ii,i,iflux,ifreq
	double precision r,s,diff,eps
	double precision xp,yp
	double precision alpha
	double precision xx(3),yy(3)
	double precision xip(3)
	double precision xis(3)
	double precision xie(3)
	double precision xi(3)

	n = 1000
	ifreq = n / 10
	ifreq = max(ifreq,1)
	eps = 1.e-5

	write(6,*) 'running lagxi_path_test: ',n

	do i=1,n
	  call make_triangle(ie,xx,yy)
	  call make_test_point(xip)
	  call random_number(alpha)
	  call random_number(r)
	  iflux = 3*r + 1
	  iflux = min(3,iflux)
	  call xit_start_end(iflux,alpha,xip,xis,xie,s)
	  diff = 0.
	  do ii=1,3
	    xi(ii) = (1.-s)*xis(ii) + s*xie(ii)
	    diff = diff + abs(xi(ii)-xip(ii))
	  end do
	  if( mod(i,ifreq) == 0 ) write(6,*) i,ie,s,diff
	  if( diff > eps ) goto 99
	end do

	return
   99	continue
	  write(6,*) 'error in parameter s: ',s,diff
	  write(6,*) ie,nel
	  write(6,*) xis
	  write(6,*) xie
	  write(6,*) xip
	  write(6,*) xi
	stop 'error stop: distances'
	end

c*****************************************************************

	subroutine lagxi_path_plot_test

c tests the internal coordinates and paths in one triangle and plot

	implicit none

	integer n,ie,i,iflux
	double precision r,s
	double precision xp,yp
	double precision alpha
	double precision xx(3),yy(3)
	double precision xip(3)

	call qopen

	xx(1) = 5.
	yy(1) = 5.
	xx(2) = 9.
	yy(2) = 6.
	xx(3) = 7.
	yy(3) = 9.

	n = 15
	write(6,*) 'running lagxi_path_plot_test: ',n

	do i=1,n
	  call make_triangle(ie,xx,yy)
	  call make_test_point(xip)
	  call random_number(alpha)
	  call random_number(r)
	  iflux = 3*r + 1
	  iflux = min(3,iflux)
	  !call xit_start_end(iflux,alpha,xip,xis,xie,s)
	  call follow_path(iflux,xx,yy,xip,alpha)
	end do

	call qclose

	end

c*****************************************************************

	subroutine make_triangle(ie,xx,yy)

	use evgeom
	use basin

	implicit none

	integer ie
	double precision xx(3),yy(3)
	
	integer ii,k
	double precision r

	call random_number(r)
	ie = r*nel + 1
	ie = min(ie,nel)

	do ii=1,3
	  k = nen3v(ii,ie)
	  xx(ii) = xgv(k)
	  yy(ii) = ygv(k)
	end do

	end

c*****************************************************************

	subroutine make_test_point(xip)

	implicit none

	double precision xip(3)
	double precision a,b

	call random_number(a)
	call random_number(b)

	xip(1) = a
	xip(2) = (1.-a)*b
	xip(3) = 1. - xip(1) - xip(2)
	
	end

c*****************************************************************

	subroutine follow_path(iflux,xx,yy,xip,alpha)

	implicit none

	integer iflux
	double precision xx(3),yy(3)
	double precision xip(3)
	double precision alpha

	logical bdebug
	real x,y
	double precision xp,yp
	double precision beta,ap,gamma

	bdebug = .false.

	call xit2xy(xx,yy,xp,yp,xip)

	write(6,*) '== ',iflux,alpha
	if( bdebug ) then
	  write(6,*) xx
	  write(6,*) yy
	  write(6,*) xip
	  write(6,*) xp,yp
	end if

	call qstart

	call plot_triang(xx,yy)
	call plot_flux(iflux,xx,yy,alpha)
	x = xp
	y = yp
	call qpsize(0.1)
	call qpoint(x,y)
	call plot_line(iflux,xx,yy,alpha,xip)

	call qend

	end

c*****************************************************************

	subroutine plot_line(iflux,xx,yy,alpha,xip)

	implicit none

	integer iflux
	double precision xx(3),yy(3)
	double precision alpha
	double precision xip(3)

	double precision xp,yp
	double precision as,ae
	double precision xi(3)
	double precision xis(3)
	double precision xie(3)

	integer i,n,i1
	real x,y
	double precision s

	call xit_start_end(iflux,alpha,xip,xis,xie,s)

	call xit2xy(xx,yy,xp,yp,xis)
	x = xp
	y = yp
	call qmove(x,y)

	call xit2xy(xx,yy,xp,yp,xie)
	x = xp
	y = yp
	call qplot(x,y)

	end

c*****************************************************************

	subroutine plot_flux(iflux,xx,yy,alpha)

	implicit none

	integer iflux
	double precision xx(3),yy(3)
	double precision alpha

	double precision xip(3)
	double precision xd,yd

	integer i1,i2,i3
	real x,y

	i1 = iflux
        i2 = mod(i1,3) + 1
        i3 = mod(i2,3) + 1

	xip(i1) = 0.
	xip(i2) = 1.-alpha
	xip(i3) = alpha

	!write(6,*) 'alpha = ',alpha
	!write(6,*) xip

	x = xx(i1)
	y = yy(i1)
	call qmove(x,y)

	call xit2xy(xx,yy,xd,yd,xip)
	x = xd
	y = yd
	call qplot(x,y)

	end

c*****************************************************************

	subroutine plot_triang(xx,yy)

	implicit none

	double precision xx(3),yy(3)

	integer ii
	real xmin,xmax,ymin,ymax
	real x,y,fact,dx,dy

	xmin = minval(xx)
	xmax = maxval(xx)
	ymin = minval(yy)
	ymax = maxval(yy)

	fact = 0.1
	dx = fact*(xmax - xmin)
	dy = fact*(ymax - ymin)

	xmin = xmin - dx
	xmax = xmax + dx
	ymin = ymin - dy
	ymax = ymax + dy

	call qworld(xmin,ymin,xmax,ymax)

	x=xx(3)
	y=yy(3)
	call qmove(x,y)
	do ii=1,3
	  x=xx(ii)
	  y=yy(ii)
	  call qplot(x,y)
	end do

	end

c*****************************************************************

	subroutine lagxi_conv_test

c tests conversion routines from and to internal coordinates

	double precision xx(3),yy(3)
	double precision xi(3),xip(3)
	double precision x,y
	double precision xis,xips,xidiff
	double precision eps,area,reg
	integer i,n,ii

	n = 500000
	ifreq = n / 10
	ifreq = max(ifreq,1)
	eps = 1.e-7
	i = 0

	write(6,*) 'running lagxi_conv_test: ',n

	do while( i < n )
	  do ii=1,3
	    call random_number(xx(ii))
	    call random_number(yy(ii))
	    call random_number(xi(ii))
	  end do

	  call xit_check(xx,yy,area,reg)
	  if( area <= 0. ) cycle
	  if( area <= 1.e-4 ) cycle
	  if( reg <= 1.e-3 ) cycle

	  i = i + 1
	  !xi(1) = 1.
	  xi(2) = xi(2)*(1.-xi(1))
	  xi(3) = 1. - xi(1) - xi(2)

	  call xit2xy(xx,yy,x,y,xi)
	  call xy2xit(xx,yy,x,y,xip)

	  xidiff = 0.
	  do ii=1,3
	    xidiff = xidiff + abs(xi(ii)-xip(ii))
	  end do

	  xis = xi(1) + xi(2) + xi(3)
	  xips = xip(1) + xip(2) + xip(3)

	  if( mod(i,ifreq) == 0 ) write(6,*) i,xidiff
	  if( xidiff .gt. eps ) then
	    write(6,*) x,y
	    write(6,*) xx
	    write(6,*) yy
	    write(6,*) xi
	    write(6,*) xip
	    call xit_check(xx,yy,area,reg)
	    write(6,*) area,reg
	    stop 'error stop xi_test'
	  end if
	end do

	end

c*****************************************************************

