
c*****************************************************************

	program lagxitest

c tests the internal coordinates and paths in one triangle

	implicit none

	integer i,iflux
	double precision r
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
	alpha = 0.4

	do i=1,15
	  call make_test_point(xip)
	  call random_number(alpha)
	  call random_number(r)
	  iflux = 3*r + 1
	  iflux = min(3,iflux)
	  call follow_path(iflux,xx,yy,xip,alpha)
	end do

	!call xi_test

	call qclose

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

	real x,y
	double precision xp,yp
	double precision beta,ap,gamma

	call xit2xy(xx,yy,xp,yp,xip)

	write(6,*) '========================= ',iflux
	write(6,*) xx
	write(6,*) yy
	write(6,*) xip
	write(6,*) xp,yp

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
	double precision beta,s

	call xit_start_end(iflux,alpha,xip,xis,xie)

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

	write(6,*) 'alpha = ',alpha
	write(6,*) xip

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

	subroutine xi_test

	double precision xx(3),yy(3)
	double precision xi(3),xip(3)
	double precision x,y
	double precision xis,xips,xidiff
	double precision eps,area,reg
	integer n,ii

	eps = 1.e-7

	do n=1,500000
	  do ii=1,3
	    call random_number(xx(ii))
	    call random_number(yy(ii))
	    call random_number(xi(ii))
	  end do

	  call xit_check(xx,yy,area,reg)
	  if( area <= 0. ) cycle
	  if( area <= 1.e-5 ) cycle
	  if( reg <= 1.e-4 ) cycle

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

	  !write(6,*) n,xis,xips
	  write(6,*) n,xidiff
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

