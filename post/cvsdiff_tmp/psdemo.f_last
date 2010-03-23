
c************************************************************************\ 
c*									*
c* xhcbs.f - test for HCBS simulation routines for FORTRAN under X11	*
c*									*
c* Copyright (c) 1992-1994 by Georg Umgiesser				*
c*									*
c* Permission to use, copy, modify, and distribute this software	*
c* and its documentation for any purpose and without fee is hereby	*
c* granted, provided that the above copyright notice appear in all	*
c* copies and that both that copyright notice and this permission	*
c* notice appear in supporting documentation.				*
c*									*
c* This file is provided AS IS with no warranties of any kind.		*
c* The author shall have no liability with respect to the		*
c* infringement of copyrights, trade secrets or any patents by		*
c* this file or any part thereof.  In no event will the author		*
c* be liable for any lost revenue or profits or other special,		*
c* indirect and consequential damages.					*
c*									*
c* Comments and additions should be sent to the author:			*
c*									*
c*			Georg Umgiesser					*
c*			ISDGM/CNR					*
c*			S. Polo 1364					*
c*			30125 Venezia					*
c*			Italy						*
c*									*
c*			Tel.   : ++39-41-5216875			*
c*			Fax    : ++39-41-2602340			*
c*			E-Mail : georg@lagoon.isdgm.ve.cnr.it		*
c*									*
c* Revision History:							*
c* 11-Feb-94: copyright notice added to all files			*
c* ..-...-92: routines written from scratch				*
c*									*
c************************************************************************/

	program psdemo

	real xx(4),yy(4)

c data for filling routine

	data xx /50.,50.,60.,90./
	data yy /50.,90.,60.,50./

c tests HCBS graphic routines for fortran under X11

c getxwd returns size of window, that will be opened
c setxwd sets size of window to open (only befor qopen is called)
c (these units are in pixel, all subsequent units are real values)

c qopen opens the window - has to be called once at the beginning

	call qopen

c qgetvp returns the dimension of the viewport (terminal) in viewport coord.
c (units of viewport coordinates are approx. in cm)

	call qgetvp(xvmin,yvmin,xvmax,yvmax)
	write(6,*) xvmin,yvmin,xvmax,yvmax

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c this page shows the fundamental plotting commands

c qstart - qend enclose one plot -- note that qend has to be called,
c ...otherwise the plot will not be flushed

	call qstart

c use qline, qmove, qplot or qpoint for plotting
c change color(grayscale) with qgray [0-1]

	call qline(2.,2.,7.,7.)
	call qplot(7.,3.)
	call qgray(0.5)
	call qmove(4.,5.)
	call qplot(6.,7.)

c qworld sets new world (window) coordinates - note that the corners
c ...of the viewport and window are coincident - if you want to avoid
c ...distortion, choose the dimensions of the world window and the
c ...viewport so that that the x/y is the same in both systems
c ...or you can also call subroutine qrcfy, that adjusts the scale

	call qworld(100.,100.,400.,300.)
	call qgray(0.8)
	call qmove(150.,250.)
	call qplot(350.,50.)

	call qend

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c this page shows how clipping is working

	call qstart

c qviewp defines a window where clipping is enabled
c clipping happens by default at the outer edge of the viewport

	call qworld(100.,100.,300.,300.)
	call qsetvp(3.,3.,13.,13.)
	call qgray(0.3)
	call qmove(200.,200.)
	dr=10./36.
	dalf=10.*3.14159/180.
	r=0.
	alf=0.
	do i=1,900
	  x=200.+r*cos(alf)
	  y=200.+r*sin(alf)
	  call qplot(x,y)
	  r=r+dr
	  alf=alf+dalf
	end do

	call qend

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c this example shows how to fill an arbitrary shape

	call qstart

c we compute x in order to not have distorsion
c ...or we can simply call qrcfy

	call qworld(0.,0.,100.,100.)
	call qrcfy
	
c qafill fills the shape defined by xx,yy with the actual foreground color
c the first argument gives the number of points in xx,yy
c the vectors xx,yy are set in data statement above

	call qgray(0.4)
	call qafill(4,xx,yy)

c now translate the triangle and fill with arbitrary color

	do i=1,4
	  xx(i) = xx(i) - 25.
	  yy(i) = yy(i) - 25.
	end do

	call qrgb(0.3,0.8,0.5)
	call qafill(4,xx,yy)

	call qend

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c here we plot colorbars - one in shades of gray, and one in color
c ... for the color we use the call to qhue()

	call qstart

	xymin = 0.
	xymax = 100.

	call qworld(xymin,xymin,xymax,xymax)
	call qsetvp(0.,0.,xvmax,xvmax)
	
	ncbars = 5		!number of color bars
	n = 100			!number of colors
	xmarg = 10.		!margin in x
	hmin = 0.0		!min hue
	hmax = 1.0		!max hue
	dy = (xymax - xymin) / (2*ncbars)	!height of space for color bar
	dx = (xymax-xymin-2.*xmarg)/n
	ddx = dx/2.
	ddy = dy/2.

	y = xymin - ddy		!start in y

c	first we do the gray scale ---------------------------------------

	xm = xmarg		!start from here in x
	y = y + 2*dy

	do ix=1,n
	    x = xm + ix*dx - dx/2.
	    hue = hmin + ix * (hmax - hmin ) / n
	    call qgray(hue)
	    call qrfill(x-ddx,y-ddy,x+ddx,y+ddy)
	end do

	call xtick(xm,y-ddy,xm+n*dx,y+ddy,dy/10.,10)

c	now we do the color ----------------------------------------------

	xm = xmarg		!start from here in x
	y = y + 2*dy

	do ix=1,n
	    x = xm + ix*dx - dx/2.
	    hue = hmin + ix * (hmax - hmin ) / n
	    call qhue(hue)
	    call qrfill(x-ddx,y-ddy,x+ddx,y+ddy)
	end do

	call xtick(xm,y-ddy,xm+n*dx,y+ddy,dy/10.,10)

c other color composition ---------------------------------------------------

	n3 = n/3
	y = y + 2*dy

	do ix=1,n
	    x = xm + ix*dx - dx/2.
	    fx = ix
	    if( ix .le. n3 ) then
		fx = ix
		red = 1. - fx/n3
		green = fx/n3
		blue = 0.
		call wheel(fx,n3,red,green)
	    else if( ix .gt. 2*n3 ) then
		fx = ix - 2*n3
		red = fx/n3
		green = 0.
		blue = 1. - fx/n3
		call wheel(fx,n3,blue,red)
	    else
		fx = ix - n3
		red = 0.
		green = 1. - fx/n3
		blue = fx/n3
		call wheel(fx,n3,green,blue)
	    end if
	    write(6,'(2i5,4f10.4)') ix,n3,fx,red,green,blue
	    call qrgb(red,green,blue)
	    call qrfill(x-ddx,y-ddy,x+ddx,y+ddy)
	end do

	call xtick(xm,y-ddy,xm+n*dx,y+ddy,dy/10.,10)

c	hsb ^ 0.5 -------------------------------------------------------

	xm = xmarg		!start from here in x
	y = y + 2*dy
	hmin = 0.
	hmax = 0.666

	do ix=1,n
	    x = xm + ix*dx - dx/2.
	    hue = hmin + ix * (hmax - hmin ) / n
	    hue = hue ** 0.5
	    call qhue(hue)
	    call qrfill(x-ddx,y-ddy,x+ddx,y+ddy)
	end do

	call xtick(xm,y-ddy,xm+n*dx,y+ddy,dy/10.,10)

c	special Maxwell ----------------------------------------------

	xm = xmarg		!start from here in x
	y = y + 2*dy
	n2 = n/2

	do ix=1,n2
	    x = xm + ix*dx - dx/2.
	    t = ix/float(n2)
	    call rglin(0.,0.,1./3.,1./3.,t,r,g)
	    call max2hsv(r,g,hue,sat,bri)
	hue = 2./3.
	sat = 1 - ix/float(n2)
	write(66,*) ix,r,g,hue,sat,bri
	    call qhsb(hue,sat,bri)
	    call qrfill(x-ddx,y-ddy,x+ddx,y+ddy)
	end do

	do ix=n2+1,n
	    x = xm + ix*dx - dx/2.
	    t = (ix-n2)/float(n2)
	    call rglin(1./3.,1./3.,1.,0.,t,r,g)
	    call max2hsv(r,g,hue,sat,bri)
	write(66,*) ix,r,g,hue,sat,bri
	    call qhsb(hue,sat,bri)
	    call qrfill(x-ddx,y-ddy,x+ddx,y+ddy)
	end do

	call xtick(xm,y-ddy,xm+n*dx,y+ddy,dy/10.,10)

c	done -----------------------------------------------------------

c	call qrgb(0.,1.,1.)
c	x = xm
c	y = xymax*4./5.
c	call qrfill(x-ddx,y-ddy,x+ddx,y+ddy)

	call qend

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c some test for text

        call qstart

        call qtext(3.,3.,'This is a test in Courier, 12 pt')

        call qtxts(20)
        call qfont('Times-Roman')
        call qtext(3.,5.,'This is Times-Roman, 20 pt')

c        call qtxts(12)
c        call qtext(1.,5.,'Just a test')
c        call qtext(999.,999.,' and more')
c        call qtext(999.,6.,' and still more')
c        call qtext(12.,999.,' and again more')

        call qtxts(10)
        call qfont('Times-Roman')
        call qtext(1.,1.,'Just a test')
        call qtext(999.,999.,' continue after the old string')
        call qtext(999.,2.,' continue in x, but one cm higher')
        call qtext(15.,999.,' same y, but at x=15')

        call qend

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c call qclose at the end of the session

	call qclose

	end

c********************************************************

	subroutine wheel(fx,n,col1,col2)

c makes color wheel (or part of it)

	implicit none

	real fx
	real col1,col2
	integer n

	real pi2,rad

	pi2 = 3.14159 / 2.
	rad = pi2 * fx / n

	col2 = sin(rad)
	col1 = cos(rad)

	end

c********************************************************

	subroutine xtick(x1,y1,x2,y2,dy,n)

c writes ticks on x axis

	implicit none

	real x1,y1,x2,y2
	real dy
	integer n

	integer i
	real dx,x

	call qgray(0.)
	call qmove(x1,y1)
	call qplot(x2,y1)
	call qplot(x2,y2)
	call qplot(x1,y2)
	call qplot(x1,y1)

	dx = (x2-x1) / n

	do i=0,n
	  x = x1 + i * dx
	  call qline(x,y1,x,y1-dy)
	end do

	end

c********************************************************

	subroutine rglin(r1,g1,r2,g2,t,r,g)

c interpolates linearly in Maxwell triangle

	implicit none

	real r1,g1,r2,g2
	real t			!t in [0-1]
	real r,g

	r = r1 + t*(r2-r1)
	g = g1 + t*(g2-g1)

	end

c********************************************************

	subroutine max2hsv(r,g,h,s,v)

c Maxwell values to HSV color space

	implicit none

	real r,g
	real h,s,v

	real b,p
	real cmax,cmin

	b = 1. - r - g

	cmax = max(r,g,b)
	cmin = min(r,g,b)

	if( cmax .gt. 1. .or. cmin .lt. 0. ) then
	  write(6,*) r,g,b
	  write(6,*) cmin,cmax
	  stop 'error stop max2hsv: internal error (1)'
	end if

	if( r .eq. cmin ) then		!between g and b
	  p = g/(1.-cmin)
	  !q = b/(1.-cmin)
	  h = 1./3. + (1.-p)/3.
	else if( g .eq. cmin ) then		!between b and r
	  p = b/(1.-cmin)
	  !q = r/(1.-cmin)
	  h = 2./3. + (1.-p)/3.
	else if( b .eq. cmin ) then		!between r and g
	  p = r/(1.-cmin)
	  !q = g/(1.-cmin)
	  h = (1.-p)/3.
	else
	  stop 'error stop max2hsv: internal error (2)'
	end if

	s = 1. - cmin/cmax
	v = 1.

	end

c********************************************************






