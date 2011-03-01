
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

	program xhcbs

	real xx(4),yy(4)

c tests HCBS graphic routines for fortran under X11

c getxwd returns size of window, that will be opened
c setxwd sets size of window to open (only befor qopen is called)
c (these units are in pixel, all subsequent units are real values)

	call getxwd(ixmin,iymin,ixmax,iymax)
	write(6,*) ixmin,iymin,ixmax,iymax
	call setxwd(ixmin,iymin,ixmax+200,iymax+100)

c qopen opens the window - has to be called once at the beginning

	call qopen

c qgetvp returns the dimension of the viewport (terminal) in viewport coord.
c (units of viewport coordinates are approx. in cm)

	call qgetvp(xvmin,yvmin,xvmax,yvmax)
	write(6,*) xvmin,yvmin,xvmax,yvmax

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c qstart - qend enclose one plot -- note that qend has to be called,
c ...otherwise the plot will not be flushed

	call qstart

c use qline, qmove, qplot or qpoint to plot

	call qline(2.,2.,7.,7.)
	call qplot(7.,3.)
	call qnewp(5)
	call qmove(4.,5.)
	call qplot(6.,7.)

	call qend
	call wait

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	call qstart

c qworld sets new world (window) coordinates - note that the corners
c ...of the viewport and window are coincident - if you want to avoid
c ...distortion, choose the dimensions of the world window and the
c ...viewport so that that the x/y is the same in both systems

c change color with qnewp [1-16]

	call qworld(100.,100.,400.,300.)
	call qnewp(8)
	call qmove(150.,250.)
	call qplot(350.,50.)

	call qend
	call wait

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	call qstart

c qviewp defines a window where clipping is enabled
c clipping happens by default at the outer edge of the viewport

	call qworld(100.,100.,300.,300.)
	call qviewp(3.,3.,13.,13.)
	call qnewp(13)
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
	call wait

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c this example shows how to fill an arbitrary shape (here triangle)

	call qstart

c we compute x in order to not have distorsion
c qviewp with all arguments 0. gives the whole window

	call qworld(0.,0.,100.*xvmax/yvmax,100.)
	call qviewp(0.,0.,0.,0.)
	
	xx(1)=50.
	yy(1)=50.
	xx(2)=50.
	yy(2)=90.
	xx(3)=60.
	yy(3)=60.
	xx(4)=90.
	yy(4)=50.

c qafill fills the shape defined by xx,yy with the actual foreground color
c the first argument gives the number of points in xx,yy

	call qnewp(4)
	call qafill(4,xx,yy)

	call qend
	call wait

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c this example gives all available colors
c qrfill fills a rectangle with the color passed

	call qstart

	call qworld(0.,0.,100.,100.)
	call qviewp(0.,0.,yvmax,yvmax)
	
c xm is size ideal, dm is plotting size
c in x,y the center of the ... are computed

	xm=100./8.
	dm=0.7*xm/2.

	icol=0
	do iy=1,2
	  y=(iy+3)*xm-xm/2.
	  y=100.-y
	  do ix=1,8
	    x=ix*xm-xm/2.
	    icol=icol+1
	    call qrfill(x-dm,y-dm,x+dm,y+dm,icol)
	  end do
	end do

	call qend
	call wait

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c call qclose at the end of the session

	call qclose

	end

	subroutine wait

	read(5,'(i10)') i

	end
