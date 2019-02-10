
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

c************************************************************************\ 
c*									*
c* colorscale.f		plots various color scales			*
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

	program colorscale

	character*50 text

	call qopen

	call qgetvp(xvmin,yvmin,xvmax,yvmax)
	write(6,*) xvmin,yvmin,xvmax,yvmax

	call qstart

	xymin = 0.
	xymax = 100.

	call qworld(xymin,xymin,xymax,xymax)
	call qsetvp(0.,0.,xvmax,xvmax)
	
	call qmove(xymin,xymin)
	call qplot(xymin,xymax)
	call qplot(xymax,xymax)
	call qplot(xymax,xymin)
	call qplot(xymin,xymin)

c-------------------------------------------------------------------
c plot color bars
c-------------------------------------------------------------------

c	----------------------------------------------------------
c	initialize
c	----------------------------------------------------------

	ncbars = 5		!number of color bars
	n = 100			!number of colors
	xmarg = 10.		!margin in x

	dy = (xymax - xymin) / (2*ncbars)	!height of space for color bar
	dx = (xymax-xymin-2.*xmarg)/n
	ddx = dx/2.
	ddy = dy/2.
	dty = 3.*dy/4.

	hmin = 0.0		!min hue
	hmax = 1.0		!max hue

	x = xmarg		!start from here in x
	y = xymin - dy		!start in y

c	----------------------------------------------------------
c	gray scale
c	----------------------------------------------------------

	y = y + 2*dy
	text = 'gray scale'
	call gray_bar(n)
	call color_bar(n,x,y,dx,dy,text)

c	----------------------------------------------------------
c	hue scale
c	----------------------------------------------------------

	y = y + 2*dy
	text = 'hue scale'
	call hue_bar(n)
	call color_bar(n,x,y,dx,dy,text)

c	----------------------------------------------------------
c	hue**0.5 scale
c	----------------------------------------------------------

	y = y + 2*dy
	text = 'hue**0.5 scale'
	call hue_squareroot_bar(n)
	call color_bar(n,x,y,dx,dy,text)

c	----------------------------------------------------------
c	red white blue scale
c	----------------------------------------------------------

	y = y + 2*dy
	text = 'red white blue scale'
	call red_white_blue_bar(n)
	call color_bar(n,x,y,dx,dy,text)

c	----------------------------------------------------------
c	red black blue scale
c	----------------------------------------------------------

	y = y + 2*dy
	text = 'red black blue scale'
	call red_black_blue_bar(n)
	call color_bar(n,x,y,dx,dy,text)

c-------------------------------------------------------------------
c end of plotting
c-------------------------------------------------------------------

	call qend

	call qclose

c-------------------------------------------------------------------
c end of routines
c-------------------------------------------------------------------

	end

c********************************************************

	subroutine color_bar(n,x0,y0,dx,dy,text)

c color bar

	implicit none

	integer n
	real x0,y0,dx,dy
	character*(*) text

	integer i
	real ddx,ddy,dty
	real x,y

	ddx = dx/2.
	ddy = dy/2.
	dty = 3.*dy/4.

	x = x0
	y = y0

	write(6,*) 'plotting ',text
	call qtext(x,y+dty,text)
	call make_color_bar(n,x,y,dx,dy)
	call xtick(x,y-ddy,x+n*dx,y+ddy,dy/10.,10)

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

	subroutine make_color_bar(n,x0,y0,dx,dy)

c makes color bar

	implicit none

	integer n
	real x0,y0,dx,dy

	integer i
	real ddx,ddy
	real x,y

	ddx = dx/2.
	ddy = dy/2.
	y = y0

	do i=1,n
	    x = x0 + i*dx - ddx
	    call qsetcpen(i)
	    call qrfill(x-ddx,y-ddy,x+ddx,y+ddy)
	end do

	end

c********************************************************

	subroutine intplin(i0,i1,i,r0,r1,r)

c linear interpolation

	implicit none

	integer i0,i1,i
	real r0,r1,r

	r = r0 + (r1-r0) * float(i-i0) / float(i1-i0)

	end

c********************************************************
c********************************************************
c********************************************************

	subroutine gray_bar(n)

c makes gray color scale 

	implicit none

	integer n
	integer i
	real g

	call qinitct(n)

	g=0.
	do i=1,n
	  call intplin(1,n,i,0.,1.,g)
	  call qsetct(i,1,g,g,g)
	end do

	end

c********************************************************

	subroutine hue_bar(n)

c makes hue color scale 

	implicit none

	integer n
	integer i
	real h

	call qinitct(n)

	do i=1,n
	  call intplin(1,n,i,0.,1.,h)
	  call qsetct(i,0,h,1.,1.)
	end do

	end

c********************************************************

	subroutine hue_squareroot_bar(n)

c makes hue color scale 

	implicit none

	integer n
	integer i
	real h,hmin,hmax

	hmin = 0.
	hmax = 0.666

	call qinitct(n)

	do i=1,n
	  call intplin(1,n,i,hmin,hmax,h)
	  h = h ** 0.5
	  call qsetct(i,0,h,1.,1.)
	end do

	end

c********************************************************

	subroutine red_black_blue_bar(n)

c makes color scale from red to black to blue

	implicit none

	integer n

	integer n2,i
	real r,g,b

	n2 = n / 2

	call qinitct(n)

	r=0.
	g=0.
	b=0.
	do i=1,n2
	  call intplin(1,n2,i,1.,0.,r)
	  call qsetct(i,1,r,g,b)
	end do

	r=0.
	g=0.
	b=0.
	do i=n2+1,n
	  call intplin(n2+1,n,i,0.,1.,b)
	  call qsetct(i,1,r,g,b)
	end do

	end

c********************************************************

	subroutine red_white_blue_bar(n)

c makes color scale from red to white to blue

	implicit none

	integer n

	integer n2,i
	real r,g,b,c
	real cmin,cmax

	n2 = n / 2
	cmin = 0.3
	cmax = 1.0

	call qinitct(n)

	r=1.
	do i=1,n2
	  call intplin(1,n2,i,cmin,cmax,c)
	  c = c*c
	  call qsetct(i,1,r,c,c)
	end do

	b=1.
	do i=n2+1,n
	  call intplin(n2+1,n,i,cmax,cmin,c)
	  c = c*c
	  call qsetct(i,1,c,c,b)
	end do

	end

c********************************************************


