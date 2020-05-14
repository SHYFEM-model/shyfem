
!--------------------------------------------------------------------------
!
!    Copyright (C) 1992,1994,2018  Georg Umgiesser
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

!--------------------------------------------------------------------------
!
! xhcbs.f - test for HCBS simulation routines for FORTRAN under X11
!
! revision log :
!
! 01.01.1992	ggu	routines written from scratch
! 11.02.1994	ggu	copyright notice added to all files
! 29.10.2018	ggu	new copyright added
!
!--------------------------------------------------------------------------

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
