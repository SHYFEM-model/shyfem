
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

!--------------------------------------------------------------------------
!
! colorbar.f		plots various color scales
!
! revision log :
!
! 01.01.1992	ggu	routines written from scratch
! 11.02.1994	ggu	copyright notice added to all files
! 20.01.2009	ggu	new color routines used
! 29.10.2018	ggu	new copyright added
!
!--------------------------------------------------------------------------

	program colorbar

	implicit none

	integer ndim
	parameter (ndim=7)

	character*50 texts(0:ndim)

	data texts	/
     +			 'gray scale'
     +			,'hue scale'
     +			,'white blue scale'
     +			,'white red scale'
     +			,'blue white red scale'
     +			,'blue black red scale'
     +			,'hue**0.5 scale'
     +			,'hue**2 scale'
     +			/

	integer ncol,ncbars,i
	real xvmin,yvmin,xvmax,yvmax
	real xymin,xymax
	real xmarg,dy,dx,ddx,ddy,dty,x,y
	character*50 text

c-------------------------------------------------------------------
c preliminaries
c-------------------------------------------------------------------

	call qopen

	call qgetvp(xvmin,yvmin,xvmax,yvmax)
	write(6,*) xvmin,yvmin,xvmax,yvmax

	call qstart

	xymin = 0.
	xymax = 100.

	call qworld(xymin,xymin,xymax,xymax)
	!call qsetvp(0.,0.,xvmax,xvmax)
	call qsetvp(0.,0.,xvmax,yvmax)
	
	call make_rect(xymin,xymin,xymax,xymax)

c-------------------------------------------------------------------
c plot color bars
c-------------------------------------------------------------------

c	----------------------------------------------------------
c	initialize
c	----------------------------------------------------------

	ncol = 100			!number of colors

	ncbars = ndim+1		!number of color bars
	xmarg = xymax/10.	!margin in x

	dy = (xymax - xymin) / (3*ncbars)	!height of space for color bar
	dx = (xymax-xymin-2.*xmarg)/ncol
	ddx = dx/2.
	ddy = dy/2.
	dty = 3.*dy/4.

	x = xmarg		!start from here in x
	y = xymin - 1.5*dy		!start in y

c-------------------------------------------------------------------
c loop over color bars
c-------------------------------------------------------------------

	do i=0,ndim
	  y = y + 3*dy
	  text = texts(i)
	  call color_bar(i,ncol,x,y,dx,dy,text)
	end do

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

	subroutine color_bar(icol,n,x0,y0,dx,dy,text)

c color bar

	implicit none

	integer icol
	integer n
	real x0,y0,dx,dy
	character*(*) text

	real ddx,ddy,dty
	real x,y
	character*60 line

	ddx = dx/2.
	ddy = dy/2.
	dty = 3.*dy/4.

	x = x0
	y = y0

	write(line,'(i2)') icol
	line = line(1:2) // " " // text

	write(6,*) 'plotting ',line
	call qtext(x,y+dty,line)
	call make_color_bar(icol,n,x,y,dx,dy)
	call xtick(x,y-ddy,x+n*dx,y+ddy,dy/10.,10)

	end

c********************************************************

	subroutine xtick(x1,y1,x2,y2,dy,n)

c writes ticks on x axis

	implicit none

	real x1,y1,x2,y2
	real dy
	real w,h
	integer n

	integer i
	real dx,x
	character*3 line

	call qgray(0.)
	call make_rect(x1,y1,x2,y2)
	call qtxtcc(0,1)

	dx = (x2-x1) / n

	do i=0,n
	  x = x1 + i * dx
	  call qline(x,y1,x,y1-dy)
	  write(line,'(f3.1)') 0.1*i
	  call qtsize(line,w,h)
	  call qtext(x,y1-dy,line)
	  !call qtext(x-0.5*w,y1-2*h,line)
	  !call make_rect(x-0.5*w,y1-2*h,x+0.5*w,y1-h)
	end do

	call qtxtcc(-1,-1)

	end

c********************************************************

	subroutine make_color_bar(icol,n,x0,y0,dx,dy)

c makes color bar

	implicit none

	integer icol
	integer n
	real x0,y0,dx,dy

	integer i
	real color
	real ddx,ddy
	real x,y

	ddx = dx/2.
	ddy = dy/2.
	y = y0

	do i=1,n
	    color = float(i) / n
	    x = x0 + i*dx - ddx
            call qsetc( icol , color )
	    call qrfill(x-ddx,y-ddy,x+ddx,y+ddy)
	end do

	end

c********************************************************

	subroutine make_rect(x1,y1,x2,y2)

	implicit none

	real x1,y1,x2,y2

	call qmove(x1,y1)
	call qplot(x2,y1)
	call qplot(x2,y2)
	call qplot(x1,y2)
	call qplot(x1,y1)

	end

c********************************************************
c********************************************************
c********************************************************

        subroutine qsetc( icol , color )

c changes color using the actual color table

        implicit none

	integer icol	!color table
        real color	!color

        if( icol .eq. 0 ) then
          call qgray(color)
        else if( icol .eq. 1 ) then
          call qhue(color)
        else if( icol .eq. 2 ) then
          call white_blue( color )
        else if( icol .eq. 3 ) then
          call white_red( color )
        else if( icol .eq. 4 ) then
          call blue_white_red( color )
        else if( icol .eq. 5 ) then
          call blue_black_red( color )
        else if( icol .eq. 6 ) then
          call hue_sqrt( color )
        else if( icol .eq. 7 ) then
          call hue_pow2( color )
        else
          write(6,*) 'icol = ',icol
          stop 'error stop qsetc: no such color table'
        end if

        end

c********************************************************
c********************************************************
c********************************************************

        subroutine white_blue( color )
	implicit none
        real color
        call qhsb(0.666,color,1.)
        end

        subroutine white_red( color )
	implicit none
        real color
        call qhsb(1.,color,1.)
        end

        subroutine blue_white_red( color )
	implicit none
        real color
        if( color .le. 0.5 ) then
          call qhsb(0.666,1.-2.*color,1.)
        else
          call qhsb(1.,2.*(color-0.5),1.)
        end if
        end

        subroutine blue_black_red( color )
	implicit none
        real color
        if( color .le. 0.5 ) then
          call qhsb(0.666,1.,1.-2.*color)
        else
          call qhsb(1.,1.,2.*(color-0.5))
        end if
        end

        subroutine hue_sqrt( color )
	implicit none
        real color
        call qhue(sqrt(color))
        end

        subroutine hue_pow2( color )
	implicit none
        real color
        call qhue(color*color)
        end

c********************************************************


