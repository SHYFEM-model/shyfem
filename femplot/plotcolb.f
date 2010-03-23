c
c $Id: plotcolb.f,v 1.5 2010-02-26 15:29:19 georg Exp $

	program plotcolb

c plots color bar to chose color

	implicit none

	real xmin,ymin,xmax,ymax
	real dy
	integer ncol,ntics,i
	character*2 line

	ncol = 100
	ntics = 10

	dy = 1.

	xmin = 2.
	xmax = xmin + 15.
	ymin = 24.

	call qopen
	call qstart

	call qtxts(20)
	call qfont('Times-Roman')
	call qtext(7.,25.,'Available Colortables')

	do i=0,6
	  ymin = ymin - 3.
	  call set_color_table(i)
	  call colb(xmin,ymin,xmax,ymin+dy,ncol,ntics)
	  write(line,'(i2)') i
	  call qtext(1.,ymin+0.3,line)
	end do

	call qend
	call qclose

	end

c*************************************************

	subroutine colb(xmin,ymin,xmax,ymax,ncol,ntics)

	implicit none

	real xmin,ymin,xmax,ymax
	integer ncol,ntics
	integer i
	real x0,col
	real db,dbt,dx
	character*4 line

	dx = xmax - xmin
	db = dx / ncol
	dbt = dx / ntics

	do i=1,ncol
	  col = i/float(ncol)
	  call qsetc(col)
	  x0 = xmin + (i-1) * db
	  call qrfill(x0,ymin,x0+db,ymax)
	end do

	call qgray(0.)
	do i=1,ntics
	  x0 = xmin + (i-1) * dbt
	  call pbox(x0,ymin,x0+dbt,ymax)
	end do

	call qtxts(12)
	call qfont('Times-Roman')

	do i=0,ntics
	  x0 = xmin + i * dbt
	  write(line,'(f4.2)') i/float(ntics)
	  call qtext(x0-0.5,ymin-0.5,line)
	end do

	end

c*************************************************

