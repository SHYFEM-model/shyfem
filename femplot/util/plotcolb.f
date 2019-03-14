
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2019  Georg Umgiesser
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

	program plotcolb

c plots color bar to chose color

	implicit none

	real xmin,ymin,xmax,ymax
	real dy,dty,xtmin
	logical berr
	integer ncol,ntics,i,ii,nmax
	integer imap
	character*80 line
	character*80 cname

	ncol = 100
	ntics = 10
	nmax = 9	!max number of color tables per page

	dy = 1.
	dty = 0.3

	xtmin = 0.2
	xmin = 4.
	xmax = xmin + 15.

	cname = ' '
	imap = 0
	call color_table_file_init(' ')
	call read_color_table(cname,imap,berr)	!this shows available CTs
	if( berr ) stop

	call qopen
	call qstart

	call qtxts(20)
	call qfont('Times-Roman')
	call qtext(7.,25.,'Available Colortables')

	ymin = 24.
	do i=0,7
	  ymin = ymin - 2.
	  call set_color_table(i)
	  call colb(xmin,ymin,xmax,ymin+dy,ncol,ntics)
	  write(line,'(i2)') i
	  call qtext(1.,ymin+dty,line)
	end do

	call qend

	call set_color_table(8)	!this indicates to used named color tables

	do ii=1,3
	  ymin = 26.
	  call qstart
	  do i=1,nmax
	    imap = nmax*(ii-1) + i
	    cname = ' '
	    call read_color_table(cname,imap,berr)
	    if( berr ) exit
	    write(6,*) 'plotting color table: ',imap,trim(cname)
	    ymin = ymin - 2.
	    call colb(xmin,ymin,xmax,ymin+dy,ncol,ntics)
	    write(line,'(i2,a,a)') imap,'  ',trim(cname)
	    call qtext(xtmin,ymin+dty,trim(line))
	  end do
	  call qend
	end do

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

