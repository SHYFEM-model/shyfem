
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

!*******************************************************************
! examples to write formatted fem files
!*******************************************************************

	program fem_test

	call profile		!one point profile
	call regular_zeta	!regular zeta (water level)
	call regular_wind	!regular wind
	call regular_conz_2d	!regular 2d grid
	call regular_conz_3d	!regular 3d grid

	end

!*******************************************************************

	subroutine profile

! shows how to write a profile for a fem file
! only one point is written
!
! more information can be found in file subfemfile.f

	implicit none

	integer, parameter :: lmax = 5

	real temp1(lmax)
	real temp2(lmax)
	real hlv(lmax)
	real regpar(7)		!not used
	integer ilhkv(1)	!just one point
	integer datetime(2)

	integer nvers,np,nvar,ntype,nlvddi,iformat
	integer iunit
	double precision dtime,dtime1,dtime2
	real hd
	character*20 string

	nvers = 0	!version of femfile, 0 for latest
	np = 1		!just one horizontal point
	nvar = 1	!just one variable
	ntype = 1	!we want with date
	nlvddi = lmax	!vertical dimension
	ilhkv(1) = lmax	!number of layers per point (only one point)
	hd = 15.	!total depth, not important for this application
	iformat = 1	!file should be formatted
	regpar = 0.	!no regular grid for values
	string = 'temperature'

! in hlv is the layer structure, values indicate the bottom of the layer
! in datetime is the reference date
! in temp1/2 are the temperature values for each layer
! the temperature values refer to the center of the layer

	hlv = (/2.,4.,7.,10.,15./)	!depth structure
	datetime = (/20150101,0/)	!reference date is 1.1.2015
	dtime1 = 0.			!time of first time stamp
	dtime2 = 86400.			!time of second time stamp
	temp1 = (/18.,17.,16.,14.,12./)	!temp of first time stamp
	temp2 = (/19.,17.,15.,13.,11./)	!temp of second time stamp

	iunit = 1
	open(iunit,file='profile.fem',form='formatted',status='unknown')

	dtime = dtime1
        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)
        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvddi,temp1)

	dtime = dtime2
        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)
        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvddi,temp2)

	close(iunit)

	end

!*******************************************************************

	subroutine regular_zeta

! shows how to write a regular water level fem file - use as initial condition
!
! more information can be found in file subfemfile.f

	implicit none

	integer, parameter :: lmax = 1
	integer, parameter :: nx = 7		!points in x direction
	integer, parameter :: ny = 8		!points in y direction

	real zeta(nx,ny)
	real hlv(lmax)		!not needed - water level is 2d
	real regpar(7)
	integer ilhkv(1)	!not needed, only one level
	integer datetime(2)

	integer nvers,np,nvar,ntype,nlvddi,iformat
	integer iunit
	integer ix,iy,ind
	double precision dtime
	real dxy,dx,dy,x0,y0,x1,y1
	real hd,flag
	character*20 string

	dxy = 0.1
	x0 = 20.4
	y0 = 54.7
	x1 = x0 + (nx-1)*dxy
	y1 = y0 + (ny-1)*dxy
	dx = dxy
	dy = dxy

	nvers = 0	!version of femfile, 0 for latest
	np = nx*ny	!total number of points
	ntype = 11	!we want with date and regular information
	nvar = 1	!zeta file has only water level as variable
	nlvddi = lmax	!vertical dimension, must be 1 for 2D
	ilhkv = lmax	!number of layers per point, not needed
	hd = 1.		!total depth, not needed
	iformat = 1	!file should be formatted
	flag = -999.
	string = 'water level'

! set array with information on regular grid
! points run from x0 to x0+(nx-1)*dx and y0 to y0+(ny-1)*dy
! values are stored row-wise from left to right starting from the lowest row
! first point is [x0,y0], second [x0+dx,y0], etc..

	regpar(1) = nx	!number of points in x direction (nx)
	regpar(2) = ny	!number of points in y direction (ny)
	regpar(3) = x0	!coordinate of first point in in x direction (x0)
	regpar(4) = y0	!coordinate of first point in in y direction (y0)
	regpar(5) = dx	!space step in x direction (dx)
	regpar(6) = dy	!space step in y direction (dy)
	regpar(7) = flag!flag for invalid value (flag)

! in hlv is the layer structure, values indicate the bottom of the layer
! in datetime is the reference date
! in zeta is the water level

	hlv = 0.			!no depth structure for wind file
	datetime = (/19970101,0/)	!reference date is 1.1.1997 - not used

! here we construct a simple 2d matrix water level field

	zeta = 0.
	do iy=1,ny
	  do ix=1,nx
	    ind = ix+iy
	    zeta(ix,iy) = ind / float(nx+ny)
	  end do
	end do

	iunit = 1
	open(iunit,file='zeta.fem',form='formatted',status='unknown')

	dtime = 0.
        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)
        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvddi,zeta)

	close(iunit)

	end

!*******************************************************************

	subroutine regular_wind

! shows how to write a regular wind fem file
!
! more information can be found in file subfemfile.f

	implicit none

	integer, parameter :: lmax = 1
	integer, parameter :: nx = 7		!points in x direction
	integer, parameter :: ny = 8		!points in y direction

	real uw(nx,ny)
	real vw(nx,ny)
	real pr(nx,ny)
	real hlv(lmax)
	real regpar(7)
	integer ilhkv(1)	!not needed, only one level
	integer datetime(2)

	integer nvers,np,nvar,ntype,nlvddi,iformat
	integer iunit
	double precision dtime
	real dxy,dx,dy,x0,y0,x1,y1
	real hd,flag
	character*20 string

	dxy = 0.1
	x0 = 20.4
	y0 = 54.7
	x1 = x0 + (nx-1)*dxy
	y1 = y0 + (ny-1)*dxy
	dx = dxy
	dy = dxy

	nvers = 0	!version of femfile, 0 for latest
	np = nx*ny	!total number of points
	ntype = 11	!we want with date and regular information
	nvar = 3	!wind file has 3 variables: wind-x, wind-y, pressure
	nlvddi = lmax	!vertical dimension, must be 1 for 2D
	ilhkv = lmax	!number of layers per point, not needed
	hd = 15.	!total depth, not needed
	iformat = 1	!file should be formatted
	flag = -999.

! set array with information on regular grid
! points run from x0 to x0+(nx-1)*dx and y0 to y0+(ny-1)*dy
! values are stored row-wise from left to right starting from the lowest row
! first point is [x0,y0], second [x0+dx,y0], etc..

	regpar(1) = nx	!number of points in x direction (nx)
	regpar(2) = ny	!number of points in y direction (ny)
	regpar(3) = x0	!coordinate of first point in in x direction (x0)
	regpar(4) = y0	!coordinate of first point in in y direction (y0)
	regpar(5) = dx	!space step in x direction (dx)
	regpar(6) = dy	!space step in y direction (dy)
	regpar(7) = flag!flag for invalid value (flag)

! in hlv is the layer structure, values indicate the bottom of the layer
! in datetime is the reference date

	hlv = 0.			!no depth structure for wind file
	datetime = (/19970101,0/)	!reference date is 1.1.1997

	uw = -10.	!spatially constant wind from SE (Scirocco)
	vw = +10.
	pr = 101325	!constant pressure of 1013.25 mbar

	iunit = 1
	open(iunit,file='wind.fem',form='formatted',status='unknown')

	dtime = 0.
        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)
        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,'wind velocity in x [m/s]'
     +                          ,ilhkv,hd
     +                          ,nlvddi,uw)
        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,'wind velocity in y [m/s]'
     +                          ,ilhkv,hd
     +                          ,nlvddi,vw)
        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,'pressure (atmospheric) [Pa]'
     +                          ,ilhkv,hd
     +                          ,nlvddi,pr)

	dtime = 86400.
        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)
        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,'wind velocity in x [m/s]'
     +                          ,ilhkv,hd
     +                          ,nlvddi,uw)
        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,'wind velocity in y [m/s]'
     +                          ,ilhkv,hd
     +                          ,nlvddi,vw)
        call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,'pressure (atmospheric) [Pa]'
     +                          ,ilhkv,hd
     +                          ,nlvddi,pr)

	close(iunit)

	end

!*******************************************************************

	subroutine regular_conz_2d

! shows how to write a regular conz fem file in 2d
!
! only one time record is written -> can be used for initial condition
!
! more information can be found in file subfemfile.f

	implicit none

	integer, parameter :: lmax = 1
	integer, parameter :: nx = 5		!points in x direction
	integer, parameter :: ny = 17		!points in y direction

	integer, parameter :: nvar = 1		!number of variables
	character*20 :: file = 'conz2d_1.fem'

	!integer, parameter :: nvar = 3		!number of variables
	!character*20 :: file = 'conz2d_3.fem'

	real val(nx,ny,nvar)
	real val1(nx,ny)
	real hlv(lmax)
	real regpar(7)
	integer ilhkv(1)	!not needed, only one level
	integer datetime(2)

	integer nvers,np,ntype,nlvddi,iformat
	integer iunit,iv,ix,iy,ind
	double precision dtime
	real dxy,dx,dy,x0,y0,x1,y1
	real hd,flag
	real y,ycrit
	character*20 string

	dxy = 1000.
	x0 = 764000.
	y0 = 4132000.
	x1 = x0 + (nx-1)*dxy
	y1 = y0 + (ny-1)*dxy
	dx = dxy
	dy = dxy

	ycrit = 4140000.

	write(6,*) x0,y0,x1,y1,dx,dy

	nvers = 0	!version of femfile, 0 for latest
	np = nx*ny	!total number of points
	ntype = 11	!we want with date and regular information
	nlvddi = lmax	!vertical dimension, must be 1 for 2D
	ilhkv = lmax	!number of layers per point, not needed
	hd = 15.	!total depth, not needed
	iformat = 1	!file should be formatted
	string = 'concentration'
	flag = -999.

! set array with information on regular grid
! points run from x0 to x0+(nx-1)*dx and y0 to y0+(ny-1)*dy
! values are stored row-wise from left to right starting from the lowest row
! first point is [x0,y0], second [x0+dx,y0], etc..

	regpar(1) = nx	!number of points in x direction (nx)
	regpar(2) = ny	!number of points in y direction (ny)
	regpar(3) = x0	!coordinate of first point in in x direction (x0)
	regpar(4) = y0	!coordinate of first point in in y direction (y0)
	regpar(5) = dx	!space step in x direction (dx)
	regpar(6) = dy	!space step in y direction (dy)
	regpar(7) = flag!flag for invalid value (flag)

! in hlv is the layer structure, values indicate the bottom of the layer
! in datetime is the reference date
! in val are the concentration values for each layer
! the concentration values refer to the center of the layer

	hlv = 0.			!no depth structure for 2d file
	datetime = (/19970101,0/)	!reference date is 1.1.1997

! here we construct the concentration (both multi and single)

	val = 0.
	do iv=1,nvar
	  do iy=1,ny
	    do ix=1,nx
	      ind = ix+iy
	      if( iv == 2 ) ind = ix + ny - iy + 1
	      if( iv == 3 ) ind = iy + nx - ix + 1
	      val(ix,iy,iv) = ind / float(nx+ny)
	    end do
	  end do
	end do

	val = 0.
	do iy=1,ny
	  y = y0 + (iy-1)*dxy
	  do ix=1,nx
	    if( y > ycrit ) val(ix,iy,1) = 1.
	  end do
	end do

	iunit = 1
	open(iunit,file=file,form='formatted',status='unknown')

	dtime = 0.
        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)
	do iv=1,nvar
          call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvddi,val(:,:,iv))
	end do

	close(iunit)

	end

!*******************************************************************

	subroutine regular_conz_3d

! shows how to write a regular conz fem file in 3d
!
! only one time record is written -> can be used for initial condition
!
! more information can be found in file subfemfile.f

	implicit none

	integer, parameter :: lmax = 5
	integer, parameter :: nx = 13		!points in x direction
	integer, parameter :: ny = 16		!points in y direction

	integer, parameter :: nxy = nx * ny

	!integer, parameter :: nvar = 1		!number of variables
	!character*20 :: file = 'conz2d_1.fem'

	integer, parameter :: nvar = 3		!number of variables
	character*20 :: file = 'conz3d_3.fem'

	real val(lmax,nx,ny,nvar)
	real hlv(lmax)
	real regpar(7)
	real depth(nxy)
	integer ilhkv(nxy)	!level structure
	integer datetime(2)

	integer nvers,ntype,nlvddi,iformat
	integer iunit,np,iv,ix,iy,ind,l
	double precision dtime
	real dxy,dx,dy,x0,y0,x1,y1
	real hd,flag
	real rnvar,rlmax,rnxy
	character*20 string

	dxy = 0.1
	x0 = 20.4
	y0 = 54.7
	x1 = x0 + (nx-1)*dxy
	y1 = y0 + (ny-1)*dxy
	dx = dxy
	dy = dxy
	flag = -999.

	!write(6,*) x0,y0,x1,y1,dx,dy

	nvers = 0	!version of femfile, 0 for latest
	np = nxy	!total number of points
	ntype = 11	!we want with date and regular information
	nlvddi = lmax	!vertical dimension, must be 1 for 2D
	ilhkv = lmax	!number of layers per point
	hd = flag	!total depth, not needed
	iformat = 1	!file should be formatted
	string = 'concentration'

! set array with information on regular grid
! points run from x0 to x0+(nx-1)*dx and y0 to y0+(ny-1)*dy
! values are stored row-wise from left to right starting from the lowest row
! first point is [x0,y0], second [x0+dx,y0], etc..

	regpar(1) = nx	!number of points in x direction (nx)
	regpar(2) = ny	!number of points in y direction (ny)
	regpar(3) = x0	!coordinate of first point in in x direction (x0)
	regpar(4) = y0	!coordinate of first point in in y direction (y0)
	regpar(5) = dx	!space step in x direction (dx)
	regpar(6) = dy	!space step in y direction (dy)
	regpar(7) = flag!flag for invalid value (flag)

! in hlv is the layer structure, values indicate the bottom of the layer
! in datetime is the reference date
! in val are the concentration values for each layer
! the concentration values refer to the center of the layer

	hlv = 0.
	if( lmax > 1 ) then
	  do l=1,lmax
	    hlv(l) = -float(l)/lmax	!sigma levels
	  end do
	end if
	datetime = (/19970101,0/)	!reference date is 1.1.1997
	depth = hd

	rnvar = nvar
	rlmax = lmax
	rnxy = nx + ny

! here we construct a 3d matrix that can be verified - no real world example

	val = 0.
	do iv=1,nvar
	  do iy=1,ny
	    do ix=1,nx
	      ind = ix+iy
	      if( iv == 2 ) ind = ix + ny - iy + 1
	      if( iv == 3 ) ind = iy + nx - ix + 1
	      do l=1,lmax
	        !val(l,ix,iy,iv) =  ind / rnxy
	        val(l,ix,iy,iv) = (l/rlmax) * ind / rnxy
	        !val(l,ix,iy,iv) = ((iv-1)/rnvar) + (l/rlmax) / rnvar
	      end do
	    end do
	  end do
	end do

	iunit = 1
	open(iunit,file=file,form='formatted',status='unknown')

	dtime = 0.
        call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)
	do iv=1,nvar
          call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,depth
     +                          ,nlvddi,val(:,:,:,iv))
	end do

	close(iunit)

	end

!*******************************************************************

