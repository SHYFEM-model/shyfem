
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017,2019  Georg Umgiesser
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

! revision log :
!
! 05.12.2017	ggu	changed VERS_7_5_39
! 14.02.2019	ggu	changed VERS_7_5_56

	program mkgeomr

c------------------------------------------------------------------
c
c makes geometry of regular basins (equilateral triangles)
c
c------------------------------------------------------------------

	implicit none

	integer nx,ny
	real dl,h

c--------------------------------------------------------------
c regular grid
c--------------------------------------------------------------

	nx = 80			!points in x
	ny = 20			!points in y
	dl = 100.		!side length of triangle
	h = 10.			!depth

	call equi(nx,ny,dl,h)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c******************************************************************

	subroutine equi(nx,ny,dl,h)

	implicit none

	integer nx,ny
	real dl,h

	real dx,dy
	real x0,y0,x,y
	integer ix,iy
	integer n,ii,i,ie,nel
	integer ndim

	real xv(nx+1,ny)
	real yv(nx+1,ny)
	integer node(nx+1,ny)
	integer nen3v(3,2*nx*ny)

	ndim = 2*nx*ny

	dx = dl
	dy = dl*sqrt(3.)/2.

	y0 = -dy
	n = 0
	ie = 0

	do iy = 1,ny
	  y0 = y0 + dy
	  x0 = 0.
	  ii = mod(iy,2)
	  if( ii .eq. 1 ) then
	    do ix=1,nx
		n = n + 1
		node(ix,iy) = n
		xv(ix,iy) = x0
		yv(ix,iy) = y0
		x0 = x0 + dx
	    end do
	    ix = nx + 1
	    node(ix,iy) = 0
	    if( iy > 1 ) then
	      do i=1,nx-1
	        call mkt(ie,nen3v,node(i,iy),node(i,iy-1),node(i+1,iy-1))
	        call mkt(ie,nen3v,node(i,iy),node(i+1,iy-1),node(i+1,iy))
	      end do
	      i = nx
	      call mkt(ie,nen3v,node(i,iy),node(i,iy-1),node(i+1,iy-1))
	      if( ie .gt. ndim ) stop 'error stop: ndim'
	    end if
	  else
	    ix = 1
	    n = n + 1
	    node(ix,iy) = n
	    xv(ix,iy) = x0
	    yv(ix,iy) = y0
	    x0 = x0 + dx/2.
	    do ix=2,nx
		n = n + 1
		node(ix,iy) = n
		xv(ix,iy) = x0
		yv(ix,iy) = y0
		x0 = x0 + dx
	    end do
	    ix = nx + 1
	    n = n + 1
	    node(ix,iy) = n
	    x0 = (nx-1)*dx
	    xv(ix,iy) = x0
	    yv(ix,iy) = y0
	    if( iy > 1 ) then
	      i = 1
	      call mkt(ie,nen3v,node(i,iy),node(i,iy-1),node(i+1,iy))
	      do i=2,nx
	        call mkt(ie,nen3v,node(i,iy),node(i-1,iy-1),node(i,iy-1))
	        call mkt(ie,nen3v,node(i,iy),node(i,iy-1),node(i+1,iy))
	      end do
	      if( ie .gt. ndim ) stop 'error stop: ndim'
	    end if
	  end if
	end do
	  
	do iy=1,ny
	  do ix=1,nx+1
	    if( node(ix,iy) .gt. 0 ) then
		x = xv(ix,iy)
		y = yv(ix,iy)
		write(6,'(i1,2i7,2f14.4)') 1,node(ix,iy),0,x,y
	    end if
	  end do
	end do

	nel = ie
	do ie=1,nel
	  write(6,'(i1,6i7,f10.2)') 2,ie,0,3,(nen3v(ii,ie),ii=1,3),h
	end do

	end

c******************************************************************

	subroutine mkt(ie,nen3v,n1,n2,n3)

	implicit none

	integer ie,n1,n2,n3
	integer nen3v(3,1)

	ie = ie + 1
	nen3v(1,ie) = n1
	nen3v(2,ie) = n2
	nen3v(3,ie) = n3

	end

c******************************************************************

