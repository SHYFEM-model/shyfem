
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

	program mkgeom

c------------------------------------------------------------------
c
c makes geometry of regular basins
c
c nxdim,nydim	number of cells in x/y dir
c xyfact	factor for x/y coordinates
c zfact		inverse factor for z coordinates 
c			-> zfact=1000 => idep=100 -> z=0.1 m
c idepc		constant depth for basin (may be modified in subs)
c
c inumb controls numbering of nodes
c
c	inumb > 0  x first      
c	inumb < 0  y first
c	abs(inumb) > 1	   every col/row has inumb more
c
c etype controls type of triangles
c
c	etype = 1		upward triangles
c	etype = -1		downward triangles
c	etype = 0		alternating triangles
c
c------------------------------------------------------------------

	implicit none

c nxdim,nydim are number of triangles

	integer nxdim,nydim
	real xyfact,zfact
	integer idepc
	integer inumb,etype
	character*72 title

c--------------------------------------------------------------
c customize parameters of regular grid
c--------------------------------------------------------------

c	parameter(title='0   (FEM-TITLE)   regular basin'
c	parameter(nxdim=50,nydim=50)
c	parameter(nxdim=40,nydim=40)
c	parameter(nxdim=8,nydim=48)
c	parameter(nxdim=5,nydim=5)
c	parameter(idepc=200,xyfact=4000.)
c	parameter(idepc=1,xyfact=1.)
c	parameter(inumb=1,etype=1)
c	parameter(zfact=1.)

c	parameter(title='0   (FEM-TITLE)   basin for chao simulations')
c	parameter(nxdim=37,nydim=45)
c	parameter(idepc=15,xyfact=3000.,inumb=100,etype=0)
c	parameter(zfact=1.)

c	parameter(title='0   (FEM-TITLE)   channel with step')
c	parameter(nxdim=100,nydim=20)
c	parameter(idepc=10,xyfact=250.,inumb=-1,etype=0)
c	parameter(zfact=1.)

c	parameter(title='0   (FEM-TITLE)   idealized inlet')
c	parameter(nxdim=144,nydim=188)
c	parameter(idepc=10,xyfact=100.,inumb=-1000,etype=0,zfact=0.001)

c	parameter(title='0   (FEM-TITLE)   regular basin for diffus')
c	parameter(nxdim=100,nydim=100)
c	parameter(idepc=10,xyfact=100.,inumb=1000,etype=0,zfact=1.)

	parameter(title='0   (FEM-TITLE)   regular basin for mpi')
	parameter(nxdim=500,nydim=500)
	parameter(idepc=10,xyfact=100.,inumb=-1,etype=0,zfact=1.)

c--------------------------------------------------------------
c do not change anything below here
c--------------------------------------------------------------

	integer idep(0:nxdim+1,0:nydim+1)
	integer node(0:nxdim,0:nydim)

	integer ix,iy
	integer nnode,nelem
	integer efact
	integer k,n,nn
	real x,y
	character*80 outfile
	integer, parameter :: nfreq = 10000

c--------------------------------------------------------------
c statement functions -> use one of them
c--------------------------------------------------------------

	integer eltype,i

c	eltype(i) = 1				!upward triangles
c	eltype(i) = -1				!downward triangles
c	eltype(i) = 1 - 2*mod(i,2)		!alternating triangles
c
c all of the above
c
c	etype = 1		upward triangles
c	etype = -1		downward triangles
c	etype = 0		alternating triangles

	eltype(i) = 1 - (1-abs(etype))*2*mod(i,2) + 2*min(0,etype)

c	1-abs(etype)	=>	1 for etype=0, 0 else
c	min(0,etype)	=>	0 for etype=0 or 1, -1 else

c--------------------------------------------------------------
c initialize arrays
c--------------------------------------------------------------

	do iy=0,nydim+1
	  do ix=0,nxdim+1
	    idep(ix,iy) = idepc
	  end do
	end do

	do iy=0,nydim
	  do ix=0,nxdim
	    node(ix,iy) = 0
	  end do
	end do

c--------------------------------------------------------------
c set geometry
c--------------------------------------------------------------

	do iy=1,nydim
	  do ix=1,nxdim
	    !call mkgeo0(ix,iy,idep,nxdim,nydim)
	    !call mkgeo1(ix,iy,idep,nxdim,nydim)
	    !call mkgeo_chao(ix,iy,idep,nxdim,nydim)
	    !call mkgeo_step(ix,iy,idep,nxdim,nydim)
	    !call mkgeo_idinlet(ix,iy,idep,nxdim,nydim)
	  end do
	end do
	
c--------------------------------------------------------------
c set nodes
c--------------------------------------------------------------

	nnode = 0

	do iy=1,nydim
	  do ix=1,nxdim
	    call mknode(ix,iy,node,idep,nxdim,nydim,nnode)
	  end do
	end do

	call renumber_node(inumb,node,nxdim,nydim)

c--------------------------------------------------------------
c open file for write
c--------------------------------------------------------------

	outfile = 'regular.grd'
	open(1,file=outfile,status='unknown',form='formatted')

c--------------------------------------------------------------
c write title
c--------------------------------------------------------------

	write(1,*) 
	write(1,'(a)') title
	write(1,*) 

c--------------------------------------------------------------
c write nodes
c--------------------------------------------------------------

	nn = 0

	do ix=0,nxdim
	  do iy=0,nydim
	    if( node(ix,iy) .gt. 0 ) then
		nn = nn + 1
		x = ix * xyfact
		y = iy * xyfact
		k = node(ix,iy)
		!call mkxy1(ix,iy,node,nxdim,nydim,x,y)
		write(1,'(i1,i9,i5,2f16.4)') 1,k,0,x,y
		if( mod(nn,nfreq) == 0 ) write(6,*) 'nodes: ',nn
	    end if
	  end do
	end do

	write(1,*) 

c--------------------------------------------------------------
c write elements
c--------------------------------------------------------------

	nelem = 0

	do iy=1,nydim
	  do ix=1,nxdim
	    if( idep(ix,iy) .gt. 0 ) then
		efact = eltype(ix+iy)
		call wrtri(ix,iy,node,idep,nxdim,nydim,1*efact,nelem,zfact)
		call wrtri(ix,iy,node,idep,nxdim,nydim,2*efact,nelem,zfact)
		if( mod(nelem,nfreq) == 0 ) write(6,*) 'elems: ',nelem
	    end if
	  end do
	end do

	write(1,*) 

c--------------------------------------------------------------
c write final message
c--------------------------------------------------------------

	write(6,*) 'total number of nodes: ',nn
	write(6,*) 'total number of elements: ',nelem
	write(6,*) 'regular grid has been written to ',trim(outfile)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c******************************************************************

	subroutine wrtri(ix,iy,node,idep,nxdim,nydim,ietype,nelem,zfact)

c writes triangle

	implicit none

	integer ix,iy
	integer nxdim,nydim
	integer ietype
	integer nelem
	real zfact
	integer node(0:nxdim,0:nydim)
	integer idep(0:nxdim+1,0:nydim+1)

	integer ielem,i
	integer nodtri(3)
	real depth

	nelem = nelem + 1
	ielem = 100*iy + 2*ix - 2 + abs(ietype)
	ielem = nelem
	depth = zfact*idep(ix,iy)

	call mktri(ix,iy,node,nxdim,nydim,ietype,nodtri)
	write(1,'(i1,6i9,f16.4)') 2,ielem,0,3
     +		,(nodtri(i),i=1,3),depth

	end
  
c******************************************************************

	subroutine mktri(ix,iy,node,nxdim,nydim,ietype,nodtri)

c makes nodes for triangle

	implicit none

	integer ix,iy
	integer nxdim,nydim
	integer ietype
	integer nodtri(3)
	integer node(0:nxdim,0:nydim)

	if( ietype .eq. 1 ) then		!make lower triangle
	  nodtri(1) = node(ix-1,iy-1)
	  nodtri(2) = node(ix,iy-1)
	  nodtri(3) = node(ix,iy)
	else if( ietype .eq. 2 ) then		!make upper triangle
	  nodtri(1) = node(ix-1,iy-1)
	  nodtri(2) = node(ix,iy)
	  nodtri(3) = node(ix-1,iy)
	else if( ietype .eq. -1 ) then		!make lower triangle
	  nodtri(1) = node(ix-1,iy-1)
	  nodtri(2) = node(ix,iy-1)
	  nodtri(3) = node(ix-1,iy)
	else if( ietype .eq. -2 ) then		!make upper triangle
	  nodtri(1) = node(ix,iy-1)
	  nodtri(2) = node(ix,iy)
	  nodtri(3) = node(ix-1,iy)
	end if

	end

c******************************************************************

	subroutine mknode(ix,iy,node,idep,nxdim,nydim,nnode)

c makes node

	implicit none

	integer ix,iy
	integer nxdim,nydim
	integer nnode
	integer node(0:nxdim,0:nydim)
	integer idep(0:nxdim+1,0:nydim+1)

	if( idep(ix,iy) .le. 0 ) return

	call new_node(ix,iy,node,nxdim,nydim,nnode)
	call new_node(ix-1,iy,node,nxdim,nydim,nnode)
	call new_node(ix,iy-1,node,nxdim,nydim,nnode)
	call new_node(ix-1,iy-1,node,nxdim,nydim,nnode)

	end

c******************************************************************

	subroutine new_node(ix,iy,node,nxdim,nydim,nnode)

c adds new node

	implicit none

	integer ix,iy
	integer nxdim,nydim
	integer nnode
	integer node(0:nxdim,0:nydim)

	if( node(ix,iy) .gt. 0 ) return		!already there

	nnode = nnode + 1
	node(ix,iy) = 100*iy + ix + 1
	node(ix,iy) = nnode

	end

c******************************************************************

	subroutine renumber_node(ihow,node,nxdim,nydim)

c renumber nodes in an intelligent way

	implicit none

	integer ihow
	integer nxdim,nydim
	integer node(0:nxdim,0:nydim)

	integer istep,inode
	integer ix,iy,n

	inode = 0
	istep = abs(ihow)

	if( ihow .gt. 0 ) then		!in x direction
	  do iy=0,nydim
	    do ix=0,nxdim
	      if( node(ix,iy) .gt. 0 ) then
	        inode = inode + 1
	        n = inode
	        if( istep .gt. 1 ) n = istep*(iy+1) + ix + 1
	        node(ix,iy) = n
	      end if
	    end do
	  end do
	else if( ihow .lt. 0 ) then
	  do ix=0,nxdim
	    do iy=0,nydim
	      if( node(ix,iy) .gt. 0 ) then
	        inode = inode + 1
	        n = inode
	        if( istep .gt. 1 ) n = istep*(ix+1) + iy + 1
	        node(ix,iy) = n
	      end if
	    end do
	  end do
	else
	  stop 'error stop renumber_node: ihow'
	end if
	
	end

c******************************************************************

	function linear(i,istart,iend,rstart,rend)

c linear interpolation

	implicit none

	real linear
	integer i,istart,iend
	real rstart,rend

	real dr,di

	dr = rend - rstart
	di = iend - istart

	linear = rstart + ( dr * ( i - istart ) ) / di

	end

c******************************************************************

	function quadrat(i,istart,iend,rstart,rend)

c quadratic interpolation

	implicit none

	real quadrat
	integer i,istart,iend
	real rstart,rend

	real dr,di,a

	dr = rend - rstart
	di = iend - istart

	a = dr / (di*di)

	quadrat = rstart + a * (i-istart)**2

	end

c******************************************************************
c******************************************************************
c
c special bathymetry and geometry
c
c******************************************************************
c******************************************************************

	subroutine mkgeo0(ix,iy,idep,nxdim,nydim)

	implicit none

	integer ix,iy
	integer nxdim,nydim
	integer idep(0:nxdim+1,0:nydim+1)

	real linear

	if( iy .ge. 24 ) then
	    idep(ix,iy) = 100
	else
	    idep(ix,iy) = linear(iy,24,1,100.,1000.)
	end if

	end

c******************************************************************

	subroutine mkgeo1(ix,iy,idep,nxdim,nydim)

c set special bathymetry (layer depth is 3 meters)

	implicit none

	integer ix,iy
	integer nxdim,nydim
	integer idep(0:nxdim+1,0:nydim+1)

	integer i

	i = ix - 19
c	i = i / 5
	i = i / 3
	if( i .gt. 0 ) then
	      idep(ix,iy) = idep(ix,iy) + i * 3
	end if

	end

c******************************************************************

	subroutine mkgeo_chao(ix,iy,idep,nxdim,nydim)

c set special geometry for chao

	implicit none

	integer ix,iy
	integer nxdim,nydim
	integer idep(0:nxdim+1,0:nydim+1)

	if( ix .le. 17 ) then
	  if( iy .le. 25 .or. iy .ge. 31 ) then
	      idep(ix,iy) = 0
	  end if
	end if

	end

c******************************************************************

	subroutine mkgeo_step(ix,iy,idep,nxdim,nydim)

c set special geometry for channel step

	implicit none

	integer ix,iy
	integer nxdim,nydim
	integer idep(0:nxdim+1,0:nydim+1)

	if( ix .le. 40 ) then
	  if( iy .le. 10 ) then
	      idep(ix,iy) = 0
	  end if
	end if


	end

c******************************************************************

	subroutine mkgeo_idinlet(ix,iy,idep,nxdim,nydim)

c set special geometry for idealized inlet

	implicit none

	integer ix,iy
	integer nxdim,nydim
	integer idep(0:nxdim+1,0:nydim+1)

	integer ix1,ix2
	real z1,z2

	if( ix .le. 20 ) then
	  idep(ix,iy) = 0
	else if( ix .le. 65 ) then
	  z1 = 0.
	  z2 = 12.5
	  ix1 = 20
	  ix2 = 65
	  idep(ix,iy) = 1000*( z1 + (ix-ix1)*(z2-z1)/(ix2-ix1) )
	else if( ix .le. 130 ) then
	  z1 = 12.5
	  z2 = 20.
	  ix1 = 65
	  ix2 = 130
	  idep(ix,iy) = 1000*( z1 + (ix-ix1)*(z2-z1)/(ix2-ix1) )
	else
	  z1 = 20.
	  z2 = 24.
	  ix1 = 130
	  ix2 = 144
	  idep(ix,iy) = 1000*( z1 + (ix-ix1)*(z2-z1)/(ix2-ix1) )
	end if
	  
	if( ix .le. 60 ) then

	  if( ix .le. 40 ) then
	    if( iy .eq. 90 .or. iy .eq. 96 ) then		!dikes
	      idep(ix,iy) = 0
	    end if
	  end if

	  if( iy .gt. 90 .and. iy .lt. 96 ) then	!inlet
	    idep(ix,iy) = 1000 * 12
	  end if

	end if

	end

c******************************************************************

	subroutine mkxy1(ix,iy,node,idep,nxdim,nydim,x,y)

c makes coordinates

	implicit none

	integer ix,iy
	integer nxdim,nydim
	integer node(0:nxdim,0:nydim)
	integer idep(0:nxdim+1,0:nydim+1)
	real x,y

	real fact,rfact
	real linear,quadrat

	rfact = 1000.

	x = ix - 4
	y = iy - 24
	if( y .ge. 0. ) then
	   x = x * rfact
	   y = y * rfact
	else
	   fact = linear(iy,24,0,1.,10.)
	   y = y * rfact * fact
	   fact = quadrat(iy,24,0,1.,10.)
	   x = x * rfact * fact
	end if

	end

c******************************************************************

