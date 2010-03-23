c
c $Id: mkgeom.f,v 1.3 2001/11/16 07:35:43 georg Exp $
c
	program mkgeom

c makes geometry

	implicit none

	integer nxdim,nydim,ibas
c	parameter(nxdim=50,nydim=50,ibas=1)
c	parameter(nxdim=40,nydim=40,ibas=1)
c	parameter(nxdim=8,nydim=48,ibas=1)
	parameter(nxdim=5,nydim=5,ibas=1)

	integer idep(0:nxdim+1,0:nydim+1)
	integer node(0:nxdim,0:nydim)

	character*72 title
	integer ix,iy
	integer nnode,nelem
	real x,y
	real fact
	integer idepc

c--------------------------------------------------------------
c set constants
c--------------------------------------------------------------

	if( ibas .eq. 1 ) then
	  title = '0   (FEM-TITLE)   regular basin'
	  idepc = 200
	  fact = 4000.
	else
	  idepc = 1
	  fact = 1
	end if

c--------------------------------------------------------------
c initialize arrays
c--------------------------------------------------------------

	do iy=0,nydim+1
	  do ix=0,nxdim+1
	    idep(ix,iy) = 0
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
	    idep(ix,iy) = idepc
	  end do
	end do
	
c--------------------------------------------------------------
c set nodes
c--------------------------------------------------------------

	nnode = 0

	do iy=0,nydim
	  do ix=0,nxdim
	    call mknode(ix,iy,node,idep,nxdim,nydim,nnode)
	  end do
	end do

c--------------------------------------------------------------
c write title
c--------------------------------------------------------------

	write(6,*) 
	write(6,*) title
	write(6,*) 

c--------------------------------------------------------------
c write nodes
c--------------------------------------------------------------

	do iy=0,nydim
	  do ix=0,nxdim
	    if( node(ix,iy) .gt. 0 ) then
		x = ix * fact
		y = iy * fact
		!call mkxy1(ix,iy,node,nxdim,nydim,x,y)
		write(6,'(3i6,2f14.4)') 1,node(ix,iy),0,x,y
	    end if
	  end do
	end do

	write(6,*) 

c--------------------------------------------------------------
c write elements
c--------------------------------------------------------------

	nelem = 0

	do iy=1,nydim
	  do ix=1,nxdim
	    if( idep(ix,iy) .gt. 0 ) then
		call wrtri(ix,iy,node,idep,nxdim,nydim,0,nelem)
		call wrtri(ix,iy,node,idep,nxdim,nydim,1,nelem)
	    end if
	  end do
	end do

	write(6,*) 

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c******************************************************************

	subroutine wrtri(ix,iy,node,idep,nxdim,nydim,iupper,nelem)

c writes triangle

	implicit none

	integer ix,iy
	integer nxdim,nydim
	integer iupper
	integer nelem
	integer node(0:nxdim,0:nydim)
	integer idep(0:nxdim+1,0:nydim+1)

	integer ielem,i
	integer nodtri(3)

	nelem = nelem + 1
	ielem = 100*iy + 2*ix - 1 + iupper
	call mktri(ix,iy,node,nxdim,nydim,iupper,nodtri)
	write(6,'(7i6,i9)') 2,ielem,0,3
     +		,(nodtri(i),i=1,3),idep(ix,iy)

	end
  
c******************************************************************

	subroutine mktri(ix,iy,node,nxdim,nydim,iupper,nodtri)

c makes triangle nodes

	implicit none

	integer ix,iy
	integer nxdim,nydim
	integer iupper
	integer nodtri(3)
	integer node(0:nxdim,0:nydim)

	if( iupper .eq. 0 ) then	!make lower triangle
	  nodtri(1) = node(ix-1,iy-1)
	  nodtri(2) = node(ix,iy-1)
	  nodtri(3) = node(ix,iy)
	else				!make upper triangle
	  nodtri(1) = node(ix-1,iy-1)
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

	logical bnode
	integer inode

	bnode = .false.
	if( idep(ix,iy) .gt. 0 ) bnode = .true.
	if( idep(ix+1,iy) .gt. 0 ) bnode = .true.
	if( idep(ix,iy+1) .gt. 0 ) bnode = .true.
	if( idep(ix+1,iy+1) .gt. 0 ) bnode = .true.

	if( bnode ) then
		nnode = nnode + 1
		inode = 100*iy + ix + 1
		node(ix,iy) = inode
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

