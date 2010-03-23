c
c reads and writes unit 47 file
c
c***********************************************************
c
c               |                   |
c               |                   |
c               |                   | i,j+1
c               |                   |
c               |                   |
c               |        i,j        |       i+1,j
c-----------------------------------------------------
c               |                   |
c               |                   |
c               |                   |
c         i-1,j |        i,j        | i,j
c               |                   |
c               |                   |
c               |        i,j-1      | 
c-----------------------------------------------------
c               |                   |
c               |                   |
c               |                   |
c               |                   |
c               |                   |
c               |                   |
c
c**********************************************************************

	subroutine wf47(nx,ny,x0,y0,dx,dy)

c writes header of matrix file (unit 47)

	implicit none

	integer nx,ny
	real x0,y0,dx,dy

	write(47) nx,ny,x0,y0,dx,dy

	end

c**********************************************************************

	subroutine rf47(nx,ny,x0,y0,dx,dy)

c reads header of matrix file (unit 47)

	implicit none

	integer nx,ny
	real x0,y0,dx,dy

	rewind(47)
	read(47) nx,ny,x0,y0,dx,dy

	end

c**********************************************************************

	subroutine wr47(ndim,mmax,nbox,nsect,i0,j0,i1,j1,ie,id
     +		,tx,ty,ar2,alen,aa,ipm,ipam)

c writes data record of matrix file (unit 47)

	implicit none

	integer ndim,mmax,nbox,nsect,i0,j0,i1,j1,ie,id
	real tx(3),ty(3)
	real ar2(4,1)
	real alen(1)
	real aa(4,1)
	integer ipm(3,1)
	integer ipam(2,1) 

	integer i,j

	write(47) mmax,nbox,nsect,i0,j0,i1,j1,ie,id
	write(47) (tx(i),ty(i),i=1,3)			!triangle vertices
	write(47) ((ar2(i,j),i=1,4),j=1,nsect)
	write(47) (alen(j),j=1,nsect)			!length of intersect.
	write(47) ((aa(i,j),i=1,4),j=1,nbox)
	write(47) ((ipm(i,j),i=1,3),j=1,nsect)		!pointers sect -> boxes
	write(47) ((ipam(i,j),i=1,2),j=1,nbox)		!area pointer

	end

c**********************************************************************

	subroutine rd47(ndim,mmax,nbox,nsect,i0,j0,i1,j1,ie,id
     +		,tx,ty,ar2,alen,aa,ipm,ipam,bstop)

c reads data record of matrix file (unit 47)

	implicit none

	integer ndim,mmax,nbox,nsect,i0,j0,i1,j1,ie,id
	real tx(3),ty(3)
	real ar2(4,1)
	real alen(1)
	real aa(4,1)
	integer ipm(3,1)
	integer ipam(2,1) 
	logical bstop

	integer i,j

	bstop = .false.

        read(47,end=99) mmax,nbox,nsect,i0,j0,i1,j1,ie,id

	if( ndim .lt. nsect ) stop 'error stop rd47: ndim'
	if( ndim .lt. nbox ) stop 'error stop rd47: ndim'

        read(47) (tx(i),ty(i),i=1,3)
        read(47) ((ar2(i,j),i=1,4),j=1,nsect)
	read(47) (alen(j),j=1,nsect)
        read(47) ((aa(i,j),i=1,4),j=1,nbox)
        read(47) ((ipm(i,j),i=1,3),j=1,nsect)
        read(47) ((ipam(i,j),i=1,2),j=1,nbox)

	return
   99	continue
	bstop = .true.
	return
	end

c**********************************************************************

