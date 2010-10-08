c
c $Id: windrose.f,v 1.1 2000/05/26 12:27:50 georg Exp $
c
c revision log :
c
c 13.04.1999	ggu	written from scratch
c
c*************************************************************

	program windrose

c plots simulation

	implicit none

c parameters
	integer ndim
	parameter (ndim=40)

	real wx(ndim)
	real wy(ndim)

	integer n
	real rmax

	n = 1
	rmax = 20.

c-----------------------
	call qopen
c-----------------------

	do while( n .gt. 0 )

	  call rdw(n,ndim,wx,wy)
	  call rose(n,wx,wy,rmax)

	end do

c-----------------------
	call qclose
c-----------------------

	end

c*****************************************************************

	subroutine rdw(n,ndim,wx,wy)

	implicit none

	integer n,ndim
	real wx(1), wy(1)

	integer i,it

	do i=1,ndim
	  read(5,*,end=1) it,wx(i),wy(i)
	end do

    1	continue

	n = i - 1

	end

c*****************************************************************

	subroutine rose(n,wx,wy,rmax)

	implicit none

	integer n
	real wx(1), wy(1)
	real rmax

	integer i

	if( n .le. 0 ) return

	write(6,*) n,rmax

	call qstart

	call qsetvp(1.,1.,15.,15.)
	call qworld(-rmax,-rmax,rmax,rmax)

	call qline(-rmax,-rmax,rmax,-rmax)
	call qline(rmax,-rmax,rmax,rmax)
	call qline(rmax,rmax,-rmax,rmax)
	call qline(-rmax,rmax,-rmax,-rmax)

	call qline(0.,-rmax,0.,rmax)
	call qline(-rmax,0.,rmax,0.)

	do i=1,n
	  call qline(0.,0.,wx(i),wy(i))
	end do

	call qend

	end

c*****************************************************************

