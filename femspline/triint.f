c
c $Id: triint.f,v 1.1 1998/08/24 10:25:46 georg Exp $
c
c subroutines for 2 (bi) and 3-dimensional (tri) linear interpolation
c
c contents :
c
c function triint(u,x,y,z)	tri-linear interpolation in cube
c function bilint(u,x,y)	bi-linear interpolation in square
c subroutine exbil		exercise bi-linear interpolation
c subroutine extri		exercise tri-linear interpolation
c function grand()		random numbers
c
c revision log :
c
c 21.08.1998	ggu	tri and bi linear interpolation
c
c***********************************************************************

	function triint(u,x,y,z)

c tri-linear interpolation in cube
c
c (x,y,z) are in the range [0...1]
c u is array of 8 values distributed in the following way
c
c
c
c                  z                      y
c                  ^                     ^
c                  |                    /
c                  |                   /
c                  |     7            /            8
c                  |       *---------------------*
c                  |      /|        /           /|
c                  |     / |       /           / |
c                  |    /  |      /           /  |
c                  |   /   |     /           /	 |                   z = 1
c                  |  /    |    /           /    |
c                  | /     |   /           /     |
c                  |/      |  /           /      |
c                5 *---------------------* 6     |
c                  |       |/            |       |
c                  |     3 *-------------|-------* 4
c                  |      /              |      /
c                  |     /               |     /
c                  |    /                |    /
c                  |   /                 |   /	                     z = 0
c                  |  /                  |  /
c                  | /                   | /
c                  |/                    |/
c                  *---------------------*-----------------------> x
c                1                         2
c

	implicit none

	real triint
	real u(8)
	real x,y,z	![0...1]

	integer i,j
	real acu,val,xy
	real a(8)

	real atri(8,8)
	save atri
	data atri /
     +			 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     +			-1 , 1 , 0 , 0 , 0 , 0 , 0 , 0 ,
     +			-1 , 0 , 1 , 0 , 0 , 0 , 0 , 0 ,
     +			 1 ,-1 ,-1 , 1 , 0 , 0 , 0 , 0 ,
     +			-1 , 0 , 0 , 0 , 1 , 0 , 0 , 0 ,
     +			 1 ,-1 , 0 , 0 ,-1 , 1 , 0 , 0 ,
     +			 1 , 0 ,-1 , 0 ,-1 , 0 , 1 , 0 ,
     +			-1 , 1 , 1 ,-1 , 1 ,-1 ,-1 , 1 
     +		 /

	do i=1,8
	  acu = 0.
	  do j=1,8
	    acu = acu + atri(j,i) * u(j)
	  end do
	  a(i) = acu
	end do

	xy = x * y

	val = a(1) + a(2) * x + a(3) * y + a(4) * xy
     +		+ z * ( a(5) + a(6) * x + a(7) * y + a(8) * xy )

	triint = val

	end

c*************************************************************

	function bilint(u,x,y)

c bi-linear interpolation in square
c
c (x,y) are in the range [0...1]
c u is array of 4 values distributed in the following way
c
c
c
c                  y
c                  ^
c                  |
c                  |
c                  |
c                  |
c                3 |                       4
c                  *---------------------*
c                  |                     |
c                  |                     |
c                  |                     |
c                  |                     |
c                  |                     |
c                  |                     |
c                  |                     |
c                  |                     |
c                  |                     |
c                  *---------------------*-----------------------> x
c                1                         2
c
c

	implicit none

	real bilint
	real u(4)
	real x,y	![0...1]

	integer i,j
	real acu,val
	real a(4)

	real atri(4,4)
	save atri
	data atri /
     +			 1 , 0 , 0 , 0 ,
     +			-1 , 1 , 0 , 0 ,
     +			-1 , 0 , 1 , 0 ,
     +			 1 ,-1 ,-1 , 1
     +		 /

	do i=1,4
	  acu = 0.
	  do j=1,4
	    acu = acu + atri(j,i) * u(j)
	  end do
	  a(i) = acu
	end do

	val = a(1) + a(2) * x + a(3) * y + a(4) * x * y

	bilint = val

	end

c*************************************************************

	subroutine xy2nat(xs,ys,x,y,xyn)

c transforms (x/y) to natural coordinates in triangle

	implicit none

	real xs(3), ys(3)	!coordinates of triangle
	real x,y		!coordinates of point
	real xyn(3)		!natural coordinates of point (return)

	integer i,i2,i3
	real area2,rarea2

	area2 = (xs(2)-xs(1))*(ys(3)-ys(1)) 
     +			- (xs(3)-xs(1))*(ys(2)-ys(1))
	rarea2 = 1. / area2

	do i=1,3
	  i2 = mod(i,3) + 1
	  i3 = mod(i2,3) + 1
	  xyn(i) = (xs(i2)*ys(i3)-xs(i3)*ys(i2)) +
     +		x * (ys(i2)-ys(i3)) + y * (xs(i3)-xs(i2))
	  xyn(i) = xyn(i) * rarea2
	end do

	end

c*************************************************************

	subroutine nat2xy(xs,ys,x,y,xyn)

c transforms natural coordinates in triangle to (x/y)

	implicit none

	real xs(3), ys(3)	!coordinates of triangle
	real x,y		!coordinates of point (return)
	real xyn(3)		!natural coordinates of point

	integer i

	x = 0.
	y = 0.
	do i=1,3
	  x = x + xs(i) * xyn(i)
	  y = y + ys(i) * xyn(i)
	end do

	end

c*************************************************************

	subroutine exbil

c exercise bi-linear interpolation

	implicit none

	real u(4)
	real x,y
	real r
	real bilint

	write(6,*) 'Enter values at 4 vertices :'
	write(6,*) '(0,0)  (1,0)  (0,1)  (1,1)'
	read(5,*) u
	write(6,*) u

	x=0.
	y=0.

	do while( x .ge. 0. .and. y .ge. 0. )
	  write(6,*) 'Enter x,y [0..1] (negative to end) :'
	  read(5,*) x,y
	  r = bilint(u,x,y)
	  write(6,*) x,y,r
	end do

	end
	  
c*************************************************************

	subroutine extri

c exercise tri-linear interpolation

	implicit none

	integer ndim
	parameter ( ndim = 30 )

	logical berror
	integer n
	integer i,j
	integer ix,iy,iz
	real u(8)
	real uvw(0:2)
	real x,y,z
	real val,compare
	real triint,grand

	real xyz(8,3)
	data xyz /
     +			 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 ,
     +			 0 , 0 , 1 , 1 , 0 , 0 , 1 , 1 ,
     +			 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1  
     +		 /

c interpolation property

	berror = .false.
	write(6,*) 'Checking interpolation property...'

	do i=1,8
	  u(i) = 0.
	end do

	do i=1,8
	  u(i) = 1.
	  do j=1,8
	    x = xyz(j,1)
	    y = xyz(j,2)
	    z = xyz(j,3)
	    val = triint(u,x,y,z)
	    compare = 0.
	    if( i .eq. j ) compare = 1.
	    if( val .ne. compare ) then
		write(6,*) i,x,y,z,val
		berror = .true.
	    end if
	  end do
	  u(i) = 0.
	end do

	if( berror ) then
	  write(6,*) 'There have been errors in check...'
	  stop 'error stop'
	else
	  write(6,*) '   ...passed'
	end if

c linear interpolation

	do j=0,2

	do i=1,8
	  if( xyz(i,j) .ne. 0. ) u(i) = 1.
	end do

	uvw(j) = 0.5
	do n=1,ndim
	  uvw( mod(j+1,3) ) = grand()
	  uvw( mod(j+2,3) ) = grand()
	  x = uvw(0)
	  y = uvw(1)
	  z = uvw(2)
	  val = triint(u,x,y,z)
	  write(6,*) n,x,y,z,val
	end do

	end do

	end

c*************************************************************

	function grand()

c random numbers

	implicit none

	real grand
	real rand

	integer iseed
	save iseed
	data iseed /37582629/

	grand = rand(iseed)

	end 

c*************************************************************
c
c uncomment next lines to run test routines

c	program hptri
c	call extri
c	call exbil
c	end

c*************************************************************

