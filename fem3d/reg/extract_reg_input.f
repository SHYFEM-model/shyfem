
	program extract_reg_input

c extracts node from regular input file

	implicit none

	integer ndim,nparam
	parameter (ndim=1000)
	parameter (nparam=5)

	real array(ndim,nparam)
	character*50 cparam(nparam)

	integer i
	integer it,n,np,itold,idtold,nx,ny
	integer iunit
	integer ix,iy
	real amin,amax,amed
	real rx,ry
	real geo(5)
	real vals(nparam)

	itold = -1
	idtold = -1
	ix = -1
	iy = -1
	rx = 19.31
	ry = 42.17

    1	continue
	  call read_record(ndim,nparam,array,cparam,it,n,np,nx,ny,geo)
	  if( ix .lt. 0 ) call get_geo(geo,rx,ry,ix,iy)
	  if( np .le. 0 ) goto 2
	  call get_val(ndim,np,nx,ny,ix,iy,array,vals)
	  write(6,*) it,(vals(i),i=1,np)
	goto 1
    2	continue

	stop
 1000	format(i12,5(a15,2x))
	end

c******************************************************************

	subroutine read_record(ndim,nparam,array,cparam,it
     +				,n,np,nx,ny,geo)

	implicit none

	integer ndim,nparam
	integer it,n,np,nx,ny
	real array(ndim,nparam)
	character*(*) cparam(nparam)
	real geo(5)

	integer i,j

	n = 0
	np = 0

	read(5,*,end=1) it,np,nx,ny,(geo(i),i=1,5)
	!write(6,*) it,np,nx,ny
	n = nx*ny
	if( np .gt. nparam ) goto 99
	if( n .gt. ndim ) goto 99

	do i=1,np
	  read(5,'(a)') cparam(i)
	  !write(6,*) cparam(i)
	  read(5,*) (array(j,i),j=1,n)
	end do

    1	continue

	return
   99	continue
	write(6,*) ndim,nparam,n,np
	stop 'error stop read_record: dimensions'
	end

c******************************************************************

	subroutine get_val(ndim,np,nx,ny,ix,iy,array,vals)

	implicit none

	integer ndim,np,nx,ny,ix,iy
	real array(ndim,np)
	real vals(np)

	integer ip,ia

	do ip=1,np
	  ia = nx*(iy-1) + ix
	  vals(ip) = array(ia,ip)
	end do

	end

c******************************************************************
	
	subroutine get_geo(geo,rx,ry,ix,iy)

	implicit none

	real geo(5)
	real rx,ry
	integer ix,iy

	real x0,y0,dx,dy

	x0 = geo(1)
	y0 = geo(2)
	dx = geo(3)
	dy = geo(4)

	ix = 1 + nint((rx-x0)/dx)
	iy = 1 + nint((ry-y0)/dy)

	write(6,*) 'coordinates: ',ix,iy

	end

c******************************************************************

