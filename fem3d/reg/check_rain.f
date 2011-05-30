
	program check_rain

	implicit none

	integer it,itold,idt
	real amin,amax,amed
	double precision aamin,aamax,aamed

	aamin = 0.
	aamed = 0.
	aamax = 0.

	read(5,*) it,amin,amed,amax

    1	continue
	  itold = it
	  read(5,*,end=2) it,amin,amed,amax
	  idt = it-itold
	  aamin = aamin + idt * amin
	  aamed = aamed + idt * amed
	  aamax = aamax + idt * amax
	goto 1
    2	continue

	aamin = aamin / 86400
	aamed = aamed / 86400
	aamax = aamax / 86400

	write(6,*) aamin,aamed,aamax

	end

c******************************************************************

	subroutine read_record(ndim,nparam,array,cparam,it,n,np)

	implicit none

	integer ndim,nparam
	integer it,n,np
	real array(ndim,nparam)
	character*(*) cparam(nparam)

	integer nx,ny,i,j
	real x0,y0,dx,dy,flag

	n = 0
	np = 0

	read(5,*,end=1) it,np,nx,ny,x0,y0,dx,dy,flag
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

	subroutine get_stats(n,array,amin,amax,amed)

	implicit none

	integer n
	real array(n)
	real amin,amax,amed

	integer j
	real a

	amin = array(1)
	amax = array(1)
	amed = 0.

	do j=1,n
	  a = array(j)
	  amin = min(amin,a)
	  amax = max(amax,a)
	  amed = amed + a
	end do

	amed = amed / n

	end

c******************************************************************

