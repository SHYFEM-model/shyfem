
	program check_reg_input

c checks regular input file

	implicit none

	integer ndim,nparam
	parameter (ndim=1000)
	parameter (nparam=5)

	real array(ndim,nparam)
	character*50 cparam(nparam)

	integer i
	integer it,n,np
	integer iunit
	integer itold
	real amin,amax,amed

	itold = -1

    1	continue
	  call read_record(ndim,nparam,array,cparam,it,n,np)
	  if( np .le. 0 ) goto 2
	  if( itold .ne. -1 .and. it .le. itold ) goto 99
	  write(6,1000) it,(cparam(i),i=1,np)
	  do i=1,np
	    call get_stats(n,array(1,i),amin,amax,amed)
	    iunit = 100 + i
	    write(iunit,*) it,amin,amed,amax
	  end do
	  itold = it
	goto 1
    2	continue

	stop
   99	continue
	write(6,*) 'itold,it,it-itold: ',itold,it,it-itold
	stop 'error stop check_reg_input: it'
 1000	format(i12,5(a15,2x))
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

