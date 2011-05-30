
	program clean_reg_input

c concatenates regular input files into one
c
c files must have been cleaned with clean_reg_input (time records)

	implicit none

	integer ndim,nparam
	parameter (ndim=1000)
	parameter (nparam=5)

	real array(ndim,nparam)
	character*50 cparam(nparam)

	integer i
	integer it,n,np,itold,idtold,nx,ny
	integer iunit
	integer itold
	real amin,amax,amed
	real geo(5)

	itold = -1

    1	continue
	  call read_record(ndim,nparam,array,cparam,it,n,np,nx,ny,geo)
	  if( np .le. 0 ) goto 2
	  !write(6,1000) it,(cparam(i),i=1,np)
	  write(6,*) it
	  if( itold .eq. -1 .or. it .gt. itold ) then
	    call write_record(ndim,nparam,array,cparam,it,n,np,nx,ny,geo)
	  else
	    write(6,*) 'eliminating record: ',it,itold,it-itold
	  end if
	  itold = it
	goto 1
    2	continue

	write(6,*) 'output written to file unit 66'

	stop
 1000	format(i12,5(a15,2x))
	end

c******************************************************************

	subroutine write_record(ndim,nparam,array,cparam,it
     +				,n,np,nx,ny,geo)

	implicit none

	integer ndim,nparam
	integer it,n,np,nx,ny
	real array(ndim,nparam)
	character*(*) cparam(nparam)
	real geo(5)

	integer i,j

	n = nx*ny

	write(66,*) it,np,nx,ny,(geo(i),i=1,5)

	do i=1,np
	  write(66,'(a)') cparam(i)
	  write(66,*) (array(j,i),j=1,n)
	end do

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

