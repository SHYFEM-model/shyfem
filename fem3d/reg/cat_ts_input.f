
	program cat_reg_input

c concatenates regular input files into one

	implicit none

	integer ndim,nparam
	parameter (ndim=1000)
	parameter (nparam=5)

	real array(ndim,nparam)
	character*50 cparam(nparam)

	integer i
	integer it,n,np,itold,idtold,nx,ny
	integer iunit
	real amin,amax,amed,val
	real geo(5)

	itold = -1
	idtold = -1

    1	continue
	  read(5,*,end=2) it,val
	  call adjust_it(it)
	  write(6,*) it,val
	  if( it .ge. 0 ) then
	    write(66,*) it,val
	  end if
	goto 1
    2	continue

	stop
 1000	format(i12,5(a15,2x))
	end

c******************************************************************

	subroutine adjust_it(it)

c adjusts time

	implicit none

	integer it

	integer idt
	integer itold,idtold,itadd
	integer it0,itorg
	save itold,idtold,itadd
	data itold,idtold,itadd /-1,-1,0/

	it0 = 86400

	if( itold .eq. -1 ) then
	  itold = it
	  return
	end if

	if( idtold .eq. -1 ) then
	  idtold = it - itold
	  itold = it
	  return
	end if

	itorg = it
	it = itorg + itadd

	if( it .lt. itold ) then
	  itadd = itold + idtold - itorg
	  write(6,*) 'time adjusted... ',it,itold,itadd,itorg+itadd
	end if

	it = itorg + itadd
	idt = it - itold
	if( idt .le. 0 ) then	!adjust
	  write(6,*) 'non-positive time step:  ',it,itold,idt
	  it = -it
	  return
	end if

	idt = it - itold
	if( idt .ne. idtold ) then
	  write(6,*) 'different time step:  ',it,itold,idt,idtold
	end if

	itold = it

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

