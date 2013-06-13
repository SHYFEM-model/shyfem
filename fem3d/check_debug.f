
	program check_debug

c checks two files written with check_debug from ht

	implicit none

	integer ndim
	parameter (ndim=2000000)

	character*60 name_on,name_off
	logical bcheck
	integer it1,it2,it
	integer nt1,nt2,nt
	integer nf1,nf2,nf
	integer nrec
	integer i,idiff

	real val1(ndim)
	real val2(ndim)

	name_on = 'debug_on.dat'
	name_off = 'debug_off.dat'

	open(1,file=name_on,status='old',form='unformatted')
	open(2,file=name_off,status='old',form='unformatted')

	do while(.true.)

	  read(1,end=9) it1
	  read(2,end=9) it2
	  if( it1 .ne. it2 ) goto 99
	  it = it1
	  write(6,*) 'time = ',it

	  nrec = 0
	  do while(.true.)
	    read(1) nt1,nf1
	    read(2) nt2,nf2
	    nrec = nrec + 1
	    if( nt1 .ne. nt2 ) goto 98
	    if( nf1 .ne. nf2 ) goto 98
	    nt = nt1
	    nf = nf1

	    if( nt .eq. 0 ) exit
	    if( nt .gt. ndim ) goto 97

	    read(1) (val1(i),i=1,nt)
	    read(2) (val2(i),i=1,nt)

	    idiff = 0
	    bcheck = .true.
	    !bcheck = it .ge. 87300
	    !bcheck = bcheck .and. nrec .ne. 5
	    if( bcheck ) then
	      call check_val(it,nrec,nt,nf,val1,val2,idiff)
	      write(6,*) nrec,nt,nf,idiff
	    end if

	  end do

	end do

    9	continue

	close(1)
	close(2)

	stop
   99	continue
	write(6,*) it1,it2
	stop 'error stop check_debug: time mismatch'
   98	continue
	write(6,*) nt1,nt2,nf1,nf2
	stop 'error stop check_debug: size mismatch'
   97	continue
	write(6,*) it,nrec,nt,ndim
	stop 'error stop check_debug: dimension'
	end

c*******************************************************************

	subroutine check_val(it,nrec,nt,nf,val1,val2,idiff)

	implicit none

	integer it,nrec
	integer nt,nf,idiff
	real val1(nt)
	real val2(nt)

	integer i,k,l

	idiff = 0

	do i=1,nt
	  if( val1(i) .ne. val2(i) ) then
	    k = 1 + (i-1)/nf
	    l = 1 + mod(i-1,nf)
	    !write(77,*) it,nrec,k,l,val1(i),val2(i)
	    idiff = idiff + 1
	  end if
	end do

	end

c*******************************************************************

