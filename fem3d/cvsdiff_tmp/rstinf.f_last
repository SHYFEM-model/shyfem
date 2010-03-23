
c info on restart file

	implicit none

	integer iunit,it,nvers,nrec,nkn,nel,nlv,ierr
	integer nread
	character*60 file

	nread = 0
	iunit = 1
	file = ' '

	write(6,*) 'Enter file name: '
	read(5,'(a)') file
	if( file .eq. ' ' ) stop

	open(iunit,file=file,status='old',form='unformatted')

    1	continue

	call skip_rst(iunit,it,nvers,nrec,nkn,nel,nlv,ierr)
	if( ierr .ne. 0 ) goto 2
	nread = nread + 1
	write(6,'(7i10)') nread,it,nvers,nrec,nkn,nel,nlv

	goto 1

    2	continue
	
	write(6,*) 'Number of records read: ',nread

	end
