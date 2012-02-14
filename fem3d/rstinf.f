
c info on restart file

	program rstinf

	implicit none

	integer iunit,it,nvers,nrec,nkn,nel,nlv,ierr
	integer nread
	character*60 file
	character*40 title1
	character*30 title2

	nread = 0
	iunit = 1
	file = ' '

	write(6,*) 'Enter file name: '
	read(5,'(a)') file
	if( file .eq. ' ' ) stop

	open(iunit,file=file,status='old',form='unformatted')

	title1 = '      irec      time   version      nrec'
	title2 = '       nkn       nel       nlv'
	write(6,1000) title1,title2
    1	continue

	call skip_rst(iunit,it,nvers,nrec,nkn,nel,nlv,ierr)
	if( ierr .ne. 0 ) goto 2
	nread = nread + 1
	write(6,'(7i10)') nread,it,nvers,nrec,nkn,nel,nlv

	goto 1

    2	continue
	
	write(6,1000) title1,title2
	write(6,*)
	write(6,*) 'Number of records read: ',nread

	stop
 1000	format(a40,a30)
	end
