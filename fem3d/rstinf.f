
c info on restart file

	program rstinf

	implicit none

	integer iunit,it,nvers,nrec,nkn,nel,nlv,iflag,ierr
	integer nread
	character*60 file
	character*35 title1
	character*40 title2

	nread = 0
	iunit = 1
	file = ' '

	write(6,*) 'Enter file name: '
	read(5,'(a)') file
	if( file .eq. ' ' ) stop

	open(iunit,file=file,status='old',form='unformatted')

	title1 = '      irec        time version nrec'
	title2 = '       nkn       nel       nlv     iflag'
c                 1234567890123456789012345678901234567890
	write(6,1000) title1,title2
    1	continue

	call skip_rst(iunit,it,nvers,nrec,nkn,nel,nlv,iflag,ierr)
	if( ierr .ne. 0 ) goto 2
	nread = nread + 1
	write(6,1010) nread,it,nvers,nrec,nkn,nel,nlv,iflag

	goto 1

    2	continue
	
	write(6,1000) title1,title2
	write(6,*)
	write(6,*) 'Number of records read: ',nread
	write(6,*)
	write(6,*) 'Meaning of iflag:'
	write(6,*) '         1          hm3v'
	write(6,*) '        10          ibarcl (T/S/rho)'
	write(6,*) '       100          iconz (cnv)'
	write(6,*) '      1000          iwvert (wlnv)'
	write(6,*) '     10000          ieco (ecological variables)'

	stop
 1000	format(a35,a40)
 1010	format(i10,i12,i8,i5,4i10)
	end
