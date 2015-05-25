
c info on restart file

	program rstinf

	implicit none

	integer iunit,it,nvers,nrec,nkn,nel,nlv,iflag,ierr
	integer nread
	double precision atime
	double precision atime_anf
	double precision atime_end
	character*20 line
	character*60 file
	character*52 title1
	character*48 title2

	nread = 0
	iunit = 1
	file = ' '

	write(6,*) 'Enter file name: '
	read(5,'(a)') file
	if( file .eq. ' ' ) stop

	open(iunit,file=file,status='old',form='unformatted')

	title1 = 'version nrec       nkn       nel       nlv     iflag'
c                 12345678901234567890123456789012345678901234567890123
	title2 = '   irec               atime                 date'

    1	continue

	call skip_rst(iunit,atime,it,nvers,nrec,nkn,nel,nlv,iflag,ierr)
	if( ierr .ne. 0 ) goto 2

	if( nread == 0 ) then
	  write(6,1000) title1
	  write(6,1010) nvers,nrec,nkn,nel,nlv,iflag
	  write(6,1001) title2
	  atime_anf = atime
	end if

	nread = nread + 1
	call dts_format_abs_time(atime,line)
	write(6,1011) nread,atime,it,line
	atime_end = atime

	goto 1

    2	continue
	
	write(6,1001) title2
	write(6,*)
	write(6,1001) title1
	write(6,1010) nvers,nrec,nkn,nel,nlv,iflag
	write(6,*)
	write(6,*) 'Number of records read: ',nread
	call dts_format_abs_time(atime_anf,line)
	write(6,*) 'Initial time in file: ',atime_anf,line
	call dts_format_abs_time(atime_end,line)
	write(6,*) 'Final time in file: ',atime_end,line
	write(6,*)
	write(6,*) 'Meaning of iflag:'
	write(6,*) '         1          hm3v'
	write(6,*) '        10          ibarcl (T/S/rho)'
	write(6,*) '       100          iconz (cnv)'
	write(6,*) '      1000          iwvert (wlnv)'
	write(6,*) '     10000          ieco (ecological variables)'

	stop
 1000	format(a52)
 1001	format(a48)
 1010	format(i7,i5,4i10)
 1011	format(i7,d20.1,i20,1x,a20)
	end
