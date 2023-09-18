!
! routines for reading time series
!
! revision log :
!
! 20.10.2014    ggu     integrating datetime into time series
! 10.02.2015    ggu     length of line set to 2048
!
! notes :
!
! keyword example setting date:		"#date: 20071001 0"
!
!*************************************************************
!-------------------------------------------------------------
        module tsfile
!-------------------------------------------------------------
        contains
!-------------------------------------------------------------

	subroutine ts_get_file_info(file,nvar)

! get info on time series file (TS)
!
! nvar is returned as the number of available columns (except ttime)
! nvar <= 0 for error or no TS file

        use fil
        use convert
        use utility
        use timing
        use shympi

	implicit none

	character*(*) file	!file name
	integer nvar		!variables (columns) in file (except time)

	character*2048 line
	integer iunit,i,iend,nrec,iline,nvar0
	double precision d(1),time1

!------------------------------------------------------
! find unit to open file
!------------------------------------------------------

        if(ln_timing) time1 = shympi_wtime()

        iunit = 90
        call find_unit(iunit)

!------------------------------------------------------
! open file
!------------------------------------------------------

	nvar = 0	!this means no such file etc..
	nrec = 0
	iline = 0
	nvar0 = 0

        open(iunit,file=file,form='formatted',status='old',err=2)

!------------------------------------------------------
! we try to read 3 records
!------------------------------------------------------

	do while( iline < 3 )
	  nrec = nrec + 1
	  read(iunit,'(a)',err=3,end=2) line
	  i = ichafs(line)
	  iend = len_trim(line)
	  if( iend > 2000 ) goto 95
	  if( i > 0 .and. line(i:i) /= '#' ) then	!valid record
	    iline = iline + 1
	    nvar = iscand(line,d,0)		!count values on line
	    if( nvar < 0 ) nvar = -nvar-1	!read error in number -nvar
	    nvar = nvar - 1			!do not count time column
	    if( nvar <= 0 ) exit		!no data found
	    if( nvar0 == 0 ) nvar0 = nvar
	    if( nvar /= nvar0 ) exit		!varying number of data
	  end if
	end do

	if( nvar /= nvar0 ) nvar = 0

        if(ln_timing) io_time = io_time + shympi_wtime() - time1

	!write(6,*) 'ts_get_file_info: ',iline,nrec,nvar,nvar0

!------------------------------------------------------
! close file
!------------------------------------------------------

    2	continue
	close(iunit)

!------------------------------------------------------
! end of routine
!------------------------------------------------------

	return
   95	continue
	write(6,*) 'read error in line ',nrec,' of file ',file
	write(6,*) 'last character read: ',iend
	stop 'error stop ts_get_file_info: line too long'
    3	continue
	write(6,*) 'read error in line ',nrec,' of file ',file
	stop 'error stop ts_get_file_info: read error'
	end

!*************************************************************

	subroutine ts_open_file(file,nvar,datetime,iunit)

        use fil
        use convert
        use utility
        use timing
        use shympi

	implicit none

	character*(*) file	!file name
	integer nvar		!variables (columns) in file (except time)
	integer datetime(2)
	integer iunit

	character*2048 line
	character*10 key
	character*60 info	!still have to use this data
	integer i,j,nrec,ioff,iend
	integer date,time
	double precision d(3),time1

!------------------------------------------------------
! open file
!------------------------------------------------------

        if(ln_timing) time1 = shympi_wtime()

	nrec = 0
	nvar = 0
	date = 0
	time = 0
	datetime = 0
	info = ' '

	iunit = ifileo(0,file,'formatted','old')

	if( iunit <= 0 ) return

!------------------------------------------------------
! check if header exists and get nvar
!------------------------------------------------------

	do while(.true.)
	  nrec = nrec + 1
	  read(iunit,'(a)',err=3,end=3) line
	  if( ichafs(line) <= 0 ) cycle		!empty line
	  if( ts_has_keyword(line) ) then
	    call ts_parse_keyword(iunit,line,datetime,info)
	  else
	    exit
	  end if
	  !if( info .ne. ' ' ) write(6,*) 'info line: ',info
	end do

!------------------------------------------------------
! set date, nvar and backspace file
!------------------------------------------------------

	iend = len_trim(line)
	if( iend > 2000 ) goto 95
	nvar = iscand(line,d,0)		!count values on line
	if( nvar < 0 ) nvar = -nvar-1	!read error in number -nvar
	nvar = nvar - 1			!do not count time column

	backspace(iunit)

        if(ln_timing) io_time = io_time + shympi_wtime() - time1

!------------------------------------------------------
! end of routine
!------------------------------------------------------

	return
   95	continue
	write(6,*) 'read error in line ',nrec,' of file ',file
	write(6,*) 'last character read: ',iend
	stop 'error stop ts_open_file: line too long'
    3	continue
	write(6,*) 'read error in line ',nrec,' of file ',file
	stop 'error stop ts_open_file: read error'
	end

!*************************************************************

	subroutine ts_peek_next_record(iunit,nvar,it,f,datetime,ierr)

! peeks into one record of time series file

	implicit none

	integer iunit
	integer nvar		!variables (columns) in file (except time)
	double precision it
	double precision f(nvar)
	integer datetime(2)
	integer ierr

	call ts_read_next_record(iunit,nvar,it,f,datetime,ierr)

	backspace(iunit)

	end

!*************************************************************

	subroutine ts_read_next_record(iunit,nvar,it,f,datetime,ierr)

! reads one record of time series file

        use timing
        use shympi

	implicit none

	integer iunit
	integer nvar		!variables (columns) in file (except time)
	double precision it
	double precision f(nvar)
	integer datetime(2)
	integer ierr

	integer i
	double precision t,time1
	character*132 line
	character*80 info

	datetime = 0
	info = ' '

!------------------------------------------------------
! check if file is open
!------------------------------------------------------

	if( iunit <= 0 ) then
	  ierr = 1
	  return
	end if

!------------------------------------------------------
! read record
!------------------------------------------------------

    1	continue

        if(ln_timing) time1 = shympi_wtime()

	read(iunit,*,end=2,err=3) t,(f(i),i=1,nvar)
	it = t
	ierr = 0
	return

        if(ln_timing) io_time = io_time + shympi_wtime() - time1

!------------------------------------------------------
! error handling
!------------------------------------------------------

    2	continue
	ierr = -1
	return

    3	continue
	backspace(iunit)	!see if it was a keyword line
	read(iunit,'(a)',end=2,err=4) line
	
	if( ts_has_keyword(line) ) then
	  call ts_parse_keyword(iunit,line,datetime,info)
	  !if( info .ne. ' ' ) write(6,*) 'info line: ',info
	  goto 1		!read next line
	end if

    4	continue
	ierr = 3
	return

!------------------------------------------------------
! end of routine
!------------------------------------------------------

	end

!*************************************************************
!*************************************************************
!*************************************************************

	function ts_has_keyword(line)

        use utility
	implicit none

	logical ts_has_keyword
	character*(*) line

	integer i

	ts_has_keyword = .false.

	i = ichafs(line)
	if( i <= 0 ) return
	if( line(i:i) == '#' ) ts_has_keyword = .true.

	end

!*************************************************************

	subroutine ts_get_keyword(line,key,ioff)

! gets keyword from line - rest of line after ioff
!
! keyword looks like:   "#key:"
!
! example:		"#date: 20071001 0"

        use utility
	implicit none

	character*(*) line
	character*(*) key
	integer ioff

	integer i,j


	key = ' '
	ioff = 1

	i = ichafs(line)
	if( i <= 0 ) return
	if( line(i:i) /= '#' ) return	!no keyword

	do j=i+1,len(line)
	  if( line(j:j) == ' ' ) return	!no white space allowed
	  if( line(j:j) == ':' ) then	!end of keyword found
	    if( j-i-1 <= 0 ) return	!no keyword present
	    key = line(i+1:j-1)
	    ioff = j+1
	    return
	  end if
	end do

	return
	end

!*************************************************************

	subroutine ts_parse_keyword(iunit,line,datetime,info)

        use fil

	implicit none

	integer iunit
	character*(*) line
	integer datetime(2)
	character*(*) info

	character*10 key
	character*75 file
	integer ioff

	call ts_get_keyword(line,key,ioff)

	if( key == 'date' ) then		!date
	  call ts_parse_datetime(line(ioff:),datetime)
	else if( key == 'info' ) then		!info
	  call ts_parse_info(line(ioff:),info)
	else if( key == ' ' ) then		!nothing
	else
	  call filna(iunit,file)
	  write(6,*) 'not recognized keyword: ',key
	  write(6,*) 'file open at unit: ',iunit
	  write(6,*) 'file name: ',file
	  stop 'error stop ts_parse_keyword: not recognized keyword'
	end if

	end

!*************************************************************

	subroutine ts_parse_datetime(line,datetime)

        use convert

	implicit none

	character*(*) line
	integer datetime(2)

	integer j
	double precision d(3)

	datetime = 0

	j = iscand(line,d,2)

	if( j >= 1 ) datetime(1) = nint(d(1))
	if( j >= 2 ) datetime(2) = nint(d(2))

	if( datetime(1) > 0 .and. datetime(1) < 10000 ) then
	  datetime(1) = 10000*datetime(1) + 101
	end if

	end

!*************************************************************

	subroutine ts_parse_info(line,info)

        use utility
	implicit none

	character*(*) line
	character*(*) info

	integer i

	info = ' '

	i = ichafs(line)
	if( i <= 0 ) return

	info = line(i:)

	end

!*************************************************************
!*************************************************************
!*************************************************************

	subroutine ts_test

	implicit none

	integer nvar
	integer iunit
	integer datetime(2)
	integer ierr
	integer nrec,i
	double precision it
	double precision f(10)
	character*60 file
	character*20 line

	line = ' '
	nrec = 0
	file = 'ts_test.txt'

	call ts_get_file_info(file,nvar)
	write(6,*) 'file info: nvar = ',nvar

	call ts_open_file(file,nvar,datetime,iunit)
	write(6,*) 'file open: nvar = ',nvar
	write(6,*) 'file open: iunit = ',iunit
	write(6,*) 'file open: datetime = ',datetime
	if( datetime(1) > 0 ) then
	  !call dtsini(datetime(1),datetime(2))
	end if

	do
	  call ts_read_next_record(iunit,nvar,it,f,datetime,ierr)
	  if( ierr .ne. 0 ) exit
	  if( datetime(1) > 0 ) then
	    write(6,*) 'datetime: ',datetime
	    !call dtsini(datetime(1),datetime(2))
	  end if
	  nrec = nrec + 1
	  !call dtsgf(nint(it),line)
	  !write(6,*) it,(f(i),i=1,nvar)
	  if( mod(nrec,1) .eq. 0 ) write(6,*) it,f(1),line
	end do

	write(6,*) nrec,' records read'

	end

!*************************************************************
!	program ts_test_main
!	call ts_test
!	end
!*************************************************************
!-------------------------------------------------------------
        end module tsfile
!-------------------------------------------------------------
