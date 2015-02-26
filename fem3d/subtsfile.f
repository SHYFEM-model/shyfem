c
c routines for reading time series
c
c revision log :
c
c 20.10.2014    ggu     integrating datetime into time series
c 10.02.2015    ggu     length of line set to 2048
c
c notes :
c
c keyword example setting date:		"#date: 20071001 0"
c
c*************************************************************

	subroutine ts_get_file_info(file,nvar)

c get info on time series file (TS)
c
c nvar is returned as the number of available columns (except ttime)
c nvar <= 0 for error or no TS file

	implicit none

	character*(*) file	!file name
	integer nvar		!variables (columns) in file (except time)

	character*2048 line
	integer iunit,i,iend,nrec,iline,nvar0
	double precision d(1)

	integer iscand,ichafs

c------------------------------------------------------
c find unit to open file
c------------------------------------------------------

        iunit = 90
        call find_unit(iunit)

c------------------------------------------------------
c open file
c------------------------------------------------------

	nvar = 0	!this means no such file etc..
	nrec = 0
	iline = 0
	nvar0 = 0

        open(iunit,file=file,form='formatted',status='old',err=2)

c------------------------------------------------------
c we try to read 3 records
c------------------------------------------------------

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

	!write(6,*) 'ts_get_file_info: ',iline,nrec,nvar,nvar0

c------------------------------------------------------
c close file
c------------------------------------------------------

    2	continue
	close(iunit)

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	return
   95	continue
	write(6,*) 'read error in line ',nrec,' of file ',file
	write(6,*) 'last character read: ',iend
	stop 'error stop ts_get_file_info: line too long'
    3	continue
	write(6,*) 'read error in line ',nrec,' of file ',file
	stop 'error stop ts_get_file_info: read error'
	end

c*************************************************************

	subroutine ts_open_file(file,nvar,datetime,iunit)

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
	double precision d(3)

	integer iscand,ichafs
	integer ifileo
	logical ts_has_keyword

c------------------------------------------------------
c open file
c------------------------------------------------------

	nrec = 0
	nvar = 0
	date = 0
	time = 0
	datetime = 0
	info = ' '

	iunit = ifileo(0,file,'formatted','old')

	if( iunit <= 0 ) return

c------------------------------------------------------
c check if header exists and get nvar
c------------------------------------------------------

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

c------------------------------------------------------
c set date, nvar and backspace file
c------------------------------------------------------

	iend = len_trim(line)
	if( iend > 2000 ) goto 95
	nvar = iscand(line,d,0)		!count values on line
	if( nvar < 0 ) nvar = -nvar-1	!read error in number -nvar
	nvar = nvar - 1			!do not count time column

	backspace(iunit)

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	return
   95	continue
	write(6,*) 'read error in line ',nrec,' of file ',file
	write(6,*) 'last character read: ',iend
	stop 'error stop ts_open_file: line too long'
    3	continue
	write(6,*) 'read error in line ',nrec,' of file ',file
	stop 'error stop ts_open_file: read error'
	end

c*************************************************************

	subroutine ts_peek_next_record(iunit,nvar,it,f,datetime,ierr)

c peeks into one record of time series file

	implicit none

	integer iunit
	integer nvar		!variables (columns) in file (except time)
	double precision it
	real f(nvar)
	integer datetime(2)
	integer ierr

	call ts_read_next_record(iunit,nvar,it,f,datetime,ierr)

	backspace(iunit)

	end

c*************************************************************

	subroutine ts_read_next_record(iunit,nvar,it,f,datetime,ierr)

c reads one record of time series file

	implicit none

	integer iunit
	integer nvar		!variables (columns) in file (except time)
	double precision it
	real f(nvar)
	integer datetime(2)
	integer ierr

	integer i
	double precision t
	character*132 line
	character*80 info

	logical ts_has_keyword

	datetime = 0
	info = ' '

c------------------------------------------------------
c check if file is open
c------------------------------------------------------

	if( iunit <= 0 ) then
	  ierr = 1
	  return
	end if

c------------------------------------------------------
c read record
c------------------------------------------------------

    1	continue

	read(iunit,*,end=2,err=3) t,(f(i),i=1,nvar)
	it = t
	ierr = 0
	return

c------------------------------------------------------
c error handling
c------------------------------------------------------

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

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	end

c*************************************************************
c*************************************************************
c*************************************************************

	function ts_has_keyword(line)

	implicit none

	logical ts_has_keyword
	character*(*) line

	integer i
	integer ichafs

	ts_has_keyword = .false.

	i = ichafs(line)
	if( i <= 0 ) return
	if( line(i:i) == '#' ) ts_has_keyword = .true.

	end

c*************************************************************

	subroutine ts_get_keyword(line,key,ioff)

c gets keyword from line - rest of line after ioff
c
c keyword looks like:   "#key:"
c
c example:		"#date: 20071001 0"

	implicit none

	character*(*) line
	character*(*) key
	integer ioff

	integer i,j

	integer ichafs

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

c*************************************************************

	subroutine ts_parse_keyword(iunit,line,datetime,info)

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

c*************************************************************

	subroutine ts_parse_datetime(line,datetime)

	implicit none

	character*(*) line
	integer datetime(2)

	integer j
	double precision d(3)

	integer iscand

	datetime = 0

	j = iscand(line,d,2)

	if( j >= 1 ) datetime(1) = nint(d(1))
	if( j >= 2 ) datetime(2) = nint(d(2))

	if( datetime(1) > 0 .and. datetime(1) < 10000 ) then
	  datetime(1) = 10000*datetime(1) + 101
	end if

	end

c*************************************************************

	subroutine ts_parse_info(line,info)

	implicit none

	character*(*) line
	character*(*) info

	integer i
	integer ichafs

	info = ' '

	i = ichafs(line)
	if( i <= 0 ) return

	info = line(i:)

	end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine ts_test

	implicit none

	integer nvar
	integer iunit
	integer datetime(2)
	integer ierr
	integer nrec,i
	double precision it
	real f(10)
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

c*************************************************************
c	program ts_test_main
c	call ts_test
c	end
c*************************************************************
