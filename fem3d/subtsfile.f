c
c routines for reading time series
c
c revision log :
c
c 20.10.2014    ggu     integrating datetime into time series
c 10.02.2015    ggu     length of line set to 2048
c 15.05.2017    ggu     new version to also read time string
c
c notes :
c
c keyword example setting date:		"#date: 20071001 0"
c time column as string:		"2007-10-01::00:00:00"
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

	integer iunit,i,ierr
	integer datetime(2)
	double precision dtime
	real, allocatable :: f(:)

	call ts_open_file(file,nvar,datetime,iunit)
	!write(6,*) 'ggguuu: ',trim(file),nvar,datetime,iunit
	if( nvar <= 0 ) return
	if( iunit <= 0 ) return

	allocate(f(nvar))

	do i=1,3
	  call ts_read_next_record(iunit,nvar,dtime,f,datetime,ierr)
	  if( ierr /= 0 ) exit
	end do

	if( ierr > 0 ) nvar = 0

	close(iunit)

	end

c*************************************************************

	subroutine ts_get_extra_time(file,dtime,datetime)

	implicit none

	character*(*) file	!file name
	double precision dtime
	integer datetime(2)

	integer nvar		!variables (columns) in file (except time)
	integer iunit
	integer ios,i,is,ierr
	character*2048 line,dummy

	integer istot,istod

	dtime = 0.
	datetime = 0

	call ts_open_file(file,nvar,datetime,iunit)
	if( nvar <= 0 .or. iunit <= 0 ) return

	line = ' '
	do while( line == ' ' )
	  read(iunit,'(a)',iostat=ios) line
	  if( ios /= 0 ) return
	  if( line == ' ' ) line = ' '
	  if( line(1:1) == '#' ) line = ' '
	end do

	is = 1
	ios = istod(line,dtime,is)		!read time column
	if( ios /= 1 ) goto 1

	do i=1,nvar
	  ios = istot(line,dummy,is)
	  if( ios /= 1 ) goto 1
	end do

	ios = istot(line,dummy,is)
	if( ios == 1 ) then
	  call string2datetime(dummy(1:20),datetime,ierr)
	end if

    1	continue
	close(iunit)

	end

c*************************************************************

	subroutine ts_open_file(file,nvar,datetime,iunit)

	implicit none

	character*(*) file	!file name
	integer nvar		!variables (columns) in file (except time)
	integer datetime(2)
	integer iunit

	integer ierr
	real f(nvar)
	double precision dtime

	integer ifileo

c------------------------------------------------------
c open file
c------------------------------------------------------

	nvar = 0
	datetime = 0

	iunit = ifileo(0,file,'formatted','old')

	if( iunit <= 0 ) return

	call ts_read_next_record(iunit,nvar,dtime,f,datetime,ierr)

	if( ierr == 0 ) then
	  backspace(iunit)
	else
	  close(iunit)
	  nvar = -1
	  iunit = 0
	end if

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	end

c*************************************************************

	subroutine ts_peek_next_record(iunit,nvar,dtime,f,datetime,ierr)

c peeks into one record of time series file

	implicit none

	integer iunit
	integer nvar		!variables (columns) in file (except time)
	double precision dtime
	real f(nvar)
	integer datetime(2)
	integer ierr

	call ts_read_next_record(iunit,nvar,dtime,f,datetime,ierr)

	backspace(iunit)

	end

c*************************************************************

	subroutine ts_read_next_record(iunit,nvar,dtime,f,datetime,ierr)

! reads one record of time series file

	implicit none

	integer iunit
	integer nvar		!variables (columns) in file (except time)
	double precision dtime
	real f(nvar)
	integer datetime(2)
	integer ierr

	integer i,iend,ios,is,nvar0
	double precision t
	character*2048 line
	character*80 info,stime

	logical ts_has_keyword
	integer ichafs,istod,istot,iscanf

!------------------------------------------------------
! check if file is open
!------------------------------------------------------

	ierr = 1
	if( iunit <= 0 ) return

!------------------------------------------------------
! read record
!------------------------------------------------------

	nvar0 = nvar

	!------------------------------------------------------
	! read keywords and skip empty lines
	!------------------------------------------------------

	line = ' '
	do while( line == ' ' )
	  read(iunit,'(a)',iostat=ierr) line
	  if( ierr /= 0 ) return

	  if( line == ' ' ) cycle
	  if( ts_has_keyword(line) ) then
	    call ts_parse_keyword(iunit,line,datetime,info)
	    line = ' '
	  end if
	end do

	!------------------------------------------------------
	! see if line is too long
	!------------------------------------------------------

	ierr = 2
	is = ichafs(line)
	iend = len_trim(line)
	if( iend > 2000 ) return

	!------------------------------------------------------
	! read time column
	!------------------------------------------------------

	ierr = 3
	ios = istod(line,dtime,is)		!read time column
	if( ios == -1 ) then			!time colum may be string
	  ios = istot(line,stime,is)		!read time column as string
	  if( ios /= 1 ) return
	  call string2datetime(stime,datetime,ierr)
	  if( ierr /= 0 ) return
	end if

	!------------------------------------------------------
	! read rest of variables
	!------------------------------------------------------

	ierr = 4
	nvar = iscanf(line(is:),f,nvar0)	!get values

	if( nvar < 0 ) nvar = -nvar-1		!read error in number -nvar
	if( nvar <= 0 ) return			!no data found
	if( nvar0 > 0 .and. nvar /= nvar0 ) return	!varying number of data

	!write(6,*) 'fffff: ',ierr,nvar0,nvar,datetime,iunit

c------------------------------------------------------
c set error code
c------------------------------------------------------

	ierr = 0

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
