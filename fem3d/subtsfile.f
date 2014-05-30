
c routines for reading time series

c*************************************************************

	subroutine ts_get_file_info(file,nvar)

	implicit none

	character*(*) file	!file name
	integer nvar		!variables (columns) in file (except time)

	character*132 line
	integer iunit,i,nrec
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

        open(iunit,file=file,form='formatted',status='old',err=2)

c------------------------------------------------------
c loop until first valid record
c------------------------------------------------------

	do while(.true.)
	  nrec = nrec + 1
	  read(iunit,'(a)',err=3,end=3) line
	  i = ichafs(line)
	  if( i > 0 .and. line(i:i) /= '#' ) exit
	end do

c------------------------------------------------------
c set nvar and close file
c------------------------------------------------------

	nvar = iscand(line,d,0)		!count values on line
	if( nvar < 0 ) nvar = -nvar-1	!read error in number -nvar
	nvar = nvar - 1			!do not count time column

	close(iunit)

c------------------------------------------------------
c end of routine
c------------------------------------------------------

    2	continue
	return
    3	continue
	write(6,*) 'read error in line ',nrec,' of file ',file
	stop 'error stop ts_get_file_info: read error'
	end

c*************************************************************

	subroutine ts_open_file(file,nvar,date,time,iunit)

	implicit none

	character*(*) file	!file name
	integer nvar		!variables (columns) in file (except time)
	integer date,time
	integer iunit

	character*132 line
	character*10 key
	integer i,j,nrec,ioff
	double precision d(3)

	integer iscand,ichafs
	integer ifileo

c------------------------------------------------------
c open file
c------------------------------------------------------

	nrec = 0
	nvar = 0
	date = 0
	time = 0

	iunit = ifileo(0,file,'formatted','old')

	if( iunit <= 0 ) return

c------------------------------------------------------
c check if header exists and get nvar
c------------------------------------------------------

	do while(.true.)
	  nrec = nrec + 1
	  read(iunit,'(a)',err=3,end=3) line
	  i = ichafs(line)
	  if( i <= 0 ) cycle
	  if( line(i:i) == '#' ) then
	    call ts_get_keyword(line,key,ioff)
	    if( key == 'date' ) then		!date
	      j = iscand(line(ioff:),d,2)
	      if( j >= 1 ) date = d(1)
	      if( j >= 2 ) time = d(2)
	    else if( key == ' ' ) then		!nothing
	    else
	      write(6,*) 'not recognized keyword: ',key,' in file ',file
	      stop 'error stop ts_open_file: not recognized keyword'
	    end if
	  else
	    exit
	  end if
	end do

c------------------------------------------------------
c set nvar and backspace file
c------------------------------------------------------

	nvar = iscand(line,d,0)		!count values on line
	if( nvar < 0 ) nvar = -nvar-1	!read error in number -nvar

	backspace(iunit)

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	return
    3	continue
	write(6,*) 'read error in line ',nrec,' of file ',file
	stop 'error stop ts_get_file_info: read error'
	end

c*************************************************************

	subroutine ts_peek_next_record(iunit,nvar,it,f,ierr)

c peeks into one record of time series file

	implicit none

	integer iunit
	integer nvar		!variables (columns) in file (except time)
	integer it
	real f(nvar)
	integer ierr

	call ts_read_next_record(iunit,nvar,it,f,ierr)

	backspace(iunit)

	end

c*************************************************************

	subroutine ts_read_next_record(iunit,nvar,it,f,ierr)

c reads one record of time series file

	implicit none

	integer iunit
	integer nvar		!variables (columns) in file (except time)
	integer it
	real f(nvar)
	integer ierr

	integer i
	double precision t

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

	read(iunit,*,end=2,err=3) t,(f(i),i=1,nvar)
	it = nint(t)
	ierr = 0
	return

c------------------------------------------------------
c error handling
c------------------------------------------------------

    2	continue
	ierr = -1
	return

    3	continue
	ierr = 3
	return

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	end

c*************************************************************

	subroutine ts_get_keyword(line,key,ioff)

c gets keyword from line - rest of line after ioff
c
c keyword looks like:   "#key:"

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

