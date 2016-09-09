
!==============================================================
	module custom_dates
!==============================================================

	implicit none

	integer, save :: idate = 0
	integer, save :: ndate = 0
	integer, save, allocatable :: restime(:)

!==============================================================
	contains
!==============================================================

	subroutine custom_dates_init(it,file)

	integer it
	character*(*) file

	integer ndim

	!write(6,*) 'custom: ',it,file

	if( idate /= 0 ) return		!called more than once

	idate = -1
	if( file == ' ' ) return			!no file given

	call get_custom_dates(file,-1,ndim,restime)
	allocate(restime(ndim))
	call get_custom_dates(file,ndim,ndate,restime)

	write(6,*) 'custom dates used: ',ndate,'  ',trim(file)

	idate = 0

	do
	  idate = idate + 1
	  if( idate > ndate ) exit
	  if( it < restime(idate) ) exit
	end do

	end subroutine custom_dates_init

c**************************************************************

	subroutine custom_dates_over(it,bover)

	integer it
	logical bover

	character*80 file

c---------------------------------------------------------------
c initialize - convert date to relative time
c---------------------------------------------------------------

	bover = .false.

	if( idate == -1 ) return

c---------------------------------------------------------------
c see if we have to reset
c---------------------------------------------------------------

	if( idate > ndate ) return
	if( it < restime(idate) ) return

c---------------------------------------------------------------
c ok, reset needed - advance to next reset time
c---------------------------------------------------------------

	do
	  idate = idate + 1
	  if( idate > ndate ) exit
	  if( it < restime(idate) ) exit
	end do

	bover = .true.

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end subroutine custom_dates_over

c**************************************************************

	subroutine get_custom_dates(file,ndim,n,restime)

c gets custom reset time from file

	character*(*) file
	integer ndim		!ndim==0 => check how many dates are given
	integer n
	integer restime(n)

	integer ianz,ios,nline,i
	integer date,time
	integer year,month,day,hour,min,sec
	integer itres,itold
	double precision d(2)
	character*80 line
	logical bdebug

	integer iscand

	n = 0
	nline = 0
	bdebug = .true.
	if( ndim == -1 ) bdebug = .false.

	open(1,file=file,status='old',form='formatted',iostat=ios)

	if( bdebug ) then
	  if( ios /= 0 ) then
	    write(6,*) 'cannot open custom reset file: ',trim(file)
	    stop 'error stop get_custom_dates: opening file'
	  else
	    write(6,*) 'reading custom reset file: ',trim(file)
	  end if
	end if
	if( ios /= 0 ) return

	do
	  read(1,'(a)',iostat=ios) line
	  nline = nline + 1
	  if( ios /= 0 ) exit
	  ianz = iscand(line,d,2)
	  if( ianz == 0 ) then
	    cycle
	  else if( ianz == 1 ) then
	    date = nint(d(1))
	    time = 0
	  else if( ianz == 2 ) then
	    date = nint(d(1))
	    time = nint(d(2))
	  else
	    write(6,*) 'parse error: ',ianz
	    write(6,*) 'line: ',trim(line)
	    write(6,*) 'file: ',trim(file)
	    stop 'error stop get_custom_dates: parse error'
	  end if

	  n = n + 1
	  if( ndim == -1 ) cycle
	  if( n > ndim ) then
	    write(6,*) 'n,ndim: ',n,ndim
	    stop 'error stop get_custom_dates: dimension error ndim'
	  end if

	  call unpacktime(time,hour,min,sec)
	  call unpackdate(date,year,month,day)
	  call dts2it(itres,year,month,day,hour,min,sec)

	  restime(n) = itres		!insert relative time
	end do

	if( ios > 0 ) then
	  write(6,*) 'read error...'
	  write(6,*) 'file: ',trim(file)
	  write(6,*) 'line number: ',nline
	  stop 'error stop get_custom_dates: read error'
	end if

	if( bdebug ) then
	  write(6,*) 'custom reset times: ',n
	  itold = restime(1) - 1
	  do i=1,n
	    itres = restime(i)
	    call dts2dt(itres,year,month,day,hour,min,sec)
	    write(6,1000) i,itres,year,month,day,hour,min,sec
 1000	    format(i5,i12,6i5)
	    if( itres <= itold ) then
	      write(6,*) 'times in custom reset must be ascending...'
	      stop 'error stop get_custom_dates: wrong order'
	    end if
	  end do
	end if

	close(1)

	end subroutine get_custom_dates

!==============================================================
	end module custom_dates
!==============================================================

