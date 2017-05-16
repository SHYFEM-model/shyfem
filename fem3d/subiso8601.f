!
! ISO 8601 for date and time specification
!
! revision log :
!
! 10.05.2017    ggu     started
! 15.05.2017    ggu     finished
!
! notes :
!
! implements partly ISO 8601 for date and time specification
!
! implements both extended form and basic form
!
! extended: YYYY-MM-DDThh:mm:ss or YYYY-MM-DD::hh:mm:ss or YYYY-MM-DD hh:mm:ss
! basic:    YYYYMMDDThhmmss     or YYYYMMDD::hhmmss     or YYYYMMDD hhmmss
!
! separator can be T, ::, or blanks
!
! mixed representation (date in extend and time in basic etc..) is not allowed
! date must always given fully (until day)
! time can be abbreviated (hh, hh:mm, hhmm)
!
!*********************************************************************

	subroutine string2dt(string,dt,ierr)

! converts date string to integer representation

	implicit none

	character*(*) string	!date string
	integer dt(8)		!year,month,day,hour,min,sec,msec,tz (return)
	integer ierr		!error if /= 0 (return)

	logical, parameter :: bdebug = .false.
	logical bextend
	integer n,nl
	character(len=max(20,len(string))) ll,time

	dt = 0

	ll = adjustl(string)
	n = len_trim(ll)

!	-------------------------------------------------------
!	parse date
!	-------------------------------------------------------

        ierr = 1
        if( n .lt. 8 ) goto 9  !we insist having at least the full date

        bextend = ( ll(5:5) == '-' .and. ll(8:8) == '-' )

        ierr = 2
        if( bextend ) then
          read(ll(1:10) ,'(i4,1x,i2,1x,i2)',err=9) dt(1:3)
	  nl = 10
	else						!try basic
          read(ll(1:8) ,'(i4,i2,i2)',err=9) dt(1:3)
	  nl = 8
	end if

        if( n .le. nl ) goto 1

!	-------------------------------------------------------
!	parse separator
!	-------------------------------------------------------

	ierr = 4
	nl = nl + 1
	if( ll(nl:nl+1) == '::' ) then
	  time = ll(nl+2:)
	else if( ll(nl:nl) == 'T' ) then
	  time = ll(nl+1:)
	else if( ll(nl:nl) /= ' ' ) then
	  goto 9
	else
	  time = adjustl(ll(nl:))
	end if

!	-------------------------------------------------------
!	parse time
!	-------------------------------------------------------

	n = len_trim(time)
	if( n == 0 ) goto 1

	ierr = 5
	if( bextend ) then
	  if( n >=6 .and. time(6:6) /= ':' ) goto 9
	  if( n >=3 .and. time(3:3) /= ':' ) goto 9
	  if( n > 6 ) then
            read(time(1:n) ,'(i2,1x,i2,1x,i2)',err=9) dt(4:6)
	    nl = 8
	  else if( n > 3 ) then
            read(time(1:n) ,'(i2,1x,i2)',err=9) dt(4:5)
	    nl = 5
	    if( time(nl+1:nl+1) == ':' ) nl = nl + 1
	  else
            read(time(1:n) ,'(i2)',err=9) dt(4)
	    nl = 2
	    if( time(nl+1:nl+1) == ':' ) nl = nl + 1
	  end if
	else
	  if( n > 4 ) then
            read(time(1:n) ,'(i2,i2,i2)',err=9) dt(4:6)
	    nl = 6
	  else if( n > 2 ) then
            read(time(1:n) ,'(i2,i2)',err=9) dt(4:5)
	    nl = 4
	  else
            read(time(1:n) ,'(i2)',err=9) dt(4)
	    nl = 2
	  end if
	end if

!	-------------------------------------------------------
!	parse rest
!	-------------------------------------------------------

	if( nl == n ) goto 1

    2   continue
 
	goto 9
!	not yet ready for milliseconds and time zone

!	-------------------------------------------------------
!	end of routine
!	-------------------------------------------------------

	return
    1   continue
	ierr = 0
	return
    9   continue
        if( bdebug ) then
          write(6,*) '*** cannot parse date: ',ierr,trim(string)
          write(6,*) '    format should be YYYY-MM-DD::hh:mm:ss'
          write(6,*) '    or iso8601 format YYYY-MM-DDThh:mm:ss'
          write(6,*) '    possible also YYYY-MM-DD[::[hh[:mm[:ss]]]]'
          write(6,*) '    or YYYY-MM-DD[T[hh[:mm[:ss]]]]'
        end if
        return
	end

!*********************************************************************

	subroutine dt2string(dt,string)

	integer dt(8)		!year,month,day,hour,min,sec,msec,tz (return)
	character*(*) string	!date string

	integer ilen,i
	character*1, save :: sep = ':'
	character*20 lineaux

	lineaux = ' '

	if( sep == 'T' ) then
	  write(lineaux,1100) dt(1:6)
	  ilen = 19
	else
	  write(lineaux,1000) dt(1:6)
	  ilen = 20
	end if

        do i=1,ilen
          if( lineaux(i:i) .eq. ' ' ) lineaux(i:i) = '0'
        end do

	string = lineaux

	return
 1000   format(i4,1h-,i2,1h-,i2,2h::,i2,1h:,i2,1h:,i2)
 1100   format(i4,1h-,i2,1h-,i2,1hT,i2,1h:,i2,1h:,i2)
	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine dt2datetime(dt,datetime)

	implicit none

	integer dt(8)		!year,month,day,hour,min,sec,msec,tz (return)
	integer datetime(2)	!date,time

	datetime(1) = 10000*dt(1) + 100*dt(2) + dt(3)
	datetime(2) = 10000*dt(4) + 100*dt(5) + dt(6)

	end
	
!*********************************************************************

	subroutine datetime2dt(datetime,dt)

	implicit none

	integer datetime(2)	!date,time
	integer dt(8)		!year,month,day,hour,min,sec,msec,tz (return)

	integer date,time,iaux

	date = datetime(1)
        iaux = date / 100
        dt(3) = date - iaux * 100
        dt(1) = iaux / 100
        dt(2) = iaux - dt(1) * 100

	time = datetime(2)
        iaux = time / 100
        dt(6) = time - iaux * 100
        dt(4) = iaux / 100
        dt(5) = iaux - dt(4) * 100

	if( date /= 10000*dt(1) + 100*dt(2) + dt(3) ) then
	  write(6,*) date,dt(1:3)
	  stop 'error stop datetime2dt: internal error (1)'
	end if
	if( time /= 10000*dt(4) + 100*dt(5) + dt(6) ) then
	  write(6,*) time,dt(4:6)
	  stop 'error stop datetime2dt: internal error (2)'
	end if

	end
	
!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine string2datetime(string,datetime,ierr)

	implicit none

	character*(*) string	!date string
	integer datetime(2)
	integer ierr

	integer dt(8)		!year,month,day,hour,min,sec,msec,tz (return)

	call string2dt(string,dt,ierr)
	if( ierr /= 0 ) return
	call dt2datetime(dt,datetime)

	end

!*********************************************************************

	subroutine datetime2string(datetime,string)

	implicit none

	integer datetime(2)
	character*(*) string	!date string

	integer dt(8)		!year,month,day,hour,min,sec,msec,tz (return)

	call datetime2dt(datetime,dt)
	call dt2string(dt,string)

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine test_iso8601

	implicit none

	integer ierr,dt(8)

	call test_iso8601_check('2017-04-23')
	call test_iso8601_check('2017-04-23::')
	call test_iso8601_check('2017-04-23::12')
	call test_iso8601_check('2017-04-23::12:')
	call test_iso8601_check('2017-04-23T12:30')
	call test_iso8601_check('2017-04-23T12:30:')
	call test_iso8601_check('2017-04-23T12:30:45')
	call test_iso8601_check('2017-04-23T')
	call test_iso8601_check('20170423T')
	call test_iso8601_check('20170423')
	call test_iso8601_check('20170423T12')
	call test_iso8601_check('20170423T1230')
	call test_iso8601_check('20170423T123045')

	call test_iso8601_check('2017-04-23T12:30:4')
	call test_iso8601_check('2017-04-23T12:3')
	call test_iso8601_check('2017-04-23T1230')
	call test_iso8601_check('2017-04-23T12:30:45.3')
	call test_iso8601_check('2017-04-23T12:30:45Z')
	call test_iso8601_check('2017-04-23T12:30:45+01')

	end

!*********************************************************************

	subroutine test_iso8601_check(string)

	implicit none

	character*(*) string
	integer ierr,dt(8)

	call string2dt(string,dt,ierr)
	write(6,'(9i5,2a)') ierr,dt,'  ',trim(string)

	end

!*********************************************************************
!	program main_iso8601
!	call test_iso8601
!	end
!*********************************************************************

