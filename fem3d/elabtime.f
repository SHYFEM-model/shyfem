
	module timeutil

	implicit none

	private

	logical, save :: bdebug = .false.

	logical, save :: bdate
	integer, save :: date_elab
	integer, save :: time_elab
	integer, save :: datetime_elab(2)

        logical, save :: btmin
        logical, save :: btmax
        logical, save :: binclusive
        double precision, save :: atmin
        double precision, save :: atmax

        INTERFACE timeutil_check_time
        MODULE PROCEDURE timeutil_check_time_i,timeutil_check_time_d
        END INTERFACE

	contains

!************************************************************

	subroutine timeutil_date_and_time(date,time)

	integer date,time

        bdate = date .gt. 0
        if( bdate ) call dtsini(date,time)
	date_elab = date
	time_elab = time
        datetime_elab = (/date,time/)

	end subroutine timeutil_date_and_time

!************************************************************

	subroutine timeutil_minmax(stmin,stmax)

	character*(*) stmin,stmax

        atmin = 0.
        atmax = 0.
        btmin = stmin .ne. ' '
        btmax = stmax .ne. ' '
        if( btmin ) call dts_string2time(stmin,atmin)
        if( btmax ) call dts_string2time(stmax,atmax)

        if( bdebug ) then
          write(6,*) 'time limits: '
          write(6,*) stmin(1:len_trim(stmin)),btmin,atmin
          write(6,*) stmax(1:len_trim(stmax)),btmax,atmax
        end if

	end subroutine timeutil_minmax

!************************************************************

	function timeutil_check_time_i(it,itnew,itold)

! integer version

	logical timeutil_check_time_i
	integer it,itnew,itold

	double precision dtime,dtimenew,dtimeold

        dtime = it
	dtimenew = itnew
	dtimeold = itold

	timeutil_check_time_i = 
     +		timeutil_check_time_d(dtime,dtimenew,dtimeold)

	end function timeutil_check_time_i

!************************************************************

	function timeutil_check_time_d(dtime,dtimenew,dtimeold)

! double version (relativ)

	logical timeutil_check_time_d
	double precision dtime,dtimenew,dtimeold

	double precision atime,atimenew,atimeold

        call dts_convert_to_atime(datetime_elab,dtime,atime)
        call dts_convert_to_atime(datetime_elab,dtimenew,atimenew)
        call dts_convert_to_atime(datetime_elab,dtimeold,atimeold)

	timeutil_check_time_d =
     +		timeutil_check_time_a(atime,atimenew,atimeold)

	end function timeutil_check_time_d

!************************************************************

	function timeutil_check_time_a(atime,atimenew,atimeold)

! double version (absolute)

	logical timeutil_check_time_a
	double precision atime,atimenew,atimeold

	logical btimew

        btimew = .true.

        if( btmin ) btimew = btimew .and. atime >= atmin
        if( btmax ) btimew = btimew .and. atime <= atmax

	timeutil_check_time_a = btimew

	if( bdebug ) then
	  write(6,*) 'exclusive..........',btimew,binclusive
	  write(6,*) 'exclusive..........',atmin,atime,atmax
	end if

	if( .not. binclusive ) return

        if( btmin ) then
	  btimew = btimew .or. (atime < atmin .and. atmin < atimenew)
	end if
        if( btmax ) then
	  btimew = btimew .or. (atimeold < atmax .and. atmax < atime)
	end if

	timeutil_check_time_a = btimew

	end function timeutil_check_time_a

!************************************************************

	function timeutil_over_time_a(atime,atimenew,atimeold)

! double version (absolute)

	logical timeutil_over_time_a
	double precision atime,atimenew,atimeold

	logical btimew

        btimew = .true.

        if( btmax ) btimew = btimew .and. atime <= atmax

	timeutil_over_time_a = .not. btimew

	if( bdebug ) then
	  write(6,*) 'exclusive..........',btimew,binclusive
	  write(6,*) 'exclusive..........',atmin,atime,atmax
	end if

	if( .not. binclusive ) return

        if( btmax ) then
	  btimew = btimew .or. (atimeold < atmax .and. atmax < atime)
	end if

	timeutil_over_time_a = .not. btimew

	end function timeutil_over_time_a

!************************************************************

	end module timeutil
