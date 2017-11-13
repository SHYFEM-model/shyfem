
!================================================================
	module elabtime
!================================================================

	implicit none

	!private

	logical, save :: bdebugtime = .false.

	logical, save :: bdate_elab
	integer, save :: date_elab
	integer, save :: time_elab
	integer, save :: datetime_elab(2)
	double precision, save :: atime00

        logical, save :: btmin = .false.
        logical, save :: btmax = .false.
        logical, save :: binclusive_elab = .false.

        double precision, save :: atmin = 0.
        double precision, save :: atmax = 0.

        INTERFACE elabtime_check_time
        MODULE PROCEDURE elabtime_check_time_i,elabtime_check_time_d
        END INTERFACE

        INTERFACE elabtime_in_time
	MODULE PROCEDURE elabtime_check_time_in,elabtime_check_time_ni
        END INTERFACE

        INTERFACE elabtime_over_time
	MODULE PROCEDURE elabtime_over_time_in,elabtime_over_time_ni
        END INTERFACE

!================================================================
	contains
!================================================================

!************************************************************

	subroutine elabtime_set_inclusive(binclusive)

! sets inclusive flag

	logical binclusive

	binclusive_elab = binclusive

	end subroutine elabtime_set_inclusive

!************************************************************

	subroutine elabtime_date_and_time(date,time)

! initializes elabtime module

	integer date,time

        bdate_elab = date .gt. 0
        if( bdate_elab ) call dtsini(date,time)

	date_elab = date
	time_elab = time
        datetime_elab = (/date,time/)
	call dts_to_abs_time(date,time,atime00)

	end subroutine elabtime_date_and_time

!************************************************************
!************************************************************
!************************************************************

	subroutine elabtime_set_minmax(stmin,stmax)

! converts stmin/stmax to absolute time

	character*(*) stmin,stmax

        btmin = stmin .ne. ' '
        btmax = stmax .ne. ' '
        if( btmin ) call dts_string2time(stmin,atmin)
        if( btmax ) call dts_string2time(stmax,atmax)

        if( bdebugtime ) then
          write(6,*) 'time limits: '
          write(6,*) stmin(1:len_trim(stmin)),btmin,atmin
          write(6,*) stmax(1:len_trim(stmax)),btmax,atmax
        end if

	end subroutine elabtime_set_minmax

!************************************************************
!************************************************************
!************************************************************

	function elabtime_check_time_i(it,itnew,itold)

! check if it is inside time window - integer version (relative)

	logical elabtime_check_time_i
	integer it,itnew,itold

	double precision dtime,dtimenew,dtimeold

        dtime = it
	dtimenew = itnew
	dtimeold = itold

	elabtime_check_time_i = 
     +		elabtime_check_time_d(dtime,dtimenew,dtimeold)

	end function elabtime_check_time_i

!************************************************************

	function elabtime_check_time_d(dtime,dtimenew,dtimeold)

! check if dtime is inside time window - double version (relative)

	logical elabtime_check_time_d
	double precision dtime,dtimenew,dtimeold

	double precision atime,atimenew,atimeold

        call dts_convert_to_atime(datetime_elab,dtime,atime)
        call dts_convert_to_atime(datetime_elab,dtimenew,atimenew)
        call dts_convert_to_atime(datetime_elab,dtimeold,atimeold)

	elabtime_check_time_d =
     +		elabtime_check_time_in(atime,atimenew,atimeold)

	end function elabtime_check_time_d

!************************************************************
!************************************************************
!************************************************************

	function elabtime_check_time_in(atime,atimenew,atimeold)

! check if atime is inside time window - double version (absolute)

	logical elabtime_check_time_in
	double precision atime,atimenew,atimeold

	logical btimew

        btimew = .true.

        if( btmin ) btimew = btimew .and. atime >= atmin
        if( btmax ) btimew = btimew .and. atime <= atmax

	elabtime_check_time_in = btimew

	if( bdebugtime ) then
	  write(6,*) 'exclusive..........',btimew,binclusive_elab
	  write(6,*) 'exclusive..........',atmin,atime,atmax
	end if

	if( btimew ) return			!already true
	if( .not. binclusive_elab ) return

        if( btmin ) then
	  btimew = btimew .or. (atime < atmin .and. atmin < atimenew)
	end if
        if( btmax ) then
	  btimew = btimew .or. (atimeold < atmax .and. atmax < atime)
	end if

	elabtime_check_time_in = btimew

	end function elabtime_check_time_in

!************************************************************

	function elabtime_check_time_ni(atime)

! check if atime is inside time window - double version (absolute)

	logical elabtime_check_time_ni
	double precision atime

	elabtime_check_time_ni 
     +			= elabtime_check_time_in(atime,atime,atime)

	end function elabtime_check_time_ni

!************************************************************

	function elabtime_over_time_in(atime,atimenew,atimeold)

! check if atime is beyond time window - double version (absolute)

	logical elabtime_over_time_in
	double precision atime,atimenew,atimeold

	logical btimew

        btimew = .true.

        if( btmax ) btimew = btimew .and. atime <= atmax

	elabtime_over_time_in = .not. btimew

	if( bdebugtime ) then
	  write(6,*) 'exclusive..........',btimew,binclusive_elab
	  write(6,*) 'exclusive..........',atmin,atime,atmax
	end if

	if( .not. binclusive_elab ) return

        if( btmax ) then
	  btimew = btimew .or. (atimeold < atmax .and. atmax < atime)
	end if

	elabtime_over_time_in = .not. btimew

	end function elabtime_over_time_in

!************************************************************

	function elabtime_over_time_ni(atime)

! check if atime is beyond time window - double version (absolute)

	logical elabtime_over_time_ni
	double precision atime

	elabtime_over_time_ni 
     +			= elabtime_over_time_in(atime,atime,atime)

	end function elabtime_over_time_ni

!************************************************************

!================================================================
	end module elabtime
!================================================================

