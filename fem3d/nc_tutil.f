
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! convert nc files to fem files: time utilities
!
! contents :
!
!
! revision log :
!
! 16.05.2017	ggu	changed VERS_7_5_27
! 25.05.2017	ggu	changed VERS_7_5_28
! 14.11.2017	ggu	changed VERS_7_5_36
! 17.11.2017	ggu	changed VERS_7_5_37
! 05.12.2017	ggu	changed VERS_7_5_39
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.07.2018	ggu	revision control introduced
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
!
!*****************************************************************
!*****************************************************************
!*****************************************************************

!=================================================================
	module nc_time
!=================================================================

	implicit none

	integer, save :: time_type
	integer, save :: date0
	integer, save :: time0
	integer, save :: datetime0(2)
	double precision, save :: atime0
	double precision, save :: time_fact

!=================================================================
	end module nc_time
!=================================================================

c*****************************************************************

        subroutine setup_nc_time(ncid,bverb)

	use nc_time

	implicit none

	integer ncid
	logical bverb

	logical, save :: bdebug = .false.
	integer var_id
	integer ifact
	character*80 atext,tstring
	character*80 time_d,time_v

        call nc_get_time_name(time_d,time_v)
	call nc_get_var_id(ncid,time_v,var_id)
	call nc_get_var_attr(ncid,var_id,'units',atext)

	if( bdebug ) then
	  write(6,*) 'setup_nc_time: debug'
	  write(6,*) trim(time_d)
	  write(6,*) trim(time_v)
	  write(6,*) var_id
	  write(6,*) trim(atext)
	end if

	call parse_time_units(bverb,atext,time_type,datetime0,time_fact)

	date0 = datetime0(1)
	time0 = datetime0(2)
	call dtsini(date0,time0)
	call dts_to_abs_time(date0,time0,atime0)
	call dts_format_abs_time(atime0,tstring)

	ifact = time_fact
	if( bverb ) write(6,*) 'setup_nc_time: ',ifact
     +					,'  ',trim(atext)
     +					,'  ',trim(tstring)

	end

c*****************************************************************

        subroutine handle_nc_time(ncid,n,atime)

	use nc_time

	implicit none

	integer ncid
	integer n			!time record to be requested
        double precision atime		!parsed absolute time (return)

        double precision t

	call nc_get_time_rec(ncid,n,t)

	if( time_type .eq. 1 ) then
          call handle_general_time(t,atime)
	else if( time_type .eq. 2 ) then
          call handle_warf_time(t,atime)
	else
	  write(6,*) 'time_type: ',time_type
	  stop 'error stop handle_nc_time: cannot handle'
	end if

	end

c*****************************************************************

	subroutine parse_time_units(bverb,atext,itype,datetime0,fact)

	use iso8601

	implicit none

	logical bverb
	character*(*) atext
	integer itype			!type of time specification
	integer datetime0(2)		!reference time
	double precision fact		!factor to be used for time conversion

	integer ierr,off
	character*80 string

	itype = 1
	off = 1
	fact = 1.
	datetime0 = 0

	if( atext(1:10) .eq. 'days since' ) then
	  off = 12
	  fact = 86400.
	else if( atext(1:10) .eq. 'Days since' ) then
	  off = 12
	  fact = 86400.
	else if( atext(1:16) .eq. 'day as %Y%m%d.%f' ) then
	  itype = 2
	  fact = 86400.
	else if( atext(1:13) .eq. 'seconds since' ) then
	  off = 15
	  fact = 1.
	else if( atext(1:13) .eq. 'minutes since' ) then
	  off = 15
	  fact = 60.
	else if( atext(1:11) .eq. 'hours since' ) then
	  off = 13
	  fact = 3600.
	else if( atext(1:11) .eq. 'Hours since' ) then
	  off = 13
	  fact = 3600.
	else if( atext .eq. ' ' ) then	!no time coordinate
	  itype = 0
	  fact = 1.
	else
	  write(6,*) 'atext: ',trim(atext)
	  stop 'error stop parse_time_units: cannot parse'
	end if

	!write(6,*) itype,off,fact
	!write(6,*) trim(atext)

	if( itype == 1 ) then
	  string = atext(off:)
	  call clean_time(string)
	  call string2date(string,datetime0,ierr)
	  if( ierr /= 0 ) then
	    write(6,*) 'error parsing time reference: ',trim(atext)
	    stop 'error stop parse_time_units: parsing time'
	  end if
	end if

	if( bverb ) then
	  call date2string(datetime0,string)
	  write(6,*) 'parsing date0: ',trim(atext)
	  write(6,*) 'parsed date:   ',trim(string)
	end if

	end

c*****************************************************************

	subroutine clean_time(string)

	implicit none

	character*(*) string

	integer len

	len = len_trim(string)

	if( string(len-1:len) == '.0' ) then
	  string(len-1:) = ' ' 
	  len = len - 2
	end if
	if( string(len-1:len) == ':0' ) then
	  string(len-1:) = ':00' 
	  len = len + 1
	end if

	end

c*****************************************************************

	subroutine parse_date_time(atext,year,month,day)

	implicit none

	character*(*) atext
	integer year,month,day

	read(atext(1:4),'(i4)') year
	read(atext(6:7),'(i2)') month
	read(atext(9:10),'(i2)') day

	end

c*****************************************************************

        subroutine handle_warf_time(time,atime)

	use nc_time

        implicit none

        double precision time
        double precision atime

	double precision secs
        integer date

        date = time				!converts to full days
        secs = (time-date) * time_fact

	call dts_to_abs_time(date,0,atime)
	atime = atime + secs

	end

c*****************************************************************

        subroutine handle_general_time(time,atime)

	use nc_time

        implicit none

        double precision time
        double precision atime

	double precision secs

	secs = time * time_fact
	atime = atime0 + secs

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine print_time_records(ncid)

c print time and date

        implicit none

        integer ncid

        integer nit,n
        double precision atime
        character*20 line

        call nc_get_time_recs(ncid,nit)
        write(6,*) 'time records found: ',nit

        do n=1,nit
          call handle_nc_time(ncid,n,atime)
          call dts_format_abs_time(atime,line)
          write(6,*) n,atime,line
        end do

        end

c*****************************************************************

        subroutine print_minmax_time_records(ncid)

c print time and date

        implicit none

        integer ncid

        integer nit
        double precision atime
        character*20 line

        call nc_get_time_recs(ncid,nit)

        if( nit == 0 ) then
          write(6,*) 'no time record found'
        else if( nit == 1 ) then
          call handle_nc_time(ncid,1,atime)
          call dts_format_abs_time(atime,line)
          write(6,*) 'one time record found: ',atime,line
        else
          write(6,*) 'time records found: ',nit
          call handle_nc_time(ncid,1,atime)
          call dts_format_abs_time(atime,line)
          write(6,*) 'first time record:  ',atime,line
          call handle_nc_time(ncid,nit,atime)
          call dts_format_abs_time(atime,line)
          write(6,*) 'last time record:   ',atime,line
        end if

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine create_date_string(ncid,it,datetime)

	use iso8601

        implicit none

        integer ncid,it
        integer datetime(2)

        integer ierr
        double precision atime
        character*80 line

        datetime = 0
        if( it == 0 ) return

        call handle_nc_time(ncid,it,atime)
        call dts_format_abs_time(atime,line)
        call string2date(line,datetime,ierr)

        if( ierr /= 0 ) then
          write(6,*) 'error converting date string: ',trim(line)
          stop 'error stop write_variables: date string'
        end if

        end

c*****************************************************************

