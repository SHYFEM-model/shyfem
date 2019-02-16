
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c time management routines
c
c contents :
c
c
c revision log :
c
c 23.09.1997	ggu	boundn deleted -> no access to data structure
c 20.03.1998	ggu	minor changes to priout
c 29.04.1998	ggu	new module for semi-implicit time-step
c 07.05.1998	ggu	check for error on return of nrdvecr
c 19.06.1998	ggu	version number is character
c 22.01.1999	ggu	oxygen section added
c 26.01.1999	ggu	new comp3d added
c 11.08.1999	ggu	new compatibility array hlhv initialized
c 19.11.1999	ggu	new routines for section vol
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 04.02.2000    ggu     no priout, dobefor/after, pritime, endtime
c 15.05.2000    ggu     hm3v substituted
c 26.05.2000    ggu     copright statement adjourned
c 21.11.2001    ggu     routines to handle advective index (aix)
c 27.11.2001    ggu     routine to handle info file (getinfo)
c 11.10.2002    ggu     aix routines deleted
c 07.02.2003    ggu     routine added: changeimp, getaz; deleted getaza
c 10.08.2003    ggu     call adjust_chezy instead sp135r
c 14.08.2003    ggu     femver transfered to subver, not called in nlsh2d
c 20.08.2003    ggu     tsmed substituted by ts_shell
c 01.09.2003    ggu     call wrousa
c 03.09.2004    ggu     call admrst, comp3d renamed to init_3d (not used)
c 03.09.2004    ggu     nlv, hlv initialized in nlsh2d (FIXME)
c 28.09.2004    ggu     read lagrangian section
c 01.12.2004    ggu     new routine set_timestep for variable time step
c 17.01.2005    ggu     get_stab_index to newcon.f, error stop in set_timestep
c 14.03.2005    ggu     syncronize idt with end of simulation (set_timestep)
c 07.11.2005    ggu     handle new section sedtr for sediments
c 23.03.2006    ggu     changed time step to real
c 23.05.2007    ggu     recall variable time step pars at every time step
c 02.10.2007    ggu     bug fix in set_timestep for very small rindex
c 10.04.2008    ccf     output in netcdf format
c 28.04.2008    ggu     in set_timestep new call to advect_stability()
c 03.09.2008    ggu     in nlsh2d different error message
c 20.11.2008    ggu     init_3d deleted, nlv initialized to 0
c 18.11.2009    ggu     new format in pritime (write also time step)
c 22.02.2010    ggu     new call to hydro_stability to compute time step
c 22.02.2010    ccf     new routine for tidal pot. (tideforc), locaus deleted
c 26.02.2010    ggu     in set_timestep compute and write ri with old dt
c 22.03.2010    ggu     some comments for better readability
c 29.04.2010    ggu     new routine set_output_frequency() ... not finished
c 04.05.2010    ggu     shell to compute energy
c 22.02.2011    ggu     in pritime() new write to terminal
c 20.05.2011    ggu     changes in set_timestep(), element removal, idtmin
c 31.05.2011    ggu     changes for BFM
c 01.06.2011    ggu     idtmin introduced
c 12.07.2011    ggu     new routine next_output(), revised set_output_frequency
c 14.07.2011    ggu     new routines for original time step
c 13.09.2011    ggu     better error check, rdtitl() more robust
c 23.01.2012    ggu     new section "proj"
c 24.01.2012    ggu     new routine setup_parallel()
c 10.02.2012    ggu     new routines to initialize and access time common block
c 05.03.2014    ggu     code prepared to repeat time step (irepeat) - not ready
c 05.03.2014    ggu     new routines get_last/first_time()
c 10.04.2014    ccf     new section "wrt" for water renewal time
c 29.10.2014    ggu     do_() routines transfered from newpri.f
c 10.11.2014    ggu     time management routines transfered to this file
c 23.09.2015    ggu     new routine convert_time_d() for double
c 24.09.2015    ggu     routines re-written for double precision
c 20.10.2015    ggu     new routines to set/get id
c 04.11.2015    ggu     allow for initial output in adjust_itmidt()
c 04.11.2017    ggu     new routine init_output_i()
c 03.10.2018    ggu     some instances of itanf and itend eliminated
c
c info :
c
c       da_out(1) = idtout      ! time step of output
c       da_out(2) = itmout      ! first output
c       da_out(3) = itout       ! next output
c       da_out(4) = 0           ! unit or shyfem file id (optional)
c
c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine adjust_itmidt_d(itmout,idtout,itout)

c sets-up output frequency and first output

	implicit none

	include 'femtime.h'

	double precision itmout		!minimum time for output
	double precision idtout		!time step for output
	double precision itout		!first output

        logical binit

        binit = ( itmout /= -1. )                !output initial time?

	if( itmout .eq. -1. ) itmout = dtanf
	if( itmout .lt. dtanf ) itmout = dtanf

	itout = itmout
	if( itmout .eq. dtanf .and. .not. binit ) itout = itout + idtout
	if( itout .gt. dtend .and. idtout .gt. 0 ) idtout = 0

	end

c********************************************************************

	subroutine set_output_frequency_d(itmout,idtout,da_out)

c sets-up array for output frequency

	implicit none

	double precision itmout		!minimum time for output
	double precision idtout		!time step for output
	double precision da_out(4)	!array where info is stored

	double precision itout

	call adjust_itmidt_d(itmout,idtout,itout)

	da_out(1) = idtout	! time step of output
	da_out(2) = itmout	! first output
	da_out(3) = itout	! next output
	da_out(4) = 0		! unit (optional)

	end

c********************************************************************

	subroutine assure_initial_output_d(da_out)

c makes sure that output will be done also for it == itanf

	implicit none

	double precision da_out(4)

	double precision itmout,itout

	itmout = da_out(2)
	itout  = da_out(3)

	itout = itmout

	da_out(3) = itout

	end

c********************************************************************

	subroutine increase_output_d(da_out)

c makes sure that itout > itmout

	implicit none

	double precision da_out(4)

	include 'femtime.h'

	double precision idtout,itmout,itout

	idtout = da_out(1)
	itmout = da_out(2)
	itout  = da_out(3)

	if( itout > itmout ) return

	itout = itmout + idtout
	if( itout .gt. itend ) idtout = 0

	da_out(1) = idtout	! time step of output
	da_out(3) = itout	! next output

	end

c********************************************************************

	function is_over_output_d(da_out)

c checks if output phase has started (it > itmout)

	implicit none

	logical is_over_output_d
	double precision da_out(4)

	include 'femtime.h'

	is_over_output_d = t_act > da_out(2)

	end

c********************************************************************

	function is_in_output_d(da_out)

c checks if we arrived at output phase (it >= itmout)

	implicit none

	logical is_in_output_d
	double precision da_out(4)

	include 'femtime.h'

	is_in_output_d = t_act >= da_out(2)

	end

c********************************************************************

	function has_output_d(da_out)

c checks if variable has any output at all

	implicit none

	logical has_output_d
	double precision da_out(4)

	has_output_d = da_out(1) > 0	!idtout > 0

	end

c********************************************************************

	function next_output_d(da_out)

c checks if time has come for output

	implicit none

	logical next_output_d
	double precision da_out(4)

	include 'femtime.h'

	double precision idtout,itout

	next_output_d = .false.
	idtout = da_out(1)
	itout  = da_out(3)

	if( idtout .le. 0 ) return
	if( itout .gt. t_act ) return

	do while( itout .le. t_act )
	  itout = itout + idtout
	end do

	da_out(3) = itout
	next_output_d = .true.

	end

c********************************************************************

	subroutine info_output_d(da_out)

c writes info on da_output

	implicit none

	double precision da_out(4)

	include 'femtime.h'

	logical has_output_d,next_output_d
	double precision dtime,atime
	character*20 aline

	write(6,*) '------ info_output start------'
	write(6,*) da_out
	dtime = da_out(1)
	write(6,*) 'idtout = ',dtime
	dtime = da_out(2)
	call dts_format_abs_time(atime0+dtime,aline)
	write(6,*) 'itmout = ',dtime,aline
	dtime = da_out(3)
	call dts_format_abs_time(atime0+dtime,aline)
	write(6,*) 'itnext = ',dtime,aline
	dtime = t_act
	call dts_format_abs_time(atime0+dtime,aline)
	write(6,*) 't_act  = ',dtime,aline
	write(6,*) 'has output =  ',has_output_d(da_out)
	write(6,*) 'next output = ',next_output_d(da_out)
	write(6,*) '------ info_output end ------'

	end

c********************************************************************

        subroutine set_id_output_d(da_out,id)

        implicit none

        integer da_out(4)
        integer id

        da_out(4) = id

        end

c********************************************************************

        subroutine get_id_output_d(da_out,id)

        implicit none

        integer da_out(4)
        integer id

        id = da_out(4)

        end

c********************************************************************

	subroutine init_output_d(itmname,idtname,da_out)

c gets time values and transforms them

	implicit none

	character*(*) itmname,idtname	!names to parse
	double precision da_out(4)	!array with time information

	double precision itmout,idtout

	call convert_date_d(itmname,itmout)
	call convert_time_d(idtname,idtout)

	call set_output_frequency_d(itmout,idtout,da_out)

	end

c********************************************************************

	subroutine convert_date_d(name,dit)

c converts date to relative time

	implicit none

	character*(*) name
	double precision dit

	integer ierr
	integer date,time
	double precision atime,atime0,dtime
	character*30 text
	logical bdebug

	double precision dgetpar

	bdebug = .true.
	bdebug = .false.

	call getfnm(name,text)

	if( text .ne. ' ' ) then
	  !call dtsgunf(it,text,ierr)
	  call dts_get_date(date,time)
	  call dts_to_abs_time(date,time,atime0)
	  call dts_string2time(text,atime,ierr)
	  if( ierr .ne. 0 ) goto 99
	  dtime = atime - atime0
	  dit = dtime
	  call dputpar(name,dit)
	else
	  dit = dgetpar(name)
	end if

	if( .not. bdebug ) return

	write(6,*) 'convert_date_d: '
	write(6,*) trim(name)
	write(6,*) text
	write(6,*) dit

	return
   99	continue
	write(6,*) 'name: ',trim(name)
	write(6,*) 'text: ',text
        write(6,*) '*** cannot parse date: ',ierr,text
        write(6,*) '    format should be YYYY-MM-DD::hh:mm:ss'
        write(6,*) '    possible also YYYY-MM-DD[::hh[:mm[:ss]]]'
        write(6,*) '    or it should be an integer (relative time)'
	stop 'error stop convert_date: cannot parse'
	end

c********************************************************************

	subroutine convert_time_d(name,didt)

c converts time period to relative time difference

	implicit none

	character*(*) name
	double precision didt

	integer ierr
	character*40 text
	logical bdebug

	double precision dgetpar

	bdebug = .true.
	bdebug = .false.

	call getfnm(name,text)

	if( bdebug ) write(6,*) 'converting time for ',name

	if( text .ne. ' ' ) then
	  call dtstimespand(didt,text,ierr)	!still in integer
	  if( ierr .ne. 0 ) goto 99
	  if( bdebug ) then
	    write(6,*) 'time span as string found'
	    write(6,*) name
	    write(6,*) text
	    write(6,*) didt
	  end if
	  call dputpar(name,didt)
	else
	  didt = dgetpar(name)
	end if

	if( bdebug ) write(6,*) 'finished converting time: ',didt

	return
   99	continue
	write(6,*) 'name: ',name
	stop 'error stop convert_time: cannot parse'
	end

c********************************************************************

