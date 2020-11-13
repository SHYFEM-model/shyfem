
!--------------------------------------------------------------------------
!
!    Copyright (C) 1997-2012,2014-2020  Georg Umgiesser
!    Copyright (C) 2008,2010,2014  Christian Ferrarin
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
c 20.01.2000	ggu	common block /dimdim/ eliminated
c 04.02.2000	ggu	no priout, dobefor/after, pritime, endtime
c 15.05.2000	ggu	hm3v substituted
c 26.05.2000	ggu	copright statement adjourned
c 21.11.2001	ggu	routines to handle advective index (aix)
c 27.11.2001	ggu	routine to handle info file (getinfo)
c 11.10.2002	ggu	aix routines deleted
c 07.02.2003	ggu	routine added: changeimp, getaz; deleted getaza
c 10.08.2003	ggu	call adjust_chezy instead sp135r
c 14.08.2003	ggu	femver transfered to subver, not called in nlsh2d
c 20.08.2003	ggu	tsmed substituted by ts_shell
c 01.09.2003	ggu	call wrousa
c 03.09.2004	ggu	call admrst, comp3d renamed to init_3d (not used)
c 03.09.2004	ggu	nlv, hlv initialized in nlsh2d (FIXME)
c 28.09.2004	ggu	read lagrangian section
c 01.12.2004	ggu	new routine set_timestep for variable time step
c 17.01.2005	ggu	get_stab_index to newcon.f, error stop in set_timestep
c 14.03.2005	ggu	syncronize idt with end of simulation (set_timestep)
c 07.11.2005	ggu	handle new section sedtr for sediments
c 23.03.2006	ggu	changed time step to real
c 23.05.2007	ggu	recall variable time step pars at every time step
c 02.10.2007	ggu	bug fix in set_timestep for very small rindex
c 10.04.2008	ccf	output in netcdf format
c 28.04.2008	ggu	in set_timestep new call to advect_stability()
c 03.09.2008	ggu	in nlsh2d different error message
c 20.11.2008	ggu	init_3d deleted, nlv initialized to 0
c 18.11.2009	ggu	new format in pritime (write also time step)
c 22.02.2010	ggu	new call to hydro_stability to compute time step
c 22.02.2010	ccf	new routine for tidal pot. (tideforc), locaus deleted
c 26.02.2010	ggu	in set_timestep compute and write ri with old dt
c 22.03.2010	ggu	some comments for better readability
c 29.04.2010	ggu	new routine set_output_frequency() ... not finished
c 04.05.2010	ggu	shell to compute energy
c 22.02.2011	ggu	in pritime() new write to terminal
c 20.05.2011	ggu	changes in set_timestep(), element removal, idtmin
c 31.05.2011	ggu	changes for BFM
c 01.06.2011	ggu	idtmin introduced
c 12.07.2011	ggu	new routine next_output(), revised set_output_frequency
c 14.07.2011	ggu	new routines for original time step
c 13.09.2011	ggu	better error check, rdtitl() more robust
c 23.01.2012	ggu	new section "proj"
c 24.01.2012	ggu	new routine setup_parallel()
c 10.02.2012	ggu	new routines to initialize and access time common block
c 05.03.2014	ggu	code prepared to repeat time step (irepeat) - not ready
c 05.03.2014	ggu	new routines get_last/first_time()
c 10.04.2014	ccf	new section "wrt" for water renewal time
c 29.10.2014	ggu	do_() routines transfered from newpri.f
c 10.11.2014	ggu	time management routines transfered to this file
c 26.11.2014	ggu	changed VERS_7_0_7
c 05.12.2014	ggu	changed VERS_7_0_8
c 12.12.2014	ggu	changed VERS_7_0_9
c 19.12.2014	ggu	accept date also as string
c 23.12.2014	ggu	fractional time step introduced
c 07.01.2015	ggu	fractional time step without rounding (itsplt=3)
c 26.02.2015	ggu	changed VERS_7_1_5
c 30.04.2015	ggu	changed VERS_7_1_9
c 23.09.2015	ggu	time step is now working with dt as double
c 10.10.2015	ggu	use bsync as global to check for syncronization
c 23.09.2016	ggu	cleaned set_timestep()
c 30.09.2016	ggu	changed VERS_7_5_18
c 20.10.2017	ggu	new get_absolute_act_time(),get_absolute_ref_time()
c 04.11.2017	ggu	changed VERS_7_5_34
c 05.12.2017	ggu	changed VERS_7_5_39
c 22.02.2018	ggu	changed VERS_7_5_42
c 23.02.2018	ggu	most parts converted from int to double
c 29.03.2018	ggu	bug fix for syncronization step
c 03.04.2018	ggu	changed VERS_7_5_43
c 03.04.2018	ggu	changed VERS_7_5_44
c 13.04.2018	ggu	hydro_stability includes explicit gravity wave
c 13.04.2018	ggu	set_timestep now is working with mpi
c 16.04.2018	ggu	write warning if time step is over recommended one
c 09.02.2019	ggu	bug fix for syncronization of last time step
c 16.02.2019	ggu	changed VERS_7_5_60
c 21.05.2019	ggu	changed VERS_7_5_62
c 15.09.2019	ggu	small changes to account for synchorization time step
c 08.02.2020	ggu	utilities in this new file
c 16.02.2020	ggu	itunit eliminated
c 16.02.2020	ggu	new routines get_time_iterations(), get_ddt()
c 03.04.2020	ggu	new routine get_real_time()
c
c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine is_time_first(bfirst)

c true if in initialization phase

	implicit none

	logical bfirst

	include 'femtime.h'

	bfirst = t_act .eq. dtanf

	end

c**********************************************************************

	subroutine is_time_last(blast)

c true if in last time step

	implicit none

	logical blast

	include 'femtime.h'

	blast = t_act .eq. dtend

	end

c**********************************************************************

        subroutine get_act_dtime(dtact)

c returns actual time (double)

        implicit none

	double precision dtact

	include 'femtime.h'

	dtact = t_act

	end

c**********************************************************************

        subroutine get_timeline(dtime,aline)

c returns time as string

        implicit none

	double precision dtime
	character*(*) aline

	include 'femtime.h'
	double precision atime

	atime = atime0 + dtime
	call dts_format_abs_time(atime,aline)

	end

c**********************************************************************

        subroutine get_act_timeline(aline)

c returns actual time as string

        implicit none

	character*(*) aline

	include 'femtime.h'

	aline = aline_act

	end

c**********************************************************************

        subroutine get_absolute_act_time(atime)

c returns actual time

        implicit none

	double precision atime

	include 'femtime.h'

	atime = t_act + atime0

	end

c**********************************************************************

        subroutine get_absolute_ref_time(atime_ref)

c returns actual time

        implicit none

	double precision atime_ref

	include 'femtime.h'

	atime_ref = atime0

	end

c**********************************************************************

        subroutine get_passed_dtime(dtime)

c returns time passed since start of simulation

        implicit none

	double precision dtime

	include 'femtime.h'

	dtime = t_act - dtanf

	end

c**********************************************************************

        subroutine get_first_dtime(dtime)

c returns first (initial) time

        implicit none

	double precision dtime

	include 'femtime.h'

	dtime = dtanf

	end

c**********************************************************************

        subroutine get_last_dtime(dtime)

c returns end time

        implicit none

	double precision dtime

	include 'femtime.h'

	dtime = dtend

	end

c**********************************************************************

        subroutine get_ddt(ddt)

c returns time step (in seconds - double version)

        implicit none

	double precision ddt		!time step (return)

	include 'femtime.h'

	ddt = dt_act

	end

c**********************************************************************

        subroutine get_timestep(dt)

c returns time step (in seconds - real version)

        implicit none

	real dt		!time step (return)

	include 'femtime.h'

	dt = dt_act

	end

c**********************************************************************

        subroutine get_orig_timestep(dt)

c returns original real time step (in real seconds)

        implicit none

	real dt		!time step (return)

	include 'femtime.h'

	dt = idtorig

	end

c**********************************************************************

        subroutine get_time_iterations(nit_done,nit_todo)

c returns iterations, already done and still to do

        implicit none

	integer nit_done
	integer nit_todo

	include 'femtime.h'

	nit_done = niter
	nit_todo = nits

	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine convert_to_dtime(aline,dtime)

	implicit none

	character*20 aline
	double precision dtime

	double precision atime,atime0

	call convert_to_atime(aline,atime)

        call get_absolute_ref_time(atime0)
        dtime = atime - atime0

	end subroutine convert_to_dtime

c********************************************************************

	subroutine convert_to_atime(aline,atime)

	use iso8601

	implicit none

	character*20 aline
	double precision atime

	integer date,time,ierr

        call string2date(aline,date,time,ierr)
        if( ierr /= 0 ) stop 'error converting date'
        call dts_to_abs_time(date,time,atime)

	end subroutine convert_to_atime

c********************************************************************
c********************************************************************
c********************************************************************

        subroutine get_real_time(atime,aline)

! returns real time as string and absolute time

        implicit none

        double precision atime
        character*20 aline

        character*20 date,time,zone
        integer values(8)

        call date_and_time(date,time,zone,values)

        date = date(1:4)//'-'//date(5:6)//'-'//date(7:8)
        time = time(1:2)//':'//time(3:4)//':'//time(5:6)

        aline = trim(date)//'::'//trim(time)
	call convert_to_atime(aline,atime)

        end

c********************************************************************
c********************************************************************
c********************************************************************

