!
! $Id: subtime.f,v 1.54 2010-03-22 15:29:31 georg Exp $
!
! time management routines
!
! contents :
!
!
! revision log :
!
! 23.09.1997	ggu	boundn deleted -> no access to data structure
! 20.03.1998	ggu	minor changes to priout
! 29.04.1998	ggu	new module for semi-implicit time-step
! 07.05.1998	ggu	check for error on return of nrdvecr
! 19.06.1998	ggu	version number is character
! 22.01.1999	ggu	oxygen section added
! 26.01.1999	ggu	new comp3d added
! 11.08.1999	ggu	new compatibility array hlhv initialized
! 19.11.1999	ggu	new routines for section vol
! 20.01.2000    ggu     common block /dimdim/ eliminated
! 04.02.2000    ggu     no priout, dobefor/after, pritime, endtime
! 15.05.2000    ggu     hm3v substituted
! 26.05.2000    ggu     copright statement adjourned
! 21.11.2001    ggu     routines to handle advective index (aix)
! 27.11.2001    ggu     routine to handle info file (getinfo)
! 11.10.2002    ggu     aix routines deleted
! 07.02.2003    ggu     routine added: changeimp, getaz; deleted getaza
! 10.08.2003    ggu     call adjust_chezy instead sp135r
! 14.08.2003    ggu     femver transfered to subver, not called in nlsh2d
! 20.08.2003    ggu     tsmed substituted by ts_shell
! 01.09.2003    ggu     call wrousa
! 03.09.2004    ggu     call admrst, comp3d renamed to init_3d (not used)
! 03.09.2004    ggu     nlv, hlv initialized in nlsh2d (FIXME)
! 28.09.2004    ggu     read lagrangian section
! 01.12.2004    ggu     new routine set_timestep for variable time step
! 17.01.2005    ggu     get_stab_index to newcon.f, error stop in set_timestep
! 14.03.2005    ggu     syncronize idt with end of simulation (set_timestep)
! 07.11.2005    ggu     handle new section sedtr for sediments
! 23.03.2006    ggu     changed time step to double precision
! 23.05.2007    ggu     recall variable time step pars at every time step
! 02.10.2007    ggu     bug fix in set_timestep for very small rindex
! 10.04.2008    ccf     output in netcdf format
! 28.04.2008    ggu     in set_timestep new call to advect_stability()
! 03.09.2008    ggu     in nlsh2d different error message
! 20.11.2008    ggu     init_3d deleted, nlv initialized to 0
! 18.11.2009    ggu     new format in pritime (write also time step)
! 22.02.2010    ggu     new call to hydro_stability to compute time step
! 22.02.2010    ccf     new routine for tidal pot. (tideforc), locaus deleted
! 26.02.2010    ggu     in set_timestep compute and write ri with old dt
! 22.03.2010    ggu     some comments for better readability
! 29.04.2010    ggu     new routine set_output_frequency() ... not finished
! 04.05.2010    ggu     shell to compute energy
! 22.02.2011    ggu     in pritime() new write to terminal
! 20.05.2011    ggu     changes in set_timestep(), element removal, idtmin
! 31.05.2011    ggu     changes for BFM
! 01.06.2011    ggu     idtmin introduced
! 12.07.2011    ggu     new routine next_output(), revised set_output_frequency
! 14.07.2011    ggu     new routines for original time step
! 13.09.2011    ggu     better error check, rdtitl() more robust
! 23.01.2012    ggu     new section "proj"
! 24.01.2012    ggu     new routine setup_parallel()
! 10.02.2012    ggu     new routines to initialize and access time common block
! 05.03.2014    ggu     code prepared to repeat time step (irepeat) - not ready
! 05.03.2014    ggu     new routines get_last/first_time()
! 10.04.2014    ccf     new section "wrt" for water renewal time
! 29.10.2014    ggu     do_() routines transfered from newpri.f
! 10.11.2014    ggu     time management routines transfered to this file
! 23.09.2015    ggu     new routine convert_time_d() for double
! 24.09.2015    ggu     new file for double precision version
! 20.10.2015    ggu     new routines to set/get unit
! 04.11.2015    ggu     allow for initial output in adjust_itmidt()
!
!************************************************************
!
!**********************************************************************
!**********************************************************************
!**********************************************************************
!---------------------------------------------------------------------
        module output
!---------------------------------------------------------------------
        contains
!---------------------------------------------------------------------

	subroutine adjust_itmidt(itmout,idtout,itout)

! sets-up output frequency and first output

	implicit none

	include 'femtime.h'

	integer itmout		!minimum time for output
	integer idtout		!time step for output
	integer itout		!first output

	logical binit

	binit = ( itmout /= -1 )		!output initial time?

	if( itmout .eq. -1 ) itmout = itanf
	if( itmout .lt. itanf ) itmout = itanf

	itout = itmout
	if( itmout .eq. itanf .and. .not. binit ) itout = itout + idtout
	if( itout .gt. itend .and. idtout .gt. 0 ) idtout = 0

	end

!********************************************************************

	subroutine set_output_frequency(itmout,idtout,ia_out)

! sets-up array for output frequency

	implicit none

	integer itmout		!minimum time for output
	integer idtout		!time step for output
	integer ia_out(4)	!array where info is stored

	integer itout

	call adjust_itmidt(itmout,idtout,itout)

	ia_out(1) = idtout	! time step of output
	ia_out(2) = itmout	! first output
	ia_out(3) = itout	! next output
	ia_out(4) = 0		! unit (optional)

	end

!********************************************************************

	subroutine assure_initial_output(ia_out)

! makes sure that output will be done also for it == itanf

	implicit none

	integer ia_out(4)

	integer itmout,itout

	itmout = ia_out(2)
	itout  = ia_out(3)

	itout = itmout

	ia_out(3) = itout

	end

!********************************************************************

	subroutine increase_output(ia_out)

! makes sure that itout > itmout

	implicit none

	integer ia_out(4)

	include 'femtime.h'

	integer idtout,itmout,itout

	idtout = ia_out(1)
	itmout = ia_out(2)
	itout  = ia_out(3)

	if( itout > itmout ) return

	itout = itmout + idtout
	if( itout .gt. itend ) idtout = 0

	ia_out(1) = idtout	! time step of output
	ia_out(3) = itout	! next output

	end

!********************************************************************

	function is_over_output(ia_out)

! checks if output phase has started (it > itmout)

	implicit none

	logical is_over_output
	integer ia_out(4)

	include 'femtime.h'

	is_over_output = it > ia_out(2)

	end

!********************************************************************

	function is_in_output(ia_out)

! checks if we arrived at output phase (it >= itmout)

	implicit none

	logical is_in_output
	integer ia_out(4)

	include 'femtime.h'

	is_in_output = it >= ia_out(2)

	end

!********************************************************************

	function has_output(ia_out)

! checks if variable has any output at all

	implicit none

	logical has_output
	integer ia_out(4)

	has_output = ia_out(1) > 0	!idtout > 0

	end

!********************************************************************

	function next_output(ia_out)

! checks if time has come for output

	implicit none

	logical next_output
	integer ia_out(4)

	include 'femtime.h'

	integer idtout,itout

	next_output = .false.
	idtout = ia_out(1)
	itout  = ia_out(3)

	if( idtout .le. 0 ) return
	if( itout .gt. it ) return

	do while( itout .le. it )
	  itout = itout + idtout
	end do

	ia_out(3) = itout
	next_output = .true.

	end

!********************************************************************

	subroutine info_output(ia_out)

! writes info on ia_output

	implicit none

	integer ia_out(4)

	include 'femtime.h'

	write(6,*) '------ info_output start------'
	write(6,*) ia_out
	write(6,*) it
	write(6,*) has_output(ia_out)
	write(6,*) next_output(ia_out)
	write(6,*) '------ info_output end ------'

	end

!********************************************************************

	subroutine set_unit_output(ia_out,iunit)

	implicit none

	integer ia_out(4)
	integer iunit

	ia_out(4) = iunit

	end

!********************************************************************

	subroutine get_unit_output(ia_out,iunit)

	implicit none

	integer ia_out(4)
	integer iunit

	iunit = ia_out(4)

	end

!********************************************************************

	subroutine init_output(itmname,idtname,ia_out)

! gets time values and transforms them

	implicit none

	character*(*) itmname,idtname	!names to parse
	integer ia_out(4)		!array with time information

	integer itmout,idtout

	call convert_date(itmname,itmout)
	call convert_time(idtname,idtout)

	call set_output_frequency(itmout,idtout,ia_out)

	end

!********************************************************************

	subroutine convert_date(name,it)

! converts date to relative time
        use para
        use dts

	implicit none

	character*(*) name
	integer it

	integer ierr
	double precision dit
	character*30 text
	logical bdebug

	bdebug = .true.
	bdebug = .false.

	call getfnm(name,text)

	if( text .ne. ' ' ) then
	  call dtsgunf(it,text,ierr)
	  if( ierr .ne. 0 ) goto 99
	  if( bdebug ) then
	    write(6,*) 'time as string found'
	    write(6,*) name
	    write(6,*) text
	    write(6,*) it
	  end if
	  dit = it
	  call dputpar(name,dit)
	else
	  it = nint(dgetpar(name))
	end if

	return
   99	continue
	write(6,*) 'name: ',name
	stop 'error stop convert_date: cannot parse'
	end

!********************************************************************

	subroutine convert_time(name,idt)

! converts time period to relative time difference

        use para
        use dts

	implicit none

	character*(*) name
	integer idt

	integer ierr
	double precision didt
	character*40 text
	logical bdebug

	bdebug = .true.
	bdebug = .false.

	call getfnm(name,text)

	if( bdebug ) write(6,*) 'converting time for ',name

	if( text .ne. ' ' ) then
	  call dtstimespan(idt,text,ierr)
	  if( ierr .ne. 0 ) goto 99
	  if( bdebug ) then
	    write(6,*) 'time span as string found'
	    write(6,*) name
	    write(6,*) text
	    write(6,*) idt
	  end if
	  didt = idt
	  call dputpar(name,didt)
	else
	  idt = nint(dgetpar(name))
	end if

	if( bdebug ) write(6,*) 'finished converting time: ',idt

	return
   99	continue
	write(6,*) 'name: ',name
	stop 'error stop convert_time: cannot parse'
	end

!********************************************************************

!---------------------------------------------------------------------
      end module output
!---------------------------------------------------------------------
