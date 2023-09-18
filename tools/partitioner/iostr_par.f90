!
! $Id: subnsh.f,v 1.54 2010-03-22 15:29:31 georg Exp $
!
! utility routines for shyfem main routine
!
! contents :
!
! subroutine prilog                     prints parameters for main
! subroutine priout(mode)               writes output files
! subroutine pritst(id)                 test output of constants and variables
! subroutine init_3d			sets up 3D vertical arrays
! subroutine nlsh2d(iunit)		read STR  parameter file for FE model
! subroutine rdtitl			reads title section
!
! subroutine impini		initializes parameters for semi-implicit time
! function bimpli(it)		checks if semi-implicit time-step is active
! function getimp		gets weight for semi-implicit time-step
! subroutine setimp(it,aweigh)	sets parameters for semi-implicit time-step
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
! 23.03.2006    ggu     changed time step to real
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
! 10.11.2014    ggu     shyfem time management routines to new file subtime.f
! 01.12.2014    ccf     handle new section waves for wave module
! 24.09.2015    ggu     call initialization for irv before reading STR file
!
!************************************************************
!-----------------------------------------------------------------
        module iostr_par
!-----------------------------------------------------------------
        contains
!-----------------------------------------------------------------

!********************************************************************
!********************************************************************
!********************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine nlsh2d(iunit)
!
! read STR  parameter file for FE model
!
! iunit		unit number of file
!
! revised 20.01.94 by ggu !$$conz - impl. of concentration in bnd(12,.)
! revised 07.04.95 by ggu !$$baroc - impl. of baroclinic salt/temp (21/22)
! revised ...06.97 by ggu !complete revision
! 18.03.1998	ggu	use variable section instead name

	use levels
        use nls
        use bnd_admin_par

	implicit none

	integer iunit

	include 'param.h'

!---------------------------------------------------------------
!---------------------------------------------------------------

	character*6 section,extr,last
	logical bdebug
	integer nsc,num
	character*80 vers

	bdebug = .true.
	bdebug = .false.

        nlv = 0         !is initialized really only in adjust_levels

	nsc = 0
	if(iunit.le.0) goto 63
	last = ' '

	if( bdebug ) write(6,*) 'start reading STR file'

	call nrdini(iunit)

! read loop over sections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do while( nrdsec(section,num,extr) .ne. 0 )

		if( bdebug ) write(6,*) 'new section: ',section,num

		call setsec(section,num)		!remember section

		nsc = nsc + 1

		if(section.eq.'title') then
			call rdtitl
			if( nsc .ne. 1 ) goto 66
		else if(section.eq.'end') then
			goto 69
		else if(section.eq.'para') then
			call nrdins(section)
		else if(section.eq.'name') then
			call nrdins(section)
		else if(section.eq.'bound') then
			call rdbnds(num)
		else if(section.eq.'levels') then
			call read_hlv
		else					!try modules
			!call modules(M_READ)
			if( .not. hasreadsec() ) then	!sec has been handled?
				goto 97			! -> no
			end if
		end if

		last = section
	end do		!loop over sections

	if( bdebug ) write(6,*) 'finished reading STR file'

! end of read %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	return
   63	continue
	write(6,*) 'Cannot read STR file on unit : ',iunit
	stop 'error stop : nlsh2d'
   66	continue
	write(6,*) 'section $title must be first section'
	stop 'error stop : nlsh2d'
   69	continue
	write(6,*) 'section $end is not allowed on its own'
	write(6,*) 'last section found was: ',last
	stop 'error stop : nlsh2d'
   77	continue
	if( nlv .eq. -1 ) then
	  write(6,*) 'read error in section $levels'
	  write(6,*) 'nlv,nlvdi: ',nlv,nlvdi
	  stop 'error stop nlsh2d: read error'
	else
	  write(6,*) 'dimension error in section $levels'
	  write(6,*) 'nlvdi = ',nlvdi,'   number of data read = ',-nlv
	  stop 'error stop nlsh2d: dimension error'
	end if
   97	continue
	write(6,*) 'Cannot handle section : ',section
	stop 'error stop nlsh2d: no such section'
	end

!************************************************************************

        subroutine read_hlv

	use levels
	use basin, only : nkn,nel,ngr,mbw
        use nls

        implicit none

	include 'param.h'

        integer n,nlvddi

        n = nls_read_vector()
	call levels_hlv_init(n)
	nlv = n
	call levels_get_dimension(nlvddi)
	if( n > nlvddi ) then
	  write(6,*) 'nlv,nlvddi: ',nlv,nlvddi
	  stop 'error stop read_hlv: dimension error'
	end if
        call nls_copy_real_vect(n,hlv)

        end subroutine read_hlv

!************************************************************************

	subroutine rdtitl

! reads title section

        use para
        use utility
        use nls

	implicit none

	include 'param.h'
	include 'simul.h'

	character*80 line,extr

	integer num

	if( nrdlin(line) .eq. 0 ) goto 65
	call triml(line)
	descrp=line
	call putfnm('title',line)

	if( nrdlin(line) .eq. 0 ) goto 65
	call triml(line)
	call putfnm('runnam',line)

	if( nrdlin(line) .eq. 0 ) goto 65
	call triml(line)
	call putfnm('basnam',line)

	!if( nrdlin(line) .gt. 0 ) goto 65
	if( nrdsec(line,num,extr) .eq. 0 ) goto 65
	if( line .ne. 'end' ) goto 64

	return
   64	continue
	write(6,*) 'error in section $title'
	stop 'error stop rdtitl: no end found'
   65	continue
	write(6,*) 'error in section $title'
	stop 'error stop rdtitl: cannot read title section'
	end

!**********************************************************************
!**********************************************************************
!**********************************************************************
!
! routines for handling semi-implicit time-step
!
! the only routines that should be necessary to be called are
! setimp(it,weight) and getazam(az,am)
!
! setimp sets the implicit parameter until time it to weight
! getazam returns az,am with the actual weight
!
! usage: call setimp in a program that would like to change the
!	 weight for a limited time (closing sections etc...)
!	 call getazam when az,am are needed
!	 (getaz if only az is needed)
!
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!********************************************************************

!-----------------------------------------------------------------
        end module iostr_par
!-----------------------------------------------------------------
