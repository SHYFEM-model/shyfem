c
c $Id: subnsh.f,v 1.54 2010-03-22 15:29:31 georg Exp $
c
c utility routines for shyfem main routine
c
c contents :
c
c subroutine prilog                     prints parameters for main
c subroutine priout(mode)               writes output files
c subroutine pritst(id)                 test output of constants and variables
c subroutine init_3d			sets up 3D vertical arrays
c subroutine nlsh2d(iunit)		read STR  parameter file for FE model
c subroutine rdtitl			reads title section
c
c subroutine impini		initializes parameters for semi-implicit time
c function bimpli(it)		checks if semi-implicit time-step is active
c function getimp		gets weight for semi-implicit time-step
c subroutine setimp(it,aweigh)	sets parameters for semi-implicit time-step
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
c 10.11.2014    ggu     shyfem time management routines to new file subtime.f
c 01.12.2014    ccf     handle new section waves for wave module
c 24.09.2015    ggu     call initialization for irv before reading STR file
c 26.05.2016    ggu     new check for sections: count_sections()
c 16.06.2016    wjm     added check for section nonhyd 
c 11.05.2018    ggu     semi.h deleted and substituted with module
c
c************************************************************

	subroutine prilog

c writes output to terminal or log file

	use shympi

	implicit none

	include 'modules.h'
	include 'femtime.h'
	include 'simul.h'

	character*80 name
        integer nrb,nbc
        integer nkbnd,nbnds
	character*20 aline

	if( .not. shympi_is_master() ) return

        nrb = nkbnd()
        nbc = nbnds()

	call getfnm('runnam',name)
	write(6,*)
	write(6,*) '     Name of run :'
	write(6,*) trim(name)

	write(6,*)
	write(6,*) '     Description of run :'
	write(6,*) trim(descrp)
	write(6,*)

	call dts_format_abs_time(atime0+dtanf,aline)
	write(6,*) '     itanf = ',aline
	call dts_format_abs_time(atime0+dtend,aline)
	write(6,*) '     itend = ',aline
	write(6,*) '     idt =   ',idt
	write(6,*) '     Iterations to go :',nits

	call getfnm('basnam',name)
	write(6,*)
	write(6,*) '     Name of basin :'
	write(6,*) trim(name)

	write(6,*)
	write(6,*) '     Description of basin :'
	call bas_info

	write(6,*) '     Description of boundary values :'
	write(6,*)
	write(6,*) '     nbc,nrb     :',nbc,nrb

	write(6,*)
	write(6,*) '     Values from parameter file :'
	write(6,*)

	call pripar(6)
	call check_parameter_values('prilog')

	call prbnds		!prints boundary info

	!call prexta		!prints extra points

	!call pretsa		!prints extra time series points

	call prflxa		!prints flux sections

	call prvola		!prints flux sections

	call prarea		!prints chezy values

	call prclos		!prints closing sections

c	call proxy		!prints oxygen section

c	call prlgr		!prints float coordinates

	call modules(M_PRINT)

	write(6,*)
	write(6,1030)
	write(6,*)

	return
 1030   format(1x,78('='))
	end

c********************************************************************

	subroutine pritst(id)

c test output of all constants and variables
c
c id    identifier

	use basin

	implicit none

	integer id

	include 'modules.h'

	include 'simul.h'

	write(6,*)
	write(6,*) '============================================='
	write(6,*) '================ test output ================'
	write(6,*) '============================================='
	write(6,1) '================ id = ',id,' ================'
	write(6,*) '============================================='
	write(6,*)

	call bas_info

	write(6,*) '/descrp/'
	write(6,*) descrp

	!call tsexta

	!call tsetsa

	call tsflxa

	call tsvola

	call tsarea

	call tsbnds

	call tsclos

c	call tsoxy	!oxygen

c	write(6,*) '/close/'	!deleted on 28.05.97

	call modules(M_TEST)
	
	write(6,*)
	write(6,*) '============================================='
	write(6,*) '============================================='
	write(6,*) '============================================='
	write(6,*)

	return
    1   format(1x,a,i6,a)
	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine do_init

c to do before time loop

	implicit none

	include 'modules.h'

	call wrboxa
	call wrousa

	end

c********************************************************************

	subroutine do_befor

c to do in time loop before time step

	implicit none

	include 'modules.h'

	call modules(M_BEFOR)

        call tideforce       !tidal potential !ccf
	call adjust_chezy

	end

c********************************************************************

	subroutine do_after

c to do in time loop after time step

	implicit none

	include 'modules.h'

	double precision dtime

	call modules(M_AFTER)

	call get_act_dtime(dtime)

c	call wrouta
	call wrousa
c	call wrexta
	call wrflxa
	call wrvola(dtime)
	call wrboxa

        call resid
        call rmsvel

        call admrst             !restart

c        call tsmed
	call ts_shell

c	call wrnetcdf		!output in netcdf format - not supported

	call custom(dtime)

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine nlsh2d(iunit)
c
c read STR  parameter file for FE model
c
c iunit		unit number of file
c
c revised 20.01.94 by ggu !$$conz - impl. of concentration in bnd(12,.)
c revised 07.04.95 by ggu !$$baroc - impl. of baroclinic salt/temp (21/22)
c revised ...06.97 by ggu !complete revision
c 18.03.1998	ggu	use variable section instead name

	use levels
	use nls

	implicit none

	integer iunit

	include 'modules.h'

c---------------------------------------------------------------
c---------------------------------------------------------------

	character*6 section,extra,last
	logical bdebug
	integer nsc,num,iline
c	integer nrdsec,nrdveci,nrdvecr
	integer nrdsec,nrdvecr
	character*80 vers,aline

	logical hasreadsec

	bdebug = .true.
	bdebug = .false.

        nlv = 0         !is initialized really only in adjust_levels

	nsc = 0
	if(iunit.le.0) goto 63
	last = ' '

	if( bdebug ) write(6,*) 'start reading STR file'

	call nrdini(iunit)
	call mod_irv_initialize

c read loop over sections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do while( nrdsec(section,num,extra) .ne. 0 )

		if( bdebug ) write(6,*) 'new section: ',section,num

		call setsec(section,num)		!remember section
		call count_sections(section)

		nsc = nsc + 1

		if(section.eq.'title') then
			call rdtitl
			if( nsc .ne. 1 ) goto 66
		else if(section.eq.'end') then
			goto 69
		else if(section.eq.'para') then
			call nrdins(section)
		else if(section.eq.'proj') then
			call nrdins(section)
		else if(section.eq.'extra') then
			call rdexta
			!call section_deleted(section,'use section $extts')
		else if(section.eq.'extts') then
			call rdetsa
		else if(section.eq.'area') then
			call rdarea
		else if(section.eq.'name') then
			call nrdins(section)
		else if(section.eq.'bound') then
			call rdbnds(num)
		else if(section.eq.'float') then
			!call rdfloa(nfldin)
			call section_deleted(section,'use section $lagrg')
		else if(section.eq.'close') then
			call rdclos(num)
		else if(section.eq.'flux') then
			call rdflxa
		else if(section.eq.'vol') then
			call rdvola
		else if(section.eq.'wrt') then		!water renewal time
                        call nrdins(section)
		else if(section.eq.'wind') then
			call section_deleted(section,'use wind file')
		else if(section.eq.'oxypar') then	!oxygen
			call nrdins(section)
		else if(section.eq.'oxyarr') then	!oxygen
			!call rdoxy
			call section_deleted(section,'use ecological model')
		else if(section.eq.'bfmsc')then        ! BFM ECO TOOL
                        call nrdins(section)
		else if(section.eq.'levels') then
			call read_hlv
                else if(section.eq.'lagrg')then
                        call nrdins(section)
                else if(section.eq.'sedtr')then         !sediment
                        call readsed
                else if(section.eq.'waves')then         !wave
                        call nrdins(section)
                else if(section.eq.'mudsec')then        !fluid mud
                        call readmud			!ARON
                else if(section.eq.'nonhyd')then        !NH model
                        call nrdins(section)	
                else if(section.eq.'connec')then       !connectivity
                        call nrdins(section)	
		else					!try modules
			call modules(M_READ)
			if( .not. hasreadsec() ) then	!sec has been handled?
				goto 97			! -> no
			end if
		end if

		last = section
	end do		!loop over sections

	call count_sections(' ')	!check section count

	if( bdebug ) write(6,*) 'finished reading STR file'

c end of read %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
	write(6,*) 'no section has been opened'
	call nls_return_line(aline,iline)
	write(6,*) trim(aline),'     line number = ',iline
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

c************************************************************************

	subroutine section_deleted(section,compat)

	implicit none

	character*(*) section,compat

	write(6,*) 'the following section has been removed: '
	write(6,*) section
	write(6,*) 'please remove section from parameter file'
	write(6,*) 'and substitute with a compatible solution'
        if( compat .ne. ' ' ) then
          write(6,*) 'a compatible solution could be: ',compat
        end if

	stop 'error stop section_deleted: section has been removed'
	end

c************************************************************************

        subroutine read_hlv

	use levels
	use basin, only : nkn,nel,ngr,mbw
        use nls

        implicit none

        integer n,nlvddi

	call nls_init_section

        n = nls_read_vector()
        !call levels_init(nkn,nel,n)
	call levels_hlv_init(n)
	nlv = n
	call levels_get_dimension(nlvddi)
	if( n > nlvddi ) then
	  write(6,*) 'nlv,nlvddi: ',nlv,nlvddi
	  stop 'error stop read_hlv: dimension error'
	end if
        call nls_copy_real_vect(n,hlv)

	call nls_finish_section

        end subroutine read_hlv

c************************************************************************

	subroutine rdtitl

c reads title section

	implicit none

	include 'simul.h'

	character*80 line,extra

	integer nrdlin,nrdsec
	integer num

	if( nrdlin(line) .eq. 0 ) goto 65
	line = adjustl(line)
	descrp=line
	call putfnm('title',line)

	if( nrdlin(line) .eq. 0 ) goto 65
	line = adjustl(line)
	call putfnm('runnam',line)

	if( nrdlin(line) .eq. 0 ) goto 65
	line = adjustl(line)
	call putfnm('basnam',line)

	!if( nrdlin(line) .gt. 0 ) goto 65
	if( nrdsec(line,num,extra) .eq. 0 ) goto 65
	if( line .ne. 'end' ) goto 64

	return
   64	continue
	write(6,*) 'error in section $title'
	stop 'error stop rdtitl: no end found'
   65	continue
	write(6,*) 'error in section $title'
	stop 'error stop rdtitl: cannot read title section'
	end

c**********************************************************************
c**********************************************************************
c**********************************************************************
c
c routines for handling semi-implicit time-step
c
c the only routines that should be necessary to be called are
c setimp(it,weight) and getazam(az,am)
c
c setimp sets the implicit parameter until time it to weight
c getazam returns az,am with the actual weight
c
c usage: call setimp in a program that would like to change the
c	 weight for a limited time (closing sections etc...)
c	 call getazam when az,am are needed
c	 (getaz if only az is needed)
c
c**********************************************************************
c**********************************************************************
c**********************************************************************

!======================================================================
	module semi_implicit
!======================================================================

	double precision, save :: 	dtimpl
	real, save :: 			weight = 0.5
	logical, save :: 		binit = .false.

!======================================================================
	end module semi_implicit
!======================================================================

c**********************************************************************

	subroutine getaz(azpar)

c returns actual az

	use semi_implicit

	implicit none

	real azpar

	include 'femtime.h'

	real ampar
	real getpar

	azpar=getpar('azpar')
	ampar=0.			!dummy

	call changeimp(it,azpar,ampar)

	end

c**********************************************************************

	subroutine getazam(azpar,ampar)

c returns actual az,am

	use semi_implicit

	implicit none

	real azpar
	real ampar

	include 'femtime.h'

	real getpar

	azpar=getpar('azpar')
	ampar=getpar('ampar')

	call changeimp(it,azpar,ampar)

	end

c**********************************************************************

	subroutine changeimp(it,azpar,ampar)

c changes parameters for semi-implicit time-step if necessary

	use semi_implicit

	implicit none

	integer it
	real azpar,ampar

	if( binit .and. it .le. dtimpl ) then
	  azpar = weight
	  ampar = weight
	end if

	end

c**********************************************************************

	subroutine setimp(it,aweigh)

c sets parameters for semi-implicit time-step

	use semi_implicit

	implicit none

	integer it
	real aweigh

	dtimpl = it
	weight = aweigh

	write(6,*) 'implicit parameters changed: ',dtimpl,weight

	end

c**********************************************************************

	function getimp()

c gets weight for semi-implicit time-step

	use semi_implicit

	implicit none

	real getimp

	getimp = weight

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

        subroutine getinfo(iunit)

c gets unit of info file

        implicit none

        integer iunit

        integer ifemop

        integer, save :: iu = 0

        if( iu .le. 0 ) then
          iu = ifemop('.inf','formatted','new')
          if( iu .le. 0 ) then
            write(6,*) 'error in opening info file'
            stop 'error stop getinfo'
          end if
        end if

        iunit = iu

        end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine setup_omp_parallel

	use shympi

	implicit none

	logical bm
	integer n,nomp
	real getpar
	logical openmp_is_parallel

	call openmp_get_max_threads(n)
	nomp = nint(getpar('nomp'))
	if( nomp > 0 ) then
	  nomp = min(nomp,n)
	else if( nomp == 0 ) then
	  nomp = 1
	else	!nomp < 0 ... use all threads
	  nomp = n
	end if
	call openmp_set_num_threads(nomp)
	call putpar('nomp',float(nomp))

	if( .not. shympi_is_master() ) return

	write(6,*) 'start of setup of parallel OMP threads'

	if( openmp_is_parallel() ) then
	  write(6,*) 'the program can run in OMP parallel mode'
	else
	  write(6,*) 'the program cannot run in OMP parallel mode'
	end if
	  
	write(6,*) 'maximum available OMP threads: ',n
	write(6,*) 'for simulation used OMP threads: ',nomp
	write(6,*) 'end of setup of parallel OMP threads'

	end

c********************************************************************
c********************************************************************
c********************************************************************

        subroutine total_energy

c writes info on total energy to info file

	use shympi

	implicit none

	real kenergy,penergy,tenergy,ksurf,paux
	real energy(3)
	character*20 aline
	logical debug

	integer, save :: iuinfo = 0

	debug = .false.

	if( iuinfo .eq. 0 ) then
          call getinfo(iuinfo)  !unit number of info file
	end if

	call energ3d(kenergy,penergy,ksurf,-1)
	!call energ3d(kenergy,penergy,ksurf,0)

	if( debug ) write(6,*) 'penergy: ',my_id,penergy
	kenergy = shympi_sum(kenergy)
	penergy = shympi_sum(penergy)
	ksurf = shympi_sum(ksurf)
	if( debug ) write(6,*) 'penergy total: ',my_id,penergy

	tenergy = kenergy + penergy

	if(shympi_is_master()) then
	  call get_act_timeline(aline)
	  write(iuinfo,1000) ' energy: ',aline
     +				,kenergy,penergy,tenergy,ksurf
 1000	  format(a,a20,4e12.4)
	end if

	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine count_sections(section)

	implicit none

	character*(*) section

	integer, parameter :: ndim = 100
	character*80, save :: sections(ndim) = ' '
	integer, save :: count(ndim) = 0
	integer, save :: ntot = 0

	integer i,ierr,ic

!----------------------------------------------------------
! check sections
!----------------------------------------------------------

	if( section == ' ' ) then	!check sections
	  ierr = 0
	  do i=1,ntot
	    ic = count(i)
	    if( ic > 1 ) then
	      if( sections(i) /= 'bound' ) cycle
	      if( sections(i) /= 'close' ) cycle
	      write(6,*) 'section ',trim(sections(i)),' - count = ',ic
	      ierr = 1
	    end if
	  end do
	  if( ierr > 0 ) then
	    stop 'error stop count_sections: not unique sections'
	  end if
	  return
	end if

!----------------------------------------------------------
! count section
!----------------------------------------------------------

	do i=1,ntot
	  if( section == sections(i) ) then
	    count(i) = count(i) + 1
	    return
	  end if
	end do

	ntot = ntot + 1
	if( ntot > ndim ) then
	  stop 'error stop count_sections: ndim'
	end if

!----------------------------------------------------------
! insert new section
!----------------------------------------------------------

	sections(ntot) = section
	count(ntot) = 1

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	end

c********************************************************************
