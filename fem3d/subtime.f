c
c $Id: subtime.f,v 1.54 2010-03-22 15:29:31 georg Exp $
c
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
c 19.12.2014    ggu     accept date also as string
c 23.12.2014    ggu     fractional time step introduced
c 07.01.2015    ggu     fractional time step without rounding (itsplt=3)
c 23.09.2015    ggu     time step is now working with dt as double
c 10.10.2015    ggu     use bsync as global to check for syncronization
c 23.09.2016    ggu     cleaned set_timestep()
c
c************************************************************
c
c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine print_time

c prints time after time step

	implicit none

	include 'femtime.h'

        integer nit1,nit2,naver
	integer idtfrac,i
        real perc,dt

	integer year,month,day,hour,min,sec
	integer, save :: isplit
	double precision daux

	character*20 line
	character*9 frac
	double precision dgetpar
	logical dts_has_date

	integer, save :: icall = 0

c---------------------------------------------------------------
c set parameters and compute percentage of simulation
c---------------------------------------------------------------

	if( icall .eq. 0 ) then
          isplit = nint(dgetpar('itsplt'))
	end if

        naver = 20
        naver = 0

        !perc = (100.*(it-itanf))/(itend-itanf)
        perc = (100.*(t_act-itanf))/(itend-itanf)

c---------------------------------------------------------------
c compute total number of iterations
c---------------------------------------------------------------

	nit1 = 0
	idtfrac = 0

	if( bsync ) then	!syncronization - do not count
	  !
	else if( idt .gt. 0 ) then
          nit1 = niter + (itend-it)/idt
	else if( dt_act > 0 ) then
          daux = (itend-t_act)/dt_act
	  if( daux > 1000000000. ) then
	    write(6,*) '******************************************'
	    write(6,*) '******************************************'
	    write(6,*) '******************************************'
	    write(6,*) t_act,dt_act,itanf,itend
	    write(6,*) daux
	    write(6,*) niter,bsync
	    write(6,*) '******************************************'
	    write(6,*) '******************************************'
	    write(6,*) '******************************************'
	    stop 'error stop print_time: internal error'
	  else
            nit1 = niter + nint((itend-t_act)/dt_act)
	    idtfrac = nint(1./dt_act)
	  end if
	else
	  write(6,*) 'idt,dt_act: ',idt,dt_act
	  write(6,*) 'warning: time step was 0'
	  stop 'error stop print_time: 0 time step'
	end if

	nit2 = nit1
	if( it .gt. itanf ) then
          nit2 = nint(niter*( 1 + (itend-t_act)/(t_act-itanf)))
	end if

        nits = nit2
        if( naver .gt. 0 ) nits = ( 1*nit1 + (naver-1)*nit2 ) / naver

c---------------------------------------------------------------
c write to terminal
c---------------------------------------------------------------

	if( dts_has_date() ) then
	  call dtsgf(it,line)
	else
	  line = ' '
	end if
	if( mod(icall,50) .eq. 0 ) write(6,1003)
	if( .true. ) then
	  !write(6,*) isplit,idtorig,dt_act
	  if( isplit == 3 .or. idtorig == 0 ) then
	    dt = dt_act
            write(6,1007) it,line,dt,niter,nits,perc
	  else if( idtfrac == 0 ) then
            write(6,1005) it,line,idt,niter,nits,perc
	  else
	    frac = ' '
	    write(frac,'(i9)') idtfrac
	    do i=9,1,-1
	      if( frac(i:i) == ' ' ) exit
	    end do
	    frac(i-1:i) = '1/'
            write(6,1006) it,line,frac,niter,nits,perc
	  end if
!          write(6,1002) it,year,month,day,hour,min,sec
!     +			,idt,niter,nits,perc
	else
          write(6,1001) it,idt,niter,nits,perc
	end if

	icall = icall + 1

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	return
 1000   format(' time =',i10,'   iterations =',i6,' / ',i6,f9.2,' %')
 1001   format(' time =',i12,'    dt =',i5,'    iterations ='
     +                 ,i8,' /',i8,f10.2,' %')
 1002   format(i12,i9,5i3,i9,i8,' /',i8,f10.2,' %')
 1003   format(8x,'time',13x,'date',14x,'dt',8x,'iterations'
     +              ,5x,'percent')
 1005   format(i12,3x,a20,1x,i9,i8,' /',i8,f10.2,' %')
 1006   format(i12,3x,a20,1x,a9,i8,' /',i8,f10.2,' %')
 1007   format(i12,3x,a20,1x,f9.2,i8,' /',i8,f10.2,' %')
	end

c********************************************************************

	subroutine print_end_time

c prints stats after last time step

	implicit none

	include 'femtime.h'

	write(6,1035) it,niter
 1035   format(' program stop at time =',i10,' seconds'/
     +         ' iterations = ',i10)

	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine setup_time

c setup and check time parameters

	implicit none

	include 'femtime.h'

	double precision didt

	double precision dgetpar

	call convert_date('itanf',itanf)
	call convert_date('itend',itend)
	call convert_time_d('idt',didt)

	if( didt .le. 0 .or. itanf+didt .gt. itend ) then
	   write(6,*) 'Error in compulsory time parameters'
	   write(6,*) 'itanf,itend,idt :',itanf,itend,didt
	   stop 'error stop : cktime'
	end if

	niter = 0
	it = itanf
	nits = (itend-itanf) / didt

	t_act = it
	dt_act = didt
	dt_orig = didt

	idt = didt
	itunit = nint(dgetpar('itunit'))
	idtorig = idt

	end

c********************************************************************

	subroutine setup_date

c setup and check date parameter

	implicit none

	integer date,time
	double precision ddate,dtime
	integer year,month,day,hour,min,sec
	integer ierr
	character*20 text
	double precision dgetpar

	call getfnm('date',text)
	if( text .ne. ' ' ) then	!string value given
	  call dtsunform(year,month,day,hour,min,sec,text,ierr)
	  if( ierr .ne. 0 ) then
	    write(6,*) 'date: ',text
	    stop 'error stop setup_date: cannot parse date'
	  end if
	  date = 10000*year + 100*month + day
	  time = 10000*hour + 100*min + sec
	  ddate = date
	  dtime = time
	  call dputpar('date',ddate)
	  call dputpar('time',dtime)
	  !write(6,*) '===================================='
	  !write(6,*) 'date as string: ',date,time
	  !write(6,*) '===================================='
	end if

	date = nint(dgetpar('date'))
	time = nint(dgetpar('time'))

	call dtsini(date,time)

	if( date >= 0 ) return

	write(6,*) 'The date parameter has not been set.'
	write(6,*) 'The date parameter indicates the absolute date'
	write(6,*) 'to which FEM time is relative to.'
	write(6,*) 'The format for the date parameter is YYYY[MMDD].'
	write(6,*) 'If you want to do without this faeture, please'
	write(6,*) 'explicitly set the date parameter to 0'
	write(6,*) 'in the parameter file.'
	write(6,*) 'Alternatively, you can also specify date as a string.'
	write(6,*) 'The format is ''YYYY-MM-DD[::hh:mm:ss]'''
	write(6,*) 'Examples:'
	write(6,*) '   date = 20120101'
	write(6,*) '   date = 2012                 # same as above'
	write(6,*) '   date = ''2012-01-01''         # same as above'
	write(6,*) '   date = 0                    # no date feature'

	stop 'error stop ckdate: date'
	end

c********************************************************************

	subroutine get_date_time(date,time)

	implicit none

	integer date,time

	double precision dgetpar

	date = nint(dgetpar('date'))
	time = nint(dgetpar('time'))

	end

c********************************************************************
c********************************************************************
c********************************************************************

        subroutine check_timestep(irepeat)

c still to be implemented

	implicit none

	integer irepeat		!on return 1 if time step has to be repeated

	irepeat = 0

	end

c********************************************************************
c********************************************************************
c********************************************************************

        subroutine set_timestep

c controls time step

        implicit none

	include 'femtime.h'

	logical bdebug
        integer idtdone,idtrest,idts
	integer idtfrac,itnext
        integer istot
        double precision dt
	real dtr
        real ri,rindex,rindex1
	real perc,rmax

        real, save :: cmax,tfact,dtmin
        integer, save :: idtsync,isplit,idtmin
        integer, save :: iuinfo
        integer, save :: icall = 0

	double precision dgetpar

	bdebug = .false.

        if( icall .eq. 0 ) then
          istot = 0

          isplit = nint(dgetpar('itsplt'))
          cmax   = dgetpar('coumax')
          tfact  = dgetpar('tfact')	!still to be commented

	  call convert_time('idtsyn',idtsync)
	  call convert_time('idtmin',idtmin)
	  dtmin = dgetpar('idtmin')

          call getinfo(iuinfo)  !unit number of info file

	  if( isplit .ge. 0 .and. itunit .ne. 1 ) then
	    write(6,*) 'isplit, itunit: ',isplit,itunit
	    stop 'error stop set_timestep: itunit /= 1 not allowed here'
	  end if

          icall = 1
        end if

c----------------------------------------------------------------------
c        idtsync = 0             !time step for syncronization
c        cmax = 1.0              !maximal Courant number permitted
c        isplit = -1             !mode for variable time step:
c                                ! -1:  time step fixed, 
c                                !      no computation of rindex
c                                !  0:  time step fixed
c                                !  1:  split time step
c                                !  2:  optimize time step (no multiple)
c                                !  3:  optimize time step (fractional)
c	 idtmin = 1		 !minimum time step allowed
c	 tfact = 0		 !factor of maximum decrease of time step
c----------------------------------------------------------------------

        if( isplit .ge. 0 ) then
          dtr = 1.
          call hydro_stability(dtr,rindex)
        else
          rindex = 0.
        end if

	istot = 0

        if( isplit .le. 0 ) then
          idts = 0
	  dt = dt_act
        else if( isplit .eq. 1 ) then
	  if( idtorig <= 0 ) then
	    stop 'error stop set_timestep: idtorig==0 and isplit==1'
	  end if
          idts = idtorig
	  dt = dt_act
	  call split_equal(rindex,cmax,dt,istot)
        else if( isplit .eq. 2 .or. isplit .eq. 3 ) then
          idts = idtsync
	  idtfrac = 0
	  dt = dt_orig
	  if( rindex > 0 ) dt = cmax / rindex	! maximum allowed time step
	  if( dt >= dt_orig ) then
	    dt = dt_orig
	  else if( dt >= 1. ) then
	    if( isplit .eq. 2 ) dt = int(dt)
	  else
	    idtfrac = ceiling(1./dt)
	    if( isplit .eq. 2 ) dt = 1. / idtfrac
	  end if
        else
          write(6,*) 'isplit = ',isplit
          stop 'error stop set_timestep: value for isplit not allowed'
        end if

c----------------------------------------------------------------------
c dt	 is proposed new time step
c idts   is time step with which to syncronize
c istot  is number of internal time steps (only for isplit = 1)
c rindex is computed stability index (refers to time step == 1)
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c	syncronize time step
c----------------------------------------------------------------------

	itnext = 0
	bsync = .false.		!true if time step has been syncronized

	if( t_act + dt .gt. itend ) then	!sync with end of sim
	  dt = itend - t_act
	  bsync = .true.
        else if( idts .gt. 0 ) then               !syncronize time step
	  itnext = itanf + idts * ceiling( (t_act-itanf)/idts )	!is integer
	  if( itnext == it ) itnext = itnext + idts
	  if( t_act + dt > itnext ) then
	    dt = itnext - t_act
	    bsync = .true.
	  end if
        end if

c----------------------------------------------------------------------
c ri     is stability index for computed time step
c dtmin  is minimum time step allowed
c----------------------------------------------------------------------

        ri = dt*rindex

        if( dt .lt. dtmin .and. .not. bsync ) then    !should never happen
	  dtr = dt
          call error_stability(dtr,rindex)
          write(6,*) 'dt is less than dtmin'
          !write(6,*) it,itanf,mod(it-itanf,idtorig)
          !write(6,*) idtnew,idtdone,idtrest,idtorig
          write(6,*) idtorig,itnext
          write(6,*) idts,idtsync
          write(6,*) t_act,dt,dtmin
          write(6,*) isplit,istot
          write(6,*) cmax,rindex,ri
	  write(6,*) 'possible computed time step:  dt = ',dt
	  write(6,*) 'minimum time step allowed: dtmin = ',dtmin
	  write(6,*) 'please lower dtmin in parameter input file'
          stop 'error stop set_timestep: time step too small'
        end if

        niter=niter+1

	dt_act = dt
	t_act = t_act + dt
	idt = dt
	it = t_act

	if( bdebug ) then
	  write(107,*) '========================'
	  write(107,*) it,idt,t_act,dt_act,bsync
	  write(107,*) itanf,itend
	  write(107,*) rindex,ri
	  write(107,*) '========================'
	end if

	!perc = (100.*(it-itanf))/(itend-itanf)
	perc = (100.*(t_act-itanf))/(itend-itanf)

        write(iuinfo,1001) '----- new timestep: ',it,idt,perc
        write(iuinfo,1002) 'set_timestep: ',it,ri,rindex,istot,idt

        return
 1001   format(a,i12,i8,f8.2)
 1002   format(a,i12,2f12.4,2i8)
        end

c**********************************************************************

	subroutine split_equal(rindex,cmax,dt,istot)

c split time step (macro timestep) into equal parts
c
c dt is the new time step, all the rest can be computed from dt

	implicit none

	real rindex,cmax
	double precision dt
	integer istot

	include 'femtime.h'

	integer idtnew
	real ri

        if( mod(it-itanf,idtorig) .ne. 0 ) return	!inside macro time step

        istot = idtorig*rindex/cmax	!starting point for istot
        if( istot < 1 ) istot = 1

        do 
          idtnew = idtorig/istot
	  if( idtnew == 0 ) exit	!istot too big
          ri = idtnew*rindex
          if( ri .le. cmax .and. istot*idtnew == idtorig ) exit
          istot = istot + 1
        end do

	if( idtnew == 0 ) then
	  write(6,*) rindex,cmax,istot,idtnew
	  stop 'error stop split_equal: cannot find new time step'
	end if

	dt = idtnew

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine is_time_first(bfirst)

c true if in initialization phase

	implicit none

	logical bfirst

	include 'femtime.h'

	bfirst = it .eq. itanf

	end

c**********************************************************************

	subroutine is_time_last(blast)

c true if in last time step

	implicit none

	logical blast

	include 'femtime.h'

	blast = it .eq. itend

	end

c**********************************************************************

        subroutine get_act_time(itact)

c returns actual time

        implicit none

	integer itact

	include 'femtime.h'

	itact = it

	end

c**********************************************************************

        subroutine get_first_time(itfirst)

c returns first (initial) time

        implicit none

	integer itfirst

	include 'femtime.h'

	itfirst = itanf

	end

c**********************************************************************

        subroutine get_last_time(itlast)

c returns end time

        implicit none

	integer itlast

	include 'femtime.h'

	itlast = itend

	end

c**********************************************************************

        subroutine get_timestep(dt)

c returns real time step (in real seconds)

        implicit none

	real dt		!time step (return)

	include 'femtime.h'

	!dt = idt/float(itunit)
	dt = dt_act

	end

c**********************************************************************

        subroutine get_orig_timestep(dt)

c returns original real time step (in real seconds)

        implicit none

	real dt		!time step (return)

	include 'femtime.h'

	dt = idtorig/float(itunit)

	end

c********************************************************************
c********************************************************************
c********************************************************************

