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
        real perc

	integer year,month,day,hour,min,sec

	character*20 line
	double precision dgetpar
	logical dts_has_date

	integer icall
	save icall
	data icall /0/

c---------------------------------------------------------------
c set parameters and compute percentage of simulation
c---------------------------------------------------------------

        naver = 20
        naver = 0

        perc = (100.*(it-itanf))/(itend-itanf)

c---------------------------------------------------------------
c compute total number of iterations
c---------------------------------------------------------------

        nit1 = niter + (itend-it)/idt
	nit2 = nit1
	if( it .gt. itanf ) then
          nit2 = nint( niter * ( 1 + float(itend-it)/(it-itanf) ) )
	end if

        nits = nit2
        if( naver .gt. 0 ) nits = ( 1*nit1 + (naver-1)*nit2 ) / naver

c---------------------------------------------------------------
c write to terminal
c---------------------------------------------------------------

	if( dts_has_date() ) then
	  if( mod(icall,50) .eq. 0 ) write(6,1003)
!	  call dts2dt(it,year,month,day,hour,min,sec)
	  call dtsgf(it,line)
          write(6,1005) it,line,idt,niter,nits,perc
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

	double precision dgetpar

	call convert_date('itanf',itanf)
	call convert_date('itend',itend)
	call convert_time('idt',idt)

	if( idt .le. 0 .or. itanf+idt .gt. itend ) then
	   write(6,*) 'Error in compulsory time parameters'
	   write(6,*) 'itanf,itend,idt :',itanf,itend,idt
	   stop 'error stop : cktime'
	end if

	niter = 0
	it = itanf
	nits = (itend-itanf) / idt

	itunit = nint(dgetpar('itunit'))
	idtorig = idt

	end

c********************************************************************

	subroutine setup_date

c setup and check date parameter

	implicit none

	integer date,time
	double precision dgetpar

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
	write(6,*) 'Examples:'
	write(6,*) '   date = 20120101'
	write(6,*) '   date = 2012       # same as above'
	write(6,*) '   date = 0          # no date feature'

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

        integer idtdone,idtrest,idts
        real dt
        real ri,rindex,rindex1,riold
	real perc
        real    cmax,rmax,tfact
	save cmax,tfact
        integer iloop,itloop
        integer idtsync,isplit,idtmin
	save idtsync,isplit,idtmin
	integer irepeat
	logical bsync

        real getpar

        integer istot,idtold,idtnew
        save    istot,idtold,idtnew
        integer iuinfo
        save    iuinfo

        integer icall
        save icall
        data icall / 0 /

        if( icall .eq. 0 ) then
          istot = 0

          isplit  = getpar('itsplt')
          cmax    = getpar('coumax')
          tfact  = getpar('tfact')	!still to be commented

	  call convert_time('idtsync',idtsync)
	  call convert_time('idtmin',idtmin)

          call getinfo(iuinfo)  !unit number of info file

	  if( isplit .ge. 0 .and. itunit .ne. 1 ) then
	    write(6,*) 'isplit, itunit: ',isplit,itunit
	    stop 'error stop set_timestep: itunit != 1 not allowed here'
	  end if

	  idtold = 0
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
c	 idtmin = 0		 !minimum time step allowed
c	 tfact = 0		 !factor of maximum decrease of time step
c----------------------------------------------------------------------

	itloop = 0

    1	continue
	itloop = itloop + 1

        if( isplit .ge. 0 ) then
          dt = idtorig
          call hydro_stability(dt,rindex)
        else
          rindex = 0.
        end if

	ri = 0.
	istot = 0
	idtnew = idt

        if( isplit .le. 0 ) then
          idts = 0
        else if( isplit .eq. 1 ) then
          idts = idtorig
          if( mod(it-itanf,idtorig) .eq. 0 ) then      !end of macro timestep
            istot = rindex/cmax
            ri = 1. + cmax
            iloop = 0
            do while( ri .gt. cmax )
              istot = istot + 1
              idtnew = idtorig/istot
              if( istot*idtnew .ne. idtorig ) idtnew = idtnew + 1
              ri = idtnew*rindex/idtorig
              iloop = iloop + 1
	      if( iloop .gt. 100 ) then
		stop 'error stop set_timestep: too many iterations'
	      end if
            end do
          end if
        else if( isplit .eq. 2 ) then
          idts = idtsync
          idtnew = idtorig
	  if( rindex / cmax .gt. 1. ) then	!BUGFIX
	    idtnew = min(idtnew,int(dt*cmax/rindex))
	  end if
        else
          write(6,*) 'isplit = ',isplit
          stop 'error stop set_timestep: value for isplit not allowed'
        end if

c----------------------------------------------------------------------
c idtnew is proposed new time step
c idts   is time step with which to syncronize
c nits   is total number of time steps
c rindex is computed stability index
c ri     is used stability index
c istot  is number of internal time steps (only for isplit = 1)
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c	syncronize time step
c----------------------------------------------------------------------

	bsync = .false.		!true if time step has been syncronized

	if( it + idtnew .gt. itend ) then	!sync with end of sim
	  idtnew = itend - it
	  bsync = .true.
        else if( idts .gt. 0 ) then               !syncronize time step
            idtdone = mod(it-itanf,idts)     !already done
            idtrest = idts - idtdone         !still to do
            if( idtnew .gt. idtrest ) then
	      idtnew = idtrest
	      bsync = .true.
	    end if
        end if

        ri = idtnew*rindex/idtorig

	if( itloop .gt. 10 ) then
	  stop 'error stop set_timestep: too many loops'
	end if

        if( idtnew .lt. 1 ) then           !should never happen
          call error_stability(dt,rindex)
          write(6,*) 'idt is less equal 0 !!!'
          write(6,*) it,itanf,mod(it-itanf,idtorig)
          write(6,*) idtnew,idtdone,idtrest,idtorig
          write(6,*) idts,idtsync,iloop
          write(6,*) isplit,istot
          write(6,*) cmax,rindex,ri
          stop 'error stop set_timestep: non positive time step'
        end if

	irepeat = 0

        if( .not. bsync .and. idtnew .lt. idtmin ) then
	  rmax = (1./idtmin)*cmax
	  call eliminate_stability(rmax)
	  write(6,*) 'repeating time step for low idt (1): ',idtnew,rmax
	  irepeat = 1
	  goto 1	!eventually this should be commented
	end if

        if( .not. bsync .and. idtnew .lt. tfact*idtold ) then
	  write(6,*) 'repeating time step for low idt (2): ',idtnew,idtold
	  irepeat = 1
	end if

	if( irepeat .gt. 0 ) then
	  idtold = idtnew
	  write(6,*) 'We really should repeat this time step'
	  write(6,*) 'code not yet ready...'
	end if

	riold = idtold*ri/idtnew	!ri with old time step
	idtold = idtnew

        niter=niter+1
        it=it+idtnew
	idt=idtnew

	perc = (100.*(it-itanf))/(itend-itanf)

        write(iuinfo,1001) '----- new timestep: ',it,idt,perc
        write(iuinfo,1002) 'set_timestep: ',it,ri,riold,rindex,istot,idt

        return
 1001   format(a,i12,i8,f8.2)
 1002   format(a,i12,3f12.4,2i8)
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

	dt = idt/float(itunit)

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

