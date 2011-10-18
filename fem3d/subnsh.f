c
c $Id: subnsh.f,v 1.54 2010-03-22 15:29:31 georg Exp $
c
c utility routines for hp
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
c
c************************************************************
c
	subroutine prilog
c
c writes output to terminal or log file

	implicit none

	include 'modules.h'

	character*80 name

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real grav,fcor,dcor,dirn,rowass,roluft
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it
	character*80 descrp,descrr
	common /descrp/ descrp
	common /descrr/ descrr

	call getfnm('runnam',name)
	write(6,*)
	write(6,*) '     Name of run :'
	write(6,*) name

	write(6,*)
	write(6,*) '     Description of run :'
	write(6,*) descrp
	write(6,*)

	write(6,*) '     itanf,itend,idt :',itanf,itend,idt
	write(6,*) '     Iterations to go :',nits

	call getfnm('basnam',name)
	write(6,*)
	write(6,*) '     Name of basin :'
	write(6,*) name

	write(6,*)
	write(6,*) '     Description of basin :'
	write(6,*) descrr
	write(6,*)

	write(6,*) '     nkn,nel,ngr,mbw :',nkn,nel,ngr,mbw
	write(6,*) '     nbc,nrz,nrq,nrb     :',nbc,nrz,nrq,nrb
	write(6,*) '     dcor,dirn       :',dcor,dirn

	write(6,*)
	write(6,*) '     Values from parameter file :'
	write(6,*)

	call pripar(6)
	call check_parameter_values('prilog')

	call prbnds		!prints boundary info

	call prwnds		!prints wind info

	call prexta		!prints extra points

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

	subroutine dobefor

c to do before time step

	implicit none

	include 'modules.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	call modules(M_BEFOR)

	call adjust_chezy

	end

c********************************************************************

	subroutine doafter

c to do after time step

	implicit none

	include 'modules.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	call modules(M_AFTER)

c	call wrouta
	call wrousa
c	call wrexta(it)
	call wrflxa(it)
	call wrvola(it)

        call resid
        call rmsvel

        call admrst             !restart

c        call tsmed
	call ts_shell

	call wrnetcdf		!output in netcdf format

	call custom(it)

	end

c********************************************************************

	subroutine pritst(id)
c
c test output of all constants and variables
c
c id    identifier
c
	implicit none

	integer id

	include 'modules.h'

	character*80 descrp,descrr
	common /descrp/ descrp
	common /descrr/ descrr

	write(6,*)
	write(6,*) '============================================='
	write(6,*) '================ test output ================'
	write(6,*) '============================================='
	write(6,1) '================ id = ',id,' ================'
	write(6,*) '============================================='
	write(6,*)
c	write(6,*) '/nkonst/'
c	write(6,*) nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
c	write(6,*) '/mkonst/'
c        write(6,*) eps1,eps2,pi,flag,high,higi
c	write(6,*) '/pkonst'
c	write(6,*) grav,fcor,dcor,dirn,rowass,roluft
c	write(6,*) '/femtim/'
c	write(6,*) itanf,itend,idt,nits,niter,it
c
c	call tslgr

	write(6,*) '/descrp/'
	write(6,*) descrp
	write(6,*) '/descrr/'
	write(6,*) descrr

	call tsexta

	call tsflxa

	call tsvola

	call tsarea

	call tsbnds

	call tswnds

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

	implicit none

	integer iunit

	include 'modules.h'

c---------------------------------------------------------------
	integer nlvdi,nlv	!for 3d model
	common /level/ nlvdi,nlv
	real hlv(1)
	common /hlv/hlv
c---------------------------------------------------------------

	character*6 section
	integer nsc,num
c	integer nrdsec,nrdveci,nrdvecr
	integer nrdsec,nrdvecr
	character*80 vers

	logical hasreadsec

        nlv = 0         !is initialized really only in adjust_levels

	nsc = 0
	if(iunit.le.0) goto 63

	call nrdini(iunit)

c read loop over sections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do while( nrdsec(section,num) .ne. 0 )

c		write(6,*) 'new section: ',section,num

		call setsec(section,num)		!remember section

		nsc = nsc + 1

		if(section.eq.'title') then
			call rdtitl
			if( nsc .ne. 1 ) goto 66
		else if(section.eq.'end') then
			goto 69
		else if(section.eq.'para') then
			call nrdins(section)
		else if(section.eq.'extra') then
			call rdexta
		else if(section.eq.'area') then
			call rdarea
		else if(section.eq.'name') then
			call nrdins(section)
		else if(section.eq.'bound') then
			call rdbnds(num)
		else if(section.eq.'float') then
c			call rdfloa(nfldin)
			stop 'error stop: float section not supported'
		else if(section.eq.'close') then
			call rdclos(num)
		else if(section.eq.'flux') then
			call rdflxa
		else if(section.eq.'vol') then
			call rdvola
		else if(section.eq.'wind') then
			call rdwnds
		else if(section.eq.'oxypar') then	!oxygen
			call nrdins(section)
		else if(section.eq.'oxyarr') then	!oxygen
c			call rdoxy
			call nrdskp
		else if(section.eq.'bfmsc')then        ! BFM ECO TOOL
                       call nrdins(section)
		else if(section.eq.'levels') then
			nlv = nrdvecr(hlv,nlvdi)
			if( nlv .lt. 0 ) goto 77
                else if(section.eq.'lagrg')then
                        call nrdins(section)
                else if(section.eq.'sedtr')then         !sediment
                        call readsed
		else					!try modules
			call modules(M_READ)
			if( .not. hasreadsec() ) then	!sec has been handled?
				goto 97			! -> no
			end if
		end if

	end do		!loop over sections

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
	stop 'error stop : nlsh2d'
   77	continue
	if( nlv .eq. -1 ) then
	  write(6,*) 'read error in section $levels'
	  stop 'error stop nlsh2d: read error'
	else
	  write(6,*) 'dimension error in section $levels'
	  write(6,*) 'nlvdim = ',nlvdi,'   number of data read = ',-nlv
	  stop 'error stop nlsh2d: dimension error'
	end if
   97	continue
	write(6,*) 'Cannot handle section : ',section
	stop 'error stop nlsh2d: no such section'
	end

c************************************************************************

	subroutine rdtitl

c reads title section

	implicit none

	character*80 descrp
	common /descrp/ descrp

	character*80 line

	integer nrdlin,nrdsec
	integer num

	if( nrdlin(line) .eq. 0 ) goto 65
	call triml(line)
	descrp=line
	if( nrdlin(line) .eq. 0 ) goto 65
	call triml(line)
	call putfnm('runnam',line)
	if( nrdlin(line) .eq. 0 ) goto 65
	call triml(line)
	call putfnm('basnam',line)

	!if( nrdlin(line) .gt. 0 ) goto 65
	if( nrdsec(line,num) .eq. 0 ) goto 65
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

	subroutine impini

c initializes parameters for semi-implicit time-step

	implicit none

	integer itimpl
	common /semimi/ itimpl
	real weight
	common /semimr/ weight

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer icall
	save icall
	data icall / 0 /

	if( icall .gt. 0 ) return

	icall = 1
	weight = 0.5
	itimpl = itanf

	end

c**********************************************************************

	function bimpli(it)

c checks if semi-implicit time-step is active

	implicit none

	logical bimpli
	integer it

	integer itimpl
	common /semimi/ itimpl

	call impini

	if( it .le. itimpl ) then
	   bimpli = .true.
	else
	   bimpli = .false.
	end if

	end

c**********************************************************************

	subroutine getaz(azpar)

c returns actual az

	implicit none

	real azpar

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real ampar
	real getpar

	azpar=getpar('azpar')
	ampar=0.			!dummy

	call changeimp(it,azpar,ampar)

	end

c**********************************************************************

	subroutine getazam(azpar,ampar)

c returns actual az,am

	implicit none

	real azpar
	real ampar

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real getpar

	azpar=getpar('azpar')
	ampar=getpar('ampar')

	call changeimp(it,azpar,ampar)

	end

c**********************************************************************

	subroutine changeimp(it,azpar,ampar)

c changes parameters for semi-implicit time-step if necessary

	implicit none

	integer it
	real azpar,ampar

	integer itimpl
	common /semimi/ itimpl
        real weight
        common /semimr/ weight

	call impini

	if( it .le. itimpl ) then
	  azpar = weight
	  ampar = weight
	end if

	end

c**********************************************************************

	subroutine setimp(it,aweigh)

c sets parameters for semi-implicit time-step

	implicit none

	integer it
	real aweigh

	integer itimpl
	common /semimi/ itimpl
        real weight
        common /semimr/ weight

	call impini

	itimpl = it
	weight = aweigh

	write(6,*) 'implicit parameters changed: ',itimpl,weight

	end

c**********************************************************************

	function getimp()

c gets weight for semi-implicit time-step

	implicit none

	real getimp

        real weight
        common /semimr/ weight

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

        integer iu
        save iu
        data iu / 0 /

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

	subroutine pritime

c prints time after time step

	implicit none

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        integer nit1,nit2,naver
        real perc

	integer time,date
	integer year,month,day,hour,min,sec
	double precision dgetpar

	save date

	integer icall
	save icall
	data icall / 0 /

c---------------------------------------------------------------
c initialize date
c---------------------------------------------------------------

	if( icall .eq. 0 ) then
	  date = nint(dgetpar('date'))
	  time = nint(dgetpar('time'))
	  call dtsini(date,time)
	end if

	if( date .ne. 0 ) then
	  call dts2dt(it,year,month,day,hour,min,sec)
	end if

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
        nit2 = nint( niter * ( 1 + float(itend-it)/(it-itanf) ) )

        nits = nit2
        if( naver .gt. 0 ) nits = ( 1*nit1 + (naver-1)*nit2 ) / naver

c---------------------------------------------------------------
c write to terminal
c---------------------------------------------------------------

	if( date .eq. 0 ) then
          write(6,1001) it,idt,niter,nits,perc
	else
	  if( mod(icall,50) .eq. 0 ) then
            write(6,1003)
	  end if
          write(6,1002) it,year,month,day,hour,min,sec
     +			,idt,niter,nits,perc
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
	end

c********************************************************************

	subroutine endtime

c prints stats after last time step

	implicit none

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	write(6,1035) it,niter
 1035   format(' program stop at time =',i10,' seconds'/
     +         ' iterations = ',i10)

	end

c********************************************************************

        subroutine set_timestep

c controls time step

        implicit none

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        integer idtdone,idtrest,idts
        real dt
        real ri,rindex,rindex1,riold
	real perc
        real    cmax,rmax
        integer iloop,itloop
        integer idtsync,isplit,idtmin
	integer itunit
	logical bsync

        real getpar

        integer idtorig,istot,idtold
        save    idtorig,istot,idtold
        integer iuinfo
        save    iuinfo

        integer icall
        save icall
        data icall / 0 /

        if( icall .eq. 0 ) then
          idtorig = idt
          istot = 0
          call set_orig_timestep(idtorig)

          isplit  = getpar('itsplt')
          cmax    = getpar('coumax')
          idtsync = getpar('idtsyn')

          call getinfo(iuinfo)  !unit number of info file

          call get_timeunit(itunit)
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
c----------------------------------------------------------------------

        isplit  = nint(getpar('itsplt'))
        cmax    = getpar('coumax')
        idtsync = nint(getpar('idtsyn'))
        idtmin  = nint(getpar('idtmin'))

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
	iloop = 0
	istot = 0

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
              idt = idtorig/istot
              if( istot*idt .ne. idtorig ) idt = idt + 1
              ri = idt*rindex/idtorig
              iloop = iloop + 1
            end do
          end if
        else if( isplit .eq. 2 ) then
          idts = idtsync
          idt = idtorig
	  if( rindex / cmax .gt. 1. ) then	!BUGFIX
	    idt = min(idt,int(dt*cmax/rindex))
	  end if
	  !write(6,*) '++++++++ ',idts,idtorig,rindex,idt
        else
          write(6,*) 'isplit = ',isplit
          stop 'error stop set_timestep: value for isplit not allowed'
        end if

	!write(6,*) 'stability: ',it,rindex
	!write(6,*) idtorig,idts,idt,istot,ri,cmax

c	syncronize at least with end of simulation

	if( it + idt .gt. itend ) idt = itend - it

c----------------------------------------------------------------------
c idt    is proposed new time step
c idts   is time step with which to syncronize
c nits   is total number of time steps
c rindex is computed stability index
c ri     is used stability index
c istot  is number of internal time steps (only for isplit = 1)
c----------------------------------------------------------------------

	bsync = .false.		!true if time step has been syncronized
        if( idts .gt. 0 ) then               !syncronize time step
            idtdone = mod(it-itanf,idts)     !already done
            idtrest = idts - idtdone         !still to do
            if( idt .gt. idtrest ) then
	      idt = idtrest
	      bsync = .true.
	    end if
        end if

        ri = idt*rindex/idtorig

	if( itloop .gt. 10 ) then
	  stop 'error stop set_timestep: too many loops'
	end if

        if( idt .lt. 1 ) then           !should never happen
          call error_stability(dt,rindex)
          write(6,*) 'idt is less equal 0 !!!'
          write(6,*) it,itanf,mod(it-itanf,idtorig)
          write(6,*) idt,idtdone,idtrest,idtorig
          write(6,*) idts,idtsync,iloop
          write(6,*) isplit,istot
          write(6,*) cmax,rindex,ri
          stop 'error stop set_timestep: non positive time step'
        end if

        if( .not. bsync .and. idt .lt. idtmin ) then
	  rmax = (1./idtmin)*cmax
	  call eliminate_stability(rmax)
	  write(6,*) 'repeating time step for low idt: ',idt,rmax
	  goto 1
	end if

	riold = idtold*ri/idt	!ri with old time step
	idtold = idt

        niter=niter+1
        it=it+idt

	perc = (100.*(it-itanf))/(itend-itanf)

        write(iuinfo,1001) '----- new timestep: ',it,idt,perc
        write(iuinfo,1002) 'set_timestep: ',it,ri,riold,rindex,istot,idt

        return
 1001   format(a,i12,i8,f8.2)
 1002   format(a,i12,3f12.4,2i8)
        end

c**********************************************************************

        subroutine get_time(itact)

c returns actual time

        implicit none

	integer itact

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	itact = it

	end

c**********************************************************************

        subroutine get_timestep(dt)

c returns real time step (in real seconds)

        implicit none

	real dt		!time step (return)

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        integer itunit,idtorig
        common /femtimu/ itunit,idtorig
	save /femtimu/

	dt = idt/float(itunit)

	end

c**********************************************************************

        subroutine set_orig_timestep(idt)

c returns original real time step (in real seconds)

        implicit none

	integer idt		!time step

        integer itunit,idtorig
        common /femtimu/ itunit,idtorig
	save /femtimu/

	idtorig = idt

	end

c**********************************************************************

        subroutine get_orig_timestep(dt)

c returns original real time step (in real seconds)

        implicit none

	real dt		!time step (return)

        integer itunit,idtorig
        common /femtimu/ itunit,idtorig
	save /femtimu/

	dt = idtorig/float(itunit)

	end

c********************************************************************

        subroutine get_timeunit(itu)

c returns time unit

        implicit none

	integer itu

        integer itunit,idtorig
        common /femtimu/ itunit,idtorig
	save /femtimu/

	itu = itunit

	end

c**********************************************************************

        subroutine set_timeunit(itu)

c sets time unit

        implicit none

	integer itu

        integer itunit,idtorig
        common /femtimu/ itunit,idtorig
	save /femtimu/

	itunit = itu

	end

c********************************************************************

	subroutine set_output_frequency(itm_out,idt_out,ia_out)

c sets-up array for output frequency

	implicit none

	integer itm_out		!minimum time for output
	integer idt_out		!time step for output
	integer ia_out(4)	!array where info is stored

	integer itanf,itend,idt
	integer itmout,idtout,itout
	real getpar

	itanf = nint(getpar('itanf'))
	itend = nint(getpar('itend'))
	idt = nint(getpar('idt'))

	itmout = itm_out
	idtout = idt_out

	if( itmout .eq. -1 ) itmout = itanf
	if( itmout .lt. itanf ) itmout = itanf
	if( idtout .lt. idt .and. idtout .gt. 0 ) idtout = idt

	itout = itmout
	if( itmout .eq. itanf ) itout = itout + idtout

	if( itout .gt. itend ) idtout = 0

	ia_out(1) = idtout	! time step of output
	ia_out(2) = itmout	! first output
	ia_out(3) = itout	! next output
	ia_out(4) = 0		! unit (optional)

	end

c********************************************************************

	function next_output(ia_out)

c checks if time has come for output

	implicit none

	logical next_output
	integer ia_out(4)

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer idtout,itnext

	idtout = ia_out(1)
	itnext = ia_out(3)
	next_output = .false.

	if( idtout .le. 0 ) return
	if( itnext .gt. it ) return

	do while( itnext .le. it )
	  itnext = itnext + idtout
	end do

	ia_out(3) = itnext
	next_output = .true.

	end

c********************************************************************
c********************************************************************
c********************************************************************

        subroutine total_energy

c writes info on total energy to info file

	implicit none

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real kenergy,penergy,tenergy

	integer iuinfo
	save iuinfo
	data iuinfo / 0 /

	if( iuinfo .eq. 0 ) then
          call getinfo(iuinfo)  !unit number of info file
	end if

	call energ3d(kenergy,penergy)
	tenergy = kenergy + penergy

	write(iuinfo,*) 'energy: ',it,kenergy,penergy,tenergy

	end

c********************************************************************

