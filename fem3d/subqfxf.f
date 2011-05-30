c
c $Id: subqfxf.f,v 1.10 2009-09-14 08:20:58 georg Exp $
c
c heat flux module (file administration)
c
c contents :
c
c subroutine qfinit(file)                       initializes heat flux module
c subroutine qfget(qs,ta,tb,uw,cc,ur,p,e,r,q)   gets actual values for air
c subroutine qfunit(iu)                         gets unit of heat flux file
c subroutine qfmake(it)                         reads from heat flux file
c subroutine qfnext(it,qs,ta,tb,uw,cc,ur,p,e,r,q,ier) reads next record
c subroutine qfcheck(it,qs,ta,tb,ur,uw,cc,p)	checks for unrealistic values
c subroutine qfcheck_file(file)			checks whole file
c subroutine qfperiodic(itp)			sets itperiod
c subroutine qfrhumid(irh)			sets irhumid
c subroutine qftest                             test drives qfxf routines
c
c revision log :
c
c 01.06.1998	ggu&lz	written from scratch (nearly)
c 24.06.1998	ggu&lz	subroutines from lucia integrated
c 24.05.2000	ggu	no details in element treatment, 3D algorithm
c 01.02.2002	ggu	only file handling routines here
c 01.02.2006	ggu	distinguish between humidity and wet bulb temp - bhumid
c 20.05.2007	ggu	itperiod may be 0 -> no periodic BC
c 12.11.2008	ggu	checks for unrealistic values in qfcheck
c 10.03.2009	ggu	call to meteo_set_matrix() for 2D arrays
c 27.08.2009	ggu	set itperiod and irhumid from outside
c 17.05.2011	ggu	compiler warnings in qfcheck()
c
c notes :
c
c qs	net solar radiation (reflection already subtracted)
c ta	air temperature
c tb	wet bulb temperature
c uw	wind speed
c cc	cloud cover (0 clear sky, 1 totally covered)
c
c to initialize call qfinit
c at every time step call qfmake -> reads if necessary new record
c to get data call qfget
c to see if a file is opened call qfunit
c
c***********************************************************

	subroutine qfinit(file)

c initializes heat flux module

	implicit none

	character*(*) file

        include 'subqfx.h'

	integer ier
	integer ifileo

	call qfinit_internal

	if( file .eq. " " ) then
	  if( ifunit .gt. 0 ) close(ifunit)
	  ifunit = -1
	  return
	else if( ifunit .gt. 0 ) then
	  write(6,*) 'heat flux file already open ... close first'
	  stop 'error stop qfinit'
	end if

	ifunit = ifileo(0,file,'form','old')
	if( ifunit .le. 0 ) then
	    write(6,'(a,a)') 'filename: ',file(1:65)
	    stop 'error stop qfinit: Cannot open file'
	end if

	write(6,*) 'heat flux file opened (version 2.0): '
	write(6,*) file

	call qfnext(itfold,qsold,taold,tbold,uwold,ccold
     +                  ,urold,pold,eold,rold,qold                              
     +                  ,ier)
	if( ier .ne. 0 ) goto 99

	call qfnext(itfnew,qsnew,tanew,tbnew,uwnew,ccnew
     +                  ,urnew,pnew,enew,rnew,qold                              
     +                  ,ier)
	if( ier .ne. 0 ) goto 99

	return
   99	continue
	write(6,*) ier
	stop 'error stop qfinit: error reading file'
	end

c***********************************************************

	subroutine qfinit_internal

c internal initialization routine

	implicit none

        include 'subqfx.h'

	integer icall
	save    icall
	data    icall / 0 /

	if( icall .eq. 0 ) then
	  ifunit = 0
	  itperiod = 0
	  irhumid = 1
	  icall = 1
	end if

	end

c***********************************************************

	subroutine qfget(qs,ta,tb,uw,cc,ur,p,e,r,q)

c gets actual values for air

c qs	solar radiation [W/m**2]
c ta	air temperature [C]
c tb	wet bulb temperature [C]
c ur	relative humidity [%], ([0-100])
c uw	wind speed [m/s]
c cc	cloud cover [0-1], 0 no clouds
c p     pressure [mb]
c e     vapor pressure [mb]
c r     mixing ratio [0-1]
c q     specific humidity [0-1]

        implicit none

	real qs,ta,tb,uw,cc,ur,p,e,r,q

        include 'subqfx.h'

	qs = qsact
	ta = taact
	tb = tbact
	uw = uwact
	cc = ccact
	ur = uract
	p  = pact
	e  = eact
	r  = ract
	q  = qact

	end

c***********************************************************

        subroutine qfunit(iu)

c gets unit of heat flux file

        implicit none

        integer iu

        include 'subqfx.h'

        iu = ifunit

        end

c***********************************************************

	subroutine qfmake(it)

c reads from heat flux file and makes vars

	implicit none

	integer it

        include 'subqfx.h'

	character*80 file
	logical bdebug,bwrite
	integer itact,ier
	integer ifileo
	real rt,rtt

	integer itmany
	save    itmany
	data    itmany / 0 /

c itperiod - periodic radiation conditions - see subqfx.h
c
c the file must contain values for time 0 and itperiod
c these records should be equal
c
c year : 31536000      day : 86400
c
c set itperiod to 0 if no periodic conditions are needed

	if( ifunit .le. 0 ) return

	bdebug = .true.
	bdebug = .false.
	bwrite = .true.

c handle periodic conditions

	itact = it - itmany * itperiod

	if( itperiod .gt. 0 .and. itact .gt. itperiod ) then
	  rewind(ifunit)
	  itmany = itmany + 1
	  itact  = itact  - itperiod
	  itfold = itfold - itperiod
	  itfnew = itfnew - itperiod
	  write(6,*) 'period in qflux: ',it,itact,itfold,itfnew
	end if

c iterate until new record has time greater than actual time

	do while( itact .gt. itfnew )

	  itfold = itfnew
	  qsold = qsnew
	  taold = tanew
	  tbold = tbnew
	  uwold = uwnew
	  ccold = ccnew
	  urold = urnew
	  pold  = pnew
	  eold  = enew
	  rold  = rnew
	  qold  = qnew

	  call qfnext(itfnew,qsnew,tanew,tbnew,uwnew,ccnew
     +                  ,urnew,pnew,enew,rnew,qnew                              
     +                  ,ier)

	  if( ier .ne. 0 ) goto 99

	end do

c do interpolation

	rt = (itact-itfold) / float(itfnew-itfold)
        rtt = 1. - rt

        itfact = itact
	 qsact = rt * qsnew + rtt * qsold
	 taact = rt * tanew + rtt * taold
	 tbact = rt * tbnew + rtt * tbold
	 uwact = rt * uwnew + rtt * uwold
	 ccact = rt * ccnew + rtt * ccold
	 uract = rt * urnew + rtt * urold
	  pact = rt * pnew  + rtt * pold
	  eact = rt * enew  + rtt * eold
	  ract = rt * rnew  + rtt * rold
	  qact = rt * qnew  + rtt * qold

	if( bdebug ) then
	  write(6,'(a,i10,5f12.2)') 'qfmake: ',itact,qsact,taact
     +                  ,tbact,uwact,ccact
	end if

	return
   99	continue
	write(6,*) it,itact,ier
	stop 'error stop qfmake: error reading file'
	end

c***********************************************************

	subroutine qfnext(it,qs,ta,tb,uw,cc,ur,p,e,r,q,ier)

c reads next record of flux file
c
c it	time [s]
c qs	solar radiation [W/m**2]
c ta	air temperature [C]
c tb	wet bulb temperature [C]
c ur	relative humidity [%], ([0-100])
c uw	wind speed [m/s]
c cc	cloud cover [0-1], 0 no clouds
c p     pressure [mb]
c e     vapor pressure [mb]
c r     mixing ratio [0-1]
c q     specific humidity [0-1]
c ier   error code
c
c format is one of the two following:
c
c it qs ta tb uw cc
c it qs ta ur uw cc

	implicit none

	integer it
	real qs,ta,tb,uw,cc
	real ur,p,e,r,q
	integer ier

        include 'subqfx.h'

        real pstd
        parameter ( pstd = 1013.25 )

	logical bdebug
	logical bhumid
	integer ios
	integer iunit
	real raux

	integer itold
	save itold
	data itold / 0 /

	bhumid = .true.		!relative humidity is given in file - else wbt
	bhumid = irhumid .gt. 0

	bdebug = .true.
	bdebug = .false.

	ier = 0

c-------------------------------------------------------
c read next record
c-------------------------------------------------------

        call qfunit(iunit)
	if( iunit .le. 0 ) return

	read(iunit,*,iostat=ios) it,qs,ta,raux,uw,cc

c-------------------------------------------------------
c error handling and debug
c-------------------------------------------------------

	if( ios .ne. 0 ) then
	  if( ios .gt. 0 ) then
	    write(6,*) 'error reading heat flux file'
	    write(6,*) 'on line close to it = ',itold
	    stop 'error stop qfnext: read error'
	  else if( ios .lt. 0 ) then
	    write(6,*) 'end of file reading heat flux file'
	    write(6,*) 'on line close to it = ',itold
	    ier = ios
	    return
	  end if
	end if

c-------------------------------------------------------
c convert relative humidity or wet bulb temperature to each other
c-------------------------------------------------------

	if( bhumid ) then		!compute wet bulb temperature
	  ur = raux
	  call rh2twb(ta,ur,tb) 
	else				!compute relative humidity
	  tb = raux
	  call twb2rh(ta,tb,ur)
	end if

        p = pstd
        call vapor(ta,p,ur,e,r,q)

	call meteo_set_matrix(qs,ta,ur,tb,uw,cc)	!set 2D matrices

c-------------------------------------------------------
c check data
c-------------------------------------------------------

	if( bdebug ) then
	  write(6,1000) 'qfnext: ',it,qs,ta,tb,uw,cc,ur,p,e,r,q
	end if

	call qfcheck(it,qs,ta,tb,ur,uw,cc,p)

	itold = it

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

        return
 1000   format(a,i10,f7.1,2f6.1,f5.1,f5.2,f6.1,f7.1,f5.1,2f6.3)
	end

c*****************************************************************************

	subroutine qfcheck(it,qs,ta,tb,ur,uw,cc,p)

c checks heat flux values for unrealistic values

	implicit none

	integer it
	real qs,ta,tb,ur,uw,cc,p

	if( qs .lt. 0.    .or. qs .gt. 1500. ) goto 99
	if( ta .lt. -100. .or. ta .gt. 100.  ) goto 99
	if( tb .lt. -100. .or. tb .gt. 100.  ) goto 99
	if( ur .lt. 0.    .or. ur .gt. 100.  ) goto 99
	if( uw .lt. 0.    .or. uw .gt. 100.  ) goto 99
	if( cc .lt. 0.    .or. cc .gt. 1.    ) goto 99
	if(  p .lt. 900.  .or.  p .gt. 1100. ) goto 99

	return
   99	continue
	write(6,*) 'Unrealistic values in heat flux file for time: ',it
	write(6,*) 'qs,ta,tb,ur,uw,cc,p'
	write(6,*) qs,ta,tb,ur,uw,cc,p
	stop 'error stop qfcheck: unrealistic values'
	end

c*****************************************************************************

	subroutine qfcheck_file(file)

c checks whole heat flux file for unusual values

	implicit none

	character*(*) file

	integer ier,it
	real qs,ta,tb,uw,cc,ur,p,e,r,q

	write(6,*) 'checking heat flux file...'

        ier = 0
        call qfinit(file)
        do while( ier .eq. 0 )
          call qfnext(it,qs,ta,tb,uw,cc,ur,p,e,r,q,ier)
        end do
        call qfinit(" ")        !close file

	write(6,*) 'heat flux file is ok.'

	end

c*****************************************************************************

	subroutine qfperiodic(itp)

c sets itperiod - can be called at any time

	implicit none

	integer itp

        include 'subqfx.h'

	call qfinit_internal
	itperiod = itp

	end

c*****************************************************************************

	subroutine qfrhumid(irh)

c sets irhumid - can be called at any time
c
c irh==1  =>  relative humidity is given (default)

	implicit none

	integer irh

        include 'subqfx.h'

	call qfinit_internal
	irhumid = irh

	end
	
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************

	subroutine qftest

c test drives qfxf routines

	implicit none

	real qs,ta,tb,uw,cc,ur
        real p,e,r,q
	character*80 file
	integer it,iddt,ityear,nmax
	integer i

	file = '9899_80.qfx'
	iddt = 300
	file = 'qflux.in'
	iddt = 3600

	ityear = 3600 * 24 * 365
	nmax = ityear / iddt

	call qfinit(file)

	do i=1,nmax
	  it = i * iddt
	  call qfmake(it)
	  call qfget(qs,ta,tb,uw,cc,ur,p,e,r,q)
	  write(6,1000) 'qfnext: ',it,qs,ta,tb,uw,cc,ur,p,e,r,q
 1000     format(a,i10,f7.1,2f6.1,f5.1,f5.2,f6.1,f7.1,f5.1,2f6.3)
	end do

	end

c*****************************************************************************
c*****************************************************************************
c*****************************************************************************

c        call qftest
c        end

c*****************************************************************************
c*****************************************************************************
c*****************************************************************************

