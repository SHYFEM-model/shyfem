!
! $Id: subqfxf.f,v 1.10 2009-09-14 08:20:58 georg Exp $
!
! heat flux module (file administration)
!
! contents :
!
! subroutine qfinit(file)                       initializes heat flux module
! subroutine qfget(qs,ta,tb,uw,cc,ur,p,e,r,q)   gets actual values for air
! subroutine qfunit(iu)                         gets unit of heat flux file
! subroutine qfmake(it)                         reads from heat flux file
! subroutine qfnext(it,qs,ta,tb,uw,cc,ur,p,e,r,q,ier) reads next record
! subroutine qfcheck(it,qs,ta,tb,ur,uw,cc,p)	checks for unrealistic values
! subroutine qfcheck_file(file)			checks whole file
! subroutine qfperiodic(itp)			sets itperiod
! subroutine qfrhumid(irh)			sets irhumid
! subroutine qftest                             test drives qfxf routines
!
! revision log :
!
! 01.06.1998	ggu&lz	written from scratch (nearly)
! 24.06.1998	ggu&lz	subroutines from lucia integrated
! 24.05.2000	ggu	no details in element treatment, 3D algorithm
! 01.02.2002	ggu	only file handling routines here
! 01.02.2006	ggu	distinguish between humidity and wet bulb temp - bhumid
! 20.05.2007	ggu	itperiod may be 0 -> no periodic BC
! 12.11.2008	ggu	checks for unrealistic values in qfcheck
! 10.03.2009	ggu	call to meteo_set_matrix() for 2D arrays
! 27.08.2009	ggu	set itperiod and irhumid from outside
! 17.05.2011	ggu	compiler warnings in qfcheck()
!
! notes :
!
! qs	net solar radiation (reflection already subtracted)
! ta	air temperature
! tb	wet bulb temperature
! uw	wind speed
! cc	cloud cover (0 clear sky, 1 totally covered)
!
! to initialize call qfinit
! at every time step call qfmake -> reads if necessary new record
! to get data call qfget
! to see if a file is opened call qfunit
!
!***********************************************************
!-----------------------------------------------------------
        module heat_admin
!-----------------------------------------------------------
        contains
!-----------------------------------------------------------

	subroutine qfinit(file)

! initializes heat flux module

        use fil

	implicit none

	character*(*) file

        include 'subqfx.h'

	integer ier

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

	call qfnext(itfold,qsold,taold,tbold,uwold,ccold,urold,pold,eold,rold,qold,ier)
	if( ier .ne. 0 ) goto 99

	call qfnext(itfnew,qsnew,tanew,tbnew,uwnew,ccnew,urnew,pnew,enew,rnew,qold,ier)
	if( ier .ne. 0 ) goto 99

	return
   99	continue
	write(6,*) ier
	stop 'error stop qfinit: error reading file'
	end

!***********************************************************

	subroutine qfinit_internal

! internal initialization routine

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

!***********************************************************

	subroutine qfget(qs,ta,tb,uw,cc,ur,p,e,r,q)

! gets actual values for air

! qs	solar radiation [W/m**2]
! ta	air temperature [C]
! tb	wet bulb temperature [C]
! ur	relative humidity [%], ([0-100])
! uw	wind speed [m/s]
! cc	cloud cover [0-1], 0 no clouds
! p     pressure [mb]
! e     vapor pressure [mb]
! r     mixing ratio [0-1]
! q     specific humidity [0-1]

        implicit none

	double precision qs,ta,tb,uw,cc,ur,p,e,r,q

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

!***********************************************************

        subroutine qfunit(iu)

! gets unit of heat flux file

        implicit none

        integer iu

        include 'subqfx.h'

        iu = ifunit

        end

!***********************************************************

	subroutine qfmake(it)

! reads from heat flux file and makes vars

	implicit none

	integer it

        include 'subqfx.h'

	character*80 file
	logical bdebug,bwrite
	integer itact,ier
	double precision rt,rtt

	integer itmany
	save    itmany
	data    itmany / 0 /

! itperiod - periodic radiation conditions - see subqfx.h
!
! the file must contain values for time 0 and itperiod
! these records should be equal
!
! year : 31536000      day : 86400
!
! set itperiod to 0 if no periodic conditions are needed

	if( ifunit .le. 0 ) return

	bdebug = .true.
	bdebug = .false.
	bwrite = .true.

! handle periodic conditions

	itact = it - itmany * itperiod

	if( itperiod .gt. 0 .and. itact .gt. itperiod ) then
	  rewind(ifunit)
	  itmany = itmany + 1
	  itact  = itact  - itperiod
	  itfold = itfold - itperiod
	  itfnew = itfnew - itperiod
	  write(6,*) 'period in qflux: ',it,itact,itfold,itfnew
	end if

! iterate until new record has time greater than actual time

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

	  call qfnext(itfnew,qsnew,tanew,tbnew,uwnew,ccnew,urnew,pnew,enew,rnew,qnew,ier)

	  if( ier .ne. 0 ) goto 99

	end do

! do interpolation

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
	  write(6,'(a,i10,5f12.2)') 'qfmake: ',itact,qsact,taact,tbact,uwact,ccact
	end if

	return
   99	continue
	write(6,*) it,itact,ier
	stop 'error stop qfmake: error reading file'
	end

!***********************************************************

	subroutine qfnext(it,qs,ta,tb,uw,cc,ur,p,e,r,q,ier)

! reads next record of flux file
!
! it	time [s]
! qs	solar radiation [W/m**2]
! ta	air temperature [C]
! tb	wet bulb temperature [C]
! ur	relative humidity [%], ([0-100])
! uw	wind speed [m/s]
! cc	cloud cover [0-1], 0 no clouds
! p     pressure [mb]
! e     vapor pressure [mb]
! r     mixing ratio [0-1]
! q     specific humidity [0-1]
! ier   error code
!
! format is one of the two following:
!
! it qs ta tb uw cc
! it qs ta ur uw cc

        use heat_util
        use meteo_forcing,      only: meteo_set_matrix
        use heat_util2

	implicit none

	integer it
	double precision qs,ta,tb,uw,cc
	double precision ur,p,e,r,q
	integer ier

        include 'subqfx.h'

        double precision pstd
        parameter ( pstd = 1013.25 )

	logical bdebug
	logical bhumid
	integer ios
	integer iunit
	double precision raux

	integer itold
	save itold
	data itold / 0 /

	bhumid = .true.		!relative humidity is given in file - else wbt
	bhumid = irhumid .gt. 0

	bdebug = .true.
	bdebug = .false.

	ier = 0

!-------------------------------------------------------
! read next record
!-------------------------------------------------------

        call qfunit(iunit)
	if( iunit .le. 0 ) return

	read(iunit,*,iostat=ios) it,qs,ta,raux,uw,cc

!-------------------------------------------------------
! error handling and debug
!-------------------------------------------------------

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

!-------------------------------------------------------
! convert relative humidity or wet bulb temperature to each other
!-------------------------------------------------------

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

!-------------------------------------------------------
! check data
!-------------------------------------------------------

	if( bdebug ) then
	  write(6,1000) 'qfnext: ',it,qs,ta,tb,uw,cc,ur,p,e,r,q
	end if

	call qfcheck(it,qs,ta,tb,ur,uw,cc,p)

	itold = it

!-------------------------------------------------------
! end of routine
!-------------------------------------------------------

        return
 1000   format(a,i10,f7.1,2f6.1,f5.1,f5.2,f6.1,f7.1,f5.1,2f6.3)
	end

!*****************************************************************************

	subroutine qfcheck(it,qs,ta,tb,ur,uw,cc,p)

! checks heat flux values for unrealistic values

	implicit none

	integer it
	double precision qs,ta,tb,ur,uw,cc,p

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

!*****************************************************************************

	subroutine qfcheck_file(file)

! checks whole heat flux file for unusual values

	implicit none

	character*(*) file

	integer ier,it
	double precision qs,ta,tb,uw,cc,ur,p,e,r,q

	write(6,*) 'checking heat flux file...'

        ier = 0
        call qfinit(file)
        do while( ier .eq. 0 )
          call qfnext(it,qs,ta,tb,uw,cc,ur,p,e,r,q,ier)
        end do
        call qfinit(" ")        !close file

	write(6,*) 'heat flux file is ok.'

	end

!*****************************************************************************

	subroutine qfperiodic(itp)

! sets itperiod - can be called at any time

	implicit none

	integer itp

        include 'subqfx.h'

	call qfinit_internal
	itperiod = itp

	end

!*****************************************************************************

	subroutine qfrhumid(irh)

! sets irhumid - can be called at any time
!
! irh==1  =>  relative humidity is given (default)

	implicit none

	integer irh

        include 'subqfx.h'

	call qfinit_internal
	irhumid = irh

	end
	
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************

	subroutine qftest

! test drives qfxf routines

	implicit none

	double precision qs,ta,tb,uw,cc,ur
        double precision p,e,r,q
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

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************

!        call qftest
!        end

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************

!-----------------------------------------------------------
        end module heat_admin
!-----------------------------------------------------------
