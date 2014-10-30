c
c $Id: suptim.f,v 1.8 2008-12-09 11:45:14 georg Exp $
c
c revision log :
c
c 12.02.1999  ggu     adapted to auto mode
c 27.05.2005  ggu     increase nrec always in oktime (even when it is the same)
c 13.11.2008  ggu     in oktime() increase irec only for new time
c 06.12.2008  ggu     in oktime() set itact to actual time
c 09.10.2010  ggu     in oktime() handle negative itfreq
c 05.09.2013  ggu     new routine endtime()
c 20.10.2014  ggu     completely restructured, old routines deleted
c
c******************************************************

c******************************************************
c******************************************************
c******************************************************
c
c itime		relative time in integer (old it)
c dtime		relative time in double precision (in fem files)
c atime		absolute time (seconds from 1.1.1)
c
c atime0	absolute time for fem time 0
c
c if no date/time is available then:
c
c		itime == dtime == atime
c		atime0 = 0
c
c******************************************************
c******************************************************
c******************************************************

        subroutine ptime_init

c initializes ptime

        implicit none

	integer ihigh
	parameter(ihigh=1000000000)
	double precision ahigh
	parameter(ahigh=4000.*365.*86400.)

	include 'timlim.h'

c the first 4 variables might be useless

        itmin = -ihigh
        itmax =  ihigh
	idto  = 0
	itact = itmin

c these are still used

        itfreq = 1
	nrec = 0

c these are the new used variables

	atime0 = 0.
	atimeact = 0.
	atimemin = -ahigh
	atimemax =  ahigh

        end

c******************************************************

        subroutine ptime_min_max

c sets time limits

        implicit none

	include 'timlim.h'

	integer iauto
	integer itanf,itend
	double precision atanf,atend

	double precision dgetpar

	iauto = nint(dgetpar('iauto'))
	if( iauto .le. 0 ) then
	  stop 'error stop  ptime_min_max: not ready for iauto=0'
	end if

	itanf = nint(dgetpar('itanf'))
	itend = nint(dgetpar('itend'))

	if( itanf .ne. -1 ) then
	  itmin = itanf
	  call ptime_i2a(itmin,atimemin)
	end if
	if( itend .ne. -1 ) then
	  itmax = itend
	  call ptime_i2a(itmax,atimemax)
	end if

	!write(6,*) '++++++++++++++++++++++++++++++++++++++++++'
	!write(6,*) atime0
	!write(6,*) itmin,itmax,atimemin,atimemax
	!write(6,*) '++++++++++++++++++++++++++++++++++++++++++'

	itfreq = nint(dgetpar('nout'))
        if( itfreq .eq. 0 ) itfreq = 1

	atanf = dgetpar('atanf')
	atend = dgetpar('atend')

	if( atanf .gt. 0. ) atimemin = atanf
	if( atend .gt. 0. ) atimemax = atend

	write(6,*)
	write(6,*) 'Using time parameters: ',itmin,itmax,itfreq
	write(6,*) 'absolute time: ',atimemin,atimemax
	write(6,*)

        end

c******************************************************

	subroutine ptime_info

	implicit none

	include 'timlim.h'

	write(6,*) '---------- ptime_info -----------'
	write(6,*) 'itmin,itmax: ',itmin,itmax
	write(6,*) 'itact: ',itact
	write(6,*) 'itfreq: ',itfreq
	write(6,*) 'nrec: ',nrec
	write(6,*) 
	write(6,*) 'atime0: ',atime0
	write(6,*) 'atimeact: ',atimeact
	write(6,*) 'atmin,atmax: ',atimemin,atimemax
	write(6,*) '---------------------------------'

	end

c******************************************************
c******************************************************
c******************************************************

	subroutine ptime_set_date_time(date,time)

	implicit none

	integer date,time

	include 'timlim.h'

	atime0 = 0.

	if( date > 0 ) then
          call dtsini(date,time)
	  call dts_to_abs_time(date,time,atime0)
	end if

	end

c******************************************************
c******************************************************
c******************************************************

	subroutine ptime_set_itime(it)

	implicit none

	integer it

	include 'timlim.h'

	atimeact = atime0 + it

	end

c******************************************************

	subroutine ptime_set_atime(atime)

	implicit none

	double precision atime

	include 'timlim.h'

	atimeact = atime

	end

c******************************************************

	subroutine ptime_set_dtime(dtime)

	implicit none

	double precision dtime

	include 'timlim.h'

	atimeact = atime0 + dtime

	end

c******************************************************

	subroutine ptime_get_itime(it)

	implicit none

	integer it

	include 'timlim.h'

	it = nint( atimeact - atime0 )

	end

c******************************************************

	subroutine ptime_get_atime(atime)

	implicit none

	double precision atime

	include 'timlim.h'

	atime = atimeact

	end

c******************************************************

	subroutine ptime_i2a(it,atime)

	implicit none

	integer it
	double precision atime

	include 'timlim.h'

	atime = atime0 + it

	end

c******************************************************
c******************************************************
c******************************************************

        function ptime_ok()

c is time ok?

        implicit none

	logical ptime_ok

	include 'timlim.h'

	integer it

        double precision atimeold
	save atimeold

	integer icall
	save icall
	data icall /0/

	if( icall .eq. 0 ) then
	  icall = icall + 1
	  nrec = nrec + 1
	  atimeold = atimeact
	end if

	if( atimeact .ne. atimeold ) then !increase for new time
	  nrec = nrec + 1
	  atimeold = atimeact
	end if

	it = atimeact - atime0
c        write(6,*) 'ptime_ok: ',nrec,itfreq,it,atimeact

	ptime_ok = .false.

	if( atimeact < atimemin ) return
	if( atimeact > atimemax ) return
	
	if( itfreq .gt. 0 .and. mod(nrec,itfreq) .eq. 0 ) then
	  ptime_ok = .true.
	else if( itfreq .eq. 0 ) then
	  ptime_ok = .true.
	else if( itfreq .lt. 0 .and. mod(nrec-1,-itfreq) .eq. 0 ) then
	  ptime_ok = .true.
	end if

	end

c******************************************************

        function ptime_end()

c is time over max limit?

        implicit none

	logical ptime_end

	include 'timlim.h'

	ptime_end = atimeact .gt. atimemax

	end

c******************************************************
c******************************************************
c******************************************************

