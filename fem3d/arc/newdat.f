
c*************************************************************
c
c absolute time is given in days and secs
c from 1/1/1 AD. 
c
c*************************************************************

	subroutine aform(days,secs,length,line)

c format time and date for output

	implicit none

	integer days,secs
	integer length
	character*(*) line

	integer l

	line = ' '

	call sform(secs,l,line)
	length = l + 1
	call dform(days,l,line(length+1:))

	length = length + l

	end

c*************************************************************

	subroutine dform(days,length,line)

c format date

	implicit none

	integer days
	integer length
	character*(*) line

	integer year,month,day
	integer i

	call upkdat(year,month,day,days)

	length=10
	write(line,'(i2,a1,i2,a1,i4)') day,'-',month,'-',year
	do i=1,length
	  if(line(i:i).eq.' ') line(i:i)='0'
	end do

	end

c*************************************************************

	subroutine sform(secs,length,line)

c format time

	implicit none

	integer secs
	integer length
	character*(*) line

	integer hour,min,sec
	integer i

	call upktim(hour,min,sec,secs)

	length=8
	write(line,'(i2,a1,i2,a1,i2)') hour,':',min,':',sec
	do i=1,length
	  if(line(i:i).eq.' ') line(i:i)='0'
	end do

	end

c*************************************************************
c*************************************************************
c*************************************************************

	function jdmon(year,month)

	implicit none

	integer jdmon
	integer year,month

	logical bises

	integer jmm(0:12)
	save jmm
	data jmm /0,31,59,90,120,151,181,212,243,273,304,334,365/

	jdmon=jmm(month)

	if( month .ge. 2 .and. bises(year) ) then
	  jdmon = jdmon + 1
	end if

	end

c*************************************************************

	function idmon(year,month)

	implicit none

	integer idmon
	integer year,month

	logical bises

	integer imm(0:12)
	save imm
	data imm /0,31,28,31,30,31,30,31,31,30,31,30,31/

	idmon=imm(month)

	if( month .eq. 2 .and. bises(year) ) then
	  idmon = idmon + 1
	end if

	end

c*************************************************************

	function bises(year)

c true if year is bisestial

	logical bises
	integer year

	if(  ( mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0 )
     +                  .or. mod(year,400) .eq. 0 ) then
		bises = .true.
	else
		bises = .false.
	end if

	end
 
c*************************************************************
c*************************************************************
c*************************************************************

	subroutine adjdat(days,secs)

c adjusts days and secs [0,86400]

	implicit none

	integer days,secs

	days = days + secs/86400
	secs = mod(secs,86400)

	if( secs .lt. 0 ) then
	  secs = secs + 86400
	  days = days - 1
	end if

	end

c*************************************************************

	subroutine addsec(days,secs,seconds)

c adds seconds to days and secs

	implicit none

	integer days,secs,seconds

	secs = secs + seconds

	call adjdat(days,secs)

	end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine j2date(year,month,day,jdays)

c computes month and day from year and julian days

	implicit none

	integer year,month,day
	integer jdays

	integer i
	integer jdmon

	do i=1,12
	  if( jdmon(year,i) .ge. jdays ) goto 1
	end do
	write(6,*) 'jdays,year: ',jdays,year
	stop 'error stop j2date: jdays'
    1	continue

	month = i
	day = jdays - jdmon(year,month-1)

	end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine nextday(year,month,day)

c adds one day to date

	implicit none

	integer year,month,day
	integer days

	call setday(year,month,day,days)
	days = days+1
	call getday(year,month,day,days)

	end

c*************************************************************

	subroutine getsecs(hour,min,sec,secs)

c splits secs into hour, min, sec

	implicit none

	integer hour,min,sec
	integer secs

	min = secs / 60
	sec = secs - min * 60

	hour = min / 60
	min = min - hour * 60

	end

c*************************************************************

	subroutine getday(year,month,day,days)

c splits days into year, month, day

	implicit none

	integer year,month,day
	integer days

	integer iy,auxdays
	integer idmon

c	estimate year (using average days in 400 years)

	iy = (days*400) / (365*400+97) + 2

c	correct year

	auxdays = days + 1
	do while( auxdays .gt. days )
	  iy = iy - 1
	  auxdays = 365*iy + iy/4 - iy/100 + iy/400
	end do

c	year found, auxdays days left

	year = iy + 1
	auxdays = days - auxdays

c	find month

	month = 0
	do while( auxdays .ge. 0 )
	  month = month + 1
	  auxdays = auxdays - idmon(year,month)
	end do

c	find day

	day = auxdays + idmon(year,month) + 1
	
	end

c*************************************************************

	subroutine setsecs(hour,min,sec,secs)

c packs hour, min, sec into secs

	implicit none

	integer hour,min,sec
	integer secs

	secs = sec + 60 * ( min + 60 * hour )

	end

c*************************************************************

	subroutine setday(year,month,day,days)

c packs year, month, day into days

	implicit none

	integer year,month,day
	integer days

	integer iy
	integer it,i
	integer jdmon
	logical bises

	iy = year - 1
	days = 365*iy + iy/4 - iy/100 + iy/400

	days = days + jdmon(year,month-1)

	days = days + day - 1

	end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine upktim(hour,min,sec,time)

c splits time [HHMMSS] into hour, min, sec

	implicit none

	integer hour,min,sec
	integer time

	min = time / 100
	sec = time - min * 100

	hour = min / 100
	min = min - hour * 100

	end

c*************************************************************

	subroutine upkdat(year,month,day,date)

c splits date [YYYYMMDD] into year, month, day

	implicit none

	integer year,month,day
	integer date

	month = date / 100
	day = date - month * 100

	year = month / 100
	month = month - year * 100
	
	end

c*************************************************************

	subroutine pktim(hour,min,sec,time)

c packs hour, min, sec into time [HHMMSS]

	implicit none

	integer hour,min,sec
	integer time

	time = sec + 100 * ( min + 100 * hour )

	end

c*************************************************************

	subroutine pkdat(year,month,day,date)

c packs year, month, day into date [YYYYMMDD]

	implicit none

	integer year,month,day
	integer date

	date = day + 100 * ( month + 100 * year )

	end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine tstdat

c tests routines

	implicit none

	integer days,newdays
	integer daymax
	integer year,month,day

	daymax = 2000 * 366

	do days = 1,daymax

	  call getday(year,month,day,days)
	  call setday(year,month,day,newdays)

	  write(6,*) days,newdays,year,month,day

	  if( days .ne. newdays ) stop 'error stop'

	end do

	end

c*************************************************************

c	program mtstdat
c	call tstdat
c	end

c*************************************************************
c*************************************************************
c*************************************************************

