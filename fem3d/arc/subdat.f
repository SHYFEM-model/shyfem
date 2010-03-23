c
c $Id: subdat.f,v 1.4 2000/03/02 09:31:23 georg Exp $
c
c date and time routines
c
c contents :
c
c subroutine absdat(iday,imon,iyear,it)
c               computes date from total days from 1.1.1
c subroutine adddat(idday,idmon,idyear,it)
c               adds days/months/years to total days it
c integer function itday(iday,imon,iyear)
c               computes total days from 1.1.1 to iday.imon.iyear
c integer function idy(imon,iyear)
c               returns number of days in month imon of year iyear
c integer function ibises(iyear)
c               decides if a year is bisestial or not
c real function ryear(it)
c               computes years from relative date it
c subroutine gdate(thour,iday,imon,iyear,imino,iorao,idayo,imono,iyearo)
c               computes absolute date given an initial date and the
c               ...number of hours from this date
c subroutine giorno(idayto,iyear,imonth,iday)
c               computes month and day from a given number of 
c               ...total days in a year
c real function hfroms(isec)
c		computes hours from seconds
c subroutine sectod(is,iday,ihour,imin,isec)
c		computes days from seconds
c
c revision log :
c
c 23.03.1998	ggu	hfroms, sectod integrated
c
c notes :
c
c this file is not used at all and can be removed
c
c***********************************************************
c
	subroutine absdat(iday,imon,iyear,it)
c
c computes date from total days from 1.1.1
c
c iday          day             (return value)
c imon          month           (return value)
c iyear         year [19..]     (return value)
c it            total number of days from 1.1.1 (input value)
c
c arguments
	integer iday,imon,iyear,it
c local variables
	integer iy,itt
c functions
	integer idy
c
c error check
	if(it.lt.0) goto 99
c find years
	iy=(it*400)/(365*400+97)+2
c correct years
    1   continue
	iy=iy-1
	itt=365*iy+iy/4-iy/100+iy/400
	if(itt.gt.it) goto 1
c years found, itt days left
	iyear=iy+1
	itt=it-itt
c find months
	imon=0
    2   continue
	imon=imon+1
	itt=itt-idy(imon,iyear)
	if(itt.ge.0) goto 2
c find days
	iday=itt+idy(imon,iyear)+1
c
	return
   99   continue
	write(6,*) 'error in total number of days : ',it
	stop 'error stop : absdat'
	end
c
c***********************************************************
c
	subroutine adddat(idday,idmon,idyear,it)
c
c adds days/months/years to total days it
c
c idday         delta day
c idmon         delta month
c idyear        delta year [19..]
c it            total number of days from 1.1.1 (input + output value)
c
c arguments
	integer idday,idmon,idyear,it
c local variables
	integer iday,imon,iyear
c functions
	integer itday
c
c get absolute date
	call absdat(iday,imon,iyear,it)
c add delta values
	iyear=iyear+idyear
	imon=imon+idmon
	iday=iday+idday
c adjust year and month
	iyear=iyear+(imon-1)/12
	imon=mod(imon-1,12)+1
c sum days
	it=itday(1,imon,iyear)+iday-1
c
	return
	end
c
c***********************************************************
c
	function itday(iday,imon,iyear)
c
c computes total days from 1.1.1 to iday.imon.iyear
c
c iday          day
c imon          month
c iyear         year [19..]
c itday         number of days (return value)
c
	integer itday
c arguments
	integer iday,imon,iyear
c local variables
	integer it,iy,i
c functions
c	integer ibises,idy
	integer idy
c
c error check
	if(iyear.le.0) goto 99
	if(imon.lt.1.or.imon.gt.12) goto 99
	if(iday.gt.idy(imon,iyear)) goto 99
	if(iday.lt.1) goto 99
c sum up years
	iy=iyear-1
	it=365*iy+iy/4-iy/100+iy/400
c other possibility to sum up years (saver but slower)
c       it=0
c       do 11 i=1,iyear-1
c       it=it+365+ibises(i)
c   11  continue
c sum up months
	do 12 i=1,imon-1
	it=it+idy(i,iyear)
   12   continue
c sum up days
	itday=it+iday-1
c
	return
   99   continue
	write(6,*) 'error in date : ',iday,imon,iyear
	stop 'error stop : itday'
	end
c
c**********************************************************
c
	function idy(imon,iyear)
c
c returns number of days in month imon of year iyear
c
c imon          month [1..12]
c iyear         year [1...]
c idy           number of days in month (return value)
c
	integer idy
c arguments
	integer imon,iyear
c local variables
	integer idyy(12)
c functions
	integer ibises
c
	save idyy
	data idyy /31,28,31,30,31,30,31,31,30,31,30,31/
c
c error check
	if(imon.lt.1.or.imon.gt.12) goto 99
c find days
	if(imon.eq.2) then
		idy=idyy(2)+ibises(iyear)
	else
		idy=idyy(imon)
	end if
c
	return
   99   continue
	write(6,*) 'error in month : ',imon
	stop 'error stop : idy'
	end
c
c**********************************************************
c
	function ibises(iyear)
c
c decides if a year is bisestial or not
c
c iyear         year [19..]
c ibises        1 if bisestial   0 if not
c
	integer ibises
c arguments
	integer iyear
c
	if(  ( mod(iyear,4).eq.0 .and. mod(iyear,100).ne.0 )
     +                  .or. mod(iyear,400).eq.0) then
		ibises=1
	else
		ibises=0
	end if
c
	return
	end
c
c********************************************************
c
	function ryear(it)
c
c computes years from relative date it
c
c it            days from 1.1.1
c ryear         number of years (return value)
c
	real ryear
c arguments
	integer it
c local variables
	integer iday,imon,iyear
	integer it1,it2
c function
	integer itday
c
	call absdat(iday,imon,iyear,it)
	it1=itday(1,1,iyear)
	it2=itday(31,12,iyear)
c
	ryear=iyear+(it-it1)/float(it2-it1+1)
c
	return
	end
c
c********************************************************************
c
	subroutine gdate(thour,iday,imon,iyear
     +                          ,imino,iorao,idayo,imono,iyearo)
c
c computes absolute date given an initial date and the
c number of hours from this date
c
c input data :
c thour         number of hours from initial date (real)
c iday          initial day
c imon          initial month
c iyear         initial year
c output data :
c imino         minute 
c iorao         hour
c idayo         day
c imono         month
c iyearo        year
c
	implicit none
c
c arguments
	real thour
	integer iday,imon,iyear
	integer imino,iorao,idayo,imono,iyearo
c local
	integer it
c functions
	integer idy
c
	it=thour
	imino=(thour-it)*60.001
	iorao=mod(it,24)
c
	idayo=iday+it/24
	imono=imon
	iyearo=iyear
c
	do while(idayo.gt.idy(imono,iyearo))
		idayo=idayo-idy(imono,iyearo)
		imono=imono+1
		if(imono.gt.12) then
			imono=1
			iyearo=iyearo+1
		end if
	end do
c
	return
	end
c
c*********************************************************
c
	subroutine giorno(idayto,iyear,imonth,iday)
c
c computes month and day from a given number of total
c ...days in a year
c
c idayto        number of total days in year [1...366]
c iyear         year
c imonth        month [1...12]
c iday          day [1...31]
c
	implicit none
c
c arguments
	integer idayto,iyear,imonth,iday
c local
	integer id,i
c functions
	integer idy
c
	id=idayto
c
	i=1
	do while(i.le.12.and.id.gt.idy(i,iyear))
	  id=id-idy(i,iyear)
	  i=i+1
	end do
c
	if(i.gt.12) then
		write(6,*) 'error in total number of days : ',idayto
		stop 'error stop : giorno'
	end if
c
	imonth=i
	iday=id
c
	return
	end
c
c***********************************************************
c
	function hfroms(isec)
c
c computes hours from seconds
c
c isec		seconds
c hfroms	hours (real)
c
	implicit none

	real hfroms
c
c arguments
	integer isec
c
	if(mod(isec,3600).eq.0) then
		hfroms=(isec/3600)
	else
		hfroms=isec/3600.
	end if
c
	return
	end
c
c***********************************************************
c
	subroutine sectod(is,iday,ihour,imin,isec)
c
c computes days from seconds
c (seconds to days)
c
c is		seconds (input)
c iday		days (output)
c ihour		hours (output)
c imin		minutes (output)
c isec		seconds (output)
c
	implicit none
c
c arguments
	integer is,iday,ihour,imin,isec
c local
	integer iss
c
	iss=is
	isec=iss-(iss/60)*60
	iss=iss/60
	imin=iss-(iss/60)*60
	iss=iss/60
	ihour=iss-(iss/24)*24
	iday=iss/24
c
	return
	end
c
c***********************************************************
c
