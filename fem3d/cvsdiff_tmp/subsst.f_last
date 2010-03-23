c
c $Id: subsst.f,v 1.3 2000/03/02 09:31:29 georg Exp $
c
c time scaling routines
c
c contents :
c
c function tconvs(itscal)		shell for tconv
c function tconv(scale1,scale2)		returns factor to scale time value
c
c*************************************************************
c
        function tconvs(itscal)
c
c shell for tconv
c
c itscal	switch to find out in what time unit
c		...has to be converted; it is assumed, that
c		...the original time unit is seconds
c tconv		factor time values have to be multiplied with
c		...(on error 1. is returned)
c
c possible values for itscal are :
c			0: sec   1: min    2:hour
c			3: day   4: month  5:year
c
        implicit none
c
	real tconvs
        integer itscal
        real tconv
c
	character*3 unit
c
	if(itscal.eq.0) then
		unit='sec'
	else if(itscal.eq.1) then
		unit='min'
	else if(itscal.eq.2) then
		unit='hou'
	else if(itscal.eq.3) then
		unit='day'
	else if(itscal.eq.4) then
		unit='mon'
	else if(itscal.eq.5) then
		unit='yea'
	end if
c
	tconvs=tconv('sec',unit)
c
	return
	end
c
c*************************************************************
c
        function tconv(scale1,scale2)
c
c returns factor to scale time value to an other unit
c
c scale1	unit of time
c scale2	unit to be converted in
c tconv		factor time values have to be multiplied with
c		...(on error 1. is returned)
c
c possible values for scale1/2 are :
c		sec,min,hour,day,month,year
c		sek,min,stunde,tag,monat,jahr
c		sec,min,ora,giorno,mese,anno
c abbreviations are possible (m=min,s=sec)
c
        implicit none
c
	real tconv
	character*(*) scale1,scale2
c
	character*2 scale
        integer i,iscale(2)
        real rscale(6)
c
	data rscale /1.,60.,3600.,86400.,2592000.,933120000./
c
	do i=1,2
c
	if(i.eq.1) then
		scale=scale1(1:2)
	else
		scale=scale2(1:2)
	end if
c
	call uplow(scale,'low')
c
	if(scale(1:1).eq.'s') then
		if(scale(2:2).eq.'t') then
			iscale(i)=3
		else
			iscale(i)=1
		end if
	else if(scale(1:1).eq.'m') then
		if(scale(2:2).eq.'e'.or.scale(2:2).eq.'o') then
			iscale(i)=5
		else
			iscale(i)=2
		end if
	else if(scale(1:1).eq.'y'.or.scale(1:1).eq.'j'
     +			.or.scale(1:1).eq.'a') then
		iscale(i)=6
	else if(scale(1:1).eq.'d'.or.scale(1:1).eq.'t'
     +			.or.scale(1:1).eq.'g') then
		iscale(i)=4
	else if(scale(1:1).eq.'h'.or.scale(1:1).eq.'o') then
		iscale(i)=3
	else
		iscale(i)=0
	end if
c
	end do
c
	if(iscale(1).ne.0.and.iscale(2).ne.0
     +				.and.iscale(1).ne.iscale(2)) then
		tconv=rscale(iscale(1))/rscale(iscale(2))
	else
		tconv=1.
	end if
c
	return
	end
