!
! $Id: subsss.f,v 1.11 2009-02-04 15:26:54 georg Exp $
!
! general utility routines
!
! contents :
!
! function izahl(z,ndec)		computes ciphers of number
! function laezi(l,ioff)		computes length of number
! function iscan(line,ioff,f)		converts alphanumeric text to numbers
! function iantw(l)			gets answer from terminal
! function inquir(l,f)			gets numbers from terminal
! function getnum(prompt)		gets number from terminal
! function igetxt(prompt,text)		gets text from terminal
! function ialfa(zahl,lh,ndec,mode)	converts number into alphanumeric chars
! subroutine uplow(text,type)		translates text to upper/lower case 
! subroutine tablnc(linold,linnew)	converts tabs into blanks
! subroutine prilin(line,iunit)		writes a line without trailing blanks
! function ichanm(line)			computes length of line
! function ichafs(line)			first occurrence of non blank character 
! function icompr(line)			eliminates blank characters from line
! function ifndch(line,iocc,cha)	finds the iocc occurrence of cha 
! function ifndoc(line,string)		finds number of occurences of string
! function xzentr(x,height,f,ndec)	funtion centers number
! function ideflt(k,text)		writes default and gets value (integer)
! function fdeflt(r,text)		writes default and gets value (double precision)
! function iround(f)			rounds double precision value
! function itypch(char)			returns typ of character
! subroutine tabula(tab)		creates tab character
! function rnext(r,mode)		finds closest value to r
! function rnexta(r,mode)		finds closest value to r (absolute)
! function istell(r)			computes exponent of r
! subroutine scplo(xmin,xmax,ymin,ymax)	scales plot
! subroutine swapr(a,b)			swaps values
! function rlen(x1,y1,x2,y2)		length of line
!
! uplow,itypch,tabula depend on the ASCII character set
!
! revision log :
!
! 01.06.1998	ggu	new routines iscan... into own file (subscn.f)
! 17.06.1998	ggu	old routines with new name (iscan0,ialfa0)
! 12.02.1999	ggu	new routine triml
! 26.01.2009	ggu	minor changes to avoid compiler warnings
! 16.02.2011	ggu	new routine trimline()
! 30.05.2014	ggu	new routine rnextsub()
! 06.05.2015	ggu	new routine logvals() for logarithmic values
! 14.09.2015	ggu	new routine i2s0()
!
!***********************************************************
        module utility
!***********************************************************

        contains

!***********************************************************
	function izahl(z,ndec)
!
! computes ciphers of number
!
! z		number for which ciphers have to be determined
! ndec		positions after the decimal point
!		(-1 : no decimal point)
! izahl		ciphers of number z with minus sign and
!		...decimal point if necessary (return value)
!
	if(ndec.lt.-1) then
		ndech=iabs(ndec)-1
		zh=z
	else if(ndec.ge.1) then
		ndech=ndec
		zh=z
	else
		ndech=0
		zh=z
	end if
!
	iz=abs(zh)+0.01
!
	if(iz.eq.0) then
		istel=1
	else
		istel=0
	end if
!
	do while (iz.gt.0)
	iz=iz/10
	istel=istel+1
	end do
!
	if(z.lt.0.) istel=istel+1
	if(ndec.ge.0) then
		istel=istel+1+ndec
	else if(ndec.lt.-1) then
		if(istel+1.lt.ndech) then
			istel=istel+1
		else
			istel=ndech
		end if
	end if
!
	izahl=istel
!
	return
	end
!
!**********************************
!
	function laezi(l,ioff)
!
! computes length of number
!
! l		text where the value is stored
! ioff		offset in l from where on the value is stored
! laezi		length of number (return value)
!
	character*(*) l
!
	character*1 blank,comma,tab,cha
	data blank,comma /' ',','/
!
	call tabula(tab)
!
	do j=ioff,len(l)
	cha=l(j:j)
	if(cha.eq.blank.or.cha.eq.comma.or.cha.eq.tab) goto 1
	end do
!
    1	continue
	laezi=j-ioff
!
	return
	end
!
!**************************************
!
	function iscan0(line,ioff,f)
!
! converts alphanumeric text to numbers
!
! line		text to be translated
! ioff		offset in text from where on values have to be read
! f		array with translated values (return value)
! iscan		total number of values read (return value)
!		-1 : error reading line
!
! as separators  blank, comma and tab can be used
! inbetween two commas, a value of 0. is assumed
!
	character*(*) line
!
	character*10 lm
	character*1 lh,blank,tab,comma,plus,minus,dot,eee
	double precision f(1)
	data lm /'1234567890'/
	data blank,comma,plus,minus,dot,eee /' ',',','+','-','.','e'/
!
	call tabula(tab)
!
	ianz=0		!total number read
	imod=0		!1:in number
	ikom=1		!1:comma read
	ipun=0		!1:decimal point read
	iexp=0		!1:exponential part
	izei=1		!sign of number
	izeiex=0	!sign of exponential
!
	ff=1.
	ffac=1.
	kexp=0
!
	length=len(line)
!
	do j=ioff,length
!
	lh=line(j:j)
	call uplow(lh,'low')
!
	if(lh.eq.blank.or.lh.eq.tab) then	!blank
		if(imod.eq.1) then	!end of number found
			f(ianz)=ff*izei*(10.**(kexp*izeiex))
			imod=0
		end if
	else if(lh.eq.comma) then		!comma
		if(imod.eq.1) then	!end of number found
			f(ianz)=ff*izei*(10.**(kexp*izeiex))
		else if(ikom.eq.1) then	!assume 0 between two commas
			ianz=ianz+1
			f(ianz)=0.
		end if
		ikom=1
		imod=0
	else if(lh.eq.eee) then			!exponential
		if(imod.eq.0) ianz=ianz+1	!number starts with 'e'
		iexp=1
		imod=1
		ikom=0
	else					!number
		if(imod.eq.0) then		!start reading
			ianz=ianz+1
			imod=1
			ipun=0
			ikom=0
			iexp=0
			izeiex=0
			kexp=0
			ff=0.
			ffac=1.
			if(lh.eq.plus) then
				izei=+1
				goto 1
			else if(lh.eq.minus) then
				izei=-1
				goto 1
			else
				izei=+1
			end if
		end if
!
		if(iexp.eq.1.and.izeiex.eq.0) then	!exponential
			if(lh.eq.plus) then
				izeiex=+1
				goto 1
			else if(lh.eq.minus) then
				izeiex=-1
				goto 1
			else
				izeiex=+1
			end if
		end if
			
		if(lh.eq.dot) then		!point
			if(ipun.eq.1) goto 99
			if(iexp.eq.1) goto 99
			ipun=1
			ffac=1.		!prepare for decimal part
			fh=0.		!aux variable
			goto 1
		end if
!
		fh=0.				!find cipher
		do i=1,10
		if(lm(i:i).eq.lh) fh=i
		end do
		if(fh.eq.0.) goto 99
		if(fh.eq.10.) fh=0.
!
		if(iexp.eq.1) then	!cypher for exponential
			kexp=10*kexp+int(fh)
		else			!cypher for number
			if(ipun.eq.0) then	
				ff=10.*ff+fh
			else
				ffac=ffac/10.
				ff=ff+ffac*fh
			end if
		end if
!
    1		continue
	end if
!
	end do
!
	if(imod.eq.1) then
		f(ianz)=ff*izei*(10.**(kexp*izeiex))
	end if
!
	iscan0=ianz
!
	return
!
   99	continue
	iscan0=-1
!
	return
	end
!
!******************************************
!
	function iantw(l)
!
! gets answer from terminal
!
! l		text written to terminal
! iantw		1 : true
!		0 : false
!
! the following answers can be given as true or false :
!
!		y j s t 1	==>	true
!		n f 0 'blank'	==>	false
!
	character*(*) l
!
	character l1*1
	data net,nat /5,6/
!
    1	continue
!
	write(nat,*) l
	read(net,1000) l1
 1000	format(a)
!
	call uplow(l1,'low')
!
	if(l1.eq.'y'.or.l1.eq.'j'.or.l1.eq.'s' &
     &			.or.l1.eq.'t'.or.l1.eq.'1') then
		iantw=1
	else if(l1.eq.'n'.or.l1.eq.'f'.or.l1.eq.'0' &
     &			.or.l1.eq.' ') then
		iantw=0
	else
		write(nat,*) 'Incorrect answer. Try again.'
		goto 1
	end if
!
	return
	end

!********************************************

	function inquire_numbers(l,f)

! gets numbers from terminal (see also iscan)
!
! do not use -> no way to know if we are out of bounds
!
! l		text written to terminal
! f		array in which the values are stored (return value)
! inquire_numbers	total number of values read in
        use convert, only: iscan

	implicit none

	integer inquire_numbers
	character*(*) l
	double precision f(1)

	character*80 lh
	integer net,nat
	data net,nat /5,6/

	lh=' '

	write(nat,1000) l
 1000	format(1x,a)

	read(net,2000) lh
 2000	format(a)

	inquire_numbers=iscan(lh,1,f)

	end

!********************************************
!
	function getnum(prompt)
!
! gets number from terminal (see also iscan)
!
! prompt	text written to terminal
! getnum	value of number 
!
        use convert, only: iscan

        implicit none
        double precision getnum
	character*(*) prompt
!
	character*80 lh
	double precision f(100)
        integer net,nat
	data net,nat /5,6/
!
    1	continue
!
	write(nat,1000) prompt
 1000	format(1x,a)
!
	lh=' '
	read(net,2000) lh
 2000	format(a)
!
	if(iscan(lh,1,f).ne.1) then
		write(nat,1000) 'Erroneous input. Repeat.'
		goto 1
	end if
!
	getnum=f(1)
!
	return
	end
!
!********************************************
!
	function igetxt(prompt,text)
!
! gets text from terminal
!
! prompt	text written to terminal
! text		text returned
! igetxt	number of characters read
!
	character*(*) prompt,text
!
	data net,nat /5,6/
!
	write(nat,1000) prompt
 1000	format(1x,a)
!
	read(net,2000) text
 2000	format(a)
!
	igetxt=ichanm(text)
!
	return
	end
!
!********************************************
!
	function ialfa0(zahl,lh,ndec,mode)
!
! converts double precision number into alphanumeric characters
! fills zahl with '*', if dimension of lh is
! ...to small for zahl
!
! zahl		number to be converted
! lh		character string where zahl is written to
!		...(is filled with '*' if dimension of lh
!		...is to small)
! ndec		ciphers after decimal point (-1 for no decimal point)
! mode		-1	left justified
!		 0	centred
!		 1	right justified
! ialfa		total number of characters written to lh
!
	character*(*) lh
!
	character*10 lm
	data lm /'1234567890'/
!
	lentxt=len(lh)
!
	istel=izahl(zahl,ndec)
!
	ialfa=istel
	zahlh=zahl
!
	if(lentxt.lt.istel) then
		do i=1,lentxt
		lh(i:i)='*'
		end do
		ialfa=lentxt
		ialfa0=lentxt
		return
	end if
!
	ncont=0
	if(zahlh.lt.0) then
		ncont=1
		lh(1:1)='-'
		istel=istel-1
		zahlh=-zahlh
	end if
!
! integer und fraction
!
	izi=istel-ndec-1
	if(ndec.gt.0) then
		izf=ndec
		ifact=10**ndec
		zahlt=zahlh*ifact+.5
		izahli=zahlt/ifact
		izahlf=(zahlh-izahli)*ifact
!		ifact=10**ndec
!		izf=ndec
!		izahlt=zahlh*ifact+.5
!		izahlf=mod(izahlt,ifact)
!		izahli=izahlt/ifact
	else
		izahlf = 0
		ifact=1
		izf=0
		izahli=zahlh+.5
	end if
!
! integer
!
	do i=izi,1,-1
	izv=mod(izahli,10)
	if(izv.eq.0) izv=10
	lh(ncont+i:ncont+i)=lm(izv:izv)
	izahli=izahli/10
	end do
!
	ncont=ncont+izi
!
! fraction
!
	if(ndec.ge.0) then
		ncont=ncont+1
		lh(ncont:ncont)='.'
		do i=izf,1,-1
		izv=mod(izahlf,10)
		if(izv.eq.0) izv=10
		lh(ncont+i:ncont+i)=lm(izv:izv)
		izahlf=izahlf/10
		end do
!
		ncont=ncont+izf
	end if
!
	do i=ncont+1,lentxt
	lh(i:i)=' '
	end do
!
! justify text
!
	lenlh=lentxt
	ial=ialfa
!
	if(mode.gt.0) then
		ndif=lenlh-ial
	else if(mode.eq.0) then
		ndif=(lenlh-ial)/2
	else
		ndif=0
	end if
!
	do i=ial,1,-1
	   lh(i+ndif:i+ndif)=lh(i:i)
	end do
	if(ndif.ge.1) lh(1:ndif)=' '
	if(ial+ndif.lt.lenlh) lh(ial+ndif+1:lenlh)=' '
!
	ialfa0 = ialfa
!
	return
	end
!
!***************************************************
!
	subroutine uplow(text,type)
!
! translates text to upper/lower case characters
! ...depending on type
!
! works only for ASCII set of characters
!
! text		text to be transformed
! type		'up'  :  transform to uppercase
!		'low' :  transform to lowercase
!
	character*(*) text,type
!
	logical bupper
	character*1 char
!
	iascii=ichar(type(1:1))
!
	if(iascii.eq.85.or.iascii.eq.117) then
		bupper=.true.
	else if(iascii.eq.76.or.iascii.eq.108) then
		bupper=.false.
	else
		return
	end if
!
	do i=1,len(text)
!
	iascii=ichar(text(i:i))
	if(bupper) then
		if(iascii.ge.97.and.iascii.le.122) then
			iascii=iascii-32
			text(i:i)=char(iascii)
		end if
	else 
		if(iascii.ge.65.and.iascii.le.90) then
			iascii=iascii+32
			text(i:i)=char(iascii)
		end if
	end if
!
	end do
!
	return
	end
!
!***********************************************************
!
	subroutine tablnc(linold,linnew)
!
! converts tabs into blanks
!
! linold	old line with tabs
! linnew	new converted line
!
! itab		number of blanks equals one tab (for fortran minimum 6)
!
	character*(*) linold,linnew
!
	character*1 tab
	data itab   / 8 /
!
	call tabula(tab)
	lengo=len(linold)
	lengn=len(linnew)
!
	iold=1
	inew=1
!
	do while(iold.le.lengo.and.inew.le.lengn)
!
	if(linold(iold:iold).ne.tab) then
		linnew(inew:inew)=linold(iold:iold)
		iold=iold+1
		inew=inew+1
	else
		ianf=inew
		iend=itab*(1+((inew-1)/itab))
		if(iend.gt.lengn) iend=lengn
		do i=ianf,iend
		   linnew(inew:inew)=' '
		   inew=inew+1
		end do		
		iold=iold+1
	end if
!
	end do
!
	return
	end
!
!*******************************************************
!
	subroutine prilin(line,iunit)
!
! writes a line of text without trailing blanks
!
! line		line of text
! iunit		unit number of file
!
	character*(*) line
!
	iend=ichanm(line)
	write(iunit,'(a)') line(1:iend)
!
	return
	end
!
!*********************************************
!
	function ichanm(line)
!
! computes length of line without trailing blanks
!
! line		line of text
! ichanm	length of line (return value)
!		... 0 : line is all blank
!
	character*(*) line
!
	character*1 blank,tab
	data blank /' '/
!
	call tabula(tab)
!
	ndim=len(line)
!
	do i=ndim,1,-1
	if(line(i:i).ne.blank.and.line(i:i).ne.tab) goto 1
	end do
!
    1	continue
	ichanm=i
!
	return
	end
!
!*********************************************
!
	function ichafs(line)
!
! computes first occurrence of a non-blank character in line
!
! line		line of text
! ichafs	position of first non-blank character (return value)
!		... 0 : no non-blank character found
!
	character*(*) line
!
	character*1 blank,tab
	data blank /' '/
!
	call tabula(tab)
!
	ndim=len(line)
!
	do i=1,ndim
	if(line(i:i).ne.blank.and.line(i:i).ne.tab) goto 1
	end do
!
	i=0
!
    1	continue
	ichafs=i
!
	return
	end
!
!*********************************************
!
	function icompr(line)
!
! eliminates all blank characters from a character string
!
! line		line of text
! icompr	length of compressed character string (return value)
!
	character*(*) line
!
	character*1 blank,tab
	data blank /' '/
!
	call tabula(tab)
!
	ndim=len(line)
!
	ih=0
	do i=1,ndim
	if(line(i:i).ne.blank.and.line(i:i).ne.tab) then
		ih=ih+1
		line(ih:ih)=line(i:i)
	end if
	end do
!
	do i=ih+1,ndim
	line(i:i)=blank
	end do
!
	icompr=ih
!
	return
	end
!
!***************************************************
!
	function ifndch(line,iocc,cha)
!
! finds the iocc occurrence of cha in line from offset ioff on
!
! line		line of text
! iocc		occurrence of character
!		... >0 : iocc occurrence is found
!		...  0 : last occurrence is found
!		... -1 : last but one occurrence is found...
! cha		character(s) to be found
! ifndch	position of cha in line (return value)
!		...  0 : cha not found
!		... -1 : error in ioff or ndim
!
	character*(*) line,cha
!
	n=0
!
	ncha=len(cha)-1
	ndim=len(line)
!
	if(iocc.gt.0) then
		iend=ndim-ncha
		do i=1,iend
		if(line(i:i+ncha).eq.cha(1:1+ncha)) then
			n=n+1
			if(n.eq.iocc) then
				ifndch=i
				return
			end if
		end if
		end do
	else
		ianf=ndim-ncha
		do i=ianf,1,-1
		if(line(i:i+ncha).eq.cha(1:1+ncha)) then
			n=n+1
			if(1-n.eq.iocc) then
				ifndch=i
				return
			end if
		end if
		end do
	end if
!
	ifndch=0
!
	return
	end
!
!*************************************************************
!
	function ifndoc(line,string)
!
! finds number of occurences of string in line
!
	character*(*) line,string
!
	ilen=len(line)
	istr=len(string)-1
!
	n=0
	iend=ilen-istr
	do i=1,iend
	if(line(i:i+istr).eq.string(1:1+istr)) n=n+1
	end do
!
	ifndoc=n
!
	return
	end
!
!*************************************************************
!
	function xzentr(x,height,f,ndec)
!
! funtion centers number
!
! x		coordinate around which number has to be centered
! height	height of number (also width)
! f		number to be centered (double precision)
! ndec		cifres after decimal point (-1 : no decimal point)
! xzentr	starting coordinate of centered number
!
	istel=izahl(f,ndec)
!
	xzentr=x-istel*height*0.5
!
	return
	end
!
!************************************************************
!
	function ideflt(k,text)
!
! writes default and gets value (integer)
!
! k		default value (integer)
! text		text written to terminal
! ideflt	integer value chosen
!
        use convert, only: iscan

	character*(*) text
!
	character line*80
	double precision f(10)
!
    1	continue
!
	write(6,*) text,' (default = ',k,' )'
	read(5,'(a)') line
!
	ianz=iscan(line,1,f)
!
	if(ianz.eq.1) then
		ideflt=iround(f(1))
	else if(ianz.eq.0) then
		ideflt=k
	else if(ianz.gt.1) then
		write(6,*) 'Only one number at a time.'
		goto 1
	else
		write(6,*) 'Read error'
		goto 1
	end if
!
	return
	end
!
!************************************************************
!
	function fdeflt(r,text)
!
! writes default and gets value (double precision)
!
! r		default value (double precision)
! text		text written to terminal
! fdeflt	double precision value chosen
!
        use convert, only: iscan
	character*(*) text
        double precision fdeflt,r
!
	character line*80
	double precision f(10)
!
    1	continue
!
	write(6,*) text,' (default = ',r,' )'
	read(5,'(a)') line
!
	ianz=iscan(line,1,f)
!
	if(ianz.eq.1) then
		fdeflt=f(1)
	else if(ianz.eq.0) then
		fdeflt=r
	else if(ianz.gt.1) then
		write(6,*) 'Only one number at a time.'
		goto 1
	else
		write(6,*) 'Read error'
		goto 1
	end if
!
	return
	end
!
!***********************************************************
!
	function iround(f)
!
! rounds double precision value
!
! f		double precision value
! iround	rounded integer value
!
        integer iround
        double precision f

	if(f.ge.0) then
		iround=f+0.5
	else
		iround=f-0.5
	end if
!
	return
	end
!
!****************************************************
!
	function itypch(char)
!
! returns typ of character
!
! works only for ASCII set of characters
!
! char		character (*1) to be looked at
! itypch	typ of character
!		1 = numeric
!		2 = letter
!		3 = special character 	=+-*/(),.'"$_!:<>%&	!'
!		4 = non-fortran character
!
	character*1 char
!
	iascii=ichar(char)
!
! number
!
	if(iascii.ge.48.and.iascii.le.57) then
		itypch=1
		return
	end if
!
! letter
!
	if(iascii.ge.65.and.iascii.le.90) then
		itypch=2
		return
	end if
!
	if(iascii.ge.97.and.iascii.le.122) then
		itypch=2
		return
	end if
!
! special character
!
	if(iascii.ge.32.and.iascii.le.47.and.iascii.ne.35) then
		itypch=3
		return
	end if
!
	if(iascii.ge.58.and.iascii.le.62.and.iascii.ne.59) then
		itypch=3
		return
	end if
!
	if(iascii.eq.9.or.iascii.eq.95) then
		itypch=3
		return
	end if
!
! no fortran character
!
	itypch=4
!
	return
	end
!
!****************************************************************
!
	subroutine tabula(tab)
!
! creates tab character
!
! works only for ASCII set of characters
!
! tab		tab character (return value)
!
	character*1 tab,char
!
	tab=char(9)
!
	return
	end

!*********************************************************

	function rnext(r,mode)
!
! finds closest value to r
!
! r		to r the closest value has to be found
! mode		mode of operation
!		...  0  : r is returned (no change)
!		... >0  : the higer value is found (further from 0)
!		... <0  : the lower value is found (closer to 0)
!		... |1| : 1. 2. 2.5 5. 8.
!		... |2| : 1. 2. 5.
!		... |3| : 1. 2. 3. 4. 5. 8.
!		... |4| : 1. 2. 3. 4. 5. 6. 7. 8. 9.
! rnext		closest value found to r
!
! val		matrix containing the closest values to be used
!		...for each mode (in ascending order)
! nval		number of values to be used for each mode
!
! if r is too small, 0 is returned
! for negative r the lower value is the value closer to zero
!
	parameter (nmodim=4,nvadim=9)
	logical bhigh
	double precision val(nvadim,nmodim)
	integer nval(nmodim)
        double precision rnext,r
!
	data val / 1. , 2. , 2.5 , 5. , 8. ,           0.,0.,0.,0.       &
     &		,  1. , 2. , 5. ,                      0.,0.,0.,0.,0.,0. &
     &		,  1. , 2. , 3. , 4. , 5. , 8. ,       0.,0.,0.          &
     &		,  1. , 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9.            &
     &		 /
	data nval / 5 , 3 , 6 , 9 /
!
	if(mode.lt.0) then	!upper or lower value
		bhigh=.false.
		m=-mode
	else
		bhigh=.true.
		m=mode
	end if
!
	if(m.gt.nmodim.or.m.eq.0) then	!wrong mode, return r
		rnext=r
		return
	end if
!
	if(abs(r).lt.1.e-5) then	!r too small, return 0
		rnext=0.
		return
	end if
!
	if(r.lt.0.) then	!r negative --> make positive
		rin=-r
		sign=-1.
	else
		rin=r
		sign=1.
	end if
!
! find next value
!
	n=nval(m)
	iexpo=int(log10(rin))-2
	expo=10.**iexpo
	nn=1
	rr=val(nn,m)*expo
	rrold = rr

	do while(rr.lt.rin)
		rrold=rr
		nn=nn+1
		if(nn.gt.n) then
			nn=1
			expo=expo*10.
		end if
		rr=val(nn,m)*expo
	end do
!
	if(bhigh.or.rr.eq.rin) then	!return next value
		rnext=rr*sign
	else
		rnext=rrold*sign
	end if
!
	return
	end

!*********************************************************

	function rnexta(r,mode)

! finds closest value to r (absolute, i.e., respects negative values)

	implicit none

	double precision rnexta
	double precision r
	integer mode

	double precision rabs
	integer m

	rabs = abs(r)
	m = mode
	if( r .lt. 0. ) m = -mode

	rabs = rnext(rabs,m)

	if( r .lt. 0. ) rabs = -rabs

	rnexta = rabs

	end

!*****************************************************************

	function rnextsub(r)

! finds best subdivision for value r
!
!		... |1| : 1. 2. 2.5 5. 8.
!		... |2| : 1. 2. 5.
!		... |3| : 1. 2. 3. 4. 5. 8.
!		... |4| : 1. 2. 3. 4. 5. 6. 7. 8. 9.

	implicit none

	double precision rnextsub
	double precision r

	integer i
	double precision eps,fact,rr,rsub
	double precision rdata(9)
	save rdata
	data rdata /0.25,0.5,1.,1.,1.,2.,1.,2.,3./

	eps = 1.e-5

	fact = 1.
	rr = r

	if( rr > 1 ) then
	  do while( rr/10. > 1 )
	    fact = fact*10.
	    rr = rr / 10.
	  end do
	else
	  do while( rr < 1 )
	    fact = fact/10.
	    rr = rr * 10.
	  end do
	end if

	if( abs(rr-2.5) < eps ) then
	  rsub = 0.5
	else
	  i = nint(rr)
	  if( i < 1 .or. i > 9 ) goto 99
	  rsub = rdata(i)
	end if

	rnextsub = rsub * fact

	return
   99	continue
	write(6,*) r,rr,fact,i
	stop 'error stop rnextsub: internal error'
	end

!*****************************************************************
!
	function istell(r)
!
! computes exponent of r
!
! -2     0.01 <= |r| < 0.1
! -1      0.1 <= |r| < 1
!  0        1 <= |r| < 10
!  1       10 <= |r| < 100
!  2      100 <= |r| < 1000
!
	implicit none

	integer istell
	double precision r

	integer i
	double precision eps,ar

	ar = abs(r)
	eps = 0.0001*ar

	i = log10(ar+eps)

	if(ar+eps.lt.1.) i=i-1

	istell = i

	return
	end
!
!*******************************************************
!
	subroutine scplo(xmin,xmax,ymin,ymax)
!
! scales plot
!
! xmin,xmax	min/max values on x axis
! ymin,ymax	min/max values on y axis
!
! inquires on terminal if another scale is desired
!
        use convert, only: ialfa

	logical bfirst,bold,bnew,bchang,bdif
	character*80 line
        double precision ymax
        double precision ymin
        double precision xmin
        double precision xmax

	save xminh,xmaxh,yminh,ymaxh
	save bfirst
	data xminh,xmaxh,yminh,ymaxh /0.,0.,0.,0./
	data bfirst /.true./
!
	bold=.not.bfirst	!ask for old values
	bchang=.true.		!values changed
!
!	bdif = xmin.ne.xminh .or. xmax.ne.xmaxh .or.	!different values
!     +		ymin.ne.yminh .or. ymax.ne.ymaxh 	!...from last call
	bdif = .true.	!always ask for old dimensions (11.10.95)
!
	do while(bchang)
	   write(6,*) 'Dimensions of plot :'
	   write(6,*)
!
	   idum=ialfa(ymax,line(1:33),4,0)
	   write(6,1001) line(1:33)
 1001	   format(1x,a)
	   write(6,1002)
 1002	   format(17x,'|'/17x,'|')
	   idum=ialfa(xmin,line(1:12),2,0)
	   idum=ialfa(xmax,line(22:33),2,0)
	   line(13:21)='----+----'
	   write(6,1001) line(1:33)
	   write(6,1002)
	   idum=ialfa(ymin,line(1:33),4,0)
	   write(6,1001) line(1:33)
	   write(6,*)
!
	   bnew=.true.			!ask for new values
	   bchang=.false.		!values changed
!
	   if(bold.and.bdif) then
	      bold=.false.
	      ianz=iantw('You want last dimensions ?')
	      if(ianz.eq.1) then
		xmin=xminh
		xmax=xmaxh
		ymin=yminh
		ymax=ymaxh
		bnew=.false.
		bchang=.true.
	      end if
	   end if
!
	   if(bnew) then
	      ianz=iantw('You want other dimensions ?')
	      if(ianz.eq.1) then
		xmin=fdeflt(xmin,'Enter xmin')
		xmax=fdeflt(xmax,'Enter xmax')
		ymin=fdeflt(ymin,'Enter ymin')
		ymax=fdeflt(ymax,'Enter ymax')
		bchang=.true.
	      end if
	   end if
	end do
!
	xminh=xmin
	xmaxh=xmax
	yminh=ymin
	ymaxh=ymax
	bfirst=.false.
!
	idum=2*idum
!
	return
	end
!
!************************************************************
!
	subroutine swapr(a,b)
!
! swaps variables a,b
!
	tmp=a
	a=b
	b=tmp
!
	return
	end
!
!************************************************************
!
	function rlen(x1,y1,x2,y2)
!
! length of line (x1,y1) --- (x2,y2)
!
	rlen = sqrt ( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) )
!
	return
	end
!
!************************************************************

	subroutine triml(line)

! trims line (deleting leading spaces)

	implicit none

	character*(*) line

	integer il

	call trimline(line,il)

	end

!************************************************************

	subroutine trimline(line,il)

! trims line (deleting leading spaces) and gives back length

	implicit none

	character*(*) line
	integer il

	integer i,j,n,length

	length = len(line)

	do i=1,length
	  if( line(i:i) .ne. ' ' ) goto 1
	end do

    1	continue

	if( i .eq. 1 ) return		!already trimmed

	n = 0
	do j=i,length
	  n = n + 1
	  line(n:n) = line(j:j)
	end do

	line(n+1:length) = ' '
	il = n

	end

!************************************************************

        subroutine logvals(amin,amax,idiv,ntk,rval,aval)

        implicit none

        double precision amin,amax          !min/max of scale
        integer idiv            !division of log scale (see below) [1-3]
        integer ntk             !dimension of aval (in), values in aval (out)
        double precision rval(ntk)          !relative x values (return)
        double precision aval(ntk)          !log scale values (return)

        double precision a1,a2,aa1,aa2,val
        double precision a,aaux,fact,r
        double precision aamin,aamax
        integer ia1,ia2,i
        integer ndim,ip

        double precision eps
        parameter (eps=0.1)

! idiv must be in [1-3]
! idiv = 1      1 10 100
! idiv = 2      1 2 10 20 100
! idiv = 3      1 2 5 10 20 50 100

        if( idiv .lt. 1 .or. idiv .gt. 3 ) goto 99
        if( amin .le. 0. .or. amin .ge. amax ) goto 97

        ndim = ntk

        a1 = alog10(real(amin))
        a2 = alog10(real(amax))
        ia1 = a1
        ia2 = a2

        if( abs(a1-ia1) .gt. eps ) then
          aa1 = ceiling(a1)
        else
          aa1 = nint(a1)
        end if
        ia1 = nint(aa1)

        if( abs(a2-ia2) .gt. eps ) then
          aa2 = floor(a2)
        else
          aa2 = nint(a2)
        end if
        ia2 = nint(aa2)

        aamin = 10.**ia1
        aamax = 10.**ia2

        write(6,*) 'log general: ',idiv,ndim,ia2-ia1+1
        write(6,*) 'log min: ',amin,a1,aa1,ia1
        write(6,*) 'log max: ',amax,a2,aa2,ia2
        write(6,*) 'log aa: ',aamin,aamax

        ip = 0
        a = aamin
        do while( a <= aamax )
          write(6,*) ip,a
          do i=1,idiv
            fact = i
            if( i .eq. 3 ) fact = 5
            aaux = a*fact
            if( aaux >= amin .and. aaux <= amax ) then
              ip = ip + 1
              if( ip .gt. ndim ) goto 98
              aval(ip) = aaux
              r = alog10(real(aaux))
              r = (r-a1)/(a2-a1)
              rval(ip) = r
            end if
          end do
          a = 10. * a
        end do
        ntk = ip

        write(6,*) 'log scale: ',ntk
        write(6,*) (aval(i),i=1,ntk)
        write(6,*) (rval(i),i=1,ntk)

        return
   97   continue
        write(6,*) 'amin,amax: ',amin,amax
        stop 'error stop logval_adjust: amin,amax'
   98   continue
        write(6,*) 'ndim = ',ndim
        write(6,*) (aval(i),i=1,ndim)
        stop 'error stop logval_adjust: ndim'
   99   continue
        write(6,*) 'idiv = ',idiv
        stop 'error stop logval_adjust: idiv'
        end

!************************************************************

	subroutine i2s0(number,string)

! converts integer to string with blanks substituted by zeros

	implicit none

	integer number
	character*(*) string

	integer i,l
	character*10 format

	l = len(string)
	write(format,'(a,i2,a)') '(i',l,')'
	write(string,format) number

	do i=1,l
	  if( string(i:i) == ' ' ) string(i:i) = '0'
	end do

	end

!************************************************************

!***********************************************************
        end module utility
!***********************************************************
