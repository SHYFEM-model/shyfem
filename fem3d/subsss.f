c
c $Id: subsss.f,v 1.11 2009-02-04 15:26:54 georg Exp $
c
c general utility routines
c
c contents :
c
c function izahl(z,ndec)		computes ciphers of number
c function laezi(l,ioff)		computes length of number
c function iscan(line,ioff,f)		converts alphanumeric text to numbers
c function iantw(l)			gets answer from terminal
c function inquir(l,f)			gets numbers from terminal
c function getnum(prompt)		gets number from terminal
c function igetxt(prompt,text)		gets text from terminal
c function ialfa(zahl,lh,ndec,mode)	converts number into alphanumeric chars
c subroutine uplow(text,type)		translates text to upper/lower case 
c subroutine tablnc(linold,linnew)	converts tabs into blanks
c subroutine prilin(line,iunit)		writes a line without trailing blanks
c function ichanm(line)			computes length of line
c function ichafs(line)			first occurrence of non blank character 
c function icompr(line)			eliminates blank characters from line
c function ifndch(line,iocc,cha)	finds the iocc occurrence of cha 
c function ifndoc(line,string)		finds number of occurences of string
c function xzentr(x,height,f,ndec)	funtion centers number
c function ideflt(k,text)		writes default and gets value (integer)
c function fdeflt(r,text)		writes default and gets value (real)
c function iround(f)			rounds real value
c function itypch(char)			returns typ of character
c subroutine tabula(tab)		creates tab character
c function rnext(r,mode)		finds closest value to r
c function rnexta(r,mode)		finds closest value to r (absolute)
c function istell(r)			computes exponent of r
c subroutine scplo(xmin,xmax,ymin,ymax)	scales plot
c subroutine swapr(a,b)			swaps values
c function rlen(x1,y1,x2,y2)		length of line
c
c uplow,itypch,tabula depend on the ASCII character set
c
c revision log :
c
c 01.06.1998	ggu	new routines iscan... into own file (subscn.f)
c 17.06.1998	ggu	old routines with new name (iscan0,ialfa0)
c 12.02.1999	ggu	new routine triml
c 26.01.2009	ggu	minor changes to avoid compiler warnings
c 16.02.2011	ggu	new routine trimline()
c 30.05.2014	ggu	new routine rnextsub()
c 06.05.2015	ggu	new routine logvals() for logarithmic values
c
c***********************************************************
c
	function izahl(z,ndec)
c
c computes ciphers of number
c
c z		number for which ciphers have to be determined
c ndec		positions after the decimal point
c		(-1 : no decimal point)
c izahl		ciphers of number z with minus sign and
c		...decimal point if necessary (return value)
c
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
c
	iz=abs(zh)+0.01
c
	if(iz.eq.0) then
		istel=1
	else
		istel=0
	end if
c
	do while (iz.gt.0)
	iz=iz/10
	istel=istel+1
	end do
c
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
c
	izahl=istel
c
	return
	end
c
c**********************************
c
	function laezi(l,ioff)
c
c computes length of number
c
c l		text where the value is stored
c ioff		offset in l from where on the value is stored
c laezi		length of number (return value)
c
	character*(*) l
c
	character*1 blank,comma,tab,cha
	data blank,comma /' ',','/
c
	call tabula(tab)
c
	do j=ioff,len(l)
	cha=l(j:j)
	if(cha.eq.blank.or.cha.eq.comma.or.cha.eq.tab) goto 1
	end do
c
    1	continue
	laezi=j-ioff
c
	return
	end
c
c**************************************
c
	function iscan0(line,ioff,f)
c
c converts alphanumeric text to numbers
c
c line		text to be translated
c ioff		offset in text from where on values have to be read
c f		array with translated values (return value)
c iscan		total number of values read (return value)
c		-1 : error reading line
c
c as separators  blank, comma and tab can be used
c inbetween two commas, a value of 0. is assumed
c
	character*(*) line
c
	character*10 lm
	character*1 lh,blank,tab,comma,plus,minus,dot,eee
	real f(1)
	data lm /'1234567890'/
	data blank,comma,plus,minus,dot,eee /' ',',','+','-','.','e'/
c
	call tabula(tab)
c
	ianz=0		!total number read
	imod=0		!1:in number
	ikom=1		!1:comma read
	ipun=0		!1:decimal point read
	iexp=0		!1:exponential part
	izei=1		!sign of number
	izeiex=0	!sign of exponential
c
	ff=1.
	ffac=1.
	kexp=0
c
	length=len(line)
c
	do j=ioff,length
c
	lh=line(j:j)
	call uplow(lh,'low')
c
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
c
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
c
		fh=0.				!find cipher
		do i=1,10
		if(lm(i:i).eq.lh) fh=i
		end do
		if(fh.eq.0.) goto 99
		if(fh.eq.10.) fh=0.
c
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
c
    1		continue
	end if
c
	end do
c
	if(imod.eq.1) then
		f(ianz)=ff*izei*(10.**(kexp*izeiex))
	end if
c
	iscan0=ianz
c
	return
c
   99	continue
	iscan0=-1
c
	return
	end
c
c******************************************
c
	function iantw(l)
c
c gets answer from terminal
c
c l		text written to terminal
c iantw		1 : true
c		0 : false
c
c the following answers can be given as true or false :
c
c		y j s t 1	==>	true
c		n f 0 'blank'	==>	false
c
	character*(*) l
c
	character l1*1
	data net,nat /5,6/
c
    1	continue
c
	write(nat,*) l
	read(net,1000) l1
 1000	format(a)
c
	call uplow(l1,'low')
c
	if(l1.eq.'y'.or.l1.eq.'j'.or.l1.eq.'s'
     +			.or.l1.eq.'t'.or.l1.eq.'1') then
		iantw=1
	else if(l1.eq.'n'.or.l1.eq.'f'.or.l1.eq.'0'
     +			.or.l1.eq.' ') then
		iantw=0
	else
		write(nat,*) 'Incorrect answer. Try again.'
		goto 1
	end if
c
	return
	end

c********************************************

	function inquire_numbers(l,f)

c gets numbers from terminal (see also iscan)
c
c do not use -> no way to know if we are out of bounds
c
c l		text written to terminal
c f		array in which the values are stored (return value)
c inquire_numbers	total number of values read in

	implicit none

	integer inquire_numbers
	character*(*) l
	real f(1)

	character*80 lh
	integer net,nat
	integer iscan
	data net,nat /5,6/

	lh=' '

	write(nat,1000) l
 1000	format(1x,a)

	read(net,2000) lh
 2000	format(a)

	inquire_numbers=iscan(lh,1,f)

	end

c********************************************
c
	function getnum(prompt)
c
c gets number from terminal (see also iscan)
c
c prompt	text written to terminal
c getnum	value of number 
c
	character*(*) prompt
c
	character*80 lh
	real f(100)
	data net,nat /5,6/
c
    1	continue
c
	write(nat,1000) prompt
 1000	format(1x,a)
c
	lh=' '
	read(net,2000) lh
 2000	format(a)
c
	if(iscan(lh,1,f).ne.1) then
		write(nat,1000) 'Erroneous input. Repeat.'
		goto 1
	end if
c
	getnum=f(1)
c
	return
	end
c
c********************************************
c
	function igetxt(prompt,text)
c
c gets text from terminal
c
c prompt	text written to terminal
c text		text returned
c igetxt	number of characters read
c
	character*(*) prompt,text
c
	data net,nat /5,6/
c
	write(nat,1000) prompt
 1000	format(1x,a)
c
	read(net,2000) text
 2000	format(a)
c
	igetxt=ichanm(text)
c
	return
	end
c
c********************************************
c
	function ialfa0(zahl,lh,ndec,mode)
c
c converts real number into alphanumeric characters
c fills zahl with '*', if dimension of lh is
c ...to small for zahl
c
c zahl		number to be converted
c lh		character string where zahl is written to
c		...(is filled with '*' if dimension of lh
c		...is to small)
c ndec		ciphers after decimal point (-1 for no decimal point)
c mode		-1	left justified
c		 0	centred
c		 1	right justified
c ialfa		total number of characters written to lh
c
	character*(*) lh
c
	character*10 lm
	data lm /'1234567890'/
c
	lentxt=len(lh)
c
	istell=izahl(zahl,ndec)
c
	ialfa=istell
	zahlh=zahl
c
	if(lentxt.lt.istell) then
		do i=1,lentxt
		lh(i:i)='*'
		end do
		ialfa=lentxt
		ialfa0=lentxt
		return
	end if
c
	ncont=0
	if(zahlh.lt.0) then
		ncont=1
		lh(1:1)='-'
		istell=istell-1
		zahlh=-zahlh
	end if
c
c integer und fraction
c
	izi=istell-ndec-1
	if(ndec.gt.0) then
		izf=ndec
		ifact=10**ndec
		zahlt=zahlh*ifact+.5
		izahli=zahlt/ifact
		izahlf=(zahlh-izahli)*ifact
c		ifact=10**ndec
c		izf=ndec
c		izahlt=zahlh*ifact+.5
c		izahlf=mod(izahlt,ifact)
c		izahli=izahlt/ifact
	else
		izahlf = 0
		ifact=1
		izf=0
		izahli=zahlh+.5
	end if
c
c integer
c
	do i=izi,1,-1
	izv=mod(izahli,10)
	if(izv.eq.0) izv=10
	lh(ncont+i:ncont+i)=lm(izv:izv)
	izahli=izahli/10
	end do
c
	ncont=ncont+izi
c
c fraction
c
	if(ndec.ge.0) then
		ncont=ncont+1
		lh(ncont:ncont)='.'
		do i=izf,1,-1
		izv=mod(izahlf,10)
		if(izv.eq.0) izv=10
		lh(ncont+i:ncont+i)=lm(izv:izv)
		izahlf=izahlf/10
		end do
c
		ncont=ncont+izf
	end if
c
	do i=ncont+1,lentxt
	lh(i:i)=' '
	end do
c
c justify text
c
	lenlh=lentxt
	ial=ialfa
c
	if(mode.gt.0) then
		ndif=lenlh-ial
	else if(mode.eq.0) then
		ndif=(lenlh-ial)/2
	else
		ndif=0
	end if
c
	do i=ial,1,-1
	   lh(i+ndif:i+ndif)=lh(i:i)
	end do
	if(ndif.ge.1) lh(1:ndif)=' '
	if(ial+ndif.lt.lenlh) lh(ial+ndif+1:lenlh)=' '
c
	ialfa0 = ialfa
c
	return
	end
c
c***************************************************
c
	subroutine uplow(text,type)
c
c translates text to upper/lower case characters
c ...depending on type
c
c works only for ASCII set of characters
c
c text		text to be transformed
c type		'up'  :  transform to uppercase
c		'low' :  transform to lowercase
c
	character*(*) text,type
c
	logical bupper
	character*1 char
c
	iascii=ichar(type(1:1))
c
	if(iascii.eq.85.or.iascii.eq.117) then
		bupper=.true.
	else if(iascii.eq.76.or.iascii.eq.108) then
		bupper=.false.
	else
		return
	end if
c
	do i=1,len(text)
c
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
c
	end do
c
	return
	end
c
c***********************************************************
c
	subroutine tablnc(linold,linnew)
c
c converts tabs into blanks
c
c linold	old line with tabs
c linnew	new converted line
c
c itab		number of blanks equals one tab (for fortran minimum 6)
c
	character*(*) linold,linnew
c
	character*1 tab
	data itab   / 8 /
c
	call tabula(tab)
	lengo=len(linold)
	lengn=len(linnew)
c
	iold=1
	inew=1
c
	do while(iold.le.lengo.and.inew.le.lengn)
c
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
c
	end do
c
	return
	end
c
c*******************************************************
c
	subroutine prilin(line,iunit)
c
c writes a line of text without trailing blanks
c
c line		line of text
c iunit		unit number of file
c
	character*(*) line
c
	iend=ichanm(line)
	write(iunit,'(a)') line(1:iend)
c
	return
	end
c
c*********************************************
c
	function ichanm(line)
c
c computes length of line without trailing blanks
c
c line		line of text
c ichanm	length of line (return value)
c		... 0 : line is all blank
c
	character*(*) line
c
	character*1 blank,tab
	data blank /' '/
c
	call tabula(tab)
c
	ndim=len(line)
c
	do i=ndim,1,-1
	if(line(i:i).ne.blank.and.line(i:i).ne.tab) goto 1
	end do
c
    1	continue
	ichanm=i
c
	return
	end
c
c*********************************************
c
	function ichafs(line)
c
c computes first occurrence of a non-blank character in line
c
c line		line of text
c ichafs	position of first non-blank character (return value)
c		... 0 : no non-blank character found
c
	character*(*) line
c
	character*1 blank,tab
	data blank /' '/
c
	call tabula(tab)
c
	ndim=len(line)
c
	do i=1,ndim
	if(line(i:i).ne.blank.and.line(i:i).ne.tab) goto 1
	end do
c
	i=0
c
    1	continue
	ichafs=i
c
	return
	end
c
c*********************************************
c
	function icompr(line)
c
c eliminates all blank characters from a character string
c
c line		line of text
c icompr	length of compressed character string (return value)
c
	character*(*) line
c
	character*1 blank,tab
	data blank /' '/
c
	call tabula(tab)
c
	ndim=len(line)
c
	ih=0
	do i=1,ndim
	if(line(i:i).ne.blank.and.line(i:i).ne.tab) then
		ih=ih+1
		line(ih:ih)=line(i:i)
	end if
	end do
c
	do i=ih+1,ndim
	line(i:i)=blank
	end do
c
	icompr=ih
c
	return
	end
c
c***************************************************
c
	function ifndch(line,iocc,cha)
c
c finds the iocc occurrence of cha in line from offset ioff on
c
c line		line of text
c iocc		occurrence of character
c		... >0 : iocc occurrence is found
c		...  0 : last occurrence is found
c		... -1 : last but one occurrence is found...
c cha		character(s) to be found
c ifndch	position of cha in line (return value)
c		...  0 : cha not found
c		... -1 : error in ioff or ndim
c
	character*(*) line,cha
c
	n=0
c
	ncha=len(cha)-1
	ndim=len(line)
c
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
c
	ifndch=0
c
	return
	end
c
c*************************************************************
c
	function ifndoc(line,string)
c
c finds number of occurences of string in line
c
	character*(*) line,string
c
	ilen=len(line)
	istr=len(string)-1
c
	n=0
	iend=ilen-istr
	do i=1,iend
	if(line(i:i+istr).eq.string(1:1+istr)) n=n+1
	end do
c
	ifndoc=n
c
	return
	end
c
c*************************************************************
c
	function xzentr(x,height,f,ndec)
c
c funtion centers number
c
c x		coordinate around which number has to be centered
c height	height of number (also width)
c f		number to be centered (real)
c ndec		cifres after decimal point (-1 : no decimal point)
c xzentr	starting coordinate of centered number
c
	istell=izahl(f,ndec)
c
	xzentr=x-istell*height*0.5
c
	return
	end
c
c************************************************************
c
	function ideflt(k,text)
c
c writes default and gets value (integer)
c
c k		default value (integer)
c text		text written to terminal
c ideflt	integer value chosen
c
	character*(*) text
c
	character line*80
	real f(10)
c
    1	continue
c
	write(6,*) text,' (default = ',k,' )'
	read(5,'(a)') line
c
	ianz=iscan(line,1,f)
c
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
c
	return
	end
c
c************************************************************
c
	function fdeflt(r,text)
c
c writes default and gets value (real)
c
c r		default value (real)
c text		text written to terminal
c fdeflt	real value chosen
c
	character*(*) text
c
	character line*80
	real f(10)
c
    1	continue
c
	write(6,*) text,' (default = ',r,' )'
	read(5,'(a)') line
c
	ianz=iscan(line,1,f)
c
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
c
	return
	end
c
c***********************************************************
c
	function iround(f)
c
c rounds real value
c
c f		real value
c iround	rounded integer value
c
	if(f.ge.0) then
		iround=f+0.5
	else
		iround=f-0.5
	end if
c
	return
	end
c
c****************************************************
c
	function itypch(char)
c
c returns typ of character
c
c works only for ASCII set of characters
c
c char		character (*1) to be looked at
c itypch	typ of character
c		1 = numeric
c		2 = letter
c		3 = special character 	=+-*/(),.'"$_!:<>%&	!'
c		4 = non-fortran character
c
	character*1 char
c
	iascii=ichar(char)
c
c number
c
	if(iascii.ge.48.and.iascii.le.57) then
		itypch=1
		return
	end if
c
c letter
c
	if(iascii.ge.65.and.iascii.le.90) then
		itypch=2
		return
	end if
c
	if(iascii.ge.97.and.iascii.le.122) then
		itypch=2
		return
	end if
c
c special character
c
	if(iascii.ge.32.and.iascii.le.47.and.iascii.ne.35) then
		itypch=3
		return
	end if
c
	if(iascii.ge.58.and.iascii.le.62.and.iascii.ne.59) then
		itypch=3
		return
	end if
c
	if(iascii.eq.9.or.iascii.eq.95) then
		itypch=3
		return
	end if
c
c no fortran character
c
	itypch=4
c
	return
	end
c
c****************************************************************
c
	subroutine tabula(tab)
c
c creates tab character
c
c works only for ASCII set of characters
c
c tab		tab character (return value)
c
	character*1 tab,char
c
	tab=char(9)
c
	return
	end

c*********************************************************

	function rnext(r,mode)
c
c finds closest value to r
c
c r		to r the closest value has to be found
c mode		mode of operation
c		...  0  : r is returned (no change)
c		... >0  : the higer value is found (further from 0)
c		... <0  : the lower value is found (closer to 0)
c		... |1| : 1. 2. 2.5 5. 8.
c		... |2| : 1. 2. 5.
c		... |3| : 1. 2. 3. 4. 5. 8.
c		... |4| : 1. 2. 3. 4. 5. 6. 7. 8. 9.
c rnext		closest value found to r
c
c val		matrix containing the closest values to be used
c		...for each mode (in ascending order)
c nval		number of values to be used for each mode
c
c if r is too small, 0 is returned
c for negative r the lower value is the value closer to zero
c
	parameter (nmodim=4,nvadim=9)
	logical bhigh
	real val(nvadim,nmodim)
	integer nval(nmodim)
c
	data val / 1. , 2. , 2.5 , 5. , 8. ,           0.,0.,0.,0.
     +		,  1. , 2. , 5. ,                      0.,0.,0.,0.,0.,0.
     +		,  1. , 2. , 3. , 4. , 5. , 8. ,       0.,0.,0.
     +		,  1. , 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9.
     +		 /
	data nval / 5 , 3 , 6 , 9 /
c
	if(mode.lt.0) then	!upper or lower value
		bhigh=.false.
		m=-mode
	else
		bhigh=.true.
		m=mode
	end if
c
	if(m.gt.nmodim.or.m.eq.0) then	!wrong mode, return r
		rnext=r
		return
	end if
c
	if(abs(r).lt.1.e-5) then	!r too small, return 0
		rnext=0.
		return
	end if
c
	if(r.lt.0.) then	!r negative --> make positive
		rin=-r
		sign=-1.
	else
		rin=r
		sign=1.
	end if
c
c find next value
c
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
c
	if(bhigh.or.rr.eq.rin) then	!return next value
		rnext=rr*sign
	else
		rnext=rrold*sign
	end if
c
	return
	end

c*********************************************************

	function rnexta(r,mode)

c finds closest value to r (absolute, i.e., respects negative values)

	implicit none

	real rnexta
	real r
	integer mode

	real rabs,rnext
	integer m

	rabs = abs(r)
	m = mode
	if( r .lt. 0. ) m = -mode

	rabs = rnext(rabs,m)

	if( r .lt. 0. ) rabs = -rabs

	rnexta = rabs

	end

c*****************************************************************

	function rnextsub(r)

c finds best subdivision for value r
c
c		... |1| : 1. 2. 2.5 5. 8.
c		... |2| : 1. 2. 5.
c		... |3| : 1. 2. 3. 4. 5. 8.
c		... |4| : 1. 2. 3. 4. 5. 6. 7. 8. 9.

	implicit none

	real rnextsub
	real r

	integer i
	real eps,fact,rr,rsub
	real rdata(9)
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

c*****************************************************************
c
	function istell(r)
c
c computes exponent of r
c
c -2     0.01 <= |r| < 0.1
c -1      0.1 <= |r| < 1
c  0        1 <= |r| < 10
c  1       10 <= |r| < 100
c  2      100 <= |r| < 1000
c
	implicit none

	integer istell
	real r

	integer i
	real eps,ar

	ar = abs(r)
	eps = 0.0001*ar

	i = log10(ar+eps)

	if(ar+eps.lt.1.) i=i-1

	istell = i

	return
	end
c
c*******************************************************
c
	subroutine scplo(xmin,xmax,ymin,ymax)
c
c scales plot
c
c xmin,xmax	min/max values on x axis
c ymin,ymax	min/max values on y axis
c
c inquires on terminal if another scale is desired
c
	logical bfirst,bold,bnew,bchang,bdif
	character*80 line
	save xminh,xmaxh,yminh,ymaxh
	save bfirst
	data xminh,xmaxh,yminh,ymaxh /0.,0.,0.,0./
	data bfirst /.true./
c
	bold=.not.bfirst	!ask for old values
	bchang=.true.		!values changed
c
c	bdif = xmin.ne.xminh .or. xmax.ne.xmaxh .or.	!different values
c     +		ymin.ne.yminh .or. ymax.ne.ymaxh 	!...from last call
	bdif = .true.	!always ask for old dimensions (11.10.95)
c
	do while(bchang)
	   write(6,*) 'Dimensions of plot :'
	   write(6,*)
c
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
c
	   bnew=.true.			!ask for new values
	   bchang=.false.		!values changed
c
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
c
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
c
	xminh=xmin
	xmaxh=xmax
	yminh=ymin
	ymaxh=ymax
	bfirst=.false.
c
	idum=2*idum
c
	return
	end
c
c************************************************************
c
	subroutine swapr(a,b)
c
c swaps variables a,b
c
	tmp=a
	a=b
	b=tmp
c
	return
	end
c
c************************************************************
c
	function rlen(x1,y1,x2,y2)
c
c length of line (x1,y1) --- (x2,y2)
c
	rlen = sqrt ( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) )
c
	return
	end
c
c************************************************************

	subroutine triml(line)

c trims line (deleting leading spaces)

	implicit none

	character*(*) line

	integer il

	call trimline(line,il)

	end

c************************************************************

	subroutine trimline(line,il)

c trims line (deleting leading spaces) and gives back length

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

c************************************************************

        subroutine logvals(amin,amax,idiv,ntk,rval,aval)

        implicit none

        real amin,amax          !min/max of scale
        integer idiv            !division of log scale (see below) [1-3]
        integer ntk             !dimension of aval (in), values in aval (out)
        real rval(ntk)          !relative x values (return)
        real aval(ntk)          !log scale values (return)

        real a1,a2,aa1,aa2,val
        real a,aaux,fact,r
        real aamin,aamax
        integer ia1,ia2,i
        integer ndim,ip

        real eps
        parameter (eps=0.1)

c idiv must be in [1-3]
c idiv = 1      1 10 100
c idiv = 2      1 2 10 20 100
c idiv = 3      1 2 5 10 20 50 100

        if( idiv .lt. 1 .or. idiv .gt. 3 ) goto 99
        if( amin .le. 0. .or. amin .ge. amax ) goto 97

        ndim = ntk

        a1 = alog10(amin)
        a2 = alog10(amax)
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
              r = alog10(aaux)
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

c************************************************************

