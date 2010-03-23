c
c $Id: subnls.f,v 1.19 2006/09/22 12:02:35 georg Exp $
c
c namelist read routines
c
c contents :
c
c subroutine nrdini(iunit)			initializes unit number
c
c subroutine setsec(bset,name,num)		memorizes section name
c
c function nrdsec(name,num)			finds next section
c
c subroutine nrdskp				skips over data in section
c function nrdlin(line)				reads next line in section
c
c function nrdnxt(name,value,text)		returns next item in section
c function nrdval(value)			returns next value
c function nrdpar(sect,name,value,text)		reads & inserts next parameter
c subroutine nrdins(sect)			reads & inserts parameters
c
c function nrdveci(ivect,ndim)			reads integer vector in section
c function nrdvecr(rvect,ndim)			reads real vector in section
c
c function nrdnls(name,value,text,line,ioff)	name list read
c function nrdvar(name,line,ioff)		reads variable name
c function nrdtxt(text,line,ioff)		reads character string
c function nrdnum(value,line,ioff)		converts text to number
c
c subroutine nrdtst				subroutine to test nrd...
c subroutine errpnt(iunit,line,ioff)		writes a pointer to character
c
c revision log :
c
c 01.06.1997	ggu	restructured (localizing nrd functions)
c 17.06.1997	ggu	nlsh, nlsa out of file
c 07.05.1998	ggu	nrdveci, nrdvecr return -1 on error
c 06.08.1998	ggu	nrdtst changed due to compiler warning of AIX
c 01.02.2000	ggu	new routine setsec
c 06.11.2000	ggu	bug found in reading sections ($NRDLIN)
c 05.08.2003	ggu	accept string enclosed in '..' and ".."
c 11.03.2005	ggu	work in double precision in nrdnum()
c 07.11.2005	ggu	better debugging output in nrdpar
c
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine setsec(name,num)

c memorizes section name

	implicit none

	character*(*) name	!section name
	integer num		!number of section

	character*80 sname	!name of section
	integer snum		!number of section
	logical sread		!section has been read
	common /secsec/ snum,sread,sname
	save /secsec/

	sname = name
	snum = num
	sread = .false.

	end

c******************************************

	subroutine getsec(name,num)

c gets section name

	implicit none

	character*(*) name	!section name
	integer num		!number of section

	character*80 sname	!name of section
	integer snum		!number of section
	logical sread		!section has been read
	common /secsec/ snum,sread,sname
	save /secsec/

	name = sname
	num = snum 

	end

c******************************************

	function handlesec(name)

c checks if can handle section name

	implicit none

	logical handlesec	!true if can handle section
	character*(*) name	!section name

	character*80 sname	!name of section
	integer snum		!number of section
	logical sread		!section has been read
	common /secsec/ snum,sread,sname
	save /secsec/

	if( name .eq. sname ) then
	  sread = .true.
	  handlesec = .true.
	else
	  handlesec = .false.
	end if

	end

c******************************************

	function hasreadsec()

c actual section has been read ?

	implicit none

	logical hasreadsec	!true if section has been read

	character*80 sname	!name of section
	integer snum		!number of section
	logical sread		!section has been read
	common /secsec/ snum,sread,sname
	save /secsec/

	hasreadsec = sread

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine nrdini(iunit)

c initializes unit number for name list read

	implicit none
	integer iunit
	integer unit
	common /nrdcom/ unit

	unit = iunit

	end

c******************************************************************

	function nrdsec(name,num)

c finds next section
c
c	finds next section and returns name and number of section
c	if section is not numbered, num = 0
c	if section is found, nrdsec = 1, else nrdsec = 0
c
c revision history :
c
c 08.09.1997	ggu	!$IREAD - bug in internal read -> use iscan

	implicit none

	integer nrdsec,num
	character*(*) name
	integer unit
	common /nrdcom/ unit

	character*80 linaux,line
	character*1 c
	real f(5)
	integer i,iend,ioff,ios
	integer nrdvar,itypch,ichafs,iscan

	nrdsec = 0
	num = 0
	name = ' '

c read until '&' or '$' as first non white space char of line found

    1	continue
	  read(unit,'(a)',iostat=ios) linaux
	  if( ios .gt. 0 ) goto 98			!read error
	  if( ios .lt. 0 ) return			!end of file

	  call tablnc(linaux,line)
	  ioff=ichafs(line)
	  if(ioff.eq.0) goto 1				!nothing on line
	  c = line(ioff:ioff)
	  if( c .ne. '&' .and. c .ne. '$' ) goto 1	!not section

c start of section found -> find name and number

	ioff=ioff+1
	iend=nrdvar(name,line,ioff)
	if( iend .le. 0 ) goto 99
	call uplow(name,'low')

c name found -> look if there is a number at end of name

	i = iend
	do while( i .gt. 0 .and. itypch(name(i:i)) .eq. 1 )	!number
	  i = i - 1
	end do
	i = i + 1

c if there is a number, strip it and put it in num

	if( i .le. iend ) then			!number to read
	  ioff = iscan(name(i:iend),1,f)  	!$IREAD
	  if( ioff .ne. 1 ) goto 97
	  num = nint(f(1))
	  name(i:iend) = ' '
	end if

c ok, new section found

	nrdsec = 1

	return
   97	continue
	write(6,*) 'error reading section number in following line'
	write(6,*) line
	stop 'error stop : nrdsec'
   98	continue
	write(6,*) 'error reading from unit ',unit
	stop 'error stop : nrdsec'
   99	continue
	write(6,*) 'error in following line'
	write(6,*) line
	stop 'error stop : nrdsec'
	end

c******************************************************************

	subroutine nrdskp

c skips over data in section

	implicit none

	character*80 line
	integer nrdlin

	do while( nrdlin(line) .ne. 0 )
	end do

	end 

c******************************************************************

	function nrdlin(line)

c reads next line in section

	implicit none
	integer nrdlin		!1 ok,  0 end of section
	character*(*) line
	integer unit
	common /nrdcom/ unit

	character*80 linaux
	character*1 c
	logical bdebug
	integer ios,ioff
	integer ichafs
	integer num

	bdebug = .false.
	bdebug = .true.

	nrdlin = 0

c loop until non empty line found

    1	continue
	  read(unit,'(a)',iostat=ios) linaux
	  if( ios .gt. 0 ) goto 98			!read error
	  if( ios .lt. 0 ) goto 97			!EOF

	  if( bdebug ) write(6,'(a)') linaux(1:79)

	  call tablnc(linaux,line)
	  ioff=ichafs(line)
	  if(ioff.eq.0) goto 1				!nothing on line

c see if end of section found - no new section may be opened

	c = line(ioff:ioff)
	if( c .eq. '&' .or. c .eq. '$' ) then
	  if( line(ioff+1:) .ne. 'end' ) goto 99
	  return					!end of section
	end if

c ok, new line read

	nrdlin = 1

	return
   97	continue
	call getsec(line,num)
	write(6,*) 'EOF found - section not ended: ',line
	stop 'error stop : nrdlin'
   98	continue
	write(6,*) 'error reading from unit ',unit,' iostat = ',ios
	stop 'error stop : nrdlin'
   99	continue
	write(6,*) 'error in following line'
	write(6,*) line
	call getsec(line,num)
	write(6,*) 'end of section expected: ',line
	stop 'error stop : nrdlin'
	end

c******************************************************************

	function nrdnxt(name,value,text)

c returns next item in current section

	implicit none

	integer nrdnxt		!return value -> type of item (>0) or 0 for
				!...end of section
	character*(*) name	!name of item read
	real value		!value of item if numeric
	character*(*) text	!text of item if string

	character*80 line
	logical bnew
	integer ioff,iweich
	integer nrdlin,nrdnls

	save bnew,line,ioff
	data bnew /.true./

	if( bnew ) then
		line = ' '
		ioff = 1
		bnew = .false.
	end if

c next section is a bug - $NRDLIN -> call to nrdlin done even if first false
c
c	iweich=nrdnls(name,value,text,line,ioff)
c	do while( iweich .eq. 0 .and. nrdlin(line) .ne. 0 )
c	    ioff = 1
c	    iweich=nrdnls(name,value,text,line,ioff)
c	end do

    1	continue
	  iweich=nrdnls(name,value,text,line,ioff)	!read next item
	  if( iweich .ne. 0 ) goto 2			!ok or error
	  if( nrdlin(line) .eq. 0 ) goto 2   !read new line -> end of section
	  ioff = 1
	  goto 1
    2	continue

	if( iweich .eq. 0 ) then	!end of section
		bnew = .true.		!prepare for next section
		name = ' '
	else if( iweich .lt. 0 ) then	!error
		write(6,*) 'error reading line :'
		write(6,*) line
		stop 'error stop nrdnxt'
	end if

	nrdnxt = iweich

	end

c******************************************************************

	function nrdval(value)

c returns next value in a value only section

	implicit none

	integer nrdval		!return value -> 1 if value read or 0 for
				!...end of section
	real value		!value of item

	character*80 line
	logical bnew
	integer ioff,istell
	integer nrdlin,nrdnum

	save bnew,line,ioff
	data bnew /.true./

	if( bnew ) then
		line = ' '
		ioff = 1
		bnew = .false.
	end if

c next section is a bug - $NRDLIN -> call to nrdlin done even if first false
c
c	istell=nrdnum(value,line,ioff)
c	do while( istell .eq. 0 .and. nrdlin(line) .ne. 0 )
c	    ioff = 1
c	    istell=nrdnum(value,line,ioff)
c	end do

    1	continue
	  istell=nrdnum(value,line,ioff)		!read next number
	  if( istell .ne. 0 ) goto 2			!ok or error
	  if( nrdlin(line) .eq. 0 ) goto 2   !read new line -> end of section
	  ioff = 1
	  goto 1
    2	continue

	if( istell .eq. 0 ) then	!end of section
		bnew = .true.		!prepare for next section
	else if( istell .lt. 0 ) then	!error
		write(6,*) 'error reading line :'
		write(6,*) line
		stop 'error stop nrdval'
	end if

	nrdval = istell

	end

c******************************************************************

	function nrdpar(sect,name,value,text)

c reads next parameter in section and inserts value

	implicit none

	integer nrdpar
	character*(*) sect,name,text
	real value

	integer iweich

	integer nrdnxt,itspar,iscpar,itsfnm,iscfnm

	iweich = nrdnxt(name,value,text)
	call uplow(name,'low')
	if( iweich .eq. 1 .or. iweich .eq. 2 ) then
		text = ' '
		if( itspar(name) .eq. 0 ) goto 93
		if( iscpar(name,sect) .eq. 0 ) goto 94
		call putpar(name,value)
	else if ( iweich .eq. 3 .or. iweich .eq. 4 ) then
		value = 0.
		if( itsfnm(name) .eq. 0 ) goto 93
		if( iscfnm(name,sect) .eq. 0 ) goto 94
		call putfnm(name,text)
	else if ( iweich .lt. 0 ) then
		goto 98
	end if

	nrdpar = iweich

	return
   93	continue
	write(6,*) 'no parameter with this name:'
	write(6,*) name
	stop 'error stop : nrdpar'
   94	continue
	write(6,*) 'parameter is in wrong section:'
	write(6,*) 'parameter type:  ',iweich
	write(6,*) 'parameter name:    ',name
	write(6,*) 'section: ',sect
	write(6,*) 'text:    ',text
	call get_sect_of(name,sect)
	write(6,*) 'section found: ',sect
	if( iweich .ge. 3 ) call prifnm
	if( iweich .le. 2 ) call pripar
	stop 'error stop : nrdpar'
   98	continue
	stop 'error stop : nrdpar (internal error (1))'
	end

c******************************************************************

	subroutine nrdins(sect)

c reads parameter section and inserts values automatically
c
c does not handle vectors

	implicit none

	character*(*) sect

	character*80 name,text
	real value
	integer iweich

	integer nrdpar

	iweich = 1
	do while( iweich .gt. 0 )
		iweich = nrdpar(sect,name,value,text)
		if( iweich .eq. 2 .or. iweich .eq. 4 ) goto 99
	end do

	return
   99	continue
	write(6,*) 'cannot read vector in nrdins for variable:'
	write(6,*) name
	write(6,*) sect
	stop 'error stop : nrdins'
	end

c******************************************************************

	function nrdveci(ivect,ndim)

c reads integer vector in section
c
c returns total number of values read
c returns -1 in case of dimension error
c returns -2 in case of read error

	implicit none

	integer nrdveci		!total number of elements read
	integer ndim		!dimension of vector
	integer ivect(ndim)	!vector

	character*80 line
	integer n,ioff,ianz
	integer nrdlin,nrdnum,ichafs
	real f

	n = 0

	do while( nrdlin(line) .eq. 1 )
		ioff=ichafs(line)
		ianz=nrdnum(f,line,ioff)
		do while(ianz.gt.0)
		   n=n+1
		   if(n.gt.ndim) goto 99
		   ivect(n)=nint(f)
		   ianz=nrdnum(f,line,ioff)
		end do
		if(ianz.lt.0) goto 98
	end do

	nrdveci = n

	return
   98	continue
	write(6,*) 'read error in following line'
	write(6,*) line
	nrdveci = -2
	return
c	stop 'error stop : nrdveci'
   99	continue
c	write(6,*) 'dimension error : ',ndim
	nrdveci = -1
	return
c	stop 'error stop : nrdveci'
	end

c******************************************************************

	function nrdvecr(rvect,ndim)

c reads real vector in section
c
c returns total number of values read
c returns -1 in case of dimension error
c returns -2 in case of read error

	implicit none

	integer nrdvecr		!total number of elements read
	integer ndim		!dimension of vector
	real rvect(ndim)	!vector

	character*80 line
	integer n,ioff,ianz
	integer nrdlin,nrdnum,ichafs
	real f

	n = 0

	do while( nrdlin(line) .eq. 1 )
		ioff=ichafs(line)
		ianz=nrdnum(f,line,ioff)
		do while(ianz.gt.0)
		   n=n+1
		   if(n.gt.ndim) goto 99
		   rvect(n)=f
		   ianz=nrdnum(f,line,ioff)
		end do
		if(ianz.lt.0) goto 98
	end do

	nrdvecr = n

	return
   98	continue
	write(6,*) 'read error in following line'
	write(6,*) line
	nrdvecr = -2
	return
c	stop 'error stop : nrdvecr'
   99	continue
c	write(6,*) 'dimension error : ',ndim
	nrdvecr = -1
	return
c	stop 'error stop : nrdvecr'
	end

c******************************************************************

	function nrdnls(name,value,text,line,ioff)
c
c name list read
c
c on every call routine reads one variable or returns 0
c ...if nothing has been read
c as seperators blank, tab and comma can be used
c
c reads lines of the form :
c		name = value
c		name = 'text'
c		name = value1,value2,...,valuen
c more than one variable in one line is permitted :
c		name1 = value , name2 = 'text' , ...
c vectors can be read from more than one line :
c		name = value1  value2
c		   value3 , ... , valuen
c
c name		name of variable (return value)
c		...for vector read name is not changed
c value		value of number (return value)
c text		text of character string (return value)
c line		text line from where names, values
c		...or text are to be read
c ioff		offset from where on line is scanned
c		...on return ioff points to the first character after
c		...the last character read (input and return value)
c nrdnls	type of variable read :
c			-1 : error
c			 0 : end of line, nothing read
c			 1 : number variable with name
c			 2 : number variable without name
c			 3 : character variable with name
c			 4 : character variable without name
c
	character*(*) name,text,line
	character*1 ll
	character*1 caux
	logical brdvar,brdequ,brdnum
c
	namlen=len(line)
	if( ioff .lt. 1 ) ioff = 1	!$$ggu 12.06.1997
c
	logtyp=0	!type of variable read (local nrdnls)
	brdvar=.true.	!variable-name to be read
	brdequ=.true.	!= sign to be read
	brdnum=.true.	!number to be read
c
c get variable name %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	i=ioff
	do while(i.le.namlen.and.brdvar)
	ll=line(i:i)
	ityp=itypch(ll)
c	write(6,*) ityp,ioff,i,ll,ichar(ll)
	if(ityp.eq.1.or.ll.eq.'+'.or.ll.eq.'-'
     +			.or.ll.eq.'.') then	!number
		logtyp=2			!...without name
		brdvar=.false.
		brdequ=.false.
	else if(ityp.eq.2) then			!letter
		logtyp=1
		ianz=nrdvar(name,line,i)
		if(ianz.eq.-1) logtyp=-1
		brdvar=.false.
	else if(line(i:i).eq.'''') then		!character string
		logtyp=4			!...without name
		brdvar=.false.
		brdequ=.false.
	else if(line(i:i).eq.' ') then		!blank
		i=i+1
	else if(line(i:i).eq.',') then		!comma
		i=i+1
	else					!not recognized
		logtyp=-1
		brdvar=.false.
	end if
	end do
c
	if(brdvar) then			!end of line
		brdequ=.false.		!...stop reading
		brdnum=.false.
	end if
c
	if(logtyp.eq.-1) brdequ=.false.
c
c get equal sign %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	do while(i.le.namlen.and.brdequ)
	if(line(i:i).eq.' ') then
		i=i+1
	else if(line(i:i).eq.'=') then
		i=i+1
		brdequ=.false.
	else
		logtyp=-1
		brdequ=.false.
	end if
	end do
c
	if(brdequ) logtyp=-1		!looking for = but not found
	if(logtyp.eq.-1) brdnum=.false.
c
c get number (or character) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	do while(i.le.namlen.and.brdnum)
	    caux = line(i:i)
	    if( caux .eq. '''' .or. caux .eq. '"' ) then   !character string
		if(logtyp.ne.4) logtyp=3		   !with name ==> 3
		ianz=nrdtxt(text,line,i)
		if(ianz.eq.-1) logtyp=-1
		brdnum=.false.
	    else if(caux.eq.' ') then
		i=i+1
	    else
		ianz=nrdnum(value,line,i)
		if(ianz.eq.-1) logtyp=-1
		brdnum=.false.
	    end if
	end do
c
	if(logtyp.eq.1.and.brdnum) then		!variable name read
		logtyp=-1			!...but no value found
	end if
c
	ioff=i
	nrdnls=logtyp
c
	return
	end
c
c********************************************************
c
	function nrdvar(name,line,ioff)
c
c reads variable name
c
c variable has to start with a letter
c ...the rest of the variable can consist of letters or numbers
c ...program stops reading when a non permitted character is found
c ...leading blank characters are discarted
c
c name		name of variable (return value)
c line		line of text from where varaible name is to be read
c ioff		offset from where on line is scanned
c		...on return ioff points to the first character after
c		...the variable name read (input and return value)
c nrdvar	number of characters in variable name (return value)
c		...(0 if no name has been read, -1 if error)
c
	character*(*) name,line
c
	character*1 cha
c
	length=len(line)
	invar=0
c
	do i=ioff,length
		cha=line(i:i)
		ityp=itypch(cha)
		if(invar.eq.0) then
			if(ityp.eq.2) then
				invar=invar+1
				name=cha
			else if(cha.eq.' ') then
c				nothing
			else
				invar=-1
				goto 1
			end if
		else
			if(ityp.eq.2.or.ityp.eq.1) then
				invar=invar+1
				name(invar:invar)=cha
			else
				goto 1
			end if
		end if
	end do
c
    1	continue
c
	nrdvar=invar
	ioff=i
c
	return
	end
c
c********************************************************
c
	function nrdtxt(text,line,ioff)
c
c reads character string
c
c reads character string enclosed in '...' or "..."
c ...leading blank characters are discarted
c
c text		character string (return value)
c line		line of text from where varaible name is to be read
c ioff		offset from where on line is scanned
c		...on return ioff points to the first character after
c		...the variable name read (input and return value)
c nrdtxt	number of characters in character string (return value)
c		...(0 if none has been read, -1 if error)
c
	character*(*) text,line
c
	character*1 cha
	character*1 delim
c
	length=len(line)
	intxt=0		!1:reading character string
	ktext=0		!number of characters in character string
c
	i=ioff
	do while(i.le.length)
		cha=line(i:i)
		if(intxt.eq.0) then
			if(cha.eq.''''.or.cha.eq.'"') then
				delim=cha
				intxt=1
				text=' '
				i=i+1
			else if(cha.eq.' ') then
				i=i+1
			else
				ktext=-1
				goto 1
			end if
		else
			if(cha.eq.delim) then
				i=i+1
				if(i.gt.length) then
					intxt=0
					goto 1
				else if(line(i:i).eq.delim) then !e.g. it''s
					ktext=ktext+1
					text(ktext:ktext)=delim
					i=i+1
				else
					intxt=0
					goto 1
				end if
			else
				ktext=ktext+1
				text(ktext:ktext)=cha
				i=i+1
			end if
		end if
	end do
c
    1	continue
c
	if(intxt.eq.1) ktext=-1
c
	nrdtxt=ktext
	ioff=i
c
	return
	end
c
c*******************************************************************
c
	function nrdnum(value,line,ioff)
c
c converts alphanumeric text to number, only one number is read
c
c as separators  blank, tab and comma can be used
c ...end of line works as comma
c ...inbetween two commas, a value of 0. is assumed
c ...in this case a length of 1 is returned in nrdnum
c
c line		text to be translated
c ioff		offset in text from where on the number has to be read
c		...on return ioff points to the first character after
c		...the number read (input and return value)
c value		translated number (return value)
c nrdnum	ciphers in number (return value)
c		 0 : nothing read, end of line
c		-1 : error reading line

	implicit none

	integer nrdnum
	real value
	character*(*) line
	integer ioff

	integer ianz,imod,ikom,ipun,iexp,izei,izeiex,ieol
	integer length,j,i
	!real ff,kexp,ffac,fh
	double precision ff,kexp,ffac,fh,val

	character*10 lm
	character*1 lh,blank,tab,comma,plus,minus,dot,eee
	save ikom,ieol
	data lm /'1234567890'/
	data blank,comma,plus,minus,dot,eee /' ',',','+','-','.','e'/
	data ikom /1/	!1:comma read
	data ieol /1/	!1:end of line
c
	call tabula(tab)
c
	ianz=0		!number consists of ianz ciphers
	imod=0		!1:in number
c	ikom=1		!1:comma read
	ipun=0		!1:decimal point read
	iexp=0		!1:exponential part
	izei=1		!sign of number
	izeiex=0	!sign of exponential
c
	ff=1.		!number
	kexp=0.		!exponential
c
	length=len(line)
c
	do j=ioff,length
c
	lh=line(j:j)
	call uplow(lh,'low')
c
	if(lh.eq.blank.or.lh.eq.tab) then		!blank
		if(imod.eq.1) then	!end of number found
			goto 2
		end if
	else if(lh.eq.comma) then	!comma
		ieol=0			!no end of line any more
		if(imod.eq.1) then	!end of number found
			goto 2
		else if(ikom.eq.1) then	!second comma found
			ff=0.		!assume 0 inbetween
			imod=1		!number read
			ikom=0		!assume no comma for next read
			ianz=1		!assume number of one cipher
			goto 2
		end if
		ikom=1			!comma read
	else if(lh.eq.eee) then		!exponential
			iexp=1
			imod=1
			ikom=0
			ieol=0
			ianz=ianz+1
	else					!number
		if(imod.eq.0) then		!start reading
			imod=1			!reading number
			ipun=0
			ikom=0
			ieol=0			!no eol any more
			ff=0.
			if(lh.eq.plus) then	! + sign
				izei=+1
				ianz=ianz+1
				goto 1
			else if(lh.eq.minus) then	! - sign
				izei=-1
				ianz=ianz+1
				goto 1
			end if
		end if
c
		if(iexp.eq.1.and.izeiex.eq.0) then	!exponential
			if(lh.eq.plus) then		! + sign
				izeiex=+1
				ianz=ianz+1
				goto 1
			else if(lh.eq.minus) then	! - sign
				izeiex=-1
				ianz=ianz+1
				goto 1
			else				!no sign, assume +
				izeiex=+1
			end if
		end if
c
		if(lh.eq.dot) then		!point
			if(ipun.eq.1) goto 99
			if(iexp.eq.1) goto 99
			ipun=1
			ffac=1.			!prepare for decimal part
			fh=0.			!aux variable
			ianz=ianz+1
			goto 1
		end if
c						!find cipher
		fh=0.
		do i=1,10
		if(lm(i:i).eq.lh) fh=i
		end do
		if(fh.eq.0.) goto 99
		if(fh.eq.10.) fh=0.
c
		if(iexp.eq.1) then		!cipher for exponential
			kexp=10*kexp+int(fh)
		else 				!cipher for number
			if(ipun.eq.0) then
				ff=10.*ff+fh
			else
				ffac=ffac/10.
				ff=ff+ffac*fh
			end if
		end if
		ianz=ianz+1
	end if
c
    1	continue
	end do
c
c exit for end of line
c
	if(ikom.eq.1.and.ieol.eq.0) then	!comma read but no
		ff=0.				!...empty line ==>
		imod=1				!...assume 0 has been read
		ianz=1
	else if(imod.eq.0) then			!nothing read
		ff=0.
	end if
c
	ikom=1			!for next call assume comma
	ieol=1			!...and end of line
c
c exit for end of number
c
    2	continue
c
	fh=1.
	do i=1,kexp
	fh=fh*10.
	end do
c
	if(izeiex.ge.0) then
		val=ff*izei*fh
	else
		val=ff*izei/fh
	end if

	value = val
	!write(6,*) 'nrdnum: ',ff,izei,fh,value,val
c
	nrdnum=ianz
	ioff=j
c
	return
c
   99	continue
	nrdnum=-1
	ioff=j
c
	return
	end

c************************************************************

	subroutine nrdtst

c subroutine to test nrd... routines
c
c to use write main as follows :
c
c--------------------------
c	call nrdtst
c	end
c--------------------------

	logical bloop
	character*80 line

	bloop = .true.

	do while( bloop )

	  line = ' '
	  write(6,*) 'Enter line (<CR> to end) :'
	  read(5,'(a)') line
	  write(6,*) line

	  if( line .eq. ' ' ) bloop = .false.

	  ianz=1
	  ioff=1
	  do while(ianz.gt.0)
	    ianz=nrdnum(f,line,ioff)
	    write(6,*) f,ioff,ianz
	  end do

	end do

	end

c*******************************************************
c
	subroutine errpnt(iunit,line,ioff)
c
c writes a pointer to character
c
c iunit		unit number pointer is written to
c ioff		pointer is written in column ioff
c
	character*(*) line
c
	character form*80, pointr*1
	data pointr /'^'/
c
	length=len(line)
c
	write(iunit,*) line
c
	if(ioff.lt.1.or.ioff.gt.length) return
c
	write(form,1000) ioff
 1000	format('(',i3,'x,a1)')
c
	write(iunit,form) pointr
c
	return
	end
