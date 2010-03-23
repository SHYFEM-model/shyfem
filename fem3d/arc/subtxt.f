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

c*************************************************************
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
c		3 = special character 	=+-*/(),.'"$_!:<>%&
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
c*************************************************************
c
