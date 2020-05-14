
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-1999,2009-2012,2014-2015,2014-2015  Georg Umgiesser
!    Copyright (C) 2017-2020  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c general utility routines
c
c contents :
c
c function izahl(z,ndec)		computes ciphers of number
c function laezi(l,ioff)		computes length of number
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
c 23.03.2010	ggu	changed v6.1.1
c 16.02.2011	ggu	new routine trimline()
c 30.03.2012	ggu	changed VERS_6_1_51
c 30.05.2014	ggu	new routine rnextsub()
c 12.12.2014	ggu	changed VERS_7_0_9
c 06.05.2015	ggu	new routine logvals() for logarithmic values
c 21.05.2015	ggu	changed VERS_7_1_11
c 10.07.2015	ggu	changed VERS_7_1_50
c 14.09.2015	ggu	new routine i2s0()
c 25.05.2017	ggu	changed VERS_7_5_28
c 14.11.2017	ggu	changed VERS_7_5_36
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.02.2020	ggu	rounding routines into new file subround.f
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
	call to_lower(l1)
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

	function inquire_numbers(l,f,max)

c gets numbers from terminal
c
c do not use -> no way to know if we are out of bounds
c
c l		text written to terminal
c f		array in which the values are stored (return value)
c inquire_numbers	total number of values read in

	implicit none

	integer inquire_numbers
	character*(*) l
	real f(max)
	integer max

	character*80 lh
	integer net,nat
	integer iscanf
	data net,nat /5,6/

	lh=' '

	write(nat,1000) l
 1000	format(1x,a)

	read(net,2000) lh
 2000	format(a)

	inquire_numbers=iscanf(lh,f,max)

	end

c********************************************
c
	function getnum(prompt)
c
c gets number from terminal
c
c prompt	text written to terminal
c getnum	value of number 
c
	character*(*) prompt
c
	character*80 lh
	real f(1)
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
	if(iscanf(lh,f,1).ne.1) then
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
	real f(1)
c
    1	continue
c
	write(6,*) text,' (default = ',r,' )'
	read(5,'(a)') line
c
	ianz=iscanf(line,f,1)
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

c************************************************************

	subroutine i2s0(number,string)

c converts integer to string with blanks substituted by zeros

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

c************************************************************

