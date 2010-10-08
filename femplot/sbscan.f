!
! $Id: sbscan.f,v 1.4 2008-10-10 09:47:22 georg Exp $
!
! scanning procedures : scan string for numbers and convert them
!
! contents :
!
! integer function istof(line,f,ioff)	converts string to number
! integer function iscanf(line,f,max)	converts string to numbers
! integer function icindx(string,c)	finds c in string
!

!****************************************************************

	function istof(line,f,ioff)

! converts string to number (reads exactly one number)
!
! a comma is treated like a blank (,, does not denote a value of 0)
!
! line		string to convert
! f		converted number (out)
! ioff		offset in string to start (in) 
!		position of first non blank char after number (out)
! istof		1: valid number in f    0: EOL    -1: read error

	implicit none

	integer istof
	character*(*) line
	real f
	integer ioff
c	integer, optional :: ioff

	logical berr,beol,bnumb,bexp,bsign,besign,bdot
	integer iof
	integer i,n,istart,ic
	integer iesign,kexp
	real ff,fh,fact,sign
	character*1 c

	integer icindx

	character*1 blank,tab,comma,plus,minus,dot
	character*4 expont
	character*10 number
	save blank,tab,comma,plus,minus,dot,expont,number
	data blank,tab,comma /' ','	',','/
	data plus,minus,dot /'+','-','.'/
	data expont /'eEdD'/
	data number /'1234567890'/

c	if( .not. present(ioff) ) ioff=1

	n=len(line)

! skip leading blanks

	do i=ioff,n
	  c=line(i:i)
	  if( c.ne.blank .and. c.ne.tab .and. c.ne.comma ) goto 3
	end do
    3	continue

! start reading number

	istart=i
	berr=.false.		!read error
	beol=i.gt.n		!end of line

	bnumb=.false.		!some number read
	bexp=.false.		!reading exponent
	bsign=.false.		!sign of number read
	besign=.false.		!sign of exponent read
	bdot=.false.		!decimal point read

	sign=1.0		!sign of number
	iesign=1		!sign of exponent

	ff=0.0			!number
	fact=1.			!factor for decimal part
	kexp=0			!exponent of number

	do i=istart,n

	c=line(i:i)

	if( berr .or. c.eq.comma .or. c.eq.blank .or. c.eq.tab ) goto 4

	if( c.eq.plus ) then
		if( .not.bexp .and. .not.bsign ) then
			bsign=.true.
		else if( bexp .and. .not.besign ) then
			besign=.true.
		else
			berr=.true.
		end if
	else if( c.eq.minus ) then
		if( .not.bexp .and. .not.bsign ) then
			bsign=.true.
			sign=-1.
		else if( bexp .and. .not.besign ) then
			besign=.true.
			iesign=-1
		else
			berr=.true.
		end if
	else if( c.eq.dot ) then
		if( .not.bdot .and. .not.bexp ) then
			bdot=.true.
		else
			berr=.true.
		end if
	else if( icindx(expont,c) .gt. 0 ) then
		if( .not.bexp .and. bnumb ) then
			bexp=.true.
		else
			berr=.true.
		end if
	else
		ic=icindx(number,c)
		if( ic.eq.0 ) then
			berr=.true.
		else if( ic.eq.10 ) then
			ic=0
		end if
		fh=ic
		bnumb=.true.

		if( bexp ) then
			kexp = 10*kexp + ic
		else if( bdot ) then
			fact = 0.1*fact
			ff = ff + fact*fh
		else
			ff = 10.*ff + fh
		end if
	end if

	end do

    4	continue

! this is the only exeption to the rule that there must be at least
! one digit read -> " , , " will give 0. inbetween the two commas
!
! -> not any more supported
!
!	if( istart.eq.i ) bnumb=.true.		!only comma read -> 0.

! skip trailing blanks

	istart=i
	do i=istart,n
	  c=line(i:i)
	  if( c.ne.blank .and. c.ne.tab .and. c.ne.comma ) goto 5
	end do
    5	continue

! if we have read a comma we step one char ahead -> not necessary anymore
!
!	if( i.le.n .and. c.eq.comma ) i=i+1

	if( beol ) then
		istof=0
		ioff=n+1
	else if( berr .or. .not.bnumb ) then
		istof=-1
		ioff=istart-1
	else
		istof=1
		ioff=i
		f = ff * sign * ( 10.**(kexp*iesign) )
	end if

	end 

!****************************************************************

	integer function iscanf(line,f,max)

! converts string to numbers (reads at most max numbers)
!
! line		string to convert
! f		array of converted numbers (return)
! max		how many numbers to convert at most
!		 0: count only numbers in line
!		-1: convert all numbers found (default)
! iscanf	total number of numbers converted ( >0 )
!		0: blank line   <0: read error in |iscanf|'th number

	implicit none

	character*(*) line
	real f(1)
	integer max
c	integer,optional :: max

	integer i,inumb,iret,maxf
	real ff

	integer istof

	maxf = max
c	if( present(max) ) then
c		maxf = max
c	else
c		maxf = size(f)
c	end if

	i=1
	inumb=0

	do while(.true.)
	  iret=istof(line,ff,i)
	  if( iret.le.0 ) goto 1
	  inumb=inumb+1
	  if( maxf.ne.0 ) f(inumb)=ff
	  if( inumb.eq.maxf ) goto 1
	end do
    1	continue
	if( iret.lt.0 ) inumb=-inumb-1

	iscanf=inumb

	end

!****************************************************************

	integer function icindx(string,c)

! finds c in string
!
! string	string where to look up c
! c		char to find
! icindx	position of c in string ( 0 if not found )

	implicit none

	character*(*) string
	character*1 c

	integer i,n

	n=len(string)

	do i=1,n
	  if( string(i:i).eq.c ) goto 1
	end do
    1	continue

	if( i.gt.n ) then
	  icindx=0
	else
	  icindx=i
	end if

	end

!****************************************************************

