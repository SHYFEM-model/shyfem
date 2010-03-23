c
c $Id: subisd.f,v 1.1 2003/09/04 11:07:35 georg Exp $
c
c is_ routines for chars
c
c contents :
c
c revision log :
c
c 03.09.2003	ggu	new routines written
c
!****************************************************************

	function is_digit(c)

c true if c is digit

	implicit none

	logical is_digit
	character*1 c

	integer i0,i9
	parameter(i0=48,i9=57)

	integer ic

	ic = ichar(c)

	if( i0 .le. ic .and. ic .le. i9 ) then
	  is_digit = .true.
	else
	  is_digit = .false.
	end if

	end

!****************************************************************

	function is_upper(c)

c true if c is upper case letter

	implicit none

	logical is_upper
	character*1 c

	integer i0,i9
	parameter(i0=65,i9=90)

	integer ic

	ic = ichar(c)

	if( i0 .le. ic .and. ic .le. i9 ) then
	  is_upper = .true.
	else
	  is_upper = .false.
	end if

	end

!****************************************************************

	function is_lower(c)

c true if c is lower case letter

	implicit none

	logical is_lower
	character*1 c

	integer i0,i9
	parameter(i0=97,i9=122)

	integer ic

	ic = ichar(c)

	if( i0 .le. ic .and. ic .le. i9 ) then
	  is_lower = .true.
	else
	  is_lower = .false.
	end if

	end

!****************************************************************

	function is_letter(c)

c true if c is letter

	implicit none

	logical is_letter
	character*1 c
	logical is_lower,is_upper

	is_letter = is_lower(c) .or. is_upper(c)

	end

!****************************************************************

	function is_alpha(c)

c true if c is alphanumeric char

	implicit none

	logical is_alpha
	character*1 c
	logical is_letter,is_digit

	is_alpha = is_letter(c) .or. is_digit(c) .or. c .eq. '_'

	end

!****************************************************************

	subroutine test_is_digit_etc

c tests routines

	implicit none

	character*1 c
	logical is_digit
	logical is_lower
	logical is_upper
	logical is_letter
	logical is_alpha

	write(6,*) 'testing is_digit ',ichar('0'),ichar('9')

	c = '0'
	if( .not. is_digit(c) ) goto 99
	c = '5'
	if( .not. is_digit(c) ) goto 99
	c = '9'
	if( .not. is_digit(c) ) goto 99
	c = 'a'
	if( is_digit(c) ) goto 99

	write(6,*) 'testing is_lower ',ichar('a'),ichar('z')

	c = 'a'
	if( .not. is_lower(c) ) goto 99
	c = 'm'
	if( .not. is_lower(c) ) goto 99
	c = 'z'
	if( .not. is_lower(c) ) goto 99
	c = 'A'
	if( is_lower(c) ) goto 99

	write(6,*) 'testing is_upper ',ichar('A'),ichar('Z')

	c = 'A'
	if( .not. is_upper(c) ) goto 99
	c = 'M'
	if( .not. is_upper(c) ) goto 99
	c = 'Z'
	if( .not. is_upper(c) ) goto 99
	c = 'a'
	if( is_upper(c) ) goto 99

	write(6,*) 'testing is_letter '

	c = 'A'
	if( .not. is_letter(c) ) goto 99
	c = 'm'
	if( .not. is_letter(c) ) goto 99
	c = '5'
	if( is_letter(c) ) goto 99
	c = '_'
	if( is_letter(c) ) goto 99

	write(6,*) 'testing is_alpha '

	c = 'A'
	if( .not. is_alpha(c) ) goto 99
	c = 'm'
	if( .not. is_alpha(c) ) goto 99
	c = '7'
	if( .not. is_alpha(c) ) goto 99
	c = '_'
	if( .not. is_alpha(c) ) goto 99
	c = '.'
	if( is_alpha(c) ) goto 99
	c = ' '
	if( is_alpha(c) ) goto 99

	write(6,*) 'test passed'

	return
   99	continue
	write(6,*) 'error for char ',c
	stop 'error stop test_is_digit_etc'
	end

!****************************************************************

c	call test_is_digit_etc
c	end

!****************************************************************

