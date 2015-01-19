c
c $Id: debug.f,v 1.10 2008-12-09 11:38:33 georg Exp $
c
c revision log :
c
c 27.03.2003    ggu     new routines for NaN check and value test
c 03.09.2003    ggu     check routines customized
c 05.12.2003    ggu     in check[12]Dr only check val if vmin!=vmax
c 06.12.2008    ggu     check for NaN changed (Portland compiler)
c 15.07.2011    ggu     new routines for checksum (CRC)
c
c notes :
c
c use of debug routines:
c
c in the calling routine you must set or unset the debug:
c
c call setdebug
c call your_routine
c call cleardebug
c
c in the routine where you want to debug:
c
c logical isdebug
c ...
c if( isdebug() ) then
c    here write your variables
c end if
c
c***************************************************************

	subroutine assigndebug(bdebug)

	implicit none

	logical bdebug

	include 'debug_aux1.h'

	debug = bdebug

	end

c***************************************************************

	subroutine setdebug

	implicit none

	include 'debug_aux1.h'

	debug = .true.

	end

c***************************************************************

	subroutine cleardebug

	implicit none

	include 'debug_aux1.h'

	debug = .false.

	end

c***************************************************************

	function isdebug()

	implicit none

	logical isdebug

	include 'debug_aux1.h'

	isdebug = debug

	end

c***************************************************************

	block data datadebug

	implicit none

	include 'debug_aux1.h'

	data debug / .false. /

	end

c***************************************************************
c***************************************************************
c***************************************************************
c***************************************************************
c***************************************************************

	function is_d_nan(val)

c tests val for NaN

	implicit none

	double precision val
	logical is_d_nan
        integer itot

        itot = 0

	if( val .gt. 0. ) itot = itot + 1
	if( val .lt. 0. ) itot = itot + 1
	if( val .eq. 0. ) itot = itot + 1

	is_d_nan = itot .ne. 1

	end

c***************************************************************

	function is_r_nan(val)

c tests val for NaN

	implicit none

	real val
	logical is_r_nan
        integer itot

        itot = 0

	if( val .gt. 0. ) itot = itot + 1
	if( val .lt. 0. ) itot = itot + 1
	if( val .eq. 0. ) itot = itot + 1

	is_r_nan = itot .ne. 1

	end

c***************************************************************

	function is_i_nan(val)

c tests val for NaN

	implicit none

	integer val
	logical is_i_nan
        integer itot

        itot = 0

	if( val .gt. 0 ) itot = itot + 1
	if( val .lt. 0 ) itot = itot + 1
	if( val .eq. 0 ) itot = itot + 1

	is_i_nan = itot .ne. 1

	end

c***************************************************************

	subroutine nantest(n,a,text)

c tests array for nan

	implicit none

	integer n
	real a(1)
	character*(*) text

	integer i
        logical is_r_nan

	do i=1,n
          if( is_r_nan( a(i) ) ) goto 99
	end do

	return
   99	continue
	write(6,*) 'NAN found while testing ',text
	write(6,*) n,i,a(i)
	stop 'error stop nantest: NAN found'
	end

c***************************************************************

	subroutine valtest(nlv,n,valmin,valmax,a,textgen,text)

c tests array for nan - use other routines below

	implicit none

	integer nlv,n
	real valmin,valmax
	real a(1)
	character*(*) textgen,text

	integer i
	integer ntot,j1,j2
        logical is_r_nan

	ntot = nlv*n

c	write(6,*) text,'   ',ntot,valmin,valmax

	do i=1,ntot
          if( is_r_nan( a(i) ) ) goto 99
	end do

	return
   99	continue
	write(6,*) textgen
	write(6,*) 'value out of range for ',text
	write(6,*) 'range: ',valmin,valmax
	write(6,*) nlv,n,i
	j1 = mod(i-1,nlv) + 1
	j2 = 1 + (i-1) / nlv
	write(6,*) j1,j2,a(i)
	stop 'error stop valtest: value out of range'
	end

c***************************************************************

	subroutine check1Dr(n,a,vmin,vmax,textgen,text)

c tests array for nan and strange values

	implicit none

	integer n
	real a(1)
	real vmin,vmax
	character*(*) textgen,text

	logical debug,bval
	integer inan,iout,i
	real val

	logical is_r_nan

        bval = vmin .lt. vmax
	debug = .true.
	debug = .false.
	inan = 0
	iout = 0

	do i=1,n
	  val = a(i)
	  if( is_r_nan(val) ) then
	    inan = inan + 1
	    if( debug ) write(6,*) 'nan ',inan,i,val
	    write(999,*) 'nan: ',inan,i,val
	  else if( bval .and. (val .lt. vmin .or. val .gt. vmax) ) then
	    iout = iout + 1
	    if( debug ) write(6,*) 'out ',iout,i,val
	    write(999,*) 'out of range: ',iout,i,val
	  end if
	end do

	if( inan .gt. 0 .or. iout .gt. 0 ) then
	  write(6,*) '*** check1Dr: ',textgen," (",text,") "
	  write(6,*) 'total number of Nan found:       ',inan
	  write(6,*) 'total number out of range found: ',iout
	  write(6,*) 'full list can be found in file fort.999'
          stop 'error stop check1Dr'
	end if

	end
	  
c***************************************************************

	subroutine check2Dr(ndim,nlv,n,a,vmin,vmax,textgen,text)

c tests array for nan and strange values

	implicit none

	integer ndim,nlv,n
	real a(ndim,1)
	real vmin,vmax
	character*(*) textgen,text

	logical debug,bval
	integer inan,iout,i,l
	real val

	logical is_r_nan

        bval = vmin .lt. vmax
	debug = .true.
	debug = .false.
	inan = 0
	iout = 0

	do i=1,n
	  do l=1,nlv
	    val = a(l,i)
	    if( is_r_nan(val) ) then
	      inan = inan + 1
	      if( debug ) write(6,*) 'nan: ',inan,l,i,val
	      write(999,*) 'nan: ',inan,l,i,val
	    else if( bval .and. (val .lt. vmin .or. val .gt. vmax) ) then
	      iout = iout + 1
	      if( debug ) write(6,*) 'out of range: ',iout,l,i,val
	      write(999,*) 'out of range: ',iout,l,i,val
	    end if
	  end do
	end do

	if( inan .gt. 0 .or. iout .gt. 0 ) then
	  write(6,*) '*** check2Dr: ',textgen," (",text,") "
	  write(6,*) 'total number of Nan found:       ',inan
	  write(6,*) 'total number out of range found: ',iout
	  write(6,*) 'full list can be found in file fort.999'
          stop 'error stop check2Dr'
	end if

	end
	  
c***************************************************************
c***************************************************************
c***************************************************************

	subroutine checksum_2d(n1dim,n,nlv,levels,data,crc)

	integer n1dim			! dimension of first index
	integer n
	integer nlv			!nlv>0 use this  nlv<=0 use levels
	integer levels(n)
	real data(n1dim,n)
	integer crc

	integer a,b

	integer ivalue
	real value
	equivalence(value,ivalue)

	a = 1
	b = 0

	do i=1,n
	  lmax = nlv
	  if( lmax .le. 0 ) lmax = levels(i)
	  do l=1,lmax
	    value = data(l,i)
	    call checksum_i(a,b,ivalue)
	  end do
	end do

	crc = 256*256 * b + a

	end

c***************************************************************

	subroutine checksum_1d(n,data,crc)

	integer n
	real data(n)
	integer crc

	integer a,b

	integer ivalue
	real value
	equivalence(value,ivalue)

	a = 1
	b = 0

	do i=1,n
	  value = data(i)
	  call checksum_i(a,b,ivalue)
	end do

	crc = 256*256 * b + a

	end

c***************************************************************

	subroutine checksum_i(a,b,data)

	integer a,b,data

	integer MOD_ADLER
	parameter ( MOD_ADLER = 65521 )

	a = mod(a+data,MOD_ADLER)
	b = mod(a+b,MOD_ADLER)

	end

c***************************************************************

