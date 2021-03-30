
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2008,2010-2012,2015-2016,2015-2016  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
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

! revision log :
!
! 27.03.2003	ggu	new routines for NaN check and value test
! 03.09.2003	ggu	check routines customized
! 05.12.2003	ggu	in check[12]Dr only check val if vmin!=vmax
! 06.12.2008	ggu	check for NaN changed (Portland compiler)
! 23.03.2010	ggu	changed v6.1.1
! 08.10.2010	ggu	changed VERS_6_1_13
! 23.03.2011	ggu	changed VERS_6_1_21
! 15.07.2011	ggu	new routines for checksum (CRC)
! 30.03.2012	ggu	changed VERS_6_1_51
! 19.01.2015	ggu	changed VERS_7_1_3
! 20.11.2015	ggu	changed VERS_7_3_15
! 18.12.2015	ggu	changed VERS_7_3_17
! 25.05.2016	ggu	changed VERS_7_5_10
! 05.10.2018	ggu	eliminated equivalent statement
! 16.10.2018	ggu	changed VERS_7_5_50
! 02.02.2019	ggu	new routines is_inf, is_nonumber, adjourned is_nan
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 22.09.2020	ggu	added test for arrays to check for Nan
! 27.03.2021	ggu	debug output routines added
! 30.03.2021	ggu	new routine set_debug_unit()
!
! notes :
!
! use of debug routines:
!
! in the calling routine you must set or unset the debug:
!
! use mod_debug
! ...
! call set_debug
!   here call your_routine
! call clear_debug
!
! in the routine where you want to debug:
!
! use mod_debug
! ...
! if( is_debug() ) then
!   here insert your debug statements
! end if
!
! to check multi-dimensional arrays for Nan, use the reshape intrinsic:
!
! is_nan( reshape(array,(/size(array)/)) )
!
!***************************************************************

!==================================================================
        module mod_debug
!==================================================================

	implicit none

	logical, save, private :: debug_intern = .false.

	integer, save, private :: iunit = 666

	integer, parameter, private :: integer_type = 1
	integer, parameter, private :: real_type = 2
	integer, parameter, private :: double_type = 3

        INTERFACE is_nan
        MODULE PROCEDURE is_d_nan, is_r_nan, is_i_nan
     +			,is_d_array1_nan
        END INTERFACE

        INTERFACE is_inf
        MODULE PROCEDURE is_d_inf, is_r_inf, is_i_inf
        END INTERFACE

        INTERFACE is_nonumber
        MODULE PROCEDURE is_d_nonumber, is_r_nonumber, is_i_nonumber
        END INTERFACE

        INTERFACE write_debug_record
        MODULE PROCEDURE write_debug_record_i1
     +			,write_debug_record_r1
     +			,write_debug_record_d1
     +			,write_debug_record_r2
     +			,write_debug_record_d2
        END INTERFACE

        INTERFACE write_debug_time
        MODULE PROCEDURE write_debug_time_intern
        END INTERFACE

        INTERFACE write_debug_final
        MODULE PROCEDURE write_debug_final_intern
        END INTERFACE

!==================================================================
        contains
!==================================================================

	subroutine assign_debug(bdebug)

	implicit none

	logical bdebug

	debug_intern = bdebug

	end subroutine

!***************************************************************

	subroutine set_debug

	implicit none

	debug_intern = .true.

	end subroutine

!***************************************************************

	subroutine clear_debug

	implicit none

	debug_intern = .false.

	end subroutine

!***************************************************************

	function is_debug()

	implicit none

	logical is_debug

	is_debug = debug_intern

	end function

!***************************************************************
!***************************************************************
!***************************************************************
! NaN routines
!***************************************************************
!***************************************************************
!***************************************************************

	function is_d_array1_nan(val)

! tests 1D array for NaN

	implicit none

	double precision val(:)
	logical is_d_array1_nan

	integer i

	is_d_array1_nan = .true.

	do i=1,size(val)
	  if( is_d_nan(val(i)) ) return
	end do

	is_d_array1_nan = .false.

	end function

!***************************************************************

	function is_d_nan(val)

! tests val for NaN

	implicit none

	double precision val
	logical is_d_nan
        integer itot

        itot = 0

	if( val .gt. 0. ) itot = itot + 1
	if( val .lt. 0. ) itot = itot + 1
	if( val .eq. 0. ) itot = itot + 1

	is_d_nan = itot .ne. 1

	end function

!***************************************************************

	function is_r_nan(val)

! tests val for NaN

	implicit none

	real val
	logical is_r_nan
        integer itot

        itot = 0

	if( val .gt. 0. ) itot = itot + 1
	if( val .lt. 0. ) itot = itot + 1
	if( val .eq. 0. ) itot = itot + 1

	is_r_nan = itot .ne. 1

	end function

!***************************************************************

	function is_i_nan(val)

! tests val for NaN

	implicit none

	integer val
	logical is_i_nan
        integer itot

        itot = 0

	if( val .gt. 0 ) itot = itot + 1
	if( val .lt. 0 ) itot = itot + 1
	if( val .eq. 0 ) itot = itot + 1

	is_i_nan = itot .ne. 1

	end function

!***************************************************************
!***************************************************************
!***************************************************************
! infinity routines
!***************************************************************
!***************************************************************
!***************************************************************

	function is_i_inf(val)

! tests val for infinity

	implicit none

	integer val
	logical is_i_inf

	is_i_inf = ( val-1 == val )

	end function

!***************************************************************

	function is_r_inf(val)

! tests val for infinity

	implicit none

	real val
	logical is_r_inf

	is_r_inf = ( val-1. == val )

	end function

!***************************************************************

	function is_d_inf(val)

! tests val for infinity

	implicit none

	double precision val
	logical is_d_inf

	is_d_inf = ( val-1. == val )

	end function

!***************************************************************
!***************************************************************
!***************************************************************
! no number routines
!***************************************************************
!***************************************************************
!***************************************************************

	function is_i_nonumber(val)

! tests val for infinity

	implicit none

	integer val
	logical is_i_nonumber

	is_i_nonumber = ( is_nan(val) .or. is_inf(val) )

	end function

!***************************************************************

	function is_r_nonumber(val)

! tests val for infinity

	implicit none

	real val
	logical is_r_nonumber

	is_r_nonumber = ( is_nan(val) .or. is_inf(val) )

	end function

!***************************************************************

	function is_d_nonumber(val)

! tests val for infinity

	implicit none

	double precision val
	logical is_d_nonumber

	is_d_nonumber = ( is_nan(val) .or. is_inf(val) )

	end function

!***************************************************************
!***************************************************************
!***************************************************************
! debug output routines -> read with check_debug
!***************************************************************
!***************************************************************
!***************************************************************

	subroutine set_debug_unit(iu)
	implicit none
	integer iu
	iunit = iu
	end

	subroutine write_debug_time_intern(dtime)
	implicit none
	double precision dtime
        write(iunit) dtime
	end

	subroutine write_debug_final_intern
	implicit none
        write(iunit) 0,0,0
	end

        subroutine write_debug_record_i1(val,text)
        implicit none
        integer val(:)
        character*(*) text
	integer ntot,nfirst
        character*80 text1
	ntot = size(val)
	nfirst = 1
        text1=text
        write(iunit) ntot,nfirst,integer_type
        write(iunit) text1
        write(iunit) val
        end

        subroutine write_debug_record_r1(val,text)
        implicit none
        real val(:)
        character*(*) text
	integer ntot,nfirst
        character*80 text1
	ntot = size(val)
	nfirst = 1
        text1=text
        write(iunit) ntot,nfirst,real_type
        write(iunit) text1
        write(iunit) val
        end

        subroutine write_debug_record_d1(val,text)
        implicit none
        double precision val(:)
        character*(*) text
	integer ntot,nfirst
        character*80 text1
	ntot = size(val)
	nfirst = 1
        text1=text
        write(iunit) ntot,nfirst,double_type
        write(iunit) text1
        write(iunit) val
        end

        subroutine write_debug_record_r2(val,text)
        implicit none
        real val(:,:)
        character*(*) text
	integer ntot,nfirst
        character*80 text1
	ntot = size(val)
	nfirst = size(val,1)
        text1=text
        write(iunit) ntot,nfirst,real_type
        write(iunit) text1
        write(iunit) val
        end

        subroutine write_debug_record_d2(val,text)
        implicit none
        double precision val(:,:)
        character*(*) text
	integer ntot,nfirst
        character*80 text1
	ntot = size(val)
	nfirst = size(val,1)
        text1=text
        write(iunit) ntot,nfirst,double_type
        write(iunit) text1
        write(iunit) val
        end

!==================================================================
        end module mod_debug
!==================================================================

!***************************************************************
!***************************************************************
!***************************************************************
! various check routines -> use check1Dr(), check2Dr()
!***************************************************************
!***************************************************************
!***************************************************************

	subroutine nantest(n,a,text)

! tests array for nan

	use mod_debug

	implicit none

	integer n
	real a(n)
	character*(*) text

	integer i

	do i=1,n
          if( is_r_nan( a(i) ) ) goto 99
	end do

	return
   99	continue
	write(6,*) 'NAN found while testing ',text
	write(6,*) n,i,a(i)
	stop 'error stop nantest: NAN found'
	end

!***************************************************************

	subroutine valtest(nlv,n,valmin,valmax,a,textgen,text)

! tests array for nan - use other routines below

	use mod_debug

	implicit none

	integer nlv,n
	real valmin,valmax
	real a(nlv*n)
	character*(*) textgen,text

	call check2Dr(nlv,nlv,n,a,valmin,valmax,textgen,text)

	end

!***************************************************************

	subroutine check1Dr(n,a,vmin,vmax,textgen,text)

! tests array for nan and strange values

	use mod_debug

	implicit none

	integer n
	real a(n)
	real vmin,vmax
	character*(*) textgen,text

	logical debug,bval
	integer inan,iout,i
	real val

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
	  
!***************************************************************

	subroutine check2Dr(nlvdi,nlv,n,a,vmin,vmax,textgen,text)

! tests array for nan and strange values

	use mod_debug

	implicit none

	integer nlvdi,nlv,n
	real a(nlvdi,n)
	real vmin,vmax
	character*(*) textgen,text

	logical debug,bval
	integer inan,iout,i,l
	real val

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
	  
!***************************************************************
!***************************************************************
!***************************************************************
! checksum routines
!***************************************************************
!***************************************************************
!***************************************************************

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
	!equivalence(value,ivalue)

	a = 1
	b = 0

	do i=1,n
	  lmax = nlv
	  if( lmax .le. 0 ) lmax = levels(i)
	  do l=1,lmax
	    value = data(l,i)
	    ivalue = transfer(value,1)
	    call checksum_i(a,b,ivalue)
	  end do
	end do

	crc = 256*256 * b + a

	end

!***************************************************************

	subroutine checksum_1d(n,data,crc)

	integer n
	real data(n)
	integer crc

	integer a,b

	integer ivalue
	real value
	!equivalence(value,ivalue)

	a = 1
	b = 0

	do i=1,n
	  value = data(i)
	  ivalue = transfer(value,1)
	  call checksum_i(a,b,ivalue)
	end do

	crc = 256*256 * b + a

	end

!***************************************************************

	subroutine checksum_i(a,b,data)

	integer a,b,data

	integer MOD_ADLER
	parameter ( MOD_ADLER = 65521 )

	a = mod(a+data,MOD_ADLER)
	b = mod(a+b,MOD_ADLER)

	end

!***************************************************************
!***************************************************************
!***************************************************************
! other routines
!***************************************************************
!***************************************************************
!***************************************************************

	subroutine divide_by_zero(iaux)

! raises a division by zero error

	implicit none

	integer iaux

        iaux = 0
        iaux = 1/iaux

	end

!***************************************************************
!***************************************************************
!***************************************************************
! test routines
!***************************************************************
!***************************************************************
!***************************************************************

	subroutine test_debug

	end

!***************************************************************
!	program main_test_debug
!	call test_debug
!	end
!***************************************************************

