
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2019-2020  Georg Umgiesser
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
! 05.11.2015	ggu	changed VERS_7_3_12
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 20.03.2020	ggu	new function bit10_return_value()

!*******************************************************************
!*******************************************************************
!*******************************************************************

        function bit10_internal(iflag,ipos,inew)

! internal routine - do not call directly

        implicit none

        integer bit10_internal
        integer iflag           !flag containing features desired
	integer ipos		!position wanted
	integer inew		!new value (-1: do not change)

	integer ipow,ifact,if,ifh,ibit

        ipow = ipos - 1
        ifact = 10**ipow
        if = iflag / ifact
        ifh = 10 * ( if / 10 )
        ibit = if - ifh

	bit10_internal = ibit

	if( inew < 0 ) return
	if( inew > 9 ) stop 'error stop bit10_internal: inew'

	iflag = iflag + ifact * ( inew - ibit )

        end

!*******************************************************************
!*******************************************************************
!*******************************************************************

        function bit10_extract_value(iflag,ipos)

! extracts value in iflag at position ipos

        implicit none

        integer bit10_extract_value
        integer iflag           !flag containing features desired
	integer ipos		!position wanted

        integer bit10_internal

        bit10_extract_value = bit10_internal(iflag,ipos,-1)

	end

!*******************************************************************

        function bit10_return_value(ipos)

! returns value for position ipos

        implicit none

	integer bit10_return_value
	integer ipos		!position wanted

        integer iflag           !flag containing features desired
	integer iaux
        integer bit10_internal

	iflag = 0
        iaux = bit10_internal(iflag,ipos,1)

	bit10_return_value = iflag

	end

!*******************************************************************

        subroutine bit10_set_value(iflag,ipos,inew)

! sets value in iflag at position ipos to inew

        implicit none

        integer iflag           !flag containing features desired
	integer ipos		!position wanted
	integer inew		!new value

	integer iaux
        integer bit10_internal

        iaux = bit10_internal(iflag,ipos,inew)

	end

!*******************************************************************

        subroutine bit10_set_pos(iflag,ipos)

! sets value in iflag at position ipos

        implicit none

        integer iflag           !flag containing features desired
	integer ipos		!position wanted

	integer iaux
        integer bit10_internal

        iaux = bit10_internal(iflag,ipos,1)

	end

!*******************************************************************

        subroutine bit10_clear_pos(iflag,ipos)

! clears value in iflag at position ipos

        implicit none

        integer iflag           !flag containing features desired
	integer ipos		!position wanted

	integer iaux
        integer bit10_internal

        iaux = bit10_internal(iflag,ipos,0)

	end

!*******************************************************************

        function bit10_is_set(iflag,ipos)

! checks if ipos in iflag is set

        implicit none

        logical bit10_is_set
        integer iflag           !flag containing features desired
	integer ipos		!position wanted

        integer bit10_extract_value

	bit10_is_set = ( bit10_extract_value(iflag,ipos) /= 0 )
	
	end

!*******************************************************************
!*******************************************************************
!*******************************************************************
! testing routines
!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine bit10_assert(text,bflag)

	implicit none

	character*(*) text
	logical bflag

	if( bflag ) then
	  write(6,*) trim(text),' : ','ok'
	else
	  write(6,*) trim(text),' : ','error'
	end if

	end

!*******************************************************************

	subroutine bit10_test

	implicit none

	integer iflag,ibit
	logical bit10_is_set
        integer bit10_extract_value

	iflag = 101010

	write(6,*) 'testing routines bit10'
	write(6,*) 'iflag = ',iflag
	call bit10_assert('bit 1 is clear',.not.bit10_is_set(iflag,1))
	call bit10_assert('bit 2 is set',bit10_is_set(iflag,2))
	call bit10_clear_pos(iflag,2)
	call bit10_assert('bit 2 cleared',iflag==101000)
	call bit10_set_pos(iflag,5)
	call bit10_assert('bit 5 set',iflag==111000)
        ibit = bit10_extract_value(iflag,3)
	call bit10_assert('bit 3 is clear',ibit==0)
        ibit = bit10_extract_value(iflag,4)
	call bit10_assert('bit 4 is set',ibit==1)
	write(6,*) 'iflag = ',iflag

	end

!*******************************************************************

!	program bit10_main
!	call bit10_test
!	end

!*******************************************************************

