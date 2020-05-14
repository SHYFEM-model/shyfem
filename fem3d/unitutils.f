
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2020  Georg Umgiesser
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
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 18.03.2020	ggu	handle case if no node number is given (j=0)

!**************************************************************************

	module shyelab_unit

	integer, save :: iunit = 100

	end module shyelab_unit

!***************************************************************

	subroutine get_new_unit(iu)

	use shyelab_unit

	implicit none

	integer iu
	logical bopen

	do
	  iunit = iunit + 1
	  inquire(unit=iunit,opened=bopen)
	  if( .not. bopen ) exit
	end do

	iu = iunit

	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine make_iunit_name(short,modi,dim,j,iu)

	implicit none

	character*(*) short,modi,dim
	integer j
	integer iu

	logical bopen
	character*5 numb
	character*80 dimen
	character*80 name

	if( j <= 0 ) then	!no node number given
	  numb = 'txt'
	else
          write(numb,'(i5)') j
          numb = adjustl(numb)
	end if

	dimen = '.' // trim(dim) // '.'
	name = trim(short)//trim(modi)//trim(dimen)//numb

	call get_new_unit(iu)

	inquire(file=name,opened=bopen)
	if( bopen ) goto 99
	inquire(unit=iu,opened=bopen)
	if( bopen ) goto 98

        !write(6,*) 'opening file : ',iu,'  ',trim(name)
        open(iu,file=name,form='formatted',status='unknown')

	return
   99	continue
	write(6,*) 'file already open: ',trim(name)
	stop 'error stop make_iunit_name: internal error (1)'
   98	continue
	write(6,*) 'unit already open: ',iu
	stop 'error stop make_iunit_name: internal error (2)'
	end

!***************************************************************

