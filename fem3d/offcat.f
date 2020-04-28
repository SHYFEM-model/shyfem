
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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
! 13.06.2013	ggu	changed VERS_6_1_65
! 05.05.2014	ggu	changed VERS_6_1_74
! 05.11.2015	ggu	changed VERS_7_3_12
! 31.03.2017	ggu	changed VERS_7_5_24
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 28.04.2020	ggu	written for Michol

c**********************************************************

	program offcat

c shows content of offline data file

	use mod_offline

	implicit none

	integer it,nkn,nel,nlv,nrec,iu,i,type,irecs
	integer ierr,ios
	character*60 name
	!double precision buffer(5)
	real buffer(1)

        call off_init(name)

	nrec = 0
	iu = 1
	irecs = 9

	open(iu,file=name,status='old',form='unformatted',iostat=ios)
	if( ios /= 0 ) stop 'error stop offinf: opening file'

	call off_read_header(iu,nkn,nel,nlv,ierr)
	rewind(iu)
	if( ierr /= 0 ) stop 'cannot read file'
	call mod_offline_init(nkn,nel,nlv)
	write(6,*) nkn,nel,nlv

	do

	  call off_read_record(iu,1,it,ierr)
	  if( ierr /= 0 ) exit
	  write(6,*) it

	end do

	if( ierr > 0 ) stop 'error stop offinf: reading file'

	close(iu)

	end

c**********************************************************

        subroutine off_init(offfile)

        use clo

        implicit none

        character*(*) offfile

        call shyfem_copyright('offinf - info on offline file')

        call clo_init('offinf','offfile','1.0')

	call clo_add_info('returns info on records of offline file')

        call clo_parse_options

        call clo_check_files(1)
        call clo_get_file(1,offfile)

        end

c**********************************************************

