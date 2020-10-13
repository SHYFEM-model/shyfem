
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013-2015,2017,2019  Georg Umgiesser
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
! 09.10.2020    ggu     adjusted using routines from module

!**********************************************************

	program offinf

	implicit none

	integer iu,ios
	character*80 filename

        call off_init(filename)

	iu = 1
	open(iu,file=filename,status='old',form='unformatted',iostat=ios)
	if( ios /= 0 ) stop 'error stop offinf: opening file'

	!call offinf_old(iu)
	!rewind(iu)

	call offinf_new(iu)

	close(iu)

	end

!**********************************************************

	subroutine offinf_old(iu)

! shows content of offline data file

	implicit none

	integer iu

	integer it,nkn,nel,nrec,i,type,irecs
	integer ios
	real buffer(1)

	nrec = 0
	irecs = 9

	do
          read(iu,iostat=ios) it,nkn,nel,type
	  if( ios /= 0 ) exit
	  nrec = nrec + 1
	  write(6,*) nrec,it,nkn,nel,type
	  if( type /= 3 ) stop 'error stop offinf: type /= 3'

	  do i=1,irecs
	    read(iu) buffer
	  end do
	end do

	if( ios > 0 ) stop 'error stop offinf: reading file'

	write(6,*) 'finished reading: ',nrec

	end

!**********************************************************

	subroutine offinf_new(iu)

! shows content of offline data file

	use mod_offline

	implicit none

	integer iu

	integer it,nkn,nel,nlv,nrec
	integer ierr,ig
	integer itfirst,itlast
	character*80 header

	header = '      nread        time'
	nrec = 0
	ig = 1

	call off_peek_header(iu,it,nkn,nel,nlv,ierr)
	if( ierr /= 0 ) stop 'error stop offinf: reading header'
	call mod_offline_init(nkn,nel,nlv)
	write(6,*) 'parameters: ',nkn,nel,nlv
	write(6,*) trim(header)
	itfirst = it

	do
	  call off_read_record(iu,ig,it,ierr)
	  if( ierr /= 0 ) exit
	  nrec = nrec + 1
	  write(6,*) nrec,it
	  itlast = it
	end do

	if( ierr > 0 ) stop 'error stop offinf: reading file'

        write(6,*) 'records read =    ',nrec
        write(6,*) 'itfirst =         ',itfirst
        write(6,*) 'itlast =          ',itlast

	end

!**********************************************************

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

!**********************************************************

