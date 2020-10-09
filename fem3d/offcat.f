
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013-2015,2017,2019-2020  Georg Umgiesser
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
! 09.10.2020	ggu	adjusted using routines from module

!**********************************************************

	program offcat

! shows content of offline data file

	use mod_offline

	implicit none

	logical bwrite
	integer nkn,nel,nlv,type
	integer ifile,nrec,nout
	integer ierr,ios,ig
	integer iu,iuout
	integer it,itfirst,itlast
	character*80 filename

        call off_init(filename)

	nrec = 0
	nout = 0
	ifile = 0
	iu = 1
	iuout=2
	ig = 1
	itfirst = 0
	itlast = 0

	open(iuout,file='out.off',status='unknown',form='unformatted')

	call off_open_next_file(iu,ierr)
	if( ierr /= 0 ) stop 'error stop offcat: no files'

	call off_peak_header(iu,it,nkn,nel,nlv,ierr)
	if( ierr /= 0 ) stop 'cannot read file'
	call mod_offline_init(nkn,nel,nlv)
	write(6,*) 'parameters: ',nkn,nel,nlv
	itfirst = it
	itlast = itfirst - 1

	do	!loop on files

	  ifile = ifile + 1

	  do	!loop on records

	    call off_read_record(iu,ig,it,ierr)
	    if( ierr /= 0 ) exit
	    nrec = nrec + 1
	    bwrite = ( it > itlast )
	    write(6,*) ifile,nrec,nout,it,bwrite
	    if( .not. bwrite ) cycle
	    nout = nout + 1
	    call off_write_record(iuout,it)
	    itlast = it

	  end do

	  if( ierr > 0 ) stop 'error stop offcat: reading file'

	  call off_open_next_file(iu,ierr)
	  if( ierr /= 0 ) exit
	end do

	write(6,*) 'records read = ',nrec
	write(6,*) 'records written = ',nout
	write(6,*) 'itfirst = ',itfirst
	write(6,*) 'itlast = ',itlast

	close(iu)
	close(iuout)

	end

!**********************************************************

        subroutine off_init(offfile)

        use clo

        implicit none

        character*(*) offfile

        call shyfem_copyright('offcat - concatenates offline files')

        call clo_init('offcat','offfile','1.0')

	call clo_add_info('concatenates records of offline files')

        call clo_parse_options

        end

!**********************************************************

	subroutine off_open_next_file(iu,ierr)

        use clo

        implicit none

	integer iu,ierr

	logical bopen
        character*80 offfile

	ierr = -1

        call clo_get_next_file(offfile)
	if( offfile == ' ' ) return

	inquire(unit=iu,opened=bopen)
	if( bopen ) close(iu)

	write(6,*) 'opening file ',trim(offfile)
	open(iu,file=offfile,status='old',form='unformatted',iostat=ierr)
	if( ierr /= 0 ) stop 'error stop offcat: opening file'

	end

!**********************************************************

