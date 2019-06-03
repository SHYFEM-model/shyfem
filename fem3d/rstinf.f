
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

! info on restart file

! revision log :
!
! 01.05.2005	ggu	written from scratch
! 23.03.2010	ggu	changed v6.1.1
! 28.09.2010	ggu	changed VERS_6_1_11
! 14.02.2012	ggu	changed VERS_6_1_44
! 29.08.2012	ggu	changed VERS_6_1_56
! 26.11.2014	ggu	changed VERS_7_0_7
! 05.06.2015	ggu	changed VERS_7_1_12
! 05.11.2015	ggu	changed VERS_7_3_12
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.02.2019	ggu	changed VERS_7_5_60
! 03.05.2019	ggu	adapted
! 21.05.2019	ggu	changed VERS_7_5_62

!******************************************************************

	program rstinf

	implicit none

	integer iunit,it,nvers,nrec,nknr,nelr,nlvr,iflag,ierr,ic
	integer nread
	double precision atime
	double precision atime_anf
	double precision atime_end
	character*20 aline
	character*80 file
	character*80 title1
	character*80 title2

!-------------------------------------------------------------------
! create title strings
!-------------------------------------------------------------------

	title1 = 'version nrec       nkn       nel       nlv' //
     +                  '     iconz     iflag'
!                 12345678901234567890123456789012345678901234567890123
	title2 = '   irec                         atime     date'

!-------------------------------------------------------------------
! initialize and open file
!-------------------------------------------------------------------

	nread = 0
	iunit = 1
	file = ' '

        call rst_init(file)

	open(iunit,file=file,status='old',form='unformatted')

!-------------------------------------------------------------------
! loop on records
!-------------------------------------------------------------------

	do

	call rst_skip_record(iunit,atime,nvers,nrec
     +					,nknr,nelr,nlvr,ic,iflag,ierr)
	if( ierr .ne. 0 ) exit

	if( nread == 0 ) then
          write(6,1000) trim(title1)
          write(6,1010) nvers,nrec,nknr,nelr,nlvr,ic,iflag
          write(6,*)
          write(6,1001) trim(title2)
	  atime_anf = atime
	end if

	nread = nread + 1
	call dts_format_abs_time(atime,aline)
	write(6,1011) nread,atime,aline
	atime_end = atime

	end do

	if( ierr > 0 ) stop 'error stop rstinf: error reading record'

!-------------------------------------------------------------------
! final message
!-------------------------------------------------------------------

        write(6,1000) trim(title1)
        write(6,1010) nvers,nrec,nknr,nelr,nlvr,ic,iflag
        write(6,*)
        write(6,1001) trim(title2)
	write(6,*)
	write(6,*) 'Number of records read: ',nread
	call dts_format_abs_time(atime_anf,aline)
	write(6,*) 'Initial time in file:   ',atime_anf,aline
	call dts_format_abs_time(atime_end,aline)
	write(6,*) 'Final time in file:     ',atime_end,aline
	write(6,*)
	write(6,*) 'Meaning of iflag:'
	write(6,*) '         1          hydro'
	write(6,*) '        10          depth'
	write(6,*) '       100          ibarcl (T/S/rho)'
	write(6,*) '      1000          iconz (cnv/conzv)'
	write(6,*) '     10000          iwvert (wlnv)'
	write(6,*) '    100000          ieco (ecological variables)'

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

	stop
 1000	format(a)
 1001	format(a)
 1010	format(i7,i5,5i10)
 1011	format(i7,f30.2,5x,a20)
	end

!******************************************************************

        subroutine rst_init(rstfile)

        use clo

        implicit none

        character*(*) rstfile

        call shyfem_copyright('rstinf - info on restart file')

        call clo_init('rstinf','rstfile','1.2')

        call clo_add_info('returns info on records of restart file')

        call clo_parse_options

        call clo_check_files(1)
        call clo_get_file(1,rstfile)

        end

!******************************************************************

