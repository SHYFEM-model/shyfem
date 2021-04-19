
!--------------------------------------------------------------------------
!
!    Copyright (C) 2005,2010,2012,2014-2015,2018-2020  Georg Umgiesser
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
! 09.03.2020	ggu	prepared for mercury restart
! 20.03.2020    ggu     adjusted for new routine calls
! 27.03.2021    ggu     new option -checkval
! 14.04.2021    ggu     bug fix - atime was integer

!******************************************************************

	program rstinf

	implicit none

	logical bread,bcheckval
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
     +                  '     iflag'
!                 12345678901234567890123456789012345678901234567890123
	title2 = '   irec                         atime     date'

!-------------------------------------------------------------------
! initialize and open file
!-------------------------------------------------------------------

	nread = 0
	iunit = 1
	file = ' '

        call rst_init(file,bcheckval)

	bread = bcheckval

	if( bread ) call open_for_read(file)

	open(iunit,file=file,status='old',form='unformatted')

!-------------------------------------------------------------------
! loop on records
!-------------------------------------------------------------------

	do

	  call rst_skip_record(iunit,atime,nvers,nrec
     +					,nknr,nelr,nlvr,iflag,ierr)
	  if( ierr .ne. 0 ) exit

	  if( nread == 0 ) then
            write(6,1000) trim(title1)
            write(6,1010) nvers,nrec,nknr,nelr,nlvr,iflag
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

        write(6,*)
	write(6,*) 'Number of records read: ',nread
	call dts_format_abs_time(atime_anf,aline)
	write(6,*) 'Initial time in file:   ',atime_anf,aline
	call dts_format_abs_time(atime_end,aline)
	write(6,*) 'Final time in file:     ',atime_end,aline
	write(6,*)
        write(6,1000) trim(title1)
        write(6,1010) nvers,nrec,nknr,nelr,nlvr,iflag
        write(6,*)

	call write_flags(iflag)

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

	subroutine open_for_read(file)

	use basin
	use levels
	use mod_hydro
	use mod_geom_dynamic
	use mod_ts
	use mod_hydro_vel

	implicit none

	character*(*) file

	integer iunit
	double precision atime
	integer nvers,nrec
	integer nk,ne,nl,iflag,ierr

	iunit = 1
	open(iunit,file=file,status='old',form='unformatted')
	call rst_skip_record(iunit,atime,nvers,nrec
     +					,nk,ne,nl,iflag,ierr)

	call basin_init(nk,ne)
	call levels_init(nk,ne,nl)
	call mod_hydro_init(nk,ne,nl)
	call mod_geom_dynamic_init(nk,ne)
	call mod_ts_init(nk,nl)
	call mod_hydro_vel_init(nk,ne,nl)

	close(iunit)

	end

!******************************************************************

        subroutine rst_init(rstfile,bcheckval)

        use clo

        implicit none

        character*(*) rstfile
	logical bcheckval

        call shyfem_copyright('rstinf - info on restart file')

        call clo_init('rstinf','rstfile','1.2')

        call clo_add_info('returns info on records of restart file')

        call clo_add_option('checkval',.false.,'check NaNs in file')
 
        call clo_parse_options

        call clo_get_option('checkval',bcheckval)
 
        call clo_check_files(1)
        call clo_get_file(1,rstfile)

        end

!******************************************************************

