
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2016,2018-2019  Georg Umgiesser
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
! 23.01.2015	ggu	changed VERS_7_1_4
! 10.06.2016	ggu	changed VERS_7_5_13
! 24.01.2018	ggu	changed VERS_7_5_41
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

	program filetype

c checks type of file

	use clo
	use basin

	implicit none

	integer iformat
	integer iunit,nvers,nvar
	integer nfile
	character*70 file

	integer ifileo

c--------------------------------------------------------------
c parameters and command line options
c--------------------------------------------------------------

        call clo_init('filetype','fem-file','1.0')

        call clo_add_info('returns type of file')

        call clo_parse_options(1)  !expecting (at least) 1 file after options

        nfile = clo_number_of_files()
	if( nfile < 1 ) stop

	call clo_get_file(1,file)

	nvar = 0
	iunit = 0
	nvers = 0
	iformat = -1

c-------------------------------------------------------------
c FEM file
c-------------------------------------------------------------

	call fem_file_test_fem_file(file,iformat)
	if( iformat >= 0 ) then
	  write(6,*) 'file format is FEM'
	  stop
	end if

c-------------------------------------------------------------
c NOS file
c-------------------------------------------------------------

	iunit = ifileo(0,file,'unformatted','old')
	call nos_is_nos_file(iunit,nvers)
	if( nvers > 0 ) then
	  write(6,*) 'file format is NOS'
	  stop
	end if
	if( iunit > 0 ) close(iunit)

c-------------------------------------------------------------
c OUS file
c-------------------------------------------------------------

	iunit = ifileo(0,file,'unformatted','old')
	call ous_is_ous_file(iunit,nvers)
	if( nvers > 0 ) then
	  write(6,*) 'file format is OUS'
	  stop
	end if
	if( iunit > 0 ) close(iunit)

c-------------------------------------------------------------
c BAS file
c-------------------------------------------------------------

	iunit = ifileo(0,file,'unformatted','old')
	call basin_test(iunit,nvers)
	if( nvers > 0 ) then
	  write(6,*) 'file format is BAS'
	  stop
	end if
	if( iunit > 0 ) close(iunit)

c-------------------------------------------------------------
c TS file
c-------------------------------------------------------------

	call ts_get_file_info(file,nvar)
	if( nvar > 0 ) then
	  write(6,*) 'file format is TS'
	  stop
	end if
	if( iunit > 0 ) close(iunit)

c-------------------------------------------------------------
c ETS file
c-------------------------------------------------------------

	iunit = ifileo(0,file,'formatted','old')
	call ets_is_ets_file(iunit,nvers)
	if( nvers > 0 ) then
	  write(6,*) 'file format is ETS'
	  stop
	end if
	if( iunit > 0 ) close(iunit)

c-------------------------------------------------------------
c unknown
c-------------------------------------------------------------

	write(6,*) 'file format is unknown'

	end

c*****************************************************************

