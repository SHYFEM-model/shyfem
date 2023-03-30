
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2019  Georg Umgiesser
!    Copyright (C) 2018  Christian Ferrarin
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

c revision log :
c
c 06.05.2015	ggu	noselab started
c 05.06.2015	ggu	many more features added
c 30.07.2015	ggu	shyelab started
c 14.09.2015	ggu	support for ext files added
c 05.10.2015	ggu	support for flx files added
c 09.10.2015	ggu	use last file to determine file type, call this routine
c 19.10.2015	ggu	changed VERS_7_3_6
c 25.05.2016	ggu	changed VERS_7_5_10
c 17.06.2016	ggu	changed VERS_7_5_15
c 09.10.2017	ggu	changed VERS_7_5_33
c 04.11.2017	ggu	new functionality tselab
c 14.11.2017	ggu	changed VERS_7_5_36
c 06.07.2018	ggu	changed VERS_7_5_48
c 30.08.2018	ccf	new functionality lgrelab
c 25.10.2018	ggu	changed VERS_7_5_51
c 16.02.2019	ggu	changed VERS_7_5_60
c
c**************************************************************

	program shyfile

	use clo
	use elabutil

c elaborates output file

	implicit none

	character*80 file,type

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

	call clo_get_last_file(file)
	call check_file_type(file,type)

	write(6,'(a)') trim(type)

        end

c***************************************************************

	subroutine check_file_type(file,type)

	use shyfile

	implicit none

	character*(*) file,type

	logical check_nos_file,check_ous_file
	logical check_ext_file,check_flx_file
	logical check_ts_file,fem_file_is_fem_file
	logical filex

	if( file == ' ') then
	  type = 'NONE'
	else if( .not. filex(file) ) then
	  type = 'NOTEXIST'
	else if( shy_is_shy_file(file) ) then
	  type = 'SHY'
	else if( shy_is_lgr_file(file) ) then
	  type = 'LGR'
	else if( check_nos_file(file) ) then
	  type = 'NOS'
	else if( check_ous_file(file) ) then
	  type = 'OUS'
	else if( check_ext_file(file) ) then
	  type = 'EXT'
	else if( check_flx_file(file) ) then
	  type = 'FLX'
	else if( fem_file_is_fem_file(file) ) then
	  type = 'FEM'
	else if( check_ts_file(file) ) then
	  type = 'TS'
	else
	  type = 'UNKNOWN'
	end if
	
        end

c***************************************************************

