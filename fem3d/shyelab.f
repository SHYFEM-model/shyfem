c
c revision log :
c
c 06.05.2015    ggu     noselab started
c 05.06.2015    ggu     many more features added
c 30.07.2015    ggu     shyelab started
c 14.09.2015    ggu     support for ext files added
c 05.10.2015    ggu     support for flx files added
c 09.10.2015    ggu     use last file to determine file type, call this routine
c 04.11.2017    ggu     new functionality tselab
c
c**************************************************************

	program shyelab

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

	if( type == 'NONE' ) then
	  call elabutil_init('NONE','shyelab')
	else if( type == 'NOTEXIST' ) then
	  write(6,*) 'file does not exists: ',trim(file)
	else if( type == 'SHY' ) then
	  call shyelab1
	else if( type == 'NOS' ) then
	  write(6,*) 'file is of NOS type'
	  write(6,*) 'please convert to SHY format or use noselab'
	else if( type == 'OUS' ) then
	  write(6,*) 'file is of OUS type'
	  write(6,*) 'please convert to SHY format or use ouselab'
	else if( type == 'EXT' ) then
	  call extelab
	else if( type == 'FLX' ) then
	  call flxelab
	else if( type == 'FEM' ) then
	  call femelab
	else if( type == 'TS' ) then
	  call tselab
	else if( type == 'UNKNOWN' ) then
	  write(6,*) 'unknown file type: ',trim(file)
	else
	  write(6,*) 'inconsistent file type: ',trim(type)
	  write(6,*) 'internal error (1)'
	end if
	
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

