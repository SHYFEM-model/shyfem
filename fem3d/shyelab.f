c
c revision log :
c
c 06.05.2015    ggu     noselab started
c 05.06.2015    ggu     many more features added
c 30.07.2015    ggu     shyelab started
c
c**************************************************************

	program shyelab

	use clo
	use elabutil

c elaborates output file

	implicit none

	integer nc
	character*80 file

	logical check_nos_file,check_ous_file

c--------------------------------------------------------------

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

	call elabutil_init('SHY')

        nc = command_argument_count()
        if( nc .le. 0 ) then
          write(6,*) 'Usage: shyelab file'
          stop 'error stop shyelab: no files given'
        end if

        call clo_get_file(1,file)

	if( check_nos_file(file) ) then
	  write(6,*) 'file is of NOS type'
	  call noselab
	else if( check_ous_file(file) ) then
	  write(6,*) 'file is of OUS type'
	  call ouselab
	else
	  write(6,*) 'cannot yet handle this file type'
	end if
	
        end

c***************************************************************

