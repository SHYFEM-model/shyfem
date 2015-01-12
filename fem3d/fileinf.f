
	program fileinf

c checks type of file

	implicit none

	integer iformat
	integer iunit,nvers,nvar
	character*70 file

	integer ifileo

        call parse_command_line(file)

	nvar = 0
	iunit = 0
	nvers = 0
	iformat = -1

c-------------------------------------------------------------
c FEM file
c-------------------------------------------------------------

	call fem_file_is_fem_file(file,iformat)
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
	call sp13test(iunit,nvers)
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

        subroutine parse_command_line(infile)

        implicit none

        character*(*) infile

        integer i,nc
        character*50 aux

        infile = ' '

        nc = command_argument_count()

        if( nc > 0 ) then
          call get_command_argument(1,infile)
          return
        end if

        write(6,*) 'Usage: fileinf [options] file'

	stop
	end

c*****************************************************************

