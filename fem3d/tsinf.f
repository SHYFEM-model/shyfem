
	program tsinf

c writes info on ts file

	use clo

	implicit none

	character*50 infile
	integer nfile
	integer i
	logical bdebug,bwrite,bout,btmin,btmax
	logical bquiet
	double precision tmin,tmax

	bdebug = .true.
	bdebug = .false.

c--------------------------------------------------------------
c parameters and command line options
c--------------------------------------------------------------

	call clo_init('tsinf','ts-file(s)','1.1')

	call clo_add_option('write',.false.,'write min/max of values')
	call clo_add_option('out',.false.,'create output file')
	call clo_add_option('quiet',.false.,'do not be verbose')
	call clo_add_option('tmin time',-1
     +				,'only process starting from time')
	call clo_add_option('tmax time',-1
     +				,'only process up to time')

	call clo_parse_options(1)  !expecting (at least) 1 file after options

	call clo_get_option('write',bwrite)
	call clo_get_option('out',bout)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('tmin',tmin)
	call clo_get_option('tmax',tmax)

	btmin = tmin .ne. -1.
	btmax = tmax .ne. -1.

	if( bout ) stop 'error stop: -out not yet implemented'
	if( btmin ) stop 'error stop: -tmin not yet implemented'
	if( btmax ) stop 'error stop: -tmax not yet implemented'

	nfile = clo_number_of_files()

	if( bdebug ) then
	  write(6,*) nfile
	  write(6,*) bwrite,bout,btmin,btmax
	  write(6,*) tmin,tmax
	end if

c--------------------------------------------------------------
c loop on files
c--------------------------------------------------------------

	do i=1,nfile
	  call clo_get_file(i,infile)
	  if( infile .ne. ' ' ) call tsinf_file(infile)
	end do

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine tsinf_file(infile)

c writes info on ts file

	use clo

	implicit none

	character*(*) infile

	character*50 name,string
	integer np,iunit,iout
	integer nvers,lmax,nvar,ntype
	integer nvar0,lmax0,np0
	integer idt,idtact
	double precision dtime,tmin,tmax,dtime0
	double precision atime,atimeold,atimeanf,atimeend
	real dmin,dmax
	integer ierr
	integer nfile
	integer irec,i,ich,nrecs
	integer datetime(2)
	logical bdebug,bfirst,bskip,bwrite,bout,btmin,btmax,boutput
	logical bquiet
	character*20 line
	real,allocatable :: data(:)
	real,allocatable :: data_minmax(:,:)

	bdebug = .true.
	bdebug = .false.

	datetime = 0
	irec = 0

	call clo_get_option('write',bwrite)
	call clo_get_option('quiet',bquiet)

c--------------------------------------------------------------
c open file
c--------------------------------------------------------------

	if( infile .eq. ' ' ) stop

	nvar = 0
	call ts_open_file(infile,nvar,datetime,iunit)
	if( iunit .le. 0 ) stop

	write(6,*) 'file name: ',infile
	write(6,*) 'columns: ',nvar

	allocate(data(nvar))
	allocate(data_minmax(2,nvar))

c--------------------------------------------------------------
c read first record
c--------------------------------------------------------------

	call ts_read_next_record(iunit,nvar,dtime,data,datetime,ierr)
	if( ierr .ne. 0 ) goto 97

	if( datetime(1) .gt. 0 .and. .not. bquiet ) then
	  write(6,*) 'date and time: ',datetime
	end if

	call dts_convert_to_atime(datetime,dtime,atime)

c--------------------------------------------------------------
c close and re-open file
c--------------------------------------------------------------

	close(iunit)

	nvar = 0
	call ts_open_file(infile,nvar,datetime,iunit)
	if( iunit .le. 0 ) stop

c--------------------------------------------------------------
c loop on all records
c--------------------------------------------------------------

	nvar0 = nvar
	irec = 0
	idt = 0
	ich = 0
	atimeanf = atime
	atimeend = atime

	do 
	  irec = irec + 1
	  atimeold = atime

	  call ts_read_next_record(iunit,nvar,dtime,data,datetime,ierr)
	  if( ierr < 0 ) exit
	  if( ierr > 0 ) goto 97

	  if( nvar .ne. nvar0 ) goto 96

	  call dts_convert_to_atime(datetime,dtime,atime)
	  call dts_format_abs_time(atime,line)

	  if( bdebug ) write(6,*) irec,atime,line

	  if( bwrite ) then
            call minmax_ts(nvar,data,data_minmax)
	  end if

	  if( irec > 1 ) then
	    if( irec == 2 ) idt = nint(atime-atimeold)
	    idtact = nint(atime-atimeold)
	    if( idtact .ne. idt ) then
	      ich = ich + 1
	      write(6,*) '* change in time step: ',irec,idt,idtact
	      idt = idtact
	    end if
	    if( idt <= 0 ) then
	      write(6,*) '*** zero or negative time step: ',irec,idt
     +				,atime,atimeold
	    end if
	  end if
	  atimeend = atime
	end do

c--------------------------------------------------------------
c finish loop - info on time records
c--------------------------------------------------------------

	if( bwrite ) then
	  write(6,*) 'min/max of data: '
	  do i=1,nvar
	    write(6,*) i,data_minmax(1,i),data_minmax(2,i)
	  end do
	end if

	nrecs = irec - 1
	write(6,*) 'nrecs:  ',nrecs
	call dts_format_abs_time(atimeanf,line)
	write(6,*) 'start time: ',atimeanf,line
	call dts_format_abs_time(atimeend,line)
	write(6,*) 'end time:   ',atimeend,line
	write(6,*) 'idt:    ',idt

	if( ich .gt. 0 ) then
	  write(6,*) '* warning: time step changed: ',ich
	end if

	close(iunit)

	deallocate(data)
	deallocate(data_minmax)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	return
   96	continue
	write(6,*) 'nvar,nvar0: ',nvar,nvar0
	write(6,*) 'cannot change total number of columns'
	stop 'error stop tsinf'
   97	continue
	write(6,*) 'record: ',irec
	write(6,*) 'cannot read record of file'
	stop 'error stop tsinf'
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine minmax_ts(nvar,data,data_minmax)

        implicit none

        integer nvar
        real data(nvar)
        real data_minmax(2,nvar)

	integer i

        integer icall
	save icall
	data icall /0/

	if( icall == 0 ) then
	  icall = 1
          do i=1,nvar
	    data_minmax(1,i) = data(i)
	    data_minmax(2,i) = data(i)
	  end do
	end if

        do i=1,nvar
	  data_minmax(1,i) = min(data_minmax(1,i),data(i))
	  data_minmax(2,i) = max(data_minmax(2,i),data(i))
        end do

        end

c*****************************************************************

