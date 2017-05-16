
	program tsinf

c writes info on ts file

	use clo

	implicit none

	character*50 infile
	integer nfile
	integer i
	logical bdebug,bwrite,bout,btmin,btmax
	logical bquiet,bconvert,bcheck
	double precision tmin,tmax

	bdebug = .true.
	bdebug = .false.

c--------------------------------------------------------------
c parameters and command line options
c--------------------------------------------------------------

	call clo_init('tsinf','ts-file(s)','1.1')

	call clo_add_option('write',.false.,'write min/max of values')
	call clo_add_option('out',.false.,'create output file')
	call clo_add_option('convert',.false.,'convert to ISO time')
	call clo_add_option('check',.false.,'check regular time step')
	call clo_add_option('quiet',.false.,'do not be verbose')
	call clo_add_option('tmin time',-1
     +				,'only process starting from time')
	call clo_add_option('tmax time',-1
     +				,'only process up to time')

	call clo_parse_options(1)  !expecting (at least) 1 file after options

	call clo_get_option('write',bwrite)
	call clo_get_option('out',bout)
	call clo_get_option('convert',bconvert)
	call clo_get_option('check',bcheck)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('tmin',tmin)
	call clo_get_option('tmax',tmax)

	btmin = tmin .ne. -1.
	btmax = tmax .ne. -1.

	if( bconvert ) bout = .true.
	!if( bout ) stop 'error stop: -out not yet implemented'
	if( btmin ) stop 'error stop: -tmin not yet implemented'
	if( btmax ) stop 'error stop: -tmax not yet implemented'

	nfile = clo_number_of_files()

	if( bdebug ) then
	  write(6,*) nfile
	  write(6,*) bwrite,bout,btmin,btmax,bconvert,bcheck
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
	integer np,iunit
	integer nvers,lmax,nvar,ntype
	integer lmax0,np0
	integer idt,idtact
	double precision dtime,tmin,tmax,dtime0
	double precision atime,atimeold,atimeanf,atimeend
	double precision atime0,atime0e
	real dmin,dmax
	integer ierr
	integer nfile
	integer irec,i,ich,nrecs
	integer datetime(2)
	logical bdebug,bfirst,bskip,bwrite,bout,btmin,btmax,boutput
	logical bquiet,bconvert,bcheck
	character*20 line
	character*20 format
	real,allocatable :: data(:)
	real,allocatable :: data_minmax(:,:)
	integer, save :: iout = 0
	integer, save :: nvar0 = 0

	bdebug = .true.
	bdebug = .false.

	datetime = 0
	irec = 0

	call clo_get_option('write',bwrite)
	call clo_get_option('out',bout)
	call clo_get_option('convert',bconvert)
	call clo_get_option('check',bcheck)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('tmin',tmin)
	call clo_get_option('tmax',tmax)

	if( bconvert ) bout = .true.

c--------------------------------------------------------------
c open file
c--------------------------------------------------------------

	if( infile .eq. ' ' ) stop

	atime0e = 0.
	if( bout ) then		!see if we get extra information on time
	  call ts_get_extra_time(infile,dtime,datetime)
	  if( datetime(1) > 0 ) then
	    dtime0 = 0.
	    call dts_convert_to_atime(datetime,dtime0,atime)
	    atime0e = atime - dtime
	  end if
	end if

	nvar = 0
	call ts_open_file(infile,nvar,datetime,iunit)
	if( iunit .le. 0 ) stop

	if( .not. bquiet ) then
	  write(6,*) 'file name: ',infile
	  write(6,*) 'columns: ',nvar
	  if( datetime(1) > 0 ) write(6,*) 'datetime: ',datetime
	end if

	if( nvar0 == 0 ) nvar0 = nvar
	if( nvar /= nvar0 ) then
	  write(6,*) 'nvar,nvar0: ',nvar,nvar0
	  stop 'error stop: nvar /= nvar0'
	end if

	allocate(data(nvar))
	allocate(data_minmax(2,nvar))

	if( bout .and. iout == 0 ) then
	  iout = 1
	  open(iout,file='out.txt',form='formatted',status='unknown')
	  !write(format,'(a,i3,a)') '(f18.2,',nvar,'g14.6)'
	  write(format,'(a,i3,a)') '(a20,',nvar,'g14.6)'
	  write(6,*) 'used format: ',trim(format)
	end if

	if( bout ) then
	  if( datetime(1) /= 0 ) then
	    dtime = 0.
	    call dts_convert_to_atime(datetime,dtime,atime)
	    atime0 = atime
	    write(6,*) 'using absolute date from file'
	  else if( atime0e /= 0 ) then
	    atime0 = atime0e
	    write(6,*) 'using absolute date from extra information'
	  else
	    write(6,*) 'no absolute time... cannot convert'
	    stop 'error stop: missing absolute time'
	  end if
	  call dts_format_abs_time(atime0,line)
	  write(6,*) 'absolute date reference: ',line
	end if

c--------------------------------------------------------------
c read first record
c--------------------------------------------------------------

	call ts_read_next_record(iunit,nvar,dtime,data,datetime,ierr)
	if( ierr .ne. 0 ) goto 97

	call dts_convert_to_atime(datetime,dtime,atime)

c--------------------------------------------------------------
c close and re-open file
c--------------------------------------------------------------

	close(iunit)

	nvar = 0
	call ts_open_file(infile,nvar,datetime,iunit)
	if( iunit .le. 0 ) stop
	!write(6,*) datetime

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

	  !write(6,*) 'ggguuu: ',datetime,dtime
	  call dts_convert_to_atime(datetime,dtime,atime)
	  call dts_format_abs_time(atime,line)

	  if( bdebug ) write(6,*) irec,atime,line

	  if( bwrite ) then
            call minmax_ts(nvar,data,data_minmax)
	    if( .not. bquiet ) write(6,*) irec,atime,line
	  end if

	  bskip = .false.
	  if( irec > 1 ) then
	    if( irec == 2 ) idt = nint(atime-atimeold)
	    idtact = nint(atime-atimeold)
	    if( idtact .ne. idt ) then
	      ich = ich + 1
	      if( bcheck ) then
	        write(6,*) '* change in time step: ',irec,idt,idtact
	      end if
	      idt = idtact
	    end if
	    bskip = idt <= 0
	    if( bcheck .and. bskip ) then
	      write(6,*) '*** zero or negative time step: ',irec,idt
     +				,atime,atimeold
	    end if
	  end if

	  if( bout .and. .not. bskip ) then
	    atime = atime0 + dtime
	    call dts_format_abs_time(atime,line)
	    write(1,format) line,data(1:nvar)
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

