!
! elaborates fem files
!
! revision log :
!
! 14.01.2015    ggu     adapted from feminf
! 20.05.2015    ggu     use bhuman to convert to human readable time
! 05.06.2015    ggu     iextract to extract nodal value
! 05.11.2015    ggu     new option chform to change format
!
!******************************************************************

	program femelab

c writes info on fem file

	use clo

	implicit none

	character*80 name,string,infile
	integer nfile
	double precision tmin,tmax
	character*80 stmin,stmax
	logical bdebug,bskip,bwrite,bout,btmin,btmax,bquiet
	logical bchform

	bdebug = .true.
	bdebug = .false.

c--------------------------------------------------------------
c parameters and command line options
c--------------------------------------------------------------

	call clo_init('femelab','fem-file','1.2')

	call clo_add_info('elaborates and rewrites a fem file')
	call clo_add_option('write',.false.,'write min/max of values')
	call clo_add_option('out',.false.,'create output file out.fem')
	call clo_add_option('chform',.false.,'change output format')
	call clo_add_option('quiet',.false.,'do not be verbose')
	call clo_add_option('tmin time',' '
     +				,'only process starting from time')
	call clo_add_option('tmax time',' '
     +				,'only process up to time')
	call clo_add_extra('format for time is YYYY-MM-DD[::hh:mm:ss]')
	call clo_add_extra('time may be integer for relative time')

	call clo_parse_options(1)  !expecting (at least) 1 file after options

	call clo_get_option('write',bwrite)
	call clo_get_option('out',bout)
	call clo_get_option('chform',bchform)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('tmin',stmin)
	call clo_get_option('tmax',stmax)

	bskip = .not. bwrite
	if( bout ) bskip = .false.
	btmin = tmin .ne. -1.
	btmax = tmax .ne. -1.

	nfile = clo_number_of_files()

	if( bdebug ) then
	  write(6,*) nfile
	  write(6,*) bwrite,bskip,bout,btmin,btmax
	  write(6,*) tmin,tmax
	end if

c--------------------------------------------------------------
c loop on files
c--------------------------------------------------------------

	if( nfile > 1 ) then
	  write(6,*) 'Can only handle one file at a time'
	  stop 'error stop femelab: too many files'
	end if

        call clo_get_file(1,infile)
        if( infile .ne. ' ' ) call femelab_file(infile)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine femelab_file(infile)

c writes info on fem file

	use clo

	implicit none

	character*(*) infile

	character*80 name,string
	integer np,iunit,iout
	integer nvers,lmax,nvar,ntype
	integer nvar0,lmax0,np0
	integer idt,idtact
	double precision dtime,atmin,atmax,atime0,atime1997
	double precision atime,atimeold,atimeanf,atimeend
	real dmin,dmax
	integer ierr
	integer nfile
	integer irec,i,ich,nrecs
	integer itype(2)
	integer iformat,iformout
	integer datetime(2),dateanf(2),dateend(2)
	integer iextract,it
	real regpar(7)
	logical bdebug,bfirst,bskip,bwrite,bout,btmin,btmax,boutput
	logical bquiet,bhuman
	logical bchform
	character*80, allocatable :: strings(:)
	character*20 line
	character*80 stmin,stmax
	real,allocatable :: data(:,:,:)
	real,allocatable :: dext(:)
	real,allocatable :: hd(:)
	real,allocatable :: hlv(:)
	integer,allocatable :: ilhkv(:)

	bdebug = .true.
	bdebug = .false.
	bhuman = .true.

	iextract = 5
	iextract = 0
	datetime = 0
	datetime(1) = 19970101
	dtime = 0.
	call fem_file_convert_time(datetime,dtime,atime)
	atime1997 = atime

        datetime = 0
        irec = 0

	call clo_get_option('write',bwrite)
	call clo_get_option('out',bout)
	call clo_get_option('chform',bchform)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('tmin',stmin)
	call clo_get_option('tmax',stmax)

	bskip = .not. bwrite
	if( bchform ) bout = .true.
	if( bout ) bskip = .false.

	atmin = 0.
	atmax = 0.
	btmin = stmin .ne. ' '
	btmax = stmax .ne. ' '
	if( btmin ) call fem_file_string2time(stmin,atmin)
	if( btmax ) call fem_file_string2time(stmax,atmax)

	!write(6,*) stmin(1:len_trim(stmin)),btmin,atmin
	!write(6,*) stmax(1:len_trim(stmax)),btmax,atmax
	
c--------------------------------------------------------------
c open file
c--------------------------------------------------------------

	if( infile .eq. ' ' ) stop

	np = 0
	call fem_file_read_open(infile,np,iunit,iformat)
	if( iunit .le. 0 ) stop

	write(6,*) 'file name: ',infile(1:len_trim(infile))
	call fem_file_get_format_description(iformat,line)
	write(6,*) 'format: ',iformat,"  (",line(1:len_trim(line)),")"

c--------------------------------------------------------------
c prepare for output if needed
c--------------------------------------------------------------

        iout = 0
	iformout = iformat
	if( bchform ) iformout = 1 - iformat
	if( iformout < 0 ) iformout = iformat

        if( bout ) then
          iout = iunit + 1
          if( iformout .eq. 1 ) then
	    open(iout,file='out.fem',status='unknown',form='formatted')
          else
	    open(iout,file='out.fem',status='unknown',form='unformatted')
          end if
        end if

c--------------------------------------------------------------
c read first record
c--------------------------------------------------------------

        call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)

	if( ierr .ne. 0 ) goto 99

	if( .not. bquiet ) then
	  write(6,*) 'nvers:  ',nvers
	  write(6,*) 'np:     ',np
	  write(6,*) 'lmax:   ',lmax
	  write(6,*) 'nvar:   ',nvar
	  write(6,*) 'ntype:  ',ntype
	end if

	allocate(hlv(lmax))
	call fem_file_make_type(ntype,2,itype)

	call fem_file_read_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar,ierr)
	if( ierr .ne. 0 ) goto 98

	if( lmax > 1 .and. .not. bquiet ) then
	  write(6,*) 'vertical layers: ',lmax
	  write(6,*) hlv
	end if
	if( itype(1) .gt. 0 .and. .not. bquiet ) then
	  write(6,*) 'date and time: ',datetime
	end if
	if( itype(2) .gt. 0 .and. .not. bquiet ) then
	  write(6,*) 'regpar: ',regpar
	end if

	call fem_file_convert_time(datetime,dtime,atime)
	atime0 = atime		!absolute time of first record

	if( iextract > 0 ) then
	  bskip = .false.
	  write(88,'(a,2i10)') '#date: ',datetime
	end if

	nvar0 = nvar
	lmax0 = lmax
	np0 = np
	allocate(strings(nvar))
	allocate(dext(nvar))
	allocate(data(lmax,np,nvar))
	allocate(hd(np))
	allocate(ilhkv(np))

	do i=1,nvar
	  call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
	  if( ierr .ne. 0 ) goto 97
	  if( .not. bquiet ) write(6,*) 'data:   ',i,'  ',trim(string)
	  strings(i) = string
	end do

c--------------------------------------------------------------
c close and re-open file
c--------------------------------------------------------------

	close(iunit)

	np = 0
	call fem_file_read_open(infile,np,iunit,iformat)
	if( iunit .le. 0 ) stop

c--------------------------------------------------------------
c loop on all records
c--------------------------------------------------------------

	irec = 0
	idt = 0
	ich = 0
	atimeanf = atime
	atimeend = atime
	atimeold = atime - 1

	do 
	  irec = irec + 1
          call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	  if( ierr .lt. 0 ) exit
	  if( ierr .gt. 0 ) goto 99
	  if( nvar .ne. nvar0 ) goto 96
	  if( lmax .ne. lmax0 ) goto 96
	  if( np .ne. np0 ) goto 96

	  call fem_file_convert_time(datetime,dtime,atime)
	  call dts_format_abs_time(atime,line)

	  if( bdebug ) write(6,*) irec,atime,line

	  call fem_file_read_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar,ierr)
	  if( ierr .ne. 0 ) goto 98

          boutput = bout .and. atime > atimeold
	  if( btmin ) boutput = boutput .and. atime >= atmin
	  if( btmax ) boutput = boutput .and. atime <= atmax

          if( boutput ) then
	    if( bhuman ) then
	      call fem_file_convert_atime(datetime,dtime,atime)
	    end if
            call fem_file_write_header(iformout,iout,dtime
     +                          ,nvers,np,lmax,nvar,ntype,lmax
     +                          ,hlv,datetime,regpar)
          end if

	  do i=1,nvar
	    if( bskip ) then
	      call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
	    else
              call fem_file_read_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,lmax,data(1,1,i)
     +                          ,ierr)
	    end if
	    if( ierr .ne. 0 ) goto 97
	    if( string .ne. strings(i) ) goto 95
            if( boutput ) then
              call fem_file_write_data(iformout,iout
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,lmax,data(1,1,i))
            end if
	    if( bwrite ) then
              call minmax_data(lmax,np,ilhkv,data(1,1,i),dmin,dmax)
	      write(6,1100) irec,i,atime,dmin,dmax,line
 1100	      format(i6,i3,f15.2,2g16.5,1x,a20)
	    end if
	    if( iextract > 0 ) then
	      dext(i) = data(1,iextract,i)
	    end if
	  end do

	  if( iextract > 0 ) then
	    it = nint(atime-atime0)
	    it = nint(atime-atime1997)
	    write(88,*) it,dext
	  end if

	  if( atime <= atimeold ) then
	    ich = ich + 1
	    write(6,*) '*** time of record is before last one... skipping'
	  else
	    atimeold = atime
	  end if
	  atimeend = atime
	end do

c--------------------------------------------------------------
c finish loop - info on time records
c--------------------------------------------------------------

	nrecs = irec - 1
	write(6,*) 'nrecs:  ',nrecs
	call dts_format_abs_time(atimeanf,line)
	write(6,*) 'start time: ',atimeanf,line
	call dts_format_abs_time(atimeend,line)
	write(6,*) 'end time:   ',atimeend,line

	if( ich .gt. 0 ) then
	  write(6,*) '* warning: records eliminated: ',ich
	end if

	close(iunit)
	if( iout > 0 ) close(iout)

	if( bout ) then
	  write(6,*) 'output written to file out.fem'
	end if

	if( iextract > 0 ) then
	  write(6,*) '********** warning **************'
	  write(6,*) 'iextract = ',iextract
	  write(6,*) 'data written to fort.88'
	  write(6,*) '********** warning **************'
	end if
c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	return
   95	continue
	write(6,*) 'variable ',i
	write(6,*) trim(string)
	write(6,*) trim(strings(i))
	write(6,*) 'cannot change description of variables'
	stop 'error stop femelab'
   96	continue
	write(6,*) 'nvar,nvar0: ',nvar,nvar0
	write(6,*) 'lmax,lmax0: ',lmax,lmax0	!this might be relaxed
	write(6,*) 'np,np0:     ',np,np0	!this might be relaxed
	write(6,*) 'cannot change number of variables'
	stop 'error stop femelab'
   97	continue
	write(6,*) 'record: ',irec
	write(6,*) 'cannot read data record of file'
	stop 'error stop femelab'
   98	continue
	write(6,*) 'record: ',irec
	write(6,*) 'cannot read second header of file'
	stop 'error stop femelab'
   99	continue
	write(6,*) 'record: ',irec
	write(6,*) 'cannot read header of file'
	stop 'error stop femelab'
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine minmax_data(nlvddi,np,ilhkv,data,vmin,vmax)

        implicit none

        integer nlvddi,np
        integer ilhkv(1)
        real data(nlvddi,1)
	real vmin,vmax

        integer k,l,lmax
        real v

        vmin = data(1,1)
        vmax = data(1,1)

        do k=1,np
          lmax = ilhkv(k)
          do l=1,lmax
            v = data(l,k)
            vmax = max(vmax,v)
            vmin = min(vmin,v)
          end do
        end do

        !write(86,*) 'min/max: ',it,vmin,vmax

        end

c*****************************************************************

