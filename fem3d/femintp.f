!
! elaborates fem files
!
! revision log :
!
! 14.01.2015    ggu     adapted from feminf
! 20.05.2015    ggu     use bhuman to convert to human readable time
! 05.06.2015    ggu     iextract to extract nodal value
! 05.11.2015    ggu     new option chform to change format
! 04.10.2016    ggu     output flags now similar to shyelab
! 05.10.2016    ggu     allow for expansion of regular grid
! 11.10.2016    ggu     introduced flag for min/max/med computation
! 31.10.2016    ggu     new flag condense (bcondense)
! 16.05.2017    ggu&mbj better handling of points to extract
!
!******************************************************************

	program femelab

! interpolates between fem records

	use clo

	implicit none

	character*80 name,string,infile
	integer nfile
	double precision tmin,tmax
	character*80 stmin,stmax
	logical bdebug,bskip,bout,btmin,btmax
	logical bverb,bwrite,bquiet,binfo
	logical bchform,bcheckdt

	bdebug = .true.
	bdebug = .false.

!--------------------------------------------------------------
! set command line options
!--------------------------------------------------------------

	call clo_init('femelab','fem-file','1.2')

	call clo_add_info('elaborates and rewrites a fem file')

        call clo_add_sep('what to do (only one of these may be given)')

	call clo_add_option('out',.false.,'create output file out.fem')
        call clo_add_option('node node',' ','extract value for node')
        call clo_add_option('coord coord',' ','extract coordinate')
	call clo_add_option('split',.false.,'splits to single variables')
        call clo_add_option('regexpand iexp',-1,'expand regular grid')

        call clo_add_sep('options in/output')

        call clo_add_option('verb',.false.,'be more verbose')
        call clo_add_option('write',.false.,'write min/max of records')
        call clo_add_option('quiet',.false.,'do not write time records')
        call clo_add_option('info',.false.,'only give info on header')
	call clo_add_option('condense',.false.
     +				,'condense data to one node')

	call clo_add_sep('additional options')

	call clo_add_option('chform',.false.,'change output format')
        call clo_add_option('checkdt',.false.
     +                          ,'check for change of time step')
	call clo_add_option('tmin time',' '
     +				,'only process starting from time')
	call clo_add_option('tmax time',' '
     +				,'only process up to time')

	call clo_add_sep('additional information')

	call clo_add_extra('format for time is YYYY-MM-DD[::hh:mm:ss]')
	call clo_add_extra('time may be integer for relative time')
	call clo_add_extra('node is internal numbering in fem file')
	call clo_add_extra('   or ix,iy of regular grid')
	call clo_add_extra('coord is x,y of point to extract')

!--------------------------------------------------------------
! parse command line options
!--------------------------------------------------------------

	call clo_parse_options(1)  !expecting (at least) 1 file after options

!--------------------------------------------------------------
! get command line options
!--------------------------------------------------------------

	call clo_get_option('verb',bverb)
	call clo_get_option('write',bwrite)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('info',binfo)

	call clo_get_option('out',bout)
	call clo_get_option('chform',bchform)
        call clo_get_option('checkdt',bcheckdt)
	call clo_get_option('tmin',stmin)
	call clo_get_option('tmax',stmax)

!--------------------------------------------------------------
! set parameters
!--------------------------------------------------------------

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

!--------------------------------------------------------------
! loop on files
!--------------------------------------------------------------

	if( nfile > 1 ) then
	  write(6,*) 'Can only handle one file at a time'
	  stop 'error stop femelab: too many files'
	end if

        call clo_get_file(1,infile)
        if( infile .ne. ' ' ) call femintp_file(infile)

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

        end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine femintp_file(infile)

! writes info on fem file

	use clo
	use fem_util

	implicit none

	character*(*) infile

	character*80 name,string
	integer np,iunit,iout
	integer nvers,lmax,nvar,ntype,nlvdi
	integer nvar0,lmax0,np0
	integer idt,idtact
	double precision dtime,atmin,atmax,atime0,atime1997
	double precision atime,atime1,atime2
	double precision rit
	real dmin,dmax,dmed
	integer ierr
	integer nfile,nintp
	integer irec,iv,ich,isk,nrecs,iu88,l,i,ivar,nout
	integer itype(2)
	integer iformat,iformout
	integer datetime(2),dateanf(2),dateend(2)
	integer iextract,it
	integer next,idte
	integer ie,nx,ny,ix,iy
	integer regexpand
	integer np_out
	real regpar0(7)
	real regpar(7)
	real flag
	real xp,yp
	logical bdebug,bfirst,bskip,bout,btmin,btmax,boutput
	logical bhuman,blayer,bcondense
	logical bverb,bwrite,bquiet,binfo
	logical bchform,bcheckdt,bdtok,bextract,breg,bintime,bexpand
	logical bsplit,bread
	integer, allocatable :: ivars(:)
	character*80, allocatable :: strings(:)
	character*20 dline,aline,fline
	character*40 eformat
	character*80 stmin,stmax
	character*80 snode,scoord
	real,allocatable :: data(:,:,:)
	real,allocatable :: data_profile(:)
	real,allocatable :: dext(:)
	real,allocatable :: d3dext(:,:)
	real,allocatable :: hd(:)
	real,allocatable :: hlv(:)
	integer,allocatable :: ilhkv(:)
	integer,allocatable :: ius(:)

        type(femfile_type) :: ffinfo_in,ffinfo_out
        type(fem_type) :: finfo1,finfo2,finfo

	integer ifileo

	bdebug = .true.
	bdebug = .false.

	call clo_get_option('verb',bverb)
	call clo_get_option('write',bwrite)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('info',binfo)
	call clo_get_option('condense',bcondense)

	call clo_get_option('out',bout)
	call clo_get_option('regexpand',regexpand)
        call clo_get_option('split',bsplit)
	call clo_get_option('chform',bchform)
        call clo_get_option('checkdt',bcheckdt)
        call clo_get_option('node',snode)
        call clo_get_option('coord',scoord)
	call clo_get_option('tmin',stmin)
	call clo_get_option('tmax',stmax)

	if( bchform ) bout = .true.
	if( bcondense ) bout = .true.
	bexpand = regexpand > -1
	boutput = .true.

	nintp = 4	!how many records to interpolate (-1)
	next = 2	!extend records before and after read records
	idte = 6*3600	!time step for records to extend

!--------------------------------------------------------------
! open file
!--------------------------------------------------------------

	if( infile .eq. ' ' ) stop

	call femutil_open_for_read(infile,0,ffinfo_in,ierr)
	if( ierr /= 0 ) stop

	write(6,*) 'file name: ',trim(infile)

!--------------------------------------------------------------
! prepare for output if needed
!--------------------------------------------------------------

	call femutil_open_for_write('out.fem'
     +			,ffinfo_in%iformat,ffinfo_out)

!--------------------------------------------------------------
! read first record
!--------------------------------------------------------------

	call femutil_init_record(finfo)
	call femutil_init_record(finfo1)
	call femutil_init_record(finfo2)

	irec = 0
	irec = irec + 1
	call femutil_read_record(ffinfo_in,finfo1,ierr)
	if( ierr .ne. 0 ) goto 99

!--------------------------------------------------------------
! write info to terminal
!--------------------------------------------------------------

	call femutil_get_time(finfo1,atime1)
	call dts_format_abs_time(atime1,dline)
	write(6,*) trim(dline)

!--------------------------------------------------------------
! write first records
!--------------------------------------------------------------

	nout = 0
	if( boutput ) then
	  finfo2 = finfo1
	  atime = atime1 - (next+1)*idte
	  do it=1,next			!extend records
	      atime = atime + idte
	      call femutil_set_time(finfo2,atime)
	      call dts_format_abs_time(atime,dline)
	      write(6,*) '   extp: ',trim(dline)
	      nout = nout + 1
	      call femutil_write_record(ffinfo_out,finfo2)
	  end do
	  nout = nout + 1
	  call femutil_write_record(ffinfo_out,finfo1)
	end if


!--------------------------------------------------------------
! loop on all records
!--------------------------------------------------------------

	do 
	  irec = irec + 1
	  if( irec > 10 ) exit
	  call femutil_read_record(ffinfo_in,finfo2,ierr)
	  if( ierr .lt. 0 ) exit
	  if( ierr .gt. 0 ) goto 99

	  call femutil_get_time(finfo2,atime2)
	  call dts_format_abs_time(atime2,dline)
	  write(6,*) trim(dline)

	  if( .not. femutil_is_compatible(finfo1,finfo2) ) goto 92

          if( boutput ) then
	    do it=1,nintp-1
	      nout = nout + 1
	      rit = float(it)/nintp
	      call record_interpolate(finfo1,finfo2,finfo,rit)
	      call femutil_get_time(finfo,atime)
	      call dts_format_abs_time(atime,dline)
	      write(6,*) '   intp: ',trim(dline)
	
	      call femutil_write_record(ffinfo_out,finfo)
	    end do
	    nout = nout + 1
	    call femutil_write_record(ffinfo_out,finfo2)
          end if

	  atime1 = atime2
	  finfo1 = finfo2
	end do

!--------------------------------------------------------------
! finish loop - last records
!--------------------------------------------------------------

	if( boutput ) then
	  finfo2 = finfo1
	  do it=1,next			!extend records
	      atime = atime + idte
	      call femutil_set_time(finfo2,atime)
	      call dts_format_abs_time(atime,dline)
	      write(6,*) '   extp: ',trim(dline)
	      call femutil_write_record(ffinfo_out,finfo2)
	  end do
	end if

!--------------------------------------------------------------
! info on time records
!--------------------------------------------------------------

	nrecs = irec - 1
	write(6,*) 'nrecs:  ',nrecs
	if( nout > 0 ) then
	  write(6,*) 'output records written: ',nout
	  write(6,*) 'output written to file out.fem'
	end if

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	return
   91	continue
	write(6,*) 'iectract,np: ',iextract,np
	stop 'error stop femelab: no such node'
   92	continue
	write(6,*) 'records are not compatible...'
	stop 'error stop femelab: compatibility'
   95	continue
	write(6,*) 'strings not in same sequence: ',iv
        write(6,*) string
        write(6,*) strings(iv)
	stop 'error stop femelab: strings'
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

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine record_interpolate(finfo1,finfo2,finfo,rit)

	use fem_util

	implicit none

        type(fem_type) :: finfo1
        type(fem_type) :: finfo2
        type(fem_type) :: finfo
	double precision rit

	integer iv,ip,lmax,l
	real val1,val2,val
	double precision atime1,atime2,atime

	finfo = finfo1

	call femutil_get_time(finfo1,atime1)
	call femutil_get_time(finfo2,atime2)
	atime = atime1 + rit*(atime2-atime1)
	call femutil_set_time(finfo,atime)

	do iv=1,finfo%nvar
	  do ip=1,finfo%np
	    lmax = finfo%ilhkv(ip)
	    do l=1,lmax
	      val1 = finfo1%data(l,ip,iv)
	      val2 = finfo2%data(l,ip,iv)
	      val = val1 + rit*(val2-val1)
	      finfo%data(l,ip,iv) = val
	    end do
	  end do
	end do

	end

!*****************************************************************

