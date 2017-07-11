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

c writes info on fem file

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

c--------------------------------------------------------------
c set command line options
c--------------------------------------------------------------

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

c--------------------------------------------------------------
c parse command line options
c--------------------------------------------------------------

	call clo_parse_options(1)  !expecting (at least) 1 file after options

c--------------------------------------------------------------
c get command line options
c--------------------------------------------------------------

	call clo_get_option('verb',bverb)
	call clo_get_option('write',bwrite)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('info',binfo)

	call clo_get_option('out',bout)
	call clo_get_option('chform',bchform)
        call clo_get_option('checkdt',bcheckdt)
	call clo_get_option('tmin',stmin)
	call clo_get_option('tmax',stmax)

c--------------------------------------------------------------
c set parameters
c--------------------------------------------------------------

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
	integer nvers,lmax,nvar,ntype,nlvdi
	integer nvar0,lmax0,np0
	integer idt,idtact
	double precision dtime,atmin,atmax,atime0,atime1997
	double precision atime,atimeold,atimeanf,atimeend
	real dmin,dmax,dmed
	integer ierr
	integer nfile
	integer irec,iv,ich,isk,nrecs,iu88,l,i,ivar
	integer itype(2)
	integer iformat,iformout
	integer datetime(2),dateanf(2),dateend(2)
	integer iextract,it
	integer ie,nx,ny,ix,iy
	integer regexpand
	integer np_out
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

	integer ifileo

	bdebug = .true.
	bdebug = .false.
	bhuman = .true.		!convert time in written fem file to dtime=0
	blayer = .false.
	blayer = .true.		!write layer structure - should be given by CLO

	iextract = 0
	iu88 = 0
	datetime = 0
	datetime(1) = 19970101
	dtime = 0.
	call dts_convert_to_atime(datetime,dtime,atime)
	atime1997 = atime

        datetime = 0
        irec = 0
	regpar = 0.

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

	bextract = ( snode /= ' ' .or. scoord /= ' ' )

	atmin = 0.
	atmax = 0.
	btmin = stmin .ne. ' '
	btmax = stmax .ne. ' '
	if( btmin ) call dts_string2time(stmin,atmin)
	if( btmax ) call dts_string2time(stmax,atmax)

	!write(6,*) stmin(1:len_trim(stmin)),btmin,atmin
	!write(6,*) stmax(1:len_trim(stmax)),btmax,atmax
	
c--------------------------------------------------------------
c open file
c--------------------------------------------------------------

	if( infile .eq. ' ' ) stop

	np = 0
	call fem_file_read_open(infile,np,iformat,iunit)
	if( iunit .le. 0 ) stop

	write(6,*) 'file name: ',infile(1:len_trim(infile))
	call fem_file_get_format_description(iformat,fline)
	write(6,*) 'format: ',iformat,"  (",trim(fline),")"

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

	write(6,*) 'nvers:  ',nvers
	write(6,*) 'np:     ',np
	write(6,*) 'lmax:   ',lmax
	write(6,*) 'nvar:   ',nvar
	write(6,*) 'ntype:  ',ntype

	allocate(hlv(lmax))
	call fem_file_make_type(ntype,2,itype)

	call fem_file_read_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar,ierr)
	if( ierr .ne. 0 ) goto 98

	if( lmax > 1 ) then
	  write(6,*) 'levels: ',lmax
	  write(6,'(5f12.4)') hlv
	end if
	if( itype(1) .gt. 0 ) then
	  write(6,*) 'date and time: ',datetime
	end if
	breg = .false.
	if( itype(2) .gt. 0 ) then
	  breg = .true.
	  write(6,*) 'regpar: '
	  !write(6,'(4f12.4)') regpar
	  write(6,'(4x,a,2i12)') 'nx,ny: ',nint(regpar(1)),nint(regpar(2))
	  write(6,'(4x,a,2f12.4)') 'x0,y0: ',regpar(3),regpar(4)
	  write(6,'(4x,a,2f12.4)') 'dx,dy: ',regpar(5),regpar(6)
	  write(6,'(4x,a,2f12.4)') 'flag : ',regpar(7)
	end if

	call dts_convert_to_atime(datetime,dtime,atime)
	atime0 = atime		!absolute time of first record

	if( bextract ) then
	  bskip = .false.
	  call handle_extract(breg,np,regpar,iextract)
	  iu88 = ifileo(88,'out.txt','form','new')
	  write(iu88,'(a,2i10)') '#date: ',datetime
	  write(iu88,'(a,2f12.5)') '#coords: ',xp,yp
	  write(eformat,'(a,i3,a)') '(i12,',nvar,'g14.6,a2,a20)'
	  write(6,*) 'using format: ',trim(eformat)
	end if

	nvar0 = nvar
	lmax0 = lmax
	nlvdi = lmax
	np0 = np
	allocate(ivars(nvar))
	allocate(strings(nvar))
	allocate(dext(nvar))
	allocate(d3dext(nlvdi,nvar))
	allocate(data(nlvdi,np,nvar))
	allocate(data_profile(nlvdi))
	allocate(hd(np))
	allocate(ilhkv(np))
	allocate(ius(nvar))
	ius = 0

c--------------------------------------------------------------
c write info to terminal
c--------------------------------------------------------------

        write(6,*) 'available variables contained in file: '
        write(6,*) 'total number of variables: ',nvar
        write(6,*) '   varnum     varid    varname'

	do iv=1,nvar
	  call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
	  if( ierr .ne. 0 ) goto 97
	  call string2ivar_n(string,ivar)
	  write(6,'(2i10,4x,a)') iv,ivar,trim(string)
	  ivars(iv) = ivar
	  strings(iv) = string
	end do
	if( binfo ) return

c--------------------------------------------------------------
c close and re-open file
c--------------------------------------------------------------

	close(iunit)

	np = 0
	call fem_file_read_open(infile,np,iformat,iunit)
	if( iunit .le. 0 ) stop

c--------------------------------------------------------------
c loop on all records
c--------------------------------------------------------------

	irec = 0
	idt = 0
	ich = 0
	isk = 0
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

	  call dts_convert_to_atime(datetime,dtime,atime)
	  call dts_format_abs_time(atime,dline)

	  if( bdebug ) write(6,*) irec,atime,dline
	  if( .not. bquiet ) then
            write(6,'(a,2f20.2,3x,a20)') 'time : ',atime,dtime,dline
	  end if

	  call fem_file_read_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar,ierr)
	  if( ierr .ne. 0 ) goto 98

	  flag = regpar(7)
	  call fem_file_make_type(ntype,2,itype)
	  breg = ( itype(2) .gt. 0 )

	  bdtok = atime > atimeold
          boutput = bout .and. bdtok
	  bread = bwrite .or. bextract .or. boutput
	  bread = bread .or. bsplit
	  bskip = .not. bread
	  bintime = .true.
	  if( btmin ) bintime = bintime .and. atime >= atmin
	  if( btmax ) bintime = bintime .and. atime <= atmax
	  boutput = boutput .and. bintime
	  if( .not. bintime ) bskip = .true.

          if( boutput ) then
	    if( bhuman ) then
	      call dts_convert_from_atime(datetime,dtime,atime)
	    end if
	    np_out = np
	    if( bcondense ) np_out = 1
            call fem_file_write_header(iformout,iout,dtime
     +                          ,0,np_out,lmax,nvar,ntype,lmax
     +                          ,hlv,datetime,regpar)
          end if

	  do iv=1,nvar
	    if( bskip ) then
	      call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
	    else
              call fem_file_read_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,data(1,1,iv)
     +                          ,ierr)
	    end if
	    if( ierr .ne. 0 ) goto 97
	    if( string .ne. strings(iv) ) goto 95
	    !write(6,*) iv,'  ',trim(string)
            if( boutput ) then
	      !call custom_elab(nlvdi,np,string,iv,data(1,1,iv))
	      if( breg .and. bexpand ) then
		call reg_set_flag(nlvdi,np,ilhkv,regpar,data(1,1,iv))
		call reg_expand_shell(nlvdi,np,lmax,regexpand
     +					,regpar,ilhkv,data(1,1,iv))
	      end if
	      if( bcondense ) then
		call fem_condense(np,lmax,data(1,1,iv),flag,data_profile)
                call fem_file_write_data(iformout,iout
     +                          ,0,np_out,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,data_profile)
	      else
                call fem_file_write_data(iformout,iout
     +                          ,0,np_out,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvdi,data(1,1,iv))
	      end if
            end if
	    if( bwrite .and. .not. bskip ) then
	      write(6,*) irec,iv,ivars(iv),trim(strings(iv))
	      do l=1,lmax
                call minmax_data(l,lmax,np,flag,ilhkv,data(1,1,iv)
     +					,dmin,dmax,dmed)
	        write(6,1000) 'l,min,aver,max : ',l,dmin,dmed,dmax
 1000	        format(a,i5,3g16.6)
	      end do
	    end if
	    if( bextract ) then
	      dext(iv) = data(1,iextract,iv)
	      d3dext(:,iv) = data(:,iextract,iv)
	    end if
	    if( bsplit ) then
	      call femsplit(iformout,ius(iv),dtime,nvers,np
     +			,lmax,nlvdi,ntype
     +			,hlv,datetime,regpar,string
     +			,ilhkv,hd,data(:,:,iv))
	    end if
	  end do

	  if( bextract .and. bdtok .and. bintime ) then
	    it = nint(atime-atime0)
	    it = nint(atime-atime1997)
	    call dts_format_abs_time(atime,aline)
	    write(iu88,eformat) it,dext,'  ',aline
	  end if

	  call check_dt(atime,atimeold,bcheckdt,irec,idt,ich,isk)
	  atimeold = atime
	  if( .not. bdtok ) cycle

	  atimeend = atime
	end do

c--------------------------------------------------------------
c finish loop - info on time records
c--------------------------------------------------------------

	nrecs = irec - 1
	write(6,*) 'nrecs:  ',nrecs
	call dts_format_abs_time(atimeanf,aline)
	write(6,*) 'start time: ',atimeanf,aline
	call dts_format_abs_time(atimeend,aline)
	write(6,*) 'end time:   ',atimeend,aline

        if( ich == 0 ) then
          write(6,*) 'idt:    ',idt
        else
          write(6,*) 'idt:     irregular ',ich,isk
        end if

	if( isk .gt. 0 ) then
	  write(6,*) '*** warning: records eliminated: ',isk
	end if

	close(iunit)
	if( iout > 0 ) close(iout)

	if( bout ) then
	  write(6,*) 'output written to file out.fem'
	end if

	if( bextract ) then
	  write(6,*) 'iextract = ',iextract
	  write(6,*) 'data written to out.txt'
	end if

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	return
   91	continue
	write(6,*) 'iectract,np: ',iextract,np
	stop 'error stop femelab: no such node'
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

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine minmax_data(level,nlvddi,np,flag,ilhkv,data
     +				,vmin,vmax,vmed)

        implicit none

	integer level		!level for which minmax to compute (0 for all)
        integer nlvddi,np
	real flag
        integer ilhkv(1)
        real data(nlvddi,1)
	real vmin,vmax,vmed

        integer k,l,lmin,lmax,lm,ntot
        real v,high
	double precision vtot

	lmin = max(1,level)
	lmax = level
	if( level == 0 ) lmax = nlvddi

	ntot = 0
	vtot = 0.
	high = 1.e+30
        vmin = high
        vmax = -high

        do k=1,np
          lm = min(ilhkv(k),lmax)
          do l=lmin,lm
            v = data(l,k)
	    if( v == flag ) cycle
	    ntot = ntot + 1
	    vtot = vtot + v
            vmax = max(vmax,v)
            vmin = min(vmin,v)
          end do
        end do

	if( ntot > 0 ) then
	  vmed = vtot / ntot
	else
	  vmin = 0.
	  vmax = 0.
	  vmed = 0.
	end if

        !write(86,*) 'min/max: ',it,vmin,vmax,vmed

        end

c*****************************************************************

	subroutine check_dt(atime,atimeold,bcheckdt,irec,idt,ich,isk)

	implicit none

	double precision atime,atimeold
	logical bcheckdt
	integer irec,idt,ich,isk

	integer idtact
	character*20 aline

          if( irec > 1 ) then
            if( irec == 2 ) idt = nint(atime-atimeold)
            idtact = nint(atime-atimeold)
            if( idtact .ne. idt ) then
              ich = ich + 1
              if( bcheckdt ) then
		call dts_format_abs_time(atime,aline)
                write(6,'(a,3i10,a,a)') '* change in time step: '
     +                          ,irec,idt,idtact,'  ',aline
              end if
              idt = idtact
            end if
            if( idt <= 0 ) then
	      isk = isk + 1
	      call dts_format_abs_time(atime,aline)
              write(6,*) '*** zero or negative time step: ',irec,idt
     +                          ,atime,atimeold,'  ',aline
            end if
          end if

	end

c*****************************************************************

	subroutine write_extract(atime,atime0,datetime
     +					,nvar,lmax,dext,d3dext)

	implicit none

	double precision atime,atime0
	integer datetime(2)
	integer nvar,lmax
	real dext(nvar)
	real d3dext(lmax,nvar)

	integer, save :: iu2d = 0
	integer, save :: iu3d = 0

	integer it
	double precision dtime
	character*80 eformat
	character*20 aline

	integer ifileo

	if( iu2d == 0 ) then
	  iu2d = ifileo(88,'out.txt','form','new')
	  write(iu2d,'(a,2i10)') '#date: ',datetime
	  write(eformat,'(a,i3,a)') '(i12,',nvar,'g14.6,a2,a20)'
	  write(6,*) 'using format: ',trim(eformat)
	end if
	if( iu2d == 0 ) then
	  iu3d = ifileo(89,'out.fem','form','new')
	end if

	dtime = atime-atime0
	it = nint(dtime)
	!it = nint(atime-atime1997)
	call dts_format_abs_time(atime,aline)
	write(iu2d,eformat) it,dext,'  ',aline

!	nvers
!        call fem_file_write_header(iformat,iu3d,dtime
!     +                          ,nvers,np,lmax
!     +                          ,nvar,ntype
!     +                          ,nlvddi,hlv,datetime,regpar)
!        call fem_file_write_data(iformat,iu3d
!     +                          ,nvers,np,lmax
!     +                          ,string
!     +                          ,ilhkv,hd
!     +                          ,nlvddi,temp1)

	end

c*****************************************************************

	subroutine femsplit(iformout,ius,dtime,nvers,np
     +			,lmax,nlvddi,ntype
     +			,hlv,datetime,regpar,string
     +			,ilhkv,hd,data)

	implicit none

	integer iformout,ius
	double precision dtime
	integer nvers,np,lmax,nlvddi,ntype
	real hlv(lmax)
	integer datetime(2)
	real regpar(7)
	character(*) string
	integer ilhkv(np)
	real hd(np)
	real data(nlvddi,np)

	integer ivar,nvar
	character*80 file,extra
	character*1 dir
	integer, save :: iusold

	if( ius == 0 ) then
	  call string2ivar_n(string,ivar)
	  call string_direction(string,dir)
	  if( dir == 'y' ) then		!is second part of vector
	    ius = iusold
	  else
	    call alpha(ivar,extra)
	    file = 'out.' // trim(extra) // '.fem'
	    call fem_file_write_open(file,iformout,ius)
	    if( ius <= 0 ) goto 99
	    iusold = ius
	  end if
	end if

	nvar = 1
	if( dir /= ' ' ) nvar = 2

	if( dir /= 'y' ) then
          call fem_file_write_header(iformout,ius,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)
	end if

        call fem_file_write_data(iformout,ius
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvddi,data)

	return
   99	continue
	write(6,*) 'cannot open file ',trim(file)
	stop 'error stop femsplit: cannot open output file'
	end

c*****************************************************************

	subroutine custom_elab(nlvdi,np,string,iv,data)

	implicit none

	integer nlvdi,np,iv
	character*(*) string
	real data(nlvdi,np)

	real fact

	return

	if( string(1:13) /= 'wind velocity' ) return
	if( iv < 1 .or. iv > 2 ) return

	fact = 2.

	!write(6,*) iv,'  ',trim(string)
	write(6,*) 'attention: wind speed changed by a factor of ',fact

	data = fact * data

	end

c*****************************************************************

	subroutine reg_expand_shell(nlvddi,np,lmax,regexpand
     +					,regpar,il,data)

c shell to call expansion routine

	implicit none

	integer nlvddi,np,lmax,regexpand
	real regpar(7)
	integer il(np)
	real data(nlvddi,np)

	integer nx,ny
	real flag

        nx = nint(regpar(1))
        ny = nint(regpar(2))
        flag = regpar(7)

	if( nx*ny /= np ) then
	  write(6,*) 'np,nx,ny,nx*ny: ',np,nx,ny,nx*ny
	  stop 'error stop reg_expand_shell: incompatible params'
	end if

	write(6,*) 'expanding regular grid: ',nx,ny,regexpand

	call reg_expand_3d(nlvddi,nx,ny,lmax,regexpand,flag,data)

	call adjust_reg_vertical(nlvddi,nx,ny,flag,data,il)

	end

c*****************************************************************

	subroutine reg_set_flag(nlvddi,np,il,regpar,data)

	implicit none

	integer nlvddi,np
	integer il(np)
	real regpar(7)
	real data(nlvddi,np)

	integer k,lmax
	real flag

	flag = regpar(7)

	do k=1,np
	  lmax = il(k)
	  data(lmax+1:nlvddi,k) = flag
	end do

	end

c*****************************************************************

	subroutine fem_condense(np,lmax,data,flag,data_profile)

	implicit none

	integer np,lmax
	real data(lmax,np)
	real flag
	real data_profile(lmax)

	integer l,i,nacu
	double precision val,acu

	do l=1,lmax
	  nacu = 0
	  acu = 0.
	  do i=1,np
	    val = data(l,i)
	    if( val /= flag ) then
	      nacu = nacu + 1
	      acu = acu + val
	    end if
	  end do
	  if( nacu == 0 ) then
	    val = flag
	  else
	    val = acu / nacu
	  end if
	  data_profile(l) = val
	end do

	end

c*****************************************************************

	subroutine handle_extract(breg,np,regpar,iextract)

	use clo

	implicit none

	logical breg
	integer np
	real regpar(7)
	integer iextract

	logical berror
	integer inode,icoord,ie
	integer nx,ny,ix,iy,ixx,iyy
	real x0,y0,dx,dy,xp,yp,x,y
	real dist,xydist
	real f(3)
	character*80 snode,scoord

	integer iscanf

        call clo_get_option('node',snode)
        call clo_get_option('coord',scoord)

	inode = iscanf(snode,f,3)
	icoord = iscanf(scoord,f,3)

	iextract = 0

	if( inode == 0 .and. icoord == 0 ) return
	if( inode /= 0 .and. icoord /= 0 ) then
	  write(6,*) 'only one of -node and -coord can be given'
	  stop 'error stop handle_extract: cannot extract'
	end if
	if( inode > 0 .and. inode > 2 ) then
	  write(6,*) 'parse error for -node: need two numbers'
	  write(6,*) '-node:  ',trim(snode)
	  stop 'error stop handle_extract: need 2 values'
	end if
	if( icoord > 0 .and. icoord /= 2 ) then
	  write(6,*) 'parse error for -coord: need two numbers'
	  write(6,*) '-coord:  ',trim(scoord)
	  stop 'error stop handle_extract: need 2 values'
	end if

	if( inode == 1 ) then		!continuous numbering
	  iextract = nint(f(1))
	  if( iextract > np .or. iextract < 1 ) then
	    write(6,*) 'cannot extract this node: ',iextract
	    write(6,*) 'min,max possible: ',1,np
	    stop 'error stop handle_extract: iextract out of range'
	  end if
	end if

	if( ( inode == 2 .or. icoord == 2 ) .and. .not. breg ) then
	  write(6,*) 'need regular file to extract with these values:'
	  write(6,*) '-node:  ',trim(snode)
	  write(6,*) '-coord: ',trim(scoord)
	  stop 'error stop handle_extract: need regular file'
	end if

	nx = nint(regpar(1))
	ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)

	if( inode == 2 ) then
	  ix = nint(f(1))
	  iy = nint(f(2))
	  berror = .false.
	  if( ix < 1 .or. ix > nx ) berror = .true.
	  if( iy < 1 .or. iy > ny ) berror = .true.
	  if( berror ) then
	    write(6,*) 'ix,iy out of range'
	    write(6,*) 'ix,iy,nx,ny: ',ix,iy,nx,ny
	    stop 'error stop handle_extract: out of range'
	  end if
	  iextract = (iy-1)*nx + ix
	end if

	if( icoord == 2 ) then
	  xp = f(1)
	  yp = f(2)
	  ixx = 1
	  iyy = 1
	  xydist = (x0-xp)**2 + (y0-yp)**2
	  do iy=1,ny
	    y = y0 + (iy-1)*dy
	    do ix=1,nx
	      x = x0 + (ix-1)*dx
	      dist = (x-xp)**2 + (y-yp)**2
	      if( dist < xydist ) then
		xydist = dist
	        ixx = ix
	        iyy = iy
	      end if
	    end do
	  end do
	  iextract = (iyy-1)*nx + ixx
	end if

	ie = iextract
	iy =  1 + (ie-1)/nx
	ix = ie - (iy-1)*nx
	xp = x0 + (ix - 1) * dx
	yp = y0 + (iy - 1) * dy
	write(6,*) 'regular grid:          ',nx,ny
	write(6,*) 'extracting point:      ',ix,iy
	write(6,*) 'extracting coords:     ',xp,yp
	write(6,*) 'extracting point (1d): ',iextract

	end

c*****************************************************************

