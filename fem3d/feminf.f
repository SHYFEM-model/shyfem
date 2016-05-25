!
! info on fem files
!
! revision log :
!
! 14.01.2015    ggu     finished and cleaned
! 10.02.2015    ggu     option for not complaining about changed time step
! 07.04.2016    ggu     option node added to extract value at node
!
!******************************************************************

	program feminf

c writes info on fem file

	use clo

	implicit none

	character*80 infile
	integer i,nfile
	logical bdebug,bskip,bwrite,bquiet

	bdebug = .true.
	bdebug = .false.

c--------------------------------------------------------------
c parameters and command line options
c--------------------------------------------------------------

	call clo_init('feminf','fem-file(s)','1.3')

	call clo_add_info('returns info on a fem file')
	call clo_add_option('node node',-1,'extract value for node')
	call clo_add_option('write',.false.,'write min/max of values')
	call clo_add_option('quiet',.false.,'do not be verbose')
	call clo_add_option('no_dt',.false.
     +				,'do not check for change of time step')

	call clo_parse_options(1)  !expecting (at least) 1 file after options

	call clo_get_option('write',bwrite)
	call clo_get_option('quiet',bquiet)

	bskip = .not. bwrite

	nfile = clo_number_of_files()

	if( bdebug ) then
	  write(6,*) nfile
	  write(6,*) bwrite,bskip,bquiet
	end if

c--------------------------------------------------------------
c loop on files
c--------------------------------------------------------------

        do i=1,nfile
          call clo_get_file(i,infile)
          if( infile .ne. ' ' ) call feminf_file(infile)
        end do

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine feminf_file(infile)

c writes info on fem file

	use clo

	implicit none

	character*(*) infile

	character*50 name,string
	integer np,iunit,iout
	integer nvers,lmax,nvar,ntype
	integer nvar0,lmax0,np0
	integer idt,idtact
	double precision dtime,tmin,tmax
	double precision atime,atimeold,atimeanf,atimeend
	real dmin,dmax
	integer node
	integer ierr
	integer nfile
	integer irec,i,ich,nrecs
	integer itype(2)
	integer iformat
	integer datetime(2),dateanf(2),dateend(2)
	real regpar(7)
	real value
	logical bdebug,bfirst,bskip,bwrite
	logical bquiet,bnodt,bcheckdt
	character*50, allocatable :: strings(:)
	character*20 line,aline
	real,allocatable :: data(:,:,:)
	real,allocatable :: hd(:)
	real,allocatable :: hlv(:)
	integer,allocatable :: ilhkv(:)

	bdebug = .true.
	bdebug = .false.

        datetime = 0
        irec = 0

	call clo_get_option('node',node)
	call clo_get_option('write',bwrite)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('no_dt',bnodt)

	bskip = .not. bwrite .and. node < 0
	bcheckdt = .not. bnodt

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

	call dts_convert_to_atime(datetime,dtime,atime)

	nvar0 = nvar
	lmax0 = lmax
	np0 = np
	allocate(strings(nvar))
	allocate(data(lmax,np,nvar))
	allocate(hd(np))
	allocate(ilhkv(np))

	do i=1,nvar
	  call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
	  if( ierr .ne. 0 ) goto 97
	  write(6,*) 'data:   ',i,'  ',string
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

	do 
	  irec = irec + 1
	  atimeold = atime
          call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	  if( ierr .lt. 0 ) exit
	  if( ierr .gt. 0 ) goto 99
	  if( nvar .ne. nvar0 ) goto 96
	  if( lmax .ne. lmax0 ) goto 96
	  if( np .ne. np0 ) goto 96

	  call dts_convert_to_atime(datetime,dtime,atime)
	  call dts_format_abs_time(atime,aline)

	  if( bdebug ) write(6,*) irec,atime,aline

	  call fem_file_read_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar,ierr)
	  if( ierr .ne. 0 ) goto 98

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
	    if( bwrite ) then
              call minmax_data(lmax,np,ilhkv,data(1,1,i),dmin,dmax)
	      write(6,1100) irec,i,atime,dmin,dmax,aline
 1100	      format(i6,i3,f15.2,2g16.5,1x,a20)
	    end if
	    if( node >= 0 .and. node <= np ) then
	      if( node > 0 ) then
		value = data(1,node,i)
	      else
		value = sum(data(1,:,i))/np
	      end if
	      write(33,*) atime,i,value
	    end if
	  end do

	  if( irec > 1 ) then
	    if( irec == 2 ) idt = nint(atime-atimeold)
	    idtact = nint(atime-atimeold)
	    if( idtact .ne. idt ) then
	      ich = ich + 1
	      if( bcheckdt ) then
	        write(6,*) '* change in time step: '
     +				,irec,idt,idtact,aline
	      end if
	      idt = idtact
	    end if
	    if( idt <= 0 ) then
	      write(6,*) '*** zero or negative time step: ',irec,idt
     +				,atime,atimeold,aline
	    end if
	  end if
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
	  write(6,*) 'idt:     irregular'
	end if

	if( ich .gt. 0 ) then
	  write(6,*) '* warning: time step changed: ',ich
	end if

	close(iunit)

	deallocate(strings)
	deallocate(data)
	deallocate(hd)
	deallocate(ilhkv)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	return
   95	continue
	write(6,*) 'variable ',i
	write(6,*) string
	write(6,*) strings(i)
	write(6,*) 'cannot change description of variables'
	stop 'error stop feminf'
   96	continue
	write(6,*) 'nvar,nvar0: ',nvar,nvar0
	write(6,*) 'lmax,lmax0: ',lmax,lmax0	!this might be relaxed
	write(6,*) 'np,np0:     ',np,np0	!this might be relaxed
	write(6,*) 'cannot change number of variables'
	stop 'error stop feminf'
   97	continue
	write(6,*) 'record: ',irec
	write(6,*) 'cannot read data record of file'
	stop 'error stop feminf'
   98	continue
	write(6,*) 'record: ',irec
	write(6,*) 'cannot read second header of file'
	stop 'error stop feminf'
   99	continue
	write(6,*) 'record: ',irec
	write(6,*) 'cannot read header of file'
	stop 'error stop feminf'
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

