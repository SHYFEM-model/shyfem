
	program feminf

c writes info on fem file

	implicit none

	character*50 name,string,infile
	integer np,iunit,iout
	integer nvers,lmax,nvar,ntype,nlvdim
	integer nvar0,lmax0,np0
	integer idt,idtact
	double precision dtime,tmin,tmax,dtime0
	double precision atime,atimeold,atimeanf,atimeend
	real dmin,dmax
	integer ierr
	integer irec,i,ich
	integer itype(2)
	integer iformat
	integer datetime(2),dateanf(2),dateend(2)
	real regpar(7)
	logical bdebug,bfirst,bskip,bwrite,bout,btmin,btmax,boutput
	logical bdate		!is date given?
	character*50, allocatable :: strings(:)
	character*20 line
	real,allocatable :: data(:,:,:)
	real,allocatable :: hd(:)
	real,allocatable :: hlv(:)
	integer,allocatable :: ilhkv(:)

c--------------------------------------------------------------
c parameters and command line options
c--------------------------------------------------------------

	bdebug = .false.

        call parse_command_line(infile,bwrite,bout,tmin,tmax)
	bskip = .not. bwrite
	if( bout ) bskip = .false.
	btmin = tmin .ne. -1.
	btmax = tmax .ne. -1.

c--------------------------------------------------------------
c open file
c--------------------------------------------------------------

	if( infile .eq. ' ' ) stop

	np = 0
	call fem_file_read_open(infile,np,iunit,iformat)
	if( iunit .le. 0 ) stop

	write(6,*) 'file name: ',infile
	write(6,*) 'iunit:     ',iunit
	write(6,*) 'format:    ',iformat

c--------------------------------------------------------------
c read first record
c--------------------------------------------------------------

        call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)

	if( ierr .ne. 0 ) goto 99

	write(6,*) 'nvers: ',nvers
	write(6,*) 'np:    ',np
	write(6,*) 'lmax:  ',lmax
	write(6,*) 'nvar:  ',nvar
	write(6,*) 'ntype: ',ntype

	allocate(hlv(lmax))
	call fem_file_make_type(ntype,2,itype)

	call fem_file_read_2header(iformat,iunit,ntype,lmax
     +			,hlv,regpar,ierr)
	if( ierr .ne. 0 ) goto 98

	if( lmax > 1 ) then
	  write(6,*) 'vertical layers: ',lmax
	  write(6,*) hlv
	end if
	if( itype(1) .gt. 0 ) then
	  write(6,*) 'date and time: ',datetime
	end if
	if( itype(2) .gt. 0 ) then
	  write(6,*) 'regpar: ',regpar
	end if

	call fem_file_convert_time(datetime,dtime,atime)

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
	  write(6,*) 'data:  ',i,'  ',string
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

	  call fem_file_convert_time(datetime,dtime,atime)
	  call dts_format_abs_time(atime,line)

	  if( bdebug ) write(6,*) irec,atime,line

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
              call minmax(lmax,np,ilhkv,data(1,1,i),dmin,dmax)
	      write(6,1100) irec,i,atime,dmin,dmax,line
	    end if
	  end do

	  if( irec > 1 ) then
	    if( irec == 2 ) idt = nint(atime-atimeold)
	    idtact = nint(atime-atimeold)
	    if( idtact .ne. idt ) then
	      ich = ich + 1
	      write(6,*) 'change in time step: ',irec,idt,idtact
	      idt = idtact
	    end if
	    if( idt <= 0 ) then
	      write(6,*) 'zero or negative time step: ',irec,idt
     +				,atime,atimeold
	    end if
	  end if
	  atimeend = atime
	end do

c--------------------------------------------------------------
c finish loop - info on time records
c--------------------------------------------------------------

	irec = irec - 1
	write(6,*) 'irec:  ',irec
	call dts_format_abs_time(atimeanf,line)
	write(6,*) 'start time: ',atimeanf,line
	call dts_format_abs_time(atimeend,line)
	write(6,*) 'end time  : ',atimeend,line
	write(6,*) 'idt:   ',idt

	if( ich .gt. 0 ) then
	  write(6,*) ' * warning: time step changed: ',ich
	end if

	close(iunit)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

 1100	format(i6,i3,f15.2,2g16.5,1x,a20)
	stop
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

        subroutine make_time(it,year0,line)

        implicit none

        integer it
        integer year0
        character*(*) line

        integer year,month,day,hour,min,sec

        line = ' '
        if( year0 .le. 0 ) return

        call dts2dt(it,year,month,day,hour,min,sec)
        call dtsform(year,month,day,hour,min,sec,line)

        end

c*****************************************************************

        subroutine write_node(it,nlvdim,np,data)

        implicit none

        integer it
        integer nlvdim,np
        real data(nlvdim,1)

        integer nnodes
        parameter(nnodes=4)
        integer nodes(nnodes)
        save nodes
        data nodes /9442,10770,13210,14219/

        integer n,i

        n = nnodes
        write(90,'(i10,10i6)') it,(ifix(data(1,nodes(i))),i=1,n)

        end

c*****************************************************************

        subroutine write_value(it,nlvdim,np,data)

        implicit none

        integer it
        integer nlvdim,np
        real data(nlvdim,1)

        integer n,nskip,i

        n = 10
        nskip = np/n

        !write(89,*) np,n,nskip,n*nskip
        write(89,'(i10,10i6)') it,(ifix(data(1,i*nskip)),i=1,n)

        end

c*****************************************************************

        subroutine minmax(nlvdim,np,ilhkv,data,vmin,vmax)

        implicit none

        integer nlvdim,np
        integer ilhkv(1)
        real data(nlvdim,1)
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

        subroutine parse_command_line(infile,bwrite,bout,tmin,tmax)

        implicit none

        character*(*) infile
	logical bwrite
	logical bout
	double precision tmin,tmax

        integer i,nc
        character*50 aux

        infile = ' '
	bwrite = .false.
	bout = .false.
	tmin = -1.
	tmax = -1.

        nc = command_argument_count()

        if( nc > 0 ) then
          i = 0
          do
            i = i + 1
            if( i > nc ) exit
            call get_command_argument(i,aux)
            if( aux(1:1) .ne. '-' ) exit
	    if( aux .eq. '-write' ) then
	      bwrite = .true.
	    else if( aux .eq. '-help' .or. aux .eq. '-h' ) then
	      i = nc + 1
              exit
	    else if( aux .eq. '-out' ) then
	      bout = .true.
	    else if( aux .eq. '-tmin' ) then
	      i = i + 1
              if( i > nc ) exit
              call get_command_argument(i,aux)
	      read(aux,'(f14.0)') tmin
	    else if( aux .eq. '-tmax' ) then
	      i = i + 1
              if( i > nc ) exit
              call get_command_argument(i,aux)
	      read(aux,'(f14.0)') tmax
            else
              write(6,*) '*** unknown option: ',aux
	      i = nc + 1
              exit
            end if
          end do
          call get_command_argument(nc,infile)
          if( i .eq. nc ) return
        end if

        write(6,*) 'Usage: feminf [options] fem-file'
        write(6,*) '   options:'
        write(6,*) '      -help|-h    this help'
        write(6,*) '      -write      write min/max of values'
        write(6,*) '      -out        create output file'
        write(6,*) '      -tmax time  only process up to time'
        write(6,*) '      -tmin time  only process starting from time'
        write(6,*) '   With no option gives general information on file'
        write(6,*) '   -write gives detailed info on every record'
        write(6,*) '   -out re-writes output file'

        stop
        end

c*****************************************************************

