
	program feminf

c writes info on fem file

	implicit none

	character*50 name,string
	integer np,iunit
	integer nvers,lmax,nvar,ntype,nlvdim
	integer it,itanf,itend,idt,itold
	double precision dtime
	real dmin,dmax
	integer ierr
	integer irec,i,nvar0,ich
	integer iformat
	logical bdebug,bfirst,bskip
	character*50, allocatable :: strings(:)
	real,allocatable :: data(:,:)
	real,allocatable :: hd(:)
	integer,allocatable :: ilhkv(:)

	bdebug = .false.
	bskip = .true.
	bskip = .false.

c--------------------------------------------------------------
c open file
c--------------------------------------------------------------

	write(6,*) 'Enter fem-file name: '
	read(5,*) name

	np = 0
	call fem_file_read_open(name,np,iunit,iformat)

	if( iunit .le. 0 ) stop

	write(6,*) 'file name: ',name
	write(6,*) 'iunit:     ',iunit
	write(6,*) 'format:    ',iformat

c--------------------------------------------------------------
c read first record
c--------------------------------------------------------------

	irec = 1

        call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,ierr)

	if( ierr .ne. 0 ) goto 99

	write(6,*) 'nvers: ',nvers
	write(6,*) 'np:    ',np
	write(6,*) 'lmax:  ',lmax
	write(6,*) 'nvar:  ',nvar
	write(6,*) 'ntype: ',ntype

	nlvdim = lmax

	call fem_file_skip_2header(iformat,iunit,lmax,ntype,ierr)
	if( ierr .ne. 0 ) goto 98

	nvar0 = nvar
	allocate(strings(nvar))
	allocate(data(nlvdim,np))
	allocate(hd(np))
	allocate(ilhkv(np))

	do i=1,nvar
	  if( bskip ) then
	    call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
	  else
            call fem_file_read_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,ilhkv
     +                          ,hd
     +                          ,string,nlvdim
     +                          ,data
     +                          ,ierr)
	  end if
	  if( ierr .ne. 0 ) goto 97
	  write(6,*) 'data:  ',i,'  ',string
	  strings(i) = string
	  if( .not. bskip ) then
            call minmax(nlvdim,np,ilhkv,data,dmin,dmax)
	    write(6,*) i,dtime,dmin,dmax
	  end if
	end do

c--------------------------------------------------------------
c loop on other data records
c--------------------------------------------------------------

	bfirst = .true.
	it = nint(dtime)
	itanf = it
	itend = it
	idt = -1
	ich = 0

	do while(.true.)
	  irec = irec + 1
	  itold = itend
          call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,ierr)
	  if( ierr .lt. 0 ) exit
	  if( ierr .gt. 0 ) goto 99
	  if( nvar .ne. nvar0 ) goto 96
	  if( bdebug ) write(6,*) irec,dtime
	  call fem_file_skip_2header(iformat,iunit,lmax,ntype,ierr)
	  if( ierr .ne. 0 ) goto 98
	  do i=1,nvar
	    if( bskip ) then
	      call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
	    else
              call fem_file_read_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,ilhkv
     +                          ,hd
     +                          ,string,nlvdim
     +                          ,data
     +                          ,ierr)
	    end if
	    if( ierr .ne. 0 ) goto 97
	    if( string .ne. strings(i) ) goto 95
	    if( .not. bskip ) then
              call minmax(nlvdim,np,ilhkv,data,dmin,dmax)
	      write(6,*) i,dtime,dmin,dmax
	    end if
	  end do
	  it = nint(dtime)
	  if( bfirst ) then
	    bfirst = .false.
	    idt = it - itold
	  end if
	  if( idt <= 0 ) then
	    write(6,*) 'zero or negative time step: ',irec,it,itold
	  end if
	  if( it-itold .ne. idt ) then
	    ich = ich + 1
	    write(6,*) 'change in time step: ',irec,idt,it-itold
	    idt = it-itold
	  end if
	  itend = it
	end do

c--------------------------------------------------------------
c finish loop - info on time records
c--------------------------------------------------------------

	irec = irec - 1
	write(6,*) 'irec:  ',irec
	write(6,*) 'itanf: ',itanf
	write(6,*) 'itend: ',itend
	write(6,*) 'idt:   ',idt

	if( ich .gt. 0 ) then
	  write(6,*) ' * warning: time step changed: ',ich
	end if

	close(iunit)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	stop
   95	continue
	write(6,*) 'variable ',i
	write(6,*) string
	write(6,*) strings(i)
	write(6,*) 'cannot change description of variables'
	stop 'error stop feminf'
   96	continue
	write(6,*) 'nvar,nvar0: ',nvar,nvar0
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

