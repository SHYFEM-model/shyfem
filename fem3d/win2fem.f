
	program win2fem

! converts files from win-format to fem-format

	implicit none

	character*50 infile
	logical bnew,bpres,bformat
	integer ios,it,id,n,n0,nvar,nvar0,itanf,nkn,i,j,itend
	integer iunit,nvers,ntype,lmax,np,nlvdim,iformat
	integer irec,ifreq,nlen
	double precision dtime
	real hlv(1)
	real hd(1)
	integer ilhkv(1)
	real, allocatable :: data(:,:)
	character*50, allocatable :: strings(:)

	bformat = .true.
	bformat = .false.

	iformat = 0
	if( bformat ) iformat = 1

!-------------------------------------------------------------
! open file
!-------------------------------------------------------------

	call parse_command_line(infile)

	if( infile .eq. ' ' ) then
          write(6,*) 'Enter name of wind file :'
          read(5,'(a)') infile
	end if

        if(infile.eq.' ') stop

        open(1,file=infile,form='unformatted',status='old')

!-------------------------------------------------------------
! read header and see what file it is
!-------------------------------------------------------------

        bnew = .false.
        read(1,iostat=ios) it,id,n,nvar
        if( ios .eq. 0 ) then
	  if( id .ge. 1001 .and. id .le. 1003 ) then
	    bnew = .true.
	  else
	    write(6,*) 'cannot read file format...'
	    stop 'error stop win2fem'
	  end if
	else
	  id = 0
	end if
        backspace(1)

	if( .not. bnew ) then
          read(1,iostat=ios) it,n
          if( ios .ne. 0 ) then
	    write(6,*) 'cannot read file format...'
	    stop 'error stop win2fem'
	  end if
	  nvar = 3
	  id = 1001
          backspace(1)
	end if

!-------------------------------------------------------------
! set parameters and allocate space
!-------------------------------------------------------------

	irec = 0
	itanf = it
	n0 = abs(n)
	nvar0 = nvar
	bpres = bnew .or. n < 0

	!write(6,*) id,itanf,n0,nvar,bnew
	write(6,*) 'file id:    ',id
	write(6,*) 'points:     ',n0
	write(6,*) 'nvar:       ',nvar
	write(6,*) 'new format: ',bnew
	if( id .eq. 1001 ) then
	  write(6,*) 'pressure  : ',bpres
	end if

	allocate(data(n0,nvar))
	allocate(strings(nvar))
	if( id .eq. 1001 ) data(:,3) = 1013.25
	if( id .eq. 1001 ) then
          strings(1) = 'wind velocity in x [m/s]'
          strings(2) = 'wind velocity in y [m/s]'
          strings(3) = 'pressure (atmospheric) [Pa]'
	else if( id .eq. 1002 ) then
	  strings(1) = 'solar radiation [W/m**2]'
          strings(2) = 'air temperature [C]'
          strings(3) = 'humidity [%]'
          strings(4) = 'cloud cover [0-1]'
	else if( id .eq. 1003 ) then
	  strings(1) = 'rain [mm/day]'
	else
	  write(6,*) 'id: ',id
	  stop 'error stop: id'
	end if

	write(6,*) 'content: '
	do i=1,nvar
	  write(6,*) '   ',strings(i)
	end do
	write(6,*) 'working...'

!-------------------------------------------------------------
! open output file and prepare for writing
!-------------------------------------------------------------

	if( bformat ) then
	  open(2,file='out.fem',status='unknown',form='formatted')
	else
	  open(2,file='out.fem',status='unknown',form='unformatted')
	end if

	iunit = 2
	nvers = 0
	ntype = 0
	lmax = 1
	np = n0
	nlvdim = 1

!-------------------------------------------------------------
! loop on input and write
!-------------------------------------------------------------

	do

          if( bnew ) then
            read(1,iostat=ios) it,id,n,nvar
	    nkn = n
          else
            read(1,iostat=ios) it,nkn
	    nvar = 3
          end if
	  if( ios .lt. 0 ) exit
	  if( ios .gt. 0 ) goto 98
          bpres = .false.
          if( nkn .lt. 0 ) then
            nkn = -nkn
            bpres = .true.
          end if

          if( nkn .ne. n0 ) goto 99
	  if( nvar .ne. nvar0 ) goto 99

	  itend = it

          if( bnew ) then
            !read(1) (data(i,1),data(i,2),data(i,3),i=1,nkn)
            read(1) ((data(i,j),j=1,nvar),i=1,nkn)
          else if( bpres ) then
            read(1) (data(i,1),data(i,2),i=1,nkn),(data(i,3),i=1,nkn)
          else
            read(1) (data(i,1),data(i,2),i=1,nkn)
	    do i=1,nkn
	      data(i,3) = 1013.25
	    end do
          end if

	  dtime = it
	  call fem_file_write_header(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,nlvdim,hlv)
	  do i=1,nvar
            call fem_file_write_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,nlvdim,ilhkv
     +                          ,strings(i),hd,data(1,i))
	  end do

	  call progress(irec,2,60)
	end do

!-------------------------------------------------------------
! write final message
!-------------------------------------------------------------

	write(6,*) 
	write(6,*) 'start/end time: ',itanf,itend
	write(6,*) 'output has been written to out.fem'

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	stop
   98	continue
	write(6,*) 'read error: ',ios
	stop 'error stop win2fem'
   99	continue
	write(6,*) 'parameter mismatch: ',nkn,n0,nvar,nvar0
	stop 'error stop win2fem'
	end

c*****************************************************************

	subroutine parse_command_line(infile)

	implicit none

	character*(*) infile

	infile = ' '

	if( command_argument_count() > 0 ) then
	  call get_command_argument(1,infile)
	end if

	if( infile .eq. ' ' ) then
	  write(6,*) 'Usage: win2fem meteo-file'
	  stop
	  !stop 'no file given'
	end if

	end

c*****************************************************************

	subroutine progress(irec,ifreq,nlen)

	implicit none

	integer irec,ifreq,nlen

	irec = irec + 1
	if( mod(irec,ifreq) .eq. 0 ) write(6,'(a1)',advance='no') '.'
	if( mod(irec,ifreq*nlen) .eq. 0 ) write(*,*)	!new line

	end

c*****************************************************************

