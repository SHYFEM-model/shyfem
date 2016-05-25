
!***************************************************************

	subroutine shy_write_time(bdate,dtime,atime,ivar)

	implicit none

	logical bdate
	double precision dtime,atime
	integer ivar

	character*20 dline

	dline = ' '
	!if( bdate ) call dtsgf(it,dline)
	if( bdate ) call dts_format_abs_time(atime,dline)
	write(6,*) 'time : ',dtime,' ',dline,'  ivar : ',ivar

	end

!***************************************************************

	subroutine shy_write_min_max(nlvdi,nn,lmax,cv3)

	implicit none

	integer nlvdi,nn,lmax
	real cv3(nlvdi,nn)

	integer l
	real rnull
	real cmin,cmax,cmed
	real cv2(nn)

	do l=1,lmax
	  cv2=cv3(l,:)
	  call mimar(cv2,nn,cmin,cmax,rnull)
          call aver(cv2,nn,cmed,rnull)
          call check1Dr(nn,cv2,0.,-1.,"NaN check","cv2")
	  write(6,1000) 'l,min,max,aver : ',l,cmin,cmax,cmed
	end do

 1000	format(a,i5,3g16.6)
	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine open_new_file(ifile,id,atstart)

	use shyfile
	use elabtime

! opens a new file, closing the old one and checking the next file
!
! ifile gives number of old file on command line
!	should be 0 for first call, is incremented at return
! id is the id of the old opened file
!	should be 0 for first call, is id of new file
! atstart is absolute time of next file that should be opened
!	is -1 if no next file exists

	implicit none

	integer ifile			!number of old file on command line
	integer id			!id of old (in) and new (out) file
	double precision atstart	!absolute time of start of next file

	integer idold
	integer date,time

	ifile = ifile + 1
	idold = id

        call open_next_file(ifile,idold,id)
        call shy_close(idold)

        call get_start_of_next_file(ifile+1,atstart)

        call shy_get_date(id,date,time)
        call elabtime_date_and_time(date,time) 

	end

!***************************************************************

	subroutine open_next_file(ifile,idold,id)

	use clo

	implicit none

	integer ifile
	integer idold,id

	character*80 file

	call clo_get_file(ifile,file)
	call open_next_file_by_name(file,idold,id)

	end

!***************************************************************

	subroutine open_next_file_by_name(file,idold,id)

	use shyfile

	implicit none

	character*80 file
	integer idold,id

	integer ierr

	if( .not. shy_exist_file(file) ) then
	  write(6,*) 'file: ',trim(file)
	  stop 'error stop shyelab: file not existing'
	end if
	if( .not. shy_is_shy_file(file) ) then
	  write(6,*) 'file: ',trim(file)
	  stop 'error stop shyelab: not a valid shy file'
	end if

	id = shy_init(file)
	if( id == 0 ) then
	  stop 'error stop open_next_file: internal error (1)'
	end if

	write(6,*) '================================'
	write(6,*) 'reading file: ',trim(file)
	write(6,*) '================================'

	call shy_read_header(id,ierr)
	if( ierr /= 0 ) then
	  write(6,*) 'error reading header, ierr = ',ierr
	  write(6,*) 'file = ',trim(file)
	  stop 'error stop open_next_file: reading header'
	end if

	if( idold > 0 ) then
	  if( .not. shy_are_compatible(idold,id) ) then
	    write(6,*) 'files are not compatible'
	    write(6,*) 'info from first file: '
	    call shy_info(idold)
	    write(6,*) 'info from second file: '
	    call shy_info(id)
	    stop 'error stop open_next_file: not compatible'
	  end if
	  !call shy_close(idold)
	end if

	!call shy_info(id)

	end

!***************************************************************

	!subroutine open_shy_file(file)

	!end

!***************************************************************

	subroutine get_start_of_next_file(ifile,atstart)

	use clo

	implicit none

	integer ifile
	double precision atstart

	logical bok
	integer datetime(2)
	double precision dtime
	character*80 file

	atstart = -1.

	call clo_get_file(ifile,file)
	if( file == ' ' ) return

	call shy_get_tstart(file,datetime,dtime,bok)
	if( .not. bok ) return

	call dts_convert_to_atime(datetime,dtime,atstart)

	end 

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine read_records(id,dtime,nvar,nndim,nlvddi
     +				,ivars,cv3,cv3all,ierr)

	use elabutil
	use shyfile
	use shyutil

	implicit none

	integer id
	double precision dtime
	integer nvar,nndim
	integer nlvddi
	integer nkn,nel
	integer ivars(4,nvar)
	real cv3(nlvddi,nndim)
	real cv3all(nlvddi,nndim,0:nvar)
	integer ierr

	integer nexp

	integer iv
	integer ivar,n,m,lmax
	logical bfirst
	double precision dtvar

	shy_znv = 0.
	shy_zenv = 0.

	iv = 0
	bfirst = .true.
	bzeta = .false.

	do
	  iv = iv + 1
	  if( iv > nvar ) exit
	  call shy_read_record(id,dtime,ivar,n,m,lmax,nlvddi,cv3,ierr)
          if( ierr .gt. 0 ) goto 75
          if( ierr .ne. 0 ) exit
	  nexp = n * m
	  bzeta = ivar == -1
	  if( bzeta ) iv = iv - 1
	  ivars(1,iv) = n
	  ivars(2,iv) = m
	  ivars(3,iv) = lmax
	  if( nexp > nndim ) goto 74
	  if( lmax > nlvddi ) goto 74
	  if( bfirst ) dtvar = dtime
	  if( dtvar /= dtime ) goto 85
	  ivars(4,iv) = ivar
	  cv3all(:,:,iv) = cv3(:,:) * fact
	  if( abs(ivar) == 1 ) then		! water level
	  !write(6,*) 'water level: ',iv,n,m
	    if( ivar == -1 .or. iv == 1 ) then
	      shy_znv = cv3(1,1:n)
	    else
	      !zenv = cv3(1,1:3*n)
	      shy_zenv = reshape(cv3(1,1:3*n),(/3,n/))	!FIXME
	    end if
	  end if
	end do

	if( ierr /= 0 .and. iv .ne. 1 ) goto 76

	return
   74	continue
	write(6,*) n,m,nexp,nndim,lmax,nlvddi
	stop 'error stop shyelab: dimension error...'
   75	continue
        write(6,*) 'error in reading file : ',ierr
	stop 'error stop shyelab: reading file'
   76	continue
        write(6,*) 'end of file between variables'
	stop 'error stop shyelab: EOF unexpected'
   85	continue
	write(6,*) 'dtime,dtvar,iv,ivar,nvar: ',dtime,dtvar,iv,ivar,nvar
	stop 'error stop shyelab: time mismatch'
	end

!***************************************************************
!***************************************************************
!***************************************************************

