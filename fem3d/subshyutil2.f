
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016-2020  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! revision log :
!
! 28.04.2016	ggu	changed VERS_7_5_9
! 25.05.2016	ggu	changed VERS_7_5_10
! 17.06.2016	ggu	changed VERS_7_5_15
! 27.06.2016	ggu	changed VERS_7_5_16
! 05.10.2016	ggu	changed VERS_7_5_19
! 09.10.2017	ggu	changed VERS_7_5_33
! 14.11.2017	ggu	changed VERS_7_5_36
! 25.10.2018	ggu	changed VERS_7_5_51
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 28.01.2020	ggu	minor changes in error message
! 28.04.2023    ggu     update function calls for belem

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

	subroutine shy_write_min_max(nlvdi,nn,il,lmax,ivar,cv3)

	implicit none

	integer nlvdi,nn,lmax,ivar
	integer il(nn)
	real cv3(nlvdi,nn)

	integer l
	real rnull
	real cmin,cmax,cmed
	real cv2(nn)

	rnull = -999.

	do l=1,lmax
	  cv2=cv3(l,:)
	  where( il < l ) cv2 = rnull
	  call mimar(cv2,nn,cmin,cmax,rnull)
          call aver(cv2,nn,cmed,rnull)
          call check1Dr(nn,cv2,0.,-1.,"NaN check","cv2")
	  write(6,1000) 'l,ivar,min,aver,max : ',l,ivar,cmin,cmed,cmax
	end do

 1000	format(a,2i5,3g16.6)
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
	if( .not. shy_is_shy_file(file) .and.
     +      .not. shy_is_lgr_file(file) ) then
	  write(6,*) 'file: ',trim(file)
	  stop 'error stop shyelab: not a valid shy file'
	end if

	id = shy_init(file)
	if( id == 0 ) then
	  stop 'error stop open_next_file: internal error (1)'
	end if

	!write(6,*) '================================'
	!write(6,*) 'reading file: ',trim(file)
	!write(6,*) '================================'

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

        subroutine open_shy_file_0(file,status,nunit)

c open SHY file
c
c nunit is 0 if no other file exists

        use clo

        implicit none

        character*(*) status
        integer nunit

        character*80 file
        integer ifileo

        nunit = 0
        if( file == ' ' ) return

        nunit = ifileo(0,file,'unform',status)

        if( nunit .le. 0 ) then
          write(6,*) 'file: ',trim(file)
          stop 'error stop open_shy_file: opening file'
        end if

        end

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

	subroutine read_records(id,dtime,bhydro,nvar,nndim,nlvddi
     +				,ivars,cv3,cv3all,ierr)

	use elabutil
	use shyfile
	use shyutil

	implicit none

	integer id
	double precision dtime
	logical bhydro
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
	logical bfirst,belem
	double precision dtvar

	shy_znv = 0.
	shy_zenv = 0.

	iv = 0
	bfirst = .true.
	bzeta = .false.
	belem = .false.

	do
	  iv = iv + 1
	  if( iv > nvar ) exit
	  belem = ( bhydro .and. iv >= 2 )
	  cv3 = 0.
	  call shy_read_record(id,dtime,belem
     +				,ivar,n,m,lmax,nlvddi,cv3,ierr)
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
	stop 'error stop read_records: dimension error...'
   75	continue
        write(6,*) 'error in reading file : ',ierr
	stop 'error stop read_records: reading file'
   76	continue
        write(6,*) 'end of file between variables'
	stop 'error stop read_records: EOF unexpected'
   85	continue
	write(6,*) 'dtime,dtvar,iv,ivar,nvar: ',dtime,dtvar,iv,ivar,nvar
	stop 'error stop read_records: time mismatch'
	end

!***************************************************************
!***************************************************************
!***************************************************************

