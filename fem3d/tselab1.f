
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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
! 14.11.2017	ggu	changed VERS_7_5_36
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 15.05.2019	ggu	new option -date0
! 16.05.2019	ggu	use sdate0 for date string

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine tselab

c writes info on ts file

	use clo
        use elabutil
        use elabtime
        use elabtime

	implicit none

	character*50 name,string
	integer np,iunit
	integer nvers,lmax,nvar,ntype
	integer lmax0,np0
	integer idt,idtact
	double precision dtime,dtime0
	double precision atime,atnew,atold,atfirst,atlast
	double precision atime0,atime0e
	real dmin,dmax
	integer ierr
	integer nfile
	integer nrec,i,ich,nrecs,iv
	integer date,time
	integer datetime(2)
	logical bfirst,bskip,debug
	character*20 dline
	character*20 format
	character*80 varline
	integer,allocatable :: ivars(:)
	character*80,allocatable :: strings(:)
	real,allocatable :: data(:)
	real,allocatable :: data_minmax(:,:)
	integer, save :: iout = 0
	integer, save :: nvar0 = 0

	logical check_ts_file

	debug = .false.
	datetime = 0
	nrec = 0

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

        call elabutil_init('TS','tselab')

        call clo_reset_files
        call clo_get_next_file(infile)
	if( infile .eq. ' ' ) stop

        if( .not. check_ts_file(infile) ) then
          stop 'error stop tselab: not a valid time series file'
        end if

	if( bconvert ) bout = .true.

c--------------------------------------------------------------
c open file
c--------------------------------------------------------------

	atime0e = 0.
	if( bout ) then		!see if we get extra information on time
	  call ts_get_extra_time(infile,dtime,datetime)
	  if( datetime(1) > 0 ) then
	    dtime0 = 0.
	    call dts_convert_to_atime(datetime,dtime0,atime)
	    atime0e = atime - dtime
	  else if( sdate0 /= ' ' ) then
	    call dts_string2time(sdate0,atime0e,ierr)
	    if( ierr /= 0 ) then
	      write(6,*) 'cannot parse date...'
	      write(6,*) 'date0 = ',trim(sdate0)
              stop 'error stop tselab: date0 invalid'
	    end if
	    write(6,*) 'using date for conversion: ',trim(sdate0)
	  end if
	end if

	nvar = 0
	call ts_open_file(infile,nvar,datetime,varline,iunit)
	if( iunit .le. 0 ) stop

	allocate(data(nvar))
	allocate(data_minmax(2,nvar))
	allocate(strings(nvar))
	allocate(ivars(nvar))

	strings = ' '
	call parse_varline(varline,nvar,strings,ivars)

	if( .not. bquiet ) then
	  write(6,*) 'file name : ',trim(infile)
	  write(6,*) 'columns   : ',nvar
	  if( datetime(1) > 0 ) then
	    write(6,*) 'datetime  : ',datetime
	  else if( atime0e > 0 ) then
	    call dts_format_abs_time(atime0e,dline)
	    write(6,*) 'extra time information: ',dline
	  end if
	  if( varline /= ' ' ) then
	    write(6,*) 'variables contained in file: '
            write(6,*) '   varnum     varid    varname'
            do iv=1,nvar
              write(6,'(2i10,4x,a)') iv,ivars(iv),trim(strings(iv))
	    end do
	  end if
	end if

	if( nvar0 == 0 ) nvar0 = nvar
	if( nvar /= nvar0 ) then
	  write(6,*) 'nvar,nvar0: ',nvar,nvar0
	  stop 'error stop: nvar /= nvar0'
	end if

	if( bout .and. iout == 0 ) then
	  iout = 1
	  open(iout,file='out.txt',form='formatted',status='unknown')
	  !write(format,'(a,i3,a)') '(f18.2,',nvar,'g14.6)'
	  write(format,'(a,i3,a)') '(a20,',nvar,'g14.6)'
	  if( debug ) write(6,*) 'used format: ',trim(format)
	end if

	atime0 = 0
	if( bout ) then
	  if( datetime(1) /= 0 ) then
	    dtime = 0.
	    call dts_convert_to_atime(datetime,dtime,atime)
	    atime0 = atime
	    if( debug ) write(6,*) 'using absolute date from file'
	  else if( atime0e /= 0 ) then
	    atime0 = atime0e
	    if( debug ) write(6,*) 
     +			'using absolute date from extra information'
	  else
	    write(6,*) 'no absolute time... cannot convert'
	    stop 'error stop: missing absolute time'
	  end if
	  call dts_format_abs_time(atime0,dline)
	  if( debug ) write(6,*) 'absolute date reference: ',dline
	end if

        date = 0
        time = 0
        call elabtime_date_and_time(date,time)  !we work with absolute time
	call elabtime_set_minmax(stmin,stmax)

c--------------------------------------------------------------
c read first record
c--------------------------------------------------------------

	call ts_read_next_record(iunit,nvar,dtime,data,datetime,ierr)
	if( ierr .ne. 0 ) goto 97

	call dts_convert_to_atime(datetime,dtime,atime)
	if( datetime(1) == 0 ) atime = atime0 + dtime

c--------------------------------------------------------------
c close and re-open file
c--------------------------------------------------------------

	close(iunit)

	if( binfo ) return

	nvar = 0
	call ts_open_file(infile,nvar,datetime,varline,iunit)
	if( iunit .le. 0 ) stop

c--------------------------------------------------------------
c loop on all records
c--------------------------------------------------------------

	nvar0 = nvar
	nrec = 0
	idt = 0
	ich = 0
	atfirst = atime
	atlast = atime

	if( .not. bquiet ) write(6,*)

	do 
	  atold = atime

	  call ts_read_next_record(iunit,nvar,dtime,data,datetime,ierr)
	  if( ierr < 0 ) exit
	  if( ierr > 0 ) goto 97
	  nrec = nrec + 1

	  if( nvar .ne. nvar0 ) goto 96

	  call dts_convert_to_atime(datetime,dtime,atime)
	  if( datetime(1) == 0 ) atime = atime0 + dtime
	  call dts_format_abs_time(atime,dline)

	  atnew = atime
          if( elabtime_over_time(atime,atnew,atold) ) exit
          if( .not. elabtime_in_time(atime,atnew,atold) ) cycle
          !if( .not. elabtime_check_time(atime,atnew,atold) ) cycle

	  if( bsdebug ) write(6,*) nrec,atime,dline

	  if( bwrite ) then
            call minmax_ts(nvar,data,data_minmax)
	  end if
	  if( bverb ) write(6,*) nrec,atime,dline

          if( bcheck ) then
            call fem_check(atime,1,1,nvar,data,flag
     +                          ,strings,scheck,bquiet)
          end if

	  bskip = .false.
	  if( nrec > 1 ) then
	    if( nrec == 2 ) idt = nint(atime-atold)
	    idtact = nint(atime-atold)
	    if( idtact .ne. idt ) then
	      ich = ich + 1
	      if( bcheckdt ) then
	        call dts_format_abs_time(atime,dline)
	        write(6,'(a,a,3i8)') '* change in time step: '
     +     					,dline,nrec,idt,idtact
	      end if
	      idt = idtact
	    end if
	    bskip = idt <= 0
	    if( bcheckdt .and. bskip ) then
	      write(6,*) '*** zero or negative time step: ',nrec,idt
	      call dts_format_abs_time(atold,dline)
	      write(6,*) '    old time: ',dline
	      call dts_format_abs_time(atime,dline)
	      write(6,*) '    new time: ',dline
	    end if
	  end if

	  if( bout .and. .not. bskip ) then
	    call dts_format_abs_time(atime,dline)
	    write(1,format) dline,data(1:nvar)
	  end if

	  atlast = atime
	end do

        if( bcheck ) then       !write final data
          atime = -1.
          call fem_check(atime,1,1,nvar,data,flag
     +                          ,strings,scheck,bquiet)
        end if

c--------------------------------------------------------------
c finish loop - info on time records
c--------------------------------------------------------------

	if( bwrite ) then
	  write(6,*) 'min/max of data: '
	  do i=1,nvar
	    write(6,*) i,data_minmax(1,i),data_minmax(2,i)
	  end do
	end if

        if( .not. bsilent ) then
          write(6,*)
	  call dts_format_abs_time(atfirst,dline)
	  write(6,*) 'first time record: ',dline
	  call dts_format_abs_time(atlast,dline)
	  write(6,*) 'last time record:  ',dline

          write(6,*)
          write(6,*) nrec ,' time records read'
          !write(6,*) nelab,' time records elaborated'
          !write(6,*) nout ,' time records written to file'
          write(6,*)
	end if

	if( bcheckdt ) then
	 if( ich .gt. 0 ) then
	  write(6,*) '* warning: changes in time step: ',ich
	 else if( .not. bquiet ) then
	  write(6,*) 'idt:    ',idt
	 end if
	end if
	if( bout .and. .not. bquiet ) then
	  write(6,*) 'output written to file out.txt'
	end if

        if( bcheck ) then       !write final message
          atime = -2.
          call fem_check(atime,1,1,nvar,data,flag
     +                          ,strings,scheck,bsilent)
        end if

	close(iunit)

	deallocate(strings)
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
	write(6,*) 'error code ',ierr
	write(6,*) 'error reading record ',nrec+1
	write(6,*) 'of file ',trim(infile)
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

	subroutine parse_varline(varline,nvar,strings,ivars)

	use shyfem_strings

	implicit none

	character*(*) varline
	integer nvar
	character*(*) strings(nvar)
	integer ivars(nvar)

	integer nvars,iv,ivar
	character*80 name,full
	character*80 s(nvar+1)
	integer iscans

	strings = ' '
	if( varline == ' ' ) return

	nvars = iscans(varline,s,0)
	nvars = nvars - 1	!subtract time column

	if( nvars /= nvar ) then
	  write(6,*) 'number of variables on varline /= '//
     +			'number of data columns'
	  write(6,*) 'varline: ',trim(varline)
	  write(6,*) 'nvars,nvar: ',nvars,nvar
	  write(6,*) 'not using varline...'
	  varline = ' '
	  return
	end if

	nvars = iscans(varline,s,nvar+1)

	do iv=1,nvar
	  strings(iv) = s(iv+1)		!get rid of time column
	  name = strings(iv)
	  call strings_get_ivar(name,ivar)
	  if( ivar >= 0 ) then
	    call strings_get_full_name(ivar,full)
	    strings(iv) = full
	  end if
	  ivars(iv) = ivar
	end do

	end

c*****************************************************************

