
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013-2015,2019  Georg Umgiesser
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
! 13.06.2013	ggu	changed VERS_6_1_65
! 05.11.2014	ggu	changed VERS_7_0_5
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_52
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 27.03.2021	ggu	adapted for other than real output
! 30.03.2021	ggu	some more enhancements

!******************************************************************

!==================================================================
	module mod_check_debug
!==================================================================

	integer, parameter :: ndim = 2000000

	integer, parameter :: iuout = 777

        integer, parameter :: integer_type = 1
        integer, parameter :: real_type = 2
        integer, parameter :: double_type = 3

	real ival1(ndim),ival2(ndim)
	real rval1(ndim),rval2(ndim)
	real dval1(ndim),dval2(ndim)

!==================================================================
	end module mod_check_debug
!==================================================================

	program check_debug

c checks two files written with check_debug from shyfem

	use mod_check_debug

	implicit none

	character*60 name_one,name_two
	character*80 text
	character*80 caption
	logical bcheck,bstop,bverbose
	logical bfirst
	double precision dtime
	integer nt,nf,ny
	integer nrec,ntrec,ierr
	integer i,idiff,idiff_tot
	integer nc

	caption=' variable        nrec     tot dim   first dim'//
     +				'       idiff'

	bstop = .false.			!stop on error
	bstop = .true.			!stop on error
	bverbose = .true.		!only write error
	bverbose = .false.		!only write error
	bcheck = .true.			!check for differences

	nc = command_argument_count()
	if( nc < 1 .or. nc > 2 ) then
	  write(6,*) 'Usage: check_debug file1 [file2]'
	  stop 'error stop check_debug: no files given'
	end if

	if( nc == 1 ) then
	  call get_command_argument(1,name_one)
	  open(1,file=name_one,status='old',form='unformatted')
	  write(6,*) 'file 1: ',trim(name_one)
	  call show_file(1)
	  stop
	end if

	call get_command_argument(1,name_one)
	call get_command_argument(2,name_two)

	open(1,file=name_one,status='old',form='unformatted')
	open(2,file=name_two,status='old',form='unformatted')

	write(6,*) 'file 1: ',trim(name_one)
	write(6,*) 'file 2: ',trim(name_two)

	ntrec = 0
	idiff_tot = 0

	call find_common_start_time

	do while(.true.)

	  call read_time(ntrec,dtime,ierr)
	  if( ierr /= 0 ) exit

	  nrec = 0
	  bfirst = .true.
	  do while(.true.)
	    call read_header(nrec,nt,nf,ny,ierr)
	    if( ierr /= 0 ) exit

	    call read_text(text)
	    call read_val(ny,nt)

	    idiff = 0
	    if( bcheck ) then
	      call check_val(text,nrec,nt,nf,ny,idiff)
	      if( idiff > 0 .or. bverbose ) then
		if( bfirst ) write(6,*) trim(caption)
		bfirst = .false.
	        write(6,1000) trim(text),nrec,nt,nf,idiff
	      end if
	      idiff_tot = idiff_tot + idiff
	    end if

	  end do

	  if( bstop .and. idiff_tot > 0 ) exit
	  if( bverbose ) write(6,*) 'nrecs checked: ',nrec
	end do

	close(1)
	close(2)

	if( idiff_tot > 0 ) then
	  write(6,*) 'detailed errors can be found in file ',iuout
	end if
	write(6,*) 'total time records compared: ',ntrec
	write(6,*) 'total differences found: ',idiff_tot

	if( idiff_tot > 0 ) then
	  write(6,*) '*** errors running check_debug...'
	  call exit(77)
	else
	  write(6,*) 'successful completion of check_debug...'
	end if

	stop
 1000	format(a10,4i12)
	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine find_common_start_time

	implicit none

	integer iu
	double precision atime1,atime2,atime,atime_high

	read(1,end=8) atime1
	read(2,end=8) atime2

	if( atime1 > atime2 ) then
	  iu = 2
	  atime_high = atime1
	  atime = atime2
	else
	  iu = 1
	  atime_high = atime2
	  atime = atime1
	end if

	do
	  if( atime >= atime_high ) exit
	  call skip_time_record(iu)
	  read(iu,end=9) atime
	end do

	if( atime /= atime_high ) goto 9

	backspace(1)
	backspace(2)

	write(6,*) 'common start time: ',atime

	return
    8	continue
	write(6,*) 'atime,atime_high'
	write(6,*) 'cannot read start times'
	stop 'error stop find_common_start_time. no starttime'
    9	continue
	write(6,*) 'atime,atime_high'
	write(6,*) 'cannot find common start time'
	stop 'error stop find_common_start_time. no time found'
	end

!******************************************************************

	subroutine skip_time_record(iu)

	implicit none

	integer iu

	integer ntot

	do
	  read(iu) ntot
	  if( ntot == 0 ) exit
	  read(iu)
	  read(iu)
	end do

	end

!******************************************************************

	subroutine show_file(iu)

	implicit none

	integer iu

	integer ntrec,ios
	double precision dtime

	ntrec = 0

	do
	  read(1,iostat=ios) dtime
	  if( ios /= 0 ) exit
	  ntrec = ntrec + 1
	  write(6,*) 'ntrec = ',ntrec,'  time = ',dtime
	  call skip_time_record(iu)
	end do

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine read_time(ntrec,dtime,ierr)

	use mod_check_debug

	implicit none

	integer ntrec
	double precision dtime
	integer ierr

	double precision dtime1,dtime2

	read(1,end=9) dtime1
	read(2,end=9) dtime2
	if( dtime1 .ne. dtime2 ) goto 99
	dtime = dtime1

	ntrec = ntrec + 1
	write(6,*) 'ntrec = ',ntrec,'  time = ',dtime
	write(iuout,*) 'ntrec = ',ntrec,'  time = ',dtime

	ierr = 0
	return
    9	continue
	ierr = -1
	return
   99	continue
	write(6,*) dtime1,dtime2
	stop 'error stop check_debug: time mismatch'
	end

!******************************************************************

	subroutine read_header(nrec,nt,nf,ny,ierr)

	use mod_check_debug

	implicit none

	double precision dtime
	integer nrec
	integer nt,nf,ny
	integer ierr

        integer nt1,nf1,ny1
        integer nt2,nf2,ny2

        read(1) nt1,nf1,ny1
        read(2) nt2,nf2,ny2

	nrec = nrec + 1

	if( nt1 .ne. nt2 ) goto 98		!total dimension
	if( nf1 .ne. nf2 ) goto 98		!first dimension
	if( ny1 .ne. ny2 ) goto 96		!value type
	nt = nt1
	nf = nf1
	ny = ny1

	ierr = -1
	if( nt .eq. 0 ) return
	if( nt .gt. ndim ) goto 97
	ierr = 0

	return
   96	continue
	write(6,*) 'nrec = ',nrec
	write(6,*) ny1,ny2
	stop 'error stop check_debug: type mismatch'
   97	continue
	write(6,*) 'nrec = ',nrec
	write(6,*) nt,ndim
	stop 'error stop check_debug: dimension'
   98	continue
	write(6,*) 'nrec = ',nrec
	write(6,*) nt1,nt2,nf1,nf2
	stop 'error stop check_debug: size mismatch'
	end

!******************************************************************

	subroutine read_text(text)

	implicit none

	character*(*) text

	character*80 text1,text2

	read(1) text1
	read(2) text2
	if( text1 .ne. text2 ) goto 96
	text = text1

	return
   96	continue
	write(6,*) trim(text1),' - ',trim(text2)
	stop 'error stop check_debug: text mismatch'
	end

!******************************************************************

	subroutine read_val(ny,nt)

	use mod_check_debug

	implicit none

	integer ny,nt

	integer i

	if( ny == integer_type) then
	  read(1) (ival1(i),i=1,nt)
	  read(2) (ival2(i),i=1,nt)
	else if( ny == real_type) then
	  read(1) (rval1(i),i=1,nt)
	  read(2) (rval2(i),i=1,nt)
	else if( ny == double_type) then
	  read(1) (dval1(i),i=1,nt)
	  read(2) (dval2(i),i=1,nt)
	end if

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine check_val(text,nrec,nt,nf,ny,idiff)

	use mod_check_debug

	implicit none

	character*(*) text
	integer nrec
	integer nt,nf,ny
	integer idiff

	write(iuout,*) nrec,nt,nf,ny,trim(text)

	if( ny == integer_type) then
	  call check_val_i(nt,nf,idiff)
	else if( ny == real_type) then
	  call check_val_r(nt,nf,idiff)
	else if( ny == double_type) then
	  call check_val_d(nt,nf,idiff)
	end if

	end

!******************************************************************

	subroutine check_val_i(nt,nf,idiff)

	use mod_check_debug

	implicit none

	integer nt,nf,idiff

	integer i,k,l

	idiff = 0

	do i=1,nt
	  if( ival1(i) .ne. ival2(i) ) then
	    k = 1 + (i-1)/nf
	    l = 1 + mod(i-1,nf)
	    write(iuout,*) k,l,ival1(i),ival2(i)
	    idiff = idiff + 1
	  end if
	end do

	end

!******************************************************************

	subroutine check_val_r(nt,nf,idiff)

	use mod_check_debug

	implicit none

	integer nt,nf,idiff

	integer i,k,l

	idiff = 0

	do i=1,nt
	  if( rval1(i) .ne. rval2(i) ) then
	    k = 1 + (i-1)/nf
	    l = 1 + mod(i-1,nf)
	    write(iuout,*) k,l,rval1(i),rval2(i)
	    idiff = idiff + 1
	  end if
	end do

	end

!******************************************************************

	subroutine check_val_d(nt,nf,idiff)

	use mod_check_debug

	implicit none

	integer nt,nf,idiff

	integer i,k,l

	idiff = 0

	do i=1,nt
	  if( dval1(i) .ne. dval2(i) ) then
	    k = 1 + (i-1)/nf
	    l = 1 + mod(i-1,nf)
	    write(iuout,*) k,l,dval1(i),dval2(i)
	    idiff = idiff + 1
	  end if
	end do

	end

!******************************************************************

