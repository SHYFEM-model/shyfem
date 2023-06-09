
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

! compares two debug files or writes info on one

! revision log :
!
! 03.06.2020	ggu	newly written based on check_debug.f
! 10.11.2021    ggu     avoid warning for stack size
! 30.03.2022    ggu     ntime was not initialized
! 02.04.2022    ggu     new routine i_info()
! 06.04.2022    ggu     new code for handling double variables
! 07.04.2022    ggu     command line options introduced
! 09.04.2022    ggu     minor changes
! 02.05.2022    ggu     new option -nodiff
! 18.05.2022    ggu     new option -maxdiff
! 05.10.2022    ggu     handle different initial time and header
! 28.03.2023    ggu     code refactorying, new options -summary
! 13.04.2023    ggu     avoid compiler errors for d/i/r_info
! 05.06.2023    ggu     handle exceptions, show records with highest errrors
! 09.06.2023    ggu     handle non existing files

! note :
!
! debug_file:
!	debug_time_record
!	debug_time_record
!	...
! debug_time_record
!	dtime
!	debug_data_record
!	debug_data_record
!	...
!	0,0,0
! debug_data_record
!	nh,nv,nt
!	text
!	val(nh*nv)
!
! first debug_time_record is header and has dtime == -1.
!
!**************************************************************************

!==========================================================================
	module mod_shympi_debug
!==========================================================================

	implicit none

	logical, save :: bstop = .true.		!stop on error
	logical, save :: bcheck = .true.	!check for differences

	logical, save :: bsilent = .false.	!be verbose
	logical, save :: bquiet = .false.	!be verbose
	logical, save :: bnodiff = .false.	!do not show differences
	logical, save :: bsummary = .false.	!do only show summary
	logical, save :: bverbose = .false.	!be verbose
	logical, save :: bbalance = .false.	!balance time records

	double precision, save :: maxdiff = 0.	!maximum difference for error

	integer, parameter :: imax = 20		!number of errors shown

	integer :: ifill
	integer :: index(imax)
	double precision :: diffmin
	double precision :: diffs(imax)

	integer, allocatable :: ival1(:),ival2(:)
	real, allocatable :: rval1(:),rval2(:)
	double precision, allocatable :: dval1(:),dval2(:)

	integer, save :: nipv = 0
	integer, save :: nipev = 0

!==========================================================================
	end module mod_shympi_debug
!==========================================================================

	program check_shympi_debug

        use clo
	use mod_shympi_debug

	implicit none

	integer nc,ierr

	call check_shympi_debug_init

	nc = clo_number_of_files()

	if( .not. bquiet ) write(6,*) 'running check_shympi_debug...'

	if( nc == 0 ) then
	  call clo_usage
	else if( nc == 1 ) then
	  call read_dbg_file
	else if( nc == 2 ) then
	  call compare_files(ierr)
	else
	  write(6,*) 'nc = ',nc
	  stop 'error stop check_shympi_debug: wrong number of files'
	end if

	if( ierr > 0 ) then
	  if( ierr == 99 ) ierr = 100	!terrible hack - FIXME
	  call exit(ierr)
	else
	  call exit(99)
	end if

	end

!**************************************************************************

	subroutine read_dbg_file

! reads one file and outputs info

        use clo
	use mod_shympi_debug

	implicit none

	integer nc
	integer ntime,nrec
	integer nh,nv,nt
	integer ios
	double precision dtime
	character*60 name_one,text

	call clo_get_file(1,name_one)

	open(1,file=name_one,status='old',form='unformatted',iostat=ios)

	if( ios /= 0 ) then
	  write(6,*) 'no such file: ',trim(name_one)
	  stop 'error opening file'
	end if

	if( .not. bquiet ) write(6,*) 'file 1: ',trim(name_one)

	ntime = 0

	do while(.true.)

	  read(1,end=9) dtime
	  if( .not. bquiet ) write(6,*) 'time = ',dtime
	  ntime = ntime + 1

	  if( bverbose ) then
	    write(6,*) '       irec          nh          nv' //
     +			'        type name'
	  end if

	  nrec = 0
	  do while(.true.)
	    read(1,end=9) nh,nv,nt
	    if( nt == 0 ) exit
	    nrec = nrec + 1
	    read(1,end=9) text
	    read(1)
	    if( bverbose ) write(6,*) nrec,nh,nv,nt,trim(text)
	  end do
	end do

    9	continue
	if( .not. bsilent ) write(6,*) 'time records read: ',ntime

	end

!**************************************************************************

	subroutine compare_files(idiff_end)

! checks two files written with check_debug from shyfem

	use clo
	use mod_shympi_debug

	implicit none

	INTERFACE 
	  subroutine allocate_arrays(nsize,ndim
     +			,ival1,ival2,rval1,rval2,dval1,dval2)
	  integer nsize,ndim
	  integer, allocatable :: ival1(:),ival2(:)
	  real, allocatable :: rval1(:),rval2(:)
	  double precision, allocatable :: dval1(:),dval2(:)
	  END subroutine
	END INTERFACE

	INTERFACE 
	  subroutine save_int(n,ival,niv,iv)
	  integer n,niv
	  integer ival(n)
	  integer, allocatable :: iv(:)
	  END subroutine
	END INTERFACE

	INTERFACE 
	  subroutine read_header(iunit,nipv,nipev,ipv,ipev)
	  integer iunit
	  integer nipv,nipev
	  integer, allocatable :: ipv(:),ipev(:)
	  END subroutine
	END INTERFACE

	INTERFACE 
	  subroutine d_info(nh,nv,val1,val2,ipv,ipev,text)
	  use mod_shympi_debug
	  integer nh,nv
	  double precision val1(nh*nv)
	  double precision val2(nh*nv)
	  integer ipv(nipv),ipev(nipev)
	  !integer ipv(:),ipev(:)
	  character*(*) text
	  END subroutine
	END INTERFACE

	INTERFACE 
	  subroutine r_info(nh,nv,val1,val2,ipv,ipev,text)
	  use mod_shympi_debug
	  integer nh,nv
	  real val1(nh*nv)
	  real val2(nh*nv)
	  integer ipv(nipv),ipev(nipev)
	  !integer ipv(:),ipev(:)
	  character*(*) text
	  END subroutine
	END INTERFACE

	INTERFACE 
	  subroutine i_info(nh,nv,val1,val2,ipv,ipev,text)
	  use mod_shympi_debug
	  integer nh,nv
	  integer val1(nh*nv)
	  integer val2(nh*nv)
	  integer ipv(nipv),ipev(nipev)
	  !integer ipv(:),ipev(:)
	  character*(*) text
	  END subroutine
	END INTERFACE

	integer idiff_end

	integer, save :: ndim = 0

	character*60 name_one,name_two
	character*80 text1,text2,text
	logical bshowdiff,bheader,blcheck
	integer nt1,nt2,nt
	integer nh1,nh2,nh
	integer nv1,nv2,nv
	integer nrec,ntot,ntime
	integer i,idiff,idiff_tot,idiff_rec,ntrerr
	integer nc
	!integer nipv,nipev
	integer ios,ierr
	integer, allocatable :: ipv(:),ipev(:)
	real rdiff
	double precision dtime,dtime1,dtime2
	double precision diff,rdiff_max

	call clo_get_file(1,name_one)
	call clo_get_file(2,name_two)

	open(1,file=name_one,status='old',form='unformatted',iostat=ios)
	if( ios /= 0 ) then
	  write(6,*) 'no such file: ',trim(name_one)
	  stop 'error opening file 1'
	end if
	open(2,file=name_two,status='old',form='unformatted',iostat=ios)
	if( ios /= 0 ) then
	  write(6,*) 'no such file: ',trim(name_two)
	  stop 'error opening file 2'
	end if

	if( .not. bquiet ) then
	  write(6,*) 'file 1: ',trim(name_one)
	  write(6,*) 'file 2: ',trim(name_two)
	end if

	if( bbalance ) then
	  call read_header(1,nipv,nipev,ipv,ipev)
	  call read_header(2,nipv,nipev,ipv,ipev)
	  call time_balance
	end if

	idiff_tot = 0
	idiff_rec = 0
	rdiff_max = 0
	ntrerr = 0		!number of time records with error
	ntime = 0
	nipv = 0
	nipev = 0

	do while(.true.)

	  if( bsummary .and. idiff_rec > 0 ) then
	    write(6,*) '  nerr,maxerr: ',idiff_rec,rdiff_max
	  end if

	  call read_time_header(dtime,ntime,ierr)
	  if( ierr == -1 ) exit

	  idiff_rec = 0
	  rdiff_max = 0
	  bheader = .true.		!write header for differences

	  nrec = 0
	  do while(.true.)

	    call read_data_header(nh,nv,nt,ntot,nrec,text,ierr)
	    if( nt == 0 .or. ierr == -1 ) exit

	    call allocate_arrays(ntot,ndim
     +			,ival1,ival2,rval1,rval2,dval1,dval2)
	    if( ntot .gt. ndim ) goto 97

	    idiff = 0

	    blcheck = bcheck		!local check
	    call handle_exception(text,ntime,blcheck)

	    if( .not. blcheck ) then
	      read(1)
	      read(2)
	      if( bverbose ) write(6,*) nrec,nh,nv,nt,idiff,trim(text)
	      cycle
	    end if
	
	    bshowdiff = .not. bnodiff

	    call read_data_record(nt,ntot)

	    if( bcheck ) then
	      call check_val(dtime,nrec,nt,nh,nv,idiff,diff)
	      bshowdiff = bshowdiff .and. idiff > 0
	    end if

	    if( text == 'ipv' ) call save_int(nh,ival1,nipv,ipv)
	    if( text == 'ipev' ) call save_int(nh,ival1,nipev,ipev)

	    if( bshowdiff .or. bverbose ) then
	      call info(nt,nh,nv,ipv,ipev,text)
	    end if

	    if( bsummary ) then
	      rdiff_max = max(rdiff_max,diff)
	    else if( idiff > 0 .or. bverbose ) then
	      if( bheader ) then
		bheader = .false.
	        write(6,'(a)') '    nrec      nh      nv      nt' //
     +				'     idiff          max-diff var'
	      end if
	      write(6,2000) nrec,nh,nv,nt,idiff,diff,' ',trim(text)
	    end if
	    idiff_tot = idiff_tot + idiff
	    idiff_rec = idiff_rec + idiff

	  end do

	  if( idiff_rec > 0 ) ntrerr = ntrerr + 1
	  if( bstop .and. idiff_tot > 0 ) exit
	  if( bverbose ) write(6,*) 'nrecs checked: ',nrec
	end do

    9	continue

	close(1)
	close(2)

	if( .not. bsilent ) then
	  write(6,*) 'time records read: ',ntime
     +			,'  differences found: ',idiff_tot
	end if

	idiff_end = idiff_tot
	idiff_end = ntrerr

 2000	format(4i8,i10,f18.6,2a)
	return
   99	continue
	write(6,*) 'times: ',dtime1,dtime2
	stop 'error stop check_debug: time mismatch'
   98	continue
	if( nh1 == 0 .and. nv1 == 0 .and. nt1 == 0 ) then
	  write(6,*) 'file 1 has finished data records for time step'
	  stop 'error stop check_debug: not enough records'
	end if
	if( nh2 == 0 .and. nv2 == 0 .and. nt2 == 0 ) then
	  write(6,*) 'file 2 has finished data records for time step'
	  stop 'error stop check_debug: not enough records'
	end if
	write(6,*) 'params nh1/nh2: ',nh1,nh2
	write(6,*) 'params nv1/nv2: ',nv1,nv2
	write(6,*) 'params nt1/nt2: ',nt1,nt2
	stop 'error stop check_debug: size or type mismatch'
   97	continue
	write(6,*) 'params: ',dtime,nrec,ntot,ndim
	stop 'error stop check_debug: dimension'
   96	continue
	write(6,*) 'text: ',trim(text1),' - ',trim(text2)
	stop 'error stop check_debug: text mismatch'
	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine read_time_header(dtime,ntime,ierr)

	use mod_shympi_debug

	implicit none

	double precision dtime
	integer ntime
	integer ierr

	double precision dtime1,dtime2

	  read(1,end=9) dtime1
	  read(2,end=9) dtime2
	  if( dtime1 .ne. dtime2 ) goto 99
	  dtime = dtime1
	  if( .not. bquiet ) write(6,*) 'time = ',dtime
	  ntime = ntime + 1

	  if( bverbose ) then
	    write(6,*) '       irec          nh          nv' //
     +			'        type        diff name'
	  end if

	ierr = 0
	return

    9	continue
	ierr = -1
	return
   99	continue
	write(6,*) 'times: ',dtime1,dtime2
	stop 'error stop check_debug: time mismatch'
	end 

!*******************************************************************

	subroutine read_data_header(nh,nv,nt,ntot,nrec,text,ierr)

	use mod_shympi_debug

	implicit none

	integer nh,nv,nt
	integer ntot,nrec
	character*(*) text
	integer ierr

	integer nh1,nv1,nt1
	integer nh2,nv2,nt2
	character*80 text1,text2

	ierr = 0

	    read(1) nh1,nv1,nt1
	    read(2) nh2,nv2,nt2
	    if( nh1 .ne. nh2 ) goto 98
	    if( nv1 .ne. nv2 ) goto 98
	    if( nt1 .ne. nt2 ) goto 98
	    nh = nh1
	    nv = nv1
	    nt = nt1
	    ntot = nh*nv

	    if( nt .eq. 0 ) goto 9
	    nrec = nrec + 1

	    !read(1,end=9) text1
	    !read(2,end=9) text2
	    read(1) text1
	    read(2) text2
	    if( text1 .ne. text2 ) goto 96
	    text = text1

	return
   9	continue
	ierr = -1
	return
   98	continue
	if( nh1 == 0 .and. nv1 == 0 .and. nt1 == 0 ) then
	  write(6,*) 'file 1 has finished data records for time step'
	  stop 'error stop check_debug: not enough records'
	end if
	if( nh2 == 0 .and. nv2 == 0 .and. nt2 == 0 ) then
	  write(6,*) 'file 2 has finished data records for time step'
	  stop 'error stop check_debug: not enough records'
	end if
	write(6,*) 'params nh1/nh2: ',nh1,nh2
	write(6,*) 'params nv1/nv2: ',nv1,nv2
	write(6,*) 'params nt1/nt2: ',nt1,nt2
	stop 'error stop check_debug: size or type mismatch'
   96	continue
	write(6,*) 'text: ',trim(text1),' - ',trim(text2)
	stop 'error stop check_debug: text mismatch'
	end

!*******************************************************************

	subroutine read_data_record(nt,ntot)

	use mod_shympi_debug

	implicit none

	integer nt,ntot

	integer i

	if( nt == 1 ) then			!integer
	  read(1) (ival1(i),i=1,ntot)
	  read(2) (ival2(i),i=1,ntot)
	else if( nt == 2 ) then			!real
	  read(1) (rval1(i),i=1,ntot)
	  read(2) (rval2(i),i=1,ntot)
	else if( nt == 3 ) then			!double
	  read(1) (dval1(i),i=1,ntot)
	  read(2) (dval2(i),i=1,ntot)
	else
	  write(6,*) 'cannot handle nt = ',nt
	  stop 'error stop: nt'
	end if

	end

!*******************************************************************

	subroutine check_val(dtime,nrec,nt,nh,nv,idiff,diff)

	use mod_shympi_debug

	implicit none

	double precision dtime
	integer nrec,nt
	integer nh,nv,idiff
	double precision diff

	if( nt == 1 ) then			!integer
	  call check_ival(dtime,nrec,nh,nv,ival1,ival2,idiff,diff)
	else if( nt == 2 ) then			!real
	  call check_rval(dtime,nrec,nh,nv,rval1,rval2,idiff,diff)
	else if( nt == 3 ) then			!double
	  call check_dval(dtime,nrec,nh,nv,dval1,dval2,idiff,diff)
	else
	  write(6,*) 'cannot handle nt = ',nt
	  stop 'error stop: nt'
	end if

	end

!*******************************************************************

	subroutine info(nt,nh,nv,ipv,ipev,text)

	use mod_shympi_debug

	implicit none

	integer nt
	integer nh,nv
	integer ipv(nipv),ipev(nipev)
	character*(*) text

	if( nt == 1 ) then			!integer
	  call i_info(nh,nv,ival1,ival2,ipv,ipev,text)
	else if( nt == 2 ) then			!real
	  call r_info(nh,nv,rval1,rval2,ipv,ipev,text)
	else if( nt == 3 ) then			!double
	  call d_info(nh,nv,dval1,dval2,ipv,ipev,text)
	else
	  write(6,*) 'cannot handle nt = ',nt
	  stop 'error stop: nt'
	end if

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine check_dval(dtime,nrec,nh,nv,val1,val2,idiff,diff)

	use mod_shympi_debug

	implicit none

	double precision dtime
	integer nrec
	integer nh,nv,idiff
	double precision val1(nh*nv)
	double precision val2(nh*nv)
	double precision diff

	integer i,k,l,ntot

	idiff = 0
	ntot = nh*nv

	do i=1,ntot
	  if( abs(val1(i)-val2(i)) > maxdiff ) then
	    k = 1 + (i-1)/nv
	    l = 1 + mod(i-1,nv)
	    !if( idiff == 0 ) write(77,*) 'check_dval...'
	    !write(77,*) dtime,nrec,k,l,val1(i),val2(i)
	    idiff = idiff + 1
	  end if
	end do

	diff = 0.
	if( idiff > 0 ) then
	  diff = maxval(abs(val1-val2))
	end if

	end

!*******************************************************************

	subroutine check_rval(dtime,nrec,nh,nv,val1,val2,idiff,diff)

	use mod_shympi_debug

	implicit none

	double precision dtime
	integer nrec
	integer nh,nv,idiff
	real val1(nh*nv)
	real val2(nh*nv)
	double precision diff

	integer i,k,l,ntot

	idiff = 0
	ntot = nh*nv

	do i=1,ntot
	  if( abs(val1(i)-val2(i)) > maxdiff ) then
	    k = 1 + (i-1)/nv
	    l = 1 + mod(i-1,nv)
	    !if( idiff == 0 ) write(77,*) 'check_rval...'
	    !write(77,*) dtime,nrec,k,l,val1(i),val2(i)
	    idiff = idiff + 1
	  end if
	end do

	diff = 0.
	if( idiff > 0 ) then
	  diff = maxval(abs(val1-val2))
	end if

	end

!*******************************************************************

	subroutine check_ival(dtime,nrec,nh,nv,val1,val2,idiff,diff)

	implicit none

	double precision dtime
	integer nrec
	integer nh,nv,idiff
	integer val1(nh*nv)
	integer val2(nh*nv)
	double precision diff

	integer i,k,l,ntot

	idiff = 0
	ntot = nh*nv

	do i=1,ntot
	  if( val1(i) .ne. val2(i) ) then
	    k = 1 + (i-1)/nv
	    l = 1 + mod(i-1,nv)
	    !if( idiff == 0 ) write(77,*) 'check_ival...'
	    !write(77,*) dtime,nrec,k,l,val1(i),val2(i)
	    idiff = idiff + 1
	  end if
	end do

	diff = 0.
	if( idiff > 0 ) then
	  diff = maxval(abs(val1-val2))
	end if

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine d_info(nh,nv,val1,val2,ipv,ipev,text)

	use mod_shympi_debug

	implicit none

	integer nh,nv
	double precision val1(nh*nv)
	double precision val2(nh*nv)
	!integer ipv(:),ipev(:)
	integer ipv(nipv),ipev(nipev)
	character*(*) text

	logical belem
	integer i,ih,iv,ierr,j
	integer ipvv(nh)
	double precision maxdif,diff
	character*80 textk,texte,text1,text2

	if( nh == size(ipv) ) then
	  belem = .false.
	  ipvv = ipv
	else if( nh == size(ipev) ) then
	  belem = .true.
	  ipvv = ipev
	else
	  write(6,*) nh,size(ipv),size(ipev)
	  stop 'error stop r_info: unknown nh'
	end if

	textk = '         irec       k       l    kext'
	texte = '         irec      ie       l   ieext'
	text2 = '              val1              val2'
	text1 = textk
	if( belem ) text1 = texte

	write(6,*) 'differences found reading ',trim(text)
	write(6,*) trim(text1) // trim(text2)

	ierr = 0
	maxdif = 0.
	ifill = 0
	do i=1,nh*nv
	  if( val1(i) /= val2(i) ) then
	    diff = abs(val1(i)-val2(i))
	    call insert_diffs(i,diff)
	    maxdif = max(maxdif,diff)
	    ierr = ierr + 1
	  end if
	end do

	do j=1,ifill
	  i = index(j)
	  iv = 1 + mod((i-1),nv)
	  ih = 1 + (i-1)/nv
	  write(6,1000) 'diff: ',i,ih,iv,ipvv(ih),val1(i),val2(i)
	end do
	write(6,*) 'maximum difference: ',maxdif

 1000	format(a,4i8,2f18.6)
	end

!*******************************************************************

	subroutine r_info(nh,nv,val1,val2,ipv,ipev,text)

	use mod_shympi_debug

	implicit none

	integer nh,nv
	real val1(nh*nv)
	real val2(nh*nv)
	!integer ipv(:),ipev(:)
	integer ipv(nipv),ipev(nipev)
	character*(*) text

	logical belem
	integer i,ih,iv,ierr,j
	integer ipvv(nh)
	double precision maxdif,diff
	character*80 textk,texte,text1,text2

	if( nh == size(ipv) ) then
	  belem = .false.
	  ipvv = ipv
	else if( nh == size(ipev) ) then
	  belem = .true.
	  ipvv = ipev
	else
	  write(6,*) nh,size(ipv),size(ipev)
	  stop 'error stop r_info: unknown nh'
	end if

	textk = '         irec       k       l    kext'
	texte = '         irec      ie       l   ieext'
	text2 = '              val1              val2'
	text1 = textk
	if( belem ) text1 = texte

	write(6,*) 'differences found reading ',trim(text)
	write(6,*) trim(text1) // trim(text2)

	ierr = 0
	maxdif = 0.
	ifill = 0
	do i=1,nh*nv
	  if( val1(i) /= val2(i) ) then
	    diff = abs(val1(i)-val2(i))
	    call insert_diffs(i,diff)
	    maxdif = max(maxdif,diff)
	    ierr = ierr + 1
	  end if
	end do

	do j=1,ifill
	  i = index(j)
	  iv = 1 + mod((i-1),nv)
	  ih = 1 + (i-1)/nv
	  write(6,1000) 'diff: ',i,ih,iv,ipvv(ih),val1(i),val2(i)
	end do
	write(6,*) 'maximum difference: ',maxdif

 1000	format(a,4i8,2f18.6)
	end

!*******************************************************************

	subroutine i_info(nh,nv,val1,val2,ipv,ipev,text)

	use mod_shympi_debug

	implicit none

	integer nh,nv
	integer val1(nh*nv)
	integer val2(nh*nv)
	!integer ipv(:),ipev(:)
	integer ipv(nipv),ipev(nipev)
	character*(*) text
	!integer, optional :: iunit

	logical belem
	integer i,ih,iv,iu,ierr
	integer ipvv(nh)
	character*80 textk,texte,text1,text2

	iu = 0
	!if( present(iunit) ) iu = iunit

	if( nh == size(ipv) ) then
	  ipvv = ipv
	else if( nh == size(ipev) ) then
	  ipvv = ipev
	else
	  write(6,*) nh,size(ipv),size(ipev)
	  stop 'error stop i_info: unknown nh'
	end if

	textk = '         irec       k       l    kext'
	texte = '         irec      ie       l   ieext'
	text2 = '              val1              val2'
	text1 = textk
	if( belem ) text1 = texte

	write(6,*) 'differences found reading ',trim(text)
	write(6,*) trim(text1) // trim(text2)
	if( iu > 0 ) write(iu,*) trim(text1) // trim(text2)

	do i=1,nh*nv
	  iv = 1 + mod((i-1),nv)
	  ih = 1 + (i-1)/nv
	  if( iu > 0 ) then
	    write(iu,1000) 'diff: ',i,ih,iv,ipvv(ih),val1(i),val2(i)
	  else if( val1(i) /= val2(i) ) then
	    ierr = ierr + 1
	    if( imax > 0 .and. ierr > imax ) cycle
	    write(6,1000) 'diff: ',i,ih,iv,ipvv(ih),val1(i),val2(i)
	  end if
	end do

 1000	format(a,4i8,2i18)
 2000	format(4i8,2i18)
	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine time_balance

! finds first time record equal in both files

	implicit none

	integer iu,nrec
	integer ierr,ierr1,ierr2
	double precision dtime,dtime1,dtime2,dtime_high

	call peek_time_record(1,dtime1,ierr1)
	call peek_time_record(2,dtime2,ierr2)
	if( ierr1 /= 0 .or. ierr2 /= 0 ) goto 99

	if( dtime1 == dtime2 ) then	!same time in initial record
	  write(6,*) 'time_balance: starting time found ',dtime1
	  return
	else if( dtime1 < dtime2 ) then
	  iu = 1
	  dtime_high = dtime2
	else
	  iu = 2
	  dtime_high = dtime1
	end if

	nrec = 0

	do while(.true.)
	  call skip_time_record(iu)
	  call peek_time_record(iu,dtime,ierr)
	  nrec = nrec + 1
	  if( ierr /= 0 ) goto 98
	  if( dtime > dtime_high ) goto 98
	  if( dtime == dtime_high ) exit
	end do

	write(6,*) 'time_balance: starting time found ',nrec,dtime

	return
   99	continue
	write(6,*) 'cannot read first time record...'
	stop 'error stop time_balance: first record'
   98	continue
	write(6,*) 'cannot find corresponding time records...'
	write(6,*) iu,nrec
	write(6,*) dtime1,dtime2
	write(6,*) dtime,dtime_high
	stop 'error stop time_balance: no corresponding time record'
	end

!*******************************************************************

	subroutine peek_time_record(iunit,dtime,ierr)

	implicit none

	integer iunit
	double precision dtime
	integer ierr

	read(iunit,end=9) dtime
	backspace(iunit)

	ierr = 0
	!write(6,*) ierr,dtime
	return
    9	continue
	ierr = -1
	return
	end

!*******************************************************************

	subroutine skip_time_record(iunit)

	implicit none

	integer iunit

	double precision dtime
	integer nh,nv,nt
	integer nrec

	nrec = 0
	read(iunit) dtime

	do while(.true.)
	    read(iunit) nh,nv,nt
	    if( nt == 0 ) exit
	    nrec = nrec + 1
	    read(iunit)
	    read(iunit)
	end do

	!write(6,*) 'skipping records: ',nrec,dtime

	end

!*******************************************************************

	subroutine read_header(iunit,nipv,nipev,ipv,ipev)

	implicit none

	integer iunit
	integer nipv,nipev
	integer, allocatable :: ipv(:),ipev(:)

	integer nrec
	integer nh,nv,nt
	character*80 text
	double precision dtime
	integer, allocatable :: ival(:)

        INTERFACE
          subroutine save_int(n,ival,niv,iv)
          integer n,niv
          integer ival(n)
          integer, allocatable :: iv(:)
          END subroutine
        END INTERFACE

	nrec = 0
	read(iunit,end=98) dtime
	if( dtime /= -1. ) goto 99

	do while(.true.)
	  read(iunit) nh,nv,nt
	  if( nt == 0 ) exit
	  nrec = nrec + 1
	  read(iunit) text
	  if( nt == 1 ) then
	    if( allocated(ival) ) deallocate(ival)
	    allocate(ival(nh*nv))
	    read(iunit) ival
	    if( text == 'ipv' ) then
	      if( allocated(ipv) ) then
		if( any(ival/=ipv) ) goto 99
	      else
		call save_int(nh,ival,nipv,ipv)
	      end if
	    end if
	    if( text == 'ipev' ) then
	      if( allocated(ipev) ) then
		if( any(ival/=ipev) ) goto 99
	      else
		call save_int(nh,ival,nipev,ipev)
	      end if
	    end if
	  else
	    read(iunit)
	  end if
	end do

	return
   98	continue
	stop 'error stop read_header: no header found'
	return
   99	continue
	stop 'error stop read_header: headers are different'
	end

!*******************************************************************

	subroutine show_extra(what,k,text,nh,nv,ntot,rval1,rval2)

! debug section - normally commented out

	implicit none

	character*(*) what,text
	integer k,nh,nv,ntot
	real rval1(ntot),rval2(ntot)

	integer iu
	integer is,ie,n,i

	iu = 666

	if( what /= text ) return

	write(iu,*) what,k,nh,nv
	is=1+(k-1)*nv
	ie=is+nv-1

	n=0
	do i=is,ie
	  n = n + 1
	  write(iu,*) n,i,rval1(i),rval2(i)
	end do

	end

!*******************************************************************

	subroutine handle_exception(text,ntime,blcheck)

	implicit none

	character*(*) text
	integer ntime
	logical blcheck

	logical bignore

	if( ntime > 2 ) return		!only do for first data record

	bignore = .false.

	if( trim(text) == 'zov' ) bignore = .true.
	if( trim(text) == 'fxv' ) bignore = .true.
	if( trim(text) == 'fyv' ) bignore = .true.

	if( bignore ) then
	  blcheck = .false.
	  write(6,*) 'ignoring ',trim(text),' in time record ',ntime
	end if

	end

!*******************************************************************

	subroutine insert_diffs(i,diff)

	use mod_shympi_debug

	implicit none

	integer i
	double precision diff

	integer jmin,j
	!integer :: ifill
	!integer :: index(imax)
	!double precision :: diffs(imax)

	if( diff == 0 ) return

	if( ifill < imax ) then
	  ifill = ifill + 1
	  diffs(ifill) = diff
	  index(ifill) = i
	  diffmin = minval(diffs(1:ifill))
	else if( diff > diffmin ) then
	  jmin = 0
	  do j=1,imax
	    if( diffs(j) == diffmin ) jmin = j
	  end do
	  if( jmin == 0 ) stop 'error stop: jmin == 0'
	  diffs(jmin) = diff
	  index(jmin) = i
	  diffmin = minval(diffs)
	end if
	
	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine save_int(n,ival,niv,iv)

	implicit none

	integer n,niv
	integer ival(n)
	integer, allocatable :: iv(:)

	allocate(iv(n))
	iv = ival
	niv = n

	end

!*******************************************************************

	subroutine allocate_arrays(nsize,ndim
     +			,ival1,ival2,rval1,rval2,dval1,dval2)

	implicit none

	integer nsize,ndim
	integer, allocatable :: ival1(:),ival2(:)
	real, allocatable :: rval1(:),rval2(:)
	double precision, allocatable :: dval1(:),dval2(:)

	if( nsize <= ndim ) return

	!write(6,*) 'allocating arrays: ',nsize,ndim

	if( allocated(ival1) ) deallocate(ival1)
	if( allocated(ival2) ) deallocate(ival2)
	if( allocated(rval1) ) deallocate(rval1)
	if( allocated(rval2) ) deallocate(rval2)
	if( allocated(dval1) ) deallocate(dval1)
	if( allocated(dval2) ) deallocate(dval2)

	ndim = nsize

	allocate(ival1(ndim))
	allocate(ival2(ndim))
	allocate(rval1(ndim))
	allocate(rval2(ndim))
	allocate(dval1(ndim))
	allocate(dval2(ndim))

	end

!*******************************************************************

        subroutine check_shympi_debug_init

        use clo
	use mod_shympi_debug

        implicit none

	logical baux
        character*80 version

	version = '2.0'

	call clo_init('check_shympi_debug','dbg-file(s)',trim(version))

        call clo_add_info('checks shyfem debug files')

        call clo_add_sep('general options:')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('quiet',.false.,'be quiet')
        call clo_add_option('nodiff',.false.,'do not show differences')
        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('nostop',.false.,'do not stop at error')
        call clo_add_option('summary',.false.,'do only summary')
        call clo_add_option('balance',.false.,'balance time records')
        call clo_add_option('maxdiff',0.,'maximum tolerated difference')

	call clo_parse_options

        call clo_get_option('silent',bsilent)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('nodiff',bnodiff)
        call clo_get_option('verbose',bverbose)
        call clo_get_option('nostop',baux)
        call clo_get_option('summary',bsummary)
        call clo_get_option('balance',bbalance)
        call clo_get_option('maxdiff',maxdiff)

        if( baux ) bstop = .false.
        if( bsilent ) bquiet = .true.
        if( bquiet ) bverbose = .false.

        end

!*******************************************************************

