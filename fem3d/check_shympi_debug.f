
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

!**************************************************************************

!==========================================================================
	module mod_shympi_debug
!==========================================================================

	implicit none

	logical, save :: bstop = .true.		!stop on error
	logical, save :: bcheck = .true.	!check for differences

	logical, save :: bverbose = .false.	!be verbose
	logical, save :: bquiet = .false.	!be verbose
	logical, save :: bsilent = .false.	!be verbose

	integer, parameter :: imax = 20		!number of errors shown

!==========================================================================
	end module mod_shympi_debug
!==========================================================================

	program check_shympi_debug

        use clo

	implicit none

	integer nc,ierr

	call check_shympi_debug_init

	nc = clo_number_of_files()

	if( nc == 0 ) then
	  call clo_usage
	else if( nc == 1 ) then
	  call read_file
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

	subroutine read_file

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

	if( ios /= 0 ) stop 'error opening file'

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
	    if( nh == 0 ) exit
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
	  subroutine d_info(nh,nv,val1,val2,ipv,ipev,text)
	  integer nh,nv
	  double precision val1(nh*nv)
	  double precision val2(nh*nv)
	  integer ipv(:),ipev(:)
	  character*(*) text
	  END subroutine
	END INTERFACE

	INTERFACE 
	  subroutine r_info(nh,nv,val1,val2,ipv,ipev,text)
	  integer nh,nv
	  real val1(nh*nv)
	  real val2(nh*nv)
	  integer ipv(:),ipev(:)
	  character*(*) text
	  END subroutine
	END INTERFACE

	INTERFACE 
	  subroutine i_info(nh,nv,val1,val2,ipv,ipev,text,iunit)
	  integer nh,nv
	  integer val1(nh*nv)
	  integer val2(nh*nv)
	  integer ipv(:),ipev(:)
	  character*(*) text
	  integer, optional :: iunit
	  END subroutine
	END INTERFACE

	integer idiff_end

	integer, save :: ndim = 0

	character*60 name_one,name_two
	character*80 text1,text2,text
	integer nt1,nt2,nt
	integer nh1,nh2,nh
	integer nv1,nv2,nv
	integer nrec,ntot,ntime
	integer i,idiff,idiff_tot
	integer nc
	integer nipv,nipev
	integer ios
	integer, allocatable :: ipv(:),ipev(:)
	real rdiff
	double precision dtime,dtime1,dtime2
	double precision diff

	integer, allocatable :: ival1(:),ival2(:)
	real, allocatable :: rval1(:),rval2(:)
	double precision, allocatable :: dval1(:),dval2(:)

	call clo_get_file(1,name_one)
	call clo_get_file(2,name_two)

	open(1,file=name_one,status='old',form='unformatted',iostat=ios)
	if( ios /= 0 ) stop 'error opening file 1'
	open(2,file=name_two,status='old',form='unformatted',iostat=ios)
	if( ios /= 0 ) stop 'error opening file 2'

	if( .not. bquiet ) then
	  write(6,*) 'file 1: ',trim(name_one)
	  write(6,*) 'file 2: ',trim(name_two)
	end if

	idiff_tot = 0
	ntime = 0
	nipv = 0
	nipev = 0

	do while(.true.)

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

	  nrec = 0
	  do while(.true.)
	    read(1) nh1,nv1,nt1
	    read(2) nh2,nv2,nt2
	    if( nh1 .ne. nh2 ) goto 98
	    if( nv1 .ne. nv2 ) goto 98
	    if( nt1 .ne. nt2 ) goto 98
	    nh = nh1
	    nv = nv1
	    nt = nt1
	    ntot = nh*nv

	    if( nt .eq. 0 ) exit
	    nrec = nrec + 1

	    read(1,end=9) text1
	    read(2,end=9) text2
	    if( text1 .ne. text2 ) goto 96
	    text = text1

	    call allocate_arrays(ntot,ndim
     +			,ival1,ival2,rval1,rval2,dval1,dval2)
	    if( ntot .gt. ndim ) goto 97

	    idiff = 0

	    if( .not. bcheck ) then
	      read(1)
	      read(2)
	      if( bverbose ) write(6,*) nrec,nh,nv,nt,idiff,trim(text)
	      cycle
	    end if

	    if( nt == 1 ) then			!integer
	      read(1) (ival1(i),i=1,ntot)
	      read(2) (ival2(i),i=1,ntot)
	      if( bcheck ) then
	        call check_ival(dtime,nrec,nh,nv,ival1,ival2,idiff,diff)
	      end if
	      if( text == 'ipv' ) call save_int(nh,ival1,nipv,ipv)
	      if( text == 'ipev' ) call save_int(nh,ival1,nipev,ipev)
	      if( idiff > 0 .and. bverbose ) then
	        call i_info(nh,nv,ival1,ival2,ipv,ipev,text)
	      end if
	    else if( nt == 2 ) then		!real
	      read(1) (rval1(i),i=1,ntot)
	      read(2) (rval2(i),i=1,ntot)
	      if( bcheck ) then
	        call check_rval(dtime,nrec,nh,nv,rval1,rval2,idiff,diff)
	      end if
	      if( idiff > 0 .or. bverbose ) then
	        call r_info(nh,nv,rval1,rval2,ipv,ipev,text)
	      end if
	    else if( nt == 3 ) then		!double
	      read(1) (dval1(i),i=1,ntot)
	      read(2) (dval2(i),i=1,ntot)
	      if( bcheck ) then
	        call check_dval(dtime,nrec,nh,nv,dval1,dval2,idiff,diff)
	      end if
	      if( idiff > 0 .or. bverbose ) then
		call d_info(nh,nv,dval1,dval2,ipv,ipev,text)
	      end if
	    else
	      write(6,*) 'cannot handle nt = ',nt
	      stop 'error stop: nt'
	    end if

	    if( idiff > 0 .or. bverbose ) then
	      write(6,2000) nrec,nh,nv,nt,idiff,diff,' ',trim(text)
	    end if
	    idiff_tot = idiff_tot + idiff

	  end do

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

	subroutine check_dval(dtime,nrec,nh,nv,val1,val2,idiff,diff)

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
	  if( val1(i) .ne. val2(i) ) then
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
	  if( val1(i) .ne. val2(i) ) then
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
	integer ipv(:),ipev(:)
	character*(*) text

	logical belem
	integer i,ih,iv,ierr
	integer ipvv(nh)
	double precision maxdif
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
	do i=1,nh*nv
	  if( val1(i) /= val2(i) ) then
	    maxdif = max(maxdif,abs(val1(i)-val2(i)))
	    ierr = ierr + 1
	    if( imax > 0 .and. ierr > imax ) cycle
	    iv = 1 + mod((i-1),nv)
	    ih = 1 + (i-1)/nv
	    write(6,1000) 'diff: ',i,ih,iv,ipvv(ih),val1(i),val2(i)
	  end if
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
	integer ipv(:),ipev(:)
	character*(*) text

	logical belem
	integer i,ih,iv,ierr
	integer ipvv(nh)
	double precision maxdif
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
	do i=1,nh*nv
	  if( val1(i) /= val2(i) ) then
	    maxdif = max(maxdif,abs(val1(i)-val2(i)))
	    ierr = ierr + 1
	    if( imax > 0 .and. ierr > imax ) cycle
	    iv = 1 + mod((i-1),nv)
	    ih = 1 + (i-1)/nv
	    write(6,1000) 'diff: ',i,ih,iv,ipvv(ih),val1(i),val2(i)
	  end if
	end do
	write(6,*) 'maximum difference: ',maxdif

 1000	format(a,4i8,2f18.6)
	end

!*******************************************************************

	subroutine i_info(nh,nv,val1,val2,ipv,ipev,text,iunit)

	use mod_shympi_debug

	implicit none

	integer nh,nv
	integer val1(nh*nv)
	integer val2(nh*nv)
	integer ipv(:),ipev(:)
	character*(*) text
	integer, optional :: iunit

	logical belem
	integer i,ih,iv,iu,ierr
	integer ipvv(nh)
	character*80 textk,texte,text1,text2

	iu = 0
	if( present(iunit) ) iu = iunit

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
        call clo_add_option('quiet',.false.,'be quiet')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('nostop',.false.,'do not stop at error')

	call clo_parse_options

        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)
        call clo_get_option('verbose',bverbose)
        call clo_get_option('nostop',baux)

        if( baux ) bstop = .false.
        if( bsilent ) bquiet = .true.

        end

!*******************************************************************

