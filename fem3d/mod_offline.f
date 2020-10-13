
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013-2015,2017-2020  Georg Umgiesser
!    Copyright (C) 2015  Christian Ferrarin
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

! routines for offline data handling
!
! revision log :
!
! 13.06.2013	ggu	new routines written from scratch
! 17.06.2013	ggu	eliminated compiler warnings
! 25.03.2014	ggu	new offline (for T/S)
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 01.04.2015	ggu	changed VERS_7_1_7
! 06.05.2015	ccf	write offline to .off file
! 06.05.2015	ccf	read offline from offlin file in section name
! 21.05.2015	ggu	changed VERS_7_1_11
! 10.07.2015	ggu	changed VERS_7_1_50
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 05.11.2015	ggu	revisited and checked
! 09.11.2015	ggu	changed VERS_7_3_13
! 16.11.2015	ggu	changed VERS_7_3_14
! 29.03.2017	ggu	bug fix - input file opened on unit 1
! 05.12.2017	ggu	changed VERS_7_5_39
! 12.11.2018	ggu	linear arrays introduced
! 18.12.2018	ggu	changed VERS_7_5_52
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 04.07.2019	ggu	solved problem for vertical velocity (1 index off)
! 17.02.2020	ggu	femtime eliminated
! 28.04.2020	ggu	restructured and taken out of suboff.f
! 09.10.2020	ggu	added comments and error checking
! 13.10.2020	ggu	default values for T/S introduced
!
! contents :
!
! mod_offline_init()	initializes offline module
! off_write_record()	writes one record to offline file
! off_read_record()	reads one record from offline file
! off_read_header()	reads header and sets nlv
! off_peek_header()	peeks into header and sets nlv
! off_init_vertical()	initializes vertical indices
! off_check_vertical()	checks if vertical indices have been initialized
! off_error()		exits with error
! off_next_record()	checks info on next record
!
! calling sequence for writing:
!
!	mod_offline_init()
!	do
!	  off_write_record()
!	end do
!
! calling sequence for reading:
!
!	off_peek_header()
!	mod_offline_init()
!	do
!	  off_read_record()
!	end do
!
!****************************************************************

!==================================================================
	module mod_offline
!==================================================================

	implicit none

	integer, parameter :: nintp = 4		!2 (linear) or 4 (cubic)

	integer, save :: ioffline = 0
	integer, save :: idtoff,itmoff,itoff
	integer, save :: iwhat			!0 (none), 1 (write), 2 (read)
	integer, save :: iread			!what is in file
	integer, save :: iuoff			!unit to read/write
	integer, save :: icall = 0
	logical, save :: bfirst = .true.
	logical, save :: bdebug = .false.

	integer, save :: nkn_off = 0
	integer, save :: nel_off = 0
	integer, save :: nlv_off = 0

	double precision, save :: dtr = 0.
	integer, save :: time(nintp)

	integer, save :: idef = 0		!use default values for T/S
	real, save :: tdef = 0
	real, save :: sdef = 0

	double precision, save, allocatable :: ut(:,:,:)
	double precision, save, allocatable :: vt(:,:,:)
	double precision, save, allocatable :: ze(:,:,:)
	double precision, save, allocatable :: wn(:,:,:)
	double precision, save, allocatable :: zn(:,:)
	double precision, save, allocatable :: sn(:,:,:)
	double precision, save, allocatable :: tn(:,:,:)

	integer, save, allocatable :: ile(:)
	integer, save, allocatable :: ilk(:)
	logical, save :: bvinit = .false.	!il files initialized

!==================================================================
	contains
!==================================================================

	subroutine mod_offline_init(nk,ne,nl)

! initializes offline module

	integer nk,ne,nl

        if( nk == nkn_off .and. ne == nel_off
     +          .and. nl == nlv_off ) return

        if( nk > 0 .or. ne > 0 .or. nl > 0 ) then
          if( nk == 0 .or. ne == 0 .or. nl == 0 ) then
            write(6,*) 'nk,ne,nl: ',nk,ne,nl
            stop 'error stop mod_offline_init: incompatible parameters'
          end if
        end if

	if( nkn_off > 0 ) then
	  deallocate(ut)
	  deallocate(vt)
	  deallocate(ze)
	  deallocate(wn)
	  deallocate(zn)
	  deallocate(sn)
	  deallocate(tn)
	  deallocate(ile)
	  deallocate(ilk)
	end if

	nkn_off = nk
	nel_off = ne
	nlv_off = nl

	allocate(ut(nl,ne,nintp))
	allocate(vt(nl,ne,nintp))
	allocate(ze(3,ne,nintp))
	allocate(wn(0:nl,nk,nintp))
	allocate(zn(nk,nintp))
	allocate(sn(nl,nk,nintp))
	allocate(tn(nl,nk,nintp))
	allocate(ile(ne))
	allocate(ilk(nk))

	wn = 0.

	end subroutine mod_offline_init

!==================================================================
	end module mod_offline
!==================================================================

!****************************************************************
	
	subroutine off_write_record(iu,it)

! writes one record to offline file

	use mod_offline

	implicit none

	integer iu,it

	integer ie,ii,k,i
	integer nlin,nlink,nline
	integer iunit
	integer nkn,nel,nlv,nlvdi
	
	double precision, save, allocatable :: wnaux(:,:)
	double precision, save, allocatable :: rlin(:)

	if( .not. bvinit ) then
	  stop 'error stop off_write_record: bvinit is false'
	else if( nkn_off <= 0  ) then
	  stop 'error stop off_write_record: offline not initialized'
	end if

!----------------------------------------------------------
! initialize
!----------------------------------------------------------

	iunit = iu
	nkn = nkn_off
	nel = nel_off
	nlv = nlv_off
	nlvdi = nlv

!----------------------------------------------------------
! set up auxiliary arrays
!----------------------------------------------------------

        call count_linear(nlvdi,nkn,1,ilk,nlink)
        call count_linear(nlvdi,nel,1,ile,nline)
	nlin = max(nlink,nline)
        if( .not. allocated(rlin) ) allocate(rlin(nlin))
        if( .not. allocated(wnaux) ) allocate(wnaux(nlvdi,nkn))

!----------------------------------------------------------
! write header
!----------------------------------------------------------

	write(iunit) it,nkn,nel,3
	write(iunit) (ile(ie),ie=1,nel)
	write(iunit) (ilk(k),k=1,nkn)

!----------------------------------------------------------
! write currents
!----------------------------------------------------------

        nlin = nline
        call dvals2linear(nlvdi,nel,1,ile,ut,rlin,nlin)
        write(iunit) (rlin(i),i=1,nlin)
        call dvals2linear(nlvdi,nel,1,ile,vt,rlin,nlin)
        write(iunit) (rlin(i),i=1,nlin)

	!write(iunit) ((ut(l,ie,1),l=1,ile(ie)),ie=1,nel)
	!write(iunit) ((vt(l,ie,1),l=1,ile(ie)),ie=1,nel)

!----------------------------------------------------------
! write water levels and vertical velocities
!----------------------------------------------------------

        nlin = nlink
	write(iunit) ((ze(ii,ie,1),ii=1,3),ie=1,nel)
        wnaux(1:nlvdi,:) = wn(1:nlvdi,:,1)
        call dvals2linear(nlvdi,nkn,1,ilk,wnaux,rlin,nlin)
        write(iunit) (rlin(i),i=1,nlin)
	!write(iunit) ((wn(l,k,1),l=1,ilk(k)),k=1,nkn)
	write(iunit) (zn(k,1),k=1,nkn)

!----------------------------------------------------------
! write T/S
!----------------------------------------------------------

        nlin = nlink
        call dvals2linear(nlvdi,nkn,1,ilk,sn,rlin,nlin)
        write(iunit) (rlin(i),i=1,nlin)
        call dvals2linear(nlvdi,nkn,1,ilk,tn,rlin,nlin)
        write(iunit) (rlin(i),i=1,nlin)

	!write(iunit) ((sn(l,k,1),l=1,ilk(k)),k=1,nkn)
	!write(iunit) ((tn(l,k,1),l=1,ilk(k)),k=1,nkn)

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	end

!****************************************************************

	subroutine off_read_record(iu,ig,it,ierr)

! reads one record from offline file

	use mod_offline

	implicit none

	integer iu,ig,it
	integer ierr

	integer ie,ii,k,i
	integer nlin,nlink,nline
	integer iunit
	integer nknaux,nelaux
	integer nkn,nel,nlv,nlvdi
	integer itype
	
	double precision, save, allocatable :: wnaux(:,:)
	double precision, save, allocatable :: rlin(:)
	integer, save, allocatable :: ileaux(:)
	integer, save, allocatable :: ilkaux(:)

	if( nkn_off <= 0  ) then
	  stop 'error stop off_read_record: offline not initialized'
	end if

!----------------------------------------------------------
! initialize
!----------------------------------------------------------

	nkn = nkn_off
	nel = nel_off
	nlv = nlv_off
	nlvdi = nlv

	iunit = iu

!----------------------------------------------------------
! read header
!----------------------------------------------------------

	read(iunit,err=99,end=98) it,nknaux,nelaux,itype
	if( nkn .ne. nknaux .or. nel .ne. nelaux ) goto 97
	if( itype .ne. 3 ) goto 96
	iread = itype

	time(ig) = it

!----------------------------------------------------------
! read vertical indices
!----------------------------------------------------------

	if( .not. allocated(ileaux) ) allocate(ileaux(nel))
	if( .not. allocated(ilkaux) ) allocate(ilkaux(nkn))
	read(iunit,iostat=ierr) (ileaux(ie),ie=1,nel)
	call off_error(ierr,it,'reading ile')
	read(iunit,iostat=ierr) (ilkaux(k),k=1,nkn)
	call off_error(ierr,it,'reading ilk')

	if( .not. bvinit ) then
	  call off_init_vertical(nkn,nel,ileaux,ilkaux)
	end if

	call off_check_vertical(nel,ileaux,ile)
	call off_check_vertical(nkn,ilkaux,ilk)

!----------------------------------------------------------
! set up auxiliary arrays
!----------------------------------------------------------

        call count_linear(nlvdi,nkn,1,ilk,nlink)
        call count_linear(nlvdi,nel,1,ile,nline)
	nlin = max(nlink,nline)
        if( .not. allocated(rlin) ) allocate(rlin(nlin))
        if( .not. allocated(wnaux) ) allocate(wnaux(nlvdi,nkn))

!----------------------------------------------------------
! read currents
!----------------------------------------------------------

	nlin = nline
        read(iunit,iostat=ierr) (rlin(i),i=1,nlin)
	call off_error(ierr,it,'reading ut')
        call dlinear2vals(nlvdi,nel,1,ile,ut(1,1,ig),rlin,nlin)
        read(iunit,iostat=ierr) (rlin(i),i=1,nlin)
	call off_error(ierr,it,'reading vt')
        call dlinear2vals(nlvdi,nel,1,ile,vt(1,1,ig),rlin,nlin)

	!read(iunit) ((ut(l,ie,ig),l=1,ile(ie)),ie=1,nel)
	!read(iunit) ((vt(l,ie,ig),l=1,ile(ie)),ie=1,nel)

!----------------------------------------------------------
! read water levels and vertical velocities
!----------------------------------------------------------

	nlin = nlink
	read(iunit,iostat=ierr) ((ze(ii,ie,ig),ii=1,3),ie=1,nel)
	call off_error(ierr,it,'reading ze')
        read(iunit,iostat=ierr) (rlin(i),i=1,nlin)
	call off_error(ierr,it,'reading wn')
        call dlinear2vals(nlvdi,nkn,1,ilk,wnaux,rlin,nlin)
	wn(1:nlvdi,:,ig) = wnaux(1:nlvdi,:)
	wn(0,:,ig) = 0.
	!read(iunit) ((wn(l,k,ig),l=1,ilk(k)),k=1,nkn)
	read(iunit,iostat=ierr) (zn(k,ig),k=1,nkn)
	call off_error(ierr,it,'reading zn')

!----------------------------------------------------------
! read T/S
!----------------------------------------------------------

	nlin = nlink
        read(iunit,iostat=ierr) (rlin(i),i=1,nlin)
	call off_error(ierr,it,'reading sn')
        call dlinear2vals(nlvdi,nkn,1,ilk,sn(1,1,ig),rlin,nlin)
	nlin = nlink
        read(iunit,iostat=ierr) (rlin(i),i=1,nlin)
	call off_error(ierr,it,'reading tn')
        call dlinear2vals(nlvdi,nkn,1,ilk,tn(1,1,ig),rlin,nlin)

	!read(iunit) ((sn(l,k,ig),l=1,ilk(k)),k=1,nkn)
	!read(iunit) ((tn(l,k,ig),l=1,ilk(k)),k=1,nkn)

!----------------------------------------------------------
! if needed set default values for T/S
!----------------------------------------------------------

	if( idef /= 0 ) then
	  sn = sdef
	  tn = tdef
	end if

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	ierr = 0

	return
   96	continue
	write(6,*) 'type: ',itype
	stop 'error stop off_read: we must have type == 3'
   97	continue
	write(6,*) 'nkn,nknaux: ',nkn,nknaux
	write(6,*) 'nel,nelaux: ',nel,nelaux
	stop 'error stop off_read: parameter mismatch'
   98	continue
	!write(6,*) 'EOF encountered: ',iu,ig
	ierr = -1
	return
	!stop 'error stop off_read: EOF encountered'
   99	continue
	write(6,*) iu,ig
	stop 'error stop off_read: error reading record'
	end

!****************************************************************

	subroutine off_read_header(iu,it,nkn,nel,nlv,ierr)

! reads header and sets nlv

	use mod_offline

	implicit none

	integer iu
	integer it,nkn,nel,nlv
	integer ierr

	integer itype
	integer ie,k,n
	integer nlve,nlvk

	integer, allocatable :: il(:)

	read(iu,err=99,end=98) it,nkn,nel,itype
	if( itype .ne. 3 ) goto 96

! determine nlv - this should be deleted once we write nlv to file

	n = max(nkn,nel)
	allocate(il(n))
	read(iu,iostat=ierr) (il(ie),ie=1,nel)
	call off_error(ierr,it,'off_read_header: reading file')
	nlve = maxval(il(1:nel))
	read(iu,iostat=ierr) (il(k),k=1,nkn)
	call off_error(ierr,it,'off_read_header: reading file')
	nlvk = maxval(il(1:nkn))

	nlv = max(nlve,nlvk)

	ierr = 0

	return
   96	continue
	write(6,*) 'type: ',itype
	stop 'error stop off_read_header: we must have type == 3'
   98	continue
	write(6,*) 'EOF encountered: ',iu
	ierr = -1
	return
   99	continue
	write(6,*) iu
	stop 'error stop off_read_header: error reading record'
	end

!****************************************************************

	subroutine off_peek_header(iu,it,nkn,nel,nlv,ierr)

! peeks into header and sets nlv

	use mod_offline

	implicit none

	integer iu
	integer it,nkn,nel,nlv
	integer ierr

	call off_read_header(iu,it,nkn,nel,nlv,ierr)

	rewind(iu)

	end

!****************************************************************

	subroutine off_init_vertical(nkn,nel,ilhv,ilhkv)

! initializes vertical indices

	use mod_offline

	implicit none

	integer nkn,nel
	integer ilhv(nel)
	integer ilhkv(nkn)

	if( nkn /= nkn_off ) goto 99
	if( nel /= nel_off ) goto 99

	ile = ilhv
	ilk = ilhkv

	bvinit = .true.

	return
   99	continue
	write(6,*) 'nkn,nel: ',nkn,nel
	write(6,*) 'nkn_off,nel_off: ',nkn_off,nel_off
	stop 'error stop off_init_vertical: non compatible'
	end

!****************************************************************

	subroutine off_check_vertical(n,ilaux,il)

! checks if vertical indices have been initialized

	implicit none

	integer n
	integer ilaux(n)
	integer il(n)

	integer i

	do i=1,n
	  if( il(i) .le. 0 ) il(i) = ilaux(i)
	  if( il(i) .ne. ilaux(i) ) then
	    write(6,*) i,il(i),ilaux(i)
	    stop 'error stop off_check_vertical: not compatible'
	  end if
	end do

	end 

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine off_error(ierr,it,text)

! exits with error

	implicit none

	integer ierr,it
	character*(*) text

	if( ierr == 0 ) return 

	write(6,*) '*** error in module offline'
	write(6,*) 'ierr = ',ierr
	write(6,*) 'it = ',it
	write(6,*) trim(text)
	stop 'error stop off_error'

	end

!****************************************************************

	subroutine off_next_record(iu,it,ierr)

! checks info on next record

	implicit none

	integer iu,it,ierr

	integer nknaux,nelaux

	read(iu,err=99,end=98) it,nknaux,nelaux
	backspace(iu)
	ierr = 0

	return
   98	continue
	it = 0
	ierr = -1
	return
   99	continue
	write(6,*) iu
	stop 'error stop off_next_record: error reading record'
	end

!****************************************************************

	subroutine off_set_default_ts(t,s)

	use mod_offline

	real t,s

	idef = 1
	tdef = t
	sdef = s

	end

!****************************************************************
!****************************************************************
!****************************************************************
!
!	program off_main
!	call off_test
!	end
!
!****************************************************************

