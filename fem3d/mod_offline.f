
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

c routines for offline data handling
c
c revision log :
c
c 13.06.2013	ggu	new routines written from scratch
c 17.06.2013	ggu	eliminated compiler warnings
c 25.03.2014	ggu	new offline (for T/S)
c 19.12.2014	ggu	changed VERS_7_0_10
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 01.04.2015	ggu	changed VERS_7_1_7
c 06.05.2015	ccf	write offline to .off file
c 06.05.2015	ccf	read offline from offlin file in section name
c 21.05.2015	ggu	changed VERS_7_1_11
c 10.07.2015	ggu	changed VERS_7_1_50
c 13.07.2015	ggu	changed VERS_7_1_51
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 05.11.2015	ggu	revisited and checked
c 09.11.2015	ggu	changed VERS_7_3_13
c 16.11.2015	ggu	changed VERS_7_3_14
c 29.03.2017	ggu	bug fix - input file opened on unit 1
c 05.12.2017	ggu	changed VERS_7_5_39
c 12.11.2018	ggu	linear arrays introduced
c 18.12.2018	ggu	changed VERS_7_5_52
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 04.07.2019	ggu	solved problem for vertical velocity (1 index off)
c 17.02.2020	ggu	femtime eliminated
c
c****************************************************************

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

c****************************************************************
	
	subroutine off_write_record(iu,it)

	use mod_offline

	implicit none

	integer iu,it

	integer ie,ii,k,i
	integer nlin,nlink,nline
	integer idout
	integer nkn,nel,nlv,nlvdi
	
	double precision, save, allocatable :: wnaux(:,:)
	double precision, save, allocatable :: rlin(:)

	if( .not. bvinit ) then
	  stop 'error stop off_write_record: bvinit is false'
	else if( nkn_off <= 0  ) then
	  stop 'error stop off_write_record: offline not initialized'
	end if

	idout = iu
	nkn = nkn_off
	nel = nel_off
	nlv = nlv_off
	nlvdi = nlv

        call count_linear(nlvdi,nkn,1,ilk,nlink)
        call count_linear(nlvdi,nel,1,ile,nline)
	nlin = max(nlink,nline)
        if( .not. allocated(rlin) ) allocate(rlin(nlin))
        if( .not. allocated(wnaux) ) allocate(wnaux(nlvdi,nkn))

        wnaux(1:nlvdi,:) = wn(1:nlvdi,:,1)

	write(iu) it,nkn,nel,3
	write(iu) (ile(ie),ie=1,nel)
	write(iu) (ilk(k),k=1,nkn)

        nlin = nline
        call dvals2linear(nlvdi,nel,1,ile,ut,rlin,nlin)
        write(idout) (rlin(i),i=1,nlin)
        nlin = nline
        call dvals2linear(nlvdi,nel,1,ile,vt,rlin,nlin)
        write(idout) (rlin(i),i=1,nlin)
	!write(iu) ((ut(l,ie,1),l=1,ile(ie)),ie=1,nel)
	!write(iu) ((vt(l,ie,1),l=1,ile(ie)),ie=1,nel)

	write(iu) ((ze(ii,ie,1),ii=1,3),ie=1,nel)
        nlin = nlink
        call dvals2linear(nlvdi,nkn,1,ilk,wnaux,rlin,nlin)
        write(idout) (rlin(i),i=1,nlin)
	!write(iu) ((wn(l,k,1),l=1,ilk(k)),k=1,nkn)
	write(iu) (zn(k,1),k=1,nkn)

        nlin = nlink
        call dvals2linear(nlvdi,nkn,1,ilk,sn,rlin,nlin)
        write(idout) (rlin(i),i=1,nlin)
        nlin = nlink
        call dvals2linear(nlvdi,nkn,1,ilk,tn,rlin,nlin)
        write(idout) (rlin(i),i=1,nlin)

	!write(iu) ((sn(l,k,1),l=1,ilk(k)),k=1,nkn)
	!write(iu) ((tn(l,k,1),l=1,ilk(k)),k=1,nkn)

	end

c****************************************************************

	subroutine off_read_record(iu,ig,it,ierr)

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

	nkn = nkn_off
	nel = nel_off
	nlv = nlv_off
	nlvdi = nlv

	read(iu,err=99,end=98) it,nknaux,nelaux,itype
	if( nkn .ne. nknaux .or. nel .ne. nelaux ) goto 97
	if( itype .ne. 3 ) goto 96
	iread = itype

	time(ig) = it

	if( .not. allocated(ileaux) ) allocate(ileaux(nel))
	if( .not. allocated(ilkaux) ) allocate(ilkaux(nel))
	read(iu) (ileaux(ie),ie=1,nel)
	read(iu) (ilkaux(k),k=1,nkn)

	if( .not. bvinit ) then
	  ile = ileaux
	  ilk = ilkaux
	end if

	call off_check_vertical(nel,ileaux,ile)
	call off_check_vertical(nkn,ilkaux,ilk)

	iunit = iu
        call count_linear(nlvdi,nkn,1,ilk,nlink)
        call count_linear(nlvdi,nel,1,ile,nline)
	nlin = max(nlink,nline)
        if( .not. allocated(rlin) ) allocate(rlin(max(nlink,nlin)))
        if( .not. allocated(wnaux) ) allocate(wnaux(nlvdi,nkn))

	nlin = nline
        read(iunit,iostat=ierr) (rlin(i),i=1,nlin)
        call dlinear2vals(nlvdi,nel,1,ile,ut(1,1,ig),rlin,nlin)
	nlin = nline
        read(iunit,iostat=ierr) (rlin(i),i=1,nlin)
        call dlinear2vals(nlvdi,nel,1,ile,vt(1,1,ig),rlin,nlin)
	!read(iu) ((ut(l,ie,ig),l=1,ile(ie)),ie=1,nel)
	!read(iu) ((vt(l,ie,ig),l=1,ile(ie)),ie=1,nel)

	read(iu) ((ze(ii,ie,ig),ii=1,3),ie=1,nel)
	nlin = nlink
        read(iunit,iostat=ierr) (rlin(i),i=1,nlin)
        call dlinear2vals(nlvdi,nkn,1,ilk,wnaux,rlin,nlin)
	wn(1:nlvdi,:,ig) = wnaux(1:nlvdi,:)
	wn(0,:,ig) = 0.
	!read(iu) ((wn(l,k,ig),l=1,ilk(k)),k=1,nkn)
	read(iu) (zn(k,ig),k=1,nkn)

	nlin = nlink
        read(iunit,iostat=ierr) (rlin(i),i=1,nlin)
        call dlinear2vals(nlvdi,nkn,1,ilk,sn(1,1,ig),rlin,nlin)
	nlin = nlink
        read(iunit,iostat=ierr) (rlin(i),i=1,nlin)
        call dlinear2vals(nlvdi,nkn,1,ilk,tn(1,1,ig),rlin,nlin)
	!read(iu) ((sn(l,k,ig),l=1,ilk(k)),k=1,nkn)
	!read(iu) ((tn(l,k,ig),l=1,ilk(k)),k=1,nkn)

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

c****************************************************************

	subroutine off_read_header(iu,nkn,nel,nlv,ierr)

	use mod_offline

	implicit none

	integer iu
	integer nkn,nel,nlv
	integer ierr

	integer itype,it
	integer ie,k,n
	integer nlve,nlvk

	integer, allocatable :: il(:)

	read(iu,err=99,end=98) it,nkn,nel,itype
	if( itype .ne. 3 ) goto 96

! determine nlv - this should be deleted once we write nlv to file

	n = max(nkn,nel)
	allocate(il(n))
	read(iu) (il(ie),ie=1,nel)
	nlve = maxval(il(1:nel))
	read(iu) (il(k),k=1,nkn)
	nlvk = maxval(il(1:nkn))

	nlv = max(nlve,nlvk)

	ierr = 0

	return
   96	continue
	write(6,*) 'type: ',itype
	stop 'error stop off_read: we must have type == 3'
   98	continue
	write(6,*) 'EOF encountered: ',iu
	ierr = -1
	return
   99	continue
	write(6,*) iu
	stop 'error stop off_read: error reading record'
	end

c****************************************************************

	subroutine off_init_vertical(nkn,nel,ilhv,ilhkv)

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

c****************************************************************

	subroutine off_check_vertical(n,ilaux,il)

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

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine off_next_record(iu,it,ierr)

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

c****************************************************************
c****************************************************************
c****************************************************************
c
c	program off_main
c	call off_test
c	end
c
c****************************************************************

