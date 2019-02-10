
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

! hash routines
!
! revision log :
!
! 30.06.2000	ggu	error check for hash table half full
! 02.12.2011	ggu	initialize hashin to zero
! 22.11.2017	ggu	new framework finished
!
!******************************************************************

!===============================================================
        module hash
!===============================================================

        implicit none

        private

        type :: entry
          integer :: nsize
          integer :: nfull
          integer :: navail
          integer :: ipos			!for visiting
          integer, allocatable :: key(:)
          integer, allocatable :: info(:)
        end type entry

        integer, parameter ::  empty_flag = -2147483601
        integer, parameter :: remove_flag = -2147483602

        logical, parameter :: bdebug = .true.
        integer, parameter :: ndim_first = 10
        integer, parameter :: size_default = 100
        logical, parameter :: bremove = .false.	!use slot of removed items

        integer, parameter :: plist(19) = (/
     +				 7,13,23,47,97,181,379
     +				,739,1559,3541,7817,13513
     +				,21973,54563,99991,223469
     +				,413827,745741,1299533
     +				/)

        integer, save :: idlast = 0
        integer, save :: ndim = 0
        type(entry), save, allocatable :: pentry(:)

        public :: hash_init       !call hash_init(id), call hash_init(id,size)
        public :: hash_delete     !call hash_delete(id)

        public :: hash_insert     !logical hash_insert(id,key,info)
        public :: hash_substitute !logical hash_substitute(id,key,info)
        public :: hash_retrieve   !logical hash_retrieve(id,key,info)
        public :: hash_remove     !logical hash_remove(id,key,info)

        public :: hash_reset      !call hash_reset(id)
        public :: hash_visit      !logical hash_visit(id,key,info)

        public :: hash_info       !call hash_info(id)
        public :: hash_check      !call hash_check(id)

        INTERFACE hash_init
        MODULE PROCEDURE         hash_init_def,hash_init_num
        END INTERFACE

!===============================================================
        contains
!===============================================================

        subroutine hash_init_alloc

	implicit none

        type(entry), allocatable :: paux(:)

        if( ndim == 0 ) then
          ndim = ndim_first
          allocate(pentry(ndim))
        else
          ndim = ndim*2
          allocate(paux(ndim))
          paux(1:ndim/2) = pentry(1:ndim/2)
          call move_alloc(paux,pentry)
        end if

        end subroutine hash_init_alloc

!******************************************************************

        subroutine hash_init_new_id(id)

	implicit none
        integer id

        idlast = idlast + 1
        if( idlast > ndim ) then
          call hash_init_alloc
        end if
        id = idlast

        call hash_init_id(id)

        end subroutine hash_init_new_id

!******************************************************************

        subroutine hash_init_id(id)

	implicit none
        integer id

        if( id > ndim ) then
          stop 'error stop hash_init_id: ndim'
        end if

        pentry(id)%nsize = 0
        pentry(id)%nfull = 0
        pentry(id)%navail = 0
        pentry(id)%ipos = 0

        if( allocated(pentry(id)%key) ) deallocate(pentry(id)%key)
        if( allocated(pentry(id)%info) ) deallocate(pentry(id)%info)

        end subroutine hash_init_id

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine hash_info(id)

	implicit none
	integer id
	integer nempty,nremove,nuse,nsize

	if( 1 == 2 ) then
	  write(6,*) 'hash_info: ',id
     +				,pentry(id)%nsize
     +				,pentry(id)%nfull
     +				,pentry(id)%navail
	else
	  write(6,*) 'hash_info: ',id
      	  write(6,*) 'nsize:  ',pentry(id)%nsize
      	  write(6,*) 'nfull:  ',pentry(id)%nfull
      	  write(6,*) 'navail: ',pentry(id)%navail
      	  nempty = count(pentry(id)%key==empty_flag)
      	  nremove = count(pentry(id)%key==remove_flag)
      	  nsize = pentry(id)%nsize
      	  write(6,*) 'empty:  ',nempty
      	  write(6,*) 'remove: ',nremove
      	  write(6,*) 'use:    ',nsize-nempty-nremove
	end if

	end subroutine hash_info

!********************************************

	subroutine hash_info_content(id)

	implicit none
	integer id
	integer i,key,info

	i = 0
	write(6,*) 'info content: ',id
	call hash_reset(id)
	do while( hash_visit(id,key,info) )
	  i = i + 1
	  write(6,*) i,key,info
	end do

	end subroutine hash_info_content

!********************************************

	subroutine hash_check(id)

	implicit none
	integer id
	integer nsize,norig,navail
	integer nfree,nrem,nkey,i

        nsize = pentry(id)%nsize
        norig = (nsize-1)/2
        navail = pentry(id)%navail
	nfree = 0
	nrem = 0
	nkey = 0
	do i=1,nsize
          if( pentry(id)%key(i) == empty_flag ) then
	    nfree = nfree + 1
          else if( pentry(id)%key(i) == remove_flag ) then
	    nrem = nrem + 1
	  else
	    nkey = nkey + 1
	  end if
	end do
	!write(6,*) 'check: ',nfree,nrem,nkey,norig
	!write(6,*) 'check: ',nfree,nrem,nkey,norig,nsize
	if( nrem+nfree+nkey /= nsize ) then
	  write(6,*) nfree,nrem,nkey,norig
	  stop 'error stop hash_check: internal error (1)'
	end if

	end subroutine hash_check

!********************************************

	subroutine hash_delete(id)

	implicit none
	integer id

	call hash_init_id(id)
	if( id == idlast ) idlast = idlast - 1

	end subroutine hash_delete

!********************************************

	subroutine hash_reallocate(id)

	implicit none
	integer id
	integer nsize,idnew,key,info,nsold

      	nsold = pentry(id)%nsize
	nsize = nsold
	!write(6,*) 'realloc start: ',id,nsize
	call hash_check(id)
	!call hash_info(id)
	!call hash_info_content(id)
	call hash_init(idnew,nsize)
      	nsize = pentry(idnew)%nsize
	!write(6,*) 'realloc new: ',idnew,nsize
	if( bdebug ) write(6,*) 'hash realloc: ',nsold,nsize
	call hash_reset(id)
	do
	  if( .not. hash_visit(id,key,info) ) exit
	  if( .not. hash_insert(idnew,key,info) ) then
	    stop 'error stop hash_reallocate: internal error (1)'
	  end if
	end do
	pentry(id) = pentry(idnew)
	!call hash_info(id)
	!call hash_info_content(id)
	call hash_delete(idnew)

	end subroutine hash_reallocate

!********************************************

	function hash_get_prime(nsize)

	implicit none
	integer hash_get_prime
	integer nsize
	integer n,i

	n = size(plist)
	do i=1,n
	  if( nsize < plist(i) ) exit
	end do
	if( i > n ) stop 'error stop hash_get_prime: no prime available'
	hash_get_prime = plist(i)

	end function hash_get_prime

!********************************************

	subroutine hash_init_def(id)

	implicit none
	integer id

	call hash_init_num(id,size_default)

	end subroutine hash_init_def
	
!********************************************

	subroutine hash_init_num(id,size)

	implicit none
	integer id,size
	integer nsize

        call hash_init_new_id(id)
	nsize = hash_get_prime(2*size)
        pentry(id)%nsize = nsize
        pentry(id)%navail = (nsize-1)/2
        allocate(pentry(id)%key(nsize))
        allocate(pentry(id)%info(nsize))
        pentry(id)%key = empty_flag
        pentry(id)%info = 0

	end subroutine hash_init_num

!***********************************************************

	function hash_lookup(id,key,brem)

	integer hash_lookup
	integer id,key
	logical brem			!return key with delete flag

	integer icount,inc,ndim,ndim2,i,lookup,pkey
	integer idel
	logical debug

	debug = ( key == 33533 )
	debug = .true.
	debug = .false.

	icount=0
	inc=1
	ndim = pentry(id)%nsize
	ndim2 = (ndim-1) / 2
	i = mod(key,ndim) + 1
	lookup = 0
	idel = 0

	if( debug ) write(6,*) 'looking for: ',key,i,brem

	do while( icount < ndim2 )
	  pkey = pentry(id)%key(i)
	  if( debug ) write(6,*) 'look:  ',i,inc,icount,ndim,pkey
	  if( pkey == empty_flag ) then
	    lookup=-i
	    if( idel > 0 ) lookup = -idel
	    exit
	  else if( brem .and. pkey == remove_flag ) then
	    if( idel == 0 ) idel = i
	    !lookup=-i
	    !exit
	  else if( pkey == key ) then
	    lookup=i
	    exit
	  end if
	  i=i+inc
	  if( i .gt. ndim ) i = mod(i,ndim) + 1
	  inc=inc+2
	  icount=icount+1
	end do

	if( icount == ndim2 .and. idel > 0 ) lookup = -idel

	if( debug ) write(6,*) 'found: ',lookup,inc,icount,ndim,pkey

	if( lookup == 0 ) then
	  stop 'error stop hash_lookup: internal error (1)'
	end if

	hash_lookup = lookup

	end function hash_lookup

!******************************************************************

	function hash_insert(id,key,info)

	implicit none
	logical hash_insert
	integer id,key,info,i

	i=hash_lookup(id,key,bremove)
	if(i.gt.0) then	!found - cannot insert
	  hash_insert = .false.
	else
	  pentry(id)%key(-i) = key
	  pentry(id)%info(-i) = info
	  pentry(id)%nfull = pentry(id)%nfull + 1
	  pentry(id)%navail = pentry(id)%navail - 1
	  hash_insert = .true.
	end if
	if( pentry(id)%navail < 0 ) call hash_reallocate(id)

	end function hash_insert

!********************************************

	function hash_substitute(id,key,info)

	implicit none
	logical hash_substitute
	integer id,key,info,i

	i=hash_lookup(id,key,bremove)
	if(i.gt.0) then	!found - substitute
	  pentry(id)%info(i) = info
	  hash_substitute = .true.
	else
	  pentry(id)%key(-i) = key
	  pentry(id)%info(-i) = info
	  pentry(id)%nfull = pentry(id)%nfull + 1
	  pentry(id)%navail = pentry(id)%navail - 1
	  hash_substitute = .false.
	end if
	if( pentry(id)%navail < 0 ) call hash_reallocate(id)

	end function hash_substitute

!********************************************

	function hash_retrieve(id,key,info)

	implicit none
	logical hash_retrieve
	integer id,key,info,i

	i=hash_lookup(id,key,.false.)
	!write(6,*) 'after lookup: ',key,i
	if(i.gt.0) then	!found
	  info = pentry(id)%info(i)
	  hash_retrieve = .true.
	else
	  info = 0
	  hash_retrieve = .false.
	end if

	end function hash_retrieve

!********************************************

	function hash_remove(id,key,info)

	implicit none
	logical hash_remove
	integer id,key,info,i

	i=hash_lookup(id,key,.false.)
	if(i.gt.0) then	!found
	  info = pentry(id)%info(i)
	  pentry(id)%key(i) = remove_flag
	  if( bremove ) then
	    pentry(id)%nfull = pentry(id)%nfull - 1
	    pentry(id)%navail = pentry(id)%navail + 1
	  end if
	  hash_remove = .true.
	else
	  info = 0
	  hash_remove = .false.
	end if

	end function hash_remove

!********************************************

	function hash_visit(id,key,info)

	implicit none
	logical hash_visit
	integer id,key,info
	integer ipos,ndim,i

	ipos = pentry(id)%ipos + 1
	ndim = pentry(id)%nsize
	do i=ipos,ndim
	  !write(6,*) 'visit: ',id,i,pentry(id)%key(i)
	  if( pentry(id)%key(i) /= empty_flag .and.
     +			pentry(id)%key(i) /= remove_flag ) exit
	end do
	if(i.le.ndim) then
	  key = pentry(id)%key(i)
	  info = pentry(id)%info(i)
	  pentry(id)%ipos = i
	  hash_visit = .true.
	else
	  key = 0
	  info = 0
	  hash_visit = .false.
	end if

	end function hash_visit

!********************************************

	subroutine hash_reset(id)

	implicit none
	integer id

	pentry(id)%ipos = 0

	end subroutine hash_reset

!===============================================================
        end module hash
!===============================================================


!*************************************************************

	subroutine test_hash_inter

! tests hash routines

	use hash

	implicit none

	integer id,mode
	integer key,info
	logical bok

	call hash_init(id)

	mode=1
	do while(mode.ne.0)
          write(6,*) '0 exit  1 insert  2 retriv  3 remove  '//
     +				'4 visit  5 info'
	  read(5,'(i10)') mode

	  if(mode.eq.1) then
		write(6,*) 'key,info :'
		read(5,*) key,info
		bok=hash_insert(id,key,info)
		write(6,*) bok,key,info
	  else if(mode.eq.2) then
		write(6,*) 'key :'
		read(5,*) key
		bok=hash_retrieve(id,key,info)
		write(6,*) bok,key,info
	  else if(mode.eq.3) then
		write(6,*) 'key :'
		read(5,*) key
		bok=hash_remove(id,key,info)
		write(6,*) bok,key,info
	  else if(mode.eq.4) then
		call hash_reset(id)
		do while( hash_visit(id,key,info) )
		  write(6,*) key,info
		end do
	  else if(mode.eq.5) then
		call hash_info(id)
	  end if
	end do

	end

!********************************************************************

	subroutine test_hash_auto

	use hash

	implicit none

	integer, parameter :: ndim = 50000
	integer, parameter :: maxnum = 50000

	logical bexist,bwrite,bcheck
	integer nmax,id,what,key,info,i,ic
	integer keys(ndim),visit(ndim)

	integer hash_irand
	logical hash_prob

!       ' 1 insert  2 retriv  3 remove  4 visit '

	bwrite = .true.
	bwrite = .false.
	bcheck = .true.		!this is slow
	bcheck = .false.

	nmax = 0
	ic = 0
	call hash_rand_init
	call hash_init(id,2)
	!call hash_info(id)

	do

	  what = hash_irand(4)
	  ic = ic + 1
	  if( ic > 300000 ) exit
	  if( nmax >= ndim ) exit
	  if( what == 1 ) then		!insert
	!if( key == 33533 ) then
	!  write(6,*) 'ggggggggg ',what,key
	!  do i=1,nmax
	!    if( keys(i) == key ) write(6,*) 'www: ',i,key,nmax
	!  end do
	!end if
	    key = hash_irand(maxnum)
	    bexist = any( keys(1:nmax) == key )
	    if( bwrite ) write(6,*) 'what: ',what,key,bexist,nmax,ic
	    if( hash_insert(id,key,key*2) ) then
	      call hash_assert(.not.bexist,'insert yes',id)
	      nmax = nmax + 1
	      call hash_assert(nmax <= ndim,'ndim',id)
	      keys(nmax) = key
	    else
	      call hash_assert(bexist,'insert no',id)
	    end if
	  else if( what == 2 ) then		!retrieve
	    if( hash_prob(0.7) .and. nmax > 0 ) then
	      key = keys(hash_irand(nmax))
	    else
	      key = hash_irand(maxnum)
	    end if
	    bexist = any( keys(1:nmax) == key )
	    if( bwrite ) write(6,*) 'what: ',what,key,bexist,nmax,ic
	    if( hash_retrieve(id,key,info) ) then
	      call hash_assert(bexist,'retrieve yes',id)
	      call hash_assert(info==key*2,'retrieve info',id)
	    else
	      call hash_assert(.not. bexist,'retrieve no',id)
	    end if
	  else if( what == 3 ) then		!remove
	    if( hash_prob(0.7) .and. nmax > 0 ) then
	      key = keys(hash_irand(nmax))
	    else
	      key = hash_irand(maxnum)
	    end if
	    bexist = any( keys(1:nmax) == key )
	    if( bwrite ) write(6,*) 'what: ',what,key,bexist,nmax,ic
	    if( hash_remove(id,key,info) ) then
	      call hash_assert(bexist,'remove yes',id)
	      call hash_assert(info==key*2,'remove info',id)
	      do i=1,nmax
	        if( keys(i) == key ) exit
	      end do
	      call hash_assert(i<=nmax,'remove found',id)
	      keys(i) = keys(nmax)
	      nmax = nmax - 1
	    else
	      call hash_assert(.not. bexist,'remove no',id)
	    end if
	  else if( what == 4 ) then		!substitute
	    if( hash_prob(0.7) .and. nmax > 0 ) then
	      key = keys(hash_irand(nmax))
	    else
	      key = hash_irand(maxnum)
	    end if
	    bexist = any( keys(1:nmax) == key )
	    if( bwrite ) write(6,*) 'what: ',what,key,bexist,nmax,ic
	    if( hash_substitute(id,key,key*2) ) then
	      call hash_assert(bexist,'subst yes',id)
	    else
	      call hash_assert(.not.bexist,'subst no',id)
	      nmax = nmax + 1
	      call hash_assert(nmax <= ndim,'ndim',id)
	      keys(nmax) = key
	    end if
	  end if
	  if( bcheck ) call hash_check(id)
	end do

	write(6,*) 'hash test successfully finished: ',nmax,ic

	end

!********************************************************************

	subroutine hash_assert(bcheck,text,id)
	use hash
	implicit none
	logical bcheck
	character*(*) text
	integer id
	if( .not. bcheck ) then
	  write(6,*) 'hash_assertion: ',trim(text)
	  call hash_info(id)
	  stop 'assertion failed'
	end if
	end

	subroutine hash_rand_init
	call random_seed
	end

	function hash_irand(max)
	implicit none
	integer hash_irand
	integer max
	real r
	call random_number(r)
	hash_irand = 1 + floor(r*max)
	!write(6,*) max,hash_irand,r
	end

	function hash_prob(ptrue)
	implicit none
	logical hash_prob
	real ptrue
	real r
	call random_number(r)
	hash_prob = ( r < ptrue )
	end

!********************************************************************

	program test_hash_main
	call test_hash_auto
	!call test_hash_inter
	end

!********************************************************************
