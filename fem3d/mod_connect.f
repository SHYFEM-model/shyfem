
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
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

!
! connectivity routines
!
! revision log :
!
! 19.05.2020	ccf	started from scratch
! 12.04.2022	ggu	adapted
!
!****************************************************************

!==================================================================
	module mod_connect
!==================================================================

	implicit none

	integer, save, private :: nkn = 0
	integer, save, private :: nel = 0
	integer, save, private :: ngr = 0

	integer, save, private, allocatable :: nen3v(:,:)

	integer, save, allocatable :: nlist(:,:)
	integer, save, allocatable :: elist(:,:)
	integer, save, allocatable :: ecv(:,:)
	integer, save, allocatable :: bound(:)

!==================================================================
	contains
!==================================================================

	subroutine connect_internal_allocate(nkn_l,nel_l,ngr_l,nen3v_l)

	integer nkn_l,nel_l,ngr_l
	integer nen3v_l(3,nel_l)

	logical brealloc

	brealloc = .true.
	if( ngr > 0 ) then
	  if( nkn == nkn_l .and. nel == nel_l ) brealloc = .false.
	end if

	if( .not. brealloc ) return

call connect_internal_deallocate

	nkn = nkn_l
	nel = nel_l
	ngr = ngr_l

	allocate(nen3v(3,nel))
	allocate(nlist(0:ngr,nkn))
	allocate(elist(0:ngr,nkn))
	allocate(bound(nkn))
	allocate(ecv(3,nel))

	nen3v = nen3v_l

	end

!******************************************************************

	subroutine connect_internal_deallocate

	if( ngr == 0 ) return

	ngr = 0
	nkn = 0
	nel = 0

	deallocate(nen3v)
	deallocate(nlist)
	deallocate(elist)
	deallocate(bound)
	deallocate(ecv)

	end

!******************************************************************

	subroutine connect_internal_init(nkn_l,nel_l,nen3v_l)

	integer nkn_l,nel_l
	integer ngr_l
	integer nen3v_l(3,nel)

	if( ngr > 0 ) call connect_internal_deallocate

	call make_grade(nkn_l,nel_l,nen3v_l,ngr_l)
	call connect_internal_allocate(nkn_l,nel_l,ngr_l,nen3v_l)

	call make_ne_list
	call sort_ne_lists
	call make_bound
	call make_ecv

	end

!******************************************************************

	subroutine connect_internal_parameters(nkn_l,nel_l,ngr_l)

	integer nkn_l,nel_l,ngr_l

	nkn_l = nkn
	nel_l = nel
	ngr_l = ngr

	end

!******************************************************************

	subroutine connect_internal_check

	if( ngr == 0 ) then
	  stop 'error stop connect_internal_check: no initialization'
	end if

	end

!******************************************************************

	subroutine connect_get_element_index(nen3v_l)

	integer nen3v_l(3,nel)

	nen3v_l = nen3v

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine make_grade(nkn,nel,nen3v,ngr)

! makes grade ngr - also computes ngrade and egrade

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	integer ngr		!maximum grade (return)

	integer ie,ii,k,ii1,k1,ngre
	integer idiff,nbnd
	integer, allocatable :: ngrade(:)
	integer, allocatable :: egrade(:)
	integer, allocatable :: nalist(:,:)

	allocate(ngrade(nkn))
	allocate(egrade(nkn))

	egrade = 0

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    egrade(k) = egrade(k) + 1
	  end do
	end do

	ngr = 1 + maxval(egrade)
	allocate(nalist(2*ngr,nkn))
	!write(6,*) 'ngr estimate: ',ngr
	
	ngrade = 0
	nalist = 0

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ii1 = 1 + mod(ii,3)
	    k1 = nen3v(ii1,ie)
	    call insert_in_list(k,k1)
	    call insert_in_list(k1,k)
	  end do
	end do

	ngr = maxval(ngrade)
	ngre = maxval(egrade)

	nbnd = 0
	do k=1,nkn
	  idiff = abs( ngrade(k) - egrade(k) )
	  if( idiff == 1 ) then
	    nbnd = nbnd + 1
	  else if( idiff /= 0 ) then
	    write(6,*) k,idiff,ngrade(k),egrade(k)
	    stop 'error stop: error in grade'
	  end if
	end do

	!write(6,*) 'nkn,nint,nbnd: ',nkn,nkn-nbnd,nbnd

	contains

	subroutine insert_in_list(k1,k2)
	integer k1,k2
	integer n,i
	n = ngrade(k1)
	do i=1,n
	  if( nalist(i,k1) == k2 ) return
	end do
	n = n + 1
	if( n > 2*ngr ) stop 'error stop: ngr'
	nalist(n,k1) = k2
	ngrade(k1) = n
	end subroutine insert_in_list

	end

!******************************************************************

	subroutine make_ne_list!(nkn,nel,nen3v,ngr)

	implicit none

	!integer nkn,nel
	!integer nen3v(3,nel)
	!integer ngr

	integer ie,ii,ii1,k,k1

	nlist = 0
	elist = 0

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ii1 = 1 + mod(ii,3)
	    k1 = nen3v(ii1,ie)
	    call insert_item_in_list(k,k1,nlist)
	    call insert_item_in_list(k1,k,nlist)
	    call insert_item_in_list(k,ie,elist)
	  end do
	end do

	end

!*******************************************************************

	subroutine sort_ne_lists

	implicit none

	logical bbound
	integer ie,k,kn,kb,iee,ii
	integer i,j,nn,ne
	integer in,ib
	integer el(ngr)
	integer nl(ngr)
	character*80 text

	do k=1,nkn
	  ne = elist(0,k)
	  nn = nlist(0,k)
	  el(1:ne) = elist(1:ne,k)
	  nl(1:nn) = nlist(1:nn,k)
	  bbound = ( nn /= ne )
	  if( bbound ) then	!boundary node, first find boundary node
	    text = 'looking for first node'
	    ie = el(1)
	    do
	      iee = ebhnd(ie,k)
	      if( iee == 0 ) exit
	      ie = iee
	    end do
	    kn = knext(ie,k)
	    call exchange_item(k,ie,1,el)
	    call exchange_item(k,kn,1,nl)
	    !call info('after finding first node',k)
	  end if
	  text = 'looking for next nodes'
	  do i=1,ne-1
	    ie = el(i)
	    iee = enext(ie,k)
	    if( iee == 0 ) goto 99
	    kn = knext(iee,k)
	    call exchange_item(k,iee,i+1,el)
	    call exchange_item(k,kn,i+1,nl)
	  end do
	  elist(1:ne,k) = el(1:ne)
	  nlist(1:nn,k) = nl(1:nn)
	  !call info('after sorting all nodes',k)
	end do

	text = 'checking lists'

	do k=1,nkn
	  ne = elist(0,k)
	  nn = nlist(0,k)
	  el(1:ne) = elist(1:ne,k)
	  nl(1:nn) = nlist(1:nn,k)
	  if( any(el(1:ne)==0) ) goto 97
	  if( any(nl(1:nn)==0) ) goto 97
	  bbound = ( nn /= ne )
	  do i=1,ne-1
	    ie = el(i)
	    iee = enext(ie,k)
	    if( iee == 0 ) goto 99
	    kb = kbhnd(ie,k)
	    kn = knext(iee,k)
	    if( kn /= kb ) goto 98
	  end do
	  if( bbound ) cycle
	  kn = knext(el(1),k)
	  kb = kbhnd(el(ne),k)
	  if( kn /= kb ) goto 98
	end do

	return
   97	continue
	write(6,*) 'phase: ',trim(text)
	write(6,*) 'zeros in list: ',k
	call info('checking nodes',k)
	stop 'error stop sort_ne_lists: zeros'
   98	continue
	write(6,*) 'phase: ',trim(text)
	write(6,*) 'k,i: ',k,i
	write(6,*) 'kn,kb not consistent: ',ie,iee,kn,kb
	call info('checking nodes',k)
	stop 'error stop sort_ne_lists: inconsistency'
   99	continue
	write(6,*) 'phase: ',trim(text)
	write(6,*) 'error in list'
	write(6,*)  'ie,iee: ',ie,iee
	call info('on error',k)
	stop 'error stop: generic'

	end

!*******************************************************************

	subroutine make_bound

! makes bound - indicator of boundary nodes

	implicit none

	integer k,n

	bound = 0

	do k=1,nkn
	  n = nlist(0,k) - elist(0,k)
	  if( n == 1 ) then
	    bound(k) = 1
	  else if( n /= 0 ) then
	    write(6,*) k,n
	    stop 'error stop make_bound: error in connectivity'
	  end if
	end do

	end

!*******************************************************************

	subroutine make_ecv

! makes neighbor information of elements

	implicit none

	logical bbound
	integer k,ne,nn,i
	integer ie,ii,iee
	integer el(ngr)

	ecv = 0

	do k=1,nkn
	  ne = elist(0,k)
	  nn = nlist(0,k)
	  el(1:ne) = elist(1:ne,k)
	  do i=1,ne-1
	    call ecv_insert(el(i),el(i+1))
	  end do
	  if( nn == ne ) call ecv_insert(el(ne),el(1))
	end do

! check neighbor list

	do ie=1,nel
	  do ii=1,3
	    iee = ecv(ii,ie)
	    bbound = is_edge_on_boundary(ii,ie)
	    if( iee == 0 ) then
	      if( bbound ) cycle	!edge on boundary... ok
	    else if( any(ecv(:,iee)==ie) ) then
	      cycle				!ie in neighbor list... ok
	    end if
	    write(6,*) ie,ii,iee,bbound
	    write(6,*) ie,ecv(:,ie)
	    if( iee > 0 ) write(6,*) iee,ecv(:,iee)
	    stop 'error stop: checking ecv'
	  end do
	end do

	contains

	subroutine get_edge_nodes(ii,ie,k1,k2)
	integer ii,ie,k1,k2
	integer ii1,ii2
	ii1 = 1+mod(ii,3)
	ii2 = 1+mod(ii1,3)
	k1 = nen3v(ii1,ie)
	k2 = nen3v(ii2,ie)
	end

	function is_edge_on_boundary(ii,ie)
	integer ii,ie
	logical is_edge_on_boundary
	integer ii1,ii2,k1,k2
	call get_edge_nodes(ii,ie,k1,k2)
	is_edge_on_boundary = .false.
	if( elist(0,k1) == nlist(0,k1) ) return
	if( elist(0,k2) == nlist(0,k2) ) return
	is_edge_on_boundary = .true.
	end

	subroutine ecv_insert(ie1,ie2)
	integer ie1,ie2
	integer ib,in
	ib = ibhnd(ie2,k)
	in = inext(ie1,k)
	ecv(in,ie1) = ie2
	ecv(ib,ie2) = ie1
	end

	end

!*******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
! local private subroutines
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine info(text,k)
	character*(*) text
	integer k
	integer i,ie,nn,ne
	nn = nlist(0,k)
	ne = elist(0,k)
	write(6,*) '----------------------'
	write(6,*) trim(text)
	write(6,*)  'k',k
	write(6,*)  'nn,ne',nn,ne
	write(6,*)  'norig',nlist(1:nn,k)
	write(6,*)  'eorig',elist(1:ne,k)
	!write(6,*)  'nl',nl(1:nn)
	!write(6,*)  'el',el(1:ne)
	write(6,*)  'element list:'
	do i=1,ne
	  ie = elist(i,k)
	  write(6,*) ie,nen3v(:,ie)
	end do
	write(6,*) '----------------------'
	end

	function ebhnd(ie,k)
	integer ebhnd
	integer ie,k
	integer i,iee,kn,ne
	kn = knext(ie,k)
	ne = elist(0,k)
	do i=1,ne
	  iee = elist(i,k)
	  if( kbhnd(iee,k) == kn ) exit
	end do
	if( i > ne ) iee = 0
	ebhnd = iee
	end

	function enext(ie,k)
	integer enext
	integer ie,k
	integer i,iee,kn,ne
	kn = kbhnd(ie,k)
	ne = elist(0,k)
	do i=1,ne
	  iee = elist(i,k)
	  if( knext(iee,k) == kn ) exit
	end do
	if( i > ne ) iee = 0
	enext = iee
	end

	function kbhnd(ie,k)
	integer kbhnd
	integer ie,k
	integer ii
	if( ie == 0 ) stop 'error stop kbhnd: iee==0'
	do ii=1,3
	  if( nen3v(ii,ie) == k ) exit
	end do
	if( ii > 3 ) stop 'error stop kbhnd: not found'
	ii = 1 + mod(ii+1,3)
	kbhnd = nen3v(ii,ie)
	end

	function knext(ie,k)
	integer knext
	integer ie,k
	integer ii
	if( ie == 0 ) stop 'error stop knext: iee==0'
	do ii=1,3
	  if( nen3v(ii,ie) == k ) exit
	end do
	if( ii > 3 ) stop 'error stop knext: not found'
	ii = 1 + mod(ii,3)
	knext = nen3v(ii,ie)
	end
	
	function ibhnd(ie,k)
	integer ibhnd
	integer ie,k
	integer ii
	if( ie == 0 ) stop 'error stop ibhnd: iee==0'
	do ii=1,3
	  if( nen3v(ii,ie) == k ) exit
	end do
	if( ii > 3 ) stop 'error stop ibhnd: not found'
	ii = 1 + mod(ii+1,3)
	ibhnd = ii
	end

	function inext(ie,k)
	integer inext
	integer ie,k
	integer ii
	if( ie == 0 ) stop 'error stop inext: iee==0'
	do ii=1,3
	  if( nen3v(ii,ie) == k ) exit
	end do
	if( ii > 3 ) stop 'error stop inext: not found'
	ii = 1 + mod(ii,3)
	inext = ii
	end

	subroutine exchange_item(k,item,jp,list)
	integer k,item,jp
	integer list(ngr)
	integer j,iaux
	j = locate_in_array(item,list)
	if( j == 0 ) then
	  write(6,*) 'looking for: ',item
	  call info('looking for item in exchange_item',k)
	  stop 'error stop exchange_item: cannot find item'
	end if
	if( j == jp ) return
	iaux = list(j)
	list(j) = list(jp)
	list(jp) = iaux
	end

	function locate_in_array(value,array)
	integer value,array(:)
	integer locate_in_array
	integer n,i
	n = size(array)
	do i=1,n
	  if( array(i) == value ) exit
	end do
	locate_in_array = i
	if( i > n ) locate_in_array = 0
	end

	subroutine insert_item_in_list(k,item,list)
	integer k,item
	integer list(0:ngr,nkn)
	integer n,i
	n = list(0,k)
	if(k==0) write(6,1000) 1,item,n,list(1:n,k)
	do i=1,n
	  if( list(i,k) == item ) return
	end do
	!if(k==6) write(6,1000) 2,item,n,list(1:n,k)
	n = n + 1
	if( n > ngr ) then
	  write(6,*) n,ngr,k,item
	  write(6,*) list(0:ngr,k)
	  stop 'error stop insert_item_in_list: n>ngr'
	end if
	list(n,k) = item
	list(0,k) = n
	if(k==0) write(6,1000) 3,item,n,list(1:n,k)
 1000	format(15i5)
	end subroutine insert_item_in_list

!==================================================================
	end module mod_connect
!==================================================================

!*******************************************************************
!*******************************************************************
!*******************************************************************
! utility routines - to be called from outside
!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine connect_init(nkn_l,nel_l,nen3v_l)

	use mod_connect

	implicit none

	integer nkn_l,nel_l
	integer nen3v_l(3,nel_l)

	call connect_internal_init(nkn_l,nel_l,nen3v_l)

	end

!*******************************************************************

	subroutine connect_release

	use mod_connect

	implicit none

	call connect_internal_deallocate

	end

!*******************************************************************

	subroutine connect_get_grades(nkn_l,ngrade,egrade,ngr_l)

	use mod_connect

	implicit none

	integer nkn_l
	integer ngrade(nkn_l)
	integer egrade(nkn_l)
	integer ngr_l

	integer nkn,nel,ngr

	call connect_internal_check

	call connect_internal_parameters(nkn,nel,ngr)

	ngrade(:) = nlist(0,:)
	egrade(:) = elist(0,:)
	ngr_l = ngr

	end

!*******************************************************************

	subroutine connect_get_bound(nkn_l,bound_l)

	use mod_connect

	implicit none

	integer nkn_l
	integer bound_l(nkn_l)

	call connect_internal_check

	bound_l = bound

	end

!*******************************************************************

	subroutine connect_get_ecv(nel_l,ecv_l)

	use mod_connect

	implicit none

	integer nel_l
	integer ecv_l(3,nel_l)

	call connect_internal_check

	ecv_l = ecv

	end

!*******************************************************************

