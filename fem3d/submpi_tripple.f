
!--------------------------------------------------------------------------
!
!    Copyright (C) 2007-2008,2010-2011,2013-2015,2013-2015  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
!    Copyright (C) 2008,2011,2013,2016  Debora Bellafiore
!    Copyright (C) 2008  Christian Ferrarin
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

! tripple point routines
!
! revision log :
!
! 18.03.2023	ggu	adjusted horizontal diffusion for mpi
! 27.03.2023	ggu	tripple point routines copied here
! 28.03.2023	ggu	insert lmax into list
! 29.03.2023	ggu	finished, new routines for getting vals
! 02.04.2023	ggu	bug fix looking for internal element
! 03.04.2023	ggu	new bug fix looking for internal element
! 12.04.2023	ggu	fill buffer_in only if my_id == id_from
! 07.06.2023    ggu     in exchange_elem_info() use utlov not utlnv (bug-fix)
! 09.06.2023    ggu     bug fix for tripple point on external boundary
!
!******************************************************************

! notes :
!
!	iexch(nexch,itmax)
!
!	  iexch(1,i) = ie	!internal element number
!	  iexch(2,i) = ipev(ie)	!external element number
!	  iexch(3,i) = my_id	!id of domain
!	  iexch(8,i) = i0	!neibor element opposite to this vertex
!
!	  iexch(4,i) = iint	!internal number of neighbor element
!	  iexch(5,i) = iext	!external number of neighbor element
!	  iexch(6,i) = ide	!domain of neighbor element
!	  iexch(7,i) = lmax	!max level of neighbor element
! 
!==================================================================
        module shympi_tripple
!==================================================================

	logical, parameter :: btripple = .true. !handles tripple points

	integer, parameter :: nexch = 8
	integer, save :: nmax_tripple = 0
	integer, save :: itrtot = -1
	integer, save, allocatable :: ielist(:,:)
	real, save, allocatable :: buffer_tripple(:,:)
	integer, save, allocatable :: ietrp(:,:)

!==================================================================
        end module shympi_tripple
!==================================================================

	subroutine tripple_points_handle

	use shympi
	use shympi_tripple

	implicit none

	if( itrtot == -1 ) call tripple_points_init

	call shympi_syncronize

	if( .not. btripple ) return	!we do not handle tripple points

	call tripple_points_exchange

	end

!******************************************************************

	subroutine tripple_points_init

! constructs list of tripple points and info on how to exchange

	use basin
	use shympi
	use shympi_tripple

	implicit none

	integer ie,itr,idn,ide,iee
	integer ii,i1,i2,i0,ip
	integer k,k1,k2,kext1,kext2
	integer i,j,iei,iext,iint,ia,itmax
	integer iu
	integer, allocatable :: ies(:)

	integer, allocatable :: itrs(:)
	integer, allocatable :: iexch(:,:)
	integer, allocatable :: iexchs(:,:)
	integer, allocatable :: ips(:,:)
	integer, allocatable :: id_elem_g(:,:)

	integer ieext,ieint,ipext

	if( itrtot >= 0 ) return	!already set up

	call shympi_syncronize
	write(6,*) 'starting tripple_points_init: ',my_id

	iu = 0
	if( bmpi_debug ) iu = 300 + my_id

	!--------------------------------------------------
	! get general information on tripple points
	!--------------------------------------------------

	itr = 0
	do ie=1,nel_unique
	  idn = id_elem(0,ie)
	  if( idn == 3 ) then
	    write(6,1000) 'tripple point found: ',my_id,ie,id_elem(:,ie)
	    itr = itr + 1
	  end if
	end do

	call shympi_syncronize

	itrtot = shympi_sum(itr)
	if( my_id == 0 ) then
	  write(6,*) 'summary for tripple points:'
	  write(6,*) '                                  '//
     +			'     domain         itr      itrtot'
	end if
	call shympi_syncronize
	write(6,*) 'total numbers of tripple points: ',my_id,itr,itrtot

	nmax_tripple = 2 * (nlv_global+1)
	allocate(buffer_tripple(nmax_tripple,itrtot))
	allocate(ielist(nexch,itrtot))
	ielist = 0

	allocate(ietrp(3,nel))
	ietrp = 0

	if( itrtot == 0 ) return

	if( .not. btripple ) then
	  write(6,*) 'not handling tripple points: ',my_id,itrtot
	  return
	end if

	call shympi_syncronize

	!--------------------------------------------------
	! allocate arrays
	!--------------------------------------------------

	itmax = shympi_max(itr)

	allocate(ies(itmax))
	allocate(iexch(nexch,itmax))
	allocate(iexchs(nexch,n_threads))
	allocate(itrs(n_threads))
	allocate(ips(itmax,n_threads))
	iexch = 0
	ies = 0
	itr = 0
	itrs = 0

	!--------------------------------------------------
	! insert local tripple elements in ies
	!--------------------------------------------------

	do ie=1,nel_unique
	  idn = id_elem(0,ie)
	  if( idn == 3 ) then
	    itr = itr + 1
	    ies(itr) = ie
	  end if
	end do

	call shympi_gather(itr,itrs)

	!--------------------------------------------------
	! find neighbor elements to the tripple elements
	!--------------------------------------------------

	do i=1,itr
	  ie = ies(i)
	  iexch(1,i) = ie		!internal element number
	  iexch(2,i) = ipev(ie)		!external element number
	  iexch(3,i) = my_id		!id of domain
	  call find_elem_neib(ie,i0,iext)
	  iexch(5,i) = iext		!external number of neibor eleement
	  iexch(8,i) = i0		!neibor element opposite to this vertex
	end do

	call shympi_syncronize

	!--------------------------------------------------
	! find id of neighbor element and remember in iexch
	!--------------------------------------------------

	call find_elem_id(nexch,itr,iexch)

	!--------------------------------------------------
	! construct pointer where to insert item in list
	!--------------------------------------------------

	ip = 0
	ips = 0
	do ia=1,n_threads
	  itr = itrs(ia)
	  do i=1,itr
	    ip = ip + 1
	    ips(i,ia) = ip	!pointer where to insert info
	  end do
	end do

	if( ip /= itrtot ) stop 'error stop tripple_points: internal (3)'

	!--------------------------------------------------
	! insert items in list (list will be identical in all domains)
	!--------------------------------------------------

	do i=1,itmax
	  call shympi_gather(iexch(:,i),iexchs)
	  do ia=1,n_threads
	    ip = ips(i,ia)
	    if( ip > 0 ) then
	      ielist(:,ip) = iexchs(:,ia)
	    end if
	  end do
	end do

	call shympi_syncronize

	!--------------------------------------------------
	! create element pointer for tripple points
	!--------------------------------------------------

	do i=1,itrtot
	  ide = ielist(3,i)
	  iext = ielist(5,i)		!only handle if iext > 0
	  if( ide /= my_id ) cycle	!only handle this domain
	  ie = ielist(1,i)
	  ii = ielist(8,i)
	  if( ii < 1 .or. ii > 3 ) then
	    stop 'error stop tripple_points: internal (21)'
	  end if
	  if( iext > 0 ) ietrp(ii,ie) = i	!pointer into tripple list
	end do

	!--------------------------------------------------
	! debug output
	!--------------------------------------------------

	call shympi_syncronize

	if( shympi_is_master() ) then
	  write(6,*) 'list of tripple points: ',itrtot
	  call ielist_info
	end if

	if(iu>0) write(iu,*) 'list of tripple points:'
	do i=1,itrtot
	  write(6,'(10i8)') my_id,ielist(:,i)
	  if(iu>0) write(iu,'(10i8)') my_id,ielist(:,i)
	  call shympi_gather(ielist(:,i),iexchs)
	  do j=1,nexch
	    if( any( iexchs(j,:) /= iexchs(j,1) ) ) then
	      write(6,*) 'values are different: ',iexchs(j,:)
	      stop 'error stop tripple_points: internal (11)'
	    end if
	  end do
	  iext = ielist(5,i)
	  iint = ieint(iext)
	  if( iint > 0 .and. iint <= nel_unique ) then
	    if( iint /= ielist(4,i) ) then
	      write(6,'(a,10i8)') 'error: ',my_id,iext,iint,nel_unique
	      stop 'error stop tripple_points: internal (8)'
	    end if
	  end if
	end do

	call shympi_syncronize

	write(6,*) 'finished tripple_points_init: ',my_id

	!stop	!debug stop

	!--------------------------------------------------
	! end of routine
	!--------------------------------------------------

	return
 1000	format(a,6i8)
	end

!******************************************************************

	subroutine tripple_points_exchange

! exchanges information on tripple points

	use basin
	use levels
	use shympi
	use shympi_tripple

	implicit none

	integer i,itr
	integer iint,iext,id
	integer id_from,id_to
	integer nmax,lmax

	integer ieint

	if( itrtot == 0 ) return
	if( itrtot < 0 ) then
	  stop 'error stop tripple_points_exchange: no init'
	end if

	do i=1,itrtot
	  !write(6,'(a,10i8)') 'tr_exchange: ',my_id,ielist(:,i)
	  iint = ielist(4,i)
	  iext = ielist(5,i)
	  id_to = ielist(3,i)
	  id_from = ielist(6,i)
	  lmax = ielist(7,i)
	  itr = i
	  call exchange_elem_info(id_from,id_to,iint,iext,lmax,itr)
	end do

	call shympi_syncronize
	
	end

!******************************************************************

	subroutine exchange_elem_info(id_from,id_to,iint,iext,lmax,itr)

	use levels
	use evgeom
	use mod_hydro
	use shympi
	use shympi_tripple

	implicit none

	integer id_from,id_to,iint,iext,lmax,itr

	logical bdebug
	integer iu
	integer l,n
	integer nmax
	real buffer_in(nmax_tripple)
	real buffer_out(nmax_tripple)

	integer ieext

	bdebug = .false.

	nmax = nmax_tripple
	n = 2*lmax + 2

        buffer_in = 0.
        if( id_from == my_id ) then
	  buffer_in(1) = lmax
	  buffer_in(2) = 12. * ev(10,iint)
	  buffer_in(3:lmax+2) = utlov(1:lmax,iint)
	  buffer_in(lmax+3:2*lmax+2) = vtlov(1:lmax,iint)
	end if

	if( n > nmax ) then
	  write(6,*) 'n,nmax: ',n,nmax
	  stop 'error stop exchange_elem_info: n>nmax'
	end if

	if( bmpi ) then
	  iu = 300 + my_id
	  buffer_out = 0.
	  call shympi_receive(id_from,id_to,n,buffer_in,buffer_out)
	  buffer_tripple(:,itr) = buffer_out(:)
	else
	  iu = 400
	  buffer_out = buffer_in
	end if

	if( bdebug ) then
	  write(iu,*) iext
	  write(iu,*) lmax,n
	  write(iu,*) buffer_out(1:n)
	end if

	end

!******************************************************************

	subroutine tripple_point_get_la(itr,lmax,area)

! gets lmax and area from neighbor element

	use shympi_tripple

	implicit none

	integer itr,lmax
	real area

	lmax = nint(buffer_tripple(1,itr))
	area = buffer_tripple(2,itr)

	end

!******************************************************************

	subroutine tripple_point_get_values(itr,nl,lmax,area,u,v)

! gets lmax, area, and u/v from neighbor element

	use shympi_tripple

	implicit none

	integer itr,nl,lmax
	real area
	real u(nl),v(nl)

	u = 0
	v = 0

	lmax = nint(buffer_tripple(1,itr))
	area = buffer_tripple(2,itr)

	if( lmax > nl ) then
	  write(6,*) lmax,nl
	  stop 'error stop tripple_point_get_values: lmax>nl'
	end if

	u(1:lmax) = buffer_tripple(3:lmax+2,itr)
	v(1:lmax) = buffer_tripple(lmax+3:2*lmax+2,itr)

	!write(555,*) lmax,area
	!write(555,*) u(1:lmax)
	!write(555,*) v(1:lmax)

	end 

!******************************************************************
!******************************************************************
!******************************************************************
! utility routines
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine find_elem_id(nexch,itr,iexch)

! computes id and internal ie from list in iexch and inserts it there

	use levels
	use mod_geom
	use shympi

	implicit none

	integer nexch,itr
	integer iexch(nexch,itr)

	integer i,iext,ide,ie,iee,iint,lmax
	integer ii,ineib,iemax,ia
	integer, allocatable :: id_elem_g(:,:)
	integer, allocatable :: ilhv_g(:)

!------------------------------------------------------------
! allocate global arrays
!------------------------------------------------------------

	allocate(id_elem_g(0:3,nel_global))
	allocate(ilhv_g(nel_global))
	call shympi_l2g_array(4,id_elem,id_elem_g)
	call shympi_l2g_array(ilhv,ilhv_g)

!------------------------------------------------------------
! find domain (ide) of neibor element
!------------------------------------------------------------

	do i=1,itr
	  iext = iexch(5,i)
	  ide = -1
	  iexch(6,i) = ide	!domain of neighbor element
	  if( iext == 0 ) cycle	!tripple point on boundary
	  do ie=1,nel_global
	    if( ip_ext_elem(ie) == iext ) then
	      iee = ie		!global internal element number of neibor
	      ide = id_elem_g(1,iee)
	      lmax = ilhv_g(iee)
	    end if
	  end do
	  if( ide == -1 ) stop 'error stop tripple_points: internal (4)'
	  if( ide < -1 .or. ide >= n_threads ) then
	    write(6,*) '*** impossible ide: ',ide
	    stop 'error stop tripple_points: internal (5)'
	  end if
	  !now look for local internal element number in domain ide
	  iint = 0
	  ia = ide + 1
	  iemax = nel_domains(ia)
	  !iemax = ne_max
	  do ie=1,iemax
	    if( iee == ip_int_elems(ie,ia) ) then !scan local element index
	      iint = ie
	    end if
	  end do
	  if( iint == 0 ) then
	    write(6,*) '*** cannot find internal element number: ',my_id
	    write(6,*) i,iext,ide,iee,iint
	    !write(6,'(a,10i7)') 'list ',iexch(:,i)
	    call iexch_info(nexch,iexch)
	    stop 'error stop tripple_points: internal (6)'
	  end if
	  write(6,*) 'ide found: ',iext,iint,ide
	  iexch(4,i) = iint	!internal element number of neighbor
	  iexch(5,i) = iext	!external element number of neighbor
	  iexch(6,i) = ide	!domain of neighbor element
	  iexch(7,i) = lmax	!ma level of neighbor element
	end do

!------------------------------------------------------------
! insert domain pointer into ieltv
!------------------------------------------------------------

	do i=1,itr
	  ie = iexch(1,i)
	  iext = iexch(5,i)
	  ide = iexch(6,i)	!domain of neibor element
	  ii = iexch(8,i)
	  ineib = ieltv(ii,ie)
	  if( ineib == 0 .and. iext == 0 ) then
	    write(6,*) 'not real tripple point... ignoring ',i
	    cycle
	  else if( ineib /= -1000 ) then
	    write(6,*) i,ineib,iext
	    write(6,*) '*** neigbor not flagged as external: ',ineib
	    call error_stop('tripple_points','internal (14)')
	  end if
	  ieltv(ii,ie) = -1000 - ide
	end do

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	end

!******************************************************************

	subroutine find_elem_neib(ie,i0,iext)

! find neighbor element opposite to node with my_id and returns external number

	use basin
	use shympi
	use mod_geom

	implicit none

	integer ie	!internal element number of tripple point element
	integer i0	!vertex opposite of neibor element (return)
	integer iext	!external number of neibor element (return)

	integer ii,i1,i2,iei,iee,ieint,ineigh
	integer k,k1,k2,kg1,kg2
	integer ids(3),ides(3)

	if( ie == 0 ) then
	  stop 'error stop find_elem_neib: ie == 0'
	end if

	ieint = ie
	iee = ipev(ie)

!------------------------------------------------------------
! find vertex opposite to neibor element and nodes on border
!------------------------------------------------------------

	k1 = 0
	k2 = 0
	i0 = 0
	do ii=1,3
	  k = nen3v(ii,ie)
	  if( id_node(k) == my_id ) then
	    i0 = ii
	    ineigh = ieltv(i0,ie)
	    i1 = mod(ii,3) + 1
	    i2 = mod(i1,3) + 1
	    k1 = nen3v(i1,ie)
	    k2 = nen3v(i2,ie)
	  end if
	  ids(ii) = id_node(k)
	  ides(ii) = id_elem(ii,ie)
	end do

	if( i0 == 0 ) then
	  stop 'error stop find_elem_neib: i0 == 0'
	end if

	if( ineigh == 0 ) then
	  write(6,*) 'tripple point is on boundary...'
	  write(6,*) 'no need to get information...'
	  write(6,*) ieint,iee,i0
	  iext = 0
	  return
	end if

!------------------------------------------------------------
! get external node numbers
!------------------------------------------------------------

	kg1 = ip_int_node(k1)
	kg2 = ip_int_node(k2)

!------------------------------------------------------------
! look for external nodes in global domain
!------------------------------------------------------------

	iei = 0
	do ie=1,nel_global
	  do ii=1,3
	    i1 = mod(ii,3) + 1
	    if( nen3v_global(i1,ie) == kg1 
     +			.and. nen3v_global(ii,ie) == kg2 ) then
	      if( iei /= 0 ) stop 'error stop find_elem_neib: int (1)'
	      iei = ie
	    end if
	  end do
	end do

	if( iei == 0 ) then
	  write(6,*) 'cannot find element to nodes ',my_id
	  write(6,*) my_id,iee,kg1,kg2
	  write(6,*) id_elem(:,ieint)
	  write(6,*) ids
	  write(6,*) ides
	  write(6,*) i0,ieltv(:,ieint)
	  call error_stop('find_elem_neib','internal (2)')
	end if

!------------------------------------------------------------
! assign external element number of neibor element
!------------------------------------------------------------

	iext = ip_ext_elem(iei)

	!write(6,1000) 'element to nodes found',my_id,iee,iext,kg1,kg2
 1000	format(a,10i6)

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine ielist_info

	use shympi_tripple

	implicit none

	integer i
	character*80 header

	header = '             ie      iee       id      ' //
     +			'ie1     iee1      id1    lmax1       i0'

	write(6,'(a)') trim(header)

	do i=1,itrtot
	  write(6,'(a,8i9)') 'info: ',ielist(:,i)
	end do

	end

!******************************************************************

	subroutine iexch_info(nexch,iexch)

	implicit none

	integer nexch
	integer iexch(nexch)

	character*80 header

	header = '             ie      iee       id      ' //
     +			'ie1     iee1      id1    lmax1       i0'
	write(6,'(a)') trim(header)
	write(6,'(a,8i9)') 'info: ',iexch(:)

	!  iexch(1,i) = ie
	!  iexch(2,i) = ipev(ie)
	!  iexch(3,i) = my_id
	!  iexch(5,i) = iext
	!  iexch(4,i) = iint
	!  iexch(5,i) = iext
	!  iexch(6,i) = ide
	!  iexch(7,i) = lmax
	!  iexch(8,i) = i0

	end

!******************************************************************
!******************************************************************
!******************************************************************
! old routines (not used)
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine test_exchange_3d

	use basin
	use levels
	use shympi

	implicit none

	integer ie,iee,l,lmax,idiff,ierr
	real val,diff
	real vals(nlvdi,nel)
	real valsaux(nlvdi,nel)
	real vals2d(nel)
	real vals2daux(nel)
	integer ivals2d(nel)
	integer ivals2daux(nel)

	integer ieext

	!write(6,*) 'running test_exchange_3d'

	ivals2d = 0
	ivals2daux = 0

	do ie=1,nel
	  iee = ieext(ie)
	  ivals2daux(ie) = iee
	  if( ie > nel_unique ) cycle
	  ivals2d(ie) = iee
	end do

	call shympi_exchange_2d_elem(ivals2d)

	write(550+my_id,*) nel,nel_unique,nel_local,nel_local-nel_unique
	write(540+my_id,*) nel,nel_unique,nel_local,nel_local-nel_unique

	ierr = 0
	do ie=1,nel
	  idiff = abs(ivals2d(ie)-ivals2daux(ie))
	  if( idiff > 0 ) then
	    ierr = ierr + 1
	    write(550+my_id,*) ierr,ie,ivals2d(ie),ivals2daux(ie),idiff
	  end if
	  write(540+my_id,*) ie,id_elem(:,ie)
	end do

	if( ierr > 0 ) then
	  write(6,*) 'running test_exchange_3d'
	  write(6,*) 'my_id = ',my_id,'  ierr = ',ierr
	  write(6,*) 'output in files 550+ and 540+'
	  write(6,*) '*** errors in test_exchange_3d ***'
	  stop 'error stop test_exchange_3d: internal (30)'
	end if

!-------------------------------------------------------------

	vals = 0.

	do ie=1,nel_unique
	  iee = ieext(ie)
	  lmax = ilhv(ie)
	  do l=1,lmax
	    val = 100*iee + l
	    vals(l,ie) = val
	  end do
	  vals2d(ie) = iee
	end do

	valsaux = vals
	vals2daux = vals2d
	call shympi_exchange_3d_elem(vals)
	call shympi_exchange_2d_elem(vals2d)

	write(670+my_id,*) nel,nel_unique,nel_local
	write(680+my_id,*) nel,nel_unique,nel_local
	write(690+my_id,*) nel,nel_unique,nel_local

	ierr = 0
	do ie=1,nel_unique
	  lmax = ilhv(ie)
	  do l=1,lmax
	    diff = abs(vals(l,ie)-valsaux(l,ie))
	    if( diff > 0 ) then
	      ierr = ierr + 1
	      write(670+my_id,*) ie,l,vals(l,ie),valsaux(l,ie),diff
	      write(680+my_id,*) ie,l,vals(l,ie),diff
	      write(690+my_id,*) ie,l,valsaux(l,ie),diff
	    end if
	  end do
	  diff = abs(vals2d(ie)-vals2daux(ie))
	  if( diff > 0 ) then
	    ierr = ierr + 1
	    write(660+my_id,*) ie,vals2d(ie),vals2daux(ie),diff
	  end if
	end do

	call shympi_check_3d_elem_r(vals,'tesssssssst')

	if( ierr > 0 ) then
	  stop 'error stop test_exchange_3d: internal (31)'
	end if

	end

!******************************************************************

	subroutine exchange_areas

	use mod_geom
	!use mod_internal
	use mod_hydro
	use evgeom
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

	integer ie,ii,iei,ia,ia0,iee
	real, parameter :: flag = -777.
	real, allocatable :: areat(:,:,:)
	real, allocatable :: areas(:,:)

	integer ieext

	allocate(areat(3,nel,n_threads))
	allocate(areas(3,nel))

	areas = flag

	do ie=1,nel
	  do ii=1,3
            iei = ieltv(ii,ie)
	    if( iei > 0 ) areas(ii,ie) = 12. * ev(10,iei)
	  end do
	end do

	call shympi_barrier
	call shympi_gather(3,areas,areat)

	ia0 = my_id + 1

	do ie=1,nel_unique
	  iee = ieext(ie)
	  do ii=1,3
	    if( areat(ii,ie,ia0) == flag ) then
	      do ia=1,n_threads
		if( ia == ia0 ) cycle
		if( areat(ii,ie,ia) /= flag ) then
		  write(6,*) ie,ii,ia0,ia,areat(ii,ie,ia)
		end if
	      end do
	    end if
	  end do
	end do

	stop 'error stop exchange_areas: generic'

	end

!******************************************************************

	subroutine check_id_elem(text)

	use basin
	use shympi

	implicit none

	character*(*) text

	integer ie,ii,id,n

	do ie=1,nel
	  n = id_elem(0,ie)
	  do ii=1,n
	    id = id_elem(ii,ie)
	    if( id < 0 .or. id >= n_threads ) then
	      write(6,*) trim(text)
	      write(6,*) 'error in id_elem check: ',ie
	      write(6,*) id_elem(:,ie)
	      call error_stop('check_id_elem','id_elem corrupt')
	    end if
	  end do
	end do

	end

!******************************************************************

