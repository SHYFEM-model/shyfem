
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2003-2004,2008-2009,2011,2013  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
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
! 06.04.1999	ggu	completely restructured
! 04.06.1999	ggu	new statistics are computed
! 08.09.2003	ggu	mode 5 -> write depth values from elements
! 23.09.2004	ggu	interpolq() changed for bathy interpolation
! 02.10.2004	ggu	interpole() for exponential interpolation
! 01.11.2004	ggu	whole program simplyfied
! 06.12.2008	ggu	smoothing introduced
! 06.04.2009	ggu	read param.h
! 29.05.2009	ggu	does only depth limiting and smoothing
! 20.11.2009	ggu	possibility to smooth only on specific areas
! 30.03.2011	ggu	new routines to delete elements
! 13.06.2013	ggu	copy_depth() renamed to transfer_depth()
! 16.04.2018	ggu	partitioning finished
! 20.04.2018	ggu	code added to check partition integrity
! 16.02.2019	ggu	changed VERS_7_5_60
! 28.05.2020	ggu	new checks for connection
!
!****************************************************************

!================================================================
	module mod_save_index
!================================================================

	implicit none

	integer, save, private :: nkn_save
	integer, save, private :: nel_save

	integer, allocatable :: save_ipv(:)
	integer, allocatable :: save_ipev(:)
	integer, allocatable :: save_nen3v(:,:)
	integer, allocatable :: aux_ipv(:)
	integer, allocatable :: aux_ipev(:)

!================================================================
	contains
!================================================================

	subroutine mod_save_index_init(nkn,nel)

	implicit none

        integer nkn, nel

        if( nkn == nkn_save .and. nel == nel_save ) return

        if( nel > 0 .or. nkn > 0 ) then
          if( nel == 0 .or. nkn == 0 ) then
            write(6,*) 'nel,nkn: ',nel,nkn
            stop 'error stop mod_save_index: incompatible parameters'
          end if
        end if

        if( nkn_save > 0 ) then
	  deallocate(save_ipv)
	  deallocate(save_ipev)
	  deallocate(save_nen3v)
	  deallocate(aux_ipv)
	  deallocate(aux_ipev)
        end if

	nkn_save = nkn
	nel_save = nel

	if( nkn == 0 ) return

	allocate(save_ipv(nkn))
	allocate(save_ipev(nel))
	allocate(save_nen3v(3,nel))
	allocate(aux_ipv(nkn))
	allocate(aux_ipev(nel))

	end subroutine mod_save_index_init

!****************************************************************

	subroutine make_new_index(nk,ne,nenv,nodep,elemp)

	use basin

	implicit none

	integer nk,ne
	integer nenv(3,nel)
	integer nodep(nk)
	integer elemp(ne)

	integer k,ie

	if( nk > nkn .or. ne > nel ) then
	  write(6,*) 'can only make smaller index: ',nk,nkn,ne,nel
	  stop 'error stop make_new_index: nk/ne > nkn/nel'
	end if

	call mod_save_index_init(nkn,nel)

	save_ipv = ipv
	save_ipev = ipev
	save_nen3v = nen3v

	nen3v(:,1:ne) = nenv(:,1:ne)

	do k=1,nk
	  aux_ipv(k) = ipv(nodep(k))
	end do
	ipv(1:nk) = aux_ipv(1:nk)

	do ie=1,ne
	  aux_ipev(ie) = ipev(elemp(ie))
	end do
	ipev(1:ne) = aux_ipev(1:ne)

	end subroutine

!****************************************************************

	subroutine restore_old_index

	use basin

	implicit none

	nen3v = save_nen3v
	ipv = save_ipv
	ipev = save_ipev

	end subroutine

!================================================================
	end module mod_save_index
!================================================================

        subroutine bas_partition

! performs partition on basin

	use mod_geom
	use mod_depth
	use evgeom
	use basin
	use grd
	use basutil

	implicit none

	integer ierr
	integer k,i,nl,il,n,ib,in,node,ilext,np,ic
	integer, allocatable :: nc(:)
	real x,y,perc
	real xx(nkn)
	real yy(nkn)

	logical inconvex,inpoly,filex

!-----------------------------------------------------------------
! open and read file containing lines
!-----------------------------------------------------------------

	if( .not. filex(lfile) ) then
	  write(6,*) 'Cannot open file ',trim(lfile)
	  stop 'error stop bas_partition: no such file'
	end if

	call grd_read(lfile)

	iarnv = 0
	nl = nl_grd

!-----------------------------------------------------------------
! loop over lines
!-----------------------------------------------------------------

        do i=1,nl
          il = i
          ilext = ipplv(il)
          n = ipntlv(il) - ipntlv(il-1)
          ib = ipntlv(il-1)
	  !write(6,*) 'line: ',i,ilext,n
	  do in=1,n
	    node = inodlv(ib+in)
	    x = xv(node)
	    y = yv(node)
	    xx(in) = x
	    yy(in) = y
	    !write(6,*) in,node,x,y
	  end do
	  if( xx(1) == xx(n) .and. yy(1) == yy(n) ) n = n - 1
	  np = 0
	  do k=1,nkn
	    x = xgv(k)
	    y = ygv(k)
	    !write(6,*) 'looking for... ',k,x,y
	    if( inpoly(n,xx,yy,x,y) ) then
	      iarnv(k) = il
	      np = np + 1
	      !write(6,*) 'inside... ',il,n,k,x,y
	    end if
	  end do
	  perc = (100.*np)/nkn
	  !write(6,*) 'inside points found... ',i,il,n,np,perc
        end do

	if( count( iarnv == 0 ) > 0 ) then	!not handled nodes
	  nl = nl + 1
	  where( iarnv == 0 ) iarnv = nl
	end if

!-----------------------------------------------------------------
! write information to terminal
!-----------------------------------------------------------------

	allocate(nc(0:nl))
	nc = 0
	do k=1,nkn
	  ic = iarnv(k)
	  if( ic < 1 .or. ic > nl ) then
	    write(6,*) 'ic,nl: ',ic,nl
	    stop 'error stop bas_partition: internal error (1)'
	  end if
	  nc(ic) = nc(ic) + 1
	end do
	write(6,*) 'Information on domains: '
	write(6,*) '   domain     nodes   percent'
	do ic=1,nl
	  write(6,'(2i10,f10.2)') ic,nc(ic),(100.*nc(ic))/nkn
	end do

!-----------------------------------------------------------------
! write file
!-----------------------------------------------------------------

        call basin_to_grd
        call grd_write('bas_partition.grd')
        write(6,*) 'The partition has been written to bas_partition.grd'

	call check_connectivity(ierr)
	call check_connections(ierr)

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
	end

!*******************************************************************

	subroutine check_connectivity(ierr)

	use basin

	implicit none

	integer ierr

	integer ier
	integer ic,nc
	integer icolor(nkn)

	ierr = 0
	icolor = iarnv
	nc = maxval(icolor)

	do ic=0,nc
	  call check_color(ic,nkn,icolor,ier)
	  ierr = ierr + ier
	end do

	end

!*******************************************************************

	subroutine check_color(ic,n,icolor,ierr)

	use basin

	implicit none

	integer ic
	integer n
	integer icolor(n)
	integer ierr

	integer cc,i,nfound

	ierr = 0
	cc = count( icolor == ic )

	do
	  do i=1,n
	    if( icolor(i) == ic ) exit		!start from this node
	  end do
	  if( i > n ) exit
	  call flood_fill_bas(i,n,icolor,nfound)
	  !write(6,*) ic,nfound,(100.*nfound)/n,' %'
	  if( nfound /= cc ) then
	    write(6,*) '  *** area is not connected...'
	    write(6,*) '      area code:     ',ic
	    write(6,*) '      contains node: ',ipv(i)
	    ierr = 1
	  end if
	end do

	end

!*******************************************************************

	subroutine flood_fill_bas(i,n,icolor,nfound)

	use basin

	implicit none

	integer i
	integer n
	integer icolor(n)
	integer nfound

	integer ic,nf,ie,ii,k
	logical bcol,bdone

	ic = icolor(i)
	icolor(i) = -1
	nfound = 1

	do
	  nf = 0
	  do ie=1,nel
	    bcol = .false.
	    bdone = .false.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( icolor(k) == ic ) bcol = .true.
	      if( icolor(k) == -1 ) bdone = .true.
	    end do
	    if( bcol .and. bdone ) then
	      do ii=1,3
	        k = nen3v(ii,ie)
	        if( icolor(k) == ic ) then
	          nf = nf + 1
		  icolor(k) = -1
		end if
	      end do
	    end if
	  end do
	  if( nf == 0 ) exit
	  nfound = nfound + nf
	end do
	
	where( icolor == -1 ) icolor = -2

	end

!*******************************************************************

	subroutine check_connections(kerr)

! this checks connections with link/lenk data structure

	use basin
	use mod_save_index

	implicit none

	integer kerr

	logical bloop
	logical bwrite
	integer nloop
	integer ic,nc,ncol,kext
	integer nk,ne
	integer nenv(3,nel)
	integer icolor(nkn)
	integer icol(nkn)
	integer nodep(nkn)
	integer elemp(nel)

	integer ipint,ipext

	bwrite = .false.
	bloop = .true.
	nloop = 0
	icolor = iarnv
	nc = maxval(icolor)

	write(6,*) '========================================'
	write(6,*) 'checking total domain... total domains = ',nc
	write(6,*) '========================================'

	call check_elem_index(nkn,nel,nen3v,kerr)

!---------------------------------------------
! loop on domains
!---------------------------------------------

	do while( bloop )

	nloop = nloop + 1
	if( nloop > 10 ) exit
	!if( nloop > 1 ) exit

	do ic=1,nc
	  ncol = count( icolor == ic )
	  if( ncol == 0 ) cycle

	  !write(6,*) '========================================'
	  if( bwrite ) write(6,*) 'checking domain ',ic,ncol
	  !write(6,*) '========================================'

	  call make_elem_index(.true.,ic,icolor
     +				,nk,ne,nenv,nodep,elemp)
	  call make_new_index(nk,ne,nenv,nodep,elemp)
	  call check_elem_index(nk,ne,nenv,kerr)

	  !-------------------------------------------
	  ! handle errors
	  !-------------------------------------------

	  if( kerr /= 0 ) then
	    kext = ipext(kerr)
	    !write(6,*) 'adjusting node kerr = ',kerr,kext,ic
	    call restore_old_index
	    call adjust_domain(ic,nkn,icolor,kext)
	    exit
	  end if

	  call restore_old_index
	end do

	  if( ic > nc ) bloop = .false.
	end do

!---------------------------------------------
! end of loop on domains
!---------------------------------------------

	if( kerr == 0 ) then
	  if( nloop > 1 ) then
	    write(6,*) 'all domains have been corrected...',nloop
	  else
	    write(6,*) 'no problems found in domains'
	  end if
	else
	  write(6,*) 'could not correct error in connections...',nloop
	end if

	iarnv = icolor

!---------------------------------------------
! end of routine
!---------------------------------------------

	end

!*******************************************************************

	subroutine make_elem_index(bghost,ic,icolor
     +				,nk,ne,nenv,nodep,elemp)

! makes element index for domain with color ic

	use basin

	implicit none

	logical bghost		!include ghost nodes (and elements)
	integer ic
	integer icolor(nkn)
	integer nk,ne
	integer nenv(3,nel)
	integer nodep(nkn)
	integer elemp(nel)

	integer ie,ii,k,nic
	integer node_aux(nkn)

	ne = 0
	node_aux = 0
	nodep = 0
	elemp = 0

	do ie=1,nel
	  nic = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( icolor(k) == ic ) nic = nic + 1
	  end do
	  if( .not. bghost ) nic = nic - 2	!only good if nic == 3
	  if( nic > 0 ) then	!color found
	    ne = ne + 1
	    nenv(:,ne) = nen3v(:,ie)
	    elemp(ne) = ie
	  end if
	end do
	    
	do ie=1,ne
	  do ii=1,3
	    k = nenv(ii,ie)
	    node_aux(k) = 1
	  end do
	end do

	nk = 0

	do k=1,nkn
	  if( node_aux(k) > 0 ) then
	    nk = nk + 1
	    node_aux(k) = nk
	    nodep(nk) = k
	  end if
	end do

	do ie=1,ne
	  do ii=1,3
	    k = nenv(ii,ie)
	    nenv(ii,ie) = node_aux(k)
	  end do
	end do

	end

!*******************************************************************

	subroutine check_elem_index(nk,ne,nenv,kerr)

! checks the element structure
!
! in order to get correct error feed back, the original arrays
! nen3v, ipv, ipev have to be altered
! they are saved at the beginning and then restored at the end

	use basin
	use mod_geom

	implicit none

	integer nk,ne
	integer nenv(3,ne)
	integer kerr

	integer k,ie,ngrm

	call estimate_max_grade(nk,ne,nenv,ngrm)
	call mod_geom_init(nk,ne,ngrm)

	!call make_links_old(nk,ne,nenv)
	call make_links(nk,ne,nenv,kerr)

	end

!*******************************************************************

	subroutine translate_color(nkn,nk,nodep,icolor,icol)

	implicit none

	integer nkn,nk
	integer nodep(nkn)
	integer icolor(nkn)
	integer icol(nk)

	integer k

	do k=1,nk
	  icol(k) = icolor(nodep(k))
	end do

	end

!*******************************************************************

	subroutine adjust_domain(ic,nkn,icolor,kerr)

	implicit none

	integer ic
	integer nkn
	integer icolor(nkn)
	integer kerr

	integer k,nc,icc,kint
	integer, allocatable :: count(:)

	integer ipext,ipint

	kint = ipint(kerr)

	nc = maxval(icolor)
	allocate(count(0:nc))
	count = 0

	do k=1,nkn
	  icc = icolor(k)
	  count(icc) = count(ic) + 1
	end do

	!write(6,*) 'colors for domain ',ic,icolor(kint)
	!do icc=1,nc
	!  write(6,*) icc,count(icc)
	!end do

	write(6,*) 're-coloring: ',kerr,kint,icolor(kint),ic

	icolor(kint) = ic

	end

!*******************************************************************

	subroutine count_elements(nkn,nel,nen3v,ic,icolor,netot,neint)

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	integer ic
	integer icolor(nkn)
	integer netot,neint

	integer k,ii,ie,icc,icount

	netot = 0
	neint = 0

	do ie=1,nel
	  icount = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    icc = icolor(k)
	    if( icc == ic ) icount = icount + 1
	  end do
	  if( icount > 0 )  netot = netot + 1
	  if( icount == 3 ) neint = neint + 1
	end do

	end

!*******************************************************************

