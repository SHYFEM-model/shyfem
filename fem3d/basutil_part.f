
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
!
!****************************************************************

        subroutine bas_partition

! performs partition on basin

	use mod_geom
	use mod_depth
	use evgeom
	use basin
	use grd
	use basutil

	implicit none

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

	call check_connectivity
	call check_connections

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
	end

!*******************************************************************

	subroutine check_connectivity

	use basin

	implicit none

	integer ic,nc
	integer icolor(nkn)

	icolor = iarnv
	nc = maxval(icolor)

	do ic=0,nc
	  call check_color(ic,nkn,icolor)
	end do

	end

!*******************************************************************

	subroutine check_color(ic,n,icolor)

	use basin

	implicit none

	integer ic
	integer n
	integer icolor(n)

	integer cc,i,nfound

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

	subroutine check_connections

! this checks connections with link/lenk data structure

	use basin

	implicit none

	integer ic,nc,ncol
	integer nk,ne
	integer nenv(3,nel)
	integer icolor(nkn)
	integer nodep(nkn)
	integer elemp(nel)

	icolor = iarnv
	nc = maxval(icolor)

	do ic=0,nc
	  ncol = count( icolor == ic )
	  if( ncol == 0 ) cycle
	  write(6,*) 'checking domain ',ic,ncol
	  call make_elem_index(.true.,ic,icolor
     +				,nk,ne,nenv,nodep,elemp)
	  call check_elem_index(nk,ne,nenv,nodep,elemp)
!	  next check is probably too restrictive
!	  call make_elem_index(.false.,ic,icolor
!     +				,nk,ne,nenv,nodep,elemp)
!	  call check_elem_index(nk,ne,nenv,nodep,elemp)
	end do

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

	subroutine check_elem_index(nk,ne,nenv,nodep,elemp)

! checks the element structure
!
! in order to get correct error feed back, the original arrays
! nen3v, ipv, ipev have to be altered
! they are saved at the beginning and then restored at the end

	use basin

	implicit none

	integer nk,ne
	integer nenv(3,ne)
	integer nodep(nk)
	integer elemp(ne)

	integer k,ie
	integer nlkdi
	integer, allocatable :: ilinkv(:)
	integer, allocatable :: lenkv(:)
	integer, allocatable :: lenkiiv(:)
	integer, allocatable :: linkv(:)
	integer, allocatable :: kantv(:,:)
	integer, allocatable :: ieltv(:,:)
	integer, allocatable :: save_ipv(:)
	integer, allocatable :: save_ipev(:)
	integer, allocatable :: save_nen3v(:,:)
	integer, allocatable :: aux_ipv(:)
	integer, allocatable :: aux_ipev(:)

        nlkdi = 3*ne+2*nk

	allocate(ilinkv(0:nk))
	allocate(lenkv(nlkdi))
	allocate(lenkiiv(nlkdi))
	allocate(linkv(nlkdi))
	allocate(kantv(2,nk))
	allocate(ieltv(3,ne))
	allocate(save_ipv(nkn))
	allocate(save_ipev(nel))
	allocate(save_nen3v(3,nel))
	allocate(aux_ipv(nk))
	allocate(aux_ipev(ne))

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

	!write(6,*) 'check_elem_index: ',nlkdi,nk,ne
        call mklenk(nlkdi,nk,ne,nenv,ilinkv,lenkv)
        call mklenkii(nlkdi,nk,ne,nenv,ilinkv,lenkv,lenkiiv)
        call mklink(nk,ilinkv,lenkv,linkv)

        call mkkant(nk,ilinkv,lenkv,linkv,kantv)
        call mkielt(nk,ne,ilinkv,lenkv,linkv,ieltv)

	ipv = save_ipv
	ipev = save_ipev
	nen3v = save_nen3v

	end

!*******************************************************************



















