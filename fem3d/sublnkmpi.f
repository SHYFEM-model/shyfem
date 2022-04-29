
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2008-2015,2019  Georg Umgiesser
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

! topological set up routines - correct for mpi domains
!
! revision log :
!
! 13.04.2022	ggu	newly written
! 26.04.2022	ggu	write grd file if needed
!
!*****************************************************************

	subroutine set_geom_mpi

! sets up geometrical arrays
!
! geom structures have already been setup on local domain
! with this call they are corrected with gloval boundary node information

	use mod_geom
	use basin
	use shympi

        implicit none

        integer i,n
	integer k,kk,kn,kb,k1,k2
	integer ie,ii
	integer id,id_neigh
	integer iunit
	integer, allocatable :: ibound(:)
        integer kerr

	integer knext,kbhnd

	if( .not. bmpi ) return

	iunit = 660 + my_id

!-------------------------------------------------------------
! make global ibound array
!-------------------------------------------------------------

	call write_grd_domain

	!return
	!write(iunit,*) '-------- set_geom_mpi ----------'
	!write(iunit,*) 'global: ',nkn_global,nel_global
	!write(iunit,*) 'local: ',nkn,nel
	allocate(ibound(nkn_global))
	!call link_set_wmpi(.true.)
	call make_links(nkn_global,nel_global,nen3v_global,ibound,kerr)
	!call link_set_wmpi(.false.)
	!stop

!-------------------------------------------------------------
! reduce to local
!-------------------------------------------------------------

! ibound are global boundary nodes
! iboundv are local boundary nodes

	call shympi_barrier
	call shympi_g2l_array(ibound,iboundv)

!-------------------------------------------------------------
! now we adjust indices on local domain
!-------------------------------------------------------------

!-------------------------------------------------------------
! adjust kantv
!-------------------------------------------------------------

	do k=1,nkn
	  if( iboundv(k) == 0 ) then
	    kantv(:,k) = 0
	  else if( .not. shympi_is_inner_node(k) ) then
	    do i=1,2
	      kk = kantv(i,k)
	      id = id_node(kk)
	      if( id /= my_id ) kantv(i,k) = 0
	    end do
	  end if
	end do

!-------------------------------------------------------------
! adjust ieltv (-1 is reserved for OB)
!-------------------------------------------------------------

	do ie=1,nel
	  do ii=1,3
	    if( ieltv(ii,ie) /= 0 ) cycle	! has neighbor
	    k = nen3v(ii,ie)
	    kn = knext(k,ie)
	    kb = kbhnd(k,ie)
	    if( iboundv(kn) == 0 .or. iboundv(kb) == 0 ) then
	      id_neigh = id_elem(1,ie)
	      if( id_neigh == -1 ) goto 99
	      ieltv(ii,ie) = -1000 - id_neigh	!internal (not real) border
	    end if
	  end do
	end do

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	return
   99	continue
	stop 'error stop set_geom_mpi: internal error (1)'
	end

!*****************************************************************

	subroutine write_grd_domain

	use shympi
	use basin

	implicit none

	integer nout
	integer k,ie,kext,ieext,itype,n
	real x,y,depth
	real, parameter :: flag = -999.
	real, allocatable :: xg(:),yg(:)
	integer, allocatable :: intype(:),ietype(:)
	integer, allocatable :: ieaux(:)

	write(6,*) 'write_grd_domain:',my_id,nkn_global,size(id_node)

        allocate(xg(nkn_global))
        allocate(yg(nkn_global))
        allocate(intype(nkn_global))
        allocate(ietype(nel_global))
        allocate(ieaux(nel))

	ieaux(:) = id_elem(0,:)

	call shympi_exchange_array(xgv,xg)
	call shympi_exchange_array(ygv,yg)
	call shympi_exchange_array(id_node,intype)
	call shympi_exchange_array(ieaux,ietype)

	if( shympi_is_master() ) then

	nout = 1
	open(nout,file='domain.grd',status='unknown',form='formatted')

	depth = flag

	do k=1,nkn_global
	  kext = ip_ext_node(k)
	  itype = intype(k)
	  x = xg(k)
	  y = yg(k)
	  call grd_write_node(nout,kext,itype,x,y,depth)
	end do

	n = 3
	do ie=1,nel_global
	  ieext = ip_ext_elem(ie)
	  itype = ietype(ie)
          call grd_write_item(nout,2,ieext,itype,n,
     +                          nen3v_global(1,ie),ip_ext_node,depth)
	end do

	close(1)

	end if

	call shympi_barrier

	end

!*****************************************************************

