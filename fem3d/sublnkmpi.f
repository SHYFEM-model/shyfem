
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

	return
	!write(iunit,*) '-------- set_geom_mpi ----------'
	!write(iunit,*) 'global: ',nkn_global,nel_global
	!write(iunit,*) 'local: ',nkn,nel
	allocate(ibound(nkn_global))
	call make_links(nkn_global,nel_global,nen3v_global,ibound,kerr)

!-------------------------------------------------------------
! reduce to local
!-------------------------------------------------------------

	call shympi_barrier
	call shympi_g2l_array(ibound,iboundv)

!-------------------------------------------------------------
! adjust kantv
!-------------------------------------------------------------

	do k=1,nkn
	  if( ibound(k) == 0 ) then
	    kantv(:,k) = 0
	  else if( .not. shympi_is_inner_node(k) ) then
	    do i=1,2
	      kk = kantv(i,k)
	      id = id_node(k)
	      if( id /= my_id ) kantv(i,k) = -id - 1	!-ia
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
	    if( ibound(kn) == 0 .or. ibound(kb) == 0 ) then
	      id_neigh = id_elem(1,ie)
	      if( id_neigh == -1 ) goto 99
	      ieltv(ii,ie) = -1000 - id_neigh	!internal (not real) border
	    end if
	  end do
	end do

!-------------------------------------------------------------
! compute dxv,dyv
!-------------------------------------------------------------

	dxv = 0.
	dyv = 0.
	do k=1,nkn
	  if( ibound(k) /= 0 .and. shympi_is_inner_node(k) ) then
	    k1 = kantv(1,k)
	    k2 = kantv(1,k)
	    dxv(k) = xgv(k1) - xgv(k2)
	    dyv(k) = ygv(k1) - ygv(k2)
	  end if
	end do

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	return
   99	continue
	stop 'error stop set_geom_mpi: internal error (1)'
	end

!*****************************************************************

