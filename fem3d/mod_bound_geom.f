
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

!==================================================================
        module mod_bound_geom
!==================================================================

c----------------------------------------------------------------------
c iopbnd(k) = 0         no open BC
c iopbnd(k) > 0         external open BC (ibtyp=1,2)
c iopbnd(k) < 0         internal open BC (ibtyp=3)
c----------------------------------------------------------------------

        implicit none

        integer, private, save :: nkn_bound_geom = 0
        integer, private, save :: nrb_bound_geom = 0

        integer, allocatable, save :: irv(:)
        integer, allocatable, save :: ierv(:,:)
        real, allocatable, save :: rhv(:)
        real, allocatable, save :: rlv(:)
        real, allocatable, save :: rrv(:)
        integer, allocatable, save :: iopbnd(:)

!==================================================================
	contains
!==================================================================

        subroutine mod_bound_geom_init(nkn,nrb)

        integer nkn
        integer nrb

	integer naux

	if( nkn == nkn_bound_geom .and. nrb == nrb_bound_geom ) return

        if( nkn_bound_geom > 0 ) then
          deallocate(irv)
          deallocate(ierv)
          deallocate(rhv)
          deallocate(rlv)
          deallocate(rrv)
          deallocate(iopbnd)
        end if

        nkn_bound_geom = nkn
        nrb_bound_geom = nrb

        if( nkn == 0 ) return

	naux = max(1,nrb)

        allocate(irv(naux))
        allocate(ierv(2,naux))
        allocate(rhv(naux))
        allocate(rlv(naux))
        allocate(rrv(naux))
        allocate(iopbnd(nkn))

	irv = 0
	ierv = 0
	rhv = 0.
	rlv = 0.
	rrv = 0.
	iopbnd = 0

        end subroutine mod_bound_geom_init

!***************************************************************

        subroutine mod_irv_init(nrb)

	integer nrb

	integer ndim,nnew
        integer, allocatable :: irv_aux(:)

	ndim = nrb_bound_geom
	nnew = ndim

        if( ndim == 0 ) then
          nnew = 10
          allocate(irv(nnew))
        else if( nrb > ndim ) then
          nnew = max(ndim*2,nrb)
          allocate(irv_aux(nnew))
          irv_aux(1:ndim) = irv(1:ndim)
          call move_alloc(irv_aux,irv)
        end if

	nrb_bound_geom = nnew

        end subroutine mod_irv_init

!***************************************************************

        subroutine mod_bound_geom_reinit(nkn,nrb)

        integer nkn
        integer nrb

	integer, allocatable :: irv_aux(:)

	if( nrb > nrb_bound_geom ) then
	  write(6,*) 'nrb,nrb_bound_geom: ',nrb,nrb_bound_geom
	  stop 'error stop mod_bound_geom_reinit: nrb > nrb_bound_geom'
	end if

	allocate(irv_aux(nrb))
	irv_aux(1:nrb) = irv(1:nrb)

	write(6,*) 'mod_bound_geom_reinit: ',nkn,nrb,nrb_bound_geom

	nrb_bound_geom = 0
	deallocate(irv)

        call mod_bound_geom_init(nkn,nrb)

	irv = 0
	irv(1:nrb) = irv_aux(1:nrb)
	deallocate(irv_aux)

	!write(6,*) 'mod_bound_geom_reinit: ',nrb
	!write(6,'(8i9)') irv

        end subroutine mod_bound_geom_reinit

!***************************************************************

        subroutine mod_bound_geom_info

	integer iu,i
	integer nkn,nrb

	iu = 88
	nkn = nkn_bound_geom
	nrb = nrb_bound_geom

	write(iu,*) 'mod_bound_geom_info: ',nkn,nrb
        write(iu,*) 'irv: ',nrb,(irv(i),i=1,nrb)
        write(iu,*) 'ierv: ',(ierv(1,i),i=1,nrb)
        write(iu,*) 'ierv: ',(ierv(2,i),i=1,nrb)
        write(iu,*) 'rhv: ',(rhv(i),i=1,nrb)
        write(iu,*) 'rlv: ',(rlv(i),i=1,nrb)
        write(iu,*) 'rrv: ',(rrv(i),i=1,nrb)
        write(iu,*) 'iopbnd: '
        do i=1,nkn
          if( iopbnd(i) .ne. 0 ) write(iu,*) i,iopbnd(i)
        end do
	write(iu,*) 'mod_bound_geom_info end'

        end subroutine mod_bound_geom_info

!***************************************************************

	function is_boundary(k)

	logical is_boundary
	integer k

	is_boundary = iopbnd(k) .ne. 0

	end function is_boundary

!***************************************************************

	function is_external_boundary(k)

	logical is_external_boundary
	integer k

	is_external_boundary = iopbnd(k) .gt. 0

	end function is_external_boundary

!***************************************************************

	function is_internal_boundary(k)

	logical is_internal_boundary
	integer k

	is_internal_boundary = iopbnd(k) .lt. 0

	end function is_internal_boundary

!***************************************************************

	function is_inner(k)

	logical is_inner
	integer k

	is_inner = iopbnd(k) .eq. 0

	end function is_inner

!==================================================================
        end module mod_bound_geom
!==================================================================

        subroutine mod_irv_initialize

	use mod_bound_geom

	implicit none

        call mod_irv_init(0)

	end subroutine mod_irv_initialize

!***************************************************************

