
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
! 19.05.2020    ccf     started from scratch
! 12.04.2022    ggu     adapted
!
!******************************************************************

!==================================================================
        module mod_flood
!==================================================================

        implicit none

!==================================================================
        contains
!==================================================================

!==================================================================
        end module mod_flood
!==================================================================

	subroutine floodfill_elem(nel,acode,ac,iareas)

	implicit none

	integer nel		!total number of elements
	integer acode(nel)	!area code of elements
	integer ac		!area code to look for
	integer iareas		!number of connected domains with area code ac

	integer ie
	integer aaux(nel)	!aux array to remember what has been filled

	iareas = 0
	aaux = 0

	do

	  do ie=1,nel
	    if( acode(ie) == ac .and. aaux(ie) == 0 ) exit
	  end do

	  if( ie > nel ) return
	  iareas = iareas + 1

	  call fill_elem_area(nel,acode,ie,aaux)

	end do

	end

!******************************************************************

	subroutine fill_elem_area(nel,acode,ie,aaux)

! flags aaux with 1 if can be reached from element ie and has acode of ie

	use queue
	use mod_connect

	implicit none

	integer nel		!total number of elements
	integer acode(nel)	!area code of elements
	integer ie		!element to start with
	integer aaux(nel)	!aux array to show what has been filled

	integer ac,ii,ieo,ien,id

	ac = acode(ie)

	call queue_init(id)
	call queue_enqueue(id,ie)

	do

	  if( .not. queue_dequeue(id,ieo) ) exit

	  aaux(ieo) = aaux(ieo) + 1
	  if( aaux(ieo) > 1 ) cycle

	  do ii=1,3
	    ien = ecv(ii,ieo)
	    if( ien > 0 ) then
	      if( acode(ien) == ac .and. aaux(ien) == 0 ) then
		call queue_enqueue(id,ien)
	      end if
	    end if
	  end do

	end do

	call queue_delete(id)

	end

!******************************************************************

	subroutine floodfill_node(nkn,acode,ac,iareas)

	implicit none

	integer nkn		!total number of nodes
	integer acode(nkn)	!area code of nodes
	integer ac		!area code to look for
	integer iareas		!number of connected domains with area code ac

	integer k
	integer aaux(nkn)	!aux array to remember what has been filled

	iareas = 0
	aaux = 0

	do

	  do k=1,nkn
	    if( acode(k) == ac .and. aaux(k) == 0 ) exit
	  end do

	  if( k > nkn ) return
	  iareas = iareas + 1

	  call fill_node_area(nkn,acode,k,aaux)

	end do

	end

!******************************************************************

	subroutine fill_node_area(nkn,acode,k,aaux)

! flags aaux with 1 if can be reached from element ie and has acode of ie

	use queue
	use mod_connect

	implicit none

	integer nkn		!total number of nodes
	integer acode(nkn)	!area code of nodes
	integer k		!node to start with
	integer aaux(nkn)	!aux array to show what has been filled

	integer ac,i,n,ko,kn,id
	integer, allocatable :: nen3v(:,:)

	ac = acode(k)

	call queue_init(id)
	call queue_enqueue(id,k)

	do

	  if( .not. queue_dequeue(id,ko) ) exit

	  aaux(ko) = aaux(ko) + 1
	  if( aaux(ko) > 1 ) cycle

	  n = nlist(0,ko)
	  do i=1,n
	    kn = nlist(i,ko)
	    if( acode(kn) == ac .and. aaux(kn) == 0 ) then
	      call queue_enqueue(id,kn)
	    end if
	  end do

	end do

	call queue_delete(id)

	end


