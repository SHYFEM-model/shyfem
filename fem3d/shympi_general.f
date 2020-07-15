
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
! general routines for mpi
!
! revision log :
!
! 07.12.2017	ggu	changed VERS_7_5_40
! 19.04.2018	ggu	changed VERS_7_5_45
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!***************************************************************

	subroutine mpi_sort_index(nk,ne)

	use basin
	use shympi

	implicit none

	integer nk,ne

	integer i,nold,nnew
	integer nks,nes,nki,nei
	logical, parameter :: bdebug = .false.

	nki = size(ipv)
	nei = size(ipev)
	nks = size(ip_sort_node)
	nes = size(ip_sort_elem)

	if( nk /= nki ) goto 97
	if( nk /= nks ) goto 97
	if( ne /= nei ) goto 97
	if( ne /= nes ) goto 97

        call isort(nk,ipv,ip_sort_node)
        call isort(ne,ipev,ip_sort_elem)

	if( bdebug ) then
	  write(6,*) 'mpi_sort_index: ',my_id
	  write(6,*) nk,ne
	  write(6,*) size(ip_sort_node),size(ip_sort_elem)
	  write(6,*) size(ipv),size(ipev)
	  write(6,*) 'ipv:',(ipv(i),i=1,nk,nk/10)
	  write(6,*) 'ip_sort_node:',(ip_sort_node(i),i=1,nk,nk/10)
	  write(6,*) 'mpi_sort_index: finished message'
	  flush(6)
	end if

        nold = ipv(ip_sort_node(1))
        do i=2,nk
          nnew = ipv(ip_sort_node(i))
          if( nnew <= nold ) goto 99
          nold = nnew
        end do
	if( any( ip_sort_node <= 0 ) ) goto 98

        nold = ipev(ip_sort_elem(1))
        do i=2,ne
          nnew = ipev(ip_sort_elem(i))
          if( nnew <= nold ) goto 99
          nold = nnew
        end do
	if( any( ip_sort_elem <= 0 ) ) goto 98

        return
   97   continue
        write(6,*) 'error in array dimensions:'
	write(6,*) nk,nki,nks
	write(6,*) ne,nei,nes
        stop 'error stop mpi_sort_index: array dimension'
   98   continue
        write(6,*) 'array has zero index: ip_sort_node,ip_sort_elem'
        stop 'error stop mpi_sort_index: zero index'
   99   continue
        write(6,*) 'array is not sorted: ',i,nold,nnew
        stop 'error stop mpi_sort_index: not sorted'
	end

!***************************************************************

