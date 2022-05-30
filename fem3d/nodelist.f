
!--------------------------------------------------------------------------
!
!    Copyright (C) 2022  Georg Umgiesser
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

! subroutines to split node list into single blocks
! blocks are separated by 0
!
! revision log :
!
! 22.05.2022    ggu     written from scratch

!***************************************************************

	subroutine list_get_number_of_blocks(n,llist,ns,nn)

! counts number of blocks in llist

	implicit none

	integer n		!total number of nodes in llist
	integer llist(n)	!original linear list
	integer ns		!total number of blocks found (return)
	integer nn		!maximum number of nodes in one block (return)

	integer i,is,nm

	is = 0
        i = 1
	nm = 0
	nn = 0

        do while( i .le. n )

          do while( i <= n .and. llist(i) == 0 )	!skip leading zeros
            i = i + 1
          end do
	  if( i > n ) exit

	  nm = 0 
          do while( i <= n .and. llist(i) /= 0 )	!look for next zeros
            i = i + 1
	    nm = nm + 1
          end do
	  is = is + 1
	  nn = max(nn,nm)
	end do

	ns = is

	end

!***************************************************************

	subroutine list_split_blocks(ns,n,nn,llist,slist)

! separates linear list into single blocks
!
! in slist(0,is) is total number of nodes of block is

	implicit none

	integer ns		!total number of sections
	integer n		!total number of nodes in llist
	integer nn		!maximum number of nodes in one block
	integer llist(n)	!original linear list
	integer slist(0:nn,ns)	!single blocks (return)

	integer i,is,nt,istart,iend

	is = 0
        i = 1

        do while( i .le. n )

          do while( i <= n .and. llist(i) == 0 )	!skip leading zeros
            i = i + 1
          end do
	  if( i > n ) exit
	  istart = i

          do while( i <= n .and. llist(i) /= 0 )	!look for next zeros
            i = i + 1
          end do
	  iend = i - 1
	  is = is + 1
	  nt = iend - istart + 1
	  if( nt > nn ) goto 99
	  slist(0,is) = nt
	  slist(1:nt,is) = llist(istart:iend)

        end do

	return
   99	continue
	write(6,*) 'nodes in block higher than dimension: ',nt,nn
	stop 'error stop list_split_blocks: nt>nn'
	end

!***************************************************************

	subroutine list_join_blocks(ns,n,nn,slist,llist)

! joins single blocks into linear list
!
! n is dimension of llist on entry, actual filling on return

	implicit none

	integer ns		!total number of sections
	integer n		!total number of nodes in llist (with zeros)
	integer nn		!maximum number of nodes in one block
	integer slist(0:nn,ns)	!single blocks
	integer llist(n)	!linear  list

	integer i,is,nt

	llist = 0

	i = 0
	do is=1,ns
	  nt = slist(0,is)
	  if( i + nt > n ) goto 99
	  llist(i+1:i+nt) = slist(1:nt,is)
	  i = i + nt + 1			!leave one space for zero
	end do

	n = i - 1
	if( n < 0 ) n = 0			!no list

	return
   99	continue
	write(6,*) 'nodes in block higher than dimension: ',i+nt,n
	stop 'error stop list_join_blocks: number of nodes too high'
	end

!***************************************************************
!***************************************************************
!***************************************************************
! test routines
!***************************************************************
!***************************************************************
!***************************************************************

	subroutine test_single_list(n,llist)

	implicit none

	integer n
	integer llist(n)
	integer, allocatable :: slist(:,:)

	integer ns,is,nt,nn,nr

	call list_get_number_of_blocks(n,llist,ns,nn)

	write(6,*) '================================='
	write(6,*) 'list has sections: ',ns,nn
	write(6,*) 'orig list:'
	write(6,1000) n
	write(6,1000) llist(1:n)

	allocate(slist(0:nn,ns))
	call list_split_blocks(ns,n,nn,llist,slist)

	write(6,*) 'separated list:'
	do is=1,ns
	  nt = slist(0,is)
	  write(6,1000) is,nt,slist(1:nt,is)
	end do

	write(6,*) 'final list:'
	nr = n
	call list_join_blocks(ns,nr,nn,slist,llist)
	write(6,1000) nr
	write(6,1000) llist(1:nr)

 1000	format(20i4)
	end

!***************************************************************

	subroutine test_list

	implicit none

	integer, parameter :: ndim = 10
	integer llist(ndim)
	integer slist(0:ndim,ndim)

	integer n,ns,is,nt

	n = ndim
	llist = (/1,2,3,0,4,5,6,0,7,8/)
	call test_single_list(n,llist)
	llist = (/0,1,2,0,0,4,5,6,7,0/)
	call test_single_list(n,llist)

	end

!***************************************************************
!	program main_test_list
!	call test_list
!	end
!***************************************************************

