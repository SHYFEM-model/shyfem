
!--------------------------------------------------------------------------
!
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

c**************************************************************

! revision log :
!
! 24.01.2018	ggu	changed VERS_7_5_41
! 14.02.2019	ggu	changed VERS_7_5_56

c**************************************************************


	subroutine heap_sort(n,ra)

	implicit none

	integer n
	real ra(1)

	call heap_make(n,ra)
	call heap_check(n,ra)
	call heap_retire(n,ra)

	end

c******************************************************************

	subroutine heap_make(n,ra)

	implicit none

	integer n
	real ra(1)

	integer l

	if( n .lt. 2 ) return

	do l=n/2,1,-1
	  call heap_demote(l,n,ra)
	end do

	end

c******************************************************************

	subroutine heap_print(n,ra,text)

	implicit none

	integer n
	real ra(1)
	character*(*) text

	integer i

	write(6,*) 'heap print: ',text

	do i=1,n
	  write(6,*) i,n,ra(i)
	end do

	end

c******************************************************************

	subroutine heap_swap(i1,i2,ra)

	implicit none

	integer i1,i2
	real ra(1)

	real rra

	rra = ra(i1)
	ra(i1) = ra(i2)
	ra(i2) = rra

	end

c******************************************************************

	subroutine heap_check(n,ra)

	implicit none

	integer n
	real ra(1)

	integer i,j

	do i=1,n/2
	  j = i+i
	  if( ra(i) .lt. ra(j) ) goto 99
	  j = j + 1
	  if( j .le. n .and. ra(i) .lt. ra(j) ) goto 99
	end do

	return
   99	continue
	write(6,*) n,i,j,ra(i),ra(j)
	stop 'error stop heap_check: heap property violated'
	end

c******************************************************************

	subroutine heap_insert(n,ra,raa)

	implicit none

	integer n
	real ra(1)
	real raa

	n = n + 1
	ra(n) = raa

	call heap_promote(n,n,ra)

	end

c******************************************************************

	subroutine heap_retire(n,ra)

	implicit none

	integer n
	real ra(1)

	integer ir

	do ir=n,2,-1
	  call heap_swap(1,ir,ra)
	  call heap_demote(1,ir-1,ra)
	  !call heap_print(n,ra,'retire')
	end do

	end

c******************************************************************

	subroutine heap_adjust(l,n,ra)

c adjusts entry l to right place

	implicit none

	integer l,n
	real ra(1)

	call heap_promote(l,n,ra)
	call heap_demote(l,n,ra)

	end

c******************************************************************

	subroutine heap_promote(l,n,ra)

c promotes entry l to right place

	implicit none

	integer l,n
	real ra(1)

	integer i,j
	real rra

	i = l		!adjust this index
	j = l/2		!node to compare with
	rra = ra(i)	!value to promote

	do while( j .ge. 1 )

	  if( rra .gt. ra(j) ) then	!demote rra
	    ra(i) = ra(j)
	    i = j
	    j = j/2
	  else				!finished
	    j = 0
	  end if

	end do

	ra(i) = rra

	end

c******************************************************************

	subroutine heap_demote(l,n,ra)

c demotes entry l to right place

	implicit none

	integer l,n
	real ra(1)

	integer i,j
	real rra

	i = l		!adjust this index
	j = l+l		!node to compare with
	rra = ra(i)	!value to demote

	do while( j .le. n )

	  if( j .lt. n .and. ra(j) .lt. ra(j+1) ) j = j + 1

	  if( rra .lt. ra(j) ) then	!demote rra
	    ra(i) = ra(j)
	    i = j
	    j = j + j
	  else				!finished
	    j = n + 1
	  end if

	end do

	ra(i) = rra

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine heap_test

	implicit none

	integer ndim
	!parameter (ndim=200000)
	!parameter (ndim=20)
	parameter (ndim=200)

	real ra(ndim)
	real ra1(ndim)
	real ra2(ndim)
	real ra3(ndim)
	real raa
	integer i,n,j
	logical berror

	real heap_random

	j = 0
	n = ndim
	berror = .false.

	do i=1,n
	  call random_number(raa)
	  ra(i) = raa
	  ra1(i) = ra(i)
	  ra2(i) = ra(i)
	  !write(6,*) i,ra(i)
	  call heap_insert(j,ra3,raa)
	  call heap_check(j,ra3)
	end do

	call heap_sort(n,ra1)
	!call sort(n,ra2)
	call heap_check(n,ra3)
	call heap_retire(n,ra3)

	!write(6,*) 1,ra1(1),ra2(1)
	do i=2,n
	  if( ra1(i) .lt. ra1(i-1) ) then
	    !write(6,*) 'not sorted... ',i,ra1(i),ra1(i-1)
	    berror = .true.
	  end if
	  if( ra1(i) .ne. ra2(i) ) then
	    !write(6,*) 'not equal... ',i,ra1(i),ra2(i)
	    berror = .true.
	  end if
	  if( ra1(i) .ne. ra3(i) ) then
	    !write(6,*) 'not equal... ',i,ra1(i),ra2(i)
	    berror = .true.
	  end if
	  !write(6,*) i,ra1(i),ra2(i),ra3(i)
	end do

	if( berror ) then
	  write(6,*) 'there have been errors...   n = ',n
	else
	  write(6,*) 'test passed.   n = ',n
	end if

	end
	  
c******************************************************************

        subroutine heap_rand_int(min,max,irand)

        implicit none
        integer min,max
        integer irand
        real r
        call random_number(r)
        irand = min + (max-min+1)*r
        end

c******************************************************************
c******************************************************************
c******************************************************************
	program heap_main
	call heap_test
	end
c******************************************************************

