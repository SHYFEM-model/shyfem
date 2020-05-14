
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017,2019  Georg Umgiesser
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

! stack utility routines
!
! still to do: implement stack for string
!
! revision log :
!
! 26.09.2017	ggu	changed VERS_7_5_32
! 05.12.2017	ggu	changed VERS_7_5_39
! 14.02.2019	ggu	changed VERS_7_5_56

!===============================================================
	module stack
!===============================================================

	implicit none

	private

        type :: entry
          integer :: top
          integer :: max
          integer :: type
          double precision, allocatable :: array(:)
          character*80, allocatable :: string(:)
        end type entry

        integer, parameter ::     no_type = 0
        integer, parameter ::  value_type = 1
        integer, parameter :: string_type = 2

        integer, parameter :: empty_error = 1
        integer, parameter ::  type_error = 2

        integer, save :: idlast = 0
        integer, save :: ndim = 0
        integer, parameter :: ndim_first = 10
        type(entry), save, allocatable :: pentry(:)

	public :: stack_init		!call stack_init(id)
	public :: stack_delete		!call stack_delete(id)
	public :: stack_push		!call stack_push(id,value)
	public :: stack_pop		!logical stack_pop(id,value)
	public :: stack_peek		!logical stack_peek(id,value)
	public :: stack_is_empty	!logical stack_is_empty(id)
	public :: stack_info		!call stack_info(id)

        INTERFACE stack_push
        MODULE PROCEDURE         stack_push_d
     +                          ,stack_push_r
     +                          ,stack_push_i
        END INTERFACE

        INTERFACE stack_pop
        MODULE PROCEDURE         stack_pop_d
     +                          ,stack_pop_r
     +                          ,stack_pop_i
        END INTERFACE

        INTERFACE stack_peek
        MODULE PROCEDURE         stack_peek_d
     +                          ,stack_peek_r
     +                          ,stack_peek_i
        END INTERFACE

!===============================================================
	contains
!===============================================================

        subroutine stack_init_alloc

        type(entry), allocatable :: paux(:)

        if( ndim == 0 ) then
          ndim = ndim_first
          allocate(pentry(ndim))
        else
          ndim = ndim*2
          allocate(paux(ndim))
          paux(1:ndim/2) = pentry(1:ndim/2)
          call move_alloc(paux,pentry)
        end if

        end subroutine stack_init_alloc

!******************************************************************

        subroutine stack_init_new_id(id)

        integer id

        idlast = idlast + 1
        if( idlast > ndim ) then
          call stack_init_alloc
        end if
        id = idlast

        call stack_init_id(id)

        end subroutine stack_init_new_id

!******************************************************************

        subroutine stack_init_id(id)

        integer id

        if( id > ndim ) then
          stop 'error stop stack_init_id: ndim'
        end if

        pentry(id)%top = 0
        pentry(id)%max = 0
        pentry(id)%type = 0

	if( allocated(pentry(id)%array) ) deallocate(pentry(id)%array)
	if( allocated(pentry(id)%string) ) deallocate(pentry(id)%string)

        end subroutine stack_init_id

!******************************************************************

	subroutine stack_error(id,error)

	integer id,error

	if( error == empty_error ) then
	  write(6,*) 'stack: ',id
	  stop 'error stop stack: stack is empty'
	else if( error == type_error ) then
	  write(6,*) 'stack: ',id
	  write(6,*) 'type: ',pentry(id)%type
	  stop 'error stop stack: variable is of wrong type'
	else
	  stop 'error stop stack: internal error (1)'
	end if

	end subroutine stack_error

!******************************************************************

	subroutine realloc_double(n,value)

	integer n
	double precision, allocatable :: value(:)

	integer nsize
	double precision, allocatable :: daux(:)

	if( n == 0 ) then
	  n = 10
	  allocate(value(n))
	else
	  nsize = min(n,size(value))
          allocate(daux(n))
          daux(1:nsize) = value(1:nsize)
          call move_alloc(daux,value)
	end if

	end subroutine realloc_double

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine stack_init(id)
	integer id
        call stack_init_new_id(id)
	end subroutine stack_init

	subroutine stack_delete(id)
	integer id
        call stack_init_id(id)
	if( id == idlast ) idlast = idlast - 1
	end subroutine stack_delete

!--------------------

	subroutine stack_push_i(id,value)
	integer id
	integer value
	call stack_push_d(id,dble(value))
	end subroutine stack_push_i

	subroutine stack_push_r(id,value)
	integer id
	real value
	call stack_push_d(id,dble(value))
	end subroutine stack_push_r

	subroutine stack_push_d(id,value)
	integer id
	double precision value
	integer n
	if( pentry(id)%top >= pentry(id)%max ) then
	  n = 2 * pentry(id)%max
	  call realloc_double(n,pentry(id)%array)
	  pentry(id)%max = n
	end if
	if( pentry(id)%type == no_type ) then
	  pentry(id)%type = value_type
	end if
	if( pentry(id)%type /= value_type ) then
	  call stack_error(id,type_error)
	end if
	pentry(id)%top = pentry(id)%top + 1
	pentry(id)%array(pentry(id)%top) = value
	end subroutine stack_push_d

!--------------------

	logical function stack_pop_i(id,value)
	integer id
	integer value
	double precision dvalue
	stack_pop_i = stack_pop_d(id,dvalue)
	value = nint(dvalue)
	end function stack_pop_i

	logical function stack_pop_r(id,value)
	integer id
	real value
	double precision dvalue
	stack_pop_r = stack_pop_d(id,dvalue)
	value = real(dvalue)
	end function stack_pop_r

	logical function stack_pop_d(id,value)
	integer id
	double precision value
	stack_pop_d = stack_peek_d(id,value)
	if( stack_pop_d ) pentry(id)%top = pentry(id)%top - 1
	end function stack_pop_d

!--------------------

	logical function stack_peek_i(id,value)
	integer id
	integer value
	double precision dvalue
	stack_peek_i = stack_peek_d(id,dvalue)
	value = nint(dvalue)
	end function stack_peek_i

	logical function stack_peek_r(id,value)
	integer id
	real value
	double precision dvalue
	stack_peek_r = stack_peek_d(id,dvalue)
	value = real(dvalue)
	end function stack_peek_r

	logical function stack_peek_d(id,value)
	integer id
	double precision value
	stack_peek_d = .false.
	if( pentry(id)%top == 0 ) return
	if( pentry(id)%type /= value_type ) then
	  call stack_error(id,type_error)
	end if
	value = pentry(id)%array(pentry(id)%top)
	stack_peek_d = .true.
	end function stack_peek_d

!--------------------

	logical function stack_is_empty(id)
	integer id
	stack_is_empty = ( pentry(id)%top == 0 )
	end function stack_is_empty

!--------------------

	subroutine stack_info(id)
	integer id
	write(6,*) 'stack_info: ',id,pentry(id)%top
     +			,pentry(id)%max,pentry(id)%type
	end subroutine stack_info

!===============================================================
	end module
!===============================================================

	subroutine stack_test

	use stack

	implicit none

	integer, parameter :: ndim = 100
	integer, parameter :: nloop = 10000
	integer, allocatable :: vals(:)
	integer val,value,nl,n,i,id,ind,nop
	logical bdebug
	real r

	bdebug = .true.
	bdebug = .false.

	call stack_init(id)
	allocate(vals(ndim))

	call random_seed
	val = 0
	ind = 0
	nop = 0

	do nl=1,nloop
	  call stack_rand_int(1,10,n)
	  if( bdebug ) write(6,*) 'push values: ',n
	  do i = 1,n
	    val = val + 1
	    !write(6,*) 'push: ',val
	    call stack_push(id,val)
	    ind = ind + 1
	    call stack_assert(ind <= ndim,'push',id)
	    vals(ind) = val
	  end do
	  nop = nop + n
	  call stack_rand_int(1,15,n)
	  if( bdebug ) write(6,*) 'pop values: ',n
	  do i = 1,n
	    if( stack_pop(id,value) ) then
	      !write(6,*) 'pop: ',value
	      call stack_assert(ind > 0,'pop empty',id)
	      call stack_assert(value==vals(ind),'pop value',id)
	      ind = ind - 1
	      nop = nop + 1
	    else
	      if( bdebug ) write(6,*) 'nothing to pop'
	      call stack_assert(ind == 0,'pop not empty',id)
	      exit
	    end if
	  end do
	end do

	call stack_delete(id)

	write(6,*) 'stack test successfully finished: ',nloop,nop,val

	end

!******************************************************************

        subroutine stack_assert(bcheck,text,id)
        use stack
        implicit none
        logical bcheck
        character*(*) text
        integer id
        if( .not. bcheck ) then
          write(6,*) 'stack_assertion: ',trim(text)
          call stack_info(id)
          stop 'assertion failed'
        end if
        end

	subroutine stack_rand_int(min,max,irand)

	implicit none
	integer min,max
	integer irand
	real r

	call random_number(r)
	irand = min + (max-min+1)*r

	end

!******************************************************************

	programme stack_main
	call stack_test
	end programme stack_main

!******************************************************************
