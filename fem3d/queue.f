
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

! queue utility routines
!
! still to do: implement queue for string
!
! revision log :
!
! 26.09.2017	ggu	changed VERS_7_5_32
! 24.01.2018	ggu	changed VERS_7_5_41
! 14.02.2019	ggu	changed VERS_7_5_56

!===============================================================
	module queue
!===============================================================

	implicit none

	private

        type :: entry
          integer :: front
          integer :: rear
          integer :: fill
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

	public :: queue_init		!call queue_init(id)
	public :: queue_delete		!call queue_delete(id)
	public :: queue_enqueue		!call queue_enqueue(id,value)
	public :: queue_dequeue		!logical queue_dequeue(id,value)
	public :: queue_peek		!logical queue_peek(id,value)
	public :: queue_is_empty	!logical queue_is_empty(id)
	public :: queue_filling 	!integer queue_filling(id)
	public :: queue_info		!call queue_info(id)

        INTERFACE queue_enqueue
        MODULE PROCEDURE         queue_enqueue_d
     +                          ,queue_enqueue_r
     +                          ,queue_enqueue_i
        END INTERFACE

        INTERFACE queue_dequeue
        MODULE PROCEDURE         queue_dequeue_d
     +                          ,queue_dequeue_r
     +                          ,queue_dequeue_i
        END INTERFACE

        INTERFACE queue_peek
        MODULE PROCEDURE         queue_peek_d
     +                          ,queue_peek_r
     +                          ,queue_peek_i
        END INTERFACE

!===============================================================
	contains
!===============================================================

        subroutine queue_init_alloc

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

        end subroutine queue_init_alloc

!******************************************************************

        subroutine queue_init_new_id(id)

        integer id

        idlast = idlast + 1
        if( idlast > ndim ) then
          call queue_init_alloc
        end if
        id = idlast

        call queue_init_id(id)

        end subroutine queue_init_new_id

!******************************************************************

        subroutine queue_init_id(id)

        integer id

        if( id > ndim ) then
          stop 'error stop queue_init_id: ndim'
        end if

        pentry(id)%front = 1
        pentry(id)%rear = 0
        pentry(id)%fill = 0
        pentry(id)%max = 0
        pentry(id)%type = 0

	if( allocated(pentry(id)%array) ) deallocate(pentry(id)%array)
	if( allocated(pentry(id)%string) ) deallocate(pentry(id)%string)

        end subroutine queue_init_id

!******************************************************************

	subroutine queue_error(id,error)

	integer id,error

	if( error == empty_error ) then
	  write(6,*) 'queue: ',id
	  stop 'error stop queue: queue is empty'
	else if( error == type_error ) then
	  write(6,*) 'queue: ',id
	  write(6,*) 'type: ',pentry(id)%type
	  stop 'error stop queue: variable is of wrong type'
	else
	  stop 'error stop queue: internal error (1)'
	end if

	end subroutine queue_error

!******************************************************************

	subroutine queue_rotate(id)

! rotates queue so that front is on first position
! after this we can make array bigger or smaller (depending on rear)

	integer id

	integer iact,max,fill,i
	double precision daux(pentry(id)%max)

        iact = pentry(id)%front
        fill = pentry(id)%fill
        max = pentry(id)%max

	if( fill == 0 ) return

	do i=1,fill
	  daux(i) = pentry(id)%array(iact)
	  iact = mod(iact,max) + 1
	end do

        pentry(id)%front = 1
        pentry(id)%rear = pentry(id)%fill
	pentry(id)%array(1:fill) = daux(1:fill)

	end subroutine queue_rotate

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

	subroutine queue_init(id)
	integer id
        call queue_init_new_id(id)
	end subroutine queue_init

	subroutine queue_delete(id)
	integer id
        call queue_init_id(id)
	if( id == idlast ) idlast = idlast - 1
	end subroutine queue_delete

!--------------------

	subroutine queue_enqueue_i(id,value)
	integer id
	integer value
	call queue_enqueue_d(id,dble(value))
	end subroutine queue_enqueue_i

	subroutine queue_enqueue_r(id,value)
	integer id
	real value
	call queue_enqueue_d(id,dble(value))
	end subroutine queue_enqueue_r

	subroutine queue_enqueue_d(id,value)
	integer id
	double precision value
	integer max,rear,fill
	max = pentry(id)%max
	if( pentry(id)%fill >= max ) then
	  call queue_rotate(id)
	  max = 2 * max
	  call realloc_double(max,pentry(id)%array)
	  pentry(id)%max = max
	end if
	if( pentry(id)%type == no_type ) then
	  pentry(id)%type = value_type
	end if
	if( pentry(id)%type /= value_type ) then
	  call queue_error(id,type_error)
	end if
	rear = pentry(id)%rear
	rear = mod(rear,max) + 1
	pentry(id)%rear = rear
	pentry(id)%fill = pentry(id)%fill + 1
	pentry(id)%array(pentry(id)%rear) = value
	end subroutine queue_enqueue_d

!--------------------

	logical function queue_dequeue_i(id,value)
	integer id
	integer value
	double precision dvalue
	queue_dequeue_i = queue_dequeue_d(id,dvalue)
	value = nint(dvalue)
	end function queue_dequeue_i

	logical function queue_dequeue_r(id,value)
	integer id
	real value
	double precision dvalue
	queue_dequeue_r = queue_dequeue_d(id,dvalue)
	value = real(dvalue)
	end function queue_dequeue_r

	logical function queue_dequeue_d(id,value)
	integer id
	double precision value
	integer front,max
	queue_dequeue_d = queue_peek_d(id,value)
	if( queue_dequeue_d ) then
	  front = pentry(id)%front
	  max = pentry(id)%max
	  front = mod(front,max) + 1
	  pentry(id)%front = front
	  pentry(id)%fill = pentry(id)%fill - 1
	else
	  value = 0.
	end if
	end function queue_dequeue_d

!--------------------

	logical function queue_peek_i(id,value)
	integer id
	integer value
	double precision dvalue
	queue_peek_i = queue_peek_d(id,dvalue)
	value = nint(dvalue)
	end function queue_peek_i

	logical function queue_peek_r(id,value)
	integer id
	real value
	double precision dvalue
	queue_peek_r = queue_peek_d(id,dvalue)
	value = real(dvalue)
	end function queue_peek_r

	logical function queue_peek_d(id,value)
	integer id
	double precision value
	queue_peek_d = .false.
	if( pentry(id)%fill == 0 ) return
	if( pentry(id)%type /= value_type ) then
	  call queue_error(id,type_error)
	end if
	value = pentry(id)%array(pentry(id)%front)
	queue_peek_d = .true.
	end function queue_peek_d

!--------------------

	logical function queue_is_empty(id)
	integer id
	queue_is_empty = ( pentry(id)%fill == 0 )
	end function queue_is_empty

!--------------------

	integer function queue_filling(id)
	integer id
	queue_filling = pentry(id)%fill 
	end function queue_filling

!--------------------

	subroutine queue_info(id)
	integer id
	write(6,*) '----------------------------'
	write(6,*) 'front: ',pentry(id)%front
	write(6,*) 'rear:  ',pentry(id)%rear
	write(6,*) 'fill:  ',pentry(id)%fill
	write(6,*) 'nmax:  ',pentry(id)%max
	write(6,*) 'type:  ',pentry(id)%type
        !write(6,*) pentry(id)%array(1:pentry(id)%max)
	call queue_print_good(id)
	write(6,*) '----------------------------'
	end subroutine queue_info

	subroutine queue_print_good(id)
	integer id
	integer front,rear,max,i
	if( pentry(id)%fill == 0 ) return
	front=pentry(id)%front
	rear=pentry(id)%rear
	max=pentry(id)%max
	if( front <= rear ) then
	  write(6,*) (pentry(id)%array(i),i=front,rear)
	else
	  write(6,*) (pentry(id)%array(mod(i-1,max)+1),i=front,rear+max)
	end if
	end subroutine queue_print_good

!===============================================================
	end module
!===============================================================

	subroutine queue_test

	use queue

	implicit none

	logical bdebug
	integer nloop,val,value,nl,n,i,id,valold
	real r

	bdebug = .true.
	bdebug = .false.

	call queue_init(id)

	call random_seed
	nloop = 1000
	val = 0
	valold = 0

	do nl=1,nloop
	  call rand_int(1,12,n)
	  if( bdebug ) write(6,*) 'enqueue values: ',n
	  do i = 1,n
	    val = val + 1
	    if( bdebug ) write(6,*) 'enqueue: ',val
	    call queue_enqueue(id,val)
	  end do
	  !call queue_info(id)
	  call rand_int(1,10,n)
	  if( bdebug ) write(6,*) 'dequeue values: ',n
	  do i = 1,n
	    if( queue_dequeue(id,value) ) then
	      valold = valold + 1
	      if( bdebug ) write(6,*) 'dequeue: ',value
	      if( value /= valold ) then
	        write(6,*) 'error dequeuing: ',value,valold
	        stop 'error stop'
	      end if
	    else
	      if( bdebug ) write(6,*) 'nothing to dequeue'
	      exit
	    end if
	  end do
	end do

	call queue_delete(id)

	write(6,*) 'queue test successfully finished: ',nloop

	contains

	subroutine rand_int(min,max,irand)
	implicit none
	integer min,max
	integer irand
	real r
	call random_number(r)
	irand = min + (max-min+1)*r
	end

	end

!******************************************************************


!******************************************************************
!	program queue_main
!	call queue_test
!	end program queue_main
!******************************************************************

