
!--------------------------------------------------------------------------
!
!    Copyright (C) 2019  Georg Umgiesser
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

! routines for fifo pipes
!
! revision log :
!
! 04.07.2019    ggu     written from scratch

!===========================================================
	module fifo_shyfem
!===========================================================

	implicit none

	character*80, save :: fifo_myname = ' '

	character*80, save :: fifo_w2s = 'fifo_w2s.fifo'
	character*80, save :: fifo_s2w = 'fifo_s2w.fifo'

	integer, save :: iu_w2s = 601
	integer, save :: iu_s2w = 602

	integer, save :: iu_send = 0
	integer, save :: iu_rec = 0

	character*80, save :: a_w2s = ' '
	character*80, save :: a_s2w = ' '

	logical, save :: bshyfem = .false.

	logical, save :: bdebug = .false.

        INTERFACE send_array
        MODULE PROCEDURE send_array_i,send_array_r,send_array_d
        END INTERFACE

        INTERFACE rec_array
        MODULE PROCEDURE rec_array_i,rec_array_r,rec_array_d
        END INTERFACE

!===========================================================
	contains
!===========================================================

!***********************************************************

	subroutine fifo_setname(name)

	implicit none

	character*(*) name

	fifo_myname = name

	end subroutine

!***********************************************************

	subroutine fifo_setup

	implicit none

	write(6,*) 'executable running: ',trim(fifo_myname)

	if( fifo_myname == 'shyfem' ) then
	  write(6,*) 'fifo_setup: ','shyfem'
	  iu_send = iu_s2w
	  iu_rec = iu_w2s
	  a_s2w = 'write'
	  a_w2s = 'read'
	  bshyfem = .true.
	else if( fifo_myname /= ' ' ) then
	  write(6,*) 'fifo_setup: ','other'
	  iu_w2s = iu_w2s + 10
	  iu_s2w = iu_s2w + 10
	  iu_send = iu_w2s
	  iu_rec = iu_s2w
	  a_s2w = 'read'
	  a_w2s = 'write'
	else
	  write(6,*) 'fifo_myname is not set'
	  stop 'error stop fifo_setup: no name'
	end if

	end subroutine

!***********************************************************

	subroutine fifo_open

	implicit none

	character*80 name,action

        write(6,*) 'fifo_open: ',trim(fifo_myname),iu_w2s,iu_s2w

	!open(iu_w2s,file=fifo_w2s,form='unformatted',action=a_w2s)
	!open(iu_s2w,file=fifo_s2w,form='unformatted',action=a_s2w)
	open(iu_w2s,file=fifo_w2s,access='stream',action=a_w2s)
	open(iu_s2w,file=fifo_s2w,access='stream',action=a_s2w)

        inquire(iu_w2s,name=name,action=action)
        write(6,*) 'fifo_open: ',trim(fifo_myname),iu_w2s                &
     &			,trim(name),'  ',trim(action)

        inquire(iu_s2w,name=name,action=action)
        write(6,*) 'fifo_open: ',trim(fifo_myname),iu_s2w                &
     &			,trim(name),'  ',trim(action)

	end subroutine

!***********************************************************

	subroutine fifo_close

	implicit none

	close(iu_send)
	close(iu_rec)

	end subroutine

!***********************************************************
!***********************************************************
!***********************************************************

	subroutine rec_array_i(n,array,ier)

	implicit none

	integer n
	integer array(n)
	integer ier

	integer n1,n2

        character*80 name,action
        inquire(iu_rec,name=name,action=action)
        write(6,*) 'receiving: ',iu_rec,trim(name),'  ',trim(action)

	read(iu_rec,iostat=ier) n1,array(1:n1),n2
	write(6,*) 'gggggggggggguuuuu: ',ier
	if( ier /= 0 ) return

	if( n1 /= n2 ) then
	  write(6,*) n,n1,n2
	  stop 'error stop rec_array_i: mismatch n1/=n2'
	else if( n1 /= n ) then
	  write(6,*) n,n1,n2
	  stop 'error stop rec_array_i: not of expected size n1/=n'
	end if

	end subroutine

!***********************************************************

	subroutine rec_array_r(n,array,ier)

	implicit none

	integer n
	real array(n)
	integer ier

	integer n1,n2

	read(iu_rec,iostat=ier) n1,array(1:n1),n2
	if( ier /= 0 ) return

	if( n1 /= n2 ) then
	  write(6,*) n,n1,n2
	  stop 'error stop rec_array_r: mismatch n1/=n2'
	else if( n1 /= n ) then
	  write(6,*) n,n1,n2
	  stop 'error stop rec_array_r: not of expected size n1/=n'
	end if

	end subroutine

!***********************************************************

	subroutine rec_array_d(n,array,ier)

	implicit none

	integer n
	double precision array(n)
	integer ier

	integer n1,n2

	read(iu_rec,iostat=ier) n1,array(1:n1),n2
	if( ier /= 0 ) return

	if( n1 /= n2 ) then
	  write(6,*) n,n1,n2
	  stop 'error stop rec_array_d: mismatch n1/=n2'
	else if( n1 /= n ) then
	  write(6,*) n,n1,n2
	  stop 'error stop rec_array_d: not of expected size n1/=n'
	end if

	end subroutine

!***********************************************************

	subroutine send_array_i(n,array)

	implicit none

	integer n
	integer array(n)

        character*80 name,action
        inquire(iu_send,name=name,action=action)
        write(6,*) 'sending: ',iu_send,trim(name),'  ',trim(action),n

	write(iu_send) n,array(1:n),n
	flush(iu_send)

	end subroutine

!***********************************************************

	subroutine send_array_r(n,array)

	implicit none

	integer n
	real array(n)

	write(iu_send) n,array(1:n),n
	flush(iu_send)

	end subroutine

!***********************************************************

	subroutine send_array_d(n,array)

	implicit none

	integer n
	double precision array(n)

	write(iu_send) n,array(1:n),n
	flush(iu_send)

	end subroutine

!===========================================================
	end module fifo_shyfem
!===========================================================

