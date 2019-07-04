
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

! routines for coupling WW3 model
!
! revision log :
!
! 04.07.2019    ggu     written from scratch

!===========================================================
        module ww3_shyfem
!===========================================================

	implicit none

	logical, save :: bww3 = .false.

!===========================================================
        end module ww3_shyfem
!===========================================================

!***********************************************************

	subroutine ww3_init

	use fifo_shyfem
	use ww3_shyfem
	use basin

	implicit none

	integer ier
	integer iwave
	integer itest(2)
	real raux(nkn)

	real getpar

	iwave = nint(getpar('iwave'))
	bww3 = ( iwave == 11 )

	if( .not. bww3 ) return

	call fifo_setname('shyfem')
	call fifo_setup
	call fifo_open

	call rec_array(2,itest,ier)
	if( ier /= 0 ) goto 99
	write(6,*) 'reading shyfem: ',itest

	call rec_array(nkn,raux,ier)
	if( ier /= 0 ) goto 99
	write(6,*) 'x ',nkn,raux(1),raux(nkn)
	call rec_array(nkn,raux,ier)
	if( ier /= 0 ) goto 99
	write(6,*) 'y ',nkn,raux(1),raux(nkn)

	return
   99	continue
	write(6,*) 'error receiving data ',ier
	stop 'error stop ww3_init: error receiving'
	end

!***********************************************************

	subroutine ww3_loop

	use fifo_shyfem
	use ww3_shyfem
	use basin
	use mod_hydro

	implicit none

	integer it
	integer itest(2)
	real raux(nkn)
	double precision dtime

	if( .not. bww3 ) return

	call get_act_dtime(dtime)
	it = nint(dtime)
	call send_array(1,(/it/))
	call send_array(nkn,znv)

	end

!***********************************************************

