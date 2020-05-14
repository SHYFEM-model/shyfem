
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
	real xaux(nkn)
	real yaux(nkn)
	logical bex

	real getpar

	iwave = nint(getpar('iwave'))
	bww3 = ( iwave == 11 )

	if( .not. bww3 ) return

	do
	  exit
	  inquire(file='started.fifo',exist=bex)
	  if( bex ) exit
	  write(6,*) 'file not existing...'
	  call sleep(1)
	end do

	call fifo_setname('shyfem')
	call fifo_setup
	call fifo_open

	write(6,*) '### shyfem-ww3 ### ready for reading from ww3...'
	call rec_array(2,itest,ier)
	write(6,*) '### shyfem-ww3 ### ier = ',ier
	if( ier /= 0 ) goto 99
	write(6,*) '### shyfem-ww3 ### reading shyfem: ',itest
	flush(6)

	!call rec_array(nkn,xaux,ier)
	!if( ier /= 0 ) goto 99
	!write(6,*) 'x ',nkn,xaux(1),xaux(nkn)
	!call rec_array(nkn,yaux,ier)
	!if( ier /= 0 ) goto 99
	!write(6,*) 'y ',nkn,yaux(1),yaux(nkn)

	!call create_index(xaux,yaux)

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

	return
	if( .not. bww3 ) return

	call get_act_dtime(dtime)
	it = nint(dtime)
	call send_array(1,(/it/))
	call send_array(nkn,znv)

	end

!***********************************************************

	subroutine create_index(xaux,yaux)

! try sort on basin points - only for testing

	use basin

	implicit none

	real xaux(nkn),yaux(nkn)

	integer k
	integer ip1,ip2,ia,ju1,ju2,ju,iaa
	integer nerror,nequal
	real x,y,ra,xx,yy
	logical bwrite
	real rk(nkn)
	integer index(nkn)
	integer is2w(nkn)

	integer locater

	write(6,*) 'create_index...',nkn

	bwrite = .true.
	bwrite = .false.
	nerror = 0
	nequal = 0

!--------------------------------------------------------
! find one array for sorting
!--------------------------------------------------------

	do k=1,nkn
	  x = xaux(k)
	  y = yaux(k)
	  rk(k) = x + y
	end do

!--------------------------------------------------------
! sort index
!--------------------------------------------------------

	call isortr(nkn,rk,index)

!--------------------------------------------------------
! check if error - compute how many entries are non-unique
!--------------------------------------------------------

	do k=2,nkn
	  ip1 = index(k-1)
	  ip2 = index(k)
	  if( rk(ip1) > rk(ip2) ) then
	    write(6,*) '*** error ',k,rk(ip1),rk(ip2)
	    write(6,*) xgv(ip1),xgv(ip2),ygv(ip1),ygv(ip2)
	    nerror = nerror + 1
	  else if( rk(ip1) == rk(ip2) ) then
	    if( bwrite ) then
	      write(6,*) 'equal ',k,rk(ip1),rk(ip2)
	      write(6,*) xgv(ip1),xgv(ip2),ygv(ip1),ygv(ip2)
	    end if
	    nequal = nequal + 1
	  end if
	end do

	write(6,*) 'sorted basin tested...',nkn,nerror,nequal
	if( nerror>0 ) stop 'error stop create_index: sort coords'

!--------------------------------------------------------
! find index for copying
!--------------------------------------------------------

	nerror = 0
	nequal = 0

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  ra = x + y
	  call locater_all(nkn,rk,index,ra,ju1,ju2)
	  if( ju1 == 0 ) then
	    ia = 0
	    write(6,*) '*** no entry found: ',k,ra
	  else if( ju1 /= ju2 ) then
	    if( bwrite ) write(6,*) k,ju1,ju2,ju2-ju1+1
	    nequal = nequal + 1
	    ia = 0
	    do ju=ju1,ju2
	      iaa = index(ju)
	      if( xaux(iaa) == x .and. yaux(iaa) == y ) ia = iaa
	    end do
	  else
	    ia = index(ju1)
	  end if
	  if( ia == 0 ) then
	    write(6,*) 'cannot find entry: ',k,x,y
	    nerror = nerror + 1
	  end if
	  is2w(k) = ia
	end do

	write(6,*) 'all nodes found...',nkn,nerror,nequal

	if( nerror > 0 ) then
	  write(6,*) 'error creating index: ',nerror
	  stop 'error stop create_index: creating index'
	end if

!--------------------------------------------------------
! final check
!--------------------------------------------------------

	nerror = 0

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  xx = xaux(is2w(k))
	  yy = yaux(is2w(k))
	  if( x /= xx .or. y /= yy ) then
	    nerror = nerror + 1
	    write(6,*) k,x,xx,y,yy
	  end if
	end do

	if( nerror > 0 ) then
	  write(6,*) 'error checking index: ',nerror
	  stop 'error stop create_index: final check of index'
	end if

!--------------------------------------------------------
! end of routine
!--------------------------------------------------------

	end

!***********************************************************

