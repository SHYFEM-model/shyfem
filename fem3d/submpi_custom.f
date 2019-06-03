
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

! revision log :
!
! 07.12.2017	ggu	changed VERS_7_5_40
! 24.01.2018	ggu	changed VERS_7_5_41
! 19.04.2018	ggu	changed VERS_7_5_45
! 26.04.2018	ggu	changed VERS_7_5_46
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine make_custom_domain_area(area_node)

	use basin
	use shympi

	implicit none

	integer area_node(nkn)

	integer ie,ii
	real r,h

	if( nkn /= 225 .or. nel /= 384 ) return

	if( n_threads == 1 ) then
	  return
	else if( n_threads == 2 ) then
	  call make_domain_area_2(area_node)
	else if( n_threads == 3 ) then
	  call make_domain_area_3(area_node)
	else if( n_threads == 4 ) then
	  call make_domain_area_4(area_node)
	else
	  write(6,*) 'n_threads = ',n_threads
	  stop 'error stop make_domain_area: cannot handle'
	end if

	do ie=1,nel
	  do ii=1,3
	    call random_number(r)
	    h = 10.* r
	    h = 2.* r
	    !hm3v(ii,ie) = h
	    !write(6,*) ie,ii,r,h
	  end do
	end do

	end

!*****************************************************************

	subroutine make_domain_area_2(area_node)

	use basin
	use shympi

	implicit none

	integer k
	integer area_node(nkn)

	do k=1,nkn
	  if( ygv(k) > 3100.  ) then
	    area_node(k) = 1
	  end if
	end do

	end

!*****************************************************************

	subroutine make_domain_area_3(area_node)

	use basin
	use shympi

	implicit none

	integer k
	integer area_node(nkn)

	do k=1,nkn
	  if( ygv(k) > 4100.  ) then
	    area_node(k) = 2
	  else if( ygv(k) > 2100.  ) then
	    area_node(k) = 1
	  end if
	end do

	end

!*****************************************************************

	subroutine make_domain_area_4(area_node)

	use basin
	use shympi

	implicit none

	integer k
	integer area_node(nkn)

	do k=1,nkn
	  if( ygv(k) > 3100.  ) then
	    area_node(k) = 2
	  end if
	  if( xgv(k) > 100.  ) then
	    area_node(k) = area_node(k) + 1
	  end if
	end do

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

