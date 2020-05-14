
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

! revision log :
!
! 19.04.2018	ggu	changed VERS_7_5_45
! 14.02.2019	ggu	changed VERS_7_5_56

	program nudge

	implicit none

	integer niter
	real c,co,cn,cobs
	real tau,dt,alfa
	real eps,high

	eps = 1.e-4
	high = 1.e+10

	niter = 0
	c = 100
	cobs = 50

	tau = 100
	dt = 750
	alfa = dt/tau

	do
	  niter = niter + 1
	  if( abs(c-cobs) < eps ) exit
	  if( abs(c) > high ) exit
	  !if( abs(cn-co) < eps ) exit
	  !cn = alfa*cobs + (1.-alfa)*c		!explicit
	  cn = (c + alfa*cobs)/(1.+alfa)	!implicit
	  write(6,*) niter,c,cn,cobs
	  c = cn
	end do

	write(6,*) 'program finished'

	end
