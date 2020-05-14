
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

! revision log :
!
! 13.03.2019	ggu	changed VERS_7_5_61

	integer ndim
	parameter(ndim=10000)

	logical bperiod
	real t(ndim)
	real v(ndim)

	bperiod = .false.
	sigma = 1.
	sigma = 0.25
	i = 0

    1	continue
	  read(5,*,end=2) taux,vaux
	  i = i + 1
	  if( i .gt. ndim ) stop 'dimension...'
	  t(i) = taux
	  v(i) = vaux
	goto 1
    2	continue
	n = i

	call gsmooth(n,t,v,sigma,bperiod)

	do i=1,n
	  write(6,*) t(i),v(i)
	end do

	end

c********************************************************

