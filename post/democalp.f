
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

	program democalp

	call plots(0,0,0)

	call plot(1.,1.,3)
	call plot(9.,9.,-2)

	call factor(0.5)

	call plot(1.,9.,2)
	call symbol( 999.,999.,1.,'text',0.,4)

	call factor(1.0)

	call symbol( 5.,4.,1.,'hans-georg',0.,10)

	call shade(3.,3.,1.,1.,9)

	call plots(0,0,0)

	call shade(2.,2.,5.,1.,3)
	call symbol(2.,2.,1.,'GEORG',0.,5)

	call plot(2.,4.,-3)
	call factor(0.5)
	call shade(2.,2.,5.,1.,3)
	call symbol(2.,2.,1.,'GEORG',0.,5)
	call shade(9.,2.,10.,2.,3)
	call symbol(9.,2.,1./0.5,'GEORG',0.,5)

	call plots(0,0,0)

	h=0.1
	dh=0.01
	x=0.5
	y=4.
	do i=1,100
	h=h+dh
	call shade(x,y,h,h,3)
	call symbol(x,y,h,'R',0.,1)
	x=x+h
	end do

	call plot(0.,0.,999)

	end

