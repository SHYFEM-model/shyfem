
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

c simulates continuous release over open boundaries
c
c revision log :
c
c 12.12.2007	ggu	written from scratch
c 23.03.2010	ggu	changed v6.1.1
c 21.05.2015	ggu	changed VERS_7_1_11
c 16.11.2015	ggu	changed VERS_7_3_14
c 25.10.2018	ggu	changed VERS_7_5_51
c 16.02.2019	ggu	changed VERS_7_5_60
c
c*******************************************************************

	subroutine lgr_larvae

c manages larvae

	use mod_lagrange

	implicit none

	include 'param.h'

	integer i,ie,ii
	real rlinit,x,y,z,rl
        double precision xx,yy
        double precision xi(3)

	real getpar

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. 0 ) then
	  rlinit = 0.
	  do i=1,nbdy
	    lgr_ar(i)%custom(1) = rlinit
	  end do
	end if

	icall = icall + 1

	do i=1,nbdy
          do ii=1,3
            xi(ii) = lgr_ar(i)%actual%xi(ii)
          end do
	  ie = lgr_ar(i)%actual%ie
          z  = lgr_ar(i)%actual%z	!rel. depth   0=surface  1=bottom
	  rl = lgr_ar(i)%custom(1)
          call xi2xy(abs(ie),xx,yy,xi)
          x   = xx
          y   = yy
	  
	  if( ie .gt. 0 ) then
	    call treat_larva(x,y,z,ie,rl)
	  end if

	  lgr_ar(i)%custom(1) = rl
          lgr_ar(i)%actual%z = z
	end do

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*******************************************************************

	subroutine treat_larva(x,y,z,ie,rl)

	implicit none

	real x,y,z
	integer ie
	real rl,r,perc

	real ggrand

	perc = 0.1
	perc = 0.
	r = ggrand(0)

	if( r .lt. perc ) z = 1.0

	end

c*******************************************************************

