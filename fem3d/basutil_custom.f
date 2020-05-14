
!--------------------------------------------------------------------------
!
!    Copyright (C) 2020  Georg Umgiesser
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
! 01.04.2020	ggu	written from scratch

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine bas_custom

	!call sort_basin
	call black_sea_nudge

	end

c***************************************************************

        subroutine black_sea_nudge

        use basin

        implicit none

        integer k,l
        real xc1,yc1,xc2,yc2
        real xm1,ym1,xm2,ym2
        real tc,tm
        real scalc,scalm
        real ncx,ncy,nmx,nmy
        real x,y
	real dreg
        real tau
	character*80 file,string

        real tauv(nkn)

        xc1 = 29.
        yc1 = 43.5
        xc2 = 30.5
        yc2 = 45.5

        xm1 = 29.5
        ym1 = 43.5
        xm2 = 31.
        ym2 = 45.5

        tc = 259200.
        tm = 86400.

        ncx = -(yc2-yc1)
        ncy = +(xc2-xc1)
        nmx = -(ym2-ym1)
        nmy = +(xm2-xm1)

        do k=1,nkn
          x = xgv(k)
          y = ygv(k)

          scalc = ncx * (x-xc1) + ncy * (y-yc1)
          scalm = nmx * (x-xm1) + nmy * (y-ym1)

          if( scalc .gt. 0. ) then              !coast
            tau = tc
          else if( scalm .lt. 0. ) then         !open sea
            tau = tm
          else                                  !intermediate
            tau = 0.5 * (tc+tm)
          end if

          tauv(k) = tau
        end do

	file = 'taunudge.fem'
	string = 'relaxation time scale'
	dreg = 0.2
	call write_regular_2d_nodes(file,string,dreg,tauv)

        write(6,*) 'relaxation time for nudging set up'

        end

c***************************************************************

