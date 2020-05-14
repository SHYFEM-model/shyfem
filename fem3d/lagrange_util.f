
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009-2010,2014-2015,2017-2019  Georg Umgiesser
!    Copyright (C) 2018  Christian Ferrarin
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

c general utilities for lagrangian model
c
c revision log :
c
c 05.02.2009	ggu	copied from other files
c 23.03.2010	ggu	changed v6.1.1
c 28.03.2014	ggu	compute total length of open boundary
c 05.05.2014	ggu	changed VERS_6_1_74
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.05.2015	ggu	changed VERS_7_1_10
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 31.03.2017	ggu	changed VERS_7_5_24
c 23.08.2018	ccf	include particle_on_side routine
c 25.10.2018	ggu	changed VERS_7_5_51
c 16.02.2019	ggu	changed VERS_7_5_60
c
c*******************************************************************

        function dist_node(k1,k2)

c returns distance between two nodes

	use basin

        implicit none

        real dist_node
        integer k1,k2

        real x1,y1,x2,y2,dx,dy

        x1 = xgv(k1)
        y1 = ygv(k1)
        x2 = xgv(k2)
        y2 = ygv(k2)

        dx = x1 - x2
        dy = y1 - y2

        dist_node = sqrt( dx*dx + dy*dy )

        end

c*******************************************************************

	subroutine dist_total(ibc,totdist)

c returns total length of open boundary

	implicit none

	integer ibc
	real totdist

	integer nk,i,k1,k2
	real dxy

	integer nkbnds,kbnds
	real dist_node

	nk = nkbnds(ibc)

        totdist = 0.
        do i=2,nk
          k1 = kbnds(ibc,i-1)
          k2 = kbnds(ibc,i)
          dxy = dist_node(k1,k2)
          totdist = totdist + dxy
        end do

	end

c*******************************************************************

	subroutine basin_center(xm,ym)

c returns center of gravity of total basin

	use basin

	implicit none

	real xm,ym

	call xy_center(nkn,xgv,ygv,xm,ym)

	end
	
c*******************************************************************

	subroutine compute_total_area(area)

	use evgeom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real area

	integer ie
	double precision a,tot_area

	tot_area = 0.

	do ie=1,nel
	  a = ev(10,ie)
	  tot_area = tot_area + a
	end do

	area = 12. * tot_area

	end

!******************************************************************
! Intersect coordinate with sides of element. Only called when
! find_elem_from_old return ie=0. Return new coordinates 
! on the element side and if is an open boundary set ie=-ie

        subroutine particle_on_side(ie,xo,yo,xn,yn)

        use basin
        use mod_geom

        implicit none

        integer, intent(inout)  :: ie           !element number
        real, intent(in)        :: xo,yo        !original coordinates
        real, intent(inout)     :: xn,yn        !new coordinates

        integer                 :: ii,i1,i2
        real                    :: x1,y1,x2,y2
        real                    :: xi,yi

        integer                 :: iint
        integer                 :: segsegint    !intersection function

        do ii=1,3
          i1 = mod(ii,3) + 1
          i2 = mod(i1,3) + 1
          x1 = xgv(i1)
          y1 = ygv(i1)
          x2 = xgv(i2)
          y2 = ygv(i2)

          iint = segsegint(xo,yo,xn,yn,x1,y1,x2,y2,xi,yi)

          if (iint > 0 ) then
            if( ieltv(ii,ie) == 0 ) then        !material boundary
              xn = xi
              yn = yi
            else if ( ieltv(ii,ie) < 0 ) then   !open boundary
              ie = -ie                          !particle exits
              xn = xi
              yn = yi
            else
              !crossing more than 1 element. Stay on the side
              xn = xi
              yn = yi
            end if
            return
          end if

        end do

        end subroutine particle_on_side

!******************************************************************
