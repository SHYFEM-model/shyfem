
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010-2011,2014-2015,2019  Georg Umgiesser
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
! 23.03.2010	ggu	changed v6.1.1
! 14.04.2011	ggu	changed VERS_6_1_22
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 23.09.2015	ggu	changed VERS_7_2_4
! 16.11.2015	ggu	changed VERS_7_3_14
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

c*******************************************************************

	subroutine init_diff_oil

	use mod_lagrange
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        include 'param.h'


	integer ie

	do ie=1,nel
	  rwhvar(ie) = rwhpar
	end do

	end

c*******************************************************************

	subroutine set_diff_oil

	use mod_lagrange
	use evgeom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'femtime.h'

	real oil(nel)
	real zfilm(nel)
	real v1v(nkn)
	real v2v(nkn)

	integer ie,i
	integer ioil,ifreq
	real qoil,rhooil,ztresh,fact
	real area,roil,zoil,dz,rfact
	real rwhnew

	integer icall
	save icall
	data icall /0/

	ioil = 1
	ifreq = 6	!frequency of output

	qoil = 10.	!quantity of oil for one particle
	rhooil = 1000.
	ztresh = 0.001	!in meters
	fact = 1.

	if( ioil .le. 0 ) return

	icall = icall + 1

	do ie=1,nel
	  oil(ie) = 0.
	end do

        do i=1,nbdy
          ie=ie_body(i)
	  oil(ie) = oil(ie) + 1.
        end do

	do ie=1,nel
	  area = 12. * ev(10,ie)
	  roil = qoil * oil(ie)
	  zoil = roil / ( rhooil * area )

	  zfilm(ie) = zoil

	  rwhnew = rwhpar
	  if( zoil .gt. ztresh ) then
	    dz = zoil - ztresh
	    rfact = dz / ztresh
	    rwhnew = ( 1. + rfact*fact ) * rwhpar
	  end if
	  rwhvar(ie) = rwhnew

	  if( zoil .gt. 0. ) then
	    write(6,*) ie,zoil,ztresh,rwhnew
	  end if

	end do

        if( ifreq .gt. 0 .and. mod(icall,ifreq) .eq. 0 ) then
	  call e2n2d(zfilm,v1v,v2v)
          call wrnos2d_index(it,icall,'film','film thickness',v1v)
        end if

	end

c*******************************************************************

