
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2013-2015,2019  Georg Umgiesser
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

c include file if basin (BAS) is read

        !integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        !common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        !real grav,fcor,dcor,dirn,rowass,roluft
        !common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 25.01.2013	ggu	changed VERS_6_1_62
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_2
! 01.04.2015	ggu	changed VERS_7_1_7
! 17.07.2015	ggu	changed VERS_7_1_53
! 16.02.2019	ggu	changed VERS_7_5_60

        integer nkn,nel,ngr,mbw
        common /nbasin/ nkn,nel,ngr,mbw
        save /nbasin/

        integer nkndi,neldi,ngrdi,mbwdi
	parameter ( nkndi = nkndim )
	parameter ( neldi = neldim )
	parameter ( ngrdi = ngrdim )
	parameter ( mbwdi = mbwdim )

        real dcorbas,dirnbas
        common /bkonst/ dcorbas,dirnbas
	save /bkonst/

        character*80 descrr
        common /descrr/ descrr
	save /descrr/

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v
	save /xgv/,/ygv/,/hm3v/

        integer nen3v(3,neldim)
        common /nen3v/nen3v
        integer ipev(neldim), ipv(nkndim)
        common /ipev/ipev, /ipv/ipv
        integer iarv(neldim)
        common /iarv/iarv
        integer iarnv(nkndim)
        common /iarnv/iarnv
	save /nen3v/,/ipev/,/ipv/,/iarv/,/iarnv/

