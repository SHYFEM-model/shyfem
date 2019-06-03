
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
! 27.01.2011	ggu	changed VERS_6_1_17
! 14.02.2012	ggu	changed VERS_6_1_44
! 16.03.2012	ggu	changed VERS_6_1_48
! 12.09.2013	ggu	changed VERS_6_1_67
! 12.12.2014	ggu	changed VERS_7_0_9
! 19.01.2015	ggu	changed VERS_7_1_3
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_52
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 30.07.2015	ggu	changed VERS_7_1_83
! 12.10.2015	ggu	changed VERS_7_3_3
! 25.05.2016	ggu	changed VERS_7_5_10
! 14.11.2017	ggu	changed VERS_7_5_36
! 22.02.2018	ggu	changed VERS_7_5_42
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

!***************************************************************

	subroutine grd_to_basin

	use grd
	use basin

	implicit none

	integer nk,ne,nl,nne,nnl
	integer k,ie,ii,ibase,n,i
	integer nhn,nhe
	real flag
	logical bdebug,bwrite

!--------------------------------------------------------
! initialize
!--------------------------------------------------------

	flag = -999.
	bdebug = .false.
	bwrite = .false.

	call grd_get_params(nk,ne,nl,nne,nnl)
	call basin_init(nk,ne)

!--------------------------------------------------------
! copy params
!--------------------------------------------------------

	dcorbas = dcor_grd
	dirnbas = dirn_grd
	descrr = title_grd

!--------------------------------------------------------
! copy arrays
!--------------------------------------------------------

	xgv = xv
	ygv = yv
	ipev = ippev
	ipv = ippnv
	iarnv = ianv
	iarv = iaev

!--------------------------------------------------------
! copy element index
!--------------------------------------------------------

	do ie=1,ne
	  ibase = ipntev(ie-1)
	  n = ipntev(ie) - ipntev(ie-1)
	  if( n .ne. 3 ) then
	    !write(6,*) (ipntev(i),i=0,10)
	    write(6,*) 'element is not triangle: ',ie,ippev(ie)
	    stop 'error stop grd_to_basin: not a triangle'
	  end if
	  do ii=1,3
	    nen3v(ii,ie) = inodev(ibase+ii)
	  end do
	end do

!--------------------------------------------------------
! check depth and copy
!--------------------------------------------------------

	nhe = 0
	do ie=1,ne
	  if( hhev(ie) .ne. flag ) nhe = nhe + 1
	end do

	nhn = 0
	do k=1,nk
	  if( hhnv(k) .ne. flag ) nhn = nhn + 1
	end do

	if( nhe > 0 .and. nhn > 0 ) then
	  write(6,*) 'nhe,nhn: ',nhe,nhn
	  if( nhe == nel .and. nhn ==nkn ) then
	    write(6,*) 'can use both depths...'
	    write(6,*) '... using element values'
	    nhn = 0
	  else if( nhe == nel ) then
	    write(6,*) 'depth on nodes incomplete...'
	    write(6,*) '... using element values'
	    nhn = 0
	  else if( nhn == nkn ) then
	    write(6,*) 'depth on elements incomplete...'
	    write(6,*) '... using nodal values'
	    nhe = 0
	  else
	    write(6,*) 'depth values are incomplete...'
	    stop 'error stop grd_to_basin: depth incomplete'
	  end if
	end if

	if( bwrite ) write(6,*) 'nhe,nhn: ',nhe,nhn

	if( nhn == 0 ) then
	  do ie=1,ne
	    do ii=1,3
	      hm3v(ii,ie) = hhev(ie)
	    end do
	  end do
	else
	  do ie=1,ne
	    do ii=1,3
	      k = nen3v(ii,ie)
	      hm3v(ii,ie) = hhnv(k)
	    end do
	  end do
	end if

!--------------------------------------------------------
! debug output
!--------------------------------------------------------

	if( bdebug ) then

	do ie=1,ne,ne/10
	  write(6,*) ie,(nen3v(ii,ie),ii=1,3)
	end do

	do k=1,nk,nk/10
	  write(6,*) k,ippnv(k),xv(k),yv(k)
	end do

	end if

!--------------------------------------------------------
! end of routine
!--------------------------------------------------------

	end

!***************************************************************

	subroutine basin_to_grd

	use grd
	use basin

	implicit none

	integer nk,ne,nl,nne,nnl
	integer k,ie,ii,ibase,n
	integer nhn,nhe
	real flag
	logical bdebug
	logical bconst,bunique
	double precision h,he,hm

!--------------------------------------------------------
! initialize
!--------------------------------------------------------

	flag = -999.
	bdebug = .false.

	nl = 0
	nne = 3*nel
	nnl = 0

	write(6,*) 'basin_to_grd: ',nkn,nel
	call grd_init(nkn,nel,nl,nne,nnl)

!--------------------------------------------------------
! copy params
!--------------------------------------------------------

	dcor_grd = dcorbas
	dirn_grd = dirnbas
	title_grd = descrr

!--------------------------------------------------------
! copy arrays
!--------------------------------------------------------

	xv(1:nkn) = xgv(1:nkn)
	yv(1:nkn) = ygv(1:nkn)
	ippev(1:nel) = ipev(1:nel)
	ippnv(1:nkn) = ipv(1:nkn)
	ianv(1:nkn) = iarnv(1:nkn)
	iaev(1:nel) = iarv(1:nel)

!--------------------------------------------------------
! copy element index
!--------------------------------------------------------

	ibase = 0
	ipntev(0) = 0
	do ie=1,nel
	  ibase = ipntev(ie-1)
	  do ii=1,3
	    ibase = ibase + 1
	    inodev(ibase) = nen3v(ii,ie)
	  end do
	  ipntev(ie) = ibase
	end do

!--------------------------------------------------------
! check depth and copy
!--------------------------------------------------------

	hhnv = flag
	hhev = flag
	h = flag

	bunique = .true.
	bconst = .true.

	do ie=1,nel
	  he = hm3v(1,ie)
	  hm = 0.
	  do ii=1,3
	    h = hm3v(ii,ie)
	    k = nen3v(ii,ie)
	    if( hhnv(k) == flag ) hhnv(k) = h
	    if( hhnv(k) /= h ) bunique = .false.
	    if( he /= h ) bconst = .false.
	    hm = hm + h
	  end do
	  hm = hm / 3.
	  hhev(ie) = hm
	end do

	if( bconst .and. bunique ) then		!constant depth
	  hhev = h
	  hhnv = flag
	else if( bunique ) then			!unique depth at nodes
	  hhev = flag
	else
	  hhnv = flag
	end if
	  
!--------------------------------------------------------
! end of routine
!--------------------------------------------------------

	end

!***************************************************************

	subroutine grd_get_nodal_depth(hk)

	use grd

	implicit none

	real hk(nk_grd)

	hk = hhnv

	end

!***************************************************************

	subroutine grd_set_unused_node_depth(hk)

	use basin
	use grd

	implicit none

	real hk(nkn)

	integer ie,ii,k
	integer iuse(nkn)

	iuse = 0

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    iuse(k) = iuse(k) + 1
	  end do
	end do

	do k=1,nkn
	  if( iuse(k) == 0 ) hhnv(k) = hk(k)
	end do

	end

!***************************************************************

	subroutine grd_set_coords(n,x,y)

	use grd

	implicit none

	integer n
	real x(n)
	real y(n)

	if( nk_grd /= n ) then
	  write(6,*) nk_grd,n
	  stop 'error stop grd_set_coords: nkn id different'
	end if

	xv = x
	yv = y

	end

!***************************************************************

