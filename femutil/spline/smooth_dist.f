
!--------------------------------------------------------------------------
!
!    Copyright (C) 2002,2004-2005,2011-2012  Georg Umgiesser
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

!**********************************************************
!
! distance computing routines
!
! revision log :
! 
! 30.01.2022    ggu     newly created from gridr.f
!
!**********************************************************

	subroutine distxy(n,xt,yt,ht,dxy,bperiod)

! sets up distance between nodes for smoothing

! dxy(i) -> distance to next node (i,i+1)

	implicit none

	integer n
	real xt(n)
	real yt(n)
	real ht(n)
	real dxy(-n:2*n)
	logical bperiod

	integer i
	real xo,yo,xn,yn
	real dx,dy,dist
	real ho,hn,h

	xn = xt(n)
	yn = yt(n)
	hn = ht(n)

	do i=1,n
	  xo = xn
	  yo = yn
	  ho = hn
	  xn = xt(i)
	  yn = yt(i)
	  hn = ht(i)

	  dx = xn - xo
	  dy = yn - yo
	  dist = sqrt( dx*dx + dy*dy )
	  if( dist == 0. ) goto 99
! if hn is set, it gives resolution factor, so give less weight to neib points
	  if( hn > 0. .and. ho > 0 ) then
	    dist = dist * (hn+ho)*0.5
	  else if( hn > 0. ) then
	    dist = dist * hn
	  else if( ho > 0. ) then
	    dist = dist * ho
	  end if

	  dxy(i-1) = dist
	  dxy(i-1+n) = dist
	  dxy(i-1-n) = dist
	end do

	dxy(2*n) = dxy(0)

	return
   99	continue
	write(6,*) 'zero distance found: ',i,i-1
	stop 'error stop distxy: zero distance'
	end

!********************************************************

	subroutine make_dist(n,xt,yt,kt,dist,bperiod)

! sets up distance on node - dist is average distance between adjacent points

	implicit none

	integer n
	real xt(n),yt(n)
	integer kt(n)
	real dist(n)
	logical bperiod

	integer k,k1,k2,id
	real dd,daver,dsigma,dmin,dmax

	id = 0
	do k=1,n
	  k1 = k-1
	  if( k == 1 ) then
	    if( bperiod ) then
	      k1 = n
	    else
	      k1 = k+1		!fake node
	    end if
	  end if
	  k2 = k+1
	  if( k == n ) then
	    if( bperiod ) then
	      k2 = 1
	    else
	      k2 = k-1		!fake node
	    end if
	  end if
	  dd = 0.
	  dd = dd + sqrt( (xt(k1)-xt(k))**2 + (yt(k1)-yt(k))**2 )
	  dd = dd + sqrt( (xt(k2)-xt(k))**2 + (yt(k2)-yt(k))**2 )
	  dist(k) = dd * 0.5
	  !write(6,*) k1,k,k2,dd,xt(k),yt(k)
	  if( dd == 0. ) then
	    id = id + 1
	    write(6,*) 'zero distance for node ',kt(k),k,k1,k2
	  end if
	end do

	if( id > 0 ) then
	  write(6,*) 'zero distances found: ',id,n
	  stop 'error stop make_dist: zero distance'
	end if

	daver = 0.
	do k=1,n
	  daver = daver + dist(k)
	end do
	daver = daver / n

	dsigma = 0.
	do k=1,n
	  dsigma = dsigma + (dist(k)-daver)**2
	end do
	dsigma = sqrt(dsigma/n)

	dmin = minval(dist)
	dmax = maxval(dist)

	write(6,*) 'dist stats: ',dmin,daver,dmax,dsigma

	end

!********************************************************

