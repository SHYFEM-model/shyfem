
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

c resolution routines
c
c contents :
c
c revision log :
c 
c 30.01.2022    ggu     newly written
c
c********************************************************

	module resol

	implicit none

	logical, save :: busebck = .false.
	integer, save :: nkbg,nebg
	real, save, allocatable :: xbg(:),ybg(:),zbg(:)
	integer, save, allocatable :: eindex(:,:)

	end module resol

c********************************************************

	subroutine init_resolution(iabg,nkn,xgv,ygv,nel,nen3v,iarv,hm3v)

	use resol

	implicit none

	integer iabg
	integer nkn
	real xgv(nkn)
	real ygv(nkn)
	integer nel
	integer nen3v(3,nel)
	integer iarv(nel)
	real hm3v(3,nel)

	integer ie,ii,k,ia,ik,iie
	integer, allocatable :: nuse(:)
	real, allocatable :: huse(:)

	allocate(nuse(nkn),huse(nkn))

	nuse = 0
	huse = 0.

	nebg = 0
	nkbg = 0

	do ie=1,nel
	  ia = iarv(ie)
	  if( ia /= iabg ) cycle	!not element with background info
	  nebg = nebg + 1
	  do ii=1,3
	    k = nen3v(ii,ie)
	    nuse(k) = nuse(k) + 1
	    huse(k) = huse(k) + hm3v(ii,ie)
	  end do
	end do

	where( nuse > 0 ) huse = huse / nuse

	nkbg = count( nuse > 0 )

	write(6,*) 'total nodes: ',nkn
	write(6,*) 'total elems: ',nel
	write(6,*) 'background nodes: ',nkbg
	write(6,*) 'background elems: ',nebg

	allocate(xbg(nkbg),ybg(nkbg),zbg(nkbg))
	allocate(eindex(3,nebg))

	ik = 0
	do k=1,nkn
	  if( nuse(k) > 0 ) then
	    ik = ik + 1
	    nuse(k) = ik
	    xbg(ik) = xgv(k)
	    ybg(ik) = ygv(k)
	    zbg(ik) = huse(k)
	    !write(6,*) ik,k,xbg(ik),ybg(ik),zbg(ik)
	  end if
	end do

	iie = 0
	do ie=1,nel
	  ia = iarv(ie)
	  if( ia /= iabg ) cycle	!not element with background info
	  iie = iie + 1
	  do ii=1,3
	    k = nen3v(ii,ie)
	    eindex(ii,iie) = nuse(k)
	  end do
	  !write(6,*) iie,eindex(:,iie)
	end do

	if( iie /= nebg ) stop 'error stop: internal error'

	busebck = nebg > 0

	end

c********************************************************

	subroutine setup_resolution(inodes,area,sigma,reduct)

! computes default resolution if inodes is given, else uses reduct/sigma

	implicit none

	integer inodes		!number of desired points in outer polygon
	real area		!area of outer polygon
	real sigma		!standard deviation for smoothing
	real reduct		!minimum distance between points

	real areascale,lengthscale

! OpResolution = (4./sqrt(27)) * area/OpIntern;

	if( inodes > 0 ) then
	  areascale = (4./sqrt(27.)) * area / inodes	!this is area
	  lengthscale = sqrt(4.*areascale/sqrt(3.))	!this is length
	else
	  areascale = 1.
	  lengthscale = 1.
	end if

	if( reduct > 0. ) reduct = lengthscale / reduct
	if( sigma > 0. ) sigma = lengthscale / sigma

	write(6,*) 'using lengthscale = ',lengthscale
	write(6,*) 'using sigma       = ',sigma
	write(6,*) 'using reduct      = ',reduct

	end

c********************************************************

	subroutine compute_weights(n,xt,yt,ht,wt)

! computes weights of points on line
!
! all points will have a weight or -1

	implicit none

	integer n
	real xt(n)
	real yt(n)
	real ht(n)
	real wt(n)

	integer i
	real rval

	wt = 0.

	do i=1,n
	  call get_weight_from_backg(xt(i),yt(i),rval)
	  if( rval == 0. ) rval = ht(i)
	  if( ht(i) < 0. ) rval = -1.
	  wt(i) = rval
	end do

	end

c********************************************************

	subroutine get_weight_from_backg(x0,y0,rval)

	use resol

	implicit none

	real x0,y0
	real rval

	integer ie,ii,k
	real xt(3),yt(3),zt(3)

	logical inconvex
	real rintrz

	do ie=1,nebg
	  do ii=1,3
	    k = eindex(ii,ie)
	    xt(ii) = xbg(k)
	    yt(ii) = ybg(k)
	    zt(ii) = zbg(k)
	  end do
	  if( inconvex(3,xt,yt,x0,y0) ) then
	    rval = rintrz(xt,yt,x0,y0,zt)
	    return
	  end if
	end do

	rval = 0.

	end

c********************************************************

