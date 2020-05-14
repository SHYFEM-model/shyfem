
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2012,2019  Georg Umgiesser
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
! 30.03.2012	ggu	changed VERS_6_1_51
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

c curve manipulating routines
c
c contents :
c
c function minmax(x,y,n,xm,ym)		gets local min/max from a curve
c function integm(x,y,ndim,xzero,rint)	computes integral from a curve (portion)
c function rinteg(x,y,ndim,rpos,rneg)	computes integral from a curve
c
c******************************************************************
c
	function minmax(x,y,n,xm,ym)
c
c gets local min/max and returns them in vector
c
c x,y		array containing curve values
c n		number of elements in x,y
c xm,ym		array containing local min/max found
c minmxa	number of local min/max found
c
        dimension x(n),y(n)
        dimension xm(n),ym(n)
c
	nm=1
	xm(nm)=x(1)
	ym(nm)=y(1)
c
	i=0
	yo=y(1)
   10	i=i+1
	if(i.gt.n) goto 80
	if(y(i).eq.yo) goto 10
	yn=y(i)
c
	idir=1
	if(yn.lt.yo) idir=-1
	yo=yn
	jo=i
c
	do 30 j=i+1,n
	yn=y(j)
	if(yn.gt.yo) then
		if(idir.eq.1) then
			yo=yn
			jo=j
		else
			nm=nm+1
			xm(nm)=(x(jo)+x(j-1))*0.5
			ym(nm)=yo
			yo=yn
			jo=j
			idir=1
		end if
	else if(yn.lt.yo) then
		if(idir.eq.1) then
			nm=nm+1
			xm(nm)=(x(jo)+x(j-1))*0.5
			ym(nm)=yo
			yo=yn
			jo=j
			idir=-1
		else
			yo=yn
			jo=j
		end if
	end if
   30	continue
c
   80	continue
c
	nm=nm+1
	xm(nm)=x(n)
	ym(nm)=y(n)
c
	minmax=nm
c
	return
	end
c
c*****************************************************************
c
	function integm(x,y,ndim,xzero,rint)
c
c computes integral from a curve and returns an array with
c ...positive/negative portion of the integral
c
c x,y		array containing x,y coordinate of curve
c ndim		total number of elements in x,y
c xzero		array containing x value of new integral portion
c rint		array containing positive/negative portion of integral
c integm	total number of pos/neg portions found
c
        dimension x(ndim),y(ndim)
        dimension xzero(ndim),rint(ndim)
c
	nm=0
	yport=0.
c
	xo=x(1)
	yo=y(1)
c
	do i=2,ndim
c
	xn=x(i)
	yn=y(i)
	dx=xn-xo
	if(yn*yo.ge.0.) then
		ym=(yn+yo)/2.
		yport=yport+ym*dx
	else
		ddx=abs(dx*yo/(yo-yn))
		yport=yport+yo*ddx/2.
		nm=nm+1
		xzero(nm)=xo+ddx
		rint(nm)=yport
		yport=yn*(dx-ddx)/2.
	end if
	xo=xn
	yo=yn
c
	end do
c
	if(yport.ne.0.) then
		nm=nm+1
		xzero(nm)=xn
		rint(nm)=yport
	end if
c
	integm=nm
c
	return
	end
c
c*****************************************************************
c
	function rinteg(x,y,ndim,rpos,rneg)
c
c computes integral from a curve
c
c x,y		array containing x,y coordinate of curve
c ndim		total number of elements in x,y
c rpos,rneg	positive/negative portion of integral
c rinteg	total value of integral (rpos+rneg)
c
        dimension x(ndim),y(ndim)
c
	yportp=0.	!positive discharge
	yportn=0.	!negative discharge
c
	xo=x(1)
	yo=y(1)
c
	do i=2,ndim
c
	xn=x(i)
	yn=y(i)
	dx=xn-xo
	if(yn*yo.ge.0.) then
		ym=(yn+yo)/2.
		if(yn.gt.0.) then
			yportp=yportp+ym*dx
		else if(yn.lt.0.) then
			yportn=yportn+ym*dx
		else
			if(yo.ge.0.) then
				yportp=yportp+ym*dx
			else
				yportn=yportn+ym*dx
			end if
		end if
	else
		ddx=abs(dx*yo/(yo-yn))
		if(yn.ge.0.) then
			yportp=yportp+yn*(dx-ddx)/2.
			yportn=yportn+yo*ddx/2.
		else
			yportn=yportn+yn*(dx-ddx)/2.
			yportp=yportp+yo*ddx/2.
		end if
	end if
	xo=xn
	yo=yn
c
	end do
c
	rpos=yportp
	rneg=yportn
	rinteg=rpos+rneg
c
	return
	end
