
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

c*********************************************************************
c
	subroutine setlst(ip,rkey,n,rflag)
c
	implicit none
c
c arguments
	integer ip(1),n
	real rkey(1),rflag
c common
	integer nmax,nins
	real rlast
	common /lstloc/ nmax,nins,rlast
c local
	integer i
c
	nmax=n
	rlast=rflag
	nins=0
c
	do i=1,n
	  ip(i)=0
	  rkey(i)=rflag
	end do
c
	return
	end
c
c*********************************************************************
c
	subroutine inslst(ip,rkey,ipact,rkact)
c
	implicit none
c
c arguments
	integer ip(1),ipact
	real rkey(1),rkact
c common
	integer nmax,nins
	real rlast
	common /lstloc/ nmax,nins,rlast
c local
	integer i
c
	if(rkact.le.rlast) return
c
	nins=nins+1
c
	do i=nmax-1,1,-1
		if(rkact.le.rkey(i)) goto 1
		rkey(i+1)=rkey(i)
		ip(i+1)=ip(i)
	end do
    1	continue
	rkey(i+1)=rkact
	ip(i+1)=ipact
c
	rlast=rkey(nmax)
c
	return
	end
c
	
c
c*********************************************************************
c
	subroutine maxlst(n)
c
	implicit none
c
c arguments
	integer n
c common
	integer nmax,nins
	real rlast
	common /lstloc/ nmax,nins,rlast
c
	n=nins
	if(nins.gt.nmax) n=nmax
c
	return
	end
