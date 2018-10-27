
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

c utility routines for operations on vectors
c
c
c subroutine mima(xx,n,xmin,xmax)	min/max of vector
c subroutine addvec(a1,a2,n)		adds two vectors
c subroutine mulvec(a1,a2,n)		multiplies two vectors
c subroutine zervec(a1,n)		initializes vector
c
c*******************************************
c
	subroutine mima(xx,n,xmin,xmax)
c
c computes min/max of vector
c
c xx		vector
c n		dimension of vector
c xmin,xmax	min/max value in vector
c
        implicit none
c
        integer n,i
        real xx(n)
        real xmin,xmax,x
c
	xmax=xx(1)
	xmin=xmax
c
	do i=1,n
          x=xx(i)
          if(x.gt.xmax) xmax=x
          if(x.lt.xmin) xmin=x
	end do
c
	return
	end
c
c*******************************************
c
	subroutine addvec(a1,a2,n)
c
c adds two vectors
c
c results are passed back in a1
c
c a1		first array
c a2		second array
c n		dimension of vectors
c
        implicit none
c
        integer n,i
        real a1(n),a2(n)
c
	do i=1,n
          a1(i)=a1(i)+a2(i)
	end do
c
	return
	end
c
c*******************************************
c
	subroutine mulvec(a1,a2,n)
c
c multiplies two vectors
c
c results are passed back in a1
c
c a1		first array
c a2		second array
c n		dimension of vectors
c
        implicit none
c
        integer n,i
        real a1(n),a2(n)
c
	do i=1,n
          a1(i)=a1(i)*a2(i)
	end do
c
	return
	end
c
c*******************************************
c
	subroutine zervec(a1,n)
c
c initializes vector
c
c a1		vector
c n		dimension of vector
c
        implicit none
c
        integer n,i
        real a1(n)
c
	do i=1,n
          a1(i)=0.
	end do
c
	return
	end
