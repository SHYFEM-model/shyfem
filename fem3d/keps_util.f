
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

c routines for k-epsilon model
c
c revision log :
c
c 25.02.1999	ggu	finished writing routines
c
c notes :
c
c u,v,rho are evaluated at the center of the layer
c other variables (k,eps,rkm,rkt) at the top and bottom of the layers
c
c********************************************************************

	subroutine tridagd(a,b,c,r,u,ud,gam,n)

c solves tridiagonal system

	implicit none

	integer n
	double precision a(n),b(n),c(n),r(n),ud(n),gam(n)
	real u(n)

	integer j
	double precision bet

	do j=1,n
	  ud(j) = u(j)
	end do

	if(b(1).eq.0.) stop 'error stop tridag: matrix singular'
	bet=b(1)
	ud(1)=r(1)/bet

	do j=2,n
	  gam(j)=c(j-1)/bet
	  bet=b(j)-a(j)*gam(j)
	  if(bet.eq.0.) stop 'error stop tridag: matrix singular'
	  ud(j)=(r(j)-a(j)*ud(j-1))/bet
	end do

	do j=n-1,1,-1
	  ud(j)=ud(j)-gam(j+1)*ud(j+1)
	end do

	do j=1,n
	  u(j) = ud(j)
	end do

	end

c************************************************************************

	subroutine tridag(a,b,c,r,u,gam,n)

c solves tridiagonal system

	implicit none

	integer n
	real a(n),b(n),c(n),r(n),u(n),gam(n)

	integer j
	real bet

	if(b(1).eq.0.) stop 'error stop tridag: matrix singular'
	bet=b(1)
	u(1)=r(1)/bet

	do j=2,n
	  gam(j)=c(j-1)/bet
	  bet=b(j)-a(j)*gam(j)
	  if(bet.eq.0.) stop 'error stop tridag: matrix singular'
	  u(j)=(r(j)-a(j)*u(j-1))/bet
	end do

	do j=n-1,1,-1
	  u(j)=u(j)-gam(j+1)*u(j+1)
	end do

	end

c*********************************************************************

	subroutine makedz(lmax,dl,dzk,dzr)

c makes level differences from dl (bottom depth of single levels)

	implicit none

	integer lmax		!number of levels [1...lmax]
	real dl(1:lmax)		!bottom depth of single levels
	real dzk(lmax)		!z difference between k (eps) levels
	real dzr(0:lmax)	!z difference between rho (u,v) levels

	integer l
	real z0

c z0 is height of surface -> if highly variable -> pass into routine
c
c dzk == dz        dzr == dZ  (Langland and Liou, MWR, May 1996, 905-918)

	z0 = 0.

	dzk(1) = dl(1) - z0
	dzr(0) = 1.			!dummy
	do l=2,lmax-1
	  dzk(l) = dl(l) - dl(l-1)
	  dzr(l-1) = 0.5 * ( dzk(l) + dzk(l-1) )
	end do
	dzk(l) = dl(l) - dl(l-1)	!last level -> l = lmax
	dzr(l-1) = 0.5 * ( dzk(l) + dzk(l-1) )
	dzr(l) = 1.			!dummy

	do l=1,lmax-1
	  dzk(l) = 1. / dzk(l)
	  dzr(l) = 1. / dzr(l)
	end do
	dzk(l) = 1. / dzk(l)		!last level -> l = lmax

	end

c*********************************************************************

	subroutine ksmooth(lmax,k,eps,rkm,rkt,kold,epsold,rkmold,rktold)

c smoothing of keps quantities

	implicit none

	integer lmax
	real k(0:1),eps(0:1)
	real rkm(0:1),rkt(0:1)
	real kold(0:1),epsold(0:1)
	real rkmold(0:1),rktold(0:1)

	call smooth(lmax,k,kold)
	call smooth(lmax,eps,epsold)
	call smooth(lmax,rkm,rkmold)
	call smooth(lmax,rkt,rktold)

	end

c************************************************************************

	subroutine smooth(lmax,var,varold)

c smoothing in time

	implicit none

	integer lmax
	real var(0:lmax)
	real varold(0:lmax)

	integer l
	real alpha,alpha1

	alpha = 0.6		!weighting of old time level
	alpha1 = 1. - alpha

	do l=1,lmax-1
	  var(l) = alpha1 * var(l) + alpha * varold(l)
	  varold(l) = var(l)
	end do

	end

c********************************************************************

	subroutine hmixk(lmax,k,lmix)

c finds depth of mixed layer

	implicit none

	integer lmax,lmix
	real k(0:lmax)

	integer l
	real d1,d2

	do l=lmax,1,-1
	  if( k(l) .gt. 1.e-5 ) then
	    lmix = l
	    return
	  end if
	end do

	lmix = 0

	end

c********************************************************************

	subroutine hmix(lmax,rho,rhoin,lmix)

c finds depth of mixed layer

	implicit none

	integer lmax,lmix
	real rho(25)
	real rhoin(25)

	integer l
	real d1,d2

	write(91,'(25i3)') (nint(1000.*(rho(l)-rhoin(l))),l=1,25)

	do l=1,lmax
	  d1 = abs( rho(l) - rhoin(l) )
	  d2 = abs( rho(l) - rhoin(1) )
	  d2 = abs( rho(l) - rho(1) )
	  if( d1 .le. 0.01 * d2 ) then
	    lmix = l
	    return
	  end if
	end do

	lmix = lmax

	end

c********************************************************************

	subroutine prvar(iunit,lmax,t,val,fact,const)

c prints val to file

	implicit none

	integer iunit
	integer lmax
	real t
	real val(1)
	real fact,const

	integer i

	write(iunit,*) t,lmax,fact
	write(iunit,*) (const+val(i),i=1,lmax)

	end

c********************************************************************

