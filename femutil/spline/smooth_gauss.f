
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

c smoothing program
c
c contents :
c
c revision log :
c 
c 01.02.2002    ggu     handle periodic line
c 06.08.2002    ggu     read also depth value
c 23.09.2004    ggu     adapted for malta, bug fix
c 19.10.2005    ggu     documentation and description
c 02.12.2011    ggu     bug fix in intpdep() and reduce()
c 02.12.2011    ggu     use depth also for smoothing (change in distxy())
c 27.05.2012    ggu     renamed ndim to nsdim in smooth()
c
c********************************************************

	subroutine gkernel(sigma,ndim,ngk,gk)

c gaussian smoothing - creates kernel to be used with gsmooth

	real sigma
	integer ndim
	integer ngk
	real gk(-ndim:ndim)

	integer i
	real pi,eps,a,b,x,val

        pi = 4.*atan(1.)
        eps = 1.e-7

        a = 1. / ( sqrt(2.*pi) * sigma )
        b = 1. / ( 2 * sigma * sigma )

	do i=0,ndim
          x = i
          val = a * exp( -x**2 * b )
	  gk(i) = val
	  gk(-i) = val
	  if( val .lt. eps ) exit
        end do
	ngk = i-1

	end

c********************************************************

	subroutine gsmooth(ndim,nl,rt,raux,ngk,gk)

c gaussian smoothing - only good for regular point spacing

	integer ndim
	integer nl
	real rt(1)
	real raux(-ndim:2*ndim)
	integer ngk
	real gk(-ndim:ndim)

	integer i,j
	real val

c set up circular array

	do i=1,nl
	  val = rt(i)
	  raux(i) = val
	  raux(i+nl) = val
	  raux(i-nl) = val
	end do

c convolution

	do i=1,nl
	  val = 0.
	  do j=-ngk,ngk
	    val = val + raux(i+j) * gk(j)
	  end do
	  rt(i) = val
	end do

	end

c********************************************************

	subroutine grsmooth(nl,sigma,rt,raux,dxy,ht,bperiod)

c gaussian smoothing - good for irregular point spacing

	integer nl
	real sigma
	real rt(nl)
	real raux(-nl:2*nl)
	real dxy(-nl:2*nl)
	real ht(nl)
	logical bperiod

	integer i,j
	real val
	real rkern,tkern
	real rintv

c set up circular array

	do i=1,nl
	  val = rt(i)
	  raux(i) = val
	  raux(i+nl) = val
	  raux(i-nl) = val
	end do

c set up gaussian kernel parameters

        pi = 4.*atan(1.)
        eps = 1.e-7

        a = 1. / ( sqrt(2.*pi) * sigma )
        b = - 1. / ( 2 * sigma * sigma )

c convolution

	do i=1,nl

	  rintv = 0.5 * ( dxy(i-1) + dxy(i) )
	  rkern = a
	  weight = rkern * rintv
	  tkern = weight
	  val = raux(i) * weight
	  jmin = i - nl
	  jmax = i + nl
	  if( .not. bperiod ) then
	    jmin = max(2,jmin)
	    jmax = min(nl-1,jmax)
	  end if

	  rkern = 1.
	  dist = 0.
	  j = i
	  do while( rkern .gt. eps .and. j .lt. jmax )
	    j = j + 1
	    dist = dist + dxy(j-1)
	    rintv = 0.5 * ( dxy(j-1) + dxy(j) )
	    rkern = a * exp( b * dist * dist )
	    weight = rkern * rintv
	    tkern = tkern + weight
	    val = val + raux(j) * weight
	  end do

	  rkern = 1.
	  dist = 0.
	  j = i
	  do while( rkern .gt. eps .and. j .gt. jmin )
	    j = j - 1
	    dist = dist + dxy(j)
	    rintv = 0.5 * ( dxy(j-1) + dxy(j) )
	    rkern = a * exp( b * dist * dist )
	    weight = rkern * rintv
	    tkern = tkern + weight
	    val = val + raux(j) * weight
	  end do

	  if( tkern .gt. 0. .and. ht(i) .ge. 0. ) then
	    rt(i) = val / tkern
	  else
	    rt(i) = raux(i)
	  end if

	end do

	end

c********************************************************

