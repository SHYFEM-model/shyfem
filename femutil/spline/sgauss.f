
!--------------------------------------------------------------------------
!
!    Copyright (C) 2019  Georg Umgiesser
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

c********************************************************

! revision log :
!
! 13.03.2019	ggu	changed VERS_7_5_61

c********************************************************


	subroutine gsmooth(nl,t,v,sigma,period)

c gaussian smoothing

	implicit none

	integer nl		!length of input arrays
	real t(nl)		!time values (independent)
	real v(nl)		!dipendent values
	real sigma		!standard deviation for smoothing
	real period		!if > 0 signal is periodic with period

c if signal is periodic, period is period of signal
c t(n) - t(1) < period must be true

	integer ndim
	parameter(ndim=10000)
	
	real taux(-ndim:2*ndim)
	real vaux(-ndim:2*ndim)

	logical bperiod
	integer i,j
	integer jmin,jmax
	real a,b
	real t0,tnew,dt,weight
	real pi,eps
	real val
	real rkern,rknew,tkern

c----------------------------------------------------
c check dimensions
c----------------------------------------------------

	if( nl .gt. ndim ) stop 'error stop gsmooth: dimension error'

c----------------------------------------------------
c set up circular array
c----------------------------------------------------

	bperiod = period .gt. 0.
	dt = t(ndim) - t(1)
	if( bperiod .and. period .le. dt ) then
	  write(6,*) 'signal is flagged periodic, but period'
	  write(6,*) 'is smaller or equal to time series'
	  write(6,*) 'period, dt : ',period,dt
	  stop 'error stop gsmooth: period'
	endi f

	do i=1,nl
	  val = t(i)
	  taux(i) = val
	  taux(i+nl) = val + period
	  taux(i-nl) = val - period
	  val = v(i)
	  vaux(i) = val
	  vaux(i+nl) = val
	  vaux(i-nl) = val
	end do

c----------------------------------------------------
c set up gaussian kernel parameters
c----------------------------------------------------

        pi = 4.*atan(1.)
        eps = 1.e-7

        a = 1. / ( sqrt(2.*pi) * sigma )
        b = - 1. / ( 2 * sigma * sigma )

c----------------------------------------------------
c convolution
c----------------------------------------------------

	do i=1,nl

c	  -------------------------------------------
c	  initialize
c	  -------------------------------------------

	  tkern = 0.		!total weight
	  val = 0.		!new value
	  t0 = taux(i)		!t value of central point

	  jmin = i - nl
	  jmax = i + nl

	  if( .not. bperiod ) then
	    jmin = max(1,jmin)
	    jmax = min(nl,jmax)
	  end if

c	  write(6,*) i,i,weight,rkern,tkern,val

c	  -------------------------------------------
c	  accumulate forward
c	  -------------------------------------------

	  rkern = a
	  j = i
	  do while( rkern .gt. eps .and. j .lt. jmax )
	    j = j + 1
	    dt = taux(j) - taux(j-1)
	    tnew = taux(j) - t0
	    rknew = a * exp( b * tnew * tnew )
	    weight = dt * ( rkern + rknew )
	    tkern = tkern + weight
	    val = val + dt * ( rkern * vaux(j-1) + rknew * vaux(j) )
	    rkern = rknew
c	    write(6,*) i,j,weight,rkern,tkern,val
	  end do

c	  write(6,*) i,j-i,nl,rkern

c	  -------------------------------------------
c	  accumulate backward
c	  -------------------------------------------

	  rkern = a
	  j = i
	  do while( rkern .gt. eps .and. j .gt. jmin )
	    j = j - 1
	    dt = taux(j+1) - taux(j)
	    tnew = taux(j) - t0
	    rknew = a * exp( b * tnew * tnew )
	    weight = dt * ( rkern + rknew )
	    tkern = tkern + weight
	    val = val + dt * ( rkern * vaux(j+1) + rknew * vaux(j) )
	    rkern = rknew
c	    write(6,*) i,j,weight,rkern,tkern,val
	  end do

c	  write(6,*) i,i-j,nl,rkern
c	  write(6,*) '*** ',i,val,tkern

c	  -------------------------------------------
c	  end of accumulation -> set value
c	  -------------------------------------------

	  if( tkern .gt. 0. ) then
	    v(i) = val / tkern
	  else
	    v(i) = vaux(i)
	  end if

	end do

c----------------------------------------------------
c end of routine
c----------------------------------------------------

	end

c********************************************************

