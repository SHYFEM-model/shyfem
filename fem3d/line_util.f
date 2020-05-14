
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
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
! 31.03.2017	ggu	changed VERS_7_5_24
! 06.07.2018	ggu	changed VERS_7_5_48
! 31.08.2018	ggu	changed VERS_7_5_49
! 25.10.2018	ggu	changed VERS_7_5_51
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

c*******************************************************************

	subroutine check_elements(n,x,y,ieflag,ikflag)

c flags elements which are in/out-side or at border of line

	use basin

	implicit none

	integer n
	real x(n),y(n)
	integer ieflag(nel)
	integer ikflag(nkn)

	logical inpoly
	logical debug
	integer ie,k,ii,iin,ikin
	integer iout,icheck
	real xlmin,xlmax,ylmin,ylmax
	real xmin,xmax,ymin,ymax
	real xa,ya,xe(3),ye(3)

	debug = .false.

c------------------------------------------------------------
c if n < 0 flag all nodes as inside
c------------------------------------------------------------

	if( n .lt. 0 ) then
	  ieflag = 1
	  ikflag = 1
	  return
	end if

c------------------------------------------------------------
c compute min/max of line
c------------------------------------------------------------

	call xy_minmax(n,x,y,xlmin,xlmax,ylmin,ylmax)

	if( debug ) then
	  write(6,*) 'in check_elements ',n
	  do ii=1,n
	    write(6,*) ii,x(ii),y(ii)
	  end do
	end if

c------------------------------------------------------------
c flag nodes that are inside line (ikflag(k) = 1)
c------------------------------------------------------------

	do k=1,nkn
	  ikflag(k) = 0

	  xa = xgv(k)
	  ya = ygv(k)

	  if( xa .lt. xlmin ) cycle
	  if( xa .gt. xlmax ) cycle
	  if( ya .lt. ylmin ) cycle
	  if( ya .gt. ylmax ) cycle

	  if( inpoly(n,x,y,xa,ya) ) ikflag(k) = 1
	end do

	ikin = 0
	do k=1,nkn
	  if( ikflag(k) .ne. 0. ) ikin = ikin + 1
	end do

c------------------------------------------------------------
c flag elements that are inside (1), outside (-1), or on border (0)
c------------------------------------------------------------

	do ie=1,nel

	  call xy_and_minmax(ie,xe,ye,xmin,xmax,ymin,ymax)	!xe,ye not used

	  ieflag(ie) = -1

	  if( xmax .lt. xlmin ) cycle
	  if( xmin .gt. xlmax ) cycle
	  if( ymax .lt. ylmin ) cycle
	  if( ymin .gt. ylmax ) cycle

	  iin = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( ikflag(k) .ne. 0. ) iin = iin + 1
	  end do

	  !if( iin .ne. 0 ) write(6,*) 'element: ',ie,iin

	  if( iin .eq. 0 ) then			!no vertex inside line
	    ieflag(ie) = -1
	  else if( iin .eq. 3 ) then		!all vertices inside line
	    ieflag(ie) = 1
	  else
	    ieflag(ie) = 0
	  end if

c         ---------------------------------------------------
c	  if ieflag != 0 we should still check if no line point is in element
c         ---------------------------------------------------

	  if( ieflag(ie) .ne. 0 ) then
	    call linepoint_in_element(n,x,y,xe,ye,iin)
	    if( iin .gt. 0 ) ieflag(ie) = 0
	  end if
	end do

c------------------------------------------------------------
c final statistics
c------------------------------------------------------------

	iin = 0
	iout = 0
	icheck = 0
	do ie=1,nel
	  if( ieflag(ie) .eq. -1 ) iout = iout + 1
	  if( ieflag(ie) .eq.  0 ) icheck = icheck + 1
	  if( ieflag(ie) .eq. +1 ) iin = iin + 1
	end do

	write(6,*) 'min/max of line: ',xlmin,xlmax,ylmin,ylmax
	write(6,*) 'total number of nodes:              ',nkn
	write(6,*) 'nodes inside line:                  ',ikin
	write(6,*) 'total number of elements:           ',nel
	write(6,*) 'elements containing no line points: ',iout
	write(6,*) 'elements fully in line:             ',iin
	write(6,*) 'elements to be checked:             ',icheck

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c*******************************************************************

	subroutine linepoint_in_element(n,x,y,xe,ye,inside)

c checks if one of the line points is in element given by xe,ye

	implicit none

	integer n
	real x(n),y(n)
	real xe(3),ye(3)
	integer inside		!0: no point inside (return)

	integer i
	integer intri

	inside = 1

	do i=1,n
          if( intri(xe,ye,x(i),y(i)) .gt. 0 ) return
	end do

	inside = 0	!no point inside found

	end

c*******************************************************************

	subroutine xy_and_minmax(ie,x,y,xmin,xmax,ymin,ymax)

c returns x,y and min/max coordinates of vertices of element ie

        use basin

        implicit none

        integer ie
        real x(3),y(3)
        real xmin,xmax
        real ymin,ymax

        include 'param.h'

        integer ii,k

        do ii=1,3
          k=nen3v(ii,ie)
          x(ii)=xgv(k)
          y(ii)=ygv(k)
        end do

        xmin=min(x(1),x(2),x(3))
        xmax=max(x(1),x(2),x(3))
        ymin=min(y(1),y(2),y(3))
        ymax=max(y(1),y(2),y(3))

        end

c*******************************************************************

        subroutine xy_minmax(n,x,y,xmin,xmax,ymin,ymax)

c returns min/max coordinates of polygon

        implicit none

        integer n
        real x(n),y(n)
        real xmin,xmax
        real ymin,ymax

        integer i

        xmin = x(1)
        xmax = x(1)
        ymin = y(1)
        ymax = y(1)

        do i=2,n
          xmin = min(xmin,x(i))
          xmax = max(xmax,x(i))
          ymin = min(ymin,y(i))
          ymax = max(ymax,y(i))
        end do

        end

c*******************************************************************

        subroutine xy_center(n,x,y,xm,ym)

c returns center of gravity of polygon

        implicit none

        integer n
        real x(n),y(n)
        real xm,ym

        integer i

        xm = 0.
        ym = 0.

        do i=1,n
          xm = xm + x(i)
          ym = ym + y(i)
        end do

        xm = xm / n
        ym = ym / n

        end

c*******************************************************************

