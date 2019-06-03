
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

c utility routines for fem model
c
c contents :
c
c subroutine baric(ie,x,y)		finds baricentre of element
c subroutine coord(k,x,y)		returns coordinates of node k (internal)
c function area_element(ie)		area for element ie (internal)
c
c function ipext(k)             returns extern node number
c function ieext(k)             returns extern element number
c function ipint(k)             returns intern node number
c function ieint(k)             returns intern element number
c
c subroutine n2int(n,nnodes,berror)	converts external nodes to internal
c subroutine e2int(n,nelems,berror)	converts external elements to internal
c
c revision log :
c
c 02.06.1997	ggu	eliminated extint,extinw,exinel (not needed)
c 24.02.1999	ggu	new subroutine n2int
c 19.11.1999	ggu	new subroutine e2int
c 27.01.2009	ggu	new subroutine coord
c 23.03.2010	ggu	changed v6.1.1
c 07.07.2011	ggu	new subroutine area_element()
c 14.07.2011	ggu	changed VERS_6_1_27
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.05.2015	ggu	changed VERS_7_1_10
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 28.04.2016	ggu	changed VERS_7_5_9
c 06.07.2018	ggu	changed VERS_7_5_48
c 16.02.2019	ggu	changed VERS_7_5_60
c
c*******************************************************

	subroutine baric(ie,x,y)

c finds baricentre of element
c
c ie		number of element
c x,y		coordinates of baricentre (return value)

	use basin

	implicit none

c arguments
	integer ie
	real x,y
c local variables
	integer i,kkk
	real xb,yb

	xb=0.
	yb=0.
	do i=1,3
	   kkk=nen3v(i,ie)
	   xb=xb+xgv(kkk)
	   yb=yb+ygv(kkk)
	end do

	x=xb/3.
	y=yb/3.

	end

c***************************************************************

        function area_element(ie)

c area for element ie

	use basin

        implicit none

c arguments
        real area_element
        integer ie
c local
        integer kn1,kn2,kn3
        real*8 x1,x2,x3,y1,y2,y3
        real*8 half

        half = 0.5

        kn1=nen3v(1,ie)
        kn2=nen3v(2,ie)
        kn3=nen3v(3,ie)

        x1=xgv(kn1)
        y1=ygv(kn1)
        x2=xgv(kn2)
        y2=ygv(kn2)
        x3=xgv(kn3)
        y3=ygv(kn3)

        area_element = half * ( (x2-x1) * (y3-y1) - (x3-x1) * (y2-y1) )

        end

c***************************************************************

	subroutine coord(k,x,y)

c returns coordinates of node k (internal)

	use basin

	implicit none

	integer k
	real x,y

	x = xgv(k)
	y = ygv(k)

	end

c***************************************************************
c
	function ipext(k)
c
c returns external node number
c
c k     internal node number
c ipext external node number, 0 if error
c
	use basin

	implicit none
	integer ipext
	integer k

	if(k.lt.1.or.k.gt.nkn) then
		ipext=0
	else
		ipext=ipv(k)
	end if
c
	return
	end
c
c***************************************************************
c
	function ieext(k)
c
c returns external element number
c
c k     internal element number
c ieext external element number, 0 if error
c
	use basin

	implicit none
	integer ieext
	integer k
c
	if(k.lt.1.or.k.gt.nel) then
		ieext=0
	else
		ieext=ipev(k)
	end if
c
	return
	end
c
c***************************************************************
c
	function ipint(k)
c
c returns internal node number
c
c k     external node number
c ipint internal node number, 0 if error
c
	use basin

	implicit none
	integer ipint
	integer k,i
c
	do i=1,nkn
	if(ipv(i).eq.k) goto 1
	end do
	i=0
    1   continue
	ipint=i
c
	return
	end
c
c***************************************************************
c
	function ieint(k)
c
c returns internal element number
c
c k     external element number
c ieint internal element number, 0 if error
c
	use basin

	implicit none
	integer ieint
	integer k,i
c
	do i=1,nel
	if(ipev(i).eq.k) goto 1
	end do
	i=0
    1   continue
	ieint=i
c
	return
	end

c***************************************************************

	subroutine n2int(n,nnodes,berror)

c converts external nodes numbers to internal ones

	implicit none

	integer n		!total number of nodes
	integer nnodes(n)	!node numbers
	logical berror		!true on return if error

	integer i,ke,ki
	integer ipint

	berror = .false.

	do i=1,n
	  ke = nnodes(i)
	  if( ke > 0 ) then
	    ki = ipint(ke)
	    if( ki <= 0 ) then
	      write(6,*) 'No such node : ',ke
	      berror = .true.
	    end if
	  else
	    !ki = ke
	    ki = 0
	  end if
	  nnodes(i) = ki
	end do

	end

c***************************************************************

	subroutine e2int(n,nelems,berror)

c converts external element numbers to internal ones

	implicit none

	integer n		!total number of elements
	integer nelems(n)	!element numbers
	logical berror		!true on return if error

	integer i,ke,ki
	integer ieint

	berror = .false.

	do i=1,n
	  ke = nelems(i)
	  ki = ieint(ke)
	  if( ki .le. 0 .and. ke .gt. 0 ) then
	    write(6,*) 'No such element : ',ke
	    berror = .true.
	  end if
	  nelems(i) = ki
	end do

	end

c***************************************************************

