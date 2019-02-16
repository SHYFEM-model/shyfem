
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

c routines for vertical interpolation
c
c contents :
c
c revision log :
c
c 04.10.2012    ggu     copied from newbsig.f
c 04.06.2013    ggu     bug fix: restore altered values zb1_save,var1_save
c 28.06.2016    ggu     bug fix: handel situation with dh == 0
c
c*****************************************************************

	subroutine intp_vert(bcons,nl1,zb1,var1,nl2,zb2,var2)

c vertical interpolation of variables from one grid to another

	implicit none

	logical bcons		!conserve total quantity
	integer nl1		!number of levels of first grid
	real zb1(0:nl1+1)	!depth of bottom of vertical boxes
	real var1(nl1+1)	!value of variable at center of box
	integer nl2		!number of levels of second grid
	real zb2(0:nl2)		!depth of bottom of vertical boxes
	real var2(nl2)		!value of variable at center of box (out)

c values are interpolated from first grid to second grid
c z levels refer to bottom of each grid (zb1(nl1) is total depth of column)
c zb(0) is surface ... normally 0
c variables are considered at center of box
c values are considered to be constant for every box
c bcons controlls if total quantity is conserved or not
c
c ATTENTION: arrays zb1, var1 MUST be dimensioned one element bigger
c	then the available data. The last element is accessed and altered 
c	by this subroutine, but at return is restored to its original value
c
c output is var2, all other variables are input values

	logical bmiss,bdebug
	integer l,j,ltop1,lbot1
	real ztop2,ztop1,zbot2,zbot1
	real ztop,zbot
	real vint1,vint2,fact
	real val
	real zb1_save,var1_save
	real dz,dh

c---------------------------------------------------------
c handle special situations
c---------------------------------------------------------

	if( nl1 < 1 .or. nl2 < 1 ) then
	  write(6,*) 'nl1,nl2: ',nl1,nl2
	  stop 'error stop intp_vert: nl1,nl2 < 1'
	else if( nl1 == 1 ) then	!only one layer in input
	  var2 = var1(1)
	  return
	end if

c---------------------------------------------------------
c initialize variables
c---------------------------------------------------------

	zb1_save = zb1(nl1+1)
	var1_save = var1(nl1+1)

	zb1(nl1+1) = zb2(nl2)
	var1(nl1+1) = var1(nl1)

	ltop1 = 0
	vint2 = 0.	!total content in second array

	bdebug = .false.

c---------------------------------------------------------
c loop over second array and interpolate onto it
c---------------------------------------------------------

	if( bdebug ) write(66,*) '----------------------------------'

	do l=1,nl2
	  ztop2 = zb2(l-1)
	  zbot2 = zb2(l)

	  if( bdebug ) write(66,*) 'l,ztop2,zbot2: ',l,ztop2,zbot2

	  do while( ltop1 .lt. nl1 .and. zb1(ltop1+1) .le. ztop2 )
	    ltop1 = ltop1 + 1
	  end do

	  bmiss = .false.
	  do lbot1=ltop1+1,nl1+1
	    if( zb1(lbot1) .ge. zbot2 ) goto 1
	  end do
	  lbot1 = nl1
	  bmiss = .true.
    1	  continue

	  ztop1 = zb1(ltop1)
	  zbot1 = zb1(lbot1)

	  if( bdebug ) write(66,*) 'l,ltop1,lbot1: ',l,ltop1,lbot1
	  if( bdebug ) write(66,*) 'l,ztop1,zbot1: ',l,ztop1,zbot1

	  if( ztop1 .gt. ztop2 .or. zbot1 .lt. zbot2 ) goto 99
	  if( bmiss ) goto 98

	  val = 0.
	  dh = 0.
	  do j=ltop1+1,lbot1
	    ztop = max(zb1(j-1),ztop2)
	    zbot = min(zb1(j),zbot2)
	    dz = max(0.,zbot-ztop)
	    dh = dh + dz
	    val = val + var1(j) * dz
	    if( bdebug ) write(66,*) 'j,ztop,zbot,val: ',j,ztop,zbot,val
	  end do

	  vint2 = vint2 + val			!integrated value
	  if( dh > 0. ) then
	    !var2(l) = val / (zbot2-ztop2)
	    var2(l) = val / dh
	  else
	    var2(l) = var1(lbot1)	!anything valid
	  end if
	  ltop1 = lbot1 - 1

	  if( bdebug ) write(66,*) 'l,var2: ',l,var2(l)

	end do

	if( bdebug ) write(66,*) '----------------------------------'

c---------------------------------------------------------
c reset modified values
c---------------------------------------------------------

	zb1(nl1+1) = zb1_save
	var1(nl1+1) = var1_save

c---------------------------------------------------------
c if bcons we have to adapt the interpolated values to conserve total value
c---------------------------------------------------------

	if( bcons ) then	!must conserve total content of scalar
	  vint1 = 0.
	  do l=1,nl1
	    ztop = zb1(l-1)
	    zbot = zb1(l)
	    vint1 = vint1 + var1(l) * ( zbot - ztop )
	  end do

	  if( vint1 .eq. vint2 ) return

	  fact = vint1 / vint2
	  do l=1,nl2
	    var2(l) = fact * var2(l)
	  end do
	end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	return
   98	continue
	stop 'error stop intp_vert: missing value not possible'
   99	continue
	write(6,*) ztop1,ztop2
	write(6,*) zbot1,zbot2
	stop 'error stop intp_vert: interval in (1) must include (2)'
	end

c*****************************************************************

        subroutine intp_aver(n,horig,vorig,femval)

        implicit none

        integer n
        real horig(0:n)
        real vorig(n)
        real femval

        integer i
        real hacu,acu,h

        acu = 0.
        hacu = 0.
        do i=1,n
          h = horig(i) - horig(i-1)
          hacu = hacu + h
          acu = acu + h * vorig(i)
        end do

        femval = acu / hacu

        end

c*****************************************************************

	subroutine intp_vert_test

	implicit none

	integer ndim
	parameter (ndim=100)
	real zb1(0:ndim+1)
	real var1(ndim+1)
	real zb2(0:ndim)
	real var2(ndim)
	logical bcons
	integer nl1,nl2,i

	nl1 = 1
	zb1(0) = 0.
	zb1(1) = 10.

	nl2 = 1
	zb2(0) = 0.
	zb2(1) = 8.

	var1(1) = 30.

	bcons = .false.
	call intp_vert(bcons,nl1,zb1,var1,nl2,zb2,var2)
	write(6,*) bcons,var1(1),var2(1)

	bcons = .true.
	call intp_vert(bcons,nl1,zb1,var1,nl2,zb2,var2)
	write(6,*) bcons,var1(1),var2(1)

	nl1 = 2
	zb1(2) = 20.
	var1(2) = 40.

	bcons = .false.
	call intp_vert(bcons,nl1,zb1,var1,nl2,zb2,var2)
	write(6,*) bcons,(var1(i),i=1,nl1),(var2(i),i=1,nl2)

	bcons = .true.
	call intp_vert(bcons,nl1,zb1,var1,nl2,zb2,var2)
	write(6,*) bcons,(var1(i),i=1,nl1),(var2(i),i=1,nl2)

	end

c*****************************************************************

	!program main_intp_vert_test
	!call intp_vert_test
	!end

c*****************************************************************

