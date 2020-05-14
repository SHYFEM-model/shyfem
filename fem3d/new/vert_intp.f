
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2019  Georg Umgiesser
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

c*************************************************************************

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2019	ggu	changed VERS_7_5_60

c*************************************************************************

	subroutine intp_vert(nl1,zb1,var1,nl2,zb2,var2)

c vertical interpolation of variables from one grid to another

	implicit none

	integer nl1		!number of levels of first grid
	real zb1(0:nl1+1)	!depth of bottom of vertical boxes
	real var1(nl1+1)	!value of variable at center of box
	integer nl2		!number of levels of second grid
	real zb2(0:nl2)		!depth of bottom of vertical boxes
	real var2(nl2)		!value of variable at center of box

c values are interpolated from first grid to second grid
c z levels refer to bottom of each grid (zb1(nl1) is total depth of column)
c zb(0) is surface ... normally 0
c variables are considered at center of box
c values are considered to be constant for every box
c
c ATTENTION: arrays zb1, var1 MUST be dimensioned one element bigger
c	then the available data. The last element is altered by this
c	subroutine, but should be of no concern for the calling program
c
c output is var1, all other variables are input values

	logical bmiss
	integer l,j,ltop1,lbot1
	real ztop2,ztop1,zbot2,zbot1
	real ztop,zbot
	real vint1,vint2,fact
	real val
	logical bcons

	zb1(nl1+1) = zb2(nl2)
	var1(nl1+1) = var1(nl1)
	bcons = .true.			!conserve total quantity
	bcons = .false.			!do not conserve total quantity

	ltop1 = 0
	vint2 = 0.

	do l=1,nl2
	  ztop2 = zb2(l-1)
	  zbot2 = zb2(l)

	  do while( ltop1 .lt. nl1 .and. zb1(ltop1+1) .le. ztop2 )
	    !write(6,*) 'adjourning top depth: ',ltop1,zb1(ltop1+1),ztop2
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

	  if( ztop1 .gt. ztop2 .or. zbot1 .lt. zbot2 ) goto 99
	  if( bmiss ) goto 98

	  !write(6,*) l,ltop1,lbot1,ztop2,zbot2,zb1(lbot1)

	  val = 0.
	  do j=ltop1+1,lbot1
	      ztop = max(zb1(j-1),ztop2)
	      zbot = min(zb1(j),zbot2)
	      val = val + var1(j) * ( zbot - ztop )
	  end do

	  vint2 = vint2 + val			!integrated value
	  var2(l) = val / (zbot2-ztop2)
	  ltop1 = lbot1 - 1

	end do

	zb1(nl1+1) = 0.
	var1(nl1+1) = 0.

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

	return
   98	continue
	stop 'error stop intp_vert: missing value not possible'
   99	continue
	write(6,*) ztop1,ztop2
	write(6,*) zbot1,zbot2
	stop 'error stop intp_vert: interval in (1) must include (2)'
	end

c*************************************************************************

	subroutine intp_test

	implicit none

	integer ndim
	parameter (ndim=100)

	integer l,nl1,nl2
	real ztot1,ztot2

	real zb1(0:ndim)
	real zb2(0:ndim)
	real var1(ndim)
	real var2(ndim)

	call intp_make_test(0,ndim,8,5.,10,4.,zb1,var1,zb2,var2)
	call intp_make_test(0,ndim,8,5.,40,1.,zb1,var1,zb2,var2)
	call intp_make_test(0,ndim,12,5.,20,3.,zb1,var1,zb2,var2)
	call intp_make_test(30,ndim,12,5.,10,6.,zb1,var1,zb2,var2)
	call intp_make_test(0,ndim,14,5.,10,7.,zb1,var1,zb2,var2)
	call intp_make_test(0,ndim,14,5.,5,14.,zb1,var1,zb2,var2)
	call intp_make_test(40,ndim,14,5.,2,35.,zb1,var1,zb2,var2)
	call intp_make_test(0,ndim,5,14.,14,5.,zb1,var1,zb2,var2)
	call intp_make_test(0,ndim,2,35.,14,5.,zb1,var1,zb2,var2)

	call intp_make_test(50,ndim,8,5.,9,4.,zb1,var1,zb2,var2)
	call intp_make_test(60,ndim,8,5.,11,4.,zb1,var1,zb2,var2)

	end
	
c*************************************************************************

	subroutine intp_make_test(iunit,ndim
     +			,nl1,dz1,nl2,dz2,zb1,var1,zb2,var2)

	implicit none

	integer iunit		!output unit (0 -> monitor)
	integer ndim		!dimensions of arrays
	integer nl1		!number of levels of first grid
	real dz1		!thickness of levels of first grid
	integer nl2		!number of levels of second grid
	real dz2		!thickness of levels of second grid
	real zb1(0:nl1)		!depth of bottom of vertical boxes
	real var1(nl1)		!value of variable at center of box
	real zb2(0:nl2)		!depth of bottom of vertical boxes
	real var2(nl2)		!value of variable at center of box

	call intp_make_array(ndim,nl1,dz1,zb1,var1)
	call intp_make_array(ndim,nl2,dz2,zb2,var2)

	call intp_vert(nl1,zb1,var1,nl2,zb2,var2)
	call intp_total(nl1,zb1,var1,nl2,zb2,var2)

	if( iunit .gt. 0 ) then
	  call intp_print(iunit,nl1,zb1,var1,nl2,zb2,var2)
	end if

	end

c*************************************************************************

	subroutine intp_make_array(ndim,nl,dz,zb,var)

	implicit none

	integer ndim,nl
	real dz
	real zb(0:ndim)
	real var(ndim)

	integer l

	do l=1,nl
	  zb(l) = l*dz
	  var(l) = 10.*sqrt(float(l))
	end do
	zb(0) = 0.

	do l=nl+1,ndim
	  zb(l) = 0.
	  var(l) = 0.
	end do

	end

c*************************************************************************

	subroutine intp_print(iunit,nl1,zb1,var1,nl2,zb2,var2)

	implicit none

	integer iunit
	integer nl1		!number of levels of first grid
	real zb1(0:nl1)		!depth of bottom of vertical boxes
	real var1(nl1)		!value of variable at center of box
	integer nl2		!number of levels of second grid
	real zb2(0:nl2)		!depth of bottom of vertical boxes
	real var2(nl2)		!value of variable at center of box

	integer l

	do l=1,nl1
	  write(iunit,*) zb1(l-1),var1(l)
	  write(iunit,*) zb1(l),var1(l)
	end do

	do l=1,nl2
	  write(iunit+1,*) zb2(l-1),var2(l)
	  write(iunit+1,*) zb2(l),var2(l)
	end do

	end

c*************************************************************************

	subroutine intp_total(nl1,zb1,var1,nl2,zb2,var2)

	implicit none

	integer nl1		!number of levels of first grid
	real zb1(0:nl1)		!depth of bottom of vertical boxes
	real var1(nl1)		!value of variable at center of box
	integer nl2		!number of levels of second grid
	real zb2(0:nl2)		!depth of bottom of vertical boxes
	real var2(nl2)		!value of variable at center of box

	integer l,n
	real ztot1,ztot2

	n = max(nl1,nl2)

	ztot1 = 0.
	ztot2 = 0.
	do l=1,n
	  !write(6,*) zb1(l),zb2(l),var1(l),var2(l)
	  ztot1 = ztot1 + var1(l) * (zb1(l)-zb1(l-1))
	  ztot2 = ztot2 + var2(l) * (zb2(l)-zb2(l-1))
	end do

	write(6,*) 'total: ',nl1,nl2,ztot1,ztot2,ztot1/zb1(nl1),ztot2/zb2(nl2)

	end

c*************************************************************************

	subroutine intp_debora

	implicit none

	integer ndim
	parameter (ndim=30)

	real z1(0:ndim+1)
	real u1(ndim+1)
	real v1(ndim+1)
	real z2(0:ndim+1)
	real u2(ndim+1)
	real v2(ndim+1)

	integer n1,n2,i,j
	character*40 file

	n1 = ndim
	file='uvin.dat'
	open(1,file=file)
	read(1,*)
	do i=1,n1
	  read(1,*) j,z1(i),u1(i),v1(i)
	  if( i .ne. j ) stop 'i != j'
	end do
	close(1)

	n2 = 16
	file='zfemout.dat'
	open(1,file=file)
	read(1,*)
	do i=1,n2
	  read(1,*) z2(i)
	end do
	close(1)

	z1(0) = 0.
	z2(0) = 0.

	do i=1,n2
	  write(6,*) i,z2(i)
	end do

	do i=1,n1
	  write(6,*) i,z1(i),u1(i),v1(i)
	end do
	call intp_vert(n1,z1,u1,n2,z2,u2)
	call intp_vert(n1,z1,v1,n2,z2,v2)
	do i=1,n2
	  write(6,*) i,z2(i),u2(i),v2(i)
	end do

	end

c*************************************************************************
c*************************************************************************
c	program intp_main
c	!call intp_test
c	call intp_debora
c	end
c*************************************************************************

