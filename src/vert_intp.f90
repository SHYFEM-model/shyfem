!
! $Id: subvintp.f,v 1.3 2009-03-24 17:49:19 georg Exp $
!
! routines for vertical interpolation
!
! contents :
!
! revision log :
!
! 04.10.2012    ggu     copied from newbsig.f
! 04.06.2013    ggu     bug fix: restore altered values zb1_save,var1_save
! 28.06.2016    ggu     bug fix: handel situation with dh == 0
!
!*****************************************************************
!------------------------------------------------------------------
        module vert_intp
!------------------------------------------------------------------
        contains
!------------------------------------------------------------------

	subroutine intp_vert(bcons,nl1,zb1,var1,nl2,zb2,var2)

! vertical interpolation of variables from one grid to another

	implicit none

	logical bcons		!conserve total quantity
	integer nl1		!number of levels of first grid
	double precision zb1(0:nl1+1)	!depth of bottom of vertical boxes
	double precision var1(nl1+1)	!value of variable at center of box
	integer nl2		!number of levels of second grid
	double precision zb2(0:nl2)		!depth of bottom of vertical boxes
	double precision var2(nl2)		!value of variable at center of box (out)

! values are interpolated from first grid to second grid
! z levels refer to bottom of each grid (zb1(nl1) is total depth of column)
! zb(0) is surface ... normally 0
! variables are considered at center of box
! values are considered to be constant for every box
! bcons controlls if total quantity is conserved or not
!
! ATTENTION: arrays zb1, var1 MUST be dimensioned one element bigger
!	then the available data. The last element is accessed and altered 
!	by this subroutine, but at return is restored to its original value
!
! output is var2, all other variables are input values

	logical bmiss,bdebug
	integer l,j,ltop1,lbot1
	double precision ztop2,ztop1,zbot2,zbot1
	double precision ztop,zbot
	double precision vint1,vint2,fact
	double precision val
	double precision zb1_save,var1_save
	double precision dz,dh

!---------------------------------------------------------
! handle special situations
!---------------------------------------------------------

	if( nl1 < 1 .or. nl2 < 1 ) then
	  write(6,*) 'nl1,nl2: ',nl1,nl2
	  stop 'error stop intp_vert: nl1,nl2 < 1'
	else if( nl1 == 1 ) then	!only one layer in input
	  var2 = var1(1)
	  return
	end if

!---------------------------------------------------------
! initialize variables
!---------------------------------------------------------

	zb1_save = zb1(nl1+1)
	var1_save = var1(nl1+1)

	zb1(nl1+1) = zb2(nl2)
	var1(nl1+1) = var1(nl1)

	ltop1 = 0
	vint2 = 0.	!total content in second array

	bdebug = .false.

!---------------------------------------------------------
! loop over second array and interpolate onto it
!---------------------------------------------------------

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

!---------------------------------------------------------
! reset modified values
!---------------------------------------------------------

	zb1(nl1+1) = zb1_save
	var1(nl1+1) = var1_save

!---------------------------------------------------------
! if bcons we have to adapt the interpolated values to conserve total value
!---------------------------------------------------------

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

!---------------------------------------------------------
! end of routine
!---------------------------------------------------------

	return
   98	continue
	stop 'error stop intp_vert: missing value not possible'
   99	continue
	write(6,*) ztop1,ztop2
	write(6,*) zbot1,zbot2
	stop 'error stop intp_vert: interval in (1) must include (2)'
	end

!*****************************************************************

        subroutine intp_aver(n,horig,vorig,femval)

        implicit none

        integer n
        double precision horig(0:n)
        double precision vorig(n)
        double precision femval

        integer i
        double precision hacu,acu,h

        acu = 0.
        hacu = 0.
        do i=1,n
          h = horig(i) - horig(i-1)
          hacu = hacu + h
          acu = acu + h * vorig(i)
        end do

        femval = acu / hacu

        end

!*****************************************************************

	subroutine intp_vert_test

	implicit none

	integer ndim
	parameter (ndim=100)
	double precision zb1(0:ndim+1)
	double precision var1(ndim+1)
	double precision zb2(0:ndim)
	double precision var2(ndim)
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

!*****************************************************************

	!program main_intp_vert_test
	!call intp_vert_test
	!end

!*****************************************************************

!------------------------------------------------------------------
        end module vert_intp
!------------------------------------------------------------------
