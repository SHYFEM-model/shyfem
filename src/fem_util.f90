!
! $Id: subnsu.f,v 1.9 2009-02-04 15:26:54 georg Exp $
!
! utility routines for fem model
!
! contents :
!
! subroutine baric(ie,x,y)		finds baricentre of element
! subroutine coord(k,x,y)		returns coordinates of node k (internal)
! function area_element(ie)		area for element ie (internal)
!
! function ipext(k)             returns extern node number
! function ieext(k)             returns extern element number
! function ipint(k)             returns intern node number
! function ieint(k)             returns intern element number
!
! subroutine n2int(n,nnodes,berror)	converts external nodes to internal
! subroutine e2int(n,nelems,berror)	converts external elements to internal
!
! revision log :
!
! 02.06.1997	ggu	eliminated extint,extinw,exinel (not needed)
! 24.02.1999	ggu	new subroutine n2int
! 19.11.1999	ggu	new subroutine e2int
! 27.01.2009	ggu	new subroutine coord
! 07.07.2011	ggu	new subroutine area_element()
!
!*******************************************************
!-------------------------------------------------------
        module fem_util
!-------------------------------------------------------
        contains
!-------------------------------------------------------

	subroutine baric(ie,x,y)

! finds baricentre of element
!
! ie		number of element
! x,y		coordinates of baricentre (return value)

	use basin

	implicit none

! arguments
	integer ie
	double precision x,y
! common blocks
! local variables
	integer i,kkk
	double precision xb,yb

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

!***************************************************************

        function area_element(ie)

! area for element ie

	use basin

        implicit none

! arguments
        double precision area_element
        integer ie
! local
        integer kn1,kn2,kn3
        double precision x1,x2,x3,y1,y2,y3
        double precision half

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

!***************************************************************

	subroutine coord(k,x,y)

! returns coordinates of node k (internal)

	use basin

	implicit none

	integer k
	double precision x,y

	x = xgv(k)
	y = ygv(k)

	end

!***************************************************************
!
	function ipext(k)
!
! returns external node number
!
! k     internal node number
! ipext external node number, 0 if error
!
	use basin

	implicit none
	integer ipext
	integer k

	if(k.lt.1.or.k.gt.nkn) then
		ipext=0
	else
		ipext=ipv(k)
	end if
!
	return
	end
!
!***************************************************************
!
	function ieext(k)
!
! returns external element number
!
! k     internal element number
! ieext external element number, 0 if error
!
	use basin

	implicit none
	integer ieext
	integer k
!
	if(k.lt.1.or.k.gt.nel) then
		ieext=0
	else
		ieext=ipev(k)
	end if
!
	return
	end
!
!***************************************************************
!
	function ipint(k)
!
! returns internal node number
!
! k     external node number
! ipint internal node number, 0 if error
!
	use basin

	implicit none
	integer ipint
	integer k,i
!
	do i=1,nkn
	if(ipv(i).eq.k) goto 1
	end do
	i=0
    1   continue
	ipint=i
!
	return
	end
!
!***************************************************************
!
	function ieint(k)
!
! returns internal element number
!
! k     external element number
! ieint internal element number, 0 if error
!
	use basin

	implicit none
	integer ieint
	integer k,i
!
	do i=1,nel
	if(ipev(i).eq.k) goto 1
	end do
	i=0
    1   continue
	ieint=i
!
	return
	end

!***************************************************************

	subroutine n2int(n,nnodes,berror)

! converts external nodes numbers to internal ones

	implicit none

	integer n		!total number of nodes
	integer nnodes(n)	!node numbers
	logical berror		!true on return if error

	integer i,ke,ki

	berror = .false.

	do i=1,n
	  ke = nnodes(i)
	  ki = ipint(ke)
	  if( ki .le. 0 .and. ke .gt. 0 ) then
	    write(6,*) 'No such node : ',ke
	    berror = .true.
	  end if
	  nnodes(i) = ki
	end do

	end

!***************************************************************

	subroutine e2int(n,nelems,berror)

! converts external element numbers to internal ones

	implicit none

	integer n		!total number of elements
	integer nelems(n)	!element numbers
	logical berror		!true on return if error

	integer i,ke,ki

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

!***************************************************************

!-------------------------------------------------------
        end module fem_util
!-------------------------------------------------------
