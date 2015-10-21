
!***************************************************************

	subroutine grd_to_basin

	use grd
	use basin

	implicit none

	integer nk,ne,nl,nne,nnl
	integer k,ie,ii,ibase,n,i
	integer nhn,nhe
	real flag
	logical bdebug

!--------------------------------------------------------
! initialize
!--------------------------------------------------------

	flag = -999.
	bdebug = .false.

	call grd_get_params(nk,ne,nl,nne,nnl)
	call basin_init(nk,ne)

!--------------------------------------------------------
! copy params
!--------------------------------------------------------

	dcorbas = dcor_grd
	dirnbas = dirn_grd
	descrr = title_grd

!--------------------------------------------------------
! copy arrays
!--------------------------------------------------------

	xgv = xv
	ygv = yv
	ipev = ippev
	ipv = ippnv
	iarnv = ianv
	iarv = iaev

!--------------------------------------------------------
! copy element index
!--------------------------------------------------------

	do ie=1,ne
	  ibase = ipntev(ie-1)
	  n = ipntev(ie) - ipntev(ie-1)
	  if( n .ne. 3 ) then
	    !write(6,*) (ipntev(i),i=0,10)
	    write(6,*) 'element is not triangle: ',ie,ippev(ie)
	    stop 'error stop grd_to_basin: not a triangle'
	  end if
	  do ii=1,3
	    nen3v(ii,ie) = inodev(ibase+ii)
	  end do
	end do

!--------------------------------------------------------
! check depth and copy
!--------------------------------------------------------

	nhe = 0
	do ie=1,ne
	  if( hhev(ie) .ne. flag ) nhe = nhe + 1
	end do

	nhn = 0
	do k=1,nk
	  if( hhnv(k) .ne. flag ) nhn = nhn + 1
	end do

	if( nhe > 0 .and. nhn > 0 ) then
	  write(6,*) 'nhe,nhn: ',nhe,nhn
	  if( nhe == nel .and. nhn ==nkn ) then
	    write(6,*) 'can use both depths...'
	    write(6,*) '... using element values'
	    nhn = 0
	  else if( nhe == nel ) then
	    write(6,*) 'depth on nodes incomplete...'
	    write(6,*) '... using element values'
	    nhn = 0
	  else if( nhn == nkn ) then
	    write(6,*) 'depth on elements incomplete...'
	    write(6,*) '... using nodal values'
	    nhe = 0
	  else
	    write(6,*) 'depth values are incomplete...'
	    stop 'error stop grd_to_basin: depth incomplete'
	  end if
	end if

	write(6,*) 'nhe,nhn: ',nhe,nhn

	if( nhn == 0 ) then
	  do ie=1,ne
	    do ii=1,3
	      hm3v(ii,ie) = hhev(ie)
	    end do
	  end do
	else
	  do ie=1,ne
	    do ii=1,3
	      k = nen3v(ii,ie)
	      hm3v(ii,ie) = hhnv(k)
	    end do
	  end do
	end if

!--------------------------------------------------------
! debug output
!--------------------------------------------------------

	if( bdebug ) then

	do ie=1,ne,ne/10
	  write(6,*) ie,(nen3v(ii,ie),ii=1,3)
	end do

	do k=1,nk,nk/10
	  write(6,*) k,ippnv(k),xv(k),yv(k)
	end do

	end if

!--------------------------------------------------------
! end of routine
!--------------------------------------------------------

	end

!***************************************************************

	subroutine basin_to_grd

	use grd
	use basin

	implicit none

	integer nk,ne,nl,nne,nnl
	integer k,ie,ii,ibase,n
	integer nhn,nhe
	real flag
	logical bdebug
	logical bconst,bunique
	double precision h,he,hm

!--------------------------------------------------------
! initialize
!--------------------------------------------------------

	flag = -999.
	bdebug = .false.

	nl = 0
	nne = 3*nel
	nnl = 0

	write(6,*) 'basin_to_grd: ',nkn,nel
	call grd_init(nkn,nel,nl,nne,nnl)

!--------------------------------------------------------
! copy params
!--------------------------------------------------------

	dcor_grd = dcorbas
	dirn_grd = dirnbas
	title_grd = descrr

!--------------------------------------------------------
! copy arrays
!--------------------------------------------------------

	xv(1:nkn) = xgv(1:nkn)
	yv(1:nkn) = ygv(1:nkn)
	ippev(1:nel) = ipev(1:nel)
	ippnv(1:nkn) = ipv(1:nkn)
	ianv(1:nkn) = iarnv(1:nkn)
	iaev(1:nel) = iarv(1:nel)

!--------------------------------------------------------
! copy element index
!--------------------------------------------------------

	ibase = 0
	ipntev(0) = 0
	do ie=1,nel
	  ibase = ipntev(ie-1)
	  do ii=1,3
	    ibase = ibase + 1
	    inodev(ibase) = nen3v(ii,ie)
	  end do
	  ipntev(ie) = ibase
	end do

!--------------------------------------------------------
! check depth and copy
!--------------------------------------------------------

	hhnv = flag
	hhev = flag

	bunique = .true.
	bconst = .true.

	do ie=1,nel
	  he = hm3v(1,ie)
	  hm = 0.
	  do ii=1,3
	    h = hm3v(ii,ie)
	    k = nen3v(ii,ie)
	    if( hhnv(k) == flag ) hhnv(k) = h
	    if( hhnv(k) /= h ) bunique = .false.
	    if( he /= h ) bconst = .false.
	    hm = hm + h
	  end do
	  hm = hm / 3.
	  hhev(ie) = hm
	end do

	if( bconst .and. bunique ) then		!constant depth
	  hhev = hm3v(1,1)
	  hhnv = flag
	else if( bunique ) then			!uniqie depth at nodes
	  hhev = flag
	else
	  hhnv = flag
	end if
	  
!--------------------------------------------------------
! end of routine
!--------------------------------------------------------

	end

!***************************************************************

