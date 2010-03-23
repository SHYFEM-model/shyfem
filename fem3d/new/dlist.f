
c*********************************************************************

	subroutine dlist_init(n,ix,bcirc)

	implicit none

	integer n
	integer ix(2,1)
	logical bcirc

	integer i,is,ie

	if( n .eq. 0 ) return

	do i=2,n-1
	  ix(1,i) = i-1
	  ix(2,i) = i+1
	end do

	if( bcirc ) then
	  is = n
	  ie = 1
	else
	  is = 0
	  ie = 0
	end if

	if( n .eq. 1 ) then
	  ix(1,1) = is
	  ix(2,1) = ie
	else
	  ix(1,1) = is
	  ix(2,1) = 2
	  ix(1,n) = n-1
	  ix(2,n) = ie
	end if

	end

c*********************************************************************

	subroutine dlist_delete(n,ix,i)

	implicit none

	integer n
	integer ix(2,1)
	integer i

	integer ia,ib

	if( n .eq. 0 ) return

	if( n .gt. 1 ) then
	  ib = ix(1,i)
	  ia = ix(2,i)
	  if( ib .gt. 0 ) ix(2,ib) = ia
	  if( ia .gt. 0 ) ix(1,ia) = ib
	end if

	ix(1,i) = 0
	ix(2,i) = 0

	end

c*********************************************************************

