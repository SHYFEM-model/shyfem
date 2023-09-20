
!==================================================================
        module mod_color
!==================================================================

!==================================================================
        contains
!==================================================================

	subroutine make_node_matrix(nel,nen3v,nkn,ids,nmax,matrix)

	implicit none

	integer nel			!total number of elements
	integer nen3v(3,nel)		!element index
	integer nkn			!total number of nodes
	integer ids(nkn)		!id (domain) of nodes
	integer nmax
	integer, allocatable :: matrix(:,:)

	integer ie,ii1,ii2,k1,k2,id1,id2

	nmax = maxval(ids)
	allocate( matrix(0:nmax,0:nmax) )
	matrix = 0

	do ie=1,nel
	  do ii1=1,3
	    k1 = nen3v(ii1,ie)
	    ii2 = 1 + mod(ii1,3)
	    k2 = nen3v(ii2,ie)
	    id1 = ids(k1)
	    id2 = ids(k2)
	    if( id1 /= id2 ) then
	      matrix(id1,id2) = matrix(id1,id2) + 1
	      matrix(id2,id1) = matrix(id2,id1) + 1
	    end if
	  end do
	end do
	
	end

!*****************************************************************

	subroutine make_elem_matrix(nel,ecv,ids,nmax,matrix)

	implicit none

	integer nel			!total number of elements
	integer ecv(3,nel)		!element index
	integer ids(nel)		!id (domain) of nodes
	integer nmax
	integer, allocatable :: matrix(:,:)

	integer ie1,ie2,ii,id1,id2

	nmax = maxval(ids)
	allocate( matrix(0:nmax,0:nmax) )
	matrix = 0

	do ie1=1,nel
	  do ii=1,3
	    ie2 = ecv(ii,ie1)
	    id1 = ids(ie1)
	    id2 = ids(ie2)
	    if( id1 /= id2 ) then
	      !matrix(id1,id2) = matrix(id1,id2) + 1
	      !matrix(id2,id1) = matrix(id2,id1) + 1
	      matrix(id1,id2) = 1
	      matrix(id2,id1) = 1
	    end if
	  end do
	end do
	
	end

!*****************************************************************

	subroutine release_color_matrix(matrix)

	integer, allocatable :: matrix(:,:)

	if( allocated(matrix) ) deallocate(matrix)

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine color_graph(n,ids,nmax,matrix,ncol,color)

! from partitioning on nodes create colors
!
! adjacent domains will have different colors
!
! ids can be starting from 0 or 1 up to nmax
! color will be 1 to ncol where ncol is maximum color needed

	implicit none

	integer n			!total number of nodes
	integer ids(n)			!id (domain) of nodes
	integer nmax			!maximum id in whole domain
	integer matrix(0:nmax,0:nmax)	!graph connection
	integer ncol			!total number of colors needed (return)
	integer color(n)		!colors for nodes (return)

	integer nmin
	integer ie,ii1,ii2,k1,k2,id1,id2
	integer id,nc,col,k,i
	integer, allocatable :: colors(:)	!id to color table

	nmin = minval(ids)

	allocate(colors(0:nmax))
	colors = 0
	ncol = 1
	colors(nmin) = ncol

	do id=nmin+1,nmax
	  call get_color(id,col)
	  colors(id) = col
	  if( col > ncol ) ncol = col
	end do

	do i=1,n
	  id = ids(i)
	  color(i) = colors(id) 
	end do

	contains

	subroutine get_color(id1,col)
	integer id1,col
	integer cols(1:nmax)
	integer id2,c
	cols = 0
	do id2=nmin,nmax
	  if(matrix(id1,id2)>0) then
	    c = colors(id2)
	    if( c > 0 ) cols(c) = 1
	  end if
	end do
	do c=1,nmax
	  if( cols(c) == 0 ) then
	    col = c
	    return
	  end if
	end do
	stop 'error stop get_color: no color found'
	end
	    
	end

!*****************************************************************
!*****************************************************************
!*****************************************************************
! old routines
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine color_nodes(nel,nen3v,nkn,ids,ncol,color)

! from partitioning on nodes create colors
!
! adjacent domains will have different colors
!
! ids can be starting from 0 or 1 up to nmax
! color will be 1 to ncol where ncol is maximum color needed

	implicit none

	integer nel			!total number of elements
	integer nen3v(3,nel)		!element index
	integer nkn			!total number of nodes
	integer ids(nkn)		!id (domain) of nodes
	integer ncol			!total number of colors needed (return)
	integer color(nkn)		!colors for nodes (return)

	integer nmin,nmax
	integer ie,ii1,ii2,k1,k2,id1,id2
	integer id,nc,col,k
	integer, allocatable :: list(:,:)	!matrix of domain connections
	integer, allocatable :: colors(:)	!id to color table

	nmin = minval(ids)
	nmax = maxval(ids)

	allocate( list(0:nmax,0:nmax) )
	list = 0

	do ie=1,nel
	  do ii1=1,3
	    k1 = nen3v(ii1,ie)
	    ii2 = 1 + mod(ii1,3)
	    k2 = nen3v(ii2,ie)
	    id1 = ids(k1)
	    id2 = ids(k2)
	    if( id1 /= id2 ) then
	      list(id1,id2) = list(id1,id2) + 1
	      list(id2,id1) = list(id2,id1) + 1
	    end if
	  end do
	end do
	
	allocate(colors(0:nmax))
	colors = 0
	ncol = 1
	colors(nmin) = ncol

	do id=nmin+1,nmax
	  call get_color(id,col)
	  colors(id) = col
	  if( col > ncol ) ncol = col
	end do

	do k=1,nkn
	  id = ids(k)
	  color(k) = colors(id) 
	end do

	contains

	subroutine get_color(id1,col)
	integer id1,col
	integer cols(1:nmax)
	integer id2,c
	cols = 0
	do id2=nmin,nmax
	  if(list(id1,id2)>0) then
	    c = colors(id2)
	    if( c > 0 ) cols(c) = 1
	  end if
	end do
	do c=1,nmax
	  if( cols(c) == 0 ) then
	    col = c
	    return
	  end if
	end do
	stop 'error stop get_color: no color found'
	end
	    
	end

!*****************************************************************

	subroutine color_elems(nel,nen3v,ecv,ids,ncol,color)

! from partitioning on elems create colors
!
! adjacent domains will have different colors
!
! ids can be starting from 0 or 1 up to nmax
! color will be 1 to ncol where ncol is maximum color needed

	implicit none

	integer nel			!total number of elements
	integer nen3v(3,nel)		!element index
	integer ecv(3,nel)		!element neighbors
	integer ids(nel)		!id (domain) of elements
	integer ncol			!total number of colors needed (return)
	integer color(nel)		!colors for elements (return)

	integer nmin,nmax
	integer ie,ii1,ii2,k1,k2,id1,id2
	integer id,nc,col,k,ie1,ie2,ii
	integer, allocatable :: list(:,:)
	integer, allocatable :: colors(:)

	nmin = minval(ids)
	nmax = maxval(ids)

	allocate( list(0:nmax,0:nmax) )
	list = 0

	do ie1=1,nel
	  do ii=1,3
	    ie2 = ecv(ii,ie)
	    id1 = ids(ie1)
	    id2 = ids(ie2)
	    if( id1 /= id2 ) then
	      list(id1,id2) = list(id1,id2) + 1
	      list(id2,id1) = list(id2,id1) + 1
	    end if
	  end do
	end do
	
	allocate(colors(0:nmax))
	colors = 0
	ncol = 1
	colors(nmin) = ncol

	do id=nmin+1,nmax
	  call get_color(id,col)
	  colors(id) = col
	  if( col > ncol ) ncol = col
	end do

	do ie=1,nel
	  id = ids(ie)
	  color(id) = colors(id) 
	end do

	contains

	subroutine get_color(id1,col)
	integer id1,col
	integer cols(1:nmax)
	integer id2,c
	cols = 0
	do id2=nmin,nmax
	  if(list(id1,id2)>0) then
	    c = colors(id2)
	    if( c > 0 ) cols(c) = 1
	  end if
	end do
	do c=1,nmax
	  if( cols(c) == 0 ) then
	    col = c
	    return
	  end if
	end do
	stop 'error stop get_color: no color found'
	end
	    
	end

!*****************************************************************

	subroutine make_translation_table(n,ids,color,nmax,tt)

	implicit none

	integer n
	integer ids(n)
	integer color(n)
	integer nmax
	integer tt(0:nmax)

	integer i,id,c

	tt = -1

	do i=1,n
	  id = ids(i)
	  c = color(i)
	  if( tt(id) < 0 ) tt(id) = c
	  if( tt(id) /= c ) then
	    write(6,*) i,id,c
	    stop 'error stop make_translation_table: inconsistency'
	  end if
	end do

	end

!*****************************************************************

!==================================================================
        end module mod_color
!==================================================================

