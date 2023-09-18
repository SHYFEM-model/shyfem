!
! $Id: sublnku.f,v 1.7 2009-09-14 08:20:58 georg Exp $
!
! topological set up routines
!
! contents :
!
! function kthis(i,ie)	        gets node at position i in ie
! function knext(k,ie) 	        gets node after k in ie
! function kbhnd(k,ie)          gets node before k in ie
! function ithis(k,ie) 	        gets i of k in ie
! function inext(k,ie) 	        gets i after k in ie
! function ibhnd(k,ie)	        gets i before k in ie
!
! subroutine link_fill(n)       returns filling of linkv
!
! subroutine get_elem_linkp(k,ipf,ipl)	gets pointer to elements around k
! subroutine get_node_linkp(k,ipf,ipl)	gets pointer to nodes around k
! subroutine get_elem_links(k,n,ibase)	gets pointer to elements around k
! subroutine get_node_links(k,n,ibase)	gets pointer to nodes around k
!
! subroutine get_elems_around(k,ndim,n,elems) returns all elems around node k
! subroutine get_nodes_around(k,ndim,n,nodes) returns all nodes around node k
!
! subroutine find_elems_to_segment(k1,k2,ie1,ie2) finds elements to segment
!
! notes :
!
! the preferred way to get nodes/elems around node k is through
!	get_elems_around()
!	get_nodes_around()
! 
! revision log :
!
! 01.08.2003	ggu	created from sublnk.f
! 16.03.2004	ggu	new routine node_links() and link_fill()
! 10.11.2007	ggu	new routine line_elems
! 28.08.2009	ggu	routine line_elems renamed to find_elems_to_segment
! 28.08.2009	ggu	new routines get_elems_around, get_nodes_around
! 09.09.2009	ggu	bug fix in find_elems_to_segment (BUGip2)
! 20.10.2011	ggu	check dimension in set_elem_links(), set_node_links()
! 16.12.2011	ggu	in lnk_elems at boundary set last value to 0
! 02.12.2015	ggu	lnk_elems and lnk_nodes eliminated
!
!****************************************************************

      module lnku


        contains

      integer function kthis(i,ie)

! gets node at position i in ie
!
! i     position
! ie    element

	use basin

        implicit none

! arguments
	!integer kthis
        !integer :: kthis
        integer i,ie
! common
	include 'param.h'

	kthis = nen3v(i,ie)

	end

!****************************************************************

       function knext(k,ie)

! gets node after k in ie
!
! k     actual node
! ie    element

	use basin

        implicit none
!
! arguments
	integer knext
        integer k,ie
! common
	include 'param.h'
! local
        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            knext=nen3v(mod(i,3)+1,ie)
            return
          end if
        end do

        knext=0

        return
        end

!****************************************************************

        function kbhnd(k,ie)

! gets node before k in ie

! k     actual node
! ie    element

	use basin

        implicit none

! arguments
	integer kbhnd
        integer k,ie
! common
	include 'param.h'
! local
        integer i
!
        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            kbhnd=nen3v(mod(i+1,3)+1,ie)
            return
          end if
        end do

        kbhnd=0

        return
        end

!****************************************************************

        function ithis(k,ie)

! gets i of k in ie (might also use lenkiiv for this)

! k     actual node
! ie    element

	use basin

        implicit none

! arguments
	integer ithis
        integer k,ie
! common
	include 'param.h'
! local
        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            ithis=i
            return
          end if
        end do

        ithis=0

        return
        end

!****************************************************************

        function inext(k,ie)

! gets i after k in ie

! k     actual node
! ie    element

	use basin

        implicit none

! arguments
	integer inext
        integer k,ie
! common
	include 'param.h'
! local
        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            inext=mod(i,3)+1
            return
          end if
        end do

        inext=0

        return
        end

!****************************************************************

        function ibhnd(k,ie)

! gets i before k in ie

! k     actual node
! ie    element

	use basin

        implicit none

! arguments
	integer ibhnd
        integer k,ie
! common
	include 'param.h'
! local
        integer i

        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            ibhnd=mod(i+1,3)+1
            return
          end if
        end do

        ibhnd=0

        return
        end

!****************************************************************

        subroutine link_fill(n)

! returns filling of linkv

	use geom
	use basin, only : nkn,nel,ngr,mbw

        implicit none

! arguments
        integer n       !filling of linkv (return)
! common
	include 'param.h'

        n = ilinkv(nkn+1)

        end

!****************************************************************

	subroutine get_elem_linkp(k,ipf,ipl)

! gets pointer to first and last element around k
!
! to loop over the neibor elements, use similar:
!
!       call get_elem_linkp(k,ipf,ipl)
!       do ip=ipf,ipl
!         ien = lenkv(ip)          !ien is number of neibor element
!       end do

	use geom

	implicit none

	integer k,ipf,ipl

	include 'param.h'

        ipf = ilinkv(k)+1
        ipl = ilinkv(k+1)

	if( lenkv(ipl) .eq. 0 ) ipl = ipl - 1	!FIXME

	end

!****************************************************************

	subroutine get_node_linkp(k,ipf,ipl)

! gets pointer to first and last node around k
!
! to loop over the neibor nodes, use similar:
!
!       call get_node_linkp(k,ipf,ipl)
!       do ip=ipf,ipl
!         kn = linkv(ip)          !kn is number of neibor node
!       end do

	use geom

	implicit none

	integer k,ipf,ipl

	include 'param.h'

        ipf = ilinkv(k)+1
        ipl = ilinkv(k+1)

	end

!****************************************************************

	subroutine get_elem_links(k,n,ibase)

! gets pointer and total number of elements around k
!
! to loop over the neibor elements, use similar:
!
!       call get_elem_links(k,n,ibase)
!       do i=1,n
!         ien = lenkv(ibase+i)          !ien is number of neibor element
!       end do

	use geom

	implicit none

	integer k,n,ibase

	include 'param.h'

	n = ilinkv(k+1)-ilinkv(k)
	ibase = ilinkv(k)

	if( lenkv(ibase+n) .eq. 0 ) n = n - 1

        end

!****************************************************************

	subroutine get_node_links(k,n,ibase)

! gets pointer and total number of nodes around k
!
! to loop over the neibor nodes, use similar:
!
!       call get_node_links(k,n,ibase)
!       do i=1,n
!         kn = linkv(ibase+i)          !kn is number of neibor node
!       end do

	use geom

	implicit none

	integer k,n,ibase

	include 'param.h'

	n = ilinkv(k+1)-ilinkv(k)
	ibase = ilinkv(k)

        end

!****************************************************************

        subroutine get_elems_around(k,ndim,n,elems)

! returns all elems around node k

	use geom

        implicit none

        integer k               !central node
        integer ndim            !dimension of elems()
        integer n               !total number of elems around k (return)
        integer elems(ndim)     !elems around k (return)

	integer i,ibase

	include 'param.h'

	n = ilinkv(k+1)-ilinkv(k)
	ibase = ilinkv(k)

	if( lenkv(ibase+n) .eq. 0 ) n = n - 1

	do i=1,n
	  elems(i) = lenkv(ibase+i)
	end do

	end

!****************************************************************

        subroutine get_nodes_around(k,ndim,n,nodes)

! returns all nodes around node k

	use geom

        implicit none

        integer k               !central node
        integer ndim            !dimension of nodes()
        integer n               !total number of nodes around k (return)
        integer nodes(ndim)     !nodes around k (return)

	integer i,ibase

	include 'param.h'

	n = ilinkv(k+1)-ilinkv(k)
	ibase = ilinkv(k)

	do i=1,n
	  nodes(i) = linkv(ibase+i)
	end do

	end

!****************************************************************

	subroutine find_elems_to_segment(k1,k2,ie1,ie2)

! finds elements to segment between nodes k1 and k2
!
! returns elements in ie1 and ie2
!
! ie1 is to the left of segment k1-k2, ie2 to the right
! if boundary segment only one ie is set, the other is zero
! if no such segment, both ie are zero

	use geom

	implicit none

! arguments
        integer k1,k2,ie1,ie2
! common
	include 'param.h'

        integer k,ipf,ipl,ip,ip2

	k = k1
        ipf=ilinkv(k)+1
        ipl=ilinkv(k+1)

	ie1 = 0
	ie2 = 0

	do ip=ipf,ipl
	  k = linkv(ip)
	  if( k .eq. k2 ) then
	    ie1 = lenkv(ip)
	    if( ip .eq. ipf ) then
		ip2 = ipl		!this sets it to 0
	    else
		ip2 = ip - 1		!previous element	!BUGip2
	    end if
	    ie2 = lenkv(ip2)
	    return
	  end if
	end do

	end

!****************************************************************
! obsolete routines
!****************************************************************

        subroutine pntfla(k,ipf,ipl)

! gets pointer to first and last element in lenkv
!
! superseeded by get_elem_linkp() - do not use anymore
!
! k     actual node
! ipf   first element (return)
! ipl   last  element (return)

	use geom

        implicit none

! arguments
        integer k,ipf,ipl
! common
	include 'param.h'

        ipf=ilinkv(k)+1
        ipl=ilinkv(k+1)

	if( lenkv(ipl) .eq. 0 ) ipl = ipl - 1	!FIXME

        end

!****************************************************************

        subroutine node_links(ie,ip)

! gets pointer to linkv for element ie
!
! attention - this is really CPU intensive

	use geom
	use basin

        implicit none

! arguments
        integer ie              !element
        integer ip(3,3)         !pointer into linkv
! common
	include 'param.h'

        integer ii,iii,k,kn,i
        integer ipf,ipl

        do ii=1,3
          k = nen3v(ii,ie)
          ipf=ilinkv(k)+1
          ipl=ilinkv(k+1)
          do iii=1,3
            kn = nen3v(iii,ie)
            if( k .eq. kn ) then
              ip(ii,iii) = 0
            else
              do i=ipf,ipl
                if( linkv(i) .eq. kn ) goto 1
              end do
              goto 99
    1         continue
              ip(ii,iii) = i
            end if
          end do
        end do

        return
   99   continue
        write(6,*) ie,ii,iii,k,kn,ipf,ipl
        write(6,*) (linkv(i),i=ipf,ipl)
        stop 'error stop node_links: internal error (1)'
        end

  
      end module lnku
