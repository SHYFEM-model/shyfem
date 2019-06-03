
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

c topological set up routines
c
c contents :
c
c function kthis(i,ie)	        gets node at position i in ie
c function knext(k,ie) 	        gets node after k in ie
c function kbhnd(k,ie)          gets node before k in ie
c function ithis(k,ie) 	        gets i of k in ie
c function inext(k,ie) 	        gets i after k in ie
c function ibhnd(k,ie)	        gets i before k in ie
c
c subroutine link_fill(n)       returns filling of linkv
c
c subroutine get_elem_linkp(k,ipf,ipl)	gets pointer to elements around k
c subroutine get_node_linkp(k,ipf,ipl)	gets pointer to nodes around k
c subroutine get_elem_links(k,n,ibase)	gets pointer to elements around k
c subroutine get_node_links(k,n,ibase)	gets pointer to nodes around k
c
c subroutine get_elems_around(k,ndim,n,elems) returns all elems around node k
c subroutine get_nodes_around(k,ndim,n,nodes) returns all nodes around node k
c
c subroutine find_elems_to_segment(k1,k2,ie1,ie2) finds elements to segment
c
c notes :
c
c the preferred way to get nodes/elems around node k is through
c	get_elems_around()
c	get_nodes_around()
c 
c revision log :
c
c 01.08.2003	ggu	created from sublnk.f
c 16.03.2004	ggu	new routine node_links() and link_fill()
c 10.11.2007	ggu	new routine line_elems
c 28.08.2009	ggu	routine line_elems renamed to find_elems_to_segment
c 28.08.2009	ggu	new routines get_elems_around, get_nodes_around
c 09.09.2009	ggu	bug fix in find_elems_to_segment (BUGip2)
c 23.03.2010	ggu	changed v6.1.1
c 20.10.2011	ggu	check dimension in set_elem_links(), set_node_links()
c 16.12.2011	ggu	in lnk_elems at boundary set last value to 0
c 24.01.2012	ggu	changed VERS_6_1_41
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.05.2015	ggu	changed VERS_7_1_10
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 02.12.2015	ggu	lnk_elems and lnk_nodes eliminated
c 16.12.2015	ggu	changed VERS_7_3_16
c 28.04.2016	ggu	changed VERS_7_5_9
c 16.02.2019	ggu	changed VERS_7_5_60
c
c****************************************************************

        function kthis(i,ie)

c gets node at position i in ie
c
c i     position
c ie    element

	use basin

        implicit none

c arguments
	integer kthis
        integer i,ie

	kthis = nen3v(i,ie)

	end

c****************************************************************
c
        function knext(k,ie)
c
c gets node after k in ie
c
c k     actual node
c ie    element
c
	use basin

        implicit none
c
c arguments
	integer knext
        integer k,ie
c local
        integer i
c
        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            knext=nen3v(mod(i,3)+1,ie)
            return
          end if
        end do
c
        knext=0
c
        return
        end
c
c****************************************************************
c
        function kbhnd(k,ie)
c
c gets node before k in ie
c
c k     actual node
c ie    element
c
	use basin

        implicit none
c
c arguments
	integer kbhnd
        integer k,ie
c local
        integer i
c
        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            kbhnd=nen3v(mod(i+1,3)+1,ie)
            return
          end if
        end do
c
        kbhnd=0
c
        return
        end
c
c****************************************************************
c
        function ithis(k,ie)
c
c gets i of k in ie (might also use lenkiiv for this)
c
c k     actual node
c ie    element
c
	use basin

        implicit none
c
c arguments
	integer ithis
        integer k,ie
c local
        integer i
c
        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            ithis=i
            return
          end if
        end do
c
        ithis=0
c
        return
        end
c
c****************************************************************
c
        function inext(k,ie)
c
c gets i after k in ie
c
c k     actual node
c ie    element
c
	use basin

        implicit none
c
c arguments
	integer inext
        integer k,ie
c local
        integer i
c
        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            inext=mod(i,3)+1
            return
          end if
        end do
c
        inext=0
c
        return
        end
c
c****************************************************************
c
        function ibhnd(k,ie)
c
c gets i before k in ie
c
c k     actual node
c ie    element
c
	use basin

        implicit none
c
c arguments
	integer ibhnd
        integer k,ie
c local
        integer i
c
        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            ibhnd=mod(i+1,3)+1
            return
          end if
        end do
c
        ibhnd=0
c
        return
        end

c****************************************************************
c****************************************************************
c****************************************************************

        subroutine link_fill(n)

c returns filling of linkv

	use mod_geom
	use basin, only : nkn,nel,ngr,mbw

        implicit none

c arguments
        integer n       !filling of linkv (return)

        n = ilinkv(nkn+1)

        end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine get_elem_linkp(k,ipf,ipl)

c gets pointer to first and last element around k
c
c to loop over the neibor elements, use similar:
c
c       call get_elem_linkp(k,ipf,ipl)
c       do ip=ipf,ipl
c         ien = lenkv(ip)          !ien is number of neibor element
c       end do

	use mod_geom

	implicit none

	integer k,ipf,ipl

        ipf = ilinkv(k)+1
        ipl = ilinkv(k+1)

	if( lenkv(ipl) .eq. 0 ) ipl = ipl - 1	!FIXME

	end

c****************************************************************

	subroutine get_node_linkp(k,ipf,ipl)

c gets pointer to first and last node around k
c
c to loop over the neibor nodes, use similar:
c
c       call get_node_linkp(k,ipf,ipl)
c       do ip=ipf,ipl
c         kn = linkv(ip)          !kn is number of neibor node
c       end do

	use mod_geom

	implicit none

	integer k,ipf,ipl

        ipf = ilinkv(k)+1
        ipl = ilinkv(k+1)

	end

c****************************************************************

	subroutine get_elem_links(k,n,ibase)

c gets pointer and total number of elements around k
c
c to loop over the neibor elements, use similar:
c
c       call get_elem_links(k,n,ibase)
c       do i=1,n
c         ien = lenkv(ibase+i)          !ien is number of neibor element
c       end do

	use mod_geom

	implicit none

	integer k,n,ibase

	n = ilinkv(k+1)-ilinkv(k)
	ibase = ilinkv(k)

	if( lenkv(ibase+n) .eq. 0 ) n = n - 1

        end

c****************************************************************

	subroutine get_node_links(k,n,ibase)

c gets pointer and total number of nodes around k
c
c to loop over the neibor nodes, use similar:
c
c       call get_node_links(k,n,ibase)
c       do i=1,n
c         kn = linkv(ibase+i)          !kn is number of neibor node
c       end do

	use mod_geom

	implicit none

	integer k,n,ibase

	n = ilinkv(k+1)-ilinkv(k)
	ibase = ilinkv(k)

        end

c****************************************************************

        subroutine get_elems_around(k,ndim,n,elems)

c returns all elems around node k

	use mod_geom

        implicit none

        integer k               !central node
        integer ndim            !dimension of elems()
        integer n               !total number of elems around k (return)
        integer elems(ndim)     !elems around k (return)

	integer i,ibase

	n = ilinkv(k+1)-ilinkv(k)
	ibase = ilinkv(k)

	if( lenkv(ibase+n) .eq. 0 ) n = n - 1

	do i=1,n
	  elems(i) = lenkv(ibase+i)
	end do

	end

c****************************************************************

        subroutine get_nodes_around(k,ndim,n,nodes)

c returns all nodes around node k

	use mod_geom

        implicit none

        integer k               !central node
        integer ndim            !dimension of nodes()
        integer n               !total number of nodes around k (return)
        integer nodes(ndim)     !nodes around k (return)

	integer i,ibase

	n = ilinkv(k+1)-ilinkv(k)
	ibase = ilinkv(k)

	do i=1,n
	  nodes(i) = linkv(ibase+i)
	end do

	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine find_elems_to_segment(k1,k2,ie1,ie2)

c finds elements to segment between nodes k1 and k2
c
c returns elements in ie1 and ie2
c
c ie1 is to the left of segment k1-k2, ie2 to the right
c if boundary segment only one ie is set, the other is zero
c if no such segment, both ie are zero

	use mod_geom

	implicit none

c arguments
        integer k1,k2,ie1,ie2

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

c****************************************************************
c****************************************************************
c****************************************************************
c obsolete routines
c****************************************************************
c****************************************************************
c****************************************************************

        subroutine pntfla(k,ipf,ipl)

c gets pointer to first and last element in lenkv
c
c superseeded by get_elem_linkp() - do not use anymore
c
c k     actual node
c ipf   first element (return)
c ipl   last  element (return)

	use mod_geom

        implicit none

c arguments
        integer k,ipf,ipl

        ipf=ilinkv(k)+1
        ipl=ilinkv(k+1)

	if( lenkv(ipl) .eq. 0 ) ipl = ipl - 1	!FIXME

        end

c****************************************************************

        subroutine node_links(ie,ip)

c gets pointer to linkv for element ie
c
c attention - this is really CPU intensive

	use mod_geom
	use basin

        implicit none

c arguments
        integer ie              !element
        integer ip(3,3)         !pointer into linkv

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

c****************************************************************

