c
c $Id: sublnka.f,v 1.7 2009-05-21 09:24:00 georg Exp $
c
c topological set up routines
c
c contents :
c
c subroutine setnod
c			 sets array inodv
c
c revision log :
c
c 01.08.2003    ggu     created from sublnk.f
c 13.08.2003    ggu     in update_geom do not call setweg and setnod
c 06.11.2008    ggu     better error handling
c 06.04.2009    ggu     nlidim -> nlkdim
c 02.12.2011    ggu     print_bound_nodes() for writing out boundary nodes
c 04.05.2015    ggu     remove equivalence - use winkv as local array
c 20.05.2015    ggu     new call to mklenkii to set up lenkiiv
c
c*****************************************************************

	subroutine update_geom

c updates geometrical array (ieltv)

	use mod_geom
	use mod_geom_dynamic
	use basin, only : nkn,nel,ngr,mbw

        implicit none

c common
	include 'param.h'

c local
        integer n

c-------------------------------------------------------------
c update ieltv
c-------------------------------------------------------------

        call update_ielt(nel,inodv,ieltv)

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c*****************************************************************

        subroutine setnod

c sets (dynamic) array inodv
c
c inodv 
c	 0: internal node  
c	>0: open boundary node
c       -1: boundary node  
c	-2: out of system
c
c if open boundary node, inodv(k) is number of boundary (ggu 15.11.2001)

	use mod_geom_dynamic
	use evgeom
	use basin

        implicit none

	double precision, parameter :: winmax = 359.99
        integer ie,ii,k,n
        integer ibc,ibtyp
	integer nbc
        double precision winkv(nkn)

        integer ipext
	integer itybnd,nkbnds,kbnds,nbnds

c initialize array to hold angles

	winkv = 0.

c sum angles

        do ie=1,nel
          if(iwegv(ie).eq.0) then !element is in system
            do ii=1,3
              k=nen3v(ii,ie)
              winkv(k)=winkv(k)+ev(10+ii,ie)
            end do
          end if
	end do

c set up inodv

        do k=1,nkn
          if(winkv(k).gt.winmax) then     !internal node
            inodv(k)=0
          else if(winkv(k).eq.0.) then  !out of system
            inodv(k)=-2
          else                          !boundary node
            inodv(k)=-1
          end if
        end do

c now mark open boundary nodes

	nbc = nbnds()

	do ibc=1,nbc
          ibtyp = itybnd(ibc)
	  n = nkbnds(ibc)
          if(ibtyp.ge.3) then       !$$ibtyp3	!$$ibtyp4
            do ii=1,n
              k=kbnds(ibc,ii)
              if(inodv(k).eq.-1) inodv(k)=ibc
            end do
          else if(ibtyp.gt.0) then
            do ii=1,n
              k=kbnds(ibc,ii)
              if(inodv(k).ne.-1) goto 99
              inodv(k)=ibc
            end do
          end if
        end do

        return
   99   continue
        write(6,*) 'error for open boundary node'
	write(6,*) ibc,n,ii,k,ipext(k),inodv(k)
        if( inodv(k) .eq. 0 ) then
          write(6,*) 'the node is an inner node'
        else if( inodv(k) .eq. -2 ) then
          write(6,*) 'the node is not in the system (dry)'
        else
          write(6,*) 'the node has already been flagged as an'
          write(6,*) 'open boundary node (listed twice?)'
        end if
        write(6,*) 'external node number : ',ipext(k)
        write(6,*) 'boundary flag : ',inodv(k)
        stop 'error stop setnod : open boundary node'
        end

c****************************************************************
c****************************************************************
c****************************************************************
c****************************************************************

	function is_internal_node(k)

	use mod_geom_dynamic

	implicit none

	logical is_internal_node
	integer k

	include 'param.h'

	is_internal_node = inodv(k) .eq. 0

	end

c****************************************************************

	function is_boundary_node(k)

	use mod_geom_dynamic

	implicit none

	logical is_boundary_node
	integer k

	include 'param.h'

	is_boundary_node = inodv(k) .ne. 0 .and. inodv(k) .ne. -2

	end

c****************************************************************

	function is_open_boundary_node(k)

	use mod_geom_dynamic

	implicit none

	logical is_open_boundary_node
	integer k

	include 'param.h'

	is_open_boundary_node = inodv(k) .gt. 0

	end

c****************************************************************

	function is_dry_node(k)

	use mod_geom_dynamic

	implicit none

	logical is_dry_node
	integer k

	include 'param.h'

	is_dry_node = inodv(k) .eq. -2

	end

c****************************************************************
c****************************************************************
c****************************************************************

