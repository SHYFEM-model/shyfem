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
	use basin

        implicit none

c local
        integer n

c-------------------------------------------------------------
c update ieltv
c-------------------------------------------------------------

        call update_ielt(nel,inodv,ieltv)

	!call exchange_ieltv

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c*****************************************************************

	subroutine exchange_ieltv

c exchanges ieltv structure

	use mod_geom
	use mod_geom_dynamic
	use basin
	use shympi

	implicit none

        integer ie,ii
	integer iaux(3,nel)
	integer iiaux(nel)
	integer i,ia,ic,nc

c-------------------------------------------------------------
c start exchanging
c-------------------------------------------------------------

	call shympi_comment('exchanging ieltv')
	iiaux(:) = ieltv(1,:)
	call shympi_exchange_2d_elem(iiaux)
	iaux(1,:) = iiaux(:)
	iiaux(:) = ieltv(2,:)
	call shympi_exchange_2d_elem(iiaux)
	iaux(2,:) = iiaux(:)
	iiaux(:) = ieltv(3,:)
	call shympi_exchange_2d_elem(iiaux)
	iaux(3,:) = iiaux(:)
	call shympi_barrier

c-------------------------------------------------------------
c print info
c-------------------------------------------------------------

        write(my_unit,*) 'printing ghost elems: ' // 'ieltv'
        write(my_unit,*) 'n_ghost_areas = ',n_ghost_areas,my_id

	do ia=1,n_ghost_areas
          ic = ghost_areas(1,ia)
          nc = ghost_areas(4,ia)
          write(my_unit,*) 'elems outer: ',ic,nc
          do i=1,nc
            ie = ghost_elems_out(i,ia)
            write(my_unit,*) ie,ipev(ie),(ieltv(ii,ie),ii=1,3)
          end do
          nc = ghost_areas(5,ia)
          write(my_unit,*) 'elems inner: ',ic,nc
          do i=1,nc
            ie = ghost_elems_in(i,ia)
            write(my_unit,*) ie,ipev(ie),(ieltv(ii,ie),ii=1,3)
          end do
        end do

	do ie=1,nel
	  do ii=1,3
	    if( ieltv(ii,ie) /= iaux(ii,ie) ) then
	      write(my_unit,*) 'ieltv: ',ie,ieltv(ii,ie),iaux(ii,ie)
	    end if
	  end do
	end do

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c****************************************************************

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
	use shympi

        implicit none

	double precision, parameter :: winmax = 359.99
        integer ie,ii,k,n
        integer ibc,ibtyp
	integer nbc,ndry
        double precision winkv(nkn)

        integer ipext
	integer itybnd,nkbnds,kbnds,nbnds

c initialize array to hold angles

	ndry = 0
	winkv = 0.

c sum angles

        do ie=1,nel
          if(iwegv(ie).eq.0) then !element is in system
            do ii=1,3
              k=nen3v(ii,ie)
              winkv(k)=winkv(k)+ev(10+ii,ie)
            end do
	  else
	    ndry = ndry + 1
          end if
	end do

	!if( ndry > 0 ) then
	!  write(6,*) 'dry elements: ',ndry,' / ',nel
	!end if

        !call shympi_comment('shympi_elem: exchange winkv')
        call shympi_exchange_and_sum_2d_nodes(winkv)

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
	      if( k <= 0 ) cycle
              if(inodv(k).eq.-1) inodv(k)=ibc
            end do
          else if(ibtyp.gt.0) then
            do ii=1,n
              k=kbnds(ibc,ii)
	      if( k <= 0 ) cycle
              if(inodv(k).ne.-1) goto 99
              inodv(k)=ibc
            end do
          end if
        end do

	!call shympi_comment('exchanging inodv')
	call shympi_exchange_2d_node(inodv)
	!call shympi_barrier

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
        stop 'error stop setnod: open boundary node'
        end

c****************************************************************
c****************************************************************
c****************************************************************

	function is_internal_node(k)

	use mod_geom_dynamic

	implicit none

	logical is_internal_node
	integer k

	is_internal_node = inodv(k) .eq. 0

	end

c****************************************************************

	function is_boundary_node(k)

	use mod_geom_dynamic

	implicit none

	logical is_boundary_node
	integer k

	is_boundary_node = inodv(k) .ne. 0 .and. inodv(k) .ne. -2

	end

c****************************************************************

	function is_open_boundary_node(k)

	use mod_geom_dynamic

	implicit none

	logical is_open_boundary_node
	integer k

	is_open_boundary_node = inodv(k) .gt. 0

	end

c****************************************************************

	function is_dry_node(k)

	use mod_geom_dynamic

	implicit none

	logical is_dry_node
	integer k

	is_dry_node = inodv(k) .eq. -2

	end

c****************************************************************
c****************************************************************
c****************************************************************

