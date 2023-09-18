!
! $Id: sublnka.f,v 1.7 2009-05-21 09:24:00 georg Exp $
!
! topological set up routines
!
! contents :
!
! subroutine setnod
!			 sets array inodv
!
! revision log :
!
! 01.08.2003    ggu     created from sublnk.f
! 13.08.2003    ggu     in update_geom do not call setweg and setnod
! 06.11.2008    ggu     better error handling
! 06.04.2009    ggu     nlidim -> nlkdim
! 02.12.2011    ggu     print_bound_nodes() for writing out boundary nodes
! 04.05.2015    ggu     remove equivalence - use winkv as local array
! 20.05.2015    ggu     new call to mklenkii to set up lenkiiv
!
!*****************************************************************
!-----------------------------------------------------------------
        module topological
!-----------------------------------------------------------------
        contains
!-----------------------------------------------------------------

	subroutine update_geom

! updates geometrical array (ieltv)

	use geom
	use geom_dynamic
	!use basin, only : nkn,nel,ngr,mbw
	use basin
	use shympi
        use mpi_io_admin
        use links

        implicit none

! local
        integer ie,ii
	integer iaux(3,nel)
	integer iiaux(nel)
	integer i,ia,ic,nc

!-------------------------------------------------------------
! update ieltv
!-------------------------------------------------------------

        if(bmpi) call up_ielt_mpi(inodv)
        call update_ielt(nel,inodv,ieltv)
        call make_aux(auxv_iei,ieltv)

	return

	!call shympi_comment('exchanging ieltv')
	iiaux(:) = ieltv(1,:)
	call shympi_exchange_2d_elem(iiaux)
	iaux(1,:) = iiaux(:)
	iiaux(:) = ieltv(2,:)
	call shympi_exchange_2d_elem(iiaux)
	iaux(2,:) = iiaux(:)
	iiaux(:) = ieltv(3,:)
	call shympi_exchange_2d_elem(iiaux)
	iaux(3,:) = iiaux(:)
	!call shympi_barrier

        write(my_unit,*) 'printing ghost elems: ' // 'ieltv'
        write(my_unit,*) 'n_ghost_areas = ',n_ghost_areas,my_id

        do ia=1,n_ghost_areas
          ic = ghost_areas(1,ia)
          nc = ghost_areas(4,ia)
          write(my_unit,*) 'elems: ',ic,nc
          do i=1,nc
            ie = ghost_elems(i,ia)
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

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!*****************************************************************

        subroutine setnod

! sets (dynamic) array inodv
!
! inodv 
!	 0: internal node  
!	>0: open boundary node
!       -1: boundary node  
!	-2: out of system
!
! if open boundary node, inodv(k) is number of boundary (ggu 15.11.2001)

	use geom_dynamic
	use evgeom
	use basin
	use shympi
        use fem_util
        use bnd_admin

        implicit none

	double precision, parameter :: winmax = 359.99
        integer ie,ii,k,n
        integer ibc,ibtyp
	integer nbc
        double precision winkv(nkn)

! initialize array to hold angles

	winkv = 0.d0

! sum angles

        do ie=1,nel
          if(iwegv(ie).eq.0) then !element is in system
            do ii=1,3
              k=nen3v(ii,ie)
              winkv(k)=winkv(k)+ev(10+ii,ie)
            end do
          end if
	end do

        call shympi_exchange_and_sum_2D_nodes(winkv)
        !call shympi_comment('shympi_elem: exchange winkv')

! set up inodv

        do k=1,nkn
          if(winkv(k).gt.winmax) then     !internal node
            inodv(k)=0
          else if(winkv(k).eq.0.) then  !out of system
            inodv(k)=-2
          else                          !boundary node
            inodv(k)=-1
          end if
        end do

! now mark open boundary nodes

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
              if(k .gt. 0) then
                if(inodv(k).ne.-1) goto 99
              end if
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
	write(6,*) ibc,n,ii,k,ipext(k),inodv(k),my_id
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

!****************************************************************
!****************************************************************
!****************************************************************

	function is_internal_node(k)

	use geom_dynamic

	implicit none

	logical is_internal_node
	integer k

	include 'param.h'

	is_internal_node = inodv(k) .eq. 0

	end

!****************************************************************

	function is_boundary_node(k)

	use geom_dynamic

	implicit none

	logical is_boundary_node
	integer k

	include 'param.h'

	is_boundary_node = inodv(k) .ne. 0 .and. inodv(k) .ne. -2

	end

!****************************************************************

	function is_open_boundary_node(k)

	use geom_dynamic

	implicit none

	logical is_open_boundary_node
	integer k

	include 'param.h'

	is_open_boundary_node = inodv(k) .gt. 0

	end

!****************************************************************

	function is_dry_node(k)

	use geom_dynamic

	implicit none

	logical is_dry_node
	integer k

	include 'param.h'

	is_dry_node = inodv(k) .eq. -2

	end

!*****************************************************************

	subroutine set_geom

! sets up geometrical arrays

	use geom
	use basin
	use shympi
        use links

        implicit none

	logical bverbose
        integer i,n
        integer nli,nbn,nin,nis,ngrd,ngrd1,ngrd2
        integer nlkdi
	integer nkbnds,kbnds,nodes,ibc,k

	bverbose = .true.
	bverbose = .false.

!-------------------------------------------------------------
! get dimension for link array
!-------------------------------------------------------------

	nlkdi = 3*nel+2*nkn
	if( ngr .gt. maxlnk ) goto 98

!-------------------------------------------------------------
! make static arrays
!-------------------------------------------------------------

        call mklenk(nlkdi,nkn,nel,nen3v,ilinkv,lenkv)
        call mklenkii(nlkdi,nkn,nel,nen3v,ilinkv,lenkv,lenkiiv)
        call mklink(nkn,ilinkv,lenkv,linkv)

        call mkkant(nkn,ilinkv,lenkv,linkv,kantv)

!        call shympi_exchange_halo_3d_nodes(2,kantv)

        call mkielt(nkn,nel,ilinkv,lenkv,linkv,ieltv)

!-------------------------------------------------------------
! write some statistics
!-------------------------------------------------------------

        ngrd=ilinkv(nkn+1)
        n=0
        do i=1,ngrd
          if(lenkv(i).eq.0) n=n+1
        end do

	nbn = n
	nin = nkn - nbn
	nis = (nel-2*nkn+nbn+2)/2

	if( bverbose ) then
          write(6,*) 'nel                      : ',nel
          write(6,*) 'nkn                      : ',nkn
          write(6,*) 'nbn                      : ',nbn
          write(6,*) 'nin                      : ',nin
          write(6,*) 'nis                      : ',nis
	end if

	ngrd1 = 4*nbn+6*nin+6*(nis-1)
        ngrd2 = 3*nel+nbn
	nli = ngrd/2

	if( bverbose ) then
          write(6,*) 'dimension                : ',nlkdi
          write(6,*) 'grades                   : ',ngrd
          write(6,*) 'formula (nis)            : ',ngrd1
          write(6,*) 'formula (nel)            : ',ngrd2
          write(6,*) 'links                    : ',nli
	end if

	if( ngrd .ne. ngrd1 ) goto 99
	if( ngrd .ne. ngrd2 ) goto 99

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	return
   98	continue
	write(6,*) 'ngr,maxlnk: ',ngr,maxlnk
	write(6,*) 'Please adjust maxlnk in links.h'
	stop 'error stop set_geom: maxlnk'
   99	continue
	stop 'error stop set_geom: ngrade'
	end

!*****************************************************************

	subroutine check_geom

! checks geometrical arrays

	use geom
	use basin
	use shympi
        use links

        implicit none

!-------------------------------------------------------------
! check static arrays
!-------------------------------------------------------------

        call checklenk(nkn,ilinkv,lenkv,nen3v)
        call checklink(nkn,ilinkv,linkv)

	if(.not.shympi_partition_on_elements()) call checkkant(nkn,kantv)
        call checkielt(nel,ieltv)

	call check_subs

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!*****************************************************************

        subroutine check_subs

! checks various subroutines

	use geom
	use basin, only : nkn,nel,ngr,mbw
        use lnku

        implicit none

	integer k,n,ne,ipf,ipl,ibase
	integer nl,i,ip,ie1,ie2,nn,k1,k2
	integer nodes(maxlnk)
	integer elems(maxlnk)

!-------------------------------------------------------------
! check element index
!-------------------------------------------------------------

	do k=1,nkn
	  call get_elem_linkp(k,ipf,ipl)
	  call get_elem_links(k,n,ibase)
	  call get_elems_around(k,maxlnk,ne,elems)
	  nl = ipl-ipf+1
	  if( n .ne. ne .or. nl .ne. ne ) goto 99
	  if( ipf .ne. ibase + 1 ) goto 96
	  if( ipl .ne. ibase + n ) goto 96
	  do i=1,n
	    ip = ipf + i - 1
	    if( ibase + i .ne. ip ) goto 98
	    ie1 = lenkv(ibase+i)
	    ie2 = elems(i)
	    if( ie1 .ne. ie2 ) goto 97
	  end do
	end do

!-------------------------------------------------------------
! check node index
!-------------------------------------------------------------

	do k=1,nkn
	  call get_node_linkp(k,ipf,ipl)
	  call get_node_links(k,n,ibase)
	  call get_nodes_around(k,maxlnk,nn,nodes)
	  nl = ipl-ipf+1
	  if( n .ne. nn .or. nl .ne. nn ) goto 89
	  if( ipf .ne. ibase + 1 ) goto 86
	  if( ipl .ne. ibase + n ) goto 86
	  do i=1,n
	    ip = ipf + i - 1
	    if( ibase + i .ne. ip ) goto 88
	    k1 = linkv(ibase+i)
	    k2 = nodes(i)
	    if( k1 .ne. k2 ) goto 87
	  end do
	end do

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	return
   86	continue
	write(6,*) k
	write(6,*) 'ipf,ipl,ibase,n: ',ipf,ipl,ibase,n
	stop 'error stop check_subs: extreme base pointers'
   87	continue
	write(6,*) k
	write(6,*) 'k1,k2: ',k1,k2
	stop 'error stop check_subs: node numbers'
   88	continue
	write(6,*) k
	write(6,*) 'i: ',i,ip,ibase + i
	stop 'error stop check_subs: base pointer'
   89	continue
	write(6,*) k
	write(6,*) 'n,nn,nl: ',n,nn,nl
	stop 'error stop check_subs: total number of nodes'
   96	continue
	write(6,*) k
	write(6,*) 'ipf,ipl,ibase,n: ',ipf,ipl,ibase,n
	stop 'error stop check_subs: extreme base pointers'
   97	continue
	write(6,*) k
	write(6,*) 'ie1,ie2: ',ie1,ie2
	stop 'error stop check_subs: element numbers'
   98	continue
	write(6,*) k
	write(6,*) 'i: ',i,ip,ibase + i
	stop 'error stop check_subs: base pointer'
   99	continue
	write(6,*) k
	write(6,*) 'n,ne,nl: ',n,ne,nl
	stop 'error stop check_subs: total number of elements'
	end

!****************************************************************
!****************************************************************
!****************************************************************

        subroutine print_bound_nodes(nkn,aux)

	use geom
        use fem_util

        implicit none

        integer nkn
        double precision aux(nkn)

        integer ib,k,kn,kstart

        ib = 0
        do k=1,nkn
          aux(k) = 0.
        end do

        open(1,file='bnd_nodes.dat',form='formatted',status='unknown')
        do k=1,nkn
          if( aux(k) .eq. 0. ) then
            if( kantv(1,k) .gt. 0 ) then
              kstart = k
              kn = kantv(1,kstart)
              ib = ib + 1
              write(1,*) 'bound ',ib
              write(1,*) ipext(kstart)
              do while( kn .ne. kstart )
                aux(kn) = 1.
                write(1,*) ipext(kn)
                kn = kantv(1,kn)
              end do
            end if
          end if
          aux(k) = 1.
        end do

        close(1)

        end 

!****************************************************************
!****************************************************************
!****************************************************************

!-----------------------------------------------------------------
        end module topological
!-----------------------------------------------------------------
