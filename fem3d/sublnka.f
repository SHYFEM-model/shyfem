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
c
c*****************************************************************

	subroutine set_geom

c sets up geometrical arrays

        implicit none

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v
	include 'links.h'
        integer kantv(2,1)
        common /kantv/kantv
        integer ieltv(3,1)
        common /ieltv/ieltv
c local
        integer i,n
        integer nli,nbn,nin,nis,ngrd,ngrd1,ngrd2
        integer nlkdim,ngrdim

c-------------------------------------------------------------
c get dimension for link array
c-------------------------------------------------------------

	call getdim('nlkdim',nlkdim)
	write(6,*) 'set_geom: nlkdim = ',nlkdim
	if( ngr .gt. maxlnk ) goto 98

c-------------------------------------------------------------
c make static arrays
c-------------------------------------------------------------

        call mklenk(nlkdim,nkn,nel,nen3v,ilinkv,lenkv)
        call mklink(nkn,ilinkv,lenkv,linkv)

        call mkkant(nkn,ilinkv,lenkv,linkv,kantv)
        call mkielt(nkn,nel,ilinkv,lenkv,linkv,ieltv)

c-------------------------------------------------------------
c write some statistics
c-------------------------------------------------------------

        ngrd=ilinkv(nkn+1)
        n=0
        do i=1,ngrd
          if(lenkv(i).eq.0) n=n+1
        end do

	nbn = n
	nin = nkn - nbn
	nis = (nel-2*nkn+nbn+2)/2

        write(6,*) 'nel                      : ',nel
        write(6,*) 'nkn                      : ',nkn
        write(6,*) 'nbn                      : ',nbn
        write(6,*) 'nin                      : ',nin
        write(6,*) 'nis                      : ',nis

	ngrd1 = 4*nbn+6*nin+6*(nis-1)
        ngrd2 = 3*nel+nbn
	nli = ngrd/2

        write(6,*) 'dimension                : ',nlkdim
        write(6,*) 'grades                   : ',ngrd
        write(6,*) 'formula (nis)            : ',ngrd1
        write(6,*) 'formula (nel)            : ',ngrd2
        write(6,*) 'links                    : ',nli

	if( ngrd .ne. ngrd1 ) goto 99
	if( ngrd .ne. ngrd2 ) goto 99

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	return
   98	continue
	write(6,*) 'ngr,maxlnk: ',ngr,maxlnk
	write(6,*) 'Please adjust maxlnk in links.h'
	stop 'error stop set_geom: maxlnk'
   99	continue
	stop 'error stop set_geom: ngrade'
	end

c*****************************************************************

	subroutine update_geom

c updates geometrical array (ieltv)

        implicit none

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer inodv(1)
        common /inodv/inodv
        integer ieltv(3,1)
        common /ieltv/ieltv
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

	subroutine check_geom

c checks geometrical arrays

        implicit none

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	include 'links.h'
        integer kantv(2,1)
        common /kantv/kantv
        integer ieltv(3,1)
        common /ieltv/ieltv

c-------------------------------------------------------------
c check static arrays
c-------------------------------------------------------------

        call checklenk(nkn,ilinkv,lenkv)
        call checklink(nkn,ilinkv,linkv)

        call checkkant(nkn,kantv)
        call checkielt(nel,ieltv)

	call check_subs

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c*****************************************************************

        subroutine check_subs

c checks various subroutines

        implicit none

        include 'links.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer k,n,ne,ipf,ipl,ibase
	integer nl,i,ip,ie1,ie2,nn,k1,k2

c-------------------------------------------------------------
c check element index
c-------------------------------------------------------------

	do k=1,nkn
	  call get_elem_linkp(k,ipf,ipl)
	  call get_elem_links(k,n,ibase)
	  call set_elem_links(k,ne)
	  nl = ipl-ipf+1
	  if( n .ne. ne .or. nl .ne. ne ) goto 99
	  if( ipf .ne. ibase + 1 ) goto 96
	  if( ipl .ne. ibase + n ) goto 96
	  do i=1,n
	    ip = ipf + i - 1
	    if( ibase + i .ne. ip ) goto 98
	    ie1 = lenkv(ibase+i)
	    ie2 = lnk_elems(i)
	    if( ie1 .ne. ie2 ) goto 97
	  end do
	end do

c-------------------------------------------------------------
c check node index
c-------------------------------------------------------------

	do k=1,nkn
	  call get_node_linkp(k,ipf,ipl)
	  call get_node_links(k,n,ibase)
	  call set_node_links(k,nn)
	  nl = ipl-ipf+1
	  if( n .ne. nn .or. nl .ne. nn ) goto 89
	  if( ipf .ne. ibase + 1 ) goto 86
	  if( ipl .ne. ibase + n ) goto 86
	  do i=1,n
	    ip = ipf + i - 1
	    if( ibase + i .ne. ip ) goto 88
	    k1 = linkv(ibase+i)
	    k2 = lnk_nodes(i)
	    if( k1 .ne. k2 ) goto 87
	  end do
	end do

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

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

c*****************************************************************
c
        subroutine setnod
c
c sets (dynamic) array inodv
c
c inodv 
c	 0: internal node  
c	>0: open boundary node
c       -1: boundary node  
c	-2: out of system
c
c if open boundary node, inodv(k) is number of boundary (ggu 15.11.2001)
c
        implicit none
c
c parameter
	real winmax
	parameter(winmax=359.8)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer iwegv(1),nen3v(3,1),inodv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /iwegv/iwegv, /nen3v/nen3v, /inodv/inodv
	include 'ev.h'
c local
        integer ie,ii,k,n
        integer ibc,ibtyp
c functions
        integer ipext
	integer itybnd,nkbnds,kbnds
c equivalence
        real winkv(1)
        equivalence(inodv(1),winkv(1))
c
c initialize array to hold angles
c
        do k=1,nkn
          winkv(k)=0.
        end do
c
c sum angles
c
        do ie=1,nel
          if(iwegv(ie).eq.0) then !element is in system
            do ii=1,3
              k=nen3v(ii,ie)
              winkv(k)=winkv(k)+ev(10+ii,ie)
            end do
          end if
	end do
c
c set up inodv
c
        do k=1,nkn
          if(winkv(k).gt.winmax) then     !internal node
            inodv(k)=0
          else if(winkv(k).eq.0.) then  !out of system
            inodv(k)=-2
          else                          !boundary node
            inodv(k)=-1
          end if
        end do
c
c now mark open boundary nodes
c
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
c
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

	implicit none

	logical is_internal_node
	integer k

        integer inodv(1)
        common /inodv/inodv

	is_internal_node = inodv(k) .eq. 0

	end

c****************************************************************

	function is_boundary_node(k)

	implicit none

	logical is_boundary_node
	integer k

        integer inodv(1)
        common /inodv/inodv

	is_boundary_node = inodv(k) .ne. 0 .and. inodv(k) .ne. -2

	end

c****************************************************************

	function is_open_boundary_node(k)

	implicit none

	logical is_open_boundary_node
	integer k

        integer inodv(1)
        common /inodv/inodv

	is_open_boundary_node = inodv(k) .gt. 0

	end

c****************************************************************

	function is_dry_node(k)

	implicit none

	logical is_dry_node
	integer k

        integer inodv(1)
        common /inodv/inodv

	is_dry_node = inodv(k) .eq. -2

	end

c****************************************************************
c****************************************************************
c****************************************************************

