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

	subroutine set_geom

c sets up geometrical arrays

	use mod_geom
	use basin

        implicit none

c common
	include 'param.h'

c local
	logical bverbose
        integer i,n
        integer nli,nbn,nin,nis,ngrd,ngrd1,ngrd2
        integer nlkdi

	bverbose = .false.

c-------------------------------------------------------------
c get dimension for link array
c-------------------------------------------------------------

	nlkdi = 3*nel+2*nkn
	if( ngr .gt. maxlnk ) goto 98

c-------------------------------------------------------------
c make static arrays
c-------------------------------------------------------------

        call mklenk(nlkdi,nkn,nel,nen3v,ilinkv,lenkv)
        call mklenkii(nlkdi,nkn,nel,nen3v,ilinkv,lenkv,lenkiiv)
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

	subroutine check_geom

c checks geometrical arrays

	use mod_geom
	use basin

        implicit none

c common
	include 'param.h'


c-------------------------------------------------------------
c check static arrays
c-------------------------------------------------------------

        call checklenk(nkn,ilinkv,lenkv,nen3v)
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

	use mod_geom
	use basin, only : nkn,nel,ngr,mbw

        implicit none

	include 'param.h'


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

c****************************************************************
c****************************************************************
c****************************************************************

        subroutine print_bound_nodes(nkn,aux)

	use mod_geom

        implicit none

        integer nkn
        real aux(nkn)

	include 'param.h'

        integer ib,k,kn,kstart

        integer ipext

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

c****************************************************************

