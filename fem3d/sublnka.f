
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2008-2015,2019  Georg Umgiesser
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
c subroutine setnod
c			 sets array inodv
c
c revision log :
c
c 01.08.2003	ggu	created from sublnk.f
c 13.08.2003	ggu	in update_geom do not call setweg and setnod
c 06.11.2008	ggu	better error handling
c 06.04.2009	ggu	nlidim -> nlkdim
c 23.03.2010	ggu	changed v6.1.1
c 02.12.2011	ggu	print_bound_nodes() for writing out boundary nodes
c 09.12.2011	ggu	changed VERS_6_1_38
c 30.03.2012	ggu	changed VERS_6_1_51
c 12.09.2013	ggu	changed VERS_6_1_67
c 18.06.2014	ggu	changed VERS_6_1_77
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 04.05.2015	ggu	remove equivalence - use winkv as local array
c 20.05.2015	ggu	new call to mklenkii to set up lenkiiv
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 24.07.2015	ggu	changed VERS_7_1_82
c 10.10.2015	ggu	changed VERS_7_3_2
c 16.12.2015	ggu	changed VERS_7_3_16
c 16.02.2019	ggu	changed VERS_7_5_60
c 20.05.2020	ggu	new way to compute link structure (still experimental)
c 28.05.2020	ggu	some more changes in constructing link structure
c
c*****************************************************************

	subroutine set_geom

c sets up geometrical arrays

	use mod_geom
	use basin

        implicit none

	logical bverbose
        integer i,n
        integer nli,nbn,nin,nis,ngrd,ngrd1,ngrd2
        integer kerr

	bverbose = .true.
	bverbose = .false.

c-------------------------------------------------------------
c check maxlnk
c-------------------------------------------------------------

	if( ngr .gt. maxlnk ) goto 98

c-------------------------------------------------------------
c make static arrays
c-------------------------------------------------------------

	call make_links_old(nkn,nel,nen3v)
	call make_links(nkn,nel,nen3v,kerr)

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

	integer k,n,ne,ipf,ipl,ibase
	integer nl,i,ip,ie1,ie2,nn,k1,k2
	integer nodes(maxlnk)
	integer elems(maxlnk)

c-------------------------------------------------------------
c check element index
c-------------------------------------------------------------

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

c-------------------------------------------------------------
c check node index
c-------------------------------------------------------------

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

c****************************************************************
c****************************************************************
c****************************************************************

        subroutine print_bound_nodes(nkn,aux)

	use mod_geom

        implicit none

        integer nkn
        real aux(nkn)

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

