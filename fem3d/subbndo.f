
!--------------------------------------------------------------------------
!
!    Copyright (C) 2001,2003-2004,2007-2012,2014-2015  Georg Umgiesser
!    Copyright (C) 2017-2019  Georg Umgiesser
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

c routines for open boundary conditions
c
c contents :
c
c subroutine bndo_init
c       sets up bndo data structure
c subroutine bndinsert(ib,area,kn)
c       inserts node kn with weight area into list (internal routine)
c
c subroutine bndo_info(iunit)
c       writes info on open boundary nodes to terminal
c
c function is_zeta_bound(k)
c	checks if node k is a zeta boundary
c
c subroutine bndo_setbc(what,nlvddi,cv,rbc,uprv,vprv)
c	sets open boundary condition for level boundaries
c subroutine bndo_impbc(what,nlvddi,cv,rbc)
c       imposes boundary conditions on open boundary
c subroutine bndo_adjbc(nlvddi,cv,uprv,vprv)
c       adjusts boundary conditions on open boundary (values on bnd already set)
c
c subroutine bndo_radiat(it,rzv)
c	imposes radiation condition for levels
c
c notes :
c
c	integer kbcdim
c	integer kopdim
c	parameter ( kbcdim = 100 )	!total number of open boundary nodes
c	parameter ( kopdim = 10 )	!maximum number of nodes close to OB
c
c       integer nbndo                   !total number of OB nodes
c
c	real xynorm(2,kbcdim)		!normal direction for OB node
c
c	integer iopbnd(nkn)		!if >0 pointer into array irv
c					!if <0 internal boundary (= -ibc)
c
c	integer ibcnod(kbcdim)		!number of boundary
c	integer kbcnod(kbcdim)		!number of boundary node
c	integer itynod(kbcdim)          !type of boundary
c
c	integer nopnod(kbcdim)		!number of internal nodes close to OB
c	integer nopnodes(kopdim,kbcdim)	!nodes close to OB
c
c	real wopnodes(kopdim,kbcdim)	!weights of nodes close to OB
c
c	iopbnd(k) = 0            no open BC
c	iopbnd(k) > 0            external open BC (ibtyp=1,2)
c	iopbnd(k) < 0            internal open BC (ibtyp=3)
c
c all dimensions are dynamical and not static
c they are setup in mod_bndo_init()
c
c the initialization routines should be called only after the
c ... arrays kantv and ieltv have been setup
c
c--------------------------------------------------------------------------
c
c revision log :
c
c 15.01.2001	ggu	written from scratch
c 03.12.2001	ggu	LEVMX - look out for missing level of near node
c 05.12.2001	ggu	NTBC - BUG -> has not been set before
c 27.03.2003	ggu	in bndo_adjbc use ambient value (bamb)
c 13.03.2004	ggu	in bndo_adjbc only for level BC (LEVELBC)
c 05.10.2004	ggu	new routine bndo_radiat, ibtyp=31
c 31.05.2007	ggu	reset BC for flux to old type (DEBHELP)
c 23.08.2007	ggu	use iopbnd as indicator for ext/int boundary nodes
c 08.04.2008	ggu	file cleaned (new bndo_setbc, bndo_impbc)
c 17.04.2008	ggu	calls to infobnd deleted (subst by get_bnd_ipar)
c 03.09.2008	ggu	new routine bndo_info_file()
c 06.11.2008	ggu	better error handling
c 12.11.2009	ggu	new array itynod and is_zeta_bound()
c 23.03.2010	ggu	changed v6.1.1
c 17.02.2011	ggu	changed VERS_6_1_18
c 25.03.2011	ggu	bug fix in bndo_impbc() -> ibcold not initialized
c 14.04.2011	ggu	changed VERS_6_1_22
c 30.03.2012	ggu	changed VERS_6_1_51
c 18.06.2014	ggu	changed VERS_6_1_77
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.06.2015	ggu	changed VERS_7_1_12
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 18.12.2015	ggu	changed VERS_7_3_17
c 05.12.2017	ggu	changed VERS_7_5_39
c 03.04.2018	ggu	changed VERS_7_5_43
c 19.04.2018	ggu	changed VERS_7_5_45
c 16.02.2019	ggu	changed VERS_7_5_60
c 26.04.2022	ggu	implementing OB in more than one domain
c 03.05.2022	ggu	lots of debug code integrated
c 22.03.2023	ggu	relax error condition on one node open boundary
c 09.05.2023    lrp     introduce top layer index variable
c
c***********************************************************************

	subroutine bndo_init

c sets up bndo data structure

	use mod_bndo
	use mod_bound_geom
	use mod_geom
	use basin
	use shympi

	implicit none

	logical bexternal
	logical berror,bdebug
	integer k,nodes,itype
	integer i,ibc
	integer inext,ilast,knext,klast
	integer ie,n,ngood,ie_mpi,id
	integer ii,iii,ib,in,kn,nb,j
	integer nbc
	integer iunit,kint,kext
	real area
	real dx,dy

	integer nkbnds,itybnd,kbnds,ipext,ipint,nbnds
	real areaele

c----------------------------------------------------------
c set up array iopbnd
c----------------------------------------------------------

	iopbnd(:) = 0

	nbndo = 0
        ndebug = 0              !unit number for debug (in common block)

	nbc = nbnds()

	do ibc = 1,nbc
	  nodes = nkbnds(ibc)
	  itype = itybnd(ibc)

	  ngood = 0
	  bexternal = ( itype .ge. 1 .and. itype .le. 2 
     +      		    .or. itype .ge. 31 .and. itype .le. 39 )

	  do i=1,nodes
	    k = kbnds(ibc,i)
	    if( k <= 0 ) cycle
	    !if( bexternal ) then
	    if( .true. ) then
	      ngood = ngood + 1
	      nbndo = nbndo + 1
	      if( nbndo .gt. kbcdim ) goto 99
	      iopbnd(k) = nbndo
	      ibcnod(nbndo) = ibc
	      kbcnod(nbndo) = k
	      itynod(nbndo) = itype
	      bexnod(nbndo) = bexternal		!flag this node as external
	    else
	      iopbnd(k) = -ibc
	    end if
	  end do

	  if( ngood == 0 ) then
	    !boundary not in domain
	  else if( ngood == nodes ) then
	    !boundary fully in domain
	  else if( itype == 2 ) then	!boundary only partially in domain
	    write(6,*) 'ngood,nodes: ',ngood,nodes
	    write(6,*) 'boundary is only partially in domain'
	    write(6,*) 'cannot handle flux OB in different domains yet'
	    stop 'error stop bndo_init: boundary between domains'
	  else			!zeta boundary - should be able to handle
	    !
	  end if
	end do

c----------------------------------------------------------
c set up normal direction
c----------------------------------------------------------

	do i=1,nbndo
	  k = kbcnod(i)
	  ibc = ibcnod(i)

	  itype = itybnd(ibc)
	  bexternal = bexnod(i)

	  if( .not. bexternal ) cycle

	  knext = kantv(1,k)
	  klast = kantv(2,k)

	  inext = 0
	  ilast = 0
	  if( knext > 0 ) inext = iopbnd(knext)
	  if( klast > 0 ) ilast = iopbnd(klast)

c	  -------------------------------
c	  internal consistency check
c	  -------------------------------

	  if( iopbnd(k) .ne. i ) then
	    stop 'error stop bndo_init: internal error (0)'
	  end if
	  if( inext .gt. 0 ) then
	   if( kbcnod(inext) .ne. knext ) then
	    stop 'error stop bndo_init: internal error (1)'
	   end if
	  end if
	  if( ilast .gt. 0 ) then
	   if( kbcnod(ilast) .ne. klast ) then
	    stop 'error stop bndo_init: internal error (2)'
	   end if
	  end if

c	  -------------------------------
c	  adjacent boundary nodes must be of same boundary
c	  -------------------------------

	  if( inext .gt. 0 ) then
	    if( ibcnod(inext) .ne. ibc ) goto 98
	  end if
	  if( ilast .gt. 0 ) then
	    if( ibcnod(ilast) .ne. ibc ) goto 98
	  end if

c	  -------------------------------
c	  get normal direction (might be wrong in case of domain border)
c	  -------------------------------

	  if( inext .gt. 0 .and. ilast .gt. 0 ) then	!inner node in OB
	    dx = xgv(knext) - xgv(klast)
	    dy = ygv(knext) - ygv(klast)
	  else if( inext .gt. 0 ) then			!first node in OB
	    dx = xgv(knext) - xgv(k)
	    dy = ygv(knext) - ygv(k)
	  else if( ilast .gt. 0 ) then			!last node in OB
	    dx = xgv(k) - xgv(klast)
	    dy = ygv(k) - ygv(klast)
	  else
	    id = id_node(k)
	    if( id == my_id ) then
	      write(6,*) 'One node open boundary not permitted'
	      write(6,*) i,k,ipext(k)
	      write(6,*) 'node number (external): ',ipext(k)
	      write(6,*) 'domain: ',my_id
	      write(6,*) 'boundary: ',ibc
	      stop 'error stop bndo'
	    else
	      dx = 0.
	      dy = 0.
	    end if
	  end if

	  xynorm(1,i) = -dy			!x-component
	  xynorm(2,i) = dx			!y-component
	end do
  
c----------------------------------------------------------
c set up weights and node list
c----------------------------------------------------------

	do i=1,nbndo
	  nopnod(i) = 0
	end do

	do ie_mpi=1,nel

	  ie = ip_sort_elem(ie_mpi)
	  n = 3
	  area = areaele(ie)

	  do ii=1,n
	    k = nen3v(ii,ie)
	    ib = iopbnd(k)
	    bexternal = bexnod(ib)
	    !if( ib .gt. 0 ) then		!insert inner nodes
	    if( bexternal ) then		!insert inner nodes
	      do iii=1,n
		kn = nen3v(iii,ie)
	        in = iopbnd(kn)
		if( in .le. 0 ) then		!only inner nodes
		  call bndinsert(ib,area,kn)
		end if
	      end do
	    end if
	  end do

	end do

c----------------------------------------------------------
c scale weights to unit
c----------------------------------------------------------

	iunit = 730 + my_id
	kext = 6651
	kext = 0
	kint = ipint(kext)

	berror = .false.

	do i=1,nbndo

	  bexternal = bexnod(i)
	  if( .not. bexternal ) cycle

	  nb = nopnod(i)
	  if( nb .le. 0 ) then
	    k = kbcnod(i)
	    ibc = ibcnod(i)
	    !write(6,*) i,k,ipext(k)
	    write(6,*) '*** No inner nodes for boundary node'
	    write(6,*) '      boundary = ',ibc,'   node = ',ipext(k)
	    berror = .true.
	  end if

	  area = 0.
	  do j=1,nb
	    area = area + wopnodes(j,i)
	  end do
	  do j=1,nb
	    if( area .gt. 0. ) then
	      wopnodes(j,i) = wopnodes(j,i) / area
	    end if
	  end do
	  
	  bdebug = ( kint == kbcnod(i) )
	  if( bdebug ) then
	    iunit = 730 + my_id
	    write(iunit,*) '--------- bndo_init ------------'
	    write(iunit,*) i,kint,kext,nb,id_node(k)
	    write(iunit,*) (wopnodes(j,i),j=1,nb)
	    write(iunit,*) '--------- end bndo_init ------------'
	  end if
	end do

	if( berror ) stop 'error stop bndo'

	write(6,*) 'finished setting up bndo_init, nbndo = ',nbndo

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	return
   98	continue
	write(6,*) 'different boundary : ',ibc
	stop 'error stop bndo'
   99	continue
	write(6,*) 'dimension error kbcdim : ',kbcdim
	stop 'error stop bndo'
	end

c***********************************************************************

	subroutine bndinsert(ib,area,kn)

c inserts node kn with weight area into list (internal routine)

	use mod_bndo

	implicit none

	integer ib		!nodal index
	real area		!weight
	integer kn		!node number to insert

	integer nb,j,k

	nb = nopnod(ib)

	do j=1,nb
	  k = nopnodes(j,ib)
	  if( k .eq. kn ) then		!already there -> add weight
	    wopnodes(j,ib) = wopnodes(j,ib) + area
	    return			!done
	  end if
	end do

	nb = nb + 1
	if( nb .gt. kopdim ) then
	  write(6,*) 'Too much inner neighbors: ',nb
	  write(6,*) 'Please raise kopdim in subbndo.h'
	  stop 'error stop bndo'
	end if

	nopnodes(nb,ib) = kn
	wopnodes(nb,ib) = area
	nopnod(ib) = nb

	end

c***********************************************************************

	subroutine bndo_info_file(file)

c writes bndo info to file

	use shympi

	implicit none

	character*(*) file

	if( file .eq. ' ' ) return

        if(shympi_is_master()) then
	  open(1,file=file,status='unknown',form='formatted')
	  call bndo_info(1)
	  close(1)
	end if

	end

c***********************************************************************

	subroutine bndo_info(iunit)

c writes info on open boundary nodes to terminal

	use mod_bndo
	use mod_bound_geom
	use levels

	implicit none

        integer iunit

	integer i,k,ibc,nb,j,iu
	integer itybnd,ipext

        iu = iunit
        if( iu .le. 0 ) iu = 6

	write(iu,*) '--------------------------------'
	write(iu,*) 'Information on open boundary nodes'
	write(iu,*) 'Total number of open boundary nodes: ',nbndo

	do i=1,nbndo

	  k = kbcnod(i)
	  ibc = ibcnod(i)
	  nb = nopnod(i)

	  if( iopbnd(k) .ne. i ) then
	    stop 'error stop bndo_info: internal error (11)'
	  end if

	  write(iu,*) '-------------------------------- bndo_info'
	  write(iu,*) i,k,ipext(k),ibc,itybnd(ibc),ilhkv(k)
	  write(iu,*) nb
	  write(iu,*) (ipext(nopnodes(j,i)),j=1,nb)
	  write(iu,*) (wopnodes(j,i),j=1,nb)
	  write(iu,*) (ilhkv(nopnodes(j,i)),j=1,nb)
	  write(iu,*) (xynorm(j,i),j=1,2)
	end do

	write(iu,*) '--------------------------------'

	end

c***********************************************************************

	function is_zeta_bound(k)

c checks if node k is a zeta boundary

	use mod_bndo
	use mod_bound_geom

	implicit none

	logical is_zeta_bound
	integer k

	integer ip

	is_zeta_bound = .false.
	ip = iopbnd(k)

	if( ip .gt. 0 ) then
	  if( itynod(ip) .eq. 1 ) then
	    is_zeta_bound = .true.
	  end if
	end if
	
	end

c***********************************************************************
c***********************************************************************
c***********************************************************************

        subroutine bndo_setbc(what,nlvddi,cv,rbc,uprv,vprv)

c sets open boundary condition for level boundaries
c
c simply calls bndo_impbc() and bndo_adjbc()

	use basin
	use shympi

        implicit none

        character*(*) what      !conz/temp/salt or else
        integer nlvddi
        real cv(nlvddi,nkn)
        real rbc(nlvddi,nkn)	!boundary condition (3D)
	real uprv(nlvddi,nkn)
	real vprv(nlvddi,nkn)

	logical bdebug
	integer k,iunit,kext
	integer ipint

	iunit = 730 + my_id
	kext = 6651
	kext = 0
	k = ipint(kext)
	bdebug = ( k > 0 )

c----------------------------------------------------------
c simply imposes whatever is in rbc
c----------------------------------------------------------

	if( bdebug ) then
	write(iunit,*) '------------ bndo_setbc 1 -------'
	write(iunit,*) cv(:,k)
	write(iunit,*) rbc(:,k)
	end if
	
        call bndo_impbc(what,nlvddi,cv,rbc)

c----------------------------------------------------------
c adjusts for ambient value, no gradient or outgoing flow
c----------------------------------------------------------

	if( bdebug ) then
	write(iunit,*) '------------ bndo_setbc 2 -------'
	write(iunit,*) cv(:,k)
	write(iunit,*) rbc(:,k)
	end if
	
	call bndo_adjbc(what,nlvddi,cv,uprv,vprv)

	if( bdebug ) then
	write(iunit,*) '------------ bndo_setbc 3 -------'
	write(iunit,*) cv(:,k)
	write(iunit,*) rbc(:,k)
	end if
	
c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	end

c***********************************************************************

        subroutine bndo_impbc(what,nlvddi,cv,rbc)

c imposes boundary conditions on open boundary

	use mod_bndo
	use mod_bound_geom
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        character*(*) what      !conz/temp/salt or else
        integer nlvddi
        real cv(nlvddi,nkn)
        real rbc(nlvddi,nkn)	!boundary condition (3D)

        logical blevel,bfix
        logical bdebug
        integer i,j,k,l
        integer ibc,ibcold
        integer nb,lmax,lmin
        integer ibtyp
        real value,rb
	real, parameter :: flag = -999.

        integer ifemopa

        bdebug = .true.
        bdebug = .false.

	ibcold = 0

        do i=1,nbndo

          k = kbcnod(i)
          nb = nopnod(i)
          ibc = ibcnod(i)

	  if( ibc .ne. ibcold ) then
	    call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	    ibcold = ibc
	  end if

	  blevel = ibtyp .eq. 1
	  bfix = ibtyp .eq. 5
          lmax = ilhkv(k)
          lmin = jlhkv(k)

	!if( bfix ) then
	!write(6,*) 'bfix: ',ibc,ibtyp
	!end if

          if( iopbnd(k) .ne. i ) then
	    stop 'error stop bndo_impbc: internal error (11)'
	  end if

	  if( blevel .or. bfix ) then
	if( bfix ) then
	!write(6,*) 'fixing bound: ',k,ibtyp,ibc,rbc(1,k)
	end if
            do l=lmin,lmax
	      rb = rbc(l,k)
              if( rb .ne. flag ) cv(l,k) = rb
            end do
	  end if

        end do

        end

c***********************************************************************

	subroutine bndo_adjbc(what,nlvddi,cv,uprv,vprv)

c adjusts boundary conditions on open boundary (values on bnd already set)
c
c adjusts for ambient value, no gradient or outgoing flow
c
c this is only done one level boundaries ( ibtyp == 1 )

	use mod_bndo
	use mod_bound_geom
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

        character*(*) what	!conz/temp/salt or else
	integer nlvddi
	real cv(nlvddi,nkn)
	real uprv(nlvddi,nkn)
	real vprv(nlvddi,nkn)

	logical bgrad0
	logical blevel
	logical bdebug
	logical bdggu
	logical bout,bamb
	logical binside
	integer i,j,k,l
	integer ibtyp,igrad0
	integer ibc,ibcold
	integer nb,nlev,flev,ko
        integer ntbc,nlevko
	real dx,dy
	real scal,bc,weight,tweight
	real value
	character*20 aline

	integer kint,kext,iunit
	integer ipext,ipint

	integer ifemopa

	bdebug = .true.
	bdebug = .false.

	ibcold = 0

	bgrad0 = .false.
	blevel = .false.

	iunit = 730 + my_id
	kext = 6651
	kext = 0
	kint = ipint(kext)

	if( bdebug ) then
	  if( ndebug .eq. 0 ) then
	    ndebug = ifemopa('bndo_adjbc (91)','.bndo','form','unknown')
            call bndo_info(ndebug)
	  end if
	  call get_act_timeline(aline)
	  write(ndebug,*) 'bndo_adjbc ........... ',what,aline
	end if

	do i=1,nbndo

	  k = kbcnod(i)
	  nb = nopnod(i)
	  ibc = ibcnod(i)
	  bdebug = ( k == kint )

	  if( iopbnd(k) .ne. i ) then
	    stop 'error stop bndo_adjbc: internal error (11)'
	  end if

	  if( ibc .ne. ibcold ) then
	    call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	    call get_bnd_ipar(ibc,'igrad0',igrad0)
	    bgrad0 = igrad0 .gt. 0
	    blevel = ibtyp .eq. 1
	    ibcold = ibc
	  end if

          if( blevel ) then       !LEVELBC !DEBHELP

	  dx = xynorm(1,i)
	  dy = xynorm(2,i)
	  nlev = ilhkv(k)
	  flev = jlhkv(k)

	  do l=flev,nlev
	    scal = dx * uprv(l,k) + dy * vprv(l,k)
	    bout = scal .le. 0.				!outgoing flow
	    bamb = cv(l,k) .le. -990.			!make ambient value
	    binside = bgrad0 .or. bout .or. bamb
	    if( binside ) then				!take from inside
	      bc = 0.
              ntbc = 0
              tweight = 0.
	      do j=1,nb
		ko = nopnodes(j,i)
                nlevko = ilhkv(ko)
                if( l .le. nlevko ) then        !LEVMX  !only if level exists
                  ntbc = ntbc + 1
		  weight = wopnodes(j,i)
                  tweight = tweight + weight
		  bc = bc + weight * cv(l,ko)
                end if
	      end do
	    else				!impose boundary value
              ntbc = nb                         !NTBC - BUG -> has not been set
	      tweight = 1.			!prob useless
	      bc = cv(l,k)
	    end if

            if( ntbc .gt. 0 ) then              !LEVMX  !at least 1 node found
	      value = bc / tweight
	      if( cv(l,k) .le. -5555. ) then	!differential value
		value = value - 10000. - cv(l,k)
	      end if
	      cv(l,k) = value
            else if( l .gt. 1 ) then            !take from above
              cv(l,k) = cv(l-1,k)
            else
	      write(6,*) i,k,ipext(k),nb,ntbc,ibc,binside
	      write(6,*) nlev,ntbc
              stop 'error stop bndo_adjbc: internal error (5)'
            end if

	  end do

          end if

	end do

	end

c***********************************************************************

	subroutine bndo_radiat(it,rzv)

c imposes radiation condition for levels

	use mod_bndo
	use mod_bound_geom
	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer it
	real rzv(nkn)

	logical bdebug
	integer i,j,k,l
	integer ibc,ibcold
	integer nb,nlev,ko
        integer ntbc,nlevko
	integer ibtyp
	real bc,weight,tweight

	integer ipext

	integer ifemopa

	bdebug = .true.
	bdebug = .false.
	ibcold = 0

	if( bdebug ) then
	  if( ndebug .eq. 0 ) then
	    ndebug = ifemopa('bndo_adjbc (91)','.bndo','form','unknown')
            call bndo_info(ndebug)
	  end if
	  write(ndebug,*) 'bndo_adjbc ........... ','radiat',it
	end if

	do i=1,nbndo

	  k = kbcnod(i)
	  nb = nopnod(i)
	  ibc = ibcnod(i)

	  if( ibc .ne. ibcold ) then
	    call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	  end if
          !write(78,*) i,k,nb,ibc,ibtyp

	  if( iopbnd(k) .ne. i ) then
	    stop 'error stop bndo_radiat: internal error (11)'
	  end if

          if( ibtyp .eq. 31 ) then       !radiation condition only

	    bc = 0.
            ntbc = 0
            tweight = 0.
	    do j=1,nb
		ko = nopnodes(j,i)
                ntbc = ntbc + 1
		weight = wopnodes(j,i)
                tweight = tweight + weight
		bc = bc + weight * znv(ko)
	    end do

            if( ntbc .gt. 0 ) then         !LEVMX  !at least 1 node found
	      rzv(k) = bc / tweight
            else
              stop 'error stop bndo_radiat: internal error (5)'
            end if

          end if

	end do

	end

c***********************************************************************

