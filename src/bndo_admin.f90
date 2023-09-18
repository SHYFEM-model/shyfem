!
! $Id: subbndo.f,v 1.18 2010-03-11 15:36:39 georg Exp $
!
! routines for open boundary conditions
!
! contents :
!
! subroutine bndo_init
!       sets up bndo data structure
! subroutine bndinsert(ib,area,kn)
!       inserts node kn with weight area into list (internal routine)
!
! subroutine bndo_info(iunit)
!       writes info on open boundary nodes to terminal
!
! function is_zeta_bound(k)
!	checks if node k is a zeta boundary
!
! subroutine bndo_setbc(it,what,nlvddi,cv,rbc,uprv,vprv)
!	sets open boundary condition for level boundaries
! subroutine bndo_impbc(it,what,nlvddi,cv,rbc)
!       imposes boundary conditions on open boundary
! subroutine bndo_adjbc(it,nlvddi,cv,uprv,vprv)
!       adjusts boundary conditions on open boundary (values on bnd already set)
!
! subroutine bndo_radiat(it,rzv)
!	imposes radiation condition for levels
!
! notes :
!
!	integer kbcdim
!	integer kopdim
!	parameter ( kbcdim = 100 )	!total number of open boundary nodes
!	parameter ( kopdim = 10 )	!maximum number of nodes close to OB
!
!       integer nbndo                   !total number of OB nodes
!
!	double precision xynorm(2,kbcdim)		!normal direction for OB node
!
!	integer iopbnd(nknddi)		!if >0 pointer into array irv
!					!if <0 internal boundary (= -ibc)
!
!	integer ibcnod(kbcdim)		!number of boundary
!	integer kbcnod(kbcdim)		!number of boundary node
!	integer itynod(kbcdim)          !type of boundary
!
!	integer nopnod(kbcdim)		!number of internal nodes close to OB
!	integer nopnodes(kopdim,kbcdim)	!nodes close to OB
!
!	double precision wopnodes(kopdim,kbcdim)	!weights of nodes close to OB
!
!	iopbnd(k) = 0            no open BC
!	iopbnd(k) > 0            external open BC (ibtyp=1,2)
!	iopbnd(k) < 0            internal open BC (ibtyp=3)
!
! the initialization routines should be called only after the
! ... arrays kantv and ieltv have been setup
!
! revision log :
!
! 15.01.2001    ggu     written from scratch
! 03.12.2001    ggu     LEVMX - look out for missing level of near node
! 05.12.2001    ggu     NTBC - BUG -> has not been set before
! 27.03.2003    ggu     in bndo_adjbc use ambient value (bamb)
! 13.03.2004    ggu     in bndo_adjbc only for level BC (LEVELBC)
! 05.10.2004    ggu     new routine bndo_radiat, ibtyp=31
! 31.05.2007    ggu     reset BC for flux to old type (DEBHELP)
! 23.08.2007    ggu     use iopbnd as indicator for ext/int boundary nodes
! 08.04.2008    ggu     file cleaned (new bndo_setbc, bndo_impbc)
! 17.04.2008    ggu     calls to infobnd deleted (subst by get_bnd_ipar)
! 03.09.2008    ggu     new routine bndo_info_file()
! 06.11.2008    ggu     better error handling
! 12.11.2009    ggu     new array itynod and is_zeta_bound()
! 25.03.2011    ggu     bug fix in bndo_impbc() -> ibcold not initialized
!
!***********************************************************************
!-----------------------------------------------------------------------
        module bndo_admin
!-----------------------------------------------------------------------
        contains
!-----------------------------------------------------------------------

	subroutine bndo_init

! sets up bndo data structure

	use bndo
	use bnd_geom
	use geom
	use basin
	use shympi
        use mpi_graph_elem
        use fem_util
        use bnd_admin
        use elems_dealing

	implicit none

	logical bexternal
	logical berror
	integer k,nodes,itype
	integer i,ibc
	integer inext,ilast,knext,klast
	integer ie,n
	integer ii,iii,ib,in,kn,nb,j
	integer nbc
	double precision area
	double precision dx,dy

	integer nelex,e
	logical bmpielems

!----------------------------------------------------------
! set up array iopbnd
!----------------------------------------------------------

	bmpielems = shympi_partition_on_elements()

	do k=1,nkn_local
	  iopbnd(k) = 0
	end do

	nbndo = 0
        ndebug = 0              !unit number for debug (in common block)

	nbc = nbnds()

	do ibc = 1,nbc
	  nodes = nkbnds(ibc)
	  itype = itybnd(ibc)

	  bexternal = ( itype .ge. 1 .and. itype .le. 2 .or. itype .ge. 31 .and. itype .le. 39 )

          if((bexternal).and. (nodes .eq. 1)) cycle

	  do i=1,nodes
	    k = kbnds(ibc,i)
	    if( bexternal ) then
	      nbndo = nbndo + 1
	      if( nbndo .gt. kbcdim ) goto 99
	      iopbnd(k) = nbndo
	      ibcnod(nbndo) = ibc
	      kbcnod(nbndo) = k
	      itynod(nbndo) = itype
	    else
	      iopbnd(k) = -ibc
	    end if
	  end do

	end do

!----------------------------------------------------------
! set up normal direction
!----------------------------------------------------------

	do i=1,nbndo
	  k = kbcnod(i)
	  ibc = ibcnod(i)


          if(bmpielems.and.bounds%last(i).ne.0) then
            knext=bounds%last(i)
            inext=bounds%last(i)
            kantv(1,k)=knext
          else
	    knext = kantv(1,k)
	    inext = iopbnd(knext)
          end if

          if(bmpielems.and.bounds%first(i).ne.0) then
            klast=bounds%first(i)
            ilast=bounds%first(i)
            kantv(2,k)=klast
          else
	    klast = kantv(2,k)
	    ilast = iopbnd(klast)
          end if


!	  -------------------------------
!	  internal consistency check
!	  -------------------------------

!          gmica
!          if(iopbnd(k).eq.1.and.my_id.eq.10)then
!            write(6,*)'inext',inext,knext,ipv(knext),ipv(k)
!            write(6,*)'ilast',ilast,klast,ipv(klast)
!          end if
!
!          if(iopbnd(k).eq.12.and.my_id.eq.9)then
!            write(6,*)'inext',inext,knext,ipv(knext)
!            write(6,*)'ilast',ilast,klast,ipv(klast)
!          end if
!
!          if(iopbnd(k).eq.70.and.my_id.eq.0)then
!            write(6,*)'inext',inext,knext,ipv(knext)
!            write(6,*)'ilast',ilast,klast,ipv(klast)
!          end if

          if(.not.bmpielems) then
	    if( iopbnd(k) .ne. i ) then
	      stop 'internal error bndo (0)'
	    end if
	    if( inext .gt. 0 ) then
	      if( kbcnod(inext) .ne. knext ) then
	        stop 'internal error bndo (1)'
	      end if
	    end if
	    if( ilast .gt. 0 ) then
	      if( kbcnod(ilast) .ne. klast ) then
	        stop 'internal error bndo (2)'
	      end if
	    end if
	  end if

!	  -------------------------------
!	  adjacent boundary nodes must be of same boundary
!	  -------------------------------

          if(.not.bmpielems) then
	    if( inext .gt. 0 ) then
	      if( ibcnod(inext) .ne. ibc ) then
	        goto 98
	      end if
	    end if
	    if( ilast .gt. 0 ) then
	      if( ibcnod(ilast) .ne. ibc ) then
	        goto 98
	      end if
	    end if
	  end if

!	  -------------------------------
!	  get normal direction
!	  -------------------------------

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
	    write(6,*) 'One node open boundary not permitted'
	    write(6,*) i,k,inext,ilast
	    write(6,*) 'node number (external): ',ipext(k)
	    write(6,*) 'boundary: ',ibc,nkbnds(ibc)
	    stop 'error stop bndo'
	  end if

	  xynorm(1,i) = -dy			!x-component
	  xynorm(2,i) = dx			!y-component
	end do

!        gmica
!        call shympi_barrier
!        if(bmpielems) then
!          if(my_id.eq.10 ) write(6,*)'xynorm:',xynorm(:,1)
!          if(my_id.eq.9 )write(6,*)'xynorm:',xynorm(:,12)
!        else
!          write(6,*)'xynorm:',xynorm(:,70)
!        end if

!----------------------------------------------------------
! set up weights and node list
!----------------------------------------------------------

	do i=1,nbndo
	  k = kbcnod(i)
	  if(k.gt.nkn) cycle
	  nopnod(i) = 0
	end do

        if(bmpielems) then
          nelex=nel_local
        else
          nelex=nel
        end if

        wopnodes=0.
        nopnodes=0.

#ifdef DEBUGON                                  !gmica: for reproducibility of sum of wopnodes
        do e=1,nel_local
           if(bmpi) then
              ie=domain%elems%mapID(e)
           else
              ie=e
           end if
#else
	do ie=1,nelex
#endif
	  !call elebase(ie,n,ibase)
	  n = 3
	  area = areaele(ie)

	  do ii=1,n
	    k = nen3v(ii,ie)
	    ib = iopbnd(k)
	    if( ib .gt. 0 ) then		!insert inner nodes
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

!              call shympi_barrier
!        if(bmpielems) then
!          call setup_bnd_neighbors(nbndo,nopnod,nopnodes)
!        end if
!              call shympi_barrier

!----------------------------------------------------------
! scale weights to unit
!----------------------------------------------------------

	berror = .false.

	do i=1,nbndo

	  k = kbcnod(i)
	  if(k.gt.nkn) cycle
	  nb = nopnod(i)
	  if( nb .le. 0 ) then
	    ibc = ibcnod(i)
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
	  
	end do

	if( berror ) stop 'error stop bndo'

	write(6,*) 'finished setting up bndo_init, nbndo = ',nbndo

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	return
   98	continue
	write(6,*) 'different boundary : ',ibc
	stop 'error stop bndo'
   99	continue
	write(6,*) 'dimension error kbcdim : ',kbcdim
	stop 'error stop bndo'
	end

!***********************************************************************

	subroutine bndinsert(ib,area,kn)

! inserts node kn with weight area into list (internal routine)

	use bndo
        use basin
        use shympi

	implicit none

	integer ib		!nodal index
	double precision area		!weight
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

        if(shympi_partition_on_elements()) then
          do j=1,bounds%nneighbors
            if(kn.eq.bounds%neighbors(j)) return
          end do
        end if

        if(my_id.eq.10.and.ib.eq.1) write(6,*)'ib:',kn,ipv(kn)

	nopnodes(nb,ib) = kn
	wopnodes(nb,ib) = area
	nopnod(ib) = nb

	end

!***********************************************************************

	subroutine bndo_info_file(file)

! writes bndo info to file

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

!***********************************************************************

	subroutine bndo_info(iunit)

! writes info on open boundary nodes to terminal

	use bndo
	use bnd_geom
	use levels
	use basin
        use fem_util
        use bnd_admin

	implicit none

        integer iunit

	integer i,k,ibc,nb,j,iu

        iu = iunit
        if( iu .le. 0 ) iu = 6

	write(iu,*) '--------------------------------'
	write(iu,*) 'Information on open boundary nodes'
	write(iu,*) 'Total number of open boundary nodes: ',nbndo

	do i=1,nbndo

	  k = kbcnod(i)
	  if(k.gt.nkn) cycle
	  ibc = ibcnod(i)
	  nb = nopnod(i)

	  if( iopbnd(k) .ne. i ) then
	    stop 'internal error bndo: (11)'
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

!***********************************************************************

	function is_zeta_bound(k)

! checks if node k is a zeta boundary

	use bndo
	use bnd_geom

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

!***********************************************************************

        subroutine bndo_setbc(it,what,nlvddi,cv,rbc,uprv,vprv)

! sets open boundary condition for level boundaries
!
! simply calls bndo_impbc() and bndo_adjbc()

	use basin
        use shympi

        implicit none

        integer it
        character*(*) what      !conz/temp/salt or else
        integer nlvddi
        double precision cv(nlvddi,nkn_local)
        double precision rbc(nlvddi,nkn)	!boundary condition (3D)
	double precision uprv(nlvddi,nkn)
	double precision vprv(nlvddi,nkn)

!----------------------------------------------------------
! simply imposes whatever is in rbc
!----------------------------------------------------------

        call bndo_impbc(it,what,nlvddi,cv,rbc)

!----------------------------------------------------------
! adjusts for ambient value, no gradient or outgoing flow
!----------------------------------------------------------

	call bndo_adjbc(it,what,nlvddi,cv,uprv,vprv)

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	end

!***********************************************************************

        subroutine bndo_impbc(it,what,nlvddi,cv,rbc)

! imposes boundary conditions on open boundary

	use bndo
	use bnd_geom
	use levels
        use bnd_admin
	use basin, only : nkn,nel,ngr,mbw
        use defnames

        implicit none

        integer it
        character*(*) what      !conz/temp/salt or else
        integer nlvddi
        double precision cv(nlvddi,nkn)
        double precision rbc(nlvddi,nkn)	!boundary condition (3D)

        logical bgrad0
        logical bdebug
        integer i,j,k,l
        integer ibc,ibcold
        integer nb,nlev
        integer ibtyp
        double precision value

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

          if( iopbnd(k) .ne. i ) stop 'internal error bndo: (11)'

	  if( ibtyp .eq. 1 ) then
            nlev = ilhkv(k)

            do l=1,nlev
              cv(l,k) = rbc(l,k)
            end do
	  end if

        end do

        end

!***********************************************************************

	subroutine bndo_adjbc(it,what,nlvddi,cv,uprv,vprv)

! adjusts boundary conditions on open boundary (values on bnd already set)
!
! adjusts for ambient value, no gradient or outgoing flow

	use bndo
	use bnd_geom
	use levels
	use basin, only : nkn,nel,ngr,mbw,ipv
        use shympi
        use fem_util
        use bnd_admin
        use defnames

	implicit none

	integer it
        character*(*) what	!conz/temp/salt or else
	integer nlvddi
	double precision cv(nlvddi,nkn_local)
	double precision uprv(nlvddi,nkn)
	double precision vprv(nlvddi,nkn)

	logical bgrad0
	logical blevel
	logical bdebug
	logical bdggu
	logical bout,bamb
	logical binside
	integer i,j,k,l
	integer ibtyp,igrad0
	integer ibc,ibcold
	integer nb,nlev,ko
        integer ntbc,nlevko
	double precision dx,dy
	double precision scal,bc,weight,tweight
	double precision value

	bdebug = .true.
	bdebug = .false.

	ibcold = 0

	bgrad0 = .false.
	blevel = .false.

	if( bdebug ) then
	  if( ndebug .eq. 0 ) then
	    ndebug = ifemopa('bndo_adjbc (91)','.bndo','form','unknown')
            call bndo_info(ndebug)
	  end if
	  write(ndebug,*) 'bndo_adjbc ........... ',what,it
	end if

	do i=1,nbndo

	  k = kbcnod(i)
	  nb = nopnod(i)
	  ibc = ibcnod(i)

	  if( iopbnd(k) .ne. i ) stop 'internal error bndo: (11)'

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

	  do l=1,nlev
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

!***********************************************************************

	subroutine bndo_radiat(it,rzv)

! imposes radiation condition for levels

	use bndo
	use bnd_geom
	use hydro_admin
	use basin, only : nkn,nel,ngr,mbw
        use fem_util
        use bnd_admin
        use defnames

	implicit none

	integer it
	double precision rzv(nkn)

	logical bdebug
	integer i,j,k,l
	integer ibc,ibcold
	integer nb,nlev,ko
        integer ntbc,nlevko
	integer ibtyp
	double precision bc,weight,tweight

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

	  if( iopbnd(k) .ne. i ) stop 'internal error bndo: (11)'

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

!***********************************************************************

!-----------------------------------------------------------------------
        end module bndo_admin
!-----------------------------------------------------------------------
