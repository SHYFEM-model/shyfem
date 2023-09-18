!
! $Id: subele.f,v 1.15 2010-03-11 15:36:39 georg Exp $
!
! routines for dealing with elements
!
! contents :
!
! function depele(ie,mode)
!       computes total (average) depth of element ie
! subroutine dep3dele(ie,mode,nlev,h)
!       computes (average) depth of element ie for all layers
! subroutine dep3dnod(k,mode,nlev,h)
!       computes depth of node k for all layers
! function zetaele(ie,mode)
!       computes (average) zeta of element ie
! function volele(ie,mode)
!       computes volume of element ie
! function areaele(ie)
!       computes area of element ie
! subroutine nindex(ie,n,kn)
!       returns node index for one element
! subroutine setdepele(ie,h)
!       sets depth for element ie
! subroutine depvele(ie,mode,n,h)
!       computes depth of vertices in element ie
! function depfvol(ie,ii,mode)
!       computes total depth in vertex ii of element ie
! function scalele(ie,sev)
!       returns average of scalar sev in element ie
! subroutine scalvele(ie,sv,n,s)
!       returns scalar value s of vertices in element ie
! subroutine tempvele(ie,n,s)
! subroutine saltvele(ie,n,s)
! subroutine conzvele(ie,n,s)
! subroutine setsvele(ie,sev,s)
!       sets scalar values sev to s in element ie
! subroutine addsvele(ie,sev,s)
!       adds scalar value s to sev in element ie
! subroutine elenmax(text,nmax)
!       checks maximum of vertices of element
! subroutine ele2node(sev,sv)
!       computes nodal values from element values (scalar)
! subroutine elebase(ie,n,ibase)
!       returns number of vertices and pointer into element array
! subroutine set_area
!       sets up area for nodes
! subroutine dvanode(l,k,mode,dep,vol,area)
!       returns depth, volume and area of node k on level l
! function depnode(l,k,mode)
!       returns depth of node k on level l
! function volnode(l,k,mode)
!       returns volume of node k on level l
! function areanode(l,k)
!       returns area of node k on level l
! subroutine setdepth(levdim,hdkn,hden,zenv,area)
!       sets up depth array for nodes
! function scalcont(mode,scal)
!       computes content of scalar in total domain
! function scalcontk(mode,scal,k)
!       computes content of scalar at node k
!
! revision log :
!
! 23.01.2004    ggu     new routines scalcontkh and scalmass
! 23.01.2004	ggu	-> FIXME: scalmass should go to other file
! 22.02.2005    ggu     routines deleted: tempvele saltvele conzvele
! 25.08.2005    ggu     bug in dvanode fixed (use * instead of +)
! 29.11.2006    ggu     in copydepth do not set to 0 old depths
! 07.11.2008    ggu     new helper routine make_new_depth()
! 16.12.2010    ggu     setdepth() changed for sigma levels
! 25.10.2011    ggu     hlhv eliminated
! 04.11.2011    ggu     adapted for hybrid coordinates
! 08.11.2011    dbf&ggu bug in setdepth(): 1 -> l
! 11.11.2011    ggu     error message for negative last layer
! 29.03.2013    ggu     avoid call to areaele -> ev(10,ie)
! 13.06.2013    ggu     new helper functions make_old_depth and copy_depth
! 29.11.2013    ggu     new subroutine masscont()
!
!****************************************************************
!----------------------------------------------------------------
        module elems_dealing
!----------------------------------------------------------------
        contains
!----------------------------------------------------------------

	function depele(ie,mode)

! computes total (average) depth of element ie

	use hydro_admin
	use basin

	implicit none

	include 'param.h'

	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	double precision depele	!depth (return)


	integer ii,n
	double precision hmed

	!call elebase(ie,n,ibase)
	n = 3

	hmed = 0.

	if( mode .lt. 0 ) then
	  do ii=1,n
	    hmed = hmed + hm3v(ii,ie) + zeov(ii,ie)
	  end do
	else if( mode .gt. 0 ) then
	  do ii=1,n
	    hmed = hmed + hm3v(ii,ie) + zenv(ii,ie)
	  end do
	else
	  do ii=1,n
	    hmed = hmed + hm3v(ii,ie)
	  end do
	end if

	depele = hmed / n

	end

!****************************************************************

	subroutine dep3dele(ie,mode,nlev,h)

! computes (average) depth of element ie for all layers

	use layer_thickness
	use levels

	implicit none

	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	integer nlev	!total number of levels
	double precision h(nlev)	!depth of layers

	integer l

	nlev = ilhv(ie)

	if( mode .gt. 0 ) then
	  do l=1,nlev
	    h(l) = hdenv(l,ie)
	  end do
	else if( mode .lt. 0 ) then
	  do l=1,nlev
	    h(l) = hdeov(l,ie)
	  end do
	else
	  stop 'error stop dep3dele: Cannot use mode = 0'
	end if

	end

!****************************************************************

	subroutine dep3dnod(k,mode,nlev,h)

! computes depth of node k for all layers

	use layer_thickness
	use levels

	implicit none

	integer k	!number of node
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	integer nlev	!total number of levels
	double precision h(nlev)	!depth of layers

	integer l

	nlev = ilhkv(k)

	if( mode .gt. 0 ) then
	  do l=1,nlev
	    h(l) = hdknv(l,k)
	  end do
	else if( mode .lt. 0 ) then
	  do l=1,nlev
	    h(l) = hdkov(l,k)
	  end do
	else
	  stop 'error stop dep3dnod: Cannot use mode = 0'
	end if

	end
	
!****************************************************************

	function zetaele(ie,mode)

! computes (average) zeta of element ie

	use hydro_admin

	implicit none

	double precision zetaele	!depth (return)
	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta

	integer ii,n
	double precision hmed

	!call elebase(ie,n,ibase)
	n = 3

	hmed = 0.

	if( mode .lt. 0 ) then
	  do ii=1,n
	    hmed = hmed + zeov(ii,ie)
	  end do
	else if( mode .gt. 0 ) then
	  do ii=1,n
	    hmed = hmed + zenv(ii,ie)
	  end do
	else
	  hmed = 0.
	end if

	zetaele = hmed / n

	end

!****************************************************************

	function volele(ie,mode)

! computes volume of element ie

	implicit none

	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	double precision volele	!volume (return)

	double precision areael, depth

	areael = areaele(ie)
	depth = depele(ie,mode)

	volele = areael * depth

	end

!****************************************************************

	function areaele(ie)

! computes area of element ie

	use evgeom

        implicit none

        integer ie      !number of element
        double precision areaele    !volume (return)

	areaele = 12. * ev(10,ie)

	end

!****************************************************************

	subroutine nindex(ie,n,kn)

! returns node index for one element

	use basin

	implicit none

	integer ie      !number of element
	integer n	!total number of nodes in element (3 for triangle)
	integer kn(1)	!node index

	integer ii

	!call elebase(ie,n,ibase)
	n = 3

	do ii=1,n
	  kn(ii) = nen3v(ii,ie)
	end do

	end

!****************************************************************

	subroutine setdepele(ie,h)

! sets depth for element ie

	use basin

	implicit none

	integer ie      !number of element
	double precision h		!depth to use

	integer ii,n

	!call elebase(ie,n,ibase)
	n = 3

	do ii=1,n
	  hm3v(ii,ie) = h
	end do

	end

!****************************************************************

	subroutine depvele(ie,mode,n,h)

! computes depth of vertices in element ie

	use hydro_admin
	use basin

	implicit none

	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	integer n	!total number of vertices (return)
	double precision h(3)	!depth at vertices (return)

	integer ii

	!call elebase(ie,n,ibase)
	n = 3

	if( mode .lt. 0 ) then
	  do ii=1,n
	    h(ii) = hm3v(ii,ie) + zeov(ii,ie)
	  end do
	else if( mode .gt. 0 ) then
	  do ii=1,n
	    h(ii) = hm3v(ii,ie) + zenv(ii,ie)
	  end do
	else
	  do ii=1,n
	    h(ii) = hm3v(ii,ie)
	  end do
	end if

	end

!****************************************************************

	function depfvol(ie,ii,mode)

! computes total depth in vertex ii of element ie

	use hydro_admin
	use basin

	implicit none

	double precision depfvol
	integer ie	!number of element
	integer ii	!number of vertex
	integer mode	!-1: old zeta   0: no zeta   1: new zeta

	integer n
	double precision h

	!call elebase(ie,n,ibase)
	n = 3

	if( mode .lt. 0 ) then
	    h = hm3v(ii,ie) + zeov(ii,ie)
	else if( mode .gt. 0 ) then
	    h = hm3v(ii,ie) + zenv(ii,ie)
	else
	    h = hm3v(ii,ie)
	end if

	depfvol = h

	end

!****************************************************************
!****************************************************************
!****************************************************************

	function scalele(ie,sev)

! returns average of scalar sev in element ie

	implicit none

	double precision scalele
	integer ie
	double precision sev(3,*)

	integer ii
	double precision sm

	sm = 0.
	do ii=1,3
	  sm = sm + sev(ii,ie)
	end do

	scalele = sm / 3.

	end

!****************************************************************

	subroutine scalvele(ie,sv,n,s)

! returns scalar value s of vertices in element ie

	implicit none

	integer ie	!number of element
	double precision sv(3,*)	!array of scalar value
	integer n	!total number of vertices (return)
	double precision s(3)	!scalar values at vertices (return)

	integer ii

	n = 3
	do ii=1,3
	  s(ii) = sv(ii,ie)
	end do

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine setsvele(ie,sev,s)

! sets scalar values sev to s in element ie

	implicit none

	integer ie
	double precision sev(3,*)
	double precision s

	integer ii

	do ii=1,3
	  sev(ii,ie) = s
	end do

	end

!****************************************************************

	subroutine addsvele(ie,sev,s)

! adds scalar value s to sev in element ie

	implicit none

	integer ie
	double precision sev(3,*)
	double precision s

	integer ii

	do ii=1,3
	  sev(ii,ie) = sev(ii,ie) + s
	end do

	end

!****************************************************************

!***********************************************************

	subroutine elenmax(text,nmax)

! checks maximum of vertices of element

	implicit none

	character*(*) text
	integer nmax

	if( nmax .lt. 3 ) then
	  write(6,*) 'Dimension of vertices too small:'
	  write(6,*) text
	  stop 'error stop elenmax'
	end if

	end

!***********************************************************

	subroutine ele2node(sev,sv)

! computes nodal values from element values (scalar)

	use hydro_admin
	use basin

	implicit none

	double precision sev(3,nel)
	double precision sv(nkn)

	integer ie,k,ii,n
	integer ip
	integer kn(10)
	double precision h(10), s(10)
	double precision areael,vol
	logical btype
	double precision v1v(nkn)
	
	btype = .false.

	call elenmax('ele2node',10)

	do k=1,nkn
	  sv(k) = 0.
	  v1v(k) = 0.
	end do

	if( btype ) then

	do ie=1,nel

	  areael = areaele(ie)
	  call nindex(ie,n,kn)
	  call depvele(ie,+1,n,h)
	  call scalvele(ie,sev,n,s)

	  do ii=1,n
	    k = kn(ii)
	    vol = areael * h(ii)
	    sv(k) = sv(k) + s(ii) * vol
	    v1v(k) = v1v(k) + vol
	  end do

	end do

	else

	do ie=1,nel

	  areael = areaele(ie)

	  do ii=1,3
	    k = nen3v(ii,ie)
	    vol = areael * ( hm3v(ii,ie) + zenv(ii,ie) )
	    sv(k) = sv(k) + sev(ii,ie) * vol
	    v1v(k) = v1v(k) + vol
	  end do

	end do

	end if

	do k=1,nkn
	  sv(k) = sv(k) / v1v(k)
	end do

	end

!***********************************************************

	subroutine elebase(ie,n,ibase)

! returns number of vertices and pointer into element array

	implicit none

	integer ie	!element number
	integer n	!total number of vertices (return)
	integer ibase	!pointer to start of values

	n = 3
	ibase = (ie-1) * 3

	end

!***********************************************************
!***********************************************************
!***********************************************************

	subroutine set_area

! sets up area for nodes

	use area
	use levels
	use basin
	use shympi
        use mpi_io_admin

	implicit none

	integer k,l,ie,ii,e
	integer nlev,n
	double precision areael,areafv

        areakv= 0.d0

#ifdef DEBUGON
        call shympi_exchange_halo_2d_elems(ilhv)

        do e=1,nel_local
          if(bmpi) then
            ie=domain%elems%mapID(e)
          else
            ie=e
          end if
#else
	do ie=1,nel
#endif
	  n = 3
	  areael = areaele(ie)
	  areafv = areael / n
	  nlev = ilhv(ie)
	  do ii=1,n
	    k = nen3v(ii,ie)
	    do l=1,nlev
	      areakv(l,k) = areakv(l,k) + areafv
	    end do
	  end do
	end do

        !call shympi_comment('exchanging areakv')
	call shympi_exchange_3d_node(areakv)
	!call shympi_barrier

        !call shympi_comment('shympi_elem: exchange areakv')
#ifndef DEBUGON
        call shympi_exchange_and_sum_3D_nodes(areakv)
#endif

	end

!***********************************************************

	subroutine dvanode(l,k,mode,dep,vol,areak)

! returns depth, volume and area of node k on level l

	use layer_thickness
	use area

	implicit none

	integer l	!layer
	integer k	!node
	integer mode	!-1: old    1: new 
	double precision dep	!layer thickness of node k and layer l
	double precision vol	!volume of node k and layer l
	double precision areak	!area of node k and layer l

	if( mode .gt. 0 ) then
	  dep = hdknv(l,k)
	else if( mode .lt. 0 ) then
	  dep = hdkov(l,k)
	else
	  write(6,*) 'mode = 0 not implemented'
	  stop 'error stop dvanode'
	end if

	areak = areakv(l,k)
        vol = dep * areak

	end

!***********************************************************

	function depnode(l,k,mode)

! returns depth of node k on level l

	use layer_thickness

	implicit none

	double precision depnode
	integer l
	integer k
	integer mode

	if( mode .gt. 0 ) then
	  depnode = hdknv(l,k)
	else if( mode .lt. 0 ) then
	  depnode = hdkov(l,k)
	else
	  write(6,*) 'mode = 0 not implemented'
	  stop 'error stop depnode'
	end if

	end

!***********************************************************

	function volnode(l,k,mode)

! returns volume of node k on level l

	implicit none

	double precision volnode
	integer l
	integer k
	integer mode

	volnode = depnode(l,k,mode) * areanode(l,k)

	end

!***********************************************************

	function areanode(l,k)

! returns area of node k on level l

	use area

	implicit none

	double precision areanode
	integer l
	integer k

	areanode = areakv(l,k)

	end

!***********************************************************

        subroutine copydepth(levdim,hdkn,hdko,hden,hdeo)

! sets up depth array for nodes

        use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer levdim
        double precision hdkn(levdim,nkn)   !depth at node, new time level
        double precision hdko(levdim,nkn)   !depth at node, old time level
        double precision hden(levdim,nel)   !depth at element, new time level
        double precision hdeo(levdim,nel)   !depth at element, old time level

        integer k,ie,l

        do k=1,nkn
          do l=1,levdim
            hdko(l,k) = hdkn(l,k)
            !hdkn(l,k) = 0.d0
          end do
        end do

        do ie=1,nel
          do l=1,levdim
            hdeo(l,ie) = hden(l,ie)
            !hden(l,ie) = 0.d0
          end do
        end do

        end

!***********************************************************

	subroutine copy_depth

! shell (helper) for copydepth

	use layer_thickness
	use levels, only : nlvdi,nlv

	implicit none

	call copydepth(nlvdi,hdknv,hdkov,hdenv,hdeov)

	end

!***********************************************************

	subroutine make_new_depth

! shell (helper) for setdepth

	use layer_thickness
	use area
	use hydro_admin
	use levels, only : nlvdi,nlv

	implicit none

	call setdepth(nlvdi,hdknv,hdenv,zenv,areakv)

	end

!***********************************************************

	subroutine make_old_depth

! shell (helper) for setdepth

	use layer_thickness
	use area
	use hydro_admin
	use levels, only : nlvdi,nlv

	implicit none

	call setdepth(nlvdi,hdkov,hdeov,zeov,areakv)

	end

!***********************************************************

	subroutine check_diff_depth

! checks differences between old and new depth values (debug)

	use layer_thickness
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,ie,l,idiff

	idiff = 0

	do k=1,nkn
	  do l=1,nlvdi
	    if( hdkov(l,k) .ne. hdknv(l,k) ) idiff = idiff + 1
	  end do
	end do

	do ie=1,nel
	  do l=1,nlvdi
	    if( hdeov(l,ie) .ne. hdenv(l,ie) ) idiff = idiff + 1
	  end do
	end do

	write(6,*) 'check_diff_depth: ',idiff

	end
	  
!***********************************************************

        subroutine setdepth(levdim,hdkn,hden,zenv,areak)

! sets up depth array for nodes

        use depth
        use evgeom
        use levels
        use basin
        use shympi
        use sigma
        use mpi_io_admin

        implicit none

        integer levdim
#ifdef DEBUGON
        integer e
        double precision hdkn(levdim,nkn_local)	!depth at node, new time level
        double precision hden(levdim,nel_local)	!depth at element, new time level
        double precision zenv(3,nel_local)    	!water level at new time level
        double precision areak(levdim,nkn_local)
#else
        double precision hdkn(levdim,nkn)	!depth at node, new time level
        double precision hden(levdim,nel)	!depth at element, new time level
        double precision zenv(3,nel)    	!water level at new time level
        double precision areak(levdim,nkn)
#endif

        logical bdebug
        logical bsigma
        integer k,l,ie,ii
        integer lmax,n,nlev,nsigma,levmin
        double precision hfirst,hlast,h,htot,z,zmed,hm
        double precision hacu,hlevel,hsigma,hsig

        double precision areael,areafv

        bdebug = .false.

!----------------------------------------------------------------
! initialize and copy
!----------------------------------------------------------------

        hden = 0.d0
        hdkn = 0.d0

!----------------------------------------------------------------
! compute volumes at node
!----------------------------------------------------------------

        call get_sigma_info(nlev,nsigma,hsigma)
        bsigma = nsigma .gt. 0
        hfirst = hldv(1)

#ifdef DEBUGON
        call shympi_exchange_halo_2d_elems(ilhv)
        call shympi_exchange_halo_2d_elems(hev)
        call shympi_exchange_halo_3d_elems3(zenv)
        call shympi_exchange_halo_3d_elems3(hm3v)

        do e=1,nel_local
          if(bmpi) then
            ie=domain%elems%mapID(e)
          else
            ie=e
          end if
#else
        do ie=1,nel
#endif
          n = 3
          areael = 12 * ev(10,ie)
          areafv = areael / n
          lmax = ilhv(ie)
          hm = hev(ie)
          zmed = 0.d0

!	  -------------------------------------------------------
!	  nodal values
!	  -------------------------------------------------------

          do ii=1,n
            k = nen3v(ii,ie)
            htot = hm3v(ii,ie)
            z = zenv(ii,ie)
            hsig = min(htot,hsigma) + z
            zmed = zmed + z
            do l=1,nsigma
              hdkn(l,k) = - hsig * hldv(l)
            end do
            if( lmax .gt. nsigma ) then
              if( lmax .eq. 1 ) then
                h = htot + z
                hdkn(1,k) = hdkn(1,k) + areafv * h
              else
                levmin = nsigma + 1
                do l=levmin,lmax-1
                  hdkn(l,k) = hdkn(l,k) + areafv * hldv(l)
                end do
                if( levmin .eq. 1 ) hdkn(1,k) = hdkn(1,k) + areafv*z
                hlast = htot - hlv(lmax-1)
                if( hlast .lt. 0. ) goto 77
                hdkn(lmax,k) = hdkn(lmax,k) + areafv * hlast
              end if
            end if
          end do

!         in hdkn is volume of finite volume

!	  -------------------------------------------------------
!	  element values
!	  -------------------------------------------------------

          k = 0
          ii = 0
          zmed = zmed / n
          htot = hev(ie)
          hsig = min(htot,hsigma) + zmed
          do l=1,nsigma
            hden(l,ie) = - hsig * hldv(l)
          end do
          if( lmax .gt. nsigma ) then
            if( lmax .eq. 1 ) then
              hden(1,ie) = htot + zmed
            else
              levmin = nsigma + 1
              do l=levmin,lmax-1
                hden(l,ie) = hldv(l)
              end do
              if( levmin .eq. 1 ) hden(1,ie) = hden(1,ie) + zmed
              hlast = htot - hlv(lmax-1)
              if( hlast .lt. 0. ) goto 77
              hden(lmax,ie) = hlast
            end if
          end if
        end do

        !call shympi_comment('shympi_elem: exchange hdkn')
#ifndef DEBUGON
        call shympi_exchange_and_sum_3D_nodes(hdkn)
#endif

!----------------------------------------------------------------
! compute depth at nodes
!----------------------------------------------------------------

        levmin = nsigma + 1
        do k=1,nkn
          do l=levmin,levdim
            areafv = areak(l,k)
            if( areafv .gt. 0. ) then
              hdkn(l,k) = hdkn(l,k) / areafv
            else
              exit
            end if
          end do
        end do

!----------------------------------------------------------------
! echange nodal values
!----------------------------------------------------------------

        !call shympi_comment('exchanging hdkn')
        call shympi_exchange_3d_node(hdkn)

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

        return
   77	continue
        write(6,*) 'last layer negative: ',hlast
        write(6,*) 'levdim,lmax: ',levdim,lmax
        write(6,*) 'ie,ii,k: ',ie,ii,k
        write(6,*) 'htot,hm: ',htot,hm
        write(6,*) 'hm3v  ',(hm3v(ii,ie),ii=1,3)
        write(6,*) 'hlv   ',(hlv(l),l=1,levdim)
        write(6,*) 'hldv  ',(hldv(l),l=1,levdim)
        stop 'error stop setdepth: last layer negative'
        end

!***********************************************************

	function masscont(mode)

! computes content of water mass in total domain

	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

        implicit none

	double precision masscont
	integer mode

	integer k,l,nlev,ntot
	double precision total

	total = 0.

        if(shympi_partition_on_elements()) then
          ntot = nkn_inner         !SHYMPI_ELEM - should be total nodes to use
        else
	  ntot = nkn      !SHYMPI_ELEM - should be total nodes to use
        end if

	do k=1,ntot
	  nlev = ilhkv(k)
	  do l=1,nlev
	    total = total + volnode(l,k,mode)
	  end do
	end do

        total = shympi_sum(total)

	masscont = total

	end

!***********************************************************

	function scalcont(mode,scal)

! computes content of scalar in total domain

	use levels
	use basin, only : nkn,nel,ngr,mbw
        use shympi

        implicit none

	double precision scalcont
	integer mode
	double precision scal(nlvdi,nkn)

	logical bdebug
	integer k,l,nlev,ntot
	double precision total

	bdebug = .false.
	total = 0.

        if(shympi_partition_on_elements()) then
          ntot = nkn_inner    !SHYMPI_ELEM - should be total nodes to use
        else
	  ntot = nkn      !SHYMPI_ELEM - should be total nodes to use
        end if

	do k=1,ntot
	  nlev = ilhkv(k)
	  do l=1,nlev
	    total = total + volnode(l,k,mode) * scal(l,k)
	    if( bdebug .and. scal(l,k) .gt. 0. ) then
		write(6,*) 'scalcont: ',l,k,scal(l,k)
	    end if
	  end do
	end do

        !call shympi_comment('shympi_sum: scalcont')
	scalcont = shympi_sum(total)

	end

!***********************************************************

	function scalcontk(mode,scal,k)
 
! computes content of scalar at node k
 
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        double precision scalcontk
        integer mode,k
        double precision scal(nlvdi,nkn)
 
        integer l,nlev
        double precision total
 
        total = 0.
 
          nlev = ilhkv(k)
          do l=1,nlev
            total = total + volnode(l,k,mode) * scal(l,k)
          end do
 
        scalcontk = total
 
        end
 
!***********************************************************

	function scalcontkh(scal,k,depth)
 
! computes content of scalar at node k (with given depth)
 
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        double precision scalcontkh
        integer k
        double precision scal(nlvdi,nkn)
        double precision depth
 
        integer l,nlev
        double precision total
 
        total = 0.
 
          nlev = ilhkv(k)
          do l=1,nlev
            total = total + depth * areanode(l,k) * scal(l,k)
          end do
 
        scalcontkh = total
 
        end
 
!***********************************************************

        subroutine scalmass(scal,depth,tstot)

! this routine into another file... FIXME

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        double precision scal(nlvdi,nkn)        !scalar field for which to compute mass
        double precision depth                      !depth of column
        double precision tstot                      !total mass (return)

        integer k
        double precision total

        total = 0.

        do k=1,nkn
          total = total + scalcontkh(scal,k,depth)
        end do

        tstot = total

        end

!***********************************************************

!----------------------------------------------------------------
        end module elems_dealing
!----------------------------------------------------------------
