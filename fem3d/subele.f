
!--------------------------------------------------------------------------
!
!    Copyright (C) 2004-2006,2008,2010-2015,2017,2019  Georg Umgiesser
!    Copyright (C) 2011  Debora Bellafiore
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

c routines for dealing with elements
c
c contents :
c
c function depele(ie,mode)
c       computes total (average) depth of element ie
c subroutine dep3dele(ie,mode,nlev,h)
c       computes (average) depth of element ie for all layers
c subroutine dep3dnod(k,mode,nlev,h)
c       computes depth of node k for all layers
c function zetaele(ie,mode)
c       computes (average) zeta of element ie
c function volele(ie,mode)
c       computes volume of element ie
c function areaele(ie)
c       computes area of element ie
c subroutine nindex(ie,n,kn)
c       returns node index for one element
c subroutine setdepele(ie,h)
c       sets depth for element ie
c subroutine depvele(ie,mode,n,h)
c       computes depth of vertices in element ie
c function depfvol(ie,ii,mode)
c       computes total depth in vertex ii of element ie
c function scalele(ie,sev)
c       returns average of scalar sev in element ie
c subroutine scalvele(ie,sv,n,s)
c       returns scalar value s of vertices in element ie
c subroutine tempvele(ie,n,s)
c subroutine saltvele(ie,n,s)
c subroutine conzvele(ie,n,s)
c subroutine setsvele(ie,sev,s)
c       sets scalar values sev to s in element ie
c subroutine addsvele(ie,sev,s)
c       adds scalar value s to sev in element ie
c subroutine elenmax(text,nmax)
c       checks maximum of vertices of element
c subroutine ele2node(sev,sv)
c       computes nodal values from element values (scalar)
c subroutine elebase(ie,n,ibase)
c       returns number of vertices and pointer into element array
c subroutine set_area
c       sets up area for nodes
c subroutine dvanode(l,k,mode,dep,vol,area)
c       returns depth, volume and area of node k on level l
c function depnode(l,k,mode)
c       returns depth of node k on level l
c function volnode(l,k,mode)
c       returns volume of node k on level l
c function areanode(l,k)
c       returns area of node k on level l
c subroutine setdepth(levdim,hdkn,hden,zenv,area)
c       sets up depth array for nodes
c function scalcont(mode,scal)
c       computes content of scalar in total domain
c function scalcontk(mode,scal,k)
c       computes content of scalar at node k
c
c revision log :
c
c 23.01.2004	ggu	new routines scalcontkh and scalmass
c 23.01.2004	ggu	-> FIXME: scalmass should go to other file
c 22.02.2005	ggu	routines deleted: tempvele saltvele conzvele
c 25.08.2005	ggu	bug in dvanode fixed (use * instead of +)
c 29.11.2006	ggu	in copydepth do not set to 0 old depths
c 07.11.2008	ggu	new helper routine make_new_depth()
c 23.03.2010	ggu	changed v6.1.1
c 16.12.2010	ggu	setdepth() changed for sigma levels
c 25.10.2011	ggu	hlhv eliminated
c 04.11.2011	ggu	adapted for hybrid coordinates
c 08.11.2011	dbf&ggu	bug in setdepth(): 1 -> l
c 11.11.2011	ggu	error message for negative last layer
c 22.11.2011	ggu	changed VERS_6_1_37
c 09.12.2011	ggu	changed VERS_6_1_38
c 30.03.2012	ggu	changed VERS_6_1_51
c 29.03.2013	ggu	avoid call to areaele -> ev(10,ie)
c 13.06.2013	ggu	new helper functions make_old_depth and copy_depth
c 12.09.2013	ggu	changed VERS_6_1_67
c 29.11.2013	ggu	new subroutine masscont()
c 05.12.2013	ggu	changed VERS_6_1_70
c 18.06.2014	ggu	changed VERS_6_1_77
c 13.10.2014	ggu	changed VERS_7_0_2
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 05.06.2015	ggu	changed VERS_7_1_12
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 24.07.2015	ggu	changed VERS_7_1_82
c 23.09.2015	ggu	changed VERS_7_2_4
c 05.11.2015	ggu	changed VERS_7_3_12
c 18.12.2015	ggu	changed VERS_7_3_17
c 05.12.2017	ggu	changed VERS_7_5_39
c 03.02.2019	ggu	in setdepth check for zero layer thickness
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c 01.04.2021	ggu	dimension check in dep3dnod()
c 02.06.2021	ggu	computing depth at node run over nkn_inner
c
c****************************************************************

	function depele(ie,mode)

c computes total (average) depth of element ie

	use mod_hydro
	use basin

	implicit none

	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	real depele	!depth (return)

	integer ii,n
	real hmed

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

c****************************************************************

	subroutine dep3dele(ie,mode,nlev,h)

c computes (average) depth of element ie for all layers

	use mod_layer_thickness
	use levels

	implicit none

	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	integer nlev	!total number of levels
	real h(nlev)	!depth of layers

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

c****************************************************************

	subroutine dep3dnod(k,mode,nlev,h)

c computes depth of node k for all layers

	use mod_layer_thickness
	use levels

	implicit none

	integer k	!number of node
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	integer nlev	!total number of levels
	real h(nlev)	!depth of layers

	integer ndim
	integer l

	ndim = nlev
	nlev = ilhkv(k)
	if( nlev > ndim ) goto 99

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

	return
   99	continue
	write(6,*) 'k,ndim,nlev: ',k,ndim,nlev
	stop 'error stop dep3dnod: nlev>ndim'
	end
	
c****************************************************************

	function zetaele(ie,mode)

c computes (average) zeta of element ie

	use mod_hydro

	implicit none

	real zetaele	!depth (return)
	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta

	integer ii,n
	real hmed

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

c****************************************************************

	function volele(ie,mode)

c computes volume of element ie

	implicit none

	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	real volele	!volume (return)

	real area, depth
	real areaele, depele

	area = areaele(ie)
	depth = depele(ie,mode)

	volele = area * depth

	end

c****************************************************************

	function areaele(ie)

c computes area of element ie

	use evgeom

        implicit none

        integer ie      !number of element
        real areaele    !volume (return)

	areaele = 12. * ev(10,ie)

	end

c****************************************************************

	subroutine nindex(ie,n,kn)

c returns node index for one element

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

c****************************************************************

	subroutine setdepele(ie,h)

c sets depth for element ie

	use basin

	implicit none

	integer ie      !number of element
	real h		!depth to use

	integer ii,n

	!call elebase(ie,n,ibase)
	n = 3

	do ii=1,n
	  hm3v(ii,ie) = h
	end do

	end

c****************************************************************

	subroutine depvele(ie,mode,n,h)

c computes depth of vertices in element ie

	use mod_hydro
	use basin

	implicit none

	integer ie	!number of element
	integer mode	!-1: old zeta   0: no zeta   1: new zeta
	integer n	!total number of vertices (return)
	real h(3)	!depth at vertices (return)

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

c****************************************************************

	function depfvol(ie,ii,mode)

c computes total depth in vertex ii of element ie

	use mod_hydro
	use basin

	implicit none

	real depfvol
	integer ie	!number of element
	integer ii	!number of vertex
	integer mode	!-1: old zeta   0: no zeta   1: new zeta

	integer n
	real h

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

c****************************************************************
c****************************************************************
c****************************************************************

	function scalele(ie,sev)

c returns average of scalar sev in element ie

	implicit none

	real scalele
	integer ie
	real sev(3,*)

	integer ii
	real sm

	sm = 0.
	do ii=1,3
	  sm = sm + sev(ii,ie)
	end do

	scalele = sm / 3.

	end

c****************************************************************

	subroutine scalvele(ie,sv,n,s)

c returns scalar value s of vertices in element ie

	implicit none

	integer ie	!number of element
	real sv(3,*)	!array of scalar value
	integer n	!total number of vertices (return)
	real s(3)	!scalar values at vertices (return)

	integer ii

	n = 3
	do ii=1,3
	  s(ii) = sv(ii,ie)
	end do

	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine setsvele(ie,sev,s)

c sets scalar values sev to s in element ie

	implicit none

	integer ie
	real sev(3,*)
	real s

	integer ii

	do ii=1,3
	  sev(ii,ie) = s
	end do

	end

c****************************************************************

	subroutine addsvele(ie,sev,s)

c adds scalar value s to sev in element ie

	implicit none

	integer ie
	real sev(3,*)
	real s

	integer ii

	do ii=1,3
	  sev(ii,ie) = sev(ii,ie) + s
	end do

	end

c****************************************************************

c***********************************************************

	subroutine elenmax(text,nmax)

c checks maximum of vertices of element

	implicit none

	character*(*) text
	integer nmax

	if( nmax .lt. 3 ) then
	  write(6,*) 'Dimension of vertices too small:'
	  write(6,*) text
	  stop 'error stop elenmax'
	end if

	end

c***********************************************************

	subroutine ele2node(sev,sv)

c computes nodal values from element values (scalar)

	use mod_hydro
	use basin

	implicit none

	real sev(3,nel)
	real sv(nkn)

	integer ie,k,ii,n
	integer ip
	integer kn(10)
	real h(10), s(10)
	real area,vol
	logical btype
	real v1v(nkn)
	
	real areaele

	btype = .false.

	call elenmax('ele2node',10)

	do k=1,nkn
	  sv(k) = 0.
	  v1v(k) = 0.
	end do

	if( btype ) then

	do ie=1,nel

	  area = areaele(ie)
	  call nindex(ie,n,kn)
	  call depvele(ie,+1,n,h)
	  call scalvele(ie,sev,n,s)

	  do ii=1,n
	    k = kn(ii)
	    vol = area * h(ii)
	    sv(k) = sv(k) + s(ii) * vol
	    v1v(k) = v1v(k) + vol
	  end do

	end do

	else

	do ie=1,nel

	  area = areaele(ie)

	  do ii=1,3
	    k = nen3v(ii,ie)
	    vol = area * ( hm3v(ii,ie) + zenv(ii,ie) )
	    sv(k) = sv(k) + sev(ii,ie) * vol
	    v1v(k) = v1v(k) + vol
	  end do

	end do

	end if

	do k=1,nkn
	  sv(k) = sv(k) / v1v(k)
	end do

	end

c***********************************************************

	subroutine elebase(ie,n,ibase)

c returns number of vertices and pointer into element array

	implicit none

	integer ie	!element number
	integer n	!total number of vertices (return)
	integer ibase	!pointer to start of values

	n = 3
	ibase = (ie-1) * 3

	end

c***********************************************************
c***********************************************************
c***********************************************************

	subroutine set_area

c sets up area for nodes

	use levels
	use basin
	use mod_area
	use shympi

	implicit none

	integer k,l,ie,ii
	integer nlev,n
	real areael,areafv
	real areaele

	areakv = 0.

	do ie=1,nel
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

	if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange areakv')
          call shympi_exchange_and_sum_3d_nodes(areakv)
	else
	  !call shympi_comment('exchanging areakv')
	  call shympi_exchange_3d_node(areakv)
	end if

	!call shympi_barrier

	end

c***********************************************************

	subroutine dvanode(l,k,mode,dep,vol,area)

c returns depth, volume and area of node k on level l

	use mod_layer_thickness
	use mod_area

	implicit none

	integer l	!layer
	integer k	!node
	integer mode	!-1: old    1: new 
	real dep	!layer thickness of node k and layer l
	real vol	!volume of node k and layer l
	real area	!area of node k and layer l

	if( mode .gt. 0 ) then
	  dep = hdknv(l,k)
	else if( mode .lt. 0 ) then
	  dep = hdkov(l,k)
	else
	  write(6,*) 'mode = 0 not implemented'
	  stop 'error stop dvanode'
	end if

	area = areakv(l,k)
        vol = dep * area

	end

c***********************************************************

	function depnode(l,k,mode)

c returns depth of node k on level l

	use mod_layer_thickness

	implicit none

	real depnode
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

c***********************************************************

	function volnode(l,k,mode)

c returns volume of node k on level l

	implicit none

	real volnode
	integer l
	integer k
	integer mode

	real depnode,areanode

	volnode = depnode(l,k,mode) * areanode(l,k)

	end

c***********************************************************

	function areanode(l,k)

c returns area of node k on level l

	use mod_area

	implicit none

	real areanode
	integer l
	integer k

	areanode = areakv(l,k)

	end

c***********************************************************

	subroutine copydepth(levdim,hdkn,hdko,hden,hdeo)

c sets up depth array for nodes

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer levdim
	real hdkn(levdim,nkn)	!depth at node, new time level
	real hdko(levdim,nkn)	!depth at node, old time level
	real hden(levdim,nel)	!depth at element, new time level
	real hdeo(levdim,nel)	!depth at element, old time level

	integer k,ie,l

	do k=1,nkn
	  do l=1,levdim
	    hdko(l,k) = hdkn(l,k)
	    !hdkn(l,k) = 0.
	  end do
	end do

	do ie=1,nel
	  do l=1,levdim
	    hdeo(l,ie) = hden(l,ie)
	    !hden(l,ie) = 0.
	  end do
	end do

	end

c***********************************************************

	subroutine copy_depth

c shell (helper) for copydepth

	use mod_layer_thickness
	use levels, only : nlvdi,nlv

	implicit none

	call copydepth(nlvdi,hdknv,hdkov,hdenv,hdeov)

	end

c***********************************************************

	subroutine make_new_depth

c shell (helper) for setdepth

	use mod_layer_thickness
	use mod_area
	use mod_hydro
	use levels, only : nlvdi,nlv

	implicit none

	call setdepth(nlvdi,hdknv,hdenv,zenv,areakv)

	end

c***********************************************************

	subroutine make_old_depth

c shell (helper) for setdepth

	use mod_layer_thickness
	use mod_area
	use mod_hydro
	use levels, only : nlvdi,nlv

	implicit none

	call setdepth(nlvdi,hdkov,hdeov,zeov,areakv)

	end

c***********************************************************

	subroutine check_diff_depth

c checks differences between old and new depth values (debug)

	use mod_layer_thickness
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
	  
c***********************************************************

	subroutine setdepth(levdim,hdkn,hden,zenv,area)

c sets up depth array for nodes

	use mod_depth
	use evgeom
	use levels
	use basin
	use shympi

	implicit none

	integer levdim
	real hdkn(levdim,nkn)	!depth at node, new time level
	real hden(levdim,nel)	!depth at element, new time level
	real zenv(3,nel)    	!water level at new time level
	real area(levdim,nkn)

        logical bdebug
        logical bsigma
	integer k,l,ie,ii
	integer lmax,n,nlev,nsigma,levmin
	real hfirst,hlast,h,htot,z,zmed,hm
	real hacu,hlevel,hsigma,hsig
	real hmin

	real areael,areafv
	real areaele

        bdebug = .false.
	hmin = -99999.
	hmin = 0.

c----------------------------------------------------------------
c initialize and copy
c----------------------------------------------------------------

	hdkn = 0.
	hden = 0.

c----------------------------------------------------------------
c compute volumes at node
c----------------------------------------------------------------

	call get_sigma_info(nlev,nsigma,hsigma)
	bsigma = nsigma .gt. 0
	hfirst = hldv(1)

	do ie=1,nel

	  !call elebase(ie,n,ibase)
	  n = 3
	  areael = 12 * ev(10,ie)
	  areafv = areael / n

	  lmax = ilhv(ie)
	  hm = hev(ie)
	  zmed = 0.

c	  -------------------------------------------------------
c	  nodal values
c	  -------------------------------------------------------

	  do ii=1,n

	    k = nen3v(ii,ie)
	    htot = hm3v(ii,ie)
	    z = zenv(ii,ie)
	    hsig = min(htot,hsigma) + z		!total depth of sigma layers
	    zmed = zmed + z

	    do l=1,nsigma
	      hdkn(l,k) = - hsig * hldv(l)	!these have already depth
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
	        if( levmin .eq. 1 ) hdkn(1,k) = hdkn(1,k) + areafv * z
	        hlast = htot - hlv(lmax-1)
		if( hlast .lt. 0. ) goto 77
	        hdkn(lmax,k) = hdkn(lmax,k) + areafv * hlast
	      end if
	    end if

	    !do l=1,lmax
	    !  if( hdkn(l,k) <= 0. ) then
	    !    write(6,*) 'no volume in node ',k
	    !    write(6,*) 'depth: ',h
	    !    write(6,*) 'nlv,lmax: ',nlv,lmax
	    !    write(6,*) 'nsigma: ',nsigma
	    !    write(6,*) 'hlv: ',hlv
	    !    write(6,*) 'hldv: ',hldv
	    !    write(6,*) 'hdkn: ',hdkn(:,k)
	    !    stop 'error stop setdepth: no volume in node'
	    !  end if
	    !end do

	  end do

!	  in hdkn is volume of finite volume around k
!	  in sigma layers hdkn already contains nodal depth

c	  -------------------------------------------------------
c	  element values
c	  -------------------------------------------------------

	  k = 0
	  ii = 0
	  zmed = zmed / n
	  htot = hev(ie)
	  hsig = min(htot,hsigma) + zmed	!total depth of sigma layers

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

	  do l=1,lmax
	    h = hden(l,ie)
	    if( h <= hmin ) then
	      write(6,*) 'error computing layer thickness'
	      write(6,*) 'no layer depth in element: ',ie,l,lmax
	      write(6,*) 'depth: ',h,htot,zmed
	      write(6,*) 'additional information available in fort.666'
	      call check_set_unit(666)
	      call check_elem(ie)
	      call check_nodes_in_elem(ie)
	      stop 'error stop setdepth: no layer in element'
	    end if
	  end do

	end do

	if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange hdkn')
          call shympi_exchange_and_sum_3d_nodes(hdkn)
	end if

c----------------------------------------------------------------
c compute depth at nodes
c----------------------------------------------------------------

	levmin = nsigma + 1
	do k=1,nkn_inner
	  lmax = ilhkv(k)
	  do l=levmin,lmax
	    areafv = area(l,k)
	    if( areafv .gt. 0. ) then
	      hdkn(l,k) = hdkn(l,k) / areafv
	    else
	      exit
	    end if
	  end do
	  do l=1,lmax
	    h = hdkn(l,k)
	    if( h <= hmin ) then
	      write(6,*) 'error computing layer thickness'
	      write(6,*) 'no layer depth in node: ',k,l,lmax
	      write(6,*) 'depth: ',h
	      write(6,*) 'nlv,lmax: ',nlv,lmax
	      write(6,*) 'nsigma: ',nsigma
	      write(6,*) 'hlv: ',hlv
	      write(6,*) 'hldv: ',hldv
	      write(6,*) 'additional information available in fort.666'
	      call check_set_unit(666)
	      call check_node(k)
	      call check_elems_around_node(k)
	      stop 'error stop setdepth: no layer in node'
	    end if
	  end do
	end do

c----------------------------------------------------------------
c echange nodal values
c----------------------------------------------------------------

	if( shympi_partition_on_nodes() ) then
	  !call shympi_comment('exchanging hdkn')
	  call shympi_exchange_3d_node(hdkn)
	  !call shympi_barrier
	end if

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

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

c***********************************************************

	function masscont(mode)

c computes content of water mass in total domain

	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

        implicit none

	double precision masscont
	integer mode

	integer k,l,nlev
	real mass
	double precision total
        real volnode

!SHYMPI_ELEM - should be total nodes to use - FIXME shympi

	total = 0.

	do k=1,nkn
	  nlev = ilhkv(k)
	  do l=1,nlev
	    total = total + volnode(l,k,mode)
	  end do
	end do

	mass = total
	mass = shympi_sum(mass)
	masscont = mass

	end

c***********************************************************

	function scalcont(mode,scal)

c computes content of scalar in total domain

	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

        implicit none

	double precision scalcont
	integer mode
	real scal(nlvdi,nkn)

	logical, parameter :: bdebug = .false.
	integer k,l,nlev
	real mass
	double precision total
        real volnode

!SHYMPI_ELEM - should be total nodes to use - FIXME shympi

	total = 0.

	do k=1,nkn
	  nlev = ilhkv(k)
	  do l=1,nlev
	    total = total + volnode(l,k,mode) * scal(l,k)
	    if( bdebug .and. scal(l,k) .gt. 0. ) then
		write(6,*) 'scalcont: ',l,k,scal(l,k)
	    end if
	  end do
	end do

	mass = total
	mass = shympi_sum(mass)
	scalcont = mass

	end

c***********************************************************

	function scalcontk(mode,scal,k)
 
c computes content of scalar at node k
 
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        double precision scalcontk
        integer mode,k
        real scal(nlvdi,nkn)
 
        integer l,nlev
        double precision total
        real volnode
 
        total = 0.
 
          nlev = ilhkv(k)
          do l=1,nlev
            total = total + volnode(l,k,mode) * scal(l,k)
          end do
 
        scalcontk = total
 
        end
 
c***********************************************************

	function scalcontkh(scal,k,depth)
 
c computes content of scalar at node k (with given depth)
 
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        double precision scalcontkh
        integer k
        real scal(nlvdi,nkn)
        real depth
 
        integer l,nlev
        double precision total
        real areanode
 
        total = 0.
 
          nlev = ilhkv(k)
          do l=1,nlev
            total = total + depth * areanode(l,k) * scal(l,k)
          end do
 
        scalcontkh = total
 
        end
 
c***********************************************************

        subroutine scalmass(scal,depth,tstot)

c this routine into another file... FIXME

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        real scal(nlvdi,nkn)        !scalar field for which to compute mass
        real depth                      !depth of column
        real tstot                      !total mass (return)

        integer k
        double precision scalcontkh
        double precision total

        total = 0.

        do k=1,nkn
          total = total + scalcontkh(scal,k,depth)
        end do

        tstot = total

        end

c***********************************************************

