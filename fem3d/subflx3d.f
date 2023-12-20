
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-2000,2002-2003,2006-2007,2009  Georg Umgiesser
!    Copyright (C) 2011-2015,2018-2020  Georg Umgiesser
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

c subroutines for computing discharge / flux
c
c contents :
c
c function flxnov(k,ibefor,iafter,istype,az)	flux through volume k
c subroutine mkweig(n,istype,is,weight)	 	computes weight
c
c function flxtype(k)				determines type of node k (1-5)
c subroutine make_fluxes(k,itype,n,rflux,tflux)	computes fluxes over sides
c
c revision log :
c
c 30.04.1998	ggu	newly written routines (subpor deleted)
c 07.05.1998	ggu	check nrdveci on return for error
c 08.05.1998	ggu	restructured with new comodity routines
c 13.09.1999	ggu	type of node computed in own routine flxtype
c 19.11.1999	ggu	iskadj into sublin
c 20.01.2000	ggu	old routines substituted, new routine extrsect
c 20.01.2000	ggu	common block /dimdim/ eliminated
c 20.01.2000	ggu	common block /dimdim/ eliminated
c 01.09.2002	ggu	ggu99 -> bug in flx routines (how to reproduce?)
c 26.05.2003	ggu	in flxnov substituted a,b with b,c
c 26.05.2003	ggu	new routine make_fluxes (for lagrangian)
c 10.08.2003	ggu	do not call setweg, setnod, setkan
c 23.03.2006	ggu	changed time step to real
c 28.09.2007	ggu	use testbndo to determine boundary node in flxtype
c 28.04.2009	ggu	links re-structured
c 23.02.2011	ggu	new routine call write_node_fluxes() for special output
c 01.06.2011	ggu	documentation to flxscs() changed
c 21.09.2011	ggu	low-level routines copied from subflxa.f
c 07.10.2011	ggu	implemented 3d flux routines
c 18.10.2011	ggu	changed VERS_6_1_33
c 20.10.2011	ggu	restructured, flx3d_k(), make_fluxes_3d()
c 16.12.2011	ggu	bug fix: in make_fluxes_2/3d() r/tflux was integer
c 24.01.2012	ggu	changed VERS_6_1_41
c 04.05.2012	ggu	bug fix: in flx3d_k correct for flux boundary
c 01.06.2012	ggu	changed VERS_6_1_53
c 13.06.2013	ggu	changed VERS_6_1_65
c 19.12.2014	ggu	changed VERS_7_0_10
c 19.01.2015	ggu	changed VERS_7_1_2
c 19.01.2015	ggu	changed VERS_7_1_3
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 16.12.2015	ggu	changed VERS_7_3_16
c 18.12.2015	ggu	changed VERS_7_3_17
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.02.2019	ggu	changed VERS_7_5_60
c 16.02.2020	ggu	femtime eliminated
c 05.03.2020	ggu	do not print flux divergence
c 21.04.2023	ggu	avoid access of node 0 in flxtype()
c 11.12.2023	ggu	prepared for new 3d flux computation
c 18.12.2023	ggu	finished for new 3d flux computation
c 20.12.2023	ggu	debugged 3d flux computation, still to be cleaned
c
c info :
c
c istype is set in flxtype():
c
c	0	not in domain
c	1	inner node
c	2	material boundary
c	3-5	open boundary
c
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine flx2d_k(k,istype,az,n,transp,ne,elems)

c computes fluxes through finite volume k (2D version)
c
c returns n and flux corrected fluxes transp
c on return transp(i) contains fluxes into finite volume (FV)
c n is total number of internal sections of FV
c ne is the total number of elements attached to the FV
c in case of an internal nodes n == ne
c for a boundary node n == ne + 1
c in any case it is the number of internal sections that is returned
c for a boundary node, transp(n) = 0

	use mod_geom
	use mod_hydro_baro
	use mod_hydro
	use evgeom

	implicit none

	integer k		!node number of finite volume
	real az			!time weighting parameter
	integer istype		!type of node (>1 => boundary node)
	integer n		!dimension/size of transp (entry/return)
	real transp(n)		!fluxes into elements (flux corrected, return) 
	integer ne		!total number of elements in elems
	integer elems(ne)	!elements around k

	logical bdebug
	integer i,ie,ii,ndim
	real aj,area,dz,dt
	real dvdt,div
	real uv,uvn,uvo
	real b,c
	real azt,tt

	integer ithis

c---------------------------------------------------------
c get parameters
c---------------------------------------------------------

	bdebug = k .eq. 6615
	bdebug = k .eq. 0

	ndim = n
	call get_timestep(dt)
	azt = 1. - az

c---------------------------------------------------------
c initialize variables
c---------------------------------------------------------

	n = ne
	if( istype .gt. 1 ) n = n + 1		!boundary
	if( n .gt. ndim ) stop 'error stop flx2d_k: n > ndim'

        transp(n) = 0.                !BUG FIX 29.5.2004

c---------------------------------------------------------
c compute transports into finite volume of node k -> transp
c computed transports are divergence corrected
c---------------------------------------------------------

	do i=1,ne
	  ie = elems(i)
	  ii = ithis(k,ie)
	  b = ev(3+ii,ie)
	  c = ev(6+ii,ie)
	  aj = ev(10,ie)
	  area = 4. * aj

	  uvn = unv(ie) * b + vnv(ie) * c
	  uvo = uov(ie) * b + vov(ie) * c
	  uv = 12. * aj * ( az * uvn + azt * uvo )

	  dz = zenv(ii,ie) - zeov(ii,ie)
	  dvdt = dz * area / dt

	  div = -dvdt
	  transp(i) = uv + div

	  !transp(i) = uv - dz * area * rdt
	  !if( bdebug ) write(88,*) i,1,uv
	  !if( bdebug ) write(88,*) i,dvdt
	end do

	end

c******************************************************************

	subroutine flx2d(k,ibefor,iafter,istype,az,flux)

c computes flux through section of finite volume k (2D version)
c
c internal section is defined by:  kbefor - k - kafter
c passed in are pointers to these section in lnk structure

	use mod_geom
	use evgeom

	implicit none

	integer k		!node number of finite volume
	integer ibefor,iafter	!pointer to pre/post node
	integer istype		!type of node (see flxtype)
	real az			!time weighting parameter
	real flux		!flux computed (return)

	logical bdebug
	integer i,n,ne
	double precision tt

	integer elems(maxlnk)
	real transp(maxlnk)
	double precision weight(maxlnk)
	double precision weight1(maxlnk)

	bdebug = k .eq. 6615
	bdebug = k .eq. 0

c---------------------------------------------------------
c compute transport through finite volume k
c---------------------------------------------------------

	n = maxlnk
	call get_elems_around(k,maxlnk,ne,elems)
	call flx2d_k(k,istype,az,n,transp,ne,elems)

	if( bdebug ) then
	  write(88,*) '2d ',k,n
	  write(88,*) (transp(i),i=1,n)
	  write(88,*) '---------------------------------'
	end if

c---------------------------------------------------------
c compute transport through section in finite volume k
c---------------------------------------------------------

	tt = 0.

	call mkweig(n,istype,ibefor,weight)

	do i=1,n
	  tt = tt - weight(i) * transp(i)	!this flux has negative sign
	end do

	call mkweig(n,istype,iafter,weight1)

	do i=1,n
	  tt = tt + weight1(i) * transp(i)
	end do

	flux = tt

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end
	
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine flx3d_k(k,istype,azr,lkmax,n,transp,ne,elems)

c computes fluxes through finite volume k (3D version)
c
c if we are on a flux boundary this is not working
c we should exclude mfluxv from the divergence computation
c -> has been done - other sources of mfluxv (rain, etc.) are also eliminated

	use mod_bound_geom
	use mod_geom
	use mod_bound_dynamic
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels

	implicit none

	integer k		!node number of finite volume
	integer istype		!type of node (see flxtype)
	real azr		!time weighting parameter
	integer lkmax		!maximum layer in finite volume k (return)
	integer n		!dimension/size of transp (entry/return)
	real transp(nlvdi,n)	!computed fluxes (return)
	integer ne		!total number of elements around k
	integer elems(ne)	!elements around k

	logical bdebug
	integer i,ie,ii,ndim
	integer l,lmax
	real dt
	double precision aj,area
	double precision ddt
	double precision uv,uvn,uvo
	double precision b,c
	double precision az,azt
	double precision div,dvdt,q,qw_top,qw_bot

	double precision dvol(nlvdi)
	double precision areal(nlvdi+1)

	integer ithis
	real areanode,volnode

c---------------------------------------------------------
c get parameters
c---------------------------------------------------------

	bdebug = k .eq. 6615
	bdebug = k .eq. 0

	ndim = n
	call get_timestep(dt)

	ddt = dt
	az = azr
	azt = 1. - az

c---------------------------------------------------------
c initialize variables
c---------------------------------------------------------

	n = ne
	if( istype .gt. 1 ) n = n + 1		!boundary
	if( n .gt. ndim ) stop 'error stop flx3d_k: n > ndim'

	lkmax = ilhkv(k)
	do l=1,lkmax
	  areal(l) = areanode(l,k)
	  dvol(l) = volnode(l,k,+1) - volnode(l,k,-1)
          transp(l,n) = 0.				!if on boundary
	end do
	areal(lkmax+1) = 0.

c---------------------------------------------------------
c compute transports into finite volume of node k -> transp
c computed transports are divergence corrected
c note: lkmax is always greater or equal than lmax
c---------------------------------------------------------

	do i=1,ne
	  ie = elems(i)
	  ii = ithis(k,ie)
	  b = ev(3+ii,ie)
	  c = ev(6+ii,ie)
	  aj = ev(10,ie)
	  area = 4. * aj
	  lmax = ilhv(ie)
	  do l=1,lmax
	    uvn = utlnv(l,ie) * b + vtlnv(l,ie) * c
	    uvo = utlov(l,ie) * b + vtlov(l,ie) * c
	    uv = 12. * aj * ( az * uvn + azt * uvo )

            dvdt = dvol(l)/ddt
	    q = mfluxv(l,k)
	    if( is_external_boundary(k) ) q = 0.
	    qw_top = area * wlnv(l-1,k)
	    qw_bot = area * wlnv(l,k)
	    if( l .eq. lmax ) qw_bot = 0.

	    div = (area/areal(l))*(q-dvdt) + qw_bot - qw_top
	    transp(l,i) = uv + div

	    !if( bdebug ) write(88,*) i,l,uv
	    !if( bdebug ) write(88,*) i,dvdt,q,area,areal(l),qw_bot,qw_top
	  end do
	  do l=lmax+1,lkmax
	    transp(l,i) = 0.
	  end do
	end do

	end

c******************************************************************

	subroutine flx3d(k,ibefor,iafter,istype,az,lkmax,flux)

c computes flux through section of finite volume k (3D version)
c
c internal section is defined by:  kbefor - k - kafter
c passed in are pointers to these section in lnk structure

	use mod_geom
	use mod_debug
	use evgeom
	use levels, only : nlvdi,nlv
	use basin

	implicit none

	integer k		!node number of finite volume
	integer ibefor,iafter	!pointer to pre/post node
	integer istype		!type of node (see flxtype)
	real az			!time weighting parameter
	integer lkmax		!maximum layer in finite volume k (return)
	real flux(nlvdi)	!computed fluxes (return)

	logical bdebug,bcheck
	logical bw
	logical bdiverg		!write error on flux divergence
	logical, parameter ::  bold = .true.
	integer i,n,ne,nepres
	integer l

	integer elems(maxlnk)
	real transp(nlvdi,ngr)

	double precision dtransp(nlvdi,ngr)

	double precision ttot,tabs,tlmax,ttmax,rtot,trmax
	double precision dflux,drflux,dabs
	double precision, save :: tmax = 0.
	double precision, save :: tomax = 0.
	double precision, save :: dmax = 0.
	double precision, save :: eps = 1.e-4
	double precision, save :: deps = 1.e-2
	double precision, save :: reps = 1.e-3
	double precision, save :: epsdiv = 1.e-2
	!double precision, save :: epsdiv = 1.e-1

	double precision tt(nlvdi)
	double precision ttnew(nlvdi)

	double precision fflux(maxlnk),qflux(maxlnk)
	double precision weight(ngr)
	double precision weight1(ngr)

	integer ipext

	bdebug = k .eq. 6615
	bdebug = k .eq. 0
	bw = .true.
	bw = .false.

	bdiverg = .false.
	bdiverg = .true.
	bcheck = .true.

c---------------------------------------------------------
c compute transport through finite volume k
c---------------------------------------------------------

	n = maxlnk
	call get_elems_around(k,maxlnk,ne,elems)
	call flx3d_k(k,istype,az,lkmax,n,transp,ne,elems)

	if( n - ne > 1 ) stop 'error stop flx3d: n-ne>1'
	if( n - ne < 0 ) stop 'error stop flx3d: n-ne<0'
	if( n /= ne ) elems(n) = 0

	if( bdebug ) then
	  write(88,*) '3d ',k,lkmax,n
	  write(88,*) (transp(1,i),i=1,n)
	end if

c---------------------------------------------------------
c check if transport is really divergence free
c---------------------------------------------------------

	dtransp = transp

	if( istype .le. 2 ) then	!no open boundary
	  tlmax = 0.
	  ttmax = 0.
	  do l=1,lkmax
	    ttot = 0.
	    tabs = 0.
	    do i=1,n
	      ttot = ttot + dtransp(l,i)
	      tabs = tabs + abs(dtransp(l,i))
	    end do
	    rtot = 0.
	    if( tabs > 0 ) rtot = abs(ttot)/tabs
	    if( bdiverg ) then
	     if( rtot .gt. epsdiv .and. abs(ttot) > 1. ) then
	      write(6,*) '*** flx3d (divergence): '
	      write(6,*) k,l,n,lkmax,istype
	      write(6,*) ttot,tabs,tlmax,rtot
	     end if
	    end if
	    tlmax = max(tlmax,abs(ttot))
	    ttmax = max(ttmax,abs(tabs))
	  end do
	  !if( ttmax > tomax ) then
	  !  tomax = ttmax
	  !  if( bw ) write(6,*) 'tomax: ',k,istype,tomax
	  !end if
	  if( tlmax > tmax ) then
	    tmax = tlmax
	    trmax = 0.
	    if( ttmax > 0. ) trmax = tlmax / ttmax
	    if( bw ) write(6,*) 'tmax: ',k,istype,real(tmax),real(trmax)
	  end if
	end if

	if( bdebug ) then
	  write(88,*) 'ttot ',ttot
	  write(88,*) '---------------------------------'
	end if

c---------------------------------------------------------
c compute transport through section in finite volume k
c---------------------------------------------------------

	tt(1:lkmax) = 0.

	if( bold ) then

	  call mkweig(n,istype,ibefor,weight)

	  do i=1,n
	  do l=1,lkmax
	    tt(l) = tt(l) - weight(i) * dtransp(l,i)	!this flux is negative
	  end do
	  end do

	  call mkweig(n,istype,iafter,weight1)

	  do i=1,n
	  do l=1,lkmax
	    tt(l) = tt(l) + weight1(i) * dtransp(l,i)
	  end do
	  end do

	end if

	  do l=1,lkmax
	    fflux(:) = dtransp(l,:)
	    call flux_with_bounds(n,l,istype,elems,fflux,qflux)
	    ttnew(l) = 0.
	    if( iafter > 0 ) ttnew(l) = ttnew(l) + qflux(iafter)
	    if( ibefor > 0 ) ttnew(l) = ttnew(l) - qflux(ibefor)
	  end do

	!if( bdebug ) then
	if( .false. ) then
	write(665,*) '--------'
	write(665,*) k,ipext(k),lkmax,istype
	write(665,*) ibefor,iafter
	do l=1,lkmax
	write(665,*) l,ttnew(l),tt(l)
	end do
	end if

	if( bcheck ) then	!checks difference between old nad new flux
	do l=1,lkmax
	  dflux = abs( tt(l) - ttnew(l) )
	  !if( dflux /= 0. .and. istype /= 1 ) then
	  !if( dflux > eps .and. istype /= 1 ) then
	  if( dflux > dmax ) then
	    dmax = dflux
	    if( bw ) write(6,*) 'dmax ',k,l,istype,dflux,dmax
	  end if
	  call count_present_elements(l,ne,elems,nepres)
	  if( dflux > deps .and. ne == nepres ) then
	    dabs = 0.5 * (abs(tt(l))+abs(ttnew(l)))
	    drflux = 0.
	    if( dabs > 0 ) drflux = dflux/dabs
	    if( drflux > reps ) then
	    write(6,*) '*** dflux ----------'
	    write(6,*) k,ipext(k),l,n,nepres
	    write(6,*) istype,ibefor,iafter
	    write(6,*) dflux,drflux,dmax
	    write(6,*) tt(l),ttnew(l)
	    end if
	  end if
	end do
	end if

	!flux(1:lkmax) = tt(1:lkmax)
	flux(1:lkmax) = ttnew(1:lkmax)

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end
	
c******************************************************************

	subroutine flux_with_bounds(n,l,istype,elems,fflux,qflux)

	use basin
	use levels
	use mod_geom
	use mod_debug

	implicit none

	integer n,l,istype
	integer elems(n)
	double precision fflux(n)
	double precision qflux(n)

	logical bdebug
	integer npres,i,ie,iu
	integer j,jback,jforw,jnext,jstart
	logical present(n)
	double precision fact,facum,faver

	npres = 0
	present = .false.
	qflux = 0.
	jstart = 0

	faver = 0.
	do i=1,n
	  ie = elems(i)
	  if( ie <= 0 ) cycle
	  if( l > ilhv(ie) ) cycle
	  present(i) = .true.
	  faver = faver + fflux(i)
	  npres = npres + 1
	end do
	if( npres > 0 ) faver = faver / npres

	if( istype <= 2 ) then
	  where( present ) fflux = fflux - faver
	end if

	if( istype > 2 ) then		!open boundary
	  if( istype == 3 ) then	!BOO
	    qflux(n) = 0.
	    jforw = n+1
	    jback = n-1
	  else if( istype == 4 ) then	!OOB
	    qflux(1) = 0.
	    jforw = 2
	    jback = 0
	  else if( istype == 5 ) then	!OOO
	    j = n/2
	    if( 2*j == n ) then	!even
	      qflux(j) = -0.5*fflux(j)
	      qflux(j+1) = 0.5*fflux(j)
	      jforw = j+2
	      jback = j-1
	    else
	      j = j + 1
	      qflux(j) = 0.
	      jforw = j+1
	      jback = j-1
	    end if
	  else
	    write(6,*) 'istype = ',istype
	    stop 'error stop flux_with_bounds: wrong istype'
	  end if
	  do j=jforw,n
	    qflux(j) = qflux(j-1) + fflux(j-1)
	  end do
	  do j=jback,1,-1
	    qflux(j) = qflux(j+1) - fflux(j)
	  end do
	else				!internal or boundary element
	  if( istype == 1 .and. npres == n ) then  !all elements are present
	    fact = 1.
	    facum = 0.
	    do i=n-1,1,-1
	      facum = facum + fact*fflux(i)
	      fact = fact + 1.
	    end do
	    qflux(1) = -facum / n
	    do i=2,n
	      qflux(i) = qflux(i-1) + fflux(i-1)
	    end do
	  else				!missing elements or material boundary
	    if( istype == 2 ) then
	      jstart = n
	    else
	      do i=1,n
	        if( .not. present(i) ) exit
	      end do
	      if( i > n ) stop 'error stop flux_with_bounds: i>n'
	      jstart = i
	    end if
	    j = jstart
	    do i=1,n
	      jnext = 1 + mod(j,n)
	      if( present(j) ) then
	        qflux(jnext) = qflux(j) + fflux(j)
	      else
	        qflux(j) = 0.
	        qflux(jnext) = 0.
	      end if
	      j = jnext
	    end do
	    qflux(jstart) = 0.
	  end if
	end if

	return

	bdebug = l > 51 .and. is_debug()
	iu = 662
	iu = 665

	!if( bdebug .and. any( qflux(1:n) /= 0. ) ) then
	if( l == 1 ) then
	  write(iu,*) ' ========================= qflux: ',istype
	  write(iu,*) n,l,istype,npres,n-npres
	  write(iu,*) jstart
	  write(iu,*) qflux(1:n)
	  write(iu,*) fflux(1:n)
	  write(iu,*) present(1:n)
	end if

	end

c******************************************************************

	subroutine count_present_elements(l,ne,elems,nepres)

	use levels

	implicit none

	integer l
	integer ne
	integer elems(ne)
	integer nepres

	integer i,ie

	nepres = 0
	do i=1,ne
	  ie = elems(i)
	  if( ie <= 0 ) stop 'count_present_elements: ie==0'
	  if( l > ilhv(ie) ) cycle
	  nepres = nepres + 1
	end do

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine mkweig(n,istype,is,weight)

c computes weight over internal section in one finite volume
c
c n		dimension (grade of node)
c istype	type of node (inner, border, ...)
c is		number of internal section to compute
c weigth	weights (return)

	implicit none

	integer n,istype,is
	double precision weight(n)

	integer it,i
	double precision start,fact,dw

	logical, save :: debug = .false.

	it = is

	do i=1,n
	  weight(i) = 0.
	end do

	if( is .eq. 0 ) then		!no section given
c	  nothing
	else if( istype .eq. 1 ) then
	  start = -0.5 * (n-1)
	  fact = 1./n
	  do i=1,n
	    weight(it) = start * fact
	    it = it + 1
	    if( it .gt. n ) it = 1
	    start = start + 1.
	  end do
	else if( istype .eq. 2 ) then
	  if( it .eq. 1 ) it = n
	  do i=it,n-1
	    weight(i) = -1.
	  end do
	else if( istype .eq. 3 ) then
	  do i=it,n-1
	    weight(i) = -1.
	  end do
	else if( istype .eq. 4 ) then
	  do i=1,it-1
	    weight(i) = +1.
	  end do
	else if( istype .eq. 5 ) then
	  if( n .eq. 2 ) then	!handle exception
	    weight(1) = -0.5 + (it-1)		! -0.5/+0.5
	  else
	    fact = 1./(n-1)
	    dw = fact * ( 1. + 1./(n-2) )
	    start = 1.
	    do i=1,n-1
	      weight(i) = (i-1) * dw
	      if( i .ge. it ) then	!upper right triangle
	        weight(i) = weight(i) - start
	      end if
	    end do
	  end if
	else
	  write(6,*) n,istype,is
	  stop 'error stop mkweig : internal error (1)'
	end if
	  
	if( debug .and. istype .ge. 3 ) then
	  write(6,*) n,istype,is
	  write(6,*) (weight(i),i=1,n)
	end if

	end

c******************************************************************

	function flxtype(k)

c determines type of node k (1-5)
c
c uses inodv only for determination of open bounday nodes
c else uses kantv

	use mod_bound_geom
	use mod_geom
	use mod_geom_dynamic

	implicit none

	integer flxtype		!type of node (1=int,2=bnd,3=BOO,4=OOB,5=OOO)
	integer k		!node number of finite volume

	logical bbefor,bafter
	integer ktype
	integer kafter,kbefor

	integer ipext

	flxtype = 0
	if( k <= 0 ) return

	if( kantv(1,k) .eq. 0 ) then			!inner point
	   ktype = 1
	else if( .not. is_external_boundary(k) ) then	!material boundary
	   ktype = 2
	else if( is_external_boundary(k) ) then		!open boundary
	   kafter = kantv(1,k)
	   kbefor = kantv(2,k)
	   bbefor = .false.
	   bafter = .false.
	   if( kbefor > 0 ) bbefor = is_external_boundary(kbefor)
	   if( kafter > 0 ) bafter = is_external_boundary(kafter)
	   if( bbefor .and. bafter ) then		!OOO
		ktype = 5
	   else if( bbefor ) then			!OOB
		ktype = 4
	   else if( bafter ) then			!BOO
		ktype = 3
	   else
		write(6,*) 'error at open boundary node ',ipext(k)
		write(6,*) 'boundary consisting of one node'
		stop 'error stop flxtype'
	   end if
	else
	   stop 'error stop flxtype: internal error (1)'
	end if

	flxtype = ktype

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine make_fluxes_3d(k,itype,lkmax,n,rflux,tflux)

c computes fluxes over sides (tflux) from fluxes into node (rflux)
c
c 3d version
c
c if on boundary (itype>1) rflux(n) is not used (because not defined)

	use levels, only : nlvdi,nlv
	use basin

	implicit none

	integer k		!node
	integer itype		!type of node (1=int,2=bnd,3=BOO,4=OOB,5=OOO)
	integer lkmax		!number of layers
	integer n		!number of sides (tfluxes)
	real rflux(nlvdi,n)	!fluxes into node (element)
	real tflux(nlvdi,n)	!fluxes through sides (return value)

	integer i,l
	real rf(ngr)
	real tf(ngr)

	do l=1,lkmax

	  do i=1,n
	    rf(i) = rflux(l,i)
	  end do

	  call make_fluxes_2d(k,itype,n,rf,tf)

	  do i=1,n
	    tflux(l,i) = tf(i)
	  end do

	end do

	end

c**********************************************************************

	subroutine make_fluxes_2d(k,itype,n,rflux,tflux)

c computes fluxes over sides (tflux) from fluxes into node (rflux)
c
c if on boundary (itype>1) rflux(n) is not used (because not defined)

	implicit none

	integer k		!node
	integer itype		!type of node (1=int,2=bnd,3=BOO,4=OOB,5=OOO)
	integer n		!number of sides (tfluxes)
	real rflux(n)	!fluxes into node (element)
	real tflux(n)	!fluxes through sides (return value)

	integer i
	real rr

	!write(6,*) 'make_fluxes_2d ',k,itype,n
	!write(6,*) 'make_fluxes_2d ',(rflux(i),i=1,n)

	if( itype .eq. 1 ) then		!internal node
		rr = 0.
		do i=1,n-1
		  rr = rr + i * rflux(i)
		end do
		rr = rr / n

		tflux(n) = rr
		do i=n-1,1,-1
		  tflux(i) = tflux(i+1) - rflux(i)
		end do
	else if( itype .eq. 2 ) then	!node on material boundary
		tflux(1) = 0.
		do i=2,n-1
		  tflux(i) = tflux(i-1) + rflux(i-1)
		  !write(6,*) 'make_fluxes_2d ',i,tflux(i)
		end do
		tflux(n) = 0.
	else if( itype .eq. 3 ) then	!BOO - boundary on left
		tflux(n) = 0.
		do i=n-1,1,-1
		  tflux(i) = tflux(i+1) - rflux(i)
		end do
	else if( itype .eq. 4 ) then	!OOB - boundary on right
		tflux(1) = 0.
		do i=2,n
		  tflux(i) = tflux(i-1) + rflux(i-1)
		end do
	else if( itype .eq. 5 ) then	!node totaly on open boundary
		rr = 0.
		do i=1,n-1
		  rr = rr + i * rflux(i)
		end do
		rr = rr / n

		tflux(n) = rr
		do i=n-1,1,-1
		  tflux(i) = tflux(i+1) - rflux(i)
		end do
	else
		stop 'error stop make_fluxes_2d: internal error (1)'
	end if

	end

c**********************************************************************

