c
c $Id: subflx.f,v 1.25 2009-05-21 09:24:00 georg Exp $
c
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
c 30.04.1998    ggu	newly written routines (subpor deleted)
c 07.05.1998    ggu	check nrdveci on return for error
c 08.05.1998    ggu	restructured with new comodity routines
c 13.09.1999    ggu	type of node computed in own routine flxtype
c 19.11.1999    ggu	iskadj into sublin
c 20.01.2000	ggu	old routines substituted, new routine extrsect
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 01.09.2002	ggu	ggu99 -> bug in flx routines (how to reproduce?)
c 26.05.2003	ggu	in flxnov substituted a,b with b,c
c 26.05.2003	ggu	new routine make_fluxes (for lagrangian)
c 10.08.2003	ggu	do not call setweg, setnod, setkan
c 23.03.2006    ggu     changed time step to real
c 28.09.2007    ggu     use testbndo to determine boundary node in flxtype
c 28.04.2009    ggu     links re-structured
c 23.02.2011    ggu     new routine call write_node_fluxes() for special output
c 01.06.2011    ggu     documentation to flxscs() changed
c 21.09.2011    ggu     low-level routines copied from subflxa.f
c 07.10.2011    ggu     implemented 3d flux routines
c 20.10.2011    ggu     restructured, flx3d_k(), make_fluxes_3d()
c 16.12.2011    ggu     bug fix: in make_fluxes_2/3d() r/tflux was integer
c 04.05.2012    ggu     bug fix: in flx3d_k correct for flux boundary
c
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine flx2d_k(k,istype,az,n,transp)

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
	real transp(1)		!fluxes into elements (flux corrected, return) 

	include 'param.h'


	logical bdebug
	integer i,ie,ii,ne,ndim
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
c get pointer into link structure
c---------------------------------------------------------

	call set_elem_links(k,ne)

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
	  ie = lnk_elems(i)
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

	integer ndim		!must be at least ngr
	parameter (ndim=100)

	include 'param.h'

	logical bdebug
	integer i,n
	real tt

	real transp(ndim)
	real weight(ndim)
	real weight1(ndim)

	bdebug = k .eq. 6615
	bdebug = k .eq. 0

c---------------------------------------------------------
c compute transport through finite volume k
c---------------------------------------------------------

	n = ndim
	call flx2d_k(k,istype,az,n,transp)

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

	subroutine flx3d_k(k,istype,az,lkmax,n,transp)

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

	include 'param.h'

	integer k		!node number of finite volume
	integer istype		!type of node (see flxtype)
	real az			!time weighting parameter
	integer lkmax		!maximum layer in finite volume k (return)
	integer n		!dimension/size of transp (entry/return)
	real transp(nlvdi,1)	!computed fluxes (return)

	logical bdebug
	integer i,ie,ii,ne,ndim
	integer l,lmax
	real aj,area,dz,dt
	real uv,uvn,uvo
	real b,c
	real azt
	real div,dvdt,q,qw_top,qw_bot

	real dvol(nlvdi)
	real areal(nlvdi+1)

	integer ithis
	real areanode,volnode

c---------------------------------------------------------
c get parameters
c---------------------------------------------------------

	bdebug = k .eq. 6615
	bdebug = k .eq. 0

	ndim = n
	call get_timestep(dt)
	azt = 1. - az

c---------------------------------------------------------
c get pointer into link structure
c---------------------------------------------------------

	call set_elem_links(k,ne)

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
	  ie = lnk_elems(i)
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

            dvdt = dvol(l)/dt
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
	use evgeom
	use levels, only : nlvdi,nlv
	use basin

	implicit none

	include 'param.h'
	include 'femtime.h'

	integer k		!node number of finite volume
	integer ibefor,iafter	!pointer to pre/post node
	integer istype		!type of node (see flxtype)
	real az			!time weighting parameter
	integer lkmax		!maximum layer in finite volume k (return)
	real flux(nlvdi)	!computed fluxes (return)

	logical bdebug
	integer i,n
	integer l
	real ttot

	real tt(nlvdi)

	real transp(nlvdi,ngr)
	real weight(ngr)
	real weight1(ngr)

	bdebug = k .eq. 6615
	bdebug = k .eq. 0

c---------------------------------------------------------
c compute transport through finite volume k
c---------------------------------------------------------

	n = ngr
	call flx3d_k(k,istype,az,lkmax,n,transp)

	if( bdebug ) then
	  write(88,*) '3d ',k,lkmax,n
	  write(88,*) (transp(1,i),i=1,n)
	end if

c---------------------------------------------------------
c check if transport is really divergence free
c---------------------------------------------------------

	if( istype .le. 2 ) then	!no open boundary
	  do l=1,lkmax
	    ttot = 0.
	    do i=1,n
	      ttot = ttot + transp(l,i)
	    end do
	    if( abs(ttot) .gt. 1. ) then
	      write(6,*) '******** flx3d (divergence): ',ttot
	      write(6,*) '     ',k,l,lkmax,istype
	    end if
	  end do
	end if

	if( bdebug ) then
	  write(88,*) 'ttot ',ttot
	  write(88,*) '---------------------------------'
	end if

c---------------------------------------------------------
c compute transport through section in finite volume k
c---------------------------------------------------------

	do l=1,lkmax
	  tt(l) = 0.
	end do

	call mkweig(n,istype,ibefor,weight)

	do i=1,n
	  do l=1,lkmax
	    tt(l) = tt(l) - weight(i) * transp(l,i)	!this flux is negative
	  end do
	end do

	call mkweig(n,istype,iafter,weight1)

	do i=1,n
	  do l=1,lkmax
	    tt(l) = tt(l) + weight1(i) * transp(l,i)
	  end do
	end do

	do l=1,lkmax
	  flux(l) = tt(l)
	end do

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

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
	real weight(1)

	integer it,i
	real start,fact,dw

	logical debug
	save debug
	data debug /.false./

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

	include 'param.h'

	integer ktype
	integer kafter,kbefor

	integer ipext


	if( kantv(1,k) .eq. 0 ) then			!inner point
	   ktype = 1
	else if( .not. is_external_boundary(k) ) then	!material boundary
	   ktype = 2
	else if( is_external_boundary(k) ) then		!open boundary
	   kafter = kantv(1,k)
	   kbefor = kantv(2,k)
	   if( is_external_boundary(kbefor) 
     +			.and. is_external_boundary(kafter) ) then	!OOO
		ktype = 5
	   else if( is_external_boundary(kbefor) ) then			!OOB
		ktype = 4
	   else if( is_external_boundary(kafter) ) then			!BOO
		ktype = 3
	   else
		write(6,*) 'error at open boundary node ',ipext(k)
		write(6,*) 'bounadry consisting of one node'
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

	include 'param.h'

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

