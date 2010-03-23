c
c $Id: flux.f,v 1.8 2009-05-21 09:24:00 georg Exp $
c
c flux routines
c
c revision log :
c
c 23.03.2006    ggu     changed time step to real
c 28.04.2009    ggu     links re-structured
c
c*******************************************************************
c
c reads prepared data from unit 47
c writes output data to unit 49
c
c*******************************************************************

	subroutine setflux

	implicit none

	integer ndim
	parameter (ndim = 100)

c	integer nxdim,nydim
c	parameter(nxdim=100,nydim=100)
c	parameter(nxdim=140,nydim=170)

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real r(ndim)		!diverg. corrected fluxes out of element
	real t(ndim)		!original fluxes into element
	real a(ndim)		!area of element
	real z(ndim)		!dz/dt in element
	real sflux(ndim)	!computed fluxes in finite volume
	integer ielem(2,ndim)	!pointer to elements and internal node numbers

	real eflux(7,3,neldim)		!fluxes for all elements

	include 'regflux.h'

	integer k,ne
	integer ktype
	integer zero
	real az,azpar
	real dt

	integer flxtype
	real getpar

	integer idtflx,itmflx
	save idtflx,itmflx

	integer icall
	save icall
	data icall / 0 /

c	-------------------------------------------------------------
c	first call and time check
c	-------------------------------------------------------------

	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
	  if( nint(getpar('regflx')) .ne. 0 ) then
	    icall = 1
	  else
	    icall = -1
	  end if

	  idtflx = nint(getpar('idtflx'))
	  itmflx = nint(getpar('itmflx'))
	  if( idtflx .le. 0 ) icall = -1
	  if( itmflx .ge. itend ) icall = -1

	  if( icall .eq. -1 ) return
	  write(6,*) 'computing regular fluxes...'

	  call rdsind		!reads index from Susanne
	end if

	if( it .lt. itmflx ) return

c	-------------------------------------------------------------
c	dimension check and initialization of variables
c	-------------------------------------------------------------

	if( ngr .gt. ndim ) stop 'error stop setflux: ndim'

	zero = 0
	call getaz(azpar)
	az = azpar
	call get_timestep(dt)

c	-------------------------------------------------------------
c	write info to GRD file
c	-------------------------------------------------------------

c	call flxwrtype		!writes ktype to ktype.grd file

c	-------------------------------------------------------------
c	loop on nodes to compute fluxes over element sides -> eflux
c	-------------------------------------------------------------

	do k=1,nkn
	  ktype = flxtype(k)
	  call flxelem(k,ne,ielem,dt,az,r,t,a,z)
	  call flxnode(k,ktype,ne,r,sflux)
	  call flxnchk(k,ktype,ne,r,sflux)
	  call flxstore(k,ktype,ne,ielem,sflux,t,a,z,eflux)
	end do

c	-------------------------------------------------------------
c	compute missing fluxes internal to elements
c	-------------------------------------------------------------

	call flxinner(nel,eflux)

c	-------------------------------------------------------------
c	do some checks on computed fluxes
c	-------------------------------------------------------------

	call flxnchk(zero,ktype,ne,r,sflux)
	call flxcheck(eflux,nel,dt)

c	-------------------------------------------------------------
c	fill global array and writes to file (unit 49)
c	-------------------------------------------------------------

	call flxfill(it,nel,eflux)

c	-------------------------------------------------------------
c	end of routine
c	-------------------------------------------------------------

	end

c*******************************************************************

	subroutine flxnode(k,ktype,ne,r,sflux)

c computes fluxes through all sections of finite volume

	implicit none

	integer k		!node
	integer ktype		!type of node(internal, boundary, ...)
	integer ne		!total number of elements attached to node
	real r(1)		!divergence corrected fluxes out of element
	real sflux(1)		!flux through sections (return) [m**3/s]

	integer i,ngrade
	real rtot

c all fluxes have unit [m**3/s]
c
c r(ie) = ( dz/dt ) * ( A(ie)/3 ) - T(ie)
c
c where 
c
c T is flux into finite volume through element ie
c A(ie) is area of element ie
c dz/dt is water level rise
c
c in element 1  the total inflow is: S1 - S2 + T1
c in element ne the total inflow is: S6 - S1 + T6
c
c on return ne+1 sections are set
c for inner point sflux(ne+1) = sflux(1)

c if node is internal		=> ngrade = ne
c if node is on boundary	=> ngrade = ne + 1

	ngrade = ne
	if( ktype .gt. 1 ) ngrade = ne + 1

c-------------------------------------------------------------
c handle different cases
c-------------------------------------------------------------

	if( ktype .eq. 1 ) then		!internal node
	  rtot = 0.
	  do i=1,ne-1
	    rtot = rtot + i * r(i)
	  end do
	  sflux(ne) = -rtot / ne
	  do i=ne-1,1,-1
	    sflux(i) = r(i) + sflux(i+1)
	  end do
	  sflux(ne+1) = sflux(1)		!these sections are the same
	else if( ktype .eq. 2 ) then		!node on material boundary
	  sflux(ne+1) = 0.
	  sflux(1) = 0.
	  do i=2,ne
	    sflux(i) = sflux(i-1) - r(i-1)
	  end do
	else if( ktype .eq. 3 ) then		!open boundary node: BOO
	  sflux(ne+1) = 0.
	  do i=ne,1,-1
	    sflux(i) = sflux(i+1) + r(i)
	  end do
	else if( ktype .eq. 4 ) then		!open boundary node: OOB
	  sflux(1) = 0.
	  do i=2,ne+1
	    sflux(i) = sflux(i-1) - r(i-1)
	  end do
	else if( ktype .eq. 5 ) then		!open boundary node: OOO
	  rtot = 0.
	  do i=1,ne
	    rtot = rtot + i * r(i)
	  end do
	  sflux(ne+1) = -rtot / (ne+1)
	  do i=ne,1,-1
	    sflux(i) = r(i) + sflux(i+1)
	  end do
	else
	  write(6,*) 'ktype: ',ktype
	  stop 'error stop flxnode: unknown type'
	end if

c	if( k .eq. 266 .or. k .eq. 257 ) then
c	  write(6,*) 'flxnode: ',k,ne
c	  write(6,*) (r(i),i=1,ne)
c	  write(6,*) (sflux(i),i=1,ne+1)
c	end if

	end

c*******************************************************************

	subroutine flxelem(k,ne,ielem,dt,az,r,t,a,z)

c computes fluxes through all elements of finite volume

	implicit none

	integer k		!node (finite volume)
	integer ne		!total number of elements at node (return)
	integer ielem(2,1)	!index of element numbers (return)
	real dt			!time step
	real az			!time weighting for continuity equation
	real r(1)		!div. corrected fluxes out of element (return)
	real t(1)		!original fluxes into element (return)
	real a(1)		!area of element (return)
	real z(1)		!dz/dt in element

        real uov(1),vov(1),unv(1),vnv(1)
        common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv
        real zenv(3,1), zeov(3,1)
        common /zenv/zenv, /zeov/zeov
	include 'ev.h'
	include 'links.h'

	integer i,ie,ii
	real aj,area,dz,uv,rdt
	real aa,bb
	real azt

	integer ithis

c set some constants

	rdt = 1./dt
	azt = 1. - az

c get pointer into link structure

	call set_elem_links(k,ne)

c compute transports out of finite volume k -> r
c computed transports are divergence corrected

	do i=1,ne
	  ie = lnk_elems(i)
	  ii = ithis(k,ie)
	  ielem(1,i) = ii
	  ielem(2,i) = ie
	  aj = ev(10,ie)
	  area = 4. * aj		!area/3
	  dz = zenv(ii,ie) - zeov(ii,ie)
	  aa = ev(3+ii,ie)
	  bb = ev(6+ii,ie)
	  uv = az * ( unv(ie) * aa + vnv(ie) * bb )
	  uv = uv + azt * ( uov(ie) * aa + vov(ie) * bb )
	  uv = 12. * aj * uv		!uv is flux into k through ie
	  r(i) = dz * rdt * area - uv	!r is diverg. corrected flux out of k
	  t(i) = uv
	  a(i) = area
	  z(i) = dz * rdt
	  if( ie .eq. -1  ) then
	    write(6,*) 'flxelem: ',ie,ii,i
	    write(6,*) r(i), t(i),a(i),z(i)
	  end if
	end do

	end

c*******************************************************************

	subroutine flxstore(k,ktype,ne,ielem,sflux,t,a,z,eflux)

c stores computed fluxes in array

	implicit none

	integer k		!node (finite volume)
	integer ktype		!type of node(internal, boundary, ...)
	integer ne		!total number of elements at node (return)
	integer ielem(2,1)	!index of element numbers (return)
	real sflux(1)		!flux through sections (return) [m**3/s]
	real t(1)		!original fluxes into element 
	real a(1)		!area of element
	real z(1)		!dz/dt in element
	real eflux(7,3,1)	!flux through element sections

	integer nen3v(3,1)
	common /nen3v/nen3v

	integer i,ie,ii
	integer inext

	do i=1,ne
	  ii = ielem(1,i)
	  ie = ielem(2,i)
	  if( nen3v(ii,ie) .ne. k ) then
	    write(6,'(6i8)') ii,ie,k,i,ne,nen3v(ii,ie)
	    stop 'error stop flxstore: internal error (1)'
	  end if
	  eflux(1,ii,ie) = sflux(i)
	  eflux(2,ii,ie) = - sflux(i+1)
	  eflux(3,ii,ie) = t(i)
	  eflux(4,ii,ie) = z(i)
	  eflux(5,ii,ie) = a(i)

	  if( ie .eq. -1 ) then
	    write(6,*) 'flxstore: ',k,ie,ii
	    write(6,*) eflux(1,ii,ie) + eflux(2,ii,ie) + eflux(3,ii,ie)
	    write(6,*) eflux(4,ii,ie) * eflux(5,ii,ie)
	  end if
	end do

	end

c*******************************************************************

	subroutine flxcheck(eflux,nel,dt)

c checks consistency of array eflux

	implicit none

	real eflux(7,3,1)	!flux through element sections
	integer nel		!total number of elements
	real dt			!time step

	include 'ev.h'
        real zenv(3,1), zeov(3,1)
        common /zenv/zenv, /zeov/zeov
	integer nen3v(3,1)
	common /nen3v/nen3v

	integer ie,ii,k,iii
	real fiitot(3)
	real fiitot1(3)
	real fiitot2(3)
	real daz
	real area,dz,rdt
	real ftot,ztot
	real err,tot,rerr
	real rin,rtin,rttin
	integer ilast

	real erra,errt

	integer flxtype

	rdt = 1. / dt
	erra = 0.
	errt = 0.

	do ie=1,nel
	  area = 4. * ev(10,ie)		!area/3
	  ftot = 0.
	  ztot = 0.
	  rin = 0.
	  rtin = 0.
	  rttin = 0.
	  do ii=1,3
	    ilast = mod(ii+1,3) + 1
	    daz = eflux(4,ii,ie) * eflux(5,ii,ie)
	    fiitot(ii) = eflux(1,ii,ie) + eflux(2,ii,ie) + eflux(3,ii,ie)
	    fiitot(ii) = fiitot(ii) - daz
	    ftot = ftot + eflux(1,ii,ie) + eflux(2,ii,ie)
	    dz = zenv(ii,ie) - zeov(ii,ie)
	    ztot = ztot + dz
	    rin = rin + eflux(3,ii,ie)
	    rtin = rtin + eflux(6,ii,ie)
	    rttin = rttin + eflux(3,ii,ie) - eflux(6,ii,ie) 
     +				+ eflux(6,ilast,ie)
	    fiitot1(ii) = eflux(1,ii,ie) + eflux(6,ii,ie) - eflux(7,ii,ie)
	    fiitot1(ii) = fiitot1(ii) - daz * 0.5
	    fiitot2(ii) = eflux(7,ii,ie) + eflux(2,ii,ie) - eflux(6,ilast,ie)
	    fiitot2(ii) = fiitot2(ii) - daz * 0.5
	  end do
	  ztot = area * ztot * rdt
	  err = abs(ftot-ztot)
	  erra = max(erra,err)
	  do ii=1,3
	    errt = max(errt,fiitot(ii))
	  end do

c	  write(67,'(i8,4e14.6)') ie,err,(fiitot(ii),ii=1,3)
	  if( ie .eq. -1 ) then
		write(6,*) 'flxcheck... '
		write(6,*) ie,err,rin,rtin,rttin
	        write(6,'(i8,4e14.6)') ie,err,(fiitot(ii),ii=1,3)
	        write(6,'(i8,4e14.6)') ie,err,(fiitot1(ii),ii=1,3)
	        write(6,'(i8,4e14.6)') ie,err,(fiitot2(ii),ii=1,3)
	  	do ii=1,3
	  	  do iii=1,6
	   	   write(6,*) ii,iii,eflux(iii,ii,ie)
	  	  end do
		end do
	  end if

	end do

	write(68,*) 'flxcheck: ',erra,errt

	end

c*******************************************************************

	subroutine flxwrtype

c writes ktype to GRD file

	implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	integer nen3v(3,1)
	common /nen3v/nen3v

	integer k,ktype
	integer ie,ii
	integer inext,knext,ilast
	integer flxtype
	real xc,yc
	real xh,yh
	real xhal(3),yhal(3)
	real x(3),y(3)
	real a1,a2

	real areat

	integer icall
	save icall
	data icall / 0 /

	if( icall .ne. 0 ) return

	icall = 1

	open(15,file='ktype.grd',status='unknown',form='formatted')
	do k=1,nkn
	  ktype = flxtype(k)
	  write(15,'(i1,2i8,2f15.4)') 1,k,ktype,xgv(k),ygv(k)
	end do
	close(15)

	open(15,file='elemc.grd',status='unknown',form='formatted')
	do ie=1,nel
	  xc = 0.
	  yc = 0.
	  do ii=1,3
	    inext = mod(ii,3) + 1
	    k = nen3v(ii,ie)
	    knext = nen3v(inext,ie)
	    xc = xc + xgv(k)
	    yc = yc + ygv(k)
	    xh = 0.5 * ( xgv(k) + xgv(knext) )
	    yh = 0.5 * ( ygv(k) + ygv(knext) )
	    xhal(ii) = xh
	    yhal(ii) = yh
	    x(ii) = xgv(k)
	    y(ii) = ygv(k)
	    write(15,'(i1,2i8,2f15.4)') 1,nel+3*ie+ii,67,xh,yh
	  end do
	  xc = xc / 3.
	  yc = yc / 3.
	  write(15,'(i1,2i8,2f15.4)') 1,ie,66,xc,yc
	end do
	close(15)

	end

c*******************************************************************

	subroutine flxnchk(k,ktype,ne,r,sflux)

c checks consistency of fluxes for finite volume

	implicit none

	integer k		!node
	integer ktype		!type of node(internal, boundary, ...)
	integer ne		!total number of elements attached to node
	real r(1)		!divergence corrected fluxes out of element
	real sflux(1)		!flux through sections (return) [m**3/s]

	integer i
	real rtot,stot,ttot,aux

	real rtota,stota,ttota
	save rtota,stota,ttota
	data rtota,stota,ttota /0.,0.,0./

	if( k .le. 0 ) then	!write and initialize
c	  write(6,'(a,3e16.8)')  ' flxnchk: ',rtota,stota,ttota
	  write(66,'(a,3e16.8)') ' flxnchk: ',rtota,stota,ttota
	  rtota = 0.
	  rtota = 0.
	  rtota = 0.
	  return
	end if

	rtot = 0.
	stot = 0.
	ttot = 0.
	do i=1,ne
	  rtot = rtot + r(i)
	  stot = stot + sflux(i)
	  aux = sflux(i) - sflux(i+1) - r(i)
	  ttot = max(ttot,abs(aux))
	end do
	rtot = abs(rtot)
	stot = abs(stot)

	if( ktype .gt. 1 ) then
	  stot = 0.
	  if( ktype .gt. 2 ) rtot = 0.
	end if

	rtota = max(rtota,rtot)
	stota = max(stota,stot)
	ttota = max(ttota,ttot)

c	write(77,'(a,2i7,3e16.8)') ' flxnchk: ',k,ktype,rtot,stot,ttot

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine flxfill(it,nel,eflux)

c fills global array with fluxes

	implicit none

	integer it
	integer nel
	real eflux(7,3,1)		!flux through element sections

	include 'regflux.h'

	integer ndim
	parameter(ndim=100)

        real ar2(4,ndim)
        real alen(ndim)
        real aa(4,ndim)
        integer ipm(3,ndim)
        integer ipam(2,ndim)

	integer node(0:nxdim,0:nydim)	!FIXME

	logical bstop,bdebug
	integer nx,ny
	integer mmax,nbox,nsect
	integer i0,j0,i1,j1
	integer ie,id
	real x0,y0,dx,dy
	real tx(3),ty(3)
	real area

	integer i,j

	integer icall
	save icall
	data icall / 0 /

	bdebug = .true.
	bdebug = .false.

c	---------------------------------------------------
c	read header
c	---------------------------------------------------

	call rf47(nx,ny,x0,y0,dx,dy)

	if( bdebug ) then
          write(6,*) 'flxfill...'
          write(6,*) 'nxdim,nydim: ',nxdim,nydim
          write(6,*) 'nx,ny:       ',nx,ny
          write(6,*) 'nel:         ',nel
	end if

	if( nx .gt. nxdim .or. ny .gt. nydim ) then
          write(6,*) 'nxdim,nydim: ',nxdim,nydim
          write(6,*) 'nx,ny:       ',nx,ny
          write(6,*) 'nel:         ',nel
	  stop 'error stop flxfill: nxdim,nydim'
	end if

c	---------------------------------------------------
c	initialize global matrices
c	---------------------------------------------------

	do j=1,ny
	  do i=1,nx
	    hmed(i,j) = 0.
	    areav(i,j) = 0.
	    zlev(i,j) = 0.
	  end do
	end do

	do j=1,ny
	  do i=1,nx
	    aflux(i,j) = 0.
	    icount(i,j) = 0
	  end do
	end do

	do j=0,ny
	  do i=0,nx
	    fflux(i,j,1) = 0.
	    fflux(i,j,2) = 0.
	    alength(i,j,1) = 0.
	    alength(i,j,2) = 0.
	  end do
	end do

c	---------------------------------------------------
c	start of loop reading elements
c	---------------------------------------------------

    1	continue

	  call rd47(ndim,mmax,nbox,nsect,i0,j0,i1,j1,ie,id
     +  	        ,tx,ty,ar2,alen,aa,ipm,ipam,bstop)

	  if( bstop ) goto 2
	  if( id .eq. 0 ) goto 1	!only for subelements

	  if( i1 .gt. nx .or. j1 .gt. ny ) goto 99

	  call addmx(mmax,nbox,nsect,i0,j0,i1,j1,ie,id
     +  	        ,ar2,alen,aa,ipm,ipam,eflux)

	  call addsmx(mmax,nbox,nsect,i0,j0,i1,j1,ie,id
     +  	        ,ar2,aa,ipm,ipam)

	  goto 1
    2	continue

c	---------------------------------------------------
c	check mass conservation in regular grid
c	---------------------------------------------------

	call chckmx(nx,ny)

c	---------------------------------------------------
c	adjust scalar values
c	---------------------------------------------------

	do j=1,ny
	  do i=1,nx
	    area = areav(i,j)
	    if( area .gt. 0. ) then
	      zlev(i,j) = zlev(i,j) / area
	      hmed(i,j) = hmed(i,j) / area
	    end if
	  end do
	end do

c	---------------------------------------------------
c	test some values
c	---------------------------------------------------

	if( icall .eq. 0 ) then
	  write(6,*) 'writing side index...'
	  call wrsidereg('side.grd',nx,ny,nxdim,nydim
     +			,areav,alength,node
     +                  ,x0,y0,dx,dy)
	end if

c	---------------------------------------------------
c	write results to file
c	---------------------------------------------------

	if( icall .eq. 0 ) then			!header, only done once
		call side2dxdy(nx,ny)
		call wf49(nx,ny)
	end if
	call wr49(it,nx,ny)			!regular data, results

	call alexb(it,nx,ny)			!for Alessandro Bergamasco

c	---------------------------------------------------
c	end of routine
c	---------------------------------------------------

	icall = 1

	return
   99	continue
	write(6,*) 'i1,j1: ',i1,j1
	write(6,*) 'nx,ny: ',nx,ny
	write(6,*) 'nxdim,nydim: ',nxdim,nydim
	stop 'error stop flxfill: i1,j1'
	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine addsmx(mmax,nbox,nsect,i0,j0,i1,j1,ie,id
     +  	        ,ar2,aa,ipm,ipam)

c fills global array with scalars of one element

	implicit none

	integer mmax,nbox,nsect
	integer i0,j0,i1,j1
	integer ie,id
        real ar2(4,1)
        real aa(4,1)
        integer ipm(3,1)
        integer ipam(2,1)

	include 'regflux.h'

        real zenv(3,1)
        common /zenv/zenv
        real hm3v(3,1)
        common /hm3v/hm3v

	integer i,j,k,ii
	integer ix,iy
	real area,z,h

c	--------------------------------------------------------
c	initialize variables
c	--------------------------------------------------------

	ii = (id+1)/2

	if( ii .le. 0 ) return

	z = zenv(ii,ie)
	h = hm3v(ii,ie)

c	--------------------------------------------------------
c	store computed scalar variables in global array
c	--------------------------------------------------------

        do j=1,nbox
          ix = i0 + ipam(1,j)
          iy = j0 + ipam(2,j)

	  area = aa(4,j)

          zlev(ix,iy) = zlev(ix,iy) + area * z
          hmed(ix,iy) = hmed(ix,iy) + area * h
          areav(ix,iy) = areav(ix,iy) + area
        end do

c	--------------------------------------------------------
c	end of routine
c	--------------------------------------------------------

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine addmx(mmax,nbox,nsect,i0,j0,i1,j1,ie,id
     +  	        ,ar2,alen,aa,ipm,ipam,eflux)

c fills global array with fluxes of one element

	implicit none

	integer mmax,nbox,nsect
	integer i0,j0,i1,j1
	integer ie,id
        real ar2(4,1)
	real alen(1)
        real aa(4,1)
        integer ipm(3,1)
        integer ipam(2,1)

	real eflux(7,3,1)		!flux through element sections

	include 'regflux.h'

	integer ndim,mdim
	parameter(ndim=50,mdim=10)

	integer i,j,k
	integer ix,iy
	integer iespec,idspec
	real acu
	real area,dzdt
	real r(5)
	real asect(ndim)

c	--------------------------------------------------------
c	initialize variables and do dimension check
c	--------------------------------------------------------

	iespec = 392
	idspec = 6
	iespec = -1
	idspec = -1

	if( ndim .lt. nsect ) goto 99
	if( ndim .lt. nbox ) goto 99
        if( mmax .gt. mdim ) goto 99

c	--------------------------------------------------------
c	extract fluxes of sub-element and put into r
c	--------------------------------------------------------

	call selects(ie,id,eflux,r)

c	--------------------------------------------------------
c	multiply matrix with r -> fluxes are in asect
c	--------------------------------------------------------

        do j=1,nsect
            acu = 0.
            do i=1,4
              acu = acu + ar2(i,j) * r(i)
            end do
            asect(j) = acu
        end do

c	--------------------------------------------------------
c	store computed fluxes in global matrix
c	--------------------------------------------------------

        do j=1,nsect
          k  =      ipm(1,j)
          ix = i0 + ipm(2,j)
          iy = j0 + ipm(3,j)
          fflux(ix,iy,k) = fflux(ix,iy,k) + asect(j)
          alength(ix,iy,k) = alength(ix,iy,k) + alen(j)
        end do

c	--------------------------------------------------------
c	store computed scalar variables in global array
c	--------------------------------------------------------

        do j=1,nbox
          ix = i0 + ipam(1,j)
          iy = j0 + ipam(2,j)

c	if( ix .eq. 52 .and. iy .eq. 8 ) then
c	  write(6,*) ix,iy,ie,id
c	end if
c	if( ix .eq. 51 .and. iy .eq. 8 ) then
c	  write(6,*) ix,iy,ie,id
c	end if

	  area = aa(4,j)
	  dzdt = r(4)

          aflux(ix,iy) = aflux(ix,iy) + area * dzdt
	  icount(ix,iy) = icount(ix,iy) + 1
        end do

c	--------------------------------------------------------
c	debug output
c	--------------------------------------------------------

	if( ie. eq. iespec .and. id .eq. idspec ) then
	  write(6,*) 'addmx... : ',ie,id,nsect,nbox
	  write(6,*) i0,j0
	  write(6,*) (r(i),i=1,5)
	  write(6,*) r(1) + r(2) + r(3) , r(4) * r(5)
	  do j=1,nsect
	    write(6,*) ipm(1,j),ipm(2,j),ipm(3,j),asect(j)
	  end do
	  do j=1,nbox
	    write(6,*) ipam(1,j),ipam(2,j),aa(4,j),aa(4,j)*r(4)
	  end do
	end if

c	--------------------------------------------------------
c	end of routine
c	--------------------------------------------------------

	return
   99	continue
	write(6,*) 'ie,id:           ',ie,id
	write(6,*) 'i0,j0,i1,j1:     ',i0,j0,i1,j1
	write(6,*) 'mmax,mdim:       ',mmax,mdim
	write(6,*) 'nsect,nbox,ndim: ',nsect,nbox,ndim
        stop 'error stop addmx: dimensions'
	end

c**********************************************************************

	subroutine selects(ie,id,eflux,r)

c stores fluxes in sections

	implicit none

	integer ie,id
	real eflux(7,3,1)
	real r(5)

	logical beven
	integer ii
	real fhalf

	ii = (id+1)/2			!node
	beven = (id/2)*2 .eq. id

	fhalf = eflux(7,ii,ie)

	if( beven ) then
	  r(1) = eflux(6,ii,ie)
	  r(2) = -fhalf
	  r(3) = eflux(1,ii,ie)
	else
	  r(1) = +fhalf
	  r(2) = eflux(3,ii,ie) - eflux(6,ii,ie)
	  r(3) = eflux(2,ii,ie)
	end if

	r(4) = eflux(4,ii,ie)	!dz/dt
	r(5) = 0.5 * eflux(5,ii,ie)	!area

	end

c**********************************************************************

	subroutine chckmx(nx,ny)

c checks mass conservation on regular grid

	implicit none

	integer nx,ny

	include 'regflux.h'

	integer i,j,l
	integer ic
	integer itot
	integer nxx
	real fact,az
	real fl,fr,fb,ft
	real ftot,ftotal
	character*130 line,iline
	character*12 chars
	save chars
	data chars /'.0123456789*'/

	fact = 10.
	ftotal = 0.
	itot = 0
	nxx = min(130,nx)

	do j=ny,1,-1
	  line = ' '
	  iline = ' '
	  do i=1,nxx

	    ic = icount(i,j)
	    az = aflux(i,j)

	    fl = fflux(i-1,j,2)
	    fr = fflux(i,j,2)
	    fb = fflux(i,j-1,1)
	    ft = fflux(i,j,1)

	    ftot = fl - fr + fb - ft - az
	    if( j .ne. 1 ) ftotal = max(ftotal,abs(ftot))

	    if( ic .gt. 0 ) itot = itot + 1

	    l = 2 + ic
	    if( l .eq. 2 ) l = 1
	    if( l .gt. 12 ) l = 12
	    iline(i:i) = chars(l:l)

	    l = 2 + nint(fact*abs(ftot))
	    if( ftot .eq. 0. ) l = 1
	    if( l .gt. 12 ) l = 12
	    line(i:i) = chars(l:l)

c	    if( ic .gt. 0 ) then
c	      write(73,*) i,j,l,ic,ftot
c	    end if
	  end do
c	  write(71,'(a)') line(1:nxx)
c	  write(72,'(a)') iline(1:nxx)
	end do

	write(74,*) 'chckmx: ',ftotal,itot,itot/float(nx*ny)

	end

c**********************************************************************

	subroutine flxinner(nel,eflux)

c stores computed fluxes in array

	implicit none

	integer nel
	real eflux(7,3,1)	!flux through element sections

	integer ie,ii
	integer inext
	real tt,drittl

	drittl = 1. / 3.

	do ie=1,nel
	  do ii=1,3
	    inext = mod(ii,3) + 1
	    tt = eflux(3,ii,ie) - eflux(3,inext,ie) 
	    eflux(6,ii,ie) = drittl * tt
	    eflux(7,ii,ie) = eflux(1,ii,ie) + eflux(6,ii,ie)
     +		    - eflux(4,ii,ie) * eflux(5,ii,ie) / 2.
	  end do
	end do

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

        subroutine wrsidereg(file,nx,ny,nxdim,nydim
     +			,am,alength,node
     +                  ,x0reg,y0reg,dx,dy)

c writes grid file

        implicit none

        character*(*) file
        integer nx,ny,nxdim,nydim
        real am(nxdim,nydim)
        real alength(0:nxdim,0:nydim,2)
        integer node(0:nxdim,0:nydim)
        real x0reg,y0reg,dx,dy

        integer ntot,nside
        integer iarea
        integer ix,iy
        real x,y

c open file

        open(1,file=file,status='unknown',form='formatted')

c create nodes

        do iy=0,ny
          do ix=0,nx
            node(ix,iy) = 0
          end do
        end do

        do iy=1,ny
          do ix=1,nx
            if( am(ix,iy) .gt. 0. ) then
                node(ix,iy) = -1
                node(ix-1,iy) = -1
                node(ix,iy-1) = -1
                node(ix-1,iy-1) = -1
            end if
          end do
        end do

c number and write nodes

        ntot = 0
        do iy=0,ny
          do ix=0,nx
            if( node(ix,iy) .eq. -1 ) then
              ntot = ntot + 1
              node(ix,iy) = ntot
              x = x0reg + ix * dx
              y = y0reg + iy * dy
              write(1,1) ntot,0,x,y
            end if
          end do
        end do

c write sides

        nside = 0

c ... horizontal

        do iy=1,ny
          do ix=1,nx
            if( alength(ix,iy,1) .gt. 0. ) then
                nside = nside + 1
                iarea = nint(100.* alength(ix,iy,1)/dx)
                write(1,2) nside,iarea,2
     +                  ,node(ix-1,iy)
     +                  ,node(ix,iy)
            end if
          end do
        end do

c ... vertical

        do iy=1,ny
          do ix=1,nx
            if( alength(ix,iy,2) .gt. 0. ) then
                nside = nside + 1
                iarea = nint(100.* alength(ix,iy,2)/dy)
                write(1,2) nside,iarea,2
     +                  ,node(ix,iy-1)
     +                  ,node(ix,iy)
            end if
          end do
        end do

c ... done

        close(1)

        return
    1   format('1 ',i8,i5,2e14.6)
    2   format('3 ',i8,2i5,2i8)
        end

c**************************************************************

	subroutine alexb(it,nx,ny)

c for Alessandro Bergamasco
c
c---------------------------------------------------------------------------
c
c Format:
c
c Header (written once)
c
c        write(iunit,*) nx,ny				!dimension
c        call wrmx(iunit,nxdim,nydim,nx,ny,hmed)	!depth of boxes
c        call wrmx(iunit,nxdim,nydim,nx,ny,areav)	!area of boxes
c
c Data (written at every output step)
c
c        write(iunit,*) it				!time
c        call wrmx(iunit,nxdim,nydim,nx,ny,zlev)	!water level in boxes
c        call wrmx(iunit,nxdim,nydim,nx,ny,uaux)	!x-velocity in boxes
c        call wrmx(iunit,nxdim,nydim,nx,ny,vaux)	!y-velocity in boxes
c
c Matrix format (subroutine wrmx) -> written line by line
c
c        write(iunit,*) ((matrix(ix,iy),ix=1,nx),iy=1,ny)
c
c---------------------------------------------------------------------------

	implicit none

	integer it
	integer nx,ny

	include 'regflux.h'

	real uaux(nxdim,nydim)
	real vaux(nxdim,nydim)

	integer ix,iy
	integer iunit
	real area,h
	real ul,ur,vb,vt
	real sl,sr,sb,st

	integer icall
	save icall
	data icall / 0 /

	iunit = 51

	if( mod(it,3600) .ne. 0 ) return
	if( it .lt. 0 ) return
c	if( it .gt. 86400 ) return
c	if( it .le. 86400 ) return
c	if( it .ge. 0 ) return			!do not write

	if( icall .eq. 0 ) then
	  icall = 1
	  write(6,*) 'writing header for alexb...'
	  write(iunit,*) nx,ny
	  call wrmx(iunit,nxdim,nydim,nx,ny,hmed)
	  call wrmx(iunit,nxdim,nydim,nx,ny,areav)
	end if

c	ix = 92
c	iy = 126
c	write(6,*) '**********************'
c	write(6,*) ix,iy,hmed(ix,iy),areav(ix,iy)
c	write(6,*) '**********************'

	do iy=1,ny
	  do ix=1,nx
	    area = areav(ix,iy)
	    if( area .gt. 0. ) then

	      h = hmed(ix,iy) + zlev(ix,iy)
	      if( h .le. 0. ) then
		write(6,*) it,ix,iy
		write(6,*) areav(ix,iy),hmed(ix,iy),zlev(ix,iy)
		write(6,*) alength(ix-1,iy,2),alength(ix,iy,2)
		write(6,*) alength(ix,iy-1,1),alength(ix,iy,1)
		write(6,*) fflux(ix-1,iy,2),fflux(ix,iy,2)
		write(6,*) fflux(ix,iy-1,1),fflux(ix,iy,1)
		stop 'error stop alexb: internal error'
	      end if

	      sl = alength(ix-1,iy,2)
	      sr = alength(ix,iy,2)
	      sb = alength(ix,iy-1,1)
	      st = alength(ix,iy,1)

	      ul = fflux(ix-1,iy,2)
	      ur = fflux(ix,iy,2)
	      vb = fflux(ix,iy-1,1)
	      vt = fflux(ix,iy,1)

	      if( sl .gt. 0. ) then
		ul = ul / sl
	      else
		ul = 0.
	      end if
	      if( sr .gt. 0. ) then
		ur = ur / sr
	      else
		ur = 0.
	      end if
	      if( sb .gt. 0. ) then
		vb = vb / sb
	      else
		vb = 0.
	      end if
	      if( st .gt. 0. ) then
		vt = vt / st
	      else
		vt = 0.
	      end if

	      uaux(ix,iy) = 0.5 * ( ul + ur ) / h
	      vaux(ix,iy) = 0.5 * ( vb + vt ) / h

	    else
	      uaux(ix,iy) = 0.
	      vaux(ix,iy) = 0.
	    end if
	  end do
	end do

	write(6,*) 'writing data for alexb... ',it
	write(iunit,*) it
	call wrmx(iunit,nxdim,nydim,nx,ny,zlev)
	call wrmx(iunit,nxdim,nydim,nx,ny,uaux)
	call wrmx(iunit,nxdim,nydim,nx,ny,vaux)

	end

c**************************************************************

	subroutine wrmx(iunit,nxdim,nydim,nx,ny,matrix)

c writes matrix

	implicit none

	integer iunit
	integer nxdim,nydim
	integer nx,ny
	real matrix(nxdim,nydim)

	integer ix,iy

	write(iunit,*) ((matrix(ix,iy),ix=1,nx),iy=1,ny)

	end

c**************************************************************

	subroutine side2dxdy(nx,ny)

c computes dx/dy from side lenghts and areas

	implicit none

	integer nx,ny

	include 'regflux.h'

	real xfact,yfact
	parameter(xfact = 100./dxreg)	!transform sidelength (250m) in percent
	parameter(yfact = 100./dyreg)	!transform sidelength (250m) in percent

	integer i,j
	real area
	real sr,sl,su,sd
	real dx,dy,dxy

	open(44,file='dxdy.dat',status='unknown',form='formatted')

	do j=1,ny
	  do i=1,nx
	    area = areav(i,j)
	    if( area .gt. 0. ) then

	      sr = alength(i,j,2)
	      sl = alength(i-1,j,2)
	      su = alength(i,j,1)
	      sd = alength(i,j-1,1)

	      dx = su + sd
	      if( su * sd .gt. 0. ) dx = 0.5 * dx	!both sides present
	      dy = sr + sl
	      if( sr * sl .gt. 0. ) dy = 0.5 * dy	!both sides present

	      if( dx .eq. 0. .and. dy .eq. 0. ) then
		stop 'error stop side2dxdy: internal error (1)'
	      else if( dx .eq. 0. .or. dy .eq. 0. ) then
		dx = sqrt(area)
		dy = dx
	      else
		dxy = dx / dy
		if( dxy .gt. 1. ) then
		  dx = sqrt( area * dxy )
		  if( dx .gt. dxreg ) dx = dxreg
		  dy = area / dx
		else
		  dy = sqrt( area / dxy )
		  if( dy .gt. dyreg ) dy = dyreg
		  dx = area / dy
		end if
	      end if

	      dxside(i,j) = dx
	      dyside(i,j) = dy

	      write(44,'(2i5,f7.2)') i,j,100.*area/(dxreg*dyreg)
	      write(44,'(7f7.2)') xfact*su,xfact*sd,yfact*sl,yfact*sr,
     +			xfact*dx,yfact*dy,100.*dxy*dy/dx
	    else

	      dxside(i,j) = 0.
	      dyside(i,j) = 0.

	    end if
	  end do
	end do

	close(44)

	end

c**************************************************************

