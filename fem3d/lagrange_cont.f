c
c $Id: lagrange_cont.f,v 1.5 2009-09-14 08:20:57 georg Exp $
c
c simulates continuous release over open boundaries
c
c revision log :
c
c 12.12.2007    ggu	written from scratch
c 12.06.2008    ggu	initialize also z
c 28.08.2009    ggu	new call to find_elems_to_segment (before line_elems)
c 16.12.2011    ggu	new routine lagr_continuous_release_ppv()
c 23.01.2012    ggu	new routine for release in point, connectivity
c
c*******************************************************************

	subroutine lagr_continuous_release_shell

c manages continuous release of particles

	implicit none

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real pps,ppv
	integer itmin,itmax

	pps = 0.5
	pps = 0.		!parts per second (per boundary)
	ppv = 1.e-5		
	ppv = 0.		!parts per volume
	itmin = 86400
	itmax = 3*86400

	!if( it .ge. itmin .and. it .le. itmax ) then
	  call lagr_continuous_release_pps(pps)
	  call lagr_continuous_release_ppv(ppv)
	!end if

	call lagr_connect_continuous_points

	end

c*******************************************************************

	subroutine lagr_continuous_release_ppv(ppv)

c continuous release - number of particles depends on volume flux

	implicit none

        include 'param.h'
        include 'lagrange.h'

	real ppv	!particles per second

	integer k1,k2
	integer ibc,nk,i,ibtyp,np
	real dt
	real q

	logical bdebug
	integer nbnds,nkbnds,kbnds,nbc
	integer nptot,iptot
	real rp,pdens
	real dist_node

	real ggrand
	real get_bflux_ppv

	if( ppv .le. 0. ) return

	call get_timestep(dt)

	bdebug = .false.

	nbc = nbnds()

	do ibc=1,nbc

	  nk = nkbnds(ibc)
	  call get_bnd_ipar(ibc,'ibtyp',ibtyp)

	  if( ibtyp .eq. 1 ) then	!only for level boundaries

	    iptot = 0
	    rp = ggrand(77)		! vary starting point of particles
	    do i=2,nk
	      k1 = kbnds(ibc,i-1)
	      k2 = kbnds(ibc,i)
	      !if( k1 .eq. 6935 ) bdebug = .true.
	      q = get_bflux_ppv(k1,k2)
	      q = max(q,0.)
	      rp = rp + q*ppv*dt
	      np = rp
	      !if( bdebug ) write(6,*) k1,k2,q,rp,np
	      if( np .gt. 0 ) then
		rp = rp - np
		iptot = iptot + np
	        call create_parts(np,k1,k2)
	      end if
	    end do

	    if( iptot .ne. 0 ) then
	      write(lunit,*) 'number of particles released: ',ibc,iptot
	    end if

	  end if

	end do

	end

c*******************************************************************

	subroutine lagr_continuous_release_pps(pps)

c continuous release - number of particles is independent of boundary length

	implicit none

        include 'param.h'
        include 'lagrange.h'

	real pps	!particles per second

	integer k1,k2
	integer ibc,nk,i,ibtyp,np
	real totdist,dxy,part,dt

	integer nbnds,nkbnds,kbnds,nbc
	integer nptot,iptot
	real rp,pdens
	real dist_node

	real ggrand

	if( pps .le. 0. ) return

	call get_timestep(dt)

	nptot = pps*dt + ggrand(77)	! transform real into statistical np
	part = nptot
	if( nptot .le. 0. ) return

	nbc = nbnds()

	do ibc=1,nbc

	  nk = nkbnds(ibc)
	  call get_bnd_ipar(ibc,'ibtyp',ibtyp)

	  if( ibtyp .eq. 1 ) then	!only for level boundaries

	    totdist = 0.
	    do i=2,nk
	      k1 = kbnds(ibc,i-1)
	      k2 = kbnds(ibc,i)
	      dxy = dist_node(k1,k2)
	      totdist = totdist + dxy
	    end do

	    iptot = 0
	    pdens = part/totdist	! particles / m
	    rp = ggrand(77)		! vary starting point of particles
	    do i=2,nk
	      k1 = kbnds(ibc,i-1)
	      k2 = kbnds(ibc,i)
	      dxy = dist_node(k1,k2)
	      !np = nint(part*dxy/totdist)
	      rp = rp + dxy * pdens
	      np = rp
	      if( np .gt. 0 ) then
		rp = rp - np
		iptot = iptot + np
	        call create_parts(np,k1,k2)
	      end if
	    end do

	    if( iptot .ne. nptot ) then
	      write(lunit,*) 'number of particles differing: ',iptot,nptot
	    end if
	    if( iptot .ne. 0 ) then
	      write(lunit,*) 'number of particles released: ',iptot,nptot
	    end if

	  end if

	end do

	end

c*******************************************************************

	subroutine create_parts(np,k1,k2)

	implicit none

	integer np,k1,k2

	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv

	integer i,ie1,ie2
	real x1,y1,x2,y2,dx,dy
	real rl,rt,x,y

	real ggrand

	if( np .le. 0 ) return

	x1 = xgv(k1)
	y1 = ygv(k1)
	x2 = xgv(k2)
	y2 = ygv(k2)

	dx = x2 - x1
	dy = y2 - y1

	call find_elems_to_segment(k1,k2,ie1,ie2)
	if( ie1 .eq. 0 .or. ie2 .ne. 0 ) then
	  write(6,*) k1,k2,ie1,ie2
	  stop 'error stop create_parts: error in boundary'
	end if
	  
	!write(6,*) 'create_parts: ',np,k1,k2,ie1

	do i=1,np
	  rl = ggrand(77)
	  rt = ggrand(77)
	  x = x1 + rl*dx
	  y = y1 + rl*dy
	  call insert_particle(ie1,rt,x,y)
	  !write(6,*) i,rl,x,y,ie1
	end do

	end

c*******************************************************************

	function get_bflux_ppv(k1,k2)

	implicit none

	real get_bflux_ppv
	integer k1,k2

        include 'param.h'
        include 'lagrange.h'
        include 'links.h'

	integer ie1,ie2,ii

	integer inext

	call find_elems_to_segment(k1,k2,ie1,ie2)

	if( ie1 .eq. 0 .or. ie2 .ne. 0 ) then
	  write(6,*) 'k1,k2,ie1,ie2: ',k1,k2,ie1,ie2
	  stop 'error stop get_bflux_ppv: nodes not at boundary'
	end if

	ii = inext(k2,ie1)

	get_bflux_ppv = flux2d(ii,ie1)
	
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine release_on_point(ppts,ie,x,y,n)

	implicit none

	real ppts		!particles to be released per time step
	integer ie		!element where particle is released
	real x,y		!coordinates where to release
	integer n		!number of particles released (return)

	integer i
	real rt
	real ggrand

	n = ppts + ggrand(77)			!particles to release

	do i=1,n
	  rt = ggrand(77)			!vary time
	  call insert_particle(ie,rt,x,y)
	end do

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

