!
! $Id: lagrange_cont.f,v 1.5 2009-09-14 08:20:57 georg Exp $
!
! simulates continuous release over open boundaries
!
! revision log :
!
! 12.12.2007    ggu	written from scratch
! 12.06.2008    ggu	initialize also z
! 28.08.2009    ggu	new call to find_elems_to_segment (before line_elems)
! 16.12.2011    ggu	new routine lagr_continuous_release_ppv()
! 23.01.2012    ggu	new routine for release in point, connectivity
! 23.10.2012    ggu	do not call connectivity here anymore
! 28.03.2014    ggu	new routine lagr_continuous_release_pps_ppv()
!
!*******************************************************************
!-------------------------------------------------------------------------
        module lagrange_cont
!-------------------------------------------------------------------------
        contains
!-------------------------------------------------------------------------

	subroutine lagr_continuous_release_shell

	implicit none

	!call lagr_continuous_release_pps
	!call lagr_continuous_release_ppv
	call lagr_continuous_release_pps_ppv

	end

!*******************************************************************

	subroutine lagr_continuous_release_ppv

! continuous release - number of particles depends on volume flux

	use lagrange_data
        use bnd_admin
        use random_gen
        use time_util

	implicit none

        include 'param.h'

	integer k1,k2
	integer ibc,nk,i,ibtyp,np
	double precision dt
	double precision q,ppv

	logical bdebug
	integer nbc
	integer nptot,iptot
	double precision rp

	call get_timestep(dt)

	bdebug = .false.

	nbc = nbnds()

	do ibc=1,nbc

	  nk = nkbnds(ibc)
	  call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	  call get_bnd_par(ibc,'lgrppv',ppv)

	  if( ibtyp .eq. 1 ) then	!only for level boundaries

	    iptot = 0
	    rp = ggrand(77)		! vary starting point of particles
	    do i=2,nk
	      k1 = kbnds(ibc,i-1)
	      k2 = kbnds(ibc,i)
	      q = get_bflux_ppv(k1,k2)
	      q = max(q,0.)
	      rp = rp + q*ppv*dt
	      np = rp
	      if( np .gt. 0 ) then
		rp = rp - np
		iptot = iptot + np
	        call create_parts(ibc,np,k1,k2)
	      end if
	    end do

	    if( iptot .ne. 0 ) then
	      write(lunit,*) 'number of particles released: ',ibc,iptot
	    end if

	  end if

	end do

	end

!*******************************************************************

	subroutine lagr_continuous_release_pps

! continuous release - number of particles is independent of boundary length

	use lagrange_data
	use lagrange_util
        use bnd_admin
        use random_gen
        use time_util

	implicit none

        include 'param.h'

	integer k1,k2
	integer ibc,nk,i,ibtyp,np
	double precision totdist,dxy,part,dt

	integer nbc
	integer iptot
	double precision rp,pps,q

	call get_timestep(dt)

	nbc = nbnds()

	do ibc=1,nbc

	  nk = nkbnds(ibc)
	  call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	  call get_bnd_par(ibc,'lgrpps',pps)

	  if( ibtyp .eq. 1 ) then	!only for level boundaries

	    call dist_total(ibc,totdist)

	    iptot = 0
	    rp = ggrand(77)		! vary starting point of particles
	    do i=2,nk
	      k1 = kbnds(ibc,i-1)
	      k2 = kbnds(ibc,i)
	      dxy = dist_node(k1,k2)
	      q = dxy/totdist
	      rp = rp + q*pps*dt
	      np = rp
	      if( np .gt. 0 ) then
		rp = rp - np
		iptot = iptot + np
	        call create_parts(ibc,np,k1,k2)
	      end if
	    end do

	    if( iptot .ne. 0 ) then
	      write(lunit,*) 'number of particles released: ',iptot
	    end if

	  end if

	end do

	end

!*******************************************************************

	subroutine lagr_continuous_release_pps_ppv

! continuous release - works both for pps and ppv
!
! replaces the routines above

	use lagrange_data
	use lagrange_util
        use random_gen
        use bnd_admin
        use time_util

	implicit none

        include 'param.h'

	integer k,k1,k2
	integer ibc,nk,i,ibtyp,np
	double precision totdist,dxy,part,dt

	logical bflux
	integer nbc
	integer iptot
	integer it
	double precision rp,pps,q

	call get_timestep(dt)
	call get_act_time(it)

	nbc = nbnds()

	do ibc=1,nbc

	  nk = nkbnds(ibc)
	  call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	  call get_bnd_par(ibc,'lgrpps',pps)

	  bflux = pps .lt. 0.	!if pps is negative it stands for parts per vol
	  pps = abs(pps)

	  if( pps .gt. 0. ) then
	   if( ibtyp .eq. 1 .or. ibtyp .eq. 2 ) then	!level or flux bnds

	    call dist_total(ibc,totdist)

	    iptot = 0
	    rp = ggrand(77)		! vary starting point of particles
	    do i=2,nk
	      k1 = kbnds(ibc,i-1)
	      k2 = kbnds(ibc,i)
	      if( bflux ) then
		q = get_bflux_ppv(k1,k2)
	      else
	        q = dist_node(k1,k2) / totdist
	      end if
	      q = max(q,0.)
	      rp = rp + q*pps*dt
	      np = rp
	      if( np .gt. 0 ) then
		rp = rp - np
		iptot = iptot + np
	        call create_parts(ibc,np,k1,k2)
	      end if
	    end do

	   else if( ibtyp .eq. 3 ) then		!release on nodes

	    iptot = 0
	    rp = ggrand(77)		! vary starting point of particles
	    do i=1,nk
	      k = kbnds(ibc,i)
	      if( bflux ) then
		call get_bnd_par(ibc,'zval',q)
	      else
	        q = 1
	      end if
	      rp = rp + q*pps*dt 
	      call release_on_node(ibc,rp,k,np)
	      rp = rp - np
	      iptot = iptot + np
	    end do

	   end if
	  end if

	  if( pps .ne. 0 ) then
	    write(lunit,*) 'particles released: ',bflux,ibc,pps,iptot
	  end if

	end do

	end

!*******************************************************************

	subroutine create_parts(ity,np,k1,k2)

	use basin
        use lnku
        use random_gen
        use lagrange_inout

	implicit none

	integer ity
	integer np,k1,k2

	include 'param.h'

	integer i,ie1,ie2
	double precision x1,y1,x2,y2,dx,dy
	double precision rl,rt,x,y

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
	  call insert_particle_3d(ie1,ity,rt,x,y)
	  !write(6,*) i,rl,x,y,ie1
	end do

	end

!*******************************************************************

	function get_bflux_ppv(k1,k2)

	use lagrange_data
	use geom
        use lnku

	implicit none

	double precision get_bflux_ppv
	integer k1,k2

        include 'param.h'

	integer ie1,ie2,ii

	call find_elems_to_segment(k1,k2,ie1,ie2)

	if( ie1 .eq. 0 .or. ie2 .ne. 0 ) then
	  write(6,*) 'k1,k2,ie1,ie2: ',k1,k2,ie1,ie2
	  stop 'error stop get_bflux_ppv: nodes not at boundary'
	end if

	ii = inext(k2,ie1)

	get_bflux_ppv = flux2d(ii,ie1)
	
	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine release_on_node(ity,ppts,k,n)

! release on node

	use basin

	implicit none

	integer ity		!type of particle
	double precision ppts		!particles to be released per time step
	integer k		!node where particle is released
	integer n

	include 'param.h'

	integer ie
	double precision x,y

	call find_elem_to_node(k,ie,x,y)
	!write(6,*) 'inserting... ',k,ie,x,y

	call release_on_point(ity,ppts,ie,x,y,n)

	end

!*******************************************************************

	subroutine find_elem_to_node(k,iee,x,y)

	use basin

! be sure particle is in an element

        use random_gen

	implicit none

	integer k
	integer iee
	double precision x,y

	integer ie,ii,kk,n,i
	double precision xs,ys,alpha
	integer ielist(ngr)

	n = 0
	alpha = 0.99

	do ie=1,nel
	  do ii=1,3
	    kk = nen3v(ii,ie)
	    if( kk .eq. k ) then
	      n = n + 1
	      if( n > ngr ) stop 'error stop find_elem_to_node: internal'
	      ielist(n) = ie
	    end if
	  end do
	end do

	if( n == 0 ) stop 'error stop find_elem_to_node: internal 2'

	i = 1 + n*ggrand(77)
	if( i > n ) i = n

	iee = ielist(i)
	ie = iee

	xs = 0.
	ys = 0.
	n = 0
	do ii=1,3
	  kk = nen3v(ii,ie)
	  if( k /= kk ) then
	    xs = xs + xgv(kk)
	    ys = ys + ygv(kk)
	  else
	    n = n + 1
	  end if
	end do

	if( n /= 1 ) stop 'error stop find_elem_to_node: internal 3'

	xs = xs / 2.
	ys = ys / 2.

	x = alpha*xgv(k) + (1.-alpha)*xs
	y = alpha*ygv(k) + (1.-alpha)*ys
	
	end

!*******************************************************************

	subroutine release_on_point(ity,ppts,ie,x,y,n)

! release from one point - works for 2D

        use random_gen
        use lagrange_inout

	implicit none

	integer ity		!type of particle
	double precision ppts		!particles to be released per time step
	integer ie		!element where particle is released
	double precision x,y		!coordinates where to release
	integer n		!number of particles released (return)

	integer i
	double precision rt

	n = ppts + ggrand(77)			!particles to release

	do i=1,n
	  rt = ggrand(77)			!vary time
	  call insert_particle_3d(ie,ity,rt,x,y)
	  !write(55,*) 'gguuyy particle: ',ie,rt,x,y
	end do

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

!-------------------------------------------------------------------------
        end module lagrange_cont
!-------------------------------------------------------------------------
