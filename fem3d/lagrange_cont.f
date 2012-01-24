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

	call lagr_continuous_points

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
	      if( k1 .eq. 6935 ) bdebug = .true.
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

	subroutine lagr_continuous_points

c continuous release from points

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer np
	integer iep(nconnect_dim)
	real xp(nconnect_dim),yp(nconnect_dim)
	save np,iep,xp,yp

	integer ip,n
	real pps,ppts,dt

	integer icall
	save icall
	data icall / 0 /

	if( icall .lt. 0 ) return

	pps = 1./200.					!particles per second

	if( icall .eq. 0 ) then				!initialize
	  if( pps .le. 0. ) then
	    icall = -1
	    return
	  end if
	  !call lagr_connect_get_coords_1(nconnect_dim,np,xp,yp)
	  call lagr_connect_get_coords(nconnect_dim,np,xp,yp)
	  call lagr_connect_find_elems(np,xp,yp,iep)

	  call lagr_connect_init(np,iep)
	end if

	icall = icall + 1

	call get_timestep(dt)
	ppts = dt * pps				!paricles per time step

	do ip=1,np
	  call release_on_point(ppts,iep(ip),xp(ip),yp(ip),n)
	  call lagr_connect_released(ip,n)	!count released particles
	end do

	if( it .eq. itend ) then
	  call lagr_connect_write(np)
	end if

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lagr_connect_init(np,iep)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer np
	integer iep(np)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer ie,ip,jp

	real areaele

	np_station = np

	do ie=1,nel
	  i_connect_elems(ie) = 0
	  i_connect_total(ie) = 0
	  t_connect_total(ie) = 0.
	end do

	do ip=1,np
	  ie = iep(ip)
	  i_connect_elems(ie) = ip
	  a_connect_area(ip) = areaele(ie)
	end do

	do ip=1,np
	  do jp=1,np
	    i_connect(ip,jp) = 0
	    t_connect(ip,jp) = 0.
	    if_connect(ip,jp) = 0
	    tf_connect(ip,jp) = 0.
	  end do
	  i_connect_released(ip) = 0
	end do

	end

c*******************************************************************

	subroutine lagr_connect_released(ip,n)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer ip,n

	i_connect_released(ip) = i_connect_released(ip)  + n

	end

c*******************************************************************

	subroutine lagr_connect_count(ibdy,ie,time,ic)

	implicit none

	include 'param.h'
	include 'lagrange.h'
	include 'lagrange_connect.h'

	integer ibdy,ie
	real time
	integer ic

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer ie_from,ip_to,ip_from
	real tarrive

	if( ie .le. 0 ) return

	i_connect_total(ie) = i_connect_total(ie) + ic
	t_connect_total(ie) = t_connect_total(ie) + time

	ip_to = i_connect_elems(ie)
	if( ip_to .le. 0 ) return

	ie_from = est(ibdy)
	ip_from = i_connect_elems(ie_from)
	if( ip_from .le. 0 ) then
	  write(6,*) 'particle not coming from source: ',ibdy,ie_from
	  stop 'error stop lagr_connect_count: no source'
	end if

	i_connect(ip_from,ip_to) = i_connect(ip_from,ip_to) + ic
	t_connect(ip_from,ip_to) = t_connect(ip_from,ip_to) + time

	if( ic .le. 0 ) return

	if( ip_to .gt. 63 ) then
	  write(6,*) 'ip_to = ',ip_to
	  stop 'error stop lagr_connect_count: internal error bitmap'
	end if

	if( btest(lgr_bitmap(ibdy),ip_to) ) return
	lgr_bitmap(ibdy) = ibset(lgr_bitmap(ibdy),ip_to)

	tarrive = it - tin(ibdy)
	if_connect(ip_from,ip_to) = if_connect(ip_from,ip_to) + ic
	tf_connect(ip_from,ip_to) = tf_connect(ip_from,ip_to) + tarrive

	end

c*******************************************************************

	subroutine lagr_connect_write(np)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer np

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
	real hev(1)
	common /hev/hev

	integer nlv,nvar,nvers,ierr,iunit,ivar
	integer ilhv(1)
	real hlv(1)
	character*80 file,title

	integer ip,jp,ie
	integer ic,icf,itot
	real tc,tcf,pab
	real area_i,area_j

	integer ifileo

	write(127,*) np

	do ip=1,np
	  write(127,*) ip,i_connect_released(ip),a_connect_area(ip)
	end do

	write(127,*) 'total matrix'

	do ip=1,np				!from
	  do jp=1,np				!to
	    call lagr_connect_write_entry(np,ip,jp)
	  end do
	end do

	write(127,*) 'lower matrix (column-wise)'

	do jp=1,np				!to
	  do ip=jp,np				!from
	    call lagr_connect_write_entry(np,ip,jp)
	  end do
	end do

	write(127,*) 'upper matrix (row-wise)'

	do ip=1,np				!from
	  do jp=ip,np				!to
	    call lagr_connect_write_entry(np,ip,jp)
	  end do
	end do

	file = 'connectivity.eos'
	nvers = 3
	nvar = 1
	nlv = 1
	ivar = 300
	title = 'connectivity'

	iunit = ifileo(0,file,'unformatted','new')
	call wfeos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
	call wseos(iunit,ilhv,hlv,hev,ierr)
	call wreos(iunit,it,ivar+1,1,ilhv,t_connect_total,ierr)
	!call wreos(iunit,it,ivar+2,nlvdim,ilhv,c,ierr)
	close(iunit)

	do ie=1,nel
	  write(128,*) ie,i_connect_total(ie),t_connect_total(ie)
	end do

	end

c*******************************************************************

	subroutine lagr_connect_write_entry(np,ip,jp)

	implicit none

	include 'param.h'
	include 'lagrange_connect.h'

	integer np,ip,jp

	integer ic,icf,itot
	real tc,tcf,pab
	real area_i,area_j

c ip is from
c jp is to

	area_i = a_connect_area(ip)
	itot = i_connect_released(ip)
	area_j = a_connect_area(jp)
	ic = i_connect(ip,jp)
	tc = t_connect(ip,jp)
	icf = if_connect(ip,jp)
	tcf = tf_connect(ip,jp)
	pab = (i_connect(ip,jp)/area_j) / (itot/area_i)
	if( icf .gt. 0 ) tcf = tcf / icf
	write(127,1000) ip,jp,ic,tc,icf,tcf,pab
 1000	format(2i3,i8,e15.7,i8,2e15.7)

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lagr_connect_find_elems(np,xp,yp,iep)

	implicit none

	integer np
	real xp(np),yp(np)
	integer iep(np)

	integer ip,ie

	do ip=1,np
	  call find_element(xp(ip),yp(ip),ie)
	  if( ie .le. 0 ) then
	    write(6,*) 'no such element: ',ip,xp(ip),yp(ip)
	    stop 'error stop find_elems_to_points: no such element'
	  end if
	  iep(ip) = ie
	end do

	end

c*******************************************************************

	subroutine lagr_connect_get_coords_1(ndim,np,xp,yp)

	implicit none

	integer ndim,np
	real xp(ndim),yp(ndim)

	np = 1
	xp(np) = 25000.
	yp(np) = 80000.

	end

c*******************************************************************

        subroutine lagr_connect_get_coords(ndim,np,xp,yp)

        implicit none

        integer ndim
        integer np,i
        real xp(ndim),yp(ndim)
        real xx,yy

        np=0

        open(1,file='coord_menor.dat',status='old',err=99)

10      continue
        read(1,*,err=98,end=12)xx,yy
        np = np +1
        if (np.gt.ndim)then
          write(6,*)'lagr_get_coords np ERROR',np,ndim
          stop
        endif
        xp(np)=xx
        yp(np)=yy
        goto 10

12      continue
        close(1)
         write(6,*)'connectivity stations',np
        do i=1,np
         write(6,*)xp(i),yp(i)
        enddo

        RETURN

99      continue
        write(6,*)'lagr_get_coords OPEN FILE ERROR'
        STOP

98      continue
        write(6,*)'lagr_get_coords READING FILE ERROR'
        STOP
        end

c******************************************************************

	subroutine lagr_connect_get_station(ie,ip_station,r_station)

	implicit none

	include 'param.h'
	include 'lagrange.h'
	include 'lagrange_connect.h'

	integer ie,ip_station
	real r_station

	ip_station = i_connect_elems(ie)
	r_station = ip_station / float(np_station)

	end

c******************************************************************
