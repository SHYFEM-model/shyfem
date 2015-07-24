c
c $Id: baswork.f,v 1.20 2010-03-22 15:29:31 georg Exp $
c
c general framework to work on files needing basin
c
c revision log :
c
c 06.04.1999	ggu	completely restructured
c 04.06.1999	ggu	new statistics are computed
c 28.11.2005	ggu	new call to makehkv
c 31.05.2007	ggu	added area and volume frequency curve
c 24.08.2007	ggu	added new routine write_grd_from_bas
c 06.04.2009    ggu     read param.h
c 12.06.2009    ggu     areatr in double precision - new algorithm
c 01.03.2010    ggu     new routine basqual() to compute grid quality
c 22.03.2010    ggu     write external element number in basqual()
c 17.05.2011    ggu     changes in freqdep()
c 12.07.2011    ggu     better treatment of freqdep()
c 16.11.2011    ggu     basin.h introduced
c 23.01.2012    ggu     new from basinf
c
c****************************************************************

        program baswork

c works on basin data structure

	use mod_geom
	use mod_depth
	use evgeom
	use basin

	implicit none

	include 'param.h'




	real haux(nkndim)

	logical bnode,belem
	integer iapini

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        if( iapini(1,0,0,0) .le. 0 ) stop

	call bas_info

	call set_ev
	call set_geom

c-----------------------------------------------------------------
c specific info
c-----------------------------------------------------------------

	call basstat
        call makehev(hev)
        call makehkv(hkv,haux)

c-----------------------------------------------------------------
c specific applications
c-----------------------------------------------------------------

	call flume
	call distance

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

	subroutine flume

c creates boundary condition for circular flume

	use mod_depth
	use basin

	implicit none

	logical bnode

	include 'param.h'


	!call bnd_nodes(1563,551,-11)	!circular_flume_7680
	call bnd_nodes(1024,386,-11)	!Hafenmodell_7680

	end

c*****************************************************************

	subroutine bnd_nodes(kstart,kend,kd)

	use mod_geom
	use basin

	implicit none

	include 'param.h'

	integer kstart,kend,kd

	integer ndim
	parameter (ndim=100)

	integer kext(ndim)
	integer kint(ndim)
	real tang(2,ndim)


	integer ip,k,i,lmax
	integer ka,kb
	integer it0,it1
	real speed
	real x,y,xy
	logical berror

	it0 = 0
	it1 = 10000000
	lmax = 18
	speed = 0.1
	lmax = 14
	speed = 0.55

	ip = 0
	do k=kstart,kend,kd
	  ip = ip + 1
	  if( ip .gt. ndim ) stop 'error stop bnd_nodes: ndim'
	  kext(ip) = k
	  kint(ip) = k
	end do

	if( kext(ip) .ne. kend ) then
	  write(6,*) kstart,kend,kd,ip,kext(ip)
	  stop 'error stop bnd_nodes: end node'
	end if

	call n2int(ip,kint,berror)
	if( berror ) stop 'error stop bnd_nodes: internal nodes'

	do i=2,ip-1
	  k = kint(i)
	  if( kint(i-1) .ne. kantv(2,k) ) goto 99
	  if( kint(i+1) .ne. kantv(1,k) ) goto 99
	end do

	write(6,*) 'flume set-up: ',kstart,kend,kd,ip
	write(6,*) 'kbound = ',(kext(i),i=1,ip)

	do i=1,ip
	  k = kint(i)
	  kb = kantv(2,k)
	  ka = kantv(1,k)
	  x = xgv(ka) - xgv(kb)
	  y = ygv(ka) - ygv(kb)
	  xy = sqrt(x*x+y*y)
	  x = speed * x / xy
	  y = speed * y / xy
	  tang(1,i) = x
	  tang(2,i) = y
	end do

	call write_bound(it0,lmax,ip,tang)
	call write_bound(it1,lmax,ip,tang)

	return
   99	continue
	write(6,*) i,k,kint(i-1),kint(i+1),kantv(1,k),kantv(2,k)
	stop 'error stop bnd_nodes: node structure'
	end

c*****************************************************************

	subroutine write_bound(it0,lmax,ip,tang)

	implicit none

	integer it0,lmax,ip
	real tang(2,ip)
	real speed

	integer i,j
	real u,v

	write(77,*) it0,lmax,ip,2

	do i=1,ip
	  u = tang(1,i)
	  write(77,*) i,(u,j=1,lmax)
	end do

	do i=1,ip
	  v = tang(2,i)
	  write(77,*) i,(v,j=1,lmax)
	end do

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine basstat

c writes statistics on basin

	use basin

	implicit none

	include 'param.h'

	integer ie,ii,k
	integer imin,imax
	real area,amin,amax,atot
        real vtot
	real aptot,vptot
	real x(3),y(3)
	real xmin,xmax,ymin,ymax
	real dxmax,dymax
	real h
        real xx,yy
        real dist,distmin
        integer i,k1,k2

	real areatr

c-----------------------------------------------------------------
c area code
c-----------------------------------------------------------------

	imin = iarv(1)
	imax = imin

	do ie=1,nel
	  imin = min(imin,iarv(ie))
	  imax = max(imax,iarv(ie))
	end do

	write(6,*) 'Area code min/max:      ',imin,imax

c-----------------------------------------------------------------
c node numbers
c-----------------------------------------------------------------

	imin = ipv(1)
	imax = imin

	do k=1,nkn
	  imin = min(imin,ipv(k))
	  imax = max(imax,ipv(k))
	end do

	write(6,*) 'Node number min/max:    ',imin,imax

c-----------------------------------------------------------------
c element numbers
c-----------------------------------------------------------------

	imin = ipev(1)
	imax = imin

	do ie=1,nel
	  imin = min(imin,ipev(ie))
	  imax = max(imax,ipev(ie))
	end do

	write(6,*) 'Element number min/max: ',imin,imax

c-----------------------------------------------------------------
c area
c-----------------------------------------------------------------

	amin = areatr(1)
	amax = amin
	atot = 0.
        vtot = 0.
	aptot = 0.
        vptot = 0.

	do ie=1,nel
	  area = areatr(ie)
	  atot = atot + area
	  amin = min(amin,area)
	  amax = max(amax,area)
          h = 0.
          do ii=1,3
            h = h + hm3v(ii,ie)
          end do
          vtot = vtot + area * h / 3.
	  if( h .gt. 0. ) then		!only positive depths
	    aptot = aptot + area
            vptot = vptot + area * h / 3.
	  end if
	end do

	write(6,*) 'Area min/max:           ',amin,amax
	write(6,*) 'Total area:             ',atot,aptot
	write(6,*) 'Total volume:           ',vtot,vptot

c-----------------------------------------------------------------
c coordinates
c-----------------------------------------------------------------

	xmin = xgv(1)
	ymin = ygv(1)
	xmax = xgv(1)
	ymax = ygv(1)

	do k=1,nkn
	  xmin = min(xmin,xgv(k))
	  ymin = min(ymin,ygv(k))
	  xmax = max(xmax,xgv(k))
	  ymax = max(ymax,ygv(k))
	end do

	write(6,*) 'X-Coordinates min/max:  ',xmin,xmax
	write(6,*) 'Y-Coordinates min/max:  ',ymin,ymax

c-----------------------------------------------------------------
c size of elements
c-----------------------------------------------------------------

	dxmax = 0.
	dymax = 0.

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    x(ii) = xgv(k)
	    y(ii) = ygv(k)
	  end do
	  xmin = min(x(1),x(2),x(3))
	  xmax = max(x(1),x(2),x(3))
	  ymin = min(y(1),y(2),y(3))
	  ymax = max(y(1),y(2),y(3))
	  dxmax = max(dxmax,xmax-xmin)
	  dymax = max(dymax,ymax-ymin)
	end do

	write(6,*) 'Element dxmax/dymax:    ',dxmax,dymax

c-----------------------------------------------------------------
c depth
c-----------------------------------------------------------------

	amin = 999999.
	amax = -amin

	do ie=1,nel
	  h = 0
	  do ii=1,3
	    h = h + hm3v(ii,ie)
	  end do
	  h = h / 3.
	  amin = min(amin,h)
	  amax = max(amax,h)
	end do

	write(6,*) 'Depth min/max:          ',amin,amax
	write(6,*) 'Depth average:          ',vtot/atot,vptot/aptot

c-----------------------------------------------------------------
c minimum distance of nodes
c-----------------------------------------------------------------

        distmin = (xmax-xmin)**2 + (ymax-ymin)**2
        distmin = 2*distmin
        k1 = 0
        k2 = 0

        do k=1,nkn
          xx = xgv(k)
          yy = ygv(k)
          do i=k+1,nkn
            dist = (xx-xgv(i))**2 + (yy-ygv(i))**2
            if( dist .lt. distmin ) then
              k1 = k
              k2 = i
              distmin = dist
            end if
          end do
        end do

        distmin = sqrt(distmin)
	write(6,*) 'min node distance:      ',distmin,ipv(k1),ipv(k2)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	write(6,*)

	end

c*******************************************************************

        function areatr(ie)

c determination of area of element
c
c ie            number of element (internal)
c areatr        element area (return value)

	use basin

	real areatr
	integer ie

	include 'param.h'

	integer ii,i1,i2,k1,k2
	double precision f,x(3),y(3)

        do ii=1,3
          k=nen3v(ii,ie)
	  x(ii) = xgv(k)
	  y(ii) = ygv(k)
        end do

	f = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))

        areatr = f / 2.D0

        end

c*******************************************************************

	subroutine distance

	use basin

	implicit none

	include 'param.h'

	real d(nkndim)
	real dd(nkndim)

	integer k,ie,ii,it,is
	integer ngood
	real high

	high = 1.e+30

	do k=1,nkn
	  d(k) = -1.
	end do

	d(100) = 0.
	d(400) = 0.
	d(1000) = 0.

	ngood = 0
	do k=1,nkn
	  if( d(k) .ge. 0. ) ngood = ngood + 1
	end do

	do while( ngood .lt. nkn )

	  do k=1,nkn
	    if( d(k) .ge. 0. ) ngood = ngood + 1
	  end do

	  do ie=1,nel

	    do k=1,nkn
	      dd(k) = high
	    end do

	    it = 0
	    is = 0
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( d(k) .ge. 0 ) then
	        it = it + 1
	        is = is + ii
	      end if
	    end do
	    
	    if( it .eq. 1 ) then
	    else if( it .eq. 2 ) then
	    end if

	  end do

	end do

	end

c*******************************************************************


