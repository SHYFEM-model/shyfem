c
c $Id: bascomp.f,v 1.3 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 01.12.2005    ggu     created from basinf
c 06.04.2009    ggu     read param.h
c
c****************************************************************

        program bascomp

c compares depth of two basins

	use evgeom
	use basin

	implicit none

	include 'param.h'

	real hm3v1(3,neldim)

        character*20 file,dfile
	integer node,nit
	integer mode
	integer ie,ii
	logical bwrite
	integer iapini

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        if( iapini(1,nkndim,neldim,0) .le. 0 ) stop

c-----------------------------------------------------------------
c copy depth
c-----------------------------------------------------------------

	do ie=1,nel
	  do ii=1,3
	    hm3v1(ii,ie) = hm3v(ii,ie)
	  end do
	end do

c-----------------------------------------------------------------
c read in second basin
c-----------------------------------------------------------------

        if( iapini(1,nkndim,neldim,0) .le. 0 ) stop

	call bas_info

c-----------------------------------------------------------------
c test
c-----------------------------------------------------------------

	call set_ev

c-----------------------------------------------------------------
c difference
c-----------------------------------------------------------------

	do ie=1,nel
	  do ii=1,3
	    hm3v(ii,ie) = hm3v(ii,ie) - hm3v1(ii,ie)
	  end do
	end do

c-----------------------------------------------------------------
c write
c-----------------------------------------------------------------

	open(1,file='new.grd',status='unknown',form='formatted')
	call wrgrd(1)
	close(1)
        write(6,*) 'file has been written to new.grd'

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

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

	real aj
	integer ii,i1,i2,k1,k2

        aj=0.
        do ii=1,3
          i1=mod(ii,3)+1
          i2=mod(i1,3)+1
          k1=nen3v(i1,ie)
          k2=nen3v(i2,ie)
          aj=aj+xgv(k1)*ygv(k2)-xgv(k2)*ygv(k1)
        end do

        areatr = aj / 2.

        end

c*******************************************************************


c*******************************************************************

	subroutine wrgrd(iunit)

c writes grd file from bas

	use basin

	implicit none

	integer iunit

	include 'param.h'

	integer k,ie,ii

	do k=1,nkn
	  write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k)
	end do

	write(iunit,*)

	do ie=1,nel
	  write(iunit,1100) 2,ipev(ie),iarv(ie)
     +		,3,(ipv(nen3v(ii,ie)),ii=1,3),hm3v(1,ie)
	end do

	return
 1000	format(i1,2i10,2e16.8)
 1100	format(i1,2i10,i4,3i10,e16.8)
	end

c*******************************************************************

c*******************************************************************

	subroutine neibors(node,nit)

c gets neibors of node node (nit iterations)

	use basin

	implicit none

	integer node
	integer nit

	include 'param.h'

	integer nodes(nkndim)
	integer icol(nkndim)
	integer iaux(nkndim)
	integer nintern,icolor
	integer k,n,ie,ii,i
	logical bcol

	integer ipext,ipint

	nintern = ipint(node)
	if( nintern .le. 0 ) then
	  write(6,*) 'node,nintern: ',node,nintern
	  stop 'error stop neibors: no such node'
	end if

	do k=1,nkn
	  icol(k) = -1
	  iaux(k) = 0
	end do

	icolor = 0
	icol(nintern) = icolor

	do n=1,nit

	  icolor = icolor + 1

	  do ie=1,nel
	    bcol = .false.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( icol(k) .ge. 0 ) bcol = .true.
	    end do

	    if( bcol ) then
	      do ii=1,3
	        k = nen3v(ii,ie)
	        if( icol(k) .lt. 0 ) iaux(k) = icolor
	      end do
	    end if
	  end do

	  do k=1,nkn
	    if( iaux(k) .gt. 0 ) then
	      icol(k) = iaux(k)
	      iaux(k) = 0
	    end if
	  end do

	end do

	n = 0
	do k=1,nkn
	  if( icol(k) .ge. 0 ) then
	    !write(6,'(3i10)') k,ipext(k),icol(k)
	    n = n + 1
	    nodes(n) = ipext(k)
	  end if
	end do

	write(6,*)
	write(6,*) (nodes(i),i=1,n)
	write(6,*)

	end

c*******************************************************************

	subroutine wnodes(file)

	use basin

	implicit none

	character*(*) file

	include 'param.h'

	integer node,it,k
	integer ipint

	open(1,file=file,status='old',form='formatted')

    1	continue
	  read(1,*,end=2) node,it
	  k = ipint(node)
	  if( k .le. 0 ) goto 99
	  write(6,1000) 1,node,it,xgv(k),ygv(k)
	  goto 1
    2	continue

	close(1)

	return
 1000	format(i1,2i10,2f14.4)
   99	continue
	write(6,*) 'node: ',node,' not available'
	stop 'error stop wnodes: no such node'
	end

c*******************************************************************

	subroutine wnodeef(file)

	use basin

	implicit none

	character*(*) file

	include 'param.h'

	integer ndim
	parameter (ndim=1000)

	integer nodes(nkndim)

	integer node,it,k
	integer ipint
	integer line(ndim)
	integer inode,iline,i

	do k=1,nkn
	  nodes(k) = 0
	end do

	open(1,file=file,status='old',form='formatted')
	open(2,file='extra.grd',status='unknown',form='formatted')

	inode = 0	!number of nodes in list
	iline = 0	!number of line
    1	continue
	  read(1,*,end=2) node,it
	  k = ipint(node)
	  if( k .le. 0 .and. node .gt. 0 ) goto 99
	  if( it .ge. 0 ) then
	    if( nodes(node) .eq. 0 ) then
	      write(2,1000) 1,node,it,xgv(k),ygv(k)
	      nodes(node) = 1
	    end if
	  else
	    if( node .gt. 0 ) then
	      if( nodes(node) .eq. 0 ) then
	        write(2,1000) 1,node,0,xgv(k),ygv(k)
	        nodes(node) = 1
	      end if
	      inode = inode + 1
	      if( inode .gt. ndim ) stop 'error stop wnodeef: ndim'
	      line(inode) = node
	    else
	      if( inode .gt. 0 ) then
		iline = iline + 1
	        write(2,2000) 3,iline,0,inode
		do i=1,inode
		  write(2,3000) line(i)
		end do
		inode = 0
	      end if
	    end if
	  end if
	  goto 1
    2	continue

	close(1)
	close(2)

	return
 1000	format(i1,2i10,2f14.4)
 2000	format(i1,3i10)
 3000	format(i12)
   99	continue
	write(6,*) 'node: ',node,' not available'
	stop 'error stop wnodes: no such node'
	end

c*******************************************************************

	subroutine wrdepth(ntot)

c writes depth values distributed in elements

	use basin

	implicit none

	include 'param.h'

	integer ntot	!number of subdivisions in element

	integer ie,ii,k,i,j,itot
	integer node,itype
	real x(3),y(3)
	real x1,x2,y1,y2,xp,yp
	real r,depth

c-----------------------------------------------------------------
c loop over elements
c-----------------------------------------------------------------

	node = 0
	itype = 0

	open(2,file='depth.grd',status='unknown',form='formatted')
	open(1,file='depth.dat',status='unknown',form='formatted')

	do ie=1,nel

	  depth = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    x(ii) = xgv(k)
	    y(ii) = ygv(k)
	    depth = depth + hm3v(ii,ie)
	  end do
	  depth = depth / 3.

	  do i=1,ntot
	    itot = i*2 - 1
	    r = itot/(2.*ntot)
	    x1 = x(3) + r * ( x(1) - x(3) )
	    x2 = x(3) + r * ( x(2) - x(3) )
	    y1 = y(3) + r * ( y(1) - y(3) )
	    y2 = y(3) + r * ( y(2) - y(3) )
	    do j=1,itot
	      r = (2*j-1)/(2.*itot)
	      xp = x1 + r * ( x2 - x1 )
	      yp = y1 + r * ( y2 - y1 )
	      node = node + 1
	      write(2,1000) 1,node,itype,xp,yp,depth
	      write(1,2000) xp,yp,depth
	    end do
	  end do

	end do

	close(2)
	close(1)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
 1000	format(i1,2i10,3f14.4)
 2000	format(3f14.4)
	end

c*******************************************************************

	subroutine readbat(file,np,xp,yp,dp)

	implicit none

	character*(*) file
	integer np
	real xp(1)
	real yp(1)
	real dp(1)

	integer ndim,n
	integer iaux,inum,ityp
	real x,y,d

	ndim = np

	n = 0
	write(6,*) 'reading file : ',file
	open(1,err=99,file=file,status='old',form='formatted')

    1	continue
	  read(1,*,end=2) iaux,inum,ityp,x,y,d
	  n = n + 1
	  if( n .gt. ndim ) stop 'error stop: ndim'
	  xp(n) = x
	  yp(n) = y
	  dp(n) = d
	  goto 1
    2	continue

	close(1)
	np = n
	write(6,*) np,' data points read'

        return
   99   continue
        write(6,*) 'Cannot read file: ',file
        stop 'error stop readbat: no such file'
	end

c*******************************************************************

	subroutine triab(x,y,area,x0,y0)

c computes area and center point of triangle

	implicit none

	real x(3)
	real y(3)
	real area,x0,y0

	area = 0.5 * ( (x(2)-x(1))*(y(3)-y(1)) 
     +			- (x(3)-x(1))*(y(2)-y(1)) )

	x0 = (x(1)+x(2)+x(3))/3.
	y0 = (y(1)+y(2)+y(3))/3.

	end

c*******************************************************************

	subroutine interpole(file)

c interpolates depth values

	use basin

	implicit none

	character*(*) file

	include 'param.h'

	integer ie,ii,k
	integer netot
	real x,y,d
	real hmin
	real depth
	integer iaux,inum,ityp
	integer n,np,i
	integer ihmed
        integer nintp
	real xt(3), yt(3)
	real xmin,ymin,xmax,ymax
	real weight,r2,w
        real r2max,sigma2
	real hmed
	real fact,dx,dy
	real area,x0,y0
	real hev(neldim)
	integer ihev(neldim)
	logical bmin
	logical ok(neldim)
	logical inconvex,inquad

	integer ndim
	parameter(ndim=200000)
	real xp(ndim), yp(ndim)
	real dp(ndim)

	bmin = .true.
	hmin = 1.
	hmin = 0.5

c-----------------------------------------------------------------
c read points
c-----------------------------------------------------------------

	np = ndim
	call readbat(file,np,xp,yp,dp)

c-----------------------------------------------------------------
c initialize
c-----------------------------------------------------------------

	netot = 0

	do ie=1,nel
	  hev(ie) = 0.
	  ihev(ie) = 0
	  ok(ie) = .false.
	  depth = 0.
	  do ii=1,3
	    depth = depth + hm3v(ii,ie)
	  end do
	  depth = depth / 3.
	  if( depth .gt. -990 ) then
	    hev(ie) = depth
	    ok(ie) = .true.
	    netot = netot + 1
	  end if
	end do

	write(6,*) 'Elements without depth (start): ',nel-netot,nel

c-----------------------------------------------------------------
c initial interpolation -> point in element
c-----------------------------------------------------------------

	fact = 2.

	do while( netot .lt. nel )

	netot = 0

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    xt(ii) = xgv(k)
	    yt(ii) = ygv(k)
	  end do

	  call triab(xt,yt,area,x0,y0)

          r2max = fact*area     !maximum radius to look for points
          sigma2 = fact*area    !standard deviation also grows

	  if( ok(ie) ) then
	    netot = netot + 1
	  else
	    depth = 0.
	    weight = 0.
            nintp = 0
	    do n=1,np
	      x = xp(n)
	      y = yp(n)
	      d = dp(n)
	      r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0)
	      if( r2 .le. r2max ) then
	        !w = exp(-r2/area)
	        w = exp(-r2/sigma2)
	        depth = depth + d * w
	        weight = weight + w
                nintp = nintp + 1
	      end if
	    end do
	    if( nintp .gt. 0 ) then
	      if( weight .le. 0. ) then
                write(6,*) nintp,weight,r2max,ie
                stop 'error stop interpole: zero weight from points'
              end if
	      depth = depth / weight
	      hev(ie) = depth
	      ok(ie) = .true.
	      netot = netot + 1
	    end if
	  end if
	end do

	write(6,*) 'Elements without depth : ',fact,nel-netot,nel

	fact = fact * 2.

	end do

c-----------------------------------------------------------------
c end up
c-----------------------------------------------------------------

	netot = 0

	do ie=1,nel
	  if( ok(ie) ) then
	    netot = netot + 1
	    hmed = hev(ie)
	    !if( hmed .lt. hmin ) hmed = hmin
	  else
	    hmed = hmin
	  end if
	  if( bmin .and. hmed .lt. hmin ) hmed = hmin
	  do ii=1,3
	    hm3v(ii,ie) = hmed
	  end do
	end do

	write(6,*) 'Elements without depth (end): ',nel-netot,nel

	if( bmin ) then
	  write(6,*) 'minimum depth used with hmin = ',hmin
	end if

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************


c*******************************************************************
