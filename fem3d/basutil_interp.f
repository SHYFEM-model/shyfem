
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c revision log :
c
c 21.02.2014	ggu	subroutines copied from basbathy
c 06.05.2015	ggu	set number of max iterations (imax)
c 10.10.2015	ggu	changed VERS_7_3_2
c 16.12.2015	ggu	depth ht is now passed in for square interpol
c 11.04.2016	ggu	meaning of ufact has changed, expo interp working
c 24.01.2018	ggu	changed VERS_7_5_41
c 18.12.2018	ggu	changed VERS_7_5_52
c 16.02.2019	ggu	changed VERS_7_5_60
c
c****************************************************************

	!program interpol
	!end

c*******************************************************************

	function dist2(x1,y1,x2,y2)

c computes squared distance (also for latlon coordinates)

	implicit none

	real dist2
	real x1,y1,x2,y2

	real pi,rad
	parameter ( pi = 3.14159 , rad = pi / 180. )

	integer latlon
	common /latlon/latlon

	real fact,y

	if( latlon .eq. 1 ) then

	  y = 0.5 * (y1+y2)
	  fact = cos( rad * y )
	  dist2 = (fact*(x2-x1))**2 + (y2-y1)**2

	else

	  dist2 = (x2-x1)**2 + (y2-y1)**2

	end if

	end
	
c*******************************************************************

	subroutine set_dist(isphe)

	implicit none

	integer isphe

	integer latlon
	common /latlon/latlon
	save /latlon/

	latlon = isphe

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine interpolq(np,xp,yp,dp,ht)

c interpolates depth values - works only on elements
c
c interpolation on squares

	use mod_depth
	use basin

	implicit none

	integer np
	real xp(np)
	real yp(np)
	real dp(np)
	real ht(max(nkn,nel))

	integer ie,ii,k
	integer netot,nnone
	integer imax
	real x,y,d
	real depth
	integer iaux,inum,ityp
	integer n,i
	integer ihmed
	real xt(3), yt(3)
	real xmin,ymin,xmax,ymax
	real hmed
	real fact,dx,dy
	real flag
	integer ihev(nel)
	logical bmin
	logical ok(nel)
	logical inconvex,inquad

c-----------------------------------------------------------------
c set maximum iteration imax
c
c imax very high makes all necessary iterations to find depth
c imax < 0 only looks in element
c imax = 0 looks also in bordering rectangle
c imax > 0 makes imax iterations increasing rectangle size
c-----------------------------------------------------------------

	imax = 0
	imax = 9999999
	!imax = -1		!do not compute over square

	flag = -999.

c-----------------------------------------------------------------
c initialize
c-----------------------------------------------------------------

	netot = 0

	do ie=1,nel
	  depth = hev(ie)
	  ihev(ie) = 0
	  if( depth .gt. -990 ) then
	    ok(ie) = .true.
	    netot = netot + 1
          else
	    ok(ie) = .false.
            hev(ie) = 0.
	  end if
	end do

	nnone = nel-netot
	write(6,*) 'Elements without depth (start): ',nnone,'/',nel

c-----------------------------------------------------------------
c initial interpolation -> point in element
c-----------------------------------------------------------------

	do n=1,np
	  x = xp(n)
	  y = yp(n)
	  d = dp(n)

	  if( mod(n,100) .eq. 0 ) then
	    write(6,*) n,(100.*n)/np,x,y,d
	  end if

	  do ie=1,nel
	    if( .not. ok(ie) ) then
	      do ii=1,3
	        k = nen3v(ii,ie)
	        xt(ii) = xgv(k)
	        yt(ii) = ygv(k)
	      end do
	      if( inconvex(3,xt,yt,x,y) ) then
	        hev(ie) = hev(ie) + d
	        ihev(ie) = ihev(ie) + 1
	      end if
	    end if
	  end do
	end do

	netot = 0

	do ie=1,nel
	  if( ihev(ie) .gt. 0 ) then
	    hev(ie) = hev(ie) / ihev(ie)
	    ok(ie) = .true.
	    netot = netot + 1
	  else if( ok(ie) ) then
	    netot = netot + 1
	  end if
	end do

	write(6,*) 'Elements without depth (start):  ',nnone,'/',nel
	write(6,*) 'Elements without depth (convex): ',nel-netot,'/',nel
	write(6,'(a)') 'interpol      iter   without      with     total'
	write(6,1000) 'convex: ',-1,nel-netot,netot,nel

c-----------------------------------------------------------------
c next interpolation -> point in square
c-----------------------------------------------------------------

	i = 0
	do while( netot .lt. nel .and. i .le. imax )
	 netot = 0
	 do ie=1,nel
	  if( .not. ok(ie) ) then
	    do ii=1,3
	      k = nen3v(ii,ie)
	      xt(ii) = xgv(k)
	      yt(ii) = ygv(k)
	    end do

	    xmin = min(xt(1),xt(2),xt(3))
	    ymin = min(yt(1),yt(2),yt(3))
	    xmax = max(xt(1),xt(2),xt(3))
	    ymax = max(yt(1),yt(2),yt(3))

	    fact = 0.5*i
	    dx = xmax - xmin
	    dy = ymax - ymin
	    xmin = xmin - fact*dx
	    ymin = ymin - fact*dy
	    xmax = xmax + fact*dx
	    ymax = ymax + fact*dy

	    hmed = 0.
	    ihmed = 0
	    do n=1,np
	      if( inquad(xp(n),yp(n),xmin,ymin,xmax,ymax) ) then
	        hmed = hmed + dp(n)
	        ihmed = ihmed + 1
	      end if
	    end do

	    if( ihmed .gt. 0 ) then
	      hmed = hmed / ihmed
	      ok(ie) = .true.
	      hev(ie) = hmed
	    end if
	  else
	    netot = netot + 1
	  end if
	 end do
	 write(6,1000) 'square: ',i,nel-netot,netot,nel
	 i = i + 1
	end do

c-----------------------------------------------------------------
c end up
c-----------------------------------------------------------------

	netot = 0

	do ie=1,nel
	  if( ok(ie) ) then
	    netot = netot + 1
	  else
	    write(6,*) '*** warning interpolq: no depth for ie = ',ie
	    !stop 'error stop interpolq: no depth for element'
	    hev(ie) = flag
	  end if
	  hm3v(:,ie) = hev(ie)
	end do

	write(6,*) 'Elements without depth (end): ',nel-netot,nel

c-----------------------------------------------------------------
c copy interpolated depths to ht
c-----------------------------------------------------------------

	ht(1:nel) = hev(1:nel)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
 1000	format(a,4i10)
	end

c*******************************************************************

	function inquad(xp,yp,xmin,ymin,xmax,ymax)

c is (xp/yp) in quad

	implicit none

	logical inquad
	real xp,yp
	real xmin,ymin,xmax,ymax

	inquad = .false.

	if( xp .lt. xmin ) return
	if( yp .lt. ymin ) return
	if( xp .gt. xmax ) return
	if( yp .gt. ymax ) return

	inquad = .true.

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine interpole(np,xp,yp,dp,nt,xt,yt,at,ht,umfact,nminimum)

c interpolates depth values - exponential interpolation

	use mod_depth
	use basin

	implicit none

	integer np		!single bathymetry points
	real xp(np)
	real yp(np)
	real dp(np)
	integer nt		!points on which to interpolate
	real xt(np)
	real yt(np)
	real at(np)		!area - used as standard deviation
	real ht(np)		!interpolated depth values (return)
	real umfact		!factor for maximum radius
	integer nminimum	!how many points needed for interpolation

	integer k
	integer netot
	real x,y,d
	real depth
	integer iaux,inum,ityp
	integer n,i,ntot
	integer ihmed
        integer nintp,nmin
	integer iloop,nintpol
	real xmin,ymin,xmax,ymax
	real weight,r2,w
        real r2max,sigma2
	real hmed
	real fact,dx,dy
	real area,x0,y0,sig2
	real pi

	logical ok(nel)

	real dist2
	logical inconvex,inquad

c-----------------------------------------------------------------
c initialize
c-----------------------------------------------------------------

	ntot = 0
	do n=1,nt
	  if( ht(n) .gt. -990. ) then
	    ok(n) = .true.
	    ntot = ntot + 1
	  else
	    ok(n) = .false.
	    ht(n) = 0.
	  end if
	end do

	write(6,*) 'Points without value (start): ',nt-ntot,nt
	write(6,*)
	write(6,*) '       loop        done       total       fact'

c-----------------------------------------------------------------
c initial interpolation -> point in element
c-----------------------------------------------------------------

        nmin = nminimum
	if( nmin .le. 0 ) nmin = 1	!at least one point is needed

	nintpol = 0
	iloop = 0
	pi = 4.*atan(1.)

	fact = 2./pi	!start with a radius slightly greater than area

	do while( ntot .lt. nt )

	ntot = 0
	iloop = iloop + 1
	!write(6,*) 'starting new loop on items: ',iloop,fact

	do i=1,nt

	  x0 = xt(i)
	  y0 = yt(i)
	  sig2 = at(i)

	  if( sig2 .le. 0. ) goto 98

          sigma2 = fact*sig2    	 !standard deviation grows
          r2max = (umfact**2)*sigma2     !maximum radius to look for points

	  !write(6,*) i,sig2,sigma2,r2max

	  if( ok(i) ) then
	    ntot = ntot + 1
	  else
	    depth = 0.
	    weight = 0.
            nintp = 0
	    do n=1,np
	      x = xp(n)
	      y = yp(n)
	      d = dp(n)
	      r2 = dist2(x0,y0,x,y)
	      if( r2 .le. r2max ) then
	        w = exp(-r2/(2.*sigma2))
	        depth = depth + d * w
	        weight = weight + w
                nintp = nintp + 1
	      end if
	    end do
	    if( nintp .ge. nmin ) then
	      if( weight .le. 0. ) then
                write(6,*) nintp,weight,r2max,i
                stop 'error stop interpole: zero weight from points'
              end if
	      ht(i) = depth / weight
	      ok(i) = .true.
	      ntot = ntot + 1
	      nintpol = nintpol + 1
	      if( mod(nintpol,100) .eq. 0 ) then
	        write(6,*) iloop,nintpol,nt,fact
	      end if
	    end if
	  end if

	end do

	write(6,*) iloop,nintpol,nt,fact
	!write(6,*) 'Items without depth : ',fact,nt-ntot,nt

	fact = fact * 2.

	end do

c-----------------------------------------------------------------
c end up
c-----------------------------------------------------------------

	ntot = 0

	do i=1,nt
	  if( ok(i) ) then
	    ntot = ntot + 1
	  else
	    write(6,*) i
	    stop 'error stop interpole: no depth for item'
	  end if
	end do

	write(6,*) 'Items without depth (end): ',nt-ntot,nt

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
   98	continue
	write(6,*) i,area
	stop 'error stop interpole: negative area'
	end

c*******************************************************************

	subroutine prepare_on_elem(nt,xt,yt,at,ht,ufact)

c prepares xt,yt,at,ht on nodes

	use mod_depth
	use basin

	implicit none

	integer nt	!total number of points prepared
	real xt(1)	!x-coordinate
	real yt(1)	!y-coordinate
	real at(1)	!area
	real ht(1)	!depth
	real ufact	!sigma or factor to be used

c if ufact > 0: use as factor times area of elem/node
c if ufact < 0: use it directly as sigma

	integer ie,ii,k
	real area,x0,y0,fact,sigma2
	real x(3),y(3)

	real areat

	nt = nel
	fact = ufact*ufact
	sigma2 = fact

	do ie=1,nel

	  do ii=1,3
	    k = nen3v(ii,ie)
	    x(ii) = xgv(k)
	    y(ii) = ygv(k)
	  end do

	  area = areat(x(1),y(1),x(2),y(2),x(3),y(3))
	  call centert(x(1),y(1),x(2),y(2),x(3),y(3),x0,y0)

	  xt(ie) = x0
	  yt(ie) = y0
	  at(ie) = fact*area
	  if( ufact < 0. ) at(ie) = sigma2
	  ht(ie) = hev(ie)

	end do

	end

c*******************************************************************

	subroutine prepare_on_node(nt,xt,yt,at,ht,ufact)

c prepares xt,yt,at,ht on nodes

	use mod_depth
	use basin

	implicit none

	integer nt
	real xt(1)
	real yt(1)
	real at(1)
	real ht(1)
	real ufact	!sigma or factor to be used

c if ufact > 0: use as factor times area of elem/node
c if ufact < 0: use it directly as sigma

	integer ie,ii,k
	real area,fact,sigma2
	real x(3),y(3)

	real areat

	nt = nkn
	fact = ufact*ufact
	sigma2 = fact

	do k=1,nkn
	  at(k) = 0.
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    x(ii) = xgv(k)
	    y(ii) = ygv(k)
	  end do
	  area = areat(x(1),y(1),x(2),y(2),x(3),y(3))
	  do ii=1,3
	    k = nen3v(ii,ie)
	    at(k) = at(k) + area / 3.
	  end do
	end do

	do k=1,nkn
	  xt(k) = xgv(k)
	  yt(k) = ygv(k)
	  at(k) = fact*at(k)
	  if( ufact < 0. ) at(k) = sigma2
	  ht(k) = hkv(k)
	end do

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine interpola(np,xp,yp,dp,ap,nt,xt,yt,at,ht)

c interpolates depth values with auto-correlation

	use mod_depth
	use basin

	implicit none

	integer np
	real xp(np)
	real yp(np)
	real dp(np)
	real ap(np)		!sigma**2 for bathy points
	integer nt
	real xt(np)
	real yt(np)
	real at(np)
	real ht(np)

	include 'param.h'




	integer k
	integer netot
	real x,y,d
	real depth
	integer iaux,inum,ityp
	integer n,i,ntot
	integer ihmed
        integer nintp,nmin
	integer nintpol
	real xmin,ymin,xmax,ymax
	real weight,r2,w
        real r2max,sigma2
	real hmed
	real fact,dx,dy
	real area,x0,y0
	real pi
	real ufact,umfact,a

	logical ok(nel)

	real dist2
	logical inconvex,inquad

c-----------------------------------------------------------------
c initialize
c-----------------------------------------------------------------

	ntot = 0
	do n=1,nt
	  if( ht(n) .gt. -990. ) then
	    ok(n) = .true.
	    ntot = ntot + 1
	  else
	    ok(n) = .false.
	    ht(n) = 0.
	  end if
	end do

	write(6,*) 'Points without depth (start): ',nt-ntot,nt

c-----------------------------------------------------------------
c initial interpolation -> point in element
c-----------------------------------------------------------------

	nintpol = 0
	ntot = 0

	do i=1,nt

	  x0 = xt(i)
	  y0 = yt(i)

	  if( ok(i) ) then
	    ntot = ntot + 1
	  else
	    depth = 0.
	    weight = 0.
            nintp = 0
	    do n=1,np
	      x = xp(n)
	      y = yp(n)
	      d = dp(n)
	      a = ap(n)
	      r2 = dist2(x0,y0,x,y)
	      !if (a == 0) a = 0.0001  	!mbj
	      cycle			!do not use this point
	      w = exp(-r2/(2.*a))
	      depth = depth + d * w
	      weight = weight + w
              nintp = nintp + 1
	    end do
	    if( weight .le. 0. ) then
              write(6,*) nintp,weight,r2max,i
              stop 'error stop interpole: zero weight from points'
            end if
	    ht(i) = depth / weight
	    ok(i) = .true.
	    ntot = ntot + 1
	    nintpol = nintpol + 1
	    if( mod(nintpol,100) .eq. 0 ) then
	      write(6,*) 'items interpolated: ',nintpol,nt
	    end if
	  end if

	end do

	write(6,*) 'Items without depth : ',nt-ntot,nt

c-----------------------------------------------------------------
c end up
c-----------------------------------------------------------------

	ntot = 0

	do i=1,nt
	  if( ok(i) ) then
	    ntot = ntot + 1
	  else
	    write(6,*) i
	    stop 'error stop interpole: no depth for item'
	  end if
	end do

	write(6,*) 'Items without depth (end): ',nt-ntot,nt

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
   98	continue
	write(6,*) i,area
	stop 'error stop interpole: negative area'
	end

c*******************************************************************

	subroutine make_auto_corr(np,xp,yp,dp,ap,ufact)

c computes minimum distance to next point
c
c this is used as sigma for the interpolation
c there are margins for improving this

	implicit none

	integer np
	real xp(np)
	real yp(np)
	real dp(np)
	real ap(np)		!computed sigma**2 values for each data point
	real ufact		!extra factor for size of sigma

	integer i,j
	real d
	double precision dm

	real dist2 

	if( ufact <= 0. ) then
	  write(6,*) 'ufact = ',ufact
	  stop 'error stop make_auto_corr: ufact'
	end if

	write(6,*) 'computing auto correlation... ',np

	do i=1,np
	  dm = 1.e+30
	  do j=1,np
	    if( i .ne. j ) then
	      d = dist2(xp(i),yp(i),xp(j),yp(j))
	      if( d .lt. dm ) dm = d
	    end if
	  end do
	  ap(i) = dm * ufact*ufact
	  !write(6,*) i,ap(i),dp(i)
	end do

	write(6,*) 'auto correlation finished ',np

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

