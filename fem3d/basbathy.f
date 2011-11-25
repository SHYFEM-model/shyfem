c
c $Id: basbathy.f,v 1.7 2010-02-26 17:35:06 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	completely restructured
c 04.06.1999	ggu	new statistics are computed
c 08.09.2003	ggu	mode 5 -> write depth values from elements
c 23.09.2004    ggu     interpolq() changed for bathy interpolation
c 02.10.2004    ggu     interpole() for exponential interpolation
c 12.05.2005    ggu     pass hmin to interpolation functions
c 06.04.2005    ggu     read param.h
c 24.04.2009	ggu	new call to rdgrd()
c 21.05.2009	ggu	restructured to allow for nodal interpolation
c 16.12.2010	ggu	bug fix in copy_depth()
c
c****************************************************************

        program basbathy

c performs bathymetry interpolation in basin
c
c takes care of lat/lon coordinates

	implicit none

	include 'param.h'

	integer ndim
	parameter(ndim=1200000)
	real xp(ndim)
	real yp(ndim)
	real dp(ndim)

        character*80 descrp
        common /descrp/ descrp

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        real xgv(nkndim), ygv(nkndim)
        real hm3v(3,neldim)
        common /xgv/xgv, /ygv/ygv
        common /hm3v/hm3v

        integer nen3v(3,neldim)
        integer ipev(neldim), ipv(nkndim)
        integer iarv(neldim)
        common /nen3v/nen3v
        common /ipev/ipev, /ipv/ipv
        common /iarv/iarv

        real hev(neldim)
        common /hev/hev
        real hkv(nkndim)
        common /hkv/hkv

        real raux(neldim)
        integer iaux(neldim)
        integer ipaux(nkndim)

	include 'evmain.h'

        character*40 bfile,gfile,nfile
        character*60 line
	integer node,nit
	integer mode,np,n,i
        integer ner,nco,nknh,nelh,nli
	integer nlidim,nlndim
	integer ike,idepth
	real ufact,umfact
	real f(5)
	logical bstop
	integer iscanf

c-----------------------------------------------------------------
c what to do
c-----------------------------------------------------------------

        write(6,*)
        write(6,*) 'I need the name of the grid file '
        write(6,*)
	write(6,*) 'Enter file name: '
	read(5,'(a)') gfile
        if( gfile .eq. ' ' ) stop
	write(6,*) 'grid is read from file : ', gfile
        write(6,*)

        write(6,*)
        write(6,*) 'I need the name of the bathymetry data file '
        write(6,*)
	write(6,*) 'Enter file name: '
	read(5,'(a)') bfile
        if( bfile .eq. ' ' ) stop
	write(6,*) 'Bathymetry is read from file : ', bfile
        write(6,*)

        write(6,*)
        write(6,*) 'Two different algorithms are available:'
        write(6,*) '  1   exponential interpolation (default)'
        write(6,*) '  2   uniform interpolation on squares'
        write(6,*)
	write(6,*) 'Enter choice: '
	read(5,'(i10)') mode
	if( mode .ne. 2 ) mode = 1
	write(6,*) 'Mode is : ', mode

        write(6,*)
        write(6,*) 'If there are some values missing you can:'
        write(6,*) '  1   interpolate on missing depth values (default)'
        write(6,*) '  2   interpolate on all elements/nodes'
        write(6,*)
	write(6,*) 'Enter choice: '
	read(5,'(i10)') idepth
	if( idepth .ne. 2 ) idepth = 1
	write(6,*) 'Choice is : ', idepth

	ike = 1
	ufact = 1.
	umfact = 2.

	if( mode .eq. 1 ) then

        write(6,*)
        write(6,*) 'For the exponential algorithm you can:'
        write(6,*) '  1   interpolate on elements (default)'
        write(6,*) '  2   interpolate on nodes'
        write(6,*)
	write(6,*) 'Enter choice: '
	read(5,'(i10)') ike
	if( ike .ne. 2 ) ike = 1
	write(6,*) 'Choice is : ', ike

        write(6,*)
	write(6,*) 'Enter parameters for expontential interpolation:'
        write(6,*)
	write(6,*) 'The std deviation is about the size of the elements'
	write(6,*) 'With ufact you can ultimately correct it (default=1)'
	write(6,*) 'The maximum radius is twice the standard deviation'
	write(6,*) 'With umfact you can correct it (default=2)'
        write(6,*)
	write(6,*) 'Enter params ufact and umfact (<CR> for default): '
        read(5,'(a)') line
        n = iscanf(line,f,2)
	if( n .lt. 0 .or. n .gt. 2 ) goto 95
        if( n .gt. 0 ) ufact = f(1)
        if( n .gt. 1 ) umfact = f(2)
        write(6,*) 'ufact,umfact :',ufact,umfact
        write(6,*)

	end if

c-----------------------------------------------------------------
c read in bathymetry file
c-----------------------------------------------------------------

	np = ndim
	call readbat(bfile,np,xp,yp,dp)

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        ner = 6
        bstop = .false.

        nlidim = 0
        nlndim = 0
        call rdgrd(
     +                   gfile
     +                  ,bstop
     +                  ,nco,nkn,nel,nli
     +                  ,nkndim,neldim,nlidim,nlndim
     +                  ,ipv,ipev,iaux
     +                  ,iaux,iarv,iaux
     +                  ,hkv,hev,raux
     +                  ,xgv,ygv
     +                  ,nen3v
     +                  ,iaux,iaux
     +                  )

        if( bstop ) stop 'error stop rdgrd'

c        call rdgrd(gfile,ner,bstop,nco,nkn,nknh,nel,nelh,nli
c     +                  ,nkndim,neldim,nliread
c     +                  ,ipv,ipev,ianv,iarv,nen3v,xgv,ygv,hev,hkv)
c
c        call ex2in(nkn,nel,ipv,ipaux,nen3v,ner,bstop)

        call ex2in(nkn,3*nel,nlidim,ipv,ipaux,nen3v,iaux,bstop)
        if( bstop ) stop 'error stop ex2in'

c-----------------------------------------------------------------
c handling depth and coordinates
c
c the program automatically checks if we have lat/lon coordinates
c if for some ???ggu???
c-----------------------------------------------------------------

	call set_dist(0)
	call set_coords_ev(0)

	call check_coords			!sets lat/lon flag
	call set_depth(idepth,nknh,nelh)

c-----------------------------------------------------------------
c general info
c-----------------------------------------------------------------

        write(6,*)
        write(6,*) ' nkn  = ',nkn, '  nel  = ',nel
        write(6,*) ' nknh = ',nknh,'  nelh = ',nelh
        write(6,*)

c-----------------------------------------------------------------
c node_test
c-----------------------------------------------------------------

	call node_test
	call set_ev

c-----------------------------------------------------------------
c interpolate
c-----------------------------------------------------------------

	if( mode .eq. 1 ) then
	  call interpole(np,xp,yp,dp,ike,ufact,umfact)
        else if( mode .eq. 2 ) then
	  call interpolq(np,xp,yp,dp)
        else
          write(6,*) 'wrong choice for mode : ',mode
          stop 'error stop'
	end if

	call copy_depth(ike)	!copyt to nodes/elements

c-----------------------------------------------------------------
c write
c-----------------------------------------------------------------

	nfile = 'basbathy.grd'
	open(1,file=nfile,status='unknown',form='formatted')
	call wrgrd(1,ike)
	close(1)
        write(6,*) 'file has been written to ',nfile

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	stop
   95	continue
	write(6,*) n,(f(i),i=1,n)
	write(6,*) line
	stop 'error stop basbathy: error in parameters'
	end

c*******************************************************************

	subroutine copy_depth(ike)

c copies depth values from elems/nodes to nodes/elems

	implicit none

	integer ike

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,neldim)
        common /nen3v/nen3v

        real hev(neldim)
        common /hev/hev
        real hkv(nkndim)
        common /hkv/hkv

	integer k,ie,ii
	real depth
	integer ic(nkndim)

	if( ike .eq. 1 ) then		!elementwise

	  do k=1,nkn
	    ic(k) = 0
	    hkv(k) = 0.
	  end do

	  do ie=1,nel
	    do ii=1,3
	      k = nen3v(ii,ie)
	      hkv(k) = hkv(k) + hev(ie)
	      ic(k) = ic(k) + 1			!BUG - this was missing
	    end do
	  end do

	  do k=1,nkn
	    hkv(k) = hkv(k) / ic(k)
	  end do

	else				!nodewise

	  do ie=1,nel
	    depth = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      depth = depth + hkv(k)
	    end do
	    hev(ie) = depth / 3.
	  end do

	end if

	end

c*******************************************************************

	subroutine set_depth(idepth,nknh,nelh)

c handles depth values

	implicit none

	integer idepth
	integer nknh,nelh

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real hev(neldim)
        common /hev/hev
        real hkv(nkndim)
        common /hkv/hkv

	integer k,ie
	real flag

	flag = -999.

	if( idepth .eq. 2 ) then
	  do k=1,nkn
	    hkv(k) = flag
	  end do
	  do ie=1,nel
	    hev(ie) = flag
	  end do
	end if
	
	nknh = 0
	nelh = 0

	do k=1,nkn
	  if( hkv(k) .gt. -990 ) nknh = nknh + 1
	end do
	do ie=1,nel
	  if( hev(ie) .gt. -990 ) nelh = nelh + 1
	end do

	end

c*******************************************************************

	subroutine check_coords

c checks if coordinates are lat/lon

	implicit none

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

	integer k,isphe
	real xmin,xmax,ymin,ymax

	xmin = xgv(1)
	xmax = xgv(1)
	ymin = ygv(1)
	ymax = ygv(1)

	do k=1,nkn
	  xmin = min(xmin,xgv(k))
	  xmax = max(xmax,xgv(k))
	  ymin = min(ymin,ygv(k))
	  ymax = max(ymax,ygv(k))
	end do

	isphe = 1
	if( xmin .lt. -180. .or. xmax .gt. 360. ) isphe = 0
	if( ymin .lt. -180. .or. ymax .gt. 180. ) isphe = 0
	call set_dist(isphe)

	call set_coords_ev(isphe)
	write(6,*) 'setting for coordinates: ',isphe
	if( isphe .ne. 0 ) write(6,*) 'using lat/lon coordinates'

	end

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

        function areatr(ie)

c determination of area of element
c
c ie            number of element (internal)
c areatr        element area (return value)

	real areatr
	integer ie

	include 'param.h'

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

        integer nen3v(3,neldim)
        common /nen3v/nen3v

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

	subroutine wrgrd(iunit,ike)

c writes grd file from bas

	implicit none

	integer iunit
	integer ike

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

        integer nen3v(3,neldim)
        integer ipev(neldim), ipv(nkndim)
        integer iarv(neldim)
        common /nen3v/nen3v
        common /ipev/ipev, /ipv/ipv
        common /iarv/iarv

        real hev(neldim)
        common /hev/hev
        real hkv(nkndim)
        common /hkv/hkv

	integer k,ie,ii

	do k=1,nkn
	  if( ike .eq. 1 ) then
	    write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k)
	  else
	    write(iunit,1000) 1,ipv(k),0,xgv(k),ygv(k),hkv(k)
	  end if
	end do

	write(iunit,*)

	do ie=1,nel
	  if( ike .eq. 1 ) then
	    write(iunit,1100) 2,ipev(ie),iarv(ie)
     +		,3,(ipv(nen3v(ii,ie)),ii=1,3),hev(ie)
	  else
	    write(iunit,1100) 2,ipev(ie),iarv(ie)
     +		,3,(ipv(nen3v(ii,ie)),ii=1,3)
	  end if
	end do

	return
 1000	format(i1,2i10,3e16.8)
 1100	format(i1,2i10,i4,3i10,e16.8)
	end

c*******************************************************************

	subroutine interpolq(np,xp,yp,dp)

c interpolates depth values

	implicit none

	integer np
	real xp(np)
	real yp(np)
	real dp(np)

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v
	real hev(neldim)
        common /hev/hev

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer ie,ii,k
	integer netot
	real x,y,d
	real depth
	integer iaux,inum,ityp
	integer n,i
	integer ihmed
	real xt(3), yt(3)
	real xmin,ymin,xmax,ymax
	real hmed
	real fact,dx,dy
	integer ihev(neldim)
	logical bmin
	logical ok(neldim)
	logical inconvex,inquad

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

	write(6,*) 'Elements without depth (start): ',nel-netot,nel

c-----------------------------------------------------------------
c initial interpolation -> point in element
c-----------------------------------------------------------------

	do n=1,np
	  x = xp(n)
	  y = yp(n)
	  d = dp(n)

	  if( mod(n,100) .eq. 0 ) then
	    !write(6,*) n,(100.*n)/np,x,y,d
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

	write(6,*) 'Elements without depth (convex): ',nel-netot,nel

c-----------------------------------------------------------------
c next interpolation -> point in square
c-----------------------------------------------------------------

	i = 0
	do while( netot .lt. nel )
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
	 write(6,*) 'Elements without depth (quad): ',i,nel-netot,nel
	 i = i + 1
	end do

c-----------------------------------------------------------------
c end up
c-----------------------------------------------------------------

	netot = 0

	do ie=1,nel
	  if( ok(ie) ) then
	    netot = netot + 1
	    hmed = hev(ie)
	  else
	    stop 'error stop interpolq: no depth for element'
	  end if
	  do ii=1,3
	    hm3v(ii,ie) = hmed
	  end do
	end do

	write(6,*) 'Elements without depth (end): ',nel-netot,nel

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

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

	subroutine node_test

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer ie,ii,k

	write(6,*) 'node_testing ... ',nel,nkn
	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( k .le. 0 ) write(6,*) ie,ii,k
	    iii = mod(ii,3) + 1
	    k1 = nen3v(iii,ie)
	    if( k .eq. k1 ) write(6,*) ie,(nen3v(iii,ie),iii=1,3)
	  end do
	end do
	write(6,*) 'end of node_testing ... '

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

	subroutine interpole(np,xp,yp,dp,ike,ufact,umfact)

c interpolates depth values

	implicit none

	integer np
	real xp(np)
	real yp(np)
	real dp(np)
	integer ike
	real ufact,umfact

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v
	real hev(neldim)
        common /hev/hev
	real hkv(nkndim)
        common /hkv/hkv

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer k
	integer netot
	real x,y,d
	real depth
	integer iaux,inum,ityp
	integer n,i,nt,ntot
	integer ihmed
        integer nintp
	integer iloop,nintpol
	real xmin,ymin,xmax,ymax
	real weight,r2,w
        real r2max,sigma2
	real hmed
	real fact,dx,dy
	real area,x0,y0
	real pi

	real xt(neldim)
	real yt(neldim)
	real at(neldim)
	real ht(neldim)
	integer ic(neldim)
	logical ok(neldim)

	real dist2
	logical inconvex,inquad

c-----------------------------------------------------------------
c initialize
c-----------------------------------------------------------------


	if( ike .eq. 1 ) then		!elementwise
	  call prepare_on_element(nt,xt,yt,at,ht)
	else
	  call prepare_on_node(nt,xt,yt,at,ht)
	end if

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

	!ufact = 1.
	!umfact = 2.

	nintpol = 0
	iloop = 0
	pi = 4.*atan(1.)
	fact = 2.
	fact = 2./pi	!start with a radius slightly greater than area
	fact = fact * ufact

	do while( ntot .lt. nt )

	ntot = 0
	iloop = iloop + 1
	write(6,*) 'starting new loop on items: ',iloop,fact

	do i=1,nt

	  x0 = xt(i)
	  y0 = yt(i)
	  area = at(i)

	  if( area .le. 0. ) goto 98

          sigma2 = fact*area    	 !standard deviation grows
          r2max = (umfact**2)*sigma2     !maximum radius to look for points

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
	      !r2 = (x-x0)*(x-x0)+(y-y0)*(y-y0)
	      if( r2 .le. r2max ) then
	        w = exp(-r2/(2.*sigma2))
	        depth = depth + d * w
	        weight = weight + w
                nintp = nintp + 1
	      end if
	    end do
	    if( nintp .gt. 0 ) then
	      if( weight .le. 0. ) then
                write(6,*) nintp,weight,r2max,i
                stop 'error stop interpole: zero weight from points'
              end if
	      ht(i) = depth / weight
	      ok(i) = .true.
	      ntot = ntot + 1
	      nintpol = nintpol + 1
	      if( mod(nintpol,100) .eq. 0 ) then
	        write(6,*) 'items interpolated: ',iloop,nintpol,nt
	      end if
	    end if
	  end if

	end do

	write(6,*) 'Items without depth : ',fact,nt-ntot,nt

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
	  if( ike .eq. 1 ) then
	    hev(i) = ht(i)
	  else
	    hkv(i) = ht(i)
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

	subroutine prepare_on_element(nt,xt,yt,at,ht)

	implicit none

	integer nt
	real xt(1)
	real yt(1)
	real at(1)
	real ht(1)

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v
	real hev(neldim)
        common /hev/hev
	real hkv(nkndim)
        common /hkv/hkv

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer ie,ii,k
	real area,x0,y0
	real xaux(3),yaux(3)

	nt = nel

	do ie=1,nel

	  do ii=1,3
	    k = nen3v(ii,ie)
	    xaux(ii) = xgv(k)
	    yaux(ii) = ygv(k)
	  end do
	  call triab(xaux,yaux,area,x0,y0)

	  xt(ie) = x0
	  yt(ie) = y0
	  at(ie) = area
	  ht(ie) = hev(ie)

	end do

	end

c*******************************************************************

	subroutine prepare_on_node(nt,xt,yt,at,ht)

	implicit none

	integer nt
	real xt(1)
	real yt(1)
	real at(1)
	real ht(1)

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v
	real hev(neldim)
        common /hev/hev
	real hkv(nkndim)
        common /hkv/hkv

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer ie,ii,k
	real area,x0,y0
	real xaux(3),yaux(3)

	nt = nkn

	do k=1,nkn
	  at(k) = 0.
	  ht(k) = 0.
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    xaux(ii) = xgv(k)
	    yaux(ii) = ygv(k)
	  end do
	  call triab(xaux,yaux,area,x0,y0)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    at(k) = at(k) + area / 3.
	    ht(k) = ht(k) + 1.
	  end do
	end do

	do k=1,nkn
	  xt(k) = xgv(k)
	  yt(k) = ygv(k)
	  at(k) = at(k) / ht(k)
	  ht(k) = hkv(k)
	end do

	end

c*******************************************************************

