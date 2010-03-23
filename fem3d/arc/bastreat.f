c
c $Id: bastreat.f,v 1.10 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	completely restructured
c 04.06.1999	ggu	new statistics are computed
c 08.09.2003	ggu	mode 5 -> write depth values from elements
c 23.09.2004    ggu     interpolq() changed for bathy interpolation
c 02.10.2004    ggu     interpole() for exponential interpolation
c 01.11.2004    ggu     whole program simplyfied
c 06.12.2008    ggu     smoothing introduced
c 06.04.2005    ggu     read param.h
c
c****************************************************************

        program bastreat

c performs modifications on basin

	implicit none

	include 'param.h'

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

	include 'evmain.h'

        character*20 file,dfile
	integer node,nit
	integer mode,idep
	logical bwrite
	integer iapini

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        if( iapini(1,nkndim,neldim,0) .le. 0 ) stop

c-----------------------------------------------------------------
c what to do
c-----------------------------------------------------------------

	bwrite = .false.
        dfile = 'depth.grd'
        file = ' '
	mode = 1

        write(6,*) 'What type of interpolation?'
        write(6,*)
        write(6,*) '1   uniform interpolation on squares'
        write(6,*) '2   exponential interpolation with max radius'
        write(6,*) '3   smoothing with default values'
        write(6,*)
	write(6,*) 'Enter choice for interpolation: '
	read(5,'(i10)') mode
	write(6,*) 'Interpolation is : ', mode

	if( mode .le. 0 ) stop

	if( mode .le. 2 ) then

        write(6,*) 'Where to interpolate?'
        write(6,*)
        write(6,*) '1   interpolate in the whole basin'
        write(6,*) '2   interpolate only where no depth data available'
        write(6,*)
	write(6,*) 'Enter choice: '
	read(5,'(i10)') idep
	write(6,*) 'Your choice is : ', idep

	if( idep .eq. 1 ) call delete_depth(nel,hm3v)

        write(6,*)
        write(6,*) 'I need the name of the bathymetry data file '
        write(6,*) '   default: ',dfile
        write(6,*)
	write(6,*) 'Enter file name: '
	read(5,'(a)') file
        if( file .eq. ' ' ) file = dfile
	write(6,*) 'Bathymetry is read from file : ', file
        write(6,*)

	end if

c-----------------------------------------------------------------
c general info
c-----------------------------------------------------------------

        write(6,*)
        write(6,*) ' nkn = ',nkn,'  nel = ',nel
        write(6,*) ' mbw = ',mbw,'  ngr = ',ngr
        write(6,*)
        write(6,*) ' dcor = ',dcor,'  dirn = ',dirn
        write(6,*)

c-----------------------------------------------------------------
c test
c-----------------------------------------------------------------

	call test
	call sp110a

c-----------------------------------------------------------------
c interpolate
c-----------------------------------------------------------------

	if( mode .eq. 1 ) then
	  bwrite = .true.
	  call interpolq(file)
        else if( mode .eq. 2 ) then
	  bwrite = .true.
	  call interpole(file)
        else if( mode .ge. 3 .and. mode .le. 3 ) then
	  bwrite = .true.
          call smooth_bathy
        else if( mode .ge. 92 .and. mode .le. 95 ) then
	  bwrite = .false.
          call other_modes(mode)
        else
          write(6,*) 'wrong choice for mode : ',mode
          stop 'error stop'
	end if

c-----------------------------------------------------------------
c write
c-----------------------------------------------------------------

	if( bwrite ) then
	  open(1,file='new.grd',status='unknown',form='formatted')
	  call wrgrd(1)
	  close(1)
          write(6,*) 'file has been written to new.grd'
	end if

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

	subroutine delete_depth(nel,hm3v)

c delete depth values from basin

        implicit none

	integer nel
        real hm3v(3,1)

	integer ie,ii
	real flag

	flag = -999.

	write(6,*) 'Deleting available depth values...'

	do ie=1,nel
	  do ii=1,3
	    hm3v(ii,ie) = flag
	  end do
	end do
	
	end

c*******************************************************************

        subroutine other_modes(mode)

c some other modes used long time ago

        implicit none

        integer mode

	if( mode .eq. 92 ) then
	  call neibors(43,2)
	  call neibors(112,2)
	  call neibors(145,2)
	  call neibors(196,2)
	  call neibors(786,2)
	  call neibors(970,2)
	  call neibors(1096,2)
	  call neibors(1497,2)
	  call neibors(1773,2)
        else if( mode .eq. 93 ) then
	  call wnodes('nodes_s.lst')
        else if( mode .eq. 94 ) then
	  call wnodeef('nodes_s.lst')
        else if( mode .eq. 95 ) then
	  call wrdepth(4)
	end if

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

	subroutine interpol(file)

c interpolates depth values

	implicit none

	character*(*) file

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer ie,ii,k
	real x,y,d
	real hmin
	integer iaux,inum,ityp
	real xt(3), yt(3)
	real hev(neldim)
	integer ihev(neldim)
	logical inconvex

	hmin = 1.

	do ie=1,nel
	  hev(ie) = 0.
	  ihev(ie) = 0
	end do

	open(1,file=file,status='old',form='formatted')

    1	continue
	  read(1,*,end=99) iaux,inum,ityp,x,y,d
	  write(6,*) inum,x,y,d

	  do ie=1,nel
	    do ii=1,3
	      k = nen3v(ii,ie)
	      xt(ii) = xgv(k)
	      yt(ii) = ygv(k)
	    end do
	    if( inconvex(3,xt,yt,x,y) ) then
	      hev(ie) = hev(ie) + d
	      ihev(ie) = ihev(ie) + 1
	    end if
	  end do

	  goto 1
   99	continue
	close(1)

	do ie=1,nel
	  if( ihev(ie) .gt. 0 ) then
	    hev(ie) = hev(ie) / ihev(ie)
	    if( hev(ie) .lt. hmin ) hev(ie) = hmin
	  end if
	  do ii=1,3
	    hm3v(ii,ie) = hev(ie)
	  end do
	end do

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

	subroutine wrgrd(iunit)

c writes grd file from bas

	implicit none

	integer iunit

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

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

	subroutine interpolq(file)

c interpolates depth values
c
c uniform interpolation on squares

	implicit none

	character*(*) file

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer ie,ii,k
	integer netot
	real x,y,d
	real hmin
	real depth
	integer iaux,inum,ityp
	integer n,np,i
	integer ihmed
	real xt(3), yt(3)
	real xmin,ymin,xmax,ymax
	real hmed
	real fact,dx,dy
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
	bmin = .false.
	hmin = 1.
	hmin = 0.5

	write(6,*) 'starting uniform interpolation on squares'

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

	do n=1,np
	  x = xp(n)
	  y = yp(n)
	  d = dp(n)

	  if( mod(n,1000) .eq. 0 ) then
	    !write(6,*) n,x,y,d,(100.*n)/np,' %'
	    write(6,3000) n,x,y,d,(100.*n)/np,' %'
 3000	    format(i10,2e14.6,f10.2,f10.1,a)
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

	subroutine test

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer ie,ii,k

	write(6,*) 'testing ... ',nel,nkn
	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( k .le. 0 ) write(6,*) ie,ii,k
	    iii = mod(ii,3) + 1
	    k1 = nen3v(iii,ie)
	    if( k .eq. k1 ) write(6,*) ie,(nen3v(iii,ie),iii=1,3)
	  end do
	end do
	write(6,*) 'end of testing ... '

	end

c*******************************************************************

	subroutine neibors(node,nit)

c gets neibors of node node (nit iterations)

	implicit none

	integer node
	integer nit

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,neldim)
        common /nen3v/nen3v

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

	implicit none

	character*(*) file

	include 'param.h'

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

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

	implicit none

	character*(*) file

	include 'param.h'

	integer ndim
	parameter (ndim=1000)

	integer nodes(nkndim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

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

	implicit none

	include 'param.h'

	integer ntot	!number of subdivisions in element

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        real hm3v(3,neldim)
        common /xgv/xgv, /ygv/ygv
        common /hm3v/hm3v

        integer nen3v(3,neldim)
        common /nen3v/nen3v

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
c
c exponential interpolation with max radius

	implicit none

	character*(*) file

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v

        integer nen3v(3,neldim)
        common /nen3v/nen3v

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
	bmin = .false.
	hmin = 1.
	hmin = 0.5

        write(6,*) 'starting exponential interpolation with max radius'

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

	subroutine smooth_bathy

c interpolates depth values

	implicit none

	include 'param.h'
	include 'evmain.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv
        real hm3v(3,neldim)
        common /hm3v/hm3v

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	integer ie,ii,k,i
	integer niter
	real x,y,d,h,hold,ao
	real hmin
	real alpha,beta
        real h1,h2,a1,a2
	integer iaux,inum,ityp
	real xt(3), yt(3)
	real hev(neldim), hkv(nkndim), v1v(nkndim)
	integer ihev(neldim)
	logical inconvex

c--------------------------------------------------------------
c set parameters
c--------------------------------------------------------------

	niter = 50
	alpha = 0.5
        h1 = 500
        a1 = 0.5
        h2 = 50
        a2 = 0.1

	write(6,*) 'smoothing bathymetry: ',niter,alpha

c--------------------------------------------------------------
c compute hev
c--------------------------------------------------------------

	do ie=1,nel
	  h = 0.
	  do ii=1,3
	    h = h + hm3v(ii,ie)
	  end do
	  hev(ie) = h/3.
	end do

c--------------------------------------------------------------
c iterate over smoother
c--------------------------------------------------------------

	do i=1,niter

	  write(6,*) 'pass ',i,' of ',niter

c	  -----------------------------------------------
c	  compute values at nodes (averages of element)
c	  -----------------------------------------------

	  do k=1,nkn
	    hkv(k) = 0.
	    v1v(k) = 0.
	  end do

	  do ie=1,nel
	    ao = ev(10,ie)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      hkv(k) = hkv(k) + ao * hev(ie)
	      v1v(k) = v1v(k) + ao
	    end do
	  end do

	  do k=1,nkn
	    hkv(k) = hkv(k) / v1v(k)
	  end do

c	  -----------------------------------------------
c	  average to element
c	  -----------------------------------------------

	  do ie=1,nel
	    h = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      h = h + hkv(k)
	    end do
	    h = h / 3.
	    hold = hev(ie)

            beta = (hold-h1)/(h2-h1)
            beta = max(beta,0.)
            beta = min(beta,1.)
            beta = a1 + beta * (a2-a1)
            !beta = alpha
            !write(6,*) ie,hold,beta

	    hev(ie) = (1.-beta) * hold + beta * h
	  end do

	end do

c-----------------------------------------------------------------
c put back into hm3v
c-----------------------------------------------------------------

	do ie=1,nel
	  do ii=1,3
	    hm3v(ii,ie) = hev(ie)
	  end do
	end do

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

