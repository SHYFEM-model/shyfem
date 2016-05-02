c
c $Id: subnev.f,v 1.15 2010-02-16 16:21:37 georg Exp $
c
c ev routines
c
c contents :
c
c subroutine set_ev				set up ev vector
c subroutine check_ev				tests if ev is set up
c subroutine cpp(x,y,rlambda,phi,rlambda0,phi0)	transforms (lon,lat) to cart
c subroutine adjust_bc(v1,v2,v3)		adjusts b/c so that sum = 0
c function area_elem(ie)			returns area of element ie
c function aomega_elem(ie)			returns aomega of element ie
c
c revision log :
c
c 31.05.1997	ggu	unnecessary routines deleted
c 27.06.1997	ggu	ev routines into own file
c 12.11.2001	ggu	cosmetic changes
c 10.08.2003	ggu	new routine check_ev
c 14.08.2003	ggu	sp110a -> set_ev
c 27.08.2007	ccf	isphe for spherical coordinate system
c 24.06.2008	ggu	compute and store also distances beteween nodes
c 18.11.2008	ggu	new routine adjust_bc() to adjust sum to 0
c 19.11.2008	ggu	helper routines to compute area/aomega
c 22.05.2009	ggu	new routine set_coords_ev() and ev_blockdata
c 12.06.2009	ggu	more stable computation of area (bug_f_64bit)
c 05.02.2010	ggu	bug fix for aj (division with 24 too early)
c 14.04.2010	ggu	new routines get_coords_ev() and check_spheric_ev()
c 07.05.2010	ggu	initialization of ev routines
c 25.01.2011	ggu	default to lat/lon if small coordinates are given
c 28.01.2011	ggu	new entry in ev for distance of nodes (17-19)
c 23.03.2011	ggu	better set-up for isphe_ev
c 24.11.2011	ggu	better routines to handle spherical coordinates
c 30.08.2012	ggu	new routine is_spherical()
c 30.09.2013	ggu	new routines gradient_normal/tangent
c 30.10.2014	ggu	new routine compute_distance()
c 21.01.2015	ggu	new routine compute_cartesian_coords()
c 31.03.2015	ggu	new routines for internal coordinates
c 19.04.2015	ggu	new routines for internal coordinates (distance)
c 22.04.2015	ggu	bug fix in xi2xy - x/y were exchanged
c 10.10.2015	ggu	bug fix in adjust_bc() - brown paper bag bug
c 16.11.2015	ggu	new routine adjust_xi()
c 15.02.2016	ggu	more debug code, assure xi is in bounds
c
c***********************************************************

c isphe_ev:
c
c -1	value not given -> try to determine automatically
c  0	cartesian
c  1	spherical (lat/lon)

!==================================================================
        module evgeom
!==================================================================

        implicit none

	integer, parameter :: evdim = 19

	integer, save, private :: nel_alloc = 0
	integer, save :: isphe_ev = -1
	logical, save :: init_ev = .false.

	logical, save :: bdebug_internal = .false.

	double precision, allocatable :: ev(:,:)

!==================================================================
        contains
!==================================================================

	subroutine ev_init(nel)

	integer nel

	if( nel == nel_alloc ) return

	if( nel_alloc > 0 ) then
	  deallocate(ev)
	end if

	nel_alloc = nel

	if( nel > 0 ) then
	  allocate(ev(evdim,nel))
	end if

	end subroutine ev_init

!==================================================================
        end module evgeom
!==================================================================

	subroutine set_ev

c set up ev vector
c
c double precision version
c
c revised on 31.08.88 by ggu (czv containes real chezy)
c revised on 25.11.88 by ggu (czv eliminated)
c revised on 28.01.92 by ggu (double precision, implicit none)

	use basin
	use evgeom

	implicit none

	include 'param.h'

	logical bdebug
	integer iedebug
	integer ie,i,kn1,kn2,kn3
	integer isphe
        double precision one,two,four,twofour,rad
	double precision x1,x2,x3,y1,y2,y3
	double precision a1,a2,a3,b1,b2,b3,c1,c2,c3
	double precision aj,aji,pi
	double precision s1,s2,s3,ss1,ss2,ss3
	double precision d1,d2,d3
	double precision dd1,dd2,dd3

	double precision xm,ym,maxmax
	double precision xx1,xx2,xx3,yy1,yy2,yy3

	double precision xlon1,ylat1,xlon2,ylat2,xlon3,ylat3	!lat/long [rad]
	double precision dlat0,dlon0			!center of projection

	if( .not. basin_has_basin() ) goto 97

	call check_spheric_ev	!checks and sets isphe_ev and init_ev
	call get_coords_ev(isphe)

	iedebug = 77
	iedebug = 0

        one = 1.
        two = 2.
        four = 4.
        twofour = 24.

	maxmax = 0.

	pi=four*atan(one)
        rad = 180./pi

	do ie=1,nel

	bdebug = ( iedebug == ie )

	kn1=nen3v(1,ie)
	kn2=nen3v(2,ie)
	kn3=nen3v(3,ie)

	if ( isphe .eq. 1 ) then		!spherical
	  call ev_make_center(ie,dlon0,dlat0)
  	  xlon1=xgv(kn1)
	  ylat1=ygv(kn1)
	  xlon2=xgv(kn2)
	  ylat2=ygv(kn2)
	  xlon3=xgv(kn3)
	  ylat3=ygv(kn3)
          call ev_g2c(x1,y1,xlon1,ylat1,dlon0,dlat0)
          call ev_g2c(x2,y2,xlon2,ylat2,dlon0,dlat0)
          call ev_g2c(x3,y3,xlon3,ylat3,dlon0,dlat0)
        else					!cartesian
  	  x1=xgv(kn1)
	  y1=ygv(kn1)
	  x2=xgv(kn2)
	  y2=ygv(kn2)
	  x3=xgv(kn3)
	  y3=ygv(kn3)
	end if

	if( bdebug ) then
	  write(6,*) '============ ev debug =================='
	  write(6,*) ie,isphe
	  write(6,*) x1,y1
	  write(6,*) x2,y2
	  write(6,*) x3,y3
	end if

	a1=x2*y3-x3*y2
	a2=x3*y1-x1*y3
	a3=x1*y2-x2*y1
	!aj=a1+a2+a3
	aj = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)		!bug_f_64bit
	aji=one/aj
	b1=(y2-y3)*aji
	c1=(x3-x2)*aji
	b2=(y3-y1)*aji
	c2=(x1-x3)*aji
	b3=(y1-y2)*aji
	c3=(x2-x1)*aji
	a1=a1*aji
	a2=a2*aji
	a3=a3*aji
	!aj=aj/twofour		!bug 5.2.2010

	if( bdebug ) then
	  write(6,*) aji
	  write(6,*) b1,b2,b3
	  write(6,*) c1,c2,c3
	end if

c natural coordinates in triangle:   xi(i) = a(i) + b(i)*x + c(i)*y    i=1,3

	call adjust_bc(b1,b2,b3)
	call adjust_bc(c1,c2,c3)

	if( bdebug ) then
	  write(6,*) aji
	  write(6,*) b1,b2,b3
	  write(6,*) c1,c2,c3
	  write(6,*) '============ ev debug end =================='
	end if

	if(aj.le.0.) goto 99

	ss1=(x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)
	s1=sqrt(ss1)
	ss2=(x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)
	s2=sqrt(ss2)
	ss3=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)
	s3=sqrt(ss3)

	d1=acos((ss2+ss3-ss1)/(two*s2*s3))*rad
	d2=acos((ss1+ss3-ss2)/(two*s1*s3))*rad
	d3=acos((ss1+ss2-ss3)/(two*s1*s2))*rad

	dd1 = aj * sqrt( b1*b1 + c1*c1 )
	dd2 = aj * sqrt( b2*b2 + c2*c2 )
	dd3 = aj * sqrt( b3*b3 + c3*c3 )

	ev(1,ie)=a1		!a values for interpolation
	ev(2,ie)=a2
	ev(3,ie)=a3
	ev(4,ie)=b1		!b values (gradient in x)
	ev(5,ie)=b2
	ev(6,ie)=b3
	ev(7,ie)=c1		!c values (gradient in y)
	ev(8,ie)=c2
	ev(9,ie)=c3
	ev(10,ie)=aj/twofour	!aera = 12 * ev(10,ie)
	ev(11,ie)=d1		!angle on vertex
	ev(12,ie)=d2
	ev(13,ie)=d3
	ev(14,ie)=dd1		!for horizontal diffusion (?)
	ev(15,ie)=dd2
	ev(16,ie)=dd3
	ev(17,ie)=s1		!distance between vertices
	ev(18,ie)=s2
	ev(19,ie)=s3

        !write(96,*) ie,(ev(i,ie),i=1,evdim)

	end do

	!write(68,*) 'maxmax: ',maxmax

	return
   97	continue
        write(6,*) 'no basin has been read'
	stop 'error stop set_ev: no basin'
   99	continue
        write(6,*) 'set_ev : nodes not in anticlockwise sense'
        write(6,*) 'elem = ',ie,'  area = ',aj/2.,'  aj = ',aj
        write(6,*) 'nodes  x  y'
        write(6,*) kn1,x1,y1
        write(6,*) kn2,x2,y2
        write(6,*) kn3,x3,y3
	stop 'error stop set_ev: clockwise'
	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine is_init_ev(binit)

c checks if ev module has been initialized

	use evgeom

	implicit none

	logical binit

	binit = init_ev

	end

c****************************************************************

	subroutine set_coords_ev(isphe)

c sets type of coordinates to use with ev module

	use evgeom

	implicit none

	integer isphe

	isphe_ev = isphe

	end

c****************************************************************

	subroutine get_coords_ev(isphe)

c gets type of coordinates that is used with ev module

	use evgeom

	implicit none

	integer isphe

	isphe = isphe_ev

	end

c****************************************************************

	function is_spherical()

c checks if coordinates are spherical (lat/lon)

	use evgeom

	implicit none

	logical is_spherical

	is_spherical = isphe_ev .ne. 0

	end

c****************************************************************

	subroutine check_spheric_ev

c checks if coordinates are lat/lon

	use basin
	use evgeom

	implicit none

	include 'param.h'

	logical bverbose
	integer k,isphe
	real xmin,xmax,ymin,ymax

	bverbose = .true.
	bverbose = .false.

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
	if( ymin .lt.  -90. .or. ymax .gt.  90. ) isphe = 0

	if( isphe_ev .eq. -1 ) then	!determine automatically
	  if( isphe .eq. 1 ) then
	    if( bverbose ) then
	      write(6,*) 'Unsure about type of coordinates.'
	      write(6,*) 'Coodinates seem as lat/lon'
	      write(6,*) 'but are not flagged as such.'
	      write(6,*) 'Using lat/lon coordinates.'
	      write(6,*) 'If this is an error, then please set'
	      write(6,*) 'parameter isphe to the desired value.'
	    end if
	  end if
	else if( isphe_ev .ne. isphe ) then	!not the same
	  if( isphe .eq. 0 ) then	!not possible -> coords out of range
	    write(6,*) 'coodinates are flaged as lat/lon'
	    write(6,*) 'but are out of range.'
	    write(6,*) 'Please set parameter isphe = 0'
	    write(6,*) 'or adjust the grid'
	    write(6,*) 'xmin,xmax: ',xmin,xmax
	    write(6,*) 'ymin,ymax: ',ymin,ymax
	    stop 'error stop check_spheric_ev: not lat/lon coordinates'
	  end if
	  isphe = isphe_ev	!take desired value
	end if

	isphe_ev = isphe
	init_ev = .true.

	if( isphe == 1 ) then
	  write(6,*) 'using lat/lon coordinates'
	else
	  write(6,*) 'using cartesian coordinates'
	end if

	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine check_ev

c tests if ev is set up correctly

	use basin, only : nkn,nel,ngr,mbw
	use evgeom

	implicit none


	integer ie,ii,i,ip
	real bmax,cmax,bs,cs,b,c
	real alpha,atria

	real eps
	save eps
	data eps /1.e-5/

	write(6,*) 'testing setup of ev...'

	atria = 180.

	do ie=1,nel

	  bmax=0.
	  cmax=0.
	  bs=0.
	  cs=0.
	  alpha = 0.

	  do ii=1,3
		b = ev(3+ii,ie)
		c = ev(6+ii,ie)
		bs = bs + b
		cs = cs + c
		if(abs(b).gt.bmax) bmax=abs(b) 
		if(abs(c).gt.cmax) cmax=abs(c) 
		alpha = alpha + ev(10+ii,ie)
	  end do

	  ip=1
	  if(ev(10,ie) .le. 0.) goto 99
	  ip=2
	  if(bmax.le.0. .or. abs(bs)/bmax.gt.eps) goto 99
	  if(cmax.le.0. .or. abs(cs)/cmax.gt.eps) goto 99
	  ip=3
	  if( abs(alpha-atria)/180. .gt. eps ) goto 99

	end do

	write(6,*) 'test of ev passed...'

	return
   99	continue
	write(6,*) 'error: ',ip
	write(6,*) 'ie,area,alpha: ',ie,ev(10,ie),alpha
	write(6,*) 'bmax,cmax,bs,cs: ',bmax,cmax,bs,cs
	write(6,*) 'ev: ',(ev(i,ie),i=1,evdim)
	stop 'error stop check_ev: errors in array ev'
	end

c***********************************************************
c***********************************************************
c***********************************************************
c internal coordinates
c***********************************************************
c***********************************************************
c***********************************************************

	subroutine xy2xi(ie,x,y,xi)

c given x/y returns internal coordinates xi

	implicit none

	integer ie
	double precision x,y
	double precision xi(3)

	integer ii
	double precision a(3),b(3),c(3)

	call xi_abc(ie,a,b,c)

	do ii=1,3
	  xi(ii) = a(ii) + b(ii)*x + c(ii)*y
	end do

	call adjust_xi(xi)

	end
	
c***********************************************************

	subroutine adjust_xi(xi)

c adjusts internal coodinates xi

	use evgeom

	implicit none

	double precision xi(3)

	logical bdebug
	integer ii,it,is
	double precision xisum,xiso,xiadj
	double precision xiorig(3)

	bdebug = bdebug_internal

	xisum = 0.
	xiso = 0.
	is = 0
	it = 0

	do ii=1,3
	  xiorig(ii) = xi(ii)
	  xiso = xiso + xi(ii)
	  if( xi(ii) < 0. ) xi(ii) = 0.
	  if( xi(ii) > 1. ) xi(ii) = 1.
	  if( xi(ii) == 0. ) then
	    is = is + ii
	    it = it + 1
	  end if
	  xisum = xisum + xi(ii)
	end do

	if( it == 3 .or. abs(xiso-1.) > 1.e-8 ) then
	  write(6,*) 'xi is wrong... cannot adjust'
	  write(6,*) xiorig
	  stop 'error stop adjust_xi'
	end if

	xiadj = (xisum-1.)/(3-it)

	if( bdebug ) then
	  write(6,*) it,is,xisum,xiadj
	  write(6,*) xi
	end if

	if( it == 2 ) then
	  is = 6 - is
	  xi(is) = 1.
	else			!both for it == 1 and 2
	  do ii=1,3
	    if( ii /= is ) then
	      xi(ii) = xi(ii) - xiadj
	      if( xi(ii) < 0. ) xi(ii) = 0.
	      if( xi(ii) > 1. ) xi(ii) = 1.
	    end if
	  end do
	end if

	xiso = 0.
	do ii=1,3
	  xiso = xiso + xi(ii)
	end do

	if( bdebug ) then
	  write(6,*) it,is,xiso,xiso-1.
	  write(6,*) xi
	end if

	if( abs(xiso-1.) > 1.e-8 ) then
	  write(6,*) 'xi is still wrong... cannot adjust'
	  write(6,*) xiorig
	  stop 'error stop adjust_xi'
	end if

	!write(6,*) 'xi_adjust start...'
	!write(6,*) xiorig
	!write(6,*) xi
	!write(6,*) 'xi_adjust end...'

	end

c***********************************************************

	subroutine xi2xy(ie,x,y,xi)

c given internal coordinates xi returns x/y

	implicit none

	integer ie
	double precision x,y
	double precision xi(3)

	double precision eps
	parameter (eps = 1.e-8)

	double precision a1,a2,det
	double precision a(3),b(3),c(3)

	call xi_abc(ie,a,b,c)

	a1 = (xi(1)-a(1))
	a2 = (xi(2)-a(2))

	det = b(2)*c(1) - b(1)*c(2)
	if( abs(det) < eps ) goto 99

	y = (a1*b(2) - a2*b(1)) / det

	det = b(1)*c(2) - b(2)*c(1)
	if( abs(det) < eps ) goto 99

	x = (a1*c(2) - a2*c(1)) / det

	return
   99	continue
	write(6,*) 'det too small: ',ie,det
	write(6,*) x,y
	write(6,*) xi
	write(6,*) a
	write(6,*) b
	write(6,*) c
	stop 'error stop xi2xy: det == 0'
	end
	
c***********************************************************

	subroutine xi_dist(ie,xia,xib,l)

c computes distance between two points given in internal coordinates
c
c the two points must be lying on the sides of the triangle

	use evgeom

	implicit none

	integer ie			!element number
	double precision xia(3)		!internal coordinates of point a
	double precision xib(3)		!internal coordinates of point b
	double precision l		!length of segment (return)

	double precision pi,rad
	parameter (pi=4.D+0*atan(1.D+0),rad=pi/180.D+0)

	integer na0,nb0,ia0,ib0
	integer ia1,ib1,ia2,ib2,ia3,ib3
	integer il,il1,il2,il3
	double precision dla,dlb,cc

	integer ilbase,icbase,ii
	parameter (ilbase=16,icbase=10)
	double precision length,cosine
	length(ii) = ev(ilbase+ii,ie)
	cosine(ii) = cos(rad*ev(icbase+ii,ie))

	ia0 = 0
	ib0 = 0
	na0 = 0
	nb0 = 0
	do ii=1,3
	  if( xia(ii) == 0. ) then
	    na0 = na0 + 1
	    ia0 = ia0 + ii
	  end if
	  if( xib(ii) == 0. ) then
	    nb0 = nb0 + 1
	    ib0 = ib0 + ii
	  end if
	end do

	!write(6,*) ia0,ib0,na0,nb0

	if( na0 == 0 .or. nb0 == 0 ) then	!points must lie on side
	  write(6,*) 'points not on side of triangle...'
	  goto 99
	end if

	if( na0 == 2 .and. nb0 == 2 ) then	!line through two vertices
	  ia1 = 6 - ia0
	  ib1 = 6 - ib0
	  dla = 0.
	  cc = 0.
	  if( ia0 == ib0 ) then
	    dlb = 0.
	  else
	    il = ia0 + ib0 - 6
	    dlb = length(il)
	  end if
	else if( na0 == 1 ) then		!point a on real side ia1
	  ia1 = ia0
	  ia2 = mod(ia1,3) + 1
	  ia3 = mod(ia2,3) + 1
	  if( xib(ia1) == 0. ) then		!a and b on same side
	    dla = abs( xia(ia2)-xib(ia2) ) * length(ia1)
	    dlb = 0.
	    cc = 0.
	  else if( xib(ia2) == 0. ) then
	    dla = xia(ia2) * length(ia1)
	    dlb = xib(ia1) * length(ia2)
	    cc = cosine(ia3)
	  else
	    dla = xia(ia3) * length(ia1)
	    dlb = xib(ia1) * length(ia3)
	    cc = cosine(ia2)
	  end if
	else if( nb0 == 1 ) then		!point b on real side ib1
	  ib1 = ib0
	  ib2 = mod(ib1,3) + 1
	  ib3 = mod(ib2,3) + 1
	  if( xia(ib1) == 0. ) then		!a and b on same side
	    dla = abs( xia(ib2)-xib(ib2) ) * length(ib1)
	    dlb = 0.
	    cc = 0.
	  else if( xia(ib2) == 0. ) then
	    dlb = xib(ib2) * length(ib1)
	    dla = xia(ib1) * length(ib2)
	    cc = cosine(ib3)
	  else
	    dlb = xib(ib3) * length(ib1)
	    dla = xia(ib1) * length(ib3)
	    cc = cosine(ib2)
	  end if
	else
	  write(6,*) 'impossible... internal error'
	  goto 99
	end if

	l = sqrt( dla**2 + dlb**2 - 2.*dla*dlb*cc )

	return
   99	continue
	write(6,*) 'error in computing distance:'
	write(6,*) ia0,ib0,na0,nb0
	write(6,*) xia
	write(6,*) xib
	stop 'error stop xi_dist: generic error'
	end

c***********************************************************

	subroutine xi_abc(ie,a,b,c)

c returns a,b,c to compute xi
c
c natural coordinates in triangle:   xi(i) = a(i) + b(i)*x + c(i)*y    i=1,3

	use basin
	use evgeom

	implicit none

	integer ie

	double precision a(3),b(3),c(3)

	include 'param.h'

	integer ii
	integer kn1,kn2,kn3
	double precision x1,y1,x2,y2,x3,y3
	double precision a1,a2,a3,aj,aji

	!write(6,*) 'ggu_xi ',ie,isphe_ev

	if( isphe_ev .eq. 0 ) then	!cartesian
	  do ii=1,3
	    a(ii) = ev(ii,ie)
	    b(ii) = ev(3+ii,ie)
	    c(ii) = ev(6+ii,ie)
	  end do
	  return
	end if

	kn1=nen3v(1,ie)
	kn2=nen3v(2,ie)
	kn3=nen3v(3,ie)

  	x1=xgv(kn1)
	y1=ygv(kn1)
	x2=xgv(kn2)
	y2=ygv(kn2)
	x3=xgv(kn3)
	y3=ygv(kn3)

	a1=x2*y3-x3*y2
	a2=x3*y1-x1*y3
	a3=x1*y2-x2*y1
	aj = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
	aji=1./aj

	a(1)=a1*aji
	a(2)=a2*aji
	a(3)=a3*aji
	b(1)=(y2-y3)*aji
	c(1)=(x3-x2)*aji
	b(2)=(y3-y1)*aji
	c(2)=(x1-x3)*aji
	b(3)=(y1-y2)*aji
	c(3)=(x2-x1)*aji

	end

c***********************************************************
c***********************************************************
c***********************************************************

        subroutine cpp(x,y,rlambda,phi,rlambda0,phi0)

c transforms (lon,lat) into cartesian coordinates (x,y) (lon,lat in radians)

        implicit none

        double precision x,y		!cartesian x,y [m]
        double precision rlambda,phi	!lon,lat [rad]
        double precision rlambda0,phi0	!center of projection [rad]

        double precision r		!earth radius [m]
	parameter ( r = 6378206.4D0 )

        x = r*(rlambda - rlambda0)*dcos(phi0)
        y = (phi-phi0)*r
        !y = phi*r

        end

c***********************************************************

	subroutine ev_make_center(ie,xm,ym)

	use basin

	implicit none

	integer ie
	double precision xm,ym

	include 'param.h'

	integer ii,k

	xm = 0.
	ym = 0.
	do ii=1,3
	  k = nen3v(ii,ie)
	  xm = xm + xgv(k)
	  ym = ym + ygv(k)
	end do
	xm = xm / 3.
	ym = ym / 3.

	end

c***********************************************************

        subroutine ev_g2c(x,y,lambda,phi,lambda0,phi0)

c transforms geographical (lon,lat) into cartesian (x,y) coordinates 
c
c (lon,lat in degrees)

        implicit none

        double precision x,y		!cartesian x,y [m]
        double precision lambda,phi	!lon,lat [deg]
        double precision lambda0,phi0	!center of projection [deg]

        double precision r		!earth radius [m]
	parameter ( r = 6378206.4D0 )
        double precision pi		!pi
	parameter ( pi = 3.14159265358979D0 )
        double precision rad		!convert degrees to radian
	parameter ( rad = pi / 180.0 )

	double precision dlambda,dphi

	dlambda = rad * (lambda - lambda0)
	dphi    = rad * (phi - phi0)
	!dphi    = rad * phi

        x = r*dlambda*dcos(rad*phi0)
        y = r*dphi

        end

c***********************************************************

        subroutine ev_c2g(x,y,lambda,phi,lambda0,phi0)

c transforms cartesian (x,y) into geographical (lon,lat) coordinates 
c
c (lon,lat in degrees)

        implicit none

        double precision x,y		!cartesian x,y [m]
        double precision lambda,phi	!lon,lat [deg]
        double precision lambda0,phi0	!center of projection [deg]

        double precision r		!earth radius [m]
	parameter ( r = 6378206.4D0 )
        double precision pi		!pi
	parameter ( pi = 3.14159265358979D0 )
        double precision rad,rrad	!convert degrees to radian
	parameter ( rad = pi / 180.0 , rrad = 1. / rad )

	lambda = lambda0 + rrad * x / (r*dcos(rad*phi0))
	phi = phi0 + rrad * y / r
	!phi = rrad * y / r

        end

c***********************************************************
c***********************************************************
c***********************************************************

	subroutine adjust_bc(v1,v2,v3)

c adjusts b/c so that sum = 0

	implicit none

	double precision v1,v2,v3

	double precision vv

	vv = v1 + v2 + v3

	if( vv .eq. 0. ) return

	if( v1 .eq. 0. .and. v2 .eq. 0. ) goto 99
	if( v1 .eq. 0. .and. v3 .eq. 0. ) goto 99
	if( v2 .eq. 0. .and. v3 .eq. 0. ) goto 99

	if( v1 .eq. 0. ) then
	  v2 = v2 - vv/2.
	  v3 = -v2
	else if( v2 .eq. 0. ) then
	  v3 = v3 - vv/2.
	  v1 = -v3			!BUGFIX
	else if( v3 .eq. 0. ) then
	  v1 = v1 - vv/2.
	  v2 = -v1			!BUGFIX
	else
	  v1 = v1 - vv/3.
	  v2 = v2 - vv/3.
	  v3 = v3 - vv/3.
	  v3 = - v1 - v2
	end if

	return
   99	continue
	write(6,*) v1,v2,v3
	stop 'error stop adjust_bc: internal error'
	end

c***********************************************************
c***********************************************************
c***********************************************************

	function area_elem(ie)

c returns area of element ie

	use evgeom

	implicit none

	real area_elem
	integer ie

	if( .not. init_ev ) then
	  stop 'error stop area_elem: ev not set up'
	end if

	area_elem = 12. * ev(10,ie)

	end

c***********************************************************

	function aomega_elem(ie)

c returns aomega of element ie

	use evgeom

	implicit none

	real aomega_elem
	integer ie

	if( .not. init_ev ) then
	  stop 'error stop aomega_elem: ev not set up'
	end if

	aomega_elem = ev(10,ie)

	end

c***********************************************************

	function weight_elem(ie)

c returns weight for element ie - if ev is not setup, return 1

	use evgeom

	implicit none

	real weight_elem
	integer ie

	if( init_ev ) then
	  weight_elem = ev(10,ie)
	else
	  weight_elem = 1.
	end if

	end

c***********************************************************
c***********************************************************
c***********************************************************

	subroutine gradient_normal(ie,ii,gx,gy,gxn,gyn)

c computes gradient normal to given direction
c
c direction is given by the two nodes not identified by ii

	use evgeom

	implicit none

	integer ie	!element number [1-nel]
	integer ii	!node that is not part of the direction [1-3]
	real gx,gy	!original vector (gradient) 
	real gxn,gyn	!projected vector (gradient) (return)

	double precision bii,cii,gn

	bii = ev(3+ii,ie)
	cii = ev(6+ii,ie)

	gn = (gx*bii + gy*cii) / (bii*bii + cii*cii)

	gxn = gn * bii
	gyn = gn * cii

	end

c***********************************************************

	subroutine gradient_tangent(ie,ii,gx,gy,gxt,gyt)

c computes gradient tangent to given direction
c
c direction is given by the two nodes not identified by ii

	use evgeom

	implicit none

	integer ie	!element number [1-nel]
	integer ii	!node that is not part of the direction [1-3]
	real gx,gy	!original vector (gradient) 
	real gxt,gyt	!projected vector (gradient) (return)

	double precision bii,cii,gn

	bii = ev(3+ii,ie)
	cii = ev(6+ii,ie)

	gn = (gx*cii - gy*bii) / (bii*bii + cii*cii)

	gxt =  gn * cii
	gyt = -gn * bii

	end

c***********************************************************

	subroutine gradient_normal_d(ie,ii,gx,gy,gxn,gyn)

c computes gradient normal to given direction
c
c direction is given by the two nodes not identified by ii

	use evgeom

	implicit none

	integer ie		 !element number [1-nel]
	integer ii		 !node that is not part of the direction [1-3]
	double precision gx,gy	 !original vector (gradient) 
	double precision gxn,gyn !projected vector (gradient) (return)

	double precision bii,cii,gn

	bii = ev(3+ii,ie)
	cii = ev(6+ii,ie)

	gn = (gx*bii + gy*cii) / (bii*bii + cii*cii)

	gxn = gn * bii
	gyn = gn * cii

	end

c***********************************************************

	subroutine gradient_tangent_d(ie,ii,gx,gy,gxt,gyt)

c computes gradient tangent to given direction
c
c direction is given by the two nodes not identified by ii

	use evgeom

	implicit none

	integer ie		 !element number [1-nel]
	integer ii		 !node that is not part of the direction [1-3]
	double precision gx,gy	 !original vector (gradient) 
	double precision gxt,gyt !projected vector (gradient) (return)

	double precision bii,cii,gn

	bii = ev(3+ii,ie)
	cii = ev(6+ii,ie)

	gn = (gx*cii - gy*bii) / (bii*bii + cii*cii)

	gxt =  gn * cii
	gyt = -gn * bii

	end

c***********************************************************

	subroutine compute_distance(xin1,yin1,xin2,yin2,dx,dy)

c computes distance between (x1,y1) and (x2,y2)
c
c takes care of lat/lon coordinates

	use evgeom

	implicit none

	real xin1,yin1
	real xin2,yin2
	real dx,dy

	double precision xl1,yl1,xl2,yl2
	double precision x1,y1,x2,y2
	double precision dlon0,dlat0

	if ( isphe_ev .eq. 1 ) then		!spherical
  	  xl1 = xin1
	  yl1 = yin1
	  xl2 = xin2
	  yl2 = yin2
	  dlon0 = 0.5*(xl1+xl2)
	  dlat0 = 0.5*(yl1+yl2)
          call ev_g2c(x1,y1,xl1,yl1,dlon0,dlat0)
          call ev_g2c(x2,y2,xl2,yl2,dlon0,dlat0)
        else					!cartesian
  	  x1 = xin1
	  y1 = yin1
	  x2 = xin2
	  y2 = yin2
	end if

	dx = x2 - x1
	dy = y2 - y1

	end

c***********************************************************

        subroutine compute_cartesian_coords(n,lon,lat,x,y)

	use evgeom

	implicit none

	integer n		!total number of points to convert
	real lon(n),lat(n)	!coordinates of points to convert
	real x(n),y(n)		!converted cartesian coordinates (return)

	integer i
	double precision lambda0,phi0
	double precision lambda,phi
	double precision xx,yy

c---------------------------------------------------------
c see if spherical coordinates
c---------------------------------------------------------

	if ( isphe_ev .eq. 0 ) then		!already cartesian
	  do i=1,n
	    x(i) = lon(i)
	    y(i) = lat(i)
	  end do
	  return
	end if

c---------------------------------------------------------
c compute reference point
c---------------------------------------------------------

	lambda0 = 0.
	phi0 = 0.
	do i=1,n
	  lambda0 = lambda0 + lon(i)
	  phi0 = phi0 + lat(i)
	end do
	lambda0 = lambda0 / n
	phi0 = phi0 / n
	
c---------------------------------------------------------
c do conversion
c---------------------------------------------------------

	do i=1,n
	  lambda = lon(i)
	  phi = lat(i)
          call ev_g2c(xx,yy,lambda,phi,lambda0,phi0)
	  x(i) = xx
	  y(i) = yy
	end do

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c***********************************************************

	subroutine ev_set_debug(bdebug)

	use evgeom

	logical bdebug

	bdebug_internal = bdebug

	end

c***********************************************************

