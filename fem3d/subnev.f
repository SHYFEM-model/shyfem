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
c old routines
c
c subroutine sp110a				set up ev vector
c subroutine testev				tests if ev is set up
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
c
c***********************************************************

	subroutine set_ev

c set up ev vector
c
c double precision version
c
c revised on 31.08.88 by ggu (czv containes real chezy)
c revised on 25.11.88 by ggu (czv eliminated)
c revised on 28.01.92 by ggu (double precision, implicit none)

	implicit none

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	include 'ev.h'
	integer isphe_ev,init_ev
	common /evcommon/ isphe_ev,init_ev

	integer nen3v(3,1)
	common /nen3v/nen3v
	real xgv(1),ygv(1)
	common /xgv/xgv,/ygv/ygv

	integer ie,i,ii,kn1,kn2,kn3
	integer isphe
        double precision one,two,four,twofour,rad
	double precision x1,x2,x3,y1,y2,y3
	double precision a1,a2,a3,b1,b2,b3,c1,c2,c3
	double precision aj,aji,pi
	double precision s1,s2,s3,ss1,ss2,ss3
	double precision d1,d2,d3
	double precision dd1,dd2,dd3

	double precision xlon1,ylat1,xlon2,ylat2,xlon3,ylat3	!lat/long [rad]
	double precision dlat0,dlon0			!center of projection

	call check_spheric_ev	!checks and sets isphe_ev

	isphe = isphe_ev
	init_ev = 1

c        write(6,*)
c        do ie=1,10
c          write(6,*) ie,(nen3v(ii,ie),ii=1,3)
c        end do

        one = 1.
        two = 2.
        four = 4.
        twofour = 24.

	pi=four*atan(one)
        rad = 180./pi

	do i=1,nel

	kn1=nen3v(1,i)
	kn2=nen3v(2,i)
	kn3=nen3v(3,i)

	if ( isphe .eq. 1 ) then		!spherical
  	  xlon1=xgv(kn1)/rad
	  ylat1=ygv(kn1)/rad
	  xlon2=xgv(kn2)/rad
	  ylat2=ygv(kn2)/rad
	  xlon3=xgv(kn3)/rad
	  ylat3=ygv(kn3)/rad
	  dlon0 = (xlon1+xlon2+xlon3) / 3.
	  dlat0 = (ylat1+ylat2+ylat3) / 3.
          call cpp(x1,y1,xlon1,ylat1,dlon0,dlat0)
          call cpp(x2,y2,xlon2,ylat2,dlon0,dlat0)
          call cpp(x3,y3,xlon3,ylat3,dlon0,dlat0)
        else					!cartesian
  	  x1=xgv(kn1)
	  y1=ygv(kn1)
	  x2=xgv(kn2)
	  y2=ygv(kn2)
	  x3=xgv(kn3)
	  y3=ygv(kn3)
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

	call adjust_bc(b1,b2,b3)
	call adjust_bc(c1,c2,c3)

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

	ev(1,i)=a1		!a values for interpolation
	ev(2,i)=a2
	ev(3,i)=a3
	ev(4,i)=b1		!b values (gradient in x)
	ev(5,i)=b2
	ev(6,i)=b3
	ev(7,i)=c1		!c values (gradient in y)
	ev(8,i)=c2
	ev(9,i)=c3
	ev(10,i)=aj/twofour	!aera * 12
	ev(11,i)=d1		!angle on vertex
	ev(12,i)=d2
	ev(13,i)=d3
	ev(14,i)=dd1		!for horizontal diffusion (?)
	ev(15,i)=dd2
	ev(16,i)=dd3
	ev(17,i)=s1		!distance between vertices
	ev(18,i)=s2
	ev(19,i)=s3

        !write(96,*) i,(ev(ii,i),ii=1,evdim)

	end do

	return
   99	continue
        write(6,*) 'set_ev : nodes not in anticlockwise sense'
        write(6,*) 'elem = ',i,'  area = ',aj
        write(6,*) 'nodes  x  y'
        write(6,*) kn1,x1,y1
        write(6,*) kn2,x2,y2
        write(6,*) kn3,x3,y3
	stop 'error stop set_ev'
	end

c****************************************************************
c****************************************************************
c****************************************************************

	blockdata ev_blockdata

c initializes isphe_ev
c
c isphe_ev:
c
c -1	value not given -> try to determine automatically
c  0	cartesian
c  1	spherical (lat/lon)

	implicit none

	integer isphe_ev,init_ev
	common /evcommon/ isphe_ev,init_ev
	save /evcommon/
	data isphe_ev,init_ev / -1 , 0 /

	end

c****************************************************************

	subroutine is_init_ev(binit)

c checks if ev module has been initialized

	implicit none

	logical binit

	integer isphe_ev,init_ev
	common /evcommon/ isphe_ev,init_ev
	
	binit = init_ev .gt. 0

	end

c****************************************************************

	subroutine set_coords_ev(isphe)

c sets type of coordinates to use with ev module

	implicit none

	integer isphe

	integer isphe_ev,init_ev
	common /evcommon/ isphe_ev,init_ev
	
	isphe_ev = isphe

	end

c****************************************************************

	subroutine get_coords_ev(isphe)

c gets type of coordinates that is used with ev module

	implicit none

	integer isphe

	integer isphe_ev,init_ev
	common /evcommon/ isphe_ev,init_ev
	
	isphe = isphe_ev

	end

c****************************************************************

	subroutine check_spheric_ev

c checks if coordinates are lat/lon

	implicit none

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real xgv(1), ygv(1)
        common /xgv/xgv, /ygv/ygv

	integer isphe_ev,init_ev
	common /evcommon/ isphe_ev,init_ev

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
	if( ymin .lt.  -90. .or. ymax .gt.  90. ) isphe = 0

	if( isphe_ev .eq. -1 ) then	!determine automatically
	  if( isphe .eq. 1 ) then
	    write(6,*) 'Unsure about type of coordinates.'
	    write(6,*) 'Coodinates seem lat/lon but are not flagged as such.'
	    write(6,*) 'Using lat/lon coordinates.'
	    write(6,*) 'If this is an error, then please set'
	    write(6,*) 'parameter isphe to the desired value.'
	    !stop 'error stop check_spheric_ev: coord type'
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

	write(6,*) 'setting for coordinates: ',isphe
	if( isphe .ne. 0 ) write(6,*) 'using lat/lon coordinates'

	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine check_ev

c tests if ev is set up correctly

	implicit none

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	include 'ev.h'
c local
	integer ie,i,ip
	real bmax,cmax,bs,cs,b,c
	real alpha,atria
c save&data
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

	  do i=1,3
		b = ev(3+i,ie)
		c = ev(6+i,ie)
		bs = bs + b
		cs = cs + c
		if(abs(b).gt.bmax) bmax=abs(b) 
		if(abs(c).gt.cmax) cmax=abs(c) 
		alpha = alpha + ev(10+i,ie)
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

c****************************************************************

	subroutine testev

c old test

	implicit none

	call check_ev

	end

c***********************************************************

	subroutine sp110a

c old set

	implicit none

	call set_ev

	end

c***********************************************************

        subroutine cpp(x,y,rlambda,phi,rlambda0,phi0)

c transforms (lon,lat) into cartesian coordinates (x,y) (lon,lat in radians)

        implicit none

        double precision x,y		!cartesian x,y [m]
        double precision rlambda,phi	!lon,lat [rad]
        double precision rlambda0,phi0	!center of projection [rad]

        double precision r		!earth radius [m]
	parameter ( r = 6378206.4 )

        x = r*(rlambda - rlambda0)*dcos(phi0)
        y = phi*r

        end

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
	  v1 = -v1
	else if( v3 .eq. 0. ) then
	  v1 = v1 - vv/2.
	  v2 = -v2
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

	function area_elem(ie)

c returns area of element ie

	implicit none

	real area_elem
	integer ie

	include 'ev.h'

	area_elem = 12. * ev(10,ie)

	end

c***********************************************************

	function aomega_elem(ie)

c returns aomega of element ie

	implicit none

	real aomega_elem
	integer ie

	include 'ev.h'

	aomega_elem = ev(10,ie)

	end

c***********************************************************

