c
c $Id: ouswork.f,v 1.15 2009-11-18 16:50:37 georg Exp $
c
c info on OUS files
c
c revision log :
c
c 02.09.2003	ggu	adapted to new OUS format
c 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
c 16.10.2007	ggu	new debug routine
c 27.10.2009    ggu     include evmain.h, compute volume
c 23.03.2011    ggu     compute real u/v-min/max of first level
c 16.12.2011    ggu     bug fix: call to init_sigma_info and makehev (common hev)
c
c***************************************************************

	program ouswork

c reads ous file and elaborates it for altimeter trace (version 0)

	implicit none

        include 'param.h'
	include 'evmain.h'

	character*80 descrr,descrp
	common /descrr/ descrr
	common /descrp/ descrp
	include 'nbasin.h'

	real xgv(nkndim), ygv(nkndim)
	real hm3v(3,neldim)
	integer nen3v(3,neldim)
	integer ipev(neldim), ipv(nkndim)
	integer iarv(neldim)
	common /xgv/xgv, /ygv/ygv
	common /hm3v/hm3v
	common /nen3v/nen3v
	common /ipev/ipev, /ipv/ipv
	common /iarv/iarv

	integer ilhv(neldim)
	real hlv(nlvdim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)
	common /ilhv/ilhv
	common /hlv/hlv
        common /utlnv/utlnv
        common /vtlnv/vtlnv

	real hev(neldim)
	common /hev/hev

	real znv(nkndim)
	real zenv(3,neldim)

        integer nvers,nin,nlv
        integer itanf,itend,idt,idtous
	integer it_record
	integer it,ie,i,ks
        integer ierr,nread,ndry
        integer nknous,nelous,nlvous
        real href,hzoff,hlvmin
	real volume
	real zmin,zmax
	real umin,umax
	real vmin,vmax

c	integer rdous,rfous
	integer iapini,ideffi

c-----------------------------------------------------------------
c initialize basin and simulation
c-----------------------------------------------------------------

	nread=0
	ks = 0

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	call set_ev

	nin=ideffi('datdir','runnam','.ous','unform','old')
	if(nin.le.0) goto 100

c-----------------------------------------------------------------
c read header of simulation
c-----------------------------------------------------------------

	nvers=1
        call rfous(nin
     +			,nvers
     +			,nknous,nelous,nlvous
     +			,href,hzoff
     +			,descrp
     +			,ierr)

	nlv=nlvous
	call dimous(nin,nkndim,neldim,nlvdim)

        write(6,*)
        write(6,*)   descrp
        write(6,*)
        write(6,*) ' nvers        : ',nvers
        write(6,*) ' href,hzoff   : ',href,hzoff
        write(6,*) ' nkn,nel      : ',nknous,nelous
        write(6,*) ' nlv          : ',nlvous
        write(6,*)

	call rsous(nin,ilhv,hlv,hev,ierr)

        call init_sigma_info(nlv,hlv)
	call makehev(hev)

	call init_date(it_record)

c-----------------------------------------------------------------
c loop on data of simulation
c-----------------------------------------------------------------

  300   continue

        call rdous(nin,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)

        if(ierr.gt.0) then
		write(6,*) 'error in reading file : ',ierr
		goto 100
        else if(ierr.lt.0) then
		goto 100
	end if

	nread=nread+1

	call mima(znv,nknous,zmin,zmax)
        call comp_vel(1,nel,hev,zenv,nlvdim,utlnv,vtlnv
     +			,umin,vmin,umax,vmax)
	call compute_volume(nel,zenv,hev,volume)

c        call debug_write_node(ks,it,nread,nkndim,neldim,nlvdim,nkn,nel,nlv
c     +          ,nen3v,zenv,znv,utlnv,vtlnv)

	!write(6,*) 
	!write(6,*) 'time : ',it
	!call get_date(it)
	!write(6,*) 
	!write(6,*) 'zmin/zmax : ',zmin,zmax
	!write(6,*) 'umin/umax : ',umin,umax
	!write(6,*) 'vmin/vmax : ',vmin,vmax
	!write(6,*) 'volume    : ',volume

	if( it .eq. it_record ) then
	  call elab_z(it,nkn,znv)
	end if

	goto 300

  100	continue

c-----------------------------------------------------------------
c end of loop
c-----------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	stop
	end

c******************************************************************

        subroutine comp_vel(level,nel,hev,zenv,nlvdim,utlnv,vtlnv
     +			,umin,vmin,umax,vmax)

        implicit none

        integer level
        integer nel
        real hev(1)
        real zenv(3,1)
        integer nlvdim
        real utlnv(nlvdim,1)
        real vtlnv(nlvdim,1)
        real umin,vmin
        real umax,vmax

        integer ie,ii
        real zmed,hmed,u,v

	umin =  1.e+30
	vmin =  1.e+30
        umax = -1.e+30
        vmax = -1.e+30

        do ie=1,nel
          zmed = 0.
          do ii=1,3
            zmed = zmed + zenv(ii,ie)
          end do
          zmed = zmed / 3.
          hmed = hev(ie) + zmed
          if( hmed .le. 0. ) stop 'error stop hmed...'

          u = utlnv(level,ie) / hmed
          v = vtlnv(level,ie) / hmed

          umin = min(umin,u)
          vmin = min(vmin,v)
          umax = max(umax,u)
          vmax = max(vmax,v)
        end do

        end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine get_date(it)

	implicit none

	integer it
	integer year,month,day,hour,min,sec
	character*20 line

	call dts2dt(it,year,month,day,hour,min,sec)
	call dtsform(year,month,day,hour,min,sec,line)

	write(6,*) line

	end

c******************************************************************

	subroutine init_date(it_record)

	implicit none

	integer it_record
	character*60 name
	integer date,time

	name='date0'

	open(1,file=name,status='old',form='formatted')
	read(1,*) date,it_record
	close(1)

	write(6,*) 'init_date: ',date,it_record

	time = 0
	call dtsini(date,time)

	end

c******************************************************************

	subroutine elab_z(it,nkn,znv)

	implicit none

	integer it
	integer nkn
	real znv(nkn)

	character*60 name
	real x,y,aux,ssl
	real zeta
	real z(3)
	real dist,dm,dmin,dhole
	integer ie,i
	logical binside,berror,bhole

	real compute_dist

	!if( it .ne. 410400 ) return

	write(6,*)
	write(6,*) '****** elaborating data for it = ',it
	write(6,*)

	call get_min_distance(dmin)
	dhole = 1.5 * dmin
	write(6,*) 'min distance is ',dmin,dhole

	name='sev.txt'
	open(1,file=name,status='old',form='formatted')

	i = 0
	ie = 0
	dm = 0.
    1	continue
	  read(1,*,end=2) y,x,aux,ssl
	  call find_elem_from_old(ie,x,y,ie)
	  zeta = -999.
	  call find_adriatico(x,y,binside)
	  !if( ie .gt. 0 .and. ssl .lt. 10. .and. binside ) then
	  if( ie .gt. 0 .and. binside ) then
	    berror = ssl .ge. 10.
	    !if( berror ) ssl = 0.
	    if( berror ) goto 3
	    call get_scal_elem(ie,znv,z)
	    call femintp(ie,z,x,y,zeta)
	    dist = compute_dist(x,y,dm)
	    if( dm .gt. 0. ) dmin = min(dmin,dm)
	    bhole = dm .gt. dhole
	    i = i + 1
	    if( berror .or. bhole ) then
	      write(65,*) 
	      write(66,*) 
	    end if
	    if( .not. berror ) then
	      write(65,*) dist,ssl,zeta
	      write(66,*) i,ssl,zeta
	    end if
	    write(67,*) i,x,y,ssl,zeta
	    write(68,1000) 1,i,0,x,y,ssl
	  end if
    3	  continue
	  !write(6,*) x,y,ssl,zeta,dm
	goto 1
    2	continue

	close(1)

	write(6,*)
	write(6,*) 'exiting...  dmin = ',dmin
	write(6,*)
	stop
 1000	format(i1,2i10,3f14.8)
	end

c******************************************************************

	subroutine find_adriatico(x,y,binside)

c finds points in Adriatic Sea

	implicit none

	real x,y
	logical binside

	real x1,y1,x2,y2
	real dx,dy,dnx,dny
	real scal

	binside = .false.

	if( y .lt. 40. ) return

	x1 = 9.84
	y1 = 44.55
	x2 = 18.41
	y2 = 40.15

	dx = x2 - x1
	dy = y2 - y1

	dny = dx
	dnx = -dy

	scal = dnx*(x-x1) + dny*(y-y1)

	if( scal .lt. 0. ) return

	binside = .true.

	end

c******************************************************************

	function compute_dist(x,y,dm)

	implicit none

	real x,y
	real dm
	real compute_dist

	real phi_0,pi,rad
	real xc,yc

	real x0,y0,xold,yold
	save x0,y0,xold,yold
	real xfact,yfact,lon_0,lat_0
	save xfact,yfact,lon_0,lat_0

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. 0 ) then
          lon_0 = 16.
          lat_0 = 43.
          phi_0 = 43.

          pi = 4. * atan(1.d+0)
          rad = pi / 180.

          yfact = 60 * 1852               !distance of 1 min = 1 nautical mile
          xfact = yfact * cos(phi_0*rad)
	end if

        xc = (x-lon_0) * xfact
        yc = (y-lat_0) * yfact

	if( icall .eq. 0 ) then
	  x0 = xc
	  y0 = yc
	  xold = xc
	  yold = yc
	end if

	icall = icall + 1

	compute_dist = sqrt( (xc-x0)**2 + (yc-y0)**2 )
	dm = sqrt( (xc-xold)**2 + (yc-yold)**2 )

	xold = xc
	yold = yc

	end

c******************************************************************

	subroutine get_min_distance(dmin)

c computes minimum distance between points from track

	implicit none

	real dmin

	real dm,x,y,dist
	real compute_dist
	logical binside
	character*80 name

	name='sev.txt'
	open(1,file=name,status='old',form='formatted')

	dm = 0.
	dmin = 1.e+30
    1	continue
	  read(1,*,end=2) y,x
	  call find_adriatico(x,y,binside)
	  if( binside ) then
	    dist = compute_dist(x,y,dm)
	    if( dm .gt. 0. ) dmin = min(dmin,dm)
	  end if
	goto 1
    2	continue

	close(1)

	end

c******************************************************************

