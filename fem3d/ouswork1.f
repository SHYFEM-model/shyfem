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

	program ouswork1

c reads ous file and elaborates it for altimeter trace (version 1)

	use mod_depth
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

        include 'param.h'

	include 'simul.h'




	real znv1(nkndim)
	real znv2(nkndim)
	real znv(nkndim)
	real zenv(3,neldim)

        integer nvers,nin,nlv
        integer itanf,itend,idt,idtous
	integer it_record
	integer it,ie,i
	integer it1,it2,itdata,id
        integer ierr,nread,ndry
        integer nknous,nelous,nlvous
	integer itmin,itmax
	logical btwin,btsim
        real href,hzoff,hlvmin
	real volume
	real zmin,zmax
	real umin,umax
	real vmin,vmax
	real x,y,z
	character*40 name

c	integer rdous,rfous
	integer iapini,ideffi

c-----------------------------------------------------------------
c initialize basin and simulation
c-----------------------------------------------------------------

	itmin = 0
	itmax = 86400		!set to 0 for no time window
	btwin = itmin .lt. itmax
	btsim = .true.		!error if time record not found in sim

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

c-----------------------------------------------------------------
c initialize date and and simulation
c-----------------------------------------------------------------

	name='data_in'
	open(11,file=name,status='old',form='formatted')
	name='data_out'
	open(12,file=name,status='unknown',form='formatted')

	call init_date		!reads file date0 and initializes date

        call rdous(nin,it1,nlvdim,ilhv,znv1,zenv,utlnv,vtlnv,ierr)
	if( ierr .ne. 0 ) goto 100
        call rdous(nin,it2,nlvdim,ilhv,znv2,zenv,utlnv,vtlnv,ierr)
	if( ierr .ne. 0 ) goto 100

	itdata = it1 - 1
	do while( itdata .lt. it1 )	!read until time in sim window
	  call read_data(id,itdata,x,y,ierr)
	  if( ierr .ne. 0 ) goto 100
	  if( btsim .and. itdata .lt. it1 ) goto 99  !no time record for data
	  if( btwin .and. itdata .lt. itmin ) goto 98		
	  if( btwin .and. itdata .gt. itmax ) goto 98	
	end do

c-----------------------------------------------------------------
c loop on data of simulation
c-----------------------------------------------------------------

  300   continue

	if( itdata .gt. it2 ) then
	  call copy_znv(nkn,it1,it2,znv1,znv2)
          call rdous(nin,it2,nlvdim,ilhv,znv2,zenv,utlnv,vtlnv,ierr)
	  if( btsim .and. ierr .ne. 0 ) goto 99	!no time records for data
	  if( ierr .ne. 0 ) goto 100		!no time records for data
	else
	  call intp_znv(nkn,it1,it2,znv1,znv2,itdata,znv)
	  call intp_data(nkn,znv,x,y,z)
	  call write_data(id,itdata,x,y,z)
	  call read_data(id,itdata,x,y,ierr)
	  if( ierr .ne. 0 ) goto 100		!no more data to handle
	  if( btwin .and. itdata .lt. 0 ) goto 98		
	  if( btwin .and. itdata .gt. 86400 ) goto 98	
	end if

	goto 300

c-----------------------------------------------------------------
c end of loop
c-----------------------------------------------------------------

  100	continue

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	stop
   98	continue
	write(6,*) 'it2,it2,itdata: ',it1,it2,itdata
	stop 'error stop: out of desired time window'
   99	continue
	write(6,*) 'it2,it2,itdata: ',it1,it2,itdata
	stop 'error stop: out of simulation time window'
	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine intp_data(nkn,znv,x,y,z)

	implicit none

	integer nkn
	real znv(nkn)
	real x,y,z

	logical binside
	real zeta
	real z3(3)

	integer ie
	save ie
	data ie /0/

	call find_elem_from_old(ie,x,y,ie)
	call find_adriatico(x,y,binside)

	zeta = -999.
	if( ie .gt. 0 .and. binside ) then
	  call get_scal_elem(ie,znv,z3)
	  call femintp(ie,z3,x,y,zeta)
	end if

	z = zeta

	end

c******************************************************************

	subroutine copy_znv(nkn,it1,it2,znv1,znv2)

	implicit none

	integer nkn,it1,it2
	real znv1(nkn),znv2(nkn)

	integer k

	it1 = it2
	do k=1,nkn
	  znv1(k) = znv2(k)
	end do

	end

c******************************************************************

	subroutine intp_znv(nkn,it1,it2,znv1,znv2,itdata,znv)

	implicit none

	integer nkn,it1,it2
	integer itdata
	real znv1(nkn),znv2(nkn)
	real znv(nkn)

	integer k
	real r

	r = float(itdata-it1)/float(it2-it1)

	do k=1,nkn
	  znv(k) = znv1(k) + r * ( znv2(k) - znv1(k) )
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

	subroutine init_date

	implicit none

	character*60 name
	integer date,time

	name='date0'

	open(1,file=name,status='old',form='formatted')
	read(1,*) date
	close(1)

	write(6,*) 'init_date: ',date

	time = 0
	call dtsini(date,time)

	end

c******************************************************************

	subroutine read_data(id,it,x,y,ierr)

c reads data file, converts date and time and returns itdate,x,y

	implicit none

	integer id,it
	real x,y
	integer ierr

	integer date,time
	integer year,month,day
	integer hour,min,sec
	real z

    9	continue
	read(11,*,end=1) id,date,time,x,y

	if( date .eq. -999 .or. time .eq. -999 ) then
	  z = -999.
	  write(12,'(3i10,3f12.5)') id,date,time,x,y,z
	  goto 9
	end if

	call unpackdate(date,year,month,day)
	call unpacktime(time,hour,min,sec)
	call dts2it(it,year,month,day,hour,min,sec)
	ierr = 0

	return
    1	continue
	ierr = -1

	end

c******************************************************************

	subroutine write_data(id,it,x,y,z)

c reads data file, converts date and time and returns itdate,x,y

	implicit none

	integer id,it
	real x,y,z

	integer date,time
	integer year,month,day
	integer hour,min,sec

	call dts2dt(it,year,month,day,hour,min,sec)
	call packdate(date,year,month,day)
	call packtime(time,hour,min,sec)

	write(12,'(3i10,3f12.5)') id,date,time,x,y,z

	end

c******************************************************************
c******************************************************************
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

