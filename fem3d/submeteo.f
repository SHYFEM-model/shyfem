c
c $Id: submeteo.f,v 1.7 2010-02-26 17:35:06 georg Exp $
c
c handle regular meteo files
c
c revision log :
c
c 10.03.2009    ggu     finished coding
c 24.03.2009    ggu     use metrain instead rqdsv
c 07.05.2009    ggu     new call to init_coords()
c 18.06.2009    ggu&dbf bug fix -> wind speed not interpolated on metws
c 23.02.2010    ggu	call to wstress changed (new wxv,wyv)
c
c*********************************************************************

	subroutine meteo_regular

c administers regular meteo files
c
c in order to use these routines, please set imreg = 1 in the STR file

	implicit none

	integer ndim
	parameter (ndim=400)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	real xgeov(1), ygeov(1)
	common /xgeov/xgeov, /ygeov/ygeov

        real tauxnv(1),tauynv(1)
        common /tauxnv/tauxnv,/tauynv/tauynv
        real wxv(1),wyv(1)
        common /wxv/wxv,/wyv/wyv
        real ppv(1)
        common /ppv/ppv

c arrays containing meteo forcings on FEM grid

        real metrad(1),methum(1)
        real mettair(1),metcc(1)
        real metwbt(1),metws(1)
	real metrain(1)
	common /metrad/metrad, /methum/methum
	common /mettair/mettair, /metcc/metcc
	common /metwbt/metwbt, /metws/metws
	common /metrain/metrain

c arrays containing meteo forcings on regular grid
c
c i-arrays contain info on array
c	1: iunit  2:itold  3: itnew  4: itact  5: nvar  6: nx  7: ny
c
c d-arrays contain data
c 2-dim arrays (.,3) contain old and new and actual time information
c
c dheat		radiation (1), Tair (2), aux (3), cloud cover (4)
c		where aux is either rel. humidity or wetbulb temperature
c dwind		windx (1), windy (2), atm. pressure (3)
c drain		rain (1)
c dextra	humidity (1) and wetbulb temperature (2)
c dws		wind speed (1)

	integer iheat(7)
	integer iwind(7)
	integer irain(7)
	real dheat(4*ndim,3)
	real dwind(3*ndim,3)
	real drain(1*ndim,3)
	real dextra(2*ndim)
	real dws(ndim)
	save iheat,dheat,iwind,dwind,irain,drain,dextra,dws

	character*60 windfile,heatfile,rainfile
	save windfile,heatfile,rainfile

	integer nvar,nx,ny
	integer mode,i
	integer nvarm,nid,nlev
	integer fuse
	real x0,y0,dx,dy
	double precision xtrans,ytrans
	double precision lon0,lat0,phi
	double precision c_param(5)
	real flag
	real val0
	real valmin,valmax

c values for debug: set nfreq and iumet > 0 for debug output
c	nfreq	debug output frequency
c	iumet	debug output unit

	integer iumet,nfreq
	save iumet,nfreq
	data iumet,nfreq / 0 , 0 /

	integer icall
	save icall
	data icall / 0 /

	if( icall .lt. 0 ) return

c------------------------------------------------------------------
c initialization
c------------------------------------------------------------------

	if( icall .eq. 0 ) then

	  write(6,*) 'initialization of regular surface forcing'

c	  ---------------------------------------------------------
c	  initialization of geometry - please customize
c	  ---------------------------------------------------------

c	  ---------------------------------------------------------
c	  geometry of regular matrix: lower left point and step size
c	  ---------------------------------------------------------

	  x0 = 11.
	  y0 = 38.
	  dx = 0.5
	  dy = 0.5
	  flag = -999.

	  call setgeo(x0,y0,dx,dy,flag)

c	  ---------------------------------------------------------
c	  projection to be used: 0: none  1: Gauss-Boaga  2: UTM  3: basic CPP
c	  ---------------------------------------------------------

	  mode = 0	!no projection

	  mode = 3	!basic
          lon0 = 14.047
          lat0 = 35.672
          phi = 36.
          c_param(1) = lon0
          c_param(2) = lat0
          c_param(3) = phi

	  mode = 1	!gauss-boaga
	  fuse = 2
	  xtrans = 2798930.
	  ytrans = 4617760.
          c_param(1) = fuse
          c_param(2) = xtrans
          c_param(3) = ytrans

	  call init_coords(mode,c_param)

	  mode = 1	!from cartesian to lat/lon
	  call convert_coords(mode,nkn,xgv,ygv,xgeov,ygeov)

c	  ---------------------------------------------------------
c	  initialization of data files
c	  ---------------------------------------------------------

	  call getfnm('qflux',heatfile)
	  call getfnm('wind',windfile)
	  call getfnm('rain',rainfile)
	  iwind(1) = -1				!flag as not open
	  iheat(1) = -1
	  irain(1) = -1

	  val0 = 0.
	  call meteo_init_array(nkn,val0,tauxnv)
	  call meteo_init_array(nkn,val0,tauynv)
	  call meteo_init_array(nkn,val0,ppv)
	  call meteo_init_array(nkn,val0,metrad)
	  call meteo_init_array(nkn,val0,mettair)
	  call meteo_init_array(nkn,val0,metcc)
	  call meteo_init_array(nkn,val0,methum)
	  call meteo_init_array(nkn,val0,metwbt)
	  call meteo_init_array(nkn,val0,metrain)
	end if

c------------------------------------------------------------------
c end of initialization - do not change anything beyond this point
c------------------------------------------------------------------

	icall = icall + 1

c------------------------------------------------------------------
c time interpolation
c------------------------------------------------------------------

        call rgf_admin(it,windfile,3,3*ndim,iwind,dwind)
        call rgf_admin(it,heatfile,4,4*ndim,iheat,dheat)
        call rgf_admin(it,rainfile,1,1*ndim,irain,drain)

c------------------------------------------------------------------
c convert some variables to needed ones (humidity, wet bulb, wind speed)
c------------------------------------------------------------------

	mode = 1	! 1: hum -> wetbulb   2: wetbulb -> hum
	call meteo_convert_hum(mode,flag,iheat,dheat,dextra)

	call meteo_convert_ws(flag,iwind,dwind,dws)

c------------------------------------------------------------------
c spatial interpolation of wind data
c------------------------------------------------------------------

        nvar  = iwind(5)
        nx    = iwind(6)
        ny    = iwind(7)

	call meteo_intp(1,nvar,nx,ny,dwind(1,3),nkn,xgeov,ygeov,wxv)
	call meteo_intp(2,nvar,nx,ny,dwind(1,3),nkn,xgeov,ygeov,wyv)
	call meteo_intp(3,nvar,nx,ny,dwind(1,3),nkn,xgeov,ygeov,ppv)

	nvar = 1
	call meteo_intp(1,nvar,nx,ny,dws,nkn,xgeov,ygeov,metws)

	call wstress(nkn,wxv,wyv,tauxnv,tauynv)

c------------------------------------------------------------------
c spatial interpolation of heat data
c------------------------------------------------------------------

        nvar  = iheat(5)
        nx    = iheat(6)
        ny    = iheat(7)

	call meteo_intp(1,nvar,nx,ny,dheat(1,3),nkn,xgeov,ygeov,metrad)
	call meteo_intp(2,nvar,nx,ny,dheat(1,3),nkn,xgeov,ygeov,mettair)
	!call meteo_intp(3,nvar,nx,ny,dheat(1,3),nkn,xgeov,ygeov,methum)
	call meteo_intp(4,nvar,nx,ny,dheat(1,3),nkn,xgeov,ygeov,metcc)

	nvar = 2
	call meteo_intp(1,nvar,nx,ny,dextra,nkn,xgeov,ygeov,methum)
	call meteo_intp(2,nvar,nx,ny,dextra,nkn,xgeov,ygeov,metwbt)

c------------------------------------------------------------------
c spatial interpolation of rain data
c------------------------------------------------------------------

        nvar  = irain(5)
        nx    = irain(6)
        ny    = irain(7)

	call meteo_intp(1,nvar,nx,ny,drain(1,3),nkn,xgeov,ygeov,metrain)

c------------------------------------------------------------------
c debug output
c------------------------------------------------------------------

	if( nfreq .gt. 0 .and. mod(icall,nfreq) .eq. 0 ) then
	  nvarm = 7	!total number of vars that are written
	  nid = 600
	  nlev = 1
	  call conwrite(iumet,'.met',nvarm,nid+1,nlev,ppv)
	  call conwrite(iumet,'.met',nvarm,nid+2,nlev,metws)
	  call conwrite(iumet,'.met',nvarm,nid+3,nlev,metrad)
	  call conwrite(iumet,'.met',nvarm,nid+4,nlev,mettair)
	  call conwrite(iumet,'.met',nvarm,nid+5,nlev,methum)
	  call conwrite(iumet,'.met',nvarm,nid+6,nlev,metcc)
	  call conwrite(iumet,'.met',nvarm,nid+7,nlev,metrain)
	end if

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	end

c*********************************************************************

	subroutine meteo_intp(ivar,nvar,nx,ny,dfile,nkn,xv,yv,valfem)

c interpolates files spatially

	implicit none

	integer ivar
	integer nvar
	integer nx,ny
	real dfile(nx,ny,nvar)
	integer nkn
	real xv(1), yv(1)
	real valfem(1)

	integer k
	real val,flag
	real x,y
	real am2val

	if( nvar .le. 0 ) return	!no data has been read

	!write(6,*) 'meteo_intp: ',nvar,ivar,nx,ny

	if( ivar .le. 0 .or. ivar .gt. nvar ) then
	  write(6,*) 'ivar,nvar: ',ivar,nvar
	  stop 'error stop meteo_intp: ivar out of range'
	end if

	call getgeoflag(flag)

	do k=1,nkn
	  x = xv(k)
	  y = yv(k)
	  val = am2val(dfile(1,1,ivar),nx,ny,x,y)
	  !if( mod(k,100) .eq. 0 ) write(6,*) x,y,valfem(k)
	  if( val .eq. flag ) then
	    write(6,*) nx,ny,x,y,val
	    stop 'error stop meteo_intp: flag found'
	  end if
	  valfem(k) = val
	end do

	end

c*********************************************************************

	subroutine meteo_init_array(n,val0,array)

c initializes 2D array

	implicit none

	integer n
	real val0
	real array(1)

	integer i

	do i=1,n
	  array(i) = val0
	end do

	end

c*********************************************************************

	subroutine meteo_convert_ws(flag,ifile,dfile,dws)

c interpolates files spatially

	implicit none

	real flag
	integer ifile(7)
	real dfile(1)
	real dws(1)

	integer n,nx,ny
	integer iuw,ivw,iws
	integer i
	real uw,vw,ws

        nx    = ifile(6)
        ny    = ifile(7)
	n = nx * ny

	iuw = (1-1) * n		!1. variable in dfile
	ivw = (2-1) * n		!2. variable in dfile
	iws = (1-1) * n		!1. variable in dws

	do i=1,n
	  uw = dfile(iuw+i)
	  vw = dfile(ivw+i)
	  if( uw .eq. flag .or. vw .eq. flag ) then		!not valid
	      ws = flag
	  else
	      ws = sqrt( uw*uw + vw*vw )
	  end if
	  dws(iws+i) = ws
	end do

	end

c*********************************************************************

	subroutine meteo_convert_hum(mode,flag,ifile,dfile,dextra)

c interpolates files spatially

	implicit none

	integer mode
	real flag
	integer ifile(7)
	real dfile(1)
	real dextra(1)

	integer n,nx,ny
	integer idb,ival,irh,iwb
	integer i
	real db,val,rh,wb

        nx    = ifile(6)
        ny    = ifile(7)
	n = nx * ny

	idb = (2-1) * n		!2. variable in dfile
	ival = (3-1) * n	!3. variable in dfile
	irh = (1-1) * n		!1. variable in dextra
	iwb = (2-1) * n		!2. variable in dextra

	do i=1,n
	  db = dfile(idb+i)
	  val = dfile(ival+i)
	  if( val .eq. flag ) then		!not valid
	      rh = flag
	      wb = flag
	  else if( mode .eq. 1 ) then		!val is humidity
	      rh = val
	      call rh2twb(db,rh,wb)
	  else if( mode .eq. 2 ) then		!val is wet bulb
	      wb = val
	      call twb2rh(db,wb,rh)
	  else
	      write(6,*) 'mode = ',mode
	      stop 'error stop meteo_convert_hum: mode'
	  end if
	  dextra(irh+i) = rh
	  dextra(iwb+i) = wb
	end do

	end

c*********************************************************************

	subroutine meteo_set_matrix(qs,ta,rh,wb,uw,cc)

c interpolates files spatially

	implicit none

	real qs,ta,rh,wb,uw,cc

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real metrad(1),methum(1)
        real mettair(1),metcc(1)
        real metwbt(1),metws(1)
	common /metrad/metrad, /methum/methum
	common /mettair/mettair, /metcc/metcc
	common /metwbt/metwbt, /metws/metws

	integer k

	do k=1,nkn
	  metrad(k) = qs
	  mettair(k) = ta
	  methum(k) = rh
	  metwbt(k) = wb
	  metws(k) = uw
	  metcc(k) = cc
	end do

	end

c*********************************************************************

	subroutine meteo_get_values(k,qs,ta,rh,wb,uw,cc,p)

c interpolates files spatially

	implicit none

	integer k
	real qs,ta,rh,wb,uw,cc,p

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real metrad(1),methum(1)
        real mettair(1),metcc(1)
        real metwbt(1),metws(1)
	common /metrad/metrad, /methum/methum
	common /mettair/mettair, /metcc/metcc
	common /metwbt/metwbt, /metws/metws

        real ppv(1)
        common /ppv/ppv

	qs = metrad(k)
	ta = mettair(k)
	rh = methum(k)
	wb = metwbt(k)
	uw = metws(k)
	cc = metcc(k)

	p = ppv(k)
	if( p .lt. 800. .or. p .gt. 1100. ) p = 1013.	!850.0 - 1085.6 mb

	end

c*********************************************************************

