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
c 26.01.2011    ggu	write wind field to debug output (iumetw)
c 05.02.2011    ggu	changed order of records in dfile (actual is first)
c 16.02.2011    ggu	pass idata to files, use geo info from files
c 18.11.2011    ggu	deleted projection code from subroutines
c 10.02.2012    ggu	limit cc and rh to acceptable range
c 16.02.2012    ggu	new routine meteo_get_solar_radiation()
c 22.02.2012    ggu	new routines for regular and ts reading of meteo
c 23.02.2012    ggu&ccf	bug fix meteo_copy_to_old and meteo_interpolate_in_time
c
c notes :
c
c info on file format can be found in subrgf.f
c
c to do: implement zdist for rain (constant rain)
c
c*********************************************************************

	subroutine meteo_forcing

c administers meteo files (still to be cleaned)
c
c in order to use these routines, please set imreg = 1 in the STR file

	implicit none

	integer ndim
	parameter (ndim=4650)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	include 'meteo.h'

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

	integer ifidim
	parameter( ifidim = 30 )

	integer iheat(ifidim)
	integer iwind(ifidim)
	integer irain(ifidim)

	real dheat(4*ndim,3)
	real dwind(3*ndim,3)
	real drain(1*ndim,3)
	real dextra(2*ndim)
	real dws(ndim)
	save iheat,dheat,iwind,dwind,irain,drain,dextra,dws

	character*60 windfile,heatfile,rainfile
	save windfile,heatfile,rainfile

	integer nvar,nx,ny,i,iunit
	integer mode,iproj,modehum
	integer nvarm,nid,nlev
	real x0,y0,dx,dy
	real flag
	real val0
	real valmin,valmax

	integer ifileo

c values for debug: set nfreq and iumet > 0 for debug output
c	nfreq	debug output frequency
c	iumet	debug output unit

	integer iumet,iumetw,nfreq
	save iumet,iumetw,nfreq
	data iumet,iumetw,nfreq / 0 , 0 , 0 /

	integer icall
	save icall
	data icall / 0 /

	if( icall .lt. 0 ) return

c------------------------------------------------------------------
c initialization
c------------------------------------------------------------------

	if( icall .eq. 0 ) then

	  write(6,*) 'initialization of meteo forcing'

c	  ---------------------------------------------------------
c	  initialization of data files
c	  ---------------------------------------------------------

	  val0 = 0.
	  call meteo_set_array(3*nkn,val0,wxv)
	  call meteo_set_array(3*nkn,val0,wyv)
	  call meteo_set_array(3*nkn,val0,ppv)
	  call meteo_set_array(3*nkn,val0,metrad)
	  call meteo_set_array(3*nkn,val0,mettair)
	  call meteo_set_array(3*nkn,val0,metcc)
	  call meteo_set_array(3*nkn,val0,methum)
	  call meteo_set_array(3*nkn,val0,metrain)

	  call meteo_set_array(nkn,val0,tauxnv)
	  call meteo_set_array(nkn,val0,tauynv)
	  call meteo_set_array(nkn,val0,metwbt)
	  call meteo_set_array(nkn,val0,metws)

	  call getfnm('qflux',heatfile)
	  call getfnm('wind',windfile)
	  call getfnm('rain',rainfile)

	  do i=1,ifidim
	    iwind(i) = 0
	    iheat(i) = 0
	    irain(i) = 0
	  end do

	  call meteo_wind_init(windfile,iwind,it)
	  call meteo_heat_init(heatfile,iheat,it)
	  call meteo_rain_init(rainfile,irain,it)

	end if

c------------------------------------------------------------------
c end of initialization
c------------------------------------------------------------------

	icall = icall + 1

c------------------------------------------------------------------
c time interpolation
c------------------------------------------------------------------

	call meteo_wind_admin(iwind,it)
	call meteo_heat_admin(iheat,it)
	call meteo_rain_admin(irain,it)

c------------------------------------------------------------------
c extra treatment of wind data
c------------------------------------------------------------------

	iunit = iwind(1)
	if( iunit .gt. 0 ) then
	  call meteo_compute_ws(iwind)
	  call wstress(nkn,wxv,wyv,tauxnv,tauynv)
	end if

c------------------------------------------------------------------
c extra treatment of heat data
c------------------------------------------------------------------

	iunit = iheat(1)
	if( iunit .gt. 0 ) then
	  !if( iwind(1) .le. 0 ) goto 99	!we need wind data for heat
	  modehum = 1	! 1: hum -> wetbulb   2: wetbulb -> hum
	  call meteo_compute_wbt(modehum,iheat)
	end if

c------------------------------------------------------------------
c extra treatment of rain data
c------------------------------------------------------------------

c------------------------------------------------------------------
c debug output
c------------------------------------------------------------------

	if( nfreq .gt. 0 .and. icall .eq. 1 ) then
	  iumetw = ifileo(0,'wind_debug.win','unform','new')
	end if
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
	  call wrwin(iumetw,it,nkn,wxv,wyv,ppv)
	end if

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	return
   99	continue
	stop 'error stop meteo_regular: need wind data for heat module'
	end

c*********************************************************************
c*********************************************************************
c*********************************************************************

	subroutine meteo_rain_init(file,ifile,it)

	implicit none

	character*(*) file
	integer ifile(1)
	integer it

	call meteo_init_file(file,ifile,1)
	call meteo_rain_next_record(ifile)	!read first record
	call meteo_rain_copy(ifile)
	call meteo_rain_admin(ifile,it)

	end

c*********************************************************************

	subroutine meteo_rain_admin(ifile,it)

	implicit none

	integer ifile(*)
	integer it

	integer itnew
	logical meteo_is_open

	if( .not. meteo_is_open(ifile) ) return

	itnew = ifile(3)

	do while( it .gt. itnew )
	  call meteo_rain_copy(ifile)
	  call meteo_rain_next_record(ifile)
	  itnew = ifile(3)
	end do

	call meteo_rain_interpolate(ifile,it)

	end

c*********************************************************************

	subroutine meteo_rain_copy(ifile)

	implicit none

	include 'meteo.h'

	integer ifile(*)

	integer n

	n = ifile(10)

	call meteo_copy_to_old(ifile,n,metrain)

	end

c*********************************************************************

	subroutine meteo_rain_next_record(ifile)

	implicit none

	include 'meteo.h'

	integer ifile(*)

	integer mode,iunit,ne,it,ierr,ip
	logical meteo_is_open
	logical bdata
	real data(9)

	if( .not. meteo_is_open(ifile) ) return

	ierr = 0
	iunit = ifile(1)
	mode = ifile(9)
	ne = ifile(10)
	ip = 2*ne+1

	if( mode .eq. 1 ) then
	  call read_rain_unformatted(iunit,it,ne,metrain(ip),ierr)
	else if( mode .eq. 3 ) then
	  call rgf_read(iunit,nmtdim,it,ifile,mdata,bdata)
	  if( bdata ) then
	    call meteo_fem_interpolate(1,ifile,mdata,metrain(ip))
	  else
	    ierr = -1
	  end if
	else if( mode .eq. 4 ) then
	  call ts_read(it,ifile,data,bdata)
	  if( bdata ) then
	    call meteo_set_array(ne,data(1),metrain(ip))
	  else
	    ierr = -1
	  end if
	else
	  goto 99
	end if

	if( ierr .lt. 0 ) then		!EOF -> close file
	  call meteo_close_file(ifile)
	  write(6,*) 'no more records for rain file... closing'
	else
	  ifile(3) = it
	  write(6,*) 'new record for rain file: ',it
	end if

	return
   99	continue
	write(6,*) 'mode = ',mode
	stop 'error stop meteo_rain_next_record: mode not supported...'
	end

c*********************************************************************

	subroutine meteo_rain_interpolate(ifile,it)

	implicit none

	include 'meteo.h'

	integer ifile(*)
	integer it

	integer n

	n = ifile(10)

	call meteo_interpolate_in_time(ifile,it,n,metrain)

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************

	subroutine meteo_heat_init(file,ifile,it)

	implicit none

	character*(*) file
	integer ifile(*)
	integer it

	call meteo_init_file(file,ifile,5)
	call meteo_heat_next_record(ifile)	!read first record
	call meteo_heat_copy(ifile)
	call meteo_heat_admin(ifile,it)

	end

c*********************************************************************

	subroutine meteo_heat_admin(ifile,it)

	implicit none

	integer ifile(*)
	integer it

	integer itnew
	logical meteo_is_open

	if( .not. meteo_is_open(ifile) ) return

	itnew = ifile(3)

	do while( it .gt. itnew )
	  call meteo_heat_copy(ifile)
	  call meteo_heat_next_record(ifile)
	  itnew = ifile(3)
	end do

	call meteo_heat_interpolate(ifile,it)

	end

c*********************************************************************

	subroutine meteo_heat_copy(ifile)

	implicit none

	include 'meteo.h'

	integer ifile(*)

	integer n

	n = ifile(10)

	call meteo_copy_to_old(ifile,n,metrad)
	call meteo_copy_to_old(ifile,n,mettair)
	call meteo_copy_to_old(ifile,n,methum)
	call meteo_copy_to_old(ifile,n,metcc)

	end

c*********************************************************************

	subroutine meteo_heat_next_record(ifile)

	implicit none

	include 'meteo.h'

	integer ifile(*)

	integer mode,iunit,ne,it,ierr,ip
	logical meteo_is_open
	logical bdata
	real data(9)

	if( .not. meteo_is_open(ifile) ) return

	ierr = 0
	iunit = ifile(1)
	mode = ifile(9)
	ne = ifile(10)
	ip = 2*ne+1

	if( mode .eq. 1 ) then
	  call read_heat_unformatted(iunit,it,ne
     +		,metrad(ip),mettair(ip),methum(ip),metcc(ip),ierr)
	else if( mode .eq. 3 ) then
	  call rgf_read(iunit,nmtdim,it,ifile,mdata,bdata)
	  if( bdata ) then
	    call meteo_fem_interpolate(1,ifile,mdata,metrad(ip))
	    call meteo_fem_interpolate(2,ifile,mdata,mettair(ip))
	    call meteo_fem_interpolate(3,ifile,mdata,methum(ip))
	    call meteo_fem_interpolate(4,ifile,mdata,metcc(ip))
	  else
	    ierr = -1
	  end if
	else if( mode .eq. 4 ) then
	  call ts_read(it,ifile,data,bdata)
	  if( bdata ) then
	    call meteo_set_array(ne,data(1),metrad(ip))
	    call meteo_set_array(ne,data(2),mettair(ip))
	    call meteo_set_array(ne,data(3),methum(ip))
	    call meteo_set_array(ne,data(3),metcc(ip))
	  else
	    ierr = -1
	  end if
	else
	  goto 99
	end if

	if( ierr .lt. 0 ) then		!EOF -> close file
	  call meteo_close_file(ifile)
	  write(6,*) 'no more records for heat file... closing'
	else
	  ifile(3) = it
	  write(6,*) 'new record for heat file: ',it
	end if

	return
   99	continue
	write(6,*) 'mode = ',mode
	stop 'error stop meteo_heat_next_record: mode not supported...'
	end

c*********************************************************************

	subroutine meteo_heat_interpolate(ifile,it)

	implicit none

	include 'meteo.h'

	integer ifile(*)
	integer it

	integer n

	n = ifile(10)

	call meteo_interpolate_in_time(ifile,it,n,metrad)
	call meteo_interpolate_in_time(ifile,it,n,mettair)
	call meteo_interpolate_in_time(ifile,it,n,methum)
	call meteo_interpolate_in_time(ifile,it,n,metcc)

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************

	subroutine meteo_wind_init(file,ifile,it)

	implicit none

	character*(*) file
	integer ifile(*)
	integer it

	call meteo_init_file(file,ifile,2)
	call meteo_wind_next_record(ifile)	!read first record
	call meteo_wind_copy(ifile)
	call meteo_wind_admin(ifile,it)

	end

c*********************************************************************

	subroutine meteo_wind_admin(ifile,it)

	implicit none

	integer ifile(*)
	integer it

	integer itnew
	logical meteo_is_open

	if( .not. meteo_is_open(ifile) ) return

	itnew = ifile(3)

	do while( it .gt. itnew )
	  call meteo_wind_copy(ifile)
	  call meteo_wind_next_record(ifile)
	  itnew = ifile(3)
	end do

	call meteo_wind_interpolate(ifile,it)

	end

c*********************************************************************

	subroutine meteo_wind_copy(ifile)

	implicit none

	include 'meteo.h'

	integer ifile(*)

	integer n

	n = ifile(10)

	call meteo_copy_to_old(ifile,n,wxv)
	call meteo_copy_to_old(ifile,n,wyv)
	call meteo_copy_to_old(ifile,n,ppv)

	end

c*********************************************************************

	subroutine meteo_wind_next_record(ifile)

	implicit none

	include 'meteo.h'

	integer ifile(*)

	integer mode,iunit,ne,it,ierr,ip
	logical meteo_is_open
	logical bdata
	real data(9)

	if( .not. meteo_is_open(ifile) ) return

	ierr = 0
	iunit = ifile(1)
	mode = ifile(9)
	ne = ifile(10)
	ip = 2*ne+1

	if( mode .eq. 1 ) then
	  call read_wind_unformatted(iunit,it,ne
     +		,wxv(ip),wyv(ip),ppv(ip),ierr)
	else if( mode .eq. 3 ) then
	  call rgf_read(iunit,nmtdim,it,ifile,mdata,bdata)
	  if( bdata ) then
	    call meteo_fem_interpolate(1,ifile,mdata,wxv(ip))
	    call meteo_fem_interpolate(2,ifile,mdata,wyv(ip))
	    call meteo_fem_interpolate(3,ifile,mdata,ppv(ip))
	  else
	    ierr = -1
	  end if
	else if( mode .eq. 4 ) then
	  call ts_read(it,ifile,data,bdata)
	  if( bdata ) then
	    call meteo_set_array(ne,data(1),wxv(ip))
	    call meteo_set_array(ne,data(2),wyv(ip))
	    call meteo_set_array(ne,data(3),ppv(ip))
	  else
	    ierr = -1
	  end if
	else
	  goto 99
	end if

	if( ierr .lt. 0 ) then		!EOF -> close file
	  call meteo_close_file(ifile)
	  write(6,*) 'no more records for wind file... closing'
	else
	  ifile(3) = it
	  write(6,*) 'new record for wind file: ',it
	end if

	return
   99	continue
	write(6,*) 'mode = ',mode
	stop 'error stop meteo_wind_next_record: mode not supported...'
	end

c*********************************************************************

	subroutine meteo_wind_interpolate(ifile,it)

	implicit none

	include 'meteo.h'

	integer ifile(*)
	integer it

	integer n

	n = ifile(10)

	call meteo_interpolate_in_time(ifile,it,n,wxv)
	call meteo_interpolate_in_time(ifile,it,n,wyv)
	call meteo_interpolate_in_time(ifile,it,n,ppv)

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************

	function meteo_is_open(ifile)

	implicit none

	logical meteo_is_open
	integer ifile(*)

	meteo_is_open = ifile(1) .gt. 0

	end

c*********************************************************************

	subroutine meteo_close_file(ifile)

	implicit none

	integer ifile(*)

	integer iunit

	iunit = ifile(1)

	ifile(1) = 0
	ifile(9) = 0

	if( iunit .gt. 0 ) close(iunit)

	end

c*********************************************************************

	subroutine meteo_init_file(file,ifile,nvar)

	implicit none

	character*(*) file
	integer ifile(*)
	integer it
	integer nvar		!expected variables (only needed for timeseries)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer iunit,mode
	integer ifileo
	logical is_meteo_unformatted
	logical is_meteo_regular
	logical is_meteo_ts

	ifile(9) = -1
	if( file .eq. ' ' ) return

	if( is_meteo_unformatted(file,0) ) then
	  iunit = ifileo(0,file,'unform','old')
	  mode = 1
	else if( is_meteo_regular(file,0) ) then
	  iunit = ifileo(0,file,'form','old')
	  mode = 3
	else if( is_meteo_ts(file,nvar) ) then
	  iunit = ifileo(0,file,'form','old')
	  mode = 4
	else
	  stop 'error stop meteo_init_file: not a known format'
	end if

	if( iunit .le. 0 ) goto 99

	ifile(1) = iunit
	ifile(5) = nvar		!works only for ts
	ifile(9) = mode
	ifile(10) = nkn

	write(6,*) 'meteo file opened: ',mode,file

	return
   99	continue
	write(6,*) file
	stop 'error stop meteo_init_file: error opening file'
	end

c*********************************************************************

	subroutine meteo_copy_to_old(ifile,n,array)

	implicit none

	integer ifile(*)
	integer n
	real array(n,3)

	integer i

	ifile(2) = ifile(3)

	do i=1,n
	  array(i,2)=array(i,3)
	end do

	end

c*********************************************************************

	subroutine meteo_interpolate_in_time(ifile,it,n,array)

	implicit none

	integer ifile(*)
	integer it
	integer n
	real array(n,3)

	integer itact,itold,itnew
	integer i
	real rit

	itold = ifile(2)
	itnew = ifile(3)
	itact = it

        if( itnew .gt. itold ) then
          rit=float(itact-itold)/float(itnew-itold)
        else
	  rit = 0
	end if

	do i=1,n
	  array(i,1)=rit*(array(i,3)-array(i,2))+array(i,2)
	end do

	ifile(4) = itact

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************

	subroutine meteo_set_array(n,val0,array)

c initializes 2D array

	implicit none

	integer n
	real val0
	real array(*)

	integer i

	do i=1,n
	  array(i) = val0
	end do

	end

c*********************************************************************

	subroutine meteo_compute_ws(ifile)

c computes wind speed

	implicit none

	include 'meteo.h'

	integer ifile(*)

	integer n,i
	real uw,vw,ws

	n = ifile(10)

	do i=1,n
	  uw = wxv(i)
	  vw = wyv(i)
	  metws(i) = sqrt( uw*uw + vw*vw )
	end do

	end

c*********************************************************************

	subroutine meteo_compute_wbt(mode,ifile)

c computes wet bulb temperature

	implicit none

	include 'meteo.h'

	integer mode
	integer ifile(*)

	integer n,i
	real db,rh,wb

	n = ifile(10)

	do i=1,n
	  db = mettair(i)
	  rh = methum(i)
	  wb = metwbt(i)
	  if( mode .eq. 1 ) then		!val is humidity
	      call rh2twb(db,rh,wb)
	  else if( mode .eq. 2 ) then		!val is wet bulb
	      call twb2rh(db,wb,rh)
	  else
	      write(6,*) 'mode = ',mode
	      stop 'error stop meteo_convert_hum: mode'
	  end if
	  methum(i) = rh
	  metwbt(i) = wb
	end do

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************

	subroutine meteo_set_heat_values(qs,ta,rh,wb,uw,cc)

c interpolates files spatially

	implicit none

	include 'meteo.h'

	real qs,ta,rh,wb,uw,cc

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

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

	subroutine meteo_get_heat_values(k,qs,ta,rh,wb,uw,cc,p)

c returns meteo parameters for one node
c
c pressure is returned in [mb]

	implicit none

	include 'meteo.h'

        integer k                       !node number
        real qs                         !solar radiation [W/m**2]
        real ta                         !air temperature [Celsius]
        real rh                         !relative humidity [%, 0-100]
        real wb                         !wet bulb temperature [Celsius]
        real uw                         !wind speed [m/s]
        real cc                         !cloud cover [0-1]
        real p                          !atmospheric pressure [mbar, hPa]

	real pstd
	parameter ( pstd = 1013.25 )

	qs = metrad(k)
	ta = mettair(k)
	rh = methum(k)
	wb = metwbt(k)
	uw = metws(k)
	cc = metcc(k)

	cc = max(0.,cc)
	cc = min(1.,cc)
	rh = max(0.,rh)
	rh = min(100.,rh)

	p = 0.01 * ppv(k)				!Pascal to mb
	if( p .lt. 800. .or. p .gt. 1100. ) p = pstd	!850.0 - 1085.6 mb

	end

c*********************************************************************

	subroutine meteo_get_solar_radiation0(k,qs)

	implicit none

	include 'meteo.h'

        integer k                       !node number
        real qs                         !solar radiation [W/m**2]

	qs = metrad(k)

	end

c*********************************************************************

