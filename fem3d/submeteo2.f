c
c $Id: submeteo.f,v 1.7 2010-02-26 17:35:06 georg Exp $
c
c handle meteo files with new fem format
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
c 20.05.2014    ggu	new routines for new file format
c
c notes :
c
c this routine is called with imreg = 3
c
c info on file format can be found in subfemfile.f
c
c*********************************************************************

!================================================================
        module meteo_forcing_module
!================================================================

	use intp_fem_file

	implicit none

	real, parameter :: pstd = 101325.	!standard pressure
	real, parameter :: dstd = 2.5e-3	!standard drag coefficient
	real, parameter :: nmile = 1852.	!nautical mile in m

	integer, save :: idwind,idheat,idrain

	integer, parameter :: nfreq = 0		!debug output
	integer, save :: iumet = 0

	integer, save :: iwtype,itdrag
	integer, save :: irtype
	integer, save :: ihtype
	real, save :: wsmax,dragco,roluft,rowass
	real, save :: pfact = 1.
	real, save :: wfact = 1.
	real, save :: sfact = 1.

	logical, save, private :: bdebug = .true.

	integer, save :: icall = 0

	character*80, save :: wxss = 'wind stress in x [N/m**2]'
	character*80, save :: wyss = 'wind stress in y [N/m**2]'
	character*80, save :: wxms = 'wind velocity in x [m/s]'
	character*80, save :: wyms = 'wind velocity in y [m/s]'
	character*80, save :: wsms = 'wind speed [m/s]'
	character*80, save :: wskn = 'wind speed [knots]'
	character*80, save :: wdir = 'wind direction [deg]'

	character*80, save :: papa = 'pressure (atmospheric) [Pa]'
	character*80, save :: pamb = 'pressure (atmospheric) [mbar]'

	character*80, save :: rain = 'rain [mm/day]'

        character*80, save :: srad = 'solar radiation [W/m**2]'
        character*80, save :: tair = 'air temperature [C]'
        character*80, save :: rhum = 'humidity [%]'
        character*80, save :: ccov = 'cloud cover [0-1]'
        character*80, save :: wbtm = 'wet bulb temperature [C]'

!================================================================
        contains
!================================================================

!*********************************************************************

	subroutine meteo_forcing_fem

! administers meteo files (still to be cleaned)
!
! in order to use these routines, please set imreg = 3 in the STR file

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer ilhkv(1)	!delete once iff_init_global is in main
        common /ilhkv/ilhkv
        real hlv(1)
        common /hlv/hlv
        real hkv(1)
        common /hkv/hkv

	include 'meteo.h'

	character*60 windfile,heatfile,rainfile
	character*4 what

	integer nvar,lmax
	integer nintp
	integer modehum
	integer nvarm,nid,nlev
	!integer it0
	double precision dtime0,dtime
	real flag
	real val0,val1,val2

	real vconst(4)
	integer nodes(1)

	if( icall .lt. 0 ) return

!------------------------------------------------------------------
! initialization
!------------------------------------------------------------------

	if( icall .eq. 0 ) then

	  write(6,*) 'initialization of meteo forcing fem'

!	  ---------------------------------------------------------
!	  initialization of data files
!	  ---------------------------------------------------------

	  !call iff_init_global(nkn,nlv,ilhkv,hkv,hlv)	!should go to main

	  call getfnm('wind',windfile)
	  call getfnm('qflux',heatfile)
	  call getfnm('rain',rainfile)

	  !it0 = itanf
	  dtime0 = itanf

	  nvar = 3
	  nintp = 2
	  what = 'wind'
	  vconst = (/ 0., 0., pstd, 0. /)
	  call iff_init(dtime0,windfile,nvar,nkn,0,nintp
     +				,nodes,vconst,idwind)

	  call meteo_set_wind_data(idwind,nvar)

	  nvar = 1
	  nintp = 2
	  what = 'rain'
	  vconst = (/ 0., 0., 0., 0. /)
	  call iff_init(dtime0,rainfile,nvar,nkn,0,nintp
     +				,nodes,vconst,idrain)

	  call meteo_set_rain_data(idrain,nvar)

	  nvar = 4
	  nintp = 2
	  what = 'heat'
	  vconst = (/ 0., 0., 50., 0. /)
	  call iff_init(dtime0,heatfile,nvar,nkn,0,nintp
     +				,nodes,vconst,idheat)

	  call meteo_set_heat_data(idheat,nvar)

	end if

!------------------------------------------------------------------
! end of initialization
!------------------------------------------------------------------

	icall = icall + 1

!------------------------------------------------------------------
! time interpolation
!------------------------------------------------------------------

	dtime = it
	lmax = 1

	if( .not. iff_is_constant(idwind) .or. icall == 1 ) then
	  call iff_time_interpolate(idwind,dtime,1,nkn,lmax,wxv)
	  call iff_time_interpolate(idwind,dtime,2,nkn,lmax,wyv)
	  if( iff_get_nvar(idwind) == 3 ) then
	    call iff_time_interpolate(idwind,dtime,3,nkn,lmax,ppv)
	!if( bdebug) then
	!  call iff_get_value(idwind,3,1,1,1000,val1)
	!  call iff_get_value(idwind,3,2,1,1000,val2)
	!  call iff_get_file_value(idwind,3,1,1000,val0)
	!end if
	  else
	    call meteo_set_array(nkn,pstd,ppv)
	  end if
	end if

	if( .not. iff_is_constant(idrain) .or. icall == 1 ) then
	  call iff_time_interpolate(idrain,dtime,1,nkn,lmax,metrain)
	end if

	if( .not. iff_is_constant(idheat) .or. icall == 1 ) then
	  call iff_time_interpolate(idheat,dtime,1,nkn,lmax,metrad)
	  call iff_time_interpolate(idheat,dtime,2,nkn,lmax,mettair)
	  call iff_time_interpolate(idheat,dtime,3,nkn,lmax,methum)
	  call iff_time_interpolate(idheat,dtime,4,nkn,lmax,metcc)
	end if

!------------------------------------------------------------------
! extra treatment of data
!------------------------------------------------------------------

	if( .not. iff_is_constant(idwind) .or. icall == 1 ) then
	  call meteo_convert_wind_data(idwind,nkn,wxv,wyv
     +			,tauxnv,tauynv,metws,ppv)
	end if

	if( .not. iff_is_constant(idheat) .or. icall == 1 ) then
	  call meteo_convert_heat_data(idheat,nkn
     +			,mettair,methum,metwbt)
	end if

!------------------------------------------------------------------
! debug output
!------------------------------------------------------------------

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

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

	return
   99	continue
	stop 'error stop meteo_regular: need wind data for heat module'
	end subroutine meteo_forcing_fem

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_set_wind_data(id,nvar)

! iwtype = 0 no wind data processed
! iwtype = 1 wind data in format (wx,wy)
! iwtype = 2 wind data in format (tx,ty) (stress)
! iwtype = 3 wind data in format (speed,dir) (speed in m/s)
!          dir = 0  -> wind from north   dir = 90 -> wind from east
! iwtype = 4 wind data in format (speed,dir) (speed in nodes)

	integer id
	integer nvar

	character*60 string,string1,string2

	real, external :: getpar

!	---------------------------------------------------------
!	check nvar and get parameters
!	---------------------------------------------------------

	if( nvar /= 2 .and. nvar /= 3 ) then
	  write(6,*) 'no support for nvar = ',nvar
	  stop 'error stop meteo_set_wind_data: wind'
	end if

        iwtype = nint(getpar('iwtype'))
        itdrag = nint(getpar('itdrag'))
        wsmax = getpar('wsmax')
        dragco = getpar('dragco')
        roluft = getpar('roluft')
        rowass = getpar('rowass')

	if( dragco .lt. 0 ) dragco = dstd

!	---------------------------------------------------------
!	handle wind
!	---------------------------------------------------------

	call iff_get_var_description(id,1,string1)
	call iff_get_var_description(id,2,string2)

	if( string1 == ' ' ) then	!TS file or constant

	  if( iff_has_file(id) ) then
	    if( iwtype .le. 0 ) then
	      write(6,*) 'wind file given but no wind type available'
	      write(6,*) 'either set wind description or set iwtype'
	      stop 'error stop meteo_set_wind_data: no wind type'
	    else if( iwtype > 4 ) then
	      write(6,*) 'not supported value for iwtype = ',iwtype
	      stop 'error stop meteo_set_wind_data: iwtype'
	    end if
	  else				!no file opened
	    iwtype = 0
	  end if

	  if( iwtype == 1 ) then
	    call iff_set_var_description(id,1,wxms)
	    call iff_set_var_description(id,2,wyms)
	  else if( iwtype == 2 ) then
	    call iff_set_var_description(id,1,wxss)
	    call iff_set_var_description(id,2,wyss)
	  else if( iwtype == 3 ) then
	    call iff_set_var_description(id,1,wsms)
	    call iff_set_var_description(id,2,wdir)
	  else if( iwtype == 4 ) then
	    call iff_set_var_description(id,1,wskn)
	    call iff_set_var_description(id,2,wdir)
	  end if

	else	!description given -> check and set iwtype

	  if( string1 == wxms .and. string2 == wyms ) then
	    iwtype = 1
	  else if( string1 == wxss .and. string2 == wyss ) then
	    iwtype = 2
	  else if( string1 == wsms .and. string2 == wdir ) then
	    iwtype = 3
	  else if( string1 == wskn .and. string2 == wdir ) then
	    iwtype = 4
	  else
	    write(6,*) 'description string for wind not recognized: '
	    write(6,*) string1
	    write(6,*) string2
	    stop 'error stop meteo_set_wind_data: wind description'
	  end if

	end if

	wfact = 1. / rowass
	if( iwtype /= 2 ) wfact = roluft / rowass
	sfact = 1.
	if( iwtype == 4 ) sfact = nmile / 3600.

!	---------------------------------------------------------
!	handle pressure
!	---------------------------------------------------------

	if( nvar == 3 ) then
	  call iff_get_var_description(id,3,string)
	  if( string == ' ' ) string = papa
	  if( string == papa ) then
	    pfact = 1.
	  else if( string == pamb ) then
	    pfact = 100.
	  else
	    write(6,*) 'description string for pressure not recognized: '
	    write(6,*) string
	    stop 'error stop meteo_set_wind_data: pressure description'
	  end if
	  if( bdebug ) then
	    write(6,*) 'pressure initialized: ',pfact,string
	    write(6,*) 'wind type: ',iwtype
	  end if
	end if

!	---------------------------------------------------------
!	remember values and write to monitor
!	---------------------------------------------------------

        call putpar('iwtype',real(iwtype))

	if( iwtype == 0 ) then
	  write(6,*) 'no wind file opened'
	else
	  write(6,*) 'wind file opened: ',iwtype
	  call iff_get_var_description(id,1,string1)
	  call iff_get_var_description(id,2,string2)
	  write(6,*) 'content: '
	  write(6,*) ' 1    ',string1
	  write(6,*) ' 2    ',string2
	  if( nvar == 3 ) then
	    call iff_get_var_description(id,3,string)
	    write(6,*) ' 3    ',string
	  end if
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end subroutine meteo_set_wind_data

!*********************************************************************

	subroutine meteo_convert_wind_data(id,n,wx,wy,tx,ty,ws,pp)

	integer id
	integer n
	real wx(n),wy(n)
	real tx(n),ty(n)
	real ws(n)
	real pp(n)

	logical bnowind,bstress,bspeed
	integer k
	real cd,wxymax,txy,wspeed,wdir

	bnowind = iwtype == 0
	bstress = iwtype == 2
	bspeed = iwtype > 2
	cd = dragco
	wxymax = 0.
	
        if( bnowind ) then              !no wind
	  !nothing to be done
        else if( bstress ) then         !data is stress -> normalize it
          if( cd .le. 0 ) cd = dstd
          do k=1,n
            tx(k) = wfact * wx(k)
            ty(k) = wfact * wy(k)
            txy = sqrt( tx(k)**2 + ty(k)**2 )
            wspeed = sqrt(txy/cd)
            wxymax = max(wxymax,wspeed)
            wx(k) = tx(k) / (cd*wspeed)
            wy(k) = ty(k) / (cd*wspeed)
	    ws(k) = wspeed
          end do
        else 
	  if( bspeed ) then             !data is speed/dir
            do k=1,n
	      wspeed = sfact * wx(k)
	      wdir = wy(k)
              call convert_wind(wspeed,wdir,wx(k),wy(k))
	      ws(k) = wspeed
	    end do
	  else				!data is wind velocity [m/s]
            do k=1,n
              wspeed = sqrt( wx(k)**2 + wy(k)**2 )
	      ws(k) = wspeed
	    end do
	  end if

          do k=1,n
            wspeed = ws(k)
            wxymax = max(wxymax,wspeed)
            if( itdrag .gt. 0 ) call get_drag(itdrag,wspeed,cd)
            tx(k) = wfact * cd * wspeed * wx(k)
            ty(k) = wfact * cd * wspeed * wy(k)
          end do
        end if

        if( wxymax .gt. wsmax ) then
          write(6,*) 'maximum wind speed: ',wxymax
          write(6,*) 'maximum allowed wind speed: ',wsmax
          write(6,*) 'Are you sure the wind is in the correct format?'
          write(6,*) 'If no, please set iwtype, else increase wsmax.'
          stop 'error stop meteo_convert_wind_data: wind speed too high'
        end if

	if( pfact /= 1. ) then
          do k=1,n
	    pp(k) = pfact * pp(k)
	  end do
	end if

	end subroutine meteo_convert_wind_data

!*********************************************************************

	subroutine meteo_set_rain_data(id,nvar)

	integer id
	integer nvar

	character*60 string

!	---------------------------------------------------------
!	check nvar and get parameters
!	---------------------------------------------------------

	if( nvar /= 1 ) then
	  write(6,*) 'no support for nvar = ',nvar
	  stop 'error stop meteo_set_rain_data: rain'
	end if

!	---------------------------------------------------------
!	handle rain
!	---------------------------------------------------------

	call iff_get_var_description(id,1,string)

	if( string == ' ' ) then	!TS file or constant
	  irtype = 0
	  if( iff_has_file(id) ) irtype = 1

	  if( irtype == 1 ) then
	    call iff_set_var_description(id,1,rain)
	  end if
	else
	  if( string == rain ) then
	    irtype = 1
	  else
	    write(6,*) 'description string for rain not recognized: '
	    write(6,*) string
	    stop 'error stop meteo_set_rain_data: rain description'
	  end if
	end if

!	---------------------------------------------------------
!	remember values and write to monitor
!	---------------------------------------------------------

	if( irtype == 0 ) then
	  write(6,*) 'no rain file opened'
	else
	  write(6,*) 'rain file opened: ',irtype
	  call iff_get_var_description(id,1,string)
	  write(6,*) 'content: '
	  write(6,*) ' 1    ',string
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end subroutine meteo_set_rain_data

!*********************************************************************

	subroutine meteo_set_heat_data(id,nvar)

	integer id
	integer nvar

	character*60 string,strings(4)
	integer i

!	---------------------------------------------------------
!	check nvar and get parameters
!	---------------------------------------------------------

	if( nvar == 5 ) then
	  write(6,*) 'old format for heatflux not supported any more'
	  write(6,*) 'time series should have only 4 columns: '
	  write(6,*) 'sol_rad t_air humidity cloud_cover'
	  write(6,*) 'no support for nvar = ',nvar
	  stop 'error stop meteo_forcing_fem: heat'
	else if( nvar /= 4 ) then
	  write(6,*) 'no support for nvar = ',nvar
	  stop 'error stop meteo_set_heat_data: heat'
	end if

!	---------------------------------------------------------
!	handle heat
!	---------------------------------------------------------

	do i=1,nvar
	  call iff_get_var_description(id,i,strings(i))
	end do

	if( strings(1) == ' ' ) then	!TS file or constant
	  ihtype = 0
	  if( iff_has_file(id) ) ihtype = 1

	  if( ihtype == 1 ) then
	    call iff_set_var_description(id,1,srad)
	    call iff_set_var_description(id,2,tair)
	    call iff_set_var_description(id,3,rhum)
	    call iff_set_var_description(id,4,ccov)
	  end if
	else
	  do i=1,nvar
	    call iff_get_var_description(id,i,strings(i))
	  end do
	  if( strings(1) == srad ) then
	    if( strings(2) == tair .and. strings(4) == ccov ) then
	      if( strings(3) == rhum ) then
	        ihtype = 1
	      else if( strings(3) == wbtm ) then
	        ihtype = 2
	      else
	        ihtype = -3
	      end if
	    else
	      ihtype = -2
	    end if
	  else
	    ihtype = -1
	  end if
	end if

	if( ihtype < 0 ) then
	  write(6,*) 'description string for heat not recognized: '
	  write(6,*) 'possible number of string with error: ',-ihtype
	  do i=1,nvar
	    call iff_get_var_description(id,i,string)
	    write(6,*) i,'    ',string
	  end do
	  stop 'error stop meteo_set_heat_data: heat description'
	end if

!	---------------------------------------------------------
!	remember values and write to monitor
!	---------------------------------------------------------

	if( ihtype == 0 ) then
	  write(6,*) 'no heat file opened'
	else
	  write(6,*) 'heat file opened: ',ihtype
	  write(6,*) 'content: '
	  do i=1,nvar
	    call iff_get_var_description(id,i,string)
	    write(6,*) i,'    ',string
	  end do
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end subroutine meteo_set_heat_data

!*********************************************************************

	subroutine meteo_convert_heat_data(id,n
     +			,mettair,methum,metwbt)

	integer id
	integer n
	real mettair(n)
	real methum(n)
	real metwbt(n)

	logical bnowind

	bnowind = ihtype == 0
	
        if( bnowind ) then              !no wind
	  !nothing to be done
        else
	  call meteo_compute_wbt(ihtype,n,mettair,methum,metwbt)
	end if

	end subroutine meteo_convert_heat_data

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_set_array(n,val0,array)

! initializes 2D array

	implicit none

	integer n
	real val0
	real array(*)

	integer i

	do i=1,n
	  array(i) = val0
	end do

	end subroutine meteo_set_array

!*********************************************************************

	subroutine meteo_compute_wbt(mode,n,tav,rhv,wbv)

! computes wet bulb temperature

	implicit none

	include 'meteo.h'

	integer mode
	integer n
	real tav(n)
	real rhv(n)
	real wbv(n)

	integer i
	real db,rh,wb

	do i=1,n
	  db = tav(i)
	  rh = rhv(i)
	  wb = wbv(i)
	  if( mode .eq. 1 ) then		!val is humidity
	      call rh2twb(db,rh,wb)
	  else if( mode .eq. 2 ) then		!val is wet bulb
	      call twb2rh(db,wb,rh)
	  else
	      write(6,*) 'mode = ',mode
	      stop 'error stop meteo_convert_hum: mode'
	  end if
	  rhv(i) = rh
	  wbv(i) = wb
	end do

	end subroutine meteo_compute_wbt

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_get_heat_values(k,qs,ta,rh,wb,uw,cc,p)

! returns meteo parameters for one node
!
! pressure is returned in [mb]

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

	p = ppv(k)
	p = 0.01 * p					  !Pascal to mb

	end subroutine meteo_get_heat_values

!*********************************************************************

	subroutine meteo_get_solar_radiation0(k,qs)

	implicit none

	include 'meteo.h'

        integer k                       !node number
        real qs                         !solar radiation [W/m**2]

	qs = metrad(k)

	end subroutine meteo_get_solar_radiation0

!*********************************************************************

        subroutine get_wind(k,wx,wy)

! helper function -> return wind for node k

        implicit none

	include 'meteo.h'

        integer k
        real wx,wy

        wx = wxv(k)
        wy = wyv(k)

        end subroutine get_wind

!*********************************************************************

        subroutine convert_wind(s,d,u,v)

        implicit none

        real s,d,u,v

        real dir
        real pi,rad
        parameter(pi=3.14159,rad=pi/180.)

        dir = d
        dir = 90. - dir + 180.
        do while( dir .lt. 0. )
          dir = dir + 360.
        end do
        dir = mod(dir,360.)

        u = s*cos(rad*dir)
        v = s*sin(rad*dir)

        end subroutine convert_wind

!*********************************************************************

        subroutine get_drag(itdrag,wxy,dragco)

! computes drag coefficient

        implicit none

        integer itdrag          !type of formula
        real wxy                !wind speed
        real dragco             !computed drag coefficient

        if( itdrag .le. 0 ) then
          !nothing
        else if( itdrag .eq. 1 ) then   !Smith and Banke (1975)
          dragco = 0.001 * (0.63 + 0.066*wxy)
        else if( itdrag .eq. 2 ) then   !Large and Pond (1981)
          if ( wxy .gt. 11. ) then
            dragco = 0.001 * (0.49 + 0.066*wxy)
          else
            dragco = 0.001 * 1.2
          end if
        else
          write(6,*) 'erroneous value for itdrag = ',itdrag
          stop 'error stop get_drag: itdrag'
        end if

        end subroutine get_drag

!*********************************************************************

!================================================================
        end module meteo_forcing_module
!================================================================

