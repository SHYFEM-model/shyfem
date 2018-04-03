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
c 30.04.2015    ggu	ice integrated
c 04.05.2015    ggu	bug in ice eliminated
c 12.05.2015    ggu	introduced ia_icefree for icefree elements
c 08.01.2016    ggu	bug fix in meteo_convert_wind_data() - no wind bug
c 10.03.2016    ggu	check for pressure to be in reasonable bounds
c 23.07.2016    ivf	new heat formulation for iheat==8
c 09.09.2016    ggu	new variable ihtype to choose between rh, wbt, dpt
c 16.09.2016    ggu	allow for Pa and mbar in pressure
c 12.01.2017    ccf	bug fix in determining pressure units
c
c notes :
c
c info on file format can be found in subfemfile.f
c
c*********************************************************************
c
c DOCS  START   S_wind
c
c In this section the wind data can be given directly without
c the creation of an external file. Note, however, that
c a wind file specified in the |name| section takes precedence
c over this section. E.g., if both a section |wind| and a
c wind file in |name| is given, the wind data from the file is used.
c
c The format of the wind data in this section is the same as the
c format in the ASCII wind file, i.e., three columns, with
c the first specifying the time in seconds and the other two columns
c giving the wind data. The interpretation of the wind data
c depends on the value of |iwtype|. For more information please
c see the description of |iwtype| in section |para|.
c
c DOCS  END

c DOCS  START   P_wind
c
c DOCS  WIND            Wind parameters
c
c The next two parameters deal with the wind stress to be
c prescribed at the surface of the basin.
c
c The wind data can either be specified in an external file (ASCII
c or binary) or directly in the parameter file in section |wind|.
c The ASCII file or the wind section contain three columns, the first
c giving the time in seconds, and the others the components of
c the wind speed. Please see below how the last two columns are
c interpreted depending on the value of |iwtype|. For the format
c of the binary file please see the relative section.
c If both a wind file and section |wind| are given, data from the
c file is used.
c
c The wind stress is normally computed with the following formula
c \beq
c \tau^x = \rho_a c_D \vert u \vert u^x \quad
c \tau^y = \rho_a c_D \vert u \vert u^y
c \eeq
c where $\rho_a,\rho_0$ is the density of air and water respectively,
c $u$ the modulus of wind speed and $u^x,u^y$ the components of
c wind speed in $x,y$ direction. In this formulation $c_D$ is a
c dimensionless drag coefficient that varies between 1.5 \ten{-3} and
c 3.2 \ten{-3}. The wind speed is normally the wind speed measured
c at a height of 10 m.
c
c |iwtype|      The type of wind data given (default 1):
c               \begin{description}
c               \item[0] No wind data is processed
c               \item[1] The components of the wind is given in [m/s]
c               \item[2] The stress ($\tau^x,\tau^y$) is directly specified
c               \item[3] The wind is given in speed [m/s] and direction
c                        [degrees]. A direction of 0\degrees{} specifies
c                        a wind from the north, 90\degrees{} a wind
c                        from the east etc.
c               \item[4] As in 3 but the speed is given in knots
c               \end{description}
c |itdrag|	Formula to compute the drag coefficient. 
c		\begin{description}
c		\item[0] constant value given in |dragco|. 
c		\item[1] Smith and Banke (1975) formula
c		\item[2] Large and Pond (1981) formula
c		\item[3] Spatio/temporally varing in function of wave. Need
c		the coupling with WWMIII.
c		\item[4] Spatio/temporally varing in function of heat flux. 
c		Only for |iheat| = 6. 
c		\end{description}
c		(Default 0)
c |dragco|      Drag coefficient used in the above formula. 
c               Please note that in case
c               of |iwtype| = 2 this value is of no interest, since the
c               stress is specified directly. (Default 2.5E-3)
c |wsmax|       Maximum wind speed allowed in [m/s]. This is in order to avoid
c               errors if the wind data is given in a different format
c               from the one spwecified by |iwtype|. (Default 50)
c |wslim|	Limit maximum wind speed to this value [m/s]. This provides
c		an easy way to exclude strong wind gusts that might
c		blow up the simulation. Use with caution. 
c		(Default -1, no limitation)
c
c DOCS  END

!================================================================
        module meteo_forcing_module
!================================================================

	use intp_fem_file

	implicit none

	real, parameter :: pstd = 101325.	!standard pressure
	real, parameter :: tkelvin = 273.15	!0C in Kelvin
	real, parameter :: dstd = 2.5e-3	!standard drag coefficient
	real, parameter :: nmile = 1852.	!nautical mile in m

	real, parameter :: amice = 1.	!use momentum reduction due to ice
					!1: use  0: do not reduce momentum

	integer, save :: idwind,idheat,idrain,idice

	integer, parameter :: nfreq = 0		!debug output
	double precision, save, private :: da_out(4) = 0

	integer, save :: iwtype,itdrag
	integer, save :: irtype
	integer, save :: ihtype
	integer, save :: ictype
	integer, save :: ia_icefree		!area type which is ice free
	real, save :: wsmax,wslim,dragco,roluft,rowass
	real, save :: pfact = 1.
	real, save :: wfact = 1.
	real, save :: sfact = 1.

	logical, save :: has_pressure = .false.

	logical, save, private :: bdebug = .true.

	integer, save, private :: icall = 0

	character*80, save :: wxss = 'wind stress - x [N/m**2]'
	character*80, save :: wyss = 'wind stress - y [N/m**2]'
	character*80, save :: wxms = 'wind velocity - x [m/s]'
	character*80, save :: wyms = 'wind velocity - y [m/s]'
	character*80, save :: wsms = 'wind speed [m/s]'
	character*80, save :: wskn = 'wind speed [knots]'
	character*80, save :: wdir = 'wind direction [deg]'

	character*80, save :: papa = 'pressure (atmospheric) [Pa]'
	character*80, save :: pamb = 'pressure (atmospheric) [mbar]'

	character*80, save :: rain = 'rain [mm/day]'
	character*80, save :: ice  = 'ice cover [0-1]'

        character*80, save :: srad = 'solar radiation [W/m**2]'
        character*80, save :: tair = 'air temperature [C]'
        character*80, save :: rhum = 'humidity (relative) [%]'
        character*80, save :: shum = 'humidity (specific)'
        character*80, save :: ccov = 'cloud cover [0-1]'
        character*80, save :: wbtm = 'wet bulb temperature [C]'
        character*80, save :: dewp = 'dew point temperature [C]'    

!================================================================
        contains
!================================================================

!*********************************************************************

	subroutine meteo_forcing_fem

! administers meteo files (still to be cleaned)

	use mod_meteo
	use mod_depth
	use levels
	use basin, only : nkn,nel,ngr,mbw

	character*60 windfile,heatfile,rainfile,icefile
	character*4 what

	integer nvar,lmax
	integer nintp
	integer modehum
	integer nvarm,nid,nlev
	integer i
	!integer it0
	double precision dtime0,dtime
	real flag
	real val0,val1,val2
	real metaux(nkn)

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
	  call getfnm('ice',icefile)

	  call get_first_dtime(dtime0)

	  write(6,'(a)') 'opening wind file...'
          nvar = 0      !not sure if 2 or 3
          call iff_get_file_nvar(windfile,nvar)
          if( nvar <= 0 ) nvar = 3      !if no file fake 3
          nintp = 2
          what = 'wind'
          vconst = (/ 0., 0., pstd, 0. /)
	  call iff_init(dtime0,windfile,nvar,nkn,0,nintp
     +				,nodes,vconst,idwind)
	  call iff_set_description(idwind,0,'meteo wind')

	  call meteo_set_wind_data(idwind,nvar)

	  write(6,'(a)') 'opening ice file...'
	  nvar = 1
	  nintp = 2
	  what = 'ice'
	  vconst = (/ 0., 0., 0., 0. /)
	  call iff_init(dtime0,icefile,nvar,nkn,0,nintp
     +				,nodes,vconst,idice)
	  call iff_set_description(idice,0,'meteo ice')
	  call iff_flag_ok(idice)	!we do not need all icd data

	  call meteo_set_ice_data(idice,nvar)

	  write(6,'(a)') 'opening rain file...'
	  nvar = 1
	  nintp = 2
	  what = 'rain'
	  vconst = (/ 0., 0., 0., 0. /)
	  call iff_init(dtime0,rainfile,nvar,nkn,0,nintp
     +				,nodes,vconst,idrain)
	  call iff_set_description(idrain,0,'meteo rain')

	  call meteo_set_rain_data(idrain,nvar)

	  write(6,'(a)') 'opening heat flux file...'
	  nvar = 4
	  nintp = 2
	  what = 'heat'
	  vconst = (/ 0., 0., 50., 0. /)
	  call iff_init(dtime0,heatfile,nvar,nkn,0,nintp
     +				,nodes,vconst,idheat)
	  call iff_set_description(idheat,0,'meteo heat')

	  call meteo_set_heat_data(idheat,nvar)

	end if

!------------------------------------------------------------------
! end of initialization
!------------------------------------------------------------------

	icall = icall + 1

!------------------------------------------------------------------
! time interpolation
!------------------------------------------------------------------

	call get_act_dtime(dtime)
	lmax = 1

	if( .not. iff_is_constant(idice) .or. icall == 1 ) then
	  call iff_read_and_interpolate(idice,dtime)
	  metice = 0.	!assume no ice cover in areas not covered by data
	  call iff_time_interpolate(idice,dtime,1,nkn,lmax,metice)
	end if
	!write(6,*) (metice(i),i=1,nkn,50)
	!write(6,*) 'constant: ',iff_is_constant(idice)

	if( .not. iff_is_constant(idwind) .or. icall == 1 ) then
	  call iff_read_and_interpolate(idwind,dtime)
	  call iff_time_interpolate(idwind,dtime,1,nkn,lmax,wxv)
	  call iff_time_interpolate(idwind,dtime,2,nkn,lmax,wyv)
	  if( iff_get_nvar(idwind) == 3 ) then
	    call iff_time_interpolate(idwind,dtime,3,nkn,lmax,ppv)
	  else
	    ppv = pstd
	  end if
	end if

	!call iff_print_info(idwind,0,.true.)

        if( .not. iff_is_constant(idrain) .or. icall == 1 ) then
          call iff_read_and_interpolate(idrain,dtime)
          call iff_time_interpolate(idrain,dtime,1,nkn,lmax,metrain)
        end if

        if( .not. iff_is_constant(idheat) .or. icall == 1 ) then
          call iff_read_and_interpolate(idheat,dtime)
          call iff_time_interpolate(idheat,dtime,1,nkn,lmax,metrad)
          call iff_time_interpolate(idheat,dtime,2,nkn,lmax,mettair)
          call iff_time_interpolate(idheat,dtime,3,nkn,lmax,metaux)
          call iff_time_interpolate(idheat,dtime,4,nkn,lmax,metcc)
        end if

!------------------------------------------------------------------
! extra treatment of data
!------------------------------------------------------------------

	if( .not. iff_is_constant(idice) .or. icall == 1 ) then
	  call meteo_convert_ice_data(idice,nkn,metice)
	end if

	if( .not. iff_is_constant(idwind) .or. icall == 1 ) then
	  call meteo_convert_wind_data(idwind,nkn,wxv,wyv
     +			,windcd,tauxnv,tauynv,metws,ppv,metice)
	end if

!	write(166,*) (wxv(i),wyv(i),windcd(i),tauxnv(i),tauynv(i)
!     +			,i=1,nkn,nkn/20)

        if( .not. iff_is_constant(idheat) .or. icall == 1 ) then
          call meteo_convert_heat_data(idheat,nkn
     +                       ,metaux,mettair,ppv,methum)
        end if

	if( .not. iff_is_constant(idrain) .or. icall == 1 ) then
	  call meteo_convert_rain_data(idrain,nkn,metrain)
	end if

!------------------------------------------------------------------
! debug output
!------------------------------------------------------------------

	if( nfreq .gt. 0 .and. mod(icall,nfreq) .eq. 0 ) then
	  nvarm = 4		!total number of vars that are written
	  nlev = 1
	  call scalar_output_file(da_out,'meteo',nvarm,20,nlev,ppv)
	  call scalar_output_file(da_out,'meteo',nvarm,28,nlev,metws)
	  call scalar_output_file(da_out,'meteo',nvarm,23,nlev,mettair)
	  call scalar_output_file(da_out,'meteo',nvarm,85,nlev,metice)
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

	integer il
	character*80 string,string1,string2
	character*10 dir,unit

	logical string_is_this_short
	real, external :: getpar

!	---------------------------------------------------------
!	check nvar and get parameters
!	---------------------------------------------------------

	if( iff_has_file(id) ) then
	 if( nvar /= 2 .and. nvar /= 3 ) then
	  write(6,*) 'no support for nvar = ',nvar
	  stop 'error stop meteo_set_wind_data: wind'
	 end if
	end if

        iwtype = nint(getpar('iwtype'))
        itdrag = nint(getpar('itdrag'))
        wsmax = getpar('wsmax')
        wslim = getpar('wslim')
        dragco = getpar('dragco')
        roluft = getpar('roluft')
        rowass = getpar('rowass')

	if( dragco .lt. 0 ) dragco = dstd

!	---------------------------------------------------------
!	handle wind
!	---------------------------------------------------------

	call iff_get_var_description(id,1,string1)
	call iff_get_var_description(id,2,string2)

	if( .not. iff_has_file(id) ) then

	  iwtype = 0

	else if( string1 == ' ' ) then	!TS file or constant

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

          if( string_is_this_short('wind',string1) .and.
     +		string_is_this_short('wind',string2) ) then
	    iwtype = 1
          else if( string_is_this_short('wstress',string1) .and.
     +		string_is_this_short('wstress',string2) ) then
	    iwtype = 2
          else if( string_is_this_short('windspeed',string1) .and.
     +		string_is_this_short('winddir',string2) ) then
	    iwtype = 3
	  else
	    write(6,*) 'description string for wind not recognized: '
	    write(6,*) trim(string1)
	    write(6,*) trim(string2)
	    stop 'error stop meteo_set_wind_data: wind description'
	  end if

	  if( iwtype == 3 ) then
	    call string_direction_and_unit(string1,dir,unit)
	    if( unit == 'knots' ) iwtype = 4
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
	  has_pressure = .true.
	  call iff_get_var_description(id,3,string)
          if( string_is_this_short('airp',string) ) then
	    call string_direction_and_unit(string1,dir,unit)
	    if( unit == 'Pa' ) then
	      pfact = 1.
	    else if( unit == 'mbar' ) then
	      pfact = 100.
	    else if( unit == ' ' ) then		!must determine later
	      pfact = 0.			!check later
	    end if
	  else if( string == ' ' ) then
	    pfact = 1.
	  else
	    write(6,*) 'description for pressure not recognized: '
	    write(6,*) trim(string)
	    stop 'error stop meteo_set_wind_data: press description'
	  end if
	  if( bdebug ) then
	    write(6,*) 'pressure initialized: ',pfact,trim(string)
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
	  write(6,*) ' 1    ',trim(string1)
	  write(6,*) ' 2    ',trim(string2)
	  if( nvar == 3 ) then
	    call iff_get_var_description(id,3,string)
	    write(6,*) ' 3    ',trim(string)
	  end if
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end subroutine meteo_set_wind_data

!*********************************************************************

	subroutine meteo_convert_wind_data(id,n,wx,wy,cdv,tx,ty,ws
     +						,pp,cice)

	integer id
	integer n
	real wx(n),wy(n)
	real cdv(n)
	real tx(n),ty(n)
	real ws(n)
	real pp(n)
	real cice(n)

	logical bnowind,bstress,bspeed
	integer k
	real cd,wxymax,txy,wspeed,wdir,fact,fice,aice
	real pmin,pmax
	character*80 string

	bnowind = iwtype == 0
	bstress = iwtype == 2
	bspeed = iwtype > 2
	cd = dragco
	wxymax = 0.
	aice = amice       !ice cover for momentum: 1: use  0: do not use
	
!	---------------------------------------------------------
!	convert wind
!	---------------------------------------------------------

        if( bnowind ) then              !no wind
	  ws = 0.
	  wx = 0.
	  wy = 0.
	  tx = 0.
	  ty = 0.
	  pp = pstd
        else if( bstress ) then         !data is stress -> normalize it
          if( cd .le. 0 ) cd = dstd
          do k=1,n
	    fice = 1. - aice*cice(k)
            tx(k) = fice * wfact * wx(k)
            ty(k) = fice * wfact * wy(k)
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
	    fice = 1. - aice*cice(k)
	    !if( k .eq. 1000 ) write(6,*) 'ice: ',k,fice
            wspeed = ws(k)
            wxymax = max(wxymax,wspeed)
	    cd = cdv(k)
            if( itdrag .gt. 0 .and. itdrag .le. 2 ) then
		call get_drag(itdrag,wspeed,cd)
	    end if
            tx(k) = fice * wfact * cd * wspeed * wx(k)
            ty(k) = fice * wfact * cd * wspeed * wy(k)
          end do
        end if

!	---------------------------------------------------------
!	limit wind
!	---------------------------------------------------------

	if( wslim > 0 .and. wxymax > wslim ) then !artificially limit wind speed
          do k=1,n
	    fice = 1. - aice*cice(k)
            wspeed = ws(k)
	    if( wspeed <= wslim ) cycle
	    ws(k) = wslim
	    fact = wslim/wspeed
	    wx(k) = fact*wx(k)
	    wy(k) = fact*wy(k)
            wspeed = ws(k)
	    cd = cdv(k)
            if( itdrag .gt. 0 .and. itdrag .le. 2 ) then
		call get_drag(itdrag,wspeed,cd)
	    end if
            tx(k) = fice * wfact * cd * wspeed * wx(k)
            ty(k) = fice * wfact * cd * wspeed * wy(k)
          end do
	  wxymax = wslim
	end if

!	---------------------------------------------------------
!	check wind speed
!	---------------------------------------------------------

        if( wxymax .gt. wsmax ) then
          write(6,*) 'maximum wind speed: ',wxymax
          write(6,*) 'maximum allowed wind speed: ',wsmax
          write(6,*) 'Are you sure the wind is in the correct format?'
          write(6,*) 'If no, please set iwtype, else increase wsmax.'
          stop 'error stop meteo_convert_wind_data: wind speed too high'
        end if

!	---------------------------------------------------------
!	convert pressure
!	---------------------------------------------------------

	if( .not. has_pressure ) return

	call iff_get_var_description(id,3,string)

	if( string == ' ' .or. pfact == 0. ) then !only if not yet determined
	  pmin = minval(pp)
	  pmax = maxval(pp)
	  if( pmin /= 0 .and. pmax /= 0. ) then
	    if( pmin > 85000 .and. pmax < 110000 ) then
	      pfact = 1.
	      string = papa
	    else if( pmin > 850 .and. pmax < 1100 ) then
	      pfact = 100.
	      string = pamb
	    else
	      pfact = 1.
	      string = 'unknown'
	    end if
	    !call iff_set_var_description(id,3,string)
	  end if
	end if

	if( pfact /= 1. ) pp = pfact * pp

	pmin = minval(pp)
	pmax = maxval(pp)

	if( pmin /= 0 .or. pmax /= 0. ) then
	  if( pmin < 85000 .or. pmax > 110000 ) then
	    write(6,*) 'pmin,pmax: ',pmin,pmax
	    write(6,*) 'pressure values out of range'
	    write(6,*) 'pressure should be given in Pascal'
	    stop 'error stop meteo_convert_wind_data: pressure'
	  end if
	end if

!	---------------------------------------------------------
!	end of routine wind speed
!	---------------------------------------------------------

	end subroutine meteo_convert_wind_data

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_set_rain_data(id,nvar)

	integer id
	integer nvar

	character*60 string

        logical string_is_this_short

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
          if(string_is_this_short('rain',string)) then
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

	subroutine meteo_convert_rain_data(id,n,r)

c convert rain from mm/day to m/s

	integer id
	integer n
	real r(n)

	integer i
        real zconv
        parameter( zconv = 1.e-3 / 86400. )

	do i=1,n
	  r(i) = r(i) * zconv
	end do

	end subroutine meteo_convert_rain_data

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_set_ice_data(id,nvar)

	integer id
	integer nvar

	character*60 string

        logical string_is_this_short
	real getpar

!	---------------------------------------------------------
!	check nvar and get parameters
!	---------------------------------------------------------

	if( nvar /= 1 ) then
	  write(6,*) 'no support for nvar = ',nvar
	  stop 'error stop meteo_set_ice_data: ice'
	end if

!	---------------------------------------------------------
!	handle ice
!	---------------------------------------------------------

	call iff_get_var_description(id,1,string)
	string = adjustl(string)

	if( string == ' ' ) then	!TS file or constant
	  ictype = 0
	  if( iff_has_file(id) ) ictype = 1

	  if( ictype == 1 ) then
	    call iff_set_var_description(id,1,ice)
	  end if
	else
          if(string_is_this_short('ice',string)) then
	    ictype = 1
	  else
	    write(6,*) 'description string for ice not recognized: '
	    write(6,*) string
	    write(6,*) 'expecting: ',trim(ice)
	    stop 'error stop meteo_set_ice_data: ice description'
	  end if
	end if

!	---------------------------------------------------------
!	handle ice free areas
!	---------------------------------------------------------

	ia_icefree = -1		!this does not change ice cover

	if( ictype /= 0 ) then	!ice file has been opened
	  ia_icefree = nint(getpar('iaicef'))
	  if( ia_icefree == -99 ) then
	    write(6,*) 'ice file has been opened'
	    write(6,*) 'but parameter iaicef has not been set'
	    write(6,*) 'please set iaicef to the area code that'
	    write(6,*) 'indicates ice free conditions'
	    write(6,*) 'if no such areas exist please set iaicef=-1'
	    stop 'error stop meteo_set_ice_data: iaicef'
	  end if
	end if

!	---------------------------------------------------------
!	remember values and write to monitor
!	---------------------------------------------------------

	if( ictype == 0 ) then
	  write(6,*) 'no ice file opened'
	else
	  write(6,*) 'ice file opened: ',ictype
	  call iff_get_var_description(id,1,string)
	  write(6,*) 'content: '
	  write(6,*) ' 1    ',string
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end subroutine meteo_set_ice_data

!*********************************************************************

	subroutine meteo_convert_ice_data(id,n,r)

c convert ice data (delete ice in ice free areas, compute statistics)

	use evgeom
	use basin

	integer id
	integer n
	real r(n)	!ice concentration

	integer k,ie,ii,ia,nflag
	real rarea,rnodes,rorig
	double precision dacu,dice,darea,area
	character*20 aline

	integer, save :: ninfo = 0
	real, parameter :: flag = -999.

	if( ninfo == 0 ) call getinfo(ninfo)

	dacu = 0.
	do k=1,n
	  dacu = dacu + r(k)
	end do
	rorig = dacu / n

	nflag = 0
	dice = 0.
	darea = 0.

	do k=1,n
	  if( r(k) == flag ) then
	    nflag = nflag + 1
	    r(k) = 0
	  end if
	end do

	do ie=1,nel
	  ia = iarv(ie)
	  area = 4. * ev(10,ie)
	  if( ia == ia_icefree ) then
	    do ii=1,3
	      k = nen3v(ii,ie)
	      r(k) = 0.
	    end do
	  else
	    do ii=1,3
	      k = nen3v(ii,ie)
	      dice = dice + area*r(k)
	      darea = darea + area
	    end do
	  end if
	end do
	dice = dice / darea

	dacu = 0.
	do k=1,n
	  dacu = dacu + r(k)
	end do
	rnodes = dacu / n

	rarea = dice
	call get_act_timeline(aline)
	write(ninfo,*) 'ice: ',aline,rarea,rnodes,nflag

	end subroutine meteo_convert_ice_data

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_set_heat_data(id,nvar)

	use shyfem_strings

	integer id
	integer nvar

	character*80 string,strings(4)
	integer i,ierr
	character*80 vapor,vshort

	logical string_is_this_short
        real getpar  

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
	call adjust_humidity_string(strings(3))		!FIXME

        ihtype = nint(getpar('ihtype'))  
	if( ihtype == 1 ) then
	  vapor = rhum
	else if( ihtype == 2 ) then
	  vapor = wbtm
	else if( ihtype == 3 ) then
	  vapor = dewp
	else if( ihtype == 4 ) then
	  vapor = shum
	else
	  write(6,*) 'ihtype = ',ihtype
	  stop 'error stop meteo_set_heat_data: erroneous ihtype'
	end if

	call strings_get_short_name(vapor,vshort)

	if( strings(1) == ' ' ) then	!TS file or constant
	  ierr = 0
	  if( .not. iff_has_file(id) ) ihtype = 0	!no heat

          if( ihtype > 0 ) then
            call iff_set_var_description(id,1,srad)
            call iff_set_var_description(id,2,tair)
            call iff_set_var_description(id,3,vapor)
            call iff_set_var_description(id,4,ccov)
          end if
        else
	  ierr = 0
          if(.not.string_is_this_short('srad',strings(1))) then
            ierr = 1
          else if(.not.string_is_this_short('airt',strings(2))) then
            ierr = 2
	  else if(.not.string_is_this_short(vshort,strings(3))) then
            ierr = 3
          else if(.not.string_is_this_short('cc',strings(4))) then
            ierr = 4
          end if
        end if

	if( ierr /= 0 ) then
	  write(6,*) 'description string for heat not recognized: '
	  write(6,*) 'possible number of string with error: ',ierr
	  if( ierr == 3 ) then
	    write(6,*) 'ihtype = ',ihtype,'  ',trim(vapor)
	    write(6,*) 'description in file: ',trim(strings(3))
	  end if
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
	    if( i == 3 ) call adjust_humidity_string(string)		!FIXME
	    write(6,*) i,'    ',trim(string)
	  end do
	end if

!	---------------------------------------------------------
!	end of routine
!	---------------------------------------------------------

	end subroutine meteo_set_heat_data

!*********************************************************************
	subroutine adjust_humidity_string(string)
	implicit none
	character*(*) string
	if( string == 'humidity [%]' ) string = rhum
	end
!*********************************************************************

        subroutine meteo_convert_heat_data(id,n
     +                  ,metaux,mettair,ppv,methum)

	integer id
	integer n
	real metaux(n)		!this is the vapor information read
	real mettair(n)
	real ppv(n)
	real methum(n)		!return

	logical bnoheat

	bnoheat = ihtype == 0
	
        if( bnoheat ) then              !no heat
	  !nothing to be done
        else
	  call meteo_convert_temperature(n,mettair)
	  call meteo_convert_vapor(ihtype,n
     +			,metaux,mettair,ppv,methum)
	end if

	end subroutine meteo_convert_heat_data

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_convert_temperature(n,tav)

! computes temperture to Celsius if in Kelvin

	use mod_meteo

	implicit none

	integer n
	real tav(n)	!dry air temperature

	where( tav > 200. ) tav = tav - tkelvin

	end subroutine meteo_convert_temperature

!*********************************************************************

	subroutine meteo_convert_vapor(mode,n,aux,tav,pav,rhv)

! computes relative humidity from other vapor values

	use mod_meteo

	implicit none

	integer mode
	integer n
	real aux(n)	!value read
	real tav(n)	!dry air temperature
	real pav(n)	!atmospheric pressure
	real rhv(n)	!relative humidity (return)

	integer i
	real ta,pp,rh,wb,dp,sh,val

	if( mode < 1 .or. mode > 4 ) then
	      write(6,*) 'mode = ',mode
	      stop 'error stop meteo_convert_hum: mode'
	end if

	do i=1,n
	  val = aux(i)
	  ta = tav(i)
	  pp = pav(i)/100.			!pressure in mbar
	  if( mode .eq. 1 ) then		!val is humidity
	      rh = val
	  else if( mode .eq. 2 ) then		!val is wet bulb
	    call wb2rh(ta,pp,val,rh)
	  else if( mode .eq. 3 ) then		!val is dew point
	    call dp2rh(ta,pp,val,rh)
	  else if( mode .eq. 4 ) then		!val is specific humidity
	    call sh2rh(ta,pp,val,rh)
	  end if
	  rhv(i) = rh
	end do

	end subroutine meteo_convert_vapor

!*********************************************************************
!*********************************************************************
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

!================================================================
        end module meteo_forcing_module
!================================================================

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

	subroutine meteo_get_heat_values(k,qs,ta,rh,twb,uw,cc,p)

! returns meteo parameters for one node
!
! pressure is returned in [mb]

	use mod_meteo

	implicit none

        integer k                       !node number
        real qs                         !solar radiation [W/m**2]
        real ta                         !air temperature [Celsius]
        real rh                         !relative humidity [%, 0-100]
        real twb                        !wet bulb temperature
        real uw                         !wind speed [m/s]
        real cc                         !cloud cover [0-1]
        real p                          !atmospheric pressure [mbar, hPa]

	qs = metrad(k)
	ta = mettair(k)
	rh = methum(k)
	uw = metws(k)
	cc = metcc(k)

	cc = max(0.,cc)
	cc = min(1.,cc)
	rh = max(0.,rh)
	rh = min(100.,rh)

	p = ppv(k)
	p = 0.01 * p					  !Pascal to mb

	call rh2wb(ta,p,rh,twb)

	end subroutine meteo_get_heat_values

!*********************************************************************

        subroutine meteo_get_heat_extra(k,dp,uuw,vvw)
 
! returns iextra meteo parameters for one node (iheat == 8)

        use mod_meteo

        implicit none

        integer k                       !node number
        real dp                         !dew point temperature [Celsius]  
        real uuw                        !u-component wind speed [m/s]
        real vvw                        !v-component wind speed [m/s]

        !dp  = metdew(k)  
	call rh2dp(mettair(k),ppv(k),methum(k),dp)
        uuw = wxv(k)
        vvw = wyv(k)

        end subroutine meteo_get_heat_extra

!*********************************************************************

	subroutine get_pe_values(k,r,e,eeff)

c returns precipitation and evaporation values
c
c eeff is evaporation used in model, if ievap==0 => eeff==0.

	use mod_meteo

	implicit none

	integer k
	real r			!rain [m/s]
	real e			!evaporation [m/s]
	real eeff		!effective evaporation [m/s] - used in model

	integer ievap
	save ievap
	data ievap /-1/

	real getpar

        if( ievap .eq. -1 ) ievap = nint(getpar('ievap'))

	r = metrain(k)
	e = evapv(k)
	eeff = e*ievap

	end subroutine get_pe_values

!*********************************************************************

	subroutine set_evap(k,e)

c sets evaporation

	use mod_meteo

	implicit none

	integer k
	real e			!evaporation [m/s]

	evapv(k) = e

	end subroutine set_evap

!*********************************************************************

	subroutine meteo_get_solar_radiation(k,qs)

	use mod_meteo

	implicit none

        integer k                       !node number
        real qs                         !solar radiation [W/m**2]

	qs = metrad(k)

	end subroutine meteo_get_solar_radiation

!*********************************************************************

        subroutine get_wind(k,wx,wy)

! helper function -> return wind for node k

	use mod_meteo

        implicit none

        integer k
        real wx,wy

        wx = wxv(k)
        wy = wyv(k)

        end subroutine get_wind

!*********************************************************************

        subroutine meteo_set_matrix(qs,ta,rh,uw,cc)

c interpolates files spatially - to be deleted

	use mod_meteo
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        real qs,ta,rh,wb,uw,cc

        integer k

        do k=1,nkn
          metrad(k) = qs
          mettair(k) = ta
          methum(k) = rh
          metws(k) = uw
          metcc(k) = cc
        end do

        end

!*********************************************************************

