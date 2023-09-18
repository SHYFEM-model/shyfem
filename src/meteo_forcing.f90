!
! $Id: submeteo.f,v 1.7 2010-02-26 17:35:06 georg Exp $
!
! handle meteo files with new fem format
!
! revision log :
!
! 10.03.2009    ggu     finished coding
! 24.03.2009    ggu     use metrain instead rqdsv
! 07.05.2009    ggu     new call to init_coords()
! 18.06.2009    ggu&dbf bug fix -> wind speed not interpolated on metws
! 23.02.2010    ggu	call to wstress changed (new wxv,wyv)
! 26.01.2011    ggu	write wind field to debug output (iumetw)
! 05.02.2011    ggu	changed order of records in dfile (actual is first)
! 16.02.2011    ggu	pass idata to files, use geo info from files
! 18.11.2011    ggu	deleted projection code from subroutines
! 10.02.2012    ggu	limit cc and rh to acceptable range
! 16.02.2012    ggu	new routine meteo_get_solar_radiation()
! 22.02.2012    ggu	new routines for regular and ts reading of meteo
! 23.02.2012    ggu&ccf	bug fix meteo_copy_to_old and meteo_interpolate_in_time
! 20.05.2014    ggu	new routines for new file format
! 30.04.2015    ggu	ice integrated
! 04.05.2015    ggu	bug in ice eliminated
! 12.05.2015    ggu	introduced ia_icefree for icefree elements
! 08.01.2016    ggu	bug fix in meteo_convert_wind_data() - no wind bug
! 10.03.2016    ggu	check for pressure to be in reasonable bounds
! 23.07.2016    ivf	new heat formulation for iheat==8
! 09.09.2016    ggu	new variable ihtype to choose between rh, wbt, dpt
! 16.09.2016    ggu	allow for Pa and mbar in pressure
! 12.01.2017    ccf	bug fix in determining pressure units
!
! notes :
!
! info on file format can be found in subfemfile.f
!
!*********************************************************************
!
! DOCS  START   S_wind
!
! In this section the wind data can be given directly without
! the creation of an external file. Note, however, that
! a wind file specified in the |name| section takes precedence
! over this section. E.g., if both a section |wind| and a
! wind file in |name| is given, the wind data from the file is used.
!
! The format of the wind data in this section is the same as the
! format in the ASCII wind file, i.e., three columns, with
! the first specifying the time in seconds and the other two columns
! giving the wind data. The interpretation of the wind data
! depends on the value of |iwtype|. For more information please
! see the description of |iwtype| in section |para|.
!
! DOCS  END

! DOCS  START   P_wind
!
! DOCS  WIND            Wind parameters
!
! The next two parameters deal with the wind stress to be
! prescribed at the surface of the basin.
!
! The wind data can either be specified in an external file (ASCII
! or binary) or directly in the parameter file in section |wind|.
! The ASCII file or the wind section contain three columns, the first
! giving the time in seconds, and the others the components of
! the wind speed. Please see below how the last two columns are
! interpreted depending on the value of |iwtype|. For the format
! of the binary file please see the relative section.
! If both a wind file and section |wind| are given, data from the
! file is used.
!
! The wind stress is normally computed with the following formula
! \beq
! \tau^x = \rho_a c_D \vert u \vert u^x \quad
! \tau^y = \rho_a c_D \vert u \vert u^y
! \eeq
! where $\rho_a,\rho_0$ is the density of air and water respectively,
! $u$ the modulus of wind speed and $u^x,u^y$ the components of
! wind speed in $x,y$ direction. In this formulation $c_D$ is a
! dimensionless drag coefficient that varies between 1.5 \ten{-3} and
! 3.2 \ten{-3}. The wind speed is normally the wind speed measured
! at a height of 10 m.
!
! |iwtype|      The type of wind data given (default 1):
!               \begin{description}
!               \item[0] No wind data is processed
!               \item[1] The components of the wind is given in [m/s]
!               \item[2] The stress ($\tau^x,\tau^y$) is directly specified
!               \item[3] The wind is given in speed [m/s] and direction
!                        [degrees]. A direction of 0\degrees{} specifies
!                        a wind from the north, 90\degrees{} a wind
!                        from the east etc.
!               \item[4] As in 3 but the speed is given in knots
!               \end{description}
! |itdrag|	Formula to compute the drag coefficient. 
!		\begin{description}
!		\item[0] constant value given in |dragco|. 
!		\item[1] Smith and Banke (1975) formula
!		\item[2] Large and Pond (1981) formula
!		\item[3] Spatio/temporally varing in function of wave. Need
!		the coupling with WWMIII.
!		\item[4] Spatio/temporally varing in function of heat flux. 
!		Only for |iheat| = 6. 
!		\end{description}
!		(Default 0)
! |dragco|      Drag coefficient used in the above formula. 
!               Please note that in case
!               of |iwtype| = 2 this value is of no interest, since the
!               stress is specified directly. (Default 2.5E-3)
! |wsmax|       Maximum wind speed allowed in [m/s]. This is in order to avoid
!               errors if the wind data is given in a different format
!               from the one spwecified by |iwtype|. (Default 50)
! |wslim|	Limit maximum wind speed to this value [m/s]. This provides
!		an easy way to exclude strong wind gusts that might
!		blow up the simulation. Use with caution. 
!		(Default -1, no limitation)
!
! DOCS  END

!================================================================
        module meteo_forcing
!================================================================

        use intp_fem_file

	implicit none

	double precision, parameter :: pstd = 101325.d0     !standard pressure
	double precision, parameter :: dstd = 2.5e-3        !standard drag coefficient
	double precision, parameter :: nmile = 1852.d0      !nautical mile in m

	integer, save :: idwind,idheat,idrain,idice

	integer, save :: nfreq = 0		!debug output
	integer, save :: iumet = 0

	integer, save :: iwtype,itdrag
	integer, save :: irtype
	integer, save :: ihtype
	integer, save :: ictype
	integer, save :: ia_icefree		!area type which is ice free
	double precision, save :: wsmax,wslim,dragco,roluft,rowass
        double precision, save :: pfact = 1.d0
        double precision, save :: wfact = 1.d0
        double precision, save :: sfact = 1.d0

	logical, save :: has_pressure = .false.

	logical, save, private :: bdebug = .true.

	integer, save, private :: icall = 0

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
	character*80, save :: ice = 'ice cover [0-1]'

        character*80, save :: srad = 'solar radiation [W/m**2]'
        character*80, save :: tair = 'air temperature [C]'
        character*80, save :: rhum = 'humidity [%]'
        character*80, save :: ccov = 'cloud cover [0-1]'
        character*80, save :: wbtm = 'wet bulb temperature [C]'
        character*80, save :: dewp = 'dew point temperature [C]'    

!================================================================
        contains
!================================================================

!*********************************************************************

        subroutine meteo_forcing_fem

! administers meteo files (still to be cleaned)

        use meteo
        use depth
        use levels
        use basin, only : nkn,nel,ngr,mbw
        use shympi
        use para
        use conz_util

        include 'femtime.h'

        character*60 windfile,heatfile,rainfile,icefile
        character*4 what

        integer nvar,lmax
        integer nintp
        integer modehum
        integer nvarm,nid,nlev
        integer i
	!integer it0
        double precision dtime0,dtime
        double precision flag
        double precision val0,val1,val2
        double precision metaux(nkn)

        double precision vconst(4)
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

	  !it0 = itanf
          dtime0 = itanf

          write(6,'(a)') 'opening wind file...'
          nvar = 0      !not sure if 2 or 3
          call iff_get_file_nvar(windfile,nvar)
          if( nvar <= 0 ) nvar = 3      !if no file fake 3
          nintp = 2
          what = 'wind'
          vconst = (/ 0.d0, 0.d0, pstd, 0.d0 /)
          if(bmpi) then
            call iff_init_mpi(dtime0,windfile,nvar,nkn,0,0,nintp,nodes,vconst,idwind,0)
          else
            call iff_init(dtime0,windfile,nvar,nkn,0,nintp,nodes,vconst,idwind)
          end if
          call iff_set_description(idwind,0,'meteo wind')

          call meteo_set_wind_data(idwind,nvar)

          write(6,'(a)') 'opening ice file...'
          nvar = 1
          nintp = 2
          what = 'ice'
          vconst = (/ 0., 0., 0., 0. /)
          if(bmpi) then
            call iff_init_mpi(dtime0,icefile,nvar,nkn,0,0,nintp,nodes,vconst,idice,0)
          else
            call iff_init(dtime0,icefile,nvar,nkn,0,nintp,nodes,vconst,idice)
          end if
          call iff_set_description(idice,0,'meteo ice')

          call meteo_set_ice_data(idice,nvar)

          write(6,'(a)') 'opening rain file...'
          nvar = 1
          nintp = 2
          what = 'rain'
          vconst = (/ 0., 0., 0., 0. /)
          if(bmpi) then
            call iff_init_mpi(dtime0,rainfile,nvar,nkn,0,0,nintp,nodes,vconst,idrain,0)
          else
            call iff_init(dtime0,rainfile,nvar,nkn,0,nintp,nodes,vconst,idrain)
          end if
          call iff_set_description(idrain,0,'meteo rain')

          call meteo_set_rain_data(idrain,nvar)

          write(6,'(a)') 'opening heat flux file...'
          nvar = 4
          nintp = 2
          what = 'heat'
          vconst = (/ 0., 0., 50., 0. /)
          if(bmpi) then
            call iff_init_mpi(dtime0,heatfile,nvar,nkn,0,0,nintp,nodes,vconst,idheat,0)
          else
            call iff_init(dtime0,heatfile,nvar,nkn,0,nintp,nodes,vconst,idheat)
          end if
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

        dtime = it
        lmax = 1

        if( .not. iff_is_constant(idice) .or. icall == 1 ) then
          call iff_read_and_interpolate(idice,dtime)
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
            call meteo_set_array(nkn,pstd,ppv)
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
          call meteo_convert_wind_data(idwind,nkn,wxv,wyv,windcd,tauxnv,tauynv,metws,ppv,metice)
        end if

!	write(166,*) (wxv(i),wyv(i),windcd(i),tauxnv(i),tauynv(i)
!     +			,i=1,nkn,nkn/20)

        if( .not. iff_is_constant(idheat) .or. icall == 1 ) then
          call meteo_convert_heat_data(idheat,nkn,metaux,mettair,methum,metwbt,metdew)
        end if

        if( .not. iff_is_constant(idrain) .or. icall == 1 ) then
          call meteo_convert_rain_data(idrain,nkn,metrain)
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

        use para

	integer id
	integer nvar

	integer il
	character*80 string,string1,string2

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
	    write(6,*) trim(string1)
	    write(6,*) trim(string2)
	    stop 'error stop meteo_set_wind_data: wind description'
	  end if

	end if

        wfact = 1.d0 / rowass
        if( iwtype /= 2 ) wfact = roluft / rowass
        sfact = 1.
        if( iwtype == 4 ) sfact = nmile / 3600.

!	---------------------------------------------------------
!	handle pressure
!	---------------------------------------------------------

	if( nvar == 3 ) then
	  has_pressure = .true.
	  call iff_get_var_description(id,3,string)
	  if( string == papa ) then
	    pfact = 1.
	  else if( string == pamb ) then
	    pfact = 100.
	  else if( string == ' ' ) then		!must determine later
	    pfact = 1.
	  else
	    write(6,*) 'description string for pressure not recognized: '
	    write(6,*) trim(string)
	    stop 'error stop meteo_set_wind_data: pressure description'
	  end if
	  if( bdebug ) then
	    write(6,*) 'pressure initialized: ',pfact,trim(string)
	    write(6,*) 'wind type: ',iwtype
	  end if
	end if

!	---------------------------------------------------------
!	remember values and write to monitor
!	---------------------------------------------------------

        call putpar('iwtype',dble(iwtype))

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

        subroutine meteo_convert_wind_data(id,n,wx,wy,cdv,tx,ty,ws,pp,cice)

        include 'femtime.h'

        integer id
        integer n
        double precision wx(n),wy(n)
        double precision cdv(n)
        double precision tx(n),ty(n)
        double precision ws(n)
        double precision pp(n)
        double precision cice(n)

        logical bnowind,bstress,bspeed
        integer k
        integer itact
        double precision cd,wxymax,txy,wspeed,wdir,fact,fice,aice
        double precision pmin,pmax
        character*80 string

        bnowind = iwtype == 0
        bstress = iwtype == 2
        bspeed = iwtype > 2
        cd = dragco
        wxymax = 0.d0
        aice = 1.       !ice cover for momentum: 1: use  0: do not use
        
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
          else                          !data is wind velocity [m/s]
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

	!call get_act_time(itact)
	!write(112,*) itact,wxymax

        if( wslim > 0 .and. wxymax > wslim ) then !artificially limit wind speed
          !call get_act_time(itact)
          !itact = it
	  !write(111,*) 'limiting wind speed: ',itact,wxymax
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

	if( string == ' ' ) then		!only if not yet determined
	  pmin = minval(pp)
	  pmax = maxval(pp)
	  if( pmin /= 0 .and. pmax /= 0. ) then
	    if( pmin > 85000 .or. pmax < 110000 ) then
	      pfact = 1.
	      string = papa
	    else if( pmin > 850 .and. pmax < 1100 ) then
	      pfact = 100.
	      string = pamb
	    else
	      pfact = 1.
	      string = 'unknown'
	    end if
	    call iff_set_var_description(id,3,string)
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

	subroutine meteo_convert_rain_data(id,n,r)

! convert rain from mm/day to m/s

	integer id
	integer n
	double precision r(n)

	integer i
        double precision zconv
        parameter( zconv = 1.e-3 / 86400. )

	do i=1,n
	  r(i) = r(i) * zconv
	end do

	end subroutine meteo_convert_rain_data

!*********************************************************************

	subroutine meteo_set_ice_data(id,nvar)

        use para

	integer id
	integer nvar

	character*60 string

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
	  if( string == ice ) then
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

! convert ice data (nothing to do)

	use evgeom
	use basin

	integer id
	integer n
	double precision r(n)	!ice concentration

	include 'femtime.h'

	integer k,ie,ii,ia
	double precision racu,rorig
	double precision rice,rarea,area

	racu = 0.
	do k=1,n
	  racu = racu + r(k)
	end do
	rorig = racu / n

	rice = 0.
	rarea = 0.

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
	      rice = rice + area*r(k)
	      rarea = rarea + area
	    end do
	  end if
	end do

	racu = 0.
	do k=1,n
	  racu = racu + r(k)
	end do
	racu = racu / n

	!write(6,*) rorig,racu
	!write(166,*) it,rice/rarea

	end subroutine meteo_convert_ice_data

!*********************************************************************
!*********************************************************************

	subroutine meteo_set_heat_data(id,nvar)

        use para

	integer id
	integer nvar

	character*60 string,strings(4)
	integer i
	character*80 vapor

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

        ihtype = nint(getpar('ihtype'))  
	if( ihtype == 1 ) then
	  vapor = rhum
	else if( ihtype == 2 ) then
	  vapor = wbtm
	else if( ihtype == 3 ) then
	  vapor = dewp
	else
	  write(6,*) 'ihtype = ',ihtype
	  stop 'error stop meteo_set_heat_data: erroneous ivapor'
	end if

	if( strings(1) == ' ' ) then	!TS file or constant
	  if( .not. iff_has_file(id) ) ihtype = 0	!no heat

          if( ihtype > 0 ) then
            call iff_set_var_description(id,1,srad)
            call iff_set_var_description(id,2,tair)
            call iff_set_var_description(id,3,vapor)
            call iff_set_var_description(id,4,ccov)
          end if
        else
          if( strings(1) /= srad ) then
            ihtype = -1
          else if( strings(2) /= tair ) then
            ihtype = -2
          else if( strings(4) /= ccov ) then
            ihtype = -4
	  else
            if( strings(3) == rhum ) then
              ihtype = 1
            else if( strings(3) == wbtm ) then  
              ihtype = 2
            else if( strings(3) == dewp ) then  
              ihtype = 3
            else
              ihtype = -3
            end if
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

        subroutine meteo_convert_heat_data(id,n,metaux,mettair,methum,metwbt,metdew)

	integer id
	integer n
	double precision metaux(n)
	double precision mettair(n)
	double precision methum(n)
	double precision metwbt(n)
	double precision metdew(n)

	logical bnoheat

	bnoheat = ihtype == 0
	
        if( bnoheat ) then              !no heat
	  !nothing to be done
        else
	  call meteo_convert_vapor(ihtype,n,metaux,mettair,methum,metwbt,metdew)
	end if

	end subroutine meteo_convert_heat_data

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine meteo_set_array(n,val0,array)

! initializes 2D array

	implicit none

	integer n
	double precision val0
	double precision array(*)

	integer i

	do i=1,n
	  array(i) = val0
	end do

	end subroutine meteo_set_array

!*********************************************************************

	subroutine meteo_convert_vapor(mode,n,aux,tav,rhv,wbv,dpv)

! computes wet bulb temperature

	use meteo
        use heat_util

	implicit none

	integer mode
	integer n
	double precision aux(n)	!value read
	double precision tav(n)	!dry air temperature
	double precision rhv(n)	!relative humidity
	double precision wbv(n)	!wet bulb temperature
	double precision dpv(n)	!dew point temperature

	integer i
	double precision db,rh,wb,dp

	do i=1,n
	  db = tav(i)
	  rh = rhv(i)
	  wb = wbv(i)
	  dp = dpv(i)
	  if( mode .eq. 1 ) then		!val is humidity
	      rh = aux(i)
	  else if( mode .eq. 2 ) then		!val is wet bulb
	      wb = aux(i)
	  else if( mode .eq. 3 ) then		!val is dew point
	      dp = aux(i)
	  else
	      write(6,*) 'mode = ',mode
	      stop 'error stop meteo_convert_hum: mode'
	  end if
	  call convert_vapor_content(mode,db,rh,wb,dp)
	  rhv(i) = rh
	  wbv(i) = wb
	  dpv(i) = dp
	end do

	end subroutine meteo_convert_vapor

!*********************************************************************
!*********************************************************************
!*********************************************************************

        subroutine get_drag(itdrag,wxy,dragco)

! computes drag coefficient

        implicit none

        integer itdrag          !type of formula
        double precision wxy                !wind speed
        double precision dragco             !computed drag coefficient

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

        subroutine convert_wind(s,d,u,v)

        implicit none

        double precision s,d,u,v

        double precision dir
        double precision pi,rad
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

	subroutine meteo_get_heat_values(k,qs,ta,rh,wb,uw,cc,p)

! returns meteo parameters for one node
!
! pressure is returned in [mb]

	use meteo

	implicit none

        integer k                       !node number
        double precision qs                         !solar radiation [W/m**2]
        double precision ta                         !air temperature [Celsius]
        double precision rh                         !relative humidity [%, 0-100]
        double precision wb                         !wet bulb temperature [Celsius]
        double precision uw                         !wind speed [m/s]
        double precision cc                         !cloud cover [0-1]
        double precision p                          !atmospheric pressure [mbar, hPa]

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

        subroutine meteo_get_heat_extra(k,dp,uuw,vvw)
 
! returns iextra meteo parameters for one node (iheat == 8)
!
! pressure is returned in [mb]

        use meteo

        implicit none

        integer k                       !node number
        double precision dp                         !dew point temperature [Celsius]  
        double precision uuw                        !u-component wind speed [m/s]
        double precision vvw                        !v-component wind speed [m/s]

        dp  = metdew(k)  
        uuw = wxv(k)
        vvw = wyv(k)

        end subroutine meteo_get_heat_extra

!*********************************************************************


	subroutine get_pe_values(k,r,e,eeff)

! returns precipitation and evaporation values
!
! eeff is evaporation used in model, if ievap==0 => eeff==0.

	use meteo
        use para

	implicit none

	integer k
	double precision r			!rain [m/s]
	double precision e			!evaporation [m/s]
	double precision eeff		!effective evaporation [m/s] - used in model

	integer ievap
	save ievap
	data ievap /-1/

        if( ievap .eq. -1 ) ievap = nint(getpar('ievap'))

	r = metrain(k)
	e = evapv(k)
	eeff = e*ievap

	end subroutine get_pe_values

!*********************************************************************

	subroutine set_evap(k,e)

! sets evaporation

	use meteo

	implicit none

	integer k
	double precision e			!evaporation [m/s]

	evapv(k) = e

	end subroutine set_evap

!*********************************************************************

	subroutine meteo_get_solar_radiation(k,qs)

	use meteo

	implicit none

        integer k                       !node number
        double precision qs                         !solar radiation [W/m**2]

	qs = metrad(k)

	end subroutine meteo_get_solar_radiation

!*********************************************************************

        subroutine get_wind(k,wx,wy)

! helper function -> return wind for node k

	use meteo

        implicit none

        integer k
        double precision wx,wy

        wx = wxv(k)
        wy = wyv(k)

        end subroutine get_wind

!*********************************************************************

        subroutine meteo_set_matrix(qs,ta,rh,wb,uw,cc)

! interpolates files spatially - to be deleted

	use meteo
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        double precision qs,ta,rh,wb,uw,cc

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

!*********************************************************************

!================================================================
        end module meteo_forcing
!================================================================
