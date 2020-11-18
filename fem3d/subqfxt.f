
!--------------------------------------------------------------------------
!
!    Copyright (C) 2002-2012,2014-2016,2018-2019  Georg Umgiesser
!    Copyright (C) 2014-2016  Christian Ferrarin
!    Copyright (C) 2016  Ivan Federico
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c heat flux module (file administration)
c
c contents :
c
c subroutine qflux_init
c		initializes routines
c subroutine qflux3d(it,dt,nkn,nlvdi,temp,dq)
c		computes new temperature (forced by heat flux) - 3d version
c
c revision log :
c
c 01.02.2002	ggu	new as administrative routines
c 11.04.2003	ggu	albedo introduced in heat2t
c 16.08.2004	ggu	some comments on usage, heat2t copied to subqfxu4.f
c 22.02.2005	ggu	subroutine qflux2d deleted
c 23.03.2006	ggu	changed time step to real
c 20.08.2007	ggu	prepared for evaporation - salinity feedback
c 17.04.2008	ggu	save evaporation in evapv for later use
c 17.11.2008	ggu	new call to qfcheck_file() before opening
c 10.03.2009	ggu	new qflux_read(), new call to meteo_get_values()
c 11.03.2009	ggu	new routine qflux_compute()
c 27.08.2009	ggu	new call to heatgill, new routine heatareg
c 11.11.2009	ggu	handle abosrption of short rad in more than one layer
c 23.03.2010	ggu	changed v6.1.1
c 04.03.2011	ggu	new routine heatgotm
c 23.03.2011	ggu	new routine check_heat() to check for Nan, new iheat
c 25.03.2011	ggu	new parameters iheat,hdecay,botabs
c 14.04.2011	ggu	changed VERS_6_1_22
c 14.02.2012	ggu	changed VERS_6_1_44
c 30.03.2012	ggu	changed VERS_6_1_51
c 29.04.2014	ccf	read qsens, qlat, long from file
c 16.06.2014	ccf	new routine heatcoare, which also update wind stress
c 20.06.2014	ccf	new routine for computing sea surface skin temperature
c 18.07.2014	ggu	changed VERS_7_0_1
c 05.11.2014	ggu	changed VERS_7_0_5
c 26.11.2014	ggu	changed VERS_7_0_7
c 05.12.2014	ggu	changed VERS_7_0_8
c 19.01.2015	ggu	changed VERS_7_1_3
c 26.02.2015	ggu	changed VERS_7_1_5
c 30.04.2015	ggu	changed VERS_7_1_9
c 05.06.2015	ggu	changed VERS_7_1_12
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 18.09.2015	ccf	do not compute heat fluxes in dry nodes
c 18.09.2015	ccf	checks only for levdbg > 2
c 26.10.2015	ggu	critical omp sections introduced (eliminated data race)
c 18.12.2015	ggu	changed VERS_7_3_17
c 07.04.2016	ggu	compute total evaporation
c 15.04.2016	ggu	changed VERS_7_5_8
c 04.05.2016	ccf	do not pass albedo into heat2t
c 04.05.2016	ggu	include effect of ice cover
c 25.05.2016	ggu	changed VERS_7_5_10
c 21.07.2016	ivn	isolp = 1, 2 length scale solar penetration (Jerlov) 
c 09.09.2016	ggu	changed VERS_7_5_17
c 29.09.2016	ivn	bug fix for isolp = 1
c 24.01.2018	ggu	changed VERS_7_5_41
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.10.2018	ggu	changed VERS_7_5_50
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 10.12.2019	ggu	ice cover introduced, handle ice-water sensible heat
c 11.11.2020	ggu	new ice model integrated, cleaned, documented
c 13.11.2020	ggu	bug fix: correct rain from [m/s] to [mm/d]
c 14.11.2020	ggu	allow for ice and other heat fluxes to coexist
c
c notes :
c
c qflux_init is called in subn11.f
c qflux_read is called in subn11.f
c qflux3d is called in newbcl.f
c qfget is called in qflux3d (here)
c
c to use heat flux module insert in section "name" of STR file:
c
c $name
c     qflux = 'qflux.dat'
c $end
c
c if such a file is given the heat flux module is used and works
c   without any other intervention (ibarcl>0)
c
c other notes in subqfxf.f
c
c*****************************************************************************

	subroutine qflux_compute(byes)

c returnes flag -> compute heat flux or not

	implicit none

	logical byes

	character*80 file

        call getfnm('qflux',file)

	byes = ( file /= ' ' )

	end

c*****************************************************************************

	subroutine qflux_init

c initializes routines

	implicit none

	character*80 file

        call getfnm('qflux',file)

	if( file .ne. ' ' ) then
	  call qfcheck_file(file)
	  call qfinit(file)
	end if

	end

c*****************************************************************************

	subroutine qflux_read(it)

c reads new meteo data

	implicit none

	integer it

	call qfmake(it)

	end

c*****************************************************************************

	subroutine qflux3d(dtime,dt,nkn,nlvddi,temp,dq)

c computes new temperature (forced by heat flux) - 3d version

	use mod_meteo
	use mod_ts
	use levels
        use basin, only : xgv,ygv  
        use mod_hydro_print  

	implicit none

        include 'subqfxm.h'

	double precision dtime
	real dt
	integer nkn
	integer nlvddi
        real temp(nlvddi,nkn) 
	double precision dq	!total energy introduced [(W/m**2)*dt*area = J]

        logical byes
        logical buseice,bicecover
        integer levdbg
	integer k
	integer l,lmax,kspec
	integer mode
	integer days,im,ih
	integer ys(8)
	real hm,sm,tm,tnew
	real salt,tfreeze
	real albedo
	real adecay,qsbottom
	real qtot
        real qs,ta,tb,uw,cc,ur,p,e,r,q
        real ddlon,ddlat  
        real dp,uuw,vvw  
	real qsens,qlat,qlong,evap,qrad,qsdown
        real qswa,qice,qsurface
	real cice,fice_free,fice_cover
	real ev,eeff
	real area
	real evaver
        real uub,vvb        
	real ik1,ik2
	real tice,u,v,uv

	double precision ddq
	double precision atime
	character*20 aline

	real, save, allocatable :: dtw(:)	!Warm layer temp. diff
	real, save, allocatable :: tws(:)	!Skin temperature (deg C)

	real hb			!depth of modelled T in first layer [m]
	real usw		!surface friction velocity [m/s]
	real qss		!Net shortwave flux 
	real cd			!wind drag coefficient

c     coefficient after Paul e Simon (1977) up to type IV
c     for coastal water (1-3-5-7-9) coefficient fitting data with Jerlov(1968)
      real,parameter, dimension(10) :: Rt =
     +          (/ 0.58, 0.62, 0.67, 0.77, 0.78,
     +         0.605, 0.61, 0.37, 0.338, 0.277 /)
      real,parameter, dimension(10) :: k1 =
     +         (/ 0.35, 0.6,  1.0,  1.5,  1.4, 
     +         0.387, 0.353, 2.395, 1.8, 1.505 /)
      real,parameter, dimension(10) :: k2 =
     +       (/ 23.0, 20.0, 17.0, 14.0, 7.9,
     +          5.05, 3.55, 0.335, 0.329, 0.325 /)

	integer itdrag

c functions
	real depnode,areanode,getpar
	integer ifemopa
	logical is_dry_node
c save
	integer, save :: icall = 0
	integer, save :: n93 = 0

	integer, save :: iheat
	integer, save :: isolp  
	integer, save :: iwtyp
	real, save :: hdecay
	real, save :: botabs

	real, save :: aice = 1.
	real, save :: fice_pen = 0.3		!penetration through ice

	integer, save :: kdebug = 0		!debug for this node
	logical, save :: bdebug = .false.
	logical, save :: baverevap = .false.
	logical, save :: bwind = .false.
	logical, save :: bice = .false.
	logical, save :: bheat = .false.
	logical, save :: bqflux = .false.

c---------------------------------------------------------
c start of routine
c---------------------------------------------------------

	dq = 0.
	if( icall < 0 ) return

c---------------------------------------------------------
c iheat		1=areg  2=pom  3=gill  4=dejak  5=gotm  
c		6=COARE 3.0
c               7=read flux from file
c		8=MFS routines (Pettenuzzo et al., 2010)
c hdecay	depth of e-folding decay of radiation
c		0. ->	everything is absorbed in first layer
c botabs	1. ->	bottom absorbs remaining radiation
c		0. ->	everything is absorbed in last layer
c---------------------------------------------------------
c format of heat file containing time series (4 data columns):
c    time srad airt rhum cc
c in case of iheat==7 the columns are:
c    time srad qsens qlat qlong
c---------------------------------------------------------

c---------------------------------------------------------
c initialize on first call
c---------------------------------------------------------

	if( icall .eq. 0 ) then
	  iheat = nint(getpar('iheat'))
          isolp = nint(getpar('isolp'))   
	  iwtyp  = nint(getpar('iwtyp'))+1
	  hdecay = getpar('hdecay')
	  botabs = getpar('botabs')

	  call shyice_init
	  call shyice_is_active(bice)		!ice model is active
	  bheat = iheat > 0			!we have to compute heat fluxes
	  call qflux_compute(bqflux)		!have qflux file

	  if( .not. bqflux ) icall = -1
	  if( .not. bheat .and. .not. bice ) icall = -1
	  if( icall < 0 ) return

!$OMP CRITICAL
	  write(6,*) 'qflux3d routines are active'
	  write(6,*) 'qflux3d: bqflux,bheat,bice: '
     +				,bqflux,bheat,bice
	  write(6,*) 'qflux3d: iheat,hdecay,botabs: '
     +				,iheat,hdecay,botabs  
!$OMP END CRITICAL

	  allocate(dtw(nkn))
	  allocate(tws(nkn))
	  dtw = 0.
	  tws(:) = temp(1,:)

          itdrag = nint(getpar('itdrag'))
	  bwind = itdrag .eq. 4
          if ( bwind ) then   
            if ( iheat .eq. 6 .or. iheat .eq. 8 ) then   
              write(6,*) 'itdrag = ',itdrag
              write(6,*) 'iheat = ',iheat
            else
              write(6,*) 'Erroneous value for itdrag = ',itdrag
              write(6,*) 'Use itdrag = 4 only with iheat = 6 or 8'
              stop 'error stop qflux3d: itdrag'
            end if
          end if

	  call shyice_init_output
	  call meteo_get_ice_usage(aice)	!use ice (1) or not (0)
	  if( bice ) fice_pen = 0.

          levdbg = nint(getpar('levdbg'))
	  bdebug = levdbg .ge. 2 

	end if

	icall = icall + 1

c---------------------------------------------------------
c set other parameters
c---------------------------------------------------------

	mode = +1	!use new time step for depth

	adecay = 0.
	if( hdecay .gt. 0. ) adecay = 1. / hdecay

	ik1 = 0.
	ik2 = 0.
	if ( isolp .eq. 1 ) then
          if (iwtyp > 10 ) then
            write(6,*) 'not valid value for iwtyp = ',iwtyp-1
            write(6,*) 'use values from 0 (Jrv t-I) to 10 (Jrv t-9)'
            stop 'error stop qflux3d: iwtyp'
          end if
          ik1 = 1. / k1(iwtyp)
          ik2 = 1. / k2(iwtyp)
	end if

c---------------------------------------------------------
c set date parameters
c---------------------------------------------------------

	call get_act_timeline(aline)
	call get_absolute_act_time(atime)
	call dts_from_abs_time_to_ys(atime,ys)
	call dts_from_abs_time_to_days_in_year(atime,days)
	days = days - 1
	im = ys(2)
	ih = ys(4)

c---------------------------------------------------------
c loop over nodes
c---------------------------------------------------------

        ddq = 0.

	do k=1,nkn

	  if (is_dry_node(k)) then	!do not compute if node is dry
	    dtw(k)   = 0.
	    tws(k)   = temp(1,k)
	    evapv(k) = 0.
	    cycle
	  end if

	  tm = temp(1,k)
	  salt = saltv(1,k)
	  call get_ice_cover(k,cice)
	  fice_cover = aice*cice
	  fice_free  = 1. - fice_cover
	  bicecover = fice_cover > 0.
	  buseice = bice .and. bicecover	!we use ice model
	  tfreeze = -0.0575*salt
	  qsens = 0.
	  qlat = 0.
	  qlong = 0.
	  evap = 0.
	  area = areanode(1,k)
	  lmax = ilhkv(k)
          if( isolp .eq. 0 .and. hdecay .le. 0. ) lmax = 1   

	  !------------------------------------------------
	  ! solar radiation - positive from air to sea
	  !------------------------------------------------

          call meteo_get_heat_values(k,qs,ta,ur,tb,uw,cc,p)
	  call make_albedo(tm,albedo)
	  qsdown = qs * (1. - albedo)
	  qss = fice_free*qsdown + fice_cover*fice_pen*qsdown
	  if( .not. bheat ) qss = 0.

	  !------------------------------------------------
	  ! heat exchange - qlong, qlat, qsens positive from sea to air
	  !------------------------------------------------

	  if( iheat .eq. 0 ) then
	    buseice = bice			!we use ice model heat fluxes
	  else if( iheat .eq. 1 ) then
	    call heatareg (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 2 ) then
	    call heatpom  (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 3 ) then
	    call heatgill (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 4 ) then
	    !call rh2wb(ta,p,ur,tb)	!FIXME
	    call heatlucia(ta,p,uw,tb,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 5 ) then
	    call heatgotm (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 6 ) then
	    call get_pe_values(k,r,ev,eeff)
	    call heatcoare(ta,p,uw,ur,cc,tws(k),r,qss,qsens,qlat,
     +                     qlong,evap,cd)
	    if ( bwind ) windcd(k) = cd
	  else if( iheat .eq. 7 ) then
	    qsens = ta
	    qlat  = ur
	    qlong = -cc	  !change sign of long wave radiation given by ISAC
	    evap  = qlat / (2.5008e6 - 2.3e3 * tm)	!pom, gill, gotm
          else if( iheat .eq. 8 ) then
            ddlon = xgv(k)    
            ddlat = ygv(k)   
            uub = uprv(1,k)  
            vvb = vprv(1,k)  
	    call meteo_get_heat_extra(k,dp,uuw,vvw)
            call heatmfsbulk(days,im,ih,ddlon,ddlat,ta,p,uuw,vvw,dp,
     +                   cc,tm,uub,vvb,qsens,qlat,qlong,evap,qswa,cd)   
            qss= fice_free*qswa  !albedo (monthly) already in qshort1 -> qswa  
            if ( bwind ) windcd(k) = cd   
          else
            write(6,*) 'iheat = ',iheat
            stop 'error stop qflux3d: value for iheat not allowed'
          end if

	  if (bdebug) call check_heat(k,tm,qsens,qlat,qlong,evap)

	  !------------------------------------------------
	  ! heat transfer from water to ice - qice positive from water to ice
	  !------------------------------------------------

	  qice = 0.
	  if( bicecover .and. .not. bice ) then	!ice on node but no ice model
	    tice = 0.
	    call getuv(1,k,u,v)
	    uv = sqrt(u*u+v*v)		!current speed
	    call ice_water_exchange(tm,tice,uv,qice)
	    if( .not. bheat ) qice = 0.
	  end if

	  !------------------------------------------------
	  ! final heat exchange - qrad, qss, qtot positive from air to sea
	  !------------------------------------------------

	  qlong = fice_free * qlong
	  qlat = fice_free * qlat
	  qsens = fice_free * qsens + fice_cover * qice
          qrad =  - ( qlong + qlat + qsens )
	  qtot = qss + qrad

	  !------------------------------------------------
	  ! distribute short wave radiation to layers and compute heat budget
	  !------------------------------------------------

	  qsurface = qrad

          do l=1,lmax
	    if( .not. bheat ) cycle	!no heat flux computed
	    if( buseice ) cycle		!we use heat flux from ice model
            hm = depnode(l,k,mode)
            tm = temp(l,k)
            if (isolp == 0) then  
              if( hdecay .le. 0. ) then
                qsbottom = 0.
              else
                qsbottom = qss * exp( -adecay*hm )
                if( l .eq. lmax ) qsbottom = botabs * qsbottom
              end if
            else if (isolp == 1) then
              !Jerlov classification (1968) from Type I to 9
               qsbottom = qss * (Rt(iwtyp) * exp(- ik1 * hm )
     +                            +  (1-Rt(iwtyp)) * exp(- ik2 * hm))
              if( l .eq. lmax ) qsbottom = botabs * qsbottom
            else
              write(6,*) 'Erroneous value for isolp = ',isolp
              write(6,*) 'Use isolp = 0 (hdecay) or 1 (Jerlov t-I)'
              stop 'error stop qflux3d: isolp'
            end if
            call heat2t(dt,hm,qss-qsbottom,qsurface,tm,tnew)
            if (bdebug) call check_heat2(k,l,qss,qsbottom,qsurface,
     +                                   albedo,tm,tnew)
            tnew = max(tnew,tfreeze)
            temp(l,k) = tnew
            albedo = 0.
            qsurface = 0.	!different from 0 only in first layer
            qss = qsbottom
          end do

c         ---------------------------------------------------------
c         compute sea surface skin temperature
c         ---------------------------------------------------------

	  tm   = temp(1,k)
	  hb   = depnode(1,k,mode) * 0.5
          usw  = max(1.e-5, sqrt(sqrt(tauxnv(k)**2 + tauynv(k)**2)))
	  call tw_skin(qss,qrad,tm,hb,usw,dt,dtw(k),tws(k))

c         ---------------------------------------------------------
c         handle ice model
c         ---------------------------------------------------------

	  if( k == kdebug ) write(444,*) icall,bice,bicecover,buseice
	  if( bice ) then
            hm = depnode(1,k,mode)
            tm = temp(1,k)
            sm = saltv(1,k)
	    call get_pe_values(k,r,e,eeff)
	    r = r * 1000 * 86400	!r is in [m/s] -> convert to [mm/d]
	    call shyice_run(k,qs,ta,ur,uw,cc,p,r,hm,tm,sm,dt)
	    if( buseice ) then	!ice on water or bheat == .false.
              temp(1,k) = tm
              saltv(1,k) = sm
	      evap = 0.
	      dtw(k) = 0.
	      call shyice_get_tsurf(k,tws(k))
	    end if
	  end if

c	  ---------------------------------------------------------
c	  evap is in [kg/(m**2 s)] -> convert it to [m/s]
c	  evap is normally positive -> we are loosing mass
c	  ---------------------------------------------------------

	  evap = evap / rhow			!evaporation in m/s
	  evapv(k) = evap			!positive if loosing mass

          ddq = ddq + qtot * dt * area

	end do

	dq = ddq

c---------------------------------------------------------
c compute total evaporation
c---------------------------------------------------------

	if( baverevap ) then
	  call aver_nodal(evapv,evaver)	!in evaver is average of evaporation m/s
	  write(678,*) dtime,evaver
	end if

c---------------------------------------------------------
c special output
c---------------------------------------------------------

	call shyice_write_output

	k = min(nkn,1000)
	k = -1
	if( k .gt. 0 ) then
	  if( n93 .eq. 0 ) then
	    n93 = ifemopa('opening file 93','.93','form','new')
	  end if
	  write(n93,*) 'qflux3d: ',dtime,temp(1,k)
	end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c*****************************************************************************

	subroutine ice_water_exchange(tw,ti,uv,qice)

! Li at al.,
! Heat transfer at ice-water interface under conditions of low flow velocities
! Journal of Hydrodynamics,2016,28(4):603-609
! DOI: 10.1016/S1001-6058(16)60664-9

! qice is heat flux [W/m**2] from water to ice (positive)

	implicit none

	real tw,ti
	real uv
	real qice

	real, parameter :: rhow = 1000.	!density of water
	real, parameter :: cw = 4179.6	!specific heat of water
	real, parameter :: ch = 1.1E-3	!heat transfer coefficient

	qice = rhow*cw*ch*uv*(tw-ti)

	end

c*****************************************************************************

	subroutine check_heat(k,tm,qsens,qlat,qlong,evap)

	use mod_debug

	implicit none

	integer k
	real tm,qsens,qlat,qlong,evap

	if( is_nan(tm) ) goto 98
	if( is_nan(qsens) ) goto 98
	if( is_nan(qlat) ) goto 98
	if( is_nan(qlong) ) goto 98
	if( is_nan(evap) ) goto 98

	return
   98	continue
	write(6,*) k,tm,qsens,qlat,qlong,evap
	stop 'error stop check_heat: Nan found'
	end

c*****************************************************************************

	subroutine check_heat2(k,l,qs,qsbottom,qrad,albedo,tm,tnew)

	use mod_debug

	implicit none

	integer k,l
	real qs,qsbottom,qrad,albedo,tm,tnew

	if( is_nan(qs) ) goto 98
	if( is_nan(qsbottom) ) goto 98
	if( is_nan(qrad) ) goto 98
	if( is_nan(albedo) ) goto 98
	if( is_nan(tm) ) goto 98
	if( is_nan(tnew) ) goto 98

	if( abs(tm) .gt. 100. ) goto 97
	if( abs(tnew) .gt. 100. ) goto 97

	return
   97	continue
	write(6,*) k,l,tm,tnew
	stop 'error stop check_heat2: t out of range'
   98	continue
	write(6,*) k,l,qs,qsbottom,qrad,albedo,tm,tnew
	stop 'error stop check_heat2: Nan found'
	end

c*****************************************************************************

