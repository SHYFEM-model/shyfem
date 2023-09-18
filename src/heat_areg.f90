!
! $Id: subqfxm3.f,v 1.1 2009-09-14 08:20:58 georg Exp $
!
! heat flux module (AREG)
!
! contents :
!
! subroutine heatareg(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)
!	computes heat fluxes from formulas in AREG
! subroutine heatpom(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)
!	computes heat fluxes from formulas in AREG (different implementation)
! subroutine vapor1(tt,p,u,e,r,q)
!	computes values for water vapor given relative humidity
! subroutine areg_fluxes(sst,tair,rhnow,speed,cldnow,qsw,qbw,ha,elat,qsurf,evap)
!	this is an adaptation of AREG heat flux routines
!
! revision log :
!
! 27.08.2009    dbf&ggu	integrated into model
! 16.02.2011    ggu	pstd introduced
!
!*******************************************************************
!-------------------------------------------------------------------
        module heat_areg
!-------------------------------------------------------------------
        contains
!-------------------------------------------------------------------

	subroutine heatareg(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)

! computes heat fluxes from formulas in AREG
!
! heat fluxes are positive upward (from sea to atmosphere)

	implicit none

	double precision t		!air temperature [C]			- in
	double precision p		!pressure [mb]				- in
	double precision w		!wind speed [m/s]			- in
	double precision ur		!relative humidity [%] ([0-100])	- in
	double precision cc		!cloud cover [0-1]			- in
	double precision ts		!sea temperature [C]			- in
	double precision qsens	!sensible heat flux [W/m**2]		- out
	double precision qlat	!latent heat flux [W/m**2]		- out
	double precision qlong	!long wave radiation [W/m**2]		- out
	double precision evap	!evaporation [kg/m**2/s]		- out

	double precision sst,tair,rhnow,speed,cldnow
        double precision qsw,qbw,ha,elat,qsurf,evap1

! qsw,qsurf,evap1 are dummy variables

	sst = ts
	tair = t
	rhnow = ur
	speed = w
	cldnow = cc

	call areg_fluxes(sst,tair,rhnow,speed,cldnow,qsw,qbw,ha,elat,qsurf,evap1)

	qlong = qbw
	qsens = ha
	qlat = elat
	evap  = qlat / 2.5e+6

	end

!*******************************************************************

	subroutine heatpom(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)

! computes heat fluxes from formulas in AREG (different implementation)
!
! heat fluxes are positive upward (from sea to atmosphere)

        use heat_util2

	implicit none

	include 'subqfxm.h'

	double precision t		!air temperature [C]			- in
	double precision p		!pressure [mb]				- in
	double precision w		!wind speed [m/s]			- in
	double precision ur		!relative humidity [%] ([0-100])	- in
	double precision cc		!cloud cover [0-1]			- in
	double precision ts		!sea temperature [C]			- in
	double precision qsens	!sensible heat flux [W/m**2]		- out
	double precision qlat	!latent heat flux [W/m**2]		- out
	double precision qlong	!long wave radiation [W/m**2]		- out
	double precision evap	!evaporation [kg/m**2/s]		- out

! other variables:
!
!       double precision q          !specific humidity [0-1]
!       double precision e          !vapor pressure [mb]
!       double precision r          !mixing ratio [0-1]

	double precision rdry,lv0,pascal
	parameter(rdry=287.04,lv0=2.5008e6,pascal=100.)

	double precision rho,cp,lv,tv
	double precision rhod	!density of moist air DEB
	double precision pw		!partial pres of water vap in moist air
	double precision pa		!atm pres of moist air
	double precision s0,s
	double precision e,r,q
	double precision es,rs,qs
	double precision dt,ch,ce
	double precision theta,theta2,theta4,thetas
	double precision cloud
	logical bignami,may

	bignami=.false.
	may=.true.

!	------------------------------------------------
!	initialization
!	------------------------------------------------

	call vapor1(t,p,ur,e,r,q)        !compute e,r,q

	tv = ( t + kelv ) * ( 1. + 0.6078 * q )
	lv = lv0 - 2.3e3 * ts
	rho = pascal * p / ( rdry * tv )
	rhod =rho * ( 1 + r )/( 1 + 1.609 * r )
	cp = 1004.6 * ( 1. + 0.8375 * q )

	call satur(ts,p,es,rs,qs)	!compute saturation values for sea

!	------------------------------------------------
!	compute Kondo related parameters
!	------------------------------------------------

	dt = ts - t
	
!	here kondo formulation for ch has to be used

	if(w.ne.0)then	
		s0= dt/w*w !this formula is incomplete (wind at 10 m)
	else
		s0= dt
	endif
	s= s0 * abs(s0)/(abs(s0)+0.01)

	ce = 0.
	ch = 0.
	if( dt .lt. 0 ) then				!stable
	   if(s .lt. 0 .and. s .gt. -3.3) then
		ch= 1.24 * (0.1 + 0.03 * s + 0.9 * exp(4.8 * s))/1000.
		ce= 1.28 * (0.1 + 0.03 * s + 0.9 *  exp(4.8 * s))/1000.
	   else if( s .le. -3.3) then
		ch = 0.
		ce = 0.
	   endif
	else						!unstable
		ch = 1.24 * (1+0.63 * s**(0.5))/1000.
		ce = 1.28 * (1+0.63 * s**(0.5))/1000.
	endif

!	------------------------------------------------
!	sensible and latent heat flux
!	------------------------------------------------

	qsens = rhod * cp *ch * w * dt

	evap = rhod * ce * w * (qs - q)
	qlat = lv * evap

!	------------------------------------------------
!	long wave radiation
!	------------------------------------------------

	theta = t + kelv
	thetas =  ts + kelv
	theta2 = theta * theta
	theta4 = theta2 * theta2
	cloud = 1. - 0.75 * cc**(3.4)

	if(may) then
	  qlong = (bolz * theta4 * (0.4 - 0.05 * sqrt(e) ) + 4 * bolz * &
     &	  theta2 * theta *( thetas - theta)) * cloud
	else if(bignami) then
	  qlong=0.98*bolz*thetas**4-bolz*theta4*(0.653+0.00535*e)*      &
     &	  (1+0.1762*cc*cc)
	else
	  stop 'error stop heatpom: must choose either may or bignami'
     	endif

!	------------------------------------------------
!	end of routine
!	------------------------------------------------

        if( qlong .le. 0. ) then
	  write(6,*) 'qlong <= 0'
          write(6,*) qsens,qlat,qlong
          write(6,*) ts,theta,theta4,thetas
          write(6,*) e,sqrt(e),0.05*sqrt(e),cc,cloud
          write(6,*) t,w,q,qs
          write(6,*) tv,lv,cp,rho
          !stop 'error stop heatpom: qlong'
        end if

	end

!*******************************************************************

	subroutine vapor1(tt,p,u,e,r,q)

! computes values for water vapor given relative humidity

	implicit none

	double precision tt         !air temperature [C]			- in
	double precision p          !pressure [mb]				- in
	double precision u          !relative humidity [%] ([0-100])	- in
	double precision e          !vapor pressure [mb]			- out
	double precision r		!mixing ratio [0-1]			- out
	double precision q		!specific humidity [0-1]		- out

	double precision ew,rw,qw			!saturation vapor
	double precision eps,fw			!molecular weight ratio
	parameter(eps=0.62197)
	double precision fice
	parameter (fice=0.00422)		!over
	double precision fsalt
	parameter (fsalt=0.98)
	double precision aux
	
	double precision a(7)
!	data a / 6984.505294,-188.9039310,2.133357675,-1.288580973e-2,
!     +		4.393587233e-5,-8.023923082e-8,6.136820929e-11  /
     	data a / 6.107799961, 4.436518521e-1, 1.428945805e-2,   &
     &      2.650648471e-4, 3.031240396e-6, 2.034080948e-8,6.136820929e-11 /

     	ew=a(1)+tt*(a(2)+tt*(a(3)+tt*(a(4)+tt*(a(5)+tt*(a(6)+tt*a(7))))))
	aux = (0.7859+0.03477*tt)/(1.+0.00412*tt)
	aux = aux + fice * tt
	fw = 1. + 1.e-6 * p * (4.5+0.0006*tt*tt)
	rw = eps / ( p/ew - 1. )
	qw = rw / ( rw +1. )
	r = 0.01 * u * rw
	q = r / ( r + 1. )
	e=0.01*u*ew

	end
	
!************************************************************************

!************************************************************************
! ADAPTATION OF AREG2 SUBROUTINE
!       subroutine fluxes(alat,alon,imonth,iday,hour,
!     *  sst,tnow,rhnow,speed,cldnow,qsw,qbw,ha,elat,
!     *  qsurf,taux,tauy,rmod,evap)
!************************************************************************

       subroutine areg_fluxes(sst,tair,rhnow,speed_in,cldnow,qsw,qbw,ha,elat,qsurf,evap)

	double precision sst,tair,rhnow,speed_in,cldnow
        double precision qsw,qbw,ha,elat,qsurf,evap

! this is an adaptation of AREG heat flux routines

!	sst		sea surface temperature [C]		- in
!	tair		air temperature [C]			- in
!	rhnow		relative humidity [%] ([0-100])		- in
!	speed_in	wind speed [m/s]			- in
!	cldnow		cloud cover [0-1]			- in
!	qsw		short wave radiation [W/m**2]		- out
!	qwb		long wave radiation [W/m**2]		- out
!	ha		sensible heat flux [W/m**2]		- out
!	elat		latent heat flux [W/m**2]		- out
!	qsurf		net surface heat flux [W/m**2]		- out
!	evap		evaporation rate [m/s]			- out
!
!	qwb,ha,elat are positive from sea to atmosphere
!	qsw is positive from atmosphere to sea
!	qsurf is positive from atmosphere to sea (qsurf=qsw-qwb-ha-elat)
!
!	originally qsw is computed inside this subroutine
!	...however, we use value computed outside
!	...therefore qsw is set to 0 inside this subroutine
!
!  1)   Water flux (WFLUX)                 [ watt/m*m ]
!  2)   Short wave flux (QSW)              [ watt/m*m ]
!  3)   Long wave flux backward (QBW)      [ watt/m*m ]
!  4)   Latent heat of evaporation (ELAT)  [ watt/m*m ]
!  5)   Sensible heat flux   (HA)          [ watt/m*m ]
!  6)   Net Upward flux (QU)               [ watt/m*m ]
!  7)   Net surface heat flux (QSURF)      [ watt/m*m ]
!  8)   Wind stress x-component   (TAUX)   [ newton/m*m ]               
!  9)   Wind stress y-component   (TAUY)   [ newton/m*m ]               
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! ^                                                                  ^
! ^  All the calculations are made in MKS                            ^
! ^  The unit conversion in CGS is made at the end !                 ^
! ^                                                                  ^
! ^  1 Newton      = 10**5 dyn                                       ^
! ^  1 Joule       = 10**7 erg                                       ^
! ^  1 Newton/m*m  = 10 dyn/cm*cm                                    ^
! ^  1 W/m*m       = 10**3 erg/sec*cm*cm                             ^
! ^  1 W/m*m = 1 J/(m*m*s) = 1/(4.19*10**4) cal/(cm*cm*s)            ^
! ^                                                                  ^
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
! bug fixes (by ggu):
!
!	if speed == 0 => speed = 0.01
!	do not change tair -> tnow
!	initialize qsw
!	cavaleri??????
!
!---------------------------------------------------------------------

      include 'subqfxm.h'

      parameter ( pstd = 1013.25 )
      parameter (const = 0.622/pstd)
      parameter (cgs = 1.e-4/4.19)
      parameter (cavaleri=1.00)
! ---- DENSITY OF PURE WATER ---------
      parameter (rho0 = 1000.e0)
!---------------------------------------------------------------------
!     coefficients ( in MKS )  :
!---------------------------------------------------------------------
!
! --- turbulent exchange coefficients
      data  ce,ch  / 1.1e-3, 1.1e-3/    
!  
! --- surface air pressure
      data ps / pstd/
!
! --- air density
      data airden   /1.2/
!     
! --- conversion factor for Langley                                            
!      data  watlan /1.434e-3/       
!
! --- conversion factor for CGS
!      data cgs /1.e+3/      
       data cgstau /10./    
!
       data precip  /0.0/
!	if(speed.gt.0)then !DEB
!
!	print*, 'fluxes : hour =', hour
!
! --- compute speed for bulk coeff.
!                                              
! DEB      speed = sqrt(unow*unow + vnow*vnow)                       
	speed = speed_in
	if( speed .lt. 0.01 ) speed = 0.01
      speedC=speed*cavaleri
! --- precipitation at the moment is supposed to be absent                        
!     precip = 0.0                                                              
! --- SST data converted in Kelvin degrees
!
	ckelv = kelv
      sstk = sst + ckelv 
!
!     TAIR=TNOW is already in Kelvin deg.
!
      tnow = tair + ckelv !DEB
!
! --- calculates the Saturation Vapor Pressure at air temp. and at sea temp.
! ---      esat(Ta) , esat(Ts)
!
      esatair = esk(tnow)
      esatoce = esk(sstk)
!
! --- calculates the saturation mixing ratios at air temp. and sea temp.
! --- wsat(Ta) , wsat(Ts)
!
      wsatair = (const06/ps) * esatair
      wsatoce = (const06/ps) * esatoce
!
! --- calculates the mixing ratio of the air
! --- w(Ta)
!
      wair = 0.01 * rhnow * wsatair
!
! --- calculates the virtual temperature of air
!
      vtnow = (tnow*(const06+wair))/(const06*(1.+wair))
!
! --- calculates the density of the moist air
!
      rhom = 100.*(ps/rgas)/vtnow
!
! ------------------ (i)      short wave
!
	qsw = 0.	!we must give a value to this

!       if (hour .eq. 0.)then
! --- at hour=0 qsw results always 0
!         qsw = 0.
!       else
! DEB           qsw=qshort(imonth,iday,hour,alat,alon,cldnow)
!       endif
!
! ------------------ (ii)      long  wave
!
! --- calculates the term :      Ea**1/2 = [r*Esat(Ta)]**1/2
      ea = sqrt(0.01*rhnow*esatair)
!
! --- calculates the NET LONGWAVE FLUX ( watt/m*m )
! --- WE USE THE qb4 OPTION
!
! --- May formula from May (1986) :
!
       qbw = (1.-0.75*(cldnow**3.4)) * (bolz*(tnow**4)*(0.4 -0.05*ea)   &
     &             + 4.*bolz*(tnow**3)*(sstk-tnow))
!
!
! ------------------ (iii)   sensible heat
!
! --- calculates the term :      ( Ts - Ta )
!
      deltemp = sstk - tnow
!

!
! --- WE USE THE ceh3 OPTION
!
! --- variable turbulent exchange coefficients ( from Kondo 1975 )
!
! --- calculate S :
!
             s=deltemp/(speed**2.)
!
! --- calculate the Stability Parameter :
!
             stp=s*abs(s)/(abs(s)+0.01)
!
! --- for stable condition :
!
	fe = 0.
	fh = 0.

             IF(s.lt.0.) THEN
                if((stp.gt.-3.3).and.(stp.lt.0.)) then
                  fh=0.1+0.03*stp+0.9*exp(4.8*stp)
                  fe=fh
                else
                 if(stp.le.-3.3) then
                  fh=0.
                  fe=fh
!
                endif
                endif
!
! --- for unstable condition :
!
             ELSE
                fh=1.0+0.63*sqrt(stp)
                fe=fh
             ENDIF
!
! --- calculate the coefficient CH,CE,CD
!
             ch=1.3e-03*fh
             ce=1.5e-03*fe
!
! --- calculates the SENSIBLE HEAT FLUX in CGS ( watt/m*m )
!
!      HA = airden*cpw*ch*speed*deltemp
      !ha = rhom*cpw*ch*speed*deltemp
      ha = rhom*cpa*ch*speed*deltemp
!	print*,'ha  ',rhom,'cpw',cpw,'ch ',ch,'sp ',speed,'delt ',deltemp
!
! ------------------ (iv)  latent heat
!
!     const = .622/pstd

!
! --- calculates the term :      esat(Ts)-r*esat(Ta)
!
      evap  = esatoce - rhnow*0.01*esatair
      dsatp = evap
!
! --- calculates the term :      Ce*|V|*[esat(Ts)-r*esat(Ta)]0.622/pstd
!
! --- evap is not the evaporation 
      evap = ce*speed*evap*const

! --- calculates the water flux
!
!      WFLUX =  evap - precip
!
! --- calculates the LATENT HEAT FLUX  LEa ( watt/m*m )
! ---     LEa = L*rho*Ce*|V|*[esat(Ts)-r*esat(Ta)]0.622/pstd
!
!      ELAT = airden*evap*heatlat(sst)
      elat = rhom*evap*heatlat(sst)
!      print*,'deb elat','rhom ',rhom,'ev ',evap,heatlat(sst),'sst ',sst
! --- evap is now (and ONLY NOW) the EVAPORATION RATE (m/s)
      evap = elat/(rho0*heatlat(sst))
!
! ------------------ net upward flux
!
! --- calculates : Qu = Qb + Ha + LEa  ( Rosati,Miyakoda 1988 ; eq. 3.2 )
!
      qu = qbw + ha + elat
!
! ------------------ net surface heat flux (total budget)
!
! --- calculates : Q = Qs - Qu  ( Rosati,Miyakoda 1988 ; eq. 3.2 )
!
      qsurf = qsw - qu
! --- convert in CGS --
!      QSURF = QSURF * cgs
!
! ------------------ wind stresses
!
! --- calculates the Drag Coefficient
!
!       deltem = -deltemp   ! WRONG

        cdx = cd(speedC,deltemp)
!
! --- calculates the wind stresses in MKS ( newton/m*m )

! ---            taux= rho*Cd*|V|u      tauy= rho*Cd*|V|v
!
!       TAUX= airden*cdx*speed*unow
!       TAUY= airden*cdx*speed*vnow
        
!        taux = rhom*cdx*speedC*unow*cavaleri
!        tauy = rhom*cdx*speedC*vnow*cavaleri
!
! --- convert in CGS ---
!        taux = cgstau*taux
!        tauy = cgstau*tauy
!
! --- calculates wind stress modulus ---
!      rmod=sqrt(taux**2+tauy**2)
!
      return
!     else
!	      print*,'zero wind',w,k
!      endif !DEB
      end
!
!=====================================================================
!
      function CD(speed,delt)
!
! --- calculates the Drag Coefficient as a function of the abs. value of the
! --- wind velocity
! --- ( Hellermann and  Rosenstein )
!
      dimension a(6)
      data a /0.934e-3,0.788e-4,0.868e-4,-0.616e-6,-.120e-5,-.214e-5/
!
      cd = a(1) + a(2)*speed + a(3)*delt + a(4)*speed*speed
     !   + a(5)*delt*delt  + a(6)*speed*delt
!
      return
      end
!
!=====================================================================
!
      double precision function ESK(t)
!
! --- compute the saturation water vapor pressure from
! --- temperature in kelvin  (Lowe,1977; JAM,16,100,1977)
!
      dimension a(7)
      data  a /6984.505294,-188.9039310,2.133357675,-1.288580973e-2,    &
     &    4.393587233e-5,-8.023923082e-8,6.136820929e-11  /
!
      esk = a(1) +t*(a(2)+t*(a(3)+t*(a(4)+t*(a(5)+t*(a(6)+t*a(7))))))
!
      return
      end
!
!======================================================================
!
!======================================================================
!
      double precision function ESC(t)
!
! --- compute the saturation water vapor pressure from
! --- temperature in celsius  (Lowe,1977; JAM,16,100,1977)
!
      dimension a(7)
      data a /6.107799961,4.436518521e-1,1.428945805e-2,        &
     &      2.650648471e-4,3.031240396e-6,2.034080948e-8,6.136820929e-11 /
!
      esc = a(1) +t*(a(2)+t*(a(3)+t*(a(4)+t*(a(5)+t*(a(6)+t*a(7))))))
!
      return
      end
!
!======================================================================
!
      function HEATLAT(t)
!
! --- calculates the Latent Heat of Vaporization ( J/kg ) as function of
! --- the sea surface temperature ( Celsius degrees )
! --- ( from A. Gill  pag. 607 )
!
! --- Constant Latent Heat of Vaporization ( Rosati,Miyakoda 1988 )
!     L = 2.501e+6  (MKS)
!
      double precision t

      heatlat = 2.5008e+6 -2.3e+3*t
!
      return
      end
!
!===================================================================
!
        function qshort(imonth,iday,hour,dlat,dlon,cldnow)
!
        parameter(pi=3.1415927,degrad=pi/180.,degradr=180./pi, eclips=23.439*degrad)
!
        dimension alpham(12),alb1(20),za(20),dza(19)
       
!
! ---   alat,alon - (lat, lon)  in radians !!
!
        data solar/1350./
        data tau /0.7/
        data aozone /0.09/
        data yrdays /360./
!        data yrdays /365./
        data sunalpha /0.03/
!
        data alb1/.719, .656, .603, .480, .385, .300, .250, .193, .164  &
     &  ,.131 , .103, .084, .071, .061, .054, .039, .036, .032, .031 ,.030 /
!
        data za/ 90., 88., 86., 84., 82., 80., 78., 76., 74., 70., 66.  &
     &  ,62., 58., 54., 50., 40., 30., 20., 10., 0.0 /
!
        data dza/8*2.0, 6*4.0, 5*10.0/
!
! --- albedo monthly values from Payne (1972) as means of the values
! --- at 40N and 30N for the Atlantic Ocean ( hence the same latitudinal
! --- band of the Mediterranean Sea ) :
!
        data alpham /0.09,0.08,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.07,0.09,0.10/
!
!----------------- calculations start -------------------------------
!
        alat=dlat*degrad
        alon=dlon*degrad
!        write(*,*)'hour in qshort',hour
        day=iday
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! *** BE CAREFUL: in the MOM version, the simple calendar was used
! *** and the input day was a day of the calendar, not a progressive
! *** day.
!        if(iday.gt.30) day=31.
!
! --- number of days in a year :
!
!        days = 30.*float(imonth-1) + day -1.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   days is the number of days elapsed until the day=iday
        days = day -1.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        th0 = 2.*pi*days/yrdays
        th02 = 2.*th0
        th03 = 3.*th0
!

!
! --- sun declination :
!
        sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0) -     &
     &  0.006758*cos(th02) + 0.000907*sin(th02) - 0.002697*cos(th03) + 0.001480*sin(th03)
!
! --- sun hour angle :
!
!        print*, 'qshort : hour =', hour
! --- 15 are the degrees of longitude the sun covers in 1 hour

        thsun = (hour -12.)*15.*degrad + alon
!
! --- cosine of the solar zenith angle :
!
        coszen =sin(alat)*sin(sundec)+cos(alat)*cos(sundec)*cos(thsun)
        if (coszen .le. 5.035D-04) then
!oddo        if(coszen .le. 0.0) then
          coszen = 0.0
          qatten = 0.0
        else
          qatten = tau**(1./coszen)
        end if
!
        R=(1.0+1.67E-2*cos(pi*2.*(days-3.0)/365.0))**2
        qzer  = coszen * solar *R
        qdir  = qzer * qatten
        qdiff = ((1.-aozone)*qzer - qdir) * 0.5
        qtot  =  qdir + qdiff
!
! --- In the sunbet formula enters days=julian days of the
! --- current year (1,365); we subtract 81 and not 82 because
! --- days = day -1.
! --- conversion of (days - 81) to radians
!
        tjul = (days -81.)*degrad
!
! --- sin of the solar noon altitude in radians :
!
       sunbet=sin(alat)*sin(eclips*sin(tjul)) + cos(alat)*cos(eclips*sin(tjul))
!
! --- solar noon altitude in degrees :
!
      sunbet = degradr*asin(sunbet)
!
!------------------------------------------------------------------------
! --- calculates the albedo as a function of the solar zenith angle :
! --- ( after Payne jas 1972 )
!------------------------------------------------------------------------
!


!------------------------------------------------------------------------
!
! --- solar zenith angle in degrees :
!
        zen=degradr*acos(coszen)
!
        if(zen.ge.74.)then
          jab=.5*(90.-zen)+1.
        elseif(zen.ge.50.)then
          jab=.23*(74.-zen)+9.
        else
          jab=.10*(50.-zen)+15.
        endif
!
        dzen=(za(jab)-zen)/dza(jab)
!
        albedo=alb1(jab)+dzen*(alb1(jab+1)-alb1(jab))
!
!------------------------------------------------------------------------
! --- ( radiation as from Reed(1977), Simpson and Paulson(1979) )
! --- calculates SHORT WAVE FLUX ( watt/m*m )
! --- ( Rosati,Miyakoda 1988 ; eq. 3.8 )
! --- clouds from COADS perpetual data set
!------------------------------------------------------------------------
!
!       qshort = qtot*(1-0.62*cldnow + .0019*sunbet)*(1.-sunalpha)
!       qshort = qtot*(1-0.62*cldnow + 0.0019*sunbet)*(1.-alpham(imonth))
!
!  correction for cloud less than 0.3 according to Reed
      if (cldnow.lt.0.3) then
        qshort = qtot*(1.-albedo)
      else
        qshort  = qtot*(1-0.62*cldnow + .0019*sunbet)*(1.-albedo)
      endif
!
        return
        end

!======================================================================

!-------------------------------------------------------------------
        end module heat_areg
!-------------------------------------------------------------------
