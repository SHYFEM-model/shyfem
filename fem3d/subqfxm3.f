
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c heat flux module (AREG)
c
c contents :
c
c subroutine heatareg(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)
c	computes heat fluxes from formulas in AREG
c subroutine heatpom(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)
c	computes heat fluxes from formulas in AREG (different implementation)
c subroutine vapor1(tt,p,u,e,r,q)
c	computes values for water vapor given relative humidity
c subroutine areg_fluxes(sst,tair,rhnow,speed,cldnow,qsw,qbw,ha,elat,qsurf,evap)
c	this is an adaptation of AREG heat flux routines
c
c revision log :
c
c 27.08.2009    dbf&ggu	integrated into model
c 16.02.2011    ggu	pstd introduced
c
c*******************************************************************

	subroutine heatareg(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)

c computes heat fluxes from formulas in AREG
c
c heat fluxes are positive upward (from sea to atmosphere)

	implicit none

	real t		!air temperature [C]			- in
	real p		!pressure [mb]				- in
	real w		!wind speed [m/s]			- in
	real ur		!relative humidity [%] ([0-100])	- in
	real cc		!cloud cover [0-1]			- in
	real ts		!sea temperature [C]			- in
	real qsens	!sensible heat flux [W/m**2]		- out
	real qlat	!latent heat flux [W/m**2]		- out
	real qlong	!long wave radiation [W/m**2]		- out
	real evap	!evaporation [kg/m**2/s]		- out

	real sst,tair,rhnow,speed,cldnow
        real qsw,qbw,ha,elat,qsurf,evap1

	if( p > 10000 ) stop 'error stop heatareg: p not in mbar'

c qsw,qsurf,evap1 are dummy variables

	sst = ts
	tair = t
	rhnow = ur
	speed = w
	cldnow = cc

	call areg_fluxes(sst,tair,rhnow,speed,cldnow
     +				,qsw,qbw,ha,elat,qsurf,evap1)

	qlong = qbw
	qsens = ha
	qlat = elat
	evap  = qlat / 2.5e+6

	end

c*******************************************************************

	subroutine heatpom(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)

c computes heat fluxes from formulas in AREG (different implementation)
c
c heat fluxes are positive upward (from sea to atmosphere)

	implicit none

	include 'subqfxm.h'

	real t		!air temperature [C]			- in
	real p		!pressure [mb]				- in
	real w		!wind speed [m/s]			- in
	real ur		!relative humidity [%] ([0-100])	- in
	real cc		!cloud cover [0-1]			- in
	real ts		!sea temperature [C]			- in
	real qsens	!sensible heat flux [W/m**2]		- out
	real qlat	!latent heat flux [W/m**2]		- out
	real qlong	!long wave radiation [W/m**2]		- out
	real evap	!evaporation [kg/m**2/s]		- out

c other variables:
c
c       real q          !specific humidity [0-1]
c       real e          !vapor pressure [mb]
c       real r          !mixing ratio [0-1]

	real rdry,lv0,pascal
	parameter(rdry=287.04,lv0=2.5008e6,pascal=100.)

	real rho,cp,lv,tv
	real rhod	!density of moist air DEB
	real pw		!partial pres of water vap in moist air
	real pa		!atm pres of moist air
	real s0,s
	real e,r,q
	real es,rs,qs
	real dt,ch,ce
	real theta,theta2,theta4,thetas
	real cloud
	logical bignami,may

	bignami=.false.
	may=.true.

	if( p > 10000 ) stop 'error stop heatpom: p not in mbar'

c	------------------------------------------------
c	initialization
c	------------------------------------------------

	call vapor1(t,p,ur,e,r,q)        !compute e,r,q

	tv = ( t + kelv ) * ( 1. + 0.6078 * q )
	lv = lv0 - 2.3e3 * ts
	rho = pascal * p / ( rdry * tv )
	rhod =rho * ( 1 + r )/( 1 + 1.609 * r )
	cp = 1004.6 * ( 1. + 0.8375 * q )

	call satur(ts,p,es,rs,qs)	!compute saturation values for sea

c	------------------------------------------------
c	compute Kondo related parameters
c	------------------------------------------------

	dt = ts - t
	
c	here kondo formulation for ch has to be used

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

c	------------------------------------------------
c	sensible and latent heat flux
c	------------------------------------------------

	qsens = rhod * cp *ch * w * dt

	evap = rhod * ce * w * (qs - q)
	qlat = lv * evap

c	------------------------------------------------
c	long wave radiation
c	------------------------------------------------

	theta = t + kelv
	thetas =  ts + kelv
	theta2 = theta * theta
	theta4 = theta2 * theta2
	cloud = 1. - 0.75 * cc**(3.4)

	if(may) then
	  qlong = (bolz * theta4 * (0.4 - 0.05 * sqrt(e) ) + 4 * bolz *
     +	  theta2 * theta *( thetas - theta)) * cloud
	else if(bignami) then
	  qlong=0.98*bolz*thetas**4-bolz*theta4*(0.653+0.00535*e)*
     +	  (1+0.1762*cc*cc)
	else
	  stop 'error stop heatpom: must choose either may or bignami'
     	endif

c	------------------------------------------------
c	end of routine
c	------------------------------------------------

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

c*******************************************************************

	subroutine vapor1(tt,p,u,e,r,q)

c computes values for water vapor given relative humidity

	implicit none

	real tt         !air temperature [C]			- in
	real p          !pressure [mb]				- in
	real u          !relative humidity [%] ([0-100])	- in
	real e          !vapor pressure [mb]			- out
	real r		!mixing ratio [0-1]			- out
	real q		!specific humidity [0-1]		- out

	real ew,rw,qw			!saturation vapor
	real eps,fw			!molecular weight ratio
	parameter(eps=0.62197)
	real fice
	parameter (fice=0.00422)		!over
	real fsalt
	parameter (fsalt=0.98)
	real aux
	
	real a(7)
c	data a / 6984.505294,-188.9039310,2.133357675,-1.288580973e-2,
c     +		4.393587233e-5,-8.023923082e-8,6.136820929e-11  /
     	data a / 6.107799961, 4.436518521e-1, 1.428945805e-2,
     +      2.650648471e-4, 3.031240396e-6, 2.034080948e-8,
     +	    6.136820929e-11 /

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
	
c************************************************************************

c************************************************************************
c ADAPTATION OF AREG2 SUBROUTINE
c       subroutine fluxes(alat,alon,imonth,iday,hour,
c     *  sst,tnow,rhnow,speed,cldnow,qsw,qbw,ha,elat,
c     *  qsurf,taux,tauy,rmod,evap)
c************************************************************************

       subroutine areg_fluxes(sst,tair,rhnow,speed_in,cldnow
     +				,qsw,qbw,ha,elat,qsurf,evap)

c this is an adaptation of AREG heat flux routines

c	sst		sea surface temperature [C]		- in
c	tair		air temperature [C]			- in
c	rhnow		relative humidity [%] ([0-100])		- in
c	speed_in	wind speed [m/s]			- in
c	cldnow		cloud cover [0-1]			- in
c	qsw		short wave radiation [W/m**2]		- out
c	qwb		long wave radiation [W/m**2]		- out
c	ha		sensible heat flux [W/m**2]		- out
c	elat		latent heat flux [W/m**2]		- out
c	qsurf		net surface heat flux [W/m**2]		- out
c	evap		evaporation rate [m/s]			- out
c
c	qwb,ha,elat are positive from sea to atmosphere
c	qsw is positive from atmosphere to sea
c	qsurf is positive from atmosphere to sea (qsurf=qsw-qwb-ha-elat)
c
c	originally qsw is computed inside this subroutine
c	...however, we use value computed outside
c	...therefore qsw is set to 0 inside this subroutine
c
c  1)   Water flux (WFLUX)                 [ watt/m*m ]
c  2)   Short wave flux (QSW)              [ watt/m*m ]
c  3)   Long wave flux backward (QBW)      [ watt/m*m ]
c  4)   Latent heat of evaporation (ELAT)  [ watt/m*m ]
c  5)   Sensible heat flux   (HA)          [ watt/m*m ]
c  6)   Net Upward flux (QU)               [ watt/m*m ]
c  7)   Net surface heat flux (QSURF)      [ watt/m*m ]
c  8)   Wind stress x-component   (TAUX)   [ newton/m*m ]               
c  9)   Wind stress y-component   (TAUY)   [ newton/m*m ]               
c
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c ^                                                                  ^
c ^  All the calculations are made in MKS                            ^
c ^  The unit conversion in CGS is made at the end !                 ^
c ^                                                                  ^
c ^  1 Newton      = 10**5 dyn                                       ^
c ^  1 Joule       = 10**7 erg                                       ^
c ^  1 Newton/m*m  = 10 dyn/cm*cm                                    ^
c ^  1 W/m*m       = 10**3 erg/sec*cm*cm                             ^
c ^  1 W/m*m = 1 J/(m*m*s) = 1/(4.19*10**4) cal/(cm*cm*s)            ^
c ^                                                                  ^
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c bug fixes (by ggu):
c
c	if speed == 0 => speed = 0.01
c	do not change tair -> tnow
c	initialize qsw
c	cavaleri??????
c
c---------------------------------------------------------------------

      include 'subqfxm.h'

      parameter ( pstd = 1013.25 )
      parameter (const = 0.622/pstd)
      parameter (cgs = 1.e-4/4.19)
      parameter (cavaleri=1.00)
c ---- DENSITY OF PURE WATER ---------
      parameter (rho0 = 1000.e0)
c---------------------------------------------------------------------
c     coefficients ( in MKS )  :
c---------------------------------------------------------------------
c
c --- turbulent exchange coefficients
      data  ce,ch  / 1.1e-3, 1.1e-3/    
c  
c --- surface air pressure
      data ps / pstd/
c
c --- air density
      data airden   /1.2/
c     
c --- conversion factor for Langley                                            
c      data  watlan /1.434e-3/       
c
c --- conversion factor for CGS
c      data cgs /1.e+3/      
       data cgstau /10./    
c
       data precip  /0.0/
c	if(speed.gt.0)then !DEB
c
c	print*, 'fluxes : hour =', hour
c
c --- compute speed for bulk coeff.
c                                              
c DEB      speed = sqrt(unow*unow + vnow*vnow)                       
	speed = speed_in
	if( speed .lt. 0.01 ) speed = 0.01
      speedC=speed*cavaleri
c --- precipitation at the moment is supposed to be absent                        
c     precip = 0.0                                                              
c --- SST data converted in Kelvin degrees
c
	ckelv = kelv
      sstk = sst + ckelv 
c
c     TAIR=TNOW is already in Kelvin deg.
c
      tnow = tair + ckelv !DEB
c
c --- calculates the Saturation Vapor Pressure at air temp. and at sea temp.
c ---      esat(Ta) , esat(Ts)
c
      esatair = esk(tnow)
      esatoce = esk(sstk)
c
c --- calculates the saturation mixing ratios at air temp. and sea temp.
c --- wsat(Ta) , wsat(Ts)
c
      wsatair = (const06/ps) * esatair
      wsatoce = (const06/ps) * esatoce
c
c --- calculates the mixing ratio of the air
c --- w(Ta)
c
      wair = 0.01 * rhnow * wsatair
c
c --- calculates the virtual temperature of air
c
      vtnow = (tnow*(const06+wair))/(const06*(1.+wair))
c
c --- calculates the density of the moist air
c
      rhom = 100.*(ps/rgas)/vtnow
c
c ------------------ (i)      short wave
c
	qsw = 0.	!we must give a value to this

c       if (hour .eq. 0.)then
c --- at hour=0 qsw results always 0
c         qsw = 0.
c       else
c DEB           qsw=qshort(imonth,iday,hour,alat,alon,cldnow)
c       endif
c
c ------------------ (ii)      long  wave
c
c --- calculates the term :      Ea**1/2 = [r*Esat(Ta)]**1/2
      ea = sqrt(0.01*rhnow*esatair)
c
c --- calculates the NET LONGWAVE FLUX ( watt/m*m )
c --- WE USE THE qb4 OPTION
c
c --- May formula from May (1986) :
c
       qbw = (1.-0.75*(cldnow**3.4))
     &             * (bolz*(tnow**4)*(0.4 -0.05*ea)
     &             + 4.*bolz*(tnow**3)*(sstk-tnow))
c
c
c ------------------ (iii)   sensible heat
c
c --- calculates the term :      ( Ts - Ta )
c
      deltemp = sstk - tnow
c

c
c --- WE USE THE ceh3 OPTION
c
c --- variable turbulent exchange coefficients ( from Kondo 1975 )
c
c --- calculate S :
c
             s=deltemp/(speed**2.)
c
c --- calculate the Stability Parameter :
c
             stp=s*abs(s)/(abs(s)+0.01)
c
c --- for stable condition :
c
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
c
                endif
                endif
c
c --- for unstable condition :
c
             ELSE
                fh=1.0+0.63*sqrt(stp)
                fe=fh
             ENDIF
c
c --- calculate the coefficient CH,CE,CD
c
             ch=1.3e-03*fh
             ce=1.5e-03*fe
c
c --- calculates the SENSIBLE HEAT FLUX in CGS ( watt/m*m )
c
c      HA = airden*cpw*ch*speed*deltemp
      !ha = rhom*cpw*ch*speed*deltemp
      ha = rhom*cpa*ch*speed*deltemp
c	print*,'ha  ',rhom,'cpw',cpw,'ch ',ch,'sp ',speed,'delt ',deltemp
c
c ------------------ (iv)  latent heat
c
c     const = .622/pstd

c
c --- calculates the term :      esat(Ts)-r*esat(Ta)
c
      evap  = esatoce - rhnow*0.01*esatair
      dsatp = evap
c
c --- calculates the term :      Ce*|V|*[esat(Ts)-r*esat(Ta)]0.622/pstd
c
c --- evap is not the evaporation 
      evap = ce*speed*evap*const

c --- calculates the water flux
c
c      WFLUX =  evap - precip
c
c --- calculates the LATENT HEAT FLUX  LEa ( watt/m*m )
c ---     LEa = L*rho*Ce*|V|*[esat(Ts)-r*esat(Ta)]0.622/pstd
c
c      ELAT = airden*evap*heatlat(sst)
      elat = rhom*evap*heatlat(sst)
c      print*,'dbf elat','rhom ',rhom,'ev ',evap,heatlat(sst),'sst ',sst
c --- evap is now (and ONLY NOW) the EVAPORATION RATE (m/s)
      evap = elat/(rho0*heatlat(sst))
c
c ------------------ net upward flux
c
c --- calculates : Qu = Qb + Ha + LEa  ( Rosati,Miyakoda 1988 ; eq. 3.2 )
c
      qu = qbw + ha + elat
c
c ------------------ net surface heat flux (total budget)
c
c --- calculates : Q = Qs - Qu  ( Rosati,Miyakoda 1988 ; eq. 3.2 )
c
      qsurf = qsw - qu
c --- convert in CGS --
c      QSURF = QSURF * cgs
c
c ------------------ wind stresses
c
c --- calculates the Drag Coefficient
c
c       deltem = -deltemp   ! WRONG

        cdx = cd(speedC,deltemp)
c
c --- calculates the wind stresses in MKS ( newton/m*m )

c ---            taux= rho*Cd*|V|u      tauy= rho*Cd*|V|v
c
c       TAUX= airden*cdx*speed*unow
c       TAUY= airden*cdx*speed*vnow
        
c        taux = rhom*cdx*speedC*unow*cavaleri
c        tauy = rhom*cdx*speedC*vnow*cavaleri
c
c --- convert in CGS ---
c        taux = cgstau*taux
c        tauy = cgstau*tauy
c
c --- calculates wind stress modulus ---
c      rmod=sqrt(taux**2+tauy**2)
c
      return
c     else
c	      print*,'zero wind',w,k
c      endif !DEB
      end
c
c=====================================================================
c
      function CD(speed,delt)
c
c --- calculates the Drag Coefficient as a function of the abs. value of the
c --- wind velocity
c --- ( Hellermann and  Rosenstein )
c
      dimension a(6)
      data a /0.934e-3,0.788e-4,0.868e-4,-0.616e-6,-.120e-5,-.214e-5/
c
      cd = a(1) + a(2)*speed + a(3)*delt + a(4)*speed*speed
     !   + a(5)*delt*delt  + a(6)*speed*delt
c
      return
      end
c
c=====================================================================
c
      real function ESK(t)
c
c --- compute the saturation water vapor pressure from
c --- temperature in kelvin  (Lowe,1977; JAM,16,100,1977)
c
      dimension a(7)
      data  a /6984.505294,-188.9039310,2.133357675,-1.288580973e-2,
     +    4.393587233e-5,-8.023923082e-8,6.136820929e-11  /
c
      esk = a(1) +t*(a(2)+t*(a(3)+t*(a(4)+t*(a(5)+t*(a(6)+t*a(7))))))
c
      return
      end
c
c======================================================================
c
c======================================================================
c
      real function ESC(t)
c
c --- compute the saturation water vapor pressure from
c --- temperature in celsius  (Lowe,1977; JAM,16,100,1977)
c
      dimension a(7)
      data a /6.107799961,4.436518521e-1,1.428945805e-2,
     +        2.650648471e-4,3.031240396e-6,2.034080948e-8,
     +        6.136820929e-11 /
c
      esc = a(1) +t*(a(2)+t*(a(3)+t*(a(4)+t*(a(5)+t*(a(6)+t*a(7))))))
c
      return
      end
c
c======================================================================
c
      function HEATLAT(t)
c
c --- calculates the Latent Heat of Vaporization ( J/kg ) as function of
c --- the sea surface temperature ( Celsius degrees )
c --- ( from A. Gill  pag. 607 )
c
c --- Constant Latent Heat of Vaporization ( Rosati,Miyakoda 1988 )
c     L = 2.501e+6  (MKS)
c
      heatlat = 2.5008e+6 -2.3e+3*t
c
      return
      end
c
c===================================================================
c
        function qshort(imonth,iday,hour,dlat,dlon,cldnow)
c
        parameter(pi=3.1415927,degrad=pi/180.,
     *  degradr=180./pi, eclips=23.439*degrad)
c
        dimension alpham(12),alb1(20),za(20),dza(19)
       
c
c ---   alat,alon - (lat, lon)  in radians !!
c
        data solar/1350./
        data tau /0.7/
        data aozone /0.09/
        data yrdays /360./
c        data yrdays /365./
        data sunalpha /0.03/
c
        data alb1/.719, .656, .603, .480, .385, .300, .250, .193, .164
     $  ,.131 , .103, .084, .071, .061, .054, .039, .036, .032, .031
     $  ,.030 /
c
        data za/ 90., 88., 86., 84., 82., 80., 78., 76., 74., 70., 66.
     $  ,62., 58., 54., 50., 40., 30., 20., 10., 0.0 /
c
        data dza/8*2.0, 6*4.0, 5*10.0/
c
c --- albedo monthly values from Payne (1972) as means of the values
c --- at 40N and 30N for the Atlantic Ocean ( hence the same latitudinal
c --- band of the Mediterranean Sea ) :
c
        data alpham /0.09,0.08,0.06,0.06,0.06,0.06,0.06,0.06,
     &               0.06,0.07,0.09,0.10/
c
c----------------- calculations start -------------------------------
c
        alat=dlat*degrad
        alon=dlon*degrad
c        write(*,*)'hour in qshort',hour
        day=iday
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c *** BE CAREFUL: in the MOM version, the simple calendar was used
c *** and the input day was a day of the calendar, not a progressive
c *** day.
c        if(iday.gt.30) day=31.
c
c --- number of days in a year :
c
c        days = 30.*float(imonth-1) + day -1.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   days is the number of days elapsed until the day=iday
        days = day -1.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        th0 = 2.*pi*days/yrdays
        th02 = 2.*th0
        th03 = 3.*th0
c

c
c --- sun declination :
c
        sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0) -
     *  0.006758*cos(th02) + 0.000907*sin(th02) -
     *  0.002697*cos(th03) + 0.001480*sin(th03)
c
c --- sun hour angle :
c
c        print*, 'qshort : hour =', hour
c --- 15 are the degrees of longitude the sun covers in 1 hour

        thsun = (hour -12.)*15.*degrad + alon
c
c --- cosine of the solar zenith angle :
c
        coszen =sin(alat)*sin(sundec)+cos(alat)*cos(sundec)*cos(thsun)
        if (coszen .le. 5.035D-04) then
Coddo        if(coszen .le. 0.0) then
          coszen = 0.0
          qatten = 0.0
        else
          qatten = tau**(1./coszen)
        end if
c
        R=(1.0+1.67E-2*cos(pi*2.*(days-3.0)/365.0))**2
        qzer  = coszen * solar *R
        qdir  = qzer * qatten
        qdiff = ((1.-aozone)*qzer - qdir) * 0.5
        qtot  =  qdir + qdiff
c
c --- In the sunbet formula enters days=julian days of the
c --- current year (1,365); we subtract 81 and not 82 because
c --- days = day -1.
c --- conversion of (days - 81) to radians
c
        tjul = (days -81.)*degrad
c
c --- sin of the solar noon altitude in radians :
c
       sunbet=sin(alat)*sin(eclips*sin(tjul)) +
     *        cos(alat)*cos(eclips*sin(tjul))
c
c --- solar noon altitude in degrees :
c
      sunbet = degradr*asin(sunbet)
c
c------------------------------------------------------------------------
c --- calculates the albedo as a function of the solar zenith angle :
c --- ( after Payne jas 1972 )
c------------------------------------------------------------------------
c


c------------------------------------------------------------------------
c
c --- solar zenith angle in degrees :
c
        zen=degradr*acos(coszen)
c
        if(zen.ge.74.)then
          jab=.5*(90.-zen)+1.
        elseif(zen.ge.50.)then
          jab=.23*(74.-zen)+9.
        else
          jab=.10*(50.-zen)+15.
        endif
c
        dzen=(za(jab)-zen)/dza(jab)
c
        albedo=alb1(jab)+dzen*(alb1(jab+1)-alb1(jab))
c
c------------------------------------------------------------------------
c --- ( radiation as from Reed(1977), Simpson and Paulson(1979) )
c --- calculates SHORT WAVE FLUX ( watt/m*m )
c --- ( Rosati,Miyakoda 1988 ; eq. 3.8 )
c --- clouds from COADS perpetual data set
c------------------------------------------------------------------------
c
c       qshort = qtot*(1-0.62*cldnow + .0019*sunbet)*(1.-sunalpha)
c       qshort = qtot*(1-0.62*cldnow + 0.0019*sunbet)*(1.-alpham(imonth))
c
c  correction for cloud less than 0.3 according to Reed
      if (cldnow.lt.0.3) then
        qshort = qtot*(1.-albedo)
      else
        qshort  = qtot*(1-0.62*cldnow + .0019*sunbet)*(1.-albedo)
      endif
c
        return
        end

c======================================================================

