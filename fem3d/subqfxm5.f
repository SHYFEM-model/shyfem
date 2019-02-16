
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

! heat flux module (COARE 2.5 and 3.0)
! It also computes momentum fluxes due to wind and rainfall.
!
! contents :
!
! subroutine heatcoare(airt,airp,ws,rh,cloud,sst,prec,rad,
!     &                  qsens,qlat,qlong,evap)
!       handles heat fluxes from COARE heat module
! subroutine humidity(hum_method,hum,pres,tw,ta,es,ea,qs,qa,rhoa)
!	computes saturation vapour pressure and specific humidity
! function psi(iflag, ZoL)
!	evaluates the stability function
! subroutine back_radiation(method,tw_c,ta_c,ea,qa,cloud,qb)
!	calculates the long-wave back radiation
! subroutine coare25(tw,ta,ws,rain,qs,qa,rhoa,evp,qe,qh,cd)
!       computes heat fluxes according COARE2.5 formulation
! subroutine coare30(tw,ta,ws,rain,rad,-qb,qs,qa,rhoa,evp,qe,qh,cd)
!       computes heat fluxes according COARE3.0 formulation
! subroutine wcstress(ws,rhoa,rain,wx,wy,taux,tauy)
!	computes momentum fluxes due to wind and rain
! subroutine tw_skin(rad,qrad,tw,hb,usw,dt,dtw,tws)
!	computes sea surface skin temperature
!
! revision log :
!
! 11.06.2014    ccf     heat module from COARE integrated
! 22.01.2016    ggu     bug fix for COARE integrated (L gets infinite)
! 05.10.2018    ggu     do not compute Ch and Ce (not needed, divide by 0)
!
!***********************************************************************

	subroutine heatcoare(airt,airp,ws,rh,cloud,sst,prec,rad,
     +			qsens,qlat,qlong,evap,cd)

	implicit none

	include 'subqfxm.h'

! input
        real, intent(in)    :: airt       !air temperature [C] 
        real, intent(in)    :: airp       !pressure [mb]               
        real, intent(in)    :: ws         !wind speed [m/s] 
        real, intent(in)    :: rh         !relative humidity [%] ([0-100])
        real, intent(in)    :: cloud      !cloud cover [0-1] 
        real, intent(in)    :: sst        !sea temperature [C]                 
        real, intent(in)    :: prec       !precipitation [m/s]
	real, intent(in)    :: rad	  !net solar flux (w/m^2)
! inout
! output
        real,intent(out)    :: qsens      !sensible heat flux [W/m**2]
        real,intent(out)    :: qlat       !latent heat flux [W/m**2] 
        real,intent(out)    :: qlong      !long wave radiation [W/m**2] 
        real,intent(out)    :: evap       !evaporation [kg/m**2/s] 
	real,intent(out)    :: cd	  !drag_coefficient
! local
	integer, parameter  :: lon_method=2
	real 		    :: ta,tw,rain
	real 		    :: es,ea,qs,qa,rhoa
	real 		    :: qe,qh,qb

        if( airp > 10000 ) stop 'error stop heatcoare: p not in mbar'

!       ---------------------------------------------------------
!       Initialize
!       ---------------------------------------------------------

	ta   = airt
	tw   = sst
	rain = prec*rhow	!from m/s to kg/m2/s

!       ---------------------------------------------------------
!       Compute saturation vapour pressure and specific humidity
!       ---------------------------------------------------------

	call humidity(rh,airp,tw,ta,es,ea,qs,qa,rhoa)

!       ---------------------------------------------------------
!       Compute long wave radiation flux
!       ---------------------------------------------------------

	call back_radiation(lon_method,tw,ta,ea,qa,cloud,qb)

	qlong = - qb

!       ---------------------------------------------------------
!       Compute sensible and latent heat fluxes and dragg coeff.
!       ---------------------------------------------------------

	!call coare25(tw,ta,ws,rain,qs,qa,rhoa,evap,qe,qh,cd)
	call coare30(tw,ta,ws,rain,rad,-qb,qs,qa,rhoa,evap,qe,qh,cd)

	qsens = - qe
	qlat  = - qh

!       ---------------------------------------------------------
!       End of routine
!       ---------------------------------------------------------

	end subroutine 

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! ADAPTED FROM GOTM CODE
!
! This routine calculated the saturation vapour pressure at SST and at
! air temperature, as well as the saturation specific humidty and the
! specific humidity.
!
	subroutine humidity(hum,airp,tw,ta,es,ea,qs,qa,rhoa)

	implicit none

	include 'subqfxm.h'

! input
	real, intent(in)        :: hum		!humidity [%]
	real, intent(in)        :: airp		!air pressure [mb]
	real, intent(in)        :: tw		!water temp [C]
	real, intent(in)        :: ta		!air temp [C]
! output
	real, intent(out)       :: qs		!bulk water spec hum (kg/kg)
	real, intent(out)	:: qa		!bulk air spec hum (kg/kg)
	real, intent(out)       :: rhoa		!moist air density (kg/m3)
	real, intent(out)       :: es
	real, intent(out)       :: ea
! local
	real        		:: rh,pres

!	------------------------------------------------
!       Initialize
!	------------------------------------------------
	pres = airp * 100.	!Conversion mb --> Pascal

!	------------------------------------------------
!       Saturation vapor pressure - using SST and AIRP
!	------------------------------------------------
	call saturation(tw,airp,es)
	es = es * 100.0		! Conversion mb --> Pascal

!	------------------------------------------------
!       Correction for seawater, following Kraus 1972
!       Correcting for salt water assuming 98% RH
!	------------------------------------------------
	es = 0.98 * es

!	------------------------------------------------
!       Saturation specific humidity
!	------------------------------------------------
	qs = const06*es/(pres-0.378*es)

	rh = 0.01 * hum

!	------------------------------------------------
!       Saturation vapor pressure at that air temperature
!	------------------------------------------------
	call saturation(ta,airp,ea)
        ea = ea * 100.0 	! Conversion mb --> Pascal

!	------------------------------------------------
!       Get actual vapor pressure and convert to specific humidity
!	------------------------------------------------
	ea = rh * ea
	qa = const06*ea/(pres-0.378*ea)

!	------------------------------------------------
!       Compute moist air density (kg/m3)
!	------------------------------------------------
	rhoa = pres/(rgas*(ta+kelv)*(1.0+const06*qa))

	return

	end subroutine humidity

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Tetens' formula for saturation vp Buck(1981) JAM 20, 1527-1532

	subroutine saturation(t,p,esat)

	implicit none
     
! input
	real, intent(in)	:: t		!temperature [C]
	real, intent(in)	:: p		!air pressure [mb]
! output
	real, intent(out)	:: esat 	!saturation vp [mb]

	esat = 6.112*exp(17.502*t/(t+240.97))*.98*(1.0007+3.46e-6*p)
	!esat = (1.0007+3.46e-6*p)*6.1121*exp(17.502*t/(240.97+t))

	return

	end subroutine saturation

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! ADAPTED FROM GOTM CODE                                               !
!								       !
!  This function evaluates the stability function, PSI, for wind       !
!  speed (iflag=1) or for air temperature and moisture (iflag=2)       !
!  profiles as function of the stability parameter, ZoL (z/L).         !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Liu, W.T., K.B. Katsaros, and J.A. Businger, 1979:  Bulk          !
!        parameterization of the air-sea exchange of heat and          !
!        water vapor including the molecular constraints at            !
!        the interface, J. Atmos. Sci, 36, 1722-1735.                  !

	function psi(iflag, ZoL)

	implicit none

	include 'subqfxm.h'

	real  :: psi
! input
	integer, intent(in) :: iflag
	real, intent(in)    :: ZoL
! parameters
	real, parameter     :: r3 = 1.0/3.0
	real, parameter     :: sqr3 = 1.7320508
! local
	real                :: Fw,chic,chik,psic,psik,c

!	------------------------------------------------
!       Initialize for the zero "ZoL" case.
!	------------------------------------------------

	psi = 0.0
	psik = 0.0
	psic = 0.0

!	------------------------------------------------
!       Unstable conditions.
!	------------------------------------------------

	if (ZoL .lt. 0.0) then
	  chik=(1.0-16.0*ZoL)**0.25
	  if (iflag .eq. 1) then
	     psik=2.0*LOG(0.5*(1.0+chik))+LOG(0.5*(1.0+chik*chik))-
     &	          2.0*ATAN(chik)+ 0.5*pi
	  else if (iflag .eq. 2) then
	        psik=2.0*LOG(0.5*(1.0+chik*chik))
	  end if
!
!  	  ------------------------------------------------
!         For very unstable conditions, use free-convection (Fairall).
!	  ------------------------------------------------
!
	  chic=(1.0-12.87*ZoL)**r3
	  psic=1.5*LOG(r3*(1.0+chic+chic*chic))-
     &	        sqr3*ATAN((1.0+2.0*chic)/sqr3)+ pi/sqr3
!
!	  ------------------------------------------------
!         Match Kansas and free-convection forms with weighting Fw.
!	  ------------------------------------------------
!
	  Fw=1.0/(1.0+ZoL*ZoL)
	  psi=Fw*psik+(1.0-Fw)*psic
!
!	------------------------------------------------
!       Stable conditions.
!	------------------------------------------------
!
	else if (ZoL .gt. 0.0) then
          c=min(50.,.35*ZoL)
	  if (iflag .eq. 1) then
            psi=-((1+1.0*ZoL)**1.0+.667*(ZoL-14.28)/exp(c)+8.525)
	  else if (iflag .eq. 2) then
            psi=-((1.+2./3.*ZoL)**1.5+.6667*(ZoL-14.28)/exp(c)+8.525)
	  end if
	  !psi=-4.7*ZoL
	end if

	return

	end function psi

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! ADAPTED FROM GOTM CODE
!
! Calculate the long-wave back radiation \label{sec:back-rad}
! with one of the four methos:
!   method = 1: Clarketal74
!   method = 2: HastenrathLamb78
!   method = 3: Bignamietal95
!   method = 4: BerliandBerliand52
! It should ne noted that the latitude must here be given in degrees.
!
	subroutine back_radiation(method,tw_c,ta_c,ea,qa,cloud,qb)

	implicit none

	include 'subqfxm.h'

! input
	integer, intent(in)             :: method
	real, intent(in)                :: tw_c,ta_c,ea,qa,cloud
! output
	real, intent(out)               :: qb
! parameters
	integer, parameter   :: clark=1      ! Clark et al, 1974
	integer, parameter   :: hastenrath=2 ! Hastenrath and Lamb, 1978
	integer, parameter   :: bignami=3    ! Bignami et al., 1995 - Medsea
	integer, parameter   :: berliand=4   ! Berliand and Berliand, 1952 - ROMS

	real, parameter, dimension(91)  :: cloud_correction_factor = (/ 
     &    0.497202,  0.501885,  0.506568,  0.511250,  0.515933, 
     &    0.520616,  0.525299,  0.529982,  0.534665,  0.539348, 
     &    0.544031,  0.548714,  0.553397,  0.558080,  0.562763, 
     &    0.567446,  0.572129,  0.576812,  0.581495,  0.586178, 
     &    0.590861,  0.595544,  0.600227,  0.604910,  0.609593, 
     &    0.614276,  0.618959,  0.623641,  0.628324,  0.633007, 
     &    0.637690,  0.642373,  0.647056,  0.651739,  0.656422, 
     &    0.661105,  0.665788,  0.670471,  0.675154,  0.679837, 
     &    0.684520,  0.689203,  0.693886,  0.698569,  0.703252, 
     &    0.707935,  0.712618,  0.717301,  0.721984,  0.726667, 
     &    0.731350,  0.736032,  0.740715,  0.745398,  0.750081, 
     &    0.754764,  0.759447,  0.764130,  0.768813,  0.773496, 
     &    0.778179,  0.782862,  0.787545,  0.792228,  0.796911, 
     &    0.801594,  0.806277,  0.810960,  0.815643,  0.820326, 
     &    0.825009,  0.829692,  0.834375,  0.839058,  0.843741, 
     &    0.848423,  0.853106,  0.857789,  0.862472,  0.867155, 
     &    0.871838,  0.876521,  0.881204,  0.885887,  0.890570, 
     &    0.895253,  0.899936,  0.904619,  0.909302,  0.913985, 
     &    0.918668 /)

! local
	real                  :: ccf
	real                  :: x1,x2,x3
	real 		      :: lat
	real 		      :: ta, tw

!	------------------------------------------------
!       Initializate
!	------------------------------------------------
	lat = 40.
	ta = ta_c + kelv
	tw = tw_c + kelv

!	------------------------------------------------
!       Calculate cloud correction factor,fortran counts from 1 !
!	Not implemented - CCF 
!	------------------------------------------------
	ccf= cloud_correction_factor(nint(abs(lat))+1)

	select case(method)

	  case(clark)
!	    ------------------------------------------------
!	    Clark et al. (1974) formula.
!	    unit of ea is Pascal, must hPa
!	    Black body defect term, clouds, water vapor correction
!	    ------------------------------------------------
	    x1 = (1.0-ccf*cloud*cloud)*(tw**4)
	    x2 = (0.39-0.05*sqrt(ea*0.01))
!	    temperature jump term
	    x3 = 4.0*(tw**3)*(tw-ta)
	    qb = -emiss*bolz*(x1*x2+x3)

	  case(hastenrath) 
!  	    ------------------------------------------------
!	    Hastenrath and Lamb (1978) formula.
!           qa in g(water)/kg(wet air)
!	    ------------------------------------------------
	    x1 = (1.0-ccf*cloud*cloud)*(tw**4)
	    x2 = (0.39-0.056*sqrt(1000.0*qa))
	    x3 = 4.0*(tw**3)*(tw-ta)
	    qb = -emiss*bolz*(x1*x2+x3)

	  case(bignami)
!	    ------------------------------------------------
!	    Bignami et al. (1995) formula (Med Sea).
!	    unit of ea is Pascal, must hPa
!	    ------------------------------------------------
	     ccf = 0.1762
	     x1 = (1.0+ccf*cloud*cloud)*ta**4
	     x2 = (0.653+0.00535*(ea*0.01))
	     x3 = emiss*(tw**4)
	     qb = -bolz*(-x1*x2+x3)

	  case(berliand)
!  	    ------------------------------------------------
!	    Berliand & Berliand (1952) formula (ROMS).
!	    ------------------------------------------------
	    x1 = (1.0-0.6823*cloud*cloud)*ta**4
	    x2 = (0.39-0.05*sqrt(0.01*ea))
	    x3 = 4.0*ta**3*(tw-ta)
	    qb = -emiss*bolz*(x1*x2+x3)

	  case default
	end select

	return
	end subroutine back_radiation

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! ADAPTED FROM GOTM CODE
!
! Heat and momentum fluxes according to Fairall et al 96a.
! Cool skin and warm layer effects are considered according to the
! suggestions of \cite{Fairalletal96b}.
! This piece of code has been adapted from the COARE code originally
! written by David Rutgers and Frank Bradley - see
! http://www.coaps.fsu.edu/COARE/flux_algor/flux.html.

	subroutine coare25(sst,airt,w,rain,qs,qa,rhoa,evap,qe,qh,Cd)

	implicit none

	include 'subqfxm.h'

! input
	real, intent(in)      :: sst		!water temperature [C]
	real, intent(in)      :: airt		!air temperature [C]
	real, intent(in)      :: w		!wind speed [m/s]
	real, intent(in)      :: rain  		!precipitation [kg/m2/s] 
        real, intent(in)      :: qs		!bulk water spec hum [kg/kg]
        real, intent(in)      :: qa		!bulk air spec hum [kg/kg]
        real, intent(in)      :: rhoa		!moist air density [kg/m3]
! output
	real, intent(out)     :: evap		!evaporation [kg/(m**2 s)]
	real, intent(out)     :: qe		!sensible heat flux [W/m**2]
	real, intent(out)     :: qh		!latent heat flux [W/m**2]
	real, intent(out)     :: Cd		!wind dragg coeff.

! parameters
!  Fairall LKB roughness Reynolds number to Von Karman
	real, parameter       :: fdg = 1.0          ! non-dimensional

!  Beta parameter evaluated from Fairall low windspeed turbulence data.
	real, parameter       :: beta = 1.2         ! non-dimensional

!  Zabl      Height (m) of atmospheric boundary layer.
	real, parameter       :: Zabl = 600.0       ! in [m]
	real, parameter       :: r3 = 1.0/3.0
!
!  Liu et al. (1979) look-up table coefficients to compute roughness
!  Reynolds number for temperature (rt) and moisture (rq) as a
!  function of wind Reynolds number (rr):
!
!       rt = Liu_a(:,1) * Rr   ** Liu_b(:,1)    temperature
!       rq = Liu_a(:,2) * Rr   ** Liu_b(:,2)    moisture
!
	real,parameter, dimension(8,2) :: Liu_a = reshape ( 
     &           (/ 0.177,  1.376,    1.026,      1.625,   
     &              4.661, 34.904, 1667.190, 588000.0,     
     &              0.292,  1.808,    1.393,      1.956,   
     &              4.994, 30.709, 1448.680, 298000.0 /),  
     &           (/ 8, 2 /) )

	real,parameter, dimension(8,2) :: Liu_b = reshape (
     &           (/  0.0,    0.929, -0.599, -1.018,        
     &              -1.475, -2.067, -2.907, -3.935,        
     &               0.0,    0.826, -0.528, -0.870,        
     &              -1.297, -1.845, -2.682, -3.616 /),     
     &           (/ 8, 2 /) )

	real,parameter, dimension(9) :: Liu_Rr = 
     &           (/    0.0,  0.11,   0.825,   3.0,         
     &                10.0, 30.0,  100.0,   300.0,        
     &              1000.0 /)

        real, parameter	      :: charn= 0.011        !Charnock
	real, parameter       :: wgust= 0.2
	real, parameter       :: zt = 2.0  !Height (m) of air temp. measurement.
	real, parameter       :: zq = 2.0  !Height (m) of air hum. measurement
	real, parameter       :: zw = 10.0 !Height (m) of winds measurement
	integer,  parameter   :: itermax = 20
! local
	real        :: cff,wgus
	real        :: L
	real        :: ta,ta_k,tw
	integer     :: iter,k
	real        :: vis_air
	real        :: tpsi,qpsi,wpsi,ZWoL,oL,ZToL,ZQoL,ZoW,ZoT, ZoQ
	real        :: Wstar,Tstar, Qstar, delQ, delT, rr,rt,rq
	real        :: TVstar,Bf, upvel,delw,Wspeed
	real        :: cd_rain
	real        :: x1,x2,x3
! function
	real        :: psi

!	------------------------------------------------
!	Initialize
!	------------------------------------------------

	cd   = 2.5e-3
	evap = 0.
	qe   = 0.
	qh   = 0.
	tw   = sst
	ta_k = airt + kelv
	ta   = airt
	delw = sqrt(w*w + wgust*wgust)

!	------------------------------------------------
!	Compute Monin-Obukhov similarity parameters for wind (Wstar),
!	heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!	------------------------------------------------

!	------------------------------------------------
!	Kinematic viscosity of dry air (m2/s), Andreas (1989).
!	------------------------------------------------
	vis_air=1.326e-5*(1.0+ta*(6.542e-3+ta*
     &	      (8.301e-6-4.84e-9*ta)))

!	------------------------------------------------
!	Compute latent heat of vaporization (J/kg) at sea surface
!	------------------------------------------------
	L = (2.501-0.00237*tw)*1.e6
!
!	------------------------------------------------
!	Assume that wind is measured relative to sea surface and include
!	gustiness.
!	Initialize.
!	------------------------------------------------
	delq=qa-qs
	delt=ta-tw

!	------------------------------------------------
!	Initial guesses for Monin-Obukhov similarity scales.
!	------------------------------------------------
	ZWoL=0.0
	ZoW=0.0005
	Wstar=0.04*delw
	Tstar=0.04*delt
	Qstar=0.04*delq
	TVstar=Tstar*(1.0+0.61*qa)+0.61*ta_k*Qstar
	rt = 0.
	rq = 0.

!	------------------------------------------------
!	Iterate until convergence.  It usually converges 
!	within four iterations.
!	------------------------------------------------
	do iter=1,itermax

!	  ------------------------------------------------
!	  Compute Monin-Obukhov stability parameter, Z/L.
!	  ------------------------------------------------
	  oL = grav*von*TVstar/(ta_k*(1.+0.61*qa)*Wstar*Wstar)
	  ZWoL = zw*oL
	  ZToL = zt*oL
	  ZQoL = zq*oL

!	  ------------------------------------------------
!	  Evaluate stability functions at Z/L.
!	  ------------------------------------------------
	  wpsi = psi(1,ZWoL)
	  tpsi = psi(2,ZToL)
	  qpsi = psi(2,ZQoL)

!	  ------------------------------------------------
!	  Compute wind scaling parameters, Wstar.
!	  ------------------------------------------------
	  ZoW=charn*Wstar*Wstar/grav+0.11*vis_air/Wstar
	  Wstar=delw*von/(log(zw/ZoW)-wpsi)

!	  ------------------------------------------------
!	  Computes roughness Reynolds number for wind (Rr), heat (Rt),
!	  and moisture (Rq). Use Liu et al. (1976) look-up table to
!	  compute "Rt" and "Rq" as function of "Rr".
!	  ------------------------------------------------
	  rr=ZoW*Wstar/vis_air
	  rr=max(min(1000.0,rr),0.0)

	  do k=1,8
	    if ((liu_rr(k).le.rr) .and. (rr.lt.liu_rr(k+1))) then
	      rt=liu_a(k,1)*rr**liu_b(k,1)
	      rq=liu_a(k,2)*rr**liu_b(k,2)
	    end if
	  end do

!	  ------------------------------------------------
!         Compute heat and moisture scaling parameters,
!         Tstar and Qstar.
!	  ------------------------------------------------
	  cff=vis_air/Wstar
	  ZoT=rt*cff
	  ZoQ=rq*cff
	  cff=von*fdg
	  Tstar=(delt)*cff/(log(zt/ZoT)-tpsi)
	  Qstar=(delq)*cff/(log(zq/ZoQ)-qpsi)

!	  ------------------------------------------------
!         Compute gustiness in wind speed.
!	  ------------------------------------------------
	  TVstar=Tstar*(1.0+0.61*qa)+0.61*ta_k*Qstar
	  bf=-grav/ta_k*Wstar*TVstar
	  if (bf .gt. 0) then
	     wgus = beta*(bf*Zabl)**r3
	  else
	     wgus = 0.
	  end if
	  delw=sqrt(w*w+wgus*wgus)

	end do

!	------------------------------------------------
!	Compute transfer coefficients for momentun (Cd), heat (Ch),
!	and moisture (Ce).
!	------------------------------------------------
	Wspeed=sqrt(w*w+wgus*wgus)
	Cd=Wstar*Wstar/(Wspeed*Wspeed)

!	------------------------------------------------
!	Compute turbulent sensible heat flux (W/m2), qe.
!	out of ocean is negative
!	------------------------------------------------
	qe=cpa*rhoa*Wstar*Tstar

!	------------------------------------------------
!	Compute sensible heatflux due to rain fall 
!	(units of qs and qa - should be kg/kg)
!	------------------------------------------------
	x1 = 2.11e-5*(ta_k/kelv)**1.94
	x2 = 0.02411*(1.0+ta*(3.309e-3-1.44e-6*ta))/(rhoa*cpa)
	x3 = qa * L /(rgas * ta_K * ta_K)
	cd_rain = 1.0/(1.0+const06*(x3*L*x1)/(cpa*x2))
	cd_rain = cd_rain*cpw*((tw-ta) + (qs-qa)*L/cpa)
	qe = qe - rain * cd_rain

!	------------------------------------------------
!	Compute turbulent latent heat flux (W/m2), qh.
!	------------------------------------------------
	qh=L*rhoa*Wstar*Qstar

!	------------------------------------------------
!	Compute Webb correction (Webb effect) to latent heat flux
!	------------------------------------------------
	upvel=-1.61*Wstar*Qstar-(1.0+1.61*qa)*Wstar*Tstar/ta_k
	qh=qh-rhoa*L*upvel*qa

!	------------------------------------------------
!	Calculation of evaporation/condensation in kg/(m**2 s)
!	------------------------------------------------
	evap = rhoa*Wstar*Qstar

	return

	end subroutine coare25

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Heat and momentum fluxes according to COARE 3.0
! This piece of code has been adapted from the COARE 3.0 code
! retrieved at https://coaps.fsu.edu/COARE/flux_algor/
! jcool = 0, since cool skin and warm layer are already 
! considered by computing SST from model temperature in the 
! first layer (see tw_skin)

	subroutine coare30(sst,airt,ws,rain,Rns,Rnl,Qs,Q,rhoa,
     +			   evap,qsens,qlat,Cd)

	implicit none

	include 'subqfxm.h'

! input
        real, intent(in) 	:: airt         !air temperature [C]
        real, intent(in) 	:: sst          !sea surface temperature [C]
        real, intent(in) 	:: ws           !wind speed [m/s]                    
        real, intent(in) 	:: rain         !precipitation [kg/m2/s]
	real, intent(in) 	:: Rns	        !net solar flux (W/m^2)
	real, intent(in) 	:: Rnl	        !net longwave flux (W/m^2)
	real, intent(in) 	:: Qs,Q	        !Spec. humidity	[kg/kg]		   
	real, intent(in) 	:: rhoa	        !moist air density (kg/m3)	   
! output
        real, intent(out) 	:: qsens        !sensible heat flux [W/m**2]         
        real, intent(out) 	:: qlat		!latent heat flux [W/m**2]           
        real, intent(out) 	:: evap		!evaporation [kg/m**2/s]             
	real, intent(out)	:: Cd		!dragg coefficient		   
! parameters
	real, parameter		:: Beta = 1.2 
	real, parameter		:: fdg  = 1.00 
	real, parameter		:: be   = 0.026 
	real, parameter		:: tcw  = 0.6 
	real, parameter		:: zu   = 10.  !wind speed measurement height (m)
	real, parameter		:: zt   = 2.   !air T measurement height (m)
	real, parameter		:: zq   = 2.   !air q measurement height (m)
	real, parameter		:: zi   = 600. !PBL depth (m)
! local
	real 	:: us		!surface current speed in the wind direction (m/s)
	integer	:: jcool	!implement cool calculation skin switch, 0=no, 1=yes
	integer	:: jwave	!implement wave dependent roughness model
	real 	:: twave, hwave	!wave period (s) and height (m)
	integer	:: i,nits	!iteration
	real 	:: Le,visa	!air constants
	real 	:: Al,bigc,wetc	!cool skin constants
	real 	:: lwave,cwave	!wave parameters
	real 	:: du,dt,dq,ta,ug,dter,dqer,ut
	real 	:: u10,usr,zo10,Cd10,Ch10,Ct10,zot,Ct,CC
	real 	:: Ribcu,Ribu,zetu,L10
	real 	:: tsr,qsr,tkt,charn
	real 	:: zet,L,zoq,Bf
	real 	:: zot10,zo,rr
	real 	:: hsb,hlb,qout,dels,qcol,alq,xlamx
	real 	:: tau,dwat,dtmp,alfac,RF
	real 	:: Ch,Ce
	real	:: aux1,aux2,aux
	character*20 :: aline

	integer :: iunit
	integer :: itact
	real	:: zetaux

! function
	real    :: psi

!	------------------------------------------------
!	Initialize
!	------------------------------------------------

	evap  = 0.
	qsens = 0.
	qlat  = 0.
	us    = 0.
	jcool = 0		!implement cool calculation skin switch, 0=no, 1=yes
	cd    = 2.5e-3

	jwave = 0 		!implement wave dependent roughness model
	twave = 0. 		!wave period (s)
	hwave = 0. 		!wave height (m)

!	------------------------------------------------
!	Set parameters
!	------------------------------------------------

	lwave = grav/2/pi*twave**2 
	cwave = grav/2/pi*twave 

	Le   = (2.501-.00237*sst)*1.e6 
	visa = 1.326e-5*(1+6.542e-3*airt+8.301e-6*airt*airt - 
     +		4.84e-9*airt*airt*airt) 
	aux = sst+3.2
	if( aux < 0. ) aux = 0.
	Al   = 2.1e-5*(aux)**0.79 	!GGUZ0 FIXME
	bigc = 16.*grav*cpw*(rhow*visw)**3/(tcw*tcw*rhoa*rhoa)
	wetc = 0.622*Le*Qs/(rgas*(sst+kelv)**2) 
	
	!Rns = Rs*.945		!oceanic albedo=0.055 daily average
	!Rnl = 0.97*(5.67e-8*(sst-0.3*jcool+kelv)**4-Rl) 
	
!	------------------------------------------------
!	Begin bulk loop
!	------------------------------------------------
!	------------------------------------------------
!	first guess
!	------------------------------------------------

	du = ws-us 
	dt = sst-airt-.0098*zt 
	dq = Qs-Q 
	ta = airt+kelv 
	ug = .5 
	dter = 0.3  
	dqer = wetc*dter 
	ut = sqrt(du*du+ug*ug) 
	u10 = ut*log(10/1e-4)/log(zu/1e-4) 
	usr = .035*u10 
	zo10 = 0.011*usr*usr/grav+0.11*visa/usr 
	Cd10 = (von/log(10/zo10))**2 
	Ch10 = 0.00115 
	Ct10 = Ch10/sqrt(Cd10) 
	zot10 = 10/exp(von/Ct10) 
	Cd = (von/log(zu/zo10))**2 
	Ct = von/log(zt/zot10) 
	CC = von*Ct/Cd 
	Ribcu = -zu/zi/.004/Beta**3 
	Ribu = -grav*zu/ta*((dt-dter*jcool)+.61*ta*dq)/ut**2 
	nits = 3 
	if (Ribu .LT. 0) then 
	  zetu=CC*Ribu/(1+Ribu/Ribcu) 
	else 
	  zetu=CC*Ribu*(1+27/9*Ribu/CC)
	endif 

	zetaux = zetu
	if( abs(zetaux) < eps ) then
	  call get_act_timeline(aline)
	  call getinfo(iunit)
	  write(iunit,*) '====================== COARE error 1'
	  write(iunit,*) aline,zetaux,zu
	  write(iunit,*) sst,airt,ws,rain
	  write(iunit,*) Rns,Rnl,Qs,Q,rhoa
	  write(iunit,*) '======================'
	  flush(iunit)
	end if

	if( abs(zetu) < eps ) zetu = sign(eps,zetu)

	L10=zu/zetu 
	if (zetu .GT. 50) then 
	  nits=1 
	endif 
	usr=ut*von/(log(zu/zo10)-psi(1,zu/L10))
	tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psi(2,zt/L10)) 
	qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zot10)-
     +	    psi(2,zq/L10)) 
	tkt=.001
	charn=0.011 
	if (ut .GT. 10) then
	  charn=0.011+(ut-10)/(18-10)*(0.018-0.011) 
	endif 
	if (ut .GT. 18) then
	  charn=0.018 
	endif 

!	------------------------------------------------
!	Bulk loop
!	------------------------------------------------

	do i=1, nits 

          zet=von*grav*zu/ta*(tsr*(1+0.61*Q)+.61*ta*qsr)/(usr*usr)/
     &        (1+0.61*Q) 
	  !disp(usr)
	  !disp(zet) 
	  if (jwave .EQ. 0) zo=charn*usr*usr/grav+0.11*visa/usr  
	  if (jwave .EQ. 1) zo=50/2/pi*lwave*(usr/cwave)**4.5+0.11*
     &                         visa/usr 	!Oost et al
	  if (jwave .EQ. 2) zo=1200*hwave*(hwave/lwave)**4.5+0.11*
     &                         visa/usr 	!Taylor and Yelland
	  rr=zo*usr/visa 

	zetaux = zet
	if( abs(zetaux) < eps ) then
	  call get_act_timeline(aline)
	  call getinfo(iunit)
	  write(iunit,*) '====================== COARE error 2'
	  write(iunit,*) aline,zetaux,zu
	  write(iunit,*) sst,airt,ws,rain
	  write(iunit,*) Rns,Rnl,Qs,Q,rhoa
	  write(iunit,*) '======================'
	  flush(iunit)
	end if

	if( abs(zet) < eps ) zet = sign(eps,zet)

	  L=zu/zet
	  zoq=min(1.15e-4,5.5e-5/rr**.6) 
	  zot=zoq 
	  usr=ut*von/(log(zu/zo)-psi(1,zu/L)) 
	  tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot)-psi(2,zt/L)) 
	  qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zoq)-psi(2,zq/L)) 
	  Bf=-grav/ta*usr*(tsr+.61*ta*qsr) 
	  if (Bf .GT. 0) then
	    ug=Beta*(Bf*zi)**.333 
	  else
	    ug=.2 
	  endif
	  ut=sqrt(du*du+ug*ug) 
	  !Rnl=0.97*(5.67e-8*(sst-dter*jcool+kelv)**4-Rl) 
	  hsb=-rhoa*cpa*usr*tsr 
	  hlb=-rhoa*Le*usr*qsr 
	  qout=Rnl+hsb+hlb 
	  dels=Rns*(.065+11*tkt-6.6e-5/tkt*(1.-exp(-tkt/8.0e-4))) !Eq.16 Shortwave
	  qcol=qout-dels 
	  alq=Al*qcol+be*hlb*cpw/Le  				 !Eq. 7 Buoy flux water

	  if (alq .GT. 0) then 
	     xlamx=6./(1.+(bigc*alq/usr**4)**.75)**.333          !Eq 13 Saunders
	     tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr)                !Eq.11 Sub. thk
	  else
            xlamx=6.0 
            tkt=min(.01,xlamx*visw/(sqrt(rhoa/rhow)*usr))        !Eq.11 Sub. thk
	  endif 
     
	  dter=qcol*tkt/tcw 					 !Eq.12 Cool skin
	  dqer=wetc*dter 
   
	enddo !bulk iter loop

!	------------------------------------------------
!	Compute stress, sensible and latent
!	formulation in hlb already includes webb
!	------------------------------------------------

	tau=rhoa*usr*usr*du/ut                 !stress
	hsb=-rhoa*cpa*usr*tsr 		       !sensible
	hlb=-rhoa*Le*usr*qsr 		       !latent
     
!	------------------------------------------------
!	Rain heat flux
!	------------------------------------------------
     
	dwat=2.11e-5*((airt+kelv)/kelv)**1.94 		!water vapour diffusivity
	dtmp=(1.+3.309e-3*airt-1.44e-6*airt*airt)*	!heat diffusivity
     &        0.02411/(rhoa*cpa)
	alfac= 1./(1.+(wetc*Le*dwat)/(cpa*dtmp))    	!wet bulb factor
	RF= rain*alfac*cpw*((sst-airt-dter*jcool)+
     &	    (Qs-Q-dqer*jcool)*Le/cpa)

!	------------------------------------------------
!	compute transfer coeffs relative to ut @meas. ht
!	------------------------------------------------

	!write(177,*) 'Ce0: ',usr,qsr,dq-dqer*jcool,ut
	!call flush(177)

	Cd=tau/rhoa/ut/max(.1,du) 
	!Ch=-usr*tsr/ut/(dt-dter*jcool) 		!possible divide by 0
	!Ce=-usr*qsr/(dq-dqer*jcool)/ut 		!possible divide by 0

!       ------------------------------------------------
!       prepare output variables
!       ------------------------------------------------

        qsens = -hsb - RF
	qlat  = -hlb
        evap  = hlb / Le

!       ------------------------------------------------
!       end of routine
!       ------------------------------------------------

        end subroutine coare30

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Compute momentum fluxes due to rainfall

	subroutine wcstress(ws,rhoa,rhow,rain,
     +			    wx,wy,taux,tauy)

	implicit none

! input
	real, intent(in)	:: ws		!wind speed [m/s]
	real, intent(in)	:: rhoa		!moist air density (kg/m3)
	real, intent(in)	:: rhow		!water density (kg/m3)
	real, intent(in)	:: rain		!precipitation (kg/m2/s)
        real, intent(in)	:: wx,wy      	!wind components [m/s]
! output
	real, intent(out)	:: taux,tauy	!wind stress components
! local
	real 			:: cff

!       ------------------------------------------------
!       Compute momentum flux (N/m2) due to rainfall (kg/m2/s)
!       ------------------------------------------------

        cff  = 0.85 * rain * ws
        taux  = taux + cff*sign(1.0,wx)
        tauy  = tauy + cff*sign(1.0,wy)

	end subroutine wcstress

!----------------------------------------------------------------------
!-----------------------------------------------------------------------
! This subroutine compute sea surface skin temperature according to
! Zeng and Beljaars (2005) to be used in the bulk flux computations

	subroutine tw_skin(rad,qrad,tw,hb,usw,dt,dtw,tws)

	implicit none

	include 'subqfxm.h'

! input
	real, intent(in)	:: rad	!Net shortwave flux [W/m2] (+ down)
	real, intent(in)	:: qrad	!LH + SH + LW [W/m2] (+ down)
	real, intent(in)	:: tw 	!sea temperature at depth hb [C]
	real, intent(in)	:: hb	!depth of modelled tw [m]
	real, intent(in)	:: usw	!surface friction velocity [m/s]
	real, intent(in)	:: dt 	!time step [s]
! input-output
	real, intent(inout)	:: dtw	!Warm layer temp. diff. [C]
! output
	real, intent(out)	:: tws	!skin temperature [C]
! parameters
	real, parameter		:: num  = 0.3    !exponent of warm layer T profile
! local
	real			:: sw	!Net shortwave flux (degC m/s), + down
	real			:: q	!LH + SH + LW [degC m/s] (+ down)
	real			:: alw	!thermal expansion coefficient
	real			:: nu	!num computed following Lee et al, 2013
	real 			:: zeta	!db / L  (L is Monin-Obukhov length)
	real 			:: qn	!Q + R_s - R(-d) heat flux in warm layer
	real			:: ds	!cool skin layer depth (m; positive)
	real			:: fs	!frac. of solar rad. absorbed in the cool sublayer
	real			:: f1	!frac. of solar rad. absorbed in the warm layer
	real			:: phi	!stability function
	real			:: dtc	!Cool skin temp. diff. (deg C)
	real			:: dtwo	!dtw frmo previous time step
	real			:: q2,qn1
	real			:: cff,cff1,cff2,cff3,cff4,cff5

	!write(6,*) rad,qrad,tw,hb,usw,dt,dtw,tws

!       ---------------------------------------------------------
!       Set parameters
!       ---------------------------------------------------------

	cff  = 1.0 / (rhow * cpw)	!convert [W/m2] to [C m/s]
        sw   = rad * cff
	q    = qrad * cff 
	dtwo = dtw			!old dtw
        alw  = 1.e-5*max(tw,1.)

!	------------------------------------------------
!	Cool skin/sublayer (eq. 4 of Z&B-2005)
!	------------------------------------------------

        cff4 = 16.*grav*alw*visw**3/diff**2
        cff5 = cff4/usw**4
        q2 = max(cff,-q)	                ! avoid iterative procedure 
                                                ! between zs and fs (impact <0.03C)
        ds = 6./(1.+(cff5*q2)**0.75)**0.333     ! sublayer thickness (m)
        ds = ds*visw/usw                        ! --> eq 6 of Z&B-2005
        fs = 0.065+11.*ds-(6.6e-5/ds)           ! fract. of solar rad. absorbed 
     &            *(1.-exp(-ds/8.e-4))          ! in sublayer
        fs = max(fs,0.01)
        dtc = ds*(q+sw*fs)/diff                 ! cool skin temp. diff
        dtc = min(dtc,0.)                       ! --> eq. 4 of Z&B-2005

!	------------------------------------------------
!	Warm layer  (eq. 11 of Z&B-2005)
!	------------------------------------------------

        f1  = 1.-0.27*exp(-2.8*hb) - 0.45*exp(-0.07*hb)
        qn  = q + sw*f1

        if(qn .lt. 0.) then                          
          cff1 = sqrt(5.*hb*grav*alw/num)   ! allow warm layer to subsist after
          qn1  = sqrt(dtwo)*usw**2/cff1     ! sunset by changing heat flux in
          qn   = max(qn,qn1)                ! Monin-Obukhov length (eq. 12 Z&B-2005)
        endif

        cff2 = von*grav*alw
        zeta = hb*cff2*qn/usw**3            ! hb/LMO

!	------------------------------------------------
!	The similarity function as given in Takaya et al, 2010
!	------------------------------------------------

        if(zeta .ge. 0.) then
          phi = 1.+((5.*zeta+4.*zeta**2)/(1.+3.*zeta+0.25*zeta**2))
          !phi=1.+5.*zeta		!Z&B-2005
        else
          phi=1./sqrt(1.-16.*zeta)
        endif
        cff3 = von*usw/(hb*phi)

!	------------------------------------------------
!	Compute nu following Lee et al, 2013
!	still not tested, use nu = num
!	------------------------------------------------

 	!nu = exp(-0.02*zeta)
	!nu = max(num,nu)
	!nu = min(nu,1.)
	nu = num

!	------------------------------------------------
!	Implicit time stepping (eq. 11 of Z&B-2005)
!	dtw(n+1)-dtw(n) = dt*RHS(dtw(n+1))
!	------------------------------------------------

        dtw = (dtwo + (nu+1.)/nu*(q+sw*f1)*dt/hb)
     &                      /(1.+(nu+1.)*cff3*dt)
        dtw = max(0.,dtw)

!	------------------------------------------------
!	Get skin temperature
!	------------------------------------------------

        tws = tw + dtw + dtc

	end subroutine tw_skin

!----------------------------------------------------------------------
!-----------------------------------------------------------------------
