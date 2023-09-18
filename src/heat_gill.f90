!
! $Id: subqfxm2.f,v 1.2 2009-09-14 08:20:58 georg Exp $
!
! heat flux module (Gill)
!
! contents :
!
! subroutine heatgill(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)
!	computes heat fluxes from formulas in Gill
! subroutine tempgill(dt,dh,qsol,t,p,w,ur,cc,ts,tsnew,rtot,evap)
!	computes heat fluxes from formulas in Gill and adjusts water temp
!
! revision log :
!
! 27.08.2009    ggu     call to heatgill changed (pass ur, and not e,q)
!
!*******************************************************************
!-------------------------------------------------------------------
        module heat_gill
!-------------------------------------------------------------------
        contains
!-------------------------------------------------------------------

	subroutine heatgill(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)

! computes heat fluxes from formulas in Gill
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
!	double precision q		!specific humidity [0-1]
!	double precision e		!vapor pressure [mb]
!	double precision r		!mixing ratio [0-1]

	double precision rdry,lv0,pascal
	parameter(rdry=287.04,lv0=2.5008e6,pascal=100.)
	double precision sigma,epsbbb,sigma0
	parameter(sigma=5.67e-8,epsbbb=0.985,sigma0=sigma*epsbbb)

	double precision rho,cp,lv,tv
	double precision e,r,q
	double precision es,rs,qs
	double precision dt,ch,ce
	double precision theta,theta2,theta4
	double precision cloud

!	------------------------------------------------
!	initialization
!	------------------------------------------------

	call vapor(t,p,ur,e,r,q)	!compute e,r,q

	tv = ( t + kelv ) * ( 1. + 0.6078 * q )
	rho = pascal * p / ( rdry * tv )
	lv = lv0 - 2.3e3 * ts
	cp = 1004.6 * ( 1. + 0.8375 * q )

	call satur(ts,p,es,rs,qs)	!compute saturation values for sea

!	------------------------------------------------
!	sensible heat flux
!	------------------------------------------------

	dt = ts - t
	if( dt .lt. 0 ) then	!stable
	  ch = 0.83e-3
	else			!unstable
	  ch = 1.10e-3
	end if

	qsens = rho * cp * ch * w * dt

!	------------------------------------------------
!	latent heat flux
!	------------------------------------------------

	ce = 1.5e-3

	evap  = rho * ce * w * ( qs - q )
	qlat = lv * evap

!	------------------------------------------------
!	long wave radiation
!	------------------------------------------------

	theta = ts + 273.15
	theta2 = theta * theta
	theta4 = theta2 * theta2
	cloud = 1. - 0.6 * cc * cc

	qlong = sigma0 * theta4 * ( 0.39 - 0.05 * sqrt(e) ) * cloud

!	------------------------------------------------
!	end of routine
!	------------------------------------------------

        if( qlong .le. 0. ) then
          write(6,*) qsens,qlat,qlong
          write(6,*) ts,theta,theta4
          write(6,*) e,sqrt(e),0.05*sqrt(e),cc,cloud
          write(6,*) t,w,q,qs
          write(6,*) tv,lv,cp,rho
          stop 'error stop heatgill: qlong'
        end if

	end

!*******************************************************************

	subroutine tempgill(dt,dh,qsol,t,p,w,ur,cc,ts,tsnew,rtot,evap)

! computes heat fluxes from formulas in Gill and adjusts water temperature
!
! heat fluxes are positive upward (from sea to atmosphere)

	implicit none

	include 'subqfxm.h'

	double precision dt		!time step [s]				- in
	double precision dh		!depth of layer [m]			- in
	double precision qsol	!solar radiation [W/m**2]		- in
	double precision t		!temperature [C]			- in
	double precision p		!pressure [mb]				- in
	double precision w		!wind speed [m/s]			- in
	double precision ur		!relative humidity [%] ([0-100])	- in
	double precision cc		!cloud cover [0-1]			- in
	double precision ts		!sea temperature [C]			- in
	double precision tsnew	!new sea temperature [C]		- out
	double precision rtot	!total heat input [W/m**2]		- out
	double precision evap	!evaporation [kg/m**2/s]		- out

        double precision ct
        double precision qsens,qlat,qlong

!	------------------------------------------------
!	constants
!	------------------------------------------------

	ct   = cpw * rhow * dh	!heat capacity / area

!	------------------------------------------------
!	compute total radiation - positive if into water
!	------------------------------------------------

	call heatgill(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)
        write(67,*) qsol,-qsens,-qlat,-qlong

	rtot = qsol - ( qsens + qlat + qlong )

!	------------------------------------------------
!	compute new temperature: dQ = dT * rho * cpw * dh / dt
!	------------------------------------------------

	tsnew = ts + rtot * dt / ct

!	------------------------------------------------
!	end of routine
!	------------------------------------------------

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

!-------------------------------------------------------------------
        end module heat_gill
!-------------------------------------------------------------------
