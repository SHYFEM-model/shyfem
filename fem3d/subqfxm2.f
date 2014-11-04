c
c $Id: subqfxm2.f,v 1.2 2009-09-14 08:20:58 georg Exp $
c
c heat flux module (Gill)
c
c contents :
c
c subroutine heatgill(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)
c	computes heat fluxes from formulas in Gill
c subroutine tempgill(dt,dh,qsol,t,p,w,ur,cc,ts,tsnew,rtot,evap)
c	computes heat fluxes from formulas in Gill and adjusts water temp
c
c revision log :
c
c 27.08.2009    ggu     call to heatgill changed (pass ur, and not e,q)
c
c*******************************************************************

	subroutine heatgill(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)

c computes heat fluxes from formulas in Gill
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
c	real q		!specific humidity [0-1]
c	real e		!vapor pressure [mb]
c	real r		!mixing ratio [0-1]

	real rdry,lv0,pascal
	parameter(rdry=287.04,lv0=2.5008e6,pascal=100.)
	real sigma,epsbbb,sigma0
	parameter(sigma=5.67e-8,epsbbb=0.985,sigma0=sigma*epsbbb)

	real rho,cp,lv,tv
	real e,r,q
	real es,rs,qs
	real dt,ch,ce
	real theta,theta2,theta4
	real cloud

c	------------------------------------------------
c	initialization
c	------------------------------------------------

	call vapor(t,p,ur,e,r,q)	!compute e,r,q

	tv = ( t + kelv ) * ( 1. + 0.6078 * q )
	rho = pascal * p / ( rdry * tv )
	lv = lv0 - 2.3e3 * ts
	cp = 1004.6 * ( 1. + 0.8375 * q )

	call satur(ts,p,es,rs,qs)	!compute saturation values for sea

c	------------------------------------------------
c	sensible heat flux
c	------------------------------------------------

	dt = ts - t
	if( dt .lt. 0 ) then	!stable
	  ch = 0.83e-3
	else			!unstable
	  ch = 1.10e-3
	end if

	qsens = rho * cp * ch * w * dt

c	------------------------------------------------
c	latent heat flux
c	------------------------------------------------

	ce = 1.5e-3

	evap  = rho * ce * w * ( qs - q )
	qlat = lv * evap

c	------------------------------------------------
c	long wave radiation
c	------------------------------------------------

	theta = ts + 273.15
	theta2 = theta * theta
	theta4 = theta2 * theta2
	cloud = 1. - 0.6 * cc * cc

	qlong = sigma0 * theta4 * ( 0.39 - 0.05 * sqrt(e) ) * cloud

c	------------------------------------------------
c	end of routine
c	------------------------------------------------

        if( qlong .le. 0. ) then
          write(6,*) qsens,qlat,qlong
          write(6,*) ts,theta,theta4
          write(6,*) e,sqrt(e),0.05*sqrt(e),cc,cloud
          write(6,*) t,w,q,qs
          write(6,*) tv,lv,cp,rho
          stop 'error stop heatgill: qlong'
        end if

	end

c*******************************************************************

	subroutine tempgill(dt,dh,qsol,t,p,w,ur,cc,ts,tsnew,rtot,evap)

c computes heat fluxes from formulas in Gill and adjusts water temperature
c
c heat fluxes are positive upward (from sea to atmosphere)

	implicit none

	include 'subqfxm.h'

	real dt		!time step [s]				- in
	real dh		!depth of layer [m]			- in
	real qsol	!solar radiation [W/m**2]		- in
	real t		!temperature [C]			- in
	real p		!pressure [mb]				- in
	real w		!wind speed [m/s]			- in
	real ur		!relative humidity [%] ([0-100])	- in
	real cc		!cloud cover [0-1]			- in
	real ts		!sea temperature [C]			- in
	real tsnew	!new sea temperature [C]		- out
	real rtot	!total heat input [W/m**2]		- out
	real evap	!evaporation [kg/m**2/s]		- out

        real ct
        real qsens,qlat,qlong

c	------------------------------------------------
c	constants
c	------------------------------------------------

	ct   = cpw * rhow * dh	!heat capacity / area

c	------------------------------------------------
c	compute total radiation - positive if into water
c	------------------------------------------------

	call heatgill(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)
        write(67,*) qsol,-qsens,-qlat,-qlong

	rtot = qsol - ( qsens + qlat + qlong )

c	------------------------------------------------
c	compute new temperature: dQ = dT * rho * cpw * dh / dt
c	------------------------------------------------

	tsnew = ts + rtot * dt / ct

c	------------------------------------------------
c	end of routine
c	------------------------------------------------

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

