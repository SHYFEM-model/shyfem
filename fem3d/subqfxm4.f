
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

c heat flux module (GOTM)
c
c contents :
c
c subroutine heatgotm(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap,cdd)
c       computes heat fluxes from GOTM heat module
c subroutine heatgotm_exchange_coeff(w,cdd,chd,ced)
c	computes exchange coefficients
c subroutine tempgotm(dt,dh,qsol,t,p,w,ur,cc,ts,tsnew,rtot,evap)
c       computes heat fluxes from formulas in GOTM and adjusts water temp
c
c revision log :
c
c 04.03.2011	ggu	heat module from gotm extracted
c 23.03.2011	ggu	bug in heatgotm() -> avoid wind speed == 0
c 29.03.2011	ggu	bug in heatgotm() -> convert pressure to Pascal
c 14.04.2011	ggu	changed VERS_6_1_22
c 16.06.2014	ccf	bug in heatgotm() -> airt instead of sst in qb
c 05.11.2014	ggu	changed VERS_7_0_5
c 24.01.2018	ggu	changed VERS_7_5_41
c 16.02.2019	ggu	changed VERS_7_5_60
c
c***********************************************************************

	subroutine heatgotm(airt,airp,ws,rh,cloud,sst
     +					,qh,qe,qb,evap)

	implicit none

	include 'subqfxm.h'

        real airt       !air temperature [C]                    - in
        real airp       !pressure [mb]                          - in
        real ws         !wind speed [m/s]                       - in
        real rh         !relative humidity [%] ([0-100])        - in
        real cloud      !cloud cover [0-1]                      - in
        real sst        !sea temperature [C]                    - in
        real qh         !sensible heat flux [W/m**2]            - out
        real qe         !latent heat flux [W/m**2]              - out
        real qb         !long wave radiation [W/m**2]           - out
        real evap       !evaporation [kg/m**2/s]                - out

	logical bdebug
	integer ihback
	real tmp
	real L,es,qs
	real ea,qa
	real tvirt,rho_air,s0,s
	real cdd,chd,ced
	real x,heat
	real cddrag     !wind drag coefficient [-]
	real w,airppa

        if( airp > 10000 ) stop 'error stop heatgotm: p not in mbar'

	ihback = 1

	!bdebug = ws .eq. 0.
	bdebug = .false.

	w = max(ws,0.01)
	airppa = 100. * airp	!pressure in Pascal

c------------------------------------------------------------
c compute exchange coefficients
c------------------------------------------------------------

	!w = sqrt(wx*wx+wy*wy)
	L = (2.5-0.00234*sst)*1.e6
	es = a1 +sst*(a2+sst*(a3+sst*(a4+sst*(a5+sst*(a6+sst*a7)))))
	es = es * 100.0 ! Conversion millibar --> Pascal
	qs = const06*es/(airppa-0.377*es) ! specific humidity at sea surface

        ea = a1 +airt*(a2+airt*(a3+airt*(a4+airt*
     +					(a5+airt*(a6+airt*a7)))))
        ea = rh*ea ! millibar --> Pascal and saturation --> actual vap. press.
        qa = const06*ea/(airppa-0.377*ea) ! specific humidity at 2m

	tvirt = (airt+kelv)*(1+qa/const06)/(1+qa)
	rho_air = airppa/(287.05*Tvirt)

!	Stability
	s0=0.25*(sst-airt)/(w+1.0e-10)**2
	s=s0*abs(s0)/(abs(s0)+0.01)

	if( bdebug ) then
	  write(6,*) 'heatgotm 1: ',es,qs,ea,qa,tvirt,s0,s
	end if

c------------------------------------------------------------
c compute transfer coefficients
c------------------------------------------------------------

	call heatgotm_exchange_coeff(w,cdd,chd,ced)

	if(s .lt. 0.) then
          if (s .gt. -3.3) then
            x = 0.1+0.03*s+0.9*exp(4.8*s)
          else
            x = 0.0
          end if
          cdd=x*cdd
          chd=x*chd
          ced=x*ced
	else
          cdd=cdd*(1.0+0.47*sqrt(s))
          chd=chd*(1.0+0.63*sqrt(s))
          ced=ced*(1.0+0.63*sqrt(s))
	end if

	if( bdebug ) then
	  write(6,*) 'heatgotm 2: ',cdd,chd,ced
	end if

c------------------------------------------------------------
c compute heat fluxes
c------------------------------------------------------------

	qe=ced*L*rho_air*w*(qs-qa)            ! latent
	qh=chd*cpa*rho_air*w*(sst-airt)       ! sensible

	tmp=airt+kelv
	if( ihback .eq. 1 ) then		!clark
	  qb=(1.0-.8*cloud*cloud)
     +            *emiss*bolz*(tmp**4)*(0.39-0.05*sqrt(ea/100.0))
     +            +4.0*emiss*bolz*(tmp**3)*(sst-airt)
	else if(ihback .eq. 2 ) then		!hastenrath
          qb=(1.0-.8*cloud*cloud)
     +            *emiss*bolz*(tmp**4)*(0.39-0.056*sqrt(1000*qa))
     +            +4.0*emiss*bolz*(tmp**3)*(sst-airt)
	else
	  stop 'error stop: ihback not allowed'
	end if

	evap = qe / L

	heat = -(qe+qh+qb)
	cddrag = cdd

	!tmp = -cdd*rho_air*w	!negative????
	!taux  = tmp*wx
	!tauy  = tmp*wy

	end

c***********************************************************************

	subroutine heatgotm_exchange_coeff(w,cdd,chd,ced)

	implicit none

	include 'subqfxm.h'

	real w,cdd,chd,ced

	real ae_d,be_d,pe_d
	real ae_h,be_h,pe_h
	real ae_e,be_e,pe_e
	real ce_h,ce_e

	if (w .lt. 2.2) then

          ae_d=0.0
	  be_d=1.08
	  pe_d=-0.15
	
          ae_h=0.0
	  be_h=1.185
	  ce_h=0.0
	  pe_h=-0.157
	
          ae_e=0.0
	  be_e=1.23
	  ce_e=0.0
	  pe_e=-0.16
	
	else if (w .lt. 5.0) then

          ae_d=0.771
	  be_d=0.0858
	  pe_d=1.0
	
          ae_h=0.927
	  be_h=0.0546
	  ce_h=0.0
	  pe_h=1.0
	
          ae_e=0.969
	  be_e=0.0521
	  ce_e=0.0
	  pe_e=1.0
	
	else if (w .lt. 8.0) then

          ae_d=0.867
	  be_d=0.0667
	  pe_d=1.0
	
          ae_h=1.15
	  be_h=0.01
	  ce_h=0.0
	  pe_h=1.0
	
          ae_e=1.18
	  be_e=0.01
	  ce_e=0.0
	  pe_e=1.0
	
	else if (w .lt. 25.0) then

          ae_d=1.2
	  be_d=0.025
	  pe_d=1.0
	
          ae_h=1.17
	  be_h=0.0075
	  ce_h=-0.00045
	  pe_h=1.0
	
          ae_e=1.196
	  be_e=0.008
	  ce_e=-0.0004
	  pe_e=1.0

	else

          ae_d=0.0
	  be_d=0.073
	  pe_d=1.0
	
          ae_h=1.652
	  be_h=-0.017
	  ce_h=0.0
	  pe_h=1.0
	
          ae_e=1.68
	  be_e=-0.016
	  ce_e=0
	  pe_e=1.0
	
	end if

	cdd=(ae_d+be_d*exp(pe_d*log(w+eps)))*1.0e-3
	chd=(ae_h+be_h*exp(pe_h*log(w+eps))+ce_h*(w-8.0)**2)*1.0e-3
	ced=(ae_e+be_e*exp(pe_e*log(w+eps))+ce_e*(w-8.0)**2)*1.0e-3

	end

c***********************************************************************

        subroutine tempgotm(dt,dh,qsol,t,p,w,ur,cc,ts,tsnew,rtot,evap)

c computes heat fluxes from formulas in Gill and adjusts water temperature
c
c heat fluxes are positive upward (from sea to atmosphere)

        implicit none

        real dt         !time step [s]                          - in
        real dh         !depth of layer [m]                     - in
        real qsol       !solar radiation [W/m**2]               - in
        real t          !temperature [C]                        - in
        real p          !pressure [mb]                          - in
        real w          !wind speed [m/s]                       - in
        real ur         !relative humidity [%] ([0-100])        - in
        real cc         !cloud cover [0-1]                      - in
        real ts         !sea temperature [C]                    - in
        real tsnew      !new sea temperature [C]                - out
        real rtot       !total heat input [W/m**2]              - out
        real evap       !evaporation [kg/m**2/s]                - out

        real cw,rhow,ct
        real qsens,qlat,qlong
	real cddrag	!drag coefficient computed by heat module

c       ------------------------------------------------
c       constants
c       ------------------------------------------------

        cw   = 3991.            !heat capacity of water
        rhow = 1026.            !density of water
        ct   = cw * rhow * dh   !heat capacity / area

c       ------------------------------------------------
c       compute total radiation - positive if into water
c       ------------------------------------------------

        call heatgotm(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)
        write(67,*) qsol,-qsens,-qlat,-qlong

        rtot = qsol - ( qsens + qlat + qlong )

c       ------------------------------------------------
c       compute new temperature: dQ = dT * rho * cw * dh / dt
c       ------------------------------------------------

        tsnew = ts + rtot * dt / ct

c       ------------------------------------------------
c       end of routine
c       ------------------------------------------------

        end

c***********************************************************************

