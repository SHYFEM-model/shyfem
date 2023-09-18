! $Id: subqfxm2.f,v 1.2 2009-09-14 08:20:58 georg Exp $
!
! heat flux module (GOTM)
!
! contents :
!
! subroutine heatgotm(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap,cdd)
!       computes heat fluxes from GOTM heat module
! subroutine heatgotm_exchange_coeff(w,cdd,chd,ced)
!	computes exchange coefficients
! subroutine tempgotm(dt,dh,qsol,t,p,w,ur,cc,ts,tsnew,rtot,evap)
!       computes heat fluxes from formulas in GOTM and adjusts water temp
!
! revision log :
!
! 04.03.2011    ggu     heat module from gotm extracted
! 23.03.2011    ggu     bug in heatgotm() -> avoid wind speed == 0
! 29.03.2011    ggu     bug in heatgotm() -> convert pressure to Pascal
! 16.06.2014    ccf     bug in heatgotm() -> airt instead of sst in qb
!
!***********************************************************************
!-----------------------------------------------------------------------
        module heat_gotm
!-----------------------------------------------------------------------
        contains
!-----------------------------------------------------------------------

	subroutine heatgotm(airt,airp,ws,rh,cloud,sst,qh,qe,qb,evap)

	implicit none

	include 'subqfxm.h'

        double precision airt       !air temperature [C]                    - in
        double precision airp       !pressure [mb]                          - in
        double precision ws         !wind speed [m/s]                       - in
        double precision rh         !relative humidity [%] ([0-100])        - in
        double precision cloud      !cloud cover [0-1]                      - in
        double precision sst        !sea temperature [C]                    - in
        double precision qh         !sensible heat flux [W/m**2]            - out
        double precision qe         !latent heat flux [W/m**2]              - out
        double precision qb         !long wave radiation [W/m**2]           - out
        double precision evap       !evaporation [kg/m**2/s]                - out

	logical bdebug
	integer ihback
	double precision tmp
	double precision L,es,qs
	double precision ea,qa
	double precision tvirt,rho_air,s0,s
	double precision cdd,chd,ced
	double precision x,heat
	double precision cddrag     !wind drag coefficient [-]
	double precision w,airppa

	ihback = 1

	!bdebug = ws .eq. 0.
	bdebug = .false.

	w = max(ws,0.01)
	airppa = 100. * airp	!pressure in Pascal

!------------------------------------------------------------
! compute exchange coefficients
!------------------------------------------------------------

	!w = sqrt(wx*wx+wy*wy)
	L = (2.5-0.00234*sst)*1.e6
	es = a1 +sst*(a2+sst*(a3+sst*(a4+sst*(a5+sst*(a6+sst*a7)))))
	es = es * 100.0 ! Conversion millibar --> Pascal
	qs = const06*es/(airppa-0.377*es) ! specific humidity at sea surface

        ea = a1 +airt*(a2+airt*(a3+airt*(a4+airt*(a5+airt*(a6+airt*a7)))))
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

!------------------------------------------------------------
! compute transfer coefficients
!------------------------------------------------------------

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

!------------------------------------------------------------
! compute heat fluxes
!------------------------------------------------------------

	qe=ced*L*rho_air*w*(qs-qa)            ! latent
	qh=chd*cpa*rho_air*w*(sst-airt)       ! sensible

	tmp=airt+kelv
	if( ihback .eq. 1 ) then		!clark
	  qb=(1.0-.8*cloud*cloud)*emiss*bolz*(tmp**4)*(0.39-0.05*sqrt(ea/100.0))        &
     &            +4.0*emiss*bolz*(tmp**3)*(sst-airt)
	else if(ihback .eq. 2 ) then		!hastenrath
          qb=(1.0-.8*cloud*cloud)*emiss*bolz*(tmp**4)*(0.39-0.056*sqrt(1000*qa))        &
     &            +4.0*emiss*bolz*(tmp**3)*(sst-airt)
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

!***********************************************************************

	subroutine heatgotm_exchange_coeff(w,cdd,chd,ced)

	implicit none

	include 'subqfxm.h'

	double precision w,cdd,chd,ced

	double precision ae_d,be_d,pe_d
	double precision ae_h,be_h,pe_h
	double precision ae_e,be_e,pe_e
	double precision ce_h,ce_e

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

!***********************************************************************

        subroutine tempgotm(dt,dh,qsol,t,p,w,ur,cc,ts,tsnew,rtot,evap)

! computes heat fluxes from formulas in Gill and adjusts water temperature
!
! heat fluxes are positive upward (from sea to atmosphere)

        implicit none

        double precision dt         !time step [s]                          - in
        double precision dh         !depth of layer [m]                     - in
        double precision qsol       !solar radiation [W/m**2]               - in
        double precision t          !temperature [C]                        - in
        double precision p          !pressure [mb]                          - in
        double precision w          !wind speed [m/s]                       - in
        double precision ur         !relative humidity [%] ([0-100])        - in
        double precision cc         !cloud cover [0-1]                      - in
        double precision ts         !sea temperature [C]                    - in
        double precision tsnew      !new sea temperature [C]                - out
        double precision rtot       !total heat input [W/m**2]              - out
        double precision evap       !evaporation [kg/m**2/s]                - out

        double precision cw,rhow,ct
        double precision qsens,qlat,qlong
	double precision cddrag	!drag coefficient computed by heat module

!       ------------------------------------------------
!       constants
!       ------------------------------------------------

        cw   = 3991.            !heat capacity of water
        rhow = 1026.            !density of water
        ct   = cw * rhow * dh   !heat capacity / area

!       ------------------------------------------------
!       compute total radiation - positive if into water
!       ------------------------------------------------

        call heatgotm(t,p,w,ur,cc,ts,qsens,qlat,qlong,evap)
        write(67,*) qsol,-qsens,-qlat,-qlong

        rtot = qsol - ( qsens + qlat + qlong )

!       ------------------------------------------------
!       compute new temperature: dQ = dT * rho * cw * dh / dt
!       ------------------------------------------------

        tsnew = ts + rtot * dt / ct

!       ------------------------------------------------
!       end of routine
!       ------------------------------------------------

        end

!***********************************************************************

!-----------------------------------------------------------------------
        end module heat_gotm
!-----------------------------------------------------------------------
