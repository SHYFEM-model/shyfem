!
! $Id: subqfxm1.f,v 1.4 2008-10-10 09:29:54 georg Exp $
!
! heat flux module default
!
! contents :
!
! subroutine subtem(dt,dh,qs,uw,cc,ta,tb,ts,tsnew,rtot)
!			compute water temperature in element of depth dh
! subroutine subtem0(dt,dh,qs,uw,cc,ta,tb,ts,tsnew,evap)
!			same as subtem, but use heatlucia
! subroutine heatlucia(t,p,w,tb,cc,ts,qsens,qlat,qlong,evap)
!			computes heat fluxes from Lucia modules
! subroutine longwave(ts,ta,tb,cc,rb)
!			long-wave radiation
! subroutine evcon(ts,ta,tb,uw,p,re,rc)
!			latent heat and convective heat
!
! revision log :
!
! 01.06.1998	ggu&lz	written from scratch (nearly)
! 24.06.1998	ggu&lz	subroutines from lucia integrated
! 30.04.2001	ggu	new routine qtotal_tb
! 09.12.2002	ggu	cleaned and re-arranged
! 23.03.2006	ggu	changed time step to double precision
! 16.02.2011	ggu	pstd introduced
!
! notes :
!
! qs	net solar radiation (reflection already subtracted)
! ta	air temperature
! tb	wet bulb temperature
! uw	wind speed
! cc	cloud cover (0 clear sky, 1 totally covered)
!
!*****************************************************************************
!-----------------------------------------------------------------------------
      module heat_default
!-----------------------------------------------------------------------------
      contains
!-----------------------------------------------------------------------------

      subroutine subtem(dt,dh,qs,uw,cc,ta,tb,ts,tsnew,evap)

!  Calcola la temperatura dell'acqua in un bacino di profondita' dh,
!  risolvendo l'equazione del bilancio termico all'interfaccia acqua-aria.
!
!  dt = timestep del modello FEM
!  dh = profondita' dello strato
!  qs = radiazione solare (short wave radiation) entrante nella superficie 
!       dell'acqua (radiazione incidente al netto della riflessione)
!  uw = velocita' del vento
!  cc = copertura nuvolosa (tra 0 e 1)
!  ts,ta,tb = temperatura dell'acqua, dell'aria, temp. dell'aria a bulbo
!             bagnato
!  tsnew = nuova temperatura dell'acqua
!
!  pa = pressione atmosferica in mbar in subroutine evcon 
!	(!!!!! verificare la sensibilita' a pa)
!
!  cpw = calore specifico dell'acqua di mare (J/kg C) da Gill, per ts=16 C
!  rhow = densita' dell'acqua di mare (kg/m3) da Gill, per ts=16 C
!
!  rtot = total radiation
!  qs = solar radiation (short wave)
!  rb = long wave radiation
!  re = evaporation term
!  rc = convection term 
!
!  Formulazione del bilancio termico da:
!  Mariotti M. e C. Dejak, 1982: 'Valutazione dei flussi energetici
!  all'interfaccia acqua-aria nella Laguna di Venezia', Ambiente Risorse 10.
!
!  Vengono usate le unita' SI.
!
!  L. Zampato - Dicembre 1997

      implicit none

      include 'subqfxm.h'

      double precision dt
      double precision dh
      double precision p
      double precision qs, uw, cc
      double precision ta, tb, ts, tsnew
      double precision evap

      double precision pstd
      parameter ( pstd = 1013.25 )

      double precision rb, re, rc, rtot
      double precision ct

! constants

      p = pstd

!  long-wave term
!   input: ts,ta,tb,cc
!   output: rb = long wave radiation

      call longwave(ts,ta,tb,cc,rb)

!  evaporation and convection terms
!   input: ts,ta,tb,uw
!   output: re,rc = evaporation and convection terms

      call evcon(ts,ta,tb,uw,p,re,rc)

!  total radiation: rtot (W/m2) (positive if into water)

      rtot = qs+rb+re+rc

!  evaporation [kg/(m**2 s)]

      evap = re / 2.5e+6                !divide by latent heat of evaporation
      evap = evap / rhow                 !in [m/s]
      evap = evap * 1000. * 86400.      !in [mm/day]

!  heat capacity/area

      ct = cpw*rhow*dh

!  new temperature      formula:  dQ = dT * rho * cpw * dh / dt

      tsnew = ts + rtot*dt/ct

      end

!*****************************************************************************

      subroutine subtem0(dt,dh,qs,uw,cc,ta,tb,ts,tsnew,evap)

! same as subtem, but use heatlucia

      implicit none

      include 'subqfxm.h'

      double precision dt
      double precision dh
      double precision qs, uw, cc
      double precision ta, tb, ts, tsnew
      double precision evap

      double precision pstd
      parameter ( pstd = 1013.25 )

      double precision rtot
      double precision qsens,qlat,qlong
      double precision ct, p

! constants

      ct = cpw*rhow*dh	!heat capacity	[ J / (m**2 K) ]
      p = pstd

      call heatlucia(ta,p,uw,tb,cc,ts,qsens,qlat,qlong,evap)

      rtot = qs - ( qlong + qlat + qsens )

!  new temperature      formula:  dQ = dT * rho * cpw * dh / dt

      tsnew = ts + rtot*dt/ct

      evap = -evap			!in [kg/(m**2 s)]
      evap = evap / rhow                 !in [m/s]
      evap = evap * 1000. * 86400.      !in [mm/day]

      end

!*****************************************************************************

        subroutine heatlucia(t,p,w,tb,cc,ts,qsens,qlat,qlong,evap)

! computes heat fluxes from Lucia modules
!
! heat fluxes are positive upward (from sea to atmosphere)

        implicit none

        double precision t          !air temperature [C]                    - in
        double precision p          !pressure [mb]                          - in
        double precision w          !wind speed [m/s]                       - in
        double precision tb         !wet bulb temperature [C]               - in
        double precision cc         !cloud cover [0-1]                      - in
        double precision ts         !sea temperature [C]                    - in
        double precision qsens      !sensible heat flux [W/m**2]            - out
        double precision qlat       !latent heat flux [W/m**2]              - out
        double precision qlong      !long wave radiation [W/m**2]           - out
        double precision evap       !evaporation [kg/(m**2 s)]              - out

	double precision ta,uw
	double precision rb,re,rc

	ta = t		!t air
	uw = w		!wind speed

	call longwave(ts,ta,tb,cc,rb)
	call evcon(ts,ta,tb,uw,p,re,rc)

	qlong = -rb
	qsens = -rc
	qlat  = -re
	evap  = -re / 2.5e+6

	end

!*****************************************************************************

      subroutine longwave(ts,ta,tb,cc,rb)

!  termine di radiazione long-wave nel bilancio termico atmosfera-mare
!
!  es, ea = emissivita' del mare e dell'aria
!  v = tensione di vapore alla temperatura t
!  ra, rb = radiazione emessa dall'atmosfera e dal mare
!  rb = termine long-wave in W/m2 da inserire nel bilancio termico
!  rbkj = termine long-wave in kJ/(m2 ora)
!
      implicit none
!
      double precision ts, ta, tb
      double precision cc
      double precision sigma, es, ea
      double precision alpha, beta
      double precision va, vd, psi, rr
      double precision ra, rs, rbkj, rb
!
      sigma=5.67/(10**(8))
      es=0.97
      ea=0.92
!     
      if (ta.ge.0.) then
         alpha=17.27
         beta=237.3
      else
         alpha=21.88
         beta=265.5
      endif
      va=6.11*exp(alpha*ta/(ta+beta))
      vd=6.11*exp(17.27*tb/(tb+237.3))
      psi=0.707 + vd/158.
      rr=1 + (0.25-0.005*(va-vd))*cc**2
!
      ra=ea*sigma*rr*psi*(ta+273)**4
      rs=es*sigma*(ts+273)**4
      rb=ra-rs
      rbkj=rb*3.6
!
      end

!*****************************************************************************

      subroutine evcon(ts,ta,tb,uw,p,re,rc)

!  termini di evaporazione-convezione nel bilancio termico atmosfera-mare
!
!  gam = costante psicrometrica in mbar/K
!  fu = funzione della velocita' del vento in W/(m2 mbar)
!  vpa = pressione parziale di vapore nell'aria  in mbar
!  re, rc = termini di evaporazione e convezione in W/m2
!  rekj, rckj = termine di evaporazione e convezione in kj/(m2 ora)

      implicit none

      double precision ts, ta, tb, uw, p, re, rc

      double precision pa, fu
      double precision gam
      double precision vd, vs, vpa
      double precision rekj, rckj

      pa = p
      gam = 0.66

      vd=6.11*exp(17.27*tb/(tb+237.3))
      vs=6.11*exp(17.27*ts/(ts+237.3))
      vpa= vd - (pa-vd)*(ta-tb)/(1540-1.3*tb)

      fu = 4.4 + 1.82*uw

      re = - fu*(vs-vpa)
      rekj = re*3.6

      rc = - fu*gam*(ts-ta)
      rckj = rc*3.6

      end

!*****************************************************************************

!-----------------------------------------------------------------------------
      end module heat_default
!-----------------------------------------------------------------------------
