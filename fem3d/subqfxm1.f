
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

c heat flux module
c
c contents :
c
c subroutine subtem(dt,dh,qs,uw,cc,ta,tb,ts,tsnew,rtot)
c			compute water temperature in element of depth dh
c subroutine subtem0(dt,dh,qs,uw,cc,ta,tb,ts,tsnew,evap)
c			same as subtem, but use heatlucia
c subroutine heatlucia(t,p,w,tb,cc,ts,qsens,qlat,qlong,evap)
c			computes heat fluxes from Lucia modules
c subroutine longwave(ts,ta,tb,cc,rb)
c			long-wave radiation
c subroutine evcon(ts,ta,tb,uw,p,re,rc)
c			latent heat and convective heat
c
c revision log :
c
c 01.06.1998	ggu&lcz	written from scratch (nearly)
c 24.06.1998	ggu&lcz	subroutines from lucia integrated
c 30.04.2001	ggu	new routine qtotal_tb
c 09.12.2002	ggu	cleaned and re-arranged
c 23.03.2006	ggu	changed time step to real
c 23.03.2010	ggu	changed v6.1.1
c 08.10.2010	ggu	changed VERS_6_1_13
c 16.02.2011	ggu	pstd introduced
c 05.11.2014	ggu	changed VERS_7_0_5
c 24.01.2018	ggu	changed VERS_7_5_41
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c
c notes :
c
c qs	net solar radiation (reflection already subtracted)
c ta	air temperature
c tb	wet bulb temperature
c uw	wind speed
c cc	cloud cover (0 clear sky, 1 totally covered)
c
c*****************************************************************************

      subroutine subtem(dt,dh,qs,uw,cc,ta,tb,ts,tsnew,evap)

c  Calcola la temperatura dell'acqua in un bacino di profondita' dh,
c  risolvendo l'equazione del bilancio termico all'interfaccia acqua-aria.
c
c  dt = timestep del modello FEM
c  dh = profondita' dello strato
c  qs = radiazione solare (short wave radiation) entrante nella superficie 
c       dell'acqua (radiazione incidente al netto della riflessione)
c  uw = velocita' del vento
c  cc = copertura nuvolosa (tra 0 e 1)
c  ts,ta,tb = temperatura dell'acqua, dell'aria, temp. dell'aria a bulbo
c             bagnato
c  tsnew = nuova temperatura dell'acqua
c
c  pa = pressione atmosferica in mbar in subroutine evcon 
c	(!!!!! verificare la sensibilita' a pa)
c
c  cpw = calore specifico dell'acqua di mare (J/kg C) da Gill, per ts=16 C
c  rhow = densita' dell'acqua di mare (kg/m3) da Gill, per ts=16 C
c
c  rtot = total radiation
c  qs = solar radiation (short wave)
c  rb = long wave radiation
c  re = evaporation term
c  rc = convection term 
c
c  Formulazione del bilancio termico da:
c  Mariotti M. e C. Dejak, 1982: 'Valutazione dei flussi energetici
c  all'interfaccia acqua-aria nella Laguna di Venezia', Ambiente Risorse 10.
c
c  Vengono usate le unita' SI.
c
c  L. Zampato - Dicembre 1997

      implicit none

      include 'subqfxm.h'

      real dt
      real dh
      real p
      real qs, uw, cc
      real ta, tb, ts, tsnew
      real evap

      real pstd
      parameter ( pstd = 1013.25 )

      real rb, re, rc, rtot
      real ct

c constants

      p = pstd

c  long-wave term
c   input: ts,ta,tb,cc
c   output: rb = long wave radiation

      call longwave(ts,ta,tb,cc,rb)

c  evaporation and convection terms
c   input: ts,ta,tb,uw
c   output: re,rc = evaporation and convection terms

      call evcon(ts,ta,tb,uw,p,re,rc)

c  total radiation: rtot (W/m2) (positive if into water)

      rtot = qs+rb+re+rc

c  evaporation [kg/(m**2 s)]

      evap = re / 2.5e+6                !divide by latent heat of evaporation
      evap = evap / rhow                 !in [m/s]
      evap = evap * 1000. * 86400.      !in [mm/day]

c  heat capacity/area

      ct = cpw*rhow*dh

c  new temperature      formula:  dQ = dT * rho * cpw * dh / dt

      tsnew = ts + rtot*dt/ct

      end

c*****************************************************************************

      subroutine subtem0(dt,dh,qs,uw,cc,ta,tb,ts,tsnew,evap)

c same as subtem, but use heatlucia

      implicit none

      include 'subqfxm.h'

      real dt
      real dh
      real qs, uw, cc
      real ta, tb, ts, tsnew
      real evap

      real pstd
      parameter ( pstd = 1013.25 )

      real rtot
      real qsens,qlat,qlong
      real ct, p

c constants

      ct = cpw*rhow*dh	!heat capacity	[ J / (m**2 K) ]
      p = pstd

      call heatlucia(ta,p,uw,tb,cc,ts,qsens,qlat,qlong,evap)

      rtot = qs - ( qlong + qlat + qsens )

c  new temperature      formula:  dQ = dT * rho * cpw * dh / dt

      tsnew = ts + rtot*dt/ct

      evap = -evap			!in [kg/(m**2 s)]
      evap = evap / rhow                 !in [m/s]
      evap = evap * 1000. * 86400.      !in [mm/day]

      end

c*****************************************************************************

        subroutine heatlucia(t,p,w,tb,cc,ts,qsens,qlat,qlong,evap)

c computes heat fluxes from Lucia modules
c
c heat fluxes are positive upward (from sea to atmosphere)

        implicit none

        real t          !air temperature [C]                    - in
        real p          !pressure [mb]                          - in
        real w          !wind speed [m/s]                       - in
        real tb         !wet bulb temperature [C]               - in
        real cc         !cloud cover [0-1]                      - in
        real ts         !sea temperature [C]                    - in
        real qsens      !sensible heat flux [W/m**2]            - out
        real qlat       !latent heat flux [W/m**2]              - out
        real qlong      !long wave radiation [W/m**2]           - out
        real evap       !evaporation [kg/(m**2 s)]              - out

	real ta,uw
	real rb,re,rc

	ta = t		!t air
	uw = w		!wind speed

        if( p > 10000 ) stop 'error stop heatlucia: p not in mbar'

	call longwave(ts,ta,tb,cc,rb)
	call evcon(ts,ta,tb,uw,p,re,rc)

	qlong = -rb
	qsens = -rc
	qlat  = -re
	evap  = -re / 2.5e+6

	end

c*****************************************************************************

      subroutine longwave(ts,ta,tb,cc,rb)

c  termine di radiazione long-wave nel bilancio termico atmosfera-mare
c
c  es, ea = emissivita' del mare e dell'aria
c  v = tensione di vapore alla temperatura t
c  ra, rb = radiazione emessa dall'atmosfera e dal mare
c  rb = termine long-wave in W/m2 da inserire nel bilancio termico
c  rbkj = termine long-wave in kJ/(m2 ora)
c
      implicit none
c
      real ts, ta, tb
      real cc
      real sigma, es, ea
      real alpha, beta
      real va, vd, psi, rr
      real ra, rs, rbkj, rb
c
      sigma=5.67/(10**(8))
      es=0.97
      ea=0.92
c     
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
c
      ra=ea*sigma*rr*psi*(ta+273)**4
      rs=es*sigma*(ts+273)**4
      rb=ra-rs
      rbkj=rb*3.6
c
      end

c*****************************************************************************

      subroutine evcon(ts,ta,tb,uw,p,re,rc)

c  termini di evaporazione-convezione nel bilancio termico atmosfera-mare
c
c  gam = costante psicrometrica in mbar/K
c  fu = funzione della velocita' del vento in W/(m2 mbar)
c  vpa = pressione parziale di vapore nell'aria  in mbar
c  re, rc = termini di evaporazione e convezione in W/m2
c  rekj, rckj = termine di evaporazione e convezione in kj/(m2 ora)

      implicit none

      real ts, ta, tb, uw, p, re, rc

      real pa, fu
      real gam
      real vd, vs, vpa
      real rekj, rckj

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

c*****************************************************************************

