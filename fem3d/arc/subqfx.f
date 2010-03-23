c
c $Id: subqfx.f,v 1.7 2001/11/16 07:35:43 georg Exp $
c
c heat flux module
c
c contents :
c
c subroutine subtem(idt,dh,qs,uw,cc,ta,tb,ts,tsnew,rtot)
c			compute water temperature in element of depth dh
c subroutine longwave(ts,ta,tb,cc,rb)
c			long-wave radiation
c subroutine evcon(ts,ta,tb,uw,re,rc)
c			latent heat and convective heat
c
c revision log :
c
c 01.06.1998	ggu&lz	written from scratch (nearly)
c 24.06.1998	ggu&lz	subroutines from lucia integrated
c 30.04.2001	ggu	new routine qtotal_tb
c
c comments :
c
c qs	net solar radiation (reflection already subtracted)
c ta	air temperature
c tb	wet bulb temperature
c uw	wind speed
c cc	cloud cover (0 clear sky, 1 totally covered)
c
c*****************************************************************************

	subroutine qsun(year,month,day,hour,cc,qs)

c computes solar radiation incident to sea water
c
c uses only date, time and cloud cover
c
c Gill, p.34

	implicit none

	integer year	!year (fully qualified, e.g., 1995)
	integer month	!month [1,12]
	integer day	!day [1,28-31]
	real hour	!hour [0-24)		!ATTENTION!! -> REAL
	real cc		!cload cover [fractional] (0-1)
	real qs		!solar radiation into water [W/m**2]

	real qtopat	!solar radiation at top of atmosphere [W/m**2]
	real albedo	!mean albedo [fractional]
	real absorb	!average absorption of atmosphere [fractional]

	real ad		!declination of sun [rad]
	real rr		!?? (correction??)
	real ai		!angle of incidence (with respect to zenith) [rad]
	real cosi	!cos(ai)

	albedo = 0.06
	absorb = 0.25

	call decl(year,month,day,ad,rr)

	call incid(hour,ad,ai,cosi)
	call radmax(rr,cosi,qtopat)

	qs = qtopat * (1.-absorb) * (1.-albedo) * (1.-0.7*cc)

	qs = max(qs,0.)

c	write(6,*) month,day,hour,qtopat,qs

	end

c*****************************************************************************

        subroutine qtotal_tb(qs,uw,cc,twb,ta,ts,qtot)

c computes total heat flux [W/m**2] -> heat fluxes positive if into water

	implicit none

	real qs		!solar heat flux incident into water [W/m**2]
	real uw		!wind velocity [m/s]
	real cc		!cload cover [fraction] (0-1)
	real twb	!wet bulb temperature [C]
	real ta		!air temperature [C]
	real ts		!sea temperature [C]
	real qtot	!total flux into water [W/m**2] (return)

	real rb,re,rc

        call longwave(ts,ta,twb,cc,rb)	!long wave radiation rb

        call evcon(ts,ta,twb,uw,re,rc)	!evaporation and convection (re,rc)

        qtot = qs + rb + re + rc

	end

c*****************************************************************************

        subroutine qtotal(qs,uw,cc,rh,ta,ts,qtot)

c computes total heat flux [W/m**2] -> heat fluxes positive if into water

	implicit none

	real qs		!solar heat flux incident into water [W/m**2]
	real uw		!wind velocity [m/s]
	real cc		!cload cover [fraction] (0-1)
	real rh		!relative humidity [%] (0-100)
	real ta		!air temperature [C]
	real ts		!sea temperature [C]
	real qtot	!total flux into water [W/m**2] (return)

	real twb	!wet bulb temperature [C]
	real rb,re,rc

	call rh2twb(ta,rh,twb)		!compute wet bulb temperature

        call longwave(ts,ta,twb,cc,rb)	!long wave radiation rb

        call evcon(ts,ta,twb,uw,re,rc)	!evaporation and convection (re,rc)

        qtot = qs + rb + re + rc

	end

c*****************************************************************************

        subroutine theat(dt,dh,qtot,ts,tsnew)

c computes new temperature due to heating

	implicit none

	real dt		!time interval for heating (time step) [sec]
	real dh		!depth of fluid layer [m]
	real qtot	!total flux into water [W/m**2] (return)
	real ts		!sea temperature [C]
	real tsnew	!new sea temperature after heating [C]

	real cw		!specific heat of water [J/(kg C)]
	real rhow	!density of sea water [kg/m**3]

        cw = 3991
        rhow = 1026.

c  new temperature      formula:  dQ = dT * rho * cw * dh / dt

        tsnew = ts + qtot * dt / ( cw * rhow * dh )

	end

c*****************************************************************************
c*****************************************************************************
c*****************************************************************************

c-----------------------------------------------------------------------------

      subroutine subtem(idt,dh,qs,uw,cc,ta,tb,ts,tsnew,rtot)

c-----------------------------------------------------------------------------
c
c  Calcola la temperatura dell'acqua in un bacino di profondita' dh,
c  risolvendo l'equazione del bilancio termico all'interfaccia acqua-aria.
c
c  idt = timestep del modello FEM
c  dh = profondita' dello strato
c  qs = radiazione solare (short wave radiation) entrante nella superficie 
c       dell'acqua (radiazione incidente al netto della riflessione)
c  uw = velocita' del vento
c  cc = copertura nuvolosa (tra 0 e 1)
c  ts,ta,tb = temperatura dell'acqua, dell'aria, temp. dell'aria a bulbo
c             bagnato
c  tsnew = nuova temperatura dell'acqua
c
c  pa = pressione atmosferica in mbar: imposta costante (pa = 1013) in
c       subroutine evcon (!!!!! verificare la sensibilita' a pa)
c
c  cw = calore specifico dell'acqua di mare (J/kg C) da Gill, per ts=16 C
c  row = densita' dell'acqua di mare (kg/m3) da Gill, per ts=16 C
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
c
c-----------------------------------------------------------------------------

      implicit none

      integer idt
      real dh
      real qs, uw, cc
      real ta, tb, ts, tsnew

      real rb, re, rc, rtot
      real cw, row, ct

c constants

      cw = 3991
      row = 1026.

c  long-wave term
c   input: ts,ta,tb,cc
c   output: rb = long wave radiation

      call longwave(ts,ta,tb,cc,rb)

c  evaporation and convection terms
c   input: ts,ta,tb,uw
c   output: re,rc = evaporation and convection terms

      call evcon(ts,ta,tb,uw,re,rc)

c  total radiation: rtot (W/m2) (positive if into water)

      rtot = qs+rb+re+rc

c  heat capacity/area

      ct = cw*row*dh

c  new temperature      formula:  dQ = dT * rho * cw * dh / idt

      tsnew = ts + rtot*idt/ct

      end

c-----------------------------------------------------------------------------

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
c     write(*,*) 'ts,ta,tb', ts,ta,tb
c     write(*,*) 'sigma ', sigma
c     write(*,*) 'psi, rr ', psi,rr
c     write(*,*) 'va, vd ', va,vd
c     write(*,*) 'ra, rs, ', ra,rs
c     write(*,*) 'rb, rbkj, ', rb,rbkj
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine evcon(ts,ta,tb,uw,re,rc)
c
c  termini di evaporazione-convezione nel bilancio termico atmosfera-mare
c
c  gam = costante psicrometrica in mbar/K
c  fu = funzione della velocita' del vento in W/(m2 mbar)
c  vpa = pressione parziale di vapore nell'aria  in mbar
c  re, rc = termini di evaporazione e convezione in W/m2
c  rekj, rckj = termine di evaporazione e convezione in kj/(m2 ora)
c
      implicit none
c
      real ts, ta, tb
      real pa, uw, fu
      real gam
      real vd, vs, vpa
      real re, rekj, rc, rckj

      pa = 1013
c
c     write(*,*) 'uw, pa, ', uw, pa
      gam = 0.66
c
      vd=6.11*exp(17.27*tb/(tb+237.3))
      vs=6.11*exp(17.27*ts/(ts+237.3))
      vpa= vd - (pa-vd)*(ta-tb)/(1540-1.3*tb)
c
      fu = 4.4 + 1.82*uw
c
      re = - fu*(vs-vpa)
      rekj = re*3.6
c
      rc = - fu*gam*(ts-ta)
      rckj = rc*3.6
c
c     write(*,*) 're, rc, ', re,rc
c
      return 
      end
c
c-------------------------------------------------------------------------------
