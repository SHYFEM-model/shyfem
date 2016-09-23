      !!------------------------------------------------------------------------
      !! ivan          
      !! Heat fluxes calculated through air-sea parametrization 
      !! with bulk formulae (Pettenuzzo et al., 2010)
      !!       
      !! Drag coefficient parametrization in heat fluxes 
      !! with Hellermann and Rosenstein (1983)
      !!       
      !!------------------------------------------------------------------------
      !!
      !! Written starting from subqfxm3.f --> heatareg
      !!
      !! subroutine heatmfsbulk --> convert variables to needed units
      !! subroutine mfs_qfluxes --> calculate qsens, qlat, qlong, evap, cd
      !! subroutine qshort1     --> calculate qshort = qswa = qs
      !! 
      !!----------------------------------------------------------------------
      !!
      !!    Method: 
      !!    CALL heatmfsbulk with:
      !!                          iheat=8
      !!                          itdrag=4 
      !!    READ Atmospheric data from qflux and wind files:
      !!    - qflux [4 vars]: QS (not-used) + T2M (C) + D2M (C) + TCC (0-1)
      !!      (we keep 4vars but QS not called)
      !!    - wind  [3 vars]: U10M (m/s) + V10M (m/s) + MSL (Pa)
      !!
      !!   Needed conversion:
      !!   - sstk: sea surface temperature   (C)-->(Kelvin)
      !!   - tairk: 2m air temperature       (C)-->(Kelvin)
      !!   - dew: 2m Dew point Temperature   (C)-->(Kelvin)
      !!   Other variables:
      !!   - pa: Mean Sea Level Pressure     (coverted in hPa in submeteo2.f) 
      !!   - uuv,vvw: 10m Wind Velocities    (m/s)
      !!   - uub,vvb: sea surface velocities (m/s) 
      !!   - cc: total cloud cover           (0-1)
      !!
      !!   Computes:
      !!   - Solar Radiation using Reed formula (1975, 1977)
      !!   - Net Long wave radiation using Bignami et al. (1995)
      !!   - Latent and Sensible heat using Kondo (1975)
      !!   - Drag coeff using Hellerman and Rosenstein (1983)
      !!   - Solar Radiation according Astronomical formulae (Reed, 1975, 1977)
      !!
      !!-----------------------------------------------------------------------



        subroutine heatmfsbulk(days,im,ih,ddlon,ddlat,ta,pa,uuw,vvw,dp,
     +                   cc,tm,uub,vvb,qsens,qlat,qlong,evap,qswa,
     +                   cd)  

        implicit none

        real ddlon,ddlat 
        real ta 
        real p   
        real uuw,vvw 
        real cc 
        real dp 
        real tm  
        real uub,vvb 
        real qsens,qlat,qlong,evap
        real qswa  
        real cd 

        real dew,sstk,tairk,norspeed,pa
        integer days,im,ih
        real, parameter  :: ckelv = 273.15


        ! convert in units needed by MFS-bulk-formulae!  
        sstk = tm + ckelv      
        tairk = ta + ckelv   
        dew = dp + ckelv    
        norspeed =sqrt (uub*uub + vvb*vvb)   

        call mfs_qfluxes(ddlon,ddlat,
     +                   days,ih,im,
     +                   sstk,uuw,vvw,
     +                   tairk,dew,cc,
     +                   uub,vvb,norspeed,pa,
     +                   qswa,
     +                   qsens,qlat,qlong,evap,
     +                   cd)


        end

        !---------------------------------------------------

        subroutine mfs_qfluxes(ddlon,ddlat,
     +                   days,ih,im,
     +                   sstk,uuw,vvw,
     +                   tairk,dew,cc,
     +                   uub,vvb,norspeed,pa,
     +                   qswa,
     +                   qsens,qlat,qlong,evap,
     +                   cd)

        implicit none

       ! surface air pressure, expsi, dry air gas constant
       real, parameter  :: ps = 1013.25
       real, parameter  :: expsi = 0.622
       real, parameter  :: rd = 287.

       ! specific heat capacity
       real, parameter  :: cp = 1005.
       real, parameter  :: ckelv = 273.15

       integer days,ih,im
       real ddlon,ddlat
       integer kku
       real uuw,vvw 
       real, parameter  :: stefan = 5.67e-8
       real, parameter  :: emic = 0.97

       real,parameter, dimension(6) :: ahr =
     +           (/ 0.934e-3, 0.788e-4, 0.868e-4,
     +              -0.616e-6,-.120e-5,-.214e-5 / )

       real pa         
       real ea
       real ch,ce,fh,fe
       real sh_now,shms,cseep
       real vtnow,sstk,dew
       real cd
       real uub,vvb 
       real rhom
       real cc 
       real wair
       real deltemp,s,stp
       real esre
       real qswa,qsens,qlong,qlat,evap
       real rel_u,rel_v,rspeed,norspeed
       real tairk

       real,parameter, dimension(5) :: a_h =
     +           (/0.0,0.927,1.15,1.17,1.652/)
       real,parameter, dimension(5) :: a_e =
     +           (/0.0,0.969,1.18,1.196,1.68/)
       real,parameter, dimension(5) :: b_h =
     +           (/1.185,0.0546,0.01,0.0075,-0.017/)
       real,parameter, dimension(5) :: b_e =
     +           (/1.23,0.0521,0.01,0.008,-0.016/)
       real,parameter, dimension(5) :: c_h =
     +           (/0.0,0.0,0.0,-0.00045,0.0/)
       real,parameter, dimension(5) :: c_e =
     +           (/0.0,0.0,0.0,-0.0004,0.0/)
       real,parameter, dimension(5) :: p_h =
     +           (/-0.157,1.0,1.0,1.0,1.0/)
       real,parameter, dimension(5) :: p_e =
     +           (/-0.16,1.0,1.0,1.0,1.0/)
       real,parameter  :: onsea = 0.98
       real,parameter  :: par1 = 640380.
       real,parameter  :: par2 = -5107.4



        !sst = sstk - ckelv !we need also sst in [C]

        ! Calculate Specific Humidity   -----------
        !------------------------------------------

        sh_now = (1/1.22)*onsea*par1*EXP(par2/dew)

        !---------------   SHORT WAVE   ----------
        qswa = 0.        
        call qshort1(im,days,ih,ddlat,ddlon,cc,qswa)

        ! WIND INPUT ---------------------------------

        !norspeed = speed
        rel_u = uuw - uub !true stress - u
        rel_v = vvw - vvb !true stress - v
        rspeed = sqrt( rel_u*rel_u + rel_v *rel_v)

        ! QLONG   -------------------------------
        wair = sh_now / (1 - sh_now)

        !calculates the virtual temperature of air
        vtnow = (tairk*(expsi+wair))/(expsi*(1.+wair))

        !calculates the density of the moist air
        rhom = 100.*(ps/rd)/vtnow

        ea = (wair / (wair+0.622) ) * pa

        qlong = emic*stefan*( sstk**4. )
     +       - ( stefan*( tairk**4. ) * ( 0.653 + 0.00535*ea ) )
     +           * ( 1. + 0.1762*( cc * cc ) )


       ! QSENS   -------------------------------

        !calculates the term :      ( Ts - Ta )
        deltemp = sstk - tairk

       !variable turbulent exchange coefficients ( from Kondo 1975 )
       !calculate S :

             s=deltemp/(norspeed**2.)

       !calculate the Stability Parameter :

             stp=s*abs(s)/(abs(s)+0.01)

         fe = 0.
         fh = 0.

       !calculate fe,fh :
            IF (s.lt.0. .and. ((stp.gt.-3.3).and.(stp.lt.0.))) THEN
                fh = 0.1+0.03*stp+0.9*exp(4.8*stp)
                fe = fh
            ELSE IF (s.lt.0. .and. stp.le.-3.3) THEN
                fh = 0.
                fe = fh
            ELSE                       ! --- for unstable condition 
                fh = 1.0+0.63*sqrt(stp)
                fe = fh
            ENDIF

         !calculate kku 
            IF (norspeed >= 0. .AND. norspeed <= 2.2)       THEN
                kku=1
            ELSE IF (norspeed > 2.2 .AND. norspeed <= 5.0)  THEN
                kku=2
            ELSE IF (norspeed > 5.0 .AND. norspeed <= 8.0)  THEN
                kku=3
            ELSE IF (norspeed > 8.0 .AND. norspeed <= 25.0) THEN
                kku=4
            ELSE IF (norspeed > 25.0 )                   THEN
                kku=5
	    ELSE
                kku=5
            ENDIF


       !calculate ce,ch :
            ch = ( a_h(kku) + b_h(kku) * norspeed ** p_h(kku)
     +                        + c_h(kku) * (norspeed - 8 )**2) * fh

            ce = ( a_e(kku) + b_e(kku) * norspeed ** p_e(kku)
     +                        + c_e(kku) * (norspeed - 8 )**2) * fe

            ch = ch / 1000.0
            ce = ce / 1000.0

            IF (norspeed < 0.3) THEN
               ch = 1.3e-03 * fh
               ce = 1.5e-03 * fe
            ELSE IF(norspeed > 50.0) THEN
               ch = 1.25e-03 * fh
               ce = 1.30e-03 * fe
            ENDIF
       

       !calculate qsens :
          qsens = rhom*cp*ch*norspeed*deltemp

       ! EVAP / QLAT   ---------------------------

       !(iv)  latent heat
       !calculates the LATENT HEAT FLUX  ( watt/m*m )
       !ELAT = L*rho*Ce*|V|*[qs(Ts)-qa(t2d)]

            shms = (1/1.22)*onsea*par1*EXP(par2/sstk) !- Saturation Specific Humidity            

            esre  = shms - sh_now  ! --- calculates the term : qs(Ta)-qa(t2d)

            cseep = ce * norspeed * esre     ! --- calculates the term: Ce*|V|*[qs(Ts)-qa(t2d)]
            evap = (cseep * rhom)  ! in [kg/m2/sec] !! ---calculates the EVAPORATION RATE [m/yr]

            qlat = rhom * cseep *
     +             (2.5008e+6 -2.3e+3 * (sstk - ckelv))


        ! WIND STRESSES  --------------------------
        ! HellerRose
        cd = ahr(1) + ahr(2)*norspeed + ahr(3)*deltemp 
     +       + ahr(4) *norspeed*norspeed 
     +       + ahr(5)*deltemp*deltemp + ahr(6)*norspeed*deltemp

        return

        end

        !**************************************************************************** 
        !**************************************************************************** 
        !**************************************************************************** 
        !**************************************************************************** 
         
        subroutine qshort1(im,days,ih,ddlat,ddlon,cc,qswa) 

        implicit none

        real, parameter  :: pi = 3.1415927
        real, parameter  :: degrad = pi/180. 
        real, parameter  :: degradr = 180./pi 
        real, parameter  :: eclips = 23.439*degrad 
        real, parameter  :: solar = 1350.
        real, parameter  :: tau = 0.7
        real, parameter  :: yrdays = 360.

        !Albedo monthly values from Payne (1972) as means of the values
        !at 40N and 30N for the Atlantic Ocean (hence the same latitudinal
        !band of the Mediterranean Sea) :

        real,parameter, dimension(12) :: alpham =
     +            (/0.095,0.08,0.065,0.065,0.06,0.06,0.06,0.06,
     +              0.065,0.075,0.09,0.10/)

        integer im,ih,days
        real ddlat,ddlon,cc
        real th0,th02,th03
        real sundec,thsun 
        real alat,alon 
        real coszen,qatten 
        real rtk
        real qzer,qdir,qdiff,qtot
        real tjul
        real sunbet,sunbetd 
        real albedo 
        real qswa 
        real, parameter  :: thco = 5.035D-04
        real, parameter  :: aozone = 0.09  

        !-------- calculations start -------------------------------

        th0 = 2.*pi*days/yrdays
        th02 = 2.*th0
        th03 = 3.*th0

        !sun declination:
        sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0) -
     +  0.006758*cos(th02) + 0.000907*sin(th02) -
     +  0.002697*cos(th03) + 0.001480*sin(th03)

        alat=ddlat*degrad
        alon=ddlon*degrad

        thsun = (ih -12.)*15.*degrad + alon

        !cosine of the solar zenith angle :

        coszen =sin(alat)*sin(sundec)+cos(alat)*cos(sundec)*cos(thsun)

        if (coszen .le. thco) then
          coszen = 0.0
          qatten = 0.0
        else
          qatten = tau**(1./coszen)
        end if

        rtk=(1.0+1.67E-2*cos(pi*2.*(days-3.0)/365.0))**2
        qzer  = coszen * solar * rtk
        qdir  = qzer * qatten
        qdiff = ((1.-aozone)*qzer - qdir) * 0.5
        qtot  =  qdir + qdiff

       !In the sunbet formula enters days=julian days of the
       !current year (1,365); we subtract 81 and not 82 because
       !days = day -1.
       !conversion of (days - 81) to radians

        tjul = (days-81.)*degrad

       !sin of the solar noon altitude in radians :

       sunbet=sin(alat)*sin(eclips*sin(tjul)) +
     +        cos(alat)*cos(eclips*sin(tjul))

       !solar noon altitude in degrees :

       sunbetd = degradr*asin(sunbet)

       !correction for cloud less than 0.3 according to Reed

       albedo = alpham(im)

       if (cc.lt.0.3) then
          qswa = qtot*(1.-albedo)
       else
          qswa  = qtot*(1-0.62*cc + .0019*sunbet)*(1.-albedo)
       endif


       return
        end



