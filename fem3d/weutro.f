
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

c weutro - EUTRO box model
c
c contents :
c
c revision log :
c
c 01.01.2000	ggu	original version until April 2002
c 19.04.2002	ggu&dmk	BUGFIX1, print_init, changes in ZOOP
c 20.06.2003	ggu&dmk	new routine for sediments
c 19.12.2003	ggu	new routines param_taranto and param_venezia
c 19.12.2003	isa	new routine rintens_tar()
c 18.02.2004	dmk	dtl non viene passato nel reattore dtl=300/86400
c 26.02.2004	dmk	BUGFIX 2D 3D for depth > 1 (???) (LIGHTFIX)
c 03.03.2004	ggu	in rlim new var btaranto (FIX)
c 04.03.2004	ggu	changes from donata integrated
c 10.06.2004	ggu	new routine param_read for Michol
c 09.08.2004	ggu	marked changes from Donata with GGU
c 17.08.2004	ggu	new routines for debug (eutro_replay)
c 21.08.2004	ggu	rearangments, renaming ($LGA)
c 24.08.2004	ggu	all changes incorporated (see check history 1-11)
c 30.08.2004	ggu	cleanup of settopseg, setbotseg
c 15.02.2006	ggu&mcg	some comments inserted for denitrification (SK18D,SK14D)
c 23.03.2010	ggu	changed v6.1.1
c 08.10.2010	ggu	changed VERS_6_1_13
c 16.03.2012	ggu	dummy restart routines added
c 24.09.2013	dmk	insert direct call to qflux input file from shyfem
c 28.09.2014	dmk	PNO3G1 cancelled
c 10.07.2015	ggu	changed VERS_7_1_50
c 01.04.2016	ggu	changed VERS_7_5_7
c 07.06.2016	ggu	changed VERS_7_5_12
c 17.06.2016	dmk	deleted rlux, link for shyfem 7_5_13
c 27.06.2016	ggu	changed VERS_7_5_16
c 31.08.2018	ggu	changed VERS_7_5_49
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c
c notes :
c
c 1
c       debug routines
c       nutlim
c       in steele temp1 = ke and not ke*depth
c       zood ?
c       wsedim, nempar (?)
c
c 2
c       changes signed with GGU
c       ddin1, ddin2, prod, cons
c       loicz()
c       wsedim
c
c 3
c       no denit (?)
c
c 5
c       denit set again
c
c 6
c       new array cold
c       eutro_check
c       param_read
c       wsedim -> cold
c
c 7
c       changes to eutro_check
c       deleted rdtemp()
c       eutro_replay, read/write_debug_param
c       ample re-formatting
c
c 8
c      new default of SOD1D
c      deleted rlim(), rewritten steele
c      changed rintens()
c      deleted rintens_tar()
c
c 9
c       renamed param_init, param_print
c       in source new rlux in arguments
c          set rlghts, call rintens()
c       cleaned up reaeration part
c       consider light attenuation in ditoro()
c       new weutro_test as a test drive for weutro
c       revised steele, rintens
c       moved luxlen, it2n to weutro_light.f
c       moved loicz, wsedim, tempcoef to weutro_sedim.f
c       new comments, reformatting
c
c 10
c       no changes
c
c 11
c       use id, rlux internally (error_check etc.)
c       new routines get/set_internal_param
c       new routine weutro_debug
c       eliminated ITYPE and segtype()
c       call BENTFLUX only if botseg is true
c       topseg used for reaeration
c       debug -> wdebug : weutro_debug(.true.), weutro_debug(.false.)
c
c
c
c State variables used: (Wasp)
c
c nh3           71      1
c no3           72      2
c opo4          73      3
c phyto         74      4
c cbod          75      5
c do            76      6
c on            77      7
c op            78      8
c zoo           79      9
c
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************

      subroutine param_init

      implicit none
      INCLUDE 'weutro.h'

c initialization of parameters

C       BYPASS OPTIONS FOR SYSTEMS 1-9. 1=BYPASS 0=SIMULATE
      

        PI = 3.14159

c--------------------------------------------------------------
c       please do not change anything in this subroutine below this point
c       all changes to parameters should be done in a custom subroutine
c       see param_venezia or param_taranto for example
c--------------------------------------------------------------

c       start parameter values

c--------------------------------------------------------------
c       general parameters
c--------------------------------------------------------------

      graztype=1        !if graztype=1 simulate zoop. variable
                  !if graztype=0  use wasp formulation
                  !Set initial value of variable zoo=0
      
c      if CCHLX is variable in the segments, introduce each value else

        CCHL=50      !range 20-50
      CCHLX(1) = CCHL

        NUTLIM=1.       !nutrient limitation option 0=minimum, 1=multiplicative.
                  ! Default=0
                        ! GGU new value (was 1) ?

c--------------------------------------------------------------
c       subroutine phyto
c--------------------------------------------------------------

c------ growth and respiration

c      K1RT=1.045       !wasp orig      
c      K1RC=0.125      !wasp orig
      K1RT=1.068      !adapting the curve to dejak model
      K1RC=0.096      !adapting the curve to dejak model
        K1T=1.068      
        K1C=2.               !GGU new value (was 1.5) ?

c      KMPHYT=1.
      KMPHYT=0.

c------ decomposition

        KPZDC=0.02      !verify the value
        KPZDT=1.08      !Default=1.0

c------ nutrients limitation

c        KMNG1=0.025      !for standard model application
                  !use a large KMNG1 (when KMNG1=0 PNH3G1=1.0)
        KMNG1=0.05      !from dejak model Venice lagoon
        KMPG1=0.01       

c      F(10,10,1) !spatially variable dissolved fraction of inorganic P
      FOPO4=0.9 !spatially variable dissolved fraction of inorganic P
         
c      attenzione: nel benthic layer asume valori diversi (0.045-0.001)

c------ grazing - zooplankton as a forcing

        K1G=0.            
c        K1G=0.08            
c        K1D=0.02       !wasp orig
      K1D=0.12      !from dejak model Venice Lagoon 
c        K1D=0.2         !prova 12sett
      ZOO=0.7      
c      if fraction is variable in segments, input each value, else:
      ZOOSG(1)=1.

c------ grazing - zooplankton as a variable

c      KGRZ=1.44      !max grazing rate day-1
        KGRZ=2.         !prova 12sett
      KPHYZ=0.5      !
c      EFF=0.5            !zoo-phyto digestion efficiency
      EFF=0.7
      KDZ=0.192      !zoo death (with excrection)

c--------------------------------------------------------------
c       subroutine organop, inorganp (P=phosphorous)
c--------------------------------------------------------------

        PCRB=0.025 !mgP/mgC
        FOP=0.5            !unitless
        K58T=1.08      !unitless      
        K58C=0.22      !day-1 
        KOPDC=0.0004      !day-1
        KOPDT=1.08      !unitless
      KPO4=1.0      !microg. P/L

c--------------------------------------------------------------
c       subroutine organicn, ammonia, nitrate
c--------------------------------------------------------------

      NCRB=0.115      !mg N/mg C
             FON=0.5     
c            K1320C=0.1      !0.09-0.13 day-1 wasp orig 
      K1320C=0.05      !from dejak model lagoon of Venice
        K1320T=1.08      !unitless
        K140C=0.09      !day-1
        K140T=1.045      !unitless
        KNIT=2.            !mgO2/L 
        KNO3=0.1      !mgO2/L       
        K1013C=0.075    !day-1
        K1013T=1.08    !unitless
        KONDC=0.0004    !day-1     
        KONDT=1.08      !unitles

c--------------------------------------------------------------
c       subroutine CBOD
c--------------------------------------------------------------

      OCRB=32./12.      !mg O2/mg C
        KDC=0.18       !0.16-0.21 day-1
        KDT=1.047      !unitless
        KDSC=0.000      !day-1
        KDST=1.08      !unitless
        KBOD=0.5      !mg N/L

c--------------------------------------------------------------
c       subroutine dissoxyg
c--------------------------------------------------------------

      WIND=3.             !m/s
        AIRTMP=22.           !C
        WTYPE=3.      !1=small 2=medim 3=large
      XICECVR=1.      !default: 1. no ice

        K2=0   !4.3      !4.1-4.7 day-1 if k2=0 then use kawind or kahydra      

c        REARSG(i)=0       !Segment specific reareation rate constant
                        !REARSG  used when  rear is not calc from wind or hydro

                  !FIXME
       SOD1D(1)=0.0      !g/m2-day 0.2-4.0 sediment oxygen demand for segment
       SODTA(1)=1.08 

c--------------------------------------------------------------
c       light limitation
c--------------------------------------------------------------

        LGHTSW=3        !LIGHT SWITCH: 1=Di Toro 2=Smith  3=Steele
c      IS2=50000.      !optimum light intensity for Steele lux/h
      IS2=1200000.      !optimum light intensity for Steele lux/day
      IS2= 20            !prova W/m2/timestep 300 sec!
c--------------------------------------------------------------
c       subroutine ditoro
c--------------------------------------------------------------

        FDAY=0.5
        IS1=300.      !(Ly/day) langleys/day
        IS1X(1)=IS1      !FIXME

c--------------------------------------------------------------
c       subroutine smith
c--------------------------------------------------------------

        PHIMX=720.      !mg C/mole photon
        XKC=0.017      !m2/mg chl a
        ITOT=500.       !ly/day
        iav=0.9 *itot/FDAY
        DTDAY=0.5           !??????verificare, serve a modulare il seno: day
        NEWDAY=1      !unknown????if G.E. 1 calcolo smith

c--------------------------------------------------------------
c       light attenuation
c--------------------------------------------------------------

c you can have up to 5 KE time functions here default is 1

      KE(1)=1.
      KEFN(1)=1
      KESG(1)=1.5      !m-1 range: 0.1-5

      end

c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************


c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************

        subroutine eutro0d(id,t,dt,vol,depth,vel,uws,stp,sal,qss
     &				,c,loads)

c EUTRO 0-Dimensional

        implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )

        integer id              !id of box
        real t                  !actual time [day]
        real dt                 !time step [day]
        real vol                !volume [m**3]
        real depth              !depth of box [m]
        real vel                !velocity [m/s]
        real uws                !wind velocity [m/s]
        real stp                !temperature [C]
        real sal                !salinity [psu] == [per mille]
        real qss                !solar radiation [W/m**2]
        real c(nstate)          !state variable [mg/L] == [g/m**3]
        real loads(nstate)      !loading for c [g/(m**3 day)]

        real cold(nstate)       !old state variable (for diagnostic purpose)
        real cds(nstate)        !source term (right hand side) [g/day]


        call source(id,t,vel,uws,stp,sal,vol,depth,qss,c,cds)
        call load0d(cds,loads,vol)
        call euler(nstate,dt,vol,c,cold,cds)
        call eutro_check(id,nstate,t,dt,vel,stp,sal,qss,vol,depth
     +          ,c,cold,cds,loads)

        end

c********************************************************************

      subroutine load0d(cds,loads,vol)

c integrate loadings

      implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )

      real cds(nstate)      !source term [g/day]
      real loads(nstate)      !loading for c [g/(m**3 day)]
      real vol            !volume of box [m**3]

      integer i
        
      do i=1,nstate
        cds(i) = cds(i) + vol * loads(i)
      end do

      end

c********************************************************************

      subroutine euler(nstate,dt,vol,c,cold,cds)

c new c is computed
c cold is returned which is just c before call

      implicit none

      integer nstate            !number of state variables
      real dt                  !time step [day]
      real vol            !volume [m**3]
      real c(nstate)            !state variable [mg/L] == [g/m**3]
      real cold(nstate)      !old state variable (return)
      real cds(nstate)      !source term [g/day]

      integer i
      real volold,volnew      !volume [m**3]
      real mass            !mass [g]
      real mder            !derivative [g/day]

      volold = vol
      volnew = vol

      do i=1,nstate
          cold(i) = c(i)
        mass = c(i) * volold
        mder = cds(i)
        c(i) = ( mass + dt * mder ) / volnew
        !if(c(i).lt.0.00001) c(i)=0.00001
      end do

      end
      
c********************************************************************

        subroutine eutro_check(id,nstate,t,dt,vel,stp,sal,qss
     +          ,vol,depth
     +          ,c,cold,cds,loads)

c checks result of integration for negatives

        implicit none

      integer id          !id of box
      integer nstate      !number of state variables
      real t              ![day]
      real dt             !time step [day]
      real vel            ![m/s]
      real stp            ![Celsius]
      real sal            ![g/L], e.g., 30.
      real qss            !solar rad from shyfem
      real vol            ![m**3]
      real depth          ![m]
      real c(nstate)      ![mg/L]
      real cold(nstate)   !old state variable (return)
      real cds(nstate)    !source term [g/day]
      real loads(nstate)  !loading for c [g/(m**3 day)]

        integer i
        logical berror

        berror = .false.

        if( vol .lt. 0. ) berror = .true.
        if( depth .lt. 0. ) berror = .true.

       do i=1,nstate
          if( c(i) .lt. 0. ) berror = .true.
        end do

        if( berror ) then
          write(6,*) '*** eutro_check: error in eutro0d'
          call write_debug_param(id,nstate,t,dt
     +                ,vol,depth,vel,stp,sal,qss,c,cold,cds,loads)
          write(6,*) '--- eutro_check: end of error message'

          write(6,*) 'id=', id
          write(6,*) 't=', t
          write(6,*) 'dt=', dt
          write(6,*) 'vol=', vol
          write(6,*) 'depth=', depth
          write(6,*) 'vel=', vel
          write(6,*) 'stp=', stp
          write(6,*) 'sal=', sal
          write(6,*) 'qss=', qss
          write(6,*) 'c=', c
          write(6,*) 'cold=', cold
          write(6,*) 'cds=', cds
          write(6,*) 'loads=', loads
          write(6,*)
          stop 'error stop eutro_check: negative values'
          do i=1,nstate
            if( c(i) .lt. 0. ) c(i) = 0.00001
          end do
        end if

        end

c********************************************************************

      subroutine eutroini

c initializes eutro routines

      implicit none

        call param_init
        !call param_taranto
        call param_venezia

        !call param_read        !GGU new for Michol

        call EUTROINT
      call param_print

      end

c********************************************************************

      subroutine param_print

c prints parameters for weutro

      implicit none

      INCLUDE 'weutro.h'

      write (6,*) 'graztype', graztype 
      write (6,*) 'CCHL' ,CCHL
      write (6,*) 'NUTLIM',NUTLIM
      write (6,*)  'growth and respiration'                                
      write (6,*) 'K1RT',K1RT
      write (6,*) 'K1RC',K1RC
      write (6,*) 'K1T',K1T
      write (6,*) 'K1C',K1C
      write (6,*) 'KMPHYT',KMPHYT
      write (6,*) 'decomposition'
      write (6,*) 'KPZDC',KPZDC
      write (6,*) 'KPZDT',KPZDT
      write (6,*) 'nutrients limitation' 
      write (6,*) 'KMNG1', KMNG1
      write (6,*) 'KMPG1', KMPG1
      write (6,*) 'FOPO4',FOPO4
      write (6,*) 'grazing'
      write (6,*) 'graztype',graztype
      write (6,*) 'wasp orig, K1G', K1G
      write (6,*) 'wasp orig, K1D', K1D
      write (6,*) 'wasp orig, ZOO', ZOO
      write (6,*)  'wasp orig, ZOOSG',ZOOSG(1) 
      write (6,*) 'KGRZ',KGRZ
      write (6,*) 'KPHYZ',KPHYZ
      write (6,*) 'EFF',EFF
      write (6,*) 'KDZ',KDZ
      write (6,*) 'subroutine organop, inorganp'
      write (6,*) 'PCRB',PCRB
      write (6,*) 'FOP',FOP
      write (6,*) 'K58T',K58T
      write (6,*) 'K58C',K58C
      write (6,*) 'KOPDC',KOPDC
      write (6,*) 'KOPDT',KOPDT
      write (6,*) 'KPO4',KPO4
      write (6,*) 'subroutine organicn, ammonia,nitrate' 
      write (6,*) 'NCRB',NCRB
      write (6,*) 'FON',FON
      write (6,*) 'K1320C',K1320C
      write (6,*) 'K1320T',K1320T
      write (6,*) 'K140C',K140C
      write (6,*) 'K140T',K140T
      write (6,*) 'KNIT',KNIT
      write (6,*) 'KNO3',KNO3
      write (6,*) 'K1013C',K1013C
      write (6,*) 'K1013T',K1013T
      write (6,*) 'KONDC',KONDC
      write (6,*) 'KONDT',KONDT
      write (6,*) 'subroutine CBOD'
      write (6,*) 'OCRB',OCRB
      write (6,*) 'KDC',KDC
      write (6,*) 'KDT',KDT 
      write (6,*) 'KDSC',KDSC
      write (6,*) 'KDST',KDST
      write (6,*) 'KBOD',KBOD
      write (6,*) ' subroutine dissoxyg'
      write (6,*) 'SOD1D(1)',SOD1D(1)
      write (6,*) 'SODTA(1)',SODTA(1)
      write (6,*) ' light limitation'
      write (6,*) 'LGHTSW',LGHTSW
      write (6,*) 'IS2',IS2
      write (6,*) 'PHIMX',PHIMX 
      write (6,*) 'XKC',XKC
      write (6,*) 'ITOT',ITOT
      write (6,*) 'iav',iav
      write (6,*) 'DTDAY',DTDAY

      end

c********************************************************************

      subroutine source(id,t,vel0,uws0,stp0,sal0,vol0,depth0,qss0,c,cds)

c vel      velocity
c stp      temperature
      
      implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )

      integer id            !id of box
      real t                  ![day]
      real vel0            ![m/s]
      real uws0		   !wind velocity [m/s]
      real stp0            ![Celsius]
      real sal0            ![g/L], e.g., 30.
      real qss0            !solar radiation from shyfem
      real vol0            ![m**3]
      real depth0            ![m]
      real c(nstate)            ![mg/L]            !FIXME
      real cds(nstate)      ![g/day]      !FIXME

c if rlux is > 0, at return it will be replaced with the attenuation
c factor at the bottom of the box -> therefore it can be used
c with the next call in 3D models looping from surface to bottom
c for 2D models the value of rlux should be 1 (value replaced at return)
c or <= 0 which leaves the value untouched
c for the surface box the value of rlux should be 1

        INCLUDE 'weutro.h'
      
      integer i,n

        idbox = id

      iseg = 1
      sedseg = .false.
      ito = 0
      
      daytime = t

c deal with light climate


c      call rintens(daytime,itot,fday,iinst)      !compute instant light

         NH3 = C (1)
         NO3 = C (2)
         OPO4 = C (3)
         PHYT = C (4)
         CBOD = C (5)
         DO = C (6)
         ON = C (7)
         OP = C (8)
       ZOO = C (9)

      stp = stp0
      stp20 = stp-20.
      sal = sal0
      vel = vel0
      wind = uws0
      h = depth0
      vol = vol0

        SA = vol0/depth0

      iinst=qss0
c      write(6,*)iinst

C                        Compute derivatives
C              For PHYT, OP, OPO4, ON, NH3, NO3, CBOD, DO

      do i=1,nstate
        cd(i,iseg) = 0.
      end do

       CALL ZOOP
         CALL PHYTO
         CALL ORGANOP
         CALL INORGANP
         CALL ORGANICN
         CALL AMMONIA
         CALL NITRATE
         CALL CBODSV
         CALL DISSOXYG

         if( botseg ) CALL BENTFLUX

c adesso sono definiti i cd

      do i=1,nstate
        cds(i) = cd(i,iseg)
      end do

c save light climate for next call


      end

c********************************************************************
c********************************************************************
c********************************************************************

      subroutine wmeteo(tempair,windspeed)

c sets meteo parameters
c
c tempair is the air temperature in degrees C
c windspeed is the wind speed in m/s

      implicit none

      include 'weutro.h'

      real tempair
      real windspeed

      airtmp = tempair
      wind   = windspeed

      end

c********************************************************************

      subroutine wturbid(turbid)

c sets turbidity parameter (units [1/m])

      implicit none

      include 'weutro.h'

      real turbid

      !kefn(1) = 1
      !ke(1)   = 1.

      kesg(1) = turbid

      end

c********************************************************************

      subroutine wlight(fracday,itotal)

c sets light parameters
c
c fracday      the fraction of the day that is light [0-1]
c itotal      average incident light intensity during whole day [ly/day]

      implicit none

      include 'weutro.h'

      real fracday
      real itotal

      fday = fracday
      itot = itotal
      itotmp = itotal            !not used

      iav=0.9 *itot/FDAY
c       IAV = 0.9*ITOTmp*RLGHTS (ISEG, 2)/FDAY !is multiplied in ditoro

      end

c********************************************************************

      subroutine icecover(ice)

c sets ice cover parameter
c ice=1: ice covers all the surface; ice=0 no ice cover      

      implicit none

      include 'weutro.h'

      real ice

      xicecvr = 1. - ice

      end

c********************************************************************

      subroutine grazing(zoopl)

c sets zooplankton state variable
c
c set in zoo the C/l value of grazing zooplankton
c may be used to implement measured zoo

      implicit none

      include 'weutro.h'

      real zoopl

      zoo = zoopl

      end

c********************************************************************

      subroutine sedflux(flnh4,flpo4)

c sets sediment fluxes
c fluxes are in [mg/(m**2 day)]

      implicit none

      include 'weutro.h'

      real flnh4, flpo4

      fnh4(1) = flnh4
      fpo4(1) = flpo4

c      write(6,*) 'sedflux : ',flnh4,flpo4

      end

c********************************************************************
c********************************************************************
c********************************************************************

        subroutine settopseg(bstate)

c sets segment as surface or not

        implicit none
        include 'weutro.h'
        logical bstate

        topseg = bstate

        end

c********************************************************************

        subroutine setbotseg(bstate)

c sets segment as bottom or not

        implicit none
        include 'weutro.h'
        logical bstate

        botseg = bstate

        end

c********************************************************************

        subroutine setsedseg(bstate)

c sets segment as sediment or not

        implicit none
        include 'weutro.h'
        logical bstate

        sedseg = bstate

        end

c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************

      SUBROUTINE EUTROINT

C     Initialize Values
C
C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:53.

      implicit none

CDDD      changed PHIMAX ---> PHIMX

      INCLUDE 'weutro.h'
      real R0MIN
      parameter( R0MIN = 1.e-15 )


      real xarg,phimax
      integer i,j,initb,mxdmp,mxseg,noseg

        IDBOX = 0

      topseg = .true.
      botseg = .true.
        wdebug = .false.

c      DUMMY = 0.
c      CHLA2 = 0.
      GP1 = 0.
      DP1 = 0.
c      GP2 = 0.
c      DP2 = 0.
c      GZ1 = 0.
c      DZ1 = 0.
c      GZ2 = 0.
c      DZ2 = 0.
c      CFOREA = 1.0
C
      noseg=segmax

      DO 1000 J = 0,noseg 
         RLGHTS (J, 1) = 0.0
         RLGHTS (J, 2) = 1.0
 1000 CONTINUE
C
      INITB = 1
      MXDMP = 4
C
C             Check to see if Michalis Menton constants are Zero
C             and readjust values to prevent floating zero divide
C
CRBA--Date: Tuesday, 1 June 1993.  Time: 09:01:20.
      XARG = ABS(NCRB)
      IF (XARG .LT. R0MIN) NCRB = 0.25
      XARG = ABS(PCRB)
      IF (XARG .LT. R0MIN) PCRB = 0.025
c      XARG = ABS(LGHTSW)
c      IF (XARG .LT. R0MIN) LGHTSW = 1.0
      if( LGHTSW .lt. 0 ) LGHTSW = 0
      if( LGHTSW .gt. 3 ) LGHTSW = 3
      
CCSC
      XARG = ABS(KBOD)
      IF (XARG .LT. R0MIN) KBOD = 1.00E-20
CCSC
      XARG = ABS(KNO3)
      IF (XARG .LT. R0MIN) KNO3 = 1.00E-20
CCSC
      XARG = ABS(KPO4)
      IF (XARG .LT. R0MIN) KPO4 = 1.00E-20
CCSC
      XARG = ABS(KNIT)
      IF (XARG .LT. R0MIN) KNIT = 1.00E-20
CCSC
      XARG = ABS(KMNG1)
      IF (XARG .LT. R0MIN) KMNG1 = 1.00E-20
CCSC
      XARG = ABS(KMPG1)
      IF (XARG .LT. R0MIN) KMPG1 = 1.00E-20
CCSC
      XARG = ABS(KMPHYT)
      IF (XARG .LT. R0MIN) KMPHYT = 1.00E-20
C
CCSC
      XARG = ABS(OCRB)
      IF (XARG .LT. R0MIN) OCRB = 32./12.
CCSC
      XARG = ABS(IS1)
      IF (XARG .LT. R0MIN) IS1 = 300.
CCSC
      XARG = ABS(CCHL)
      IF (XARG .LT. R0MIN) CCHL = 30.
CCSC
      XARG = ABS(FON)
      IF (XARG .LT. R0MIN) FON = 1.0
CCSC
      XARG = ABS(FOP)
      IF (XARG .LT. R0MIN) FOP = 1.0
CCSC
      XARG = ABS(PHIMX)
      IF (XARG .LT. R0MIN) PHIMX = 720.
CCSC
      XARG = ABS(XKC)
      IF (XARG .LT. R0MIN) XKC = 0.017
C
C  Check for Zero Temperature Correction Factors and readjust to 1.0
C
CCSC
      XARG = ABS(K1320T)
      IF (XARG .LT. R0MIN) K1320T = 1.0
CCSC
      XARG = ABS(K140T)
      IF (XARG .LT. R0MIN) K140T = 1.0
CCSC
      XARG = ABS(K1T)
      IF (XARG .LT. R0MIN) K1T = 1.0
CCSC
      XARG = ABS(K1RT)
      IF (XARG .LT. R0MIN) K1RT = 1.0
CCSC
      XARG = ABS(KDT)
      IF (XARG .LT. R0MIN) KDT = 1.0
CCSC
      XARG = ABS(K1013T)
      IF (XARG .LT. R0MIN) K1013T = 1.0
CCSC
      XARG = ABS(KONDT)
      IF (XARG .LT. R0MIN) KONDT = 1.0
CCSC
      XARG = ABS(K58T)
      IF (XARG .LT. R0MIN) K58T = 1.0
CCSC
      XARG = ABS(KOPDT)
      IF (XARG .LT. R0MIN) KOPDT = 1.0
CCSC
      XARG = ABS(KPZDT)
      IF (XARG .LT. R0MIN) KPZDT = 1.0
CCSC
      XARG = ABS(KDST)
      IF (XARG .LT. R0MIN) KDST = 1.0

c sediment fluxes

      do i=1,segmax
          FNH4 (i) = 0.
          FPO4 (i) = 0.
      end do

c      write(6,*) segmax,FNH4(1),FPO4(1)

      RETURN
      END

c********************************************************************

      SUBROUTINE AMMONIA

C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:55.

      implicit none

      real SR13ON,SR13P,SK13P1

      INCLUDE 'weutro.h'
      include 'donata.h'      !GGU (and others)
C
C       *-*-*-*-*  SYSTEM 1 - AMMONIA (NH3-N)  *-*-*-*-*
C
C                        Sources
C               Mineralization of organic nitrogen
C
      SR13ON = SK1013
      denit = 0.
      denit = SK14D
C
C                  Phytoplankton Death
C
      SR13P = NCRB*DPP*(1.0 - FON)
C
C                        Sinks
C                    Algal Uptake
C
      SK13P1 = PNH3G1*NCRB*GP1*PHYT
C
C                      Nitrification
C
      IF (DO .GT. 1.0E-10) THEN
         SK1314 = (K1320C*K1320T**STP20)*NH3*DO/(KNIT + DO)
      ELSE
         SK1314 = 0.0
      END IF
      IF (STP .LT. 7.) SK1314 = 0.0
       ddin1=sk1013+SR13P-SK13P1        !GGU new
c      write(6,*)ddin1

C
C                   Formulate Derivative
C
c      write(6,*) 'ammonia debug :'
c      write(6,*) SR13P,SR13ON,SK13P1,SK1314,VOL,CD (1, ISEG)
c      write(6,*) PNH3G1,NCRB,GP1,PHYT

      if (wdebug) then
      write(99,*)'ammondebug:',SR13P,denit,SR13ON,SK13P1,SK1314
      end if

c      CD (1, ISEG) = (SR13P +denit+ SR13ON - SK13P1 - SK1314)*VOL      !just for tests
      CD (1, ISEG) = (SR13P + SR13ON - SK13P1 - SK1314)*VOL     !ORIG

c      write(6,*) CD (1, ISEG),iseg
C
      RETURN
      END

c********************************************************************

      SUBROUTINE NITRATE

C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:56.

      implicit none
      INCLUDE 'weutro.h'
      include 'donata.h'
      real SR1413,SK14P1
C
C    *-*-*-*-* SYSTEM 2 - Nitrate (+Nitrite)  (NO3-N+NO2-N)  *-*-*-*-*
C
C                         Sources
C                      Nitrification
C
      SR1413 = SK1314
C
C                          Sinks
C                      Algal Uptake
C
      SK14P1 = (1. - PNH3G1)*NCRB*GP1*PHYT
C
C                     Denitrification
C
      SK14D = (K140C*K140T**STP20)*NO3
      IF (DO .GT. 0) SK14D = SK14D*KNO3/(KNO3 + DO)
C
c for no denitrification set SK14D to 0.

c      IF (SK14D .LT. 1.00E-24) SK14D = 1.00E-24

      ddin2=-SK14P1-SK14D       !GGU new
c      write (6,*) ddin2


C
C                   Formulate Derivative
C
      if (wdebug) then
        write(99,*)'NO2debug:',SR1413,SK14P1,SK14D,DO,NO3
      end if


      CD (2, ISEG) = (SR1413 - SK14P1 - SK14D)*VOL
C
      RETURN
      END

c********************************************************************

      SUBROUTINE INORGANP

C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:56.

      implicit none
      INCLUDE 'weutro.h'
      include 'donata.h'
      real SR8P,SR8OP,SK8P
C
C   *-*-*-*-*  SYSTEM 3 - Dissolved Inorganic Phosphorous  *-*-*-*-*
C
C                            Sources
C               Mineralization of Organic Phosphorous
C
      SR8OP = SK58
C
C                         Phytoplankton Death
C
      SR8P = PCRB*DPP*(1. - FOP)
C
C                            Sinks
C                        Algal uptake
C
      SK8P = PCRB*GPP

        prod = SR8OP+SR8P       !GGU new
      cons= SK8P              !GGU new
c      write(6,*) 'prod phosph',prod
C
C                   Formulate Derivative
C
      if (wdebug) then
      write(99,*)'inPdebug:',SR8P,SR8OP,SK8P
      end if

      CD (3, ISEG) = (SR8P + SR8OP - SK8P)*VOL
C
      RETURN
      END

c********************************************************************

      SUBROUTINE PHYTO

C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:57.

      implicit none

      INCLUDE 'weutro.h'
      real R0MIN
      parameter( R0MIN = 1.e-15 )
      REAL XARG,XDIFF,GIT1,CN,DOPO4
      real ttr, ttr1
C
C        *-*-*-*-*-*  System 4 - Phytoplankton  *-*-*-*-*-*
C
c         open (9, file = 'outgp.dat')      !ggu


c      write(6,*) 'phyto, sedseg,iseg : ',SEDSEG,iseg
      IF (SEDSEG) THEN

c      write(6,*) 'phyto,  ctrl: ',K1C,PHYT,K1T,STP20
c      write(6,*) 'phyto, ctrl: ',GITMX1,KEFN(1)
            GP1 = K1C
            GPP = GP1*PHYT
            PNH3G1 = 0.

         ELSE
c      write(6,*) 'phyto, ctrl : ',K1C,PHYT,K1T,STP20
            GITMX1 = K1C*K1T**STP20
            GITMAX = GITMX1
            IKE = KEFN (ISEG)
      
c         write(6,*) 'phyto, ctrl1 : ',GITMAX
C
C               Compute growth rate reduction due to light conditions
C               using either Dick Smith's or Di Toro's formulation
C
            RLIGHT = 1.
            if( LGHTSW .eq. 1 ) call ditoro
            if( LGHTSW .eq. 2 ) call smith
            if( LGHTSW .eq. 3 ) call steele

             GIT1 = RLIGHT*GITMAX
c            write(6,*) 'phyto, ctrl2 : ',GITMAX,GIT1

C                   Compute ammonia preference
            PNH3G1 = 0.0
C
C
            IF (NH3 .GE. 1.0E-5)
     1         PNH3G1 = NH3*NO3/((KMNG1 + NH3)*(KMNG1 + NO3))
     2         + NH3*KMNG1/((NH3 + NO3)*(KMNG1 + NO3))
C
C              Compute Michaelis Limitations
C
            CN = NH3 + NO3
c      write(6,*)KMNG1,PNH3G1,NH3,NO3,CN
            XEMP1 = CN/(KMNG1 + CN)

               DOPO4 = OPO4*FOPO4
            XEMP2 = DOPO4/(KMPG1 + DOPO4)
      
C
C       Compute Growth Rate Reduction due to Nutrient Limitation
C
CCSC
            XARG = ABS(NUTLIM)
            IF (XARG .LT. R0MIN) RNUTR = AMIN1 (XEMP1, XEMP2)

c      write(6,*) '-------------- ggu --------------'
c      write(6,*) NUTLIM,XARG,R0MIN
c      write(6,*) XEMP1,XEMP2,RNUTR
c      write(6,*) '-------------- ggu --------------'
CCSC
            IF (NUTLIM .EQ. 1.) RNUTR = XEMP1*XEMP2
c            write(6,*) RNUTR
            GP1 = RNUTR*GIT1
            GPP = GP1*PHYT
      
c
C       ********************************************
C                     Respiration Rate
C       ********************************************
C
         RESP = K1RC*K1RT**STP20
c      
C
C       ALGAL RESPIRATION + DEATH + GRAZING
C
c         DP1 = RESP + K1D + K1G*ZOO*ZOOSG(ISEG)      !old

         DP1 = RESP + K1D             !FIXMED 

         DPP = DP1*PHYT      !BUGFIX1   LAA
c         DPP = DPP + DP1*PHYT       
c         RESP = RESP*PHYT      !BUGFIX1

c      write(6,*) DP1,RESP,K1D,K1G,ZOO,ZOOSG(iseg)


C         IF (PHTY .GT. 1.0E-6)THEN
      XEMPRC=1.
Cdoni            XEMPRC = PHYT/(KMPHYT + PHYT)
C         ELSE
C            XEMPRC = 1.0E-6/(KMPHYT + 1.0E-6)
C         ENDIF

C
      END IF      !sedsed
C
c      write(6,*) 'phyto debug ggu'
      if (wdebug) then
      write(99,*)'phytodebug:',GP1,DPP,GRZ
      end if

c        write(9,'(3(f8.4,2x))') RLIGHT,DP1, RESP       !ggu

c      CD (4, ISEG) = (GP1 - DP1-GRZ)*PHYT*VOL      !BUGFIX1
      CD (4, ISEG) = (GPP - DPP - GRZ)*VOL
      RETURN
      END

c********************************************************************

      subroutine zoop

c zooplankton

      implicit none
      include 'weutro.h'
      real gra

c      *-*-*-*- system 9 zooplankton *-*-*-
c
c      Sources: zooplankton growth

c zoosk      Sink term: zooplankton death -->source term for organop organicn cbodsv
c          grazing inefficiency -> source term for  organop organicn cbodsv

      if( graztype.eq.1 ) then
        DPP = 0.
          GRZ=KGRZ*ZOO*PHYT/(PHYT+KPHYZ)
        zoosk = (1-EFF)*GRZ + KDZ*ZOO
        CD (9, ISEG) = (EFF*GRZ - KDZ*ZOO)*VOL
      else      !BUGFIX1
        DPP = K1G*ZOO*PHYT*ZOOSG(ISEG)      !original formulation
        GRZ = 0.
        zoosk = 0. 
        CD (9, ISEG) = 0.
      end if
      !zood=KDZ*ZOO           !GGU ???

      if (wdebug) then
      write(99,*)'zoodebug:',GRZ
      end if

      return
      end

c********************************************************************

      SUBROUTINE CBODSV

C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:57.

      implicit none
      INCLUDE 'weutro.h'
      real SR18P,SK18D,GRC
C
C                *-*-*-*-* SYSTEM 5 - CBOD  *-*-*-*-*
C                             Sources
C
      IF ( .NOT. SEDSEG) GO TO 1000
      SR18P = OCRB*DPP
      GO TO 1010
C
C                        Phytoplankton 'DEATH'
C
 1000 CONTINUE
      DEATH = K1D*PHYT
      SR18P = OCRB*DEATH
      GRC= ZOOSK * OCRB
C
C                              Sinks
C                            Oxidation
C
 1010 CONTINUE
      IF ( .NOT. SEDSEG) THEN
c         IF (DO .GT. 1.0E-15) THEN
         IF (DO .GT. 1.0E-15 .AND. CBOD .GT. 1.0E-5) THEN
            SK180 = (KDC*KDT**STP20)*((CBOD*DO)/(KBOD + DO))
         ELSE
            SK180 = 0.0
         END IF
      ELSE
         IF (CBOD .GT. 1.0E-15) THEN
            SK180 = (KDSC*KDST**STP20)*CBOD
         ELSE
            SK180 = 0.0
         END IF
      END IF
C
C                         Denitrification
c
c ggu & mcg: next formula is indepenedent of BOD -> insert dependence

       SK18D = (5./4.)*(32./14.)*SK14D
C
C                      Formulate Derivative
      if (wdebug) then
      write(99,*)'cboddebug:',sr18p,grc,sk180,sk18d
      end if
C
      CD (5, ISEG) = (SR18P+GRC - SK180 - SK18D)*VOL
C
      RETURN
      END

c********************************************************************

      SUBROUTINE DISSOXYG

C     Last Revised:  Date: Monday, 26 August 1991.  Time: 10:37:53.

      implicit none
      INCLUDE 'weutro.h'
      real R0MIN
      parameter( R0MIN = 1.e-15 )
      REAL XARG1, XARG2
      REAL K2WIND, K2HYDRA
      real SR190,SR19PA,SR19PB,SK1913,SK1918,SK19S
      real TK,RLNCS,CL
C
C              *-*-*-*-*  SYSTEM 6 - OXYGEN  *-*-*-*-*
C
C                          Sources
C
      K20 = 0.0
      IF ( .NOT. SEDSEG) GO TO 1000
      SR190 = 0.0
      SR19PA = 0.0
      SR19PB = 0.0
      SK19P = 0.0
      GO TO 1010
C
C                         Reaeration
C
 1000 CONTINUE
      IF ( topseg .AND. XICECVR .GT. 0.0) THEN

c if surface element and not completly covered with ice

c elimino rear e rearsg perch\E9 non ci interessano, sono una time function 
c e una segment specific reareation rate, ma noi usiamo la reareation from
c  kawind o kahydra o imponendo un valore della costante di reareazione

         IF (K2 .EQ. 0.0) THEN
            CALL KAWIND (WIND, STP, AIRTMP, H, WTYPE,  K2WIND) 
            CALL KAHYDRA (K2HYDRA)
            IF (K2WIND .GT. K2HYDRA) THEN
               KA = K2WIND * XICECVR
            ELSE
               KA = K2HYDRA * XICECVR
            END IF
         ELSE
            IF (K2 .GT. 0)THEN
                KA = ((K2*1.028**STP20)* XICECVR)
            ELSE
            KA = 0.
            ENDIF
         ENDIF
      ELSE
         KA = 0.0
      END IF

C       Calculate oxygen saturation level for current water
C       temperature; DOSAT is expressed as mg oxygen per liter

      CL = SAL/1.80655
      TK = STP + 273.
      RLNCS = - 139.34411 + (1.575701E05/TK) - (6.642308E07/TK**2) +
     1   (1.243800E10/TK**3) - (8.621949E11/TK**4) -
     2   (CL*(3.1929E-02 - (19.428/TK) + (3.8673E03/TK**2)))
C
      CS = EXP (RLNCS)
C
      SR190 = KA*(CS - DO)

c      write(88,*) KA,CS,SR190
C
C                 Evolution by phytoplankton
C          growth of phytoplankton using CO2 and NH3
C
      SR19PA = PNH3G1*GP1*PHYT*32./12.
C
C       Growth of phytoplankton using CO2 and NO3 (2NO3 = 2NH3 + 302)
C
      SR19PB = (1. - PNH3G1)*GP1*PHYT*32.*(1./12. + 1.5*NCRB/14.)

C
      SR19P = SR19PA + SR19PB
C
C      NOTE: SR19P = GPP*(32/12 + (1.5*NCRB/14)*(1-PNH3G1))
C
C                             Sinks
C                       Algal Respiration

c      SK19P = OCRB*RESP                  !BUGFIX1
      SK19P = OCRB*RESP*PHYT
C
C        Nitrification (NH3-N + 2O2 = NO3-N + H2O + H)
C
 1010 CONTINUE
      SK1913 = 64./14.*SK1314
C
C                      Oxidation of CBOD
C
      SK1918 = SK180
C
C             Sediment Oxygen Demand (1-D Networks)
C
      SK19S = SOD1D (ISEG)*SODTA(ISEG)**stp20/H
C
C======================================================================
C                     Formulate Derivative
C
      if (wdebug) then
      write(99,*)'doxdebug:',SR190,SR19PA,SR19PB,SK19P
     +                        ,SK1913,SK1918,SK19S
      end if

      CD (6, ISEG) = (SR190 + SR19PA + SR19PB - SK19P
     1   - SK1913 - SK1918 - SK19S)*VOL

c        write(54,'(8(f8.4,2x))') SR190,SR19PA,SR19PB,SK19P,SK1913,SK1918
c     & ,SK19S 
c        write(54,'(3(f8.4,2x))') SR19PA,SR19PB,GP1      !ggu
       
      RETURN
      END

c********************************************************************

      SUBROUTINE ORGANICN

C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:59.

      implicit none
      INCLUDE 'weutro.h'
      real SR10P,GRN
C
C
C             *-*-*-*-*  SYSTEM 7 Organic Nitrogen      *-*-*-*-*
C
C                          Sources
C             Phytoplankton Respiration and 'DEATH'
C            GRAZING and zooplankton death
C                  GRN=  ZOO --> ORGANICN

C
      SR10P = NCRB*DPP*FON
      GRN   = ZOOSK*NCRB      !BUGFIX1
C
C                         Sinks
C         Mineralization of Dissolved Organic Nitrogen
C
      IF ( .NOT. SEDSEG) SK1013 = (K1013C*K1013T**STP20)*ON*XEMPRC
      IF (SEDSEG) SK1013 = (KONDC*KONDT**STP20)*ON
C
C                   Formulate Derivative
C
      if (wdebug) then
      write(99,*)'orgNdebug:',SR10P,GRN,SK1013
      end if
      CD (7, ISEG) = (SR10P+ GRN - SK1013)*VOL
C
      RETURN
      END

c********************************************************************

      SUBROUTINE ORGANOP

C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:59.

      implicit none
      real SR5P,GRP

      INCLUDE 'weutro.h'
C
C        *-*-*-*-*-*  SYSTEM 8 Organic Phosphorous  *-*-*-*-*-*
C
C                            Sources
C                 Phytoplankton respiration and 'death'
C               GRAZING and zooplankton death
C                  GRP=  ZOO --> ORGANOP
C
      SR5P = PCRB*DPP*FOP
      GRP  = ZOOSK * PCRB      !BUGFIX1
C
C                             Sinks
C
C           Mineralization of dissolved organic phosphorus and
C           Phytoplankton respiration and 'death'
C
      IF ( .NOT. SEDSEG) SK58 = (K58C*K58T**STP20)*OP*XEMPRC
      IF (SEDSEG) SK58 = (KOPDC*KOPDT**STP20)*OP
C
C                        Formulate Derivative
C
      if (wdebug) then
      write(99,*)'orgPdebug:',SR5P,GRP,SK58
      end if
      CD (8, ISEG) = (SR5P + GRP - SK58)*VOL
C
      RETURN
      END

c********************************************************************
c
c IAVBOT                bottom light                            not used
c IAVBOTX(ISEG)         bottom light (variable in seg)          not used
c daytime               = t
c KESG(ISEG)            given turbidity [1/m]
c KE(IKE)               time function for turbidity
c SKE                   turbidity (given and phyto) [1/m]
c
c RLIGHT                limiting factor
c RLGHTS(ISEG, 1)       = RLIGHT
c RLGHTS(ISEG, 2)       fraction of surface light in segment (0-1)
c                          is 1 in surface layer
c ITO                   pointer to lower segment
c IS1,IS1X(ISEG)        saturation light intensity
c IS2                   optimum light intensity for Steele
c FDAY                  fraction of day that is daylight (unitless)
c IAV                   average incident light intensity during daylight hours
c                          (iav=0.9*itot/fday)
c ITOT                  total incident light intensity during one day [ly/day]
c ITOTMP                ...                                     not used
c
c smith and steele use ITOT
c ditoro uses IAV
c
c summary: we need IAV, ITOT, IINST, FDAY, IS1, IS2, KESG
c
c primary:   ITOT, FDAY
c derived:   IAV, IINST
c constants: IS1, IS2, KESG
c
c set primary, compute derived, use constants
c
c moreover: 
c      RLGHTS(ISEG,2) is attenuation coefficient (enter) [0-1]
c      RLGHTS(ITO,2) is attenuation coefficient at bottom (return) [0-1]
c
c********************************************************************

      SUBROUTINE DITORO

C*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
C        Di Toro et al Light Formulation
C*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
c
c needs IAV, IS1 and FDAY (plus KESG)

      implicit none
      INCLUDE 'weutro.h'
      real CCHL1,TCHLA,KESHD,SKE
      real temp1, temp2, temp3

      CCHL1 = CCHL
      TCHLA = PHYT/CCHL1
      KESHD = (0.0088*1000.*TCHLA + 0.054*(1000.*TCHLA)**0.6667)
      SKE = KESG (ISEG)
      IF (IKE .GT. 0 .AND. IKE .LE. 5) SKE = SKE*KE (IKE)
      SKE = SKE + KESHD
      TEMP1 = SKE*H

C         Get average solar radiation during daylight hours

      TEMP2 = IAV/IS1
      TEMP2 = RLGHTS(ISEG,2)*IAV/IS1      !consider light attenuation - $LGA
      TEMP3 = EXP ( - TEMP1)
      IAVBOT=IAV * TEMP3
      RLGHTS (ITO, 2) = RLGHTS (ISEG, 2)*TEMP3
      RLIGHT = 2.718*FDAY/TEMP1*(EXP ( - TEMP2*
     1   TEMP3) - EXP ( - TEMP2))
      RLGHTS (ISEG, 1) = RLIGHT

c      write(6,*) 'ditoro debug ggu'
c      write(6,*) H,SKE,TEMP1,TEMP2,TEMP3,IAVBOT
c      write(6,*) FDAY,RLGHTS (ISEG, 1),RLIGHT
c      write(6,*) ISEG,ITO,IS1,IAV, 2.718*FDAY,temp1
c      write(6,*) H
c       write(10,'(6(f8.5,2x))')KESHD,SKE,PHYT,RLIGHT

      RETURN
      END

c********************************************************************

      SUBROUTINE SMITH

C*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
C          Dick Smith variable carbon/chlorophyll ratio
C*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
c
c needs ITOT and FDAY (plus KESG)
c IS1 is not needed (optimal light value is computed)
c RLGHTS (ISEG, 2) is always 1 as implemented now

      implicit none
      INCLUDE 'weutro.h'

      real I0,IMAX,IAVSG,SUM,KESHD,SKE
c      real IAVBOTX(segmax)
      real CCHL1,TCHLA
       real temp1, temp2, temp3
      integer I
C
C            IAV = IAV/1440.
C
C               (Average solar radiation during daylight hours)
C
C/cm2 day*[0.43 visible/total*10,000 cm2/m2*E/52,000]=E visible/m2-day
C             0.083  is FU (E/M2 - Ly or E/10 Kcal)
C
CRBA--Date: Tuesday, 1 June 1993.  Time: 09:05:36.
      CCHL1 = CCHLX(ISEG)
      IS1 = IS1X(ISEG)
      IAVBOT=IAVBOTX(ISEG)
C      CCHL1 = 0.3*0.083*PHIMX*XKC*IAV/(GITMX1*2.718)
C      IF (CCHL1 .LT. 20.) CCHL1 = 20.
      TCHLA = PHYT/CCHL1
      RLIGHT = RLGHTS (ISEG, 1)
C
      IF (NEWDAY .GE. 1) THEN
C
C           Dick Smith formulation integrated every day
C
         KESHD = XKC*1000.*TCHLA      !FIXME
         SKE = KESG (ISEG)
         IF (IKE .GT. 0 .AND. IKE .LE. 5) SKE = SKE*KE (IKE)
         SKE = SKE + KESHD
         TEMP1 = SKE*H
         TEMP2 = 0.083*PHIMX*XKC/(GITMAX*CCHL1*2.718)
         TEMP3 = EXP ( - TEMP1)
         RLGHTS (ITO, 2) = RLGHTS (ISEG, 2)*TEMP3
         IMAX = PI*ITOT*RLGHTS (ISEG, 2)/(2.*FDAY)
         SUM = 0.0
         DO 1000 I = 1, 25
            DTDAY = FLOAT (I - 1)/24.
            IF (DTDAY .GT. FDAY) GO TO 1010
            I0 = IMAX*SIN (PI*DTDAY/FDAY)
            SUM = SUM + 2.7183/TEMP1*
     1         (EXP ( - TEMP2*I0*TEMP3) - EXP ( - TEMP2*I0))
 1000    CONTINUE
 1010    CONTINUE
         RLIGHT = SUM/24.
         RLGHTS (ISEG, 1) = RLIGHT
CRBA--Date: Tuesday, 1 June 1993.  Time: 09:06:22.
C        Adapt carbon to chlorophyll ratio:
         IAVSG=IAV*(1.0-TEMP3)/TEMP1
         CCHLX(ISEG)=0.3*0.083*PHIMX*XKC*IAVSG/(GITMX1*2.718)
         IF(CCHLX(ISEG).LT.20.0) CCHLX(ISEG)=20.0
         IS1X(ISEG)=1/TEMP2
         IAVBOTX(ISEG)=IAV*TEMP3
      END IF
      RETURN
      END

c***************************************************************

      subroutine steele

c computes light limitation with steele formula
c
c uses instantaneous light intensity IINST, and not IAV or ITOT
c
c needs IINST and IS2 (plus KESG)

      implicit none

      include 'weutro.h'

      real expon
      parameter ( expon = 2.718 )

      real TCHLA,KESHD,SKE
        real temp1,temp2,temp3

        TCHLA = PHYT/CCHL
        KESHD = (0.0088*1000.*TCHLA + 0.054*(1000.*TCHLA)**0.6667)

        SKE = KESG (ISEG)
        IF (IKE .GT. 0 .AND. IKE .LE. 5) SKE = SKE*KE (IKE)
        SKE = SKE + KESHD

c       write(66,*) daytime,itot,fday,iinst,ske

      temp1 = ske * h
      temp1 = ske                                    !$WRONG_FORMULA
      temp2 = iinst / is2
      temp2 = RLGHTS(ISEG,2) * iinst / is2      !light attenuation $LGA
      temp3 = exp ( - temp1 )

        RLGHTS (ITO, 2) = RLGHTS (ISEG, 2)*TEMP3

      RLIGHT = (expon/temp1) * ( exp(-temp2*temp3) - exp(-temp2) )
        RLIGHT = temp2*temp3*exp(1-temp2*temp3)         !2D_FORMULA

        RLGHTS (ISEG, 1) = RLIGHT
c        write(66,*) ,'2', temp1,temp2,temp3,'fine2'
c        write(66,*) ,'3', RLGHTS,'fine3'
c        write(66,*) ,'4', RLIGHT,'fine4'

      end

c***************************************************************

      subroutine rintens(t,itot,fday,iinst)

c computes light intensity during a day given itot and fday
c
c uses formula 5.7 in WASP user manual
c
c needs ITOT (average incident light intensity over one whole day)
c and FDAY (fraction of day length)

      implicit none

      real t            !day [0-365]
      real itot      !average incident light intensity over one whole day
      real fday      !fraction of day length [0-1]
        real iinst      !instantaneous light intensity at time t (return)
      
      real pi
      parameter( pi = 3.14159 )

      real tday,aux

      real t_old,iinst_old
      save t_old,iinst_old
      data t_old,iinst_old /0.,0./      !hope we are not in Antarktic

      if( t .eq. t_old ) then            !if for same daytime we are done
        iinst = iinst_old
        return
      else
        t_old = t
      end if

c      write(6,*) ,'rintens1',iinst

      tday = t - int(t)                  !fraction in day
      if( tday .lt. 0. ) tday = tday + 1.      !negative days
      tday = tday - 0.5                  !maximum at noon

c        intens = itot * (300./86400.)         !LIGHTFIX !FIXME

      aux = pi / fday

      if( abs(tday) .le. fday/2 ) then
        iinst = 0.5 * itot * aux * cos( tday * aux )
      else
        iinst = 0.
      end if

      iinst_old = iinst

c        write(66,*) ,'intens' ,t, tday,itot,iinst
      end

c********************************************************************

      SUBROUTINE KAWIND (WS, TW, TA, depth, WTYPE,  RK)

C*     THIS SUBROUTINE CALCULATES:
C*
C*              RK = REAERATION COEFFICIENT (RK) (M/DAY)*
C*
C*     GIVEN:
C*              WS = Wind Speed (WS) (m/s)
C*              TA = Temperature of the Air  (Degrees C)
C*              TW = Water Temperature (Degrees C)
C*
C* Using the method presented in:
C*
C*           Jour. of Env Eng, Vol. 109, NO.3,PP.731-752,
C*           June 1983, Author: D.J.O'Connor, TITLE: "Wind Effects on
C*           Gas- Liquid Transfer Coefficients"
C*
C*====================================================================
C*
C* THIS SUBROUTINE WAS WRITTEN BY:
C*     SANDRA BIRD
C*     USAE WATERWAYS EXPERIMENT STATION (WES-ES-Q)
C*     VICKSBURG, MISSISSIPPI
C* AND MODIFIED BY JAMES L. MARTIN
C*
C*
C*====================================================================
C*
C*   Parameters used in the model include:
C*
C*        Transitional Shear Velocity - UT(cm/sec)
C*        Critical Shear Velocity - UC (cm/sec)
C*        Vonkarman's Constant (KARMAN)
C*        Equilibrium Roughness - ZE (cm)
C*        1/LAM Is a Reynold's Number
C*        GAM is a a Nondimensional Coefficient Dependent on
C*        Water Body Size (WTYPE).
C*        LAM, GAM, UT, UC and ZE are Dependent on Water Body
C*        Size (See O'Conners Paper for Table of Values).
C*
C*       UT       UC      ZE    LAM     GAM
C*       10.0     11.    .35    3.0     5.          Large Scale
C*       10.0     11.    .25    3.0     6.5         Intermediate
C*        9.      22.    .25   10.     10.          Small Scale
C*
CC******************************************************************

      implicit none

      real R0MIN
      parameter( R0MIN = 1.e-15 )
c      INCLUDE 'weutro.h'
C*
C*  Declarations:
C*

        real WS         !wind speed, m/s
        real TW         !water temperature C
        real TA         !air temperature C
        real depth      !defined depth(segmax) in geometrical
        real WTYPE      !type of water body
        real RK         !reareation term calculated in kawind

      REAL XARG1,XARG2,XWDIFF
      Character*1 cont1
      REAL*4 KARMAN, LAM, KA3, TMPAIR
       real DIFF       !diffusivity in water, from TW CM**2/SE
        real VW         !auxiliar variable to calculate RK=rearsg
        real VA         !viscosity of air, CM**2/SEC  function of TA
        real PA         !density of air g/cm**3
        real PW         !density of water g/cm**3
        real WH         !constant in kawind: =1000
        real SRCD       !constant in kawin: =0.04
        real EF         !auxiliar variable
        real F1         !auxiliar variable
        real F2         !auxiliar variable
        real FP1        !auxiliar variable
        real FP2        !auxiliar variable
        real FP3        !auxiliar variable
        real FP4        !auxiliar variable
        integer N
        real SRCD2      !auxiliar variable
        real ERR        !auxiliar variable
        real  US        !auxiliar variable
        real Z0         !auxiliar variable
        real RK1        !auxiliar variable
        real GAMU       !!auxiliar variable
        real RK2        !!auxiliar variable
        real RK3        !!auxiliar variable

c       !!!!???!!!decidere      quale tenere depth(iseg) o depth. main1 o param?
        real ut !Transitional Shear Velocity, cm/s defined in kawind
        real uc !Critical Shear Velocity, cm/sec  defined in kawind
        real ze !Equilibrium Roughness, cm
        real gam!Nondimensional Coefficient Dependent on water body type in kaw.

      integer IWTYPE

      real*4 cddrag
      COMMON /KAHOLD/ ut, uc, ze, lam, gam, TMPAIR
      data cont1/'$'/
C*
C*   Determine Water Body Type, if WTYPE=0., then default is large
C*   Water Body:
C*
CCSC
CRBA      XARG1 = ABS(WTYPE)
cdd      XWDIFF= WTYPE - 3.0
cdd      XARG2 = ABS(XWDIFF)
      IWTYPE = NINT(WTYPE)
      IF (IWTYPE .eq. 3) then
CRBA      IF (XARG1.LT.R0MIN.OR.XARG2.LT.R0MIN) THEN
cddd      IF (XARG2.LT.R0MIN) THEN
         UT = 10.0
         UC = 11.0
         ZE = 0.35
         LAM = 3.0
         GAM = 5.0
CCSC
cddd         XWDIFF = WTYPE - 1.0
cddd         XARG1 = ABS(XWDIFF)
      ELSE IF (IWTYPE .eq. 1 ) THEN
cddd         IF (XARG1 .LT. R0MIN) THEN
            UT = 9.0
            UC = 22.0
            ZE = 0.25
            LAM = 10.0
            GAM = 10.0
      ELSE if( WTYPE .eq. 2 ) then
            UT = 10.0
            UC = 11.0
            ZE = 0.25
            LAM = 3.0
            GAM = 6.5
      else
c            write(6,*) 'iwtype: ',iwtype
            stop 'error stop'
      END IF
C*
C CALCULATE DIFFUSIVITY OF OXYGEN IN WATER (DIFF) (CM**2/SEC), VISCOSIT
C OF WATER (VW) (CM**2/SEC),VISCOSITY OF AIR (VA) (CM**2/SEC),DENSITY
C OF WATER (PW) (G/CM**3), DENSITY OF AIR (PA) (G/CM**3)
      DIFF = 4.58E-07*TW + 1.2E-05
C  NOTE: IF OTHER CHEMICALS WERE USED, THEN THEIR DIFFUSIVITIES
C  MAY VARY. FOR EXAMPLE FOR TCDD:   (JLM)
C      DIFF=4.83E-6
C
      VW = 0.0164 - .00024514*TW
      VA = .133 + .0009*TA
      PA = .00129 - .0000040*TA
      PW = 1.00
      WS = WS*100.
      RK = 1.
C USE NEWTON RAPHSON METHOD TO CALCULATE THE SQUARE ROOT OF THE DRAG
C COEFFICIENT
      N = 0
C PARAMETERS USED IN THE MODEL INCLUDE TRANSITIONAL SHEAR VELOCITY - UT(
C CRITICAL SHEAR VELOCITY - UC (CM/SEC); VONKARMAN'S CONSTANT (KARMAN);
C EQUILIBRIUM ROUGHNESS - ZE (CM); 1/LAM IS A REYNOLD'S NUMBER; GAM IS
C NONDIMENSIONAL COEFFICIENT DEPENDENT ON WATER BODY SIZE.  LAM, GAM, UT
C UC AND ZE ARE DEPENDENT ON WATER BODY SIZE
      KARMAN = 0.4
      KA3 = KARMAN**.3333
      WH = 1000.
C MAKE INITIAL GUESS FOR SQUARE ROOT OF THE DRAG COEFFICIENT
      SRCD = 0.04
 1000 CONTINUE
      N = N + 1
C CALCULATE VALUE OF FUNCTION(F2) AND DERIVATIVE OF FUNCTION(FP)
      EF = EXP ( - SRCD*WS/UT)
      F1 = LOG ((WH/ZE) + (WH*LAM/VA)*SRCD*WS*EF)
      F2 = F1 - KARMAN/SRCD
      FP1 = 1./((WH/ZE) + (LAM*WH/VA)*SRCD*WS*EF)
      FP2 = ((WH*LAM)/(VA*UT))*SRCD*WS**2*EF
      FP3 = (WH*LAM/VA)*WS*EF
      FP4 = FP1*(FP2 + FP3) + (KARMAN/(SRCD**2))
C CALCULATE A NEW GUESS FOR SQUARE ROOT OF DRAG AND COMPARE TO
C PREVIOUS GUESS AND LOOP BACK THROUGH N-R WITH NEW GUESS IF
C APPROPRIATE
      SRCD2 = SRCD - F2/FP4
      ERR = ABS (SRCD - SRCD2)
      IF (ERR .GT. 0.0005 .AND. N .LT. 8) THEN
         SRCD = SRCD2
         GO TO 1000
      END IF
      IF (ERR .GT. 0.005 .AND. N .GE. 8) GO TO 1010
      CDDRAG = SRCD**2
      US = SRCD*WS
      Z0 = 1./((1./ZE) + LAM*US*EXP ( - US/UT)/VA)
      WS = WS/100.
      IF (WS .LT. 6.0) GO TO 1020
      IF (WS .GE. 6.0 .AND. WS .LE. 20.0) GO TO 1030
      IF (WS .GT. 20.0) GO TO 1040
C CALC 1050S VALUES FOR WINDSPEEDS LESS THAN 6.0 M/SEC
 1020 CONTINUE
      RK1 = (DIFF/VW)**0.666667*SRCD*(PA/PW)**0.5
      RK = RK1*KA3*WS/GAM
      RK = RK*3600.*24.
      GO TO 1050
C CALC 1050S VALUES FOR WINDSPEED GREATER THAN 20 M/S
 1040 CONTINUE
CRBA--Date: Tuesday, 1 June 1993.  Time: 09:15:47.
C      RK = (DIFF*PA*VA*US/(0.1*PW*VW))**0.5
      RK = (DIFF*PA*VA*US/(KARMAN*ZE*PW*VW))**0.5
      RK = RK*3600.*24./100.
      GO TO 1050
C CALC 1050S VALUES FOR WINDSPEED BETWEEN 6 AND 20 M/S
 1030 CONTINUE
      GAMU = GAM*US*EXP ( - US/UC + 1.)/UC
      RK1 = (DIFF/VW)**.6667*KA3*(PA/PW)**0.5*US/GAMU
      RK2 = (DIFF*US*PA*VA/(KARMAN*Z0*PW*VW))**0.5
      RK3 = (1./RK1) + (1./RK2)
      RK = 1./RK3
      RK = RK*3600.*24./100.
      GO TO 1050
 1050 continue
      GO TO 1060
 1010 CONTINUE
      WRITE (6, 6000)
 6000 FORMAT(5X,'SOLUTION DID NOT CONVERGE')
 1060 CONTINUE
      rk = rk/depth
cddd
      RETURN
      END

c********************************************************************

      SUBROUTINE KAHYDRA (K2HYDRA) 

C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:58.

      implicit none
      INCLUDE 'weutro.h'
      REAL K2HYDRA
      real CFOREA,AVDEPE,AVVELE,REAK,EXPREV,EXPRED,TRANDP,DIF
C
C
C                      Calculate Oxygen Reaeration
C
      CFOREA = 1.0
      AVDEPE = H
      AVVELE = ABS(VEL)
C
C
C         Calculate reaeration coefficient for free-flowing reach
C         Calculate reaeration coefficient as a power function of
C         average hydraulic depth and velocity; determine exponents
C         to depth and velocity terms and assign value to REAK
C
      IF (AVDEPE .GT. 0.61) GO TO 1000
C
C          Use Owen's formulation for reaeration
C
      REAK = 5.349
      EXPREV = 0.67
      EXPRED = - 1.85
      GO TO 1010
C
 1000 CONTINUE
C
C       Calculate transition depth; transition depth determines
C       which method of calculation is used given the current
C       velocity
C
      IF (AVVELE .GE. 0.518) GO TO 1020
      TRANDP = 0.0
      GO TO 1030
 1020 CONTINUE
      TRANDP = 4.411*(AVVELE**2.9135)
 1030 CONTINUE
C
      DIF = AVDEPE - TRANDP
      IF (DIF .GT. 0.0) GO TO 1040
C
C                 Use Churchill's formulation for reaeration
C
      REAK = 5.049
      EXPREV = .969
      EXPRED = - 1.673
      GO TO 1050
C
 1040 CONTINUE
C
C                 Use O'Connor-Dobbins formulation for reaeration
C
      REAK = 3.93
      EXPREV = 0.5
      EXPRED = - 1.5
C
 1050 CONTINUE
C
 1010 CONTINUE
C
C                               Calculate reaeration coefficient
C
      K20 = REAK*(AVVELE**EXPREV)*(AVDEPE**EXPRED)
      IF (K20 .GT. 24.) K20 = 24.
      K2HYDRA = K20*1.028**STP20
 1060 CONTINUE
      RETURN
      END

c****************************************************

        function oxysat(temp,salt)

c computes oxygen level at saturation

        implicit none

        real oxysat
        real temp,salt

        oxysat = (14.6244-0.367134*temp+4.4972E-3*temp**2
     :     -0.0966*salt+0.00205*salt*temp+2.739E-4*salt**2)
c     :    /32.

        end

c****************************************************

      SUBROUTINE BENTFLUX

C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:54.

      implicit none

      INCLUDE 'weutro.h'

      REAL AUX, FLUXP, FLUXN
C
C              Benthic ammonium and phosphate fluxes
C
       AUX = SA * 0.001   !conversion to g/day, fluxes are [mg/(m**2 day)]

         FLUXN = FNH4 (ISEG)*AUX
         CD (1, ISEG) = CD (1, ISEG) + FLUXN
         FLUXP = FPO4 (ISEG)*AUX
         CD (3, ISEG) = CD (3, ISEG) + FLUXP

c      write(77,*) FLUXN,ISEG,AUX,FNH4(ISEG)      !ggu
C
      RETURN
      END

c********************************************************************

        subroutine param_read

        implicit none
        INCLUDE 'weutro.h'

        integer ndim
        parameter(ndim=10)
        integer n_params
        real v_params(ndim)

        integer i

        open(1,file='michol.dat',status='old',form='formatted')
        read(1,*) n_params
        if( n_params .gt. ndim ) stop 'error stop: ndim - n_params'
        read(1,*) (v_params(i),i=1,n_params)
        close(1)

        K1C = v_params(1)
        K1D = v_params(2)

        write(6,*) 'using params in weutro from file michol.dat'

        end


c********************************************************************

      subroutine param_venezia

      implicit none
      INCLUDE 'weutro.h'

c        IS2 = 63       !optimum for Steele input Kj/m2/day, (or Kj/m2/300sec?dmk 20/9/2013)timestep 300 sec
        IS2 = 20        !optimum for Steele input Watt*h/mq/300sec (coherently with the previous line-dmk 20/9/2013)FixME
        IS2 = 200       !optimum for Steele input Watt*h/mq/300sec (coherently with the previous line-dmk 20/9/2013)FixME
        KESG(1)=0.85     !m-1 range: 0.1-5
c        KGRZ=0.8        !run07KGRZ= 1.2 primi 2 run 2013 mgl (inclusi nuovi_input)
        KGRZ=0.80
        KDZ = 0.150       ! 0.168 run08
c        K1C = 2.        !primi 2 run 2013
c        K1C = 2.5       ! test15 e test16 (squentin) 
c         K1C = 1.025      ! test17 (squentin)
         K1C = 0.45       ! test19 (squintin) originale 0.9

c        K58C = 0.44     ! test17 (squentin) increase(*2) remineralization of OP
c        k58C = 0.66     ! test19 (squintin) increase 10% remineralization of P
         k58C = 0.60     ! ultimo valore, provare 0.20
         KMNG1 = 0.125  ! test22 (squintin) decrease KMNG1(0.072),increase ammonia preference (PNH3G1)
         KMPG1= 0.096
      
        EFF = 0.5
        K1013C=0.037    !day-1 prova 4 agosto (primi 2 run)
c       K1013C =0.075   !test21(squintin) 0.035 descrease 50% remineralization DON (0.075 original)
         FON=1.0         !prova 4 agosto
c        FON=0.5        !prova 22 luglio 2015
         FOP=1.          !prova 4 agosto
c       K140C=0.18      !prova 11 agosto !che sia per questo che va a valori negativi? 24 sett 2013
        K140C=1.6      !increase denitrification at 1 d-1. San quintin value
     
        K1320C=0.010    !rate for nitrification 0.09-0.3
        KNO3=0.020          !mgO2/L       

        SOD1D(1)=0.1   

        write(6,*) 'using params in weutro for Venezia'

        end

c********************************************************************

      subroutine param_taranto

      implicit none
      INCLUDE 'weutro.h'

        K1C=2.88
        K1G=1.
        K1D=0.05
        KGRZ=1.8
        KDZ=0.1
        SOD1D(1)=0.2
        IS2=300/2.06
        KE(1)=0.5
        KESG(1)=1

        write(6,*) 'using params in weutro for Taranto'

        end

c********************************************************************

        subroutine weutro_debug(debug)

        implicit none

        logical debug

        include 'weutro.h'

        wdebug = debug

        end

c********************************************************************

        subroutine eutro_replay

        implicit none

        integer nstate
        parameter (nstate=9)

        integer id
        real t                  ![day]
        real dt                 !time step [day]
        real vol                ![m**3]
        real depth              ![m]
        real vel                ![m/s]
        real uws                ![m/s]
        real stp                ![Celsius]
        real sal                ![g/L], e.g., 30.
        real qss                !solar rad from shyfem
        real c(nstate)          ![mg/L]
        real cold(nstate)       !old state variable (return)
        real cds(nstate)        !source term [g/day]
        real loads(nstate)      !loading for c [g/(m**3 day)]

        integer i

        write(6,*) 'initializing...'
        call eutroini

        write(6,*) 'reading...'
        call read_debug_param(id,nstate,t,dt
     +                ,vol,depth,vel,stp,sal,qss,c,cold,cds,loads)

        call write_debug_param(id,nstate,t,dt
     +                ,vol,depth,vel,stp,sal,qss,c,cold,cds,loads)

        do i=1,nstate
          c(i) = cold(i)
        end do

        write(6,*) 'running...'
        call eutro0d(id,t,dt,vol,depth,vel,uws,stp,sal,qss,c,loads)
        write(6,*) 'finished...'

        call write_debug_param(id,nstate,t,dt
     +                ,vol,depth,vel,stp,sal,qss,c,cold,cds,loads)

        end

c********************************************************************

        subroutine read_debug_param(id,nstate,t,dt
     +                ,vol,depth,vel,stp,sal,qss,c,cold,cds,loads)

        implicit none

        integer id              !id of box
        integer nstate          !number of state variables
        real t                  ![day]
        real dt                 !time step [day]
        real vol                ![m**3]
        real depth              ![m]
        real vel                ![m/s]
        real stp                ![Celsius]
        real sal                ![g/L], e.g., 30.
        real qss                !solar rad from shyfem
        real c(nstate)          ![mg/L]
        real cold(nstate)       !old state variable (return)
        real cds(nstate)        !source term [g/day]
        real loads(nstate)      !loading for c [g/(m**3 day)]

        real tair,wspeed,turbid,fracday,itotal
        real ice,flnh4,flpo4
        logical tops,bots,seds

        integer ns,i

        read(5,*) id,ns,t,dt
        if( ns .ne. nstate ) stop 'error stop read_debug_param: nstate'

        read(5,*) vel,stp,sal,vol,depth,qss
        read(5,*) tair,wspeed,turbid,fracday,itotal
        read(5,*) ice,flnh4,flpo4
        read(5,*) tops,bots,seds
        read(5,*) (c(i),i=1,nstate)
        read(5,*) (cold(i),i=1,nstate)
        read(5,*) (cds(i),i=1,nstate)
        read(5,*) (loads(i),i=1,nstate)

        call set_internal_param(tair,wspeed,turbid,fracday,itotal
     +                  ,ice,flnh4,flpo4
     +                  ,tops,bots,seds)

        end

c***************************************************************

        subroutine write_debug_param(id,nstate,t,dt
     +                ,vol,depth,vel,stp,sal,qss,c,cold,cds,loads)

        implicit none

        integer id              !id of box
        integer nstate          !number of state variables
        real t                  ![day]
        real dt                 !time step [day]
        real vol                ![m**3]
        real depth              ![m]
        real vel                ![m/s]
        real stp                ![Celsius]
        real sal                ![g/L], e.g., 30.
        real qss                !solar rad from shyfem
        real c(nstate)          ![mg/L]
        real cold(nstate)       !old state variable (return)
        real cds(nstate)        !source term [g/day]
        real loads(nstate)      !loading for c [g/(m**3 day)]

        real tair,wspeed,turbid,fracday,itotal
        real ice,flnh4,flpo4
        logical tops,bots,seds

        integer i

        call get_internal_param(tair,wspeed,turbid,fracday,itotal
     +                  ,ice,flnh4,flpo4
     +                  ,tops,bots,seds)

        write(6,*) id,nstate,t,dt
        write(6,*) vel,stp,sal,vol,depth,qss
        write(6,*) tair,wspeed,turbid,fracday,itotal
        write(6,*) ice,flnh4,flpo4
        write(6,*) tops,bots,seds
        write(6,*) (c(i),i=1,nstate)
        write(6,*) (cold(i),i=1,nstate)
        write(6,*) (cds(i),i=1,nstate)
        write(6,*) (loads(i),i=1,nstate)

        end

c***************************************************************

        subroutine get_internal_param(tair,wspeed,turbid,fracday,itotal
     +                  ,ice,flnh4,flpo4
     +                  ,tops,bots,seds)

c gets internal parameters (that can be altered through functions)

        implicit none

        include 'weutro.h'

        real tair,wspeed,turbid,fracday,itotal
        real ice,flnh4,flpo4
        logical tops,bots,seds

        tair = airtmp
        wspeed = wind
        turbid = kesg(1)
        fracday = fday
        itotal = itot
        ice = xicecvr
        flnh4 = fnh4(1)
        flpo4 = fpo4(1)
        tops = topseg
        bots = botseg
        seds = sedseg

        end

c***************************************************************

        subroutine set_internal_param(tair,wspeed,turbid,fracday,itotal
     +                  ,ice,flnh4,flpo4
     +                  ,tops,bots,seds)

c sets internal parameters (that can be altered through functions)

        implicit none

        include 'weutro.h'

        real tair,wspeed,turbid,fracday,itotal
        real ice,flnh4,flpo4
        logical tops,bots,seds

        airtmp = tair
        wind = wspeed
        kesg(1) = turbid
        fday = fracday
        itot = itotal
        xicecvr = ice
        fnh4(1) = flnh4
        fpo4(1) = flpo4
        topseg = tops
        botseg = bots
        sedseg = seds

        end

c***************************************************************

c      subroutine weutro_test

c test drives weutro

c      integer nstate
c      parameter (nstate=9)

c       integer iddt,nt,iday,n,it
c        integer id
c        real vol,depth,vel,stp,sal
c        real dt,t
c        real rlux
c        real pi
c        real itot,fday,iinst,imax

c      real c(nstate)
c      real loads(nstate)

c       data c /0.05, 0.4, 0.01, 0.05, 2., 11., 0.2, 0.01, 0.015/
c        data loads /0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

c        id = 1234
c      iddt = 300
c        rlat = 45.              !latitude

c        pi = 4.*atan(1.)
c      vol = 1.e+6
c      depth = 1.
c      vel = 1.
c      stp = 20.
c      sal = 18.

c      nt = 86400 / iddt      !time steps per day
c        if( nt*iddt .ne. 86400 ) then
c            stop 'error stop: iddt must devide 86400'
c      end if
c      dt = iddt/86400.
c      t = 0.
c        it = 0

c       call eutroini
c       call luxlen_init('lux.dat')      

c      do iday=1,365

c          stp = 17. + 13. * cos(2.*pi*(iday-172)/365.)   !simple temp curve
c      write(6,*) 'ma qui ci arrivo?'

c        call luxlen(t,itot,fday)
c          itot = itot * (300./86400.)   !for data Donata
c          call get_radiation(iday,rlat,fday,itot,imax)
c          itot = itot / 8.              !just to have same magnitude

c        call wlight(fday,itot)

c        do n=1,nt
c            it = it + iddt
c          t = it / 86400.
c            rlux = 1.
c            call eutro0d(id,t,dt,vol,depth,vel,stp,sal,qss,c,loads)
c        end do

c          write(6,*) t,stp,c
c          write(75,'(11e14.6)') t,stp,c
c      end do

c      end


c**************************************************************
c**************************************************************
c**************************************************************
c restart files -> fill in real routines
c**************************************************************
c**************************************************************
c**************************************************************

        subroutine write_restart_eco(iunit)
        implicit none
        integer iunit
        integer nstate,nkn,i
        nstate = 0
        nkn = 0
        write(iunit) nstate,nkn
        end
        subroutine skip_restart_eco(iunit)
        implicit none
        integer iunit
        integer nstate,nkn,i
        read(iunit) nstate,nkn
        do i=1,nstate
          read(iunit)
        end do
        end
        subroutine read_restart_eco(iunit)
        implicit none
        integer iunit
        call skip_restart_eco(iunit)
        end





c***************************************************************

c        program main
c        call eutro_replay
c      call weutro_test
c        end

c***************************************************************

