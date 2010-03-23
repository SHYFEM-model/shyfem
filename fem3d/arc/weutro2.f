c
c weutro - EUTRO for SHYFEM
c
c********************************************************************
c
c original version until April 2002
c
c 19-04-2002	ggu&dmk	BUGFIX1, print_init, changes in ZOOP
c
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************

	subroutine param

	implicit none
	INCLUDE 'weutro.h'

c initialization of parameters

C       BYPASS OPTIONS FOR SYSTEMS 1-9. 1=BYPASS 0=SIMULATE
	
        SYSBY(1)=0
        SYSBY(2)=0
        SYSBY(3)=0
        SYSBY(4)=0
        SYSBY(5)=0
        SYSBY(6)=0
        SYSBY(7)=0
        SYSBY(8)=0
        SYSBY(9)=0


        PI = 3.14159

	ITYPE(1) = 1
	graztype=1  	!if graztype=1 simulate zoop. variable
			!if graztype=0  use wasp formulation
			!Set initial value of variable zoo=0

c--------------------------------------------------------------
c
c--------------------------------------------------------------

	
c	if CCHLX is variable in the segments, introduce each value else
        CCHL=30	!range 20-50
	CCHLX(1) = CCHL

        NUTLIM=1.       !nutrient limitation option 0=minimum, 1=multiplicative.
			! Default=0

c       parameters values
c
c--------------------------------------------------------------
c       subroutine phyto
c--------------------------------------------------------------
c------ growth and respiration

c	K1RT=1.045i	!wasp orig	
c	K1RC=0.125	!wasp orig
	K1RT=1.068	!adapting the curve to dejak model
	K1RC=0.096	!adapting the curve to dejak model
        K1T=1.068      
        K1C=1.5 	

c	KMPHYT=1.
	KMPHYT=0.
c
c------ decomposition
c
        KPZDC=0.02      !verify the value
        KPZDT=1.08	!Default=1.0

c
c--------nutrients limitation-
c
c        KMNG1=0.025	!for standard model application
			!use a large KMNG1 (when KMNG1=0 PNH3G1=1.0)
        KMNG1=0.05	!from dejak model Venice lagoon
        KMPG1=0.01       
c
c	F(10,10,1) !spatially variable dissolved fraction of inorganic P
	FOPO4=0.9 !spatially variable dissolved fraction of inorganic P
	   
c	attenzione: nel benthic layer asume valori diversi (0.045-0.001)
c--------grazing-- zooplankton is a forcing----------
c
        K1G=0.		
c        K1G=0.08		
c        K1D=0.02       !wasp orig
	K1D=0.12	!from dejak model Venice Lagoon 
c        K1D=0.2         !prova 12sett
	ZOO=0.7	
c	if fraction is variable in segments, input each value, else:
	ZOOSG(1)=1.
c--------grazing--- zooplankton as a variable---------
c	KGRZ=1.44	!max grazing rate day-1
        KGRZ=2.         !prova 12sett
	KPHYZ=0.5	!
c	EFF=0.5		!zoo-phyto digestion efficiency
	EFF=0.7
	KDZ=0.192	!zoo death (with excrection)
c--------------------------------------------------------------
c       subroutine organop, inorganp (P=phosphorous)
c--------------------------------------------------------------

        PCRB=0.025 !mgP/mgC
        FOP=0.5      	!unitless
        K58T=1.08	!unitless      
        K58C=0.22	!day-1 
        KOPDC=0.0004	!day-1
        KOPDT=1.08	!unitless
	KPO4=1.0	!microg. P/L

c--------------------------------------------------------------
c       subroutine organicn, ammonia,nitrate
c--------------------------------------------------------------

	NCRB=0.115	!mg N/mg C
       	FON=0.5     
c      	K1320C=0.1      !0.09-0.13 day-1 wasp orig 
	K1320C=0.05	!from dejak model lagoon of Venice
        K1320T=1.08	!unitless
        K140C=0.09	!day-1
        K140T=1.045	!unitless
        KNIT=2.		!mgO2/L 
        KNO3=0.1	!mgO2/L       
        K1013C=0.075    !day-1
        K1013T=1.08    !unitless
        KONDC=0.0004    !day-1     
        KONDT=1.08      !unitles

c-------------------------------------------------------------
c       subroutine CBOD
c--------------------------------------------------------------

	OCRB=32./12.	!mg O2/md C
        KDC=0.18       !0.16-0.21 day-1
        KDT=1.047      !unitless
        KDSC=0.000	!day-1
        KDST=1.08	!unitless
        KBOD=0.5	!mg N/L

c--------------------------------------------------------------
c       subroutine dissoxyg
c--------------------------------------------------------------

	WIND=3.       	!m/s
        AIRTMP=22.     	!C
        WTYPE=3.	!1=small 2=medim 3=large
	XICECVR=1.	!default: 1. no ice

        K2=0   !4.3		!4.1-4.7 day-1 if k2=0 then use kawind or kahydra      

c        REARSG(i)=0 	!Segment specific reareation rate constant
                        !REARSG  used when  rear is not calc from wind or hydro

			!FIXME
       SOD1D(1)=2.0	!g/m2-day 0.2-4.0 sediment oxygen demand for segment
       SODTA(1)=1.08 
c
c--------------------------------------------------------------
c       light limitation
c--------------------------------------------------------------
c
        LGHTSW=3        !LIGHT SWITCH: 1=Di Toro 2=Smith  3=Steele
c	IS2=50000.	!optimum light intensity for Steele lux/h
	IS2=1200000.	!optimum light intensity for Steele lux/day
c			!
c--------------------------------------------------------------
c       subroutine ditoro
c--------------------------------------------------------------
c
        FDAY=0.5
        IS1=300.	!(Ly/day) langleys/day
        IS1X(1)=IS1	!FIXME
c
c--------------------------------------------------------------
c       subroutine smith
c--------------------------------------------------------------c
c
        PHIMX=720.	!mg C/mole photon
        XKC=0.017	!m2/mg chl a
        ITOT=500.       !ly/day
        iav=0.9 *itot/FDAY
        DTDAY=0.5     	!??????verificare, serve a modulare il seno: day
        NEWDAY=1	!unknown????if G.E. 1 calcolo smith

c you can have up to 5 KE time functions here default is 1

	KE(1)=1.
	KEFN(1)=1
	KESG(1)=1.5	!m-1 range: 0.1-5

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

        subroutine eutro0d(t,dt,vol,depth,vel,stp,sal,c,loads)

c EUTRO 0-Dimensional

        implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )

	real t			!actual time [day]
        real dt                 !time step [day]
        real vol                !volume [m**3]
        real depth              !depth of box [m]
        real vel                !velocity [m/s]
        real stp                !temperature [C]
        real sal                !salinity [psu] == [per mille]
        real c(nstate)          !state variable [mg/L] == [g/m**3]
	real loads(nstate)	!loading for c [g/(m**3 day)]

        real cds(nstate)        !source term (right hand side) [g/day]
        integer icall
        save icall
        data icall / 0 /

        if( icall .eq. 0 ) then
          call settopseg(.true.)        !marks segment as surface
          call setbotseg(.true.)        !marks segment as bottom
          icall = 1
        end if

        call source(t,vel,stp,sal,vol,depth,c,cds)
	call load0d(cds,loads,vol)
        call euler(dt,vol,c,cds)

        end

c********************************************************************

	subroutine load0d(cds,loads,vol)

c integrate loadings

	implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )

	real cds(nstate)	!source term [g/day]
	real loads(nstate)	!loading for c [g/(m**3 day)]
	real vol		!volume of box [m**3]

	integer i

	do i=1,nstate
	  cds(i) = cds(i) + vol * loads(i)
	end do

	end

c********************************************************************

	subroutine euler(dt,vol,c,cds)

	implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )

	real dt			!time step [day]
	real vol		!volume [m**3]
	real c(nstate)		!state variable [mg/L] == [g/m**3]
	real cds(nstate)	!source term [g/day]

	integer i
	real volold,volnew	!volume [m**3]
	real mass		!mass [g]
	real mder		!derivative [g/day]

	volold = vol
	volnew = vol

	do i=1,nstate
	  mass = c(i) * volold
	  mder = cds(i)
	  c(i) = ( mass + dt * mder ) / volnew
	end do

c	 write (7, 21) c
c   21 format (9f8.2)

	end
	
c********************************************************************

	subroutine eutroini

c initializes eutro routines

	implicit none

        call param
        call EUTROINT
	call print_init

	end

c********************************************************************

	subroutine print_init

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

	subroutine source(t,vel0,stp0,sal0,vol0,depth0,c,cds)

c vel	velocity
c stp	temperature
	
	implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )

	real t			![day]
	real vel0		![m/s]
	real stp0		![Celsius]
	real sal0		![g/L], e.g., 30.
	real vol0		![m**3]
	real depth0		![m]
	real c(nstate)		![mg/L]		!FIXME
	real cds(nstate)	![g/day]	!FIXME


        INCLUDE 'weutro.h'
	
	integer i,n

	iseg = 1
	sedseg = .false.
	ito = 0
	
	daytime = t

cdd         IF (ITYPE (ISEG) .GE. 3.) SEDSEG = .TRUE.
cddd         ITO = IBOTSG (ISEG)

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
	h = depth0
	vol = vol0
        SA = vol0/depth0

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

	 CALL BENTFLUX

c adesso sono definiti i cd

	do i=1,nstate
	  cds(i) = cd(i,iseg)
	end do

	end

c********************************************************************
c********************************************************************
c********************************************************************

	subroutine wmeteo(tempair,windspeed)

c sets meteo parameters
c tempair is the air temperature C degrees
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

c sets turbidity parameter (units ?)

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
c fracday is the fraction of the day that is light
c itotal incident light intensity during the daylight langley/m2

	implicit none

	include 'weutro.h'

	real fracday
	real itotal

	fday = fracday
	itot = itotal
	itotmp = itotal

	iav=0.9 *itot/FDAY
c      IAV = 0.9*ITOTmp*RLGHTS (ISEG, 2)/FDAY
	write(6,*)'eutro.f',iav

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

c sets grazing parameter
c set in zoo the C/l value of grazing zooplankton

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

	write(6,*) 'sedflux : ',flnh4,flpo4

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

	subroutine segtype(segty)

c sets type of segment
c
c                        !1=surface water segment
c                        !2=subsurface water segment
c                        !3=upper bed segment
c                        !4=lower bed segment

	implicit none
	include 'weutro.h'
	integer segty

	itype(1)=segty

	end

c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************

      SUBROUTINE EUTROINT
C
C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:53.
C
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*
	implicit none
C
CDDD	changed PHIMAX ---> PHIMX
      INCLUDE 'weutro.h'
	real R0MIN
	parameter( R0MIN = 1.e-15 )


	real xarg,phimax
	integer i,j,initb,mxdmp,mxseg,noseg


C
C                             Initialize Values
C
	topseg = .true.
	botseg = .true.

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

      DO 1000 J = 1,noseg 
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
C      IF (KBOD .EQ. 0.0) KBOD = 1.00E-20
      IF (XARG .LT. R0MIN) KBOD = 1.00E-20
CCSC
      XARG = ABS(KNO3)
C      IF (KNO3 .EQ. 0.0) KNO3 = 1.00E-20
      IF (XARG .LT. R0MIN) KNO3 = 1.00E-20
CCSC
      XARG = ABS(KPO4)
C      IF (KPO4 .EQ. 0.0) KPO4 = 1.00E-20
      IF (XARG .LT. R0MIN) KPO4 = 1.00E-20
CCSC
      XARG = ABS(KNIT)
C      IF (KNIT .EQ. 0.0) KNIT = 1.00E-20
      IF (XARG .LT. R0MIN) KNIT = 1.00E-20
CCSC
      XARG = ABS(KMNG1)
C      IF (KMNG1 .EQ. 0.0) KMNG1 = 1.00E-20
      IF (XARG .LT. R0MIN) KMNG1 = 1.00E-20
CCSC
      XARG = ABS(KMPG1)
C      IF (KMPG1 .EQ. 0.0) KMPG1 = 1.00E-20
      IF (XARG .LT. R0MIN) KMPG1 = 1.00E-20
CCSC
      XARG = ABS(KMPHYT)
C      IF (KMPHYT .EQ. 0.0) KMPHYT = 1.00E-20
      IF (XARG .LT. R0MIN) KMPHYT = 1.00E-20
C
CCSC
      XARG = ABS(OCRB)
C      IF (OCRB .EQ. 0.) OCRB = 32./12.
      IF (XARG .LT. R0MIN) OCRB = 32./12.
CCSC
      XARG = ABS(IS1)
C      IF (IS1 .EQ. 0.) IS1 = 300.
      IF (XARG .LT. R0MIN) IS1 = 300.
CCSC
      XARG = ABS(CCHL)
C      IF (CCHL .EQ. 0.) CCHL = 30.
      IF (XARG .LT. R0MIN) CCHL = 30.
CCSC
      XARG = ABS(FON)
C      IF (FON .EQ. 0.) FON = 1.0
      IF (XARG .LT. R0MIN) FON = 1.0
CCSC
      XARG = ABS(FOP)
C      IF (FOP .EQ. 0.) FOP = 1.0
      IF (XARG .LT. R0MIN) FOP = 1.0
CCSC
      XARG = ABS(PHIMX)
C      IF (PHIMX .EQ. 0.) PHIMX = 720.
      IF (XARG .LT. R0MIN) PHIMX = 720.
CCSC
      XARG = ABS(XKC)
C      IF (XKC .EQ. 0.) XKC = 0.017
      IF (XARG .LT. R0MIN) XKC = 0.017
C
C  Check for Zero Temperature Correction Factors and readjust to 1.0
C
CCSC
      XARG = ABS(K1320T)
C      IF (K1320T .EQ. 0.) K1320T = 1.0
      IF (XARG .LT. R0MIN) K1320T = 1.0
CCSC
      XARG = ABS(K140T)
C      IF (K140T .EQ. 0.) K140T = 1.0
      IF (XARG .LT. R0MIN) K140T = 1.0
CCSC
      XARG = ABS(K1T)
C      IF (K1T .EQ. 0.) K1T = 1.0
      IF (XARG .LT. R0MIN) K1T = 1.0
CCSC
      XARG = ABS(K1RT)
C      IF (K1RT .EQ. 0.) K1RT = 1.0
      IF (XARG .LT. R0MIN) K1RT = 1.0
CCSC
      XARG = ABS(KDT)
C      IF (KDT .EQ. 0.) KDT = 1.0
      IF (XARG .LT. R0MIN) KDT = 1.0
CCSC
      XARG = ABS(K1013T)
C      IF (K1013T .EQ. 0.) K1013T = 1.0
      IF (XARG .LT. R0MIN) K1013T = 1.0
CCSC
      XARG = ABS(KONDT)
C      IF (KONDT .EQ. 0.) KONDT = 1.0
      IF (XARG .LT. R0MIN) KONDT = 1.0
CCSC
      XARG = ABS(K58T)
C      IF (K58T .EQ. 0.) K58T = 1.0
      IF (XARG .LT. R0MIN) K58T = 1.0
CCSC
      XARG = ABS(KOPDT)
C      IF (KOPDT .EQ. 0.) KOPDT = 1.0
      IF (XARG .LT. R0MIN) KOPDT = 1.0
CCSC
      XARG = ABS(KPZDT)
C      IF (KPZDT .EQ. 0.) KPZDT = 1.0
      IF (XARG .LT. R0MIN) KPZDT = 1.0
CCSC
      XARG = ABS(KDST)
C      IF (KDST .EQ. 0.) KDST = 1.0
      IF (XARG .LT. R0MIN) KDST = 1.0
C
C
C        Initialize internal clock for Dick Smith Light Formulation
C
C             Convert IC'S for Phyt from ug/L to mg C/L
C             Set all SOD Thetas to 1.0 if equal 0
C
c      DO 1010 ISEG = 1, NOSEG
CRBA--Date: Tuesday, 1 June 1993.  Time: 09:02:09.
cd         CCHLX(ISEG)=CCHL
c         IF (C (4, ISEG) .LT. 1.0e-24) C (4, ISEG) = 1.0e-24
c         C (4, ISEG) = C (4, ISEG)*CCHL/1000.
CCSC
c         XARG = ABS(SODTA(ISEG))
C         IF (SODTA (ISEG) .EQ. 0.) SODTA (ISEG) = 1.0
cd         IF (XARG .LT. R0MIN) SODTA (ISEG) = 1.0
c 1010 CONTINUE
C
cdd         SA = BVOL (ISEG)/DEPTHG (ISEG)
c	SA=bvol/depth
cd         FNH4 (ISEG) = FNH4 (ISEG)*SA*0.001
cd         FPO4 (ISEG) = FPO4 (ISEG)*SA*0.001
Cdd         ENDIF
 1050 CONTINUE

c sediment fluxes

	do i=1,segmax
          FNH4 (i) = 0.
          FPO4 (i) = 0.
	end do

c	write(6,*) segmax,FNH4(1),FPO4(1)

      RETURN
      END
      SUBROUTINE AMMONIA
C
C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:55.
C
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*
C
	implicit none

	real SR13ON,SR13P,SK13P1, denit

      INCLUDE 'weutro.h'
C
C       *-*-*-*-*  SYSTEM 1 - AMMONIA (NH3-N)  *-*-*-*-*
C
C                        Sources
C               Mineralization of organic nitrogen
C
      SR13ON = SK1013
c	denit = SK14D
	denit = 0.
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
         SK1314 = 0.
      END IF
      IF (STP .LT. 7.) SK1314 = 0.0
C
C                   Formulate Derivative
C
c	write(6,*) 'ammonia debug :'
c	write(6,*) SR13P,SR13ON,SK13P1,SK1314,VOL,CD (1, ISEG)
c	write(6,*) PNH3G1,NCRB,GP1,PHYT

      CD (1, ISEG) = (SR13P +denit+ SR13ON - SK13P1 - SK1314)*VOL

c	write(6,*) CD (1, ISEG),iseg
C
      RETURN
      END
      SUBROUTINE NITRATE
C
C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:56.
C
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*
C
	implicit none
      INCLUDE 'weutro.h'
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
      IF (SK14D .LT. 1.00E-24) SK14D = 1.00E-24
C
C                   Formulate Derivative
C
      CD (2, ISEG) = (SR1413 - SK14P1 - SK14D)*VOL
C
      RETURN
      END
      SUBROUTINE INORGANP
C
C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:56.
C
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*
C
	implicit none
      INCLUDE 'weutro.h'
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
C
C                   Formulate Derivative
C
      CD (3, ISEG) = (SR8P + SR8OP - SK8P)*VOL
C
      RETURN
      END
      SUBROUTINE PHYTO
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*
C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:57.
C
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*
C
	implicit none

      INCLUDE 'weutro.h'
	real R0MIN
	parameter( R0MIN = 1.e-15 )
      REAL XARG,XDIFF,GIT1,CN,DOPO4
	real ttr, ttr1
c	write (44,*), phyt,zoo	!ggu
C
C        *-*-*-*-*-*  System 4 - Phytoplankton  *-*-*-*-*-*
C
c         open (9, file = 'outgp.dat')	!ggu


c      write(6,*) 'phyto, sedseg,iseg : ',SEDSEG,iseg
      IF (SEDSEG) THEN


C
         GP1 = 0.0
         GPP = 0.0
         RESP = 0.0
         DP1 = KPZDC*KPZDT**STP20
         DPP = DP1*PHYT
C
C                      Growth Rate
C
      ELSE


c      write(6,*) 'phyto,  ctrl: ',K1C,PHYT,K1T,STP20
c      write(6,*) 'phyto, ctrl: ',GITMX1,KEFN(1)
         IF (SYSBY (4) .EQ. 1) THEN
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
c	write(6,*)KMNG1,PNH3G1,NH3,NO3,CN
            XEMP1 = CN/(KMNG1 + CN)
cd            DOPO4 = OPO4*F (2, 3, ISEG)
		DOPO4 = OPO4*FOPO4
            XEMP2 = DOPO4/(KMPG1 + DOPO4)
	
C
C       Compute Growth Rate Reduction due to Nutrient Limitation
C
CCSC
            XARG = ABS(NUTLIM)
c            IF (NUTLIM .EQ. 0.) RNUTR = AMIN1 (XEMP1, XEMP2)
            IF (XARG .LT. R0MIN) RNUTR = AMIN1 (XEMP1, XEMP2)

c	write(6,*) '-------------- ggu --------------'
c	write(6,*) NUTLIM,XARG,R0MIN
c	write(6,*) XEMP1,XEMP2,RNUTR
c	write(6,*) '-------------- ggu --------------'
CCSC
cd            XDIFF = NUTLIM - 1.0
cd            XARG = ABS(XDIFF)
            IF (NUTLIM .EQ. 1.) RNUTR = XEMP1*XEMP2
cd            IF (XARG .LT. R0MIN) RNUTR = XEMP1*XEMP2
c		write(6,*) XARG 
c		write(6,*) RNUTR
            GP1 = RNUTR*GIT1
            GPP = GP1*PHYT
         END IF		!SYSBY
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
c	if (PHYT .eq. 0) then
c	write(6,*)'phyto ctrl, phyt=' phyt
c               stop 'error stop'
C
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C********************************************
C                     Respiration Rate
C********************************************
C
         RESP = K1RC*K1RT**STP20
c	
C
C       ALGAL RESPIRATION + DEATH + GRAZING
C
c         DP1 = RESP + K1D + K1G*ZOO*ZOOSG(ISEG)	!old

         DP1 = RESP + K1D 		!FIXMED 

c         DPP = DP1*PHYT	!BUGFIX1
         DPP = DPP + DP1*PHYT
c         RESP = RESP*PHYT	!BUGFIX1

c	write(6,*) DP1,RESP,K1D,K1G,ZOO,ZOOSG(iseg)


C         IF (PHTY .GT. 1.0E-6)THEN
	XEMPRC=1.
Cdoni            XEMPRC = PHYT/(KMPHYT + PHYT)
C         ELSE
C            XEMPRC = 1.0E-6/(KMPHYT + 1.0E-6)
C         ENDIF

C
      END IF	!sedsed
C
c	write(6,*) 'phyto debug ggu'
c	write(6,*)RNUTR,XEMP1,XEMP2 

c	write(36,'(8(f8.4,2x))')K1C,K1T,STP20,GITMX1,RLIGHT,RNUTR,GP1,DP1
c        write(9,'(3(f8.4,2x))') GP1,DP1, RESP 	!ggu

c      CD (4, ISEG) = (GP1 - DP1-GRZ)*PHYT*VOL	!BUGFIX1
      CD (4, ISEG) = (GPP - DPP - GRZ)*VOL
      RETURN
      END
	subroutine zoop
	implicit none
	include 'weutro.h'
	real gra
c
c
c
c	*-*-*-*- system 9 zooplankton *-*-*-
c
c	Sources: zooplankton growth
c

c
c zoosk	Sink term: zooplankton death -->source term for organop organicn cbodsv
c	    grazing inefficiency -> source term for  organop organicn cbodsv
c

	if( graztype.eq.1 ) then
	  DPP = 0.
          GRZ=KGRZ*ZOO*PHYT/(PHYT+KPHYZ)
	  zoosk = (1-EFF)*GRZ + KDZ*ZOO
	  CD (9, ISEG) = (EFF*GRZ - KDZ*ZOO)*VOL
	else	!BUGFIX1
	  DPP = K1G*ZOO*PHYT*ZOOSG(ISEG)	!original formulation
	  GRZ = 0.
	  zoosk = 0. 
	  CD (9, ISEG) = 0.
	end if

	return
	end


      SUBROUTINE CBODSV
C
C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:57.
C
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*
C
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
         IF (DO .GT. 1.0E-15) THEN
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
C
cddd      IF (CBOD .GT. 0.)TEMPBOD= SK180/CBOD
      SK18D = (5./4.)*(32./14.)*SK14D
C
C                      Formulate Derivative
C
      CD (5, ISEG) = (SR18P+GRC - SK180 - SK18D)*VOL
C
      RETURN
      END
      SUBROUTINE DISSOXYG
C
C     Last Revised:  Date: Monday, 26 August 1991.  Time: 10:37:53.
C
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*
C
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
      IF (ITYPE (ISEG) .EQ. 1 .AND. XICECVR .GT. 0.0) THEN

c if surface element and not completly covered with ice

CCSC
         XARG1 = ABS(K2)
cd         XARG2 = ABS(REARSG(ISEG))
cd	XARG2=0.
cd         IF (K2 .EQ. 0. .AND. REARSG(ISEG) .EQ. 0.0) THEN
         IF (K2 .EQ. 0.0) THEN
c         IF (XARG1.LT.R0MIN.AND.XARG2.LT.R0MIN) THEN
            K2WSAVE = 0.
            K2HSAVE = 0.
            CALL KAWIND (WIND, STP, AIRTMP, H, WTYPE,  K2WIND) 
            CALL KAHYDRA (K2HYDRA)
ccdd	write(6,*)'dissoxyg:',K2HYDRA,K2WIND
            K2WSAVE = K2WIND
            K2HSAVE = K2HYDRA
            IF (K2WIND .GT. K2HYDRA) THEN
               KA = K2WIND * XICECVR
            ELSE
               KA = K2HYDRA * XICECVR
            END IF
         ELSE
            IF (K2 .GT. 0)THEN
                KA = ((K2*1.028**STP20)* XICECVR)
            ELSE
cd                KA = REARSG(ISEG) * REAR
		KA = 0.
            ENDIF
c elimino rear e rearsg perché non ci interessano, sono una time function 
c e una segment specific reareation rate, ma noi usiamo la reareation from
c  kawind o kahydra o imponendo un valore della costante di reareazione
         ENDIF
      ELSE
         KA = 0.0
      END IF
C
C       Calculate oxygen saturation level for current water
C       temperature; DOSAT is expressed as mg oxygen per liter
C

      CL = SAL/1.80655
      TK = STP + 273.
      RLNCS = - 139.34411 + (1.575701E05/TK) - (6.642308E07/TK**2) +
     1   (1.243800E10/TK**3) - (8.621949E11/TK**4) -
     2   (CL*(3.1929E-02 - (19.428/TK) + (3.8673E03/TK**2)))
C
      CS = EXP (RLNCS)
C
      SR190 = KA*(CS - DO)
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
C
c      SK19P = OCRB*RESP			!BUGFIX1
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
C=======================================================================
C                     Formulate Derivative
C
      CD (6, ISEG) = (SR190 + SR19PA + SR19PB - SK19P
     1   - SK1913 - SK1918 - SK19S)*VOL
c        write(54,'(8(f8.4,2x))') SR190,SR19PA,SR19PB,SK19P,SK1913,SK1918
c     & ,SK19S 
c        write(54,'(3(f8.4,2x))') SR19PA,SR19PB,GP1	!ggu
       

C
C
      RETURN
      END
      SUBROUTINE ORGANICN
C
C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:59.
C
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*
C
	implicit none
      INCLUDE 'weutro.h'
	real SR10P,GRN
C
C
C             *-*-*-*-*  SYSTEM 7 Organic Nitrogen      *-*-*-*-*
C
C                          Sources
C             Phytoplankton Respiration and 'DEATH'
C		GRAZING and zooplankton death
C                  GRN=  ZOO --> ORGANICN

C
      SR10P = NCRB*DPP*FON
      GRN   = ZOOSK*NCRB	!BUGFIX1
C
C                         Sinks
C         Mineralization of Dissolved Organic Nitrogen
C
      IF ( .NOT. SEDSEG) SK1013 = (K1013C*K1013T**STP20)*ON*XEMPRC
      IF (SEDSEG) SK1013 = (KONDC*KONDT**STP20)*ON
C
C                   Formulate Derivative
C
      CD (7, ISEG) = (SR10P+ GRN - SK1013)*VOL
C
      RETURN
      END
      SUBROUTINE ORGANOP
C
C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:59.
C
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*
C
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
      GRP  = ZOOSK * PCRB	!BUGFIX1
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
      CD (8, ISEG) = (SR5P + GRP - SK58)*VOL
C
      RETURN
      END
      SUBROUTINE DITORO
C*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
C        Di Toro et al Light Formulation
C*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
C
	implicit none
      INCLUDE 'weutro.h'
	real CCHL1,TCHLA,KESHD,SKE
	real temp1, temp2, temp3
         open (10, file = 'outtoro.dat')

C
      CCHL1 = CCHL
      TCHLA = PHYT/CCHL1
      KESHD = (0.0088*1000.*TCHLA + 0.054*(1000.*TCHLA)**0.6667)
      SKE = KESG (ISEG)
      IF (IKE .GT. 0 .AND. IKE .LE. 5) SKE = SKE*KE (IKE)
      SKE = SKE + KESHD
      TEMP1 = SKE*H
C
C         Get average solar radiation during daylight hours
C
ccd
      TEMP2 = IAV/IS1
      TEMP3 = EXP ( - TEMP1)
      IAVBOT=IAV * TEMP3
      RLGHTS (ITO, 2) = RLGHTS (ISEG, 2)*TEMP3
      RLIGHT = 2.718*FDAY/TEMP1*(EXP ( - TEMP2*
     1   TEMP3) - EXP ( - TEMP2))
      RLGHTS (ISEG, 1) = RLIGHT

c	write(6,*) 'ditoro debug ggu'
c	write(6,*) H,SKE,TEMP1,TEMP2,TEMP3,IAVBOT
c	write(6,*) FDAY,RLGHTS (ISEG, 1),RLIGHT
c	write(6,*) ISEG,ITO,IS1,IAV, 2.718*FDAY,temp1
c	write(6,*) H
        write(10,'(6(f8.5,2x))')KESHD,SKE,PHYT,RLIGHT

C
C
      RETURN
      END
      SUBROUTINE SMITH
C*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
C          Dick Smith variable carbon/chlorophyll ratio
C*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
C
	implicit none
      INCLUDE 'weutro.h'

	real I0,IMAX,IAVSG,SUM,KESHD,SKE
c	real IAVBOTX(segmax)
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
         KESHD = XKC*1000.*TCHLA	!FIXME
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
      SUBROUTINE KAWIND (WS, TW, TA, depth, WTYPE,  RK)
C*******************************************************************
C*
C*     SUBROUTINE REAERK CALCULATES:
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
C*=====================================================================
C*
C* THIS SUBROUTINE WAS WRITTEN BY:
C*     SANDRA BIRD
C*     USAE WATERWAYS EXPERIMENT STATION (WES-ES-Q)
C*     VICKSBURG, MISSISSIPPI
C* AND MODIFIED BY JAMES L. MARTIN
C*
C*
C*=====================================================================
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
C*

	implicit none

	real R0MIN
	parameter( R0MIN = 1.e-15 )
c	INCLUDE 'weutro.h'
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
		write(6,*) 'iwtype: ',iwtype
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
      SUBROUTINE KAHYDRA (K2HYDRA) 
C
C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:58.
C
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*
C
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
C
C     Last Revised:  Date: Thursday, 1 February 1990.  Time: 16:32:54.
C
C*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*X*
C
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

c	write(77,*) FLUXN,ISEG,AUX,FNH4(ISEG)	!ggu
C
      RETURN
      END

c***************************************************************

	subroutine steele

	implicit none

	include 'weutro.h'

	real TCHLA,KESHD,SKE
	real rlim

        TCHLA = PHYT/CCHL
        KESHD = (0.0088*1000.*TCHLA + 0.054*(1000.*TCHLA)**0.6667)

        SKE = KESG (ISEG)
        IF (IKE .GT. 0 .AND. IKE .LE. 5) SKE = SKE*KE (IKE)
        SKE = SKE + KESHD

	RLIGHT = rlim(daytime,IS2,SKE,H)

c	write(6,*) 'steele: ',daytime,IS2,SKE,H
c	write(6,*) PHYT,CCHL,KESHD,RLIGHT

	end

c***************************************************************

	function rlim(t,is,ke,depth)

c computes light limitation (formula 5.6 in WASP user manual)

	implicit none

	real rlim	!computed light limitation
	real t		!day [0-365]
	real is		!saturation light intensity
	real ke		!light extinction coefficient [1/m]
	real depth	!depth [m]

	real expon
	parameter ( expon = 2.718 )

	real lux,fday
	real i0,iobottom
	real temp1,temp2,temp3
	real rintens

	call luxlen(t,lux,fday)
	i0 = rintens(t,lux,fday)

c	write(6,*) 'rlim: ',t,lux,fday,ke,depth,is,i0

	temp1 = ke * depth
	temp2 = i0 / is
	temp3 = exp ( - temp1 )

	iobottom = i0 * temp3		!for later ...

c	rlim = (expon/temp1) * ( exp(-temp2*temp3) - exp(-temp2) )
c	steele formulation for rlim:

	rlim=temp2*temp3*exp(1-temp2*temp3)
c	write(6,*) 'rlim(1): ',rlim,temp1,temp2,temp3
c            write(72,'(9(f13.4,2x))')I0,rlim	!ggu


	end
	
c***************************************************************

	subroutine luxlen(t,lux,fday)

c given a day t returns light intensity lux and day length fday

	implicit none

	real t		!day [0-365]
	real lux	!light intensity [units??]
	real fday	!day length [day]

	real luxv(365)
	real fdayv(365)
	save luxv, fdayv

	integer n
	real tdummy
	integer it2n

	integer icall
	save icall
	data icall / 0 /

c--------------------------------------------------------------
c read in file
c--------------------------------------------------------------

	if( icall .eq. 0 ) then

          open(2,file='../input/lux1.dat',status='old',form='formatted')
	  do n=1,365
            read(2,*) tdummy,luxv(n),fdayv(n)
c	write(72,*) tdummy,luxv(n),fdayv(n)	!ggu
	  end do
	  close(2)

c	  do n=1,365
c            fdayv(n) = 18./24.	!FIXME
c            fdayv(n) = 12./24.	!FIXME
c	  end do
	  icall = 1
	end if

c--------------------------------------------------------------
c compute values
c--------------------------------------------------------------

	n = it2n(t)

	lux = luxv(n)
	fday = fdayv(n)

	end

c***************************************************************

	function rintens(t,intens,fday)

c computes light intensity during a day (formula 5.7 in WASP user manual)

	implicit none

	real rintens	!light intensity at time t
	real t		!day [0-365]
	real intens	!total light over one day
	real fday	!day length [day]
	
	real pi
	parameter( pi = 3.14159 )

	real tday,aux,i0,aux2

	tday = t - int(t)			!fraction in day
	if( tday .lt. 0. ) tday = tday + 1.	!negative days
	tday = tday - 0.5			!maximum at noon

	aux = pi / fday
	aux2=intens/(24*fday)

	if( abs(tday) .le. fday/2 ) then
	  i0 = (intens/2.) * aux * cos( tday * aux )
	else
	  i0 = 0.
	end if

	rintens = i0

c	rintens = max( i0 , 0. )		!only positive values
c	write (87,*)I0
	end

c***************************************************************

	function it2n(t)

c converts time [day] to pointer into array [1-365]

	implicit none

	integer it2n
	real t

	integer n

	n = t
	if( t .lt. 0. ) n = n - 1
	n = mod(n,365)
	if( n .lt. 0 ) n = n + 365	!handle negative days
	n = n + 1			!n is in [1,365]

	it2n = n

	end

c***************************************************************


c********************************************************************

	subroutine rdtemp(t,temp)

c reads temperature file

	implicit none

	real t		!actual time [day]
	real temp	!temperature [C] (return)

	real tempv(365)
	save tempv

	integer i,n
	integer it2n

        integer icall,nold
        save icall,nold
        data icall,nold / 0 , 0 /

	if( icall .eq. 0 ) then

c         open(1,file='temperature.con',status='old',form='formatted')
         open(1,file='../input/temperature1.FIX',status='old',
     &form='formatted')

         write(6,*) '../temperature file opened...'
         do i=1,365
           read(1,*) temp

	    tempv(i) = temp
	  end do
	  close(1)
	  write(6,*) 'temperature file read and closed...'
	  icall = 1
	end if

	n = it2n(t)

	temp = tempv(n)

        if( n .ne. nold ) then
          nold = n
          write(6,*) 'New temperature read for temp ',n,t,temp
        end if

c	write(6,*) 'rdtemp: ',t,n,temp

	end

c********************************************************************

