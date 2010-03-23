C****************************************************************************
C
C                  SED96:  A SEDIMENT TRANSPORT MODEL
C                           FOR CONTINENTAL SHELVES
C
C                  GEOLOGICAL SURVEY OF CANADA (ATLANTIC)
C                     BEDFORD INSTITUTE OF OCEANOGRAPHY
C
C****************************************************************************
C
C Revision Log:
C
C Jan, 2003	ccf	(sedtg1.f)
C			based on sedg6g.f
C Jan, 2003	ccf	(sedtg2.f)
C			made it easy, deleted old routine, call SEDISDGM
C			made subroutine ZVEL and GETIOPT1
C Jan, 2003	ccf	(sedtg3.f)
C			deleted not used variable RK
C Apr, 2003	ccf	(sedi3d_2.f)
C			bh,sutr,betr variable only in nodes
C			grainsize reade from .str file
C			variable introduced in subsys.f
C Apr, 2003	ccf	(sedi3d_3.f)
C			use only bh, deleted sutr,betr
C			use SED instead of SEDM
C Apr, 2003	ccf 	(sedi3d_4.f)
C			different computaion for suspended mean
C			concentration
C May, 2003	ccf	(sedi3d_5.f)
C			deposition calculated in subroutine massconc
C Jul, 2003	ccf	(sedi3d_6.f)
C			implemented wave module (subwaves.f)
C Aug, 2003	ccf	(sedi3d_7.f)
C			included cohesive sediment
C Sep, 2003	ccf	(sedi3d_8.f)
C			rearrange wave call out
C Nov, 2003	ccf	(sedi3d_9.f)
C			new method using the deposition - erosion term
C			rearrange wave call out
C Feb, 2004	ccf	(sedi3d_10.f)
C			bug fix in calcoulus of bh
C			bug fix in n2e3d
C			compute from ceq from getceq
C Feb, 2004	ccf	(sedi3d_11.f)
C			read CONC0 (initial concentration from .str file)
C			correct definition of the erosion-deposition term
C Mar, 2004	ccf	(sedi3d_12.f)
C			read physical parameter only once
C Mar, 2004	ccf	(sedi3d_13.f)
C			if UZ=0 then Z=0. to avoid problem in FRICFAC
C Mar, 2004	ccf	(sedi3d_14.f)
C			add mass balence check
C			CONC0 in mg/l
C Mar, 2004	ccf	(sedi3d_15.f)
C			bug fix in declared physical parameters
C Apr, 2004     ccf     (sedi3d_16.f)
C			adjust wave common, add ed and wws common
C Apr, 2004     ccf     (sedi3d_17.f)
C                       deleted routine getceq
C Ago, 2004     ccf     (sedi3d_18.f)
C			get C0A from sedtrans (the detph averaged reference 
C			concentration)
C Ago, 2004     ccf     (sedi3d_19.f)
C                       first advection and diffusion, then sedtrans, then settling
C			and bottom flux
C Ago, 2004     ccf     (sedi3d_20.f)
C                       some changes
c 05.03.2004    ggu     integrated changes into main model
c
c 14.01.2004    ggu     fixed error in calling n2e3d
c 14.01.2004    ggu     new call to scal3sh -> made 3D arrau difhv
C
C****************************************************************************
C
C AVAILABLE OPTIONS SEDIMENT TRANSPORT PREDICTOR (IOPT1) ARE:
C   1 - ENGELUND-HANSEN (1967) TOTAL LOAD EQUATION (GD > 0.15 mm)
C   2 - EINSTEIN-BROWN (1950) BEDLOAD EQUATION (0.3 < GD < 29 mm)
C   3 - BAGNOLD (1963) TOTAL LOAD EQUATION (0.18 < GD <0.45 mm)
C   4 - YALIN (1963) BEDLOAD EQUATION (GD > 0.2 mm)
C   7 - AMOS AND GREENBERG (1980) COHESIVE SEDIMENTS ( GD < 0.063 mm)


      subroutine sedi(it,idt)

      implicit none

      include 'param.h'

      INTEGER IOPT1             !SEDIMENT TRANSPORT FORMULA OPTION NUMBER
      DOUBLE PRECISION D        !WATER DEPTH (M)
      DOUBLE PRECISION UZ       !AMBIENT CURRENT AT HEIGHT Z ABOVE THE SEAFLOOR (M/S)
      DOUBLE PRECISION Z        !HEIGHT OF UZ ABOVE SEAFLOOR
      DOUBLE PRECISION CDIR     !DIRECTION OF AMBIENT CURRENT (DEGREES TRUE)
      DOUBLE PRECISION HT       !WAVE HEIGHT (M)
      DOUBLE PRECISION PER      !WAVE PERIOD (S)
      DOUBLE PRECISION WDIR     !WAVE PROPOGATION DIRECTION (DEGREES TRUE)
      DOUBLE PRECISION GD       !SEDIMENT GRAIN DIAMETER (M)
      DOUBLE PRECISION RLINP    !INPUT RIPPLE LENGTH (M)
      DOUBLE PRECISION RHINP    !INPUT RIPPLE HIGHT (M)
      DOUBLE PRECISION BETA     !BED SLOPE (DEGREE)
      DOUBLE PRECISION RHOS     !DENSITY OF SEDIMENT MINERAL(S) (KG/M**3)
      DOUBLE PRECISION RHOW     !DENSITY OF FLUID (WATER)  (KG/M**3)
      DOUBLE PRECISION FRACT    !FRACTION OF THE TOTAL SEDIMENT WITH GRAIN SIZE GD
      DOUBLE PRECISION CONC0    !INITIAL ESTIMATE OF SEDIMENT CONCENTRATION (ie mg/l)
      DOUBLE PRECISION TAOCE    !CRITICAL STRESS FOR EROSION (Pa)
      DOUBLE PRECISION TAOCD    !CRITICAL STRESS FOR DEPOSITION (Pa)
      DOUBLE PRECISION TIMEDR   !DEPOSITION OR EROSION DURATION (minutes)
      DOUBLE PRECISION WS       !SETTLING VELOCITY FOR COHESIVE SEDIMENT (m/s)
      DOUBLE PRECISION PRS      !PROBABILITY OF RESUSPENSION (NORMALLY ASSUMED = 0)
      DOUBLE PRECISION RKERO    !PROPORTIONALITY COEFF FOR EROSION RATE (DEFAULT = 1.62)
      DOUBLE PRECISION ZB0      !INITIALE DEPTH OF THE LAYER EXPOSED TO THE FLOW
      DOUBLE PRECISION AULVA    !PERCENTAGE OF AREA COVERED BY THE ALGAE 'ULVA' (%)
      DOUBLE PRECISION WSULVA   !SETTLING VELOCITY OF ULVA (M/S)

C------------- OUTPUT VARIABLES -----------------

      DOUBLE PRECISION UB       !MAX WAVE INDUCED ORBITAL VELOCITY AT THE BOTTOM (M/S)
      DOUBLE PRECISION AB       !EXCERSION LENGTH OF BOTTOM WAVE ORBIT (M)
      DOUBLE PRECISION WL       !WAVE LENGTH (M)
      DOUBLE PRECISION FCW      !BOTTOM (SKIN) FRICTION FACTOR
      DOUBLE PRECISION DELTACW  !HEIGHT OF THE WAVE-CURRENT BOUNDARY LAYER
      DOUBLE PRECISION UA       !CURRENT SPEED USED IN BOTTOM STRESS CALC. (M/SEC)
      DOUBLE PRECISION PHIB     !ANGLE BETWEEN WAVE AND CURRENT DIRECTIONS (RADIANS)
      DOUBLE PRECISION U100     !CURRENT SPEED AT 1 M. ABOVE SEABED (M/SEC)
      DOUBLE PRECISION USTCS    !CURRENT SKIN-FRICTION SHEAR VELOCITY OF GM
      DOUBLE PRECISION USTWS    !WAVE SKIN-FRICTION SHEAR VELOCITY OF GM
      DOUBLE PRECISION USTCWS   !COMBINED SKIN-FRICTION SHEAR VELOCITY OF GM
      DOUBLE PRECISION USTCWSE  !EFFECTIVE COMBINED SKIN-FRICTION SHEAR VELOCITY
      DOUBLE PRECISION USTC     !TOTAL CURRENT SHEAR VELOCITY OF GM
      DOUBLE PRECISION USTW     !TOTAL WAVE SHEAR VELOCITY OF GM
      DOUBLE PRECISION USTCW    !COMBINED TOTAL SHEAR VELOCITY OF GM
      DOUBLE PRECISION Z0       !BED ROUGHNESS LENGTH (M)
      DOUBLE PRECISION Z0C      !APPARENT BED ROUGHNESS LENGTH (M)
      DOUBLE PRECISION RHEIGHT  !PREDICTED RIPPLE HEIGHT
      DOUBLE PRECISION RLENGTH  !PREDICTED RIPPLE LENGTH
      DOUBLE PRECISION USTCRB   !CRITICAL SHEAR VEL OF BEDLOAD TRANS (M/SEC)
      DOUBLE PRECISION USTCRS   !CRITICAL SHEAR VEL OF SUSPENDED LOAD TRANSPORT (M/SEC)
      DOUBLE PRECISION TB1      !TIME AT WHICH BEDLOAD TRANSPORT CEASES (SEC) 
      DOUBLE PRECISION TB2      !TIME AT WHICH BEDLOAD TRANSPORT RECOMMENCES (SEC)
      DOUBLE PRECISION TS1      !TIME AT WHICH SUSPENDED LOAD TRANSPORT CEASES (SEC)
      DOUBLE PRECISION TS2      !TIME AT WHICH SUSPENDED LOAD TRANSPORT RECOMMENCES (SEC)
      DOUBLE PRECISION PERBED   !% OF TIME SPENT IN ONLY BEDLOAD TRANSPORT PHASE
      DOUBLE PRECISION PERSUSP  !% OF TIME SPENT IN SUSPENDED LOAD TRANSPORT PHASE
      DOUBLE PRECISION QS       !SUSPENDED SEDIMENT TRANSPORT RATE (KG/M/S)
      DOUBLE PRECISION QSDIR    !DIRECTION OF SUSPENDED SEDIMENT TRANSPORT (DEGREE)
      DOUBLE PRECISION C0	!REFERENCE CONCENTRATION AT Z0 (KG/M^3)
      DOUBLE PRECISION C0A	!DEPTH AVERAGED REFERENCE CONCENTRATION AT Z0 (KG/M^3)
      DOUBLE PRECISION SED      !TIME-AVERAGED NET SEDIMENT TRANSPORT AS VOLUME (M**3/S/M)
      DOUBLE PRECISION SEDM     !TIME-AVERAGED NET SEDIMENT TRANSPORT AS MASS (KG/S/M)
      DOUBLE PRECISION SEDDIR   !DIRECTION OF NET SEDIMENT TRANSPORT (AZIMUTH,DEGREES)
      DOUBLE PRECISION CONC     !FINAL CALCULATED SEDIMENT CONCENTRATION (mg/l)
      DOUBLE PRECISION RD0      !INITIAL DEPOSITION RATE (kg/m^2/s)
      DOUBLE PRECISION RD       !DEPOSITION RATE (kg/m^2/s)
      DOUBLE PRECISION RE0      !INITIAL EROSION RATE (kg/m^2/s)
      DOUBLE PRECISION RE       !EROSION RATE (kg/m^2/s)
      DOUBLE PRECISION TIME0    !CALCULATED TIME (minutes) WHEN CONC. BE LESS THAN 1 MG/L
      DOUBLE PRECISION QS0      !INITIAL MASS SEDIMENT TRANSPORT RATE (kg/m/s)
      DOUBLE PRECISION PHIM     !INTERNAL FRICTION ANGLE
      DOUBLE PRECISION TAOCWS   !AVERAGED EFFECTIVE SHEAR STRESS
      DOUBLE PRECISION TAOS     !SOLID TRANSMITTED STRESS DUE TO ULVA
      DOUBLE PRECISION RES      !ULVA RESUSPENSION PARAMETER
      DOUBLE PRECISION USTCWSB  !TRANSPORT-RELATED COMBINED SHEAR VELOCITY
      DOUBLE PRECISION USTUP    !CRITICAL SHEAR VEL OF SHEET FLOW TRANSPORT (M/SEC)
      DOUBLE PRECISION FALL     !SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
      real USTAR		!SHARE VELOCITY [M/S]

c -------------------------------------------------------------
c fem variables
c -------------------------------------------------------------

      integer ndim		!number of sediment grain size class
      parameter (ndim=1)
      logical cohesive

      integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
      common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
      real eps1,eps2,pi,flag,high,higi
      common /mkonst/ eps1,eps2,pi,flag,high,higi

      integer nen3v(3,1)	!node number
      common /nen3v/nen3v
      integer ilhv(1),ilhkv(1)	!number of element and node level
      common /ilhv/ilhv, /ilhkv/ilhkv
      real ev(13,1)		!element information vector ??
      common /ev/ev
      real rcv(1)               !boundary ??
      common /rcv/rcv
      
      real ulnv(nlvdim,1),vlnv(nlvdim,1)!x and y integrated velocity 3d
      common /ulnv/ulnv, /vlnv/vlnv
      real hdenv(nlvdim,neldim) 	!element depth 3d
      common /hdenv/hdenv

      integer idt,dt		!time step in seconds
      integer it		!time in seconds
      integer isedi		!call parameter
      integer ius,itmcon,idtcon !output parameter
      integer ie,ii,k,l		!counter
      integer levmax		!number of element level
      integer icall

      real wws			!settling velocity [m/s]
      real u,v			!x and y current components
      real b,c			!x and y derivated form function
      real bflux		!flux of bedload mass [m/s]
      real sflux		!flux of suspended mass [m/s]
      real dep,area,vol		!depth,area and volume of box [m,m2,m3]
      real por			!porosity
      real bh(nkndim)    	!bottom height variation [m]
      real deparea		!deposition area [m2]
      real ceq			!depth-averaged equilibrium suspended concentration [kg/m3]

      real waveh(neldim)        !wave height [m]
      real wavep(neldim)        !wave period [s]
      real waved(neldim)        !wave direction (same as wind direction)
      common  /waveh/waveh, /wavep/wavep, /waved/waved

      double precision ede		!deposition-erosion per element [kg/m2s]
      double precision edn(nkndim)	!deposition-erosion per node [kg/m2s]

      real sce(nlvdim,neldim)	!suspended concentration per element [kg/m3]
      real scn(nlvdim,nkndim)	!suspended concentration per node [kg/m3]
      real rhe(neldim)          !ripple height [m]
      real rle(neldim)          !ripple length [m]

      real aaux(nlvdim,neldim)
      real difv(0:nlvdim,1)
      common /difv/difv
      real difhv(nlvdim,1)
      common /difhv/difhv
      real rkpar,difmol
      character*10 what
      real uuz,ccdir,ssed,sedir
      real aux1,aux

C function
      integer iround
      real getpar
      real areaele
      real getceq

C save and data
      save ius,itmcon,idtcon
      save rkpar,difmol
      save dt,gd,cohesive,IOPT1

      save icall
      data icall /0/

c ----------------------------------------------------------
c Initialization
c ----------------------------------------------------------

        if( icall .le. -1 ) return

        if( icall .eq. 0 ) then
          isedi = iround(getpar('isedi'))
          if( isedi .le. 0 ) icall = -1
          if( icall .le. -1 ) return
          icall = 1

c         --------------------------------------------------
c         Sediment parameter
c         --------------------------------------------------

          GD = getpar('sedgrs')			!sediment grain size
          CONC0 = iround(getpar('sedref'))	!initial concentration

          call getiopt1(GD,IOPT1)		!get IOPT1
          cohesive = iopt1 .eq. 7		!true=cohesive; false=non cohesive

c         --------------------------------------------------
c	  Initialize state variables
c         --------------------------------------------------

          wws = 0.
          ede = 0.

          do k = 1,nkn
            bh(k) = 0.
            edn(k) = 0.
            levmax = ilhkv(k)
            do l=1,levmax
              scn(l,k) = CONC0
            end do
          end do

          do ie = 1,nel
            levmax = ilhv(ie)
            do l = 1,levmax
              sce(l,ie) = 0.
              aaux(l,ie) = 0.
            end do
            rhe(ie) = 0.
            rle(ie) = 0.
          end do      
             
c         --------------------------------------------------
c         Parameters for transport/diffusion resolution
c         --------------------------------------------------

          rkpar=getpar('chpar')
          difmol=getpar('difmol')

	  dt = idt

c         --------------------------------------------------
c	  Initialize output
c         --------------------------------------------------

          ius = 61
          itmcon = iround(getpar('itmcon'))
          idtcon = iround(getpar('idtcon'))

          call confop(ius,itmcon,idtcon,ndim,'sed')

          write(6,*) 'sediment model initialized...'

	endif

c -------------------------------------------------------------------
c normal call
c -------------------------------------------------------------------

        what = 'sedtrans'

c       -------------------------------------------------------------
c       Phisical parameter
c       -------------------------------------------------------------

        BETA  = 0.			! [degree]
        FRACT = 1.			! 

        ZB0    = 0.001			! [m]
        AULVA  = 0.			! [%]
        WSULVA = 0.			! [m/s]

        CONC0  = 0.			! [mg/l]
        TAOCE  = 0.005			! [Pa]
        TAOCD  = 0.002			! [Pa]
        TIMEDR = 5.			! [minutes]
        WS     = 0.0002			! [m/s]
        PRS    = 0.		 
        RKERO  = 1.62

        RHOW = 1023.			! [kg/m**3]

        if ( cohesive ) then
          RHOS = 1150.			! [kg/m**3]
        else
          RHOS = 2650.
        endif

        por = 0.
        aux1 = dt/(1-por)

c       -------------------------------------------------------------------
c       transport and diffusion
c       -------------------------------------------------------------------

        call scal3sh(what,scn(1,1),1,rcv,rkpar,difhv,difv,difmol)
        
        call n2e3d(nlvdim,sce,scn,aaux)		!from node to element conc.

c       -------------------------------------------------------------------
c       start loop on elements
c       -------------------------------------------------------------------
	
	do ie = 1,nel

          HT = waveh(ie)			!wave height
          PER = wavep(ie)			!wave period
          WDIR = waved(ie)			!wave direction

          l = ilhv(ie)				!bottom element level
          D = hdenv(l,ie)			!depth of element
          area = areaele(ie) 			!area of element
          vol = area * D			!volume of element
          u = ulnv(l,ie)			!u current velocity component
          v = vlnv(l,ie)			!v current velocity component
          RHINP = rhe(ie)			!ripple height
          RLINP = rhe(ie)			!ripple length

          call c2p(u,v,uuz,ccdir)		!get UZ and CDIR
          UZ = uuz
          CDIR = ccdir
          call zvel(UZ,D,GD,Z,USTAR)		!get Z

          if (cohesive) CONC0 = sce(l,ie)*1000.	!cohesive initial concentration [mg/l]

          call sedisdgm(D,UZ,Z,CDIR,HT,PER,WDIR,GD,RHINP,RLINP,
     @BETA,RHOS,RHOW,IOPT1,FRACT,CONC0,TAOCE,TAOCD,TIMEDR,
     @WS,PRS,RKERO,ZB0,AULVA,WSULVA,               !**INPUT VARIABLES**

     @UB,AB,WL,FCW,DELTACW,UA,U100,PHIB,USTCS,USTWS,USTCWS,USTCWSB,
     @USTC,USTW,USTCW,Z0,Z0C,RHEIGHT,RLENGTH,USTCRB,USTCRS,TS1,TB1,
     @TS2,TB2,PERBED,PERSUSP,QS,QSDIR,C0,C0A,SEDM,SED,SEDDIR,CONC,RD0,
     @RD,USTUP,FALL,USTCWSE,RE0,RE,TIME0,QS0,PHIM,TAOCWS,TAOS,RES)!**OUTPUT VARIABLES**

          rhe(ie) = RHEIGHT
          rle(ie) = RLENGTH

c         -----------------------------------------------------------------
c         compute suspended erosion-deposition rate "ede" [kg/m2s]
c         -----------------------------------------------------------------
c         if ede > 0 erosion occurs, if ede < 0 deposition occurs

          if ( cohesive ) then
            wws = WS
            ede = RE - RD
            SED = 0.
          else
            wws = FALL
            ede = wws*(C0A - sce(l,ie))
          endif

c         -----------------------------------------------------------------
c         compute the sediment transport on each element
c         -----------------------------------------------------------------
c         Using the sediment continuity equation:
c             dh/dt + 1/(1-p)*((dqx(b)/dx + dqy(b)/dy) + ED) = 0
c         h = bottom hight [m]
c         p = porosity
c         qx(b),qy(b) = x and y volumetric bedload transport rate [m3/ms]
c         ED = (erosion - deposition) rate (suspended sediment) [kg/m2s]

          SEDDIR = SEDDIR / (45. / atan (1.))			!rad

          do ii = 1,3
            k = nen3v(ii,ie)
            b = ev(3+ii,ie)                             	!1/m
            c = ev(6+ii,ie)                             	!1/m
            bflux = 3*SED*(-b*sin(SEDDIR)-c*cos(SEDDIR))	!m/s
            sflux = ede*area/3.					!kg/m2s

c           ---------------------------------------------------------------
c           compute bottom elevation change due to bedload
c           ---------------------------------------------------------------

            bh(k) = bh(k) + bflux*aux1 				!m
            edn(k) = edn(k) + sflux		         	!kg/s
          end do

c       -------------------------------------------------------------------
c       end of element loop
c       -------------------------------------------------------------------

        end do

c	-------------------------------------------------------------------
c	commmon for cohesive and non-cohesive sediments
c	-------------------------------------------------------------------

        call settl(dt,scn,wws,edn)		!settling & sediment flux

        do k = 1,nkn
          bh(k) = bh(k) - edn(k)*aux1/rhos	!bottom height variation due
        end do					!to suspended transport 

        call massbalance(it,rhos,scn,bh)	!check sed mass conservation

c       -------------------------------------------------------------------
c       write of results (file SED)
c       -------------------------------------------------------------------

        call confil(ius,itmcon,idtcon,80,1,bh(1))
        call confil(ius,itmcon,idtcon,81,nlvdim,scn(1,1))

c -------------------------------------------------------------------
c end of routine
c -------------------------------------------------------------------

	end

C ********************************************************************
C *******************************************************
C SUBROUTINE ZVEL
C *******************************************************
C THIS SUBROUTINE COMPUTES Z, THE HEIGHT OF UZ (THE MAIN CURRENT 
C VELOCITY) ABOVE SEAFLOOR, ASSUMING A LOG PROFILE OF THE VELOCITY 
C ON THE DEPTH, BASED ON THE LAW: u = (ustar/k)*ln(z/zo)

	SUBROUTINE ZVEL(UZ,D,GD,Z,USTAR)

	IMPLICIT NONE

	DOUBLE PRECISION D	!WATER DEPTH (M)
	DOUBLE PRECISION UZ	!AMBIENT CURRENT AT HEIGHT Z ABOVE THE SEAFLOOR (M/S)
	DOUBLE PRECISION GD	!SEDIMENT GRAIN DIAMETER (M)
	DOUBLE PRECISION Z	!HEIGHT OF UZ ABOVE SEAFLOOR

	REAL USTAR		!SHEAR VELOCITY
	REAL BR			!BOTTOM ROUGHNESS HEIGHT
	REAL Z0			!BOTTOM ROUGHNESS LENGTH
	REAL K			!VON KARMAN COSTANT
	PARAMETER(K=0.4)

        BR = 2.5*GD		!SEE AMOS ARTICLE
        Z0 = BR/30.
        USTAR = (UZ*K*D/Z0)/(D/Z0*(LOG(D/Z0)-1)+1)
        Z = Z0*EXP(UZ*K/USTAR)
        if (UZ .eq. 0) Z=D/2.
	
	END

C ********************************************************************
C *******************************************************
C SUBROUTINE GETIOPT1
C *******************************************************
C THIS SUBROUTINE HANDLE THE SEDIMENT TRANSPORT FORMULA IN
C FUNCTION OF THE GRAIN SIZE. AVAILABLE OPTIONS FOR IOPT1 ARE:
C   1 - ENGELUND-HANSEN (1967) TOTAL LOAD EQUATION (GD > 0.15 mm)
C   2 - EINSTEIN-BROWN (1950) BEDLOAD EQUATION (0.3 < GD < 29 mm)
C   3 - BAGNOLD (1963) TOTAL LOAD EQUATION (0.18 < GD <0.45 mm)
C   4 - YALIN (1963) BEDLOAD EQUATION (GD > 0.2 mm) 
C   7 - AMOS AND GREENBERG (1980) COHESIVE SEDIMENTS ( GD < 0.063 mm)

	SUBROUTINE GETIOPT1(GD,IOPT1)

	IMPLICIT NONE

	INTEGER IOPT1		!SEDIMENT TRANSPORT FORMULA OPTION NUMBER	
	DOUBLE PRECISION GD	!SEDIMENT GRAIN DIAMETER (M)

	IF ( GD .GE. 0.001) THEN			!GRANULO
	   IOPT1 = 2
	ELSEIF (GD.LT.0.001.AND.GD.GE.0.0005) THEN	!SABBIA GROSSA
	   IOPT1 = 2
	ELSEIF (GD.LT.0.0005.AND.GD.GE.0.00025) THEN	!SABBIA MEDIA
	   IOPT1 = 3
	ELSEIF (GD.LT.0.00025.AND.GD.GE.0.000125) THEN	!SABBIA FINE
	   IOPT1 = 2
	ELSEIF (GD.LT.0.000125.AND.GD.GT.0.000063) THEN	!SABBIA M. FINE
	   IOPT1 = 2
	ELSEIF (GD .LE. 0.000063 ) THEN			!SILT E ARGILLA
	   IOPT1 = 7
	ENDIF

	END

C ********************************************************************
C *******************************************************
C SUBROUTINE SETTL
C *******************************************************
C THIS SUBROUTINE COMPUTES THE CONTRIBUTE OF THE SEDIMENT SETTLING AND 
C OF THE SEDIMENT FLUX BEETWEN WATER COLUMN AND BOTTOM

        subroutine settl(dt,scn,wws,edn)

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhv(1),ilhkv(1)  	!number of element and node level
        common /ilhv/ilhv, /ilhkv/ilhkv

        integer dt			!time step [s]
        real scn(nlvdim,nkndim)   	!suspended concentration per node [kg/m3]
        real wws                  	!settling velocity [m/s]
        double precision edn(nkndim)  	!deposition-erosion per node [kg/m2s]

        integer k,l
        real depnode,dep
        real areanode,area

        do k = 1,nkn
          l = ilhkv(k)
          dep = depnode(l,k,1)
          area = areanode(l,k)                !area of node
          aux = dt/dep

          edn(k) = edn(k)/area
          if (scn(l,k).lt.0.001 .and. edn(k).lt.0.) then
            edn(k) = -scn(l,k)/aux
            scn(l,k) = 0.
          else
            scn(l,k) = scn(l,k) + edn(k)*aux
          end if
        end do

        end

C ********************************************************************
C compute the sediment mass balance between suspended material and what
C is erosded/deposited (for closed basin)

        subroutine massbalance(it,rhos,scn,bh)

        implicit none

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real hdenv(nlvdim,neldim)         !element depth 3d
        common /hdenv/hdenv
        integer ilhv(1),ilhkv(1)  !number of element and node level
        common /ilhv/ilhv, /ilhkv/ilhkv
	
        real scn(nlvdim,nkndim)   !suspended concentration per node [kg/m3]
        real bh(nkndim)           !
        double precision rhos     !density of sediment mineral(s) (kg/m**3)

        real areanode,volnode
        real vol,area
        real somb,soms,relat

        integer k,l,nlev,it

        soms = 0.
        somb = 0. 

        do k = 1,nkn
          nlev = ilhkv(k)
          do l = 1,nlev
            area = areanode(l,k)                !area of node
            vol = volnode(l,k,+1)               !volume of node
            soms = soms + scn(l,k)*vol		!suspended contribute
          end do 
          somb = somb + bh(k)*rhos*area         !bottom contribute
        end do

        if (somb .eq. 0.) then
          relat = 0.
        else 
          relat = abs( (soms + somb)/somb )
        end if
        write(22,*)it,relat,soms,somb,soms + somb

c        if (relat .gt. 1.e-3) then
c          write(*,*)'sediment mass balance error:',relat
c        end if

        end
