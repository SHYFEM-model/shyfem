c
c $Id: subqfxu1.f,v 1.7 2009-03-11 16:25:59 georg Exp $
c
c heat flux module (utility 1)
c
c contents :
c
c subroutine rh2twb(db,rh,wb)           computes wb from db and rh
c subroutine twb2rh(db,wb,rh)		computes rh from db and wb
c subroutine psy(isel,db,rh,wb,dp)      convert between values
c ...                                   others
c subroutine wbt(t,rh,td,tw)		computes td, wb from ta, rh
c
c revision log :
c
c 01.02.2002    ggu     new import into file
c 09.12.2002    ggu     radiation routines taken out
c 28.03.2003    ggu     new routine wbt substituted in rh2twb
c 30.07.2003    ggu     new routine test_wbt added
c 01.02.2006    ggu     new routine twb2rh(db,wb,rh)
c 21.08.2007    ggu     bug fix for wbt -> account for rh = 0
c 18.02.2009    ggu     dew point introduced into psy()
c 16.02.2011    ggu     conversion from specific humidity (not working)
c
c*****************************************************************************

        subroutine rh2twb(db,rh,wb)

c computes wet bulb temperature from temperature and relative humidity

	implicit none

	real db		!(dry) air temperature [C]
	real rh		!relative humidity [%] (0-100)
	real wb		!wet bulb temperature [C] (output)

	integer isel	!select the mode
        real dp         !dew point temperature (not used)

        !isel = 3
        !call psy(isel,db,rh,wb,dp)

        call wbt(db,rh,dp,wb)

	end

c*****************************************************************************

        subroutine twb2rh(db,wb,rh)

c computes relative humidity from temperature and wet bulb temperature

	implicit none

	real db		!(dry) air temperature [C]
	real wb		!wet bulb temperature [C]
	real rh		!relative humidity [%] (0-100) (output)

	integer isel	!select the mode
	real dp		!dew point

        isel = 1
        call psy(isel,db,rh,wb,dp)

	end

c*****************************************************************************

      subroutine psy(isel,db,rh,wb,dp)

c grandezze psicrometriche 
c
c temperatures in Celsius, pressure in Pascal
c
C     DB Temp. a bulbo asciutto 0C
C     WB Temp. a bulbo umido 0C,
C     DP Temp. di rugiada 0C
C     PB Press. atmosferica N/m2,
C     PV Press. parz. di vapore N/m2
C     W Umidita' specifica Kg/Kg,
C     H Entalpia KJ/Kg
C     V Volume specifico m3/Kg
C     RH Umidita' relativa % [0-100]

      REAL H
  
      real pstd
      parameter ( pstd = 1013.0 )

      PB=pstd*100.			!std pressure [pa]
      
      IF(ISEL .EQ. 0) GO TO 9	!nothing

      IF(ISEL .EQ. 1) GO TO 1	!input db,wb
      IF(ISEL .EQ. 2) GO TO 2	!input db,dp
      IF(ISEL .EQ. 3) GO TO 3	!input db,rh
      IF(ISEL .EQ. 4) GO TO 4	!input dp,h
      IF(ISEL .EQ. 5) GO TO 5	!input db,w
   
    1 CALL PSDBWB(DB,WB,DP,PB,PV,W,H,V,RH)
      GO TO 9
    
    2 CALL PSDBDP(DB,WB,DP,PB,PV,W,H,V,RH)
      GO TO 9
    
    3 CALL PSDBRH(DB,WB,DP,PB,PV,W,H,V,RH)
      GO TO 9
    
    4 CALL PSDPH(DB,WB,DP,PB,PV,W,H,V,RH)
      GO TO 9
    
    5 CALL PSDBW(DB,WB,DP,PB,PV,W,H,V,RH)
      GO TO 9
      
    9 CONTINUE
      
      END
       
c*****************************************************************************

	subroutine convert_humidity(mode,rh,w,tair)

c converts specific from/to relative humidity

	implicit none

	integer mode	! 0: relative to specific  1: specific to relative
	real rh		! relative humidity [%]
	real w		! specific humidity
	real tair	! air temperature [Celsius]

	real DB,WB,DP,PB,PV,H,V

	db = tair

	if( mode .eq. 0 ) then
	  CALL PSDBRH(DB,WB,DP,PB,PV,W,H,V,RH)
	else if( mode .eq. 1 ) then
	  CALL PSDBW(DB,WB,DP,PB,PV,W,H,V,RH)
	else
	  stop 'error stop humidity: mode value not allowed'
	end if

	end

c*****************************************************************************
c
c	subroutine rh2w(rh,w)
c
c	w = 0.622 * rh * rhows / ( rho - rhows) * 100.
c
c    w = specific humidity of air vapor mixture (kg/kg)
c    rh = relative humidity (%)
c    rhows = density of water vapor (kg/m3)
c    rho = density of the moist or humid air (kg/m3)
c
c	end
c
c*****************************************************************************
       
c the next routine compute a series of values from two input values (and PB)
c
c example: PSDBWB uses DB and WB to compute DP,PV,W,H,V,RH

      SUBROUTINE PSDBWB(DB,WB,DP,PB,PV,W,H,V,RH)
      PVP=PVSF(WB)
      IF (DB-WB)1,4,7
    1 WRITE(*,10)DB,WB
   10 FORMAT(5X,'ERROR IN PSDBWB : DB=',F6.2,'  WB=',F6.2
     *                  ,'   WB greater than DB')
      STOP
    7 HDB0=1.006*DB
      WB1=WBFF(0.0,DB)
      IF (WB .GE. WB1) GO TO 5
      WRITE(*,11)DB,WB
   11 FORMAT(5X,'ERROR IN PSDBWB : DB=',F6.2,'  WB=',F6.2
     *                  ,'   temperatures not compatible')
      STOP
    4 PV=PVP
      GO TO 3
    5 WSTAR=0.62198*PVP/(PB-PVP)
      HL=2501+1.83*DB-4.194*WB
      CH=1.006*(WB-DB)+WSTAR*(2501.-2.364*WB)
      EX=CH/(HL*0.62198)
      PV=ABS(PB*EX/(1.+EX))
    3 W=0.62198*PV/(PB-PV)
      V=8314.32/28.9645*(DB+273.16)*(1.+1.6078*W)/PB
      RH=PV/PVSF(DB)*100.
      DP=DPF(PV)
      H=DB*1.006+(2501+1.83*DP)*W
      RETURN
      END
      
      SUBROUTINE PSDBDP(DB,WB,DP,PB,PV,W,H,V,RH)
      IF(DP .LE. DB) GO TO 1
      WRITE(*,10)DB,DP
   10 FORMAT(5X,'ERROR IN PSDBDP : DB=',F6.2,'  DP=',F6.2
     *                  ,'   DP greater than DB')
      STOP
    1 PV=PVSF(DP)
      PVS=PVSF(DB)
      RH=PV/PVS*100.
      W=0.62198*PV/(PB-PV)
      V=8314.32/28.9645*(DB+273.16)*(1.+1.6078*W)/PB
      H=DB*1.006+(2501+1.83*DB)*W
      WB=WBFF(W,DB)
      RETURN
      END
      
      SUBROUTINE PSDBRH(DB,WB,DP,PB,PV,W,H,V,RH)
      IF (RH .EQ. 0.0)RH=0.001
      IF (RH .GE. 0.001 .AND. RH .LE. 100) GO TO 1
      WRITE(*,11)DB,RH
   11 FORMAT(5X,'ERROR IN PSDBRH : DB=',F6.2,'  RH=',F6.2
     *                  ,'   relative humidity not compatible')
      STOP
    1 PVS=PVSF(DB)
      X=RH-100.
      IF (ABS(X)-0.09)10,10,20
   10 RH=100.
   20 PV=RH/100.*PVS
      W=0.62198*PV/(PB-PV)
      V=8314.32/28.9645*(DB+273.16)*(1.+1.6078*W)/PB
      H=DB*1.006+(2501+1.83*DB)*W
      DP=DPF(PV)
      WB=WBFF(W,DB)
      RETURN
      END
      
      SUBROUTINE PSDPH(DB,WB,DP,PB,PV,W,H,V,RH)
      PV=PVSF(DP)
      W=0.62198*PV/(PB-PV)
      H2=(1.006*DP+(2501.+1.83*DP)*W)-0.0001
      IF(H .GE. H2)GO TO 1
      WRITE(*,10)DP,H
   10 FORMAT(5X,'ERROR IN PSDPH : DP=',F6.2,'  H=',F6.3
     *                  ,'   data not compatible')
      STOP
    1 DB=(H-2501*W)/(1.006+1.83*W)
      WB=WBFF(W,DB)
      V=8314.32/28.9645*(DB+273.16)*(1.+1.6078*W)/PB
      PWS=PVSF(DB)
      RH=PV/PWS*100.
      RETURN
      END
      
      SUBROUTINE PSDBW(DB,WB,DP,PB,PV,W,H,V,RH)
      IF(W .EQ. 0.0)W=1.E-07
      IF(W .GE. 1.E-07)GO TO 2
      WRITE(*,11)DB,W
   11 FORMAT(5X,'ERROR IN PSDBW : DB=',F6.2,' W=',E10.3
     *                  ,'   specific humidity not compatible')
      STOP
    2 XS=0.62198*PVSF(DB)/(PB-PVSF(DB))
      IF (W .LE. XS)GO TO 3
      WRITE(*,11)DB,W
      STOP
    3 H=1.006*DB+(2501.+1.83*DB)*W
      V=8314.32/28.9645*(DB+273.16)*(1.+1.6078*W)/PB
      PV=W*PB/(0.62198+W)
      RH=PV/PVSF(DB)*100.
      WB=WBFF(W,DB)
      DP=DPF(PV)
      RETURN
      END
      
      FUNCTION DPF(PV)
      PU=PV/3386.389
      Y=ALOG(PU)
      IF(PU .GT. 0.18036)GO TO 1
      DPF=22.21111+13.81833*Y+0.49594*Y*Y
      GO TO 2
    1 DPF=26.13722+16.98833*Y+1.04961*Y*Y
    2 RETURN
      END
      
      FUNCTION PVSF(TEMP)
      real pstd
      parameter ( pstd = 1013.0 )
      PRIF=pstd*100.
      T=TEMP+273.16
      Z=273.16/T
      IF (T .LT. 273.16)GO TO 3
      P1=10.79586*(1.-Z)-2.2195983
      P2=5.02808*ALOG10(Z)
      A1=-8.29692*((1./Z)-1.)
      P3=1.5047E-04*(1.-10.**A1)
      A2=4.76955*(1.-Z)
      P4=0.42873E-03*((10.**A2)-1.)
      GO TO 4
    3 P1=-9.096936*(Z-1.)
      P2=-3.56654*ALOG10(Z)
      P3=0.876817*(1.-1./Z)
      P4=-2.2195983
    4 SUM=P1+P2+P3+P4
      PVSF=PRIF*10.**SUM
      RETURN
      END
      
      FUNCTION WBFF(W,DB)
      WB1=DB
   15 WS1=XSAT(WB1)
      W1=(WS1*(2501-2.364*WB1)+1.006*(WB1-DB))/(2501+1.83*DB-4.194*WB1)
      Y1=W-W1
      IF(ABS(Y1)-0.00003)11,11,16
   16 IF(Y1)9,11,14
   14 WB1=WB1+0.5
      GO TO 15
    9 WB2=WB1-1
      WS2=XSAT(WB2)
      W2=(WS2*(2501-2.364*WB2)+1.006*(WB2-DB))/(2501+1.83*DB-4.194*WB2)
      Y2=W-W2
      IF(ABS(Y2)-0.00003)10,10,20
   20 IF(Y1*Y2)6,7,8
    8 WB1=WB2
      Y1=Y2
      GO TO 9
    7 IF(Y1)10,11,10
   11 WBFF=WB1
      GO TO 4
   10 WBFF=WB2
      GO TO 4
    6 Z=ABS(Y1/Y2)
      WBFF=(WB2*Z+WB1)/(1.+Z)
    4 RETURN
      END
      
      FUNCTION XSAT(TEMP)
      DIMENSION C(8)
      DATA C/0.15036368E-02,0.87755907E-04,0.17497041E-04,-0.11118891E-5
     *,0.56661491E-07,-0.13317505E-08,0.16671702E-10,-0.80607734E-13/
      TEMP=TEMP+11.
      SUM=0.0
      DO 10 II=1,8
      I=9-II
      SUM=C(I)+SUM*TEMP
   10 CONTINUE
      XSAT=SUM
      TEMP=TEMP-11
      RETURN
      END

c*****************************************************************************

        subroutine wbt(t,rh,td,tw)

c computes dew point and wet bulb temperature from dry temp and rel. hum.

        implicit none

        real t          !dry air temperature
        real rh         !relative humidity
        real td         !dew point air temperature
        real tw         !wet bulb temperature (out)

        real es0
        parameter(es0=0.611)    !reference saturation vapor pressure in kPa
        real t0
        parameter(t0=237.3)     !reference temperature ?

        real p,taux,es,e,lne,gamma,delta

        p = 100                 !pressure in kPa

        es = es0 * exp(17.27*t/(t+t0))
        e = 0.01 * rh * es      !e in kPa

c	--------------------------------------------------------
c	next if to handle rh = 0.
c	--------------------------------------------------------

	if( e .lt. 1.e-2 ) then
	  td = -t0
	  tw = t
	  return
	end if

        lne = log(e)
        td = ( 116.9 + t0 * lne ) / ( 16.78 - lne )

        gamma = 0.00066 * p

        taux = td
        delta = 4098 * e / (taux + t0)**2
        tw = ( gamma * t + delta * td ) / ( gamma + delta )

        taux = 0.5 * ( td + tw )
        delta = 4098 * e / (taux + t0)**2
        tw = ( gamma * t + delta * td ) / ( gamma + delta )

        end

c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************

	subroutine test_wbt2

	implicit none

	integer i
	real db,rh,wb,dp,td,wb2,wb3

        db = 12.

	do i=100,0,-1
	  rh = i
          call rh2twb(db,rh,wb)
          call wbt(db,rh,td,wb2)
          call psy(3,db,rh,wb3,dp)
	  write(6,*) rh,db,wb,wb2,wb3,td
	end do

	end

c***********************************************************************

	subroutine test_wbt

c -14.70000       76.00000

        db = -12.
        rh = 81.

        db = -14.7
        rh = 76.

        isel = 3

    1   continue
        read(5,*,end=2) it,qs,ta,rh,wind,cc

        db = ta
        call rh2twb(db,rh,wb)
        call wbt(db,rh,td,tw)
        call psy(isel,db,rh,wb,dp)

        write(6,*) db,rh,td,wb,tw
        goto 1
    2   continue

        end

c***********************************************************************

	subroutine test_specific

	tair = 25.
	rh = 70.
	call convert_humidity(0,rh,w,tair)
	write(6,*) tair,w,rh
	call convert_humidity(1,rh1,w,tair)
	write(6,*) tair,w,rh,rh1

	tair = -1.
	tair = 25.
	w = 2.e-3
	call convert_humidity(1,rh,w,tair)
	call convert_humidity(0,rh,w1,tair)
	write(6,*) tair,rh,w,w1

	rh = 70.
	tair = 20.
	call convert_humidity(0,rh,w,tair)
	call convert_humidity(1,rh1,w,tair)
	write(6,*) tair,w,rh,rh1

	end

c***********************************************************************

c	program main_test_wbt
c	call test_wbt
c	call test_wbt2
c	call test_specific
c	end

c***********************************************************************

