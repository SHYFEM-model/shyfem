c
c $Id: subqfx1.f,v 1.2 2003/06/20 15:42:22 georg Exp $
c
c helper functions for heat fluxes
c
c contents :
c
c revision log :
c
c 28.03.2003    ggu     new function wbt for wet bulb temperature (converges)
c
******************************************************************************

        subroutine rh2twb(db,rh,wb)

c computes wet bulb temperature from temperature and relative humidity

	implicit none

	real db		!(dry) air temperature [C]
	real rh		!relative humidity [%] (0-100)
	real wb		!wet bulb temperature [C] (output)

	integer isel	!select the mode
	real td		!dew point temperature (not used)

	!isel = 3
        !call psy(isel,db,rh,wb)

        call wbt(db,rh,td,wb)

	end

******************************************************************************

      subroutine psy(isel,db,rh,wb)

c grandezze psicrometriche 
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
  
      PB=101325			!std pressure [pa]
      
      IF(ISEL .EQ. 0) GO TO 4	!nothing
      IF(ISEL .EQ. 1) GO TO 5	!db,wb -> rh
      IF(ISEL .EQ. 2) GO TO 6
      IF(ISEL .EQ. 3) GO TO 7	!db,rh -> wb
      IF(ISEL .EQ. 4) GO TO 8
      IF(ISEL .EQ. 5) GO TO 9
   
    5 CALL PSDBWB(DB,WB,DP,PB,PV,W,H,V,RH)
      GO TO 4
    
    6 CALL PSDBDP(DB,WB,DP,PB,PV,W,H,V,RH)
      GO TO 4
    
    7 CALL PSDBRH(DB,WB,DP,PB,PV,W,H,V,RH)
      GO TO 4
    
    8 CALL PSDPH(DB,WB,DP,PB,PV,W,H,V,RH)
      GO TO 4
    
    9 CALL PSDBW(DB,WB,DP,PB,PV,W,H,V,RH)
      GO TO 4
      
    4 CONTINUE
      
      END
       
******************************************************************************
       
c the next routine compute a series of values from two input values (and PB)
c
c example: PSDBWB uses DB and WB to compute DP,PV,W,H,V,RH

      SUBROUTINE PSDBWB(DB,WB,DP,PB,PV,W,H,V,RH)
      PVP=PVSF(WB)
      IF (DB-WB)1,4,7
    1 WRITE(*,10)DB,WB
   10 FORMAT(5X,'ERRORE IN PSDBWB : DB=',F6.2,'  WB=',F6.2,'   WB maggi
     *ore di DB')
      STOP
    7 HDB0=1.006*DB
      WB1=WBFF(0.0,DB)
      IF (WB .GE. WB1) GO TO 5
      WRITE(*,11)DB,WB
   11 FORMAT(5X,'ERRORE IN PSDBWB : DB=',F6.2,'  WB=',F6.2,'   temperatu
     *re non compatibili')
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
   10 FORMAT(5X,'ERRORE IN PSDBDP : DB=',F6.2,'  DP=',F6.2,'   DP maggio
     *re di DB')
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
   11 FORMAT(5X,'ERRORE IN PSDBRH : DB=',F6.2,'  RH=',F6.2,'   umidita
     *relativa non compatibile')
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
   10 FORMAT(5X,'ERRORE IN PSDPH : DP=',F6.2,'  H=',F6.3,'   dati non c
     *ompatibili')
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
   11 FORMAT(5X,'ERRORE IN PSDBW : DB=',F6.2,' W=',E10.3,' umidita spe
     *cifica non compatibile')
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
      
c***********************************************************************
c
c helper functions
c
c DPF: compute dew point temperature from partial vapor pressure

      FUNCTION DPF(PV)
      PU=PV/3386.389
      Y=ALOG(PU)
      IF(PU .GT. 0.18036)GO TO 1
      DPF=22.21111+13.81833*Y+0.49594*Y*Y
      GO TO 2
    1 DPF=26.13722+16.98833*Y+1.04961*Y*Y
    2 RETURN
      END
      
c PVSF: compute saturation partial vapor pressure from dry air temperature

      FUNCTION PVSF(TEMP)
      PRIF=101325.
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
      
c WBFF: compute wet bulb temperature from specific humidity and dry air temp.

      FUNCTION WBFF(W,DB)
      WB1=DB
	icount=0
   15 WS1=XSAT(WB1)
      W1=(WS1*(2501-2.364*WB1)+1.006*(WB1-DB))/(2501+1.83*DB-4.194*WB1)
      Y1=W-W1
	!write(6,*) 'WBFF : ',WB1,WS1,Y1,ABS(Y1)-0.00003
      IF(ABS(Y1)-0.00003)11,11,16
   16 IF(Y1)9,11,14
   14 WB1=WB1+0.5
      GO TO 15
    9 WB2=WB1-1
	!write(6,*) 'WBFF9 : ',WB1,WS1,Y1,ABS(Y1)-0.00003
      WS2=XSAT(WB2)
      W2=(WS2*(2501-2.364*WB2)+1.006*(WB2-DB))/(2501+1.83*DB-4.194*WB2)
      Y2=W-W2
	!write(6,*) 'WBFFy : ',y2,ABS(Y2)-0.00003,Y1*y2
	icount=icount+1
	if(icount .gt. 50 ) stop 'error stop WBFF: count exceeded'
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
      
c WBFF: compute saturation value (?) from temperature

      FUNCTION XSAT(TEMP)
      DIMENSION C(8)
      DATA C/0.15036368E-02,0.87755907E-04,0.17497041E-04,-0.11118891E-5
     *,0.56661491E-07,-0.13317505E-08,0.16671702E-10,-0.80607734E-13/
      TEMP=TEMP+11.
      T=TEMP+11.
      SUM=0.0
      DO 10 II=1,8
      I=9-II
      SUM=C(I)+SUM*TEMP
   10 CONTINUE
      XSAT=SUM
      TEMP=TEMP-11
      RETURN
      END

******************************************************************************

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

c-----------------------------------------------------------------------------

      subroutine srad(ia,im,id,ah,qrt,qrdif,qstot,cosi,qsmax,cc)

c-----------------------------------------------------------------------------
c
c  input:   ia, im, id, ah, qrt, qrdif
c  output:  qstot, qsmax, cc
c
c    ia,im,id,ah = year, month, day, hour
c    qrt = total solar radiation flux [W/m**2]
c    qrdif = diffuse solar radiation flux [W/m**2]
c    qstot = net solar radiation flux (reflection already subtracted) [W/m**2]
c    qsmax = maximum solar radiation flux, theoretically calculated [W/m**2]
c    cc = cloud cover
c
c    cosi = cos angle of incidence
c
c cloud cover cc: it is computed only during the day, when cosi>0
c                 during the night its value is the same as the last 
c                 hour (with cosi > 0) of that day.

	real ccold
	save ccold
	data ccold /1./

        qidir = qrt - qrdif

        call decl(ia,im,id,ad,rr)
        call incid(ah,ad,ai,cosi)

        if(cosi.ge.0.) then			!day
          call rifl(ai,ro)
          call totq(qidir,qrdif,ro,qstot)
          call radmax(rr,cosi,qsmax)
          call cloud(qstot,qsmax,cc)
	  ccold = cc
        else					!night
          qstot=0.
          qsmax=0.
	  cc = ccold
        endif

        end

c-----------------------------------------------------------------------------

      subroutine decl(iy,im,id,ad,rr)

c input: id,im,iy (day,month,year)
c output:
c  ad = declinazione solare [rad]
c  rr = ??
c
c 1 Jan ==> nd = 1

        integer idmon

      pi=4*atan(1.)

      nd=id+jdmon(iy,im-1)

      ad=(23.45*pi/180.)*sin(2.*pi*(284+nd)/365.)
      rr = 1. + 0.033*cos(2.*pi*nd/365.)

      end

c-----------------------------------------------------------------------------

      subroutine incid(ah,ad,ai,cosi)

c input:
c  ah = hour
c  ad = declinazione solare
c output:
c  ai = angolo di incidenza della luce solare (in rad rispetto allo zenith)
c  cosi = cos(ai)

      pi=4*atan(1.)
      omega=2*pi/24.
      alat=45.5
      alatr=alat*pi/180.
      t12=ah-12.

      cosi=sin(ad)*sin(alatr)+cos(ad)*cos(alatr)*cos(omega*t12)
      ai=acos(cosi)
     
      end

c-----------------------------------------------------------------------------

      subroutine rifl(ai,ro)

c input:  ai = angolo di incidenza della luce solare [rad]
c output: ro = percentage of reflected radiation [0-1]

      an2=1.33

      aj = asin(sin(ai)/an2)
      ap = aj + ai
      am = aj - ai

      ro = 0.5*( (sin(am)/sin(ap))**2 + (tan(am)/tan(ap))**2 ) 

      end

c-----------------------------------------------------------------------------

      subroutine totq(qidir,qidif,ro,qitot)

      rodif=0.0589
      qitot = (1-ro)*qidir + (1-rodif)*qidif

      end

c-----------------------------------------------------------------------------

         subroutine radmax(rr,cosi,qsmax)

c input:
c   rr = ??
c   cosi = cos(ai)
c output:
c   qsmax = theoretical maximum radiation [W/m**2]
c
c  Lazzarin R., 1981: Sistemi solari attivi
c  sc = costante solare [W/m**2]

      sc = 1353.
      qsmax = rr*sc*cosi

      end

c-----------------------------------------------------------------------------

         subroutine cloud(qstot,qsmax,cc)

         cctmp = 1. - qstot/qsmax
         if (cctmp.ge.0.) cc=cctmp

      end

c-----------------------------------------------------------------------------

        function jdmon(year,month)
 
        implicit none
 
        integer jdmon
        integer year,month
 
        logical bises
 
        integer jmm(0:12)
        save jmm
        data jmm /0,31,59,90,120,151,181,212,243,273,304,334,365/
 
        jdmon=jmm(month)
 
        if( month .ge. 2 .and. bises(year) ) then
          jdmon = jdmon + 1
        end if
 
        end                                                                     

c-----------------------------------------------------------------------------

        function bises(year)
 
c true if year is bisestial
 
        logical bises
        integer year
 
        if(  ( mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0 )
     +                  .or. mod(year,400) .eq. 0 ) then
                bises = .true.
        else
                bises = .false.
        end if
 
        end                                                                     

c-----------------------------------------------------------------------------

