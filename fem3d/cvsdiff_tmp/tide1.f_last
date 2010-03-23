
c****************************************************************

        subroutine tidesc(iy,im,id,ih,imn,isec)

ccc scenario astronomico ad un certo tempo, dato anno (quattro cifre),
ccc  mese, giorno, ora, minuto.
ccc Il tempo si riferisce al Europa centrale, facilm. riportabile a UT.
ccc   Nel common /ASTRO/ sono registrati i parametri importanti:
ccc    coordinate nella volta celeste di sole e luna, nonche` le loro
ccc     distanze dal centro della terra.

        implicit real*8(a-h,o-z)
        common /astro/cosz, sinz, cosps, sinps, rr,
     +                coszm,sinzm,cospsm,sinpsm,rrl
       DIMENSION IMONTH(12)
       DATA EC,SM,SIN23,COS23,SINI,COSI,C23CI,S23SI/.0549,.074804,.39798
     $06547d0,.9173938077d0,.0896765581d0,.9959709407d0,.9136975736d0,
     $.0356895353d0/
      DATA  IMONTH/0,31,59,90,120,151,181,212,243,273,304,334/
       DATA INIT/0/
       DAYT=IH+IMN/60.d0+isec/3600.d0
 133       YEART=DAYT+24*(IMONTH(IM)+ID-1)
 2       ST=YEART+((IY-1900)*365+FLOAT((IY-1900)/4))*24  -1
       IF(((IY/4)*4 .EQ.IY).AND.(IY.NE.1900).AND.(IM.LT.3))ST=ST-24.d0
          T=(ST +12.)/876600.
       H=4.8816280d0+628.3319509d0*T+52.d-7*(T**2)
       P=4.9082295d0+0.0300053d0*T+79.d-7*(T**2)
       E=0.01675104d0-4.18d-5*T-1.26d-7*(T**2)
       RS=1.+E*DCOS(H-P)+(E**2)*DCOS(2.*(H-P))
       SL=H+2.*E*DSIN(H-P)+1.25*(E**2)*DSIN(2.*(H-P))
       COSZS=DSIN(SL)*SIN23
       Z=dacos(COSZS)
       AA=.5*(SL+Z-1.57079)
       TPS=DTAN(AA)*1.52386101d0
       PS1=2.*dATAN(TPS)
       PSG=DMOD(.2617993879d0*(ST-12.)+H,6.2831853071796d0)
       PSS=PS1-PSG
       HS=H
       ES=E
       PSOL=P
       H=4.7200089d0+8399.7092745d0*T+0.0000346d0*(T**2)
       P=5.8351526d0+71.0180412d0*T-0.0001801d0*(T**2)
       SN=4.5236016d0-33.7571463d0*T+0.0000363d0*(T**2)
         SSN=DSIN(SN)
       CSN=DCOS(SN)
      COSOC=C23CI-S23SI*CSN
       OC=dacos(COSOC)
       SINOC=dSQRT(1.-COSOC**2)
       SINNU=SINI*SSN/SINOC
         ANU=dasin(SINNU)
        COSNU=dSQRT(1.-SINNU**2)
       SINAO=SSN     *SIN23/SINOC
       COSAO=CSN     *COSNU+COS23*SSN     *SINNU
       TAN2AO=SINAO/(1.+COSAO)
        AO=2.*dATAN(TAN2AO)
        H=DMOD(H,6.28318530717796d0)
        P=DMOD(P,6.28318530717796d0)
      SN=DMOD(SN,6.28318530717796d0)
       HMP=H-P
        SIG=H-SN+AO
C    $*SM))*DSIN(2.*(H-HS))-3.*SM*ES*DSIN(HS-PSOL)+2.125*(SM**2)*EC*DSIN
       SL=SIG+2.*EC*dSIN(HMP)+1.25*(EC**2)*dSIN(2.*HMP)+SM*EC*(3.75+263.
     $*SM/16.)*DSIN(H-2.*HS+P)+SM**2*(11./8.+59.*SM/12.+75.*(EC**2)/(16.
     $*SM))*DSIN(2.*(H-HS))                       +2.125*(SM**2)*EC*DSIN
     $(3.*H-2.*HS-P)+(77./16.)*(SM**2)*ES*DSIN(2.*H-3.*HS+PSOL)
       EW=EC**2
      RR=1.+(1./(1.+SM**2/6.))*(EC*dCOS(HMP)+EC**2*dCOS(2.*HMP)+SM*EC*(
     $15.
     $/8.+329.*SM/64.)*DCOS(H-2.*HS+P)+(SM**2)*(1.+19.*SM/6.+15.*EW/(4.*
     $SM))*DCOS(2.*(H-HS))-1.5*(SM**2) *ES*DCOS(HS-PSOL)+(33./16.)*(SM**
     $2)*EC*DCOS(3.*H-2.*HS-P)+(7./2.)*(SM**2)*ES*DCOS(2.*H-3.*HS+PSOL))
        rrl=rr
        rr=rs
       SW=SL
       COSZ=dSIN(SW)*SINOC
       Z=dacos(COSZ)
      TAN2P=DTAN(.5*(SL+Z-1.57079))*dSIN(.5*(1.57079+OC))/dSIN(.5*
     $(1.57079
     $-OC))
       PSI=2.*dATAN(TAN2P)
        PS=ANU+PSI-PSG
        coszm=cosz
        cosz=coszs
        psm=ps
        ps=pss
        SINZ=dSQRT(1.-COSZ**2)
       SINZM=dSQRT(1.-COSZM**2)
       SINPS=dSIN(PS)
       COSPS=dCOS(PS)
       SINPSM=dSIN(PSM)
       COSPSM=dCOS(PSM)
        return
        end

c****************************************************************

        subroutine tidelc(lon,lat,zequi)

c lon	longitudine
c lat	latitudine
c zequi	equilibrium tide in meters

c        subroutine local(coszv,sinzv,cospv,sinpv,zeq)
ccc data la situazione astronomica desunta dal common /ASTRO/,
ccc  date le coordinate della localita` di interesse, sotto forma
ccc  di coseno e seno della colatitudine (angolo dal polo) e della
ccc  longitudine (positiva a est), si ottiene la marea di equilibrio (cm).
ccc  Tutti i numeri reali sono in doppia precisione

        implicit real*8(a-h,o-z)
	real lon,lat,zequi
	parameter(pi2=6.28318530717796d0,rad=pi2/360.d0)
        common /astro/cosz, sinz, cosps, sinps, rr,
     +                coszm,sinzm,cospsm,sinpsm,rrl

	colat = rad * (90.-lat)
	coszv = dcos( colat )
	sinzv = dsin( colat )

	alon = rad * lon
	cospv = dcos( alon )
	sinpv = dsin( alon )

       COSOL=COSZ*COSZV+SINZ*SINZV              *(COSPV*COSPS  +SINPV*
     *SINPS  )
       COLUN=COSZM*COSZV+SINZV*SINZM              *(COSPV*COSPSM
     *+SINPV*SINPSM  )
       ZEQ=(16.427*(1.5*COSOL**2-.5)+7.00365E-4*RR*(1.6666*COSOL**3-1.5*
     *COSOL))*RR**3+(35.785*(1.5*COLUN**2-.5)+.59378*RRL*(1.6666*COLUN**
     .   3 -1.5*COLUN))*RRL**3

	zequi = 0.01 * zeq

        return
        end

c****************************************************************

