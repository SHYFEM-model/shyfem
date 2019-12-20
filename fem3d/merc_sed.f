
         subroutine mercury_sed_react(dtday,
     +                          k,temp,area,C,Cw,
     +             Shgsil, Shgpom, Smhgsil, Smhgpom,
     +             fdiss1w,fdiss2w,fdoc1w,fdoc2w,
     +             silt,pom,vr,bvels,bvelp)

c	
************************************************************************
* 
* !! reaction mass yied coefficients: 1.07 (met); 0.93 (demet) !!
* !! oppure calcolare in mol e poi convertire (?) !!

c -----------------------------------------------------------------------
c ---- Call from bentic sediment model: ---------------------------------
c 	 Vr, Bvels, Bvelp
c      silt, POM, por    
c-----------------------------------------------------------------------
c
c INITIAL CONDITIONS: concentrations of hg and mehg [ mg(hg)/kg(sed) ]
c that are converted to hgit and mehgt [ ug(hg)/m3(s+w) ] 
c
************************************************************************
	implicit none
c	include 'mercury.h'
	integer m,nvmerc !water or sediment box
	parameter (nvmerc=2)

      real dtday, time, col
	real C(nvmerc)		!mercury variables C(1):hgit, C(2):mehgt
	real Cw(3)		!mercury in water var : Hg0, Hg2, MeHg
	real CD(nvmerc)		!derivatives CD(1), CD(2)
	real CDw(3)		!mercury in water var : Hg0, Hg2, MeHg
	real temp  ! for t corrected diffusion coefficient
	real depth, area, vol	! box property
        real dt    !time step, day
        real cold(nvmerc)
        real coldw(3)
        integer k

c------------------------------------------------------------------
c----- variables - initial COND. - ug(hg)/g(sed)  
c	------------------------------------------------------------------
	real Hg2sed, MeHgsed	   ![ug(hg)/g(sed)]
	real hgit, hgp, hgd,hgp1, hgp2     ! [ug m-3] or [ng(hg)/l(w+s)] 
	real mehgt, mehgd_ngl, mehgp1, mehgp2, mehgpi
      real fdiss1w, fdoc1w
      real  fdiss2w, fdoc2w
c
	! hgp = Hg2sed*Csilt/10**3 !   
	! kd/por =  Hg2sed/ hgd 
	! hgit = hgp + hgd	        [ug m-3]or [ng(hg) l(w+s)-1] 
c ----------------------------------------------------------------------
c     ! Hg2x  = hgit   * fx1
c     ! MeHgx = mehgt * fx1     [ug m-3]or [ng(hg) l(w+s)-1] 
c ----------------------------------------------------------------------
	real Hg2Cl, Hg2DOC, Hgpsilt, Hg2POM, Hg2P      ! [ug m-3]or [ng(hg) l(w+s)-1] 
	real MeHgCl, MeHgDOC,MeHgPOM, MeHgsilt, MeHgP  ! [ug m-3]or [ng(hg) l(w+s)-1] 
	real Hg2D, MeHgD ![ng l(w+s)-1]  
	real Hg2dpw, MeHgdpw, Hg2pw, MeHgpw ![ng l(w)-1] hg species in pore-water


        real hgitw, mehgtw,HgDw,MeHgDw

c 	Hg2dpw, MeHgdpw sono usate x calcolo fllussi diff; Hg2pw, MeHgpw calcolate dopo integrazione
c	real Hg0sed    ! ? ---- reductive demethylation produces hg0 ---- 
c
c----- partition coefficients and fractions *****************************
c	------------------------------------------------------------------
      real K1silt, K1doc,K1POM,K1tp, num1  !Kd of Hg2 to silt,doc !K1sand, 
      real K2silt, K2doc, K2POM, K2tp,num2 !Kd of MeHg to silt, sand, doc, poc
      real partden1, partden2	!auxiliar variable in computation
      real d1,si1,sa1,P1, d2,si2,sa2,P2
      real faq1, fsilt1, fsand1, fPOM1, fdoc1, ftot1 ! fraction of hg sp. in sed
      real faq2, fsilt2, fsand2, fPOM2, fdoc2, ftot2 ! fraction of mehg sp. in sed
c --------------------------------------------------------------------
c----- Hg and MeHg resuspension and burial -----------------------------
c   --------------------------------------------------------------------
      real silt, POM, DOC                 ! [g  m-3]
	real Cdoc, Csilt,  Cpom, Ctot       ! [g cm-3] or [kg l-1]
	real Vr, Rhgsil, Rhgpom, Rmhgsil, Rmhgpom !resuspension vel [m s-1] and fluxes [..]         
	real Bvelp, Bvels, bursHg, burpHg, bursMHg, burpMHg !burial vel [m s-1] and fluxes
c ----------------------------------------------------------------------
c----- Pore-water diffusion --------------------------------------------
c   --------------------------------------------------------------------	
	real pw_m3,pw_L, por, tor !pore-w volume,liters, porosity and tortuosity 
	real hgd_ugl, hgd_ngl     !conversions for diss hg in pore-water
c	
 	real Dchg_t, Dchg_25 !diffusion coefficient [m2 s-1]
	real ft1, JHgD,   num3, num4, Jngm2d !difffusion flux
      real ft2, JMHgD,  num5, JMngm2d !difffusion flux 	
c --------------------------------------------------------------------
c----- Hg methylation and MeHg demethylation ---------------------------
c   --------------------------------------------------------------------    
      real ksmet, ksdem, Qbac
      real Eadem, tkel, tkref, R, Rcal, deltat
      real cksmet,sksme, cksdem, sksdem
c
c       deposition terms
c
        real Shgsil,Shgpom,Smhgpom,Smhgsil

c	------------------------------------------------------
c------  BOX features --------------------------------------------------
c	------------------------------------------------------
	depth = 0.02             ! m
	vol= depth*area 	       ! m3
      
        DOC = 15.                  ! [mg/l] di DOC in sediment
        por = 0.7
c	-----------------------------------------------------------------------------
c	---------REACTIONS--rate constants-------------------------------------------

      ksmet = 0.09      ! d-1
      ksdem = 0.12       ! d-1
      Qbac = 1.5	

c	------------------------------------------------------	
c -- partition coefficients [L/kg] for mercury into silt,sand,DOC,POC**
	K1silt= 10.**(5.84) 	!Hg2 in silt   
        K1doc=  1.0*10.**4.	!Hg2 in doc
	K1POM=  10.**(6.72) 	!Kd Hg2 to POM particles 
	 
	K2silt= 1.0*10.**4.	!MeHg in silt  
 	K2doc=  2.0*10.**3.	!MeHg  in doc
	K2POM=  10.**6          !MeHg in POC particles 


c ---------------Assign old variables-----------------------

        hgitw=Cw(2)
        mehgtw=Cw(3)
        hgit =C(1)
        mehgt=C(2)

c -----------------------------------------------------------
c ----- Lamborg et al., 2016 --------------------------------
c ------- log(kdPOM) = 6.72+-2%; kdCacO3 = log(6.71)+-4% ----
c -----------------------------------------------------------
c --------log(kd-lithogenic part) =5.84+-9% -----------------
c -----------------------------------------------------------
c --- Soerensen Baltic sea: one Kd proportional to LOI ------
c ---- LOI = 4.2854*p_OC+0.859-------------------------------
c -----kd hg = 2.97 + 0.15*LOI ---(= 10**3 -10**4)----------- 
c -----kd mehg = 1.98 + 0.18*LOI ---(= 10**2.2 -10**3.1)-----

***** call from solids water module********************************
c	
c	Sw, POMw, Vds, Vdp   ! solids water and sink velocities
c--------------------------------------------------------------
****** call from Hg water module***********************************
c 
c     hgitw, MeHgw ! hgit in water [ng l-1] or [μg m-3]    
c     fsilt1w,fpom1w, fdiss1w  !fraction Hg-silt -pom -diss water
c      fsilt2w, fpom2w, fdiss2w  ! fraction MeHg-silt -pom -diss water
c
****** call from solids module***********************************
****** CONCENTRATIONS of silt, POM, DOC in sed*********part 04**********
******* converted from [g m-3] to [g cm-3] or [kg l-1]***********
 
      Cdoc  = por*DOC/10.**6. ! g(DOC)/m3(w) * 0.8 m3(w)/m3(s+w) -> mg(DOC)/L(s+w)/10**6 -> kgl-1
      Csilt = silt/10.**6.  ! from [g m-3] to [g cm-3] or [kg l-1]      
      Cpom  = POM /10.**6.  ! [kg l-1]   
      Ctot  = Csilt+Cpom      
 
	num1 = (K1silt *Csilt) + (K1POM*Cpom) !+(K1DOC*Cdoc) ! L/kg* kg/L
      num2 = (K2silt *Csilt) + (K2POM*Cpom) !+(K2DOC*Cdoc)	
	K1tp = num1/Ctot!  KD to all solids
	K2tp = num2/Ctot
 	
	pw_m3 = por*vol    ! [m-3(w)] 
	pw_L  = pw_m3*1000.               ! g m-3 * m3 -> g of pore water
 
c------------------------------------------------------------------------
c ---------- Hg and MeHg Fractions ------------- ------------------------
c------------------------------------------------------part 06-----------
c ---------- kD  [g(hg)/kg(s)] / [g(hg)/l(w)] --> [l(w)/kg(s)] ----------
c ---------- por [l(w)/l(s+w)] ------------------------------------------ 
c ---------- kD/por = [l(w)/kg(s)]*[l(s+w)/l(w)]*Csolid [kg(s) l-1(s+w)]- 
c------------------------------------------------------------------------

      d1  = (K1doc/por) * Cdoc  !   
      si1 = (K1silt/por)* Csilt
      P1  = (K1POM/por) * Cpom   ! kd'*m
           
      d2  = (K2doc/por)  * Cdoc    
      si2 = (K2silt/por) * Csilt   
      P2  = (K2POM/por)  * Cpom    
      
      partden1 = 1. + si1 + P1 + d1  ! partden1=1.0+(d1+si1+sa1+P1)   
      partden2 = 1. + si2 + P2 + d2  ! partden2 =1+(si1+sa2+ P2 + d2)	
 
         faq1= 1. /partden1     !fraction of freely dissolved  Hg2	
      faq2 = 1./partden2	!fraction of freely dissolved  MeHg
 	
      fsilt1 = si1/partden1	!(0.000001*k1silt*silt)/partden1 !silt sorbed Hg2
      fsilt2 = si2/partden2 !silt sorbed MeHg
 	
	fPOM1 = P1/partden1  !(0.000001*K1POM*POC)/partden1	!POC matter sorbed Hg2
 	fPOM2 = P2/partden2	!POC matter sorbed MeHg
 	
 	fdoc1 = d1/partden1
	fdoc2 = d2/partden2	!doc sorbed MeHg
 	
	ftot1 = faq1 + fsilt1 + fPOM1 +fdoc1    !fsand1	
	ftot2 = faq2 + fsilt2 + fPOM2 +fdoc2     !+fsand2+
 
c     write(*,*) '    ' 
      write(*,*) ':::::::::::::: PARTITIONING :::::::::::::::::::::::::' 
c      write(*,*) 'K1doc :',K1doc , 'K2doc :', K2doc	!---------- 
c      write(*,*) 'K1silt', K1silt, 'K2silt:', K2silt	!------------
c      write(*,*) 'K1POM',  K1POM,  'K2POM:', K2POM	!------------
c      write(*,*) '    '           
c      write(*,*) 'por', por, 'Doc', DOC	!-----------------------
c      write(*,*) 'partden1', partden1, 'partden2', partden2	   
c      write(*,*) 'fsilt1', fsilt1, 'fsilt2', fsilt2	!------------- 
      write(*,*) 'CDOC', Cdoc, 'Cpom', Cpom	!------------------- 
      write(*,*) 'fPOM1', fPOM1, 'fPOM2', fPOM2	!------------------- 
      write(*,*) 'fdoc1', fdoc1, 'fdoc2', fdoc2	!-------------------       
c      write(*,*) 'ftot1', ftot1, 'ftot2', ftot2	!----------------                
      write(*,*) ':::::::::::::::::::::::::::::::::::::::::::::::::::::'      
      write(*,*) ' '       
c------------------------------------------------------------------------
c------------------------------------------------------------------------
c ---------- Benthic dynamics   -----------------------------------------   
c------------------------------------------------------part 07-----------
c ----------Hg MeHg RESUSPENSION  --------------------------------------- 
c        
      Rhgsil = Vr * hgit * fsilt1 *area  !m s-1 * ug m-3 * m2-> ug s-1  
      Rhgpom = Vr * hgit * fPOM1 *area   !m s-1 * ug m-3 * m2-> ug s-1  
c        
      Rmhgsil = Vr * mehgt * fsilt2 *area  ! m s-1 * ug m-3 -> ug s-1  
      Rmhgpom = Vr * mehgt * fPOM2 *area   ! m s-1 * ug m-3 -> ug s-1   
c 
      write(*,*) '    ' 
      write(*,*) ':::::::::::::: RESUSPENSION :::::::::::::::::::::::::'                       
      write(*,*)  'res [ug s-1] hg-silt:', Rhgsil,'hg-pom:', Rhgpom   
      write(*,*)  'res [ug s-1] mhg-silt:', Rmhgsil,'mhg-pom:', Rmhgpom  
      write(*,*) '    '  
      write(*,*)  'Vr:',Vr,'hgit:',hgit,'fsilt1:', fsilt1       
      write(*,*) ':::::::::::::::::::::::::::::::::::::::::::::::::::::'      
      write(*,*) ' '  
c------------------------------------------------------------------------    
c --------- Hg MeHg BURIAL  ---------------------------------------------
c------------------------------------------------------------------------
      bursHg  = (hgit*fsilt1 *Bvels*area) ![ug s-1]=[ug m-3]*[m s-1]*[m2]
      burpHg  = (hgit*fpom1  *Bvelp*area)  
       
      bursMHg = (mehgt*fsilt2*Bvels*area) ![ug s-1]=[ug m-3]*[m s-1]*[m2]
      burpMHg = (mehgt*fpom2*Bvelp*area)          
        
      write(*,*) '    ' 
      write(*,*) '::::::::::::::: BURIAL ::::::::::::::::::::::::'                       
      write(*,*)  'bur [ug s-1] hg-silt:', bursHg, 'hg-pom',burpHg
      write(*,*)  'bur [ug s-1] mhg-silt:', bursMHg, 'mhg-pom',burpMHg
      write(*,*) '    '  
      write(*,*)  'Bvels:', Bvels,'fsilt1:', fsilt1,'fpom1:', fpom1              
      write(*,*) '::::::::::::::::::::::::::::::::::::::::::::::'      
      write(*,*) '    ' 
c-----------------------------------------------------------------------
c --------- Hg MeHg DIFFUSION ------------------------------------------
c-------------------------------------------------------part 08---------
c---------- leggere Oxygen conc from biogeochem. model ? ---------------
c --Soerensen et al., 2016 --- if Oxy<0 :  Dchg_25_ANOX = 10.0/10**10 --
c ----------------------------------------------------------------------
c ---Sunderland et al., 2010 -------------------------------------------
c ---------- Dchg_25 = 2.0/10.**10.- Hg aggregated with macromolecular --
c -------------------------------------- colloidal organic matter -------
c-------- Dchg_25 = 9.5/10.**10.- HgCl4----------------------------------
c ----------------------------------------------------------------------- 
c ---------- Dchg_25 = 1.2/10.**9. -- MeHg sulfides (CH3HgSH0)----------- 
c -----------------------------------------------------------------------
      HgDw  = hgitw*(fdiss1w+fdoc1w)  ! Hg diss acque [ng l-1(w+s)] - da model water
      MeHgDw  = mehgtw*(fdiss2w+fdoc2w)  ! Hg diss acque [ng l-1(w+s)] - da model water
      Dchg_25 = 2.0/(10.**10.)                ! [m2 s-1]
      Dchg_t = Dchg_25/(1.+0.048*(25.-temp))  ! temperature dependent difffusion coefficient
      tor = 1. - log(por**2.)
     
      Hg2dpw  =  hgit  * (fdoc1+faq1)/por    !Soerensen et al., 2016  
      MeHgdpw =  mehgt * (fdoc1+faq1)/por

      num3 = (por*Dchg_t)/tor        
      num4 = (Hg2dpw - HgDw)/(depth/2.)  ![ug m-3] * [m-1] -> [ug m-4]  Hg2dpw = Hg2D/por?     
      num5 = (MeHgdpw - MeHgDw)/(depth/2.)  ![ug m-3] * [m-1] -> [ug m-4]  Hg2dpw = Hg2D/por? 

      ft1 =  - num3*num4    ![ug m-4] * [m2 s-1] -> [ug m-2 s-1]
      ft2 =  - num3*num5     
    
      JHgD     = ft1*area  ![ug s-1] < 0: da sedimenti ad acque
      Jngm2d     = ft1*1000.*86400. ![ng m-2] < 0: da sedimenti ad acque

    ! calculated:  7 - 570       [ng m-2 d-1] Emili et al. 2012
    ! measured  : 2.000 - 70.000 [ng m-2 d-1] Emili et al. 2012
      
      JMHgD     = ft2*area  ![ug s-1] < 0: da sedimenti ad acque
      JMngm2d   = ft2*1000.*86400. ![ug m-2] < 0: da sedimenti ad acque
c     
c      JHgD_kgy     = JHgD *365/10**9! kg y-1
      write(*,*) '    ' 
      write(*,*) '::::::::::::::: DIFFUSION FLUX ::::::::::::::::::::::'         
      write(*,*) 'Jngm2d:', Jngm2d,num3,num4,Hg2dpw,MeHgdpw
      write(*,*) 'hgit,mehgt', hgit,mehgt,fdoc1,faq1
      write(*,*) 'hgdw,depth', hgdw,depth,Hg2dpw
      write(*,*) '7 - 570       [ng m-2 d-1] Emili et al. 2012'
      if (JHgD .LT. 0.0) then 
      write(*,*) 'diff flux from sediment to water'
      else if (JHgD .GT. 0.0) then 
      write(*,*) 'diff flux from water to sediment'
      else 
      write(*,*) 'diff flux = 0'
      end if
      write(*,*) 'JHgD', JHgD   
      write(*,*) '::::::::::::::: DIFFUSION FLUX ::::::::::::::::::::::'               
      write(*,*) '    ' 
      write(*,*) 'hgit' , hgit, 'mehgt', mehgt

c-----------------------------------------------------------------------
c --------- Hg methylation and MeHg demethylation ---------------------- 
c -----------------------------------------------------------------------
      cksmet=(ksmet*Qbac**((temp-20.)/10.))/86400.       ! [d-1] to [s-1]
      sksme=cksmet*(hgit*(faq1+fdoc1))*vol    ! [s-1] * [ug m-3] = ug s-1

c -----------------------------------------------------------------------
c --------- Soerensen et al. 2016 FOR WATER COLUMN ----------------------
c -- OXIC COND.----- kmet = rmr/100 ! rmr = remineralization rate--------
c -- ANOXIC COND.--- kmet = (PO4 - PO4subox)*0.0005----------------------
c -----------------------------------------------------------------------
c ---met - demet in sedimentfrom observations and ratio MeHg/Hg----------
c -----------------------------------------------------------------------

	Eadem=10.	      !WASP manual around 10 kcal/mol
      tkel=temp+273.
      tkref=293.  	!reference temperature, K
   	R=8.314472        ![J K-1 mol-1] universal gas constant
	Rcal=R/4.184
      deltat = tkel-tkref

      cksdem=(ksdem*exp(Eadem*1000*(deltat/(Rcal*tkel*tkref))))/86400
      sksdem=(cksdem*mehgt*(faq2+fdoc2))*vol ! [s-1] * [ug m-3] = ug s-1

c -------------------------------------------------------------------------------- 
c --------- part 09--------------------------------------------------------------- 
c ------------------------integration--------------------------------------------- 
	
c       Check Ginevra, 
 	C(1) = hgit              ! ug m-3
 	C(2) = mehgt	
 	
      if (bursHg .GE. 0.) then       
      CD(1) = +Shgsil +Shgpom -Rhgsil -Rhgpom +JHgD -bursHg -burpHg 
     & +sksdem-sksme
      else 
      CD(1) = +Shgsil +Shgpom -Rhgsil +JHgD -Rhgpom +sksdem-sksme
      end if
      if (bursMHg .GE. 0.) then       
      CD(2) = +Smhgsil+Smhgpom-Rmhgsil-Rmhgpom+JMHgD-bursMHg-burpMHg
     & -sksdem+sksme
      else 
      CD(2) = +Smhgsil+Smhgpom-Rmhgsil-Rmhgpom+JMHgD -sksdem+sksme
      end if    

        CDw=0
        CDw(2)=Rhgsil+Rhgpom            !Check Ginevra
        CDw(3)=Rmhgsil+Rhgpom           !Check Ginevra
 
	call merc_euler (nvmerc,dt,vol,vol,c,cold,cd) !claurent-OGS:add second volume
	call merc_euler (3,dt,vol,vol,cw,coldw,cdw) !claurent-OGS:add second volume
c	write(*,*) 'Hg2sed' , Hg2sed, 'MeHgsed', MeHgsed
	
      write(*,*) 'hgit' , hgit, 'mehgt', mehgt
     
!	------------------------------------------------------------------
	Hgpsilt= hgit * fsilt1  ! [ug m-3] or [ng(hg) l(w+s)-1]
	Hg2POM = hgit * fPOM1   
	Hg2doc = hgit * fdoc1
      Hg2Cl =  hgit * faq1    ! hgD in the tot volume (s+w)

      write (*,*) '' 
      write (*,*) 'Hgpsilt',Hgpsilt, 'Hg2POM',Hg2POM
      write (*,*) ''
      write (*,*) 'Hg2doc',Hg2doc, 'Hg2Cl',Hg2Cl 
      write (*,*) ''
      write (*,*) 'hgit', hgit

        
	MeHgsilt= mehgt * fsilt2  ![ug(hg) m-3(w+s)] or [ng l-1]
   	MeHgPOM = mehgt * fPOM2
   	MeHgdoc = mehgt * fdoc2
      MeHgCl  = mehgt * faq2  ! mehgD in the tot volume (s+w)   

	Hg2P  = (Hgpsilt + Hg2POM)       ![ug m-3(w+s)] or [ng l-1]  
	MeHgP = (MeHgsilt + MeHgPOM)        

      Hg2D = (Hg2Cl + Hg2doc) 
      Hg2pw = Hg2D/por
     
  	MeHgD = (MeHgCl + MeHgdoc) ! [ug m-3(w+s)] or [ng l-1] 
  	MeHgpw  = MeHgD/por

c     add1 = Hgpsilt/silt -- NO, viene molto diverso
c     add2 = Hg2POM/POM

      Hg2sed = Hg2P/(silt+POM)
      MeHgsed = MeHgP/(silt+POM)

c     Cdoc  = por*DOC/10.**6. ! g(DOC)/m3(w) * 0.8 m3(w)/m3(s+w) -> mg(DOC)/L(s+w)/10**6 -> kgl-1
c     Csilt = silt/10.**6.  ! from [g m-3] to [g cm-3] or [kg l-1]      
c     Cpom  = POM /10.**6.  ! [kg l-1]   
c     Ctot  = Csilt+Cpom 

      write (*,*) '' 
      write (*,*) 'Hg2P   ', Hg2P
      write (*,*) ''
      write (*,*) 'Hg2D   ', Hg2D
      write (*,*) ''  
 	
      write (88,*) (C(m), m=1,nvmerc)
      write (89,*) Hg2sed, Hg2pw 
      write (91,*) Hgpsilt, Hg2POM,Hg2doc,Hg2Cl 
      write (92,*) Shgsil, Shgpom,Rhgsil, Rhgpom,JHgD, bursHg,burpHg
      write (93,*) Smhgsil,Smhgpom,Rmhgsil,Rmhgpom,JMHgD,bursMHg,burpMHg
      write (94,*) MeHgsed, MeHgpw 
      write (95,*) MeHgsilt, MeHgPOM,MeHgdoc,MeHgCl 
      write (96,*) sksme, sksdem

c------ dissolved hg and mehg in porewater ------------------------------
c     hgd = (faq/por)*hgit    !Thomann and Di Toro, 1983
c
      write(*,*) 'Hg2sed', Hg2sed ,'MeHgsed', MeHgsed
      write (*,*) ''

c     stop
c	Hg0sed=C(1)
c     Hg2sed=C(1)+DepHg2-burHg
c     MeHgsed=C(3)+ DepMeHg-burMeHg

c     write (89,*) time, i,(C(m), m=1,nvmerc) !, ' variables'
c     write (91,*) DepHg2,DepMeHg,burHg,burMeHg,dt !, ' variables'        
c --------------------------------------------------------
c--------------------------------------------------
      end 
c	
