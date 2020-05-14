!--------------------------------------------------------------------------
!
!    Copyright (C) 2017  Ginevra Rosat
!    Copyright (C) 2017  Donata Melaku Canu
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main
!    directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------
!
! mercury dynamics in water subroutine
!
! revision log :
!
! 04.05.2017  dmc        updated mercury routine with merc_water
! 17.05.2017	dmc&gir	add volatilization subroutine gas merc_gas_exchange
!
!--------------------------------------------------------------------------
!
      subroutine mercury_react(id,bsurf,bbottom,boxtype,dtday,vol
     +                         ,depth,k,temp,uwind10,area,sal,qrad
     +                         ,C,loads,vds,vdp,conz1,conz2,    !tday,
     +                         Shgsil,Shgpom,Smhgsil,Smhgpom,
     +                         faq1,faq2,fdoc1,fdoc2)
   
       implicit none

        integer nstate
        parameter (nstate=3)

        real C(nstate)  !mercury variables C(1):Hg0, C(2):Hg2,C(3):MeHg
        real CD(nstate) !derivatives CD(1):difHg0,CD(2):difHg2,CD(3):difMeHg
        integer m       !time indicator to write outputs, debug
        integer k       !node
        integer id      !id of node and level
        logical bsurf   !is on surface?
        logical volatnew!enable volatilisation formula (Sorensen et al.,2010)
        logical bbottom !is on bottom?
        logical boxtype !is water box?

        real cold(nstate) !old state variables
        real dtday                 !time step                     [day]
        real qrad                  !irradiance                   [W/m2]
        real depth, area, vol      !dept, area, volume    [m],[m2],[m3]
        real temp,sal                    !temperature, salinity [C],[-]
        real tkel                            !temperature - T       [K]
        real tkref                           !ref temperature, 293  [K]
        real uwind10  !wind velocity from hydrodynamic module     [m/s]
        real vds,vdp  !velocity of solids deposition to sediment  [m/s]
        real conz1,conz2 !silt and pom in water from sed4merc    [mg/l]
        real flux        !Hg0 air-sea exchange flux           [ng/m2 h]
c
        real loads(nstate)      !atmospheric loadings

c	variables
        real Hg0d, Hg0a     !Hg0 in water and atmosphere [ng/l], [ug/m3]
        real Hg2, MeHg      !total HgII and MeHg in water         [ng/l]
        real Hg2d, Hg2DOC   !Dissolved phases of HgII in water    [ng/L]
        real MeHgd,MeHgDOC  !Dissolved phases of MeHg in water    [ng/L]
        real Hg2silt,Hg2sand,Hg2org,Hg2sorb!Solid phases of HgII in water [ng/L]
        real MeHgsilt,MeHgsand,MeHgorg,MeHgsorb !Solid phases of MeHg in water [ng/L]
        real Hg2dsed, MeHgsed !total HgII and MeHg in sediment
c
c	partition coefficients

      real K1silt, K1sand, K1doc, K1org  !Partition coefficients for HgII [l/kg]
      real K2silt, K2sand, K2doc, K2org  !Partition coefficients for MeHg [l/kg]
      real partden1, partden2   !auxiliar variables in computation
      real faq1,fdoc1 !ionic and DOC-complexed fraction of Hg2  [-]
      real faq2,fdoc2 !ionic and DOC-complexed fraction of MeHg [-]
      real fsilt1, fsand1, forg1 !fraction of Hg2 in silt, sand ond POM [-]
      real fsilt2, fsand2, forg2 !fraction of MeHg in silt, sand ond POM [- ]
      real ftot1, ftot2 !sum of all the fraction for Hg2 and MeHg
c
c     solids in water
      real DOC  !concentrations of dissolved organic carbon in water [mg/l]
      real silt !concentrations of silt in water [mg/l]
      real sand !concentrations of sand in water [mg/l]
      real org  !concentrations of POM in water  [mg/l]
c
c     Hg0 sink reactions, parameters
c
      real skvo        !Hg0 air-sea exchange flux            [ng/m2 h]
      real He          !Henry 's Law constant                [atm/mol]
      real R, Rcal     !universal gas constant  [J/mol K], [cal/mol K]
      real kox         !dark oxidation rate                      [1/d]
      real OxAct   !oxydation activation energy            [kcal/mol]
      real ckox, skox

c     HgII sink reactions, parameters
c
       real skme           !bacterial methylation               [ug/m3d]
       real Qbac           !T correction factor for methylation      [-]
       real kmeth, ckmeth  !bacterial meth rate, T-adjusted rate   [1/d]
       real skph           !photoreduction                      [ug/m3d]
       real kphr,ckph      !photored rate, light-adjusted rate     [1/d]
       real lref, ke  !reference light,extintion coefficient[w/m2],[1/m]
       real ladj           !normalised irradiance
       real xHg2d,xHg2DOC,xHg2sorb !activate Hg species in transformations
c
c      MeHg sink reactions, parameters
       real skdem        !bacterial demethylation               [ug/m3d]
       real kdem, ckdem  !bacterial demeth rate, T-adjusted rate   [1/d]
       real cordem       !T correction factor for demethylation      [-]
       real Eadem        !activation energy for demethylation [kcal/mol]

c      photoreductive demethylation
       real skphdem         !photoreductive demethylation       [ug/m3d]
       real kphdem, ckphdem !photodem rate,light-adjusted rate     [1/d]
       real xMeHgd,xMeHgDOC,xMeHgsorb!activate MeHg species in transformations

c      deposition fluxes and rates
       real Dhgsil, Dhgpom, Dmhgsil, Dmhgpom ! Deposition of HgP and MeHgp
       real Shgsil, Shgpom, Smhgsil, Smhgpom ! Deposition rates   

       integer ipext,ipint,kext     !nodes external and internal numbers
       integer fortfilenum, iter
       real tday
                       
c       --------------------------------------------------
c       volatilization formula
        volatnew=.true.
c	---------------------------------------------------

        silt=conz1  !FIXME      !mg/L
        org=conz2  ! FIXME! mg/L
        sand=0.!
c        silt=conz        !mg/L
        DOC=3. !mg/L (INPUT)

        Hg0a=0.0016     ![ug/m3] atmospheric Hg concentration (INPUT)

c	InHg2 = 0.001 !0.000005
c	InMeHg = 0.001 !0.0000005
c	-----------------------------------------------------
c	parameters:set and check before running a simulation
c	------------------------------------------------------
c	FIXME
c	------------------------------------------------------
c	availability of  Hg species
c	----------------------
        xHg2d=1.        !FIXME WASp uses a dimensional value expressed in [L/kg] 
        xHg2DOC=1.      !FIXME 	
        xHg2sorb=0.     !FIXME WASP:
        xMeHgd=1.
        xMeHgDOC=1.
        xMeHgsorb=0.
c	---------------------
c	partition coefficients for mercury into  silt,sand,DOC,ORG-sediment sorbed
c	part coefficients are [L/kg]
c	---------------------
        k1silt=200000.   !Hg2 in silt     [L/kg]
        k1sand=0.       !Hg2 in sand     [L/kg]
        k1doc=20000.    !Hg2 in doc      [L/kg]
        k1org=400000.    !Hg2 in organic particles [L/kg]
        k2silt=100000.   !MeHg in silt  [L/kg]
        k2sand=0.       !MeHg in sand  [L/kg]   
        k2doc=10000.    !MeHg  in doc  [L/kg]
        k2org=500000.    !MeHg in organic particles [L/kg]
c       -----------------------------------------------------
c	------------------------------------------------------
c	initial conditions: mercury species concentration in [mg/L] or  [mg/kg]  (sed)
c	------------------------------------------------------
c	----------------------------------------------------
c	global constants
c	----------------
        tkref=293       !reference temperature, K
c  	------------------------------------------------------
c       volatilization
c	--------------
        He=0.0071       ![atm/mol] Henry 's Law constant
        R=8.314472      ![J K-1 mol-1] universal gas constant
c        write(6,*) kvol
c
c       oxydation of Hg0d
c	-----------------
        OxAct=2 !activation energy for oxydation, WASP impl. FIXME =2 k[kcal mol-1]
        kox=0.00001    !0001  !FIXME oxydation parameter 
c
c	methylation
c	----------
        kmeth=0.006 ![day-1] wq  methylation rate at 20°C Monperrus et al., 2007
        Qbac = 1.5      !FIXME 
c
c       photoreduction of Hg2d to Hg0
c	----------------------------
        kphr=0.05       ![day-1] photoreduction rate constant of HgII to Hg0 (default WASP) 
        lref=240 !950                !reference light intensity for kph [watt/m2] FIXME (default WASP)
        ke=1.05             !light extintion coefficient FIXME

c	------------------------------------
c	photoreductive demethylation
c	--------------------------
        kphdem=0.00015    ![day-1] source of this parameter? FIXME
c	------------------------------------------------------------------------
c	bacterial demethylation
c	---------------------
        Eadem=2 !WASP manual around 10 kcal/mol
c ksdem lignano:0.159 s.andrea:0.093 buso:0.139 morgo:0.139 grado:0.064 primero:0.064 Hines et al. 2012
        kdem=0.15     !069	!bacterial demethylation in water	!FIXME Monperrus?

c	temp=22.	!FIXME

c	call rddepth (depth)	!depth of the element
c        depth=1.        !FIXME 1 m
c        vol=1.        !FIXME 1 m3
        
        skvo=0
        skox=0
        skph=0 
        skphdem=0
        skox=0
        skph=0
        skdem=0
        skme=0
        skdem=0
        skphdem=0

c       _______________________________________________________
c       assigne old value to mercury variables

        if(boxtype) then

       Hg0d=C(1)
       Hg2=C(2)        !+ InHg2
       MeHg=C(3)       !+ InMeHg

        else

        Hg0d=0
        Hg2=C(2)        !+ InHg2
        MeHg=C(3)       !+ InMeHg

        end if

c --------------------------------------------------------------
c
        tkel=temp+273
        Rcal=R/4.184
c	------------------------------------
c	partition of mercury spp (HG2 and MeHg) into solid phases
c	--------------------------------------

        partden1=(k1silt*silt+K1sand*sand+K1org*org+K1doc*DOC)   
        partden1=1+0.000001*partden1   
          !1+(k1doc*DOC)+((k1org*org)+(k1silt*silt)+(k1sand*sand))   
        faq1=1./partden1!fraction of freely dissolved  Hg2
        fsilt1=(0.000001*k1silt*silt)/partden1
        !(0.000001*k1silt*silt)/partden1 !silt sorbed Hg2
        fsand1=(0.000001*k1sand*sand)/partden1
        !(0.000001**k1sand*sand)/partden1!sand sorbed Hg2
        forg1=(0.000001*k1org*org)/partden1
        !(0.000001*k1org*org)/partden1	!organic matter sorbed Hg2
        fdoc1=(0.000001*k1doc*DOC)/partden1
        !(0.000001**k1doc*DOC)/partden1 !doc sorbed Hg2
        ftot1=faq1+fsilt1+fsand1+ forg1+ fdoc1

        partden2=(k2silt*silt+K2sand*sand+K2org*org+K2doc*DOC)
        partden2=1+0.000001*partden2
        faq2=1./partden2        !fraction of freely dissolved  MeHg
        fsilt2=(0.000001*k2silt*silt)/partden2 !silt sorbed MeHg
        fsand2=(0.000001*k2sand*sand)/partden2  !sand sorbed MeHg
        forg2=(0.000001*k2org*org)/partden2     !organic matter sorbed MeHg
        fdoc2=(0.000001*k2doc*DOC)/partden2     !doc sorbed MeHg
        ftot2=faq2+fsilt2+fsand2+ forg2+ fdoc2

c         write(*,*) ' '
c         write(*,*) 'ftot HgWat:',ftot1,'ftot MeHgWat:',ftot2
c         write(*,*) ' '        
  
c        write(69,*) faq2, fsilt2, fsand2, forg2, fdoc2, ftot2
c	-----------------------------------
c	calculation of mercury fractions amongst the different phases
c	------------------------------------
c       write(6,*) Hg2, 'Hg2?'
        Hg2d=Hg2*faq1
        Hg2silt=Hg2*fsilt1
        Hg2sand=Hg2*fsand1
        Hg2org=Hg2*forg1
        Hg2doc=Hg2*fdoc1
        MeHgd=MeHg*faq2
        MeHgsilt=MeHg*fsilt2
        MeHgsand=MeHg*fsand2
        MeHgorg=MeHg*forg2
        MeHgdoc=MeHg*fdoc2
        Hg2sorb=Hg2silt+Hg2sand+Hg2org
        MeHgsorb=MeHgsilt+MeHgsand+MeHgorg

c		-----------------------------
c		Hg0d--> Hg0a volatilization
c
        if(bsurf) then
c        if(volatnew)then
        call merc_gas_exchange(sal,temp,area,uwind10,Hg0d,Hg0a,skvo)
c        else

c        ckvol=kvol*(ktvol**(temp-20))   !calculated volatilization rate at temperature temp
c        write(6,*) ckvol, 'ckvol'
c        skvo=ckvol*(Hg0d-(Hg0a/(He/R*tkel))) ![g/m3/day] of Hg0 exchanged with the atmosphere.
c        end if
        else
        skvo=0
        end if

c       write(6,*) skvo,'skvo',bsurf, 'bsurf'
c			--------------------
c	Hg0d --> Hg2d oxydation
c
c	conversion of R, gas constant, from [J K-1 mol-1] to [cal k-1 mol-1]
      ckox=kox*exp(OxAct*1000*((tkel-tkref)/(Rcal*tkel*tkref))) !1000: conversion from kcal to cal
      skox=ckox*Hg0d
c       write(81,*) tkel,skvo,skox !'volatilization,oxydation'
c       ok

c	----------------------------------------
c		Hg2d--> Hg0d photoreduction 
c
        ladj=qrad/lref*((1-exp(-ke*depth))/ke*depth)	
        ckph=kphr*ladj   !FIXME unità di misura
	skph=ckph*(Hg2d*xHg2d+Hg2DOC*xHg2DOC+Hg2sorb*xHg2sorb)	!FIXME unità di misura

c	---------------------------------------------
c		Hg2d --> MeHgd methylation

	ckmeth=kmeth*Qbac**((temp-20)/10)
	skme=ckmeth*(Hg2d*xHg2d+Hg2DOC*xHg2DOC+Hg2sorb*xHg2sorb)
c	write(83,*) temp, kmeth, Qbac, ckmeth, skme !, ' kmeth Qbac ckmeth skme'
c        write(83,*) Hg2d,xHg2d,Hg2DOC,xHg2DOC,Hg2sorb,xHg2sorb,'test'
c
c	----------------------------------------------
c	MeHg --> Hg0 photoreductive demethylation in the water column
c
	ckphdem=kphdem*ladj	!FIXME insert a different ladj?
	skphdem=ckphdem*(MeHgd*xMeHgd+MeHgDOC*xMeHgDOC+MeHgsorb*xMeHgsorb)
c	write (85,*) skphdem,ckphdem !, ' photochemical demethylation'
c
c       ______________________________________________
c	MeHg --> HgII bacterial demethylation in water 
c	
c	bacterial demethylation
c	real skdem,kdem
c	real ckdem, Eadem
	
	cordem=exp(Eadem*1000*((tkel-tkref)/(Rcal*tkel*tkref)))
	ckdem=kdem*cordem
	skdem=ckdem*(MeHgd*xMeHgd+MeHgDOC*xMeHgDOC+MeHgsorb*xMeHgsorb)
c	write (86,*) tkel,cordem,ckdem,skdem	!, ' water column bacterial demethylation'
c
c       Compute Compute deposition fluxes and rates 

        Dhgsil = vds*Hg2silt   ! [m s-1]*[g m-3]*[-]=[g m2s-1]
        Dhgpom = vdp*Hg2org    ! [m s-1]*[g m-3]*[-]=[g m2s-1]
        Dmhgsil= vds*MeHgsilt  ! [m s-1]*[g m-3]*[-]=[g m2s-1]
        Dmhgpom= vdp*MeHgorg   ! [m s-1]*[g m-3]*[-]=[g m2s-1]
        Shgsil = Dhgsil * area *86400  ! [g m-2 day-1] * [m2] = [g day-1] 
        Shgpom = Dhgpom * area*86400 ! [g/day]
        Smhgsil = Dmhgsil* area*86400 ![g/day]
        Smhgpom = Dmhgpom* area*86400 ![g/day]

c       __________________________________________________

        C(1)=Hg0d
        C(2)=Hg2
        C(3)=MeHg
c      write(6,*) C(1),Hg0d,C(2),Hg2,C(3),MeHg,'merc var'
c	
c	CD= transformations 1:Hg0 2:Hg2 3:MeHg	
        CD(1) = -skvo*area+vol*(- skox + skph + skphdem) !g/day
        CD(2) = -Shgsil-Shgpom+(skox - skph + skdem - skme)*vol        !g/day
        CD(3) = -Smhgsil-Smhgpom+ (skme - skdem - skphdem)*vol    !mass, g/day

      kext=ipext(k)
c      if (kext .EQ. 1372) then 
  
c      write (487,*) vdp, vds, 'merc_water'
c      write (488,*) skvo/area
c      write (489,*) Shgsil/area     
c      write (490,*) Shgpom/area
c      write (491,*) skme/area
c      write (492,*) skdem/area
c      write (493,*) skphdem/area ! ' Transformations'
c      write (494,*) skox/area ! ' Transformations'
c      write (495,*) (C(m), m=1,nstate),kext    !' HgW vars old'
c      write (496,*) vol    !' HgW vars old'
c      write (497,*) area    !' HgW vars old'
c      write (498,*) vds
c      write (499,*) vdp
c      end if
     
c       write (444,*) (C(m), m=1,nstate),k    !' HgW vars old'
c       write(88,*) skox,skph,skphdem,skdem,skme
c       write(88,*)areaivol,'area and vol'

c      integration



       if(silt.LT.0) then
       write(445,*) 'SILTW<=0',silt,'node',ipext(k),depth,'m_watbef'
       end if

       if(org .LT.0) then
       write(448,*) 'POMW<=0',org,'node',ipext(k),depth,'m_watbef'
       end if


       if(C(1) .LT.0) then
c       write(555,*) 'hg0W<=0',C(1),'node',ipext(k),'z=',depth,'m_watbef'
       write(*,*) 'hg0W<=0',C(1),'node',ipext(k),'z=',depth,'m_watbef'
       end if
       
       if(C(2) .LT.0) then
c       write(666,*) 'hg2W<=0',C(2),'node',ipext(k),'z=',depth,'m_watbef'
       write(*,*) 'hg2W<=0',C(2),'node',ipext(k),'z=',depth,'m_watbef'
       end if

       if(C(3) .LT.0) then
       write(*,*) 'mhgw<=0',C(3),'node',ipext(k),'z=',depth,'m_watbef'
c       write(777,*) 'mhgw<=0',C(3),'node',ipext(k),'z=',depth,'m_watbef'
       end if

      call load0d_merc(dtday,cd,loads,vol)
c
c      call merc_euler (3,dt,vol,vol,c,cold,cd)  ! claurent-OGS here volold=volnew=vol
      call merc_euler (3,dtday,vol,c,cold,cd)  ! claurent-OGS here volold=volnew=vol
      
       if(C(1) .LT.0) then
       write(*,*) 'hg0W<=0',C(1),'node',ipext(k),'z=',depth,'m_wataft'
       end if

       if(C(2) .LT.0) then
       write(*,*) 'hg2W<=0',C(2),'node',ipext(k),'z=',depth,'m_wataft'
       end if

       if(C(3) .LT.0) then
       write(*,*) 'mhgw<=0',C(3),'node',ipext(k),'z=',depth,'m_wataft'
       end if


        if (C(1).le.0.) C(1)=0.00000001  !FIXME controllo mostruoso
c le trasf. di C1 la portano negativa! skvo e skox prevalgono	
c       write(*,*) 'Hgod',Hg0d,'Hg2',Hg2,'MeHg',MeHg, 'nodo',k

c        write (444,*) (C(m), m=1,nstate),'k',k,' HgW vars new' !

c     iter=nint(tday*86400.)

c       if (MOD (iter,1800) .EQ. 0) then
       if (k .GE. 1) then
           kext=ipext(k)
           fortfilenum=-1
           if(kext==3985)then
               fortfilenum=350
           elseif(kext==3986) then
               fortfilenum=351
           elseif(kext==3982) then
               fortfilenum=352
           elseif(kext==4007) then
               fortfilenum=353
           elseif(kext==3763) then
               fortfilenum=354
           elseif(kext==3765) then
               fortfilenum=355
           elseif(kext==3764)then
               fortfilenum=356
           elseif(kext==3762) then
               fortfilenum=357
           else if(kext==2150)then
               fortfilenum=358
           elseif(kext==2009) then
               fortfilenum=359
           elseif(kext==2359) then
               fortfilenum=360
           elseif(kext==2358) then
               fortfilenum=361
           elseif(kext==2341) then
               fortfilenum=362
           elseif(kext==2408) then
               fortfilenum=363
           elseif(kext==2191)then
               fortfilenum=364
           elseif(kext==2192) then
               fortfilenum=365
           elseif(kext==2654)then
               fortfilenum=366
           elseif(kext==2505) then
                fortfilenum=367
           elseif(kext==2655) then
               fortfilenum=368
          elseif(kext==2653) then
               fortfilenum=369
           else if(kext==1372)then
               fortfilenum=370
           elseif(kext==1375) then
               fortfilenum=371
           elseif(kext==1331) then
               fortfilenum=372
           elseif(kext==1378) then
               fortfilenum=373
           elseif(kext==4347) then
               fortfilenum=374
           elseif(kext==3216) then
               fortfilenum=375
           elseif(kext==3057)then
               fortfilenum=376
           elseif(kext==2953) then
               fortfilenum=377
           elseif(kext==3217) then
               fortfilenum=378
           elseif(kext==2405) then
               fortfilenum=379
           elseif(kext==2407)then
               fortfilenum=380
           elseif(kext==2284) then
               fortfilenum=381
           elseif(kext==2404) then
               fortfilenum=382
           endif
           if(fortfilenum.ge.0)then
c               if(fortfilenum==350)
c     +             write(*,*) 'stamp to file 350... at iter=',iter,
c     +             ', tday=', tday
               write(fortfilenum,"(2(i10,','),4(f15.7,','))")
     +         iter,kext,depth, Hg0d, Hg2, MeHg
           endif
         endif



       end
c---------------------------------------------------------
c       set atmospheric loading on the surface layer
c--------------------------------------------------------
c********************************************************************

      subroutine load0d_merc(dtday,cds,loads,vol)

c integrate loadings

      implicit none

      integer nstate          !total number of state parameters
      parameter( nstate =     3 )

      real cds(nstate)      !source term [g]
      real loads(nstate)      !loading for c [g/ day)on the element]
      real vol            !volume of box [m**3]
      real dtday
      integer i

      loads(1)=0.
      loads(2)=0.
      loads(3)=0.

c      write(*,*) loads(nstate),'loads'
c        write(6,*) cds(1),cds(2),cds(3),dt,'cds before load'
      do i=1,nstate
        cds(i) = cds(i) +  loads(i)*dtday
      end do

c       write(6,*) loads(1),loads(2),loads(3),dtday,'loads'
c        write(6,*) cds(1),cds(2),cds(3),dt,'cds after load'
      end

c --------------------------------------------------------
c	end do !END OF TIME CYCLE
c---------------------------------------------------------
c---------------------------------------------------
c********************************************************************

      subroutine  merc_gas_exchange(salin,temp,area,uwind10,Hg0,Hg0atm,
     &               flux2)
c
c 4.05.2017	dmc	

c	computing mercury exchange between air and the water column
c	Sorensen et al., 2010 An Improved Global Mercury Model
c	for Air Sea Exchange of Mercury: High Concentrations over the North
c	Atlantic


	implicit none
        logical bvis1		!select viscosity calculation bvis1 true
	save bvis1		!bvis true use model output else use Soerensen

c	general constant s

c	 parameters 
        real AW		!Constant based on the Weibull distribution of 
c			!wind speeds over oceans
        real mw		!molecular weight of water [g mol-1]\
        real molHg	!molal volume of mercury at its normal boiling
                        !temperature [cm3 mol-1]
        real phiw	!solvent association factor introduced to define 
			!the effective molecular weight of the solvent
			! with respect to the diffusion process 
	
c	auxiliar variables	
         real a,b, p,c,d
         real e,g,h,m,n,o
         real v1,v2,v3,v4,v5,v6


c	from the hydrodynamic model
     	real temp	!water temperature °C
     	real salin	!water salinity
       	real uwind10	!wind speed normalised at 10 m above sea surface
     	real area	!element surface

c  	from mercury module
      	real Hg0	!Hg0 concentration in water [g/ m3] or mg/L
       	real Hg0atm	!Hg0 concentration in air   [g/ m3] or mg/L

c	variables and parameters calculated in this routine
	
     	real mex	!Air_sea Exchange of mercury Hg0 at each time step
c			![kg*s-1]
     	real mex2	!alternative formulation for test
     	real flux2	!alternative formulation for test
     	real flux	!Hg0 air-sea exchange flux [ng m-2 h-1]
     	real Hlw	!dimensionless Henry's law constant
     	real kw		!water side mass transfer coefficient for steady winds
      	real kwb	!water side mass transfer coefficient Borgest et al.2004
			!kwb Borges et al., 2004 formulation for microtidal syst
      	real SchHg	!Schmidt number for mercury
      	real ScCO2	!Schmidt number for CO2
     	real kvis	!kinematic viscosity [cm2 s-1]
     	real diff	!diffusivity (Wilke-Cang method) [cm s-1]
     	real visc	 ! water viscosity [cP]
        real rhow 	! density of the (sea)water  (KG/M**3)
      	real tempk	!temperature [K]
    	
        mw=18.0		!molecular weight of water [g mol-1]
        molHg= 12.74	!molal volume of mercury at its normal boiling
     			!temperature [cm3 mol-1]
        phiw=2.26	!solvent association factor 
       	AW=0.25		!Weibull Constant based on wind distribution

c       from the hydrodynamic model

	tempk=temp+273.15	!temperature Kelvin

C ======================================================================
C ======================================================================
C
C Compute the density and the dynamic viscosity of water from the temperature
C and the salinity

C compute the dynamic/molecular viscosity
c      VISC0=1.802863d-3 - 6.1086d-5*TEMP + 1.31419d-06*TEMP**2 -
c       &1.35576d-08*TEMP**3 + 2.15123d-06*SALIN + 3.59406d-11*SALIN**2

        a=0.0001529
        b=0.000016826
        p=1.013253
        c=8.3885*(10**(-8))
        d=p**(2)
        e=0.0024727
        g=4.8429*(10**(-5))
        h=4.7172*(10**(-6))
        m=7.5986*(10**(-8))
        n=6.0574*(10**(-6))
        o= 2.676*(10**(-9))
        v1= temp*(0.06144-temp*(0.001451-temp*b))
        v2=a*p
        v3=c*d
        v4=e*salin
        v5=(n*p-o*d)*temp
        v6=((temp*g)-temp*(h-temp*m))*salin
 
        visc=(1.791- v1-v2+v3+v4+ v5+v6)/1000     ![kg m-1* s-1]

c mpute the water density according to Brydon et al. 1999, J. Geoph. Res.
C 104/C1, 1537-1540, equation 2 with Coefficient of Table 4, without pressure
C component. Ranges TEMP -2 - 40øC, S 0-42, surface water.
C      RHOW=9.20601d-2 + 5.10768d-2*TEMP + 8.05999d-1*SALIN
C     &     -7.40849d-3*TEMP**2 - 3.01036d-3*SALIN*TEMP +
C     %     3.32267d-5*TEMP**3 + 3.21931d-5*SALIN*TEMP**2
C      RHOW=RHOW+1000d0

C compute the water density according to EOS80, Fofonoff 198599,
C J. Geoph. Res. 90/C2, 3332-3342, without pressure component.
c	[kg * m-2]

      RHOW=999.842594d0 +6.793952d-2*TEMP -9.095290d-3*TEMP**2
     &   +1.00168d-4*TEMP**3 -1.120083d-6*TEMP**4 +6.536332d-9*TEMP**5
     & +(8.24493d-1 -4.0899d-3*TEMP +7.6438d-5*TEMP**2
     &   -8.2467d-7*TEMP**3 +5.3875d-9*TEMP**4) * SALIN
     & +(-5.72466d-3 +1.0227d-4*TEMP -1.6546d-6*TEMP**2) * SALIN**1.5d0
     & +4.8314d-4*SALIN**2

	bvis1=.true.
	if(bvis1)then

	kvis=(VISC/RHOW)*10000	![cm2 s-1]
	else

	kvis=0.017*exp(-0.025*tempk)
	
	end if


C ======================================================================
C ======================================================================

        diff=((7.4*0.00000001)*((phiw*mw)**0.5)*tempk) 
        diff=diff/(visc*1000*(molHg**0.6))      !diffusivity

	Hlw=exp((-2403.3/tempk)+6.92)

	ScCO2=0.11*temp**2-6.16*temp+644.7	!Soerensen da Poissant et al 2000
c       ScCO2=2073.1+125.62*temp+3.6276*temp*temp-0.043219*temp**3 !Wanninkhof salt
c       ScCO2=1911.1+118.11*temp+3.4527*temp*temp-0.04132*temp**3 !Wanninkhof fresh
        SchHg=kvis/diff
c        write(*,*) Aw,uwind10, SchHg, ScCO2
        kw=Aw/100*24*uwind10**2*(SchHg/ScCO2)**(-0.5)  ![m day-1]
        flux=(kw)*(Hg0-Hg0atm/Hlw)       !g m-2 day-1] FIXMME
        !mex=area*flux/1000    ![kg s-1]

c	alternative equation of 

        kwb=(0.1+2.26*uwind10*(SchHg/ScCO2)**(-0.5))/100*24  ! [m day-1] Borges et al., 2004
        flux2=(kwb)*(Hg0-Hg0atm/Hlw)        ![g m-2 day-1]
        !mex2=area*flux2/1000            ![kg day-1]


c        write(97,*)visc,RHOW,diff,temp
c        write(98,*)ScCO2,SchHg, Hlw
c	write(82,*)kw,kwb
c	write(71,*)uwind10,Hg0,Hg0atm,flux,flux2 


	end			!end of routine






