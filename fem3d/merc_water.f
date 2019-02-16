
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

c	program merc.f

        subroutine mercury_react(id,bsurf,bbottom,boxtype,dt,vol,depth,k
     +                          ,temp,qrad,C,loads)
        
c

        implicit none

        integer nstate
	parameter (nstate=3)

        real C(nstate)		!mercury variables C(1):Hg1, C(2):Hg2, C(3):MeHg
        real CD(nstate)		!derivatives CD(1):diffHg0, CD(2):diffHg2, CD(3):diffMeHg
        integer m       !time indicator for write outputs, debug
        integer id              !id of node and level
        logical bsurf           !is on surface?
        logical bbottom         !is on bottom?
        logical boxtype
        real cold(nstate)       !old state variables
        real dt         !time step, [day]
        real qrad
        real  temp,light
        real depth              ! box property
	integer k
        real vol                !box property
        real tkel       !temperature K
        real tkref      !reference temperature, 293 K
c
       real loads(nstate)      !atmospheric loadings

c	variables
c
        real Hg0d, Hg0a
        real Hg2,Hg2d, Hg2DOC,Hg2sorb
        real MeHg,MeHgd,MeHgDOC,MeHgsorb
        real MeHgsed, Hg2dsed
        real Hg2DOCs
        real sksme
c
c	partition coefficients
c
      	real K1silt !           !Hg2 in silt
        real K1sand ! !Hg2 in sand
        real K1doc  ! !Hg2 in doc
        real K1org  ! !Hg2 in organic particles
        real K2silt ! !MeHg in silt
        real K2sand ! !MeHg in sand
        real K2doc  ! !MeHg  in doc
        real K2org  ! !MeHg in organic particles
        real partden1   !auxiliar variable in computation
        real partden2   !auxiliar varible in computation
        real Hg2silt
        real Hg2sand
        real Hg2org
        real MeHgsilt
        real MeHgsand
        real MeHgorg
        real faq1
        real fsilt1
        real fsand1
        real forg1
        real fdoc1
        real ftot1
        real ftot2
        real faq2
        real fsilt2
        real fsand2
        real forg2
        real fdoc2

c	solids in water
c
        real DOC
        real silt  
        real sand  
        real org  
        
c       

c	Hg0 sink reactions, parameters
c
        real kvol, ktvol                !volatilization
        real ckvol, skvo
        real He, R, OxAct, kox          !oxydation
        real Rcal, ckox, skox
c
c	HgII sink reactions, parameters
c
        real ksmeth, kmeth, Qbac,ckmeth,skme    !bacterial methylation
        real cksmeth
        real kphr, lref, ke             !photoreduction
        real ladj,skph,ckph             !photoreduction
        real xHg2d,xHg2DOC,xHg2sorb
c
c	MeHg sink reactions, parameters
c	bacterial demethylation
        real skdem,kdem
        real ckdem, Eadem, cordem
c	
c	sediment demethylation
c
        real cksdem,ksdem
        real sksdem

c	photoreductive demethylation
c
        real skphdem
        real ckphdem
        real kphdem
        real xMeHgd,xMeHgDOC,xMeHgsorb


c	---------------------------------------------------

        DOC=3   !mg/L
        silt=2.5        !mg/L
        sand=0.5!mg/L
        org=2.5 !mg/L (INPUT)

        Hg0a=0.0016     ![ng/L] atmospheric Hg concentration (INPUT)

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
        k1silt=200000   !Hg2 in silt
        k1sand=100      !Hg2 in sand
        k1doc=200000    !Hg2 in doc
        k1org=400000    !Hg2 in organic particles
        k2silt=10000    !MeHg in silt  230000
        k2sand=0        !MeHg in sand
        k2doc=100000    !MeHg  in doc
        k2org=500000    !MeHg in organic particles
c       -----------------------------------------------------

        
c
c	------------------------------------------------------
c	initial conditions: mercury species concentration in [ng/L] or  [mg/kg]  (sed)
c	------------------------------------------------------
c	-----------------------------------------------------
c        Hg0d=0.05 !1.42	!initial concentration of dissolved Hg0d
c        Hg2=4   !FIXME from serafm simulation
c        MeHg=0.25 !0.01	!FIXME 	from SERAFM simulation
c FIXME forse togliere	Hg2dsed= 1.23	!FIXME  from SERAFM simulation
c FIXME forse togliere	MeHgsed= 0.34	!FIXME  from SERAFM simulation
c	----------------------------------------------------
c	global constants
c	----------------
        tkref=293       !reference temperature, K
c  	------------------------------------------------------
c       volatilization
c	--------------
        kvol=0.1        ![day-1] volatilization rate default value at 20°C
        ktvol=1.04      !volatilization temperature correction Theta for elemental mercury
        He=0.0071       ![atm/mol] Henry 's Law constant
        R=8.314472      ![J K-1 mol-1] universal gas constant
c        write(6,*) kvol
c
c       oxydation of Hg0d
c	-----------------
        OxAct=2 !activation energy for oxydation, WASP impl. FIXME =2 k[kcal mol-1]
        kox=0.01    !0001  !FIXME oxydation parameter 
c
c	methylation
c	----------
        ksmeth=0.006 ![day-1] sed methylation rate at 20°C Lignano:0.006 Hines et al., 2012
c       ksmeth=0.013 ![day-1] sed methylation rate at 20°C S.Andrea:0.013 Hines et al., 2012
c       ksmeth=0.006 ![day-1] sed methylation rate at 20°C Buso:0.006 Hines et al., 2012
c       ksmeth=0.006 ![day-1] sed methylation rate at 20°C Morgo:0.006 Hines et al., 2012
c       ksmeth=0.009 ![day-1] sed methylation rate at 20°C Grado:0.009 Hines et al., 2012
c       ksmeth=0.009 ![day-1] sed methylation rate at 20°C Primero:0.009 Hines et al., 2012

        kmeth=0.006 ![day-1] wq  methylation rate at 20°C Monperrus et al., 2007
        Qbac = 1.5      !FIXME 
c
c       photoreduction of Hg2d to Hg0
c	----------------------------
        kphr=0.05       ![day-1] photoreduction rate constant of HgII to Hg0 (default WASP) 
        lref=240 !950                !reference light intensity for kph [watt/m2] FIXME (default WASP)
        ke=1.05             !light extintion coefficient FIXME
c	
c	ke=0.0525	!lignano
c	ke=2.25		!s. andrea
c	ke=0.0525	!buso, morgo
c	ke=1.05		!grado, primero

c	------------------------------------
c	photoreductive demethylation
c	--------------------------
        kphdem=0.015    ![day-1] source of this parameter? FIXME
c	------------------------------------------------------------------------
c	bacterial demethylation
c	---------------------
        ksdem=0.159     !bacterial sediment demethylation
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

c        write(6,*) 'partition fractions: ',faq1,faq2,fsilt1,fsilt2
c        write(70,*) faq1, fsilt1, fsand1, forg1, fdoc1, ftot1!test
c        write(69,*) faq2, fsilt2, fsand2, forg2, fdoc2, ftot2
c	-----------------------------------
c	calculation of mercury fractions amongst the different phases
c	------------------------------------
c        write(6,*) Hg2, 'Hg2?'
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

        if(boxtype) then

c
c		-----------------------------
c		Hg0d--> Hg0a volatilization
c
        if(bsurf) then

        ckvol=kvol*(ktvol**(temp-20))   !calculated volatilization rate at temperature temp
c        write(6,*) ckvol, 'ckvol'
        skvo=ckvol*(Hg0d-(Hg0a/(He/R*tkel))) ![g/m3/day] of Hg0 exchanged with the atmosphere.
        else
        skvo=0
        end if

c        write(6,*) skvo,'skvo',bsurf, 'bsurf'
c			--------------------
c	Hg0d --> Hg2d oxydation
c
c	conversion of R, gas constant, from [J K-1 mol-1] to [cal k-1 mol-1]
        ckox=kox*exp(OxAct*1000*((tkel-tkref)/(Rcal*tkel*tkref)))	!1000: conversion from kcal to cal
        skox=ckox*Hg0d
c        write(81,*) ckvol,tkel,skvo,skox !'volatilization,oxydation'
c       ok

c	----------------------------------------
c		Hg2d--> Hg0d photoreduction 
c
        light=qrad
        ladj=light/lref*((1-exp(-ke*depth))/ke*depth)	
        ckph=kphr*ladj   !FIXME unità di misura
	skph=ckph*(Hg2d*xHg2d+Hg2DOC*xHg2DOC+Hg2sorb*xHg2sorb)	!FIXME unità di misura
c	write (82,*) skph, ladj,light, ckph,  ' photoreduction'


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
	write (85,*) skphdem,ckphdem !, ' photochemical demethylation'
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
	write (86,*) tkel,cordem,ckdem,skdem	!, ' water column bacterial demethylation'
	
c       __________________________________________________
	

	C(1)=Hg0d
	C(2)=Hg2
	C(3)=MeHg
c	write(6,*) C(1),Hg0d,C(2),Hg2,C(3),MeHg,'merc var'
c	
c	CD= transformations 1:Hg0 2:Hg2 3:MeHg	
 	CD(1) = - skvo - skox + skph + skphdem
 	CD(2) = skox - skph + skdem - skme
 	CD(3) = skme - skdem - skphdem
        write (88,*) skvo,(CD(m), m=1,nstate),'f.88' !, ' Transformations'
        write (88,*) (C(m), m=1,nstate),boxtype,'variables' !, ' Transformations'

c	integration

c        write(6,*) vol,c,cold,cd        !
c
        call load0d_merc(dt,cd,loads,vol)
	call merc_euler (3,dt,vol,c,cold,cd)
        write (88,*) (C(m), m=1,nstate),boxtype,'variables after' !, ' Transformations'
	
        else
        
c       ______________________________________________
c       Sediment
c       Sediment methylation: Hg2d --> MeHgd
c
        ckmeth=ksmeth*Qbac**((temp-20)/10)
        skme=ckmeth*(Hg2d*xHg2d+Hg2DOC*xHg2DOC+Hg2sorb*xHg2sorb)
c        write(84,*) kmeth, Qbac, ckmeth, skme,'mt' !,' ksmeth Qbac ckmeth
c        write(84,*) Hg2d,xHg2d,Hg2DOC,xHg2DOC,Hg2sorb,xHg2sorb,'test'
c
c       MeHg ---> Hg II sediment demethylation

        cordem=exp(Eadem*1000*((tkel-tkref)/(Rcal*tkel*tkref)))
        cksdem=ksdem*cordem
        sksdem=cksdem*(MeHgd*xMeHgd+MeHgDOC*xMeHgDOC+MeHgsorb*xMeHgsorb)
	write (84,*) tkel,ksdem,cksdem,sksdem,'dm'!, ' water column bacterial demethylation'


        C(1)=0
        C(2)=Hg2
        C(3)=MeHg

c         write(6,*) C(1),Hg0d,C(2),Hg2,C(3),MeHg,'merc var_sed'
c       CD= transformations 1:Hg0 2:Hg2 3:MeHg  

        CD(2) = skdem - skme
        CD(3) = skme - sksdem 

c        write (88,*) skvo,(CD(m), m=1,nstate),'f.88 sed' !, 'sed tr'
c        write (88,*) (C(m), m=1,nstate),boxtype,'var sed before' !, 'sed tr'

c       sarebbe da sistemare, mettendo solo due var anche in mercury)
          call merc_euler (2,dt,vol,c(2),cold(2),cd(2))
c        call merc_euler (3,dt,vol,c,cold,cd)
c        write (88,*) (C(m), m=1,nstate),boxtype,'var sed after' !, 'sed tr'
        
        end if

        end
c---------------------------------------------------------
c       set atmospheric loading on the surface layer
c--------------------------------------------------------
c********************************************************************

      subroutine load0d_merc(dt,cds,loads,vol)

c integrate loadings

      implicit none

        integer nstate          !total number of state parameters
        parameter( nstate =     3 )

      real cds(nstate)      !source term [g/day]
      real loads(nstate)      !loading for c [g/(m**3 day)]
      real vol            !volume of box [m**3]
        real dt

      integer i

c        write(6,*) cds(1),cds(2),cds(3),dt,'cds before load'
      do i=1,nstate
        cds(i) = cds(i) +  loads(i)*dt
      end do

c        write(6,*) loads(1),loads(2),loads(3),dt,'loads'
c        write(6,*) cds(1),cds(2),cds(3),dt,'cds after load'

      end



c --------------------------------------------------------
c	end do !END OF TIME CYCLE
c---------------------------------------------------------
c---------------------------------------------------
c********************************************************************



