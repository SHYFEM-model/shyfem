
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
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

! revision log :
!
! 03.04.2018	ggu	changed VERS_7_5_43
! 31.08.2018	ggu	changed VERS_7_5_49
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

c weutro_ulva
c notes:
c State variables used: (Wasp)
c
c nh3           71      1   -
c no3 , nox     72      2   -
c opo4          73      3   -
c phyto         74      4
c cbod          75      5
c do , oxy      76      6   -
c on , onsed    77      7
c op , opsed    78      8
c zoo           79      9
c
c********************************************************************
c********************************************************************
c********************************************************************
c
        subroutine wulva(k,t,dt,vol,depth,vel,stp,qrad,c,culva)

c EUTRO 0-Dimensional (Sediments)

        implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )
        integer nulstate         !total number of state parameters
        parameter( nulstate = 2 )

        integer k               !nodi
        real t                  !actual time [day]
        real dt                 !time step [day]
        real vol                !volume [m**3]
        real depth              !depth of box [m]
        real vel                !velocity [m/s]
        real stp                !temperature [C]
        real qrad               !radiazione
        real c(nstate)          !state variable [mg/L] == [g/m**3]
        real cold(nstate)       !old state variable [mg/L] == [g/m**3]
        real culva(nulstate)     !state variable (ulva) [mg/L] == [g/m**3]

        include 'donata.h'

        logical bulva
        integer i
        real Wdt                !Wdt = dry weight [gC/l] DryW
        real Quota              !Quota = intratissual Nitrogen [gN/gC dw]
        real Biom               !Biomassa

        real fi0,fugl,fugt,furt,fupl,mtu 
        real tsnh,tsno,gusp,gu,gquo,ru,fu,tsrp
        real opo4,nh3,nox,oxy,sedC
        real onsed,opsed,tsn,Wdtnew,Quonew,WetWnew,gUlva,dic

        real p_kmtu1,p_kmtu2,p_kmtdo,p_vru,p_Vnhu,p_Knhu,p_Vnou,p_Knou 
        real p_Qmax,p_Qmin,p_Vgud,p_Qlc,p_upc,p_ruc,p_ugt1,p_ugt2 
        real p_urt1,p_urt2,p_Pmax,p_ugl1,p_kest,p_kpar,p_ud2fw,p_kpo4 
        real p_roc,p_ron 
        real cds(nstate)        !source term (right hand side) [g/day]
        real ca(nstate)         !auxiliary state vector [g/m**3]
        real caold(nstate)      !old auxiliary state vector [g/m**3]

        integer icall
        save icall
        data icall / 0 /

        bulva = .false.
        bulva = .true.

        if( .not. bulva ) return

        if( icall .eq. 0 ) then
          call settopseg(.true.)        !marks segment as surface
          call setbotseg(.true.)        !marks segment as bottom
          icall = 1
        end if


!----------------------------------------------------------
!        Model parameters @ March-2015
!----------------------------------------------------------

         p_kmtu1 = 0.84                   !coefficiente calcolo mortalita' Ulva
         p_kmtu2 = 0.0012                 !coefficiente calcolo mortalita' Ulva
         p_kmtdo = 1.0                    !mortalita' oraria massima Ulva(.05)
         p_vru = 2.54                     !velocita' respirazione ulva
         p_Vnhu = 5.2                     !velocita' NH4 per Ulva
         p_Knhu = 0.7                     !costante di semisaturazio per NH4 (0.2-0.6_BSQ)
         p_Vnou = 0.9                     !velocita' NOx per Ulva
         p_Knou = 0.07                    !costante di semisaturazio per NOx (0.05-0.2_BSQ)
         p_Qmax = 45.0                    !quota massima
         p_Qmin = 10.0                    !quota minima
         p_Vgud = 0.0187                  !velocita 'crescita' giornaliera 0.45 oraria (0.36-0.050)
         p_Qlc = 8.0                      !livello critico di quota
         p_upc = 0.5                      !rapporto mgP / g peso secco per ulva  
         p_ruc = 0.25                     !rapporto carbonio/peso secco in ulva (0.30_BSQ)
         p_ugt1 = 0.2                     !coeff T° growth 1
         p_ugt2 = 12.5                    !coeff T° growth 2
         p_urt1 = 0.3                     !coeff T° respiration 1
         p_urt2 = 10.0                    !coeff T° respiration 2
         p_Pmax = 27.5                    !produzione max di ossigeno [mgO2/gdw h] Solidoro (EM94)
         p_ugl1 = 8.67                    !corrisponde a 5800 * 0.001496
         p_kest = 1.0                     !light estinction coefficient 0.751  
         p_kpar = 0.46                    !stefO8:fattore conversione radsol tot->PAR (Ciavatta 2008)
         p_ud2fw = 7.5                    !rapporto gdw/l kg fw m2 !ulva gdwl-1=7.5 Kgfw m-2 per h=1m
         p_kpo4 = 0.01                     !costante di semisaturazione per ulva per PO4
         p_roc = 2.66                     !oxygen/carbon ratio (2.89_BSQ)           
         p_ron = 4.5                      !oxygen/nitrogen ratio (0.30_BSQ)         

!----------------------------------------------------------
! variabili WASP, iniziali
!----------------------------------------------------------

        nh3   = c(1)
        nox   = c(2)
        opo4  = c(3)
        oxy   = c(4)
        Wdt   = culva(1)
        Quota = culva(2)

!----------------------------------------------------------
!      Functions
!----------------------------------------------------------
!...     Light limiting function

        fi0 = qrad * p_kpar * exp(-p_kest * depth) !converto radsol tot (W/m2)-PAR ULVA!!
        fugl = (1 - EXP(-fi0 / p_ugl1) )
c
!...    Temperature limiting function for growth and respiration
c
        fugt = (1/ (1 + exp(-p_ugt1 * (stp-p_ugt2)) ) )   !funz. crescita
        furt = (1/ (1 + exp(-p_urt1 * (stp-p_urt2)) ) )   !funz. respirazione
c
c...	Funzione limitante fosforo pb(kpo4) = 0.01
c
        fupl = opo4 /(p_kpo4 + opo4) 
c
!...    Mortality function : physiological + DO induced

        mtu = min(p_kmtu2 * (Wdt**p_kmtu1), p_kmtdo*Wdt)
!
!...    Nitrogen Uptake
!
        tsnh = ( (p_Vnhu * nh3) / (nh3+p_Knhu) ) 
     +           * ( (p_Qmax - Quota) / (p_Qmax - p_Qmin) )

        tsno = ( (p_Vnou * nox) / (nox+p_Knou) )  
     +           * ( (p_Qmax - Quota) / (p_Qmax - p_Qmin) )

!...    growth function
!
        gusp = p_Vgud * ( (Quota-p_Qmin) / (Quota-p_Qlc) )  
     +          * fugt * fugl * fupl
!
!..     Compute growth effect on dry weight and internal quota
!
        gu = gusp * Wdt
        gquo = gusp * Quota
!
!...    Respiration [mgO/L]

        ru = p_Vru * furt * Wdt
!
!...    Photosinthesis [mgO/L]
!
        fu = (p_Pmax / p_Vgud) * gusp * Wdt    ! da Solid. EM94
!
!...    Uptake rate of P
!
        tsrp = gusp * p_upc 
!
!...    Uptake from the water column
c
        ca(1)   =  nh3  
        ca(2)   =  nox  
        ca(3)   =  opo4 
        ca(4)   =  oxy  
        cds(1)  = (-1.0)*  tsrp * Wdt  ! D_opo4        
        cds(2)  = (-1.0)*  tsnh * Wdt  ! D_nh3   
        cds(3)  = (-1.0)*  tsno * Wdt  ! D_nox   
        cds(4)  = (+1.0)*  (fu - ru)   ! D_oxy          !!!proportional to growth

!...    Dead matter to sediment
        ca(7)   = c(7) ! onsed
        ca(8)   = c(8) ! opsed
        !ca(9)  = cs(3)
        cds(7)  = mtu * Quota ! D_onsed 
        cds(8)  = mtu * p_upc ! D_opsed 
        !cds(9)  = mtu * p_ruc ! D_sedC  
!
!
!...    Integration of State Variables
!
        tsn=(tsnh+tsno)           !Nitrogen cumulative uptake
        ca(5)=Wdt
        ca(6)=Quota
        cds(5)=(gu-mtu)              ! Wdtnew = Wdt + (gu-mtu) * dt  !tal vez es la biomasa
        cds(6)=((tsn-gquo)/(1+gusp)) !Quonew = Quota + ((tsn-gquo)/(1+gusp)) * dt   !è una sorta di luxury uptake

        call euler(2,dt,vol,ca(5),caold(5),cds(5)) 

        Wdtnew=ca(5)
        Quonew=ca(6)
!
!...    Computed on the basis of the gdw/l to Kgfw/m2 ratio
        WetWnew = Wdtnew * p_ud2fw
!
!...    Net growth of Ulva r.
        gUlva= max(0.,(Wdtnew-Wdt)) 
!
!...    Uptake/release of DIC        
        dic = gUlva * p_ruc        
!
!--------------------------------------------------------------
!...    Control for intratissual N quota concentration thresholds
!--------------------------------------------------------------
!...    MINIMUM QUOTA - p_Qmin
!
        if (Quonew. lt. p_Qmin) then

        Wdtnew = Wdt - mtu * dt
        Quonew  = Quota              ! è una sorta di luxury uptake
        WetWnew = Wdtnew * p_ud2fw
!
!...     Uptake from the water column
        ca(1)   = nh3   ! D_opo4 = 0.    
        ca(2)   = nox   ! D_nh3 = 0.     
        ca(3)   = opo4  ! D_nox = 0.     
        ca(4)   = oxy   ! D_oxy = - ru   
             ! D_dic = 0.    
!
!...    Dead matter to sediment 
!        D_sedC  = mtu * p_ruc   
!        D_onsed = mtu * Quota
!        D_opsed = mtu * p_upc

        endif
!
!...    MAXIMUM QUOTA - p_Qmax
!
        if (Quonew. gt. p_Qmax) then
                Quonew  = p_Qmax
        endif
!
!...    euler integration of the state variables
!
        call euler(4,dt,vol,ca(1),caold(1),cds(1)) 
        call euler(2,dt,vol,ca(7),caold(7),cds(7)) 
!
!...    Transfer to state vector
!
        c(1) = ca(1)  
        c(2) = ca(2)  
        c(3) = ca(3)  
        c(4) = ca(4)  
        c(7)= ca(7)
        c(8)= ca(8)

        if (Wdtnew.le.0.) then
        culva(1)  = 0.       !VPhyto(i,j,DryW) = 0.
        culva(2)  = 0.       !VPhyto(i,j,quo)  = 0.      
                             !VPhyto(i,j,Biom) = 0.
        else
        culva(1)  = Wdtnew   ! VPhyto(i,j,DryW) = Wdtnew
        culva(2)  = Quonew   ! VPhyto(i,j,quo)  = Quonew       
                             ! VPhyto(i,j,Biom) = WetWnew
        endif
        
        return
        end

c****************************************************************************

