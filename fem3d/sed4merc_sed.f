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
! sediment dynamics in the seabed for mercury subroutines (sed4merc_sed.f)
! the sediment bed is initialized from %OC in mercury.f (L245 - L257)

! revision log :
!			                	dmc= donata canu; gr= ginevra rosati; cl= célia laurent
! 30.05.2018	dmc&gr integration of the 0d module into SHYFEM
! 01.07.2019	gr&cl  sediment active layer added (Carniello et al., 2012)

!*****************************************************************

      subroutine sed4merc_sed(k,dt,area,C,wat_vol,
     +                 taub,Cs,Dssink_sum,Dpsink_sum,
     +                 Sres, Pres,Vr,Bvels,Bvelp,
     +                 ds_gm2s, dp_gm2s,
     +                 dZbedk,dZactivk)


       implicit none

      integer nstate
	parameter (nstate=2)
      real C(nstate)	     !variables in water: C(1):Silt, C(2):POM            [mg/L]
      real CD(nstate)	     !derivatives in water CD(1) Dsilt, CD(2) DPOM       [g/d]
      real CS(nstate)        !sed  variables in sediment: C(1):Silt, C(2):POM  [g/m3]
      real CDS(nstate)	     !deriv in sediment CD(1) Dsilt, CD(2) DPOM
      real cold(nstate)      !old state variables
      real dt                !time step                                        [day]

      real area              !elements area                                    [m2]
      real wat_vol           !elements water volume                            [m3]
      real sed_vol_old       !volume of sediment elements at time t-1          [m3]
      real sed_vol_new       !volume of sediment elements at time t            [m3]

      real sw, POMw              !silt and POM in water                        [mg/L]/[g/m3]
      real silt, POM             !silt and POM in sediment                     [g/m3]
      real siltm, POMm           !silt and POM mass in sediment                [g]
      real silt_s0,POM_s0        !silt and POM in layer0 (below active layer)  [g/m3]

      real p_POC, p_POM, p_silt  !OC, POM and silt in sediment                 [%]
      real OC_mg_g, OM_mg_g      !OC and POM in sediment                       [mg/g]
      real Pdens                 !sediment particle weighted density           [g/m3]
      real DryD, BulkD           !Dry and Bulk density                         [g/cm3]
      real por                   !sediment porosity                            [-]

      real Dssink_sum, Dpsink_sum      !deposition for silt and POM            [g/sec]
      real ds_gm2s, dp_gm2s, dep_gm2s  !deposition flux for silt, POM and tot  [g/m2s]

      real taub,tCE                 !Bottom Shear Stress and Critical Shear Stress for Erosion [Pa]
      real MagR, Er                 !Magnitude of resuspension and Erosion flux       [-] and  [g/m2s]
      real dme, logdme, es, dme2     !terms for computing resuspension multiplier dme2
      real Vr                       !resuspension velocity                                    [m s-1]
      real Pres, Sres, Tres         !resuspension of silt, POM, and total                     [g/s]
      real rs_gm2s,rp_gm2s,res_gm2s !resuspension flux for silt, POM, and total               [g/m2s]

      real Bvels, Bvelp      !burial velocity of silt and POM                        [m/s]
      real Bvel      !burial velocity of total solids                                [m/s]
      real Bsflux, Bpflux    !burial of silt and POM                                 [g/s]
      real dZactivk          !thickness of the active layer                          [m]
      real dZact0            !critical thickness of the active layer                 [m]
      real dZdig             !thickness of sediment layer eroded                     [m]
      real dZit              !active layer variation during current iteration        [m]
      real dZbedk            !total variation of the sea bed depth since the beginning of the simulation
      real prct_0,prct_c     !fraction of solids from layer0 and active layer        [-]

      integer m              !time indicator for write outputs, debug
      integer k              !node
      integer ipext          !node external number
      logical constant_parameters

        constant_parameters=.False.
        if(constant_parameters)then
          area=2000.
          wat_vol= 1.8*area
        endif

        dZact0=0.02
        silt_s0 = 513767.
        POM_s0 =  27613.
c       _______________________________________________________
c       assigne old value to water (sw,POMw) and sediment (silt, POM) variables

        sed_vol_old=area*dZactivk

        sw=C(1)
        POMw=C(2)

        silt=CS(1)
        POM=CS(2)

c       __________________________________________________
C       processes sed4merc_sed
c       __________________________________________________

              p_POM = POM/(POM+silt)*100.
              DryD = (POM+silt)/10.**6.
              OM_mg_g = 10.0 * p_POM
              OC_mg_g = OM_mg_g/1.7
              p_silt=100-p_POM

c___________ Compute weigthed particle density [g cm-3], porosity [-], Bulk density [g(s+w) cm-3]

       Pdens = ((1.25*p_POM)+(2.65*(100.0 - p_POM)))/100.0  ![g(s)/cm3(s)]

       por = 1.0 - (DryD/Pdens)        ![L(w)/L(s+w)]

       BulkD = por + DryD              ![g/cm3]
c
c____________ Input critical shear for erosion tCE ____________
c
       tCE  = 0.44                     ![Pa]

c___________ Compute resuspension multiplier(dme/dt) from BulkD _________
c___________ from Hwang and Metha (1989) in Tetra Tech (2002) _________
c
        es = 0.198/(BulkD-1.0023)
        logdme = 0.23*exp(es)
        dme = 10.**logdme                      ![mg cm-2 hr-1]
        dme2 = ((dme/1000.)/3600.)*10.**4.     ![g m-2 s-1]
c
c___________ Resuspension Occurrence _________________________________

        if (taub > tCE) then
          Vr = dme2/(POM+silt)
          MagR = (taub/tCE - 1.) !**2
        else
          MagR = 0.
          Vr = dme2/(POM+silt)
        end if

c___________ Compute Resuspension Fluxes and rates _________________________

       Er = dme2 * MagR                  ![g m-2 s-1]
       Sres  = (Er*area)*(p_silt/100.)  ![g s-1] = [g m-2 s-1] * [m2]
       Pres  = (Er*area)*(p_POM/100.)
       Tres   =  Sres + Pres
c
       rs_gm2s = Sres/area
       rp_gm2s = Pres/area
       res_gm2s= Tres/area

c___________ Compute Burial Velocities and Fluxes ____________________________

       Bvels = (ds_gm2s-rs_gm2s)/silt    ![m s-1]
       Bvelp = (dp_gm2s-rp_gm2s)/POM
c
       Bvel = (dep_gm2s-res_gm2s)/(silt+POM) ![[m s-1]
c
       Bsflux = Bvels*silt*area
       Bpflux = Bvelp*POM*area

       siltm = silt*sed_vol_old         ![g]
       POMm = POM*sed_vol_old
c
        CS(1)=silt            ![g/m3]
        CS(2)=POM             ![g/m3]

        C(1) = sw
        C(2) = POMw
c        Bsflux=0.0
c        Bpflux=0.0

c___________ Compute sediment thickness variations ____________________________

        dZit=( (Dpsink_sum-Pres)/POM + (Dssink_sum-Sres)/silt)
     +         *dt*86400./area ! [m] = [g s-1]*[g-1 m3] *[day]*[s day-1]*[m-2]
        if(dZactivk+dZit<dZact0)then
          dZdig=dZact0-(dZactivk+dZit)
          prct_0=(dZdig/dZact0)
          prct_c=(dZact0-dZdig)/dZact0
        else
          dZdig=0
          prct_0=0.
          prct_c=1.0
        endif
        cs(1)=cs(1)*prct_c+silt_s0*prct_0
        cs(2)=cs(2)*prct_c+POM_s0*prct_0
        sed_vol_new=area*(dZactivk+dZit+dZdig)
        dZactivk=dZactivk+dZit+dZdig
c_______________________________________________________________________________
c
c____________Positive burial push sediment below the active layer_______________
c____________Negative burial is equal to net resuspension ______________________

        if (Bsflux .GE. 0.) then
              CDS(1) =(Dssink_sum-Sres-Bsflux) *86400  ![g/day] silt sed
        else
              CDS(1) =(Dssink_sum-Sres) *86400
        end if

        if (Bpflux .GE. 0.) then
              CDS(2) =(Dpsink_sum-Pres-Bpflux) *86400 ![g/day] pom sed
        else
              CDS(2) =(Dpsink_sum-Pres) *86400
        end if
c
        CD(1) = (+Sres) *86400.    !variation in days
        CD(2) = (+Pres) *86400.

        write (88,*) (CS(m), m=1,nstate),'cs variables beforer'
        write (88,*) (C(m), m=1,nstate),'c variables beforer'

        call merc_euler (2,dt,sed_vol_old,sed_vol_new,cs,cold,cds) !claurent-OGS: vol old and new are different for sea bed sediment
c       call merc_euler (2,dt,sed_vol_old,cs,cold,cds)  ! c(i)=( c(i)*vol+dt*cd(i) )/vol

         call merc_euler (2,dt,wat_vol+dZbedk*area,
     +                        wat_vol+(dZbedk+dZit)*area, c,cold,cd)  ! c(i)=( c(i)*vol+dt*cd(i) )/vol
c       call merc_euler (2,dt,wat_vol+dZbedk*area, c,cold,cd)  ! c(i)=( c(i)*vol+dt*cd(i) )/vol

        dZbedk = dZbedk+dZit  ! update sea bed depth

        c(1)=max(0.001,c(1))
        c(2)=max(0.001,c(2))  ! Celia: allow POM_w to be NULL (for drying nodes)

        ! Celia: if this trick is kept, unsinked silt and pom should be
        ! substracted from the bottom sediment (doing an inversed euler)
        ! in order to maintain mass conservation


        write (88,*) (CS(m), m=1,nstate),'cs variables after' !, ' Transformations'
        write (88,*) (C(m), m=1,nstate),'c variables after' !, ' Transformations'
        write (88,*) wat_vol,cds,cd, 'wat_vol and derivative cds,cd'

       silt   = CS(1)
       POM    = CS(2)

        sw = C(1)
        POMw=C(2)

       if (POMw .LE. 0.0) then
       write(6,*) 'instability - negative POM_w in sed4merc_sed kext=',
     +   ipext(k)
       stop
       else if (sw .LE. 0.0) then
       write(6,*) 'instability - negative silt_w in sed4merc_sed kext=',
     +   ipext(k)
       stop
       end if
       end
c********************************************************************
