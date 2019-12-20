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
! sediment dynamics in water for mercury subroutines (sed4merc_water.f)
! revision log :
!                     dmc= donata canu; gr= ginevra rosati; cl= célia laurent
!
! 30.05.2018    dmc&gr integration of the 0d module into SHYFEM
! 10.07.2019    cl     3d sinking (sum from all water levels), debug

!*****************************************************************
        subroutine sed4merc_water(bbottom,dtday,tday,wat_vol,
     +                          wdepth,k,temp,sal,taub,area,
     +                          C,Dssink,Dpsink,Vds,Vdp,
     +                          ds_gm2s, dp_gm2s)

        implicit none

        integer nstate
	parameter (nstate=2)

        real C(nstate)	   !sed variables: C(1):Silt, C(2):POM  [mg/L]
        real CD(nstate)	   !derivatives CD(1) Dsilt, CD(2) DPOM [g/d]
        integer m              !time indicator for write outputs, debug
        integer k              !node
        logical bbottom        !is on bottom?
        real cold(nstate)      !old state variables
        real dtday,tday        !time step [day], time [day]
        real temp,sal          !temperature [C] and salinity [-]
        real tkel              !temperature [K]
        real tkref             !reference temperature, 293 [K]
        real loads(nstate)     !atmospheric loadings

       real Sw, POMw           !State variables, silt and POM [mg/L]
       real dsilt,sdens, spd   !diameter [m], particle density [g/cm3], particle density [kg/m3] for silt
       real dPOM, podens,ppd   !diameter [m], particle density [g/cm3], particle density [kg/m3] for POM
       real vis, swd           !viscosity [Pa s-1] and density [kg/m3] for seawater
       real g                  !acceleration gravity [m/sec2]

       real wdepth, area, wat_vol   !depth [m], area [m2] and volume [m3] of water elements
       real Vss,Vsp                 !Stoke's settling vel for silt and POM               [m/s]
       real taub                    !bottom stress from subssed.f                         [Pa]
       real tCDs                    !critical shear stress for deposition                 [Pa]
       real Swm,POMwm               !masses of silt and POM in water                       [g]
       real ter1                    !intermediate term for calculation                   [m/s]
       real Pd                      !Pd, probability of deposition                         [-]
       real Vds, Vdp                !Pd x Stoke's velocity for silt and POM              [m/s]
       real ds_gm2s, dp_gm2s        !deposition Flux for silt and POM                 [g/m2 s]
       real dep_gm2s                !total deposition flux (silt+POM)                 [g/m2 s]
       real Dssink,Dpsink,Dsink     !Sink of silt (Dssink), POM (Dpsink) and sum (Dsink) [g/s]

       integer ipext,ipint,kext     !nodes external and internal numbers
       integer fortfilenum, iter
       logical constant_parameters
c       ____________________________________________________
c _____ Assigne old value to variables________________________

       Sw=C(1)       ![mg/L]
       POMw=C(2)
c      _______________________________________________________
c
        tkel=temp+273

c _____ Enable the use of constant parameters for debug ________________________

        constant_parameters=.False.

        if(constant_parameters)then
          vis = 0.0012
          swd = 1030.0
          area=2000.
          wdepth = 1.8
          temp = 15.
          wat_vol = wdepth*area

        else
          call sed4merc_gas_exchange(sal,temp,area,vis,swd)
        endif
c ______________________________________________________________________________
c _____ silt and POM particle properties ________________________
                                            ! move to INITI diameters FIXME
        dsilt   = 2./10.**5.        ![m]
        dPOM    = 5./10.**5.
        sdens   = 2.65              ![g/cm3]
        podens  = 1.25
        spd = sdens*10.**3.         ![kg/m3]
        ppd = podens*10.**3.
c      _________________________________________________________________________

c _____ Input critical shear for deposition ____________________________________
                                            ! move to INITI (?) FIXME
       tCDs = 0.5                   ![Pa]

c _____ Compute Stoke''s settling velocities for silt and POM ____________________

       g  = 9.81                       ![m/s2]
       ter1 = g/(18.*vis)
       Vss= ter1*(spd-swd)*(dsilt**2.) ![m/s]
       Vsp= ter1*(ppd-swd)*(dPOM**2.)

c _____ Deposition Occurrence ____________________
        if (taub < tCDs) then
          Pd = (1. - taub/tCDs)
        else
          Pd = 0.
        end if

       Vds = Pd*Vss                    ![m/s]
       Vdp = Pd*Vsp
       ds_gm2s = Vds*Sw
       dp_gm2s = Vdp*POMw              ![g/m2s]
c
       Dssink = ds_gm2s *area          ![g/s]
       Dpsink = dp_gm2s *area
       Dsink  = Dssink + Dpsink
c
       dep_gm2s = Dsink/area

c _____  Mass of Silt and POM in water ____________________

       Swm   = Sw*wat_vol            ![g]
       POMwm = POMw *wat_vol

        if (POMw .LE. 0.0) then
        write(6,*)'instability - negative POMw in sed4merc_wat kext=',
     +  ipext(k)
        stop
        else if (Sw .LE. 0.0) then
        write(6,*)'instability - negative siltw in sed4merc_wat kext=',
     +  ipext(k)
        stop
        end if

C       _________________________________________________________

	C(1)=Sw
	C(2)=POMw
c       __________________________________________________________

 	CD(1) = -Dssink *86400
 	CD(2) = -Dpsink *86400
	call merc_euler (2,dtday,wat_vol,wat_vol,c,cold,cd)    ! c(i)=( c(i)*vol+dt*cd(i) )/vol

c _____  Extract SPM concentrations at nodes close to sampling stations _______

      iter=nint(tday*86400.)

      if (MOD (iter,1800) .EQ. 0) then
           kext=ipext(k)
           fortfilenum=-1
           if(kext==3985)then
               fortfilenum=250
           elseif(kext==3986) then
               fortfilenum=251
           elseif(kext==3982) then
               fortfilenum=252
           elseif(kext==4007) then
               fortfilenum=253
           elseif(kext==3763) then
               fortfilenum=254
           elseif(kext==3765) then
               fortfilenum=255
           elseif(kext==3764)then
               fortfilenum=256
           elseif(kext==3762) then
               fortfilenum=257
           else if(kext==2150)then
               fortfilenum=258
           elseif(kext==2009) then
               fortfilenum=259
           elseif(kext==2359) then
               fortfilenum=260
           elseif(kext==2358) then
               fortfilenum=261
           elseif(kext==2341) then
               fortfilenum=262
           elseif(kext==2408) then
               fortfilenum=263
           elseif(kext==2191)then
               fortfilenum=264
           elseif(kext==2192) then
               fortfilenum=265
           elseif(kext==2654)then
               fortfilenum=266
           elseif(kext==2505) then
               fortfilenum=267
           elseif(kext==2655) then
               fortfilenum=268
           elseif(kext==2653) then
               fortfilenum=269
           else if(kext==1372)then
               fortfilenum=270
           elseif(kext==1375) then
               fortfilenum=271
           elseif(kext==1331) then
               fortfilenum=272
           elseif(kext==1378) then
               fortfilenum=273
           elseif(kext==3216) then
               fortfilenum=275
           elseif(kext==3057)then
               fortfilenum=276
           elseif(kext==2953) then
               fortfilenum=277
           elseif(kext==3217) then
               fortfilenum=278
           elseif(kext==2405) then
               fortfilenum=279
           elseif(kext==2407)then
               fortfilenum=280
           elseif(kext==2284) then
               fortfilenum=281
           elseif(kext==2404) then
               fortfilenum=282
           endif
           if(fortfilenum.ge.0)then
               if(fortfilenum==250)
     +             write(*,*) 'stamp to file 250... at iter=',iter,
     +             ', tday=', tday
               write(fortfilenum,"(2(i10,','),4(f15.7,','))")
     +         iter,kext,wdepth, Sw, POMw, taub
           endif
         endif
        end
c---------------------------------------------------------
c       set atmospheric loading on the surface layer
c--------------------------------------------------------
c********************************************************************

      subroutine load0ds(dt,cds,loads,vol)
c integrate loadings

      implicit none

      integer nstate          !total number of state parameters
      parameter(nstate = 2)

      real cds(nstate)     !source term [g]
      real loads(nstate)   !loading for c [g/ day)on the element]
      real vol             !volume of box [m**3]
      real dt

      integer i

      do i=1,nstate
        cds(i) = cds(i) +  loads(i)*dt
      end do
      end

c********************************************************************
	subroutine  sed4merc_gas_exchange(salin,temp,area,visc,rhow)

C Compute the density and the dynamic viscosity of water from the temperature
C and the salinity

	implicit none
        logical bvis1	!select viscosity calculation bvis1 true
	 save bvis1		!bvis true use model output else use Soerensen

c	auxiliar variables
         real a,b, p,c,d
         real e,g,h,m,n,o
         real v1,v2,v3,v4,v5,v6

c	from the hydrodynamic model
     	real temp	!water temperature [°C]
     	real salin	!water salinity
     	real area	!element surface

     	real visc	 ! water viscosity [cP]
       real rhow 	! density of the (sea)water  [kg/m**3]
      	real tempk	!temperature [K]

	tempk=temp+273.1

C ======================================================================
C ======================================================================
C compute the dynamic/molecular viscosity
c      VISC0=1.802863d-3 - 6.1086d-5*TEM + 1.31419d-06*TEM**2 -
c       &1.35576d-08*TEM**3 + 2.15123d-06*SALIN + 3.59406d-11*SALIN**2

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

c compute the water density according to Brydon et al. 1999, J. Geoph. Res.
C 104/C1, 1537-1540, equation 2 with Coefficient of Table 4, without pressure
C component. Ranges TEM -2 - 40øC, S 0-42, surface water.
C      RHOW=9.20601d-2 + 5.10768d-2*TEM + 8.05999d-1*SALIN
C     &     -7.40849d-3*TEM**2 - 3.01036d-3*SALIN*TEM +
C     %     3.32267d-5*TEM**3 + 3.21931d-5*SALIN*TEM**2
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

	end			!end of routine
