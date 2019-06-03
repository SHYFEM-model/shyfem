
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

! revision log :
!
! 23.03.2010	ggu	changed v6.1.1
! 12.12.2014	ggu	changed VERS_7_0_9
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

c********************************************************************

        subroutine haka0d(time,dt,vol,h,temp,aI,e,eload)

c box ecological model for hakata bay

        implicit none

        integer nstate          !total number of state variables
        parameter( nstate = 5 )

        real time          !actual time [sec]
        real dt            !time step [sec]
        real vol           !volume of box [m**3]
        real h             !depth of box [m]
        real temp          !temerature [C]
        real aI            !light at surface of box [??]
        real e(nstate)     !state variables [g/m**3]
        real eload(nstate) !loading for state variables [g/sec]
        
        real es(nstate)    !source for equation [g/sec]

        call hakaeco(time,vol,h,temp,aI,e,eload,es)
        call euler(dt,vol,e,es)

        end

c********************************************************************

        subroutine hakaeco(time,vol,h,temp,aI,e,el,es)

c handles biochemical reactor
c
c state variables are in g/m**3
c source term is in g/sec
c
c new concentration has to be computed as:
c
c ( (CV)_new - (CV)_old ) / dt = source
c
c =>  C_new = ( C_old V_old + dt * source ) / V_new

        implicit none

        integer nstate
        parameter( nstate = 5 )

        real time       !actual time [sec]
        real vol        !volume of box [m**3]
        real h          !depth of box [m]
        real temp       !temerature [C]
        real aI         !light at surface of box [??]
        real e(nstate)  !state variables [g/m**3]
        real el(nstate) !loading for state variables [g/sec]
        real es(nstate) !source for equation [g/sec]

        real Kx,Kz,Ks,Iopt,L,Kzp,kmp,kmz,kg,kk
        parameter(Kx=140,Kz=0.0002,Ks=0.03,
     +             Iopt=6,
     +             L=10000,Kzp=0.00001, kmp=0.069,
     +             kmz=0.0693, kg=0.0693,kk=0.063)
c     +             Iopt=500000,
      
        real A1,A2,A1A2,A3,B1,B2,B3,B4,C1,C2,D1
        real aIr,T
        real PHP,ZOO,DET,DOP,DIP
        real PHP2,ZOO2
        real dPHP,dZOO,dDET,dDOP,dDIP

        real al,b,gr,f,ram,th,Vmax
        real sec
        real cmp,cmz,vpio,vpit,vpdo,vpdt,vdio,vdit

c----------------------------------------------------------------  

      cmp=1.67e-4
c       mortarity of phytoplankton at 0 degree ; /sec
      cmz=3.5e-4
c       mortarity of zooplankton at 0 degree ; /sec
c     kmp=0.069
c       temperature dependency of mortarity of phytoplankton 
c       ; /degreee 
c     kmz=0.0693
c       temperature dependency of mortarity of zooplankton 
c       ; /degreee 
      vpio=3.47e-7
c       decomposition speed of detritus to DIP at 0 degree
c       ; /sec      
      vpit=0.0693
c       temperature dependency of decomposition of detritus to DIP
c       ; /degree
      vpdo=2.89e-7
c       decomposition speed of detritus to DOP at 0 degree
c       ; /sec      
      vpdt=0.0693
c       temperature dependency of decomposition of detritus to DOP
c       ; /degree
      vdio=3.47e-7
c       decomposition speed of DOP to DIP at 0 degree
c       ; /sec      
      vdit=0.0693
c       temperature dependency of decomposition of DOP to DIP
c       ; /degree
      b=0.3
c       constant for fecal pellet generation ; 
      GR=0.2
c       maximum grazing rate of zooplankton ; /day
c      kg=0.0693
c       temperature dependency of grazing ; /degree       
      f=1.05
c       sedimentation of detritus ; mg/sec
c      Iopt=500000.
c       optimum light intensity ; cal*day/m**2
c      Ks=0.03
c       half saturation concentration ; 
c      I1=313000.    
c       average light intensity ;

       ram=723.
       th=8.33e-5
       al=0.4
       Vmax=1.2
       sec=86400.
       A2=0.135
       
c----------------------------------------------------------------  
  
        T = temp

        PHP = e(1)
        ZOO = e(2)
        DET = e(3)
        DOP = e(4)
        DIP = e(5)

        PHP2 = PHP*PHP
        ZOO2 = ZOO*ZOO
      
c----------------------------------------------------------------  

        aIr = (aI/Iopt)**2
      
        A1=(Vmax/sec)*exp(kk*T)*(aIr*exp(1.0-aIr))*(DIP/(Ks+DIP))
        A1A2=A2*A1
        A3=cmp*exp(kmp*T)
        
        B1=(GR/sec)*exp(kg*T)*(1-exp(ram*(-PHP+th)))
        B2=al*B1
        B3=b*B1
        !B2=cmz*exp(kmz*T)      !??
        B4=cmz*exp(kmz*T)
        
        C1=vpio*exp(vpit*T)
        C2=vpdo*exp(vpdt*T)
        
        D1=vdio*exp(vdit*T)
        
c----------------------------------------------------------------  

        dPHP = vol * ( A1*PHP - B1*ZOO  - A1A2*PHP - A3*PHP2 )
        dZOO = vol * ( B1*ZOO - B2*ZOO  - B3*ZOO - B4*ZOO2 )
        dDET = vol * ( A3*PHP2 + B3*ZOO + B4*ZOO2 - C1*DET - C2*DET )
        dDOP = vol * ( A1A2*PHP + C2*DET - D1*DOP )
        dDIP = vol * ( -A1*PHP + B2*ZOO  + C1*DET + D1*DOP )

c----------------------------------------------------------------  

        es(1) = el(1) + dPHP
        es(2) = el(2) + dZOO
        es(3) = el(3) + dDET
        es(4) = el(4) + dDOP
        es(5) = el(5) + dDIP

c----------------------------------------------------------------  

        end

c********************************************************************

	subroutine euler(dt,vol,c,cds)

	implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 5 )

	real dt			!time step [sec]
	real vol		!volume [m**3]
	real c(nstate)		!state variable [mg/L] == [g/m**3]
	real cds(nstate)	!source term [g/sec]

	integer i

	do i=1,nstate
	  c(i) = ( vol * c(i) + dt * cds(i) ) / vol
	end do

	end
	
c********************************************************************

