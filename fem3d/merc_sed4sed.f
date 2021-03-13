!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Ginevra Rosat
!    Copyright (C) 2018,2020  Donata Melaku Canu
!    Copyright (C) 2019  Celia Laurent
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
!
! sediment dynamics in the seabed for mercury subroutines (sed4merc_sed.f)
! the sediment bed is initialized in mercury.f (L245 - L257)
!
! revision log :
!
! 30.05.2018	dmc&gir	integration of the 0d module into SHYFEM
! 15.07.2019	gir&clr	sediment active layer added (Carniello et al., 2012)
! 24.02.2020    dmc	read variable tCE(k) according to the element type
!
! notes :
!
! dmc= donata canu; gir= ginevra rosati; clr= célia laurent
!
!*****************************************************************

        subroutine sed4merc_sed(k,dt,area,C,wat_vol,
     +                 taub,Cs,Dssink_sum,Dpsink_sum,
     +                 Sres, Pres,Vr,Bvels,Bvelp,
     +                 ds_gm2s, dp_gm2s,tCE,        
     +                 dZbedk,dZactivk)         


        implicit none

      integer nstate
      parameter (nstate=2)
      real C(nstate)         !variables in water: C(1):Silt, C(2):POM          [mg/L]
      real CD(nstate)        !derivatives in water CD(1) Dsilt, CD(2) DPOM     [g/d]
      real CS(nstate)        !sed  variables in sediment: C(1):Silt, C(2):POM  [g/m3]
      real CDS(nstate)       !deriv in sediment CD(1) Dsilt, CD(2) DPOM
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
      real dme, logdme, es,dme2     !terms for computing resuspension multiplier dme2           
      real Vr, Vrs, Vrp             !resuspension velocity                                    [m s-1]
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
      real dZbed
      real wdep
	real, save :: esmax = 0.
	real, save :: esmax0 = 4.

      integer m              !time indicator for write outputs, debug
      integer k              !node     
      integer ipext          !node external number
      logical constant_parameters
	integer nits,nita

	call get_time_iterations(nits,nita)

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

        sed_vol_old=area*dZactivk   !GINEVRA CHECK   ???

	if( c(1) < 0 ) write(6,*) '***** 1 ggu ',c
	if( c(1) < 0 ) call mdebug('ggu 1',k,2,c)

        sw=C(1)      ! [mg/l]
        POMw=C(2)

        silt=CS(1)   ! [g/m3]
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
       por = 1.0 - (DryD/Pdens)        ! [L(w)/L(s+w)]
       BulkD = por + DryD              ! [g/cm3] 
c
c____________ Input critical shear for erosion tCE ____________
c
c       tCE  = 0.65                     ![Pa]

c        write(*,*) tCE, 'tCE read from mercury.f'
        

c___________ Compute resuspension multiplier(dme/dt) from BulkD _________  
c___________ from Hwang and Metha (1989) in Tetra Tech (2002) _________         
c

	!write(666,*) aline,'  ',dtime,'  -----------------------'
	!write(666,*) Pdens,por,p_POM
	!write(666,*) por,DryD
	!write(666,*) BulkD,por,DryD
	!call flush(666)

        es = 0.198/(BulkD-1.0023)
        logdme = 0.23*exp(es)

	call massert('es',k,es,0.,5.)

        !if( es > 1.1*esmax ) then
        !  call mdebug('es>esmax',k,2,(/es,logdme/))
        !  esmax = es
        !end if
        if( es > esmax0 ) then
          esmax = max(es,esmax)
          !logdme = 0.23*exp(es)
          !call mdebug('es>esmax0',k,2,(/es,logdme/))
          call mdebug('es adjusted',k,3,(/es,esmax/))
          es = esmax0
          logdme = 0.23*exp(es)
        end if

        dme = 10.**logdme                      ![mg cm-2 hr-1]
        dme2 = ((dme/1000.)/3600.)*10.**4.     ![g m-2 s-1]       
c
c       write(665,*) k,es,logdme,BulkD,OC_mg_g
c       write(666,*) k,taub,dme2, por 
c       write(667,*) k,dme2, OC_mg_g 

c___________ Resuspension Occurrence _________________________________ 

        if (taub>1) then 
         taub=1
        end if  
 
c        write(*,*) '---------------------'
c        write(*,*) k,'nodo', taub, 'taub'
c        write(*,*) '---------------------'

        if (taub > tCE) then   
         MagR = (taub/tCE - 1.)       
        else                             
         MagR = 0.
        end if

	if( c(1) < 0 ) write(6,*) '***** 2 ggu ',c
	if( c(1) < 0 ) call mdebug('ggu 2',k,2,c)
c___________ Compute Resuspension Fluxes and rates _________________________ 

       Er = dme2 * MagR                  ![g m-2 s-1]
       Vr = Er/(POM+silt)                ![g m-2 s-1]*[m3 g]-->  [m s-1]

       Vrs=Er/silt
       Vrp=Er/POM

       Sres  = (Er*area)*(p_silt/100.)  ![g s-1] = [g m-2 s-1] * [m2]     
       Pres  = (Er*area)*(p_POM/100.)
       Tres  =  Sres + Pres  
c  
       rs_gm2s = Er*(p_silt/100)
       rp_gm2s = Er*(p_POM/100)
       res_gm2s= Tres/area         
     
c       if (ipext(k) .EQ. 2284) then 
c       write(888,*) Vr, Er, Sres, Pres, taub
c       end if 

c       if (taub .GT. 0.4) then
c       write (779,*) ipext(k),Vr, Er, Sres, Pres, taub
c       end if
c___________ Compute Burial Velocities and Fluxes ____________________________
       
       Bvels = (ds_gm2s-rs_gm2s)/silt    ![m s-1]
       Bvelp = (dp_gm2s-rp_gm2s)/POM 
c        
c       write(*,*) 'dep', ds_gm2s,'res', rs_gm2s, k       

c       Bvel = (dep_gm2s-res_gm2s)/(silt+POM) ![[m s-1]
c       
       Bsflux = Bvels*silt*area    ![m s-1]*[g m-3]*[m2] = [g s-1]
       Bpflux = Bvelp*POM*area 

       siltm = silt*sed_vol_old    ![g]   
       POMm = POM*sed_vol_old                                 
c       

       if (ipext(k) .EQ. 1372) then
       write(665,*) Bvels,Bvelp,Bsflux,Bpflux, 'Sed_4Merc_SED' 
       write(664,*) silt, POM, 'sed4Merc_SED'
c       write(486,*) Dssink_sum, Dpsink_sum, 'sed4Merc_Sed'
       end if 
      
        CS(1)=silt            ![g/m3]
        CS(2)=POM             ![g/m3]

        C(1) = sw
        C(2) = POMw
	if( c(1) < 0 ) write(6,*) '***** 3 ggu ',c
	if( c(1) < 0 ) call mdebug('ggu 3',k,2,c)

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
c________________________________________________________________________________    
c       
c____________Positive burial push sediment below the active layer__________________
c____________Negative burial is equal to net resuspension thus not included________
        
	!if( nits >= 6498 .and. nits < 6530 ) then
	if( abs(Sres) > 1.e+6 ) then
	  call mdebug('burial',k,3,(/Dssink_sum,Sres,Bsflux/))
	  Sres = 0.
	  Bsflux = 0.
	end if
	if( abs(Pres) > 1.e+6 ) then
	  call mdebug('burial',k,3,(/Dpsink_sum,Pres,Bpflux/))
	  Pres = 0.
	  Bpflux = 0.
	end if

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
        
c        write (988,*) (CS(m), m=1,nstate),'k',k,'SEdSolid var beforer'
c        write (987,*) (C(m), m=1,nstate),'k',k, 'SEdwater var beforer' 

       wdep=wat_vol/area


       if (cs(1) .LE. 0.0) then
       write(*,*) 'SiltSed<0',cs(1),'s4m_sedaft kext=',ipext(k)
       stop
       else if (cs(2) .LE. 0.0) then
       write(*,*) 'POMwSed<0',cs(2),'s4m_sedaft kext=',ipext(k) 
       stop
       end if

	if( c(1) < 0 ) write(6,*) '***** 4 ggu ',c
	if( c(1) < 0 ) call mdebug('ggu 4',k,2,c)
        call merc_euler_sed(2,dt,sed_vol_old,sed_vol_new,cs,cold,cds) 
        !claurent-OGS: vol old and new are different for sea bed sediment 

       if (cs(1) .LE. 0.0) then
       write(*,*) 'SiltSed<0',cs(1),'s4m_sedaft kext=',ipext(k)
       stop
       else if (cs(2) .LE. 0.0) then
       write(*,*) 'POMSed<0',cs(2),'s4m_sedaft kext=',ipext(k) 
       stop
       end if

       if (ipext(k) .EQ. 1372) then
         write(999,*) ipext(k), Sres,Dssink_sum,Bvels
       elseif(ipext(k) .EQ. 3057) then
         write(998,*) ipext(k),Sres,Dssink_sum,Bvels
      elseif(ipext(k) .EQ. 2654) then
         write(997,*) ipext(k),Sres,Dssink_sum,Bvels
       end if


c        call merc_euler (2,dt,sed_vol_old,cs,cold,cds)

c         call merc_euler_sed(2,dt,wat_vol+dZbedk*area,
c     +                        wat_vol+(dZbedk+dZit)*area, c,cold,cd)  
c    

       !write(222,*) 'Vold', sed_vol_old, 'Vnew', sed_vol_new, k   

	if( c(1) < 0 ) write(6,*) '***** 5 ggu ',c
	if( c(1) < 0 ) call mdebug('ggu 5',k,2,c)
       if (c(1) .LE. 0.0) then
       write(*,*) 'Siltw<0',c(1),wdep,'s4m_sedbef kext=',ipext(k)
c       write(443,*) 'Siltw<0',c(1), 's4m_sedbef kext=',ipext(k)
       else if (c(2) .LE. 0.0) then
c       write(447,*) 'POMw<0',c(2), 's4m_sedbef kext=',ipext(k) 
       write(*,*) 'POMw<0',c(2),wdep, 's4m_sedaft kext=',ipext(k)
	stop 'error stop... conz<0 (ggu)'
       end if

c      call merc_euler (2,dt,
c     +                        wat_vol, c,cold,cd)

         call merc_euler_sed(2,dt,wat_vol+dZbedk*area,
     +                        wat_vol+(dZbedk+dZit)*area, c,cold,cd)  
    

        dZbed  = dZbedk+dZit  ! Celia: update sea bed depth

	if( c(1) < 0 ) write(6,*) '***** 6 ggu ',c
	if( c(1) < 0 ) call mdebug('ggu 6a',k,2,c)
	if( c(2) < 0 ) call mdebug('ggu 6b',k,2,c)

       if (c(1) .LE. 0.0) then
       write(*,*) 'Siltw<0',c(1),wdep,'s4m_sedaft kext=',ipext(k)
c       write(443,*) 'Siltw<0',c(1), 's4m_sedaft kext=',ipext(k)
       else if (c(2) .LE. 0.0) then
c       write(447,*) 'POMw<0',c(2), 's4m_sedaft kext=',ipext(k)       
       write(*,*) 'POMw<0',c(2),wdep,'s4m_sedaft kext=',ipext(k) 
      end if

        c(1)=max(0.001,c(1))  
        c(2)=max(0.001,c(2))  ! Celia: allow POM_w to be NULL (for drying nodes)

        ! Celia: if this trick is kept, unsinked silt and pom should be
        ! substracted from the bottom sediment (doing an inversed euler)
        ! in order to maintain mass conservation

c        write (988,*) (CS(m), m=1,nstate),'k',k,'SEdSolid var after' !' Transformations'
c        write (987,*) (C(m), m=1,nstate),'k',k,'SEdwater var after' 
     
c        write (777,*) wat_vol,cd,k, 'wat_vol and water derivative cd'
c        write (888,*) wat_vol,cds,k, 'wat_vol and sed derivative cds'

        silt   = CS(1)  
        POM    = CS(2)

        sw = C(1)
        POMw=C(2)

       if (POMw .LE. 0.0) then 
       write(*,*), 'instability - negative POM_w in sed4merc_sed kext=',
     +   ipext(k)
       stop
       else if (sw .LE. 0.0) then  
       write(*,*), 'instability - negative slt_w in sed4merc_sed kext=',
     +   ipext(k)
       stop 
       end if 
       end

!*****************************************************************

	subroutine mdebug(text,k,n,f)

	implicit none

	character*(*) text
	integer k
	integer n
	real f(n)

	integer iunit
	double precision dtime
	character*20 aline

	integer ipext

	iunit = 777

	call get_act_dtime(dtime)
	call get_timeline(dtime,aline)

	write(iunit,*) '------------------------------',k,ipext(k)
	write(iunit,*) trim(text),dtime,'  ',aline
	write(iunit,*) n,f
	call flush(iunit)

	end

!*****************************************************************

	subroutine massert(text,k,c,cmin,cmax)

	implicit none

	character*(*) text
	integer k
	real c,cmin,cmax

	double precision dtime
	character*20 aline

	integer ipext

	if( c >= cmin .and. c <= cmax ) return

	call get_act_dtime(dtime)
	call get_timeline(dtime,aline)

	write(6,*) 'assertion violated: ',trim(text)
	write(6,*) k,ipext(k)
	write(6,*) dtime,'  ',aline
	write(6,*) c,cmin,cmax

	stop 'error stop massert: assertion violated'

	end

!*****************************************************************

