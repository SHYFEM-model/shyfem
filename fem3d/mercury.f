
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016  Georg Umgiesser
!    Copyright (C) 2017-2018,2020  Donata Melaku Canu
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

c mercury routines
c
c contents :
c
c revision log :
c
c 15.05.2016    ggu     started mercury from bio3d
c 17.05.2017    dmc     updated with merc_water
c 17.05.2017    dmc     and volatilization routine wind dependent
c 02.08.2018    dmc     add solids in water and at the bottom
c 05.09.2018    dmc	Hg0atm è immesso cost in mercury_react,
c ...                   da inserire come file esterno attraverso main
c 24.02.2020    dmc     tcek reads from subroutine init_crit_thre_erosio
c ...                   according to element type--> node
c 22.11.2020    ggu     some more on restart
c
c
c********************************************************************
c
c notes :
c
c State variables used: (mercury) 
c
c State variables used: (mercury)
c Hg0           81      1
c Hg2           82      2
c Hg3           83      3
c msed1         84      4
c msed2         85      5 
c
c********************************************************************

!====================================================================
      module mercury
!====================================================================

      implicit none

      integer, parameter :: npstate = 3   !pelagic state variables
      integer, parameter :: nsstate = 2 !bottom state variables

      integer, parameter :: nsolwst = 2  !solids in w state variables
      integer, parameter :: nsolsst = 2  !solids at bottom state variables

      real, save, allocatable :: emp(:,:,:) !Hg pelagic state vector
      real, save, allocatable :: ems(:,:)   !Hg sediment state vector

      real, save, allocatable :: emsolw(:,:,:) !solids in water state vector
      real, save, allocatable :: emsols(:,:)   !solids at bottom state vector

      real, save, allocatable :: etau(:)    ! bottom friction vector           !claurent-OGS: created to produce outputs fields
      real, save, allocatable :: dZbed(:)   ! variations of the sea bed depth  !claurent-OGS: created to produce outputs fields
                                            ! due to erosion or deposition     !claurent-OGS: created to produce outputs fields 
      real, save, allocatable :: dZactiv(:) ! thickness of the active layer    !claurent-OGS: created to produce outputs fields

!       real, save, allocatable :: eload(:,:,:)   !atmospheric loading

	integer, save :: ia_out(4)
	double precision, save :: da_out(4)

	integer, save :: iubp,iubs,iubsolw,iubsols

!====================================================================
	end module mercury
!====================================================================

        subroutine mercury_module

c general interface to mercury module

        implicit none

	include 'femtime.h'

        real dt


	dt = dt_act
        call mercury3d(dt)

        end

c********************************************************************

	subroutine mercury3d(dt)

c eco-model cosimo

	use mod_diff_visc_fric
	use levels
	use basin
	use mercury
 
	implicit none

!	integer it	!time in seconds
	real dt		!time step in seconds


        character*10 what
        character*10 what2

	integer k,i,l,lmax
	integer imerc
        integer id,idc
	integer nintp,ivar
	integer nbc
	real t,s
	real u,v

	real epela(npstate)
	real esedi(nsstate)


	real esolw(nsolwst)
	real esols(nsolsst)

        real epload(npstate) !loading for each node
        real loadsup(npstate) 

! next array specifies boundary conditions if non are given
! mercury concentration are in mg/L equivalent to  g/m3

      real, save :: epbound(npstate) = (/0.0,0.0,0.0/)  !default bound cond. Hg in water
      real, save :: epinit(npstate) = (/0.08,5.9,0.02/) !default in. cond. Hg0,HgII,MeHg in water
      real, save :: esinit(nsstate) = (/3.0,0.06/)       !default in. cond. HgII, MeHg in sediment [mg/kg]
      real, save :: eploadin(npstate) = (/0.0,0.01,0.0/) !atm load kg/day

      real, save :: esolbound(nsolwst) = (/5.0,0.1/)   !default bound cond.solids in water  !grosati-OGS:calibration 
      real, save :: esolwinit(nsolwst) = (/3.0,1./)       !default in. cond. solids in water
      real, save :: esolsinit(nsolsst) = (/3.,0./)       !initial OC%, dummy var.  
c       esolsinit: initial value is the OC% in sediment. From this value
c       we compute the % of POM and weighted particle density

        real tpstot(npstate)              !for mass test
        real tsstot(nsstate)
        real tsoltot(nsolwst)
        real tsolstot(nsolsst)

	integer, save, allocatable :: idmerc(:)
	integer, save, allocatable :: ids4merc(:)

	logical bsurf,bbottom
	integer iunit
	integer j
        real qrad
        real dtday
        real area,vol
        real getpar
        logical has_output,next_output
        logical has_output_d,next_output_d

        integer mode
        real ai,lsurf

        logical boxtype
        logical bcheck
        logical bresi,breact,bdecay
        integer ie,ii
        integer kspec
        integer nvar
        double precision dtime0,dtime
        real d          !element tickness 
        real depth      !element depth
        !real cbod,nh3,krear,sod
        real vel
        real windspeed,tempair
        real uws       !qss,tas,urs,tbs,ccs,ps !meteo_get_values
        real wx,wy        !get_wind
        real tday,t0,tsec
        real stp
        real mass
        real wsink
        real flux,conz(nlvdi)
        real conz1,conz2
        real Dpsink, Dssink
        real Dpsink_sum, Dssink_sum   !claurent-OGS: stores sum of sinks for multiple water levels above sea bed 

        logical t1
        real t2,t3,t4
        integer t5
        real t6, t7,t8,t9
        real t10(2)
        real t11,t12
        real vds,vdp

        real taubot(nkn)        !ottom stress from subssed.f simple_sed
        real tcek(nkn)          !node crit. threshold for erosion [N/m**2]

        real tau
        real Sres, Pres

        real Shgsil, Shgpom, Smhgsil, Smhgpom ! Deposition rates 
        real faq1,faq2,fdoc1,fdoc2      ! frazioni di merc in water

        real Hg2sed,MeHgsed,siltin,POMin
        real hgp1,hgp2,hgp,mehgp,por,k1tp,k2tp
        real mehgp1,mehgp2,hgd,mehgd,hgit,mehgt

        real silt,pom,vr,bvels,bvelp    !FIXME VARIABILI TEMPORANEE 
        real p_poc,p_pom,p_silt,OC_mg_g,pdens,DryD
        real ds_gm2s, dp_gm2s   !claurent-OGS: transfer of values between sed4merc routines to avoid double calculation 
        real Rhgsil,Rhgpom,Rmhgsil,Rmhgpom 

        integer nbnds

        real, save :: rkpar,difmol
        integer, save :: icall = 0
        integer fortfilenum
        integer itype
         

        include 'femtime.h'

c------------------------------------------------------------------
c	initial and boundary conditions  [mg/l]			??
c	initial loadings                 [g/m**2/sec]		??
c------------------------------------------------------------------

	breact = .true.		!use reactor

        what = 'mercury'
        what2 = 's4mercury'

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

	if( icall .le. -1 ) return

	if( icall .eq. 0 ) then
	  imerc = nint(getpar('imerc'))
	  if( imerc .le. 0 ) icall = -1
	  if( icall .le. -1 ) return
	  icall = 1

           call get_first_dtime(dtime0)
        dtime0=itanf

c       ___________________________________________________
c        initialization data transformed to variables units
c       ----------------------------------------------------
c       solids in sediment initialization

        p_poc=esolsinit(1)
        p_pom=p_poc*1.7
        p_silt=100.-p_pom
        OC_mg_g=10.0*p_poc
        pdens = ((1.25*p_POM)+(2.65*(100.0 - p_POM)))/100.0
        DryD=1.776-0.363*log(OC_mg_g)

        por=1.-(DryD/pdens)
        siltin=DryD*10.**6*p_silt/100.
        POMin=DryD*10.**6*p_pom/100.

        write(*,*) 'siltin', siltin, 'DryD', DryD

        esolsinit(1)=siltin     !GINEVRA FIXME, unità di misura
        esolsinit(2)=POMin

       write(*,*) 'esols', esolsinit

c       mercury in sediment initialization

        k1tp=10.**5. !FIXME parametro che deriva dalla kd, mettere 
c                   una routine unica per settare tutti i parametri        
        k2tp=14590.  !FIXME

        Hg2sed=esinit(1)   !ug(hg)/g(sed)
        MeHgsed=esinit(2)
        
        hgp1   = Hg2sed  * siltin  ! mg(hg)/kg(s)]
        hgp2   = Hg2sed  * POMin  
        mehgp1 = MeHgsed * siltin
        mehgp2 = MeHgsed * POMin
        
        hgp   = hgp1 + hgp2 
        mehgp = mehgp1 + mehgp2

        hgd   = (Hg2sed *(por/k1tp))*10.**6. ![ug(hg) m-3(w+s)]
        mehgd = (MeHgsed*(por/k2tp))*10.**6.


c ------------- total hg and mehg -----------------------

      hgit =  hgp + hgd      ! [ug m-3]or [ng(hg) l(w+s)-1]  Hg in particulate AND dissolved
      mehgt = mehgd + mehgp  ! [ug m-3]or [ng(hg) l(w+s)-1] MeHg in particulate AND dissolved
        
c      write(*,*) 'hgit',hgit, k, 'mercury.f'  
c      write(*,*) 'mehgt',mehgt,  'mercury.f'         
     
 
        esinit(1)=hgit
        esinit(2)=mehgt

c         --------------------------------------------------
c	  initialize state variables
c         --------------------------------------------------

	  allocate(emp(nlvdi,nkndi,npstate))
	  allocate(ems(nkndi,nsstate))
	  allocate(emsolw(nlvdi,nkndi,nsolwst))
	  allocate(emsols(nkndi,nsolsst))
	  allocate(etau(nkndi))      !claurent-OGS: created to produce outputs fields
	  allocate(dZbed(nkndi))     !claurent-OGS: created to produce outputs fields
	  allocate(dZactiv(nkndi))   !claurent-OGS: created to produce outputs fields 


           do i=1,npstate
            emp(:,:,i) = epinit(i)
          end do

           do i=1,nsstate
            ems(:,i) = esinit(i)
          end do

           do i=1,nsolwst
            emsolw(:,:,i) = esolwinit(i)
          end do

           do i=1,nsolsst
            emsols(:,i) = esolsinit(i)
          end do
          
          etau(:)=0.0                  !claurent-OGS: created to produce outputs fields
          dZbed(:)=0.0    ! [meters]   !claurent-OGS: created to produce outputs fields
          dZactiv(:)=0.05 ! [meters]   !claurent-OGS: created to produce outputs fields 
c         --------------------------------------------------
c	  initial conditions from file (only for pelagic part)
c         --------------------------------------------------

c this is still not working FIXME
c	  nvar = npstate
c	  call mercury_init_file(dtime0,nvar,nlvdi,nlv,nkn,epinit,emp)
c	  nvar = nsolw
c	  call mercury_init_file(dtime0,nvar,nlvdi,nlv,nkn,esolwinit,esolw)

c         --------------------------------------------------
c	  set boundary conditions for all state variables
c         --------------------------------------------------

          nbc = nbnds()
          allocate(idmerc(nbc))
          idmerc = 0
	  call get_first_dtime(dtime0)
          nintp = 2
	  nvar = npstate
        !write(6,*) 'npstate', nvar,epbound,idmerc
          call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +                          ,epbound,idmerc)

          nbc = nbnds()
          allocate(ids4merc(nbc))
          ids4merc = 0
	  call get_first_dtime(dtime0)
          nintp = 2
	  nvar = nsolwst
          call bnds_init_new(what2,dtime0,nintp,nvar,nkn,nlv
     +                          ,esolbound,ids4merc)
c         --------------------------------------------------
c	  initialize eco model
c         --------------------------------------------------

	  call mercury_init

c         --------------------------------------------------
c	  parameters for transport/diffusion resolution
c         --------------------------------------------------

          rkpar=getpar('chpar')
          difmol=getpar('difmol')

c         --------------------------------------------------
c	  initialize output 
c         --------------------------------------------------

	  call mercury_init_file_output
	  call mercury_write_file_output(dtime0)

	  write(6,*) 'mercury model initialized...'

	end if

c-------------------------------------------------------------------
c normal call
c-------------------------------------------------------------------

	wsink = 0.
        Shgsil=0
        Shgpom=0
        Smhgsil=0 
        Smhgpom=0
        Vdp=0
        Vds=0
        ds_gm2s=0
        dp_gm2s=0

c       FIX ME
        Sres=0
        Pres=0
        Vr=0
        Bvels=0
        Bvelp=0
        Rhgsil=0
        Rhgpom=0
        Rmhgsil=0
        Rmhgpom=0        
c------------------------------------------------------------------
c       compute loading/unit surface	!FIXME -> this could be done only once
c------------------------------------------------------------------

	call mercury_surface_loading(npstate,eploadin,loadsup)

c-------------------------------------------------------------------
c       time management
c-------------------------------------------------------------------

        t0 = 0.
        dtday = dt / 86400.

        tsec = it
        tday = it / 86400. + t0         !time in days, FEM 0 is day t0


       if( it .le.dtime0+dt ) then
         do fortfilenum=250,282
            write(fortfilenum,
     +       "(2(a10,','),4(a15,','))")
     +       'it','kext','wdepth','silt_w','POM_w','tauB'
         end do
       endif


       if( it .le.dtime0+dt ) then
         do fortfilenum=350,382
            write(fortfilenum,
     +       "(2(a10,','),4(a15,','))")
     +       'it','kext','wdepth','Hgod','Hg2','MeHg'
         end do
       endif

c	-------------------------------------------------------------------
c	loop on elements for reactor
c	-------------------------------------------------------------------

        mode = +1               !new time level for volume and depth

	if( breact ) then	!use reactor ?

        !call simple_sedi_bottom_stress(taubot)  !claurent-OGS: get friction coefficient array only once for all nodes
        call bottom_stress(taubot)  !claurent-OGS: get friction coefficient array only once for all nodes
c       questa call è nel loop, il taub viene letto ad ogni passo.
c       mettere fuori? FIXME
        
        call init_crit_thre_erosion(tcek)

	do k=1,nkn		!loop on nodes

          lmax = ilhkv(k)
	  call get_light(k,qrad)
          tau=taubot(k)    !claurent-OGS: get friction coeff at current k node
          etau(k)=tau      !claurent-OGS: store value to write 2Dfield in output
c       FIXME conz è letta fuori dal ciclo sui livelli, è la conz(lmax)
 
          call get_wind(k,wx,wy)
          uws=(wx*wx+wy*wy)**(1./2.)
                

          depth=0               !reinitialize the depth at each node
          Dssink_sum=0.0   !claurent-OGS
          Dpsink_sum=0.0   !claurent-OGS
          do l=1,lmax
        
            call dvanode(l,k,mode,d,vol,area)   !gets depth, volume and area
            call getts(l,k,t,s)                 !gets temp and salt
            call getuv(l,k,u,v)                 !gets velocities u/v
            vel = sqrt(u*u+v*v)

            id = 1000*k+l
	    bsurf = (l==1)
	    bbottom = (l==lmax)
                
            depth=depth+d   !compute the depth at the centre of the element

            !write(6,*) bsurf,bbottom, l,d,k,depth,area,vol, 'mercury'
            !write(6,*) epload(:),'epload' !
            !write(6,*) esolw,k, 'flux and conz,k' !

            epload(:) = loadsup(:)*area         !(g/m2)*m2 area of element
            epela(:) = emp(l,k,:)
            esolw(:) = emsolw(l,k,:)

            boxtype=.true.                
c      FIXME conz(l) da leggere nel ciclo sui livelli


        call sed4merc_water(bbottom,dtday,tday,vol,d,k,t,s,tau
     +                          ,area,esolw,
     +                          Dssink,Dpsink,Vds,Vdp,   
     +                          ds_gm2s, dp_gm2s) !claurent-OGS:get values then send them to sed4merc_sed 
           
              Dssink_sum=Dssink_sum+Dssink !claurent-OGS: stores sum of sinks for multiple ... 
              Dpsink_sum=Dpsink_sum+Dpsink !claurent-OGS: ... water levels above sea bed 
           
c        write(86,*) Dssink, Dpsink, 'dssink and dpsink in mercury.f'
                emsolw(l,k,:) = esolw(:)

               conz1=esolw(1)
               conz2=esolw(2)
      call mercury_react(id,bsurf,bbottom,boxtype,dtday,vol
     +                  ,d,k,t,uws,area,s,qrad,epela,epload
     +                  ,Vds,Vdp,conz1,conz2,             !tday,
     +                  Shgsil,Shgpom,Smhgsil,Smhgpom,
     +                  faq1,faq2,fdoc1,fdoc2)

              emp(l,k,:) = epela(:)
           
          
        if (bbottom) then
        
          esedi(:)=ems(k,:)
          esols(:)=emsols(k,:)
        
c         write(6,*) esols, 'esols' 
          call sed4merc_sed(k,dtday,area,esolw,vol,
     +                         tau,esols,Dssink_sum,Dpsink_sum,  !claurent-OGS: send sum instead of single sinks
     +                           Sres,Pres,Vr,Bvels,Bvelp,
     +                        ds_gm2s, dp_gm2s,tcek(k),       !claurent-OGS: values required by sed4merc_sed
     +                        dZbed(k),dZactiv(k))    !claurent-OGS: get thicknesses for extraction of the fields in output                
          
          emsolw(l,k,:)=esolw(:)
          emsols(k,:)=esols(:)
        
      
          silt=esols(1)
          pom= esols(2)

c               write(*,*) 'silt_dopo_sed4merc', silt
         
          epload(:)=0
          esedi(:)=ems(k,:)          !Hg in sed
          epela(:)=emp(l,k,:)        !Hg in water
        
!         write(*,*) 'ems_before',ems(k,:)
!         write(*,*) 'esedi_before',esedi(:)
 
          call mercury_sed_react(dtday,
     +                         k,t,area,esedi,epela,
     +                  Shgsil, Shgpom, Smhgsil, Smhgpom,
     +             faq1,faq2,fdoc1,fdoc2,
     +             silt,pom,Vr,Bvels,Bvelp)
         

          ems(k,:) = esedi(:)
          emp(l,k,:) = epela(:)
 
c         write(*,*) 'silt_dopo_mercury_sed', silt
c         write(*,*) 'ems_after',ems(k,:)
c         write(*,*) 'esedi_after',esedi(:)

        end if ! bbottom
c           end if ! d>0.01
          end do
	end do

	end if	!breact
c	-------------------------------------------------------------------
c	advection and diffusion
c	-------------------------------------------------------------------

	dtime = it
	call bnds_read_new(what,idmerc,dtime)

!$OMP PARALLEL PRIVATE(i)
!$OMP DO SCHEDULE(DYNAMIC)	

	do i=1,npstate

          call scal_adv(what,i
     +                          ,emp(1,1,i),idmerc
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

          call tsmass (emp(1,1,i),1,nlvdi,tpstot(i)) !mass control

	end do
!$OMP END DO NOWAIT   ! claurent-OGS: added end of parallel region as the PARALLEL...
!$OMP END PARALLEL    ! claurent-OGS: ...  region must be associated to only one loop		

	call bnds_read_new(what2,ids4merc,dtime)   ! claurent-OGS: read boundary conditions for solids 

!$OMP PARALLEL PRIVATE(i)  ! claurent-OGS: added new parallel region as the PARALLEL...
!$OMP DO SCHEDULE(DYNAMIC) ! claurent-OGS: ...  region must be associated to only one loop		     
	do i=1,nsolwst

          call scal_adv(what2,i                          ! claurent-OGS: sends solid boundary conditions 
     +                          ,emsolw(1,1,i),ids4merc  ! claurent-OGS: sends solid boundary conditions
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

          call tsmass (emsolw(1,1,i),1,nlvdi,tpstot(i)) !mass control

	end do

!$OMP END DO NOWAIT	
!$OMP END PARALLEL	

        do i=1,nsstate
          call scalmass(ems(1,i),0.1,tsstot(i))   !mass ctrl sed
        end do

        do i=1,nsolsst
          call scalmass(emsols(1,i),0.1,tsstot(i))   !mass ctrl sed
        end do

c 6.8. dmc fin qui

c	-------------------------------------------------------------------
c	write of results
c	-------------------------------------------------------------------

	call mercury_write_file_output(dtime)

c	-------------------------------------------------------------------
c	debug output
c	-------------------------------------------------------------------

c	-------------------------------------------------------------------
c	end of routine
c	-------------------------------------------------------------------

	end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine mercury_init

! initializes mercury routines

	implicit none

	end

c*************************************************************

	subroutine mercury_init_file(dtime,nvar,nlvddi,nlv,nkn,val0,val)

c initialization of mercury from file

        implicit none

        double precision dtime
        integer nvar    !npstate+nsstate, mercury in water & sed
        integer nlvddi
        integer nlv
        integer nkn
        real val0(nvar)
        real val(nlvddi,nkn,nvar)

        call tracer_file_init('mercury init','mercin',dtime
     +                          ,nvar,nlvddi,nlv,nkn,val0,val)

        end

c*************************************************************

        subroutine mercury_init_file_output

        use basin
        use levels
        use mercury

        implicit none

        integer ishyff,nvar,id
        logical has_output,has_output_d
        real getpar

        ishyff = nint(getpar('ishyff'))

        call init_output('itmcon','idtcon',ia_out)
        if( ishyff == 1 ) ia_out = 0
        if( has_output(ia_out) ) then
          call open_scalar_file(ia_out,nlv,npstate,'mer')
          iubp = ia_out(4)
          call open_scalar_file(ia_out,1,nsstate,'mes')
          iubs = ia_out(4)
        end if

        call init_output_d('itmcon','idtcon',da_out)
        if( ishyff == 0 ) da_out = 0
        if( has_output_d(da_out) ) then
          nvar = npstate + nsstate+nsolwst+nsolsst
          call shyfem_init_scalar_file('merc',nvar,.false.,id)
          da_out(4) = id
        end if

        write(6,*) 'merc init file output done'

        end 

c*************************************************************

	subroutine mercury_write_file_output(dtime)

	use basin
	use levels
	use mercury

	implicit none

	double precision dtime

	integer nvar,id,idc,i
	logical next_output,next_output_d

        if( next_output(ia_out) ) then
	  ia_out(4) = iubp
	  do i=1,npstate
	    idc = 250 + i
	    call write_scalar_file(ia_out,idc,nlvdi,emp(1,1,i))
	  end do

	  ia_out(4) = iubs
	  do i=1,nsstate
	    idc = 270 + i
	    call write_scalar_file(ia_out,idc,1,ems(1,i))
	  end do

	  ia_out(4) = iubsolw
	  do i=1,nsolwst
	    idc = 280 + i
	    call write_scalar_file(ia_out,idc,1,emsolw(1,1,i))
	  end do

	  ia_out(4) = iubsols
	  do i=1,nsolsst
	    idc = 290 + i
	    call write_scalar_file(ia_out,idc,1,emsols(1,i))
	  end do

	  ia_out(4) = 0                                    !claurent-OGS: created to produce outputs fields
          idc = 297                                        !claurent-OGS: created to produce outputs fields
          call write_scalar_file(ia_out,idc,1,dZbed(1))    !claurent-OGS: created to produce outputs fields 
	  ia_out(4) = 0                                    !claurent-OGS: created to produce outputs fields    
          idc = 298                                        !claurent-OGS: created to produce outputs fields
          call write_scalar_file(ia_out,idc,1,dZactiv(1))  !claurent-OGS: created to produce outputs fields 
	  ia_out(4) = 0                                    !claurent-OGS: created to produce outputs fields   
          idc = 299                                        !claurent-OGS: created to produce outputs fields
          call write_scalar_file(ia_out,idc,1,etau(1))     !claurent-OGS: created to produce outputs fields 
        end if                                             
                                                           
        if( next_output_d(da_out) ) then                   
          id = nint(da_out(4))
          do i=1,npstate
            idc = 250 + i
            call shy_write_scalar_record(id,dtime,idc,nlvdi
     +                                         ,emp(1,1,i))
          end do
          do i=1,nsstate
            idc = 270 + i
            call shy_write_scalar_record(id,dtime,idc,1
     +                                        ,ems(1,i))
          end do
          do i=1,nsolwst
            idc = 280 + i                                   ! claurent-OGS: was 270 + i
            call shy_write_scalar_record(id,dtime,idc,1
     +                                        ,emsolw(1,1,i))
          end do
          do i=1,nsolsst
            idc = 290 + i                                   ! claurent-OGS: was 270 + i
            call shy_write_scalar_record(id,dtime,idc,1
     +                                       ,emsols(1,i))
          end do

          idc = 297                                                !claurent-OGS: created to produce outputs fields
          call shy_write_scalar_record(id,dtime,idc,1,dZbed(1))    !claurent-OGS: created to produce outputs fields
          idc = 298                                                !claurent-OGS: created to produce outputs fields 
          call shy_write_scalar_record(id,dtime,idc,1,dZactiv(1))  !claurent-OGS: created to produce outputs fields
          idc = 299                                                !claurent-OGS: created to produce outputs fields
          call shy_write_scalar_record(id,dtime,idc,1,etau(1))     !claurent-OGS: created to produce outputs fields 

        end if

	end

c*************************************************************
c*************************************************************
c*************************************************************

       subroutine merc_euler(nstate,dt,vol,c,cold,cds) !

! new c is computed
! cold is returned which is just c before call

      implicit none

      integer nstate            !number of state variables
      real dt                  !time step [day]
      real vol            !volume [m**3]
      real volold,volnew      !volume [m**3]   ! claurent-OGS: differenciates volold and volnew
      real c(nstate)            !state variable [mg/L] == [g/m**3]
      real cold(nstate)      !old state variable (return)
      real cds(nstate)      !source term [g/day]

      integer i
      real mass            !mass [g]
      real mder            !derivative [g/day]

      volold = vol
      volnew = vol

      do i=1,nstate
        cold(i) = c(i)
        mass = c(i) * volold
        mder = cds(i)
        c(i) = ( mass + dt * mder ) / volnew
      !  write(88,*)dt,'dt-day'
        !write(88,*) mder,'mder'
        !write(88,*) cds(i),'cds'
        !write(88,*) i, c(i),nstate,'c(i)'
        !if(c(i).lt.0.00001) c(i)=0.00001
      end do
      end
     
c*************************************************************

      subroutine merc_euler_sed(nstate,dt,volold,volnew,c,cold,cds) 

claurent-OGS: differenciates volold and volnew
c      subroutine merc_euler(nstate,dt,vol,c,cold,cds) !

! new c is computed
! cold is returned which is just c before call

      implicit none

      integer nstate            !number of state variables
      real dt                  !time step [day]
      real vol            !volume [m**3]
      real volold,volnew      !volume [m**3]   ! claurent-OGS:
      real c(nstate)            !state variable [mg/L] == [g/m**3]
      real cold(nstate)      !old state variable (return)
      real cds(nstate)      !source term [g/day]

      integer i
      real mass            !mass [g]
      real mder            !derivative [g/day]

c      volold = vol
c      volnew = vol

      do i=1,nstate
        cold(i) = c(i)
        mass = c(i) * volold
        mder = cds(i)
        c(i) = ( mass + dt * mder ) / volnew
        !write(88,*)dt,'dt-day' 
        !write(88,*) 'old', volold, 'new', volnew, nstate
        !write(88,*) mder,'mder'
        !write(88,*) cds(i),'cds'
        !write(88,*) i, c(i),nstate,'c(i)'
        !if(c(i).lt.0.00001) c(i)=0.00001
      end do
      end
 
c*************************************************************
c*************************************************************
c*************************************************************

	subroutine mercury_surface_loading(npstate,eploadin,loadsup)

	use basin

	implicit none

	integer npstate
	real eploadin(npstate)
	real loadsup(npstate)

	integer mode
	integer k,l,i
	real d,vol,area,areatot

        loadsup=0
        mode=+1
        l=1

        areatot=0
        do k=1,nkn
          call dvanode(l,k,mode,d,vol,area)       !gets depth, volume, area
          areatot=areatot+area
        end do

        do i=1,npstate
          loadsup(i)=eploadin(i)*1000/areatot          !loading g/m2
        end do

        !write(6,*) loadsup(1),areatot,'loadsup1,areatot'
        !write(6,*) loadsup(2),areatot,'loadsup2,areatot'
        !write(6,*) loadsup(3),areatot,'loadsup3,areatot'

	end

c*************************************************************

        subroutine init_crit_thre_erosion(tcek)

        use basin

        implicit none

        integer ie, ii,k,ia
        real tce,tcek(nkn),tceaux
        
        tceaux=0.001
        tce=0.001

        do k=1,nkn
             tcek(k)=0.001
        end do 


        do ie=1,nel
          ia = iarv(ie)
          if( ia == 0 ) tce = 0.3      !FIXME
          if( ia== 1 )  tce = 0.3
          if( ia== 2 )  tce = 0.86
          if( ia== 3 )  tce = 0.86
          if( ia== 4 )  tce = 0.86
          if( ia== 5 )  tce = 0.86
          if( ia== 6 )  tce = 0.86
          if( ia== 7 )  tce = 0.86
          if( ia== 8 )  tce = 0.86
          if( ia== 9 )  tce = 0.86

c       mettere tutti gli if
          do ii=1,3
            k = nen3v(ii,ie) !3 nodi --> 1 elemento 
            tceaux=tcek(k)
c            tau = taubot(k)
c                val minimo              tcek(k)=minimo tra tce e tcek(k)
                tceaux=max(tceaux,tce)
                 tcek(k)=tceaux
c                 write(*,*) tce,k,tceaux,ia,ie
          end do
        end do

        end

!*************************************************************

        subroutine write_restart_mercury(iunit)
        implicit none
        integer iunit
	integer, save :: nvers = 1
	integer, save :: narray = 1
	write(iunit) nvers
	write(iunit) narray
        end

        subroutine skip_restart_mercury(iunit)
        implicit none
        integer iunit
	integer, save :: nvers = 1
	integer, save :: narray = 1
	read(iunit) nvers
	read(iunit) narray
	do i=1,narray
	  read(iunit)
	end do
        end

        subroutine read_restart_mercury(iunit)
        implicit none
        integer iunit
        end

!*************************************************************

