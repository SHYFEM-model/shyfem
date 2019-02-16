
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

! **************************************************
! -------------Fluid Mud Module -------------------
! **************************************************
!
! This routine manages the fluid mud computation.
!
! Written by : Christian Ferrarin, Deborah Bellafiore, Georg Umgiesser & Aron Roland
!
! Fundet by 03KIS065 by the Federal Ministry for Education and Research (BMBF) through the KFKI
! 
! Project Leader: Prof. Dr.-Ing. Prof. h.c. Ulrich Zanke
! Responsible   : Dr.-Ing. Aron Roland, IWW, Technische Universitat Darmstadt
!
!
        subroutine submud

	use mod_bound_geom
	use mod_sinking
	use mod_fluidmud
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use levels
	use basin

          implicit none

          include 'param.h'
  

        real difvm(0:nlvdim,nkndim)  !vertical diffusivity for transport of mud
        

        integer, save :: testnode, icycle
        logical, save :: ldebug,lsink,lhind,llambda,ladvlam
        logical, save :: lfloc,ladvfloc


        integer, allocatable :: innode(:)          !in boundary node number array (bound1)
        integer, allocatable :: ounode(:)          !out boundary node number array (bound2)

        integer ius,id,itmcon,idtcon	!output parameter

        integer kranf,krend,ibc
        integer kin,kou,nnode,knode,nk
          integer it                        !time in seconds
          real dt                           !time step in seconds
  

        integer iround
        real getpar

        save nk, innode, ounode


! -------------------------------------------------------------
! local mud variables
! -------------------------------------------------------------
        integer icall			!Initialization parameter
        integer imud                    !Fluid mud call parameter
        integer, parameter :: nsmud = 1                  !Number of grainsize classes
        integer nslam                   !Number of grainsize classes
        integer nsdmf                   !Number of grainsize classes
        integer nintp			!Interpolation type for boundary
        integer ivar			!This indicates only where to take boundary condition
        integer k,l,lmax		!Counters
        character*10 what
        double precision rhow		!Fluid density (kg/m3)
        double precision cgel(nlvdim),phigel(nlvdim)
        real mudref		!Initial fluid mud concentration (kg/m3)
        double precision visk		!molecular viscosity
        double precision phi(nlvdim),phip(nlvdim)
        integer lthick			!Initial fluid mud layer thickness (layer number)
        real mudhpar			!Fluid mud diffusion coefficient [m**2/s]
        real difmol			!Molecolar diffusion coefficient [m**2/s]
        real mwsink			!Fluid mud floc settling velocity [m/s]
	    real wsink,sinkaccel
	    real fact, dm1 
        real tsec			!Simulation time, real [s]
        real bnd3_mud(nb3dim,0:nbcdim)  !Array containing boundary state
        real bnd3_lam(nb3dim,0:nbcdim)  !Array containing boundary state
        real bnd3_dmf(nb3dim,0:nbcdim)  !Array containing boundary state
        real fmbnd(nsmud),fak,smooth
        real t_now,t_start		        !Boundary vector [kg/m3]
        logical, save :: circle,ldumpmud,linitmud,lsetbound

        save what,mudref,lthick,sinkaccel
        save fmbnd
        save bnd3_mud
        save mudhpar,difmol
        save mwsink,visk
        save ius,id,itmcon,idtcon
        save icall
        save t_now, t_start, smooth
        data icall /0/
        data t_now/0./

!       ----------------------------------------------------------
!       Initialization; This section is called only the first time step when ICALL = 0
!       ----------------------------------------------------------

        t_now = t_now + 1.

	call get_timestep(dt)
	call get_act_dtime(dtime)
	it = dtime

        imud = iround(getpar('imud'))
        rhow = getpar('rowass')

        if( imud .le. 0 ) then
          icall = -1
          return
        end if

	stop 'error stop submud: not adapted yet for new framework'

        if( icall .eq. 0 ) then

          itmcon = iround(getpar('itmcon'))
          idtcon = iround(getpar('idtcon'))

!         --------------------------------------------------
!         Initializes state variables and constants
!         --------------------------------------------------

          what     = 'mud'
          nslam    = 1
          nsdmf    = 1
          linitmud = getpar('linitmud')
          t_start  = getpar('t_start')
          smooth   = getpar('smooth')
          ldumpmud = getpar('ldumpmud')
          testnode = getpar('testnode')
          ldebug   = getpar('ldebug') 
          circle   = getpar('lcircle')
          mudref   = getpar('mudref')
          mudhpar  = getpar('mudhpar')
          difmol   = getpar('difmol')
          icycle   = getpar('icycle')
          mwsink   = getpar('mwsink')
          lthick   = getpar('lthick')
          rhosed   = getpar('rhosed')
          rhow     = getpar('rowass')
          dm0      = getpar('dm0')
          dm1      = getpar('dm1')
          visk     = getpar('vismol')
          lsink    = getpar('lsink')
          lhind    = getpar('lhind')
          sinkaccel= getpar('sinkaccel')
          llambda  = getpar('llambda')
          ladvlam  = getpar('ladvlam')
          lfloc    = getpar('lfloc')
          ladvfloc  = getpar('ladvfloc')
          lsetbound = getpar('lsetbound')

          dmf_mud  = dm1 
          z0bkmud  = 0.
          nf       = 2.

          do k = 1,nkn

!           --------------------------------------------------
!           add mud at first time step ... for Malcharek exp. 
!           --------------------------------------------------
            if (linitmud) then
              lmax = ilhkv(k)
              if (xgv(k) .gt. 7.) then
                do l = 1,lmax
                  if (l.lt.lthick) then
                    mudc(l,k) = 0.
                  elseif(l.ge.lthick) then
                    mudc(l,k) = mudref
                  end if
                end do
              elseif (xgv(k) .lt. 0.) then
                do l = 1,lmax
                  if (l.lt.lthick) then
                    mudc(l,k) = 0.
                  elseif(l.ge.lthick) then
                    mudc(l,k) = mudref
                  end if
                end do
              endif 
            endif
!           --------------------------------------------------
!           Sets sink vel. ... this is overwritten ... 
!           --------------------------------------------------
            call get_mudrho(k,rhow,dm0,ldebug,testnode,icycle)
!           --------------------------------------------------
!           Sets sink vel. ... this is overwritten ... 
!           --------------------------------------------------
            call get_cgel(k,dm0,cgel,phigel,ldebug,testnode,icycle)
!           --------------------------------------------------
!           Sets sink vel. ... this is overwritten ... 
!           --------------------------------------------------
            call get_mudc(k,rhow,dm0,phi,phip,ldebug,testnode,icycle)
!           --------------------------------------------------
!           Sets sink vel. ... this is overwritten ... 
!           --------------------------------------------------
            if (lsink) call set_hind_wsink(k,nlvdi,nkn,phi,phip,
     &                         phigel,ldebug,testnode,icycle,
     &                         sinkaccel,circle,lhind)
!           --------------------------------------------------
!           Initialize Lambda 
!           --------------------------------------------------
            lambda = 1.
          end do
!         --------------------------------------------------
!         Sets boundary conditions for all state variables
!         --------------------------------------------------
          fmbnd = 0.
          nintp = 2

          call bnds_init(what,mud2dn,nintp,nsmud,nb3dim,bnd3_mud,fmbnd)

          if (llambda) call bnds_init(what,lam2dn,nintp,nslam,
     &                                nb3dim,bnd3_lam,fmbnd)
          if (lfloc)  call bnds_init(what,dmf2dn,nintp,nsdmf,
     &                                nb3dim,bnd3_dmf,fmbnd)
!         -------------------------------------------------------------------
!         Store boundary 1 (in) nodes number into array
!         -------------------------------------------------------------------
          ibc = 1
          nk = 0

          if (circle) then
           call kanfend(ibc,kranf,krend)
           nnode = krend-kranf+1
           allocate(innode(nnode));innode=0
           do k = kranf,krend
             nk = nk + 1
             knode = irv(k)
             write(*,*) k, kranf, krend, knode, 'in nodes'
             innode(nk) = knode
           end do
!         -------------------------------------------------------------------
!         Store boundary 2 (out) nodes number into array
!         -------------------------------------------------------------------
           ibc = 2
           nk = 0
           call kanfend(ibc,kranf,krend)
           nnode = krend-kranf+1
           allocate(ounode(nnode));ounode=0
           do k = kranf,krend
             nk = nk + 1
             knode = irv(k)
             write(*,*) k, kranf, krend, knode, 'out nodes'
             ounode(nk) = knode
           end do
          endif
!         --------------------------------------------------
!         Initializes output
!         --------------------------------------------------
          id = 21
          ius = 0
          call confop(ius,itmcon,idtcon,nlv,1,'mud')
          write(6,*) 'fluid mud model initialized...'
          icall = 1
          return
        endif ! 1st call 

        if (ldumpmud) then
	  write(6,*)'ldump ha un errore'
          do k = 1,nkn
c            if(k.eq.349.or.k.eq.357.or.k.eq.369.or.k.eq.368.or.!DEB
c     +	k.eq.365.or.k.eq.360.or.k.eq.353)then!DEB
            if(k.eq.5004.or.k.eq.5005.or.k.eq.5003.or.k.eq.5002.or.!DEB
     +	k.eq.5000)then!DEB
            lmax = ilhkv(k)
            !if (xgv(k) .gt. 7. .or. xgv(k) .lt. 0.) then
            do l = 1,lmax
              if (l.lt.lthick) then
                mudc(l,k) = 0.
              elseif(l.ge.lthick) then
                !fak = 0.5d0+0.5d0*tanh((t_now-t_start)/smooth)
                !if (fak .gt. 0.9999999) cycle
                !mudc(l,k) = fak * mudref + (1.-fak)*0.
                !!if(k==452)write(*,*) fak, mudc(l,k)
                !!if(k==57)write(680,*) l,fak, mudc(l,k)
                mudc(l,k) = mudref !DEB
              end if
            end do
          endif !DEB
            !endif
          enddo
        endif

        if (.false.) then
          do k = 1,nkn
            lmax = ilhkv(k)
            if (xgv(k) .gt. 199000) then
              do l = 1,lmax
                if (l.lt.lthick) then
                  mudc(l,k) = 0.
                elseif(l.ge.lthick) then
                  fak = 0.5d0+0.5d0*tanh((t_now-t_start)/smooth)
                  if (fak .gt. 0.9999999) cycle
                  mudc(l,k) = fak * mudref + (1.-fak)*0.
                  !if(k==103)write(*,*) fak, mudc(l,k)
                end if
              end do
            endif
          enddo
        endif
! -------------------------------------------------------------------
! Normal call
! -------------------------------------------------------------------
        tsec = it
        ivar = 1
!       -------------------------------------------------------------
!       Boundary condition for fluid mud
!       -------------------------------------------------------------
        call scal_bnd(what,tsec,bnd3_mud)
        if (ladvlam)  call scal_bnd(what,tsec,bnd3_lam)
        if (ladvfloc) call scal_bnd(what,tsec,bnd3_dmf)
!       -------------------------------------------------------------
!       Transport and diffusion of fluid mud
!       -------------------------------------------------------------
        wsink = 1.
        ivar = 0
        fact = 1.
	!wsinkv = 0.0001 !DEB for test winterverp const acc
        call scal_adv_mud(what,ivar,fact
     &               ,mudc,bnd3_mud
     &               ,mudhpar,wsink,wsinkv
     &               ,difhv,difv,difmol,mudref,circle,lsetbound)
!       -------------------------------------------------------------
!       Boundary condition for fluid mud
!       -------------------------------------------------------------
        if (ladvlam) call scal_bnd(what,tsec,bnd3_lam)
!       -------------------------------------------------------------
!       Transport of lambda
!       -------------------------------------------------------------
       if (ladvlam)  call scal_adv_mud(what,ivar,fact
     &              ,lambda,bnd3_lam
     &              ,mudhpar,wsink,wprvs
     &              ,difhv,difv,difmol,mudref,circle,.false.)
!       -------------------------------------------------------------
!       Boundary condition for flocs
!       -------------------------------------------------------------
        if (ladvfloc) call scal_bnd(what,tsec,bnd3_dmf)
!       -------------------------------------------------------------
!       Transport and diffusion of floc diameter
!       -------------------------------------------------------------
       if (ladvfloc)  call scal_adv_mud(what,ivar,fact
     &              ,real(dmf_mud),bnd3_dmf
     &              ,mudhpar,wsink,wprvs
     &              ,difhv,difv,difmol,mudref,circle,.false.)
!
        do k = 1, nkn
!       ------------------------------------------------------
!       compute new floc diameter
!        ------------------------------------------------------
          call get_floc_diam(k,nlvdim,dt,dm0,ldebug,testnode,icycle) ! Hollandisches Betrugerpack, bei gerinen Scherraten (0.5Hz-1Hz) und hohen Konzentration geht das Flokkulationsmodell von Winterwerp (1998) gegen einen unendlichen Korndurchmesser!
!       ------------------------------------------------------
!       Get mud density  
!       -----------------------------------------------------
          call get_mudrho(k,rhow,dm0,ldebug,testnode,icycle)
!       ------------------------------------------------------
!       Get gelling concentration  
!       -----------------------------------------------------
          call get_cgel(k,dm0,cgel,phigel,ldebug,testnode,icycle)
!       ------------------------------------------------------
!       Get mud vol. concentration  
!       ------------------------------------------------------
          call get_mudc(k,rhow,dm0,phi,phip,ldebug,testnode,icycle)
!       -------------------------------------------------------
!       Update settling velocity
!       -------------------------------------------------------
          if (lsink) call set_hind_wsink(k,nlvdim,nkn,phi,phip,phigel,
     &                        ldebug,testnode,icycle,
     &                        sinkaccel,circle,lhind)
          !wsinkv(:,k) = 0.
!       -------------------------------------------------------
!       Update structural parameter lambda 
!       -------------------------------------------------------
          if (llambda) call get_lambda(k,nlvdim,dt,
     &                               ldebug,testnode,icycle)

        end do ! nodes
!       -------------------------------------------------------
!       Check mass conservation
!       -------------------------------------------------------
        call masscons(it,mudc)
!       -------------------------------------------------------
!       Fill output arrays
!       -------------------------------------------------------
        call confil(ius,itmcon,idtcon,id,nlvdi,mudc)
!        call confil(ius,itmcon,idtcon,id,nlvdi,lambda)
!        call confil(ius,itmcon,idtcon,id,nlvdi,vts)
!        call confil(ius,itmcon,idtcon,id,nlvdi,real(dmf_mud))

        !write(*,*) sum(wsinkv), 'from mud'

        end subroutine submud
! -------------------------------------------------------------
!
! -------------------------------------------------------------
       subroutine get_mudc(k,rhow,dm0,phi,phip,ldebug,testnode,icycle)

	use mod_fluidmud
	use levels

        implicit none

        include 'param.h'

        integer k,l,nlev
        double precision phip(nlvdim),phi(nlvdim)
        double precision rhow          !Sediment mineral density (kg/m3)
        double precision dm0           !Actual and primary diameter


	integer lmax
        integer testnode,icycle,iwrite
        logical ldebug
        save iwrite
        data iwrite/0/


        lmax=ilhkv(k)

        do l = 1, lmax
          phip(l)=mudc(l,k)/rhosed
          phi(l) =mudc(l,k)/rhosed*(dmf_mud(l,k)/dm0)**(3-nf(l,k))
          !phi(l)=((rhosed-rhow)/(rhomud(l,k)-rhow))*(mudc(l,k)/rhosed)
          if (ldebug .and. k == testnode .and. 
     &                  mod(iwrite,icycle) .eq. 0) then
            write(1111,'(A10,I10,10F16.8)') 'get_mudc',l, 
     &                                mudc(l,k), phip(l),phi(l),nf(l,k),
     &                                rhosed, rhow, 
     &                                rhomud(l,k),
     &                                mudc(l,k) / rhosed,
     &                       (dmf_mud(l,k)/dm0)**(3-nf(l,k))
          endif
        end do

        if (k == testnode) iwrite = iwrite + 1

       end subroutine get_mudc
! -------------------------------------------------------------
!
! -------------------------------------------------------------
       subroutine get_mudrho(k,rhow,dm0,ldebug,testnode,icycle)

	use mod_fluidmud
	use levels

        implicit none
! Purpose:
! Computes mud density based on fractal dimension and ratio between actual and primary diameter
!
        include 'param.h'

        integer k,l,nlev
        double precision rhow          !Sediment mineral density (kg/m3)
        double precision dm0           !Actual and primary diameter
	integer lmax
        integer testnode,icycle,iwrite
        logical ldebug
        save iwrite
        data iwrite/0/

        lmax=ilhkv(k)

        do l = 1, lmax
          rhomud(l,k) = rhow + (rhosed - rhow) * 
     &                  (dm0/dmf_mud(l,k))**(3.-nf(l,k))
            if (ldebug .and. k == testnode .and.
     &                  mod(iwrite,icycle) .eq. 0) then
             write(1112,'(A10,I10,7F15.8)') 'get_mudrho', l, 
     &                  rhow, rhosed, rhomud(l,k),
     &                  dm0, dmf_mud(l,k),nf(l,k)
           endif
        enddo

        if (k == testnode) iwrite = iwrite + 1

       end subroutine
! -------------------------------------------------------------
!
! -------------------------------------------------------------
       subroutine set_rhomud(k,l,rhop)

	use mod_fluidmud

        implicit none
! Purpose:
! Computes computes baroclinic forcing due to mud 
!
        include 'param.h'

        real*8 drho, maxdrho

        integer, intent(in) :: l,k
        real, intent(inout) :: rhop 

        maxdrho = 10.

!        write(*,*) rhop, mudc(l,k), (1.-(rhop+1000.)/rhomud(l,k)), 
!     &             mudc(l,k)*(1.-(rhop+1000.)/rhomud(l,k))

        drho = mudc(l,k)*(1.-(rhop+1000.)/rhosed) ! using the mud concentration is very tricky ...

        if (drho .gt. maxdrho) then
           drho = MIN(maxdrho,drho)
        endif

        rhop = rhop + drho ! using the mud concentration is very tricky ...

      end subroutine
! -------------------------------------------------------------
!
! -------------------------------------------------------------
       subroutine get_cgel(k,dm0,cgel,phigel,ldebug,testnode,icycle)

	use mod_fluidmud
	use levels

        implicit none
! Purpose:
! Computes the gelling concentration based on fractal dimension 
! The gelling concentration defines the state of mud where already a certain network is formed and the sink-velocity is strongly hindered

        include 'param.h'

        integer k,l,nlev


        double precision rhow          !Sediment mineral density (kg/m3)
        double precision dm0           !Actual and primary diameter
        double precision cgel(nlvdim)  ! gel, concentration
        double precision phigel(nlvdim)! gel, concentration
	integer lmax
        integer testnode,icycle, iwrite
        logical ldebug
        save iwrite
        data iwrite/0/

        lmax = ilhkv(k)

        do l = 1, lmax
          cgel(l) = rhosed * (dm0/dmf_mud(l,k))**(3.-nf(l,k))
          phigel(l) = 1. 
          if (ldebug .and. k == testnode .and.
     &                  mod(iwrite,icycle) .eq. 0) then
          write(1113,'(A10,I10,8F15.8)') 'get_cgel',l,
     &            rhomud(l,k),dmf_mud(l,k),dm0,nf(l,k),
     &           (dmf_mud(l,k)/dm0)**(3.-nf(l,k)),
     &            cgel(l),phigel(l)
           endif
        end do

        if (k == testnode) iwrite = iwrite + 1
        !write(*,*) 'entering get_cgel'
       end subroutine
! -------------------------------------------------------------
!
! -------------------------------------------------------------
        subroutine set_hind_wsink(k,nlvdi,nkn,phi,phip,phigel,
     &                            ldebug,testnode,icycle,sinkaccel,
     &                            lcircle,lhind)

	use mod_bound_geom
	use mod_sinking
	use mod_fluidmud
	use mod_ts
	use levels

        implicit none

! Purpose:
! Computes computes the settling vel. 
!

        include 'param.h'

        integer nlvdi,nkn               !number of level and nodes
        double precision d              !grain size
        double precision visk           !molecular viskosity 
        real getpar,vismol

        integer l,k,lmax


        double precision :: rhost,nu,dst,ws0,ws,rhow,fak,dmi
        double precision, intent(in):: phi(nlvdim)         !Fluid density (kg/m3)
        double precision, intent(in):: phip(nlvdim)         !Fluid density (kg/m3)
        double precision, intent(in):: phigel(nlvdim)          !Fluid density (kg/m3)
        integer testnode,icycle,iwrite
        integer kranf,krend
        integer kin,kou,nnode,ibc,knode,icall,nk,kk
        logical ldebug, lcircle, lhind
        real smooth,fak2,sinkaccel
        save iwrite,icall,nk,nnode

        data iwrite/0/
        data icall/0/

        rhow = getpar('rowass')
        vismol = getpar('vismol')
        smooth = 0.1

        lmax = ilhkv(k)
 
        wsinkv(0,k) = 0.

        do l=1,lmax
          rhost = (rhomud(l,k)-rhow)/rhow ! 1.65
          dmi   = dmf_mud(l,k)
          nu    = 1.73E-6 !vismol!max(vismol,visv(l,k)) ! 1.E-6
          dst   = (rhost*9.81/nu**2)**(1./3.)*dmi ! approx. 3.5 for particles of 200micros
          ws0   = 11*nu/dmi*(SQRT(1+.01*dst**3)-1) 
          fak   = (1.-min(1.d0,phi(l))*(1.-phip(l)))/(1+2.5*phi(l)) ! hindered settling 
          if (lhind) then
            ws = ws0 * fak * sinkaccel
          else 
            ws = ws0 * sinkaccel
          endif
          fak2   = 0.5d0+0.5d0*tanh((phi(l)-phigel(l))/smooth)
          wsinkv(l,k) = fak2 * 0. + (1.-fak2)*ws
          if(l.gt.1) wsinkv(l-1,k) = fak2*wsinkv(l,k)+ ! ????
     &                             (1.-fak2)*wsinkv(l-1,k)
          if((phi(l)-phigel(l)).gt.0.01 .and.l.gt.1)then ! no consolidation model yet ...
            wsinkv(l,k) = 0.d0
            wsinkv(l-1,k) = 0.d0
          endif
          if (l==lmax) then 
            wsinkv(l,k) = 0.
            ws          = 0.
          endif 
          if (wsinkv(l,k) .ne. wsinkv(l,k)) then
            write(*,*) dst, rhost, nu, dmi 
            stop 'nan in wsinkv'
          endif
        end do

        do l = 1, lmax
          if (ldebug .and. k == testnode .and.
     &                  mod(iwrite,icycle) .eq. 0) then
             !write(*,*) testnode, l, k, wsinkv(l,k),'from mud'
             write(1114,'(A10,I10,11F16.8)') 'set_wsink',
     &         l,
     &         dmf_mud(l,k),
     &         rhomud(l,k),
     &         phi(l),
     &         phip(l),
     &         phigel(l),
     &         -wsinkv(l,k)*1000. 
           endif
        end do

        if (k == testnode) iwrite = iwrite + 1

        end subroutine
! -------------------------------------------------------------
!
! -------------------------------------------------------------
        subroutine stress_mud(nldim,k,ldebug,testnode,icycle)

	use mod_fluidmud
	use mod_gotm_aux
	use mod_ts
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer nldim

        include 'param.h'

! Purpose:
! Computes the rheoligical stresses and the viscosity of the mud
!

	include 'pkonst.h'


        integer k,l,nlev
        real*8 aux,dh,du,dv,m2,dbuoy,tau_test,lambda_e,visk_new
        real*8 rho1, rho2, visk_bar,mu8,mu0,beta,c,tau0,visklim
        real h(nlvdim)
        real*8 cnpar,stress_x,stress_y, rhobar, g_dot, rhop, tau,visk
        real*8 g2_dot,g_dot_thr,smooth,viskmax,tf,tau_new,li, dvisk
        real rhop_r,tau0_r

        integer testnode,icycle,iwrite
        logical ldebug
        save iwrite
        data iwrite/0/

        visklim = 1.

        if( nldim .ne. nlvdim ) stop 'error stop stress_mud: dimension'
          call dep3dnod(k,+1,nlev,h)
          vts(0,k) = 0.
          do l=1,nlev
            rhobar = rhov(l,k)/1000.
            g_dot = sqrt(shearf2(l,k))
            rhop = min(0.2,max(0.d0,rhobar))
            call set_toorman_constants(rhop,mu8,mu0,beta)
            c  = mu0 - mu8
            rhop_r = rhop
            call set_yieldstress(rhop_r,tau0_r)
            tau0 = tau0_r
            g_dot_thr = 0.000001
            smooth    = 2.
            viskmax   = 100.
            visk      = 0.
            lambda_e  = 0.
            tau       = 0.
            tau_test  = 0.
!Toormann mit Oberrecht - equilibrium flow curve; toorman eq. 6
            lambda_e = 1.d0/(1.d0+beta*g_dot)
            tau     = tau0+(mu8+c+beta*tau0*lambda_e)*g_dot
!Toorman mit Oberrecht - full solution based on toorman eq. 11
            li      = lambda(l,k)
            tau_new = li*tau0+(mu8+li*c+beta*tau0*lambda_e)*g_dot
!dtau/dg_dot for eq. 6
            visk=-beta**2*tau0*g_dot/(1.d0+beta*g_dot)**2+
     &            mu8+c*+beta*tau0/(1.+beta*g_dot)
!dtau/dg_dot for eq. 11 is missing ... maple!
            visk_new=-beta**2*tau0*g_dot/(1.d0+beta*g_dot)**2+
     &            mu8+c*li*+beta*tau0/(1.+beta*g_dot)
            visk = visk/1000.
            dvisk = min(viskmax,visk)-vts(l,k)
            if (dvisk .gt. visklim) then
              dvisk = SIGN(MIN(visklim,ABS(dvisk)),dvisk)
            endif
            vts(l,k) = min(viskmax,vts(l,k)+dvisk) 
            if (ldebug .and. k == testnode .and.
     &                  mod(iwrite,icycle) .eq. 0) then            
                write(1115,'(A10,I5,17F14.8)') 'mudvisc',l,rhop,
     &                   mu0,mu8,beta,tau,tau_new,visk,visk_new,  
     &                   g_dot,tau/g_dot,tau_new/g_dot
            endif
         end do

         if (k == testnode) iwrite = iwrite + 1

      end subroutine
! -------------------------------------------------------------
!
! -------------------------------------------------------------
        subroutine set_yield(nldim,k,tstress,ldebug,testnode,icycle)

	use mod_fluidmud
	use mod_gotm_aux
	use mod_ts
	use mod_roughness
	use mod_diff_visc_fric
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer nldim

        include 'param.h'
! Purpose:
! Computes the yield stress and sets the certain visc. to model the yield stress regime
!
	include 'pkonst.h'

        real tstress(nlvdim)

        integer k,l,nlev
        real aux,dh,du,dv,m2,dbuoy,tau_test,lambda_e
        real rho1, rho2, visk_bar,mu8,mu0,beta,c,tau0
        real h(nlvdim)
        real cnpar,stress_x,stress_y, rhobar, g_dot, rhop, tau,visk
        real g2_dot,g_dot_thr,smooth,viskmax,tf
        real fak

        integer testnode,icycle,iwrite
        logical ldebug
        save iwrite
        data iwrite/0/

        smooth = 0.001

        if( nldim .ne. nlvdim ) stop 'error stop stress_mud: dimension'

          call dep3dnod(k,+1,nlev,h)

          do l=1,nlev
            rhop = max(0.d0,dble(rhov(l,k)))/1000.
            call set_yieldstress(rhop,tau0)
            if (tau0 .gt. 0. .and. rhop .gt. 0.) then
              fak = 0.5d0+0.5d0*tanh((tstress(l)-tau0)/smooth)
              visv_yield(l,k) = fak * visv(l,k) + (1.-fak) * 99.
              difv_yield(l,k) = fak * difv(l,k) + (1.-fak) *  0.
              if (l==nlev) then 
                z0bkmud(k)  = fak * z0bk(k)   + (1.-fak) * 99.
              else
                z0bkmud(k)  = 0.
              endif
            else
              fak = 999.
              visv_yield(l,k) = visv(l,k)
              difv_yield(l,k) = difv(l,k)
              if (l==nlev) then
                z0bkmud(k) = z0bk(k)
              else
                z0bkmud(k) = 0.
              endif
            endif
            if (ldebug .and. k == testnode .and.
     &                  mod(iwrite,icycle) .eq. 0) then
              write(1118,'(A10,I5,17F14.8)') 'set yield',
     &        l, rhop,tstress(l)-tau0,tstress(l),tau0,
     &        visv_yield(l,k),difv_yield(l,k),
     &        visv(l,k),difv(l,k),z0bk(k),z0bkmud(k),fak, 
     &        tanh((tstress(l)-tau0)/smooth)
            endif
          end do
          
          if (k == testnode) iwrite = iwrite + 1

      end subroutine
! -------------------------------------------------------------
!
! -------------------------------------------------------------
      SUBROUTINE set_yieldstress(rhop,tau_y)

        IMPLICIT NONE
!
        REAL, INTENT(IN)  :: rhop
        REAL, INTENT(OUT) :: tau_y

        tau_y = 5000.d0*rhop**3-340.d0*rhop**2+10.d0*rhop

      END SUBROUTINE
! -------------------------------------------------------------
!
! -------------------------------------------------------------
      SUBROUTINE set_toorman_constants(rhop,mu8,mu0,beta)

        IMPLICIT NONE
!
        REAL*8, INTENT(IN)  :: rhop
        REAL*8, INTENT(OUT) :: mu8,mu0,beta

        REAL*8  :: rhop2

        rhop2= rhop*rhop

        mu8  = 1.16*rhop2+0.1*rhop
        mu0  = 0.012*exp(63*rhop)-0.012
        beta = 72*rhop2+6*rhop

      END SUBROUTINE
! -------------------------------------------------------------
!
! -------------------------------------------------------------
        subroutine set_mud_roughness(k,lmax,alpha)
!
!Purpose:
!Computes the roughness reduction due to the presence of mud

	use mod_sinking
	use mod_fluidmud
	use mod_gotm_aux
	use mod_ts
	use mod_diff_visc_fric
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer nldim

        include 'param.h'

	integer lmax
	real alpha
	include 'pkonst.h'


        integer k,nlev

        real*8 aux,dh,du,dv,m2,dbuoy,dh1,dh2,rho3
        real*8 rho1, rho2, visk_bar,tstress,ufric
        real*8 cnpar,stress_x,stress_y, rhobar
        real*8  g_dot, rhop, tau,visk, drho, ri

        real, parameter :: bpar = 1.
        real, parameter :: beta = 0.7
        real, parameter :: mpar = 1.


	integer i
        real h(nlvdim)

        call dep3dnod(k,+1,nlev,h)

        lmax = ilhkv(k)
        rho1 = rhov(lmax  ,k)+1000.
        rho2 = rhov(lmax-1,k)+1000.
        rho3 = rhov(lmax-2,k)+1000.
        dh1 = 0.5 * ( h(lmax) + h(lmax-1) )
        dh2 = 0.5 * ( h(lmax-1) + h(lmax-2) )
        drho = 0.5*(abs(rho2-rho1)/dh1 + abs(rho3-rho2)/dh2)
        rhobar = 1./3. * (rho1+rho2+rho3)
        g_dot = sqrt(shearf2(lmax,k)) 
        if (g_dot .gt. 0.000001) then
!AR: take care with the mud!
          visk_bar = visv(lmax,k) + vts(lmax,k)
          tstress  = visk_bar * g_dot
          ufric = sqrt(tstress)
        else
          ufric = 0.
        endif 
        if (ufric .gt. 0.00001) then
          ri = 9.81/rhobar*drho/g_dot**2
          alpha = exp(-(1+beta*wsinkv(lmax-1,k)/ufric)*
     &         (1-exp(-bpar*ri**mpar)))
          !write(*,*) k, alpha, drho
        else 
          alpha = 1.
        endif

        alpha = 1.

      end subroutine
! -------------------------------------------------------------
!
! -------------------------------------------------------------
        subroutine set_bottom_stress(k,tstress)

	use mod_sinking
	use mod_fluidmud
	use mod_gotm_aux
	use mod_ts
	use mod_diff_visc_fric
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none
!Purpose: Computes discrete bottom stress 

        integer nldim

        include 'param.h'

	integer lmax
	real alpha
	include 'pkonst.h'


        integer k,l,nlev
        real tstress
        real*8 aux,dh,du,dv,m2,dbuoy
        real*8 rho1,rho2,visk_bar,ufric
        real*8 cnpar,stress_x,stress_y, rhobar
        real*8  g_dot, rhop, tau,visk, drho, ri


        lmax = ilhkv(k)
        g_dot = sqrt(shearf2(lmax,k))
        visk_bar = visv(lmax,k) !+ vts(lmax,k)
        tstress  = visk_bar * g_dot 
        !write(*,*) k, visk_bar, tstress, du, dv, dh,lmax
      end subroutine
! -------------------------------------------------------------
!
! -------------------------------------------------------------
        subroutine set_bottom_visk(k,ufric,depth,visk1,visk2)

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer nldim

!Purpose: computes the turb. visc. based on mixing length. Only for the 1st two layers.

        include 'param.h'

	integer lmax

        integer k,l,nlev
        real*8 ufric,depth
        real visk1,visk2
        real h(nlvdim)

        call dep3dnod(k,+1,nlev,h)

        lmax = ilhkv(k)

!AR: 2do check dep3nod with respect to nlev
        visk1 = 0.41 * ufric * h(lmax-1)*(1 - h(lmax-1)/depth)
        visk2 = 0.41 * ufric * h(lmax-2)*(1 - h(lmax-2)/depth)
!        write(*,'(I10,5F15.10)')k,visk1,visk2,h(nlev-1),depth,ufric
      end subroutine
! -------------------------------------------------------------
!
! -------------------------------------------------------------
      subroutine get_floc_diam(k,nldim,dt,dm0,ldebug,testnode,icycle)
	use mod_fluidmud
	use mod_gotm_aux
	use mod_ts
	use mod_diff_visc_fric
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer :: k,l,nldim,nlev

!Purpose: Compute Flocculation

        include 'param.h'

	include 'pkonst.h'



        real, intent(in)    :: dt                   ! time step
        double precision, intent(in)    :: dm0                  ! Size of primary floc

        double precision, parameter     :: pi = 3.14159265

        double precision, parameter     :: delta = 3.357          ! Coefficient
        double precision, parameter     :: psi = -0.093
        double precision, parameter     :: p  = 1.0
        double precision, parameter     :: q  = 0.5 


        double precision, parameter     :: Fc  = 2.0      ! Characteristic fractal dimension ! formel4 seite 59 khelifa and hill 
        double precision, parameter     :: dfc = 2.E-6    ! Characteristic size of floc ! 2microMeter laut formel4 seite 59 khelifa and hill
        double precision, parameter     :: Fy  = 1.E-10   ! 10^-10 Yield strength of floc ( N )

        double precision, parameter     :: ka = 0.98      ! emp. coeff.
        double precision, parameter     :: kb = 3.3E-5    ! emp. coeff.

        double precision, parameter     :: ka2 = 14.6!0.98            ! emp. coeff.
        double precision, parameter     :: kb2 = 14.E3!3.3E-5      ! emp. coeff.

        double precision, parameter     :: dequi = 300.E-6

        double precision, parameter     :: one_third = 1.d0/3.d0
        double precision, parameter     :: two_third = 2.d0/3.d0

        double precision :: dddt, dddt1, dddt2, dddt3, conz
        double precision :: alpha, beta, dv, du, dh, cnpar, g_dot_thr
        double precision :: g_dot,g2_dot,dold,dnew,mu,coeff,dinit
        double precision :: growth, decay, rk(3), diffd

        integer testnode,iwrite,icycle
        logical ldebug

        real h(nlvdim)
        save iwrite
        data iwrite/0/

        g_dot_thr = 0.0000001
!AR: tak care with visv should be mu
        nlev = ilhkv(k)

        mu = 1.8d-3!*(rhov(l,k)+1000.)
        beta = log(Fc/3.d0)/log(dfc/dm0)

        do l=1,nlev

           g_dot  = sqrt(shearf2(l,k))

           dold = dmf_mud(l,k)
           dinit = dold 

           conz = max(0.,mudc(l,k))/rhomud(l,k)
           !conz = 0.65/2650.
           !g_dot = 7.3

!Fractal dimension is a serious problem ... here Miki, JGR 2007
           nf(l,k)= delta*(dinit/dm0)**psi
           !nf(l,k)= 2.

           coeff=(g_dot*dm0**beta)/(beta*log(dold/dm0)+1.d0)
           growth=one_third*ka*conz*
     &            dm0**(nf(l,k)-3)*dold**(-nf(l,k)+4.-beta)
           decay=one_third*kb*(mu*g_dot/fy)**0.5*dm0**(-1)*
     &        dold**(-beta+2)*(dold-dm0)
           dddt1 = coeff*(growth-decay)
           dnew  = dold + dddt1 * dt
           dold  = max(dm0,min(dequi,dnew))

           if (dnew .ne. dnew ) then
             write(*,*) l, nf(l,k),dold/dm0,
     &                     g_dot,nf(l,k),
     &                     beta, dm0, visv(l,k),
     &                     dt,growth,decay,dold*1E6,dnew*1E6
             stop 'NaN in dmf_mud rk1'
           else if (dnew .lt. 0.) then
             write(*,*) l, nf(l,k),dold/dm0,
     &                     g_dot,nf(l,k),
     &                     beta, dm0, visv(l,k),
     &                     dt,growth,decay,dold*1E6,dnew*1E6
             stop 'neg. diam in dmf_mud rk1'
           endif

!           goto 999

!Fractal dimension is a serious problem ... here Miki, JGR 2007
           nf(l,k)= delta*(dinit/dm0)**psi
           !nf(l,k)= 2.

           coeff=(g_dot*dm0**beta)/(beta*log(dold/dm0)+1.d0)
           growth=one_third*ka*conz*
     &            dm0**(nf(l,k)-3)*dold**(-nf(l,k)+4.-beta)
           decay=one_third*kb*(mu*g_dot/fy)**0.5*dm0**(-1)*
     &        dold**(-beta+2)*(dold-dm0)
           dddt2 = coeff*(growth-decay)
           dnew = 0.75d0*dmf_mud(l,k)+0.25d0*dnew+0.25d0*dddt2*dt
           dold  = max(dm0,min(dequi,dnew))

!Fractal dimension is a serious problem ... here Miki, JGR 2007
           nf(l,k)= delta*(dnew/dm0)**psi
           !nf(l,k)= 2.

           if (dnew .ne. dnew ) then
             write(*,*) l, nf(l,k),dold/dm0,
     &                     g_dot,nf(l,k),
     &                     beta, dm0, visv(l,k),
     &                     dt,growth,decay,dold,dnew
             stop 'NaN in dmf_mud rk2'
           else if (dnew .lt. 0.) then
             write(*,*) l, nf(l,k),dold/dm0,
     &                     g_dot,nf(l,k),
     &                     beta, dm0, visv(l,k),
     &                     dt,growth,decay,dold,dnew

             stop 'neg. diam in dmf_mud rk2'
           endif

           coeff=(g_dot*dm0**beta)/(beta*log(dold/dm0)+1.d0)
           growth=one_third*ka*conz*
     &            dm0**(nf(l,k)-3)*dold**(-nf(l,k)+4.-beta)
           decay=one_third*kb*(mu*g_dot/fy)**0.5*dm0**(-1)*
     &        dold**(-beta+2)*(dold-dm0)
!           dddt3 = sign(min(0.1*dequi,abs(coeff*(growth-decay))),
!     &             coeff*(growth-decay))
           dddt3 = coeff*(growth-decay)
           dnew = one_third*dmf_mud(l,k)+
     &            two_third*dnew+two_third*dddt3*dt
           dold  = max(dm0,min(dequi,dnew))

!Fractal dimension is a serious problem ... here Miki, JGR 2007
           nf(l,k)= delta*(dnew/dm0)**psi
           !nf(l,k)= 2.

           if (dnew .ne. dnew ) then
             write(*,*) l, nf(l,k),dold/dm0,
     &                     g_dot,nf(l,k),
     &                     beta, dm0, visv(l,k),
     &                     dt,growth,decay,dold,dnew
             stop 'NaN in dmf_mud rk3'
           else if (dnew .lt. 0.) then
             stop 'neg. diam in dmf_mud rk3'
           endif

999        continue

            if (ldebug .and. k == testnode .and.
     &                  mod(iwrite,icycle) .eq. 0) then
             write(1116,'(A10,I10,10F20.10)') 'get_floc', l, 
     &        mudc(l,k)/rhosed, nf(l,k), g_dot, dinit*1E6, dnew*1E6,
     &        dddt1*dt, coeff, growth*1E6, decay*1E6
           endif                                                                

        end do

        nf(nlev,k) = nf(nlev-1,k)
        dmf_mud(nlev,k) = dmf_mud(nlev-1,k)

        if (k == testnode) iwrite = iwrite + 1

        end subroutine get_floc_diam
! -------------------------------------------------------------
!
! -------------------------------------------------------------       
        subroutine get_lambda(k,nldim,dt,ldebug,testnode,icycle)
	use mod_fluidmud
	use mod_gotm_aux
	use mod_ts
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer :: k,l,nldim,nlev
!
!Purpose: Compute memory effects ... 
!
        include 'param.h'

	include 'pkonst.h'

        logical ldebug
        integer testnode, icycle



        real, intent(in)    :: dt                   ! time step

        double precision, parameter     :: pi = 3.14159265
        double precision, parameter     :: gdot_thr = 10.E-5 

        double precision :: dddt, dddt1, dddt2, dddt3
        double precision :: lambda_new, lambda_old, lambda_e
        double precision :: growth, decay, rk(3), rhobar, g_dot
        double precision :: rhop,mu8,mu0,beta,b_lambda

        double precision, parameter :: a_lambda = 0.03 ! for the ems this has a range from 0.03-0.07
        double precision, parameter :: one_third = 1.d0/3.d0
        double precision, parameter :: two_third = 2.d0/3.d0

        integer iwrite
        real h(nlvdim),li
        save iwrite
        data iwrite/0/

        nlev = ilhkv(k)

        do l=1,nlev
           rhobar = rhov(l,k)/1000.
           g_dot = sqrt(shearf2(l,k))
           rhop = max(0.d0,rhobar)
           call set_toorman_constants(rhop,mu8,mu0,beta)
           b_lambda = beta * a_lambda
           li       = lambda(l,k)

           lambda_e = 1.d0/(1.d0+beta*g_dot) ! eq. struct. parameter

             dddt1 = -(a_lambda+b_lambda*g_dot)*(li-lambda_e) ! source term
           lambda_new = li + dddt1 * dt ! rk3 ... 1step 
           lambda_old = lambda_new ! copy

             dddt2 = -(a_lambda+b_lambda*g_dot)*(lambda_new-lambda_e) ! source term
           lambda_new = 0.75d0*li+0.25d0*lambda_new+0.25d0*dddt2*dt ! rk3 ... 2step 
           lambda_old = lambda_new ! copy

             dddt3 = -(a_lambda+b_lambda*g_dot)*(lambda_new-lambda_e) !rk3 ... 3step 
           lambda_new = one_third*li+
     &            two_third*lambda_new+two_third*dddt3*dt
           !if (lambda(l,k) .gt. 1.) lambda(l,k) = 1. 
           lambda(l,k) = min(1.,max(0.,lambda_new))
           if (ldebug .and. k == testnode .and.
     &                 mod(iwrite,icycle) .eq. 0) then
             write(1117,'(A10,I10,10F20.10)') 'get_lambda', l, 
     &       dddt1,dddt2,dddt3,g_dot,rhop,lambda(l,k)
           endif                                                                
        end do

        if (k == testnode) iwrite = iwrite + 1

        end subroutine get_lambda
! -------------------------------------------------------------
!
! -------------------------------------------------------------  
      SUBROUTINE RISK(Ustar,R)
      IMPLICIT NONE

      real(kind=8), intent(in)  :: ustar
      real(kind=8), intent(out) :: R

      IF ( Ustar >1.5 ) THEN
         R = 1.
      ELSEIF ( Ustar <=0.7 ) THEN
         R = 0.
      ELSE
!         R = 1./(10.*Ustar**(-18)+1.)
         R = SQRT((Ustar-0.7)/0.8)
      ENDIF

      END SUBROUTINE RISK
! -------------------------------------------------------------
!
! -------------------------------------------------------------
      SUBROUTINE SHIELDS(d,rhost,rhow,nu,usc)
      IMPLICIT NONE

      
      real(kind=8), intent(in)  :: d, nu, rhow 
      real(kind=8), intent(out) :: rhost, usc
      real(kind=8) :: rgn , frstc, rhosusp, dst

      rhost=(rhosusp-rhow)/rhow ! 
      rgn = (Rhost*9.81/nu**2)**(1./3.)
      dst = rgn*D
      IF ( dst<6. ) THEN
        frstc = .109*dst**(-.5)
      ELSEIF ( dst<10. ) THEN
        frstc = .14*dst**(-.64)
      ELSEIF ( dst<20. ) THEN
        frstc = .04*dst**(-.1)
      ELSEIF ( dst<150. ) THEN
        frstc = .013*dst**(.29)
      ELSEIF ( dst>=150. ) THEN
        frstc = .055
      ENDIF
      usc = SQRT(frstc*(Rhost*9.81*D))

      END SUBROUTINE SHIELDS
! -------------------------------------------------------------
!
! -------------------------------------------------------------
! SUBROUTINE READMUD
! This subroutine reads the simulation fluid mud parameter from the 
! input parameter .str file

        subroutine readmud

        implicit none

        character*80 name               !name of item read
        character*80 text               !text of item if string
        real value                      !value of item if numeric
        double precision dvalue         !value of item if numeric
        integer iweich
        integer nrdpar

! DOCS  START   S_mudpar_h
!
! DOCS  COMPULS        Compulsory fluid mud parameters
!
! These parameters are compulsory parameters that define if the
! fluid mud transport module is run and what kind of output is written.
!
! |mudsec|       Fluid mud section name.

        call sctpar('mudsec')             !sets default section
        call sctfnm('mudsec')

! |imud|        Flag if the computation on the fluid mud is done:
!               \begin{description}
!               \item[0] Do nothing (default)
!               \item[1] Compute fluid mud
!               \end{description}

        call addpar('imud',0.)

! |dm0|         Primary particle diameter [m]
!               (Default 0).

        call addpar('dm0',4e-06)

! |dm1|         Initial particle diameter [m]
!               (Default 0).

        call addpar('dm1',4e-06)

! |rhosed|      Primary particle density [kg/m3]
!               (Default 0).

        call addpar('rhosed',2650.)

! |mudref|      Fluid mud initial reference concentration [kg/m3]
!               (Default 0).

        call addpar('mudref',0.)

! |mudhpar|     Fluid mud diffusion coefficient (Default 0).

        call addpar('mudhpar',0.)

! |mwsink|      Fluid mud floc sinking velocity [m/s] (Default 0).

        call addpar('mwsink',0.)

! |lthick|      Initial fluid mud layer thickness [layer number] 
!               (Default 1 = whole water column).

        call addpar('lthick',1.)

! |ldebug|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('ldebug',0.)

! |tnode|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('testnode',0.)

! |linitmud|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('linitmud',0.)

! |ldumpmud|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('ldumpmud',0.)

! |lcircle|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('lcircle',0.)

! |t_start|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('t_start',0.)

! |smooth|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('smooth',0.)

! |icycle|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('icycle',1.)

! |lmudvisc|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('lmudvisc',0.)

! |lbaro|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('lbaro',0.)

! |lbaro|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('lsink',0.)

! |lbaro|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('lhind',0.)

! |lbaro|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('llambda',0.)

! |lbaro|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('ladvlam',0.)

! |lsinkaccel|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('sinkaccel',1.)

! |lsinkaccel|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('ladvfloc',0.)

! |lsinkaccel|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('lfloc',0.)
! |lsinkaccel|      Circular boundary conditions for Transport  
!               (Default 0 = No circular boundary conditions).

        call addpar('lsetbound',0.)

! DOCS  FILENAME        Boundary conditions
!
! Boundary conditions have to be given in a file in every
! section |bound|.
!
! |mud2dn|      File name that contains boundary conditions for
!               concentration of the fluid mud.

! DOCS  END

!       --------------------------------------------------
!       Starts the reading loop
!       --------------------------------------------------

        iweich = 1
        do while(iweich.ne.0)
          iweich=nrdpar('mudsec',name,dvalue,text)
          value = real(dvalue)
        end do

        return

        end
! -------------------------------------------------------------
!
! -------------------------------------------------------------
        subroutine scal_adv_mud(what,ivar,fact
     +                          ,scal,bnd3
     +                          ,rkpar,wsink,wsinkv
     +                          ,difhv,difv,difmol
     +                          ,mudref,lcircle,lsetbound)

!--------------------------------------------------------------
! shell for scalar (for parallel version)
! special version with factor for BC and variable sinking velocity
!--------------------------------------------------------------
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

        implicit none 
        include 'param.h'

	integer k,l !DEB
	character*(*) what
        integer ivar
        real fact           !factor for boundary condition
        real scal(nlvdim,nkndim)
        real bnd3(nb3dim,0:nbcdim)
        real mudref

        real rkpar
        real wsink
        real wsinkv(0:nlvdim,nkndim)
        real difhv(nlvdim,1)
        real difv(0:nlvdim,1)
        real difmol

        logical lcircle,lsetbound

        real bnd3_aux(nb3dim)
        real r3v(nlvdim,nkndim)
        real robs
        integer iwhat,ichanm
        character*20 whatvar,whataux

        robs = 0.
!--------------------------------------------------------------
! make identifier for variable
!--------------------------------------------------------------
        whatvar = what
        if( ivar .ne. 0 ) then
          write(whataux,'(i2)') ivar
          whatvar = what // whataux
        end if
        iwhat = ichanm(whatvar)
!--------------------------------------------------------------
! transfer boundary conditions of var ivar to 3d matrix r3v
!--------------------------------------------------------------
	call bnds_trans(whatvar(1:iwhat)
     +              ,nb3dim,bnd3,bnd3_aux
     +                          ,ivar,nlvdim,r3v)
        if(lcircle)then
          call circflume(r3v,scal) 
        else if (lsetbound) then
          call mudboundary(r3v,mudref)
        else if (lsetbound .and. lcircle) then
          stop 'invalid combination for lcircle, lsetbound'
	endif
	
	
	!write(679,*)'bordo',(r3v(l,154),l=1,nlv)!DEB
!--------------------------------------------------------------
! multiply boundary condition with factor
!--------------------------------------------------------------
        if( fact .ne. 1. ) then
          call mult_scal_bc(r3v,fact)
        end if
!--------------------------------------------------------------
! do advection and diffusion
!--------------------------------------------------------------
        call scal3sh(whatvar(1:iwhat)
     +              ,scal,nlvdim
     +                          ,r3v,scal,robs
     +              ,rkpar,wsink,wsinkv
     +                          ,difhv,difv,difmol)

       end subroutine
! -------------------------------------------------------------
!
! -------------------------------------------------------------
! This routine sets the boundary condition for the fluid mud 
! concentration for the case of anular flume in the way that conc
! going out from the system is imposed in the influx nodes

        subroutine circflume(r3v,scal)

	use mod_bound_geom
	use mod_bound_dynamic
	use mod_hydro_print
	use levels

        implicit none

        include 'param.h'

        real r3v(nlvdim,1)              !boundary conc array

        real scal(nlvdim,1)             !Fluid mud concentration array (kg/m3)      !ccf

        integer, allocatable :: innode(:)          !in boundary node number array (bound1)
        integer, allocatable :: ounode(:)          !out boundary node number array (bound2)

        integer ibc,knode
        integer kranf,krend
        integer kin,kou,nnode
        integer k,n,l,nk
        integer icall                   !Initialization parameter
        save icall
        save innode,ounode,nk,nnode
        data icall /0/


        integer lmax
        real diff

! ----------------------------------------------------------
! Initialization
! ----------------------------------------------------------

        if( icall .eq. 0 ) then
!         -------------------------------------------------------------------
!         Store boundary 1 (in) nodes number into array
!         -------------------------------------------------------------------
          ibc = 1
          nk = 0
          call kanfend(ibc,kranf,krend)
          nnode = krend-kranf+1
          allocate(innode(nnode));innode=0
          do k = kranf,krend
             nk = nk + 1
             knode = irv(k) 
             write(6,*) k, kranf, krend, knode, 'in nodes'
             innode(nk) = knode
          end do
!         -------------------------------------------------------------------
!         Store boundary 2 (out) nodes number into array
!         -------------------------------------------------------------------
          ibc = 2
          nk = 0
          call kanfend(ibc,kranf,krend)
          nnode = krend-kranf+1
          allocate(ounode(nnode));ounode=0
          do k = kranf,krend
             nk = nk + 1
             knode = irv(k)
             write(6,*) k, kranf, krend, knode, 'out nodes'
             ounode(nk) = knode
          end do
          icall = 1
        end if
!       -------------------------------------------------------------------
!       Normal call:
!       Set boundary conc for bound 1 equal to conc in bound 2 (invert order)
!       -------------------------------------------------------------------

        do k = 1,nk
          kin = innode(k)
          kou = ounode(nk-k+1)
          !write(*,*) k, kin, kou
          lmax = ilhkv(kin)
          do l = 1,lmax
            diff = scal(l,kin)-scal(l,kou)
!AR: here is the bug with the sink vel. it would be really great to fix this 
!... diff is not .eq. 0 when wsinkv is .ne. 0
!... so it must have something to do with the vertical flux 
!... since this worked before it has something to do with the corrections made by georg. 
            r3v(l,kou) = scal(l,kou)  
            r3v(l,kin) = scal(l,kou)  
            !r3v(l,kin) = scal(l,kou)*abs(mfluxv(l,kou)/mfluxv(l,kin)) !corrected
            !write(*,*) l, diff, scal(l,kin), scal(l,kou)
          end do
        end do
      end
! -------------------------------------------------------------
!
! -------------------------------------------------------------
        subroutine mudboundary(r3v,mudref)

	use mod_bound_geom
	use mod_fluidmud
	use levels, only : nlvdi,nlv

        implicit none

        include 'param.h'

        real,intent(in) :: mudref       !mudref  
 
        real r3v(nlvdim,1)              !boundary conc array


        integer, allocatable :: innode(:)          !in boundary node number array (bound1)
        integer, allocatable :: ounode(:)          !out boundary node number array (bound2)

        integer ibc,knode
        integer kranf,krend
        integer kin,kou,nnode
        integer k,n,l,nk
        integer icall, lthick                   !Initialization parameter
        save icall, lthick
        save innode,ounode,nk,nnode
        real getpar
        data icall /0/

! ----------------------------------------------------------
! Initialization
! ----------------------------------------------------------
! This section is called only the first time step when ICALL = 0

        if( icall .eq. 0 ) then
!         -------------------------------------------------------------------
!         Store boundary 1 (in) nodes number into array
!         -------------------------------------------------------------------
          ibc = 1
          nk = 0
          call kanfend(ibc,kranf,krend)
          nnode = krend-kranf+1
          allocate(innode(nnode),ounode(nnode));innode=0;ounode=0
          do k = kranf,krend
             nk = nk + 1
             knode = irv(k)
             innode(nk) = knode
          end do
!         -------------------------------------------------------------------
!         Store boundary 2 (out) nodes number into array
!         -------------------------------------------------------------------
          ibc = 2
          nk = 0
          call kanfend(ibc,kranf,krend)
          do k = kranf,krend
             nk = nk + 1
             knode = irv(k)
             ounode(nk) = knode
          end do
          icall = 1
          lthick   = getpar('lthick')
        end if
!       -------------------------------------------------------------------
!       Normal call:
!       Set boundary conc for bound 1 equal to conc in bound 2 (invert order)
!       -------------------------------------------------------------------

        do k = 1,nk
          kin = innode(k)
          kou = ounode(nk-k+1)
          do l = 1,nlvdi
            if (l.lt.lthick) then
              !r3v(l,kin)  = mudref 
              r3v(l,kou)  = 0.
            else
              !r3v(l,kin)  = mudref 
              r3v(l,kou)  = mudref
            endif
          end do
        end do

      end subroutine
! -------------------------------------------------------------
!
! -------------------------------------------------------------
      subroutine write_xfn
	use mod_fluidmud
	use mod_depth
	use mod_gotm_aux
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use basin

        implicit none

        integer :: k,l,nldim

        include 'param.h'

	include 'pkonst.h'



        logical, save :: lfirst
        data lfirst/.true./
 
        real h(nlvdim) , r
        real hh(nkndim),dep(nkndim)
        integer ie,ip,il,iwrite,icycle
        real,save :: time
        save iwrite
        data time/0./
        data iwrite/0/

        if (lfirst) then
          write(5001,*) 0
          write(5001,*) nkn
          do ip = 1, nkn
            call dep3dnod(ip,+1,nlvdi,h)
            dep(ip)= hkv(ip) 
            hh(ip) = h(nlvdi)
            write(5001,*) ip-1,xgv(ip),ygv(ip),dep(ip)
          end do
          write(5001,*) nel
          do ie = 1, nel
            write(5001,*)nen3v(1,ie)-1,nen3v(2,ie)-1,nen3v(3,ie)-1,0,ie
          enddo
        endif

        if (lfirst) then 
          OPEN(5002, FILE  = 'ergzus.bin'  , FORM = 'UNFORMATTED')
          OPEN(5003, FILE  = 'ergrho.bin'  , FORM = 'UNFORMATTED')
          OPEN(5004, FILE  = 'ergmudc.bin'  , FORM = 'UNFORMATTED')
          OPEN(5005, FILE  = 'erglam.bin'  , FORM = 'UNFORMATTED')
          OPEN(5006, FILE  = 'ergdmf.bin'  , FORM = 'UNFORMATTED')
          lfirst = .false.
        endif
 
        r = 1000.

        !write(*,*) mod(iwrite,icycle), iwrite, icycle
	if(lfirst)then !DEB
        do il = 1, nlvdi
            time = time + 1.
            WRITE(5002) TIME
            WRITE(5002) (uprv(il,ip), vprv(il,ip), hh(ip),ip = 1, nkn)
            CALL FLUSH(5002)
            WRITE(5003) TIME
            WRITE(5003) (uprv(il,ip),vprv(il,ip),rhov(il,ip)+r,ip=1,nkn)
            CALL FLUSH(5003)
            WRITE(5004) TIME
            WRITE(5004) (uprv(il,ip),vprv(il,ip),mudc(il,ip),ip=1,nkn)
            CALL FLUSH(5005)
            WRITE(5005) TIME
            WRITE(5005) (uprv(il,ip),vprv(il,ip),lambda(il,ip),ip=1,nkn)
            CALL FLUSH(5005)
            WRITE(5006) TIME
            WRITE(5006)(uprv(il,ip),vprv(il,ip),
     &                  real(1.E6*dmf_mud(il,ip)),ip=1,nkn)
            CALL FLUSH(5006)
        enddo
  	endif !DEB
      end 
! -------------------------------------------------------------
!
! -------------------------------------------------------------
        subroutine masscons(it,mudc)

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'


        integer it
        real mudc(nlvdim,nkndim)    !Fluid mud concentration array (kg/m3)
        real totm
        real volnode,vol

        integer k,l,lmax
        logical is_r_nan

        totm = 0.
        do k = 1,nkn
           lmax = ilhkv(k)
           do l = 1,lmax
            vol = volnode(l,k,1)
             totm = totm + mudc(l,k)*vol
             if( is_r_nan(mudc(l,k)) ) then
               write(6,*) '*** nan in mudc...'
               write(6,*) l,k,mudc(l,k)
              end if
          end do
        end do

        write(74,*)it,totm

        end
! -------------------------------------------------------------
!
! -------------------------------------------------------------
