!****************************************************************
!
! $Id: submud.f,v 1.0 2009-11-29 19:38:34 aron Exp $
!
!****************************************************************************
! ***********************************
! ------ SUBROUTINE SUBMUD ---------
! ***********************************
!
! This routine manages the fluid mud computation.

! Revision Log :
!
! Jun, 2010	ccf	routines written from scratch

!****************************************************************************

        subroutine submud(it,dt)

          implicit none

          integer it                        !time in seconds
          real dt                           !time step in seconds
  
        include 'param.h'
  
! -------------------------------------------------------------
! fem variables
! -------------------------------------------------------------

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer ilhkv(nkndim)           !number of node levels
        common /ilhkv/ilhkv
        real difv(0:nlvdim,nkndim)	!vertical diffusivity
        real difvm(0:nlvdim,nkndim)  !vertical diffusivity for transport of mud
        common /difv/difv
        real difhv(nlvdim,1)		!horizontal diffusivity
        common /difhv/difhv
        integer nlvdi,nlv               !total number of levels
        common /level/ nlvdi,nlv
        real rhov(nlvdim,nkndim)	!water density
        common /rhov/rhov
        double precision dmf_mud(nlvdim,nkndim)  !floc diameter
        common /dmf_mud/dmf_mud
        double precision rhomud(nlvdim,nkndim) !Mud floc particle density (kg/m3)
        common /rhomud/rhomud
        real lambda(nlvdim,nkndim) !Mud floc particle density (kg/m3)
        common /lambda/lambda
        integer testnode 
        logical ldebug
        real z0bkmud(nkndim)                   !bottom roughenss on nodes
        common /z0bkmud/z0bkmud

!        real dmf_mud1000(nlvdim,nkndim)

        real wprv(nlvdim,nkndim)
        common /wprv/wprv

        double precision numv(0:nlvdim,nkndim)  !viscosity (momentum)
        double precision nuhv(0:nlvdim,nkndim)  !diffusivity (scalars)
        double precision tken(0:nlvdim,nkndim)  !turbulent kinetic energy
        double precision eps (0:nlvdim,nkndim)  !dissipation rate
        double precision rls (0:nlvdim,nkndim)  !length scale

        common /numv/numv
        common /nuhv/nuhv
        common /tken_gotm/tken
        common /eps_gotm/eps
        common /rls_gotm/rls

        integer ius,id,itmcon,idtcon	!output parameter
        integer iround
        real getpar
  
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
        double precision rhosed		!Mud primary particle density (kg/m3)
        common /rhosed/ rhosed
        double precision cgel(nlvdim),phigel(nlvdim)
        real rowass     !Mud floc particle density (kg/m3)
        common /rowass/ rowass
        double precision mudref		!Initial fluid mud concentration (kg/m3)
        double precision dm0		!Primary particle diameter [m]
        double precision nf(nlvdim)		!Fractal dimension
        common /dm0/ dm0
        double precision visk		!molecular viscosity
        double precision phi(nlvdim),phip(nlvdim)
        integer lthick			!Initial fluid mud layer thickness (layer number)
        real mudc(nlvdim,nkndim)	!Fluid mud concentration array (kg/m3)
        common /mudc/mudc
        real mudhpar			!Fluid mud diffusion coefficient [m**2/s]
        real difmol			!Molecolar diffusion coefficient [m**2/s]
        real mwsink			!Fluid mud floc settling velocity [m/s]
        real wsink(0:nlvdim,nkndim)	!Settling velocity array
        common /wsink/wsink
        real tsec			!Simulation time, real [s]
        real bnd3_mud(nb3dim,0:nbcdim)  !Array containing boundary state
        real bnd3_lam(nb3dim,0:nbcdim)  !Array containing boundary state
        real bnd3_dmf(nb3dim,0:nbcdim)  !Array containing boundary state
        real fmbnd(nsmud),fak,smooth
        real t_now,t_start		        !Boundary vector [kg/m3]
        character*80 mud2dn(1)		!File for boundary condition
        common /mud2dn/mud2dn
        character*80 lam2dn(1)          !File for boundary condition
        common /lam2dn/lam2dn
        character*80 dmf2dn(1)          !File for boundary condition
        common /dmf2dn/dmf2dn
        logical circle,ldumpmud,linitmud
        include 'testbndo.h'

        save /mud2dn/
        save /lam2dn/
        save /dmf2dn/
        save what,mudref,lthick
        save fmbnd
        save bnd3_mud
        save mudhpar,difmol
        save mwsink,visk
        save ius,id,itmcon,idtcon
        save circle
        save icall
        save t_now
        data icall /0/
        data t_now/0./

        ldebug = .true.
        testnode = 2449 
        ldumpmud = .true.
        linitmud = .false.

! ----------------------------------------------------------
! Initialization
! ----------------------------------------------------------
! This section is called only the first time step when ICALL = 0

        t_now = t_now + 1.
        t_start = 0. 
        smooth  = 10.

        if( icall .le. -1 ) return

        imud = iround(getpar('imud'))
        rhow = getpar('rowass')
        if( imud .le. 0 ) then
          icall = -1
          return
        end if

        if( icall .eq. 0 ) then

          itmcon = iround(getpar('itmcon'))
          idtcon = iround(getpar('idtcon'))
          if( it .lt. itmcon ) return
          icall = 1

!         --------------------------------------------------
!         Initializes state variables and constants
!         --------------------------------------------------

          circle = .false.
          !circle = .true.

          what = 'mud'
          nslam = 1
          nsdmf = 1
          mudref = getpar('mudref')
          mudhpar = getpar('mudhpar')
          difmol = getpar('difmol')
          mwsink = getpar('mwsink')
          lthick = getpar('lthick')
          rhosed = getpar('rhosed')
          rhow = getpar('rowass')
          dm0 = getpar('dm0')
          visk = getpar('vismol')

! init mud 
          dmf_mud = 100./10**6.
          nf = 2.
          z0bkmud = 0.
          
          do k = 1,nkn
!         --------------------------------------------------
!         add mud at first time step 
!         --------------------------------------------------

          if (linitmud) then
            lmax = ilhkv(k)
            do l = 1,lmax
              if (l.lt.lthick) then
                mudc(l,k) = 0.
              elseif(l.ge.lthick) then
                mudc(l,k) = mudref
              end if
            end do
          endif
!         --------------------------------------------------
!         Sets sink vel. ... this is overwritten ... 
!         --------------------------------------------------
!            call get_mudrho(k,rhow,nf,dm0,ldebug,testnode)
!         --------------------------------------------------
!         Sets sink vel. ... this is overwritten ... 
!         --------------------------------------------------
!            call get_cgel(k,nf,dm0,cgel,phigel,ldebug,testnode)
!         --------------------------------------------------
!         Sets sink vel. ... this is overwritten ... 
!         --------------------------------------------------
!            call get_mudc(k,rhow,nf,dm0,phi,phip,ldebug,testnode)
!         --------------------------------------------------
!         Sets sink vel. ... this is overwritten ... 
!         --------------------------------------------------
!            call set_hind_wsink(k,nlvdi,nkn,phi,phip,
!     &                         phigel,ldebug,testnode)
            wsink(:,k) = 0.
          end do

          !stop 'check init'

          !wsink = 0.
!         --------------------------------------------------
!         Sets boundary conditions for all state variables
!         --------------------------------------------------
          fmbnd = 0.
          nintp = 2
          call bnds_init(what,mud2dn,nintp,nsmud,nb3dim,bnd3_mud,fmbnd)
               !bnds_init(text,file,nintp,nvar,ndim,array,aconst)
!          call bnds_init(what,lam2dn,nintp,nslam,nb3dim,bnd3_lam,fmbnd)
!          call bnds_init(what,dmf2dn,nintp,nsdmf,nb3dim,bnd3_dmf,fmbnd)
!         --------------------------------------------------
!         Initializes output
!         --------------------------------------------------
          id = 21
          ius = 0
          call confop(ius,itmcon,idtcon,nlv,1,'mud')
          write(6,*) 'fluid mud model initialized...'

        endif ! 1st call 

!        write(*,*) t_now, t_start

        if (ldumpmud) then
          do k = 1,nkn
            lmax = ilhkv(k)
!         --------------------------------------------------
!         add mud at first time step 
!         --------------------------------------------------
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
          end do
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
        do k = 1, nkn 
!          lmax = ilhkv(k)
!          difv(lmax,k) = 0. ! No Diffusion 
!          difv(1,k)=0.
!          mudc(lmax,k)=mudc(lmax-1,k) ! Zero flux ...
!          mudc(1,k)=mudc(2,k)
        end do
!       -------------------------------------------------------------
!       Transport and diffusion of fluid mud
!       -------------------------------------------------------------

        difhv = 0.
        difv  = 0.
        mudhpar = 0.
        wsink = 0. 

          if (circle) then
            call scal_adv_circ(what,ivar
     &                      ,mudc,bnd3_mud
     &                      ,mudhpar,wsink
     &                      ,difhv,difv,0.)
          else
!            call scal_adv(what,ivar
!     &                 ,mudc,bnd3_mud
!     &                 ,mudhpar,wsink
!     &                 ,difhv,difv,0.)

!    subroutine scal_adv_fact(what,ivar,fact
!     +              ,scal,bnd3
!     +              ,rkpar,wsink,wsinkv
!     +                          ,difhv,difv,difmol)

            call scal_adv_fact(what,ivar
     &                 ,mudc,bnd3_mud
     &                 ,mudhpar,wsink
     &                 ,difhv,difv,0.)

          endif
!       -------------------------------------------------------------
!       Boundary condition for fluid mud
!       -------------------------------------------------------------
!        call scal_bnd(what,tsec,bnd3_lam)
!       -------------------------------------------------------------
!       Transport of lambda
!       -------------------------------------------------------------
!          if (circle) then
!            call scal_adv_circ(what,ivar
!     &                      ,lambda,bnd3_mud
!     &                      ,0.,wprv
!     &                      ,0.,0.,0.)
!          else
!            call scal_adv(what,ivar
!     &                 ,lambda,bnd3_mud
!     &                 ,0.,wprv
!     &                 ,0.,0.,0.)
!
!          endif
!       -------------------------------------------------------------
!       Boundary condition for flocs
!       -------------------------------------------------------------
!       call scal_bnd(what,tsec,bnd3_dmf)
!       -------------------------------------------------------------
!       Transport and diffusion of floc diameter
!       -------------------------------------------------------------
!          dmf_mud1000 = real(dmf_mud*1000000.d0)
!          if (circle) then
!            call scal_adv_circ(what,ivar
!     &                      ,dmf_mud,bnd3_mud
!     &                      ,0.,wprv
!     &                      ,0.,0.,0.)
!          else
!            call scal_adv(what,ivar
!     &                 ,dmf_mud,bnd3_mud
!     &                 ,0.,wprv
!     &                 ,0.,0.,0.)
!          endif
!          dmf_mud = real(dmf_mud1000/1000000.d0)
!
        do k = 1, nkn
!       ------------------------------------------------------
!       compute new floc diameter
!        ------------------------------------------------------
!          call get_floc_diam(k,nlvdim,dt,dm0,nf,ldebug,testnode)
!       ------------------------------------------------------
!       Get mud density  
!       -----------------------------------------------------
!          call get_mudrho(k,rhow,nf,dm0,ldebug,testnode)
!       ------------------------------------------------------
!       Get gelling concentration  
!       -----------------------------------------------------
!          call get_cgel(k,nf,dm0,cgel,phigel,ldebug,testnode)
!       ------------------------------------------------------
!       Get mud vol. concentration  
!       ------------------------------------------------------
!          call get_mudc(k,rhow,nf,dm0,phi,phip,ldebug,testnode)
!       -------------------------------------------------------
!       Update settling velocity
!       -------------------------------------------------------
!          call set_hind_wsink(k,nlvdim,nkn,phi,phip,phigel,
!     &                        ldebug,testnode)
          !wsink(:,k) = 0.
        end do ! nodes

!       -------------------------------------------------------------
!       Check mass conservation
!       -------------------------------------------------------------
        call masscons(it,mudc)
        !write(*,*) sum(mudc)
!       -------------------------------------------------------
!       Fill output arrays
!       -------------------------------------------------------
        !call confil(ius,itmcon,idtcon,id,nlvdi,mudc)
! -------------------------------------------------------------
! End of routine
! -------------------------------------------------------------
        end subroutine submud
! -------------------------------------------------------------
       subroutine get_mudc(k,rhow,nf,dm0,phi,phip,ldebug,testnode)

        implicit none

        include 'param.h'

        integer k,l,nlev
        double precision phip(nlvdim),phi(nlvdim)
        double precision rhomud(nlvdim,nkndim)        !Mud floc particle density (kg/m3)
        common /rhomud/ rhomud
        double precision dmf_mud(nlvdim,nkndim)  !floc diameter
        common /dmf_mud/dmf_mud
        real mudc(nlvdim,nkndim)  !floc diameter
        common /mudc/mudc
        double precision rhosed                !Mud primary particle density (kg/m3)
        common /rhosed/ rhosed
        double precision rhow          !Sediment mineral density (kg/m3)
        double precision nf(nlvdim)    !Fractal dimension
        double precision dm0           !Actual and primary diameter

        integer ilhkv(nkndim),lmax           !number of node levels
        common /ilhkv/ilhkv
        integer testnode
        logical ldebug


        lmax=ilhkv(k)

        do l = 1, lmax
          phip(l) = mudc(l,k) / rhosed
          phi(l)  = mudc(l,k) / rhosed * (dmf_mud(l,k)/dm0)**(3-nf(l))
          !phi(l)=((rhosed-rhow)/(rhomud(l,k)-rhow))*(mudc(l,k)/rhosed)
          if (ldebug .and. k == testnode) then
            write(1111,'(A10,I10,6F16.8)') 'get_mudc',l, 
     &                                phip(l),phi(l),
     &                                rhosed, rhow, 
     &                                rhomud(l,k), mudc(l,k)
          endif
        end do

       end subroutine get_mudc

!*********************************************************************
! Get volume concentration in the node. 
! Inverting order to be used inside GOTM
       subroutine get_mudrho(k,rhow,nf,dm0,ldebug,testnode)

        implicit none

        include 'param.h'

        integer k,l,nlev
        double precision rhomud(nlvdim,nkndim)        !Mud floc particle density (kg/m3)
        common /rhomud/ rhomud
        double precision dmf_mud(nlvdim,nkndim)  !floc diameter
        common /dmf_mud/dmf_mud
        double precision rhosed                !Mud primary particle density (kg/m3)
        common /rhosed/ rhosed
        double precision rhow          !Sediment mineral density (kg/m3)
        double precision nf(nlvdim)    !Fractal dimension
        double precision dm0           !Actual and primary diameter
        integer ilhkv(nkndim),lmax           !number of node levels
        common /ilhkv/ilhkv
        integer testnode
        logical ldebug

        lmax=ilhkv(k)

        do l = 1, lmax
          rhomud(l,k) = rhow + (rhosed - rhow) * 
     &                  (dm0/dmf_mud(l,k))**(3.-nf(l))
           if (ldebug .and. k == testnode) then
             write(1112,'(A10,I10,7F15.8)') 'get_mudrho', l, 
     &                  rhow, rhosed, rhomud(l,k),
     &                  dm0, dmf_mud(l,k),nf(l)
           endif
        enddo
       end subroutine
! Get volume concentration in the node. 
! Inverting order to be used inside GOTM
!*********************************************************************
       subroutine get_cgel(k,nf,dm0,cgel,phigel,ldebug,testnode)

        implicit none

        include 'param.h'

        integer k,l,nlev

        double precision rhomud(nlvdim,nkndim)        !Mud floc particle density (kg/m3)
        common /rhomud/ rhomud
        double precision dmf_mud(nlvdim,nkndim)  !floc diameter
        common /dmf_mud/dmf_mud

        double precision rhow          !Sediment mineral density (kg/m3)
        double precision nf(nlvdim)    !Fractal dimension
        double precision dm0           !Actual and primary diameter
        double precision cgel(nlvdim)  ! gel, concentration
        double precision phigel(nlvdim)! gel, concentration
        integer ilhkv(nkndim),lmax           !number of node levels
        common /ilhkv/ilhkv
        double precision rhosed                !Mud primary particle density (kg/m3)
        common /rhosed/ rhosed
        integer testnode
        logical ldebug

        lmax = ilhkv(k)

        do l = 1, lmax
          cgel(l) = rhosed * (dm0/dmf_mud(l,k))**(3.-nf(l))
          phigel(l) = 1. 
           if (ldebug .and. k == testnode) then
          write(1113,'(A10,I10,8F15.8)') 'get_cgel',l,
     &            rhosed,dmf_mud(l,k),dm0,nf(l),
     &           (dmf_mud(l,k)/dm0)**(3.-nf(l)),
     &            cgel(l),phigel(l)
           endif
        end do
        !write(*,*) 'entering get_cgel'
       end
!*********************************************************************
! Check mass conservation of fluid mud in suspension

        subroutine masscons(it,mudc)

        implicit none

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhkv(nkndim)           !number of node levels
        common /ilhkv/ilhkv

        integer it
        real mudc(nlvdim,nkndim)	!Fluid mud concentration array (kg/m3)
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
            write(6,*) "*** nan in mudc..."
            write(6,*) l,k,mudc(l,k)
          end if

	  end do
	end do

        write(74,*)it,totm

	end 

!*********************************************************************
        subroutine set_hind_wsink(k,nlvdi,nkn,phi,phip,phigel,
     &                            ldebug,testnode)

        implicit none

        include 'param.h'

        integer nlvdi,nkn               !number of level and nodes
        double precision d              !grain size
        double precision visk           !molecular viskosity 
        real mudc(nlvdim,nkndim)        !Fluid mud concentration array (kg/m3)
        real getpar,vismol
        double precision rhosed                !Mud primary particle density (kg/m3)
        common /rhosed/ rhosed
        common /mudc/mudc
        real wsink(0:nlvdim,nkndim) !Fluid mud concentration array (kg/m3)
        common /wsink/wsink
        double precision rhomud(nlvdim,nkndim)        !Density of mixture 
        common /rhomud/rhomud
        real rhov(nlvdim,nkndim)        !Density of mixture 
        common /rhov/rhov
        integer ilhkv(nkndim)           !number of node levels
        common /ilhkv/ilhkv
        double precision dmf_mud(nlvdim,nkndim)  !floc diameter
        common /dmf_mud/dmf_mud
        integer l,k,lmax
        double precision :: rhost,nu,dst,ws0,ws,rhow,fak,dmi
        double precision, intent(in):: phi(nlvdim)         !Fluid density (kg/m3)
        double precision, intent(in):: phip(nlvdim)         !Fluid density (kg/m3)
        double precision, intent(in):: phigel(nlvdim)          !Fluid density (kg/m3)
        integer testnode
        logical ldebug
        real smooth,fak2
        rhow = getpar('rowass')
        vismol = getpar('vismol')

        smooth = 0.1

         lmax = ilhkv(k)

         wsink(0,k) = 0.
 
         do l=1,lmax
           rhost = (rhomud(l,k)-rhow)/rhow ! 1.65
           dmi   = dmf_mud(l,k)
           nu    = 1.73E-6 !vismol!max(vismol,visv(l,k)) ! 1.E-6
           dst   = (rhost*9.81/nu**2)**(1./3.)*dmi ! approx. 3.5 for particles of 200micros
           ws0   = 11*nu/dmi*(SQRT(1+.01*dst**3)-1)
           fak   = (1.-min(1.d0,phi(l))*(1.-phip(l)))/(1+2.5*phi(l))
           ws    = ws0 * fak
           fak2   = 0.5d0+0.5d0*tanh((phi(l)-phigel(l))/smooth)
           wsink(l,k) = fak2 * 0. + (1.-fak2)*ws
           if(l.gt.1) wsink(l-1,k) = fak2*wsink(l,k)+
     &                             (1.-fak2)*wsink(l-1,k)
           !if (k==198) write(*,*) (phi(l)-phigel(l)),l
           if((phi(l)-phigel(l)).gt.0.01 .and.l.gt.1)then ! no consolidation model yet ...
             wsink(l,k) = 0.d0
             wsink(l-1,k) = 0.d0
           endif
           if (l==lmax) then 
             wsink(l,k) = 0.
             ws             = 0.
           endif
         end do

         do l = 1, lmax
           if (ldebug .and. k==testnode)then
             write(1114,'(A10,I10,11F16.8)') 'set_wsink',
     &         l,
     &         dmf_mud(l,k),
     &         (((rhosed-rhow)/rhow)*9.81/nu**2)**(1./3.)*dmf_mud(l,k),
     &         11*nu/dmi*(SQRT(1+.01*((((rhosed-rhow)/rhow)*
     &         9.81/nu**2)**(1./3.)*dmf_mud(l,k))**3)-1),
     &         wsink(l,k),
     &         (1.-min(1.d0,phi(l))*(1.-phip(l)))/(1+2.5*phi(l)),
     &         (11*nu/dmi*(SQRT(1+.01*((((rhosed-rhow)/rhow)*
     &         9.81/nu**2)**(1./3.)*dmf_mud(l,k))**3)-1))*
     &         (1.-min(1.d0,phi(l))*(1.-phip(l)))/(1+2.5*phi(l)),
     &         wsink(l,k)*1000.,
     &         phi(l),
     &         phip(l),
     &         phigel(l)
             endif
         end do

        end subroutine
!*****************************************************************

        subroutine stress_mud(nldim,k,ldebug,testnode)

        implicit none

        integer nldim

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        real vts(0:nlvdim,nkndim)
        common /vts/vts
        real rhov(nlvdim,nkndim)
        common /rhov/rhov
        double precision rhomud         ! Mud floc particle density (kg/m3)
        common /rhomud/ rhomud
        real mudc(nlvdim,1)        !Fluid mud concentration array (kg/m3)      !ccf
        common /mudc/mudc
        real shearf2(nlvdim,nkndim)
        common /saux1/shearf2

        integer k,l,nlev
        real*8 aux,dh,du,dv,m2,dbuoy,tau_test,lambda_e
        real*8 rho1, rho2, visk_bar,mu8,mu0,beta,c,tau0
        real h(nlvdim)
        real*8 cnpar,stress_x,stress_y, rhobar, g_dot, rhop, tau,visk
        real*8 g2_dot,g_dot_thr,smooth,viskmax,tf

        integer testnode
        logical ldebug

        if( nldim .ne. nlvdim ) stop 'error stop stress_mud: dimension'

        aux = -grav / rowass
        cnpar = 1

!tau = 786.3684 (d-1)^2.7
!       + 1781.741 (d-1)^2.7 tanh(0.2867 g)
!       + 10.9993 (d-1)^2.7 g

!        do k=1,nkn
          call dep3dnod(k,+1,nlev,h)
          vts(0,k) = 0.
          do l=1,nlev
            rhobar = rhov(l,k)/1000.
            g_dot = sqrt(shearf2(l,k))
            rhop = max(0.d0,rhobar)
            call set_toorman_constants(rhop,mu8,mu0,beta)
            c  = mu0 - mu8
            call set_yieldstress(rhop,tau0)
            g_dot_thr = 0.000001
            smooth    = 2.
            viskmax   = 1000.
            visk      = 0.
            lambda_e  = 0.
            tau       = 0.
            tau_test  = 0.
!Mathieu ... tanh .. taug nix zu große werte bei kleinen g_dot
            tau_test = 786.3684d0*(rhop)**2.7+
     &         1781.741d0*(rhop)**2.7*tanh(max(20.d0,0.2867d0*g_dot))+
     &           10.9993d0*(rhop)**2.7*g_dot 
!Toormann mit Oberrecht 
            lambda_e = 1.d0/(1.d0+beta*g_dot)
            tau  = tau0 + (mu8 + c + beta*tau0*lambda_e) * g_dot
!take care with the 1000
            visk=-beta**2*tau0*g_dot/(1.d0+beta*g_dot)**2+
     &            mu8+c+beta*tau0/(1.+beta*g_dot)
            visk = visk/1000.
            vts(l,k) = min(1000.,visk) 
            
            if (ldebug .and. k == testnode) then
              write(1115,'(A10,I5,17F14.8)') 'mudvisc',l,rhop,
     &                 rhov(l,k),g_dot-g_dot_thr,
     &                 mu0,mu8,beta,c,lambda_e,
     &                 tau0,tau,tau_test,vts(l,k)
            endif
         end do
!          vts(nlev,k) = 0. 
!          vts(1,k)    = 0.
!        end do
      end subroutine
!*****************************************************************
        subroutine set_yield(nldim,k,tstress,ldebug,testnode)

        implicit none

        integer nldim

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        real vts(0:nlvdim,nkndim)
        common /vts/vts
        real rhov(nlvdim,nkndim)
        common /rhov/rhov
        double precision rhomud         ! Mud floc particle density (kg/m3)
        common /rhomud/ rhomud
        real mudc(nlvdim,1)        !Fluid mud concentration array (kg/m3)      !ccf
        common /mudc/mudc
        real shearf2(nlvdim,nkndim)
        common /saux1/shearf2
        real tstress(nlvdim)
        real z0bkmud(nkndim)                   !bottom roughenss on nodes
        common /z0bkmud/z0bkmud
        real z0bk(nkndim)                   !bottom roughenss on nodes
        common /z0bk/z0bk
        real visv(0:nlvdim,nkndim)
        common /visv/visv
        real disv(0:nlvdim,nkndim)
        common /disv/disv
        real visv_yield(0:nlvdim,nkndim)
        common /visv_yield/visv_yield
        real disv_yield(0:nlvdim,nkndim)
        common /disv_yield/disv_yield

        integer k,l,nlev
        real aux,dh,du,dv,m2,dbuoy,tau_test,lambda_e
        real rho1, rho2, visk_bar,mu8,mu0,beta,c,tau0
        real h(nlvdim)
        real cnpar,stress_x,stress_y, rhobar, g_dot, rhop, tau,visk
        real g2_dot,g_dot_thr,smooth,viskmax,tf
        real fak

        integer testnode
        logical ldebug

        smooth = 0.000001

        if( nldim .ne. nlvdim ) stop 'error stop stress_mud: dimension'

          call dep3dnod(k,+1,nlev,h)

          do l=1,nlev
            rhop = max(0.d0,dble(rhov(l,k)))/1000.
            call set_yieldstress(rhop,tau0)
            if (tau0 .gt. 0. .and. rhop .gt. 0.) then
              fak = 0.5d0+0.5d0*tanh((tstress(l)-tau0)/smooth)
              visv_yield(l,k) = fak * visv(l,k) + (1.-fak) * 99.
              disv_yield(l,k) = fak * disv(l,k) + (1.-fak) *  0.
              if (l==nlev) z0bkmud(k)  = fak * z0bk(k)   + (1.-fak) * 99.
            else
              visv_yield(l,k) = visv(l,k)
              disv_yield(l,k) = disv(l,k)
              if (l==nlev) z0bkmud(k) = z0bk(k)
            endif
            if (ldebug .and. k == testnode) then
              write(1117,'(A10,I5,17F14.8)') 'set yield',
     &        l, rhop,tstress(l)-tau0,tstress(l),tau0,
     &        visv_yield(l,k),disv_yield(l,k),
     &        visv(l,k),disv(l,k),z0bk(k),z0bkmud(k),fak, 
     &        tanh((tstress(l)-tau0)/smooth)
            endif
          end do

      end subroutine
!---------------------------------------------------------
!=5000*rhop^3-340*rhop^2+10*rhop
      SUBROUTINE set_yieldstress(rhop,tau_y)

        IMPLICIT NONE
!
        REAL, INTENT(IN)  :: rhop
        REAL, INTENT(OUT) :: tau_y

        tau_y = 5000.d0*rhop**3-340.d0*rhop**2+10.d0*rhop

!        write(*,*) rhopp, tau_y

      END SUBROUTINE
!---------------------------------------------------------
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
        !write(*,'(5F20.10)') rhop, mu8, 1.16*rhop**2.+0.1*rhop

      END SUBROUTINE
!---------------------------------------------------------
        subroutine set_mud_roughness(k,lmax,alpha)

        implicit none

        integer nldim

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw,lmax
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft,alpha
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        real vts(0:nlvdim,nkndim)
        common /vts/vts
        real wsink(0:nlvdim,nkndim)
        common /wsink/wsink
        real rhov(nlvdim,nkndim)
        common /rhov/rhov
        double precision rhomud         ! Mud floc particle density (kg/m3)
        common /rhomud/ rhomud
        real mudc(nlvdim,1)        !Fluid mud concentration array (kg/m3)      !ccf
        common /mudc/mudc
        real shearf2(nlvdim,nkndim)
        common /saux1/shearf2

        integer k,nlev

        real*8 aux,dh,du,dv,m2,dbuoy
        real*8 rho1, rho2, visk_bar,tstress,ufric
        real*8 cnpar,stress_x,stress_y, rhobar
        real*8  g_dot, rhop, tau,visk, drho, ri

        real, parameter :: bpar = 1.
        real, parameter :: beta = 0.7
        real, parameter :: mpar = 1.

        double precision visv(0:nlvdim,nkndim)  !viscosity (momentum)

        common /visv/visv
        integer ilhkv(nkndim)           !number of node levels
        common /ilhkv/ilhkv
        real h(nlvdim)

        call dep3dnod(k,+1,nlev,h)

        lmax = ilhkv(k)
        rho1 = rhov(lmax  ,k)+1000.
        rho2 = rhov(lmax-1,k)+1000.
        rhobar = 0.5 * (rho1 + rho2)
        rhobar = rhov(lmax,k)+1000.
        dh = 0.5 * ( h(lmax) + h(lmax-1) )
        drho = abs((rho2-rho1))/dh
        g_dot = sqrt(shearf2(lmax,k)) 
        if (g_dot .gt. 0.0000001) then
!AR: take care with the mud!
          visk_bar = visv(lmax,k) !+ vts(lmax,k)
          tstress  = visk_bar * g_dot
          ufric = sqrt(tstress)
          ri = 9.81/rhobar*drho/g_dot**2
          alpha = exp(-(1+beta*wsink(lmax,k)/ufric)*
     &         (1-exp(-bpar*ri**mpar)))
        else 
          alpha = 1.
        endif
        alpha = 1.
        !write(*,*) k, alpha, wsink(lmax,k), ri, g_dot**2, ufric
      end subroutine

!---------------------------------------------------------------------------
        subroutine set_bottom_stress(k,tstress)

        implicit none

        integer nldim

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw,lmax
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft,alpha
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        real vts(0:nlvdim,nkndim)
        common /vts/vts
        real wsink(0:nlvdim,nkndim)
        common /wsink/wsink
        real rhov(nlvdim,nkndim)
        common /rhov/rhov
        double precision rhomud         ! Mud floc particle density (kg/m3)
        common /rhomud/ rhomud
        real mudc(nlvdim,1)        !Fluid mud concentration array (kg/m3)      !ccf
        common /mudc/mudc
        real shearf2(nlvdim,nkndim)
        common /saux1/shearf2

        integer k,l,nlev
        real tstress
        real*8 aux,dh,du,dv,m2,dbuoy
        real*8 rho1,rho2,visk_bar,ufric
        real*8 cnpar,stress_x,stress_y, rhobar
        real*8  g_dot, rhop, tau,visk, drho, ri

        real visv(0:nlvdim,nkndim)      !viscosity (momentum)
        common /visv/visv      !viscosity (momentum)
        integer ilhkv(nkndim)           !number of node levels
        common /ilhkv/ilhkv

        lmax = ilhkv(k)
        g_dot = sqrt(shearf2(lmax,k))
        visk_bar = visv(lmax,k) !+ vts(lmax,k)
        tstress  = visk_bar * g_dot 
        !write(*,*) k, visk_bar, tstress, du, dv, dh,lmax
      end subroutine
!---------------------------------------------------------------------------
        subroutine set_bottom_visk(k,ufric,depth,visk1,visk2)

        implicit none

        integer nldim

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw,lmax
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer k,l,nlev
        real*8 ufric,depth
        real visk1,visk2
        integer ilhkv(nkndim)           !number of node levels
        common /ilhkv/ilhkv
        real h(nlvdim)

        call dep3dnod(k,+1,nlev,h)

        lmax = ilhkv(k)

!AR: 2do check dep3nod with respect to nlev
        visk1 = 0.41 * ufric * h(lmax-1)*(1 - h(lmax-1)/depth)
        visk2 = 0.41 * ufric * h(lmax-2)*(1 - h(lmax-2)/depth)
!        write(*,'(I10,5F15.10)')k,visk1,visk2,h(nlev-1),depth,ufric
      end subroutine
!
!---------------------------------------------------------------------------
!
      subroutine get_floc_diam(k,nldim,dt,dm0,nf,ldebug,testnode)
        implicit none

        integer :: k,l,nldim,nlev

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        double precision rhosed                !Mud primary particle density (kg/m3)
        common /rhosed/ rhosed

        integer ilhkv(nkndim)           !number of node levels
        common /ilhkv/ilhkv
        real vts(0:nlvdim,nkndim)
        common /vts/vts
        real visv(0:nlvdim,nkndim)
        common /visv/visv
        real rhov(nlvdim,nkndim)
        common /rhov/rhov
        double precision rhomud         ! Mud floc particle density (kg/m3)
        common /rhomud/ rhomud
        real mudc(nlvdim,1)             !Fluid mud concentration array (kg/m3)      !ccf
        common /mudc/mudc
        double precision dmf_mud(nlvdim,nkndim)  !floc diameter
        common /dmf_mud/dmf_mud
        real shearf2(nlvdim,nkndim)
        common /saux1/shearf2
!        double precision eps (0:nlvdim,nkndim)  !dissipation rate
!        common /eps_gotm/eps

        real, intent(in)    :: dt                   ! time step
        double precision, intent(in)    :: dm0                  ! Size of primary floc
        double precision, intent(out)   :: nf(nlvdim)           ! 3D fractal dimension of floc

        double precision, parameter     :: pi = 3.14159265

        double precision, parameter     :: alpha = 3.0          ! Coefficient
        double precision, parameter     :: p  = 1.0
        double precision, parameter     :: q  = 0.5 


        double precision, parameter     :: Fc  = 2.0      ! Characteristic fractal dimension ! formel4 seite 59 khelifa and hill 
        double precision, parameter     :: dfc = 2000.E-6    ! Characteristic size of floc ! 2microMeter laut formel4 seite 59 khelifa and hill
        double precision, parameter     :: Fy  = 1.E-10   ! 10^-10 Yield strength of floc ( N )

        double precision, parameter     :: ka = 0.98            ! emp. coeff.
        double precision, parameter     :: kb = 3.3E-5      ! emp. coeff.

        double precision, parameter     :: ka2 = 14.6!0.98            ! emp. coeff.
        double precision, parameter     :: kb2 = 14.E3!3.3E-5      ! emp. coeff.

        double precision, parameter     :: dequi = 300E-6

        double precision :: dddt, dddt1, dddt2, dddt3
        double precision :: beta, dv, du, dh, cnpar, g_dot_thr
        double precision :: g_dot,g2_dot,dold,dnew,mu,coeff
        double precision :: growth, decay, rk(3)

        integer testnode
        logical ldebug

        real h(nlvdim)

        g_dot_thr = 0.0000001
!AR: tak care with visv should be mu
        nlev = ilhkv(k)

        do l=1,nlev

           g_dot  = sqrt(shearf2(l,k))
           dold = dmf_mud(l,k)

           g_dot = 0.2

           if (dmf_mud(l,k) .lt. 0.) then
             write(*,'(I10,10F20.10)') l,
     &                                g_dot,
     &                                beta, dm0, visv(l,k),
     &                                dt,dold,dnew
             write(*,*) nf(l), dold

             stop 'neg. diam in dmf_mud rk0'
           endif

           if (dmf_mud(l,k) .ne. dmf_mud(l,k) ) then
             write(*,'(I10,10F15.10)') l, nf(l),dold/dm0,
     &                                g_dot,
     &                                beta, dm0, visv(l,k),
     &                                dt,dold,dnew
             stop 'NaN in dmf_mud rk0'
           endif

           mu = 1.E-3!*(rhov(l,k)+1000.)
           beta = log(Fc/3.d0)/log(dfc/dm0)

           nf(l)= alpha * (dold/dm0)**beta
           coeff=(g_dot*dm0**beta)/(beta*log(dold/dm0)+1.d0)
           growth=0.3333*ka*mudc(l,k)/rhosed*
     &            dm0**(nf(l)-3)*dold**(-nf(l)+4.-beta)
           decay=0.3333*kb*(mu*g_dot/fy)**0.5*
     &        dm0**(-1)*dold**(-beta+2)*(dold-dm0)
!           dddt1 = sign(min(0.1*dequi,abs(coeff*(growth-decay))),
!     &             coeff*(growth-decay)) 
           dddt1 = coeff*(growth-decay)
           dnew = dold + dddt1 * dt
           dold = dnew
           !dmf_mud(l,k) = dnew

           if (dnew .ne. dnew ) then
             write(*,'(I10,10F15.10)') l, nf(l),dold/dm0,
     &                                g_dot,
     &                                beta, dm0, visv(l,k),
     &                                dt,dold,dnew
             stop 'NaN in dmf_mud rk1'
           endif

           if (dnew .lt. 0.) then
             write(*,'(I10,10F20.10)') l,
     &                                g_dot,
     &                                beta, dm0, visv(l,k),
     &                                dt,dold,dnew
             write(*,*) nf(l), dold

             stop 'neg. diam in dmf_mud rk1'
           endif

        
           nf(l)= alpha * (dold/dm0)**beta
           coeff=(g_dot*dm0**beta)/(beta*log(dold/dm0)+1.d0)
           growth=0.3333*ka*mudc(l,k)/rhosed*
     &            dm0**(nf(l)-3)*dold**(-nf(l)+4.-beta)
           decay=0.3333*kb*(mu*g_dot/fy)**0.5*dm0**(-1)*
     &        dold**(-beta+2)*(dold-dm0)
!           dddt2 = sign(min(0.1*dequi,abs(coeff*(growth-decay))),
!     &             coeff*(growth-decay))
           dddt2 = coeff*(growth-decay)
           dnew = 0.75*dmf_mud(l,k)+0.25*dnew+0.25*dddt2*dt
           dold = dnew

           if (dnew .ne. dnew) then
             write(*,'(I10,10F15.10)') l, nf(l),dold/dm0,
     &                                g_dot,
     &                                beta, dm0, visv(l,k),
     &                                dt,dold,dnew

             stop 'NaN in dmf_mud rk2'
           endif

           if (dnew .lt. 0.) then
             write(*,'(I10,10F20.10)') l,
     &                                g_dot,
     &                                beta, dm0, visv(l,k),
     &                                dt,dold,dnew
             write(*,*) nf(l), dold

             stop 'neg. diam in dmf_mud rk2'
           endif

           nf(l)= alpha * (dold/dm0)**beta
           coeff=(g_dot*dm0**beta)/(beta*log(dold/dm0)+1.d0)
           growth=0.3333*ka*mudc(l,k)/rhosed*
     &            dm0**(nf(l)-3)*dold**(-nf(l)+4.-beta)
           decay=0.3333*kb*(mu*g_dot/fy)**0.5*dm0**(-1)*
     &        dold**(-beta+2)*(dold-dm0)
!           dddt3 = sign(min(0.1*dequi,abs(coeff*(growth-decay))),
!     &             coeff*(growth-decay))
           dddt3 = coeff*(growth-decay)
           dnew = 0.3333*dmf_mud(l,k)+0.6666*dnew+0.6666*dddt3*dt
           dmf_mud(l,k) = dnew
           if (dmf_mud(l,k) .ne. dmf_mud(l,k)) then
             write(*,'(I10,10F20.10)') l,
     &                                g_dot, 
     &                                beta, dm0, visv(l,k),
     &                                dt,dold,dnew
             write(*,*) nf(l), dold
           
             stop 'NaN in dmf_mud rk3'
           endif

           if (dmf_mud(l,k) .lt. 0.) then
             write(*,'(I10,10F20.10)') l,
     &                                g_dot,
     &                                beta, dm0, visv(l,k),
     &                                dt,dold,dnew
             write(*,*) nf(l), dold

             stop 'neg. diam in dmf_mud rk3'
           endif


           if (ldebug .and. k == testnode) then
             write(1116,'(A10,I10,10F20.10)') 'get_floc', l, 
     &                  dddt2*1.E6,dm0*1.E6,dold*1.E6,dnew*1.E6, 
     &                  growth*1.E6,decay*1.E6,(growth-decay)*1.E6,
     &                  beta, nf(l), g_dot
           endif                                                                

        end do

        nf(nlev) = nf(nlev-1)
        dmf_mud(nlev,k) = dmf_mud(nlev-1,k)

        end subroutine get_floc_diam

!******************************************************************
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
!******************************************************************
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
!**********************************************************************
!*                                                                    *
!**********************************************************************
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

!*********************************************************************

        subroutine scal_adv_circ(what,ivar
     +                          ,scal,bnd3
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

c shell for scalar (for parallel version)

        implicit none

        include 'param.h'

        character*(*) what
        integer ivar
        real scal(nlvdim,nkndim)
        real bnd3(nb3dim,0:nbcdim)

        real rkpar
        real wsink(0:nlvdim,nkndim)
        real difhv(nlvdim,1)
        real difv(0:nlvdim,1)
        real difmol

        real cobs(nlvdim,1) !observations (for nudging)
        real robs       !use nudging


        real bnd3_aux(nb3dim)
        real r3v(nlvdim,nkndim)

        integer iwhat,ichanm
        character*10 whatvar,whataux

        cobs = 0.
        robs = 0.

c--------------------------------------------------------------
c make identifier for variable
c--------------------------------------------------------------

        whatvar = what
        if( ivar .ne. 0 ) then
          write(whataux,'(i2)') ivar
          whatvar = what // whataux
        end if
        iwhat = ichanm(whatvar)

c--------------------------------------------------------------
c transfer boundary conditions of var ivar to 3d matrix r3v
c--------------------------------------------------------------

        call bnds_trans(whatvar(1:iwhat)
     +                          ,nb3dim,bnd3,bnd3_aux
     +                          ,ivar,nlvdim,r3v)

c--------------------------------------------------------------
c adjust boundary conditions for the flume case in order that 
c concentration going out from the system is imposed in the 
c corresponding influx nodes
c--------------------------------------------------------------

        call circflume(r3v,scal)

c--------------------------------------------------------------
c do advection and diffusion
c--------------------------------------------------------------

        call scal3sh(whatvar(1:iwhat)
     +                          ,scal,nlvdim
     +                          ,r3v,cobs,robs
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

        end

!*********************************************************************
! This routine sets the boundary condition for the fluid mud 
! concentration for the case of anular flume in the way that conc
! going out from the system is imposed in the influx nodes

        subroutine circflume(r3v,scal)

        implicit none

        include 'param.h'

        real r3v(nlvdim,1)              !boundary conc array
        real scal(nlvdim,1)             !conc array

        integer irv(1)                  !boundary node number array
        common /irv/irv
        integer nlvdi,nlv               !total number of levels
        common /level/ nlvdi,nlv
        real mfluxv(nlvdim,1)
        common /mfluxv/mfluxv

        integer nbnode                  !Number of boundary nodes
        !parameter (nbnode=17)          !circular flume
        parameter (nbnode=21)           !straight flume

        integer innode(nbnode)          !in boundary node number array (bound1)
        integer ounode(nbnode)          !out boundary node number array (bound2)
        integer ibc,knode
        integer kranf,krend
        integer kin,kou
        integer k,n,l,nk
        integer icall                   !Initialization parameter
        save icall
        save innode,ounode,nk
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

        end if

!       -------------------------------------------------------------------
!       Normal call:
!       Set boundary conc for bound 1 equal to conc in bound 2 (invert order)
!       -------------------------------------------------------------------
        do k = 1,nk
          kin = innode(k)
          kou = ounode(nk-k+1)
          do l = 1,nlvdi
            r3v(l,kou) = scal(l,kou)
            r3v(l,kin) = scal(l,kou)
            !r3v(l,kin) = scal(l,kou)*abs(mfluxv(l,kou)/mfluxv(l,kin)) !corrected
          end do
        end do

      end
!*********************************************************************
      subroutine write_xfn()
        implicit none

        integer :: k,l,nldim

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw,nlvdi
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw,nlvdi
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        double precision rhosed                !Mud primary particle density (kg/m3)
        common /rhosed/ rhosed

        real vts(0:nlvdim,nkndim)
        common /vts/vts
        real visv(0:nlvdim,nkndim)
        common /visv/visv
        real rhov(nlvdim,nkndim)
        common /rhov/rhov
        double precision rhomud         ! Mud floc particle density (kg/m3)
        common /rhomud/ rhomud
        real mudc(nlvdim,1)             !Fluid mud concentration array (kg/m3)      !ccf
        common /mudc/mudc
        double precision dmf_mud(nlvdim,nkndim)  !floc diameter
        common /dmf_mud/dmf_mud
        real shearf2(nlvdim,nkndim)
        common /saux1/shearf2
        integer nen3v(3,1)
        common /nen3v/nen3v
        real xgv(1),ygv(1) ! coordinate nodi
        common /xgv/xgv, /ygv/ygv
	    real uprv(nlvdim,nkndim)
        common /uprv/uprv
   	    real vprv(nlvdim,nkndim)
        common /vprv/vprv
        real hkv(nkndim)
        common /hkv/hkv

        logical, save :: lfirst
        data lfirst/.true./
 
        real h(nlvdim) 
        real hh(nkndim),dep(nkndim)
        integer ie,ip,il,iwrite,icycle
        real,save :: time
        save iwrite
        data time/0./
        data iwrite/100/

        icycle = 1

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
          OPEN(5003, FILE  = 'ergmud.bin'  , FORM = 'UNFORMATTED')
          OPEN(5004, FILE  = 'ergmudc.bin'  , FORM = 'UNFORMATTED')
          lfirst = .false.
        endif

        !write(*,*) mod(iwrite,icycle), iwrite, icycle

        if (mod(iwrite,icycle) .eq. 0) then       
          do il = 1, nlvdi
            time = time + 1.
            WRITE(5002) TIME
            WRITE(5002) (uprv(il,ip), vprv(il,ip), hh(ip),ip = 1, nkn)
            CALL FLUSH(5002)
            WRITE(5003) TIME
            WRITE(5003) (uprv(il,ip),vprv(il,ip),rhov(il,ip),ip=1,nkn)
            CALL FLUSH(5003)
            WRITE(5004) TIME
            WRITE(5004) (uprv(il,ip),vprv(il,ip),mudc(il,ip),ip=1,nkn)
            CALL FLUSH(5004)
!            write(*,*) time, nlvdi
          enddo
        endif
        iwrite = iwrite + 1
  
      end 
!*********************************************************************
