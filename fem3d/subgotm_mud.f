
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

c gotm module
c
c contents :
c
c revision log :
c
c 10.08.2003	ggu	call gotm_init
c 05.10.2004	ggu	new administration routine turb_closure, Munk-And.
c 23.03.2006	ggu	changed time step to real
c 20.11.2006	ggu	new version of keps, for gotm changed common blocks
c 20.10.2007	ccf	new version of gotm (4.0.0)
c 10.04.2008	ggu	integrated in main branch
c 18.09.2008	ccf	bug fix for m2 in setm2n2
c 02.12.2008	ggu	bug in gotm_init: no limiting values for initialization
c 18.12.2008	ggu	bug in GOTM module and setm2n2() corrected
c 16.02.2011	ggu	write n2max to info file, profiles in special node
c 26.06.2012	ggu	changed VERS_6_1_55
c 29.03.2013	ggu	avoid call to areaele -> ev(10,ie)
c 25.10.2013	ggu	changed VERS_6_1_68
c 18.06.2014	ggu	changed VERS_6_1_77
c 19.12.2014	ggu	changed VERS_7_0_10
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 25.05.2016	ggu	changed VERS_7_5_10
c 16.02.2019	ggu	changed VERS_7_5_60
c
c**************************************************************

	subroutine turb_closure

c administers turbulence closure

	implicit none

	real getpar,areaele

	integer iturb
	save iturb
	data iturb / 0 /

	if( iturb .lt. 0 ) return

	if( iturb .eq. 0 ) then
	  iturb = nint(getpar('iturb'))
	  if( iturb .le. 0 ) iturb = -1
	  if( iturb .lt. 0 ) return
	  write(*,*) 'starting turbulence model: ',iturb
	end if

	if( iturb .eq. 1 ) then		!Gotm
	  call gotm_shell
	else if( iturb .eq. 2 ) then	!keps
	  call keps_shell
	else if( iturb .eq. 3 ) then	!Munk Anderson
	  call munk_anderson_shell
	else
	  write(6,*) 'Value iturb not possible: ',iturb
	  stop 'error stop turb_closure: iturb'
	end if

	end

c**************************************************************

	subroutine munk_anderson_shell

c computes turbulent quantities with Munk - Anderson model

	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	include 'femtime.h'
	include 'pkonst.h'

c---------------------------------------------------------------
c aux arrays superposed onto other aux arrays
c---------------------------------------------------------------

	real shearf2(nlvdim,nkndim)
	real buoyf2(nlvdim,nkndim)
	real richard(nlvdim,nkndim)


	integer k,l
	integer nlev
	integer mode
	real h(nlvdim)
	real ri,vis,dif
	real diftur,vistur
	real a,b,alpha,beta

	real getpar

	integer icall
	save icall
	data icall / 0 /

c------------------------------------------------------
c initialization
c------------------------------------------------------

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  write(*,*) 'starting Munk Anderson turbulence model'
	end if

	icall = icall + 1

c------------------------------------------------------
c set up buoyancy frequency and shear frequency
c------------------------------------------------------

	call setm2n2(nlvdim,buoyf2,shearf2)

c------------------------------------------------------
c set up parameters
c------------------------------------------------------

	mode = +1
	a = 10.
	b = 3.33
	alpha = -1./2.
	beta = -3./2.
	vistur = getpar('vistur')
	diftur = getpar('diftur')

c------------------------------------------------------
c compute richardson number for each water column
c------------------------------------------------------

	do k=1,nkn

	    call dep3dnod(k,mode,nlev,h)

	    do l=1,nlev-1
	      ri = buoyf2(l,k) / shearf2(l,k)
	      vis = vistur*(1.+a*ri)**alpha
	      dif = diftur*(1.+b*ri)**beta

	      visv(l,k) = vis
	      difv(l,k) = dif
	      richard(l,k) = ri
	    end do

            if( nlev .eq. 1 ) then
	      visv(0,k) = vistur
	      difv(0,k) = diftur
	      visv(1,k) = vistur
	      difv(1,k) = diftur
	      richard(1,k) = 0.
            else
	      visv(0,k) = visv(1,k)
	      difv(0,k) = difv(1,k)
	      visv(nlev,k) = visv(nlev-1,k)
	      difv(nlev,k) = difv(nlev-1,k)
	      richard(nlev,k) = richard(nlev-1,k)
            end if

	end do

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	end

c**************************************************************

	subroutine gotm_shell

c computes turbulent quantities with GOTM model

	use mod_meteo
	use mod_sinking
	use mod_fluidmud
	use mod_gotm_aux
	use mod_ts
	use mod_roughness
	use mod_diff_visc_fric
	use mod_hydro_print
	use levels
	use basin

	implicit none

	include 'param.h'

	integer ndim
	parameter(ndim=nlvdim)

	double precision dt
	double precision u_taus,u_taub

	double precision hh(0:ndim)
	double precision ddepth(nlvdim)
	double precision plots(nlvdim,2)
	double precision nn(0:ndim), ss(0:ndim)
	double precision richgrad(0:ndim), richflux(0:ndim)

	double precision num(0:ndim), nuh(0:ndim)
	double precision ken(0:ndim), dis(0:ndim)
	double precision len(0:ndim)

	double precision num_old(0:ndim), nuh_old(0:ndim)
	double precision ken_old(0:ndim), dis_old(0:ndim)
	double precision len_old(0:ndim)

 
 

	include 'femtime.h'
	include 'pkonst.h'


c---------------------------------------------------------------
c aux arrays superposed onto other aux arrays
c---------------------------------------------------------------

	real taub(nkndim)
	real areaac(nkndim)



	integer ioutfreq,ks
	integer i, k,l
	integer laux
	integer nlev
	real g
	real czdef,taubot
	save czdef
	real vismol 
	save vismol
	double precision nu

	real h(nlvdim)
	double precision depth		!total depth [m]
	double precision z0s,z0b	!surface/bottom roughness length [m]
	double precision rlmax
	
	integer nltot
	logical bwrite

	real charnock_val		!emp. Charnok constant (1955)
	parameter(charnock_val=1400.)	!default value = 1400.
	double precision z0s_min	!minimum value of z0s
	parameter(z0s_min=0.02)
	real ubot,vbot,rr


	double precision tmp(nlvdim)
	real dtreal
	real getpar,areaele
	real*8 rhop,tau0


! fluid mud variables
	real rig(nlvdim), rif(nlvdim)
      real taux(nlvdim), tauy(nlvdim), tautot(nlvdim), visk_bottom1
      real visk_bottom2
	double precision tau_y
	real tstress(nlvdim)
	real scn			! Schmidt number (=0.5)
	integer imud
	logical bmud
	save bmud

	integer namlst
	character*80 fn	
	integer ius,id,itmcon,idtcon    !output parameter
	save ius,id,itmcon,idtcon
	integer icall, iwrite, icycle, testnode,lmudvisc
	save icall, iwrite, icycle, testnode,lmudvisc
	data icall / 0 /
	data iwrite / 0 /
	logical, save :: ldebug


c------------------------------------------------------
c documentation
c------------------------------------------------------

c total stress :  tau_x = a |u| u   tau_y = a |u| v
c tau_tot^2 = a^2 |u|^2 ( u^2 + v^2 ) a^2 |u|^4
c tau_tot = sqrt(tau_x^2+tau_y^2) = a |u|^2
c
c u_taus=sqrt(sqrt(tx*tx+ty*ty))	(friction velocity)

c------------------------------------------------------
c initialization
c------------------------------------------------------

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  czdef = getpar('czdef')
	  vismol = getpar('vismol')
	  lmudvisc = getpar('lmudvisc')
	  write(*,*) 'starting GOTM turbulence model'
      imud = getpar('imud')
	  bmud = imud .gt. 0
      ldebug   = nint(getpar('ldebug')) .ne. 0
      testnode = getpar('testnode')
      icycle   = getpar('icycle')
c         --------------------------------------------------------
c         Initializes gotm arrays 
c         --------------------------------------------------------

	  call gotm_init

c         --------------------------------------------------------
c         Get file name containing GOTM turbulence model parameters
c         --------------------------------------------------------

          call getfnm('gotmpa',fn)

	  call init_gotm_turb(10,fn,ndim)

	  id = 22
          ius = 0
          call confop(ius,itmcon,idtcon,nlv,1,'lam')

      call write_xfn

	end if

	icall = icall + 1

	call get_timestep(dtreal)
	dt = dtreal
	g = grav

c------------------------------------------------------
c set up bottom stress on nodes
c------------------------------------------------------

	call bnstress(czdef,taub,areaac)

c------------------------------------------------------
c set up buoyancy frequency and shear frequency
c------------------------------------------------------

	call setm2n2(nlvdim,buoyf2,shearf2)

c------------------------------------------------------
c call gotm for each water column
c------------------------------------------------------

	rlmax = 0.
	nltot = 0
    
	do k=1,nkn

	    call dep3dnod(k,+1,nlev,h)

            if( nlev .eq. 1 ) goto 1

c           ------------------------------------------------------
c           update boyancy and shear-frequency vectors
c           ------------------------------------------------------

	    do l=1,nlev-1
	      laux = nlev - l
	      nn(laux) = buoyf2(l,k)
	      ss(laux) = shearf2(l,k)
	    end do
	    nn(0) = 0.
	    nn(nlev) = 0.
	    ss(0) = ss(1)
	    ss(nlev) = ss(nlev-1)

c       ------------------------------------------------------
c       compute layer thickness and total depth
c       ------------------------------------------------------

            depth = 0.
            ddepth = 0.
            do l=1,nlev
              hh(nlev-l+1) = h(l)
              depth = depth + h(l)
              ddepth(l) = depth 
            end do
            do l =1,nlev
              ddepth(l) = depth - ddepth(l) 
            end do

c           ------------------------------------------------------
c           compute surface friction velocity (m/s)
c           ------------------------------------------------------

	    u_taus = sqrt( sqrt( tauxnv(k)**2 + tauynv(k)**2 ) )

	    z0s = charnock_val*u_taus**2/g
	    z0s = max(z0s,z0s_min)

c           ------------------------------------------------------
c           compute bottom friction velocity (m/s)
c           ------------------------------------------------------

	    z0b = z0bk(k)
	    u_taub = sqrt( taub(k) )

	    ubot = uprv(nlev,k)
  	    vbot = vprv(nlev,k)
	    rr = 0.4/(log((z0b+hh(1)/2)/z0b))
            u_taub = rr*sqrt( ubot*ubot + vbot*vbot )

c           ------------------------------------------------------
c           update 1-dimensional vectors
c           ------------------------------------------------------

	    do l=0,nlev
	      num(l) = numv_gotm(l,k)
	      nuh(l) = nuhv_gotm(l,k)
	      ken(l) = tken_gotm(l,k)
	      dis(l) = eps_gotm(l,k)
	      len(l) = rls_gotm(l,k)
	      num_old(l) = numv_gotm(l,k)
	      nuh_old(l) = nuhv_gotm(l,k)
	      ken_old(l) = tken_gotm(l,k)
	      dis_old(l) = eps_gotm(l,k)
	      len_old(l) = rls_gotm(l,k)
	    end do

	    !call save_gotm_init

c           ------------------------------------------------------
c           call GOTM turbulence routine
c           ------------------------------------------------------

 	    call do_gotm_turb   (
     &				  nlev,dt,depth
     &				 ,u_taus,u_taub
     &				 ,z0s,z0b,hh
     &	                         ,nn,ss
     &				 ,num,nuh,ken,dis,len
     &				 )

!           ------------------------------------------------------
!           update viscosity from fluid mud
!           ------------------------------------------------------
        !bmud = .false.
        !bmud = .true. !DEB
        if ( bmud ) then
!         ------------------------------------------------------
!         compute effective viscosity for sediment rheological stresses
!         ------------------------------------------------------
!AR: take care nlev choosen maybe inconsistent with moving bottom
          call set_bottom_visk(k,u_taub,depth,visk_bottom1,visk_bottom2) 
          !write(*,*) visk_bottom1,visk_bottom2,num(0),num(1)
          num(0) = max(num(0),dble(visk_bottom1))
          num(1) = max(num(1),dble(visk_bottom2))
          nuh(0) = num(0)/0.7
          nuh(1) = num(1)/0.7
          if ( lmudvisc .ne. 0 ) then
		call stress_mud(nlvdim,k,ldebug,testnode,icycle)
	  end if
        end if
c           ------------------------------------------------------
c           copy back to node vectors
c           ------------------------------------------------------

	    bwrite = .false.

	    do l=0,nlev
	      numv_gotm(l,k) = num(l)
	      nuhv_gotm(l,k) = nuh(l)
	      tken_gotm(l,k) = ken(l)
	      eps_gotm(l,k)  = dis(l)
	      rls_gotm(l,k)  = len(l)
	      rlmax = max(rlmax,len(l))	!ggu
	      if( len(l) .gt. 100. ) then
		    nltot = nltot + 1
		    bwrite = .true.
	      end if
	    end do

	    if( bwrite ) then

	      !write(45,*) it,k
	      !call save_gotm_write

	      write(89,*) '==========================='
	      write(89,*) 'time___: ',it
	      write(89,*) 'node___: ',k
	      write(89,*) 'nlev___: ',nlev
	      write(89,*) 'dt_____: ',dt
	      write(89,*) 'depth__: ',depth
	      write(89,*) 'utaus_b: ',u_taus,u_taub
	      write(89,*) 'z0s_z0b: ',z0s,z0b
	      write(89,*) 'hh_____: ',(hh     (l),l=0,nlev)
	      write(89,*) 'nn_____: ',(nn     (l),l=0,nlev)
	      write(89,*) 'ss_____: ',(ss     (l),l=0,nlev)
	      write(89,*) 'num_old: ',(num_old(l),l=0,nlev)
	      write(89,*) 'nuh_old: ',(nuh_old(l),l=0,nlev)
	      write(89,*) 'ken_old: ',(ken_old(l),l=0,nlev)
	      write(89,*) 'dis_old: ',(dis_old(l),l=0,nlev)
	      write(89,*) 'len_old: ',(len_old(l),l=0,nlev)
	      write(89,*) 'num____: ',(num    (l),l=0,nlev)
	      write(89,*) 'nuh____: ',(nuh    (l),l=0,nlev)
	      write(89,*) 'ken____: ',(ken    (l),l=0,nlev)
	      write(89,*) 'dis____: ',(dis    (l),l=0,nlev)
	      write(89,*) 'len____: ',(len    (l),l=0,nlev)
	      write(89,*) '==========================='
	    end if

c           ------------------------------------------------------
c           update viscosity and diffusivity
c           ------------------------------------------------------

	    do l=0,nlev
	      laux = nlev - l
	      visv(l,k) = num(laux)
	      difv(l,k) = nuh(laux)
	    end do

        if ( bmud ) then
          call stress3d(nlvdim,k,tstress,ldebug,testnode,icycle)
          call richardson(k,nlvdim,rif,rig,ldebug,testnode,icycle)
!          call set_yield(nlvdim,k,tstress,ldebug,testnode,icycle)
        end if

        if (k == testnode .and. ldebug .and. 
     &            mod(iwrite,icycle) .eq. 0) then
          do i = 1, nlev
            if (.false.) then
          write(*,'(I10,7F15.8,7F15.10)') i,
     &                             ddepth(i),
     &                             ken(nlev-i+1),
     &                             dis(nlev-i+1),
     &                             visv(i,k),
     &                             difv(i,k),
     &                             rhov(i,k),
     &                             -wsinkv(i,k)*1000,
     &                             vts(i,k),
     &                             tstress(i),
     &                             mudc(i,k),
     &                             uprv(i,k),
     &                             dmf_mud(i,k)*1000,
     &                             rig(i),
     &                             rif(i)
            endif
          end do
        end if

          !write(*,*) '----------------------------------'
          !write(*,*) testnode,ldebug,iwrite,icycle,mod(iwrite,icycle) 

          if (k == testnode .and. ldebug .and. 
     &            mod(iwrite,icycle) .eq. 0) then
             do i = 1, nlev 
               !write(*,*) testnode, i, k, wsinkv(i,k), 'from gotm'
          write(1020,'(I10,7F15.8,5F15.8,2E15.6)') i,
     &                             ddepth(i),
     &                             ken(nlev-i+1),
     &                             dis(nlev-i+1),
     &                             num(nlev-i+1),
     &                             nuh(nlev-i+1),
     &                             rhov(i,k),
     &                             -wsinkv(i,k)*1000,
     &                             vts(i,k),
     &                             tstress(i),
     &                             mudc(i,k),
     &                             uprv(i,k),
     &                             dmf_mud(i,k)*1000,
     &                             rig(i),
     &                             rif(i)
              end do
          end if! testnode
          if (k == testnode) iwrite = iwrite + 1

    1     continue
	end do

	call checka(nlvdim,shearf2,buoyf2,taub)

	if (mod(iwrite,icycle) .eq. 0) then
		call write_xfn
		write(*,*) '********* XFN OUTPUT **********',iwrite,icycle
	endif

	!write(*,*) sum(wsinkv), 'from gotm'
 
	!write(70,*) 'rlmax: ',it,rlmax,nltot

!	ks = 0			!internal node number
!	ioutfreq = 3600		!output frequency
!	if( ks .gt. 0 .and. mod(it,ioutfreq) .eq. 0 ) then
!	  write(188,*) it,nlev,(visv(l,ks),l=1,nlev)
!	  write(189,*) it,nlev,(difv(l,ks),l=1,nlev)
!	end if

c------------------------------------------------------
c end of routine
c------------------------------------------------------

	end

c**************************************************************

	subroutine gotm_init

c initializes gotm arrays

	use mod_gotm_aux
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'





	integer l,k
        double precision num_min, nuh_min
        double precision tken_min, eps_min, rls_min

        num_min = 1.e-6
        nuh_min = 1.e-6
        tken_min = 1.e-10
        eps_min = 1.e-12
        rls_min = 1.e-10

        do k=1,nkn
          do l=0,nlv
            numv_gotm(l,k) = num_min
            nuhv_gotm(l,k) = nuh_min
            tken_gotm(l,k) = tken_min
            eps_gotm(l,k)  = eps_min
            rls_gotm(l,k)  = rls_min
          end do
        end do

	end

c**************************************************************

	subroutine gotm_get(k,nlev,num,nuh,tk,ep,rl)

c returns internal parameters from turbulence closure

	use mod_gotm_aux

	implicit none

	integer k		!node
	integer nlev		!number of levels to return
	real num(0:nlev)	!viscosity
	real nuh(0:nlev)	!diffusivity
	real tk (0:nlev)	!kinetic energy
	real ep (0:nlev)	!dissipation
	real rl (0:nlev)	!length scale

	include 'param.h'




	integer l,laux

        do l=0,nlev
	  laux = nlev - l
          num(l) = numv_gotm(laux,k)
          nuh(l) = nuhv_gotm(laux,k)
          tk(l)  = tken_gotm(laux,k)
          ep(l)  = eps_gotm (laux,k)
          rl(l)  = rls_gotm (laux,k)
        end do

	end

c**************************************************************

      subroutine Yevol(Nmx,dt,h,nuh,difmol,Qsour,Yvar,Yd,
     &                 Bcup,Yup,Bcdw,Ydw,Taur)

      implicit none

	integer Nmx		!number of levels
	real dt			!time step [s]
	real h(Nmx)		!level thickness [m]
	real nuh(0:Nmx)		!turbulent diffusivity [m**2/s]
	real difmol		!molecular diffusivity [m**2/s]
	real Qsour(Nmx)		!internal source term
	real Yvar(Nmx)		!variable to diffuse
	real Yd(Nmx)		!relaxation value for variable
	integer Bcup		!type of boundary condition (upper)
	real Yup		!value for boundary condition (upper)
	integer Bcdw		!type of boundary condition (lower)
	real Ydw		!value for boundary condition (lower)
	real Taur		!time scale to use with Yd [s]

c if Taur <= 0	-> do not use Yd
c Bcup/Bcdw :   1 = Neuman condition   2 = Dirichlet condition
c Yup,Ydw   :   value to use for boundary condition
c
c example: Bcup = 1  and Yup = 0  --> no flux condition across surface

!     This subroutine computes the evolution of any 
!     state variable "Y" after Mixing/Advection/Source
!     - with Neuman    Boundary Condition: Bc=1
!     - with Dirichlet Boundary Condition: Bc=2
!     Convention: fluxex are taken positive upward

	integer MaxNmx
	parameter(MaxNmx = 500)

	integer i
	double precision au(0:MaxNmx),bu(0:MaxNmx)
	double precision cu(0:MaxNmx),du(0:MaxNmx)
	double precision avh(0:MaxNmx),Y(0:MaxNmx)
	double precision cnpar,a,c

	if( Nmx .gt. MaxNmx ) stop 'error stop Tridiagonal: MaxNmx'

      cnpar = 1.      		!Crank Nicholson parameter (implicit)

!
! copy to auxiliary values
!

	 Y(0) = 0.
	 avh(0) = 0.
	 do i=1,Nmx
	   Y(i) = Yvar(i)
	   avh(i) = nuh(i) + difmol
	 end do

! Water column 
! => [a b c].X=[d] where X=[2, ...,i,i+1,...Nmx-1]

         do i=2,Nmx-1
          c    =2*dt*avh(i)  /(h(i)+h(i+1))/h(i) 
          a    =2*dt*avh(i-1)/(h(i)+h(i-1))/h(i)
          cu(i)=-cnpar*c                                 !i+1,n+1
          au(i)=-cnpar*a                                 !i-1,n+1
          bu(i)=1-au(i)-cu(i)                            !i  ,n+1
          du(i)=Y(i)+dt*Qsour(i)                         !i  ,n
     &          +(1-cnpar)*(a*Y(i-1)-(a+c)*Y(i)+c*Y(i+1))
         end do

! Surface 
! => [a b /].X=[d] where X=[Nmx-1,Nmx,/]

         if (Bcup.eq.1) then                       !BC Neuman
          a      =2*dt*avh(Nmx-1)/(h(Nmx)+h(Nmx-1))/h(Nmx)
          au(Nmx)=-cnpar*a
          bu(Nmx)=1-au(Nmx)
          du(Nmx)=Y(Nmx)
     &            +dt*(Qsour(Nmx)-Yup/h(Nmx))   
     &            +(1-cnpar)*a*(Y(Nmx-1)-Y(Nmx))
         else if (Bcup.eq.2) then                    !BC Dirichlet
          au(Nmx)=0.
          bu(Nmx)=1.
          du(Nmx)=Yup
         end if

! Bottom  
! => [/ b c].X=[d] where X=[/,1,2]

         if (Bcdw.eq.1) then                        !BC Neuman              
          c    =2*dt*avh(1)/(h(1)+h(2))/h(1)
          cu(1)=-cnpar*c 
          bu(1)=1-cu(1)
          du(1)=Y(1)             
     &          +dt*(Qsour(1)+Ydw/h(1))
     &          +(1-cnpar)*c*(Y(2)-Y(1))
         else if (Bcdw.eq.2) then                    !BC Dirichlet
          cu(1)=0.
          bu(1)=1.
          du(1)=Ydw
         end if
!
! Implicit internal restoring
!
         if ( Taur .gt. 0. .and. Taur .lt. 1.E10 ) then
          do i=1,Nmx
           bu(i)=bu(i)+dt/Taur
           du(i)=du(i)+dt/Taur*Yd(i)
          end do
         end if
!
! Implicit vertical mixing
!

         call Tridiagonal(Nmx,1,Nmx,au,bu,cu,du,Y)

	 do i=1,Nmx
	   Yvar(i) = Y(i)
	 end do

      return
      end

c**************************************************************

      subroutine Tridiagonal(Nmx,fi,lt,au,bu,cu,du,value)
      implicit none

	integer MaxNmx
	parameter(MaxNmx=500)
!
! !USES:
!
! !INPUT PARAMETERS:
      integer Nmx
      integer fi,lt
      double precision au(0:Nmx),bu(0:Nmx),cu(0:Nmx),du(0:Nmx)
!
! !OUTPUT PARAMETERS:
      double precision value(0:Nmx)
!
! !DESCRIPTION:
!     Used trough out the program.
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  25Nov98   The GOTM developers  Initial code.
!
! !LOCAL VARIABLES:
      double precision ru(0:MaxNmx),qu(0:MaxNmx)
      integer i
!
!EOP
!BOC
	if( Nmx .gt. MaxNmx ) stop 'error stop Tridiagonal: MaxNmx'

      ru(lt)=au(lt)/bu(lt)
      qu(lt)=du(lt)/bu(lt)

      do i=lt-1,fi+1,-1
         ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
         qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
      end do

      qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

      value(fi)=qu(fi)
      do i=fi+1,lt
         value(i)=qu(i)-ru(i)*value(i-1)
      end do

      return
      end
!EOC

c**************************************************************

	subroutine setm2n2(nldim,buoyf2,shearf2)

c sets buoyancy and shear frequency
c
c this is done directly for each node
c in rhov is already rho^prime = rho - rho0 (deviation)
!  Discretisation of vertical shear squared according to Burchard (2002)
!  in order to guarantee conservation of kinetic energy when transformed
!  from mean kinetic energy to turbulent kinetic energy.
c
c bug fix in computation of shearf2 -> abs() statements to avoid negative vals

	use mod_ts
	use mod_hydro_print
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nldim
	real buoyf2(nldim,1)
	real shearf2(nldim,1)

	include 'param.h'

	include 'pkonst.h'
	include 'femtime.h'



	integer k,l,nlev
	real aux,dh,du,dv,m2,dbuoy
	real h(nlvdim)
	real cnpar			!numerical "implicitness" parameter
	real n2max,n2
	real nfreq,nperiod

	integer iuinfo
	save iuinfo
	data iuinfo / 0 /
 
	if( nldim .ne. nlvdim ) stop 'error stop setbuoyf: dimension'

        aux = -grav / rowass
	cnpar = 1
	n2max = 0.
 
        do k=1,nkn
          call dep3dnod(k,+1,nlev,h)

          do l=1,nlev-1
            dh = 0.5 * ( h(l) + h(l+1) )
            dbuoy = aux * ( rhov(l,k) - rhov(l+1,k) )
            n2 = dbuoy / dh
			n2max = max(n2max,n2)
            buoyf2(l,k) = n2

            du = 0.5*(
     &       (cnpar*abs((uprv(l+1,k)-uprv(l,k))
     &		*(uprv(l+1,k)-upro(l,k)))+
     &       (1.-cnpar)*abs((upro(l+1,k)-upro(l,k))
     &		*(upro(l+1,k)-uprv(l,k))))
     &       /dh/h(l)
     &      +(cnpar*abs((uprv(l+1,k)-uprv(l,k))
     &		*(upro(l+1,k)-uprv(l,k)))+
     &       (1.-cnpar)*abs((upro(l+1,k)-upro(l,k))
     &		*(uprv(l+1,k)-upro(l,k))))
     &       /dh/h(l+1)
     &                )

            dv = 0.5*(
     &       (cnpar*abs((vprv(l+1,k)-vprv(l,k))
     &		*(vprv(l+1,k)-vpro(l,k)))+
     &       (1.-cnpar)*abs((vpro(l+1,k)-vpro(l,k))
     &		*(vpro(l+1,k)-vprv(l,k))))
     &       /dh/h(l)
     &      +(cnpar*abs((vprv(l+1,k)-vprv(l,k))
     &		*(vpro(l+1,k)-vprv(l,k)))+
     &       (1.-cnpar)*abs((vpro(l+1,k)-vpro(l,k))
     &		*(vprv(l+1,k)-vpro(l,k))))
     &       /dh/h(l+1)
     &                )

            !m2 = du**2 + dv**2
            m2 = du + dv
            shearf2(l,k) = m2
          end do
        end do

	nfreq = sqrt(n2max)
	nperiod = 0.
	if( nfreq .gt. 0. ) nperiod = 1. / nfreq
	if( iuinfo .le. 0 ) call getinfo(iuinfo)
	write(iuinfo,*) 'n2max: ',it,n2max,nfreq,nperiod

	end

c**************************************************************

      subroutine richardson(k,nldim,rif,rig,ldebug,testnode,icycle)

	use mod_fluidmud
	use mod_gotm_aux
	use mod_ts
	use mod_hydro_print
	use basin, only : nkn,nel,ngr,mbw

        implicit none
 
        integer nldim

        include 'param.h'

	include 'pkonst.h'

        real :: rif(nlvdim), rig(nlvdim)
        real*8 :: dis(nlvdim)





        logical ldebug
        integer k,l,nlev,iwrite,icycle,testnode
        real aux,dh,du,dv,m2,g_dot
        real h(nlvdim), rho1, rho2
        save iwrite
        data iwrite/0/

        if( nldim .ne. nlvdim ) stop 'error stop setbuoyf: dimension'

        call dep3dnod(k,+1,nlev,h)
 
        !write(*,*) sum(eps(:,k))

        do l=1,nlev
          if (buoyf2(l,k) .gt. 0.0000001) then
            rig(l) = shearf2(l,k)/buoyf2(l,k)
            rif(l) = 1./(1.-eps_gotm(l,k)/buoyf2(l,k))
          endif
          if (ldebug .and. k == testnode .and.
     &                mod(iwrite,icycle) .eq. 0) then
            write(1120,'(A10,I5,10F20.12)') 'richardson',
     &      l, shearf2(l,k), buoyf2(l,k), eps_gotm(l,k), 
     &      tken_gotm(l,k), rig(l), rif(l) 
          endif
        end do

        if (k==testnode) iwrite = iwrite + 1
      end subroutine
c**************************************************************

      subroutine stress3d(nldim,k,tstress,ldebug,testnode,icycle)

! m * g = kg * m/s2 --- force 
! kg * m/s2 / m2 = kg / (m * s2) --- stress
! du/dz * nu * rho = 1/s * m2/s * kg/m3 = kg / (m * s2)
!
	use mod_fluidmud
	use mod_gotm_aux
	use mod_ts
	use mod_diff_visc_fric
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer nldim

        include 'param.h'

	include 'pkonst.h'

        real tstress(nlvdim)


        logical ldebug

        integer k,l,nlev,testnode,icycle,iwrite
        real aux,dh,du,dv,m2,dbuoy
        real h(nlvdim), rho1, rho2, visk_bar, sr1, sr2
        real cnpar,stress_x,stress_y, rhobar, g_dot
        save iwrite
        data iwrite/0/

        if( nldim .ne. nlvdim ) stop 'error stop setbuoyf: dimension'
        call dep3dnod(k,+1,nlev,h)
        do l=1,nlev
          g_dot = sqrt(shearf2(l,k))
          rhobar = rhov(l,k)+1000.
          visk_bar = visv(l,k) + vts(l,k)
          tstress(l) = g_dot * visk_bar * rhobar
          if (ldebug .and. k == testnode .and.
     &                  mod(iwrite,icycle) .eq. 0) then
            write(1119,'(I10,12F16.8)')l,g_dot,
     &      vts(l,k),visv(l,k),
     &      rhobar, tstress(l)
          endif
        end do
        if (k == testnode) iwrite = iwrite + 1
        tstress(nlev) = tstress(nlev-1)
      end subroutine stress3d

c**************************************************************

	subroutine bnstress(czdef,taub,areaac)

c computes bottom stress at nodes
c
c this is evaluated for every element and then averaged for each node
c taub (stress at bottom) is accumulated and weighted by area
 
	use mod_hydro_vel
	use evgeom
	use levels
	use basin

	implicit none

	real czdef
	real taub(1)
	real areaac(1)

	include 'param.h'



	integer k,ie,ii,n,nlev,imud
	real aj,taubot,getpar

	real areaele

	imud = getpar('imud')

c	---------------------------------------------------
c	initialize arrays
c	---------------------------------------------------

        do k=1,nkn
          taub(k) = 0.
          areaac(k) = 0.
        end do
 
c	---------------------------------------------------
c	accumulate
c	---------------------------------------------------
!AR: mud 
        if (imud == 0) then
          do ie=1,nel
 
            !call elebase(ie,n,ibase)
	    n = 3
            aj = ev(10,ie)
	    nlev = ilhv(ie)

            taubot = czdef * ( ulnv(nlev,ie)**2 + vlnv(nlev,ie)**2 )
            do ii=1,n
              k = nen3v(ii,ie)
              taub(k) = taub(k) + taubot * aj
              areaac(k) = areaac(k) + aj
            end do
          end do
        else
          do ie=1,nel
            !call elebase(ie,n,ibase)
	    n = 3
            aj = ev(10,ie)
            nlev = ilhv(ie)
            do ii=1,n
              k = nen3v(ii,ie)
!AR: take care with nlev assumed to be constant for the nodes
              call set_bottom_stress(k,taubot)
              taub(k) = taub(k) + taubot * aj
              areaac(k) = areaac(k) + aj
            end do
          end do
        endif

c	---------------------------------------------------
c	compute bottom stress
c	---------------------------------------------------

        do k=1,nkn
          if( areaac(k) .le. 0. ) stop 'error stop bnstress: (2)'
          taub(k) = taub(k) / areaac(k)
        end do

	end

c**************************************************************

	subroutine notused
 
c this is evaluated for every element and then averaged for each node
c the average is weigthed with the volume of each element
c level l (for node) refers to interface between layers l and l+1
c
c taub (stress at bottom) is also accumulated and weighted by area

	!implicit none

	use mod_ts
	use mod_hydro_vel
	use levels
	use basin

	include 'param.h'
 

 
        real shearf2(nlvdim,nkndim)
        real buoyf2(nlvdim,nkndim)
        real volf2(nlvdim,nkndim)
        real taub(nkndim)
        real areaac(nkndim)
 

	real h(nlvdim)

	  czdef = getpar('czdef')

        do k=1,nkn
	  nlev = ilhkv(k)
	  do l=1,nlev
	    shearf2(l,k) = 0.
	    volf2(l,k) = 0.
	  end do
          taub(k) = 0.
          areaac(k) = 0.
        end do
 
        do ie=1,nel
 
          call dep3dele(ie,+1,nlev,h)
          !call elebase(ie,n,ibase)
	  n = 3
          area = areaele(ie)
          arean = area/n
 
          do l=1,nlev-1
            dh = 0.5 * (h(l)+h(l+1))
            du = ulnv(l,ie) - ulnv(l+1,ie)
            dv = vlnv(l,ie) - vlnv(l+1,ie)
            m2 = (du**2 + dv**2) / dh**2
            vol = dh * arean
            do ii=1,n
              k = nen3v(ii,ie)
              shearf2(l,k) = shearf2(l,k) + m2 * vol
              volf2(l,k) = volf2(l,k) + vol
            end do
          end do
 
          taubot = czdef * ( ulnv(nlev,ie)**2 + vlnv(nlev,ie)**2 )
          do ii=1,n
            k = nen3v(ii,ie)
            taub(k) = taub(k) + taubot * arean
            areaac(k) = areaac(k) + arean
          end do
 
        end do
 
        do k=1,nkn
          nlev = ilhkv(k)
          do l=1,nlev-1
            if( volf2(l,k) .le. 0. ) stop 'error stop ... (1)'
            shearf2(l,k) = shearf2(l,k) / volf2(l,k)
          end do
          if( areaac(k) .le. 0. ) stop 'error stop ... (2)'
          taub(k) = taub(k) / areaac(k)
        end do

	end

c**************************************************************

	subroutine checka(nldim,shearf2,buoyf2,taub)

c checks arrays for nan or other strange values

	use mod_meteo
	use mod_diff_visc_fric
	use mod_hydro_print
	use mod_hydro_vel
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nldim
	real buoyf2(nldim,1)
	real shearf2(nldim,1)
	real taub(1)

	include 'param.h'




	if( nlvdim .ne. nldim ) stop 'error stop checka'

	call nantest(nkn*nlvdim,shearf2,'shearf2')
	call nantest(nkn*nlvdim,buoyf2,'buoyf2')
	call nantest(nkn,taub,'taub')

	call nantest(nkn*(nlvdim+1),visv,'visv')
	call nantest(nkn*(nlvdim+1),difv,'difv')

	call nantest(nkn*nlvdim,uprv,'uprv')
	call nantest(nkn*nlvdim,vprv,'vprv')
	call nantest(nel*nlvdim,ulnv,'ulnv')
	call nantest(nel*nlvdim,vlnv,'vlnv')

	call nantest(nkn,tauxnv,'tauxnv')
	call nantest(nkn,tauynv,'tauynv')
 
	end

c**************************************************************

	subroutine checkb(text)

c checks arrays for strange values

	use mod_meteo
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) text

	include 'param.h'





	integer one,three
	real zero,valmax
	real amin,amax

	one = 1
	three = 3
	zero = 0.

	valmax = 5.
	call valtest(one,nkn,-valmax,valmax,znv,text,'znv')
	call valtest(three,nel,-valmax,valmax,zenv,text,'zenv')

	valmax = 10000.
	call valtest(nlvdim+1,nkn,zero,valmax,visv,text,'visv')
	call valtest(nlvdim+1,nkn,zero,valmax,difv,text,'difv')

	valmax = 100.
	call valtest(nlvdim,nkn,-valmax,valmax,rhov,text,'rhov')

	valmax = 3.
	call valtest(nlvdim,nel,-valmax,valmax,ulnv,text,'ulnv')
	call valtest(nlvdim,nel,-valmax,valmax,vlnv,text,'vlnv')

	valmax = 100.
	call valtest(nlvdim,nel,-valmax,valmax,utlnv,text,'utlnv')
	call valtest(nlvdim,nel,-valmax,valmax,vtlnv,text,'vtlnv')

	end

c**************************************************************

	subroutine keps_shell

	use mod_turbulence
	use mod_depth
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'









	integer k,lmax,l
	real rho0,rhoair
	real cd,cb
	real wx,wy
	real ubot,vbot
	real dtreal,dt
	real taus,taub

	integer icall
	save icall
	data icall /0/

	if( icall .eq. 0 ) then
	  call keps_init
	  icall = 1
	end if

	call get_timestep(dtreal)
	dt = dtreal
	rho0 = 1024.
	rhoair = 1.225
	cd = 2.5e-3
	cb = 2.5e-3

	do k=1,nkn

	  lmax = ilhkv(k)
	  ubot = uprv(lmax,k)
	  vbot = vprv(lmax,k)
	  call get_wind(k,wx,wy)
	  taus = rhoair * cd * (wx*wx+wy*wy)
	  taub = rho0 * cb * (ubot*ubot+vbot*vbot)

          call keps(lmax,dt,rho0
     +		,taus,taub
     +		,hdknv(1,k),uprv(1,k),vprv(1,k)
     +		,rhov(1,k),visv(0,k),difv(0,k),tken(0,k),eps(0,k))

	end do

	end

c**************************************************************

        subroutine keps_init

c initializes arrays for keps routine

	use mod_turbulence
	use mod_diff_visc_fric
	use basin, only : nkn,nel,ngr,mbw

        implicit none

	real kmin,epsmin,lenmin,avumol,avtmol,avsmol
c       parameter(kmin=1.e-10,epsmin=1.e-12,lenmin=0.01)
        parameter(kmin=3.e-6,epsmin=5.e-10,lenmin=0.01)
        parameter(avumol=1.3e-6,avtmol=1.4e-7,avsmol=1.1e-9)

	include 'param.h'






        integer k,l

	do k=1,nkn
          do l=0,nlvdim
            tken(l,k) = kmin
            eps(l,k) = epsmin
            rls(l,k) = lenmin
            visv(l,k) = avumol
            difv(l,k) = avsmol
	  end do
        end do

        end

c**************************************************************
c**************************************************************
c debug routines
c**************************************************************
c**************************************************************

	subroutine write_node_vel_info(iunit,it,k,ndim,depth
     +			,u_taus,u_taub,z0s,z0b
     +			,num,nuh,ken,dis,len,nn,ss,hh)

	implicit none

	integer iunit,it,k,ndim
	double precision depth,u_taus,u_taub,z0s,z0b

	double precision hh(0:ndim)
	double precision nn(0:ndim), ss(0:ndim)

	double precision num(0:ndim), nuh(0:ndim)
	double precision ken(0:ndim), dis(0:ndim)
	double precision len(0:ndim)

	integer l

	write(iunit,*) '----------------------------------'
	write(iunit,*) 'time,node,nlev: ',it,k,ndim
	write(iunit,*) 'depth: ',depth
	write(iunit,*) 'u_taus,u_taub: ',u_taus,u_taub
	write(iunit,*) 'z0s,z0b: ',z0s,z0b
	write(iunit,*) 'num:',(num(l),l=0,ndim)
	write(iunit,*) 'nuh:',(nuh(l),l=0,ndim)
	write(iunit,*) 'ken:',(ken(l),l=0,ndim)
	write(iunit,*) 'dis:',(dis(l),l=0,ndim)
	write(iunit,*) 'len:',(len(l),l=0,ndim)
	write(iunit,*) 'nn :',(nn(l),l=0,ndim)
	write(iunit,*) 'ss :',(ss(l),l=0,ndim)
	write(iunit,*) 'hh :',(hh(l),l=0,ndim)

	end

c**************************************************************

	subroutine write_node_bin_info(iunit,it,k,ndim,depth
     +			,u_taus,u_taub,z0s,z0b
     +			,num,nuh,ken,dis,len,nn,ss,hh)

	implicit none

	integer iunit,it,k,ndim
	double precision depth,u_taus,u_taub,z0s,z0b

	double precision hh(0:ndim)
	double precision nn(0:ndim), ss(0:ndim)

	double precision num(0:ndim), nuh(0:ndim)
	double precision ken(0:ndim), dis(0:ndim)
	double precision len(0:ndim)

	integer l

	write(iunit) it,k,ndim
	write(iunit) depth
	write(iunit) u_taus,u_taub
	write(iunit) z0s,z0b
	write(iunit) (num(l),l=0,ndim)
	write(iunit) (nuh(l),l=0,ndim)
	write(iunit) (ken(l),l=0,ndim)
	write(iunit) (dis(l),l=0,ndim)
	write(iunit) (len(l),l=0,ndim)
	write(iunit) (nn(l),l=0,ndim)
	write(iunit) (ss(l),l=0,ndim)
	write(iunit) (hh(l),l=0,ndim)

	end

c**************************************************************

