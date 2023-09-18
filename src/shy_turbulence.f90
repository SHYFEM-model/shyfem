
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

! gotm module
!
! contents :
!
! revision log :
!
! 10.08.2003	ggu	call gotm_init
! 05.10.2004	ggu	new administration routine turb_closure, Munk-And.
! 23.03.2006	ggu	changed time step to double precision
! 20.11.2006	ggu	new version of keps, for gotm changed common blocks
! 20.10.2007	ccf	new version of gotm (4.0.0)
! 10.04.2008	ggu	integrated in main branch
! 18.09.2008	ccf	bug fix for m2 in setm2n2
! 02.12.2008	ggu	bug in gotm_init: no limiting values for initialization
! 18.12.2008	ggu	bug in GOTM module and setm2n2() corrected
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2011	ggu	write n2max to info file, profiles in special node
! 30.03.2012	ggu	changed VERS_6_1_51
! 29.03.2013	ggu	avoid call to areaele -> ev(10,ie)
! 13.06.2013	ggu	changed VERS_6_1_65
! 25.03.2014	ggu	new offline
! 18.06.2014	ggu	changed VERS_6_1_77
! 27.06.2014	ggu	changed VERS_6_1_78
! 05.11.2014	ggu	changed VERS_7_0_5
! 05.12.2014	ggu	changed VERS_7_0_8
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_2
! 19.01.2015	ggu	changed VERS_7_1_3
! 10.07.2015	ggu	changed VERS_7_1_50
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 30.07.2015	ggu	changed VERS_7_1_83
! 03.12.2015	ccf	levdbg introduced for checka
! 16.12.2015	ggu	changed VERS_7_3_16
! 18.12.2015	ggu	changed VERS_7_3_17
! 25.05.2016	ggu	changed VERS_7_5_10
! 05.12.2017	ggu	changed VERS_7_5_39
! 03.02.2019	ggu	in gotm_shell check for 0layer and z0s/bmin (GGUZ0)
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
!
!**************************************************************
!--------------------------------------------------------------------
        module shy_turbulence
!--------------------------------------------------------------------
        contains
!--------------------------------------------------------------------

        subroutine turb_closure

! administers turbulence closure

        use para
        use offline_data

        implicit none

        logical boff

        integer iturb
        save iturb
        data iturb / 0 /

        if( iturb .lt. 0 ) return

        call is_offline(4,boff)
        if( boff ) return

        if( iturb .eq. 0 ) then
          iturb = nint(getpar('iturb'))
          if( iturb .le. 0 ) iturb = -1
          if( iturb .lt. 0 ) return
          write(*,*) 'starting turbulence model: iturb = ',iturb
        end if

        if( iturb .eq. 1 ) then         !Gotm
          call gotm_shell
        else if( iturb .eq. 2 ) then    !keps
          call keps_shell
        else if( iturb .eq. 3 ) then    !Munk Anderson
          call munk_anderson_shell
        else
          write(6,*) 'Value iturb not possible: ',iturb
          stop 'error stop turb_closure: iturb'
        end if

        end

!**************************************************************

	subroutine munk_anderson_shell

! computes turbulent quantities with Munk - Anderson model

	use diffusion
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
        use para
        use elems_dealing

	implicit none

	include 'femtime.h'
	include 'pkonst.h'

!---------------------------------------------------------------
! aux arrays superposed onto other aux arrays
!---------------------------------------------------------------

	double precision shearf2(nlvdi,nkn)
	double precision buoyf2(nlvdi,nkn)
	double precision richard(nlvdi,nkn)
	double precision h(nlvdi)


	integer k,l
	integer nlev
	integer mode
	double precision ri,vis,dif
	double precision diftur,vistur
	double precision a,b,alpha,beta

	integer icall
	save icall
	data icall / 0 /

!------------------------------------------------------
! initialization
!------------------------------------------------------

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  write(*,*) 'starting Munk Anderson turbulence model'
	end if

	icall = icall + 1

!------------------------------------------------------
! set up buoyancy frequency and shear frequency
!------------------------------------------------------

	call setm2n2(nlvdi,buoyf2,shearf2)

!------------------------------------------------------
! set up parameters
!------------------------------------------------------

	mode = +1
	a = 10.
	b = 3.33
	alpha = -1./2.
	beta = -3./2.
        vistur = getpar('vistur')
        diftur = getpar('diftur')

!------------------------------------------------------
! compute richardson number for each water column
!------------------------------------------------------

	do k=1,nkn

	    nlev = nlvdi
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

!------------------------------------------------------
! end of routine
!------------------------------------------------------

	end

!**************************************************************

        subroutine gotm_shell

! computes turbulent quantities with GOTM model

        use meteo
        use gotm_aux
        use ts
        use roughness
        use diffusion
        use hydro_print
        use levels, only : nlvdi,nlv
        use basin
        use shympi
        use para
        use elems_dealing
        use restart
        use waves_admin
        use time_util
        use turbulence_util

        implicit none

        include 'femtime.h'
        include 'pkonst.h'

        double precision dt
        double precision u_taus,u_taub

        double precision hh(0:nlvdi)
        double precision nn(0:nlvdi), ss(0:nlvdi)

        double precision num(0:nlvdi), nuh(0:nlvdi)
        double precision ken(0:nlvdi), dis(0:nlvdi)
        double precision len(0:nlvdi)

        double precision num_old(0:nlvdi), nuh_old(0:nlvdi)
        double precision ken_old(0:nlvdi), dis_old(0:nlvdi)
        double precision len_old(0:nlvdi)
 
!---------------------------------------------------------------
! aux arrays superposed onto other aux arrays
!---------------------------------------------------------------

#ifdef DEBUGON
        double precision taub(nkn_local)
#else
        double precision taub(nkn)
#endif

        integer ioutfreq,ks
        integer iunit,id
        integer k,l
        integer laux
        integer nlev
        double precision g
        double precision czdef
        save czdef

        double precision h(nlvdi)
        double precision depth          !total depth [m]
        double precision z0s,z0b        !surface/bottom roughness length [m]
        double precision rlmax,dz0
        integer nltot
        logical bwrite

        double precision, parameter :: dz0min = 1.1     !min value for dz0=d/z0
        double precision, parameter :: charnock_val=1400.           !emp. Charnock constant
        double precision ubot,vbot,rr

        double precision dtreal

        logical bwave,bgotm
        save bwave

        character*80 fn
        integer icall
        save icall
        data icall / 0 /

        integer, save   :: levdbg

!------------------------------------------------------
! documentation
!------------------------------------------------------

! total stress :  tau_x = a |u| u   tau_y = a |u| v
! tau_tot^2 = a^2 |u|^2 ( u^2 + v^2 ) a^2 |u|^4
! tau_tot = sqrt(tau_x^2+tau_y^2) = a |u|^2
!
! u_taus=sqrt(sqrt(tx*tx+ty*ty))	(friction velocity)

!------------------------------------------------------
! initialization
!------------------------------------------------------

        if( icall .lt. 0 ) return

        if( icall .eq. 0 ) then

          call has_gotm(bgotm)
          if( .not. bgotm ) then
           write(6,*) 'the model has been compiled without GOTM support'
           write(6,*) 'please eneable GOTM=true in the Rules.make file'
           write(6,*) 'or set iturb to another value'
           stop 'error stop gotm_shell: no GOTM support'
          end if

          czdef = getpar('czdef')
          bwave = has_waves()
           
          !! ivb - DON'T DO IN CASE OF RESTART  
          if (.not. rst_use_restart(7)) then  
            visv = 0.d0
            difv = 0.d0
          end if  

          if( nlvdi <= 1 ) then
            icall = -1
            return
          end if

!         --------------------------------------------------------
!         Initializes gotm arrays 
!         --------------------------------------------------------

          write(*,*) 'starting initializing GOTM turbulence model'

          !! ivb - DON'T DO IN CASE OF RESTART  
          if (.not. rst_use_restart(7)) then   
            call gotm_init
          end if  

!         --------------------------------------------------------
!         Get file name containing GOTM turbulence model parameters
!         --------------------------------------------------------

          call getfnm('gotmpa',fn)

          iunit = 10
          call init_gotm_turb(iunit,fn,nlvdi)

          levdbg = nint(getpar('levdbg'))

          write(*,*) 'finished initializing GOTM turbulence model'

          icall = 1
        end if

        call get_timestep(dtreal)
        dt = dtreal
        g = grav

!------------------------------------------------------
! set up bottom stress on nodes
!------------------------------------------------------

        call bnstress(czdef,taub)

        !call shympi_comment('exchanging taub')
        call shympi_exchange_2d_node(taub)
        !call shympi_barrier

!------------------------------------------------------
! set up buoyancy frequency and shear frequency
!------------------------------------------------------

        shearf2 = 0.d0
        buoyf2  = 0.d0
        call setm2n2(nlvdi,buoyf2,shearf2)

!------------------------------------------------------
! call gotm for each water column
!------------------------------------------------------

        rlmax = 0.
        nltot = 0

        do k=1,nkn

            nlev = nlvdi
            call dep3dnod(k,+1,nlev,h)

            if( count( h(1:nlev) <= 0. ) > 0 ) goto 97
            if( nlev .eq. 1 ) goto 1

!           ------------------------------------------------------
!           update boyancy and shear-frequency vectors
!           ------------------------------------------------------

            do l=1,nlev-1
              laux = nlev - l
              nn(laux) = buoyf2(l,k)
              ss(laux) = shearf2(l,k)
            end do
            nn(0) = 0.
            nn(nlev) = 0.
            ss(0) = ss(1)
            ss(nlev) = ss(nlev-1)

!           ------------------------------------------------------
!           compute layer thickness and total depth
!           ------------------------------------------------------

            depth = 0.d0
            do l=1,nlev
              hh(nlev-l+1) = h(l)
              depth = depth + h(l)
            end do

!           ------------------------------------------------------
!           compute surface friction velocity (m/s)
!           ------------------------------------------------------

            u_taus = sqrt( sqrt( tauxnv(k)**2 + tauynv(k)**2 ) )

            if ( bwave ) then
              z0s = z0sk(k)
            else
              z0s = charnock_val*u_taus**2/g
            end if
            z0s = max(z0s,z0smin)                       !GGUZ0

!           ------------------------------------------------------
!           compute bottom friction velocity (m/s)
!           ------------------------------------------------------

            z0b = z0bk(k)
            z0b = max(z0bmin,z0b)                       !GGUZ0
            u_taub = sqrt( taub(k) )

            ubot = uprv(nlev,k)
            vbot = vprv(nlev,k)
            dz0 = hh(1)/z0b
            if( dz0 < dz0min ) dz0 = dz0min             !GGUZ0
            rr = 0.4/( log( 0.5*(1.+dz0) ) )
            !rr = 0.4/(log((z0b+hh(1)/2)/z0b))
            u_taub = rr*sqrt( ubot*ubot + vbot*vbot )

!           ------------------------------------------------------
!           update 1-dimensional vectors
!           ------------------------------------------------------

            do l=0,nlev !ivb    !needed for restart
              laux = nlev - l
              numv_gotm(laux,k) = visv(l,k)
              nuhv_gotm(laux,k) = difv(l,k)
              tken_gotm(laux,k) = tken(l,k)
              eps_gotm(laux,k) = eps(l,k)
              rls_gotm(laux,k) = rls(l,k)
            end do

            num(:) = numv_gotm(:,k)
            nuh(:) = nuhv_gotm(:,k)
            ken(:) = tken_gotm(:,k)
            dis(:) = eps_gotm(:,k)
            len(:) = rls_gotm(:,k)
            num_old(:) = numv_gotm(:,k)
            nuh_old(:) = nuhv_gotm(:,k)
            ken_old(:) = tken_gotm(:,k)
            dis_old(:) = eps_gotm(:,k)
            len_old(:) = rls_gotm(:,k)

	    !call save_gotm_init

!           ------------------------------------------------------
!           call GOTM turbulence routine
!           ------------------------------------------------------

            
            bwrite = ipv(k) .eq. 44041
            bwrite = .false.
            if (bwrite) then !ivb
              write(6,*) '1==========================='
              write(6,*) '1time___: ',it
              write(6,*) '1node___: ',k
              write(6,*) '1nlev___: ',nlev
              write(6,*) '1dt_____: ',dt
              write(6,*) '1depth__: ',depth
              write(6,*) '1utaus_b: ',u_taus,u_taub
              write(6,*) '1z0s_z0b: ',z0s,z0b
              write(6,*) '1hh_____: ',(hh     (l),l=0,nlev)
              write(6,*) '1nn_____: ',(nn     (l),l=0,nlev)
              write(6,*) '1ss_____: ',(ss     (l),l=0,nlev)
              write(6,*) '1num____: ',(num    (l),l=0,nlev)
              write(6,*) '1nuh____: ',(nuh    (l),l=0,nlev)
              write(6,*) '1ken____: ',(ken    (l),l=0,nlev)
              write(6,*) '1dis____: ',(dis    (l),l=0,nlev)
              write(6,*) '1len____: ',(len    (l),l=0,nlev)
              write(6,*) '1==========================='
            end if

            call do_gotm_turb (nlev,dt,depth,u_taus,u_taub,z0s,z0b,hh,nn,ss,num,nuh,ken,dis,len)

            bwrite = ipv(k) .eq. 44041
            bwrite = .false.
            if (bwrite) then !ivb
              write(6,*) '2==========================='
              write(6,*) '2time___: ',it
              write(6,*) '2node___: ',k
              write(6,*) '2nlev___: ',nlev
              write(6,*) '2dt_____: ',dt
              write(6,*) '2depth__: ',depth
              write(6,*) '2utaus_b: ',u_taus,u_taub
              write(6,*) '2z0s_z0b: ',z0s,z0b
              write(6,*) '2hh_____: ',(hh     (l),l=0,nlev)
              write(6,*) '2nn_____: ',(nn     (l),l=0,nlev)
              write(6,*) '2ss_____: ',(ss     (l),l=0,nlev)
              write(6,*) '2num____: ',(num    (l),l=0,nlev)
              write(6,*) '2nuh____: ',(nuh    (l),l=0,nlev)
              write(6,*) '2ken____: ',(ken    (l),l=0,nlev)
              write(6,*) '2dis____: ',(dis    (l),l=0,nlev)
              write(6,*) '2len____: ',(len    (l),l=0,nlev)
              write(6,*) '2==========================='
            end if
!           ------------------------------------------------------
!           copy back to node vectors
!           ------------------------------------------------------

            numv_gotm(:,k) = num(:)
            nuhv_gotm(:,k) = nuh(:)
            tken_gotm(:,k) = ken(:)
            eps_gotm(:,k)  = dis(:)
            rls_gotm(:,k)  = len(:)

            bwrite = .false.
            do l=0,nlev
              rlmax = max(rlmax,len(l))	!ggu
              if( len(l) .gt. 100. ) then
                nltot = nltot + 1
                bwrite = .true.
              end if

            end do

            bwrite = k == 34
            bwrite = .false.
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

!           ------------------------------------------------------
!           update viscosity and diffusivity
!           ------------------------------------------------------

            do l=0,nlev !ivb
              laux = nlev - l
              visv(l,k) = num(laux)
              difv(l,k) = nuh(laux)
              tken(l,k) = ken(laux)     !ivb needed for restart              
              eps(l,k) = dis(laux)      !ivb needed for restart
              rls(l,k) = len(laux)      !ivb needed for restart
            end do

    1     continue
        end do

        if( levdbg >= 2 ) call checka(nlvdi,shearf2,buoyf2,taub)
 
	ks = 0			!internal node number
	ioutfreq = 3600		!output frequency
	if( ks .gt. 0 .and. mod(it,ioutfreq) .eq. 0 ) then
	  write(188,*) it,nlev,(visv(l,ks),l=1,nlev)
	  write(189,*) it,nlev,(difv(l,ks),l=1,nlev)
	end if

!------------------------------------------------------
! end of routine
!------------------------------------------------------

	return
   97	continue
	write(6,*) 'layers without depth...'
	write(6,*) k,nlev
	write(6,*) h(1:nlev)
	stop 'error stop gotm_shell: no layer'
	end

!**************************************************************

        subroutine gotm_init

! initializes gotm arrays

        use gotm_aux

        implicit none

        double precision, parameter    :: num_min  = 1.e-6
        double precision, parameter    :: nuh_min  = 1.e-6
        double precision, parameter    :: tken_min = 1.e-10
        double precision, parameter    :: eps_min  = 1.e-12
        double precision, parameter    :: rls_min  = 1.e-10

        numv_gotm = num_min
        nuhv_gotm = nuh_min
        tken_gotm = tken_min
        eps_gotm  = eps_min
        rls_gotm  = rls_min

        end

!**************************************************************

	subroutine gotm_get(k,nlev,num,nuh,tk,ep,rl)

! returns internal parameters from turbulence closure

	use gotm_aux

	implicit none

	integer k		!node
	integer nlev		!number of levels to return
	double precision num(0:nlev)	!viscosity
	double precision nuh(0:nlev)	!diffusivity
	double precision tk (0:nlev)	!kinetic energy
	double precision ep (0:nlev)	!dissipation
	double precision rl (0:nlev)	!length scale

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

!**************************************************************

      subroutine Yevol(Nmx,dt,h,nuh,difmol,Qsour,Yvar,Yd,Bcup,Yup,Bcdw,Ydw,Taur)

      implicit none

	integer Nmx		!number of levels
	double precision dt			!time step [s]
	double precision h(Nmx)		!level thickness [m]
	double precision nuh(0:Nmx)		!turbulent diffusivity [m**2/s]
	double precision difmol		!molecular diffusivity [m**2/s]
	double precision Qsour(Nmx)		!internal source term
	double precision Yvar(Nmx)		!variable to diffuse
	double precision Yd(Nmx)		!relaxation value for variable
	integer Bcup		!type of boundary condition (upper)
	double precision Yup		!value for boundary condition (upper)
	integer Bcdw		!type of boundary condition (lower)
	double precision Ydw		!value for boundary condition (lower)
	double precision Taur		!time scale to use with Yd [s]

! if Taur <= 0	-> do not use Yd
! Bcup/Bcdw :   1 = Neuman condition   2 = Dirichlet condition
! Yup,Ydw   :   value to use for boundary condition
!
! example: Bcup = 1  and Yup = 0  --> no flux condition across surface

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
          du(i)=Y(i)+dt*Qsour(i)+(1-cnpar)*(a*Y(i-1)-(a+c)*Y(i)+c*Y(i+1)) !i,n
         end do

! Surface 
! => [a b /].X=[d] where X=[Nmx-1,Nmx,/]

         if (Bcup.eq.1) then                       !BC Neuman
          a      =2*dt*avh(Nmx-1)/(h(Nmx)+h(Nmx-1))/h(Nmx)
          au(Nmx)=-cnpar*a
          bu(Nmx)=1-au(Nmx)
          du(Nmx)=Y(Nmx)+dt*(Qsour(Nmx)-Yup/h(Nmx))+(1-cnpar)*a*(Y(Nmx-1)-Y(Nmx))
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
          du(1)=Y(1)+dt*(Qsour(1)+Ydw/h(1))+(1-cnpar)*c*(Y(2)-Y(1))
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

!**************************************************************

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

!**************************************************************

	subroutine setm2n2(nldim,buoyf2,shearf2)

! sets buoyancy and shear frequency
!
! this is done directly for each node
! in rhov is already rho^prime = rho - rho0 (deviation)
!  Discretisation of vertical shear squared according to Burchard (2002)
!  in order to guarantee conservation of kinetic energy when transformed
!  from mean kinetic energy to turbulent kinetic energy.
!
! bug fix in computation of shearf2 -> abs() statements to avoid negative vals

	use ts
	use hydro_print
	use basin, only : nkn,nel,ngr,mbw
	use shympi
        use elems_dealing
        use defnames

	implicit none

	integer nldim
	double precision buoyf2(nldim,nkn)
	double precision shearf2(nldim,nkn)

	include 'pkonst.h'
	include 'femtime.h'

	integer k,l,nlev
	double precision aux,dh,du,dv,m2,dbuoy
	double precision h(nldim)
	double precision cnpar			!numerical "implicitness" parameter
	double precision n2max,n2
	double precision nfreq,nperiod

	integer iuinfo
	save iuinfo
	data iuinfo / 0 /
 
        aux = -grav / rowass
	cnpar = 1
	n2max = 0.
 
        do k=1,nkn
	  nlev = nldim
          call dep3dnod(k,+1,nlev,h)
          do l=1,nlev-1
            dh = 0.5 * ( h(l) + h(l+1) )
            dbuoy = aux * ( rhov(l,k) - rhov(l+1,k) )
            n2 = dbuoy / dh
	    n2max = max(n2max,n2)
            buoyf2(l,k) = n2

            du = 0.5*((cnpar*abs((uprv(l+1,k)-uprv(l,k))*(uprv(l+1,k)-upro(l,k)))+      &
     &       (1.-cnpar)*abs((upro(l+1,k)-upro(l,k))*(upro(l+1,k)-uprv(l,k))))/dh/h(l)   &
     &      +(cnpar*abs((uprv(l+1,k)-uprv(l,k))*(upro(l+1,k)-uprv(l,k)))+               &
     &       (1.-cnpar)*abs((upro(l+1,k)-upro(l,k))*(uprv(l+1,k)-upro(l,k))))           &
     &       /dh/h(l+1))

            dv = 0.5*((cnpar*abs((vprv(l+1,k)-vprv(l,k))*(vprv(l+1,k)-vpro(l,k)))+      &
     &       (1.-cnpar)*abs((vpro(l+1,k)-vpro(l,k))*(vpro(l+1,k)-vprv(l,k))))/dh/h(l)   &
     &      +(cnpar*abs((vprv(l+1,k)-vprv(l,k))*(vpro(l+1,k)-vprv(l,k)))+               &
     &       (1.-cnpar)*abs((vpro(l+1,k)-vpro(l,k))*(vprv(l+1,k)-vpro(l,k))))           &
     &       /dh/h(l+1))

            !m2 = du**2 + dv**2
            m2 = du + dv
            shearf2(l,k) = m2
          end do
        end do

        n2max = shympi_max(n2max)

	nfreq = sqrt(n2max)
	nperiod = 0.
	if( nfreq .gt. 0. ) nperiod = 1. / nfreq
	if( iuinfo .le. 0 ) call getinfo(iuinfo)
	if(shympi_is_master()) then
	  write(iuinfo,*) 'n2max: ',it,n2max,nfreq,nperiod
	end if

	end

!**************************************************************

        subroutine bnstress(czdef,taub)

! computes bottom stress at nodes
!
! this is evaluated for every element and then averaged for each node
! taub (stress at bottom) is accumulated and weighted by area
 
        use hydro_vel
        use evgeom
        use levels
        use basin
        use shympi
        use area

        implicit none

        double precision czdef
#ifdef DEBUGON
        integer e
        double precision taub(nkn_local)
#else
        double precision taub(nkn)
#endif

        integer k,ie,ii,n,nlev
        double precision aj,taubot

!	---------------------------------------------------
!	initialize arrays
!	---------------------------------------------------

        do k=1,nkn
          taub(k) = 0.d0
        end do
 
!	---------------------------------------------------
!	accumulate
!	---------------------------------------------------
#ifdef DEBUGON
        call shympi_exchange_halo_2d_elems(ilhv)
        call shympi_exchange_halo_3d_elems(ulnv)
        call shympi_exchange_halo_3d_elems(vlnv)

        do e=1,nel_local
            if(bmpi) then
               ie=domain%elems%mapID(e)
            else
               ie=e
            end if
#else
        do ie=1,nel
#endif
          !call elebase(ie,n,ibase)
          n = 3
          aj = 4. * ev(10,ie)
          nlev = ilhv(ie)

          taubot = czdef * ( ulnv(nlev,ie)**2 + vlnv(nlev,ie)**2 )
          do ii=1,n
            k = nen3v(ii,ie)
            taub(k) = taub(k) + taubot * aj
          end do

        end do

!       shympi_elem: exchange taub
!       for area as weight use surface value areakv(1,k)

#ifndef DEBUGON
        !call shympi_comment('shympi_elem: exchange taub')
        call shympi_exchange_and_sum_2d_nodes(taub)
#endif

!	---------------------------------------------------
!	compute bottom stress
!	---------------------------------------------------

	if( any( areakv(1,:) <= 0. ) ) stop 'error stop bnstress: (2)'
        taub = taub / areakv(1,:)

	end

!**************************************************************

	subroutine checka(nldim,shearf2,buoyf2,taub)

! checks arrays for nan or other strange values

	use meteo
	use diffusion
	use hydro_print
	use hydro_vel
        use chk_NaN
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nldim
	double precision buoyf2(nldim,nkn)
	double precision shearf2(nldim,nkn)
	double precision taub(nkn)

	call nantest(nkn*nldim,shearf2,'shearf2')
	call nantest(nkn*nldim,buoyf2,'buoyf2')
	call nantest(nkn,taub,'taub')

	call nantest(nkn*(nldim+1),visv,'visv')
	call nantest(nkn*(nldim+1),difv,'difv')

	call nantest(nkn*nldim,uprv,'uprv')
	call nantest(nkn*nldim,vprv,'vprv')
	call nantest(nel*nldim,ulnv,'ulnv')
	call nantest(nel*nldim,vlnv,'vlnv')

	call nantest(nkn,tauxnv,'tauxnv')
	call nantest(nkn,tauynv,'tauynv')
 
	end

!**************************************************************

	subroutine checkb(text)

! checks arrays for strange values

	use meteo
	use ts
	use diffusion
	use hydro_print
	use hydro_vel
	use hydro_admin
        use chk_NaN
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) text

	integer one,three
	double precision zero,valmax
	double precision amin,amax

	one = 1
	three = 3
	zero = 0.

	valmax = 5.
	call valtest(one,nkn,-valmax,valmax,znv,text,'znv')
	call valtest(three,nel,-valmax,valmax,zenv,text,'zenv')

	valmax = 10000.
	call valtest(nlvdi+1,nkn,zero,valmax,visv,text,'visv')
	call valtest(nlvdi+1,nkn,zero,valmax,difv,text,'difv')

	valmax = 100.
	call valtest(nlvdi,nkn,-valmax,valmax,rhov,text,'rhov')

	valmax = 3.
	call valtest(nlvdi,nel,-valmax,valmax,ulnv,text,'ulnv')
	call valtest(nlvdi,nel,-valmax,valmax,vlnv,text,'vlnv')

	valmax = 100.
	call valtest(nlvdi,nel,-valmax,valmax,utlnv,text,'utlnv')
	call valtest(nlvdi,nel,-valmax,valmax,vtlnv,text,'vtlnv')

	end

!**************************************************************

        subroutine keps_shell

        use basin, only : nkn,nel,ngr,mbw
        use turbulence_util
        use layer_thickness
        use ts
        use diffusion
        use hydro_print
        use levels
        use meteo_forcing
        use restart
        use kepsilon
        use time_util

        implicit none

        include 'femtime.h'

        integer k,lmax,l
        double precision rho0,rhoair
        double precision cd,cb
        double precision wx,wy
        double precision ubot,vbot
        double precision dtreal,dt
        double precision taus,taub

        integer icall
        save icall
        data icall /0/

        if( icall .eq. 0 ) then
          if (.not. rst_use_restart(7)) then  !ivb
            call keps_init
          end if  
          icall = 1
        end if

        call get_timestep(dtreal)
        dt      = dtreal
        rho0    = 1024.d0
        rhoair  = 1.225d0
        cd      = 2.5e-3
        cb      = 2.5e-3

        do k=1,nkn

          lmax = ilhkv(k)
          ubot = uprv(lmax,k)
          vbot = vprv(lmax,k)
          call get_wind(k,wx,wy)
          taus = rhoair * cd * (wx*wx+wy*wy)
          taub = rho0 * cb * (ubot*ubot+vbot*vbot)

          call keps(lmax,dt,rho0,taus,taub                      &
!     +		,hdknv(1,k),uprv(1,k),vprv(1,k)
!     +		,rhov(1,k),visv(0,k),difv(0,k),tken(0,k),eps(0,k))
     &		,hdknv(1:lmax,k),uprv(1:lmax,k),vprv(1:lmax,k)  &  !ivb
     &		,rhov(1:lmax,k),visv(0:lmax,k),difv(0:lmax,k)   &   !ivb
     &          ,tken(0:lmax,k),eps(0:lmax,k))                  !ivb         

        end do

        end

!**************************************************************

        subroutine keps_init

! initializes arrays for keps routine

        use turbulence_util
        use diffusion
        use levels, only : nlvdi,nlv
        use basin, only : nkn,nel,ngr,mbw

        implicit none

        double precision kmin,epsmin,lenmin,avumol,avtmol,avsmol
!       parameter(kmin=1.e-10,epsmin=1.e-12,lenmin=0.01)
        parameter(kmin=3.e-6,epsmin=5.e-10,lenmin=0.01)
        parameter(avumol=1.3e-6,avtmol=1.4e-7,avsmol=1.1e-9)

        integer k,l

        do k=1,nkn
          do l=0,nlvdi
            tken(l,k)   = kmin
            eps(l,k)    = epsmin
            rls(l,k)    = lenmin
            visv(l,k)   = avumol
            difv(l,k)   = avsmol
          end do
        end do

        end

!**************************************************************
!**************************************************************
! debug routines
!**************************************************************
!**************************************************************

	subroutine write_node_vel_info(iunit,it,k,ndim,depth,u_taus,u_taub,z0s,z0b      &
     &			,num,nuh,ken,dis,len,nn,ss,hh)

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

!**************************************************************

	subroutine write_node_bin_info(iunit,it,k,ndim,depth,u_taus,u_taub,z0s,z0b      &
     &			,num,nuh,ken,dis,len,nn,ss,hh)

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

!**************************************************************

!--------------------------------------------------------------------
        end module shy_turbulence
!--------------------------------------------------------------------
