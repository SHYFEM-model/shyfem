
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    This file has been adapted from GOTM.
!    Please see also the original copyright notice in the GOTM distribution.

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
! 07.06.2011	ggu	changed VERS_6_1_25
! 12.12.2014	ggu	changed VERS_7_0_9
! 19.01.2015	ggu	changed VERS_7_1_3
! 25.05.2016	ggu	changed VERS_7_5_10
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

c internal gotm routines
c
c these are routines that should not be used anymore
c please use the routines directly from the GOTM library
c
c********************************************************************

	subroutine has_gotm(bgotm)

	implicit none

	logical bgotm

	bgotm = .false.

	end

c********************************************************************

	subroutine init_gotm_turb(iunit,fn,ndim)

	implicit none

	integer iunit
	character*(*) fn
	integer ndim

	call gotmturb_init 

	end

c********************************************************************

        subroutine do_gotm_turb   (
     &                            Nmx
     &                           ,dt,depth
     &                           ,u_taus,u_taub
     &                           ,z0s,z0b
     &                           ,hh,nn,ss
     &                           ,num,nuh
     &                           ,ken,dis,len
     &                           )

	implicit none

        integer Nmx 
        double precision dt,depth
        double precision u_taus,u_taub
        double precision z0s,z0b
        double precision hh(0:Nmx),NN(0:Nmx),SS(0:Nmx)
        double precision num(0:Nmx),nuh(0:Nmx)
        double precision ken(0:Nmx),dis(0:Nmx),len(0:Nmx)

        call gotmturb   (
     &                            Nmx,dt,hh
     &                           ,nn,ss
     &                           ,num,nuh,ken,dis,len
     &                           ,u_taus,u_taub
     &                  )

	end

c********************************************************************
c********************************************************************
c********************************************************************

!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GOTMTURB() - the turbulence part of GOTM 
!
! !INTERFACE:
      subroutine gotmturb(Nmx,dt,h,NN,SS,num,nuh,k,eps,L,
     &                    u_taus,u_taub)

      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
      integer Nmx 
      double precision h(0:Nmx),NN(0:Nmx),SS(0:Nmx),num(0:Nmx),
     &                 nuh(0:Nmx),k(0:Nmx),eps(0:Nmx),L(0:Nmx)
      double precision u_taus,u_taub,dt
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
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
!  14Nov00   Changes due to separation of turbulent part 
!
! !LOCAL VARIABLES:
      integer i
      logical FIRST 
      save    FIRST
      data    FIRST / .true./
      double precision xRf(0:MaxN),ko(0:MaxN),depth    
      double precision P(0:MaxN),B(0:MaxN),cmue1(0:MaxN),cmue2(0:MaxN)    

	logical istdebug!ggu
	logical bdebug	!ggu
	integer j	!ggu
	integer nmlast	!ggu
	integer istep	!ggu
	save istep
	data istep / 0 /
!
!EOP
!-------------------------------------------------------------------------
!BOC

	nmlast = 30	!ggu
	bdebug = .true.		!ggu
	bdebug = istdebug()	!ggu
	!bdebug = .false.

      if (FIRST) then 
         call gotmturb_init
         FIRST=.false.
      end if

	if( Nmx .le. 1 ) return	!ggu

	z0b = 0.03 * 0.05	!ggu
	z0s = 0.03 * 0.05	!ggu

      if (Nmx.gt.MaxN) then
         write(*,*) 'MaxN is smaller than Nmx.'
         write(*,*) 'Increase MaxN in gotmturb.i or decrease Nmx.'
         write(*,*) 'Program interrupted in gotmturb.f'
         stop
      end if 

      do i=1,Nmx-1
         if (k(i).lt.k_min) k(i)=k_min
         if (eps(i).lt.epsmin) eps(i)=epsmin
         if (L(i).lt.L_min) L(i)=L_min
      end do 
      
      depth=0. 
      do i=1,Nmx
         depth=depth+h(i)
      end do

	if( bdebug ) then	!ggu

c	if( istep .eq. 0 ) then
c	write(95,*) 'depth...' 
c	doi=nmx,1,-1
c	  j = Nmx - i + 1
c	  write(95,'(i3,7e11.3)') j,h(i)
c	end do
c	end if

c	call write_common(99)
c	call write_dt(91,dt)

c	write(95,*) istep,depth,u_taus,u_taub
c	write(95,*) z0b,z0s
c	write(95,*) k_min,epsmin,L_min
c	write(95,*) Nmx,dt
c	do i=nmx,nmlast,-1
c	  j = Nmx - i + 1
c	  write(95,'(i3,7e11.3)') j,num(i),nuh(i)
c     +			,k(i),eps(i),L(i),ss(i),nn(i)
c	end do
	istep = istep + 1
	end if

      call StabilityFunctions(Nmx,k,L,SS,NN,
     &                        cmue1,cmue2,P,B,num,nuh,xRf)

	call pr_info('Stab',91,istep,Nmx,nmlast
     +			,num,nuh,k,eps,L,ss,nn)

      call TKE(Nmx,dt,cmue1,cmue2,num,P,B,eps,L,k,ko,h,u_taus,u_taub,
     &          NN,SS)

	call pr_info('TKE',91,istep,Nmx,nmlast
     +			,num,nuh,k,eps,L,ss,nn)

      call LengthScale(Nmx,dt,k,ko,eps,L,num,nuh,NN,SS,h,
     &                  cmue1,cmue2,P,B,u_taus,u_taub,depth,xRf)

	call pr_info('Len',91,istep,Nmx,nmlast
     +			,num,nuh,k,eps,L,ss,nn)

      call KolPran(Nmx,cmue1,cmue2,k,L,num,nuh,u_taus,u_taub)

	call pr_info('KolP',91,istep,Nmx,nmlast
     +			,num,nuh,k,eps,L,ss,nn)

      call InternalWave(Nmx,num,nuh,NN,SS,k)

	call pr_info('IWave',91,istep,Nmx,nmlast
     +			,num,nuh,k,eps,L,ss,nn)

c	if( bdebug ) then	!ggu
c	write(96,*) istep,depth,u_taus,u_taub
c	write(96,*) z0b,z0s
c	write(96,*) k_min,epsmin,L_min
c	write(96,*) Nmx,dt
c	do i=nmx,nmlast,-1
c	  j = Nmx - i + 1
c	  write(96,'(i3,7e11.3)') j,num(i),nuh(i)
c     +			,k(i),eps(i),L(i),ss(i),nn(i)
c	end do
c	end if

      return
      end
!EOC

c**********************************************************************

	subroutine write_vars(iunit,it,low,high,fact,var)

c writes var to file

	implicit none

	integer iunit,it,low,high
	real fact
	double precision var(0:1)

	integer l

        write(iunit,*) it,high-low+1,fact
        write(iunit,*) (var(l),l=high,low,-1)                                        

	end

c**********************************************************************

	subroutine write_dt(iunit,dt)

	implicit none

	integer iunit
	double precision dt
	logical bdebug
	logical istdebug

	bdebug = istdebug()

	if( .not. bdebug ) return

	write(iunit,*) 'time step : ',dt

	end

c**********************************************************************

	subroutine write_common(iunit)

	implicit none

	integer iunit
	logical bdebug
	logical istdebug

      double precision kappa,z0b,
     &              cm0,z0s,ce1,ce2,ce3minus,ce3plus,
     &              k_min,epsmin,L_min,
     &              sig_k,sig_e,Prandtl0,cde,cdL,cmucst,
     &              sl,e1,e2,e3,a1,a2,b1,b2,c1,
     &              galp,qeghmin,qeghmax,qeghcrit,
     &              alfa,c2,c3,klimiw,rich_cr,numiw,nuhiw,numshear   
      integer MaxN,Stab,MYLength,TKEMeth,LengthMeth,Iwmodel
      logical fluxcond,lengthlim  
      logical qesmooth

      common /turbconstants/  kappa,z0b,
     &                    cm0,z0s,ce1,ce2,ce3minus,ce3plus,
     &                    k_min,epsmin,L_min,
     &                    sig_k,sig_e,Prandtl0,cde,cdL,
     &                    sl,e1,e2,e3,a1,a2,b1,b2,c1,
     &                    galp,cmucst,qeghmin,qeghmax,
     &                    qeghcrit,
     &                    alfa,c2,c3,klimiw,rich_cr,numiw,nuhiw,
     &                    numshear,Stab,TKEMeth,LengthMeth,IwModel,
     &                    MYLength,
     &                    fluxcond,lengthlim,  
     &                    qesmooth 

	bdebug = istdebug()

	if( .not. bdebug ) return

	write(iunit,*) '-------------------------------'
	write(iunit,'(e14.6)') kappa,z0b,
     &                    cm0,z0s,ce1,ce2,ce3minus,ce3plus,
     &                    k_min,epsmin,L_min,
     &                    sig_k,sig_e,Prandtl0,cde,cdL,
     &                    sl,e1,e2,e3,a1,a2,b1,b2,c1,
     &                    galp,cmucst,qeghmin,qeghmax,
     &                    qeghcrit,
     &                    alfa,c2,c3,klimiw,rich_cr,numiw,nuhiw,
     &                    numshear
	write(iunit,'(i10)') Stab,TKEMeth,LengthMeth,IwModel,
     &                    MYLength
	write(iunit,*) fluxcond,lengthlim,qesmooth 
	write(iunit,*) '-------------------------------'

	end

c**********************************************************************

	subroutine pr_info(text,iunit,istep,Nmx,nmlast
     +				,num,nuh,k,eps,L,ss,nn)

	implicit none

	character*(*) text
	logical bdebug
	integer iunit,istep,nmx,nmlast
      double precision NN(0:Nmx),SS(0:Nmx),num(0:Nmx),
     &                 nuh(0:Nmx),k(0:Nmx),eps(0:Nmx),L(0:Nmx)

	integer i,j
	integer last
	logical istdebug

	bdebug = istdebug()

	if( .not. bdebug ) return

	last = Nmx - nmlast
	last = max(0,last)

	write(iunit,*) text
	write(iunit,*) text
	write(iunit,*) text
	write(iunit,*) istep,nmx,last
	do i=nmx,last,-1
	  j = Nmx - i + 1
	  write(iunit,'(i3,7e11.3)') j,num(i),nuh(i)
     +			,k(i),eps(i),L(i),ss(i),nn(i)
	end do

	end

c***************************************************************
 
        function istdebug()
 
        implicit none
 
        logical istdebug
 
        istdebug = .true.
        istdebug = .false.
 
        end
 
c***************************************************************                
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GOTMTURB_INIT() - initialisation of the turbulence part of GOTM  
!
! !INTERFACE:
      subroutine gotmturb_init 

      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
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
!  14Nov00   Changes due to separation of turbulent part
!  22Feb00   Hardcoded namelist parameters introduced (ggu)
!
! !LOCAL VARIABLES:
      integer NAMLST
      parameter(NAMLST=29) 
      logical NMLIST
      data NMLIST / .false. /	!read from namelist (true) or hardcoded (false)
!
!EOP
!-------------------------------------------------------------------------
!BOC

!     +---------------------------------------------------------------+
!     |  Namelist specifications                                      |
!     +---------------------------------------------------------------+
!     General turbulence settings
      namelist /turb/
     &          TKEMeth,LengthMeth,lengthlim,k_min,epsmin,L_min, 
     &          Stab,Prandtl0,cmucst,galp,kappa,cm0    

!     Parameters for k-eps Model
      namelist /keps/   ce1,ce2,ce3minus,ce3plus,sig_k,fluxcond

!     Parameters for MY Model
      namelist /my/     sl,e1,e2,e3,MYLength

!     Parameters quasi-equilibrium stability function
      namelist /stabfunc/
     &         a1,a2,b1,b2,c2,c3,qesmooth,qeghmax,qeghmin,qeghcrit

!     Parameters internal wave model.
      namelist /iw/ IwModel,alfa,klimiw,rich_cr,numiw,nuhiw,numshear

!     +---------------------------------------------------------------+
!     |  Open and read name list parameters                           |
!     +---------------------------------------------------------------+

      if( NMLIST ) then
         open(NAMLST,status='unknown',file='gotmturb.inp')
         read(NAMLST,turb)
         read(NAMLST,keps)
         read(NAMLST,my)
         read(NAMLST,stabfunc)
         read(NAMLST,iw)
         close(NAMLST) 
      else
         call hardlist
      end if

         cde   = cm0*cm0*cm0
         sig_e = kappa*kappa*cm0/(ce2-ce1)/cde
         c1=(1.-b1**(-1./3.)/a1-6*a1/b1)/3. !See Kantha & Clayson 1994, eq. (23)
         if ((LengthMeth.eq.Ispramix).and.(Stab.ne.FluxRich)) then
             write(*,*) 'You have chosen LenthMeth = ',Ispramix,'.'
             write(*,*) 'Please, choose then Stab =',FluxRich,'.'
             write(*,*) 'Program aborted.'
             stop
          end if 
          if ((Stab.lt.1).or.(Stab.gt.12)) then
             write(*,*) 'Stab should be between 1 and 12.'
             write(*,*) 'Program aborted.'
             stop 
          end if 

      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HARDLIST() - hardcoded namelist parameters
!
! !INTERFACE:
      subroutine hardlist

      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
! hardcoded namelist inputs
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
!  22Feb00   newly introduced (ggu)
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
!BONL
!-------------------------------------------------------------------------
!         JRC/SAI/ME Generel Ocean Turbulence Model (GOTM)               !
!-------------------------------------------------------------------------
!
! ! namelist.inp: Model configuration file 
!
! TKEMeth    : Method for computing TKE (Turbulent Kinetic Energy)
!     : 1, Algebraic equation
!     : 2, Dynamic equation for k-epsilon model
!     : 3, Dynamic equation for Mellor-Yamada model
! LengthMethod : Method for computing turbulent length scale L
!     : 1, Parabolic shape
!     : 2, Triangle shape
!     : 3, Xing and Davies [1995]
!     : 4, Robert and Ouellet [1987] 
!     : 5, Blackadar (two boundaries) [1962] 
!     : 6, Bougeault and Andre [1986] 
!     : 7, Eifler and Schrimpf (ISPRAMIX) [1992] 
!     : 8, Dynamic dissipation rate equation
!     : 9, Dynamic Mellor-Yamada kL equation
! lengthlim : Limitation of L for stable stratification, Galperin et al. [1988]
!     : .true.,  Length limitation on  
!     : .false., Length limitation off   
! k_min : Minimum value of TKE [J/kg]
! epsmin: Minimum value of Dissipation rate [W/kg] 
! L_min : Minimum value of turbulent length scale [m] 
! Stab: Method for computing stability functions (SF)
! Note that the given values for cm0,cmust,Prandtl0 are recommendations
! For values for ce3minus, see below.
!     : 1, Kantha and Clayson [1994],      full version, cm0 = 0.5544
!     : 2, Burchard and Baumert [1995],    full version, cm0 = 0.5900
!     : 3, Canuto et al. [2000] version A, full version, cm0 = 0.5270
!     : 4, Canuto et al. [2000] version B, full version, cm0 = 0.5540
!     : 5, Kantha and Clayson [1994],      quasi-eq. version, cm0 = 0.5544
!     : 6, Burchard and Baumert [1995],    quasi-eq. version, cm0 = 0.5900
!     : 7, Canuto et al. [2000] version A, quasi-eq. version, cm0 = 0.5270
!     : 8, Canuto et al. [2000] version B, quasi-eq. version, cm0 = 0.5540
!     : 9, Constant stability functions,   cm0 = cmust = 0.5477, Prandtl0=0.74
!     :10, Munk and Anderson [1954],       cm0 = cmust = 0.5477, Prandtl0=0.74
!     :11, Schumann and Gerz [1995],       cm0 = cmust = 0.5477, Prandtl0=0.74
!     :12, Eifler and Schrimpf [1992],     cm0 = cmust = 0.5477, Prandtl0=0.74
! cmucst   : Value of SF cmue1 if Stab>=9, for LengthMeth=6 set cmucst=0.1
! Prandtl0 : Neutral turbulent Prandtl number, only relevant if Stab>=9
! galp     : Coefficient for length limitation, should be 0.53  
! kappa : von Karman constant
! cm0   : Stability function for momentum for unstratified flow, 
!         Galperin et al. [1988], important for relation between k, L and eps.
!
!&turb
      TKEMeth    = 2 
      LengthMeth = 8 
      lengthlim  = .false.
      k_min      = 3e-6
      epsmin     = 5e-10
      L_min      = 0.01
      Stab       = 3 
      Prandtl0   = 0.714
      cmucst     = 0.1
      galp       = 0.53
      kappa      = 0.4
      cm0        = 0.527 
!&end
! 
! ce1     : Empirical coefficient in dissipation equation
! ce2     : Empirical coefficient in dissipation equation
! ce3minus: Empirical coefficient in dissipation equation (stable strat.)
!   Recommended values for ce3minus (steady-state Richardson number=0.25) are:
!   Stab = 1 --> ce3minus = -0.404
!   Stab = 2 --> ce3minus = -0.444
!   Stab = 3 --> ce3minus = -0.629
!   Stab = 4 --> ce3minus = -0.566
!   Stab = 5 --> ce3minus = -0.404
!   Stab = 6 --> ce3minus = -0.444
!   Stab = 7 --> ce3minus = -0.629
!   Stab = 8 --> ce3minus = -0.566
!   Stab = 9 --> ce3minus = +0.499
!   Stab =10 --> ce3minus = +0.035
!   Stab =11 --> ce3minus = -0.368
!   Stab =12 --> ce3minus = +0.239
! ce3plus : Empirical coefficient in dissipation equation (unstable strat.)
! sig_k   : Schmidt number for TKE eddy diffusivity
! fluxcond: Switch for type of boundary condtions for k and eps 
!         : .true.,  flux boundary conditions   
!         : .false., Dirichlet boundary conditions   
!&keps
      ce1        = 1.44
      ce2        = 1.92
      ce3minus   = -0.629
      ce3plus    = 1.0  
      sig_k      = 1.
      fluxcond   = .true. 
!&end
!
! sl : Parameter for calculating eddy diffusivities of k and kL (sl=cl/sqrt(2))
! e1 : Coefficient in Mellor-Yamada kL equation 
! e2 : Coefficient in Mellor-Yamada kL equation 
! e3 : Coefficient in Mellor-Yamada kL equation, see Burchard 2000, JPO.  
! MYLength: Prescribed barotropic lengthscale in Mellor-Yamada kL equation
!     : 1, Parabolic shape
!     : 2, Triangle shape
!     : 3, Linear from surface, infinit depth  
!
!&my
      sl         = 0.2
      e1         = 1.8
      e2         = 1.33
      e3         = 5.1
      MYLength   = 3
!&end
!
! a1,a2,b1,b2,c1 : Coefficients in Galperin quasi-equilibrium SF for Stab=2  
! qesmooth: .true. smoothing of SF for unstable stratification for Stab=2
! qeghmax : Maximum value of stratification parameter gh for Stab=2
! qeghmin : Minimum value of stratification parameter gh for Stab=2
! qeghcrit: Critical value of gh to start smoothing for Stab=2
!
!&stabfunc
      a1         = 0.92     
      a2         = 0.74 
      b1         = 16.6
      b2         = 10.1
      c2         = 0.0
      c3         = 0.0
      qesmooth   = .true. 
      qeghmax    = 0.0233 
      qeghmin    = -0.28 
      qeghcrit   = 0.02
!&end 
!
! IWModel : Method for modelling internal wave effects on mixing
!    : 0, No IW Model (except from possible limitations of TKE and epsilon)
!    : 1, Mellor [1989] 
!    : 2, Large et al. [1994] 
! alfa    : Coefficient for Mellor [1989] IW model
! klimiw  : Critical value of TKE in Large et al. [1994] IW model
! rich_cr : Critical Richardson number for shear instability       
! numshear: Background diffusivity for shear instability           
! numiw   : Background viscosity for internal wave breaking        
! nuhiw   : Background viscosity for internal wave breaking       
!
!&iw
      IWModel    = 2
      alfa       = 0.7
      klimiw     = 1e-6
      rich_cr    = 0.7
      numiw      = 1.e-4
      nuhiw      = 1.e-5
      numshear   = 5.e-3
!&end
!EONL

	z0b = 0.03 * 0.05
	z0s = 0.03 * 0.05

	end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: AlgebraicLength() - diff. algebraic length scales.
!
! !INTERFACE:
      subroutine AlgebraicLength(Nmx,k,ko,eps,L,NN,SS,h,
     &                     u_taus,u_taub,depth,xRf)
      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
      integer Nmx
      double precision k(0:Nmx),ko(0:Nmx),
     &                 NN(0:Nmx),SS(0:Nmx),
     &                 h(0:Nmx),
     &                 u_taus,u_taub,depth,xRf(0:Nmx)
!
!     kappa is read from namelist.inp and declared in const.i
! !OUTPUT PARAMETERS:
      double precision L(0:Nmx),eps(0:Nmx)
!
! !DESCRIPTION:
!     These subroutine computes the vertical profile of the mixing lengthscale
!     from an algebraic relation.
!
!     1) Parabola
!     \begin{equation}
!     l=\kappa \frac {l_b l_s} {l_b+l_s}
!     \end{equation}
!     2) triangle
!     \begin{equation}
!     l=\kappa max(l_b,l_s)
!     \end{equation}
!     3) Distorted Parabola. Robert-Ouellet
!     \begin{equation}
!     l=\kappa l_b (1-\frac {l_b} {D})
!     \end{equation}
!     4) Xing, parabolic with exponential $d_b$
!     \begin{equation}
!     l_b=\kappa l_b e^{-\beta l_b}
!     \end{equation}
!     5) Blackadar (in the presence of two boundaries)
!     \begin{equation}
!     l=\kappa \frac {1} {\frac {1} {l_b}+\frac {1} {l_s}+\frac {1} {l_a}} \\
!     l_a=\gamma_0 \frac {\int_{0}^{D} k^{1/2} zdz} {\int_{0}^{D} k^{1/2} dz}
!     \end{equation}
!     6) see ispralength.f
!
!     At the end of the subroutine, the dissipation rate is calculated using:
!
!     \begin{equation}\label{DefDissip}
!     \varepsilon = (c_{\mu}^0)^3 \frac{k^{3/2}}{L}. 
!     \end{equation}
!
! !BUGS:
!
! !SEE ALSO:
!     potentialml()  
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  25Nov98   The GOTM developers  Initial code.
!  14Nov00   Changes due to separation of turbulent part
!
! !LOCAL VARIABLES:
      integer i
      double precision ds,db,dbxing
      double precision beta,gamma,La,La_up,La_down   
      double precision Lcrit
!
!EOP
!-------------------------------------------------------------------------
!BOC
      db=0.0
!
! Parabola shape
!       
      if (LengthMeth.eq.Parabola) then
        do i=1,Nmx-1
          db=db+h(i)
          ds=depth-db
          L(i)=kappa*(ds+z0s)*(db+z0b)/(ds+db+z0b+z0s)
        end do
      end if
!
! Triangle shape
!
      if (LengthMeth.eq.Triangle) then 
        do i=1,Nmx-1
          db=db+h(i)
          ds=depth-db
          L(i)=kappa*min(ds+z0s,db+z0b)
        end do
      end if 
!  
! Xing and Davies (1995). 
! Modification of parabolic mixing length. db changes:
!
      if (LengthMeth.eq.Xing) then
        beta=-2.
        do i=1,Nmx-1
          db=db+h(i)
          ds=depth-db
          dbxing=db*dexp(-beta*db)
          L(i)=kappa*(ds+z0s)*(dbxing+z0b)/(ds+dbxing+z0s+z0b)
        end do
      end if
!
! Robert and Ouellet(1987). Similar to parabolic
!
      if (LengthMeth.eq.RobertOuellet) then
        do i=1,Nmx-1
          db=db+h(i)
          ds=depth-db
          L(i)= kappa*(db+z0b)*sqrt((ds+z0s)/(ds+db+z0b+z0s)) 
        end do
      end if
!
! Blackadar (1962). 
! In the form suggested by Luyten et al. (1996) for two boundary layers.
!
      if (LengthMeth.eq.Blackadar) then
        La_up=0.
        La_down=0.
        do i=1,Nmx-1
          db=db+h(i) 
          La_up=La_up+sqrt(ko(i))*(db+z0b)*h(i)
          La_down=La_down+sqrt(ko(i))*h(i) 
        end do
        gamma=0.2
        La=gamma*La_up/La_down
        db=0.0
        do i=1,Nmx-1
          db=db+h(i)
          ds=depth-db
          L(i)=1/(1/(kappa*(ds+z0s))+1/(kappa*(db+z0b))+1/La)
        end do
      end if  
!
!  Ispramix
!
      if (LengthMeth.eq.Ispramix) then
        call  Ispralength(Nmx,k,ko,eps,L,NN,SS,h,
     &                     u_taus,u_taub,depth,xRf) 
      end if

! Boundary conditions for L
!
      L(0)=kappa*z0b
      L(Nmx)=kappa*z0s 
 
      do i=0,Nmx
        if ((NN(i).gt.0).and.(lengthlim)) then
          Lcrit=sqrt(2*galp*galp*k(i)/NN(i))
          if (L(i).gt.Lcrit) L(i)=Lcrit
        end if
        if (L(i).lt.L_min) L(i)=L_min
        eps(i)=cde*sqrt(k(i)*k(i)*k(i))/L(i)
        if (eps(i).lt.epsmin) eps(i)=epsmin
      end do  
      return 
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DissipationEq() - the dissipation part of the k-$\varepsilon$ model 
!
! !INTERFACE:
      subroutine DissipationEq(Nmx,dt,k,ko,eps,L,num,nuh,NN,h,
     &                    cmue1,cmue2,P,B,u_taus,u_taub)

      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!     Tridiagonal()
!
! !INPUT PARAMETERS:
      integer Nmx 
      double precision cmue1(0:Nmx),P(0:Nmx),B(0:Nmx),k(0:Nmx),
     &                 h(0:Nmx),
     &                 num(0:Nmx),nuh(0:Nmx),cmue2(0:Nmx),ko(0:Nmx),
     &                 NN(0:Nmx)
      double precision u_taus,u_taub,dt 
!
! !OUTPUT PARAMETERS:
      double precision eps(0:Nmx),L(0:Nmx)
!
! !DESCRIPTION:
!     This subroutine calculates the dissipation rate in the framework
!     of the k-epsilon model:
!
!     \begin{equation}\label{eps_eq}
!     \partial_t \varepsilon -
!     \partial_z(\nu_{\varepsilon}\partial_z \varepsilon) =
!     \frac{\varepsilon}{k} \left(c_{\varepsilon 1}P + c_{\varepsilon 3}B - c_{\varepsilon 2}\varepsilon \right).
!     \end{equation}
!
!     As boundary condtions a choice between Dirichlet (fluxcond=.false.)
!     and Neumann flux conditions (fluxcond=.true.) has to be made.   
!
!     Dirichlet conditions:
!
!     \begin{equation}\label{Standardeps}
!     \varepsilon =
!     \left( c_{\mu}^0 \right)^3 \frac{k^{3/2}}{\kappa (\tilde z + z_0)}.
!     \end{equation}
!
!     Neumann flux conditions:
!
!     \begin{equation}\label{Fluxeps}
!     \frac{\nu_t}{\sigma_{\varepsilon}} \partial_{\tilde z} \varepsilon =
!     -\left( c_{\mu}^0 \right)^3
!     \frac{\nu_t}{\sigma_{\varepsilon}}
!     \frac{k^{3/2}}{\kappa (\tilde z + z_0)^2}. 
!     \end{equation}
!
!     At the end of the subroutine, the Galperin et al. [1988] limitation
!     and the calculation of the macro length scale is carried out. 
!
! !BUGS:
!
! !SEE ALSO:
!     algebraiclength(), ispralength(), potentialml()   
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  25Nov98   The GOTM developers  Initial code.
!  14Nov00   Changes due to separation of turbulent part
!
! !LOCAL VARIABLES:
      double precision avh(0:MaxN),au(0:MaxN),
     &                 bu(0:MaxN),cu(0:MaxN),
     &                 du(0:MaxN),flux(0:MaxN),
     &                 pminus(0:MaxN),pplus(0:MaxN)
      double precision Prod,Buoy,Diss,cee3,epslim 
      integer i
!
!EOP
!-------------------------------------------------------------------------
!BOC

      do i=1,Nmx 
         avh(i)=0.5/sig_e*(num(i-1)+num(i))
      end do 
      do i=1,Nmx-1 
         flux(i)=0 
      end do 

      if (fluxcond) then
         flux(1    )=avh(1  )*cde*(ko(1    )**1.5)
     &              /(kappa*(z0b+0.5*h(1  ))**2.)
         flux(Nmx-1)=avh(Nmx)*cde*(ko(Nmx-1)**1.5)
     &              /(kappa*(z0s+0.5*h(Nmx))**2.)
! A bug in the previous two lines has been found
! by Patrick Luyten, MUMM, Belgium. kappa had been squared as well before.
! See the GOTM report, 1999 for the correct mathematical formulation.
         avh(1  )=0
         avh(Nmx)=0
      else
         avh(1  )=u_taub**4*2/sig_e/(eps(0)+eps(1))
         avh(Nmx)=u_taus**4*2/sig_e/(eps(Nmx)+eps(Nmx-1))
      end if 
 
      do 202 i=1,Nmx-1
         if (B(i).gt.0) then
            cee3=ce3plus 
         else
            cee3=ce3minus 
         end if
         Prod=ce1*eps(i)/ko(i)*P(i)
         Buoy=cee3*eps(i)/ko(i)*B(i)
         Diss=ce2*eps(i)*eps(i)/ko(i)
         if (Prod+Buoy.gt.0) then
            pplus(i)=Prod+Buoy
            pminus(i)=Diss
         else
            pplus(i)=Prod
            pminus(i)=Diss-Buoy
         end if
202   continue
 
      do 203 i=1,Nmx-1
         au(i)=-2.*dt*avh(i)/(h(i)+h(i+1))/h(i)
         cu(i)=-2.*dt*avh(i+1)/(h(i)+h(i+1))/h(i+1)
         bu(i)=1.-au(i)-cu(i)+pminus(i)*dt/eps(i)
         du(i)=(1+pplus(i)*dt/eps(i))*eps(i)
     &         +flux(i)*dt/(0.5*(h(i)+h(i+1)))
203   continue

      if (fluxcond) then
         call Tridiagonalx(Nmx,1,Nmx-1,au,bu,cu,du,eps)
         eps(0  ) = cde*sqrt(k(0  )*k(0  )*k(0  ))/kappa/z0b 
         eps(Nmx) = cde*sqrt(k(Nmx)*k(Nmx)*k(Nmx))/kappa/z0s 
      else

         cu(0)=0
         bu(0)=1.
         du(0)=cde*sqrt(k(0)*k(0)*k(0))/kappa/z0b 
 
         bu(Nmx)=1.
         au(Nmx)=0
         du(Nmx)=cde*sqrt(k(Nmx)*k(Nmx)*k(Nmx))/kappa/z0s 

         call Tridiagonalx(Nmx,0,Nmx,au,bu,cu,du,eps)
      end if 

 
      do i=0,Nmx
         if ((NN(i).gt.0).and.(lengthlim)) then
            epslim=cm0*cm0*cm0/sqrt(2.)/galp*k(i)*sqrt(NN(i)) 
         else
            epslim=epsmin 
         end if 
         if (eps(i).lt.epslim) eps(i)=epslim
         L(i)=cde*sqrt(k(i)*k(i)*k(i))/eps(i)
      end do

      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cmue_bb() - Burchard and Baumert [1995] stability functions.
!
! !INTERFACE:
      subroutine cmue_bb(aas,aan,cmue1,cmue2)
      implicit none
 
!
! !USES:
!
! !INPUT PARAMETERS:
      double precision aas,aan
!
! !OUTPUT PARAMETERS:
      double precision cmue1,cmue2
!
! !DESCRIPTION:
!     This subroutine computes Burchard and Baumert [1995] stability functions.
!
!
! !BUGS:
!
! !SEE ALSO:
!     cmue_kc(),cmue_ca(),cmue_cb(),cmue_kcqe(),cmue_bbqe(),
!     cmue_caqe(),cmue_cbqe(),cmue_cn(),cmue_ma(),cmue_sg(),cmue_rf()
! 
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  18Nov99  Introduced for GOTM2.0.0 
!
! !LOCAL VARIABLES:
      double precision an,as
      double precision k(1:5),kk(1:4),X,V,D,A,E,H,ww 
      integer i,j,count
      double precision c1,c2,c3,c1t,c2t,c3t,ct
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,cm0
!
!EOP
!-------------------------------------------------------------------------
!BOC

      c1=1.8
      c2=0.6
      c3=0.6
      c1t=3.0
      c2t=0.33
      c3t=0.33
      ct=1.6

      a1=c1-1.0
      a2=c1-1.0
      a3=c1t-0.5
      a4=c1t-0.5
      a5=c2
      a6=c2
      a7=1-c3
      a8=1-c3
      a9=1-c2t
      a10=1-c3t
      a11=1.5-c3

      cm0 = 0.5900 

      as = 1./cm0**6 * aas
      an = 1./cm0**6 * aan


      k(1)= a1+a2+2*(a3+a4+a10*ct*an)
      k(2)=4.*a3*(a4+a10*ct*an)
     &     +2.*(a1+a2)*(a3+a4+a10*ct*an)+a1*a2
     &     +8./3.*a11*an+2*a7*an-2./3.*(1.-a5)*a6*as
      k(3)=16./3.*a11*(0.5*a1+a3)*an
     &     +4.*a7*(0.5*a2+a4+a10*ct*an)*an
     &     -4./3.*a6*(1.-a5)*(a3+a4+a10*ct*an)*as
     &     +4.*a3*(a1+a2)*(a4+a10*ct*an)   
     &     +2.*a1*a2*(a3+a4+a10*ct*an)
     &     -8./3.*(c1-1.)*(0.25*as*(1.-a5)-0.5*an)
      k(4)=16./3.*a1*a3*a11*an+4.*a2*a7*(a4+a10*ct*an)*an
     &     -8./3.*a3*a6*(1.-a5)*(a4+a10*ct*an)*as
     &     +16./3.*a7*a11*an*an
     &     +8./3.*a6*a7*a9*as*an+4.*a1*a2*a3*(a4+a10*ct*an)
     &     -8./3.*(c1-1.)*(0.5*(1.-a5)*(a3+a4+a10*ct*an)*as
     &     -(0.5*a1+a3)*an)
      k(5)=-8./3.*(c1-1.)*(((1.-a5)*a3*(a4+a10*ct*an)
     &     -a7*a9*an)*as-a7*an*an-a1*a3*an)



c      k(1) = 11.60 +  2.14*an
c      k(2) = 41.64 + 17.35*an - 0.160*as
c      k(3) = 46.40 + 38.15*an - 1.813*as + 1.715*an*an - 0.343*an*as
c      k(4) = 16.00 + 25.85*an - 6.133*as + 3.292*an*an - 1.744*an*as
c      k(5) =          4.27*an - 5.333*as + 0.853*an*an - 1.715*an*as

      kk(1) = 4*k(1)      
      kk(2) = 3*k(2)      
      kk(3) = 2*k(3)      
      kk(4) =   k(4)      

      
      X=0.25*as 
      if (an.lt.0.) X=X-2.0*an 

      count =0
111   V=1.0
      count=count+1
      do i=1,5
         V=V*X+k(i) 
      end do  

      D=5.0
      do i=1,4
         D=D*X+kk(i) 
      end do  

      X=X-V/D

      if (abs(V/D).gt.1e-4) goto 111

      X=X-1.

      A=1./(1.+1.072*an/(3.0+0.5*X))
      D=1.+0.4*an/(1.8+X)/(3.0+0.5*X)
      E=1.+1.2*A*an/(1.8+X)/(3.0+0.5*X)
      H=1.-0.67*A*an/(3.0+0.5*X)/(3.0+0.5*X)
      ww=2./3.*0.8/(E*(1.8+X)-0.16/(1.8+X)*H/D*as)

      cmue1 = 0.4/(1.8+X)*H/D*ww
      cmue2 = A/(3.0+0.5*X)*ww

      cmue1 = 1./cm0**3 * cmue1
      cmue2 = 1./cm0**3 * cmue2


      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cmue_bbqe() - Equilibrium version of the Burchard and Baumert [1995]!                         stability functions.
!
! !INTERFACE:
      subroutine cmue_bbqe(aan,cmue1,cmue2)
      implicit none
 
!
! !USES:
!
! !INPUT PARAMETERS:
      double precision aan
!
! !OUTPUT PARAMETERS:
      double precision cmue1,cmue2
!
! !DESCRIPTION:
!     This subroutine computes the equilibrium version of the 
!     Burchard and Baumert [1995] stability functions.
!
! !BUGS:
!
! !SEE ALSO:
!     cmue_kc(),cmue_bb(),cmue_ca(),cmue_cb(),cmue_kcQe(),
!     cmue_caQe(),cmue_cbQe(),cmue_cn(),cmue_ma(),cmue_sg(),cmue_rf()
! 
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  18Nov99  Introduced for GOTM2.0.0 
!
! !LOCAL VARIABLES:
      double precision an,as 
      double precision k(1:5),kk(1:4),X,V,D,A,E,H,ww 
      integer i,j,count
      double precision c1,c2,c3,c1t,c2t,c3t,ct
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,cm0
      double precision Prinv,xxx,err 
!
!EOP
!-------------------------------------------------------------------------
!BOC



      c1=1.8
      c2=0.6
      c3=0.5
      c1t=3.0
      c2t=0.33
      c3t=0.33
      ct=1.6

      a1=c1-1.0
      a2=c1-1.0
      a3=c1t-0.5
      a4=c1t-0.5
      a5=c2
      a6=c2
      a7=1-c3
      a8=1-c3
      a9=1-c2t
      a10=1-c3t
      a11=1.5-c3


      cm0 = 0.5900   

      an = 1./cm0**6 * aan

      X=0.0

      A=1./(1.+1.072*an/(3.0+0.5*X))
      D=1.+0.4*an/(1.8+X)/(3.0+0.5*X)
      E=1.+1.2*A*an/(1.8+X)/(3.0+0.5*X)
      H=1.-0.67*A*an/(3.0+0.5*X)/(3.0+0.5*X)

      Prinv = A*D/H*c1/(c1t*(1-c2))

      cmue1=0.4/(1.8+X)*H/D*0.33

      do i=1,10         ! Loop for iteration to calculate as 

      as = 1. + an*Prinv + 1/cmue1

      ww=2./3.*0.8/(E*(1.8+X)-0.16/(1.8+X)*H/D*as)
      xxx   = 0.4/(1.8+X)*H/D*ww
      err = abs(cmue1-xxx)
      cmue1 = xxx
      end do


      ww=2./3.*0.8/(E*(1.8+X)-0.16/(1.8+X)*H/D*as)

      cmue1 = 0.4/(1.8+X)*H/D*ww
      cmue2 = A/(3.0+0.5*X)*ww

      cmue1 = 1./cm0**3 * cmue1
      cmue2 = 1./cm0**3 * cmue2 


      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cmue_ca() - Canuto et al. [2000] version A non-equilibrium
!                       stability functions.
!
! !INTERFACE:
      subroutine cmue_ca(as,an,cmue1,cmue2)
      implicit none
 
!
! !USES:
!
! !INPUT PARAMETERS:
      double precision as,an
!
! !OUTPUT PARAMETERS:
      double precision cmue1,cmue2
!
! !DESCRIPTION:
!     This subroutine computes Canuto et al. [2000] version A non-equilibrium
!     stability functions.
!
!
! !BUGS:
!
! !SEE ALSO:
!     cmue_kc(),cmue_bb(),cmue_cb(),cmue_kcQe(),cmue_bbQe(),
!     cmue_caQe(),cmue_cbQe(),cmue_cn(),cmue_ma(),cmue_sg(),cmue_rf()
! 
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  18Nov99  Introduced for GOTM2.0.0 
!
! !LOCAL VARIABLES:
      double precision tn,ts,sm,sh,cmy1,cmy2                               
      double precision L1,L2,L3,L4,L5,L6,L7,L8
      double precision s0,s1,s2,s4,s5,s6
      double precision d0,d1,d2,d3,d4,d5
      double precision d 
      double precision tnmin,tnmax,tsmax,cm0 
!
!EOP
!-------------------------------------------------------------------------
!BOC


      parameter (L1 =  0.1070) 
      parameter (L2 =  0.0032) 
      parameter (L3 =  0.0864) 
      parameter (L4 =  0.1200) 
      parameter (L5 = 11.9000) 
      parameter (L6 =  0.4000) 
      parameter (L7 =  0.0000) 
      parameter (L8 =  0.4800) 
      cm0 = 0.5270

      ts = 4./cm0**6 * as
      tn = 4./cm0**6 * an

      s0 = 1.5*L1*L5*L5
      s1 = -L4*(L6+L7)+2.*L4*L5*(L1-1./3.*L2-L3)+1.5*L1*L5*L8
      s2 = -0.375*L1*(L6*L6-L7*L7)
      s4 = 2.*L5
      s5 = 2.*L4
      s6 = 2./3.*L5*(3.*L3*L3-L2*L2)-0.5*L5*L1*(3.*L3-L2)
     &     +0.75*L1*(L6-L7)

      d0 = 3.*L5*L5
      d1 = L5*(7.*L4+3.*L8)
      d2 = L5*L5*(3.*L3*L3-L2*L2)-0.75*(L6*L6-L7*L7)
      d3 = L4*(4.*L4+3.*L8)
      d4 = L4*(L2*L6-3.*L3*L7-L5*(L2*L2-L3*L3))+L5*L8*(3.*L3*L3-L2*L2)
      d5 = 0.25*(L2*L2-3*L3*L3)*(L6*L6-L7*L7)

      tnmin = -12.27
      
      if (tn.lt.tnmin) tn = tnmin

      tsmax = (d0+d1*tn+d3*tn*tn)/(d2+d4*tn)

      if (ts.gt.tsmax) ts = tsmax

      d = d0 + d1*tn + d2*ts + d3*tn*tn + d4*tn*ts + d5*ts*ts

      sm = (s0 + s1*tn + s2*ts) / d
      sh = (s4 + s5*tn + s6*ts) / d

      cmue1 = 2./cm0**3 * sm
      cmue2 = 2./cm0**3 * sh


      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cmue_caqe() - Equilibrium version of the Canuto et al. [2000]  
!                         version A stability functions.
!
! !INTERFACE:
      subroutine cmue_caqe(an,cmue1,cmue2)
      implicit none
 
!
! !USES:
!
! !INPUT PARAMETERS:
      double precision an
!
! !OUTPUT PARAMETERS:
      double precision cmue1,cmue2
!
! !DESCRIPTION:
!     This subroutine computes Equilibrium version of the Canuto et al. [2000]
!     version A stability functions.
!
! !BUGS:
!
! !SEE ALSO:
!     cmue_kc(),cmue_bb(),cmue_ca(),cmue_cb(),cmue_kcQe(),cmue_bbQe(),
!     cmue_cbQe(),cmue_cn(),cmue_ma(),cmue_sg(),cmue_rf()
! 
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  18Nov99  Introduced for GOTM2.0.0 
!
! !LOCAL VARIABLES:
      double precision tn,ts,sm,sh,cmy1,cmy2  
      double precision L1,L2,L3,L4,L5,L6,L7,L8
      double precision s0,s1,s2,s4,s5,s6
      double precision d0,d1,d2,d3,d4,d5
      double precision d 
      double precision tnmin,tnmax,tsmax,cm0
      double precision PP,QQ
!
!EOP
!-------------------------------------------------------------------------
!BOC



      parameter (L1 =  0.1070) 
      parameter (L2 =  0.0032) 
      parameter (L3 =  0.0864) 
      parameter (L4 =  0.1200) 
      parameter (L5 = 11.9000) 
      parameter (L6 =  0.4000) 
      parameter (L7 =  0.0000) 
      parameter (L8 =  0.4800) 
      cm0 = 0.5270

      tn = 4./cm0**6 * an

      s0 = 1.5*L1*L5*L5
      s1 = -L4*(L6+L7)+2.*L4*L5*(L1-1./3.*L2-L3)+1.5*L1*L5*L8
      s2 = -0.375*L1*(L6*L6-L7*L7)
      s4 = 2.*L5
      s5 = 2.*L4
      s6 = 2./3.*L5*(3.*L3*L3-L2*L2)-0.5*L5*L1*(3.*L3-L2)
     &     +0.75*L1*(L6-L7)

      d0 = 3.*L5*L5
      d1 = L5*(7.*L4+3.*L8)
      d2 = L5*L5*(3.*L3*L3-L2*L2)-0.75*(L6*L6-L7*L7)
      d3 = L4*(4.*L4+3.*L8)
      d4 = L4*(L2*L6-3.*L3*L7-L5*(L2*L2-L3*L3))+L5*L8*(3.*L3*L3-L2*L2)
      d5 = 0.25*(L2*L2-3*L3*L3)*(L6*L6-L7*L7)


      tnmin = -12.27
      
      if (tn.lt.tnmin) tn = tnmin

      PP=(s0+(s1-s6)*tn-2.*(d2+d4*tn))/(s2-2.*d5)
      QQ=-(2.*(d0+d1*tn+d3*tn*tn)+(s4+s5*tn)*tn)/(s2-2.*d5)

      ts=-0.5*PP-sqrt(PP**2/4.-QQ)

      tsmax = (d0+d1*tn+d3*tn*tn)/(d2+d4*tn)

      d = d0 + d1*tn + d2*ts + d3*tn*tn + d4*tn*ts + d5*ts*ts

      sm = (s0 + s1*tn + s2*ts) / d
      sh = (s4 + s5*tn + s6*ts) / d

      cmue1 = 2./cm0**3 * sm
      cmue2 = 2./cm0**3 * sh


      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cmue_cb() - Canuto et al. [2000] version B non-equilibrium
!                       stability functions.
!
! !INTERFACE:
      subroutine cmue_cb(as,an,cmue1,cmue2)
      implicit none
 
!
! !USES:
!
! !INPUT PARAMETERS:
      double precision as,an
!
! !OUTPUT PARAMETERS:
      double precision cmue1,cmue2
!
! !DESCRIPTION:
!     This subroutine computes Canuto et al. [2000] version B non-equilibrium 
!     stability functions.
!
! !BUGS:
!
! !SEE ALSO:
!     cmue_kc(),cmue_bb(),cmue_ca(),cmue_kcQe(),cmue_bbQe(),
!     cmue_caQe(),cmue_cbQe(),cmue_cn(),cmue_ma(),cmue_sg(),cmue_rf()
! 
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  18Nov99  Introduced for GOTM2.0.0 
!
! !LOCAL VARIABLES:
      double precision tn,ts,sm,sh,cmy1,cmy2                               
      double precision L1,L2,L3,L4,L5,L6,L7,L8
      double precision s0,s1,s2,s4,s5,s6
      double precision d0,d1,d2,d3,d4,d5
      double precision d 
      double precision tnmin,tnmax,tsmax,cm0 
!
!EOP
!-------------------------------------------------------------------------
!BOC

      parameter (L1 =  0.1270) 
      parameter (L2 =  0.00336) 
      parameter (L3 =  0.0906) 
      parameter (L4 =  0.1010) 
      parameter (L5 = 11.2000) 
      parameter (L6 =  0.4000) 
      parameter (L7 =  0.0000) 
      parameter (L8 =  0.3180) 
      cm0 = 0.5540

      ts = 4./cm0**6 * as
      tn = 4./cm0**6 * an

      s0 = 1.5*L1*L5*L5
      s1 = -L4*(L6+L7)+2.*L4*L5*(L1-1./3.*L2-L3)+1.5*L1*L5*L8
      s2 = -0.375*L1*(L6*L6-L7*L7)
      s4 = 2.*L5
      s5 = 2.*L4
      s6 = 2./3.*L5*(3.*L3*L3-L2*L2)-0.5*L5*L1*(3.*L3-L2)
     &     +0.75*L1*(L6-L7)

      d0 = 3.*L5*L5
      d1 = L5*(7.*L4+3.*L8)
      d2 = L5*L5*(3.*L3*L3-L2*L2)-0.75*(L6*L6-L7*L7)
      d3 = L4*(4.*L4+3.*L8)
      d4 = L4*(L2*L6-3.*L3*L7-L5*(L2*L2-L3*L3))+L5*L8*(3.*L3*L3-L2*L2)
      d5 = 0.25*(L2*L2-3*L3*L3)*(L6*L6-L7*L7)

      tnmin = -12.27
      
      if (tn.lt.tnmin) tn = tnmin

      tsmax = (d0+d1*tn+d3*tn*tn)/(d2+d4*tn)

      if (ts.gt.tsmax) ts = tsmax

      d = d0 + d1*tn + d2*ts + d3*tn*tn + d4*tn*ts + d5*ts*ts

      sm = (s0 + s1*tn + s2*ts) / d
      sh = (s4 + s5*tn + s6*ts) / d

      cmue1 = 2./cm0**3 * sm
      cmue2 = 2./cm0**3 * sh

      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cmue_cbqe() - Equilibrium version of the Canuto et al. [2000] 
!                         version B stability functions.
!
! !INTERFACE:
      subroutine cmue_cbqe(an,cmue1,cmue2)
      implicit none
 
!
! !USES:
!
! !INPUT PARAMETERS:
      double precision an
!
! !OUTPUT PARAMETERS:
      double precision cmue1,cmue2
!
! !DESCRIPTION:
!     This subroutine computes Equilibrium version of the Canuto et al. [2000]
!     version B stability functions.
!
! !BUGS:
!
! !SEE ALSO:
!     cmue_kc(),cmue_bb(),cmue_ca(),cmue_cb(),cmue_kcQe(),cmue_bbQe(),
!     cmue_caQe(),cmue_cn(),cmue_ma(),cmue_sg(),cmue_rf()
! 
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  18Nov99  Introduced for GOTM2.0.0 
!
! !LOCAL VARIABLES:
      double precision tn,ts,sm,sh,cmy1,cmy2  
      double precision L1,L2,L3,L4,L5,L6,L7,L8
      double precision s0,s1,s2,s4,s5,s6
      double precision d0,d1,d2,d3,d4,d5
      double precision d 
      double precision tnmin,tnmax,tsmax,cm0
      double precision PP,QQ
!
!EOP
!-------------------------------------------------------------------------
!BOC


      parameter (L1 =  0.1270) 
      parameter (L2 =  0.00336) 
      parameter (L3 =  0.0906) 
      parameter (L4 =  0.1010) 
      parameter (L5 = 11.2000) 
      parameter (L6 =  0.4000) 
      parameter (L7 =  0.0000) 
      parameter (L8 =  0.3180) 
      cm0 = 0.5540


      tn = 4./cm0**6 * an

      s0 = 1.5*L1*L5*L5
      s1 = -L4*(L6+L7)+2.*L4*L5*(L1-1./3.*L2-L3)+1.5*L1*L5*L8
      s2 = -0.375*L1*(L6*L6-L7*L7)
      s4 = 2.*L5
      s5 = 2.*L4
      s6 = 2./3.*L5*(3.*L3*L3-L2*L2)-0.5*L5*L1*(3.*L3-L2)
     &     +0.75*L1*(L6-L7)

      d0 = 3.*L5*L5
      d1 = L5*(7.*L4+3.*L8)
      d2 = L5*L5*(3.*L3*L3-L2*L2)-0.75*(L6*L6-L7*L7)
      d3 = L4*(4.*L4+3.*L8)
      d4 = L4*(L2*L6-3.*L3*L7-L5*(L2*L2-L3*L3))+L5*L8*(3.*L3*L3-L2*L2)
      d5 = 0.25*(L2*L2-3*L3*L3)*(L6*L6-L7*L7)


      tnmin = -12.27
      
      if (tn.lt.tnmin) tn = tnmin

      PP=(s0+(s1-s6)*tn-2.*(d2+d4*tn))/(s2-2.*d5)
      QQ=-(2.*(d0+d1*tn+d3*tn*tn)+(s4+s5*tn)*tn)/(s2-2.*d5)

      ts=-0.5*PP-sqrt(PP**2/4.-QQ)

      tsmax = (d0+d1*tn+d3*tn*tn)/(d2+d4*tn)

      d = d0 + d1*tn + d2*ts + d3*tn*tn + d4*tn*ts + d5*ts*ts

      sm = (s0 + s1*tn + s2*ts) / d
      sh = (s4 + s5*tn + s6*ts) / d
      
      cmue1 = 2./cm0**3 * sm
      cmue2 = 2./cm0**3 * sh


      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cmue_cn() - const. stability functions.
!
! !INTERFACE:
      subroutine cmue_cn(cmue1,cmue2)
      implicit none

      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
!
!  Prandtl0, cmust are read from gotmturb.inp and declared in gotmturb.i
!
! !OUTPUT PARAMETERS:
      double precision cmue1,cmue2
!
! !DESCRIPTION:
!     This subroutine computes constant stability functions.  
!
! !BUGS:
!
! !SEE ALSO:
!
!     cmue_kc(),cmue_ca(),cmue_cb(),cmue_kcqe(),cmue_bbqe(),
!     cmue_caqe(),cmue_cbqe(),cmue_bb(),cmue_ma(),cmue_sg(),cmue_rf()
!
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
!
!EOP
!-------------------------------------------------------------------------
!BOC
       
      cmue1=cmucst
      cmue2=cmucst/Prandtl0 

      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cmue_kc() - Kantha and Clayson [1994] non-equilibrium 
!                       stability functions.
!
! !INTERFACE:
      subroutine cmue_kc(as,an,cmue1,cmue2)
      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
      double precision as,an
!
! !OUTPUT PARAMETERS:
      double precision cmue1,cmue2
!
! !DESCRIPTION:
!     This subroutine computes Mellor and Yamada [1982] stability functions
!     with the extension by Kantha and Clayson [1994] to non-zero values 
!     for $c_2$ and $c_3$.
!
! !BUGS:
!
! !SEE ALSO:
!     cmue_bb(),cmue_ca(),cmue_cb(),cmue_kcQe(),cmue_bbQe(),
!     cmue_caQe(),cmue_cbQe(),cmue_cn(),cmue_ma(),cmue_sg(),cmue_rf()
! 
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  18Nov99  Introduced for GOTM2.0.0 
!
! !LOCAL VARIABLES:
      double precision gm,gh,sm,sh
      double precision c11,c12,c13,c21,c22,c23

!
!EOP
!-------------------------------------------------------------------------
!BOC

      gm=0.5*as             !Transformation to MY notation
      gh=-0.5*an            !Transformation to MY notation
      if (gh.gt.0.029) gh=0.029
      if (gm.gt.0.825-25.0*gh) gm=0.825-25.0*gh
      c11=6*a1*a2*gm
      c12=1-3*a2*b2*(1-c3)*gh-12*a1*a2*gh
      c13=a2
      c21=1+6*a1*a1*gm-9*a1*a2*gh
      c22=-12*a1*a1*gh-9*a1*a2*(1-c2)*gh
      c23=a1*(1-3*c1)
      sm=(c12*c23-c22*c13)/(c12*c21-c22*c11)
      sh=(c21*c13-c11*c23)/(c12*c21-c22*c11)
      cmue1=sqrt(2.)*sm     !Retransformation to GOTM notation
      cmue2=sqrt(2.)*sh     !Retransformation to GOTM notation

     
      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cmue_kcqe() - Kantha and Clayson [19995] stability functions.
!
! !INTERFACE:
      subroutine cmue_kcqe(an,cmue1,cmue2)
      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
      double precision an
!
! !OUTPUT PARAMETERS:
      double precision cmue1,cmue2
!
! !DESCRIPTION:
!     This subroutine computes Kantha and Clayson [19995] stability functions
!     including smoothing for convective conditions as discussed by Burchard
!     and Petersen [1999].
!
!
! !BUGS:
!
! !SEE ALSO:
!     cmue_kc(),cmue_bb(),cmue_ca(),cmue_cb(),cmue_bbQe(),
!     cmue_caQe(),cmue_cbQe(),cmue_cn(),cmue_ma(),cmue_sg(),cmue_rf()
! 
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  18Nov99  Introduced for GOTM2.0.0 
!
! !LOCAL VARIABLES:
      double precision gh,sm,sh
!
!
!EOP
!-------------------------------------------------------------------------
!BOC

      gh=-0.5*an           ! Transformation to MY notation
      qeghmax=1/(a2*(b1+12*a1+3*b2*(1-c3)))
      if (qesmooth) then
         if (gh.gt.qeghcrit)
     &         gh=gh-(gh-qeghcrit)**2/(gh+qeghmax-2*qeghcrit)
      else
         if (gh.gt.qeghmax) gh=qeghmax
      end if
      if (gh.lt.qeghmin)  gh = qeghmin
      sm=1-3*c1-6*a1/b1-3*a2*gh*((b2*(1-c3)-3*a2*(1-c2))*(1-6*a1/b1)
     &   -3*c1*(b2*(1-c3)+6.*a1))
      sm=a1*sm/((1-3*a2*gh*(6*a1+b2*(1-c3)))*(1-9*a1*a2*gh))
      sh=a2*(1-6*a1/b1)/(1-3.*a2*gh*(6*a1+b2*(1-c3)))
      cmue1=sqrt(2.)*sm    ! Retransformation to GOTM notation
      cmue2=sqrt(2.)*sh    ! Retransformation to GOTM notation

      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cmue_ma() - Munk-Anderson stability functions.
!
! !INTERFACE:
      subroutine cmue_ma(cmue1,cmue2,NN,SS)
      implicit none

      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
      double precision NN,SS  
!
!     Prandtl0,cmucst are read from gotmturb.inp and declared in gotmturb.i
!
! !OUTPUT PARAMETERS:
      double precision cmue1,cmue2
!
! !DESCRIPTION:
!     This subroutine computes stability functions 
!     according to Munk \& Anderson [1948]:  
!
!     \begin{equation}
!     \begin{array}{ll}
!     c_{\mu} = cmucst \\
!     c_{\mu}'= \frac{c_{\mu}}{P_r^0}
!     \frac{(1+10 R_i)^{1/2}}{(1+3.33 R_i)^{3/2}}, &  R_i \geq 0,\\
!     c_{\mu}'= c_{\mu}, &  R_i<0,
!     \end{array}
!     \end{equation}
!
! !BUGS:
!
! !SEE ALSO:
!
!     cmue_kc(),cmue_ca(),cmue_cb(),cmue_kcqe(),cmue_bbqe(),
!     cmue_caqe(),cmue_cbqe(),cmue_cn(),cmue_bb(),cmue_sg(),cmue_rf()
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
      double precision Ri,Prandtl
!
!EOP
!-------------------------------------------------------------------------
!BOC
   
      Ri=NN/(SS+1e-8)   ! Gradient Richardson number 
      if (Ri.ge.1e-10) then 
         Prandtl=Prandtl0*(1.+3.33*Ri)**1.5/sqrt(1.+10.0*Ri)
      else 
         Prandtl=Prandtl0
      end if 
   
      cmue1=cmucst
      cmue2=cmucst/Prandtl 

      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cmue_rf() - simple formula stability functions.
!
! !INTERFACE:
      subroutine cmue_rf(cmue1,cmue2,NN,SS,xRf)
      implicit none

      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
      double precision NN,SS  
!
!     Prandtl0, cmucst are read from gotmturb.inp and declared in gotmturb.i
!     cmucst is also read from gotmturb.inp and declared in gotmturb.i
!
! !OUTPUT PARAMETERS:
      double precision cmue1,cmue2,xRf
!
! !DESCRIPTION:
!     This subroutine computes stability functions with some "simple" formulas
!
!     Stab=Isprastab
!     \begin{equation}
!     \begin{array}{ll}
!     c_{\mu} = cmucst \\
!     c_{\mu}'= \frac {1} {Prandtl0} (1-R_f)^{1/2}
!     \end{array}
!     \end{equation}
!     \begin{equation} 
!     1-R_f=(\sqrt{R_i^2+1}-R_i)^2 
!     \end{equation}
!     \begin{equation}
!     R_i=\frac {1} {2 Prandtl0} \frac {N^2} {S^2}
!     \end{equation}
!
!     $R_f$ is limited for supercritically stable stratification $1.8<R_f$.  
!
!     For details, see the GOTM report.
!
! !BUGS:
!
! !SEE ALSO:
!
!     cmue_kc(),cmue_ca(),cmue_cb(),cmue_kcqe(),cmue_bbqe(),
!     cmue_caqe(),cmue_cbqe(),cmue_cn(),cmue_ma(),cmue_sg(),cmue_bb()
!
! 
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
!
! !REVISION HISTORY:
!
!  25Nov98   The GOTM developers  Initial code.
!
! !LOCAL VARIABLES:
      double precision Ri,Prandtl_inv
!
!EOP
!-------------------------------------------------------------------------
!BOC

!   Calculation of xRf=(1-Rf), where Rf is the flux Richardson number
       
      xRf=0.
      Ri=0.5/Prandtl0*NN/(SS+1e-8)
      xRf=(sqrt(Ri*Ri+1)-Ri)**2 
      if (xRf.gt.2.) xRf=2.
      Prandtl_inv=1/Prandtl0*sqrt(xRf)

      if (Prandtl_inv.lt.0.18) Prandtl_inv=0.18
      if (Prandtl_inv.gt.2.0)  Prandtl_inv=2.0
   
      cmue1=cmucst
      cmue2=cmucst*Prandtl_inv

      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cmue_sg() - Simple stability functions with Schumann and Gerz [1995]
!                       parameterization of the Prandtl number.
!
! !INTERFACE:
      subroutine cmue_sg(cmue1,cmue2,NN,SS)
      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
      double precision NN,SS
!
!     Prandtl0,cmucst are read from gotmturb.inp and declared in gotmturb.i
!
! !OUTPUT PARAMETERS:
      double precision cmue1,cmue2
!
! !DESCRIPTION:
!     This subroutine computes Schumann and Gerz [1995] stability functions.
!
!    \begin{equation}
!    c_{\mu}=c_{\mu}^0,\qquad c'_{\mu}=\frac{c_{\mu}^0}{P_r}
!    \end{equation}
!
!    with constant $c_{\mu}^0$. We choose here for the Prandtl number $P_r$ a
!    formulation suggested by {\it Schumann and Gerz} [1995]:
!
!    \begin{equation}
!    P_r=P_r^0\exp\left(-\frac{R_i}{P_r^0R_i^{\infty}}\right)
!    -\frac{R_i}{R_i^{\infty}}
!    \end{equation}
!
!    with the neutral Prandtl number $P_r^0=0.74$ and $R_i^{\infty}=0.25$.
!
!
! !BUGS:
!
! !SEE ALSO:
!     cmue_kc(),cmue_bb(),cmue_ca(),cmue_cb(),cmue_kcQe(),cmue_bbQe(),
!     cmue_caQe(),cmue_cbQe(),cmue_cn(),cmue_ma(),cmue_rf()
! 
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  18Nov99  Introduced for GOTM2.0.0 
!
! !LOCAL VARIABLES:
      double precision Ri,Prandtl
!
!EOP
!-------------------------------------------------------------------------
!BOC

      Ri=NN/(SS+1e-8)   ! Gradient Richardson number
      if (Ri.ge.1e-10) then
         Prandtl=Prandtl0*exp(-Ri/(Prandtl0*0.25))+Ri/0.25
      else
         Prandtl=Prandtl0
      end if

      cmue1=cmucst
      cmue2=cmucst/Prandtl

     
      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: InternalWave() - Internal wave and shear instability  
!                            parameterisation according to
!                            Large et al. [1994]. 
!
! !INTERFACE:
      subroutine InternalWave(Nmx,num,nuh,NN,SS,k) 

      implicit none 
  
      include 'gotmturb.i' 
!
! !USES:
!
! !INPUT PARAMETERS:
      integer Nmx
      double precision NN(0:Nmx),SS(0:Nmx),k(0:Nmx)
!     klimiw is read from gotmturb.inp and declared in gotmturb.i
!
! !OUTPUT PARAMETERS:
      double precision num(0:Nmx),nuh(0:Nmx)
!
! !DESCRIPTION:
!    Imposes eddy viscosity and diffisivity characteristic 
!    of internal wave activity and shear instability when there is extinction 
!    of turbulence as suggested by Large et al. [1994]. 
!    In this case, these new num and nuh 
!    are used instead of those computed with the model.
!
!    When k is small (extinction of turbulence, diagnosed by $k<klimiw$), 
!    $\nu_t$ and $\nu'_t$ are set to empirical values typical 
!    in the presence of internal wave activity (IW) and shear 
!    instability (SI). 
!    {\large
!    \begin{equation}
!    \nu_t=(\nu_t)^{IW}+(\nu_t)^{SI}, \quad
!    \nu_t'=(\nu_t')^{IW}+(\nu'_t)^{SI}
!    \end{equation}
!    \vfill
!    \begin{equation}
!    (\nu_t)^{IW}=10^{-4}, \quad         
!    (\nu'_t)^{IW}=5 10^{-5}
!    \end{equation} 
!    \vfill
!    \begin{eqnarray}
!    (\nu_t)^{SI}=(\nu_t')^{SI}=0, & R_i>0.7 \\
!    (\nu_t)^{SI}=(\nu_t')^{SI}=5 10^{-3} \left[1-\left(\frac {R_i} 
!    {0.7}\right)^2\right]^3, & 0<R_i<0.7 \\
!    (\nu_t)^{SI}= (\nu_t')^{SI}=5 10^{-3}, & R_i < 0
!    \end{eqnarray}
! 
!    The unit of all diffusivities is $m^2 s^{-1}$
! 
! !BUGS:
!
! !SEE ALSO:
!     cmue_qe(), stabilityfunctions() 
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  25Nov98   The GOTM developers  Initial code.
!  14Feb00   Initial damping deleted, Iwmodel test included.

!
! !LOCAL VARIABLES:
      double precision rich(0:MaxN)	 
      double precision rich2,pot,x      
      integer i
!
!EOP
!-------------------------------------------------------------------------
!BOC
!

      if (Iwmodel.eq.2) then
         rich2 = rich_cr*rich_cr
         do i=1,Nmx-1 
            if (k(i).le.klimiw) then
               rich(i)=NN(i)/(SS(i)+1.e-10)
               if (rich(i).lt.rich_cr) then
                  if (rich(i).gt.0) then
                     pot=1-rich(i)*rich(i)/rich2 
                     x=numshear*pot*pot*pot
                     num(i)=numiw+x 
                     nuh(i)=nuhiw+x  
                  else
                     num(i)=numiw+numshear
                     nuh(i)=nuhiw+numshear
                  endif          
               else
                  num(i)=numiw
                  nuh(i)=nuhiw
               endif
            endif   
         end do 
      end if 

       return 
       end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Ispralength() - length scale from ISPRAMIX
!
! !INTERFACE:
      subroutine Ispralength(Nmx,k,ko,eps,L,NN,SS,h,
     &                     u_taus,u_taub,depth,xRf)
      implicit none

      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
      integer Nmx
      double precision k(0:Nmx),ko(0:Nmx),eps(0:Nmx),
     &                 NN(0:Nmx),SS(0:Nmx),
     &                 h(0:Nmx),u_taus,u_taub,depth,xRf(0:Nmx)

!
! !OUTPUT PARAMETERS:
      double precision L(0:Nmx)  
!
! !DESCRIPTION:
!
!     This subroutine calculates the 
!     lengthscale used in the ISPRAMIX model of JRC, Ispra, Italy
!     (Eifler and Schrimpf [1992] and Demirov et al. [1998]).
!     $L$ in both mixed layers is obtained from a {\it Blackadar} [1962]
!     type formula:
!     \begin{equation}\label {Lmixed}
!     L=\frac {\kappa \tilde z} {1+\frac {\kappa \tilde z} {c_2 \cdot h_m}}
!     (1-R_f)^e
!     \end{equation}
!     where $\tilde z$
!     is the distance from the interface (surface or bottom). The
!     fraction in (\ref{Lmixed})
!     predicts an approximation to a linear behavior of $L$ near boundaries 
!     and a value proportional to the thickness of the mixed
!     layer far from the interface, $L=c_2 h_m$, where $c_2=0.065$
!     is estimated from experimental data as discussed in
!     {\it Eifler and Schrimpf} [1992].
!     The factor $(1-R_f)$, with the flux Richardson
!     number $R_f=-B/P$, accounts for the effect
!     of stratification on the length scale.
!     The parameter $e$ is here a tuning parameter
!     (pers.\ comm.\ Walter Eifler, JRC, Ispra, Italy)
!     which is usually set to $e=1$.
!
! !BUGS:
!
! !SEE ALSO: 
!     potentialml(), algebraiclength()   
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  25Nov98   The GOTM developers  Initial code.
!  14Nov00   Changes due to separation of turbulent part
!
! !LOCAL VARIABLES:
      integer i,SLind,BLind,Index,Index2
      double precision hms,hmb,db,ds
      double precision kml,c2_i,cd,c3_i
!
!EOP
!-------------------------------------------------------------------------
!BOC
      kml   = 1.e-5
      c2_i  = 0.065
      cd    = 0.125

!  Calculation of surface mixed layer depth 
      hms=0.
      SLind=1
      do i=Nmx,1,-1
         hms=hms+h(i)
         if (k(i).le.kml) then 
            SLind=i
            goto 500
         end if   
      end do
  500  continue
!    Calculation of bottom mixed layer depth
      hmb=0.
      BLind=Nmx
      do i=1,Nmx
         hmb=hmb+h(i)
         if (k(i).le.kml) then
            BLind=i
            goto 501 
         end if
      end do
  501  Continue

! If there is no point where k < kml, the water column is assumed to be mixed.
      if (BLind.gt.SLind) then 
         hms=0.5*depth 
         hmb=0.5*depth
         BLind=int(Nmx/2)
         SLind=int(Nmx/2)+1 
      endif

! Calculation of mixing length in bottom layer 
      db=0.
      do i=1,BLind 
         db=db+h(i)
         L(i)=kappa*db/(1.+kappa*db/(c2_i*hmb+L_min))*xRf(i)**3 
         if (L(i).lt.L_min) L(i)=L_min
      end do 

! Calculation of mixing length in surface layer
      ds=h(Nmx)
      do i=Nmx-1,SLind,-1
         ds=ds+h(i)
         L(i)=kappa*ds/(1.+kappa*ds/(c2_i*hms+L_min))*xRf(i)**3
         if (L(i).lt.L_min) L(i)=L_min
      end do

! Calculation of mixing length in the intermediate region
       
      c3_i=L(SLind)*sqrt(NN(SLind)/k(SLind))
      if (c3_i.lt.1e-10) c3_i=0.
      Index=Slind-1
      do i=SLind-1,BLind+1,-1
         if (NN(i).le.0.) then
            L(i)=L_min
         else
            L(i)=max(c3_i*sqrt(k(i)/NN(i)),L_min)
            if (L(i).gt.L(SLind)) L(i)=L(SLind) 
         endif
         if (L(i).eq.L_min) then
            Index=i 
            goto 503
         end if
      end do
 503  continue
      c3_i=L(BLind)*sqrt(NN(BLind)/k(BLind))
      if (c3_i.lt.1e-10) c3_i=0.
      Index2=BLind+1  
      do i=BLind+1,Index
         if (NN(i).le.0.) then
            L(i)=L_min
         else
            L(i)=max(c3_i*sqrt(k(i)/NN(i)),L_min)
            if(L(i).gt.L(BLind)) L(i)=L(BLind) 
         endif
         if (L(i).eq.L_min) then
            Index2=i 
            goto 504
         end if
      end do 
  504  continue 
      do i=Index2+1,Index-1
         L(i)=L_min
      end do
     
      return
      end
!EOC

!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: KolPran() - calculates turbulent eddy viscosities
!
! !INTERFACE:
      subroutine KolPran(Nmx,cmue1,cmue2,k,L,num,nuh,u_taus,u_taub)

      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
      integer Nmx 
      double precision cmue1(0:Nmx),cmue2(0:Nmx),k(0:Nmx),L(0:Nmx),
     &                 u_taus,u_taub
!
!     z0b and z0s are read and declared in gotmturb.i 
! !OUTPUT PARAMETERS:
      double precision num(0:Nmx),nuh(0:Nmx)
!     
! !DESCRIPTION:
!     Eddy viscosity/diffusivity are calculated by means of the relation of 
!     Kolmogorov and Prandtl from the computed values of k, L and 
!     stability functions. 
!     \begin{equation}
!     \nu_t = c_{\mu} \sqrt{k}L;\quad \nu'_t = c'_{\mu} \sqrt{k}L,
!     \end{equation}
!
! !BUGS:
!
! !SEE ALSO:
!     stabilityfunctions(),tke(),lengthscale()
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
      integer i
!
!EOP
!-------------------------------------------------------------------------
!BOC
      do i=1,Nmx-1
         num(i)=cmue1(i)*sqrt(k(i))*L(i)
         nuh(i)=cmue2(i)*sqrt(k(i))*L(i)
      end do

      num(0  )=kappa*u_taub*z0b
      num(Nmx)=kappa*u_taus*z0s 
      nuh(0  )=kappa*u_taub*z0b
      nuh(Nmx)=kappa*u_taus*z0s 
 
      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: LengthScaleEq() - use an diff. eq. to get the lengthscale.
!
! !INTERFACE:
      subroutine LengthScaleEq(Nmx,dt,k,ko,eps,L,num,nuh,NN,h,
     &                    cmue1,cmue2,P,B,u_taus,u_taub,depth)
      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!     Tridiagonal()
!
! !INPUT PARAMETERS:
      integer Nmx
      double precision k(0:Nmx),ko(0:Nmx),
     &                 num(0:Nmx),NN(0:Nmx),h(0:Nmx),
     &                 cmue1(0:Nmx),cmue2(0:Nmx),P(0:Nmx),B(0:Nmx),
     &                 nuh(0:Nmx),u_taus,u_taub,depth,dt 
!
! !OUTPUT PARAMETERS:
      double precision eps(0:Nmx),L(0:Nmx)
!
! !DESCRIPTION:
!     This subroutine calculates the lengthscale equation according to
!     Mellor and Yamada [1982]:
!
!     \begin{equation}\label{kL_eq}
!     \partial_t (kL) - \partial_z\left(\nu_L\partial_z (kL)\right) =
!     L \left(c_{L1}P + c_{L3}B
!     -  \left(1 + E_2\left(\frac{L}{L_z}\right)^2\right)\varepsilon \right).
!     \end{equation}
!
!     At the end of the subroutine, the Galperin et al. [1988] length
!     limitation is applied and the disispation rate calculated. 
!
! !BUGS:
!
! !SEE ALSO:
!     algebraiclength(), ispralength(), potentialml(), dissipationeq() 
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  25Nov98   The GOTM developers  Initial code.
!  14Nov00   Changes due to separation of turbulent part
!
! !LOCAL VARIABLES:
      double precision avh(0:MaxN),au(0:MaxN),bu(0:MaxN),
     &                 cu(0:MaxN),du(0:MaxN),q2l(0:MaxN),
     &                 q3(0:MaxN),phi_minus(0:MaxN),
     &                 phi_plus(0:MaxN),Lz(0:MaxN),
     &                 ds,db,Prod,Buoy,
     &                 Diss,Lcrit 
      integer i 
!
!EOP
!-------------------------------------------------------------------------
!BOC
      do i=1,Nmx-1
         q2l(i)=2.*ko(i)*L(i)
         q3 (i)=sqrt(8.*k(i)*k(i)*k(i))
      end do
 
      do i=1,Nmx
         avh(i)=0.5*(num(i-1)+num(i))
      end do
     
      db=0.0    ! Diagnostic Length Scale  
      ds=0.0 
      do i=1,Nmx-1   
         db=db+h(i) 
         ds=depth-db  
         if (MYLength.eq.1) 
     &      Lz(i)=kappa*(ds+z0s)*(db+z0b)/(ds+z0s+db+z0b) ! Parabola shape  
         if (MYLength.eq.2) 
     &      Lz(i)=kappa*min(ds+z0s,db+z0b)                ! Triangle shape  
         if (MYlength.eq.3) 
     &      Lz(i)=kappa*(ds+z0s)                          ! For infinite depth
      end do
       
      do i=1,Nmx-1
         Prod=e1*L(i)*P(i)
         Buoy=e3*L(i)*B(i) ! See Burchard 2000 for the value of E3
         Diss=-q3(i)/b1*
     &         (1.+e2*(L(i)/Lz(i))*(L(i)/Lz(i)))
         if (Prod+Buoy.gt.0) then
            phi_plus(i)=Prod+Buoy
            phi_minus(i)=-Diss
         else
            phi_plus(i)=Prod
            phi_minus(i)=-Buoy-Diss
         end if
      end do
 
      do i=1,Nmx-1
         au(i)=-2.*dt*avh(i)/(h(i)+h(i+1))/h(i)
         cu(i)=-2.*dt*avh(i+1)/(h(i)+h(i+1))/h(i+1)
         bu(i)=1.-au(i)-cu(i)+phi_minus(i)*dt/q2l(i)
         du(i)=(1+phi_plus(i)*dt/q2l(i))*q2l(i)
      end do
 
      cu(0)=0
      bu(0)=1.
      du(0)=2.*k(0)*kappa*z0b 
 
      bu(Nmx)=1.
      au(Nmx)=0
      du(Nmx)=2.*k(Nmx)*kappa*z0s 
 
      call Tridiagonalx(Nmx,0,Nmx,au,bu,cu,du,q2l)
 
      do i=0,Nmx
        L(i)=q2l(i)/(2.*k(i))
        if ((NN(i).gt.0).and.(lengthlim)) then
          Lcrit=sqrt(2*galp*galp*k(i)/NN(i)) 
          if (L(i).gt.Lcrit) L(i)=Lcrit  
        end if 
        if (L(i).lt.L_min) L(i)=L_min
        eps(i)=cde*sqrt(k(i)*k(i)*k(i))/L(i)
        if (eps(i).lt.epsmin) eps(i)=epsmin
      end do
 
      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: LengthScale() - a routine to get the lengthscale 
!
! !INTERFACE:
       subroutine LengthScale(Nmx,dt,k,ko,eps,L,num,nuh,NN,SS,h,
     &                     cmue1,cmue2,P,B,u_taus,u_taub,depth,
     &                     xRf)

      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!     DissipationEq(), LengthScaleEq(), AlgebraicLength(), PotentialMl()
!
! !INPUT PARAMETERS:
      integer Nmx
      double precision k(0:Nmx),ko(0:Nmx),
     &                 num(0:Nmx),NN(0:Nmx),SS(0:Nmx),h(0:Nmx),
     &                 depth,
     &                 cmue1(0:Nmx),cmue2(0:Nmx),P(0:Nmx),B(0:Nmx),
     &                 nuh(0:Nmx),u_taus,u_taub,xRf(0:Nmx),dt
!
! !OUTPUT PARAMETERS:
      double precision eps(0:Nmx),L(0:Nmx)
!
! !DESCRIPTION:
!     Calls different subroutines that calculate the lengthscale $L$  
!     and the dissipation rate $\epsilon$ with different methods.
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
!  14Nov00   Changes due to separation of turbulent part
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC

      if (LengthMeth.eq.DissEq)
     &   call DissipationEq(Nmx,dt,k,ko,eps,L,num,nuh,NN,h,
     &                    cmue1,cmue2,P,B,u_taus,u_taub)
 
      if (LengthMeth.eq.LengthEq)
     &   call LengthScaleEq(Nmx,dt,k,ko,eps,L,num,nuh,NN,h,
     &                    cmue1,cmue2,P,B,u_taus,u_taub,depth)
!
!  Bougeault and Andre (1986)
!
      if (LengthMeth.eq.BougeaultAndre)  
     &   call  potentialml(Nmx,k,ko,eps,L,h,u_taus,u_taub,depth,NN)
!
      if ((LengthMeth.ne.DissEq).and.(LengthMeth.ne.LengthEq).and.
     &    (LengthMeth.ne.BougeaultAndre)) 
     &       call algebraiclength(Nmx,k,ko,eps,L,NN,SS,h,
     &                            u_taus,u_taub,depth,xRf)

 
      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: PotentialML - length sacle using 2 master lenghth scales
!
! !INTERFACE:
      subroutine potentialml(Nmx,k,ko,eps,L,h,u_taus,u_taub,depth,NN)

      implicit none
 
      include 'gotmturb.i'
!
! !USES:

!
! !INPUT PARAMETERS:
      integer Nmx
      double precision k(0:Nmx),ko(0:Nmx),h(0:Nmx),
     &                 u_taus,u_taub,depth,NN(0:Nmx)
!
! !OUTPUT PARAMETERS:
      double precision L(0:Nmx),eps(0:Nmx)
!
! !DESCRIPTION:
!     Computes the length scale by defining two master 
!     length scales $l_u$ and $l_d$
!     \begin{equation}
!     \begin{array}{l}
!     \int_{z_0}^{z_0+l_u(z_0)} [b(z_0)-b(z)] dz =k(z_0) \\
!     \int_{z_0-l_d(z_0)}^{z_0} [b(z)-b(z_0)] dz =k(z_0)
!     \end{array}
!     \end{equation}
! 
!      From $l_u$ and $l_d$ two length scales are defined $l_k$ 
!      (characteristic mixing length)
!      and $l_\epsilon$ (characteristic dissipation length):
!      \begin{equation}
!      \begin{array}{l}
!      l_k(z_0)= min[l_d(z_0),l_u(z_0)] \\
!      l_{\epsilon}(z_0)={[l_d(z_0)l_u(z_0)]}^{1/2}
!      \end{array}
!      \end{equation}
! 
!      $l_k$ is used in kolpran() to compute eddy viscosity/difussivity  
!      (is transported as L()). $l_{\epsilon}$ is ed to compute $\epsilon$:
!      \begin{equation}
!      \epsilon=C_{\epsilon}k^{3/2}l_{\epsilon}^{-1}, with C_{\epsilon}=0.7
!      \end{equation}
 
! !BUGS:
!
! !SEE ALSO:
!     algebraiclength(), ispralength(), dissipationeq(), lengthscaleeq() 
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  25Nov98   The GOTM developers  Initial code.
!  14Feb00   Completely rewritten by Manuel Ruiz Villarreal, Henrique Coelho,
!            and Hans Burchard 
!
! !LOCAL VARIABLES:
!
      integer i,j,jj
      double precision ds(0:MaxN),db(0:MaxN)
      double precision lu(0:MaxN),ld(0:MaxN)
      double precision lk(0:MaxN),leps(0:MaxN)
      double precision Lcrit,NNmin,z,buoydiff,integral,ceps
!EOP
!-------------------------------------------------------------------------
!BOC

      parameter(NNmin=1.e-8)

      db(0)=0.
      ds(Nmx)=0.
      do i=1,Nmx-1
         db(i)=db(i-1)+h(i)      ! distance of intercace i from bottom 
         ds(i)=depth-db(i)       ! distance of intercace i from surface 
      end do
!
! Calculation of lu and ld by solving the integral equation following 
! Gaspar (1990). Some other approximations of the integral equation 
! are possible.
!
! Computation of lupward
!

         do i=1,Nmx-1
            lu(i)=0.
            integral=0.
            buoydiff=0.
            do j=i+1,Nmx
                buoydiff=buoydiff+NN(j-1)*0.5*(h(j)+h(j-1))
                integral=integral+buoydiff*h(j)
                if (integral.ge.ko(i)) then
	           if(j.ne.Nmx) then
                      if(j.ne.i+1) then
                         lu(i)=lu(i)-(integral-ko(i))/buoydiff
                      else 
!      To avoid lu(i) from becoming too large if NN(i) is too small
	                 if(NN(i).gt.NNmin) then
	                    lu(i)=sqrt(2.)*sqrt(ko(i))/sqrt(NN(i))
                         else
	                    lu(i)=h(i)
                         end if
                      end if 
                      goto 600
                   end if
                end if
                lu(i)=lu(i)+h(j)
             end do 
600          continue
!     Implicitely done in the do loop: if (lu(i).gt.ds(i)) lu(i)=ds(i) 
!     lu limited by distance to surface 
          end do

!      Computation of ldownward
          do i=Nmx-1,1,-1
             ld(i)=0.
	     integral=0.
             buoydiff=0.
             do j=i-1,1,-1 
                buoydiff=buoydiff+NN(j)*0.5*(h(j+1)+h(j))
                integral=integral-buoydiff*h(j)
                if (integral.ge.ko(i)) then
                   if(j.ne.0) then
                      if(j.ne.i-1) then
                         ld(i)=ld(i)-(integral-ko(i))/buoydiff
                      else
!      To avoid ld(i) from becoming too large if NN(i) is too small
                         if(NN(i).gt.NNmin) then
                            ld(i)=sqrt(2.)*sqrt(ko(i))/sqrt(NN(i))
                         else
                            ld(i)=h(i)
                         end if
                      end if 
                      goto 610
                   end if
                end if
                ld(i)=ld(i)+h(j)
             end do
610	     continue
!            if (ld(i).gt.db(i)) ld(i)=db(i) !ld limited by distance to bottom
          end do         


!   Calculation of lk and leps, mixing and dissipation lengths
      do i=Nmx-1,1,-1 
!   Suggested by Gaspar:        lk(i)   = min(lu(i),ld(i))
         lk(i)=sqrt(lu(i)*ld(i))
         leps(i) = sqrt(lu(i)*ld(i)) 
      end do

!    We set L=lk because it is the one we use to calculate num and nuh
 
      ceps=0.7
      do i=1,Nmx-1
         L(i)=lk(i)
      end do      
     
      L(0)=kappa*z0b
      L(Nmx)=kappa*z0s
! Gaspar uses null gradient
      do i=0,Nmx
         if ((NN(i).gt.0).and.(lengthlim)) then
            Lcrit=sqrt(2*galp*galp*k(i)/NN(i))
            if (L(i).gt.Lcrit) L(i)=Lcrit
         end if
         if (L(i).lt.L_min) L(i)=L_min
         eps(i)=cde*sqrt(k(i)*k(i)*k(i))/L(i)
         if(eps(i).lt.epsmin) eps(i)=epsmin
      end do  
 
      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: StabilityFunctions() - calculates the stability functions - NOT Finished.
!
! !INTERFACE:
      subroutine StabilityFunctions(Nmx,k,L,SS,NN,
     &            cmue1,cmue2,P,B,num,nuh,xRf)

      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!     Buoyancy(), cmue_cn(), cmue_qe(), function eqstate()
!
! !INPUT PARAMETERS:
      integer Nmx 
      double precision L(0:Nmx),num(0:Nmx),nuh(0:Nmx),k(0:Nmx) 
!
! !OUTPUT PARAMETERS:
      double precision NN(0:Nmx),SS(0:Nmx),cmue1(0:Nmx),cmue2(0:Nmx),
     &                 xRf(0:Nmx),P(0:Nmx),B(0:Nmx)
!
! !DESCRIPTION:
!     This subroutine calculates different parameters needed for the
!     turbulence equations such as
!
!     NN: square of Brunt-Vaisala frequency
!     SS: square of shear frequency    
!     P : shear production of turbulent kinetic energy
!     B : buoyancy production of turbulent kinetic energy
!
!     Furthermore, the stability functions are called here, and the
!     needed non-dimensional parameters as and an calculated. 
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
!  14Nov00   Changes due to separation of turbulent part
!
! !LOCAL VARIABLES:
      integer i
      double precision LLk,dudz,dvdz
      double precision as(0:MaxN),an(0:MaxN)
!
!EOP
!-------------------------------------------------------------------------
!BOC

      do i=Nmx-1,1,-1 
         LLk=L(i)*L(i)/k(i)
         as(i)=LLk*SS(i)
         an(i)=LLk*NN(i)
         if (Iwmodel.ne.1) then
            P(i)=SS(i)*num(i) 
         else 
            P(i)=(SS(i)+alfa*NN(i))*num(i) 
         end if 
         B(i)=-NN(i)*nuh(i) 
 
      if (Stab.eq.KanClay)
     &      call cmue_kc(as(i),an(i),cmue1(i),cmue2(i))
      if (Stab.eq.BurBaum)
     &      call cmue_bb(as(i),an(i),cmue1(i),cmue2(i))
      if (Stab.eq.CanutoA)
     &      call cmue_ca(as(i),an(i),cmue1(i),cmue2(i))
      if (Stab.eq.CanutoB)
     &      call cmue_cb(as(i),an(i),cmue1(i),cmue2(i))
      if (Stab.eq.KanClayQe)
     &      call cmue_kcqe(an(i),cmue1(i),cmue2(i))
      if (Stab.eq.BurBaumQe)
     &      call cmue_bbqe(an(i),cmue1(i),cmue2(i))
      if (Stab.eq.CanutoAQe)
     &      call cmue_caqe(an(i),cmue1(i),cmue2(i))
      if (Stab.eq.CanutoBQe)
     &      call cmue_cbqe(an(i),cmue1(i),cmue2(i))
      if (Stab.eq.Constan)
     &      call cmue_cn(cmue1(i),cmue2(i))
      if (Stab.eq.MunkAnd)
     &      call cmue_ma(cmue1(i),cmue2(i),NN(i),SS(i))
      if (Stab.eq.SchumGerz)
     &      call cmue_sg(cmue1(i),cmue2(i),NN(i),SS(i))
      if (Stab.eq.FluxRich)
     &      call cmue_rf(cmue1(i),cmue2(i),NN(i),SS(i),xRf(i))

       end do

      cmue1(0)=cmue1(1)
      cmue1(Nmx)=cmue1(Nmx-1)
      cmue2(0)=cmue2(1)
      cmue2(Nmx)=cmue2(Nmx-1)
 
      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: TKEalgebraic() - local equilibrium assumption to get TKE  
!
! !INTERFACE:
      subroutine tkealgebraic(Nmx,cmue1,cmue2,k,ko,eps,L,h,NN,SS,
     &                        u_taus,u_taub)
      implicit none

      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
      integer Nmx
      double precision cmue1(0:Nmx),cmue2(0:Nmx),eps(0:Nmx),L(0:Nmx),
     &                 h(0:Nmx),NN(0:Nmx),SS(0:Nmx),
     &                 u_taus,u_taub
!
! !OUTPUT PARAMETERS:
      double precision k(0:Nmx),ko(0:Nmx)
!
! !DESCRIPTION:
!     This subroutine calculates the turbulent kinetic energy based
!     on the local equilibrium assumption
!
!     \begin{equation}
!     P+B-\varepsilon=0.
!     \end{equation}  
!
! !BUGS:
!
! !SEE ALSO:
!     theeq() 
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  25Nov98   The GOTM developers  Initial code.
!  14Nov00   Changes due to separation of turbulent part
!
! !LOCAL VARIABLES:
      integer i
!
!EOP
!-------------------------------------------------------------------------
!BOC

      do i=1,Nmx-1   
         ko(i)=k(i)
         k(i)=L(i)*L(i)/cde*(cmue1(i)*SS(i)-cmue2(i)*NN(i))
      end do 
      k(0  )=u_taub*u_taub/sqrt(cm0*cde) 
      k(Nmx)=u_taus*u_taus/sqrt(cm0*cde)
      do i=0,Nmx   
         if (k(i).lt.k_min) k(i)=k_min 
      end do 
 
      return
     
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: TKEeq() - use a diff. eq. to get the TKE
!
! !INTERFACE:
      subroutine TKEeq(Nmx,dt,num,P,B,eps,L,k,ko,h,u_taus,u_taub)
      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!     Tridiagonal()
!
! !INPUT PARAMETERS:
      integer Nmx 
      double precision num(0:Nmx),P(0:Nmx),B(0:Nmx),eps(0:Nmx),
     &        L(0:Nmx),k(0:Nmx),ko(0:Nmx),h(0:Nmx)
      double precision u_taus,u_taub,dt
!
! !OUTPUT PARAMETERS:
!!      double precision num(0:Nmx),P(0:Nmx),B(0:Nmx),eps(0:Nmx),
!!     &        L(0:Nmx),k(0:Nmx),ko(0:Nmx),h(0:Nmx)
!!      double precision u_taus,u_taub
!
! !DESCRIPTION:
!     This subroutine calculates the turbulent kinetic energy as
!     needed for one- or two-equation models:
!
!     \begin{equation}\label{k_eq}
!     \partial_t k - \partial_z(\nu_k\partial_z k) =  P + B -\varepsilon,
!     \end{equation}
!
!     The diffusion coefficient depends on the type of model (k-epsilon
!     or Mellor-Yamada). 
!
!     As boundary conditions a choice between Dirichlet (fluxcond=.false.)
!     and Neumann no-flux conditions (fluxcond=.true.) has to be made.
!
!     Dirichlet condition:
!
!     \begin{equation}
!     k=\left(\frac{u_*}{c_{\mu}^0}\right)^2. 
!     \end{equation}
!
!     The sink terms are treated quasi-implicitely in oder to guarantee
!     positivity. 
    
!
! !BUGS:
!
! !SEE ALSO:
!     tkealgebraic() 
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  25Nov98   The GOTM developers  Initial code.
!  14Nov00   Changes due to separation of turbulent part
!
! !LOCAL VARIABLES:
      integer i
      double precision avh(0:MaxN),au(0:MaxN),bu(0:MaxN),
     &                 cu(0:MaxN),du(0:MaxN),
     &                 numtke(0:MaxN)  
      double precision pminus(0:MaxN),pplus(0:MaxN),Prod,Buoy,Diss
!
	logical istdebug
!EOP
!-------------------------------------------------------------------------
!BOC
!     +----------------------------------------------------------------+
!     | Special treatment of diffusivity for MY model                  |
!     +----------------------------------------------------------------+
      if (TKEMeth.eq.TKEMY) then
         do i=1,Nmx-1
            numtke(i)=Sl*sqrt(2.*k(i))*L(i)
         end do
      else 
         do i=1,Nmx-1
            numtke(i)=num(i)/sig_k 
         end do
      end if
 
      do i=0,Nmx
          ko(i)=k(i)
      end do
 
      do i=2,Nmx-1
         avh(i)=0.5*(numtke(i-1)+numtke(i))
      end do

      if (fluxcond) then 
         avh(1  )=0
         avh(Nmx)=0
      else 
         if (TKEMeth.eq.TKE_keps) then  
            avh(1  )=u_taub**4*2/(eps(0  )+eps(1    ))
            avh(Nmx)=u_taus**4*2/(eps(Nmx)+eps(Nmx-1))
         else   
            avh(1  )= 0.5*(num(0  )+numtke(1    ))
            avh(Nmx)= 0.5*(num(Nmx)+numtke(Nmx-1))
         end if
      end if
      do i=1,Nmx-1
         Prod=P(i)
         Buoy=B(i)
         Diss=eps(i)
         !if (Prod+Buoy.gt.0) then
         if (Buoy.gt.0) then		!ggu
            pplus(i)=Prod+Buoy
            pminus(i)=Diss
         else
            pplus(i)=Prod
            pminus(i)=Diss-Buoy
         end if
      end do
 
      do i=1,Nmx-1
         au(i)=-2.*dt*avh(i)/(h(i)+h(i+1))/h(i)
         cu(i)=-2.*dt*avh(i+1)/(h(i)+h(i+1))/h(i+1)
         bu(i)=1.-au(i)-cu(i)+pminus(i)*dt/k(i)
         du(i)=(1+pplus(i)*dt/k(i))*k(i)
         du(i)=k(i)+pplus(i)*dt
c	 if( istdebug() ) then		!ggu ggg
c	   write(91,*) i,du(i),du(199)
c	 end if
      end do

c ggu

        call pr_info('in after loop abcd',91,-1,Nmx,30
     +                  ,pplus,pminus,k,au,bu,cu,du) 

	if( istdebug() ) then
	  i= Nmx - 1
	  write(91,*) 'tkeeq extra...'
	  write(91,*) Nmx,i,Nmx - i + 1,dt
	  write(91,*) k(i),pplus(i), k(i)+pplus(i)*dt,du(i)
	end if
        call pr_info('in TKEeq',91,-1,Nmx,30
     +                  ,num,avh,k,eps,L,P,B)                                 

        call pr_info('in TKEeq...',91,-1,Nmx,30
     +                  ,pplus,pminus,k,eps,L,P,B)                                 

      if (fluxcond) then
!     +-------------------------------------------------------------+
!     | No-flux conditions for TKE                                  | 
!     +-------------------------------------------------------------+
        call pr_info('befor tridiag',91,-1,Nmx,30
     +                  ,num,avh,k,au,bu,cu,du)                                 

	if( istdebug() ) then
	  write(91,*) 'tridiag extra...'
	  write(91,*) Nmx,1,Nmx-1
	end if
         call Tridiagonalx(Nmx,1,Nmx-1,au,bu,cu,du,k)
         k(0  ) = u_taub*u_taub/sqrt(cm0*cde) 
         k(Nmx) = u_taus*u_taus/sqrt(cm0*cde) 
      else
!     +-------------------------------------------------------------+
!     | Dirichlet conditions for TKE                                | 
!     +-------------------------------------------------------------+
         cu(0)=0
         bu(0)=1.
         du(0)=u_taub*u_taub/sqrt(cm0*cde)
 
         bu(Nmx)=1.
         au(Nmx)=0
         du(Nmx)=u_taus*u_taus/sqrt(cm0*cde)

         call Tridiagonalx(Nmx,0,Nmx,au,bu,cu,du,k)
      end if 

        call pr_info('after tridiag',91,-1,Nmx,30
     +                  ,num,avh,k,au,bu,cu,du)                                 

      do 104 i=0,Nmx
         if (k(i).lt.k_min) k(i)=k_min
104   continue
222   format (5(F15.10,1x)) 
        call pr_info('after linit',91,-1,Nmx,30
     +                  ,num,avh,k,au,bu,cu,du)                                 

 
      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: TKE() - a wrapper routine to get the TKE
!
! !INTERFACE:
      subroutine TKE(Nmx,dt,cmue1,cmue2,num,P,B,eps,L,k,ko,h,
     &              u_taus,u_taub,NN,SS)
      implicit none
 
      include 'gotmturb.i'
!
! !USES:
!     TKEeq()
!     TKEalgebraic()
!
! !INPUT PARAMETERS:
      integer Nmx 
      double precision cmue1(0:Nmx),cmue2(0:Nmx),num(0:Nmx),P(0:Nmx),
     &        B(0:Nmx),eps(0:Nmx),
     &        L(0:Nmx),h(0:Nmx),
     &        NN(0:Nmx),SS(0:Nmx)
      double precision u_taus,u_taub,dt
!
! !OUTPUT PARAMETERS:
       double precision k(0:Nmx),ko(0:Nmx)
!
! !DESCRIPTION:
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
!  14Nov00   Changes due to separation of turbulent part
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
	call write_common(91)	!ggu	only if istdebug() true

      if ((TKEMeth.eq.TKE_keps).or.(TKEMeth.eq.TKEMY)) 
     &  call TKEeq(Nmx,dt,num,P,B,eps,L,k,ko,h,u_taus,u_taub)

      if (TKEMeth.eq.TKEloceq) 
     &  call TKEalgebraic(Nmx,cmue1,cmue2,k,ko,eps,L,h,NN,SS,
     &                    u_taus,u_taub)
      
      return
      end
!EOC
!-------------------------------------------------------------------------
!             General Ocean Turbulence Model (GOTM)                      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: TridiagonalX() - solves a linear set or eqs. 
!
! !INTERFACE:
      subroutine TridiagonalX(Nmx,fi,lt,au,bu,cu,du,value)
      implicit none

      include 'gotmturb.i'
!
! !USES:
!
! !INPUT PARAMETERS:
      integer Nmx
      integer fi,lt
      double precision au(0:Nmx),bu(0:Nmx),cu(0:Nmx),
     &                 du(0:Nmx)
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
!  14Nov00   Changes due to separation of turbulent part
!
! !LOCAL VARIABLES:
      double precision ru(0:MaxN),qu(0:MaxN)
      integer i
!
!EOP
!-------------------------------------------------------------------------
!BOC

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
