#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenAmmonium
!
! DESCRIPTION
!   Description of the diagenetic ammonium processes in the sediment
!   Details on the equations and the method used to calculate
!   the equilibrium and transient profiles can be found in
!   Ruardij et al., 1995. Neth. J. Sea Res. 33(3/4):453-483
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenAmmoniumDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: K4n, K3n, G2o
  ! The following Benthic-states are used (NOT in fluxes): D1m, K14n, D2m, K24n, &
  ! D7m
  ! The following global vars are modified: dummy
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY, LocalDelta
  ! The following Benthic 1-d global boxvars are modified : M4n, KNH4, jG2K3o, &
  ! jK4N4n
  ! The following Benthic 1-d global boxvars got a value: M14n, M24n
  ! The following Benthic 1-d global boxvars are used: reBTn, reATn, &
  ! irrenh, ETW_Ben, N4n_Ben, Depth_Ben
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global box parametes are used: p_d_tot, p_clDxm, &
  ! p_q10diff, p_qon_nitri
  ! The following global constants are used: RLEN
  ! The following constants are used: LAYERS, &
  ! LAYER1, LAYER2, DIFFUSION, FOR_ALL_LAYERS, POROSITY, &
  ! ADSORPTION, DOUBLE_DEFINE, EXPONENTIAL_TERM, DEFINE, CONSTANT_TERM, &
  ! ZERO_EXPONENTIAL_TERM, LINEAR_TERM, SET_CONTINUITY, FLAG, MASS, &
  ! SET_BOUNDARY, EQUATION, SET_LAYER_INTEGRAL, SET_LAYER_INTEGRAL_UNTIL, &
  ! LAYER3, INPUT_TERM, STANDARD, ADD, INTEGRAL, RFLUX, DERIVATIVE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: K4n, K3n, G2o, D1m, K14n, D2m, K24n, D7m, D2STATE
  use mem, ONLY: ppK4n, ppK3n, ppG2o, ppD1m, ppK14n, ppD2m, &
    ppK24n, ppD7m, dummy, BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, &
    BoxNumberY, NO_BOXES_Y, BoxNumber, BoxNumberXY, LocalDelta, M4n, KNH4, jG2K3o, &
    jK4N4n, M14n, M24n, reBTn, reATn, irrenh, ETW_Ben, N4n_Ben, Depth_Ben, iiBen, &
    iiPel, flux
  use constants, ONLY: LAYERS, LAYER1, LAYER2, &
    DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, DOUBLE_DEFINE, &
    EXPONENTIAL_TERM, DEFINE, CONSTANT_TERM, ZERO_EXPONENTIAL_TERM, LINEAR_TERM, &
    SET_CONTINUITY, FLAG, MASS, SET_BOUNDARY, EQUATION, SET_LAYER_INTEGRAL, &
    SET_LAYER_INTEGRAL_UNTIL, LAYER3, INPUT_TERM, STANDARD, ADD, INTEGRAL, RFLUX, &
    DERIVATIVE
  use mem_Param,  ONLY: p_poro, p_d_tot, p_clDxm, p_q10diff, p_qon_nitri
  use mem_BenAmmonium


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:InitializeSet, &
  ! DefineSet, CompleteSet, CalculateSet, CalculateTau, CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: InitializeSet, DefineSet, CompleteSet, &
    CalculateSet, CalculateTau, CalculateFromSet


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:eTq
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: eTq


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: IntegralExp



!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
!
! !REVISION_HISTORY
!   September 1999 by M. Vichi
!               Commented version
!
!
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: zuBT
  real(RLEN)  :: zuD1
  real(RLEN)  :: alpha
  real(RLEN)  :: diff
  real(RLEN)  :: Tau
  real(RLEN)  :: labda
  real(RLEN)  :: a15
  real(RLEN)  :: cK4n
  real(RLEN)  :: jK4K3n
  real(RLEN)  :: sK4K3
  real(RLEN)  :: cO2
  real(RLEN)  :: r

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! user defined external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, external  :: D3toD1
  integer, external  :: D2toD1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  BoxNumberZ = 1
  DO BoxNumberY=1,NO_BOXES_Y
    DO BoxNumberX=1,NO_BOXES_X
      BoxNumber=D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ)
      BoxNumberXY=D2toD1(BoxNumberX,BoxNumberY)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore-water average concentrations from the state variables
      ! (Diagnostic variables, not used in calculations)
      ! Calculate pore-water oxygen concentration in the oxic layer
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M4n(BoxNumberXY) = K4n(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_p+ &
        1.0D+00)/( D1m(BoxNumberXY))
      M14n(BoxNumberXY) = K14n(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_p+ &
        1.0D+00)/( D2m(BoxNumberXY)- D1m(BoxNumberXY))
      M24n(BoxNumberXY) = K24n(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_p+ &
        1.0D+00)/( p_d_tot- D2m(BoxNumberXY))

      cO2  =   G2o(BoxNumberXY)/ p_poro(BoxNumberXY)/ D1m(BoxNumberXY)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D7.m is the average penetration depth for N-detritus
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      alpha  =   1.0D+00/ max(  p_clDxm,  D7m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate Mineralization m2 --> m3 porewater
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Average in the oxic layer:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      zuBT  =   reBTn(BoxNumberXY)/ p_poro(BoxNumberXY)/ D1m(BoxNumberXY)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Anoxic Mineralization at D1.m, using the exponential distribution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      zuD1 = max( 1.D-20, reATn(BoxNumberXY))/ p_poro(BoxNumberXY)/ IntegralExp( &
        - alpha, p_d_tot- D1m(BoxNumberXY))


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      ! nitrification rate: temperature and oxygen
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* p_poro(BoxNumberXY)* irrenh(BoxNumberXY)* &
        eTq( ETW_Ben(BoxNumberXY), p_q10diff)
      sK4K3 = p_sK4K3* eTq( ETW_Ben(BoxNumberXY), p_q10)* cO2/( &
        cO2+ p_clO2)* M4n(BoxNumberXY)/( M4n(BoxNumberXY)+ p_clM4)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! if availability of carbon for degradation is low , the nitrification &
      ! will be hampered
      ! by the lack of carbon for nitrification bacteria. As proxy for &
      ! the degrdability
      ! of the carbon the respiration per mass of bacteria is calulated.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the exponential terms of the solution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      labda  =   sqrt(  sK4K3/ diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient of the zero order term
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      a15  =   zuBT/ sK4K3

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layer depths
      ! - environmental conditions (diffusion, p_porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      KNH4(BoxNumberXY) = InitializeSet( KNH4(BoxNumberXY), N_layers, &
        N_coeff)
      call DefineSet( KNH4(BoxNumberXY), LAYERS, LAYER1, &
        LAYER2, D1m(BoxNumberXY), D2m(BoxNumberXY))

      call DefineSet( KNH4(BoxNumberXY), DIFFUSION, FOR_ALL_LAYERS, 0, diff, &
        dummy)

      call DefineSet( KNH4(BoxNumberXY), POROSITY, FOR_ALL_LAYERS, 0, &
        p_poro(BoxNumberXY), dummy)

      call DefineSet( KNH4(BoxNumberXY), ADSORPTION, FOR_ALL_LAYERS, 0, p_p, &
        dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! A(z) = a11*exp(labda*z) + a12*exp(-labda*z) + a13*z^2 + a14*z + a15
      ! 2nd layer:
      ! A(z) = a21*exp[-alpha*(z-D1.m)] + a24*z + a25
      ! 3rd layer:
      ! A(z) = a31*exp[-alpha*(z-D2.m)] + a34*z + a35
      !    a34 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call DefineSet( KNH4(BoxNumberXY), DOUBLE_DEFINE, 11, EXPONENTIAL_TERM, &
        labda, sK4K3)

      call DefineSet( KNH4(BoxNumberXY), DOUBLE_DEFINE, 12, EXPONENTIAL_TERM, - &
        labda, sK4K3)

      call DefineSet( KNH4(BoxNumberXY), DEFINE, 15, CONSTANT_TERM, dummy, &
        dummy)

      call DefineSet( KNH4(BoxNumberXY), DEFINE, 21, ZERO_EXPONENTIAL_TERM, - &
        alpha, dummy)

      call DefineSet( KNH4(BoxNumberXY), DEFINE, 24, LINEAR_TERM, dummy, &
        dummy)

      call DefineSet( KNH4(BoxNumberXY), DEFINE, 25, CONSTANT_TERM, dummy, &
        dummy)

      call DefineSet( KNH4(BoxNumberXY), DEFINE, 31, ZERO_EXPONENTIAL_TERM, - &
        alpha, dummy)

      call DefineSet( KNH4(BoxNumberXY), DEFINE, 35, CONSTANT_TERM, dummy, &
        dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call CompleteSet( KNH4(BoxNumberXY), SET_CONTINUITY, FLAG, MASS, &
        dummy, dummy)

      call CompleteSet( KNH4(BoxNumberXY), SET_BOUNDARY, LAYER1, &
        EQUATION, 0.0D+00, N4n_Ben(BoxNumberXY))

      call CompleteSet( KNH4(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER2, &
        LAYER2, dummy, K14n(BoxNumberXY))

      call CompleteSet( KNH4(BoxNumberXY), SET_LAYER_INTEGRAL_UNTIL, LAYER3, &
        LAYER3, p_d_tot, K24n(BoxNumberXY))

      call CompleteSet( KNH4(BoxNumberXY), INPUT_TERM, 15, STANDARD, dummy, &
        a15)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate for the above defined set of boundary conditions
      ! the steady-state profiles and return the vertically integrated
      ! concentration in the first layer.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cK4n = CalculateSet( KNH4(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER1, &
        LAYER1, dummy, 0.0D+00)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the adaptation time to the steady-state profile
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Tau  =   CalculateTau(  sK4K3,  diff,  p_p,  D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate the average value of K4n over the actual time step
      ! (transient value).
      ! This value depends on the adaptation time, the actual time step,
      ! the ''old'' value and the ''equilibrium value''
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cK4n = cK4n+( K4n(BoxNumberXY)- cK4n)* IntegralExp( - LocalDelta/ &
        Tau, 1.0D+00)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Derive the equations for the transient profiles, assuming the same
      ! solution as for the steady-state case and using cK4n as new constraint.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      dummy  =   CalculateSet(  KNH4(BoxNumberXY),  ADD,  0,  0,  dummy,  cK4n)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Start calculation of fluxes:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Next calculation is done in BenNitrogenShiftin which
      ! include all fluxes between layers!
      ! All the nutrient mineralization source term in the anoxic layer
      ! has been added to K14.n in BenBacDynamics
      ! However in the model this layer is subdivided and hence a partition
      ! flux is here calculated according to the exponential distribution.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate Anoxic Mineralization at D2.m
      !
      !  zuD2 = zuD1 * exp( -alpha * (D2.m - D1.m));
      !
      ! K14.n -> K24.n = zuD2 * p_poro * IntegralExp(-alpha, p_d_tot - D2.m);
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate Nitrification flux in the first layer and the related
      ! oxygen consumption flux:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK4K3n = sK4K3* CalculateFromSet( KNH4(BoxNumberXY), INTEGRAL, &
        RFLUX, 0.0D+00, D1m(BoxNumberXY))

      call flux(BoxNumberXY, iiBen, ppK4n, ppK3n, jK4K3n )

      jG2K3o(BoxNumberXY)  =   jK4K3n* p_qon_nitri
      call flux(BoxNumberXY, iiBen, ppG2o, ppG2o, -( jG2K3o(BoxNumberXY)) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimation of the Vertical fluxes from the set of transient solutions:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK4N4n(BoxNumberXY) = max( CalculateFromSet( KNH4(BoxNumberXY), &
        DERIVATIVE, RFLUX, 0.0D+00, 0.0D+00), - 0.5D+00* &
        N4n_Ben(BoxNumberXY)* Depth_Ben(BoxNumberXY))
      call flux(BoxNumberXY, iiBen, ppK4n, ppK4n, -( jK4N4n(BoxNumberXY)) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !  All transport between layers are done in BenNitrogenShifting:
      !
      !  jK14K4n = CalculateFromSet(KNH4, DERIVATIVE, RFLUX, D1.m, 0.0);
      !  jK24K14n = CalculateFromSet(KNH4, DERIVATIVE, RFLUX, D2.m, 0.0);
      !  jK34K24n = CalculateFromSet(KNH4, DERIVATIVE, RFLUX, p_d_tot, 0.0);
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    end DO



  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
