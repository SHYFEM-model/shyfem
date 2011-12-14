#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenNitrate
!
! DESCRIPTION
!   Description of the diagenetic nitrate processes in the sediment
!       Details on the equations and the method used to calculate
!       the equilibrium and transient profiles can be found in
!       Ruardij et al., 1995. Neth. J. Sea Res. 33(3/4):453-483   
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenNitrateDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: K3n, G4n
  ! The following Benthic-states are used (NOT in fluxes): D2m, D6m, D1m, K6r
  ! The following global vars are modified: dummy
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY, LocalDelta
  ! The following Benthic 1-d global boxvars are modified : KNO3, jK3N3n
  ! The following Benthic 1-d global boxvars got a value: M3n
  ! The following Benthic 1-d global boxvars are used: rrATo, irrenh, ETW_Ben, &
  ! KNH4, N3n_Ben
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global box parametes are used: p_d_tot, p_q10diff, &
  ! p_qro, p_qon_dentri
  ! The following global constants are used: RLEN
  ! The following constants are used: GET, &
  ! LABDA_1, LABDA_2, COEFFICIENT, LAYERS, LAYER1, &
  ! DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, DOUBLE_DEFINE, &
  ! ZERO_EXPONENTIAL_TERM, DEFINE, QUADRATIC_TERM, LINEAR_TERM, CONSTANT_TERM, &
  ! EXPONENTIAL_TERM, SET_CONTINUITY, FLAG, MASS, SET_BOUNDARY, &
  ! EQUATION, INPUT_TERM, PARAMETER, START_ADD_TERM, INPUT_SUBTRACT_TERM, &
  ! SET_LAYER_INTEGRAL_UNTIL, LAYER2, ADD, DERIVATIVE, RFLUX, INTEGRAL

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: K3n, G4n, D2m, D6m, D1m, K6r, D2STATE
  use mem, ONLY: ppK3n, ppG4n, ppD2m, ppD6m, ppD1m, ppK6r, &
    dummy, BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
    BoxNumber, BoxNumberXY, LocalDelta, KNO3, jK3N3n, M3n, rrATo, irrenh, ETW_Ben, &
    KNH4, N3n_Ben, iiBen, iiPel, flux
  use constants, ONLY: GET, LABDA_1, LABDA_2, &
    COEFFICIENT, LAYERS, LAYER1, DIFFUSION, FOR_ALL_LAYERS, POROSITY, &
    ADSORPTION, DOUBLE_DEFINE, ZERO_EXPONENTIAL_TERM, DEFINE, QUADRATIC_TERM, &
    LINEAR_TERM, CONSTANT_TERM, EXPONENTIAL_TERM, SET_CONTINUITY, FLAG, MASS, &
    SET_BOUNDARY, EQUATION, INPUT_TERM, PARAMETER, START_ADD_TERM, &
    INPUT_SUBTRACT_TERM, SET_LAYER_INTEGRAL_UNTIL, LAYER2, ADD, DERIVATIVE, RFLUX, &
    INTEGRAL
  use mem_Param,  ONLY: p_poro, p_d_tot, p_q10diff, p_qro, p_qon_dentri
  use mem_BenNitrate


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:GetInfoFromSet, &
  ! InitializeSet, DefineSet, CompleteSet, CalculateSet, CalculateTau, &
  ! CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: GetInfoFromSet, InitializeSet, DefineSet, &
    CompleteSet, CalculateSet, CalculateTau, CalculateFromSet


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
  real(RLEN)  :: sK4K3
  real(RLEN)  :: sK3G4
  real(RLEN)  :: diff
  real(RLEN)  :: gamma
  real(RLEN)  :: labda
  real(RLEN)  :: a11
  real(RLEN)  :: a12
  real(RLEN)  :: a15
  real(RLEN)  :: n12
  real(RLEN)  :: cK3n
  real(RLEN)  :: Tau
  real(RLEN)  :: zATo
  real(RLEN)  :: alpha
  real(RLEN)  :: jK3G4n
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
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M3n(BoxNumberXY) = K3n(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_p+ &
        1.0D+00)/( D2m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D6.m is the average penetration depth for C-detritus
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      alpha  =   1.0D+00/ D6m(BoxNumberXY)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate Mineralization m2 --> m3 porewater
      ! Calculate the total anoxic mineralization in mmol O/m3/d
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      zATo = rrATo(BoxNumberXY)/ p_poro(BoxNumberXY)/ IntegralExp( - &
        alpha, p_d_tot- D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      ! denitrification: temperature and coupling with anoxic mineralization
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY)* p_poro(BoxNumberXY)* &
        eTq( ETW_Ben(BoxNumberXY), p_q10diff)
      sK3G4  =   p_sK3G4* eTq(  ETW_Ben(BoxNumberXY),  p_q10)* zATo/ p_zATo

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Technical correction when K6.r is nearly zero. This denitrification
      ! is limited because not enough red. equiv. are present to oxidize &
      ! material
      ! Calculate net consumption of reduction equivalents and limit &
      ! denitrifaction rate
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      r = max( 0.0D+00, sK3G4* p_qro* p_qon_dentri* &
        K3n(BoxNumberXY)- rrATo(BoxNumberXY)* p_qro)
      sK3G4  =   sK3G4* K6r(BoxNumberXY)/( r+ K6r(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the exponential terms of the solution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      gamma  =   sqrt(  sK3G4/ diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Get coefficients describing ammonium in the oxic layer :
      ! 1. labda of the exponential curve
      ! 2. parameter of the nitrification term
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      labda = GetInfoFromSet( KNH4(BoxNumberXY), GET, LABDA_1, 11, dummy, &
        dummy)
      sK4K3 = GetInfoFromSet( KNH4(BoxNumberXY), GET, LABDA_2, 12, dummy, &
        dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Get coefficients of all terms of equation valid for the
      ! first layer of ammonium (integration constants)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      a11 = GetInfoFromSet( KNH4(BoxNumberXY), GET, COEFFICIENT, 11, &
        dummy, dummy)
      a12 = GetInfoFromSet( KNH4(BoxNumberXY), GET, COEFFICIENT, 12, &
        dummy, dummy)
      a15 = GetInfoFromSet( KNH4(BoxNumberXY), GET, COEFFICIENT, 15, &
        dummy, dummy)
      n12  =   sK4K3* a12

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layer depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      KNO3(BoxNumberXY) = InitializeSet( KNO3(BoxNumberXY), N_layers, &
        N_coeff)
      call DefineSet( KNO3(BoxNumberXY), LAYERS, LAYER1, 0, &
        D1m(BoxNumberXY), dummy)

      call DefineSet( KNO3(BoxNumberXY), DIFFUSION, FOR_ALL_LAYERS, 0, diff, &
        dummy)

      call DefineSet( KNO3(BoxNumberXY), POROSITY, FOR_ALL_LAYERS, 0, &
        p_poro(BoxNumberXY), dummy)

      call DefineSet( KNO3(BoxNumberXY), ADSORPTION, FOR_ALL_LAYERS, 0, p_p, &
        dummy)



      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! N(z) = n11*exp(labda*z) + n12*exp(-labda*z) + n13*z^2 + n14*z + n15
      ! 2nd layer:
      ! N(z) = n21*exp(gamma*z) + n22*exp(-gamma*z)
      !    n22 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call DefineSet( KNO3(BoxNumberXY), DOUBLE_DEFINE, 11, &
        ZERO_EXPONENTIAL_TERM, labda, sK4K3)

      call DefineSet( KNO3(BoxNumberXY), DOUBLE_DEFINE, 12, &
        ZERO_EXPONENTIAL_TERM, - labda, sK4K3)

      call DefineSet( KNO3(BoxNumberXY), DEFINE, 13, QUADRATIC_TERM, dummy, &
        dummy)

      call DefineSet( KNO3(BoxNumberXY), DEFINE, 14, LINEAR_TERM, dummy, &
        dummy)

      call DefineSet( KNO3(BoxNumberXY), DEFINE, 15, CONSTANT_TERM, dummy, &
        dummy)


      call DefineSet( KNO3(BoxNumberXY), DOUBLE_DEFINE, 21, EXPONENTIAL_TERM, - &
        gamma, sK3G4)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call CompleteSet( KNO3(BoxNumberXY), SET_CONTINUITY, FLAG, MASS, &
        dummy, dummy)

      call CompleteSet( KNO3(BoxNumberXY), SET_BOUNDARY, LAYER1, &
        EQUATION, 0.0D+00, N3n_Ben(BoxNumberXY))


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! a11 / (labda * labda * diff) = a12 / (labda * labda * diff)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      select case ( abs(a11)< 1.0D-30)

        case( .TRUE. )
          call CompleteSet( KNO3(BoxNumberXY), INPUT_TERM, 11, PARAMETER, &
            dummy, 0.0D+00)




        case( .FALSE. )
          call CompleteSet( KNO3(BoxNumberXY), START_ADD_TERM, 11, PARAMETER, &
            dummy, a11)

          call CompleteSet( KNO3(BoxNumberXY), INPUT_SUBTRACT_TERM, 12, &
            PARAMETER, dummy, a12)




      end select


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! a11 / (labda * labda * diff) = a15 / (2 * diff)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      select case ( abs(a15)< 1.0D-30)

        case( .TRUE. )
          call CompleteSet( KNO3(BoxNumberXY), INPUT_TERM, 13, PARAMETER, &
            dummy, 0.0D+00)




        case( .FALSE. )
          call CompleteSet( KNO3(BoxNumberXY), START_ADD_TERM, 12, PARAMETER, &
            dummy, a12)

          call CompleteSet( KNO3(BoxNumberXY), INPUT_SUBTRACT_TERM, 13, &
            PARAMETER, dummy, a15)



      end select


      call CompleteSet( KNO3(BoxNumberXY), INPUT_TERM, 12, PARAMETER, dummy, &
        n12)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate for the above defined set of boundary conditions
      ! the steady-state profiles and return the vertically integrated
      ! concentration
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cK3n = CalculateSet( KNO3(BoxNumberXY), SET_LAYER_INTEGRAL_UNTIL, &
        LAYER1, LAYER2, D2m(BoxNumberXY), 0.0D+00)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the adaptation time to the steady-state profile
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Tau  =   CalculateTau(  0.0D+00,  diff,  p_p,  D2m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate the average value of K3n over the actual time step
      ! (transient value).
      ! This value depends on the adaptation time, the actual time step,
      ! the ''old'' value and the ''equilibrium value''
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cK3n = cK3n+( K3n(BoxNumberXY)- cK3n)* IntegralExp( - LocalDelta/ &
        Tau, 1.0D+00)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Derive the equations for the transient profiles, assuming the same
      ! solution as for the steady-state case and using cK3n as new constraint.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      dummy  =   CalculateSet(  KNO3(BoxNumberXY),  ADD,  0,  0,  dummy,  cK3n)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Start calculation of fluxes:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Flux at the water/sediment interface
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK3N3n(BoxNumberXY) = CalculateFromSet( KNO3(BoxNumberXY), DERIVATIVE, &
        RFLUX, 0.0D+00, dummy)
      call flux(BoxNumberXY, iiBen, ppK3n, ppK3n, -( jK3N3n(BoxNumberXY)) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Flux at the lower boundary
      ! All transport between layers are done in BenNitrogenShifting:
      !
      ! jK13K3n = CalculateFromSet(KNO3, DERIVATIVE, RFLUX, D2.m, dummy);
      !
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Nitrification is already calculated by BenAmmonium:
      ! K4.n -> K3.n= CalculateFromSet(KM4n, INTEGRAL, RFLUX, 0.0, D1.m)*sK4K3;
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Denitrification:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK3G4n = CalculateFromSet( KNO3(BoxNumberXY), INTEGRAL, &
        RFLUX, D1m(BoxNumberXY), D2m(BoxNumberXY))* sK3G4
      call flux(BoxNumberXY, iiBen, ppK3n, ppG4n, jK3G4n )


    end DO


  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
