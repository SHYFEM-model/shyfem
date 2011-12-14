#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenQ1Transport
!
! DESCRIPTION
!   Description of the DOM dynamics in the sediments
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenQ1TransportDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: Q11c, Q1c, Q11n, Q1n, &
  ! Q11p, Q1p
  ! The following Benthic-states are used (NOT in fluxes): D1m, D6m, D2m
  ! The following global vars are modified: dummy
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY, LocalDelta
  ! The following Benthic 1-d global boxvars are modified : KQ1
  ! The following Benthic 1-d global boxvars are used: irrenh, ETW_Ben, &
  ! shiftD1m
  ! The following Benthic 2-d global boxvars  are used: ruHI, reHI
  ! The following groupmember vars  are used: iiH1, iiH2
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global box parametes are used: p_q10diff, p_clDxm, p_d_tot
  ! The following global constants are used: RLEN
  ! The following constants are used: &
  ! LAYERS, LAYER1, LAYER2, DIFFUSION, FOR_ALL_LAYERS, &
  ! POROSITY, ADSORPTION, DEFINE, EXPONENTIAL_TERM, CONSTANT_TERM, &
  ! ZERO_EXPONENTIAL_TERM, LINEAR_TERM, SET_CONTINUITY, FLAG, MASS, &
  ! SET_BOUNDARY, DERIVATIVE, SET_LAYER_INTEGRAL, LAYER3, INPUT_TERM, &
  ! STANDARD, START_ADD_TERM, INPUT_SUBTRACT_TERM, SET_LAYER_INTEGRAL_UNTIL, &
  ! ADD, SHIFT, RFLUX, ONE_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: Q11c, Q1c, Q11n, Q1n, Q11p, Q1p, D1m, D6m, D2m, D2STATE
  use mem, ONLY: ppQ11c, ppQ1c, ppQ11n, ppQ1n, ppQ11p, ppQ1p, &
    ppD1m, ppD6m, ppD2m, dummy, BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, &
    BoxNumberY, NO_BOXES_Y, BoxNumber, BoxNumberXY, LocalDelta, KQ1, irrenh, &
    ETW_Ben, shiftD1m, ruHI, reHI, iiH1, iiH2, iiBen, iiPel, flux
  use constants, ONLY: LAYERS, LAYER1, LAYER2, &
    DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, DEFINE, EXPONENTIAL_TERM, &
    CONSTANT_TERM, ZERO_EXPONENTIAL_TERM, LINEAR_TERM, SET_CONTINUITY, FLAG, &
    MASS, SET_BOUNDARY, DERIVATIVE, SET_LAYER_INTEGRAL, LAYER3, INPUT_TERM, &
    STANDARD, START_ADD_TERM, INPUT_SUBTRACT_TERM, SET_LAYER_INTEGRAL_UNTIL, ADD, &
    SHIFT, RFLUX, ONE_PER_DAY
  use mem_Param,  ONLY: p_poro, p_q10diff, p_clDxm, p_d_tot
  use mem_BenQ1Transport


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
  ! The following sesame functions are used:IntegralExp, insw
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: IntegralExp, insw



!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
!
! !REVISION_HISTORY
!   by Piet Ruardij  *:0 at Sun Dec 04 23:09:55 CET 2005
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
  real(RLEN)  :: M
  real(RLEN)  :: a15
  real(RLEN)  :: r
  real(RLEN)  :: alpha
  real(RLEN)  :: gamma
  real(RLEN)  :: diff
  real(RLEN)  :: Tau
  real(RLEN)  :: cQ1c
  real(RLEN)  :: zu
  real(RLEN)  :: sQ1
  real(RLEN)  :: zuD1
  real(RLEN)  :: rQ11
  real(RLEN)  :: jQ11Q1c
  real(RLEN)  :: jQ1Q11c
  real(RLEN)  :: flow
  real(RLEN)  :: Dnew

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
      ! Get total Net Benthic DOC (Q1.c)
      ! production/consumption in the oxic layer (m2 --> m3 porewater)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      sQ1 = max( 0.001D+00, ruHI(iiH1, BoxNumberXY)/ &
        D1m(BoxNumberXY)/ p_poro(BoxNumberXY))/( 1.0D-80+ Q1c(BoxNumberXY))
      M  =   reHI(iiH1,BoxNumberXY)

      a15  =   M/ sQ1

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY)* p_poro(BoxNumberXY)* &
        eTq( ETW_Ben(BoxNumberXY), p_q10diff)
      gamma  =   sqrt(  sQ1/ diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate specific bacterial consumption rate in anoxic layers (limited)
      ! Calculate coefficient for the exponential terms of the solution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      alpha  =   1.0D+00/ max(  p_clDxm,  D6m(BoxNumberXY))
      rQ11 = ( reHI(iiH2,BoxNumberXY)- ruHI(iiH2,BoxNumberXY))/ &
        p_poro(BoxNumberXY)/( p_d_tot- D1m(BoxNumberXY))
      zuD1 = max( 1.D-20, rQ11)/ IntegralExp( - alpha, p_d_tot- &
        D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layers depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      KQ1(BoxNumberXY) = InitializeSet( KQ1(BoxNumberXY), N_layers, &
        N_coeff)

      call DefineSet( KQ1(BoxNumberXY), LAYERS, LAYER1, &
        LAYER2, D1m(BoxNumberXY), D2m(BoxNumberXY))

      call DefineSet( KQ1(BoxNumberXY), DIFFUSION, FOR_ALL_LAYERS, 0, diff, &
        dummy)

      call DefineSet( KQ1(BoxNumberXY), POROSITY, FOR_ALL_LAYERS, 0, &
        p_poro(BoxNumberXY), dummy)

      call DefineSet( KQ1(BoxNumberXY), ADSORPTION, FOR_ALL_LAYERS, 0, p_p, &
        dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! Q(z) = q13*z^2 + q14*z + q15
      ! 2nd layer:
      ! Q(z) = q21*exp(gamma*z) + q22*exp(-gamma*z)
      ! 3rd layer:
      ! Q(z) = q31*exp(gamma*z) + q32*exp(-gamma*z) `
      !    q32 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call DefineSet( KQ1(BoxNumberXY), DEFINE, 11, EXPONENTIAL_TERM, gamma, &
        dummy)

      call DefineSet( KQ1(BoxNumberXY), DEFINE, 12, EXPONENTIAL_TERM, - &
        gamma, dummy)

      call DefineSet( KQ1(BoxNumberXY), DEFINE, 15, CONSTANT_TERM, dummy, &
        dummy)


      call DefineSet( KQ1(BoxNumberXY), DEFINE, 21, ZERO_EXPONENTIAL_TERM, - &
        alpha, dummy)

      call DefineSet( KQ1(BoxNumberXY), DEFINE, 24, LINEAR_TERM, dummy, &
        dummy)

      call DefineSet( KQ1(BoxNumberXY), DEFINE, 25, CONSTANT_TERM, dummy, &
        dummy)


      call DefineSet( KQ1(BoxNumberXY), DEFINE, 31, ZERO_EXPONENTIAL_TERM, - &
        alpha, dummy)

      call DefineSet( KQ1(BoxNumberXY), DEFINE, 35, CONSTANT_TERM, dummy, &
        dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !4
      call CompleteSet( KQ1(BoxNumberXY), SET_CONTINUITY, FLAG, MASS, dummy, &
        dummy)

      !5
      call CompleteSet( KQ1(BoxNumberXY), SET_BOUNDARY, LAYER1, DERIVATIVE, &
        0.0D+00, 0.0D+00)

      !6:
      call CompleteSet( KQ1(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER2, &
        LAYER3, dummy, Q11c(BoxNumberXY))

      !7:
      r  =   exp( - alpha*( D1m(BoxNumberXY)- D2m(BoxNumberXY)))

      select case ( r> 1.0D-20)

        case( .FALSE. )
          call CompleteSet( KQ1(BoxNumberXY), INPUT_TERM, 31, STANDARD, &
            dummy, 0.0D+00)




        case( .TRUE. )
          call CompleteSet( KQ1(BoxNumberXY), START_ADD_TERM, 31, STANDARD, &
            dummy, r)

          call CompleteSet( KQ1(BoxNumberXY), INPUT_SUBTRACT_TERM, 21, STANDARD, &
            dummy, 1.0D+00)



      end select

      !8:
      call CompleteSet( KQ1(BoxNumberXY), INPUT_TERM, 15, STANDARD, dummy, &
        a15)



      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate for the above defined set of boundary conditions
      ! the steady-state profiles and return the vertically integrated
      ! concentration
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cQ1c = CalculateSet( KQ1(BoxNumberXY), SET_LAYER_INTEGRAL_UNTIL, &
        LAYER1, LAYER1, D1m(BoxNumberXY), 0.0D+00)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the adaptation time to the steady-state profile
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Tau  =   CalculateTau(  sQ1,  diff,  p_p,  D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate the average value of Q11 over the actual time step
      ! (transient value).
      ! This value depends on the adaptation time, the actual time step,
      ! the ''old'' value and the ''equilibrium value''
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cQ1c = cQ1c+( Q1c(BoxNumberXY)- cQ1c)* IntegralExp( - LocalDelta/ &
        Tau, 1.0D+00)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Derive the equations for the transient profiles, assuming the same
      ! solution as for the steady-state case and using cQ11c as new constraint.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      dummy  =   CalculateSet(  KQ1(BoxNumberXY),  ADD,  0,  0,  dummy,  cQ1c)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Start calculation of fluxes:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate new depth of the oxygen horizon and the flux related
      ! to the shifting
      ! Add the flux at D1.m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Dnew  =   D1m(BoxNumberXY)+ LocalDelta* shiftD1m(BoxNumberXY)
      flow = CalculateFromSet( KQ1(BoxNumberXY), SHIFT, LAYER1, &
        D1m(BoxNumberXY), Dnew)/ LocalDelta

      flow = flow+ CalculateFromSet( KQ1(BoxNumberXY), DERIVATIVE, &
        RFLUX, D1m(BoxNumberXY), dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Limit for too large fluxes
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      select case ( flow< 0.0D+00)

        case( .TRUE. )
          jQ1Q11c = (- flow* Q1c(BoxNumberXY)/ ONE_PER_DAY/( &
            Q1c(BoxNumberXY)/ ONE_PER_DAY- flow))* insw( - flow)
          jQ11Q1c  =   0.0D+00



        case( .FALSE. )
          jQ11Q1c = ( flow* Q11c(BoxNumberXY)/ ONE_PER_DAY/( Q11c(BoxNumberXY)/ &
            ONE_PER_DAY+ flow))* insw( flow)
          jQ1Q11c  =   0.0D+00



      end select

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! One of the 2 fluxes between the some constituents is 0!
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call flux(BoxNumberXY, iiBen, ppQ1c, ppQ11c, jQ1Q11c )
      call flux(BoxNumberXY, iiBen, ppQ11c, ppQ1c, jQ11Q1c )

      call flux(BoxNumberXY, iiBen, ppQ1n, ppQ11n, jQ1Q11c/ Q1c(BoxNumberXY)* &
        Q1n(BoxNumberXY) )
      call flux(BoxNumberXY, iiBen, ppQ11n, ppQ1n, jQ11Q1c/ Q11c(BoxNumberXY)* &
        Q11n(BoxNumberXY) )

      call flux(BoxNumberXY, iiBen, ppQ1p, ppQ11p, jQ1Q11c/ Q1c(BoxNumberXY)* &
        Q1p(BoxNumberXY) )
      call flux(BoxNumberXY, iiBen, ppQ11p, ppQ1p, jQ11Q1c/ Q11c(BoxNumberXY)* &
        Q11p(BoxNumberXY) )




    end DO


  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
