#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenAnoxic
!
! DESCRIPTION
!   Description of the anoxic diagenetic processes in the sediment
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
  subroutine BenAnoxicDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: K6r, G2o
  ! The following Benthic-states are used (NOT in fluxes): D6m, D1m
  ! The following global vars are modified: dummy
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY, LocalDelta
  ! The following Benthic 1-d global boxvars are modified : M6r, KRED, jK6N6r, &
  ! jG2K7o
  ! The following Benthic 1-d global boxvars are used: rrATo, irrenh, ETW_Ben, &
  ! KNO3, N6r_Ben
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global box parametes are used: p_d_tot, p_clDxm, &
  ! p_qro, p_q10diff, p_qon_dentri
  ! The following global constants are used: RLEN
  ! The following constants are used: GET, &
  ! LABDA_1, LABDA_2, COEFFICIENT, LAYERS, LAYER1, &
  ! DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, DEFINE, &
  ! EXPONENTIAL_TERM, ZERO_EXPONENTIAL_TERM, DOUBLE_DEFINE, CONSTANT_TERM, &
  ! SET_CONTINUITY, FLAG, MASS, SET_BOUNDARY, EQUATION, INPUT_TERM, &
  ! PARAMETER, SET_LAYER_INTEGRAL_UNTIL, LAYER2, ADD, DERIVATIVE, RFLUX, &
  ! INTEGRAL, ONE_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: K6r, G2o, D6m, D1m, D2STATE
  use mem, ONLY: ppK6r, ppG2o, ppD6m, ppD1m, dummy, BoxNumberZ, &
    NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, BoxNumber, &
    BoxNumberXY, LocalDelta, M6r, KRED, jK6N6r, jG2K7o, rrATo, irrenh, ETW_Ben, &
    KNO3, N6r_Ben, iiBen, iiPel, flux
  use constants, ONLY: GET, LABDA_1, LABDA_2, COEFFICIENT, &
    LAYERS, LAYER1, DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, &
    DEFINE, EXPONENTIAL_TERM, ZERO_EXPONENTIAL_TERM, DOUBLE_DEFINE, &
    CONSTANT_TERM, SET_CONTINUITY, FLAG, MASS, SET_BOUNDARY, EQUATION, &
    INPUT_TERM, PARAMETER, SET_LAYER_INTEGRAL_UNTIL, LAYER2, ADD, DERIVATIVE, &
    RFLUX, INTEGRAL, ONE_PER_DAY
  use mem_Param,  ONLY: p_poro, p_d_tot, p_clDxm, p_qro, p_q10diff, p_qon_dentri
  use mem_BenAnoxic


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
  real(RLEN)  :: gamma
  real(RLEN)  :: alpha
  real(RLEN)  :: zuD1
  real(RLEN)  :: diff
  real(RLEN)  :: labda
  real(RLEN)  :: sK3G4
  real(RLEN)  :: n21
  real(RLEN)  :: Tau
  real(RLEN)  :: cK6r
  real(RLEN)  :: jATK6r
  real(RLEN)  :: jK6BTr
  real(RLEN)  :: jK6G4r
  real(RLEN)  :: jK16K6r

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

      M6r(BoxNumberXY) = K6r(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_p+ 1.0D+00)/( &
        p_d_tot)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D6.m is the average penetration depth for C-detritus
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      alpha  =   1.0D+00/ max(  p_clDxm,  D6m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Convert anoxic mineralization (mmol S/m2/d)
      ! This rate is already assigned to the dynamical equation for K6.r
      ! in BenBacDyanmics for H2:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jATK6r  =   p_qro* rrATo(BoxNumberXY)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate Mineralization m2 --> m3 porewater
      ! Anoxic mineralization at D1.m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      zuD1 = jATK6r/ p_poro(BoxNumberXY)/ IntegralExp( - alpha, &
        p_d_tot- D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY)* p_poro(BoxNumberXY)* &
        eTq( ETW_Ben(BoxNumberXY), p_q10diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the exponential terms of the solution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      gamma  =   sqrt(  p_rOS/ diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Get coefficients describing Nitrate in anoxic layer :
      ! 1. labda of the exponential curve, and the denitrification rate
      ! 2. parameter of the denitrification term (integration constant)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      labda = GetInfoFromSet( KNO3(BoxNumberXY), GET, LABDA_1, 21, dummy, &
        dummy)
      sK3G4 = GetInfoFromSet( KNO3(BoxNumberXY), GET, LABDA_2, 21, &
        dummy, dummy)* p_qro* p_qon_dentri
      n21 = - GetInfoFromSet( KNO3(BoxNumberXY), GET, COEFFICIENT, 21, &
        dummy, dummy)* sK3G4

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layers depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      KRED(BoxNumberXY) = InitializeSet( KRED(BoxNumberXY), N_layers, &
        N_coeff)
      call DefineSet( KRED(BoxNumberXY), LAYERS, LAYER1, 0, &
        D1m(BoxNumberXY), dummy)

      call DefineSet( KRED(BoxNumberXY), DIFFUSION, FOR_ALL_LAYERS, 0, diff, &
        dummy)

      call DefineSet( KRED(BoxNumberXY), POROSITY, FOR_ALL_LAYERS, 0, &
        p_poro(BoxNumberXY), dummy)

      call DefineSet( KRED(BoxNumberXY), ADSORPTION, FOR_ALL_LAYERS, 0, p_p, &
        dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! R(z) = r11*exp(gamma*z) + r12*exp(-gamma*z)
      ! 2nd layer:
      ! R(z) = r21*exp[-alpha*(z-D1.m)] + r22*exp(-labda*z) + r23*z^2 + r24*z + &
      ! r25
      !    r23, r24 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call DefineSet( KRED(BoxNumberXY), DEFINE, 11, EXPONENTIAL_TERM, - &
        gamma, dummy)

      call DefineSet( KRED(BoxNumberXY), DEFINE, 12, EXPONENTIAL_TERM, &
        gamma, dummy)

      call DefineSet( KRED(BoxNumberXY), DEFINE, 21, ZERO_EXPONENTIAL_TERM, - &
        alpha, dummy)

      call DefineSet( KRED(BoxNumberXY), DOUBLE_DEFINE, 22, &
        ZERO_EXPONENTIAL_TERM, labda, sK3G4)

      call DefineSet( KRED(BoxNumberXY), DEFINE, 25, CONSTANT_TERM, dummy, &
        dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call CompleteSet( KRED(BoxNumberXY), SET_CONTINUITY, FLAG, MASS, &
        dummy, dummy)

      call CompleteSet( KRED(BoxNumberXY), SET_BOUNDARY, LAYER1, &
        EQUATION, 0.0D+00, N6r_Ben(BoxNumberXY))


      call CompleteSet( KRED(BoxNumberXY), INPUT_TERM, 22, PARAMETER, dummy, &
        n21)

      call CompleteSet( KRED(BoxNumberXY), INPUT_TERM, 21, PARAMETER, dummy, &
        zuD1)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate for the above defined set of boundary conditions
      ! the steady-state profiles and return the vertically integrated
      ! concentration.
      !
      ! Technical improvements: in case of utlimate low mineralization rates and
      ! a nearly empty reduction equivalent pool. There is a chance the
      ! estimated equilibrium value is negative. Therfore cK6r is limited to &
      ! values >=0
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cK6r = max( 0.0D+00, CalculateSet( KRED(BoxNumberXY), &
        SET_LAYER_INTEGRAL_UNTIL, LAYER1, LAYER2, p_d_tot, 0.0D+00))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the adaptation time to the steady-state profile
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Tau  =   CalculateTau(  sK3G4,  diff,  p_p,  p_d_tot)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate the average value of K6r over the actual time step
      ! (transient value).
      ! This value depends on the adaptation time, the actual time step,
      ! the ''old'' value and the ''equilibrium value''
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cK6r = cK6r+( K6r(BoxNumberXY)- cK6r)* IntegralExp( - LocalDelta/ &
        Tau, 1.0D+00)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Derive the equations for the transient profiles, assuming the same
      ! solution as for the steady-state case and using cK6r as new constraint.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      dummy  =   CalculateSet(  KRED(BoxNumberXY),  ADD,  0,  0,  dummy,  cK6r)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimation of the Vertical flux at surface
      ! from the set of transient solutions:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK6N6r(BoxNumberXY) = CalculateFromSet( KRED(BoxNumberXY), DERIVATIVE, &
        RFLUX, 0.0D+00, dummy)
      jK6N6r(BoxNumberXY) = jK6N6r(BoxNumberXY)* insw( ( &
        M6r(BoxNumberXY)- N6r_Ben(BoxNumberXY))* jK6N6r(BoxNumberXY))


      call flux(BoxNumberXY, iiBen, ppK6r, ppK6r, -( jK6N6r(BoxNumberXY)) )


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Reoxidation Flux from (S2-, Fe3+, Mg3+ to SO4-, FE2+, Mg2+
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK6BTr = p_rOS* CalculateFromSet( KRED(BoxNumberXY), INTEGRAL, &
        RFLUX, 0.0D+00, D1m(BoxNumberXY))

      jG2K7o(BoxNumberXY)  =   jK6BTr/ p_qro

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! loss flux due to denitrification
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK6G4r = - GetInfoFromSet( KRED(BoxNumberXY), INTEGRAL, PARAMETER, &
        22, D1m(BoxNumberXY), p_d_tot)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! flux at the lower boundary
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK16K6r = max( CalculateFromSet( KRED(BoxNumberXY), DERIVATIVE, &
        RFLUX, p_d_tot, dummy), - 0.1D+00* K6r(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! At low value of K6.r there is possibility that K6.r becomes negative
      ! The surplus of loss due to denitrification is now directly
      ! subtracted from the oxygen consumption.
      ! Make denitrification equal to anoxic mineralization
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if ( jK6N6r(BoxNumberXY)+ jK6G4r+ jK6BTr> K6r(BoxNumberXY)/ ONE_PER_DAY+ &
        jATK6r+ jK16K6r) then

        jK6BTr  =   0.0D+00
        jG2K7o(BoxNumberXY)  =  -( jK6G4r+ jK6N6r(BoxNumberXY)- jATK6r)/ p_qro

        jK16K6r  =   0.0D+00
        jK6G4r  =   jATK6r

      end if


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! execute flux statements:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call flux(BoxNumberXY, iiBen, ppG2o, ppG2o, -( jG2K7o(BoxNumberXY)) )
      call flux(BoxNumberXY, iiBen, ppK6r, ppK6r, -( jK6BTr) )
      call flux(BoxNumberXY, iiBen, ppK6r, ppK6r, -( jK6G4r) )
      call flux(BoxNumberXY, iiBen, ppK6r, ppK6r, -(- jK16K6r) )

    end DO


  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
