#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenPhosphate
!
! DESCRIPTION
!   Description of the phosphate diagenitic processes in the sediment
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
  subroutine BenPhosphateDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: K1p, K11p, K21p
  ! The following Benthic-states are used (NOT in fluxes): D1m, D2m, D8m
  ! The following global vars are modified: dummy
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY, LocalDelta
  ! The following Benthic 1-d global boxvars are modified : M1p, M11p, M21p, &
  ! KPO4, jK1N1p
  ! The following Benthic 1-d global boxvars are used: reBTp, &
  ! reATp, irrenh, ETW_Ben, N1p_Ben, Depth_Ben, shiftD1m, shiftD2m
  ! The following Benthic 1-d global boxpars  are used: p_poro, p_p_ae
  ! The following 0-d global box parametes are used: p_d_tot, p_clDxm, p_q10diff
  ! The following global constants are used: RLEN
  ! The following constants are used: QUADRATIC_TERM, &
  ! ZERO_EXPONENTIAL_TERM, LAYERS, LAYER1, LAYER2, LAYER3, &
  ! DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, LAYER4, &
  ! DEFINE, LINEAR_TERM, CONSTANT_TERM, SET_CONTINUITY, FLAG, &
  ! MASS, SET_BOUNDARY, EQUATION, SET_LAYER_INTEGRAL, SET_LAYER_INTEGRAL_UNTIL, &
  ! INPUT_TERM, PARAMETER, START_ADD_TERM, STANDARD, INPUT_SUBTRACT_TERM, ADD, &
  ! DERIVATIVE, RFLUX, SHIFT, ONE_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: K1p, K11p, K21p, D1m, D2m, D8m, D2STATE
  use mem, ONLY: ppK1p, ppK11p, ppK21p, ppD1m, ppD2m, ppD8m, &
    dummy, BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
    BoxNumber, BoxNumberXY, LocalDelta, M1p, M11p, M21p, KPO4, jK1N1p, reBTp, &
    reATp, irrenh, ETW_Ben, N1p_Ben, Depth_Ben, shiftD1m, shiftD2m, iiBen, iiPel, &
    flux
  use constants, ONLY: QUADRATIC_TERM, ZERO_EXPONENTIAL_TERM, LAYERS, &
    LAYER1, LAYER2, LAYER3, DIFFUSION, FOR_ALL_LAYERS, POROSITY, &
    ADSORPTION, LAYER4, DEFINE, LINEAR_TERM, CONSTANT_TERM, SET_CONTINUITY, &
    FLAG, MASS, SET_BOUNDARY, EQUATION, SET_LAYER_INTEGRAL, &
    SET_LAYER_INTEGRAL_UNTIL, INPUT_TERM, PARAMETER, START_ADD_TERM, &
    STANDARD, INPUT_SUBTRACT_TERM, ADD, DERIVATIVE, RFLUX, SHIFT, ONE_PER_DAY
  use mem_Param,  ONLY: p_poro, p_p_ae, p_d_tot, p_clDxm, p_q10diff
  use mem_BenPhosphate


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
!   April 15, 1994 by EGM Embsen and P Ruardij:
!               Created a new version of the this process
!               so that it can be used with OpenSESAME.
!       September 1999 by M. Vichi
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
  real(RLEN)  :: diff
  real(RLEN)  :: dx
  real(RLEN)  :: jK11K1p
  real(RLEN)  :: jK21K11p
  real(RLEN)  :: jK31K21p
  real(RLEN)  :: zuBT
  real(RLEN)  :: zuD1
  real(RLEN)  :: zuD2
  real(RLEN)  :: alpha
  real(RLEN)  :: cK1p
  real(RLEN)  :: Tau
  real(RLEN)  :: clM1p
  real(RLEN)  :: r
  real(RLEN)  :: Dnew
  integer  :: term

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

      M1p(BoxNumberXY) = K1p(BoxNumberXY)/ p_poro(BoxNumberXY)/( &
        p_p_ae(BoxNumberXY)+ 1.0D+00)/( D1m(BoxNumberXY))
      M11p(BoxNumberXY) = K11p(BoxNumberXY)/ p_poro(BoxNumberXY)/( &
        p_p_ae(BoxNumberXY)+ 1.0D+00)/( D2m(BoxNumberXY)- D1m(BoxNumberXY))
      M21p(BoxNumberXY) = K21p(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_p_an+ &
        1.0D+00)/( p_d_tot- D2m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D8.m is the average penetration depth for P-detritus
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      alpha  =   1.0D+00/ max(  p_clDxm,  D8m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate Mineralization m2 --> m3 porewater
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Average in the oxic layer:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      select case ( D1m(BoxNumberXY)< p_chD1m)

        case( .TRUE. )
          zuBT = max( 1.079D-6, reBTp(BoxNumberXY))/ &
            p_poro(BoxNumberXY)/ D1m(BoxNumberXY)
          term  =   QUADRATIC_TERM




        case( .FALSE. )
          zuBT = max( 1.D-20, reBTp(BoxNumberXY))/ p_poro(BoxNumberXY)/ &
            IntegralExp( - alpha, D1m(BoxNumberXY))
          term  =   ZERO_EXPONENTIAL_TERM


      end select


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Anoxic Mineralization at D1.m, using the exponential distribution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      zuD1 = max( 1.D-20, reATp(BoxNumberXY))/ p_poro(BoxNumberXY)/ IntegralExp( &
        - alpha, p_d_tot- D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Anoxic Mineralization at D2.m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      zuD2  =   zuD1* exp( - alpha*( D2m(BoxNumberXY)- D1m(BoxNumberXY)))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY)* p_poro(BoxNumberXY)* &
        eTq( ETW_Ben(BoxNumberXY), p_q10diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Subdivide anoxic layer in two sublayers only for the calculation.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      dx  =  ( p_d_tot+ D2m(BoxNumberXY))* 0.5D+00

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layers depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      KPO4(BoxNumberXY) = InitializeSet( KPO4(BoxNumberXY), N_layers, &
        N_coeff)

      call DefineSet( KPO4(BoxNumberXY), LAYERS, LAYER1, &
        LAYER2, D1m(BoxNumberXY), D2m(BoxNumberXY))

      call  DefineSet(  KPO4(BoxNumberXY),  LAYERS,  LAYER3,  0,  dx,  dummy)

      call DefineSet( KPO4(BoxNumberXY), DIFFUSION, FOR_ALL_LAYERS, 0, diff, &
        dummy)

      call DefineSet( KPO4(BoxNumberXY), POROSITY, FOR_ALL_LAYERS, 0, &
        p_poro(BoxNumberXY), dummy)

      call DefineSet( KPO4(BoxNumberXY), ADSORPTION, LAYER1, &
        LAYER2, p_p_ae(BoxNumberXY), p_p_ae(BoxNumberXY))

      call DefineSet( KPO4(BoxNumberXY), ADSORPTION, LAYER3, LAYER4, &
        p_p_an, p_p_an)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! P(z) = p13*z^2 + p14*z + p15
      ! 2nd layer:
      ! P(z) = p21*exp[-alpha*(z-D1.m)] + p24*z + p25
      ! 3rd layer:
      ! P(z) = p31*exp[-alpha*(z-D2.m)] + p34*z + p35
      ! 4th layer:
      ! P(z) = p41*exp[-alpha*(z-dx)] + p44*z + p45
      !    p44 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call  DefineSet(  KPO4(BoxNumberXY),  DEFINE,  13,  term, - alpha,  dummy)

      call DefineSet( KPO4(BoxNumberXY), DEFINE, 14, LINEAR_TERM, dummy, &
        dummy)

      call DefineSet( KPO4(BoxNumberXY), DEFINE, 15, CONSTANT_TERM, dummy, &
        dummy)


      call DefineSet( KPO4(BoxNumberXY), DEFINE, 21, ZERO_EXPONENTIAL_TERM, - &
        alpha, dummy)

      call DefineSet( KPO4(BoxNumberXY), DEFINE, 24, LINEAR_TERM, dummy, &
        dummy)

      call DefineSet( KPO4(BoxNumberXY), DEFINE, 25, CONSTANT_TERM, dummy, &
        dummy)


      call DefineSet( KPO4(BoxNumberXY), DEFINE, 31, ZERO_EXPONENTIAL_TERM, - &
        alpha, dummy)

      call DefineSet( KPO4(BoxNumberXY), DEFINE, 34, LINEAR_TERM, dummy, &
        dummy)

      call DefineSet( KPO4(BoxNumberXY), DEFINE, 35, CONSTANT_TERM, dummy, &
        dummy)


      call DefineSet( KPO4(BoxNumberXY), DEFINE, 41, ZERO_EXPONENTIAL_TERM, - &
        alpha, dummy)

      call DefineSet( KPO4(BoxNumberXY), DEFINE, 45, CONSTANT_TERM, dummy, &
        dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call CompleteSet( KPO4(BoxNumberXY), SET_CONTINUITY, FLAG, MASS, &
        dummy, dummy)

      call CompleteSet( KPO4(BoxNumberXY), SET_BOUNDARY, LAYER1, &
        EQUATION, 0.0D+00, N1p_Ben(BoxNumberXY))

      call CompleteSet( KPO4(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER2, &
        LAYER2, dummy, K11p(BoxNumberXY))

      call CompleteSet( KPO4(BoxNumberXY), SET_LAYER_INTEGRAL_UNTIL, LAYER3, &
        LAYER4, p_d_tot, K21p(BoxNumberXY))


      r  =   exp( - alpha*( dx- D2m(BoxNumberXY)))

      select case ( r> 1.0D-20)

        case( .FALSE. )
          call CompleteSet( KPO4(BoxNumberXY), INPUT_TERM, 41, PARAMETER, &
            dummy, 0.0D+00)




        case( .TRUE. )
          call CompleteSet( KPO4(BoxNumberXY), START_ADD_TERM, 41, STANDARD, &
            dummy, r)

          call CompleteSet( KPO4(BoxNumberXY), INPUT_SUBTRACT_TERM, 31, &
            STANDARD, dummy, 1.0D+00)



      end select


      call CompleteSet( KPO4(BoxNumberXY), INPUT_TERM, 13, PARAMETER, dummy, &
        zuBT)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate for the above defined set of boundary conditions
      ! the steady-state profiles and return the vertically integrated
      ! concentration in the first layer.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cK1p = CalculateSet( KPO4(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER1, &
        LAYER1, dummy, 0.0D+00)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the adaptation time to the steady-state profile
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Tau = CalculateTau( 0.0D+00, diff, p_p_ae(BoxNumberXY), &
        D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate the average value of K1p over the actual time step
      ! (transient value).
      ! This value depends on the adaptation time, the actual time step,
      ! the ''old'' value and the ''equilibrium value''
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cK1p = cK1p+( K1p(BoxNumberXY)- cK1p)* IntegralExp( - LocalDelta/ &
        Tau, 1.0D+00)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Derive the equations for the transient profiles, assuming the same
      ! solution as for the steady-state case and using cK1p as new constraint.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      dummy  =   CalculateSet(  KPO4(BoxNumberXY),  ADD,  0,  0,  dummy,  cK1p)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Start calculation of fluxes:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Vertical fluxes :
      ! There are 2 problems with this model in this version both connected
      ! with the shifting of the layers:
      ! 1. shifting from the oxic+denitrification layer with high
      !  adsorped fraction phosphate to the lower anoxic layer with a very
      !  low percentage of adsorped phosphate.
      ! 2. Too large changes in spring due to large change of D1.m:
      !  This lead sometimes to a calculated phosphate gradient which
      !  has at some depth negative values.
      !
      !  Solution:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate flux at the sediment/water interface:
      ! Check on; to high fluxes from pelagic and on concisteny of gradient
      ! ( only flux of M1p > N1p!)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK1N1p(BoxNumberXY) = max( CalculateFromSet( KPO4(BoxNumberXY), &
        DERIVATIVE, RFLUX, 0.0D+00, 0.0D+00), - 0.5D+00* &
        N1p_Ben(BoxNumberXY)* Depth_Ben(BoxNumberXY))
      jK1N1p(BoxNumberXY) = jK1N1p(BoxNumberXY)* insw( &
        jK1N1p(BoxNumberXY)*( M1p(BoxNumberXY)- N1p_Ben(BoxNumberXY)))

      call flux(BoxNumberXY, iiBen, ppK1p, ppK1p, -( jK1N1p(BoxNumberXY)) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate new depth of the oxygen horizon
      ! and the flux of phosphate related to this shifting
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Dnew  =   D1m(BoxNumberXY)+ shiftD1m(BoxNumberXY)* LocalDelta

      jK11K1p = CalculateFromSet( KPO4(BoxNumberXY), SHIFT, LAYER1, &
        D1m(BoxNumberXY), Dnew)/ LocalDelta

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! limit flux according to the actual phosphate content in the layer
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      select case ( jK11K1p> 0.0D+00)
        case( .TRUE. )
          jK11K1p = ( jK11K1p* K11p(BoxNumberXY)/ ONE_PER_DAY/( &
            K11p(BoxNumberXY)/ ONE_PER_DAY+ jK11K1p))


        case( .FALSE. )
          jK11K1p = -(- jK11K1p* K1p(BoxNumberXY)/ ONE_PER_DAY/( &
            K1p(BoxNumberXY)/ ONE_PER_DAY- jK11K1p))



      end select


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate diffusive flux at the oxic/denitrification interface:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      r = CalculateFromSet( KPO4(BoxNumberXY), DERIVATIVE, RFLUX, &
        D1m(BoxNumberXY), 1.0D+00)
      jK11K1p  =   jK11K1p+ r* insw( ( M11p(BoxNumberXY)- M1p(BoxNumberXY))* r)


      call flux(BoxNumberXY, iiBen, ppK11p, ppK1p, jK11K1p* insw( jK11K1p) )
      call flux(BoxNumberXY, iiBen, ppK1p, ppK11p, - jK11K1p* insw( - jK11K1p) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! All the nutrient mineralization source term in the anoxic layer
      ! has been added to K11.p in BenBacDynamics
      ! However in the model this layer is subdivided and hence a partition
      ! flux is here calculated according to the exponential distribution.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK21K11p = - zuD2* p_poro(BoxNumberXY)* IntegralExp( - alpha, &
        p_d_tot- D2m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate new depth of the sulphide horizon:
      ! and the flux of phosphate related to this shifting
      ! (this calculation involves the change of the adsorption coefficient)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Dnew  =   D2m(BoxNumberXY)+ shiftD2m(BoxNumberXY)* LocalDelta

      jK21K11p = jK21K11p+ CalculateFromSet( KPO4(BoxNumberXY), SHIFT, &
        LAYER2, D2m(BoxNumberXY), Dnew)/ LocalDelta

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate diffusive flux at the denitrification/anoxic interface:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK21K11p = jK21K11p+ CalculateFromSet( KPO4(BoxNumberXY), DERIVATIVE, &
        RFLUX, D2m(BoxNumberXY), 1.0D+00)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! limit flux according to the actual phosphate content in the layer
      ! if the flux is upwards.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      select case ( jK21K11p> 0.0D+00)

        case( .TRUE. )
          jK21K11p = ( jK21K11p* K21p(BoxNumberXY)/( K21p(BoxNumberXY)+ &
            jK21K11p))


        case( .FALSE. )
          jK21K11p = max( - M21p(BoxNumberXY)* &
            abs(shiftD2m(BoxNumberXY))* p_poro(BoxNumberXY)*( 1.0D+00+ p_p_an), &
            jK21K11p)



      end select

      call flux(BoxNumberXY, iiBen, ppK21p, ppK11p, jK21K11p* insw( jK21K11p) )
      call flux(BoxNumberXY, iiBen, ppK11p, ppK21p, - jK21K11p* insw( - &
        jK21K11p) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate flux at the lower boundary
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK31K21p = CalculateFromSet( KPO4(BoxNumberXY), DERIVATIVE, RFLUX, &
        p_d_tot, 0.0D+00)
      if ( jK31K21p< 0.0D+00) then
        jK31K21p = ( jK31K21p* K21p(BoxNumberXY)/( K21p(BoxNumberXY)- &
          jK31K21p))
      end if


      call flux(BoxNumberXY, iiBen, ppK21p, ppK21p, -(- jK31K21p) )

    end DO


  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
