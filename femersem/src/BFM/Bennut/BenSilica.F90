#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenSilica
!
! DESCRIPTION
!   Description of the diagenitic processes in the sediment
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
  subroutine BenSilicaDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: K5s, Q6s, D9m
  ! The following Benthic-states are used (NOT in fluxes): D1m, D2m
  ! The following global vars are modified: dummy
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY, idummy, LocalDelta
  ! The following Benthic 1-d global boxvars are modified : M5s, KSIO3, KSIO3E, &
  ! jK5N5s
  ! The following Benthic 1-d global boxvars are used: irrenh, ETW_Ben, &
  ! N5s_Ben, shiftD2m
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global box parametes are used: p_clD1D2m, p_q10diff, &
  ! p_clDxm, p_d_tot
  ! The following global constants are used: RLEN
  ! The following constants are used: LAYERS, &
  ! LAYER1, DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, &
  ! DEFINE, QUADRATIC_TERM, LINEAR_TERM, CONSTANT_TERM, PARAMETER_DEFINE, &
  ! BESSELI_EXP_TERM, SET_CONTINUITY, STANDARD, SET_BOUNDARY, EQUATION, &
  ! INPUT_TERM, PARAMETER, SET_LAYER_INTEGRAL_UNTIL, LAYER2, ADD, INTEGRAL, &
  ! DERIVATIVE, RFLUX, MASS, EXPONENTIAL_INTEGRAL

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: K5s, Q6s, D9m, D1m, D2m, D2STATE
  use mem, ONLY: ppK5s, ppQ6s, ppD9m, ppD1m, ppD2m, dummy, &
    BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
    BoxNumber, BoxNumberXY, idummy, LocalDelta, M5s, KSIO3, KSIO3E, jK5N5s, &
    irrenh, ETW_Ben, N5s_Ben, shiftD2m, iiBen, iiPel, flux
  use constants, ONLY: LAYERS, LAYER1, DIFFUSION, &
    FOR_ALL_LAYERS, POROSITY, ADSORPTION, DEFINE, QUADRATIC_TERM, LINEAR_TERM, &
    CONSTANT_TERM, PARAMETER_DEFINE, BESSELI_EXP_TERM, SET_CONTINUITY, STANDARD, &
    SET_BOUNDARY, EQUATION, INPUT_TERM, PARAMETER, SET_LAYER_INTEGRAL_UNTIL, &
    LAYER2, ADD, INTEGRAL, DERIVATIVE, RFLUX, MASS, EXPONENTIAL_INTEGRAL
  use mem_Param,  ONLY: p_poro, p_clD1D2m, p_q10diff, p_clDxm, p_d_tot
  use mem_BenSilica


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:InitializeSet, &
  ! DefineSet, CompleteSet, CalculateSet, CalculateTau, CopySet, &
  ! CalculateFromSet, GetInfoFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: InitializeSet, DefineSet, CompleteSet, &
    CalculateSet, CalculateTau, CopySet, CalculateFromSet, GetInfoFromSet


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
  integer  :: i
  real(RLEN)  :: cD1m
  real(RLEN)  :: cD2m
  real(RLEN)  :: cD2mNew
  real(RLEN)  :: cShiftD2m
  real(RLEN)  :: chM5s
  real(RLEN)  :: cM5s
  real(RLEN)  :: Tau
  real(RLEN)  :: alpha
  real(RLEN)  :: diff
  real(RLEN)  :: M5b0
  real(RLEN)  :: M5b_0_d1
  real(RLEN)  :: M5bD1
  real(RLEN)  :: zuBT
  real(RLEN)  :: suD1
  real(RLEN)  :: rmQ6s
  real(RLEN)  :: s
  real(RLEN)  :: jK15K5s
  real(RLEN)  :: jQ6K5s
  real(RLEN)  :: jQ6K15s
  real(RLEN)  :: smQ6

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
      ! Here M5s is used in the calculation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cD1m  =   min(  max(  D1m(BoxNumberXY),   p_clD1m),  p_chD2m- p_clD1D2m)
      cD2m  =   min(  max(  D2m(BoxNumberXY),   p_clD2m),  p_chD2m)
      M5s(BoxNumberXY) = K5s(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_p+ 1.0D+00)/ &
        cD2m

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      ! saturation value: temperature
      ! dissolution rate: temperature
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY)* p_poro(BoxNumberXY)* &
        eTq( ETW_Ben(BoxNumberXY), p_q10diff)
      chM5s = p_chM5s+ p_cvM5s*( eTq( ETW_Ben(BoxNumberXY), p_q10)- &
        1.0D+00)
      smQ6  =   p_smQ6* eTq(  ETW_Ben(BoxNumberXY),  p_q10diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D9.m is the average penetration depth for biogenic Si
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      alpha  =   1.0D+00/ max(  p_clDxm,  D9m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate total biogenic silica from m2 --> m3 porewater
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M5b0 = Q6s(BoxNumberXY)/ p_poro(BoxNumberXY)/ IntegralExp( - alpha, &
        p_d_tot)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Average content of Biogenic silica in the oxic layer
      ! and calculation of the zero-order dissolution term
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M5b_0_d1  =   M5b0* IntegralExp( - alpha,  cD1m)/ cD1m
      zuBT  =   smQ6* M5b_0_d1

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Biogenic silica at cD1m and calculation of the dissolution rate
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M5bD1  =   M5b0* exp( - alpha* cD1m)
      suD1  =   smQ6* M5bD1/ chM5s

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layer depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      KSIO3(BoxNumberXY) = InitializeSet( KSIO3(BoxNumberXY), N_layers, &
        N_coeff)

      call  DefineSet(  KSIO3(BoxNumberXY),  LAYERS,  LAYER1,  0,  cD1m,  dummy)

      call DefineSet( KSIO3(BoxNumberXY), DIFFUSION, FOR_ALL_LAYERS, idummy, &
        diff, dummy)

      call DefineSet( KSIO3(BoxNumberXY), POROSITY, FOR_ALL_LAYERS, &
        idummy, p_poro(BoxNumberXY), dummy)

      call DefineSet( KSIO3(BoxNumberXY), ADSORPTION, FOR_ALL_LAYERS, idummy, &
        p_p, dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! C = Ssat - S
      !
      ! 1st layer:
      ! C(z) = c13*z^2 + c14*z + c15
      ! 2nd layer:
      ! C(z) = c21*I0*exp[-alpha*(z-cD1m
      !    I0 = modified Bessel function of 0-order
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call DefineSet( KSIO3(BoxNumberXY), DEFINE, 13, QUADRATIC_TERM, dummy, &
        dummy)

      call DefineSet( KSIO3(BoxNumberXY), DEFINE, 14, LINEAR_TERM, dummy, &
        dummy)

      call DefineSet( KSIO3(BoxNumberXY), DEFINE, 15, CONSTANT_TERM, dummy, &
        dummy)


      call DefineSet( KSIO3(BoxNumberXY), PARAMETER_DEFINE, 21, &
        BESSELI_EXP_TERM, - alpha, suD1)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call CompleteSet( KSIO3(BoxNumberXY), SET_CONTINUITY, STANDARD, idummy, &
        dummy, dummy)

      call CompleteSet( KSIO3(BoxNumberXY), SET_BOUNDARY, LAYER1, &
        EQUATION, 0.0D+00, chM5s- N5s_Ben(BoxNumberXY))

      call CompleteSet( KSIO3(BoxNumberXY), INPUT_TERM, 13, PARAMETER, dummy, &
        - zuBT)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate for the above defined set of boundary conditions
      ! the steady-state profiles and return the average concentration
      ! in the oxic and denitrification layers
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cM5s = CalculateSet( KSIO3(BoxNumberXY), SET_LAYER_INTEGRAL_UNTIL, LAYER1, &
        LAYER2, cD2m, 0.0D+00)/ cD2m

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the adaptation time to the steady-state profile
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Tau  =   CalculateTau(  0.0D+00,  diff,  p_p,  cD2m)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate the average value of M5s over the actual time step
      ! (transient value).
      ! This value depends on the adaptation time, the actual time step,
      ! the ''old'' value and the ''equilibrium value''
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cM5s = cM5s+( chM5s- M5s(BoxNumberXY)- cM5s)* IntegralExp( - LocalDelta/ &
        Tau, 1.0D+00)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! 1.Store equilibrium profile
      ! 2.Derive the equations for the transient profiles, assuming the same
      !  solution as for the steady-state case and using cM5s as new constraint.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      KSIO3E(BoxNumberXY) = CopySet( KSIO3(BoxNumberXY), &
        KSIO3E(BoxNumberXY))
      dummy = CalculateSet( KSIO3(BoxNumberXY), ADD, 0, 0, dummy, cD2m* &
        cM5s)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate the pore-water average concentrations for the standard &
      ! ''D2.n''
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M5s(BoxNumberXY) = chM5s- CalculateFromSet( KSIO3E(BoxNumberXY), &
        INTEGRAL, STANDARD, 0.0D+00, D2m(BoxNumberXY))/ D2m(BoxNumberXY)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Start calculation of fluxes:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate flux at the sediment/water interface:
      ! Flux limitation at very low values of N5s
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jK5N5s(BoxNumberXY) = - CalculateFromSet( KSIO3(BoxNumberXY), DERIVATIVE, &
        RFLUX, 0.0D+00, dummy)
      jK5N5s(BoxNumberXY) = jK5N5s(BoxNumberXY)* insw( ( &
        M5s(BoxNumberXY)- N5s_Ben(BoxNumberXY))* jK5N5s(BoxNumberXY))




      call flux(BoxNumberXY, iiBen, ppK5s, ppK5s, -( jK5N5s(BoxNumberXY)) )


      jK15K5s = - CalculateFromSet( KSIO3(BoxNumberXY), DERIVATIVE, RFLUX, &
        cD2m, dummy)
      call flux(BoxNumberXY, iiBen, ppK5s, ppK5s, -(- jK15K5s) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate new depth of the sulphide horizon
      ! and the flux of silicate related to this shifting
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if ( abs(shiftD2m(BoxNumberXY))> 0.0D+00) then

        cD2mNew = min( p_chD2m, max( p_clD2m, &
          D2m(BoxNumberXY)+ shiftD2m(BoxNumberXY)* LocalDelta))

        cShiftD2m  =   cD2mNew- cD2m

        if ( abs(cShiftD2m)> 0.0D+00) then
          s = chM5s* cShiftD2m* p_poro(BoxNumberXY)*( 1.0D+00+ &
            p_p)- CalculateFromSet( KSIO3(BoxNumberXY), INTEGRAL, MASS, cD2m, &
            cD2mNew)
          if ( cShiftD2m< 0.0D+00) then
            !        s=max(s, cShiftD2m*M5s*p_poro*(1.0+p_p));
            s  =   max(  s, - LocalDelta* K5s(BoxNumberXY))
          end if

          ! recalculation to rates per day....
          call flux(BoxNumberXY, iiBen, ppK5s, ppK5s, -(- s/ LocalDelta) )
          jK15K5s  =   jK15K5s+ s/ LocalDelta
        end if

      end if



      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! the dissolution fluxes:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jQ6K15s = max( 0.0D+00, suD1* CalculateFromSet( &
        KSIO3E(BoxNumberXY), EXPONENTIAL_INTEGRAL, RFLUX, cD2m, p_d_tot))
      call flux(BoxNumberXY, iiBen, ppQ6s, ppQ6s, -( jQ6K15s) )

      jQ6K5s = max( 0.0D+00, suD1* CalculateFromSet( &
        KSIO3E(BoxNumberXY), EXPONENTIAL_INTEGRAL, RFLUX, cD1m, cD2m))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Determine dissolution rate in oxidized layer (mMol/m3).
      ! Zero order process:
      ! Maximalization: this important source can not cause higher values
      ! than equilibrium flux of jQ6M5s is limited in such a way that M5s can
      ! never reach a value higher than the chM5s:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      rmQ6s = - GetInfoFromSet( KSIO3E(BoxNumberXY), INTEGRAL, PARAMETER, &
        13, 0.0D+00, cD1m)

      s = max( 0.0D+00, min( ( chM5s- M5s(BoxNumberXY))* &
        cD1m* p_poro(BoxNumberXY)*( 1.0D+00+ p_p)+ jK5N5s(BoxNumberXY)- jK15K5s- &
        jQ6K5s, rmQ6s))

      call flux(BoxNumberXY, iiBen, ppQ6s, ppK5s, jQ6K5s+ s )


      call flux(BoxNumberXY, iiBen, ppD9m, ppD9m, ( 0.5D+00* cD2m- &
        D9m(BoxNumberXY))*( jQ6K5s)/( 1.0D-80+ Q6s(BoxNumberXY)) )
      call flux(BoxNumberXY, iiBen, ppD9m, ppD9m, ( 0.5D+00*( &
        p_d_tot- cD2m)- D9m(BoxNumberXY))*( jQ6K15s)/( 1.0D-80+ Q6s(BoxNumberXY)) &
        )


    end DO


  end DO


  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
