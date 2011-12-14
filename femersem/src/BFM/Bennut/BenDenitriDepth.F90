#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenDenitriDepth
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenDenitriDepthDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: D2m
  ! The following Benthic-states are used (NOT in fluxes): D1m
  ! The following global scalar vars are used: BoxNumberZ, &
  ! NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, BoxNumber, &
  ! BoxNumberXY, dummy
  ! The following Benthic 1-d global boxvars are modified : shiftD2m
  ! The following Benthic 1-d global boxvars  are used: KNO3, KNH4
  ! The following 0-d global box parametes are used: p_d_tot, p_clD1D2m
  ! The following global constants are used: RLEN
  ! The following constants are used: EQUATION, STANDARD, GET, LABDA_1, &
  ! ONE_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: D2m, D1m, D2STATE
  use mem, ONLY: ppD2m, ppD1m, BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, &
    BoxNumberY, NO_BOXES_Y, BoxNumber, BoxNumberXY, dummy, shiftD2m, KNO3, KNH4, &
    iiBen, iiPel, flux
  use constants,  ONLY: EQUATION, STANDARD, GET, LABDA_1, ONE_PER_DAY
  use mem_Param,  ONLY: p_d_tot, p_clD1D2m
  use mem_BenDenitriDepth


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:CalculateFromSet, GetInfoFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface,   ONLY: CalculateFromSet, GetInfoFromSet



!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
!
! !REVISION_HISTORY
!   !
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
  integer  :: control
  real(RLEN)  :: pmM3n
  real(RLEN)  :: ushiftD2m
  real(RLEN)  :: D2mNew
  real(RLEN)  :: M3n_D1m

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! user defined external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, external  :: D3toD1
  integer, external  :: D2toD1
  integer, external  :: PrintSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  BoxNumberZ = 1
  DO BoxNumberY=1,NO_BOXES_Y
    DO BoxNumberX=1,NO_BOXES_X
      BoxNumber=D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ)
      BoxNumberXY=D2toD1(BoxNumberX,BoxNumberY)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate concentration of nitrate in porewater in M3n:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M3n_D1m = CalculateFromSet( KNO3(BoxNumberXY), EQUATION, &
        STANDARD, D1m(BoxNumberXY), dummy)

      select case ( M3n_D1m< 0.0D+00)

        case( .TRUE. )

          control  =   PrintSet(  KNH4(BoxNumberXY))
          control  =   PrintSet(  KNO3(BoxNumberXY))




        case( .FALSE. )

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Calculate fraction with which new depth of D2.m is calculated:
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          pmM3n  =   max(  p_pmM3n* M3n_D1m,  p_clM3n_D2)/ M3n_D1m

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! According solution nitrate concentrate decreases
          ! exponentially. Use exponent to calculate depth at which
          ! concentration is pmM3n * M3n_D1m. Use this new depth as
          ! uncorrected new denitrification depth
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          D2mNew = min( p_cmD2m, ( log( pmM3n)/ &
            GetInfoFromSet( KNO3(BoxNumberXY), GET, LABDA_1, 21, dummy, &
            dummy)+ D1m(BoxNumberXY)))

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! 1. Calculate uncorrected shift of D2.m
          ! 2. limit shift incase of D2mnew moves in the direction of D1m
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          ushiftD2m = max( min( p_d_tot- p_clD1D2m, &
            D2mNew), D1m(BoxNumberXY)+ p_clD1D2m)- D2m(BoxNumberXY)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Correct by damping the change of D2m in case large changes:
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          shiftD2m(BoxNumberXY) = ushiftD2m/ ONE_PER_DAY* &
            (D2m(BoxNumberXY)/( D2m(BoxNumberXY)+ abs(ushiftD2m)))**(p_xdamping)

          call flux(BoxNumberXY, iiBen, ppD2m, ppD2m, shiftD2m(BoxNumberXY) )




      end select


    end DO


  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
