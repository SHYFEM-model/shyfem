#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenNitrogenShifting
!
! DESCRIPTION
!   Description of shifting of dissolved N (amm and nitrate) between
!       layers
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenNitrogenShiftingDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: K14n, K4n, K24n, K3n
  ! The following Benthic-states are used (NOT in fluxes): D1m, D7m, D2m
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY, LocalDelta, dummy
  ! The following Benthic 1-d global boxvars are used: shiftD1m, KNH4, reATn, &
  ! shiftD2m, KNO3
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global box parametes are used: p_clDxm, p_d_tot
  ! The following global constants are used: RLEN
  ! The following constants are used: SHIFT, LAYER1, DERIVATIVE, RFLUX, LAYER2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: K14n, K4n, K24n, K3n, D1m, D7m, D2m, D2STATE
  use mem, ONLY: ppK14n, ppK4n, ppK24n, ppK3n, ppD1m, ppD7m, &
    ppD2m, BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
    BoxNumber, BoxNumberXY, LocalDelta, dummy, shiftD1m, KNH4, reATn, shiftD2m, &
    KNO3, iiBen, iiPel, flux
  use constants,  ONLY: SHIFT, LAYER1, DERIVATIVE, RFLUX, LAYER2
  use mem_Param,  ONLY: p_poro, p_clDxm, p_d_tot


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface,   ONLY: CalculateFromSet


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:insw, IntegralExp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: insw, IntegralExp



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
  real(RLEN)  :: Dnew
  real(RLEN)  :: shiftmass
  real(RLEN)  :: zuD1
  real(RLEN)  :: zuD2
  real(RLEN)  :: jK14K4n
  real(RLEN)  :: jK24K14n
  real(RLEN)  :: jK34K24n
  real(RLEN)  :: jK13K3n
  real(RLEN)  :: alpha

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


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Ammonium Fluxes at the oxic/denitrification boundary
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Dnew  =   D1m(BoxNumberXY)+ LocalDelta* shiftD1m(BoxNumberXY)

      ! Calculate mass shifted in upwards direction:
      shiftmass = CalculateFromSet( KNH4(BoxNumberXY), SHIFT, LAYER1, &
        D1m(BoxNumberXY), Dnew)/ LocalDelta

      jK14K4n = CalculateFromSet( KNH4(BoxNumberXY), DERIVATIVE, RFLUX, &
        D1m(BoxNumberXY), 0.0D+00)+ shiftmass

      call flux(BoxNumberXY, iiBen, ppK14n, ppK4n, jK14K4n* insw( jK14K4n) )
      call flux(BoxNumberXY, iiBen, ppK4n, ppK14n, - jK14K4n* insw( - jK14K4n) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! All the nutrient mineralization source term in the anoxic layer
      ! has been added to K14.n in BenBacDynamics
      ! However in the model this layer is subdivided and hence a partition
      ! flux is here calculated according to the exponential distribution.
      !
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D7.m is the average penetration depth for N-detritus
      !                 +
      ! Anoxic Mineralization at D1.m, using the exponential distribution
      !                 +
      !          Anoxic Mineralization at D2.m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      alpha  =   1.0D+00/ max(  p_clDxm,  D7m(BoxNumberXY))
      zuD1 = max( 1.D-20, reATn(BoxNumberXY))/ p_poro(BoxNumberXY)/ IntegralExp( &
        - alpha, p_d_tot- D1m(BoxNumberXY))
      zuD2  =   zuD1* exp( - alpha*( D2m(BoxNumberXY)- D1m(BoxNumberXY)))

      jK24K14n = - zuD2* p_poro(BoxNumberXY)* IntegralExp( - alpha, &
        p_d_tot- D2m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Ammonium Fluxes at the denitrification/anoxic boundary
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Dnew  =   D2m(BoxNumberXY)+ LocalDelta* shiftD2m(BoxNumberXY)

      ! Calculate mass shifted in upwards direction:
      shiftmass = CalculateFromSet( KNH4(BoxNumberXY), SHIFT, LAYER2, &
        D2m(BoxNumberXY), Dnew)/ LocalDelta

      jK24K14n = jK24K14n+ CalculateFromSet( KNH4(BoxNumberXY), DERIVATIVE, &
        RFLUX, D2m(BoxNumberXY), 0.0D+00)+ shiftmass

      call flux(BoxNumberXY, iiBen, ppK24n, ppK14n, jK24K14n* insw( jK24K14n) )
      call flux(BoxNumberXY, iiBen, ppK14n, ppK24n, - jK24K14n* insw( - &
        jK24K14n) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Fluxes at the lower boundary
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Ammonium:
      jK34K24n = CalculateFromSet( KNH4(BoxNumberXY), DERIVATIVE, RFLUX, &
        p_d_tot, 0.0D+00)
      call flux(BoxNumberXY, iiBen, ppK24n, ppK24n, -(- jK34K24n) )

      ! Nitrate:
      shiftmass = CalculateFromSet( KNO3(BoxNumberXY), SHIFT, LAYER2, &
        D2m(BoxNumberXY), Dnew)/ LocalDelta
      jK13K3n = CalculateFromSet( KNO3(BoxNumberXY), DERIVATIVE, RFLUX, &
        D2m(BoxNumberXY), dummy)+ shiftmass

      call flux(BoxNumberXY, iiBen, ppK3n, ppK3n, -(- jK13K3n) )


    end DO


  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
