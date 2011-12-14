#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: WindOxReaeration_3
!
! DESCRIPTION
!   Model describes reaeration between air and water column.
!       as forced by temperature and wind.
!
!       The equation and correlation used in this routine
!       are found in the 
!		R. Wanninkhof (1992), Relationship between windspeed and gas
!		exchange over the oecean
!               J. GeoPhys. Res. 97, 7373-7382
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine OxygenReaerationDynamics
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: O2o
  ! The following global scalar vars are used: BoxNumberZ, &
  ! BoxNumberY, NO_BOXES_Y, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumber, &
  ! BoxNumberXY, Wind
  ! The following Pelagic 1-d global boxvars  are used: ETW, cxoO2, Depth
  ! The following Benthic 1-d global boxvars are modified : jOAO2o
  ! The following 0-d global box parametes are used: AssignAirPelFluxesInBFMFlag
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: O2o, D3STATE
  use mem, ONLY: ppO2o, BoxNumberZ, BoxNumberY, NO_BOXES_Y, NO_BOXES_Z, &
    BoxNumberX, NO_BOXES_X, BoxNumber, BoxNumberXY, Wind, ETW, cxoO2, Depth, &
    jOAO2o, iiBen, iiPel, flux
  use mem_Param,  ONLY: AssignAirPelFluxesInBFMFlag
  use mem_WindOxReaeration_3


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:CalcSchmidtNumberOx
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: CalcSchmidtNumberOx



!  
!
! !AUTHORS
!   11 March 1998 Original version by P. Ruardij
!	              JWB 1999/03/25 Corrected k 
!
!
!
! !REVISION_HISTORY
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
  real(RLEN)  :: reacon
  real(RLEN)  :: p_schmidt

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! user defined external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, external  :: D3toD1
  integer, external  :: D2toD1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  BoxNumberZ = NO_BOXES_Z
  DO BoxNumberY=1,NO_BOXES_Y
    DO BoxNumberX=1,NO_BOXES_X
      BoxNumber=D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ)
      BoxNumberXY=D2toD1(BoxNumberX,BoxNumberY)


      !
      ! The authors assumed a Schimdt number of CO2 (=reference) of 660.0
      !

      p_schmidt  =   CalcSchmidtNumberOx(  ETW(BoxNumber))/ 660.0D+00

      !
      ! Calculate wind dependency:
      !`

      reacon  =   k* (Wind)**(2.0D+00)/ sqrt(  p_schmidt)

      jOAO2o(BoxNumberXY)  =   reacon*( cxoO2(BoxNumber)- O2o(BoxNumber))

      if ( AssignAirPelFluxesInBFMFlag) then
        call flux(BoxNumber, iiPel, ppO2o, ppO2o, jOAO2o(BoxNumberXY)/ &
          Depth(BoxNumber) )
      end if



    end DO


  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
