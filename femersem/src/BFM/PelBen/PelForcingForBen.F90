#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelForcingForBen
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PelForcingForBenDynamics
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): R6c, R6n, R6p, R6s, &
  ! N1p, N3n, N4n, N5s, N6r, O2o
  ! The following box states are used (NOT in fluxes): PhytoPlankton
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY
  ! The following Pelagic 1-d global boxvars  are used: ETW, Depth
  ! The following Benthic 1-d global boxvars are modified : PIc, PIn, PIp, PIs
  ! The following Benthic 1-d global boxvars got a value: RIc, &
  ! RIn, RIp, RIs, N1p_Ben, N3n_Ben, N4n_Ben, N5s_Ben, N6r_Ben, O2o_Ben, &
  ! ETW_Ben, Depth_Ben
  ! The following groupmember vars  are used: iiPhytoPlankton, iiP1
  ! The following constituent constants  are used: iiC, iiN, iiP, iiS
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem, ONLY: R6c, R6n, R6p, R6s, N1p, N3n, N4n, N5s, N6r, O2o, &
    PhytoPlankton, D2STATE
  use mem, ONLY: ppR6c, ppR6n, ppR6p, ppR6s, ppN1p, ppN3n, &
    ppN4n, ppN5s, ppN6r, ppO2o, ppPhytoPlankton, BoxNumberZ, NO_BOXES_Z, &
    BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, BoxNumber, BoxNumberXY, ETW, &
    Depth, PIc, PIn, PIp, PIs, RIc, RIn, RIp, RIs, N1p_Ben, N3n_Ben, &
    N4n_Ben, N5s_Ben, N6r_Ben, O2o_Ben, ETW_Ben, Depth_Ben, iiPhytoPlankton, iiP1, &
    iiC, iiN, iiP, iiS, iiBen, iiPel, flux



!  
!
! !AUTHORS
!   Piet Ruardij
!
!
!
! !REVISION_HISTORY
!   Created at Wed Jun 16 02:04:44 PM CEST 2004
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
  ! Local Vectors used  of group vectors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN), dimension(:), pointer  ::lcl_PhytoPlankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i

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


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Compute total phytoplankton conc. used as food for filtereeders
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      PIc(BoxNumberXY)  =   0.0D+00
      PIn(BoxNumberXY)  =   0.0D+00
      PIp(BoxNumberXY)  =   0.0D+00
      PIs(BoxNumberXY)  =   0.0D+00

      do i = 1 , ( iiPhytoPlankton)

        lcl_PhytoPlankton => PhytoPlankton(i,iiC)
        PIc(BoxNumberXY)  =   PIc(BoxNumberXY)+ lcl_PhytoPlankton(BoxNumber)
        lcl_PhytoPlankton => PhytoPlankton(i,iiN)
        PIn(BoxNumberXY)  =   PIn(BoxNumberXY)+ lcl_PhytoPlankton(BoxNumber)
        lcl_PhytoPlankton => PhytoPlankton(i,iiP)
        PIp(BoxNumberXY)  =   PIp(BoxNumberXY)+ lcl_PhytoPlankton(BoxNumber)
        if ( i== iiP1) then
          lcl_PhytoPlankton => PhytoPlankton(i,iiS)
          PIs(BoxNumberXY)  =   PIs(BoxNumberXY)+ lcl_PhytoPlankton(BoxNumber)
        end if

      end do


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Compute total detritus conc. used as food for filtereeders
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      RIc(BoxNumberXY)  =   R6c(BoxNumber)
      RIn(BoxNumberXY)  =   R6n(BoxNumber)
      RIp(BoxNumberXY)  =   R6p(BoxNumber)
      RIs(BoxNumberXY)  =   R6s(BoxNumber)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Derive Forcing for benthos
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !Nutrient Forcing:
      N1p_Ben(BoxNumberXY)  =   N1p(BoxNumber)
      N3n_Ben(BoxNumberXY)  =   N3n(BoxNumber)
      N4n_Ben(BoxNumberXY)  =   N4n(BoxNumber)
      N5s_Ben(BoxNumberXY)  =   N5s(BoxNumber)
      N6r_Ben(BoxNumberXY)  =   N6r(BoxNumber)

      !Oxygen Forcing:
      O2o_Ben(BoxNumberXY)  =   O2o(BoxNumber)

      ! Temperature in the benthos is made equal tothe temperature of the
      ! adjacent level (layer) of the pelagic
      ETW_Ben(BoxNumberXY)  =   ETW(BoxNumber)

      ! depth of the level aboce the sediment
      Depth_Ben(BoxNumberXY)  =   Depth(BoxNumber)

    end DO



  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
