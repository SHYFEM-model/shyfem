#include "DEBUG.h"
#include "INCLUDE.h"

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcOxygenSaturation_3
!
! DESCRIPTION
!   !       Function for calculation of the Oxygen Saturation
!
!	based on:
!
!	WEISS 1970 DEEP SEA RES 17, 721-735.
!	units of ln(ml(STP)/l)
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  SUBROUTINE CalcOxygenSaturation()
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): O2o
  ! The following Pelagic 1-d global boxvars are modified : cxoO2
  ! The following Pelagic 1-d global boxvars got a value: eO2mO2
  ! The following Pelagic 1-d global boxvars  are used: ETW, ESW
  ! The following global constants are used: RLEN,ZERO_KELVIN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO_KELVIN
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem,  ONLY: O2o
#endif
  use mem, ONLY: ppO2o, cxoO2, eO2mO2, ETW, ESW, NO_BOXES, iiBen, iiPel, &
    flux_vector



!  
!
! !AUTHORS
!   Piet Ruardij
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
  real(RLEN),dimension(NO_BOXES)  :: h
  real(RLEN),dimension(NO_BOXES)  :: abt

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  ! calc absolute temperature divided by 100.0;
  ! input for the next empirical equation.

  abt  =  ( ETW(:)- ZERO_KELVIN)/ 100.0D+00

  !   calc theoretical oxygen saturation for temp + salinity
  !   From WEISS 1970 DEEP SEA RES 17, 721-735.
  !   units of ln(ml(STP)/l)

  h = - 173.4292D+00+ 249.6339D+00/ abt+ 143.3483D+00* log( &
    abt)- 21.8492D+00* abt+ ESW(:)*(- 0.033096D+00+ 0.014259D+00* abt- &
    0.0017D+00* (abt)**(2.0D+00))

  ! convert units to ml(STP)/l

  h  =   exp(  h)

  ! convert to mMol/m3

  !
  !   calc volume of an ideal gas at standard temp (25C) and
  !   pressure (1.e-3 atm)
  !   p_videal = (8.3145 * 298.15 / 101325.0) = 24.4665e-3;
  !

  cxoO2(:)  =   h/ 24.4665D-3
  eO2mO2(:)  =   O2o(:)/ cxoO2(:)






  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
