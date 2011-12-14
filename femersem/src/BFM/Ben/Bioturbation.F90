#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Bioturbation
!
! DESCRIPTION
!   Describes sediment reworking by benthic organisms and their effect on
!   vertical transport of dissolved (bioirrigation) and particulate 
!   (bioturbation) matter
! 
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BioturbationDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: D6m, D7m, D8m, D9m
  ! The following Benthic-states are used (NOT in fluxes): Y2c, Y5c, Y1c, Y4c
  ! The following Benthic 1-d global boxvars are modified : turenh
  ! The following Benthic 1-d global boxvars got a value: irrenh, rrBTo, &
  ! reBTn, reBTp, rrATo, reATn, reATp
  ! The following Benthic 1-d global boxvars  are used: ETW_Ben
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem,  ONLY: D6m, D7m, D8m, D9m, Y2c, Y5c, Y1c, Y4c
#endif
  use mem, ONLY: ppD6m, ppD7m, ppD8m, ppD9m, ppY2c, ppY5c, ppY1c, &
    ppY4c, turenh, irrenh, rrBTo, reBTn, reBTp, rrATo, reATn, reATp, &
    ETW_Ben, NO_BOXES_XY, iiBen, iiPel, flux_vector
  use mem_Bioturbation


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:MM_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: MM_vector



!  
!
! !AUTHORS
!   W. Ebenhoeh and C. Kohlmeier.
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
  real(RLEN),dimension(NO_BOXES_XY)  :: et
  real(RLEN),dimension(NO_BOXES_XY)  :: Ytur
  real(RLEN),dimension(NO_BOXES_XY)  :: Yirr
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature Response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  et  =   (p_q10)**(( ETW_Ben(:)- 10.0D+00)* 0.1D+00)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of bioturbation. ''turenh'' is a global box variable!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  Ytur  =   Y2c(:)+ Y5c(:)+ p_turY1* Y1c(:)

  turenh(:)  =   1.0D+00+ p_cmtur* MM_vector(  Ytur,  p_chtur)* et

  call flux_vector(iiBen, ppD6m,ppD6m, p_Etur* turenh(:)*( 1.0D+00- exp( - &
    p_cturm/ D6m(:)))/ D6m(:))
  call flux_vector(iiBen, ppD7m,ppD7m, p_Etur* turenh(:)*( 1.0D+00- exp( - &
    p_cturm/ D7m(:)))/ D7m(:))
  call flux_vector(iiBen, ppD8m,ppD8m, p_Etur* turenh(:)*( 1.0D+00- exp( - &
    p_cturm/ D8m(:)))/ D8m(:))
  call flux_vector(iiBen, ppD9m,ppD9m, p_Etur* turenh(:)*( 1.0D+00- exp( - &
    p_cturm/ D9m(:)))/ D9m(:))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of bioirrigation. ''irrenh'' is a global box variable!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  Yirr  =   Y2c(:)+ Y5c(:)+ p_irrY4* Y4c(:)

  irrenh(:)  =   1.0D+00+ p_cmirr* MM_vector(  Yirr,  p_chirr)* et

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! In the following processes all respirations and excretions
  ! are added to the rr???? and re??? variables.
  ! There rates are input to the Benthic Nutrient model
  ! first these variables are initialized:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrBTo(:)  =   0.0D+00  ! mgO2/m2 # Total Benthic oxic respiration
  reBTn(:)  =   0.0D+00  ! mmN/m2  # Total Benthic oxic N mineralization
  reBTp(:)  =   0.0D+00  ! mmP/m2  # Total Benthic oxic P mineralization
  rrATo(:)  =   0.0D+00  ! mgO2/m2 # Total Benthic anoxic respiration
  reATn(:)  =   0.0D+00  ! mmN/m2  # Total Benthic anoxic N mineralization
  reATp(:)  =   0.0D+00  ! mmP/m2  # Total Benthic anoxic P mineralization







  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
