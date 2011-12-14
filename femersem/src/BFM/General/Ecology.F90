#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Ecology
!
! DESCRIPTION
!   !	This submodel calls all other submodels
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine EcologyDynamics
!
! !USES:
  ! The following 0-d global box parametes are used: CalcPelagicFlag, &
  ! CalcBenthicFlag
  ! The following global constants are used: RLEN
  ! The following constants are used: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: iiBen, iiPel, iiReset, flux
  use constants,  ONLY: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL
  use mem_Param,  ONLY: CalcPelagicFlag, CalcBenthicFlag


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:ResetTotMassVar
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: ResetTotMassVar



!  
!
! !AUTHORS
!   ERSEM team	
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

#ifdef DEBUG
  call  flux(1,iiReset,1,1,0.00D+00)
#endif

  call  ResetTotMassVar( )


  if ( CalcPelagicFlag) then

    call PelagicSystemDynamics

  end if


  call SettlingDynamics

  select case ( CalcBenthicFlag)

    case ( BENTHIC_RETURN )  ! Simple benthic return
      call BenthicReturn1Dynamics


    case ( BENTHIC_BIO )  ! Intermediate benthic return
      call PelForcingForBenDynamics
      call BenthicSystemDynamics
      call BenthicNutrient2Dynamics


    case ( BENTHIC_FULL )  ! Full benthic nutrients
      call PelForcingForBenDynamics
      call BenthicSystemDynamics
      call BenthicNutrient3Dynamics



  end select


  call ControlBenPartNutrientBuffersDynamics

  call BentoPelCoupDynamics

  call SedimentationDynamics

  call BenCheckMassConservationDynamics




  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
