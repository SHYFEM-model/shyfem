#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenthicSystem
!
! DESCRIPTION
!   !   This is the top level of the benthic submodel.
!       All the biological processes affecting the benthic dynamics are called
!       in a specified sequence according to the calculation flags.
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenthicSystemDynamics
!
! !USES:
  ! The following Benthic-states are used (NOT in fluxes): Y1c, Y1n, Y1p, Y2c, &
  ! Y2n, Y2p, Y4c, Y4n, Y4p, Y5c, Y5n, Y5p, H1c, H1n, H1p, H2c, H2n, H2p
  ! The following groupmember vars  are used: iiY1, iiY2, iiY4, iiY5, iiH1, iiH2
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem, ONLY: ppY1c, ppY1n, ppY1p, ppY2c, ppY2n, ppY2p, ppY4c, ppY4n, &
    ppY4p, ppY5c, ppY5n, ppY5p, ppH1c, ppH1n, ppH1p, ppH2c, ppH2n, ppH2p, iiY1, &
    iiY2, iiY3, iiY4, iiY5, iiH1, iiH2, iiBen, iiPel, flux
  use mem_Param, ONLY: CalcBenOrganisms, CalcBenBacteria

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following group processes are used: BenOrganismDynamics, BenBacDynamics
  use global_interface, ONLY: BenOrganismDynamics, BenBacDynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



!
!
! !AUTHORS
!   ERSEM group
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


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of turenh and irrenh in Bioturbation...
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call BioturbationDynamics

  if ( CalcBenOrganisms(iiY1)) then
    call BenOrganismDynamics( iiY1, ppY1c, ppY1n, ppY1p)
  end if

  if ( CalcBenOrganisms(iiY2)) then
    call BenOrganismDynamics( iiY2, ppY2c, ppY2n, ppY2p)
  end if

  if ( CalcBenOrganisms(iiY4)) then
    call BenOrganismDynamics( iiY4, ppY4c, ppY4n, ppY4p)
  end if

  if ( CalcBenOrganisms(iiY5)) then
    call BenOrganismDynamics( iiY5, ppY5c, ppY5n, ppY5p)
  end if

  if ( CalcBenOrganisms(iiY3)) then
    call FilterFeederDynamics
  end if

  if ( CalcBenBacteria(iiH1)) then
    call BenBacDynamics( iiH1, ppH1c, ppH1n, ppH1p)
  end if

  if ( CalcBenBacteria(iiH2)) then
    call BenBacDynamics( iiH2, ppH2c, ppH2n, ppH2p)
  end if

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
