#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelagicSystem
!
! DESCRIPTION
!   This is the Pelagic Submodel. 
!   All the pelagic biogeochemical modules are called in sequence
!   according to the logical switches
!        
!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PelagicSystemDynamics
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): P1p, P2p, P3p, P4p, &
  ! B1p, Z3p, Z4p, Z5p, Z6p, R1p, R6p, N1p, P1n, P2n, P3n, P4n, B1n, Z3n, Z4n, &
  ! Z5n, Z6n, R1n, R6n, N3n, N4n, O4n, P1s, R6s, N5s
  ! The following Pelagic-states are used (NOT in fluxes): P1c, P1n, P1p, P1s, &
  ! P1l, P2c, P2n, P2p, P2s, P2l, P3c, P3n, P3p, P3s, P3l, P4c, P4n, P4p, P4s, &
  ! P4l, Z3c, Z3n, Z3p, Z4c, Z4n, Z4p, Z5c, Z5n, Z5p, Z6c, Z6n, Z6p
  ! The following Pelagic 1-d global boxvars  are used: Depth
  ! The following Benthic 1-d global boxvars are modified : totpeln, totpels
  ! The following Benthic 1-d global boxvars  are used: totpelp
  ! The following groupmember vars are used: iiP1, iiP2, iiP3, iiP4, iiZ3, iiZ4, &
  ! iiZ5, iiZ6
  ! The following 0-d global box parametes are used: &
  ! ChlLightFlag
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem, ONLY: P1p, P2p, P3p, P4p, B1p, Z3p, Z4p, Z5p, Z6p, R1p, R6p, N1p, &
    P1n, P2n, P3n, P4n, B1n, Z3n, Z4n, Z5n, Z6n, R1n, R6n, N3n, N4n, O4n, P1s, &
    R6s, N5s, D2STATE
  use mem, ONLY: ppP1p, ppP2p, ppP3p, ppP4p, ppB1p, ppZ3p, ppZ4p, ppZ5p, &
    ppZ6p, ppR1p, ppR6p, ppN1p, ppP1n, ppP2n, ppP3n, ppP4n, ppB1n, ppZ3n, ppZ4n, &
    ppZ5n, ppZ6n, ppR1n, ppR6n, ppN3n, ppN4n, ppO4n, ppP1s, ppR6s, ppN5s, ppP1c, &
    ppP1l, ppP2c, ppP2l, ppP3c, ppP3l, ppP4c,  ppP4l, ppZ3c,ppZ4c, ppZ5c, ppZ6c, &
    ppP2s, ppP3s, ppP4s,  &
    Depth, totpeln, totpels, totpelp, iiP1, iiP2, iiP3, iiP4,&
    iiZ3, iiZ4, iiZ5, iiZ6, iiBen, iiPel, flux 
   
  use mem_Param, ONLY: ChlLightFlag, CalcPhytoPlankton,CalcMicroZooPlankton, &
    CalcMesoZooPlankton, CalcBacteria, CalcPelChemistry

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:CalcChlorophylla, &
  ! CalcOxygenSaturation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: CalcChlorophylla, CalcOxygenSaturation

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following group processes are &
  ! used: PhotoAvailableRadiationDynamics, PhytoDynamics, &
  ! LightAdaptationDynamics, MesoZooDynamics, MicroZooDynamics
  use global_interface, ONLY: PhotoAvailableRadiationDynamics, &
    PhytoDynamics, LightAdaptationDynamics, MesoZooDynamics, MicroZooDynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



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

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  call  CalcChlorophylla( )


  call PelGlobalDynamics

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute oxygen variables: cxoO2 eO2mO2
  ! calculate oxygen OxygenReaeration
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call  CalcOxygenSaturation( )


  call OxygenReaerationDynamics

  if ( ChlLightFlag== 1) then
    if ( CalcPhytoPlankton(iiP1)) then
      call PhotoAvailableRadiationDynamics( iiP1, ppP1c, ppP1n, ppP1p, ppP1s, &
        ppP1l)
    end if

    if ( CalcPhytoPlankton(iiP2)) then
      call PhotoAvailableRadiationDynamics( iiP2, ppP2c, ppP2n, ppP2p, ppP2s, &
        ppP2l)
    end if

    if ( CalcPhytoPlankton(iiP3)) then
      call PhotoAvailableRadiationDynamics( iiP3, ppP3c, ppP3n, ppP3p, ppP3s, &
        ppP3l)
    end if

    if ( CalcPhytoPlankton(iiP4)) then
      call PhotoAvailableRadiationDynamics( iiP4, ppP4c, ppP4n, ppP4p, ppP4s, &
        ppP4l)
    end if

  end if




  if ( CalcPhytoPlankton(iiP1)) then
    call PhytoDynamics( iiP1, ppP1c, ppP1n, ppP1p, ppP1s, ppP1l)
  end if

  if ( CalcPhytoPlankton(iiP2)) then
    call PhytoDynamics( iiP2, ppP2c, ppP2n, ppP2p, ppP2s, ppP2l)
  end if

  if ( CalcPhytoPlankton(iiP3)) then
    call PhytoDynamics( iiP3, ppP3c, ppP3n, ppP3p, ppP3s, ppP3l)
  end if

  if ( CalcPhytoPlankton(iiP4)) then
    call PhytoDynamics( iiP4, ppP4c, ppP4n, ppP4p, ppP4s, ppP4l)
  end if


  if ( ChlLightFlag== 1) then
    if ( CalcPhytoPlankton(iiP1)) then
      call LightAdaptationDynamics( iiP1, ppP1c, ppP1n, ppP1p, ppP1s, ppP1l)
    end if

    if ( CalcPhytoPlankton(iiP2)) then
      call LightAdaptationDynamics( iiP2, ppP2c, ppP2n, ppP2p, ppP2s, ppP2l)
    end if

    if ( CalcPhytoPlankton(iiP3)) then
      call LightAdaptationDynamics( iiP3, ppP3c, ppP3n, ppP3p, ppP3s, ppP3l)
    end if

    if ( CalcPhytoPlankton(iiP4)) then
      call LightAdaptationDynamics( iiP4, ppP4c, ppP4n, ppP4p, ppP4s, ppP4l)
    end if

  end if
  
  if ( CalcMesoZooPlankton(iiZ3)) then
    call MesoZooDynamics( iiZ3, ppZ3c, ppZ3n, ppZ3p)
  end if
  if ( CalcMesoZooPlankton(iiZ3)) then
    call MesoZooDynamics( iiZ4, ppZ4c, ppZ4n, ppZ4p)
  end if
  if ( CalcMicroZooPlankton(iiZ5)) then
    call MicroZooDynamics( iiZ5, ppZ5c, ppZ5n, ppZ5p)
  end if
  if ( CalcMicroZooPlankton(iiZ6)) then
    call MicroZooDynamics( iiZ6, ppZ6c, ppZ6n, ppZ6p)
  end if
  if ( CalcBacteria) then
    call PelBacDynamics
  end if
  if ( CalcPelChemistry) then
    call PelChemDynamics
  end if

    totpelp(1) = sum((P1p+ P2p+ P3p+ P4p+ B1p+ Z3p+ Z4p+ Z5p+ Z6p+ R1p+ R6p+ &
      N1p)* Depth)
    totpeln(1) =sum( (P1n+ P2n+ P3n+ P4n+ B1n+ Z3n+ Z4n+ Z5n+ Z6n+ R1n+ R6n+ &
      N3n+ N4n+ O4n)* Depth)
    totpels(1)  = sum( (P1s+ R6s+ N5s)* Depth)



  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
