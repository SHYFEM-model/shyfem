#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelGlobal
!
! DESCRIPTION
!   !
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PelGlobalDynamics
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): R6p, R6c, R6n, R6s, &
  ! P1s, P1c, B1p, B1c, B1n
  ! The following box states are used (NOT in fluxes): &
  ! MicroZooPlankton, MesoZooPlankton, PhytoPlankton
  ! The following Pelagic 1-d global boxvars got a value: flP1R6s, flPTN6r, &
  ! qpR6c, qnR6c, qsR6c, qpB1c, qnB1c, sediR6
  ! The following Pelagic 2-d global boxvars got a value: qp_mz, qn_mz, qpZc, &
  ! qnZc, qpPc, qnPc, qlPc, qsPc, sediPI
  ! The following groupmember vars are used: iiMicroZooPlankton, &
  ! iiMesoZooPlankton, iiPhytoPlankton, iiP1
  ! The following constituent constants  are used: iiP, iiC, iiN, iiL
  ! The following 0-d global box parametes are used: p_small
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem, ONLY: R6p, R6c, R6n, R6s, P1s, P1c, B1p, &
    B1c, B1n, MicroZooPlankton, MesoZooPlankton, PhytoPlankton
#endif
  use mem, ONLY: ppR6p, ppR6c, ppR6n, ppR6s, ppP1s, ppP1c, &
    ppB1p, ppB1c, ppB1n, ppMicroZooPlankton, ppMesoZooPlankton, ppPhytoPlankton, &
    flP1R6s, flPTN6r, qpR6c, qnR6c, qsR6c, qpB1c, qnB1c, sediR6, qp_mz, &
    qn_mz, qpZc, qnZc, qpPc, qnPc, qlPc, qsPc, sediPI, iiMicroZooPlankton, &
    iiMesoZooPlankton, iiPhytoPlankton, iiP1, iiP, iiC, iiN, iiL, NO_BOXES, &
    iiBen, iiPel, flux_vector
  use mem_Param,  ONLY: p_small
  use mem_PelGlobal



!  
!
! !AUTHORS
!   Piet Ruardij
!
!
!
! !REVISION_HISTORY
!   Created at Tue Apr 20 09:11:59 AM CEST 2004
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

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Reset var in which silica fluxes are collected:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  flP1R6s(:)  =   0.0D+00
  flPTN6r(:)  =   0.0D+00
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in pelagic detritus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  qpR6c(:)  =   R6p(:)/( p_small+ R6c(:))
  qnR6c(:)  =   R6n(:)/( p_small+ R6c(:))
  qsR6c(:)  =   R6s(:)/( p_small+ R6c(:))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in microzooplankton and HNAN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , ( iiMicroZooPlankton)

    qp_mz(i,:)  =   MicroZooPlankton(i,iiP)/( p_small+ MicroZooPlankton(i,iiC))
    qn_mz(i,:)  =   MicroZooPlankton(i,iiN)/( p_small+ MicroZooPlankton(i,iiC))
  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in omnivorous and herbivorous mesozooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , ( iiMesoZooPlankton)

    qpZc(i,:)  =   MesoZooPlankton(i,iiP)/( p_small+ MesoZooPlankton(i,iiC))
    qnZc(i,:)  =   MesoZooPlankton(i,iiN)/( p_small+ MesoZooPlankton(i,iiC))
  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in phytoplankton
  ! Compute light prop.or chl. quota in phytoplankton (dep. on ChlLightFlag)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , ( iiPhytoPlankton)

    qpPc(i,:)  =   PhytoPlankton(i,iiP)/( p_small+ PhytoPlankton(i,iiC))
    qnPc(i,:)  =   PhytoPlankton(i,iiN)/( p_small+ PhytoPlankton(i,iiC))
    qlPc(i,:)  =   PhytoPlankton(i,iiL)/( p_small+ PhytoPlankton(i,iiC))
  end do

  qsPc(iiP1,:)  =   P1s(:)/( p_small+ P1c(:))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in Pelagic Bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  qpB1c(:)  =   B1p(:)/( p_small+ B1c(:))
  qnB1c(:)  =   B1n(:)/( p_small+ B1c(:))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute sedimentation velocities
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sediR6(:)  =   p_rR6m
  do i = 1 , ( iiPhytoPlankton)

    sediPI(i,:)  =   p_rPIm( i)
  end do






  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
