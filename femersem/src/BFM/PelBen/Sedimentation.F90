#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Sedimentation
!
! DESCRIPTION
!   Define all fluxes of material which arrive on the sediment:
!	a, fluxes of detritus (slow degradable and labile organic detritus)
!       b. Changes in distributions states which describes exponential distribution.
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine SedimentationDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: Q6c, Q6n, Q6p, Q6s, &
  ! Q1c, Q1n, Q1p, D6m, D7m, D8m, D9m
  ! The following Benthic 1-d global boxvars are used: rutQ6c, rutQ6n, &
  ! rutQ6p, rutQ6s, rutQ1c, rutQ1n, rutQ1p
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem,  ONLY: Q6c, Q6n, Q6p, Q6s, Q1c, Q1n, Q1p, D6m, D7m, D8m, D9m
#endif
  use mem, ONLY: ppQ6c, ppQ6n, ppQ6p, ppQ6s, ppQ1c, ppQ1n, ppQ1p, &
    ppD6m, ppD7m, ppD8m, ppD9m, rutQ6c, rutQ6n, rutQ6p, rutQ6s, rutQ1c, &
    rutQ1n, rutQ1p, NO_BOXES_XY, iiBen, iiPel, flux_vector



!  
!
! !AUTHORS
!   Piet Ruardij
!
!
!
! !REVISION_HISTORY
!   Created at Mon Nov 21 09:11:50 CET 2005
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

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculated Fluxes from Pelagic to Benthic
  ! These fluxes are the sum sedimentation flux + the flux of
  ! material excreted by filterfeeders and originating from the pelagic.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector( iiBen, ppQ6c,ppQ6c, rutQ6c(:) )
  call flux_vector( iiBen, ppQ6n,ppQ6n, rutQ6n(:) )
  call flux_vector( iiBen, ppQ6p,ppQ6p, rutQ6p(:) )
  call flux_vector( iiBen, ppQ6s,ppQ6s, rutQ6s(:) )

  call flux_vector( iiBen, ppQ1c,ppQ1c, rutQ1c(:) )
  call flux_vector( iiBen, ppQ1n,ppQ1n, rutQ1n(:) )
  call flux_vector( iiBen, ppQ1p,ppQ1p, rutQ1p(:) )

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of changes due to sedimentation of detritus in
  ! distribution state variables (Dx.m is a undetermined source).
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector(iiBen, ppD6m,ppD6m,( 0.0D+00- D6m(:))* rutQ6c(:)/ Q6c(:))
  call flux_vector(iiBen, ppD7m,ppD7m,( 0.0D+00- D7m(:))* rutQ6n(:)/ Q6n(:))
  call flux_vector(iiBen, ppD8m,ppD8m,( 0.0D+00- D8m(:))* rutQ6p(:)/ Q6p(:))
  call flux_vector(iiBen, ppD9m,ppD9m,( 0.0D+00- D9m(:))* rutQ6s(:)/ Q6s(:))





  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
