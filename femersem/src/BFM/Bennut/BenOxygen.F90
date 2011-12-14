#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenOxygen
!
! DESCRIPTION
!   Description of first order oxic processes in the sediment
!       and computation of oxygen penetration depth
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenOxygenDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: D1m, G2o, D2m
  ! The following global scalar vars are used: LocalDelta
  ! The following Benthic 1-d global boxvars are modified : shiftD1m, jG2O2o
  ! The following Benthic 1-d global boxvars are used: ETW_Ben, irrenh, &
  ! rrBTo, jG2K3o, jG2K7o, O2o_Ben
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global box parametes are used: p_small, p_d_tot, &
  ! CalcBenthicFlag
  ! The following global constants are used: RLEN
  ! The following constants are used: SEC_PER_DAY, ONE_PER_DAY, BENTHIC_BIO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem,  ONLY: D1m, G2o, D2m
#endif
  use mem, ONLY: ppD1m, ppG2o, ppD2m, LocalDelta, shiftD1m, jG2O2o, &
    ETW_Ben, irrenh, rrBTo, jG2K3o, jG2K7o, O2o_Ben, NO_BOXES_XY, iiBen, &
    iiPel, flux_vector
  use constants,  ONLY: SEC_PER_DAY, ONE_PER_DAY, BENTHIC_BIO
  use mem_Param,  ONLY: p_poro, p_small, p_d_tot, CalcBenthicFlag
  use mem_BenOxygen



!  
!
! !AUTHORS
!   P. Ruardij
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
  real(RLEN),dimension(NO_BOXES_XY)  :: r
  real(RLEN),dimension(NO_BOXES_XY)  :: diff
  real(RLEN),dimension(NO_BOXES_XY)  :: zmG2o
  real(RLEN),dimension(NO_BOXES_XY)  :: D1mNew
  real(RLEN),dimension(NO_BOXES_XY)  :: G2oNew
  real(RLEN),dimension(NO_BOXES_XY)  :: unc_shiftD1m
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Emperical equation derived from Broecker and Peng (1973)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  diff = SEC_PER_DAY* 1.0D-9* (10.0D+00)**((- 984.26D+00/( 273.0D+00+ &
    ETW_Ben(:))+ 3.672D+00))* p_poro* irrenh(:)* p_exsaf

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Recalculate total consumption from /m2 to /m3 pw:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  zmG2o  =  ( rrBTo(:)+ jG2K3o(:)+ jG2K7o(:))/( D1m(:)* p_poro)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Determine new thickness:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  D1mNew  =   sqrt(  2.0D+00* diff* O2o_Ben(:)/( p_small+ zmG2o))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate rate of change of thickness of the aerobic layer:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  unc_shiftD1m  =   max(  p_mD1m,  D1mNew)- D1m(:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Damping the change of D1m in case of large changes,
  ! in case of too thick D1m, and in case of too large negative changes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  shiftD1m(:) = unc_shiftD1m* (D1m(:)/( D1m(:)+ &
    abs(unc_shiftD1m)))**(p_xdampingD1m)*( p_chD1m/( p_chD1m+ D1m(:)))/ &
    ONE_PER_DAY
  r  =   max(  0.0D+00,  p_d_tot- p_chD1m- D1m(:))
  shiftD1m(:)  =   shiftD1m(:)* r/( r+ p_chD1m)

  call flux_vector( iiBen, ppD1m,ppD1m, shiftD1m(:) )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! recalculate the new D1mNew at the actual time step:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  D1mNew  =   D1m(:)+ shiftD1m(:)* LocalDelta

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! calculate the consumption which belongs to the corrected D1mNew
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  zmG2o  =  ( 2.0D+00* diff* O2o_Ben(:))/( D1mNew* D1mNew)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! New oxygen conc. in the sediment:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  G2oNew = D1mNew*( O2o_Ben(:)- 0.66667D+00* zmG2o* D1mNew* D1mNew/( 2.0D+00* &
    diff))* p_poro

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! flux to pelagic: correct flux for rate of change of G2o
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  jG2O2o(:) = -( rrBTo(:)+ jG2K3o(:)+ jG2K7o(:))-( G2oNew- G2o(:))/ &
    ONE_PER_DAY
  call flux_vector( iiBen, ppG2o,ppG2o,-( jG2O2o(:)) )

  if ( CalcBenthicFlag== BENTHIC_BIO) then
    ! Compute shifting of the denitrification layer here in case of running only
    ! the benthic submodel and NOT the benthic nutrient model.
    r  =   p_d_tot- D2m(:)
    call flux_vector( iiBen, ppD2m,ppD2m, shiftD1m(:)* r/( r+ 0.01D+00) )
  end if





  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
