#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenthicReturn2
!
! DESCRIPTION
!   This process is an intermediate parameterisation of benthic 
!   remineralisation processes.
!   The sub-model of benthic organisms is required to use this 
!   parameterisation (at least benthic bacteria).
!   Oxygen demand at the water-sed interface is associated to amount 
!   of reduction equivalents present in the sediments.
!   Nutrients are released to the water column at constant specific
!   rates, according to the pore-water concentration.
!   Nitrogen remineralisation is partitioned 
!   into ammonium and nitrate flux with a constant value. 
!   A constant portion of biogenic silica in the sediments is
!   released to the water column as silicate.
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenthicReturn2Dynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: G2o, K6r, K1p, K11p, &
  ! K4n, K14n, Q6s
  ! The following Benthic-states are used (NOT in fluxes): D1m, D2m
  ! The following Benthic 1-d global boxvars are modified : jK3N3n, jK4N4n
  ! The following Benthic 1-d global boxvars got a value: jG2K7o, jK1N1p, jK5N5s
  ! The following 0-d global box parametes are used: p_qro, p_d_tot
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem,  ONLY: G2o, K6r, K1p, K11p, K4n, K14n, Q6s, D1m, D2m
#endif
  use mem, ONLY: ppG2o, ppK6r, ppK1p, ppK11p, ppK4n, ppK14n, ppQ6s, &
    ppD1m, ppD2m, jK3N3n, jK4N4n, jG2K7o, jK1N1p, jK5N5s, NO_BOXES_XY, iiBen, &
    iiPel, flux_vector
  use mem_Param,  ONLY: p_qro, p_d_tot
  use mem_BenthicReturn2



!  
!
! !AUTHORS
!   Piet Ruardij 
!
!
!
! !REVISION_HISTORY
!   Created at Fri Apr 30 21:30:27 CEST 2004
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
  real(RLEN),dimension(NO_BOXES_XY)  :: rate
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  ! Oxygen consumption in the sediments
  rate = p_reminO2* max( 0.0D+00, - G2o(:)/ D1m(:)+ K6r(:)/ p_qro/( &
    p_d_tot- D1m(:)))* D1m(:)
  call flux_vector( iiBen, ppG2o,ppG2o,-( rate) )
  call flux_vector( iiBen, ppK6r,ppK6r,-( rate* p_qro) )
  jG2K7o(:)  =   rate

  !----------------------------------------------------------------------
  ! Phosphorus remineralization in the sediments
  !----------------------------------------------------------------------
  rate  =   p_reminN1* K1p(:)/ D1m(:)
  ! jK1N1p is used in BenPelCoup to define the pelagic flux
  call flux_vector( iiBen, ppK1p,ppK1p,-( rate) )
  jK1N1p(:)  =   rate
  rate  =   p_K11K1p* K11p(:)/( D2m(:)- D1m(:))
  call flux_vector( iiBen, ppK11p,ppK1p, rate )
  ! K21.p is not used in this model version

  !----------------------------------------------------------------------
  ! Nitrogen remineralization in the sediments
  !----------------------------------------------------------------------
  rate  =   p_reminN4* K4n(:)
  ! K3.n is not used in this model version
  jK3N3n(:)  =   rate* p_pQIN3
  jK4N4n(:)  =   rate*( 1.0D+00- p_pQIN3)
  ! jK3N3n is used in BenPelCoup to define the pelagic flux
  ! jK4N4n is used in BenPelCoup to define the pelagic flux
  call flux_vector( iiBen, ppK4n,ppK4n,-( jK3N3n(:)+ jK4N4n(:)) )
  rate  =   p_K14K4n* K14n(:)/( D2m(:)- D1m(:))
  call flux_vector( iiBen, ppK14n,ppK4n, rate )

  !----------------------------------------------------------------------
  ! Silicate remineralization in the sediments
  !----------------------------------------------------------------------
  rate  =   p_reminQ6s* Q6s(:)
  call flux_vector( iiBen, ppQ6s,ppQ6s,-( rate) )
  ! K5.s is not used in this model version

  jK5N5s(:)  =   rate





  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
