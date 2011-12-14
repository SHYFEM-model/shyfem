#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenthicReturn1
!
! DESCRIPTION
!   This process is a very simple parameterisation of benthic 
!   remineralisation.
!   Benthic organism sub-model cannot be used with this model.
!   A constant portion of the organic matter in the sediments is
!   released to the water column as inorganic nutrients.
!   Oxygen consumption is stoichiometrically associated to carbon 
!   remineralisation rates and nitrogen remineralisation is partitioned 
!   into ammonium and nitrate flux with a constant value. 
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenthicReturn1Dynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: Q6c, Q1c, Q6p, Q1p, &
  ! K1p, Q6n, Q1n, K3n, K4n, Q6s, K5s
  ! The following Benthic 1-d global boxvars are modified : jG2O2o, jK1N1p, &
  ! jK3N3n, jK4N4n, jK5N5s
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem,  ONLY: Q6c, Q1c, Q6p, Q1p, K1p, Q6n, Q1n, K3n, K4n, Q6s, K5s
#endif
  use mem, ONLY: ppQ6c, ppQ1c, ppQ6p, ppQ1p, ppK1p, ppQ6n, ppQ1n, &
    ppK3n, ppK4n, ppQ6s, ppK5s, jG2O2o, jK1N1p, jK3N3n, jK4N4n, jK5N5s, &
    NO_BOXES_XY, iiBen, iiPel, flux_vector
  use mem_BenthicReturn1



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

  rate  =   p_reminQ6c* Q6c(:)
  call flux_vector( iiBen, ppQ6c,ppQ6c,-( rate) )
  jG2O2o(:)  =  - rate/ 12.D+00

  rate  =   p_reminQ1c* Q1c(:)
  call flux_vector( iiBen, ppQ1c,ppQ1c,-( rate) )
  jG2O2o(:)  =   jG2O2o(:)- rate/ 12.D+00

  rate  =   p_reminQ6p* Q6p(:)
  call flux_vector( iiBen, ppQ6p,ppQ6p,-( rate) )
  jK1N1p(:)  =   rate

  rate  =   p_reminQ1p* Q1p(:)
  call flux_vector( iiBen, ppQ1p,ppQ1p,-( rate) )
  jK1N1p(:)  =   jK1N1p(:)+ rate

  ! K1.p is not used in this model version
  ! hence the amount is added which is subtracted in BentoPelCoupling
  call flux_vector( iiBen, ppK1p,ppK1p,-(- jK1N1p(:)) )

  rate  =   p_reminQ6n* Q6n(:)
  call flux_vector( iiBen, ppQ6n,ppQ6n,-( rate) )
  jK3N3n(:)  =   rate* p_pQIN3
  jK4N4n(:)  =   rate*( 1.0D+00- p_pQIN3)

  rate  =   p_reminQ1n* Q1n(:)
  call flux_vector( iiBen, ppQ1n,ppQ1n,-( rate) )
  jK3N3n(:)  =   jK3N3n(:)+ rate* p_pQIN3
  jK4N4n(:)  =   jK4N4n(:)+ rate*( 1.0D+00- p_pQIN3)

  ! K3.n and K4.n are not used in this model version
  ! hence the amount is added which is subtracted in BentoPelCoupling
  call flux_vector( iiBen, ppK3n,ppK3n,-(- jK3N3n(:)) )
  call flux_vector( iiBen, ppK4n,ppK4n,-(- jK4N4n(:)) )


  rate  =   p_reminQ6s* Q6s(:)
  call flux_vector( iiBen, ppQ6s,ppQ6s,-( rate) )
  jK5N5s(:)  =   rate

  ! K5.s is not used in this model version
  ! hence the amount is added which is subtracted in BentoPelCoupling
  call flux_vector( iiBen, ppK5s,ppK5s,-(- jK5N5s(:)) )


  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
