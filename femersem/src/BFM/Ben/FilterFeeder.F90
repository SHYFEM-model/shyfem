#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: FilterFeeder
!
! DESCRIPTION
!   This process describes the carbon dynamics and associated
!   nutrient dynamics in benthic organism Y3 (suspension feeders)
!   Y3 is handled separately because it also feeds from the water column.
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine FilterFeederDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: Y3c, Y3n, Y3p, Q6c, &
  ! Q6n, Q6p, G2o, K4n, K1p, D6m, D7m, D8m
  ! The following Benthic-states are used (NOT in fluxes): D1m
  ! The following Benthic 1-d global boxvars are modified : rrBTo, reBTn, &
  ! reBTp, rutQ6c, rutQ6n, rutQ6p, rutQ6s
  ! The following Benthic 1-d global boxvars got a value: jPIY3c, jRIY3c, &
  ! jRIY3n, jRIY3p, jRIY3s
  ! The following Benthic 1-d global boxvars are used: ETW_Ben, PIc, RIc, PIn, &
  ! PIp, PIs, RIn, RIp, RIs
  ! The following 0-d global box parametes are used: p_d_tot
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem, ONLY: Y3c, Y3n, Y3p, Q6c, Q6n, Q6p, G2o, K4n, K1p, D6m, D7m, D8m, &
    D1m
#endif
  use mem, ONLY: ppY3c, ppY3n, ppY3p, ppQ6c, ppQ6n, ppQ6p, ppG2o, ppK4n, &
    ppK1p, ppD6m, ppD7m, ppD8m, ppD1m, rrBTo, reBTn, reBTp, rutQ6c, rutQ6n, &
    rutQ6p, rutQ6s, jPIY3c, jRIY3c, jRIY3n, jRIY3p, jRIY3s, ETW_Ben, PIc, RIc, &
    PIn, PIp, PIs, RIn, RIp, RIs, NO_BOXES_XY, iiBen, iiPel, flux_vector
  use mem_Param,  ONLY: p_d_tot
  use mem_FilterFeeder


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, eramp_vector, &
  ! MM_vector, PartQ_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, eramp_vector, MM_vector, PartQ_vector



!  
!
! !AUTHORS
!   W. Ebenhoeh and C. Kohlmeier 
!
!
!
! !REVISION_HISTORY
!   !
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
  integer,dimension(NO_BOXES_XY)  :: i
  real(RLEN),dimension(NO_BOXES_XY)  :: clm
  real(RLEN),dimension(NO_BOXES_XY)  :: cmm
  real(RLEN),dimension(NO_BOXES_XY)  :: cm
  real(RLEN),dimension(NO_BOXES_XY)  :: et
  real(RLEN),dimension(NO_BOXES_XY)  :: eo
  real(RLEN),dimension(NO_BOXES_XY)  :: food
  real(RLEN),dimension(NO_BOXES_XY)  :: food_src
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_c
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_n
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_p
  real(RLEN),dimension(NO_BOXES_XY)  :: eF
  real(RLEN),dimension(NO_BOXES_XY)  :: sgu
  real(RLEN),dimension(NO_BOXES_XY)  :: rgu
  real(RLEN),dimension(NO_BOXES_XY)  :: snu
  real(RLEN),dimension(NO_BOXES_XY)  :: snuQ6
  real(RLEN),dimension(NO_BOXES_XY)  :: se_u
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uQ6
  real(RLEN),dimension(NO_BOXES_XY)  :: choice
  real(RLEN),dimension(NO_BOXES_XY)  :: rtY3c
  real(RLEN),dimension(NO_BOXES_XY)  :: rtY3n
  real(RLEN),dimension(NO_BOXES_XY)  :: rtY3p
  real(RLEN),dimension(NO_BOXES_XY)  :: rrc
  real(RLEN),dimension(NO_BOXES_XY)  :: sm
  real(RLEN),dimension(NO_BOXES_XY)  :: ren
  real(RLEN),dimension(NO_BOXES_XY)  :: rep
  real(RLEN),dimension(NO_BOXES_XY)  :: retR6c
  real(RLEN),dimension(NO_BOXES_XY)  :: retR6n
  real(RLEN),dimension(NO_BOXES_XY)  :: retR6p
  real(RLEN),dimension(NO_BOXES_XY)  :: retQ6c
  real(RLEN),dimension(NO_BOXES_XY)  :: retQ6n
  real(RLEN),dimension(NO_BOXES_XY)  :: retQ6p
  real(RLEN),dimension(NO_BOXES_XY)  :: rePIc
  real(RLEN),dimension(NO_BOXES_XY)  :: rePIn
  real(RLEN),dimension(NO_BOXES_XY)  :: rePIp
  real(RLEN),dimension(NO_BOXES_XY)  :: reR6c
  real(RLEN),dimension(NO_BOXES_XY)  :: reR6n
  real(RLEN),dimension(NO_BOXES_XY)  :: reR6p
  real(RLEN),dimension(NO_BOXES_XY)  :: reQ6c
  real(RLEN),dimension(NO_BOXES_XY)  :: reQ6n
  real(RLEN),dimension(NO_BOXES_XY)  :: reQ6p
  real(RLEN),dimension(NO_BOXES_XY)  :: ruPIc
  real(RLEN),dimension(NO_BOXES_XY)  :: ruPIn
  real(RLEN),dimension(NO_BOXES_XY)  :: ruPIp
  real(RLEN),dimension(NO_BOXES_XY)  :: ruPIs
  real(RLEN),dimension(NO_BOXES_XY)  :: ruR6c
  real(RLEN),dimension(NO_BOXES_XY)  :: ruR6n
  real(RLEN),dimension(NO_BOXES_XY)  :: ruR6p
  real(RLEN),dimension(NO_BOXES_XY)  :: ruR6s
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6c
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6n
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6p
  real(RLEN),dimension(NO_BOXES_XY)  :: ruK1p
  real(RLEN),dimension(NO_BOXES_XY)  :: ruK4n
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Physiological temperature response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  et  =   eTq_vector(  ETW_Ben(:),  p_q10)

  eo  =   eramp_vector(  D1m(:),  0.005D+00)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  food  =   0.0D+00

  ! For other benthic organisms:

  food_src  =   PIc(:)* p_dwat
  food  =   food+ p_PI* food_src* MM_vector(  food_src,  p_clu)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Detritus, (if eaten) first calculate available amount
  ! and add it to the total amount of food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  food_src  =   RIc(:)* p_dwat
  food  =   food+ p_R6* food_src* MM_vector(  food_src,  p_clu)

  clm  =   p_clm
  cm  =   p_cm
  availQ6_c  =   Q6c(:)* PartQ_vector(  D6m(:),  clm,  cm,  p_d_tot)
  availQ6_n  =   Q6n(:)* PartQ_vector(  D7m(:),  clm,  cm,  p_d_tot)
  availQ6_p  =   Q6p(:)* PartQ_vector(  D8m(:),  clm,  cm,  p_d_tot)

  food_src  =   availQ6_c
  food  =   food+ p_puQ6* food_src* MM_vector(  food_src,  p_clu)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Correct for too much food:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  eF  =   MM_vector(  food,  p_chu)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Correction of growth rate for environmental factors:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! Growth rate at actual amount:

  rgu  =   p_su* Y3c(:)* et* eo* eF

  ! Relative growth rate corrected for actual amount of food:

  sgu  =   rgu/ food

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Net uptake:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  snu  =   sgu*( 1.0D+00- p_pue)
  snuQ6  =   sgu*( 1.0D+00- p_pueQ6)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Execreted part:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  se_u  =   sgu- snu
  se_uQ6  =   sgu- snuQ6

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of uptake rate:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pelagic Phytoplankton:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  choice  =   p_PI* MM_vector(  PIc(:),  p_clu)* p_dwat

  ruPIc  =   PIc(:)* sgu* choice
  ruPIn  =   PIn(:)* sgu* choice
  ruPIp  =   PIp(:)* sgu* choice
  ruPIs  =   PIs(:)* sgu* choice

  call flux_vector( iiBen, ppY3c,ppY3c, ruPIc )
  call flux_vector( iiBen, ppY3n,ppY3n, ruPIn )
  call flux_vector( iiBen, ppY3p,ppY3p, ruPIp )


  ! flux definitions from P -> Y3 are found in BentoPelCoup
  jPIY3c(:)  =   ruPIc

  rePIc  =   PIc(:)* se_u* choice
  rePIn  =   PIn(:)* se_u* p_pudsil* choice
  rePIp  =   PIp(:)* se_u* p_pudsil* choice

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pelagic Detritus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  choice  =   p_R6* MM_vector(  RIc(:),  p_clu)* p_dwat

  ruR6c  =   RIc(:)* sgu* choice
  ruR6n  =   RIn(:)* sgu* choice
  ruR6p  =   RIp(:)* sgu* choice
  ruR6s  =   RIs(:)* sgu* choice

  call flux_vector( iiBen, ppY3c,ppY3c, ruR6c )
  call flux_vector( iiBen, ppY3n,ppY3n, ruR6n )
  call flux_vector( iiBen, ppY3p,ppY3p, ruR6p )

  reR6c  =   RIc(:)* se_uQ6* choice
  reR6n  =   RIn(:)* se_uQ6* p_pudsil* choice
  reR6p  =   RIp(:)* se_uQ6* p_pudsil* choice

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Benthic Detritus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  choice  =   p_puQ6* MM_vector(  availQ6_c,  p_clu)

  ruQ6c  =   sgu* choice* availQ6_c
  ruQ6n  =   sgu* choice* availQ6_n
  ruQ6p  =   sgu* choice* availQ6_p

  call flux_vector( iiBen, ppQ6c,ppY3c, ruQ6c )
  call flux_vector( iiBen, ppQ6n,ppY3n, ruQ6n )
  call flux_vector( iiBen, ppQ6p,ppY3p, ruQ6p )

  reQ6c  =   se_uQ6* choice* availQ6_c
  reQ6n  =   se_uQ6* p_pudsil* choice* availQ6_n
  reQ6p  =   se_uQ6* p_pudsil* choice* availQ6_p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Book keeping
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rtY3c  =   ruPIc+ ruR6c+ ruQ6c
  rtY3n  =   ruPIn+ ruR6n+ ruQ6n
  rtY3p  =   ruPIp+ ruR6p+ ruQ6p

  retR6c  =   rePIc+ reR6c
  retR6n  =   rePIn+ reR6n
  retR6p  =   rePIp+ reR6p

  retQ6c  =   reQ6c
  retQ6n  =   reQ6n
  retQ6p  =   reQ6p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of respiration:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrc  =   p_sr* Y3c(:)* et+ p_pur*( food* sgu- retR6c- retQ6c)

  call flux_vector( iiBen, ppY3c,ppY3c,-( rrc) )
  call flux_vector(iiBen, ppG2o,ppG2o,-( rrc/ 12.0D+00))

  rtY3c  =   rtY3c- rrc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- =-=-=-=-=-=-=-=-=-=-=-=-=-

  sm  =   p_sd* et

  reQ6c  =   Y3c(:)* sm
  reQ6n  =   Y3n(:)* sm
  reQ6p  =   Y3p(:)* sm

  retQ6c  =   retQ6c+ reQ6c
  retQ6n  =   retQ6n+ reQ6n
  retQ6p  =   retQ6p+ reQ6p

  ! in case of a negative value of one of the following values there is a &
  ! situation
  ! of startvation and very low biomass values. Check on quota in the food is &
  ! out of order
  rtY3c  =   max(  0.0D+00,  rtY3c- retQ6c)
  rtY3n  =   max(  0.0D+00,  rtY3n- retQ6n)
  rtY3p  =   max(  0.0D+00,  rtY3p- retQ6p)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of nutrient release and correction of C:N:P
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ren  =   rtY3n- rtY3c* p_qn
  rep  =   rtY3p- rtY3c* p_qp

  where ( ren< 0.0D+00)

    reQ6c  =  - ren/ p_qn
    retQ6c  =   retQ6c+ reQ6c
    rtY3c  =   rtY3c- reQ6c

    ren  =   rtY3n- rtY3c* p_qn
    rep  =   rtY3p- rtY3c* p_qp

  end where


  where ( rep< 0.0D+00)

    reQ6c  =  - rep/ p_qp
    retQ6c  =   retQ6c+ reQ6c
    rtY3c  =   rtY3c- reQ6c

    ren  =   rtY3n- rtY3c* p_qn
    rep  =   rtY3p- rtY3c* p_qp

  end where


  ren = min( max( 0.0D+00, ren), max( 0.0D+00, ren-( p_qn* Y3c(:)- &
    Y3n(:))))
  rep = min( max( 0.0D+00, rep), max( 0.0D+00, rep-( p_qp* Y3c(:)- &
    Y3p(:))))

  call flux_vector( iiBen, ppY3n,ppK4n, ren )
  call flux_vector( iiBen, ppY3p,ppK1p, rep )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Add respiration and excretion to the benthic totals
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrBTo(:)  =   rrBTo(:)+ rrc/ 12.0D+00
  reBTn(:)  =   reBTn(:)+ ren
  reBTp(:)  =   reBTp(:)+ rep



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total flux from Suspension feeders to Q6:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector( iiBen, ppY3c,ppQ6c, retQ6c )
  call flux_vector( iiBen, ppY3n,ppQ6n, retQ6n )
  call flux_vector( iiBen, ppY3p,ppQ6p, retQ6p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of changes due to uptake of Benthic detritus in distribution of
  ! state variables (Dx.m is an undetermined source)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  cmm  =  ( clm+ cm)* 0.5D+00
  call flux_vector(iiBen, ppD6m,ppD6m,( cmm- D6m(:))*( retQ6c- ruQ6c)/ Q6c(:))
  call flux_vector(iiBen, ppD7m,ppD7m,( cmm- D7m(:))*( retQ6n- ruQ6n)/ Q6n(:))
  call flux_vector(iiBen, ppD8m,ppD8m,( cmm- D8m(:))*( retQ6p- ruQ6p)/ Q6p(:))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total flux from Suspension feeders to R6:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector( iiBen, ppY3c,ppY3c,-( retR6c) )
  call flux_vector( iiBen, ppY3n,ppY3n,-( retR6n) )
  call flux_vector( iiBen, ppY3p,ppY3p,-( retR6p) )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate NET flux from R6 to Suspension feeders :
  ! (can be negative!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  jRIY3c(:)  =   ruR6c- retR6c
  jRIY3n(:)  =   ruR6n- retR6n
  jRIY3p(:)  =   ruR6p- retR6p
  ! The ruR6s which is uptaken is directly relased back to R6: net food flux &
  ! from Y3 to/from R6 is 0
  jRIY3s(:)  =  - ruPIs

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! excretion of food orginating of Pelagic food will be sedimented and
  ! hence is considered as added to the total Pelagic->Ben flux.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rutQ6c(:)  =   rutQ6c(:)+ retR6c
  rutQ6n(:)  =   rutQ6n(:)+ retR6n
  rutQ6p(:)  =   rutQ6p(:)+ retR6p

  ! The silicate is directly transferred to Q6.s
  ! the ruPis which is put back in R6 is however sedimentating:
  rutQ6s(:)  =   rutQ6s(:)+ ruPIs+ ruR6s




  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
