#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenOrganism
!
! DESCRIPTION
!   !    This submodel describes the carbon dynamics and associated
!    nutrient dynamics in benthic organisms
! 
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenOrganismDynamics(y,  ppyc, ppyn, ppyp)
!
! !USES:

  ! For the following Benthic-states fluxes are defined: Q6c, Q6n, Q6p, G2o, &
  ! K4n, K1p, D6m, D7m, D8m
  ! The following Benthic-states are used (NOT in fluxes): D1m
  ! For the following Benthic-group-states fluxes are defined: BenOrganisms, &
  ! BenBacteria
  ! The following Benthic 1-d global boxvars are modified : rrBTo, reBTn, reBTp
  ! The following Benthic 1-d global boxvars  are used: ETW_Ben
  ! The following groupmember vars are used: iiBenOrganisms, iiBenBacteria, &
  ! iiY1, iiY2, iiY4, iiY5
  ! The following constituent constants  are used: iiC, iiN, iiP
  ! The following 0-d global box parametes are used: p_d_tot
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem, ONLY: D2STATE, Q6c, Q6n, Q6p, G2o, K4n, K1p, D6m, D7m, &
    D8m, D1m, BenOrganisms, BenBacteria
#endif
  use mem, ONLY: ppQ6c, ppQ6n, ppQ6p, ppG2o, ppK4n, ppK1p, &
    ppD6m, ppD7m, ppD8m, ppD1m, ppBenOrganisms, ppBenBacteria, rrBTo, reBTn, &
    reBTp, ETW_Ben, iiBenOrganisms, iiBenBacteria, iiY1, iiY2, iiY4, iiY5, iiC, &
    iiN, iiP, NO_BOXES_XY, iiBen, iiPel, flux_vector
  use mem_Param,  ONLY: p_d_tot
  use mem_BenOrganism


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, eramp_vector, &
  ! MM_vector, PartQ_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, eramp_vector, MM_vector, PartQ_vector


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: y
  integer,intent(IN) :: ppyc
  integer,intent(IN) :: ppyn
  integer,intent(IN) :: ppyp

!  
!
! !AUTHORS
!   W. Ebenhoh and C. Kohlmeier.
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
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY) :: yc
  real(RLEN),dimension(NO_BOXES_XY) :: yn
  real(RLEN),dimension(NO_BOXES_XY) :: yp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  real(RLEN),dimension(NO_BOXES_XY)  :: clm
  real(RLEN),dimension(NO_BOXES_XY)  :: cm
  real(RLEN),dimension(NO_BOXES_XY)  :: cmm
  real(RLEN),dimension(NO_BOXES_XY)  :: et
  real(RLEN),dimension(NO_BOXES_XY)  :: eo
  real(RLEN),dimension(NO_BOXES_XY)  :: food
  real(RLEN),dimension(NO_BOXES_XY)  :: food_src
  real(RLEN),dimension(NO_BOXES_XY)  :: eF
  real(RLEN),dimension(NO_BOXES_XY)  :: rgu
  real(RLEN),dimension(NO_BOXES_XY)  :: sgu
  real(RLEN),dimension(NO_BOXES_XY)  :: sgu_y
  real(RLEN),dimension(NO_BOXES_XY)  :: sguQ6
  real(RLEN),dimension(NO_BOXES_XY)  :: snu
  real(RLEN),dimension(NO_BOXES_XY)  :: snuQ6
  real(RLEN),dimension(NO_BOXES_XY)  :: se_u
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uQ6
  real(RLEN),dimension(NO_BOXES_XY)  :: choice
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_c
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_n
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_p
  real(RLEN),dimension(NO_BOXES_XY)  :: rtyc
  real(RLEN),dimension(NO_BOXES_XY)  :: rtyn
  real(RLEN),dimension(NO_BOXES_XY)  :: rtyp
  real(RLEN),dimension(NO_BOXES_XY)  :: rqt6c
  real(RLEN),dimension(NO_BOXES_XY)  :: rqt6n
  real(RLEN),dimension(NO_BOXES_XY)  :: rqt6p
  real(RLEN),dimension(NO_BOXES_XY)  :: rq6c
  real(RLEN),dimension(NO_BOXES_XY)  :: rq6n
  real(RLEN),dimension(NO_BOXES_XY)  :: rq6p
  real(RLEN),dimension(NO_BOXES_XY)  :: ruYIc
  real(RLEN),dimension(NO_BOXES_XY)  :: ruYIn
  real(RLEN),dimension(NO_BOXES_XY)  :: ruYIp
  real(RLEN),dimension(NO_BOXES_XY)  :: ruBIc
  real(RLEN),dimension(NO_BOXES_XY)  :: ruBIn
  real(RLEN),dimension(NO_BOXES_XY)  :: ruBIp
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6c
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6n
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6p
  real(RLEN),dimension(NO_BOXES_XY)  :: rrc
  real(RLEN),dimension(NO_BOXES_XY)  :: sm
  real(RLEN),dimension(NO_BOXES_XY)  :: ren
  real(RLEN),dimension(NO_BOXES_XY)  :: rep
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  yc = D2STATE(ppyc,:)
  yn = D2STATE(ppyn,:)
  yp = D2STATE(ppyp,:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Physiological temperature and oxygen response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  et  =   eTq_vector(  ETW_Ben(:),  p_q10(y))

  eo  =   eramp_vector(  D1m(:),  p_sdm(y))

  ! As alternative the following function can be used
  ! eo = MM(pow(O2.o, 3.0), p_chdo);
  ! and the parameter p_chdo must be defined


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  food  =   0.0D+00

  rtyc  =   0.0D+00
  rtyn  =   0.0D+00
  rtyp  =   0.0D+00

  rqt6c  =   0.0D+00
  rqt6n  =   0.0D+00
  rqt6p  =   0.0D+00

  ! For other benthic organisms:

  do i = 1 , ( iiBenOrganisms)


    food_src  =   BenOrganisms(i,iiC)* p_Yn(y,i)
    food  =   food+ food_src* MM_vector(  food_src,  p_clu(y))

  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Bacteria:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


  do i = 1 , ( iiBenBacteria)


    food_src  =   BenBacteria(i,iiC)* p_Hn(y,i)
    food  =   food+ food_src* MM_vector(  food_src,  p_clu(y))

  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Detritus (if eaten) First calculate the available portion
  ! and then add it to the total amount of food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if ( p_puQ6(y)> 0.0D+00) then

    clm  =   p_clm(y)
    cm  =   p_cm(y)
    availQ6_c  =   Q6c(:)* PartQ_vector(  D6m(:),  clm,  cm,  p_d_tot)
    availQ6_n  =   Q6n(:)* PartQ_vector(  D7m(:),  clm,  cm,  p_d_tot)
    availQ6_p  =   Q6p(:)* PartQ_vector(  D8m(:),  clm,  cm,  p_d_tot)

    food_src  =   p_puQ6(y)* availQ6_c
    food  =   food+ food_src* MM_vector(  food_src,  p_clu(y))

  end if



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Correct for too much food...
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  eF  =   MM_vector(  food,  p_chu(y))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Correction of growth rate for environmental factors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! Growth rate at actual amount:

  sgu_y  =   p_su(y)* et* eo* eF

  ! Relative growth rate corrected for actual amount of food:

  sgu  =  ( sgu_y* yc)/ food

  ! Net uptake:

  snu  =   sgu*( 1.0D+00- p_pue(y))
  snuQ6  =   sgu*( 1.0D+00- p_pueQ6(y))

  ! Execreted part:

  se_u  =   sgu- snu
  se_uQ6  =   sgu- snuQ6


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of uptake rate:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Benthic organisms:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  do i = 1 , ( iiBenOrganisms)


    choice = p_Yn(y,i)* MM_vector( p_Yn(y, i)* BenOrganisms(i, iiC), &
      p_clu(y))

    ruYIc  =   BenOrganisms(i,iiC)* sgu* choice
    ruYIn  =   BenOrganisms(i,iiN)* sgu* choice
    ruYIp  =   BenOrganisms(i,iiP)* sgu* choice

    ! In case of cannibalism rate of change in state = zero!

    if ( i/= y) then
      call flux_vector( iiBen, ppBenOrganisms(i,iiC),ppyc, ruYIc )
      call flux_vector( iiBen, ppBenOrganisms(i,iiN),ppyn, ruYIn )
      call flux_vector( iiBen, ppBenOrganisms(i,iiP),ppyp, ruYIp )
    end if


    rtyc  =   rtyc+ ruYIc
    rtyn  =   rtyn+ ruYIn
    rtyp  =   rtyp+ ruYIp

    rq6c  =   BenOrganisms(i,iiC)* se_u* choice
    rq6n  =   BenOrganisms(i,iiN)* se_u* p_pudil(y)* choice
    rq6p  =   BenOrganisms(i,iiP)* se_u* p_pudil(y)* choice

    rqt6c  =   rqt6c+ rq6c
    rqt6n  =   rqt6n+ rq6n
    rqt6p  =   rqt6p+ rq6p

  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Bacteria:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  do i = 1 , ( iiBenBacteria)


    choice = p_Hn(y,i)* MM_vector( p_Hn(y, i)* BenBacteria(i, iiC), &
      p_clu(y))

    ruBIc  =   BenBacteria(i,iiC)* sgu* choice
    ruBIn  =   BenBacteria(i,iiN)* sgu* choice
    ruBIp  =   BenBacteria(i,iiP)* sgu* choice

    call flux_vector( iiBen, ppBenBacteria(i,iiC),ppyc, ruBIc )
    call flux_vector( iiBen, ppBenBacteria(i,iiN),ppyn, ruBIn )
    call flux_vector( iiBen, ppBenBacteria(i,iiP),ppyp, ruBIp )

    rtyc  =   rtyc+ ruBIc
    rtyn  =   rtyn+ ruBIn
    rtyp  =   rtyp+ ruBIp

    rq6c  =   BenBacteria(i,iiC)* se_u* choice
    rq6n  =   BenBacteria(i,iiN)* se_u* p_pudil(y)* choice
    rq6p  =   BenBacteria(i,iiP)* se_u* p_pudil(y)* choice

    rqt6c  =   rqt6c+ rq6c
    rqt6n  =   rqt6n+ rq6n
    rqt6p  =   rqt6p+ rq6p

  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! detritus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( p_puQ6(y)> 0.0D+00)

    case( .TRUE. )
      choice  =   p_puQ6(y)* MM_vector(  p_puQ6(y)* availQ6_c,  p_clu(y))

      ruQ6c  =   sgu* choice* availQ6_c
      ruQ6n  =   sgu* choice* availQ6_n
      ruQ6p  =   sgu* choice* availQ6_p

      call flux_vector( iiBen, ppQ6c,ppyc, ruQ6c )
      call flux_vector( iiBen, ppQ6n,ppyn, ruQ6n )
      call flux_vector( iiBen, ppQ6p,ppyp, ruQ6p )

      rtyc  =   rtyc+ ruQ6c
      rtyn  =   rtyn+ ruQ6n
      rtyp  =   rtyp+ ruQ6p

      rq6c  =   se_uQ6* choice* availQ6_c
      rq6n  =   se_uQ6* p_pudil(y)* choice* availQ6_n
      rq6p  =   se_uQ6* p_pudil(y)* choice* availQ6_p

      rqt6c  =   rqt6c+ rq6c
      rqt6n  =   rqt6n+ rq6n
      rqt6p  =   rqt6p+ rq6p



    case( .FALSE. )
      ruQ6c  =   0.0D+00
      ruQ6n  =   0.0D+00
      ruQ6p  =   0.0D+00


  end select




  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of respiration:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrc  =   p_sr(y)* yc* et+ p_pur(y)*( sgu* food- rqt6c)

  call flux_vector( iiBen, ppyc,ppyc,-( rrc) )
  call flux_vector(iiBen, ppG2o,ppG2o,-( rrc/ 12.0D+00))
  rtyc  =   rtyc- rrc


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  sm  =   p_sd(y)* et

  rq6c  =   yc* sm
  rq6n  =   yn* sm
  rq6p  =   yp* sm

  rqt6c  =   rqt6c+ rq6c
  rqt6n  =   rqt6n+ rq6n
  rqt6p  =   rqt6p+ rq6p


  ! in case of a negative value of one of the following values there is a &
  ! situation
  ! of startvation and very low biomass values. Check on quota in the food is &
  ! out of order
  rtyc  =   max(  0.0D+00,  rtyc- rqt6c)
  rtyn  =   max(  0.0D+00,  rtyn- rqt6n)
  rtyp  =   max(  0.0D+00,  rtyp- rqt6p)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of nutrient release and correction of C:N:P
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ren  =   rtyn- rtyc* p_qn(y)
  rep  =   rtyp- rtyc* p_qp(y)

  where ( ren< 0.0D+00)

    rq6c  =  - ren/ p_qn(y)
    rqt6c  =   rqt6c+ rq6c
    rtyc  =   rtyc- rq6c

    ren  =   rtyn- rtyc* p_qn(y)
    rep  =   rtyp- rtyc* p_qp(y)

  end where

  where ( rep< 0.0D+00)

    rq6c  =  - rep/ p_qp(y)
    rqt6c  =   rqt6c+ rq6c
    rtyc  =   rtyc- rq6c

    ren  =   rtyn- rtyc* p_qn(y)
    rep  =   rtyp- rtyc* p_qp(y)

  end where


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Correction for cases where initial conditions deviate strongly from
  ! Redfield C:N:P. In this way the C:N:P does not become too extreme
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ren = min( max( 0.0D+00, ren), max( 0.0D+00, ren-( p_qn(y)* yc- &
    yn)))
  rep = min( max( 0.0D+00, rep), max( 0.0D+00, rep-( p_qp(y)* yc- &
    yp)))


  call flux_vector( iiBen, ppyn,ppK4n, ren )
  call flux_vector( iiBen, ppyp,ppK1p, rep )


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total flux of Y to Q6:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector( iiBen, ppyc,ppQ6c, rqt6c )
  call flux_vector( iiBen, ppyn,ppQ6n, rqt6n )
  call flux_vector( iiBen, ppyp,ppQ6p, rqt6p )


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Add respiration and excretion to the benthic totals
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrBTo(:)  =   rrBTo(:)+ rrc/ 12.0D+00
  reBTn(:)  =   reBTn(:)+ ren
  reBTp(:)  =   reBTp(:)+ rep


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Assign organism-dependent parameters
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( y)

    case ( iiY1 )
      cmm  =   0.0D+00



    case ( iiY2 )
      cmm  =   D6m(:)/ 2.0D+00



    case ( iiY4 )
      cmm  =  ( p_cm(y)+ p_clm(y))/ 2.0D+00



    case ( iiY5 )
      cmm  =   D1m(:)/ 2.0D+00



  end select



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of changes due to uptake of detritus in distribution of
  ! state variables (Dx.m is an undetermined source).
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector(iiBen, ppD6m,ppD6m,( cmm- D6m(:))*( rqt6c- ruQ6c)/ Q6c(:))
  call flux_vector(iiBen, ppD7m,ppD7m,( cmm- D7m(:))*( rqt6n- ruQ6n)/ Q6n(:))
  call flux_vector(iiBen, ppD8m,ppD8m,( cmm- D8m(:))*( rqt6p- ruQ6p)/ Q6p(:))




  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
