#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenBac
!
! DESCRIPTION
!   !    This submodel describes the carbon dynamics and associated
!    nutrient dynamics in benthic bacteria (represented
!    by state variables H1-H2) 
! 
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenBacDynamics(hx,  pphxc, pphxn, pphxp)
!
! !USES:

  ! For the following Benthic-states fluxes are defined: Q6c, Q6n, Q6p, G2o, &
  ! K6r, D6m, D7m, D8m
  ! The following Benthic-states are used (NOT in fluxes): D1m, D2m
  ! For the following Benthic-group-states fluxes are defined: &
  ! BenDetritus, BenthicAmmonium, BenthicPhosphate
  ! The following Benthic 1-d global boxvars are modified : rrBTo, reBTn, &
  ! reBTp, rrATo, reATn, reATp
  ! The following Benthic 1-d global boxvars  are used: ETW_Ben
  ! The following Benthic 2-d global boxvars got a value: ruHI
  ! The following groupmember vars  are used: iiH1, iiH2
  ! The following constituent constants  are used: iiC, iiN, iiP
  ! The following 0-d global box parametes are used: p_d_tot, p_small, &
  ! p_pe_R1c, p_pe_R1n, p_pe_R1p, p_qro
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem, ONLY: D2STATE, Q6c, Q6n, Q6p, G2o, K6r, D6m, &
    D7m, D8m, D1m, D2m, BenDetritus, BenthicAmmonium, BenthicPhosphate
#endif
  use mem, ONLY: ppQ6c, ppQ6n, ppQ6p, ppG2o, ppK6r, ppD6m, &
    ppD7m, ppD8m, ppD1m, ppD2m, ppBenDetritus, ppBenthicAmmonium, &
    ppBenthicPhosphate, rrBTo, reBTn, reBTp, rrATo, reATn, reATp, ETW_Ben, ruHI, &
    iiH1, iiH2, iiC, iiN, iiP, NO_BOXES_XY, iiBen, iiPel, flux_vector
  use mem_Param,  ONLY: p_d_tot, p_small, p_pe_R1c, p_pe_R1n, p_pe_R1p, p_qro
  use mem_BenBac


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector, &
  ! PartQ_vector, eramp_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun, ONLY: eTq_vector, MM_vector, PartQ_vector, eramp_vector, &
    insw_vector


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: hx
  integer,intent(IN) :: pphxc
  integer,intent(IN) :: pphxn
  integer,intent(IN) :: pphxp

!  
!
! !AUTHORS
!   W. Ebenhoh and C. Kohlmeier 
! 
!
!
! !REVISION_HISTORY
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
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY) :: hxc
  real(RLEN),dimension(NO_BOXES_XY) :: hxn
  real(RLEN),dimension(NO_BOXES_XY) :: hxp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY)  :: clm
  real(RLEN),dimension(NO_BOXES_XY)  :: cm
  real(RLEN),dimension(NO_BOXES_XY)  :: chm
  real(RLEN),dimension(NO_BOXES_XY)  :: cmm
  real(RLEN),dimension(NO_BOXES_XY)  :: et
  real(RLEN),dimension(NO_BOXES_XY)  :: eN
  real(RLEN),dimension(NO_BOXES_XY)  :: eo
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_c
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_p
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_n
  real(RLEN),dimension(NO_BOXES_XY)  :: suQ1
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ1c
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ1n
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ1p
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6c
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6n
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6p
  real(RLEN),dimension(NO_BOXES_XY)  :: ruKIn
  real(RLEN),dimension(NO_BOXES_XY)  :: ruKIp
  real(RLEN),dimension(NO_BOXES_XY)  :: rrc
  real(RLEN),dimension(NO_BOXES_XY)  :: sm
  real(RLEN),dimension(NO_BOXES_XY)  :: ren
  real(RLEN),dimension(NO_BOXES_XY)  :: rep
  real(RLEN),dimension(NO_BOXES_XY)  :: sHc
  real(RLEN),dimension(NO_BOXES_XY)  :: sHn
  real(RLEN),dimension(NO_BOXES_XY)  :: sHp
  real(RLEN),dimension(NO_BOXES_XY)  :: rq6c
  real(RLEN),dimension(NO_BOXES_XY)  :: rqt6c
  real(RLEN),dimension(NO_BOXES_XY)  :: rqt6n
  real(RLEN),dimension(NO_BOXES_XY)  :: rqt6p
  real(RLEN),dimension(NO_BOXES_XY)  :: qnQ6c
  real(RLEN),dimension(NO_BOXES_XY)  :: qpQ6c
  real(RLEN),dimension(NO_BOXES_XY)  :: rut
  real(RLEN),dimension(NO_BOXES_XY)  :: rum
  real(RLEN),dimension(NO_BOXES_XY)  :: run
  real(RLEN),dimension(NO_BOXES_XY)  :: rug
  real(RLEN),dimension(NO_BOXES_XY)  :: netgrowth
  real(RLEN),dimension(NO_BOXES_XY)  :: r
  real(RLEN),dimension(NO_BOXES_XY)  :: runp
  real(RLEN),dimension(NO_BOXES_XY)  :: runn
  integer,dimension(NO_BOXES_XY)  :: i
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  hxc = D2STATE(pphxc,:)
  hxn = D2STATE(pphxn,:)
  hxp = D2STATE(pphxp,:)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Assign functional group-dependent parameters:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( hx)

    case ( iiH1 )
      clm  =   0.0D+00
      cm  =   D1m(:)
      chm  =   D1m(:)
      cmm  =   D1m(:)/ 2.0D+00



    case ( iiH2 )
      clm  =   D1m(:)
      cm  =   D2m(:)
      chm  =   p_d_tot
      cmm  =   D6m(:)+ D1m(:)



  end select



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Physiological temperature response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  et  =   eTq_vector(  ETW_Ben(:),  p_q10(hx))

  eo  =   MM_vector(  cm- clm,  p_cdm(hx))



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Detritus (if eaten): calculate available amount
  ! and add it to the total amount of food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  availQ6_c  =   Q6c(:)* PartQ_vector(  D6m(:),  clm,  chm,  p_d_tot)
  availQ6_n  =   Q6n(:)* PartQ_vector(  D7m(:),  clm,  chm,  p_d_tot)
  availQ6_p  =   Q6p(:)* PartQ_vector(  D8m(:),  clm,  chm,  p_d_tot)


  qnQ6c  =   availQ6_n/( availQ6_c+ 1.0D-30)
  qpQ6c  =   availQ6_p/( availQ6_c+ 1.0D-30)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Growth is controlled by quality of detritus (N and P content):
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  eN  =   eramp_vector(  qnQ6c,  p_qnc(hx))* eramp_vector(  qpQ6c,  p_qpc(hx))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total substrate availability:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruQ6c  =   p_suhQ6(hx)* availQ6_c* eN+ p_sulQ6(hx)* availQ6_c

  ruQ1c  =   p_suQ1(hx)* BenDetritus(p_iQ1(hx),iiC)

  rut  =   ruQ6c+ ruQ1c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential uptake by bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rum  =   p_sum(hx)* et* eo* hxc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Actual uptake by bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rug  =   min(  rum,  rut)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon fluxes into bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruQ6c  =   rug* ruQ6c/ rut
  ruQ1c  =   rug* ruQ1c/ rut

  call flux_vector( iiBen, ppQ6c,pphxc, ruQ6c )
  call flux_vector( iiBen, ppBenDetritus(p_iQ1(hx),iiC),pphxc, ruQ1c )
  ruHI(hx,:)  =   ruQ1c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient fluxes into bacteria from carbon fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruQ6n  =   ruQ6c* qnQ6c
  ruQ6p  =   ruQ6c* qpQ6c
  ruQ1n  =   p_suQ1(hx)* BenDetritus(p_iQ1(hx),iiN)* rug/ rut
  ruQ1p  =   p_suQ1(hx)* BenDetritus(p_iQ1(hx),iiP)* rug/ rut


  call flux_vector( iiBen, ppQ6n,pphxn, ruQ6n )
  call flux_vector( iiBen, ppQ6p,pphxp, ruQ6p )


  call flux_vector( iiBen, ppBenDetritus(p_iQ1(hx),iiN),pphxn, ruQ1n )
  call flux_vector( iiBen, ppBenDetritus(p_iQ1(hx),iiP),pphxp, ruQ1p )


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of respiration:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrc  =   p_srr(hx)* hxc* et+( ruQ1c+ ruQ6c)* p_pur(hx)

  run  =   max(  0.0D+00,  rug- rrc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of potential nutrient uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruKIn  =   BenthicAmmonium(p_iK4(hx),iiN)* p_sumKIn(hx)
  ruKIp  =   BenthicPhosphate(p_iK1(hx),iiP)* p_sumKIp(hx)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon correction: all C which cannot be used for growth due to
  ! lack of nutrients is excreted!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  netgrowth  =   min(  run, ( p_small+ ruQ6n+ ruQ1n+ ruKIn)/ p_qlnc(hx))
  netgrowth  =   min(  netgrowth, ( ruQ6p+ ruQ1p+ ruKIp)/ p_qlpc(hx))
  netgrowth  =   max(  netgrowth,  0.0D+00)
  call flux_vector( iiBen, pphxc,pphxc,-( run- netgrowth) )
  run  =   netgrowth

  ren  =   run*(( ruQ6n+ ruQ1n)/( p_small+ run)- p_qnc(hx))
  rep  =   run*(( ruQ6p+ ruQ1p)/( p_small+ run)- p_qpc(hx))


  r  =   insw_vector(  ren)
  call flux_vector( iiBen, pphxn,ppBenthicAmmonium(p_iK4(hx),iiN), ren* r )

  ! shortage of nutrients : ren < 0 --> Nutrient uptake
  runn  =   min( - ren,  ruKIn)*( 1.0D+00- r)
  call flux_vector(iiBen, ppBenthicAmmonium(p_iK4(hx),iiN),pphxn, runn)

  r  =   insw_vector(  rep)
  call flux_vector( iiBen, pphxp,ppBenthicPhosphate(p_iK1(hx),iiP), rep* r )

  ! shortage of nutrients : rep < 0 --> Nutrient uptake
  runp  =   min( - rep,  ruKIp)*( 1.0D+00- r)
  call flux_vector(iiBen, ppBenthicPhosphate(p_iK1(hx),iiP),pphxp, runp)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  sm  =   p_sd(hx)*( 1.0D+00- eo)

  rqt6c  =   hxc* sm*( 1.0D+00- p_pe_R1c)
  rqt6n  =   hxn* sm*( 1.0D+00- p_pe_R1n)
  rqt6p  =   hxp* sm*( 1.0D+00- p_pe_R1p)

  call flux_vector( iiBen, pphxc,ppBenDetritus(p_iQ1(hx),iiC), hxc* sm* p_pe_R1c &
    )
  call flux_vector( iiBen, pphxn,ppBenDetritus(p_iQ1(hx),iiN), hxn* sm* p_pe_R1n &
    )
  call flux_vector( iiBen, pphxp,ppBenDetritus(p_iQ1(hx),iiP), hxp* sm* p_pe_R1p &
    )

  call flux_vector( iiBen, pphxc,ppQ6c, rqt6c )
  call flux_vector( iiBen, pphxn,ppQ6n, rqt6n )
  call flux_vector( iiBen, pphxp,ppQ6p, rqt6p )



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Assigning depends on type of bacteria and layers in which they occur:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( hx)

    case ( iiH1 )

      call flux_vector( iiBen, pphxc,pphxc,-( rrc) )
      call flux_vector(iiBen, ppG2o,ppG2o,-( rrc/ 12.0D+00))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Add respiration and excretion to benthic totals:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      rrBTo(:)  =   rrBTo(:)+ rrc/ 12.0D+00
      reBTn(:)  =   reBTn(:)+ ren
      reBTp(:)  =   reBTp(:)+ rep





    case ( iiH2 )

      call flux_vector( iiBen, pphxc,pphxc,-( rrc) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Respiration in anoxic circumstances produces reduced material
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call flux_vector( iiBen, ppK6r,ppK6r, rrc/ 12.0D+00* p_qro )


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Add respiration and excretion to benthic totals:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      rrATo(:)  =   rrATo(:)+ rrc/ 12.0D+00
      reATn(:)  =   reATn(:)+ ren
      reATp(:)  =   reATp(:)+ rep




  end select



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of changes due to uptake of detritus in distribution
  ! of state variables (Dx.m is a undetermined source)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector(iiBen, ppD6m,ppD6m,( cmm- D6m(:))*( ruQ6c- rqt6c)/ Q6c(:))
  call flux_vector(iiBen, ppD7m,ppD7m,( cmm- D7m(:))*( ruQ6n- rqt6n)/ Q6n(:))
  call flux_vector(iiBen, ppD8m,ppD8m,( cmm- D8m(:))*( ruQ6p- rqt6p)/ Q6p(:))





  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
