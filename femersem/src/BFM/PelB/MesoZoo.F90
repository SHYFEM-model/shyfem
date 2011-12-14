#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MesoZoo
!
! DESCRIPTION
!   This submodel describes the carbon dynamics and associated
!    nutrient dynamics in carnivorous mesozooplankton (represented
!    by the state variable Z3) and in omnivorous zooplankton (in
!    the model known as Z4).
!    
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine MesoZooDynamics(zoo,  ppzooc, ppzoon, ppzoop)
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: O2o, N1p, N4n, R6c, &
  ! R6p, R6n
  ! For the following Pelagic-group-states fluxes are &
  ! defined: PhytoPlankton, MicroZooPlankton, MesoZooPlankton
  ! The following Pelagic 1-d global boxvars are modified : flP1R6s
  ! The following Pelagic 1-d global boxvars  are used: ETW
  ! The following Pelagic 2-d global boxvars are used: qnPc, qpPc, qlPc, qsPc, &
  ! qn_mz, qp_mz, qnZc, qpZc
  ! The following groupmember vars are used: iiPhytoPlankton, &
  ! iiMicroZooPlankton, iiMesoZooPlankton, iiP1
  ! The following constituent constants  are used: iiC, iiN, iiP, iiL
  ! The following 0-d global box parametes are used: p_small
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem, ONLY: D3STATE, O2o, N1p, N4n, R6c, R6p, &
    R6n, PhytoPlankton, MicroZooPlankton, MesoZooPlankton
#endif
  use mem, ONLY: ppO2o, ppN1p, ppN4n, ppR6c, ppR6p, &
    ppR6n, ppPhytoPlankton, ppMicroZooPlankton, ppMesoZooPlankton, flP1R6s, ETW, &
    qnPc, qpPc, qlPc, qsPc, qn_mz, qp_mz, qnZc, qpZc, iiPhytoPlankton, &
    iiMicroZooPlankton, iiMesoZooPlankton, iiP1, iiC, iiN, iiP, iiL, NO_BOXES, &
    iiBen, iiPel, flux_vector
  use mem_Param,  ONLY: p_small
  use mem_MesoZoo


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, MM_vector


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: zoo
  integer,intent(IN) :: ppzooc
  integer,intent(IN) :: ppzoon
  integer,intent(IN) :: ppzoop

!  
!
! !AUTHORS
!   N. Broekhuizen and A.D. Bryant, ERSEM group
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
  real(RLEN),dimension(NO_BOXES) :: zooc
  real(RLEN),dimension(NO_BOXES) :: zoon
  real(RLEN),dimension(NO_BOXES) :: zoop
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  integer,dimension(NO_BOXES)  :: nut_lim
  real(RLEN),dimension(NO_BOXES)  :: put_u
  real(RLEN),dimension(NO_BOXES)  :: temp_p
  real(RLEN),dimension(NO_BOXES)  :: temp_n
  real(RLEN),dimension(NO_BOXES)  :: rumc
  real(RLEN),dimension(NO_BOXES)  :: rugc
  real(RLEN),dimension(NO_BOXES)  :: et
  real(RLEN),dimension(NO_BOXES)  :: rrs_c
  real(RLEN),dimension(NO_BOXES)  :: rrs_n
  real(RLEN),dimension(NO_BOXES)  :: rrs_p
  real(RLEN),dimension(NO_BOXES)  :: rra_c
  real(RLEN),dimension(NO_BOXES)  :: rra_n
  real(RLEN),dimension(NO_BOXES)  :: rra_p
  real(RLEN),dimension(NO_BOXES)  :: rut_c
  real(RLEN),dimension(NO_BOXES)  :: rut_n
  real(RLEN),dimension(NO_BOXES)  :: rut_p
  real(RLEN),dimension(NO_BOXES)  :: rd_c
  real(RLEN),dimension(NO_BOXES)  :: rd_n
  real(RLEN),dimension(NO_BOXES)  :: rd_p
  real(RLEN),dimension(NO_BOXES)  :: sdo
  real(RLEN),dimension(NO_BOXES)  :: rdo_c
  real(RLEN),dimension(NO_BOXES)  :: rdo_n
  real(RLEN),dimension(NO_BOXES)  :: rdo_p
  real(RLEN),dimension(NO_BOXES)  :: ret_c
  real(RLEN),dimension(NO_BOXES)  :: ret_n
  real(RLEN),dimension(NO_BOXES)  :: ret_p
  real(RLEN),dimension(NO_BOXES)  :: ru_c
  real(RLEN),dimension(NO_BOXES)  :: ru_n
  real(RLEN),dimension(NO_BOXES)  :: ru_p
  real(RLEN),dimension(NO_BOXES)  :: pu_e_n
  real(RLEN),dimension(NO_BOXES)  :: pu_e_p
  real(RLEN),dimension(NO_BOXES)  :: prI_R6
  real(RLEN),dimension(NO_BOXES)  :: pe_R6c
  real(RLEN),dimension(NO_BOXES)  :: pe_N1p
  real(RLEN),dimension(NO_BOXES)  :: pe_N4n
  real(RLEN),dimension(NO_BOXES,4)  :: rumPIc
  real(RLEN),dimension(NO_BOXES,2)  :: rumMIZc
  real(RLEN),dimension(NO_BOXES,2)  :: rumMEZc
  real(RLEN),dimension(NO_BOXES)  :: ruPIc
  real(RLEN),dimension(NO_BOXES)  :: ruMIZc
  real(RLEN),dimension(NO_BOXES)  :: ruMEZc
  real(RLEN),dimension(NO_BOXES)  :: rq6c
  real(RLEN),dimension(NO_BOXES)  :: rq6n
  real(RLEN),dimension(NO_BOXES)  :: rq6p
  real(RLEN),dimension(NO_BOXES)  :: rrc
  real(RLEN),dimension(NO_BOXES)  :: ren
  real(RLEN),dimension(NO_BOXES)  :: rep
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  zooc = D3STATE(ppzooc,:)
  zoon = D3STATE(ppzoon,:)
  zoop = D3STATE(ppzoop,:)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !Physiological temperature response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq_vector(  ETW(:),  p_q10(zoo))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! rua: food uptake generalized. Addition of new FG become much more simpler!


  rumc  =   0.0D+00
  do i = 1 , ( iiPhytoPlankton)

    rumPIc(:, i)  =   p_puPI(zoo,i)* PhytoPlankton(i,iiC)
    rumc  =   rumc+ rumPIc(:, i)
  end do


  do i = 1 , ( iiMicroZooPlankton)

    rumMIZc(:, i)  =   p_puMIZ(zoo,i)* MicroZooPlankton(i,iiC)
    rumc  =   rumc+ rumMIZc(:, i)
  end do


  do i = 1 , ( iiMesoZooPlankton)

    rumMEZc(:, i)  =   p_puMEZ(zoo,i)* MesoZooPlankton(i,iiC)
    rumc  =   rumc+ rumMEZc(:, i)
  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rugc  =   et* p_sum(zoo)* MM_vector(  p_vum(zoo)* rumc,  p_sum(zoo))* zooc
  put_u  =   rugc/ rumc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total Gross Uptakes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rut_c  =   0.0D+00
  rut_n  =   0.0D+00
  rut_p  =   0.0D+00
  do i = 1 , ( iiPhytoPlankton)

    ruPIc  =   put_u* rumPIc(:, i)
    call flux_vector( iiPel, ppPhytoPlankton(i,iiC),ppzooc, ruPIc )
    call flux_vector( iiPel, ppPhytoPlankton(i,iiN),ppzoon, ruPIc* qnPc(i,:) )
    call flux_vector( iiPel, ppPhytoPlankton(i,iiP),ppzoop, ruPIc* qpPc(i,:) )
    rut_c  =   rut_c+ ruPIc
    rut_n  =   rut_n+ ruPIc* qnPc(i,:)
    rut_p  =   rut_p+ ruPIc* qpPc(i,:)
    ! Chl is transferred to the sink
    call flux_vector( iiPel, ppPhytoPlankton(i,iiL),ppPhytoPlankton(i,iiL),-( &
      ruPIc* qlPc(i,:)) )
    if ( i== iiP1) then
      ! P1s is directly transferred to R6s
      ! PhytoPlankton[i].s -> R6.s = ruPIc * qsPc[i];
      flP1R6s(:)  =   flP1R6s(:)+ ruPIc* qsPc(i,:)
    end if

  end do


  do i = 1 , ( iiMicroZooPlankton)

    ruMIZc  =   put_u* rumMIZc(:, i)
    call flux_vector( iiPel, ppMicroZooPlankton(i,iiC),ppzooc, ruMIZc )
    call flux_vector( iiPel, ppMicroZooPlankton(i,iiN),ppzoon, ruMIZc* &
      qn_mz(i,:) )
    call flux_vector( iiPel, ppMicroZooPlankton(i,iiP),ppzoop, ruMIZc* &
      qp_mz(i,:) )
    rut_c  =   rut_c+ ruMIZc
    rut_n  =   rut_n+ ruMIZc* qn_mz(i,:)
    rut_p  =   rut_p+ ruMIZc* qp_mz(i,:)
  end do


  do i = 1 , ( iiMesoZooPlankton)

    ruMEZc  =   put_u* rumMEZc(:, i)
    ! intra-group predation is not computed
    if ( i/= zoo) then
      call flux_vector( iiPel, ppMesoZooPlankton(i,iiC),ppzooc, ruMEZc )
      call flux_vector( iiPel, ppMesoZooPlankton(i,iiN),ppzoon, ruMEZc* &
        qnZc(i,:) )
      call flux_vector( iiPel, ppMesoZooPlankton(i,iiP),ppzoop, ruMEZc* &
        qpZc(i,:) )
    end if

    rut_c  =   rut_c+ ruMEZc
    rut_n  =   rut_n+ ruMEZc* qnZc(i,:)
    rut_p  =   rut_p+ ruMEZc* qpZc(i,:)
  end do



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Proportion of ingested food respired by zoo = prIR6/Z4R6
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  prI_R6  =   1.0D+00- p_puI_u(zoo)- p_peI_R6(zoo)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Assimilated material
  ! Respectively Carbon, Nitrogen and Phosphorus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ru_c  =   p_puI_u(zoo)* rut_c
  ru_n  =  ( p_puI_u(zoo)+ prI_R6)* rut_n
  ru_p  =  ( p_puI_u(zoo)+ prI_R6)* rut_p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! P:C and N:C ratios in assimilate
  ! Nitrogen & Phosphorus:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pu_e_n  =   ru_n/( p_small+ ru_c)
  pu_e_p  =   ru_p/( p_small+ ru_c)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Determine whether C, P or N is the Limiting Nutrient. Variable nut_lim
  ! holds the kind of nutrient limitation (1, 2, 3).
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  nut_lim  =   1

  temp_p  =   pu_e_p/ qpZc(zoo,:)
  temp_n  =   pu_e_n/ qnZc(zoo,:)



    WHERE ( (( temp_p< temp_n) .OR.( abs(temp_p- temp_n)< p_small)) )
      where ( pu_e_p< qpZc(zoo,:))
        nut_lim  =   2
      end where




    ELSEWHERE
      where ( pu_e_n< qnZc(zoo,:))
        nut_lim  =   3
      end where




  END WHERE


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration and basal metabolism
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rra_c  =   prI_R6* rut_c
  rra_n  =   0.0D+00
  rra_p  =   0.0D+00

  rrs_c  =   p_srs(zoo)* et* zooc
  rrs_n  =   p_srs(zoo)* et* zoon
  rrs_p  =   p_srs(zoo)* et* zoop

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Defecation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ret_c  =   p_peI_R6(zoo)* rut_c
  ret_n  =   p_peI_R6(zoo)* rut_n
  ret_p  =   p_peI_R6(zoo)* rut_p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Natural mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rd_c  =   p_sd(zoo)* et* zooc
  rd_n  =   p_sd(zoo)* et* zoon
  rd_p  =   p_sd(zoo)* et* zoop

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Density dependent mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sdo  =   p_sdo(zoo)* (zooc)**(p_sds(zoo))
  rdo_c  =   sdo* zooc
  rdo_n  =   sdo* zoon
  rdo_p  =   sdo* zoop

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Eliminate excess of non-limiting nutrients
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


    WHERE (( nut_lim)==1)
      pe_R6c  =   0.0D+00
      pe_N1p = max( 0.0D+00, ( 1.0D+00- p_peI_R6(zoo))* rut_p- p_qpc(zoo)* &
        ru_c)
      pe_N1p  =   pe_N1p/( p_small+ rut_p)
      pe_N4n = max( 0.0D+00, ( 1.0D+00- p_peI_R6(zoo))* rut_n- p_qnc(zoo)* &
        ru_c)
      pe_N4n  =   pe_N4n/( p_small+ rut_n)



    ELSEWHERE (( nut_lim)==2)
      pe_N1p  =   0.0D+00
      pe_R6c  =  max(0.0,( p_qpc(zoo)* ru_c)-( 1.0D+00- p_peI_R6(zoo))* rut_p)
      pe_R6c  =   pe_R6c/( p_small+ p_qpc(zoo)* rut_c)
      pe_N4n = max( 0.0D+00, ( 1.0D+00- p_peI_R6(zoo))* rut_n- p_qnc(zoo)*( &
        ru_c- pe_R6c* rut_c))
      pe_N4n  =   pe_N4n/( p_small+ rut_n)



    ELSEWHERE (( nut_lim)==3)
      pe_N4n  =   0.0D+00
      pe_R6c  = max(0.0, ( p_qnc(zoo)* ru_c)-( 1.0D+00- p_peI_R6(zoo))* rut_n)
      pe_R6c  =   pe_R6c/( p_small+ p_qnc(zoo)* rut_c)
      pe_N1p = max( 0.0D+00, ( 1.0D+00- p_peI_R6(zoo))* rut_p- p_qpc(zoo)*( &
        ru_c- pe_R6c* rut_c))
      pe_N1p  =   pe_N1p/( p_small+ rut_p)



  END WHERE
  ! End of select(nut_lim)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes for eliminated excess nutrients
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rq6c  =   rd_c+ ret_c+ rdo_c+ pe_R6c* rut_c
  rq6p  =   rd_p+ ret_p+ rdo_p
  rq6n  =   rd_n+ ret_n+ rdo_n

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration: activity + basal metabolism
  ! Excretion: activity + basal metabolism + excess non-limiting nutrients
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrc  =   rra_c+ rrs_c
  ren  =   rra_n+ rrs_n+ pe_N4n* rut_n
  rep  =   rra_p+ rrs_p+ pe_N1p* rut_p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! flow statements
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call flux_vector( iiPel, ppO2o,ppO2o,-( rrc/ 12.0D+00) )
  call flux_vector( iiPel, ppzooc,ppzooc,-( rrc) )
  call flux_vector( iiPel, ppzoop,ppN1p, rep )
  call flux_vector( iiPel, ppzoon,ppN4n, ren )

  call flux_vector( iiPel, ppzooc,ppR6c, rq6c )
  call flux_vector( iiPel, ppzoop,ppR6p, rq6p )
  call flux_vector( iiPel, ppzoon,ppR6n, rq6n )




  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
