#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MicroZoo
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine MicroZooDynamics(zoo,  ppzooc, ppzoon, ppzoop)
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: B1c, B1n, B1p, O2o, &
  ! R1c, R6c, R1n, R6n, R1p, R6p, N4n, N1p
  ! For the following Pelagic-group-states fluxes are defined: &
  ! PhytoPlankton, MicroZooPlankton
  ! The following Pelagic 1-d global boxvars are modified : flP1R6s
  ! The following Pelagic 1-d global boxvars are used: ETW, eO2mO2, qnB1c, &
  ! qpB1c
  ! The following Pelagic 2-d global boxvars are used: qnPc, qpPc, qn_mz, qp_mz, &
  ! qlPc, qsPc
  ! The following groupmember vars are used: iiPhytoPlankton, &
  ! iiMicroZooPlankton, iiP1
  ! The following constituent constants  are used: iiC, iiN, iiP, iiL
  ! The following 0-d global box parametes are used: p_pe_R1c, p_pe_R1n, &
  ! p_pe_R1p, p_small
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem, ONLY: D3STATE, B1c, B1n, B1p, O2o, R1c, R6c, R1n, R6n, &
    R1p, R6p, N4n, N1p, PhytoPlankton, MicroZooPlankton
#endif
  use mem, ONLY: ppB1c, ppB1n, ppB1p, ppO2o, ppR1c, ppR6c, &
    ppR1n, ppR6n, ppR1p, ppR6p, ppN4n, ppN1p, ppPhytoPlankton, ppMicroZooPlankton, &
    flP1R6s, ETW, eO2mO2, qnB1c, qpB1c, qnPc, qpPc, qn_mz, qp_mz, &
    qlPc, qsPc, iiPhytoPlankton, iiMicroZooPlankton, iiP1, iiC, iiN, iiP, iiL, &
    NO_BOXES, iiBen, iiPel, flux_vector
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, p_small
  use mem_MicroZoo


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
!   ERSEM group, Hanneke Baretta-Bekker
!
!
!
! !REVISION_HISTORY
!   by Piet Ruardij at Thu Mar 16 08:34:04 CET 2006
!       s: BFMI
!       d: 
!	

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
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES) :: zooc
  real(RLEN),dimension(NO_BOXES) :: zoon
  real(RLEN),dimension(NO_BOXES) :: zoop
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  real(RLEN),dimension(NO_BOXES)  :: CORROX
  real(RLEN),dimension(NO_BOXES)  :: put_u
  real(RLEN),dimension(NO_BOXES)  :: et
  real(RLEN),dimension(NO_BOXES)  :: eO2
  real(RLEN),dimension(NO_BOXES)  :: rumc
  real(RLEN),dimension(NO_BOXES)  :: rumn
  real(RLEN),dimension(NO_BOXES)  :: rump
  real(RLEN),dimension(NO_BOXES)  :: rugc
  real(RLEN),dimension(NO_BOXES)  :: rugn
  real(RLEN),dimension(NO_BOXES)  :: rugp
  real(RLEN),dimension(NO_BOXES)  :: runc
  real(RLEN),dimension(NO_BOXES)  :: runn
  real(RLEN),dimension(NO_BOXES)  :: runp
  real(RLEN),dimension(NO_BOXES)  :: efood
  real(RLEN),dimension(NO_BOXES)  :: rrsc
  real(RLEN),dimension(NO_BOXES)  :: rrac
  real(RLEN),dimension(NO_BOXES)  :: reac
  real(RLEN),dimension(NO_BOXES)  :: rdc
  real(RLEN),dimension(NO_BOXES)  :: rrtc
  real(RLEN),dimension(NO_BOXES)  :: ruB1c
  real(RLEN),dimension(NO_BOXES)  :: ruPIc
  real(RLEN),dimension(NO_BOXES)  :: ruZIc
  real(RLEN),dimension(NO_BOXES)  :: rumB1c
  real(RLEN),dimension(NO_BOXES,4)  :: rumPIc
  real(RLEN),dimension(NO_BOXES,2)  :: rumZIc
  real(RLEN),dimension(NO_BOXES)  :: rric
  real(RLEN),dimension(NO_BOXES)  :: rr1c
  real(RLEN),dimension(NO_BOXES)  :: rr6c
  real(RLEN),dimension(NO_BOXES)  :: rr1p
  real(RLEN),dimension(NO_BOXES)  :: rr1n
  real(RLEN),dimension(NO_BOXES)  :: rrip
  real(RLEN),dimension(NO_BOXES)  :: rr6p
  real(RLEN),dimension(NO_BOXES)  :: rep
  real(RLEN),dimension(NO_BOXES)  :: rrin
  real(RLEN),dimension(NO_BOXES)  :: rr6n
  real(RLEN),dimension(NO_BOXES)  :: ren
  real(RLEN),dimension(NO_BOXES)  :: pu_ra
  real(RLEN),dimension(NO_BOXES)  :: r
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  zooc = D3STATE(ppzooc,:)
  zoon = D3STATE(ppzoon,:)
  zoop = D3STATE(ppzoop,:)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq_vector(  ETW(:),  p_q10(zoo))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Oxygen limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  CORROX  =   1.0D+00+ p_chro(zoo)
  eO2  =   min(  1.0D+00,  CORROX* MM_vector(  eO2mO2(:),   p_chro(zoo)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Available food, etc...
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rumB1c  =   p_suB1(zoo)* B1c(:)* B1c(:)/( B1c(:)+ p_minfood(zoo))
  rumc  =   rumB1c
  rumn  =   rumB1c* qnB1c(:)
  rump  =   rumB1c* qpB1c(:)

  do i = 1 , ( iiPhytoPlankton)

    rumPIc(:, i) = p_suPI(zoo,i)* PhytoPlankton(i,iiC)* &
      PhytoPlankton(i,iiC)/( PhytoPlankton(i,iiC)+ p_minfood(zoo))
    rumc  =   rumc+ rumPIc(:, i)
    rumn  =   rumn+ rumPIc(:, i)* qnPc(i,:)
    rump  =   rump+ rumPIc(:, i)* qpPc(i,:)
  end do


  do i = 1 , ( iiMicroZooPlankton)

    rumZIc(:, i) = p_suZI(zoo,i)* &
      MicroZooPlankton(i,iiC)* MicroZooPlankton(i,iiC)/( MicroZooPlankton(i,iiC)+ &
      p_minfood(zoo))
    rumc  =   rumc+ rumZIc(:, i)
    rumn  =   rumn+ rumZIc(:, i)* qn_mz(i,:)
    rump  =   rump+ rumZIc(:, i)* qp_mz(i,:)
  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  efood  =   MM_vector(  rumc,  p_chuc(zoo))
  rugc  =   p_sum(zoo)* et* zooc* efood

  r  =   min(  rumn/ p_qn_mz(zoo),  rump/ p_qp_mz(zoo))

  pu_ra  =   max(  p_pu_ra(zoo),  1.0D+00- r/ rumc)

  put_u  =   rugc/ rumc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes into microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruB1c  =   put_u* rumB1c
  call flux_vector( iiPel, ppB1c,ppzooc, ruB1c )
  call flux_vector( iiPel, ppB1n,ppzoon, ruB1c* qnB1c(:) )
  call flux_vector( iiPel, ppB1p,ppzoop, ruB1c* qpB1c(:) )
  rugn  =   ruB1c* qnB1c(:)
  rugp  =   ruB1c* qpB1c(:)

  do i = 1 , ( iiPhytoPlankton)

    ruPIc  =   put_u* rumPIc(:, i)
    call flux_vector( iiPel, ppPhytoPlankton(i,iiC),ppzooc, ruPIc )
    call flux_vector( iiPel, ppPhytoPlankton(i,iiN),ppzoon, ruPIc* qnPc(i,:) )
    call flux_vector( iiPel, ppPhytoPlankton(i,iiP),ppzoop, ruPIc* qpPc(i,:) )
    ! Chl is transferred to the sink
    call flux_vector( iiPel, ppPhytoPlankton(i,iiL),ppPhytoPlankton(i,iiL),-( &
      ruPIc* qlPc(i,:)) )
    if ( i== iiP1) then
      ! P1s is directly transferred to R6s
      ! PhytoPlankton[i].s -> R6.s = ruPIc * qsPc[i];
      flP1R6s(:)  =   flP1R6s(:)+ ruPIc* qsPc(i,:)
    end if

    rugn  =   rugn+ ruPIc* qnPc(i,:)
    rugp  =   rugp+ ruPIc* qpPc(i,:)
  end do


  do i = 1 , ( iiMicroZooPlankton)

    ruZIc  =   put_u* rumZIc(:, i)
    ! intra-group predation is not computed
    if ( i/= zoo) then
      call flux_vector( iiPel, ppMicroZooPlankton(i,iiC),ppzooc, ruZIc )
      call flux_vector( iiPel, ppMicroZooPlankton(i,iiN),ppzoon, ruZIc* &
        qn_mz(i,:) )
      call flux_vector( iiPel, ppMicroZooPlankton(i,iiP),ppzoop, ruZIc* &
        qp_mz(i,:) )
    end if

    rugn  =   rugn+ ruZIc* qn_mz(i,:)
    rugp  =   rugp+ ruZIc* qp_mz(i,:)
  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !       Fluxes from microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Rest, activity, total respiration fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrsc  =   p_srs(zoo)* et* zooc
  rrac  =   rugc* pu_ra
  rrtc  =   rrsc+ rrac

  call flux_vector( iiPel, ppzooc,ppzooc,-( rrtc) )
  call flux_vector( iiPel, ppO2o,ppO2o,-( rrtc/ 12.0D+00) )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Mortality (rdc) + Excetion (reac)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rdc  =  (( 1.0D+00- eO2)* p_sdo(zoo)+ p_sd(zoo))* zooc
  reac  =   rugc* p_pu_ea(zoo)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes due to mortality and excetion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rric  =  ( reac+ rdc)
  rr1c  =   rric* p_pe_R1c
  rr6c  =   rric*( 1.0D+00- p_pe_R1c)

  call flux_vector( iiPel, ppzooc,ppR1c, rr1c )
  call flux_vector( iiPel, ppzooc,ppR6c, rr6c )
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !     Nutrient dynamics in microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrin  =   rugn* p_pu_ea(zoo)+ rdc* qn_mz(zoo,:)
  rr1n  =   rrin* p_pe_R1n
  rr6n  =   rrin- rr1n

  call flux_vector( iiPel, ppzoon,ppR1n, rr1n )
  call flux_vector( iiPel, ppzoon,ppR6n, rr6n )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Phosphorus dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrip  =   rugp* p_pu_ea(zoo)+ rdc* qp_mz(zoo,:)
  rr1p  =   rrip* p_pe_R1p
  rr6p  =   rrip- rr1p

  call flux_vector( iiPel, ppzoop,ppR1p, rr1p )
  call flux_vector( iiPel, ppzoop,ppR6p, rr6p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Dissolved nutrient dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  runc  =   max(  0.0D+00,  rugc*( 1.0D+00- p_pu_ea(zoo))- rrac)
  runn  =   max(  0.0D+00,  rugn*( 1.0D+00- p_pu_ea(zoo))+ rrsc* qn_mz(zoo, :))
  runp  =   max(  0.0D+00,  rugp*( 1.0D+00- p_pu_ea(zoo))+ rrsc* qp_mz(zoo, :))

  ren  =   max(  0.0D+00,  runn/( p_small+ runc)- p_qn_mz(zoo))* runc
  rep  =   max(  0.0D+00,  runp/( p_small+ runc)- p_qp_mz(zoo))* runc
  call flux_vector( iiPel, ppzoon,ppN4n, ren )
  call flux_vector( iiPel, ppzoop,ppN1p, rep )







  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
