!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.50-g
!
! SUBROUTINE
!   AllocateMem
!
! FILE
!   AllocateMem
!
! DESCRIPTION
!   Allocation of memory for Global State variables and other Variables
!  
!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! AUTHORS
!   mfstep/ERSEM team
!
! CHANGE_LOG
!   ---
!
! COPYING
!   
!   Copyright (C) 2004 P. Ruardij, the mfstep group, the ERSEM team 
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
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!
  subroutine AllocateMem
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  use (import) other modules
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem
  use mem
  use mem_Param
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  integer:: status
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Start the allocation of pelagic state global
  ! matrix and pointers
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-"
  write(LOGUNIT,*) "# Allocating State Variables and Rates array ..."

#ifndef NOT_STANDALONE
 
          allocate(D3STATE(1:NO_D3_BOX_STATES,1:NO_BOXES),stat=status)
          if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D3STATE")
          D3STATE = ZERO
          allocate(D3SOURCE(1:NO_D3_BOX_STATES,1:NO_D3_BOX_STATES,1:NO_BOXES),stat=status)
          if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D3SOURCE")
          D3SOURCE = ZERO
          allocate(D3SINK(1:NO_D3_BOX_STATES,1:NO_D3_BOX_STATES,1:NO_BOXES) ,stat=status)
          if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D3SINK")
          D3SINK = ZERO
          allocate(D3STATETYPE(1:NO_D3_BOX_STATES ),stat=status)
          if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem","D3STATETYPE")
          D3STATETYPE = ZERO

 
          allocate(D2STATE(1:NO_D2_BOX_STATES,1:NO_BOXES_XY),stat=status)
          if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D2STATE")
          D2STATE = ZERO
          allocate(D2SOURCE(1:NO_D2_BOX_STATES,1:NO_D2_BOX_STATES,1:NO_BOXES_XY),stat=status)
          if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D2SOURCE")
          D2SOURCE = ZERO
          allocate(D2SINK(1:NO_D2_BOX_STATES,1:NO_D2_BOX_STATES,1:NO_BOXES_XY) ,stat=status)
          if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D2SINK")
          D2SINK = ZERO
          allocate(D2STATETYPE(1:NO_D2_BOX_STATES ),stat=status)
          if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem","D2STATETYPE")
          D2STATETYPE = ZERO

 
          allocate(D3DIAGNOS(1:NO_D3_BOX_DIAGNOSS,1:NO_BOXES),stat=status)
          if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D3DIAGNOS")
          D3DIAGNOS = ZERO

 
          allocate(D2DIAGNOS(1:NO_D2_BOX_DIAGNOSS,1:NO_BOXES_XY),stat=status)
          if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D2DIAGNOS")
          D2DIAGNOS = ZERO

#endif




  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocation of Pelagic variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        O2o => D3STATE(ppO2o,:); O2o=ZERO
        N1p => D3STATE(ppN1p,:); N1p=ZERO
        N3n => D3STATE(ppN3n,:); N3n=ZERO
        N4n => D3STATE(ppN4n,:); N4n=ZERO
        O4n => D3STATE(ppO4n,:); O4n=ZERO
        N5s => D3STATE(ppN5s,:); N5s=ZERO
        N6r => D3STATE(ppN6r,:); N6r=ZERO
        B1c => D3STATE(ppB1c,:); B1c=ZERO
        B1n => D3STATE(ppB1n,:); B1n=ZERO
        B1p => D3STATE(ppB1p,:); B1p=ZERO
        P1c => D3STATE(ppP1c,:); P1c=ZERO
        P1n => D3STATE(ppP1n,:); P1n=ZERO
        P1p => D3STATE(ppP1p,:); P1p=ZERO
        P1l => D3STATE(ppP1l,:); P1l=ZERO
        P1s => D3STATE(ppP1s,:); P1s=ZERO
        P2c => D3STATE(ppP2c,:); P2c=ZERO
        P2n => D3STATE(ppP2n,:); P2n=ZERO
        P2p => D3STATE(ppP2p,:); P2p=ZERO
        P2l => D3STATE(ppP2l,:); P2l=ZERO
        P3c => D3STATE(ppP3c,:); P3c=ZERO
        P3n => D3STATE(ppP3n,:); P3n=ZERO
        P3p => D3STATE(ppP3p,:); P3p=ZERO
        P3l => D3STATE(ppP3l,:); P3l=ZERO
        P4c => D3STATE(ppP4c,:); P4c=ZERO
        P4n => D3STATE(ppP4n,:); P4n=ZERO
        P4p => D3STATE(ppP4p,:); P4p=ZERO
        P4l => D3STATE(ppP4l,:); P4l=ZERO
        Z3c => D3STATE(ppZ3c,:); Z3c=ZERO
        Z3n => D3STATE(ppZ3n,:); Z3n=ZERO
        Z3p => D3STATE(ppZ3p,:); Z3p=ZERO
        Z4c => D3STATE(ppZ4c,:); Z4c=ZERO
        Z4n => D3STATE(ppZ4n,:); Z4n=ZERO
        Z4p => D3STATE(ppZ4p,:); Z4p=ZERO
        Z5c => D3STATE(ppZ5c,:); Z5c=ZERO
        Z5n => D3STATE(ppZ5n,:); Z5n=ZERO
        Z5p => D3STATE(ppZ5p,:); Z5p=ZERO
        Z6c => D3STATE(ppZ6c,:); Z6c=ZERO
        Z6n => D3STATE(ppZ6n,:); Z6n=ZERO
        Z6p => D3STATE(ppZ6p,:); Z6p=ZERO
        R1c => D3STATE(ppR1c,:); R1c=ZERO
        R1n => D3STATE(ppR1n,:); R1n=ZERO
        R1p => D3STATE(ppR1p,:); R1p=ZERO
        R2c => D3STATE(ppR2c,:); R2c=ZERO
        R6c => D3STATE(ppR6c,:); R6c=ZERO
        R6n => D3STATE(ppR6n,:); R6n=ZERO
        R6p => D3STATE(ppR6p,:); R6p=ZERO
        R6s => D3STATE(ppR6s,:); R6s=ZERO
        R7c => D3STATE(ppR7c,:); R7c=ZERO


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocation of Benthic variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        Y1c => D2STATE(ppY1c,:); Y1c=ZERO
        Y1n => D2STATE(ppY1n,:); Y1n=ZERO
        Y1p => D2STATE(ppY1p,:); Y1p=ZERO
        Y2c => D2STATE(ppY2c,:); Y2c=ZERO
        Y2n => D2STATE(ppY2n,:); Y2n=ZERO
        Y2p => D2STATE(ppY2p,:); Y2p=ZERO
        Y3c => D2STATE(ppY3c,:); Y3c=ZERO
        Y3n => D2STATE(ppY3n,:); Y3n=ZERO
        Y3p => D2STATE(ppY3p,:); Y3p=ZERO
        Y4c => D2STATE(ppY4c,:); Y4c=ZERO
        Y4n => D2STATE(ppY4n,:); Y4n=ZERO
        Y4p => D2STATE(ppY4p,:); Y4p=ZERO
        Y5c => D2STATE(ppY5c,:); Y5c=ZERO
        Y5n => D2STATE(ppY5n,:); Y5n=ZERO
        Y5p => D2STATE(ppY5p,:); Y5p=ZERO
        Q6c => D2STATE(ppQ6c,:); Q6c=ZERO
        Q6n => D2STATE(ppQ6n,:); Q6n=ZERO
        Q6p => D2STATE(ppQ6p,:); Q6p=ZERO
        Q6s => D2STATE(ppQ6s,:); Q6s=ZERO
        Q1c => D2STATE(ppQ1c,:); Q1c=ZERO
        Q1n => D2STATE(ppQ1n,:); Q1n=ZERO
        Q1p => D2STATE(ppQ1p,:); Q1p=ZERO
        Q11c => D2STATE(ppQ11c,:); Q11c=ZERO
        Q11n => D2STATE(ppQ11n,:); Q11n=ZERO
        Q11p => D2STATE(ppQ11p,:); Q11p=ZERO
        H1c => D2STATE(ppH1c,:); H1c=ZERO
        H1n => D2STATE(ppH1n,:); H1n=ZERO
        H1p => D2STATE(ppH1p,:); H1p=ZERO
        H2c => D2STATE(ppH2c,:); H2c=ZERO
        H2n => D2STATE(ppH2n,:); H2n=ZERO
        H2p => D2STATE(ppH2p,:); H2p=ZERO
        K1p => D2STATE(ppK1p,:); K1p=ZERO
        K11p => D2STATE(ppK11p,:); K11p=ZERO
        K21p => D2STATE(ppK21p,:); K21p=ZERO
        K4n => D2STATE(ppK4n,:); K4n=ZERO
        K14n => D2STATE(ppK14n,:); K14n=ZERO
        K24n => D2STATE(ppK24n,:); K24n=ZERO
        K3n => D2STATE(ppK3n,:); K3n=ZERO
        K5s => D2STATE(ppK5s,:); K5s=ZERO
        K6r => D2STATE(ppK6r,:); K6r=ZERO
        G2o => D2STATE(ppG2o,:); G2o=ZERO
        G4n => D2STATE(ppG4n,:); G4n=ZERO
        D1m => D2STATE(ppD1m,:); D1m=ZERO
        D2m => D2STATE(ppD2m,:); D2m=ZERO
        D6m => D2STATE(ppD6m,:); D6m=ZERO
        D7m => D2STATE(ppD7m,:); D7m=ZERO
        D8m => D2STATE(ppD8m,:); D8m=ZERO
        D9m => D2STATE(ppD9m,:); D9m=ZERO


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Start the allocation of other pelagic states
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-"
  write(LOGUNIT,*) "# Allocating Other Global Variables .."

        ppsunPI(iiP1)=19
        ppsunPI(iiP2)=20
        ppsunPI(iiP3)=21
        ppsunPI(iiP4)=22
        ppqpPc(iiP1)=23
        ppqpPc(iiP2)=24
        ppqpPc(iiP3)=25
        ppqpPc(iiP4)=26
        ppqnPc(iiP1)=27
        ppqnPc(iiP2)=28
        ppqnPc(iiP3)=29
        ppqnPc(iiP4)=30
        ppqsPc(iiP1)=31
        ppqsPc(iiP2)=32
        ppqsPc(iiP3)=33
        ppqsPc(iiP4)=34
        ppqlPc(iiP1)=35
        ppqlPc(iiP2)=36
        ppqlPc(iiP3)=37
        ppqlPc(iiP4)=38
        ppqpZc(iiZ3)=39
        ppqpZc(iiZ4)=40
        ppqnZc(iiZ3)=41
        ppqnZc(iiZ4)=42
        ppqp_mz(iiZ5)=43
        ppqp_mz(iiZ6)=44
        ppqn_mz(iiZ5)=45
        ppqn_mz(iiZ6)=46
        ppsediPI(iiP1)=47
        ppsediPI(iiP2)=48
        ppsediPI(iiP3)=49
        ppsediPI(iiP4)=50
        ppeiPI(iiP1)=51
        ppeiPI(iiP2)=52
        ppeiPI(iiP3)=53
        ppeiPI(iiP4)=54
        ppEPLi(iiP1)=55
        ppEPLi(iiP2)=56
        ppEPLi(iiP3)=57
        ppEPLi(iiP4)=58

        ETW => D3DIAGNOS(ppETW,:); ETW=ZERO
        ESW => D3DIAGNOS(ppESW,:); ESW=ZERO
        EIR => D3DIAGNOS(ppEIR,:); EIR=ZERO
        ESS => D3DIAGNOS(ppESS,:); ESS=ZERO
        cxoO2 => D3DIAGNOS(ppcxoO2,:); cxoO2=ZERO
        XO2o => D3DIAGNOS(ppXO2o,:); XO2o=ZERO
        eO2mO2 => D3DIAGNOS(ppeO2mO2,:); eO2mO2=ZERO
        Chla => D3DIAGNOS(ppChla,:); Chla=ZERO
        flP1R6s => D3DIAGNOS(ppflP1R6s,:); flP1R6s=ZERO
        flPTN6r => D3DIAGNOS(ppflPTN6r,:); flPTN6r=ZERO
        flN3O4n => D3DIAGNOS(ppflN3O4n,:); flN3O4n=ZERO
        qpR6c => D3DIAGNOS(ppqpR6c,:); qpR6c=ZERO
        qnR6c => D3DIAGNOS(ppqnR6c,:); qnR6c=ZERO
        qsR6c => D3DIAGNOS(ppqsR6c,:); qsR6c=ZERO
        qpB1c => D3DIAGNOS(ppqpB1c,:); qpB1c=ZERO
        qnB1c => D3DIAGNOS(ppqnB1c,:); qnB1c=ZERO
        sediR6 => D3DIAGNOS(ppsediR6,:); sediR6=ZERO
        xEPS => D3DIAGNOS(ppxEPS,:); xEPS=ZERO

        sunPI => D3DIAGNOS(ppsunPI(iiP1): ppsunPI(iiP4),:)
        sunPI=ZERO
        qpPc => D3DIAGNOS(ppqpPc(iiP1): ppqpPc(iiP4),:)
        qpPc=ZERO
        qnPc => D3DIAGNOS(ppqnPc(iiP1): ppqnPc(iiP4),:)
        qnPc=ZERO
        qsPc => D3DIAGNOS(ppqsPc(iiP1): ppqsPc(iiP4),:)
        qsPc=ZERO
        qlPc => D3DIAGNOS(ppqlPc(iiP1): ppqlPc(iiP4),:)
        qlPc=ZERO
        qpZc => D3DIAGNOS(ppqpZc(iiZ3): ppqpZc(iiZ4),:)
        qpZc=ZERO
        qnZc => D3DIAGNOS(ppqnZc(iiZ3): ppqnZc(iiZ4),:)
        qnZc=ZERO
        qp_mz => D3DIAGNOS(ppqp_mz(iiZ5): ppqp_mz(iiZ6),:)
        qp_mz=ZERO
        qn_mz => D3DIAGNOS(ppqn_mz(iiZ5): ppqn_mz(iiZ6),:)
        qn_mz=ZERO
        sediPI => D3DIAGNOS(ppsediPI(iiP1): ppsediPI(iiP4),:)
        sediPI=ZERO
        eiPI => D3DIAGNOS(ppeiPI(iiP1): ppeiPI(iiP4),:)
        eiPI=ZERO
        EPLi => D3DIAGNOS(ppEPLi(iiP1): ppEPLi(iiP4),:)
        EPLi=ZERO


        ppruHI(iiH1)=63
        ppruHI(iiH2)=64
        ppreHI(iiH1)=65
        ppreHI(iiH2)=66
        ppretPIc(iiP1)=67
        ppretPIc(iiP2)=68
        ppretPIc(iiP3)=69
        ppretPIc(iiP4)=70
        ppretPIn(iiP1)=71
        ppretPIn(iiP2)=72
        ppretPIn(iiP3)=73
        ppretPIn(iiP4)=74
        ppretPIp(iiP1)=75
        ppretPIp(iiP2)=76
        ppretPIp(iiP3)=77
        ppretPIp(iiP4)=78
        ppretPIl(iiP1)=79
        ppretPIl(iiP2)=80
        ppretPIl(iiP3)=81
        ppretPIl(iiP4)=82
        ppretPIs(iiP1)=83
        ppretPIs(iiP2)=84
        ppretPIs(iiP3)=85
        ppretPIs(iiP4)=86

        jOAO2o => D2DIAGNOS(ppjOAO2o,:); jOAO2o=ZERO
        rutQ6c => D2DIAGNOS(pprutQ6c,:); rutQ6c=ZERO
        rutQ6n => D2DIAGNOS(pprutQ6n,:); rutQ6n=ZERO
        rutQ6p => D2DIAGNOS(pprutQ6p,:); rutQ6p=ZERO
        rutQ6s => D2DIAGNOS(pprutQ6s,:); rutQ6s=ZERO
        rutQ1c => D2DIAGNOS(pprutQ1c,:); rutQ1c=ZERO
        rutQ1n => D2DIAGNOS(pprutQ1n,:); rutQ1n=ZERO
        rutQ1p => D2DIAGNOS(pprutQ1p,:); rutQ1p=ZERO
        rrBTo => D2DIAGNOS(pprrBTo,:); rrBTo=ZERO
        rrATo => D2DIAGNOS(pprrATo,:); rrATo=ZERO
        reBTn => D2DIAGNOS(ppreBTn,:); reBTn=ZERO
        reBTp => D2DIAGNOS(ppreBTp,:); reBTp=ZERO
        reATn => D2DIAGNOS(ppreATn,:); reATn=ZERO
        reATp => D2DIAGNOS(ppreATp,:); reATp=ZERO
        turenh => D2DIAGNOS(ppturenh,:); turenh=ZERO
        irrenh => D2DIAGNOS(ppirrenh,:); irrenh=ZERO
        shiftD1m => D2DIAGNOS(ppshiftD1m,:); shiftD1m=ZERO
        shiftD2m => D2DIAGNOS(ppshiftD2m,:); shiftD2m=ZERO
        jG2K3o => D2DIAGNOS(ppjG2K3o,:); jG2K3o=ZERO
        jG2K7o => D2DIAGNOS(ppjG2K7o,:); jG2K7o=ZERO
        M1p => D2DIAGNOS(ppM1p,:); M1p=ZERO
        M11p => D2DIAGNOS(ppM11p,:); M11p=ZERO
        M21p => D2DIAGNOS(ppM21p,:); M21p=ZERO
        M4n => D2DIAGNOS(ppM4n,:); M4n=ZERO
        M14n => D2DIAGNOS(ppM14n,:); M14n=ZERO
        M24n => D2DIAGNOS(ppM24n,:); M24n=ZERO
        M3n => D2DIAGNOS(ppM3n,:); M3n=ZERO
        M5s => D2DIAGNOS(ppM5s,:); M5s=ZERO
        M6r => D2DIAGNOS(ppM6r,:); M6r=ZERO
        RIc => D2DIAGNOS(ppRIc,:); RIc=ZERO
        RIn => D2DIAGNOS(ppRIn,:); RIn=ZERO
        RIp => D2DIAGNOS(ppRIp,:); RIp=ZERO
        RIs => D2DIAGNOS(ppRIs,:); RIs=ZERO
        PIc => D2DIAGNOS(ppPIc,:); PIc=ZERO
        PIn => D2DIAGNOS(ppPIn,:); PIn=ZERO
        PIp => D2DIAGNOS(ppPIp,:); PIp=ZERO
        PIs => D2DIAGNOS(ppPIs,:); PIs=ZERO
        jPIY3c => D2DIAGNOS(ppjPIY3c,:); jPIY3c=ZERO
        jRIY3c => D2DIAGNOS(ppjRIY3c,:); jRIY3c=ZERO
        jRIY3n => D2DIAGNOS(ppjRIY3n,:); jRIY3n=ZERO
        jRIY3p => D2DIAGNOS(ppjRIY3p,:); jRIY3p=ZERO
        jRIY3s => D2DIAGNOS(ppjRIY3s,:); jRIY3s=ZERO
        Depth_Ben => D2DIAGNOS(ppDepth_Ben,:); Depth_Ben=ZERO
        ETW_Ben => D2DIAGNOS(ppETW_Ben,:); ETW_Ben=ZERO
        O2o_Ben => D2DIAGNOS(ppO2o_Ben,:); O2o_Ben=ZERO
        N1p_Ben => D2DIAGNOS(ppN1p_Ben,:); N1p_Ben=ZERO
        N3n_Ben => D2DIAGNOS(ppN3n_Ben,:); N3n_Ben=ZERO
        N4n_Ben => D2DIAGNOS(ppN4n_Ben,:); N4n_Ben=ZERO
        N5s_Ben => D2DIAGNOS(ppN5s_Ben,:); N5s_Ben=ZERO
        N6r_Ben => D2DIAGNOS(ppN6r_Ben,:); N6r_Ben=ZERO
        jG2O2o => D2DIAGNOS(ppjG2O2o,:); jG2O2o=ZERO
        jK1N1p => D2DIAGNOS(ppjK1N1p,:); jK1N1p=ZERO
        jK3N3n => D2DIAGNOS(ppjK3N3n,:); jK3N3n=ZERO
        jK4N4n => D2DIAGNOS(ppjK4N4n,:); jK4N4n=ZERO
        jK5N5s => D2DIAGNOS(ppjK5N5s,:); jK5N5s=ZERO
        jK6N6r => D2DIAGNOS(ppjK6N6r,:); jK6N6r=ZERO
        totpeln => D2DIAGNOS(pptotpeln,:); totpeln=ZERO
        totpelp => D2DIAGNOS(pptotpelp,:); totpelp=ZERO
        totpels => D2DIAGNOS(pptotpels,:); totpels=ZERO
        totbenn => D2DIAGNOS(pptotbenn,:); totbenn=ZERO
        totbenp => D2DIAGNOS(pptotbenp,:); totbenp=ZERO
        totbens => D2DIAGNOS(pptotbens,:); totbens=ZERO

        ruHI => D2DIAGNOS(ppruHI(iiH1): ppruHI(iiH2),:)
        ruHI=ZERO
        reHI => D2DIAGNOS(ppreHI(iiH1): ppreHI(iiH2),:)
        reHI=ZERO
        retPIc => D2DIAGNOS(ppretPIc(iiP1): ppretPIc(iiP4),:)
        retPIc=ZERO
        retPIn => D2DIAGNOS(ppretPIn(iiP1): ppretPIn(iiP4),:)
        retPIn=ZERO
        retPIp => D2DIAGNOS(ppretPIp(iiP1): ppretPIp(iiP4),:)
        retPIp=ZERO
        retPIl => D2DIAGNOS(ppretPIl(iiP1): ppretPIl(iiP4),:)
        retPIl=ZERO
        retPIs => D2DIAGNOS(ppretPIs(iiP1): ppretPIs(iiP4),:)
        retPIs=ZERO


        allocate(Depth(1:NO_BOXES), stat=status);  Depth = ZERO
        allocate(ABIO_eps(1:NO_BOXES), stat=status);  ABIO_eps = ZERO



        allocate(KPO4(1:NO_BOXES_XY), stat=status);  KPO4 = ZERO
        allocate(KNH4(1:NO_BOXES_XY), stat=status);  KNH4 = ZERO
        allocate(KNO3(1:NO_BOXES_XY), stat=status);  KNO3 = ZERO
        allocate(KRED(1:NO_BOXES_XY), stat=status);  KRED = ZERO
        allocate(KSIO3(1:NO_BOXES_XY), stat=status);  KSIO3 = ZERO
        allocate(KSIO3E(1:NO_BOXES_XY), stat=status);  KSIO3E = ZERO
        allocate(KQ1(1:NO_BOXES_XY), stat=status);  KQ1 = ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Start the allocation of vars for calculation of combined fluxes for output
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        allocate(flx_calc_nr(0:10),stat=status)
        allocate(flx_CalcIn(1:10),stat=status)
        allocate(flx_option(1:10),stat=status)
        allocate(flx_t(1:116),stat=status)
        allocate(flx_SS(1:116),stat=status)
        allocate(flx_states(1:116),stat=status)
        allocate(flx_ostates(1:116),stat=status)
        flx_calc_nr(0)=0
        flx_cal_ben_start=9


        ! ruPTc=P.c <- *        (normal flux):
        flx_calc_nr(1)= 4; flx_CalcIn(1)=iiPel; flx_option(1)=0
        flx_t(1)=+1.00;flx_SS(1)=0; flx_states(1)=ppP1c;flx_ostates(1)=ppP1c
        flx_t(2)=+1.00;flx_SS(2)=0; flx_states(2)=ppP2c;flx_ostates(2)=ppP2c
        flx_t(3)=+1.00;flx_SS(3)=0; flx_states(3)=ppP3c;flx_ostates(3)=ppP3c
        flx_t(4)=+1.00;flx_SS(4)=0; flx_states(4)=ppP4c;flx_ostates(4)=ppP4c

        ! ruPTn=P.n <- N3n+N4n        (normal flux):
        flx_calc_nr(2)= 12; flx_CalcIn(2)=iiPel; flx_option(2)=0
        flx_t(5)=+1.00;flx_SS(5)=0; flx_states(5)=ppP1n;flx_ostates(5)=ppN3n
        flx_t(6)=+1.00;flx_SS(6)=0; flx_states(6)=ppP1n;flx_ostates(6)=ppN4n
        flx_t(7)=+1.00;flx_SS(7)=0; flx_states(7)=ppP2n;flx_ostates(7)=ppN3n
        flx_t(8)=+1.00;flx_SS(8)=0; flx_states(8)=ppP2n;flx_ostates(8)=ppN4n
        flx_t(9)=+1.00;flx_SS(9)=0; flx_states(9)=ppP3n;flx_ostates(9)=ppN3n
        flx_t(10)=+1.00;flx_SS(10)=0; flx_states(10)=ppP3n;flx_ostates(10)=ppN4n
        flx_t(11)=+1.00;flx_SS(11)=0; flx_states(11)=ppP4n;flx_ostates(11)=ppN3n
        flx_t(12)=+1.00;flx_SS(12)=0; flx_states(12)=ppP4n;flx_ostates(12)=ppN4n

        ! ruPTP=P.p <- N1p        (normal flux):
        flx_calc_nr(3)= 16; flx_CalcIn(3)=iiPel; flx_option(3)=0
        flx_t(13)=+1.00;flx_SS(13)=0; flx_states(13)=ppP1p;flx_ostates(13)=ppN1p
        flx_t(14)=+1.00;flx_SS(14)=0; flx_states(14)=ppP2p;flx_ostates(14)=ppN1p
        flx_t(15)=+1.00;flx_SS(15)=0; flx_states(15)=ppP3p;flx_ostates(15)=ppN1p
        flx_t(16)=+1.00;flx_SS(16)=0; flx_states(16)=ppP4p;flx_ostates(16)=ppN1p

        ! netPP=(P.c<- *)-(P.c -> *)-(P.c->R1c+R2c+R6c)        (normal flux):
        flx_calc_nr(4)= 36; flx_CalcIn(4)=iiPel; flx_option(4)=0
        flx_t(17)=+1.00;flx_SS(17)=0; flx_states(17)=ppP1c;flx_ostates(17)=ppP1c
        flx_t(18)=-1.00;flx_SS(18)=1; flx_states(18)=ppP1c;flx_ostates(18)=ppP1c
        flx_t(19)=-1.00;flx_SS(19)=1; flx_states(19)=ppP1c;flx_ostates(19)=ppR1c
        flx_t(20)=-1.00;flx_SS(20)=1; flx_states(20)=ppP1c;flx_ostates(20)=ppR2c
        flx_t(21)=-1.00;flx_SS(21)=1; flx_states(21)=ppP1c;flx_ostates(21)=ppR6c
        flx_t(22)=+1.00;flx_SS(22)=0; flx_states(22)=ppP2c;flx_ostates(22)=ppP2c
        flx_t(23)=-1.00;flx_SS(23)=1; flx_states(23)=ppP2c;flx_ostates(23)=ppP2c
        flx_t(24)=-1.00;flx_SS(24)=1; flx_states(24)=ppP2c;flx_ostates(24)=ppR1c
        flx_t(25)=-1.00;flx_SS(25)=1; flx_states(25)=ppP2c;flx_ostates(25)=ppR2c
        flx_t(26)=-1.00;flx_SS(26)=1; flx_states(26)=ppP2c;flx_ostates(26)=ppR6c
        flx_t(27)=+1.00;flx_SS(27)=0; flx_states(27)=ppP3c;flx_ostates(27)=ppP3c
        flx_t(28)=-1.00;flx_SS(28)=1; flx_states(28)=ppP3c;flx_ostates(28)=ppP3c
        flx_t(29)=-1.00;flx_SS(29)=1; flx_states(29)=ppP3c;flx_ostates(29)=ppR1c
        flx_t(30)=-1.00;flx_SS(30)=1; flx_states(30)=ppP3c;flx_ostates(30)=ppR2c
        flx_t(31)=-1.00;flx_SS(31)=1; flx_states(31)=ppP3c;flx_ostates(31)=ppR6c
        flx_t(32)=+1.00;flx_SS(32)=0; flx_states(32)=ppP4c;flx_ostates(32)=ppP4c
        flx_t(33)=-1.00;flx_SS(33)=1; flx_states(33)=ppP4c;flx_ostates(33)=ppP4c
        flx_t(34)=-1.00;flx_SS(34)=1; flx_states(34)=ppP4c;flx_ostates(34)=ppR1c
        flx_t(35)=-1.00;flx_SS(35)=1; flx_states(35)=ppP4c;flx_ostates(35)=ppR2c
        flx_t(36)=-1.00;flx_SS(36)=1; flx_states(36)=ppP4c;flx_ostates(36)=ppR6c

        ! ruZTc=(Z.c<-P.c+B1c+Z.c)-(Z.c->R2c+R6c)        (normal flux):
        flx_calc_nr(5)= 76; flx_CalcIn(5)=iiPel; flx_option(5)=0
        flx_t(37)=-1.00;flx_SS(37)=1; flx_states(37)=ppZ3c;flx_ostates(37)=ppR2c
        flx_t(38)=+1.00;flx_SS(38)=0; flx_states(38)=ppZ3c;flx_ostates(38)=ppB1c
        flx_t(39)=+1.00;flx_SS(39)=0; flx_states(39)=ppZ3c;flx_ostates(39)=ppP1c
        flx_t(40)=+1.00;flx_SS(40)=0; flx_states(40)=ppZ3c;flx_ostates(40)=ppP2c
        flx_t(41)=+1.00;flx_SS(41)=0; flx_states(41)=ppZ3c;flx_ostates(41)=ppP3c
        flx_t(42)=+1.00;flx_SS(42)=0; flx_states(42)=ppZ3c;flx_ostates(42)=ppP4c
        flx_t(43)=+1.00;flx_SS(43)=0; flx_states(43)=ppZ3c;flx_ostates(43)=ppZ3c
        flx_t(44)=+1.00;flx_SS(44)=0; flx_states(44)=ppZ3c;flx_ostates(44)=ppZ4c
        flx_t(45)=+1.00;flx_SS(45)=0; flx_states(45)=ppZ3c;flx_ostates(45)=ppZ5c
        flx_t(46)=+1.00;flx_SS(46)=0; flx_states(46)=ppZ3c;flx_ostates(46)=ppZ6c
        flx_t(47)=-1.00;flx_SS(47)=1; flx_states(47)=ppZ4c;flx_ostates(47)=ppR2c
        flx_t(48)=+1.00;flx_SS(48)=0; flx_states(48)=ppZ4c;flx_ostates(48)=ppB1c
        flx_t(49)=+1.00;flx_SS(49)=0; flx_states(49)=ppZ4c;flx_ostates(49)=ppP1c
        flx_t(50)=+1.00;flx_SS(50)=0; flx_states(50)=ppZ4c;flx_ostates(50)=ppP2c
        flx_t(51)=+1.00;flx_SS(51)=0; flx_states(51)=ppZ4c;flx_ostates(51)=ppP3c
        flx_t(52)=+1.00;flx_SS(52)=0; flx_states(52)=ppZ4c;flx_ostates(52)=ppP4c
        flx_t(53)=+1.00;flx_SS(53)=0; flx_states(53)=ppZ4c;flx_ostates(53)=ppZ3c
        flx_t(54)=+1.00;flx_SS(54)=0; flx_states(54)=ppZ4c;flx_ostates(54)=ppZ4c
        flx_t(55)=+1.00;flx_SS(55)=0; flx_states(55)=ppZ4c;flx_ostates(55)=ppZ5c
        flx_t(56)=+1.00;flx_SS(56)=0; flx_states(56)=ppZ4c;flx_ostates(56)=ppZ6c
        flx_t(57)=-1.00;flx_SS(57)=1; flx_states(57)=ppZ5c;flx_ostates(57)=ppR2c
        flx_t(58)=+1.00;flx_SS(58)=0; flx_states(58)=ppZ5c;flx_ostates(58)=ppB1c
        flx_t(59)=+1.00;flx_SS(59)=0; flx_states(59)=ppZ5c;flx_ostates(59)=ppP1c
        flx_t(60)=+1.00;flx_SS(60)=0; flx_states(60)=ppZ5c;flx_ostates(60)=ppP2c
        flx_t(61)=+1.00;flx_SS(61)=0; flx_states(61)=ppZ5c;flx_ostates(61)=ppP3c
        flx_t(62)=+1.00;flx_SS(62)=0; flx_states(62)=ppZ5c;flx_ostates(62)=ppP4c
        flx_t(63)=+1.00;flx_SS(63)=0; flx_states(63)=ppZ5c;flx_ostates(63)=ppZ3c
        flx_t(64)=+1.00;flx_SS(64)=0; flx_states(64)=ppZ5c;flx_ostates(64)=ppZ4c
        flx_t(65)=+1.00;flx_SS(65)=0; flx_states(65)=ppZ5c;flx_ostates(65)=ppZ5c
        flx_t(66)=+1.00;flx_SS(66)=0; flx_states(66)=ppZ5c;flx_ostates(66)=ppZ6c
        flx_t(67)=-1.00;flx_SS(67)=1; flx_states(67)=ppZ6c;flx_ostates(67)=ppR2c
        flx_t(68)=+1.00;flx_SS(68)=0; flx_states(68)=ppZ6c;flx_ostates(68)=ppB1c
        flx_t(69)=+1.00;flx_SS(69)=0; flx_states(69)=ppZ6c;flx_ostates(69)=ppP1c
        flx_t(70)=+1.00;flx_SS(70)=0; flx_states(70)=ppZ6c;flx_ostates(70)=ppP2c
        flx_t(71)=+1.00;flx_SS(71)=0; flx_states(71)=ppZ6c;flx_ostates(71)=ppP3c
        flx_t(72)=+1.00;flx_SS(72)=0; flx_states(72)=ppZ6c;flx_ostates(72)=ppP4c
        flx_t(73)=+1.00;flx_SS(73)=0; flx_states(73)=ppZ6c;flx_ostates(73)=ppZ3c
        flx_t(74)=+1.00;flx_SS(74)=0; flx_states(74)=ppZ6c;flx_ostates(74)=ppZ4c
        flx_t(75)=+1.00;flx_SS(75)=0; flx_states(75)=ppZ6c;flx_ostates(75)=ppZ5c
        flx_t(76)=+1.00;flx_SS(76)=0; flx_states(76)=ppZ6c;flx_ostates(76)=ppZ6c

        ! rrPTo=(O2o->*)        (normal flux):
        flx_calc_nr(6)= 77; flx_CalcIn(6)=iiPel; flx_option(6)=0
        flx_t(77)=+1.00;flx_SS(77)=1; flx_states(77)=ppO2o;flx_ostates(77)=ppO2o

        ! rePTn=(N3n+N4n<-B1n+Z.n)        (normal flux):
        flx_calc_nr(7)= 87; flx_CalcIn(7)=iiPel; flx_option(7)=0
        flx_t(78)=+1.00;flx_SS(78)=0; flx_states(78)=ppN3n;flx_ostates(78)=ppB1n
        flx_t(79)=+1.00;flx_SS(79)=0; flx_states(79)=ppN3n;flx_ostates(79)=ppZ3n
        flx_t(80)=+1.00;flx_SS(80)=0; flx_states(80)=ppN3n;flx_ostates(80)=ppZ4n
        flx_t(81)=+1.00;flx_SS(81)=0; flx_states(81)=ppN3n;flx_ostates(81)=ppZ5n
        flx_t(82)=+1.00;flx_SS(82)=0; flx_states(82)=ppN3n;flx_ostates(82)=ppZ6n
        flx_t(83)=+1.00;flx_SS(83)=0; flx_states(83)=ppN4n;flx_ostates(83)=ppB1n
        flx_t(84)=+1.00;flx_SS(84)=0; flx_states(84)=ppN4n;flx_ostates(84)=ppZ3n
        flx_t(85)=+1.00;flx_SS(85)=0; flx_states(85)=ppN4n;flx_ostates(85)=ppZ4n
        flx_t(86)=+1.00;flx_SS(86)=0; flx_states(86)=ppN4n;flx_ostates(86)=ppZ5n
        flx_t(87)=+1.00;flx_SS(87)=0; flx_states(87)=ppN4n;flx_ostates(87)=ppZ6n

        ! rePTp=(N1p<-B1p+Z.p)        (normal flux):
        flx_calc_nr(8)= 92; flx_CalcIn(8)=iiPel; flx_option(8)=0
        flx_t(88)=+1.00;flx_SS(88)=0; flx_states(88)=ppN1p;flx_ostates(88)=ppB1p
        flx_t(89)=+1.00;flx_SS(89)=0; flx_states(89)=ppN1p;flx_ostates(89)=ppZ3p
        flx_t(90)=+1.00;flx_SS(90)=0; flx_states(90)=ppN1p;flx_ostates(90)=ppZ4p
        flx_t(91)=+1.00;flx_SS(91)=0; flx_states(91)=ppN1p;flx_ostates(91)=ppZ5p
        flx_t(92)=+1.00;flx_SS(92)=0; flx_states(92)=ppN1p;flx_ostates(92)=ppZ6p

        ! suPTc=P.c <- *        (specific flux):
        flx_calc_nr(9)= 96; flx_CalcIn(9)=iiPel; flx_option(9)=1
        flx_t(93)=+1.00;flx_SS(93)=0; flx_states(93)=ppP1c;flx_ostates(93)=ppP1c
        flx_t(94)=+1.00;flx_SS(94)=0; flx_states(94)=ppP2c;flx_ostates(94)=ppP2c
        flx_t(95)=+1.00;flx_SS(95)=0; flx_states(95)=ppP3c;flx_ostates(95)=ppP3c
        flx_t(96)=+1.00;flx_SS(96)=0; flx_states(96)=ppP4c;flx_ostates(96)=ppP4c


        ! netPPm2=(P.c<- *)-(P.c -> *)-(P.c->R1c+R2c+R6c)        (perm2 flux):
        flx_calc_nr(10)= 116; flx_CalcIn(10)=iiPel; flx_option(10)=2
        flx_t(97)=+1.00;flx_SS(97)=0; flx_states(97)=ppP1c;flx_ostates(97)=ppP1c
        flx_t(98)=-1.00;flx_SS(98)=1; flx_states(98)=ppP1c;flx_ostates(98)=ppP1c
        flx_t(99)=-1.00;flx_SS(99)=1; flx_states(99)=ppP1c;flx_ostates(99)=ppR1c
        flx_t(100)=-1.00;flx_SS(100)=1; flx_states(100)=ppP1c
        flx_ostates(100)=ppR2c
        flx_t(101)=-1.00;flx_SS(101)=1; flx_states(101)=ppP1c
        flx_ostates(101)=ppR6c
        flx_t(102)=+1.00;flx_SS(102)=0; flx_states(102)=ppP2c
        flx_ostates(102)=ppP2c
        flx_t(103)=-1.00;flx_SS(103)=1; flx_states(103)=ppP2c
        flx_ostates(103)=ppP2c
        flx_t(104)=-1.00;flx_SS(104)=1; flx_states(104)=ppP2c
        flx_ostates(104)=ppR1c
        flx_t(105)=-1.00;flx_SS(105)=1; flx_states(105)=ppP2c
        flx_ostates(105)=ppR2c
        flx_t(106)=-1.00;flx_SS(106)=1; flx_states(106)=ppP2c
        flx_ostates(106)=ppR6c
        flx_t(107)=+1.00;flx_SS(107)=0; flx_states(107)=ppP3c
        flx_ostates(107)=ppP3c
        flx_t(108)=-1.00;flx_SS(108)=1; flx_states(108)=ppP3c
        flx_ostates(108)=ppP3c
        flx_t(109)=-1.00;flx_SS(109)=1; flx_states(109)=ppP3c
        flx_ostates(109)=ppR1c
        flx_t(110)=-1.00;flx_SS(110)=1; flx_states(110)=ppP3c
        flx_ostates(110)=ppR2c
        flx_t(111)=-1.00;flx_SS(111)=1; flx_states(111)=ppP3c
        flx_ostates(111)=ppR6c
        flx_t(112)=+1.00;flx_SS(112)=0; flx_states(112)=ppP4c
        flx_ostates(112)=ppP4c
        flx_t(113)=-1.00;flx_SS(113)=1; flx_states(113)=ppP4c
        flx_ostates(113)=ppP4c
        flx_t(114)=-1.00;flx_SS(114)=1; flx_states(114)=ppP4c
        flx_ostates(114)=ppR1c
        flx_t(115)=-1.00;flx_SS(115)=1; flx_states(115)=ppP4c
        flx_ostates(115)=ppR2c
        flx_t(116)=-1.00;flx_SS(116)=1; flx_states(116)=ppP4c
        flx_ostates(116)=ppR6c


  end subroutine

