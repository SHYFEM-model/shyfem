!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.50-g
!
! MODULE
!   ModuleMem
!
! FILE
!   ModuleMem
!
! DESCRIPTION
!   Definition of Global Shared Memory
!  
!   This module contains all the structural definitions of the BFM
!   and sets up the memory layout.
!   It is automatically generated from the prototype file 
!   BFM/proto/ModuleMem.proto by including the information from 
!   BFM/General/GlobalDefsBFM.model
!   Do not directly edit this code because changes will be lost at
!   any new compilation.
!
! AUTHORS
!   Piet Ruardij and Marcello Vichi
!
! CHANGE_LOG
!   ---
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the BFM team
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
!
!! IMPORTANT NOTE:
!! Do not change the lines starting with two comment characters "!" 
!! These lines are used by the parser to generate the final module file

!

#include"cppdefs.h"
#include "DEBUG.h"

      module mem
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Modules can optionally use (import) other modules
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        use global_mem, only:RLEN, ZERO
#if DEBUG && BFM_GOTM
        use gotm_error_msg, only: gotm_error
#endif
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Default all is private
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        private
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! State variables Info
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        !    3d name                                                  description            unit
! ---------- ------------------------------------------------------------ ---------------
!        O2o                                                        Oxgen      mmol O2/m3
!        N1p                                                    Phosphate       mmol P/m3
!        N3n                                                      Nitrate       mmol N/m3
!        N4n                                                     Ammonium       mmol N/m3
!        O4n                                                 NitrogenSink       mmol N/m3
!        N5s                                                     Silicate      mmol Si/m3
!        N6r                                        Reduction Equivalents     mmol S--/m3
!        B1c                                             Pelagic Bacteria         mg C/m3
!        B1n                                             Pelagic Bacteria       mmol N/m3
!        B1p                                             Pelagic Bacteria       mmol P/m3
!        P1c                                                      Diatoms         mg C/m3
!        P1n                                                      Diatoms       mmol N/m3
!        P1p                                                      Diatoms       mmol P/m3
!        P1l                                                      Diatoms       mg Chl/m3
!        P1s                                                      Diatoms     mmmol Si/m3
!        P2c                                                  Flagellates         mg C/m3
!        P2n                                                  Flagellates       mmol N/m3
!        P2p                                                  Flagellates       mmol P/m3
!        P2l                                                  Flagellates       mg Chl/m3
!        P3c                                            PicoPhytoPlankton         mg C/m3
!        P3n                                            PicoPhytoPlankton       mmol N/m3
!        P3p                                            PicoPhytoPlankton       mmol P/m3
!        P3l                                            PicoPhytoPlankton       mg Chl/m3
!        P4c                                              Dinoflagellates         mg C/m3
!        P4n                                              Dinoflagellates       mmol N/m3
!        P4p                                              Dinoflagellates       mmol P/m3
!        P4l                                              Dinoflagellates       mg Chl/m3
!        Z3c                                  Carnivorous mesozooplankton         mg C/m3
!        Z3n                                  Carnivorous mesozooplankton       mmol N/m3
!        Z3p                                  Carnivorous mesozooplankton       mmol P/m3
!        Z4c                                   Omnivorous mesozooplankton         mg C/m3
!        Z4n                                   Omnivorous mesozooplankton       mmol N/m3
!        Z4p                                   Omnivorous mesozooplankton       mmol P/m3
!        Z5c                                             Microzooplankton         mg C/m3
!        Z5n                                             Microzooplankton       mmol N/m3
!        Z5p                                             Microzooplankton       mmol P/m3
!        Z6c                         Heterotrophic nanoflagellates (HNAN)         mg C/m3
!        Z6n                         Heterotrophic nanoflagellates (HNAN)       mmol N/m3
!        Z6p                         Heterotrophic nanoflagellates (HNAN)       mmol P/m3
!        R1c                                  Labile Organic Carbon (LOC)         mg C/m3
!        R1n                                  Labile Organic Carbon (LOC)       mmol N/m3
!        R1p                                  Labile Organic Carbon (LOC)       mmol P/m3
!        R2c                                       CarboHydrates (sugars)         mg C/m3
!        R6c                             Particulate Organic Carbon (POC)         mg C/m3
!        R6n                             Particulate Organic Carbon (POC)       mmol N/m3
!        R6p                             Particulate Organic Carbon (POC)       mmol P/m3
!        R6s                             Particulate Organic Carbon (POC)     mmmol Si/m3
!        R7c                                               Refractory DOC         mg C/m3


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        !    2d name                                                  description            unit
! ---------- ------------------------------------------------------------ ---------------
!        Y1c                                                   Epibenthos         mg C/m2
!        Y1n                                                   Epibenthos       mmol N/m2
!        Y1p                                                   Epibenthos       mmol P/m2
!        Y2c                                              Deposit feeders         mg C/m2
!        Y2n                                              Deposit feeders       mmol N/m2
!        Y2p                                              Deposit feeders       mmol P/m2
!        Y3c                                           Suspension feeders         mg C/m2
!        Y3n                                           Suspension feeders       mmol N/m2
!        Y3p                                           Suspension feeders       mmol P/m2
!        Y4c                                                  Meiobenthos         mg C/m2
!        Y4n                                                  Meiobenthos       mmol N/m2
!        Y4p                                                  Meiobenthos       mmol P/m2
!        Y5c                                            Benthic predators         mg C/m2
!        Y5n                                            Benthic predators       mmol N/m2
!        Y5p                                            Benthic predators       mmol P/m2
!        Q6c                                   Particulate organic carbon         mg C/m2
!        Q6n                                   Particulate organic carbon       mmol N/m2
!        Q6p                                   Particulate organic carbon       mmol P/m2
!        Q6s                                   Particulate organic carbon       mmolSi/m2
!        Q1c                                        Labile organic carbon         mg C/m2
!        Q1n                                        Labile organic carbon       mmol N/m2
!        Q1p                                        Labile organic carbon       mmol P/m2
!       Q11c                                        Labile organic carbon         mg C/m2
!       Q11n                                        Labile organic carbon       mmol N/m2
!       Q11p                                        Labile organic carbon       mmol P/m2
!        H1c                                     Aerobic benthic bacteria         mg C/m2
!        H1n                                     Aerobic benthic bacteria       mmol N/m2
!        H1p                                     Aerobic benthic bacteria       mmol P/m2
!        H2c                                   Anaerobic benthic bacteria         mg C/m2
!        H2n                                   Anaerobic benthic bacteria       mmol N/m2
!        H2p                                   Anaerobic benthic bacteria       mmol P/m2
!        K1p                                      Phosphate in oxic layer       mmol P/m3
!       K11p                           Phosphate in denitrification layer       mmol P/m3
!       K21p                                    Phosphate in anoxic layer       mmol P/m3
!        K4n                                       Ammonium in oxic layer       mmol N/m3
!       K14n                            Ammonium in denitrification layer       mmol N/m3
!       K24n                                     Ammonium in anoxic layer       mmol N/m3
!        K3n                                         Nitrate in sediments       mmol N/m2
!        K5s                                        Silicate in sediments      mmol Si/m2
!        K6r                           Reduction equivalents in sediments      mmolS--/m2
!        G2o                                                   Benthic O2      mmol O2/m3
!        G4n                                  N2 sink for benthic system.       mmol N/m3
!        D1m                                     Oxygen penetration depth               m
!        D2m                                        Denitrification depth               m
!        D6m                                  Penetration Depth organic C               m
!        D7m                                  Penetration Depth organic N               m
!        D8m                                  Penetration Depth organic P               m
!        D9m                                 Penetration Depth organic Si               m


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! GLOBAL system CONSTANTS
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public  ::NO_D3_BOX_STATES= 48

        integer,parameter,public  ::NO_D2_BOX_STATES= 48

        integer,parameter,public  ::NO_D3_BOX_DIAGNOSS= 58

        integer,parameter,public  ::NO_D2_BOX_DIAGNOSS= 86

        integer,parameter,public  ::NO_D3_BOX_FLUX= 9

        integer,parameter,public  ::NO_D2_BOX_FLUX= 1

        integer,public  ::NO_BOXES
        integer,public  ::NO_BOXES_X
        integer,public  ::NO_BOXES_Y
        integer,public  ::NO_BOXES_Z
        integer,public  ::NO_STATES
        integer,public  ::NO_BOXES_XY
        integer,parameter,public  ::iiPel= 0
        integer,parameter,public  ::iiBen= 1000
        integer,parameter,public  ::iiRESET= -1000
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Specification of State variables
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        real(RLEN),public,pointer,dimension(:,:)  :: D3STATE
        real(RLEN),public,pointer,dimension(:,:,:)  :: D3SOURCE
        real(RLEN),public,pointer,dimension(:,:,:)  :: D3SINK
        integer,public,pointer,dimension(:)  :: D3STATETYPE

        real(RLEN),public,pointer,dimension(:,:)  :: D2STATE
        real(RLEN),public,pointer,dimension(:,:,:)  :: D2SOURCE
        real(RLEN),public,pointer,dimension(:,:,:)  :: D2SINK
        integer,public,pointer,dimension(:)  :: D2STATETYPE

        real(RLEN),public,pointer,dimension(:,:)  :: D3DIAGNOS

        real(RLEN),public,pointer,dimension(:,:)  :: D2DIAGNOS


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !! GLOBAL definition of Pelagic (D3) state variables
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        integer,parameter,public :: ppO2o=1, ppN1p=2, ppN3n=3, ppN4n=4, ppO4n=5,&
          ppN5s=6, ppN6r=7, ppB1c=8, ppB1n=9, ppB1p=10, ppP1c=11, ppP1n=12,&
          ppP1p=13, ppP1l=14, ppP1s=15, ppP2c=16, ppP2n=17, ppP2p=18, ppP2l=19,&
          ppP3c=20, ppP3n=21, ppP3p=22, ppP3l=23, ppP4c=24, ppP4n=25, ppP4p=26,&
          ppP4l=27, ppZ3c=28, ppZ3n=29, ppZ3p=30, ppZ4c=31, ppZ4n=32, ppZ4p=33,&
          ppZ5c=34, ppZ5n=35, ppZ5p=36, ppZ6c=37, ppZ6n=38, ppZ6p=39, ppR1c=40,&
          ppR1n=41, ppR1p=42, ppR2c=43, ppR6c=44, ppR6n=45, ppR6p=46, ppR6s=47,&
          ppR7c=48, ppP2s=0, ppP3s=0, ppP4s=0


        real(RLEN),public,dimension(:),pointer :: O2o, N1p, N3n, N4n, O4n, N5s,&
          N6r, B1c, B1n, B1p, P1c, P1n, P1p, P1l, P1s, P2c, P2n, P2p, P2l, P3c,&
          P3n, P3p, P3l, P4c, P4n, P4p, P4l, Z3c, Z3n, Z3p, Z4c, Z4n, Z4p, Z5c,&
          Z5n, Z5p, Z6c, Z6n, Z6p, R1c, R1n, R1p, R2c, R6c, R6n, R6p, R6s, R7c


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !! GLOBAL definition of Benthic (D2) state variables
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public :: ppY1c=1, ppY1n=2, ppY1p=3, ppY2c=4, ppY2n=5,&
          ppY2p=6, ppY3c=7, ppY3n=8, ppY3p=9, ppY4c=10, ppY4n=11, ppY4p=12,&
          ppY5c=13, ppY5n=14, ppY5p=15, ppQ6c=16, ppQ6n=17, ppQ6p=18, ppQ6s=19,&
          ppQ1c=20, ppQ1n=21, ppQ1p=22, ppQ11c=23, ppQ11n=24, ppQ11p=25, ppH1c=26,&
          ppH1n=27, ppH1p=28, ppH2c=29, ppH2n=30, ppH2p=31, ppK1p=32, ppK11p=33,&
          ppK21p=34, ppK4n=35, ppK14n=36, ppK24n=37, ppK3n=38, ppK5s=39, ppK6r=40,&
          ppG2o=41, ppG4n=42, ppD1m=43, ppD2m=44, ppD6m=45, ppD7m=46, ppD8m=47,&
          ppD9m=48


        real(RLEN),public,dimension(:),pointer :: Y1c, Y1n, Y1p, Y2c, Y2n, Y2p,&
          Y3c, Y3n, Y3p, Y4c, Y4n, Y4p, Y5c, Y5n, Y5p, Q6c, Q6n, Q6p, Q6s, Q1c,&
          Q1n, Q1p, Q11c, Q11n, Q11p, H1c, H1n, H1p, H2c, H2n, H2p, K1p, K11p,&
          K21p, K4n, K14n, K24n, K3n, K5s, K6r, G2o, G4n, D1m, D2m, D6m, D7m, D8m,&
          D9m


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Constituent parameters:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public       :: iiC=1, iiN=2, iiP=3, iiL=4, iiS=5


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Group parameters:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public :: iiPhytoPlankton=4, iiP1=1, iiP2=2, iiP3=3,&
          iiP4=4
        integer,parameter,public       :: iiMesoZooPlankton=2, iiZ3=1, iiZ4=2
        integer,parameter,public       :: iiMicroZooPlankton=2, iiZ5=1, iiZ6=2


        integer,parameter,public :: iiBenOrganisms=5, iiY1=1, iiY2=2, iiY3=3,&
          iiY4=4, iiY5=5
        integer,parameter,public       :: iiBenDetritus=2, iiQ1=1, iiQ11=2
        integer,parameter,public       :: iiBenBacteria=2, iiH1=1, iiH2=2
        integer,parameter,public :: iiBenthicPhosphate=3, iiK1=1, iiK11=2,&
          iiK21=3
        integer,parameter,public :: iiBenthicAmmonium=3, iiK4=1, iiK14=2,&
          iiK24=3


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  Global Variables
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        integer,public  :: BoxNumber
        integer,public  :: BoxNumberX
        integer,public  :: BoxNumberY
        integer,public  :: BoxNumberZ
        integer,public  :: BoxNumberXY

        real(RLEN),public                                      :: &
          dummy,  &   !
          LocalDelta,  &   !
          Wind,  &   !
          SUNQ,  &   !Daylength in hours (h)
          ThereIsLight      !Forcing for day/night
 
        integer,public                                      :: &
          idummy      !

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        !! GLOBAL definition of Pelagic (D3) variables which can be outputted in netcdf
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-

!    3d name                                                  description            unit
! ---------- ------------------------------------------------------------ ---------------
!        ETW                                                  temperature               C
!        ESW                                                     Salinity              --
!        EIR                                                   Irradiance            W/m2
!        ESS                                           Suspended sediment            g/m3
!      cxoO2                      Temperature-dependent oxygen saturation      mmol O2/m3
!       XO2o                                      Net O2 conc. = O2 - H2S      mmol O2/m3
!     eO2mO2                                   Relative Oxygen saturation               %
!       Chla                                                   Chlorophyl       mg Chl/m3
!    flP1R6s                                            part. silica flux     mmol Si /m3
!    flPTN6r                                     anaerobic mineralixation     mmol S--/m3
!    flN3O4n                                  denitrification flux (sink)     mmol N/m3/d
!      qpR6c                                             Quotum p/c in R6     mmol N/mg C
!      qnR6c                                             Quotum n/c in R6     mmol P/mg C
!      qsR6c                                            Quotum si/c in R6    mmol Si/mg C
!      qpB1c                                             Quotum p/c in B1     mmol N/mg C
!      qnB1c                                             Quotum n/c in B1     mmol P/mg C
!     sediR6                                  Detritus sedimentation rate             m/d
!       xEPS                                 Total extinction coefficient             1/m

! sunPI(iiP1)                         Specific Net Prod. P1(PhytoPlankton)     mgC/mgC/day
! sunPI(iiP2)                         Specific Net Prod. P2(PhytoPlankton)     mgC/mgC/day
! sunPI(iiP3)                         Specific Net Prod. P3(PhytoPlankton)     mgC/mgC/day
! sunPI(iiP4)                         Specific Net Prod. P4(PhytoPlankton)     mgC/mgC/day
! qpPc(iiP1)                              Quotum p/c in P1(PhytoPlankton)     mmol P/mg C
! qpPc(iiP2)                              Quotum p/c in P2(PhytoPlankton)     mmol P/mg C
! qpPc(iiP3)                              Quotum p/c in P3(PhytoPlankton)     mmol P/mg C
! qpPc(iiP4)                              Quotum p/c in P4(PhytoPlankton)     mmol P/mg C
! qnPc(iiP1)                              Quotum n/c in P1(PhytoPlankton)     mmol N/mg C
! qnPc(iiP2)                              Quotum n/c in P2(PhytoPlankton)     mmol N/mg C
! qnPc(iiP3)                              Quotum n/c in P3(PhytoPlankton)     mmol N/mg C
! qnPc(iiP4)                              Quotum n/c in P4(PhytoPlankton)     mmol N/mg C
! qsPc(iiP1)                             Quotum si/c in P1(PhytoPlankton)    mmol Si/mg C
! qsPc(iiP2)                             Quotum si/c in P2(PhytoPlankton)    mmol Si/mg C
! qsPc(iiP3)                             Quotum si/c in P3(PhytoPlankton)    mmol Si/mg C
! qsPc(iiP4)                             Quotum si/c in P4(PhytoPlankton)    mmol Si/mg C
! qlPc(iiP1)                            Quotum chl/c in P1(PhytoPlankton)    mg Chl /mg C
! qlPc(iiP2)                            Quotum chl/c in P2(PhytoPlankton)    mg Chl /mg C
! qlPc(iiP3)                            Quotum chl/c in P3(PhytoPlankton)    mg Chl /mg C
! qlPc(iiP4)                            Quotum chl/c in P4(PhytoPlankton)    mg Chl /mg C
! qpZc(iiZ3)                               Quotum p/c Z3(MesoZooPlankton)     mmol P/mg C
! qpZc(iiZ4)                               Quotum p/c Z4(MesoZooPlankton)     mmol P/mg C
! qnZc(iiZ3)                               Quotum n/c Z3(MesoZooPlankton)     mmol N/mg C
! qnZc(iiZ4)                               Quotum n/c Z4(MesoZooPlankton)     mmol N/mg C
! qp_mz(iiZ5)                              Quotum p/c Z5(MicroZooPlankton)     mmol P/mg C
! qp_mz(iiZ6)                              Quotum p/c Z6(MicroZooPlankton)     mmol P/mg C
! qn_mz(iiZ5)                              Quotum n/c Z5(MicroZooPlankton)     mmol N/mg C
! qn_mz(iiZ6)                              Quotum n/c Z6(MicroZooPlankton)     mmol N/mg C
! sediPI(iiP1)                         P1(PhytoPlankton) sedimentation rate             m/d
! sediPI(iiP2)                         P2(PhytoPlankton) sedimentation rate             m/d
! sediPI(iiP3)                         P3(PhytoPlankton) sedimentation rate             m/d
! sediPI(iiP4)                         P4(PhytoPlankton) sedimentation rate             m/d
! eiPI(iiP1)             Regulating factor for light in P1(PhytoPlankton)               -
! eiPI(iiP2)             Regulating factor for light in P2(PhytoPlankton)               -
! eiPI(iiP3)             Regulating factor for light in P3(PhytoPlankton)               -
! eiPI(iiP4)             Regulating factor for light in P4(PhytoPlankton)               -
! EPLi(iiP1)                           Optimal light in P1(PhytoPlankton)            W/m2
! EPLi(iiP2)                           Optimal light in P2(PhytoPlankton)            W/m2
! EPLi(iiP3)                           Optimal light in P3(PhytoPlankton)            W/m2
! EPLi(iiP4)                           Optimal light in P4(PhytoPlankton)            W/m2

        integer,parameter,public :: ppETW=1, ppESW=2, ppEIR=3, ppESS=4,&
          ppcxoO2=5, ppXO2o=6, ppeO2mO2=7, ppChla=8, ppflP1R6s=9, ppflPTN6r=10,&
          ppflN3O4n=11, ppqpR6c=12, ppqnR6c=13, ppqsR6c=14, ppqpB1c=15,&
          ppqnB1c=16, ppsediR6=17, ppxEPS=18

        integer,public ::&
          ppsunPI(iiPhytoPlankton), ppqpPc(iiPhytoPlankton),&
          ppqnPc(iiPhytoPlankton), ppqsPc(iiPhytoPlankton),&
          ppqlPc(iiPhytoPlankton), ppqpZc(iiMesoZooPlankton),&
          ppqnZc(iiMesoZooPlankton), ppqp_mz(iiMicroZooPlankton),&
          ppqn_mz(iiMicroZooPlankton), ppsediPI(iiPhytoPlankton),&
          ppeiPI(iiPhytoPlankton), ppEPLi(iiPhytoPlankton)

        real(RLEN),public,dimension(:),pointer :: ETW, ESW, EIR, ESS, cxoO2,&
          XO2o, eO2mO2, Chla, flP1R6s, flPTN6r, flN3O4n, qpR6c, qnR6c, qsR6c,&
          qpB1c, qnB1c, sediR6, xEPS

        real(RLEN),public,dimension(:,:),pointer :: sunPI, qpPc, qnPc, qsPc,&
          qlPc, qpZc, qnZc, qp_mz, qn_mz, sediPI, eiPI, EPLi


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !!  GLOBAL definition of Benthic (D2) variables which can be outputted in netcdf
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!    2d name                                                  description            unit
! ---------- ------------------------------------------------------------ ---------------
!     jOAO2o                                       oxygen reaeration flux  mmol O2/m2/day
!     rutQ6c                                 Total input of R6 to benthos       mg C/m2/d
!     rutQ6n                                 Total input of R6 to benthos     mmol N/m2/d
!     rutQ6p                                 Total input of R6 to benthos     mmol P/m2/d
!     rutQ6s                                 Total input of R6 to benthos    mmol Si/m2/d
!     rutQ1c                                 Total input of R1 to benthos       mg C/m2/d
!     rutQ1n                                 Total input of R1 to benthos     mmol N/m2/d
!     rutQ1p                                 Total input of R1 to benthos     mmol P/m2/d
!      rrBTo                                Total Benthic oxic respiratio    mmol O2/m2/d
!      rrATo                             Total Benthic anoxic respiration    mmol O2/m2/d
!      reBTn                                Total Benthic oxic respiratio     mmol N/m2/d
!      reBTp                                Total Benthic oxic respiratio     mmol P/m2/d
!      reATn                             Total Benthic anoxic respiration     mmol N/m2/d
!      reATp                             Total Benthic anoxic respiration     mmol P/m2/d
!     turenh                      Enhancement factor due to bioirrigation               -
!     irrenh                       Enhancement factor due to bioturbation               -
!   shiftD1m                                        Rate of change of D1m             m/d
!   shiftD2m                                        Rate of change of D2m             m/d
!     jG2K3o                         Oxygen consumption by nitrification.    mmol O2/m2/d
!     jG2K7o                     ReOxidation of Red.Equiv. in  oxic layer   mmol S--/m2/d
!        M1p                                      phosphate in oxic layer       mmol P/m3
!       M11p                           phosphate in denitrification layer       mmol P/m3
!       M21p                                    phosphate in anoxic layer       mmol P/m3
!        M4n                                       ammonium in oxic layer       mmol N/m3
!       M14n                            ammonium in denitrification layer       mmol N/m3
!       M24n                                     ammonium in anoxic layer       mmol N/m3
!        M3n                                         nitrate in porewater       mmol N/m3
!        M5s                       silicate in oxic_denitrification layer      mmol Si/m3
!        M6r                            reduction equivalent in porewater     mmol S--/m3
!        RIc                               Detritus Food for Filterfeeder         mg C/m3
!        RIn                               Detritus Food for Filterfeeder      mmol N /m3
!        RIp                               Detritus Food for Filterfeeder       mmol N/m3
!        RIs                               Detritus Food for Filterfeeder      mmol Si/m3
!        PIc                                     Total phytoplankton Food         mg C/m3
!        PIn                                     Total phytoplankton Food      mmol N /m3
!        PIp                                     Total phytoplankton Food       mmol N/m3
!        PIs                                     Total phytoplankton Food      mmol Si/m3
!     jPIY3c                     phytoplankton filtered by filter feeders       mg C/m2/d
!     jRIY3c                           detrtus filtered by filter feeders       mg C/m2/d
!     jRIY3n                           detrtus filtered by filter feeders     mmol N/m2/d
!     jRIY3p                           detrtus filtered by filter feeders     mmol P/m2/d
!     jRIY3s                           detrtus filtered by filter feeders    mmol Si/m2/d
!  Depth_Ben                           depth of the  layer above sediment               m
!    ETW_Ben                                                  temperature               C
!    O2o_Ben                                 oxygen conc. in the pelagica      mmol O2/m3
!    N1p_Ben                               phosphate conc. in the pelagic       mmol P/m3
!    N3n_Ben                                 nitrate conc. in the pelagic       mmol N/m3
!    N4n_Ben                                ammonium conc. in the pelagic       mmol N/m3
!    N5s_Ben                                silicate conc. in the pelagic      mmol Si/m3
!    N6r_Ben                             red. equiv. conc. in the pelagic     mmol S--/m3
!     jG2O2o                      flux of oxygen from sediment to pelagic       mmol/m2/d
!     jK1N1p                   flux of phosphate from sediment to pelagic       mmol/m2/d
!     jK3N3n                     flux of nitrate from sediment to pelagic       mmol/m2/d
!     jK4N4n                    flux of ammonium from sediment to pelagic       mmol/m2/d
!     jK5N5s                    flux of silicate from sediment to pelagic       mmol/m2/d
!     jK6N6r                 flux of red. equiv. from sediment to pelagic       mmol/m2/d
!    totpeln                                total mass present in pelagic       mmol N/m3
!    totpelp                                total mass present in pelagic       mmol P/m3
!    totpels                                total mass present in pelagic      mmol Si/m3
!    totbenn                                total mass present in benthos       mmol N/m3
!    totbenp                                total mass present in benthos       mmol P/m3
!    totbens                                total mass present in benthos      mmol Si/m3

! ruHI(iiH1)                          uptake of Q1/Q11 by H1(BenBacteria)       mg C/m2/d
! ruHI(iiH2)                          uptake of Q1/Q11 by H2(BenBacteria)       mg C/m2/d
! reHI(iiH1)                          excretion of Q1/Q11 H1(BenBacteria)       mg C/m2/d
! reHI(iiH2)                          excretion of Q1/Q11 H2(BenBacteria)       mg C/m2/d
! retPIc(iiP1)                  P1(PhytoPlankton) uptaken by filter feeders       mg C/m2/d
! retPIc(iiP2)                  P2(PhytoPlankton) uptaken by filter feeders       mg C/m2/d
! retPIc(iiP3)                  P3(PhytoPlankton) uptaken by filter feeders       mg C/m2/d
! retPIc(iiP4)                  P4(PhytoPlankton) uptaken by filter feeders       mg C/m2/d
! retPIn(iiP1)                  P1(PhytoPlankton) uptaken by filter feeders     mmol N/m2/d
! retPIn(iiP2)                  P2(PhytoPlankton) uptaken by filter feeders     mmol N/m2/d
! retPIn(iiP3)                  P3(PhytoPlankton) uptaken by filter feeders     mmol N/m2/d
! retPIn(iiP4)                  P4(PhytoPlankton) uptaken by filter feeders     mmol N/m2/d
! retPIp(iiP1)                  P1(PhytoPlankton) uptaken by filter feeders     mmol P/m2/d
! retPIp(iiP2)                  P2(PhytoPlankton) uptaken by filter feeders     mmol P/m2/d
! retPIp(iiP3)                  P3(PhytoPlankton) uptaken by filter feeders     mmol P/m2/d
! retPIp(iiP4)                  P4(PhytoPlankton) uptaken by filter feeders     mmol P/m2/d
! retPIl(iiP1)                  P1(PhytoPlankton) uptaken by filter feeders    mmol Si/m2/d
! retPIl(iiP2)                  P2(PhytoPlankton) uptaken by filter feeders    mmol Si/m2/d
! retPIl(iiP3)                  P3(PhytoPlankton) uptaken by filter feeders    mmol Si/m2/d
! retPIl(iiP4)                  P4(PhytoPlankton) uptaken by filter feeders    mmol Si/m2/d
! retPIs(iiP1)                  P1(PhytoPlankton) uptaken by filter feeders    mg Chl /m2/d
! retPIs(iiP2)                  P2(PhytoPlankton) uptaken by filter feeders    mg Chl /m2/d
! retPIs(iiP3)                  P3(PhytoPlankton) uptaken by filter feeders    mg Chl /m2/d
! retPIs(iiP4)                  P4(PhytoPlankton) uptaken by filter feeders    mg Chl /m2/d

        integer,parameter,public :: ppjOAO2o=1, pprutQ6c=2, pprutQ6n=3,&
          pprutQ6p=4, pprutQ6s=5, pprutQ1c=6, pprutQ1n=7, pprutQ1p=8,&
          pprrBTo=9, pprrATo=10, ppreBTn=11, ppreBTp=12, ppreATn=13,&
          ppreATp=14, ppturenh=15, ppirrenh=16, ppshiftD1m=17, ppshiftD2m=18,&
          ppjG2K3o=19, ppjG2K7o=20, ppM1p=21, ppM11p=22, ppM21p=23,&
          ppM4n=24, ppM14n=25, ppM24n=26, ppM3n=27, ppM5s=28,&
          ppM6r=29, ppRIc=30, ppRIn=31, ppRIp=32, ppRIs=33,&
          ppPIc=34, ppPIn=35, ppPIp=36, ppPIs=37, ppjPIY3c=38,&
          ppjRIY3c=39, ppjRIY3n=40, ppjRIY3p=41, ppjRIY3s=42, ppDepth_Ben=43,&
          ppETW_Ben=44, ppO2o_Ben=45, ppN1p_Ben=46, ppN3n_Ben=47, ppN4n_Ben=48,&
          ppN5s_Ben=49, ppN6r_Ben=50, ppjG2O2o=51, ppjK1N1p=52, ppjK3N3n=53,&
          ppjK4N4n=54, ppjK5N5s=55, ppjK6N6r=56, pptotpeln=57, pptotpelp=58,&
          pptotpels=59, pptotbenn=60, pptotbenp=61, pptotbens=62

        integer,public :: ppruHI(iiBenBacteria),&
          ppreHI(iiBenBacteria), ppretPIc(iiPhytoPlankton),&
          ppretPIn(iiPhytoPlankton), ppretPIp(iiPhytoPlankton),&
          ppretPIl(iiPhytoPlankton), ppretPIs(iiPhytoPlankton)

        real(RLEN),public,dimension(:),pointer :: jOAO2o, rutQ6c, rutQ6n,&
          rutQ6p, rutQ6s, rutQ1c, rutQ1n, rutQ1p, rrBTo, rrATo, reBTn, reBTp,&
          reATn, reATp, turenh, irrenh, shiftD1m, shiftD2m, jG2K3o, jG2K7o, M1p,&
          M11p, M21p, M4n, M14n, M24n, M3n, M5s, M6r, RIc,&
          RIn, RIp, RIs, PIc, PIn, PIp, PIs, jPIY3c, jRIY3c,&
          jRIY3n, jRIY3p, jRIY3s, Depth_Ben, ETW_Ben, O2o_Ben, N1p_Ben, N3n_Ben,&
          N4n_Ben, N5s_Ben, N6r_Ben, jG2O2o, jK1N1p, jK3N3n, jK4N4n, jK5N5s,&
          jK6N6r, totpeln, totpelp, totpels, totbenn, totbenp, totbens

        real(RLEN),public,dimension(:,:),pointer :: ruHI, reHI, retPIc, retPIn,&
          retPIp, retPIl, retPIs

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  Other 3d-Global Variables 
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        real(RLEN),public,dimension(:),allocatable             :: &
          Depth,  &   !
          ABIO_eps      ! the abiotic extinction coefficient calculated in sediment models (GOTM/GETM)



        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  Other 2d-Global Variables 
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,public,dimension(:),allocatable             :: &
          KPO4,  &   !
          KNH4,  &   !
          KNO3,  &   !
          KRED,  &   !
          KSIO3,  &   !
          KSIO3E,  &   !
          KQ1      !

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  variables to generate flux_output 
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
           real(RLEN), public,dimension(:),allocatable ::flx_t
           integer,    public,dimension(:),allocatable ::flx_SS
           integer,    public,dimension(:),allocatable ::flx_states
           integer,    public,dimension(:),allocatable ::flx_ostates
           integer,    public,dimension(:),allocatable ::flx_calc_nr
           integer,    public,dimension(:),allocatable ::flx_CalcIn
           integer,    public,dimension(:),allocatable ::flx_option
           integer,    public                          ::flx_cal_ben_start
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !! SHARED GLOBAL FUNCTIONS (must be below contains)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        public flux, flux_vector, Source, Source_D3_vector, Source_D2_vector , &  
               make_flux_output
        public ppPhytoPlankton, ppMesoZooPlankton, ppMicroZooPlankton,&
          PhytoPlankton, MesoZooPlankton, MicroZooPlankton

        public ppBenOrganisms, ppBenDetritus, ppBenBacteria,&
          ppBenthicPhosphate, ppBenthicAmmonium, BenOrganisms, BenDetritus,&
          BenBacteria, BenthicPhosphate, BenthicAmmonium


        contains

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !! Group Pelagic (D3) state functions
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        function ppPhytoPlankton(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          integer  ::ppPhytoPlankton
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiP1)
              i = ppP1c
            CASE (iiP2)
              i = ppP2c
            CASE (iiP3)
              i = ppP3c
            CASE (iiP4)
              i = ppP4c
          END SELECT
          i=i+ (constituent -iiC )

          ppPhytoPlankton = i
        END function

        function ppMesoZooPlankton(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          integer  ::ppMesoZooPlankton
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiZ3)
              i = ppZ3c
            CASE (iiZ4)
              i = ppZ4c
          END SELECT
          i=i+ (constituent -iiC )

          ppMesoZooPlankton = i
        END function

        function ppMicroZooPlankton(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          integer  ::ppMicroZooPlankton
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiZ5)
              i = ppZ5c
            CASE (iiZ6)
              i = ppZ6c
          END SELECT
          i=i+ (constituent -iiC )

          ppMicroZooPlankton = i
        END function

        function PhytoPlankton(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          real(RLEN),dimension(:),pointer  ::PhytoPlankton
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiP1)
              i = ppP1c
            CASE (iiP2)
              i = ppP2c
            CASE (iiP3)
              i = ppP3c
            CASE (iiP4)
              i = ppP4c
          END SELECT
          i=i+ (constituent -iiC )

          PhytoPlankton => D3STATE(i,:)
        END function

        function MesoZooPlankton(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          real(RLEN),dimension(:),pointer  ::MesoZooPlankton
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiZ3)
              i = ppZ3c
            CASE (iiZ4)
              i = ppZ4c
          END SELECT
          i=i+ (constituent -iiC )

          MesoZooPlankton => D3STATE(i,:)
        END function

        function MicroZooPlankton(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          real(RLEN),dimension(:),pointer  ::MicroZooPlankton
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiZ5)
              i = ppZ5c
            CASE (iiZ6)
              i = ppZ6c
          END SELECT
          i=i+ (constituent -iiC )

          MicroZooPlankton => D3STATE(i,:)
        END function


          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !! Group Benthic (D2) state functions
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        function ppBenOrganisms(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          integer  ::ppBenOrganisms
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiY1)
              i = ppY1c
            CASE (iiY2)
              i = ppY2c
            CASE (iiY3)
              i = ppY3c
            CASE (iiY4)
              i = ppY4c
            CASE (iiY5)
              i = ppY5c
          END SELECT
          i=i+ (constituent -iiC )

          ppBenOrganisms = i
        END function

        function ppBenDetritus(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          integer  ::ppBenDetritus
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiQ1)
              i = ppQ1c
            CASE (iiQ11)
              i = ppQ11c
          END SELECT
          i=i+ (constituent -iiC )

          ppBenDetritus = i
        END function

        function ppBenBacteria(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          integer  ::ppBenBacteria
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiH1)
              i = ppH1c
            CASE (iiH2)
              i = ppH2c
          END SELECT
          i=i+ (constituent -iiC )

          ppBenBacteria = i
        END function

        function ppBenthicPhosphate(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          integer  ::ppBenthicPhosphate
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiK1)
              i = ppK1p
            CASE (iiK11)
              i = ppK11p
            CASE (iiK21)
              i = ppK21p
          END SELECT
          i=i+ (constituent -iiP )

          ppBenthicPhosphate = i
        END function

        function ppBenthicAmmonium(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          integer  ::ppBenthicAmmonium
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiK4)
              i = ppK4n
            CASE (iiK14)
              i = ppK14n
            CASE (iiK24)
              i = ppK24n
          END SELECT
          i=i+ (constituent -iiN )

          ppBenthicAmmonium = i
        END function

        function BenOrganisms(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          real(RLEN),dimension(:),pointer  ::BenOrganisms
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiY1)
              i = ppY1c
            CASE (iiY2)
              i = ppY2c
            CASE (iiY3)
              i = ppY3c
            CASE (iiY4)
              i = ppY4c
            CASE (iiY5)
              i = ppY5c
          END SELECT
          i=i+ (constituent -iiC )

          BenOrganisms => D2STATE(i,:)
        END function

        function BenDetritus(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          real(RLEN),dimension(:),pointer  ::BenDetritus
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiQ1)
              i = ppQ1c
            CASE (iiQ11)
              i = ppQ11c
          END SELECT
          i=i+ (constituent -iiC )

          BenDetritus => D2STATE(i,:)
        END function

        function BenBacteria(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          real(RLEN),dimension(:),pointer  ::BenBacteria
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiH1)
              i = ppH1c
            CASE (iiH2)
              i = ppH2c
          END SELECT
          i=i+ (constituent -iiC )

          BenBacteria => D2STATE(i,:)
        END function

        function BenthicPhosphate(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          real(RLEN),dimension(:),pointer  ::BenthicPhosphate
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiK1)
              i = ppK1p
            CASE (iiK11)
              i = ppK11p
            CASE (iiK21)
              i = ppK21p
          END SELECT
          i=i+ (constituent -iiP )

          BenthicPhosphate => D2STATE(i,:)
        END function

        function BenthicAmmonium(n,constituent)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          real(RLEN),dimension(:),pointer  ::BenthicAmmonium
          integer, intent(IN)  ::n
          integer, intent(IN)  ::constituent

          integer i

          SELECT CASE (n)
            CASE (iiK4)
              i = ppK4n
            CASE (iiK14)
              i = ppK14n
            CASE (iiK24)
              i = ppK24n
          END SELECT
          i=i+ (constituent -iiN )

          BenthicAmmonium => D2STATE(i,:)
        END function



          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! flux functions
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          subroutine flux_vector(iiSub,origin,destination,flux)

            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            use constants, ONLY: RLEN, ZERO,  SEC_PER_DAY
            use global_mem, ONLY: LOGUNIT
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            ! Implicit typing is never allowed
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            IMPLICIT NONE

            integer,intent(IN) :: iiSub
            integer,intent(IN) :: origin
            integer,intent(IN) :: destination
            real(RLEN),intent(IN) :: flux(:)

            integer :: i
            character :: D23*8

            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            !BEGIN compute
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            TESTNANVECTOR(flux,origin,destination)
            CHECKFLUX(-1,iiSub,origin,destination)

            IF ( origin /= destination )  THEN
              IF ( minval(flux) < ZERO) then
                do i=1,size(flux)
                  if (flux(i)< 0.0D+00) write(LOGUNIT,'(''at level:'',I4)') i
                enddo
                D23="Pelagic"
                if ( iiSub == iiBen) D23="Benthic"
                write(LOGUNIT,'(''In '',A,'':origin='',i4,'' destination='',i4)') &
                  D23, origin,destination
                write(LOGUNIT,'(''In '',A,'':origin='',i4,'' &
                  destination='',i4)') D23, origin,destination
                write(LOGUNIT,'(''flux='',(G16.8))') flux
                STDERR  "Error in flux_vector function: negative flux !"
#if DEBUG && BFM_GOTM
                call gotm_error("flux_vector","negative flux")
#endif
              ENDIF ! minval<0
              SELECT CASE ( iiSub )
                case (iiPel)
                  D3SINK(origin,destination,:)  =  flux/SEC_PER_DAY
                  D3SOURCE(destination,origin,:)=  flux/SEC_PER_DAY
                case (iiBen)
                  D2SINK(origin,destination,:) =  flux/SEC_PER_DAY
                  D2SOURCE(destination,origin,:)   = flux/SEC_PER_DAY
              end select
            ELSE
              SELECT CASE ( iiSub )
                CASE (iiPel)
                  WHERE (flux > 0.0D+00 )
                    D3SOURCE(origin,destination,:) =D3SOURCE(origin,destination,:) &
                      + flux/SEC_PER_DAY
                  ELSEWHERE
                    D3SINK(destination,origin,:) =D3SINK(destination,origin,:) - &
                      flux/SEC_PER_DAY
                  END WHERE
                CASE (iiBen)
                  WHERE (flux > 0.0D+00 )
                    D2SOURCE(destination,origin,:) =D2SOURCE(destination,origin,:) &
                      + flux/SEC_PER_DAY
                  ELSEWHERE
                    D2SINK(origin,destination,:) =D2SINK(origin,destination,:) - &
                      flux/SEC_PER_DAY
                  END WHERE
              end select
            ENDIF !origin <> destination

            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            !END compute
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            RETURN
          END subroutine flux_vector

#ifdef DEBUG
          subroutine testnan_vector(array,origin,destination)
          use global_mem, ONLY: LOGUNIT

            real(RLEN),intent(IN)    :: array(:)
            integer,intent(IN) :: origin
            integer,intent(IN) :: destination
            integer:: i=0
            do i=1,size(array)
              if (isnan(array(i))== .TRUE. ) then
                write(LOGUNIT,'(''at level:'',I4)') i
                write(LOGUNIT,'(''origin='',i4,'' destination='',i4)') &
                  origin,destination
                write(*,'(''at level:'',I4)') i
                write(*,'(''origin='',i4,'' destination='',i4)') &
                  origin,destination
                STDERR 'Nan value in flux'
                stop 1002
              endif
            enddo
          END subroutine testnan_vector

          subroutine testnan(scalar,origin,destination)

            real(RLEN),intent(IN)    :: scalar
            integer,intent(IN) :: origin
            integer,intent(IN) :: destination
            if (isnan(scalar)== .TRUE. ) then
            write(*,'(''origin='',i4,'' destination='',i4)') origin,destination
            STDERR 'Nan value in scalar flux'
            stop 1003
          endif
        END subroutine testnan
#endif

        subroutine flux(grid_nr,iiSub,origin,destination,flow)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          use constants, ONLY: RLEN, ZERO, SEC_PER_DAY
          use global_mem, ONLY: LOGUNIT
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          integer,intent(IN)    :: grid_nr
          integer,intent(IN)    :: iiSub
          integer,intent(IN)    :: origin
          integer,intent(IN)    :: destination
          real(RLEN),intent(IN) :: flow

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          character    :: D23*8
          !BEGIN compute
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          TESTNAN(flow,origin,destination)
          CHECKFLUX(grid_nr,iiSub,origin,destination)

          IF ( origin /= destination ) THEN
            IF ( flow < ZERO) then
              D23="Pelagic"
              if ( iiSub == iiBen) D23="Benthic"
              write(LOGUNIT,'(''In '',A,'':origin='',i4,'' destination='',i4)') &
                D23, origin,destination
              write(LOGUNIT,*) "Error in (scalar) vector  function: negative flux!"
              write(LOGUNIT,*) "origin,destination:", origin,destination
              write(LOGUNIT,*) flow
              STDERR "Error in (scalar)flux function:negative flux !"
#if DEBUG && BFM_GOTM
              call gotm_error("flux","negative flux")
#endif
            ENDIF ! flow<0
            SELECT CASE ( iiSub )
              CASE (iiPel)
                D3SINK(origin,destination,grid_nr)=  flow/SEC_PER_DAY
                D3SOURCE(destination,origin,grid_nr)= flow/SEC_PER_DAY
              CASE (iiBen)
                D2SINK(origin,destination,grid_nr)=  flow/SEC_PER_DAY
                D2SOURCE(destination,origin,grid_nr)= flow/SEC_PER_DAY
            END SELECT
          ELSE
            SELECT CASE ( iiSub )
              CASE (iiPel)
                if (flow > 0.0 ) then
                  D3SOURCE(destination,origin,grid_nr)= & 
                    D3SOURCE(destination,origin,grid_nr) + flow/SEC_PER_DAY
                ELSE
                  D3SINK(origin,destination,grid_nr)= &
                     D3SINK(origin,destination,grid_nr) - flow/SEC_PER_DAY
                ENDIF
              CASE (iiBen)
                if (flow > 0.0 ) then
                  D2SOURCE(destination,origin,grid_nr)= &
                    D2SOURCE(destination,origin,grid_nr) + flow/SEC_PER_DAY
                ELSE
                  D2SINK(origin,destination,grid_nr)= &
                    D2SINK(origin,destination,grid_nr) - flow/SEC_PER_DAY
                ENDIF
            END SELECT
          ENDIF
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !END compute
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          return
        end subroutine flux

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! vector function to get actual rate of change in the pelagic
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        function Source_D3_vector(iistate)
          use constants, ONLY: RLEN, ZERO, SEC_PER_DAY

          IMPLICIT NONE

          integer, intent(IN) ::iistate
          real(RLEN) ::Source_D3_vector(size(D3SOURCE,DIM=3))
          ! Array in sum is by sum seen as 2D-array: DIM=1 and NOT 2
          Source_D3_vector=(sum(D3SOURCE(iistate,:,:),DIM=1)- &
             sum(D3SINK(iistate,:,:),DIM=1))*SEC_PER_DAY
        end function Source_D3_vector

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! vector function to get actual rate of change in the benthic
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        function Source_D2_vector(iistate)
          use constants, ONLY: RLEN, ZERO, SEC_PER_DAY

          IMPLICIT NONE

          integer, intent(IN) ::iistate
          real(RLEN) ::Source_D2_vector(size(D2SOURCE,DIM=3))
          ! Array in sum is by sum seen as 2D-array: DIM=1 and NOT 2
          Source_D2_vector=(sum(D2SOURCE(iistate,:,:),DIM=1)- &
            sum(D2SINK(iistate,:,:),DIM=1))*SEC_PER_DAY
        end function Source_D2_vector
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! function to get actual rate of change
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        function source(iiSub,iibox,iistate)
          use constants, ONLY: RLEN, ZERO, SEC_PER_DAY

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          real(RLEN)  ::Source
          integer, intent(IN)  ::iiSub
          integer, intent(IN)  ::iibox
          integer, intent(IN)  ::iistate
          IF ( iiSub == iiPel )  THEN
            Source = (sum(D3SOURCE(iistate,:,iibox))- &
              sum(D3SINK(iistate,:,iibox)))*SEC_PER_DAY
          ELSEIF ( iiSub == iiBen )  THEN
            Source = (sum(D2SOURCE(iistate,:,iibox))- &
              sum(D2SINK(iistate,:,iibox)))*SEC_PER_DAY
          ENDIF
        end function source

        subroutine unicflux(grid_nr,iiSub,origin,destination)
        use constants, ONLY: RLEN

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          use global_mem, ONLY: LOGUNIT
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          IMPLICIT NONE

          integer,intent(IN)    :: grid_nr
          integer,intent(IN)    :: origin
          integer,intent(IN)    :: iiSub
          integer,intent(IN)    :: destination

          real(RLEN) :: tot
          character*20:: type

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !BEGIN compute
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


          SELECT CASE ( iiSub )
            CASE (iiPel)
              type="D3"
              IF ( grid_nr <=0  ) THEN
                tot=sum(D3SINK(origin,destination,:))
              ELSE
                tot=D3SINK(origin,destination,grid_nr)
              ENDIF
            CASE (iiBen)
              type="D2"
              IF ( grid_nr <=0  ) THEN
                tot=sum(D2SINK(origin,destination,:))
              ELSE
                tot=D2SINK(origin,destination,grid_nr)
              ENDIF
            CASE (iiReset)
              D3SINK(:,:,:)=0.0D+00
              D2SINK(:,:,:)=0.0D+00
          END SELECT
          IF ( tot > 0.0D+00  ) THEN
            write(LOGUNIT,'(''Double defintion '',A2,''-flux'')')type
            write(LOGUNIT,'(''origin:'',I3,'' destination:'',I3)') origin, destination
            IF ( origin /= destination ) THEN
              STDERR 'double definition of fluxes'
              stop 1006
            ENDIF
          ENDIF
        !END compute
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        return
      end subroutine unicflux

      subroutine make_flux_output(mode, nr0,nlev, h, out)
      use constants, ONLY: RLEN, SEC_PER_DAY

      implicit none
      integer,intent(IN)                  ::mode
      integer,intent(IN)                  ::nr0
      integer,intent(IN)                  ::nlev
      real(RLEN),intent(IN),dimension(0:nlev)  :: h
      real(RLEN),intent(OUT),dimension(0:nlev)  :: out


      integer      ::nr
      integer      ::i
      integer      ::k
      integer      ::klev
      real(RLEN),dimension(0:nlev):: tot
      real(RLEN)   :: r

      nr=nr0;if ( mode == 2 ) nr=nr0+flx_cal_ben_start
      klev=nlev ; if ( flx_CalcIn(nr) == iiBen)  klev=1 
      out(:)=0.0
      if ( flx_CalcIn(nr) == iiBen) then
        do i=flx_calc_nr(nr-1)+1,flx_calc_nr(nr)
          if (flx_SS(i) ==1 ) then
             out(1:klev)= out(1:klev) + flx_t(i) * D2SINK(flx_states(i),flx_ostates(i),:)
          else
             out(1:klev)= out(1:klev) + flx_t(i) * D2SOURCE(flx_states(i),flx_ostates(i),:)
          endif
        enddo
      else
        do i=flx_calc_nr(nr-1)+1,flx_calc_nr(nr)
          if (flx_SS(i) ==1 ) then
             out(1:klev)= out(1:klev) + flx_t(i) * D3SINK(flx_states(i),flx_ostates(i),:)
          else
             out(1:klev)= out(1:klev) + flx_t(i) * D3SOURCE(flx_states(i),flx_ostates(i),:)
          endif
        enddo
      endif

      out=out*SEC_PER_DAY        

      SELECT CASE ( flx_option(nr) )
        CASE(1) !Specific rate
          tot=1.0D-80;
          k=0
          if ( flx_CalcIn(nr) == iiBen) then
             do i=flx_calc_nr(nr-1)+1,flx_calc_nr(nr)
               if ( k.ne. flx_states(i) ) then
                  k=flx_states(i)
                  tot(1:klev)=tot(1:klev) + D2STATE(k,:)
               endif
              enddo
          else 
             do i=flx_calc_nr(nr-1)+1,flx_calc_nr(nr)
               if ( k.ne. flx_states(i) ) then
                 k=flx_states(i)
                 tot(1:klev)=tot(1:klev) + D3STATE(k,:)
               endif
             enddo
          endif
          out(1:klev)=out(1:klev)/tot(1:klev)
        CASE(2) ! summing the column :perm2
          ! d3 -->d2 var.
          r=sum(out(1:klev) *h(1:klev))
          out=0.0
          out(1)=r
      END SELECT


      return
      end subroutine make_flux_output
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! end of contain section
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  end module

