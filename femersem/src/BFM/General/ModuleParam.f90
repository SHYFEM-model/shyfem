!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Param
!
! DESCRIPTION
!   List of global model parameters 
!      (global variables which can be changed during the model initialization

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  MODULE mem_Param

!
! !USES:

  USE global_mem
  USE constants
  USE mem, ONLY: iiPhytoPlankton, iiMesoZooPlankton, &
    iiMicroZooPlankton, iiBenOrganisms, iiBenDetritus, iiBenBacteria, &
    iiBenthicPhosphate, iiBenthicAmmonium

!  
!
! !AUTHORS
!   mfstep/ERSEM team
!
! !REVISION_HISTORY
!   --------
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
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  PUBLIC

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Global Model Parameters
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! 0d-parameter used in pelagic submodel
  real(RLEN)   :: &
      p_small=1.0D-80  ,  &
      p_q10diff=1.49  ,  &  ! Temperature-dependency porewater diffusion
      p_qro=0.5  ,  &  ! stoichiometry O2-->S2-
      p_qon_dentri=1.25  ,  &  ! stoichiometry O2-->N denitrification
      p_qon_nitri=2.0  ,  &  ! stoichiometry O2-->N nitrification
      p_clDxm=0.001  ! minimal value of D?.m for calculation of the alpha
  !  alpha is used in expo.func and values > 1/0.001 leadt

  ! 0d-parameter used in pelagic submodel
  logical   :: &
      CalcPelagicFlag=.TRUE.  ! Switch for pelagic system
  ! 0d-parameter used in pelagic submodel
  integer   :: &
      CalcBenthicFlag=3  ! Switch for benthic system

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Allocate the logical flags for switch on the LFG
  !  Initialize to TRUE (overwritten by the namelist values)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  logical   :: CalcPhytoPlankton(iiPhytoPlankton) = .TRUE.
  logical   :: CalcMicroZooPlankton(iiMicroZooPlankton) = .TRUE.
  logical   :: CalcMesoZooPlankton(iiMesoZooPlankton) = .TRUE.
  logical   :: CalcBenOrganisms(iiBenOrganisms) = .TRUE.
  logical   :: CalcBenBacteria(iiBenBacteria) = .TRUE.
  logical   :: CalcBacteria = .TRUE.

  logical   :: &
      CalcPelChemistry=.TRUE.  ,  &  !
      AssignPelBenFluxesInBFMFlag=.TRUE.  ,  &  ! Switches to make choice to define boundary
      AssignAirPelFluxesInBFMFlag=.TRUE.        ! fluxes in physical of biological model
  !      %dim2D%

  ! 1d-parameter used in benthic submodel
  real(RLEN),public,dimension(:),allocatable   :: &
      p_p_ae  ,  &
      p_poro
  ! 0d-parameter used in pelagic submodel
  integer   :: &
      seqnr_cloud=2
  ! 0d-parameter used in pelagic submodel
  real(RLEN)   :: &
      XLatitude=54.0  ,  &  ! Latitude
      p_PAR=0.50  ! Photosynthetically available radiation
  ! 0d-parameter used in pelagic submodel
  integer   :: &
      ChlLightFlag=2  ,  &  ! Switch between light prop.(=1) or Chla.(=2) as a state
      LightForcingFlag=1  ! Switch between instantaneous light and day light average
  ! 1d-parameter used in pelagic submodel
  real(RLEN)   :: &
      p_qchlc(iiPhytoPlankton)  ! Fixed/Maximum quotum Chla:C dependent on ChlLightFlag [mg Chla (mg C)-1]
  ! %p_qchlc%   0.05, 0.03,  0.07,  0.02

  ! 0d-parameter used in pelagic submodel
  real(RLEN)   :: &
      p_eps0=0.04  ,  &  ! Background extinction (abiotic)
      p_epsESS=0.04e-3  ,  &  ! Inorg. suspended matter extinction coeff. (abiotic)
      p_InitSink=100.0  ,  &  ! parameter to Initialize BenthicSInk var.
      p_d_tot=0.30  ,  &  ! m # Thickness of modelled benthic sediment layers
      p_clD1D2m=0.01  ,  &  ! m # minimum distancebetween D1m and D2m
      p_pe_R1c=0.60  ,  &  ! Fraction of excretion going to PLOC
      p_pe_R1n=0.72  ,  &  ! Fraction of excretion going to PLOC
      p_pe_R1p=0.832  ,  &  ! Fraction of excretion going to PLOC
      p_pe_R1s=0.06  ,  &
      p_epsChla=10.0e-3  ,  &  ! Chla-contribution to extinction
      p_epsR6=0.1e-3  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitParam
  
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitParam()
  use mem
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /Param_parameters/ p_small, p_q10diff, p_qro, p_qon_dentri,      &
    p_qon_nitri, p_clDxm, CalcPelagicFlag, CalcBenthicFlag,                 &
    CalcPhytoPlankton,CalcMicroZooPlankton,                                 &
    CalcPelChemistry,CalcMesoZooPlankton,CalcBenOrganisms,CalcBenBacteria,  &
    CalcBacteria, AssignPelBenFluxesInBFMFlag, AssignAirPelFluxesInBFMFlag, &
    p_PAR, ChlLightFlag, LightForcingFlag, p_qchlc, p_eps0, p_epsESS,       &
    p_InitSink, p_d_tot, p_clD1D2m, p_pe_R1c, p_pe_R1n, p_pe_R1p, p_pe_R1s, &
    p_epsChla, p_epsR6
   integer :: i
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !BEGIN compute

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading Param parameters.."
   open(NMLUNIT,file='Param.nml',status='old',action='read',err=100)
   read(NMLUNIT,nml=Param_parameters,err=101)
   close(NMLUNIT)
   write(LOGUNIT,*) "#  Namelist is:"
   write(LOGUNIT,nml=Param_parameters)
  

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitParam.f90","Param.nml")
101 call error_msg_prn(NML_READ,"InitParam.f90","Param_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitParam

  end module
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
