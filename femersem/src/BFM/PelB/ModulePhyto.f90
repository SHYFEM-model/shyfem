!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Phyto
!
! DESCRIPTION
!   Parameter values for the phytoplankton groups
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_Phyto
!
! !USES:

  use global_mem
  use mem,  ONLY: iiPhytoPlankton, iiP1, iiP2, iiP3, iiP4

!  
!
! !AUTHORS
!   the ERSEM group, Marcello Vichi, JWB, HBB
!
!
!
! !REVISION_HISTORY
!   !
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
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Phyto PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  !  ---------------- Physiological parameters -----------------
  !
  real(RLEN)  :: p_q10(iiPhytoPlankton)  ! Doubling temperature
  real(RLEN)  :: p_sum(iiPhytoPlankton)  ! Maximal productivity at 10 degrees C
  real(RLEN)  :: p_srs(iiPhytoPlankton)  ! Respiration rate at 10 degrees C
  real(RLEN)  :: p_sdmo(iiPhytoPlankton)  ! Max.specific nutrient-stress lysis rate
  real(RLEN)  :: p_thdo(iiPhytoPlankton)  ! Half value for nutrient stress lysis
  real(RLEN)  :: p_seo(iiPhytoPlankton)  ! Extra lysis rate for P4
  real(RLEN)  :: p_pu_ea(iiPhytoPlankton)  ! Fraction of pp excreted as PLOC/PDET
  real(RLEN)  :: p_pu_ra(iiPhytoPlankton)  ! Activity respiration rate
  !
  !  ---------------- Nutrient parameters in phytoplankton -----------------
  !
  integer  :: p_limnut(iiPhytoPlankton)  ! switch for nut. limitation (Liebig is default)
  real(RLEN)  :: p_qnlc(iiPhytoPlankton)
  real(RLEN)  :: p_qnRc(iiPhytoPlankton)
  real(RLEN)  :: p_xqn(iiPhytoPlankton)
  real(RLEN)  :: p_qplc(iiPhytoPlankton)
  real(RLEN)  :: p_qpRc(iiPhytoPlankton)
  real(RLEN)  :: p_xqp(iiPhytoPlankton)
  real(RLEN)  :: p_chPs(iiPhytoPlankton)  ! half-value of SIO4-lim (mmol Si m-3)
  real(RLEN)  :: p_qslc(iiPhytoPlankton)  ! Minimum quotum Si in PI
  real(RLEN)  :: p_qsRc(iiPhytoPlankton)  ! Reference quotum Si in PI
  real(RLEN)  :: p_xqs(iiPhytoPlankton)  ! factor for max quotum S
  real(RLEN)  :: p_qus(iiPhytoPlankton)  ! affinity of PI for Si
  real(RLEN)  :: p_qun(iiPhytoPlankton)
  real(RLEN)  :: p_qup(iiPhytoPlankton)
  real(RLEN)  :: p_lN4(iiPhytoPlankton)
  real(RLEN)  :: p_alpha_chl(iiPhytoPlankton)  ! Initial slope P-I curve
  !
  !  ------------- Chlorophyll parameters -----------
  !  skel: Skeletonema costatum pav: Pavlova lutheri
  !  syn: Synechoccus sp. (significant alpha decrease with irradiance)
  !  gyr: Gyrodinium sp. iso: Isochrysis galbana
  !              skel     iso      syn      gyr
  real(RLEN)  :: p_sdchl(iiPhytoPlankton)  ! Specific turnover rate for Chla [d-1]
  ! p_qchlc =    0.03,    0.025,   0.1,     0.02    # Maximum quotum Chla:C
  !             +-0.024  +-0.001  +-0.003  +-0.004
  !  Thalassiosira sp. [0.05+-0.01]
  !              skel     pav      syn      gyr
  real(RLEN)  :: p_esNI(iiPhytoPlankton)  ! Nutrient stress threshold for Sinking
  ! p_alpha_chl = 1.0e-5, 0.46e-5*2.0, 2.0e-5, 0.68e-5 # Initial slope P-I curve
  !  Thalassiosira sp. [0.48-0.63]
  real(RLEN)  :: p_res(iiPhytoPlankton)  ! Sinking velocity (m/d)
  !   p_qchlc  = 0.05,      0.03,      0.07,      0.02 # Maximum quotum Chla:C
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitPhyto
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitPhyto()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /Phyto_parameters/ p_q10, p_sum, p_srs, p_sdmo, p_seo, p_pu_ea, &
    p_pu_ra, p_qnlc, p_qplc, p_qslc, p_qnRc, p_qpRc, p_qsRc, p_qun, p_qup, &
    p_qus, p_xqn, p_xqp, p_xqs, p_esNI, p_thdo, p_res, p_chPs, p_lN4, &
    p_limnut, p_alpha_chl, p_sdchl
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading Phyto parameters.."
open(NMLUNIT,file='Phyto.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=Phyto_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=Phyto_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitPhyto.f90","Phyto.nml")
101 call error_msg_prn(NML_READ,"InitPhyto.f90","Phyto_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitPhyto
  end module mem_Phyto
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
