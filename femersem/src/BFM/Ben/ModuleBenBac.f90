!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenBac
!
! DESCRIPTION
!   The parameter value file for benthic bacteria (H1,H2)   
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_BenBac
!
! !USES:

  use global_mem
  use mem,  ONLY: iiBenBacteria, iiH1, iiH2

!  
!
! !AUTHORS
!   ERSEM group, HBB   
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
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! BenBac PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: p_cdm(iiBenBacteria)  ! Half mortality layer        (m)
  real(RLEN)  :: p_qnc(iiBenBacteria)  ! Optimal Internal N quota
  real(RLEN)  :: p_qpc(iiBenBacteria)  ! Optimal Internal P quota
  real(RLEN)  :: p_qlnc(iiBenBacteria)  ! Optimal Internal N quota
  real(RLEN)  :: p_qlpc(iiBenBacteria)  ! Optimal Internal P quota
  real(RLEN)  :: p_q10(iiBenBacteria)  ! Q10
  real(RLEN)  :: p_suhQ6(iiBenBacteria)  ! Specific (high) uptake rate (1/d)
  real(RLEN)  :: p_sulQ6(iiBenBacteria)  ! Specific (slow) uptake rate (1/d)
  real(RLEN)  :: p_sum(iiBenBacteria)  ! Potential uptake rate       (1/d)
  real(RLEN)  :: p_pue(iiBenBacteria)  ! Fraction of Q6 degradated as Q1
  real(RLEN)  :: p_suQ1(iiBenBacteria)  ! Specific uptakte rate of Q1 (1/d)
  real(RLEN)  :: p_pur(iiBenBacteria)  ! Fraction of uptake respired
  real(RLEN)  :: p_srr(iiBenBacteria)  ! Specific respiration        (1/d)
  real(RLEN)  :: p_sumKIn(iiBenBacteria)  ! max. uptake of KIn (mmN/m2)
  real(RLEN)  :: p_sumKIp(iiBenBacteria)  ! max. uptake of KIp (mmP/m2)
  real(RLEN)  :: p_sd(iiBenBacteria)  ! Specific mortality          (1/d)
  integer  :: p_iK4(iiBenBacteria)  ! BenthicAmmonium(p_IK4) =1,2 --> K4n K14n
  integer  :: p_iK1(iiBenBacteria)  ! BenthicPhosphate(p_IK1) =1,2 --> K1p K11p
  integer  :: p_iQ1(iiBenBacteria)  ! BenDetritus(p_IQ1) =1,2 --> Q1 Q11
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitBenBac
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitBenBac()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /BenBac_parameters/ p_cdm, p_qpc, p_qlpc, p_qnc, p_qlnc, p_q10, &
    p_suhQ6, p_sulQ6, p_sum, p_pue, p_suQ1, p_pur, p_srr, p_sumKIn, p_sumKIp, &
    p_sd, p_iK4, p_iK1, p_iQ1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading BenBac parameters.."
open(NMLUNIT,file='BenBac.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=BenBac_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=BenBac_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitBenBac.f90","BenBac.nml")
101 call error_msg_prn(NML_READ,"InitBenBac.f90","BenBac_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitBenBac
  end module mem_BenBac
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
