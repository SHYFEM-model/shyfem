!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MicroZoo
!
! DESCRIPTION
!   The parameter values for the microzooplankton submodel.
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_MicroZoo
!
! !USES:

  use global_mem
  use mem,  ONLY: iiMicroZooPlankton, iiZ5, iiZ6

!  
!
! !AUTHORS
!   ERSEM group, Hanneke Baretta-Bekker
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
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! MicroZoo PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  !  The variable Z5 represents microozooplankton and Z6 represents
  !  heterotrophic nanoflagellates (HNAN)
  !
  !
  real(RLEN)  :: p_q10(iiMicroZooPlankton)  ! Q10 value
  real(RLEN)  :: p_srs(iiMicroZooPlankton)  ! Respiration rate at 10 degrees Celsius
  real(RLEN)  :: p_sum(iiMicroZooPlankton)  ! Max. rel daily uptake as a fraction of biomass
  real(RLEN)  :: p_sdo(iiMicroZooPlankton)  ! Mortality due to oxygen limitation
  real(RLEN)  :: p_sd(iiMicroZooPlankton)  ! Temperature independent mortality
  real(RLEN)  :: p_pu_ra(iiMicroZooPlankton)  ! Activity respiration
  real(RLEN)  :: p_pu_ea(iiMicroZooPlankton)  ! Activity excretion
  real(RLEN)  :: p_chro(iiMicroZooPlankton)  ! Oxygen saturation where respiration is 0.5
  real(RLEN)  :: p_chuc(iiMicroZooPlankton)  ! Food concentration where total uptake rate is 0.5
  real(RLEN)  :: p_minfood(iiMicroZooPlankton)  ! Concentration below which feeding on a particular
  !  foodsource is depressed
  real(RLEN)  :: p_suB1(iiMicroZooPlankton)  ! /day   #relative B1 uptake by zoo
  real(RLEN)  :: p_qn_mz(iiMicroZooPlankton)  ! Maximum quotum P
  real(RLEN)  :: p_qp_mz(iiMicroZooPlankton)  ! Maximum quotum N
  real(RLEN)  :: p_suPI(iiMicroZooPlankton,4)  ! /day   #relative P1 uptake by zoo
  real(RLEN)  :: p_suZI(iiMicroZooPlankton,2)  ! /day   #relative Z5 uptake by zoo
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitMicroZoo
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitMicroZoo()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /MicroZoo_parameters/ p_q10, p_srs, p_sum, p_sdo, p_sd, p_pu_ra, &
    p_pu_ea, p_chro, p_chuc, p_minfood, p_suPI, p_suZI, p_suB1, p_qp_mz, p_qn_mz
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading MicroZoo parameters.."
open(NMLUNIT,file='MicroZoo.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=MicroZoo_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=MicroZoo_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitMicroZoo.f90","MicroZoo.nml")
101 call error_msg_prn(NML_READ,"InitMicroZoo.f90","MicroZoo_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitMicroZoo
  end module mem_MicroZoo
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
