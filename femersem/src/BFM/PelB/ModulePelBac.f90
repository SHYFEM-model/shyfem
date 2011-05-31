!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelBac
!
! DESCRIPTION
!   Here are the values for the pelagic bacteria.
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_PelBac
!
! !USES:

  use global_mem

!  
!
! !AUTHORS
!   Original version by J.W. Baretta and Piet Ruardij
!
!
!
! !REVISION_HISTORY
!   2006-03-17 Piet Ruardij : density dependent mortality
!  rua: density dependent mortality set on zero now!!
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
  ! PelBac PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: p_version  ! Switch for DOM uptake parameterization
  !  p_version=1 <LUCA> Polimenes version
  !  p_version=2 <BFM> version
  real(RLEN)  :: p_q10  ! Q10-value (temperature dependency)
  real(RLEN)  :: p_chdo  ! Michaelis const for O2 dependence (mmol/m3)
  real(RLEN)  :: p_sd  ! Independent specific mortality (1/d)
  real(RLEN)  :: p_sd2  ! Density dependent mortality (value: 0.009) (1/d)
  real(RLEN)  :: p_suR1  ! Specific potential DOM availability (1/d)
  real(RLEN)  :: p_suR2  ! Specific potential DOM availability (1/d)
  real(RLEN)  :: p_suR6  ! Availability of POM (1/d)
  real(RLEN)  :: p_sum  ! Specific potential uptake (1/d)
  real(RLEN)  :: p_pu_ra  ! Activity respiration (-)
  real(RLEN)  :: p_pu_ra_o  ! Decrease in Ass. efficiency at low O2 conc (-).
  real(RLEN)  :: p_srs  ! Specific rest respiration (1/day)
  real(RLEN)  :: p_qnc  ! Optimal N/C ratio (model units) 45:9:1
  real(RLEN)  :: p_qpc  ! Optimal P/C ratio (model units) C:N:P
  real(RLEN)  :: p_qlnc  ! Minimal N/C ratio (model units) 45:9:1 <BFM>
  real(RLEN)  :: p_qlpc  ! Minimal P/C ratio (model units) C:N:P (BFM>
  real(RLEN)  :: p_qun  ! nutrient affinity ( mmol/mgC/day) <BFM>
  real(RLEN)  :: p_qup  ! nutrient affinity ( mmol/mgC/day) <BFM>
  real(RLEN)  :: p_lN4  ! ammonium conc. at which nutrate uptake are equal (BFM)
  real(RLEN)  :: p_pu_ea_R7  ! excretion of reg. org. met. (-) <LUCA>
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitPelBac
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitPelBac()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /PelBac_parameters/ p_version, p_q10, p_chdo, p_sd, p_sd2, p_suR1, &
    p_suR2, p_suR6, p_sum, p_pu_ra, p_pu_ra_o, p_pu_ea_R7, p_srs, p_qpc, p_qlpc, &
    p_qnc, p_qlnc, p_qun, p_qup, p_lN4
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading PelBac parameters.."
open(NMLUNIT,file='PelBac.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=PelBac_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=PelBac_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitPelBac.f90","PelBac.nml")
101 call error_msg_prn(NML_READ,"InitPelBac.f90","PelBac_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitPelBac
  end module mem_PelBac
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
