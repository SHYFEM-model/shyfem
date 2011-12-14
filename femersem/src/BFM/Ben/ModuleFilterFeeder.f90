!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: FilterFeeder
!
! DESCRIPTION
!   Parameter value file for suspension feeders (Y3)
!   
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_FilterFeeder
!
! !USES:

  use global_mem

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
  ! FilterFeeder PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: p_dwat  ! Food layer in the water
  real(RLEN)  :: p_su  ! Growth rate
  real(RLEN)  :: p_q10  ! q10
  real(RLEN)  :: p_R6  ! Food matrix detritus on pelagic det. (R6)
  real(RLEN)  :: p_clu  ! Lower limit of availability of a food source
  real(RLEN)  :: p_chu  ! Upper limit of total food
  real(RLEN)  :: p_pue  ! Excreted fraction of uptake
  real(RLEN)  :: p_pueQ6  ! Excreted fraction of uptake
  real(RLEN)  :: p_pudsil  ! Selectvity for nutrients
  real(RLEN)  :: p_sr  ! Relative respiration rate
  real(RLEN)  :: p_pur  ! Part of uptake which is used for respiration
  real(RLEN)  :: p_sd  ! Specific Mortality
  real(RLEN)  :: p_qn  ! Fixed nutrient quotum N:C
  real(RLEN)  :: p_qp  ! Fixed nutrient quotum P:C
  real(RLEN)  :: p_clm  ! Upper depth of accessed sediment layer
  real(RLEN)  :: p_cm  ! Lower  depth of accessed sediment layer
  real(RLEN)  :: p_puQ6  ! Food matrix detritus on the sediment (Q6)
  real(RLEN)  :: p_PI  ! Food matrix Y3 on diatoms
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitFilterFeeder
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitFilterFeeder()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /FilterFeeder_parameters/ p_dwat, p_su, p_q10, p_R6, p_clu, p_chu, &
    p_pue, p_pueQ6, p_pudsil, p_sr, p_pur, p_sd, p_qn, p_qp, p_clm, p_cm, p_puQ6, &
    p_PI
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading FilterFeeder parameters.."
open(NMLUNIT,file='FilterFeeder.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=FilterFeeder_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=FilterFeeder_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitFilterFeeder.f90","FilterFeeder.nml")
101 call error_msg_prn(NML_READ,"InitFilterFeeder.f90","FilterFeeder_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitFilterFeeder
  end module mem_FilterFeeder
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
