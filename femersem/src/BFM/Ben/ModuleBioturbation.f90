!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Bioturbation
!
! DESCRIPTION
!   Here are the values belonging to the benthic submodel.
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_Bioturbation
!
! !USES:

  use global_mem

!  
!
! !AUTHORS
!   Original version by W. Ebenhoeh and C Kohlmeier.
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
  ! Bioturbation PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: p_q10  ! -       # temperature dependency
  real(RLEN)  :: p_chtur  ! mgC/m2  # MM-constant describing control
  !  of biomass YI on bioturbation
  real(RLEN)  :: p_chirr  ! mgC/m2  # MM-constant describing control
  !  of biomass YI on bioirrigation
  real(RLEN)  :: p_turY1  ! -       # Multiplication factor to correct
  !  influence of a group of YI on
  !  bioturbation
  real(RLEN)  :: p_irrY4  ! -       # Multiplication factor to correct
  !  influence of a group of YI on
  !  bioirrigation
  real(RLEN)  :: p_Etur  ! m2/d    # Potential turbation constant
  real(RLEN)  :: p_cturm  ! m       # Average depth where the
  !  bioturbation activity is
  !  located
  real(RLEN)  :: p_cmtur  ! -       # maximum enhancement of detritus
  !  penetration due to
  !  bioturbation
  real(RLEN)  :: p_cmirr  ! -       # maximum enhancement of diffusion
  !  due to bioirrigation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitBioturbation
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitBioturbation()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /Bioturbation_parameters/ p_q10, p_chtur, p_chirr, p_turY1, p_irrY4, &
    p_Etur, p_cturm, p_cmtur, p_cmirr
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading Bioturbation parameters.."
open(NMLUNIT,file='Bioturbation.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=Bioturbation_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=Bioturbation_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitBioturbation.f90","Bioturbation.nml")
101 call error_msg_prn(NML_READ,"InitBioturbation.f90","Bioturbation_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitBioturbation
  end module mem_Bioturbation
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
