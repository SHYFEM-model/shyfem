#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelChem
!
! DESCRIPTION
!   !    This process describes the additional dynamics of dissolved
!       compounds in the watercolumn. Parameterized processes are:
!       - nitrification
!       - denitrification
!       - reoxidation of reduction equivalents
!        - dissolution of biogenic silica
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PelChemDynamics
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: N4n, N3n, O2o, O4n, &
  ! N6r, R6s, N5s, P1s
  ! The following Pelagic 1-d global boxvars are modified : flN3O4n
  ! The following Pelagic 1-d global boxvars  are used: ETW, flPTN6r, flP1R6s
  ! The following 0-d global box parametes are used: p_qon_nitri, p_qro, &
  ! p_qon_dentri
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem,  ONLY: N4n, N3n, O2o, O4n, N6r, R6s, N5s, P1s
#endif
  use mem, ONLY: ppN4n, ppN3n, ppO2o, ppO4n, ppN6r, ppR6s, ppN5s, &
    ppP1s, flN3O4n, ETW, flPTN6r, flP1R6s, NO_BOXES, iiBen, iiPel, flux_vector
  use mem_Param,  ONLY: p_qon_nitri, p_qro, p_qon_dentri
  use mem_PelChem


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:MM_vector, eTq_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: MM_vector, eTq_vector, insw_vector



!  
!
! !AUTHORS
!   Original version by P. Ruardij and M. Vichi
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
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES)  :: fN4N3n
  real(RLEN),dimension(NO_BOXES)  :: fN6O2r
  real(RLEN),dimension(NO_BOXES)  :: eo
  real(RLEN),dimension(NO_BOXES)  :: er
  real(RLEN),dimension(NO_BOXES)  :: osat
  real(RLEN),dimension(NO_BOXES)  :: rPAo
  real(RLEN),dimension(NO_BOXES)  :: fR6N5s
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !
  ! Regulating factors
  !

  eo  =   MM_vector(  O2o(:),  p_clO2o)
  er  =   MM_vector(  N6r(:),  p_clN6r)

  !
  ! Nitrification in the water
  !

  fN4N3n  =   p_sN4N3* N4n(:)* eTq_vector(  ETW(:),  p_q10N4N3)* eo
  call flux_vector( iiPel, ppN4n,ppN3n, fN4N3n )

  call flux_vector( iiPel, ppO2o,ppO2o,-( fN4N3n* p_qon_nitri) )

  !
  ! Denitrification in the water
  !

  rPAo  =   flPTN6r(:)/ p_qro
  flN3O4n(:) = p_sN3O4n* eTq_vector( ETW(:), p_q10N4N3)* er* rPAo/ p_rPAo* &
    N3n(:)

  call flux_vector( iiPel, ppN3n,ppO4n, flN3O4n(:) )
  call flux_vector( iiPel, ppN6r,ppN6r,-( p_qro* flN3O4n(:)* p_qon_dentri* &
    insw_vector( -( O2o(:)- N6r(:)/ p_qro))) )


  !
  ! Reoxidation of reduction equivalents
  !

  fN6O2r  =   p_rOS* N6r(:)* eo

  call flux_vector( iiPel, ppN6r,ppN6r,-( fN6O2r) )
  call flux_vector( iiPel, ppO2o,ppO2o,-( fN6O2r/ p_qro) )

  !
  ! Regeneration of dissolved silica
  !

  fR6N5s  =   p_sR6N5* eTq_vector(  ETW(:),  p_q10R6N5)* R6s(:)
  call flux_vector( iiPel, ppR6s,ppN5s, fR6N5s )

  !
  !
  !

  call flux_vector( iiPel, ppP1s,ppR6s, flP1R6s(:) )




  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
