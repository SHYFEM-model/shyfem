#include "DEBUG.h"
#include "INCLUDE.h"

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcVerticalExtinction
!
! DESCRIPTION
!   Calculates the vertical extinction.
!     
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  SUBROUTINE CalcVerticalExtinction()
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): R6c
  ! The following box states are used (NOT in fluxes): PhytoPlankton
  ! The following Pelagic 1-d global boxvars are modified : xEPS
  ! The following Pelagic 1-d global boxvars  are used: ABIO_eps, ESS
  ! The following groupmember vars  are used: iiPhytoPlankton
  ! The following constituent constants  are used: iiC, iiL
  ! The following 0-d global box parametes are used: p_eps0, &
  ! p_epsR6, p_epsESS, ChlLightFlag, p_epsChla
  ! The following 1-d global parameter vars are used: p_qchlc
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem,  ONLY: R6c, PhytoPlankton
#endif
  use mem, ONLY: ppR6c, ppPhytoPlankton, xEPS, ABIO_eps, ESS, iiPhytoPlankton, &
    iiC, iiL, NO_BOXES, iiBen, iiPel, flux_vector
  use mem_Param, ONLY: p_eps0, p_epsR6, p_epsESS, ChlLightFlag, p_epsChla, &
    p_qchlc



!  
!
! !AUTHORS
!   ERSEM-team
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
  integer  :: i


  select case ( p_eps0== 0.0D+00)

    case( .TRUE. )
      xEPS(:)  =   ABIO_eps(:)+ p_epsR6* R6c(:)



    case( .FALSE. )
      xEPS(:)  =   p_eps0+ p_epsESS* ESS(:)+ p_epsR6* R6c(:)


  end select


  select case ( ChlLightFlag)

    case ( 1 )
      do i = 1 , ( iiPhytoPlankton)

        xEPS(:)  =   xEPS(:)+ p_epsChla* p_qchlc( i)* PhytoPlankton(i,iiC)
      end do




    case ( 2 )
      do i = 1 , ( iiPhytoPlankton)

        xEPS(:)  =   xEPS(:)+ p_epsChla* PhytoPlankton(i,iiL)
	
      end do



  end select
 




  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
