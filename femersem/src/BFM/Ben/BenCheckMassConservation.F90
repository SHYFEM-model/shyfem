#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenCheckMassConservation
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenCheckMassConservationDynamics
!
! !USES:
  ! The following Benthic-states are used (NOT in fluxes): Q1p, Q6p, Q1n, Q6n, &
  ! Q6s, Y1p, Y2p, Y3p, Y4p, Y5p, H1p, H2p, Q11p, K1p, K11p, Y1n, Y2n, Y3n, &
  ! Y4n, Y5n, H1n, H2n, Q11n, K4n, K14n, K21p, G4n, K3n, K24n, K5s
  ! The following Benthic 1-d global boxvars got a value: totbenp, totbenn, &
  ! totbens
  ! The following 0-d global box parametes are used: CalcBenthicFlag
  ! The following global constants are used: RLEN
  ! The following constants are used: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem, ONLY: Q1p, Q6p, Q1n, Q6n, Q6s, Y1p, Y2p, Y3p, Y4p, Y5p, H1p, H2p, &
    Q11p, K1p, K11p, Y1n, Y2n, Y3n, Y4n, Y5n, H1n, H2n, Q11n, K4n, K14n, K21p, &
    G4n, K3n, K24n, K5s
#endif
  use mem, ONLY: ppQ1p, ppQ6p, ppQ1n, ppQ6n, ppQ6s, ppY1p, ppY2p, &
    ppY3p, ppY4p, ppY5p, ppH1p, ppH2p, ppQ11p, ppK1p, ppK11p, ppY1n, ppY2n, &
    ppY3n, ppY4n, ppY5n, ppH1n, ppH2n, ppQ11n, ppK4n, ppK14n, ppK21p, &
    ppG4n, ppK3n, ppK24n, ppK5s, totbenp, totbenn, totbens, NO_BOXES_XY, iiBen, &
    iiPel, flux_vector
  use constants,  ONLY: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL
  use mem_Param,  ONLY: CalcBenthicFlag



!  
!
! !AUTHORS
!   Piet Ruardij
!
!
!
! !REVISION_HISTORY
!   Created at Mon Nov 21 09:44:23 CET 2005
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
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  select case ( CalcBenthicFlag)

    case ( 0 )
      totbenp(:)  =   0.0D+00
      totbenn(:)  =   0.0D+00
      totbens(:)  =   0.0D+00



    case ( BENTHIC_RETURN )  ! Simple benthic return
      ! Mass conservation variables
      totbenp(:)  =  ( Q1p(:)+ Q6p(:))
      totbenn(:)  =  ( Q1n(:)+ Q6n(:))
      totbens(:)  =  ( Q6s(:))



    case ( BENTHIC_BIO )  ! Intermediate benthic return
      ! Mass conservation variables
      totbenp(:) = ( Y1p(:)+ Y2p(:)+ Y3p(:)+ Y4p(:)+ Y5p(:)+ H1p(:)+ &
        H2p(:)+ Q1p(:)+ Q6p(:)+ Q11p(:)+ K1p(:)+ K11p(:))
      totbenn(:) = ( Y1n(:)+ Y2n(:)+ Y3n(:)+ Y4n(:)+ Y5n(:)+ H1n(:)+ &
        H2n(:)+ Q1n(:)+ Q6n(:)+ Q11n(:)+ K4n(:)+ K14n(:))
      totbens(:)  =  ( Q6s(:))



    case ( BENTHIC_FULL )  ! Full benthic nutrients
      totbenp(:) = ( Y1p(:)+ Y2p(:)+ Y3p(:)+ Y4p(:)+ Y5p(:)+ H1p(:)+ &
        H2p(:)+ Q1p(:)+ Q6p(:)+ Q11p(:)+ K1p(:)+ K11p(:)+ K21p(:))
      totbenn(:) = ( Y1n(:)+ Y2n(:)+ Y3n(:)+ Y4n(:)+ Y5n(:)+ H1n(:)+ &
        H2n(:)+ Q1n(:)+ Q6n(:)+ Q11n(:)+ G4n(:)+ K3n(:)+ K4n(:)+ K14n(:)+ K24n(:))
      totbens(:)  =   K5s(:)+ Q6s(:)


  end select






  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
