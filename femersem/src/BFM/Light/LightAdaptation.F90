#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: LightAdaptation
!
! DESCRIPTION
!   This routine describes the photoadaptation of the phytoplankton to 
!	the prevailing irradiance level at depth
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine LightAdaptationDynamics(phyto, ppphytoc, ppphyton, ppphytop, &
    ppphytos, ppphytol)
!
! !USES:
  ! For the following Pelagic-group-states fluxes are defined: PhytoPlankton
  ! The following Pelagic 1-d global vars are used: D3STATETYPE
  ! The following Pelagic 1-d global boxvars  are used: Depth, xEPS, EIR
  ! The following Pelagic 2-d global boxvars  are used: EPLi
  ! The following constituent constants  are used: iiL, iiC
  ! The following global constants are used: RLEN,NOTRANSPORT

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,NOTRANSPORT
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem,  ONLY: D3STATE, PhytoPlankton
#endif
  use mem, ONLY: ppPhytoPlankton, D3STATETYPE, Depth, xEPS, EIR, EPLi, &
    iiL, iiC, Source_D3_vector, NO_BOXES, iiBen, iiPel, flux_vector
  use mem_LightAdaptation


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: phyto
  integer,intent(IN) :: ppphytoc
  integer,intent(IN) :: ppphyton
  integer,intent(IN) :: ppphytop
  integer,intent(IN) :: ppphytos
  integer,intent(IN) :: ppphytol

!  
!
! !AUTHORS
!   Original version by W. Ebenhoeh, Oldenburg University
!                           Hanneke Baretta-Bekker, VKI
!       Translated to OpenSesame by Piet Ruardij
!	Phytoplankton species dependency added by M. Vichi, INGV
!
!
!
! !REVISION_HISTORY
!   File created on 8 feb. 1997
!	Modified by Daji and JWB, 19/6/1998
!	Checked by D.Mills and JWB 030429
!	Horrible error removed in addepth
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
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES) :: phytoc
  real(RLEN),dimension(NO_BOXES) :: phyton
  real(RLEN),dimension(NO_BOXES) :: phytop
  real(RLEN),dimension(NO_BOXES) :: phytos
  real(RLEN),dimension(NO_BOXES) :: phytol
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  real(RLEN),dimension(NO_BOXES)  :: addepth
  real(RLEN),dimension(NO_BOXES)  :: adfactor
  real(RLEN),dimension(NO_BOXES)  :: sum_state_phyto
  real(RLEN),dimension(NO_BOXES)  :: sum_rate_phyto
  real(RLEN),dimension(NO_BOXES)  :: loss
  real(RLEN),dimension(NO_BOXES)  :: rate_PLi
  real(RLEN),dimension(NO_BOXES)  :: rate_EPLi
  real(RLEN),dimension(NO_BOXES)  :: new_EPLi
  real(RLEN),dimension(NO_BOXES)  :: eir_c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  phytoc = D3STATE(ppphytoc,:)
  phyton = D3STATE(ppphyton,:)
  phytop = D3STATE(ppphytop,:)
  phytos = D3STATE(ppphytos,:)
  phytol = D3STATE(ppphytol,:)


  ! EPLi[%phyto] has been already calculated in Irradiance_PI.p

  ! computation of the change in IRR_OPT
  ! adaptation to the light in the depth p_addepth below surface.
  ! p_addepth is superseded if this function is used in a
  ! high resolution vertical grid. In that case the central depth
  ! of the layer is used.


  ! addepth  = min(Depth/2.0, p_addepth);
  addepth  =   min(  Depth(:),  p_addepth(phyto))
  adfactor  =   exp( - xEPS(:)* addepth)

  eir_c  =   EIR(:)* adfactor

  select case ( p_isw(phyto))

    case ( 1 )
      new_EPLi  =   max(  eir_c,  p_clEPLi(phyto))
      new_EPLi  =   min(  new_EPLi,  p_chEPLi(phyto))



    case ( 2 )
      new_EPLi = max( 2.0D+00* eir_c* p_chEPLi(phyto)/( &
        eir_c+ p_chEPLi(phyto)), p_clEPLi(phyto))
      new_EPLi  =   min(  new_EPLi,  p_chEPLi(phyto))


  end select


  ! Speed of adaptation is controlled by p_ruPLi ( 1 maximum speed
  !                        0 no adaptation )

  i  =   ppPhytoPlankton(phyto,iiL)
  select case ( D3STATETYPE( i))

    case ( NOTRANSPORT )

      rate_EPLi  =   p_ruEPLi(phyto)*( new_EPLi- EPLi(phyto,:))
      call flux_vector( iiPel, &
        ppPhytoPlankton(phyto,iiL),ppPhytoPlankton(phyto,iiL), rate_EPLi )




    case default

      rate_PLi = Source_D3_vector(ppPhytoPlankton(phyto,iiC))* EPLi(phyto,:)+ &
        p_ruEPLi(phyto)*( new_EPLi- EPLi(phyto,:))* phytoc
      call flux_vector( iiPel, &
        ppPhytoPlankton(phyto,iiL),ppPhytoPlankton(phyto,iiL), rate_PLi )




  end select






  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
