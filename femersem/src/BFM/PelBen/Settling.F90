#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Settling
!
! DESCRIPTION
!   This process describes the dynamics of sedimentation and
!    deposition of phytoplankton (P1, P2, P3, P4, ...) and detritus (R6)
!    in the benthic system.
!    A burial velocity is defined, which controls the magnitude
!    of the inflow rate of detritus from the pelagic form R6 to
!    the benthic form Q6.
!    The processes described here taks only place in the lowest
!    boxes in the Z-direction.
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine SettlingDynamics
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: R1c, R6c, R1n, R6n, &
  ! R1p, R6p, R6s
  ! For the following Pelagic-group-states fluxes are defined: PhytoPlankton
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY
  ! The following Pelagic 1-d global boxvars  are used: Depth
  ! The following Benthic 1-d global boxvars are modified : rutQ6c, &
  ! rutQ6n, rutQ6p, rutQ6s, rutQ1c, rutQ1n, rutQ1p
  ! The following Pelagic 2-d global boxvars  are used: sediPI
  ! The following groupmember vars  are used: iiPhytoPlankton, iiP1
  ! The following constituent constants  are used: iiC, iiN, iiP, iiS
  ! The following 0-d global box parametes are used: CalcBenthicFlag, &
  ! p_pe_R1c, p_pe_R1n, p_pe_R1p
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: R1c, R6c, R1n, R6n, R1p, R6p, R6s, PhytoPlankton, D3STATE
  use mem, ONLY: ppR1c, ppR6c, ppR1n, ppR6n, ppR1p, ppR6p, &
    ppR6s, ppPhytoPlankton, BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, &
    BoxNumberY, NO_BOXES_Y, BoxNumber, BoxNumberXY, Depth, rutQ6c, rutQ6n, rutQ6p, &
    rutQ6s, rutQ1c, rutQ1n, rutQ1p, sediPI, iiPhytoPlankton, iiP1, iiC, iiN, iiP, &
    iiS, iiBen, iiPel, flux
  use mem_Param,  ONLY: CalcBenthicFlag, p_pe_R1c, p_pe_R1n, p_pe_R1p
  use mem_Settling



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
  ! Local Vectors used  of group vectors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN), dimension(:), pointer  ::lcl_PhytoPlankton
  real(RLEN), dimension(:), pointer  ::lcl_ppPhytoPlankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  integer  :: seq
  real(RLEN)  :: sedi
  real(RLEN)  :: psed
  real(RLEN)  :: ruQIc
  real(RLEN)  :: ruQIn
  real(RLEN)  :: ruQIp
  real(RLEN)  :: ruQIs
  real(RLEN)  :: ruQ1c
  real(RLEN)  :: ruQ1n
  real(RLEN)  :: ruQ1p
  real(RLEN)  :: ruQ1s
  real(RLEN)  :: ruQ6c
  real(RLEN)  :: ruQ6n
  real(RLEN)  :: ruQ6p
  real(RLEN)  :: ruQ6s
  real(RLEN)  :: ruQ6m

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! user defined external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, external  :: D3toD1
  integer, external  :: D2toD1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  BoxNumberZ = 1
  DO BoxNumberY=1,NO_BOXES_Y
    DO BoxNumberX=1,NO_BOXES_X
      BoxNumber=D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ)
      BoxNumberXY=D2toD1(BoxNumberX,BoxNumberY)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Starting real settling:
      ! phytoplankton settling are tmeporarly stalled as R1 and R6 and
      ! subsequently exported (BentoPelCoup) adn sedimentated in bnethos &
      ! (sedimentation)
      ! All fluxes are seperated in sink and sources fluxes:
      ! to avoid problems with the definitions of Pel. fluxes .
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      rutQ6c(BoxNumberXY)  =   0.0D+00
      rutQ6n(BoxNumberXY)  =   0.0D+00
      rutQ6p(BoxNumberXY)  =   0.0D+00
      rutQ6s(BoxNumberXY)  =   0.0D+00

      rutQ1c(BoxNumberXY)  =   0.0D+00
      rutQ1n(BoxNumberXY)  =   0.0D+00
      rutQ1p(BoxNumberXY)  =   0.0D+00


      if ( CalcBenthicFlag> 0) then
        do i = 1 , ( iiPhytoPlankton)

          sedi  =   sediPI(i,BoxNumber)
          if ( sedi> 0.0D+00) then


            lcl_PhytoPlankton => PhytoPlankton(i,iiC)
            ruQIc  =   sedi* lcl_PhytoPlankton(BoxNumber)
            ruQ1c  =   p_pe_R1c* ruQIc
            ruQ6c  =   ruQIc- ruQ1c
            call flux(BoxNumber, iiPel, ppPhytoPlankton(i,iiC), &
              ppPhytoPlankton(i,iiC), -( ruQIc/ Depth(BoxNumber)) )
            call flux(BoxNumber, iiPel, ppR1c, ppR1c, ruQ1c/ Depth(BoxNumber) )
            call flux(BoxNumber, iiPel, ppR6c, ppR6c, ruQ6c/ Depth(BoxNumber) )
            rutQ1c(BoxNumberXY)  =   rutQ1c(BoxNumberXY)+ ruQ1c
            rutQ6c(BoxNumberXY)  =   rutQ6c(BoxNumberXY)+ ruQ6c

            lcl_PhytoPlankton => PhytoPlankton(i,iiN)
            ruQIn  =   sedi* lcl_PhytoPlankton(BoxNumber)
            ruQ1n  =   p_pe_R1n* ruQIn
            ruQ6n  =   ruQIn- ruQ1n
            call flux(BoxNumber, iiPel, ppPhytoPlankton(i,iiN), &
              ppPhytoPlankton(i,iiN), -( ruQIn/ Depth(BoxNumber)) )
            call flux(BoxNumber, iiPel, ppR1n, ppR1n, ruQ1n/ Depth(BoxNumber) )
            call flux(BoxNumber, iiPel, ppR6n, ppR6n, ruQ6n/ Depth(BoxNumber) )
            rutQ1n(BoxNumberXY)  =   rutQ1n(BoxNumberXY)+ ruQ1n
            rutQ6n(BoxNumberXY)  =   rutQ6n(BoxNumberXY)+ ruQ6n

            lcl_PhytoPlankton => PhytoPlankton(i,iiP)
            ruQIp  =   sedi* lcl_PhytoPlankton(BoxNumber)
            ruQ1p  =   p_pe_R1p* ruQIp
            ruQ6p  =   ruQIp- ruQ1p
            call flux(BoxNumber, iiPel, ppPhytoPlankton(i,iiP), &
              ppPhytoPlankton(i,iiP), -( ruQIp/ Depth(BoxNumber)) )
            call flux(BoxNumber, iiPel, ppR1p, ppR1p, ruQ1p/ Depth(BoxNumber) )
            call flux(BoxNumber, iiPel, ppR6p, ppR6p, ruQ6p/ Depth(BoxNumber) )
            rutQ1p(BoxNumberXY)  =   rutQ1p(BoxNumberXY)+ ruQ1p
            rutQ6p(BoxNumberXY)  =   rutQ6p(BoxNumberXY)+ ruQ6p

            if ( i== iiP1) then
              lcl_PhytoPlankton => PhytoPlankton(i,iiS)
              ruQ6s  =   sedi* lcl_PhytoPlankton(BoxNumber)
              call flux(BoxNumber, iiPel, ppPhytoPlankton(i,iiS), &
                ppPhytoPlankton(i,iiS), -( ruQ6s/ Depth(BoxNumber)) )
              call flux(BoxNumber, iiPel, ppR6s, ppR6s, ruQ6s/ Depth(BoxNumber) &
                )
              rutQ6s(BoxNumberXY)  =   rutQ6s(BoxNumberXY)+ ruQ6s
            end if

          end if


        end do


        ! R6 into Q6:

        ruQ6c  =   p_burvel* R6c(BoxNumber)
        ruQ6n  =   p_burvel* R6n(BoxNumber)
        ruQ6p  =   p_burvel* R6p(BoxNumber)
        ruQ6s  =   p_burvel* R6s(BoxNumber)

        rutQ6c(BoxNumberXY)  =   rutQ6c(BoxNumberXY)+ ruQ6c
        rutQ6n(BoxNumberXY)  =   rutQ6n(BoxNumberXY)+ ruQ6n
        rutQ6p(BoxNumberXY)  =   rutQ6p(BoxNumberXY)+ ruQ6p
        rutQ6s(BoxNumberXY)  =   rutQ6s(BoxNumberXY)+ ruQ6s

      end if


    end DO


  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
