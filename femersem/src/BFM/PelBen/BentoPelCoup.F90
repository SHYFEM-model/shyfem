#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BentoPelCoup
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BentoPelCoupDynamics
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: R6c, R6n, R6p, R6s, &
  ! O2o, N1p, N3n, N4n, N5s, N6r, R1c, R1n, R1p
  ! For the following Pelagic-group-states fluxes are defined: PhytoPlankton
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY
  ! The following Pelagic 1-d global boxvars  are used: Depth
  ! The following Benthic 1-d global boxvars are used: jPIY3c, &
  ! PIc, jRIY3c, jRIY3n, jRIY3p, jRIY3s, jG2O2o, jK1N1p, jK3N3n, jK4N4n, &
  ! jK5N5s, jK6N6r, rutQ6c, rutQ6n, rutQ6p, rutQ6s, rutQ1c, rutQ1n, rutQ1p
  ! The following Benthic 2-d global boxvars are modified : retPIc, retPIn, &
  ! retPIp, retPIl, retPIs
  ! The following groupmember vars  are used: iiPhytoPlankton, iiP1
  ! The following constituent constants  are used: iiC, iiN, iiP, iiL, iiS
  ! The following 0-d global box parametes are used: &
  ! CalcBenthicFlag, AssignPelBenFluxesInBFMFlag
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem, ONLY: R6c, R6n, R6p, R6s, O2o, N1p, N3n, N4n, N5s, N6r, &
    R1c, R1n, R1p, PhytoPlankton, D2STATE
  use mem, ONLY: ppR6c, ppR6n, ppR6p, ppR6s, ppO2o, ppN1p, &
    ppN3n, ppN4n, ppN5s, ppN6r, ppR1c, ppR1n, ppR1p, ppPhytoPlankton, BoxNumberZ, &
    NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, BoxNumber, &
    BoxNumberXY, Depth, jPIY3c, PIc, jRIY3c, jRIY3n, jRIY3p, jRIY3s, jG2O2o, &
    jK1N1p, jK3N3n, jK4N4n, jK5N5s, jK6N6r, rutQ6c, rutQ6n, rutQ6p, rutQ6s, &
    rutQ1c, rutQ1n, rutQ1p, retPIc, retPIn, retPIp, retPIl, retPIs, &
    iiPhytoPlankton, iiP1, iiC, iiN, iiP, iiL, iiS, iiBen, iiPel, flux
  use mem_Param,  ONLY: CalcBenthicFlag, AssignPelBenFluxesInBFMFlag



!  
!
! !AUTHORS
!   Piet Ruardij  
!
!
!
! !REVISION_HISTORY
!   Created at Mon Apr 19 00:08:12 CEST 2004
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
  ! Local Vectors used  of group vectors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN), dimension(:), pointer  ::lcl_PhytoPlankton
  real(RLEN), dimension(:), pointer  ::lcl_ppPhytoPlankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  real(RLEN)  :: uptake
  real(RLEN)  :: Pc

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


      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate Phyto Fluxes to Filterfeeder from Pelagic for
      ! all phyt types/constituents
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      if ( CalcBenthicFlag> 0) then
        if ( jPIY3c(BoxNumberXY)> 0) then
          do i = 1 , ( iiPhytoPlankton)

            lcl_PhytoPlankton => PhytoPlankton(i,iiC)
            Pc  =   lcl_PhytoPlankton(BoxNumber)
            if ( Pc> 0.0D+00) then
              uptake  =   jPIY3c(BoxNumberXY)* Pc/ PIc(BoxNumberXY)
              retPIc(i,BoxNumberXY)  =   uptake
              lcl_PhytoPlankton => PhytoPlankton(i,iiN)
              retPIn(i,BoxNumberXY) = uptake* lcl_PhytoPlankton(BoxNumber)/ &
                Pc
              lcl_PhytoPlankton => PhytoPlankton(i,iiP)
              retPIp(i,BoxNumberXY) = uptake* lcl_PhytoPlankton(BoxNumber)/ &
                Pc
              lcl_PhytoPlankton => PhytoPlankton(i,iiL)
              retPIl(i,BoxNumberXY) = uptake* lcl_PhytoPlankton(BoxNumber)/ &
                Pc
              if ( i== iiP1) then
                lcl_PhytoPlankton => PhytoPlankton(i,iiS)
                retPIs(i,BoxNumberXY) = uptake* lcl_PhytoPlankton(BoxNumber)/ &
                  Pc
              end if

            end if

          end do

        end if


        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Net detritus Fluxes to Benthic from Pelagic by Y3
        ! net flux= uptake - excretion of food : flux may ber negative!
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call flux(BoxNumber, iiPel, ppR6c, ppR6c, -( jRIY3c(BoxNumberXY)/ &
          Depth(BoxNumber)) )
        call flux(BoxNumber, iiPel, ppR6n, ppR6n, -( jRIY3n(BoxNumberXY)/ &
          Depth(BoxNumber)) )
        call flux(BoxNumber, iiPel, ppR6p, ppR6p, -( jRIY3p(BoxNumberXY)/ &
          Depth(BoxNumber)) )
        call flux(BoxNumber, iiPel, ppR6s, ppR6s, -( jRIY3s(BoxNumberXY)/ &
          Depth(BoxNumber)) )


        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! All Fluxes to Benthic from Pelagic defined for the
        ! Pelagic State variables
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        if ( AssignPelBenFluxesInBFMFlag) then

          ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Nutrient Fluxes to Benthic from Pelagic
          ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          call flux(BoxNumber, iiPel, ppO2o, ppO2o, jG2O2o(BoxNumberXY)/ &
            Depth(BoxNumber) )
          call flux(BoxNumber, iiPel, ppN1p, ppN1p, jK1N1p(BoxNumberXY)/ &
            Depth(BoxNumber) )
          call flux(BoxNumber, iiPel, ppN3n, ppN3n, jK3N3n(BoxNumberXY)/ &
            Depth(BoxNumber) )
          call flux(BoxNumber, iiPel, ppN4n, ppN4n, jK4N4n(BoxNumberXY)/ &
            Depth(BoxNumber) )
          call flux(BoxNumber, iiPel, ppN5s, ppN5s, jK5N5s(BoxNumberXY)/ &
            Depth(BoxNumber) )
          call flux(BoxNumber, iiPel, ppN6r, ppN6r, jK6N6r(BoxNumberXY)/ &
            Depth(BoxNumber) )

          ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! PhytoPlankton Fluxes to Benthic from Pelagic by Y3
          ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          if ( jPIY3c(BoxNumberXY)> 0.0D+00) then
            do i = 1 , ( iiPhytoPlankton)

              if ( Pc> 0.0D+00) then
                call flux(BoxNumber, iiPel, ppPhytoPlankton(i,iiC), &
                  ppPhytoPlankton(i,iiC), -( retPIc(i,BoxNumberXY)/ &
                  Depth(BoxNumber)) )
                call flux(BoxNumber, iiPel, ppPhytoPlankton(i,iiN), &
                  ppPhytoPlankton(i,iiN), -( retPIn(i,BoxNumberXY)/ &
                  Depth(BoxNumber)) )
                call flux(BoxNumber, iiPel, ppPhytoPlankton(i,iiP), &
                  ppPhytoPlankton(i,iiP), -( retPIp(i,BoxNumberXY)/ &
                  Depth(BoxNumber)) )
                call flux(BoxNumber, iiPel, ppPhytoPlankton(i,iiL), &
                  ppPhytoPlankton(i,iiL), -( retPIl(i,BoxNumberXY)/ &
                  Depth(BoxNumber)) )
                if ( i== iiP1) then
                  !No Y3.s defined, all silicate uptaken is moved into sink
                  call flux(BoxNumber, iiPel, ppPhytoPlankton(i,iiS), &
                    ppPhytoPlankton(i,iiS), -( retPIs(i,BoxNumberXY)/ &
                    Depth(BoxNumber)) )
                end if

              end if

            end do

          end if




          ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Total Sedimentation flux to Benthic from Pelagic defined
          ! for the pealgic state variables.
          ! (See sedimentation for definition of the fluxes of benthic &
          ! variables)
          ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          call flux(BoxNumber, iiPel, ppR6c, ppR6c, -( rutQ6c(BoxNumberXY)/ &
            Depth(BoxNumber)) )
          call flux(BoxNumber, iiPel, ppR6n, ppR6n, -( rutQ6n(BoxNumberXY)/ &
            Depth(BoxNumber)) )
          call flux(BoxNumber, iiPel, ppR6p, ppR6p, -( rutQ6p(BoxNumberXY)/ &
            Depth(BoxNumber)) )
          call flux(BoxNumber, iiPel, ppR6s, ppR6s, -( rutQ6s(BoxNumberXY)/ &
            Depth(BoxNumber)) )

          call flux(BoxNumber, iiPel, ppR1c, ppR1c, -( rutQ1c(BoxNumberXY)/ &
            Depth(BoxNumber)) )
          call flux(BoxNumber, iiPel, ppR1n, ppR1n, -( rutQ1n(BoxNumberXY)/ &
            Depth(BoxNumber)) )
          call flux(BoxNumber, iiPel, ppR1p, ppR1p, -( rutQ1p(BoxNumberXY)/ &
            Depth(BoxNumber)) )
        end if


      end if


    end DO


  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
