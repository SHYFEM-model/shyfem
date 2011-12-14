#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ControlBenPartNutrientBuffers
!
! DESCRIPTION
!   This process forces the system to calculate a flux from
!       dissolved pelagic nutrient to the particualte nutreit pool
!       in the sediment. The reason for this is that if one of the
!       benthic pools become zero the partial differential equations 
!       used to calclate the diagenict processes and vertical fluxes
!       are not valid any more.
!       Therfor this process maximize input to benthos in
!       cases when nearly no particulate detritus enters the benthic system.
!       This is done by transferring dissolved nutrients N3.n,N4.n,N5.s
!       in the layer adjacent to the sediment into
!       benthic particulate nutrients (Q6.n,Q6.p,Q6.s).
!       A fixed fraction per day will be transferred to the benthos :
!        a. Minimized to  a fixed maximum value
!        b. Only when the settling rate is lower than the fixed maximum value
!        c. Only when het particulate nutrient buffer is really low.
!           This is determined to check if the quotient of the nutrient buffer
!           and the  fixed fraction is below a parameter, standing for the
!           specific refill rate.
!       The values of fixed fraction and the maximum value for the 
!       different nutrients are calculated according to Redfield
!       The user decides which nutrients should be handled
!          
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine ControlBenPartNutrientBuffersDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: Q6n, D9m, D7m, Q6p, &
  ! D8m, Q6s
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY
  ! The following Benthic 1-d global boxvars are modified : jK4N4n, jK1N1p, &
  ! jK5N5s
  ! The following Benthic 1-d global boxvars are used: N4n_Ben, &
  ! Depth_Ben, rutQ6n, N1p_Ben, rutQ6p, N5s_Ben, rutQ6s
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: Q6n, D9m, D7m, Q6p, D8m, Q6s, D2STATE
  use mem, ONLY: ppQ6n, ppD9m, ppD7m, ppQ6p, ppD8m, ppQ6s, &
    BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
    BoxNumber, BoxNumberXY, jK4N4n, jK1N1p, jK5N5s, N4n_Ben, Depth_Ben, rutQ6n, &
    N1p_Ben, rutQ6p, N5s_Ben, rutQ6s, iiBen, iiPel, flux
  use mem_ControlBenPartNutrientBuffers



!  
!
! !AUTHORS
!   Piet Ruardij
!
!
!
! !REVISION_HISTORY
!   Created at Fri Jan 27 08:47:37 CET 2006
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
  real(RLEN)  :: r
  real(RLEN)  :: rate

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


      !----N:----------------------------------------------------------------
      if ( p_control_n) then
        rate = min( p_rn, p_s* N4n_Ben(BoxNumberXY)* &
          Depth_Ben(BoxNumberXY))
        if ( rutQ6n(BoxNumberXY)< p_rn .AND. Q6n(BoxNumberXY)/ rate< &
          p_refill_time) then
          rate  =   max(  0.0D+00,  rate- rutQ6n(BoxNumberXY))
          call flux(BoxNumberXY, iiBen, ppQ6n, ppQ6n, -(- rate) )
          call flux(BoxNumberXY, iiBen, ppD9m, ppD9m, ( 0.0D+00- &
            D7m(BoxNumberXY))* rate/( 1.D-80+ Q6n(BoxNumberXY)) )
          jK4N4n(BoxNumberXY)  =   jK4N4n(BoxNumberXY)- rate
        end if

      end if


      !----P:----------------------------------------------------------------
      if ( p_control_p) then
        r  =   p_rn/ 16.0D+00  ! Correct according Redfield.
        rate  =   min(  r,  p_s* N1p_Ben(BoxNumberXY)* Depth_Ben(BoxNumberXY))
        if ( rutQ6p(BoxNumberXY)< r .AND. Q6p(BoxNumberXY)/ rate< p_refill_time) &
          then
          rate  =   max(  0.0D+00,  rate- rutQ6p(BoxNumberXY))
          call flux(BoxNumberXY, iiBen, ppQ6p, ppQ6p, -(- rate) )
          call flux(BoxNumberXY, iiBen, ppD9m, ppD9m, ( 0.0D+00- &
            D8m(BoxNumberXY))* rate/( 1.D-80+ Q6p(BoxNumberXY)) )
          jK1N1p(BoxNumberXY)  =   jK1N1p(BoxNumberXY)- rate
        end if

      end if


      !----Si:---------------------------------------------------------------
      if ( p_control_s) then
        r  =   p_rn* 22.0D+00/ 16.0D+00  ! Correct according Readfield.
        rate  =   min(  r,  p_s* N5s_Ben(BoxNumberXY)* Depth_Ben(BoxNumberXY))
        if ( rutQ6s(BoxNumberXY)< r .AND. Q6s(BoxNumberXY)/ rate< p_refill_time) &
          then
          rate  =   max(  0.0D+00,  rate- rutQ6s(BoxNumberXY))
          call flux(BoxNumberXY, iiBen, ppQ6s, ppQ6s, -(- rate) )
          call flux(BoxNumberXY, iiBen, ppD9m, ppD9m, ( 0.0D+00- &
            D9m(BoxNumberXY))* rate/( 1.D-80+ Q6s(BoxNumberXY)) )
          jK5N5s(BoxNumberXY)  =   jK5N5s(BoxNumberXY)- rate
        end if

      end if




    end DO



  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
