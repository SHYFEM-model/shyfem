!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   CalculateFromSet.f90
!
! FILE
!   CalculateFromSet.f90
!
! DESCRIPTION
!       NUTR: which series of results are used (sequence number)
!       mode: which kind of action is made:
!       1/INTEGRAL= calculate integral between from and to
!       0 or EQUATION= calculate equation at from
!       -1 or DERIVATIVE /-2 or SDERIVATIVE= calculate first/second derivative at from
!       from    : input (see above)
!       to      :input (see above) 
!
!   This file is generated from f77 code, using a code generator which
!   transposes from the F77 code into F90
!   the code USES module file as used in the BFM model
!   F90 code generator written by P. Ruardij.
!
! AUTHORS
!   Piet Ruardij ,NIOZ
!
! CHANGE_LOG
!
!
! COPYING
!
!   Copyright (C) 2004 P. Ruardij, the mfstep group, the ERSEM team
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
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!
      REAL(RLEN) FUNCTION CalculateFromSet(NUTR,mode,input,from,to)
        USE global_mem, ONLY:RLEN
        USE constants
        USE bennut_interface,ONLY: noutput,CalculateLayer,CalculateShift
        IMPLICIT  NONE
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::nutr ! Specification
        integer,intent(IN) ::input ! Specification
        REAL(RLEN),intent(IN) ::from ! Specification
        REAL(RLEN),intent(IN) ::to ! Specification
        !mode =1 : perfrom an integration.
        !if x_from>x_to an integration is perform between x_to and x_from but the sign is &
        ! changed
        !if x_to<x_from normal integration
        !mode =0 : calculate the equation
        !mode=-1 : perform an differentiation.

        REAL(RLEN) ::x_from
        REAL(RLEN) ::x_to
        REAL(RLEN) ::xo
        REAL(RLEN) ::xb
        REAL(RLEN) ::r
        REAL(RLEN) ::s
        integer ::j

        select case (mode)
          case (SHIFT)
            CalculateFromSet = CalculateShift(NUTR,input,from,to);
          case (INTEGRAL,EXPONENTIAL_INTEGRAL)
             s=to-from
             if (abs(s)/(from+to+.5D-6) < 0.5D-6) then
               !if difference between x and from is less than machine precision &
               ! answer will be 0.
               CalculateFromSet=0.D+00
               return
             elseif ( s < 0.0D+00 ) then
               !if integral is asked for negative interval, border of intervales &
               ! are exchanged and answer will be made negative:
                 x_to=from
                 x_from=to
             else 
                 x_from=from
                 x_to=to
             endif
             !calculation of the integral...
             !find first layer get upper border of the first layer.
             call CalculateLayer(NUTR,0,x_from,j,xb)
             !make x0 equatian to under border....
             xo=x_from
             r=0.D+00
             do while (x_to > xb)
               !calculate integration (mode=1) of equation for nutrient NUTR,
               !equation nr j between the boundaries xo and xb
               r=r+noutput(NUTR,j,mode,input,xo,xb)
               !make under border of the next equation equal to the upper border
               !of the actual equation
               xo=xb
               !goto the next equation
               !get upper border of the next equation/layer.
               call CalculateLayer(NUTR,j+1,xo,j,xb)
             enddo
             r=r+noutput(NUTR,j,mode,input,xo,x_to)
             CalculateFromSet=sign(r,s)
          !calculation the equation or 1 derivative or 2 derivative
          case (EQUATION,DERIVATIVE,SDERIVATIVE)
             !find layer/equation
             call CalculateLayer(NUTR,0,from,j,xb)
             CalculateFromSet=noutput(NUTR,j,mode,input,from,0.0D+00)
        end select
        return
      end

