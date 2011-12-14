!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   noutput.f90
!
! FILE
!   noutput.f90
!
! DESCRIPTION
!   !

!NUTR: kind of nutrient
!mode: 0 : copy result from ninput to ???s(??,NUTR) vars
!mode: -1<<1000 : fill bborder
!mode: 0<<1000 : calculate equation mode according the way input
!mode 2000<<3000: get coefficient equation
!mode 3000<<4000: get labda equation
!mode 7000<<8000: calculate one term without multiplying with coefficient
!mode =-1001 : get number of equations with unknwon values
!mode =-1002 : get number of total terms in the equatons !
!
!   This file is generated from f77 code, using a code generator which
!   transposes from the F77 code into F90
!   the code USES module file as used in the BFM model
!   F90 code generator written by P. Ruardij.
!
! AUTHORS
!    Piet Ruardij
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
      REAL(RLEN) FUNCTION noutput(NUTR,mode,option,input,xinput,yinput)
        USE global_mem, ONLY:RLEN
        USE mem, ONLY:NO_BOXES_XY
        USE bennut_type
        USE bennut_variables
        USE bennut_constants
        USE constants
        USE bennut_interface,ONLY:calculate_equation,transfer
        IMPLICIT  NONE
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::option ! Specification
        integer,intent(IN) ::input ! Specification
        integer,intent(IN) ::nutr ! Specification
        REAL(RLEN),intent(IN) ::xinput ! Specification
        REAL(RLEN),intent(IN) ::yinput ! Specification

        REAL(RLEN) ::r
        REAL(RLEN) ::bC
        REAL(RLEN) ::fC(10)
        integer ::l,m

        !calculation of the result, in case of integration calculation of the &
        ! result at
        !the under bborder
        if (option == INTEGRAL) then
            if (abs(xinput-yinput)/(.5D+00*(xinput+yinput)) < 1.D-6) then
              noutput=0.D+00
              return
            endif
        endif
        !Get equation number or term number :
        !CALCULATE EQUATION:
        !Get first term number:
        l=sets(NUTR)%lst(mode)
        !Get last equation number:
        m=sets(NUTR)%lfi(mode)-l+1
        bC=sets(NUTR)%b(mode)
        fC(1:m)=sets(NUTR)%factor(l:l+m)
        !r=calculate_equation(option,xinput,sets(NUTR)%coeffs(l), &
        !                                                     bC,fC,m)
	r = 0
        if (option == INTEGRAL.or.option == EXPONENTIAL_INTEGRAL) then
          !r= calculate_equation(option,yinput,sets(NUTR)%coeffs(l),&
          !                                                bC,fC,m) -r
	  r = 0
          if (input == RFLUX.or. input == MASS) then
            r=r*sets(NUTR)%poro(mode)
            if(input == MASS)r=r*(sets(NUTR)%ads(mode)+1.D+00)
          endif
        elseif (option == DERIVATIVE) then
          if (input == RFLUX) &
             r=r*sets(NUTR)%diff(mode)*sets(NUTR)%poro(mode)
        endif
        noutput=r
      end
