!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   transfer.f90
!
! FILE
!   transfer.f90
!
! DESCRIPTION
!   
!  
!   This file is generated from f77 code, using a code generator which
!   transposes from the F77 code into F90
!   the code USES module file as used in the BFM model
!   F90 code generator written by P. Ruardij. 
!
! AUTHORS
!   
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
      REAL(RLEN) FUNCTION transfer(mode,coeff,input,diff)
        USE bennut_constants
        USE bennut_type
        USE constants
        IMPLICIT  NONE
        integer,intent(IN) ::mode ! Specification
        type (ty_coeff),intent(IN) ::coeff ! Specification
        REAL(RLEN),intent(IN) ::diff ! Specification
        REAL(RLEN),intent(IN) ::input ! Specification

        REAL(RLEN) ::r
        integer ::term 

        term=coeff%ia
        if (input.eq.0.0D+00) then
          transfer=0.0D+00
          return
        endif

        select case (term)
          case (ZERO_EXPONENTIAL_TERM)
            r=-coeff%labda(1)**2*diff
            if (mode.eq.PARA2COEFF) then
              transfer=input/r
            else
              transfer=input*r
            endif
          case (QUADRATIC_TERM)
            r=-2.D+00*diff
            if (mode.eq.PARA2COEFF) then
              transfer=input/r
            else
              transfer=input*r
            endif
          case ( CONSTANT_TERM,LINEAR_TERM) 
            transfer=input
          case (BESSELI_EXP_TERM, BESSELK_EXP_TERM)
            if (mode.eq.PARA2COEFF) then
              transfer=sqrt(input/diff)/abs(coeff%labda(1))/2.D+00
            elseif (mode.eq.LABDA_2) then
              r=coeff%labda(1) 
              transfer = 2.0 / abs( r ) * sqrt(coeff%labda(2) / diff);
            elseif (mode.eq.LABDA_1) then
              transfer = coeff%labda(1) / 2.0 
            else
              stop 'transfer'
            endif
          case default
            stop 'transfer '
        end select

        return
      end

