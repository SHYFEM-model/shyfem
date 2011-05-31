!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   calculate_equation.f90
!
! FILE
!   calculate_equation.f90
!
! DESCRIPTION
!   mode <1000 calculate equation mode
!   mode >1000 calulate one term (term is input as a number existing
!   of equationumber*10+termnumber
!  
!   This file is generated from f77 code, using a code generator which
!   transposes from the F77 code into F90
!   the code USES module file as used in the BFM model
!   F90 code generator written by P. Ruardij. 
!
! AUTHORS
!   Piet Ruardij ( rua@nioz.nl)
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
      REAL(RLEN) FUNCTION calculate_equation(mode,x,coeffs,b,factor,nn)
        USE global_mem, ONLY:RLEN
        USE bennut_type, ONLY:ty_coeff
        USE bennut_interface,ONLY: funcalc
        IMPLICIT  NONE
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::nn ! Specification
        type (ty_coeff),intent(IN) ::coeffs(nn) ! Specification
        REAL(RLEN),intent(IN) ::x ! Specification
        REAL(RLEN),intent(IN) ::b ! Specification
        REAL(RLEN),intent(IN) ::factor(nn) ! Specification

        integer ::j
        REAL(RLEN) ::r

          r=0.0D+00
          do j=1,nn
            !find all terms equation mode for x
            r=r +funcalc(mode,0,coeffs(j),b,x)*factor(j)
          enddo

        calculate_equation=r
        return
      end

