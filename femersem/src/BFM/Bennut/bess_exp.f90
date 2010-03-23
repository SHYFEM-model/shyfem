!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   bess_exp.f90
!
! FILE
!   bess_exp.f90
!
! DESCRIPTION
!   FILE: Solve_coupled.F  nutrient dynamics in the benthos 
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
      REAL(RLEN) FUNCTION BESS_EXP(X,LABDA,FN)
        USE global_mem, ONLY:RLEN
        USE constants
        IMPLICIT  NONE
        REAL(RLEN),intent(IN)    ::x ! Specification
        REAL(RLEN),intent(IN)    ::labda(2) ! Specification

        INTERFACE                             ! Specification
          REAL(RLEN) FUNCTION FN(X )          ! Specification
          USE global_mem, ONLY:RLEN
          REAL(RLEN),INTENT(IN)   ::X              ! Specification
        END FUNCTION                         ! Specification
        END  INTERFACE                       ! Specification     

        real(RLEN)  :: r
        real(RLEN)  :: s

        r=labda(1)*x
        r=dexp(r)
        s=r*labda(2)

        BESS_EXP=FN(X)
        RETURN
      end

