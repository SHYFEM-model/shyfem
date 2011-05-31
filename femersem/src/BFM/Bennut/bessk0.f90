!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   bessk0.f90
!
! FILE
!   bessk0.f90
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
!   Numerical Recipes Software 
!
! CHANGE_LOG
!
! COPYING
!   
!   (C) Copr. 1986-92 Numerical Recipes Software .
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
      REAL(RLEN) FUNCTION BESSK0(X)
        USE global_mem, ONLY:RLEN
        USE bennut_interface,ONLY:BESSI0
        IMPLICIT  NONE
        REAL(RLEN),intent(IN) ::x ! Specification

        REAL(RLEN) ::y
        REAL(RLEN) ::p1
        REAL(RLEN) ::p2
        REAL(RLEN) ::p3
        REAL(RLEN) ::p4
        REAL(RLEN) ::p5
        REAL(RLEN) ::p6
        REAL(RLEN) ::p7
        REAL(RLEN) ::q1
        REAL(RLEN) ::q2
        REAL(RLEN) ::q3
        REAL(RLEN) ::q4
        REAL(RLEN) ::q5
        REAL(RLEN) ::q6
        REAL(RLEN) ::q7
        REAL(RLEN) ::bx
        REAL(RLEN) ::ax
        parameter( P1=-0.57721566D0,P2=0.42278420D0,P3=0.23069756D0, &
        P4=0.3488590D-1,P5=0.262698D-2,P6=0.10750D-3,P7=0.74D-5)
        parameter( Q1=1.25331414D0,Q2=-0.7832358D-1,Q3=0.2189568D-1, &
        Q4=-0.1062446D-1,Q5=0.587872D-2,Q6=-0.251540D-2,Q7=0.53208D-3)

        IF (X.LE.2.0D0) THEN
          BX=BESSI0(X)
          Y=X*X/4.0D0
          AX=X/2.0D0
          BESSK0=(-LOG(AX)*BX)+(P1+Y*(P2+Y*(P3+ &
          Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
        ELSE
          Y=(2.0D0/X)
          BESSK0=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+ &
          Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
        ENDIF
        RETURN
      END

