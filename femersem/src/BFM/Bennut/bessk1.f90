!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   bessk1.f90
!
! FILE
!   bessk1.f90
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
      REAL(RLEN) FUNCTION BESSK1(X)
       USE constants
       USE bennut_interface,ONLY:BESSI1
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
        parameter( P1=1.0D0,P2=0.15443144D0,P3=-0.67278579D0, &
        P4=-0.18156897D0,P5=-0.1919402D-1,P6=-0.110404D-2,P7=-0.4686D-4)
        parameter( Q1=1.25331414D0,Q2=0.23498619D0,Q3=-0.3655620D-1, &
        Q4=0.1504268D-1,Q5=-0.780353D-2,Q6=0.325614D-2,Q7=-0.68245D-3)

        IF (X.LE.2.0D+00) THEN
          BX=BESSI1(X)
          Y=X*X/4.0D+00
          AX=X/2.0D+00
          BESSK1=(LOG(AX)*BX)+(1.0D+00/X)*(P1+Y*(P2+ &
          Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
        ELSE
          Y=2.0D+00/X
          BESSK1=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+ &
          Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
        ENDIF
        RETURN
      END
