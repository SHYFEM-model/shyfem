!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   bessi0.f90
!
! FILE
!   bessi0.f90
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
      REAL(RLEN) FUNCTION BESSI0(X)
        USE constants
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
        REAL(RLEN) ::ax
        REAL(RLEN) ::q1
        REAL(RLEN) ::q2
        REAL(RLEN) ::q3
        REAL(RLEN) ::q4
        REAL(RLEN) ::q5
        REAL(RLEN) ::q6
        REAL(RLEN) ::q7
        REAL(RLEN) ::q8
        REAL(RLEN) ::q9
        parameter (P1=1.0D0,P2=3.5156229D0,P3=3.0899424D0, &
        P4=1.2067492D0, P5=0.2659732D0,P6=0.360768D-1,P7=0.45813D-2)
        parameter( Q1=0.39894228D0,Q2=0.1328592D-1, &
        Q3=0.2253190D-2,Q4=-0.1575650D-2,Q5=0.916281D-2,Q6=-0.2057706D-1, &
        Q7=0.2635537D-1,Q8=-0.1647633D-1,Q9=0.392377D-2)

        IF (ABS(X).LT.3.75D+00) THEN
          Y=(X/3.75D0)**2
          BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
        ELSE
          AX=ABS(X)
          Y=3.75D0/AX
          BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4 &
          +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
        ENDIF
        RETURN
      END

