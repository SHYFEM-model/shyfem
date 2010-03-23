!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! SUBROUTINE
!   svdcmp.f90
!
! FILE
!   svdcmp.f90
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
!   {!  (C) Copr. 1986-92 Numerical Recipes Software .)1:.}
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
      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
        USE global_mem, ONLY:RLEN
        parameter (NMAX=100)
        integer,intent(INOUT) ::m ! Specification
        integer,intent(INOUT) ::n ! Specification
        integer,intent(INOUT) ::mp ! Specification
        integer,intent(INOUT) ::np ! Specification
        REAL(RLEN),intent(INOUT) ::a(MP,NP) ! Specification
        REAL(RLEN),intent(INOUT) ::w(NP) ! Specification
        REAL(RLEN),intent(INOUT) ::v(NP,NP) ! Specification
        REAL(RLEN) ::rv1(NMAX)

        integer ::i
        integer ::l
        integer ::j
        integer ::k
        integer ::nm
        integer ::its
        REAL(RLEN) ::g
        REAL(RLEN) ::scale
        REAL(RLEN) ::anorm
        REAL(RLEN) ::f
        REAL(RLEN) ::s
        REAL(RLEN) ::h
        REAL(RLEN) ::c
        REAL(RLEN) ::y
        REAL(RLEN) ::z
        G=0.0D+00
        SCALE=0.0D+00
        ANORM=0.0D+00
        DO 25 I=1,N
          L=I+1
          RV1(I)=SCALE*G
          G=0.0D+00
          S=0.0D+00
          SCALE=0.0D+00
          IF (I.LE.M) THEN
            DO 11 K=I,M
              SCALE=SCALE+ABS(A(K,I))
 11         CONTINUE
            IF (SCALE.NE.0.0D+00) THEN
              DO 12 K=I,M
                A(K,I)=A(K,I)/SCALE
                S=S+A(K,I)*A(K,I)
 12           CONTINUE
              F=A(I,I)
              G=-SIGN(SQRT(S),F)
              H=F*G-S
              A(I,I)=F-G
              IF (I.NE.N) THEN
                DO 15 J=L,N
                  S=0.0D+00
                  DO 13 K=I,M
                    S=S+A(K,I)*A(K,J)
 13               CONTINUE
                  F=S/H
                  DO 14 K=I,M
                    A(K,J)=A(K,J)+F*A(K,I)
 14               CONTINUE
 15             CONTINUE
              ENDIF
              DO 16 K= I,M
                A(K,I)=SCALE*A(K,I)
 16           CONTINUE
            ENDIF
          ENDIF
          W(I)=SCALE *G
          G=0.0D+00
          S=0.0D+00
          SCALE=0.0D+00
          IF ((I.LE.M).AND.(I.NE.N)) THEN
            DO 17 K=L,N
              SCALE=SCALE+ABS(A(I,K))
 17         CONTINUE
            IF (SCALE.NE.0.0D+00) THEN
              DO 18 K=L,N
                A(I,K)=A(I,K)/SCALE
                S=S+A(I,K)*A(I,K)
 18           CONTINUE
              F=A(I,L)
              G=-SIGN(SQRT(S),F)
              H=F*G-S
              A(I,L)=F-G
              DO 19 K=L,N
                RV1(K)=A(I,K)/H
 19           CONTINUE
              IF (I.NE.M) THEN
                DO 23 J=L,M
                  S=0.0D+00
                  DO 21 K=L,N
                    S=S+A(J,K)*A(I,K)
 21               CONTINUE
                  DO 22 K=L,N
                    A(J,K)=A(J,K)+S*RV1(K)
 22               CONTINUE
 23             CONTINUE
              ENDIF
              DO 24 K=L,N
                A(I,K)=SCALE*A(I,K)
 24           CONTINUE
            ENDIF
          ENDIF
          ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
 25     CONTINUE
        DO 32 I=N,1,-1
          IF (I.LT.N) THEN
            IF (G.NE.0.0D+00) THEN
              DO 26 J=L,N
                V(J,I)=(A(I,J)/A(I,L))/G
 26           CONTINUE
              DO 29 J=L,N
                S=0.0D+00
                DO 27 K=L,N
                  S=S+A(I,K)*V(K,J)
 27             CONTINUE
                DO 28 K=L,N
                  V(K,J)=V(K,J)+S*V(K,I)
 28             CONTINUE
 29           CONTINUE
            ENDIF
            DO 31 J=L,N
              V(I,J)=0.0D+00
              V(J,I)=0.0D+00
 31         CONTINUE
          ENDIF
          V(I,I)=1.0D+00
          G=RV1(I)
          L=I
 32     CONTINUE
        DO 39 I=N,1,-1
          L=I+1
          G=W(I)
          IF (I.LT.N) THEN
            DO 33 J=L,N
              A(I,J)=0.0D+00
 33         CONTINUE
          ENDIF
          IF (G.NE.0.0D+00) THEN
            G=1.0D+00/G
            IF (I.NE.N) THEN
              DO 36 J=L,N
                S=0.0D+00
                DO 34 K=L,M
                  S=S+A(K,I)*A(K,J)
 34             CONTINUE
                F=(S/A(I,I))*G
                DO 35 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
 35             CONTINUE
 36           CONTINUE
            ENDIF
            DO 37 J=I,M
              A(J,I)=A(J,I)*G
 37         CONTINUE
          ELSE
            DO 38 J= I,M
              A(J,I)=0.0D+00
 38         CONTINUE
          ENDIF
          A(I,I)=A(I,I)+1.0D+00
 39     CONTINUE
        DO 49 K=N,1,-1
          DO 48 ITS=1,30
            DO 41 L=K,1,-1
              NM=L-1
              IF ((ABS(RV1(L))+ANORM).EQ.ANORM) GO TO 2
              IF ((ABS(W(NM))+ANORM).EQ.ANORM) GO TO 1
 41         CONTINUE
 1          C=0.0D+00
            S=1.0D+00
            DO 43 I=L,K
              F=S*RV1(I)
              IF ((ABS(F)+ANORM).NE.ANORM) THEN
                G=W(I)
                H=SQRT(F*F+G*G)
                W(I)=H
                H=1.0D+00/H
                C= (G*H)
                S=-(F*H)
                DO 42 J=1,M
                  Y=A(J,NM)
                  Z=A(J,I)
                  A(J,NM)=(Y*C)+(Z*S)
                  A(J,I)=-(Y*S)+(Z*C)
 42             CONTINUE
              ENDIF
 43         CONTINUE
 2          Z=W(K)
            IF (L.EQ.K) THEN
              IF (Z.LT.0.0D+00) THEN
                W(K)=-Z
                DO 44 J=1,N
                  V(J,K)=-V(J,K)
 44             CONTINUE
              ENDIF
              GO TO 3
            ENDIF
            IF (ITS.EQ.30) THEN
              WRITE(0,*) 'No convergence in 30 iterations'
              stop 'svdcmp'
            ENDIF
            X=W(L)
            NM=K-1
            Y=W(NM)
            G=RV1(NM)
            H=RV1(K)
            F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0D+00*H*Y)
            G=SQRT(F*F+1.0D+00)
            F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
            C=1.0D+00
            S=1.0D+00
            DO 47 J=L,NM
              I=J+1
              G=RV1(I)
              Y=W(I)
              H=S*G
              G=C*G
              Z=SQRT(F*F+H*H)
              RV1(J)=Z
              C=F/Z
              S=H/Z
              F= (X*C)+(G*S)
              G=-(X*S)+(G*C)
              H=Y*S
              Y=Y*C
              DO 45 NM=1,N
                X=V(NM,J)
                Z=V(NM,I)
                V(NM,J)= (X*C)+(Z*S)
                V(NM,I)=-(X*S)+(Z*C)
 45           CONTINUE
              Z=SQRT(F*F+H*H)
              W(J)=Z
              IF (Z.NE.0.0D+00) THEN
                Z=1.0D+00/Z
                C=F*Z
                S=H*Z
              ENDIF
              F= (C*G)+(S*Y)
              X=-(S*G)+(C*Y)
              DO 46 NM=1,M
                Y=A(NM,J)
                Z=A(NM,I)
                A(NM,J)= (Y*C)+(Z*S)
                A(NM,I)=-(Y*S)+(Z*C)
 46           CONTINUE
 47         CONTINUE
            RV1(L)=0.0D+00
            RV1(K)=F
            W(K)=X
 48       CONTINUE
 3        CONTINUE
 49     CONTINUE
        RETURN
      END

