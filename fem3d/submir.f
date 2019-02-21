
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c linear equation solvers
c
c contents :
c
c      SUBROUTINE LEQT1B (A,N,NLC,NUC,IA,B,M,IB,IJOB,XL,IER)
c      SUBROUTINE LEQT2B (A,N,NLC,NUC,IA,B,M,IB,IJOB,U,IU,XL,IER)
c      SUBROUTINE LEQ1PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,IER)
c      SUBROUTINE LEQ2PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,WK,IER)
c      SUBROUTINE LUDAPB (A,N,NC,IA,UL,IU,D1,D2,IER)
c      SUBROUTINE LUELPB (UL,B,N,NC,IA,X)
c      SUBROUTINE LUREPB (A,N,NC,IA,UL,IU,B,X,IDGT,RES,IER)
c      SUBROUTINE DEQT1B (A,N,NLC,NUC,IA,B,M,IB,IJOB,XL,IER)
c      SUBROUTINE DEQT2B (A,N,NLC,NUC,IA,B,M,IB,IJOB,U,IU,XL,IER)
c      SUBROUTINE DEQ1PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,IER)
c      SUBROUTINE DEQ2PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,WK,IER)
c      SUBROUTINE DUDAPB (A,N,NC,IA,UL,IU,D1,D2,IER)
c      SUBROUTINE DUELPB (UL,B,N,NC,IA,X)
c      SUBROUTINE DUREPB (A,N,NC,IA,UL,IU,B,X,IDGT,RES,IER)
c      SUBROUTINE UERTST (IER,NAME)
c      SUBROUTINE UGETIO (IOPT,NIN,NOUT)
c      SUBROUTINE GELB(R,A,M,N,MUD,MLD,EPS,IER)
c      SUBROUTINE DGELB(R,A,M,N,MUD,MLD,EPS,IER)
c      SUBROUTINE MCHB(R,A,M,N,MUD,IOP,EPS,IER)
c      SUBROUTINE DMCHB(R,A,M,N,MUD,IOP,EPS,IER)
c      subroutine loctst(i,j,n,m)
c      function locmy(i,j,n,m)
c      function locimm(i,j,n,m)
c      function loccer(i,j,n,m)
c      function locssp(i,j,n,m)
c      function locsps(i,j,n,m)
c
c revision log :
c
c 03.04.1997	ggu	general - compiler warnings for gfortran
c 24.05.1997	ggu	general - compiler warnings -> call to dp routines
c
c********************************************************************

      SUBROUTINE LEQT1B (A,N,NLC,NUC,IA,B,M,IB,IJOB,XL,IER)
C
      DIMENSION          A(IA,*),XL(N,*),B(IB,*)
      DATA               ZERO/0./,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JBEG = NLC+1
      NLC1 = JBEG
      IF (IJOB .EQ. 2) GO TO 80
      RN = N
C                                  RESTRUCTURE THE MATRIX
C                                  FIND RECIPROCAL OF THE LARGEST
C                                  ABSOLUTE VALUE IN ROW I
      I = 1
      NC = JBEG+NUC
      NN = NC
      JEND = NC
      IF (N .EQ. 1 .OR. NLC .EQ. 0) GO TO 25
    5 K = 1
      P = ZERO
      DO 10 J = JBEG,JEND
         A(I,K) = A(I,J)
         Q =  ABS(A(I,K))
         IF (Q .GT. P) P = Q
         K = K+1
   10 CONTINUE
      IF (P .EQ. ZERO) GO TO 135
      XL(I,NLC1) = ONE/P
      IF (K .GT. NC) GO TO 20
      DO 15 J = K,NC
         A(I,J) = ZERO
   15 CONTINUE
   20 I = I+1
      JBEG = JBEG-1
      IF (JEND-JBEG .EQ. N) JEND = JEND-1
      IF (I .LE. NLC) GO TO 5
      JBEG = I
      NN = JEND
   25 JEND = N-NUC
      DO 40 I = JBEG,N
         P = ZERO
         DO 30 J = 1,NN
            Q =  ABS(A(I,J))
            IF (Q .GT. P) P = Q
   30    CONTINUE
         IF (P .EQ. ZERO) GO TO 135
         XL(I,NLC1) = ONE/P
         IF (I .EQ. JEND) GO TO 37
         IF (I .LT. JEND) GO TO 40
         K = NN+1
         DO 35 J = K,NC
            A(I,J) = ZERO
   35    CONTINUE
   37    NN = NN-1
   40 CONTINUE
      L = NLC
C                                  L-U DECOMPOSITION
      DO 75 K = 1,N
         P =  ABS(A(K,1))*XL(K,NLC1)
         I = K
         IF (L .LT. N) L = L+1
         K1 = K+1
         IF (K1 .GT. L) GO TO 50
         DO 45 J = K1,L
            Q = ABS(A(J,1))*XL(J,NLC1)
            IF (Q .LE. P) GO TO 45
            P = Q
            I = J
   45    CONTINUE
   50    XL(I,NLC1) = XL(K,NLC1)
         XL(K,NLC1) = I
C                                  SINGULARITY FOUND
         Q = RN+P
         IF (Q .EQ. RN) GO TO 135
C                                  INTERCHANGE ROWS I AND K
         IF (K .EQ. I) GO TO 60
         DO 55 J = 1,NC
            P = A(K,J)
            A(K,J) = A(I,J)
            A(I,J) = P
   55    CONTINUE
   60    IF (K1 .GT. L) GO TO 75
         DO 70 I = K1,L
            P = A(I,1)/A(K,1)
            IK = I-K
            XL(K1,IK) = P
            DO 65 J = 2,NC
               A(I,J-1) = A(I,J)-P*A(K,J)
   65    CONTINUE
         A(I,NC) = ZERO
   70    CONTINUE
   75 CONTINUE
      IF (IJOB .EQ. 1) GO TO 9005
C                                  FORWARD SUBSTITUTION
   80 L = NLC
      DO 105 K = 1,N
         I = XL(K,NLC1)
         IF (I .EQ. K) GO TO 90
         DO 85 J = 1,M
            P = B(K,J)
            B(K,J) = B(I,J)
            B(I,J) = P
   85    CONTINUE
   90    IF (L .LT. N) L = L+1
         K1 = K+1
         IF (K1 .GT. L) GO TO 105
         DO 100 I = K1,L
            IK = I-K
            P = XL(K1,IK)
            DO 95 J = 1,M
               B(I,J) = B(I,J)-P*B(K,J)
   95       CONTINUE
  100    CONTINUE
  105 CONTINUE
C                                  BACKWARD SUBSTITUTION
      JBEG = NUC+NLC
      DO 125 J = 1,M
         L = 1
         K1 = N+1
         DO 120 I = 1,N
            K = K1-I
            P = B(K,J)
            IF (L .EQ. 1) GO TO 115
            DO 110 KK = 2,L
               IK = KK+K
               P = P-A(K,KK)*B(IK-1,J)
  110       CONTINUE
  115       B(K,J) = P/A(K,1)
            IF (L .LE. JBEG) L = L+1
  120    CONTINUE
  125 CONTINUE
      GO TO 9005
  135 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,'LEQT1B')
 9005 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE LEQT2B (A,N,NLC,NUC,IA,B,M,IB,IJOB,U,IU,XL,IER)
C
      DIMENSION          A(IA,*),U(IU,*),XL(N,*),B(IB,*)
      DOUBLE PRECISION   SUM
      DATA               ZERO/0.0/,ITMAX/50/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      NLCP1 = NLC+1
      NCC = NUC+NLCP1
      NLCP2 = NLC+2
      NMNUC = N-NUC
      NU2 = NCC+1
      NU3 = NCC+2
      IF (IJOB .EQ. 2) GO TO 15
C                                  SAVE MATRIX A
      DO 10 I = 1,N
         DO 5 J = 1,NCC
            U(I,J) = A(I,J)
    5    CONTINUE
   10 CONTINUE
C                                  FACTOR MATRIX A
      CALL LEQT1B(U,N,NLC,NUC,IU,B,M,IB,1,XL,IER)
      IF (IER .NE. 0) GO TO 9000
      IF (IJOB .EQ. 1) GO TO 9005
C                                  SAVE THE RIGHT HAND SIDES
   15 DO 60 J = 1,M
         DO 20 I = 1,N
            U(I,NU2) = B(I,J)
   20    CONTINUE
C                                  OBTAIN A SOLUTION
         CALL LEQT1B(U,N,NLC,NUC,IU,U(1,NU2),1,IU,2,XL,IER)
C                                  COMPUTE THE NORM OF THE SOLUTION
         XNORM = ZERO
         DO 25 I = 1,N
            XNORM = MAX(XNORM, ABS(U(I,NU2)))
   25    CONTINUE
         IF (XNORM .EQ. ZERO) GO TO 60
C                                  COMPUTE THE RESIDUALS
         DO 45 ITER = 1,ITMAX
            NC = NCC
            KK = 1
            DO 35 I = 1,N
               SUM = B(I,J)
               L = NLCP2-I
               IR = MAX0(L,1)
               IF (L .LE. 0) KK = KK+1
               K = KK
               DO 30 JJ = IR,NC
                  SUM = SUM-DBLE(A(I,JJ))*DBLE(U(K,NU2))
                  K = K+1
   30          CONTINUE
               U(I,NU3) = SUM
               IF (I .GE. NMNUC) NC = NC-1
   35       CONTINUE
            CALL LEQT1B(U,N,NLC,NUC,IU,U(1,NU3),1,IU,2,XL,IER)
            DXNORM = ZERO
C                                  UPDATE THE SOLUTION
            DO 40 I = 1,N
               U(I,NU2) = U(I,NU2)+U(I,NU3)
               DXNORM = MAX(DXNORM, ABS(U(I,NU3)))
   40       CONTINUE
            IF (XNORM+DXNORM .EQ. XNORM) GO TO 50
   45    CONTINUE
         IER = 130
C                                  STORE THE SOLUTION
   50    DO 55 JK = 1,N
            B(JK,J) = U(JK,NU2)
   55    CONTINUE
         IF (IER .NE. 0) GO TO 9000
   60 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'LEQT2B')
 9005 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE LEQ1PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,IER)
C
      REAL               A(IA,*),B(IB,*),D1,D2
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0*IDGT
C                                  DECOMPOSITION OF MATRIX A INTO
C                                  L*L-TRANSPOSE
      CALL LUDAPB(A,N,NC,IA,A,IA,D1,D2,IER)
      IF (IER .NE. 0) GO TO 9000
      DO 5 I = 1,M
C                                  SOLUTION OF AX = B
         CALL LUELPB(A,B(1,I),N,NC,IA,B(1,I))
    5 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'LEQ1PB')
 9005 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE LEQ2PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,WK,IER)
C
      REAL               A(IA,*),B(IB,*),WK(N,*),D1,D2
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      K = NC+2
      K1 = K+1
      CALL LUDAPB(A,N,NC,IA,WK,N,D1,D2,IER)
      IF (IER .NE. 0) GO TO 9000
         DO 10 I = 1,M
         CALL LUELPB(WK,B(1,I),N,NC,N,WK(1,K))
         CALL LUREPB(A,N,NC,IA,WK,N,B(1,I),WK(1,K),IDGT,WK(1,K1),IER)
            DO 5 J = 1,N
            B(J,I) = WK(J,K)
    5       CONTINUE
         IF (IER .NE. 0) GO TO 15
   10    CONTINUE
      GO TO 9005
   15 IER = 130
 9000 CONTINUE
      CALL UERTST(IER,'LEQ2PB')
 9005 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE LUDAPB (A,N,NC,IA,UL,IU,D1,D2,IER)
C
      REAL               ZERO,A(IA,*),UL(IU,*),D1,D2,ONE,SUM,
     +                   RN,FOUR,SIXTN,SIXTH
      DATA               ZERO/0.0/,FOUR/4.0/,SIXTN/16./,
     +                   SIXTH/.0625/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      RN = ONE/(N*SIXTN)
      D1 = ONE
      D2 = ZERO
      NCP1 = NC+1
      IF (NC .EQ. 0) GO TO 15
C                                  INITIALIZE ZERO ELEMENTS
      DO 10 I = 1,NC
         DO 5 J = I,NC
            K = NCP1-J
            UL(I,K) = ZERO
    5    CONTINUE
   10 CONTINUE
C                                  I IS ROW INDEX OF ELEMENT BEING
C                                  COMPUTED
   15 DO 60 I = 1,N
         IMNCP1 = I-NCP1
         I1 = MAX0(1,1-IMNCP1)
C                                  J IS COLUMN INDEX OF ELEMENT BEING
C                                  COMPUTED
         DO 60 J = I1,NCP1
C                                  L IS ROW INDEX OF PREVIOUSLY COMPUTED
C                                  VECTOR BEING USED TO COMPUTE INNER
C                                  PRODUCT
            L = IMNCP1+J
            I2 = NCP1-J
            SUM = A(I,J)
            JM1 = J-1
            IF (JM1) 30,30,20
   20       DO 25 K = 1,JM1
C                                  M IS COLUMN INDEX
               M = I2+K
               SUM = SUM-UL(I,K)*UL(L,M)
   25       CONTINUE
   30       IF (J .NE. NCP1) GO TO 55
            IF(A(I,J)+SUM*RN .LE. A(I,J))GO TO 65
            UL(I,J) = ONE/SQRT(SUM)
C                                  UPDATE THE DETERMINANT
            D1 = D1*SUM
   35       IF (ABS(D1)-ONE) 45,45,40
   40       D1 = D1*SIXTH
            D2 = D2+FOUR
            GO TO 35
   45       IF (ABS(D1)-SIXTH) 50,50,60
   50       D1 = D1*SIXTN
            D2 = D2-FOUR
            GO TO 45
   55       UL(I,J) = SUM*UL(L,NCP1)
   60 CONTINUE
      GO TO 9005
   65 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,'LUDAPB')
 9005 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE LUELPB (UL,B,N,NC,IA,X)
C
      REAL               UL(IA,*),B(*),X(*),ZERO,SUM
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  SOLUTION LY = B
      NC1 = NC+1
      IW = 0
      L = 0
      DO 15 I = 1,N
         SUM = B(I)
         IF (NC .LE. 0) GO TO 10
         IF (IW .EQ. 0) GO TO 9
         L = L+1
         IF (L .GT. NC) L = NC
         K = NC1-L
         KL = I-L
         DO 5 J = K,NC
            SUM = SUM -X(KL) * UL(I,J)
            KL = KL+1
    5    CONTINUE
         GO TO 10
    9    IF (SUM .NE. ZERO) IW = 1
   10    X(I) = SUM*UL(I,NC1)
   15 CONTINUE
C                                  SOLUTION UX = Y
   20 X(N) = X(N)*UL(N,NC1)
      IF (N .LE. 1) GO TO 40
      N1 = N+1
      DO 35 I = 2,N
         K = N1-I
         SUM = X(K)
         IF (NC .LE. 0) GO TO 30
         KL = K+1
         K1 = MIN0(N,K+NC)
         L = 1
         DO 25 J = KL,K1
            SUM = SUM -X(J) * UL(J,NC1-L)
            L = L+1
   25    CONTINUE
   30    X(K) = SUM*UL(K,NC1)
   35 CONTINUE
   40 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE LUREPB (A,N,NC,IA,UL,IU,B,X,IDGT,RES,IER)
C
      REAL               A(IA,*),B(*),UL(IU,*),RES(*),X(*),XNORM,
     +                   ZERO,DXNORM
      DOUBLE PRECISION   SUM
      DATA               ITMAX/50/
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      NC1 = NC+1
      XNORM = ZERO
      DO 5 I = 1,N
         XNORM = MAX(XNORM,ABS(X(I)))
    5 CONTINUE
      IF (XNORM .NE. ZERO) GO TO 10
      IDGT = 50
      GO TO 9005
C                                  ITERATION LOOP
   10 DO 45 ITER = 1,ITMAX
         L = 0
         DO 30 I = 1,N
            K = MAX0(1,NC1-L)
            LK = MAX0(1,I-NC)
            L = I
            SUM = B(I)
            DO 15 J = K,NC1
               SUM = SUM - DBLE(A(I,J))*DBLE(X(LK))
               LK = LK+1
   15       CONTINUE
            LL = MIN0(NC,N-I)
            IF (LL .LT. 1) GO TO 25
            DO 20 J = 1,LL
               SUM = SUM-DBLE(A(I+J,NC1-J))*DBLE(X(LK))
               LK = LK+1
   20       CONTINUE
   25       RES(I) = SUM
   30    CONTINUE
         CALL LUELPB(UL,RES,N,NC,IU,RES)
         DXNORM = ZERO
         XNORM = ZERO
         DO 35 I = 1,N
            X(I) = X(I)+RES(I)
            DXNORM = MAX(DXNORM,ABS(RES(I)))
            XNORM = MAX(XNORM,ABS(X(I)))
   35    CONTINUE
         IF (ITER .NE. 1) GO TO 40
         IF (DXNORM .GT. ZERO) IDGT = -LOG10(DXNORM/XNORM)
   40    IF (XNORM + DXNORM .EQ. XNORM)GO TO 9005
   45 CONTINUE
C                                  ITERATIVE IMPROVEMENT FAILED
      IER = 129
 9000 CONTINUE
      CALL UERTST(IER,'LUREPB')
 9005 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE DEQT1B (A,N,NLC,NUC,IA,B,M,IB,IJOB,XL,IER)
C
	implicit integer (i-n)
	implicit double precision (a-h)
	implicit double precision (o-z)
      DIMENSION          A(IA,*),XL(N,*),B(IB,*)
      DATA               ZERO/0./,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JBEG = NLC+1
      NLC1 = JBEG
      IF (IJOB .EQ. 2) GO TO 80
      RN = N
C                                  RESTRUCTURE THE MATRIX
C                                  FIND RECIPROCAL OF THE LARGEST
C                                  DABSOLUTE VALUE IN ROW I
      I = 1
      NC = JBEG+NUC
      NN = NC
      JEND = NC
      IF (N .EQ. 1 .OR. NLC .EQ. 0) GO TO 25
    5 K = 1
      P = ZERO
      DO 10 J = JBEG,JEND
         A(I,K) = A(I,J)
         Q =  DABS(A(I,K))
         IF (Q .GT. P) P = Q
         K = K+1
   10 CONTINUE
      IF (P .EQ. ZERO) GO TO 135
      XL(I,NLC1) = ONE/P
      IF (K .GT. NC) GO TO 20
      DO 15 J = K,NC
         A(I,J) = ZERO
   15 CONTINUE
   20 I = I+1
      JBEG = JBEG-1
      IF (JEND-JBEG .EQ. N) JEND = JEND-1
      IF (I .LE. NLC) GO TO 5
      JBEG = I
      NN = JEND
   25 JEND = N-NUC
      DO 40 I = JBEG,N
         P = ZERO
         DO 30 J = 1,NN
            Q =  DABS(A(I,J))
            IF (Q .GT. P) P = Q
   30    CONTINUE
         IF (P .EQ. ZERO) GO TO 135
         XL(I,NLC1) = ONE/P
         IF (I .EQ. JEND) GO TO 37
         IF (I .LT. JEND) GO TO 40
         K = NN+1
         DO 35 J = K,NC
            A(I,J) = ZERO
   35    CONTINUE
   37    NN = NN-1
   40 CONTINUE
      L = NLC
C                                  L-U DECOMPOSITION
      DO 75 K = 1,N
         P =  DABS(A(K,1))*XL(K,NLC1)
         I = K
         IF (L .LT. N) L = L+1
         K1 = K+1
         IF (K1 .GT. L) GO TO 50
         DO 45 J = K1,L
            Q = DABS(A(J,1))*XL(J,NLC1)
            IF (Q .LE. P) GO TO 45
            P = Q
            I = J
   45    CONTINUE
   50    XL(I,NLC1) = XL(K,NLC1)
         XL(K,NLC1) = I
C                                  SINGULARITY FOUND
         Q = RN+P
         IF (Q .EQ. RN) GO TO 135
C                                  INTERCHANGE ROWS I AND K
         IF (K .EQ. I) GO TO 60
         DO 55 J = 1,NC
            P = A(K,J)
            A(K,J) = A(I,J)
            A(I,J) = P
   55    CONTINUE
   60    IF (K1 .GT. L) GO TO 75
         DO 70 I = K1,L
            P = A(I,1)/A(K,1)
            IK = I-K
            XL(K1,IK) = P
            DO 65 J = 2,NC
               A(I,J-1) = A(I,J)-P*A(K,J)
   65    CONTINUE
         A(I,NC) = ZERO
   70    CONTINUE
   75 CONTINUE
      IF (IJOB .EQ. 1) GO TO 9005
C                                  FORWARD SUBSTITUTION
   80 L = NLC
      DO 105 K = 1,N
         I = XL(K,NLC1)
         IF (I .EQ. K) GO TO 90
         DO 85 J = 1,M
            P = B(K,J)
            B(K,J) = B(I,J)
            B(I,J) = P
   85    CONTINUE
   90    IF (L .LT. N) L = L+1
         K1 = K+1
         IF (K1 .GT. L) GO TO 105
         DO 100 I = K1,L
            IK = I-K
            P = XL(K1,IK)
            DO 95 J = 1,M
               B(I,J) = B(I,J)-P*B(K,J)
   95       CONTINUE
  100    CONTINUE
  105 CONTINUE
C                                  BACKWARD SUBSTITUTION
      JBEG = NUC+NLC
      DO 125 J = 1,M
         L = 1
         K1 = N+1
         DO 120 I = 1,N
            K = K1-I
            P = B(K,J)
            IF (L .EQ. 1) GO TO 115
            DO 110 KK = 2,L
               IK = KK+K
               P = P-A(K,KK)*B(IK-1,J)
  110       CONTINUE
  115       B(K,J) = P/A(K,1)
            IF (L .LE. JBEG) L = L+1
  120    CONTINUE
  125 CONTINUE
      GO TO 9005
  135 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,'LEQT1B')
 9005 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE DEQT2B (A,N,NLC,NUC,IA,B,M,IB,IJOB,U,IU,XL,IER)
C
	implicit integer (i-n)
	implicit double precision (a-h)
	implicit double precision (o-z)
      DIMENSION          A(IA,*),U(IU,*),XL(N,*),B(IB,*)
      DOUBLE PRECISION   SUM
      DATA               ZERO/0.0/,ITMAX/50/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      NLCP1 = NLC+1
      NCC = NUC+NLCP1
      NLCP2 = NLC+2
      NMNUC = N-NUC
      NU2 = NCC+1
      NU3 = NCC+2
      IF (IJOB .EQ. 2) GO TO 15
C                                  SAVE MATRIX A
      DO 10 I = 1,N
         DO 5 J = 1,NCC
            U(I,J) = A(I,J)
    5    CONTINUE
   10 CONTINUE
C                                  FACTOR MATRIX A
cggu  CALL LEQT1B(U,N,NLC,NUC,IU,B,M,IB,1,XL,IER)
      CALL DEQT1B(U,N,NLC,NUC,IU,B,M,IB,1,XL,IER)
      IF (IER .NE. 0) GO TO 9000
      IF (IJOB .EQ. 1) GO TO 9005
C                                  SAVE THE RIGHT HAND SIDES
   15 DO 60 J = 1,M
         DO 20 I = 1,N
            U(I,NU2) = B(I,J)
   20    CONTINUE
C                                  OBTAIN A SOLUTION
cggu     CALL LEQT1B(U,N,NLC,NUC,IU,U(1,NU2),1,IU,2,XL,IER)
         CALL DEQT1B(U,N,NLC,NUC,IU,U(1,NU2),1,IU,2,XL,IER)
C                                  COMPUTE THE NORM OF THE SOLUTION
         XNORM = ZERO
         DO 25 I = 1,N
            XNORM = DMAX1(XNORM, DABS(U(I,NU2)))
   25    CONTINUE
         IF (XNORM .EQ. ZERO) GO TO 60
C                                  COMPUTE THE RESIDUALS
         DO 45 ITER = 1,ITMAX
            NC = NCC
            KK = 1
            DO 35 I = 1,N
               SUM = B(I,J)
               L = NLCP2-I
               IR = MAX0(L,1)
               IF (L .LE. 0) KK = KK+1
               K = KK
               DO 30 JJ = IR,NC
                  SUM = SUM-DBLE(A(I,JJ))*DBLE(U(K,NU2))
                  K = K+1
   30          CONTINUE
               U(I,NU3) = SUM
               IF (I .GE. NMNUC) NC = NC-1
   35       CONTINUE
cggu        CALL LEQT1B(U,N,NLC,NUC,IU,U(1,NU3),1,IU,2,XL,IER)
            CALL DEQT1B(U,N,NLC,NUC,IU,U(1,NU3),1,IU,2,XL,IER)
            DXNORM = ZERO
C                                  UPDATE THE SOLUTION
            DO 40 I = 1,N
               U(I,NU2) = U(I,NU2)+U(I,NU3)
               DXNORM = DMAX1(DXNORM, DABS(U(I,NU3)))
   40       CONTINUE
            IF (XNORM+DXNORM .EQ. XNORM) GO TO 50
   45    CONTINUE
         IER = 130
C                                  STORE THE SOLUTION
   50    DO 55 JK = 1,N
            B(JK,J) = U(JK,NU2)
   55    CONTINUE
         IF (IER .NE. 0) GO TO 9000
   60 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'LEQT2B')
 9005 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE DEQ1PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,IER)
C
	implicit integer (i-n)
	implicit double precision (a-h)
	implicit double precision (o-z)
      double precision               A(IA,*),B(IB,*),D1,D2
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0*IDGT
C                                  DECOMPOSITION OF MATRIX A INTO
C                                  L*L-TRANSPOSE
cggu  CALL LUDAPB(A,N,NC,IA,A,IA,D1,D2,IER)
      CALL DUDAPB(A,N,NC,IA,A,IA,D1,D2,IER)
      IF (IER .NE. 0) GO TO 9000
      DO 5 I = 1,M
C                                  SOLUTION OF AX = B
cggu     CALL LUELPB(A,B(1,I),N,NC,IA,B(1,I))
         CALL DUELPB(A,B(1,I),N,NC,IA,B(1,I))
    5 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'LEQ1PB')
 9005 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE DEQ2PB (A,N,NC,IA,B,IB,M,IDGT,D1,D2,WK,IER)
C
	implicit integer (i-n)
	implicit double precision (a-h)
	implicit double precision (o-z)
      double precision          A(IA,*),B(IB,*),WK(N,*),D1,D2
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      K = NC+2
      K1 = K+1
cggu  CALL LUDAPB(A,N,NC,IA,WK,N,D1,D2,IER)
      CALL DUDAPB(A,N,NC,IA,WK,N,D1,D2,IER)
      IF (IER .NE. 0) GO TO 9000
         DO 10 I = 1,M
cggu     CALL LUELPB(WK,B(1,I),N,NC,N,WK(1,K))
cggu     CALL LUREPB(A,N,NC,IA,WK,N,B(1,I),WK(1,K),IDGT,WK(1,K1),IER)
         CALL DUELPB(WK,B(1,I),N,NC,N,WK(1,K))
         CALL DUREPB(A,N,NC,IA,WK,N,B(1,I),WK(1,K),IDGT,WK(1,K1),IER)
            DO 5 J = 1,N
            B(J,I) = WK(J,K)
    5       CONTINUE
         IF (IER .NE. 0) GO TO 15
   10    CONTINUE
      GO TO 9005
   15 IER = 130
 9000 CONTINUE
      CALL UERTST(IER,'LEQ2PB')
 9005 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE DUDAPB (A,N,NC,IA,UL,IU,D1,D2,IER)
C
	implicit integer (i-n)
	implicit double precision (a-h)
	implicit double precision (o-z)
      double precision    ZERO,A(IA,*),UL(IU,*),D1,D2,ONE,SUM,
     +                   RN,FOUR,SIXTN,SIXTH
      DATA               ZERO/0.0/,FOUR/4.0/,SIXTN/16./,
     +                   SIXTH/.0625/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      RN = ONE/(N*SIXTN)
      D1 = ONE
      D2 = ZERO
      NCP1 = NC+1
      IF (NC .EQ. 0) GO TO 15
C                                  INITIALIZE ZERO ELEMENTS
      DO 10 I = 1,NC
         DO 5 J = I,NC
            K = NCP1-J
            UL(I,K) = ZERO
    5    CONTINUE
   10 CONTINUE
C                                  I IS ROW INDEX OF ELEMENT BEING
C                                  COMPUTED
   15 DO 60 I = 1,N
         IMNCP1 = I-NCP1
         I1 = MAX0(1,1-IMNCP1)
C                                  J IS COLUMN INDEX OF ELEMENT BEING
C                                  COMPUTED
         DO 60 J = I1,NCP1
C                                  L IS ROW INDEX OF PREVIOUSLY COMPUTED
C                                  VECTOR BEING USED TO COMPUTE INNER
C                                  PRODUCT
            L = IMNCP1+J
            I2 = NCP1-J
            SUM = A(I,J)
            JM1 = J-1
            IF (JM1) 30,30,20
   20       DO 25 K = 1,JM1
C                                  M IS COLUMN INDEX
               M = I2+K
               SUM = SUM-UL(I,K)*UL(L,M)
   25       CONTINUE
   30       IF (J .NE. NCP1) GO TO 55
            IF(A(I,J)+SUM*RN .LE. A(I,J))GO TO 65
            UL(I,J) = ONE/DSQRT(SUM)
C                                  UPDATE THE DETERMINANT
            D1 = D1*SUM
   35       IF (DABS(D1)-ONE) 45,45,40
   40       D1 = D1*SIXTH
            D2 = D2+FOUR
            GO TO 35
   45       IF (DABS(D1)-SIXTH) 50,50,60
   50       D1 = D1*SIXTN
            D2 = D2-FOUR
            GO TO 45
   55       UL(I,J) = SUM*UL(L,NCP1)
   60 CONTINUE
      GO TO 9005
   65 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,'LUDAPB')
 9005 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE DUELPB (UL,B,N,NC,IA,X)
C
	implicit integer (i-n)
	implicit double precision (a-h)
	implicit double precision (o-z)
      double precision               UL(IA,*),B(*),X(*),ZERO,SUM
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  SOLUTION LY = B
      NC1 = NC+1
      IW = 0
      L = 0
      DO 15 I = 1,N
         SUM = B(I)
         IF (NC .LE. 0) GO TO 10
         IF (IW .EQ. 0) GO TO 9
         L = L+1
         IF (L .GT. NC) L = NC
         K = NC1-L
         KL = I-L
         DO 5 J = K,NC
            SUM = SUM -X(KL) * UL(I,J)
            KL = KL+1
    5    CONTINUE
         GO TO 10
    9    IF (SUM .NE. ZERO) IW = 1
   10    X(I) = SUM*UL(I,NC1)
   15 CONTINUE
C                                  SOLUTION UX = Y
   20 X(N) = X(N)*UL(N,NC1)
      IF (N .LE. 1) GO TO 40
      N1 = N+1
      DO 35 I = 2,N
         K = N1-I
         SUM = X(K)
         IF (NC .LE. 0) GO TO 30
         KL = K+1
         K1 = MIN0(N,K+NC)
         L = 1
         DO 25 J = KL,K1
            SUM = SUM -X(J) * UL(J,NC1-L)
            L = L+1
   25    CONTINUE
   30    X(K) = SUM*UL(K,NC1)
   35 CONTINUE
   40 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE DUREPB (A,N,NC,IA,UL,IU,B,X,IDGT,RES,IER)
C
	implicit integer (i-n)
	implicit double precision (a-h)
	implicit double precision (o-z)
      double precision      A(IA,*),B(*),UL(IU,*),RES(*),X(*),XNORM,
     +                   ZERO,DXNORM
      DOUBLE PRECISION   SUM
      DATA               ITMAX/50/
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      NC1 = NC+1
      XNORM = ZERO
      DO 5 I = 1,N
         XNORM = DMAX1(XNORM,DABS(X(I)))
    5 CONTINUE
      IF (XNORM .NE. ZERO) GO TO 10
      IDGT = 50
      GO TO 9005
C                                  ITERATION LOOP
   10 DO 45 ITER = 1,ITMAX
         L = 0
         DO 30 I = 1,N
            K = MAX0(1,NC1-L)
            LK = MAX0(1,I-NC)
            L = I
            SUM = B(I)
            DO 15 J = K,NC1
               SUM = SUM+DBLE(-A(I,J))*DBLE(X(LK))
               LK = LK+1
   15       CONTINUE
            LL = MIN0(NC,N-I)
            IF (LL .LT. 1) GO TO 25
            DO 20 J = 1,LL
               SUM = SUM+DBLE(-A(I+J,NC1-J))*DBLE(X(LK))
               LK = LK+1
   20       CONTINUE
   25       RES(I) = SUM
   30    CONTINUE
cggu     CALL LUELPB(UL,RES,N,NC,IU,RES)
         CALL DUELPB(UL,RES,N,NC,IU,RES)
         DXNORM = ZERO
         XNORM = ZERO
         DO 35 I = 1,N
            X(I) = X(I)+RES(I)
            DXNORM = DMAX1(DXNORM,DABS(RES(I)))
            XNORM = DMAX1(XNORM,DABS(X(I)))
   35    CONTINUE
         IF (ITER .NE. 1) GO TO 40
         IF (DXNORM .GT. ZERO) IDGT = -DLOG10(DXNORM/XNORM)
   40    IF (XNORM + DXNORM .EQ. XNORM)GO TO 9005
   45 CONTINUE
C                                  ITERATIVE IMPROVEMENT FAILED
      IER = 129
 9000 CONTINUE
      CALL UERTST(IER,'LUREPB')
 9005 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE UERTST (IER,NAME)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
	character*6	name
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
	character*6	namset,nameq
	character*1	ieq
	data	namset,nameq,ieq /'UERSET','      ','='/
C                                  FIRST EXECUTABLE STATEMENT
      DATA               LEVEL/4/,IEQDF/0/
      IF (IER.GT.999) GO TO 25
      IF (IER.LT.-32) GO TO 55
      IF (IER.LE.128) GO TO 5
      IF (LEVEL.LT.1) GO TO 30
C                                  PRINT TERMINAL MESSAGE
      CALL UGETIO(1,NIN,IOUNIT)
      IF (IEQDF.EQ.1) WRITE(IOUNIT,35) IER,NAMEQ,IEQ,NAME
      IF (IEQDF.EQ.0) WRITE(IOUNIT,35) IER,NAME
      GO TO 30
    5 IF (IER.LE.64) GO TO 10
      IF (LEVEL.LT.2) GO TO 30
C                                  PRINT WARNING WITH FIX MESSAGE
      CALL UGETIO(1,NIN,IOUNIT)
      IF (IEQDF.EQ.1) WRITE(IOUNIT,40) IER,NAMEQ,IEQ,NAME
      IF (IEQDF.EQ.0) WRITE(IOUNIT,40) IER,NAME
      GO TO 30
   10 IF (IER.LE.32) GO TO 15
C                                  PRINT WARNING MESSAGE
      IF (LEVEL.LT.3) GO TO 30
      CALL UGETIO(1,NIN,IOUNIT)
      IF (IEQDF.EQ.1) WRITE(IOUNIT,45) IER,NAMEQ,IEQ,NAME
      IF (IEQDF.EQ.0) WRITE(IOUNIT,45) IER,NAME
      GO TO 30
   15 CONTINUE
C                                  CHECK FOR UERSET CALL
	if( name .ne. namset )  goto 25
c      DO 20 I=1,3
c         IF (NAME(I).NE.NAMSET(I)) GO TO 25
c   20 CONTINUE
      LEVOLD = LEVEL
      LEVEL = IER
      IER = LEVOLD
      IF (LEVEL.LT.0) LEVEL = 4
      IF (LEVEL.GT.4) LEVEL = 4
      GO TO 30
   25 CONTINUE
      IF (LEVEL.LT.4) GO TO 30
C                                  PRINT NON-DEFINED MESSAGE
      CALL UGETIO(1,NIN,IOUNIT)
      IF (IEQDF.EQ.1) WRITE(IOUNIT,50) IER,NAMEQ,IEQ,NAME
      IF (IEQDF.EQ.0) WRITE(IOUNIT,50) IER,NAME
   30 IEQDF = 0
      RETURN
   35 FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3,
     +       20H) FROM IMSL ROUTINE ,a6,A1,a6)
   40 FORMAT(36H *** WARNING WITH FIX ERROR  (IER = ,I3,
     +       20H) FROM IMSL ROUTINE ,a6,A1,a6)
   45 FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3,
     +       20H) FROM IMSL ROUTINE ,a6,A1,a6)
   50 FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5,
     +       20H) FROM IMSL ROUTINE ,a6,A1,a6)
C                                  SAVE P FOR P = R CASE
C                                    P IS THE PAGE NAME
C                                    R IS THE ROUTINE NAME
   55 IEQDF = 1
	nameq = name
c      DO 60 I=1,3
c   60 NAMEQ(I) = NAME(I)
   65 RETURN
      END
C
c********************************************************************
C
      SUBROUTINE UGETIO (IOPT,NIN,NOUT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,NIN,NOUT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NIND,NOUTD
      DATA               NIND/5/,NOUTD/6/
C                                  FIRST EXECUTABLE STATEMENT
      IF (IOPT.EQ.3) GO TO 10
      IF (IOPT.EQ.2) GO TO 5
      IF (IOPT.NE.1) GO TO 9005
      NIN = NIND
      NOUT = NOUTD
      GO TO 9005
    5 NIND = NIN
      GO TO 9005
   10 NOUTD = NOUT
 9005 RETURN
      END
C
c********************************************************************
C
	SUBROUTINE GELB(R,A,M,N,MUD,MLD,EPS,IER)
C
	DIMENSION R(*),A(*)
	J=0
	JJ=0
	IF(MLD)47,1,1
    1	IF(MUD)47,2,2
    2	MC=1+MLD+MUD
	IF(MC+1-M-M)3,3,47
    3	IF(MC-M)5,5,4
    4	MC=M
    5	MU=MC-MUD-1
	ML=MC-MLD-1
	MR=M-ML
	MZ=(MU*(MU+1))/2
	MA=M*MC-(ML*(ML+1))/2
	NM=N*M
	IER=0
	PIV=0.
	IF(MLD)14,14,6
    6	JJ=MA
	J=MA-MZ
	KST=J
	DO 9 K=1,KST
	TB=A(J)
	A(JJ)=TB
	TB=ABS(TB)
	IF(TB-PIV)8,8,7
    7	PIV=TB
    8	J=J-1
    9	JJ=JJ-1
	IF(MZ)14,14,10
   10	JJ=1
	J=1+MZ
	IC=1+MUD
	DO 13 I=1,MU
	DO 12 K=1,MC
	A(JJ)=0.
	IF(K-IC)11,11,12
   11	A(JJ)=A(J)
	J=J+1
   12	JJ=JJ+1
   13	IC=IC+1
   14	TOL=EPS*PIV
	KST=1
	IDST=MC
	IC=MC-1
	DO 38 K=1,M
	IF(K-MR-1)16,16,15
   15	IDST=IDST-1
   16	ID=IDST
	ILR=K+MLD
	IF(ILR-M)18,18,17
   17	ILR=M
   18	II=KST
	PIV=0.
	DO 22 I=K,ILR
	TB=ABS(A(II))
	IF(TB-PIV)20,20,19
   19	PIV=TB
	J=I
	JJ=II
   20	IF(I-MR)22,22,21
   21	ID=ID-1
   22	II=II+ID
	IF(PIV)47,47,23
   23	IF(IER)26,24,26
   24	IF(PIV-TOL)25,25,26
   25	IER=K-1
   26	PIV=1./A(JJ)
	ID=J-K
	DO 27 I=K,NM,M
	II=I+ID
	TB=PIV*R(II)
	R(II)=R(I)
   27	R(I)=TB
	II=KST
	J=JJ+IC
	DO 28 I=JJ,J
	TB=PIV*A(I)
	A(I)=A(II)
	A(II)=TB
   28	II=II+1
	IF(K-ILR)29,34,34
   29	ID=KST
	II=K+1
	MU=KST+1
	MZ=KST+IC
	DO 33 I=II,ILR
	ID=ID+MC
	JJ=I-MR-1
	IF(JJ)31,31,30
   30	ID=ID-JJ
   31	PIV=-A(ID)
	J=ID+1
	DO 32 JJ=MU,MZ
	A(J-1)=A(J)+PIV*A(JJ)
   32	J=J+1
	A(J-1)=0.
	J=K
	DO 33 JJ=I,NM,M
	R(JJ)=R(JJ)+PIV*R(J)
   33	J=J+M
   34	KST=KST+MC
	IF(ILR-MR)36,35,35
   35	IC=IC-1
   36	ID=K-MR
	IF(ID)38,38,37
   37	KST=KST-ID
   38	CONTINUE
	IF(MC-1)46,46,39
   39	IC=2
	KST=MA+ML-MC+2
	II=M
	DO 45 I=2,M
	KST=KST-MC
	II=II-1
	J=II-MR
	IF(J)41,41,40
   40	KST=KST+J
   41	DO 43 J=II,NM,M
	TB=R(J)
	MZ=KST+IC-2
	ID=J
	DO 42 JJ=KST,MZ
	ID=ID+1
   42	TB=TB-A(JJ)*R(ID)
   43	R(J)=TB
	IF(IC-MC)44,45,45
   44	IC=IC+1
   45	CONTINUE
   46	RETURN
   47	IER=-1
	RETURN
	END
C
c********************************************************************
C
	SUBROUTINE DGELB(R,A,M,N,MUD,MLD,EPS,IER)
C
	DIMENSION R(*),A(*)
	DOUBLE PRECISION R,A,PIV,TB,TOL
	J=0
	JJ=0
	IF(MLD)47,1,1
    1	IF(MUD)47,2,2
    2	MC=1+MLD+MUD
	IF(MC+1-M-M)3,3,47
    3	IF(MC-M)5,5,4
    4	MC=M
    5	MU=MC-MUD-1
	ML=MC-MLD-1
	MR=M-ML
	MZ=(MU*(MU+1))/2
	MA=M*MC-(ML*(ML+1))/2
	NM=N*M
	IER=0
	PIV=0.D0
	IF(MLD)14,14,6
    6	JJ=MA
	J=MA-MZ
	KST=J
	DO 9 K=1,KST
	TB=A(J)
	A(JJ)=TB
	TB=DABS(TB)
	IF(TB-PIV)8,8,7
    7	PIV=TB
    8	J=J-1
    9	JJ=JJ-1
	IF(MZ)14,14,10
   10	JJ=1
	J=1+MZ
	IC=1+MUD
	DO 13 I=1,MU
	DO 12 K=1,MC
	A(JJ)=0.D0
	IF(K-IC)11,11,12
   11	A(JJ)=A(J)
	J=J+1
   12	JJ=JJ+1
   13	IC=IC+1
   14	TOL=EPS*PIV
	KST=1
	IDST=MC
	IC=MC-1
	DO 38 K=1,M
	IF(K-MR-1)16,16,15
   15	IDST=IDST-1
   16	ID=IDST
	ILR=K+MLD
	IF(ILR-M)18,18,17
   17	ILR=M
   18	II=KST
	PIV=0.D0
	DO 22 I=K,ILR
	TB=DABS(A(II))
	IF(TB-PIV)20,20,19
   19	PIV=TB
	J=I
	JJ=II
   20	IF(I-MR)22,22,21
   21	ID=ID-1
   22	II=II+ID
	IF(PIV)47,47,23
   23	IF(IER)26,24,26
   24	IF(PIV-TOL)25,25,26
   25	IER=K-1
   26	PIV=1.D0/A(JJ)
	ID=J-K
	DO 27 I=K,NM,M
	II=I+ID
	TB=PIV*R(II)
	R(II)=R(I)
   27	R(I)=TB
	II=KST
	J=JJ+IC
	DO 28 I=JJ,J
	TB=PIV*A(I)
	A(I)=A(II)
	A(II)=TB
   28	II=II+1
	IF(K-ILR)29,34,34
   29	ID=KST
	II=K+1
	MU=KST+1
	MZ=KST+IC
	DO 33 I=II,ILR
	ID=ID+MC
	JJ=I-MR-1
	IF(JJ)31,31,30
   30	ID=ID-JJ
   31	PIV=-A(ID)
	J=ID+1
	DO 32 JJ=MU,MZ
	A(J-1)=A(J)+PIV*A(JJ)
   32	J=J+1
	A(J-1)=0.D0
	J=K
	DO 33 JJ=I,NM,M
	R(JJ)=R(JJ)+PIV*R(J)
   33	J=J+M
   34	KST=KST+MC
	IF(ILR-MR)36,35,35
   35	IC=IC-1
   36	ID=K-MR
	IF(ID)38,38,37
   37	KST=KST-ID
   38	CONTINUE
	IF(MC-1)46,46,39
   39	IC=2
	KST=MA+ML-MC+2
	II=M
	DO 45 I=2,M
	KST=KST-MC
	II=II-1
	J=II-MR
	IF(J)41,41,40
   40	KST=KST+J
   41	DO 43 J=II,NM,M
	TB=R(J)
	MZ=KST+IC-2
	ID=J
	DO 42 JJ=KST,MZ
	ID=ID+1
   42	TB=TB-A(JJ)*R(ID)
   43	R(J)=TB
	IF(IC-MC)44,45,45
   44	IC=IC+1
   45	CONTINUE
   46	RETURN
   47	IER=-1
	RETURN
	END
c
c******************************************************************
c
	SUBROUTINE MCHB(R,A,M,N,MUD,IOP,EPS,IER)
c
        DIMENSION R(*),A(*)
	DOUBLE PRECISION TOL,SUM,PIV
	PIV=0.
	IF(IABS(IOP)-3)1,1,43
    1	IF(MUD)43,2,2
    2	MC=MUD+1
	IF(M-MC)43,3,3
    3	MR=M-MUD
	IER=0
	IF(IOP)24,4,4
    4	IEND=0
	LLDST=MUD
	DO 23 K=1,M
	IST=IEND+1
	IEND=IST+MUD
	J=K-MR
	IF(J)6,6,5
    5	IEND=IEND-J
    6	IF(J-1)8,8,7
    7	LLDST=LLDST-1
    8	LMAX=MUD
	J=MC-K
	IF(J)10,10,9
    9	LMAX=LMAX-J
   10	ID=0
	TOL=A(IST)*EPS
	DO 23 I=IST,IEND
	SUM=0.D0
	IF(LMAX)14,14,11
   11	LL=IST
	LLD=LLDST
	DO 13 L=1,LMAX
	LL=LL-LLD
	LLL=LL+ID
	SUM=SUM+A(LL)*A(LLL)
	IF(LLD-MUD)12,13,13
   12	LLD=LLD+1
   13	CONTINUE
   14	SUM=DBLE(A(I))-SUM
	IF(I-IST)15,15,20
   15	IF(SUM)43,43,16
   16	IF(SUM-TOL)17,17,19
   17	IF(IER)18,18,19
   18	IER=K-1
   19	PIV=DSQRT(SUM)
	A(I)=PIV
	PIV=1.D0/PIV
	GO TO 21
   20	A(I)=SUM*PIV
   21	ID=ID+1
	IF(ID-J)23,23,22
   22	LMAX=LMAX-1
   23	CONTINUE
	IF(IOP)24,44,24
   24	ID=N*M
	IEND=IABS(IOP)-2
	IF(IEND)25,35,25
   25	IST=1
	LMAX=0
	J=-MR
	LLDST=MUD
	DO 34 K=1,M
	PIV=A(IST)
	IF(PIV)26,43,26
   26	PIV=1.D0/PIV
	DO 30 I=K,ID,M
	SUM=0.D0
	IF(LMAX)30,30,27
   27	LL=IST
	LLL=I
	LLD=LLDST
	DO 29 L=1,LMAX
	LL=LL-LLD
	LLL=LLL-1
	SUM=SUM+A(LL)*R(LLL)
	IF(LLD-MUD)28,29,29
   28	LLD=LLD+1
   29	CONTINUE
   30	R(I)=PIV*(DBLE(R(I))-SUM)
	IF(MC-K)32,32,31
   31	LMAX=K
   32	IST=IST+MC
	J=J+1
	IF(J)34,34,33
   33	IST=IST-J
	LLDST=LLDST-1
   34	CONTINUE
	IF(IEND)35,35,44
   35	IST=M+(MUD*(M+M-MC))/2+1
	LMAX=0
	K=M
   36	IEND=IST-1
	IST=IEND-LMAX
	PIV=A(IST)
	IF(PIV)37,43,37
   37	PIV=1.D0/PIV
	L=IST+1
	DO 40 I=K,ID,M
	SUM=0.D0
	IF(LMAX)40,40,38
   38	LLL=I
	DO 39 LL=L,IEND
	LLL=LLL+1
   39	SUM=SUM+A(LL)*R(LLL)
   40	R(I)=PIV*(DBLE(R(I))-SUM)
	IF(K-MR)42,42,41
   41	LMAX=LMAX+1
   42	K=K-1
	IF(K)44,44,36
   43	IER=-1
   44	RETURN
	END
c
c******************************************************************
c
	SUBROUTINE DMCHB(R,A,M,N,MUD,IOP,EPS,IER)
c
        DIMENSION R(*),A(*)
	DOUBLE PRECISION TOL,SUM,PIV,R,A
	PIV=0.
	IF(IABS(IOP)-3)1,1,43
    1	IF(MUD)43,2,2
    2	MC=MUD+1
	IF(M-MC)43,3,3
    3	MR=M-MUD
	IER=0
	IF(IOP)24,4,4
    4	IEND=0
	LLDST=MUD
	DO 23 K=1,M
	IST=IEND+1
	IEND=IST+MUD
	J=K-MR
	IF(J)6,6,5
    5	IEND=IEND-J
    6	IF(J-1)8,8,7
    7	LLDST=LLDST-1
    8	LMAX=MUD
	J=MC-K
	IF(J)10,10,9
    9	LMAX=LMAX-J
   10	ID=0
	TOL=A(IST)*EPS
	DO 23 I=IST,IEND
	SUM=0.D0
	IF(LMAX)14,14,11
   11	LL=IST
	LLD=LLDST
	DO 13 L=1,LMAX
	LL=LL-LLD
	LLL=LL+ID
	SUM=SUM+A(LL)*A(LLL)
	IF(LLD-MUD)12,13,13
   12	LLD=LLD+1
   13	CONTINUE
   14	SUM=A(I)-SUM
	IF(I-IST)15,15,20
   15	IF(SUM)43,43,16
   16	IF(SUM-TOL)17,17,19
   17	IF(IER)18,18,19
   18	IER=K-1
   19	PIV=DSQRT(SUM)
	A(I)=PIV
	PIV=1.D0/PIV
	GO TO 21
   20	A(I)=SUM*PIV
   21	ID=ID+1
	IF(ID-J)23,23,22
   22	LMAX=LMAX-1
   23	CONTINUE
	IF(IOP)24,44,24
   24	ID=N*M
	IEND=IABS(IOP)-2
	IF(IEND)25,35,25
   25	IST=1
	LMAX=0
	J=-MR
	LLDST=MUD
	DO 34 K=1,M
	PIV=A(IST)
	IF(PIV)26,43,26
   26	PIV=1.D0/PIV
	DO 30 I=K,ID,M
	SUM=0.D0
	IF(LMAX)30,30,27
   27	LL=IST
	LLL=I
	LLD=LLDST
	DO 29 L=1,LMAX
	LL=LL-LLD
	LLL=LLL-1
	SUM=SUM+A(LL)*R(LLL)
	IF(LLD-MUD)28,29,29
   28	LLD=LLD+1
   29	CONTINUE
   30	R(I)=PIV*(R(I)-SUM)
	IF(MC-K)32,32,31
   31	LMAX=K
   32	IST=IST+MC
	J=J+1
	IF(J)34,34,33
   33	IST=IST-J
	LLDST=LLDST-1
   34	CONTINUE
	IF(IEND)35,35,44
   35	IST=M+(MUD*(M+M-MC))/2+1
	LMAX=0
	K=M
   36	IEND=IST-1
	IST=IEND-LMAX
	PIV=A(IST)
	IF(PIV)37,43,37
   37	PIV=1.D0/PIV
	L=IST+1
	DO 40 I=K,ID,M
	SUM=0.D0
	IF(LMAX)40,40,38
   38	LLL=I
	DO 39 LL=L,IEND
	LLL=LLL+1
   39	SUM=SUM+A(LL)*R(LLL)
   40	R(I)=PIV*(R(I)-SUM)
	IF(K-MR)42,42,41
   41	LMAX=LMAX+1
   42	K=K-1
	IF(K)44,44,36
   43	IER=-1
   44	RETURN
	END
c
c*************************************************************************
c
        subroutine loctst(i,j,n,m)
c
c control indices
c
c (i,j)   position of element in square matrix (row,column)
c n       dimension of square matrix
c m       band width of square matrix
c
        implicit none
c
        integer i,j,n,m
c
        if(i.lt.1.or.i.gt.n.or.j.lt.1.or.j.gt.n) then
          write(6,*) '(i,j) out of matrix'
          write(6,*) 'i,j,n : ',i,j,n
          stop 'error stop : loctst'
        end if
c
        if(iabs(i-j).gt.m) then
          write(6,*) '(i,j) out of band'
          write(6,*) 'i,j,|i-j|,m : ',i,j,iabs(i-j),m
          stop 'error stop : loctst'
        end if
c
        return
        end
c
c*************************************************************************
c
        function locmy(i,j,n,m)
c
c access my routines
c
c (i,j)   position of element in square matrix (row,column)
c n       dimension of square matrix
c m       band width of square matrix
c locmy   position of element in band matrix
c
        implicit none
c
	integer locmy
        integer i,j,n,m
c
        if(i.eq.j) then
          locmy=i
        else if(i.gt.j) then
          locmy=n+j*(m-1)+i-m-1
        else!if(i.lt.j) then
          locmy=n+n*m+i*(m-1)+j-m-1
        end if
c
        return
        end
c
c*************************************************************************
c
        function locimm(i,j,n,m)
c
c access imsl routines
c
c formula if mlo!=mup : loc = n*(mlo+j-i)+i
c
c (i,j)   position of element in square matrix (row,column)
c n       dimension of square matrix
c m       band width of square matrix
c locimm  position of element in band matrix
c
        implicit none
c
	integer locimm
        integer i,j,n,m
c
        locimm = n*(m+j-i)+i
c
        return
        end
c
c*************************************************************************
c
        function loccer(i,j,n,m)
c
c access kernlib routines
c
c formula if mlo!=mup : loc = n*(j-max(i-mlo,1))+i
c
c (i,j)   position of element in square matrix (row,column)
c n       dimension of square matrix
c m       band width of square matrix
c loccer  position of element in band matrix
c
        implicit none
c
	integer loccer
        integer i,j,n,m
c
        loccer = n*(j-max(i-m,1))+i
c
        return
        end
c
c*************************************************************************
c
        function locssp(i,j,n,m)
c
c access ssp routines
c
c (i,j)   position of element in square matrix (row,column)
c n       dimension of square matrix
c m       band width of square matrix
c locssp  position of element in band matrix
c
        implicit none
c
	integer locssp
        integer i,j,n,m
c
	if(i-j.gt.m.or.j-i.gt.m) then
	  locssp = 0
	else if(n.le.m) then	!this is for a full matrix
	  locssp = n*(i-1) + j
        else if(i.lt.m) then
          locssp = 2*m*i - m + j - m*(m+1)/2 + (m-i)*(m-i+1)/2
        else if(i.gt.n-m+1) then
          locssp = 2*m*i - m + j - m*(m+1)/2
     +                - (i-(n-m+1))*(i-(n-m))/2
        else
          locssp = 2*m*i - m + j - m*(m+1)/2
        end if
c
        return
        end
c
c*************************************************************************
c
        function locsps(i,j,n,m)
c
c access ssp routines (symmetric compressed storage mode)
c
c only main and upper diagonals - if an element in the lower
c ...diagonals is referenced, 0 is returned
c
c (i,j)   position of element in square matrix (row,column)
c n       dimension of square matrix
c m       band width of square matrix
c locsps  position of element in band matrix
c
c original formula : locsps = (1+m)*(i-1)+abs(j-i)+1
c
        implicit none
c
	integer locsps
        integer i,j,n,m
c
	if(i.gt.j.or.j-i.gt.m) then
	  locsps = 0
        else if(i.gt.n-m+1) then
          locsps = m*(i-1)+j
     +                - (i-(n-m+1))*(i-(n-m))/2
        else
          locsps = m*(i-1)+j
        end if
c
        return
        end

c*************************************************************************

