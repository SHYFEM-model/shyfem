      REAL(RLEN)FUNCTION QGAUS_EXP(A,B,LABDA,FN)
        USE global_mem, ONLY:RLEN
        USE bennut_interface,ONLY:BESS_EXP

        REAL(RLEN) ::	A,B,LABDA(2)         ! Specification

        INTERFACE                            ! Specification
           REAL(RLEN) FUNCTION FN(X)         ! Specification
            USE global_mem, ONLY:RLEN
            REAL(RLEN),INTENT(IN)   ::X      ! Specification
          END FUNCTION                       ! Specification
        END  INTERFACE                       ! Specification    

        REAL(RLEN) ::	X(5),W(5)
        REAL(RLEN) ::	XM,XR,SS,DX


        X(1)= .1488743389D0
        X(2)= .4333953941D0
        X(3)= .6794095682D0
        X(4)= .8650633666D0
        X(5)= .9739065285D0
        W(1)= .2955242247D0
        W(2)= .2692667193D0
        W(3)= .2190863625D0
        W(4)= .1494513491D0
        W(5)= .0666713443D0

        XM=0.5*(B+A)
        XR=0.5*(B-A)
        SS=0.0D+00
        DO J=1,5
          DX=XR*X(J)
          SS=SS+W(J)*(BESS_EXP(XM+DX,LABDA,FN)+BESS_EXP(XM-DX,LABDA,FN))
        ENDDO
      SS=XR*SS
      QGAUS_EXP=SS
      RETURN
      END


