!*==runrc.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
      SUBROUTINE RUNRC(N,Rhs,Sol,Ipar,Fpar,Wk,Guess,A,Ja,Ia,Au,Jau,Ju,
     &                 SOLVER)
      IMPLICIT NONE
      INTEGER N , Ipar(16) , Ia(N+1) , Ja(*) , Ju(*) , Jau(*)
      REAL*8 Fpar(16) , Rhs(N) , Sol(N) , Guess(N) , Wk(*) , A(*) , 
     &       Au(*)
      EXTERNAL SOLVER
!-----------------------------------------------------------------------
!     the actual tester. It starts the iterative linear system solvers
!     with a initial guess suppied by the user.
!
!     The structure {au, jau, ju} is assumed to have the output from
!     the ILU* routines in ilut.f.
!
!-----------------------------------------------------------------------
!     local variables
!
      INTEGER i , iou , its
      REAL*8 res , DNRM2
!     real dtime, dt(2), time
!     external dtime
      EXTERNAL DNRM2
      SAVE its , res

!
!     ipar(2) can be 0, 1, 2, please don't use 3
!
      IF ( Ipar(2)>2 ) THEN
         PRINT * , 'I can not do both left and right preconditioning.'
         RETURN
      ENDIF
!
!     normal execution
!
      its = 0
      res = 0.0D0
!
      DO i = 1 , N
         Sol(i) = Guess(i)
      ENDDO
!
      iou = 6
      Ipar(1) = 0
      DO
!     time = dtime(dt)
         CALL SOLVER(N,Rhs,Sol,Ipar,Fpar,Wk)
!
!     output the residuals
!
!         write (iou, *) its, real(res)
         IF ( Ipar(7)/=its ) its = Ipar(7)
         res = Fpar(5)
!
	!write(6,*) 'ipar(1) = ',ipar(1),its

         IF ( Ipar(1)==1 ) THEN
            CALL AMUX(N,Wk(Ipar(8)),Wk(Ipar(9)),A,Ja,Ia)
            CYCLE
         ELSEIF ( Ipar(1)==2 ) THEN
            CALL ATMUX(N,Wk(Ipar(8)),Wk(Ipar(9)),A,Ja,Ia)
            CYCLE
         ELSEIF ( Ipar(1)==3 .OR. Ipar(1)==5 ) THEN
            CALL LUSOL(N,Wk(Ipar(8)),Wk(Ipar(9)),Au,Jau,Ju)
            CYCLE
         ELSEIF ( Ipar(1)==4 .OR. Ipar(1)==6 ) THEN
            CALL LUTSOL(N,Wk(Ipar(8)),Wk(Ipar(9)),Au,Jau,Ju)
            CYCLE
         ELSEIF ( Ipar(1)<=0 ) THEN
            IF ( Ipar(1)==0 ) THEN
!            print *, 'Iterative sovler has satisfied convergence test.'
            ELSEIF ( Ipar(1)==-1 ) THEN
               PRINT * , 'Iterative solver has iterated too many times.'
            ELSEIF ( Ipar(1)==-2 ) THEN
!            print *, 'Iterative solver was not given enough work space.'
!            print *, 'The work space should at least have ', ipar(4), &
!     &           ' elements.'
            ELSEIF ( Ipar(1)/=-3 ) THEN
!            print *, 'Iterative sovler is facing a break-down.'
               PRINT * , 'Iterative solver terminated. code =' , Ipar(1)
            ENDIF
         ENDIF
!     time = dtime(dt)
!      write (iou, *) ipar(7), real(fpar(6))
!      write (iou, *) '# retrun code =', ipar(1), &
!     &     '	convergence rate =', fpar(7)
!     write (iou, *) '# total execution time (sec)', time
!
!     check the error
!
         CALL AMUX(N,Sol,Wk,A,Ja,Ia)
         DO i = 1 , N
            Wk(N+i) = Sol(i) - 1.0D0
            Wk(i) = Wk(i) - Rhs(i)
         ENDDO
!      write (iou, *) '# the actual residual norm is', dnrm2(n,wk,1)
!      write (iou, *) '# the error norm is', dnrm2(n,wk(1+n),1)
!
         IF ( iou/=6 ) CLOSE (iou)
         EXIT
      ENDDO
      END SUBROUTINE RUNRC
!*==distdot.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!-----end-of-runrc
!-----------------------------------------------------------------------
      FUNCTION DISTDOT(N,X,Ix,Y,Iy)
      INTEGER N , Ix , Iy
      REAL*8 DISTDOT , X(*) , Y(*) , DDOT
      EXTERNAL DDOT
      DISTDOT = DDOT(N,X,Ix,Y,Iy)
      END FUNCTION DISTDOT
!*==afun.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!-----end-of-distdot
!-----------------------------------------------------------------------
!
      FUNCTION AFUN(X,Y,Z)
      REAL*8 AFUN , X , Y , Z
      AFUN = -1.0D0
      END FUNCTION AFUN
!*==bfun.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      FUNCTION BFUN(X,Y,Z)
      REAL*8 BFUN , X , Y , Z
      BFUN = -1.0D0
      END FUNCTION BFUN
!*==cfun.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      FUNCTION CFUN(X,Y,Z)
      REAL*8 CFUN , X , Y , Z
      CFUN = -1.0D0
      END FUNCTION CFUN
!*==dfun.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      FUNCTION DFUN(X,Y,Z)
      REAL*8 DFUN , X , Y , Z , GAMmax , GAMmay , ALPha
      COMMON /FUNC  / GAMmax , GAMmay , ALPha
      DFUN = GAMmax*EXP(X*Y)
      END FUNCTION DFUN
!*==efun.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      FUNCTION EFUN(X,Y,Z)
      REAL*8 EFUN , X , Y , Z , GAMmax , GAMmay , ALPha
      COMMON /FUNC  / GAMmax , GAMmay , ALPha
      EFUN = GAMmay*EXP(-X*Y)
      END FUNCTION EFUN
!*==ffun.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      FUNCTION FFUN(X,Y,Z)
      REAL*8 FFUN , X , Y , Z
      FFUN = 0.0D0
      END FUNCTION FFUN
!*==gfun.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      FUNCTION GFUN(X,Y,Z)
      REAL*8 GFUN , X , Y , Z , GAMmax , GAMmay , ALPha
      COMMON /FUNC  / GAMmax , GAMmay , ALPha
      GFUN = ALPha
      END FUNCTION GFUN
!*==hfun.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      FUNCTION HFUN(X,Y,Z)
      REAL*8 HFUN , X , Y , Z , GAMmax , GAMmay , ALPha
      COMMON /FUNC  / GAMmax , GAMmay , ALPha
      HFUN = ALPha*SIN(GAMmax*X+GAMmay*Y-Z)
      END FUNCTION HFUN
!*==betfun.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
!
      FUNCTION BETFUN(Side,X,Y,Z)
      REAL*8 BETFUN , X , Y , Z
      CHARACTER*2 Side
      BETFUN = 1.0
      END FUNCTION BETFUN
!*==gamfun.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      FUNCTION GAMFUN(Side,X,Y,Z)
      REAL*8 GAMFUN , X , Y , Z
      CHARACTER*2 Side
      IF ( Side=='x2' ) THEN
         GAMFUN = 5.0
      ELSEIF ( Side=='y1' ) THEN
         GAMFUN = 2.0
      ELSEIF ( Side=='y2' ) THEN
         GAMFUN = 7.0
      ELSE
         GAMFUN = 0.0
      ENDIF
      END FUNCTION GAMFUN
!*==afunbl.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!-----------------------------------------------------------------------
!     functions for the block PDE's
!-----------------------------------------------------------------------
      SUBROUTINE AFUNBL(Nfree,X,Y,Z,Coeff)
      END SUBROUTINE AFUNBL
!*==bfunbl.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      SUBROUTINE BFUNBL(Nfree,X,Y,Z,Coeff)
      END SUBROUTINE BFUNBL
!*==cfunbl.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      SUBROUTINE CFUNBL(Nfree,X,Y,Z,Coeff)
!
      END SUBROUTINE CFUNBL
!*==dfunbl.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      SUBROUTINE DFUNBL(Nfree,X,Y,Z,Coeff)
!
      END SUBROUTINE DFUNBL
!*==efunbl.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      SUBROUTINE EFUNBL(Nfree,X,Y,Z,Coeff)
      END SUBROUTINE EFUNBL
!*==ffunbl.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      SUBROUTINE FFUNBL(Nfree,X,Y,Z,Coeff)
      END SUBROUTINE FFUNBL
!*==gfunbl.spg  processed by SPAG 6.50Rc at 22:35 on  9 Feb 2009
!
      SUBROUTINE GFUNBL(Nfree,X,Y,Z,Coeff)
      END SUBROUTINE GFUNBL
