!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! INTERFACE
!   ModuleBenNutInterface
!
! FILE
!   ModuleBenNutInterface
!
! DESCRIPTION
!   Definition of Explicit Interfaces for benthic nutrient model
!  
!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! AUTHORS
!   mfstep/ERSEM team
!
! CHANGE_LOG
!   ---
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
      MODULE bennut_interface

      USE global_mem, ONLY:RLEN

      INTERFACE

      REAL(RLEN) FUNCTION bess_exp(X,labda,FN)
      USE global_mem, ONLY:RLEN
      implicit none
      REAL(RLEN),intent(IN)    ::X
      REAL(RLEN),intent(IN)    ::labda(2)
      INTERFACE
         REAL(RLEN) FUNCTION FN(X)
         USE global_mem, ONLY:RLEN
         REAL(RLEN),INTENT(IN)   ::X
      END FUNCTION
      end INTERFACE
      end FUNCTION bess_exp
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION BESSI0(X)
      USE global_mem, ONLY:RLEN
      implicit none
      REAL(RLEN),intent(IN) ::x
      end FUNCTION BESSI0
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION BESSI1(X)
      USE global_mem, ONLY:RLEN
      implicit none
      REAL(RLEN),intent(IN) ::x
      end FUNCTION BESSI1
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION BESSK0(X)
      USE global_mem, ONLY:RLEN
      implicit none
      REAL(RLEN),intent(IN) ::x
      end FUNCTION BESSK0
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION BESSK1(X)
      USE global_mem, ONLY:RLEN
      implicit none
      REAL(RLEN),intent(IN) ::x
      end FUNCTION BESSK1
      end INTERFACE

      INTERFACE

      SUBROUTINE calcadd(irow,k,C,nn,fc)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::irow
      integer,intent(IN) ::nn
      integer,intent(IN) ::k
      REAL(RLEN),intent(INOUT) ::c(nn,nn)
      REAL(RLEN),intent(IN) ::fc
      end       SUBROUTINE calcadd
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION CalculateFromCondition(irow,C,Y,nn)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::irow
      integer,intent(IN) ::nn
      REAL(RLEN),intent(IN) ::c(nn,nn)
      REAL(RLEN),intent(IN) ::y(nn)
      end FUNCTION CalculateFromCondition
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION calculate_equation(mode,xinput,coeffs,&
                                                      b,factor,nn)
      USE global_mem, ONLY:RLEN
      USE bennut_type, ONLY:ty_coeff
      implicit none
      integer,intent(IN) ::mode
      integer,intent(IN) ::nn
!     integer,intent(IN) ::input
      type (ty_coeff),intent(IN) ::coeffs(nn)
      REAL(RLEN),intent(IN) ::xinput
      REAL(RLEN),intent(IN) ::b
      REAL(RLEN),intent(IN) ::factor(nn)
      end FUNCTION calculate_equation
      end INTERFACE

      INTERFACE

      SUBROUTINE CalculateLayer(NUTR,mode,x,layernr,xb)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::NUTR
      integer,intent(IN) ::mode
      integer,intent(INOUT) ::layernr
      REAL(RLEN),intent(IN) ::x
      REAL(RLEN),intent(INOUT) ::xb
      end       SUBROUTINE CalculateLayer
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION calculate_one_term(mode,option,xinput, &
      coeff,b)
      USE global_mem, ONLY:RLEN
      USE bennut_type
      implicit none
      integer,intent(IN) ::mode
      integer,intent(IN) ::option
      type (ty_coeff),intent(IN) ::coeff
      REAL(RLEN),intent(IN) ::xinput
      REAL(RLEN),intent(IN) ::b
      end FUNCTION calculate_one_term
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION CalculateTau(sMI,xdiffMI,ptMI,DXm)
      USE global_mem, ONLY:RLEN
      implicit none
      REAL(RLEN),intent(IN) ::smi
      REAL(RLEN),intent(IN) ::xdiffmi
      REAL(RLEN),intent(IN) ::ptmi
      REAL(RLEN),intent(IN) ::dxm
      end FUNCTION CalculateTau
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION chebev(a,b,c,m,x)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::m
      REAL(RLEN),intent(IN) ::a
      REAL(RLEN),intent(IN) ::b
      REAL(RLEN),intent(IN) ::x
      REAL(RLEN),intent(IN) ::c(m)
      end FUNCTION chebev
      end INTERFACE

      INTERFACE

      SUBROUTINE CompleteSet(NUTR,mode,option,input,x1,y1)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::NUTR
      integer,intent(IN) ::mode
      integer,intent(IN) ::input
      integer,intent(IN) ::option
      REAL(RLEN),intent(IN) ::x1
      REAL(RLEN),intent(IN) ::y1
      end SUBROUTINE CompleteSet
      end INTERFACE

      INTERFACE

      INTEGER FUNCTION InitializeSet(NUTR,option,input)
      implicit none
      integer,intent(IN) ::NUTR
      integer,intent(IN) ::input
      integer,intent(IN) ::option
      end FUNCTION InitializeSet
      end INTERFACE

      INTERFACE

      SUBROUTINE DefineSet(NUTR,mode,option,input,x1,y1)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::NUTR
      integer,intent(IN) ::mode
      integer,intent(IN) ::input
      integer,intent(IN) ::option
      REAL(RLEN),intent(IN) ::x1
      REAL(RLEN),intent(IN) ::y1
      end SUBROUTINE DefineSet
      end INTERFACE

      INTERFACE

      SUBROUTINE filly(nr,yinput,y)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(INOUT) ::nr
      REAL(RLEN),intent(INOUT) ::yinput
      REAL(RLEN),intent(INOUT) ::y(nr)
      end       SUBROUTINE filly
      end INTERFACE

      INTERFACE

      SUBROUTINE exchange_integer(a,b)
      implicit none
      integer,intent(INOUT) ::a
      integer,intent(INOUT) ::b
      end       SUBROUTINE exchange_integer
      end INTERFACE

      INTERFACE

      SUBROUTINE exchange_real(a,b)
      USE global_mem, ONLY:RLEN
      implicit none
      REAL(RLEN),intent(INOUT) ::a
      REAL(RLEN),intent(INOUT) ::b
      end       SUBROUTINE exchange_real
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION find_zero_c(k,NUTR,modx,option,bar,xnn,nn,a, &
      X1,X2,acc)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::nn
      integer,intent(IN) ::k
      REAL(RLEN),intent(IN) ::a
      REAL(RLEN),intent(IN) ::x1
      REAL(RLEN),intent(IN) ::x2
      REAL(RLEN),intent(IN) ::acc
      integer,intent(IN) ::nutr(nn)
      integer,intent(IN) ::modx(nn)
      integer,intent(IN) ::option(nn)
      REAL(RLEN),intent(IN) ::bar(nn)
      REAL(RLEN),intent(IN) ::xnn(2,nn)
      end FUNCTION find_zero_c
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION find_zero(NUTR1,NUTR2,mode,k,a,b,X1,X2,acc)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::nutr1
      integer,intent(IN) ::nutr2
      integer,intent(IN) ::mode
      integer,intent(IN) ::k
      REAL(RLEN),intent(IN) ::a
      REAL(RLEN),intent(IN) ::b
      REAL(RLEN),intent(IN) ::x1
      REAL(RLEN),intent(IN) ::x2
      REAL(RLEN),intent(IN) ::acc
      end FUNCTION find_zero
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION funcalc(mode,termx,coeff,xx,x)
      USE global_mem, ONLY:RLEN
      USE bennut_type
      implicit none
      integer,intent(IN) ::mode
      integer,intent(IN) ::termx
      type (ty_coeff),intent(IN) ::coeff
      REAL(RLEN),intent(IN) ::xx
      REAL(RLEN),intent(IN) ::x
      end FUNCTION funcalc
      end INTERFACE

      INTERFACE

      SUBROUTINE input_para(mode,option,input &
      ,xinput,yinput,para,equa,status)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::mode
      integer,intent(IN) ::equa
      integer,intent(IN) ::option
      integer,intent(IN) ::input
      integer,intent(INOUT) ::status
      REAL(RLEN),intent(IN) ::xinput
      REAL(RLEN),intent(IN) ::yinput
      REAL(RLEN),intent(INOUT) ::para(equa)
      end       SUBROUTINE input_para
      end INTERFACE

      INTERFACE

      integer FUNCTION kfind(nr,coeffs,n)
      USE global_mem, ONLY:RLEN
      USE bennut_type
      implicit none
      integer,intent(IN) ::nr
      integer,intent(IN) ::n
      type (ty_coeff),intent(IN) ::coeffs(n)
      end FUNCTION kfind
      end INTERFACE

      INTERFACE

      SUBROUTINE LUBKSB(N,A,INDX,B)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::n
      REAL(RLEN),intent(IN) ::a(N,N)
      integer,intent(IN) ::indx(N)
      REAL(RLEN),intent(INOUT) ::b(N)
      end       SUBROUTINE LUBKSB
      end INTERFACE

      INTERFACE

      SUBROUTINE LUDCMP(N,A,INDX,D)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::n
      REAL(RLEN),intent(INOUT) ::a(N,N)
      integer,intent(INOUT) ::indx(N)
      REAL(RLEN),intent(OUT) ::d
      end       SUBROUTINE LUDCMP
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION noutput(NUTR,mode,option,input,xinput,yinput)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::mode
      integer,intent(IN) ::option
      integer,intent(IN) ::input
      integer,intent(IN) ::nutr
      REAL(RLEN),intent(IN) ::xinput
      REAL(RLEN),intent(IN) ::yinput
      end FUNCTION noutput
      end INTERFACE
 
      INTERFACE

      REAL(RLEN) FUNCTION QGAUS_EXP(A,B,labda,FN)
      USE global_mem, ONLY:RLEN
      implicit none
      REAL(RLEN),intent(IN) ::a
      REAL(RLEN),intent(IN) ::b
      REAL(RLEN),intent(IN) ::labda(2)
      INTERFACE
        REAL(RLEN) FUNCTION FN(Y)
          USE global_mem, ONLY:RLEN
          REAL(RLEN),INTENT(IN)   ::Y
        END FUNCTION
      end INTERFACE
      end FUNCTION QGAUS_EXP
      end INTERFACE

      INTERFACE

      integer FUNCTION read_coeff(file)
      implicit none
      character,intent(IN) ::file*(*)
      end FUNCTION read_coeff
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION rtsafe(k,NUTR,mode,option,b,XNN,NN,a,X1,X2,acc)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::nn
      integer,intent(IN) ::nutr(nn)
      integer,intent(IN) ::mode(nn)
      integer,intent(IN) ::option(nn)
      integer,intent(IN) ::k
      REAL(RLEN),intent(IN) ::b(NN)
      REAL(RLEN),intent(IN) ::xnn(2,NN)
      REAL(RLEN),intent(IN) ::a
      REAL(RLEN),intent(IN) ::x1
      REAL(RLEN),intent(IN) ::x2
      REAL(RLEN),intent(IN) ::acc
      end FUNCTION rtsafe
      end INTERFACE

      INTERFACE

      SUBROUTINE rtsafe_func(mode,NUTR,input,option,b,XNN,NN,a,point,F,DF)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(INOUT) ::mode
      integer,intent(INOUT) ::nn
      integer,intent(INOUT) ::nutr(NN)
      integer,intent(INOUT) ::input(NN)
      integer,intent(INOUT) ::option(NN)
      REAL(RLEN),intent(INOUT) ::a
      REAL(RLEN),intent(INOUT) ::b(nn)
      REAL(RLEN),intent(INOUT) ::xnn(2,NN)
      REAL(RLEN),intent(INOUT) ::point
      REAL(RLEN),intent(INOUT) ::f
      REAL(RLEN),intent(INOUT) ::df
      end       SUBROUTINE rtsafe_func
      end INTERFACE

      INTERFACE

      integer FUNCTION select_coeff(mode,box)
      implicit none
      integer,intent(IN) ::mode
      integer,intent(IN) ::box
      end FUNCTION select_coeff
      end INTERFACE

      INTERFACE

      SUBROUTINE set_max_sing(nutrnr,Y,N)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(INOUT) ::nutrnr
      integer,intent(INOUT) ::n
      REAL(RLEN),intent(INOUT) ::y(n)
      end       SUBROUTINE set_max_sing
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION GetInfoFromSet(NUTR,option,input,term,x1,y1)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::nutr
      integer,intent(IN) ::term
      integer,intent(IN) ::option
      integer,intent(IN) ::input
      REAL(RLEN),intent(IN) ::x1
      REAL(RLEN),intent(IN) ::y1
      end FUNCTION GetInfoFromSet
      end INTERFACE

      INTERFACE

      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(INOUT) ::m
      integer,intent(INOUT) ::mp
      integer,intent(INOUT) ::n
      integer,intent(INOUT) ::np
      REAL(RLEN),intent(INOUT) ::b(mp)
      REAL(RLEN),intent(INOUT) ::u(mp,np)
      REAL(RLEN),intent(INOUT) ::v(np,np)
      REAL(RLEN),intent(INOUT) ::w(np)
      REAL(RLEN),intent(INOUT) ::x(np)
      end       SUBROUTINE svbksb
      end INTERFACE

      INTERFACE

      integer FUNCTION store_coeff(file,box)
      implicit none
      character file*(*)
      integer,intent(IN) ::box
      end FUNCTION store_coeff
      end INTERFACE

      INTERFACE

      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(INOUT) ::m
      integer,intent(INOUT) ::n
      integer,intent(INOUT) ::mp
      integer,intent(INOUT) ::np
      REAL(RLEN),intent(INOUT) ::a(MP,NP)
      REAL(RLEN),intent(INOUT) ::w(NP)
      REAL(RLEN),intent(INOUT) ::v(NP,NP)
      end       SUBROUTINE SVDCMP
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION transfer(mode,coeff,input,diff)
      USE global_mem, ONLY:RLEN
      USE bennut_type
      implicit none
      integer,intent(IN) ::mode
      type (ty_coeff),intent(IN) ::coeff
      REAL(RLEN),intent(IN) ::diff
      REAL(RLEN),intent(IN) ::input
      end FUNCTION transfer
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION zbrent(k,NUTR,mode,option,bar,xnn,nn,a1,X1,X2,acc)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::k
      integer,intent(IN) ::nn
      integer,intent(IN) ::nutr(nn)
      integer,intent(IN) ::mode(nn)
      integer,intent(IN) ::option(nn)
      REAL(RLEN),intent(IN) ::a1
      REAL(RLEN),intent(IN) ::bar(nn)
      REAL(RLEN),intent(IN) ::xnn(nn)
      REAL(RLEN),intent(IN) ::x1
      REAL(RLEN),intent(IN) ::x2
      REAL(RLEN),intent(IN) ::acc
      end FUNCTION zbrent
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION CalculateFromSet(NUTR,mode,input,xh,xh1)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::NUTR
      integer,intent(IN) ::mode
      integer,intent(IN) ::input
      REAL(RLEN),intent(IN) ::xh
      REAL(RLEN),intent(IN) ::xh1
      end FUNCTION CalculateFromSet
      end INTERFACE

      INTERFACE

      REAL(RLEN) FUNCTION CalculateShift(NUTR,input,xh,xh1)
      USE global_mem, ONLY:RLEN
      implicit none
      integer,intent(IN) ::NUTR
      integer,intent(IN) ::input
      REAL(RLEN),intent(IN) ::xh
      REAL(RLEN),intent(IN) ::xh1
      end FUNCTION CalculateShift
      end INTERFACE

      INTERFACE
        INTEGER FUNCTION CopySet(NUTR1,NUTR2)
        IMPLICIT  NONE
        integer,intent(IN) ::NUTR1 ! Specification
        integer,intent(IN) ::NUTR2 ! Specification
        end FUNCTION CopySet
      end INTERFACE

      INTERFACE

      SUBROUTINE AddEquation(nn_boundaries,option,input,sncf,s, xinput)
        USE global_mem, ONLY:RLEN
        USE bennut_type
        IMPLICIT  NONE

        integer, intent(IN) :: nn_boundaries
        integer, intent(IN) :: option
        integer, intent(IN) :: input
        type (ty_set),intent(IN) :: sncf
        real(RLEN), intent(IN) :: s
        real(RLEN), intent(IN) :: xinput
        end SUBROUTINE 
      end INTERFACE
  
      INTERFACE
       FUNCTION imaxloc(n,arr)
         USE global_mem, ONLY:RLEN
        integer, intent(IN) :: n
         REAL(RLEN), DIMENSION(n), INTENT(IN) :: arr
         INTEGER :: imaxloc
       END FUNCTION imaxloc
      end INTERFACE
 
      INTERFACE
        FUNCTION outerprod(n,a,b)
        USE global_mem, ONLY:RLEN
        integer, intent(IN) :: n
        REAL(RLEN), DIMENSION(n), INTENT(IN) :: a,b
        REAL(RLEN), DIMENSION(n,n) :: outerprod
       END FUNCTION outerprod
      end INTERFACE

      INTERFACE
       integer FUNCTION io_coeff(mode,option,input)
        IMPLICIT  NONE
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::option ! Specification
        integer,intent(IN) ::input ! Specification
       END FUNCTION io_coeff
      end INTERFACE

      INTERFACE
        REAL(RLEN) FUNCTION CalculateSet(NUTR, mode,option,input,&
                                                          xinput,yinput)
        USE global_mem, ONLY:RLEN
        IMPLICIT  NONE
        integer,intent(IN) ::NUTR ! Specification
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::input ! Specification
        integer,intent(IN) ::option ! Specification
        REAL(RLEN),intent(IN) ::xinput ! Specification
        REAL(RLEN),intent(IN) ::yinput ! Specification
       END FUNCTION CalculateSet

      end INTERFACE

      end MODULE bennut_interface
