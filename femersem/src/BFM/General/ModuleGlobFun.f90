!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ModuleGlobFun
!
! DESCRIPTION
!   List of general model functions

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  MODULE mem_globalfun
!
! !USES:
  USE global_mem, ONLY:RLEN, ZERO, BASETEMP

!  
!
! !AUTHORS
!   mfstep/ERSEM team
!
! !REVISION_HISTORY
!   --------
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
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
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  ! SHARED GLOBAL FUNCTIONS (must be below contains)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  contains
FUNCTION INSW_VECTOR(input)
        real(RLEN),intent(IN) ::input(:)
        real(RLEN) ::INSW_VECTOR(size(input))

        INSW_VECTOR =0.0D+00
        where (input > 0.0D+00 ) INSW_VECTOR=1.0D+00 

        end function INSW_VECTOR
FUNCTION MM_VECTOR(vector,param)
        real(RLEN),intent(IN) ::param
        real(RLEN)            ::vector(:)
        real(RLEN)            ::MM_VECTOR(size(vector))

        MM_VECTOR= VECTOR / ( VECTOR+  PARAM)

        end function MM_VECTOR
FUNCTION eramp_VECTOR(x, m)
        real(RLEN),intent(IN) ::x(:)
        real(RLEN),intent(IN) ::m
        real(RLEN)            ::eramp_VECTOR(size(x))

        eramp_VECTOR =0.0D+00
        where (x > 0 ) 
           eramp_VECTOR=1.0D+00 
           where (X< M) eramp_VECTOR=X/M;
        endwhere

        end function
FUNCTION PartQ_vector(p, d_a, d_b, d_m)
        real(RLEN),intent(IN) ::p(:)
        real(RLEN),intent(IN) ::d_a(:)
        real(RLEN),intent(IN) ::d_b(:)
        real(RLEN),intent(IN) ::d_m
        real(RLEN)            ::PartQ_VECTOR(size(p))

        real(RLEN),dimension(:),allocatable ::c1
        REAL(RLEN),dimension(:),allocatable ::b1 
        REAL(RLEN),dimension(:),allocatable ::a1 
        REAL(RLEN),dimension(:),allocatable ::norm 
        REAL(RLEN),dimension(:),allocatable ::r 

        ALLOCATE(c1(size(p)))
        ALLOCATE(b1(size(p)))
        ALLOCATE(a1(size(p)))
        ALLOCATE(norm(size(p)))
        ALLOCATE(r(size(p)))

        c1 = min(p * (-1.0 * log(1.0D-20)), d_m);
        b1 = min(d_b, c1);
        a1 = min(d_a, b1);
        r=0.0D+00
        where ( d_a ==0.0 ) r=1.0D+00

        where (c1 > 0.0D+0 .and. p /= 0.0D+00)
           norm = 1.0D+00 - exp((- c1) / p);
           PartQ_vector= (exp( -a1 / p) - exp(- b1 / p)) / norm;
        elsewhere 
           PartQ_VECTOR =r
        endwhere

        end function
function eTq_VECTOR(temp,p_q10)

        IMPLICIT NONE
        real(RLEN)            :: temp(:)
        real(RLEN),intent(IN) :: p_q10
        real(RLEN)            ::eTq_VECTOR(size(temp))

        eTq_VECTOR=  exp(  log(  p_q10)*( temp- BASETEMP)/ BASETEMP)
        end function eTq_VECTOR
function IntegralExp(alfa,x)
        IMPLICIT NONE
        real(RLEN),intent(IN)       :: alfa
        real(RLEN),intent(IN)       :: x
        real(RLEN)                 :: IntegralExp
        
        IntegralExp=(exp(alfa * x) -1.0D+00)/alfa
        end function IntegralExp
FUNCTION INSW(input)
        real(RLEN),intent(IN) ::input
        real(RLEN) ::INSW

        INSW =0.0D+00
        if (input > 0.0D+00 ) INSW=1.0D+00 

        end function INSW
  end module
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
