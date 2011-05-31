!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   funcalc.f90
!
! FILE
!   funcalc.f90
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
      REAL(RLEN) FUNCTION funcalc(mode,chterm,coeff,basis,x)
        USE global_mem, ONLY:RLEN
        USE bennut_type
        USE constants
        USE bennut_interface,ONLY:BESSK1, BESSK0, BESSI1, BESSI0, &
          QGAUS_EXP
        IMPLICIT  NONE
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::chterm ! Specification
        type (ty_coeff),intent(IN) ::coeff ! Specification
        REAL(RLEN),intent(IN) ::basis ! Specification
        REAL(RLEN),intent(IN) ::x ! Specification

        REAL(RLEN),parameter :: pi=3.141592D+00
        REAL(RLEN) ::r
        REAL(RLEN) ::rx
        REAL(RLEN) ::s
        integer ::i
        integer ::j
        integer ::term

        rx=x-basis
        term=coeff%ia
        if ( chterm .ne. 0 )term=max(term-2,min(term,0))
        select case (term)
        case (CONSTANT_TERM)
          select case (mode)
            case (INTEGRAL)
               funcalc=rx
            case (EQUATION,DERIVATIVE,SDERIVATIVE)
               funcalc=DBLE(max(1+mode,0))
            case default
               stop 'funcalc CONSTANT_TERM'
          end select
        case (LINEAR_TERM:)
          !integration of terms......
          select case (mode)
            case (INTEGRAL)
              funcalc=(rx**(term+1))/(term+1)
            case (EQUATION)
               funcalc=rx**(term)
            case (DERIVATIVE,SDERIVATIVE)
              !calculation for differention
              j=term+mode
              if (j.lt.0) then
                funcalc=0.0D+00
              else
                r=term
                if (mode == SDERIVATIVE .and. term > 2) r=(term-1)*r
                if (j > 0) r=r*(rx**(j))
                funcalc=r
              endif
            case default
               stop 'error funcalc FIRST_ORDER_TERM'
          end select
        case ( EXPONENTIAL_TERM,ZERO_EXPONENTIAL_TERM)
          r=exp(coeff%labda(1)*rx)

          select case (mode)
            case (INTEGRAL)
              funcalc=r/coeff%labda(1)
            case (EQUATION)
              funcalc=r
            case (DERIVATIVE,SDERIVATIVE)
               funcalc=(coeff%labda(1)**(-mode))*r
            case (EXPONENTIAL_INTEGRAL)
              funcalc=0.0D+00
            case default
              stop 'funcalc EXPONENTIAL_TERM'
          end select 
        case ( BESSELI_EXP_TERM)
          r=exp(coeff%labda(1)*rx)
          s=r*coeff%labda(2)

          select case (mode)
            case (EXPONENTIAL_INTEGRAL)
              !equation: 11.3.25 :
              if (r == 0.0D+00) then
                funcalc=0.0D+00
              else
                funcalc=r*bessi1(s)/(coeff%labda(1)*coeff%labda(2))
              endif
            case (INTEGRAL)
              r=0.0D+00
              if (rx.ne.0.0D+00)r= QGAUS_EXP(0.0D+00,rx,coeff%labda,bessi0)
              funcalc=r
            case (EQUATION)
                funcalc=bessi0(s) 
            case (DERIVATIVE)
                funcalc=bessi1(s)*coeff%labda(1)*s
            case default
              stop 'funcalc BESSELI_EXP_TERM'
            end select 
        case ( BESSELK_EXP_TERM)
          r=exp(coeff%labda(1)*rx)
          s=r*coeff%labda(2)

          select case (mode)
            case (EXPONENTIAL_INTEGRAL)
              if (r == 0.0D+00) then
                funcalc=0.0D+00
              else
                funcalc=r*bessk1(s)/(coeff%labda(1)*coeff%labda(2))
              endif
            case (INTEGRAL)
              r=0.0D+00
              if (rx.ne.0.0D+00)r= qgaus_exp(0.0D+00,rx,coeff%labda,bessk0)
              funcalc=r
            case (EQUATION)
              funcalc=bessk0(s) 
            case (DERIVATIVE)
              funcalc=-bessk1(s)*coeff%labda(1)*s
            case default
              stop 'funcalc BESSELK_EXP_TERM'
          end select 
        case default
           stop 'funcalc no such term'
        end select 
        return
      end
