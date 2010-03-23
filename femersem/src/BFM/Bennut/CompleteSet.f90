!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   CompleteSet
!   filly
!
! FILE
!   completeset.f90
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
      SUBROUTINE CompleteSet(NUTR,mode,option,input,xinput,yinput)
        USE global_mem, ONLY:RLEN,ALLOC,error_msg_prn
        USE bennut_variables
        USE constants
        USE bennut_constants
        USE bennut_interface,ONLY:kfind, AddEquation, transfer
        IMPLICIT  NONE
        integer,intent(IN) ::NUTR ! Specification ...dummy
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::input ! Specification
        integer,intent(IN) ::option ! Specification
        REAL(RLEN),intent(IN) ::xinput ! Specification
        REAL(RLEN),intent(IN) ::yinput ! Specification

        integer ::i
        integer ::j
        integer ::k
        integer ::l
        integer ::m
        integer ::n
        REAL(RLEN) ::r
        REAL(RLEN) ::s
        REAL(RLEN) ::t
        REAL(RLEN) ::u
        logical ::control

        if (NUTR /= nutr_seq) stop 'error: wrong use of routine'

        select case (ns%status)
          case(DEFINE)
              nn_boundaries=0
              call DefineSet(NUTR,SET_BOUNDARY,0,0,0.0D+00,0.0D+00)
              C=0.0D+00
              Y2=0.0D+00
              ns%status=SET_BOUNDARY
              ModeMass=.false.
              ns%imethod=0
          case (SET_BOUNDARY,ADD) 
          case default 
             STOP 'equation only defined after finishing definition'
        end select 

        select case (mode)
!         case (SET)
!           select case (option)
!              case (FLAG)
!                if (input == STANDARD) ModeMass=.false.
!                if (input == MASS) ModeMass=.true.
!              case (METHOD)
!                 ns%imethod=input
!              case(COEFFICIENT)
!                nn_boundaries=nn_boundaries+1
!                j=kfind(input,ns%coeffs,ns%nn)
!                call calcadd(nn_boundaries,j,C,ns%nn,1.0D+00)
!                call filly(nn_boundaries,xinput,Y2)
!           end select 
          case ( ADD,SUBTRACT,SET_BOUNDARY)
            if (mode == SET_BOUNDARY) nn_boundaries=nn_boundaries+1
            s=1.0D+00 
            if (mode  == SUBTRACT) s=-1.0D+00
            if ( option > 0 ) &
              call AddEquation(nn_boundaries,option,input,ns,s,xinput)
            call filly(nn_boundaries, s* yinput ,Y2)
          case (INPUT_TERM,START_ADD_TERM,INPUT_ADD_TERM,&
                                              INPUT_SUBTRACT_TERM)  
            !special continuity equations.........
            if (mode == INPUT_TERM.or. &
              mode == START_ADD_TERM) nn_boundaries=nn_boundaries+1
            s=1.0D+00
            if (mode == INPUT_SUBTRACT_TERM ) s=-1.0D+00
            if (mode == INPUT_TERM) then
              r=1.D+00
              t=yinput
            else
              if (yinput == 0.0D+00) &
                             stop 'CompleteSet:input_term yinnput=0.0'
              r=1.D+00/yinput
              t=1.D+00
            endif
            j=kfind(option,ns%coeffs,ns%nn)
            if (j > 100000) stop 'CompleteSet=40??'
            i=ns%coeffs(j)%il/10
            if (input == PARAMETER) then
              r=transfer(COEFF2PARA,ns%coeffs(j),r,ns%diff(i))
            elseif(input /= STANDARD) then
              stop 'CompleteSet:error mode=input_term'
            endif
            call calcadd(nn_boundaries,j,C,ns%nn,s*r)
            call filly(nn_boundaries,s*t,Y2)
          case (SET_CONTINUITY)
            if (option == FLAG) then
              if (input == MASS) ModeMass=.true.
            endif
            s=-1.D+00
            !For all lsts:
            if (ns%equa > 1) then
              do m=2,ns%equa
                !For the first derivative and the equation:
                do l=EQUATION,DERIVATIVE,DERIVATIVE
                  if (ns%lst(m-1) /= ns%lst(m)) then
                    nn_boundaries=nn_boundaries+1
                    !Make sum F1(x)=F2(X)
                    do k=m-1,m
                      t=1.D+00
                      if (l == DERIVATIVE ) t=ns%diff(k)*ns%poro(k)
                      s=-s
                      call AddEquation(nn_boundaries,k,l,ns,s*t,ns%b(m))
                    enddo
                  endif
                enddo
              enddo
            endif
          case (SET_LAYER_INTEGRAL,SET_LAYER_INTEGRAL_UNTIL )
            t=1.0D+00
            s=1.0D+00
            nn_boundaries=nn_boundaries+1
            do m=option,input
              if (ModeMass) t=ns%poro(m)*(ns%ads(m)+1.D+00)
              r=ns%b(m)
              u=ns%b(m+1)
              if (mode == SET_LAYER_INTEGRAL_UNTIL) then
                if (m == input) u=xinput
              endif
              do k=m,m+1
                s=-s
                call AddEquation(nn_boundaries,m,INTEGRAL,ns,s*t,r)
                r=u
              enddo
            enddo
            call filly(nn_boundaries,yinput,Y2)
        end select
        return
      end
 
!
      SUBROUTINE filly(nr,yinput,y)
        USE global_mem, ONLY:RLEN
        IMPLICIT  NONE
        integer,intent(IN) ::nr ! Specification
        REAL(RLEN),intent(IN)    ::yinput ! Specification
        REAL(RLEN),intent(INOUT) ::y(nr) ! Specification

        y(nr)= y(nr)+yinput

        return
      end

      SUBROUTINE calcadd(irow,k,C,nn,fc)
       USE global_mem, ONLY:RLEN
       IMPLICIT  NONE
        integer,intent(IN) ::irow ! Specification
        integer,intent(IN) ::nn ! Specification
        integer,intent(IN) ::k ! Specification
        REAL(RLEN),intent(INOUT) ::c(nn,nn) ! Specification
        REAL(RLEN),intent(IN) ::fc ! Specification

        C(irow,k)=C(irow,k)+fc

        return
      end




      SUBROUTINE AddEquation(nn_bound,layer,mode,nt,s, xinput)
        USE global_mem, ONLY:RLEN
        USE bennut_interface, ONLY:funcalc
        USE bennut_variables
        USE bennut_type
        IMPLICIT  NONE

        integer, intent(IN) :: nn_bound
        integer, intent(IN) :: layer
        integer, intent(IN) :: mode
        type (ty_set), intent(IN) :: nt
        real(RLEN), intent(IN) :: s
        real(RLEN), intent(IN) :: xinput

        integer ::j
        real(RLEN) ::r

        do j=nt%lst(layer),nt%lfi(layer)
           r=s*funcalc(mode,0,nt%coeffs(j),nt%b(layer),xinput)
           call calcadd(nn_bound,j,C,nt%nn,r)
        enddo

      end   
