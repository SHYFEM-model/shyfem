!Piet Ruardij (rua@nioz.nl)-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   DefineSet
!
! FILE
!   DeFineset.f90
!
! DESCRIPTION
!   !
!  
!   This file is generated from f77 code, using a code generator which
!   transposes from the F77 code into F90
!   the code USES module file as used in the BFM model
!   F90 code generator written by P. Ruardij. 
!
! AUTHORS
! Piet Ruardij (rua@nioz.nl)  
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
      SUBROUTINE DefineSet(NUTR,mode,option,input,xinput,yinput)
        USE global_mem, ONLY:RLEN
        USE bennut_variables
        USE constants
        USE bennut_constants
        USE bennut_interface,ONLY:input_para, transfer

        IMPLICIT  NONE
        integer,intent(IN) ::NUTR ! Specification
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::input ! Specification
        integer,intent(IN) ::option ! Specification
        REAL(RLEN),intent(IN) ::xinput ! Specification
        REAL(RLEN),intent(IN) ::yinput ! Specification


        integer ::i
        integer ::j
        integer ::j_1

        integer,save  ::nr=0

        REAL(RLEN) ::r
        REAL(RLEN) ::s

        if (NUTR.ne.nutr_seq) stop 'error: wrong use of routine'

        select case (mode)
        case(LAYERS)
          nr=0
          if (ns%status == LAYERS) then
            call input_para(0,option,input &
            ,xinput,yinput,ns%b(2),ns%equa-1,i)
            if (i == 1) ns%status= DIFFUSION
          else
            stop 'error ninput:input_lsts'
          endif
        case(DIFFUSION)
          if (ns%status == DIFFUSION) then
            if (ns%b(ns%equa+1).lt.0.0D+00) ns%b(ns%equa+1)=1.0D30
            call input_para(1,option,input &
            ,xinput,yinput,ns%diff,ns%equa,i)
            if (i == 1) ns%status= POROSITY
          else
            stop 'error ninput:input_lsts'
          endif
        case(POROSITY)
          if (ns%status == POROSITY) then
            call input_para(1,option,input &
            ,xinput,yinput,ns%poro,ns%equa,i)
            if (i == 1) ns%status= ADSORPTION
          else
            stop 'error ninput:input_lsts'
          endif
        case(ADSORPTION)
          if (ns%status == ADSORPTION) then
            call input_para(1,option,input &
            ,xinput,yinput,ns%ads,ns%equa,i)
            if (i == 1) ns%status= DEFINE
          else
            stop 'error ninput:input_lsts'
          endif
        case(DEFINE,DOUBLE_DEFINE,PARAMETER_DEFINE)
          if (ns%status /= DEFINE.and.ns%equa /= 1) stop &
            'defineset:mode=DEFINE/DOUBLE_DEFINE:not allowed now'
          i=option/10
          j=option-10*i
          if (i > ns%equa.or.(j < 0.or.j > 9)) stop &
            'mode=DEFINE/DOUBLE_DEFINE/DOUBLE_DEINFE:not in interval'
          nr=nr+1
          ns%coeffs(nr)%ia=input
          ns%coeffs(nr)%il=option
          ns%coeffs(nr)%labda=0.0D+00
          if ( input < 0 ) ns%coeffs(nr)%labda(1)=xinput
          if (mode /= DEFINE) ns%coeffs(nr)%labda(2)=yinput
          if (mode == PARAMETER_DEFINE) then
            r = transfer( LABDA_1,ns%coeffs(nr), xinput, ns%diff(i))
            s = transfer( LABDA_2,ns%coeffs(nr), yinput, ns%diff(i))
            ns%coeffs(nr)%labda(1) = r;
            ns%coeffs(nr)%labda(2) = s;
          endif
        case (SET_BOUNDARY)
          if (nr == ns%nn.and.(ns%status == DEFINE.or.ns%equa == 1)) then
            !wrap away all terms not necessary
            !the set of DE will consists nr coefficients (=terms)
            !  and so and we have to define at most nn equations to solve the &
            ! system.
            !lst()=0 : new definition:
            if (ns%lst(1) == 0) then
              ns%nn=nr
              do i=1,nr
                if ( any(ns%coeffs(i+1:nr)%il < ns%coeffs(i)%il)) &
                    stop 'defineset: terms notin increasing sequence'
                j=ns%coeffs(i)%il/10
                if (ns%lst(j) == 0) then
                  ns%lst(j)=i
                  j_1=j-1
                  if (j_1.gt.0) then
                    ns%lfi(j_1)=i-1
                    if (ns%lst(j_1) == 0) &
                      stop 'defineset: terms notin increasing sequence'
                  endif
                endif
              enddo
              ns%lfi(ns%equa)=ns%nn
            endif
            nn_boundaries=0
          endif
        end select
        return
      end

