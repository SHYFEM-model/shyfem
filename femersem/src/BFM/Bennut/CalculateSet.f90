!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   calculateset.f90
!
! FILE
!   calculateset.f90
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
      REAL(RLEN) FUNCTION CalculateSet(NUTR, mode,option,input,&
                                                          xinput,yinput)
        USE global_mem, ONLY:RLEN,ALLOC,error_msg_prn
        USE bennut_variables, ONLY:ns,Y2,C,nutr_seq,nn_boundaries
        USE constants
        USE bennut_constants
        USE bennut_interface,ONLY:CompleteSet, CalculateFromCondition
        IMPLICIT  NONE
        integer,intent(IN) ::NUTR ! Specification
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::input ! Specification
        integer,intent(IN) ::option ! Specification
        REAL(RLEN),intent(IN) ::xinput ! Specification
        REAL(RLEN),intent(IN) ::yinput ! Specification

        real(RLEN),dimension(NCOEFF)  :: Y3
        real(RLEN),dimension(NCOEFF*NCOEFF):: C2
        real(RLEN),dimension(NCOEFF*NCOEFF):: C3
        integer,dimension(NCOEFF) ::iindx
   
        REAL(RLEN) ::xh

        if (NUTR /= nutr_seq) then
           stop 'error: wrong use of routine'
        endif

        !calculate coefficients.......
        if (ns%nn == 0.or.nn_boundaries == 0) then
          write(0,'(''ns%nn='',i4,'' equations='',i4)') &
                                         ns%nn,nn_boundaries
          stop 'ns%nn=0 or nn_boundaries equation=0'
        elseif (ns%nn.ne.nn_boundaries) then
          write(0,'(''termsn='',i4,'' equations='',i4)')    &
                                                 ns%nn,nn_boundaries
        endif

        if ( mode == ADD .or. mode == SUBTRACT ) then
          call CompleteSet(NUTR,ADD,0,0,0.0D+00,yinput)
          ns%factor(1:nn_boundaries)=Y2(1:nn_boundaries)
          if (ns%imethod == 0) then
            call ludcmp(nn_boundaries,C,iindx,xh)
            call lubksb(nn_boundaries,C,iindx,ns%factor)
          else
            call svdcmp(C,nn_boundaries,nn_boundaries,ns%nn,   &
                                                 ns%nn,Y3,C3)
            call set_max_sing(NUTR,Y3,ns%nn)
            call svbksb(C,Y3,C3,nn_boundaries,nn_boundaries,   &
                                          ns%nn,ns%nn,Y2,ns%factor)
          endif
          ns%status=READY
          CalculateSet= yinput
        else 
          call re_Store(nn_boundaries,C,C2,ns%nn,ns%nn)
          ns%factor(1:nn_boundaries)=Y2(1:nn_boundaries)
          if (ns%imethod == 0) then
            call ludcmp(nn_boundaries,C2,iindx,xh)
            call lubksb(nn_boundaries,C2,iindx,ns%factor)
          else
            call svdcmp(C2,nn_boundaries,nn_boundaries,ns%nn,ns%nn,Y3,C3)
            call set_max_sing(NUTR,Y3,ns%nn)
            call svbksb(C2,Y3,C3,nn_boundaries,nn_boundaries,   &
                                           ns%nn,ns%nn,Y2,ns%factor)
          endif
          if ( mode ==0 ) then
            ns%status=READY
            CalculateSet=0.0D+00
          else
            !Set n'th row on zero:
            Y2(nn_boundaries:nn_boundaries)=0.0D+00
            nn_boundaries=ns%nn-1
            call re_Store(-nn_boundaries,C,C,ns%nn,ns%nn)
            call CompleteSet(NUTR,mode,option,input,xinput,yinput)
            xh=CalculateFromCondition(ns%nn,C,ns%factor,ns%nn)
            xh=max(yinput,xh)
            CalculateSet=xh
            ns%status=ADD
          endif
          return
        endif
      end

      REAL(RLEN) FUNCTION CalculateFromCondition (irow,C,Y,nn)
        USE constants,ONLY: RLEN
        IMPLICIT  NONE
        integer,intent(IN) ::irow ! Specification
        integer,intent(IN) ::nn ! Specification
        REAL(RLEN),intent(IN) ::c(nn,nn) ! Specification
        REAL(RLEN),intent(IN) ::y(nn) ! Specification
 
        CalculateFromCondition=dot_product(C(irow,1:nn),Y(1:nn))
        return
      end
