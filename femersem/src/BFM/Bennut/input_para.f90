!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! SUBROUTINE
!   input_para.f90
!
! FILE
!   input_para.f90
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
!     Piet Ruardij, rua@nioz.nl
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
      SUBROUTINE input_para(mode,option,input &
        ,xinput,yinput,para,equa,status)
        USE global_mem, ONLY:RLEN
        IMPLICIT  NONE
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::equa ! Specification
        integer,intent(IN) ::option ! Specification
        integer,intent(IN) ::input ! Specification
        integer,intent(INOUT) ::status ! Specification
        REAL(RLEN),intent(IN) ::xinput ! Specification
        REAL(RLEN),intent(IN) ::yinput ! Specification
        REAL(RLEN),intent(INOUT) ::para(equa) ! Specification

        integer ::i

        status=1
        if (option.gt.0) then
          if (option.le.equa) then
            para(option)=xinput
            if (input.eq.option+1) para(option+1)=yinput
          endif
          do i=2,equa
            if (para(i).lt.0.0D+00) status=0
          enddo
        elseif (mode.eq.1) then
          para(1:equa)=xinput
        elseif(mode.ne.0) then
          stop 'error in input_para'
        endif
        return
      end

