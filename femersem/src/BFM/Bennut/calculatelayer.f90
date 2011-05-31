!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! SUBROUTINE
!   CalculateLayer
!   check_on_bigger
!
! FILE
!   calculatelayer.f90
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
      SUBROUTINE CalculateLayer(NUTR,mode,x,layernr,xb)
        USE global_mem, ONLY:RLEN
        USE bennut_type
        USE bennut_variables, ONLY:sets
        IMPLICIT  NONE
        integer,intent(IN) ::nutr ! Specification
        integer,intent(IN) ::mode ! Specification
        REAL(RLEN),intent(IN) ::x ! Specification
        integer,intent(OUT) ::layernr ! Specification
        REAL(RLEN),intent(OUT) ::xb ! Specification

        integer ::j
        REAL(RLEN) ::dummy

        select case (mode)
          case(0)
            j=1
            !get upper border where first equation is valid..
            xb=sets(NUTR)%b(j+1)
            do while ((x-xb)/(xb+1.0D-30).gt.-1.0D-6)
              !if x > xb (=xabove) check next interval:
              j=j+1
              xb=sets(NUTR)%b(j+1)
            enddo
            layernr=j
        case(1:)
            layernr=mode
            xb=sets(NUTR)%b(mode+1)
        case default
            layernr=sets(NUTR)%equa
        end select

        return
      end

