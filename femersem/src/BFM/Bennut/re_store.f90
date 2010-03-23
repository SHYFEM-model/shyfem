!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! SUBROUTINE
!   re_store.f90
!
! FILE
!   re_store.f90
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
      SUBROUTINE re_store(nr,C,C2,column,row)
        USE global_mem, ONLY:RLEN
        IMPLICIT  NONE
        integer,intent(INOUT) ::nr ! Specification
        integer,intent(INOUT) ::column ! Specification
        integer,intent(INOUT) ::row ! Specification
        REAL(RLEN),intent(INOUT) ::c(column,row) ! Specification
        REAL(RLEN),intent(INOUT) ::c2(column,row) ! Specification

        if (nr.lt.0) then
          C2(-nr+1:column,1:row)=0.0D+00
        else
          if (nr.gt.0) C2(1:nr,1:row)=C(1:nr,1:row)
          if (nr.lt.column) C2(nr+1:column,1:row)=0.0D+00;
        endif
      end

