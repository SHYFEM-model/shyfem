!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! SUBROUTINE
!   set_max_sing.f90
!
! FILE
!   set_max_sing.f90
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

      SUBROUTINE set_max_sing(NUTRNR,Y,N)
        USE global_mem, ONLY:RLEN
        USE mem, ONLY:BoxNumberX,BoxNumberY
        IMPLICIT  NONE
        integer,intent(INOUT) ::NUTRNR ! Specification
        integer,intent(INOUT) ::n ! Specification
        REAL(RLEN),intent(INOUT) ::y(n) ! Specification

        integer ::i
        integer ::j
        REAL(RLEN) ::ymax
        REAL(RLEN) ::ymin

        j=0
        ymax=0.0D+00
        do i=1,n
          if (y(i).gt.ymax) ymax=y(i)
        enddo
        ymin=1.D-12* ymax
        do i=1,n
          if (y(i).lt.ymin) then
            y(i)=0.0D+00
            if (j.eq.0) j=i
          endif
        enddo
        if (nutrnr.gt.0.and.j.gt.0) then
          write(0,'(''Calculateset 1:ill conditioned set BoxNumberX='' &
                               ,I4,''BoxNumberY='',I4,''SetNr='',I4)') &
            BoxNumberX, BoxNUmberY,NUTRNR
        endif

        return
      end

