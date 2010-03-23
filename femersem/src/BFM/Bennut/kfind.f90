!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   kfind.f90
!
! FILE
!   kfind.f90
!
! DESCRIPTION
!   {
!}
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
      integer FUNCTION kfind(nr,coeffs,n)
        use bennut_type
        IMPLICIT  NONE
        integer,intent(IN) ::nr ! Specification
        integer,intent(IN) ::n ! Specification
        type (ty_coeff),intent(IN) ::coeffs(n) ! Specification

        integer ::i
        integer ::j

        j=100000000+n
        i=1
        do while (nr.gt.coeffs(i)%il.and.i.lt.n)
          i=i+1
        enddo
        if (coeffs(i)%il.eq.nr) j=i

        kfind= j
        return
      end

