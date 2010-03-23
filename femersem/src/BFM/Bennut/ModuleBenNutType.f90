!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! MODULE
!   bennut_type
! FILE
!   ModuleBenNutType.f90
!
! DESCRIPTION
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

        MODULE bennut_type
        USE global_mem, ONLY:RLEN
        USE constants, ONLY: NCOEFF,NLAYER 


        IMPLICIT NONE

        type ty_coeff
             integer :: il
             integer :: ia
             real(RLEN) :: labda(2)
        end type ty_coeff


        type ty_set 
            integer :: imethod
            integer :: equa, status, nn
            type (ty_coeff)::  coeffs(NCOEFF)
            real(RLEN)::  factor(NCOEFF)
            real(RLEN)::  diff(NLAYER)
            real(RLEN)::  poro(NLAYER)
            real(RLEN)::  ads(NLAYER)
            real(RLEN):: b(NLAYER)
            integer:: lfi(NLAYER)       
            integer:: lst(NLAYER)       
            real(RLEN):: other(NLAYER)
         end type ty_set

        end MODULE bennut_type
