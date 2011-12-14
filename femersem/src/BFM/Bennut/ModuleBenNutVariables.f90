!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! MODULE
!   bennut_variables
! FILE
!   ModuleBenNutVariables.f90
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

        MODULE bennut_variables
        USE global_mem, ONLY:RLEN
        USE bennut_type
        USE constants, ONLY: NCOEFF,NLAYER 

        IMPLICIT NONE

 
        integer                    :: nflag=0
        integer,public,pointer,dimension(:)        :: fflag
        type (ty_set),public,pointer,dimension(:)  :: sets
        type (ty_set),pointer,public               :: ns 

   
        logical,public                  :: ModeMass
        integer,public                  :: nutr_seq, nn_boundaries
        integer,public                  :: nutr_control=0
        real(RLEN),public,dimension(NCOEFF)  :: Y2
        real(RLEN),public,dimension(NCOEFF*NCOEFF):: C

        end MODULE bennut_variables
