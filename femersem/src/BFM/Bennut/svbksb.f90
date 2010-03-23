!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! SUBROUTINE
!   svbksb.f90
!
! FILE
!   svbksb.f90
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
!   Numerical Recipes Software 
!
! CHANGE_LOG
!   {!  (C) Copr. 1986-92 Numerical Recipes Software .)1:.}
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
      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
        USE global_mem, ONLY:RLEN
        integer,intent(INOUT) ::m ! Specification
        integer,intent(INOUT) ::mp ! Specification
        integer,intent(INOUT) ::n ! Specification
        integer,intent(INOUT) ::np ! Specification
        integer ::nmax
        REAL(RLEN),intent(INOUT) ::b(mp) ! Specification
        REAL(RLEN),intent(INOUT) ::u(mp,np) ! Specification
        REAL(RLEN),intent(INOUT) ::v(np,np) ! Specification
        REAL(RLEN),intent(INOUT) ::w(np) ! Specification
        REAL(RLEN),intent(INOUT) ::x(np) ! Specification
        parameter (NMAX=500)
        integer ::i
        integer ::j
        integer ::jj
        REAL(RLEN) ::s
        REAL(RLEN) ::tmp(NMAX)
        do 12 j=1,n
          s=0.0D0
          if(w(j).ne.0.0D0)then
            do 11 i=1,m
              s=s+u(i,j)*b(i)
 11         continue
            s=s/w(j)
          endif
          tmp(j)=s
 12     continue
        do 14 j=1,n
          s=0.0D0
          do 13 jj=1,n
            s=s+v(j,jj)*tmp(jj)
 13       continue
          x(j)=s
 14     continue
        return
      END
  ! (C) Copr. 1986-92 Numerical Recipes Software .)1:.

