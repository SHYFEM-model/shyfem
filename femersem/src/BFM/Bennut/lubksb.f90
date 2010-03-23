!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!          BFM - Biogeochemical Flux Model version 2.3
!
! SUBROUTINE
!   lubksb.f90
!
! FILE
!   lubksb.f90
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



       SUBROUTINE lubksb(n,a,indx,b)
       USE global_mem, ONLY:RLEN
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n
       REAL(RLEN), DIMENSION(n,n), INTENT(IN) :: a
       INTEGER, DIMENSION(n), INTENT(IN) :: indx
       REAL(RLEN), DIMENSION(n), INTENT(INOUT) :: b
       INTEGER :: i,ii,ll
       REAL(RLEN) :: summ
       ii=0
       do i=1,n
              ll=indx(i)
              summ=b(ll)
              b(ll)=b(i)
              if (ii /= 0) then
                     summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
              else if (summ /= 0.0) then
                     ii=i
              end if
              b(i)=summ
       end do
       do i=n,1,-1
              b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
       end do
       END SUBROUTINE lubksb

