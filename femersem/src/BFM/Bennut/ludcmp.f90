!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!          BFM - Biogeochemical Flux Model version 2.3
!
! SUBROUTINE
!   ludcmp.f90
!
! FILE
!   ludcmp.f90
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



       SUBROUTINE ludcmp(n,a,indx,d)
       USE global_mem, ONLY:RLEN
       USE bennut_interface, ONLY:outerprod,imaxloc
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n
       REAL(RLEN), DIMENSION(n,n), INTENT(INOUT) :: a
       INTEGER, DIMENSION(n), INTENT(OUT) :: indx
       REAL(RLEN), INTENT(OUT) :: d
       REAL(RLEN), DIMENSION(n) :: vv
       REAL(RLEN), PARAMETER :: TINY=1.0D-20
       INTEGER :: j,imax

       d=1.0
       vv=maxval(abs(a),dim=2)
       if (any(vv == 0.0D+00)) stop 'singular matrix in ludcmp'
       vv=1.0D+00/vv
       do j=1,n
              imax=(j-1)+imaxloc(n-j+1,vv(j:n)*abs(a(j:n,j)))
              if (j /= imax) then
                     call swap(n,a(imax,:),a(j,:))
                     d=-d
                     vv(imax)=vv(j)
              end if
              indx(j)=imax
              if (a(j,j) == 0.0D+00) a(j,j)=TINY
              a(j+1:n,j)=a(j+1:n,j)/a(j,j)
              a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(n-j,a(j+1:n,j),a(j,j+1:n))
       end do
       END SUBROUTINE ludcmp

      SUBROUTINE swap(n,a,b)
      USE global_mem, ONLY:RLEN
      integer,intent(IN) ::n ! Specification
      REAL(RLEN), DIMENSION(n), INTENT(INOUT) :: a,b
      REAL(RLEN), DIMENSION(n) :: dum
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap

       FUNCTION imaxloc(n,arr)
       USE global_mem, ONLY:RLEN
       REAL(RLEN), DIMENSION(n), INTENT(IN) :: arr
       INTEGER :: imaxloc
       INTEGER, DIMENSION(1) :: imax
       imax=maxloc(arr(:))
       imaxloc=imax(1)
       END FUNCTION imaxloc

       FUNCTION outerprod(n,a,b)
        USE global_mem, ONLY:RLEN
       REAL(RLEN), DIMENSION(n), INTENT(IN) :: a,b
       REAL(RLEN), DIMENSION(n,n) :: outerprod
       outerprod = spread(a,dim=2,ncopies=n) * &
              spread(b,dim=1,ncopies=n)
       END FUNCTION outerprod


