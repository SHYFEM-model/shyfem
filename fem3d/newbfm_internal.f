
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! routines for BFM module
!
! revision log :
!
! 22.02.2016	ggu&erp	new bfm routines created from newconz
! 22.03.2016	ggu	changed VERS_7_5_6
! 06.06.2016	ggu	initialization from file changed
! 09.06.2016 	leslie	modified (Search LAA)
! 10.06.2016	ggu	changed VERS_7_5_13
! 17.06.2016	ggu	changed VERS_7_5_15
! 28.06.2016	ggu	initialize bfmv, new routine bfm_check()
! 09.09.2016	ggu	changed VERS_7_5_17
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 01.05.2019	erp	some routines written for linkage to BFM
! 17.10.2019	ggu	these routines eliminated into own file

!==================================================================
        module mod_bfm_internal
!==================================================================

        implicit none

	!integer, save :: ibfm_state = 50	!number of state variables

!==================================================================
	contains
!==================================================================

!==================================================================
        end module mod_bfm_internal
!==================================================================

        subroutine bfm_init_internal

	implicit none

        call BFM0D_NO_BOXES(1,1,1,1,1)
        call Initialize()
        
        end subroutine

c*********************************************************************

      subroutine bfm_reactor_internal
      
      use basin
      use levels
      use mod_sinking
      use mod_bfm

      implicit none

      include 'femtime.h'

      integer :: i,l,k
      double precision, dimension(10) :: er
      double precision, dimension(4) :: wsink
      logical :: bsur,bot
      double precision, dimension(50) :: b,b_apx,mat_dbl
      double precision, dimension(114) :: d
      double precision, dimension(11) :: d2    
      double precision, dimension(nlv,nkn) :: matrix_light
      
      matrix_light(:,:) = 0.
      call light_abs(nkn,nlv,ibfm_state,bfmv,matrix_light)
      
      DO k=1,nkn
	  DO l=1,ilhkv(k)
	  b_apx = DBLE(bfmv(l,k,:))
	  !bfmv(l,k,:) =b_apx      
	  bsur = (k .eq. 1)
	  bot = .FALSE.
      
	  er(1)  = DBLE(tempv(l,k))!tn (ji,jj,jk)        ! Temperature (Celsius)
	  er(2)  = DBLE(saltv(l,k)) !(ji,jj,jk)        ! Salinity PSU
	  er(3)  = 1025 + rhov(l,k)        ! Density Kg/m3
	  er(4)  = 0                 ! from 0 to 1 adimensional
	  er(5)  = 390           ! CO2 Mixing Ratios (ppm)  390
	  er(6)  = (max(matrix_light(l,k),0.1))/0.217   ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217 (LAA)
	  er(7)  = 24    ! fotoperiod expressed in hours
	  er(8)  = DBLE(hdknv(l,k))              ! depth in meters of the given cell
	  er(9)  = DBLE(sqrt(wxv(k)**2+wxv(k)**2)) !vatm(ji,jj) * surf_mask(jk) ! wind speed (m/s)
	  er(10) = 8.         ! PH
	  
!!!#ifdef BFM_active     
	  
	  call BFM0D_Input_EcologyDynamics(bsur,bot,b_apx,ibfm,er)
	  call BFM0D_reset()

	  call EcologyDynamics()
	  
	  if (bsur) then
               call BFM0D_Output_EcologyDynamics_surf(b,wsink,d,d2)
          else
               call BFM0D_Output_EcologyDynamics(b,wsinkv,d)
          endif
!!!#endif          
          wsinkv(l,k) = REAL(wsink(1)+wsink(2)+wsink(3)+wsink(4))/4
	  
	  bfmv(l,k,:) = REAL(b_apx + b(:)*dt_act)
	  
	  ENDDO
	ENDDO

      end subroutine
      
c*********************************************************************
      
      subroutine light_abs(nkn,nlv,nst,bfmv,matrix_light)
      
	use levels, only: ilhkv
	use mod_layer_thickness

      implicit none

      integer,intent(in) :: nkn,nlv,nst
      real bfmv(nlv,nkn,nst)
      double precision,dimension(nlv,nkn),intent(inout)::matrix_light

      integer :: k,l
      real :: solrad
      double precision :: coeff,eps1,eps0
      
      eps1= 0.0088
      eps0=0.04

      do k=1,nkn
	  
	  call meteo_get_solar_radiation(k,solrad)
	  matrix_light(1,k) = solrad
      
	  do l=2,ilhkv(k)
	    coeff = -(eps1*(bfmv(l-1,k,14)+bfmv(l-1,k,19)+ 
     &	          bfmv(l-1,k,23)+bfmv(l-1,k,27))+eps0)*hdknv(l-1,k)
	    matrix_light(l,k) = matrix_light(l-1,k) * exp(coeff) 
	  enddo
	  
      
	  do l=1,ilhkv(k)
	    coeff = -0.5*(eps1*(bfmv(l,k,14)+bfmv(l,k,19)+ 
     &	          bfmv(l,k,23)+bfmv(l,k,27))+eps0)*hdknv(l,k)
	    matrix_light(l,k) = matrix_light(l,k) * exp(coeff) 
	  enddo
	  
      enddo

      end subroutine

c*********************************************************************

      SUBROUTINE write_matrix(mtw,righe,colonne,percorso)

          IMPLICIT NONE

          INTEGER :: i, j, righe,colonne
          character(LEN=*), INTENT(in) :: percorso
          double precision, INTENT(in) :: mtw(righe,colonne)
    
          OPEN(unit=115, file=percorso)
	  print *,righe,colonne
          Do i=1,righe
          WRITE(115,'(*(F14.7))') (mtw(i,j),j=1,colonne)
          END DO
          CLOSE(0)

       END SUBROUTINE write_matrix
      
c*********************************************************************

