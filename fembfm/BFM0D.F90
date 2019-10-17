#include "cppdefs.h"
!
! the next directive specifies if we have some new features available
! uncomment if you use a version higher than 5.1 or 5.1 develop
!
!#define version_ggu_5_1_dev
!
subroutine BFM0D_INIT_IO_CHANNELS()
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifdef version_ggu_5_1_dev
  use global_mem, only:RLEN,LOGUNIT,NMLUNIT,bfm_file_FirstUnit
  use api_bfm, ONLY: GetLun
#else
  use global_mem, only:RLEN,LOGUNIT,NMLUNIT
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
#ifdef version_ggu_5_1_dev
  LOGUNIT=GetLun()
  NMLUNIT=GetLun()+1
#else
  LOGUNIT=777
#endif

end subroutine BFM0D_INIT_IO_CHANNELS
!
subroutine BFM0D_Input_EcologyDynamics(sur,bot,BFM0D_trn,dim_BFM0D_trn,BFM0D_er)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN
  use api_bfm, ONLY: SRFindices, BOTindices
  use mem
  use mem_CO2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  logical   , intent(in) :: sur, bot
  integer dim_BFM0D_trn
  real(RLEN), intent(in) :: BFM0D_trn(dim_BFM0D_trn), BFM0D_er(10)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer i

    IF(sur) then
       SRFindices(1) = 1
    ELSE
       SRFindices(1) = 0
    ENDIF


    IF(bot) then
       BOTindices(1) = 1
    ELSE
       BOTindices(1) = 0
    ENDIF


  DO i=1,dim_BFM0D_trn
      D3STATE(i,1)=BFM0D_trn(i)
!       LEVEL1 'BFM0D_Input_EcologyDynamics:D3STATE(',i,')=',D3STATE(i,1) 
  END DO

  ! Environmental regulating factors
  ETW    = BFM0D_er(1)
! LEVEL1 'BFM0D_Input_EcologyDynamics:ETW', ETW
  ESW    = BFM0D_er(2)
! LEVEL1 'BFM0D_Input_EcologyDynamics:ESW', ESW
  ERHO   = BFM0D_er(3)
! LEVEL1 'BFM0D_Input_EcologyDynamics:ERHO', ERHO
  EICE   = BFM0D_er(4)
! LEVEL1 'BFM0D_Input_EcologyDynamics:EICE', EICE
#ifdef INCLUDE_PELCO2
  AtmCO2%fnow = BFM0D_er(5)
! LEVEL1 'BFM0D_Input_EcologyDynamics:EPCO2air', EPCO2air
#endif
  EIR    = BFM0D_er(6)
! LEVEL1 'BFM0D_Input_EcologyDynamics:EIR', EIR
  SUNQ   = BFM0D_er(7)
! LEVEL1 'BFM0D_Input_EcologyDynamics:DL', SUNQ
  DEPTH  = BFM0D_er(8)
! LEVEL1 'BFM0D_Input_EcologyDynamics:DEPTH', DEPTH
  EWIND  = BFM0D_er(9)
! LEVEL1 'BFM0D_Input_EcologyDynamics:EWIND,', EWIND
  ph     = BFM0D_er(10)
! LEVEL1 'BFM0D_Input_EcologyDynamics:PH,', ph

end subroutine BFM0D_Input_EcologyDynamics
!
subroutine BFM0D_NO_BOXES(N,X,Y,Z,XY)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN
  use mem, ONLY:NO_BOXES,NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z,NO_BOXES_XY, NO_STATES, &
       NO_D3_BOX_STATES
  use api_bfm, ONLY: SRFindices,BOTindices

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, intent(in) :: N,X,Y,Z,XY

  NO_BOXES   = N  
  NO_BOXES_X = X
  NO_BOXES_Y = Y
  NO_BOXES_Z = Z
  NO_BOXES_XY = XY
  NO_STATES   = NO_D3_BOX_STATES * NO_BOXES
  write(*,*) 'NO_BOXES_TOT', N
  write(*,*) 'NO_BOXES_X', X
  write(*,*) 'NO_BOXES_Y', Y
  write(*,*) 'NO_BOXES_Z', Z
  write(*,*) 'NO_BOXES_XY', XY
  write(*,*) 'allocating NO_BOXES_XY'
  allocate(SRFindices(NO_BOXES_XY))
  allocate(BOTindices(NO_BOXES_XY))
  
  write(*,*) 'NO_STATES', NO_STATES
  
  write(*,*) 'NO_D3_BOX_STATES',NO_D3_BOX_STATES 

end subroutine BFM0D_NO_BOXES
!
subroutine BFM0D_reset
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  use (import) other modules
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, only:ZERO
  use mem

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Start the allocation of pelagic state global
  ! matrix and pointers
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer i, dim_BFM0D_tra

  call ResetFluxes()

#ifdef version_ggu_5_1_dev
  SiPToutr = ZERO
  SiPToutc = ZERO
  SiPTouto = ZERO
  SiPToutn = ZERO
  SiPToutp = ZERO
  SiPTouti = ZERO
  SoO2Airo = ZERO
  SoPTinr  = ZERO
  SoPTinc  = ZERO
  SoPTino  = ZERO
  SoPTinn  = ZERO
  SoPTinp  = ZERO
  SoPTini  = ZERO
  SoRIc    = ZERO
  SoRIn    = ZERO
  SoRIp    = ZERO
  SoRIs    = ZERO
#endif

#ifdef INCLUDE_PELCO2
  jsurO3c=0
#endif

end subroutine BFM0D_reset

! the next subroutines are only available in higher than 5.1 version
! therefore we include dummies for the older versions

#ifndef version_ggu_5_1_dev

subroutine BFM0D_Output_EcologyDynamics(BFM0D_tra, BFM0D_sed, local_BFM0D_dia)

  use global_mem, ONLY:RLEN,ZERO
  use mem

  IMPLICIT NONE
  real(RLEN), intent(out) :: BFM0D_sed( iiPhytoPlankton )
  real(RLEN), intent(out) :: BFM0D_tra( NO_D3_BOX_STATES )
  INTEGER, parameter :: jptra_var = 90
  INTEGER, parameter :: jptra_flux = 23
  INTEGER, parameter :: jptra_dia = jptra_var + jptra_flux
  real(RLEN), intent(out) :: local_BFM0D_dia(jptra_dia) 

  BFM0D_tra = ZERO
  BFM0D_sed = ZERO
  local_BFM0D_dia = ZERO

end subroutine BFM0D_Output_EcologyDynamics


subroutine BFM0D_Output_EcologyDynamics_surf(BFM0D_tra, BFM0D_sed, local_BFM0D_dia,local_BFM0D_dia2D)

  use global_mem, ONLY:RLEN,ZERO
  use mem

  IMPLICIT NONE
  real(RLEN), intent(out) :: BFM0D_sed( iiPhytoPlankton )
  real(RLEN), intent(out) :: BFM0D_tra( NO_D3_BOX_STATES )
  INTEGER, parameter :: jptra_var = 90
  INTEGER, parameter :: jptra_flux = 23
  INTEGER, parameter :: jptra_dia = jptra_var + jptra_flux
  INTEGER, parameter :: jptra = 50
  INTEGER, parameter :: jptra_dia_2d =  11

  real(RLEN), intent(out) :: local_BFM0D_dia2d(jptra_dia_2d)
  real(RLEN), intent(out) :: local_BFM0D_dia(jptra_dia) 

  BFM0D_tra = ZERO
  BFM0D_sed = ZERO
  local_BFM0D_dia = ZERO
  local_BFM0D_dia2D = ZERO

end subroutine BFM0D_Output_EcologyDynamics_surf

#endif
