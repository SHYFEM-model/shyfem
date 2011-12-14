!-----------------------------------------------------------------------
!BOP
!
! !PROGRAM: main
! 
! !INTERFACE:
   PROGRAM main
!
! !DESCRIPTION: 
!
! !USES:
   use standalone
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO) and Marcello Vichi (INGV)
!
! !LOCAL VARIABLES:
! 
!EOP
!-----------------------------------------------------------------------
!BOC

      call init_standalone
      call timestepping

      END PROGRAM main
!EOC
