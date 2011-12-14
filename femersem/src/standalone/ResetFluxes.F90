#include <cppdefs.h>
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Euler-forward time-integration
!
! !INTERFACE
   subroutine ResetFluxes
!
! !DESCRIPTION
!  Ueler-forward integration with time step adjustment
! !USES
   use Mem, ONLY: NO_D2_BOX_STATES,NO_BOXES_XY,D2SOURCE,sunq, &
         NO_D3_BOX_STATES,NO_BOXES,D3SOURCE, D3SINK, D2SINK,&
         ESS,ESW,ETW,EIR,eipi,xeps,Depth,sediR6,sediPI
   implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
   integer :: i
!EOP
!-----------------------------------------------------------------------
!BOC
!
   ! Reset source term arrays 
#ifdef DEBUG
   D3SOURCE = 0.0
   D3SINK = 0.0
   D2SOURCE = 0.0
   D2SINK = 0.0
#else
   ! only the diagonal
   do i=1,NO_D3_BOX_STATES
      D3SOURCE(i,i,:) = 0.0
      D3SINK(i,i,:) = 0.0
   end do
   do i=1,NO_D2_BOX_STATES
      D2SOURCE(i,i,:) = 0.0
      D2SINK(i,i,:) = 0.0
   end do
#endif

   sunq= 0.0
   eipi = 0.0

   ESS = 0.0
   ESW = 0.0
   ETW= 0.0
   EIR = 0.0
   xeps = 0.0

   sediPI=0.
   sediR6=0.

   end subroutine ResetFluxes
!EOC
!-----------------------------------------------------------------------
