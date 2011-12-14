#include <cppdefs.h>
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Euler-forward time-integration
!
! !INTERFACE
   SUBROUTINE integrationEfw
!
! !DESCRIPTION
!  Euler-forward integration with time step adjustment
! !USES
   use global_mem, ONLY:RLEN
   use mem, ONLY:NO_D3_BOX_STATES, NO_D2_BOX_STATES, &
         NO_BOXES,D3SOURCE,D3STATE,D2SOURCE,D2STATE,NO_BOXES_XY, &
         D3STATETYPE,D2STATETYPE,D3SINK,D2SINK
   use standalone
   use api_bfm
   use time, only: update_time
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
   real(RLEN),parameter      :: eps=0.
   real(RLEN)                :: min3D,min2D
   integer                   :: i,j,ll
   integer,dimension(2,2)    :: blccc
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   LEVEL1 'integration efw: starting delt = ',delt
#endif
   bbccc3D=D3STATE
   bbccc2D=D2STATE
   TLOOP : DO
   ! Integration step:
      DO j=1,NO_D3_BOX_STATES
         IF (D3STATETYPE(j).ge.0) THEN
            D3STATE(j,:) = D3STATE(j,:) + delt*sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
         END IF
      END DO
      if (bio_setup>=2) then
         DO j=1,NO_D2_BOX_STATES
            IF (D2STATETYPE(j).ge.0) THEN
               D2STATE(j,:) = D2STATE(j,:) + delt*sum(D2SOURCE(j,:,:)-D2SINK(j,:,:),1)
            END IF
         END DO
      end if
      nmin=nmin+nstep 
   !  Check for negative concentrations
      min3D=minval(D3STATE)
      min2D=minval(D2STATE)
      IF(min3D.lt.eps.OR.min2D.lt.eps) THEN ! cut timestep
         IF (nstep.eq.1) THEN
            LEVEL1 'Necessary Time Step too small! Exiting...'
            blccc(:,1)=minloc(D3STATE)
            blccc(:,2)=minloc(D2STATE)
            LEVEL1 blccc
            LEVEL1 'Pelagic Variable:',trim(var_names(stPelStateS+blccc(1,1)-1))
            LEVEL1 'Value: ',D3STATE(blccc(1,1),blccc(2,1)),' Rate: ', &
                        bbccc3D(blccc(1,1),blccc(2,1))
            LEVEL1 'Benthic Variable:',trim(var_names(stBenStateS+blccc(1,2)-1))
            LEVEL1 'Value: ',D2STATE(blccc(1,2),blccc(2,2)),' Rate: ', &
                        bbccc2D(blccc(1,2),blccc(2,2))
            LEVEL1 'EXIT at  time ',timesec
            STOP 'integration-efw'
         END IF
         nstep=nstep/2
         nmin=0
         D3STATE=bbccc3D
         D2STATE=bbccc2D
         dtm1=delt
         delt=nstep*mindelt
         timesec=ntime*maxdelt
         LEVEL1 'Time Step cut! delt= ',delt,' nstep= ',nstep
      ELSE
#ifdef DEBUG
      LEVEL2 'Internal time step= ',nmin
#endif
         IF(nmin.eq.nmaxdelt) EXIT TLOOP
         timesec=timesec+delt
#ifdef DEBUG
      LEVEL2 'Internal time= ',timesec
#endif
      END IF
      ! Recalculate Sources:
!     call ResetFluxes
!     call envforcing_bfm
      call EcologyDynamics
   END DO TLOOP
   nstep=nmaxdelt
   nmin=0
   ntime=ntime+1
   call update_time(ntime)
   dtm1=delt
   delt=nstep*mindelt
   timesec=delt*ntime
#ifdef DEBUG
   LEVEL1 'ntime: ',ntime
   LEVEL1 'Integration time: ',timesec
#endif

   END SUBROUTINE integrationEfw
!EOC
!-----------------------------------------------------------------------
