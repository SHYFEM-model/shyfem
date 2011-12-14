#include <cppdefs.h>
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Runge-Kutta 2nd-order time-integration
!
! !INTERFACE
   SUBROUTINE integrationRK2
!
! !DESCRIPTION
!  Runge-Kutta 2nd-order integration with time step adjustment
! !USES
   use global_mem, ONLY:RLEN
   use mem, ONLY:NO_D3_BOX_STATES, NO_D2_BOX_STATES, &
         NO_BOXES,D3SOURCE,D3STATE,D2SOURCE,D2STATE,NO_BOXES_XY, &
         D3STATETYPE,D2STATETYPE,D3SINK,D2SINK
   use standalone
   use api_bfm
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
   real(RLEN),parameter    :: eps=0.
   real(RLEN)              :: min3D,min2D
   integer                 :: i,j,ll
   integer,dimension(2,2)  :: blccc
!
! !EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   LEVEL1 'integration RK: starting delt = ',delt
#endif
   bbccc3D=D3STATE
   bbccc2D=D2STATE
   TLOOP : DO
   ! Integration step:
      bccc3D=sum(D3SOURCE-D3SINK,2)
      bccc2D=sum(D2SOURCE-D2SINK,2)
      ccc_tmp3D=D3STATE
      ccc_tmp2D=D2STATE
	
      DO j=1,NO_D3_BOX_STATES
         IF (D3STATETYPE(j).ge.0) THEN
            D3STATE(j,:) = ccc_tmp3D(j,:) + delt*sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
!	    PRINT*,'3D',j,ccc_tmp3D(j,:),delt,sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
  	   
         END IF
      END DO
      DO j=1,NO_D2_BOX_STATES
         IF (D2STATETYPE(j).ge.0) THEN
            D2STATE(j,:) = ccc_tmp2D(j,:) + delt*sum(D2SOURCE(j,:,:)-D2SINK(j,:,:),1)
!	    PRINT*,'2D' ,j,ccc_tmp2D(j,:),delt,sum(D2SOURCE(j,:,:)-D2SINK(j,:,:),1)
         END IF
      END DO
      nmin=nmin+nstep 
      ! Check for negative concentrations
      min3D=minval(D3STATE)
      min2D=minval(D2STATE)
      
!      PRINT*,min3D,min2D,eps,delt
!      PRINT*,'                         '
      IF(min3D.lt.eps.OR.min2D.lt.eps) THEN ! cut timestep
         IF (nstep.eq.1) THEN
            blccc(:,1)=minloc(D3STATE)
            blccc(:,2)=minloc(D2STATE)
            LEVEL1 'Necessary Time Step too small! Exiting...'
            LEVEL1 blccc
            LEVEL1 var_names(stPelStateS+blccc(1,1)-1)
            LEVEL1 var_names(stBenStateS+blccc(1,2)-1)
            LEVEL1 ccc_tmp3D(blccc(1,1),blccc(2,1)), &
                           bccc3D(blccc(1,1),blccc(2,1))
            LEVEL1 ccc_tmp2D(blccc(1,2),blccc(2,2)), &
                           bccc2D(blccc(1,2),blccc(2,2))
            LEVEL1 'EXIT at time: ',timesec
            D3STATE=bbccc3D
            D2STATE=bbccc2D
            STOP
         END IF
         nstep=nstep/2
         nmin=0
         D3STATE=bbccc3D
         D2STATE=bbccc2D
         dtm1=maxdelt
         delt=nstep*mindelt
         timesec=ntime*maxdelt
#ifdef DEBUG
         LEVEL2 'Time Step cut! delt= ',delt/2.,' nstep= ',nstep
#endif
         ! Recalculate Sources:
         call ResetFluxes
         call envforcing_bfm
         call EcologyDynamics
      ELSE
#ifdef DEBUG
         LEVEL2 'Internal time step= ',nmin
#endif
         timesec=timesec+delt
#ifdef DEBUG
         LEVEL2 'Internal time= ',timesec
#endif   
         ! Recalculate sources:
         call ResetFluxes
         call envforcing_bfm
         call EcologyDynamics
         DO j=1,NO_D3_BOX_STATES
            IF (D3STATETYPE(j).ge.0) THEN
               D3STATE(j,:) = ccc_tmp3D(j,:) + &
                  .5*delt*(sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)+bccc3D(j,:))
            END IF
         END DO
         DO j=1,NO_D2_BOX_STATES
            IF (D2STATETYPE(j).ge.0) THEN
               D2STATE(j,:) = ccc_tmp2D(j,:) + &
                  .5*delt*(sum(D2SOURCE(j,:,:)-D2SINK(j,:,:),1)+bccc2D(j,:))
            END IF
         END DO
         min3D=minval(D3STATE)
         min2D=minval(D2STATE)
         IF (min3D.lt.eps.OR.min2D.lt.eps) THEN ! cut timestep
            IF (nstep.eq.1) THEN
               LEVEL1 'Necessary Time Step too small! Exiting...'
               blccc(1,:)=minloc(D3STATE)
               blccc(2,:)=minloc(D2STATE)
               LEVEL1 blccc
               LEVEL1 ccc_tmp3D(blccc(1,1),blccc(1,2)), &
                           bccc3D(blccc(1,1),blccc(1,2))
               LEVEL1 ccc_tmp2D(blccc(2,1),blccc(2,2)), &
                           bccc2D(blccc(2,1),blccc(2,2))
               LEVEL1 'EXIT at time: ',timesec
               D3STATE=bbccc3D
               D2STATE=bbccc2D
               STOP 'integration-RK2'
            END IF
            nstep=nstep/2
            nmin=0
            D3STATE=bbccc3D
            D2STATE=bbccc2D
            dtm1=delt
            delt=nstep*mindelt
            timesec=ntime*maxdelt
            LEVEL1 'Time Step cut at RK2! delt= ',delt,' nstep= ',nstep
         ELSE
            IF (nmin.eq.nmaxdelt) EXIT TLOOP
         ENDIF
         ! Recalculate Sources:
         call ResetFluxes
         call envforcing_bfm
         call EcologyDynamics
      END IF
   END DO TLOOP
   nstep=nmaxdelt
   nmin=0
   ntime=ntime+1
   delt=nstep*mindelt
   timesec=delt*ntime
#ifdef DEBUG
   LEVEL1 'ntime: ',ntime
   LEVEL1 'Integration time: ',time
#endif


   END SUBROUTINE integrationRK2
!EOC
!-----------------------------------------------------------------------
