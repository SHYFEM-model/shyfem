#include <cppdefs.h>
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Leap Frog time-integration
!
! !INTERFACE
  subroutine integrationLf
!
! !DESCRIPTION
!  Leap Frog time-integration with time-step adjustment
!
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
   real(RLEN),parameter    :: ass=.05
   real(RLEN)              :: min3D,min2D
   integer                 :: i,j,n
   integer,dimension(2,2)  :: blccc
   real(RLEN),dimension(NO_D3_BOX_STATES,NO_BOXES)    :: bc3D
   real(RLEN),dimension(NO_D2_BOX_STATES,NO_BOXES_XY) :: bc2D
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   LEVEL1 'integration-lf: starting delt = ',delt
#endif
   ! save initial states and additional variables
   bccc3D=D3STATE
   bccc2D=D2STATE
   bc3D=bbccc3D
   bc2D=bbccc2D

   TLOOP : DO
      ! Integration step:
      DO j=1,NO_D3_BOX_STATES
         IF(D3STATETYPE(j).ge.0) THEN
            ccc_tmp3D(j,:) = bbccc3D(j,:) + delt*sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
         END IF
      END DO
      DO j=1,NO_D2_BOX_STATES
         IF(D2STATETYPE(j).ge.0) THEN
                  ccc_tmp2D(j,:) = bbccc2D(j,:) + delt*sum(D2SOURCE(j,:,:)-D2SINK(j,:,:),1)
         END IF
      END DO
      nmin=nmin+nstep 
      ! Check for negative concentrations
      min3D=minval(ccc_tmp3D)
      min2D=minval(ccc_tmp2D)
      IF (min3D.lt.eps.OR.min2D.lt.eps) THEN ! cut timestep
         IF (nstep.eq.1) THEN
            LEVEL1 'Necessary Time Step too small! Exiting...'
               blccc(1,:)=minloc(ccc_tmp3D)
               blccc(2,:)=minloc(ccc_tmp2D)
               LEVEL2 blccc
               LEVEL2 ccc_tmp3D(blccc(1,1),blccc(1,2)), &
                           bbccc3D(blccc(1,1),blccc(1,2))
               LEVEL2 ccc_tmp2D(blccc(2,1),blccc(2,2)), &
                           bbccc2D(blccc(2,1),blccc(2,2))
               LEVEL2 'EXIT at time: ',timesec
            STOP
         END IF
         nstep=nstep/2
         nmin=0
         D3STATE=bccc3D
         D2STATE=bccc2D
         dtm1=.5*delt
         delt=2.*nstep*mindelt
         timesec=maxdelt*ntime
         n=nmaxdelt/nstep
      ELSE
         ! filtering and advancement:
         write(6,*) 'filtered'
         DO j=1,NO_D3_BOX_STATES
            IF (D3STATETYPE(j).ge.0) THEN
               bbccc3D(j,:) = D3STATE(j,:) + &
                     ass*(bbccc3D(j,:)-2.*D3STATE(j,:)+ccc_tmp3D(j,:))
            END IF
         END DO
         DO j=1,NO_D2_BOX_STATES
            IF (D2STATETYPE(j).ge.0) THEN
               bbccc2D(j,:) = D2STATE(j,:) + &
                     ass*(bbccc2D(j,:)-2.*D2STATE(j,:)+ccc_tmp2D(j,:))
            END IF
         END DO
         D3STATE=ccc_tmp3D
         D2STATE=ccc_tmp2D
!        call ResetFluxes
!        call envforcing_bfm
         call EcologyDynamics
      END IF
      IF(nmin.eq.nmaxdelt) EXIT TLOOP
      IF(nmin.eq.0) THEN
         ! 2nd order approximation of backward State from Taylor expansion:
         bbccc2D=bc2D/n**2+D2STATE*(1.-1./n**2)+ &
            sum((D2SOURCE-D2SINK),2)*.5*delt*(1./n**2-1./n)
         bbccc3D=bc3D/n**2+D3STATE*(1.-1./n**2)+ &
            sum((D3SOURCE-D3SINK),2)*.5*delt*(1./n**2-1./n)
#ifdef DEBUG
         LEVEL2 'Time Step cut! delt= ',delt/2.,' nstep= ',nstep
#endif
      ELSE
#ifdef DEBUG
         LEVEL2 'Internal time step= ',nmin
#endif
         timesec=timesec+nstep*mindelt
#ifdef DEBUG
         LEVEL2 'Internal time= ',timesec
#endif
      END IF
   END DO TLOOP

   IF (nstep.ne.nmaxdelt) THEN 
      ! filter and forwarding of central value if the time step was cut:
#ifdef DEBUG
      LEVEL2 'Full Step Filter after cutting'
#endif
      bbccc3D=bccc3D
      bbccc2D=bccc2D
      DO j=1,NO_D3_BOX_STATES
         IF (D3STATETYPE(j).ge.0) THEN
            bbccc3D(j,:) = bccc3d(j,:) + &
               ass*(bc3D(j,:)-2.*bccc3d(j,:)+D3STATE(j,:))
         END IF
      END DO
      DO j=1,NO_D2_BOX_STATES
         IF (D2STATETYPE(j).ge.0) THEN
            bbccc2d(j,:) = bccc2d(j,:) + &
               ass*(bc2D(j,:)-2.*bccc2d(j,:)+D2STATE(j,:))
         END IF
      END DO
   ENDIF
   nstep=nmaxdelt
   nmin=0
   ntime=ntime+1
   dtm1=.5*delt
   delt=2*nstep*mindelt
   timesec=delt/2.*ntime
#ifdef DEBUG
   LEVEL1 'ntime: ',ntime
   LEVEL1 'Integration time: ',timesec
#endif

 END subroutine integrationLf
!EOC
!-----------------------------------------------------------------------
