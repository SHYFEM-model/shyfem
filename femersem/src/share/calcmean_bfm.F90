#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: compute the average field
!
! !INTERFACE:
   subroutine calcmean_bfm(mode)
!
! !DESCRIPTION:
!  This is a sophisticated subroutine to accumulate instantaneous
!  values and compute averages over each output interval.
!  It stores only the variables defined in the {/tt namelist} variable
!  {/tt ave_save}.
!
! !USES:
   use api_bfm
   use mem, only: D3STATE,D2STATE,D3DIAGNOS,D2DIAGNOS
   use mem, only: Depth,NO_BOXES,NO_BOXES_XY,make_flux_output
   implicit none
!
! !INPUT PARAMETERS:
   integer,intent(IN)                  :: mode
!
! !REVISION HISTORY:
!  Original author(s): Piet Ruardij (NIOZ)
!
! !LOCAL VARIABLES:
    integer                     ::i
    integer                     ::j
    integer                     ::k
    integer                     ::rc

!EOP
!-----------------------------------------------------------------------
!BOC

   select case (mode)
      case(INIT)   ! initialization
         i=count(var_ave(stPelStateS:stPelFluxE) .eqv. .true. )
         if ( i > 0 ) then
            allocate(D3ave(1:i,1:NO_BOXES),stat=rc)
            if (rc /= 0) stop 'init_bio(): Error allocating D3ave'
            D3ave=0.0
         endif
         i=count(var_ave(stBenStateS:stBenFluxE) .eqv. .true. )
         if ( i > 0 ) then
            allocate(D2ave(1:i,1:NO_BOXES_XY),stat=rc)
            if (rc /= 0) stop 'init_bio(): Error allocating D3ave'
            D2ave=0.0
         end if
         ave_count=0.0
      case(MEAN)    ! prepare for printing
         D3ave=D3ave/ave_count
         D2ave=D2ave/ave_count
         ave_count=0.0
      case(ACCUMULATE)   ! Start of new time-step
         ave_count=ave_count+1.0
         !---------------------------------------------
         ! Compute pelagic means
         !---------------------------------------------
         if (stPelStateE /= 0 ) then
            k=0
            j=0
            do i=stPelStateS,stPelStateE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D3ave(k,:)=D3STATE(j,:)
                  else
                     D3ave(k,:)=D3ave(k,:)+D3STATE(j,:)
                  end if
               end if
            end do
            j=0
            do i=stPelDiagS,stPelDiagE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D3ave(k,:)=D3DIAGNOS(j,:)
                  else
                     D3ave(k,:)=D3ave(k,:)+D3DIAGNOS(j,:)
                  end if
               end if
            end do
            j=0
            do i=stPelFluxS,stPelFluxE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  call make_flux_output(1,j,NO_BOXES,Depth,c1dim)
                  if ( ave_count < 1.5 ) then
                     D3ave(k,:)=c1dim
                  else
                     D3ave(k,:)=D3ave(k,:)+c1dim
                  end if
               end if
            end do
         else if (stBenStateE /= 0) then
            !---------------------------------------------
            ! Compute benthic means
            !---------------------------------------------
            k=0
            j=0
            do i=stBenStateS,stBenStateE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D2ave(k,0:1)=D2STATE(j,0:1)
                  else
                     D2ave(k,0:1)=D2ave(k,0:1)+D2STATE(j,0:1)
                  end if
               end if
            end do
            j=0
            do i=stBenDiagS,stBenDiagE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  if ( ave_count < 1.5 ) then
                     D2ave(k,:)=D2DIAGNOS(j,:)
                  else
                     D2ave(k,:)=D2ave(k,:)+D2DIAGNOS(j,:)
                  end if
               end if
            end do
            j=0
            do i=stBenFluxS,stBenFluxE
               j=j+1
               if ( var_ave(i) ) then
                  k=k+1
                  call make_flux_output(2,j,NO_BOXES,Depth,c1dim)
                  if ( ave_count < 1.5 ) then
                     D2ave(k,:)=c1dim(:)
                  else
                     D2ave(k,:)=D2ave(k,:)+c1dim(:)
                  end if
               end if
            end do
         end if
   end select

end subroutine calcmean_bfm
!EOC
