!$Id: $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_bfm --- BFM bio model \label{sec:bio_bfm}
!
! !INTERFACE:
   module bio_bfm
!
! !DESCRIPTION:
!
!
! !USES:
!  default: all is private.
   use bio_var
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio_bfm, pointers_gotm_bfm,            &
          var_info_bfm, envforcing_bfm, do_bio_bfm,   &
          allocate_memory_bfm,reset_diagonal,         &
          test_on_negative_states, end_bio_bfm
!
! !PRIVATE DATA MEMBERS:
   REALTYPE,dimension(:),allocatable :: cdepth
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!
!  $Log: $
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the template bio module
!
! !INTERFACE:
   subroutine init_bio_bfm(nlev,out_unit)
!
! !DESCRIPTION:
!  Here, the main communication of array dimensions between GOTM
!  and BFM is done.
!
!
! !USES:
   use mem, only: NO_D3_BOX_STATES, NO_BOXES,          &
                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,  &
                  NO_D2_BOX_STATES, NO_BOXES_XY,       &
                  NO_D2_BOX_DIAGNOSS, NO_D3_BOX_DIAGNOSS,&
                  NO_D2_BOX_FLUX, NO_D3_BOX_FLUX,&
                  NO_STATES

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)   :: nlev
   integer,          intent(in)   :: out_unit
!
   integer :: i,rc
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_bfm'


   ! BFM  --> GOTM
   numc  = NO_D3_BOX_STATES
   numbc = NO_D2_BOX_STATES
   numc_diag  = NO_D3_BOX_DIAGNOSS
   numbc_diag = NO_D2_BOX_DIAGNOSS
   numc_flux  = NO_D3_BOX_FLUX
   numbc_flux = NO_D2_BOX_FLUX
   ! numcc is the number of transported variables
   numcc = numc

   ! GOTM --> BFM
   NO_BOXES_X  = 1
   NO_BOXES_Y  = 1
   NO_BOXES_Z  = nlev
   NO_BOXES    = NO_BOXES_X * NO_BOXES_Y * NO_BOXES_Z
   NO_BOXES_XY = NO_BOXES_X * NO_BOXES_Y
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES +   &
                 NO_D2_BOX_STATES * NO_BOXES_XY
   !LOGUNIT = out_unit

   LEVEL3 'pelagic variables =',numc
   LEVEL3 'pelagic transported variables =',numcc
   LEVEL3 'benthic variables =',numbc
   LEVEL3 'pelagic variables prepared for output',numc_diag
   LEVEL3 'benthic variables prepared for output',numbc_diag
   LEVEL3 'NO_BOXES_X=',NO_BOXES_X
   LEVEL3 'NO_BOXES_Y=',NO_BOXES_Y
   LEVEL3 'NO_BOXES_Z=',NO_BOXES_Z
   LEVEL3 'NO_BOXES=',NO_BOXES
   LEVEL3 'NO_BOXES_XY=',NO_BOXES_XY
   LEVEL3 'NO_STATES=',NO_STATES
   LEVEL3 'Step 1 of GOTM <-> BFM initialisation done ...'

   sfl=_ZERO_
   sfl_read=_ZERO_
   allocate(cdepth(1:NO_BOXES),stat=rc)
   return

   end subroutine init_bio_bfm
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialize BFM and GETM shared memory
!
! !INTERFACE:
   subroutine pointers_gotm_bfm()
!
! !DESCRIPTION:
! Allocate pointers to GOTM memory
!
! !USES:
   use mem, only: D3STATE,D3SOURCE,D3SINK,D3STATETYPE, &
                  D3DIAGNOS,D2STATE,D2SOURCE,D2SINK,   &
                  D2STATETYPE,NO_BOXES,NO_BOXES_XY,    &
                  D2DIAGNOS,NO_D2_BOX_STATES,          &
                  NO_D2_BOX_DIAGNOSS

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:

   !---------------------------------------------
   ! Pelagic pointers
   !---------------------------------------------
   D3STATE  => cc(:,1:NO_BOXES)
   D3SOURCE => pp(:,:,1:NO_BOXES)
   D3SINK   => dd(:,:,1:NO_BOXES)
   D3STATETYPE => pelvar_type
   if (numc_diag > 0) D3DIAGNOS => diag(:,1:NO_BOXES)

   !---------------------------------------------
   ! Benthic pointers
   !---------------------------------------------
   if (bio_setup >=2 ) then
      D2STATE  => ccb(:,1:NO_BOXES_XY)
      D2SOURCE => ppb(:,:,1:NO_BOXES_XY)
      D2SINK   => ddb(:,:,1:NO_BOXES_XY)
      D2STATETYPE => benvar_type
      if (numbc_diag>0) D2DIAGNOS => diagb(:,1:NO_BOXES_XY)
   else
      ! allocate memory anyhow to avoid problems with BFM allocation
      allocate(D2STATE(1:NO_D2_BOX_STATES,1:NO_BOXES_XY))
      allocate(D2SOURCE(1:NO_D2_BOX_STATES,1:NO_D2_BOX_STATES,1:NO_BOXES_XY))
      allocate(D2SINK(1:NO_D2_BOX_STATES,1:NO_D2_BOX_STATES,1:NO_BOXES_XY))
      allocate(D2STATETYPE(1:NO_D2_BOX_STATES ))
      if (numbc_diag>0)  &
         allocate(D2DIAGNOS(1:NO_D2_BOX_DIAGNOSS,1:NO_BOXES_XY))
   end if

   end subroutine pointers_gotm_bfm
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Providing info on variables
!
! !INTERFACE:
   subroutine var_info_bfm()
!
! !DESCRIPTION:
!  This subroutine provides information on the variables. To be used
!  when storing data in NetCDF files.
!
! !USES:
   use mem
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call set_var_info_bfm
   return
   end subroutine var_info_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light and other environmental forcing used in the BFM
!
! !INTERFACE
   subroutine envforcing_bfm(nlev,h,t,s,I_0,wind_gotm, &
                             bioshade_feedback,bioshade,abioshade)
!
! !DESCRIPTION
!
! !USES
! BFM modules
use constants, ONLY: E2W
use mem_Param, ONLY: p_eps0, p_epsESS, p_PAR,p_small
use mem,       ONLY: NO_BOXES, R6c, PhytoPlankton, xEPS, ESS, &
                     iiPhytoPlankton, iiL, Chla, ETW, ESW, Wind,    &
                     Depth, EIR, ABIO_eps
use mem_Param,  ONLY: p_eps0, p_epsESS,p_poro

IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                              :: nlev
   REALTYPE, intent(in)                 :: h(0:nlev)
   REALTYPE, intent(in)                 :: t(0:nlev)
   REALTYPE, intent(in)                 :: s(0:nlev)
   REALTYPE, intent(in)                 :: I_0
   REALTYPE, intent(in)                 :: wind_gotm
   logical, intent(in)                  :: bioshade_feedback
   REALTYPE, intent(in),optional        :: abioshade(0:nlev)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: bioshade(0:nlev)!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer             :: i,n
   REALTYPE            :: psilt
!EOP
!-----------------------------------------------------------------------
!BOC

!   LEVEL2 'calculating environmental forcings for the BFM'

   !---------------------------------------------
   ! Assign depths of layers
   ! temperature and salinity
   !---------------------------------------------
   Depth(:) = h(1:nlev)
    cdepth(NO_BOXES) =depth(NO_BOXES)
    Wind=wind_gotm
    do n=NO_BOXES-1,1,-1
       cdepth(n)=depth(n)+ cdepth(n+1)
    enddo
   ETW(:) = t(1:nlev)
   ESW(:) = s(1:nlev)
   psilt=(p_poro(1) - 0.38662 )/ 0.00415
   ESS(:) = 1250.0 * min(1.0,20.0/cdepth(1))* psilt/7.0

   !---------------------------------------------
   ! Compute extinction coefficient
   !---------------------------------------------

   if (p_eps0 ==0.0 ) then
     ABIO_eps(:) = abioshade(1:nlev)
   end if

   call  CalcVerticalExtinction( )

   !---------------------------------------------
   ! Notice that irradiance in the BFM is in
   ! uE/m2/s and is defined at the top of each
   ! layer (the derivation of the middle-layer
   ! EIR for production is done in the
   ! Phytoplankton routines)
   !---------------------------------------------
   EIR(nlev) = max(p_small,p_PAR*I_0/E2W)
   do i=nlev,2,-1
     EIR(i-1) = EIR(i)*exp(-xEPS(i)*Depth(i))
   end do

   !---------------------------------------------
   ! bioshade is instead derived in the
   ! middle of the layer and it's non-dimensional
   !---------------------------------------------
   if (bioshade_feedback) &
     bioshade(1:nlev) =  EIR(:)*exp(-xEPS(:)*Depth(:)*0.5)/ EIR(nlev)

   end subroutine envforcing_bfm
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of the BFM model
!
! !INTERFACE
   subroutine do_bio_bfm(first)
!
! !DESCRIPTION
!  This subroutine is a wrapper for the computing core of the BFM
!
! !USES
   use mem_param, only: AssignAirPelFluxesInBFMFlag,AssignPelBenFluxesInBFMFlag
   use mem, only: sediPI, sediR6, iiC,iiN,iiP,iiS,iiL, &
                  ppR6c, ppR6n, ppR6p, ppR6s, NO_BOXES_Z,   &
                  ppR1c, ppR1n, ppR1p,   &
                  ppO2o,ppN1p,ppN3n,ppN4n,ppN5s,ppN6r,  &
                  NO_D3_BOX_STATES, Depth, jOAO2o,                &
                  ppPhytoPlankton,iiPhytoPlankton, &
                  jG2O2o,jK1N1p,jK3N3n,jK4N4n,jK5N5s,jK6N6r, &
                  retPIc,retPIn,retPIp,retPIs,retPIl,  &
                  rutQ6c,rutQ6n,rutQ6p,rutQ6s,  &
                  rutQ1c,rutQ1n,rutQ1p
   use constants,  only: SEC_PER_DAY
   use mem_settling, only: p_burvel
   use gotm_error_msg, only:gotm_error

   IMPLICIT NONE
!
   logical                     :: first
   integer                     :: n,k,i
   REALTYPE                    :: burvel, topm3psec
!See "USE association" above
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from template by Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Compute BFM terms
   !---------------------------------------------
   call EcologyDynamics

   !---------------------------------------------
   ! Transfer sinking velocities (m/d -> m/s)
   !---------------------------------------------
   if ( bio_setup ==2 ) return
   do i=1,iiPhytoPlankton
     k=ppPhytoPlankton(i,iiC)
     ws(k,1:NO_BOXES_Z) = -sediPI(i,1:NO_BOXES_Z)/SEC_PER_DAY
     ws(k,0)= ws(k,1);
     k=ppPhytoPlankton(i,iiN)
     ws(k,1:NO_BOXES_Z) = -sediPI(i,1:NO_BOXES_Z)/SEC_PER_DAY
     ws(k,0)= ws(k,1);
     k=ppPhytoPlankton(i,iiP)
     ws(k,1:NO_BOXES_Z) = -sediPI(i,1:NO_BOXES_Z)/SEC_PER_DAY
     ws(k,0)= ws(k,1);
     k=ppPhytoPlankton(i,iiL)
     ws(k,1:NO_BOXES_Z) = -sediPI(i,1:NO_BOXES_Z)/SEC_PER_DAY
     ws(k,0)= ws(k,1);
     k=ppPhytoPlankton(i,iiS)
     if ( k > 0 ) then
       ws(k,1:NO_BOXES_Z) = -sediPI(i,1:NO_BOXES_Z)/SEC_PER_DAY
       ws(k,0)= ws(k,1);
     endif
   enddo

   ws(ppR6c,1:NO_BOXES_Z) = -sediR6/SEC_PER_DAY
   ws(ppR6c,0) = ws(ppR6c,1)
   ws(ppR6n,1:NO_BOXES_Z) = -sediR6/SEC_PER_DAY
   ws(ppR6n,0) = ws(ppR6n,1)
   ws(ppR6p,1:NO_BOXES_Z) = -sediR6/SEC_PER_DAY
   ws(ppR6p,0) = ws(ppR6p,1)
   ws(ppR6s,1:NO_BOXES_Z) = -sediR6/SEC_PER_DAY
   ws(ppR6s,0) = ws(ppR6s,1)

   !---------------------------------------------
   ! Surface fluxes
   !---------------------------------------------
   topm3psec=_ONE_/Depth(NO_BOXES_Z)/ SEC_PER_DAY
   if ( .NOT. AssignAirPelFluxesInBFMFlag ) then
     select case (surface_flux_method)
        case (-1)! absolutely nothing
        case (0) ! constant
           sfl(ppN3n) =   0.09  *topm3psec
           sfl(ppN4n) =   0.10  *topm3psec
           sfl(ppN1p) =   0.0  !0.0
           sfl(ppO2o) =   jOAO2o(1) *topm3psec
        case (2) ! from file via sfl_read
           ! fluxes are in mmol m-2 d-1
           sfl(ppN3n) =   1.0*sfl_read(1)/SEC_PER_DAY
           sfl(ppN4n) =   1.0*sfl_read(2)/SEC_PER_DAY
           sfl(ppN1p) =0.0
!          sfl(ppN1p) =   1.0*sfl_read(3)/SEC_PER_DAY
        case (3) ! sfl array filled externally - for 3D models
        case default
     end select
   endif

   !---------------------------------------------
   ! Bottom fluxes
   !---------------------------------------------
   topm3psec=1.0/Depth(1)/ SEC_PER_DAY
   if (bio_setup == 3 .and. ( .NOT.AssignPelBenFluxesInBFMFlag)) then

      bfl(ppR6c) = ( -rutQ6c(1))*topm3psec
      bfl(ppR6n) = ( -rutQ6n(1))*topm3psec
      bfl(ppR6p) = ( -rutQ6p(1))*topm3psec
      bfl(ppR6s) = ( -rutQ6s(1))*topm3psec

      bfl(ppR1c) =  -rutQ1c(1)*topm3psec
      bfl(ppR1n) =  -rutQ1n(1)*topm3psec
      bfl(ppR1p) =  -rutQ1p(1)*topm3psec

      bfl(ppO2o) = jG2O2o(1)*topm3psec
      bfl(ppN1p) = jK1N1p(1)*topm3psec
      bfl(ppN3n) = jK3N3n(1)*topm3psec
      bfl(ppN4n) = jK4N4n(1)*topm3psec
      bfl(ppN5s) = jK5N5s(1)*topm3psec
      bfl(ppN6r) = jK6N6r(1)*topm3psec

      do i=1,iiPhytoPlankton
        bfl(ppPhytoPlankton(i,iiC)) = -retPIc(i,1)*topm3psec
        bfl(ppPhytoPlankton(i,iiN)) = -retPIn(i,1)*topm3psec
        bfl(ppPhytoPlankton(i,iiP)) = -retPIp(i,1)*topm3psec
        bfl(ppPhytoPlankton(i,iiL)) = -retPIl(i,1)*topm3psec
        k=ppPhytoPlankton(i,iiS)
        if ( k > 0 ) bfl(k) = -retPIs(i,1)*topm3psec
      enddo

!MAV check this
!     burvel= - p_burvel/SEC_PER_DAY
!     k=0
!     burvel=0.08 * cdepth(1)/SEC_PER_DAY
!     do n=1,NO_D3_BOX_STATES
!        if ( ws(n,NO_BOXES) /=0.0 ) then
!           ws(n,:)=-min(burvel,-ws(n,:))
!        endif
!     enddo

   endif

   return

   end subroutine do_bio_bfm
!EOC
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Reset diagonal id a 3d array
!
! !INTERFACE:
   subroutine reset_diagonal(n,pp)
!
! !DESCRIPTION:
!    Reset of the diagonal
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                    :: n
   REALTYPE,dimension(:,:,:),intent(inout) :: pp
!
! !REVISION HISTORY:
!  Original author(s): Piet Ruardij
!
! !LOCAL VARIABLES:
   integer                   :: i
!EOP
!-----------------------------------------------------------------------
!BOC
     do i=1,n
       pp(i,i,:) = _ZERO_
     end do

   return
   end subroutine reset_diagonal
!EOC
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocate_bfm
!
! !INTERFACE:
        subroutine allocate_memory_bfm(nlev)
!
! !INPUT PARAMETERS:
        implicit none
        integer,intent(IN)            ::nlev
!
! !LOCAL VARAIBELS:
   integer                   :: rc
!
! !DESCRIPTION:
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  28-04-2006  Piet Ruardij Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

   if ( numc_diag > 0 ) then
     allocate(diag(1:numc_diag,0:nlev),stat=rc)
     if (rc /= 0) STOP 'init_bio: Error allocating (cc)'


     diag=_ZERO_                                                     !BFM
   endif


   if (bio_setup >= 2) then                                         !BFM
     ! allocate benthic state variables                             !BFM
     allocate(ccb(1:numbc,0:1),stat=rc)                             !BFM
     if (rc /= 0) STOP 'init_bio: Error allocating (ccb)'           !BFM
     allocate(ppb(1:numbc,1:numbc,0:1),stat=rc)                     !BFM
     if (rc /= 0) STOP 'init_bio: Error allocating (ppb)'           !BFM
     allocate(ddb(1:numbc,1:numbc,0:1),stat=rc)                     !BFM
     if (rc /= 0) STOP 'init_bio: Error allocating (ppb)'           !BFM

     ccb=_ZERO_                                                     !BFM
     ppb=_ZERO_                                                     !BFM
     ddb=_ZERO_                                                     !BFM

     ! allocate variable holding type and save attributes           !BFM
     allocate(benvar_type(1:numbc),stat=rc)                         !BFM
     if (rc /= 0) STOP 'init_bio: Error allocating (benvar_type)'   !BFM
     benvar_type = 0

     if ( numbc_diag > 0 ) then
       allocate(diagb(1:numbc_diag,0:1),stat=rc)
       if (rc /= 0) STOP 'init_bio: Error allocating (cc)'

       diagb=_ZERO_                                                     !BFM
     endif

   end if



 end subroutine allocate_memory_bfm
!EOC
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Test negative concentrations
!
! !INTERFACE:
       subroutine test_on_negative_states ( statenr,nlev, after, error )
!
! !USES:
       use gotm_error_msg, only:set_warning_for_getm
       IMPLICIT NONE
!
! !INPUT PARAMETERS:
       integer,intent(in)                      :: statenr
       integer,intent(in)                      :: nlev
       character(len=*),intent(IN)             :: after
!
! !OUTPUT PARAMETERS:
       integer,intent(OUT)                     :: error
!          Array cldim is modified if ncecessry
!
! !LOCAL VARAIBELS:
        integer              ::k
        integer              ::i

! !DESCRIPTION:
!   Routine to check for negative values.
!   Negative values are corrected with the aveage of neighbouring
!   grid points. A warning is given.
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!      created by P. Ruardij 21-06-2006
!
!
!EOP
!-------------------------------------------------------------------------
!BOC
       error=0
       if (minval(c1dim(:)) < 0.00D+00) then                    !BFM
          STDERR "statenr=",statenr, c1dim(:)
          STDERR "Negative value after call to ",after         !BFM
          k=0                                                  !BFM
          do i = 0,nlev                                        !BFM
              if ( c1dim(i).lt.0.0D+00) then                   !BFM
                  k=i                                          !BFM
                  STDERR i,c1dim(i)                            !BFM
                  call set_warning_for_getm()
                  if ( i == 1 ) then
                     if ( c1dim(i+1) > 0.0 )  then
                       c1dim(i)=0.5* c1dim(i+1)
                       k=0
                     endif
                  elseif ( i == nlev ) then
                     if ( c1dim(i-1) > 0.0 )  then
                       c1dim(i)=0.5* c1dim(i-1)
                       k=0
                     endif
                  else if ( c1dim(i-1) > 0.0 .and. c1dim(i+1)>0.0 ) then
                     k=0
                     c1dim(i)=(c1dim(i-1)+c1dim(i+1)) * 0.5
                  else if ( c1dim(i-1) > 0.0 ) then
                       c1dim(i)=0.5* c1dim(i-1)
                       k=0
                  else if ( c1dim(i+1) > 0.0 ) then
                       c1dim(i)=0.5* c1dim(i+1)
                       k=0
                  endif
              endif                                            !BFM
          end do                                               !BFM
          if (k.ne.0)  error=k                                 !BFM
       end if

     end subroutine test_on_negative_states


!EOC
!-------------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio_bfm
!
! !DESCRIPTION:
!  Nothing done here --- supplied for completeness
!  with GOTM bio structure.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine end_bio_bfm
!EOC

!-----------------------------------------------------------------------


   end module bio_bfm

!-----------------------------------------------------------------------
! Copyright by the GOTM-team and BFM-team under the GNU Public License
! www.gnu.org
!-----------------------------------------------------------------------
