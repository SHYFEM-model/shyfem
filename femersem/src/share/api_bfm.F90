#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bfm
!
! !INTERFACE:
   MODULE api_bfm
!
! !DESCRIPTION: 
! API for the BFM. 
! Storage of variables and diagnostics
! To be used in all the coupled applications except
! GOTM, where it actually originated from.
! The GOTM module netcdfout is needed
! Appropriate functions are already available in the GOTM structure

!
! !USE:
   use mem, only: NO_D3_BOX_STATES, NO_BOXES,            &
                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,    &
                  NO_D2_BOX_STATES, NO_BOXES_XY,         &
                  NO_D2_BOX_DIAGNOSS, NO_D3_BOX_DIAGNOSS,&
                  NO_STATES, Depth, NO_D3_BOX_FLUX,      &
                  NO_D2_BOX_FLUX

   implicit none

!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bfm
!
! !PUBLIC DATA MEMBERS:
   logical                            :: bio_calc,bioshade_feedback
   integer                            :: bio_setup  =1
   integer                            :: surface_flux_method=-1
   integer                            :: n_surface_fluxes=-1
   character(len=PATH_MAX)            :: out_dir,out_fname,out_title
   integer                            :: out_units
   integer                            :: out_delta,out_secs

   !---------------------------------------------
   ! Dimension lengths for output
   !---------------------------------------------
   integer        :: lon_len
   integer        :: lat_len
   integer        :: depth_len
   integer        :: ocepoint_len,surfpoint_len,botpoint_len

   !---------------------------------------------
   ! BFM variable information for output
   !---------------------------------------------
   integer, dimension(:), allocatable    :: var_ids
   logical, dimension(:), allocatable    :: var_ave
   REALTYPE                              :: ave_count
   REALTYPE,allocatable,dimension(:,:)   :: D3ave
   REALTYPE,allocatable,dimension(:,:)   :: D2ave
   character(len=64), dimension(:), allocatable :: var_names
   character(len=64), dimension(:), allocatable :: var_units
   character(len=64), dimension(:), allocatable :: var_long
   !---------------------------------------------
   ! Indices of the various output variables
   !---------------------------------------------
   integer,public                            :: stPelStateS=0
   integer,public                            :: stPelDiagS=0
   integer,public                            :: stPelFluxS=0
   integer,public                            :: stBenStateS=0
   integer,public                            :: stBenDiagS=0
   integer,public                            :: stBenFluxS=0
   integer,public                            :: stPelStateE=0
   integer,public                            :: stPelDiagE=0
   integer,public                            :: stPelFluxE=0
   integer,public                            :: stBenStateE=0
   integer,public                            :: stBenDiagE=0
   integer,public                            :: stBenFluxE=0

   !---------------------------------------------
   ! Additional output variables
   !---------------------------------------------
   REALTYPE, dimension(:), allocatable   :: c1dim

!
! !REVISION HISTORY:
!  Author(s): Marcello Vichi and Piet Ruardij
!  Uses functions and portions of code from GOTM
!  by Hans Burchard and Karsten Bolding
!
! !LOCAL VARIABLES:
!
! !BUGS
!
!EOP
!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bfm module
!
! !INTERFACE:
   subroutine init_bfm(namlst)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  Adapted to BFM by Marcello Vichi
!
! !LOCAL VARIABLES:
   integer                   :: rc,i,j,n
   namelist /bfm_nml/ bio_calc,bio_setup,                  &
                      out_fname,out_dir,out_units,         &
                      out_title,out_delta,out_secs,        &
                      bioshade_feedback
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'init_bfm'

   !---------------------------------------------
   ! Provide sensible values for namelist parameters
   !---------------------------------------------
   bio_calc=.TRUE.
   bio_setup=1
   out_fname='bfm'
   out_dir='.'
   out_title='Another great BFM simulation!'
   out_units=0
   out_delta=100
   out_secs=100.0
   bioshade_feedback=.FALSE.

   !---------------------------------------------
   !  Open and read the namelist
   !---------------------------------------------
   open(namlst,file='bfm.nml',action='read',status='old',err=99)
   read(namlst,nml=bfm_nml,err=98)
   close(namlst)

   LEVEL2 "Writing NetCDF output to file: ",trim(out_fname)
   LEVEL3 "Output frequency every ",out_delta,"time-steps"

   select case (bio_setup)
      case (0)
      case (1) ! Pelagic only
        LEVEL2 "Using a Pelagic setup (bio_setup=1)"
        LEVEL3 'pelagic variables =',NO_D3_BOX_STATES
        LEVEL3 'pelagic transported variables ='
        LEVEL3 'pelagic diagnostic variables =', NO_D3_BOX_DIAGNOSS
      case (2) ! Benthic only
        LEVEL2 "Using a Benthic-only setup (bio_setup=2)"
        LEVEL3 'benthic variables =',NO_D2_BOX_STATES
        LEVEL3 'benthic diagnostic variables=', NO_D2_BOX_DIAGNOSS
      case (3) ! Pelagic-Benthic coupling
        LEVEL2 "Using a Pelagic-Benthic coupled setup (bio_setup=3)"
        LEVEL3 'pelagic variables =',NO_D3_BOX_STATES
        LEVEL3 'pelagic transported variables ='
        LEVEL3 'pelagic diagnostic variables =', NO_D3_BOX_DIAGNOSS
        LEVEL3 'benthic variables =',NO_D2_BOX_STATES
        LEVEL3 'benthic diagnostic variables=', NO_D2_BOX_DIAGNOSS
   end select

   LEVEL2 'Dimensional informations:'
   LEVEL3 'NO_BOXES_X=',NO_BOXES_X
   LEVEL3 'NO_BOXES_Y=',NO_BOXES_Y
   LEVEL3 'NO_BOXES_Z=',NO_BOXES_Z
   LEVEL3 'NO_BOXES=',NO_BOXES
   LEVEL3 'NO_BOXES_XY=',NO_BOXES_XY
   LEVEL3 'NO_STATES=',NO_STATES
   LEVEL3 'Step 1 of BFM initialisation done ...'
   ! dimension lengths used in the netcdf output
   lon_len = NO_BOXES_X
   lat_len = NO_BOXES_Y
   depth_len = NO_BOXES_Z
   ocepoint_len = NO_BOXES
   surfpoint_len = NO_BOXES_XY
   botpoint_len = NO_BOXES_XY

   !---------------------------------------------
   ! Allocate arrays with attributes of state
   ! variables
   !---------------------------------------------
   ! total number of output states
   n = NO_D3_BOX_STATES+NO_D3_BOX_FLUX+NO_D3_BOX_DIAGNOSS+ &
      NO_D2_BOX_STATES+NO_D2_BOX_FLUX+NO_D2_BOX_DIAGNOSS
   allocate(var_ids(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_ids'
   var_ids=0;
   allocate(var_ave(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_ave'
   var_ave=.false.;
   allocate(var_names(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_names'
   allocate(var_units(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_units'
   allocate(var_long(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_long'
   ! temporary diagnostic variable
   allocate(c1dim(0:NO_BOXES),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating c1dim'

   return
98 FATAL 'I could not read bfm.nml'
   stop 'init_bfm'
99 LEVEL2 'I could not open bfm.nml'
   LEVEL2 'Simulation starting without the BFM'
   bio_calc = .false.

  end subroutine init_bfm
!EOC


!-----------------------------------------------------------------------

   END MODULE api_bfm

!-----------------------------------------------------------------------

