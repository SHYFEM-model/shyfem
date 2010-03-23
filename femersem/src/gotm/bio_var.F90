!$Id: bio_var.F90,v 1.7 2005-12-02 20:57:27 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_var --- declaration of biological variables
!
! !INTERFACE:
   module bio_var
!
! !DESCRIPTION:
!  Here all variables necessary for the biogeochemical models are
!  declared, mostly as allocatable variables.
!
! !USES:
!  default: all is public.
   public
!
! !PUBLIC DATA MEMBERS:
   integer                               :: bio_model
   integer                               :: numc,numcc
   REALTYPE                              :: I_0
   REALTYPE, dimension(:), allocatable   :: zlev
   REALTYPE, dimension(:), allocatable   :: par
   REALTYPE, dimension(:,:), allocatable,target :: cc,ws    !BFM
   integer                               :: surface_flux_method=-1
   integer                               :: n_surface_fluxes=-1
   REALTYPE, dimension(:), allocatable   :: sfl_read
   REALTYPE, dimension(:), allocatable   :: sfl,bfl
   integer, dimension(:), allocatable    :: posconc
   logical, dimension(:), allocatable    :: mussels_inhale
   logical, dimension(:,:), allocatable  :: particle_active
   integer, dimension(:,:), allocatable  :: particle_indx
   REALTYPE, dimension(:,:), allocatable :: particle_pos

   integer, dimension(:), allocatable    :: var_ids
   character(len=64), dimension(:), allocatable :: var_names
   character(len=64), dimension(:), allocatable :: var_units
   character(len=64), dimension(:), allocatable :: var_long

   REALTYPE, parameter                   :: secs_pr_day=86400.

   integer                               :: bio_setup =1        !BFM
   integer,public                        :: pelvar_save_all=0   !BFM
   REALTYPE, dimension(:), allocatable   :: c1dim

#ifdef BFM_GOTM
   ! additional storage variables for benthic and diagnostics
   integer              :: numbc
   integer              :: numc_diag,numbc_diag
   integer              :: numc_flux,numbc_flux
   ! Start and End markers for variable and diagnostics storage
   integer              :: stPelStateS=1
   integer              :: stPelDiagS=1
   integer              :: stPelFluxS=1
   integer              :: stBenStateS=1
   integer              :: stBenDiagS=1
   integer              :: stBenFluxS=1
   integer              :: stPelStateE=0
   integer              :: stPelDiagE=0
   integer              :: stPelFluxE=0
   integer              :: stBenStateE=0
   integer              :: stBenDiagE=0
   integer              :: stBenFluxE=0

   ! parameter values for the attributes pelvar_type and benvar_type
   integer, parameter   :: SINKSOURCE=-1
   integer, parameter   :: NOTRANSPORT=0
   integer, parameter   :: HORTRANSPORT=10
   integer, parameter   :: ALLTRANSPORT=20

   !additional BFM pelagic arrays
   REALTYPE, dimension(:),     allocatable         ::  SSt,RRa
   REALTYPE, dimension(:,:,:), allocatable, target ::  dd,pp

   !benthic BFM arrays
   REALTYPE, dimension(:,:),   allocatable, target :: ccb
   REALTYPE, dimension(:,:,:), allocatable, target :: ddb,ppb

   !diagnostic output arrays  (pelagic and benthic)
   REALTYPE, dimension(:,:),   allocatable, target :: diag
   REALTYPE, dimension(:,:),   allocatable, target :: diagb

   ! type and save attributes of pelagic and benthic variables
   integer, dimension(:),      allocatable, target :: pelvar_type
   integer, dimension(:),      allocatable, target :: benvar_type

   ! store arrays for average computations (pelagic and benthic)
   logical , dimension(:)  ,   allocatable           :: var_ave
   REALTYPE, dimension(:,:),   allocatable, target   :: cc_ave
   REALTYPE, dimension(:,:),   allocatable, target   :: ccb_ave
   REALTYPE                                          :: ave_count
#endif !BFM

!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  BFM additions:      Piet Ruardij & Marcello Vichi
!
!  $Log: bio_var.F90,v $
!  Revision 1.7  2005-12-02 20:57:27  hb
!  Documentation updated and some bugs fixed
!
!  Revision 1.6  2005-11-17 09:58:18  hb
!  explicit argument for positive definite variables in diff_center()
!
!  Revision 1.5  2004/07/30 09:22:20  hb
!  use bio_var in specific bio models - simpliefied internal interface
!
!  Revision 1.4  2004/03/30 11:32:48  kbk
!  select between eulerian or lagrangian solver
!
!  Revision 1.3  2003/10/16 15:42:16  kbk
!  simple mussesl model implemented - filter only
!
!  Revision 1.2  2003/09/16 12:11:24  hb
!  added new biological model - bio_iow
!
!  Revision 1.1  2003/07/23 12:27:31  hb
!  more generic support for different bio models
!
!
!EOP
!-----------------------------------------------------------------------

   end module bio_var

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
