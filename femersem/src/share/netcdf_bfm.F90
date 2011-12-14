#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: netcdf_bfm --- Save the BFM results in NetCDF
!
! !INTERFACE:
   module netcdf_bfm
!
! !DESCRIPTION:
!  This module provides routines for saving the results using
!  NetCDF format.
!
! !USES:
   use api_bfm, only: out_dir,out_fname,out_title,out_units,out_delta,out_secs
   use api_bfm, only: stPelStateS,stPelDiagS,stPelFluxS,stBenStateS,stBenDiagS,stBenFluxS
   use api_bfm, only: stPelStateE,stPelDiagE,stPelFluxE,stBenStateE,stBenDiagE,stBenFluxE
   use api_bfm, only: lon_len,lat_len,ocepoint_len,surfpoint_len,botpoint_len,depth_len
   use api_bfm, only: bio_setup,var_ids,var_names,var_long,var_units,c1dim
   use api_bfm, only: D3ave,D2ave,var_ave,ave_count
   use mem,     only: NO_BOXES,NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z,NO_BOXES_XY,Depth
   use mem,     only: make_flux_output
   implicit none

!
! !PUBLIC MEMBER FUNCTIONS:
   public init_netcdf_bfm,init_save_bfm,save_bfm,close_ncdf,check_err
   public define_mode,new_nc_variable,set_attributes,store_data
!
! !PUBLIC DATA MEMBERS:
   !---------------------------------------------
   ! netCDF file specifications
   !---------------------------------------------
   integer,public                :: ncid_bfm
   integer,public                :: ncdf_time_unit
   ! record counter
   integer,public                :: recnum = 0
   integer                       :: time_len
   !---------------------------------------------
   ! Dimension IDs
   !---------------------------------------------
   integer                            :: lon_dim,lat_dim,depth_dim
   integer                            :: ocepoint_dim
   integer                            :: surfpoint_dim,botpoint_dim
   integer                            :: time_dim
   integer, parameter                 :: dim1=1,dim4=4
   integer                            :: dims(dim4)
   !---------------------------------------------
   ! Coordinate variables IDs
   !---------------------------------------------
   integer          :: lon_id,lat_id,z_id,z1_id,time_id
   integer          :: zeta_id
   integer          :: depth_id,ocepoint_id,surfpoint_id,botpoint_id

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  Modifications: Marcello Vichi
!
!EOP
!
! !PRIVATE DATA MEMBERS
!  variable ids
   integer, private          :: start(4),edges(4)
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the netcdf output
!
! !INTERFACE:
   subroutine init_netcdf_bfm(title,start_time,time_unit,lat,lon,z,dz, &
                              lat2d,lon2d,oceanpoint,surfacepoint,bottompoint)
!
! !DESCRIPTION:
!  Prepare the netcdf output file which is finalized in init_save_bfm
!
! !USES:
   implicit none
!
! !INPUT/OUTPUT PARAMETERS:
   character(len=*), intent(in)                 :: title,start_time
   integer, intent(in)                          :: time_unit
   REALTYPE, intent(in),optional                :: lat,lon
   REALTYPE, intent(in),dimension(:,:),optional :: lat2d,lon2d
   REALTYPE, intent(in),dimension(:),optional   :: z,dz
   integer, intent(in),dimension(:),optional   :: oceanpoint
   integer, intent(in),dimension(:),optional   :: surfacepoint,bottompoint
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  Modifications: Marcello Vichi
!
!EOP
!
! !LOCAL VARIABLES:
   character(len=PATH_MAX)   :: ext,fname
   integer                   :: iret
   character(len=128)        :: ncdf_time_str,history
!  dimension lengths (not used yet)
   integer                   :: lon_len
   integer                   :: lat_len
   integer                   :: depth_len
   integer                   :: time_len
!!
!-------------------------------------------------------------------------
!BOC
   LEVEL1 'init_netcdf_bfm'

   !---------------------------------------------
   ! Prepare the netcdf file
   !---------------------------------------------

   end subroutine init_netcdf_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Intialise the storage of results in NetCDF
!
! !INTERFACE:
   subroutine init_save_bfm()
!
! !DESCRIPTION:
! Preparation of the netcdf output.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:

!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  Adapted to BFM: Marcello Vichi (INGV) & Piet Ruardij (NIOZ)
!
! !LOCAL VARIABLES:
   integer, save             :: nn       ! number pel.var to be saved 
   integer, save             :: nnb      ! number ben.var to be saved 
   integer                   :: iret,rc
   integer                   :: out_unit=67
   integer                   :: i,j,n
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'init_save_bfm'
   !---------------------------------------------
   ! Enter define mode
   !---------------------------------------------
   return
   end subroutine init_save_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Store the results
!
! !INTERFACE:
  subroutine save_bfm(time)
!
! !DESCRIPTION:
! output of BFM variables 
!
! !USES:
   use mem, only: D3STATE,D3DIAGNOS,D2STATE,D2DIAGNOS
   implicit none
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)     :: time
! !LOCAL VARIABLES:
   integer                   :: iret
   integer                   :: i,j,k,n
   REALTYPE                  :: temp_time

   return
   end subroutine save_bfm
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Close files used for saving model results
!
! !INTERFACE:
   subroutine close_ncdf(ncid)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Closes the NetCDF file.
!
! !INPUT PARAMETERS:
   integer, intent(in)       :: ncid
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: iret
!
!-------------------------------------------------------------------------
!BOC
   LEVEL1 'Output has been written in NetCDF'


   return
   end subroutine close_ncdf
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Begin or end define mode
!
! !INTERFACE:
   integer function define_mode(ncid,action)
!
! !DESCRIPTION:
!  Depending on the value of the argument {\tt action},
!  this routine put NetCDF in the `define' mode or not.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)       :: ncid
   logical, intent(in)       :: action
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer         :: iret
!
   define_mode = 0
!-----------------------------------------------------------------------
!BOC
   end function define_mode
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define a new NetCDF variable
!
! !INTERFACE:
   integer function new_nc_variable(ncid,name,data_type,n,dims,id)
!
! !DESCRIPTION:
!  This routine is used to define a new variable to store in a NetCDF file.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: ncid
   character(len=*), intent(in)        :: name
   integer, intent(in)                 :: data_type,n
   integer, intent(in)                 :: dims(:)
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: id
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: iret
!
!-----------------------------------------------------------------------
!BOC
   id = 0
   new_nc_variable = 0
   return
   end function new_nc_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set attributes for a NetCDF variable.
!
! !INTERFACE:
   integer function set_attributes(ncid,id,                         &
                                   units,long_name,                 &
                                   valid_min,valid_max,valid_range, &
                                   scale_factor,add_offset,         &
                                   FillValue,missing_value,         &
                                   C_format,FORTRAN_format,         &
                                   compress,formula_term)
!
! !DESCRIPTION:
!  This routine is used to set a number of attributes for
!  variables. The routine makes heavy use of the {\tt optional} keyword.
!  The list of recognized keywords is very easy to extend. 
!  The CF-1.0 convention is used.
!
! !USES:
!  IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: ncid,id
   character(len=*), optional          :: units,long_name
   REALTYPE, optional                  :: valid_min,valid_max
   REALTYPE, optional                  :: valid_range(2)
   REALTYPE, optional                  :: scale_factor,add_offset
   REALTYPE, optional                  :: FillValue,missing_value
   character(len=*), optional          :: C_format,FORTRAN_format
   character(len=*), optional          :: compress,formula_term
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!
! !LOCAL VARIABLES:
   integer                   :: len,iret
   REAL_4B                   :: vals(2)
!
!EOP
!-----------------------------------------------------------------------
!BOC
   set_attributes = 0
   return
   end function set_attributes
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Store values in a NetCDF file
!
! !INTERFACE:
   integer function store_data(ncid,id,var_shape,nbox, &
                               iscalar,iarray,scalar,array,garray)
!
! !DESCRIPTION:
!  This routine is used to store a variable in the NetCDF file.
!  The subroutine uses {\tt optional} parameters to find out which data
!  type to save.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: ncid,id,var_shape,nbox
   integer, optional                   :: iscalar
   integer, optional                   :: iarray(1:nbox)
   REALTYPE, optional                  :: scalar
   REALTYPE, optional                  :: array(1:nbox)
   REALTYPE, optional                  :: garray(1:nbox)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  Modifications: Marcello Vichi
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: iret,n=0
   integer                   :: idum(1:nbox)
   REAL_4B                   :: r4,dum(1:nbox)
!
!-----------------------------------------------------------------------
!BOC
   store_data = 0
   return
   end function store_data
!EOC


!-----------------------------------------------------------------------

   subroutine check_err(iret)
     integer iret
   end subroutine check_err
!-----------------------------------------------------------------------

   end module netcdf_bfm


!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

