
!--------------------------------------------------------------------------
!   
!    Copyright (C) 2005-2018  Petras Zemlys
!    Copyright (C) 2005-2018  Ali Erturk
!   
!    This file is part of SHYFEM. (m)
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

!module aquabc_II_vars
!module aquabc_II_layers_data
!not used module aquabc_II_wc_ini 
!module aquabc_II_sed_ini

!==================================================================
       module aquabc_II_vars
!==================================================================
!
! Initilizes all parameters from former include files
! and contains some commonly used variables
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!             Main parameters of Water quality model
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer, parameter :: nstate   =  30   ! Number of Water column state variables
      integer, parameter :: noutput  =  54   ! Number WC variables plus derived number of derived variables
      integer, parameter :: nconst   =  291  ! Number of WC constants
      integer, parameter :: nsstate  =  24   ! Number of Bottom Sediments state variables
      integer, parameter :: nsoutput =  26   ! Number of BS variables plus number of derived variables!
      integer, parameter :: nsconst  =  171  ! Number of BS constants
      integer, parameter :: noslay   =  7    ! Number of layers in bottom sediments

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Commonly used variables
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! Arrays for binary output of  state vars and derived vars
       real, save, allocatable :: wc_output(:,:,:)  ! for state vars  together with derived vars to be written to output files
                                                    ! derived vars are placed as nstate+1, ... elements
       real, save, allocatable :: sed_output(:,:,:) ! BS output
      
       double precision, save  :: da_out(4)         ! for controlling of binary output for WC
       
       double precision, save  :: da_outs(4)        ! for controlling of binary output for BS
       
       real, save              :: hlv_sed(noslay) ! For binary output
   
!==================================================================
! contains
!==================================================================


!==================================================================
      end module aquabc_II_vars
!==================================================================


module aquabc_II_layers_data

! Module for management saved outputs in pelagic model

    !use AQUABC_II_GLOBAL, only : DBL_PREC
    implicit none

    integer :: NUM_SHYFEM_LAYERS
    
    type SAVED_OUTPUT_MATRIX_DS
        double precision, dimension(:,:), allocatable :: MATRIX
    end type SAVED_OUTPUT_MATRIX_DS
    
   type(SAVED_OUTPUT_MATRIX_DS), dimension(:), allocatable :: SAVED_OUTPUTS_IN_LAYERS 
    
contains

    subroutine ALLOC_SAVED_OUTPUTS_IN_LAYERS &
               (NUM_SHYFEM_LAYERS, NUM_WC_SAVED_OUTPUTS, NUM_NODES_IN_LAYERS)
        
        ! INPUT ARGUMENTS
        integer, intent(in) :: NUM_SHYFEM_LAYERS     ! Number of layers in SHYFEM grid
        integer, intent(in) :: NUM_WC_SAVED_OUTPUTS  ! Number of saved outputs in water column
        
        ! Array that keeps the number of nodes for each layer
        integer, dimension(NUM_SHYFEM_LAYERS), intent(in) :: NUM_NODES_IN_LAYERS
        
        ! AUXILLARY VARIABLES
        integer :: SHYFEM_LAYER_NO, NUM_NODES
           
        allocate(SAVED_OUTPUTS_IN_LAYERS(NUM_SHYFEM_LAYERS))
        
        do SHYFEM_LAYER_NO = 1, NUM_SHYFEM_LAYERS
            NUM_NODES = NUM_NODES_IN_LAYERS(SHYFEM_LAYER_NO)
            allocate(SAVED_OUTPUTS_IN_LAYERS(SHYFEM_LAYER_NO) % MATRIX(NUM_NODES, NUM_WC_SAVED_OUTPUTS))
            SAVED_OUTPUTS_IN_LAYERS(SHYFEM_LAYER_NO) % MATRIX(:,:) = 0.0D0
        end do
        
    end subroutine ALLOC_SAVED_OUTPUTS_IN_LAYERS

    
    subroutine DEALLOC_SAVED_OUTPUTS_IN_LAYERS(NUM_SHYFEM_LAYERS)
        
        ! INPUT ARGUMENTS
        integer, intent(in) :: NUM_SHYFEM_LAYERS     ! Number of layers in SHYFEM grid

        ! AUXILLARY VARIABLES
        integer :: SHYFEM_LAYER_NO
        
        do SHYFEM_LAYER_NO = 1, NUM_SHYFEM_LAYERS
            deallocate(SAVED_OUTPUTS_IN_LAYERS(SHYFEM_LAYER_NO)%MATRIX)
        end do

        deallocate(SAVED_OUTPUTS_IN_LAYERS)        
    end subroutine DEALLOC_SAVED_OUTPUTS_IN_LAYERS
    
end module aquabc_II_layers_data




! !==================================================================
 ! module aquabc_II_wc_ini
! !==================================================================
! !
! ! Initilizes some variables necessary for WC calculations

   ! implicit none
   ! double precision, save :: frac_avail_DON
   
! !==================================================================
 ! contains
! !==================================================================
 ! subroutine calc_frac_avail_DON
 ! ! Calculates fraction of available DON for cyanobacteria
 ! ! it is called in pelagic model every step
 
  ! frac_avail_DON = 0.25  ! sedt
  ! !frac_avail_DON = 0.20  ! sedt not used
  ! !frac_avail_DON = 0.15  ! sedt1
  ! !frac_avail_DON = 0.10  ! sedt2
  ! !frac_avail_DON = 0.01  ! sedt3 

 ! end subroutine calc_frac_avail_DON
  

! !==================================================================
 ! end module aquabc_II_wc_ini
! !==================================================================


!==================================================================
 module aquabc_II_sed_ini
!==================================================================
!  Initialises some sediments properties and constants also
!  
!  Initialization is done using routine sed_properties_ini(bsedim)
!  called in aquabc_II_fem_interface initialisation part
!  Recalculates BS state variables
!



	  implicit none
	  

	   ! Variables used in BS model
	   ! It is assumed nkndi == nkn
	   
       ! Thickness of sediment layers	   
	   double precision, allocatable :: SED_DEPTHS      (:)       ! noslay,  
	   double precision, allocatable :: SED_DEPTHS_fast (:,:)     ! nkndi,noslay	   
	   
	   ! Variables that could get values from sedtrans model  if it is used
       double precision, allocatable  :: SED_POROSITIES_fast(:,:)        ! nkndi,noslay
       double precision, allocatable  :: SED_DENSITIES_fast (:,:)        ! nkndi,noslay
       double precision, allocatable  :: H_ERODEP_fast(:)                ! nkndi, eroded or deposited layer thickness
       double precision, allocatable  :: PART_MIXING_COEFFS_fast (:,:,:) ! nkndi,noslay, nsstate
       double precision, allocatable  :: SED_BURRIALS_fast(:,:)          ! nkndi,noslay
        
       
       ! Variables that values  are hardcoded for layers and used      
       !                  for initialisation of 2D variables       

       double precision :: SURF_MIXLEN
       double precision :: ADVECTIVE_VELOCITY
       
       ! Variables used always in WC and also BS model if it is run.  
       double precision, allocatable :: SETTLING_VELOCITIES(:)  !  nstate   
       
       double precision :: frac_avail_DON  
      
      
       real :: WATER_DENSITY     
      


!==================================================================
 contains
!==================================================================


 subroutine sed_properties_first(bsedim)
 
 use aquabc_II_vars
 use aquabc_pel_state_var_indexes
 use basin, only: nkndi
 
 ! Sediment properties initialization before the first time step
 implicit none
 
 logical bsedim
 
! ASSIGNMENT VALUES FOR HARDCODED WC VARIABLES      
! Settling velocities for WC variables

      if(nstate .ne. 30) then
       print *,'SED_PROPERTIES_FIRST:'       
       print *,'Number of statevars should be 30 but is', nstate
       stop
      end if
      
      if (.not. allocated(SETTLING_VELOCITIES)) allocate(SETTLING_VELOCITIES(nstate))   


      settling_velocities(NH4_N_INDEX         )  = 0.0D+0 !units: m/day ! AMMONIUM NITROGEN
      settling_velocities(NO3_N_INDEX         )  = 0.0D+0               ! NITRATE NITROGEN
      settling_velocities(PO4_P_INDEX         )  = 0.0D+0               ! ORTHOPHOSPHATE PHOSPHORUS
      settling_velocities(DISS_OXYGEN_INDEX   )  = 0.0D+0               ! DISSOLVED OXYGEN
      settling_velocities(DIA_C_INDEX         ) = 3.0D-1               ! DIATOMS CARBON
      settling_velocities(ZOO_C_INDEX         ) = 2.5D-1               ! ZOOPLANKTON CARBON
      settling_velocities(ZOO_N_INDEX         ) = 0.0D+0               ! ZOOPLANKTON NITROGEN
      settling_velocities(ZOO_P_INDEX         ) = 0.0D+0               ! ZOOPLANKTON PHOSPHORUS
      settling_velocities(DET_PART_ORG_C_INDEX) = 1.0D-1               ! DETRITUS PARTICULATE ORG. CARBON
      settling_velocities(DET_PART_ORG_N_INDEX) = 1.0D-1               ! DETRITUS PARTICULATE ORG. NITROGEN
      settling_velocities(DET_PART_ORG_P_INDEX) = 1.0D-1               ! DETRITUS PARTICULATE ORG. PHOSPHORUS
      settling_velocities(DISS_ORG_C_INDEX    ) = 0.0D+0               ! DISSOLVED ORGANIC CARBON
      settling_velocities(DISS_ORG_N_INDEX    ) = 0.0D+0               ! DISSOLVED ORGANIC NITROGEN
      settling_velocities(DISS_ORG_P_INDEX    ) = 0.0D+0               ! DISSOLVED ORGANIC PHOSPHORUS
      settling_velocities(CYN_C_INDEX         ) = 2.0D-1               ! CYANOBACTERIA CARBON
      settling_velocities(OPA_C_INDEX         ) = 2.0D-1               ! OTHER PHYTOPLANKTON CARBON
      settling_velocities(DISS_Si_INDEX       ) = 0.0D+0               ! DISSOLOVED SILICA
      settling_velocities(PART_Si_INDEX       ) = 2.0D-1               ! BIOGENIC SILICA
      settling_velocities(FIX_CYN_C_INDEX     ) = 2.0D-1               ! Nitrogen fixing cyano-bateria
      settling_velocities(INORG_C_INDEX       ) = 0.0D+0               ! Dissolved inorganic carbon
      settling_velocities(TOT_ALK_INDEX       ) = 0.0D+0               ! Alkalinity
      settling_velocities(FE_II_INDEX         ) = 0.5D+0               ! Iron2+
      settling_velocities(FE_III_INDEX        ) = 0.5D+0               ! Iron3+
      settling_velocities(MN_II_INDEX         ) = 0.5D+0               ! Manganese2+
      settling_velocities(MN_IV_INDEX         ) = 0.5D+0               ! Manganese4+    
      settling_velocities(CA_INDEX            ) = 0.5D+0               !Calcium total
      settling_velocities(MG_INDEX            ) = 0.5D+0               !Magnesium total
      settling_velocities(S_PLUS_6_INDEX      ) = 0.0D+0               !Sulphate sulphur total
      settling_velocities(S_MINUS_2_INDEX     ) = 0.0D+0               !Sulphide sulphur total
      settling_velocities(CH4_C_INDEX         ) = 0.0D+0               !Methane carbon total

      
      
      
      
      
! END ASSIGNMENT VALUES FOR HARDCODED WC VARIABLES

    if(bsedim) then 
    
      if (.not. allocated(SED_DEPTHS))      allocate(SED_DEPTHS(noslay))
      if (.not. allocated(SED_DEPTHS_fast)) allocate(SED_DEPTHS_fast (nkndi,noslay))

      ! Thickness of sediment layers in meters       
!       sed_depths(1) = 0.005D+0
!       sed_depths(2) = 0.005D+0
!       sed_depths(3) = 0.01D+0
!       sed_depths(4) = 0.01D+0
!       sed_depths(5) = 0.02D+0
!       sed_depths(6) = 0.02D+0
!       sed_depths(7) = 0.03D+0
      
      call read_sed_levels
      
      !  Assignment to the array for faster performance
      sed_depths_fast (:,1) = SED_DEPTHS(1) !0.005D+0
      sed_depths_fast (:,2) = SED_DEPTHS(2) !0.005D+0
      sed_depths_fast (:,3) = SED_DEPTHS(3) !0.005D+0
      sed_depths_fast (:,4) = SED_DEPTHS(4) !0.005D+0
      sed_depths_fast (:,5) = SED_DEPTHS(5) !0.03D+0
      sed_depths_fast (:,6) = SED_DEPTHS(6) !0.05D+0
      sed_depths_fast (:,7) = SED_DEPTHS(7) !0.05D+0 
       
   end if
 
 end subroutine sed_properties_first
 
!**************************************************
 
 subroutine read_sed_levels
 
! Reads BS levels from file that name is given in main *.str
! Calculates thickness of sediment layers
 
    use aquabc_II_vars
 
    !integer, parameter :: max_bs_lev  = 15
    implicit none
    
    real :: bs_lev (noslay)      ! BS levels read from file
    integer nlevf                ! number of levels obtained from file
    character*80 filen           ! file name
    integer ifile                ! number of symbols in file name 
    integer iunit                ! file unit
    integer ifileo               ! function to obtain file unit
    integer i
    
    call getfnm('bbs_lev',filen)
    !call trimline(filen,ifile)
    iunit = ifileo(0,trim(adjustl(filen)),'form','old')
    
    if( iunit .le. 0 ) then         
       write(6,*) 'READ_SED_LEVELS: Cannot open BS levels file'
       write(6,*) 'filename: ',trim(filen)
       stop
    end if
    
    write(6,*) 'READ_SED_LEVELS: Reading BS levels file '
    write(6,*) 'filename: ',trim(adjustl(filen))
    
    
    do nlevf = 1,noslay
     read(iunit,*, end = 200, err = 200) bs_lev(nlevf)     
!      if(nlevf .gt. noslay) then
!        write(6,*) 'READ_SED_LEVELS: Number of levels in file '
!        write(6,*) '        is greather than defined by NOSLAY'
!        write(6,*) 'NOSLAY: ', noslay,'nlevf: ', nlevf       
!        stop
!      end if 
    end do
    


!     if(nlevf .lt. noslay) then
!        write(6,*) 'READ_SED_LEVELS: Number of levels in file '
!        write(6,*) '        is less than defined by NOSLAY'
!        write(6,*) 'NOSLAY: ', noslay,'nlevf: ', nlevf        
!        stop
!      end if
     
     SED_DEPTHS(1) = bs_lev(1)
     i=1
     if(SED_DEPTHS(i) .le. 0.) then
      write(6,*) 'READ_SED_LEVELS: BS layer thicknes is negative or zero ' 
      write(6,*) 'SED_DEPTHS(i): ', SED_DEPTHS(i),'i: ', i        
      stop
     end if
     
     do i=2,noslay
       SED_DEPTHS(i) = bs_lev(i) - bs_lev(i-1)
       if(SED_DEPTHS(i) .le. 0.) then
        write(6,*) 'READ_SED_LEVELS: BS layer thicknes is negative or zero ' 
        write(6,*) 'SED_DEPTHS(i): ', SED_DEPTHS(i),'i: ', i        
        stop
       end if 
     
     end do
     
     write(6,*) 'READ_SED_LEVELS: BS levels file read succesfully'
      
     return
     
200  continue
     write(6,*) 'READ_SED_LEVELS: Error reading BS levels file '       
 
 end subroutine read_sed_levels
 
!**************************************************

 subroutine sed_properties_ini(bsedim)
 
!  Management of sediment properties initialisation every time step:
!
!  -Settling velocities  values are hardcoded.
!  -Sediment physical properties are obatained from sed. trans. model and prepared
!   to use in aquabc
!  -Thickness of sediment layers is hardcoded here
!  -SED_BURRIALS, PART_MIXING_COEFFS,SURF_MIXLEN, ADVECTIVE_VELOCITY, WATERDENSITY are hardcoded.
!
!   if bsedim == .false. only settling velocities and dissolved fractions are processed.
!   Called:
!        by aquabc_II_fem_interface
!  Note:
!   noslay is not equal to nlbdim in sediment transport

       use para_aqua
       use aquabc_II_vars  
       
       use levels, only: nlvdi,ilhkv  
       use basin,  only: nkndi,ipv
       
       use mod_sediment, only: nlbdim, tcn, tmsus, bdens, bleve
       ! Meaning:
       ! nlbdim            - maximum number of sediment layers
       ! tmsus(nkn)        - total suspended(eroded) sediment load,kg/s per reactor
       !                      > 0 - erosion, < 0 - deposition
       ! tcn  (nlvdi ,nkn) - total suspended(in water) concentration,kg/m3
       ! bleve(nlbdim,nkn) - bed levels(depths below sediment surface), m
       ! bdens(nlbdim,nkn) - bed dry bulk density of bottom sediments, kg/m^3
 
       
       implicit none
       
       !include 'aquabc_II.h'
       !include 'param.h'
       include 'femtime.h'
       !include 'levels.h'
       !include 'basin.h'   !ipv is used from there
       
       !integer ipv(nkndi)	!external node numbers
       !common  /ipv/ipv

      logical bsedim
      
      real getpar !function for obtaining parameters values defined in *.str file
      
      !integer, parameter :: nlbdim  = 20 ! maximum number of sediment levels in sedtrans model
                                         ! It should be the same as in sed_param.h. Temporary solution 
                      ! noslay - number of layers in aquabc from include
      !integer nlbdim1 ! maximum number of sediment levels past by sedtrans model(for errors control)              
      integer nlb     ! actual number of sediment layers in sedtrans
      
      integer i, k, j, j_int
      
      integer isedi   
                   
      double precision  sed_levels(noslay+1)  ! levels for aquabc BS
      double precision  sed_centers(noslay)   ! sedimen layer centers for aquabc model
      
      real bed_dens_int
      
      real depth,volold,area !variable to get reactor area
      
      !Variables to store all nodes obtained from sedtrans model       
      !real stress
      !real t_total_load(nkn)
      !real t_total_susp_conc(nlvdi,nkn)
      !real t_bed_lev (nlbdim,nkn) 
      !real t_bed_dens(nlbdim,nkn)
      
      !Variables for one node obtained from sedtrans model       
      !real stress
      real total_load
      real total_susp_conc(nlvdi)
      real bed_lev(nlbdim) ! sediment transport model has 
      real bed_dens(nlbdim)! the last layer with not defined thickness
                           ! It is just layer below last level.
                           
      real bed_lev_c(nlbdim)! bed layers centers (calculated localy)                    
                            ! number of elements is nlb             
                           
      double precision :: SED_DENSITIES      (noslay) 
      double precision :: SED_POROSITIES     (noslay)       
      double precision :: SED_BURRIALS       (noslay)
      double precision :: PART_MIXING_COEFFS (noslay)
      real :: SED_ERODEP
      real :: H_ERODEP !Auxilary

!     Crosschecking of parameters:
      if(nstate .ne. 30) then
       print *, 'SED_PROPERTIES_INI:' 
       print *, 'Number of WC variables is not equal 30 but',nstate
       print *,  'Change its value in module SED_PROPERTIES_INI' 
       stop 
      end if
      
      if(nsstate .ne. 24) then
       print *, 'SED_PROPERTIES_INI:' 
       print *, 'Number of BS variables is not equal 24 but',nstate
       print *,  'Change its value in module SED_PROPERTIES_INI' 
       stop 
      end if
      
      
      if (nkndi .ne. 1309) then
       print *, 'SED_PROPERTIES_INI:' 
       print *, 'Number of nodes is not equal 1309 but', nkndi
       print *, 'Change its value in module SED_PROPERTIES_INI'
       stop 
      end if
      

      
      
      if(.not. bsedim) return
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! PREPARATION FOR BS SIMULATION

      if (noslay .ne. 7) then
       print *, 'SED_PROPERTIES_INI:' 
       print *, 'Number of BS layers is not equal 7 but', noslay
       print *, 'Change its value in module SED_PROPERTIES_INI'
       stop 
      end if
      
      if(.not. allocated(SED_POROSITIES_fast))     allocate(SED_POROSITIES_fast(nkndi,noslay))            
      if(.not. allocated(SED_DENSITIES_fast))      allocate(SED_DENSITIES_fast (nkndi,noslay))            
      if(.not. allocated(H_ERODEP_fast))           allocate(H_ERODEP_fast(nkndi)) 
      if(.not. allocated(PART_MIXING_COEFFS_fast)) allocate(PART_MIXING_COEFFS_fast (nkndi,noslay, nsstate)) 
      if(.not. allocated(SED_BURRIALS_fast))       allocate(SED_BURRIALS_fast(nkndi,noslay)) 
      

      isedi = nint(getpar('isedi'))
      
      WATER_DENSITY = 1.  ! should be checked for units in wet density calculations                         

      
     
     if (isedi .ne. 0) then
     
! OBTAINING INFORMATION FROM SEDTRANS MODEL      
      sed_levels(1) = 0.       
      do i = 2, noslay + 1
       sed_levels(i)    = sed_levels(i-1) + SED_DEPTHS(i-1)
       sed_centers(i-1) = (sed_levels(i-1) + sed_levels(i))/2.
      end do
       

      
      
      
      
          
      do k=1,nkndi
      
!--------------------------------------------------------------------------      
!     This subroutine is now called inside sedtrans model obtaining values for all nodes
!           
       !call get_sedim_prop(k, nlbdim1,stress, total_load, &
                           !total_susp_conc, bed_lev, bed_dens)
!--------------------------------------------------------------------------
! Assigning values obtained from STM module mod_sediment to local variables 
      total_load               = tmsus(k)                          
      total_susp_conc(1:nlvdi) = tcn  (1:nlvdi, k) 
      bed_lev (1:nlbdim)       = bleve(1:nlbdim,k)
      bed_dens(1:nlbdim)       = bdens(1:nlbdim,k)

      ! Meaning of variables:      
      !nlbdim                                      - maximum number of levels      
      !total_load            - tmsus(nkndi)        - total suspended(eroded) sediment load,kg/s per reactor
      !                      -                      > 0 - erosion, < 0 - deposition
      !total_susp_conc(nlvdi)- tcn(nlvdi,nkndi)    - total suspended(in water) concentration,kg/m3
      !bed_lev(nlbdim)       - bleve(nlbdim,nkndi) - bed levels(depths below sediment surface), m
      !bed_dens(nlbdim)      - bdens(nlbdim,nkndi) - bed dry bulk density of bottom sediments, kg/m^3
      
      ! porosity can be calculated from bulk density: fi = 1 - (ro_bulk/ro_particle)
      ! ro_particle can be assumed approximatelly equal 2.65 g/cm3 = 2650 kg/m3
      !----------------------------------------------------------------------
 
              
       ! Erosion, deposition per current time step.
       
       !eroded/deposited mass per time step:
       call dvanode(1,k,-1,depth,volold,area)        !gets old depth, volume and area
                                                     ! Is the area the same fo all layers? fixme             
       SED_ERODEP = total_load*dt_act/area
       
       !eroded/deposited thickness per time step:              
       H_ERODEP_fast(k)   = SED_ERODEP/bed_dens(1) !> 0 - erosion, < 0 - deposition
       H_ERODEP = H_ERODEP_fast(k)
       
       ! Case when erosion is higher than first layer, fixme      
       if (abs(H_ERODEP) .gt. SED_DEPTHS_fast (k,1) .and. H_ERODEP .gt. 0.) then
       
         
         H_ERODEP_fast(k) = SED_DEPTHS_fast (k,1)

         print *, ''
         print *, 'eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee'       
         print *, 'SED_PROPERTIES_INI: Erosion thickness is higher than thicknes of 1-st layer'
         print *, 'Only erosion of whole 1-st layer is assummed' 
         print *, 'Internal node No:', k
         print *, 'External node No:', ipv(k)
         print *, 'H_ERODEP:', H_ERODEP
         print *, 'SED_DEPTHS_fast(k,1):', SED_DEPTHS_fast(k,1)
         print *, 'total_load ', total_load 
         print *, 'dt_act     ', dt_act     
         print *, 'SED_ERODEP ', SED_ERODEP 
         print *, 'bed_dens(1)', bed_dens(1)
         print *, '*********All sediment properties returned by GET_SEDIM_PROP:'
         !print *,'bottom stress=',stress
         print *,'total_bed_load=',total_load
         print *,'total_susp_conc=',total_susp_conc(ilhkv(k))       
         print *,'bed_lev=',bed_lev
         print *,'bed_dens=',bed_dens
         !print *, 'H_ERODEP_fast=', H_ERODEP_fast
         print *, 'eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee'
         print *, ''
         !stop
       
       end if
       

       
       !number of active levels in sedtrans 
       !(number of levels not equal to zero + first level)
       nlb  = count(bed_lev .gt. 0) + 1 
       j_int = 0
             
       
       ! layer centers in sedtrans
       bed_lev_c(:) = 0.
       do i=1,nlb-1       
        bed_lev_c(i) = (bed_lev(i) + bed_lev(i+1))/2.         
       end do
       
       ! interpolation        
       do i=1,noslay
          if(sed_centers(i) .gt. bed_lev_c(nlb-1)) then
             bed_dens_int = bed_dens(nlb-1)
          else if(sed_centers(i) .lt. bed_lev_c(1)) then
             bed_dens_int = bed_dens(1)
          else
            do j=1,nlb-1
              if(sed_centers(i) .ge. bed_lev_c(j) .and. &
                 sed_centers(i) .le. bed_lev_c(j+1))  j_int = j  
            end do
            bed_dens_int = bed_dens(j_int) + (bed_dens(j_int+1) - &  
                           bed_dens(j_int))/(bed_lev_c(j_int+1) - bed_lev_c(j_int)) * &
                           (sed_centers(i) - bed_lev_c(j_int))
          end if
          
          SED_POROSITIES_fast(k,i) = 1. - (bed_dens_int/2650.)
          ! from wet to dry bulk density:  kg/m3 to g/cm3
          ! Wet bulk density. Water density now is 1.
          SED_DENSITIES_fast (k,i) = 0.001*bed_dens_int + WATER_DENSITY * SED_POROSITIES_fast(k,i)
          
       end do

     if (ipv(k)	.eq. 0) then  !198 2265 (vidmares)
                         
       print *, '*********Sediment properties returned by GET_SEDIM_PROP for node 2265:******'
       print *,'node=', k 
       !print *,'stress=',stress
       print *,'total_bed_load=',total_load
       print *,'total_susp_conc=',total_susp_conc(ilhkv(k))       
       print *,'bed_lev=',bed_lev
       print *,'bed_dens=',bed_dens
       
       print *, ' '
       print *, '*********Sediment properties calculated for AQUABC, node 2265:******' 
       print *, 'H_ERODEP:', H_ERODEP_fast(k)
       print *, 'SED_LAYER_THICKNESS:', SED_DEPTHS_fast(k,:)
       print *, 'time step:     ', dt_act     
       print *, 'Bed load per time step: ', SED_ERODEP 
       print *, 'SED_POROSITIES=', SED_POROSITIES_fast(k,:)
       print *, 'SED_WET_DENSITIES', SED_DENSITIES_fast (k,:)     
       
      end if       
  
        
      end do !k=1,nkndi
! END OBTAINING INFORMATION FROM SEDTRANS MODEL



      
! CASE WHEN SEDTRANS NOT USED(sedi==0)    
     else  
      H_ERODEP_fast(:) = 0.
     end if	
      
     !write(1000, *) H_ERODEP_fast(2265)      

      

! INITIALISATION OF SED PRPPERTIES IN CASE WHEN SEDTRANS IS NOT USED
    if(isedi .eq. 0) then 
    
       if (noslay .ne. 7) then                                      
        print *, 'SED_PROPERTIES_INI:'                              
        print *, 'Number of BS layers is not equal 7 but', noslay   
        print *, 'Change its value and assignments to variables in module SED_PROPERTIES_INI'           
        stop                                                                          
      end if

!     Porosities of BS layers. For the case if sediment transport model is not used
!     Local       
      SED_POROSITIES(1) = 0.40
      SED_POROSITIES(2) = 0.40
      SED_POROSITIES(3) = 0.40
      SED_POROSITIES(4) = 0.40
      SED_POROSITIES(5) = 0.30
      SED_POROSITIES(6) = 0.25
      SED_POROSITIES(7) = 0.25
      
!     Assignment to the array for faster performance
!     Global 
      SED_POROSITIES_fast(:,1) = SED_POROSITIES(1) !0.40D+0
      SED_POROSITIES_fast(:,2) = SED_POROSITIES(2) !0.40D+0
      SED_POROSITIES_fast(:,3) = SED_POROSITIES(3) !0.40D+0
      SED_POROSITIES_fast(:,4) = SED_POROSITIES(4) !0.40D+0
      SED_POROSITIES_fast(:,5) = SED_POROSITIES(5) !0.30D+0
      SED_POROSITIES_fast(:,6) = SED_POROSITIES(6) !0.25D+0
      SED_POROSITIES_fast(:,7) = SED_POROSITIES(7) !0.25D+0
        
!     Wet sed. densities for BS layers, g/cm3 = kg/L,
!     For the case if sediment model is not used
!     Local     
      SED_DENSITIES(1) = 1.75
      SED_DENSITIES(2) = 1.75
      SED_DENSITIES(3) = 1.75
      SED_DENSITIES(4) = 1.75
      SED_DENSITIES(5) = 1.75
      SED_DENSITIES(6) = 1.75
      SED_DENSITIES(7) = 1.75
      
!     Densities for BS layers
!     Global
      SED_DENSITIES_fast(:,1) = SED_DENSITIES(1) ! 1.75D+0
      SED_DENSITIES_fast(:,2) = SED_DENSITIES(2) ! 1.75D+0
      SED_DENSITIES_fast(:,3) = SED_DENSITIES(3) ! 1.75D+0
      SED_DENSITIES_fast(:,4) = SED_DENSITIES(4) ! 1.75D+0
      SED_DENSITIES_fast(:,5) = SED_DENSITIES(5) ! 1.75D+0
      SED_DENSITIES_fast(:,6) = SED_DENSITIES(6) ! 1.75D+0
      SED_DENSITIES_fast(:,7) = SED_DENSITIES(7) ! 1.75D+0 
      
    end if       
! END INITIALISATION OF SED PRPPERTIES IN CASE WHEN SEDTRANS IS NOT USED         

! INITIALISATION OF PROCESS RATE PARAMETERS FOR BS MODEL.
! In future should be read from file. Fixme
! Only 'fast' variables are global    
      
!     Burial rate of BS state variables (around 1cm/year) !m/day      
      SED_BURRIALS(1) =  2.7397D-5
      SED_BURRIALS(2) =  2.7397D-5
      SED_BURRIALS(3) =  2.7397D-5
      SED_BURRIALS(4) =  2.7397D-5
      SED_BURRIALS(5) =  2.7397D-5
      SED_BURRIALS(6) =  2.7397D-5
      SED_BURRIALS(7) =  2.7397D-5
      
!     Burial rate of BS state variables (around 1cm/year) !m/day
      SED_BURRIALS_fast(:,1) = SED_BURRIALS(1) !2.7397D-5 !2.7397D-5    !m/day
      SED_BURRIALS_fast(:,2) = SED_BURRIALS(2) !2.7397D-5 !2.7397D-5
      SED_BURRIALS_fast(:,3) = SED_BURRIALS(3) !2.7397D-5 !2.7397D-5
      SED_BURRIALS_fast(:,4) = SED_BURRIALS(4) !2.7397D-5
      SED_BURRIALS_fast(:,5) = SED_BURRIALS(5) !2.7397D-5
      SED_BURRIALS_fast(:,6) = SED_BURRIALS(6) !2.7397D-5
      SED_BURRIALS_fast(:,7) = SED_BURRIALS(7) !2.7397D-5

      
!     Particle mixing diffusion coeff.  
!     Made zero for solutes in sediment routine!     
      PART_MIXING_COEFFS(1) = 2.64D-5
      PART_MIXING_COEFFS(2) = 2.64D-5
      PART_MIXING_COEFFS(3) = 2.64D-5
      PART_MIXING_COEFFS(4) = 2.64D-5
      PART_MIXING_COEFFS(5) = 2.64D-5
      PART_MIXING_COEFFS(6) = 2.64D-5
      PART_MIXING_COEFFS(7) = 2.64D-5
      
!     Particle mixing diffusion coeff.
!     Made zero for solutes in sediment routine! I am not sure, additional check is necessary. fixme
      PART_MIXING_COEFFS_fast (:,1, :) = PART_MIXING_COEFFS(1) ! 2.64D-5 ! 0.d+0
      PART_MIXING_COEFFS_fast (:,2, :) = PART_MIXING_COEFFS(2) ! 2.64D-5 ! 0.d+0
      PART_MIXING_COEFFS_fast (:,3, :) = PART_MIXING_COEFFS(3) ! 2.64D-5 ! 0.d+0
      PART_MIXING_COEFFS_fast (:,4, :) = PART_MIXING_COEFFS(4) ! 2.64D-5
      PART_MIXING_COEFFS_fast (:,5, :) = PART_MIXING_COEFFS(5) ! 2.64D-5
      PART_MIXING_COEFFS_fast (:,6, :) = PART_MIXING_COEFFS(6) ! 2.64D-5
      PART_MIXING_COEFFS_fast (:,7, :) = PART_MIXING_COEFFS(7) ! 2.64D-5

      
!     Surface mixing length for exchange with WC     
      SURF_MIXLEN = 0.10
      
!     advective velocity in BS  
      ADVECTIVE_VELOCITY = 0
!     
! END INITIALISATION OF PROCESS RATE PARAMETERS FOR BS MODEL.     
      
 end subroutine sed_properties_ini 
 
 
      
!**********************************************      
!**********************************************
 subroutine sed_recalc_ini(es)
  
! Recalculates initial NH4 and PO4 as  concentrations in porewater
! to concentrations in total sediments (pore water + solids)
       
        use para_aqua
        use aquabc_II_vars 
        use basin
        
        implicit none
        
        !include 'aquabc_II.h'
        !include 'nbasin.h'
        
        real, intent(inout) :: es        (noslay,nkn,nsstate)  
!        real, intent(inout) :: sed_output(noslay,nkn,nsoutput)
        real SOLUTE_FRACTIONS_NH4N(nkn,NOSLAY)
        real SOLUTE_FRACTIONS_PO4P(nkn,NOSLAY)
        real SOLID_CONCS(nkn,noslay)
        
        double precision SOLID_PART_COEFF_NH4
        double precision SOLID_PART_COEFF_PO4
        
        integer i,j,k



          call para_get_value('SOLID_PART_COEFF_NH4',SOLID_PART_COEFF_NH4)  !41
          call para_get_value('SOLID_PART_COEFF_PO4',SOLID_PART_COEFF_PO4)  !42
          

            ! dry bulk density
            SOLID_CONCS(1:nkn,1:noslay) = (SED_DENSITIES_fast(1:nkn,1:noslay) - &
                          (WATER_DENSITY * SED_POROSITIES_fast(1:nkn,1:noslay)))
                          
            ! solute fractions of conentrations per total sediment volume
            SOLUTE_FRACTIONS_NH4N(1:nkn,1:noslay) =  &              
              (1.0 / (1.0 + (SOLID_CONCS(1:nkn,1:noslay) * SOLID_PART_COEFF_NH4)))

            SOLUTE_FRACTIONS_PO4P(1:nkn,1:noslay) =  &              
              (1.0 / (1.0 + (SOLID_CONCS(1:nkn,1:noslay) * SOLID_PART_COEFF_PO4)))

            sed_output(1:noslay,1:nkn,nsstate+1) = es(1:noslay,1:nkn,1) ! saves conc. per pore water for output
            sed_output(1:noslay,1:nkn,nsstate+2) = es(1:noslay,1:nkn,5) ! saves conc. per pore water for output
            
          do k=1,nkn
           do i=1,noslay
           
            ! recalculates NH4 to concentrations per total sed. vol.:
            es(i,k,1) = SED_POROSITIES_fast(k,i)* es(i,k,1)/ &
                                                    SOLUTE_FRACTIONS_NH4N(k,i)           
            
            ! recalculates PO4 to concentrations per total sed. vol.:           
            es(i,k,5) = SED_POROSITIES_fast(k,i) * es(i,k,5)/ &
                                                    SOLUTE_FRACTIONS_PO4P(k,i)

           end do !noslay
          end do  !nkn
          
          ! saves state variables to the output 
          sed_output(1:noslay,1:nkn,1:nsstate) = es(1:noslay,1:nkn,1:nsstate)
            
 end subroutine sed_recalc_ini

!************************************************

 subroutine sed_recalc_final(es)
 
! Recalculates NH4 and PO4 from total concentrations insediments (in pore water + solids)
! to solute concentrations in porewater
 
        use para_aqua
        use aquabc_II_vars 
        use basin
                  
        implicit none
        
        !include 'aquabc_II.h'
        !include 'nbasin.h'
        
        real, intent(inout) :: es(noslay,nkn,nsstate)  
!        real sed_output(NOSLAY,nkn,nsoutput)
      
        
        double precision SOLUTE_FRACTION_NH4N(nkn,noslay)
        double precision SOLUTE_FRACTION_PO4P(nkn,noslay)
        double precision SOLID_CONCS(nkn,noslay)
        
        double precision SOLID_PART_COEFF_NH4
        double precision SOLID_PART_COEFF_PO4
        
        integer i,j,k


           call para_get_value('SOLID_PART_COEFF_NH4',SOLID_PART_COEFF_NH4)  !41
           call para_get_value('SOLID_PART_COEFF_PO4',SOLID_PART_COEFF_PO4)  !4
!          

                !dry bulk density
                SOLID_CONCS(1:nkn,1:noslay) = (SED_DENSITIES_fast(1:nkn,1:noslay) - & 
                                        (WATER_DENSITY * SED_POROSITIES_fast(1:nkn,1:noslay)))
                                        
                ! Solute fractions of concentrations per total sediment volume
                SOLUTE_FRACTION_NH4N(1:nkn,1:noslay) =      &
                      (1.0 / (1.0 + (SOLID_CONCS(1:nkn,1:noslay) * SOLID_PART_COEFF_NH4)))
                
                SOLUTE_FRACTION_PO4P(1:nkn,1:noslay) =      &                     
                      (1.0 / (1.0 + (SOLID_CONCS(1:nkn,1:noslay) * SOLID_PART_COEFF_PO4)))
                

            ! recalculates total(solute + solids) conc. per total volume to conc. of solute per pore water
            do k = 1, nkn
               do j = 1, noslay
                
                es(j,k,1) = es(j,k,1) * &
                       SOLUTE_FRACTION_NH4N(k,j)/ SED_POROSITIES_fast(k,j)
                
                es(j,k,5) = es(j,k,5) *  &
                        SOLUTE_FRACTION_PO4P(k,j)/ SED_POROSITIES_fast(k,j)
               end do
            end do
            
 end subroutine sed_recalc_final

!==================================================================
 end module aquabc_II_sed_ini
!==================================================================
