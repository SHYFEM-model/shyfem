
!--------------------------------------------------------------------------
!
!    Copyright (C) 2005-2018  Petras Zemlys
!    Copyright (C) 2005-2018  Ali Erturk   
!    Copyright (C) 2005-2018  Georg Umgiesser
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

! CONTENT:
!  module AQUABC_II_GLOBAL
!  subroutine AQUABC_II      - controls calls to pelagic and BS routines
!  subroutine cur_euler      - calculates state variables from rates using Euler
!  FUNCTION STRANGERS(VALUE) - checks for NaNs, Ifs and strange values
!  subroutine layer_light    - calculates the intensity of light after penetration of current layer
!  subroutine light_bottom   - is used by subroutine layer_light
!  subroutine derived_vars   - calculates derived variables
!
!  subroutine get_nodal_area_code(k,ncode)(placed into subn35.f) - routine to get sediment type areas
!********************************************************************
!********************************************************************

module AQUABC_II_GLOBAL
    implicit none
    integer :: DBL_PREC
    parameter(DBL_PREC = selected_real_kind(15, 307))
end module AQUABC_II_GLOBAL

!*****************************************************
!*****************************************************


subroutine AQUABC_II(nkn,lmax_fast, nlayers, &
           t, dt,                            &
           vol_fast,vol_old_fast,area_fast,  &
           depth_fast,vel_fast,              &
           wtemp_fast, sal_fast,             &
           wind_fast,atemp_fast,             &
           ITOT_fast,FDAY,                   &
           ice_cover_fast,                   &
           loads,                            &
           par,nconst,                       &
           e_fast0,e_fast1,nstate,           &
           wc_output, noutput,               &
           dg, NDIAGVAR,                     &
           bsedim,                           &
           par_sed,nsconst,                  &
           INIT_SED_STATE_VARS_fast,  FINAL_SED_STATE_VARS_fast,  &
           nsstate,noslay,                   &
           waves_fast,                       &
           output_sed_full, nsoutput,        &
           INTER_sed_full,  NDIAGVAR_sed,ifirst_it)


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! Aquatic biochemical cycling model
    ! This program serves as second level interface between SHYFEM
    ! and WC, BS kinetic models
    !
    ! Consists of WC pelagic and sediment module (on development)
    ! Calls BS module when all WC layers are processed.
    ! Arguments:
    
    ! Water column:
    !   nkn,lmax_fast, nlayers                    - number of nodes, max number of layers for nodes
    !   t,dt                                      -
    !   vol_fast,vol_old_fast,area_fast           -
    !   depth_fast,vel_fast                       -
    !   wtemp_fast, sal_fast                      -
    !   ITOT_fast,FDAY                            -
    !   loads                                     -
    !   par - WC constants file                   - not used more. Fixme
    !   nconst - number of WC constants           -
    !   e_fast0,e_fast1,nstate                    -
    !   wc_output, noutput                        -
    !   INTER, NDIAGVAR
    !                             -
    ! Bottom sediments:
    !   bsedim - logical (true if BS is processed -
    !   par_sed - BS constants                    - not used more. fixme
    !   nsconst - number of BS constants          -
    !   INIT_SED_STATE_VARS_fast                  -
    !   FINAL_SED_STATE_VARS_fast                 -
    !   INTER_sed_full, NDIAGVAR_sed              -
    !
    !   ifirst_it !first time iteration indicator
    !
    !
    ! Notes:
    !    Incomming time and step (t,dt) should be in days
    !    Program returns reals in single precision and passes to kinetic routinas in
    !    double precision
    !    All variables with suffix 'fast' and 'full' are for speedup of calculations by
    !    vectorization:
    !                   full - variable is as it is used in fem interface
    !                   fast - elements of arrays are arranged according
    !                          dimensions: 1 - nodes, 2 - state variables, 3 - layers
    !
    !    Sediment types:
    !       0 = fine sand (blue);
    !       1 = silty fine sand (green)
    !       2 = coarse silt(sky-blue)
    !       3 = clayey coarse silt(red)
    !       4 = silty clayey mud (magenta)
    !       5 = clay(braun)
    !       6 = klaipeda strait(yelow)
    !       7 = Baltic Sea(dark grey)
    !      10 = Nemunas delta(light grey)
    !
    
    !
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use aquabc_II_sed_ini
    use para_aqua       
             
    use levels, only: nlvdi
    use basin, only: ipv, nkndi
    use aquabc_pel_state_var_indexes
    use aquabc_II_layers_data

    implicit none

    ! This include is commented because it is present in aquabc_II_sed_ini
    !include 'param.h'
    !include 'basin.h' !ipv is used from there

    !integer ipv(nkndim)	!external node numbers
    !common  /ipv/ipv

    logical bsedim                      !  false - no bottom sediments processing
                                        !  true - b. sed. procesing
    integer nkn                         !  total number of nodes (reactors)
    integer nlayers                     !  total number of levels in grid
    integer nnode                       !  internal node number
    integer exn_node                    !  external node number

    integer layer                       !  processed WC layer number
    integer LAYER_max                   !  maximum number of WC layers for the node(reactor)
    integer, allocatable :: node_active_num(:) !  internal numbers of active nodes for the current layer
    integer nactive_nodes                ! number of active nodes in current layer
    integer nactive_nodes_arr(nlayers)   ! array of number of active nodes in  layers
    real state(nstate)                   !  WC state variable [mg/L] == [g/m**3]
    !real output_wc(noutput)             !  WC state variables and derived variables for one node

    ! Arrays for faster calculations by vectorization for all nodes(reactors)
    real e_fast0(nkn,nstate,nlayers)     ! initial state vector arranged                                              !

    real e_fast1(nkn,nstate,nlayers)     ! final state vector arranged

    double precision, allocatable :: STATE_fast(:,:)  ! state variables array for one layer for call to pelagic kinetics
    real     states_fast(nkn,nstate)       ! state variables array for one layer
    integer  lmax_fast  (nkn)              ! max lev. numbers
    real    depth_fast  (nkn,nlayers)      ! depth array
    real   depth_layer  (nkn)              ! depth array for layer
    real      vel_fast  (nkn,nlayers)      ! current velocity
    real    wtemp_fast  (nkn,nlayers)      ! water temperature
    real      sal_fast  (nkn,nlayers)      ! water salinity
    real     wind_fast  (nkn)              ! wind speed
    real    atemp_fast  (nkn)              ! air temperature
    real  ice_cover_fast(nkn)              ! ice cover - area fraction of reactor
    real    light_fast  (nkn)              ! incident light for the current layer
    real    ITOT_fast   (nkn)              ! incident light for the first layer
    real    light_next  (nkn)              ! incident light for the next layer

    real      vol_fast  (nkn,nlayers)      ! reactor volume
    real  vol_old_fast  (nkn,nlayers)      ! reactor volume before one step
    real     area_fast  (nkn,nlayers)      ! reactor area
    real    waves_fast  (nkn,2)            ! wind waves:  1-significant  wave height
                                             !              2 - mean wave period
    real rmass_fast(nkn)                   ! resuspension mass, units?
    real ph_full(nkn, nlayers)


    real             wc_output (nlayers,nkn,noutput)    ! array for outputs (state vars + auxilaries)
    DOUBLE PRECISION, allocatable :: PH(:)              ! the same double precision

	real             dg  (nlayers,nkn,nstate,NDIAGVAR)  ! processes values for diagnostics
	DOUBLE PRECISION,allocatable :: DGS (:,:,:)           ! processes values for diagnostics

	real             rates_fast(nkn,nstate,nlayers)        ! derivatives for state variables
	DOUBLE PRECISION,allocatable :: DERIVATIVES_fast(:,:) ! derivatives for state variables
	DOUBLE PRECISION,allocatable :: SAVED_OUTPUTS_fast(:,:) ! saved outputs for one layer

    ! sediments
	integer SEDIMENT_TYPE_fast(nkn)
    real FRACTION_OF_DEPOSITION_fast(nkn,nstate)

	integer           n_driving_functions
    parameter        (n_driving_functions = 10)
    double precision, allocatable :: DRIVING_FUNCTIONS(:,:) ! for pelagic model

    real fluxes_from_upper_fast    (nkn,nstate)
    real fluxes_to_next_fast       (nkn,nstate) ! for upper layers
    !save fluxes_to_next_fast
    real not_deposited_fluxes_fast (nkn,nstate)
    real settling_rates_fast       (nkn,nstate)    ! to the bottom from the last layer

    double precision INIT_SED_STATE_VARS_fast  (nkn,noslay,nsstate) ! initial state array for processing
    double precision FINAL_SED_STATE_VARS_fast (nkn,noslay,nsstate) ! state array after one time step
    real             output_sed_full           (noslay,nkn,nsoutput) ! derived variables as is fem interface
    real             INTER_sed_full            (noslay,nkn,nsstate,NDIAGVAR_sed)! process rates or any intermediates as is in fem interface

    double precision FLUXES_TO_SEDIMENTS_fast(nkn,nsstate)

    ! The following variables moved to subroutine sed_properties_ini in aquabc_II_ased_ini.f90
    !
    !SED_DEPTHS_fast    (nkn,noslay)
    !SED_POROSITIES_fast(nkn,noslay)
    !SED_DENSITIES_fast (nkn,noslay)
    !SED_DEPTHS (NOSLAY)                ! Thicknes of BS layers, m
    !SED_POROSITIES(NOSLAY)             ! Sediment porosities
    !SED_DENSITIES(NOSLAY)              ! Sed. density , kg/m3
    !SURF_MIXLEN                        ! Surface mixing length for exchange with WC, m
    !SED_BURRIALS(NOSLAY)               ! Burial rates (advection) for each of BS layers, m/day
    !SETTLING_VELOCITIES(nstate)        ! m/day
    double precision dissolved_fractions_full(nkn,nstate,nlayers) ! Dissolved fractions for WC variables. For some vars calculated in pelagic model
    double precision DISSOLVED_FRACTIONS(nstate)                  ! The same for current node and  layer
    !ADVECTIVE_VELOCITY                 ! For dissolved
    !PART_MIXING_COEFFS_fast (nkn,noslay, nsstate)
    !SED_BURRIALS_fast(nkn,noslay)
    !

    double precision SED_DIFFUSIONS_fast(nkn,noslay,nsstate)
    double precision SED_TEMPS_fast     (nkn,noslay)
    double precision SURF_WATER_CONCS_fast(nkn,nsstate)

    double precision FLUXES_FROM_SEDIMENTS_fast(nkn,nsstate)
    real flux  !auxilary variable
    double precision PROCESSES_sed_fast(nkn,noslay,nsstate,ndiagvar_sed)
    double precision SED_OUTPUTS_fast  (nkn,noslay,nsoutput)




    real FLUXES_TO_ALUKAS_fast(nkn,nstate)

    ! End of arrays for faster calculations by vectorization for all nodes(reactors)

    real vol,volold               !volume [m**3]
    real area, depth              !depth of box [m]
    real VEL                      !velocity [m/s]
    real TEMP                     !water temperature [C]
    real WIND,AIRTMP
    real sal                      !salinity [psu] == [per mille]
    !real pH
    real ITOT,FDAY		          !light intensity ly/day, photoperiod
                                  !(fraction of day = 1, compatibility with light lim algorithm)
    double precision KE           ! Background extinction coefficient, from parameters


    integer nstate                      !  number of state variables
    integer noutput                     ! number of state variables + number of derived variables for the output
    integer nconst                      !  number of parameters not used fixxme
    integer i,j,k,kk,l
    real t 		                  !current time [day]
    real dt                       !time step [day]


    real cold(nstate)             !old state variable (for diagnostic purpose)
    real loads(nstate)	          !loading for c [g/(m**3 day)]. Array should be 3d. fixme
    real rates(nstate)            !WC state variables derivatives [g/m**3/day]

    real par(nconst)      !(800) ! WC model constants

    integer NDIAGVAR
    real INTER(nstate,NDIAGVAR)

    double precision STATE_VARIABLES(nstate)
    double precision DERIVATIVES(nstate)   !derivatives  [g/m3/day]
    double precision MODEL_CONSTANTS(nconst)


    double precision PROCESSES(nstate,NDIAGVAR)

    integer           n_saved_outputs
    parameter        (n_saved_outputs = 5)
    double precision SAVED_OUTPUTS(nkndi,n_saved_outputs,nlvdi )  !for storage of WC saved outputs 
                                                                  !with active and inactive nodes in layers
    integer    num_sed_saved_outputs                                                             
    parameter (num_sed_saved_outputs = 5)                                                              
    double precision, allocatable :: SED_SAVED_OUTPUTS(:,:,:) 
                    !(nkn,noslay, num_sed_saved_outputs)
    save SED_SAVED_OUTPUTS                                                            

    double precision SURFL
    double precision IA

    !integer    iprod
    !parameter (iprod=nkn*nlayers)
    !integer WCKIN_CALLED_BEFORE(nkn, nlayers) ! To save indicators if node
    !                                            ! and layer are already processed. Not used yet
    !data    WCKIN_CALLED_BEFORE /iprod*0/
    !save    WCKIN_CALLED_BEFORE
    integer       CALLED_BEFORE   ! Value of  WCKIN_CALLED_BEFORE
    integer   nflags
    parameter(nflags = 5)
    INTEGER  FLAGS(nflags) !For pelagic kinetics (1) - not used; (2): 1-surface box, 0-below surface
                           ! sediments variables
    integer num_sed_flags                       
    parameter(num_sed_flags = 3)
    INTEGER  SED_FLAGS(nflags) !For BS kinetics (1) - not used; (2): 1-surface box, 0-below surface
                           ! sediments variables
                           
    character*80 file ! for indication of repeated run

    real    par_sed(nsconst) !(100)  !  BS model parmeters
    integer nsstate         !  number of botomm sediments state variables
    integer nsconst         !  number of parameters
    Integer NOSLAY          !  number of BS layers
    integer nsoutput        !  number of BS model outputs

    real state_sed(NOSLAY,nsstate)         ! BS state variable, g/m3 of pore water or solids or sediments
    real output_sed(NOSLAY,nsoutput)       ! BS model outputs (states and and auxilaries)
    integer          NUM_SED_OUTPUTS  !nsoutput


    integer   NUM_FLUX_RECEIVING_SED_LAYERS
    parameter(NUM_FLUX_RECEIVING_SED_LAYERS = 1) !Number of first layers that receives material from WC


    double precision DRIVING_FUNCTIONS_FOR_FLUXES(3) ! Fixme

    double precision SETTLING_RATES(nstate)     !g/m2/day

    integer          NUM_FLUXES_TO
    double precision FLUXES_TO_SEDIMENTS(nsstate) !g/m2/day

    !DOUBLE PRECISION SED_TEMPS(NOSLAY) ! BS Layers temperatures, C        
    double precision SALS(nkn,noslay) ! Salinity for Diffusion coeff., psu
    double precision TS               ! Temperature for Diffusion coeff., C

    INTEGER          NUM_FLUXES_FROM_SEDIMENTS
    DOUBLE PRECISION FLUXES_FROM_SEDIMENTS(nsstate) !g/m2/day


    DOUBLE PRECISION SURF_WATER_CONCS(nsstate)  ! WC concentrations to calculate difusion between WC and 1-st BS layer, g/m3

    INTEGER   NUM_SED_DRIV
    parameter(NUM_SED_DRIV = 1)
    DOUBLE PRECISION SED_DRIVING_FUNCTIONS(NOSLAY, NUM_SED_DRIV)

    DOUBLE PRECISION SED_DIFFUSIONS(NOSLAY, nsstate)

    DOUBLE PRECISION SED_MOD_1_ALUKAS_MOLDI_C   ! Function name for Diff. coeff. calculation

    DOUBLE PRECISION FLUXES_TO_ALUKAS(nstate)    ! Fluxes from BS to ALUKAS, g/m2/day



    DOUBLE PRECISION SED_MODEL_CONSTANTS(nsconst)

    INTEGER SEDIMENT_TYPE ! taken from grid file

    INTEGER NUM_NOT_DEPOSITED_FLUXES
    DOUBLE PRECISION FRACTION_OF_DEPOSITION(nstate)
    real             FRACTION_OF_DEPOSITION_FACTOR

    DOUBLE PRECISION NOT_DEPOSITED_FLUXES(nstate) !g/m2/day
    !DOUBLE PRECISION PART_MIXING_COEFFS (NOSLAY, nsstate)

    DOUBLE PRECISION PSTIME, TIME_STEP ! days

    integer NDIAGVAR_sed
    real INTER_sed(noslay,nsstate,NDIAGVAR_sed)      				!Intermediate variables for BS (processes)
    DOUBLE PRECISION PROCESSES_sed(noslay,nsstate, NDIAGVAR_sed)    !The same double precision

    !for setling fluxes management between layers
    DOUBLE PRECISION fluxes_from_upper(nstate)
    DOUBLE PRECISION fluxes_to_next(nstate)
    !save  fluxes_to_next

    integer isedi ! sediment transport on(>0), off(==0)
    real getpar !function for obtaining parameters values defined in *.str file

    ! For debugging purposes
    integer ix, iy, iz, error
    real x,y,z
    logical debug_stranger
    integer STRANGERS
    !For indication of first iteration
    integer ifirst_it
    
    ! INTERFACES    
        interface
            subroutine PELAGIC_KINETICS &
               (node_active,       nkn, &
                STATE_VARIABLES  , DERIVATIVES  , nstate, &
                MODEL_CONSTANTS  , nconst               , &
                DRIVING_FUNCTIONS, n_driving_functions  , &
                FLAGS            , nflags               , &
                PROCESS_RATES    , NDIAGVAR             , &
                SAVED_OUTPUTS    , n_saved_outputs      , &
                PH, &
                TIME, TIME_STEP  , CALLED_BEFORE, SURFACE_BOXES, ZOOP_OPTION_1, &
				ADVANCED_REDOX_OPTION)

                integer  :: nkn, nstate, nconst, n_driving_functions, nflags
                integer  :: n_saved_outputs, NDIAGVAR

                !internal numbers of nodes, used for diagnostics. Not implemented yet
                integer node_active(nkn)

                double precision, dimension(nkn,nstate),             intent(in)   :: STATE_VARIABLES
                double precision, dimension(nkn,nstate),             intent(out)  :: DERIVATIVES
                double precision, dimension(nconst),                 intent(in)   :: MODEL_CONSTANTS
                double precision, dimension(nkn,n_driving_functions),intent(in)   :: DRIVING_FUNCTIONS
                integer,          dimension(nflags),                 intent(in)   :: FLAGS
                double precision, dimension(nkn,nstate, NDIAGVAR),   intent(out)  :: PROCESS_RATES
                double precision, dimension(nkn)                                  :: PH
                double precision, dimension(nkn,n_saved_outputs),    intent(inout):: SAVED_OUTPUTS
                double precision, intent(in)                                      :: TIME
                double precision, intent(in)                                      :: TIME_STEP
                integer, intent(inout)                                            :: CALLED_BEFORE
                
                !Optional arguments
                integer, dimension(nkn), intent(in), optional                     :: SURFACE_BOXES
                integer                , intent(in), optional                     :: ZOOP_OPTION_1
                integer                , intent(in), optional                     :: ADVANCED_REDOX_OPTION
            end subroutine PELAGIC_KINETICS
        end interface

        interface
            subroutine SEDIMENT_MODEL_1 &
               (nkn,INIT_SED_STATE_VARS, SED_DEPTHS , SED_POROSITIES,  &
                SED_DENSITIES          , PART_MIXING_COEFFS         ,  &
                SED_DIFFUSIONS         , SURF_MIXLEN, SED_BURRIALS  ,  &
                SURF_WATER_CONCS       , SED_TEMPS                  ,  &
                NUM_SED_VARS           , NUM_SED_LAYERS             ,  &
                SED_MODEL_CONSTANTS    , NUM_SED_CONSTS             ,  &
                SED_DRIVING_FUNCTIONS  , NUM_SED_DRIV               ,  & ! not used yet
                SED_FLAGS              , NUM_SED_FLAGS              ,  &
                FLUXES_TO_SEDIMENTS    , NUM_FLUXES_TO_SEDIMENTS    ,  &
                NUM_FLUX_RECEIVING_SED_LAYERS, ADVECTIVE_VELOCITY   ,  &
                PSTIME, TIME_STEP                                   ,  &
                H_ERODEP                                            ,  &
                FINAL_SED_STATE_VARS                                ,  &
                FLUXES_FROM_SEDIMENTS, NUM_FLUXES_FROM_SEDIMENTS    ,  &
                PROCESSES_sed        , NDIAGVAR_sed                 ,  &
                SED_OUTPUTS          , NUM_SED_OUTPUTS              ,  &
                SED_SAVED_OUTPUTS    , NUM_SED_SAVED_OUTPUTS        ,  &
                SED_BURRIAL_RATE_OUTPUTS, ADVANCED_REDOX_OPTION)

                integer NUM_SED_VARS
                integer NUM_SED_LAYERS
                integer NUM_SED_CONSTS
                integer NUM_SED_DRIV
                integer NUM_FLUXES_TO_SEDIMENTS   ! for BS state variables
                integer NUM_FLUXES_FROM_SEDIMENTS ! for BS state variables
                integer NDIAGVAR_sed
                integer NUM_FLUX_RECEIVING_SED_LAYERS
                integer NUM_SED_OUTPUTS
                integer NUM_SED_FLAGS
                integer NUM_SED_SAVED_OUTPUTS

                !INPUT ARGUMENTS
                integer nkn      ! number of reactors (nodes)
                double precision INIT_SED_STATE_VARS  (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
                double precision SED_DEPTHS           (nkn,NUM_SED_LAYERS)
                double precision SED_POROSITIES       (nkn,NUM_SED_LAYERS)
                double precision SED_DENSITIES        (nkn,NUM_SED_LAYERS) !Bulk wet density
                double precision PART_MIXING_COEFFS   (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
                double precision SED_DIFFUSIONS       (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
                double precision SURF_MIXLEN
                double precision SED_BURRIALS         (nkn,NUM_SED_LAYERS)
                double precision SURF_WATER_CONCS     (nkn,NUM_SED_VARS)
                double precision SED_TEMPS            (nkn,NUM_SED_LAYERS) !Water temperature
                double precision SED_MODEL_CONSTANTS  (NUM_SED_CONSTS)
                double precision PROCESSES_sed        (nkn,NUM_SED_LAYERS, NUM_SED_VARS, NDIAGVAR_sed)
                double precision SED_DRIVING_FUNCTIONS(NUM_SED_LAYERS, NUM_SED_DRIV)
                double precision FLUXES_TO_SEDIMENTS  (nkn,NUM_FLUXES_TO_SEDIMENTS)
                double precision ADVECTIVE_VELOCITY
                double precision PSTIME
                double precision TIME_STEP
                double precision H_ERODEP             (nkn) ! Eroded(>0) or deposited(<0) thickness of first layer per time step, m

                !OUTPUT ARGUMENTS
                double precision FINAL_SED_STATE_VARS (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
                double precision FLUXES_FROM_SEDIMENTS(nkn,NUM_FLUXES_FROM_SEDIMENTS)
                double precision SED_OUTPUTS          (nkn,NUM_SED_LAYERS, NUM_SED_OUTPUTS) ! it is assigned to other array
                                                                                            ! with indices order required by hydrodynamics
                double precision, dimension(nkn,NUM_SED_LAYERS, NUM_SED_SAVED_OUTPUTS)   :: SED_SAVED_OUTPUTS
                integer, dimension(NUM_SED_FLAGS) :: SED_FLAGS
                
                !Optional arguments
                double precision, optional, dimension(nkn, NUM_SED_LAYERS, NUM_SED_VARS) :: SED_BURRIAL_RATE_OUTPUTS
				integer                , intent(in), optional                     :: ADVANCED_REDOX_OPTION
            end subroutine SEDIMENT_MODEL_1
        end interface 
                      
         
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    If(nlayers .ne. nlvdi) then
        print *, 'AQUABC: nlayers is not equal to nlvdi'
        stop
    end if
    
    If(nkn .ne. nkndi) then                       
        print *, 'AQUABC: nkn is not equal to nkndi'
        stop                                            
    end if 
                                           
    error = 0               !Indicator of errors for routine termination
    debug_stranger = .true. !True if debug for strangers
    !debug_stranger = .false.
    
    ! Definition of active nodes in each layer for pelagic model
    
    do layer = 1, nlayers
        nactive_nodes_arr(layer) = count(layer <= lmax_fast(1:nkn))
    end do

    !-----------------------------------------------------------
    ! Meaning of flags for pelagic kinetics :
    !  1 - safe mode indication 
    !  2 - surface layer indication
    !  3 - first iteration indication
    !  4 - ititial value for dissolved fraction FEII  calc. method
    !  5 - ititial value for dissolved fraction FEIII calc. method
    !
    !  Meaning of flags for BS kinetics :                     
    !   1 - first iteration indication                                                          
    !   2 - ititial value for dissolved fraction FEII  calc. method                             
    !   3 - ititial value for dissolved fraction FEIII calc. method                             
    !-----------------------------------------------------------
    
    ! Flags initialisation for subroutine PELAGIC_Kinetics
    FLAGS(1) = 0
    FLAGS(2) = 0
    FLAGS(3) = 0
    
    SED_FLAGS(1) = 0


    !first time iteration indication and necessary preparations
    if(ifirst_it .eq. 1) then
    
        !first iteration indication for pelagic and BS: 1-first it., 0-other iterations.         
        
         call getfnm('bioin',file)
         
         if (adjustl(trim(file)) .ne. 'INPUT/dump_wc.dat') then 
            FLAGS(3)   = 1
            SED_FLAGS(1)   = 1  
         end if                         
        
        !allocation of structure for saved outputs
        
        call alloc_saved_outputs_in_layers &
                      (nlayers, n_saved_outputs, nactive_nodes_arr)
        
        if(bsedim) then              
            allocate (SED_SAVED_OUTPUTS(nkn,noslay, num_sed_saved_outputs))
        end if          
        
        ifirst_it  = 0
    end if

    ! Flags for dissolved fraction calculation method for metals at first time iteration
    FLAGS(4) =  2  !   2: initial FE_II_DISS fraction from params,  1: from saturation 
    FLAGS(5) =  2  !   2: initial FE_III_DISS fraction from params, 1: from saturation
    
    SED_FLAGS(2) =  2  !   2: initial FE_II_DISS fraction from params,  1: from saturation 
    SED_FLAGS(3) =  2  !   2: initial FE_III_DISS fraction from params, 1: from saturation
     
    if(nstate.ne.30) then
        print *, 'AQUABC: Number of state variables is wrong', nstate
        stop
    end if
    
    ! Initialisation of arrays
    dg  (:,:,:,:)              = 0.
    e_fast1(:,:,:)             = 0.
    rates_fast(:,:,:)          = 0.
    dissolved_fractions_full(:,:,:) = 0.
    SED_DRIVING_FUNCTIONS(1:noslay,1:NUM_SED_DRIV)   = 0.
    SAVED_OUTPUTS(1:nkn,1:n_saved_outputs,1:nlayers) = 0.

    ! Asignment of some array bounds past to subroutines in order not to change their name
    NUM_FLUXES_FROM_SEDIMENTS = nstate   !nstate
    NUM_FLUXES_TO             = nsstate  !nsstate Number of fluxes from WC to BS
    NUM_NOT_DEPOSITED_FLUXES  = nstate   !nstate
    NUM_SED_OUTPUTS           = nsoutput !nsoutput

    ! To doble precision for call to pelagic model
    PSTIME = t
    TIME_STEP = dt
    MODEL_CONSTANTS(1:nconst) = par(1:nconst)



    ! Checking for NaNs in state vars
    if (debug_stranger) then
        do i = 1,nlayers
            do j = 1, nstate
                do k = 1, nkn
                    if(i .le. lmax_fast(k)) then
                        if(STRANGERS(e_fast0(k,j,i)) .eq. 1) then
                              !external node number
                              exn_node = ipv(k)
                              print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                              print *, 'AQUABC before call to pelagic:'
                              print *, 'Variable ',j,'Layer ',i,'Cell ',k, exn_node
                              print *, 'WC state variable is strange'
                              print *, 'State=',e_fast0(k,:,i)
                              print *, 'Layers=',lmax_fast(k)
                              print *, 'Layers depth=',depth_fast(k,1:lmax_fast(k))
                              error =1
                        end if
                    end if
                end do
            end do
        end do
        if (error .eq. 1) stop
    end if


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !       Processing water column and sinks to the bottom
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
    !------------------------------------------------------------
    !        Preparation of sinks from WC to BS
    !------------------------------------------------------------
    ! Settling velocities and dissolved fractions assigned in module now
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !       FRACTION OF DEPOSITED MATERIAL (TEMPORARY ASIGNMENT)
    !       Executed also in case of no botom sediments, because it is needed for
    !       FLX_ALUKAS_TO_SED_MOD_1


    ! Loop on nodes
    do nnode = 1,nkn
        call  get_nodal_area_code(nnode,SEDIMENT_TYPE) ! routine should be vectorized fixme
        exn_node = ipv(nnode)
        
        if(SEDIMENT_TYPE .lt. 0 .or. SEDIMENT_TYPE .gt. 10) then
            print *, 'aquabc: incorrect sediment type:', SEDIMENT_TYPE,'node:', exn_node
            stop
        end if
        
        SEDIMENT_TYPE_fast(nnode) = SEDIMENT_TYPE
    end do
    !End of loop on nodes
    
    ! case when sedtrans is switched off
    ! it is overwritten if sedtrans is switched on

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FRACTION_OF_DEPOSITION_FACTOR = 0.5
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !SEDIMENT_TYPE_fast(:)=0    !should be switched of when grid will be prepared with sediment types: fixme

    !Loop on state vaíables
    do i=1,nstate
        where     (SEDIMENT_TYPE_fast .eq. 0 .or. SEDIMENT_TYPE_fast .eq. 6) 
                                               
            !sand and Klaipeda strait
            FRACTION_OF_DEPOSITION_fast(:,i) = 0.06D+0
        elsewhere (SEDIMENT_TYPE_fast .eq. 1)
            !silty fine sand
            FRACTION_OF_DEPOSITION_fast(:,i) = 0.1D+0 !0.3D+0
        elsewhere (SEDIMENT_TYPE_fast .eq. 2 .or. SEDIMENT_TYPE_fast .eq. 3.or. &
                   SEDIMENT_TYPE_fast .eq. 4 .or. SEDIMENT_TYPE_fast .eq. 5)
            !coarse silt, clayey coarse silt, silty clayey mud, mud
            FRACTION_OF_DEPOSITION_fast(:,i) = 0.1D+0 !0.2D+0 !0.6D+0
        elsewhere (SEDIMENT_TYPE_fast .eq. 7 .or. SEDIMENT_TYPE_fast .eq. 10)
            ! Baltic sea and Nemunas delta
            FRACTION_OF_DEPOSITION_fast(:,i) = 0.D+0
        end where
    end do
    !End of loop on state vaíables
    
          
    ! 10 - the area type where no sediment transport takes place (Nemunas delta)
    ! Introduced to cope with situation when it is much erosion due to to low conc. in the river      


    !Multiplying by scale factor for deposition control
    FRACTION_OF_DEPOSITION_fast(:,:) = &
        FRACTION_OF_DEPOSITION_fast(:,:) * FRACTION_OF_DEPOSITION_FACTOR


    ! case when sedtrans is switched on
    isedi = nint(getpar('isedi'))

    if(isedi .gt. 0) then
        do i=1,nstate
            where     (H_ERODEP_fast .gt. 0.d0 )
                FRACTION_OF_DEPOSITION_fast(:,i) = 0.
                !> 0 - erosion, < 0 - deposition
            elsewhere (H_ERODEP_fast .le. 0.d0 )
                !FRACTION_OF_DEPOSITION_fast(:,i) = 1.
                ! for experiment taken as was without sediment transport:
                 FRACTION_OF_DEPOSITION_fast(:,i) = FRACTION_OF_DEPOSITION_fast(:,i) !fixme
            end where
        end do
    end if

    !------------------------------------------------------------
    !     End of preparation sinks from WC to BS
    !------------------------------------------------------------




    !llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
    !       Loop on WC layers
    !llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
    do layer = 1, nlayers

        if(layer .eq. 1) then
            FLAGS(2) = 1
        else
            FLAGS(2) = 0
        end if
        
        states_fast(1:nkn,1:nstate) = e_fast0(1:nkn,1:nstate,layer)
        depth_layer(1:nkn)          = depth_fast(1:nkn,layer)


        ! Light management
        !       current layer
        light_fast(:) = 0.

	    if (layer .eq. 1) then
	        light_fast(1:nkn) = ITOT_fast(1:nkn)
	    else
	        where(layer <= lmax_fast(1:nkn))
                light_fast(1:nkn) = light_next(1:nkn)
            end where
	    end if

        !next layer
        if (layer .lt. nlayers) then

            call layer_light(nkn,states_fast,    &
                             nstate,             &
                             par,                &
                             nconst,             &
                             layer+1,            &
                             lmax_fast,          &
                             depth_layer,        &
                             light_fast, light_next)
        end if

        ! Preparation of arrays for call to pelagic model
        ! without nodes that does  not exist for the current layer

        ! number of active nodes in current layer
        nactive_nodes = count(layer <= lmax_fast(1:nkn))

        allocate(node_active_num    (nactive_nodes))
        allocate(STATE_fast         (nactive_nodes,nstate))       
        allocate(DERIVATIVES_fast   (nactive_nodes,nstate))        
        allocate(PH                 (nactive_nodes))
        allocate(DGS                (nactive_nodes,nstate,ndiagvar))
        allocate(DRIVING_FUNCTIONS  (nactive_nodes,10))
        allocate(SAVED_OUTPUTS_fast (nactive_nodes,n_saved_outputs))
        !allocate(WC_OUTPUTS         (nactive_nodes,noutput))
        
        STATE_fast       (:,:)  = 0.
        DERIVATIVES_fast (:,:)  = 0.
        PH (:)                  = 0.
        DGS (:,:,:)             = 0.
        DRIVING_FUNCTIONS(:,:)  = 0.
        SAVED_OUTPUTS_fast(:,:) = 0.
        !WC_OUTPUTS       (:,:) = 0.
        
        j=1
      
        ! Loop on nodes
        do k=1,nkn

            if(layer <= lmax_fast(k)) then
                node_active_num(j)  = k
                STATE_fast(j,1:nstate)  = e_fast0(k,1:nstate,layer)
                !WC_OUTPUTS(j,1:noutput) = wc_output(layer,k,1:noutput)

                ! Driving functions management
                call para_get_value('KE', KE)
                DRIVING_FUNCTIONS(j, 1) = wtemp_fast(k,layer)  ! w. temp.
                DRIVING_FUNCTIONS(j, 2) = sal_fast  (k,layer)  ! salinity
                DRIVING_FUNCTIONS(j, 3) = light_fast(k)        ! Instantaneous light W/m2, conversion to langlays is done in pelagic kinetics
                DRIVING_FUNCTIONS(j, 4) = FDAY                 ! fraction of day with light. 1 - for instantenous light, just for compatibility with older code
                DRIVING_FUNCTIONS(j, 5) = atemp_fast(k)        ! air temp.
                DRIVING_FUNCTIONS(j, 6) = wind_fast (k)        ! wind spead
                DRIVING_FUNCTIONS(j, 7) = 0.0D0                ! elevation
                DRIVING_FUNCTIONS(j, 8) = depth_fast(k,layer)  ! water depth  or layer thickness in case of 3d
                DRIVING_FUNCTIONS(j, 9) = KE   !0.8D0          ! background light extinction coefficient. Later should be calculated. fixme
                                                               ! should be compatible with value in subroutine alight_alukas
                DRIVING_FUNCTIONS(j,10) = ice_cover_fast(k)    ! ice cover - area fraction of reactor
                                                               ! was reserved for pH. Not necessary when co2sys was introduced

                !end of driving functions section

                j=j+1

            end if ! layer <=
        end do ! k (nodes)
        ! End of loop on nodes

        SAVED_OUTPUTS_fast(1:nactive_nodes,1:n_saved_outputs) = &
                           SAVED_OUTPUTS_IN_LAYERS(layer) % MATRIX(1:nactive_nodes,1:n_saved_outputs)

        CALL INIT_PELAGIC_MODEL_CONSTANTS
											  
        CALL PELAGIC_KINETICS &
                (node_active_num, nactive_nodes,&
                 STATE_fast  , DERIVATIVES_fast, nstate,  &
                 MODEL_CONSTANTS  , nconst,               &
                 DRIVING_FUNCTIONS, n_driving_functions,  &
                 FLAGS            , nflags,               &
                 DGS              , NDIAGVAR,             &
                 SAVED_OUTPUTS_fast    , n_saved_outputs, &
                 pH ,               &
                 PSTIME, TIME_STEP, CALLED_BEFORE)



        do k = 1, nactive_nodes
            rates_fast(node_active_num(k),1:nstate,layer) = DERIVATIVES_fast(k,1:nstate)
            ph_full(node_active_num(k),layer) = PH(k)
            dg(layer,node_active_num(k),1:nstate,:)       = DGS(k,1:nstate,:)
            SAVED_OUTPUTS(node_active_num(k),1:n_saved_outputs,layer) = SAVED_OUTPUTS_fast(k,1:n_saved_outputs)            
        end do
          
        SAVED_OUTPUTS_IN_LAYERS(layer) % MATRIX(1:nactive_nodes,1:n_saved_outputs) = &
            SAVED_OUTPUTS_fast(1:nactive_nodes,1:n_saved_outputs)

        if (debug_stranger) then
            do j = 1, nstate
               do k = 1, nkn
                   if(layer <= lmax_fast(k)) then
                       if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
                           !external node number
                           exn_node = ipv(k)
                           print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                           print *, 'AQUABC after call to pelagic kinetics:'
                           print *, 'Variable ',j,'Cell ',k, exn_node
                           print *, 'WC derivative is NaN'
                           print *, 'Derivative=',rates_fast(k,j,layer)
                           print *, 'Layers=',lmax_fast(k)
                           print *, 'Layers depth=',depth_fast(k,1:lmax_fast(k))
                           error =1
                       end if
                   end if
               end do
            end do
         
            if (error .eq. 1) stop
        end if


        !  SETTLING FLUXES MANAGAMENT ALSO FOR 3D CASE
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        fluxes_from_upper_fast(1:nkn,1:nstate) = 0.

        if (layer > 1) then
            do k=1,nkn
                if(layer <= lmax_fast(k)) then
                    fluxes_from_upper_fast(k,1:nstate) = &
                        fluxes_to_next_fast(k,1:nstate)
                end if
            end do
        end if


        ! Particulate material fluxes to and from the  layers
        ! settling rates comes in g/m2/day

        fluxes_to_next_fast(1:nkn,1:nstate) = 0.


        ! Loop on nodes
        do k=1,nkn
            nnode = k
            !---------------------------------------------------------

            if(nstate .ne. 30) then
                print *, 'AQUABC: Wrong number of state variables for'
                print *, ' for DISSOLVED FRACTIONS'
                stop
            end if

            ! Dissolved fractions for WC variables
            dissolved_fractions_full(k,NH4_N_INDEX         ,layer) = 1.0D+0         ! AMMONIUM NITROGEN
            dissolved_fractions_full(k,NO3_N_INDEX         ,layer) = 1.0D+0         ! NITRATE NITROGEN
            dissolved_fractions_full(k,PO4_P_INDEX         ,layer) = 1.0D+0* SAVED_OUTPUTS(k,5,layer) ! ORTHOPHOSPHATE PHOSPHORUS
            dissolved_fractions_full(k,DISS_OXYGEN_INDEX   ,layer) = 1.0D+0         ! DISSOLVED OXYGEN
            dissolved_fractions_full(k,DIA_C_INDEX         ,layer) = 0.0D+0         ! DIATOMS CARBON
            dissolved_fractions_full(k,ZOO_C_INDEX         ,layer) = 0.0D+0         ! ZOOPLANKTON CARBON
            dissolved_fractions_full(k,ZOO_N_INDEX         ,layer) = 0.0D+0         ! ZOOPLANKTON CARBON
            dissolved_fractions_full(k,ZOO_P_INDEX         ,layer) = 0.0D+0         ! ZOOPLANKTON PHOSPHORUS
            dissolved_fractions_full(k,DET_PART_ORG_C_INDEX,layer) = 0.0D+0         ! DETRITUS PARTICULATE ORG. CARBON
            dissolved_fractions_full(k,DET_PART_ORG_N_INDEX,layer) = 0.0D+0         ! DETRITUS PARTICULATE ORG. NITROGEN
            dissolved_fractions_full(k,DET_PART_ORG_P_INDEX,layer) = 0.0D+0         ! DETRITUS PARTICULATE ORG. PHOSPHORUS
            dissolved_fractions_full(k,DISS_ORG_C_INDEX    ,layer) = 1.0D+0         ! DISSOLVED ORGANIC CARBON
            dissolved_fractions_full(k,DISS_ORG_N_INDEX    ,layer) = 1.0D+0         ! DISSOLVED ORGANIC NITROGEN
            dissolved_fractions_full(k,DISS_ORG_P_INDEX    ,layer) = 1.0D+0         ! DISSOLVED ORGANIC PHOSPHORUS
            dissolved_fractions_full(k,CYN_C_INDEX         ,layer) = 0.0D+0         ! CYANOBACTERIA CARBON
            dissolved_fractions_full(k,OPA_C_INDEX         ,layer) = 0.0D+0         ! OTHER PHYTOPLANKTON CARBON
            dissolved_fractions_full(k,DISS_Si_INDEX       ,layer) = 1.0D+0         ! DISSOLOVED SILICA
            dissolved_fractions_full(k,PART_Si_INDEX       ,layer) = 0.0D+0         ! BIOGENIC SILICA
            dissolved_fractions_full(k,FIX_CYN_C_INDEX     ,layer) = 0.0D+0         ! Nitrogen fixing cyano-bateria
            dissolved_fractions_full(k,INORG_C_INDEX       ,layer) = 1.0D+0         ! Inorganic carbon (dissolved in this version)
            dissolved_fractions_full(k,TOT_ALK_INDEX       ,layer) = 1.0D+0         ! Alcalinitity
            dissolved_fractions_full(k,FE_II_INDEX         ,layer) = 1.0D+0 * SAVED_OUTPUTS(k,1,layer) ! Iron2+      : Special treatment because of the possible dissloved and particulate feactions of ions
            dissolved_fractions_full(k,FE_III_INDEX        ,layer) = 1.0D+0 * SAVED_OUTPUTS(k,2,layer) ! Iron3+      : Special treatment because of the possible dissloved and particulate feactions of ions
            dissolved_fractions_full(k,MN_II_INDEX         ,layer) = 1.0D+0 * SAVED_OUTPUTS(k,3,layer) ! Manganese2+ : Special treatment because of the possible dissloved and particulate feactions of ions
            dissolved_fractions_full(k,MN_IV_INDEX         ,layer) = 1.0D+0 * SAVED_OUTPUTS(k,4,layer) ! Manganese4+ : Special treatment because of the possible dissloved and particulate feactions of ions
            dissolved_fractions_full(k,CA_INDEX            ,layer) = 1.D+0
            dissolved_fractions_full(k,MG_INDEX            ,layer) = 1.D+0
            dissolved_fractions_full(k,S_PLUS_6_INDEX      ,layer) = 1.D+0
            dissolved_fractions_full(k,S_MINUS_2_INDEX     ,layer) = 1.D+0
            dissolved_fractions_full(k,CH4_C_INDEX         ,layer) = 1.D+0


            if(layer .eq. lmax_fast(k)) then   !processing nodes which current layer is nearbottom layer
                STATE_VARIABLES       (1:nstate) = e_fast0                    (k,1:nstate,layer)
                FRACTION_OF_DEPOSITION(1:nstate) = FRACTION_OF_DEPOSITION_fast(k,1:nstate)
                DISSOLVED_FRACTIONS   (1:nstate) = dissolved_fractions_full   (k,1:nstate,layer)

                CALL FLX_ALUKAS_II_TO_SED_MOD_1 &
                       (STATE_VARIABLES               , nstate,              &
                        MODEL_CONSTANTS               , nconst,              &
                        DRIVING_FUNCTIONS_FOR_FLUXES  , 3,                   &
                        SETTLING_VELOCITIES, DISSOLVED_FRACTIONS,            & !inputs
                        1.0D+0 , nnode, PSTIME,                              & !bottom factor
                        SETTLING_RATES, FLUXES_TO_SEDIMENTS, NUM_FLUXES_TO,  &
                        SEDIMENT_TYPE , FRACTION_OF_DEPOSITION,              &
                        NOT_DEPOSITED_FLUXES, NUM_NOT_DEPOSITED_FLUXES)

                not_deposited_fluxes_fast(k, 1:nstate)  = NOT_DEPOSITED_FLUXES(1:nstate)
                settling_rates_fast      (k, 1:nstate)  = SETTLING_RATES      (1:nstate)
                FLUXES_TO_SEDIMENTS_fast (k, 1:nsstate) = FLUXES_TO_SEDIMENTS (1:nsstate)


                rates_fast(k,1:nstate,layer) = rates_fast(k,1:nstate,layer)                                  &
                                               - settling_rates_fast(k,1:nstate)       /depth_fast(k,layer)         &
                                               + not_deposited_fluxes_fast(k, 1:nstate)/depth_fast(k,layer)  &
                                               + fluxes_from_upper_fast(k,1:nstate)    /depth_fast(k,layer)
                
                dg(layer,k,1:nstate,24) = settling_rates_fast      (k,1:nstate)
                dg(layer,k,1:nstate,25) = not_deposited_fluxes_fast(k,1:nstate)

                ! Checking for strange values in WC derivatives
                if (debug_stranger) then
                    do j = 1, nstate
                        if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
                            !external node number
                            exn_node = ipv(k)
                            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                            print *, 'AQUABC after updating derivatives (bottom layers):'
                            print *, 'Variable ',j,' Node ',k, exn_node,' Layer ',layer
                            print *, 'Derivative is NaN'
                            print *, 'Derivative = '          ,  rates_fast(k,j,layer)
                            print *, 'State = '               ,  STATE_fast(k,j)
                            print *, 'Settling_rates:'        ,  settling_rates_fast(k,1:nstate)
                            print *, 'DISS. FRAC. for var. 25',  SAVED_OUTPUTS(k,1,layer)
                            print *, 'DISS. FRAC. for var. 26',  SAVED_OUTPUTS(k,2,layer)
                            print *, 'DISS. FRAC. for var. 27',  SAVED_OUTPUTS(k,3,layer)
                            print *, 'DISS. FRAC. for var. 28',  SAVED_OUTPUTS(k,4,layer)
                            print *, 'Not_deposited_fluxes:'  ,  not_deposited_fluxes_fast(k, j)
                            print *, 'Fluxes_from_upper: '    ,  fluxes_from_upper_fast(k,1:nstate)
                            print *, 'Depth:   '              ,  depth_fast(k,layer)
                            error =1
                        end if
                    end do
                    
                    if (error .eq. 1) stop
                end if

            end if ! end of processing nodes which current layer is nearbottom layer

            !------------------------------------------------------------
            if(layer .lt. lmax_fast(k)) then !processing nodes which current layer is not nearbottom layer
            
                STATE_VARIABLES     (1:nstate) = e_fast0                 (k,1:nstate,layer)
                DISSOLVED_FRACTIONS (1:nstate) = dissolved_fractions_full(k,1:nstate,layer)
                
                ! Settling velocities and dissolved fractions initialised
                ! using module aquabc_II_sed_ini
                
                call SETTLING_TO_THE_NEXT_LAYER &
                      (STATE_VARIABLES    , nstate   ,           &
                       SETTLING_VELOCITIES, DISSOLVED_FRACTIONS, &
                       nnode, LAYER, PSTIME,                     &
                       fluxes_to_next)
                
                fluxes_to_next_fast(k,1:nstate) = fluxes_to_next(1:nstate)
                
                rates_fast(k,1:nstate,layer)    = rates_fast(k,1:nstate,layer)                      &
                                         - fluxes_to_next_fast   (k,1:nstate)/depth_fast(k,layer) &
                                         + fluxes_from_upper_fast(k,1:nstate)/depth_fast(k,layer)
                
                dg(layer,k,1:nstate,24) = fluxes_to_next(1:nstate)

                ! Checking for strange values in WC derivatives
                if (debug_stranger) then
                    do j = 1, nstate
                        if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
                            !external node number
                            exn_node = ipv(k)
                            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                            print *, 'AQUABC after updating derivatives (upper layers):'
                            print *, 'Variable ',j,' Cell ',k, exn_node,' Layer ',layer
                            print *, 'Derivative is NaN'
                            print *, 'Derivative=',rates_fast(k,j,layer)
                            print *, 'Depth:   '  , depth_fast(k,layer)
                            print *, 'Max number of layers:', lmax_fast(k)
                            error =1
                        end if
                    end do
                    if (error .eq. 1) stop
                end if
                !-----------------------------------------------------------

            end if
            ! end of processing nodes which current layer is not nearbottom

        end do 
        ! Loop on nodes
        
        ! End of processing particulate material fluxes to and from the  layers


        if(allocated(node_active_num ))  deallocate(node_active_num)
        if(allocated(STATE_fast ))       deallocate(STATE_fast)
        !if(allocated(WC_OUTPUTS ))       deallocate(WC_OUTPUTS)
        if(allocated(DRIVING_FUNCTIONS ))deallocate(DRIVING_FUNCTIONS)
        if(allocated(DERIVATIVES_fast )) deallocate(DERIVATIVES_fast)
        if(allocated(DGS))               deallocate(DGS)
    end do !layers


    !llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
    !       End of loop on WC layers
    !llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
    
    ! Checking for NaN, Inf and strange values in WC derivatives
    if (debug_stranger) then
        do layer = 1,nlayers
            do j = 1, nstate
                do k = 1, nkn
                    if(layer <= lmax_fast(k)) then
                        if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
                            ! external node number
                            exn_node = ipv(k)
                            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                            print *, 'AQUABC after updating derivatives by settling etc:'
                            print *, 'Variable ',j,'Layer ',layer,'Cell ', k, exn_node
                            print *, 'WC Derivative is NaN'
                            print *, 'Derivative=',rates_fast(k,j,layer)
                            error =1
                        end if
                    end if
                end do
            end do
        end do
        if (error .eq. 1) stop
    end if

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     BOTTOM SEDIMENTS PROCESSING
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    
    ! In this section it is assummed that variable 'layer' is a nearbottom WC layer but has a value 1 so
    ! only 2d calculations are performed. fixme

    if(bsedim)  then

        !Checking for strange values in BS initial state
        if (debug_stranger) then
            do j = 1, nsstate
                do i = 1, noslay
                    do k=1, nkn
                        flux = INIT_SED_STATE_VARS_fast(k,i,j)
                        if(STRANGERS(flux) .eq. 1) then
                            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                            print *, 'AQUABC: Layer ', i, 'Variable ',j,'Cell ',k
                            print *, 'BS initial state is NaN'
                            print *, 'INITIAL(k,i,j)=',INIT_SED_STATE_VARS_fast(k,i,j)
                            error = 1
                        end if
                    end do
                end do
            end do  !
            if(error.eq. 1) stop
        end if
        
        ! ************************************************************************
        ! Initialisation of sediment properties is done by routine initialisation
        ! routine from module aquabc_II_sed_ini
        ! Later should be read from files. fixme

        if (noslay .ne. 7) then
            print *, 'AQUABC_II:'
            print *, 'Number of BS layers is not equal 7 but', noslay
            print *, 'Change its value in aquabc_II'
            stop
        end if

        ! Porosities and densities are processed from sediment transport model.
        
        !--------------------------------------------------------------
        
        
        !          Initilized in module aquabc_II_sed_ini now
                   !ADVECTIVE_VELOCITY = 0. !1.0D-5 !m/s
        
        !          Surface mixing length for exchange with WC
        !          Initilized in module aquabc_II_sed_ini now
                   !SURF_MIXLEN = 0.10  !m 0.1
        ! ********************************************************************
        
        
        ! Checking for NaN, Inf and strange values in WC derivatives
        if (debug_stranger) then
            do layer = 1,nlayers
                do j = 1, nstate
                    do k = 1, nkn
                        if(layer <= lmax_fast(k)) then
                            if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
                                !external node number
                                exn_node = ipv(k)
                                print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                                print *, 'AQUABC before diffusion coeff.:'
                                print *, 'Variable ',j,'Layer ',layer,'Cell ', k, exn_node
                                print *, 'Derivative is NaN'
                                print *, 'Derivative=',rates_fast(k,j,layer)
                                error =1
                            end if
                        end if
                    end do
                end do
            end do
            if (error .eq. 1) stop
        end if

        !Calculation of diffusion coefficients
        do k =1,nkn
            do i=1,noslay              
                SED_TEMPS_fast(k,i) = wtemp_fast(k,lmax_fast(k)) ! sed. temperatures for sediments subroutine.
                TS                  = SED_TEMPS_fast(k,i)        ! Assumed equal to water column temp.
                SALS(k,i)           = sal_fast  (k,lmax_fast(k)) ! sed. salinity. Assumed equal to water column salinity (last layer)
                     
                !do j=1,nsstate
                !    SED_DIFFUSIONS_fast(k,i,j) = &
                !              SED_MOD_1_ALUKAS_MOLDI_C(j, TS, SALS(k,i), 0.0D+0)*86400. ! m2/day
                !end do
                          
            end do
        end do

        
        

        call SED_MOD_1_ALUKAS_MOLDI_C_VEC(SED_TEMPS_fast,  SALS, nkn, noslay, nsstate, SED_DIFFUSIONS_fast) 
        SED_DIFFUSIONS_fast = SED_DIFFUSIONS_fast*86400. ! m2/day
!         SED_TEMPS_fast      = 12.0D0
!         SED_DIFFUSIONS_fast = 0.0864000d0

        !  WC concentrations corresponding to BS variables for molecular diffusion
        !  Change these statments if water column vars number is changed

        if(nsstate.ne.24) then
            print *, 'AQUABC: To get correctly values for Surface water concentrations'
            print *, 'number of sediment state variables should be 24 but is', nsstate
            stop
        end if

        do k=1,nkn
            SURF_WATER_CONCS_fast(k, 1)  =  e_fast0(k, NH4_N_INDEX,     lmax_fast(k)) !state(1)
            SURF_WATER_CONCS_fast(k, 2)  =  e_fast0(k, NO3_N_INDEX,     lmax_fast(k)) !state(2)
            SURF_WATER_CONCS_fast(k, 3)  =  e_fast0(k, DISS_ORG_N_INDEX,lmax_fast(k)) !state(13)
            SURF_WATER_CONCS_fast(k, 4)  =  0.0D+0  ! Should be fixed in case of particle mixing is switched on in BS
            SURF_WATER_CONCS_fast(k, 5)  =  e_fast0(k, PO4_P_INDEX,     lmax_fast(k)) !state(3)
            SURF_WATER_CONCS_fast(k, 6)  =  e_fast0(k, DISS_ORG_P_INDEX,lmax_fast(k)) !state(14)
            SURF_WATER_CONCS_fast(k, 7)  =  0.0D+0  ! Should be fixed in case of particle mixing is switched on in BS
            SURF_WATER_CONCS_fast(k, 8)  =  e_fast0(k, DISS_OXYGEN_INDEX,lmax_fast(k)) !state(4)
            SURF_WATER_CONCS_fast(k, 9)  =  e_fast0(k, DISS_ORG_C_INDEX,lmax_fast(k)) !state(12)
            SURF_WATER_CONCS_fast(k,10)  =  0.0D+0   ! Should be fixed in case of particle mixing is switched on in BS
            SURF_WATER_CONCS_fast(k,11)  =  e_fast0(k, DISS_Si_INDEX,   lmax_fast(k)) !state(20)
            SURF_WATER_CONCS_fast(k,12)  =  0.0D+0   ! Should be fixed in case of particle mixing is switched on in BS
            SURF_WATER_CONCS_fast(k,13)  =  e_fast0(k, INORG_C_INDEX,   lmax_fast(k)) !state(23)
            SURF_WATER_CONCS_fast(k,14)  =  e_fast0(k, TOT_ALK_INDEX,   lmax_fast(k)) !state(24)
            SURF_WATER_CONCS_fast(k,15)  =  sal_fast(k,lmax_fast(k))   !Salinity
            
            SURF_WATER_CONCS_fast(k,16)  =  e_fast0(k, FE_II_INDEX,     lmax_fast(k)) !state(25)
            SURF_WATER_CONCS_fast(k,17)  =  e_fast0(k, FE_III_INDEX,    lmax_fast(k)) !state(26)
            SURF_WATER_CONCS_fast(k,18)  =  e_fast0(k, MN_II_INDEX,     lmax_fast(k)) !state(27)
            SURF_WATER_CONCS_fast(k,19)  =  e_fast0(k, MN_IV_INDEX,     lmax_fast(k)) !state(28)
            
            SURF_WATER_CONCS_fast(k,20)  =  e_fast0(k, CA_INDEX,        lmax_fast(k)) !state(29)
            SURF_WATER_CONCS_fast(k,21)  =  e_fast0(k, MG_INDEX,        lmax_fast(k)) !state(30)
            SURF_WATER_CONCS_fast(k,22)  =  e_fast0(k, S_PLUS_6_INDEX,  lmax_fast(k)) !state(31)
            SURF_WATER_CONCS_fast(k,23)  =  e_fast0(k, S_MINUS_2_INDEX, lmax_fast(k)) !state(32)
            SURF_WATER_CONCS_fast(k,24)  =  e_fast0(k, CH4_C_INDEX,     lmax_fast(k)) !state(33)
        end do

        ! Checking for NaN, Inf and strange values in WC derivatives
        if (debug_stranger) then
            do layer = 1,nlayers
                do j = 1, nstate
                    do k = 1, nkn
                        if(layer <= lmax_fast(k)) then
                            if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
                                !external node number
                                exn_node = ipv(k)
                                print *, 'AQUABC after surf. conc:'
                                print *, 'Variable ',j,'Layer ',layer,'Cell ', k, exn_node
                                print *, 'WC Derivative is NaN'
                                print *, 'Derivative=',rates_fast(k,j,layer)
                                error =1
                            end if
                        end if
                    end do
                end do
            end do
            if (error .eq. 1) stop
        end if

        !BS constants
        do i=1,nsconst
            SED_MODEL_CONSTANTS(i) = par_sed(i)
        end do

        !BS state vars initial values
        !INIT_SED_STATE_VARS_fast(:,:,:) = 0.!state_sed_fast(:,:,:)

        ! Checking for NaN, Inf and strange values in BS initial state vars
        if (debug_stranger) then
            do i = 1, nsstate
                do j = 1,noslay
                    do k = 1, nkn
                      flux = INIT_SED_STATE_VARS_fast(k,j,i)
                      if(STRANGERS(flux) .eq. 1) then
                          ! external node number
                          exn_node = ipv(k)
                          print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                          print *, 'AQUABC before BS kinetics:'
                          print *, 'Variable ',i,'Layer ',j,'Cell ', k, exn_node
                          print *, 'BS state variable is NaN'
                          print *, 'Value=',INIT_SED_STATE_VARS_fast(k,j,i)
                          error =1
                      end if
                    end do
                end do
            end do
            if (error .eq. 1) stop
        end if

        ! Checking for NaN, Inf and strange values in WC derivatives
        if (debug_stranger) then
            do layer = 1,nlayers
                do j = 1, nstate
                    do k = 1, nkn
                        if(layer <= lmax_fast(k)) then
                            if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
                                ! external node number
                                exn_node = ipv(k)
                                print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                                print *, 'AQUABC befor call to BS:'
                                print *, 'Variable ',j,'Layer ',layer,'Cell ', k, exn_node
                                print *, 'WC Derivative is NaN'
                                print *, 'Derivative=',rates_fast(k,j,layer)
                                error =1
                            end if
                        end if
                    end do
                end do
            end do
            if (error .eq. 1) stop
        end if

         !call resuspension(nkn,nlayers,lmax_fast, &
         !                  waves_fast,depth_fast,rmass_fast)


        !Sediment drivig functions not used yet
        !Call to the main sediment routine

        call INIT_BSED_MODEL_CONSTANTS()

        call SEDIMENT_MODEL_1 (nkn,  &
                  INIT_SED_STATE_VARS_fast  , SED_DEPTHS_fast , SED_POROSITIES_fast,&
                   SED_DENSITIES_fast        , PART_MIXING_COEFFS_fast         ,  &
                   SED_DIFFUSIONS_fast       , SURF_MIXLEN, SED_BURRIALS_fast  ,  &
                   SURF_WATER_CONCS_fast     , SED_TEMPS_fast                  ,  &
                   nsstate                   , NOSLAY                          ,  &
                   SED_MODEL_CONSTANTS       , nsconst                         ,  &
                   SED_DRIVING_FUNCTIONS     , NUM_SED_DRIV                    ,  &
                   SED_FLAGS                 , NUM_SED_FLAGS                   ,  &
                   FLUXES_TO_SEDIMENTS_fast  , NUM_FLUXES_TO                   ,  &
                   NUM_FLUX_RECEIVING_SED_LAYERS, ADVECTIVE_VELOCITY         ,  &
                   PSTIME, TIME_STEP         ,                               &
                   H_ERODEP_fast             ,                               &
                   FINAL_SED_STATE_VARS_fast ,                               &
                   FLUXES_FROM_SEDIMENTS_fast, NUM_FLUXES_FROM_SEDIMENTS  ,  &
                   PROCESSES_sed_fast        , NDIAGVAR_sed               ,  &
                   SED_OUTPUTS_fast          , NUM_SED_OUTPUTS            ,  &
                   SED_SAVED_OUTPUTS         , NUM_SED_SAVED_OUTPUTS         &
                   )
        
        ! Correcting slightly negaative or small values              !                                          1 - CO2SYS is called with one dimensional array
        where ( FINAL_SED_STATE_VARS_fast .lt. 1.0D-10 )        
         FINAL_SED_STATE_VARS_fast = 1.0D-10        
        end where
        
        ! reshaping arrays
        do j=1,nsstate
            do i=1,noslay
                do k=1,nkn
                    do l=1,ndiagvar_sed
                        INTER_sed_full (i,k,j,l) = PROCESSES_sed_fast(k,i,j,l)
                    end do
                end do
            end do
        end do

        do j=1,nsoutput
            do i=1,noslay
                do k=1,nkn
                    output_sed_full(i,k,j)  = SED_OUTPUTS_fast(k,i,j)
                end do
            end do
        end do
   
        ! Checking for strange values in BS final state
        if (debug_stranger) then
            do j = 1, nsstate
                do k=1,nkn
                    do i = 1, noslay
                        flux = FINAL_SED_STATE_VARS_fast(k,i,j)
                        if(STRANGERS(flux) .eq. 1) then
                            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                            print *, 'AQUABC after call to BS:'
                            print *, 'Layer ', i, 'Variable ',j,'Cell ',k
                            print *, 'Final state is NaN'
                            print *, 'FINAL(k,i,j)=',FINAL_SED_STATE_VARS_fast(k,i,j)
                            error =1
                        end if
                    end do
                end do
            end do
            if(error .eq. 1) stop
        end if

        ! Checking for strange values in fluxes from BS
        if (debug_stranger) then
            do i = 1, nsstate
                do k=1,nkn
                    flux = FLUXES_FROM_SEDIMENTS_fast(k,i)
                    if(STRANGERS(flux) .eq. 1) then
                        print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                        print *, 'AQUABC: Layer ', i, 'Variable ',i,'Cell ',k
                        print *, 'Flux from BS to WC by BS vars is NaN'
                        print *, 'Flux(k,i)=',FLUXES_FROM_SEDIMENTS_fast(k,i)
                        error =1
                    end if
                end do
            end do
            if(error .eq. 1) stop
        end if
        
        !  Fluxes from BS to WC
        do k=1,nkn
            FLUXES_FROM_SEDIMENTS(:) = FLUXES_FROM_SEDIMENTS_fast(k,:)
            
            CALL FLX_SED_MOD_1_TO_ALUKAS_II                            &
                   (FLUXES_FROM_SEDIMENTS, NUM_FLUXES_FROM_SEDIMENTS, &
                    FLUXES_TO_ALUKAS     , nstate)
            
            FLUXES_TO_ALUKAS_fast(k,:) = FLUXES_TO_ALUKAS(:)
        end do
        
        ! Checking for strange values in fluxes to WC from BS
        if (debug_stranger) then
            do i = 1, nstate
                do k=1,nkn
                    flux = FLUXES_TO_ALUKAS_fast(k,i)
                    if(STRANGERS(flux) .eq. 1) then
                        print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                        print *, 'AQUABC:'
                        print *, 'Layer ', i, 'Variable ',i,'Cell ',k
                        print *, 'Flux from BS to WC by WC vars is NaN'
                        print *, 'Flux(k,i)=',FLUXES_TO_ALUKAS_fast(k,i)
                        error =1
                    end if
                end do
            end do
            if(error .eq. 1) stop
        end if
        
        ! Sending fluxes to intermediate variables
        do k=1,nkn
            do i = 1, noslay
                do j = 1, nsstate
                
                    if (i .eq. 1) then
                        INTER_sed_full(i,k,j, 8) =  FLUXES_TO_SEDIMENTS_fast(k,j)
                        INTER_sed_full(i,k,j, 9) =  FLUXES_FROM_SEDIMENTS_fast(k,j)
                    else
                        INTER_sed_full(i,k,j, 8) = 0.D+0
                        INTER_sed_full(i,k,j, 9) = 0.D+0
                    end if
                
                end do
            end do
        end do

        ! Updating WC bottom layer derivatives
        
        ! Checking for strange values in WC derivatives
        if (debug_stranger) then
            do i = 1, nstate
                do k=1,nkn
                    ! print  *, 'node=', k, 'Variable ',i,'Layer' ,lmax_fast(k),'rate ',rates_fast(k,i,lmax_fast(k))
                    if(STRANGERS(rates_fast(k,i,lmax_fast(k))) .eq. 1) then
                        print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                        print *, 'AQUABC before updating by fluxes from from BS: '
                        print *,  'Variable ',i,'Node ',k, 'Layer ',lmax_fast(k)
                        print *, 'WC bottom layer derivative is NaN'
                        print *, 'Derivative',rates_fast(k,i,lmax_fast(k))
                        print *, 'Flux from BS:', FLUXES_TO_ALUKAS_fast(k,i)
                        print *, 'Layer depth: ', depth_fast(k,lmax_fast(k))
                        error =1
                    end if
                end do
            end do
            if(error .eq. 1) stop
        end if

        ! Updating WC dervatives by fluxes from BS
        do i = 1,nstate
            do k=1,nkn
            
              rates_fast(k,i,lmax_fast(k))  = rates_fast(k,i,lmax_fast(k)) + FLUXES_TO_ALUKAS_fast(k,i)/depth_fast(k,lmax_fast(k))
              !if(STRANGERS(rates_fast(k,i,lmax_fast(k))) .eq. 1) print *,'After ','k=',k,'i=',i, 'rates=',rates_fast(k,i,lmax_fast(k))
              dg(lmax_fast(k),k,i,26) = FLUXES_TO_ALUKAS_fast(k,i)
            
            end do
        end do

        ! Updating

        ! Checking for strange values in WC bottom layer derivatives
        if (debug_stranger) then
            do i = 1, nstate
               do k=1,nkn
                   if(STRANGERS(rates_fast(k,i,lmax_fast(k))) .eq. 1) then
                        print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                        print *, 'AQUABC after call to BS: '
                        print *,  'Variable ',i,'Cell ',k, 'Layer ',lmax_fast(k)
                        print *, 'WC bottom layer derivative is NaN'
                        print *, 'Derivative',rates_fast(k,i,lmax_fast(k))
                        print *, 'Flux from BS:', FLUXES_TO_ALUKAS_fast(k,i)
                        print *, 'Layer depth: ', depth_fast(k,lmax_fast(k))
                        print *, 'Fluxes/Depth', FLUXES_TO_ALUKAS_fast(k,i)/depth_fast(k,lmax_fast(k))
                        error =1
                   end if
               end do
            end do
            if(error .eq. 1) stop
        end if
    end if  !bsedim

    !ccccccccccccccccccccccccccc
    !     End of BS processing c
    !ccccccccccccccccccccccccccc




    !      Updating derivatives by loading. Not tested: fixme
    !cccccccccccccccccccccccccccccccccccccc
    !       call load0d(rates,loads,vol) !Atention, check cds updates(*Vol!) fixme
    
    !      Updating state variables by new derivatives
    !ccccccccccccccccccccccccccccccccccccccccccccccccc
    
    
    ! Integration
    do layer = 1,nlayers
        do k=1,nkn ! fixme for vectorization
            if(layer <= lmax_fast(k)) then
                state(1:nstate) = e_fast0(k,1:nstate,layer)
             rates(1:nstate) = rates_fast(k,1:nstate,layer)
             vol             = vol_fast    (k,layer)
             volold          = vol_old_fast(k,layer)
            
             call cur_euler(nstate, dt, vol, volold, rates, state)
            
             e_fast1(k,1:nstate,layer) = state(1:nstate)
            end if
        end do
    end do

      if (debug_stranger) then
          do i = 1,nlayers
              do j = 1, nstate
                  do k = 1, nkn
                      if(layer <= lmax_fast(k)) then
                          if(STRANGERS(e_fast1(k,j,i)) .eq. 1) then
                              !external node number
                              exn_node = ipv(k)
                              print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                              print *, 'AQUABC WC after integration:'
                              print *, 'Variable ',j,'Cell ',k, exn_node
                              print *, 'State variable is NaN'
                              print *, 'State=',e_fast1(k,j,i)
                              error =1
                          end if
                      end if ! layer
                  end do !k
              end do !j
          end do !i
          if(error .eq. 1) stop
      end if

      call derived_vars(nkn, nstate, nlayers, lmax_fast, e_fast1, dissolved_fractions_full, pH_full, &
                     wtemp_fast, sal_fast, noutput, wc_output)
                        
   return
   end ! subroutine aquabc



!**************************************************************
!************routines intended to use for sed trans************
!****************************not used *************************

!  subroutine resuspension(nkn,nlayers,lmax_fast, &
!                          waves_fast, depth_fast, rmass_fast)
!
!   ! Calculates resuspended mass of sediments
!   ! not finished. Tought to be used without sediment transport model
!
!   ! Inputs:
!   integer nkn, nlayers
!   integer lmax_fast(nkn)
!   real waves_fast(nkn,2)
!   real depth_fast(nkn,nlayers)
!
!   !Output:
!   real rmass_fast (nkn)
!
!   ! Internal vars:
!   real depth(nkn)
!   real, parameter :: g = 9.8
!   integer k
!
!
!   ! Calculate total depth
!   do k=1,nkn
!     depth(k) = sum(depth_fast(k,1:lmax_fast(k)))
!   end do
!
!
!
!  end subroutine resuspension
!
! !---------------------------
!
!  subroutine wavelength(nkn,wave_height,wave_length)
!  ! Calculates wave length
!  ! ! not finished. Tought to be used without wave model
!
!   ! Inputs:
!   integer nkn
!   real wave_height(nkn)
!
!   ! Output:
!   real wave_length
!
!
!  end subroutine wavelength




!***************************************************************
!***************************************************************
!***************************************************************


            subroutine layer_light(nkn,state,&
                            nstate,          &
                            par,             &
                            nconst,          &
                            layer,           &
                            lmax_fast,       &
                            depth_layer,     &
                            light_surface, alight_bottom)


!==========================================================================
! Calculates incident light after penetration of the given layer of the node
! Only nodes that exist for current layer are processed
! Petras Zemlys 2014.10.2
!
! Inputs:
!       state   -  state variable values for the layer and node
!       nstate  -  number of state variables
!       par     -  model constants
!       nconst  -  number of model constants
!       light_surface    -  incident light to the layer surface
!       layer   -  layer number
!       lmax    -  number of layers for the node
!       node    -  internal node nmber
! Outputs:
!       light_bottom   - incident light to the surface of layer
! Uses:
!       subroutine light_bottom
! Note:
!       Constants array is not used more but still
!       is not removed from parameters list. Fixme
!=============================================================

      use para_aqua
      use aquabc_pel_state_var_indexes

      implicit none

      integer nkn

      real state(nkn,nstate)
      real par (nconst)

      real light_surface(nkn)  ! incident light to the layer surface
      real alight_bottom (nkn)  ! incident light to the layer bottom

      real DIA_C              (nkn)
      real CYN_C              (nkn)
      real OPA_C              (nkn)
      real FIX_CYN_C          (nkn)

      double precision DIA_C_TO_CHLA
      double precision CYN_C_TO_CHLA
      double precision FIX_CYN_C_TO_CHLA
      double precision OPA_C_TO_CHLA

      real DIA_C_TO_CHLA_fast      (nkn)
      real CYN_C_TO_CHLA_fast      (nkn)
      real FIX_CYN_C_TO_CHLA_fast  (nkn)
      real OPA_C_TO_CHLA_fast      (nkn)
      real CHLA                    (nkn)


      real depth_layer  (nkn)
      !real vol, area  !just for the call dvanode, not used in calculations

      integer lmax_fast(nkn)! number of layers  for the nodes
      integer nstate, nconst, nnode, layer !,   mode, l

            ! Ratios C to Chla
            call para_get_value('DIA_C_TO_CHLA'    ,DIA_C_TO_CHLA)
            call para_get_value('CYN_C_TO_CHLA'    ,CYN_C_TO_CHLA)
            call para_get_value('FIX_CYN_C_TO_CHLA',FIX_CYN_C_TO_CHLA)
            call para_get_value('OPA_C_TO_CHLA'    ,OPA_C_TO_CHLA)

! Following assignments  depend on state variables  numbering!
! Change them if numbering is changed
            where (lmax_fast(1:nkn) >= layer) !
          ! Phytoplankton concentrations mgC/l:
             DIA_C     =  state(:, DIA_C_INDEX)
             CYN_C     =  state(:, CYN_C_INDEX)
             OPA_C     =  state(:, OPA_C_INDEX)
             FIX_CYN_C =  state(:, FIX_CYN_C_INDEX)

             DIA_C_TO_CHLA_fast(1:nkn)     = DIA_C_TO_CHLA
             CYN_C_TO_CHLA_fast(1:nkn)     = CYN_C_TO_CHLA
             FIX_CYN_C_TO_CHLA_fast(1:nkn) = FIX_CYN_C_TO_CHLA
             OPA_C_TO_CHLA_fast(1:nkn)     = OPA_C_TO_CHLA

             ! Total chorophyl a
            CHLA = ((DIA_C / DIA_C_TO_CHLA_fast) +     &
                   (CYN_C / CYN_C_TO_CHLA_fast)  +     &
                   (OPA_C / OPA_C_TO_CHLA_fast)  +     &
                   (FIX_CYN_C /FIX_CYN_C_TO_CHLA_fast )) * 1.0D3
            end where

           call light_bottom(nkn,lmax_fast,layer, depth_layer,CHLA, &
                             light_surface, alight_bottom)

      end  ! routine layer_light
!***********************************************************
!***********************************************************
!***********************************************************

      subroutine light_bottom(nkn,lmax,layer, DEPTH,CHLA,LIGHT,ALIGHT_bottom)

      ! Calculates light intensity on the bootom of the
      ! layer from light intensity on the surface of the layer and layer depth
      ! Takes  into account only Chl a in this version
      ! Calculates only for nodes that exist for current layer

          use para_aqua

          IMPLICIT NONE

          integer nkn      !number of nodes
          integer layer    !layer number
          integer lmax(nkn)!number of layers for the node

          real DEPTH (nkn) ! depth of layer surface
          real CHLA  (nkn) ! Chl a concentration in the overlaying layer

          real KE_fast   (nkn)
          real XKC_fast  (nkn)
          double precision KE     ! Background light extinction coefficient
          double precision XKC    ! extinction per chlorophyl unit

          real LIGHT (nkn) ! incident light
          real EXTIN (nkn) ! Total extinction coefficient

          real ALIGHT_bottom (nkn)

          ALIGHT_bottom(:) = 0.

          call para_get_value('KE', KE)
          KE_fast = KE
          != 0.8   ! later come as input parmeters background light extinction !fixme
                      ! Should be compatible with value in aquabc_II
          call para_get_value('XKC', XKC)
          XKC_fast = XKC
          != 0.010 ! extinction per chlorophyl unit                           !fixme
                      ! should be compatible with the value in subroutine LIM_LIGHT

!          EXTIN = KE + (8.8D-3 * CHLA) +
!     *            (5.4D-2 * (CHLA ** (2.0D0 / 3.0D0)))
         where( layer <= lmax)

          EXTIN = KE_fast + (XKC_fast * CHLA) ! Extinction coefficient
          ALIGHT_bottom = LIGHT * EXP(-EXTIN * DEPTH)

         end where

      END

!********************************************************************
!********************************************************************

!**************************************************************************
!**************************************************************************
        subroutine SETTLING_TO_THE_NEXT_LAYER            &
            (STATE_VARIABLES    , NUM_VARS   ,           &
             SETTLING_VELOCITIES, DISSOLVED_FRACTIONS,   &
             CELLNO, LAYER, PSTIME,                      &
             SETTLING_RATES)

        implicit none

        !Dummy arguments
        integer, intent(in) :: NUM_VARS, CELLNO, LAYER

      double precision, intent(in)    :: STATE_VARIABLES    (NUM_VARS)
      double precision, intent(in)    :: SETTLING_VELOCITIES(NUM_VARS)
      double precision, intent(in)    :: DISSOLVED_FRACTIONS(NUM_VARS)
      double precision, intent(inout) :: SETTLING_RATES     (NUM_VARS)
      double precision, intent(in)    :: PSTIME

      integer i

        do i = 1, NUM_VARS
          SETTLING_RATES(i) = STATE_VARIABLES(i) *  &
            (1.0D+0 - DISSOLVED_FRACTIONS(i))* SETTLING_VELOCITIES(i)
        end do

      end subroutine SETTLING_TO_THE_NEXT_LAYER

!**************************************************************************
!**************************************************************************
       subroutine cur_euler(nstate,dt,vol,volold,cds,c)

! new c is computed and replaced by new values


       implicit none

       integer nstate	!number of state variables
       real dt			!time step [day]
       real vol		    !volume [m**3]
       real c(nstate)	!state variable [mg/L] == [g/m**3]
       real cds(nstate)	!source term [g/day]
       real massnew
       integer i
       real volold, volnew	!volume [m**3]
       real mass		    !mass [g]
       real mder		    !derivative [g/day]

       !volold = vol
       volnew = vol

       do i=1,nstate
         mass    = c(i) * volold
         mder    = cds(i) * volold + c(i)*(vol-volold)
         massnew = mass + dt*mder
         c(i)    = massnew/volnew
!       c(i) = c(i) + dt*cds(i)
         if(c(i).lt.0.) c(i) = 1.e-20
       end do

       end


!********************************************************************
!********************************************************************
! Cheks for NaN and Inf
! Input is single precision!
! Scalar version

      INTEGER FUNCTION STRANGERS(VALUE)

      use mod_debug

      implicit none

      REAL VALUE, BIGNUMBER, RATIO, LLIMIT, ULIMIT
      
      logical VALUE_NaN
      logical VALUE_Inf
      logical VALUE_strange
      logical VALUE_out

      !print *, 'VALUE=', VALUE

      BIGNUMBER=1.0E30
      LLIMIT = -10000.
      ULIMIT =  10000. !Values for WC
      
      STRANGERS=0
      
      VALUE_out     = (VALUE < LLIMIT) .or. (VALUE > ULIMIT)
      VALUE_NAN     = is_nan(VALUE)
      VALUE_Inf     = abs(VALUE) > BIGNUMBER
      VALUE_strange = VALUE_NAN .or. VALUE_Inf .or. VALUE_out
      
      if(VALUE_strange) STRANGERS = 1      
!       print *,value
            
      return
      end

!*********************************************************
!*********************************************************
!*********************************************************

subroutine derived_vars(nkn,nstate, nlayers, lmax_fast, STATE_VARIABLES,  &
                        dissolved_fractions, pH, temp, sal, noutput, wc_outputs)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine for calculation of derived variables (calculated from state variables
! and model) and storage them to the array for final output to files
!
! Array WC_OUTPUTS is used for the output and gets (done in aquabc) also values of
! state variables as first nstate variables. Therefore is assumed that derived
! variables are placed after state variables.
!
! This routine is called in AQUABC after integration and derived
! variables values correspond to the same time momment, except variables
! passed from pelagic model (pH) that values correspond to the time one step earlier  !
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use AQUABC_II_GLOBAL
  use para_aqua
  use aquabc_pel_state_var_indexes

  implicit none
  include 'param.h'

! Inputs:
  integer         , intent(in) :: nkn, nstate, nlayers
  integer         , intent(in) :: lmax_fast  (nkn)                      ! max lev. numbers
  real            , intent(in) :: STATE_VARIABLES(nkn,nstate,nlayers)
  double precision, intent(in) :: dissolved_fractions(nkn,nstate,nlayers)
  real            , intent(in) :: pH  (nkn, nlayers)
  real            , intent(in) :: temp(nkn, nlayers)
  real            , intent(in) :: sal (nkn, nlayers)
  integer         , intent(in) :: noutput

! Output:
  real,  intent(out) :: WC_OUTPUTS(nlayers,nkn,noutput)  ! Array for state vars + auxilaries
                                                         ! derived vars are stored as last variables
                                                         ! state vars are stored as first vars later
                                                         ! before final output to files
! Auxilaries:
  integer i, k, l, lmax, ext_node
  double precision   :: WC_OUTPUTS_temp(nkn,noutput,nlayers) ! Array for temporary storage

! State variables
    double precision :: NH4_N             (nkn, nlayers)
    double precision :: NO3_N             (nkn, nlayers)
    double precision :: PO4_P             (nkn, nlayers)
    double precision :: DISS_OXYGEN       (nkn, nlayers)
    double precision :: CHEM_AUT_BAC_C    (nkn, nlayers)
    double precision :: AER_HET_BAC_C     (nkn, nlayers)
    double precision :: FAC_AN_HET_BAC_C  (nkn, nlayers)
    double precision :: DIA_C             (nkn, nlayers)
    double precision :: ZOO_C             (nkn, nlayers)
    !double precision :: ZOO_N             (nkn, nlayers)
    !double precision :: ZOO_P             (nkn, nlayers)
    double precision :: DET_PART_ORG_C    (nkn, nlayers)
    double precision :: DET_PART_ORG_N    (nkn, nlayers)
    double precision :: DET_PART_ORG_P    (nkn, nlayers)
    double precision :: DISS_ORG_C        (nkn, nlayers)
    double precision :: DISS_ORG_N        (nkn, nlayers)
    double precision :: DISS_ORG_P        (nkn, nlayers)
    double precision :: CYN_C             (nkn, nlayers)
    double precision :: OPA_C             (nkn, nlayers)
    !double precision :: DISS_Si           (nkn, nlayers)
    !double precision :: PART_Si           (nkn, nlayers)
    double precision :: FIX_CYN_C         (nkn, nlayers)
    double precision :: FEII              (nkn, nlayers)
    double precision :: FEIII             (nkn, nlayers)
    double precision :: MNII              (nkn, nlayers)
    double precision :: MNIV              (nkn, nlayers)
    double precision :: CA                (nkn, nlayers)
    double precision :: MG                (nkn, nlayers)

    ! Passed derived variables


    ! New state variables added 22 September 2014
    !double precision(kind=DBL_PREC) :: INORG_C        (nkn, nlayers)   !Inorganic carbon
    !double precision(kind=DBL_PREC) :: TOT_ALK        (nkn, nlayers)   !Total alkalinity

! Constants

  double precision ::              CHEM_AUT_BAC_N_TO_C
  double precision ::              CHEM_AUT_BAC_P_TO_C

  double precision ::               AER_HET_BAC_N_TO_C
  double precision ::               AER_HET_BAC_P_TO_C

  double precision ::            FAC_AN_HET_BAC_N_TO_C
  double precision ::            FAC_AN_HET_BAC_P_TO_C

  double precision ::                       DIA_N_TO_C
  double precision ::                       DIA_P_TO_C

  double precision ::                       CYN_N_TO_C
  double precision ::                       CYN_P_TO_C

  double precision ::                   FIX_CYN_N_TO_C
  double precision ::                   FIX_CYN_P_TO_C

  double precision ::                       OPA_N_TO_C
  double precision ::                       OPA_P_TO_C

  double precision ::                       ZOO_N_TO_C
  double precision ::                       ZOO_P_TO_C
  
  double precision :: T ! auxiliary for temp                    
  double precision :: S ! auxiliary for sal
  
  double precision ::  DO_SATURATION !function


  WC_OUTPUTS_temp(:,:,:) = 0.

! Only necessary state variables are assigned here
    NH4_N           (:,:)      = STATE_VARIABLES(:,NH4_N_INDEX          ,:)      ! AMMONIUM NITROGEN
    NO3_N           (:,:)      = STATE_VARIABLES(:,NO3_N_INDEX          ,:)      ! NITRATE NITROGEN
    PO4_P           (:,:)      = STATE_VARIABLES(:,PO4_P_INDEX          ,:)      ! ORTHOPHOSPHATE PHOSPHORUS
    DISS_OXYGEN     (:,:)      = STATE_VARIABLES(:,DISS_OXYGEN_INDEX    ,:)      ! DISSOLVED OXYGEN
    DIA_C           (:,:)      = STATE_VARIABLES(:,DIA_C_INDEX          ,:)      ! DIATOMS CARBON
    ZOO_C           (:,:)      = STATE_VARIABLES(:,ZOO_C_INDEX          ,:)      ! ZOOPLANKTON CARBON
   !ZOO_N           (:,:)      = STATE_VARIABLES(:,ZOO_N_INDEX          ,:)     ! ZOOPLANKTON CARBON
   !ZOO_P           (:,:)      = STATE_VARIABLES(:,ZOO_P_INDEX          ,:)     ! ZOOPLANKTON PHOSPHORUS
    DET_PART_ORG_C  (:,:)      = STATE_VARIABLES(:,DET_PART_ORG_C_INDEX ,:)     ! DETRITUS PARTICULATE ORG. CARBON
    DET_PART_ORG_N  (:,:)      = STATE_VARIABLES(:,DET_PART_ORG_N_INDEX ,:)     ! DETRITUS PARTICULATE ORG. NITROGEN
    DET_PART_ORG_P  (:,:)      = STATE_VARIABLES(:,DET_PART_ORG_P_INDEX ,:)     ! DETRITUS PARTICULATE ORG. PHOSPHORUS
    DISS_ORG_C      (:,:)      = STATE_VARIABLES(:,DISS_ORG_C_INDEX     ,:)     ! DISSOLVED ORGANIC CARBON
    DISS_ORG_N      (:,:)      = STATE_VARIABLES(:,DISS_ORG_N_INDEX     ,:)     ! DISSOLVED ORGANIC NITROGEN
    DISS_ORG_P      (:,:)      = STATE_VARIABLES(:,DISS_ORG_P_INDEX     ,:)     ! DISSOLVED ORGANIC PHOSPHORUS
    CYN_C           (:,:)      = STATE_VARIABLES(:,CYN_C_INDEX          ,:)     ! NON FIXING CYANOBACTERIA CARBON
    OPA_C           (:,:)      = STATE_VARIABLES(:,OPA_C_INDEX          ,:)     ! OTHER PHYTOPLANKTON CARBON
   !DISS_Si         (:,:)      = STATE_VARIABLES(:,DISS_Si_INDEX        ,:)     ! DISSOLOVED SILICA
   !PART_Si         (:,:)      = STATE_VARIABLES(:,PART_Si_INDEX        ,:)     ! PARTICULATE SILICA
    FIX_CYN_C       (:,:)      = STATE_VARIABLES(:,FIX_CYN_C_INDEX      ,:)     ! FIXING CYANOBACTERIA CARBON
   !INORG_C         (:,:)      = STATE_VARIABLES(:,INORG_C_INDEX        ,:)     ! INORG CARBON CARBON
   !TOT_ALK         (:,:)      = STATE_VARIABLES(:,TOT_ALK_INDEX        ,:)     ! TOTAL ALKALNITY
    FEII             (:,:)     = STATE_VARIABLES(:,FE_II_INDEX          ,:)     ! Iron
    FEIII            (:,:)     = STATE_VARIABLES(:,FE_III_INDEX         ,:)     !
    MNII             (:,:)     = STATE_VARIABLES(:,MN_II_INDEX          ,:)     ! Manganese
    MNIV             (:,:)     = STATE_VARIABLES(:,MN_IV_INDEX          ,:)     !
    CA               (:,:)     = STATE_VARIABLES(:,CA_INDEX             ,:)     ! Calcium
    MG               (:,:)     = STATE_VARIABLES(:,MG_INDEX             ,:)     ! Magnesium



     ! Only necessary parameters are obtained

!       call para_get_value('CHEM_AUT_BAC_N_TO_C'             ,              CHEM_AUT_BAC_N_TO_C) ! 16 Chemoautotrophic bacteria Nitrogen to Carbon ratio
!       call para_get_value('CHEM_AUT_BAC_P_TO_C'             ,              CHEM_AUT_BAC_P_TO_C) ! 17 Chemoautotrophic bacteria Phosphorus to Carbon ratio
! 
!       call para_get_value('AER_HET_BAC_N_TO_C'              ,               AER_HET_BAC_N_TO_C) ! 42 Aerobic heterotrophic bacteria Nitrogen to Carbon ratio
!       call para_get_value('AER_HET_BAC_P_TO_C'              ,               AER_HET_BAC_P_TO_C) ! 43 Aerobic heterotrophic bacteria Phosphorus to Carbon ratio
! 
!       call para_get_value('FAC_AN_HET_BAC_N_TO_C'           ,            FAC_AN_HET_BAC_N_TO_C) ! 60 not used! Facultative anaerobic heterotrophic bacteria Nitrogen to Carbon ratio
!       call para_get_value('FAC_AN_HET_BAC_P_TO_C'           ,            FAC_AN_HET_BAC_P_TO_C) ! 61 not used! Facultative anaerobic heterotrophic bacteria Phosphorus to Carbon ratio

      call para_get_value('DIA_N_TO_C'                      ,                       DIA_N_TO_C) ! 83 Diatoms Nitrogen to Carbon ratio
      call para_get_value('DIA_P_TO_C'                      ,                       DIA_P_TO_C) ! 84 Diatoms Phosphorus to Carbon ratio

      call para_get_value('CYN_N_TO_C'                      ,                       CYN_N_TO_C) !106 Non-fixing cyanobacteria Nitrogen to Carbon ratio ,was 0.1
      call para_get_value('CYN_P_TO_C'                      ,                       CYN_P_TO_C) !107 Non-fixing cyanobacteria Phosphorus to Carbon ratio

      call para_get_value('FIX_CYN_N_TO_C'                  ,                   FIX_CYN_N_TO_C) !128 Fixing cyanobacteria Nitrogen to Carbon ratio
      call para_get_value('FIX_CYN_P_TO_C'                  ,                   FIX_CYN_P_TO_C) !129 Fixing cyanobacteria Phosphorus to Carbon ratio

      call para_get_value('OPA_N_TO_C'                      ,                       OPA_N_TO_C) !152 OtherPhyto Nitrogen to Carbon ratio
      call para_get_value('OPA_P_TO_C'                      ,                       OPA_P_TO_C) !153 OtherPhyto Phosphorus to Carbon ratio

      call para_get_value('ZOO_N_TO_C'                      ,                       ZOO_N_TO_C) !196 Zooplankton Nitrogen to Carbon ratio
      call para_get_value('ZOO_P_TO_C'                      ,                       ZOO_P_TO_C) !197 Zooplankton Phosphorus to Carbon ratio






!         Checking NOUTPUT compatibility with the code
          if (nstate+24 .gt. noutput) then
          print *, 'DERIVED_VARS: noutput is to small'
          print *, nstate+24, noutput
          stop
         end if

 ! Calculating and sending Total N, Total P and Total IN to the output array
 ! Note: for elements layer > lmax_fast(k) zeros are assumed, so division operation
 ! can lead to exceptions.

 !       Total N
         WC_OUTPUTS_temp(1:nkn,nstate+1,1:nlayers) =    NH4_N    +   NO3_N + &    !Dissolved inorganic
!                                  CHEM_AUT_BAC_C     * CHEM_AUT_BAC_N_TO_C   + &   !Bacteria
!                                  AER_HET_BAC_C      * AER_HET_BAC_N_TO_C    + &
!                                  FAC_AN_HET_BAC_C   * FAC_AN_HET_BAC_N_TO_C + &
                                 DIA_C * DIA_N_TO_C + &                           !Phyto
                                 CYN_C * CYN_N_TO_C + &
                                 OPA_C * OPA_N_TO_C + &
                                 FIX_CYN_C      * FIX_CYN_N_TO_C + &
                                 DET_PART_ORG_N + DISS_ORG_N +   &                 !Part + Diss
                                 ZOO_C          * ZOO_N_TO_C                       !Zoo



         !Total P
         WC_OUTPUTS_temp(1:nkn,nstate+2,1:nlayers) = PO4_P + &                  !Phosphates
!                                 CHEM_AUT_BAC_C   * CHEM_AUT_BAC_P_TO_C  + &     !Bacteria
!                                 AER_HET_BAC_C    * AER_HET_BAC_P_TO_C   + &
!                                 FAC_AN_HET_BAC_C * FAC_AN_HET_BAC_P_TO_C+ &
                                DIA_C * DIA_P_TO_C  + &                        !Phyto
                                CYN_C * CYN_P_TO_C + &
                                OPA_C * OPA_P_TO_C + &
                                FIX_CYN_C      * CYN_P_TO_C + &
                                DET_PART_ORG_P + DISS_ORG_P + &                !Part + Diss
                                ZOO_C          * ZOO_P_TO_C                    !Zoo


         ! Total IN
         WC_OUTPUTS_temp(1:nkn,nstate+3,1:nlayers) = NH4_N + NO3_N
         ! pH

         !print *, pH(1:3, 1)
         !stop
         WC_OUTPUTS_temp(1:nkn,nstate+4,1:nlayers) = pH
         ! Total phytoplankton
         WC_OUTPUTS_temp(1:nkn,nstate+5,1:nlayers) = DIA_C + OPA_C + CYN_C + FIX_CYN_C
         ! Total cyano bacteria
         WC_OUTPUTS_temp(1:nkn,nstate+6,1:nlayers) = CYN_C + FIX_CYN_C
         
         WC_OUTPUTS_temp(1:nkn,nstate+ 7,1:nlayers) =  dissolved_fractions(:,PO4_P_INDEX  ,:) * PO4_P
         WC_OUTPUTS_temp(1:nkn,nstate+ 9,1:nlayers) =  dissolved_fractions(:,FE_II_INDEX  ,:) * FEII
         WC_OUTPUTS_temp(1:nkn,nstate+11,1:nlayers) =  dissolved_fractions(:,FE_III_INDEX ,:) * FEIII
         WC_OUTPUTS_temp(1:nkn,nstate+13,1:nlayers) =  dissolved_fractions(:,MN_II_INDEX  ,:) * MNII
         WC_OUTPUTS_temp(1:nkn,nstate+15,1:nlayers) =  dissolved_fractions(:,MN_IV_INDEX  ,:) * MNIV
         WC_OUTPUTS_temp(1:nkn,nstate+17,1:nlayers) =  dissolved_fractions(:,CA_INDEX     ,:) * CA
         WC_OUTPUTS_temp(1:nkn,nstate+19,1:nlayers) =  dissolved_fractions(:,MG_INDEX     ,:) * MG
         
         ! Total dissolved iron
         WC_OUTPUTS_temp(1:nkn,nstate+21,1:nlayers) = & 
            WC_OUTPUTS_temp(1:nkn,nstate+ 9,1:nlayers) + WC_OUTPUTS_temp(1:nkn,nstate+11,1:nlayers)
         
         WC_OUTPUTS_temp(1:nkn,nstate+ 8,1:nlayers) =  (1 - dissolved_fractions(:,PO4_P_INDEX  ,:)) * PO4_P
         WC_OUTPUTS_temp(1:nkn,nstate+10,1:nlayers) =  (1 - dissolved_fractions(:,FE_II_INDEX  ,:)) * FEII
         WC_OUTPUTS_temp(1:nkn,nstate+12,1:nlayers) =  (1 - dissolved_fractions(:,FE_III_INDEX ,:)) * FEIII
         WC_OUTPUTS_temp(1:nkn,nstate+14,1:nlayers) =  (1 - dissolved_fractions(:,MN_II_INDEX  ,:)) * MNII
         WC_OUTPUTS_temp(1:nkn,nstate+16,1:nlayers) =  (1 - dissolved_fractions(:,MN_IV_INDEX  ,:)) * MNIV
         WC_OUTPUTS_temp(1:nkn,nstate+18,1:nlayers) =  (1 - dissolved_fractions(:,CA_INDEX     ,:)) * CA
         WC_OUTPUTS_temp(1:nkn,nstate+20,1:nlayers) =  (1 - dissolved_fractions(:,MG_INDEX     ,:)) * MG
         
         ! Total particulate iron
         WC_OUTPUTS_temp(1:nkn,nstate+22,1:nlayers) = &
            WC_OUTPUTS_temp(1:nkn,nstate+ 10,1:nlayers) + WC_OUTPUTS_temp(1:nkn,nstate+12,1:nlayers)
            
         ! Total iron   
         WC_OUTPUTS_temp(1:nkn,nstate+23,1:nlayers) = &
            WC_OUTPUTS_temp(1:nkn,nstate+ 21,1:nlayers) + WC_OUTPUTS_temp(1:nkn,nstate+22,1:nlayers)

         do k = 1,nkn
             lmax = lmax_fast(k)
             
              do l = 1,lmax               
               T = temp (k,l)
               S = sal  (k,l)
               WC_OUTPUTS_temp(k,nstate+24,l) = (DISS_OXYGEN(k,l)/DO_SATURATION(T, S, 0.D0))*100.D0
              end do
             
         end do

         ! Reshape and assign temporrary array to final
         !WC_OUTPUTS = reshape(WC_OUTPUTS_temp,shape=(/nlayers,nkn,noutput/), order=(/3,1,2/))

         do k = 1,nkn
             lmax = lmax_fast(k)
             do i = nstate+1,noutput
              do l = 1,lmax
               WC_OUTPUTS(l,k,i) = WC_OUTPUTS_temp(k,i,l)
              end do
             end do
         end do

   end subroutine derived_vars
!************************************************************************
!************************************************************************

!********************************************************************

	subroutine load0d(cds,loads,vol)

! integrate loadings

	implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )

	real cds(nstate)	!source term [g/day]
	real loads(nstate)	!loading for c [g/(m**3 day)]
	real vol		!volume of box [m**3]

	integer i
        
	do i=1,nstate
	  cds(i) = cds(i) + vol * loads(i)
	end do

	end

!********************************************************************

































