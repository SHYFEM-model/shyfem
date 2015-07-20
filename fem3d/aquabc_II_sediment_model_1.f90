!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Main routines for bottom sediment model No1 
! CONTENT:
!  subroutine SEDIMENT_MODEL_1   -  main routine
!  FUNCTION STRANGER(VALUE)      -  checks for NaNs (input type is double precision)
!  subroutine FLX_ALUKAS_II_TO_SED_MOD_1
!  FUNCTION SED_MOD_1_ALUKAS_MOLDI_C - molecular diffusion coefficients
!  subroutine SED_MOD_1_CVISC
!  subroutine FLX_SED_MOD_1_TO_ALUKAS_II 
!
! Produced by Ali Erturk 2010 June
!
! Enhancement by Ali Erturk 13th of February 2011
!           - Added advection as a new transport process
!
! Enhancement by Ali Erturk and Petras Zemlys 26th of September 2014
!           - Expanded the model by adding inorganic carbon, alkality and
!             salinity as new state variables and pH as a derived variable
!
! Change by Ali Erturk 2nd of February 2015
!           - Carbon related state variables are subroutinized
!************************************************************************

subroutine SEDIMENT_MODEL_1 &
           (nkn,INIT_SED_STATE_VARS  , SED_DEPTHS , SED_POROSITIES,  &
            SED_DENSITIES        , PART_MIXING_COEFFS         ,  &
            SED_DIFFUSIONS       , SURF_MIXLEN, SED_BURRIALS  ,  &
            SURF_WATER_CONCS     , SED_TEMPS                  ,  &
            NUM_SED_VARS         , NUM_SED_LAYERS             ,  &
            SED_MODEL_CONSTANTS  , NUM_SED_CONSTS             ,  &
            SED_DRIVING_FUNCTIONS, NUM_SED_DRIV               ,  & ! not used yet
            FLUXES_TO_SEDIMENTS  , NUM_FLUXES_TO_SEDIMENTS    ,  &
            NUM_FLUX_RECEIVING_SED_LAYERS, ADVECTIVE_VELOCITY ,  &
            PSTIME, TIME_STEP                                 ,  &
            FINAL_SED_STATE_VARS ,                               &
            FLUXES_FROM_SEDIMENTS, NUM_FLUXES_FROM_SEDIMENTS  ,  &
            PROCESSES_sed        , NDIAGVAR_sed ,                &
            SED_OUTPUTS          , NUM_SED_OUTPUTS)

    use CO2SYS_CDIAC
    use AQUABC_II_GLOBAL
    use para_aqua

    implicit none
    
    include 'param.h'
    
    integer nkn ! number of reactors (nodes) 
    integer ipv(nkndim)	!external node numbers
    common  /ipv/ipv    

    !ARGUMENTS RELATED TO ARAY SIZES

    !NUMBER OF CONSTANTS FOR CONTROL
    integer NUM_SED_MOD_CONSTS
    parameter(NUM_SED_MOD_CONSTS = 42)

    integer NUM_SED_VARS
    integer NUM_SED_LAYERS
    integer NUM_SED_CONSTS
    integer NUM_SED_DRIV
    integer NUM_FLUXES_TO_SEDIMENTS
    integer NUM_FLUXES_FROM_SEDIMENTS
    integer NDIAGVAR_sed
    integer NUM_FLUX_RECEIVING_SED_LAYERS 
    integer NUM_SED_OUTPUTS
    !PARAMETER(NUM_SED_OUTPUTS = 14)
       
    !INPUT ARGUMENTS
    double precision INIT_SED_STATE_VARS  (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision SED_DEPTHS           (nkn,NUM_SED_LAYERS)
    double precision SED_DENSITIES        (nkn,NUM_SED_LAYERS)
    double precision SED_POROSITIES       (nkn,NUM_SED_LAYERS)
    double precision PART_MIXING_COEFFS   (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision SED_DIFFUSIONS       (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision SURF_MIXLEN          
    double precision SED_BURRIALS         (nkn,NUM_SED_LAYERS)
    double precision SURF_WATER_CONCS     (nkn,NUM_SED_VARS)
    double precision SED_TEMPS            (nkn,NUM_SED_LAYERS)
    double precision SED_MODEL_CONSTANTS  (NUM_SED_CONSTS)
    double precision PROCESSES_sed        (nkn,NUM_SED_LAYERS, NUM_SED_VARS, NDIAGVAR_sed)
    double precision SED_DRIVING_FUNCTIONS(NUM_SED_LAYERS, NUM_SED_DRIV)
    double precision FLUXES_TO_SEDIMENTS  (nkn,NUM_FLUXES_TO_SEDIMENTS)

    !ADVECTION      
    double precision ADVECTIVE_VELOCITY
    double precision ADV_ENTERING_CONC(nkn)
    
    integer CELLNO
    integer LAYER
    double precision PSTIME
    double precision TIME_STEP
    
    !OUTPUT ARGUMENTS
    double precision FINAL_SED_STATE_VARS (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision FLUXES_FROM_SEDIMENTS(nkn,NUM_FLUXES_FROM_SEDIMENTS)
    double precision SED_OUTPUTS          (nkn,NUM_SED_LAYERS, NUM_SED_OUTPUTS)
    
    !COUNTERS
    integer I, J,k,l

    !SEDIMENT STATE VARIABLES
    double precision SED_NH4N(nkn,NUM_SED_LAYERS)
    double precision SED_NO3N(nkn,NUM_SED_LAYERS)
    double precision SED_DON (nkn,NUM_SED_LAYERS)
    double precision SED_PON (nkn,NUM_SED_LAYERS)
    double precision SED_PO4P(nkn,NUM_SED_LAYERS)
    double precision SED_DOP (nkn,NUM_SED_LAYERS)
    double precision SED_POP (nkn,NUM_SED_LAYERS)
    double precision SED_DOXY(nkn,NUM_SED_LAYERS)
    double precision SED_DOC (nkn,NUM_SED_LAYERS)
    double precision SED_POC (nkn,NUM_SED_LAYERS)
    double precision SED_DSi (nkn,NUM_SED_LAYERS)
    double precision SED_PSi (nkn,NUM_SED_LAYERS)
    
    ! New state variables added 27 September 2014
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: INORG_C   !Inorganic carbon
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: TOT_ALK   !Total alkalinity
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: SALT      !Total alkalinity
    ! End of new state variables added 27 September 2014
    
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: pH
       
    !SEDIMENT MODEL COEFFICIENTS
    !Dissolution of particulate organic carbon
    double precision K_OXIC_DISS_POC   
    double precision K_ANOXIC_DISS_POC 
    double precision THETA_DISS_POC    
    double precision KHS_DISS_POC

    !Dissolution of particulate organic nitrogen
    double precision K_OXIC_DISS_PON   
    double precision K_ANOXIC_DISS_PON 
    double precision THETA_DISS_PON    
    double precision KHS_DISS_PON 
    
    !Dissolution of particulate organic phosphorus
    double precision K_OXIC_DISS_POP   
    double precision K_ANOXIC_DISS_POP 
    double precision THETA_DISS_POP    
    double precision KHS_DISS_POP
    
    !Dissolution of particulate silicon
    double precision K_OXIC_DISS_PSi   
    double precision K_ANOXIC_DISS_PSi 
    double precision THETA_DISS_PSi    
    double precision KHS_DISS_PSi
    
    !Mineralization of dissolved organic carbon
    double precision K_OXIC_MINER_DOC
    double precision K_ANOXIC_MINER_DOC
    double precision THETA_MINER_DOC
    double precision KHS_MINER_DOC 
    double precision O_TO_C
    
    !Mineralization of dissolved organic nitrogen
    double precision K_OXIC_MINER_DON
    double precision K_ANOXIC_MINER_DON
    double precision THETA_MINER_DON
    double precision KHS_MINER_DON
    
    !Mineralization of dissolved organic phosphorus
    double precision K_OXIC_MINER_DOP
    double precision K_ANOXIC_MINER_DOP
    double precision THETA_MINER_DOP
    double precision KHS_MINER_DOP
    
    !Nitrification
    double precision K_NITR
    double precision THETA_NITR
    double precision KHS_NITR_NH4N
    double precision KHS_NITR_DOXY
    
    !Denitrification
    double precision K_DENITR
    double precision THETA_DENITR
    double precision KHS_DENITR_NO3N
    double precision KHS_DENITR_DOC
    double precision KHS_DENITR_DOXY
    double precision DENITR_YIELD
    
    !Anoxia     
    double precision DOXY_AT_ANOXIA
    
    !Solid partition      
    double precision SOLID_PART_COEFF_NH4
    double precision SOLID_PART_COEFF_PO4
    
    
    !SEDIMENT KINETIC PROCESS RATES
    double precision, dimension(nkn,NUM_SED_VARS) :: R_DISS_POC
    double precision, dimension(nkn,NUM_SED_VARS) :: R_DISS_PON
    double precision, dimension(nkn,NUM_SED_VARS) :: R_DISS_POP
    double precision, dimension(nkn,NUM_SED_VARS) :: R_MINER_DOC
    double precision, dimension(nkn,NUM_SED_VARS) :: R_MINER_DON
    double precision, dimension(nkn,NUM_SED_VARS) :: R_MINER_DOP
    double precision, dimension(nkn,NUM_SED_VARS) :: R_NITR
    double precision, dimension(nkn,NUM_SED_VARS) :: R_DENITR
    double precision, dimension(nkn,NUM_SED_VARS) :: R_DISS_PSi
    double precision, dimension(nkn,NUM_SED_VARS) :: DEOXYGENATION
    
    !SEDIMENT TRANSPORT PROCESS RATES
    double precision NEIGHBOUR_CONC     (nkn)
    double precision UPPER_CONC_GRADIENT(nkn)
    double precision SED_MIXLEN         (nkn)
    double precision SED_IN_ADVEC_RATES (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision SED_OUT_ADVEC_RATES(nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision SED_DIFFUSION_RATES(nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision PART_MIXING_RATES  (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision SED_BURRIAL_RATES  (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision UNIT_AREA_MASSES   (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    
    !DERIVS
    double precision DERIVS             (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision TRANSPORT_DERIVS   (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision ADVECTION_DERIVS   (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision DIFFUSION_DERIVS   (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision BURIAL_DERIVS      (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision PART_MIXING_DERIVS (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision SETTLING_DERIVS    (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision KINETIC_DERIVS     (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    
    double precision TEMP
    double precision SETTLING_AFFECTED_DEPTH(nkn)
    
    integer STRANGER  !Function checking for strange values 
    integer STRANGERSD 
    integer error     !Error indicator
    logical VALUE_strange(nkn) ! For NaN and Inf checking
    ! For strange values processing    
    real(kind = DBL_PREC),allocatable, dimension (:) :: STRANGERS     ! Strange values
    integer              ,allocatable, dimension (:) :: NODES_STRANGE ! node numbers with strange values
    integer :: index_strange (nkn)                                    ! array with values 1:nkn
    integer :: nstrange                                               ! number of nodes with strange values 
    logical debug_stranger
    
    double precision DIFF_CORRECTION_FACTOR(nkn)
    
    integer TIME_LOOP
    integer NUM_SUB_TIME_STEPS
    double precision INTERMED_RESULTS(nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision VOLUME_FRACTION(nkn)
    integer IN_WHICH_PHASE(NUM_SED_VARS)
    
    double precision SOLUTE_FRACTIONS(nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision WATER_DENSITY
    double precision SOLID_CONCS(nkn,NUM_SED_LAYERS)
    
    double precision SINK_INTENSITY, DIFF_SINK(nkn) ! for difussion boundary sink gradient definition
    
    INTEGER switch_kinetics
    INTEGER switch_burial
    INTEGER switch_partmixing
    INTEGER switch_diffusion
    INTEGER switch_advection
    INTEGER switch_settling
    double precision DIFF_DRAG

    
    !27 September 2014
    !New variables for DIC and ALK
    real(kind = DBL_PREC), allocatable, dimension (:) :: CO2SYS_PAR1
    real(kind = DBL_PREC), allocatable, dimension (:) :: CO2SYS_PAR2
    integer, allocatable, dimension (:) :: CO2SYS_PAR1TYPE
    integer, allocatable, dimension (:) :: CO2SYS_PAR2TYPE
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_SALT
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_TEMPIN
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_TEMPOUT
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_PRESIN
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_PRESOUT
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_SI
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_PO4
    integer, allocatable, dimension (:) :: CO2SYS_pHSCALEIN
    integer, allocatable, dimension (:) :: CO2SYS_K1K2CONSTANTS
    integer, allocatable, dimension (:) :: CO2SYS_KSO4CONSTANTS
    real(kind = DBL_PREC),allocatable, dimension (:,:) :: CO2SYS_OUT_DATA
    character(len=34), allocatable, dimension (:) :: CO2SYS_NICEHEADERS
    integer :: CO2SYS_NUM_SAMPLES
    integer :: RUN_CO2SYS
    parameter(RUN_CO2SYS = 1)
    integer :: CO2SYS_ntps
    
    !New variables for DIC and ALK
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: TOTAL_DIC_KINETIC_SOURCES
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: TOTAL_DIC_KINETIC_SINKS
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: T_A
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_K_H
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: K_H
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: POWER
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_CO2
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: CO2_SAT
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: K_A_CALC_CO2
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: CO2_ATM_EXHANGE
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: DIC_KINETIC_DERIVATIVE
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: ALPHA_0
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: ALPHA_1

    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: PKH
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: FRAC_NH3
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: FRAC_NH4
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_CHEM_AUT_BAC_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_AER_HET_BAC_INT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_FAC_AN_HET_BAC_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_DIA_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_CYN_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_OPA_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_FIX_CYN_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_AER_HET_BAC_N_OX
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_FAC_AN_HET_BAC_N_OX
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_ZOO_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_ABIOTIC_DON_MIN
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: ALK_GAINED_BY_AMMONIUM_GEN
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_DENITRIFICATION
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_DIA_GROWTH
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_CYN_GROWTH
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_OPA_GROWTH
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_NON_FIX_CYN_GROWTH
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_AER_HET_BAC_GROWTH
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: ALK_GAINED_BY_NITRATE_CONS
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_CHEM_AUT_BAC_GROWTH
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: ALK_LOST_BY_AMMONIUM_CONS
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_NITRIFICATION_NH4
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: N_NITRIFICATION_NH3
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: ALK_LOST_BY_NITRIFICATION
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: H_PLUS
    integer :: KP_OPTION
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: K_ONE_TIP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: K_TWO_TIP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: K_THREE_TIP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: FRACTION_DIVISOR_TIP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: ALPHA_H2PO4
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: ALPHA_HPO4
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: ALPHA_PO4
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: PHOSPHATE_EQ_CONSTANT
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: ALK_GAINED_BY_PHOSPHATE_CONS
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: ALK_LOST_BY_PHOSPHATE_GEN
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: ALK_KINETIC_DERIVATIVE
    
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_CHEM_AUT_BAC_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_AER_HET_BAC_INT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_FAC_AN_HET_BAC_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_DIA_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_CYN_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_OPA_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_FIX_CYN_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_AER_HET_BAC_N_OX
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_FAC_AN_HET_BAC_N_OX
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_ZOO_TOT_RESP
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_ABIOTIC_DON_MIN
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_DENITRIFICATION
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_DIA_GROWTH
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_CYN_GROWTH
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_OPA_GROWTH
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_NON_FIX_CYN_GROWTH
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_AER_HET_BAC_GROWTH
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_AER_HET_BAC_P_OX
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_FAC_AN_HET_BAC_P_OX
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: P_CHEM_AUT_BAC_GROWTH
    integer :: CONSIDER_ALKALNITY_DERIVATIVE
    integer :: CONSIDER_INORG_C_DERIVATIVE
    integer :: CONSIDER_CO2_REARATION
    
    !End of new variables for DIC and ALK

    
    CONSIDER_INORG_C_DERIVATIVE = 0	!avoid compiler warnings
    CONSIDER_ALKALNITY_DERIVATIVE = 0	!avoid compiler warnings

    
    !**************************************************************************
    !*                                                                        *
    !*                            EXECUTION PART                              *
    !*                                                                        *
    !**************************************************************************
    
    debug_stranger = .true.
    debug_stranger = .false. !True if check for strangers
     
    ! Multiplier for diffusion rate for the first layer (used reverse for negative)
      DIFF_DRAG = 6.D0
      
      
    !Initializations
    !SEDIMENT KINETIC PROCESS RATES
     R_DISS_POC    = 0.0D0
     R_DISS_PON    = 0.0D0
     R_DISS_POP    = 0.0D0
     R_MINER_DOC   = 0.0D0
     R_MINER_DON   = 0.0D0
     R_MINER_DOP   = 0.0D0
     R_NITR        = 0.0D0
     R_DENITR      = 0.0D0
     R_DISS_PSi    = 0.0D0
    
    !SEDIMENT TRANSPORT PROCESS RATES
     NEIGHBOUR_CONC     (:)      = 0.0D0
     UPPER_CONC_GRADIENT(:)      = 0.0D0
     SED_MIXLEN         (:)      = 0.0D0
     SED_IN_ADVEC_RATES (:,:, :) = 0.0D0
     SED_OUT_ADVEC_RATES(:,:, :) = 0.0D0
     SED_DIFFUSION_RATES(:,:, :) = 0.0D0
     PART_MIXING_RATES  (:,:, :) = 0.0D0
     SED_BURRIAL_RATES  (:,:, :) = 0.0D0
     UNIT_AREA_MASSES   (:,:, :) = 0.0D0

    !DERIVS
     DERIVS             (:,:, :) = 0.0D0
     TRANSPORT_DERIVS   (:,:, :) = 0.0D0
     ADVECTION_DERIVS   (:,:, :) = 0.0D0
     DIFFUSION_DERIVS   (:,:, :) = 0.0D0
     BURIAL_DERIVS      (:,:, :) = 0.0D0
     PART_MIXING_DERIVS (:,:, :) = 0.0D0
     SETTLING_DERIVS    (:,:, :) = 0.0D0
     KINETIC_DERIVS     (:,:, :) = 0.0D0
    
     TEMP = 0.0D0
     SETTLING_AFFECTED_DEPTH = 0.0D0

    
     DEOXYGENATION            = 0.0D0
     DIFF_CORRECTION_FACTOR(:)   = 0.0D0
    
     TIME_LOOP          = 0
     NUM_SUB_TIME_STEPS = 0
     INTERMED_RESULTS(:,:, :) = 0.0D0
     VOLUME_FRACTION(:)  = 0.0D0
    
     SOLUTE_FRACTIONS(:,:, :) = 0.0D0
     SOLID_CONCS(:,:)         = 0.0D0
    
     SINK_INTENSITY = 0.0D0
     DIFF_SINK(:)   = 0.0D0    !End of initializations
   
!  Switching derivatives: 
!     0 - switched off (is calculated but makes equal zero)
!     1 - is calculated 
    switch_kinetics   = 1
    switch_burial     = 1
    switch_partmixing = 0  !is made equal zero in the code, because of bug
    switch_diffusion  = 1
    switch_advection  = 0
    switch_settling   = 1
      
    !Checking state and output dimensions. Change them if they are changed
    if(NUM_SED_LAYERS .ne. 6) then
       print *, 'SEDIMENT_MODEL_1: Wrong number of layers'
       stop
    end if             
    
    if(NUM_SED_VARS .ne. 15) then
       print *, 'SEDIMENT_MODEL_1: Wrong number of state variables'
       stop
    end if  
    
    if(NUM_SED_OUTPUTS .ne. 17) then
       print *, 'SEDIMENT_MODEL_1: Wrong number of outputs'
       stop
    end if  
          
    WATER_DENSITY = 1.0D0
    
    do i = 1, NUM_SED_LAYERS
        SOLID_CONCS(:,i) = (SED_DENSITIES(:,i) - (WATER_DENSITY * SED_POROSITIES(:,i)))
    
        !Initialisation of solute fractions      
        do j = 1, NUM_SED_VARS
            !Fraction of variable in solute form. 
            !Not used for solids (Should be zero for solids)
            SOLUTE_FRACTIONS (:,i, j) = 1.0D0 
            PART_MIXING_RATES(:,i, j) = 0.0D0
        end do
    end do
    
    !IN_WHICH_PHASE
    !0 : SOLUTE
    !1 : SOLID
    !2 : ALL SEDIMENTS
        
    IN_WHICH_PHASE(1)  = 2 !2 
    IN_WHICH_PHASE(2)  = 0
    IN_WHICH_PHASE(3)  = 0
    IN_WHICH_PHASE(4)  = 1 !1
    IN_WHICH_PHASE(5)  = 2 !2
    IN_WHICH_PHASE(6)  = 0
    IN_WHICH_PHASE(7)  = 1 !1
    IN_WHICH_PHASE(8)  = 0
    IN_WHICH_PHASE(9)  = 0
    IN_WHICH_PHASE(10) = 1 !1
    IN_WHICH_PHASE(11) = 0
    IN_WHICH_PHASE(12) = 1 !1
    IN_WHICH_PHASE(13) = 0
    IN_WHICH_PHASE(14) = 0
    IN_WHICH_PHASE(15) = 0
    
    error = 0
     
    !call set_3d_d_array(NUM_SED_LAYERS,NUM_SED_VARS,NDIAGVAR_sed
    !		,PROCESSES_sed,0.0D+0)
    
    PROCESSES_sed(:,:,:,:) = 0.0D+0
    
   if(debug_stranger) then
    do i = 1, NUM_SED_LAYERS
        do j = 1, NUM_SED_VARS 
         do k = 1,nkn 
                  
            if(STRANGER(INIT_SED_STATE_VARS(k,I,J)) .eq. 1) then
                print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j,'Cell ',k 
                print *, 'Initial state is NaN'
                print *, 'INITIAL(i,j)=',INIT_SED_STATE_VARS(k,I,J)
                error =1
            end if
        end do
     end do
    end do

    if (error .eq. 1) stop        
   end if
    
    !INITIALIZE SEDIMENT STATE VARIABLES    
    SED_NH4N(:,:) = INIT_SED_STATE_VARS(:,:, 1)
    SED_NO3N(:,:) = INIT_SED_STATE_VARS(:,:, 2)
    SED_DON (:,:) = INIT_SED_STATE_VARS(:,:, 3)
    SED_PON (:,:) = INIT_SED_STATE_VARS(:,:, 4)
    SED_PO4P(:,:) = INIT_SED_STATE_VARS(:,:, 5)
    SED_DOP (:,:) = INIT_SED_STATE_VARS(:,:, 6)
    SED_POP (:,:) = INIT_SED_STATE_VARS(:,:, 7)
    SED_DOXY(:,:) = INIT_SED_STATE_VARS(:,:, 8)
    SED_DOC (:,:) = INIT_SED_STATE_VARS(:,:, 9)
    SED_POC (:,:) = INIT_SED_STATE_VARS(:,:, 10)
    SED_DSi (:,:) = INIT_SED_STATE_VARS(:,:, 11)
    SED_PSi (:,:) = INIT_SED_STATE_VARS(:,:, 12)
    INORG_C (:,:) = INIT_SED_STATE_VARS(:,:, 13)
    TOT_ALK (:,:) = INIT_SED_STATE_VARS(:,:, 14)
    SALT    (:,:) = INIT_SED_STATE_VARS(:,:, 15)
    
    !INITIALIZE SEDIMENT MODEL COEFFICIENTS
    if(NUM_SED_CONSTS .ne. NUM_SED_MOD_CONSTS) then
       print *, 'SEDIMENT_MODEL_1: Wrong number of constants'
    end if
    
!   Constants    
!     K_OXIC_DISS_POC      = SED_MODEL_CONSTANTS(1)      !1 ! PARTICULATE ORGANIC CARBON   Dissolution rate constant  at 20 C (aerobic) - 1/day    
!     K_ANOXIC_DISS_POC    = SED_MODEL_CONSTANTS(2)      !2 ! PARTICULATE ORGANIC CARBON   Dissolution rate constant n at 20 C (anoxic) - 1/day    
!     THETA_DISS_POC       = SED_MODEL_CONSTANTS(3)      !3 ! PARTICULATE ORGANIC CARBON   Temperature correction for dissolution                  
!     KHS_DISS_POC         = SED_MODEL_CONSTANTS(4)      !4 ! PARTICULATE ORGANIC CARBON   Half saturation concentration  for dissolution          
!     K_OXIC_DISS_PON      = SED_MODEL_CONSTANTS(5)      !5 ! PARTICULATE ORGANIC NITROGEN  Dissolution rate constant  at 20 C (aerobic) - 1/day   
!     K_ANOXIC_DISS_PON    = SED_MODEL_CONSTANTS(6)      !6 ! PARTICULATE ORGANIC NITROGEN  Dissolution rate constant  at 20 C (anoxic) - 1/day    
!     THETA_DISS_PON       = SED_MODEL_CONSTANTS(7)      !7 ! PARTICULATE ORGANIC NITROGEN  Temperature correction for dissolution                 
!     KHS_DISS_PON         = SED_MODEL_CONSTANTS(8)      !8 ! PARTICULATE ORGANIC NITROGEN  Half saturation concentration  for dissolution         
!     K_OXIC_DISS_POP      = SED_MODEL_CONSTANTS(9)      !9 ! PARTICULATE ORGANIC PHOSPHORUS  Dissolution rate constant  at 20 C (aerobic) - 1/day 
!     K_ANOXIC_DISS_POP    = SED_MODEL_CONSTANTS(10)     !10! PARTICULATE ORGANIC PHOSPHORUS  Dissolution rate constant  at 20 C (anoxic) - 1/day  
!     THETA_DISS_POP       = SED_MODEL_CONSTANTS(11)     !11! PARTICULATE ORGANIC PHOSPHORUS  Temperature correction for dissolution               
!     KHS_DISS_POP         = SED_MODEL_CONSTANTS(12)     !12! PARTICULATE ORGANIC PHOSPHORUS  Half saturation concentration   for dissolution      
!     K_OXIC_DISS_PSi      = SED_MODEL_CONSTANTS(13)     !13! PARTICULATE SILICON    Dissolution rate constant  at 20 C (aerobic) - 1/day          
!     K_ANOXIC_DISS_PSi    = SED_MODEL_CONSTANTS(14)     !14! PARTICULATE SILICON    Dissolution rate constant  at 20 C (anoxic) - 1/day           
!     THETA_DISS_PSi       = SED_MODEL_CONSTANTS(15)     !15! PARTICULATE SILICON    Temperature correction for dissolution                        
!     KHS_DISS_PSi         = SED_MODEL_CONSTANTS(16)     !16! PARTICULATE SILICON    Half saturation concentration  for dissolution                
!     K_OXIC_MINER_DOC     = SED_MODEL_CONSTANTS(17)     !17! DISSOLVED ORGANIC CARBON   Mineralization rate constant  at 20 C (aerobic) - 1/day   
!     K_ANOXIC_MINER_DOC   = SED_MODEL_CONSTANTS(18)     !18! DISSOLVED ORGANIC CARBON   Mineralization rate constant  at 20 C (anoxic) - 1/day    
!     THETA_MINER_DOC      = SED_MODEL_CONSTANTS(19)     !19! DISSOLVED ORGANIC CARBON   Temperature correction for dissolution                    
!     KHS_MINER_DOC        = SED_MODEL_CONSTANTS(20)     !20! DISSOLVED ORGANIC CARBON   Half saturation concentration  for mineralization         
!     K_OXIC_MINER_DON     = SED_MODEL_CONSTANTS(21)     !21! DISSOLVED ORGANIC NITROGEN  Mineralization rate constant  at 20 C (aerobic) - 1/day  
!     K_ANOXIC_MINER_DON   = SED_MODEL_CONSTANTS(22)     !22! DISSOLVED ORGANIC NITROGEN  Mineralization rate constant  at 20 C (anoxic) - 1/day   
!     THETA_MINER_DON      = SED_MODEL_CONSTANTS(23)     !23! DISSOLVED ORGANIC NITROGEN  Temperature correction                                   
!     KHS_MINER_DON        = SED_MODEL_CONSTANTS(24)     !24! DISSOLVED ORGANIC NITROGEN  Half saturation concentration  for mineralization        
!     K_OXIC_MINER_DOP     = SED_MODEL_CONSTANTS(25)     !25! DISSOLVED ORGANIC PHOSPHORUS Mineralization rate constant   at 20 C (aerobic) - 1/day
!     K_ANOXIC_MINER_DOP   = SED_MODEL_CONSTANTS(26)     !26! DISSOLVED ORGANIC PHOSPHORUS Mineralization rate constant   at 20 C (anoxic) - 1/day 
!     THETA_MINER_DOP      = SED_MODEL_CONSTANTS(27)     !27! DISSOLVED ORGANIC PHOSPHORUS Temperature correction for dissolution                  
!     KHS_MINER_DOP        = SED_MODEL_CONSTANTS(28)     !28! DISSOLVED ORGANIC PHOSPHORUS Half saturation concentration  for mineralization       
!     O_TO_C               = SED_MODEL_CONSTANTS(29)     !29! Oxygen to carbon ratio                                                               
!     K_NITR               = SED_MODEL_CONSTANTS(30)     !30! Nitrification rate constant at 20 C - 1/day                                          
!     THETA_NITR           = SED_MODEL_CONSTANTS(31)     !31! Temperature correction for nitrification                                             
!     KHS_NITR_NH4N        = SED_MODEL_CONSTANTS(32)     !32! Half saturation constant of nitrification for NH4N - mg/L N                          
!     KHS_NITR_DOXY        = SED_MODEL_CONSTANTS(33)     !33! Half saturation constant of nitrification for DOXY - mg/L O2                         
!     K_DENITR             = SED_MODEL_CONSTANTS(34)     !34! Denitrification rate constant at 20 C - 1/day                                        
!     THETA_DENITR         = SED_MODEL_CONSTANTS(35)     !35! Temperature correction for denitrification                                           
!     KHS_DENITR_NO3N      = SED_MODEL_CONSTANTS(36)     !36! Half saturation constant of denitrification for NO3N - mg/L N                        
!     KHS_DENITR_DOC       = SED_MODEL_CONSTANTS(37)     !37! Half saturation constant of denitrification for DOC - mg/L C                         
!     KHS_DENITR_DOXY      = SED_MODEL_CONSTANTS(38)     !38! Half saturation constant of denitrification for DOXY - mg/L O                        
!     DENITR_YIELD         = SED_MODEL_CONSTANTS(39)     !39! Denitrification yield                                                                
!     DOXY_AT_ANOXIA       = SED_MODEL_CONSTANTS(40)     !40! DOXY, under which anoxia begins - mg/L O2                                            
!     SOLID_PART_COEFF_NH4 = SED_MODEL_CONSTANTS(41)     !41! Solid part coeff for ammonium nitrogen (kg^-1)                                       
!     SOLID_PART_COEFF_PO4 = SED_MODEL_CONSTANTS(42)     !42! Solid part coeff for phosphate phosphorus (kg^-1)                                    

!   Constants    
      call para_get_value('K_OXIC_DISS_POC'       , K_OXIC_DISS_POC      ) !1 ! PARTICULATE ORGANIC CARBON   Dissolution rate constant  at 20 C (aerobic) - 1/day    
      call para_get_value('K_ANOXIC_DISS_POC'     , K_ANOXIC_DISS_POC    ) !2 ! PARTICULATE ORGANIC CARBON   Dissolution rate constant n at 20 C (anoxic) - 1/day    
      call para_get_value('THETA_DISS_POC'        , THETA_DISS_POC       ) !3 ! PARTICULATE ORGANIC CARBON   Temperature correction for dissolution                  
      call para_get_value('KHS_DISS_POC'          , KHS_DISS_POC         ) !4 ! PARTICULATE ORGANIC CARBON   Half saturation concentration  for dissolution          
      call para_get_value('K_OXIC_DISS_PON'       , K_OXIC_DISS_PON      ) !5 ! PARTICULATE ORGANIC NITROGEN  Dissolution rate constant  at 20 C (aerobic) - 1/day   
      call para_get_value('K_ANOXIC_DISS_PON'     , K_ANOXIC_DISS_PON    ) !6 ! PARTICULATE ORGANIC NITROGEN  Dissolution rate constant  at 20 C (anoxic) - 1/day    
      call para_get_value('THETA_DISS_PON'        , THETA_DISS_PON       ) !7 ! PARTICULATE ORGANIC NITROGEN  Temperature correction for dissolution                 
      call para_get_value('KHS_DISS_PON'          , KHS_DISS_PON         ) !8 ! PARTICULATE ORGANIC NITROGEN  Half saturation concentration  for dissolution         
      call para_get_value('K_OXIC_DISS_POP'       , K_OXIC_DISS_POP      ) !9 ! PARTICULATE ORGANIC PHOSPHORUS  Dissolution rate constant  at 20 C (aerobic) - 1/day 
      call para_get_value('K_ANOXIC_DISS_POP'     , K_ANOXIC_DISS_POP    ) !10! PARTICULATE ORGANIC PHOSPHORUS  Dissolution rate constant  at 20 C (anoxic) - 1/day  
      call para_get_value('THETA_DISS_POP'        , THETA_DISS_POP       ) !11! PARTICULATE ORGANIC PHOSPHORUS  Temperature correction for dissolution               
      call para_get_value('KHS_DISS_POP'          , KHS_DISS_POP         ) !12! PARTICULATE ORGANIC PHOSPHORUS  Half saturation concentration   for dissolution      
      call para_get_value('K_OXIC_DISS_PSi'       , K_OXIC_DISS_PSi      ) !13! PARTICULATE SILICON    Dissolution rate constant  at 20 C (aerobic) - 1/day          
      call para_get_value('K_ANOXIC_DISS_PSi'     , K_ANOXIC_DISS_PSi    ) !14! PARTICULATE SILICON    Dissolution rate constant  at 20 C (anoxic) - 1/day           
      call para_get_value('THETA_DISS_PSi'        , THETA_DISS_PSi       ) !15! PARTICULATE SILICON    Temperature correction for dissolution                        
      call para_get_value('KHS_DISS_PSi'          , KHS_DISS_PSi         ) !16! PARTICULATE SILICON    Half saturation concentration  for dissolution                
      call para_get_value('K_OXIC_MINER_DOC'      , K_OXIC_MINER_DOC     ) !17! DISSOLVED ORGANIC CARBON   Mineralization rate constant  at 20 C (aerobic) - 1/day   
      call para_get_value('K_ANOXIC_MINER_DOC'    , K_ANOXIC_MINER_DOC   ) !18! DISSOLVED ORGANIC CARBON   Mineralization rate constant  at 20 C (anoxic) - 1/day    
      call para_get_value('THETA_MINER_DOC'       , THETA_MINER_DOC      ) !19! DISSOLVED ORGANIC CARBON   Temperature correction for dissolution                    
      call para_get_value('KHS_MINER_DOC'         , KHS_MINER_DOC        ) !20! DISSOLVED ORGANIC CARBON   Half saturation concentration  for mineralization         
      call para_get_value('K_OXIC_MINER_DON'      , K_OXIC_MINER_DON     ) !21! DISSOLVED ORGANIC NITROGEN  Mineralization rate constant  at 20 C (aerobic) - 1/day  
      call para_get_value('K_ANOXIC_MINER_DON'    , K_ANOXIC_MINER_DON   ) !22! DISSOLVED ORGANIC NITROGEN  Mineralization rate constant  at 20 C (anoxic) - 1/day   
      call para_get_value('THETA_MINER_DON'       , THETA_MINER_DON      ) !23! DISSOLVED ORGANIC NITROGEN  Temperature correction                                   
      call para_get_value('KHS_MINER_DON'         , KHS_MINER_DON        ) !24! DISSOLVED ORGANIC NITROGEN  Half saturation concentration  for mineralization        
      call para_get_value('K_OXIC_MINER_DOP'      , K_OXIC_MINER_DOP     ) !25! DISSOLVED ORGANIC PHOSPHORUS Mineralization rate constant   at 20 C (aerobic) - 1/day
      call para_get_value('K_ANOXIC_MINER_DOP'    , K_ANOXIC_MINER_DOP   ) !26! DISSOLVED ORGANIC PHOSPHORUS Mineralization rate constant   at 20 C (anoxic) - 1/day 
      call para_get_value('THETA_MINER_DOP'       , THETA_MINER_DOP      ) !27! DISSOLVED ORGANIC PHOSPHORUS Temperature correction for dissolution                  
      call para_get_value('KHS_MINER_DOP'         , KHS_MINER_DOP        ) !28! DISSOLVED ORGANIC PHOSPHORUS Half saturation concentration  for mineralization       
      call para_get_value('O_TO_C'                , O_TO_C               ) !29! Oxygen to carbon ratio                                                               
      call para_get_value('K_NITR'                , K_NITR               ) !30! Nitrification rate constant at 20 C - 1/day                                          
      call para_get_value('THETA_NITR'            , THETA_NITR           ) !31! Temperature correction for nitrification                                             
      call para_get_value('KHS_NITR_NH4N'         , KHS_NITR_NH4N        ) !32! Half saturation constant of nitrification for NH4N - mg/L N                          
      call para_get_value('KHS_NITR_DOXY'         , KHS_NITR_DOXY        ) !33! Half saturation constant of nitrification for DOXY - mg/L O2                         
      call para_get_value('K_DENITR'              , K_DENITR             ) !34! Denitrification rate constant at 20 C - 1/day                                        
      call para_get_value('THETA_DENITR'          , THETA_DENITR         ) !35! Temperature correction for denitrification                                           
      call para_get_value('KHS_DENITR_NO3N'       , KHS_DENITR_NO3N      ) !36! Half saturation constant of denitrification for NO3N - mg/L N                        
      call para_get_value('KHS_DENITR_DOC'        , KHS_DENITR_DOC       ) !37! Half saturation constant of denitrification for DOC - mg/L C                         
      call para_get_value('KHS_DENITR_DOXY'       , KHS_DENITR_DOXY      ) !38! Half saturation constant of denitrification for DOXY - mg/L O                        
      call para_get_value('DENITR_YIELD'          , DENITR_YIELD         ) !39! Denitrification yield                                                                
      call para_get_value('DOXY_AT_ANOXIA'        , DOXY_AT_ANOXIA       ) !40! DOXY, under which anoxia begins - mg/L O2                                            
      call para_get_value('SOLID_PART_COEFF_NH4'  , SOLID_PART_COEFF_NH4 ) !41! Solid part coeff for ammonium nitrogen (kg^-1)                                       
      call para_get_value('SOLID_PART_COEFF_PO4'  , SOLID_PART_COEFF_PO4 ) !42! Solid part coeff for phosphate phosphorus (kg^-1)
    
    
    
!       print *, 'K_OXIC_DISS_POC     ',  K_OXIC_DISS_POC      
!       print *, 'K_ANOXIC_DISS_POC   ',  K_ANOXIC_DISS_POC    
!       print *, 'THETA_DISS_POC      ',  THETA_DISS_POC       
!       print *, 'KHS_DISS_POC        ',  KHS_DISS_POC         
!       print *, 'K_OXIC_DISS_PON     ',  K_OXIC_DISS_PON      
!       print *, 'K_ANOXIC_DISS_PON   ',  K_ANOXIC_DISS_PON    
!       print *, 'THETA_DISS_PON      ',  THETA_DISS_PON       
!       print *, 'KHS_DISS_PON        ',  KHS_DISS_PON         
!       print *, 'K_OXIC_DISS_POP     ',  K_OXIC_DISS_POP      
!       print *, 'K_ANOXIC_DISS_POP   ',  K_ANOXIC_DISS_POP    
!       print *, 'THETA_DISS_POP      ',  THETA_DISS_POP       
!       print *, 'KHS_DISS_POP        ',  KHS_DISS_POP         
!       print *, 'K_OXIC_DISS_PSi     ',  K_OXIC_DISS_PSi      
!       print *, 'K_ANOXIC_DISS_PSi   ',  K_ANOXIC_DISS_PSi    
!       print *, 'THETA_DISS_PSi      ',  THETA_DISS_PSi       
!       print *, 'KHS_DISS_PSi        ',  KHS_DISS_PSi         
!       print *, 'K_OXIC_MINER_DOC    ',  K_OXIC_MINER_DOC     
!       print *, 'K_ANOXIC_MINER_DOC  ',  K_ANOXIC_MINER_DOC   
!       print *, 'THETA_MINER_DOC     ',  THETA_MINER_DOC      
!       print *, 'KHS_MINER_DOC       ',  KHS_MINER_DOC        
!       print *, 'K_OXIC_MINER_DON    ',  K_OXIC_MINER_DON     
!       print *, 'K_ANOXIC_MINER_DON  ',  K_ANOXIC_MINER_DON   
!       print *, 'THETA_MINER_DON     ',  THETA_MINER_DON      
!       print *, 'KHS_MINER_DON       ',  KHS_MINER_DON        
!       print *, 'K_OXIC_MINER_DOP    ',  K_OXIC_MINER_DOP     
!       print *, 'K_ANOXIC_MINER_DOP  ',  K_ANOXIC_MINER_DOP   
!       print *, 'THETA_MINER_DOP     ',  THETA_MINER_DOP      
!       print *, 'KHS_MINER_DOP       ',  KHS_MINER_DOP        
!       print *, 'O_TO_C              ',  O_TO_C               
!       print *, 'K_NITR              ',  K_NITR               
!       print *, 'THETA_NITR          ',  THETA_NITR           
!       print *, 'KHS_NITR_NH4N       ',  KHS_NITR_NH4N        
!       print *, 'KHS_NITR_DOXY       ',  KHS_NITR_DOXY        
!       print *, 'K_DENITR            ',  K_DENITR             
!       print *, 'THETA_DENITR        ',  THETA_DENITR         
!       print *, 'KHS_DENITR_NO3N     ',  KHS_DENITR_NO3N      
!       print *, 'KHS_DENITR_DOC      ',  KHS_DENITR_DOC       
!       print *, 'KHS_DENITR_DOXY     ',  KHS_DENITR_DOXY      
!       print *, 'DENITR_YIELD        ',  DENITR_YIELD         
!       print *, 'DOXY_AT_ANOXIA      ',  DOXY_AT_ANOXIA       
!       print *, 'SOLID_PART_COEFF_NH4',  SOLID_PART_COEFF_NH4 
!       print *, 'SOLID_PART_COEFF_PO4',  SOLID_PART_COEFF_PO4       
!    stop
    !OLD(INITIAL) STATE VARIABLES ASIGNED TO ARRAY 'INTERMED_RESULTS'        
    INTERMED_RESULTS(:,:, :) = INIT_SED_STATE_VARS(:,:, :)
    
    NUM_SUB_TIME_STEPS = 1

    !LOOP ON SUBSTEPS (Substeps not used. For correct use averaging shoud be implemented) 
    do TIME_LOOP = 1, NUM_SUB_TIME_STEPS
    
        !UNIT AREA MASSES (all calculations inside on unit area masses, g/m2)
        do i = 1, NUM_SED_LAYERS
        
            do j = 1, NUM_SED_VARS
        
                if (IN_WHICH_PHASE(j).eq.0) then
                    VOLUME_FRACTION(:) = SED_POROSITIES(:,i)
                end if
                
                if (IN_WHICH_PHASE(j).eq.1) then
                    VOLUME_FRACTION(:) = 1.0D0 - SED_POROSITIES(:,i)
                end if
        
                if (IN_WHICH_PHASE(j).eq.2) then
                    VOLUME_FRACTION(:) = 1.0D0
                end if
        
                 !g/m3(~mg/l) * m =  g/m2
                UNIT_AREA_MASSES(:,i, j) = &
                     INTERMED_RESULTS(:,i, j) * SED_DEPTHS(:,i) * VOLUME_FRACTION(:)
            end do
           
        end do
        !UNIT AREA MASSES
        
        !PREPARATION FOR CALCULATION DERIVATIVES 
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                  
        do i = 1, NUM_SED_LAYERS
        
            !Calculate solute fractions of NH4 and PO4      
            SOLUTE_FRACTIONS(:,i, 1) = (1.0D0 / SED_POROSITIES(:,i)) * &
                (1.0D0 / (1.0D0 + (SOLID_CONCS(:,i) * SOLID_PART_COEFF_NH4)))
        
            SOLUTE_FRACTIONS(:,i, 5) = (1.0D0 / SED_POROSITIES(:,i)) * &
                (1.0D0 / (1.0D0 + (SOLID_CONCS(:,i) * SOLID_PART_COEFF_PO4)))
               
            do j = 1, NUM_SED_VARS
        
                if (IN_WHICH_PHASE(j).eq.0) then
                    VOLUME_FRACTION = SED_POROSITIES(:,i)
                end if
                
                if (IN_WHICH_PHASE(j).eq.1) then
                    VOLUME_FRACTION = 1.0D0 - SED_POROSITIES(:,i)
                end if
        
                if (IN_WHICH_PHASE(j).eq.2) then
                    VOLUME_FRACTION = 1.0D0
                end if
                
                ! PREPARATION FOR DifFERENT PHASES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
                !For phase 0 and 2          
                if ((IN_WHICH_PHASE(j).eq.0).or.(IN_WHICH_PHASE(j).eq.2)) then
            
                    if (i.eq.1) then  ! For the first layer
        
                        !ADVECTION
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_advection.ne.0) then
                        
                            if (ADVECTIVE_VELOCITY <= 0.0D0) then
                                ADV_ENTERING_CONC(:) = SURF_WATER_CONCS(:,j)
                            else
                                ADV_ENTERING_CONC(:) = &
                                    INTERMED_RESULTS(:,i + 1, j) * SOLUTE_FRACTIONS(:,i + 1, j)
                            end if
        
                        end if
                        
                        !DIFFUSION
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_diffusion.ne.0) then
                            NEIGHBOUR_CONC(:) = SURF_WATER_CONCS(:,j)
                            SED_MIXLEN     = SURF_MIXLEN
                        end if
                        
                    else !for  layer numbers  > 1
        
                       !ADVECTION ENTERING CONC.
                       !Entering concentration is from lower to upper layer
                       !in direction of the flow
                       
                       !24/06/2012: If clause Added by Ali for safer switch off    
                       if (switch_advection.ne.0) then
        
                           if (ADVECTIVE_VELOCITY <= 0.0D0) then
                               ADV_ENTERING_CONC(:) = INTERMED_RESULTS(:,i - 1, j) * &
                                                   SOLUTE_FRACTIONS(:,i - 1, j)
                           else
        
                               if (i.eq.NUM_SED_LAYERS) then !for the last layer
                            
                                   !Let nothing comes from below (corrected by Petras)
                                    ADV_ENTERING_CONC = 0.0D0 
                                    !      INTERMED_RESULTS(i, j) * 
                                    !      SOLUTE_FRACTIONS(i, j)
        
                                else !for the middle layers
                                    ADV_ENTERING_CONC(:) =  &
                                        INTERMED_RESULTS(:,i + 1, j) *SOLUTE_FRACTIONS(:,i + 1, j)
                                end if
        
                            end if
                        
                        end if
                        
                        !DIFFUSION
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_diffusion.ne.0) then
                            NEIGHBOUR_CONC(:) = &
                                INTERMED_RESULTS(:,i - 1, j) * SOLUTE_FRACTIONS(:,i - 1, j)
        
                            SED_MIXLEN(:)     = 0.5D0 * (SED_DEPTHS(:,i - 1) + SED_DEPTHS(:,i))
                        end if
                        
                    end if
        
                    !ADVECTION IN AND OUT RATES
                    !24/06/2012: If clause Added by Ali for safer switch off
                    if (switch_advection.ne.0) then
                        SED_IN_ADVEC_RATES(:,i, j) = &
                            ADV_ENTERING_CONC(:) * dabs(ADVECTIVE_VELOCITY)
        
                        SED_OUT_ADVEC_RATES(:,i, j) = INTERMED_RESULTS(:,i, j) * &
                            SOLUTE_FRACTIONS(:,i, j) * DABS(ADVECTIVE_VELOCITY)
                    end if
                        
                    !DIFFUSION
                    !24/06/2012: If clause Added by Ali for safer switch off
                    if (switch_diffusion.ne.0) then
                        UPPER_CONC_GRADIENT(:) = (INTERMED_RESULTS(:,i, j) * &
                               SOLUTE_FRACTIONS(:,i, j)) - NEIGHBOUR_CONC(:)
                               
                     if(debug_stranger) then         
                      do k= 1,nkn
                        if (STRANGER(UPPER_CONC_GRADIENT(k)).eq.1) then
                            print *, 'CELL NO                  : ', k
                            print *, 'LAYER NO                 : ', i
                            print *, 'STATE VARIABLE NO        : ', j
                            print *, 'UPPER_CONC_GRADIENT      : ', UPPER_CONC_GRADIENT
                            print *, '  INTERMED_RESULTS(i, j) : ',INTERMED_RESULTS(k,i, j)
                            print *, '  SOLUTE_FRACTIONS(i, j) : ',SOLUTE_FRACTIONS(k,i, j)
                            print *, '  NEIGHBOUR_CONC         : ', NEIGHBOUR_CONC(k)
                            error = 1
                        end if
                       end do                       
                       if(error .eq. 1) stop
                      end if
            
                        DIFF_CORRECTION_FACTOR(:) = 1.0D0 / &
                            (1.0D0 + (3.0D0 * (1.0D0 - SED_POROSITIES(:,i))))
            
                        SED_DIFFUSION_RATES(:,i, j) = DIFF_CORRECTION_FACTOR(:) * &
                            (UPPER_CONC_GRADIENT(:) * SED_DIFFUSIONS(:,i, j)) / SED_MIXLEN(:)
                            
                        ! Increasing diffusion rates from first layer
                        
                        if (i .eq. 1 ) then
                          where (SED_DIFFUSION_RATES(:,1, j) .ge. 0.D0) &
                             SED_DIFFUSION_RATES(:,i, j) = DIFF_DRAG * SED_DIFFUSION_RATES(:,i, j)
                          where (SED_DIFFUSION_RATES(:,1, j) .lt. 0.D0) &
                             SED_DIFFUSION_RATES(:,i, j) = (1.D0/DIFF_DRAG) * SED_DIFFUSION_RATES(:,i, j)
                        end if
                    end if
                    
                    !FLUX FROM SEDIMENTS TO WATER COLUMN
                    if (i.eq.1) then
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_diffusion.ne.0) then
                            FLUXES_FROM_SEDIMENTS(:,j) = SED_DIfFUSION_RATES(:,i, j)
                        end if
                        
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_advection.ne.0) then 
                        
                            if (ADVECTIVE_VELOCITY <= 0.0D0) then
                                FLUXES_FROM_SEDIMENTS(:,j) = &
                                    FLUXES_FROM_SEDIMENTS(:,j) - SED_IN_ADVEC_RATES (:,i, j)
                            else
                                FLUXES_FROM_SEDIMENTS(:,j) = &
                                    FLUXES_FROM_SEDIMENTS(:,j) + SED_OUT_ADVEC_RATES(:,i, j)
                            end if
                        
                        end if
                        
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if ((switch_diffusion.eq.0).and.(switch_advection.eq.0)) then
                            FLUXES_FROM_SEDIMENTS(:,j) = 0.0D0
                        end if
                    end if
        
                    !FOR PARTICLE MIXING
                    if (IN_WHICH_PHASE(j).eq.0) then
                    
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_partmixing.ne.0) then
                            PART_MIXING_RATES(:,i, j) = 0.0D0
                        end if
                    
                    end if
        
                    if (IN_WHICH_PHASE(j).eq.2) then
        
                        if (I.eq.1) then
                            
                            !24/06/2012: If clause Added by Ali for safer switch off
                            if (switch_partmixing.ne.0) then
                                PART_MIXING_RATES(:,i, j) = 0.0D0
                            end if
                            
                        else
                        
                            !24/06/2012: If clause Added by Ali for safer switch off
                            if (switch_partmixing.ne.0) then
                                NEIGHBOUR_CONC(:) = &
                                    INTERMED_RESULTS(:,I - 1, J) * (1.0D0 - SOLUTE_FRACTIONS(:,i, j))
                                
                                SED_MIXLEN(:) = 0.5D0 * (SED_DEPTHS(:,I - 1) + SED_DEPTHS(:,i))
                                
                                UPPER_CONC_GRADIENT(:) = & 
                                    (INTERMED_RESULTS(:,i, j) * (1.0D0 - SOLUTE_FRACTIONS(:,i, j))) - &
                                    NEIGHBOUR_CONC
                                
                               if(debug_stranger) then    
                                do k=1,nkn
                                    if (STRANGER(UPPER_CONC_GRADIENT(k)).eq.1) then
                                        print *, 'CELL NO                  : ', CELLNO
                                        print *, 'LAYER NO                 : ', i
                                        print *, 'STATE VARIABLE NO        : ', j
                                        print *, 'UPPER_CONC_GRADIENT      : ', UPPER_CONC_GRADIENT
                                        print *, '  INTERMED_RESULTS(i, j) : ', INTERMED_RESULTS(k,i, j)
                                        print *, '  SOLUTE_FRACTIONS(i, j) : ',SOLUTE_FRACTIONS(k,i, j)
                                        print *, '  NEIGHBOUR_CONC         : ', NEIGHBOUR_CONC(k)
                                        error=1
                                    end if
                                end do
                                if(error .eq. 1) stop
                               end if
                                
                                PART_MIXING_RATES(:,i, j) = &
                                    (UPPER_CONC_GRADIENT * PART_MIXING_COEFFS(:,i, j)) / SED_MIXLEN
                            end if
                            
                        end if
        
                    end if
        
                else !for other phases(1)
                    FLUXES_FROM_SEDIMENTS (:,j) = 0.0D0
                    SED_DIFFUSION_RATES(:,i, j) = 0.0D0
                   
                    if (I.eq.1) then
                        
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_partmixing.ne.0) then
                            PART_MIXING_RATES(:,i, j) = 0.0D0
                        end if
                    
                    else
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_partmixing.ne.0) then
                            NEIGHBOUR_CONC(:) = INTERMED_RESULTS(:,i - 1, j)
                            SED_MIXLEN(:)     = 0.5D0 * (SED_DEPTHS(:,I - 1) + SED_DEPTHS(:,i))
                            
                            UPPER_CONC_GRADIENT(:) = &
                                  (INTERMED_RESULTS(:,i, j) * SOLUTE_FRACTIONS(:,i, j)) - NEIGHBOUR_CONC(:)
                            
                           if(debug_stranger) then      
                            do k=1,nkn                            
                             if (STRANGER(UPPER_CONC_GRADIENT(k)) .eq. 1) then
                                 print *, 'CELL NO                  : ', k
                                 print *, 'LAYER NO                 : ', i
                                 print *, 'STATE VARIABLE NO        : ', j
                                 print *, 'UPPER_CONC_GRADIENT      : ', UPPER_CONC_GRADIENT(k)
                                 print *, '  INTERMED_RESULTS(i, j) : ', INTERMED_RESULTS(k,i, j)
                                 print *, '  SOLUTE_FRACTIONS(i, j) : ', SOLUTE_FRACTIONS(k,i, j)
                                 print *, '  NEIGHBOUR_CONC         : ', NEIGHBOUR_CONC(k)
                                 error = 1
                             end if                            
                            end do
                            if(error .eq. 1) stop
                           end if
                            
                            PART_MIXING_RATES(:,i, j) = &
                                  (UPPER_CONC_GRADIENT(:) * PART_MIXING_COEFFS(:,i, j)) / SED_MIXLEN(:)
                        end if !switch_partmixing
                        
                    end if
        
                end if !end PREPARATION FOR DifFERENT PHASES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               
        
                !24/06/2012: If clause Added by Ali for safer switch off
                if (switch_burial.ne.0) then
                    SED_BURRIAL_RATES(:,i, j) = INTERMED_RESULTS(:,i, j) * SED_BURRIALS(:,i) 
                end if
                
                DERIVS           (:,i, j) = 0.0D0
                TRANSPORT_DERIVS (:,i, j) = 0.0D0
                KINETIC_DERIVS   (:,i, j) = 0.0D0
                SETTLING_DERIVS  (:,i, j) = 0.0D0
            end do
        
        end do  !end PREPARATION FOR CALCULATION DERIVATIVES 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
 

!---------------------------------------------------------------------------------------------------------
       !TRANSPORT DERIVATIVES
        do i = 1, (NUM_SED_LAYERS - 1)
        
            do j = 1, NUM_SED_VARS
        
                if (I.GT.1) then 
        
                    !TRANSPORT_DERIVS(i, j) = &
                    !    SED_DIFFUSION_RATES(i + 1, j) - SED_DIFFUSION_RATES(i, j) + & 
                    !    SED_BURRIAL_RATES  (i - 1, j) - SED_BURRIAL_RATES  (i, j) + &
                    !    PART_MIXING_RATES  (i + 1, j) - PART_MIXING_RATES  (i, j)
        
                    !Introduced by Petras for diagnostic output     
                    DIFFUSION_DERIVS  (:,i, j) = &
                        SED_DIFFUSION_RATES(:,i + 1, j) - SED_DIFFUSION_RATES(:,i, j)
        
                    BURIAL_DERIVS     (:,i, j) = &
                        SED_BURRIAL_RATES  (:,i - 1, j) - SED_BURRIAL_RATES  (:,i, j)
        
                    PART_MIXING_DERIVS(:,i, j) = &
                        PART_MIXING_RATES  (:,i + 1, j) - PART_MIXING_RATES  (:,i, j)
                        
                    ! Introduced by Petras to avoid double calculation                        
                    TRANSPORT_DERIVS(:,i, j) = &
                            DIFFUSION_DERIVS  (:,i, j) + &
                            BURIAL_DERIVS     (:,i, j) + &
                            PART_MIXING_DERIVS(:,i, j)
                        
                                              
                else ! i=1
        
!                    TRANSPORT_DERIVS(i, j) = &
!                        SED_DIFFUSION_RATES(i + 1, j) - SED_DIFFUSION_RATES(i, j)     - &
!                        SED_BURRIAL_RATES  (i, j)     + PART_MIXING_RATES  (i + 1, j) - &
!                        PART_MIXING_RATES  (i, j)
        
                    !Introduced by Petras for diagnostic output     
                    DIFFUSION_DERIVS  (:,i, j) =  &
                        SED_DIFFUSION_RATES(:,i + 1, j) - SED_DIFFUSION_RATES(:,i, j)
        
                    BURIAL_DERIVS     (:,i, j) = (-1.0D0)*SED_BURRIAL_RATES  (:,i, j)
        
                    PART_MIXING_DERIVS(:,i, j) = &                 
                        PART_MIXING_RATES  (:,i + 1, j) - PART_MIXING_RATES  (:,i, j)
                        
                    ! Introduced by Petras to avoid double calculation. It should be done everywhere fixme
                    TRANSPORT_DERIVS(:,i, j) = &
                            DIFFUSION_DERIVS  (:,i, j) + &
                            BURIAL_DERIVS     (:,i, j) + &
                            PART_MIXING_DERIVS(:,i, j)
                            
                end if
        
                !CONSIDER ADVECTION
        !        TRANSPORT_DERIVS(i, j) = TRANSPORT_DERIVS(i, j) + &
        !            SED_IN_ADVEC_RATES(i, j) - SED_OUT_ADVEC_RATES(i, j)
        
                   
                ADVECTION_DERIVS(:,i, j) = &
                    SED_IN_ADVEC_RATES(:,i, j) - SED_OUT_ADVEC_RATES(:,i, j)
                    
                TRANSPORT_DERIVS(:,i, j) = TRANSPORT_DERIVS(:,i, j) + ADVECTION_DERIVS(:,i, j)
            end do
        
        end do
        
        !Last layers processing (why mixing is not processed here?)     
        do j = 1, NUM_SED_VARS
        
            if (NUM_SED_LAYERS.GT.1) then
             
                 SINK_INTENSITY = 1.D0 !fixme should not be hardcoded
                 
                 DIFF_SINK(:) = (-1)*SINK_INTENSITY*  DABS(SED_DIFFUSION_RATES(:,NUM_SED_LAYERS, j)) !Lower boundary gradient
                 DIFFUSION_DERIVS(:,NUM_SED_LAYERS, j) = &
                              DIFF_SINK(:) - SED_DIFFUSION_RATES(:,NUM_SED_LAYERS, j)
        
                TRANSPORT_DERIVS(:,NUM_SED_LAYERS, j) = &
                    DIFFUSION_DERIVS  (:,NUM_SED_LAYERS, j) + &
                    SED_BURRIAL_RATES  (:,NUM_SED_LAYERS - 1, j) - &
                    SED_BURRIAL_RATES  (:,NUM_SED_LAYERS, j) + &
                    SED_IN_ADVEC_RATES (:,NUM_SED_LAYERS, j) - &!Added by Petras(last layer was not processed for advection)
                    SED_OUT_ADVEC_RATES(:,NUM_SED_LAYERS, j) 
                    !PART_MIXING_RATES  (i + 1, j)          !Added by Petras(last layer was not processed. Rates for the last layer not calculated! fixme
                    !PART_MIXING_RATES  (i, j)     *        
           
        
                    
                 BURIAL_DERIVS(:,NUM_SED_LAYERS, j) = &              
                    SED_BURRIAL_RATES (:,NUM_SED_LAYERS - 1, j) - &
                    SED_BURRIAL_RATES (:,NUM_SED_LAYERS, j)
        
                 PART_MIXING_DERIVS(:,NUM_SED_LAYERS, j) = 0. !temporary fix 
          
        
                 ADVECTION_DERIVS(:,NUM_SED_LAYERS, j) = &
                    SED_IN_ADVEC_RATES (:,NUM_SED_LAYERS, j) - &
                    SED_OUT_ADVEC_RATES(:,NUM_SED_LAYERS, j)
            else
                TRANSPORT_DERIVS(:,NUM_SED_LAYERS, j) = &
                    ((-1.0D0) * SED_DIFFUSION_RATES(:,NUM_SED_LAYERS, j)) + &
                    SED_BURRIAL_RATES  (:,NUM_SED_LAYERS, j) + &
                    SED_IN_ADVEC_RATES (:,NUM_SED_LAYERS, j) - & !Added by Petras(last layer was not processed for advection)
                    SED_OUT_ADVEC_RATES(:,NUM_SED_LAYERS, j) 
            
                !Introduced by Petras for testing:
                DIFFUSION_DERIVS(:,NUM_SED_LAYERS, j) = &
                    ((-1.0D0) * SED_DIFFUSION_RATES(:,NUM_SED_LAYERS, J))
        
                BURIAL_DERIVS   (:,NUM_SED_LAYERS, j) = &
                    SED_BURRIAL_RATES (:,NUM_SED_LAYERS, j)
        
                ADVECTION_DERIVS(:,NUM_SED_LAYERS, j) = &
                    SED_IN_ADVEC_RATES (:,NUM_SED_LAYERS, j) - &
                    SED_OUT_ADVEC_RATES(:,NUM_SED_LAYERS, j)
        
                PART_MIXING_DERIVS(:,NUM_SED_LAYERS, J) = 0. !temporary fix                
            end if         
        
        end do
        

   
       !end TRANSPORT DERIVATIVES
!---------------------------------------------------------------------------------------------------




        
        !Particle mixing derivatives do not work correctly because of a bug and
        !are therefore removed for now.
        
        !THE FOLLOWING COMMENTED LINES ARE TO BE DELETED WHEN SAFE SWITCH
        !OFF ALGORITHMS INTRODUECED BY ALI AT 24/06/2012 IN PETRAS HOUSE ARE
        !COMPLETED AND TESTED    
        ! Switch transport derivatives to zero if necessary
        !if (switch_burial     .eq. 0) BURIAL_DERIVS     (:,:) = 0.
        !if (switch_partmixing .eq. 0) PART_MIXING_DERIVS(:,:) = 0.
        !if (switch_diffusion  .eq. 0) DIFFUSION_DERIVS  (:,:) = 0.
        !if (switch_advection  .eq. 0) ADVECTION_DERIVS  (:,:) = 0.
        !
        !TRANSPORT_DERIVS(:,:) = DIFFUSION_DERIVS(:,:) + &
        !                        BURIAL_DERIVS(:,:)  + &
        !                        ADVECTION_DERIVS(:,:) + &
        !                        BURIAL_DERIVS(:,:) + &
        !                        PART_MIXING_DERIVS(:,:)
        
        
        
        
!kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk        
                              !KINETIC DERIVATIVES
        
        !24/06/2012: If clause Added by Ali for safer switch off
        if (switch_kinetics .ne. 0) then                      
            if (RUN_CO2SYS .eq. 1) then
                        
!$OMP PARALLEL PRIVATE       &
!$OMP (CO2SYS_PAR1          ,&
!$OMP  CO2SYS_PAR2          ,&
!$OMP  CO2SYS_PAR1TYPE      ,&
!$OMP  CO2SYS_PAR2TYPE      ,&
!$OMP  CO2SYS_SALT          ,&
!$OMP  CO2SYS_TEMPIN        ,&
!$OMP  CO2SYS_TEMPOUT       ,&
!$OMP  CO2SYS_PRESIN        ,&
!$OMP  CO2SYS_PRESOUT       ,&
!$OMP  CO2SYS_SI            ,&
!$OMP  CO2SYS_PO4           ,&
!$OMP  CO2SYS_pHSCALEIN     ,&
!$OMP  CO2SYS_K1K2CONSTANTS ,&
!$OMP  CO2SYS_KSO4CONSTANTS, CO2SYS_OUT_DATA,    &
!$OMP  CO2SYS_NICEHEADERS, CO2SYS_ntps)
                                           
 
!$OMP DO SCHEDULE(STATIC)

         
                do i = 1,NUM_SED_LAYERS
                    allocate(CO2SYS_PAR1         (nkn), &
                             CO2SYS_PAR2         (nkn), &
                             CO2SYS_PAR1TYPE     (nkn), &
                             CO2SYS_PAR2TYPE     (nkn), &
                             CO2SYS_SALT         (nkn), &
                             CO2SYS_TEMPIN       (nkn), &
                             CO2SYS_TEMPOUT      (nkn), &
                             CO2SYS_PRESIN       (nkn), &
                             CO2SYS_PRESOUT      (nkn), &
                             CO2SYS_SI           (nkn), &
                             CO2SYS_PO4          (nkn), &
                             CO2SYS_pHSCALEIN    (nkn), &
                             CO2SYS_K1K2CONSTANTS(nkn), &
                             CO2SYS_KSO4CONSTANTS(nkn))        
            
                    call ASSIGN_DBL_VECTOR_CONTENT(CO2SYS_PAR1, TOT_ALK(:,i) * 1.0D6)

                    CO2SYS_PAR2         (:) = INORG_C(:,i) * 1.0D6
                    CO2SYS_PAR1TYPE     (:) = 1
                    CO2SYS_PAR2TYPE     (:) = 2
                    CO2SYS_SALT         (:) = SALT(:,i)
                    CO2SYS_TEMPIN       (:) = SED_TEMPS(:,i)
                    CO2SYS_TEMPOUT      (:) = 0.0D0  !Does not matter for this example
                    CO2SYS_PRESIN       (:) = 0.0D0  !Does not matter for this example
                    CO2SYS_PRESOUT      (:) = 0.0D0  !Does not matter for this example
                    CO2SYS_SI           (:) = (SED_DSI (:,i) / 28.0855D0) * 1.0D3
                    CO2SYS_PO4          (:) = (SED_PO4P(:,i) / 30.9737D0) * 1.0D3
                    CO2SYS_pHSCALEIN    (:) = 1
                    CO2SYS_K1K2CONSTANTS(:) = 4
                    CO2SYS_KSO4CONSTANTS(:) = 1

                    CO2SYS_ntps = nkn
            
                    call CO2SYS(CO2SYS_PAR1         , CO2SYS_PAR2  , CO2SYS_PAR1TYPE , &
                                CO2SYS_PAR2TYPE     , CO2SYS_SALT  , CO2SYS_TEMPIN   , &
                                CO2SYS_TEMPOUT      , CO2SYS_PRESIN, CO2SYS_PRESOUT  , &
                                CO2SYS_SI           , CO2SYS_PO4   , CO2SYS_pHSCALEIN, &
                                CO2SYS_K1K2CONSTANTS, CO2SYS_KSO4CONSTANTS, CO2SYS_OUT_DATA , &
                                CO2SYS_NICEHEADERS  , &
                                CO2SYS_ntps)

                    pH         (:,i) = CO2SYS_OUT_DATA(:, 18)
                    K_ONE_TIP  (:,i) = CO2SYS_OUT_DATA(:, 75)
                    K_TWO_TIP  (:,i) = CO2SYS_OUT_DATA(:, 76)
                    K_THREE_TIP(:,i) = CO2SYS_OUT_DATA(:, 77)    
                        
                    if(allocated(CO2SYS_PAR1            )) deallocate(CO2SYS_PAR1            )
                    if(allocated(CO2SYS_PAR2            )) deallocate(CO2SYS_PAR2            )
                    if(allocated(CO2SYS_PAR1TYPE        )) deallocate(CO2SYS_PAR1TYPE        )
                    if(allocated(CO2SYS_PAR2TYPE        )) deallocate(CO2SYS_PAR2TYPE        )
                    if(allocated(CO2SYS_SALT            )) deallocate(CO2SYS_SALT            )
                    if(allocated(CO2SYS_TEMPIN          )) deallocate(CO2SYS_TEMPIN          )
                    if(allocated(CO2SYS_TEMPOUT         )) deallocate(CO2SYS_TEMPOUT         )
                    if(allocated(CO2SYS_PRESIN          )) deallocate(CO2SYS_PRESIN          )
                    if(allocated(CO2SYS_PRESOUT         )) deallocate(CO2SYS_PRESOUT         )
                    if(allocated(CO2SYS_SI              )) deallocate(CO2SYS_SI              )
                    if(allocated(CO2SYS_PO4             )) deallocate(CO2SYS_PO4             )
                    if(allocated(CO2SYS_pHSCALEIN       )) deallocate(CO2SYS_pHSCALEIN       )
                    if(allocated(CO2SYS_K1K2CONSTANTS   )) deallocate(CO2SYS_K1K2CONSTANTS   )
                    if(allocated(CO2SYS_KSO4CONSTANTS   )) deallocate(CO2SYS_KSO4CONSTANTS   )
                    if(allocated(CO2SYS_OUT_DATA        )) deallocate(CO2SYS_OUT_DATA        )
                    if(allocated(CO2SYS_NICEHEADERS     )) deallocate(CO2SYS_NICEHEADERS     )
                end do !layers 
!$OMP END DO NOWAIT	
!$OMP END PARALLEL
            end if !RUN_CO2SYS

            call SED_POC_DISSOLUTION &
                 (K_OXIC_DISS_POC  , &
                  K_ANOXIC_DISS_POC, &
                  THETA_DISS_POC   , &
                  KHS_DISS_POC     , &
                  DOXY_AT_ANOXIA   , &
                  SED_POC          , &
                  SED_DOXY         , &
                  SED_TEMPS        , &
                  nkn              , &
                  NUM_SED_LAYERS   , &
                  R_DISS_POC)

           call SED_DOC_MINERALIZATION &
                (K_OXIC_MINER_DOC  , &
                 K_ANOXIC_MINER_DOC, &
                 THETA_MINER_DOC   , &
                 KHS_MINER_DOC     , &
                 DOXY_AT_ANOXIA    , &
                 SED_DOC           , &
                 SED_DOXY          , &
                 SED_TEMPS         , &
                 nkn               , &
                 NUM_SED_LAYERS    , &
                 R_MINER_DOC)

            where (SED_DOXY.GE.DOXY_AT_ANOXIA)

                R_DISS_PON  = K_OXIC_DISS_PON    * (THETA_DISS_PON  ** (SED_TEMPS - 2.0D1)) * &
                    (SED_PON / (SED_PON + KHS_DISS_PON))  * SED_PON
            
                R_DISS_POP  = K_OXIC_DISS_POP    * (THETA_DISS_POP  ** (SED_TEMPS - 2.0D1)) * &
                    (SED_POP / (SED_POP + KHS_DISS_POP))  * SED_POP
            
                R_DISS_PSi  = K_OXIC_DISS_PSi    * (THETA_DISS_PSi  ** (SED_TEMPS - 2.0D1)) * &
                    (SED_PSi / (SED_PSi + KHS_DISS_PSi))  * SED_PSi
            
                DEOXYGENATION = O_TO_C * R_MINER_DOC
            
                R_MINER_DON = K_OXIC_MINER_DON   * (THETA_MINER_DON ** (SED_TEMPS - 2.0D1)) * &
                    (SED_DON / (SED_DON + KHS_MINER_DON)) * SED_DON
            
                R_MINER_DOP = K_OXIC_MINER_DOP   * (THETA_MINER_DOP ** (SED_TEMPS - 2.0D1)) * &
                    (SED_DOP / (SED_DOP + KHS_MINER_DOP)) * SED_DOP
            end where
            
            where (SED_DOXY.LT.DOXY_AT_ANOXIA)
            
                R_DISS_PON  = K_ANOXIC_DISS_PON  * (THETA_DISS_PON  ** (SED_TEMPS - 2.0D1)) * &
                    (SED_PON / (SED_PON + KHS_DISS_PON))  * SED_PON

                R_DISS_POP  = K_ANOXIC_DISS_POP  * (THETA_DISS_POP  ** (SED_TEMPS - 2.0D1)) * &
                    (SED_POP / (SED_POP + KHS_DISS_POP))  * SED_POP

                R_DISS_PSi  = K_ANOXIC_DISS_PSi  * (THETA_DISS_PSi  ** (SED_TEMPS - 2.0D1)) * &
                    (SED_PSi / (SED_PSi + KHS_DISS_PSi))  * SED_PSi
            
                DEOXYGENATION = (O_TO_C * R_MINER_DOC) * (SED_DOXY / &
                    (SED_DOXY + (DOXY_AT_ANOXIA / 2.0D0)))
            
                R_MINER_DON = K_ANOXIC_MINER_DON * (THETA_MINER_DON ** (SED_TEMPS - 2.0D1)) * &
                    (SED_DON / (SED_DON + KHS_MINER_DON)) * SED_DON
            
                R_MINER_DOP = K_ANOXIC_MINER_DOP * (THETA_MINER_DOP ** (SED_TEMPS - 2.0D1)) * &
                    (SED_DOP / (SED_DOP + KHS_MINER_DOP)) * SED_DOP
            end where
                
            !Nitrification
            R_NITR   = K_NITR * (THETA_NITR ** (SED_TEMPS - 2.0D1)) * &
                       (SED_DOXY / (SED_DOXY + KHS_NITR_DOXY)) * &
                       (SED_NH4N / (SED_NH4N + KHS_NITR_NH4N)) * SED_NH4N
            
            R_DENITR = K_DENITR * (THETA_DENITR ** (SED_TEMPS - 2.0D1)) * &
                       (KHS_DENITR_DOXY / (SED_DOXY + KHS_DENITR_DOXY))* &
                       (SED_NO3N / (SED_NO3N + KHS_DENITR_NO3N)) * &
                       (SED_DOC  / (SED_DOC  + KHS_DENITR_DOC))  * SED_NO3N
            
            R_DISS_POC    = R_DISS_POC    * SED_DEPTHS * (1.0D0 - SED_POROSITIES)
            R_DISS_PON    = R_DISS_PON    * SED_DEPTHS * (1.0D0 - SED_POROSITIES)
            R_DISS_POP    = R_DISS_POP    * SED_DEPTHS * (1.0D0 - SED_POROSITIES)
            R_DISS_PSi    = R_DISS_PSi    * SED_DEPTHS * (1.0D0 - SED_POROSITIES)
            R_MINER_DOC   = R_MINER_DOC   * SED_DEPTHS * SED_POROSITIES
            R_MINER_DON   = R_MINER_DON   * SED_DEPTHS * SED_POROSITIES
            R_MINER_DOP   = R_MINER_DOP   * SED_DEPTHS * SED_POROSITIES
            DEOXYGENATION = DEOXYGENATION * SED_DEPTHS * SED_POROSITIES    
            R_NITR        = R_NITR        * SED_DEPTHS * SED_POROSITIES
            R_DENITR      = R_DENITR      * SED_DEPTHS * SED_POROSITIES
            
            KINETIC_DERIVS(:,:, 1)  = R_MINER_DON(:,:) - R_NITR(:,:)
            KINETIC_DERIVS(:,:, 2)  = R_NITR(:,:) - R_DENITR(:,:)
            KINETIC_DERIVS(:,:, 3)  = R_DISS_PON(:,:) - R_MINER_DON(:,:)
            KINETIC_DERIVS(:,:, 4)  = (-1.0D0) * R_DISS_PON(:,:)
            KINETIC_DERIVS(:,:, 5)  = R_MINER_DOP(:,:)
            KINETIC_DERIVS(:,:, 6)  = R_DISS_POP(:,:) - R_MINER_DOP(:,:)
            KINETIC_DERIVS(:,:, 7)  = (-1.0D0) * R_DISS_POP(:,:)
            KINETIC_DERIVS(:,:, 8)  = ((-4.57D0) * R_NITR(:,:)) - DEOXYGENATION(:,:)
            KINETIC_DERIVS(:,:, 9)  = R_DISS_POC(:,:) - R_MINER_DOC(:,:) - (R_DENITR(:,:) / DENITR_YIELD)
            KINETIC_DERIVS(:,:, 10) = (-1.0D0) * R_DISS_POC(:,:)
            KINETIC_DERIVS(:,:, 11) = R_DISS_PSi(:,:)
            KINETIC_DERIVS(:,:, 12) = (-1.0D0) * R_DISS_PSi(:,:)
            
            ! Kinetic sub model for dissolved inorganic carbon
            TOTAL_DIC_KINETIC_SOURCES(:,:) = R_MINER_DOC(:,:)
            TOTAL_DIC_KINETIC_SINKS  (:,:) = 0.0D0
            
            ! Calculate the total derivative and convert it to moles
            if (CONSIDER_INORG_C_DERIVATIVE > 0) then
                PROCESSES_sed(:,:, 13, 1) = TOTAL_DIC_KINETIC_SOURCES /12000.0D0
                PROCESSES_sed(:,:, 13, 2) = TOTAL_DIC_KINETIC_SINKS   /12000.0D0
            
                DIC_KINETIC_DERIVATIVE = &
                    ((TOTAL_DIC_KINETIC_SOURCES - TOTAL_DIC_KINETIC_SINKS) / 12000.0D0)
                
                KINETIC_DERIVS(:,:, 13) = DIC_KINETIC_DERIVATIVE
            else
                PROCESSES_sed (:,:, 13, 1) = 0.0D0
                PROCESSES_sed (:,:, 13, 2) = 0.0D0
                KINETIC_DERIVS(:,:, 13)    = 0
            end if
            
            
            ! -------------------------------------------------------------------------
            ! KINETIC SUBMODEL FOR ALKALINITY
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! PREPARE FOR ALKALINITY-NITROGEN INTERACTIONS
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! Calculate the NH4 and NH3 fractions in ammonia
            ! -------------------------------------------------------------------------
            T_A      = SED_TEMPS + 2.7316D2
            pKH      = 9.018D-2 + (2.72992D3 / T_A)
            FRAC_NH3 = 1.0D0 / (1.0D0 + (10.0D0 ** (pKH - pH)))
            FRAC_NH4 = 1.0D0 - FRAC_NH3
            ! -------------------------------------------------------------------------
            ! End of calculate NH4 and NH3 fractions in ammonia
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! END OF PREPARE FOR ALKALINITY-NITROGEN INTERACTIONS
            ! -------------------------------------------------------------------------
            
            
            ! -------------------------------------------------------------------------
            ! NITROGEN BASED SOURCES
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! Calculate the alkalinity gain by ammonium generation
            ! (1 eq alk for each ammonium generated since one positive ion is gained)
            ! -------------------------------------------------------------------------
            N_CHEM_AUT_BAC_TOT_RESP   = 0.0D0
            N_AER_HET_BAC_INT_RESP    = 0.0D0
            N_FAC_AN_HET_BAC_TOT_RESP = 0.0D0
            N_DIA_TOT_RESP            = 0.0D0
            N_CYN_TOT_RESP            = 0.0D0
            N_OPA_TOT_RESP            = 0.0D0
            N_FIX_CYN_TOT_RESP        = 0.0D0
            N_AER_HET_BAC_N_OX        = 0.0D0
            N_FAC_AN_HET_BAC_N_OX     = 0.0D0
            N_ZOO_TOT_RESP            = 0.0D0
            N_ABIOTIC_DON_MIN         = R_MINER_DON * FRAC_NH4
            
            ALK_GAINED_BY_AMMONIUM_GEN     = &
                (N_CHEM_AUT_BAC_TOT_RESP   + N_AER_HET_BAC_INT_RESP + &
                 N_FAC_AN_HET_BAC_TOT_RESP + N_DIA_TOT_RESP         + &
                 N_CYN_TOT_RESP            + N_OPA_TOT_RESP         + &
                 N_FIX_CYN_TOT_RESP        + N_AER_HET_BAC_N_OX     + &
                 N_FAC_AN_HET_BAC_N_OX     + N_ZOO_TOT_RESP         + &
                 N_ABIOTIC_DON_MIN) / 14007.0D0
            ! -------------------------------------------------------------------------
            ! End of calculate the alkalinity gain by ammonium generation
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! Calculate the alkality gain by nitrate consumption
            ! 1 eq alk for each nitrate consumed since one negative ion is lost
            ! -------------------------------------------------------------------------
            N_DENITRIFICATION     = R_DENITR
            N_DIA_GROWTH          = 0.0D0
            N_CYN_GROWTH          = 0.0D0
            N_OPA_GROWTH          = 0.0D0
            N_NON_FIX_CYN_GROWTH  = 0.0D0
            N_AER_HET_BAC_GROWTH  = 0.0D0
            
           ALK_GAINED_BY_NITRATE_CONS = &
                (N_DENITRIFICATION     + N_DIA_GROWTH          + &
                 N_CYN_GROWTH          + N_OPA_GROWTH          + &
                 N_NON_FIX_CYN_GROWTH  + N_AER_HET_BAC_GROWTH) / 14007.0D0
            ! -------------------------------------------------------------------------
            ! End of calculate the alkality gain by denitrification
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! END OF NITROGEN BASED SOURCES
            ! -------------------------------------------------------------------------
            
            
            ! -------------------------------------------------------------------------
            ! NITROGEN BASED SINKS
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! Calculate the alkalinity loss by ammonium consumption
            ! 1 eq alk for each ammonium consumed since one positive ion is lost
            ! -------------------------------------------------------------------------
            N_CHEM_AUT_BAC_GROWTH = 0.0D0
            N_DIA_GROWTH          = 0.0D0
            N_CYN_GROWTH          = 0.0D0
            N_OPA_GROWTH          = 0.0D0
            N_NON_FIX_CYN_GROWTH  = 0.0D0
            N_AER_HET_BAC_GROWTH  = 0.0D0
            
            ALK_LOST_BY_AMMONIUM_CONS  = &
                (N_CHEM_AUT_BAC_GROWTH + N_DIA_GROWTH + &
                 N_CYN_GROWTH          + N_OPA_GROWTH + &
                 N_NON_FIX_CYN_GROWTH  + N_AER_HET_BAC_GROWTH) / 14007.0D0
            ! -------------------------------------------------------------------------
            ! End of calculate the alkalinity loss by ammonium consumption
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! Calculate the alkalinity loss by nitrification
            ! - 2 eq alk for each ammonium nitrified since one positive ion is lost and
            !   one negative ion is gained
            !
            ! - 1 eq alk for each NH3 nitrified since one uncharged nitrogen is lost and
            !   one negative ion is gained
            ! -------------------------------------------------------------------------
            N_NITRIFICATION_NH4 = R_NITR * FRAC_NH4
            N_NITRIFICATION_NH3 = R_NITR * FRAC_NH3
            
            ALK_LOST_BY_NITRIFICATION = &
                ((2.0D0 * N_NITRIFICATION_NH4) + N_NITRIFICATION_NH3) / 14007.0D0
            ! -------------------------------------------------------------------------
            ! END OF NITROGEN BASED SINKS
            ! -------------------------------------------------------------------------
            
            
            ! -------------------------------------------------------------------------
            ! PREPARE FOR PHOSPHORUS-ALKALINITY COUPLING
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! Calculate the dissociation constants of H2PO3
            ! -------------------------------------------------------------------------
            H_PLUS(:,:) = 10.0D0 ** (-pH(:,:))
            
            KP_OPTION = 1
            
            select case (KP_OPTION)
                
                ! Option 0 : Use fixed values from the general chemistry book
                case (0)
                    K_ONE_TIP   = 10.0D0 ** (-2.15D0)
                    K_TWO_TIP   = 10.0D0 ** (-7.20D0)
                    K_THREE_TIP = 10.0D0 ** (-2.15D0)
            
                ! Option 1 : Use values calculated by CO2SYS. Programmed together with call to co2sys
                case (1)
                     continue
!                    K_ONE_TIP   = K_ONE_TIP  
!                    K_TWO_TIP   = K_TWO_TIP  
!                    K_THREE_TIP = K_THREE_TIP
            end select
            
            FRACTION_DIVISOR_TIP = &
                   (H_PLUS *    H_PLUS * H_PLUS) + &
                (K_ONE_TIP * K_ONE_TIP * H_PLUS) + &
                (K_ONE_TIP * K_TWO_TIP * H_PLUS) + &
                (K_ONE_TIP * K_TWO_TIP * K_THREE_TIP)
            
            ALPHA_H2PO4 = (K_ONE_TIP * H_PLUS    * H_PLUS     ) / FRACTION_DIVISOR_TIP
            ALPHA_HPO4  = (K_ONE_TIP * K_TWO_TIP * H_PLUS     ) / FRACTION_DIVISOR_TIP
            ALPHA_PO4   = (K_ONE_TIP * K_TWO_TIP * K_THREE_TIP) / FRACTION_DIVISOR_TIP
            
            PHOSPHATE_EQ_CONSTANT = &
                ALPHA_H2PO4 + (2.0D0 * ALPHA_HPO4) + (3.0D0 * ALPHA_HPO4)
            ! -------------------------------------------------------------------------
            ! End of calculate the dissociation constants of H2PO3
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! END OF PREPARE FOR PHOSPHORUS-ALKALINITY COUPLING
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! PHOSPHORUS BASED SOURCES
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! Calculate the alkality gain by phosphate consumption
            ! 1 eq alk for each H2PO4 consumed since one negative ion charge is lost
            ! 2 eq alk for each HPO4 consumed since two negative ion charges are lost
            ! 3 eq alk for each PO4 consumed since three negative ion charges are lost
            ! -------------------------------------------------------------------------
            P_CHEM_AUT_BAC_GROWTH = 0.0D0
            P_DIA_GROWTH          = 0.0D0
            P_CYN_GROWTH          = 0.0D0
            P_OPA_GROWTH          = 0.0D0
            P_NON_FIX_CYN_GROWTH  = 0.0D0
            P_AER_HET_BAC_GROWTH  = 0.0D0
            
            ALK_GAINED_BY_PHOSPHATE_CONS = &
                (P_CHEM_AUT_BAC_GROWTH + P_DIA_GROWTH + &
                 P_CYN_GROWTH          + P_OPA_GROWTH + &
                 P_NON_FIX_CYN_GROWTH  + P_AER_HET_BAC_GROWTH) / 30974.0D0
            ! -------------------------------------------------------------------------
            ! End of calculate the alkality gain by phosphate consumption
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! END OF PHOSPHORUS BASED SOURCES
            ! -------------------------------------------------------------------------
            
            
            ! -------------------------------------------------------------------------
            ! PHOSPHORUS BASED SINKS
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! Calculate the alkality loss by phosphate generation
            ! 1 eq alk for each H2PO4 consumed since one negative ion charge is lost
            ! 2 eq alk for each HPO4 consumed since two negative ion charges are lost
            ! 3 eq alk for each PO4 consumed since three negative ion charges are lost
            ! -------------------------------------------------------------------------
            P_CHEM_AUT_BAC_TOT_RESP   = 0.0D0
            P_AER_HET_BAC_INT_RESP    = 0.0D0
            P_FAC_AN_HET_BAC_TOT_RESP = 0.0D0
            P_DIA_TOT_RESP            = 0.0D0
            P_CYN_TOT_RESP            = 0.0D0
            P_OPA_TOT_RESP            = 0.0D0
            P_FIX_CYN_TOT_RESP        = 0.0D0
            P_AER_HET_BAC_P_OX        = 0.0D0
            P_FAC_AN_HET_BAC_P_OX     = 0.0D0
            P_ZOO_TOT_RESP            = 0.0D0
            P_ABIOTIC_DON_MIN         = R_MINER_DOP * PHOSPHATE_EQ_CONSTANT
            
            ALK_LOST_BY_PHOSPHATE_GEN = &
                (P_CHEM_AUT_BAC_TOT_RESP   + P_AER_HET_BAC_INT_RESP + &
                 P_FAC_AN_HET_BAC_TOT_RESP + P_DIA_TOT_RESP         + &
                 P_CYN_TOT_RESP            + P_OPA_TOT_RESP         + &
                 P_FIX_CYN_TOT_RESP        + P_AER_HET_BAC_P_OX     + &
                 P_FAC_AN_HET_BAC_P_OX     + P_ZOO_TOT_RESP         + &
                 P_ABIOTIC_DON_MIN) / 30974.0D0
            ! -------------------------------------------------------------------------
            ! End of calculate the alkality loss by phosphate generation
            ! -------------------------------------------------------------------------
            
            ! -------------------------------------------------------------------------
            ! END OF PHOSPHORUS BASED SINKS
            ! -------------------------------------------------------------------------
            
            
            
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! TO DO : Processes that should be added for nutrients
            !         - Precipitation / settling
            !         - Sediment fluxes
            !         - Adsorption
            !         - Desorption
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            
            ! Calculate the total derivative and convert in moles
            if (CONSIDER_ALKALNITY_DERIVATIVE > 0) then
                ALK_KINETIC_DERIVATIVE = &
                    ALK_GAINED_BY_AMMONIUM_GEN   + &
                    ALK_GAINED_BY_NITRATE_CONS   + &
                    ALK_GAINED_BY_PHOSPHATE_CONS - &
                    ALK_LOST_BY_AMMONIUM_CONS    - &
                    ALK_LOST_BY_NITRIFICATION    - &
                    ALK_LOST_BY_PHOSPHATE_GEN   
            else
                ALK_KINETIC_DERIVATIVE       = 0.0D0
                ALK_GAINED_BY_AMMONIUM_GEN   = 0.0D0
                ALK_GAINED_BY_NITRATE_CONS   = 0.0D0
                ALK_GAINED_BY_PHOSPHATE_CONS = 0.0D0
                ALK_LOST_BY_AMMONIUM_CONS    = 0.0D0
                ALK_LOST_BY_NITRIFICATION    = 0.0D0
                ALK_LOST_BY_PHOSPHATE_GEN    = 0.0D0
            end if
            ! -------------------------------------------------------------------------
            ! END OF KINETIC SUBMODEL FOR ALKALINITY
            ! -------------------------------------------------------------------------
            
            KINETIC_DERIVS(:,:, 14)   = ALK_KINETIC_DERIVATIVE
                
            PROCESSES_sed(:,:, 14, 1) = pH(:,:)  !pH
            PROCESSES_sed(:,:, 14, 2) = ALK_GAINED_BY_AMMONIUM_GEN  
            PROCESSES_sed(:,:, 14, 3) = ALK_GAINED_BY_NITRATE_CONS  
            PROCESSES_sed(:,:, 14, 4) = ALK_GAINED_BY_PHOSPHATE_CONS
            PROCESSES_sed(:,:, 14, 5) = ALK_LOST_BY_AMMONIUM_CONS   
            PROCESSES_sed(:,:, 14, 6) = ALK_LOST_BY_NITRIFICATION   
            PROCESSES_sed(:,:, 14, 7) = ALK_LOST_BY_PHOSPHATE_GEN   
            
                      
            !end KINETIC DERIVATIVES

            
            !KINETIC_DERIVS(SED_NH4N) = R_MINER_DON - R_NITR
            PROCESSES_sed(:,:, 1, 10) = R_MINER_DON
            PROCESSES_sed(:,:, 1, 11) = R_NITR     
            
            !KINETIC_DERIVS(SED_NO3N) = R_NITR - R_DENITR
            PROCESSES_sed(:,:, 2, 10) = R_MINER_DON
            PROCESSES_sed(:,:, 2, 11) = R_NITR     
            
            !KINETIC_DERIVS(SED_DON) = R_DISS_PON - R_MINER_DON
            PROCESSES_sed(:,:, 3, 10) = R_DISS_PON 
            PROCESSES_sed(:,:, 3, 11) = R_MINER_DON
            
            !KINETIC_DERIVS(SED_PON) = (-1.0D0) * R_DISS_PON 
            PROCESSES_sed(:,:, 4, 10) = R_DISS_PON
            
            !KINETIC_DERIVS(SED_PO4P) = R_MINER_DOP 
            PROCESSES_sed(:,:, 5, 10) = R_MINER_DOP
            
            !KINETIC_DERIVS(SED_DOP) = R_DISS_POP - R_MINER_DOP
            PROCESSES_sed(:,:, 6, 10) = R_DISS_POP 
            PROCESSES_sed(:,:, 6, 11) = R_MINER_DOP
            
            !KINETIC_DERIVS(SED_POP) = (-1.0D0) * R_DISS_POP 
            PROCESSES_sed(:,:, 7, 10) = R_DISS_POP
            
            !KINETIC_DERIVS(SED_DOXY) = ((-4.57D0) * R_NITR) - DEOXYGENATION
            PROCESSES_sed(:,:, 8, 10) = 4.57D0 * R_NITR
            PROCESSES_sed(:,:, 8, 11) = DEOXYGENATION  
            
            !KINETIC_DERIVS(SED_DOC) = R_DISS_POC - R_MINER_DOC - (R_DENITR / DENITR_YIELD)
            PROCESSES_sed(:,:, 9, 10) = R_DISS_POC 
            PROCESSES_sed(:,:, 9, 11) = R_MINER_DOC
            PROCESSES_sed(:,:, 9, 12) = (R_DENITR / DENITR_YIELD)
            
            !KINETIC_DERIVS(SED_POC) = (-1.0D0) * R_DISS_POC 
            PROCESSES_sed(:,:, 10, 10) = R_DISS_POC
            
            !KINETIC_DERIVS(SED_DSi) = R_DISS_PSi
            PROCESSES_sed(:,:, 11, 10) = R_DISS_PSi
            
            !KINETIC_DERIVS(SED_DSi) = R_DISS_PSi
            PROCESSES_sed(:,:, 11, 10) = R_DISS_PSi
            
            !KINETIC_DERIVS(SED_PSi) = (-1.0D0) * R_DISS_PSi
            PROCESSES_sed(:,:, 12, 10) = R_DISS_PSi 
                   
        end if ! switch_kinetics
        
        !Deposition(settling) derivatives:      
        !24/06/2012: If clause Added by Ali for safer switch off
        if (switch_settling.ne.0) then
            
            
            SETTLING_AFFECTED_DEPTH(:) = sum(SED_DEPTHS(:,1:NUM_FLUX_RECEIVING_SED_LAYERS),2)
            
            do i = 1, NUM_FLUX_RECEIVING_SED_LAYERS
                do j = 1, NUM_SED_VARS
                    SETTLING_DERIVS(:,i, j)  = FLUXES_TO_SEDIMENTS(:,j) * &
                        (SED_DEPTHS(:,i) / SETTLING_AFFECTED_DEPTH(:))
            
                    PROCESSES_sed(:,i, j, 7) = FLUXES_TO_SEDIMENTS(:,j) * &
                        (SED_DEPTHS(:,i) / SETTLING_AFFECTED_DEPTH(:))
                end do
            end do
            !End deposition fluxes      
        
        end if
        
        !THE FOLLOWING COMMENTED LINES ARE TO BE DELETED WHEN SAFE SWITCH
        !OFF ALGORITHMS INTRODUECED BY ALI AT 24/06/2012 IN PETRAS HOUSE ARE
        !COMPLETED AND TESTED
        !Switch to zero kinetic and settling derivatives if necessary      
        !if (switch_kinetics .eq. 0) KINETIC_DERIVS(:, :) = 0. 
        !if (switch_settling .eq. 0) SETTLING_DERIVS(:, :) = 0.
        
         do i=1,nkn
         index_strange(i) = i
         end do
         


        !Time integration
        do i = 1, NUM_SED_LAYERS
        
            do j = 1, NUM_SED_VARS  
             
               if(debug_stranger) then
                if(STRANGERSD(TRANSPORT_DERIVS(:,i, j),VALUE_strange,nkn) .eq. 1) then
                   nstrange = count(VALUE_strange)
                   allocate(STRANGERS    (nstrange))
                   allocate(NODES_STRANGE(nstrange))
             
                   l=1
                   do k=1,nkn
                    if(VALUE_strange(k)) then
                      STRANGERS    (l) = TRANSPORT_DERIVS(k,i,j)
                      NODES_STRANGE(l) = index_strange(k)
                      l=l+1
                    end if
                   end do
                       
                    print *, 'AQUABC_SEDIMENT1: Layer ', i, 'Variable ',j  !,'Cell ',k 
                    print *, 'Transport derivative is NaN'
                    
                    print *, 'NODE_NUMBERS=',NODES_STRANGE
                    print *, 'VALUES=',STRANGERS
                    
                     print *, 'DERIV(:,I,J)                  : ', TRANSPORT_DERIVS  (1:10,i, j)
                     print *, '  BURIAL_DERIVS     (:,i, j)  : ', BURIAL_DERIVS     (1:10,i, j)
                     print *, '  PART_MIXING_DERIVS(:,i, j)  : ', PART_MIXING_DERIVS(1:10,i, j)
                     print *, '  DIFFUSION_DERIVS  (:,i, j)  : ', DIFFUSION_DERIVS  (1:10,i, j)
                     print *, '  BURIAL_DERIVS     (:,i, j)  : ', BURIAL_DERIVS     (1:10,i, j)
                     print *, '  PART_MIXING_DERIVS(:,i, j)  : ', PART_MIXING_DERIVS(1:10,i, j)
                    
                    if (I.lt.NUM_SED_LAYERS) then
                        print *, '  PART_MIXING_RATES  (1:10,i + 1, j) : ', &
                                 PART_MIXING_RATES     (1:10,i + 1, j)
                    end if
                    
                    print *, '  PART_MIXING_RATES  (1:10,i, j) : ',PART_MIXING_RATES  (1:10,i, j)
                    print *, '  ADVECTION_DERIVS   (1:10,i, j) : ', ADVECTION_DERIVS (1:10,i, j)
                    error =1
                    deallocate(STRANGERS, NODES_STRANGE)
                end if
               end if
               
               if(debug_stranger) then 
                if(STRANGERSD(KINETIC_DERIVS(:,i, j),VALUE_strange,nkn) .eq. 1) then
                   nstrange = count(VALUE_strange)
                   allocate(STRANGERS    (nstrange))
                   allocate(NODES_STRANGE(nstrange))
             
                   l=1
                   do k=1,nkn
                    if(VALUE_strange(k)) then
                      STRANGERS    (l) = KINETIC_DERIVS(k,i,j)
                      NODES_STRANGE(l) = index_strange(k)
                      l=l+1
                    end if
                   end do
                   
                   print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j !, 'Cell ',CELLNO                  
                   print *, 'Kinetic derivative is NaN'                    
                   print *, 'NODE_NUMBERS=',NODES_STRANGE
                   print *, 'VALUES=',STRANGERS

                  !  print *, 'Deriv(i,j)=',KINETIC_DERIVS(:,i, j)
                  
                  deallocate(STRANGERS, NODES_STRANGE)
                  error =1
                end if
               end if
               
               if(debug_stranger) then  
                if(STRANGERSD(SETTLING_DERIVS(:,i, j),VALUE_strange,nkn) .eq. 1) then
                   nstrange = count(VALUE_strange)
                   allocate(STRANGERS    (nstrange))
                   allocate(NODES_STRANGE(nstrange))
             
                   l=1
                   do k=1,nkn
                    if(VALUE_strange(k)) then
                      STRANGERS    (l) = SETTLING_DERIVS(k,i,j)
                      NODES_STRANGE(l) = index_strange(k)
                      l=l+1
                    end if
                   end do
                   print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j !, 'Cell ',CELLNO  
                   print *, 'Settling derivative is NaN'
                   print *, 'NODE_NUMBERS=',NODES_STRANGE 
                   print *, 'VALUES=',STRANGERS           
!                    print *, 'Deriv(i,j)=',SETTLING_DERIVS(:,i, j)

                    error =1
                    deallocate(STRANGERS, NODES_STRANGE)
                end if
               end if 
        
                !sum of all derivatives                   
                DERIVS(:,i, j) = TRANSPORT_DERIVS(:,i, j) + KINETIC_DERIVS(:,i, j) + &
                                 SETTLING_DERIVS (:,i, j)
                        
                 PROCESSES_sed(:,i,j, 15) =  BURIAL_DERIVS     (:,i, j) 
                 PROCESSES_sed(:,i,j, 16) =  PART_MIXING_DERIVS(:,i, j)
                 PROCESSES_sed(:,i,j, 17) =  DIFFUSION_DERIVS  (:,i, j)
                 PROCESSES_sed(:,i,j, 18) =  ADVECTION_DERIVS  (:,i, j)
                 PROCESSES_sed(:,i,j, 19) =  TRANSPORT_DERIVS  (:,i, j)
                 PROCESSES_sed(:,i,j, 20) =  KINETIC_DERIVS    (:,i, j)
                 PROCESSES_sed(:,i,j, 21) =  SETTLING_DERIVS   (:,i, j)
                 PROCESSES_sed(:,i,j, 22) =  DERIVS            (:,i, j)
                 if(i .eq. 1) then 
                  PROCESSES_sed(:,i,j, 23) =  FLUXES_TO_SEDIMENTS(:,j)
                 else
                  PROCESSES_sed(:,i,j, 23) = 0.D0
                 end if
                 
                 ! Number 7 is reserved for settling derivatives (done in settling section)
                 ! Number 8 is reserved for flux to BS (done in aquabc)
                 ! Number 9 is reserved for flux from BS (done in aquabc)                  

        
                if (IN_WHICH_PHASE(j).eq.0) then
                    VOLUME_FRACTION = SED_POROSITIES(:,i)
                end if
                
                if (IN_WHICH_PHASE(j).eq.1) then
                    VOLUME_FRACTION = 1.0D0 - SED_POROSITIES(:,i)
                end if
        
                if (IN_WHICH_PHASE(j).eq.2) then
                    VOLUME_FRACTION = 1.0D0
                end if
        
                
                UNIT_AREA_MASSES(:,i, j) = UNIT_AREA_MASSES(:,i, j) + &
                     (DERIVS(:,i, j) * (TIME_STEP / NUM_SUB_TIME_STEPS))
        
                INTERMED_RESULTS(:,i, j) = UNIT_AREA_MASSES(:,i, j) / &
                    (SED_DEPTHS(:,i) * VOLUME_FRACTION)
                    
                    
            if(debug_stranger) then        
             if(STRANGERSD(INTERMED_RESULTS(:,i, j),VALUE_strange,nkn) .eq. 1) then
             
                    nstrange = count(VALUE_strange)
                    allocate(STRANGERS    (nstrange))
                    allocate(NODES_STRANGE(nstrange))
              
                    l=1
                    do k=1,nkn
                     if(VALUE_strange(k)) then
                       STRANGERS    (l) = INTERMED_RESULTS(k,i,j)
                       NODES_STRANGE(l) = index_strange(k)
                       l=l+1
                     end if
                    end do
                    print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j
                    print *, 'Intermediate result on substeps is NaN'
                    print *, 'NODE_NUMBERS=',NODES_STRANGE
                    print *, 'VALUES=',STRANGERS 
             
                    deallocate(STRANGERS, NODES_STRANGE)
             
!                 print *, 'FINAL(i,j)=',FINAL_SED_STATE_VARS(i, j)
!                 print *, 'time= ', PSTIME
                 stop
             end if  ! strangers checking
            end if                       
                    
            end do ! j-variables
        
        end do     ! i-layers 
        !End time integration        
        
        if (error .eq. 1) then
          print *, 'time= ', PSTIME 
          stop
        endif
    
    end do !end OF SUB TIME STEPS : CAUTION, FLUXES ARE NOT AVERAGED YET. 
           !fixme if substeps are necessary            
           
           
           
!======================================================================================           
           
           
           
         
    !Final assignments for state var. and var. for the output      
    do i = 1, NUM_SED_LAYERS
    
        do j = 1, NUM_SED_VARS
            FINAL_SED_STATE_VARS(:,i, j) = INTERMED_RESULTS(:,i, j)
            
           if(debug_stranger) then 
            if(STRANGERSD(FINAL_SED_STATE_VARS(:,i, j),VALUE_strange,nkn) .eq. 1) then
            
                   nstrange = count(VALUE_strange)
                   allocate(STRANGERS    (nstrange))
                   allocate(NODES_STRANGE(nstrange))
             
                   l=1
                   do k=1,nkn
                    if(VALUE_strange(k)) then
                      STRANGERS    (l) = FINAL_SED_STATE_VARS(k,i,j)
                      NODES_STRANGE(l) = index_strange(k)
                      l=l+1
                    end if
                   end do
                   print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j
                   print *, 'Final state is NaN'
                   print *, 'NODE_NUMBERS=',NODES_STRANGE
                   print *, 'VALUES=',STRANGERS 
            
                   deallocate(STRANGERS, NODES_STRANGE)

!                print *, 'FINAL(i,j)=',FINAL_SED_STATE_VARS(i, j)
!                print *, 'time= ', PSTIME
                stop
            end if 
           end if  ! strangers     
    
            SED_OUTPUTS(:,i, j) = INTERMED_RESULTS(:,i, j) ! State variables
        end do !j - variables
    
        SED_OUTPUTS(:,I, NUM_SED_VARS + 1) = SOLUTE_FRACTIONS(:,I, 1) * INTERMED_RESULTS(:,I, 1) !Dissolved  NH4N
        SED_OUTPUTS(:,I, NUM_SED_VARS + 2) = SOLUTE_FRACTIONS(:,I, 5) * INTERMED_RESULTS(:,I, 5) !Dissolved   PO4P
        

    end do !i - layers
    
    i=1  ! just for breakpoint in case of debugging    
            
    end ! end of sediment routine
    
    !******************************************************************* 
    !*******************************************************************
        
      INTEGER FUNCTION STRANGER(VALUE)
      
! cheks for NaN and Inf in double precision

      DOUBLE PRECISION VALUE, BIGNUMBER, RATIO, LLIMIT, ULIMIT
      LLIMIT = -1.0D4
      ULIMIT =  1.0D4
      BIGNUMBER=1.0D300
      STRANGER=0
      
      if (isnan(VALUE)) then
       STRANGER=1
      end if
      
      if (isnan(VALUE)) then
       STRANGER=1
      end if
      
      RATIO = BIGNUMBER/VALUE
      
      if(RATIO .eq. 0.D0) then
       STRANGER=1
      end if
      
      if(VALUE .le. LLIMIT .or. VALUE .ge. ULIMIT) then      
       STRANGER=1
      end if

      return
      end
         
    !************************************************************************
    !************************************************************************
      
subroutine FLX_ALUKAS_II_TO_SED_MOD_1 &
           (STATE_VARIABLES    , NUM_VARS   ,         &
            MODEL_CONSTANTS    , NUM_CONSTS ,         &
            DRIVING_FUNCTIONS  , NUM_DRIV   ,         &
            SETTLING_VELOCITIES, DISSOLVED_FRACTIONS, &
            BOTTOM_FACTOR , CELLNO, PSTIME,           &
            SETTLING_RATES, FLUXES, NUM_FLUXES,       &
            SEDIMENT_TYPE , FRACTION_OF_DEPOSITION,   & 
            NOT_DEPOSITED_FLUXES, NUM_NOT_DEPOSITED_FLUXES)
            
! Outputs:
!  FLUXES               - for each BS state variable
!  SETTLING_RATES       - in g/m2/day
!  NOT_DEPOSITED_FLUXES - for each WC state variable

    use AQUABC_II_GLOBAL
    use para_aqua
    implicit none
    
    integer NUM_VARS
    integer NUM_CONSTS
    integer NUM_FLUXES
    integer NUM_DRIV
    
    double precision STATE_VARIABLES
    DIMENSION STATE_VARIABLES(NUM_VARS)
    
    double precision MODEL_CONSTANTS
    DIMENSION MODEL_CONSTANTS(NUM_CONSTS)
    
    double precision SETTLING_VELOCITIES
    DIMENSION SETTLING_VELOCITIES(NUM_VARS)
    
    double precision DISSOLVED_FRACTIONS
    DIMENSION DISSOLVED_FRACTIONS(NUM_VARS)
    
    double precision SETTLING_RATES
    DIMENSION SETTLING_RATES(NUM_VARS)
    
    double precision FLUXES
    DIMENSION FLUXES(NUM_VARS)
    
    double precision DRIVING_FUNCTIONS
    DIMENSION DRIVING_FUNCTIONS(NUM_DRIV)
    
    integer CELLNO
    integer LAYER
    
    double precision PSTIME
    double precision BOTTOM_FACTOR
    
    integer SEDIMENT_TYPE
    integer NUM_NOT_DEPOSITED_FLUXES
    double precision FRACTION_OF_DEPOSITION(NUM_NOT_DEPOSITED_FLUXES)
    double precision NOT_DEPOSITED_FLUXES(NUM_NOT_DEPOSITED_FLUXES)
    
    double precision NH4_FLUX
    double precision NO3_FLUX
    double precision PO4_FLUX
    double precision DISS_OXYGEN_FLUX
    double precision NITR_BAC_C_FLUX
    double precision AER_HET_BAC_C_FLUX
    double precision DENITR_BAC_C_FLUX
    double precision DIA_C_FLUX
    double precision ZOO_C_FLUX
    double precision ZOO_N_FLUX
    double precision ZOO_P_FLUX
    double precision DET_PART_ORG_C_FLUX
    double precision DET_PART_ORG_N_FLUX
    double precision DET_PART_ORG_P_FLUX
    double precision DISS_ORG_C_FLUX
    double precision DISS_ORG_N_FLUX
    double precision DISS_ORG_P_FLUX
    double precision CYN_C_FLUX
    double precision OPA_C_FLUX
    double precision DISS_Si_FLUX
    double precision BIOG_Si_FLUX
    double precision FIX_CYN_C_FLUX
    real(kind = DBL_PREC) :: INORG_C_FLUX
    real(kind = DBL_PREC) :: TOT_ALK_FLUX
    
    double precision SETTLING_FACTORS
    DIMENSION SETTLING_FACTORS(NUM_VARS)
    
    double precision DIA_N_TO_C        
    double precision CYN_N_TO_C        
    double precision OPA_N_TO_C        
    double precision FIX_CYN_N_TO_C
    double precision NITR_BAC_N_TO_C   
    double precision AER_HET_BAC_N_TO_C
    double precision DENITR_BAC_N_TO_C 
    
    double precision DIA_P_TO_C        
    double precision CYN_P_TO_C        
    double precision OPA_P_TO_C
    double precision FIX_CYN_P_TO_C        
    double precision NITR_BAC_P_TO_C   
    double precision AER_HET_BAC_P_TO_C
    double precision DENITR_BAC_P_TO_C 
    
    double precision DIA_Si_TO_C
    
    integer I
! Change these lines if numbers of constants was changed    
!     NITR_BAC_N_TO_C     = MODEL_CONSTANTS(16)
!     NITR_BAC_P_TO_C     = MODEL_CONSTANTS(17)    
!     AER_HET_BAC_N_TO_C  = MODEL_CONSTANTS(42)
!     AER_HET_BAC_P_TO_C  = MODEL_CONSTANTS(43)    
!     DENITR_BAC_N_TO_C   = MODEL_CONSTANTS(60)
!     DENITR_BAC_P_TO_C   = MODEL_CONSTANTS(61)    
!     DIA_N_TO_C          = MODEL_CONSTANTS(83)
!     DIA_P_TO_C          = MODEL_CONSTANTS(84)    
!     CYN_N_TO_C          = MODEL_CONSTANTS(106)
!     CYN_P_TO_C          = MODEL_CONSTANTS(107)    
!     FIX_CYN_N_TO_C      = MODEL_CONSTANTS(128)
!     FIX_CYN_P_TO_C      = MODEL_CONSTANTS(129)    
!     OPA_N_TO_C          = MODEL_CONSTANTS(152)
!     OPA_P_TO_C          = MODEL_CONSTANTS(153)                    
!     DIA_Si_TO_C         = MODEL_CONSTANTS(85)
    
     call para_get_value( 'CHEM_AUT_BAC_N_TO_C'   ,    NITR_BAC_N_TO_C    )
     call para_get_value( 'CHEM_AUT_BAC_P_TO_C'   ,    NITR_BAC_P_TO_C    )
     call para_get_value( 'AER_HET_BAC_N_TO_C',    AER_HET_BAC_N_TO_C )
     call para_get_value( 'AER_HET_BAC_P_TO_C',    AER_HET_BAC_P_TO_C )
     call para_get_value( 'FAC_AN_HET_BAC_N_TO_C' ,    DENITR_BAC_N_TO_C  )
     call para_get_value( 'FAC_AN_HET_BAC_P_TO_C' ,    DENITR_BAC_P_TO_C  )
     call para_get_value( 'DIA_N_TO_C'        ,    DIA_N_TO_C         )
     call para_get_value( 'DIA_P_TO_C'        ,    DIA_P_TO_C         )
     call para_get_value( 'CYN_N_TO_C'        ,    CYN_N_TO_C         )
     call para_get_value( 'CYN_P_TO_C'        ,    CYN_P_TO_C         )
     call para_get_value( 'FIX_CYN_N_TO_C'    ,    FIX_CYN_N_TO_C     )
     call para_get_value( 'FIX_CYN_P_TO_C'    ,    FIX_CYN_P_TO_C     )
     call para_get_value( 'OPA_N_TO_C'        ,    OPA_N_TO_C         )
     call para_get_value( 'OPA_P_TO_C'        ,    OPA_P_TO_C         )
     call para_get_value( 'DIA_Si_TO_C'       ,    DIA_Si_TO_C        )
    
    do i = 1, NUM_VARS
        SETTLING_FACTORS(i) = BOTTOM_FACTOR * &
                              (1.0D+0 - DISSOLVED_FRACTIONS(i))* &
                              SETTLING_VELOCITIES(i)
    
        if (SETTLING_FACTORS(i).lt.0.0D0) then
            SETTLING_FACTORS(i) = 0.0D0
        end if
    
        SETTLING_RATES(i)   = STATE_VARIABLES(i) *  &
                              (1.0D+0 - DISSOLVED_FRACTIONS(i))* &
                              SETTLING_VELOCITIES(i)
    end do
    
    NH4_FLUX            = STATE_VARIABLES(1)  * SETTLING_FACTORS(1)
    NO3_FLUX            = STATE_VARIABLES(2)  * SETTLING_FACTORS(2)
    PO4_FLUX            = STATE_VARIABLES(3)  * SETTLING_FACTORS(3)
    DISS_OXYGEN_FLUX    = STATE_VARIABLES(4)  * SETTLING_FACTORS(4)
    NITR_BAC_C_FLUX     = STATE_VARIABLES(5)  * SETTLING_FACTORS(5)
    AER_HET_BAC_C_FLUX  = STATE_VARIABLES(6)  * SETTLING_FACTORS(6)
    DENITR_BAC_C_FLUX   = STATE_VARIABLES(7)  * SETTLING_FACTORS(7)
    DIA_C_FLUX          = STATE_VARIABLES(8)  * SETTLING_FACTORS(8)
    ZOO_C_FLUX          = STATE_VARIABLES(9)  * SETTLING_FACTORS(9)
    ZOO_N_FLUX          = STATE_VARIABLES(10) * SETTLING_FACTORS(10)
    ZOO_P_FLUX          = STATE_VARIABLES(11) * SETTLING_FACTORS(11)
    DET_PART_ORG_C_FLUX = STATE_VARIABLES(12) * SETTLING_FACTORS(12)
    DET_PART_ORG_N_FLUX = STATE_VARIABLES(13) * SETTLING_FACTORS(13)
    DET_PART_ORG_P_FLUX = STATE_VARIABLES(14) * SETTLING_FACTORS(14)
    DISS_ORG_C_FLUX     = STATE_VARIABLES(15) * SETTLING_FACTORS(15)
    DISS_ORG_N_FLUX     = STATE_VARIABLES(16) * SETTLING_FACTORS(16)
    DISS_ORG_P_FLUX     = STATE_VARIABLES(17) * SETTLING_FACTORS(17)
    CYN_C_FLUX          = STATE_VARIABLES(18) * SETTLING_FACTORS(18)
    OPA_C_FLUX          = STATE_VARIABLES(19) * SETTLING_FACTORS(19)
    DISS_Si_FLUX        = STATE_VARIABLES(20) * SETTLING_FACTORS(20)
    BIOG_Si_FLUX        = STATE_VARIABLES(21) * SETTLING_FACTORS(21)
    FIX_CYN_C_FLUX      = STATE_VARIABLES(22) * SETTLING_FACTORS(22)
    INORG_C_FLUX        = STATE_VARIABLES(23) * SETTLING_FACTORS(23)
    TOT_ALK_FLUX        = STATE_VARIABLES(24) * SETTLING_FACTORS(24)
    
    !NOT DEPOSITED FLUXES      
    NOT_DEPOSITED_FLUXES(1)  = NH4_FLUX            * (1.0D0 - FRACTION_OF_DEPOSITION(1))
    NOT_DEPOSITED_FLUXES(2)  = NO3_FLUX            * (1.0D0 - FRACTION_OF_DEPOSITION(2))
    NOT_DEPOSITED_FLUXES(3)  = PO4_FLUX            * (1.0D0 - FRACTION_OF_DEPOSITION(3))
    NOT_DEPOSITED_FLUXES(4)  = DISS_OXYGEN_FLUX    * (1.0D0 - FRACTION_OF_DEPOSITION(4))
    NOT_DEPOSITED_FLUXES(5)  = NITR_BAC_C_FLUX     * (1.0D0 - FRACTION_OF_DEPOSITION(5))
    NOT_DEPOSITED_FLUXES(6)  = AER_HET_BAC_C_FLUX  * (1.0D0 - FRACTION_OF_DEPOSITION(6))
    NOT_DEPOSITED_FLUXES(7)  = DENITR_BAC_C_FLUX   * (1.0D0 - FRACTION_OF_DEPOSITION(7))
    NOT_DEPOSITED_FLUXES(8)  = DIA_C_FLUX          * (1.0D0 - FRACTION_OF_DEPOSITION(8))
    NOT_DEPOSITED_FLUXES(9)  = ZOO_C_FLUX          * (1.0D0 - FRACTION_OF_DEPOSITION(9))
    NOT_DEPOSITED_FLUXES(10) = ZOO_N_FLUX          * (1.0D0 - FRACTION_OF_DEPOSITION(10))
    NOT_DEPOSITED_FLUXES(11) = ZOO_P_FLUX          * (1.0D0 - FRACTION_OF_DEPOSITION(11))
    NOT_DEPOSITED_FLUXES(12) = DET_PART_ORG_C_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(12))
    NOT_DEPOSITED_FLUXES(13) = DET_PART_ORG_N_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(13))
    NOT_DEPOSITED_FLUXES(14) = DET_PART_ORG_P_FLUX * (1.0D0 - FRACTION_OF_DEPOSITION(14))
    NOT_DEPOSITED_FLUXES(15) = DISS_ORG_C_FLUX     * (1.0D0 - FRACTION_OF_DEPOSITION(15))
    NOT_DEPOSITED_FLUXES(16) = DISS_ORG_N_FLUX     * (1.0D0 - FRACTION_OF_DEPOSITION(16))
    NOT_DEPOSITED_FLUXES(17) = DISS_ORG_P_FLUX     * (1.0D0 - FRACTION_OF_DEPOSITION(17))
    NOT_DEPOSITED_FLUXES(18) = CYN_C_FLUX          * (1.0D0 - FRACTION_OF_DEPOSITION(18))
    NOT_DEPOSITED_FLUXES(19) = OPA_C_FLUX          * (1.0D0 - FRACTION_OF_DEPOSITION(19))
    NOT_DEPOSITED_FLUXES(20) = DISS_Si_FLUX        * (1.0D0 - FRACTION_OF_DEPOSITION(20))
    NOT_DEPOSITED_FLUXES(21) = BIOG_Si_FLUX        * (1.0D0 - FRACTION_OF_DEPOSITION(21))
    NOT_DEPOSITED_FLUXES(22) = FIX_CYN_C_FLUX      * (1.0D0 - FRACTION_OF_DEPOSITION(22))
    NOT_DEPOSITED_FLUXES(23) = INORG_C_FLUX        * (1.0D0 - FRACTION_OF_DEPOSITION(23))
    NOT_DEPOSITED_FLUXES(24) = TOT_ALK_FLUX        * (1.0D0 - FRACTION_OF_DEPOSITION(24))
    
    !AMMONIA NITROGEN FLUX
    FLUXES(1) = NH4_FLUX * FRACTION_OF_DEPOSITION(1)
    
    !NITRATE NITROGEN FLUX
    FLUXES(2) = NO3_FLUX * FRACTION_OF_DEPOSITION(2)
    
    !DISSOLVED ORGANIC NITROGEN FLUX
    FLUXES(3) = DISS_ORG_N_FLUX * FRACTION_OF_DEPOSITION(16)
    
    !PARTICULATE ORGANIC NITROGEN FLUX
    FLUXES(4) = &
       (DIA_C_FLUX          * FRACTION_OF_DEPOSITION(8)  * DIA_N_TO_C)        + &
       (CYN_C_FLUX          * FRACTION_OF_DEPOSITION(18) * CYN_N_TO_C)        + &
       (OPA_C_FLUX          * FRACTION_OF_DEPOSITION(19) * OPA_N_TO_C)        + &
       (FIX_CYN_C_FLUX      * FRACTION_OF_DEPOSITION(22) * FIX_CYN_N_TO_C)    + &
       (NITR_BAC_C_FLUX     * FRACTION_OF_DEPOSITION(5)  * NITR_BAC_N_TO_C)   + &
       (AER_HET_BAC_C_FLUX  * FRACTION_OF_DEPOSITION(6)  * AER_HET_BAC_N_TO_C)+ &
       (DENITR_BAC_C_FLUX   * FRACTION_OF_DEPOSITION(7)  * DENITR_BAC_N_TO_C) + &
       (ZOO_N_FLUX          * FRACTION_OF_DEPOSITION(10))                     + &
       (DET_PART_ORG_N_FLUX * FRACTION_OF_DEPOSITION(13))
         
    !PHOSPHATE FLUX
    FLUXES(5) = PO4_FLUX * FRACTION_OF_DEPOSITION(3)
    
    !DISSOLVED ORGANIC PHOSPHORUS FLUX
    FLUXES(6) = DISS_ORG_P_FLUX * FRACTION_OF_DEPOSITION(17)
    
    !PARTICULATE ORGANIC PHOSPHORUS FLUX
    FLUXES(7) = &
       (DIA_C_FLUX          * FRACTION_OF_DEPOSITION(8)  * DIA_P_TO_C)        + &
       (CYN_C_FLUX          * FRACTION_OF_DEPOSITION(18) * CYN_P_TO_C)        + &
       (OPA_C_FLUX          * FRACTION_OF_DEPOSITION(19) * OPA_P_TO_C)        + &
       (FIX_CYN_C_FLUX      * FRACTION_OF_DEPOSITION(22) * FIX_CYN_P_TO_C)    + &
       (NITR_BAC_C_FLUX     * FRACTION_OF_DEPOSITION(5)  * NITR_BAC_P_TO_C)   + &
       (AER_HET_BAC_C_FLUX  * FRACTION_OF_DEPOSITION(6)  * AER_HET_BAC_P_TO_C)+ &
       (DENITR_BAC_C_FLUX   * FRACTION_OF_DEPOSITION(7)  * DENITR_BAC_P_TO_C) + &
       (ZOO_P_FLUX          * FRACTION_OF_DEPOSITION(11)) +                     &
       (DET_PART_ORG_P_FLUX * FRACTION_OF_DEPOSITION(14))
    
    !DISSOLVED OXYGEN
    FLUXES(8) = 0.0D0
    
    !DISSOLVED ORGANIC CARBON
    FLUXES(9) = DISS_ORG_C_FLUX * FRACTION_OF_DEPOSITION(15)
    
    !PARTICULATE ORGANIC CARBON FLUX
    FLUXES(10) = &
       (DIA_C_FLUX          * FRACTION_OF_DEPOSITION(8))  + &
       (CYN_C_FLUX          * FRACTION_OF_DEPOSITION(18)) + &
       (OPA_C_FLUX          * FRACTION_OF_DEPOSITION(19)) + &
       (FIX_CYN_C_FLUX      * FRACTION_OF_DEPOSITION(22)) + &
       (NITR_BAC_C_FLUX     * FRACTION_OF_DEPOSITION(5))  + &
       (AER_HET_BAC_C_FLUX  * FRACTION_OF_DEPOSITION(6))  + &
       (DENITR_BAC_C_FLUX   * FRACTION_OF_DEPOSITION(7))  + &
       (ZOO_C_FLUX          * FRACTION_OF_DEPOSITION(9))  + &
       (DET_PART_ORG_C_FLUX * FRACTION_OF_DEPOSITION(12))
    
    !DISSOLVED SILICON
    FLUXES(11) = DISS_Si_FLUX * FRACTION_OF_DEPOSITION(20)
    
    !PARTICULATE  SILICON FLUX
    FLUXES(12) = &
       (DIA_C_FLUX   * FRACTION_OF_DEPOSITION(8) * DIA_Si_TO_C) + &
       (BIOG_Si_FLUX * FRACTION_OF_DEPOSITION(21))
    
    !  print *,'FRACTIONS_OF_DEPOSITION:', FRACTION_OF_DEPOSITION
    !  print *,'FLUXES:', FLUXES
    !  stop

    !INORGANIC CARBON FLUX
    FLUXES(13) = INORG_C_FLUX * FRACTION_OF_DEPOSITION(13)

    !ALKALINITY FLUX
    FLUXES(14) = TOT_ALK_FLUX * FRACTION_OF_DEPOSITION(14)
    
    !Salinity flux
    FLUXES(15) = 0.D1
     
end subroutine FLX_ALUKAS_II_TO_SED_MOD_1
      
!********************************************************************
!********************************************************************
    
    
!THIS IS A USER PROGRAMMED FUNCTION. IT IS THE USERS RESPONSIBILITY
!TO PROGRAM THE EQUATIONS WHICH GIVE THE MOLECULAR DIFFUSION COEFFICIENTS
!OF STATE VARIABLES IN WATER AS A FUNCTION OF TEMPERATURE AND SALINITY
function SED_MOD_1_ALUKAS_MOLDI_C(SVARNO, T, SAL, TVAR) result(MOL_DIFF)
    
    implicit none   
    !SVARNO : Stave variable no
    !T      : Temperature in Celcisus
    !SAL    : Salinity (ppt)
    !TVAR   : Generic temporal variable
    
    integer SVARNO
    double precision T
    double precision SAL
    double precision TVAR
    
    double precision TS
    double precision SS
    double precision V25
    double precision ONE
    double precision VTK
    double precision ZERO
    double precision D
    double precision TK
    
    double precision MOL_DIFF
    
    MOL_DIFF = 0.

    TK = T + 273.16
    
    !NH4N
    !Boudreau 1997, Springer-Verlag
    if (SVARNO.eq.1) then
        MOL_DIFF = (9.5D0 + 4.13D-1 * T) * 1.0D-6
    end if
    
    !NO3N
    !Boudreau 1997, Springer-Verlag
    if (SVARNO.eq.2) then
        MOL_DIFF = (9.5D0 + 3.88D-1 * T) * 1.0D-6
    end if
    
    !Dissolved organic nitrogen
    if (SVARNO.eq.3) then
        MOL_DIFF = 1.0D-6
    end if
    
    !Particulate organic nitrogen
    if (SVARNO.eq.4) then
        MOL_DIFF = 0.0D0
    end if
    
    !PO4P
    !Boudreau 1997, Springer-Verlag
    if (SVARNO.eq.5) then
        MOL_DIFF = (2.62D0 + 1.43D-1 * T) * 1.0D-6
    end if
    
    !Dissolved organic phosphorus
    if (SVARNO.eq.6) then
        MOL_DIFF = 1.0D-6
    end if
    
    !Particulate organic phosphorus
    if (SVARNO.eq.7) then
        MOL_DIFF = 0.0D0
    end if
    
    !DOXY
    !Reference : Fossing et al., 2004 (NERI Technical report)
    if (SVARNO.eq.8) then
        MOL_DIFF = (1.17D1 + (3.44D-1 * T) + &
                  (5.05D-3 * (T ** 2.0D0))) * 1.0D-6
    end if
    
    !Dissolved organic carbon
    if (SVARNO.eq.9) then
         MOL_DIFF = 1.0D-6
    end if
    
    !Particulate organic carbon
    if (SVARNO.eq.10) then
        MOL_DIFF = 0.0D0
    end if
    
    !DSi
    !From Boudreau 1997, Springer-Verlag
    !Wollast and Garrels (1971) found D(H4SiO4) at 25 deg C
    !and 36.1 ppt S., Assume that this value can be scaled by
    !the Stokes-Einstein relationship to any other temperature.
    if (SVARNO.eq.10) then
        MOL_DIFF = 1.0E-05
        TS = 25.0
        SS = 36.1
        call SED_MOD_1_CVISC(V25, SS , TS, 1.0D0)
        call SED_MOD_1_CVISC(VTK, SAL, T , 1.0D0)
        MOL_DIFF = MOL_DIFF * V25 / 298.16 * TK / VTK
    end if
    
    !Particulate silicon
    if (SVARNO.eq.12) then
        MOL_DIFF = 0.0D0
    end if
    
    !Dissolved inorganic carbon (to be corrected with better formulation)
    if (SVARNO.eq.13) then
        MOL_DIFF = 1.0D-6
    end if

    !Alkalinity (to be corrected with better formulation)
    if (SVARNO.eq.14) then
        MOL_DIFF = 1.0D-6
    end if

    !Salinity (to be corrected with better formulation)
    if (SVARNO.eq.15) then
        MOL_DIFF = 1.0D-6
    end if

    MOL_DIFF = MOL_DIFF * 1.0D-4
    
end function SED_MOD_1_ALUKAS_MOLDI_C
    
!***********************************************************************
!***********************************************************************


subroutine SED_MOD_1_CVISC(V,S,T,P)

    !Calculates the shear viscosity of water using the equation
    !given by Kukulka et al. (1987).
    !Calculated viscosity is in centipoise.
    !
    !Valid for 0<T<30 and 0<S<36.
    
    double precision V
    double precision S
    double precision T
    double precision P
    
    V = 1.7910 - T * (6.144D-02 - T*(1.4510D-03 - T*1.6826D-05)) - &
        1.5290D-04 * P + 8.3885D-08 * P * P + 2.4727D-03 * S + &
        (6.0574D-06*P - 2.6760D-09*P*P)*T + (T * (4.8429D-05 - &
        T * (4.7172D-06 - T * 7.5986D-08))) * S
    
end subroutine SED_MOD_1_CVISC

!*******************************************************************
!*******************************************************************
    
    
subroutine FLX_SED_MOD_1_TO_ALUKAS_II &
           (FLUXES_FROM_SEDIMENT, NUM_FLUXES_FROM_SEDIMENT, &
            FLUXES_TO_ALUKAS    , NUM_FLUXES_TO_ALUKAS)
   
    integer :: NUM_FLUXES_FROM_SEDIMENT
    integer :: NUM_FLUXES_TO_ALUKAS
    
    double precision :: FLUXES_FROM_SEDIMENT(NUM_FLUXES_FROM_SEDIMENT)
    double precision :: FLUXES_TO_ALUKAS    (NUM_FLUXES_TO_ALUKAS)
       
    !NH4 FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(1) = FLUXES_FROM_SEDIMENT(1)
    
    !NO3 FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(2) = FLUXES_FROM_SEDIMENT(2)
    
    !PO4 FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(3) = FLUXES_FROM_SEDIMENT(5)
    
    !DISS_OXYGEN FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(4) = FLUXES_FROM_SEDIMENT(8)
    
    !CHEMO_AUT_BAC_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(5) = 0.0D0
    
    !AER_HET_BAC_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(6) = 0.0D0
    
    !FAC_AN_AER_HET_BAC_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(7) = 0.0D0
    
    !DIA_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(8) = 0.0D0
    
    !ZOO_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(9) = 0.0D0
    
    !ZOO_N FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(10) = 0.0D0
    
    !ZOO_P FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(11) = 0.0D0
    
    !DET_PART_ORG_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(12) = 0.0D0
    
    !DET_PART_ORG_N FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(13) = 0.0D0
    
    !DET_PART_ORG_P FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(14) = 0.0D0
    
    !DISS_ORG_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(15) = FLUXES_FROM_SEDIMENT(9)
    
    !DISS_ORG_N FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(16) = FLUXES_FROM_SEDIMENT(3)
    
    !DISS_ORG_P FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(17) = FLUXES_FROM_SEDIMENT(6)      
    
    !CYN_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(18) = 0.0D0
    
    !OPA_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(19) = 0.0D0
                
    !DISS_Si FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(20) = FLUXES_FROM_SEDIMENT(11)
    
    !PART_Si FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(21) = 0.0D0 
 
    !FIX_CYN_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(22) = 0.0D0

    !INORG_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(23) = FLUXES_FROM_SEDIMENT(13)

    !INORG_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(24) = FLUXES_FROM_SEDIMENT(14)

end subroutine FLX_SED_MOD_1_TO_ALUKAS_II
 
!************************************************************
!************************************************************
