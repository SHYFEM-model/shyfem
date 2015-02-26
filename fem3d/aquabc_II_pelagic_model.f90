! Pelagic kinetic model ALUKAS_II
! Version with variables calculated in subroutines 
! Version with dissolved inorganic carbon and alkalinity as
! state variables.

!Contains:
! subroutine derived_vars
! SUBROUTINE cur_smith
! subroutine LIM_LIGHT
! function GROWTH_AT_TEMP
! function AMMONIA_PREF
! function DO_SATURATION
! function KAWIND
! function AMMONIA_VOLATILIZATION
! function STRANGERSD



!******************************************************************
!******************************************************************
!******************************************************************
subroutine PELAGIC_KINETICS &
           (node_active,       nkn, &
            STATE_VARIABLES  , DERIVATIVES  , nstate, &
            MODEL_CONSTANTS  , nconst               , &
            DRIVING_FUNCTIONS, n_driving_functions  , &
            FLAGS            , nflags               , &
            PROCESS_RATES    , NDIAGVAR             , &
            SAVED_OUTPUTS    , n_saved_outputs      , &
            WC_OUTPUTS       , noutput              , &
            TIME, TIME_STEP  , CALLED_BEFORE)

    use CO2SYS_CDIAC
    use AQUABC_II_GLOBAL
    use para_aqua

    implicit none
    include 'param.h'

    integer ipv(nkndim)	!external node numbers
    common  /ipv/ipv

    logical VALUE_strange(nkn) ! For NaN and Inf checking

    integer  :: nkn, nstate, noutput, nconst, n_driving_functions, nflags
    integer  :: n_saved_outputs, NDIAGVAR, STRANGERSD, error

    integer node_active(nkn) ! internal numbers of nodes, used for diagnostics. Not implemented yet

    double precision, dimension(nkn,nstate),             intent(in)   :: STATE_VARIABLES
    double precision, dimension(nkn,nstate),             intent(out)  :: DERIVATIVES
    double precision, dimension(nconst),                 intent(in)   :: MODEL_CONSTANTS
    double precision, dimension(nkn,n_driving_functions),intent(in)   :: DRIVING_FUNCTIONS
    integer,          dimension(nflags),                 intent(in)   :: FLAGS
    double precision, dimension(nkn,nstate, NDIAGVAR),   intent(out)  :: PROCESS_RATES
    double precision, dimension(nkn,noutput),            intent(inout):: WC_OUTPUTS   !For saving derived vars as nstate+1, ... elements
    double precision, dimension(n_saved_outputs),        intent(inout):: SAVED_OUTPUTS !For saving some variables to be used for all nodes?
    double precision, intent(in) :: TIME
    double precision, intent(in) :: TIME_STEP

!    integer, intent(in)          :: BOX_NO
    integer, intent(inout)       :: CALLED_BEFORE

    integer :: i,j,k

    !*********************************************'
    !*                                           *'
    !* INSERT YOUR PELAGIC ECOLOGY KINETICS HERE *'
    !*                                           *'
    !*********************************************'

    !FUNCTIONS
    double precision   :: DO_SATURATION
    double precision   :: KAWIND

    !VARIABLE DEFINITIONS

    !State variables
    double precision :: NH4_N              (nkn)
    double precision :: NO3_N              (nkn)
    double precision :: PO4_P              (nkn)
    double precision :: DISS_OXYGEN        (nkn)
    double precision :: CHEM_AUT_BAC_C     (nkn)
    double precision :: AER_HET_BAC_C      (nkn)
    double precision :: FAC_AN_HET_BAC_C   (nkn)
    double precision :: DIA_C              (nkn)
    double precision :: ZOO_C              (nkn)
    double precision :: ZOO_N              (nkn)
    double precision :: ZOO_P              (nkn)
    double precision :: DET_PART_ORG_C     (nkn)
    double precision :: DET_PART_ORG_N     (nkn)
    double precision :: DET_PART_ORG_P     (nkn)
    double precision :: DISS_ORG_C         (nkn)
    double precision :: DISS_ORG_N         (nkn)
    double precision :: DISS_ORG_P         (nkn)
    double precision :: CYN_C              (nkn)
    double precision :: OPA_C              (nkn)
    double precision :: DISS_Si            (nkn)
    double precision :: PART_Si            (nkn)
    double precision :: FIX_CYN_C          (nkn)
    ! New state variables added 22 September 2014
    real(kind=DBL_PREC) :: INORG_C(nkn)   !Inorganic carbon
    real(kind=DBL_PREC) :: TOT_ALK(nkn)   !Total alkalinity
    ! End of new state variables added 22 September 2014


    double precision :: PHYT_TOT_C(nkn)

    !Driving functions
    double precision :: TEMP     (nkn)
    double precision :: pH       (nkn)
    double precision :: FDAY     (nkn) 
    double precision :: DEPTH    (nkn)
    double precision :: I_A      (nkn)
    double precision :: SALT     (nkn)
    double precision :: ELEVATION(nkn)
    double precision :: AIRTEMP  (nkn)
    double precision :: WINDS    (nkn)
    double precision :: K_B_E    (nkn)
    double precision :: ice_cover(nkn)

    !Flags
    integer :: SAFE_MODE
    integer :: SURFACE_BOX

    !Model constants
    double precision :: K_A             !prescribed aeration coefficient, if negative should be calculated
    double precision :: K_A_CALC(nkn)   !calculated aeration reactor specific coefficient
    double precision :: THETA_K_A

    double precision :: KG_CHEM_AUT_BAC_20
    double precision :: THETA_KG_CHEM_AUT_BAC
    double precision :: KR_CHEM_AUT_BAC_20
    double precision :: THETA_KR_CHEM_AUT_BAC
    double precision :: KD_CHEM_AUT_BAC_20
    double precision :: THETA_KD_CHEM_AUT_BAC
    double precision :: KHS_NH4N_CHEM_AUT_BAC
    double precision :: KHS_PO4P_CHEM_AUT_BAC
    double precision :: KHS_O2_CHEM_AUT_BAC
    double precision :: DO_STR_HYPOX_CHEM_AUT_BAC_D
    double precision :: THETA_HYPOX_CHEM_AUT_BAC_D
    double precision :: EXPON_HYPOX_CHEM_AUT_BAC_D
    double precision :: CHEM_AUT_BAC_N_TO_C
    double precision :: CHEM_AUT_BAC_P_TO_C
    double precision :: CHEM_AUT_BAC_O2_TO_C
    double precision :: YIELD_CHEM_AUT_BAC

    double precision :: KG_AER_HET_BAC_20
    double precision :: THETA_KG_AER_HET_BAC
    double precision :: KR_AER_HET_BAC_20
    double precision :: THETA_KR_AER_HET_BAC
    double precision :: KD_AER_HET_BAC_20
    double precision :: THETA_KD_AER_HET_BAC
    double precision :: KHS_ORGC_AER_HET_BAC
    double precision :: KHS_ORGN_AER_HET_BAC
    double precision :: KHS_ORGP_AER_HET_BAC
    double precision :: OX_ORGN_AER_HET_BAC
    double precision :: OX_ORGP_AER_HET_BAC
    double precision :: KHS_O2_AER_HET_BAC
    double precision :: KHS_DIN_AER_HET_BAC
    double precision :: KHS_DIP_AER_HET_BAC
    double precision :: YIELD_OC_AER_HET_BAC
    double precision :: KHS_PHYT_AER_HET_BAC
    double precision :: DO_STR_HYPOX_AER_HET_BAC_D
    double precision :: THETA_HYPOX_AER_HET_BAC_D
    double precision :: EXPON_HYPOX_AER_HET_BAC_D
    double precision :: AER_HET_BAC_N_TO_C
    double precision :: AER_HET_BAC_P_TO_C
    double precision :: AER_HET_BAC_O2_TO_C

    double precision :: KG_FAC_AN_HET_BAC_20
    double precision :: THETA_KG_FAC_AN_HET_BAC
    double precision :: KR_FAC_AN_HET_BAC_20
    double precision :: THETA_KR_FAC_AN_HET_BAC
    double precision :: KD_FAC_AN_HET_BAC_20
    double precision :: THETA_KD_FAC_AN_HET_BAC
    double precision :: KHS_NO3N_FAC_AN_HET_BAC
    double precision :: KHS_ORGC_FAC_AN_HET_BAC
    double precision :: KHS_ORGN_FAC_AN_HET_BAC
    double precision :: KHS_ORGP_FAC_AN_HET_BAC
    double precision :: REV_KHS_O2_FAC_AN_HET_BAC
    double precision :: NO3N_LACK_STR_FAC_AN_HET_BAC_D
    double precision :: THETA_NO3_LACK_FAC_AN_HET_BAC_D
    double precision :: EXP_NO3_LACK_FAC_AN_HET_BAC_D
    double precision :: FAC_AN_HET_BAC_N_TO_C
    double precision :: FAC_AN_HET_BAC_P_TO_C
    double precision :: FAC_AN_HET_BAC_O2_TO_C
    double precision :: YIELD_FAC_AN_HET_BAC

    double precision :: KG_DIA_OPT_TEMP
    double precision :: KAPPA_DIA_UNDER_OPT_TEMP
    double precision :: KAPPA_DIA_OVER_OPT_TEMP
    double precision :: KR_DIA_20
    double precision :: THETA_KR_DIA
    double precision :: KD_DIA_20
    double precision :: THETA_KD_DIA
    double precision :: KHS_DIN_DIA
    double precision :: KHS_DIP_DIA
    double precision :: KHS_O2_DIA
    double precision :: KHS_NH4N_PREF_DIA
    double precision :: I_S_DIA
    double precision :: DO_STR_HYPOX_DIA_D
    double precision :: THETA_HYPOX_DIA_D
    double precision :: EXPON_HYPOX_DIA_D
    double precision :: DIA_N_TO_C
    double precision :: DIA_P_TO_C
    double precision :: DIA_O2_TO_C
    double precision :: DIA_C_TO_CHLA
    double precision :: DIA_C_TO_CHLA_NEW

    double precision :: KG_CYN_OPT_TEMP
    double precision :: KAPPA_CYN_UNDER_OPT_TEMP
    double precision :: KAPPA_CYN_OVER_OPT_TEMP
    double precision :: KR_CYN_20
    double precision :: THETA_KR_CYN
    double precision :: KD_CYN_20
    double precision :: THETA_KD_CYN
    double precision :: KHS_DIN_CYN
    double precision :: KHS_DIP_CYN
    double precision :: KHS_O2_CYN
    double precision :: KHS_NH4N_PREF_CYN
    double precision :: I_S_CYN
    double precision :: DO_STR_HYPOX_CYN_D
    double precision :: THETA_HYPOX_CYN_D
    double precision :: EXPON_HYPOX_CYN_D
    double precision :: CYN_N_TO_C
    double precision :: CYN_P_TO_C
    double precision :: CYN_O2_TO_C
    double precision :: CYN_C_TO_CHLA
    double precision :: CYN_C_TO_CHLA_NEW

    double precision :: KG_OPA_OPT_TEMP
    double precision :: KAPPA_OPA_UNDER_OPT_TEMP
    double precision :: KAPPA_OPA_OVER_OPT_TEMP
    double precision :: KR_OPA_20
    double precision :: THETA_KR_OPA
    double precision :: KD_OPA_20
    double precision :: THETA_KD_OPA
    double precision :: KHS_DIN_OPA
    double precision :: KHS_DIP_OPA
    double precision :: KHS_O2_OPA
    double precision :: KHS_NH4N_PREF_OPA
    double precision :: I_S_OPA
    double precision :: DO_STR_HYPOX_OPA_D
    double precision :: THETA_HYPOX_OPA_D
    double precision :: EXPON_HYPOX_OPA_D
    double precision :: OPA_N_TO_C
    double precision :: OPA_P_TO_C
    double precision :: OPA_O2_TO_C
    double precision :: OPA_C_TO_CHLA
    double precision :: OPA_C_TO_CHLA_NEW

    double precision :: KG_ZOO_OPT_TEMP
    double precision :: KAPPA_ZOO_UNDER_OPT_TEMP
    double precision :: KAPPA_ZOO_OVER_OPT_TEMP
    double precision :: GRAT_ZOO_DIA
    double precision :: GRAT_ZOO_CYN
    double precision :: GRAT_ZOO_OPA
    double precision :: GRAT_ZOO_FIX_CYN
    double precision :: GRAT_ZOO_CHEM_AUT_BAC
    double precision :: GRAT_ZOO_AER_HET_BAC
    double precision :: GRAT_ZOO_FAC_AN_HET_BAC
    double precision :: GRAT_ZOO_DET_PART_ORG_C
    double precision :: PREF_ZOO_DIA
    double precision :: PREF_ZOO_CYN
    double precision :: PREF_ZOO_FIX_CYN
    double precision :: PREF_ZOO_OPA
    double precision :: PREF_ZOO_CHEM_AUT_BAC
    double precision :: PREF_ZOO_AER_HET_BAC
    double precision :: PREF_ZOO_FAC_AN_HET_BAC
    double precision :: PREF_ZOO_DET_PART_ORG_C
    double precision :: KHS_DIA_C_ZOO
    double precision :: KHS_CYN_C_ZOO
    double precision :: KHS_OPA_C_ZOO
    double precision :: KHS_FIX_CYN_C_ZOO
    double precision :: KHS_CHEM_AUT_BAC_C_ZOO
    double precision :: KHS_AER_HET_BAC_C_ZOO
    double precision :: KHS_FAC_AN_HET_BAC_C_ZOO
    double precision :: KHS_DET_PART_ORG_C_ZOO
    
    double precision :: KHS_MIN_P          
    double precision :: KHS_MIN_N          
    double precision :: KHS_NITR_NH4_N    
    double precision :: KHS_NITR_OXY
    
    double precision :: KHS_AMIN_P 
    double precision :: KHS_AMIN_N 
    double precision :: KHS_DISS_N 
    double precision :: KHS_DISS_P 
       
         
    double precision :: FOOD_MIN_ZOO
    double precision :: THETA_KE_ZOO
    double precision :: KR_ZOO_20
    double precision :: THETA_KR_ZOO
    double precision :: KD_ZOO_20
    double precision :: THETA_KD_ZOO
    double precision :: DO_STR_HYPOX_ZOO_D
    double precision :: THETA_HYPOX_ZOO_D
    double precision :: EXPON_HYPOX_ZOO_D
    double precision :: ZOO_N_TO_C
    double precision :: ZOO_P_TO_C
    double precision :: ZOO_O2_TO_C
    double precision :: KDISS_DET_PART_ORG_C_20
    double precision :: THETA_KDISS_DET_PART_ORG_C
    double precision :: KDISS_DET_PART_ORG_N_20
    double precision :: THETA_KDISS_DET_PART_ORG_N
    double precision :: KDISS_DET_PART_ORG_P_20
    double precision :: THETA_KDISS_DET_PART_ORG_P

    double precision :: KHS_DSI_DIA
    double precision :: DIA_SI_TO_C
    double precision :: KDISS_PART_SI_20
    double precision :: THETA_KDISS_PART_SI

    double precision :: DIA_OPT_TEMP_LR
    double precision :: DIA_OPT_TEMP_UR
    double precision :: CYN_OPT_TEMP_LR
    double precision :: CYN_OPT_TEMP_UR
    double precision :: FIX_CYN_OPT_TEMP_LR
    double precision :: FIX_CYN_OPT_TEMP_UR
    double precision :: OPA_OPT_TEMP_LR
    double precision :: OPA_OPT_TEMP_UR
    double precision :: ZOO_OPT_TEMP_LR
    double precision :: ZOO_OPT_TEMP_UR
    
    double precision :: KE_ZOO                       
    double precision :: FRAC_ZOO_EX_ORG

    !Abitotic mineralization and nitrification (bacteria not modelled)
    double precision :: K_MIN_DOC_20
    double precision :: THETA_K_MIN_DOC
    double precision :: K_MIN_DON_20
    double precision :: THETA_K_MIN_DON
    double precision :: K_MIN_DOP_20
    double precision :: THETA_K_MIN_DOP
    double precision :: K_NITR_20
    double precision :: THETA_K_NITR

    double precision :: EFF_CHEM_AUT_BAC_GROWTH
    double precision :: EFF_AER_HET_BAC_GROWTH
    double precision :: EFF_FAC_AN_HET_BAC_GROWTH
    double precision :: EFF_DIA_GROWTH
    double precision :: EFF_CYN_GROWTH
    double precision :: EFF_OPA_GROWTH
    double precision :: EFF_ZOO_GROWTH

    !Nitrogen fixing cyanobacteria
    double precision :: KG_FIX_CYN_OPT_TEMP
    double precision :: EFF_FIX_CYN_GROWTH
    double precision :: KAPPA_FIX_CYN_UNDER_OPT_TEMP
    double precision :: KAPPA_FIX_CYN_OVER_OPT_TEMP
    double precision :: KR_FIX_CYN_20
    double precision :: THETA_KR_FIX_CYN
    double precision :: KD_FIX_CYN_20
    double precision :: THETA_KD_FIX_CYN
    double precision :: KHS_DIN_FIX_CYN
    double precision :: KHS_DIP_FIX_CYN
    double precision :: KHS_O2_FIX_CYN
    double precision :: KHS_NH4N_PREF_FIX_CYN
    double precision :: I_S_FIX_CYN
    double precision :: DO_STR_HYPOX_FIX_CYN_D
    double precision :: THETA_HYPOX_FIX_CYN_D
    double precision :: EXPON_HYPOX_FIX_CYN_D
    double precision :: FIX_CYN_N_TO_C
    double precision :: FIX_CYN_P_TO_C
    double precision :: FIX_CYN_O2_TO_C
    double precision :: FIX_CYN_C_TO_CHLA
    double precision :: FIX_CYN_C_TO_CHLA_NEW
    double precision :: R_FIX
    double precision :: K_FIX
    double precision :: FRAC_FIX_CYN_EXCR
    
    double precision :: FRAC_CYN_EXCR                
    double precision :: FRAC_OPA_EXCR                
    double precision :: FRAC_DIA_EXCR                
    double precision :: FAC_PHYT_AMIN_DON            
    double precision :: FAC_PHYT_AMIN_DOC             
    double precision :: FAC_PHYT_AMIN_DOP            
    double precision :: FAC_PHYT_DET_PART_ORG_C         
    double precision :: FAC_PHYT_DET_PART_ORG_P      
    double precision :: FAC_PHYT_DET_PART_ORG_N
    !end of constatnts     

    !Main process rates
    double precision :: R_AERATION               (nkn)
    double precision :: R_CHEM_AUT_BAC_GROWTH    (nkn)
    double precision :: R_CHEM_AUT_BAC_RESP      (nkn)
    double precision :: R_CHEM_AUT_BAC_INT_RESP  (nkn)
    double precision :: R_CHEM_AUT_BAC_TOT_RESP  (nkn)
    double precision :: R_CHEM_AUT_BAC_DEATH     (nkn)
    double precision :: R_AER_HET_BAC_GROWTH     (nkn)
    double precision :: R_AER_HET_BAC_RESP       (nkn)
    double precision :: R_AER_HET_BAC_INT_RESP   (nkn)
    double precision :: R_AER_HET_BAC_TOT_RESP   (nkn)
    double precision :: R_AER_HET_BAC_DEATH      (nkn)
    double precision :: R_FAC_AN_HET_BAC_GROWTH  (nkn)
    double precision :: R_FAC_AN_HET_BAC_RESP    (nkn)
    double precision :: R_FAC_AN_HET_BAC_INT_RESP(nkn)
    double precision :: R_FAC_AN_HET_BAC_TOT_RESP(nkn)
    double precision :: R_FAC_AN_HET_BAC_DEATH   (nkn)

    double precision :: R_DIA_GROWTH  (nkn)
    double precision :: R_DIA_MET     (nkn)
    double precision :: R_DIA_RESP    (nkn)
    double precision :: R_DIA_EXCR    (nkn)

    double precision :: R_DIA_INT_RESP(nkn)
    double precision :: R_DIA_TOT_RESP(nkn)
    double precision :: R_DIA_DEATH   (nkn)
    
    
       
    double precision :: R_CYN_GROWTH   (nkn)
    double precision :: R_CYN_RESP     (nkn)
    double precision :: R_CYN_MET      (nkn)
    double precision :: R_CYN_EXCR     (nkn)
    
    double precision :: R_CYN_INT_RESP (nkn)
    double precision :: R_CYN_TOT_RESP (nkn)
    double precision :: R_CYN_DEATH    (nkn)

    double precision :: R_FIX_CYN_GROWTH  (nkn)
    double precision :: R_FIX_CYN_RESP    (nkn)
    double precision :: R_FIX_CYN_MET     (nkn)
    double precision :: R_FIX_CYN_EXCR    (nkn)    
    double precision :: R_FIX_CYN_INT_RESP(nkn)
    double precision :: R_FIX_CYN_TOT_RESP(nkn)
    double precision :: R_FIX_CYN_DEATH   (nkn)


    double precision :: R_OPA_GROWTH   (nkn)
    double precision :: R_OPA_RESP     (nkn)
    double precision :: R_OPA_MET      (nkn)
    double precision :: R_OPA_EXCR     (nkn)

    double precision :: R_OPA_INT_RESP (nkn)
    double precision :: R_OPA_TOT_RESP (nkn)
    double precision :: R_OPA_DEATH    (nkn)           


    double precision :: R_ZOO_GROWTH                 (nkn)
    double precision :: R_ZOO_FEEDING_DIA            (nkn)
    double precision :: R_ZOO_FEEDING_CYN            (nkn)
    double precision :: R_ZOO_FEEDING_OPA            (nkn)
    double precision :: R_ZOO_FEEDING_FIX_CYN        (nkn)
    double precision :: R_ZOO_FEEDING_CHEM_AUT_BAC   (nkn)
    double precision :: R_ZOO_FEEDING_AER_HET_BAC    (nkn)
    double precision :: R_ZOO_FEEDING_FAC_AN_HET_BAC (nkn)
    double precision :: R_ZOO_FEEDING_DET_PART_ORG_C (nkn)
    double precision :: R_ZOO_RESP                   (nkn)
    double precision :: R_ZOO_INT_RESP               (nkn)
    double precision :: R_ZOO_TOT_RESP               (nkn)
    double precision :: R_ZOO_DEATH                  (nkn)
    double precision :: R_ZOO_EX_DON                 (nkn)
    double precision :: R_ZOO_EX_DOP                 (nkn)
    double precision :: R_ZOO_EX_DOC                 (nkn)
    double precision :: R_ZOO_EX_C                   (nkn)
    
    
    double precision :: R_ZOO_EX_NH4                 (nkn)
    double precision :: R_ZOO_EX_PO4                 (nkn)
    double precision :: R_DET_PART_ORG_C_DISSOLUTION (nkn)
    double precision :: LIM_PHYT_DISS_DET_PART_ORG_C (nkn)


    double precision :: R_AER_HET_BAC_N_OX (nkn)
    double precision :: LIM_N_MIN_DON_N    (nkn)    
    double precision :: R_AER_HET_BAC_P_OX (nkn)
    double precision :: LIM_P_MIN_DOP_P    (nkn) 
        
    double precision :: R_FAC_AN_HET_BAC_N_OX (nkn)
    double precision :: R_FAC_AN_HET_BAC_P_OX (nkn)
    double precision :: R_PART_SI_DISS        (nkn)
    double precision :: R_NON_FIX_CYN_GROWTH  (nkn)
    double precision :: R_FIX_FIX_CYN_GROWTH  (nkn)

    double precision :: R_ABIOTIC_DOC_MIN (nkn)
    double precision :: R_ABIOTIC_DON_MIN (nkn)
    double precision :: R_ABIOTIC_DOP_MIN (nkn)                   
                       

    double precision :: LIM_PHYT_AMIN_DOC  (nkn)
    double precision :: LIM_N_AMIN_DON     (nkn)    
    double precision :: LIM_PHY_N_AMIN_DON (nkn)
    double precision :: LIM_P_AMIN_DOP     (nkn)    
    double precision :: LIM_PHY_P_AMIN_DOP (nkn)
    

    double precision :: R_ABIOTIC_NITR (nkn)
    double precision :: LIM_NITR_OXY   (nkn)    
    double precision :: LIM_NITR_NH4_N (nkn)


    double precision :: R_AMMONIA_VOLATIL (nkn)   !Ammonia volatilization (mgN/L)

    !Derived process rates
    double precision :: R_NITRIFICATION              (nkn)
    double precision :: R_DENITRIFICATION            (nkn)

    double precision :: R_DET_PART_ORG_N_DISSOLUTION (nkn)
    double precision :: LIM_N_DISS_DET_PART_ORG_N    (nkn)
    
    
    double precision :: LIM_PHY_N_DISS_DET_PART_ORG_N(nkn)
    double precision :: R_DET_PART_ORG_P_DISSOLUTION (nkn)
    double precision :: LIM_P_DISS_DET_PART_ORG_P    (nkn)
    
    
    double precision :: LIM_PHY_P_DISS_DET_PART_ORG_P(nkn)
    
    
    double precision :: FAC_NO3N_LACK_FAC_AN_HET_BAC_D (nkn)



    !Auxillary variables
    double precision :: DISS_OXYGEN_SAT                (nkn)
    double precision :: FAC_HYPOX_CHEM_AUT_BAC_D       (nkn)
    double precision :: KD_CHEM_AUT_BAC                (nkn)
    double precision :: FAC_HYPOX_AER_HET_BAC_D        (nkn)
    double precision :: KD_AER_HET_BAC                 (nkn)
    
    double precision :: KD_FAC_AN_HET_BAC              (nkn)
    double precision :: CHLA                           (nkn)
    double precision :: K_E                            (nkn)
    double precision :: KG_DIA                         (nkn)
    double precision :: LIM_KG_DIA                     (nkn)
    double precision :: FAC_HYPOX_DIA_D                (nkn)
    double precision :: PREF_NH4N_DIA                  (nkn)
    double precision :: KD_DIA                         (nkn)
    double precision :: KG_ZOO                         (nkn)
    double precision :: KG_ZOO_DIA                     (nkn)
    double precision :: KG_ZOO_CHEM_AUT_BAC            (nkn)
    double precision :: KG_ZOO_AER_HET_BAC             (nkn)
    double precision :: KG_ZOO_FAC_AN_HET_BAC          (nkn)
    double precision :: KG_ZOO_DET_PART_ORG_C          (nkn)
    double precision :: FOOD_AVAIL_ZOO                 (nkn)
    double precision :: FOOD_FACTOR_ZOO_DIA            (nkn)
    double precision :: FOOD_FACTOR_ZOO_CHEM_AUT_BAC   (nkn)
    double precision :: FOOD_FACTOR_ZOO_AER_HET_BAC    (nkn)
    double precision :: FOOD_FACTOR_ZOO_FAC_AN_HET_BAC (nkn)
    double precision :: FOOD_FACTOR_ZOO_DET_PART_ORG_C (nkn)
    double precision :: ACTUAL_ZOO_N_TO_C              (nkn)
    double precision :: ACTUAL_ZOO_P_TO_C              (nkn)

    ! not used
    integer          :: ZOO_N_EXCESS
    integer          :: ZOO_P_EXCESS
    integer          :: ZOO_N_DEFICIT
    integer          :: ZOO_P_DEFICIT
    double precision :: EXPECTED_ZOO_C
    double precision :: EXPECTED_ZOO_N
    double precision :: EXPECTED_ZOO_P
    double precision :: ACTUAL_N_TO_C_OVER_N_TO_C
    double precision :: ACTUAL_P_TO_C_OVER_P_TO_C
    ! end not used

    double precision :: FAC_HYPOX_ZOO_D  (nkn)
    double precision :: ACTUAL_DET_N_TO_C(nkn)
    double precision :: ACTUAL_DET_P_TO_C(nkn)
    double precision :: KD_ZOO           (nkn)

    double precision :: KD_FIX_CYN         (nkn)
    double precision :: FAC_HYPOX_FIX_CYN_D(nkn)
    double precision :: PREF_NH4N_FIX_CYN  (nkn)

    double precision :: KG_ZOO_FIX_CYN               (nkn)
    double precision :: LIM_DISS_NO3_N_FAC_AN_HET_BAC(nkn)
    double precision :: FOOD_FACTOR_ZOO_FIX_CYN      (nkn)

    double precision :: LIM_KG_DIA_DOXY (nkn)
    double precision :: LIM_KG_DIA_NUTR (nkn)
    double precision :: LIM_KG_DIA_LIGHT(nkn)

    double precision :: DIA_LIGHT_SAT    (nkn)
    double precision :: CYN_LIGHT_SAT    (nkn)
    double precision :: FIX_CYN_LIGHT_SAT(nkn)
    double precision :: OPA_LIGHT_SAT    (nkn)

    
    
    double precision :: LIM_KG_CYN_DOXY    (nkn)
    double precision :: LIM_KG_CYN_NUTR    (nkn)
    double precision :: LIM_KG_CYN_LIGHT   (nkn)
    double precision :: LIM_KG_CYN         (nkn)
    double precision :: KD_CYN             (nkn)
    double precision :: FAC_HYPOX_CYN_D    (nkn)
    double precision :: PREF_NH4N_CYN      (nkn)
    double precision :: KG_OPA             (nkn)
    double precision :: KG_CYN             (nkn)
    double precision :: LIM_KG_OPA_DOXY    (nkn)
    double precision :: LIM_KG_OPA_NUTR    (nkn)
    double precision :: LIM_KG_OPA_LIGHT   (nkn)
    double precision :: LIM_KG_OPA         (nkn)
    double precision :: KD_OPA             (nkn)
    double precision :: FAC_HYPOX_OPA_D    (nkn)
    double precision :: PREF_NH4N_OPA      (nkn)
    double precision :: KG_ZOO_CYN         (nkn)
    double precision :: KG_ZOO_OPA         (nkn)
    
    double precision :: FOOD_FACTOR_ZOO_CYN(nkn)
    double precision :: FOOD_FACTOR_ZOO_OPA(nkn)

    double precision :: PREF_NH4N_AER_HET_BAC(nkn)

    !limitation factors
    double precision :: LIM_TEMP_CHEM_AUT_BAC (nkn)
    double precision :: LIM_NH4_N_CHEM_AUT_BAC(nkn)
    double precision :: LIM_PO4_P_CHEM_AUT_BAC(nkn)
    double precision :: LIM_OXY_CHEM_AUT_BAC  (nkn)

    double precision :: LIM_TEMP_AER_HET_BAC      (nkn)
    double precision :: LIM_DISS_ORG_C_AER_HET_BAC(nkn)
    double precision :: LIM_DISS_ORG_N_AER_HET_BAC(nkn)
    double precision :: LIM_DISS_ORG_P_AER_HET_BAC(nkn)
    double precision :: LIM_OXY_AER_HET_BAC       (nkn)
    double precision :: LIM_DIN_AER_HET_BAC       (nkn)
    double precision :: LIM_DIP_AER_HET_BAC       (nkn)
    double precision :: LIM_PHYT_C_AER_HET_BAC    (nkn)

    double precision :: LIM_TEMP_FAC_AN_HET_BAC      (nkn)
    double precision :: LIM_DISS_ORG_C_FAC_AN_HET_BAC(nkn)
    double precision :: LIM_DISS_ORG_N_FAC_AN_HET_BAC(nkn)
    double precision :: LIM_DISS_ORG_P_FAC_AN_HET_BAC(nkn)
    double precision :: LIM_OXY_FAC_AN_HET_BAC       (nkn)

    double precision :: LIM_KG_DIA_TEMP   (nkn)
    double precision :: LIM_KG_DIA_N      (nkn)
    double precision :: LIM_KG_DIA_P      (nkn)
    double precision :: LIM_KG_DIA_DISS_Si(nkn)

    double precision :: LIM_KG_CYN_TEMP(nkn)
    double precision :: LIM_KG_CYN_N   (nkn)
    double precision :: LIM_KG_CYN_P   (nkn)

    double precision :: LIM_KG_OPA_TEMP(nkn)
    double precision :: LIM_KG_OPA_N   (nkn)
    double precision :: LIM_KG_OPA_P   (nkn)

    double precision KG_FIX_CYN         (nkn)
    
    double precision LIM_KG_FIX_CYN_TEMP    (nkn)
    double precision LIM_KG_FIX_CYN_DOXY    (nkn)
    double precision LIM_KG_FIX_CYN_LIGHT   (nkn)
    double precision LIM_KG_NON_FIX_CYN_N   (nkn)
    double precision LIM_KG_NON_FIX_CYN_P   (nkn)
    double precision LIM_KG_NON_FIX_CYN_NUTR(nkn)
    double precision LIM_KG_NON_FIX_CYN     (nkn)

    double precision LIM_KG_FIX_FIX_CYN_N   (nkn)
    double precision LIM_KG_FIX_FIX_CYN_P   (nkn)
    double precision LIM_KG_FIX_FIX_CYN_NUTR(nkn)
    double precision LIM_KG_FIX_FIX_CYN     (nkn)

    integer smith
    double precision FIX_CYN_DEPTH ! equivalent depth for the fixers selfshading (assumed they are on the surface)
! 1 - Smith light limitation with constant C:chla
! 0 - Ali light limitation


    !18 September 2014
  !New variables for DIC and ALK

    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_PAR1
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_PAR2
    integer,              allocatable, dimension (:) :: CO2SYS_PAR1TYPE
    integer,              allocatable, dimension (:) :: CO2SYS_PAR2TYPE
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_SALT
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_TEMPIN
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_TEMPOUT
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_PRESIN
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_PRESOUT
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_SI
    real(kind = DBL_PREC),allocatable, dimension (:) :: CO2SYS_PO4
    integer,              allocatable, dimension (:) :: CO2SYS_pHSCALEIN
    integer,              allocatable, dimension (:) :: CO2SYS_K1K2CONSTANTS
    integer,              allocatable, dimension (:) :: CO2SYS_KSO4CONSTANTS

    real(kind = DBL_PREC),allocatable, dimension (:,:) :: CO2SYS_OUT_DATA
    character(len=34),    allocatable, dimension (:) :: CO2SYS_NICEHEADERS

    integer :: CO2SYS_NUM_SAMPLES
    integer :: RUN_CO2SYS
    parameter (RUN_CO2SYS = 1)
    integer :: CO2SYS_ntps     ! number of nodes

    ! Variables and constansts for Dissolved Inorganic Carbon
    real(kind = DBL_PREC) :: TOTAL_DIC_KINETIC_SOURCES(nkn)
    real(kind = DBL_PREC) :: TOTAL_DIC_KINETIC_SINKS  (nkn)

    real(kind = DBL_PREC) :: T_A                   (nkn)
    real(kind = DBL_PREC) :: P_K_H                 (nkn)
    real(kind = DBL_PREC) :: K_H                   (nkn)
    real(kind = DBL_PREC) :: POWER                 (nkn)
    real(kind = DBL_PREC) :: P_CO2                 (nkn)
    real(kind = DBL_PREC) :: CO2_SAT               (nkn)
    real(kind = DBL_PREC) :: K_A_CALC_CO2          (nkn)
    real(kind = DBL_PREC) :: CO2_ATM_EXHANGE       (nkn)
    real(kind = DBL_PREC) :: DIC_KINETIC_DERIVATIVE(nkn)
    real(kind = DBL_PREC) :: ALPHA_0               (nkn)
    real(kind = DBL_PREC) :: ALPHA_1               (nkn)

    ! Variables and constansts for Alkalinity?
    real(kind = DBL_PREC) :: PKH                        (nkn)
    real(kind = DBL_PREC) :: FRAC_NH3                   (nkn)
    real(kind = DBL_PREC) :: FRAC_NH4                   (nkn)
    real(kind = DBL_PREC) :: N_CHEM_AUT_BAC_TOT_RESP    (nkn)
    real(kind = DBL_PREC) :: N_AER_HET_BAC_INT_RESP     (nkn)
    real(kind = DBL_PREC) :: N_FAC_AN_HET_BAC_TOT_RESP  (nkn)
    real(kind = DBL_PREC) :: N_DIA_TOT_RESP             (nkn)
    real(kind = DBL_PREC) :: N_CYN_TOT_RESP             (nkn)
    real(kind = DBL_PREC) :: N_OPA_TOT_RESP             (nkn)
    real(kind = DBL_PREC) :: N_FIX_CYN_TOT_RESP         (nkn)
    real(kind = DBL_PREC) :: N_AER_HET_BAC_N_OX         (nkn)
    real(kind = DBL_PREC) :: N_FAC_AN_HET_BAC_N_OX      (nkn)
    real(kind = DBL_PREC) :: N_ZOO_TOT_RESP             (nkn)
    real(kind = DBL_PREC) :: N_ABIOTIC_DON_MIN          (nkn)
    real(kind = DBL_PREC) :: ALK_GAINED_BY_AMMONIUM_GEN (nkn)
    real(kind = DBL_PREC) :: N_DENITRIFICATION          (nkn)
    real(kind = DBL_PREC) :: N_DIA_GROWTH               (nkn)
    real(kind = DBL_PREC) :: N_CYN_GROWTH               (nkn)
    real(kind = DBL_PREC) :: N_OPA_GROWTH               (nkn)
    real(kind = DBL_PREC) :: N_NON_FIX_CYN_GROWTH       (nkn)
    real(kind = DBL_PREC) :: N_AER_HET_BAC_GROWTH       (nkn)
    real(kind = DBL_PREC) :: ALK_GAINED_BY_NITRATE_CONS (nkn)
    real(kind = DBL_PREC) :: N_CHEM_AUT_BAC_GROWTH      (nkn)
    real(kind = DBL_PREC) :: ALK_LOST_BY_AMMONIUM_CONS  (nkn)
    real(kind = DBL_PREC) :: N_NITRIFICATION_NH4        (nkn)
    real(kind = DBL_PREC) :: N_NITRIFICATION_NH3        (nkn)
    real(kind = DBL_PREC) :: ALK_LOST_BY_NITRIFICATION  (nkn)
    real(kind = DBL_PREC) :: H_PLUS                     (nkn)

         ! Indicators
    integer :: TEST_MODE
    integer :: KP_OPTION
         ! Constansts:
    real(kind = DBL_PREC) :: A_1
    real(kind = DBL_PREC) :: A_2
    real(kind = DBL_PREC) :: A_3
    real(kind = DBL_PREC) :: B_1
    real(kind = DBL_PREC) :: B_2
    real(kind = DBL_PREC) :: C_1
    real(kind = DBL_PREC) :: C_2

    real(kind = DBL_PREC) :: K_ONE_TIP                   (nkn)
    real(kind = DBL_PREC) :: K_TWO_TIP                   (nkn)
    real(kind = DBL_PREC) :: K_THREE_TIP                 (nkn)
    real(kind = DBL_PREC) :: FRACTION_DIVISOR_TIP        (nkn)
    real(kind = DBL_PREC) :: ALPHA_H2PO4                 (nkn)
    real(kind = DBL_PREC) :: ALPHA_HPO4                  (nkn)
    real(kind = DBL_PREC) :: ALPHA_PO4                   (nkn)
    real(kind = DBL_PREC) :: PHOSPHATE_EQ_CONSTANT       (nkn)
    real(kind = DBL_PREC) :: ALK_GAINED_BY_PHOSPHATE_CONS(nkn)
    real(kind = DBL_PREC) :: ALK_LOST_BY_PHOSPHATE_GEN   (nkn)
    real(kind = DBL_PREC) :: ALK_KINETIC_DERIVATIVE      (nkn)

    real(kind = DBL_PREC) :: P_CHEM_AUT_BAC_TOT_RESP     (nkn)
    real(kind = DBL_PREC) :: P_AER_HET_BAC_INT_RESP      (nkn)
    real(kind = DBL_PREC) :: P_FAC_AN_HET_BAC_TOT_RESP   (nkn)
    real(kind = DBL_PREC) :: P_DIA_TOT_RESP              (nkn)
    real(kind = DBL_PREC) :: P_CYN_TOT_RESP              (nkn)
    real(kind = DBL_PREC) :: P_OPA_TOT_RESP              (nkn)
    real(kind = DBL_PREC) :: P_FIX_CYN_TOT_RESP          (nkn)
    real(kind = DBL_PREC) :: P_AER_HET_BAC_N_OX          (nkn)
    real(kind = DBL_PREC) :: P_FAC_AN_HET_BAC_N_OX       (nkn)
    real(kind = DBL_PREC) :: P_ZOO_TOT_RESP              (nkn)
    real(kind = DBL_PREC) :: P_ABIOTIC_DON_MIN           (nkn)
    real(kind = DBL_PREC) :: P_DENITRIFICATION           (nkn)
    real(kind = DBL_PREC) :: P_DIA_GROWTH                (nkn)
    real(kind = DBL_PREC) :: P_CYN_GROWTH                (nkn)
    real(kind = DBL_PREC) :: P_OPA_GROWTH                (nkn)
    real(kind = DBL_PREC) :: P_NON_FIX_CYN_GROWTH        (nkn)
    real(kind = DBL_PREC) :: P_AER_HET_BAC_GROWTH        (nkn)
    real(kind = DBL_PREC) :: P_AER_HET_BAC_P_OX          (nkn)
    real(kind = DBL_PREC) :: P_FAC_AN_HET_BAC_P_OX       (nkn)
    real(kind = DBL_PREC) :: P_CHEM_AUT_BAC_GROWTH       (nkn)


    integer :: CONSIDER_ALKALNITY_DERIVATIVE
    integer :: CONSIDER_INORG_C_DERIVATIVE
    integer :: CONSIDER_CO2_REARATION

    integer :: NUM_SCREEN_OUTPUT_NODES
    parameter (NUM_SCREEN_OUTPUT_NODES = 0)

    integer, dimension(NUM_SCREEN_OUTPUT_NODES) :: SCREEN_OUTPUT_NODES

    integer, dimension(nstate) :: DERIVS_FOR_SCREEN_OUTPUTS

    ! For strange values processing
    ! Strange values

    real(kind = DBL_PREC),allocatable, dimension (:) :: STRANGERS

    ! node numbers with strange values (numbers of active nodes in the layer)
    integer              ,allocatable, dimension (:) :: NODES_STRANGE

    !internal node numbers with strange values
    integer              ,allocatable, dimension (:) :: NODES_STRANGE_int

    !external node numbers with strange values
    integer              ,allocatable, dimension (:) :: NODES_STRANGE_ext

    ! number of nodes with strange values
    integer :: nstrange
    logical debug_stranger

    integer :: NODE_COUNT
    !End of new variables for DIC and ALK

    debug_stranger = .true. ! True if check for strange values
    debug_stranger = .false.

    CONSIDER_ALKALNITY_DERIVATIVE = 1
    CONSIDER_INORG_C_DERIVATIVE   = 1
    CONSIDER_CO2_REARATION        = 1

    ! indicator for light limitation  algorithm:
    ! 1 - modified Smith, 0 - Ali
    smith = 1

    PROCESS_RATES(:,:,:) = 0.0D0

    error = 0

    if(debug_stranger) then
        do i = 1, nstate
            if(STRANGERSD(STATE_VARIABLES(1:nkn,i),VALUE_strange,nkn) .eq. 1) then
                nstrange = count(VALUE_strange)
                allocate(STRANGERS    (nstrange))
                allocate(NODES_STRANGE(nstrange))
                allocate(NODES_STRANGE_int(nstrange))
                allocate(NODES_STRANGE_ext(nstrange))

                j=1
                do k=1,nkn
                  if(VALUE_strange(k)) then
                   STRANGERS    (j) = STATE_VARIABLES(k,i)
                   NODES_STRANGE(j) = k
                   NODES_STRANGE_int(j) = node_active(k)
                   NODES_STRANGE_ext(j) = ipv(node_active(k))
                   j=j+1
                  end if
                end do

                print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                print *, 'PELAGIC_KINETICS:  Variable ',i !'Cell ',k
                print *, 'Initial state is NaN or Inf:'
                print *, 'NODE_NUMBERS int.=',NODES_STRANGE_int
                print *, 'NODE_NUMBERS ext.=',NODES_STRANGE_ext
                print *, 'VALUES=',STRANGERS

                deallocate(STRANGERS)
                deallocate(NODES_STRANGE)
                deallocate(NODES_STRANGE_int)
                deallocate(NODES_STRANGE_ext)

                error =1
            end if
        end do

        if(error .eq. 1) stop
    end if

    !INITIALIZE STATE VARIABLES
    if(nstate.ne.24) then
        print *, 'PELAGIC_KINETICS: Number of state variables is wrong', nstate
        stop
    end if

    ! When changed this should be also updated in derived_vars
    NH4_N           (:)      = STATE_VARIABLES(:,1)      ! AMMONIUM NITROGEN
    NO3_N           (:)      = STATE_VARIABLES(:,2)      ! NITRATE NITROGEN
    PO4_P           (:)      = STATE_VARIABLES(:,3)      ! ORTHOPHOSPHATE PHOSPHORUS
    DISS_OXYGEN     (:)      = STATE_VARIABLES(:,4)      ! DISSOLVED OXYGEN
    CHEM_AUT_BAC_C  (:)      = STATE_VARIABLES(:,5)      ! CHEMOAUTIC BACTERIA CARBON
    AER_HET_BAC_C   (:)      = STATE_VARIABLES(:,6)      ! AEROBIC HETEROTROPHIC BACTERIA CARBON
    FAC_AN_HET_BAC_C(:)      = STATE_VARIABLES(:,7)      ! FACULTATIVE ANAEROBIC HETEROTROPHIC BACTERIA CARBON
    DIA_C           (:)      = STATE_VARIABLES(:,8)      ! DIATOMS CARBON
    ZOO_C           (:)      = STATE_VARIABLES(:,9)      ! ZOOPLANKTON CARBON
    ZOO_N           (:)      = STATE_VARIABLES(:,10)     ! ZOOPLANKTON NITROGEN
    ZOO_P           (:)      = STATE_VARIABLES(:,11)     ! ZOOPLANKTON PHOSPHORUS
    DET_PART_ORG_C  (:)      = STATE_VARIABLES(:,12)     ! DETRITUS PARTICULATE ORG. CARBON
    DET_PART_ORG_N  (:)      = STATE_VARIABLES(:,13)     ! DETRITUS PARTICULATE ORG. NITROGEN
    DET_PART_ORG_P  (:)      = STATE_VARIABLES(:,14)     ! DETRITUS PARTICULATE ORG. PHOSPHORUS
    DISS_ORG_C      (:)      = STATE_VARIABLES(:,15)     ! DISSOLVED ORGANIC CARBON
    DISS_ORG_N      (:)      = STATE_VARIABLES(:,16)     ! DISSOLVED ORGANIC NITROGEN
    DISS_ORG_P      (:)      = STATE_VARIABLES(:,17)     ! DISSOLVED ORGANIC PHOSPHORUS
    CYN_C           (:)      = STATE_VARIABLES(:,18)     ! NON FIXING CYANOBACTERIA CARBON
    OPA_C           (:)      = STATE_VARIABLES(:,19)     ! OTHER PHYTOPLANKTON CARBON
    DISS_Si         (:)      = STATE_VARIABLES(:,20)     ! DISSOLOVED SILICA
    PART_Si         (:)      = STATE_VARIABLES(:,21)     ! PARTICULATE SILICA
    FIX_CYN_C       (:)      = STATE_VARIABLES(:,22)     ! FIXING CYANOBACTERIA CARBON
    INORG_C         (:)      = STATE_VARIABLES(:,23)     ! INORG CARBON CARBON
    TOT_ALK         (:)      = STATE_VARIABLES(:,24)     ! TOTAL ALKALNITY

    !INITIALIZE DRIVING_FUNCTIONS
    if(n_driving_functions.ne.10) then
        print *, 'PELAGIC_KINETICS: Number of elements in DRIVING_FUNCTIONS is not 10', &
                 n_driving_functions
        stop
    end if

    TEMP     (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 1)
    SALT     (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 2)
    I_A      (1:nkn) =(DRIVING_FUNCTIONS(1:nkn, 3) * 5.0D-1 * 8.64D4 * 0.238846) / 1.0D4
    FDAY     (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 4)
    AIRTEMP  (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 5)
    WINDS    (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 6)
    ELEVATION(1:nkn) = DRIVING_FUNCTIONS(1:nkn, 7)
    DEPTH    (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 8)
    K_B_E    (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 9)
    ice_cover(1:nkn) = DRIVING_FUNCTIONS(1:nkn,10)

    !INITIALIZE FLAGS
    if(nflags.ne.2) then
       print *, 'PELAGIC_KINETICS: Number of elements in FLAGS is wrong', nflags
       stop
    end if

    SAFE_MODE   = FLAGS(1)
    SURFACE_BOX = FLAGS(2)


    !INITIALIZE MODEL CONSTANTS
    if(nconst.ne.229) then
        print *, 'PELAGIC_KINETICS: Number of elements in MODEL_CONSTANTS is wrong', nconst
        stop
    end if

! When changed this should be also updated in derived_vars.
! Note: this list does not correspond to existing list in parameters file.
!                               K_A   =   MODEL_CONSTANTS(  1) !  1! Aeration coefficient (if negative calculates internally)
!                         THETA_K_A   =   MODEL_CONSTANTS(  2) !  2! Temperature correction factor for aeration
!                KG_CHEM_AUT_BAC_20   =   MODEL_CONSTANTS(  3) !  3! Chemoautotrophic bacteria Growth rate
!           EFF_CHEM_AUT_BAC_GROWTH   =   MODEL_CONSTANTS(  4) !  4! Chemoautotrophic bacteria growth efficiency
!             THETA_KG_CHEM_AUT_BAC   =   MODEL_CONSTANTS(  5) !  5! Chemoautotrophic bacteria Temperature correction for growth rate
!                KR_CHEM_AUT_BAC_20   =   MODEL_CONSTANTS(  6) !  6! Chemoautotrophic bacteria Respiration rate
!             THETA_KR_CHEM_AUT_BAC   =   MODEL_CONSTANTS(  7) !  7! Chemoautotrophic bacteria Temperature correction for respiration rate
!                KD_CHEM_AUT_BAC_20   =   MODEL_CONSTANTS(  8) !  8! Chemoautotrophic bacteria Mortality rate
!             THETA_KD_CHEM_AUT_BAC   =   MODEL_CONSTANTS(  9) !  9! Chemoautotrophic bacteria Temperature correction for Mortality rate
!             KHS_NH4N_CHEM_AUT_BAC   =   MODEL_CONSTANTS( 10) ! 10! Chemoautotrophic bacteria Half saturation growth for NH4N
!             KHS_PO4P_CHEM_AUT_BAC   =   MODEL_CONSTANTS( 11) ! 11! Chemoautotrophic bacteria Half saturation growth for PO4P
!               KHS_O2_CHEM_AUT_BAC   =   MODEL_CONSTANTS( 12) ! 12! Chemoautotrophic bacteria Half saturation growth for O2
!       DO_STR_HYPOX_CHEM_AUT_BAC_D   =   MODEL_CONSTANTS( 13) ! 13! Chemoautotrophic bacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!        THETA_HYPOX_CHEM_AUT_BAC_D   =   MODEL_CONSTANTS( 14) ! 14! Chemoautotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
!        EXPON_HYPOX_CHEM_AUT_BAC_D   =   MODEL_CONSTANTS( 15) ! 15! Chemoautotrophic bacteria Exponent constant for Dissolved oxygen stress
!               CHEM_AUT_BAC_N_TO_C   =   MODEL_CONSTANTS( 16) ! 16! Chemoautotrophic bacteria Nitrogen to Carbon ratio
!               CHEM_AUT_BAC_P_TO_C   =   MODEL_CONSTANTS( 17) ! 17! Chemoautotrophic bacteria Phosphorus to Carbon ratio
!              CHEM_AUT_BAC_O2_TO_C   =   MODEL_CONSTANTS( 18) ! 18! Chemoautotrophic bacteria Oxygen to Carbon ratio
!                YIELD_CHEM_AUT_BAC   =   MODEL_CONSTANTS( 19) ! 19! Chemoautotrophic bacteria Yield of Carbon per unit amonia nitrogen
!                 KG_AER_HET_BAC_20   =   MODEL_CONSTANTS( 20) ! 20! Aerobic heterotrophic bacteria Growth rate
!            EFF_AER_HET_BAC_GROWTH   =   MODEL_CONSTANTS( 21) ! 21! Aerobic heterotrophic bacteria growth efficiency
!              THETA_KG_AER_HET_BAC   =   MODEL_CONSTANTS( 22) ! 22! Aerobic heterotrophic bacteria Temperature correction for growth rate
!                 KR_AER_HET_BAC_20   =   MODEL_CONSTANTS( 23) ! 23! Aerobic heterotrophic bacteria Respiration rate
!              THETA_KR_AER_HET_BAC   =   MODEL_CONSTANTS( 24) ! 24! Aerobic heterotrophic bacteria Temperature correction for respiration rate
!                 KD_AER_HET_BAC_20   =   MODEL_CONSTANTS( 25) ! 25! Aerobic heterotrophic bacteria Mortality rate
!              THETA_KD_AER_HET_BAC   =   MODEL_CONSTANTS( 26) ! 26! Aerobic heterotrophic bacteria Temperature correction for Mortality rate
!              KHS_ORGC_AER_HET_BAC   =   MODEL_CONSTANTS( 27) ! 27! Aerobic heterotrophic bacteria Half saturation growth for OC
! 
!              KHS_ORGN_AER_HET_BAC   =   MODEL_CONSTANTS( 28)  ! 28! Aerobic heterotrophic bacteria Half saturation growth for ON
!              KHS_ORGP_AER_HET_BAC   =   MODEL_CONSTANTS( 29)  ! 29! Aerobic heterotrophic bacteria Half saturation growth for OP
!                KHS_O2_AER_HET_BAC   =   MODEL_CONSTANTS( 30)  ! 30! Aerobic heterotrophic bacteria Half saturation growth for Oxygen
!               KHS_DIN_AER_HET_BAC   =   MODEL_CONSTANTS( 31)  ! 31! Aerobic heterotrophic bacteria Half saturation growth for inorganic nitrogen
!               KHS_DIP_AER_HET_BAC   =   MODEL_CONSTANTS( 32)  ! 32! Aerobic heterotrophic bacteria Half saturation growth for inorganic phosphorus
!              KHS_PHYT_AER_HET_BAC   =   MODEL_CONSTANTS( 33)  ! 33! Aerobic heterotrophic bacteria Half saturation growth for Phytoplankton C (not used as a resource)
!              YIELD_OC_AER_HET_BAC   =   MODEL_CONSTANTS( 34)  ! 34! Aerobic heterotrophic bacteria Yield of bacteria carbon per unit of organic carbon
!               OX_ORGN_AER_HET_BAC   =   MODEL_CONSTANTS( 35)  ! 35! Aerobic heterotrophic bacteria ON oxidation rate mg N per mg C of bacteria production
! 
!                        KHS_MIN_N    =   MODEL_CONSTANTS(36 )  ! 36! Aerobic heterotrophic bacteria ON mineralisation reverse half saturation for DIN
!              OX_ORGP_AER_HET_BAC    =   MODEL_CONSTANTS(37 )  ! 37! Aerobic heterotrophic bacteria OP mineralisation rate mg P per mg C of bacteria production
!                        KHS_MIN_P    =   MODEL_CONSTANTS(38 )  ! 38! Aerobic heterotrophic bacteria ON mineralisation reverse half saturation for DIP
!        DO_STR_HYPOX_AER_HET_BAC_D   =   MODEL_CONSTANTS(39 )  ! 39! Aerobic heterotrophic bacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!         THETA_HYPOX_AER_HET_BAC_D   =   MODEL_CONSTANTS(40 )  ! 40! Aerobic heterotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
!         EXPON_HYPOX_AER_HET_BAC_D   =   MODEL_CONSTANTS(41 )  ! 41! Aerobic heterotrophic bacteria Exponent constant for Dissolved oxygen stress
!                AER_HET_BAC_N_TO_C   =   MODEL_CONSTANTS(42 )  ! 42! Aerobic heterotrophic bacteria Nitrogen to Carbon ratio
!                AER_HET_BAC_P_TO_C   =   MODEL_CONSTANTS(43 )  ! 43! Aerobic heterotrophic bacteria Phosphorus to Carbon ratio
!               AER_HET_BAC_O2_TO_C   =   MODEL_CONSTANTS(44 )  ! 44! Aerobic heterotrophic bacteria Oxygen to Carbon ratio for respiration
!              KG_FAC_AN_HET_BAC_20   =   MODEL_CONSTANTS(45 )  ! 45! Facultative anaerobic heterotrophic bacteria Growth rate of
!         EFF_FAC_AN_HET_BAC_GROWTH   =   MODEL_CONSTANTS(46 )  ! 46! not used! Facultative anaerobic heterotrophic bacteria growth efficiency
!           THETA_KG_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(47 )  ! 47! not used! Facultative anaerobic heterotrophic bacteria Temperature correction for growth rate
!              KR_FAC_AN_HET_BAC_20   =   MODEL_CONSTANTS(48 )  ! 48! not used! Facultative anaerobic heterotrophic bacteria Respiration rate
!           THETA_KR_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(49 )  ! 49! not used! Facultative anaerobic heterotrophic bacteria Temperature correction for respiration rate
!              KD_FAC_AN_HET_BAC_20   =   MODEL_CONSTANTS(50 )  ! 50! not used! Facultative anaerobic heterotrophic bacteria Mortality rate
!           THETA_KD_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(51 )  ! 51! not used! Facultative anaerobic heterotrophic bacteria Temperature correction for Mortality rate
!           KHS_NO3N_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(52 )  ! 52! Facultative anaerobic heterotrophic bacteria Half saturation growth for NO3N
!           KHS_ORGC_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(53 )  ! 53! not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for OC
!           KHS_ORGN_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(54 )  ! 54! not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for ON
!           KHS_ORGP_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(55 )  ! 55! not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for OP
!         REV_KHS_O2_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(56 )  ! 56! not used! Facultative anaerobic heterotrophic bacteria Reverse Half saturation growth for O2
!    NO3N_LACK_STR_FAC_AN_HET_BAC_D   =   MODEL_CONSTANTS(57 )  ! 57! not used! Facultative anaerobic heterotrophic bacteria NO3N stress concentration
!   THETA_NO3_LACK_FAC_AN_HET_BAC_D   =   MODEL_CONSTANTS(58 )  ! 58! not used! Facultative anaerobic heterotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
!     EXP_NO3_LACK_FAC_AN_HET_BAC_D   =   MODEL_CONSTANTS(59 )  ! 59! not used! Facultative anaerobic heterotrophic bacteria Exponent constant for Dissolved oxygen stress
!             FAC_AN_HET_BAC_N_TO_C   =   MODEL_CONSTANTS(60 )  ! 60! not used! Facultative anaerobic heterotrophic bacteria Nitrogen to Carbon ratio
!             FAC_AN_HET_BAC_P_TO_C   =   MODEL_CONSTANTS(61 )  ! 61! not used! Facultative anaerobic heterotrophic bacteria Phosphorus to Carbon ratio
!            FAC_AN_HET_BAC_O2_TO_C   =   MODEL_CONSTANTS(62 )  ! 62! not used! Facultative anaerobic heterotrophic bacteria Oxygen to Carbon ratio for respiration
!              YIELD_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(63 )  ! 63! Facultative anaerobic heterotrophic bacteria Yield of carbon per unit nitrates nitrogen
!                   KG_DIA_OPT_TEMP   =   MODEL_CONSTANTS(64 )  ! 64! Diatoms Growth rate
!                   DIA_OPT_TEMP_LR   =   MODEL_CONSTANTS(65 )  ! 65! Diatoms optimal temperature lower range
!                   DIA_OPT_TEMP_UR   =   MODEL_CONSTANTS(66 )  ! 66! Diatoms optimal temperature upper range
!                    EFF_DIA_GROWTH   =   MODEL_CONSTANTS(67 )  ! 67! Diatoms Effective growth. (1-EG)*growth - losses for respiration and excretion
!          KAPPA_DIA_UNDER_OPT_TEMP   =   MODEL_CONSTANTS(68 )  ! 68! Diatoms Temperature correction for growth lower temperature
!           KAPPA_DIA_OVER_OPT_TEMP   =   MODEL_CONSTANTS(69 )  ! 69! Diatoms Temperature correction for growth upper temperature
!                         KR_DIA_20   =   MODEL_CONSTANTS(70 )  ! 70! Diatoms Respiration rate
!                      THETA_KR_DIA   =   MODEL_CONSTANTS(71 )  ! 71! Diatoms Temperature correction for basal respiration rate
!                         KD_DIA_20   =   MODEL_CONSTANTS(72 )  ! 72! Diatoms Mortality rate
!                      THETA_KD_DIA   =   MODEL_CONSTANTS(73 )  ! 73! Diatoms Temperature correction for Mortality rate
!                       KHS_DIN_DIA   =   MODEL_CONSTANTS(74 )  ! 74! Diatoms Half saturation growth for DIN
!                       KHS_DIP_DIA   =   MODEL_CONSTANTS(75 )  ! 75! Diatoms Half saturation growth for DIP
!                       KHS_DSi_DIA   =   MODEL_CONSTANTS(76 )  ! 76! Diatoms Half saturation growth for DSi
!                        KHS_O2_DIA   =   MODEL_CONSTANTS(77 )  ! 77! Diatoms Half saturation growth for O2
!                     FRAC_DIA_EXCR   =   MODEL_CONSTANTS(78 )  ! 78! Diatoms Fraction of excretion in metabolism rate
!                           I_S_DIA   =   MODEL_CONSTANTS(79 )  ! 79! Diatoms Light saturation (langleys)
!                DO_STR_HYPOX_DIA_D   =   MODEL_CONSTANTS(80 )  ! 80! Diatoms Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!                 THETA_HYPOX_DIA_D   =   MODEL_CONSTANTS(81 )  ! 81! Diatoms Multiplier of the exponent for Dissolved oxygen stress
!                 EXPON_HYPOX_DIA_D   =   MODEL_CONSTANTS(82 )  ! 82! Diatoms Exponent constant for Dissolved oxygen stress
!                        DIA_N_TO_C   =   MODEL_CONSTANTS(83 )  ! 83! Diatoms Nitrogen to Carbon ratio
!                        DIA_P_TO_C   =   MODEL_CONSTANTS(84 )  ! 84! Diatoms Phosphorus to Carbon ratio
!                       DIA_Si_TO_C   =   MODEL_CONSTANTS(85 )  ! 85! Diatoms Silica to Carbon ratio
!                       DIA_O2_TO_C   =   MODEL_CONSTANTS(86 )  ! 86! Diatoms Oxygen to Carbon ratio for respiration
!                     DIA_C_TO_CHLA   =   MODEL_CONSTANTS(87 )  ! 87! Diatoms Carbon to Chlorophil a ratio
!                   KG_CYN_OPT_TEMP   =   MODEL_CONSTANTS(88 )  ! 88! Non-fixing cyanobacteria Growth rate
!                   CYN_OPT_TEMP_LR   =   MODEL_CONSTANTS(89 )  ! 89! Non-fixing cyanobacteria optimal temperature lower range
!                   CYN_OPT_TEMP_UR   =   MODEL_CONSTANTS(90 )  ! 90! Non-fixing cyanobacteria optimal temperature upper range
!                    EFF_CYN_GROWTH   =   MODEL_CONSTANTS(91 )  ! 91! Non-fixing cyanobacteria Effective growth. (1-EG)*growth - losses for respiration and excretion
!          KAPPA_CYN_UNDER_OPT_TEMP   =   MODEL_CONSTANTS(92 )  ! 92! Non-fixing cyanobacteria Temperature correction for growth lower temperature
!           KAPPA_CYN_OVER_OPT_TEMP   =   MODEL_CONSTANTS(93 )  ! 93! Non-fixing cyanobacteria Temperature correction for growth upper temperature
!                         KR_CYN_20   =   MODEL_CONSTANTS(94 )  ! 94! Non-fixing cyanobacteria Respiration rate
!                      THETA_KR_CYN   =   MODEL_CONSTANTS(95 )  ! 95! Non-fixing cyanobacteria Temperature correction for respiration rate
!                         KD_CYN_20   =   MODEL_CONSTANTS(96 )  ! 96! Non-fixing cyanobacteria Mortality rate
!                      THETA_KD_CYN   =   MODEL_CONSTANTS(97 )  ! 97! Non-fixing cyanobacteria Temperature correction for Mortality rate
!                       KHS_DIN_CYN   =   MODEL_CONSTANTS(98 )  ! 98! Non-fixing cyanobacteria Half saturation growth for DIN
!                       KHS_DIP_CYN   =   MODEL_CONSTANTS(99 )  ! 99! Non-fixing cyanobacteria Half saturation growth for DIP
!                        KHS_O2_CYN   =   MODEL_CONSTANTS(100)  !100! Non-fixing cyanobacteria Half saturation growth for O2
!                     FRAC_CYN_EXCR   =   MODEL_CONSTANTS(101)  !101! Non-fixing cyanobacteria Fraction of excretion in metabolism rate
!                           I_S_CYN   =   MODEL_CONSTANTS(102)  !102! Non-fixing cyanobacteria Light saturation (langleys)
!                DO_STR_HYPOX_CYN_D   =   MODEL_CONSTANTS(103)  !103! Non-fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!                 THETA_HYPOX_CYN_D   =   MODEL_CONSTANTS(104)  !104! Non-fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
!                 EXPON_HYPOX_CYN_D   =   MODEL_CONSTANTS(105)  !105! Non-fixing cyanobacteria Exponent constant for Dissolved oxygen stress
!                        CYN_N_TO_C   =   MODEL_CONSTANTS(106)  !106! Non-fixing cyanobacteria Nitrogen to Carbon ratio ,was 0.1
!                        CYN_P_TO_C   =   MODEL_CONSTANTS(107)  !107! Non-fixing cyanobacteria Phosphorus to Carbon ratio
!                       CYN_O2_TO_C   =   MODEL_CONSTANTS(108)  !108! Non-fixing cyanobacteria Oxygen to Carbon ratio for respiration
!                     CYN_C_TO_CHLA   =   MODEL_CONSTANTS(109)  !109! Non-fixing cyanobacteria Carbon to Chlorophyl a ratio
!               KG_FIX_CYN_OPT_TEMP   =   MODEL_CONSTANTS(110)  !110! Fixing cyanobacteria Growth rate
!               FIX_CYN_OPT_TEMP_LR   =   MODEL_CONSTANTS(111)  !111! Fixing Cyanobacteria optimal temperature lower range
!               FIX_CYN_OPT_TEMP_UR   =   MODEL_CONSTANTS(112)  !112! Fixing Cyanobacteria optimal temperature upper range
!                EFF_FIX_CYN_GROWTH   =   MODEL_CONSTANTS(113)  !113! Fixing cyanobacteria Effective growth. (1-EG)*growth - losses for RESP and excretion
!      KAPPA_FIX_CYN_UNDER_OPT_TEMP   =   MODEL_CONSTANTS(114)  !114! Fixing cyanobacteria Temperature correction for growth lower temperature
!       KAPPA_FIX_CYN_OVER_OPT_TEMP   =   MODEL_CONSTANTS(115)  !115! Fixing cyanobacteria Temperature correction for growth upper temperature
!                     KR_FIX_CYN_20   =   MODEL_CONSTANTS(116)  !116! Fixing cyanobacteria RESP rate
!                  THETA_KR_FIX_CYN   =   MODEL_CONSTANTS(117)  !117! Fixing cyanobacteria Temperature correction for RESP rate
!                     KD_FIX_CYN_20   =   MODEL_CONSTANTS(118)  !118! Fixing cyanobacteria Mortality rate of nitrification bacteria
!                  THETA_KD_FIX_CYN   =   MODEL_CONSTANTS(119)  !119! Fixing cyanobacteria Temperature correction for Mortality rate
!                   KHS_DIN_FIX_CYN   =   MODEL_CONSTANTS(120)  !120! Fixing cyanobacteria Half saturation growth for DIN
!                   KHS_DIP_FIX_CYN   =   MODEL_CONSTANTS(121)  !121! Fixing cyanobacteria Half saturation growth for DIP
!                    KHS_O2_FIX_CYN   =   MODEL_CONSTANTS(122)  !122! Fixing cyanobacteria Half saturation growth for O2
!                 FRAC_FIX_CYN_EXCR   =   MODEL_CONSTANTS(123)  !123! Fixing cyanobacteria Fraction of excretion in metabolism rate
!                       I_S_FIX_CYN   =   MODEL_CONSTANTS(124)  !124! Fixing cyanobacteria Light saturation (langleys)
!            DO_STR_HYPOX_FIX_CYN_D   =   MODEL_CONSTANTS(125)  !125! Fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!             THETA_HYPOX_FIX_CYN_D   =   MODEL_CONSTANTS(126)  !126! Fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
!             EXPON_HYPOX_FIX_CYN_D   =   MODEL_CONSTANTS(127)  !127! Fixing cyanobacteria Exponent constant for Dissolved oxygen stress
!                    FIX_CYN_N_TO_C   =   MODEL_CONSTANTS(128)  !128! Fixing cyanobacteria Nitrogen to Carbon ratio
!                    FIX_CYN_P_TO_C   =   MODEL_CONSTANTS(129)  !129! Fixing cyanobacteria Phosphorus to Carbon ratio
!                   FIX_CYN_O2_TO_C   =   MODEL_CONSTANTS(130)  !130! Fixing cyanobacteria Oxygen to Carbon ratio for respiration
!                 FIX_CYN_C_TO_CHLA   =   MODEL_CONSTANTS(131)  !131! Fixing cyanobacteria Carbon to Chlorophyl a ratio
!                             R_FIX   =   MODEL_CONSTANTS(132)  !132! Fixing cyanobacteria Ratio between non-fixing and fixing fractions growth rate
!                             K_FIX   =   MODEL_CONSTANTS(133)  !133! Fixing cyanobacteria Effectivity parameter of switching to nitrogen fixation
!                   KG_OPA_OPT_TEMP   =   MODEL_CONSTANTS(134)  !134! OtherPhyto Growth rate
!                   OPA_OPT_TEMP_LR   =   MODEL_CONSTANTS(135)  !135! OtherPhyto optimal temperature lower range
!                   OPA_OPT_TEMP_UR   =   MODEL_CONSTANTS(136)  !136! OtherPhyto optimal temperature upper range
!                    EFF_OPA_GROWTH   =   MODEL_CONSTANTS(137)  !137! OtherPhyto Effective growth. (1-EG)*growth - losses for respiration and excretion
!          KAPPA_OPA_UNDER_OPT_TEMP   =   MODEL_CONSTANTS(138)  !138! OtherPhyto Temperature correction for growth lower temperature
!           KAPPA_OPA_OVER_OPT_TEMP   =   MODEL_CONSTANTS(139)  !139! OtherPhyto Temperature correction for growth upper temperature
!                         KR_OPA_20   =   MODEL_CONSTANTS(140)  !140! OtherPhyto Respiration rate
!                      THETA_KR_OPA   =   MODEL_CONSTANTS(141)  !141! OtherPhyto Temperature correction for respiration rate
!                         KD_OPA_20   =   MODEL_CONSTANTS(142)  !142! OtherPhyto Mortality rate
!                      THETA_KD_OPA   =   MODEL_CONSTANTS(143)  !143! OtherPhyto Temperature correction for Mortality rate
!                       KHS_DIN_OPA   =   MODEL_CONSTANTS(144)  !144! OtherPhyto Half saturation growth for DIN
!                       KHS_DIP_OPA   =   MODEL_CONSTANTS(145)  !145! OtherPhyto Half saturation growth for DIP
!                        KHS_O2_OPA   =   MODEL_CONSTANTS(146)  !146! OtherPhyto Half saturation growth for O2
!                     FRAC_OPA_EXCR   =   MODEL_CONSTANTS(147)  !147! OtherPhyto Fraction of excretion in metabolism rate
!                           I_S_OPA   =   MODEL_CONSTANTS(148)  !148! OtherPhyto Light saturation (langleys)
!                DO_STR_HYPOX_OPA_D   =   MODEL_CONSTANTS(149)  !149! OtherPhyto Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!                 THETA_HYPOX_OPA_D   =   MODEL_CONSTANTS(150)  !150! OtherPhyto Multiplier of the exponent for Dissolved oxygen stress
!                 EXPON_HYPOX_OPA_D   =   MODEL_CONSTANTS(151)  !151! OtherPhyto Exponent constant for Dissolved oxygen stress
!                        OPA_N_TO_C   =   MODEL_CONSTANTS(152)  !152! OtherPhyto Nitrogen to Carbon ratio
!                        OPA_P_TO_C   =   MODEL_CONSTANTS(153)  !153! OtherPhyto Phosphorus to Carbon ratio
!                       OPA_O2_TO_C   =   MODEL_CONSTANTS(154)  !154! OtherPhyto Oxygen to Carbon ratio for respiration
!                     OPA_C_TO_CHLA   =   MODEL_CONSTANTS(155)  !155! OtherPhyto Carbon to Chlorophyl a ratio
!                   KG_ZOO_OPT_TEMP   =   MODEL_CONSTANTS(156)  !156! Zooplankton Growth rate
!                   ZOO_OPT_TEMP_LR   =   MODEL_CONSTANTS(157)  !157! Zooplankton optimal temperature lower range
!                   ZOO_OPT_TEMP_UR   =   MODEL_CONSTANTS(158)  !158! Zooplankton optimal temperature upper range
!                    EFF_ZOO_GROWTH   =   MODEL_CONSTANTS(159)  !159! Zooplankton Effective growth. (1-EG)*growth - losses for respiration and excretion
!          KAPPA_ZOO_UNDER_OPT_TEMP   =   MODEL_CONSTANTS(160)  !160! Zooplankton Temperature correction for growth lower temperature
!           KAPPA_ZOO_OVER_OPT_TEMP   =   MODEL_CONSTANTS(161)  !161! Zooplankton Temperature correction for growth upper temperature
!                      GRAT_ZOO_DIA   =   MODEL_CONSTANTS(162)  !162! Zooplankton Grazing rate (growhth rate multiplier) on diatoms
!                      GRAT_ZOO_CYN   =   MODEL_CONSTANTS(163)  !163! Zooplankton Grazing rate (growhth rate multiplier) on Cyanobacteria
!                      GRAT_ZOO_OPA   =   MODEL_CONSTANTS(164)  !164! Zooplankton Grazing rate (growhth rate multiplier) on fixing Cyanobacteria
!                  GRAT_ZOO_FIX_CYN   =   MODEL_CONSTANTS(165)  !165! Zooplankton Grazing rate (growhth rate multiplier) on OtherPhyto
!             GRAT_ZOO_CHEM_AUT_BAC   =   MODEL_CONSTANTS(166)  !166! Zooplankton Grazing rate (growhth rate multiplier) on NITR_BAC
!              GRAT_ZOO_AER_HET_BAC   =   MODEL_CONSTANTS(167)  !167! Zooplankton Grazing rate (growhth rate multiplier) on AER_HET_BAC
!           GRAT_ZOO_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(168)  !168! Zooplankton Grazing rate (growhth rate multiplier) on DENITR_BAC
!           GRAT_ZOO_DET_PART_ORG_C   =   MODEL_CONSTANTS(169)  !169! Zooplankton Grazing rate (growhth rate multiplier) on part. ORG_C
!                      PREF_ZOO_DIA   =   MODEL_CONSTANTS(170)  !170! Zooplankton Preference for diatoms
!                      PREF_ZOO_CYN   =   MODEL_CONSTANTS(171)  !171! Zooplankton Preference for Cyanobacteria
!                  PREF_ZOO_FIX_CYN   =   MODEL_CONSTANTS(172)  !172! Zooplankton Preference for fixing Cyanobacteria
!                      PREF_ZOO_OPA   =   MODEL_CONSTANTS(173)  !173! Zooplankton Preference for OtherPhyto
!             PREF_ZOO_CHEM_AUT_BAC   =   MODEL_CONSTANTS(174)  !174! Zooplankton Preference for NITR_BAC
!              PREF_ZOO_AER_HET_BAC   =   MODEL_CONSTANTS(175)  !175! Zooplankton Preference for AER_HET_BAC
!           PREF_ZOO_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(176)  !176! Zooplankton Preference for DENITR_BAC
!           PREF_ZOO_DET_PART_ORG_C   =   MODEL_CONSTANTS(177)  !177! Zooplankton Preference for part. ORG_C
!                     KHS_DIA_C_ZOO   =   MODEL_CONSTANTS(178)  !178! Zooplankton Half saturation growth for diatoms
!                     KHS_CYN_C_ZOO   =   MODEL_CONSTANTS(179)  !179! Zooplankton Half saturation growth for Cyanobacteria
!                 KHS_FIX_CYN_C_ZOO   =   MODEL_CONSTANTS(180)  !180! Zooplankton Half saturation growth for fixing Cyanobacteria
!                     KHS_OPA_C_ZOO   =   MODEL_CONSTANTS(181)  !181! Zooplankton Half saturation growth for OtherPhyto
!            KHS_CHEM_AUT_BAC_C_ZOO   =   MODEL_CONSTANTS(182)  !182! Zooplankton Half saturation growth for NITR_BAC
!             KHS_AER_HET_BAC_C_ZOO   =   MODEL_CONSTANTS(183)  !183! Zooplankton Half saturation growth for AER_HET_BAC
!          KHS_FAC_AN_HET_BAC_C_ZOO   =   MODEL_CONSTANTS(184)  !184! Zooplankton Half saturation growth for DENITR_BAC
!            KHS_DET_PART_ORG_C_ZOO   =   MODEL_CONSTANTS(185)  !185! Zooplankton Half saturation growth for part. ORG_C
!                      FOOD_MIN_ZOO   =   MODEL_CONSTANTS(186)  !186! Zooplankton Minimum food conc. for feeding
!                            KE_ZOO   =   MODEL_CONSTANTS(187)  !187! not used Zooplankton Excretion rate as growth fraction
!                   FRAC_ZOO_EX_ORG   =   MODEL_CONSTANTS(188)  !188! not used Zooplankton Excretion rate organic fraction
!                         KR_ZOO_20   =   MODEL_CONSTANTS(189)  !189! Zooplankton Respiration rate
!                      THETA_KR_ZOO   =   MODEL_CONSTANTS(190)  !190! Zooplankton Respiration rate Temperature correction
!                         KD_ZOO_20   =   MODEL_CONSTANTS(191)  !191! Zooplankton Mortality rate
!                      THETA_KD_ZOO   =   MODEL_CONSTANTS(192)  !192! Zooplankton Mortality rate Temperature correction
!                DO_STR_HYPOX_ZOO_D   =   MODEL_CONSTANTS(193)  !193! Zooplankton Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!                 THETA_HYPOX_ZOO_D   =   MODEL_CONSTANTS(194)  !194! Zooplankton Multiplier of the exponent for Dissolved oxygen stress
!                 EXPON_HYPOX_ZOO_D   =   MODEL_CONSTANTS(195)  !195! Zooplankton Exponent constant for Dissolved oxygen stress
!                        ZOO_N_TO_C   =   MODEL_CONSTANTS(196)  !196! Zooplankton Nitrogen to Carbon ratio
!                        ZOO_P_TO_C   =   MODEL_CONSTANTS(197)  !197! Zooplankton Phosphorus to Carbon ratio
!                       ZOO_O2_TO_C   =   MODEL_CONSTANTS(198)  !198! Zooplankton Oxygen to Carbon ratio for respiration
!           KDISS_DET_PART_ORG_C_20   =   MODEL_CONSTANTS(199)  !199! Particulate Detritus Carbon Dissolution rate not dependent on phytoplankton
!        THETA_KDISS_DET_PART_ORG_C   =   MODEL_CONSTANTS(200)  !200! Particulate Detritus Carbon Dissolution rate Temperature correction
!           FAC_PHYT_DET_PART_ORG_C   =   MODEL_CONSTANTS(201)  !201! Particulate Detritus Carbon Phytoplankton linear factor for dissolution rate
!           KDISS_DET_PART_ORG_N_20   =   MODEL_CONSTANTS(202)  !202! Particulate Detritus Nitrogen Dissolution rate not dependent on phytoplankton
!        THETA_KDISS_DET_PART_ORG_N   =   MODEL_CONSTANTS(203)  !203! Particulate Detritus Nitrogen Dissolution rate Temperature correction
!                        KHS_DISS_N   =   MODEL_CONSTANTS(204)  !204! Particulate Detritus Nitrogen dissolution reverse half saturation for DIN
!           FAC_PHYT_DET_PART_ORG_N   =   MODEL_CONSTANTS(205)  !205! Particulate Detritus Nitrogen Phytoplankton linear factor for dissolution rate
!           KDISS_DET_PART_ORG_P_20   =   MODEL_CONSTANTS(206)  !206! Particulate Detritus Phosphorus Dissolution rate not dependent on phytoplankton
!        THETA_KDISS_DET_PART_ORG_P   =   MODEL_CONSTANTS(207)  !207! Particulate Detritus Phosphorus Dissolution rate Temperature correction
!                        KHS_DISS_P   =   MODEL_CONSTANTS(208)  !208! Particulate Detritus Phosphorus  dissolution reverse half saturation for DIP
!           FAC_PHYT_DET_PART_ORG_P   =   MODEL_CONSTANTS(209)  !209! Particulate Detritus Phosphorus  Phytoplankton linear factor for dissolution rate
!                  KDISS_PART_Si_20   =   MODEL_CONSTANTS(210)  !210! Particulate Silica Dissolution rate
!               THETA_KDISS_PART_Si   =   MODEL_CONSTANTS(211)  !211! Particulate Silica Dissolution rate Temperature correction
!                      K_MIN_DOC_20   =   MODEL_CONSTANTS(212)  !212! Dissolved carbon  mineralisation rate
!                   THETA_K_MIN_DOC   =   MODEL_CONSTANTS(213)  !213! Dissolved carbon  mineralisation rate Temperature constant
!                 FAC_PHYT_AMIN_DOC   =   MODEL_CONSTANTS(214)  !214! Dissolved carbon  Phytoplankton linear factor for mineralisation rate
!                      K_MIN_DON_20   =   MODEL_CONSTANTS(215)  !215! Dissolved nitrogen  mineralisation rate not dependent on phytoplankton
!                   THETA_K_MIN_DON   =   MODEL_CONSTANTS(216)  !216! Dissolved nitrogen  mineralisation rate Temperature constant
!                        KHS_AMIN_N   =   MODEL_CONSTANTS(217)  !217! Dissolved nitrogen  reverse half saturation for DIN
!                 FAC_PHYT_AMIN_DON   =   MODEL_CONSTANTS(218)  !218! Dissolved nitrogen Phytoplankton linear factor for mineralisation rate
!                      K_MIN_DOP_20   =   MODEL_CONSTANTS(219)  !219! Dissolved phosphorus  mineralisation rate not dependent on phytoplankton
!                   THETA_K_MIN_DOP   =   MODEL_CONSTANTS(220)  !220! Dissolved phosphorus  mineralisation rate Temperature constant
!                        KHS_AMIN_P   =   MODEL_CONSTANTS(221)  !221! Dissolved phosphorus reverse half saturation for DIP
!                 FAC_PHYT_AMIN_DOP   =   MODEL_CONSTANTS(222)  !222! Dissolved phosphorus Phytoplankton linear factor for mineralisation rate
!                         K_NITR_20   =   MODEL_CONSTANTS(223)  !223! Amonia nitrification rate
!                      KHS_NITR_OXY   =   MODEL_CONSTANTS(224)  !224! Amonia nitrification half saturation for Oxygen
!                    KHS_NITR_NH4_N   =   MODEL_CONSTANTS(225)  !225! Amonia nitrification half saturation for Amonia
!                      THETA_K_NITR   =   MODEL_CONSTANTS(226)  !226! Amonia nitrification rate Temperature constant

     ! When changed this should be also updated in derived_vars
     call para_get_value('K_A'                             ,                              K_A) !  1 Aeration coefficient (if negative calculates internally)
     call para_get_value('THETA_K_A'                       ,                        THETA_K_A) !  2 Temperature correction factor for aeration
     call para_get_value('KG_CHEM_AUT_BAC_20'              ,               KG_CHEM_AUT_BAC_20) !  3 Chemoautotrophic bacteria Growth rate
     call para_get_value('EFF_CHEM_AUT_BAC_GROWTH'         ,          EFF_CHEM_AUT_BAC_GROWTH) !  4 Chemoautotrophic bacteria growth efficiency
     call para_get_value('THETA_KG_CHEM_AUT_BAC'           ,            THETA_KG_CHEM_AUT_BAC) !  5 Chemoautotrophic bacteria Temperature correction for growth rate
     call para_get_value('KR_CHEM_AUT_BAC_20'              ,               KR_CHEM_AUT_BAC_20) !  6 Chemoautotrophic bacteria Respiration rate
     call para_get_value('THETA_KR_CHEM_AUT_BAC'           ,            THETA_KR_CHEM_AUT_BAC) !  7 Chemoautotrophic bacteria Temperature correction for respiration rate
     call para_get_value('KD_CHEM_AUT_BAC_20'              ,               KD_CHEM_AUT_BAC_20) !  8 Chemoautotrophic bacteria Mortality rate
     call para_get_value('THETA_KD_CHEM_AUT_BAC'           ,            THETA_KD_CHEM_AUT_BAC) !  9 Chemoautotrophic bacteria Temperature correction for Mortality rate
     call para_get_value('KHS_NH4N_CHEM_AUT_BAC'           ,            KHS_NH4N_CHEM_AUT_BAC) ! 10 Chemoautotrophic bacteria Half saturation growth for NH4N
     call para_get_value('KHS_PO4P_CHEM_AUT_BAC'           ,            KHS_PO4P_CHEM_AUT_BAC) ! 11 Chemoautotrophic bacteria Half saturation growth for PO4P
     call para_get_value('KHS_O2_CHEM_AUT_BAC'             ,              KHS_O2_CHEM_AUT_BAC) ! 12 Chemoautotrophic bacteria Half saturation growth for O2
     call para_get_value('DO_STR_HYPOX_CHEM_AUT_BAC_D'     ,      DO_STR_HYPOX_CHEM_AUT_BAC_D) ! 13 Chemoautotrophic bacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_CHEM_AUT_BAC_D'      ,       THETA_HYPOX_CHEM_AUT_BAC_D) ! 14 Chemoautotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_CHEM_AUT_BAC_D'      ,       EXPON_HYPOX_CHEM_AUT_BAC_D) ! 15 Chemoautotrophic bacteria Exponent constant for Dissolved oxygen stress
     call para_get_value('CHEM_AUT_BAC_N_TO_C'             ,              CHEM_AUT_BAC_N_TO_C) ! 16 Chemoautotrophic bacteria Nitrogen to Carbon ratio
     call para_get_value('CHEM_AUT_BAC_P_TO_C'             ,              CHEM_AUT_BAC_P_TO_C) ! 17 Chemoautotrophic bacteria Phosphorus to Carbon ratio
     call para_get_value('CHEM_AUT_BAC_O2_TO_C'            ,             CHEM_AUT_BAC_O2_TO_C) ! 18 Chemoautotrophic bacteria Oxygen to Carbon ratio
     call para_get_value('YIELD_CHEM_AUT_BAC'              ,               YIELD_CHEM_AUT_BAC) ! 19 Chemoautotrophic bacteria Yield of Carbon per unit amonia nitrogen
     call para_get_value('KG_AER_HET_BAC_20'               ,                KG_AER_HET_BAC_20) ! 20 Aerobic heterotrophic bacteria Growth rate
     call para_get_value('EFF_AER_HET_BAC_GROWTH'          ,           EFF_AER_HET_BAC_GROWTH) ! 21 Aerobic heterotrophic bacteria growth efficiency
     call para_get_value('THETA_KG_AER_HET_BAC'            ,             THETA_KG_AER_HET_BAC) ! 22 Aerobic heterotrophic bacteria Temperature correction for growth rate
     call para_get_value('KR_AER_HET_BAC_20'               ,                KR_AER_HET_BAC_20) ! 23 Aerobic heterotrophic bacteria Respiration rate
     call para_get_value('THETA_KR_AER_HET_BAC'            ,             THETA_KR_AER_HET_BAC) ! 24 Aerobic heterotrophic bacteria Temperature correction for respiration rate
     call para_get_value('KD_AER_HET_BAC_20'               ,                KD_AER_HET_BAC_20) ! 25 Aerobic heterotrophic bacteria Mortality rate
     call para_get_value('THETA_KD_AER_HET_BAC'            ,             THETA_KD_AER_HET_BAC) ! 26 Aerobic heterotrophic bacteria Temperature correction for Mortality rate
     call para_get_value('KHS_ORGC_AER_HET_BAC'            ,             KHS_ORGC_AER_HET_BAC) ! 27 Aerobic heterotrophic bacteria Half saturation growth for OC
     call para_get_value('KHS_ORGN_AER_HET_BAC'            ,             KHS_ORGN_AER_HET_BAC) ! 28 Aerobic heterotrophic bacteria Half saturation growth for ON
     call para_get_value('KHS_ORGP_AER_HET_BAC'            ,             KHS_ORGP_AER_HET_BAC) ! 29 Aerobic heterotrophic bacteria Half saturation growth for OP
     call para_get_value('KHS_O2_AER_HET_BAC'              ,               KHS_O2_AER_HET_BAC) ! 30 Aerobic heterotrophic bacteria Half saturation growth for Oxygen
     call para_get_value('KHS_DIN_AER_HET_BAC'             ,              KHS_DIN_AER_HET_BAC) ! 31 Aerobic heterotrophic bacteria Half saturation growth for inorganic nitrogen
     call para_get_value('KHS_DIP_AER_HET_BAC'             ,              KHS_DIP_AER_HET_BAC) ! 32 Aerobic heterotrophic bacteria Half saturation growth for inorganic phosphorus
     call para_get_value('KHS_PHYT_AER_HET_BAC'            ,             KHS_PHYT_AER_HET_BAC) ! 33 Aerobic heterotrophic bacteria Half saturation growth for Phytoplankton C (not used as a resource)
     call para_get_value('YIELD_OC_AER_HET_BAC'            ,             YIELD_OC_AER_HET_BAC) ! 34 Aerobic heterotrophic bacteria Yield of bacteria carbon per unit of organic carbon
     call para_get_value('OX_ORGN_AER_HET_BAC'             ,              OX_ORGN_AER_HET_BAC) ! 35 Aerobic heterotrophic bacteria ON oxidation rate mg N per mg C of bacteria production
     call para_get_value('KHS_MIN_N'                       ,                       KHS_MIN_N ) ! 36 Aerobic heterotrophic bacteria ON mineralisation reverse half saturation for DIN
     call para_get_value('OX_ORGP_AER_HET_BAC'             ,             OX_ORGP_AER_HET_BAC ) ! 37 Aerobic heterotrophic bacteria OP mineralisation rate mg P per mg C of bacteria production
     call para_get_value('KHS_MIN_P'                       ,                       KHS_MIN_P ) ! 38 Aerobic heterotrophic bacteria ON mineralisation reverse half saturation for DIP
     call para_get_value('DO_STR_HYPOX_AER_HET_BAC_D'      ,       DO_STR_HYPOX_AER_HET_BAC_D) ! 39 Aerobic heterotrophic bacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_AER_HET_BAC_D'       ,        THETA_HYPOX_AER_HET_BAC_D) ! 40 Aerobic heterotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_AER_HET_BAC_D'       ,        EXPON_HYPOX_AER_HET_BAC_D) ! 41 Aerobic heterotrophic bacteria Exponent constant for Dissolved oxygen stress
     call para_get_value('AER_HET_BAC_N_TO_C'              ,               AER_HET_BAC_N_TO_C) ! 42 Aerobic heterotrophic bacteria Nitrogen to Carbon ratio
     call para_get_value('AER_HET_BAC_P_TO_C'              ,               AER_HET_BAC_P_TO_C) ! 43 Aerobic heterotrophic bacteria Phosphorus to Carbon ratio
     call para_get_value('AER_HET_BAC_O2_TO_C'             ,              AER_HET_BAC_O2_TO_C) ! 44 Aerobic heterotrophic bacteria Oxygen to Carbon ratio for respiration
     call para_get_value('KG_FAC_AN_HET_BAC_20'            ,             KG_FAC_AN_HET_BAC_20) ! 45 Facultative anaerobic heterotrophic bacteria Growth rate of
     call para_get_value('EFF_FAC_AN_HET_BAC_GROWTH'       ,        EFF_FAC_AN_HET_BAC_GROWTH) ! 46 not used! Facultative anaerobic heterotrophic bacteria growth efficiency
     call para_get_value('THETA_KG_FAC_AN_HET_BAC'         ,          THETA_KG_FAC_AN_HET_BAC) ! 47 not used! Facultative anaerobic heterotrophic bacteria Temperature correction for growth rate
     call para_get_value('KR_FAC_AN_HET_BAC_20'            ,             KR_FAC_AN_HET_BAC_20) ! 48 not used! Facultative anaerobic heterotrophic bacteria Respiration rate
     call para_get_value('THETA_KR_FAC_AN_HET_BAC'         ,          THETA_KR_FAC_AN_HET_BAC) ! 49 not used! Facultative anaerobic heterotrophic bacteria Temperature correction for respiration rate
     call para_get_value('KD_FAC_AN_HET_BAC_20'            ,             KD_FAC_AN_HET_BAC_20) ! 50 not used! Facultative anaerobic heterotrophic bacteria Mortality rate
     call para_get_value('THETA_KD_FAC_AN_HET_BAC'         ,          THETA_KD_FAC_AN_HET_BAC) ! 51 not used! Facultative anaerobic heterotrophic bacteria Temperature correction for Mortality rate
     call para_get_value('KHS_NO3N_FAC_AN_HET_BAC'         ,          KHS_NO3N_FAC_AN_HET_BAC) ! 52 Facultative anaerobic heterotrophic bacteria Half saturation growth for NO3N
     call para_get_value('KHS_ORGC_FAC_AN_HET_BAC'         ,          KHS_ORGC_FAC_AN_HET_BAC) ! 53 not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for OC
     call para_get_value('KHS_ORGN_FAC_AN_HET_BAC'         ,          KHS_ORGN_FAC_AN_HET_BAC) ! 54 not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for ON
     call para_get_value('KHS_ORGP_FAC_AN_HET_BAC'         ,          KHS_ORGP_FAC_AN_HET_BAC) ! 55 not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for OP
     call para_get_value('REV_KHS_O2_FAC_AN_HET_BAC'       ,        REV_KHS_O2_FAC_AN_HET_BAC) ! 56 not used! Facultative anaerobic heterotrophic bacteria Reverse Half saturation growth for O2
     call para_get_value('NO3N_LACK_STR_FAC_AN_HET_BAC_D'  ,   NO3N_LACK_STR_FAC_AN_HET_BAC_D) ! 57 not used! Facultative anaerobic heterotrophic bacteria NO3N stress concentration
     call para_get_value('THETA_NO3_LACK_FAC_AN_HET_BAC_D' ,  THETA_NO3_LACK_FAC_AN_HET_BAC_D) ! 58 not used! Facultative anaerobic heterotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXP_NO3_LACK_FAC_AN_HET_BAC_D'   ,    EXP_NO3_LACK_FAC_AN_HET_BAC_D) ! 59 not used! Facultative anaerobic heterotrophic bacteria Exponent constant for Dissolved oxygen stress
     call para_get_value('FAC_AN_HET_BAC_N_TO_C'           ,            FAC_AN_HET_BAC_N_TO_C) ! 60 not used! Facultative anaerobic heterotrophic bacteria Nitrogen to Carbon ratio
     call para_get_value('FAC_AN_HET_BAC_P_TO_C'           ,            FAC_AN_HET_BAC_P_TO_C) ! 61 not used! Facultative anaerobic heterotrophic bacteria Phosphorus to Carbon ratio
     call para_get_value('FAC_AN_HET_BAC_O2_TO_C'          ,           FAC_AN_HET_BAC_O2_TO_C) ! 62 not used! Facultative anaerobic heterotrophic bacteria Oxygen to Carbon ratio for respiration
     call para_get_value('YIELD_FAC_AN_HET_BAC'            ,             YIELD_FAC_AN_HET_BAC) ! 63 Facultative anaerobic heterotrophic bacteria Yield of carbon per unit nitrates nitrogen
     call para_get_value('KG_DIA_OPT_TEMP'                 ,                  KG_DIA_OPT_TEMP) ! 64 Diatoms Growth rate
     call para_get_value('DIA_OPT_TEMP_LR'                 ,                  DIA_OPT_TEMP_LR) ! 65 Diatoms optimal temperature lower range
     call para_get_value('DIA_OPT_TEMP_UR'                 ,                  DIA_OPT_TEMP_UR) ! 66 Diatoms optimal temperature upper range
     call para_get_value('EFF_DIA_GROWTH'                  ,                   EFF_DIA_GROWTH) ! 67 Diatoms Effective growth. (1-EG)*growth - losses for respiration and excretion
     call para_get_value('KAPPA_DIA_UNDER_OPT_TEMP'        ,         KAPPA_DIA_UNDER_OPT_TEMP) ! 68 Diatoms Temperature correction for growth lower temperature
     call para_get_value('KAPPA_DIA_OVER_OPT_TEMP'         ,          KAPPA_DIA_OVER_OPT_TEMP) ! 69 Diatoms Temperature correction for growth upper temperature
     call para_get_value('KR_DIA_20'                       ,                        KR_DIA_20) ! 70 Diatoms Respiration rate
     call para_get_value('THETA_KR_DIA'                    ,                     THETA_KR_DIA) ! 71 Diatoms Temperature correction for basal respiration rate
     call para_get_value('KD_DIA_20'                       ,                        KD_DIA_20) ! 72 Diatoms Mortality rate
     call para_get_value('THETA_KD_DIA'                    ,                     THETA_KD_DIA) ! 73 Diatoms Temperature correction for Mortality rate
     call para_get_value('KHS_DIN_DIA'                     ,                      KHS_DIN_DIA) ! 74 Diatoms Half saturation growth for DIN
     call para_get_value('KHS_DIP_DIA'                     ,                      KHS_DIP_DIA) ! 75 Diatoms Half saturation growth for DIP
     call para_get_value('KHS_DSi_DIA'                     ,                      KHS_DSi_DIA) ! 76 Diatoms Half saturation growth for DSi
     call para_get_value('KHS_O2_DIA'                      ,                       KHS_O2_DIA) ! 77 Diatoms Half saturation growth for O2
     call para_get_value('FRAC_DIA_EXCR'                   ,                    FRAC_DIA_EXCR) ! 78 Diatoms Fraction of excretion in metabolism rate
     call para_get_value('I_S_DIA'                         ,                          I_S_DIA) ! 79 Diatoms Light saturation (langleys)
     call para_get_value('DO_STR_HYPOX_DIA_D'              ,               DO_STR_HYPOX_DIA_D) ! 80 Diatoms Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_DIA_D'               ,                THETA_HYPOX_DIA_D) ! 81 Diatoms Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_DIA_D'               ,                EXPON_HYPOX_DIA_D) ! 82 Diatoms Exponent constant for Dissolved oxygen stress
     call para_get_value('DIA_N_TO_C'                      ,                       DIA_N_TO_C) ! 83 Diatoms Nitrogen to Carbon ratio
     call para_get_value('DIA_P_TO_C'                      ,                       DIA_P_TO_C) ! 84 Diatoms Phosphorus to Carbon ratio
     call para_get_value('DIA_Si_TO_C'                     ,                      DIA_Si_TO_C) ! 85 Diatoms Silica to Carbon ratio
     call para_get_value('DIA_O2_TO_C'                     ,                      DIA_O2_TO_C) ! 86 Diatoms Oxygen to Carbon ratio for respiration
     call para_get_value('DIA_C_TO_CHLA'                   ,                    DIA_C_TO_CHLA) ! 87 Diatoms Carbon to Chlorophil a ratio
     call para_get_value('KG_CYN_OPT_TEMP'                 ,                  KG_CYN_OPT_TEMP) ! 88 Non-fixing cyanobacteria Growth rate
     call para_get_value('CYN_OPT_TEMP_LR'                 ,                  CYN_OPT_TEMP_LR) ! 89 Non-fixing cyanobacteria optimal temperature lower range
     call para_get_value('CYN_OPT_TEMP_UR'                 ,                  CYN_OPT_TEMP_UR) ! 90 Non-fixing cyanobacteria optimal temperature upper range
     call para_get_value('EFF_CYN_GROWTH'                  ,                   EFF_CYN_GROWTH) ! 91 Non-fixing cyanobacteria Effective growth. (1-EG)*growth - losses for respiration and excretion
     call para_get_value('KAPPA_CYN_UNDER_OPT_TEMP'        ,         KAPPA_CYN_UNDER_OPT_TEMP) ! 92 Non-fixing cyanobacteria Temperature correction for growth lower temperature
     call para_get_value('KAPPA_CYN_OVER_OPT_TEMP'         ,          KAPPA_CYN_OVER_OPT_TEMP) ! 93 Non-fixing cyanobacteria Temperature correction for growth upper temperature
     call para_get_value('KR_CYN_20'                       ,                        KR_CYN_20) ! 94 Non-fixing cyanobacteria Respiration rate
     call para_get_value('THETA_KR_CYN'                    ,                     THETA_KR_CYN) ! 95 Non-fixing cyanobacteria Temperature correction for respiration rate
     call para_get_value('KD_CYN_20'                       ,                        KD_CYN_20) ! 96 Non-fixing cyanobacteria Mortality rate
     call para_get_value('THETA_KD_CYN'                    ,                     THETA_KD_CYN) ! 97 Non-fixing cyanobacteria Temperature correction for Mortality rate
     call para_get_value('KHS_DIN_CYN'                     ,                      KHS_DIN_CYN) ! 98 Non-fixing cyanobacteria Half saturation growth for DIN
     call para_get_value('KHS_DIP_CYN'                     ,                      KHS_DIP_CYN) ! 99 Non-fixing cyanobacteria Half saturation growth for DIP
     call para_get_value('KHS_O2_CYN'                      ,                       KHS_O2_CYN) !100 Non-fixing cyanobacteria Half saturation growth for O2
     call para_get_value('FRAC_CYN_EXCR'                   ,                    FRAC_CYN_EXCR) !101 Non-fixing cyanobacteria Fraction of excretion in metabolism rate
     call para_get_value('I_S_CYN'                         ,                          I_S_CYN) !102 Non-fixing cyanobacteria Light saturation (langleys)
     call para_get_value('DO_STR_HYPOX_CYN_D'              ,               DO_STR_HYPOX_CYN_D) !103 Non-fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_CYN_D'               ,                THETA_HYPOX_CYN_D) !104 Non-fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_CYN_D'               ,                EXPON_HYPOX_CYN_D) !105 Non-fixing cyanobacteria Exponent constant for Dissolved oxygen stress
     call para_get_value('CYN_N_TO_C'                      ,                       CYN_N_TO_C) !106 Non-fixing cyanobacteria Nitrogen to Carbon ratio ,was 0.1
     call para_get_value('CYN_P_TO_C'                      ,                       CYN_P_TO_C) !107 Non-fixing cyanobacteria Phosphorus to Carbon ratio
     call para_get_value('CYN_O2_TO_C'                     ,                      CYN_O2_TO_C) !108 Non-fixing cyanobacteria Oxygen to Carbon ratio for respiration
     call para_get_value('CYN_C_TO_CHLA'                   ,                    CYN_C_TO_CHLA) !109 Non-fixing cyanobacteria Carbon to Chlorophyl a ratio
     call para_get_value('KG_FIX_CYN_OPT_TEMP'             ,              KG_FIX_CYN_OPT_TEMP) !110 Fixing cyanobacteria Growth rate
     call para_get_value('FIX_CYN_OPT_TEMP_LR'             ,              FIX_CYN_OPT_TEMP_LR) !111 Fixing Cyanobacteria optimal temperature lower range
     call para_get_value('FIX_CYN_OPT_TEMP_UR'             ,              FIX_CYN_OPT_TEMP_UR) !112 Fixing Cyanobacteria optimal temperature upper range
     call para_get_value('EFF_FIX_CYN_GROWTH'              ,               EFF_FIX_CYN_GROWTH) !113 Fixing cyanobacteria Effective growth. (1-EG)*growth - losses for RESP and excretion
     call para_get_value('KAPPA_FIX_CYN_UNDER_OPT_TEMP'    ,     KAPPA_FIX_CYN_UNDER_OPT_TEMP) !114 Fixing cyanobacteria Temperature correction for growth lower temperature
     call para_get_value('KAPPA_FIX_CYN_OVER_OPT_TEMP'     ,      KAPPA_FIX_CYN_OVER_OPT_TEMP) !115 Fixing cyanobacteria Temperature correction for growth upper temperature
     call para_get_value('KR_FIX_CYN_20'                   ,                    KR_FIX_CYN_20) !116 Fixing cyanobacteria RESP rate
     call para_get_value('THETA_KR_FIX_CYN'                ,                 THETA_KR_FIX_CYN) !117 Fixing cyanobacteria Temperature correction for RESP rate
     call para_get_value('KD_FIX_CYN_20'                   ,                    KD_FIX_CYN_20) !118 Fixing cyanobacteria Mortality rate of nitrification bacteria
     call para_get_value('THETA_KD_FIX_CYN'                ,                 THETA_KD_FIX_CYN) !119 Fixing cyanobacteria Temperature correction for Mortality rate
     call para_get_value('KHS_DIN_FIX_CYN'                 ,                  KHS_DIN_FIX_CYN) !120 Fixing cyanobacteria Half saturation growth for DIN
     call para_get_value('KHS_DIP_FIX_CYN'                 ,                  KHS_DIP_FIX_CYN) !121 Fixing cyanobacteria Half saturation growth for DIP
     call para_get_value('KHS_O2_FIX_CYN'                  ,                   KHS_O2_FIX_CYN) !122 Fixing cyanobacteria Half saturation growth for O2
     call para_get_value('FRAC_FIX_CYN_EXCR'               ,                FRAC_FIX_CYN_EXCR) !123 Fixing cyanobacteria Fraction of excretion in metabolism rate
     call para_get_value('I_S_FIX_CYN'                     ,                      I_S_FIX_CYN) !124 Fixing cyanobacteria Light saturation (langleys)
     call para_get_value('DO_STR_HYPOX_FIX_CYN_D'          ,           DO_STR_HYPOX_FIX_CYN_D) !125 Fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_FIX_CYN_D'           ,            THETA_HYPOX_FIX_CYN_D) !126 Fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_FIX_CYN_D'           ,            EXPON_HYPOX_FIX_CYN_D) !127 Fixing cyanobacteria Exponent constant for Dissolved oxygen stress
     call para_get_value('FIX_CYN_N_TO_C'                  ,                   FIX_CYN_N_TO_C) !128 Fixing cyanobacteria Nitrogen to Carbon ratio
     call para_get_value('FIX_CYN_P_TO_C'                  ,                   FIX_CYN_P_TO_C) !129 Fixing cyanobacteria Phosphorus to Carbon ratio
     call para_get_value('FIX_CYN_O2_TO_C'                 ,                  FIX_CYN_O2_TO_C) !130 Fixing cyanobacteria Oxygen to Carbon ratio for respiration
     call para_get_value('FIX_CYN_C_TO_CHLA'               ,                FIX_CYN_C_TO_CHLA) !131 Fixing cyanobacteria Carbon to Chlorophyl a ratio
     call para_get_value('R_FIX'                           ,                            R_FIX) !132 Fixing cyanobacteria Ratio between non-fixing and fixing fractions growth rate
     call para_get_value('K_FIX'                           ,                            K_FIX) !133 Fixing cyanobacteria Effectivity parameter of switching to nitrogen fixation
     call para_get_value('KG_OPA_OPT_TEMP'                 ,                  KG_OPA_OPT_TEMP) !134 OtherPhyto Growth rate
     call para_get_value('OPA_OPT_TEMP_LR'                 ,                  OPA_OPT_TEMP_LR) !135 OtherPhyto optimal temperature lower range
     call para_get_value('OPA_OPT_TEMP_UR'                 ,                  OPA_OPT_TEMP_UR) !136 OtherPhyto optimal temperature upper range
     call para_get_value('EFF_OPA_GROWTH'                  ,                   EFF_OPA_GROWTH) !137 OtherPhyto Effective growth. (1-EG)*growth - losses for respiration and excretion
     call para_get_value('KAPPA_OPA_UNDER_OPT_TEMP'        ,         KAPPA_OPA_UNDER_OPT_TEMP) !138 OtherPhyto Temperature correction for growth lower temperature
     call para_get_value('KAPPA_OPA_OVER_OPT_TEMP'         ,          KAPPA_OPA_OVER_OPT_TEMP) !139 OtherPhyto Temperature correction for growth upper temperature
     call para_get_value('KR_OPA_20'                       ,                        KR_OPA_20) !140 OtherPhyto Respiration rate
     call para_get_value('THETA_KR_OPA'                    ,                     THETA_KR_OPA) !141 OtherPhyto Temperature correction for respiration rate
     call para_get_value('KD_OPA_20'                       ,                        KD_OPA_20) !142 OtherPhyto Mortality rate
     call para_get_value('THETA_KD_OPA'                    ,                     THETA_KD_OPA) !143 OtherPhyto Temperature correction for Mortality rate
     call para_get_value('KHS_DIN_OPA'                     ,                      KHS_DIN_OPA) !144 OtherPhyto Half saturation growth for DIN
     call para_get_value('KHS_DIP_OPA'                     ,                      KHS_DIP_OPA) !145 OtherPhyto Half saturation growth for DIP
     call para_get_value('KHS_O2_OPA'                      ,                       KHS_O2_OPA) !146 OtherPhyto Half saturation growth for O2
     call para_get_value('FRAC_OPA_EXCR'                   ,                    FRAC_OPA_EXCR) !147 OtherPhyto Fraction of excretion in metabolism rate
     call para_get_value('I_S_OPA'                         ,                          I_S_OPA) !148 OtherPhyto Light saturation (langleys)
     call para_get_value('DO_STR_HYPOX_OPA_D'              ,               DO_STR_HYPOX_OPA_D) !149 OtherPhyto Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_OPA_D'               ,                THETA_HYPOX_OPA_D) !150 OtherPhyto Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_OPA_D'               ,                EXPON_HYPOX_OPA_D) !151 OtherPhyto Exponent constant for Dissolved oxygen stress
     call para_get_value('OPA_N_TO_C'                      ,                       OPA_N_TO_C) !152 OtherPhyto Nitrogen to Carbon ratio
     call para_get_value('OPA_P_TO_C'                      ,                       OPA_P_TO_C) !153 OtherPhyto Phosphorus to Carbon ratio
     call para_get_value('OPA_O2_TO_C'                     ,                      OPA_O2_TO_C) !154 OtherPhyto Oxygen to Carbon ratio for respiration
     call para_get_value('OPA_C_TO_CHLA'                   ,                    OPA_C_TO_CHLA) !155 OtherPhyto Carbon to Chlorophyl a ratio
     call para_get_value('KG_ZOO_OPT_TEMP'                 ,                  KG_ZOO_OPT_TEMP) !156 Zooplankton Growth rate
     call para_get_value('ZOO_OPT_TEMP_LR'                 ,                  ZOO_OPT_TEMP_LR) !157 Zooplankton optimal temperature lower range
     call para_get_value('ZOO_OPT_TEMP_UR'                 ,                  ZOO_OPT_TEMP_UR) !158 Zooplankton optimal temperature upper range
     call para_get_value('EFF_ZOO_GROWTH'                  ,                   EFF_ZOO_GROWTH) !159 Zooplankton Effective growth. (1-EG)*growth - losses for respiration and excretion
     call para_get_value('KAPPA_ZOO_UNDER_OPT_TEMP'        ,         KAPPA_ZOO_UNDER_OPT_TEMP) !160 Zooplankton Temperature correction for growth lower temperature
     call para_get_value('KAPPA_ZOO_OVER_OPT_TEMP'         ,          KAPPA_ZOO_OVER_OPT_TEMP) !161 Zooplankton Temperature correction for growth upper temperature
     call para_get_value('GRAT_ZOO_DIA'                    ,                     GRAT_ZOO_DIA) !162 Zooplankton Grazing rate (growhth rate multiplier) on diatoms
     call para_get_value('GRAT_ZOO_CYN'                    ,                     GRAT_ZOO_CYN) !163 Zooplankton Grazing rate (growhth rate multiplier) on Cyanobacteria
     call para_get_value('GRAT_ZOO_OPA'                    ,                     GRAT_ZOO_OPA) !164 Zooplankton Grazing rate (growhth rate multiplier) on fixing Cyanobacteria
     call para_get_value('GRAT_ZOO_FIX_CYN'                ,                 GRAT_ZOO_FIX_CYN) !165 Zooplankton Grazing rate (growhth rate multiplier) on OtherPhyto
     call para_get_value('GRAT_ZOO_CHEM_AUT_BAC'           ,            GRAT_ZOO_CHEM_AUT_BAC) !166 Zooplankton Grazing rate (growhth rate multiplier) on NITR_BAC
     call para_get_value('GRAT_ZOO_AER_HET_BAC'            ,             GRAT_ZOO_AER_HET_BAC) !167 Zooplankton Grazing rate (growhth rate multiplier) on AER_HET_BAC
     call para_get_value('GRAT_ZOO_FAC_AN_HET_BAC'         ,          GRAT_ZOO_FAC_AN_HET_BAC) !168 Zooplankton Grazing rate (growhth rate multiplier) on DENITR_BAC
     call para_get_value('GRAT_ZOO_DET_PART_ORG_C'         ,          GRAT_ZOO_DET_PART_ORG_C) !169 Zooplankton Grazing rate (growhth rate multiplier) on part. ORG_C
     call para_get_value('PREF_ZOO_DIA'                    ,                     PREF_ZOO_DIA) !170 Zooplankton Preference for diatoms
     call para_get_value('PREF_ZOO_CYN'                    ,                     PREF_ZOO_CYN) !171 Zooplankton Preference for Cyanobacteria
     call para_get_value('PREF_ZOO_FIX_CYN'                ,                 PREF_ZOO_FIX_CYN) !172 Zooplankton Preference for fixing Cyanobacteria
     call para_get_value('PREF_ZOO_OPA'                    ,                     PREF_ZOO_OPA) !173 Zooplankton Preference for OtherPhyto
     call para_get_value('PREF_ZOO_CHEM_AUT_BAC'           ,            PREF_ZOO_CHEM_AUT_BAC) !174 Zooplankton Preference for NITR_BAC
     call para_get_value('PREF_ZOO_AER_HET_BAC'            ,             PREF_ZOO_AER_HET_BAC) !175 Zooplankton Preference for AER_HET_BAC
     call para_get_value('PREF_ZOO_FAC_AN_HET_BAC'         ,          PREF_ZOO_FAC_AN_HET_BAC) !176 Zooplankton Preference for DENITR_BAC
     call para_get_value('PREF_ZOO_DET_PART_ORG_C'         ,          PREF_ZOO_DET_PART_ORG_C) !177 Zooplankton Preference for part. ORG_C
     call para_get_value('KHS_DIA_C_ZOO'                   ,                    KHS_DIA_C_ZOO) !178 Zooplankton Half saturation growth for diatoms
     call para_get_value('KHS_CYN_C_ZOO'                   ,                    KHS_CYN_C_ZOO) !179 Zooplankton Half saturation growth for Cyanobacteria
     call para_get_value('KHS_FIX_CYN_C_ZOO'               ,                KHS_FIX_CYN_C_ZOO) !180 Zooplankton Half saturation growth for fixing Cyanobacteria
     call para_get_value('KHS_OPA_C_ZOO'                   ,                    KHS_OPA_C_ZOO) !181 Zooplankton Half saturation growth for OtherPhyto
     call para_get_value('KHS_CHEM_AUT_BAC_C_ZOO'          ,           KHS_CHEM_AUT_BAC_C_ZOO) !182 Zooplankton Half saturation growth for NITR_BAC
     call para_get_value('KHS_AER_HET_BAC_C_ZOO'           ,            KHS_AER_HET_BAC_C_ZOO) !183 Zooplankton Half saturation growth for AER_HET_BAC
     call para_get_value('KHS_FAC_AN_HET_BAC_C_ZOO'        ,         KHS_FAC_AN_HET_BAC_C_ZOO) !184 Zooplankton Half saturation growth for DENITR_BAC
     call para_get_value('KHS_DET_PART_ORG_C_ZOO'          ,           KHS_DET_PART_ORG_C_ZOO) !185 Zooplankton Half saturation growth for part. ORG_C
     call para_get_value('FOOD_MIN_ZOO'                    ,                     FOOD_MIN_ZOO) !186 Zooplankton Minimum food conc. for feeding
     call para_get_value('KE_ZOO'                          ,                           KE_ZOO) !187 not used Zooplankton Excretion rate as growth fraction
     call para_get_value('FRAC_ZOO_EX_ORG'                 ,                  FRAC_ZOO_EX_ORG) !188 not used Zooplankton Excretion rate organic fraction
     call para_get_value('KR_ZOO_20'                       ,                        KR_ZOO_20) !189 Zooplankton Respiration rate
     call para_get_value('THETA_KR_ZOO'                    ,                     THETA_KR_ZOO) !190 Zooplankton Respiration rate Temperature correction
     call para_get_value('KD_ZOO_20'                       ,                        KD_ZOO_20) !191 Zooplankton Mortality rate
     call para_get_value('THETA_KD_ZOO'                    ,                     THETA_KD_ZOO) !192 Zooplankton Mortality rate Temperature correction
     call para_get_value('DO_STR_HYPOX_ZOO_D'              ,               DO_STR_HYPOX_ZOO_D) !193 Zooplankton Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_ZOO_D'               ,                THETA_HYPOX_ZOO_D) !194 Zooplankton Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_ZOO_D'               ,                EXPON_HYPOX_ZOO_D) !195 Zooplankton Exponent constant for Dissolved oxygen stress
     call para_get_value('ZOO_N_TO_C'                      ,                       ZOO_N_TO_C) !196 Zooplankton Nitrogen to Carbon ratio
     call para_get_value('ZOO_P_TO_C'                      ,                       ZOO_P_TO_C) !197 Zooplankton Phosphorus to Carbon ratio
     call para_get_value('ZOO_O2_TO_C'                     ,                      ZOO_O2_TO_C) !198 Zooplankton Oxygen to Carbon ratio for respiration
     call para_get_value('KDISS_DET_PART_ORG_C_20'         ,          KDISS_DET_PART_ORG_C_20) !199 Particulate Detritus Carbon Dissolution rate not dependent on phytoplankton
     call para_get_value('THETA_KDISS_DET_PART_ORG_C'      ,       THETA_KDISS_DET_PART_ORG_C) !200 Particulate Detritus Carbon Dissolution rate Temperature correction
     call para_get_value('FAC_PHYT_DET_PART_ORG_C'         ,          FAC_PHYT_DET_PART_ORG_C) !201 Particulate Detritus Carbon Phytoplankton linear factor for dissolution rate
     call para_get_value('KDISS_DET_PART_ORG_N_20'         ,          KDISS_DET_PART_ORG_N_20) !202 Particulate Detritus Nitrogen Dissolution rate not dependent on phytoplankton
     call para_get_value('THETA_KDISS_DET_PART_ORG_N'      ,       THETA_KDISS_DET_PART_ORG_N) !203 Particulate Detritus Nitrogen Dissolution rate Temperature correction
     call para_get_value('KHS_DISS_N'                      ,                       KHS_DISS_N) !204 Particulate Detritus Nitrogen dissolution reverse half saturation for DIN
     call para_get_value('FAC_PHYT_DET_PART_ORG_N'         ,          FAC_PHYT_DET_PART_ORG_N) !205 Particulate Detritus Nitrogen Phytoplankton linear factor for dissolution rate
     call para_get_value('KDISS_DET_PART_ORG_P_20'         ,          KDISS_DET_PART_ORG_P_20) !206 Particulate Detritus Phosphorus Dissolution rate not dependent on phytoplankton
     call para_get_value('THETA_KDISS_DET_PART_ORG_P'      ,       THETA_KDISS_DET_PART_ORG_P) !207 Particulate Detritus Phosphorus Dissolution rate Temperature correction
     call para_get_value('KHS_DISS_P'                      ,                       KHS_DISS_P) !208 Particulate Detritus Phosphorus  dissolution reverse half saturation for DIP
     call para_get_value('FAC_PHYT_DET_PART_ORG_P'         ,          FAC_PHYT_DET_PART_ORG_P) !209 Particulate Detritus Phosphorus  Phytoplankton linear factor for dissolution rate
     call para_get_value('KDISS_PART_Si_20'                ,                 KDISS_PART_Si_20) !210 Particulate Silica Dissolution rate
     call para_get_value('THETA_KDISS_PART_Si'             ,              THETA_KDISS_PART_Si) !211 Particulate Silica Dissolution rate Temperature correction
     call para_get_value('K_MIN_DOC_20'                    ,                     K_MIN_DOC_20) !212 Dissolved carbon  mineralisation rate
     call para_get_value('THETA_K_MIN_DOC'                 ,                  THETA_K_MIN_DOC) !213 Dissolved carbon  mineralisation rate Temperature constant
     call para_get_value('FAC_PHYT_AMIN_DOC'               ,                FAC_PHYT_AMIN_DOC) !214 Dissolved carbon  Phytoplankton linear factor for mineralisation rate
     call para_get_value('K_MIN_DON_20'                    ,                     K_MIN_DON_20) !215 Dissolved nitrogen  mineralisation rate not dependent on phytoplankton
     call para_get_value('THETA_K_MIN_DON'                 ,                  THETA_K_MIN_DON) !216 Dissolved nitrogen  mineralisation rate Temperature constant
     call para_get_value('KHS_AMIN_N'                      ,                       KHS_AMIN_N) !217 Dissolved nitrogen  reverse half saturation for DIN
     call para_get_value('FAC_PHYT_AMIN_DON'               ,                FAC_PHYT_AMIN_DON) !218 Dissolved nitrogen Phytoplankton linear factor for mineralisation rate
     call para_get_value('K_MIN_DOP_20'                    ,                     K_MIN_DOP_20) !219 Dissolved phosphorus  mineralisation rate not dependent on phytoplankton
     call para_get_value('THETA_K_MIN_DOP'                 ,                  THETA_K_MIN_DOP) !220 Dissolved phosphorus  mineralisation rate Temperature constant
     call para_get_value('KHS_AMIN_P'                      ,                       KHS_AMIN_P) !221 Dissolved phosphorus reverse half saturation for DIP
     call para_get_value('FAC_PHYT_AMIN_DOP'               ,                FAC_PHYT_AMIN_DOP) !222 Dissolved phosphorus Phytoplankton linear factor for mineralisation rate
     call para_get_value('K_NITR_20'                       ,                        K_NITR_20) !223 Amonia nitrification rate
     call para_get_value('KHS_NITR_OXY'                    ,                     KHS_NITR_OXY) !224 Amonia nitrification half saturation for Oxygen
     call para_get_value('KHS_NITR_NH4_N'                  ,                   KHS_NITR_NH4_N) !225 Amonia nitrification half saturation for Amonia
     call para_get_value('THETA_K_NITR'                    ,                     THETA_K_NITR) !226 Amonia nitrification rate Temperature constant                    
                     
                     
                     
                     
                     
                     
    ! Calling CO2SYS
    CO2SYS_NUM_SAMPLES = nkn ! number of nodes
    CO2SYS_ntps = nkn ! correction of bug: just CO2SYS_NUM_SAMPLES is not passed to co2sys

    allocate(CO2SYS_PAR1         (CO2SYS_NUM_SAMPLES), &
             CO2SYS_PAR2         (CO2SYS_NUM_SAMPLES), &
             CO2SYS_PAR1TYPE     (CO2SYS_NUM_SAMPLES), &
             CO2SYS_PAR2TYPE     (CO2SYS_NUM_SAMPLES), &
             CO2SYS_SALT         (CO2SYS_NUM_SAMPLES), &
             CO2SYS_TEMPIN       (CO2SYS_NUM_SAMPLES), &
             CO2SYS_TEMPOUT      (CO2SYS_NUM_SAMPLES), &
             CO2SYS_PRESIN       (CO2SYS_NUM_SAMPLES), &
             CO2SYS_PRESOUT      (CO2SYS_NUM_SAMPLES), &
             CO2SYS_SI           (CO2SYS_NUM_SAMPLES), &
             CO2SYS_PO4          (CO2SYS_NUM_SAMPLES), &
             CO2SYS_pHSCALEIN    (CO2SYS_NUM_SAMPLES), &
             CO2SYS_K1K2CONSTANTS(CO2SYS_NUM_SAMPLES), &
             CO2SYS_KSO4CONSTANTS(CO2SYS_NUM_SAMPLES))

    if (RUN_CO2SYS .eq. 1) then
        CO2SYS_PAR1         (1:nkn) = TOT_ALK(1:nkn) * 1.0D6
        CO2SYS_PAR2         (1:nkn) = INORG_C(1:nkn) * 1.0D6
        CO2SYS_PAR1TYPE     (1:nkn) = 1
        CO2SYS_PAR2TYPE     (1:nkn) = 2
        CO2SYS_SALT         (1:nkn) = SALT(1:nkn)
        CO2SYS_TEMPIN       (1:nkn) = TEMP(1:nkn)
        CO2SYS_TEMPOUT      (1:nkn) = 0.0D0  !Does not matter for this case
        CO2SYS_PRESIN       (1:nkn) = 0.0D0  !Does not matter for this case
        CO2SYS_PRESOUT      (1:nkn) = 0.0D0  !Does not matter for this case
        CO2SYS_SI           (1:nkn) = (Diss_SI(1:nkn) / 28.0855D0) * 1.0D3
        CO2SYS_PO4          (1:nkn) = (PO4_P  (1:nkn) / 30.9737D0) * 1.0D3
        CO2SYS_pHSCALEIN    (1:nkn) = 1
        CO2SYS_K1K2CONSTANTS(1:nkn) = 4
        CO2SYS_KSO4CONSTANTS(1:nkn) = 1

        call CO2SYS(CO2SYS_PAR1         , CO2SYS_PAR2  , CO2SYS_PAR1TYPE , &
                    CO2SYS_PAR2TYPE     , CO2SYS_SALT  , CO2SYS_TEMPIN   , &
                    CO2SYS_TEMPOUT      , CO2SYS_PRESIN, CO2SYS_PRESOUT  , &
                    CO2SYS_SI           , CO2SYS_PO4   , CO2SYS_pHSCALEIN, &
                    CO2SYS_K1K2CONSTANTS, CO2SYS_KSO4CONSTANTS, CO2SYS_OUT_DATA , &
                    CO2SYS_NICEHEADERS  , &
                    CO2SYS_ntps)

        pH         (1:nkn) = CO2SYS_OUT_DATA(1:nkn, 18)
        K_ONE_TIP  (1:nkn) = CO2SYS_OUT_DATA(1:nkn, 75)
        K_TWO_TIP  (1:nkn) = CO2SYS_OUT_DATA(1:nkn, 76)
        K_THREE_TIP(1:nkn) = CO2SYS_OUT_DATA(1:nkn, 77)
        ALPHA_0    (1:nkn) = CO2SYS_OUT_DATA(1:nkn, 23)
        H_PLUS     (1:nkn) = 10.0D0 ** (-CO2SYS_OUT_DATA(1:nkn,18))

        deallocate(CO2SYS_PAR1         )
        deallocate(CO2SYS_PAR2         )
        deallocate(CO2SYS_PAR1TYPE     )
        deallocate(CO2SYS_PAR2TYPE     )
        deallocate(CO2SYS_SALT         )
        deallocate(CO2SYS_TEMPIN       )
        deallocate(CO2SYS_TEMPOUT      )
        deallocate(CO2SYS_PRESIN       )
        deallocate(CO2SYS_PRESOUT      )
        deallocate(CO2SYS_SI           )
        deallocate(CO2SYS_PO4          )
        deallocate(CO2SYS_pHSCALEIN    )
        deallocate(CO2SYS_K1K2CONSTANTS)
        deallocate(CO2SYS_KSO4CONSTANTS)
        deallocate(CO2SYS_OUT_DATA     )
        deallocate(CO2SYS_NICEHEADERS  )
    end if ! call co2sys

    ! Calculate derived variables
    call derived_vars(nkn,pH,STATE_VARIABLES, nstate, &
                      MODEL_CONSTANTS, nconst,WC_OUTPUTS, noutput)

    !*****************************************
    !     D I S S O L V E D  O X Y G E N     !
    !*****************************************

    do k=1,nkn
        DISS_OXYGEN_SAT(k) = DO_SATURATION(TEMP(k), SALT(k), ELEVATION(k))

        if (SURFACE_BOX == 1) then

            if (K_A < 0.0D0) then
                K_A_CALC(k) = KAWIND(WINDS(k), TEMP(k), AIRTEMP(k), DEPTH(k), 3.0D0)
                R_AERATION(k) = K_A_CALC(k) * (DISS_OXYGEN_SAT(k) - DISS_OXYGEN(k))
            else
                K_A_CALC(k) = K_A
                R_AERATION(k) = K_A_CALC(k) * (DISS_OXYGEN_SAT(k) - DISS_OXYGEN(k)) * &
                &    (THETA_K_A ** (TEMP(k) - 2.0D1))
            end if

            !----------------------------------------------------------------------
            ! 2 February 2015 
            ! New code added to account the effect of ice cover.
            !----------------------------------------------------------------------
            R_AERATION(k) = (1.0D0 - ice_cover(k)) * R_AERATION(k)
            !----------------------------------------------------------------------
            ! End of new code added to account the effect of ice cover.
            !----------------------------------------------------------------------
        else
            R_AERATION(k) = 0.0D0
        end if
    end do

    ! Calculate the total phytoplankton.
    PHYT_TOT_C = DIA_C + CYN_C + OPA_C + FIX_CYN_C

    !*************!
    !*************!
    !  BACTERIA   !
    !*************!
    !*************!


    !***************************************!
    !     C H E M O A U T O T R O P H S     !
    !***************************************!

    call CHEMOAUTOTROPHS_1 &
           (KG_CHEM_AUT_BAC_20         , &
            THETA_KG_CHEM_AUT_BAC      , &
            KR_CHEM_AUT_BAC_20         , &
            THETA_KR_CHEM_AUT_BAC      , &
            KD_CHEM_AUT_BAC_20         , &
            THETA_KD_CHEM_AUT_BAC      , &
            KHS_NH4N_CHEM_AUT_BAC      , &
            KHS_PO4P_CHEM_AUT_BAC      , &
            KHS_O2_CHEM_AUT_BAC        , &
            DO_STR_HYPOX_CHEM_AUT_BAC_D, &
            THETA_HYPOX_CHEM_AUT_BAC_D , &
            EXPON_HYPOX_CHEM_AUT_BAC_D , &
            CHEM_AUT_BAC_N_TO_C        , &
            CHEM_AUT_BAC_P_TO_C        , &
            CHEM_AUT_BAC_O2_TO_C       , &
            YIELD_CHEM_AUT_BAC         , &
            EFF_CHEM_AUT_BAC_GROWTH    , &
            TIME_STEP                  , &
            nkn                        , &
            TEMP                       , &
            NH4_N                      , &
            NO3_N                      , &
            PO4_P                      , &
            DISS_OXYGEN                , &
            CHEM_AUT_BAC_C             , &
            LIM_TEMP_CHEM_AUT_BAC      , &
            LIM_NH4_N_CHEM_AUT_BAC     , &
            LIM_PO4_P_CHEM_AUT_BAC     , &
            LIM_OXY_CHEM_AUT_BAC       , &
            R_CHEM_AUT_BAC_GROWTH      , &
            R_CHEM_AUT_BAC_RESP        , &
            R_CHEM_AUT_BAC_INT_RESP    , &
            R_CHEM_AUT_BAC_DEATH       , &
            KD_CHEM_AUT_BAC            , &
            FAC_HYPOX_CHEM_AUT_BAC_D)

    !**********************************************************************!
    !     A E R O B I C and A N A E R O B I C  H E T E R O T R O P H S     !
    !**********************************************************************!

    call HETEROTROPHS_1 &
           (KG_AER_HET_BAC_20            , &
            EFF_AER_HET_BAC_GROWTH       , &
            THETA_KG_AER_HET_BAC         , &
            KR_AER_HET_BAC_20            , &
            THETA_KR_AER_HET_BAC         , &
            KD_AER_HET_BAC_20            , &
            THETA_KD_AER_HET_BAC         , &
            KHS_ORGC_AER_HET_BAC         , &
            KHS_ORGN_AER_HET_BAC         , &
            KHS_ORGP_AER_HET_BAC         , &
            KHS_O2_AER_HET_BAC           , &
            KHS_DIN_AER_HET_BAC          , &
            KHS_DIP_AER_HET_BAC          , &
            KHS_PHYT_AER_HET_BAC         , &
            YIELD_OC_AER_HET_BAC         , &
            OX_ORGN_AER_HET_BAC          , &
            KHS_MIN_N                    , &
            OX_ORGP_AER_HET_BAC          , &
            KHS_MIN_P                    , &
            DO_STR_HYPOX_AER_HET_BAC_D   , &
            THETA_HYPOX_AER_HET_BAC_D    , &
            EXPON_HYPOX_AER_HET_BAC_D    , &
            AER_HET_BAC_N_TO_C           , &
            AER_HET_BAC_P_TO_C           , &
            AER_HET_BAC_O2_TO_C          , &
            KG_FAC_AN_HET_BAC_20         , &
            KHS_NO3N_FAC_AN_HET_BAC      , &
            YIELD_FAC_AN_HET_BAC         , &
            TIME_STEP                    , &
            nkn                          , &
            TEMP                         , &
            DISS_OXYGEN                  , &
            DISS_ORG_C                   , &
            DISS_ORG_N                   , &
            DISS_ORG_P                   , &
            NH4_N                        , &
            NO3_N                        , &
            PO4_P                        , &
            AER_HET_BAC_C                , &
            PHYT_TOT_C                   , &
            LIM_TEMP_AER_HET_BAC         , &
            LIM_DISS_ORG_C_AER_HET_BAC   , &
            LIM_DISS_ORG_N_AER_HET_BAC   , &
            LIM_DISS_ORG_P_AER_HET_BAC   , &
            LIM_OXY_AER_HET_BAC          , &
            LIM_DIN_AER_HET_BAC          , &
            LIM_DIP_AER_HET_BAC          , &
            LIM_PHYT_C_AER_HET_BAC       , &
            R_AER_HET_BAC_GROWTH         , &
            R_AER_HET_BAC_INT_RESP       , &
            KD_AER_HET_BAC               , &
            FAC_HYPOX_AER_HET_BAC_D      , &
            LIM_DISS_NO3_N_FAC_AN_HET_BAC, &
            R_DENITRIFICATION            , &
            R_AER_HET_BAC_DEATH          , &
            R_AER_HET_BAC_RESP           , &
            PREF_NH4N_AER_HET_BAC)

    !*********************************************************!
    !     F A C U L T A T I V E   A N A E R O B S not used    !
    !*********************************************************!

    call HETEROTROPHS_2 &
           (KG_FAC_AN_HET_BAC_20           , &
            THETA_KG_FAC_AN_HET_BAC        , &
            KR_FAC_AN_HET_BAC_20           , &
            THETA_KR_FAC_AN_HET_BAC        , &
            KD_FAC_AN_HET_BAC_20           , &
            THETA_KD_FAC_AN_HET_BAC        , &
            KHS_NO3N_FAC_AN_HET_BAC        , &
            KHS_ORGC_FAC_AN_HET_BAC        , &
            KHS_ORGN_FAC_AN_HET_BAC        , &
            KHS_ORGP_FAC_AN_HET_BAC        , &
            REV_KHS_O2_FAC_AN_HET_BAC      , &
            NO3N_LACK_STR_FAC_AN_HET_BAC_D , &
            THETA_NO3_LACK_FAC_AN_HET_BAC_D, &
            EXP_NO3_LACK_FAC_AN_HET_BAC_D  , &
            FAC_AN_HET_BAC_N_TO_C          , &
            FAC_AN_HET_BAC_P_TO_C          , &
            FAC_AN_HET_BAC_O2_TO_C         , &
            YIELD_FAC_AN_HET_BAC           , &
            EFF_FAC_AN_HET_BAC_GROWTH      , &
            TIME_STEP                      , &
            nkn                            , &
            TEMP                           , &
            NO3_N                          , &
            DISS_ORG_C                     , &
            DISS_ORG_N                     , &
            DISS_ORG_P                     , &
            DISS_OXYGEN                    , &
            FAC_AN_HET_BAC_C               , &
            LIM_TEMP_FAC_AN_HET_BAC        , &
            LIM_DISS_NO3_N_FAC_AN_HET_BAC  , &
            LIM_DISS_ORG_C_FAC_AN_HET_BAC  , &
            LIM_DISS_ORG_N_FAC_AN_HET_BAC  , &
            LIM_DISS_ORG_P_FAC_AN_HET_BAC  , &
            LIM_OXY_FAC_AN_HET_BAC         , &
            R_FAC_AN_HET_BAC_GROWTH        , &
            R_FAC_AN_HET_BAC_RESP          , &
            R_FAC_AN_HET_BAC_INT_RESP      , &
            KD_FAC_AN_HET_BAC              , &
            R_FAC_AN_HET_BAC_DEATH)

    !**********************************!
    !**********************************!
    !     P H Y T O P L A N K T O N    !
    !**********************************!
    !**********************************!

    ! total chlorophyl in micrograms
    CHLA = ((DIA_C / DIA_C_TO_CHLA) + (CYN_C / CYN_C_TO_CHLA) + &
            (FIX_CYN_C / FIX_CYN_C_TO_CHLA) + (OPA_C / OPA_C_TO_CHLA)) * 1.0D3

    ! light extinction coeficient for smith = 0
    K_E = K_B_E + (8.8D-3 * CHLA) + (5.4D-2 * (CHLA ** (2.0D0 / 3.0D0)))


    !********************
    !      DIATOMS      !
    !********************

    !Calculations for diatom growth
    call DIATOMS(KG_DIA_OPT_TEMP         , &
                 DIA_OPT_TEMP_LR         , &
                 DIA_OPT_TEMP_UR         , &
                 EFF_DIA_GROWTH          , &
                 KAPPA_DIA_UNDER_OPT_TEMP, &
                 KAPPA_DIA_OVER_OPT_TEMP , &
                 KR_DIA_20               , &
                 THETA_KR_DIA            , &
                 KD_DIA_20               , &
                 THETA_KD_DIA            , &
                 KHS_DIN_DIA             , &
                 KHS_DIP_DIA             , &
                 KHS_DSi_DIA             , &
                 KHS_O2_DIA              , &
                 FRAC_DIA_EXCR           , &
                 I_S_DIA                 , &
                 DO_STR_HYPOX_DIA_D      , &
                 THETA_HYPOX_DIA_D       , &
                 EXPON_HYPOX_DIA_D       , &
                 DIA_N_TO_C              , &
                 DIA_P_TO_C              , &
                 DIA_Si_TO_C             , &
                 DIA_O2_TO_C             , &
                 DIA_C_TO_CHLA           , &
                 DIA_LIGHT_SAT           , &
                 NH4_N                   , &
                 NO3_N                   , &
                 PO4_P                   , &
                 DISS_OXYGEN             , &
                 DIA_C                   , &
                 ZOO_C                   , &
                 DISS_Si                 , &
                 TEMP                    , &
                 I_A                     , &
                 K_E                     , &
                 DEPTH                   , &
                 CHLA                    , &
                 K_B_E                   , &
                 FDAY                    , &
                 TIME_STEP               , &
                 SMITH                   , &
                 nkn                     , &
                 KG_DIA                  , &
                 ALPHA_0                 , &
                 ALPHA_1                 , &
                 LIM_KG_DIA_TEMP         , &
                 LIM_KG_DIA_LIGHT        , &
                 LIM_KG_DIA_DOXY         , &
                 LIM_KG_DIA_N            , &
                 LIM_KG_DIA_P            , &
                 LIM_KG_DIA_DISS_Si      , &
                 LIM_KG_DIA_NUTR         , &
                 LIM_KG_DIA              , &
                 R_DIA_GROWTH            , &
                 R_DIA_MET               , &
                 R_DIA_RESP              , &
                 R_DIA_EXCR              , &
                 R_DIA_INT_RESP          , &
                 KD_DIA                  , &
                 FAC_HYPOX_DIA_D         , &
                 R_DIA_DEATH             , &
                 PREF_NH4N_DIA)

    !**************************************
    ! NON-NITROGEN FIXING CYANOBACTERIA   !
    !*************************************

    !Calculations for non-fixing cyanobacteria growth
    call CYANOBACTERIA &
           (KG_CYN_OPT_TEMP         , &
            CYN_OPT_TEMP_LR         , &
            CYN_OPT_TEMP_UR         , &
            EFF_CYN_GROWTH          , &
            KAPPA_CYN_UNDER_OPT_TEMP, &
            KAPPA_CYN_OVER_OPT_TEMP , &
            KR_CYN_20               , &
            THETA_KR_CYN            , &
            KD_CYN_20               , &
            THETA_KD_CYN            , &
            KHS_DIN_CYN             , &
            KHS_DIP_CYN             , &
            KHS_O2_CYN              , &
            FRAC_CYN_EXCR           , &
            I_S_CYN                 , &
            DO_STR_HYPOX_CYN_D      , &
            THETA_HYPOX_CYN_D       , &
            EXPON_HYPOX_CYN_D       , &
            CYN_N_TO_C              , &
            CYN_P_TO_C              , &
            CYN_O2_TO_C             , &
            CYN_C_TO_CHLA           , &
            CYN_LIGHT_SAT           , &
            NH4_N                   , &
            NO3_N                   , &
            PO4_P                   , &
            DISS_OXYGEN             , &
            CYN_C                   , &
            ZOO_C                   , &
            TEMP                    , &
            I_A                     , &
            K_E                     , &
            DEPTH                   , &
            CHLA                    , &
            K_B_E                   , &
            FDAY                    , &
            TIME_STEP               , &
            SMITH                   , &
            nkn                     , &
            KG_CYN                  , &
            ALPHA_0                 , &
            ALPHA_1                 , &
            LIM_KG_CYN_TEMP         , &
            LIM_KG_CYN_LIGHT        , &
            LIM_KG_CYN_DOXY         , &
            LIM_KG_CYN_N            , &
            LIM_KG_CYN_P            , &
            LIM_KG_CYN_NUTR         , &
            LIM_KG_CYN              , &
            R_CYN_GROWTH            , &
            R_CYN_MET               , &
            R_CYN_RESP              , &
            R_CYN_EXCR              , &
            R_CYN_INT_RESP          , &
            KD_CYN                  , &
            FAC_HYPOX_CYN_D         , &
            R_CYN_DEATH             , &
            PREF_NH4N_CYN)

    !********************************
    ! NITROGEN FIXING CYANOBACTERIA !
    !********************************

    call FIX_CYANOBACTERIA  &
           (KG_FIX_CYN_OPT_TEMP          , &
            FIX_CYN_OPT_TEMP_LR          , &
            FIX_CYN_OPT_TEMP_UR          , &
            EFF_FIX_CYN_GROWTH           , &
            KAPPA_FIX_CYN_UNDER_OPT_TEMP , &
            KAPPA_FIX_CYN_OVER_OPT_TEMP  , &
            KR_FIX_CYN_20                , &
            THETA_KR_FIX_CYN             , &
            KD_FIX_CYN_20                , &
            THETA_KD_FIX_CYN             , &
            KHS_DIN_FIX_CYN              , &
            KHS_DIP_FIX_CYN              , &
            KHS_O2_FIX_CYN               , &
            KHS_NH4N_PREF_FIX_CYN        , &
            I_S_FIX_CYN                  , &
            DO_STR_HYPOX_FIX_CYN_D       , &
            THETA_HYPOX_FIX_CYN_D        , &
            EXPON_HYPOX_FIX_CYN_D        , &
            FIX_CYN_N_TO_C               , &
            FIX_CYN_P_TO_C               , &
            FIX_CYN_O2_TO_C              , &
            FIX_CYN_C_TO_CHLA            , &
            FIX_CYN_C_TO_CHLA_NEW        , &
            FIX_CYN_LIGHT_SAT            , &
            FRAC_FIX_CYN_EXCR            , &
            R_FIX                        , &
            K_FIX                        , &
            TIME_STEP                    , &
            SMITH                        , &
            nkn                          , &
            FDAY                         , &
            I_A                          , &
            K_E                          , &
            K_B_E                        , &
            DEPTH                        , &
            CHLA                         , &
            TEMP                         , &
            NH4_N                        , &
            NO3_N                        , &
            PO4_P                        , &
            DISS_OXYGEN                  , &
            FIX_CYN_C                    , &
            ALPHA_0                      , &
            ALPHA_1                      , &
            KG_FIX_CYN                   , &
            LIM_KG_FIX_CYN_LIGHT         , &
            LIM_KG_FIX_CYN_TEMP          , &
            LIM_KG_FIX_CYN_DOXY          , &
            LIM_KG_NON_FIX_CYN_N         , &
            LIM_KG_NON_FIX_CYN_P         , &
            LIM_KG_NON_FIX_CYN_NUTR      , &
            LIM_KG_FIX_FIX_CYN_N         , &
            LIM_KG_FIX_FIX_CYN_P         , &
            LIM_KG_FIX_FIX_CYN_NUTR      , &
            LIM_KG_NON_FIX_CYN           , &
            LIM_KG_FIX_FIX_CYN           , &
            R_NON_FIX_CYN_GROWTH         , &
            R_FIX_FIX_CYN_GROWTH         , &
            R_FIX_CYN_GROWTH             , &
            R_FIX_CYN_MET                , &
            R_FIX_CYN_RESP               , &
            R_FIX_CYN_EXCR               , &
            R_FIX_CYN_INT_RESP           , &
            KD_FIX_CYN                   , &
            FAC_HYPOX_FIX_CYN_D          , &
            R_FIX_CYN_DEATH              , &
            PREF_NH4N_FIX_CYN)
    
    !*****************************
    ! OTHER PLANKTONIC ALGAE     !
    !*****************************
    call OTHER_PLANKTONIC_ALGAE &
           (KG_OPA_OPT_TEMP         , &
            OPA_OPT_TEMP_LR         , &
            OPA_OPT_TEMP_UR         , &
            EFF_OPA_GROWTH          , &
            KAPPA_OPA_UNDER_OPT_TEMP, &
            KAPPA_OPA_OVER_OPT_TEMP , &
            KR_OPA_20               , &
            THETA_KR_OPA            , &
            KD_OPA_20               , &
            THETA_KD_OPA            , &
            KHS_DIN_OPA             , &
            KHS_DIP_OPA             , &
            KHS_O2_OPA              , &
            FRAC_OPA_EXCR           , &
            I_S_OPA                 , &
            DO_STR_HYPOX_OPA_D      , &
            THETA_HYPOX_OPA_D       , &
            EXPON_HYPOX_OPA_D       , &
            OPA_N_TO_C              , &
            OPA_P_TO_C              , &
            OPA_O2_TO_C             , &
            OPA_C_TO_CHLA           , &
            OPA_LIGHT_SAT           , &
            NH4_N                   , &
            NO3_N                   , &
            PO4_P                   , &
            DISS_OXYGEN             , &
            OPA_C                   , &
            ZOO_C                   , &
            TEMP                    , &
            I_A                     , &
            K_E                     , &
            DEPTH                   , &
            CHLA                    , &
            K_B_E                   , &
            FDAY                    , &
            TIME_STEP               , &
            SMITH                   , &
            nkn                     , &
            KG_OPA                  , &
            ALPHA_0                 , &
            ALPHA_1                 , &
            LIM_KG_OPA_TEMP         , &
            LIM_KG_OPA_LIGHT        , &
            LIM_KG_OPA_DOXY         , &
            LIM_KG_OPA_N            , &
            LIM_KG_OPA_P            , &
            LIM_KG_OPA_NUTR         , &
            LIM_KG_OPA              , &
            R_OPA_GROWTH            , &
            R_OPA_MET               , &
            R_OPA_RESP              , &
            R_OPA_EXCR              , &
            R_OPA_INT_RESP          , &
            KD_OPA                  , &
            FAC_HYPOX_OPA_D         , &
            R_OPA_DEATH             , &
            PREF_NH4N_OPA)

    !******************************!
    !******************************!
    !     Z O O P L A N K T O N    !
    !******************************!
    !******************************!
    call ZOOPLANKTON &
           (KG_ZOO_OPT_TEMP               , &
            ZOO_OPT_TEMP_LR               , &
            ZOO_OPT_TEMP_UR               , &
            EFF_ZOO_GROWTH                , &
            KAPPA_ZOO_UNDER_OPT_TEMP      , &
            KAPPA_ZOO_OVER_OPT_TEMP       , &
            GRAT_ZOO_DIA                  , &
            GRAT_ZOO_CYN                  , &
            GRAT_ZOO_OPA                  , &
            GRAT_ZOO_FIX_CYN              , &
            GRAT_ZOO_CHEM_AUT_BAC         , &
            GRAT_ZOO_AER_HET_BAC          , &
            GRAT_ZOO_FAC_AN_HET_BAC       , &
            GRAT_ZOO_DET_PART_ORG_C       , &
            PREF_ZOO_DIA                  , &
            PREF_ZOO_CYN                  , &
            PREF_ZOO_FIX_CYN              , &
            PREF_ZOO_OPA                  , &
            PREF_ZOO_CHEM_AUT_BAC         , &
            PREF_ZOO_AER_HET_BAC          , &
            PREF_ZOO_FAC_AN_HET_BAC       , &
            PREF_ZOO_DET_PART_ORG_C       , &
            KHS_DIA_C_ZOO                 , &
            KHS_CYN_C_ZOO                 , &
            KHS_FIX_CYN_C_ZOO             , &
            KHS_OPA_C_ZOO                 , &
            KHS_CHEM_AUT_BAC_C_ZOO        , &
            KHS_AER_HET_BAC_C_ZOO         , &
            KHS_FAC_AN_HET_BAC_C_ZOO      , &
            KHS_DET_PART_ORG_C_ZOO        , &
            FOOD_MIN_ZOO                  , &
            KE_ZOO                        , &
            FRAC_ZOO_EX_ORG               , &
            KR_ZOO_20                     , &
            THETA_KR_ZOO                  , &
            KD_ZOO_20                     , &
            THETA_KD_ZOO                  , &
            DO_STR_HYPOX_ZOO_D            , &
            THETA_HYPOX_ZOO_D             , &
            EXPON_HYPOX_ZOO_D             , &
            ZOO_N_TO_C                    , &
            ZOO_P_TO_C                    , &
            ZOO_O2_TO_C                   , &
            TEMP                          , &
            DISS_OXYGEN                   , &
            DIA_C                         , &
            CYN_C                         , &
            OPA_C                         , &
            FIX_CYN_C                     , &
            CHEM_AUT_BAC_C                , &
            AER_HET_BAC_C                 , &
            FAC_AN_HET_BAC_C              , &
            DET_PART_ORG_C                , &
            ZOO_C                         , &
            TIME_STEP                     , &
            nkn                           , &
            KG_ZOO                        , &
            KG_ZOO_DIA                    , &
            KG_ZOO_CYN                    , &
            KG_ZOO_OPA                    , &
            KG_ZOO_FIX_CYN                , &
            KG_ZOO_CHEM_AUT_BAC           , &
            KG_ZOO_AER_HET_BAC            , &
            KG_ZOO_FAC_AN_HET_BAC         , &
            KG_ZOO_DET_PART_ORG_C         , &
            KD_ZOO                        , &
            FOOD_FACTOR_ZOO_DIA           , &
            FOOD_FACTOR_ZOO_CYN           , &
            FOOD_FACTOR_ZOO_OPA           , &
            FOOD_FACTOR_ZOO_FIX_CYN       , &
            FOOD_FACTOR_ZOO_CHEM_AUT_BAC  , &
            FOOD_FACTOR_ZOO_AER_HET_BAC   , &
            FOOD_FACTOR_ZOO_FAC_AN_HET_BAC, &
            FOOD_FACTOR_ZOO_DET_PART_ORG_C, &
            R_ZOO_FEEDING_DIA             , &
            R_ZOO_FEEDING_CYN             , &
            R_ZOO_FEEDING_FIX_CYN         , &
            R_ZOO_FEEDING_OPA             , &
            R_ZOO_FEEDING_CHEM_AUT_BAC    , &
            R_ZOO_FEEDING_AER_HET_BAC     , &
            R_ZOO_FEEDING_FAC_AN_HET_BAC  , &
            R_ZOO_FEEDING_DET_PART_ORG_C  , &
            R_ZOO_INT_RESP                , &
            R_ZOO_RESP                    , &
            R_ZOO_EX_DON                  , &
            R_ZOO_EX_DOP                  , &
            R_ZOO_EX_DOC                  , &
            R_ZOO_DEATH                   , &
            ACTUAL_ZOO_N_TO_C             , &
            ACTUAL_ZOO_P_TO_C             , &
            R_ZOO_GROWTH                  , &
            FAC_HYPOX_ZOO_D)

    if(debug_stranger) then
        if (STRANGERSD(R_ZOO_GROWTH,VALUE_strange,nkn).eq.1) then
            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'

            nstrange = count(VALUE_strange)
            allocate(STRANGERS    (nstrange))
            allocate(NODES_STRANGE(nstrange))

            j=1

            do k=1,nkn
                if(VALUE_strange(k)) then
                    STRANGERS    (j) = R_ZOO_GROWTH(k)
                    NODES_STRANGE(j) = node_active(k)
                    j=j+1
                end if
            end do

            print *, 'PELAGIC KINETICS: '
            write(*,*) 'TIME          : ', TIME
            write(*,*) 'R_ZOO_GROWTH is NaN or Inf:'
            print *, 'NODE_NUMBERS=',NODES_STRANGE
            print *, 'VALUES=',STRANGERS

            deallocate(STRANGERS)
            deallocate(NODES_STRANGE)

            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'
            write(*,*) 'R_ZOO_FEEDING_DIA                : ', R_ZOO_FEEDING_DIA
            write(*,*) 'R_ZOO_FEEDING_CYN                : ', R_ZOO_FEEDING_CYN
            write(*,*) 'R_ZOO_FEEDING_OPA                : ', R_ZOO_FEEDING_OPA
            write(*,*) 'R_ZOO_FEEDING_FIX_CYN            : ', R_ZOO_FEEDING_CYN
            write(*,*) 'R_ZOO_FEEDING_AER_HET_BAC        : ', R_ZOO_FEEDING_AER_HET_BAC
            write(*,*) 'R_ZOO_FEEDING_FAC_AN_HET_BAC     : ', R_ZOO_FEEDING_FAC_AN_HET_BAC
            write(*,*) 'R_ZOO_FEEDING_DET_PART_ORG_C     : ', R_ZOO_FEEDING_DET_PART_ORG_C
            write(*,*) '    ZOO_C                        : ', ZOO_C
            write(*,*) '    KG_ZOO_DIA                   : ', KG_ZOO_DIA
            write(*,*) '    KG_ZOO_CYN                   : ', KG_ZOO_CYN
            write(*,*) '    KG_ZOO_OPA                   : ', KG_ZOO_OPA
            write(*,*) '    FOOD_FACTOR_ZOO_OPA          : ', FOOD_FACTOR_ZOO_OPA
            write(*,*) '    KG_ZOO_CHEM_AUT_BAC          : ', KG_ZOO_CHEM_AUT_BAC
            write(*,*) '    FOOD_FACTOR_ZOO_CHEM_AUT_BAC : ', FOOD_FACTOR_ZOO_CHEM_AUT_BAC
            write(*,*) '    KG_ZOO_AER_HET_BAC           : ', KG_ZOO_AER_HET_BAC
            write(*,*) '    FOOD_FACTOR_ZOO_AER_HET_BAC  : ', FOOD_FACTOR_ZOO_AER_HET_BAC
            write(*,*) '    KG_ZOO_DET_PART_ORG_C        : ', KG_ZOO_DET_PART_ORG_C
            write(*,*) '    KG_ZOO_DET_PART_ORG_C        : ', KG_ZOO_DET_PART_ORG_C
            write(*,*) '    KG_ZOO                       : ', KG_ZOO
            stop
        end if
    end if

    if(debug_stranger) then
        if (STRANGERSD(R_ZOO_RESP,VALUE_strange,nkn).eq.1) then
            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            print *, 'PELAGIC_KINETICS:'
            write(*,*) 'TIME   : ', TIME


            nstrange = count(VALUE_strange)
            allocate(STRANGERS    (nstrange))
            allocate(NODES_STRANGE(nstrange))

            j=1
            do k=1,nkn
                if(VALUE_strange(k)) then
                    STRANGERS    (j) = R_ZOO_RESP(k)
                    NODES_STRANGE(j) = node_active(k)
                    j=j+1
                end if
            end do

            write(*,*) 'R_ZOO_RESP is NaN or Inf:'
            print *, 'Initial state is NaN or Inf:'
            print *, 'NODE_NUMBERS=',NODES_STRANGE
            print *, 'VALUES=',STRANGERS

            deallocate(STRANGERS)
            deallocate(NODES_STRANGE)

            write(*,*) 'R_ZOO_RESP is not a number or infinite :', R_ZOO_RESP
            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'
            write(*,*) 'R_ZOO_GROWTH      : ', R_ZOO_GROWTH
            write(*,*) 'EFF_ZOO_GROWTH    : ', EFF_ZOO_GROWTH
            stop
        end if
    end if


    !************************
    !NITRIFICATION           !
    !************************
    R_NITRIFICATION   = (1.0D0 / YIELD_CHEM_AUT_BAC)   * R_CHEM_AUT_BAC_GROWTH

    !moved to merged aerobic-anaerobic heterototrophs
    !R_DENITRIFICATION = (1.0D0 / YIELD_FAC_AN_HET_BAC) * R_FAC_AN_HET_BAC_GROWTH


    !*********************************************************************!
    !     D E A T H   O R G A N I C     P A R T I C L E S    DISSOLUTION !
    !*********************************************************************!

    call ORGANIC_CARBON_DISSOLUTION &
           (FAC_PHYT_DET_PART_ORG_C     , &
            KDISS_DET_PART_ORG_C_20     , &
            THETA_KDISS_DET_PART_ORG_C  , &
            nkn                         , &
            TEMP                        , &
            DET_PART_ORG_C              , &
            PHYT_TOT_C                  , &
            LIM_PHYT_DISS_DET_PART_ORG_C, &
            R_DET_PART_ORG_C_DISSOLUTION)

    ACTUAL_DET_N_TO_C = DET_PART_ORG_N / DET_PART_ORG_C
    ACTUAL_DET_P_TO_C = DET_PART_ORG_P / DET_PART_ORG_C

    ! Nitrogen dissolution

    ! Accerelation of hydrolysis when DIN is scarce
    LIM_N_DISS_DET_PART_ORG_N = KHS_DISS_N / (KHS_DISS_N + (NH4_N + NO3_N))
    LIM_PHY_N_DISS_DET_PART_ORG_N = LIM_N_DISS_DET_PART_ORG_N * FAC_PHYT_DET_PART_ORG_N * PHYT_TOT_C

    R_DET_PART_ORG_N_DISSOLUTION = (KDISS_DET_PART_ORG_N_20 + LIM_PHY_N_DISS_DET_PART_ORG_N) * &
       (THETA_KDISS_DET_PART_ORG_N ** (TEMP - 2.0D1)) * DET_PART_ORG_N

    ! Phosphorus dissolution

    ! Accerelation of hydrolysis when DIP is scarce
    LIM_P_DISS_DET_PART_ORG_P = KHS_DISS_P / (KHS_DISS_P + PO4_P)
    LIM_PHY_P_DISS_DET_PART_ORG_P = LIM_P_DISS_DET_PART_ORG_P * FAC_PHYT_DET_PART_ORG_P * PHYT_TOT_C

    R_DET_PART_ORG_P_DISSOLUTION = (KDISS_DET_PART_ORG_P_20 + LIM_PHY_P_DISS_DET_PART_ORG_P) * &
       (THETA_KDISS_DET_PART_ORG_P ** (TEMP - 2.0D1)) * DET_PART_ORG_P

   !********************************************
   ! RESPIRATION AND EXCRETION ( fixme ) RATES  !
   !********************************************

    !Chemoautotrophic bacteria total respiration rate
    R_CHEM_AUT_BAC_TOT_RESP      = R_CHEM_AUT_BAC_INT_RESP   + R_CHEM_AUT_BAC_RESP

    !Aerobic heterotrophic bacteria total respiration rate
    R_AER_HET_BAC_TOT_RESP       = R_AER_HET_BAC_INT_RESP    + R_AER_HET_BAC_RESP

    !Facultative anaerobic heterotrophic bacteria total respiration rate
    R_FAC_AN_HET_BAC_TOT_RESP    = R_FAC_AN_HET_BAC_INT_RESP + R_FAC_AN_HET_BAC_RESP

    !Diatom total respiration rate
    R_DIA_TOT_RESP               = R_DIA_RESP                + R_DIA_INT_RESP

    !Non-fixing cyanobacteria total respiration rate
    R_CYN_TOT_RESP               = R_CYN_RESP                + R_CYN_INT_RESP

    !Other planktonic algae total respiration rate
    R_OPA_TOT_RESP               = R_OPA_RESP                + R_OPA_INT_RESP

    !Nitrogen fixing cyanobacteria total respiration rate
    R_FIX_CYN_TOT_RESP           = R_FIX_CYN_RESP            + R_FIX_CYN_INT_RESP

    !Zooplankton total respiration rate
    R_ZOO_TOT_RESP               = R_ZOO_INT_RESP            + R_ZOO_RESP

    !***************************************************
    ! DISSOLVED ORGANIC N P MINERALIZATION BY BACTERIA !
    !**********************************8**** ***********
    !Rates for recycle of nutrients by aerobic heterotrophic bacteria during
    !organic matter oxidation.

    ! Accerelation of mineralisation when NIP is scarce
    LIM_N_MIN_DON_N = KHS_MIN_N / (KHS_MIN_N + (NH4_N + NO3_N))

    R_AER_HET_BAC_N_OX    = R_AER_HET_BAC_GROWTH * OX_ORGN_AER_HET_BAC * (1.D0 + LIM_N_MIN_DON_N)

    ! Accerelation of mineralisation when DIP is scarce
    LIM_P_MIN_DOP_P = KHS_MIN_P / (KHS_MIN_P + PO4_P)
    R_AER_HET_BAC_P_OX    = R_AER_HET_BAC_GROWTH * OX_ORGP_AER_HET_BAC * (1.D0 + LIM_P_MIN_DOP_P)

    !Rates for recycle of nutrients by facultative anaerobic heterotrophic bacteria during
    !organic matter oxidation. Not used more!
    R_FAC_AN_HET_BAC_N_OX = R_FAC_AN_HET_BAC_GROWTH * (DISS_ORG_N / DISS_ORG_C)
    R_FAC_AN_HET_BAC_P_OX = R_FAC_AN_HET_BAC_GROWTH * (DISS_ORG_P / DISS_ORG_C)


    !*********************************************************************!
    !     SILICON                                                         !
    !*********************************************************************!

    !Dissolution rate of biogenic silicon
    R_PART_Si_DISS = KDISS_PART_SI_20 * &
    &       (THETA_KDISS_PART_SI ** (TEMP - 2.0D1)) * PART_SI


    !*********************************************************************!
    !     MINERALIZATION OF DOC, DON, DOP when bacteria are not modelled.
    !     Called abiotic in the sense of modelling method
    !*********************************************************************!

    !Algal dependent mineralisation rate
    call ORGANIC_CARBON_MINERALIZATION &
           (FAC_PHYT_AMIN_DOC           , &
            K_MIN_DOC_20                , &
            THETA_K_MIN_DOC             , &
            nkn                         , &
            TEMP                        , &
            DISS_ORG_C                  , &
            PHYT_TOT_C                  , &
            LIM_PHYT_AMIN_DOC           , &
            R_ABIOTIC_DOC_MIN)

    ! Accerelation of mineralisation when DIN is scarce
    LIM_N_AMIN_DON = KHS_AMIN_N / (KHS_AMIN_N + (NH4_N + NO3_N))
    LIM_PHY_N_AMIN_DON = LIM_N_AMIN_DON * FAC_PHYT_AMIN_DON * PHYT_TOT_C

    R_ABIOTIC_DON_MIN = (K_MIN_DON_20 + LIM_PHY_N_AMIN_DON)* &
                        (THETA_K_MIN_DON ** (TEMP - 2.0D1)) * DISS_ORG_N

    ! Accerelation of mineralisation when DIP is scarce
    LIM_P_AMIN_DOP = KHS_AMIN_P / (KHS_AMIN_P + PO4_P)
    LIM_PHY_P_AMIN_DOP = LIM_P_AMIN_DOP * FAC_PHYT_AMIN_DOP * PHYT_TOT_C

    R_ABIOTIC_DOP_MIN = (K_MIN_DOP_20 + LIM_PHY_P_AMIN_DOP) * &
                        (THETA_K_MIN_DOP ** (TEMP - 2.0D1)) * DISS_ORG_P

    !*******************************************************************************************!
    !     Nitrification of ammonia when bacteria are not modelled.
    !     Called abiotic in the sense of modelling method
    !*******************************************************************************************!
     LIM_NITR_OXY = DISS_OXYGEN / (KHS_NITR_OXY + DISS_OXYGEN)
     LIM_NITR_NH4_N = NH4_N / (KHS_NITR_NH4_N + NH4_N)

     R_ABIOTIC_NITR = K_NITR_20 * LIM_NITR_OXY * LIM_NITR_NH4_N * &
                      (THETA_K_NITR ** (TEMP - 2.0D1)) * NH4_N


    !*********************************************************************!
    !     VOLATILIZATION OF UNIONIZED AMMONI.
    !*********************************************************************!
    call AMMONIA_VOLATILIZATION(R_AMMONIA_VOLATIL, NH4_N, pH, TEMP, K_A_CALC,nkn)

    !----------------------------------------------------------------------
    ! 2 February 2015 
    ! New code added to account the effect of ice cover.
    !----------------------------------------------------------------------
    R_AMMONIA_VOLATIL = R_AMMONIA_VOLATIL * (1.0D0 - ice_cover)
    !----------------------------------------------------------------------
    ! End of new code added to account the effect of ice cover.
    !----------------------------------------------------------------------

    !------------------------------------------------------------------------------------------------
    !      Final calculation of derivatives
    !------------------------------------------------------------------------------------------------


    !AMMONIA NITROGEN
    PROCESS_RATES(1:nkn,1, 1)  = R_CHEM_AUT_BAC_TOT_RESP   * CHEM_AUT_BAC_N_TO_C
    PROCESS_RATES(1:nkn,1, 2)  = R_AER_HET_BAC_INT_RESP    * AER_HET_BAC_N_TO_C    !?
    PROCESS_RATES(1:nkn,1, 3)  = R_FAC_AN_HET_BAC_TOT_RESP * FAC_AN_HET_BAC_N_TO_C
    PROCESS_RATES(1:nkn,1, 4)  = R_DIA_TOT_RESP            * DIA_N_TO_C
    PROCESS_RATES(1:nkn,1, 5)  = R_CYN_TOT_RESP            * CYN_N_TO_C
    PROCESS_RATES(1:nkn,1, 6)  = R_OPA_TOT_RESP            * OPA_N_TO_C
    PROCESS_RATES(1:nkn,1, 7)  = R_FIX_CYN_TOT_RESP        * FIX_CYN_N_TO_C
    PROCESS_RATES(1:nkn,1, 8)  = R_AER_HET_BAC_N_OX
    PROCESS_RATES(1:nkn,1, 9)  = R_FAC_AN_HET_BAC_N_OX
    PROCESS_RATES(1:nkn,1, 10) = R_ZOO_TOT_RESP         * ACTUAL_ZOO_N_TO_C
    PROCESS_RATES(1:nkn,1, 11) = R_CHEM_AUT_BAC_GROWTH  * CHEM_AUT_BAC_N_TO_C
    PROCESS_RATES(1:nkn,1, 12) = R_DIA_GROWTH         * PREF_NH4N_DIA     * DIA_N_TO_C
    PROCESS_RATES(1:nkn,1, 13) = R_CYN_GROWTH         * PREF_NH4N_CYN     * CYN_N_TO_C
    PROCESS_RATES(1:nkn,1, 14) = R_OPA_GROWTH         * PREF_NH4N_OPA     * OPA_N_TO_C
    PROCESS_RATES(1:nkn,1, 15) = R_NON_FIX_CYN_GROWTH * PREF_NH4N_FIX_CYN * FIX_CYN_N_TO_C
    PROCESS_RATES(1:nkn,1, 16) = R_NITRIFICATION + R_ABIOTIC_NITR
    PROCESS_RATES(1:nkn,1, 17) = R_ABIOTIC_DON_MIN
    PROCESS_RATES(1:nkn,1, 18) = R_AER_HET_BAC_GROWTH * PREF_NH4N_AER_HET_BAC * AER_HET_BAC_N_TO_C !uptake by aerobic heterotrophic bacteria
    ! Auxiliary
    PROCESS_RATES(1:nkn,1, 19) = PREF_NH4N_DIA
    PROCESS_RATES(1:nkn,1, 20) = PREF_NH4N_CYN
    PROCESS_RATES(1:nkn,1, 21) = PREF_NH4N_OPA
    PROCESS_RATES(1:nkn,1, 22) = PREF_NH4N_FIX_CYN
    PROCESS_RATES(1:nkn,1, 23) = (-1.0D0) * R_AMMONIA_VOLATIL


    DERIVATIVES(1:nkn,1) = PROCESS_RATES(1:nkn,1, 1)  + PROCESS_RATES(1:nkn,1, 2)  + PROCESS_RATES(1:nkn,1, 3)  + &
                     PROCESS_RATES(1:nkn,1, 4)  + PROCESS_RATES(1:nkn,1, 5)  + PROCESS_RATES(1:nkn,1, 6)  + &
                     PROCESS_RATES(1:nkn,1, 7)  + PROCESS_RATES(1:nkn,1, 8)  + PROCESS_RATES(1:nkn,1, 9)  + &
                     PROCESS_RATES(1:nkn,1, 10) - PROCESS_RATES(1:nkn,1, 11) - PROCESS_RATES(1:nkn,1, 12) - &
                     PROCESS_RATES(1:nkn,1, 13) - PROCESS_RATES(1:nkn,1, 14) - PROCESS_RATES(1:nkn,1, 15) - &
                     PROCESS_RATES(1:nkn,1, 16) + PROCESS_RATES(1:nkn,1, 17) - PROCESS_RATES(1:nkn,1, 18) - &
                     PROCESS_RATES(1:nkn,1, 23)

    !i=0
    !i=STRANGERSD(DERIVATIVES(1:nkn,1),VALUE_strange,nkn)

    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,1),VALUE_strange,nkn).eq.1) then
            nstrange = count(VALUE_strange)
            allocate(STRANGERS    (nstrange))
            allocate(NODES_STRANGE(nstrange))
            allocate(NODES_STRANGE_int(nstrange))
            allocate(NODES_STRANGE_ext(nstrange))

            j=1
            do k=1,nkn
                if(VALUE_strange(k)) then
                    STRANGERS    (j) = DERIVATIVES(k,1)
                    NODES_STRANGE(j) = node_active(k)
                    NODES_STRANGE_int(j) = node_active(k)
                    NODES_STRANGE_ext(j) = ipv(node_active(k))
                    j=j+1
                end if
            end do

            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            print *, 'PELAGIC_KINETICS:'
            write(*,*) 'TIME   : ', TIME

            write(*,*) 'DERIVATIVES(1) is not a number or infinite:'
            print *, 'NODE_NUMBERS=',NODES_STRANGE
            print *, 'NODE_NUMBERS int.=',NODES_STRANGE_int
            print *, 'NODE_NUMBERS ext.=',NODES_STRANGE_ext
            print *, 'VALUES=',STRANGERS

            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'
            write(*,*) 'NH4_N                        : ',  (NH4_N                    &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'R_CHEM_AUT_BAC_TOT_RESP       : ', (R_CHEM_AUT_BAC_TOT_RESP  &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    R_CHEM_AUT_BAC_INT_RESP   : ', (R_CHEM_AUT_BAC_INT_RESP  &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    R_CHEM_AUT_BAC_RESP       : ', (R_CHEM_AUT_BAC_RESP      &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    CHEM_AUT_BAC_N_TO_C       : ',  CHEM_AUT_BAC_N_TO_C

            write(*,*) ''
            write(*,*) 'R_AER_HET_BAC_TOT_RESP        : ', (R_AER_HET_BAC_TOT_RESP    &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    R_AER_HET_BAC_INT_RESP    : ', (R_AER_HET_BAC_INT_RESP    &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    R_AER_HET_BAC_RESP        : ', (R_AER_HET_BAC_RESP        &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    AER_HET_BAC_N_TO_C        : ',  AER_HET_BAC_N_TO_C

            write(*,*) ''
            write(*,*) 'R_FAC_AN_HET_BAC_TOT_RESP     : ', (R_FAC_AN_HET_BAC_TOT_RESP  &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    R_FAC_AN_HET_BAC_INT_RESP : ', (R_FAC_AN_HET_BAC_INT_RESP  &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    R_FAC_AN_HET_BAC_RESP     : ', (R_FAC_AN_HET_BAC_RESP      &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    FAC_AN_HET_BAC_N_TO_C     : ',  FAC_AN_HET_BAC_N_TO_C

            write(*,*) ''
            write(*,*) 'R_DIA_TOT_RESP                : ', (R_DIA_TOT_RESP  &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    R_DIA_INT_RESP            : ', (R_DIA_INT_RESP  &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    R_DIA_RESP                : ', (R_DIA_RESP      &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    DIA_N_TO_C                : ',  DIA_N_TO_C

            write(*,*) ''
            write(*,*) 'R_ZOO_TOT_RESP                : ', (R_ZOO_TOT_RESP  &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    R_ZOO_INT_RESP            : ', (R_ZOO_INT_RESP  &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    R_ZOO_RESP                : ', (R_ZOO_RESP      &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'ACTUAL_ZOO_N_TO_C             : ', (ACTUAL_ZOO_N_TO_C &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    ZOO_C                     : ', (ZOO_C             &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    ZOO_N                     : ', (ZOO_N             &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'R_CHEM_AUT_BAC_GROWTH         : ', (R_CHEM_AUT_BAC_GROWTH &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    CHEM_AUT_BAC_N_TO_C       : ',  CHEM_AUT_BAC_N_TO_C

            write(*,*) ''
            write(*,*) 'R_DIA_GROWTH                  : ', (R_DIA_GROWTH   &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    KG_DIA                    : ', (KG_DIA         &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    LIM_KG_DIA                : ', (LIM_KG_DIA     &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '        FDAY                  : ',  (FDAY          &
                                                            (NODES_STRANGE(j)),j=1,nstrange)

            write(*,*) '        K_E                   : ', (K_E            &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '            CHLA              : ', (CHLA           &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '            DIA_C             : ', (DIA_C          &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '            DIA_C_TO_CHLA     : ',  DIA_C_TO_CHLA

            write(*,*) '            K_B_E             : ', (K_B_E          &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '        DEPTH                 : ', (DEPTH          &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '        DISS_OXYGEN           : ', (DISS_OXYGEN    &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '        KHS_O2_DIA            : ',  KHS_O2_DIA

            write(*,*) '        ALPHA_1               : ', (ALPHA_1        &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '        ALPHA_0               : ', (ALPHA_0        &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '        NH4_N                 : ', (NH4_N          &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '        NO3_N                 : ', (NO3_N          &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '        KHS_DIN_DIA           : ',  KHS_DIN_DIA

            write(*,*) '        PO4_P                 : ', (PO4_P          &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '        KHS_DIP_DIA           : ',  KHS_DIP_DIA

            write(*,*) '    DIA_C                     : ', (DIA_C          &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'PREF_NH4N_DIA                 : ', (PREF_NH4N_DIA  &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'DIA_N_TO_C                    : ',  DIA_N_TO_C

            write(*,*) 'R_NITRIFICATION               : ', (R_NITRIFICATION   &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'R_AMMONIA_VOLATIL             : ', (R_AMMONIA_VOLATIL &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'pH                            : ', (pH                &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'TEMP                          : ', (TEMP              &
                                                           (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'K_A_CALC                      : ', (K_A_CALC          &
                                                           (NODES_STRANGE(j)),j=1,nstrange)

            deallocate(STRANGERS)
            deallocate(NODES_STRANGE)
            stop
        end if
    end if

    !NITRATE NITROGEN
    PROCESS_RATES(1:nkn,2, 1) = R_NITRIFICATION + R_ABIOTIC_NITR
    PROCESS_RATES(1:nkn,2, 2) = R_DENITRIFICATION
    PROCESS_RATES(1:nkn,2, 3) = R_DIA_GROWTH         * (1.0D0 - PREF_NH4N_DIA)     * DIA_N_TO_C
    PROCESS_RATES(1:nkn,2, 4) = R_CYN_GROWTH         * (1.0D0 - PREF_NH4N_CYN)     * CYN_N_TO_C
    PROCESS_RATES(1:nkn,2, 5) = R_OPA_GROWTH         * (1.0D0 - PREF_NH4N_OPA)     * OPA_N_TO_C
    PROCESS_RATES(1:nkn,2, 6) = R_NON_FIX_CYN_GROWTH * (1.0D0 - PREF_NH4N_FIX_CYN) * FIX_CYN_N_TO_C
    PROCESS_RATES(1:nkn,2, 7) = R_AER_HET_BAC_GROWTH * (1.0D0 - PREF_NH4N_AER_HET_BAC) * AER_HET_BAC_N_TO_C
    ! Auxiliary
    PROCESS_RATES(1:nkn,2, 8) = PREF_NH4N_DIA
    PROCESS_RATES(1:nkn,2, 9) = PREF_NH4N_CYN
    PROCESS_RATES(1:nkn,2,10) = PREF_NH4N_OPA


    DERIVATIVES(1:nkn,2) = PROCESS_RATES(1:nkn,2, 1) - PROCESS_RATES(1:nkn,2, 2) - PROCESS_RATES(1:nkn,2, 3) - &
                     PROCESS_RATES(1:nkn,2, 4) - PROCESS_RATES(1:nkn,2, 5) - PROCESS_RATES(1:nkn,2, 6)- PROCESS_RATES(1:nkn,2,7)

    !PHOSPHATE PHOSPHORUS
    PROCESS_RATES(1:nkn,3,  1) = R_CHEM_AUT_BAC_TOT_RESP   * CHEM_AUT_BAC_P_TO_C
    PROCESS_RATES(1:nkn,3,  2) = R_AER_HET_BAC_INT_RESP    * AER_HET_BAC_P_TO_C
    PROCESS_RATES(1:nkn,3,  3) = R_FAC_AN_HET_BAC_TOT_RESP * FAC_AN_HET_BAC_P_TO_C
    PROCESS_RATES(1:nkn,3,  4) = R_DIA_TOT_RESP            * DIA_P_TO_C
    PROCESS_RATES(1:nkn,3,  5) = R_CYN_TOT_RESP            * CYN_P_TO_C
    PROCESS_RATES(1:nkn,3,  6) = R_OPA_TOT_RESP            * OPA_P_TO_C
    PROCESS_RATES(1:nkn,3,  7) = R_FIX_CYN_TOT_RESP        * FIX_CYN_P_TO_C
    PROCESS_RATES(1:nkn,3,  8) = R_AER_HET_BAC_P_OX
    PROCESS_RATES(1:nkn,3,  9) = R_FAC_AN_HET_BAC_P_OX
    PROCESS_RATES(1:nkn,3, 10) = R_ZOO_TOT_RESP            * ACTUAL_ZOO_P_TO_C
    PROCESS_RATES(1:nkn,3, 11) = R_CHEM_AUT_BAC_GROWTH     * CHEM_AUT_BAC_P_TO_C
    PROCESS_RATES(1:nkn,3, 12) = R_DIA_GROWTH              * DIA_P_TO_C
    PROCESS_RATES(1:nkn,3, 13) = R_CYN_GROWTH              * CYN_P_TO_C
    PROCESS_RATES(1:nkn,3, 14) = R_OPA_GROWTH              * OPA_P_TO_C
    PROCESS_RATES(1:nkn,3, 15) = R_FIX_CYN_GROWTH          * FIX_CYN_P_TO_C
    PROCESS_RATES(1:nkn,3, 16) = R_ABIOTIC_DOP_MIN
    PROCESS_RATES(1:nkn,3, 17) = R_AER_HET_BAC_GROWTH  * AER_HET_BAC_P_TO_C
    ! Auxiliary
    PROCESS_RATES(1:nkn,3, 18) = K_MIN_DOP_20
    PROCESS_RATES(1:nkn,3, 19) = THETA_K_MIN_DOP
    PROCESS_RATES(1:nkn,3, 20) = TEMP
    PROCESS_RATES(1:nkn,3, 21) = DISS_ORG_P


    DERIVATIVES(1:nkn,3) = PROCESS_RATES(1:nkn,3, 1)  + PROCESS_RATES(1:nkn,3, 2)  + PROCESS_RATES(1:nkn,3, 3)  + &
                     PROCESS_RATES(1:nkn,3, 4)  + PROCESS_RATES(1:nkn,3, 5)  + PROCESS_RATES(1:nkn,3, 6)  + &
                     PROCESS_RATES(1:nkn,3, 7)  + PROCESS_RATES(1:nkn,3, 8)  + PROCESS_RATES(1:nkn,3, 9)  + &
                     PROCESS_RATES(1:nkn,3, 10) - PROCESS_RATES(1:nkn,3, 11) - PROCESS_RATES(1:nkn,3, 12) - &
                     PROCESS_RATES(1:nkn,3, 13) - PROCESS_RATES(1:nkn,3, 14) - PROCESS_RATES(1:nkn,3, 15) + &
                     PROCESS_RATES(1:nkn,3, 16) - PROCESS_RATES(1:nkn,3, 17)

    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,3),VALUE_strange,nkn).eq.1) then

        nstrange = count(VALUE_strange)
        allocate(STRANGERS    (nstrange))
        allocate(NODES_STRANGE(nstrange))

        j=1
            do k=1,nkn
                if(VALUE_strange(k)) then
                    STRANGERS    (j) = DERIVATIVES(k,3)
                    NODES_STRANGE(j) = k
                    j=j+1
                end if
            end do

            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            write(*,*) 'TIME   : ', TIME
            print *, 'PELAGIC_KINETICS:  Derivative(3) '
            print *, 'is NaN or Inf:'
            print *, 'NODE_NUMBERS=',NODES_STRANGE
            print *, 'VALUES=',STRANGERS
            error=1

            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'
            write(*,*) 'R_ABIOTIC_DOP_MIN   : ', (R_ABIOTIC_DOP_MIN(NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    KMIN_DOP_20     : ', K_MIN_DOP_20, GRAT_ZOO_DET_PART_ORG_C
            write(*,*) '    THETA_KMIN_DOP  : ', THETA_K_MIN_DOP, PREF_ZOO_DIA
            write(*,*) '    TEMP            : ', (TEMP(NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    DISS_ORG_P      : ', (DISS_ORG_P(NODES_STRANGE(j)),j=1,nstrange)

            deallocate(STRANGERS)
            deallocate(NODES_STRANGE)
            stop
        end if
    end if

    !DISSOLVED SILICA SILICON
    PROCESS_RATES(1:nkn,20, 1) = R_PART_Si_DISS
    !PROCESS_RATES(1:nkn,20, 2) = R_DIA_TOT_RESP * DIA_SI_TO_C  !Should it come from respiration? Silica is in shell
    PROCESS_RATES(1:nkn,20, 2) = R_DIA_GROWTH   * DIA_SI_TO_C

    DERIVATIVES(1:nkn,20) = PROCESS_RATES(1:nkn,20, 1) - PROCESS_RATES(1:nkn,20, 2)

    !DISSOLVED OXYGEN
    PROCESS_RATES(1:nkn,4, 1)  = R_AERATION
    PROCESS_RATES(1:nkn,4, 2)  = R_DIA_GROWTH       * (1.3 - 0.3*PREF_NH4N_DIA    )* DIA_O2_TO_C        ! formulation from EFDC
    PROCESS_RATES(1:nkn,4, 3)  = R_CYN_GROWTH       * (1.3 - 0.3*PREF_NH4N_CYN    )* CYN_O2_TO_C
    PROCESS_RATES(1:nkn,4, 4)  = R_OPA_GROWTH       * (1.3 - 0.3*PREF_NH4N_OPA    )* OPA_O2_TO_C
    PROCESS_RATES(1:nkn,4, 5)  = R_FIX_CYN_GROWTH   * (1.3 - 0.3*PREF_NH4N_FIX_CYN)* FIX_CYN_O2_TO_C
    PROCESS_RATES(1:nkn,4, 6)  = R_CHEM_AUT_BAC_TOT_RESP * CHEM_AUT_BAC_O2_TO_C
    PROCESS_RATES(1:nkn,4, 7)  = R_AER_HET_BAC_TOT_RESP  * AER_HET_BAC_O2_TO_C
    PROCESS_RATES(1:nkn,4, 8)  = R_DIA_TOT_RESP          * DIA_O2_TO_C
    PROCESS_RATES(1:nkn,4, 9)  = R_CYN_TOT_RESP          * CYN_O2_TO_C
    PROCESS_RATES(1:nkn,4, 10) = R_OPA_TOT_RESP          * OPA_O2_TO_C
    PROCESS_RATES(1:nkn,4, 11) = R_ZOO_TOT_RESP          * ZOO_O2_TO_C
    PROCESS_RATES(1:nkn,4, 12) = (R_NITRIFICATION + R_ABIOTIC_NITR) * 4.57D0
    PROCESS_RATES(1:nkn,4, 13) = R_ABIOTIC_DOC_MIN * (32.D0/16.D0)
    ! Auxiliary
    PROCESS_RATES(1:nkn,4, 14) = K_A_CALC
    PROCESS_RATES(1:nkn,4, 15) = DISS_OXYGEN_SAT

    DERIVATIVES(1:nkn,4) = PROCESS_RATES(1:nkn,4, 1)  + PROCESS_RATES(1:nkn,4, 2)  + PROCESS_RATES(1:nkn,4, 3) + &
                     PROCESS_RATES(1:nkn,4, 4)  + PROCESS_RATES(1:nkn,4, 5)  - PROCESS_RATES(1:nkn,4, 6) - &
                     PROCESS_RATES(1:nkn,4, 7)  - PROCESS_RATES(1:nkn,4, 8)  - PROCESS_RATES(1:nkn,4, 9) - &
                     PROCESS_RATES(1:nkn,4, 10) - PROCESS_RATES(1:nkn,4, 11) - PROCESS_RATES(1:nkn,4, 12)- &
                     PROCESS_RATES(1:nkn,4, 13)

    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,4),VALUE_strange,nkn).eq.1) then
            error=1

            nstrange = count(VALUE_strange)
            allocate(STRANGERS    (nstrange))
            allocate(NODES_STRANGE(nstrange))
            allocate(NODES_STRANGE_int(nstrange))
            allocate(NODES_STRANGE_ext(nstrange))

            j=1

            do k=1,nkn
                if(VALUE_strange(k)) then
                    STRANGERS    (j) = DERIVATIVES(k,4)
                    NODES_STRANGE(j) = k
                    NODES_STRANGE_int(j) = node_active(k)
                    NODES_STRANGE_ext(j) = ipv(node_active(k))
                    j=j+1
                end if
            end do

            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            write(*,*) 'TIME   : ', TIME
            print *, 'PELAGIC_KINETICS:  Derivative(4) '
            print *, 'is NaN or Inf:'
            print *, 'NODE_NUMBERS=',NODES_STRANGE
            print *, 'NODE_NUMBERS int.=',NODES_STRANGE_int
            print *, 'NODE_NUMBERS ext.=',NODES_STRANGE_ext
            print *, 'VALUES=',STRANGERS

            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'
            write(*,*) 'R_AERATION               : ', (R_AERATION      &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'DISS_OXYGEN              : ', (DISS_OXYGEN     &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'DISS_OXYGEN_SAT          : ', (DISS_OXYGEN_SAT &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'K_A parameter            : ',  K_A
            write(*,*) 'K_A_CALC                 : ', (K_A_CALC        &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'R_DIA_GROWTH             : ', (R_DIA_GROWTH &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    DIA_O2_TO_C          : ', DIA_O2_TO_C

            write(*,*) ''
            write(*,*) 'R_CYN_GROWTH             : ', (R_CYN_GROWTH &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    CYN_O2_TO_C          : ', CYN_O2_TO_C
            write(*,*) ''
            write(*,*) 'R_OPA_GROWTH             : ', (R_OPA_GROWTH &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    OPA_O2_TO_C          : ', OPA_O2_TO_C
            write(*,*) ''
            write(*,*) 'R_CHEM_AUT_BAC_TOT_RESP  : ', (R_CHEM_AUT_BAC_TOT_RESP &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    CHEM_AUT_BAC_O2_TO_C : ',  CHEM_AUT_BAC_O2_TO_C
            write(*,*) ''
            write(*,*) 'R_AER_HET_BAC_TOT_RESP   : ', (R_AER_HET_BAC_TOT_RESP  &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    AER_HET_BAC_O2_TO_C  : ', AER_HET_BAC_O2_TO_C
            write(*,*) ''
            write(*,*) 'R_DIA_TOT_RESP           : ', (R_DIA_TOT_RESP          &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'R_CYN_TOT_RESP           : ', (R_CYN_TOT_RESP          &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'R_OPA_TOT_RESP           : ', (R_OPA_TOT_RESP          &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'R_ZOO_EX_DOC             : ', (R_ZOO_EX_DOC            &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'R_ZOO_TOT_RESP           : ', (R_ZOO_TOT_RESP          &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    ZOO_O2_TO_C          : ', ZOO_O2_TO_C
            write(*,*) ''
            write(*,*) 'R_NITRIFICATION          : ', (R_NITRIFICATION         &
                                                      (NODES_STRANGE(j)),j=1,nstrange)

            print *,'TEMP     =', (TEMP       (NODES_STRANGE(j)),j=1,nstrange)
            print *,'SALT     =', (SALT       (NODES_STRANGE(j)),j=1,nstrange)
            print *,'AIRTEMP  =', (AIRTEMP    (NODES_STRANGE(j)),j=1,nstrange)
            print *,'WINDS    =', (WINDS      (NODES_STRANGE(j)),j=1,nstrange)
            print *,'ELEVATION=', (ELEVATION  (NODES_STRANGE(j)),j=1,nstrange)
            print *,'DEPTH    =', (DEPTH      (NODES_STRANGE(j)),j=1,nstrange)
            print *,'SURFACE_BOX:', SURFACE_BOX
            print *,'FLAGS(2)=',FLAGS(2)

            print *,'R_AERATION '     ,(R_AERATION     (j),j=1,nkn)
            print *,'DISS_OXYGEN_SAT ',(DISS_OXYGEN_SAT(j),j=1,nkn)
            print *,'K_A_CALC '       ,(K_A_CALC       (j),j=1,nkn)
            print *,'DISS_OXYGEN'     ,(DISS_OXYGEN    (j),j=1,nkn)
            stop
        end if
    end if

    !CHEMOAUTOTROPHIC BACTERIA CARBON
    PROCESS_RATES(1:nkn,5, 1) = R_CHEM_AUT_BAC_GROWTH
    PROCESS_RATES(1:nkn,5, 2) = R_CHEM_AUT_BAC_TOT_RESP
    PROCESS_RATES(1:nkn,5, 3) = R_CHEM_AUT_BAC_DEATH
    PROCESS_RATES(1:nkn,5, 4) = R_ZOO_FEEDING_CHEM_AUT_BAC
    ! Auxiliary
    PROCESS_RATES(1:nkn,5, 5) = LIM_TEMP_CHEM_AUT_BAC
    PROCESS_RATES(1:nkn,5, 6) = LIM_NH4_N_CHEM_AUT_BAC
    PROCESS_RATES(1:nkn,5, 7) = LIM_PO4_P_CHEM_AUT_BAC
    PROCESS_RATES(1:nkn,5, 8) = LIM_OXY_CHEM_AUT_BAC

    DERIVATIVES(1:nkn,5) = PROCESS_RATES(1:nkn,5, 1) - PROCESS_RATES(1:nkn,5, 2) - PROCESS_RATES(1:nkn,5, 3) - &
                     PROCESS_RATES(1:nkn,5, 4)

    !AEROBIC HETEROTROPHIC BACTERIA CARBON
    PROCESS_RATES(1:nkn,6, 1) = R_AER_HET_BAC_GROWTH
    !    PROCESS_RATES(1:nkn,6, 2) = R_AER_HET_BAC_TOT_RESP
    PROCESS_RATES(1:nkn,6, 2) = R_AER_HET_BAC_INT_RESP
    PROCESS_RATES(1:nkn,6, 3) = R_AER_HET_BAC_DEATH
    PROCESS_RATES(1:nkn,6, 4) = R_ZOO_FEEDING_AER_HET_BAC
    ! Auxiliary
    PROCESS_RATES(1:nkn,6, 5) = LIM_TEMP_AER_HET_BAC
    PROCESS_RATES(1:nkn,6, 6) = LIM_DISS_ORG_C_AER_HET_BAC
    PROCESS_RATES(1:nkn,6, 7) = LIM_DISS_ORG_N_AER_HET_BAC
    PROCESS_RATES(1:nkn,6, 8) = LIM_DISS_ORG_P_AER_HET_BAC
    PROCESS_RATES(1:nkn,6, 9) = LIM_DIN_AER_HET_BAC
    PROCESS_RATES(1:nkn,6,10) = LIM_DIP_AER_HET_BAC
    PROCESS_RATES(1:nkn,6,11) = LIM_OXY_AER_HET_BAC
    PROCESS_RATES(1:nkn,6,12) = LIM_PHYT_C_AER_HET_BAC
    PROCESS_RATES(1:nkn,6,13) = R_AER_HET_BAC_TOT_RESP

    DERIVATIVES(1:nkn,6) = PROCESS_RATES(1:nkn,6, 1) - PROCESS_RATES(1:nkn,6, 2) - PROCESS_RATES(1:nkn,6, 3) - &
                     PROCESS_RATES(1:nkn,6, 4)

    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,6),VALUE_strange,nkn).eq.1) then
            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            print *, 'PELAGIC_KINETICS:'
            write(*,*) 'TIME   : ', TIME

            write(*,*) 'DERIVATIVES(6) is not a number or infinite : ', DERIVATIVES(1:nkn,6)
            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'
            write(*,*) 'R_AER_HET_BAC_GROWTH         : ', R_AER_HET_BAC_GROWTH
            write(*,*) 'R_AER_HET_BAC_TOT_RESP       : ', R_AER_HET_BAC_TOT_RESP
            write(*,*) 'R_AER_HET_BAC_DEATH          : ', R_AER_HET_BAC_DEATH
            write(*,*) 'R_ZOO_FEEDING_AER_HET_BAC    : ', R_ZOO_FEEDING_AER_HET_BAC
            write(*,*) 'LIM_TEMP_AER_HET_BAC         : ', LIM_TEMP_AER_HET_BAC
            write(*,*) 'LIM_DISS_ORG_C_AER_HET_BAC   : ', LIM_DISS_ORG_C_AER_HET_BAC
            write(*,*) 'LIM_DIN_AER_HET_BAC          : ', LIM_DIN_AER_HET_BAC
            write(*,*) 'LIM_DIP_AER_HET_BAC          : ', LIM_DIP_AER_HET_BAC
            write(*,*) 'LIM_OXY_AER_HET_BAC          : ', LIM_OXY_AER_HET_BAC
            stop
        end if
    end if

    !FACULTATIVE ANAEROBIC HETEROTROPHIC BACTERIA CARBON Limitation from nutrients fixme
    PROCESS_RATES(1:nkn,7, 1) = R_FAC_AN_HET_BAC_GROWTH
    PROCESS_RATES(1:nkn,7, 2) = R_FAC_AN_HET_BAC_TOT_RESP
    PROCESS_RATES(1:nkn,7, 3) = R_FAC_AN_HET_BAC_DEATH
    PROCESS_RATES(1:nkn,7, 4) = R_ZOO_FEEDING_FAC_AN_HET_BAC
    ! Auxiliary
    PROCESS_RATES(1:nkn,7, 5) = LIM_TEMP_FAC_AN_HET_BAC
    PROCESS_RATES(1:nkn,7, 6) = LIM_DISS_ORG_C_FAC_AN_HET_BAC
    PROCESS_RATES(1:nkn,7, 7) = LIM_DISS_ORG_N_FAC_AN_HET_BAC
    PROCESS_RATES(1:nkn,7, 8) = LIM_DISS_ORG_P_FAC_AN_HET_BAC
    PROCESS_RATES(1:nkn,7, 9) = LIM_OXY_FAC_AN_HET_BAC

    DERIVATIVES(1:nkn,7) = PROCESS_RATES(1:nkn,7, 1) - PROCESS_RATES(1:nkn,7, 2) - PROCESS_RATES(1:nkn,7, 3) - &
                     PROCESS_RATES(1:nkn,7, 4)

    !DIATOMS CARBON
    PROCESS_RATES(1:nkn,8, 1)  = R_DIA_GROWTH
    PROCESS_RATES(1:nkn,8, 2)  = R_DIA_TOT_RESP
    PROCESS_RATES(1:nkn,8, 3)  = R_DIA_DEATH
    PROCESS_RATES(1:nkn,8, 4)  = R_ZOO_FEEDING_DIA
    ! Auxiliary
    PROCESS_RATES(1:nkn,8, 5)  = LIM_KG_DIA_TEMP
    PROCESS_RATES(1:nkn,8, 6)  = LIM_KG_DIA_DOXY
    PROCESS_RATES(1:nkn,8, 7)  = LIM_KG_DIA_N
    PROCESS_RATES(1:nkn,8, 8)  = LIM_KG_DIA_P
    PROCESS_RATES(1:nkn,8, 9)  = LIM_KG_DIA_DISS_Si
    PROCESS_RATES(1:nkn,8, 10) = LIM_KG_DIA_LIGHT
    PROCESS_RATES(1:nkn,8, 11) = DIA_LIGHT_SAT

    DERIVATIVES(1:nkn,8) = PROCESS_RATES(1:nkn,8, 1) - PROCESS_RATES(1:nkn,8, 2) - PROCESS_RATES(1:nkn,8, 3) - &
                     PROCESS_RATES(1:nkn,8, 4)

    !NON-NITROGEN FIXING CYANOBACTERIA CARBON
    PROCESS_RATES(1:nkn,18, 1) = R_CYN_GROWTH
    PROCESS_RATES(1:nkn,18, 2) = R_CYN_TOT_RESP
    PROCESS_RATES(1:nkn,18, 3) = R_CYN_DEATH
    PROCESS_RATES(1:nkn,18, 4) = R_ZOO_FEEDING_CYN
    ! Auxiliary
    PROCESS_RATES(1:nkn,18, 5) = LIM_KG_CYN_TEMP
    PROCESS_RATES(1:nkn,18, 6) = LIM_KG_CYN_DOXY
    PROCESS_RATES(1:nkn,18, 7) = LIM_KG_CYN_N
    PROCESS_RATES(1:nkn,18, 8) = LIM_KG_CYN_P
    PROCESS_RATES(1:nkn,18, 9) = LIM_KG_CYN_LIGHT
    PROCESS_RATES(1:nkn,18,10) = I_A             !light langlays
    PROCESS_RATES(1:nkn,18, 11) = CYN_LIGHT_SAT
    PROCESS_RATES(1:nkn,18, 12) = TEMP

    DERIVATIVES(1:nkn,18) = PROCESS_RATES(1:nkn,18, 1) - PROCESS_RATES(1:nkn,18, 2) - PROCESS_RATES(1:nkn,18, 3) - &
                      PROCESS_RATES(1:nkn,18, 4)

    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,18),VALUE_strange,nkn).eq.1) then
            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            print *, 'PELAGIC_KINETICS:'
            write(*,*) 'TIME   : ', TIME
            write(*,*) 'DERIVATIVES(18) is nor a number or infinite.'
            write(*,*) DERIVATIVES(1:nkn,18)
            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'
            write(*,*) 'R_CYN_GROWTH             : ', R_CYN_GROWTH
            write(*,*) '    LIM_KG_CYN           : ', LIM_KG_CYN
            write(*,*) '        LIM_KG_CYN_DOXY  : ', LIM_KG_CYN_DOXY
            write(*,*) '        LIM_KG_CYN_NUTR  : ', LIM_KG_CYN_NUTR
            write(*,*) '        LIM_KG_CYN_LIGHT : ', LIM_KG_CYN_LIGHT
            write(*,*) '            ALPHA_0      : ', ALPHA_0
            write(*,*) '            ALPHA_1      : ', ALPHA_1
            write(*,*) '            I_A          : ', I_A
            write(*,*) '            I_S_CYN      : ', I_S_CYN
            write(*,*) '    CYN_C                : ', CYN_C
            write(*,*) '    KG_CYN               : ', KG_CYN
            write(*,*) ''
            write(*,*) 'R_CYN_TOT_RESP           : ', R_CYN_TOT_RESP
            write(*,*) ''
            write(*,*) 'R_CYN_DEATH              : ', R_CYN_DEATH
            write(*,*) ''
            write(*,*) 'R_ZOO_FEEDING_CYN        : ', R_ZOO_FEEDING_CYN
            stop
        end if
    end if

    !NITROGEN FIXING CYANOBACTERIA CARBON
    PROCESS_RATES(1:nkn,22, 1)  = R_FIX_CYN_GROWTH
    PROCESS_RATES(1:nkn,22, 2)  = R_FIX_CYN_TOT_RESP
    PROCESS_RATES(1:nkn,22, 3)  = R_FIX_CYN_DEATH
    PROCESS_RATES(1:nkn,22, 4)  = R_ZOO_FEEDING_FIX_CYN
    ! Auxiliary
    PROCESS_RATES(1:nkn,22, 5)  = R_NON_FIX_CYN_GROWTH
    PROCESS_RATES(1:nkn,22, 6)  = R_FIX_FIX_CYN_GROWTH
    PROCESS_RATES(1:nkn,22, 7)  = R_FIX_FIX_CYN_GROWTH * FIX_CYN_N_TO_C !Nitrogen fixation
    PROCESS_RATES(1:nkn,22, 8)  = LIM_KG_FIX_CYN_TEMP
    PROCESS_RATES(1:nkn,22, 9)  = LIM_KG_FIX_CYN_DOXY
    PROCESS_RATES(1:nkn,22, 10)  = LIM_KG_FIX_FIX_CYN_N
    PROCESS_RATES(1:nkn,22, 11) = LIM_KG_FIX_FIX_CYN_P
    PROCESS_RATES(1:nkn,22, 12) = LIM_KG_NON_FIX_CYN_N
    PROCESS_RATES(1:nkn,22, 13) = LIM_KG_NON_FIX_CYN_P
    PROCESS_RATES(1:nkn,22, 14) = LIM_KG_FIX_CYN_LIGHT
    PROCESS_RATES(1:nkn,22, 15) = ((NH4_N + NO3_N)/14.D0)/(PO4_P/31.D0) ! NP molar ratio
    PROCESS_RATES(1:nkn,22, 16) = NH4_N + NO3_N
    PROCESS_RATES(1:nkn,22, 17) = FIX_CYN_LIGHT_SAT

    DERIVATIVES(1:nkn,22) = PROCESS_RATES(1:nkn,22, 1) - PROCESS_RATES(1:nkn,22, 2) - PROCESS_RATES(1:nkn,22, 3) - &
                      PROCESS_RATES(1:nkn,22, 4)

    !OTHER PLANKTONIC ALGAE CARBON
    PROCESS_RATES(1:nkn,19, 1) = R_OPA_GROWTH
    PROCESS_RATES(1:nkn,19, 2) = R_OPA_TOT_RESP
    PROCESS_RATES(1:nkn,19, 3) = R_OPA_DEATH
    PROCESS_RATES(1:nkn,19, 4) = R_ZOO_FEEDING_OPA
    ! Auxiliary
    PROCESS_RATES(1:nkn,19, 5) = LIM_KG_OPA_TEMP
    PROCESS_RATES(1:nkn,19, 6) = LIM_KG_OPA_DOXY
    PROCESS_RATES(1:nkn,19, 7) = LIM_KG_OPA_N
    PROCESS_RATES(1:nkn,19, 8) = LIM_KG_OPA_P
    PROCESS_RATES(1:nkn,19, 9) = LIM_KG_OPA_LIGHT
    PROCESS_RATES(1:nkn,19, 10) = OPA_LIGHT_SAT

    DERIVATIVES(1:nkn,19) = PROCESS_RATES(1:nkn,19, 1) - PROCESS_RATES(1:nkn,19, 2) - PROCESS_RATES(1:nkn,19, 3) - &
                      PROCESS_RATES(1:nkn,19, 4)

    !ZOOPLANKTON CARBON
    PROCESS_RATES(1:nkn,9, 1)  = R_ZOO_GROWTH
    PROCESS_RATES(1:nkn,9, 2)  = 0.   !No exrection for a while
    PROCESS_RATES(1:nkn,9, 3)  = R_ZOO_TOT_RESP
    PROCESS_RATES(1:nkn,9, 4)  = R_ZOO_DEATH

    PROCESS_RATES(1:nkn,9, 5)  = R_ZOO_FEEDING_DIA
    PROCESS_RATES(1:nkn,9, 6)  = R_ZOO_FEEDING_CYN
    PROCESS_RATES(1:nkn,9, 7)  = R_ZOO_FEEDING_OPA
    PROCESS_RATES(1:nkn,9, 8)  = R_ZOO_FEEDING_FIX_CYN
    PROCESS_RATES(1:nkn,9, 9)  = R_ZOO_FEEDING_CHEM_AUT_BAC
    PROCESS_RATES(1:nkn,9, 10) = R_ZOO_FEEDING_AER_HET_BAC
    PROCESS_RATES(1:nkn,9, 11) = R_ZOO_FEEDING_FAC_AN_HET_BAC
    PROCESS_RATES(1:nkn,9, 12) = R_ZOO_FEEDING_DET_PART_ORG_C


    DERIVATIVES(1:nkn,9) = PROCESS_RATES(1:nkn,9, 1) - PROCESS_RATES(1:nkn,9, 2) - PROCESS_RATES(1:nkn,9, 3) - &
                     PROCESS_RATES(1:nkn,9, 4)

    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,9),VALUE_strange,nkn).eq.1) then
            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            print *, 'PELAGIC_KINETICS:'
            write(*,*) 'TIME   : ', TIME
            write(*,*) 'DERIVATIVES(9) is not a number or infinite'
            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'
            write(*,*) 'ZOO_C                                  : ', ZOO_C
            write(*,*) 'R_ZOO_GROWTH                           : ', R_ZOO_GROWTH
            write(*,*) '    R_ZOO_FEEDING_DIA                  : ', R_ZOO_FEEDING_DIA
            write(*,*) '    R_ZOO_FEEDING_CYN                  : ', R_ZOO_FEEDING_CYN
            write(*,*) '    R_ZOO_FEEDING_OPA                  : ', R_ZOO_FEEDING_OPA
            write(*,*) '    R_ZOO_FEEDING_FIX_CYN              : ', R_ZOO_FEEDING_FIX_CYN
            write(*,*) '    R_ZOO_FEEDING_CHEM_AUT_BAC         : ', R_ZOO_FEEDING_CHEM_AUT_BAC
            write(*,*) '    R_ZOO_FEEDING_AER_HET_BAC          : ', R_ZOO_FEEDING_AER_HET_BAC
            write(*,*) '    R_ZOO_FEEDING_FAC_AN_HET_BAC       : ', R_ZOO_FEEDING_FAC_AN_HET_BAC
            write(*,*) '    R_ZOO_FEEDING_DET_PART_ORG_C       : ', R_ZOO_FEEDING_DET_PART_ORG_C
            write(*,*) '        KG_ZOO_DET_PART_ORG_C          : ', KG_ZOO_DET_PART_ORG_C
            write(*,*) '        FOOD_FACTOR_ZOO_DET_PART_ORG_C : ', FOOD_FACTOR_ZOO_DET_PART_ORG_C
            write(*,*) '            PREF_ZOO_DET_PART_ORG_C    : ', PREF_ZOO_DET_PART_ORG_C
            write(*,*) '            DET_PART_ORG_C             : ', DET_PART_ORG_C
            write(*,*) '            FOOD_MIN_ZOO               : ', FOOD_MIN_ZOO
            write(*,*) '            FOOD_AVAIL_ZOO             : ', FOOD_AVAIL_ZOO
            write(*,*) '            KHS_DET_PART_ORG_C_ZOO     : ', KHS_DET_PART_ORG_C_ZOO
            write(*,*) 'R_ZOO_EX_DOC                           : ', R_ZOO_EX_DOC
            write(*,*) 'R_ZOO_TOT_RESP                         : ', R_ZOO_TOT_RESP
            write(*,*) 'R_ZOO_DEATH                            : ', R_ZOO_DEATH
            stop
        end if
    end if

    !ZOOPLANKTON NITROGEN
    PROCESS_RATES(1:nkn,10, 1)  = R_ZOO_FEEDING_DIA            * DIA_N_TO_C
    PROCESS_RATES(1:nkn,10, 2)  = R_ZOO_FEEDING_CYN            * CYN_N_TO_C
    PROCESS_RATES(1:nkn,10, 3)  = R_ZOO_FEEDING_OPA            * OPA_N_TO_C
    PROCESS_RATES(1:nkn,10, 4)  = R_ZOO_FEEDING_FIX_CYN        * FIX_CYN_N_TO_C
    PROCESS_RATES(1:nkn,10, 5)  = R_ZOO_FEEDING_CHEM_AUT_BAC   * CHEM_AUT_BAC_N_TO_C
    PROCESS_RATES(1:nkn,10, 6)  = R_ZOO_FEEDING_AER_HET_BAC    * AER_HET_BAC_N_TO_C
    PROCESS_RATES(1:nkn,10, 7)  = R_ZOO_FEEDING_FAC_AN_HET_BAC * FAC_AN_HET_BAC_N_TO_C
    PROCESS_RATES(1:nkn,10, 8)  = R_ZOO_FEEDING_DET_PART_ORG_C * ACTUAL_DET_N_TO_C
    PROCESS_RATES(1:nkn,10, 9)  = 0. !R_ZOO_EX_DON No excretion for a while
    PROCESS_RATES(1:nkn,10, 10) = R_ZOO_TOT_RESP               * ACTUAL_ZOO_N_TO_C
    PROCESS_RATES(1:nkn,10, 11) = R_ZOO_DEATH                  * ACTUAL_ZOO_N_TO_C
    PROCESS_RATES(1:nkn,10, 12) = ACTUAL_ZOO_N_TO_C

    !     DERIVATIVES(1:nkn,10) = PROCESS_RATES(1:nkn,10, 1)  + PROCESS_RATES(1:nkn,10, 2) + PROCESS_RATES(1:nkn,10, 3) + &
    !                       PROCESS_RATES(1:nkn,10, 4)  + PROCESS_RATES(1:nkn,10, 5) + PROCESS_RATES(1:nkn,10, 6) + &
    !                       PROCESS_RATES(1:nkn,10, 7)  + PROCESS_RATES(1:nkn,10, 8) - PROCESS_RATES(1:nkn,10, 9) - &
    !                       PROCESS_RATES(1:nkn,10, 10) - PROCESS_RATES(1:nkn,10, 11)
    DERIVATIVES(1:nkn,10) = DERIVATIVES(1:nkn,9)*ACTUAL_ZOO_N_TO_C
    !ZOOPLANKTON PHOSPHORUS
    PROCESS_RATES(1:nkn,11, 1)  = R_ZOO_FEEDING_DIA            * DIA_P_TO_C
    PROCESS_RATES(1:nkn,11, 2)  = R_ZOO_FEEDING_CYN            * CYN_P_TO_C
    PROCESS_RATES(1:nkn,11, 3)  = R_ZOO_FEEDING_OPA            * OPA_P_TO_C
    PROCESS_RATES(1:nkn,11, 4)  = R_ZOO_FEEDING_FIX_CYN        * FIX_CYN_P_TO_C
    PROCESS_RATES(1:nkn,11, 5)  = R_ZOO_FEEDING_CHEM_AUT_BAC   * CHEM_AUT_BAC_P_TO_C
    PROCESS_RATES(1:nkn,11, 6)  = R_ZOO_FEEDING_AER_HET_BAC    * AER_HET_BAC_P_TO_C
    PROCESS_RATES(1:nkn,11, 7)  = R_ZOO_FEEDING_FAC_AN_HET_BAC * FAC_AN_HET_BAC_P_TO_C
    PROCESS_RATES(1:nkn,11, 8)  = R_ZOO_FEEDING_DET_PART_ORG_C * ACTUAL_DET_P_TO_C
    PROCESS_RATES(1:nkn,11, 9)  = 0. ! R_ZOO_EX_DOP no excretion
    PROCESS_RATES(1:nkn,11, 10) = R_ZOO_TOT_RESP                * ACTUAL_ZOO_P_TO_C
    PROCESS_RATES(1:nkn,11, 11) = R_ZOO_DEATH                   * ACTUAL_ZOO_P_TO_C
    PROCESS_RATES(1:nkn,11, 12) = ACTUAL_ZOO_P_TO_C

    !     DERIVATIVES(1:nkn,11) = PROCESS_RATES(1:nkn,11, 1)  + PROCESS_RATES(1:nkn,11, 2) + PROCESS_RATES(1:nkn,11, 3) + &
    !                       PROCESS_RATES(1:nkn,11, 4)  + PROCESS_RATES(1:nkn,11, 5) + PROCESS_RATES(1:nkn,11, 6) + &
    !                       PROCESS_RATES(1:nkn,11, 7)  + PROCESS_RATES(1:nkn,11, 8) - PROCESS_RATES(1:nkn,11, 9) - &
    !                       PROCESS_RATES(1:nkn,11, 10) - PROCESS_RATES(1:nkn,11, 11)
    DERIVATIVES(1:nkn,11) =  DERIVATIVES(1:nkn,9) * ACTUAL_ZOO_P_TO_C

    !DEAD ORGANIC CARBON PARTICLES
    PROCESS_RATES(1:nkn,12, 1)  = R_DIA_DEATH
    PROCESS_RATES(1:nkn,12, 2)  = R_CYN_DEATH
    PROCESS_RATES(1:nkn,12, 3)  = R_OPA_DEATH
    PROCESS_RATES(1:nkn,12, 4)  = R_FIX_CYN_DEATH
    PROCESS_RATES(1:nkn,12, 5)  = R_CHEM_AUT_BAC_DEATH
    PROCESS_RATES(1:nkn,12, 6)  = R_AER_HET_BAC_DEATH
    PROCESS_RATES(1:nkn,12, 7)  = R_FAC_AN_HET_BAC_DEATH
    PROCESS_RATES(1:nkn,12, 8)  = R_ZOO_DEATH
    PROCESS_RATES(1:nkn,12, 9)  = R_ZOO_FEEDING_DET_PART_ORG_C
    PROCESS_RATES(1:nkn,12, 10) = R_DET_PART_ORG_C_DISSOLUTION

    DERIVATIVES(1:nkn,12) = PROCESS_RATES(1:nkn,12, 1) + PROCESS_RATES(1:nkn,12, 2) + PROCESS_RATES(1:nkn,12, 3) + &
                      PROCESS_RATES(1:nkn,12, 4) + PROCESS_RATES(1:nkn,12, 5) + PROCESS_RATES(1:nkn,12, 6) + &
                      PROCESS_RATES(1:nkn,12, 7) + PROCESS_RATES(1:nkn,12, 8) - PROCESS_RATES(1:nkn,12, 9) - &
                      PROCESS_RATES(1:nkn,12, 10)

    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,12),VALUE_strange,nkn).eq.1) then
            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            print *, 'PELAGIC_KINETICS:'
            write(*,*) 'TIME   : ', TIME
            write(*,*) 'DERIVATIVES(12) is not a number or infinite.'
            write(*,*) DERIVATIVES(1:nkn,12)
            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'
            write(*,*) 'R_DIA_DEATH                        : ', R_DIA_DEATH
            write(*,*) 'R_CYN_DEATH                        : ', R_CYN_DEATH
            write(*,*) 'R_OPA_DEATH                        : ', R_OPA_DEATH
            write(*,*) 'R_FIX_CYN_DEATH                    : ', R_FIX_CYN_DEATH
            write(*,*) 'R_CHEM_AUT_BAC_DEATH               : ', R_CHEM_AUT_BAC_DEATH
            write(*,*) 'R_AER_HET_BAC_DEATH                : ', R_AER_HET_BAC_DEATH
            write(*,*) 'R_FAC_AN_HET_BAC_DEATH             : ', R_FAC_AN_HET_BAC_DEATH
            write(*,*) '    KD_FAC_AN_HET_BAC              : ', KD_FAC_AN_HET_BAC
            write(*,*) '    NO3N_LACK_STR_FAC_AN_HET_BAC_D : ', NO3N_LACK_STR_FAC_AN_HET_BAC_D
            write(*,*) '    FAC_AN_HET_BAC_C               : ', FAC_AN_HET_BAC_C
            write(*,*) 'R_ZOO_DEATH                        : ', R_ZOO_DEATH
            write(*,*) 'R_ZOO_FEEDING_DET_PART_ORG_C       : ', R_ZOO_FEEDING_DET_PART_ORG_C
            write(*,*) 'R_DET_PART_ORG_C_DISSOLUTION       : ', R_DET_PART_ORG_C_DISSOLUTION
            write(*,*) '    KDISS_DET_PART_ORG_C_20        : ', KDISS_DET_PART_ORG_C_20
            write(*,*) '    THETA_KDISS_DET_PART_ORG_C     : ', THETA_KDISS_DET_PART_ORG_C
            write(*,*) '    DET_PART_ORG_C                 : ', DET_PART_ORG_C
            stop
        end if
    end if

    !DEAD ORGANIC NITROGEN PARTICLES
    PROCESS_RATES(1:nkn,13, 1)  = R_DIA_DEATH                  * DIA_N_TO_C
    PROCESS_RATES(1:nkn,13, 2)  = R_CYN_DEATH                  * CYN_N_TO_C
    PROCESS_RATES(1:nkn,13, 3)  = R_OPA_DEATH                  * OPA_N_TO_C
    PROCESS_RATES(1:nkn,13, 4)  = R_FIX_CYN_DEATH              * FIX_CYN_N_TO_C
    PROCESS_RATES(1:nkn,13, 5)  = R_CHEM_AUT_BAC_DEATH         * CHEM_AUT_BAC_N_TO_C
    PROCESS_RATES(1:nkn,13, 6)  = R_AER_HET_BAC_DEATH          * AER_HET_BAC_N_TO_C
    PROCESS_RATES(1:nkn,13, 7)  = R_FAC_AN_HET_BAC_DEATH       * FAC_AN_HET_BAC_N_TO_C
    PROCESS_RATES(1:nkn,13, 8)  = R_ZOO_DEATH                  * ACTUAL_ZOO_N_TO_C
    PROCESS_RATES(1:nkn,13, 9)  = R_ZOO_FEEDING_DET_PART_ORG_C * ACTUAL_DET_N_TO_C
    PROCESS_RATES(1:nkn,13, 10) = R_DET_PART_ORG_N_DISSOLUTION
    ! Auxiliary
    PROCESS_RATES(1:nkn,13, 11) = ACTUAL_DET_N_TO_C

    DERIVATIVES(1:nkn,13) = PROCESS_RATES(1:nkn,13, 1) + PROCESS_RATES(1:nkn,13, 2) + PROCESS_RATES(1:nkn,13, 3) + &
                      PROCESS_RATES(1:nkn,13, 4) + PROCESS_RATES(1:nkn,13, 5) + PROCESS_RATES(1:nkn,13, 6) + &
                      PROCESS_RATES(1:nkn,13, 7) + PROCESS_RATES(1:nkn,13, 8) - PROCESS_RATES(1:nkn,13, 9) - &
                      PROCESS_RATES(1:nkn,13, 10)

    !DEAD ORGANIC PHOSPHORUS PARTICLES
    PROCESS_RATES(1:nkn,14, 1)  = R_DIA_DEATH                  * DIA_P_TO_C
    PROCESS_RATES(1:nkn,14, 2)  = R_CYN_DEATH                  * CYN_P_TO_C
    PROCESS_RATES(1:nkn,14, 3)  = R_OPA_DEATH                  * OPA_P_TO_C
    PROCESS_RATES(1:nkn,14, 4)  = R_FIX_CYN_DEATH              * FIX_CYN_P_TO_C
    PROCESS_RATES(1:nkn,14, 5)  = R_CHEM_AUT_BAC_DEATH         * CHEM_AUT_BAC_P_TO_C
    PROCESS_RATES(1:nkn,14, 6)  = R_AER_HET_BAC_DEATH          * AER_HET_BAC_P_TO_C
    PROCESS_RATES(1:nkn,14, 7)  = R_FAC_AN_HET_BAC_DEATH       * FAC_AN_HET_BAC_P_TO_C
    PROCESS_RATES(1:nkn,14, 8)  = R_ZOO_DEATH                  * ACTUAL_ZOO_P_TO_C
    PROCESS_RATES(1:nkn,14, 9)  = R_ZOO_FEEDING_DET_PART_ORG_C * ACTUAL_DET_P_TO_C
    PROCESS_RATES(1:nkn,14, 10) = R_DET_PART_ORG_P_DISSOLUTION
    ! Auxiliary
    PROCESS_RATES(1:nkn,14, 11) = ACTUAL_DET_P_TO_C

    DERIVATIVES(1:nkn,14) = PROCESS_RATES(1:nkn,14, 1) + PROCESS_RATES(1:nkn,14, 2) + PROCESS_RATES(1:nkn,14, 3) + &
                      PROCESS_RATES(1:nkn,14, 4) + PROCESS_RATES(1:nkn,14, 5) + PROCESS_RATES(1:nkn,14, 6) + &
                      PROCESS_RATES(1:nkn,14, 7) + PROCESS_RATES(1:nkn,14, 8) - PROCESS_RATES(1:nkn,14, 9) - &
                      PROCESS_RATES(1:nkn,14, 10)

    !PARTICULATE SILICA
    PROCESS_RATES(1:nkn,21, 1) = R_DIA_DEATH
    PROCESS_RATES(1:nkn,21, 2) = R_ZOO_FEEDING_DIA * DIA_Si_TO_C
    PROCESS_RATES(1:nkn,21, 3) = R_PART_Si_DISS

    DERIVATIVES(1:nkn,21) = PROCESS_RATES(1:nkn,21, 1) + PROCESS_RATES(1:nkn,21, 2) - PROCESS_RATES(1:nkn,21, 3)

    !DISSOLVED ORGANIC CARBON
    PROCESS_RATES(1:nkn,15, 1) = R_DET_PART_ORG_C_DISSOLUTION
    PROCESS_RATES(1:nkn,15, 2) = 0.   !No excretion for a while
    PROCESS_RATES(1:nkn,15, 3) = R_AER_HET_BAC_GROWTH/YIELD_OC_AER_HET_BAC
    PROCESS_RATES(1:nkn,15, 4) = R_FAC_AN_HET_BAC_GROWTH
    PROCESS_RATES(1:nkn,15, 5) = R_ABIOTIC_DOC_MIN
    PROCESS_RATES(1:nkn,15, 6) = R_DIA_EXCR + R_CYN_EXCR + R_FIX_CYN_EXCR + R_OPA_EXCR

    DERIVATIVES(1:nkn,15) = PROCESS_RATES(1:nkn,15, 1) + PROCESS_RATES(1:nkn,15, 2) - PROCESS_RATES(1:nkn,15, 3) - &
                      PROCESS_RATES(1:nkn,15, 4) - PROCESS_RATES(1:nkn,15, 5) + PROCESS_RATES(1:nkn,15, 6)

    !DISSOLVED ORGANIC NITROGEN
    PROCESS_RATES(1:nkn,16, 1) = R_DET_PART_ORG_N_DISSOLUTION
    PROCESS_RATES(1:nkn,16, 2) = 0. ! R_ZOO_EX_DON No excretion for a while
    PROCESS_RATES(1:nkn,16, 3) = R_AER_HET_BAC_N_OX
    PROCESS_RATES(1:nkn,16, 4) = R_FAC_AN_HET_BAC_GROWTH * FAC_AN_HET_BAC_N_TO_C
    PROCESS_RATES(1:nkn,16, 5) = R_ABIOTIC_DON_MIN
    PROCESS_RATES(1:nkn,16, 6) = R_DIA_EXCR *DIA_N_TO_C + R_CYN_EXCR * CYN_N_TO_C + &
                                 R_FIX_CYN_EXCR * FIX_CYN_N_TO_C + R_OPA_EXCR * OPA_N_TO_C

    DERIVATIVES(1:nkn,16) = PROCESS_RATES(1:nkn,16, 1) + PROCESS_RATES(1:nkn,16, 2) - PROCESS_RATES(1:nkn,16, 3) - &
                      PROCESS_RATES(1:nkn,16, 4) - PROCESS_RATES(1:nkn,16, 5) + PROCESS_RATES(1:nkn,16, 6)

    !DISSOLVED ORGANIC PHOSPHORUS
    PROCESS_RATES(1:nkn,17, 1) = R_DET_PART_ORG_P_DISSOLUTION
    PROCESS_RATES(1:nkn,17, 2) = 0. !R_ZOO_EX_DOP No excretion
    PROCESS_RATES(1:nkn,17, 3) = R_AER_HET_BAC_P_OX
    PROCESS_RATES(1:nkn,17, 4) = R_FAC_AN_HET_BAC_GROWTH * FAC_AN_HET_BAC_P_TO_C
    PROCESS_RATES(1:nkn,17, 5) = R_ABIOTIC_DOP_MIN
    PROCESS_RATES(1:nkn,17, 6) = R_DIA_EXCR *DIA_P_TO_C + R_CYN_EXCR * CYN_P_TO_C + &
                                 R_FIX_CYN_EXCR * FIX_CYN_P_TO_C + R_OPA_EXCR * OPA_P_TO_C

    DERIVATIVES(1:nkn,17) = PROCESS_RATES(1:nkn,17, 1) + PROCESS_RATES(1:nkn,17, 2) - PROCESS_RATES(1:nkn,17, 3) - &
                      PROCESS_RATES(1:nkn,17, 4) - PROCESS_RATES(1:nkn,17, 5) + PROCESS_RATES(1:nkn,17, 6)


    ! Kinetic sub model for dissolved inorganic carbon

    ! Sources
    R_CHEM_AUT_BAC_TOT_RESP    = PROCESS_RATES(1:nkn,5 , 2)
    R_AER_HET_BAC_TOT_RESP     = PROCESS_RATES(1:nkn,6 , 2)
    R_FAC_AN_HET_BAC_TOT_RESP  = PROCESS_RATES(1:nkn,7 , 2)
    R_DIA_TOT_RESP             = PROCESS_RATES(1:nkn,8 , 2)
    R_CYN_TOT_RESP             = PROCESS_RATES(1:nkn,18, 2)
    R_FIX_CYN_TOT_RESP         = PROCESS_RATES(1:nkn,22, 2)
    R_OPA_TOT_RESP             = PROCESS_RATES(1:nkn,19, 2)
    R_ZOO_RESP                 = PROCESS_RATES(1:nkn,9 , 3)
    R_ABIOTIC_DOC_MIN          = PROCESS_RATES(1:nkn,15, 5)

    TOTAL_DIC_KINETIC_SOURCES = &
        R_CHEM_AUT_BAC_TOT_RESP   + R_AER_HET_BAC_TOT_RESP + &
        R_FAC_AN_HET_BAC_TOT_RESP + R_DIA_TOT_RESP         + &
        R_CYN_TOT_RESP            + R_FIX_CYN_TOT_RESP     + &
        R_OPA_TOT_RESP            + R_ZOO_RESP             + &
        R_ABIOTIC_DOC_MIN

    ! Sinks
    R_CHEM_AUT_BAC_GROWTH      = PROCESS_RATES(1:nkn,5 , 1)
    R_DIA_GROWTH               = PROCESS_RATES(1:nkn,8 , 1)
    R_CYN_GROWTH               = PROCESS_RATES(1:nkn,18, 1)
    R_FIX_CYN_GROWTH           = PROCESS_RATES(1:nkn,22, 1)
    R_OPA_GROWTH               = PROCESS_RATES(1:nkn,19, 1)

    TOTAL_DIC_KINETIC_SINKS = &
       R_CHEM_AUT_BAC_GROWTH + R_DIA_GROWTH + R_CYN_GROWTH + &
       R_FIX_CYN_GROWTH      + R_OPA_GROWTH

    ! Atmospheric exchange
    T_A = TEMP + 273.6D0

    ! Calculate the saturaion concentration of CO2 in water
    P_K_H   = -(2385.73D0 / T_A) - (0.0152642D0 * T_A) + 14.0184D0
    K_H     = 10.0D0 ** (-P_K_H)
    P_CO2   = 10.0D0 ** (-3.416D0)
    CO2_SAT = K_H * P_CO2         !In moles

    ! Calculate the rearation rate constant for CO2
    K_A_CALC_CO2 = K_A_CALC * 0.923D0

    ! Finally calculate the atmospheric exchange rate
    CO2_ATM_EXHANGE = K_A_CALC_CO2 * (CO2_SAT - (ALPHA_0/1.0D6))

    !----------------------------------------------------------------------
    ! 2 February 2015 
    ! New code added to account the effect of ice cover.
    !----------------------------------------------------------------------
    CO2_ATM_EXHANGE = CO2_ATM_EXHANGE * (1.0D0 - ice_cover)
    !----------------------------------------------------------------------
    ! End of new code added to account the effect of ice cover.
    !----------------------------------------------------------------------

    ! End of atmospheric exchange

    ! Calculate the total derivative and convert it to moles
    if (CONSIDER_INORG_C_DERIVATIVE > 0) then
        PROCESS_RATES(1:nkn,23, 1) = TOTAL_DIC_KINETIC_SOURCES /12000.0D0
        PROCESS_RATES(1:nkn,23, 2) = TOTAL_DIC_KINETIC_SINKS  /12000.0D0

        if (CONSIDER_CO2_REARATION > 0) then
            PROCESS_RATES(1:nkn,23, 3) = CO2_ATM_EXHANGE
        else
            PROCESS_RATES(1:nkn,23, 3) = 0.0D0
            CO2_ATM_EXHANGE = 0.0D0
        end if

        PROCESS_RATES(1:nkn,23, 4) = R_CHEM_AUT_BAC_TOT_RESP
        PROCESS_RATES(1:nkn,23, 5) = R_AER_HET_BAC_TOT_RESP
        PROCESS_RATES(1:nkn,23, 6) = R_FAC_AN_HET_BAC_TOT_RESP
        PROCESS_RATES(1:nkn,23, 7) = R_DIA_TOT_RESP
        PROCESS_RATES(1:nkn,23, 8) = R_CYN_TOT_RESP
        PROCESS_RATES(1:nkn,23, 9) = R_FIX_CYN_TOT_RESP
        PROCESS_RATES(1:nkn,23,10) = R_OPA_TOT_RESP
        PROCESS_RATES(1:nkn,23,11) = R_ZOO_RESP
        PROCESS_RATES(1:nkn,23,12) = R_ABIOTIC_DOC_MIN
        PROCESS_RATES(1:nkn,23,13) = R_CHEM_AUT_BAC_GROWTH
        PROCESS_RATES(1:nkn,23,14) = R_DIA_GROWTH
        PROCESS_RATES(1:nkn,23,15) = R_CYN_GROWTH
        PROCESS_RATES(1:nkn,23,16) = R_FIX_CYN_GROWTH
        PROCESS_RATES(1:nkn,23,17) = R_OPA_GROWTH

        DIC_KINETIC_DERIVATIVE = CO2_ATM_EXHANGE + &
            ((TOTAL_DIC_KINETIC_SOURCES - TOTAL_DIC_KINETIC_SINKS) / 12000.0D0)

        DERIVATIVES(1:nkn,23) = DIC_KINETIC_DERIVATIVE
    else
        PROCESS_RATES(1:nkn,23, 1) = 0.0D0
        PROCESS_RATES(1:nkn,23, 2) = 0.0D0

        if (CONSIDER_CO2_REARATION > 0) then
            PROCESS_RATES(1:nkn,23, 3) = CO2_ATM_EXHANGE
        else
            PROCESS_RATES(1:nkn,23, 3) = 0.0D0
            CO2_ATM_EXHANGE = 0.0D0
        end if

        PROCESS_RATES(1:nkn,23, 4:17) = 0.0D0

        if (CONSIDER_CO2_REARATION > 0) then
            DERIVATIVES(1:nkn,23) = CO2_ATM_EXHANGE
        else
            DERIVATIVES(1:nkn,23) = 0.0D0
        end if

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
    T_A      = TEMP + 2.7316D2
    pKH      = 9.018D-2 + (2.72992D3 / T_A)
    FRAC_NH3 = 1.0D0 / (1.0D0 + (10.0D0 ** (pKH - PH)))
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
    N_CHEM_AUT_BAC_TOT_RESP   = PROCESS_RATES(1:nkn,1, 1)  * FRAC_NH4
    N_AER_HET_BAC_INT_RESP    = PROCESS_RATES(1:nkn,1, 2)  * FRAC_NH4
    N_FAC_AN_HET_BAC_TOT_RESP = PROCESS_RATES(1:nkn,1, 3)  * FRAC_NH4
    N_DIA_TOT_RESP            = PROCESS_RATES(1:nkn,1, 4)  * FRAC_NH4
    N_CYN_TOT_RESP            = PROCESS_RATES(1:nkn,1, 5)  * FRAC_NH4
    N_OPA_TOT_RESP            = PROCESS_RATES(1:nkn,1, 6)  * FRAC_NH4
    N_FIX_CYN_TOT_RESP        = PROCESS_RATES(1:nkn,1, 7)  * FRAC_NH4
    N_AER_HET_BAC_N_OX        = PROCESS_RATES(1:nkn,1, 8)  * FRAC_NH4
    N_FAC_AN_HET_BAC_N_OX     = PROCESS_RATES(1:nkn,1, 9)  * FRAC_NH4
    N_ZOO_TOT_RESP            = PROCESS_RATES(1:nkn,1, 10) * FRAC_NH4
    N_ABIOTIC_DON_MIN         = PROCESS_RATES(1:nkn,1, 17) * FRAC_NH4

    ALK_GAINED_BY_AMMONIUM_GEN = &
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
    N_DENITRIFICATION     = PROCESS_RATES(1:nkn,2, 2)
    N_DIA_GROWTH          = PROCESS_RATES(1:nkn,2, 3)
    N_CYN_GROWTH          = PROCESS_RATES(1:nkn,2, 4)
    N_OPA_GROWTH          = PROCESS_RATES(1:nkn,2, 5)
    N_NON_FIX_CYN_GROWTH  = PROCESS_RATES(1:nkn,2, 6)
    N_AER_HET_BAC_GROWTH  = PROCESS_RATES(1:nkn,2, 7)

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
    N_CHEM_AUT_BAC_GROWTH = PROCESS_RATES(1:nkn,1, 11) * FRAC_NH4
    N_DIA_GROWTH          = PROCESS_RATES(1:nkn,1, 12) * FRAC_NH4
    N_CYN_GROWTH          = PROCESS_RATES(1:nkn,1, 13) * FRAC_NH4
    N_OPA_GROWTH          = PROCESS_RATES(1:nkn,1, 14) * FRAC_NH4
    N_NON_FIX_CYN_GROWTH  = PROCESS_RATES(1:nkn,1, 15) * FRAC_NH4
    N_AER_HET_BAC_GROWTH  = PROCESS_RATES(1:nkn,1, 18) * FRAC_NH4

    ALK_LOST_BY_AMMONIUM_CONS = &
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
    N_NITRIFICATION_NH4 = PROCESS_RATES(1:nkn,1,16) * FRAC_NH4
    N_NITRIFICATION_NH3 = PROCESS_RATES(1:nkn,1,16) * FRAC_NH3

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


    KP_OPTION = 1

    if (any(SALT < 0.0d0)) then
        write(*,*) '!!!SALT!!!   !!!SALT!!!    !!!SALT!!!'
        write(*,*) 'SALINITY IS NEGATIVE'
        write(*,*) '!!!SALT!!!   !!!SALT!!!    !!!SALT!!!'
    endif

    select case (KP_OPTION)

        ! Option 0 : Use fixed values from the general chemistry book
        case (0)
            K_ONE_TIP   = 10.0D0 ** (-2.15D0)
            K_TWO_TIP   = 10.0D0 ** (-7.20D0)
            K_THREE_TIP = 10.0D0 ** (-2.15D0)

        ! Option 1 : Use values calculated by CO2SYS
        case (1)
            continue

        ! Option 2 : Use values from the chemical ocenography book, the unsafe formula
        case (2)
            A_1 = 115.525D0
            A_2 = -4576.7525D0
            A_3 = -18.453D0
            B_1 = 0.69171D0
            B_2 = -106.736D0
            C_1 = -0.01844D0
            C_2 = -0.6543D0

            K_ONE_TIP   = exp(A_1 + (A_2 / TEMP) + (A_3 * log(TEMP))   + &
                              ((B_1 + (B_2 / TEMP)) * (SALT ** 0.5D0)) + &
                              ((C_1 + (C_2 / TEMP)) * SALT))

            A_1 = 172.0883D0
            A_2 = -8814.715D0
            A_3 = -27.927D0
            B_1 = 1.3566D0
            B_2 = -160.34D0
            C_1 = -0.05778D0
            C_2 = -0.37335D0

            K_TWO_TIP   = exp(A_1 + (A_2 / TEMP) + (A_3 * log(TEMP))   + &
                              ((B_1 + (B_2 / TEMP)) * (SALT ** 0.5D0)) + &
                              ((C_1 + (C_2 / TEMP)) * SALT))

            A_1 = -18.141D0
            A_2 = -3070.75D0
            A_3 = 0.0D0
            B_1 = 2.81197D0
            B_2 = 17.27039D0
            C_1 = -0.09984D0
            C_2 = -44.99486D0

            K_THREE_TIP = exp(A_1 + (A_2 / TEMP) + (A_3 * log(TEMP))   + &
                              ((B_1 + (B_2 / TEMP)) * (SALT ** 0.5D0)) + &
                              ((C_1 + (C_2 / TEMP)) * SALT))
    end select

    FRACTION_DIVISOR_TIP = &
        (H_PLUS * H_PLUS * H_PLUS) + (K_ONE_TIP * K_ONE_TIP * H_PLUS) + &
        (K_ONE_TIP * K_TWO_TIP * H_PLUS) + (K_ONE_TIP * K_TWO_TIP * K_THREE_TIP)

    ALPHA_H2PO4 = (K_ONE_TIP * H_PLUS    * H_PLUS)      / FRACTION_DIVISOR_TIP
    ALPHA_HPO4  = (K_ONE_TIP * K_TWO_TIP * H_PLUS)      / FRACTION_DIVISOR_TIP
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
    P_CHEM_AUT_BAC_GROWTH = PROCESS_RATES(1:nkn,3, 11) * PHOSPHATE_EQ_CONSTANT
    P_DIA_GROWTH          = PROCESS_RATES(1:nkn,3, 12) * PHOSPHATE_EQ_CONSTANT
    P_CYN_GROWTH          = PROCESS_RATES(1:nkn,3, 13) * PHOSPHATE_EQ_CONSTANT
    P_OPA_GROWTH          = PROCESS_RATES(1:nkn,3, 14) * PHOSPHATE_EQ_CONSTANT
    P_NON_FIX_CYN_GROWTH  = PROCESS_RATES(1:nkn,3, 15) * PHOSPHATE_EQ_CONSTANT
    P_AER_HET_BAC_GROWTH  = PROCESS_RATES(1:nkn,3, 17) * PHOSPHATE_EQ_CONSTANT

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
    P_CHEM_AUT_BAC_TOT_RESP   = PROCESS_RATES(1:nkn,3, 1)  * PHOSPHATE_EQ_CONSTANT
    P_AER_HET_BAC_INT_RESP    = PROCESS_RATES(1:nkn,3, 2)  * PHOSPHATE_EQ_CONSTANT
    P_FAC_AN_HET_BAC_TOT_RESP = PROCESS_RATES(1:nkn,3, 3)  * PHOSPHATE_EQ_CONSTANT
    P_DIA_TOT_RESP            = PROCESS_RATES(1:nkn,3, 4)  * PHOSPHATE_EQ_CONSTANT
    P_CYN_TOT_RESP            = PROCESS_RATES(1:nkn,3, 5)  * PHOSPHATE_EQ_CONSTANT
    P_OPA_TOT_RESP            = PROCESS_RATES(1:nkn,3, 6)  * PHOSPHATE_EQ_CONSTANT
    P_FIX_CYN_TOT_RESP        = PROCESS_RATES(1:nkn,3, 7)  * PHOSPHATE_EQ_CONSTANT
    P_AER_HET_BAC_P_OX        = PROCESS_RATES(1:nkn,3, 8)  * PHOSPHATE_EQ_CONSTANT
    P_FAC_AN_HET_BAC_P_OX     = PROCESS_RATES(1:nkn,3, 9)  * PHOSPHATE_EQ_CONSTANT
    P_ZOO_TOT_RESP            = PROCESS_RATES(1:nkn,3, 10) * PHOSPHATE_EQ_CONSTANT
    P_ABIOTIC_DON_MIN         = PROCESS_RATES(1:nkn,3, 17) * PHOSPHATE_EQ_CONSTANT

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

    DERIVATIVES(1:nkn,24) = ALK_KINETIC_DERIVATIVE

    PROCESS_RATES(1:nkn,24, 1) = pH(1:nkn)
    PROCESS_RATES(1:nkn,24, 2) = ALK_GAINED_BY_AMMONIUM_GEN
    PROCESS_RATES(1:nkn,24, 3) = ALK_GAINED_BY_NITRATE_CONS
    PROCESS_RATES(1:nkn,24, 4) = ALK_GAINED_BY_PHOSPHATE_CONS
    PROCESS_RATES(1:nkn,24, 5) = ALK_LOST_BY_AMMONIUM_CONS
    PROCESS_RATES(1:nkn,24, 6) = ALK_LOST_BY_NITRIFICATION
    PROCESS_RATES(1:nkn,24, 7) = ALK_LOST_BY_PHOSPHATE_GEN

    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,24),VALUE_strange,nkn).eq.1) then
            nstrange = count(VALUE_strange)
            allocate(STRANGERS    (nstrange))
            allocate(NODES_STRANGE(nstrange))
            allocate(NODES_STRANGE_int(nstrange))
            allocate(NODES_STRANGE_ext(nstrange))
            j=1

            do k=1,nkn
                if(VALUE_strange(k)) then
                    STRANGERS    (j) = DERIVATIVES(k,24)
                    NODES_STRANGE(j) = k
                    j=j+1
                end if
            end do

            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            print *, 'PELAGIC_KINETICS:  Derivative 24'
            print *, 'is NaN or Inf:'
            print *, 'NODE_NUMBERS=',NODES_STRANGE
            print *, 'VALUES=',STRANGERS

            print *, '----------------------------------'
            print *, 'Related variables:'
            print *,'pH                          :',(pH                           &
                                                    (NODES_STRANGE(j)),j=1,nstrange)
            print *,'ALK_GAINED_BY_AMMONIUM_GEN  :',(ALK_GAINED_BY_AMMONIUM_GEN   &
                                                    (NODES_STRANGE(j)),j=1,nstrange)
            print *,'ALK_GAINED_BY_NITRATE_CONS  :',(ALK_GAINED_BY_NITRATE_CONS   &
                                                    (NODES_STRANGE(j)),j=1,nstrange)
            print *,'ALK_GAINED_BY_PHOSPHATE_CONS:',(ALK_GAINED_BY_PHOSPHATE_CONS &
                                                    (NODES_STRANGE(j)),j=1,nstrange)
            print *,'ALK_LOST_BY_AMMONIUM_CONS   :',(ALK_LOST_BY_AMMONIUM_CONS    &
                                                    (NODES_STRANGE(j)),j=1,nstrange)
            print *,'ALK_LOST_BY_NITRIFICATION   :',(ALK_LOST_BY_NITRIFICATION    &
                                                    (NODES_STRANGE(j)),j=1,nstrange)
            print *,'ALK_LOST_BY_PHOSPHATE_GEN   :',(ALK_LOST_BY_PHOSPHATE_GEN    &
                                                    (NODES_STRANGE(j)),j=1,nstrange)

            error=1
            deallocate(STRANGERS    )
            deallocate(NODES_STRANGE)
            deallocate(NODES_STRANGE_int)
            deallocate(NODES_STRANGE_ext)
        end if

        if(error .eq. 1) stop
    end if

    ! Checking all derivatives
    if(debug_stranger) then
        do i = 1, nstate
            if (STRANGERSD(DERIVATIVES(1:nkn,i),VALUE_strange,nkn).eq.1) then
                nstrange = count(VALUE_strange)
                allocate(STRANGERS    (nstrange))
                allocate(NODES_STRANGE(nstrange))

                j=1

                do k=1,nkn
                    if(VALUE_strange(k)) then
                        STRANGERS    (j) = DERIVATIVES(k,i)
                        NODES_STRANGE(j) = k
                        j=j+1
                    end if
                end do

                print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                print *, 'PELAGIC_KINETICS:  Derivative ',i !'Cell ',k
                print *, 'is NaN or Inf:'
                print *, 'NODE_NUMBERS=',NODES_STRANGE
                print *, 'VALUES=',STRANGERS
                error=1
                deallocate(STRANGERS)
                deallocate(NODES_STRANGE)
            end if
        end do

        if(error .eq. 1) stop
    end if
    !*********************************************'
    !*                                           *'
    !*          END OF ECOLOGY KINETICS          *'
    !*                                           *'
    !*      DO NOT CHANGE THE FOLLOWING CODE     *'
    !*                                           *'
    !*********************************************'   
    
end subroutine PELAGIC_KINETICS


!*********************************************************
!*********************************************************
!*********************************************************

subroutine derived_vars(nkn,pH,STATE_VARIABLES, nstate, &
                          MODEL_CONSTANTS, nconst, &
                          WC_OUTPUTS, noutput)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine for calculation of derived variables (calculated from state variables
! and model or other internally used constants and derived variables
!
! Array WC_OUTPUTS is used for the output and gets (done in aquabc) also values of
! state variables as first nstate variables. Therefore is assumed that derived
! variables are placed after state variables.
!
! While this routine is called in pelagic model before integration,
! derived variables are calculated using input for the current step
! and will correspond to the time that is one time step earlier
! than the time after integration
!  Note:
!    MODEL_CONSTANTS are not necessary more but still are left
!    on parameters list. Fixme
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use AQUABC_II_GLOBAL
  use para_aqua
  
  implicit none
  include 'param.h'

  double precision, dimension(nkn,nstate), intent(in)     :: STATE_VARIABLES
  integer  :: nkn, nstate, nconst, noutput, ext_node
  double precision, dimension(nconst),      intent(in)    :: MODEL_CONSTANTS
  double precision, dimension(nkn,noutput), intent(inout) :: WC_OUTPUTS

! State variables
    double precision :: NH4_N             (nkn)
    double precision :: NO3_N             (nkn)
    double precision :: PO4_P             (nkn)
    double precision :: DISS_OXYGEN       (nkn)
    double precision :: CHEM_AUT_BAC_C    (nkn)
    double precision :: AER_HET_BAC_C     (nkn)
    double precision :: FAC_AN_HET_BAC_C  (nkn)
    double precision :: DIA_C             (nkn)
    double precision :: ZOO_C             (nkn)
    double precision :: ZOO_N             (nkn)
    double precision :: ZOO_P             (nkn)
    double precision :: DET_PART_ORG_C    (nkn)
    double precision :: DET_PART_ORG_N    (nkn)
    double precision :: DET_PART_ORG_P    (nkn)
    double precision :: DISS_ORG_C        (nkn)
    double precision :: DISS_ORG_N        (nkn)
    double precision :: DISS_ORG_P        (nkn)
    double precision :: CYN_C             (nkn)
    double precision :: OPA_C             (nkn)
    double precision :: DISS_Si           (nkn)
    double precision :: PART_Si           (nkn)
    double precision :: FIX_CYN_C         (nkn)
    double precision :: pH                (nkn)
    ! New state variables added 22 September 2014
    real(kind=DBL_PREC) :: INORG_C        (nkn)   !Inorganic carbon
    real(kind=DBL_PREC) :: TOT_ALK        (nkn)   !Total alkalinity

! Constants
 double precision ::                              K_A
 double precision ::                        THETA_K_A
 double precision ::               KG_CHEM_AUT_BAC_20
 double precision ::          EFF_CHEM_AUT_BAC_GROWTH
 double precision ::            THETA_KG_CHEM_AUT_BAC
 double precision ::               KR_CHEM_AUT_BAC_20
 double precision ::            THETA_KR_CHEM_AUT_BAC
 double precision ::               KD_CHEM_AUT_BAC_20
 double precision ::            THETA_KD_CHEM_AUT_BAC
 double precision ::            KHS_NH4N_CHEM_AUT_BAC
 double precision ::            KHS_PO4P_CHEM_AUT_BAC
 double precision ::              KHS_O2_CHEM_AUT_BAC
 double precision ::      DO_STR_HYPOX_CHEM_AUT_BAC_D
 double precision ::       THETA_HYPOX_CHEM_AUT_BAC_D
 double precision ::       EXPON_HYPOX_CHEM_AUT_BAC_D
 double precision ::              CHEM_AUT_BAC_N_TO_C
 double precision ::              CHEM_AUT_BAC_P_TO_C
 double precision ::             CHEM_AUT_BAC_O2_TO_C
 double precision ::               YIELD_CHEM_AUT_BAC
 double precision ::                KG_AER_HET_BAC_20
 double precision ::           EFF_AER_HET_BAC_GROWTH
 double precision ::             THETA_KG_AER_HET_BAC
 double precision ::                KR_AER_HET_BAC_20
 double precision ::             THETA_KR_AER_HET_BAC
 double precision ::                KD_AER_HET_BAC_20
 double precision ::             THETA_KD_AER_HET_BAC
 double precision ::             KHS_ORGC_AER_HET_BAC
 double precision ::             KHS_ORGN_AER_HET_BAC
 double precision ::             KHS_ORGP_AER_HET_BAC
 double precision ::               KHS_O2_AER_HET_BAC
 double precision ::              KHS_DIN_AER_HET_BAC
 double precision ::              KHS_DIP_AER_HET_BAC
 double precision ::              KHS_PHYT_AER_HET_BAC
 double precision ::            YIELD_OC_AER_HET_BAC
 double precision ::             OX_ORGN_AER_HET_BAC
 double precision ::                       KHS_MIN_N
 double precision ::             OX_ORGP_AER_HET_BAC
 double precision ::                       KHS_MIN_P
 double precision ::       DO_STR_HYPOX_AER_HET_BAC_D
 double precision ::        THETA_HYPOX_AER_HET_BAC_D
 double precision ::        EXPON_HYPOX_AER_HET_BAC_D
 double precision ::               AER_HET_BAC_N_TO_C
 double precision ::               AER_HET_BAC_P_TO_C
 double precision ::              AER_HET_BAC_O2_TO_C
 double precision ::             KG_FAC_AN_HET_BAC_20
 double precision ::        EFF_FAC_AN_HET_BAC_GROWTH
 double precision ::          THETA_KG_FAC_AN_HET_BAC
 double precision ::             KR_FAC_AN_HET_BAC_20
 double precision ::          THETA_KR_FAC_AN_HET_BAC
 double precision ::             KD_FAC_AN_HET_BAC_20
 double precision ::          THETA_KD_FAC_AN_HET_BAC
 double precision ::          KHS_NO3N_FAC_AN_HET_BAC
 double precision ::          KHS_ORGC_FAC_AN_HET_BAC
 double precision ::          KHS_ORGN_FAC_AN_HET_BAC
 double precision ::          KHS_ORGP_FAC_AN_HET_BAC
 double precision ::        REV_KHS_O2_FAC_AN_HET_BAC
 double precision ::   NO3N_LACK_STR_FAC_AN_HET_BAC_D
 double precision ::  THETA_NO3_LACK_FAC_AN_HET_BAC_D
 double precision ::    EXP_NO3_LACK_FAC_AN_HET_BAC_D
 double precision ::            FAC_AN_HET_BAC_N_TO_C
 double precision ::            FAC_AN_HET_BAC_P_TO_C
 double precision ::           FAC_AN_HET_BAC_O2_TO_C
 double precision ::             YIELD_FAC_AN_HET_BAC
 double precision ::                  KG_DIA_OPT_TEMP
 double precision ::                  DIA_OPT_TEMP_LR
 double precision ::                  DIA_OPT_TEMP_UR
 double precision ::                   EFF_DIA_GROWTH
 double precision ::         KAPPA_DIA_UNDER_OPT_TEMP
 double precision ::          KAPPA_DIA_OVER_OPT_TEMP
 double precision ::                        KR_DIA_20
 double precision ::                     THETA_KR_DIA
 double precision ::                        KD_DIA_20
 double precision ::                     THETA_KD_DIA
 double precision ::                      KHS_DIN_DIA
 double precision ::                      KHS_DIP_DIA
 double precision ::                      KHS_DSi_DIA
 double precision ::                       KHS_O2_DIA
 double precision ::                    FRAC_DIA_EXCR
 double precision ::                          I_S_DIA
 double precision ::               DO_STR_HYPOX_DIA_D
 double precision ::                THETA_HYPOX_DIA_D
 double precision ::                EXPON_HYPOX_DIA_D
 double precision ::                       DIA_N_TO_C
 double precision ::                       DIA_P_TO_C
 double precision ::                      DIA_Si_TO_C
 double precision ::                      DIA_O2_TO_C
 double precision ::                    DIA_C_TO_CHLA
 double precision ::                  KG_CYN_OPT_TEMP
 double precision ::                  CYN_OPT_TEMP_LR
 double precision ::                  CYN_OPT_TEMP_UR
 double precision ::                   EFF_CYN_GROWTH
 double precision ::         KAPPA_CYN_UNDER_OPT_TEMP
 double precision ::          KAPPA_CYN_OVER_OPT_TEMP
 double precision ::                        KR_CYN_20
 double precision ::                     THETA_KR_CYN
 double precision ::                        KD_CYN_20
 double precision ::                     THETA_KD_CYN
 double precision ::                      KHS_DIN_CYN
 double precision ::                      KHS_DIP_CYN
 double precision ::                       KHS_O2_CYN
 double precision ::                    FRAC_CYN_EXCR
 double precision ::                          I_S_CYN
 double precision ::               DO_STR_HYPOX_CYN_D
 double precision ::                THETA_HYPOX_CYN_D
 double precision ::                EXPON_HYPOX_CYN_D
 double precision ::                       CYN_N_TO_C
 double precision ::                       CYN_P_TO_C
 double precision ::                      CYN_O2_TO_C
 double precision ::                    CYN_C_TO_CHLA
 double precision ::              KG_FIX_CYN_OPT_TEMP
 double precision ::              FIX_CYN_OPT_TEMP_LR
 double precision ::              FIX_CYN_OPT_TEMP_UR
 double precision ::               EFF_FIX_CYN_GROWTH
 double precision ::     KAPPA_FIX_CYN_UNDER_OPT_TEMP
 double precision ::      KAPPA_FIX_CYN_OVER_OPT_TEMP
 double precision ::                    KR_FIX_CYN_20
 double precision ::                 THETA_KR_FIX_CYN
 double precision ::                    KD_FIX_CYN_20
 double precision ::                 THETA_KD_FIX_CYN
 double precision ::                  KHS_DIN_FIX_CYN
 double precision ::                  KHS_DIP_FIX_CYN
 double precision ::                   KHS_O2_FIX_CYN
 double precision ::                FRAC_FIX_CYN_EXCR
 double precision ::                      I_S_FIX_CYN
 double precision ::           DO_STR_HYPOX_FIX_CYN_D
 double precision ::            THETA_HYPOX_FIX_CYN_D
 double precision ::            EXPON_HYPOX_FIX_CYN_D
 double precision ::                   FIX_CYN_N_TO_C
 double precision ::                   FIX_CYN_P_TO_C
 double precision ::                  FIX_CYN_O2_TO_C
 double precision ::                FIX_CYN_C_TO_CHLA
 double precision ::                            R_FIX
 double precision ::                            K_FIX
 double precision ::                  KG_OPA_OPT_TEMP
 double precision ::                  OPA_OPT_TEMP_LR
 double precision ::                  OPA_OPT_TEMP_UR
 double precision ::                   EFF_OPA_GROWTH
 double precision ::         KAPPA_OPA_UNDER_OPT_TEMP
 double precision ::          KAPPA_OPA_OVER_OPT_TEMP
 double precision ::                        KR_OPA_20
 double precision ::                     THETA_KR_OPA
 double precision ::                        KD_OPA_20
 double precision ::                     THETA_KD_OPA
 double precision ::                      KHS_DIN_OPA
 double precision ::                      KHS_DIP_OPA
 double precision ::                       KHS_O2_OPA
 double precision ::                    FRAC_OPA_EXCR
 double precision ::                          I_S_OPA
 double precision ::               DO_STR_HYPOX_OPA_D
 double precision ::                THETA_HYPOX_OPA_D
 double precision ::                EXPON_HYPOX_OPA_D
 double precision ::                       OPA_N_TO_C
 double precision ::                       OPA_P_TO_C
 double precision ::                      OPA_O2_TO_C
 double precision ::                    OPA_C_TO_CHLA
 double precision ::                  KG_ZOO_OPT_TEMP
 double precision ::                  ZOO_OPT_TEMP_LR
 double precision ::                  ZOO_OPT_TEMP_UR
 double precision ::                   EFF_ZOO_GROWTH
 double precision ::         KAPPA_ZOO_UNDER_OPT_TEMP
 double precision ::          KAPPA_ZOO_OVER_OPT_TEMP
 double precision ::                     GRAT_ZOO_DIA
 double precision ::                     GRAT_ZOO_CYN
 double precision ::                     GRAT_ZOO_OPA
 double precision ::                 GRAT_ZOO_FIX_CYN
 double precision ::            GRAT_ZOO_CHEM_AUT_BAC
 double precision ::             GRAT_ZOO_AER_HET_BAC
 double precision ::          GRAT_ZOO_FAC_AN_HET_BAC
 double precision ::          GRAT_ZOO_DET_PART_ORG_C
 double precision ::                     PREF_ZOO_DIA
 double precision ::                     PREF_ZOO_CYN
 double precision ::                 PREF_ZOO_FIX_CYN
 double precision ::                     PREF_ZOO_OPA
 double precision ::            PREF_ZOO_CHEM_AUT_BAC
 double precision ::             PREF_ZOO_AER_HET_BAC
 double precision ::          PREF_ZOO_FAC_AN_HET_BAC
 double precision ::          PREF_ZOO_DET_PART_ORG_C
 double precision ::                    KHS_DIA_C_ZOO
 double precision ::                    KHS_CYN_C_ZOO
 double precision ::                KHS_FIX_CYN_C_ZOO
 double precision ::                    KHS_OPA_C_ZOO
 double precision ::           KHS_CHEM_AUT_BAC_C_ZOO
 double precision ::            KHS_AER_HET_BAC_C_ZOO
 double precision ::         KHS_FAC_AN_HET_BAC_C_ZOO
 double precision ::           KHS_DET_PART_ORG_C_ZOO
 double precision ::                     FOOD_MIN_ZOO
 double precision ::                           KE_ZOO
 double precision ::                  FRAC_ZOO_EX_ORG
 double precision ::                        KR_ZOO_20
 double precision ::                     THETA_KR_ZOO
 double precision ::                        KD_ZOO_20
 double precision ::                     THETA_KD_ZOO
 double precision ::               DO_STR_HYPOX_ZOO_D
 double precision ::                THETA_HYPOX_ZOO_D
 double precision ::                EXPON_HYPOX_ZOO_D
 double precision ::                       ZOO_N_TO_C
 double precision ::                       ZOO_P_TO_C
 double precision ::                      ZOO_O2_TO_C
 double precision ::          KDISS_DET_PART_ORG_C_20
 double precision ::       THETA_KDISS_DET_PART_ORG_C
 double precision ::          FAC_PHYT_DET_PART_ORG_C
 double precision ::          KDISS_DET_PART_ORG_N_20
 double precision ::       THETA_KDISS_DET_PART_ORG_N
 double precision ::                       KHS_DISS_N
 double precision ::          FAC_PHYT_DET_PART_ORG_N
 double precision ::          KDISS_DET_PART_ORG_P_20
 double precision ::       THETA_KDISS_DET_PART_ORG_P
 double precision ::                       KHS_DISS_P
 double precision ::          FAC_PHYT_DET_PART_ORG_P
 double precision ::                 KDISS_PART_Si_20
 double precision ::              THETA_KDISS_PART_Si
 double precision ::                     K_MIN_DOC_20
 double precision ::                  THETA_K_MIN_DOC
 double precision ::                FAC_PHYT_AMIN_DOC
 double precision ::                     K_MIN_DON_20
 double precision ::                  THETA_K_MIN_DON
 double precision ::                       KHS_AMIN_N
 double precision ::                FAC_PHYT_AMIN_DON
 double precision ::                     K_MIN_DOP_20
 double precision ::                  THETA_K_MIN_DOP
 double precision ::                       KHS_AMIN_P
 double precision ::                FAC_PHYT_AMIN_DOP
 double precision ::                        K_NITR_20
 double precision ::                     KHS_NITR_OXY
 double precision ::                   KHS_NITR_NH4_N
 double precision ::                     THETA_K_NITR


    NH4_N           (:)      = STATE_VARIABLES(:,1)      ! AMMONIUM NITROGEN
    NO3_N           (:)      = STATE_VARIABLES(:,2)      ! NITRATE NITROGEN
    PO4_P           (:)      = STATE_VARIABLES(:,3)      ! ORTHOPHOSPHATE PHOSPHORUS
    DISS_OXYGEN     (:)      = STATE_VARIABLES(:,4)      ! DISSOLVED OXYGEN
    CHEM_AUT_BAC_C  (:)      = STATE_VARIABLES(:,5)      ! CHEMOAUTIC BACTERIA CARBON
    AER_HET_BAC_C   (:)      = STATE_VARIABLES(:,6)      ! AEROBIC HETEROTROPHIC BACTERIA CARBON
    FAC_AN_HET_BAC_C(:)      = STATE_VARIABLES(:,7)      ! FACULTATIVE ANAEROBIC HETEROTROPHIC BACTERIA CARBON
    DIA_C           (:)      = STATE_VARIABLES(:,8)      ! DIATOMS CARBON
    ZOO_C           (:)      = STATE_VARIABLES(:,9)      ! ZOOPLANKTON CARBON
    ZOO_N           (:)      = STATE_VARIABLES(:,10)     ! ZOOPLANKTON CARBON
    ZOO_P           (:)      = STATE_VARIABLES(:,11)     ! ZOOPLANKTON PHOSPHORUS
    DET_PART_ORG_C  (:)      = STATE_VARIABLES(:,12)     ! DETRITUS PARTICULATE ORG. CARBON
    DET_PART_ORG_N  (:)      = STATE_VARIABLES(:,13)     ! DETRITUS PARTICULATE ORG. NITROGEN
    DET_PART_ORG_P  (:)      = STATE_VARIABLES(:,14)     ! DETRITUS PARTICULATE ORG. PHOSPHORUS
    DISS_ORG_C      (:)      = STATE_VARIABLES(:,15)     ! DISSOLVED ORGANIC CARBON
    DISS_ORG_N      (:)      = STATE_VARIABLES(:,16)     ! DISSOLVED ORGANIC NITROGEN
    DISS_ORG_P      (:)      = STATE_VARIABLES(:,17)     ! DISSOLVED ORGANIC PHOSPHORUS
    CYN_C           (:)      = STATE_VARIABLES(:,18)     ! NON FIXING CYANOBACTERIA CARBON
    OPA_C           (:)      = STATE_VARIABLES(:,19)     ! OTHER PHYTOPLANKTON CARBON
    DISS_Si         (:)      = STATE_VARIABLES(:,20)     ! DISSOLOVED SILICA
    PART_Si         (:)      = STATE_VARIABLES(:,21)     ! PARTICULATE SILICA
    FIX_CYN_C       (:)      = STATE_VARIABLES(:,22)     ! FIXING CYANOBACTERIA CARBON
    INORG_C         (:)      = STATE_VARIABLES(:,23)     ! INORG CARBON CARBON
    TOT_ALK         (:)      = STATE_VARIABLES(:,24)     ! TOTAL ALKALNITY


!                              K_A   =   MODEL_CONSTANTS(  1) !  1! Aeration coefficient (if negative calculates internally)
!                        THETA_K_A   =   MODEL_CONSTANTS(  2) !  2! Temperature correction factor for aeration
!               KG_CHEM_AUT_BAC_20   =   MODEL_CONSTANTS(  3) !  3! Chemoautotrophic bacteria Growth rate
!          EFF_CHEM_AUT_BAC_GROWTH   =   MODEL_CONSTANTS(  4) !  4! Chemoautotrophic bacteria growth efficiency
!            THETA_KG_CHEM_AUT_BAC   =   MODEL_CONSTANTS(  5) !  5! Chemoautotrophic bacteria Temperature correction for growth rate
!               KR_CHEM_AUT_BAC_20   =   MODEL_CONSTANTS(  6) !  6! Chemoautotrophic bacteria Respiration rate
!            THETA_KR_CHEM_AUT_BAC   =   MODEL_CONSTANTS(  7) !  7! Chemoautotrophic bacteria Temperature correction for respiration rate
!               KD_CHEM_AUT_BAC_20   =   MODEL_CONSTANTS(  8) !  8! Chemoautotrophic bacteria Mortality rate
!            THETA_KD_CHEM_AUT_BAC   =   MODEL_CONSTANTS(  9) !  9! Chemoautotrophic bacteria Temperature correction for Mortality rate
!            KHS_NH4N_CHEM_AUT_BAC   =   MODEL_CONSTANTS( 10) ! 10! Chemoautotrophic bacteria Half saturation growth for NH4N
!            KHS_PO4P_CHEM_AUT_BAC   =   MODEL_CONSTANTS( 11) ! 11! Chemoautotrophic bacteria Half saturation growth for PO4P
!              KHS_O2_CHEM_AUT_BAC   =   MODEL_CONSTANTS( 12) ! 12! Chemoautotrophic bacteria Half saturation growth for O2
!      DO_STR_HYPOX_CHEM_AUT_BAC_D   =   MODEL_CONSTANTS( 13) ! 13! Chemoautotrophic bacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!       THETA_HYPOX_CHEM_AUT_BAC_D   =   MODEL_CONSTANTS( 14) ! 14! Chemoautotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
!       EXPON_HYPOX_CHEM_AUT_BAC_D   =   MODEL_CONSTANTS( 15) ! 15! Chemoautotrophic bacteria Exponent constant for Dissolved oxygen stress
!              CHEM_AUT_BAC_N_TO_C   =   MODEL_CONSTANTS( 16) ! 16! Chemoautotrophic bacteria Nitrogen to Carbon ratio
!              CHEM_AUT_BAC_P_TO_C   =   MODEL_CONSTANTS( 17) ! 17! Chemoautotrophic bacteria Phosphorus to Carbon ratio
!             CHEM_AUT_BAC_O2_TO_C   =   MODEL_CONSTANTS( 18) ! 18! Chemoautotrophic bacteria Oxygen to Carbon ratio
!               YIELD_CHEM_AUT_BAC   =   MODEL_CONSTANTS( 19) ! 19! Chemoautotrophic bacteria Yield of Carbon per unit amonia nitrogen
!                KG_AER_HET_BAC_20   =   MODEL_CONSTANTS( 20) ! 20! Aerobic heterotrophic bacteria Growth rate
!           EFF_AER_HET_BAC_GROWTH   =   MODEL_CONSTANTS( 21) ! 21! Aerobic heterotrophic bacteria growth efficiency
!             THETA_KG_AER_HET_BAC   =   MODEL_CONSTANTS( 22) ! 22! Aerobic heterotrophic bacteria Temperature correction for growth rate
!                KR_AER_HET_BAC_20   =   MODEL_CONSTANTS( 23) ! 23! Aerobic heterotrophic bacteria Respiration rate
!             THETA_KR_AER_HET_BAC   =   MODEL_CONSTANTS( 24) ! 24! Aerobic heterotrophic bacteria Temperature correction for respiration rate
!                KD_AER_HET_BAC_20   =   MODEL_CONSTANTS( 25) ! 25! Aerobic heterotrophic bacteria Mortality rate
!             THETA_KD_AER_HET_BAC   =   MODEL_CONSTANTS( 26) ! 26! Aerobic heterotrophic bacteria Temperature correction for Mortality rate
!             KHS_ORGC_AER_HET_BAC   =   MODEL_CONSTANTS( 27) ! 27! Aerobic heterotrophic bacteria Half saturation growth for OC
!
!             KHS_ORGN_AER_HET_BAC   =   MODEL_CONSTANTS( 28)  ! 28! Aerobic heterotrophic bacteria Half saturation growth for ON
!             KHS_ORGP_AER_HET_BAC   =   MODEL_CONSTANTS( 29)  ! 29! Aerobic heterotrophic bacteria Half saturation growth for OP
!               KHS_O2_AER_HET_BAC   =   MODEL_CONSTANTS( 30)  ! 30! Aerobic heterotrophic bacteria Half saturation growth for Oxygen
!              KHS_DIN_AER_HET_BAC   =   MODEL_CONSTANTS( 31)  ! 31! Aerobic heterotrophic bacteria Half saturation growth for inorganic nitrogen
!              KHS_DIP_AER_HET_BAC   =   MODEL_CONSTANTS( 32)  ! 32! Aerobic heterotrophic bacteria Half saturation growth for inorganic phosphorus
!              KHS_PHYT_AER_HET_BAC  =   MODEL_CONSTANTS( 33)  ! 33! Aerobic heterotrophic bacteria Half saturation growth for Phytoplankton C (not used as a resource)
!            YIELD_OC_AER_HET_BAC    =   MODEL_CONSTANTS( 34)  ! 34! Aerobic heterotrophic bacteria Yield of bacteria carbon per unit of organic carbon
!             OX_ORGN_AER_HET_BAC    =   MODEL_CONSTANTS( 35)  ! 35! Aerobic heterotrophic bacteria ON oxidation rate mg N per mg C of bacteria production
!
!                       KHS_MIN_N    =   MODEL_CONSTANTS(36 )  ! 36! Aerobic heterotrophic bacteria ON mineralisation reverse half saturation for DIN
!             OX_ORGP_AER_HET_BAC    =   MODEL_CONSTANTS(37 )  ! 37! Aerobic heterotrophic bacteria OP mineralisation rate mg P per mg C of bacteria production
!                       KHS_MIN_P    =   MODEL_CONSTANTS(38 )  ! 38! Aerobic heterotrophic bacteria ON mineralisation reverse half saturation for DIP
!       DO_STR_HYPOX_AER_HET_BAC_D   =   MODEL_CONSTANTS(39 )  ! 39! Aerobic heterotrophic bacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!        THETA_HYPOX_AER_HET_BAC_D   =   MODEL_CONSTANTS(40 )  ! 40! Aerobic heterotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
!        EXPON_HYPOX_AER_HET_BAC_D   =   MODEL_CONSTANTS(41 )  ! 41! Aerobic heterotrophic bacteria Exponent constant for Dissolved oxygen stress
!               AER_HET_BAC_N_TO_C   =   MODEL_CONSTANTS(42 )  ! 42! Aerobic heterotrophic bacteria Nitrogen to Carbon ratio
!               AER_HET_BAC_P_TO_C   =   MODEL_CONSTANTS(43 )  ! 43! Aerobic heterotrophic bacteria Phosphorus to Carbon ratio
!              AER_HET_BAC_O2_TO_C   =   MODEL_CONSTANTS(44 )  ! 44! Aerobic heterotrophic bacteria Oxygen to Carbon ratio for respiration
!             KG_FAC_AN_HET_BAC_20   =   MODEL_CONSTANTS(45 )  ! 45! Facultative anaerobic heterotrophic bacteria Growth rate of
!        EFF_FAC_AN_HET_BAC_GROWTH   =   MODEL_CONSTANTS(46 )  ! 46! not used! Facultative anaerobic heterotrophic bacteria growth efficiency
!          THETA_KG_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(47 )  ! 47! not used! Facultative anaerobic heterotrophic bacteria Temperature correction for growth rate
!             KR_FAC_AN_HET_BAC_20   =   MODEL_CONSTANTS(48 )  ! 48! not used! Facultative anaerobic heterotrophic bacteria Respiration rate
!          THETA_KR_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(49 )  ! 49! not used! Facultative anaerobic heterotrophic bacteria Temperature correction for respiration rate
!             KD_FAC_AN_HET_BAC_20   =   MODEL_CONSTANTS(50 )  ! 50! not used! Facultative anaerobic heterotrophic bacteria Mortality rate
!          THETA_KD_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(51 )  ! 51! not used! Facultative anaerobic heterotrophic bacteria Temperature correction for Mortality rate
!          KHS_NO3N_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(52 )  ! 52! Facultative anaerobic heterotrophic bacteria Half saturation growth for NO3N
!          KHS_ORGC_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(53 )  ! 53! not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for OC
!          KHS_ORGN_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(54 )  ! 54! not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for ON
!          KHS_ORGP_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(55 )  ! 55! not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for OP
!        REV_KHS_O2_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(56 )  ! 56! not used! Facultative anaerobic heterotrophic bacteria Reverse Half saturation growth for O2
!   NO3N_LACK_STR_FAC_AN_HET_BAC_D   =   MODEL_CONSTANTS(57 )  ! 57! not used! Facultative anaerobic heterotrophic bacteria NO3N stress concentration
!  THETA_NO3_LACK_FAC_AN_HET_BAC_D   =   MODEL_CONSTANTS(58 )  ! 58! not used! Facultative anaerobic heterotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
!    EXP_NO3_LACK_FAC_AN_HET_BAC_D   =   MODEL_CONSTANTS(59 )  ! 59! not used! Facultative anaerobic heterotrophic bacteria Exponent constant for Dissolved oxygen stress
!            FAC_AN_HET_BAC_N_TO_C   =   MODEL_CONSTANTS(60 )  ! 60! not used! Facultative anaerobic heterotrophic bacteria Nitrogen to Carbon ratio
!            FAC_AN_HET_BAC_P_TO_C   =   MODEL_CONSTANTS(61 )  ! 61! not used! Facultative anaerobic heterotrophic bacteria Phosphorus to Carbon ratio
!           FAC_AN_HET_BAC_O2_TO_C   =   MODEL_CONSTANTS(62 )  ! 62! not used! Facultative anaerobic heterotrophic bacteria Oxygen to Carbon ratio for respiration
!             YIELD_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(63 )  ! 63! Facultative anaerobic heterotrophic bacteria Yield of carbon per unit nitrates nitrogen
!                  KG_DIA_OPT_TEMP   =   MODEL_CONSTANTS(64 )  ! 64! Diatoms Growth rate
!                  DIA_OPT_TEMP_LR   =   MODEL_CONSTANTS(65 )  ! 65! Diatoms optimal temperature lower range
!                  DIA_OPT_TEMP_UR   =   MODEL_CONSTANTS(66 )  ! 66! Diatoms optimal temperature upper range
!                   EFF_DIA_GROWTH   =   MODEL_CONSTANTS(67 )  ! 67! Diatoms Effective growth. (1-EG)*growth - losses for respiration and excretion
!         KAPPA_DIA_UNDER_OPT_TEMP   =   MODEL_CONSTANTS(68 )  ! 68! Diatoms Temperature correction for growth lower temperature
!          KAPPA_DIA_OVER_OPT_TEMP   =   MODEL_CONSTANTS(69 )  ! 69! Diatoms Temperature correction for growth upper temperature
!                        KR_DIA_20   =   MODEL_CONSTANTS(70 )  ! 70! Diatoms Respiration rate
!                     THETA_KR_DIA   =   MODEL_CONSTANTS(71 )  ! 71! Diatoms Temperature correction for basal respiration rate
!                        KD_DIA_20   =   MODEL_CONSTANTS(72 )  ! 72! Diatoms Mortality rate
!                     THETA_KD_DIA   =   MODEL_CONSTANTS(73 )  ! 73! Diatoms Temperature correction for Mortality rate
!                      KHS_DIN_DIA   =   MODEL_CONSTANTS(74 )  ! 74! Diatoms Half saturation growth for DIN
!                      KHS_DIP_DIA   =   MODEL_CONSTANTS(75 )  ! 75! Diatoms Half saturation growth for DIP
!                      KHS_DSi_DIA   =   MODEL_CONSTANTS(76 )  ! 76! Diatoms Half saturation growth for DSi
!                       KHS_O2_DIA   =   MODEL_CONSTANTS(77 )  ! 77! Diatoms Half saturation growth for O2
!                    FRAC_DIA_EXCR   =   MODEL_CONSTANTS(78 )  ! 78! Diatoms Fraction of excretion in metabolism rate
!                          I_S_DIA   =   MODEL_CONSTANTS(79 )  ! 79! Diatoms Light saturation (langleys)
!               DO_STR_HYPOX_DIA_D   =   MODEL_CONSTANTS(80 )  ! 80! Diatoms Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!                THETA_HYPOX_DIA_D   =   MODEL_CONSTANTS(81 )  ! 81! Diatoms Multiplier of the exponent for Dissolved oxygen stress
!                EXPON_HYPOX_DIA_D   =   MODEL_CONSTANTS(82 )  ! 82! Diatoms Exponent constant for Dissolved oxygen stress
!                       DIA_N_TO_C   =   MODEL_CONSTANTS(83 )  ! 83! Diatoms Nitrogen to Carbon ratio
!                       DIA_P_TO_C   =   MODEL_CONSTANTS(84 )  ! 84! Diatoms Phosphorus to Carbon ratio
!                      DIA_Si_TO_C   =   MODEL_CONSTANTS(85 )  ! 85! Diatoms Silica to Carbon ratio
!                      DIA_O2_TO_C   =   MODEL_CONSTANTS(86 )  ! 86! Diatoms Oxygen to Carbon ratio for respiration
!                    DIA_C_TO_CHLA   =   MODEL_CONSTANTS(87 )  ! 87! Diatoms Carbon to Chlorophil a ratio
!                  KG_CYN_OPT_TEMP   =   MODEL_CONSTANTS(88 )  ! 88! Non-fixing cyanobacteria Growth rate
!                  CYN_OPT_TEMP_LR   =   MODEL_CONSTANTS(89 )  ! 89! Non-fixing cyanobacteria optimal temperature lower range
!                  CYN_OPT_TEMP_UR   =   MODEL_CONSTANTS(90 )  ! 90! Non-fixing cyanobacteria optimal temperature upper range
!                   EFF_CYN_GROWTH   =   MODEL_CONSTANTS(91 )  ! 91! Non-fixing cyanobacteria Effective growth. (1-EG)*growth - losses for respiration and excretion
!         KAPPA_CYN_UNDER_OPT_TEMP   =   MODEL_CONSTANTS(92 )  ! 92! Non-fixing cyanobacteria Temperature correction for growth lower temperature
!          KAPPA_CYN_OVER_OPT_TEMP   =   MODEL_CONSTANTS(93 )  ! 93! Non-fixing cyanobacteria Temperature correction for growth upper temperature
!                        KR_CYN_20   =   MODEL_CONSTANTS(94 )  ! 94! Non-fixing cyanobacteria Respiration rate
!                     THETA_KR_CYN   =   MODEL_CONSTANTS(95 )  ! 95! Non-fixing cyanobacteria Temperature correction for respiration rate
!                        KD_CYN_20   =   MODEL_CONSTANTS(96 )  ! 96! Non-fixing cyanobacteria Mortality rate
!                     THETA_KD_CYN   =   MODEL_CONSTANTS(97 )  ! 97! Non-fixing cyanobacteria Temperature correction for Mortality rate
!                      KHS_DIN_CYN   =   MODEL_CONSTANTS(98 )  ! 98! Non-fixing cyanobacteria Half saturation growth for DIN
!                      KHS_DIP_CYN   =   MODEL_CONSTANTS(99 )  ! 99! Non-fixing cyanobacteria Half saturation growth for DIP
!                       KHS_O2_CYN   =   MODEL_CONSTANTS(100)  !100! Non-fixing cyanobacteria Half saturation growth for O2
!                    FRAC_CYN_EXCR   =   MODEL_CONSTANTS(101)  !101! Non-fixing cyanobacteria Fraction of excretion in metabolism rate
!                          I_S_CYN   =   MODEL_CONSTANTS(102)  !102! Non-fixing cyanobacteria Light saturation (langleys)
!               DO_STR_HYPOX_CYN_D   =   MODEL_CONSTANTS(103)  !103! Non-fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!                THETA_HYPOX_CYN_D   =   MODEL_CONSTANTS(104)  !104! Non-fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
!                EXPON_HYPOX_CYN_D   =   MODEL_CONSTANTS(105)  !105! Non-fixing cyanobacteria Exponent constant for Dissolved oxygen stress
!                       CYN_N_TO_C   =   MODEL_CONSTANTS(106)  !106! Non-fixing cyanobacteria Nitrogen to Carbon ratio ,was 0.1
!                       CYN_P_TO_C   =   MODEL_CONSTANTS(107)  !107! Non-fixing cyanobacteria Phosphorus to Carbon ratio
!                      CYN_O2_TO_C   =   MODEL_CONSTANTS(108)  !108! Non-fixing cyanobacteria Oxygen to Carbon ratio for respiration
!                    CYN_C_TO_CHLA   =   MODEL_CONSTANTS(109)  !109! Non-fixing cyanobacteria Carbon to Chlorophyl a ratio
!              KG_FIX_CYN_OPT_TEMP   =   MODEL_CONSTANTS(110)  !110! Fixing cyanobacteria Growth rate
!              FIX_CYN_OPT_TEMP_LR   =   MODEL_CONSTANTS(111)  !111! Fixing Cyanobacteria optimal temperature lower range
!              FIX_CYN_OPT_TEMP_UR   =   MODEL_CONSTANTS(112)  !112! Fixing Cyanobacteria optimal temperature upper range
!               EFF_FIX_CYN_GROWTH   =   MODEL_CONSTANTS(113)  !113! Fixing cyanobacteria Effective growth. (1-EG)*growth - losses for RESP and excretion
!     KAPPA_FIX_CYN_UNDER_OPT_TEMP   =   MODEL_CONSTANTS(114)  !114! Fixing cyanobacteria Temperature correction for growth lower temperature
!      KAPPA_FIX_CYN_OVER_OPT_TEMP   =   MODEL_CONSTANTS(115)  !115! Fixing cyanobacteria Temperature correction for growth upper temperature
!                    KR_FIX_CYN_20   =   MODEL_CONSTANTS(116)  !116! Fixing cyanobacteria RESP rate
!                 THETA_KR_FIX_CYN   =   MODEL_CONSTANTS(117)  !117! Fixing cyanobacteria Temperature correction for RESP rate
!                    KD_FIX_CYN_20   =   MODEL_CONSTANTS(118)  !118! Fixing cyanobacteria Mortality rate of nitrification bacteria
!                 THETA_KD_FIX_CYN   =   MODEL_CONSTANTS(119)  !119! Fixing cyanobacteria Temperature correction for Mortality rate
!                  KHS_DIN_FIX_CYN   =   MODEL_CONSTANTS(120)  !120! Fixing cyanobacteria Half saturation growth for DIN
!                  KHS_DIP_FIX_CYN   =   MODEL_CONSTANTS(121)  !121! Fixing cyanobacteria Half saturation growth for DIP
!                   KHS_O2_FIX_CYN   =   MODEL_CONSTANTS(122)  !122! Fixing cyanobacteria Half saturation growth for O2
!                FRAC_FIX_CYN_EXCR   =   MODEL_CONSTANTS(123)  !123! Fixing cyanobacteria Fraction of excretion in metabolism rate
!                      I_S_FIX_CYN   =   MODEL_CONSTANTS(124)  !124! Fixing cyanobacteria Light saturation (langleys)
!           DO_STR_HYPOX_FIX_CYN_D   =   MODEL_CONSTANTS(125)  !125! Fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!            THETA_HYPOX_FIX_CYN_D   =   MODEL_CONSTANTS(126)  !126! Fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
!            EXPON_HYPOX_FIX_CYN_D   =   MODEL_CONSTANTS(127)  !127! Fixing cyanobacteria Exponent constant for Dissolved oxygen stress
!                   FIX_CYN_N_TO_C   =   MODEL_CONSTANTS(128)  !128! Fixing cyanobacteria Nitrogen to Carbon ratio
!                   FIX_CYN_P_TO_C   =   MODEL_CONSTANTS(129)  !129! Fixing cyanobacteria Phosphorus to Carbon ratio
!                  FIX_CYN_O2_TO_C   =   MODEL_CONSTANTS(130)  !130! Fixing cyanobacteria Oxygen to Carbon ratio for respiration
!                FIX_CYN_C_TO_CHLA   =   MODEL_CONSTANTS(131)  !131! Fixing cyanobacteria Carbon to Chlorophyl a ratio
!                            R_FIX   =   MODEL_CONSTANTS(132)  !132! Fixing cyanobacteria Ratio between non-fixing and fixing fractions growth rate
!                            K_FIX   =   MODEL_CONSTANTS(133)  !133! Fixing cyanobacteria Effectivity parameter of switching to nitrogen fixation
!                  KG_OPA_OPT_TEMP   =   MODEL_CONSTANTS(134)  !134! OtherPhyto Growth rate
!                  OPA_OPT_TEMP_LR   =   MODEL_CONSTANTS(135)  !135! OtherPhyto optimal temperature lower range
!                  OPA_OPT_TEMP_UR   =   MODEL_CONSTANTS(136)  !136! OtherPhyto optimal temperature upper range
!                   EFF_OPA_GROWTH   =   MODEL_CONSTANTS(137)  !137! OtherPhyto Effective growth. (1-EG)*growth - losses for respiration and excretion
!         KAPPA_OPA_UNDER_OPT_TEMP   =   MODEL_CONSTANTS(138)  !138! OtherPhyto Temperature correction for growth lower temperature
!          KAPPA_OPA_OVER_OPT_TEMP   =   MODEL_CONSTANTS(139)  !139! OtherPhyto Temperature correction for growth upper temperature
!                        KR_OPA_20   =   MODEL_CONSTANTS(140)  !140! OtherPhyto Respiration rate
!                     THETA_KR_OPA   =   MODEL_CONSTANTS(141)  !141! OtherPhyto Temperature correction for respiration rate
!                        KD_OPA_20   =   MODEL_CONSTANTS(142)  !142! OtherPhyto Mortality rate
!                     THETA_KD_OPA   =   MODEL_CONSTANTS(143)  !143! OtherPhyto Temperature correction for Mortality rate
!                      KHS_DIN_OPA   =   MODEL_CONSTANTS(144)  !144! OtherPhyto Half saturation growth for DIN
!                      KHS_DIP_OPA   =   MODEL_CONSTANTS(145)  !145! OtherPhyto Half saturation growth for DIP
!                       KHS_O2_OPA   =   MODEL_CONSTANTS(146)  !146! OtherPhyto Half saturation growth for O2
!                    FRAC_OPA_EXCR   =   MODEL_CONSTANTS(147)  !147! OtherPhyto Fraction of excretion in metabolism rate
!                          I_S_OPA   =   MODEL_CONSTANTS(148)  !148! OtherPhyto Light saturation (langleys)
!               DO_STR_HYPOX_OPA_D   =   MODEL_CONSTANTS(149)  !149! OtherPhyto Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!                THETA_HYPOX_OPA_D   =   MODEL_CONSTANTS(150)  !150! OtherPhyto Multiplier of the exponent for Dissolved oxygen stress
!                EXPON_HYPOX_OPA_D   =   MODEL_CONSTANTS(151)  !151! OtherPhyto Exponent constant for Dissolved oxygen stress
!                       OPA_N_TO_C   =   MODEL_CONSTANTS(152)  !152! OtherPhyto Nitrogen to Carbon ratio
!                       OPA_P_TO_C   =   MODEL_CONSTANTS(153)  !153! OtherPhyto Phosphorus to Carbon ratio
!                      OPA_O2_TO_C   =   MODEL_CONSTANTS(154)  !154! OtherPhyto Oxygen to Carbon ratio for respiration
!                    OPA_C_TO_CHLA   =   MODEL_CONSTANTS(155)  !155! OtherPhyto Carbon to Chlorophyl a ratio
!                  KG_ZOO_OPT_TEMP   =   MODEL_CONSTANTS(156)  !156! Zooplankton Growth rate
!                  ZOO_OPT_TEMP_LR   =   MODEL_CONSTANTS(157)  !157! Zooplankton optimal temperature lower range
!                  ZOO_OPT_TEMP_UR   =   MODEL_CONSTANTS(158)  !158! Zooplankton optimal temperature upper range
!                   EFF_ZOO_GROWTH   =   MODEL_CONSTANTS(159)  !159! Zooplankton Effective growth. (1-EG)*growth - losses for respiration and excretion
!         KAPPA_ZOO_UNDER_OPT_TEMP   =   MODEL_CONSTANTS(160)  !160! Zooplankton Temperature correction for growth lower temperature
!          KAPPA_ZOO_OVER_OPT_TEMP   =   MODEL_CONSTANTS(161)  !161! Zooplankton Temperature correction for growth upper temperature
!                     GRAT_ZOO_DIA   =   MODEL_CONSTANTS(162)  !162! Zooplankton Grazing rate (growhth rate multiplier) on diatoms
!                     GRAT_ZOO_CYN   =   MODEL_CONSTANTS(163)  !163! Zooplankton Grazing rate (growhth rate multiplier) on Cyanobacteria
!                     GRAT_ZOO_OPA   =   MODEL_CONSTANTS(164)  !164! Zooplankton Grazing rate (growhth rate multiplier) on fixing Cyanobacteria
!                 GRAT_ZOO_FIX_CYN   =   MODEL_CONSTANTS(165)  !165! Zooplankton Grazing rate (growhth rate multiplier) on OtherPhyto
!            GRAT_ZOO_CHEM_AUT_BAC   =   MODEL_CONSTANTS(166)  !166! Zooplankton Grazing rate (growhth rate multiplier) on NITR_BAC
!             GRAT_ZOO_AER_HET_BAC   =   MODEL_CONSTANTS(167)  !167! Zooplankton Grazing rate (growhth rate multiplier) on AER_HET_BAC
!          GRAT_ZOO_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(168)  !168! Zooplankton Grazing rate (growhth rate multiplier) on DENITR_BAC
!          GRAT_ZOO_DET_PART_ORG_C   =   MODEL_CONSTANTS(169)  !169! Zooplankton Grazing rate (growhth rate multiplier) on part. ORG_C
!                     PREF_ZOO_DIA   =   MODEL_CONSTANTS(170)  !170! Zooplankton Preference for diatoms
!                     PREF_ZOO_CYN   =   MODEL_CONSTANTS(171)  !171! Zooplankton Preference for Cyanobacteria
!                 PREF_ZOO_FIX_CYN   =   MODEL_CONSTANTS(172)  !172! Zooplankton Preference for fixing Cyanobacteria
!                     PREF_ZOO_OPA   =   MODEL_CONSTANTS(173)  !173! Zooplankton Preference for OtherPhyto
!            PREF_ZOO_CHEM_AUT_BAC   =   MODEL_CONSTANTS(174)  !174! Zooplankton Preference for NITR_BAC
!             PREF_ZOO_AER_HET_BAC   =   MODEL_CONSTANTS(175)  !175! Zooplankton Preference for AER_HET_BAC
!          PREF_ZOO_FAC_AN_HET_BAC   =   MODEL_CONSTANTS(176)  !176! Zooplankton Preference for DENITR_BAC
!          PREF_ZOO_DET_PART_ORG_C   =   MODEL_CONSTANTS(177)  !177! Zooplankton Preference for part. ORG_C
!                    KHS_DIA_C_ZOO   =   MODEL_CONSTANTS(178)  !178! Zooplankton Half saturation growth for diatoms
!                    KHS_CYN_C_ZOO   =   MODEL_CONSTANTS(179)  !179! Zooplankton Half saturation growth for Cyanobacteria
!                KHS_FIX_CYN_C_ZOO   =   MODEL_CONSTANTS(180)  !180! Zooplankton Half saturation growth for fixing Cyanobacteria
!                    KHS_OPA_C_ZOO   =   MODEL_CONSTANTS(181)  !181! Zooplankton Half saturation growth for OtherPhyto
!           KHS_CHEM_AUT_BAC_C_ZOO   =   MODEL_CONSTANTS(182)  !182! Zooplankton Half saturation growth for NITR_BAC
!            KHS_AER_HET_BAC_C_ZOO   =   MODEL_CONSTANTS(183)  !183! Zooplankton Half saturation growth for AER_HET_BAC
!         KHS_FAC_AN_HET_BAC_C_ZOO   =   MODEL_CONSTANTS(184)  !184! Zooplankton Half saturation growth for DENITR_BAC
!           KHS_DET_PART_ORG_C_ZOO   =   MODEL_CONSTANTS(185)  !185! Zooplankton Half saturation growth for part. ORG_C
!                     FOOD_MIN_ZOO   =   MODEL_CONSTANTS(186)  !186! Zooplankton Minimum food conc. for feeding
!                           KE_ZOO   =   MODEL_CONSTANTS(187)  !187! not used Zooplankton Excretion rate as growth fraction
!                  FRAC_ZOO_EX_ORG   =   MODEL_CONSTANTS(188)  !188! not used Zooplankton Excretion rate organic fraction
!                        KR_ZOO_20   =   MODEL_CONSTANTS(189)  !189! Zooplankton Respiration rate
!                     THETA_KR_ZOO   =   MODEL_CONSTANTS(190)  !190! Zooplankton Respiration rate Temperature correction
!                        KD_ZOO_20   =   MODEL_CONSTANTS(191)  !191! Zooplankton Mortality rate
!                     THETA_KD_ZOO   =   MODEL_CONSTANTS(192)  !192! Zooplankton Mortality rate Temperature correction
!               DO_STR_HYPOX_ZOO_D   =   MODEL_CONSTANTS(193)  !193! Zooplankton Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!                THETA_HYPOX_ZOO_D   =   MODEL_CONSTANTS(194)  !194! Zooplankton Multiplier of the exponent for Dissolved oxygen stress
!                EXPON_HYPOX_ZOO_D   =   MODEL_CONSTANTS(195)  !195! Zooplankton Exponent constant for Dissolved oxygen stress
!                       ZOO_N_TO_C   =   MODEL_CONSTANTS(196)  !196! Zooplankton Nitrogen to Carbon ratio
!                       ZOO_P_TO_C   =   MODEL_CONSTANTS(197)  !197! Zooplankton Phosphorus to Carbon ratio
!                      ZOO_O2_TO_C   =   MODEL_CONSTANTS(198)  !198! Zooplankton Oxygen to Carbon ratio for respiration
!          KDISS_DET_PART_ORG_C_20   =   MODEL_CONSTANTS(199)  !199! Particulate Detritus Carbon Dissolution rate not dependent on phytoplankton
!       THETA_KDISS_DET_PART_ORG_C   =   MODEL_CONSTANTS(200)  !200! Particulate Detritus Carbon Dissolution rate Temperature correction
!          FAC_PHYT_DET_PART_ORG_C   =   MODEL_CONSTANTS(201)  !201! Particulate Detritus Carbon Phytoplankton linear factor for dissolution rate
!          KDISS_DET_PART_ORG_N_20   =   MODEL_CONSTANTS(202)  !202! Particulate Detritus Nitrogen Dissolution rate not dependent on phytoplankton
!       THETA_KDISS_DET_PART_ORG_N   =   MODEL_CONSTANTS(203)  !203! Particulate Detritus Nitrogen Dissolution rate Temperature correction
!                       KHS_DISS_N   =   MODEL_CONSTANTS(204)  !204! Particulate Detritus Nitrogen dissolution reverse half saturation for DIN
!          FAC_PHYT_DET_PART_ORG_N   =   MODEL_CONSTANTS(205)  !205! Particulate Detritus Nitrogen Phytoplankton linear factor for dissolution rate
!          KDISS_DET_PART_ORG_P_20   =   MODEL_CONSTANTS(206)  !206! Particulate Detritus Phosphorus Dissolution rate not dependent on phytoplankton
!       THETA_KDISS_DET_PART_ORG_P   =   MODEL_CONSTANTS(207)  !207! Particulate Detritus Phosphorus Dissolution rate Temperature correction
!                       KHS_DISS_P   =   MODEL_CONSTANTS(208)  !208! Particulate Detritus Phosphorus  dissolution reverse half saturation for DIP
!          FAC_PHYT_DET_PART_ORG_P   =   MODEL_CONSTANTS(209)  !209! Particulate Detritus Phosphorus  Phytoplankton linear factor for dissolution rate
!                 KDISS_PART_Si_20   =   MODEL_CONSTANTS(210)  !210! Particulate Silica Dissolution rate
!              THETA_KDISS_PART_Si   =   MODEL_CONSTANTS(211)  !211! Particulate Silica Dissolution rate Temperature correction
!                     K_MIN_DOC_20   =   MODEL_CONSTANTS(212)  !212! Dissolved carbon  mineralisation rate
!                  THETA_K_MIN_DOC   =   MODEL_CONSTANTS(213)  !213! Dissolved carbon  mineralisation rate Temperature constant
!                FAC_PHYT_AMIN_DOC   =   MODEL_CONSTANTS(214)  !214! Dissolved carbon  Phytoplankton linear factor for mineralisation rate
!                     K_MIN_DON_20   =   MODEL_CONSTANTS(215)  !215! Dissolved nitrogen  mineralisation rate not dependent on phytoplankton
!                  THETA_K_MIN_DON   =   MODEL_CONSTANTS(216)  !216! Dissolved nitrogen  mineralisation rate Temperature constant
!                       KHS_AMIN_N   =   MODEL_CONSTANTS(217)  !217! Dissolved nitrogen  reverse half saturation for DIN
!                FAC_PHYT_AMIN_DON   =   MODEL_CONSTANTS(218)  !218! Dissolved nitrogen Phytoplankton linear factor for mineralisation rate
!                     K_MIN_DOP_20   =   MODEL_CONSTANTS(219)  !219! Dissolved phosphorus  mineralisation rate not dependent on phytoplankton
!                  THETA_K_MIN_DOP   =   MODEL_CONSTANTS(220)  !220! Dissolved phosphorus  mineralisation rate Temperature constant
!                       KHS_AMIN_P   =   MODEL_CONSTANTS(221)  !221! Dissolved phosphorus reverse half saturation for DIP
!                FAC_PHYT_AMIN_DOP   =   MODEL_CONSTANTS(222)  !222! Dissolved phosphorus Phytoplankton linear factor for mineralisation rate
!                        K_NITR_20   =   MODEL_CONSTANTS(223)  !223! Amonia nitrification rate
!                     KHS_NITR_OXY   =   MODEL_CONSTANTS(224)  !224! Amonia nitrification half saturation for Oxygen
!                   KHS_NITR_NH4_N   =   MODEL_CONSTANTS(225)  !225! Amonia nitrification half saturation for Amonia
!                     THETA_K_NITR   =   MODEL_CONSTANTS(226)  !226! Amonia nitrification rate Temperature constant

     ! all parameters are obtained though not all are used
     call para_get_value('K_A'                             ,                              K_A) !  1 Aeration coefficient (if negative calculates internally)
     call para_get_value('THETA_K_A'                       ,                        THETA_K_A) !  2 Temperature correction factor for aeration
     call para_get_value('KG_CHEM_AUT_BAC_20'              ,               KG_CHEM_AUT_BAC_20) !  3 Chemoautotrophic bacteria Growth rate
     call para_get_value('EFF_CHEM_AUT_BAC_GROWTH'         ,          EFF_CHEM_AUT_BAC_GROWTH) !  4 Chemoautotrophic bacteria growth efficiency
     call para_get_value('THETA_KG_CHEM_AUT_BAC'           ,            THETA_KG_CHEM_AUT_BAC) !  5 Chemoautotrophic bacteria Temperature correction for growth rate
     call para_get_value('KR_CHEM_AUT_BAC_20'              ,               KR_CHEM_AUT_BAC_20) !  6 Chemoautotrophic bacteria Respiration rate
     call para_get_value('THETA_KR_CHEM_AUT_BAC'           ,            THETA_KR_CHEM_AUT_BAC) !  7 Chemoautotrophic bacteria Temperature correction for respiration rate
     call para_get_value('KD_CHEM_AUT_BAC_20'              ,               KD_CHEM_AUT_BAC_20) !  8 Chemoautotrophic bacteria Mortality rate
     call para_get_value('THETA_KD_CHEM_AUT_BAC'           ,            THETA_KD_CHEM_AUT_BAC) !  9 Chemoautotrophic bacteria Temperature correction for Mortality rate
     call para_get_value('KHS_NH4N_CHEM_AUT_BAC'           ,            KHS_NH4N_CHEM_AUT_BAC) ! 10 Chemoautotrophic bacteria Half saturation growth for NH4N
     call para_get_value('KHS_PO4P_CHEM_AUT_BAC'           ,            KHS_PO4P_CHEM_AUT_BAC) ! 11 Chemoautotrophic bacteria Half saturation growth for PO4P
     call para_get_value('KHS_O2_CHEM_AUT_BAC'             ,              KHS_O2_CHEM_AUT_BAC) ! 12 Chemoautotrophic bacteria Half saturation growth for O2
     call para_get_value('DO_STR_HYPOX_CHEM_AUT_BAC_D'     ,      DO_STR_HYPOX_CHEM_AUT_BAC_D) ! 13 Chemoautotrophic bacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_CHEM_AUT_BAC_D'      ,       THETA_HYPOX_CHEM_AUT_BAC_D) ! 14 Chemoautotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_CHEM_AUT_BAC_D'      ,       EXPON_HYPOX_CHEM_AUT_BAC_D) ! 15 Chemoautotrophic bacteria Exponent constant for Dissolved oxygen stress
     call para_get_value('CHEM_AUT_BAC_N_TO_C'             ,              CHEM_AUT_BAC_N_TO_C) ! 16 Chemoautotrophic bacteria Nitrogen to Carbon ratio
     call para_get_value('CHEM_AUT_BAC_P_TO_C'             ,              CHEM_AUT_BAC_P_TO_C) ! 17 Chemoautotrophic bacteria Phosphorus to Carbon ratio
     call para_get_value('CHEM_AUT_BAC_O2_TO_C'            ,             CHEM_AUT_BAC_O2_TO_C) ! 18 Chemoautotrophic bacteria Oxygen to Carbon ratio
     call para_get_value('YIELD_CHEM_AUT_BAC'              ,               YIELD_CHEM_AUT_BAC) ! 19 Chemoautotrophic bacteria Yield of Carbon per unit amonia nitrogen
     call para_get_value('KG_AER_HET_BAC_20'               ,                KG_AER_HET_BAC_20) ! 20 Aerobic heterotrophic bacteria Growth rate
     call para_get_value('EFF_AER_HET_BAC_GROWTH'          ,           EFF_AER_HET_BAC_GROWTH) ! 21 Aerobic heterotrophic bacteria growth efficiency
     call para_get_value('THETA_KG_AER_HET_BAC'            ,             THETA_KG_AER_HET_BAC) ! 22 Aerobic heterotrophic bacteria Temperature correction for growth rate
     call para_get_value('KR_AER_HET_BAC_20'               ,                KR_AER_HET_BAC_20) ! 23 Aerobic heterotrophic bacteria Respiration rate
     call para_get_value('THETA_KR_AER_HET_BAC'            ,             THETA_KR_AER_HET_BAC) ! 24 Aerobic heterotrophic bacteria Temperature correction for respiration rate
     call para_get_value('KD_AER_HET_BAC_20'               ,                KD_AER_HET_BAC_20) ! 25 Aerobic heterotrophic bacteria Mortality rate
     call para_get_value('THETA_KD_AER_HET_BAC'            ,             THETA_KD_AER_HET_BAC) ! 26 Aerobic heterotrophic bacteria Temperature correction for Mortality rate
     call para_get_value('KHS_ORGC_AER_HET_BAC'            ,             KHS_ORGC_AER_HET_BAC) ! 27 Aerobic heterotrophic bacteria Half saturation growth for OC
     call para_get_value('KHS_ORGN_AER_HET_BAC'            ,             KHS_ORGN_AER_HET_BAC) ! 28 Aerobic heterotrophic bacteria Half saturation growth for ON
     call para_get_value('KHS_ORGP_AER_HET_BAC'            ,             KHS_ORGP_AER_HET_BAC) ! 29 Aerobic heterotrophic bacteria Half saturation growth for OP
     call para_get_value('KHS_O2_AER_HET_BAC'              ,               KHS_O2_AER_HET_BAC) ! 30 Aerobic heterotrophic bacteria Half saturation growth for Oxygen
     call para_get_value('KHS_DIN_AER_HET_BAC'             ,              KHS_DIN_AER_HET_BAC) ! 31 Aerobic heterotrophic bacteria Half saturation growth for inorganic nitrogen
     call para_get_value('KHS_DIP_AER_HET_BAC'             ,              KHS_DIP_AER_HET_BAC) ! 32 Aerobic heterotrophic bacteria Half saturation growth for inorganic phosphorus
     call para_get_value('KHS_PHYT_AER_HET_BAC'            ,             KHS_PHYT_AER_HET_BAC) ! 33 Aerobic heterotrophic bacteria Half saturation growth for Phytoplankton C (not used as a resource)
     call para_get_value('YIELD_OC_AER_HET_BAC'            ,             YIELD_OC_AER_HET_BAC) ! 34 Aerobic heterotrophic bacteria Yield of bacteria carbon per unit of organic carbon
     call para_get_value('OX_ORGN_AER_HET_BAC'             ,              OX_ORGN_AER_HET_BAC) ! 35 Aerobic heterotrophic bacteria ON oxidation rate mg N per mg C of bacteria production
     call para_get_value('KHS_MIN_N'                       ,                       KHS_MIN_N ) ! 36 Aerobic heterotrophic bacteria ON mineralisation reverse half saturation for DIN
     call para_get_value('OX_ORGP_AER_HET_BAC'             ,             OX_ORGP_AER_HET_BAC ) ! 37 Aerobic heterotrophic bacteria OP mineralisation rate mg P per mg C of bacteria production
     call para_get_value('KHS_MIN_P'                       ,                       KHS_MIN_P ) ! 38 Aerobic heterotrophic bacteria ON mineralisation reverse half saturation for DIP
     call para_get_value('DO_STR_HYPOX_AER_HET_BAC_D'      ,       DO_STR_HYPOX_AER_HET_BAC_D) ! 39 Aerobic heterotrophic bacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_AER_HET_BAC_D'       ,        THETA_HYPOX_AER_HET_BAC_D) ! 40 Aerobic heterotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_AER_HET_BAC_D'       ,        EXPON_HYPOX_AER_HET_BAC_D) ! 41 Aerobic heterotrophic bacteria Exponent constant for Dissolved oxygen stress
     call para_get_value('AER_HET_BAC_N_TO_C'              ,               AER_HET_BAC_N_TO_C) ! 42 Aerobic heterotrophic bacteria Nitrogen to Carbon ratio
     call para_get_value('AER_HET_BAC_P_TO_C'              ,               AER_HET_BAC_P_TO_C) ! 43 Aerobic heterotrophic bacteria Phosphorus to Carbon ratio
     call para_get_value('AER_HET_BAC_O2_TO_C'             ,              AER_HET_BAC_O2_TO_C) ! 44 Aerobic heterotrophic bacteria Oxygen to Carbon ratio for respiration
     call para_get_value('KG_FAC_AN_HET_BAC_20'            ,             KG_FAC_AN_HET_BAC_20) ! 45 Facultative anaerobic heterotrophic bacteria Growth rate of
     call para_get_value('EFF_FAC_AN_HET_BAC_GROWTH'       ,        EFF_FAC_AN_HET_BAC_GROWTH) ! 46 not used! Facultative anaerobic heterotrophic bacteria growth efficiency
     call para_get_value('THETA_KG_FAC_AN_HET_BAC'         ,          THETA_KG_FAC_AN_HET_BAC) ! 47 not used! Facultative anaerobic heterotrophic bacteria Temperature correction for growth rate
     call para_get_value('KR_FAC_AN_HET_BAC_20'            ,             KR_FAC_AN_HET_BAC_20) ! 48 not used! Facultative anaerobic heterotrophic bacteria Respiration rate
     call para_get_value('THETA_KR_FAC_AN_HET_BAC'         ,          THETA_KR_FAC_AN_HET_BAC) ! 49 not used! Facultative anaerobic heterotrophic bacteria Temperature correction for respiration rate
     call para_get_value('KD_FAC_AN_HET_BAC_20'            ,             KD_FAC_AN_HET_BAC_20) ! 50 not used! Facultative anaerobic heterotrophic bacteria Mortality rate
     call para_get_value('THETA_KD_FAC_AN_HET_BAC'         ,          THETA_KD_FAC_AN_HET_BAC) ! 51 not used! Facultative anaerobic heterotrophic bacteria Temperature correction for Mortality rate
     call para_get_value('KHS_NO3N_FAC_AN_HET_BAC'         ,          KHS_NO3N_FAC_AN_HET_BAC) ! 52 Facultative anaerobic heterotrophic bacteria Half saturation growth for NO3N
     call para_get_value('KHS_ORGC_FAC_AN_HET_BAC'         ,          KHS_ORGC_FAC_AN_HET_BAC) ! 53 not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for OC
     call para_get_value('KHS_ORGN_FAC_AN_HET_BAC'         ,          KHS_ORGN_FAC_AN_HET_BAC) ! 54 not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for ON
     call para_get_value('KHS_ORGP_FAC_AN_HET_BAC'         ,          KHS_ORGP_FAC_AN_HET_BAC) ! 55 not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for OP
     call para_get_value('REV_KHS_O2_FAC_AN_HET_BAC'       ,        REV_KHS_O2_FAC_AN_HET_BAC) ! 56 not used! Facultative anaerobic heterotrophic bacteria Reverse Half saturation growth for O2
     call para_get_value('NO3N_LACK_STR_FAC_AN_HET_BAC_D'  ,   NO3N_LACK_STR_FAC_AN_HET_BAC_D) ! 57 not used! Facultative anaerobic heterotrophic bacteria NO3N stress concentration
     call para_get_value('THETA_NO3_LACK_FAC_AN_HET_BAC_D' ,  THETA_NO3_LACK_FAC_AN_HET_BAC_D) ! 58 not used! Facultative anaerobic heterotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXP_NO3_LACK_FAC_AN_HET_BAC_D'   ,    EXP_NO3_LACK_FAC_AN_HET_BAC_D) ! 59 not used! Facultative anaerobic heterotrophic bacteria Exponent constant for Dissolved oxygen stress
     call para_get_value('FAC_AN_HET_BAC_N_TO_C'           ,            FAC_AN_HET_BAC_N_TO_C) ! 60 not used! Facultative anaerobic heterotrophic bacteria Nitrogen to Carbon ratio
     call para_get_value('FAC_AN_HET_BAC_P_TO_C'           ,            FAC_AN_HET_BAC_P_TO_C) ! 61 not used! Facultative anaerobic heterotrophic bacteria Phosphorus to Carbon ratio
     call para_get_value('FAC_AN_HET_BAC_O2_TO_C'          ,           FAC_AN_HET_BAC_O2_TO_C) ! 62 not used! Facultative anaerobic heterotrophic bacteria Oxygen to Carbon ratio for respiration
     call para_get_value('YIELD_FAC_AN_HET_BAC'            ,             YIELD_FAC_AN_HET_BAC) ! 63 Facultative anaerobic heterotrophic bacteria Yield of carbon per unit nitrates nitrogen
     call para_get_value('KG_DIA_OPT_TEMP'                 ,                  KG_DIA_OPT_TEMP) ! 64 Diatoms Growth rate
     call para_get_value('DIA_OPT_TEMP_LR'                 ,                  DIA_OPT_TEMP_LR) ! 65 Diatoms optimal temperature lower range
     call para_get_value('DIA_OPT_TEMP_UR'                 ,                  DIA_OPT_TEMP_UR) ! 66 Diatoms optimal temperature upper range
     call para_get_value('EFF_DIA_GROWTH'                  ,                   EFF_DIA_GROWTH) ! 67 Diatoms Effective growth. (1-EG)*growth - losses for respiration and excretion
     call para_get_value('KAPPA_DIA_UNDER_OPT_TEMP'        ,         KAPPA_DIA_UNDER_OPT_TEMP) ! 68 Diatoms Temperature correction for growth lower temperature
     call para_get_value('KAPPA_DIA_OVER_OPT_TEMP'         ,          KAPPA_DIA_OVER_OPT_TEMP) ! 69 Diatoms Temperature correction for growth upper temperature
     call para_get_value('KR_DIA_20'                       ,                        KR_DIA_20) ! 70 Diatoms Respiration rate
     call para_get_value('THETA_KR_DIA'                    ,                     THETA_KR_DIA) ! 71 Diatoms Temperature correction for basal respiration rate
     call para_get_value('KD_DIA_20'                       ,                        KD_DIA_20) ! 72 Diatoms Mortality rate
     call para_get_value('THETA_KD_DIA'                    ,                     THETA_KD_DIA) ! 73 Diatoms Temperature correction for Mortality rate
     call para_get_value('KHS_DIN_DIA'                     ,                      KHS_DIN_DIA) ! 74 Diatoms Half saturation growth for DIN
     call para_get_value('KHS_DIP_DIA'                     ,                      KHS_DIP_DIA) ! 75 Diatoms Half saturation growth for DIP
     call para_get_value('KHS_DSi_DIA'                     ,                      KHS_DSi_DIA) ! 76 Diatoms Half saturation growth for DSi
     call para_get_value('KHS_O2_DIA'                      ,                       KHS_O2_DIA) ! 77 Diatoms Half saturation growth for O2
     call para_get_value('FRAC_DIA_EXCR'                   ,                    FRAC_DIA_EXCR) ! 78 Diatoms Fraction of excretion in metabolism rate
     call para_get_value('I_S_DIA'                         ,                          I_S_DIA) ! 79 Diatoms Light saturation (langleys)
     call para_get_value('DO_STR_HYPOX_DIA_D'              ,               DO_STR_HYPOX_DIA_D) ! 80 Diatoms Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_DIA_D'               ,                THETA_HYPOX_DIA_D) ! 81 Diatoms Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_DIA_D'               ,                EXPON_HYPOX_DIA_D) ! 82 Diatoms Exponent constant for Dissolved oxygen stress
     call para_get_value('DIA_N_TO_C'                      ,                       DIA_N_TO_C) ! 83 Diatoms Nitrogen to Carbon ratio
     call para_get_value('DIA_P_TO_C'                      ,                       DIA_P_TO_C) ! 84 Diatoms Phosphorus to Carbon ratio
     call para_get_value('DIA_Si_TO_C'                     ,                      DIA_Si_TO_C) ! 85 Diatoms Silica to Carbon ratio
     call para_get_value('DIA_O2_TO_C'                     ,                      DIA_O2_TO_C) ! 86 Diatoms Oxygen to Carbon ratio for respiration
     call para_get_value('DIA_C_TO_CHLA'                   ,                    DIA_C_TO_CHLA) ! 87 Diatoms Carbon to Chlorophil a ratio
     call para_get_value('KG_CYN_OPT_TEMP'                 ,                  KG_CYN_OPT_TEMP) ! 88 Non-fixing cyanobacteria Growth rate
     call para_get_value('CYN_OPT_TEMP_LR'                 ,                  CYN_OPT_TEMP_LR) ! 89 Non-fixing cyanobacteria optimal temperature lower range
     call para_get_value('CYN_OPT_TEMP_UR'                 ,                  CYN_OPT_TEMP_UR) ! 90 Non-fixing cyanobacteria optimal temperature upper range
     call para_get_value('EFF_CYN_GROWTH'                  ,                   EFF_CYN_GROWTH) ! 91 Non-fixing cyanobacteria Effective growth. (1-EG)*growth - losses for respiration and excretion
     call para_get_value('KAPPA_CYN_UNDER_OPT_TEMP'        ,         KAPPA_CYN_UNDER_OPT_TEMP) ! 92 Non-fixing cyanobacteria Temperature correction for growth lower temperature
     call para_get_value('KAPPA_CYN_OVER_OPT_TEMP'         ,          KAPPA_CYN_OVER_OPT_TEMP) ! 93 Non-fixing cyanobacteria Temperature correction for growth upper temperature
     call para_get_value('KR_CYN_20'                       ,                        KR_CYN_20) ! 94 Non-fixing cyanobacteria Respiration rate
     call para_get_value('THETA_KR_CYN'                    ,                     THETA_KR_CYN) ! 95 Non-fixing cyanobacteria Temperature correction for respiration rate
     call para_get_value('KD_CYN_20'                       ,                        KD_CYN_20) ! 96 Non-fixing cyanobacteria Mortality rate
     call para_get_value('THETA_KD_CYN'                    ,                     THETA_KD_CYN) ! 97 Non-fixing cyanobacteria Temperature correction for Mortality rate
     call para_get_value('KHS_DIN_CYN'                     ,                      KHS_DIN_CYN) ! 98 Non-fixing cyanobacteria Half saturation growth for DIN
     call para_get_value('KHS_DIP_CYN'                     ,                      KHS_DIP_CYN) ! 99 Non-fixing cyanobacteria Half saturation growth for DIP
     call para_get_value('KHS_O2_CYN'                      ,                       KHS_O2_CYN) !100 Non-fixing cyanobacteria Half saturation growth for O2
     call para_get_value('FRAC_CYN_EXCR'                   ,                    FRAC_CYN_EXCR) !101 Non-fixing cyanobacteria Fraction of excretion in metabolism rate
     call para_get_value('I_S_CYN'                         ,                          I_S_CYN) !102 Non-fixing cyanobacteria Light saturation (langleys)
     call para_get_value('DO_STR_HYPOX_CYN_D'              ,               DO_STR_HYPOX_CYN_D) !103 Non-fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_CYN_D'               ,                THETA_HYPOX_CYN_D) !104 Non-fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_CYN_D'               ,                EXPON_HYPOX_CYN_D) !105 Non-fixing cyanobacteria Exponent constant for Dissolved oxygen stress
     call para_get_value('CYN_N_TO_C'                      ,                       CYN_N_TO_C) !106 Non-fixing cyanobacteria Nitrogen to Carbon ratio ,was 0.1
     call para_get_value('CYN_P_TO_C'                      ,                       CYN_P_TO_C) !107 Non-fixing cyanobacteria Phosphorus to Carbon ratio
     call para_get_value('CYN_O2_TO_C'                     ,                      CYN_O2_TO_C) !108 Non-fixing cyanobacteria Oxygen to Carbon ratio for respiration
     call para_get_value('CYN_C_TO_CHLA'                   ,                    CYN_C_TO_CHLA) !109 Non-fixing cyanobacteria Carbon to Chlorophyl a ratio
     call para_get_value('KG_FIX_CYN_OPT_TEMP'             ,              KG_FIX_CYN_OPT_TEMP) !110 Fixing cyanobacteria Growth rate
     call para_get_value('FIX_CYN_OPT_TEMP_LR'             ,              FIX_CYN_OPT_TEMP_LR) !111 Fixing Cyanobacteria optimal temperature lower range
     call para_get_value('FIX_CYN_OPT_TEMP_UR'             ,              FIX_CYN_OPT_TEMP_UR) !112 Fixing Cyanobacteria optimal temperature upper range
     call para_get_value('EFF_FIX_CYN_GROWTH'              ,               EFF_FIX_CYN_GROWTH) !113 Fixing cyanobacteria Effective growth. (1-EG)*growth - losses for RESP and excretion
     call para_get_value('KAPPA_FIX_CYN_UNDER_OPT_TEMP'    ,     KAPPA_FIX_CYN_UNDER_OPT_TEMP) !114 Fixing cyanobacteria Temperature correction for growth lower temperature
     call para_get_value('KAPPA_FIX_CYN_OVER_OPT_TEMP'     ,      KAPPA_FIX_CYN_OVER_OPT_TEMP) !115 Fixing cyanobacteria Temperature correction for growth upper temperature
     call para_get_value('KR_FIX_CYN_20'                   ,                    KR_FIX_CYN_20) !116 Fixing cyanobacteria RESP rate
     call para_get_value('THETA_KR_FIX_CYN'                ,                 THETA_KR_FIX_CYN) !117 Fixing cyanobacteria Temperature correction for RESP rate
     call para_get_value('KD_FIX_CYN_20'                   ,                    KD_FIX_CYN_20) !118 Fixing cyanobacteria Mortality rate of nitrification bacteria
     call para_get_value('THETA_KD_FIX_CYN'                ,                 THETA_KD_FIX_CYN) !119 Fixing cyanobacteria Temperature correction for Mortality rate
     call para_get_value('KHS_DIN_FIX_CYN'                 ,                  KHS_DIN_FIX_CYN) !120 Fixing cyanobacteria Half saturation growth for DIN
     call para_get_value('KHS_DIP_FIX_CYN'                 ,                  KHS_DIP_FIX_CYN) !121 Fixing cyanobacteria Half saturation growth for DIP
     call para_get_value('KHS_O2_FIX_CYN'                  ,                   KHS_O2_FIX_CYN) !122 Fixing cyanobacteria Half saturation growth for O2
     call para_get_value('FRAC_FIX_CYN_EXCR'               ,                FRAC_FIX_CYN_EXCR) !123 Fixing cyanobacteria Fraction of excretion in metabolism rate
     call para_get_value('I_S_FIX_CYN'                     ,                      I_S_FIX_CYN) !124 Fixing cyanobacteria Light saturation (langleys)
     call para_get_value('DO_STR_HYPOX_FIX_CYN_D'          ,           DO_STR_HYPOX_FIX_CYN_D) !125 Fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_FIX_CYN_D'           ,            THETA_HYPOX_FIX_CYN_D) !126 Fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_FIX_CYN_D'           ,            EXPON_HYPOX_FIX_CYN_D) !127 Fixing cyanobacteria Exponent constant for Dissolved oxygen stress
     call para_get_value('FIX_CYN_N_TO_C'                  ,                   FIX_CYN_N_TO_C) !128 Fixing cyanobacteria Nitrogen to Carbon ratio
     call para_get_value('FIX_CYN_P_TO_C'                  ,                   FIX_CYN_P_TO_C) !129 Fixing cyanobacteria Phosphorus to Carbon ratio
     call para_get_value('FIX_CYN_O2_TO_C'                 ,                  FIX_CYN_O2_TO_C) !130 Fixing cyanobacteria Oxygen to Carbon ratio for respiration
     call para_get_value('FIX_CYN_C_TO_CHLA'               ,                FIX_CYN_C_TO_CHLA) !131 Fixing cyanobacteria Carbon to Chlorophyl a ratio
     call para_get_value('R_FIX'                           ,                            R_FIX) !132 Fixing cyanobacteria Ratio between non-fixing and fixing fractions growth rate
     call para_get_value('K_FIX'                           ,                            K_FIX) !133 Fixing cyanobacteria Effectivity parameter of switching to nitrogen fixation
     call para_get_value('KG_OPA_OPT_TEMP'                 ,                  KG_OPA_OPT_TEMP) !134 OtherPhyto Growth rate
     call para_get_value('OPA_OPT_TEMP_LR'                 ,                  OPA_OPT_TEMP_LR) !135 OtherPhyto optimal temperature lower range
     call para_get_value('OPA_OPT_TEMP_UR'                 ,                  OPA_OPT_TEMP_UR) !136 OtherPhyto optimal temperature upper range
     call para_get_value('EFF_OPA_GROWTH'                  ,                   EFF_OPA_GROWTH) !137 OtherPhyto Effective growth. (1-EG)*growth - losses for respiration and excretion
     call para_get_value('KAPPA_OPA_UNDER_OPT_TEMP'        ,         KAPPA_OPA_UNDER_OPT_TEMP) !138 OtherPhyto Temperature correction for growth lower temperature
     call para_get_value('KAPPA_OPA_OVER_OPT_TEMP'         ,          KAPPA_OPA_OVER_OPT_TEMP) !139 OtherPhyto Temperature correction for growth upper temperature
     call para_get_value('KR_OPA_20'                       ,                        KR_OPA_20) !140 OtherPhyto Respiration rate
     call para_get_value('THETA_KR_OPA'                    ,                     THETA_KR_OPA) !141 OtherPhyto Temperature correction for respiration rate
     call para_get_value('KD_OPA_20'                       ,                        KD_OPA_20) !142 OtherPhyto Mortality rate
     call para_get_value('THETA_KD_OPA'                    ,                     THETA_KD_OPA) !143 OtherPhyto Temperature correction for Mortality rate
     call para_get_value('KHS_DIN_OPA'                     ,                      KHS_DIN_OPA) !144 OtherPhyto Half saturation growth for DIN
     call para_get_value('KHS_DIP_OPA'                     ,                      KHS_DIP_OPA) !145 OtherPhyto Half saturation growth for DIP
     call para_get_value('KHS_O2_OPA'                      ,                       KHS_O2_OPA) !146 OtherPhyto Half saturation growth for O2
     call para_get_value('FRAC_OPA_EXCR'                   ,                    FRAC_OPA_EXCR) !147 OtherPhyto Fraction of excretion in metabolism rate
     call para_get_value('I_S_OPA'                         ,                          I_S_OPA) !148 OtherPhyto Light saturation (langleys)
     call para_get_value('DO_STR_HYPOX_OPA_D'              ,               DO_STR_HYPOX_OPA_D) !149 OtherPhyto Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_OPA_D'               ,                THETA_HYPOX_OPA_D) !150 OtherPhyto Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_OPA_D'               ,                EXPON_HYPOX_OPA_D) !151 OtherPhyto Exponent constant for Dissolved oxygen stress
     call para_get_value('OPA_N_TO_C'                      ,                       OPA_N_TO_C) !152 OtherPhyto Nitrogen to Carbon ratio
     call para_get_value('OPA_P_TO_C'                      ,                       OPA_P_TO_C) !153 OtherPhyto Phosphorus to Carbon ratio
     call para_get_value('OPA_O2_TO_C'                     ,                      OPA_O2_TO_C) !154 OtherPhyto Oxygen to Carbon ratio for respiration
     call para_get_value('OPA_C_TO_CHLA'                   ,                    OPA_C_TO_CHLA) !155 OtherPhyto Carbon to Chlorophyl a ratio
     call para_get_value('KG_ZOO_OPT_TEMP'                 ,                  KG_ZOO_OPT_TEMP) !156 Zooplankton Growth rate
     call para_get_value('ZOO_OPT_TEMP_LR'                 ,                  ZOO_OPT_TEMP_LR) !157 Zooplankton optimal temperature lower range
     call para_get_value('ZOO_OPT_TEMP_UR'                 ,                  ZOO_OPT_TEMP_UR) !158 Zooplankton optimal temperature upper range
     call para_get_value('EFF_ZOO_GROWTH'                  ,                   EFF_ZOO_GROWTH) !159 Zooplankton Effective growth. (1-EG)*growth - losses for respiration and excretion
     call para_get_value('KAPPA_ZOO_UNDER_OPT_TEMP'        ,         KAPPA_ZOO_UNDER_OPT_TEMP) !160 Zooplankton Temperature correction for growth lower temperature
     call para_get_value('KAPPA_ZOO_OVER_OPT_TEMP'         ,          KAPPA_ZOO_OVER_OPT_TEMP) !161 Zooplankton Temperature correction for growth upper temperature
     call para_get_value('GRAT_ZOO_DIA'                    ,                     GRAT_ZOO_DIA) !162 Zooplankton Grazing rate (growhth rate multiplier) on diatoms
     call para_get_value('GRAT_ZOO_CYN'                    ,                     GRAT_ZOO_CYN) !163 Zooplankton Grazing rate (growhth rate multiplier) on Cyanobacteria
     call para_get_value('GRAT_ZOO_OPA'                    ,                     GRAT_ZOO_OPA) !164 Zooplankton Grazing rate (growhth rate multiplier) on fixing Cyanobacteria
     call para_get_value('GRAT_ZOO_FIX_CYN'                ,                 GRAT_ZOO_FIX_CYN) !165 Zooplankton Grazing rate (growhth rate multiplier) on OtherPhyto
     call para_get_value('GRAT_ZOO_CHEM_AUT_BAC'           ,            GRAT_ZOO_CHEM_AUT_BAC) !166 Zooplankton Grazing rate (growhth rate multiplier) on NITR_BAC
     call para_get_value('GRAT_ZOO_AER_HET_BAC'            ,             GRAT_ZOO_AER_HET_BAC) !167 Zooplankton Grazing rate (growhth rate multiplier) on AER_HET_BAC
     call para_get_value('GRAT_ZOO_FAC_AN_HET_BAC'         ,          GRAT_ZOO_FAC_AN_HET_BAC) !168 Zooplankton Grazing rate (growhth rate multiplier) on DENITR_BAC
     call para_get_value('GRAT_ZOO_DET_PART_ORG_C'         ,          GRAT_ZOO_DET_PART_ORG_C) !169 Zooplankton Grazing rate (growhth rate multiplier) on part. ORG_C
     call para_get_value('PREF_ZOO_DIA'                    ,                     PREF_ZOO_DIA) !170 Zooplankton Preference for diatoms
     call para_get_value('PREF_ZOO_CYN'                    ,                     PREF_ZOO_CYN) !171 Zooplankton Preference for Cyanobacteria
     call para_get_value('PREF_ZOO_FIX_CYN'                ,                 PREF_ZOO_FIX_CYN) !172 Zooplankton Preference for fixing Cyanobacteria
     call para_get_value('PREF_ZOO_OPA'                    ,                     PREF_ZOO_OPA) !173 Zooplankton Preference for OtherPhyto
     call para_get_value('PREF_ZOO_CHEM_AUT_BAC'           ,            PREF_ZOO_CHEM_AUT_BAC) !174 Zooplankton Preference for NITR_BAC
     call para_get_value('PREF_ZOO_AER_HET_BAC'            ,             PREF_ZOO_AER_HET_BAC) !175 Zooplankton Preference for AER_HET_BAC
     call para_get_value('PREF_ZOO_FAC_AN_HET_BAC'         ,          PREF_ZOO_FAC_AN_HET_BAC) !176 Zooplankton Preference for DENITR_BAC
     call para_get_value('PREF_ZOO_DET_PART_ORG_C'         ,          PREF_ZOO_DET_PART_ORG_C) !177 Zooplankton Preference for part. ORG_C
     call para_get_value('KHS_DIA_C_ZOO'                   ,                    KHS_DIA_C_ZOO) !178 Zooplankton Half saturation growth for diatoms
     call para_get_value('KHS_CYN_C_ZOO'                   ,                    KHS_CYN_C_ZOO) !179 Zooplankton Half saturation growth for Cyanobacteria
     call para_get_value('KHS_FIX_CYN_C_ZOO'               ,                KHS_FIX_CYN_C_ZOO) !180 Zooplankton Half saturation growth for fixing Cyanobacteria
     call para_get_value('KHS_OPA_C_ZOO'                   ,                    KHS_OPA_C_ZOO) !181 Zooplankton Half saturation growth for OtherPhyto
     call para_get_value('KHS_CHEM_AUT_BAC_C_ZOO'          ,           KHS_CHEM_AUT_BAC_C_ZOO) !182 Zooplankton Half saturation growth for NITR_BAC
     call para_get_value('KHS_AER_HET_BAC_C_ZOO'           ,            KHS_AER_HET_BAC_C_ZOO) !183 Zooplankton Half saturation growth for AER_HET_BAC
     call para_get_value('KHS_FAC_AN_HET_BAC_C_ZOO'        ,         KHS_FAC_AN_HET_BAC_C_ZOO) !184 Zooplankton Half saturation growth for DENITR_BAC
     call para_get_value('KHS_DET_PART_ORG_C_ZOO'          ,           KHS_DET_PART_ORG_C_ZOO) !185 Zooplankton Half saturation growth for part. ORG_C
     call para_get_value('FOOD_MIN_ZOO'                    ,                     FOOD_MIN_ZOO) !186 Zooplankton Minimum food conc. for feeding
     call para_get_value('KE_ZOO'                          ,                           KE_ZOO) !187 not used Zooplankton Excretion rate as growth fraction
     call para_get_value('FRAC_ZOO_EX_ORG'                 ,                  FRAC_ZOO_EX_ORG) !188 not used Zooplankton Excretion rate organic fraction
     call para_get_value('KR_ZOO_20'                       ,                        KR_ZOO_20) !189 Zooplankton Respiration rate
     call para_get_value('THETA_KR_ZOO'                    ,                     THETA_KR_ZOO) !190 Zooplankton Respiration rate Temperature correction
     call para_get_value('KD_ZOO_20'                       ,                        KD_ZOO_20) !191 Zooplankton Mortality rate
     call para_get_value('THETA_KD_ZOO'                    ,                     THETA_KD_ZOO) !192 Zooplankton Mortality rate Temperature correction
     call para_get_value('DO_STR_HYPOX_ZOO_D'              ,               DO_STR_HYPOX_ZOO_D) !193 Zooplankton Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
     call para_get_value('THETA_HYPOX_ZOO_D'               ,                THETA_HYPOX_ZOO_D) !194 Zooplankton Multiplier of the exponent for Dissolved oxygen stress
     call para_get_value('EXPON_HYPOX_ZOO_D'               ,                EXPON_HYPOX_ZOO_D) !195 Zooplankton Exponent constant for Dissolved oxygen stress
     call para_get_value('ZOO_N_TO_C'                      ,                       ZOO_N_TO_C) !196 Zooplankton Nitrogen to Carbon ratio
     call para_get_value('ZOO_P_TO_C'                      ,                       ZOO_P_TO_C) !197 Zooplankton Phosphorus to Carbon ratio
     call para_get_value('ZOO_O2_TO_C'                     ,                      ZOO_O2_TO_C) !198 Zooplankton Oxygen to Carbon ratio for respiration
     call para_get_value('KDISS_DET_PART_ORG_C_20'         ,          KDISS_DET_PART_ORG_C_20) !199 Particulate Detritus Carbon Dissolution rate not dependent on phytoplankton
     call para_get_value('THETA_KDISS_DET_PART_ORG_C'      ,       THETA_KDISS_DET_PART_ORG_C) !200 Particulate Detritus Carbon Dissolution rate Temperature correction
     call para_get_value('FAC_PHYT_DET_PART_ORG_C'         ,          FAC_PHYT_DET_PART_ORG_C) !201 Particulate Detritus Carbon Phytoplankton linear factor for dissolution rate
     call para_get_value('KDISS_DET_PART_ORG_N_20'         ,          KDISS_DET_PART_ORG_N_20) !202 Particulate Detritus Nitrogen Dissolution rate not dependent on phytoplankton
     call para_get_value('THETA_KDISS_DET_PART_ORG_N'      ,       THETA_KDISS_DET_PART_ORG_N) !203 Particulate Detritus Nitrogen Dissolution rate Temperature correction
     call para_get_value('KHS_DISS_N'                      ,                       KHS_DISS_N) !204 Particulate Detritus Nitrogen dissolution reverse half saturation for DIN
     call para_get_value('FAC_PHYT_DET_PART_ORG_N'         ,          FAC_PHYT_DET_PART_ORG_N) !205 Particulate Detritus Nitrogen Phytoplankton linear factor for dissolution rate
     call para_get_value('KDISS_DET_PART_ORG_P_20'         ,          KDISS_DET_PART_ORG_P_20) !206 Particulate Detritus Phosphorus Dissolution rate not dependent on phytoplankton
     call para_get_value('THETA_KDISS_DET_PART_ORG_P'      ,       THETA_KDISS_DET_PART_ORG_P) !207 Particulate Detritus Phosphorus Dissolution rate Temperature correction
     call para_get_value('KHS_DISS_P'                      ,                       KHS_DISS_P) !208 Particulate Detritus Phosphorus  dissolution reverse half saturation for DIP
     call para_get_value('FAC_PHYT_DET_PART_ORG_P'         ,          FAC_PHYT_DET_PART_ORG_P) !209 Particulate Detritus Phosphorus  Phytoplankton linear factor for dissolution rate
     call para_get_value('KDISS_PART_Si_20'                ,                 KDISS_PART_Si_20) !210 Particulate Silica Dissolution rate
     call para_get_value('THETA_KDISS_PART_Si'             ,              THETA_KDISS_PART_Si) !211 Particulate Silica Dissolution rate Temperature correction
     call para_get_value('K_MIN_DOC_20'                    ,                     K_MIN_DOC_20) !212 Dissolved carbon  mineralisation rate
     call para_get_value('THETA_K_MIN_DOC'                 ,                  THETA_K_MIN_DOC) !213 Dissolved carbon  mineralisation rate Temperature constant
     call para_get_value('FAC_PHYT_AMIN_DOC'               ,                FAC_PHYT_AMIN_DOC) !214 Dissolved carbon  Phytoplankton linear factor for mineralisation rate
     call para_get_value('K_MIN_DON_20'                    ,                     K_MIN_DON_20) !215 Dissolved nitrogen  mineralisation rate not dependent on phytoplankton
     call para_get_value('THETA_K_MIN_DON'                 ,                  THETA_K_MIN_DON) !216 Dissolved nitrogen  mineralisation rate Temperature constant
     call para_get_value('KHS_AMIN_N'                      ,                       KHS_AMIN_N) !217 Dissolved nitrogen  reverse half saturation for DIN
     call para_get_value('FAC_PHYT_AMIN_DON'               ,                FAC_PHYT_AMIN_DON) !218 Dissolved nitrogen Phytoplankton linear factor for mineralisation rate
     call para_get_value('K_MIN_DOP_20'                    ,                     K_MIN_DOP_20) !219 Dissolved phosphorus  mineralisation rate not dependent on phytoplankton
     call para_get_value('THETA_K_MIN_DOP'                 ,                  THETA_K_MIN_DOP) !220 Dissolved phosphorus  mineralisation rate Temperature constant
     call para_get_value('KHS_AMIN_P'                      ,                       KHS_AMIN_P) !221 Dissolved phosphorus reverse half saturation for DIP
     call para_get_value('FAC_PHYT_AMIN_DOP'               ,                FAC_PHYT_AMIN_DOP) !222 Dissolved phosphorus Phytoplankton linear factor for mineralisation rate
     call para_get_value('K_NITR_20'                       ,                        K_NITR_20) !223 Amonia nitrification rate
     call para_get_value('KHS_NITR_OXY'                    ,                     KHS_NITR_OXY) !224 Amonia nitrification half saturation for Oxygen
     call para_get_value('KHS_NITR_NH4_N'                  ,                   KHS_NITR_NH4_N) !225 Amonia nitrification half saturation for Amonia
     call para_get_value('THETA_K_NITR'                    ,                     THETA_K_NITR) !226 Amonia nitrification rate Temperature constant                    
                     

 ! Calculating and sending Total N, Total P and Total IN to the output

 !       Total N

          if (nstate+4 .gt. noutput) then
          print *, 'DERIVED_VARS: noutput is to small'
          print *, nstate+3, noutput
          stop
         end if

         WC_OUTPUTS(1:nkn,nstate+1) =    NH4_N    +   NO3_N + &            !Dissolved inorganic
                                 CHEM_AUT_BAC_C     * CHEM_AUT_BAC_N_TO_C   + &   !Bacteria
                                 AER_HET_BAC_C      * AER_HET_BAC_N_TO_C    + &
                                 FAC_AN_HET_BAC_C   * FAC_AN_HET_BAC_N_TO_C + &
                                 DIA_C * DIA_N_TO_C + &                            !Phyto
                                 CYN_C * CYN_N_TO_C + &
                                 OPA_C * OPA_N_TO_C + &
                                 FIX_CYN_C      * FIX_CYN_N_TO_C + &
                                 DET_PART_ORG_N + DISS_ORG_N +   &                !Part + Diss
                                 ZOO_C          * ZOO_N_TO_C                             !Zoo


!        Total P


         WC_OUTPUTS(1:nkn,nstate+2) = PO4_P + &                                 !Phosphates
                                CHEM_AUT_BAC_C   * CHEM_AUT_BAC_P_TO_C  + &     !Bacteria
                                AER_HET_BAC_C    * AER_HET_BAC_P_TO_C   + &
                                FAC_AN_HET_BAC_C * FAC_AN_HET_BAC_P_TO_C+ &
                                DIA_C * DIA_P_TO_C  + &                        !Phyto
                                CYN_C * CYN_P_TO_C + &
                                OPA_C * OPA_P_TO_C + &
                                FIX_CYN_C      * CYN_P_TO_C + &
                                DET_PART_ORG_P + DISS_ORG_P + &         !Part + Diss
                                ZOO_C          * ZOO_P_TO_C                    !Zoo



         WC_OUTPUTS(1:nkn,nstate+3) = NH4_N + NO3_N                                  ! Total IN
         WC_OUTPUTS(1:nkn,nstate+4) = pH
   end
!************************************************************************
!************************************************************************

!*********************************************'
!*                                           *'
!*            AUXILLARY FUNCTIONS            *'
!*                                           *'
!*********************************************'

!Subroutine to calculate actual growth
subroutine GROWTH_AT_TEMP(TEMP,K_GROWTH, Lower_TEMP,Upper_TEMP, K_AT_OPT_TEMP,  &
                          KAPPA_UNDER_OPT_TEMP, KAPPA_OVER_OPT_TEMP,nkn) !result(K_GROWTH)
! Output:
! K_GROWTH,
!
! Inputs:
! TEMP,
! Lower_TEMP,              Constant
! Upper_TEMP,              Constant
! K_AT_OPT_TEMP,           Constant
! KAPPA_UNDER_OPT_TEMP,    Constant
! KAPPA_OVER_OPT_TEMP,     Constant
! nkn                      Constant

    implicit none

    integer nkn


    double precision, intent(in) :: TEMP          (nkn)
    double precision             :: K_GROWTH      (nkn)

    double precision, intent(in) :: Lower_TEMP,Upper_TEMP
    double precision, intent(in) :: K_AT_OPT_TEMP
    double precision, intent(in) :: KAPPA_UNDER_OPT_TEMP
    double precision, intent(in) :: KAPPA_OVER_OPT_TEMP

    where (TEMP <= Lower_TEMP)
        K_GROWTH = K_AT_OPT_TEMP * &
            exp((-1.0D0) * KAPPA_UNDER_OPT_TEMP * dabs(Lower_TEMP - TEMP))
    end where

    where ((TEMP>Lower_TEMP).and.(TEMP<Upper_TEMP))
        K_GROWTH = K_AT_OPT_TEMP
    end where

    where (TEMP >= Upper_TEMP)
        K_GROWTH = K_AT_OPT_TEMP * &
            exp((-1.0D0) * KAPPA_OVER_OPT_TEMP * dabs(Upper_TEMP - TEMP))
    end where

    !GROWTH_AT_TEMP = K_GROWTH

end subroutine GROWTH_AT_TEMP


!Function to calculate ammonia preference (WASP)
subroutine AMMONIA_PREFS(AMMONIA_PREF,NH3, NOx, kn,nkn)

    integer nkn

    double precision :: AMMONIA_PREF(nkn)
    double precision, intent(in) :: NH3(nkn)
    double precision, intent(in) :: NOx(nkn)
    double precision, intent(in) :: kn
    double precision             :: PN (nkn)

!     IF (NH3 .lt. 1.0D-6) THEN
!      PN = 0.0D0
!     ELSE
!     PN = (NH3 * NOx) / ((kn + NH3) * (kn + NOx)) + &
!         (kn * NH3) / ((NH3 + NOx) * (kn + NOx))
!     END IF

    where (NH3 .lt. 1.0D-6)
     PN = 0.0D0
    elsewhere
     PN = (NH3 * NOx) / ((kn + NH3) * (kn + NOx)) + &
        (kn * NH3) / ((NH3 + NOx) * (kn + NOx))
    end where

    AMMONIA_PREF = PN


end subroutine AMMONIA_PREFS

!Function, which returns saturation concentration of dissolved oxygen
double precision function DO_SATURATION(T, S, H) !result(CS)

    !Water temperature (in Celcius)
    double precision, intent(in) :: T

    !Salinity (in ppt)
    double precision, intent(in) :: S

    !Elevation (in m)
    double precision, intent(in) :: H


    !Water temperature (in Kelvin)
    double precision :: T_KELVIN

    !Altitude (in feet)
    double precision :: H_FEET

    double precision :: LN_CSF
    double precision :: LN_CSS
    double precision :: CSS
    double precision :: CSP

    !double precision :: CS
    double precision       :: CS
    !Pressure at altitude H (in atm)
    double precision :: P

    !Standart pressure (in mmHg)
    double precision :: P0

    double precision :: LN_PWV

    !Partial pressure of water vapor (in atm)
    double precision :: PWV

    !A constant
    double precision :: THETA

    !print *, 'TEMP',T,'SALT',S,'ELEVATION',H
    T_KELVIN = T + 273.15
    H_FEET = H / 0.3048

    !Calculate the effect of temperature on dissolved oxygen saturation
    LN_CSF = -139.34411 + (157570.1 / T_KELVIN) - &
    &        (66423080.0 / (T_KELVIN ** 2.0D0)) + &
    &        (12438000000.0 / (T_KELVIN ** 3.0D0)) - &
    &        (862194900000.0 / (T_KELVIN ** 4.0D0))

    !Calculate the effect of salinity on dissolved oxygen saturation
    LN_CSS = LN_CSF - S * &
    &        (0.017674 - (10.754 / T_KELVIN) + &
    &         (2140.7 / (T_KELVIN ** 2.0D0)))

    CSS = exp(LN_CSS)

    !Calculate the effect of altitude on dissolved oxygen saturation

    !Calculate THETA
    THETA = 0.000975 - (0.00001426 * T) + (0.00000006436 * (T ** 2.0D0))

    !Set standard pressure to mean sea level
    P0 = 760.0

    !Calculate atmospheric pressure at altitude H
    P = (P0 - (0.02667 * H_FEET)) / 760.0

    !Calculate vapour pressure of water(DIKKAT)
    LN_PWV = 11.8571 - (3840.7 / T_KELVIN) - &
    &        (216961.0 / (T_KELVIN ** 2.0D0))

    PWV = exp(LN_PWV)

    !Final calculation including altitude effect
    CSP = CSS  * P * (((1.0D0 - (PWV / P)) * (1.0D0 - (THETA * P))) &
    &   / ((1 - PWV) * (1.0D0 - THETA)))

    CS = CSP
    DO_SATURATION = CS
end function DO_SATURATION


!Function to calculate kind based reareation constant
!Borrowed from EUTRO5, Ambrose et al., 1993
double precision function KAWIND(WINDS, TW, TA, DEPTH, WTYPE) !result(KA_W)

    !WS         wind speed, m/s
    !TW         water temperature C
    !TA         air temperature C
    !DEPTH      defined depth(segmax) in geometrical
    !WTYPE      type of water body

    double precision, intent(in) :: WINDS
    double precision, intent(in) :: TW
    double precision, intent(in) :: TA
    double precision, intent(in) :: DEPTH
    double precision, intent(in) :: WTYPE

    double precision :: KA_W

    double precision :: WS

    !RK         reareation term calculated in kawind
    double precision :: RK
    double precision :: R0MIN

    integer :: IWTYPE
    integer :: N

    !UT     : SHEAR VELOCITY (CM/SEC)
    !UC     : CRITICAL SHEAR VELOCITY (CM/SEC)
    !KARMAN : VONKARMAN'S CONSTANT
    !ZE     : EQUILIBRIUM ROUGHNESS (CM)
    !1/LAM  : A REYNOLD'S NUMBER
    !GAM    : NONDIMENSIONAL COEFFICIENT DEPENDENT ON WATER BODY SIZE.
    double precision :: UT
    double precision :: UC
    double precision :: KARMAN
    double precision :: ZE
    double precision :: LAM
    double precision :: GAM

    !DifF : DifFUSIVITY OF OXYGEN IN WATER (CM**2/SEC)
    !VW   : VISCOSITY OF WATER             (CM**2/SEC)
    !VA   : VISCOSITY OF AIR               (CM**2/SEC)
    !PW   : DENSITY OF WATER               (G/CM**3)
    !PA   : DENSITY OF AIR                 (G/CM**3)
    double precision :: DifF
    double precision :: VW
    double precision :: VA
    double precision :: PA
    double precision :: PW

    double precision :: KA3
    double precision :: WH
    double precision :: SRCD
    double precision :: ERR
    double precision :: EF
    double precision :: F1
    double precision :: F2
    double precision :: FP1
    double precision :: FP2
    double precision :: FP3
    double precision :: FP4
    double precision :: SRCD2
    double precision :: CDDRAG
    double precision :: US
    double precision :: Z0
    double precision :: RK1
    double precision :: RK2
    double precision :: RK3
    double precision :: GAMU

    R0MIN = 1.e-15
    IWTYPE = int(WTYPE)
    WS = WINDS

    if (IWTYPE.EQ.3) then
        UT  = 10.0
        UC  = 11.0
        ZE  = 0.35
        LAM = 3.0
        GAM = 5.0
    else
        if (IWTYPE .EQ.1) then
            UT  = 9.0
            UC  = 22.0
            ZE  = 0.25
            LAM = 10.0
            GAM = 10.0
        else
            if (IWTYPE.EQ.2) then
                UT  = 10.0
                UC  = 11.0
                ZE  = 0.25
                LAM = 3.0
                GAM = 6.5
            else
                WRITE(*,*) 'KAWIND : WRONG VALUE FOR WATERBODY TYPE'
                STOP
            end if
        end if

    end if


    DifF = 4.58E-07 * TW + 1.2E-05
    VW   = 0.0164 - 0.00024514*TW
    VA   = 0.133 + 0.0009*TA
    PA   = 0.00129 - 0.0000040*TA
    PW   = 1.00
    WS   = WS * 100.0
    RK   = 1.0

    !NEWTON RAPHSON METHOD TO CALCULATE THE SQUARE ROOT OF THE DRAG
    !COEFFICIENT
    N = 0

    KARMAN = 0.4
    KA3    = KARMAN**0.3333
    WH     = 1000.0

    !INITIAL GUESS FOR SQUARE ROOT OF THE DRAG COEFFICIENT
    SRCD = 0.04
    ERR  = 1.0

    do N = 1, 9
        !CALCULATE VALUE OF FUNCTION(F2) AND
        !DERIVATIVE OF FUNCTION(FP)
        EF  = exp( - SRCD*WS/UT)
        F1  = LOG((WH/ZE) + (WH*LAM/VA)*SRCD*WS*EF)
        F2  = F1 - KARMAN / SRCD
        FP1 = 1.0/((WH/ZE) + (LAM*WH/VA)*SRCD*WS*EF)
        FP2 = ((WH*LAM)/(VA*UT))*SRCD*(WS**2.0)*EF
        FP3 = (WH*(LAM / VA)) * WS * EF
        FP4 = FP1*(FP2 + FP3) + (KARMAN/(SRCD**2.0))

       !A NEW GUESS FOR SQUARE ROOT OF DRAG AND COMPARE TO
       !PREVIOUS GUESS AND LOOP BACK THROUGH N-R WITH NEW GUESS
       !if APPROPRIATE
       SRCD2 = SRCD - F2/FP4
       ERR = ABS(SRCD - SRCD2)
       SRCD = SRCD2

       if (ERR.LE.0.0005) then
           EXIT
       end if

    end do

    if ((ERR.GT.0.005).AND.(N.EQ.9)) then
        WRITE(*,*) 'KAWIND : SOLUTION DID NOT CONVERGE'
        STOP
    end if

    CDDRAG = SRCD**2.0
    US     = SRCD * WS
    Z0     = 1.0 / ((1.0 / ZE) + LAM * US * exp(-US / UT) / VA)
    WS     = WS / 100.0

    if (WS.LT.6.0) then
        RK1 = ((DifF / VW)**0.666667) * SRCD * ((PA / PW)**0.5)
        RK  = RK1 * KA3 * (WS/GAM)
        RK  = RK * 3600.0 *24.0
        RK  = RK / DEPTH
    end if

    if ((WS.GE.6.0).AND.(WS.LE.20.0)) then
        GAMU = GAM * US * (exp(-(US / UC) + 1.0) / UC)
        RK1  = ((DifF/VW)**0.6667) * KA3 * ((PA/PW)**0.5) * (US / GAMU)
        RK2  = ((DifF * US * PA * VA) / (KARMAN * Z0 * PW * VW))**0.5
        RK3  = (1.0 / RK1) + (1.0 / RK2)
        RK   = 1.0 / RK3
        RK   = RK * 3600.0 * (24.0 / 100.0)
        RK   = RK / DEPTH
    end if

    if (WS.GT.20.0) then
        RK = ((DifF * PA * VA * US) / (KARMAN * ZE * PW * VW))**0.5
        RK = RK * 3600.0 * (24.0 / 100.0)
        RK = RK / DEPTH
    end if

    KA_W = RK
    KAWIND = KA_W

end function KAWIND

!************************************************************************

!********************************************************************
!********************************************************************
subroutine LIM_LIGHT(Ia, TCHLA, GITMAX, H, ke, LLIGHT, CCHL_RATIO,LIGHT_SAT,nkn)

!   Dick-Smith light limitation
!
!   Version to use in ALUKAS with instanteneous light intensity and constant C to Chla ratio
!   Vectorized version
!
!   Double precision variables
!
!   Parameters:
! Inputs:
!      Ia         -   Instanteneous light intensity (langleys), renamed to ITOT in the program
!      TCHLA      -   total chlorophyl for all phytoplankton groups
!      CCHL_RATIO -  carbon to chlorophyl a ratio
!      GITMAX     - temperature corrected maximum relative growth for phytoplankton group
!      H          - depth
!      ke         - light extinction coefficient for the water free of phytoplankton
! Outputs:
!      LLIGHT     - light limitation factor
!      LIGHT_SAT  - saturation light intensity, returned for control

! Constants:
!      KC    - HARDCODED. Chlorophyll light extinction coefficient (1/m)
!      PHIMX - HARDCODED. Max. Quantum Yield


      use para_aqua
      
      implicit none
      integer nkn  ! number of nodes in grid

      double precision Ia    (nkn)
      double precision TCHLA (nkn)
      double precision GITMAX(nkn)
      double precision H     (nkn)
      double precision ke    (nkn)

      double precision CCHL_RATIO

      double precision KC
      double precision PHIMX

      double precision LLIGHT (nkn)

      double precision LIGHT_SAT (nkn)
      double precision KESHD     (nkn)
      double precision SKE       (nkn)
      double precision TEMP1     (nkn)
      double precision TEMP2     (nkn)
      double precision TEMP3     (nkn)
      double precision PI


      PI = 3.14159D0
      call para_get_value('XKC', KC)
      !KC     =  0.010D0 !0.016D0    ! Chloroph. extinction, ( mcg Chla/l/m)
      call para_get_value('PHIMX', PHIMX)   
      !=  720.0D0    !PHY   Quantum yield const. mg C/mole photon



!      KESHD = KC * 1000.0 * TCHLA
      KESHD = KC * TCHLA
      SKE = ke
      SKE = SKE + KESHD
      TEMP1 = SKE * H

      if (any(GITMAX .lt. 1D-20) .or. CCHL_RATIO .lt. 1D-20) then
          write(6,*) 'LIM_LIGHT: TEMP2 is NaN ', 'GITMAX=', GITMAX, 'CCHL_RATIO=', CCHL_RATIO
          stop
      end if


      TEMP2 = (0.083D0 * PHIMX * KC) / (GITMAX * CCHL_RATIO * 2.718D0) !1/TEMP2 - the saturating light intensity

      LIGHT_SAT = 1.0D0/TEMP2

      TEMP3 = EXP( - TEMP1)                                        ! fraction of the light at bottom

      LLIGHT   = (2.7183D0 / TEMP1) * &
              (EXP( -TEMP2 * Ia * TEMP3) - EXP( -TEMP2 * Ia))

end subroutine LIM_LIGHT

!********************************************************************
!********************************************************************

!********************************************************************
     SUBROUTINE cur_smith(Ia,TCHLA,CCHLXI,GITMAX,H,ke, LLIGHT,CCHLX)

!    Can not be used for instanteneous light
!    while total dayllight is unknown to adjust C to Chla ratio!

!   Dick-Smith light limitation formulation(adapted from EUTRO).
!
!   Version to use in ALUKAS with instanteneous light intensity and variable(calculated inside) C to Chla ratio
!   Pytoplankton growth is modelled full day
!   Double precision variables
!
!   Parameters:
!      PHOTO - fraction of the day with light, renamed by FDAY later in the program.
!      Ia    - Instanteneous light intensity (langleys), renamed to ITOT in the program
!      TCHLA - total chlorophyl for all phytoplankton groups
!      CCHLXI- carbon to chlorophyl ratio for the current step
!      GITMAX - temperature corrected maximum relative growth for phytoplankton group
!      H      - depth
!      ke     - light extinction coefficient for the water free of phytoplankton
!      CCHLX  - carbon to chlorophyl ratio calculated by subroutine for the next time step
!      XKC    - HARDCODED. Chlorophyll light extinction coefficient (1/m)
!      PHIMAX - HARDCODED. Max. Quantum Yield


          double precision PHOTO
          double precision Ia
          double precision TCHLA
          double precision CCHLXI
          double precision GITMAX
          double precision H
          double precision ke
          double precision XKC
          double precision PHIMX

          double precision LLIGHT
          double precision CCHLX


          double precision FDAY
          double precision ITOT
          double precision CCHL1
          double precision KESHD
          double precision SKE
          double precision TEMP1
          double precision TEMP2
          double precision TEMP3
          double precision IMAX
          double precision SUM
          double precision DTDAY
          double precision I0
          double precision RLIGHT
          double precision IAV
          double precision IAVSG


          double precision PI
          INTEGER I


          PI = 3.14159D0

          XKC     =  0.016   ! Chloroph. extinction, ( mcg Chla/l/m)
                             ! 0.04 is approximatelly the value that corresponds to the angle of curve given in Chapra
          PHIMX   =  720.0   !PHY   Quantum yield const. mg C/mole photon

          FDAY=1.0D0
          ITOT = Ia
!         CCHL1 = CCHLX(ISEG)
          CCHL1 = CCHLXI

!         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!         %?IAVBOT=IAVBOTX(ISEG)
!         %PHYT=TPHY
!         %GITMAX=k1c*rtmult(TEMP, t_1, t_2, t_3, t_4, k_1, 0.98, 0.98, k_4)
!         %TCHLA = PHYT/CCHL1
!         %RLIGHT = RLGHTS (ISEG, 1)
!         %
!         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!
          KESHD = XKC * TCHLA   !Chla commes in micrograms
          SKE = ke
          SKE = SKE + KESHD
          TEMP1 = SKE * H

          if (GITMAX .lt. 1D-20 .or. CCHL1 .lt. 1D-20) then
           write(6,*) 'SMITH: TEMP2 is NaN ', 'GITMAX=', GITMAX,&
                    'CCHL1=', CCHL1
           stop
          end if

          TEMP2 = 0.083D0 * PHIMX * XKC / (GITMAX * CCHL1 * 2.718D0) !1/TEMP2 - the saturating light intensity
          TEMP3 = EXP( - TEMP1) !fraction of the light at bottom

! Light limitation varies during the day
            RLIGHT   = 2.7183D0 / TEMP1 * &
              (EXP( -TEMP2 * ITOT * TEMP3) - EXP( -TEMP2 * ITOT))
            LLIGHT=RLIGHT

!         Adapt carbon to chlorophyll ratio:
          IAV=0.9D0 * ITOT/FDAY              ! It can not be used for instantenous light because total light is unknown!
          IAVSG=IAV*(1.0D0-TEMP3)/TEMP1
          CCHLX=0.3D0 * 0.083D0 * PHIMX * XKC * IAVSG / (GITMAX * 2.718D0)
!          print *, 'CCHLX ', CCHLX

          IF (CCHLX.LT.15.0D0) THEN
              CCHLX=15.0D0
!              write(6,*)PHIMX,XKC,IAVSG,GITMAX
          END IF


      END SUBROUTINE

!********************************************************************
!********************************************************************
!Function to calculate the volatilization rate of unionized ammonia
subroutine AMMONIA_VOLATILIZATION(AMMONIA_VOLATIL_RATE,NH4N, pH, TEMP, KA,nkn)

! Output:
!    AMMONIA_VOLATIL_RATE
! Inputs:
!    NH4N,
!    pH,
!    TEMP,
!    KA,
!    nkn

    implicit none
    integer nkn

    double precision, intent(in) :: NH4N(nkn)      !NH4N (mg/L)
    double precision, intent(in) :: pH  (nkn)      !pH
    double precision, intent(in) :: TEMP(nkn)      !Temperature in Celcisus
    double precision, intent(in) :: KA  (nkn)      !Reaeration rate

    double precision :: AMMONIA_VOLATIL_RATE(nkn)  !Ammonia volatilization rate (mgN/L.day)
    double precision :: NH3N (nkn)                 !Concentration of unionized ammonia (mgN/L)
    double precision :: NH3S                       !Satuation concentration of unionized ammonia in water (mgNH3N/L)

    !Declaration of function name
    !double precision :: UNIONIZED_AMMONIA

    NH3S = 0.0D0                              !Taken zero for a while assuming that the partial pressure of
                                              !unionized ammonia is zero in the atmosphere (unpolluted air)
     call UNIONIZED_AMMONIA(NH3N, NH4N, pH, TEMP,nkn)

    !1.4D1/1.70D1 : Ratio of NH3N:NH3
    !3.2D1/1.7D1  : Ratio of molecular weight of oxygen to ammonia
    AMMONIA_VOLATIL_RATE = KA * ((NH3S * (1.4D1/1.70D1)) - NH3N) * ((3.2D1/1.7D1)**2.5D-1)
end subroutine AMMONIA_VOLATILIZATION


!Function to calculate the concentration of unionized ammonia
subroutine UNIONIZED_AMMONIA(NH3N, NH4N, pH, TEMP,nkn)
! Output:
!     NH3N
! Inputs:
!   NH4N,
!   pH,
!   TEMP,
!   nkn

    implicit none
    integer nkn

    double precision NH3N(nkn)     !Concentration of unionized ammonia (mgN/L)

    double precision, intent(in) :: NH4N(nkn)      !NH4N (mg/L)
    double precision, intent(in) :: pH  (nkn)      !pH
    double precision, intent(in) :: TEMP(nkn)      !Temperature in Celcisus

    double precision :: FRAC_NH3(nkn)              !Fraction of unionized ammonia
    double precision :: pKH     (nkn)              !Concentration of unionized ammonia (mgN/L)
    double precision :: T_KELVIN(nkn)              !Temperature in Kelvins

    T_KELVIN = TEMP + 2.7316D2
    pKH      = 9.018D-2 + (2.72992D3 / T_KELVIN)
    FRAC_NH3 = 1.0D0 / (1.0D0 + (1.0D1 ** (pKH - pH)))
    NH3N     = FRAC_NH3 * NH4N
end subroutine UNIONIZED_AMMONIA

integer function STRANGERSD(VALUE,VALUE_strange,nkn)

    ! Cheks for NaN and Inf in 1D array with nkn elements
    ! Input is double precision!
      implicit none

      integer nkn

      double precision VALUE(nkn), BIGNUMBER, RATIO
      double precision value_left, value_right
      logical VALUE_NaN(nkn)
      logical VALUE_Inf(nkn)
      logical VALUE_strange(nkn)
      logical VALUE_out(nkn)

      BIGNUMBER=1.0D30
      STRANGERSD=0
      value_left = -10000.D0
      value_right=  10000.D0

      VALUE_out     = (VALUE < value_left) .or. (VALUE > value_right)
      VALUE_NAN     = isnan(VALUE)
      VALUE_Inf     = abs(VALUE) > BIGNUMBER
      VALUE_strange = VALUE_NAN .or. VALUE_Inf .or. VALUE_out

      if(any(VALUE_strange)) STRANGERSD = 1
      return

end function STRANGERSD

!************************************************************************
!************************************************************************