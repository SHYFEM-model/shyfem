! Auxilary routines for the pelagic model

! Contents:
!subroutine DIATOMS
!subroutine CYANOBACTERIA
!subroutine OTHER_PLANKTONIC_ALGAE
!subroutine FIX_CYANOBACTERIA
!subroutine ZOOPLANKTON
!subroutine CHEMOAUTOTROPHS_1
!subroutine HETEROTROPHS_1
!subroutine HETEROTROPHS_2
!subroutine ORGANIC_CARBON_DISSOLUTION
!subroutine ORGANIC_CARBON_MINERALIZATION


subroutine DIATOMS(KG_DIA_OPT_TEMP         , &
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

    use AQUABC_II_GLOBAL
    implicit none

    ! -------------------------------------------------------------------------
    ! Ingoing variables
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), intent(in) :: KG_DIA_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: DIA_OPT_TEMP_LR
    real(kind = DBL_PREC), intent(in) :: DIA_OPT_TEMP_UR
    real(kind = DBL_PREC), intent(in) :: EFF_DIA_GROWTH
    real(kind = DBL_PREC), intent(in) :: KAPPA_DIA_UNDER_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: KAPPA_DIA_OVER_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: KR_DIA_20
    real(kind = DBL_PREC), intent(in) :: THETA_KR_DIA
    real(kind = DBL_PREC), intent(in) :: KD_DIA_20
    real(kind = DBL_PREC), intent(in) :: THETA_KD_DIA
    real(kind = DBL_PREC), intent(in) :: KHS_DIN_DIA
    real(kind = DBL_PREC), intent(in) :: KHS_DIP_DIA
    real(kind = DBL_PREC), intent(in) :: KHS_DSi_DIA
    real(kind = DBL_PREC), intent(in) :: KHS_O2_DIA
    real(kind = DBL_PREC), intent(in) :: FRAC_DIA_EXCR
    real(kind = DBL_PREC), intent(in) :: I_S_DIA
    real(kind = DBL_PREC), intent(in) :: DO_STR_HYPOX_DIA_D
    real(kind = DBL_PREC), intent(in) :: THETA_HYPOX_DIA_D
    real(kind = DBL_PREC), intent(in) :: EXPON_HYPOX_DIA_D
    real(kind = DBL_PREC), intent(in) :: DIA_N_TO_C
    real(kind = DBL_PREC), intent(in) :: DIA_P_TO_C
    real(kind = DBL_PREC), intent(in) :: DIA_Si_TO_C
    real(kind = DBL_PREC), intent(in) :: DIA_O2_TO_C
    real(kind = DBL_PREC), intent(in) :: DIA_C_TO_CHLA

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NH4_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NO3_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PO4_P
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_OXYGEN
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DIA_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: ZOO_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_Si

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: I_A
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: K_E
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DEPTH
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: CHLA
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: K_B_E

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FDAY(nkn)
    real(kind = DBL_PREC), intent(in) :: TIME_STEP

    integer, intent(in) :: SMITH
    integer, intent(in) :: nkn
    ! -------------------------------------------------------------------------
    ! End of ingoing variables
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Outgoing variables
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_DIA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: ALPHA_0
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: ALPHA_1
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_DIA_TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_DIA_LIGHT
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_DIA_DOXY
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_DIA_N
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_DIA_P
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_DIA_DISS_Si
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_DIA_NUTR
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_DIA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_DIA_GROWTH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_DIA_MET
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_DIA_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_DIA_EXCR
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_DIA_INT_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KD_DIA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FAC_HYPOX_DIA_D
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_DIA_DEATH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PREF_NH4N_DIA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: DIA_LIGHT_SAT
    ! -------------------------------------------------------------------------
    ! End of outgoing variables
    ! -------------------------------------------------------------------------

     call GROWTH_AT_TEMP &
          (TEMP,KG_DIA, DIA_OPT_TEMP_LR, DIA_OPT_TEMP_UR, KG_DIA_OPT_TEMP,  &
           KAPPA_DIA_UNDER_OPT_TEMP, KAPPA_DIA_OVER_OPT_TEMP,nkn)

     if (smith .eq. 0) then
         !May be replaced by Smith formulation
         ALPHA_0 = (I_A / I_S_DIA) * exp(-1.0D0 * K_E * 0.0D0)
         ALPHA_1 = (I_A / I_S_DIA) * exp(-1.0D0 * K_E * DEPTH)

         LIM_KG_DIA_LIGHT = &
               (((2.718 * FDAY) / (K_E * DEPTH)) * &
                (exp(-1.0D0 * ALPHA_1) - exp(-1.0D0 * ALPHA_0)))
     end if

     if (smith .eq. 1) then
         call LIM_LIGHT(I_A, CHLA, KG_DIA, DEPTH, K_B_E, LIM_KG_DIA_LIGHT, &
                        DIA_C_TO_CHLA, DIA_LIGHT_SAT,nkn)
     end if

     if( KG_DIA_OPT_TEMP .ne.  0.D0) then
         LIM_KG_DIA_TEMP = KG_DIA / KG_DIA_OPT_TEMP
     else
         LIM_KG_DIA_TEMP = 0.0D0
     end if

     LIM_KG_DIA_DOXY    = DISS_OXYGEN / (KHS_O2_DIA + DISS_OXYGEN)
     LIM_KG_DIA_N       = (NH4_N + NO3_N) / (KHS_DIN_DIA + NH4_N + NO3_N)
     LIM_KG_DIA_P       = PO4_P   / (KHS_DIP_DIA + PO4_P)
     LIM_KG_DIA_DISS_Si = DISS_Si / (KHS_DSi_DIA + DISS_Si)
     LIM_KG_DIA_NUTR    = min(LIM_KG_DIA_N, LIM_KG_DIA_P, LIM_KG_DIA_DISS_Si)
     LIM_KG_DIA         = LIM_KG_DIA_LIGHT*min(LIM_KG_DIA_DOXY, LIM_KG_DIA_NUTR)

     !Diatom photo growth rate
     R_DIA_GROWTH       = KG_DIA * LIM_KG_DIA * DIA_C

    !Diatom  metabolic rate. Formulation of division metabolic rate to excretion
    !and respiration is correct only for high oxygen conc. Fixme

    R_DIA_MET     = R_DIA_GROWTH * (1.0D0 - EFF_DIA_GROWTH)
    R_DIA_RESP = (1.D0-FRAC_DIA_EXCR) * R_DIA_MET
    R_DIA_EXCR = FRAC_DIA_EXCR * R_DIA_MET
    !Diatom basal respiration rate
    R_DIA_INT_RESP = KR_DIA_20 * (THETA_KR_DIA ** (TEMP - 2.0D1)) * LIM_KG_DIA_DOXY * DIA_C

    !Calculations for diatom death
    KD_DIA = KD_DIA_20 * (THETA_KD_DIA ** (TEMP - 2.0D1))

    where (DISS_OXYGEN <= DO_STR_HYPOX_DIA_D)

        where (DISS_OXYGEN / DO_STR_HYPOX_DIA_D > 1.0D-1)
            FAC_HYPOX_DIA_D = THETA_HYPOX_DIA_D ** &
                (EXPON_HYPOX_DIA_D * &
                    (DO_STR_HYPOX_DIA_D - DISS_OXYGEN))
        elsewhere
            FAC_HYPOX_DIA_D = TIME_STEP / (5.0D-1 * KD_DIA)
            R_DIA_INT_RESP = 0.0D0
            R_DIA_RESP     = 0.0D0
            R_DIA_GROWTH   = 0.0D0
        end where
    else where
        FAC_HYPOX_DIA_D = 1.0D0
    end where

    R_DIA_DEATH = KD_DIA * FAC_HYPOX_DIA_D * DIA_C

    call AMMONIA_PREFS(PREF_NH4N_DIA,NH4_N, NO3_N, KHS_DIN_DIA,nkn)
end subroutine DIATOMS




subroutine CYANOBACTERIA &
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

    use AQUABC_II_GLOBAL
    implicit none

    ! -------------------------------------------------------------------------
    ! Ingoing variables
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), intent(in) :: KG_CYN_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: CYN_OPT_TEMP_LR
    real(kind = DBL_PREC), intent(in) :: CYN_OPT_TEMP_UR
    real(kind = DBL_PREC), intent(in) :: EFF_CYN_GROWTH
    real(kind = DBL_PREC), intent(in) :: KAPPA_CYN_UNDER_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: KAPPA_CYN_OVER_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: KR_CYN_20
    real(kind = DBL_PREC), intent(in) :: THETA_KR_CYN
    real(kind = DBL_PREC), intent(in) :: KD_CYN_20
    real(kind = DBL_PREC), intent(in) :: THETA_KD_CYN
    real(kind = DBL_PREC), intent(in) :: KHS_DIN_CYN
    real(kind = DBL_PREC), intent(in) :: KHS_DIP_CYN
    real(kind = DBL_PREC), intent(in) :: KHS_O2_CYN
    real(kind = DBL_PREC), intent(in) :: FRAC_CYN_EXCR
    real(kind = DBL_PREC), intent(in) :: I_S_CYN
    real(kind = DBL_PREC), intent(in) :: DO_STR_HYPOX_CYN_D
    real(kind = DBL_PREC), intent(in) :: THETA_HYPOX_CYN_D
    real(kind = DBL_PREC), intent(in) :: EXPON_HYPOX_CYN_D
    real(kind = DBL_PREC), intent(in) :: CYN_N_TO_C
    real(kind = DBL_PREC), intent(in) :: CYN_P_TO_C
    real(kind = DBL_PREC), intent(in) :: CYN_O2_TO_C
    real(kind = DBL_PREC), intent(in) :: CYN_C_TO_CHLA

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NH4_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NO3_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PO4_P
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_OXYGEN
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: CYN_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: ZOO_C

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: I_A
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: K_E
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DEPTH
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: CHLA
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: K_B_E

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FDAY(nkn)
    real(kind = DBL_PREC), intent(in) :: TIME_STEP

    integer, intent(in) :: SMITH
    integer, intent(in) :: nkn
    ! -------------------------------------------------------------------------
    ! End of ingoing variables
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Outgoing variables
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: ALPHA_0
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: ALPHA_1
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_CYN_TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_CYN_LIGHT
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_CYN_DOXY
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_CYN_N
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_CYN_P
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_CYN_NUTR
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_CYN_GROWTH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_CYN_MET
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_CYN_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_CYN_EXCR
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_CYN_INT_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KD_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FAC_HYPOX_CYN_D
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_CYN_DEATH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PREF_NH4N_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: CYN_LIGHT_SAT
    ! -------------------------------------------------------------------------
    ! End of outgoing variables
    ! -------------------------------------------------------------------------

    call GROWTH_AT_TEMP &
         (TEMP,KG_CYN, CYN_OPT_TEMP_LR, CYN_OPT_TEMP_UR, KG_CYN_OPT_TEMP,  &
          KAPPA_CYN_UNDER_OPT_TEMP, KAPPA_CYN_OVER_OPT_TEMP,nkn)

    if (smith .eq. 0) then
        !May be replaced by Smith formulation
        ALPHA_0 = (I_A / I_S_CYN) * exp(-1.0D0 * K_E * 0.0D0)
        ALPHA_1 = (I_A / I_S_CYN) * exp(-1.0D0 * K_E * DEPTH)

        LIM_KG_CYN_LIGHT = &
            (((2.718 * FDAY) / (K_E * DEPTH)) * &
             (exp(-1.0D0 * ALPHA_1) - exp(-1.0D0 * ALPHA_0)))
    end if

    if (smith .eq. 1) then
        call LIM_LIGHT(I_A, CHLA, KG_CYN, DEPTH, K_B_E, LIM_KG_CYN_LIGHT, CYN_C_TO_CHLA, CYN_LIGHT_SAT,nkn)
    end if


    if( KG_CYN_OPT_TEMP .ne.  0.D0) then
       LIM_KG_CYN_TEMP = KG_CYN / KG_CYN_OPT_TEMP
    else
       LIM_KG_CYN_TEMP = 0.D0
    end if

    LIM_KG_CYN_DOXY = DISS_OXYGEN / (KHS_O2_CYN + DISS_OXYGEN)
    LIM_KG_CYN_N    = (NH4_N + NO3_N) / (KHS_DIN_CYN + NH4_N + NO3_N)
    LIM_KG_CYN_P    = PO4_P   / (KHS_DIP_CYN + PO4_P)
    LIM_KG_CYN_NUTR = min(LIM_KG_CYN_N, LIM_KG_CYN_P)
    !LIM_KG_CYN      = min(LIM_KG_CYN_DOXY, LIM_KG_CYN_NUTR, LIM_KG_CYN_LIGHT)
    LIM_KG_CYN      = LIM_KG_CYN_LIGHT*min(LIM_KG_CYN_DOXY, LIM_KG_CYN_NUTR) !changed by Petras
    !Non-fixing cyanobacteria growth rate
    R_CYN_GROWTH    = KG_CYN * LIM_KG_CYN * CYN_C

    !Non-fixing cyanobacteria metabolic rate, respiration and excretion
    R_CYN_MET = R_CYN_GROWTH * (1.0D0 - EFF_CYN_GROWTH)
    R_CYN_RESP = (1.D0-FRAC_CYN_EXCR) * R_CYN_MET
    R_CYN_EXCR = FRAC_CYN_EXCR * R_CYN_MET

    !Non-fixing cyanobacteria dark respiration rate
    R_CYN_INT_RESP = KR_CYN_20 * (THETA_KR_CYN ** (TEMP - 2.0D1)) * LIM_KG_CYN_DOXY * CYN_C

    !Calculations for non-fixing cyanobacteria death rate
    KD_CYN = KD_CYN_20 * (THETA_KD_CYN ** (TEMP - 2.0D1))

    where (DISS_OXYGEN <= DO_STR_HYPOX_CYN_D)

        where (DISS_OXYGEN / DO_STR_HYPOX_CYN_D > 1.0D-1)
            FAC_HYPOX_CYN_D = THETA_HYPOX_CYN_D ** &
            &    (EXPON_HYPOX_CYN_D * (DO_STR_HYPOX_CYN_D - DISS_OXYGEN))
        elsewhere
            FAC_HYPOX_CYN_D = TIME_STEP / (5.0D-1 * KD_CYN)
            R_CYN_INT_RESP = 0.0D0
            R_CYN_RESP     = 0.0D0
            R_CYN_GROWTH   = 0.0D0
        end where
    elsewhere
        FAC_HYPOX_CYN_D = 1.0D0
    end where

    !Non-fixing cyanobacteria death rate
    R_CYN_DEATH = KD_CYN * FAC_HYPOX_CYN_D * CYN_C

    !PREF_NH4N_CYN = NH4_N / (NH4_N + KHS_NH4N_PREF_CYN)
    call AMMONIA_PREFS(PREF_NH4N_CYN, NH4_N, NO3_N, KHS_DIN_CYN,nkn)
end subroutine CYANOBACTERIA




subroutine OTHER_PLANKTONIC_ALGAE &
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

    use AQUABC_II_GLOBAL
    implicit none

    ! -------------------------------------------------------------------------
    ! Ingoing variables
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), intent(in) :: KG_OPA_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: OPA_OPT_TEMP_LR
    real(kind = DBL_PREC), intent(in) :: OPA_OPT_TEMP_UR
    real(kind = DBL_PREC), intent(in) :: EFF_OPA_GROWTH
    real(kind = DBL_PREC), intent(in) :: KAPPA_OPA_UNDER_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: KAPPA_OPA_OVER_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: KR_OPA_20
    real(kind = DBL_PREC), intent(in) :: THETA_KR_OPA
    real(kind = DBL_PREC), intent(in) :: KD_OPA_20
    real(kind = DBL_PREC), intent(in) :: THETA_KD_OPA
    real(kind = DBL_PREC), intent(in) :: KHS_DIN_OPA
    real(kind = DBL_PREC), intent(in) :: KHS_DIP_OPA
    real(kind = DBL_PREC), intent(in) :: KHS_O2_OPA
    real(kind = DBL_PREC), intent(in) :: FRAC_OPA_EXCR
    real(kind = DBL_PREC), intent(in) :: I_S_OPA
    real(kind = DBL_PREC), intent(in) :: DO_STR_HYPOX_OPA_D
    real(kind = DBL_PREC), intent(in) :: THETA_HYPOX_OPA_D
    real(kind = DBL_PREC), intent(in) :: EXPON_HYPOX_OPA_D
    real(kind = DBL_PREC), intent(in) :: OPA_N_TO_C
    real(kind = DBL_PREC), intent(in) :: OPA_P_TO_C
    real(kind = DBL_PREC), intent(in) :: OPA_O2_TO_C
    real(kind = DBL_PREC), intent(in) :: OPA_C_TO_CHLA

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NH4_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NO3_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PO4_P
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_OXYGEN
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: OPA_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: ZOO_C

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: I_A
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: K_E
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DEPTH
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: CHLA
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: K_B_E

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FDAY(nkn)
    real(kind = DBL_PREC), intent(in) :: TIME_STEP

    integer, intent(in) :: SMITH
    integer, intent(in) :: nkn
    ! -------------------------------------------------------------------------
    ! End of ingoing variables
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Outgoing variables
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_OPA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: ALPHA_0
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: ALPHA_1
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_OPA_TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_OPA_LIGHT
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_OPA_DOXY
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_OPA_N
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_OPA_P
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_OPA_NUTR
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_OPA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_OPA_GROWTH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_OPA_MET
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_OPA_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_OPA_EXCR
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_OPA_INT_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KD_OPA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FAC_HYPOX_OPA_D
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_OPA_DEATH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PREF_NH4N_OPA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: OPA_LIGHT_SAT
    ! -------------------------------------------------------------------------
    ! End of outgoing variables
    ! -------------------------------------------------------------------------

    ALPHA_0 = (I_A / I_S_OPA) * exp(-1.0D0 * K_E * 0.0D0)
    ALPHA_1 = (I_A / I_S_OPA) * exp(-1.0D0 * K_E * DEPTH)

    call GROWTH_AT_TEMP &
         (TEMP,KG_OPA, OPA_OPT_TEMP_LR, OPA_OPT_TEMP_UR, KG_OPA_OPT_TEMP,  &
          KAPPA_OPA_UNDER_OPT_TEMP, KAPPA_OPA_OVER_OPT_TEMP,nkn)

    if (smith .eq. 0) then
        !May be replaced by Smith formulation
        LIM_KG_OPA_LIGHT = &
             (((2.718 * FDAY) / (K_E * DEPTH)) * &
              (exp(-1.0D0 * ALPHA_1) - exp(-1.0D0 * ALPHA_0)))
    end if

    if (smith .eq. 1) then
        call LIM_LIGHT(I_A, CHLA, KG_OPA, DEPTH, K_B_E, LIM_KG_OPA_LIGHT, OPA_C_TO_CHLA, OPA_LIGHT_SAT,nkn)
    end if

    if( KG_OPA_OPT_TEMP .ne.  0.D0) then
        LIM_KG_OPA_TEMP = KG_OPA / KG_OPA_OPT_TEMP
    else
        LIM_KG_OPA_TEMP = 0.D0
    end if

    LIM_KG_OPA_DOXY = DISS_OXYGEN / (KHS_O2_OPA + DISS_OXYGEN)
    LIM_KG_OPA_N    = (NH4_N + NO3_N) / (KHS_DIN_OPA + NH4_N + NO3_N)
    LIM_KG_OPA_P    = PO4_P   / (KHS_DIP_OPA + PO4_P)
    LIM_KG_OPA_NUTR = min(LIM_KG_OPA_N, LIM_KG_OPA_P)
    !    LIM_KG_OPA      = min(LIM_KG_OPA_DOXY, LIM_KG_OPA_NUTR, LIM_KG_OPA_LIGHT)
    LIM_KG_OPA      = LIM_KG_OPA_LIGHT * min(LIM_KG_OPA_DOXY, LIM_KG_OPA_NUTR) !changed by Petras

    !Other planktonic algae growth rate
    R_OPA_GROWTH    = KG_OPA * LIM_KG_OPA * OPA_C

    !Other planktonic algae metabolism, respiration, excretion rate

    R_OPA_MET    = R_OPA_GROWTH * (1.0D0 - EFF_OPA_GROWTH)
    R_OPA_RESP  =  (1.D0-FRAC_OPA_EXCR) * R_OPA_MET
    R_OPA_EXCR = FRAC_OPA_EXCR * R_OPA_MET

    !Other planktonic algae dark respiration rate
    R_OPA_INT_RESP  = KR_OPA_20 * (THETA_KR_OPA ** (TEMP - 2.0D1)) * LIM_KG_OPA_DOXY * OPA_C

    !Calculations for other planktonic algae death
    KD_OPA = KD_OPA_20 * (THETA_KD_OPA ** (TEMP - 2.0D1))

    where (DISS_OXYGEN <= DO_STR_HYPOX_OPA_D)

        where (DISS_OXYGEN / DO_STR_HYPOX_OPA_D > 1.0D-1)
            FAC_HYPOX_OPA_D = THETA_HYPOX_OPA_D ** &
                (EXPON_HYPOX_OPA_D * &
                    (DO_STR_HYPOX_OPA_D - DISS_OXYGEN))
        elsewhere
            FAC_HYPOX_OPA_D = TIME_STEP / (5.0D-1 * KD_OPA)
            R_OPA_INT_RESP  = 0.0D0
            R_OPA_RESP      = 0.0D0
            R_OPA_GROWTH    = 0.0D0
        end where
    elsewhere
        FAC_HYPOX_OPA_D = 1.0D0
    end where

    !Other planktonic algae death rate
    R_OPA_DEATH = KD_OPA * FAC_HYPOX_OPA_D * OPA_C

    !PREF_NH4N_OPA = NH4_N / (NH4_N + KHS_NH4N_PREF_OPA)
    call AMMONIA_PREFS(PREF_NH4N_OPA, NH4_N, NO3_N, KHS_DIN_OPA,nkn)
end subroutine OTHER_PLANKTONIC_ALGAE



subroutine FIX_CYANOBACTERIA  &
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

    use AQUABC_II_GLOBAL
    implicit none

    ! -------------------------------------------------------------------------
    ! Ingoing variables
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), intent(in) :: KG_FIX_CYN_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: FIX_CYN_OPT_TEMP_LR
    real(kind = DBL_PREC), intent(in) :: FIX_CYN_OPT_TEMP_UR
    real(kind = DBL_PREC), intent(in) :: EFF_FIX_CYN_GROWTH
    real(kind = DBL_PREC), intent(in) :: KAPPA_FIX_CYN_UNDER_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: KAPPA_FIX_CYN_OVER_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: KR_FIX_CYN_20
    real(kind = DBL_PREC), intent(in) :: THETA_KR_FIX_CYN
    real(kind = DBL_PREC), intent(in) :: KD_FIX_CYN_20
    real(kind = DBL_PREC), intent(in) :: THETA_KD_FIX_CYN
    real(kind = DBL_PREC), intent(in) :: KHS_DIN_FIX_CYN
    real(kind = DBL_PREC), intent(in) :: KHS_DIP_FIX_CYN
    real(kind = DBL_PREC), intent(in) :: KHS_O2_FIX_CYN
    real(kind = DBL_PREC), intent(in) :: KHS_NH4N_PREF_FIX_CYN
    real(kind = DBL_PREC), intent(in) :: I_S_FIX_CYN
    real(kind = DBL_PREC), intent(in) :: DO_STR_HYPOX_FIX_CYN_D
    real(kind = DBL_PREC), intent(in) :: THETA_HYPOX_FIX_CYN_D
    real(kind = DBL_PREC), intent(in) :: EXPON_HYPOX_FIX_CYN_D
    real(kind = DBL_PREC), intent(in) :: FIX_CYN_N_TO_C
    real(kind = DBL_PREC), intent(in) :: FIX_CYN_P_TO_C
    real(kind = DBL_PREC), intent(in) :: FIX_CYN_O2_TO_C
    real(kind = DBL_PREC), intent(in) :: FIX_CYN_C_TO_CHLA
    real(kind = DBL_PREC), intent(in) :: FIX_CYN_C_TO_CHLA_NEW
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FIX_CYN_LIGHT_SAT
    real(kind = DBL_PREC), intent(in) :: FRAC_FIX_CYN_EXCR
    real(kind = DBL_PREC), intent(in) :: R_FIX
    real(kind = DBL_PREC), intent(in) :: K_FIX

    real(kind = DBL_PREC), intent(in) :: TIME_STEP
    integer, intent(in) :: SMITH
    integer, intent(in) :: nkn
!
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FDAY
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: I_A
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: K_E
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: K_B_E
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DEPTH
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: CHLA
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NH4_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NO3_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PO4_P
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_OXYGEN
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FIX_CYN_C
    ! -------------------------------------------------------------------------
    ! End of ingoing variables
    ! -------------------------------------------------------------------------

    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: ALPHA_0
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: ALPHA_1
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_FIX_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_FIX_CYN_LIGHT
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_FIX_CYN_TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_FIX_CYN_DOXY
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_NON_FIX_CYN_N
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_NON_FIX_CYN_P
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_NON_FIX_CYN_NUTR
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_FIX_FIX_CYN_N
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_FIX_FIX_CYN_P
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_FIX_FIX_CYN_NUTR
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_NON_FIX_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_KG_FIX_FIX_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_NON_FIX_CYN_GROWTH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_FIX_FIX_CYN_GROWTH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_FIX_CYN_GROWTH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_FIX_CYN_MET
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_FIX_CYN_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_FIX_CYN_EXCR
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_FIX_CYN_INT_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KD_FIX_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FAC_HYPOX_FIX_CYN_D
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_FIX_CYN_DEATH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PREF_NH4N_FIX_CYN

    !Auxillary variable
    real(kind = DBL_PREC) :: FIX_CYN_DEPTH

    !Caculations for nitrogen fixing cyanobacteria growth
    call GROWTH_AT_TEMP &
         (TEMP,KG_FIX_CYN, FIX_CYN_OPT_TEMP_LR, FIX_CYN_OPT_TEMP_UR, KG_FIX_CYN_OPT_TEMP,  &
          KAPPA_FIX_CYN_UNDER_OPT_TEMP, KAPPA_FIX_CYN_OVER_OPT_TEMP,nkn)

    if (smith .eq. 0) then
        !May be replaced by Smith formulation
        ALPHA_0 = (I_A / I_S_FIX_CYN) * exp(-1.0D0 * K_E * 0.0D0)
        ALPHA_1 = (I_A / I_S_FIX_CYN) * exp(-1.0D0 * K_E * DEPTH)

        LIM_KG_FIX_CYN_LIGHT = &
            (((2.718 * FDAY) / (K_E * DEPTH)) * &
             (exp(-1.0D0 * ALPHA_1) - exp(-1.0D0 * ALPHA_0)))
        !WC_OUTPUTS(nstate+3) = FIX_CYN_C_TO_CHLA
    end if

    if (smith .eq. 1) then
    	 FIX_CYN_DEPTH = 1.  !1.2 is assumed that all fixers are in the layer of this depth (Introduced 2013 working with Ali)
    	                     ! Canged to 1 by Petras 2014 10 13
         call LIM_LIGHT(I_A, CHLA, KG_FIX_CYN, DEPTH, K_B_E, LIM_KG_FIX_CYN_LIGHT, FIX_CYN_C_TO_CHLA, FIX_CYN_LIGHT_SAT,nkn)
         LIM_KG_FIX_CYN_LIGHT = FIX_CYN_DEPTH*LIM_KG_FIX_CYN_LIGHT

         !WC_OUTPUTS(nstate+3) = FIX_CYN_C_TO_CHLA
    end if

    if( KG_FIX_CYN_OPT_TEMP .ne.  0.D0) then
        LIM_KG_FIX_CYN_TEMP = KG_FIX_CYN / KG_FIX_CYN_OPT_TEMP
    else
        LIM_KG_FIX_CYN_TEMP = 0.D0
    end if

    LIM_KG_FIX_CYN_DOXY     = DISS_OXYGEN / (KHS_O2_FIX_CYN + DISS_OXYGEN)

    !Nutrient limitation of fixing cyanobacteria in non-fixing fraction
    LIM_KG_NON_FIX_CYN_N    = (NH4_N + NO3_N) / (KHS_DIN_FIX_CYN + NH4_N + NO3_N)
    LIM_KG_NON_FIX_CYN_P    = PO4_P / (KHS_DIP_FIX_CYN + PO4_P)
    LIM_KG_NON_FIX_CYN_NUTR = min(LIM_KG_NON_FIX_CYN_N, LIM_KG_NON_FIX_CYN_P)

    !Nutrient limitation of fixing cyanobacteria in fixing fraction
    LIM_KG_FIX_FIX_CYN_N    = (K_FIX / (K_FIX + NH4_N + NO3_N))
    LIM_KG_FIX_FIX_CYN_P    = LIM_KG_NON_FIX_CYN_P
    LIM_KG_FIX_FIX_CYN_NUTR = min(LIM_KG_FIX_FIX_CYN_N, LIM_KG_FIX_FIX_CYN_P)

    !Growth limitation of fixing cyanobacteria in non-fixing fraction
    LIM_KG_NON_FIX_CYN   = LIM_KG_FIX_CYN_LIGHT*min(LIM_KG_FIX_CYN_DOXY, LIM_KG_NON_FIX_CYN_NUTR)

    !Growth limitation of fixing cyanobacteria in fixing fraction
    LIM_KG_FIX_FIX_CYN   = LIM_KG_FIX_CYN_LIGHT* min(LIM_KG_FIX_CYN_DOXY, LIM_KG_FIX_FIX_CYN_NUTR)

    !Growth rate of fixing cyanobacteria in non-fixing fraction
    R_NON_FIX_CYN_GROWTH = KG_FIX_CYN * LIM_KG_NON_FIX_CYN * FIX_CYN_C

    !Growth rate of fixing cyanobacteria in fixing state
    R_FIX_FIX_CYN_GROWTH = R_FIX * KG_FIX_CYN * LIM_KG_FIX_FIX_CYN * FIX_CYN_C

    !Total growth rate of fixing cyanobacteria as a sum of non-fixing and
    !fixing fractions.
    R_FIX_CYN_GROWTH = R_NON_FIX_CYN_GROWTH + R_FIX_FIX_CYN_GROWTH

    !Nitrogen fixing cyanobacteria metabolism, respiration, excretion rate
    R_FIX_CYN_MET = R_FIX_CYN_GROWTH * (1.0D0 - EFF_FIX_CYN_GROWTH)
    R_FIX_CYN_RESP = (1.D0-FRAC_FIX_CYN_EXCR) * R_FIX_CYN_MET
    R_FIX_CYN_EXCR = FRAC_FIX_CYN_EXCR * R_FIX_CYN_MET

    !Nitrogen fixing cyanobacteria dark respiration rate
    R_FIX_CYN_INT_RESP = KR_FIX_CYN_20 * (THETA_KR_FIX_CYN ** (TEMP - 2.0D1)) * &
           LIM_KG_FIX_CYN_DOXY * FIX_CYN_C

    !Nitrogen fixing cyanobacteria death rate
    KD_FIX_CYN = KD_FIX_CYN_20 * (THETA_KD_FIX_CYN ** (TEMP - 2.0D1))

    where (DISS_OXYGEN <= DO_STR_HYPOX_FIX_CYN_D)

        where (DISS_OXYGEN / DO_STR_HYPOX_FIX_CYN_D > 1.0D-1)
            FAC_HYPOX_FIX_CYN_D = THETA_HYPOX_FIX_CYN_D ** &
            &    (EXPON_HYPOX_FIX_CYN_D * (DO_STR_HYPOX_FIX_CYN_D - DISS_OXYGEN))
        elsewhere
            FAC_HYPOX_FIX_CYN_D = TIME_STEP / (5.0D-1 * KD_FIX_CYN)
            R_FIX_CYN_INT_RESP = 0.0D0
            R_FIX_CYN_RESP     = 0.0D0
            R_FIX_CYN_GROWTH   = 0.0D0
        end where
    elsewhere
        FAC_HYPOX_FIX_CYN_D = 1.0D0
    end where

    !Nitrogen fixing cyanobacteria death rate
    R_FIX_CYN_DEATH = KD_FIX_CYN * FAC_HYPOX_FIX_CYN_D * FIX_CYN_C

    !PREF_NH4N_FIX_CYN = NH4_N / (NH4_N + KHS_NH4N_PREF_FIX_CYN)
    call AMMONIA_PREFS(PREF_NH4N_FIX_CYN, NH4_N, NO3_N, KHS_DIN_FIX_CYN,nkn)
end subroutine FIX_CYANOBACTERIA



subroutine ZOOPLANKTON &
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

    use AQUABC_II_GLOBAL
    implicit none

    ! -------------------------------------------------------------------------
    ! Ingoing variables
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), intent(in) :: KG_ZOO_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: ZOO_OPT_TEMP_LR
    real(kind = DBL_PREC), intent(in) :: ZOO_OPT_TEMP_UR
    real(kind = DBL_PREC), intent(in) :: EFF_ZOO_GROWTH
    real(kind = DBL_PREC), intent(in) :: KAPPA_ZOO_UNDER_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: KAPPA_ZOO_OVER_OPT_TEMP
    real(kind = DBL_PREC), intent(in) :: GRAT_ZOO_DIA
    real(kind = DBL_PREC), intent(in) :: GRAT_ZOO_CYN
    real(kind = DBL_PREC), intent(in) :: GRAT_ZOO_OPA
    real(kind = DBL_PREC), intent(in) :: GRAT_ZOO_FIX_CYN
    real(kind = DBL_PREC), intent(in) :: GRAT_ZOO_CHEM_AUT_BAC
    real(kind = DBL_PREC), intent(in) :: GRAT_ZOO_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: GRAT_ZOO_FAC_AN_HET_BAC
    real(kind = DBL_PREC), intent(in) :: GRAT_ZOO_DET_PART_ORG_C
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_DIA
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_CYN
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_FIX_CYN
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_OPA
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_CHEM_AUT_BAC
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_FAC_AN_HET_BAC
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_DET_PART_ORG_C
    real(kind = DBL_PREC), intent(in) :: KHS_DIA_C_ZOO
    real(kind = DBL_PREC), intent(in) :: KHS_CYN_C_ZOO
    real(kind = DBL_PREC), intent(in) :: KHS_FIX_CYN_C_ZOO
    real(kind = DBL_PREC), intent(in) :: KHS_OPA_C_ZOO
    real(kind = DBL_PREC), intent(in) :: KHS_CHEM_AUT_BAC_C_ZOO
    real(kind = DBL_PREC), intent(in) :: KHS_AER_HET_BAC_C_ZOO
    real(kind = DBL_PREC), intent(in) :: KHS_FAC_AN_HET_BAC_C_ZOO
    real(kind = DBL_PREC), intent(in) :: KHS_DET_PART_ORG_C_ZOO
    real(kind = DBL_PREC), intent(in) :: FOOD_MIN_ZOO
    real(kind = DBL_PREC), intent(in) :: KE_ZOO
    real(kind = DBL_PREC), intent(in) :: FRAC_ZOO_EX_ORG
    real(kind = DBL_PREC), intent(in) :: KR_ZOO_20
    real(kind = DBL_PREC), intent(in) :: THETA_KR_ZOO
    real(kind = DBL_PREC), intent(in) :: KD_ZOO_20
    real(kind = DBL_PREC), intent(in) :: THETA_KD_ZOO
    real(kind = DBL_PREC), intent(in) :: DO_STR_HYPOX_ZOO_D
    real(kind = DBL_PREC), intent(in) :: THETA_HYPOX_ZOO_D
    real(kind = DBL_PREC), intent(in) :: EXPON_HYPOX_ZOO_D
    real(kind = DBL_PREC), intent(in) :: ZOO_N_TO_C
    real(kind = DBL_PREC), intent(in) :: ZOO_P_TO_C
    real(kind = DBL_PREC), intent(in) :: ZOO_O2_TO_C

    real(kind = DBL_PREC), intent(in) :: TIME_STEP
    integer, intent(in) :: nkn

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_OXYGEN
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DIA_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: CYN_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: OPA_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FIX_CYN_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: CHEM_AUT_BAC_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: AER_HET_BAC_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FAC_AN_HET_BAC_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DET_PART_ORG_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: ZOO_C
    ! -------------------------------------------------------------------------
    ! End of ingoing variables
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Outgoing variables
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_ZOO
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_ZOO_DIA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_ZOO_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_ZOO_OPA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_ZOO_FIX_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_ZOO_CHEM_AUT_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_ZOO_AER_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_ZOO_FAC_AN_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_ZOO_DET_PART_ORG_C
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KD_ZOO
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_DIA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_OPA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_FIX_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_CHEM_AUT_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_AER_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_FAC_AN_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_DET_PART_ORG_C
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_FEEDING_DIA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_FEEDING_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_FEEDING_FIX_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_FEEDING_OPA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_FEEDING_CHEM_AUT_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_FEEDING_AER_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_FEEDING_FAC_AN_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_FEEDING_DET_PART_ORG_C
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_INT_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_EX_DON
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_EX_DOP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_EX_DOC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_DEATH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: ACTUAL_ZOO_N_TO_C
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: ACTUAL_ZOO_P_TO_C
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_GROWTH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FAC_HYPOX_ZOO_D
    ! -------------------------------------------------------------------------
    ! End of outgoing variables
    ! -------------------------------------------------------------------------

    !Zooplankton growth rate
    call GROWTH_AT_TEMP &
         (TEMP,KG_ZOO, ZOO_OPT_TEMP_LR, ZOO_OPT_TEMP_UR, KG_ZOO_OPT_TEMP,  &
          KAPPA_ZOO_UNDER_OPT_TEMP, KAPPA_ZOO_OVER_OPT_TEMP,nkn)

    KG_ZOO_DIA            = KG_ZOO * GRAT_ZOO_DIA
    KG_ZOO_CYN            = KG_ZOO * GRAT_ZOO_CYN
    KG_ZOO_OPA            = KG_ZOO * GRAT_ZOO_OPA
    KG_ZOO_FIX_CYN        = KG_ZOO * GRAT_ZOO_FIX_CYN
    KG_ZOO_CHEM_AUT_BAC   = KG_ZOO * GRAT_ZOO_CHEM_AUT_BAC
    KG_ZOO_AER_HET_BAC    = KG_ZOO * GRAT_ZOO_AER_HET_BAC
    KG_ZOO_FAC_AN_HET_BAC = KG_ZOO * GRAT_ZOO_FAC_AN_HET_BAC
    KG_ZOO_DET_PART_ORG_C = KG_ZOO * GRAT_ZOO_DET_PART_ORG_C

    where (DIA_C > FOOD_MIN_ZOO)
        FOOD_FACTOR_ZOO_DIA = (PREF_ZOO_DIA * (DIA_C - FOOD_MIN_ZOO)) / &
            (DIA_C + KHS_DIA_C_ZOO)
    elsewhere
        FOOD_FACTOR_ZOO_DIA = 0.0D0
    end where

    where (CYN_C > FOOD_MIN_ZOO)
        FOOD_FACTOR_ZOO_CYN = (PREF_ZOO_CYN * (CYN_C - FOOD_MIN_ZOO)) / &
            (CYN_C + KHS_CYN_C_ZOO)
    elsewhere
        FOOD_FACTOR_ZOO_CYN = 0.0D0
    end where

    where (OPA_C > FOOD_MIN_ZOO)
        FOOD_FACTOR_ZOO_OPA = (PREF_ZOO_OPA * (OPA_C - FOOD_MIN_ZOO)) / &
            (OPA_C + KHS_OPA_C_ZOO)
    elsewhere
        FOOD_FACTOR_ZOO_OPA = 0.0D0
    end where

    where (FIX_CYN_C > FOOD_MIN_ZOO)
        FOOD_FACTOR_ZOO_FIX_CYN = (PREF_ZOO_FIX_CYN * (FIX_CYN_C - FOOD_MIN_ZOO)) / &
            (FIX_CYN_C + KHS_FIX_CYN_C_ZOO)
    elsewhere
        FOOD_FACTOR_ZOO_FIX_CYN = 0.0D0
    end where

    where (CHEM_AUT_BAC_C > FOOD_MIN_ZOO)
        FOOD_FACTOR_ZOO_CHEM_AUT_BAC = &
            (PREF_ZOO_CHEM_AUT_BAC * (CHEM_AUT_BAC_C - FOOD_MIN_ZOO)) / &
            (CHEM_AUT_BAC_C + KHS_CHEM_AUT_BAC_C_ZOO)
    elsewhere
        FOOD_FACTOR_ZOO_CHEM_AUT_BAC = 0.0D0
    end where

    where (AER_HET_BAC_C > FOOD_MIN_ZOO)
        FOOD_FACTOR_ZOO_AER_HET_BAC = &
            (PREF_ZOO_AER_HET_BAC * (AER_HET_BAC_C - FOOD_MIN_ZOO)) / &
            (AER_HET_BAC_C + KHS_AER_HET_BAC_C_ZOO)
    elsewhere
        FOOD_FACTOR_ZOO_AER_HET_BAC = 0.0D0
    end where

    where (FAC_AN_HET_BAC_C > FOOD_MIN_ZOO)
        FOOD_FACTOR_ZOO_FAC_AN_HET_BAC = &
            (PREF_ZOO_FAC_AN_HET_BAC * (FAC_AN_HET_BAC_C - FOOD_MIN_ZOO)) / &
            (FAC_AN_HET_BAC_C + KHS_FAC_AN_HET_BAC_C_ZOO)
    elsewhere
        FOOD_FACTOR_ZOO_FAC_AN_HET_BAC = 0.0D0
    end where

    where (DET_PART_ORG_C > FOOD_MIN_ZOO)
        FOOD_FACTOR_ZOO_DET_PART_ORG_C = &
            (PREF_ZOO_DET_PART_ORG_C * (DET_PART_ORG_C - FOOD_MIN_ZOO)) / &
            (DET_PART_ORG_C + KHS_DET_PART_ORG_C_ZOO)
    elsewhere
        FOOD_FACTOR_ZOO_DET_PART_ORG_C = 0.0D0
    end where

    R_ZOO_FEEDING_DIA            = KG_ZOO_DIA            * FOOD_FACTOR_ZOO_DIA            * ZOO_C
    R_ZOO_FEEDING_CYN            = KG_ZOO_CYN            * FOOD_FACTOR_ZOO_CYN            * ZOO_C
    R_ZOO_FEEDING_FIX_CYN        = KG_ZOO_FIX_CYN        * FOOD_FACTOR_ZOO_FIX_CYN        * ZOO_C
    R_ZOO_FEEDING_OPA            = KG_ZOO_OPA            * FOOD_FACTOR_ZOO_OPA            * ZOO_C
    R_ZOO_FEEDING_CHEM_AUT_BAC   = KG_ZOO_CHEM_AUT_BAC   * FOOD_FACTOR_ZOO_CHEM_AUT_BAC   * ZOO_C
    R_ZOO_FEEDING_AER_HET_BAC    = KG_ZOO_AER_HET_BAC    * FOOD_FACTOR_ZOO_AER_HET_BAC    * ZOO_C
    R_ZOO_FEEDING_FAC_AN_HET_BAC = KG_ZOO_FAC_AN_HET_BAC * FOOD_FACTOR_ZOO_FAC_AN_HET_BAC * ZOO_C
    R_ZOO_FEEDING_DET_PART_ORG_C = KG_ZOO_DET_PART_ORG_C * FOOD_FACTOR_ZOO_DET_PART_ORG_C * ZOO_C

    !Zooplankton excretion rate
    ACTUAL_ZOO_N_TO_C = ZOO_N_TO_C
    ACTUAL_ZOO_P_TO_C = ZOO_P_TO_C

    !Zooplankton Growth rate
    R_ZOO_GROWTH =  R_ZOO_FEEDING_DIA         + R_ZOO_FEEDING_CYN            + &
                    R_ZOO_FEEDING_OPA         + R_ZOO_FEEDING_CHEM_AUT_BAC   + &
                    R_ZOO_FEEDING_AER_HET_BAC + R_ZOO_FEEDING_FAC_AN_HET_BAC + &
                    R_ZOO_FEEDING_DET_PART_ORG_C

    R_ZOO_RESP = R_ZOO_GROWTH * (1.0D0 - EFF_ZOO_GROWTH)
    R_ZOO_INT_RESP = KR_ZOO_20 * (THETA_KR_ZOO ** (TEMP - 2.0D1)) * ZOO_C

    !Zooplankton death rate
    KD_ZOO = KD_ZOO_20 * (THETA_KD_ZOO ** (TEMP - 2.0D1))

    where (DISS_OXYGEN <= DO_STR_HYPOX_ZOO_D)

        where (DISS_OXYGEN / DO_STR_HYPOX_ZOO_D > 1.0D-1)
            FAC_HYPOX_ZOO_D = THETA_HYPOX_ZOO_D ** &
                (EXPON_HYPOX_ZOO_D * &
                    (DO_STR_HYPOX_ZOO_D - DISS_OXYGEN))
        elsewhere
            FAC_HYPOX_ZOO_D = TIME_STEP / (5.0D-1 * KD_ZOO)
            R_ZOO_FEEDING_DIA            = 0.0D0
            R_ZOO_FEEDING_CYN            = 0.0D0
            R_ZOO_FEEDING_OPA            = 0.0D0
            R_ZOO_FEEDING_CHEM_AUT_BAC   = 0.0D0
            R_ZOO_FEEDING_AER_HET_BAC    = 0.0D0
            R_ZOO_FEEDING_FAC_AN_HET_BAC = 0.0D0
            R_ZOO_FEEDING_DET_PART_ORG_C = 0.0D0
            R_ZOO_INT_RESP               = 0.0D0
            R_ZOO_RESP                   = 0.0D0
            R_ZOO_EX_DON                 = 0.0D0
            R_ZOO_EX_DOP                 = 0.0D0
            R_ZOO_EX_DOC                 = 0.0D0
        end where
    elsewhere
        FAC_HYPOX_ZOO_D = 1.0D0
    end where

    R_ZOO_DEATH = KD_ZOO * FAC_HYPOX_ZOO_D * ZOO_C
end subroutine ZOOPLANKTON



subroutine CHEMOAUTOTROPHS_1 &
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

    use AQUABC_II_GLOBAL
    implicit none

    real(kind = DBL_PREC), intent(in) :: KG_CHEM_AUT_BAC_20
    real(kind = DBL_PREC), intent(in) :: THETA_KG_CHEM_AUT_BAC
    real(kind = DBL_PREC), intent(in) :: KR_CHEM_AUT_BAC_20
    real(kind = DBL_PREC), intent(in) :: THETA_KR_CHEM_AUT_BAC
    real(kind = DBL_PREC), intent(in) :: KD_CHEM_AUT_BAC_20
    real(kind = DBL_PREC), intent(in) :: THETA_KD_CHEM_AUT_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_NH4N_CHEM_AUT_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_PO4P_CHEM_AUT_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_O2_CHEM_AUT_BAC
    real(kind = DBL_PREC), intent(in) :: DO_STR_HYPOX_CHEM_AUT_BAC_D
    real(kind = DBL_PREC), intent(in) :: THETA_HYPOX_CHEM_AUT_BAC_D
    real(kind = DBL_PREC), intent(in) :: EXPON_HYPOX_CHEM_AUT_BAC_D
    real(kind = DBL_PREC), intent(in) :: CHEM_AUT_BAC_N_TO_C
    real(kind = DBL_PREC), intent(in) :: CHEM_AUT_BAC_P_TO_C
    real(kind = DBL_PREC), intent(in) :: CHEM_AUT_BAC_O2_TO_C
    real(kind = DBL_PREC), intent(in) :: YIELD_CHEM_AUT_BAC
    real(kind = DBL_PREC), intent(in) :: EFF_CHEM_AUT_BAC_GROWTH

    real(kind = DBL_PREC), intent(in) :: TIME_STEP
    integer, intent(in) :: nkn

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NH4_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NO3_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PO4_P
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_OXYGEN
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: CHEM_AUT_BAC_C

    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_TEMP_CHEM_AUT_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_NH4_N_CHEM_AUT_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_PO4_P_CHEM_AUT_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_OXY_CHEM_AUT_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_CHEM_AUT_BAC_GROWTH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_CHEM_AUT_BAC_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_CHEM_AUT_BAC_INT_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_CHEM_AUT_BAC_DEATH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KD_CHEM_AUT_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FAC_HYPOX_CHEM_AUT_BAC_D

    !Chemoautotrophic bacteria growth rate
    LIM_TEMP_CHEM_AUT_BAC  = THETA_KG_CHEM_AUT_BAC ** (TEMP - 2.0D1)
    LIM_NH4_N_CHEM_AUT_BAC = NH4_N / (NH4_N + KHS_NH4N_CHEM_AUT_BAC)
    LIM_PO4_P_CHEM_AUT_BAC = PO4_P / (PO4_P + KHS_PO4P_CHEM_AUT_BAC)
    LIM_OXY_CHEM_AUT_BAC   = DISS_OXYGEN / (KHS_O2_CHEM_AUT_BAC + DISS_OXYGEN)

    R_CHEM_AUT_BAC_GROWTH = KG_CHEM_AUT_BAC_20 * LIM_TEMP_CHEM_AUT_BAC * &
        min(LIM_NH4_N_CHEM_AUT_BAC, LIM_PO4_P_CHEM_AUT_BAC) * &
        LIM_OXY_CHEM_AUT_BAC * CHEM_AUT_BAC_C

    !Chemoautotrophic bacteria respiration rate
    R_CHEM_AUT_BAC_RESP = R_CHEM_AUT_BAC_GROWTH * (1.0D0 - EFF_CHEM_AUT_BAC_GROWTH)

    !Chemoautotrophic bacteria internal respiration rate
    R_CHEM_AUT_BAC_INT_RESP = &
        KR_CHEM_AUT_BAC_20 * (THETA_KR_CHEM_AUT_BAC ** (TEMP - 2.0D1)) * &
        LIM_OXY_CHEM_AUT_BAC * CHEM_AUT_BAC_C

    !Chemoautotrophic bacteria death rate
    KD_CHEM_AUT_BAC = KD_CHEM_AUT_BAC_20 * (THETA_KD_CHEM_AUT_BAC ** (TEMP - 2.0D1))

    where (DISS_OXYGEN <= DO_STR_HYPOX_CHEM_AUT_BAC_D)

        where (DISS_OXYGEN / DO_STR_HYPOX_CHEM_AUT_BAC_D > 1.0D-1)
            FAC_HYPOX_CHEM_AUT_BAC_D = THETA_HYPOX_CHEM_AUT_BAC_D ** &
                (EXPON_HYPOX_CHEM_AUT_BAC_D * &
                    (DO_STR_HYPOX_CHEM_AUT_BAC_D - DISS_OXYGEN))
        elsewhere
            FAC_HYPOX_CHEM_AUT_BAC_D = TIME_STEP / (5.0D-1 * KD_CHEM_AUT_BAC)

            R_CHEM_AUT_BAC_INT_RESP = 0.0D0
            R_CHEM_AUT_BAC_RESP     = 0.0D0
            R_CHEM_AUT_BAC_GROWTH   = 0.0D0
        end where

    elsewhere
        FAC_HYPOX_CHEM_AUT_BAC_D = 1.0D0
    end where

    R_CHEM_AUT_BAC_DEATH = KD_CHEM_AUT_BAC * FAC_HYPOX_CHEM_AUT_BAC_D * CHEM_AUT_BAC_C
end subroutine CHEMOAUTOTROPHS_1



subroutine HETEROTROPHS_1 &
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

    use AQUABC_II_GLOBAL
    implicit none

    real(kind = DBL_PREC), intent(in) :: KG_AER_HET_BAC_20
    real(kind = DBL_PREC), intent(in) :: EFF_AER_HET_BAC_GROWTH
    real(kind = DBL_PREC), intent(in) :: THETA_KG_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KR_AER_HET_BAC_20
    real(kind = DBL_PREC), intent(in) :: THETA_KR_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KD_AER_HET_BAC_20
    real(kind = DBL_PREC), intent(in) :: THETA_KD_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_ORGC_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_ORGN_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_ORGP_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_O2_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_DIN_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_DIP_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_PHYT_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: YIELD_OC_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: OX_ORGN_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_MIN_N
    real(kind = DBL_PREC), intent(in) :: OX_ORGP_AER_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_MIN_P
    real(kind = DBL_PREC), intent(in) :: DO_STR_HYPOX_AER_HET_BAC_D
    real(kind = DBL_PREC), intent(in) :: THETA_HYPOX_AER_HET_BAC_D
    real(kind = DBL_PREC), intent(in) :: EXPON_HYPOX_AER_HET_BAC_D
    real(kind = DBL_PREC), intent(in) :: AER_HET_BAC_N_TO_C
    real(kind = DBL_PREC), intent(in) :: AER_HET_BAC_P_TO_C
    real(kind = DBL_PREC), intent(in) :: AER_HET_BAC_O2_TO_C
    real(kind = DBL_PREC), intent(in) :: KG_FAC_AN_HET_BAC_20
    real(kind = DBL_PREC), intent(in) :: KHS_NO3N_FAC_AN_HET_BAC
    real(kind = DBL_PREC), intent(in) :: YIELD_FAC_AN_HET_BAC

    real(kind = DBL_PREC), intent(in) :: TIME_STEP
    integer, intent(in) :: nkn

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_OXYGEN
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_ORG_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_ORG_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_ORG_P
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NH4_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NO3_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PO4_P
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: AER_HET_BAC_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PHYT_TOT_C


    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_TEMP_AER_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_DISS_ORG_C_AER_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_DISS_ORG_N_AER_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_DISS_ORG_P_AER_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_OXY_AER_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_DIN_AER_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_DIP_AER_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_PHYT_C_AER_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_AER_HET_BAC_GROWTH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_AER_HET_BAC_INT_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KD_AER_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FAC_HYPOX_AER_HET_BAC_D
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_DISS_NO3_N_FAC_AN_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_DENITRIFICATION
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_AER_HET_BAC_DEATH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_AER_HET_BAC_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PREF_NH4N_AER_HET_BAC


    !Aerobic and facultative anaerobic heterotrophic bacteria growth rate in oxic condition
    LIM_TEMP_AER_HET_BAC       = THETA_KG_AER_HET_BAC ** (TEMP - 2.0D1)
    LIM_DISS_ORG_C_AER_HET_BAC = DISS_ORG_C / (KHS_ORGC_AER_HET_BAC + DISS_ORG_C)
    LIM_DISS_ORG_N_AER_HET_BAC = DISS_ORG_N / (KHS_ORGN_AER_HET_BAC + DISS_ORG_N)
    LIM_DISS_ORG_P_AER_HET_BAC = DISS_ORG_P / (KHS_ORGP_AER_HET_BAC + DISS_ORG_P)
    LIM_OXY_AER_HET_BAC        = DISS_OXYGEN /(KHS_O2_AER_HET_BAC   + DISS_OXYGEN)
    LIM_DIN_AER_HET_BAC        = (NH4_N + NO3_N)/(KHS_DIN_AER_HET_BAC + NH4_N + NO3_N)
    LIM_DIP_AER_HET_BAC        = PO4_P/(KHS_DIP_AER_HET_BAC + PO4_P)
    LIM_PHYT_C_AER_HET_BAC     = PHYT_TOT_C/(KHS_PHYT_AER_HET_BAC + PHYT_TOT_C)

    call AMMONIA_PREFS(PREF_NH4N_AER_HET_BAC,NH4_N, NO3_N, KHS_DIN_AER_HET_BAC,nkn)

    R_AER_HET_BAC_GROWTH = KG_AER_HET_BAC_20 * LIM_TEMP_AER_HET_BAC * &
        min(LIM_DISS_ORG_C_AER_HET_BAC, LIM_DISS_ORG_N_AER_HET_BAC, LIM_DISS_ORG_P_AER_HET_BAC, &
            LIM_DIN_AER_HET_BAC, LIM_DIP_AER_HET_BAC) * &
            LIM_OXY_AER_HET_BAC * LIM_PHYT_C_AER_HET_BAC * AER_HET_BAC_C

    !Aerobic- anaerobic heterotrophic bacteria internal respiration rate in oxic condition
    R_AER_HET_BAC_INT_RESP = R_AER_HET_BAC_GROWTH * (1.0D0 - EFF_AER_HET_BAC_GROWTH) + & ! active respiration
        KR_AER_HET_BAC_20 * (THETA_KR_AER_HET_BAC ** (TEMP - 2.0D1)) * &                 ! basal respiration
        LIM_OXY_AER_HET_BAC * AER_HET_BAC_C

    !Aerobic-anaerobic heterotrophic bacteria death rate in oxic conditions
    KD_AER_HET_BAC = KD_AER_HET_BAC_20 * (THETA_KD_AER_HET_BAC ** (TEMP - 2.0D1)) !normal condition

    ! approaching to anoxic conditions:   (not tested!)
    where (DISS_OXYGEN / DO_STR_HYPOX_AER_HET_BAC_D > 2.0D0)          ! Critical value when exponential death starts should be a parameter fixme
            FAC_HYPOX_AER_HET_BAC_D = THETA_HYPOX_AER_HET_BAC_D ** &  ! Returning from low oxygen to higher should be different fixme
                (EXPON_HYPOX_AER_HET_BAC_D * &
                    (DISS_OXYGEN -  DO_STR_HYPOX_AER_HET_BAC_D))
    end where

    ! switching to nitrate respiration
    where (DISS_OXYGEN <= DO_STR_HYPOX_AER_HET_BAC_D)

            R_AER_HET_BAC_INT_RESP  = 0.0D0

            LIM_DISS_NO3_N_FAC_AN_HET_BAC = NO3_N / (NO3_N + KHS_NO3N_FAC_AN_HET_BAC)

            R_AER_HET_BAC_GROWTH    = KG_FAC_AN_HET_BAC_20 *     &
               LIM_TEMP_AER_HET_BAC * &
               LIM_DISS_NO3_N_FAC_AN_HET_BAC * &
               min(LIM_DISS_ORG_C_AER_HET_BAC, LIM_DISS_ORG_N_AER_HET_BAC, &
               LIM_DISS_ORG_P_AER_HET_BAC, LIM_DIN_AER_HET_BAC, LIM_DIP_AER_HET_BAC) &
                * LIM_PHYT_C_AER_HET_BAC * AER_HET_BAC_C
            R_DENITRIFICATION = (1.0D0 / YIELD_FAC_AN_HET_BAC) * R_AER_HET_BAC_GROWTH

    elsewhere
        FAC_HYPOX_AER_HET_BAC_D = 1.0D0
        R_DENITRIFICATION = 0.0D0
    end where

    R_AER_HET_BAC_DEATH = KD_AER_HET_BAC * FAC_HYPOX_AER_HET_BAC_D * AER_HET_BAC_C

    !Aerobic heterotrophic bacteria external respiration(oxidation) rate
    !R_AER_HET_BAC_RESP = R_AER_HET_BAC_GROWTH * (1.0D0 - EFF_AER_HET_BAC_GROWTH)
    R_AER_HET_BAC_RESP = R_AER_HET_BAC_GROWTH * (1-YIELD_OC_AER_HET_BAC)/YIELD_OC_AER_HET_BAC !no loss of bacteria carbon
end subroutine HETEROTROPHS_1



subroutine HETEROTROPHS_2 &
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

    use AQUABC_II_GLOBAL
    implicit none


    real(kind = DBL_PREC), intent(in) :: KG_FAC_AN_HET_BAC_20
    real(kind = DBL_PREC), intent(in) :: THETA_KG_FAC_AN_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KR_FAC_AN_HET_BAC_20
    real(kind = DBL_PREC), intent(in) :: THETA_KR_FAC_AN_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KD_FAC_AN_HET_BAC_20
    real(kind = DBL_PREC), intent(in) :: THETA_KD_FAC_AN_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_NO3N_FAC_AN_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_ORGC_FAC_AN_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_ORGN_FAC_AN_HET_BAC
    real(kind = DBL_PREC), intent(in) :: KHS_ORGP_FAC_AN_HET_BAC
    real(kind = DBL_PREC), intent(in) :: REV_KHS_O2_FAC_AN_HET_BAC
    real(kind = DBL_PREC), intent(in) :: NO3N_LACK_STR_FAC_AN_HET_BAC_D
    real(kind = DBL_PREC), intent(in) :: THETA_NO3_LACK_FAC_AN_HET_BAC_D
    real(kind = DBL_PREC), intent(in) :: EXP_NO3_LACK_FAC_AN_HET_BAC_D
    real(kind = DBL_PREC), intent(in) :: FAC_AN_HET_BAC_N_TO_C
    real(kind = DBL_PREC), intent(in) :: FAC_AN_HET_BAC_P_TO_C
    real(kind = DBL_PREC), intent(in) :: FAC_AN_HET_BAC_O2_TO_C
    real(kind = DBL_PREC), intent(in) :: YIELD_FAC_AN_HET_BAC
    real(kind = DBL_PREC), intent(in) :: EFF_FAC_AN_HET_BAC_GROWTH

    real(kind = DBL_PREC), intent(in) :: TIME_STEP
    integer, intent(in) :: nkn

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NO3_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_ORG_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_ORG_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_ORG_P
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_OXYGEN
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FAC_AN_HET_BAC_C

    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_TEMP_FAC_AN_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_DISS_NO3_N_FAC_AN_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_DISS_ORG_C_FAC_AN_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_DISS_ORG_N_FAC_AN_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_DISS_ORG_P_FAC_AN_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_OXY_FAC_AN_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_FAC_AN_HET_BAC_GROWTH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_FAC_AN_HET_BAC_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_FAC_AN_HET_BAC_INT_RESP
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KD_FAC_AN_HET_BAC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_FAC_AN_HET_BAC_DEATH

    !Facultative anaerobic heterotrohic bacteria growth rate. Nitrate respiration threshold 0.5mg/l O2
    ! Switched of for the Curonian lagoon making zero boundaries and initial condition
    ! Formulation should be corrected for anaerobic conditions
    LIM_TEMP_FAC_AN_HET_BAC       =  THETA_KG_FAC_AN_HET_BAC ** (TEMP - 2.0D1)
    LIM_DISS_NO3_N_FAC_AN_HET_BAC = NO3_N / (NO3_N + KHS_NO3N_FAC_AN_HET_BAC)
    LIM_DISS_ORG_C_FAC_AN_HET_BAC = DISS_ORG_C / (KHS_ORGC_FAC_AN_HET_BAC + DISS_ORG_C)
    LIM_DISS_ORG_N_FAC_AN_HET_BAC = DISS_ORG_N / (KHS_ORGN_FAC_AN_HET_BAC + DISS_ORG_N)
    LIM_DISS_ORG_P_FAC_AN_HET_BAC = DISS_ORG_P / (KHS_ORGP_FAC_AN_HET_BAC + DISS_ORG_P)

    LIM_OXY_FAC_AN_HET_BAC        = &
        REV_KHS_O2_FAC_AN_HET_BAC / (REV_KHS_O2_FAC_AN_HET_BAC + DISS_OXYGEN)

    R_FAC_AN_HET_BAC_GROWTH = &
        KG_FAC_AN_HET_BAC_20 * LIM_TEMP_FAC_AN_HET_BAC * &
        LIM_DISS_NO3_N_FAC_AN_HET_BAC * &
        min(LIM_DISS_ORG_C_FAC_AN_HET_BAC, LIM_DISS_ORG_N_FAC_AN_HET_BAC, &
            LIM_DISS_ORG_P_FAC_AN_HET_BAC) * LIM_OXY_FAC_AN_HET_BAC * FAC_AN_HET_BAC_C

    !Facultative anaerobic heterotrohic bacteria respiration rate
    R_FAC_AN_HET_BAC_RESP = R_FAC_AN_HET_BAC_GROWTH * (1.0D0 - EFF_FAC_AN_HET_BAC_GROWTH)

    !Facultative anaerobic heterotrohic bacteria internal respiration rate
    R_FAC_AN_HET_BAC_INT_RESP = &
        KR_FAC_AN_HET_BAC_20 * (THETA_KR_FAC_AN_HET_BAC ** (TEMP - 2.0D1)) * &
        LIM_DISS_NO3_N_FAC_AN_HET_BAC * LIM_OXY_FAC_AN_HET_BAC * FAC_AN_HET_BAC_C

    !Facultative anaerobic heterotrohic bacteria death rate
    KD_FAC_AN_HET_BAC = KD_FAC_AN_HET_BAC_20 * (THETA_KD_FAC_AN_HET_BAC ** (TEMP - 2.0D1))

    ! The following ic commented by Petras because it caused compilation errors
!    where (NO3_N <= NO3N_LACK_STR_FAC_AN_HET_BAC_D)


!         where (NO3_N / NO3N_LACK_STR_FAC_AN_HET_BAC_D > 1.0D-1)
!             FAC_NO3N_LACK_FAC_AN_HET_BAC_D = &
!                 THETA_NO3_LACK_FAC_AN_HET_BAC_D ** &
!                     (EXP_NO3_LACK_FAC_AN_HET_BAC_D * &
!                         (NO3N_LACK_STR_FAC_AN_HET_BAC_D - NO3_N))
!         elsewhere
!             NO3N_LACK_STR_FAC_AN_HET_BAC_D = TIME_STEP / (5.0D-1 * KD_FAC_AN_HET_BAC)
!             R_FAC_AN_HET_BAC_INT_RESP      = 0.0D0
!             R_FAC_AN_HET_BAC_RESP          = 0.0D0
!             R_FAC_AN_HET_BAC_GROWTH        = 0.0D0
!         end where
!
!     elsewhere
!
!         NO3N_LACK_STR_FAC_AN_HET_BAC_D = 1.0D0
!     end where

    R_FAC_AN_HET_BAC_DEATH = KD_FAC_AN_HET_BAC * &
                 NO3N_LACK_STR_FAC_AN_HET_BAC_D * FAC_AN_HET_BAC_C
end subroutine HETEROTROPHS_2



subroutine ORGANIC_CARBON_DISSOLUTION &
           (FAC_PHYT_DET_PART_ORG_C     , &
            KDISS_DET_PART_ORG_C_20     , &
            THETA_KDISS_DET_PART_ORG_C  , &
            nkn                         , &
            TEMP                        , &
            DET_PART_ORG_C              , &
            PHYT_TOT_C                  , &
            LIM_PHYT_DISS_DET_PART_ORG_C, &
            R_DET_PART_ORG_C_DISSOLUTION)

    use AQUABC_II_GLOBAL
    implicit none

    real(kind = DBL_PREC), intent(in) :: FAC_PHYT_DET_PART_ORG_C
    real(kind = DBL_PREC), intent(in) :: KDISS_DET_PART_ORG_C_20
    real(kind = DBL_PREC), intent(in) :: THETA_KDISS_DET_PART_ORG_C

    integer, intent(in) :: nkn

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DET_PART_ORG_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PHYT_TOT_C

    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_PHYT_DISS_DET_PART_ORG_C
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_DET_PART_ORG_C_DISSOLUTION

    !Algal dependent hydrolysis rate
    LIM_PHYT_DISS_DET_PART_ORG_C = FAC_PHYT_DET_PART_ORG_C * PHYT_TOT_C

    R_DET_PART_ORG_C_DISSOLUTION = &
           (KDISS_DET_PART_ORG_C_20 + LIM_PHYT_DISS_DET_PART_ORG_C) * &
           (THETA_KDISS_DET_PART_ORG_C ** (TEMP - 2.0D1)) * DET_PART_ORG_C

end subroutine ORGANIC_CARBON_DISSOLUTION




subroutine ORGANIC_CARBON_MINERALIZATION &
           (FAC_PHYT_AMIN_DOC           , &
            K_MIN_DOC_20                , &
            THETA_K_MIN_DOC             , &
            nkn                         , &
            TEMP                        , &
            DISS_ORG_C                  , &
            PHYT_TOT_C                  , &
            LIM_PHYT_AMIN_DOC           , &
            R_ABIOTIC_DOC_MIN)

    use AQUABC_II_GLOBAL
    implicit none

    real(kind = DBL_PREC), intent(in) :: FAC_PHYT_AMIN_DOC
    real(kind = DBL_PREC), intent(in) :: K_MIN_DOC_20
    real(kind = DBL_PREC), intent(in) :: THETA_K_MIN_DOC

    integer, intent(in) :: nkn

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_ORG_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PHYT_TOT_C

    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_PHYT_AMIN_DOC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ABIOTIC_DOC_MIN

    LIM_PHYT_AMIN_DOC = FAC_PHYT_AMIN_DOC * PHYT_TOT_C

    R_ABIOTIC_DOC_MIN = (K_MIN_DOC_20 + LIM_PHYT_AMIN_DOC) * &
                        (THETA_K_MIN_DOC ** (TEMP - 2.0D1)) * DISS_ORG_C
end subroutine ORGANIC_CARBON_MINERALIZATION
