
!--------------------------------------------------------------------------
!
!    Copyright (C) 2005-2018  Ali Erturk  
!    Copyright (C) 2005-2018  Petras Zemlys     
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
!subroutine REDOX_AND_SPECIATION
!subroutine DO_SATURATION_VEC
!subroutine IRON_II_DISSOLUTION
!subroutine IRON_II_OXIDATION
!subroutine IP_SOLUBLE_FRACTION
!subroutine CALC_DISS_ME_CONC


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

    R_DIA_MET  = R_DIA_GROWTH * (1.0D0 - EFF_DIA_GROWTH)
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
            DON                     , &
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
            frac_avail_DON          , &
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
            PREF_NH4N_DON_CYN)

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
    real(kind = DBL_PREC), intent(in) :: frac_avail_DON
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NH4_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NO3_N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DON

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
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PREF_NH4N_DON_CYN
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

    LIM_KG_CYN_N    = (NH4_N + (DON * frac_avail_DON) + NO3_N) / &
                      (KHS_DIN_CYN + NH4_N + (DON * frac_avail_DON) + NO3_N)

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
    call AMMONIA_DON_PREFS(PREF_NH4N_DON_CYN, NH4_N, DON, frac_avail_DON, NO3_N, KHS_DIN_CYN,nkn)
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
            frac_avail_DON               , &
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
            DON                          , &
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
            PREF_NH4N_DON_FIX_CYN)

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
    real(kind = DBL_PREC), intent(in) :: frac_avail_DON
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DON
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
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PREF_NH4N_DON_FIX_CYN

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
    	                     ! Changed to 1 by Petras 2014 10 13
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
    LIM_KG_NON_FIX_CYN_N    = (NH4_N + (DON * frac_avail_DON) + NO3_N) / &
                              (KHS_DIN_FIX_CYN + NH4_N +(DON * frac_avail_DON) + NO3_N)

    LIM_KG_NON_FIX_CYN_P    = PO4_P / (KHS_DIP_FIX_CYN + PO4_P)
    LIM_KG_NON_FIX_CYN_NUTR = min(LIM_KG_NON_FIX_CYN_N, LIM_KG_NON_FIX_CYN_P)

    !Nutrient limitation of fixing cyanobacteria in fixing fraction
    LIM_KG_FIX_FIX_CYN_N    = (K_FIX / (K_FIX + NH4_N +(DON * frac_avail_DON) + NO3_N))
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
    call AMMONIA_DON_PREFS(PREF_NH4N_DON_FIX_CYN, NH4_N, DON, frac_avail_DON, NO3_N, KHS_DIN_FIX_CYN,nkn)

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
            GRAT_ZOO_DET_PART_ORG_C       , &
            PREF_ZOO_DIA                  , &
            PREF_ZOO_CYN                  , &
            PREF_ZOO_FIX_CYN              , &
            PREF_ZOO_OPA                  , &
            PREF_ZOO_DET_PART_ORG_C       , &
            KHS_DIA_C_ZOO                 , &
            KHS_CYN_C_ZOO                 , &
            KHS_FIX_CYN_C_ZOO             , &
            KHS_OPA_C_ZOO                 , &
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
            DET_PART_ORG_C                , &
            ZOO_C                         , &
            TIME_STEP                     , &
            nkn                           , &
            KG_ZOO                        , &
            KG_ZOO_DIA                    , &
            KG_ZOO_CYN                    , &
            KG_ZOO_OPA                    , &
            KG_ZOO_FIX_CYN                , &
            KG_ZOO_DET_PART_ORG_C         , &
            KD_ZOO                        , &
            FOOD_FACTOR_ZOO_DIA           , &
            FOOD_FACTOR_ZOO_CYN           , &
            FOOD_FACTOR_ZOO_OPA           , &
            FOOD_FACTOR_ZOO_FIX_CYN       , &
            FOOD_FACTOR_ZOO_DET_PART_ORG_C, &
            R_ZOO_FEEDING_DIA             , &
            R_ZOO_FEEDING_CYN             , &
            R_ZOO_FEEDING_FIX_CYN         , &
            R_ZOO_FEEDING_OPA             , &
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
    real(kind = DBL_PREC), intent(in) :: GRAT_ZOO_DET_PART_ORG_C
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_DIA
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_CYN
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_FIX_CYN
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_OPA
    real(kind = DBL_PREC), intent(in) :: PREF_ZOO_DET_PART_ORG_C
    real(kind = DBL_PREC), intent(in) :: KHS_DIA_C_ZOO
    real(kind = DBL_PREC), intent(in) :: KHS_CYN_C_ZOO
    real(kind = DBL_PREC), intent(in) :: KHS_FIX_CYN_C_ZOO
    real(kind = DBL_PREC), intent(in) :: KHS_OPA_C_ZOO
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
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KG_ZOO_DET_PART_ORG_C
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: KD_ZOO
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_DIA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_OPA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_FIX_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FOOD_FACTOR_ZOO_DET_PART_ORG_C
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_FEEDING_DIA
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_FEEDING_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_FEEDING_FIX_CYN
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ZOO_FEEDING_OPA
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
    R_ZOO_FEEDING_DET_PART_ORG_C = KG_ZOO_DET_PART_ORG_C * FOOD_FACTOR_ZOO_DET_PART_ORG_C * ZOO_C

    !Zooplankton excretion rate
    ACTUAL_ZOO_N_TO_C = ZOO_N_TO_C
    ACTUAL_ZOO_P_TO_C = ZOO_P_TO_C

    !Zooplankton Growth rate
    R_ZOO_GROWTH =  R_ZOO_FEEDING_DIA         + R_ZOO_FEEDING_CYN            + &
                    R_ZOO_FEEDING_OPA         + R_ZOO_FEEDING_DET_PART_ORG_C

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
            K_MIN_DOC_DOXY_20           , &
            K_MIN_DOC_NO3N_20           , &
            K_MIN_DOC_MN_IV_20          , &
            K_MIN_DOC_FE_III_20         , &
            K_MIN_DOC_S_PLUS_6_20       , &
            K_MIN_DOC_DOC_20            , &
            THETA_K_MIN_DOC_DOXY        , &
            THETA_K_MIN_DOC_NO3N        , &
            THETA_K_MIN_DOC_MN_IV       , &
            THETA_K_MIN_DOC_FE_III      , &
            THETA_K_MIN_DOC_S_PLUS_6    , &
            THETA_K_MIN_DOC_DOC         , &
            K_HS_DOC_MIN_DOXY           , &
            K_HS_DOC_MIN_NO3N           , &
            K_HS_DOC_MIN_MN_IV          , &
            K_HS_DOC_MIN_FE_III         , &
            K_HS_DOC_MIN_S_PLUS_6       , &
            K_HS_DOC_MIN_DOC            , &
            K_HS_DOXY_RED_LIM           , &
            K_HS_NO3N_RED_LIM           , &
            K_HS_MN_IV_RED_LIM          , &
            K_HS_FE_III_RED_LIM         , &
            K_HS_S_PLUS_6_RED_LIM       , &
            K_HS_DOXY_RED_INHB          , &
            K_HS_NO3N_RED_INHB          , &
            K_HS_MN_IV_RED_INHB         , &
            K_HS_FE_III_RED_INHB        , &
            K_HS_S_PLUS_6_RED_INHB      , &
            PH_MIN_DOC_MIN_DOXY         , &  !Min. pH for the optimum pH range for DOC mineralization with DOXY     as final electron acceptor (subroutine input)
            PH_MIN_DOC_MIN_NO3N         , &  !Min. pH for the optimum pH range for DOC mineralization with NO3N     as final electron acceptor (subroutine input)
            PH_MIN_DOC_MIN_MN_IV        , &  !Min. pH for the optimum pH range for DOC mineralization with MN_IV    as final electron acceptor (subroutine input)
            PH_MIN_DOC_MIN_FE_III       , &  !Min. pH for the optimum pH range for DOC mineralization with FE_III   as final electron acceptor (subroutine input)
            PH_MIN_DOC_MIN_S_PLUS_6     , &  !Min. pH for the optimum pH range for DOC mineralization with S_PLUS_6 as final electron acceptor (subroutine input)
            PH_MIN_DOC_MIN_DOC          , &  !Min. pH for the optimum pH range for DOC mineralization with DOC      as final electron acceptor (subroutine input)
            PH_MAX_DOC_MIN_DOXY         , &  !Max. pH for the optimum pH range for DOC mineralization with DOXY     as final electron acceptor (subroutine input)
            PH_MAX_DOC_MIN_NO3N         , &  !Max. pH for the optimum pH range for DOC mineralization with NO3N     as final electron acceptor (subroutine input)
            PH_MAX_DOC_MIN_MN_IV        , &  !Max. pH for the optimum pH range for DOC mineralization with MN_IV    as final electron acceptor (subroutine input)
            PH_MAX_DOC_MIN_FE_III       , &  !Max. pH for the optimum pH range for DOC mineralization with FE_III   as final electron acceptor (subroutine input)
            PH_MAX_DOC_MIN_S_PLUS_6     , &  !Max. pH for the optimum pH range for DOC mineralization with S_PLUS_6 as final electron acceptor (subroutine input)
            PH_MAX_DOC_MIN_DOC          , &  !Max. pH for the optimum pH range for DOC mineralization with DOC      as final electron acceptor (subroutine input)
            nkn                         , &
            TEMP                        , &
            DISS_ORG_C                  , &
            PHYT_TOT_C                  , &
            DOXY                        , &
            NO3N                        , &
            MN_IV                       , &
            FE_III                      , &
            S_PLUS_6                    , &
            PH                          , &
            LIM_DOXY_RED                , &
            LIM_NO3N_RED                , &
            LIM_MN_IV_RED               , &
            LIM_FE_III_RED              , &
            LIM_S_PLUS_6_RED            , &
            LIM_DOC_RED                 , &
            LIM_PHYT_AMIN_DOC           , &
            PH_CORR_DOC_MIN_DOXY        , &  !pH correction for DOC mineralization with DOXY     as final electron acceptor (subroutine output)
            PH_CORR_DOC_MIN_NO3N        , &  !pH correction for DOC mineralization with NO3N     as final electron acceptor (subroutine output)
            PH_CORR_DOC_MIN_MN_IV       , &  !pH correction for DOC mineralization with MN_IV    as final electron acceptor (subroutine output)
            PH_CORR_DOC_MIN_FE_III      , &  !pH correction for DOC mineralization with FE_III   as final electron acceptor (subroutine output)
            PH_CORR_DOC_MIN_S_PLUS_6    , &  !pH correction for DOC mineralization with S_PLUS_6 as final electron acceptor (subroutine output)
            PH_CORR_DOC_MIN_DOC         , &  !pH correction for DOC mineralization with DOC      as final electron acceptor (subroutine output)
            K_NO3_RED                   , &
            K_MN_IV_RED                 , &
            K_FE_III_RED                , &
            K_S_PLUS_6_RED              , &
            K_DOC_RED                   , &
            R_ABIOTIC_DOC_MIN_DOXY      , &  !Process rate  for DOC mineralization with DOXY     as final electron acceptor (subroutine output)
            R_ABIOTIC_DOC_MIN_NO3N      , &  !Process rate  for DOC mineralization with NO3N     as final electron acceptor (subroutine output)
            R_ABIOTIC_DOC_MIN_MN_IV     , &  !Process rate  for DOC mineralization with MN_IV    as final electron acceptor (subroutine output)
            R_ABIOTIC_DOC_MIN_FE_III    , &  !Process rate  for DOC mineralization with FE_III   as final electron acceptor (subroutine output)
            R_ABIOTIC_DOC_MIN_S_PLUS_6  , &  !Process rate  for DOC mineralization with S_PLUS_6 as final electron acceptor (subroutine output)
            R_ABIOTIC_DOC_MIN_DOC)           !Process rate  for DOC mineralization with DOC      as final electron acceptor (subroutine output)

    ! ----------------------------------------------------------------------------------------
    ! Subroutine for organic carbon mineraliztion
    ! This subroutine is almost completely rewritten to be compitable with the redox sequences
    ! ----------------------------------------------------------------------------------------
    use AQUABC_II_GLOBAL
    implicit none

    real(kind = DBL_PREC), intent(in) :: FAC_PHYT_AMIN_DOC

    real(kind = DBL_PREC), intent(in) :: K_MIN_DOC_DOXY_20
    real(kind = DBL_PREC), intent(in) :: K_MIN_DOC_NO3N_20
    real(kind = DBL_PREC), intent(in) :: K_MIN_DOC_MN_IV_20
    real(kind = DBL_PREC), intent(in) :: K_MIN_DOC_FE_III_20
    real(kind = DBL_PREC), intent(in) :: K_MIN_DOC_S_PLUS_6_20
    real(kind = DBL_PREC), intent(in) :: K_MIN_DOC_DOC_20
    real(kind = DBL_PREC), intent(in) :: THETA_K_MIN_DOC_DOXY
    real(kind = DBL_PREC), intent(in) :: THETA_K_MIN_DOC_NO3N
    real(kind = DBL_PREC), intent(in) :: THETA_K_MIN_DOC_MN_IV
    real(kind = DBL_PREC), intent(in) :: THETA_K_MIN_DOC_FE_III
    real(kind = DBL_PREC), intent(in) :: THETA_K_MIN_DOC_S_PLUS_6
    real(kind = DBL_PREC), intent(in) :: THETA_K_MIN_DOC_DOC
    real(kind = DBL_PREC), intent(in) :: K_HS_DOC_MIN_DOXY
    real(kind = DBL_PREC), intent(in) :: K_HS_DOC_MIN_NO3N
    real(kind = DBL_PREC), intent(in) :: K_HS_DOC_MIN_MN_IV
    real(kind = DBL_PREC), intent(in) :: K_HS_DOC_MIN_FE_III
    real(kind = DBL_PREC), intent(in) :: K_HS_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC), intent(in) :: K_HS_DOC_MIN_DOC
    real(kind = DBL_PREC), intent(in) :: K_HS_DOXY_RED_LIM
    real(kind = DBL_PREC), intent(in) :: K_HS_NO3N_RED_LIM
    real(kind = DBL_PREC), intent(in) :: K_HS_MN_IV_RED_LIM
    real(kind = DBL_PREC), intent(in) :: K_HS_FE_III_RED_LIM
    real(kind = DBL_PREC), intent(in) :: K_HS_S_PLUS_6_RED_LIM
    real(kind = DBL_PREC), intent(in) :: K_HS_DOXY_RED_INHB
    real(kind = DBL_PREC), intent(in) :: K_HS_NO3N_RED_INHB
    real(kind = DBL_PREC), intent(in) :: K_HS_MN_IV_RED_INHB
    real(kind = DBL_PREC), intent(in) :: K_HS_FE_III_RED_INHB
    real(kind = DBL_PREC), intent(in) :: K_HS_S_PLUS_6_RED_INHB
    real(kind = DBL_PREC), intent(in) :: PH_MIN_DOC_MIN_DOXY
    real(kind = DBL_PREC), intent(in) :: PH_MIN_DOC_MIN_NO3N
    real(kind = DBL_PREC), intent(in) :: PH_MIN_DOC_MIN_MN_IV
    real(kind = DBL_PREC), intent(in) :: PH_MIN_DOC_MIN_FE_III
    real(kind = DBL_PREC), intent(in) :: PH_MIN_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC), intent(in) :: PH_MIN_DOC_MIN_DOC
    real(kind = DBL_PREC), intent(in) :: PH_MAX_DOC_MIN_DOXY
    real(kind = DBL_PREC), intent(in) :: PH_MAX_DOC_MIN_NO3N
    real(kind = DBL_PREC), intent(in) :: PH_MAX_DOC_MIN_MN_IV
    real(kind = DBL_PREC), intent(in) :: PH_MAX_DOC_MIN_FE_III
    real(kind = DBL_PREC), intent(in) :: PH_MAX_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC), intent(in) :: PH_MAX_DOC_MIN_DOC

    integer, intent(in) :: nkn

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_ORG_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PHYT_TOT_C
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DOXY
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NO3N
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: MN_IV
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FE_III
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PH
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: LIM_DOXY_RED
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: LIM_NO3N_RED
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: LIM_MN_IV_RED
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: LIM_FE_III_RED
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: LIM_S_PLUS_6_RED
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: LIM_DOC_RED

    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_PHYT_AMIN_DOC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PH_CORR_DOC_MIN_DOXY
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PH_CORR_DOC_MIN_NO3N
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PH_CORR_DOC_MIN_MN_IV
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PH_CORR_DOC_MIN_FE_III
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PH_CORR_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PH_CORR_DOC_MIN_DOC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ABIOTIC_DOC_MIN_DOXY
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ABIOTIC_DOC_MIN_NO3N
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ABIOTIC_DOC_MIN_MN_IV
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ABIOTIC_DOC_MIN_FE_III
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ABIOTIC_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_ABIOTIC_DOC_MIN_DOC
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: K_NO3_RED
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: K_MN_IV_RED
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: K_FE_III_RED
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: K_S_PLUS_6_RED
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: K_DOC_RED
    
   

    LIM_PHYT_AMIN_DOC = FAC_PHYT_AMIN_DOC * PHYT_TOT_C

    call CALCULATE_PH_CORR(PH_CORR_DOC_MIN_DOXY    , PH, PH_MIN_DOC_MIN_DOXY    , PH_MAX_DOC_MIN_DOXY    , nkn)
    call CALCULATE_PH_CORR(PH_CORR_DOC_MIN_NO3N    , PH, PH_MIN_DOC_MIN_NO3N    , PH_MAX_DOC_MIN_NO3N    , nkn)
    call CALCULATE_PH_CORR(PH_CORR_DOC_MIN_MN_IV   , PH, PH_MIN_DOC_MIN_MN_IV   , PH_MAX_DOC_MIN_MN_IV   , nkn)
    call CALCULATE_PH_CORR(PH_CORR_DOC_MIN_FE_III  , PH, PH_MIN_DOC_MIN_FE_III  , PH_MAX_DOC_MIN_FE_III  , nkn)
    call CALCULATE_PH_CORR(PH_CORR_DOC_MIN_S_PLUS_6, PH, PH_MIN_DOC_MIN_S_PLUS_6, PH_MAX_DOC_MIN_S_PLUS_6, nkn)
    call CALCULATE_PH_CORR(PH_CORR_DOC_MIN_DOC     , PH, PH_MIN_DOC_MIN_DOC     , PH_MAX_DOC_MIN_DOC     , nkn)

    R_ABIOTIC_DOC_MIN_DOXY = &
        (K_MIN_DOC_DOXY_20 + LIM_PHYT_AMIN_DOC) * (THETA_K_MIN_DOC_DOXY ** (TEMP - 2.0D1)) * &
        LIM_DOXY_RED * PH_CORR_DOC_MIN_DOXY * (DISS_ORG_C / (DISS_ORG_C + K_HS_DOC_MIN_DOXY)) * &
        DISS_ORG_C

    R_ABIOTIC_DOC_MIN_NO3N = &
        K_MIN_DOC_NO3N_20  * (THETA_K_MIN_DOC_NO3N ** (TEMP - 2.0D1)) * &
        LIM_NO3N_RED * PH_CORR_DOC_MIN_NO3N * (DISS_ORG_C / (DISS_ORG_C + K_HS_DOC_MIN_NO3N)) * &
        DISS_ORG_C

    R_ABIOTIC_DOC_MIN_MN_IV = &
        K_MIN_DOC_MN_IV_20  * (THETA_K_MIN_DOC_MN_IV ** (TEMP - 2.0D1)) * &
        LIM_MN_IV_RED * PH_CORR_DOC_MIN_MN_IV * (DISS_ORG_C / (DISS_ORG_C + K_HS_DOC_MIN_MN_IV)) * &
        DISS_ORG_C

    R_ABIOTIC_DOC_MIN_FE_III = &
        K_MIN_DOC_FE_III_20  * (THETA_K_MIN_DOC_FE_III ** (TEMP - 2.0D1)) * &
        LIM_FE_III_RED * PH_CORR_DOC_MIN_FE_III * (DISS_ORG_C / (DISS_ORG_C + K_HS_DOC_MIN_FE_III)) * &
        DISS_ORG_C

    R_ABIOTIC_DOC_MIN_S_PLUS_6 = &
        K_MIN_DOC_S_PLUS_6_20  * (THETA_K_MIN_DOC_S_PLUS_6 ** (TEMP - 2.0D1)) * &
        LIM_S_PLUS_6_RED * PH_CORR_DOC_MIN_S_PLUS_6 * (DISS_ORG_C / (DISS_ORG_C + K_HS_DOC_MIN_S_PLUS_6)) * &
        DISS_ORG_C

    R_ABIOTIC_DOC_MIN_DOC = &
        (K_MIN_DOC_DOC_20  * (THETA_K_MIN_DOC_DOC ** (TEMP - 2.0D1)) * &
         LIM_DOC_RED * PH_CORR_DOC_MIN_DOXY * (DISS_ORG_C / (DISS_ORG_C + K_HS_DOC_MIN_DOC)) * DISS_ORG_C)

end subroutine ORGANIC_CARBON_MINERALIZATION

!***********************************************************************
!***********************************************************************       
    
        
subroutine REDOX_AND_SPECIATION &
    (DOXY, NO3N, MN_IV, FE_III, S_PLUS_6, DISS_ORG_C, &
     S_MINUS_2 , MN_II, FE_II , HCO3    , CO3, &
     TEMP, SALT, PH, ELEVATION, &
     K_HS_DOXY_RED_LIM   , K_HS_NO3N_RED_LIM , K_HS_MN_IV_RED_LIM , &
     K_HS_FE_III_RED_LIM , K_HS_S_PLUS_6_RED_LIM, &
     K_HS_DOXY_RED_INHB  , K_HS_NO3N_RED_INHB, K_HS_MN_IV_RED_INHB, &
     K_HS_FE_III_RED_INHB, K_HS_S_PLUS_6_RED_INHB, nkn, &
     LIM_DOXY_RED        , LIM_NO3N_RED          , LIM_MN_IV_RED  , &
     LIM_FE_III_RED      , LIM_S_PLUS_6_RED      , LIM_DOC_RED, &
     PE, FE_II_DISS, FE_III_DISS, MN_II_DISS)

    use AQUABC_II_GLOBAL
    use mod_debug

    implicit none

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DOXY          ! Dissolved oxygen (mg/L)
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: NO3N          ! Nitrate nitrogen (mg/L)
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: MN_IV         ! Mn IV            (mg/L)
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FE_III        ! Fe III           (mg/L)
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: S_PLUS_6      ! S +VI            (mg/L)
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DISS_ORG_C    ! DOC              (mg/L)
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: S_MINUS_2     ! S -II            (mg/L)
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: MN_II         ! Mn II            (mg/L)
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FE_II         ! Fe II            (mg/L)
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: HCO3          ! Bicarbonates     (moles)
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: CO3           ! Carbonate        (moles)

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: SALT
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PH
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: ELEVATION

    real(kind = DBL_PREC), intent(in) :: K_HS_DOXY_RED_LIM
    real(kind = DBL_PREC), intent(in) :: K_HS_NO3N_RED_LIM
    real(kind = DBL_PREC), intent(in) :: K_HS_MN_IV_RED_LIM
    real(kind = DBL_PREC), intent(in) :: K_HS_FE_III_RED_LIM
    real(kind = DBL_PREC), intent(in) :: K_HS_S_PLUS_6_RED_LIM
    real(kind = DBL_PREC), intent(in) :: K_HS_DOXY_RED_INHB
    real(kind = DBL_PREC), intent(in) :: K_HS_NO3N_RED_INHB
    real(kind = DBL_PREC), intent(in) :: K_HS_MN_IV_RED_INHB
    real(kind = DBL_PREC), intent(in) :: K_HS_FE_III_RED_INHB
    real(kind = DBL_PREC), intent(in) :: K_HS_S_PLUS_6_RED_INHB

    integer, intent(in) :: nkn

    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_DOXY_RED
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_NO3N_RED
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_MN_IV_RED
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_FE_III_RED
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_S_PLUS_6_RED
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: LIM_DOC_RED
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PE
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FE_II_DISS
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FE_III_DISS
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: MN_II_DISS
    
    integer, dimension(nkn) :: REDUCED_AGENT_NO

    real(kind = DBL_PREC), dimension(nkn, 6) :: REDUCER_LIM_FACTORS
    real(kind = DBL_PREC), dimension(nkn) :: CS, H_PLUS, SO4_MOLAR, HS_MOLAR, S_MINUS_2_MOLAR

    real(kind = DBL_PREC), dimension(nkn) :: FE_CO3_OVER_FE_II
    real(kind = DBL_PREC), dimension(nkn) :: FE_OH_2_OVER_FE_II
    real(kind = DBL_PREC), dimension(nkn) :: FE_S_OVER_FE_II
    real(kind = DBL_PREC), dimension(nkn) :: FE_S_2_OVER_FE_II

    real(kind = DBL_PREC), dimension(nkn) :: FE_3_O_4_OVER_FE_II
    real(kind = DBL_PREC), dimension(nkn) :: FE_2_O_3_OVER_FE_II
    real(kind = DBL_PREC), dimension(nkn) :: FE_OOH_OVER_FE_II
    real(kind = DBL_PREC), dimension(nkn) :: FE_OH_3_OVER_FE_II

    real(kind = DBL_PREC), dimension(nkn) :: MN_CO3_OVER_MN_II
    real(kind = DBL_PREC), dimension(nkn) :: MN_OH_2_OVER_MN_II
    real(kind = DBL_PREC), dimension(nkn) :: MN_S_OVER_MN_II
    
    real(kind = DBL_PREC), dimension(nkn, 4) :: FE_II_ACTIVITY_RATIOS
    real(kind = DBL_PREC), dimension(nkn, 3) :: FE_III_ACTIVITY_RATIOS
    real(kind = DBL_PREC), dimension(nkn, 3) :: MN_II_ACTIVITY_RATIOS
    
    integer, dimension(nkn) :: FE_II_SALT_NO
    integer, dimension(nkn) :: FE_III_SALT_NO
    integer, dimension(nkn) :: MN_II_SALT_NO

    real(kind = DBL_PREC), dimension(nkn) :: FREE_FE_II
    real(kind = DBL_PREC), dimension(nkn) :: FREE_FE_III
    real(kind = DBL_PREC), dimension(nkn) :: FREE_MN_II

    real(kind = DBL_PREC), dimension(nkn) :: TEMP_K
    real(kind = DBL_PREC), dimension(nkn) :: LOG_KW
    real(kind = DBL_PREC), dimension(nkn) :: KW
    real(kind = DBL_PREC), dimension(nkn) :: OH_MINUS
    real(kind = DBL_PREC), dimension(nkn) :: FE_III_FREE
    real(kind = DBL_PREC), dimension(nkn) :: H_PLUS_OVER_2
    real(kind = DBL_PREC), dimension(nkn) :: H_PLUS_OVER_3
    real(kind = DBL_PREC), dimension(nkn) :: H_PLUS_OVER_4
    
    
    real(kind = DBL_PREC) :: K1_FE_OH_3      
    real(kind = DBL_PREC) :: BETA_2_FE_OH_3  
    real(kind = DBL_PREC) :: BETA_3_FE_OH_3  
    real(kind = DBL_PREC) :: BETA_4_FE_OH_3  
    real(kind = DBL_PREC) :: BETA_2_2_FE_OH_3
    real(kind = DBL_PREC) :: BETA_4_3_FE_OH_3
    
    H_PLUS = 10.0D0 ** (-PH)
    
    LIM_DOXY_RED = DOXY  / (DOXY + K_HS_DOXY_RED_LIM)

    LIM_NO3N_RED = (NO3N / (NO3N + K_HS_NO3N_RED_LIM)) * &
        (K_HS_DOXY_RED_INHB / (DOXY + K_HS_DOXY_RED_INHB))

    LIM_MN_IV_RED = (MN_IV  / (MN_IV + K_HS_MN_IV_RED_LIM)) * &
        (K_HS_DOXY_RED_INHB / (DOXY + K_HS_DOXY_RED_INHB))  * &
        (K_HS_NO3N_RED_INHB / (NO3N + K_HS_NO3N_RED_INHB))

    LIM_FE_III_RED = (FE_III  / (FE_III + K_HS_FE_III_RED_LIM)) * &
        (K_HS_DOXY_RED_INHB  / (DOXY  + K_HS_DOXY_RED_INHB))    * &
        (K_HS_NO3N_RED_INHB  / (NO3N  + K_HS_NO3N_RED_INHB))    * &
        (K_HS_MN_IV_RED_INHB / (MN_IV + K_HS_MN_IV_RED_INHB))

    LIM_S_PLUS_6_RED = (S_PLUS_6 / (S_PLUS_6 + K_HS_S_PLUS_6_RED_LIM)) * &
        (K_HS_DOXY_RED_INHB   / (DOXY   + K_HS_DOXY_RED_INHB))  * &
        (K_HS_NO3N_RED_INHB   / (NO3N   + K_HS_NO3N_RED_INHB))  * &
        (K_HS_MN_IV_RED_INHB  / (MN_IV  + K_HS_MN_IV_RED_INHB)) * &
        (K_HS_FE_III_RED_INHB / (FE_III + K_HS_FE_III_RED_INHB))

    LIM_DOC_RED = 1.0D0 - &
        (LIM_DOXY_RED + LIM_NO3N_RED + LIM_MN_IV_RED + LIM_FE_III_RED + LIM_S_PLUS_6_RED)

    where (LIM_DOC_RED < 0.0D0)
        LIM_DOC_RED = 0.0D0
    end where

    REDUCER_LIM_FACTORS(:,1) = LIM_DOXY_RED
    REDUCER_LIM_FACTORS(:,2) = LIM_NO3N_RED
    REDUCER_LIM_FACTORS(:,3) = LIM_MN_IV_RED
    REDUCER_LIM_FACTORS(:,4) = LIM_FE_III_RED
    REDUCER_LIM_FACTORS(:,5) = LIM_S_PLUS_6_RED
    REDUCER_LIM_FACTORS(:,6) = LIM_DOC_RED

    REDUCED_AGENT_NO = maxloc(REDUCER_LIM_FACTORS, dim = 2);

    call DO_SATURATION_VEC(TEMP, SALT, ELEVATION, nkn, CS);

    ! Dissolved oxygen is reduced
    where (REDUCED_AGENT_NO == 1)
        PE = 20.75D0 - log10(1.0D0 / (((0.21D0 * (DOXY/CS))**0.25D0) * H_PLUS))
    end where

    ! Nitrate is reduced
    where (REDUCED_AGENT_NO == 2)
        PE = 21.05 - log10(1.0D0 / ((NO3N/14000.0D0)*(H_PLUS**1.2D0)))
    end where

    ! Mn IV is reduced
    where (REDUCED_AGENT_NO == 3)
        PE = 20.8D0 - log10(((MN_II/54938.0D0) ** 0.5D0) / (((MN_IV/54938.0D0) ** 0.5D0)*(H_PLUS**2.0D0)))
    end where

    ! FE III is reduced
    where (REDUCED_AGENT_NO == 4)
        PE = 13.0D0 - log10(FE_II / FE_III)
    end where

    ! S VI is reduced
    where (REDUCED_AGENT_NO == 5)
        HS_MOLAR = (S_MINUS_2 / 32000.0D0) * &
            ((H_PLUS * 8.9D-8)  / ((H_PLUS * H_PLUS) + (H_PLUS * 8.9D-8) + (8.9D-8 * 1.2D-13)))

        PE = 4.25D0 - log10((HS_MOLAR**0.125D0) / (((S_PLUS_6 / 32000.0D0) **0.125D0) * (H_PLUS**1.125D0)));
    end where

    ! Methanogenesis
    where (REDUCED_AGENT_NO == 6)
        PE = -0.2D0 - log10(1.0D0 / ((DISS_ORG_C / 12000D0)**0.25D0) * H_PLUS)
    end where

    ! Since PE is known, let's have a look in FE_II and FE_III solids following
    ! the Appendix 8.1, pp 513-515 in Stumm and Morgen (1996),
    ! "Aquatic Chemistry, Chemical Equilibria and Rates in Natural Waters"

    ! FE_CO3  / DISS_FE_II    
    FE_CO3_OVER_FE_II  = 10.0D0**(-0.2D0 +  PH + log10(HCO3))

    ! FE_OH_2 / DISS_FE_II
    FE_OH_2_OVER_FE_II = 10.0D0**(-11.7D0 + (2.0D0*PH))

    ! FE_S / DISS_FE_II
    FE_S_OVER_FE_II    = 10.0D0**(38.0D0  - (8.0D0*PH) + log10((S_PLUS_6 / 32000.0D0)) - (8.0D0*PE))

    ! FE_S_2 / DISS_FE_II
    FE_S_2_OVER_FE_II  = 10.0D0**(86.8D0  - (16.0D0*PH) + (2.0D0*log10((S_PLUS_6 / 32000.0D0))) - (14.0D0*PE))

    
    ! Now find out which Fe II salt is more likely to form for each reactor
    FE_II_ACTIVITY_RATIOS (:, 1) = FE_CO3_OVER_FE_II
    FE_II_ACTIVITY_RATIOS (:, 2) = FE_OH_2_OVER_FE_II
    FE_II_ACTIVITY_RATIOS (:, 3) = FE_S_OVER_FE_II
    FE_II_ACTIVITY_RATIOS (:, 4) = FE_S_2_OVER_FE_II    
    
    FE_II_SALT_NO  = maxloc(FE_II_ACTIVITY_RATIOS , dim = 2)

    !FeCO3
    where(FE_II_SALT_NO == 1)
        FREE_FE_II = 10.0D0**(-0.3D0 - PH + log10(HCO3))
    end where

    !Fe(OH)2
    where(FE_II_SALT_NO == 2)
        FREE_FE_II = 10.0D0**(13.3D0 - (2.0D0*PH))
    end where

    S_MINUS_2_MOLAR = (S_MINUS_2 / 32000.0D0) * &
        ((8.9D-8 * 1.2D-13)  / ((H_PLUS * H_PLUS) + (H_PLUS * 8.9D-8) + (8.9D-8 * 1.2D-13)))

    !FeS
    where((FE_II_SALT_NO == 3).and.(S_MINUS_2_MOLAR > 1.0D-12) )        
        FREE_FE_II = 10.0D0**(-18.64D0) / S_MINUS_2_MOLAR
    end where

    where((FE_II_SALT_NO == 3).and.(S_MINUS_2_MOLAR <= 1.0D-12) )        
        FREE_FE_II = FE_II / 56000.0D0
    end where
    
    
    !FeS2
    where((FE_II_SALT_NO == 4).and.(S_MINUS_2_MOLAR > 1.0D-12) )        
        FREE_FE_II = (10.0D0**(-26.89D0)) / (4.0D0 * S_MINUS_2_MOLAR * S_MINUS_2_MOLAR)
    end where

    where((FE_II_SALT_NO == 4).and.(S_MINUS_2_MOLAR <= 1.0D-12) )        
        FREE_FE_II = FE_II / 56000.0D0
    end where

    ! For a while no complex formation (to be fixed)
    FE_II_DISS = FREE_FE_II
    
    ! -------------------------------------------------------------------------
    ! END OF FE_II Species
    ! -------------------------------------------------------------------------


    ! Calculate the complexation inn case of iron III Hydroxides
    K1_FE_OH_3       = 10.0D0**(-3.05D0)
    BETA_2_FE_OH_3   = 10.0D0**(-6.31D0)
    BETA_3_FE_OH_3   = 10.0D0**(-13.8D0)
    BETA_4_FE_OH_3   = 10.0D0**(-22.7D0)
    BETA_2_2_FE_OH_3 = 10.0D0**(-2.91D0)
    BETA_4_3_FE_OH_3 = 10.0D0**(-5.77D0)

    !where((FE_III_SALT_NO == 2).or.(FE_III_SALT_NO == 3))
    !for a while assume that iron 3+ solubity is only related to Fe(OH3)
        TEMP_K = TEMP + 273.15D0
        
        LOG_KW = &
            -2.839710D2 + (1.3323D4 / TEMP_K) - (5.069842D-2 * TEMP_K) + &
            (1.0224447D2 * log10(TEMP_K)) - (1.119669D6 / (TEMP_K * TEMP_K))
        
        KW = 10**(LOG_KW)
        OH_MINUS = KW / H_PLUS
        !FREE_FE_III = (6.0D0**(-38.0D0))/(OH_MINUS * OH_MINUS * OH_MINUS)
        FE_III_FREE = 10.0D0**(3.96D0 - (3.0D0 * PH))
        
        H_PLUS_OVER_2 = H_PLUS * H_PLUS
        H_PLUS_OVER_3 = H_PLUS * H_PLUS * H_PLUS
        H_PLUS_OVER_4 = H_PLUS * H_PLUS * H_PLUS * H_PLUS

        FE_III_DISS = FE_III_FREE * &
            (1.0D0 + (K1_FE_OH_3 / (H_PLUS))  + (BETA_2_FE_OH_3 / H_PLUS_OVER_2) + &
             (BETA_3_FE_OH_3 / H_PLUS_OVER_3) + (BETA_4_FE_OH_3 / H_PLUS_OVER_4) + &
             ((2.0D0*BETA_2_2_FE_OH_3*FE_III_FREE)/H_PLUS_OVER_2) + &
             ((3.0D0*BETA_4_3_FE_OH_3*FE_III_FREE*FE_III_FREE)/H_PLUS_OVER_4));
    !end where
    
    if(is_nan(FE_III_DISS)) then
    !if(any(FE_III_DISS /= FE_III_DISS)) then
    
       print *,'REDOX_AND_SPECIATION:'
       print *,'FE_III_DISS is NaN:', FE_III_DISS
       print *,'FREE_FE_III:'   , FREE_FE_III
       print *,'K1_FE_OH_3:'    , K1_FE_OH_3
       print *,'BETA_2_FE_OH_3:', BETA_2_FE_OH_3
       print *,'H_PLUS_OVER_3:' , H_PLUS_OVER_3
       print *,'BETA_4_FE_OH_3:', BETA_4_FE_OH_3
       print *,'H_PLUS_OVER_4:' , H_PLUS_OVER_4
       print *,'FE_III_FREE:'   , FE_III_FREE         
       stop
    end if 
    
    
    
    ! -------------------------------------------------------------------------
    ! END OF FE_III Species
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! MN_II Species
    ! -------------------------------------------------------------------------

    ! MN_CO3  / DISS_MN_II
    MN_CO3_OVER_MN_II  = 10.0D0**(-0.2D0 +  PH + log10(HCO3))

    ! MN_OH_2 / DISS_MN_II
    MN_OH_2_OVER_MN_II = 10.0D0**(-15.0D0 + (2.0D0*PH))

    ! MN_S / DISS_MN_II
    MN_S_OVER_MN_II    = 10.0D0**(34.0D0  - (8.0D0*PH) + log10((S_PLUS_6 / 32000.0D0)) - (8.0D0*PE))
    
    ! Now find out which Mn II salt is more likely to form for each reactor
    MN_II_ACTIVITY_RATIOS (:, 1) = MN_CO3_OVER_MN_II
    MN_II_ACTIVITY_RATIOS (:, 2) = MN_OH_2_OVER_MN_II
    MN_II_ACTIVITY_RATIOS (:, 3) = MN_S_OVER_MN_II

    MN_II_SALT_NO  = maxloc(MN_II_ACTIVITY_RATIOS , dim = 2)

    !MnCO3
    where(MN_II_SALT_NO == 1)
        FREE_MN_II = (10.0D0 ** (8.03D0)) / CO3 
    end where

    !Mn(OH)2
    where(FE_II_SALT_NO == 2)
        FREE_MN_II = (10.0D0 ** (11.14D0))/ (4 * OH_MINUS * OH_MINUS)
    end where

    !MnS
    where((MN_II_SALT_NO == 3).and.(S_MINUS_2_MOLAR > 1.0D-12))
        FREE_MN_II = 10.0D0**(-10.19D0) / S_MINUS_2_MOLAR
    end where

    where((MN_II_SALT_NO == 3).and.(S_MINUS_2_MOLAR <= 1.0D-12))
        FREE_MN_II = MN_II / 54938.0D0
    end where
    
    ! For a while no complex formation (to be fixed)
    MN_II_DISS = FREE_MN_II
    ! -------------------------------------------------------------------------
    ! END OF MN_II Species
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! MN_IV Species
    ! -------------------------------------------------------------------------
    ! We assume that MnO2 is the main MN 4+ species
    ! -------------------------------------------------------------------------
    ! END OF MN_IV Species
    ! -------------------------------------------------------------------------
end subroutine REDOX_AND_SPECIATION        
        


!Function, which returns saturation concentration of dissolved oxygen
subroutine DO_SATURATION_VEC(T, S, H, nkn, CS)

    !Water temperature (in Celcius)
    double precision, dimension(nkn), intent(in) :: T

    !Salinity (in ppt)
    double precision, dimension(nkn), intent(in) :: S

    !Elevation (in m)
    double precision, dimension(nkn), intent(in) :: H

    integer, intent(in) :: nkn

    !Water temperature (in Kelvin)
    double precision, dimension(nkn) :: T_KELVIN

    !Altitude (in feet)
    double precision, dimension(nkn) :: H_FEET

    double precision, dimension(nkn) :: LN_CSF
    double precision, dimension(nkn) :: LN_CSS
    double precision, dimension(nkn) :: CSS
    double precision, dimension(nkn) :: CSP

    double precision, dimension(nkn) :: CS
    !Pressure at altitude H (in atm)
    double precision, dimension(nkn) :: P

    !Standart pressure (in mmHg)
    double precision, dimension(nkn) :: P0

    double precision, dimension(nkn) :: LN_PWV

    !Partial pressure of water vapor (in atm)
    double precision, dimension(nkn) :: PWV

    !A constant
    double precision, dimension(nkn) :: THETA

    T_KELVIN = T + 273.15
    H_FEET = H / 0.3048D0

    !Calculate the effect of temperature on dissolved oxygen saturation
    LN_CSF = -139.34411D0 + (157570.1d0 / T_KELVIN) - &
    &        (66423080.0D0     / (T_KELVIN * T_KELVIN)) + &
    &        (12438000000.0D0  / (T_KELVIN * T_KELVIN * T_KELVIN)) - &
    &        (862194900000.0D0 / (T_KELVIN * T_KELVIN * T_KELVIN * T_KELVIN))

    !Calculate the effect of salinity on dissolved oxygen saturation
    LN_CSS = LN_CSF - S * &
    &        (0.017674D0 - (10.754 / T_KELVIN) + &
    &         (2140.7D0 / (T_KELVIN ** 2.0D0)))

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
    end subroutine DO_SATURATION_VEC

! -----------------------------------------------------------------------------
! The subroutine below calculates the Iron II dissolution according to
! Stumm and Lee (1960).    
! -----------------------------------------------------------------------------
!
!                       Initial development 2 nd of July 2016
!    
!                                by Ali Ertürk
! -----------------------------------------------------------------------------
subroutine IRON_II_DISSOLUTION(HS2_TOT, PH, TOT_ALK, nkn, FE_II_TOT)
    use AQUABC_II_GLOBAL
    implicit none

    integer, intent(in) ::  nkn

    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: HS2_TOT       ! [H2S] + [HS-]  +  [S]   (mol/L)
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PH            ! PH
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: TOT_ALK       ! Total alkalinity (mol/L)
    
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: FE_II_TOT     ! Total Fe2+ (mol/L)
    
    
    real(kind = DBL_PREC), dimension(nkn) :: K_1      ! Equilibrium constant for the reaction : Fe(OH)2    <-------->    Fe2+    +     2(OH-) 
    real(kind = DBL_PREC), dimension(nkn) :: K_2      ! Equilibrium constant for the reaction : Fe(OH)2    <-------->    [Fe(OH)]+     +    OH- 
    real(kind = DBL_PREC), dimension(nkn) :: K_3      ! Equilibrium constant for the reaction : Fe(OH)2    +    OH-    <-------->     [Fe(OH)3]- 
    real(kind = DBL_PREC), dimension(nkn) :: K_4      ! Equilibrium constant for the reaction : Fe(CO3)    <-------->    Fe++    +     (CO3)--
    real(kind = DBL_PREC), dimension(nkn) :: K_5      ! Equilibrium constant for the reaction : Fe(OH)2    +    OH-    <-------->     [Fe(OH)]+    +    CO3-- 
    real(kind = DBL_PREC), dimension(nkn) :: K_6      ! Equilibrium constant for the reaction : HCO3-      <-------->    H+      +    CO3--         
    real(kind = DBL_PREC), dimension(nkn) :: K_7      ! Equilibrium constant for the reaction : FeS        <-------->    Fe++    +    S--
    real(kind = DBL_PREC), dimension(nkn) :: K_8      ! Equilibrium constant for the reaction : FeS        +    OH-    <-------->     [Fe(OH)]+    +    S--
    real(kind = DBL_PREC), dimension(nkn) :: K_9      ! Equilibrium constant for the reaction : FeS        +  3(OH-)   <-------->     [Fe(OH3)]-   +    S--
    real(kind = DBL_PREC), dimension(nkn) :: K_10_A   ! Equilibrium constant for the reaction : H2S        +    H+     <-------->     HS-
    real(kind = DBL_PREC), dimension(nkn) :: K_10_B   ! Equilibrium constant for the reaction : HS-        +    H+     <-------->     S--
    real(kind = DBL_PREC), dimension(nkn) :: K_W      ! Equilibrium constant for water dissociation
    
    real(kind = DBL_PREC), dimension(nkn) :: H_PLUS               ! [H+]    (mol/L) 
    real(kind = DBL_PREC), dimension(nkn) :: OH_MINUS             ! [OH-]   (mol/L)
    
    real(kind = DBL_PREC), dimension(nkn,3) :: ALL_FE_II
    ! ALL_FE_II(:,1) : Fe(OH)2 solubility (mol/L)
    ! ALL_FE_II(:,2) : FeCO3   solubility (mol/L)
    ! ALL_FE_II(:,3) : FeS     solubility (mol/L) 

    ! For now, equlibrium constants are hard-coded. In the future, they could be
    ! new nodel constants. Another option to calculate more realistic values
    ! for the model constants is to use the temperature and thermodynamic constants.
    K_1    = 8.0D-16
    K_2    = 4.0D-10
    K_3    = 8.3D-6
    K_4    = 2.1D-11
    K_5    = 1.0D-5
    K_6    = 4.8D-11
    K_7    = 6.0D-18
    K_8    = 3.0D-12
    K_9    = 6.2D-8
    K_10_A = 1.0D-7
    K_10_B = 1.3D-13
    K_W    = 1.0D-14

    H_PLUS   = 10.0D0 ** (-PH)
    OH_MINUS = K_W / H_PLUS

    ALL_FE_II(:,1) = &
        ((K_1  / (K_W * K_W)) * H_PLUS * H_PLUS) + ((K_2 / K_W) * H_PLUS) + &
        ((K_3 * K_W) / H_PLUS)

    ALL_FE_II(:,2) = &
        ((H_PLUS + (2.0D0 * K_6)) / (TOT_ALK * K_6)) * (K_4 + ((K_5 * K_W) / H_PLUS))
    
    ALL_FE_II(:,3) = 1.0D0
    
    where (HS2_TOT < 1.0D-20)
        HS2_TOT   = 1.0D-20
        FE_II_TOT = minval(ALL_FE_II(:,1:2), dim=2)
    elsewhere
        ALL_FE_II(:,3)  = &
            ((K_7 / HS2_TOT) * &
             (1.0D0 + (H_PLUS / K_10_B) + ((H_PLUS * H_PLUS) / (K_10_B * K_10_A)))) + &
             (K_7 * ((OH_MINUS / K_8) + ((27.0D0 * OH_MINUS * OH_MINUS * OH_MINUS) / K_9)))

      FE_II_TOT(:) = minval(ALL_FE_II, dim=2)
    end where
    
end subroutine IRON_II_DISSOLUTION

    
    
! -----------------------------------------------------------------------------
! The subroutine below calculates the Iron II oxidation according to
! Morgen and Lahav (2007).    
! -----------------------------------------------------------------------------
!
!                       Initial development 6 th of July 2016
!    
!                                by Ali Ertürk
! -----------------------------------------------------------------------------
subroutine IRON_II_OXIDATION(FE_II_DISS, DOXY, PH, TEMP, SALT, ELEVATION, nkn, R_FE_II_OXIDATION)
    use AQUABC_II_GLOBAL
    implicit none

    integer, intent(in) ::  nkn

    real(kind = DBL_PREC), dimension(nkn), intent(in) :: FE_II_DISS    ! Total dissolved Fe2+   (mg/L)
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PH            ! PH
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: DOXY          ! Dissolved oxygen
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: SALT
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: ELEVATION
    
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: R_FE_II_OXIDATION     ! Rate of Fe2+ oxidation (mg/L/day)
    
    ! Use environmental conditions to estimate the values of following rate and equilibrium constants
    real(kind = DBL_PREC), dimension(nkn) :: K_1           ! Equilibrium constant K_1 
    real(kind = DBL_PREC), dimension(nkn) :: K_2           ! Equilibrium constant K_2
    real(kind = DBL_PREC), dimension(nkn) :: K_3           ! Equilibrium constant K_3
    real(kind = DBL_PREC), dimension(nkn) :: K_W           ! Equilibrium constant for water dissociation
    real(kind = DBL_PREC), dimension(nkn) :: k_FE_II       ! Rate constant for Fe2+ oxidation
    real(kind = DBL_PREC), dimension(nkn) :: k_FE_OH_PLUS  ! Rate constant for Fe2(OH)+ oxidation
    real(kind = DBL_PREC), dimension(nkn) :: k_FE_OH_2     ! Rate constant for Fe2(OH)2 oxidation
    
    ! Auxillary variables
    real(kind = DBL_PREC), dimension(nkn) :: H_PLUS        ! [H+]    (mol/L) 
    real(kind = DBL_PREC), dimension(nkn) :: OH_MINUS      ! [OH-]   (mol/L)
    real(kind = DBL_PREC), dimension(nkn) :: FE_II_TOT     ! Fe2+Tot (mol/L)
    real(kind = DBL_PREC), dimension(nkn) :: CS
    
    ! For now, equlibrium constants are hard-coded. In the future, they could be
    ! new nodel constants. Another option to calculate more realistic values
    ! for the model constants is to use the temperature and thermodynamic constants.
    K_1    = 10.0D0**(4.50D0)
    K_2    = 10.0D0**(2.93D0)
    K_3    = 10.0D0**(3.57D0)

    K_W    = 1.0D-14

    k_FE_II      = 6.0D-5
    k_FE_OH_PLUS = 1.7D0
    k_FE_OH_2    = 4.3D5

    H_PLUS    = 10.0D0 ** (-PH)
    OH_MINUS  = K_W / H_PLUS
    FE_II_TOT = FE_II_DISS / 56000.0D0
    
    call DO_SATURATION_VEC(TEMP, SALT, ELEVATION, nkn, CS);
    
    R_FE_II_OXIDATION = &
        ((k_FE_II      / (1.0D0 + ((K_1*K_W) / (H_PLUS))  + &
                           ((K_1*K_2*K_W*K_W) / (H_PLUS*H_PLUS)) + &
                           ((K_1*K_2*K_W*K_W*K_W*K_W) / (H_PLUS*H_PLUS*H_PLUS)))) + & 
         (k_FE_OH_PLUS / (1.0D0 + ((H_PLUS)  / (K_1*K_W)) + &
                          ((K_2*K_W) / (H_PLUS))  +  &
                          ((K_2*K_3*K_W*K_W) / (H_PLUS*H_PLUS))))  + &
         (k_FE_OH_2    / (1.0D0 + ((K_3*K_W) / (H_PLUS))  +  &
                          ((H_PLUS) / (K_2*K_W)) + &
                           ((H_PLUS*H_PLUS) / (K_1*K_2*K_W*K_W))))) * FE_II_TOT
    
    R_FE_II_OXIDATION = (R_FE_II_OXIDATION * (DOXY / CS)) * 56000.0D0
end subroutine IRON_II_OXIDATION


! -----------------------------------------------------------------------------
! The subroutine below calculates the inorganic phosphorus dissolution according to
! Veroni Snoeyink and Jenkins, 1980.    
! -----------------------------------------------------------------------------
!
!                       Initial development 6 th of July 2016
!    
!                                by Ali Ertürk
! -----------------------------------------------------------------------------
subroutine IP_SOLUBLE_FRACTION &
            (FE_III, PO4P, K_A_1, K_A_2, K_A_3, PH, nkn, nlayers, DIP_OVER_IP)
    
    use AQUABC_II_GLOBAL
    implicit none

    integer, intent(in) ::  nkn, nlayers

    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(in) :: FE_III
    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(in) :: K_A_1    ! First dissociation constant for H3PO4
    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(in) :: K_A_2    ! Second dissociation constant for H3PO4
    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(in) :: K_A_3    ! Third dissociation constant for H3PO4
    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(in) :: PO4P     ! Total inorganic phosphorus (mg/L)
    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(in) :: PH
    
    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(inout) :: DIP_OVER_IP
    
    ! Auxillary variables
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: C_T_FE_III ! [Fe3+]tot (moles/L)
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: COEFF_1    ! Intermediate ceofficient as a function of [H+], K_1, K_4
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: COEFF_2    ! Intermediate ceofficient as a function of [H+], K_A_1, K_A_2, K_A_3
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: KS0        ! Equlibrium constant for AlPO4 <-----> Al3+  +  PO4---
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: H_PLUS     ! [H+]    (mol/L)
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: C_T_PO4    ! PO4 solubility in moles
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: K_1        ! Equlibrium constant for Al+++    +     H2O    <-------->    Al(OH)++    +     H+
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: K_4        ! Equlibrium constant for Al+++    +    4H2O    <-------->    Al(OH)4-    +     4H+
    
    KS0        = 1.0D-21
    K_1        = 1.0D-5
    K_4        = 10.0D0 ** (-21.7D0)
    H_PLUS     = 10.0D0 ** (-PH)
    C_T_FE_III = FE_III / 56000.0D0
    COEFF_1    = 1.0D0 + (K_1/H_PLUS) + (K_4/(H_PLUS*H_PLUS*H_PLUS*H_PLUS))
    COEFF_2    = 1.0D0 + ((H_PLUS*H_PLUS*H_PLUS)/(K_A_1*K_A_2*K_A_3)) + ((H_PLUS*H_PLUS)/(K_A_2*K_A_3)) + (H_PLUS/K_A_3)
    
    ! Main equation
    ! KS0 = (C_T_FE_III / COEFF_1) * (C_T_PO4 / COEFF_2)
    !   ==> C_T_PO4 = (KS0 * COEFF_1 * COEFF_2) / C_T_FE_III 
    C_T_PO4 = (KS0 * COEFF_1 * COEFF_2) / C_T_FE_III    
    
    DIP_OVER_IP = C_T_PO4 / C_T_FE_III 
    
    where(DIP_OVER_IP > 1.0D0)
        DIP_OVER_IP = 1.0D0
    end where

end subroutine IP_SOLUBLE_FRACTION

!********************************************

! ---------------------------------------------------------------------------------------
! Subroutine to calculate the dissolved metal concentration at the end of the timestep
! and averaged over the timestep using some simplified dissolution kinetics, where
! dissolution is directly propotional to the concentration gradient from equlibrium
! solubility. If more metal is dissolved than the equibirium solubility, then
! the equation will work in a reverse way.
! ---------------------------------------------------------------------------------------
!                  Developer        : Ali Erturk
!                  Development date : 9 August 2016
! ---------------------------------------------------------------------------------------
!
! INPUT VARIABLES
!     TOT_ME              : Total concentration of the metal (FE_II, FE_III, CA, MG, etc.)
!     ME_DISS_INIT        : Initial dissolved metal concentration
!     ME_SOLUB_EQ         : Estimated solubility
!     k_DISS_ME           : Dissolution precipitation rate constant
!     t                   : Time (to be taken as time step)
!     nkn                 : Number of nodes
!     nlayers             : Number of layers
!
! OUTPUT VARIABLES
!     DISS_ME_CONC_TS_END : Dissolved metal concentration at the end of timestep
!     DISS_ME_CONC_TS_AVG : Dissolved metal concentration averaged over the timestep
! ---------------------------------------------------------------------------------------
subroutine CALC_DISS_ME_CONC &
           (TOT_ME             , &! Pass from the state variable representing the total metal concentration (FE_II, FE_III, CA, etc.)
            ME_DISS_INIT       , &! Pass from the last time step (SAVED_OUTPUTS) 
            ME_SOLUB_EQ        , &! Pass from equlibrum chemistry calculations (for example FE_II_DISS_EQ)
            k_DISS_ME          , &! Rate constant
            t                  , &! This is the time step
            nkn                , &
            nlayers            , &
            DISS_ME_CONC_TS_END, &
            DISS_ME_CONC_TS_AVG)

    use AQUABC_II_GLOBAL
    implicit none

    ! In going variables
    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(in) :: TOT_ME
    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(in) :: ME_DISS_INIT
    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(in) :: ME_SOLUB_EQ

    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(in) :: k_DISS_ME
    real(kind = DBL_PREC), intent(in) :: t

    integer, intent(in) :: nkn
    integer, intent(in) :: nlayers

    ! Out going variables
    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(inout) :: DISS_ME_CONC_TS_END
    real(kind = DBL_PREC), dimension(nkn, nlayers), intent(inout) :: DISS_ME_CONC_TS_AVG

    ! Auxillary variables
    !real(kind = DBL_PREC), dimension(nkn, nlayers) :: A_2
    !real(kind = DBL_PREC), dimension(nkn, nlayers) :: A_3
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: C
    !real(kind = DBL_PREC), dimension(nkn, nlayers) :: Y_0
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: INT_ME_DISS_t   
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: INT_ME_DISS_zero
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: LOG_INT_ME_DISS_t   
    real(kind = DBL_PREC), dimension(nkn, nlayers) :: LOG_INT_ME_DISS_zero
    

    ! Updated by Ali and Petras, 15 th of August 2016. 
    
    !Y_0 = ME_DISS_INIT
    
    where (ME_DISS_INIT > ME_SOLUB_EQ)
        ! This is the oversaturated case, so dissolution reaction will move oppsite way to
        ! chemical precipitation. In this case the rate is dependent on dissolved metal as
        ! well and we end up with nonlinear equation of second degree due to dissolved metal
        ! (Bernoulli type ODE) and is still analytically solveable. 
        
        ! Function              : ME_DISS(t) = -k_DISS_ME * (ME_DISS - ME_SOLUB_EQ) * (ME_DISS - ME_SOLUB_EQ)
        
        !In future more advanced model considering cyristal or particle growth will be developed.

        DISS_ME_CONC_TS_END = ME_SOLUB_EQ + (1.0D0 / (k_DISS_ME + (1.0D0/(ME_DISS_INIT - ME_SOLUB_EQ)))) 

        LOG_INT_ME_DISS_t    = log(dabs((1.0D0 / (ME_DISS_INIT - ME_SOLUB_EQ)) + (k_DISS_ME * t)))
        LOG_INT_ME_DISS_zero = log(dabs((1.0D0 / (ME_DISS_INIT - ME_SOLUB_EQ))                  ))
        
        INT_ME_DISS_t    =  (ME_SOLUB_EQ * t) + ((1.0D0 / (k_DISS_ME)) * LOG_INT_ME_DISS_t) 
        INT_ME_DISS_zero =  (1.0D0 / (k_DISS_ME)) * LOG_INT_ME_DISS_zero
        
        DISS_ME_CONC_TS_AVG = (1.0D0/t) * (INT_ME_DISS_t - INT_ME_DISS_zero)        
        
        ! Old code
        !A_2 = k_DISS_ME * ME_SOLUB_EQ
        !A_3 = k_DISS_ME
        !C   = (log(Y_0) / A_2) - (log(A_2 - (A_3 * Y_0)) / A_2)
        
        ! Calculate the dissolved metal at the end of timestep        
        !DISS_ME_CONC_TS_END = (A_2 * exp(A_2 * (t + C))) / ((A_3 * exp(A_2 * (t + C))) + 1.0D0)

        ! Calculate the definite integral for the solution and divide it by time to
        ! get dissolved metal averaged over timestep  
        !DISS_ME_CONC_TS_AVG  = (1.0D0/t) * &
        !    (((1.0D0 * A_3)*(log(dabs((1.0D0/A_3)*((A_3*exp(C*A_2)*exp(t*A_2))+1.0D0)))-(C * A_2))) - &
        !    ((1.0D0 * A_3)*(log(dabs((1.0D0/A_3)*((A_3*exp(C*A_2))+1.0D0)))-(C * A_2))))

    else where (ME_DISS_INIT < ME_SOLUB_EQ)
        ! This is the undersaturated case where more metals will be dissolved. 
        ! In this case the rate is dependent on particulate metal as
        ! well and we end up with nonlinear equation of second degree due to dissolved metal
        ! (Ricatti type ODE) and is still analytically solveable.
        
        ! Function              : ME_DISS(t)
        ! Differential equation : diff(ME_DISS) == k_DISS_ME * (ME_SOLUB_EQ - ME_DISS) * (TOT_ME - ME_DISS - ME_SOLUB_EQ)
        ! Initial condition     : ME_DISS(0) == ME_DISS_INIT)
        
        ! In future more advanced model considering dissolution related to particle
        ! size and more environmental colditions will be considered.

        ! Calculate the dissolved metal at the end of timestep 
        DISS_ME_CONC_TS_END = ME_SOLUB_EQ + &
            ((3.0D0*ME_SOLUB_EQ)/ &
             (1.0D0 -(((ME_DISS_INIT-(4.0D0*ME_SOLUB_EQ))/(ME_DISS_INIT-ME_SOLUB_EQ))* &
                      exp(3.0D0*t*k_DISS_ME*ME_SOLUB_EQ))))

        ! Calculate the definite integral for the solution and divide it by time to
        ! get dissolved metal averaged over timestep
        
        ! 19th of August 2016 updates
        C = (ME_DISS_INIT - (4.0D0*ME_SOLUB_EQ)) / (ME_DISS_INIT - ME_SOLUB_EQ)
        
        LOG_INT_ME_DISS_t    = log(dabs(C*exp(3.0D0*k_DISS_ME*ME_SOLUB_EQ*t) - 1.0D0))
        LOG_INT_ME_DISS_zero = log(dabs(C - 1.0D0))
        
        INT_ME_DISS_t    = (ME_SOLUB_EQ * t) + (3.0D0*ME_SOLUB_EQ*(t - (LOG_INT_ME_DISS_t/(3.0D0*k_DISS_ME*ME_SOLUB_EQ))))
        INT_ME_DISS_zero = (3.0D0*ME_SOLUB_EQ*(-(LOG_INT_ME_DISS_zero/(3.0D0*k_DISS_ME*ME_SOLUB_EQ))))
        
        ! Old code
        !LOG_INT_ME_DISS_t    = (1.0D0/((4.0D0*ME_SOLUB_EQ) - ME_DISS_INIT)) * &
        !    (ME_DISS_INIT - ME_SOLUB_EQ + (((4.0D0*ME_SOLUB_EQ) - ME_DISS_INIT)*exp(3.0D0*k_DISS_ME*ME_SOLUB_EQ*t)))
        !
        !LOG_INT_ME_DISS_zero = (1.0D0/((4.0D0*ME_SOLUB_EQ)-ME_DISS_INIT))* &
        !    (ME_DISS_INIT - ME_SOLUB_EQ+((4.0D0*ME_SOLUB_EQ)-ME_DISS_INIT))
        !
        !INT_ME_DISS_t    = (-1.0D0/k_DISS_ME)*(log(dabs(LOG_INT_ME_DISS_t)) - (4.0D0*t*ME_SOLUB_EQ*k_DISS_ME))
        !INT_ME_DISS_zero = (-1.0D0/k_DISS_ME)*log(dabs(LOG_INT_ME_DISS_zero))

        DISS_ME_CONC_TS_AVG = (1.0D0/t) * (INT_ME_DISS_t - INT_ME_DISS_zero)
    else where
        DISS_ME_CONC_TS_END = ME_DISS_INIT
        DISS_ME_CONC_TS_AVG = ME_DISS_INIT             
    end where
    
    ! Check for saturation case where logarithms in both over and undersaturation solutions may give Nans
    ! and assume that neither dissolution nor precipitation occurs.
    ! where (is_nan(DISS_ME_CONC_TS_END).or.is_nan(DISS_ME_CONC_TS_AVG))
    !     DISS_ME_CONC_TS_END = ME_DISS_INIT
    !     DISS_ME_CONC_TS_AVG = ME_DISS_INIT
    ! end where
    DISS_ME_CONC_TS_AVG = 0.5 * (ME_DISS_INIT + DISS_ME_CONC_TS_END)
end subroutine CALC_DISS_ME_CONC
