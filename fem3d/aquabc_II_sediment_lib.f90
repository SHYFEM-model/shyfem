
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

! Content:
!subroutine SED_POC_DISSOLUTION
!subroutine SED_DOC_MINERALIZATION
!subroutine SED_REDOX_AND_SPECIATION
!subroutine DO_SATURATION_MAT
!subroutine SED_IRON_II_DISSOLUTION

subroutine SED_POC_DISSOLUTION &
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

    use AQUABC_II_GLOBAL
    implicit none

    real(kind = DBL_PREC), intent(in) :: K_OXIC_DISS_POC
    real(kind = DBL_PREC), intent(in) :: K_ANOXIC_DISS_POC
    real(kind = DBL_PREC), intent(in) :: THETA_DISS_POC
    real(kind = DBL_PREC), intent(in) :: KHS_DISS_POC
    real(kind = DBL_PREC), intent(in) :: DOXY_AT_ANOXIA

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: SED_POC
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: SED_DOXY
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: SED_TEMPS

    integer :: nkn
    integer :: NUM_SED_LAYERS

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: R_DISS_POC

    where (SED_DOXY.GE.DOXY_AT_ANOXIA)

        R_DISS_POC  = &
            K_OXIC_DISS_POC    * (THETA_DISS_POC  ** (SED_TEMPS - 2.0D1)) * &
            (SED_POC / (SED_POC + KHS_DISS_POC)) * SED_POC

    end where

    where (SED_DOXY.LT.DOXY_AT_ANOXIA)

        R_DISS_POC  = &
            K_ANOXIC_DISS_POC  * (THETA_DISS_POC  ** (SED_TEMPS - 2.0D1)) * &
            (SED_POC / (SED_POC + KHS_DISS_POC))  * SED_POC
    end where

end subroutine SED_POC_DISSOLUTION


subroutine SED_DOC_MINERALIZATION &
           (SED_K_MIN_DOC_DOXY_20       , &
            SED_K_MIN_DOC_NO3N_20       , &
            SED_K_MIN_DOC_MN_IV_20      , &
            SED_K_MIN_DOC_FE_III_20     , &
            SED_K_MIN_DOC_S_PLUS_6_20   , &
            SED_K_MIN_DOC_DOC_20        , &
            SED_THETA_K_MIN_DOC_DOXY    , &
            SED_THETA_K_MIN_DOC_NO3N    , &
            SED_THETA_K_MIN_DOC_MN_IV   , &
            SED_THETA_K_MIN_DOC_FE_III  , &
            SED_THETA_K_MIN_DOC_S_PLUS_6, &
            SED_THETA_K_MIN_DOC_DOC     , &
            SED_K_HS_DOC_MIN_DOXY       , &
            SED_K_HS_DOC_MIN_NO3N       , &
            SED_K_HS_DOC_MIN_MN_IV      , &
            SED_K_HS_DOC_MIN_FE_III     , &
            SED_K_HS_DOC_MIN_S_PLUS_6   , &
            SED_K_HS_DOC_MIN_DOC        , &
            SED_K_HS_DOXY_RED_LIM       , &
            SED_K_HS_NO3N_RED_LIM       , &
            SED_K_HS_MN_IV_RED_LIM      , &
            SED_K_HS_FE_III_RED_LIM     , &
            SED_K_HS_S_PLUS_6_RED_LIM   , &
            SED_K_HS_DOXY_RED_INHB      , &
            SED_K_HS_NO3N_RED_INHB      , &
            SED_K_HS_MN_IV_RED_INHB     , &
            SED_K_HS_FE_III_RED_INHB    , &
            SED_K_HS_S_PLUS_6_RED_INHB  , &
            SED_PH_MIN_DOC_MIN_DOXY     , &
            SED_PH_MIN_DOC_MIN_NO3N     , &
            SED_PH_MIN_DOC_MIN_MN_IV    , &
            SED_PH_MIN_DOC_MIN_FE_III   , &
            SED_PH_MIN_DOC_MIN_S_PLUS_6 , &
            SED_PH_MIN_DOC_MIN_DOC      , &
            SED_PH_MAX_DOC_MIN_DOXY     , &
            SED_PH_MAX_DOC_MIN_NO3N     , &
            SED_PH_MAX_DOC_MIN_MN_IV    , &
            SED_PH_MAX_DOC_MIN_FE_III   , &
            SED_PH_MAX_DOC_MIN_S_PLUS_6 , &
            SED_PH_MAX_DOC_MIN_DOC      , &
            SED_TEMPS                   , &
            SED_DOC                     , &
            SED_DOXY                    , &
            SED_NO3N                    , &
            MN_IV_DISS                  , &
            FE_III_DISS                 , &
            S_PLUS_6                    , &
            PH                          , &
            nkn                         , &
            NUM_SED_LAYERS              , &
            PH_CORR_DOC_MIN_DOXY        , &
            PH_CORR_DOC_MIN_NO3N        , &
            PH_CORR_DOC_MIN_MN_IV       , &
            PH_CORR_DOC_MIN_FE_III      , &
            PH_CORR_DOC_MIN_S_PLUS_6    , &
            PH_CORR_DOC_MIN_DOC         , &
            LIM_DOXY_RED                , &
            LIM_NO3N_RED                , &
            LIM_MN_IV_RED               , &
            LIM_FE_III_RED              , &
            LIM_S_PLUS_6_RED            , &
            LIM_DOC_RED                 , &
            K_NO3_RED                   , &
            K_MN_IV_RED                 , &
            K_FE_III_RED                , &
            K_S_PLUS_6_RED              , &
            K_DOC_RED                   , &
            R_MINER_DOC_DOXY            , &
            R_MINER_DOC_NO3N            , &
            R_MINER_DOC_MN_IV           , &
            R_MINER_DOC_FE_III          , &
            R_MINER_DOC_S_PLUS_6        , &
            R_MINER_DOC_DOC)

    use AQUABC_II_GLOBAL
    implicit none

    real(kind = DBL_PREC), intent(in) :: SED_K_MIN_DOC_DOXY_20
    real(kind = DBL_PREC), intent(in) :: SED_K_MIN_DOC_NO3N_20
    real(kind = DBL_PREC), intent(in) :: SED_K_MIN_DOC_MN_IV_20
    real(kind = DBL_PREC), intent(in) :: SED_K_MIN_DOC_FE_III_20
    real(kind = DBL_PREC), intent(in) :: SED_K_MIN_DOC_S_PLUS_6_20
    real(kind = DBL_PREC), intent(in) :: SED_K_MIN_DOC_DOC_20
    real(kind = DBL_PREC), intent(in) :: SED_THETA_K_MIN_DOC_DOXY
    real(kind = DBL_PREC), intent(in) :: SED_THETA_K_MIN_DOC_NO3N
    real(kind = DBL_PREC), intent(in) :: SED_THETA_K_MIN_DOC_MN_IV
    real(kind = DBL_PREC), intent(in) :: SED_THETA_K_MIN_DOC_FE_III
    real(kind = DBL_PREC), intent(in) :: SED_THETA_K_MIN_DOC_S_PLUS_6
    real(kind = DBL_PREC), intent(in) :: SED_THETA_K_MIN_DOC_DOC
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_DOC_MIN_DOXY
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_DOC_MIN_NO3N
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_DOC_MIN_MN_IV
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_DOC_MIN_FE_III
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_DOC_MIN_DOC
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_DOXY_RED_LIM
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_NO3N_RED_LIM
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_MN_IV_RED_LIM
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_FE_III_RED_LIM
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_S_PLUS_6_RED_LIM
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_DOXY_RED_INHB
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_NO3N_RED_INHB
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_MN_IV_RED_INHB
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_FE_III_RED_INHB
    real(kind = DBL_PREC), intent(in) :: SED_K_HS_S_PLUS_6_RED_INHB
    real(kind = DBL_PREC), intent(in) :: SED_PH_MIN_DOC_MIN_DOXY
    real(kind = DBL_PREC), intent(in) :: SED_PH_MIN_DOC_MIN_NO3N
    real(kind = DBL_PREC), intent(in) :: SED_PH_MIN_DOC_MIN_MN_IV
    real(kind = DBL_PREC), intent(in) :: SED_PH_MIN_DOC_MIN_FE_III
    real(kind = DBL_PREC), intent(in) :: SED_PH_MIN_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC), intent(in) :: SED_PH_MIN_DOC_MIN_DOC
    real(kind = DBL_PREC), intent(in) :: SED_PH_MAX_DOC_MIN_DOXY
    real(kind = DBL_PREC), intent(in) :: SED_PH_MAX_DOC_MIN_NO3N
    real(kind = DBL_PREC), intent(in) :: SED_PH_MAX_DOC_MIN_MN_IV
    real(kind = DBL_PREC), intent(in) :: SED_PH_MAX_DOC_MIN_FE_III
    real(kind = DBL_PREC), intent(in) :: SED_PH_MAX_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC), intent(in) :: SED_PH_MAX_DOC_MIN_DOC

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: SED_TEMPS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: SED_DOC
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: SED_DOXY
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: SED_NO3N
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: MN_IV_DISS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: FE_III_DISS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: PH

    integer :: nkn
    integer :: NUM_SED_LAYERS

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: PH_CORR_DOC_MIN_DOXY
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: PH_CORR_DOC_MIN_NO3N
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: PH_CORR_DOC_MIN_MN_IV
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: PH_CORR_DOC_MIN_FE_III
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: PH_CORR_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: PH_CORR_DOC_MIN_DOC
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: LIM_DOXY_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: LIM_NO3N_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: LIM_MN_IV_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: LIM_FE_III_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: LIM_S_PLUS_6_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: LIM_DOC_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: K_NO3_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: K_MN_IV_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: K_FE_III_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: K_S_PLUS_6_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: K_DOC_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: R_MINER_DOC_DOXY
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: R_MINER_DOC_NO3N
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: R_MINER_DOC_MN_IV
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: R_MINER_DOC_FE_III
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: R_MINER_DOC_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: R_MINER_DOC_DOC

    !call CALCULATE_PH_CORR_SED(PH_CORR_DOC_MIN, PH, SED_PH_MIN_DOC_MIN, &
    !                       SED_PH_MIN_DOC_MAX, nkn, NUM_SED_LAYERS)

    call CALCULATE_PH_CORR_SED &
             (PH_CORR_DOC_MIN_DOXY     , PH , SED_PH_MIN_DOC_MIN_DOXY  , &
              SED_PH_MAX_DOC_MIN_DOXY  , nkn, NUM_SED_LAYERS)
    
    call CALCULATE_PH_CORR_SED &
             (PH_CORR_DOC_MIN_NO3N     , PH , SED_PH_MIN_DOC_MIN_NO3N  , &
              SED_PH_MAX_DOC_MIN_NO3N  , nkn, NUM_SED_LAYERS)
    
    call CALCULATE_PH_CORR_SED &
             (PH_CORR_DOC_MIN_MN_IV    , PH , SED_PH_MIN_DOC_MIN_MN_IV , &
              SED_PH_MAX_DOC_MIN_MN_IV , nkn, NUM_SED_LAYERS)
    
    call CALCULATE_PH_CORR_SED &
             (PH_CORR_DOC_MIN_FE_III   , PH , SED_PH_MIN_DOC_MIN_FE_III, &
              SED_PH_MAX_DOC_MIN_FE_III, nkn, NUM_SED_LAYERS)
    
    call CALCULATE_PH_CORR_SED(PH_CORR_DOC_MIN_S_PLUS_6, PH, SED_PH_MIN_DOC_MIN_S_PLUS_6, SED_PH_MAX_DOC_MIN_S_PLUS_6, nkn, &
                               NUM_SED_LAYERS)
    
    call CALCULATE_PH_CORR_SED(PH_CORR_DOC_MIN_DOC     , PH, SED_PH_MIN_DOC_MIN_DOC     , SED_PH_MAX_DOC_MIN_DOC     , nkn, &
                               NUM_SED_LAYERS)
        
    LIM_DOXY_RED = SED_DOXY  / (SED_DOXY + SED_K_HS_DOXY_RED_LIM)

    LIM_NO3N_RED = (SED_NO3N / (SED_NO3N + SED_K_HS_NO3N_RED_LIM)) * &
        (SED_K_HS_DOXY_RED_INHB / (SED_DOXY + SED_K_HS_DOXY_RED_INHB))

    LIM_MN_IV_RED = (MN_IV_DISS  / (MN_IV_DISS + SED_K_HS_MN_IV_RED_LIM)) * &
        (SED_K_HS_DOXY_RED_INHB / (SED_DOXY + SED_K_HS_DOXY_RED_INHB))            * &
        (SED_K_HS_NO3N_RED_INHB / (SED_NO3N + SED_K_HS_NO3N_RED_INHB))

    LIM_FE_III_RED = (FE_III_DISS  / (FE_III_DISS + SED_K_HS_FE_III_RED_LIM)) * &
        (SED_K_HS_DOXY_RED_INHB  / (SED_DOXY   + SED_K_HS_DOXY_RED_INHB))         * &
        (SED_K_HS_NO3N_RED_INHB  / (SED_NO3N   + SED_K_HS_NO3N_RED_INHB))         * &
        (SED_K_HS_MN_IV_RED_INHB / (MN_IV_DISS + SED_K_HS_MN_IV_RED_INHB))

    LIM_S_PLUS_6_RED = (S_PLUS_6 / (S_PLUS_6 + SED_K_HS_S_PLUS_6_RED_LIM)) * &
        (SED_K_HS_DOXY_RED_INHB   / (SED_DOXY    + SED_K_HS_DOXY_RED_INHB))    * &
        (SED_K_HS_NO3N_RED_INHB   / (SED_NO3N    + SED_K_HS_NO3N_RED_INHB))    * &
        (SED_K_HS_MN_IV_RED_INHB  / (MN_IV_DISS  + SED_K_HS_MN_IV_RED_INHB))   * &
        (SED_K_HS_FE_III_RED_INHB / (FE_III_DISS + SED_K_HS_FE_III_RED_INHB))

    LIM_DOC_RED = 1.0D0 - &
        (LIM_DOXY_RED + LIM_NO3N_RED + LIM_MN_IV_RED + LIM_FE_III_RED + LIM_S_PLUS_6_RED)

    where (LIM_DOC_RED < 0.0D0)
        LIM_DOC_RED = 0.0D0
    end where

    R_MINER_DOC_DOXY = &
        SED_K_MIN_DOC_DOXY_20 * (SED_THETA_K_MIN_DOC_DOXY ** (SED_TEMPS - 2.0D1)) * &
        LIM_DOXY_RED * PH_CORR_DOC_MIN_DOXY * (SED_DOC / (SED_DOC + SED_K_HS_DOC_MIN_DOXY)) * &
        SED_DOC

    R_MINER_DOC_NO3N = &
        SED_K_MIN_DOC_NO3N_20 * (SED_THETA_K_MIN_DOC_NO3N ** (SED_TEMPS - 2.0D1)) * &
        LIM_NO3N_RED * PH_CORR_DOC_MIN_NO3N * (SED_DOC / (SED_DOC + SED_K_HS_DOC_MIN_NO3N)) * &
        SED_DOC

    R_MINER_DOC_MN_IV = &
        SED_K_MIN_DOC_MN_IV_20 * (SED_THETA_K_MIN_DOC_MN_IV ** (SED_TEMPS - 2.0D1)) * &
        LIM_MN_IV_RED * PH_CORR_DOC_MIN_MN_IV * (SED_DOC / (SED_DOC + SED_K_HS_DOC_MIN_MN_IV)) * &
        SED_DOC

    R_MINER_DOC_FE_III = &
        SED_K_MIN_DOC_FE_III_20 * (SED_THETA_K_MIN_DOC_FE_III ** (SED_TEMPS - 2.0D1)) * &
        LIM_FE_III_RED * PH_CORR_DOC_MIN_FE_III * (SED_DOC / (SED_DOC + SED_K_HS_DOC_MIN_FE_III)) * &
        SED_DOC

    R_MINER_DOC_S_PLUS_6 = &
        SED_K_MIN_DOC_S_PLUS_6_20 * (SED_THETA_K_MIN_DOC_S_PLUS_6 ** (SED_TEMPS - 2.0D1)) * &
        LIM_S_PLUS_6_RED * PH_CORR_DOC_MIN_S_PLUS_6 * (SED_DOC / (SED_DOC + SED_K_HS_DOC_MIN_S_PLUS_6)) * &
        SED_DOC

    R_MINER_DOC_DOC = &
        (SED_K_MIN_DOC_DOC_20 * (SED_THETA_K_MIN_DOC_DOC ** (SED_TEMPS - 2.0D1)) * &
         LIM_DOC_RED * PH_CORR_DOC_MIN_DOXY * (SED_DOC / (SED_DOC + SED_K_HS_DOC_MIN_DOC)) * SED_DOC)

end subroutine SED_DOC_MINERALIZATION
    

subroutine SED_REDOX_AND_SPECIATION &
           (DOXY, NO3N, MN_IV, FE_III, S_PLUS_6, DISS_ORG_C, &
            S_MINUS_2 , MN_II, FE_II , HCO3    , CO3, &
            TEMP, SALT, PH, ELEVATION, &
            K_HS_DOXY_RED_LIM   , K_HS_NO3N_RED_LIM , K_HS_MN_IV_RED_LIM , &
            K_HS_FE_III_RED_LIM , K_HS_S_PLUS_6_RED_LIM, &
            K_HS_DOXY_RED_INHB  , K_HS_NO3N_RED_INHB, K_HS_MN_IV_RED_INHB, &
            K_HS_FE_III_RED_INHB, K_HS_S_PLUS_6_RED_INHB, nkn, NUM_SED_LAYERS, &
            LIM_DOXY_RED        , LIM_NO3N_RED          , LIM_MN_IV_RED  , &
            LIM_FE_III_RED      , LIM_S_PLUS_6_RED      , LIM_DOC_RED, &
            PE, FE_II_DISS, FE_III_DISS, MN_II_DISS)

    use AQUABC_II_GLOBAL
    implicit none

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: DOXY          ! Dissolved oxygen (mg/L)
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: NO3N          ! Nitrate nitrogen (mg/L)
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: MN_IV         ! Mn IV            (mg/L)
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: FE_III        ! Fe III           (mg/L)
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: S_PLUS_6      ! S +VI            (mg/L)
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: DISS_ORG_C    ! DOC              (mg/L)
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: S_MINUS_2     ! S -II            (mg/L)
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: MN_II         ! Mn II            (mg/L)
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: FE_II         ! Fe II            (mg/L)
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: HCO3          ! Bicarbonates     (moles)
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: CO3           ! Carbonate        (moles)

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: TEMP
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: SALT
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: PH
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(in) :: ELEVATION

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

    integer, intent(in) :: nkn, NUM_SED_LAYERS

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: LIM_DOXY_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: LIM_NO3N_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: LIM_MN_IV_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: LIM_FE_III_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: LIM_S_PLUS_6_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: LIM_DOC_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: PE
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: FE_II_DISS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: FE_III_DISS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: MN_II_DISS

    integer, dimension(nkn, NUM_SED_LAYERS) :: REDUCED_AGENT_NO

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS, 6) :: REDUCER_LIM_FACTORS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: CS, H_PLUS, SO4_MOLAR, HS_MOLAR, S_MINUS_2_MOLAR

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_CO3_OVER_FE_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_OH_2_OVER_FE_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_S_OVER_FE_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_S_2_OVER_FE_II

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_3_O_4_OVER_FE_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_2_O_3_OVER_FE_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_OOH_OVER_FE_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_OH_3_OVER_FE_II

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MN_CO3_OVER_MN_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MN_OH_2_OVER_MN_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MN_S_OVER_MN_II

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS, 4) :: FE_II_ACTIVITY_RATIOS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS, 3) :: FE_III_ACTIVITY_RATIOS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS, 3) :: MN_II_ACTIVITY_RATIOS

    integer, dimension(nkn, NUM_SED_LAYERS) :: FE_II_SALT_NO
    integer, dimension(nkn, NUM_SED_LAYERS) :: FE_III_SALT_NO
    integer, dimension(nkn, NUM_SED_LAYERS) :: MN_II_SALT_NO

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FREE_FE_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FREE_FE_III
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FREE_MN_II

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: TEMP_K
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: LOG_KW
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: KW
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: OH_MINUS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_III_FREE
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: H_PLUS_OVER_2
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: H_PLUS_OVER_3
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: H_PLUS_OVER_4


    real(kind = DBL_PREC) :: K1_FE_OH_3
    real(kind = DBL_PREC) :: BETA_2_FE_OH_3
    real(kind = DBL_PREC) :: BETA_3_FE_OH_3
    real(kind = DBL_PREC) :: BETA_4_FE_OH_3
    real(kind = DBL_PREC) :: BETA_2_2_FE_OH_3
    real(kind = DBL_PREC) :: BETA_4_3_FE_OH_3

    integer :: i, j

    !This is done to avoid exception. This a bug, variable does not get any value in this routine. fixme
    H_PLUS(:,:) = 10.0D0**(-PH)
    
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

    REDUCER_LIM_FACTORS(:,:,1) = LIM_DOXY_RED
    REDUCER_LIM_FACTORS(:,:,2) = LIM_NO3N_RED
    REDUCER_LIM_FACTORS(:,:,3) = LIM_MN_IV_RED
    REDUCER_LIM_FACTORS(:,:,4) = LIM_FE_III_RED
    REDUCER_LIM_FACTORS(:,:,5) = LIM_S_PLUS_6_RED
    REDUCER_LIM_FACTORS(:,:,6) = LIM_DOC_RED

    REDUCED_AGENT_NO = maxloc(REDUCER_LIM_FACTORS, dim = 3)
    
    call DO_SATURATION_MAT(TEMP, SALT, ELEVATION, nkn, NUM_SED_LAYERS, CS)

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

        PE = 4.25D0 - log((HS_MOLAR**0.125D0) / (((S_PLUS_6 / 32000.0D0) **0.125D0) * (H_PLUS**1.125D0)))
    end where

    ! Methanogenesis
    where (REDUCED_AGENT_NO == 6)
        PE = -0.2D0 - log10(1.0D0 / ((DISS_ORG_C / 12000D0)**0.25D0) * H_PLUS)
    end where

    ! Since PE is known, let's have a look in FE_II and FE_III solids following
    ! the Appendix 8.1, pp 513-515 in Stumm and Morgen (1996),
    ! "Aquatic Chemistry, Chemical Equilibria and Rates in Natural Waters"

    ! FE_CO3  / DISS_FE_II
    !write(*,*) '******************** HCO3 **************************** in sub'
    !write(*,*) HCO3
    !write(*,*) '******************** PH ****************************   in sub'
    !write(*,*) PH
    FE_CO3_OVER_FE_II  = 10.0D0**(-0.2D0 +  PH + log10(HCO3))

    ! FE_OH_2 / DISS_FE_II
    FE_OH_2_OVER_FE_II = 10.0D0**(-11.7D0 + (2.0D0*PH))

    ! FE_S / DISS_FE_II
    FE_S_OVER_FE_II    = &
        10.0D0**(38.0D0  - (8.0D0*PH) + log10((S_PLUS_6 / 32000.0D0)) - (8.0D0*PE))
    
    
    ! FE_S_2 / DISS_FE_II
    FE_S_2_OVER_FE_II  = &
        10.0D0**(86.8D0  - (16.0D0*PH) + (2.0D0*log10((S_PLUS_6 / 32000.0D0))) - &
                 (14.0D0*PE))

    ! Now find out which Fe II salt is more likely to form for each reactor
    FE_II_ACTIVITY_RATIOS (:, :, 1) = FE_CO3_OVER_FE_II
    FE_II_ACTIVITY_RATIOS (:, :, 2) = FE_OH_2_OVER_FE_II
    FE_II_ACTIVITY_RATIOS (:, :, 3) = FE_S_OVER_FE_II
    FE_II_ACTIVITY_RATIOS (:, :, 4) = FE_S_2_OVER_FE_II    
    
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
        ((8.9D-8 * 1.2D-13)  / ((H_PLUS * H_PLUS) + (H_PLUS * 8.9D-8) + &
         (8.9D-8 * 1.2D-13)))

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


    ! (FE_3_O_4^(1/3)) / DISS_FE_II
    !FE_3_O_4_OVER_FE_II = 10.0D0**(3.0D0*(-10.1+((2.66D0)*PH) + ((0.66D0)*PE)))

    ! -------------------------------------------------------------------------
    ! FE_III Species
    ! -------------------------------------------------------------------------

    
    !! (FE_2_O_3^(1/2)) / DISS_FE_II
    !FE_2_O_3_OVER_FE_II = 10.0D0**(2.0D0 * (-11.1D0 + (3.0D0*PH) + PE))
    !
    !! FE_OOH / DISS_FE_II
    !FE_OOH_OVER_FE_II   = 10.0D0**(-11.3D0 + (3.0D0*PH) + PE)
    !
    !! FE_OH_3 / DISS_FE_II
    !FE_OH_3_OVER_FE_II  = 10.0D0**(-17.1D0 + (3.0D0*PH) + PE)
    !
    !
    !! Now find out which Fe III salt is more likely to form for each reactor
    !FE_III_ACTIVITY_RATIOS(:, :, 1) = FE_2_O_3_OVER_FE_II
    !FE_III_ACTIVITY_RATIOS(:, :, 2) = FE_OOH_OVER_FE_II
    !FE_III_ACTIVITY_RATIOS(:, :, 3) = FE_OH_3_OVER_FE_II
    !
    !FE_III_SALT_NO = maxloc(FE_III_ACTIVITY_RATIOS, dim = 2)
    !
    !!Fe2O3
    !where(FE_III_SALT_NO == 1)
    !    FREE_FE_III = 0.0D0     !Better approximation some time
    !end where
    !
    !!FeOOH
    !where(FE_III_SALT_NO == 2)
    !    TEMP_K = TEMP + 273.15D0
    !    
    !    LOG_KW = &
    !        -2.839710D2 + (1.3323D4 / TEMP_K) - (5.069842D-2 * TEMP_K) + &
    !        (1.0224447D2 * log10(TEMP_K)) - (1.119669D6 / (TEMP_K * TEMP_K))
    !
    !    KW = 10**(LOG_KW)        
    !    OH_MINUS = H_PLUS / KW
    !    FREE_FE_III = (10.0D0**(42.97D0))/(27.0D0 * OH_MINUS * OH_MINUS * OH_MINUS)
    !end where
    !
    ! if(any(isnan(FREE_FE_III))) then
    !   print *,'REDOX_AND_SPECIATION:'
    !   print *,'1.FREE_FE_III is NaN:', FREE_FE_III
    !   stop
    ! end if
    !
    !!Fe(OH)3
    !where(FE_III_SALT_NO == 3)
    !    TEMP_K = TEMP + 273.15D0
    !    
    !    LOG_KW = &
    !        -2.839710D2 + (1.3323D4 / TEMP_K) - (5.069842D-2 * TEMP_K) + &
    !        (1.0224447D2 * log10(TEMP_K)) - (1.119669D6 / (TEMP_K * TEMP_K))
    !
    !    OH_MINUS = H_PLUS / KW
    !    FREE_FE_III = (10.0D0**(37.08D0))/(27.0D0 * OH_MINUS * OH_MINUS * OH_MINUS)
    !end where
    !
    ! if(any(isnan(FREE_FE_III))) then
    !   print *,'REDOX_AND_SPECIATION:'
    !   print *,'2.FREE_FE_III is NaN:', FREE_FE_III
    !   print *,'OH_MINUS:',OH_MINUS
    !   print *,'KW:',KW
    !   print *,'H_PLUS:',H_PLUS
    !   stop
    ! end if

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
        
        KW = 10.0D0**(LOG_KW)
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
    
    if(any(isnan(FE_III_DISS))) then
       print *,'FE_III_DISS is NaN,  SEDIMENT LIBRARY'
       
       write(unit = *, fmt = '(2A15, 13A20)') 'BOX', 'SED LAYER', 'PH' , &                       
             'FE_III_DISS', 'FE_III_FREE', 'K1_FE_OH_3', &
             'H_PLUS', 'H_PLUS_OVER_2', 'H_PLUS_OVER_3', 'H_PLUS_OVER_4', &
             'BETA_2_FE_OH_3', 'BETA_3_FE_OH_3', 'BETA_4_FE_OH_3', &
             'BETA_2_2_FE_OH_3', 'BETA_4_3_FE_OH_3'
       
       do i = 1, nkn
           do j = 1, NUM_SED_LAYERS
               write(unit = *, fmt = '(2i15, 13f20.15)') i, j, PH(i, j), &
                     FE_III_DISS   (i, j), FE_III_FREE     (i, j), K1_FE_OH_3            , &
                     H_PLUS        (i, j), H_PLUS_OVER_2   (i, j), H_PLUS_OVER_3   (i, j), &
                     H_PLUS_OVER_4 (i, j), BETA_2_FE_OH_3        , BETA_3_FE_OH_3        , &
                     BETA_4_FE_OH_3      , BETA_2_2_FE_OH_3      , BETA_4_3_FE_OH_3
           end do
       end do
       
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
    MN_II_ACTIVITY_RATIOS (:, :, 1) = MN_CO3_OVER_MN_II
    MN_II_ACTIVITY_RATIOS (:, :, 2) = MN_OH_2_OVER_MN_II
    MN_II_ACTIVITY_RATIOS (:, :, 3) = MN_S_OVER_MN_II

    MN_II_SALT_NO  = maxloc(MN_II_ACTIVITY_RATIOS , dim = 3)

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
end subroutine SED_REDOX_AND_SPECIATION
    
    
!Function, which returns saturation concentration of dissolved oxygen
subroutine DO_SATURATION_MAT(T, S, H, nkn, NUM_SED_LAYERS, CS)

    !Water temperature (in Celcius)
    double precision, dimension(nkn, NUM_SED_LAYERS), intent(in) :: T

    !Salinity (in ppt)
    double precision, dimension(nkn, NUM_SED_LAYERS), intent(in) :: S

    !Elevation (in m)
    double precision, dimension(nkn, NUM_SED_LAYERS), intent(in) :: H

    integer, intent(in) :: nkn, NUM_SED_LAYERS

    !Water temperature (in Kelvin)
    double precision, dimension(nkn, NUM_SED_LAYERS) :: T_KELVIN

    !Altitude (in feet)
    double precision, dimension(nkn, NUM_SED_LAYERS) :: H_FEET

    double precision, dimension(nkn, NUM_SED_LAYERS) :: LN_CSF
    double precision, dimension(nkn, NUM_SED_LAYERS) :: LN_CSS
    double precision, dimension(nkn, NUM_SED_LAYERS) :: CSS
    double precision, dimension(nkn, NUM_SED_LAYERS) :: CSP

    double precision, dimension(nkn, NUM_SED_LAYERS) :: CS
    !Pressure at altitude H (in atm)
    double precision, dimension(nkn, NUM_SED_LAYERS) :: P

    !Standart pressure (in mmHg)
    double precision, dimension(nkn, NUM_SED_LAYERS) :: P0

    double precision, dimension(nkn, NUM_SED_LAYERS) :: LN_PWV

    !Partial pressure of water vapor (in atm)
    double precision, dimension(nkn, NUM_SED_LAYERS) :: PWV

    !A constant
    double precision, dimension(nkn, NUM_SED_LAYERS) :: THETA

    T_KELVIN = T + 273.15
    H_FEET = H / 0.3048D0

    !Calculate the effect of temperature on dissolved oxygen saturation
    LN_CSF = -139.34411D0 + (157570.1d0 / T_KELVIN) - &
    &        (66423080.0D0     / (T_KELVIN * T_KELVIN)) + &
    &        (12438000000.0D0  / (T_KELVIN * T_KELVIN * T_KELVIN)) - &
    &        (862194900000.0D0 / (T_KELVIN * T_KELVIN * T_KELVIN * T_KELVIN))

    !Calculate the effect of salinity on dissolved oxygen saturation
    LN_CSS = LN_CSF - S * &
        (0.017674D0 - (10.754 / T_KELVIN) + (2140.7D0 / (T_KELVIN ** 2.0D0)))

    CSS = exp(LN_CSS)

    !Calculate the effect of altitude on dissolved oxygen saturation

    !Calculate THETA
    THETA = 0.000975 - (0.00001426 * T) + (0.00000006436 * (T ** 2.0D0))

    !Set standard pressure to mean sea level
    P0 = 760.0

    !Calculate atmospheric pressure at altitude H
    P = (P0 - (0.02667 * H_FEET)) / 760.0

    !Calculate vapour pressure of water(DIKKAT)
    LN_PWV = 11.8571 - (3840.7 / T_KELVIN) - (216961.0 / (T_KELVIN ** 2.0D0))

    PWV = exp(LN_PWV)

    !Final calculation including altitude effect
    CSP = CSS  * P * (((1.0D0 - (PWV / P)) * (1.0D0 - (THETA * P))) &
        / ((1 - PWV) * (1.0D0 - THETA)))

    CS = CSP
    end subroutine DO_SATURATION_MAT

! -----------------------------------------------------------------------------
! The subroutine below calculates the Iron II dissolution according to
! Stumm and Lee (1960).    
! -----------------------------------------------------------------------------
!
!                       Initial development 2nf of July 2016
!    
!                                by Ali Ertürk
! -----------------------------------------------------------------------------
subroutine SED_IRON_II_DISSOLUTION(HS2_TOT, PH, TOT_ALK, nkn, NUM_SED_LAYERS, FE_II_TOT)
    use AQUABC_II_GLOBAL
    implicit none

    integer, intent(in) ::  nkn, NUM_SED_LAYERS

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: HS2_TOT       ! [H2S] + [HS-]  +  [S]   (mol/L)
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: PH            ! PH
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: TOT_ALK       ! Total alkalinity (mol/L)
    
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS), intent(inout) :: FE_II_TOT     ! Total Fe2+ (mol/L)
    
    ! Equilibrium constant for the reaction : Fe(OH)2    <-------->    Fe2+    +     2(OH-) 
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_1      
    
    ! Equilibrium constant for the reaction : Fe(OH)2    <-------->    [Fe(OH)]+     +    OH- 
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_2      
    
    ! Equilibrium constant for the reaction : Fe(OH)2    +    OH-    <-------->     [Fe(OH)3]- 
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_3      
    
    ! Equilibrium constant for the reaction : Fe(CO3)    <-------->    Fe++    +     (CO3)--
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_4      
    
    ! Equilibrium constant for the reaction : Fe(OH)2    +    OH-    <-------->     [Fe(OH)]+    +    CO3-- 
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_5      
    
    ! Equilibrium constant for the reaction : HCO3-      <-------->    H+      +    CO3--         
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_6      
    
    ! Equilibrium constant for the reaction : FeS        <-------->    Fe++    +    S--
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_7      
    
    ! Equilibrium constant for the reaction : FeS        +    OH-    <-------->     [Fe(OH)]+    +    S--
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_8      
    
    ! Equilibrium constant for the reaction : FeS        +  3(OH-)   <-------->     [Fe(OH3)]-   +    S--
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_9      
    
    ! Equilibrium constant for the reaction : H2S        +    H+     <-------->     HS-
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_10_A   
    
    ! Equilibrium constant for the reaction : HS-        +    H+     <-------->     S--
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_10_B   
    
    ! Equilibrium constant for water dissociation
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_W      
    
    
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: H_PLUS               ! [H+]    (mol/L) 
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: OH_MINUS             ! [OH-]   (mol/L)
    
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS, 3) :: ALL_FE_II
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

    ALL_FE_II(:,:,1) = &
        ((K_1  / (K_W * K_W)) * H_PLUS * H_PLUS) + ((K_2 / K_W) * H_PLUS) + &
        ((K_3 * K_W) / H_PLUS)

    ALL_FE_II(:,:,2) = &
        ((H_PLUS + (2.0D0 * K_6)) / (TOT_ALK * K_6)) * (K_4 + ((K_5 * K_W) / H_PLUS))
    
    ALL_FE_II(:,:,3) = 1.0D0
    
    where (HS2_TOT < 1.0D-12)
        HS2_TOT  (:,:) = 1.0D-20
        FE_II_TOT(:,:) = minval(ALL_FE_II(:,:,1:2), dim=3)
    elsewhere        
        ALL_FE_II(:,:,3)  = &
            ((K_7 / HS2_TOT) * &
             (1.0D0 + (H_PLUS / K_10_B) + ((H_PLUS * H_PLUS) / (K_10_B * K_10_A)))) + &
             (K_7 * ((OH_MINUS / K_8) + ((27.0D0 * OH_MINUS * OH_MINUS * OH_MINUS) / K_9)))

      FE_II_TOT(:,:) = minval(ALL_FE_II, dim=3)
    end where

end subroutine SED_IRON_II_DISSOLUTION

    
subroutine FLX_SED_MOD_1_TO_ALUKAS_II_VEC &
           (FLUXES_FROM_SEDIMENT, NUM_FLUXES_FROM_SEDIMENT, &
            FLUXES_TO_ALUKAS    , nkn, NUM_FLUXES_TO_ALUKAS)
            
    ! Note: routine does not take into account particulate material resuspension
    ! as a flux yet. Everytthing what is resuspended is added to the same fluxes. fixme

    implicit none

    integer :: nkn
    integer :: NUM_FLUXES_FROM_SEDIMENT
    integer :: NUM_FLUXES_TO_ALUKAS

    double precision :: FLUXES_FROM_SEDIMENT(nkn, NUM_FLUXES_FROM_SEDIMENT)
    double precision :: FLUXES_TO_ALUKAS    (nkn, NUM_FLUXES_TO_ALUKAS)

    if(NUM_FLUXES_TO_ALUKAS .ne. 30) then
        print *,'FLX_SED_MOD_1_TO_ALUKAS_II:'
        print *,'To get values correctly by FLUXES_TO_ALUKAS'
        print *,'num. of statevars should be 30 but is', NUM_FLUXES_TO_ALUKAS
        stop
    end if

    !NH4 FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 1) = FLUXES_FROM_SEDIMENT(:, 1)

    !NO3 FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 2) = FLUXES_FROM_SEDIMENT(:, 2)

    !PO4 FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 3) = FLUXES_FROM_SEDIMENT(:, 5)

    !DISS_OXYGEN FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 4) = FLUXES_FROM_SEDIMENT(:, 8)

    !DIA_C FLUX TO ALUKAS
    !ZOO_C FLUX TO ALUKAS
    !ZOO_N FLUX TO ALUKAS
    !ZOO_P FLUX TO ALUKAS
    !DET_PART_ORG_C FLUX TO ALUKAS
    !DET_PART_ORG_N FLUX TO ALUKAS
    !DET_PART_ORG_P FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 5:11) = 0.0D0
    
    !DISS_ORG_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 12) = FLUXES_FROM_SEDIMENT(:, 9)

    !DISS_ORG_N FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 13) = FLUXES_FROM_SEDIMENT(:, 3)

    !DISS_ORG_P FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 14) = FLUXES_FROM_SEDIMENT(:, 6)

    !CYN_C FLUX TO ALUKAS
    !OPA_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 15:16) = 0.0D0

    !DISS_Si FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 17) = FLUXES_FROM_SEDIMENT(:, 11)

    !PART_Si FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 18) = FLUXES_FROM_SEDIMENT(:, 12)

    !FIX_CYN_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 19) = 0.0D0

    !INORG_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(:, 20) = FLUXES_FROM_SEDIMENT(:, 13)
    FLUXES_TO_ALUKAS(:, 21) = FLUXES_FROM_SEDIMENT(:, 14)

    !Metal fluxes to sediments
    FLUXES_TO_ALUKAS(:,22:30) = FLUXES_FROM_SEDIMENT(:,16:24)
end subroutine FLX_SED_MOD_1_TO_ALUKAS_II_VEC
    

    
!THIS IS A USER PROGRAMMED FUNCTION. IT IS THE USERS RESPONSIBILITY
!TO PROGRAM THE EQUATIONS WHICH GIVE THE MOLECULAR DIFFUSION COEFFICIENTS
!OF STATE VARIABLES IN WATER AS A FUNCTION OF TEMPERATURE AND SALINITY
subroutine SED_MOD_1_ALUKAS_MOLDI_C_VEC(T,  SAL, nkn, NUM_SED_LAYERS, NUM_SED_VARS, MOL_DIFF)
    implicit none
    !INPUTS:
    !SVARNO : State variable no
    !T      : Temperature in Celsius
    !SAL    : Salinity (ppt)
    !OTPUT:
    !MOL_DIFF : Molecular diffusion coefficient, units cm^2/sec inside
    !           and converted to m^2/s for the output

    integer nkn, NUM_SED_LAYERS, NUM_SED_VARS
    
    double precision T       (nkn, NUM_SED_LAYERS)
    double precision SAL     (nkn, NUM_SED_LAYERS)
    double precision MOL_DIFF(nkn, NUM_SED_LAYERS, NUM_SED_VARS)
    
    double precision TS   !Standard temperature
    double precision SS   !Standard salinity
    double precision P    !Pressure

    double precision V25(nkn, NUM_SED_LAYERS)
    double precision VTK(nkn, NUM_SED_LAYERS)
    double precision TK(nkn, NUM_SED_LAYERS)

    if(NUM_SED_VARS .ne. 24) then
        print *,'SED_MOD_1_ALUKAS_MOLDI_C:'
        print *,'To get values correctly by Molecular diffusion'
        print *,'BS statevars No should not gt than 24 but is', NUM_SED_VARS
        stop
    end if

    TK = T + 273.16
    TS = 25.0D0
    SS = 36.1D0
    P  =  1.0D0
    
    ! Shear viscosities
    !call SED_MOD_1_CVISC(V25, SS , TS, 1.0D0)
    !call SED_MOD_1_CVISC(VTK, SAL, T , 1.0D0)

    V25 = 1.7910 - TS * (6.144D-02 - TS*(1.4510D-03 - TS*1.6826D-05)) - &
          1.5290D-04 * P + 8.3885D-08 * P * P + 2.4727D-03 * SS + &
          (6.0574D-06*P - 2.6760D-09*P*P)*TS + (TS * (4.8429D-05 - &
          TS * (4.7172D-06 - TS * 7.5986D-08))) * SS    
    
    VTK = 1.7910 - T * (6.144D-02 - T*(1.4510D-03 - T*1.6826D-05)) - &
          1.5290D-04 * P + 8.3885D-08 * P * P + 2.4727D-03 * SAL + &
          (6.0574D-06*P - 2.6760D-09*P*P)*T + (T * (4.8429D-05 - &
          T * (4.7172D-06 - T * 7.5986D-08))) * SAL
    
    MOL_DIFF(:,:,1) = (9.5D0 + 4.13D-1 * T(:,:)) * 1.0D-6    !NH4N: Boudreau 1997, Springer-Verlag
    MOL_DIFF(:,:,2) = (9.5D0 + 3.88D-1 * T(:,:)) * 1.0D-6    !NO3N: Boudreau 1997, Springer-Verlag
    MOL_DIFF(:,:,3) = 1.0D-6                                 !Dissolved organic nitrogen (assumed value)
    MOL_DIFF(:,:,4) = 0.0D0                                  !Particulate organic nitrogen (not subjected to molecular diffusion)
    MOL_DIFF(:,:,5) = (2.62D0 + 1.43D-1 * T(:,:)) * 1.0D-6   !PO4P: Boudreau 1997, Springer-Verlag
    MOL_DIFF(:,:,6) = 1.0D-6                                 !Dissolved organic phosphorus (assumed value)
    MOL_DIFF(:,:,7) = 0.0D0                                  !Particulate organic phosphorus (not subjected to molecular diffusion)
    
    !DOXY: Fossing et al., 2004 (NERI Technical report)    
    MOL_DIFF(:,:,8) = (1.17D1 + (3.44D-1 * T(:,:)) + (5.05D-3 * (T(:,:) ** 2.0D0))) * 1.0D-6
    
    MOL_DIFF(:,:,9)  = 1.0D-6                                !Dissolved organic carbon (assumed value)
    MOL_DIFF(:,:,10) = 0.0D0                                 !Particulate organic carbon (not subjected to molecular diffusion)
    
    !DSi
    !From Boudreau 1997, Springer-Verlag
    !Wollast and Garrels (1971) found D(H4SiO4) at 25 deg C
    !and 36.1 ppt S., Assume that this value can be scaled by
    !the Stokes-Einstein relationship to any other temperature.    
    MOL_DIFF(:,:,11) = 1.0D-5 * (V25(:,:) / 298.16D0) * (TK(:,:) / VTK(:,:))
    
    !Particulate organic silicon (not subjected to molecular diffusion)
    MOL_DIFF(:,:,12) = 0.0D0
    
    !Dissolved inorganic carbon, alkalinity, salinity (assumed value 1.0D-6 to be corrected with better formulation)
    MOL_DIFF(:,:,13:15) = 1.0D-6                                
    
    MOL_DIFF(:,:,16) = (3.31D0 + 0.1500D0 * T(:,:)) * 1.0D-5     !FEII
    MOL_DIFF(:,:,17) = 1.0D-6                                    !FEIII
    MOL_DIFF(:,:,18) = (3.18D0 + 0.1553D0 * T(:,:)) * 1.0D-5     !MNII
    MOL_DIFF(:,:,19) = 1.0D-6                                    !MN_IV

    MOL_DIFF(:,:,20) = (3.60D0 + 0.179D0 * T(:,:)) * 1.0D-5    !Ca
    MOL_DIFF(:,:,21) = (3.43D0 + 0.144D0 * T(:,:)) * 1.0D-5    !Mg
    MOL_DIFF(:,:,22) = 4.72D-9 * TK(:,:)/(35.2D0**0.6) * (V25(:,:)/VTK(:,:))    !S_PLUS_6
    MOL_DIFF(:,:,23) = (4.88D0 + 0.232D0 * T(:,:)) * 1.0D-5    !S_MINUS_2
 
    !CH4_C
    MOL_DIFF(:,:,24) = 5.7524D-3 * exp(-(3300D0 + 8.94104D-4*(TK(:,:)-228D0)**1.5D0)/(1.9858775D0*TK(:,:)))

    ! Converting to m^2/sec
    MOL_DIFF = MOL_DIFF * 1.0D-4

end subroutine SED_MOD_1_ALUKAS_MOLDI_C_VEC

!***********************************************************************
!***********************************************************************   

