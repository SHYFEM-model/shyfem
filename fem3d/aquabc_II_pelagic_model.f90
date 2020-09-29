
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

! Pelagic kinetic model ALUKAS_II
! Version with variables calculated in subroutines
! Version with dissolved inorganic carbon and alkalinity as
! state variables.

!Contains:
! module aquabc_II_wc_ini
!     contains : subroutine calc_frac_avail_DON

! module PELAGIC_MODEL_CONSTANTS
! subroutine INIT_PELAGIC_MODEL_CONSTANTS
! subroutine PELAGIC_KINETICS
! SUBROUTINE cur_smith
! subroutine LIM_LIGHT
! function GROWTH_AT_TEMP
! subroutine AMMONIA_PREF
! function DO_SATURATION
! function KAWIND
! function AMMONIA_VOLATILIZATION
! subroutine CALCULATE_PH_CORR
! subroutine FLX_ALUKAS_II_TO_SED_MOD_1 
! function STRANGERSD


!==================================================================
 module aquabc_II_wc_ini
!==================================================================
!
! Initilizes some variables necessary for WC calculations

    implicit none
    double precision, save :: frac_avail_DON
   
!==================================================================
 contains
!==================================================================
    subroutine calc_frac_avail_DON
        ! Calculates fraction of available DON for cyanobacteria
        ! it is called in pelagic model every step
 
        frac_avail_DON = 0.25  ! sedt
     end subroutine calc_frac_avail_DON

!==================================================================
end module aquabc_II_wc_ini
!==================================



module PELAGIC_MODEL_CONSTANTS
    use AQUABC_II_GLOBAL

    !Model constants
    
    !prescribed aeration coefficient, if negative should be calculated
    real(kind = DBL_PREC) :: K_A             
    real(kind = DBL_PREC) :: THETA_K_A
    real(kind = DBL_PREC) :: XKC
    real(kind = DBL_PREC) :: PHIMX    
    real(kind = DBL_PREC) :: KG_DIA_OPT_TEMP
    real(kind = DBL_PREC) :: KAPPA_DIA_UNDER_OPT_TEMP
    real(kind = DBL_PREC) :: KAPPA_DIA_OVER_OPT_TEMP
    real(kind = DBL_PREC) :: KR_DIA_20
    real(kind = DBL_PREC) :: THETA_KR_DIA
    real(kind = DBL_PREC) :: KD_DIA_20
    real(kind = DBL_PREC) :: THETA_KD_DIA
    real(kind = DBL_PREC) :: KHS_DIN_DIA
    real(kind = DBL_PREC) :: KHS_DIP_DIA
    real(kind = DBL_PREC) :: KHS_O2_DIA
    real(kind = DBL_PREC) :: KHS_NH4N_PREF_DIA
    real(kind = DBL_PREC) :: I_S_DIA
    real(kind = DBL_PREC) :: DO_STR_HYPOX_DIA_D
    real(kind = DBL_PREC) :: THETA_HYPOX_DIA_D
    real(kind = DBL_PREC) :: EXPON_HYPOX_DIA_D
    real(kind = DBL_PREC) :: DIA_N_TO_C
    real(kind = DBL_PREC) :: DIA_P_TO_C
    real(kind = DBL_PREC) :: DIA_O2_TO_C
    real(kind = DBL_PREC) :: DIA_C_TO_CHLA
    real(kind = DBL_PREC) :: DIA_C_TO_CHLA_NEW

    real(kind = DBL_PREC) :: KG_CYN_OPT_TEMP
    real(kind = DBL_PREC) :: KAPPA_CYN_UNDER_OPT_TEMP
    real(kind = DBL_PREC) :: KAPPA_CYN_OVER_OPT_TEMP
    real(kind = DBL_PREC) :: KR_CYN_20
    real(kind = DBL_PREC) :: THETA_KR_CYN
    real(kind = DBL_PREC) :: KD_CYN_20
    real(kind = DBL_PREC) :: THETA_KD_CYN
    real(kind = DBL_PREC) :: KHS_DIN_CYN
    real(kind = DBL_PREC) :: KHS_DIP_CYN
    real(kind = DBL_PREC) :: KHS_O2_CYN
    real(kind = DBL_PREC) :: KHS_NH4N_PREF_CYN
    real(kind = DBL_PREC) :: I_S_CYN
    real(kind = DBL_PREC) :: DO_STR_HYPOX_CYN_D
    real(kind = DBL_PREC) :: THETA_HYPOX_CYN_D
    real(kind = DBL_PREC) :: EXPON_HYPOX_CYN_D
    real(kind = DBL_PREC) :: CYN_N_TO_C
    real(kind = DBL_PREC) :: CYN_P_TO_C
    real(kind = DBL_PREC) :: CYN_O2_TO_C
    real(kind = DBL_PREC) :: CYN_C_TO_CHLA
    real(kind = DBL_PREC) :: CYN_C_TO_CHLA_NEW

    real(kind = DBL_PREC) :: KG_OPA_OPT_TEMP
    real(kind = DBL_PREC) :: KAPPA_OPA_UNDER_OPT_TEMP
    real(kind = DBL_PREC) :: KAPPA_OPA_OVER_OPT_TEMP
    real(kind = DBL_PREC) :: KR_OPA_20
    real(kind = DBL_PREC) :: THETA_KR_OPA
    real(kind = DBL_PREC) :: KD_OPA_20
    real(kind = DBL_PREC) :: THETA_KD_OPA
    real(kind = DBL_PREC) :: KHS_DIN_OPA
    real(kind = DBL_PREC) :: KHS_DIP_OPA
    real(kind = DBL_PREC) :: KHS_O2_OPA
    real(kind = DBL_PREC) :: KHS_NH4N_PREF_OPA
    real(kind = DBL_PREC) :: I_S_OPA
    real(kind = DBL_PREC) :: DO_STR_HYPOX_OPA_D
    real(kind = DBL_PREC) :: THETA_HYPOX_OPA_D
    real(kind = DBL_PREC) :: EXPON_HYPOX_OPA_D
    real(kind = DBL_PREC) :: OPA_N_TO_C
    real(kind = DBL_PREC) :: OPA_P_TO_C
    real(kind = DBL_PREC) :: OPA_O2_TO_C
    real(kind = DBL_PREC) :: OPA_C_TO_CHLA
    real(kind = DBL_PREC) :: OPA_C_TO_CHLA_NEW

    real(kind = DBL_PREC) :: KG_ZOO_OPT_TEMP
    real(kind = DBL_PREC) :: KAPPA_ZOO_UNDER_OPT_TEMP
    real(kind = DBL_PREC) :: KAPPA_ZOO_OVER_OPT_TEMP
    real(kind = DBL_PREC) :: GRAT_ZOO_DIA
    real(kind = DBL_PREC) :: GRAT_ZOO_CYN
    real(kind = DBL_PREC) :: GRAT_ZOO_OPA
    real(kind = DBL_PREC) :: GRAT_ZOO_FIX_CYN
    real(kind = DBL_PREC) :: GRAT_ZOO_DET_PART_ORG_C
    real(kind = DBL_PREC) :: PREF_ZOO_DIA
    real(kind = DBL_PREC) :: PREF_ZOO_CYN
    real(kind = DBL_PREC) :: PREF_ZOO_FIX_CYN
    real(kind = DBL_PREC) :: PREF_ZOO_OPA
    real(kind = DBL_PREC) :: PREF_ZOO_DET_PART_ORG_C
    real(kind = DBL_PREC) :: KHS_DIA_C_ZOO
    real(kind = DBL_PREC) :: KHS_CYN_C_ZOO
    real(kind = DBL_PREC) :: KHS_OPA_C_ZOO
    real(kind = DBL_PREC) :: KHS_FIX_CYN_C_ZOO
    real(kind = DBL_PREC) :: KHS_DET_PART_ORG_C_ZOO

    real(kind = DBL_PREC) :: KHS_MIN_P
    real(kind = DBL_PREC) :: KHS_MIN_N
    real(kind = DBL_PREC) :: KHS_NITR_NH4_N
    real(kind = DBL_PREC) :: KHS_NITR_OXY

    real(kind = DBL_PREC) :: KHS_AMIN_P
    real(kind = DBL_PREC) :: KHS_AMIN_N
    real(kind = DBL_PREC) :: KHS_DISS_N
    real(kind = DBL_PREC) :: KHS_DISS_P


    real(kind = DBL_PREC) :: FOOD_MIN_ZOO
    real(kind = DBL_PREC) :: THETA_KE_ZOO
    real(kind = DBL_PREC) :: KR_ZOO_20
    real(kind = DBL_PREC) :: THETA_KR_ZOO
    real(kind = DBL_PREC) :: KD_ZOO_20
    real(kind = DBL_PREC) :: THETA_KD_ZOO
    real(kind = DBL_PREC) :: DO_STR_HYPOX_ZOO_D
    real(kind = DBL_PREC) :: THETA_HYPOX_ZOO_D
    real(kind = DBL_PREC) :: EXPON_HYPOX_ZOO_D
    real(kind = DBL_PREC) :: ZOO_N_TO_C
    real(kind = DBL_PREC) :: ZOO_P_TO_C
    real(kind = DBL_PREC) :: ZOO_O2_TO_C
    real(kind = DBL_PREC) :: KDISS_DET_PART_ORG_C_20
    real(kind = DBL_PREC) :: THETA_KDISS_DET_PART_ORG_C
    real(kind = DBL_PREC) :: KDISS_DET_PART_ORG_N_20
    real(kind = DBL_PREC) :: THETA_KDISS_DET_PART_ORG_N
    real(kind = DBL_PREC) :: KDISS_DET_PART_ORG_P_20
    real(kind = DBL_PREC) :: THETA_KDISS_DET_PART_ORG_P

    real(kind = DBL_PREC) :: KHS_DSI_DIA
    real(kind = DBL_PREC) :: DIA_SI_TO_C
    real(kind = DBL_PREC) :: KDISS_PART_SI_20
    real(kind = DBL_PREC) :: THETA_KDISS_PART_SI

    real(kind = DBL_PREC) :: DIA_OPT_TEMP_LR
    real(kind = DBL_PREC) :: DIA_OPT_TEMP_UR
    real(kind = DBL_PREC) :: CYN_OPT_TEMP_LR
    real(kind = DBL_PREC) :: CYN_OPT_TEMP_UR
    real(kind = DBL_PREC) :: FIX_CYN_OPT_TEMP_LR
    real(kind = DBL_PREC) :: FIX_CYN_OPT_TEMP_UR
    real(kind = DBL_PREC) :: OPA_OPT_TEMP_LR
    real(kind = DBL_PREC) :: OPA_OPT_TEMP_UR
    real(kind = DBL_PREC) :: ZOO_OPT_TEMP_LR
    real(kind = DBL_PREC) :: ZOO_OPT_TEMP_UR

    real(kind = DBL_PREC) :: KE_ZOO
    real(kind = DBL_PREC) :: FRAC_ZOO_EX_ORG

    !Abitotic mineralization and nitrification (bacteria not modelled)
    real(kind = DBL_PREC) :: K_NITR_20
    real(kind = DBL_PREC) :: THETA_K_NITR

    real(kind = DBL_PREC) :: EFF_DIA_GROWTH
    real(kind = DBL_PREC) :: EFF_CYN_GROWTH
    real(kind = DBL_PREC) :: EFF_OPA_GROWTH
    real(kind = DBL_PREC) :: EFF_ZOO_GROWTH

    !Nitrogen fixing cyanobacteria
    real(kind = DBL_PREC) :: KG_FIX_CYN_OPT_TEMP
    real(kind = DBL_PREC) :: EFF_FIX_CYN_GROWTH
    real(kind = DBL_PREC) :: KAPPA_FIX_CYN_UNDER_OPT_TEMP
    real(kind = DBL_PREC) :: KAPPA_FIX_CYN_OVER_OPT_TEMP
    real(kind = DBL_PREC) :: KR_FIX_CYN_20
    real(kind = DBL_PREC) :: THETA_KR_FIX_CYN
    real(kind = DBL_PREC) :: KD_FIX_CYN_20
    real(kind = DBL_PREC) :: THETA_KD_FIX_CYN
    real(kind = DBL_PREC) :: KHS_DIN_FIX_CYN
    real(kind = DBL_PREC) :: KHS_DIP_FIX_CYN
    real(kind = DBL_PREC) :: KHS_O2_FIX_CYN
    real(kind = DBL_PREC) :: KHS_NH4N_PREF_FIX_CYN
    real(kind = DBL_PREC) :: I_S_FIX_CYN
    real(kind = DBL_PREC) :: DO_STR_HYPOX_FIX_CYN_D
    real(kind = DBL_PREC) :: THETA_HYPOX_FIX_CYN_D
    real(kind = DBL_PREC) :: EXPON_HYPOX_FIX_CYN_D
    real(kind = DBL_PREC) :: FIX_CYN_N_TO_C
    real(kind = DBL_PREC) :: FIX_CYN_P_TO_C
    real(kind = DBL_PREC) :: FIX_CYN_O2_TO_C
    real(kind = DBL_PREC) :: FIX_CYN_C_TO_CHLA
    real(kind = DBL_PREC) :: FIX_CYN_C_TO_CHLA_NEW
    real(kind = DBL_PREC) :: R_FIX
    real(kind = DBL_PREC) :: K_FIX
    real(kind = DBL_PREC) :: FRAC_FIX_CYN_EXCR

    real(kind = DBL_PREC) :: FRAC_CYN_EXCR
    real(kind = DBL_PREC) :: FRAC_OPA_EXCR
    real(kind = DBL_PREC) :: FRAC_DIA_EXCR
    real(kind = DBL_PREC) :: FAC_PHYT_AMIN_DON
    real(kind = DBL_PREC) :: FAC_PHYT_AMIN_DOC
    real(kind = DBL_PREC) :: FAC_PHYT_AMIN_DOP
    real(kind = DBL_PREC) :: FAC_PHYT_DET_PART_ORG_C
    real(kind = DBL_PREC) :: FAC_PHYT_DET_PART_ORG_P
    real(kind = DBL_PREC) :: FAC_PHYT_DET_PART_ORG_N

    ! New model constants added 9 September 2015
    real(kind = DBL_PREC) :: k_OX_FE_II
    real(kind = DBL_PREC) :: k_RED_FE_III
    real(kind = DBL_PREC) :: k_OX_MN_II
    real(kind = DBL_PREC) :: k_RED_MN_IV
    real(kind = DBL_PREC) :: KHS_DOXY_FE_III_RED
    real(kind = DBL_PREC) :: KHS_DOXY_MN_IV_RED
    ! End of new model constants added 22 September 2015

    ! New model constants introduced 27 January 2016 for the redox sequences
    real(kind = DBL_PREC) :: K_MIN_DOC_DOXY_20
    real(kind = DBL_PREC) :: K_MIN_DOC_NO3N_20
    real(kind = DBL_PREC) :: K_MIN_DOC_MN_IV_20
    real(kind = DBL_PREC) :: K_MIN_DOC_FE_III_20
    real(kind = DBL_PREC) :: K_MIN_DOC_S_PLUS_6_20
    real(kind = DBL_PREC) :: K_MIN_DOC_DOC_20
    real(kind = DBL_PREC) :: THETA_K_MIN_DOC_DOXY
    real(kind = DBL_PREC) :: THETA_K_MIN_DOC_NO3N
    real(kind = DBL_PREC) :: THETA_K_MIN_DOC_MN_IV
    real(kind = DBL_PREC) :: THETA_K_MIN_DOC_FE_III
    real(kind = DBL_PREC) :: THETA_K_MIN_DOC_S_PLUS_6
    real(kind = DBL_PREC) :: THETA_K_MIN_DOC_DOC
    real(kind = DBL_PREC) :: K_HS_DOC_MIN_DOXY
    real(kind = DBL_PREC) :: K_HS_DOC_MIN_NO3N
    real(kind = DBL_PREC) :: K_HS_DOC_MIN_MN_IV
    real(kind = DBL_PREC) :: K_HS_DOC_MIN_FE_III
    real(kind = DBL_PREC) :: K_HS_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: K_HS_DOC_MIN_DOC
    real(kind = DBL_PREC) :: K_HS_DOXY_RED_LIM
    real(kind = DBL_PREC) :: K_HS_NO3N_RED_LIM
    real(kind = DBL_PREC) :: K_HS_MN_IV_RED_LIM
    real(kind = DBL_PREC) :: K_HS_FE_III_RED_LIM
    real(kind = DBL_PREC) :: K_HS_S_PLUS_6_RED_LIM
    real(kind = DBL_PREC) :: K_HS_DOXY_RED_INHB
    real(kind = DBL_PREC) :: K_HS_NO3N_RED_INHB
    real(kind = DBL_PREC) :: K_HS_MN_IV_RED_INHB
    real(kind = DBL_PREC) :: K_HS_FE_III_RED_INHB
    real(kind = DBL_PREC) :: K_HS_S_PLUS_6_RED_INHB
    real(kind = DBL_PREC) :: PH_MIN_DOC_MIN_DOXY
    real(kind = DBL_PREC) :: PH_MIN_DOC_MIN_NO3N
    real(kind = DBL_PREC) :: PH_MIN_DOC_MIN_MN_IV
    real(kind = DBL_PREC) :: PH_MIN_DOC_MIN_FE_III
    real(kind = DBL_PREC) :: PH_MIN_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: PH_MIN_DOC_MIN_DOC
    real(kind = DBL_PREC) :: PH_MAX_DOC_MIN_DOXY
    real(kind = DBL_PREC) :: PH_MAX_DOC_MIN_NO3N
    real(kind = DBL_PREC) :: PH_MAX_DOC_MIN_MN_IV
    real(kind = DBL_PREC) :: PH_MAX_DOC_MIN_FE_III
    real(kind = DBL_PREC) :: PH_MAX_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: PH_MAX_DOC_MIN_DOC
    ! End of New model constants introduced 27 January 2016 for the redox sequences

    ! New model constants introduced 28 January 2016 for the redox sequences
    real(kind = DBL_PREC) :: K_MIN_DON_DOXY_20
    real(kind = DBL_PREC) :: K_MIN_DON_NO3N_20
    real(kind = DBL_PREC) :: K_MIN_DON_MN_IV_20
    real(kind = DBL_PREC) :: K_MIN_DON_FE_III_20
    real(kind = DBL_PREC) :: K_MIN_DON_S_PLUS_6_20
    real(kind = DBL_PREC) :: K_MIN_DON_DOC_20
    real(kind = DBL_PREC) :: THETA_K_MIN_DON_DOXY
    real(kind = DBL_PREC) :: THETA_K_MIN_DON_NO3N
    real(kind = DBL_PREC) :: THETA_K_MIN_DON_MN_IV
    real(kind = DBL_PREC) :: THETA_K_MIN_DON_FE_III
    real(kind = DBL_PREC) :: THETA_K_MIN_DON_S_PLUS_6
    real(kind = DBL_PREC) :: THETA_K_MIN_DON_DOC
    real(kind = DBL_PREC) :: K_HS_DON_MIN_DOXY
    real(kind = DBL_PREC) :: K_HS_DON_MIN_NO3N
    real(kind = DBL_PREC) :: K_HS_DON_MIN_MN_IV
    real(kind = DBL_PREC) :: K_HS_DON_MIN_FE_III
    real(kind = DBL_PREC) :: K_HS_DON_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: K_HS_DON_MIN_DOC
    real(kind = DBL_PREC) :: PH_MIN_DON_MIN_DOXY
    real(kind = DBL_PREC) :: PH_MIN_DON_MIN_NO3N
    real(kind = DBL_PREC) :: PH_MIN_DON_MIN_MN_IV
    real(kind = DBL_PREC) :: PH_MIN_DON_MIN_FE_III
    real(kind = DBL_PREC) :: PH_MIN_DON_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: PH_MIN_DON_MIN_DOC
    real(kind = DBL_PREC) :: PH_MAX_DON_MIN_DOXY
    real(kind = DBL_PREC) :: PH_MAX_DON_MIN_NO3N
    real(kind = DBL_PREC) :: PH_MAX_DON_MIN_MN_IV
    real(kind = DBL_PREC) :: PH_MAX_DON_MIN_FE_III
    real(kind = DBL_PREC) :: PH_MAX_DON_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: PH_MAX_DON_MIN_DOC

    real(kind = DBL_PREC) :: K_MIN_DOP_DOXY_20
    real(kind = DBL_PREC) :: K_MIN_DOP_NO3N_20
    real(kind = DBL_PREC) :: K_MIN_DOP_MN_IV_20
    real(kind = DBL_PREC) :: K_MIN_DOP_FE_III_20
    real(kind = DBL_PREC) :: K_MIN_DOP_S_PLUS_6_20
    real(kind = DBL_PREC) :: K_MIN_DOP_DOC_20
    real(kind = DBL_PREC) :: THETA_K_MIN_DOP_DOXY
    real(kind = DBL_PREC) :: THETA_K_MIN_DOP_NO3N
    real(kind = DBL_PREC) :: THETA_K_MIN_DOP_MN_IV
    real(kind = DBL_PREC) :: THETA_K_MIN_DOP_FE_III
    real(kind = DBL_PREC) :: THETA_K_MIN_DOP_S_PLUS_6
    real(kind = DBL_PREC) :: THETA_K_MIN_DOP_DOC
    real(kind = DBL_PREC) :: K_HS_DOP_MIN_DOXY
    real(kind = DBL_PREC) :: K_HS_DOP_MIN_NO3N
    real(kind = DBL_PREC) :: K_HS_DOP_MIN_MN_IV
    real(kind = DBL_PREC) :: K_HS_DOP_MIN_FE_III
    real(kind = DBL_PREC) :: K_HS_DOP_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: K_HS_DOP_MIN_DOC
    real(kind = DBL_PREC) :: PH_MIN_DOP_MIN_DOXY
    real(kind = DBL_PREC) :: PH_MIN_DOP_MIN_NO3N
    real(kind = DBL_PREC) :: PH_MIN_DOP_MIN_MN_IV
    real(kind = DBL_PREC) :: PH_MIN_DOP_MIN_FE_III
    real(kind = DBL_PREC) :: PH_MIN_DOP_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: PH_MIN_DOP_MIN_DOC
    real(kind = DBL_PREC) :: PH_MAX_DOP_MIN_DOXY
    real(kind = DBL_PREC) :: PH_MAX_DOP_MIN_NO3N
    real(kind = DBL_PREC) :: PH_MAX_DOP_MIN_MN_IV
    real(kind = DBL_PREC) :: PH_MAX_DOP_MIN_FE_III
    real(kind = DBL_PREC) :: PH_MAX_DOP_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: PH_MAX_DOP_MIN_DOC
    ! End of new model constants introduced 28 January 2016 for the redox sequences

    ! New model constats added in 29 th of January 2016
    real(kind = DBL_PREC) :: k_OX_CH4
    real(kind = DBL_PREC) :: THETA_k_OX_CH4
    real(kind = DBL_PREC) :: k_HS_OX_CH4_DOXY
    real(kind = DBL_PREC) :: k_OX_H2S
    real(kind = DBL_PREC) :: THETA_k_OX_H2S
    real(kind = DBL_PREC) :: k_HS_OX_H2S_DOXY
    ! End of new model constats added in 29 th of January 2016


    ! New model constants added in 9 August 2016
    real(kind = DBL_PREC) :: k_DISS_FE_II_20
    real(kind = DBL_PREC) :: THETA_k_DISS_FE_II
    real(kind = DBL_PREC) :: INIT_MULT_FE_II_DISS
    real(kind = DBL_PREC) :: k_DISS_FE_III_20
    real(kind = DBL_PREC) :: THETA_k_DISS_FE_III
    real(kind = DBL_PREC) :: INIT_MULT_FE_III_DISS
    ! End of new model constants added in 9 August 2016

    ! constants
    real(kind = DBL_PREC) :: PH_NITR_NH4_MIN
    real(kind = DBL_PREC) :: PH_NITR_NH4_MAX
    ! -------------------------------------------------------------------------
    !
    ! -------------------------------------------------------------------------

    !end of constatnts
end module PELAGIC_MODEL_CONSTANTS



subroutine INIT_PELAGIC_MODEL_CONSTANTS
    use PELAGIC_MODEL_CONSTANTS
    use para_aqua

    call para_get_value('K_A'                              ,                              K_A) !  1 Aeration coefficient (if negative calculates internally)
    call para_get_value('THETA_K_A'                        ,                        THETA_K_A) !  2 Temperature correction factor for aeration
    call para_get_value('XKC', XKC)
    call para_get_value('PHIMX', PHIMX)

    call para_get_value('KG_DIA_OPT_TEMP'                 ,                  KG_DIA_OPT_TEMP) !67  Diatoms Growth rate
    call para_get_value('DIA_OPT_TEMP_LR'                 ,                  DIA_OPT_TEMP_LR) !68  Diatoms optimal temperature lower range
    call para_get_value('DIA_OPT_TEMP_UR'                 ,                  DIA_OPT_TEMP_UR) !69  Diatoms optimal temperature upper range
    call para_get_value('EFF_DIA_GROWTH'                  ,                   EFF_DIA_GROWTH) !70  Diatoms Effective growth. (1-EG)*growth - losses for respiration and excretion
    call para_get_value('KAPPA_DIA_UNDER_OPT_TEMP'        ,         KAPPA_DIA_UNDER_OPT_TEMP) !71  Diatoms Temperature correction for growth lower temperature
    call para_get_value('KAPPA_DIA_OVER_OPT_TEMP'         ,          KAPPA_DIA_OVER_OPT_TEMP) !72  Diatoms Temperature correction for growth upper temperature
    call para_get_value('KR_DIA_20'                       ,                        KR_DIA_20) !73  Diatoms Respiration rate
    call para_get_value('THETA_KR_DIA'                    ,                     THETA_KR_DIA) !74  Diatoms Temperature correction for basal respiration rate
    call para_get_value('KD_DIA_20'                       ,                        KD_DIA_20) !75  Diatoms Mortality rate
    call para_get_value('THETA_KD_DIA'                    ,                     THETA_KD_DIA) !76  Diatoms Temperature correction for Mortality rate
    call para_get_value('KHS_DIN_DIA'                     ,                      KHS_DIN_DIA) !77  Diatoms Half saturation growth for DIN
    call para_get_value('KHS_DIP_DIA'                     ,                      KHS_DIP_DIA) !78  Diatoms Half saturation growth for DIP
    call para_get_value('KHS_DSi_DIA'                     ,                      KHS_DSi_DIA) !79  Diatoms Half saturation growth for DSi
    call para_get_value('KHS_O2_DIA'                      ,                       KHS_O2_DIA) !80  Diatoms Half saturation growth for O2
    call para_get_value('FRAC_DIA_EXCR'                   ,                    FRAC_DIA_EXCR) !81  Diatoms Fraction of excretion in metabolism rate
    call para_get_value('I_S_DIA'                         ,                          I_S_DIA) !82  Diatoms Light saturation (langleys)
    call para_get_value('DO_STR_HYPOX_DIA_D'              ,               DO_STR_HYPOX_DIA_D) !83  Diatoms Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
    call para_get_value('THETA_HYPOX_DIA_D'               ,                THETA_HYPOX_DIA_D) !84  Diatoms Multiplier of the exponent for Dissolved oxygen stress
    call para_get_value('EXPON_HYPOX_DIA_D'               ,                EXPON_HYPOX_DIA_D) !85  Diatoms Exponent constant for Dissolved oxygen stress
    call para_get_value('DIA_N_TO_C'                      ,                       DIA_N_TO_C) !86  Diatoms Nitrogen to Carbon ratio
    call para_get_value('DIA_P_TO_C'                      ,                       DIA_P_TO_C) !87  Diatoms Phosphorus to Carbon ratio
    call para_get_value('DIA_Si_TO_C'                     ,                      DIA_Si_TO_C) !88  Diatoms Silica to Carbon ratio
    call para_get_value('DIA_O2_TO_C'                     ,                      DIA_O2_TO_C) !89  Diatoms Oxygen to Carbon ratio for respiration
    call para_get_value('DIA_C_TO_CHLA'                   ,                    DIA_C_TO_CHLA) !90  Diatoms Carbon to Chlorophil a ratio
    call para_get_value('KG_CYN_OPT_TEMP'                 ,                  KG_CYN_OPT_TEMP) !91  Non-fixing cyanobacteria Growth rate
    call para_get_value('CYN_OPT_TEMP_LR'                 ,                  CYN_OPT_TEMP_LR) !92  Non-fixing cyanobacteria optimal temperature lower range
    call para_get_value('CYN_OPT_TEMP_UR'                 ,                  CYN_OPT_TEMP_UR) !93  Non-fixing cyanobacteria optimal temperature upper range
    call para_get_value('EFF_CYN_GROWTH'                  ,                   EFF_CYN_GROWTH) !94  Non-fixing cyanobacteria Effective growth. (1-EG)*growth - losses for respiration and excretion
    call para_get_value('KAPPA_CYN_UNDER_OPT_TEMP'        ,         KAPPA_CYN_UNDER_OPT_TEMP) !95  Non-fixing cyanobacteria Temperature correction for growth lower temperature
    call para_get_value('KAPPA_CYN_OVER_OPT_TEMP'         ,          KAPPA_CYN_OVER_OPT_TEMP) !96  Non-fixing cyanobacteria Temperature correction for growth upper temperature
    call para_get_value('KR_CYN_20'                       ,                        KR_CYN_20) !97  Non-fixing cyanobacteria Respiration rate
    call para_get_value('THETA_KR_CYN'                    ,                     THETA_KR_CYN) !98  Non-fixing cyanobacteria Temperature correction for respiration rate
    call para_get_value('KD_CYN_20'                       ,                        KD_CYN_20) !99  Non-fixing cyanobacteria Mortality rate
    call para_get_value('THETA_KD_CYN'                    ,                     THETA_KD_CYN) !100 Non-fixing cyanobacteria Temperature correction for Mortality rate
    call para_get_value('KHS_DIN_CYN'                     ,                      KHS_DIN_CYN) !101 Non-fixing cyanobacteria Half saturation growth for DIN
    call para_get_value('KHS_DIP_CYN'                     ,                      KHS_DIP_CYN) !102 Non-fixing cyanobacteria Half saturation growth for DIP
    call para_get_value('KHS_O2_CYN'                      ,                       KHS_O2_CYN) !103 Non-fixing cyanobacteria Half saturation growth for O2
    call para_get_value('FRAC_CYN_EXCR'                   ,                    FRAC_CYN_EXCR) !104 Non-fixing cyanobacteria Fraction of excretion in metabolism rate
    call para_get_value('I_S_CYN'                         ,                          I_S_CYN) !105 Non-fixing cyanobacteria Light saturation (langleys)
    call para_get_value('DO_STR_HYPOX_CYN_D'              ,               DO_STR_HYPOX_CYN_D) !106 Non-fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
    call para_get_value('THETA_HYPOX_CYN_D'               ,                THETA_HYPOX_CYN_D) !107 Non-fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
    call para_get_value('EXPON_HYPOX_CYN_D'               ,                EXPON_HYPOX_CYN_D) !108 Non-fixing cyanobacteria Exponent constant for Dissolved oxygen stress
    call para_get_value('CYN_N_TO_C'                      ,                       CYN_N_TO_C) !109 Non-fixing cyanobacteria Nitrogen to Carbon ratio ,was 0.1
    call para_get_value('CYN_P_TO_C'                      ,                       CYN_P_TO_C) !110 Non-fixing cyanobacteria Phosphorus to Carbon ratio
    call para_get_value('CYN_O2_TO_C'                     ,                      CYN_O2_TO_C) !111 Non-fixing cyanobacteria Oxygen to Carbon ratio for respiration
    call para_get_value('CYN_C_TO_CHLA'                   ,                    CYN_C_TO_CHLA) !112 Non-fixing cyanobacteria Carbon to Chlorophyl a ratio
    call para_get_value('KG_FIX_CYN_OPT_TEMP'             ,              KG_FIX_CYN_OPT_TEMP) !113 Fixing cyanobacteria Growth rate
    call para_get_value('FIX_CYN_OPT_TEMP_LR'             ,              FIX_CYN_OPT_TEMP_LR) !114 Fixing Cyanobacteria optimal temperature lower range
    call para_get_value('FIX_CYN_OPT_TEMP_UR'             ,              FIX_CYN_OPT_TEMP_UR) !115 Fixing Cyanobacteria optimal temperature upper range
    call para_get_value('EFF_FIX_CYN_GROWTH'              ,               EFF_FIX_CYN_GROWTH) !116 Fixing cyanobacteria Effective growth. (1-EG)*growth - losses for RESP and excretion
    call para_get_value('KAPPA_FIX_CYN_UNDER_OPT_TEMP'    ,     KAPPA_FIX_CYN_UNDER_OPT_TEMP) !117 Fixing cyanobacteria Temperature correction for growth lower temperature
    call para_get_value('KAPPA_FIX_CYN_OVER_OPT_TEMP'     ,      KAPPA_FIX_CYN_OVER_OPT_TEMP) !118 Fixing cyanobacteria Temperature correction for growth upper temperature
    call para_get_value('KR_FIX_CYN_20'                   ,                    KR_FIX_CYN_20) !119 Fixing cyanobacteria RESP rate
    call para_get_value('THETA_KR_FIX_CYN'                ,                 THETA_KR_FIX_CYN) !120 Fixing cyanobacteria Temperature correction for RESP rate
    call para_get_value('KD_FIX_CYN_20'                   ,                    KD_FIX_CYN_20) !121 Fixing cyanobacteria Mortality rate of nitrification bacteria
    call para_get_value('THETA_KD_FIX_CYN'                ,                 THETA_KD_FIX_CYN) !122 Fixing cyanobacteria Temperature correction for Mortality rate
    call para_get_value('KHS_DIN_FIX_CYN'                 ,                  KHS_DIN_FIX_CYN) !123 Fixing cyanobacteria Half saturation growth for DIN
    call para_get_value('KHS_DIP_FIX_CYN'                 ,                  KHS_DIP_FIX_CYN) !124 Fixing cyanobacteria Half saturation growth for DIP
    call para_get_value('KHS_O2_FIX_CYN'                  ,                   KHS_O2_FIX_CYN) !125 Fixing cyanobacteria Half saturation growth for O2
    call para_get_value('FRAC_FIX_CYN_EXCR'               ,                FRAC_FIX_CYN_EXCR) !126 Fixing cyanobacteria Fraction of excretion in metabolism rate
    call para_get_value('I_S_FIX_CYN'                     ,                      I_S_FIX_CYN) !127 Fixing cyanobacteria Light saturation (langleys)
    call para_get_value('DO_STR_HYPOX_FIX_CYN_D'          ,           DO_STR_HYPOX_FIX_CYN_D) !128 Fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
    call para_get_value('THETA_HYPOX_FIX_CYN_D'           ,            THETA_HYPOX_FIX_CYN_D) !129 Fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
    call para_get_value('EXPON_HYPOX_FIX_CYN_D'           ,            EXPON_HYPOX_FIX_CYN_D) !130 Fixing cyanobacteria Exponent constant for Dissolved oxygen stress
    call para_get_value('FIX_CYN_N_TO_C'                  ,                   FIX_CYN_N_TO_C) !131 Fixing cyanobacteria Nitrogen to Carbon ratio
    call para_get_value('FIX_CYN_P_TO_C'                  ,                   FIX_CYN_P_TO_C) !132 Fixing cyanobacteria Phosphorus to Carbon ratio
    call para_get_value('FIX_CYN_O2_TO_C'                 ,                  FIX_CYN_O2_TO_C) !133 Fixing cyanobacteria Oxygen to Carbon ratio for respiration
    call para_get_value('FIX_CYN_C_TO_CHLA'               ,                FIX_CYN_C_TO_CHLA) !134 Fixing cyanobacteria Carbon to Chlorophyl a ratio
    call para_get_value('R_FIX'                           ,                            R_FIX) !135 Fixing cyanobacteria Ratio between non-fixing and fixing fractions growth rate
    call para_get_value('K_FIX'                           ,                            K_FIX) !136 Fixing cyanobacteria Effectivity parameter of switching to nitrogen fixation
    call para_get_value('KG_OPA_OPT_TEMP'                 ,                  KG_OPA_OPT_TEMP) !137 OtherPhyto Growth rate
    call para_get_value('OPA_OPT_TEMP_LR'                 ,                  OPA_OPT_TEMP_LR) !138 OtherPhyto optimal temperature lower range
    call para_get_value('OPA_OPT_TEMP_UR'                 ,                  OPA_OPT_TEMP_UR) !139 OtherPhyto optimal temperature upper range
    call para_get_value('EFF_OPA_GROWTH'                  ,                   EFF_OPA_GROWTH) !140 OtherPhyto Effective growth. (1-EG)*growth - losses for respiration and excretion
    call para_get_value('KAPPA_OPA_UNDER_OPT_TEMP'        ,         KAPPA_OPA_UNDER_OPT_TEMP) !141 OtherPhyto Temperature correction for growth lower temperature
    call para_get_value('KAPPA_OPA_OVER_OPT_TEMP'         ,          KAPPA_OPA_OVER_OPT_TEMP) !142 OtherPhyto Temperature correction for growth upper temperature
    call para_get_value('KR_OPA_20'                       ,                        KR_OPA_20) !143 OtherPhyto Respiration rate
    call para_get_value('THETA_KR_OPA'                    ,                     THETA_KR_OPA) !144 OtherPhyto Temperature correction for respiration rate
    call para_get_value('KD_OPA_20'                       ,                        KD_OPA_20) !145 OtherPhyto Mortality rate
    call para_get_value('THETA_KD_OPA'                    ,                     THETA_KD_OPA) !146 OtherPhyto Temperature correction for Mortality rate
    call para_get_value('KHS_DIN_OPA'                     ,                      KHS_DIN_OPA) !147 OtherPhyto Half saturation growth for DIN
    call para_get_value('KHS_DIP_OPA'                     ,                      KHS_DIP_OPA) !148 OtherPhyto Half saturation growth for DIP
    call para_get_value('KHS_O2_OPA'                      ,                       KHS_O2_OPA) !149 OtherPhyto Half saturation growth for O2
    call para_get_value('FRAC_OPA_EXCR'                   ,                    FRAC_OPA_EXCR) !150 OtherPhyto Fraction of excretion in metabolism rate
    call para_get_value('I_S_OPA'                         ,                          I_S_OPA) !151 OtherPhyto Light saturation (langleys)
    call para_get_value('DO_STR_HYPOX_OPA_D'              ,               DO_STR_HYPOX_OPA_D) !152 OtherPhyto Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
    call para_get_value('THETA_HYPOX_OPA_D'               ,                THETA_HYPOX_OPA_D) !153 OtherPhyto Multiplier of the exponent for Dissolved oxygen stress
    call para_get_value('EXPON_HYPOX_OPA_D'               ,                EXPON_HYPOX_OPA_D) !154 OtherPhyto Exponent constant for Dissolved oxygen stress
    call para_get_value('OPA_N_TO_C'                      ,                       OPA_N_TO_C) !155 OtherPhyto Nitrogen to Carbon ratio
    call para_get_value('OPA_P_TO_C'                      ,                       OPA_P_TO_C) !156 OtherPhyto Phosphorus to Carbon ratio
    call para_get_value('OPA_O2_TO_C'                     ,                      OPA_O2_TO_C) !157 OtherPhyto Oxygen to Carbon ratio for respiration
    call para_get_value('OPA_C_TO_CHLA'                   ,                    OPA_C_TO_CHLA) !158 OtherPhyto Carbon to Chlorophyl a ratio
    call para_get_value('KG_ZOO_OPT_TEMP'                 ,                  KG_ZOO_OPT_TEMP) !159 Zooplankton Growth rate
    call para_get_value('ZOO_OPT_TEMP_LR'                 ,                  ZOO_OPT_TEMP_LR) !160 Zooplankton optimal temperature lower range
    call para_get_value('ZOO_OPT_TEMP_UR'                 ,                  ZOO_OPT_TEMP_UR) !161 Zooplankton optimal temperature upper range
    call para_get_value('EFF_ZOO_GROWTH'                  ,                   EFF_ZOO_GROWTH) !162 Zooplankton Effective growth. (1-EG)*growth - losses for respiration and excretion
    call para_get_value('KAPPA_ZOO_UNDER_OPT_TEMP'        ,         KAPPA_ZOO_UNDER_OPT_TEMP) !163 Zooplankton Temperature correction for growth lower temperature
    call para_get_value('KAPPA_ZOO_OVER_OPT_TEMP'         ,          KAPPA_ZOO_OVER_OPT_TEMP) !164 Zooplankton Temperature correction for growth upper temperature
    call para_get_value('GRAT_ZOO_DIA'                    ,                     GRAT_ZOO_DIA) !165 Zooplankton Grazing rate (growhth rate multiplier) on diatoms
    call para_get_value('GRAT_ZOO_CYN'                    ,                     GRAT_ZOO_CYN) !166 Zooplankton Grazing rate (growhth rate multiplier) on Cyanobacteria
    call para_get_value('GRAT_ZOO_OPA'                    ,                     GRAT_ZOO_OPA) !167 Zooplankton Grazing rate (growhth rate multiplier) on fixing Cyanobacteria
    call para_get_value('GRAT_ZOO_FIX_CYN'                ,                 GRAT_ZOO_FIX_CYN) !168 Zooplankton Grazing rate (growhth rate multiplier) on OtherPhyto
    call para_get_value('GRAT_ZOO_DET_PART_ORG_C'         ,          GRAT_ZOO_DET_PART_ORG_C) !172 Zooplankton Grazing rate (growhth rate multiplier) on part. ORG_C
    call para_get_value('PREF_ZOO_DIA'                    ,                     PREF_ZOO_DIA) !173 Zooplankton Preference for diatoms
    call para_get_value('PREF_ZOO_CYN'                    ,                     PREF_ZOO_CYN) !174 Zooplankton Preference for Cyanobacteria
    call para_get_value('PREF_ZOO_FIX_CYN'                ,                 PREF_ZOO_FIX_CYN) !175 Zooplankton Preference for fixing Cyanobacteria
    call para_get_value('PREF_ZOO_OPA'                    ,                     PREF_ZOO_OPA) !176 Zooplankton Preference for OtherPhyto
    call para_get_value('PREF_ZOO_DET_PART_ORG_C'         ,          PREF_ZOO_DET_PART_ORG_C) !180 Zooplankton Preference for part. ORG_C
    call para_get_value('KHS_DIA_C_ZOO'                   ,                    KHS_DIA_C_ZOO) !181 Zooplankton Half saturation growth for diatoms
    call para_get_value('KHS_CYN_C_ZOO'                   ,                    KHS_CYN_C_ZOO) !182 Zooplankton Half saturation growth for Cyanobacteria
    call para_get_value('KHS_FIX_CYN_C_ZOO'               ,                KHS_FIX_CYN_C_ZOO) !183 Zooplankton Half saturation growth for fixing Cyanobacteria
    call para_get_value('KHS_OPA_C_ZOO'                   ,                    KHS_OPA_C_ZOO) !184 Zooplankton Half saturation growth for OtherPhyto
    call para_get_value('KHS_DET_PART_ORG_C_ZOO'          ,           KHS_DET_PART_ORG_C_ZOO) !188 Zooplankton Half saturation growth for part. ORG_C
    call para_get_value('FOOD_MIN_ZOO'                    ,                     FOOD_MIN_ZOO) !189 Zooplankton Minimum food conc. for feeding
    call para_get_value('KE_ZOO'                          ,                           KE_ZOO) !190 not used Zooplankton Excretion rate as growth fraction
    call para_get_value('FRAC_ZOO_EX_ORG'                 ,                  FRAC_ZOO_EX_ORG) !191 not used Zooplankton Excretion rate organic fraction
    call para_get_value('KR_ZOO_20'                       ,                        KR_ZOO_20) !192 Zooplankton Respiration rate
    call para_get_value('THETA_KR_ZOO'                    ,                     THETA_KR_ZOO) !193 Zooplankton Respiration rate Temperature correction
    call para_get_value('KD_ZOO_20'                       ,                        KD_ZOO_20) !194 Zooplankton Mortality rate
    call para_get_value('THETA_KD_ZOO'                    ,                     THETA_KD_ZOO) !195 Zooplankton Mortality rate Temperature correction
    call para_get_value('DO_STR_HYPOX_ZOO_D'              ,               DO_STR_HYPOX_ZOO_D) !196 Zooplankton Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
    call para_get_value('THETA_HYPOX_ZOO_D'               ,                THETA_HYPOX_ZOO_D) !197 Zooplankton Multiplier of the exponent for Dissolved oxygen stress
    call para_get_value('EXPON_HYPOX_ZOO_D'               ,                EXPON_HYPOX_ZOO_D) !198 Zooplankton Exponent constant for Dissolved oxygen stress
    call para_get_value('ZOO_N_TO_C'                      ,                       ZOO_N_TO_C) !199 Zooplankton Nitrogen to Carbon ratio
    call para_get_value('ZOO_P_TO_C'                      ,                       ZOO_P_TO_C) !200 Zooplankton Phosphorus to Carbon ratio
    call para_get_value('ZOO_O2_TO_C'                     ,                      ZOO_O2_TO_C) !201 Zooplankton Oxygen to Carbon ratio for respiration
    call para_get_value('KDISS_DET_PART_ORG_C_20'         ,          KDISS_DET_PART_ORG_C_20) !202 Particulate Detritus Carbon Dissolution rate not dependent on phytoplankton
    call para_get_value('THETA_KDISS_DET_PART_ORG_C'      ,       THETA_KDISS_DET_PART_ORG_C) !203 Particulate Detritus Carbon Dissolution rate Temperature correction
    call para_get_value('FAC_PHYT_DET_PART_ORG_C'         ,          FAC_PHYT_DET_PART_ORG_C) !204 Particulate Detritus Carbon Phytoplankton linear factor for dissolution rate
    call para_get_value('KDISS_DET_PART_ORG_N_20'         ,          KDISS_DET_PART_ORG_N_20) !205 Particulate Detritus Nitrogen Dissolution rate not dependent on phytoplankton
    call para_get_value('THETA_KDISS_DET_PART_ORG_N'      ,       THETA_KDISS_DET_PART_ORG_N) !206 Particulate Detritus Nitrogen Dissolution rate Temperature correction
    call para_get_value('KHS_DISS_N'                      ,                       KHS_DISS_N) !207 Particulate Detritus Nitrogen dissolution reverse half saturation for DIN
    call para_get_value('FAC_PHYT_DET_PART_ORG_N'         ,          FAC_PHYT_DET_PART_ORG_N) !208 Particulate Detritus Nitrogen Phytoplankton linear factor for dissolution rate
    call para_get_value('KDISS_DET_PART_ORG_P_20'         ,          KDISS_DET_PART_ORG_P_20) !209 Particulate Detritus Phosphorus Dissolution rate not dependent on phytoplankton
    call para_get_value('THETA_KDISS_DET_PART_ORG_P'      ,       THETA_KDISS_DET_PART_ORG_P) !210 Particulate Detritus Phosphorus Dissolution rate Temperature correction
    call para_get_value('KHS_DISS_P'                      ,                       KHS_DISS_P) !211 Particulate Detritus Phosphorus  dissolution reverse half saturation for DIP
    call para_get_value('FAC_PHYT_DET_PART_ORG_P'         ,          FAC_PHYT_DET_PART_ORG_P) !212 Particulate Detritus Phosphorus  Phytoplankton linear factor for dissolution rate
    call para_get_value('KDISS_PART_Si_20'                ,                 KDISS_PART_Si_20) !213 Particulate Silica Dissolution rate
    call para_get_value('THETA_KDISS_PART_Si'             ,              THETA_KDISS_PART_Si) !214 Particulate Silica Dissolution rate Temperature correction
    call para_get_value('FAC_PHYT_AMIN_DOC'               ,                FAC_PHYT_AMIN_DOC) !217 Dissolved carbon  Phytoplankton linear factor for mineralisation rate
    call para_get_value('KHS_AMIN_N'                      ,                       KHS_AMIN_N) !220 Dissolved nitrogen  reverse half saturation for DIN
    call para_get_value('FAC_PHYT_AMIN_DON'               ,                FAC_PHYT_AMIN_DON) !221 Dissolved nitrogen Phytoplankton linear factor for mineralisation rate
    call para_get_value('KHS_AMIN_P'                      ,                       KHS_AMIN_P) !224 Dissolved phosphorus reverse half saturation for DIP
    call para_get_value('FAC_PHYT_AMIN_DOP'               ,                FAC_PHYT_AMIN_DOP) !225 Dissolved phosphorus Phytoplankton linear factor for mineralisation rate
    call para_get_value('K_NITR_20'                       ,                        K_NITR_20) !226 Amonia nitrification rate
    call para_get_value('KHS_NITR_OXY'                    ,                     KHS_NITR_OXY) !227 Amonia nitrification half saturation for Oxygen
    call para_get_value('KHS_NITR_NH4_N'                  ,                   KHS_NITR_NH4_N) !228 Amonia nitrification half saturation for Amonia
    call para_get_value('THETA_K_NITR'                    ,                     THETA_K_NITR) !229 Amonia nitrification rate Temperature constant

    call para_get_value('PH_NITR_NH4_MIN'                 ,            PH_NITR_NH4_MIN  )  !236   optimum lower range for pH correction factor for nitrification
    call para_get_value('PH_NITR_NH4_MAX'                 ,            PH_NITR_NH4_MAX  )  !237   optimum upper range for pH correction factor for nitrification

    call para_get_value('k_OX_FE_II'                      ,          k_OX_FE_II         )  !244!    Oxidation rate for iron 2+
    call para_get_value('k_RED_FE_III'                    ,          k_RED_FE_III       )  !245!    reduction rate for iron 3+
    call para_get_value('k_OX_MN_II'                      ,          k_OX_MN_II         )  !246!    oxidation rate for manganese 2+
    call para_get_value('k_RED_MN_IV'                     ,          k_RED_MN_IV        )  !247!    reduction rate for manganese 4+
    call para_get_value('KHS_DOXY_FE_III_RED'             ,          KHS_DOXY_FE_III_RED)  !248!    reversed Monod half saturation of DOXY for iron 3+ reduction
    call para_get_value('KHS_DOXY_MN_IV_RED'              ,          KHS_DOXY_MN_IV_RED )  !249!    reversed Monod half saturation of DOXY for manganese 4+ reduction

    ! New model constants introduced 27 January 2016 for the redox sequences:
    call para_get_value('K_MIN_DOC_DOXY_20'       ,  K_MIN_DOC_DOXY_20       )
    call para_get_value('K_MIN_DOC_NO3N_20'       ,  K_MIN_DOC_NO3N_20       )
    call para_get_value('K_MIN_DOC_MN_IV_20'      ,  K_MIN_DOC_MN_IV_20      )
    call para_get_value('K_MIN_DOC_FE_III_20'     ,  K_MIN_DOC_FE_III_20     )
    call para_get_value('K_MIN_DOC_S_PLUS_6_20'   ,  K_MIN_DOC_S_PLUS_6_20   )
    call para_get_value('K_MIN_DOC_DOC_20'        ,  K_MIN_DOC_DOC_20        )
    call para_get_value('THETA_K_MIN_DOC_DOXY'    ,  THETA_K_MIN_DOC_DOXY    )
    call para_get_value('THETA_K_MIN_DOC_NO3N'    ,  THETA_K_MIN_DOC_NO3N    )
    call para_get_value('THETA_K_MIN_DOC_MN_IV'   ,  THETA_K_MIN_DOC_MN_IV   )
    call para_get_value('THETA_K_MIN_DOC_FE_III'  ,  THETA_K_MIN_DOC_FE_III  )
    call para_get_value('THETA_K_MIN_DOC_S_PLUS_6',  THETA_K_MIN_DOC_S_PLUS_6)
    call para_get_value('THETA_K_MIN_DOC_DOC'     ,  THETA_K_MIN_DOC_DOC     )
    call para_get_value('K_HS_DOC_MIN_DOXY'       ,  K_HS_DOC_MIN_DOXY       )
    call para_get_value('K_HS_DOC_MIN_NO3N'       ,  K_HS_DOC_MIN_NO3N       )
    call para_get_value('K_HS_DOC_MIN_MN_IV'      ,  K_HS_DOC_MIN_MN_IV      )
    call para_get_value('K_HS_DOC_MIN_FE_III'     ,  K_HS_DOC_MIN_FE_III     )
    call para_get_value('K_HS_DOC_MIN_S_PLUS_6'   ,  K_HS_DOC_MIN_S_PLUS_6   )
    call para_get_value('K_HS_DOC_MIN_DOC'        ,  K_HS_DOC_MIN_DOC        )
    call para_get_value('K_HS_DOXY_RED_LIM'       ,  K_HS_DOXY_RED_LIM       )
    call para_get_value('K_HS_NO3N_RED_LIM'       ,  K_HS_NO3N_RED_LIM       )
    call para_get_value('K_HS_MN_IV_RED_LIM'      ,  K_HS_MN_IV_RED_LIM      )
    call para_get_value('K_HS_FE_III_RED_LIM'     ,  K_HS_FE_III_RED_LIM     )
    call para_get_value('K_HS_S_PLUS_6_RED_LIM'   ,  K_HS_S_PLUS_6_RED_LIM   )
    call para_get_value('K_HS_DOXY_RED_INHB'      ,  K_HS_DOXY_RED_INHB      )
    call para_get_value('K_HS_NO3N_RED_INHB'      ,  K_HS_NO3N_RED_INHB      )
    call para_get_value('K_HS_MN_IV_RED_INHB'     ,  K_HS_MN_IV_RED_INHB     )
    call para_get_value('K_HS_FE_III_RED_INHB'    ,  K_HS_FE_III_RED_INHB    )
    call para_get_value('K_HS_S_PLUS_6_RED_INHB'  ,  K_HS_S_PLUS_6_RED_INHB  )
    call para_get_value('PH_MIN_DOC_MIN_DOXY'     ,  PH_MIN_DOC_MIN_DOXY     )
    call para_get_value('PH_MIN_DOC_MIN_NO3N'     ,  PH_MIN_DOC_MIN_NO3N     )
    call para_get_value('PH_MIN_DOC_MIN_MN_IV'    ,  PH_MIN_DOC_MIN_MN_IV    )
    call para_get_value('PH_MIN_DOC_MIN_FE_III'   ,  PH_MIN_DOC_MIN_FE_III   )
    call para_get_value('PH_MIN_DOC_MIN_S_PLUS_6' ,  PH_MIN_DOC_MIN_S_PLUS_6 )
    call para_get_value('PH_MIN_DOC_MIN_DOC'      ,  PH_MIN_DOC_MIN_DOC      )
    call para_get_value('PH_MAX_DOC_MIN_DOXY'     ,  PH_MAX_DOC_MIN_DOXY     )
    call para_get_value('PH_MAX_DOC_MIN_NO3N'     ,  PH_MAX_DOC_MIN_NO3N     )
    call para_get_value('PH_MAX_DOC_MIN_MN_IV'    ,  PH_MAX_DOC_MIN_MN_IV    )
    call para_get_value('PH_MAX_DOC_MIN_FE_III'   ,  PH_MAX_DOC_MIN_FE_III   )
    call para_get_value('PH_MAX_DOC_MIN_S_PLUS_6' ,  PH_MAX_DOC_MIN_S_PLUS_6 )
    call para_get_value('PH_MAX_DOC_MIN_DOC'      ,  PH_MAX_DOC_MIN_DOC      )

    ! New model constants introduced 28 January 2016 for the redox sequences
    call para_get_value('K_MIN_DON_DOXY_20'         ,    K_MIN_DON_DOXY_20         )
    call para_get_value('K_MIN_DON_NO3N_20'         ,    K_MIN_DON_NO3N_20         )
    call para_get_value('K_MIN_DON_MN_IV_20'        ,    K_MIN_DON_MN_IV_20        )
    call para_get_value('K_MIN_DON_FE_III_20'       ,    K_MIN_DON_FE_III_20       )
    call para_get_value('K_MIN_DON_S_PLUS_6_20'     ,    K_MIN_DON_S_PLUS_6_20     )
    call para_get_value('K_MIN_DON_DOC_20'          ,    K_MIN_DON_DOC_20          )
    call para_get_value('THETA_K_MIN_DON_DOXY'      ,    THETA_K_MIN_DON_DOXY      )
    call para_get_value('THETA_K_MIN_DON_NO3N'      ,    THETA_K_MIN_DON_NO3N      )
    call para_get_value('THETA_K_MIN_DON_MN_IV'     ,    THETA_K_MIN_DON_MN_IV     )
    call para_get_value('THETA_K_MIN_DON_FE_III'    ,    THETA_K_MIN_DON_FE_III    )
    call para_get_value('THETA_K_MIN_DON_S_PLUS_6'  ,    THETA_K_MIN_DON_S_PLUS_6  )
    call para_get_value('THETA_K_MIN_DON_DOC'       ,    THETA_K_MIN_DON_DOC       )
    call para_get_value('K_HS_DON_MIN_DOXY'         ,    K_HS_DON_MIN_DOXY         )
    call para_get_value('K_HS_DON_MIN_NO3N'         ,    K_HS_DON_MIN_NO3N         )
    call para_get_value('K_HS_DON_MIN_MN_IV'        ,    K_HS_DON_MIN_MN_IV        )
    call para_get_value('K_HS_DON_MIN_FE_III'       ,    K_HS_DON_MIN_FE_III       )
    call para_get_value('K_HS_DON_MIN_S_PLUS_6'     ,    K_HS_DON_MIN_S_PLUS_6     )
    call para_get_value('K_HS_DON_MIN_DOC'          ,    K_HS_DON_MIN_DOC          )
    call para_get_value('PH_MIN_DON_MIN_DOXY'       ,    PH_MIN_DON_MIN_DOXY       )
    call para_get_value('PH_MIN_DON_MIN_NO3N'       ,    PH_MIN_DON_MIN_NO3N       )
    call para_get_value('PH_MIN_DON_MIN_MN_IV'      ,    PH_MIN_DON_MIN_MN_IV      )
    call para_get_value('PH_MIN_DON_MIN_FE_III'     ,    PH_MIN_DON_MIN_FE_III     )
    call para_get_value('PH_MIN_DON_MIN_S_PLUS_6'   ,    PH_MIN_DON_MIN_S_PLUS_6   )
    call para_get_value('PH_MIN_DON_MIN_DOC'        ,    PH_MIN_DON_MIN_DOC        )
    call para_get_value('PH_MAX_DON_MIN_DOXY'       ,    PH_MAX_DON_MIN_DOXY       )
    call para_get_value('PH_MAX_DON_MIN_NO3N'       ,    PH_MAX_DON_MIN_NO3N       )
    call para_get_value('PH_MAX_DON_MIN_MN_IV'      ,    PH_MAX_DON_MIN_MN_IV      )
    call para_get_value('PH_MAX_DON_MIN_FE_III'     ,    PH_MAX_DON_MIN_FE_III     )
    call para_get_value('PH_MAX_DON_MIN_S_PLUS_6'   ,    PH_MAX_DON_MIN_S_PLUS_6   )
    call para_get_value('PH_MAX_DON_MIN_DOC'        ,    PH_MAX_DON_MIN_DOC        )
    call para_get_value('K_MIN_DOP_DOXY_20'         ,    K_MIN_DOP_DOXY_20         )
    call para_get_value('K_MIN_DOP_NO3N_20'         ,    K_MIN_DOP_NO3N_20         )
    call para_get_value('K_MIN_DOP_MN_IV_20'        ,    K_MIN_DOP_MN_IV_20        )
    call para_get_value('K_MIN_DOP_FE_III_20'       ,    K_MIN_DOP_FE_III_20       )
    call para_get_value('K_MIN_DOP_S_PLUS_6_20'     ,    K_MIN_DOP_S_PLUS_6_20     )
    call para_get_value('K_MIN_DOP_DOC_20'          ,    K_MIN_DOP_DOC_20          )
    call para_get_value('THETA_K_MIN_DOP_DOXY'      ,    THETA_K_MIN_DOP_DOXY      )
    call para_get_value('THETA_K_MIN_DOP_NO3N'      ,    THETA_K_MIN_DOP_NO3N      )
    call para_get_value('THETA_K_MIN_DOP_MN_IV'     ,    THETA_K_MIN_DOP_MN_IV     )
    call para_get_value('THETA_K_MIN_DOP_FE_III'    ,    THETA_K_MIN_DOP_FE_III    )
    call para_get_value('THETA_K_MIN_DOP_S_PLUS_6'  ,    THETA_K_MIN_DOP_S_PLUS_6  )
    call para_get_value('THETA_K_MIN_DOP_DOC'       ,    THETA_K_MIN_DOP_DOC       )
    call para_get_value('K_HS_DOP_MIN_DOXY'         ,    K_HS_DOP_MIN_DOXY         )
    call para_get_value('K_HS_DOP_MIN_NO3N'         ,    K_HS_DOP_MIN_NO3N         )
    call para_get_value('K_HS_DOP_MIN_MN_IV'        ,    K_HS_DOP_MIN_MN_IV        )
    call para_get_value('K_HS_DOP_MIN_FE_III'       ,    K_HS_DOP_MIN_FE_III       )
    call para_get_value('K_HS_DOP_MIN_S_PLUS_6'     ,    K_HS_DOP_MIN_S_PLUS_6     )
    call para_get_value('K_HS_DOP_MIN_DOC'          ,    K_HS_DOP_MIN_DOC          )
    call para_get_value('PH_MIN_DOP_MIN_DOXY'       ,    PH_MIN_DOP_MIN_DOXY       )
    call para_get_value('PH_MIN_DOP_MIN_NO3N'       ,    PH_MIN_DOP_MIN_NO3N       )
    call para_get_value('PH_MIN_DOP_MIN_MN_IV'      ,    PH_MIN_DOP_MIN_MN_IV      )
    call para_get_value('PH_MIN_DOP_MIN_FE_III'     ,    PH_MIN_DOP_MIN_FE_III     )
    call para_get_value('PH_MIN_DOP_MIN_S_PLUS_6'   ,    PH_MIN_DOP_MIN_S_PLUS_6   )
    call para_get_value('PH_MIN_DOP_MIN_DOC'        ,    PH_MIN_DOP_MIN_DOC        )
    call para_get_value('PH_MAX_DOP_MIN_DOXY'       ,    PH_MAX_DOP_MIN_DOXY       )
    call para_get_value('PH_MAX_DOP_MIN_NO3N'       ,    PH_MAX_DOP_MIN_NO3N       )
    call para_get_value('PH_MAX_DOP_MIN_MN_IV'      ,    PH_MAX_DOP_MIN_MN_IV      )
    call para_get_value('PH_MAX_DOP_MIN_FE_III'     ,    PH_MAX_DOP_MIN_FE_III     )
    call para_get_value('PH_MAX_DOP_MIN_S_PLUS_6'   ,    PH_MAX_DOP_MIN_S_PLUS_6   )
    call para_get_value('PH_MAX_DOP_MIN_DOC'        ,    PH_MAX_DOP_MIN_DOC        )

    ! New model constats added in 29 th of January 2016
    call para_get_value('k_OX_CH4'                  ,     k_OX_CH4                 )
    call para_get_value('THETA_k_OX_CH4'            ,     THETA_k_OX_CH4           )
    call para_get_value('k_HS_OX_CH4_DOXY'          ,     k_HS_OX_CH4_DOXY         )
    call para_get_value('k_OX_H2S'                  ,     k_OX_H2S                 )
    call para_get_value('THETA_k_OX_H2S'            ,     THETA_k_OX_H2S           )
    call para_get_value('k_HS_OX_H2S_DOXY'          ,     k_HS_OX_H2S_DOXY         )

    ! New model constants added in 9th August 2016
    call para_get_value('k_DISS_FE_II_20'           ,    k_DISS_FE_II_20           )
    call para_get_value('THETA_k_DISS_FE_II'        ,    THETA_k_DISS_FE_II        )
    call para_get_value('INIT_MULT_FE_II_DISS'      ,    INIT_MULT_FE_II_DISS      )

    call para_get_value('k_DISS_FE_III_20'          ,    k_DISS_FE_III_20          )
    call para_get_value('THETA_k_DISS_FE_III'       ,    THETA_k_DISS_FE_III       )
    call para_get_value('INIT_MULT_FE_III_DISS'     ,    INIT_MULT_FE_III_DISS     )
    ! End of new model constants added in 9th August 2016
end subroutine INIT_PELAGIC_MODEL_CONSTANTS



!******************************************************************
!******************************************************************
!******************************************************************

!******************************************************************
! 30th of Novemver 2015
!
! Initial additions for interfacing dissolved and particulate
! fractions of FE_II, FE_IV, MN_II, MN_IV
!
!******************************************************************

!******************************************************************
!
! Change by Ali Erturk 25th of November 2016
!
!          Initial version of balance beased aquatic chemistry model
!          just for iron
!
! Change by Ali Erturk and Petras Zemlys 27th of November 2016
!
!         - Redox sequence
!
!         - New state variables added
!
!              - CA        (Calcium)
!              - MG        (Magnesium)
!              - S_PLUS_6  (Sulphate sulphur)
!              - S_MINUS_2 (Sulphide sulphur)
!              - CH4_C      (Methane carbon)
!
! Changes by Ali Erturk for KUMADUBI 3 rd of July 2016
!
! 9 September 2018
! Changes by Ali Erturk related to clean the code, where all 
! useless bacteria state variables and related processes are
! removed
!
!*******************************************************************


subroutine PELAGIC_KINETICS &
           (node_active      ,           nkn, &
            STATE_VARIABLES  , DERIVATIVES  , nstate, &
            MODEL_CONSTANTS  , nconst               , &
            DRIVING_FUNCTIONS, n_driving_functions  , &
            FLAGS            , nflags               , &
            PROCESS_RATES    , NDIAGVAR             , &
            SAVED_OUTPUTS    , n_saved_outputs      , &
            PH               , &
            TIME, TIME_STEP  , CALLED_BEFORE        , &
            SURFACE_BOXES    , ZOOP_OPTION_1        , &
            ADVANCED_REDOX_OPTION)

    use CO2SYS_CDIAC
    use AQUABC_II_GLOBAL
    use PELAGIC_MODEL_CONSTANTS
	!use para_aqua
    use aquabc_II_wc_ini
    !use basin, only: ipv !0d correction

    use AQUABC_PEL_STATE_VAR_INDEXES
    
    implicit none
    !include 'param.h'

    !integer ipv(nkndim)	       !external node numbers
    !common  /ipv/ipv

    logical VALUE_strange(nkn) !For NaN and Inf checking

    integer  :: nkn, nstate, nconst, n_driving_functions, nflags
    integer  :: n_saved_outputs, NDIAGVAR, STRANGERSD, error

    integer node_active(nkn)   !internal numbers of nodes, used for diagnostics. Not implemented yet

    double precision, dimension(nkn,nstate),             intent(in)   :: STATE_VARIABLES
    double precision, dimension(nkn,nstate),             intent(out)  :: DERIVATIVES
    double precision, dimension(nconst),                 intent(in)   :: MODEL_CONSTANTS
    double precision, dimension(nkn,n_driving_functions),intent(in)   :: DRIVING_FUNCTIONS
    integer,          dimension(nflags),                 intent(in)   :: FLAGS
    double precision, dimension(nkn,nstate, NDIAGVAR),   intent(out)  :: PROCESS_RATES
    double precision, dimension(nkn) :: PH
    
    !For saving some variables to be used for all nodes?
    double precision, dimension(nkn,n_saved_outputs),    intent(inout):: SAVED_OUTPUTS
    double precision, intent(in) :: TIME
    double precision, intent(in) :: TIME_STEP

    integer, intent(inout)        :: CALLED_BEFORE

    !Optional arguments that are needed for ESTAS implementation but
    !kept optional to allow SHYFEM AQUABC compatibility
    integer, dimension(nkn), intent(in), optional :: SURFACE_BOXES
    integer                , intent(in), optional :: ZOOP_OPTION_1
    integer                , intent(in), optional :: ADVANCED_REDOX_OPTION
    
    integer :: i,j,k

    !*********************************************'
    !*                                           *'
    !* PELAGIC ECOLOGY KINETICS
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

    ! New state variables added 9 September 2015
    real(kind = DBL_PREC), dimension(nkn) :: FE_II
    real(kind = DBL_PREC), dimension(nkn) :: FE_III
    real(kind = DBL_PREC), dimension(nkn) :: MN_II
    real(kind = DBL_PREC), dimension(nkn) :: MN_IV
    ! End of new state variables added 22 September 2015

    ! New state variables added 27 January 2016
    real(kind = DBL_PREC), dimension(nkn) :: CA
    real(kind = DBL_PREC), dimension(nkn) :: MG
    real(kind = DBL_PREC), dimension(nkn) :: S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn) :: S_MINUS_2
    real(kind = DBL_PREC), dimension(nkn) :: CH4_C
    ! End of new state variables added 27 January 2016

    double precision :: PHYT_TOT_C(nkn)

    !Driving functions
    double precision :: TEMP     (nkn)
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

    ! Introduced by Ali
    integer :: FIRST_TIME_STEP
    integer :: INIT_OPTION_OF_FE_II_DISS
    integer :: INIT_OPTION_OF_FE_III_DISS
    ! End of flags

	double precision :: K_A_CALC(nkn)   !calculated aeration reactor specific coefficient
	
    !Main process rates
    ! New kinetic rates added 9 September 2015
    real(kind = DBL_PREC), dimension(nkn) :: R_FE_II_OXIDATION
    real(kind = DBL_PREC), dimension(nkn) :: R_FE_III_REDUCTION
    real(kind = DBL_PREC), dimension(nkn) :: R_MN_II_OXIDATION
    real(kind = DBL_PREC), dimension(nkn) :: R_MN_IV_REDUCTION
    ! End of new kinetic rates added 9 September 2015

    ! New kinetic rates introduced 27 January 2016 for the redox sequences
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DOC_MIN_DOXY
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DOC_MIN_NO3N
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DOC_MIN_MN_IV
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DOC_MIN_FE_III
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DOC_MIN_DOC
    ! End of new kinetic rates introduced 27 January 2016 for the redox sequences

    ! New kinetic rates introduced 28 January 2016 for the redox sequences
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DON_MIN_DOXY
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DON_MIN_NO3N
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DON_MIN_MN_IV
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DON_MIN_FE_III
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DON_MIN_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DON_MIN_DOC

    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DOP_MIN_DOXY
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DOP_MIN_NO3N
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DOP_MIN_MN_IV
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DOP_MIN_FE_III
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DOP_MIN_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn) :: R_ABIOTIC_DOP_MIN_DOC
    ! End of new kinetic rates introduced 28 January 2016 for the redox sequences

    ! New kinetic rates introduced 29 January 2016 for the redox sequences
    real(kind = DBL_PREC), dimension(nkn) :: R_SULPHATE_REDUCTION
    real(kind = DBL_PREC), dimension(nkn) :: R_SULPHIDE_OXIDATION
    real(kind = DBL_PREC), dimension(nkn) :: H2S_ATM_EXCHANGE
    real(kind = DBL_PREC), dimension(nkn) :: R_METHANOGENESIS
    real(kind = DBL_PREC), dimension(nkn) :: R_METHANE_OXIDATION
    real(kind = DBL_PREC), dimension(nkn) :: CH4_ATM_EXCHANGE
    ! End of new kinetic rates introduced 29 January 2016 for the redox sequences

    real(kind = DBL_PREC), dimension(nkn) :: R_AERATION


    real(kind = DBL_PREC), dimension(nkn) :: R_DIA_GROWTH
    real(kind = DBL_PREC), dimension(nkn) :: R_DIA_MET
    real(kind = DBL_PREC), dimension(nkn) :: R_DIA_RESP
    real(kind = DBL_PREC), dimension(nkn) :: R_DIA_EXCR

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
    double precision :: R_ZOO_FEEDING_DET_PART_ORG_C (nkn)
    double precision :: R_ZOO_RESP                   (nkn)
    double precision :: R_ZOO_INT_RESP               (nkn)
    double precision :: R_ZOO_TOT_RESP               (nkn)
    double precision :: R_ZOO_DEATH                  (nkn)
    double precision :: R_ZOO_EX_DON                 (nkn)
    double precision :: R_ZOO_EX_DOP                 (nkn)
    double precision :: R_ZOO_EX_DOC                 (nkn)

    double precision :: R_DET_PART_ORG_C_DISSOLUTION (nkn)
    double precision :: LIM_PHYT_DISS_DET_PART_ORG_C (nkn)

    double precision :: LIM_N_MIN_DON_N    (nkn)
    double precision :: LIM_P_MIN_DOP_P    (nkn)

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
    double precision :: R_DENITRIFICATION            (nkn)

    double precision :: R_DET_PART_ORG_N_DISSOLUTION (nkn)
    double precision :: LIM_N_DISS_DET_PART_ORG_N    (nkn)

    double precision :: LIM_PHY_N_DISS_DET_PART_ORG_N(nkn)
    double precision :: R_DET_PART_ORG_P_DISSOLUTION (nkn)
    double precision :: LIM_P_DISS_DET_PART_ORG_P    (nkn)

    double precision :: LIM_PHY_P_DISS_DET_PART_ORG_P(nkn)


    !Auxillary variables
    double precision :: DISS_OXYGEN_SAT                (nkn)
    double precision :: CHLA                           (nkn)
    double precision :: K_E                            (nkn)
    double precision :: KG_DIA                         (nkn)
    double precision :: LIM_KG_DIA                     (nkn)
    double precision :: FAC_HYPOX_DIA_D                (nkn)
    double precision :: PREF_NH4N_DIA                  (nkn)
    double precision :: KD_DIA                         (nkn)
    double precision :: KG_ZOO                         (nkn)
    double precision :: KG_ZOO_DIA                     (nkn)
    double precision :: KG_ZOO_DET_PART_ORG_C          (nkn)
    double precision :: FOOD_AVAIL_ZOO                 (nkn)
    double precision :: FOOD_FACTOR_ZOO_DIA            (nkn)
    double precision :: FOOD_FACTOR_ZOO_DET_PART_ORG_C (nkn)
    double precision :: ACTUAL_ZOO_N_TO_C              (nkn)
    double precision :: ACTUAL_ZOO_P_TO_C              (nkn)

    double precision :: FAC_HYPOX_ZOO_D  (nkn)
    double precision :: ACTUAL_DET_N_TO_C(nkn)
    double precision :: ACTUAL_DET_P_TO_C(nkn)
    double precision :: KD_ZOO           (nkn)

    double precision :: KD_FIX_CYN         (nkn)
    double precision :: FAC_HYPOX_FIX_CYN_D(nkn)
    double precision :: PREF_NH4N_DON_FIX_CYN  (nkn)

    double precision :: KG_ZOO_FIX_CYN               (nkn)
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
    double precision :: PREF_NH4N_DON_CYN  (nkn)
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

    !limitation factors
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

    real(kind = DBL_PREC) :: ALPHA_0               (nkn)
    real(kind = DBL_PREC) :: ALPHA_1               (nkn)

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
    character(len=34),    allocatable, dimension (:)   :: CO2SYS_NICEHEADERS

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
    real(kind = DBL_PREC) :: H2CO3                 (nkn)

    !Introduced by Petras instead of using co2sys_outdata directly
    real(kind = DBL_PREC) :: HCO3                  (nkn)
    real(kind = DBL_PREC) :: CO3                   (nkn)


    ! Variables and constansts for Alkalinity?
    real(kind = DBL_PREC) :: PKH                        (nkn)
    real(kind = DBL_PREC) :: FRAC_NH3                   (nkn)
    real(kind = DBL_PREC) :: FRAC_NH4                   (nkn)
    real(kind = DBL_PREC) :: N_DIA_TOT_RESP             (nkn)
    real(kind = DBL_PREC) :: N_CYN_TOT_RESP             (nkn)
    real(kind = DBL_PREC) :: N_OPA_TOT_RESP             (nkn)
    real(kind = DBL_PREC) :: N_FIX_CYN_TOT_RESP         (nkn)
    real(kind = DBL_PREC) :: N_ZOO_TOT_RESP             (nkn)
    real(kind = DBL_PREC) :: N_ABIOTIC_DON_MIN          (nkn)
    real(kind = DBL_PREC) :: ALK_GAINED_BY_AMMONIUM_GEN (nkn)
    real(kind = DBL_PREC) :: N_DENITRIFICATION          (nkn)
    real(kind = DBL_PREC) :: N_DIA_GROWTH               (nkn)
    real(kind = DBL_PREC) :: N_CYN_GROWTH               (nkn)
    real(kind = DBL_PREC) :: N_OPA_GROWTH               (nkn)
    real(kind = DBL_PREC) :: N_NON_FIX_CYN_GROWTH       (nkn)
    real(kind = DBL_PREC) :: ALK_GAINED_BY_NITRATE_CONS (nkn)
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

    real(kind = DBL_PREC) :: P_DIA_TOT_RESP              (nkn)
    real(kind = DBL_PREC) :: P_CYN_TOT_RESP              (nkn)
    real(kind = DBL_PREC) :: P_OPA_TOT_RESP              (nkn)
    real(kind = DBL_PREC) :: P_FIX_CYN_TOT_RESP          (nkn)
    real(kind = DBL_PREC) :: P_ZOO_TOT_RESP              (nkn)
    real(kind = DBL_PREC) :: P_ABIOTIC_DOP_MIN           (nkn)
    real(kind = DBL_PREC) :: P_DIA_GROWTH                (nkn)
    real(kind = DBL_PREC) :: P_CYN_GROWTH                (nkn)
    real(kind = DBL_PREC) :: P_OPA_GROWTH                (nkn)
    real(kind = DBL_PREC) :: P_NON_FIX_CYN_GROWTH        (nkn)


    integer :: CONSIDER_ALKALNITY_DERIVATIVE
    integer :: CONSIDER_INORG_C_DERIVATIVE
    integer :: CONSIDER_CO2_REARATION

    !integer :: NUM_SCREEN_OUTPUT_NODES
    !parameter (NUM_SCREEN_OUTPUT_NODES = 0)
    !integer, dimension(NUM_SCREEN_OUTPUT_NODES) :: SCREEN_OUTPUT_NODES

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

    ! -------------------------------------------------------------------------
    ! Variables added for the new pH correction algorithm
    ! 31st August 2015
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_NITR_NH4
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DENITR_NO3

    ! -------------------------------------------------------------------------
    ! End of variables added for the new pH correction algorithm
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Variables added 30 November 2015 for dissolved and particulate species
    ! of FE_II, FE_III, MN_II, MN_IV. These will be internally calculated by
    ! this subroutine and will be returned in some suitable way to AQUABC since
    ! they will be used in settling calculations. Job by Petras and his fellows.
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn) :: MULT_FE_II_DISS  !Dissolved fraction of FE_II
    real(kind = DBL_PREC), dimension(nkn) :: MULT_FE_II_PART  !Particulate fraction of FE_II
    real(kind = DBL_PREC), dimension(nkn) :: MULT_FE_III_DISS !Dissolved fraction of FE_III
    real(kind = DBL_PREC), dimension(nkn) :: MULT_FE_III_PART !Particulate fraction of FE_III

    real(kind = DBL_PREC), dimension(nkn) :: MULT_MN_II_DISS  !Dissolved fraction of MN_II
    real(kind = DBL_PREC), dimension(nkn) :: MULT_MN_II_PART  !Particulate fraction of MN_II
    real(kind = DBL_PREC), dimension(nkn) :: MULT_MN_IV_DISS  !Dissolved fraction of MN_IV
    real(kind = DBL_PREC), dimension(nkn) :: MULT_MN_IV_PART  !Particulate fraction of MN_IV
    ! -------------------------------------------------------------------------
    ! End of variables added 30 November 2015 for dissolved and particulate
    ! species of FE_II, FE_III, MN_II, MN_IV
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Variables added 25 January 2016 for a simple version of equlibrium
    ! based aquatic chemistry calculations for Fe2+ and Fe3+
    ! -------------------------------------------------------------------------
    !
    ! Variable list:
    ! K_EQ_S_1                    : 1st Dissociaciton constant for H2S
    ! K_EQ_S_2                    : 2nd Dissociaciton constant for H2S
    ! K_SP_FES                    : Solubility pproduct for FeS
    ! HS2_TOT                     : Total H2S in moles, (will be replaced by state variable)
    ! H2S_DIVISOR                 : Auxillary variable
    ! FRAC_H2S_IN_H2S_TOT         : Fraction of S-- in total H2S
    ! FRAC_HS_MINUS_IN_H2S_TOT    : Fraction of HS- in total H2S
    ! FRAC_S_MINUS_TWO_IN_H2S_TOT : Fraction of H2S in total H2S
    ! H2S                         : H2S in moles
    ! HS_MINUS                    : HS- in moles
    ! S_MINUS_TWO                 : S-- in moles
    ! FE_II_DISS                  : Dissolved Fe2+
    ! FE_II_PART                  : Particulate Fe2+
    ! FE_III_DISS                 : Dissolved Fe3+
    ! FE_III_PART                 : Particulate Fe3+

    real(kind = DBL_PREC), dimension(nkn) :: K_EQ_S_1
    real(kind = DBL_PREC), dimension(nkn) :: K_EQ_S_2
    real(kind = DBL_PREC), dimension(nkn) :: K_SP_FES
    real(kind = DBL_PREC), dimension(nkn) :: HS2_TOT
    real(kind = DBL_PREC), dimension(nkn) :: H2S_DIVISOR
    real(kind = DBL_PREC), dimension(nkn) :: FRAC_H2S_IN_H2S_TOT
    real(kind = DBL_PREC), dimension(nkn) :: FRAC_HS_MINUS_IN_H2S_TOT
    real(kind = DBL_PREC), dimension(nkn) :: FRAC_S_MINUS_TWO_IN_H2S_TOT
    real(kind = DBL_PREC), dimension(nkn) :: H2S
    real(kind = DBL_PREC), dimension(nkn) :: HS_MINUS
    real(kind = DBL_PREC), dimension(nkn) :: S_MINUS_TWO
    real(kind = DBL_PREC), dimension(nkn) :: FE_II_DISS

    real(kind = DBL_PREC), dimension(nkn) :: FE_II_PART
    real(kind = DBL_PREC), dimension(nkn) :: FE_III_DISS
    real(kind = DBL_PREC), dimension(nkn) :: FE_III_PART
    real(kind = DBL_PREC), dimension(nkn) :: MN_II_DISS
    real(kind = DBL_PREC), dimension(nkn) :: MN_II_PART

    ! -------------------------------------------------------------------------
    ! New variables added (9 August 2016)
    ! -------------------------------------------------------------------------
    ! FE_II_DISS_EQ               : Dissolved Fe2+ in equilibrium (solubility of Fe2+)
    ! FE_III_DISS_EQ              : Dissolved Fe2+ in equilibrium (solubility of Fe3+)
    real(kind = DBL_PREC), dimension(nkn) :: FE_II_DISS_EQ
    real(kind = DBL_PREC), dimension(nkn) :: FE_III_DISS_EQ
    ! -------------------------------------------------------------------------
    ! End of new variables added (9 August 2016)
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! End of variables added 25 January 2016 for a simple version of equlibrium
    ! based aquatic chemistry calculations for Fe2+ and Fe3+
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! New auxillary variables introduced 27 January 2016
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DOC_MIN_DOXY
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DOC_MIN_NO3N
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DOC_MIN_MN_IV
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DOC_MIN_FE_III
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DOC_MIN_DOC
    real(kind = DBL_PREC), dimension(nkn) :: LIM_DOXY_RED
    real(kind = DBL_PREC), dimension(nkn) :: LIM_NO3N_RED
    real(kind = DBL_PREC), dimension(nkn) :: LIM_MN_IV_RED
    real(kind = DBL_PREC), dimension(nkn) :: LIM_FE_III_RED
    real(kind = DBL_PREC), dimension(nkn) :: LIM_S_PLUS_6_RED
    real(kind = DBL_PREC), dimension(nkn) :: LIM_DOC_RED
    real(kind = DBL_PREC), dimension(nkn) :: K_NO3_RED
    real(kind = DBL_PREC), dimension(nkn) :: K_MN_IV_RED
    real(kind = DBL_PREC), dimension(nkn) :: K_FE_III_RED
    real(kind = DBL_PREC), dimension(nkn) :: K_S_PLUS_6_RED
    real(kind = DBL_PREC), dimension(nkn) :: K_DOC_RED
    real(kind = DBL_PREC), dimension(nkn) :: MN_IV_DISS
    ! -------------------------------------------------------------------------
    ! End of new auxillary variables introduced 27 January 2016
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! New auxillary variables introduced 28 January 2016
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DON_MIN_DOXY
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DON_MIN_NO3N
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DON_MIN_MN_IV
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DON_MIN_FE_III
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DON_MIN_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DON_MIN_DOC

    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DOP_MIN_DOXY
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DOP_MIN_NO3N
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DOP_MIN_MN_IV
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DOP_MIN_FE_III
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DOP_MIN_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn) :: PH_CORR_DOP_MIN_DOC
    ! -------------------------------------------------------------------------
    ! End of new auxillary variables introduced 28 January 2016
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! New auxillary variables introduced 29 January 2016
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn) :: K_A_CH4
    real(kind = DBL_PREC), dimension(nkn) :: K_A_H2S
    real(kind = DBL_PREC), dimension(nkn) :: CH4_SAT
    real(kind = DBL_PREC), dimension(nkn) :: H2S_SAT
    ! -------------------------------------------------------------------------
    ! End of new auxillary variables introduced 29 January 2016
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! New auxillary variables introduced / January 2016
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn) :: PE
    ! -------------------------------------------------------------------------
    ! End of new auxillary variables introduced / January 2016
    ! -------------------------------------------------------------------------


    ! -------------------------------------------------------------------------
    ! New auxillary variables introduced 6 July 2016
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn) :: DIP_OVER_IP
    ! -------------------------------------------------------------------------
    ! End of new auxillary variables introduced 6 July 2016
    ! -------------------------------------------------------------------------
    integer     :: iron_oxidation

    ! -------------------------------------------------------------------------
    ! New auxillary variables introduced 9 August 2016
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn) :: DISS_FE_II_CONC_TS_END
    real(kind = DBL_PREC), dimension(nkn) :: DISS_FE_II_CONC_TS_AVG
    real(kind = DBL_PREC), dimension(nkn) :: DISS_FE_III_CONC_TS_END
    real(kind = DBL_PREC), dimension(nkn) :: DISS_FE_III_CONC_TS_AVG
    
    
    ! -------------------------------------------------------------------------
    ! New auxillary variables introduced 9 August 2016
    ! -------------------------------------------------------------------------
    integer  :: DO_ADVANCED_REDOX_SIMULATION
    
    !debug_stranger = .true. ! True if check for strange values
    debug_stranger = .false.
    
    CONSIDER_ALKALNITY_DERIVATIVE = 1
    CONSIDER_INORG_C_DERIVATIVE   = 1
    CONSIDER_CO2_REARATION        = 1

    ! indicator for light limitation  algorithm:
    ! 1 - modified Smith, 0 - Ali
    smith = 1

    ! indicator of iron oxidation formulation
    !iron_oxidation = 1  !complex, no calibration
    iron_oxidation = 0  !simple, with calibration parameter

    PROCESS_RATES(:,:,:) = 0.0D0
    DERIVATIVES  (:,:)   = 0.0D0
    error = 0
    DO_ADVANCED_REDOX_SIMULATION = 0

    if (present(ADVANCED_REDOX_OPTION)) then
        if (ADVANCED_REDOX_OPTION > 0) then
            DO_ADVANCED_REDOX_SIMULATION = 1
        end if
    end if

    ! calculates value of available DON for cyanobacteria
    ! variable name 'frac_avail_DON'
    call calc_frac_avail_DON

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
                        NODES_STRANGE_ext(j) = (node_active(k))
                        j=j+1
                    end if
                end do

                print *, '========================================='
                print *, 'PELAGIC_KINETICS:  Variable ',i !'Cell ',k
                print *, 'Initial state is NaN or Inf:'
                print *, 'NODE_NUMBERS int. =',NODES_STRANGE_int
                print *, 'NODE_NUMBERS ext. =',NODES_STRANGE_ext
                print *, 'VALUES            =',STRANGERS

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
    if(nstate.ne.30) then
        print *, 'PELAGIC_KINETICS: Number of state variables is wrong', nstate
        stop
    end if

    NH4_N           (:)      = STATE_VARIABLES(:,NH4_N_INDEX         )     ! AMMONIUM NITROGEN
    NO3_N           (:)      = STATE_VARIABLES(:,NO3_N_INDEX         )     ! NITRATE NITROGEN
    PO4_P           (:)      = STATE_VARIABLES(:,PO4_P_INDEX         )     ! ORTHOPHOSPHATE PHOSPHORUS
    DISS_OXYGEN     (:)      = STATE_VARIABLES(:,DISS_OXYGEN_INDEX   )     ! DISSOLVED OXYGEN
    DIA_C           (:)      = STATE_VARIABLES(:,DIA_C_INDEX         )     ! DIATOMS CARBON
    ZOO_C           (:)      = STATE_VARIABLES(:,ZOO_C_INDEX         )     ! ZOOPLANKTON CARBON
    ZOO_N           (:)      = STATE_VARIABLES(:,ZOO_N_INDEX         )     ! ZOOPLANKTON NITROGEN
    ZOO_P           (:)      = STATE_VARIABLES(:,ZOO_P_INDEX         )     ! ZOOPLANKTON PHOSPHORUS
    DET_PART_ORG_C  (:)      = STATE_VARIABLES(:,DET_PART_ORG_C_INDEX)     ! DETRITUS PARTICULATE ORG. CARBON
    DET_PART_ORG_N  (:)      = STATE_VARIABLES(:,DET_PART_ORG_N_INDEX)     ! DETRITUS PARTICULATE ORG. NITROGEN
    DET_PART_ORG_P  (:)      = STATE_VARIABLES(:,DET_PART_ORG_P_INDEX)     ! DETRITUS PARTICULATE ORG. PHOSPHORUS
    DISS_ORG_C      (:)      = STATE_VARIABLES(:,DISS_ORG_C_INDEX    )     ! DISSOLVED ORGANIC CARBON
    DISS_ORG_N      (:)      = STATE_VARIABLES(:,DISS_ORG_N_INDEX    )     ! DISSOLVED ORGANIC NITROGEN
    DISS_ORG_P      (:)      = STATE_VARIABLES(:,DISS_ORG_P_INDEX    )     ! DISSOLVED ORGANIC PHOSPHORUS
    CYN_C           (:)      = STATE_VARIABLES(:,CYN_C_INDEX         )     ! NON FIXING CYANOBACTERIA CARBON
    OPA_C           (:)      = STATE_VARIABLES(:,OPA_C_INDEX         )     ! OTHER PHYTOPLANKTON CARBON
    DISS_Si         (:)      = STATE_VARIABLES(:,DISS_Si_INDEX       )     ! DISSOLOVED SILICA
    PART_Si         (:)      = STATE_VARIABLES(:,PART_Si_INDEX       )     ! PARTICULATE SILICA
    FIX_CYN_C       (:)      = STATE_VARIABLES(:,FIX_CYN_C_INDEX     )     ! FIXING CYANOBACTERIA CARBON
    INORG_C         (:)      = STATE_VARIABLES(:,INORG_C_INDEX       )     ! INORG CARBON CARBON
    TOT_ALK         (:)      = STATE_VARIABLES(:,TOT_ALK_INDEX       )     ! TOTAL ALKALNITY
    FE_II           (:)      = STATE_VARIABLES(:,FE_II_INDEX         )     ! Iron charged as plus 2
    FE_III          (:)      = STATE_VARIABLES(:,FE_III_INDEX        )     ! Iron chargen as plus 3
    MN_II           (:)      = STATE_VARIABLES(:,MN_II_INDEX         )     ! Manganese charged as plus 2
    MN_IV           (:)      = STATE_VARIABLES(:,MN_IV_INDEX         )     ! Manganese charged as plus 4
    CA              (:)      = STATE_VARIABLES(:,CA_INDEX            )
    MG              (:)      = STATE_VARIABLES(:,MG_INDEX            )
    S_PLUS_6        (:)      = STATE_VARIABLES(:,S_PLUS_6_INDEX      )
    S_MINUS_2       (:)      = STATE_VARIABLES(:,S_MINUS_2_INDEX     )
    CH4_C           (:)      = STATE_VARIABLES(:,CH4_C_INDEX         )

    !INITIALIZE DRIVING_FUNCTIONS
    if(n_driving_functions.ne.10) then
        print *, 'PELAGIC_KINETICS: Number of elements in DRIVING_FUNCTIONS is not 10', &
                 n_driving_functions
        stop
    end if

    TEMP     (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 1)
    SALT     (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 2)

    ! Conversion from W/m^2 to langleys
    I_A      (1:nkn) =(DRIVING_FUNCTIONS(1:nkn, 3) * 5.0D-1 * 8.64D4 * 0.238846) / 1.0D4

    FDAY     (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 4)
    AIRTEMP  (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 5)
    WINDS    (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 6)
    ELEVATION(1:nkn) = DRIVING_FUNCTIONS(1:nkn, 7)
    DEPTH    (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 8)
    K_B_E    (1:nkn) = DRIVING_FUNCTIONS(1:nkn, 9)
    ice_cover(1:nkn) = DRIVING_FUNCTIONS(1:nkn,10)

    !INITIALIZE FLAGS
    if(nflags.ne.5) then
       print *, 'PELAGIC_KINETICS: Number of elements in FLAGS is wrong', nflags
       stop
    end if

    SAFE_MODE                  = FLAGS(1)
    SURFACE_BOX                = FLAGS(2)
    FIRST_TIME_STEP            = FLAGS(3)
    INIT_OPTION_OF_FE_II_DISS  = FLAGS(4)
    INIT_OPTION_OF_FE_III_DISS = FLAGS(5)

    !INITIALIZE MODEL CONSTANTS
!     if(nconst.ne.291) then
!         print *, 'PELAGIC_KINETICS: Number of elements in MODEL_CONSTANTS is wrong', nconst
!         stop
!     end if

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

        call CO2SYS(CO2SYS_PAR1         , CO2SYS_PAR2         , CO2SYS_PAR1TYPE , &
                    CO2SYS_PAR2TYPE     , CO2SYS_SALT         , CO2SYS_TEMPIN   , &
                    CO2SYS_TEMPOUT      , CO2SYS_PRESIN       , CO2SYS_PRESOUT  , &
                    CO2SYS_SI           , CO2SYS_PO4          , CO2SYS_pHSCALEIN, &
                    CO2SYS_K1K2CONSTANTS, CO2SYS_KSO4CONSTANTS, CO2SYS_OUT_DATA , &
                    CO2SYS_NICEHEADERS  , &
                    CO2SYS_ntps)

        pH         (1:nkn) = CO2SYS_OUT_DATA(1:nkn, 18)
        K_ONE_TIP  (1:nkn) = CO2SYS_OUT_DATA(1:nkn, 75)
        K_TWO_TIP  (1:nkn) = CO2SYS_OUT_DATA(1:nkn, 76)
        K_THREE_TIP(1:nkn) = CO2SYS_OUT_DATA(1:nkn, 77)
        H2CO3      (1:nkn) = CO2SYS_OUT_DATA(1:nkn, 23)
        H_PLUS     (1:nkn) = 10.0D0 ** (-CO2SYS_OUT_DATA(1:nkn,18))
        HCO3       (1:nkn) = CO2SYS_OUT_DATA(1:nkn, 21)
        CO3        (1:nkn) = CO2SYS_OUT_DATA(1:nkn, 22)

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

    HCO3 = HCO3 / 1000000.0
    CO3  = CO3  / 1000000.0
    ! -----------------------------------------------------------------
    ! Additions made 30th November 2015 related to dissolved and solid
    ! iron and mangese species, soon to be replaced by the
    ! aquatic chemistry model that will be called
    !
    !
    !        EXACTLY HERE
    !
    ! -----------------------------------------------------------------


    ! -------------------------------------------------------------------------
    ! New equlibrium calculations added at 25th of January 2016
    ! This is the first version of the aquatic chemistry model that was
    ! promissed 30th November 2015
    ! -------------------------------------------------------------------------

	! CO2SYS_OUT_DATA(:, 21), CO2SYS_OUT_DATA(:, 22)can not be used because deallocated

    if (DO_ADVANCED_REDOX_SIMULATION > 0) then
	    call REDOX_AND_SPECIATION &
            (DISS_OXYGEN, NO3_N, MN_IV, FE_III, S_PLUS_6, DISS_ORG_C, &
             S_MINUS_2 , MN_II, FE_II , HCO3, CO3, & 
             TEMP, SALT, PH, ELEVATION, &
             K_HS_DOXY_RED_LIM   , K_HS_NO3N_RED_LIM , K_HS_MN_IV_RED_LIM , &
             K_HS_FE_III_RED_LIM , K_HS_S_PLUS_6_RED_LIM, &
             K_HS_DOXY_RED_INHB  , K_HS_NO3N_RED_INHB, K_HS_MN_IV_RED_INHB, &
             K_HS_FE_III_RED_INHB, K_HS_S_PLUS_6_RED_INHB, nkn, &
             LIM_DOXY_RED        , LIM_NO3N_RED          , LIM_MN_IV_RED  , &
             LIM_FE_III_RED      , LIM_S_PLUS_6_RED      , LIM_DOC_RED, &
             PE, FE_II_DISS_EQ   , FE_III_DISS_EQ, MN_II_DISS)
        
        FE_III_DISS_EQ = FE_III_DISS_EQ * 56000.0D0
        PROCESS_RATES(1:nkn,FE_III_INDEX, 4) = FE_III_DISS_EQ
        ! ---------------------------------------------------------------------------
        ! Handle Fe2+ dissolution changes by Ali Ertürk, 2 nd of july 2016
        ! ---------------------------------------------------------------------------
        
        ! For now, take fixed values for equlibrium and solubility constants
        ! In future, change these to temperature and/or other environmental
        ! variable based equations.
        K_EQ_S_1  (:) = 8.90D-8
        K_EQ_S_2  (:) = 1.20D-13
        K_SP_FES  (:) = 8.00D-19
        !K_SP_FEOH3(:) = 1.00D-36
        !K_SP_FEPO4(:) = 9.91D-16
        
        ! -------------------------------------------------------------------------
        ! Calculate H2S Species
        ! -------------------------------------------------------------------------
        
        ! Convert sulphide to moles for aquatic chemistry calculations
        HS2_TOT(:) = S_MINUS_2(:) / 32000.0D0
        
        H2S_DIVISOR(:) = (H_PLUS(:) * H_PLUS(:)) + &
            (H_PLUS(:) * K_EQ_S_1(:)) + (K_EQ_S_1(:) * K_EQ_S_2(:))
        
        FRAC_HS_MINUS_IN_H2S_TOT   (:) = (H_PLUS(:)   * K_EQ_S_1(:)) / H2S_DIVISOR(:)
        FRAC_S_MINUS_TWO_IN_H2S_TOT(:) = (K_EQ_S_1(:) * K_EQ_S_2(:)) / H2S_DIVISOR(:)
        
        FRAC_H2S_IN_H2S_TOT(:)         = &
            1.0D0 - (FRAC_HS_MINUS_IN_H2S_TOT(:) + FRAC_S_MINUS_TWO_IN_H2S_TOT(:))
        
        H2S        (:) = HS2_TOT(:) * FRAC_H2S_IN_H2S_TOT        (:)
        HS_MINUS   (:) = HS2_TOT(:) * FRAC_HS_MINUS_IN_H2S_TOT   (:)
        S_MINUS_TWO(:) = HS2_TOT(:) * FRAC_S_MINUS_TWO_IN_H2S_TOT(:)
        ! -------------------------------------------------------------------------
        ! End of calculate H2S Species
        ! -------------------------------------------------------------------------
        
        ! -------------------------------------------------------------------------
        ! Changes by Ali Ertürk (2 July 2016)
        ! Change updated by Ali (9 August 2016)
        ! -------------------------------------------------------------------------
        call IRON_II_DISSOLUTION(HS2_TOT, PH, TOT_ALK, nkn, FE_II_DISS_EQ)
        FE_II_DISS_EQ = FE_II_DISS_EQ * 56000.0D0
        
        !write(2000,*) FE_II_DISS(557)
        ! -------------------------------------------------------------------------
        ! End of changes by Ali Ertürk (2 July 2016)
        ! -------------------------------------------------------------------------
        
        if(debug_stranger) then
            if (STRANGERSD(FE_II_DISS,VALUE_strange,nkn).eq.1) then
                print *, '========================================'
                nstrange = count(VALUE_strange)
                allocate(STRANGERS    (nstrange))
                allocate(NODES_STRANGE(nstrange))
        
                j=1
        
                do k=1,nkn
                    if(VALUE_strange(k)) then
                        STRANGERS    (j) = FE_II_DISS(k)
                        NODES_STRANGE(j) = node_active(k)
                        j=j+1
                    end if
                end do
        
                print *, 'PELAGIC KINETICS: '
                write(*,*) 'TIME          : ', TIME
                write(*,*) 'FE_II_DISS is NaN or Inf:'
                print *, 'NODE_NUMBERS=',NODES_STRANGE
                print *, 'VALUES=',STRANGERS
        
                write(*,*)
                write(*,*) 'Related variables'
                write(*,*) '-----------------'
                write(*,*) 'HS2_TOT      : ', (HS2_TOT(NODES_STRANGE(j)),j=1,nstrange)
                write(*,*) 'PH           : ', (PH     (NODES_STRANGE(j)),j=1,nstrange)
                write(*,*) 'TOT_ALK      : ', (TOT_ALK(NODES_STRANGE(j)),j=1,nstrange)
        
                deallocate(STRANGERS)
                deallocate(NODES_STRANGE)
                stop
            end if
        end if
        
        ! -------------------------------------------------------------------------
        ! Updated by Ali and Petras, 9 August 2016
        !
        ! If it is the first timestep, initialize the dissolved fraction by assuming
        ! that the system is in equlibrum (option 1) or the dissolved fractions will
        ! be initialized (option 2). Otherwise calculate the dissolved fractions will
        ! be calculated by simple submodel based on the analytical solution.
        ! -------------------------------------------------------------------------
        where (FE_II_DISS_EQ > FE_II)
            FE_II_DISS_EQ = FE_II
        end where
        
        if (FIRST_TIME_STEP > 0) then
        
            select case (INIT_OPTION_OF_FE_II_DISS)
        
                case(1)
                    ! Handle the different between saturated and nonsaturated cases.
                    ! If more Fe2+ is allowed to dissolve than the total Fe2+ present then
                    !
                    ! Assume unsaturated case, ie
                    !     - Dissolved Fe2+ is equal to total Fe2+
                    !     - Multiplier for dissolved Fe2+ is equal to 1
                    !     - Multiplier for particulate Fe2+ is equal to 0
                    ! otherwise
                    !     - Dissolved Fe2+ is smaller than to total Fe2+
                    !     - Multiplier for dissolved Fe2+ is between 0 and 1
                    !     - Multiplier for particulate Fe2+ is between 0 and 1
        
                    where(FE_II_DISS_EQ >= FE_II)
                        FE_II_DISS      = FE_II
                        MULT_FE_II_DISS = 1.0D0
                        MULT_FE_II_PART = 0.0D0
                    elsewhere
                        MULT_FE_II_DISS = FE_II_DISS_EQ / FE_II
                        MULT_FE_II_PART = 1.0D0 - MULT_FE_II_DISS
                    end where
        
                case(2)
                    ! Get the initial fraction of dissolved Fe2+ from model constants and
                    ! recalculate the initial concentration accordingly.
                    MULT_FE_II_DISS = INIT_MULT_FE_II_DISS
                    MULT_FE_II_PART = 1.0D0 - MULT_FE_II_DISS
        
            end select
        
            call CALC_DISS_ME_CONC &
                 (FE_II                                                    , & ! Total Fe2+
                  (MULT_FE_II_DISS * FE_II)                                , & ! Dissolved Fe2+ from previous time step
                  FE_II_DISS_EQ                                            , & ! Equilibrium concentration for dissolved Fe2+
                  (k_DISS_FE_II_20 * (THETA_k_DISS_FE_II**(TEMP - 20.0D0))), & ! Dissolution rate constant for Fe2+
                  TIME_STEP                                                , & ! Time step in days
                  nkn                                                      , &
                  1                                                        , & ! number of layers
                  DISS_FE_II_CONC_TS_END                                   , & ! Estimated dissolved Fe2+ at the end of time step (for output)
                  DISS_FE_II_CONC_TS_AVG)                                      ! Estimated avg. dissolved Fe2+ during the timestep to be used for kinetic calculations
        else
            call CALC_DISS_ME_CONC &
                 (FE_II                                                    , & ! Total Fe2+
                  (SAVED_OUTPUTS(:,1) * FE_II)                             , & ! Dissolved Fe2+ from previous time step
                  FE_II_DISS_EQ                                            , & ! Equilibrium concentration for dissolved Fe2+
                  (k_DISS_FE_II_20 * (THETA_k_DISS_FE_II**(TEMP - 20.0D0))), & ! Dissolution rate constant for Fe2+
                  TIME_STEP                                                , & ! Time step in days
                  nkn                                                      , &
                  1                                                        , & ! number of layers
                  DISS_FE_II_CONC_TS_END                                   , & ! Estimated dissolved Fe2+ at the end of time step (for output)
                  DISS_FE_II_CONC_TS_AVG)                                      ! Estimated avg. dissolved Fe2+ during the timestep to be used for kinetic calculations
        end if
        
        
        where(DISS_FE_II_CONC_TS_AVG >= FE_II)
            FE_II_DISS      = FE_II
            MULT_FE_II_DISS = 1.0D0
            MULT_FE_II_PART = 0.0D0
        elsewhere
            MULT_FE_II_DISS = DISS_FE_II_CONC_TS_AVG / FE_II
            MULT_FE_II_PART = 1.0D0 - MULT_FE_II_DISS
        end where
        
        where(DISS_FE_II_CONC_TS_END >= FE_II)
            DISS_FE_II_CONC_TS_END = FE_II
        end where
        
        ! -------------------------------------------------------------------------
        ! End of calculate the dissolved and particulate fractions of Fe2+
        ! -------------------------------------------------------------------------
        
        ! -------------------------------------------------------------------------
        ! Calculate the dissolved and particulate fractions of Fe3+
        ! -------------------------------------------------------------------------
        
        ! -------------------------------------------------------------------------
        ! Updated by Ali and Petras, 9 August 2016
        !
        ! If it is the first timestep, initialize the dissolved fraction by assuming
        ! that the system is in equlibrum (option 1) or the dissolved fractions will
        ! be initialized (option 2). Otherwise calculate the dissolved fractions will
        ! be calculated by simple submodel based on the analytical solution.
        ! -------------------------------------------------------------------------
        PROCESS_RATES(1:nkn,FE_III_INDEX, 5) = FE_III_DISS_EQ
        PROCESS_RATES(1:nkn,FE_III_INDEX, 6) = FE_III
        
        where (FE_III_DISS_EQ > FE_III)
            FE_III_DISS_EQ = FE_III
        end where
        
        PROCESS_RATES(1:nkn,FE_III_INDEX, 7) = FE_III_DISS_EQ
        
        if (FIRST_TIME_STEP > 0) then
        
            select case (INIT_OPTION_OF_FE_III_DISS)
        
                case(1)
                    ! Handle the different between saturated and nonsaturated cases.
                    ! If more Fe3+ is allowed to dissolve than the total Fe3+ present then
                    !
                    ! Assume unsaturated case, ie
                    !     - Dissolved Fe3+ is equal to total Fe3+
                    !     - Multiplier for dissolved Fe3+ is equal to 1
                    !     - Multiplier for particulate Fe3+ is equal to 0
                    ! otherwise
                    !     - Dissolved Fe3+ is smaller than to total Fe3+
                    !     - Multiplier for dissolved Fe3+ is between 0 and 1
                    !     - Multiplier for particulate Fe3+ is between 0 and 1
        
                    where(FE_III_DISS_EQ >= FE_III)
                        FE_III_DISS      = FE_III
                        MULT_FE_III_DISS = 1.0D0
                        MULT_FE_III_PART = 0.0D0
                    elsewhere
                        MULT_FE_III_DISS = FE_III_DISS_EQ / FE_III
                        MULT_FE_III_PART = 1.0D0 - MULT_FE_III_DISS
                    end where
        
                case(2)
                    ! Get the initial fraction of dissolved Fe3+ from model constants and
                    ! recalculate the initial concentration accordingly.
                    MULT_FE_III_DISS = INIT_MULT_FE_III_DISS
                    MULT_FE_III_PART = 1.0D0 - MULT_FE_III_DISS
        
                end select
        
                
            call CALC_DISS_ME_CONC &
                 (FE_III                                                     , & ! Total Fe3+
                  (MULT_FE_III_DISS * FE_III)                                , & ! Dissolved Fe3+ from previous time step
                  FE_III_DISS_EQ                                             , & ! Equilibrium concentration for dissolved Fe3+
                  (k_DISS_FE_III_20 * (THETA_k_DISS_FE_III**(TEMP - 20.0D0))), & ! Dissolution rate constant for Fe3+
                  TIME_STEP                                                  , & ! Time step in days
                  nkn                                                        , &
                  1                                                          , & ! number of layers
                  DISS_FE_III_CONC_TS_END                                    , & ! Estimated dissolved Fe3+ at the end of time step (for output)
                  DISS_FE_III_CONC_TS_AVG)                                       ! Estimated avg. dissolved Fe3+ during the timestep to be used for kinetic calculations
        else
            call CALC_DISS_ME_CONC &
                 (FE_III                                                     , & ! Total Fe3+
                  (SAVED_OUTPUTS(:,2) * FE_III)                              , & ! Dissolved Fe3+ from previous time step
                  FE_III_DISS_EQ                                             , & ! Equilibrium concentration for dissolved Fe3+
                  (k_DISS_FE_III_20 * (THETA_k_DISS_FE_III**(TEMP - 20.0D0))), & ! Dissolution rate constant for Fe3+
                  TIME_STEP                                                  , & ! Time step in days
                  nkn                                                        , &
                  1                                                          , & ! number of layers
                  DISS_FE_III_CONC_TS_END                                    , & ! Estimated dissolved Fe3+ at the end of time step (for output)
                  DISS_FE_III_CONC_TS_AVG)                                       ! Estimated avg. dissolved Fe3+ during the timestep to be used for kinetic calculations
        end if
        
        PROCESS_RATES(1:nkn,FE_III_INDEX, 8) = DISS_FE_III_CONC_TS_AVG
        PROCESS_RATES(1:nkn,FE_III_INDEX, 9) = DISS_FE_III_CONC_TS_END
        
            
        where(DISS_FE_III_CONC_TS_AVG >= FE_III)
            FE_III_DISS      = FE_III
            MULT_FE_III_DISS = 1.0D0
            MULT_FE_III_PART = 0.0D0
        elsewhere
            MULT_FE_III_DISS = DISS_FE_III_CONC_TS_AVG / FE_III
            MULT_FE_III_PART = 1.0D0 - MULT_FE_III_DISS
        end where
        
        where(DISS_FE_III_CONC_TS_END >= FE_III)
            DISS_FE_III_CONC_TS_END = FE_III
        end where
        
        
        !MULT_FE_III_DISS = DISS_FE_III_CONC_TS_AVG / FE_III
        !MULT_FE_III_PART = 1.0D0 - MULT_FE_III_DISS
        ! -------------------------------------------------------------------------
        ! End of calculate the dissolved and particulate fractions of Fe3+
        ! -------------------------------------------------------------------------
        
        ! Handle the different between saturated and nonsaturated cases.
        ! If more Mn2+ is allowed to dissolve than the total Mn2+ present then
        !
        ! Assume unsaturated case, ie
        !     - Dissolved Mn2+ is equal to total Mn2+
        !     - Multiplier for dissolved Mn2+ is equal to 1
        !     - Multiplier for particulate Mn2+ is equal to 0
        ! otherwise
        !     - Dissolved Mn2+ is smaller than to total Mn2+
        !     - Multiplier for dissolved Mn2+ is between 0 and 1
        !     - Multiplier for particulate Mn2+ is between 0 and 1
        
        where(MN_II_DISS >= MN_II)
            MN_II_DISS      = MN_II
            MULT_MN_II_DISS = 1.0D0
            MULT_MN_II_PART = 0.0D0
        elsewhere
            MULT_MN_II_DISS = MN_II_DISS / MN_II
            MULT_MN_II_PART = 1.0D0 - MULT_MN_II_DISS
        end where
        
        ! -------------------------------------------------------------------------
        ! New equlibrium calculations added at 25th of January 2016
        ! -------------------------------------------------------------------------
        
        MULT_MN_IV_DISS(:)  = 0.0D0
        MULT_MN_IV_PART(:)  = 1.0D0
        ! -----------------------------------------------------------------------
        ! End of additions made 30th November 2015 related to dissolved and solid
        ! iron and mangese species
        ! -----------------------------------------------------------------------
        
        ! -----------------------------------------------------------------------
        ! Introduced 26th of January 2016 to allow dissolved fractions of iron
        ! (Fe2+ and Fe3+) and manganese (Mn2+ and Mn4+) to outside to be used
        ! by AQUABC settling calculations
        !
        ! Updated by Ali and Petras 9th of August 2016
        ! -----------------------------------------------------------------------
        SAVED_OUTPUTS(:,1) = DISS_FE_II_CONC_TS_END  / FE_II
        SAVED_OUTPUTS(:,2) = DISS_FE_III_CONC_TS_END / FE_III
        SAVED_OUTPUTS(:,3) = MULT_MN_II_DISS(:)
        SAVED_OUTPUTS(:,4) = MULT_MN_IV_DISS(:)
        ! -----------------------------------------------------------------------
        ! End of additions 26th of January 2016
        ! -----------------------------------------------------------------------
        
        !Calculate the dissolved MN IV
        MN_IV_DISS(:) = MN_IV(:) * MULT_MN_IV_DISS(:)
        
        ! Calculate derived variables
        !     call derived_vars(nkn,pH,STATE_VARIABLES, nstate, &
        !                       MODEL_CONSTANTS, nconst,WC_OUTPUTS, noutput)
        
        call IP_SOLUBLE_FRACTION &
             (FE_III     , &
              PO4_P      , &
              K_ONE_TIP  , &
              K_TWO_TIP  , &
              K_THREE_TIP, &
              PH         , &
              nkn        , &
              1          , &
              DIP_OVER_IP)
        
        SAVED_OUTPUTS(:,5) = DIP_OVER_IP(:)
    else
        LIM_DOXY_RED = DISS_OXYGEN  / (DISS_OXYGEN + K_HS_DOXY_RED_LIM)

        LIM_NO3N_RED = (NO3_N / (NO3_N + K_HS_NO3N_RED_LIM)) * &
            (K_HS_DOXY_RED_INHB / (DISS_OXYGEN + K_HS_DOXY_RED_INHB))
    
        DIP_OVER_IP = 1.0D0
            
        SAVED_OUTPUTS(:,1) = 0.0D0
        SAVED_OUTPUTS(:,2) = 0.0D0
        SAVED_OUTPUTS(:,3) = 0.0D0
        SAVED_OUTPUTS(:,4) = 0.0D0
        SAVED_OUTPUTS(:,5) = 1.0D0
    end if
    !*****************************************
    !     D I S S O L V E D  O X Y G E N     !
    !*****************************************
    if (present(SURFACE_BOXES)) then
        do k=1,nkn
            DISS_OXYGEN_SAT(k) = DO_SATURATION(TEMP(k), SALT(k), ELEVATION(k))

            if (SURFACE_BOXES(k) == 1) then !first layer

                if (K_A < 0.0D0) then
                    K_A_CALC(k) = KAWIND(WINDS(k), TEMP(k), AIRTEMP(k), DEPTH(k), 3.0D0)
                    R_AERATION(k) = K_A_CALC(k) * (DISS_OXYGEN_SAT(k) - DISS_OXYGEN(k))
                else
                    K_A_CALC(k) = K_A
                
                    R_AERATION(k) = K_A_CALC(k) * (DISS_OXYGEN_SAT(k) - DISS_OXYGEN(k)) * &
                         (THETA_K_A ** (TEMP(k) - 2.0D1))
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
                R_AERATION(k) = 0.0D0 ! other layers
            end if
        end do
    else
        do k=1,nkn
            DISS_OXYGEN_SAT(k) = DO_SATURATION(TEMP(k), SALT(k), ELEVATION(k))

            if (SURFACE_BOX == 1) then !first layer

                if (K_A < 0.0D0) then
                    K_A_CALC(k) = KAWIND(WINDS(k), TEMP(k), AIRTEMP(k), DEPTH(k), 3.0D0)
                    R_AERATION(k) = K_A_CALC(k) * (DISS_OXYGEN_SAT(k) - DISS_OXYGEN(k))
                else
                    K_A_CALC(k) = K_A
                
                    R_AERATION(k) = K_A_CALC(k) * (DISS_OXYGEN_SAT(k) - DISS_OXYGEN(k)) * &
                         (THETA_K_A ** (TEMP(k) - 2.0D1))
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
                R_AERATION(k) = 0.0D0 ! other layers
            end if
        end do
    end if
    
    ! Calculate the total phytoplankton.
    PHYT_TOT_C = DIA_C + CYN_C + OPA_C + FIX_CYN_C

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
                 (PO4_P * DIP_OVER_IP)   , & ! Change, 6 July 2016, original call was PO4P
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
            DISS_ORG_N              , &
            (PO4_P * DIP_OVER_IP)   , & ! Change, 6 July 2016, original call was PO4P
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
            DISS_ORG_N                   , &
            (PO4_P * DIP_OVER_IP)        , & ! Change, 6 July 2016, original call was PO4P
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
            (PO4_P * DIP_OVER_IP)   , & ! Change, 6 July 2016, original call was PO4P
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

    if (present(ZOOP_OPTION_1)) then
        if (ZOOP_OPTION_1 > 0) then
            ACTUAL_ZOO_N_TO_C = ZOO_N / ZOO_C
            ACTUAL_ZOO_P_TO_C = ZOO_P / ZOO_C
        end if
    end if

    if(debug_stranger) then
        if (STRANGERSD(R_ZOO_GROWTH,VALUE_strange,nkn).eq.1) then
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
            write(*,*) 'R_ZOO_FEEDING_DET_PART_ORG_C     : ', R_ZOO_FEEDING_DET_PART_ORG_C
            write(*,*) '    ZOO_C                        : ', ZOO_C
            write(*,*) '    KG_ZOO_DIA                   : ', KG_ZOO_DIA
            write(*,*) '    KG_ZOO_CYN                   : ', KG_ZOO_CYN
            write(*,*) '    KG_ZOO_OPA                   : ', KG_ZOO_OPA
            write(*,*) '    FOOD_FACTOR_ZOO_OPA          : ', FOOD_FACTOR_ZOO_OPA
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
    LIM_P_DISS_DET_PART_ORG_P = KHS_DISS_P / (KHS_DISS_P + DIP_OVER_IP*PO4_P)
    LIM_PHY_P_DISS_DET_PART_ORG_P = LIM_P_DISS_DET_PART_ORG_P * FAC_PHYT_DET_PART_ORG_P * PHYT_TOT_C

    R_DET_PART_ORG_P_DISSOLUTION = (KDISS_DET_PART_ORG_P_20 + LIM_PHY_P_DISS_DET_PART_ORG_P) * &
       (THETA_KDISS_DET_PART_ORG_P ** (TEMP - 2.0D1)) * DET_PART_ORG_P

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


    !*********************************************************************!
    !     SILICON                                                         !
    !*********************************************************************!

    !Dissolution rate of biogenic silicon
    R_PART_Si_DISS = KDISS_PART_SI_20 * &
            (THETA_KDISS_PART_SI ** (TEMP - 2.0D1)) * PART_SI


    !*********************************************************************!
    !     MINERALIZATION OF DOC, DON, DOP whith bacteria are not modelled.
    !     Called abiotic in the sense of modelling method
    !*********************************************************************!

    !Algal dependent mineralisation rate
    
    if (DO_ADVANCED_REDOX_SIMULATION > 0) then
        call ORGANIC_CARBON_MINERALIZATION &
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
                 DISS_OXYGEN                 , &
                 NO3_N                       , &
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
    else
        LIM_PHYT_AMIN_DOC = FAC_PHYT_AMIN_DOC * PHYT_TOT_C

        call CALCULATE_PH_CORR &
             (PH_CORR_DOC_MIN_DOXY, PH, PH_MIN_DOC_MIN_DOXY, PH_MAX_DOC_MIN_DOXY, nkn)
        
        call CALCULATE_PH_CORR &
             (PH_CORR_DOC_MIN_NO3N, PH, PH_MIN_DOC_MIN_NO3N, PH_MAX_DOC_MIN_NO3N, nkn)

        R_ABIOTIC_DOC_MIN_DOXY = &
            (K_MIN_DOC_DOXY_20 + LIM_PHYT_AMIN_DOC) * &
            (THETA_K_MIN_DOC_DOXY ** (TEMP - 2.0D1)) * LIM_DOXY_RED * &
            PH_CORR_DOC_MIN_DOXY * &
            (DISS_ORG_C / (DISS_ORG_C + K_HS_DOC_MIN_DOXY)) * DISS_ORG_C

        R_ABIOTIC_DOC_MIN_NO3N = &
            K_MIN_DOC_NO3N_20  * (THETA_K_MIN_DOC_NO3N ** (TEMP - 2.0D1)) * &
            LIM_NO3N_RED * PH_CORR_DOC_MIN_NO3N * &
            (DISS_ORG_C / (DISS_ORG_C + K_HS_DOC_MIN_NO3N)) * DISS_ORG_C
     
        R_ABIOTIC_DOC_MIN_MN_IV    = 0.0D0
        R_ABIOTIC_DOC_MIN_FE_III   = 0.0D0
        R_ABIOTIC_DOC_MIN_S_PLUS_6 = 0.0D0
        R_ABIOTIC_DOC_MIN_DOC      = 0.0D0
    end if

    ! Accerelation of mineralisation when DIN is scarce
    LIM_N_AMIN_DON = KHS_AMIN_N / (KHS_AMIN_N + (NH4_N + NO3_N))
    LIM_PHY_N_AMIN_DON = LIM_N_AMIN_DON * FAC_PHYT_AMIN_DON * PHYT_TOT_C

    ! -------------------------------------------------------------------------
    ! DON mineralization compatible with redox cycle
    ! -------------------------------------------------------------------------

    ! 28 January 2016, the following commented lines are replaced in order to be
    ! compitable with the redox sequence
    ! call CALCULATE_PH_CORR(PH_CORR_DON_MIN, PH, PH_MIN_DON_MIN, PH_MIN_DON_MAX, nkn)
    !
    !R_ABIOTIC_DON_MIN = &
    !    (K_MIN_DON_20 + LIM_PHY_N_AMIN_DON) * (THETA_K_MIN_DON ** (TEMP - 2.0D1)) * &
    !    PH_CORR_DON_MIN * DISS_ORG_N * (1 - frac_avail_DON)

    !(1 - frac_avail_DON) counts fraction available for minerasilation bybacteria
    !frac_avail_DON - fraction available for cyanobacteria

    call CALCULATE_PH_CORR &
         (PH_CORR_DON_MIN_DOXY, PH, PH_MIN_DON_MIN_DOXY, PH_MAX_DON_MIN_DOXY, nkn)
    
    call CALCULATE_PH_CORR &
         (PH_CORR_DON_MIN_NO3N, PH, PH_MIN_DON_MIN_NO3N, PH_MAX_DON_MIN_NO3N, nkn)

    R_ABIOTIC_DON_MIN_DOXY = &
        (K_MIN_DON_DOXY_20 + LIM_PHY_N_AMIN_DON) * &
        (THETA_K_MIN_DON_DOXY ** (TEMP - 2.0D1)) * &
        LIM_DOXY_RED * PH_CORR_DON_MIN_DOXY * &
        (DISS_ORG_N / (DISS_ORG_N + K_HS_DON_MIN_DOXY)) * DISS_ORG_N * &
        (1.0D0 - frac_avail_DON)

    ! No phytoplankton or cyanobacteria when there is no oxygen so mineralization
    ! rate calculation differs
    R_ABIOTIC_DON_MIN_NO3N = &
        K_MIN_DON_NO3N_20  * (THETA_K_MIN_DON_NO3N ** (TEMP - 2.0D1)) * &
        LIM_NO3N_RED * PH_CORR_DON_MIN_NO3N * &
        (DISS_ORG_N / (DISS_ORG_N + K_HS_DON_MIN_NO3N)) * DISS_ORG_N

    R_ABIOTIC_DON_MIN_MN_IV    = 0.0D0
    R_ABIOTIC_DON_MIN_FE_III   = 0.0D0
    R_ABIOTIC_DON_MIN_S_PLUS_6 = 0.0D0        
    R_ABIOTIC_DON_MIN_DOC      = 0.0D0

    if (DO_ADVANCED_REDOX_SIMULATION > 0) then
        
        call CALCULATE_PH_CORR &
             (PH_CORR_DON_MIN_MN_IV   , PH, PH_MIN_DON_MIN_MN_IV   , &
              PH_MAX_DON_MIN_MN_IV   , nkn)
        
        call CALCULATE_PH_CORR &
             (PH_CORR_DON_MIN_FE_III  , PH, PH_MIN_DON_MIN_FE_III  , &
              PH_MAX_DON_MIN_FE_III  , nkn)
        
        call CALCULATE_PH_CORR &
             (PH_CORR_DON_MIN_S_PLUS_6, PH, PH_MIN_DON_MIN_S_PLUS_6, &
              PH_MAX_DON_MIN_S_PLUS_6, nkn)
        
        call CALCULATE_PH_CORR &
             (PH_CORR_DON_MIN_DOC     , PH, PH_MIN_DON_MIN_DOC     , &
              PH_MAX_DON_MIN_DOC     , nkn)
        
        R_ABIOTIC_DON_MIN_MN_IV   = &
            K_MIN_DON_MN_IV_20     * (THETA_K_MIN_DON_MN_IV    ** (TEMP - 2.0D1)) * &
            LIM_MN_IV_RED * PH_CORR_DON_MIN_MN_IV * &
            (DISS_ORG_N / (DISS_ORG_N + K_HS_DON_MIN_MN_IV)) * DISS_ORG_N
        
        R_ABIOTIC_DON_MIN_FE_III   = &
            K_MIN_DON_FE_III_20    * (THETA_K_MIN_DON_FE_III   ** (TEMP - 2.0D1)) * &
            LIM_FE_III_RED * PH_CORR_DON_MIN_FE_III * &
            (DISS_ORG_N / (DISS_ORG_N + K_HS_DON_MIN_FE_III)) * DISS_ORG_N
        
        R_ABIOTIC_DON_MIN_S_PLUS_6 = &
            K_MIN_DON_S_PLUS_6_20  * (THETA_K_MIN_DON_S_PLUS_6 ** (TEMP - 2.0D1)) * &
            LIM_S_PLUS_6_RED * PH_CORR_DON_MIN_S_PLUS_6 * &
            (DISS_ORG_N / (DISS_ORG_N + K_HS_DON_MIN_S_PLUS_6)) * DISS_ORG_N
        
        R_ABIOTIC_DON_MIN_DOC      = &
            (K_MIN_DON_DOC_20      * (THETA_K_MIN_DON_DOC      ** (TEMP - 2.0D1)) * &
            LIM_DOC_RED * PH_CORR_DON_MIN_DOC * &
            (DISS_ORG_N / (DISS_ORG_N + K_HS_DON_MIN_DOC)) * DISS_ORG_N)
    end if


    ! -------------------------------------------------------------------------
    ! End of DON mineralization compatible with redox cycle
    ! -------------------------------------------------------------------------

    ! Accerelation of mineralisation when DIP is scarce
    LIM_P_AMIN_DOP = KHS_AMIN_P / (KHS_AMIN_P + DIP_OVER_IP*PO4_P)
    LIM_PHY_P_AMIN_DOP = LIM_P_AMIN_DOP * FAC_PHYT_AMIN_DOP * PHYT_TOT_C

    ! -------------------------------------------------------------------------
    ! DOP mineralization compatible with redox cycle
    ! -------------------------------------------------------------------------

    call CALCULATE_PH_CORR(PH_CORR_DOP_MIN_DOXY    , PH, PH_MIN_DOP_MIN_DOXY    , PH_MAX_DOP_MIN_DOXY    , nkn)
    call CALCULATE_PH_CORR(PH_CORR_DOP_MIN_NO3N    , PH, PH_MIN_DOP_MIN_NO3N    , PH_MAX_DOP_MIN_NO3N    , nkn)

    R_ABIOTIC_DOP_MIN_DOXY = &
        (K_MIN_DOP_DOXY_20 + LIM_PHY_P_AMIN_DOP) * (THETA_K_MIN_DOP_DOXY ** (TEMP - 2.0D1)) * &
        LIM_DOXY_RED * PH_CORR_DOP_MIN_DOXY * (DISS_ORG_P / (DISS_ORG_P + K_HS_DOP_MIN_DOXY)) * &
        DISS_ORG_P

    ! No phytoplankton or cyanobacteria when there is no oxygen so mineralization
    ! rate calculation differs
    R_ABIOTIC_DOP_MIN_NO3N = &
        K_MIN_DOP_NO3N_20  * (THETA_K_MIN_DOP_NO3N ** (TEMP - 2.0D1)) * &
        LIM_NO3N_RED * PH_CORR_DOP_MIN_NO3N * (DISS_ORG_P / (DISS_ORG_P + K_HS_DOP_MIN_NO3N)) * &
        DISS_ORG_P

    R_ABIOTIC_DOP_MIN_MN_IV    = 0.0D0
    R_ABIOTIC_DOP_MIN_FE_III   = 0.0D0
    R_ABIOTIC_DOP_MIN_S_PLUS_6 = 0.0D0
    R_ABIOTIC_DOP_MIN_DOC      = 0.0D0
        
    if (DO_ADVANCED_REDOX_SIMULATION > 0) then
        call CALCULATE_PH_CORR(PH_CORR_DOP_MIN_MN_IV   , PH, PH_MIN_DOP_MIN_MN_IV   , PH_MAX_DOP_MIN_MN_IV   , nkn)
        call CALCULATE_PH_CORR(PH_CORR_DOP_MIN_FE_III  , PH, PH_MIN_DOP_MIN_FE_III  , PH_MAX_DOP_MIN_FE_III  , nkn)
        call CALCULATE_PH_CORR(PH_CORR_DOP_MIN_S_PLUS_6, PH, PH_MIN_DOP_MIN_S_PLUS_6, PH_MAX_DOP_MIN_S_PLUS_6, nkn)
        call CALCULATE_PH_CORR(PH_CORR_DOP_MIN_DOC     , PH, PH_MIN_DOP_MIN_DOC     , PH_MAX_DOP_MIN_DOC     , nkn)
        
        R_ABIOTIC_DOP_MIN_MN_IV = &
            K_MIN_DOP_MN_IV_20  * (THETA_K_MIN_DOP_MN_IV ** (TEMP - 2.0D1)) * &
            LIM_MN_IV_RED * PH_CORR_DOP_MIN_MN_IV * (DISS_ORG_P / (DISS_ORG_P + K_HS_DOP_MIN_MN_IV)) * &
            DISS_ORG_P
        
        R_ABIOTIC_DOP_MIN_FE_III = &
            K_MIN_DOP_FE_III_20  * (THETA_K_MIN_DOP_FE_III ** (TEMP - 2.0D1)) * &
            LIM_FE_III_RED * PH_CORR_DOP_MIN_FE_III * (DISS_ORG_P / (DISS_ORG_P + K_HS_DOP_MIN_FE_III)) * &
            DISS_ORG_P
        
        R_ABIOTIC_DOP_MIN_S_PLUS_6 = &
            K_MIN_DOP_S_PLUS_6_20  * (THETA_K_MIN_DOP_S_PLUS_6 ** (TEMP - 2.0D1)) * &
            LIM_S_PLUS_6_RED * PH_CORR_DOP_MIN_S_PLUS_6 * (DISS_ORG_P / (DISS_ORG_P + K_HS_DOP_MIN_S_PLUS_6)) * &
            DISS_ORG_P
        
        R_ABIOTIC_DOP_MIN_DOC = &
            (K_MIN_DOP_DOC_20  * (THETA_K_MIN_DOP_DOC ** (TEMP - 2.0D1)) * &
             LIM_DOC_RED * PH_CORR_DOP_MIN_DOC * (DISS_ORG_P / (DISS_ORG_P + K_HS_DOP_MIN_DOC)) * DISS_ORG_P)
    end if
    ! -------------------------------------------------------------------------
    ! End of DOP mineralization compatible with redox cycle
    ! -------------------------------------------------------------------------

    !*******************************************************************************************!
    !     Nitrification of ammonia by bacteria are not modelled.
    !     Called abiotic in the sense of modelling method
    !*******************************************************************************************!
    LIM_NITR_OXY = DISS_OXYGEN / (KHS_NITR_OXY + DISS_OXYGEN)
    LIM_NITR_NH4_N = NH4_N / (KHS_NITR_NH4_N + NH4_N)

    call CALCULATE_PH_CORR(PH_CORR_NITR_NH4, PH, PH_NITR_NH4_MIN, PH_NITR_NH4_MAX, nkn)

    R_ABIOTIC_NITR = K_NITR_20 * LIM_NITR_OXY * LIM_NITR_NH4_N * &
                     PH_CORR_NITR_NH4 * (THETA_K_NITR ** (TEMP - 2.0D1)) * NH4_N

    ! -------------------------------------------------------------------------
    ! DENITRIFICATION
    ! -------------------------------------------------------------------------

    ! Introduced 28 January 2016 by Ali
    !
    ! [CH2O] + 4/5[NO3N-] -----> 1/2[N2] + [CO2]
    !
    ! Therefore, for each gram of DOC, that is mineralized over denitrification
    ! process,
    !
    !  - 14/(12 * 1.25) = 0.93 grams of NO3N is converted to N2
    !  - 1 gram of carbondioxide carbon is produced

    R_DENITRIFICATION = 0.93D0 * R_ABIOTIC_DOC_MIN_NO3N
    ! -------------------------------------------------------------------------
    ! END OF DENITRIFICATION
    ! -------------------------------------------------------------------------


    ! -------------------------------------------------------------------------
    ! MANGANESE REDUCTION
    ! -------------------------------------------------------------------------

    ! Introduced 29 January 2016 by Ali
    !
    ! [CH2O] + 2[MN_IV] -----> [CO2] + 2[MN_II]
    !
    ! Therefore, for each gram of DOC, that is mineralized over denitrification
    ! process,
    !
    !  - (2*52)/12 = 8.66 grams of manganese IV is reduced to manganese II
    !  - 1 gram of carbondioxide carbon is produced

    R_MN_IV_REDUCTION = 8.66D0 * R_ABIOTIC_DOC_MIN_MN_IV
    ! -------------------------------------------------------------------------
    ! END OF MANGANESE REDUCTION
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! IRON REDUCTION
    ! -------------------------------------------------------------------------

    ! Introduced 29 January 2016 by Ali
    !
    ! [CH2O] + 4[FE_III] -----> [CO2] + 4[FE_II]
    !
    ! Therefore, for each gram of DOC, that is mineralized over denitrification
    ! process,
    !
    !  - (4*56)/12 = 18.66 grams of iron III is reduced to iron II
    !  - 1 gram of carbondioxide carbon is produced

    R_FE_III_REDUCTION = 18.66D0 * R_ABIOTIC_DOC_MIN_FE_III
    ! -------------------------------------------------------------------------
    ! END OF IRON REDUCTION
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! SULPHATE REDUCTION
    ! -------------------------------------------------------------------------

    ! Introduced 28 January 2016 by Ali
    !
    ! [CH2O] + 1/2[SO4--] -----> [CO2] + 1/2[S--]
    !
    ! Therefore, for each gram of DOC, that is mineralized over denitrification
    ! process,
    !
    !  - (32/2)/12 = 1.33 grams of S_PLUS_6 is converted to S_MINUS_2
    !  - 1 gram of carbondioxide is produced

    R_SULPHATE_REDUCTION = 1.33D0 * R_ABIOTIC_DOC_MIN_S_PLUS_6
    ! -------------------------------------------------------------------------
    ! END OF SULPHATE REDUCTION
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! METHANOGENESIS
    ! -------------------------------------------------------------------------

    ! Introduced 29 January 2016 by Ali
    !
    ! [CH2O]  -----> 1/2[CH4] + 1/2[CO2]
    !
    ! Therefore, for each gram of DOC, that is mineralized over denitrification
    ! process,
    !
    !  - 0.5 grams of methane carbon is produced from DON
    !  - 0.5 grams of carbondioxide carbon is produced

    R_METHANOGENESIS = 0.5D0 * R_ABIOTIC_DOC_MIN_DOC
    ! -------------------------------------------------------------------------
    ! END OF METHANOGENESIS
    ! -------------------------------------------------------------------------

    !*********************************************************************!
    !     VOLATILIZATION OF UNIONIZED AMMONI.
    !*********************************************************************!
    call AMMONIA_VOLATILIZATION(R_AMMONIA_VOLATIL, NH4_N, pH, TEMP, K_A_CALC,nkn)

    !----------------------------------------------------------------------
    ! 2 February 2015
    ! New code added to account the effect of ice cover.
    !----------------------------------------------------------------------
    R_AMMONIA_VOLATIL = R_AMMONIA_VOLATIL * (1.0D0 - ice_cover)
    
    where (R_AMMONIA_VOLATIL < 0.0D0)
        R_AMMONIA_VOLATIL = 0.0D0
    end where
    !do k=1,nkn
    !    if (R_AMMONIA_VOLATIL(k) < 0.0D0) then
    !        write(*,*) 'Unexpected value for ammonia volatilization at node', k 
    !        write(*,*) 'R_AMMONIA_VOLATIL : ', R_AMMONIA_VOLATIL(k)
    !        write(*,*) 'ice_cover         : ', ice_cover        (k)
    !        write(*,*) 'NH4_N             : ', NH4_N            (k)
    !        write(*,*) 'pH                : ', pH               (k)
    !        write(*,*) 'TEMP              : ', TEMP             (k)
    !        write(*,*) 'K_A_CALC          : ', K_A_CALC         (k)
    !        stop
    !    end if
    !end do
    
    !----------------------------------------------------------------------
    ! End of new code added to account the effect of ice cover.
    !----------------------------------------------------------------------

    ! ---------------------------------------------------------------------
    ! Changes by Ali Ertürk, 6 th of July 2016
    !
    ! Following lines are commented
    ! ---------------------------------------------------------------------

    if (DO_ADVANCED_REDOX_SIMULATION > 0) then
        ! New kinetic rate calculations added 9 September 2015
        ! For now, no temparature corrections. Effect on temperature and other
        ! environmental conditions may be included after more detailed investigations
        
        ! Updated in 25th January 2016 where only dissolved fractions are allowed to be oxidized or reduced
        
        ! Iron
        if(iron_oxidation .eq. 0) then
            !simple formulation, k_OX_FE_II is calibration parameter
            where (DISS_OXYGEN < 1)
                R_FE_II_OXIDATION  = k_OX_FE_II * DISS_OXYGEN * (10.0D0 ** (PH - 7.0D0)) * FE_II
            elsewhere
                R_FE_II_OXIDATION  = k_OX_FE_II * (10.0D0 ** (PH - 7.0D0)) * FE_II
            end where
        end if
        
        if(iron_oxidation .eq. 1) then
            ! In the future include several options for heavy metal oxidation and possibly reduction
            ! Morgen and Lahav (2007) formulation. No calibration
            call IRON_II_OXIDATION(FE_II_DISS, DISS_OXYGEN, PH, TEMP, SALT, ELEVATION, nkn, R_FE_II_OXIDATION)
            ! ---------------------------------------------------------------------
            ! End of changes by Ali Ertürk, 6 th of July 2016
            ! ---------------------------------------------------------------------
        end if
        ! 29 January 2016
        ! Following commented lines are replaced by the new redox sequence based DOC mineralization
        ! it is the next visit of Ali and the redox sequences as described
        ! by Van Chappen and Wang 2015 and Katsev papers are now included.
        
        ! Manganese
        where (DISS_OXYGEN < 1)
            R_MN_II_OXIDATION  = k_OX_MN_II * DISS_OXYGEN * (10.0D0 ** (PH - 7.0D0))* MN_II
        elsewhere
            R_MN_II_OXIDATION  = k_OX_MN_II * (10.0D0 ** (PH - 7.0D0)) * MN_II
        end where
        
        ! 29 January 2016
        ! Following commented lines are replaced by the new redox sequence based DOC mineralization
        ! it is the next visit of Ali and the redox sequences as described
        ! by Van Chappen and Wang 2015 and Katsev papers are now included.
        
        ! End of new kinetic rate calculations added 9 September 2015
        
        ! -------------------------------------------------------------------------
        ! 29 January 2016 KINETICS OF NEW STATE VARIABLES
        ! -------------------------------------------------------------------------
        K_A_CH4 = K_A_CALC * 1.188D0
        K_A_H2S = K_A_CALC * 0.984D0
        CH4_SAT = 0.0D0 ! Assume that no methane is present in the atmosphere
        H2S_SAT = 0.0D0 ! Assume that no H2S is present in the atmosphere
        
        CH4_ATM_EXCHANGE = K_A_CH4 * (CH4_SAT - CH4_C)
        H2S_ATM_EXCHANGE = K_A_H2S * (H2S_SAT - (H2S * 32000.D0))
        
        R_METHANE_OXIDATION = &
            k_OX_CH4 * (THETA_k_OX_CH4 ** (TEMP - 20.0D0)) * CH4_C * &
            (DISS_OXYGEN / (k_HS_OX_CH4_DOXY + DISS_OXYGEN))
        
        R_SULPHIDE_OXIDATION = &
            k_OX_H2S * (THETA_k_OX_H2S ** (TEMP - 20.0D0)) * S_MINUS_2 * &
            (DISS_OXYGEN / (k_HS_OX_H2S_DOXY + DISS_OXYGEN))
    else
        R_FE_II_OXIDATION    = 0.0D0
        R_MN_II_OXIDATION    = 0.0D0
        R_SULPHIDE_OXIDATION = 0.0D0
        R_METHANE_OXIDATION  = 0.0D0
    end if
    ! -------------------------------------------------------------------------
    ! END OF KINETICS OF NEW STATE VARIABLES
    ! -------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------------
    !      Final calculation of derivatives
    !------------------------------------------------------------------------------------------------


    !AMMONIA NITROGEN
    PROCESS_RATES(1:nkn,NH4_N_INDEX, 1) = R_DIA_TOT_RESP     * DIA_N_TO_C
    PROCESS_RATES(1:nkn,NH4_N_INDEX, 2) = R_CYN_TOT_RESP     * CYN_N_TO_C
    PROCESS_RATES(1:nkn,NH4_N_INDEX, 3) = R_OPA_TOT_RESP     * OPA_N_TO_C
    PROCESS_RATES(1:nkn,NH4_N_INDEX, 4) = R_FIX_CYN_TOT_RESP * FIX_CYN_N_TO_C
    PROCESS_RATES(1:nkn,NH4_N_INDEX, 5) = R_ZOO_TOT_RESP     * ACTUAL_ZOO_N_TO_C
    
    PROCESS_RATES(1:nkn,NH4_N_INDEX, 6) = R_DIA_GROWTH       * PREF_NH4N_DIA * DIA_N_TO_C

    PROCESS_RATES(1:nkn,NH4_N_INDEX, 7) = R_CYN_GROWTH * PREF_NH4N_DON_CYN * CYN_N_TO_C * &
        (NH4_N / ((NH4_N + (DISS_ORG_N *  frac_avail_DON))))

    PROCESS_RATES(1:nkn,NH4_N_INDEX, 8) = R_OPA_GROWTH * PREF_NH4N_OPA * OPA_N_TO_C

    PROCESS_RATES(1:nkn,NH4_N_INDEX, 9) = &
        R_NON_FIX_CYN_GROWTH * PREF_NH4N_DON_FIX_CYN * FIX_CYN_N_TO_C * &
        (NH4_N / ((NH4_N + (DISS_ORG_N *  frac_avail_DON))))

    PROCESS_RATES(1:nkn,NH4_N_INDEX, 10) = R_ABIOTIC_NITR

    PROCESS_RATES(1:nkn,NH4_N_INDEX, 11) = &
        R_ABIOTIC_DON_MIN_DOXY   + R_ABIOTIC_DON_MIN_NO3N     + R_ABIOTIC_DON_MIN_MN_IV + &
        R_ABIOTIC_DON_MIN_FE_III + R_ABIOTIC_DON_MIN_S_PLUS_6 + R_ABIOTIC_DON_MIN_DOC

    PROCESS_RATES(1:nkn,NH4_N_INDEX, 12) = R_AMMONIA_VOLATIL
    
    ! Auxiliary
    PROCESS_RATES(1:nkn,NH4_N_INDEX, 13) = PREF_NH4N_DIA
    PROCESS_RATES(1:nkn,NH4_N_INDEX, 14) = PREF_NH4N_DON_CYN
    PROCESS_RATES(1:nkn,NH4_N_INDEX, 15) = PREF_NH4N_OPA
    PROCESS_RATES(1:nkn,NH4_N_INDEX, 16) = PREF_NH4N_DON_FIX_CYN


    DERIVATIVES(1:nkn,NH4_N_INDEX) = &
        PROCESS_RATES(1:nkn,NH4_N_INDEX, 1)  + PROCESS_RATES(1:nkn,NH4_N_INDEX, 2)  + &
        PROCESS_RATES(1:nkn,NH4_N_INDEX, 3)  + PROCESS_RATES(1:nkn,NH4_N_INDEX, 4)  + &
        PROCESS_RATES(1:nkn,NH4_N_INDEX, 5)  - PROCESS_RATES(1:nkn,NH4_N_INDEX, 6)  - &
        PROCESS_RATES(1:nkn,NH4_N_INDEX, 7)  - PROCESS_RATES(1:nkn,NH4_N_INDEX, 8)  - &
        PROCESS_RATES(1:nkn,NH4_N_INDEX, 9)  - PROCESS_RATES(1:nkn,NH4_N_INDEX, 10) + &
        PROCESS_RATES(1:nkn,NH4_N_INDEX, 11) - PROCESS_RATES(1:nkn,NH4_N_INDEX, 12)

    ! Code to debug NH4N
    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,NH4_N_INDEX),VALUE_strange,nkn).eq.1) then
            nstrange = count(VALUE_strange)
            allocate(STRANGERS        (nstrange))
            allocate(NODES_STRANGE    (nstrange))
            allocate(NODES_STRANGE_int(nstrange))
            allocate(NODES_STRANGE_ext(nstrange))

            j=1
            do k=1,nkn
                if(VALUE_strange(k)) then
                    STRANGERS        (j) = DERIVATIVES(k,NH4_N_INDEX)
                    NODES_STRANGE    (j) = node_active(k)
                    NODES_STRANGE_int(j) = node_active(k)
                    NODES_STRANGE_ext(j) = (node_active(k))
                    j=j+1
                end if
            end do

            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            print *, 'PELAGIC_KINETICS:'
            write(*,*) 'TIME   : ', TIME

            write(*,*) 'DERIVATIVES(NH4_N_INDEX) is not a number or infinite:'
            print *,   'NODE_NUMBERS      =', NODES_STRANGE
            print *,   'NODE_NUMBERS int. =', NODES_STRANGE_int
            print *,   'NODE_NUMBERS ext. =', NODES_STRANGE_ext
            print *,   'VALUES            =', STRANGERS

            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'
            write(*,*) 'NH4_N                        : ',  (NH4_N                    &
                                                           (NODES_STRANGE(j)),j=1,nstrange)

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
            
            write(*,*) K_B_E
            write(*,*) NODES_STRANGE
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
    PROCESS_RATES(1:nkn,NO3_N_INDEX, 1) = R_ABIOTIC_NITR
    PROCESS_RATES(1:nkn,NO3_N_INDEX, 2) = R_DENITRIFICATION
    PROCESS_RATES(1:nkn,NO3_N_INDEX, 3) = R_DIA_GROWTH         * (1.0D0 - PREF_NH4N_DIA)         * DIA_N_TO_C
    PROCESS_RATES(1:nkn,NO3_N_INDEX, 4) = R_CYN_GROWTH         * (1.0D0 - PREF_NH4N_DON_CYN)     * CYN_N_TO_C
    PROCESS_RATES(1:nkn,NO3_N_INDEX, 5) = R_OPA_GROWTH         * (1.0D0 - PREF_NH4N_OPA)         * OPA_N_TO_C
    PROCESS_RATES(1:nkn,NO3_N_INDEX, 6) = R_NON_FIX_CYN_GROWTH * (1.0D0 - PREF_NH4N_DON_FIX_CYN) * FIX_CYN_N_TO_C

    ! Auxiliary
    PROCESS_RATES(1:nkn,NO3_N_INDEX, 7) = PREF_NH4N_DIA
    PROCESS_RATES(1:nkn,NO3_N_INDEX, 8) = PREF_NH4N_DON_CYN
    PROCESS_RATES(1:nkn,NO3_N_INDEX, 9) = PREF_NH4N_OPA


    DERIVATIVES(1:nkn,NO3_N_INDEX) = &
        PROCESS_RATES(1:nkn,NO3_N_INDEX, 1) - PROCESS_RATES(1:nkn,NO3_N_INDEX, 2) - &
        PROCESS_RATES(1:nkn,NO3_N_INDEX, 3) - PROCESS_RATES(1:nkn,NO3_N_INDEX, 4) - &
        PROCESS_RATES(1:nkn,NO3_N_INDEX, 5) - PROCESS_RATES(1:nkn,NO3_N_INDEX, 6)


    !PHOSPHATE PHOSPHORUS
    PROCESS_RATES(1:nkn,PO4_P_INDEX, 1) = R_DIA_TOT_RESP     * DIA_P_TO_C
    PROCESS_RATES(1:nkn,PO4_P_INDEX, 2) = R_CYN_TOT_RESP     * CYN_P_TO_C
    PROCESS_RATES(1:nkn,PO4_P_INDEX, 3) = R_OPA_TOT_RESP     * OPA_P_TO_C
    PROCESS_RATES(1:nkn,PO4_P_INDEX, 4) = R_FIX_CYN_TOT_RESP * FIX_CYN_P_TO_C
    PROCESS_RATES(1:nkn,PO4_P_INDEX, 5) = R_ZOO_TOT_RESP     * ACTUAL_ZOO_P_TO_C
    PROCESS_RATES(1:nkn,PO4_P_INDEX, 6) = R_DIA_GROWTH       * DIA_P_TO_C
    PROCESS_RATES(1:nkn,PO4_P_INDEX, 7) = R_CYN_GROWTH       * CYN_P_TO_C
    PROCESS_RATES(1:nkn,PO4_P_INDEX, 8) = R_OPA_GROWTH       * OPA_P_TO_C
    PROCESS_RATES(1:nkn,PO4_P_INDEX, 9) = R_FIX_CYN_GROWTH   * FIX_CYN_P_TO_C

    PROCESS_RATES(1:nkn,PO4_P_INDEX, 10) = &
        R_ABIOTIC_DOP_MIN_DOXY   + R_ABIOTIC_DOP_MIN_NO3N     + R_ABIOTIC_DOP_MIN_MN_IV + &
        R_ABIOTIC_DOP_MIN_FE_III + R_ABIOTIC_DOP_MIN_S_PLUS_6 + R_ABIOTIC_DOP_MIN_DOC
    
    ! Auxiliary
    PROCESS_RATES(1:nkn,PO4_P_INDEX, 11) = TEMP
    PROCESS_RATES(1:nkn,PO4_P_INDEX, 12) = DISS_ORG_P


    DERIVATIVES(1:nkn,PO4_P_INDEX) = &
        PROCESS_RATES(1:nkn,PO4_P_INDEX, 1) + PROCESS_RATES(1:nkn,PO4_P_INDEX, 2)  + &
        PROCESS_RATES(1:nkn,PO4_P_INDEX, 3) + PROCESS_RATES(1:nkn,PO4_P_INDEX, 4)  + &
        PROCESS_RATES(1:nkn,PO4_P_INDEX, 5) - PROCESS_RATES(1:nkn,PO4_P_INDEX, 6) - &
        PROCESS_RATES(1:nkn,PO4_P_INDEX, 7) - PROCESS_RATES(1:nkn,PO4_P_INDEX, 8) - &
        PROCESS_RATES(1:nkn,PO4_P_INDEX, 9) + PROCESS_RATES(1:nkn,PO4_P_INDEX, 10)

    ! Debug for PO4P
    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,PO4_P_INDEX),VALUE_strange,nkn).eq.1) then

        nstrange = count(VALUE_strange)
        allocate(STRANGERS    (nstrange))
        allocate(NODES_STRANGE(nstrange))

        j=1
            do k=1,nkn
                if(VALUE_strange(k)) then
                    STRANGERS    (j) = DERIVATIVES(k,PO4_P_INDEX)
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
            write(*,*) '    TEMP            : ', (TEMP(NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) '    DISS_ORG_P      : ', (DISS_ORG_P(NODES_STRANGE(j)),j=1,nstrange)

            deallocate(STRANGERS)
            deallocate(NODES_STRANGE)
            stop
        end if
    end if

    !DISSOLVED SILICA SILICON
    PROCESS_RATES(1:nkn,DISS_Si_INDEX, 1) = R_PART_Si_DISS
    !PROCESS_RATES(1:nkn,DISS_Si_INDEX, 2) = R_DIA_TOT_RESP * DIA_SI_TO_C  !Should it come from respiration? Silica is in shell
    PROCESS_RATES(1:nkn,DISS_Si_INDEX, 2) = R_DIA_GROWTH   * DIA_SI_TO_C

    DERIVATIVES(1:nkn,DISS_Si_INDEX) = &
        PROCESS_RATES(1:nkn,DISS_Si_INDEX, 1) - PROCESS_RATES(1:nkn,DISS_Si_INDEX, 2)

    !DISSOLVED OXYGEN
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 1)  = R_AERATION
    ! formulation from EFDC
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 2)  = R_DIA_GROWTH       * (1.3D0 - 0.3D0*PREF_NH4N_DIA    )    * DIA_O2_TO_C
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 3)  = R_CYN_GROWTH       * (1.3D0 - 0.3D0*PREF_NH4N_DON_CYN)    * CYN_O2_TO_C
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 4)  = R_OPA_GROWTH       * (1.3D0 - 0.3D0*PREF_NH4N_OPA    )    * OPA_O2_TO_C
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 5)  = R_FIX_CYN_GROWTH   * (1.3D0 - 0.3D0*PREF_NH4N_DON_FIX_CYN)* FIX_CYN_O2_TO_C
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 6)  = R_DIA_TOT_RESP     * DIA_O2_TO_C
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 7)  = R_CYN_TOT_RESP     * CYN_O2_TO_C
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 8)  = R_OPA_TOT_RESP     * OPA_O2_TO_C
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 9)  = R_FIX_CYN_TOT_RESP * FIX_CYN_O2_TO_C
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 10) = R_ZOO_TOT_RESP     * ZOO_O2_TO_C
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 11) = R_ABIOTIC_NITR     * 4.57D0
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 12) = 2.66D0 * R_ABIOTIC_DOC_MIN_DOXY
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 13) = 0.43D0 * R_FE_II_OXIDATION
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 14) = 0.88D0 * R_MN_II_OXIDATION
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 15) = 2.00D0 * R_SULPHIDE_OXIDATION
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 16) = 2.66D0 * R_METHANE_OXIDATION
    

    ! Auxiliary
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 17) = K_A_CALC
    PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 18) = DISS_OXYGEN_SAT

    DERIVATIVES(1:nkn,DISS_OXYGEN_INDEX) = &
        PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 1)  + PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 2)  + &
        PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 3)  + PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 4)  + &
        PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 5)  - PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 6)  - &
        PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 7)  - PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 8)  - &
        PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 9)  - PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 10) - &
        PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 11) - PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 12) - &
        PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 13) - PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 14) - &
        PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 15) - PROCESS_RATES(1:nkn,DISS_OXYGEN_INDEX, 16)

    !Debug code for dissolved oxygen
    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,DISS_OXYGEN_INDEX),VALUE_strange,nkn).eq.1) then
            error=1

            nstrange = count(VALUE_strange)
            allocate(STRANGERS    (nstrange))
            allocate(NODES_STRANGE(nstrange))
            allocate(NODES_STRANGE_int(nstrange))
            allocate(NODES_STRANGE_ext(nstrange))

            j=1

            do k=1,nkn
                if(VALUE_strange(k)) then
                    STRANGERS    (j) = DERIVATIVES(k,DISS_OXYGEN_INDEX)
                    NODES_STRANGE(j) = k
                    NODES_STRANGE_int(j) = node_active(k)
                    NODES_STRANGE_ext(j) = (node_active(k))
                    j=j+1
                end if
            end do

            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            write(*,*) 'TIME   : ', TIME
            print *, 'PELAGIC_KINETICS:  Derivative(DISS_OXYGEN_INDEX) '
            print *, 'is NaN or Inf:'
            print *, 'NODE_NUMBERS=',NODES_STRANGE
            print *, 'NODE_NUMBERS int.=',NODES_STRANGE_int
            print *, 'NODE_NUMBERS ext.=',NODES_STRANGE_ext
            print *, 'VALUES=',STRANGERS

            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'

            write(*,*) '    ZOO_O2_TO_C          : ', ZOO_O2_TO_C
            write(*,*) '    DIA_O2_TO_C          : ', DIA_O2_TO_C
            write(*,*) '    CYN_O2_TO_C          : ', CYN_O2_TO_C
            write(*,*) '    OPA_O2_TO_C          : ', OPA_O2_TO_C
            print *,'SURFACE_BOX:', SURFACE_BOX
            print *,'FLAGS(2)=',FLAGS(2)

            write(*,*) '    K_A parameter        : ',  K_A

            write(*,*) ''
            write(*,*) 'R_AERATION               : ', (R_AERATION      &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'DISS_OXYGEN              : ', (DISS_OXYGEN     &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'DISS_OXYGEN_SAT          : ', (DISS_OXYGEN_SAT &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'K_A_CALC                 : ', (K_A_CALC        &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'R_DIA_GROWTH             : ', (R_DIA_GROWTH &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'R_CYN_GROWTH             : ', (R_CYN_GROWTH &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'R_OPA_GROWTH             : ', (R_OPA_GROWTH &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
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

            write(*,*) ''
            write(*,*) 'R_FE_II_OXIDATION        : ', (R_FE_II_OXIDATION        &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'R_MN_II_OXIDATION        : ', (R_MN_II_OXIDATION        &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'R_METHANE_OXIDATION      : ', (R_METHANE_OXIDATION      &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) ''
            write(*,*) 'R_SULPHIDE_OXIDATION     : ', (R_SULPHIDE_OXIDATION     &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            print *,'TEMP     =', (TEMP       (NODES_STRANGE(j)),j=1,nstrange)
            print *,'SALT     =', (SALT       (NODES_STRANGE(j)),j=1,nstrange)
            print *,'AIRTEMP  =', (AIRTEMP    (NODES_STRANGE(j)),j=1,nstrange)
            print *,'WINDS    =', (WINDS      (NODES_STRANGE(j)),j=1,nstrange)
            print *,'ELEVATION=', (ELEVATION  (NODES_STRANGE(j)),j=1,nstrange)
            print *,'DEPTH    =', (DEPTH      (NODES_STRANGE(j)),j=1,nstrange)

            stop
        end if
    end if

    !DIATOMS CARBON
    PROCESS_RATES(1:nkn,DIA_C_INDEX, 1)  = R_DIA_GROWTH
    PROCESS_RATES(1:nkn,DIA_C_INDEX, 2)  = R_DIA_TOT_RESP
	PROCESS_RATES(1:nkn,DIA_C_INDEX, 3)  = R_DIA_EXCR
    PROCESS_RATES(1:nkn,DIA_C_INDEX, 4)  = R_DIA_DEATH
    PROCESS_RATES(1:nkn,DIA_C_INDEX, 5)  = R_ZOO_FEEDING_DIA
    ! Auxiliary
    PROCESS_RATES(1:nkn,DIA_C_INDEX, 6)  = LIM_KG_DIA_TEMP
    PROCESS_RATES(1:nkn,DIA_C_INDEX, 7)  = LIM_KG_DIA_DOXY
    PROCESS_RATES(1:nkn,DIA_C_INDEX, 8)  = LIM_KG_DIA_N
    PROCESS_RATES(1:nkn,DIA_C_INDEX, 9)  = LIM_KG_DIA_P
    PROCESS_RATES(1:nkn,DIA_C_INDEX, 10) = LIM_KG_DIA_DISS_Si
    PROCESS_RATES(1:nkn,DIA_C_INDEX, 11) = LIM_KG_DIA_LIGHT
    PROCESS_RATES(1:nkn,DIA_C_INDEX, 12) = DIA_LIGHT_SAT

    DERIVATIVES(1:nkn,DIA_C_INDEX) = &
        PROCESS_RATES(1:nkn,DIA_C_INDEX, 1) - PROCESS_RATES(1:nkn,DIA_C_INDEX, 2) - &
        PROCESS_RATES(1:nkn,DIA_C_INDEX, 3) - PROCESS_RATES(1:nkn,DIA_C_INDEX, 4) - &
		PROCESS_RATES(1:nkn,DIA_C_INDEX, 5)

    !NON-NITROGEN FIXING CYANOBACTERIA CARBON
    PROCESS_RATES(1:nkn,CYN_C_INDEX, 1)  = R_CYN_GROWTH
    PROCESS_RATES(1:nkn,CYN_C_INDEX, 2)  = R_CYN_TOT_RESP
	PROCESS_RATES(1:nkn,CYN_C_INDEX, 3)  = R_CYN_EXCR
    PROCESS_RATES(1:nkn,CYN_C_INDEX, 4)  = R_CYN_DEATH
    PROCESS_RATES(1:nkn,CYN_C_INDEX, 5)  = R_ZOO_FEEDING_CYN
    ! Auxiliary
    PROCESS_RATES(1:nkn,CYN_C_INDEX, 6)  = LIM_KG_CYN_TEMP
    PROCESS_RATES(1:nkn,CYN_C_INDEX, 7)  = LIM_KG_CYN_DOXY
    PROCESS_RATES(1:nkn,CYN_C_INDEX, 8)  = LIM_KG_CYN_N
    PROCESS_RATES(1:nkn,CYN_C_INDEX, 9)  = LIM_KG_CYN_P
    PROCESS_RATES(1:nkn,CYN_C_INDEX, 10) = LIM_KG_CYN_LIGHT
    PROCESS_RATES(1:nkn,CYN_C_INDEX, 11) = I_A             !light langlays
    PROCESS_RATES(1:nkn,CYN_C_INDEX, 12) = CYN_LIGHT_SAT
    PROCESS_RATES(1:nkn,CYN_C_INDEX, 13) = TEMP

    DERIVATIVES(1:nkn,CYN_C_INDEX) = &
        PROCESS_RATES(1:nkn,CYN_C_INDEX, 1) - PROCESS_RATES(1:nkn,CYN_C_INDEX, 2) - &
        PROCESS_RATES(1:nkn,CYN_C_INDEX, 3) - PROCESS_RATES(1:nkn,CYN_C_INDEX, 4) - &
		PROCESS_RATES(1:nkn,CYN_C_INDEX, 5)

    ! Debug code for CYN_C
    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,CYN_C_INDEX),VALUE_strange,nkn).eq.1) then
            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            print *, 'PELAGIC_KINETICS:'
            write(*,*) 'TIME   : ', TIME
            write(*,*) 'DERIVATIVES(CYN_C_INDEX) is nor a number or infinite.'
            write(*,*) DERIVATIVES(1:nkn,CYN_C_INDEX)
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
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 1)  = R_FIX_CYN_GROWTH
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 2)  = R_FIX_CYN_TOT_RESP
	PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 3)  = R_FIX_CYN_EXCR
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 4)  = R_FIX_CYN_DEATH
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 5)  = R_ZOO_FEEDING_FIX_CYN
    
    ! Auxiliary
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 6)  = R_NON_FIX_CYN_GROWTH
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 7)  = R_FIX_FIX_CYN_GROWTH
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 8)  = R_FIX_FIX_CYN_GROWTH * FIX_CYN_N_TO_C !Nitrogen fixation
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 9)  = LIM_KG_FIX_CYN_TEMP
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 10) = LIM_KG_FIX_CYN_DOXY
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 11) = LIM_KG_FIX_FIX_CYN_N
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 12) = LIM_KG_FIX_FIX_CYN_P
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 13) = LIM_KG_NON_FIX_CYN_N
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 14) = LIM_KG_NON_FIX_CYN_P
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 15) = LIM_KG_FIX_CYN_LIGHT
    
	PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 16) = &
	    ((NH4_N + NO3_N)/14.D0)/(DIP_OVER_IP*PO4_P/31.D0) ! NP molar ratio
    
	PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 17) = NH4_N + NO3_N
    PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 18) = FIX_CYN_LIGHT_SAT

    DERIVATIVES(1:nkn,FIX_CYN_C_INDEX) = &
        PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 1) - PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 2) - &
        PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 3) - PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 4) - &
		PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 5)

    !OTHER PLANKTONIC ALGAE CARBON
    PROCESS_RATES(1:nkn,OPA_C_INDEX, 1) = R_OPA_GROWTH
    PROCESS_RATES(1:nkn,OPA_C_INDEX, 2) = R_OPA_TOT_RESP
	PROCESS_RATES(1:nkn,OPA_C_INDEX, 3) = R_OPA_EXCR
    PROCESS_RATES(1:nkn,OPA_C_INDEX, 4) = R_OPA_DEATH
    PROCESS_RATES(1:nkn,OPA_C_INDEX, 5) = R_ZOO_FEEDING_OPA
	
    ! Auxiliary
    PROCESS_RATES(1:nkn,OPA_C_INDEX, 6 ) = LIM_KG_OPA_TEMP
    PROCESS_RATES(1:nkn,OPA_C_INDEX, 7 ) = LIM_KG_OPA_DOXY
    PROCESS_RATES(1:nkn,OPA_C_INDEX, 8 ) = LIM_KG_OPA_N
    PROCESS_RATES(1:nkn,OPA_C_INDEX, 9 ) = LIM_KG_OPA_P
    PROCESS_RATES(1:nkn,OPA_C_INDEX, 10) = LIM_KG_OPA_LIGHT
    PROCESS_RATES(1:nkn,OPA_C_INDEX, 11) = OPA_LIGHT_SAT

    DERIVATIVES(1:nkn,OPA_C_INDEX) = &
        PROCESS_RATES(1:nkn,OPA_C_INDEX, 1) - PROCESS_RATES(1:nkn,OPA_C_INDEX, 2) - &
        PROCESS_RATES(1:nkn,OPA_C_INDEX, 3) - PROCESS_RATES(1:nkn,OPA_C_INDEX, 4) - &
		PROCESS_RATES(1:nkn,OPA_C_INDEX, 5)

    !ZOOPLANKTON CARBON
    PROCESS_RATES(1:nkn,ZOO_C_INDEX, 1) = R_ZOO_GROWTH
    PROCESS_RATES(1:nkn,ZOO_C_INDEX, 2) = 0.0
	
	if (present(ZOOP_OPTION_1)) then
        if (ZOOP_OPTION_1 > 0) then
            PROCESS_RATES(1:nkn,ZOO_C_INDEX, 2) = R_ZOO_EX_DOC
	    end if
	end if
	
    PROCESS_RATES(1:nkn,ZOO_C_INDEX, 3) = R_ZOO_TOT_RESP
    PROCESS_RATES(1:nkn,ZOO_C_INDEX, 4) = R_ZOO_DEATH

    PROCESS_RATES(1:nkn,ZOO_C_INDEX, 5) = R_ZOO_FEEDING_DIA
    PROCESS_RATES(1:nkn,ZOO_C_INDEX, 6) = R_ZOO_FEEDING_CYN
    PROCESS_RATES(1:nkn,ZOO_C_INDEX, 7) = R_ZOO_FEEDING_OPA
    PROCESS_RATES(1:nkn,ZOO_C_INDEX, 8) = R_ZOO_FEEDING_FIX_CYN
    PROCESS_RATES(1:nkn,ZOO_C_INDEX, 9) = R_ZOO_FEEDING_DET_PART_ORG_C


    DERIVATIVES(1:nkn,ZOO_C_INDEX) = &
        PROCESS_RATES(1:nkn,ZOO_C_INDEX, 1) - PROCESS_RATES(1:nkn,ZOO_C_INDEX, 2) - &
        PROCESS_RATES(1:nkn,ZOO_C_INDEX, 3) - PROCESS_RATES(1:nkn,ZOO_C_INDEX, 4)

    ! Debug code for ZOO_C
    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,ZOO_C_INDEX),VALUE_strange,nkn).eq.1) then
            print *, 'PELAGIC_KINETICS:'
            write(*,*) 'TIME   : ', TIME
            write(*,*) 'DERIVATIVES(ZOO_C_INDEX) is not a number or infinite'
            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'
            write(*,*) 'ZOO_C                                  : ', ZOO_C
            write(*,*) 'R_ZOO_GROWTH                           : ', R_ZOO_GROWTH
            write(*,*) '    R_ZOO_FEEDING_DIA                  : ', R_ZOO_FEEDING_DIA
            write(*,*) '    R_ZOO_FEEDING_CYN                  : ', R_ZOO_FEEDING_CYN
            write(*,*) '    R_ZOO_FEEDING_OPA                  : ', R_ZOO_FEEDING_OPA
            write(*,*) '    R_ZOO_FEEDING_FIX_CYN              : ', R_ZOO_FEEDING_FIX_CYN
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
    DERIVATIVES(1:nkn,ZOO_N_INDEX) = DERIVATIVES(1:nkn,ZOO_C_INDEX) * ACTUAL_ZOO_N_TO_C

    !ZOOPLANKTON PHOSPHORUS
    DERIVATIVES(1:nkn,ZOO_P_INDEX) = DERIVATIVES(1:nkn,ZOO_C_INDEX) * ACTUAL_ZOO_P_TO_C
    
    if (present(ZOOP_OPTION_1)) then
        if (ZOOP_OPTION_1 > 0) then

            PROCESS_RATES(1:nkn,ZOO_N_INDEX, 1) = R_ZOO_FEEDING_DIA            * DIA_N_TO_C
            PROCESS_RATES(1:nkn,ZOO_N_INDEX, 2) = R_ZOO_FEEDING_CYN            * CYN_N_TO_C
            PROCESS_RATES(1:nkn,ZOO_N_INDEX, 3) = R_ZOO_FEEDING_OPA            * OPA_N_TO_C
            PROCESS_RATES(1:nkn,ZOO_N_INDEX, 4) = R_ZOO_FEEDING_FIX_CYN        * FIX_CYN_N_TO_C
            PROCESS_RATES(1:nkn,ZOO_N_INDEX, 5) = R_ZOO_FEEDING_DET_PART_ORG_C * ACTUAL_DET_N_TO_C
            PROCESS_RATES(1:nkn,ZOO_N_INDEX, 6) = R_ZOO_EX_DON
            PROCESS_RATES(1:nkn,ZOO_N_INDEX, 7) = R_ZOO_TOT_RESP               * ACTUAL_ZOO_N_TO_C
            PROCESS_RATES(1:nkn,ZOO_N_INDEX, 8) = R_ZOO_DEATH                  * ACTUAL_ZOO_N_TO_C
            
            PROCESS_RATES(1:nkn,ZOO_N_INDEX, 9) = ACTUAL_ZOO_N_TO_C
    
            PROCESS_RATES(1:nkn,ZOO_P_INDEX, 1) = R_ZOO_FEEDING_DIA            * DIA_P_TO_C
            PROCESS_RATES(1:nkn,ZOO_P_INDEX, 2) = R_ZOO_FEEDING_CYN            * CYN_P_TO_C
            PROCESS_RATES(1:nkn,ZOO_P_INDEX, 3) = R_ZOO_FEEDING_OPA            * OPA_P_TO_C
            PROCESS_RATES(1:nkn,ZOO_P_INDEX, 4) = R_ZOO_FEEDING_FIX_CYN        * FIX_CYN_P_TO_C
            PROCESS_RATES(1:nkn,ZOO_P_INDEX, 5) = R_ZOO_FEEDING_DET_PART_ORG_C * ACTUAL_DET_P_TO_C
            PROCESS_RATES(1:nkn,ZOO_P_INDEX, 6) = R_ZOO_EX_DOP
            PROCESS_RATES(1:nkn,ZOO_P_INDEX, 7) = R_ZOO_TOT_RESP               * ACTUAL_ZOO_P_TO_C
            PROCESS_RATES(1:nkn,ZOO_P_INDEX, 8) = R_ZOO_DEATH                  * ACTUAL_ZOO_P_TO_C
            
            PROCESS_RATES(1:nkn,ZOO_P_INDEX, 9) = ACTUAL_ZOO_P_TO_C    
    
            DERIVATIVES(1:nkn,ZOO_N_INDEX) = &
                PROCESS_RATES(1:nkn,ZOO_N_INDEX, 1)  + PROCESS_RATES(1:nkn,ZOO_N_INDEX, 2)  + &
                PROCESS_RATES(1:nkn,ZOO_N_INDEX, 3)  + PROCESS_RATES(1:nkn,ZOO_N_INDEX, 4)  + &
                PROCESS_RATES(1:nkn,ZOO_N_INDEX, 5)  - PROCESS_RATES(1:nkn,ZOO_N_INDEX, 6)  - &
                PROCESS_RATES(1:nkn,ZOO_N_INDEX, 7)  - PROCESS_RATES(1:nkn,ZOO_N_INDEX, 8)

            DERIVATIVES(1:nkn,ZOO_P_INDEX) = &
                PROCESS_RATES(1:nkn,ZOO_P_INDEX, 1)  + PROCESS_RATES(1:nkn,ZOO_P_INDEX, 2)  + &
                PROCESS_RATES(1:nkn,ZOO_P_INDEX, 3)  + PROCESS_RATES(1:nkn,ZOO_P_INDEX, 4)  + &
                PROCESS_RATES(1:nkn,ZOO_P_INDEX, 5)  - PROCESS_RATES(1:nkn,ZOO_P_INDEX, 6)  - &
                PROCESS_RATES(1:nkn,ZOO_P_INDEX, 7)  - PROCESS_RATES(1:nkn,ZOO_P_INDEX, 8)
        end if                  
    end if
    
    !DEAD ORGANIC CARBON PARTICLES
    PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 1) = R_DIA_DEATH
    PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 2) = R_CYN_DEATH
    PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 3) = R_OPA_DEATH
    PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 4) = R_FIX_CYN_DEATH
    PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 5) = R_ZOO_DEATH
    PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 6) = R_ZOO_FEEDING_DET_PART_ORG_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 7) = R_DET_PART_ORG_C_DISSOLUTION

    DERIVATIVES(1:nkn,DET_PART_ORG_C_INDEX) = &
        PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 1) + PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 2) + &
        PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 3) + PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 4) + &
        PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 5) - PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 6) - &
        PROCESS_RATES(1:nkn,DET_PART_ORG_C_INDEX, 7)

    ! Debug code for DET_PART_ORG_C
    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,DET_PART_ORG_C_INDEX),VALUE_strange,nkn).eq.1) then
            print *, 'PELAGIC_KINETICS:'
            write(*,*) 'TIME   : ', TIME
            write(*,*) 'DERIVATIVES(DET_PART_ORG_C_INDEX) is not a number or infinite.'
            write(*,*) DERIVATIVES(1:nkn,DET_PART_ORG_C_INDEX)
            write(*,*)
            write(*,*) 'Related variables'
            write(*,*) '-----------------'
            write(*,*) 'R_DIA_DEATH                        : ', R_DIA_DEATH
            write(*,*) 'R_CYN_DEATH                        : ', R_CYN_DEATH
            write(*,*) 'R_OPA_DEATH                        : ', R_OPA_DEATH
            write(*,*) 'R_FIX_CYN_DEATH                    : ', R_FIX_CYN_DEATH
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
    PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 1) = R_DIA_DEATH                  * DIA_N_TO_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 2) = R_CYN_DEATH                  * CYN_N_TO_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 3) = R_OPA_DEATH                  * OPA_N_TO_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 4) = R_FIX_CYN_DEATH              * FIX_CYN_N_TO_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 5) = R_ZOO_DEATH                  * ACTUAL_ZOO_N_TO_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 6) = R_ZOO_FEEDING_DET_PART_ORG_C * ACTUAL_DET_N_TO_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 7) = R_DET_PART_ORG_N_DISSOLUTION
    ! Auxiliary
    PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 8) = ACTUAL_DET_N_TO_C

    DERIVATIVES(1:nkn,DET_PART_ORG_N_INDEX) = &
        PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 1) + PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 2) + &
        PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 3) + PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 4) + &
        PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 5) - PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 6) - &
        PROCESS_RATES(1:nkn,DET_PART_ORG_N_INDEX, 7)

    !DEAD ORGANIC PHOSPHORUS PARTICLES
    PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 1) = R_DIA_DEATH                  * DIA_P_TO_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 2) = R_CYN_DEATH                  * CYN_P_TO_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 3) = R_OPA_DEATH                  * OPA_P_TO_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 4) = R_FIX_CYN_DEATH              * FIX_CYN_P_TO_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 5) = R_ZOO_DEATH                  * ACTUAL_ZOO_P_TO_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 6) = R_ZOO_FEEDING_DET_PART_ORG_C * ACTUAL_DET_P_TO_C
    PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 7) = R_DET_PART_ORG_P_DISSOLUTION
    ! Auxiliary
    PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 8) = ACTUAL_DET_P_TO_C

    DERIVATIVES(1:nkn,DET_PART_ORG_P_INDEX) = &
        PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 1) + PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 2) + &
        PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 3) + PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 4) + &
        PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 5) - PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 6) - &
        PROCESS_RATES(1:nkn,DET_PART_ORG_P_INDEX, 7)

    !PARTICULATE SILICA
    PROCESS_RATES(1:nkn,PART_Si_INDEX, 1) = R_DIA_DEATH
    PROCESS_RATES(1:nkn,PART_Si_INDEX, 2) = R_ZOO_FEEDING_DIA * DIA_Si_TO_C
    PROCESS_RATES(1:nkn,PART_Si_INDEX, 3) = R_PART_Si_DISS

    DERIVATIVES(1:nkn,PART_Si_INDEX) = &
        PROCESS_RATES(1:nkn,PART_Si_INDEX, 1) + PROCESS_RATES(1:nkn,PART_Si_INDEX, 2) - &
        PROCESS_RATES(1:nkn,PART_Si_INDEX, 3)

    !DISSOLVED ORGANIC CARBON
    PROCESS_RATES(1:nkn,DISS_ORG_C_INDEX, 1) = R_DET_PART_ORG_C_DISSOLUTION
    PROCESS_RATES(1:nkn,DISS_ORG_C_INDEX, 2) = 0.   !No excretion for a while

    PROCESS_RATES(1:nkn,DISS_ORG_C_INDEX, 3)  = &
        R_ABIOTIC_DOC_MIN_DOXY   + R_ABIOTIC_DOC_MIN_NO3N     + R_ABIOTIC_DOC_MIN_MN_IV + &
        R_ABIOTIC_DOC_MIN_FE_III + R_ABIOTIC_DOC_MIN_S_PLUS_6 + R_ABIOTIC_DOC_MIN_DOC

    PROCESS_RATES(1:nkn,DISS_ORG_C_INDEX, 4) = R_DIA_EXCR + R_CYN_EXCR + R_FIX_CYN_EXCR + R_OPA_EXCR

    DERIVATIVES(1:nkn,DISS_ORG_C_INDEX) = &
        PROCESS_RATES(1:nkn,DISS_ORG_C_INDEX, 1) + PROCESS_RATES(1:nkn,DISS_ORG_C_INDEX, 2) - &
        PROCESS_RATES(1:nkn,DISS_ORG_C_INDEX, 3) + PROCESS_RATES(1:nkn,DISS_ORG_C_INDEX, 4)

    !DISSOLVED ORGANIC NITROGEN
    PROCESS_RATES(1:nkn,DISS_ORG_N_INDEX, 1) = R_DET_PART_ORG_N_DISSOLUTION
    PROCESS_RATES(1:nkn,DISS_ORG_N_INDEX, 2) = 0. ! R_ZOO_EX_DON No excretion for a while

    PROCESS_RATES(1:nkn,DISS_ORG_N_INDEX, 3)  = &
        R_ABIOTIC_DON_MIN_DOXY   + R_ABIOTIC_DON_MIN_NO3N     + R_ABIOTIC_DON_MIN_MN_IV + &
        R_ABIOTIC_DON_MIN_FE_III + R_ABIOTIC_DON_MIN_S_PLUS_6 + R_ABIOTIC_DON_MIN_DOC

    PROCESS_RATES(1:nkn,DISS_ORG_N_INDEX, 4) = &
        (R_DIA_EXCR * DIA_N_TO_C)         + (R_CYN_EXCR * CYN_N_TO_C) + &
        (R_FIX_CYN_EXCR * FIX_CYN_N_TO_C) + (R_OPA_EXCR * OPA_N_TO_C)

    PROCESS_RATES(1:nkn,DISS_ORG_N_INDEX, 5) = &
        (R_CYN_GROWTH * PREF_NH4N_DON_CYN * CYN_N_TO_C) * &
        ((DISS_ORG_N * frac_avail_DON) / ((NH4_N + (DISS_ORG_N *  frac_avail_DON))))

    PROCESS_RATES(1:nkn,DISS_ORG_N_INDEX, 6) = &
        R_NON_FIX_CYN_GROWTH * PREF_NH4N_DON_FIX_CYN * CYN_N_TO_C * &
         ((DISS_ORG_N * frac_avail_DON)  / ((NH4_N + (DISS_ORG_N *  frac_avail_DON))))

    DERIVATIVES(1:nkn,DISS_ORG_N_INDEX) = &
        PROCESS_RATES(1:nkn,DISS_ORG_N_INDEX, 1) + PROCESS_RATES(1:nkn,DISS_ORG_N_INDEX, 2) - &
        PROCESS_RATES(1:nkn,DISS_ORG_N_INDEX, 3) + PROCESS_RATES(1:nkn,DISS_ORG_N_INDEX, 4) - &
        PROCESS_RATES(1:nkn,DISS_ORG_N_INDEX, 5) - PROCESS_RATES(1:nkn,DISS_ORG_N_INDEX, 6)

     ! Debug code for  DISS_ORG_N
     if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,DISS_ORG_N_INDEX),VALUE_strange,nkn).eq.1) then
            nstrange = count(VALUE_strange)
            allocate(STRANGERS    (nstrange))
            allocate(NODES_STRANGE(nstrange))
            j=1

            do k=1,nkn
                if(VALUE_strange(k)) then
                    STRANGERS    (j) = DERIVATIVES(k,DISS_ORG_N_INDEX)
                    NODES_STRANGE(j) = node_active(k)
                    j=j+1
                end if
            end do
            print *, 'PELAGIC_KINETICS:  Derivative 16 is strange'
            write(*,*) 'TIME   : ', TIME
            print *, 'NODE_NUMBERS=',NODES_STRANGE
            print *, 'VALUES=',STRANGERS

            print *, 'RELATED VARIABLES     :'
            write(*,*) 'K_MIN_DON_DOC_20    : ',   K_MIN_DON_DOC_20
            write(*,*) 'THETA_K_MIN_DON_DOC : ',   THETA_K_MIN_DON_DOC
            write(*,*) 'K_HS_DON_MIN_DOC    : ',   K_HS_DON_MIN_DOC

            write(*,*) 'R_ABIOTIC_DON_MIN_DOXY : ',   (R_ABIOTIC_DON_MIN_DOXY    &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'R_ABIOTIC_DON_MIN_NO3N : ',   (R_ABIOTIC_DON_MIN_NO3N    &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'R_ABIOTIC_DON_MIN_MN_IV: ',   (R_ABIOTIC_DON_MIN_MN_IV    &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'R_ABIOTIC_DON_MIN_FE_III : ', (R_ABIOTIC_DON_MIN_FE_III &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'R_ABIOTIC_DON_MIN_S_PLUS_6: ',(R_ABIOTIC_DON_MIN_S_PLUS_6 &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'R_ABIOTIC_DON_MIN_DOC     : ',(R_ABIOTIC_DON_MIN_DOC &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'R_ABIOTIC_DOC_MIN_DOC     : ',(R_ABIOTIC_DOC_MIN_DOC &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'LIM_DOC_RED               : ',(LIM_DOC_RED &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'PH_CORR_DON_MIN_DOC       : ',(PH_CORR_DON_MIN_DOC &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'DISS_ORG_N                : ',(DISS_ORG_N &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'LIM_NO3N_RED              : ',(LIM_NO3N_RED &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'LIM_MN_IV_RED             : ',(LIM_MN_IV_RED &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'LIM_FE_III_RED            : ',(LIM_FE_III_RED &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            write(*,*) 'LIM_S_PLUS_6_RED          : ',(LIM_S_PLUS_6_RED &
                                                      (NODES_STRANGE(j)),j=1,nstrange)
            deallocate(STRANGERS)
            deallocate(NODES_STRANGE)

            stop
        end if
    end if


    !DISSOLVED ORGANIC PHOSPHORUS
    PROCESS_RATES(1:nkn,DISS_ORG_P_INDEX, 1) = R_DET_PART_ORG_P_DISSOLUTION
    PROCESS_RATES(1:nkn,DISS_ORG_P_INDEX, 2) = 0. !R_ZOO_EX_DOP No excretion

    PROCESS_RATES(1:nkn,DISS_ORG_P_INDEX, 3) =  &
        R_ABIOTIC_DOP_MIN_DOXY   + R_ABIOTIC_DOP_MIN_NO3N     + R_ABIOTIC_DOP_MIN_MN_IV + &
        R_ABIOTIC_DOP_MIN_FE_III + R_ABIOTIC_DOP_MIN_S_PLUS_6 + R_ABIOTIC_DOP_MIN_DOC

    PROCESS_RATES(1:nkn,DISS_ORG_P_INDEX, 4) = &
        (R_DIA_EXCR * DIA_P_TO_C)         + (R_CYN_EXCR * CYN_P_TO_C) + &
        (R_FIX_CYN_EXCR * FIX_CYN_P_TO_C) + (R_OPA_EXCR * OPA_P_TO_C)

    DERIVATIVES(1:nkn,DISS_ORG_P_INDEX) = &
        PROCESS_RATES(1:nkn,DISS_ORG_P_INDEX, 1) + PROCESS_RATES(1:nkn,DISS_ORG_P_INDEX, 2) - &
        PROCESS_RATES(1:nkn,DISS_ORG_P_INDEX, 3) + PROCESS_RATES(1:nkn,DISS_ORG_P_INDEX, 4)


    ! Kinetic sub model for dissolved inorganic carbon

    ! Sources
    R_DIA_TOT_RESP     = PROCESS_RATES(1:nkn,DIA_C_INDEX     , 2)
    R_CYN_TOT_RESP     = PROCESS_RATES(1:nkn,CYN_C_INDEX     , 2)
    R_FIX_CYN_TOT_RESP = PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX , 2)
    R_OPA_TOT_RESP     = PROCESS_RATES(1:nkn,OPA_C_INDEX     , 2)
    R_ZOO_RESP         = PROCESS_RATES(1:nkn,ZOO_C_INDEX     , 3)
    R_ABIOTIC_DOC_MIN  = PROCESS_RATES(1:nkn,DISS_ORG_C_INDEX, 3)

    TOTAL_DIC_KINETIC_SOURCES = &
        R_DIA_TOT_RESP             + R_CYN_TOT_RESP             + &
        R_FIX_CYN_TOT_RESP         + R_OPA_TOT_RESP             + &
        R_ZOO_RESP                 + R_ABIOTIC_DOC_MIN_DOXY     + &
        R_ABIOTIC_DOC_MIN_NO3N     + R_ABIOTIC_DOC_MIN_MN_IV    + &
        R_ABIOTIC_DOC_MIN_FE_III   + R_ABIOTIC_DOC_MIN_S_PLUS_6 + &
        (0.5D0 * R_METHANOGENESIS) + R_METHANE_OXIDATION

    ! Sinks
    R_DIA_GROWTH       = PROCESS_RATES(1:nkn,DIA_C_INDEX    , 1)
    R_CYN_GROWTH       = PROCESS_RATES(1:nkn,CYN_C_INDEX    , 1)
    R_FIX_CYN_GROWTH   = PROCESS_RATES(1:nkn,FIX_CYN_C_INDEX, 1)
    R_OPA_GROWTH       = PROCESS_RATES(1:nkn,OPA_C_INDEX    , 1)

    TOTAL_DIC_KINETIC_SINKS = R_DIA_GROWTH + R_CYN_GROWTH + R_FIX_CYN_GROWTH + R_OPA_GROWTH

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
    CO2_ATM_EXHANGE = K_A_CALC_CO2 * (CO2_SAT - (H2CO3/1.0D6))

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
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 1) = TOTAL_DIC_KINETIC_SOURCES / 12000.0D0
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 2) = TOTAL_DIC_KINETIC_SINKS   / 12000.0D0

        if (CONSIDER_CO2_REARATION > 0) then
            PROCESS_RATES(1:nkn,INORG_C_INDEX, 3) = CO2_ATM_EXHANGE
        else
            PROCESS_RATES(1:nkn,INORG_C_INDEX, 3) = 0.0D0
            CO2_ATM_EXHANGE = 0.0D0
        end if

        PROCESS_RATES(1:nkn,INORG_C_INDEX, 4)  = R_DIA_TOT_RESP
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 5)  = R_CYN_TOT_RESP
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 6)  = R_FIX_CYN_TOT_RESP
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 7)  = R_OPA_TOT_RESP
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 8)  = R_ZOO_RESP
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 9)  = R_ABIOTIC_DOC_MIN
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 10) = R_DIA_GROWTH
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 11) = R_CYN_GROWTH
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 12) = R_FIX_CYN_GROWTH
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 13) = R_OPA_GROWTH

        DIC_KINETIC_DERIVATIVE = CO2_ATM_EXHANGE + &
            ((TOTAL_DIC_KINETIC_SOURCES - TOTAL_DIC_KINETIC_SINKS) / 12000.0D0)

        DERIVATIVES(1:nkn,INORG_C_INDEX) = DIC_KINETIC_DERIVATIVE
    else
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 1) = 0.0D0
        PROCESS_RATES(1:nkn,INORG_C_INDEX, 2) = 0.0D0

        if (CONSIDER_CO2_REARATION > 0) then
            PROCESS_RATES(1:nkn,INORG_C_INDEX, 3) = CO2_ATM_EXHANGE
        else
            PROCESS_RATES(1:nkn,INORG_C_INDEX, 3) = 0.0D0
            CO2_ATM_EXHANGE = 0.0D0
        end if

        PROCESS_RATES(1:nkn,INORG_C_INDEX, 4:13) = 0.0D0

        if (CONSIDER_CO2_REARATION > 0) then
            DERIVATIVES(1:nkn,INORG_C_INDEX) = CO2_ATM_EXHANGE
        else
            DERIVATIVES(1:nkn,INORG_C_INDEX) = 0.0D0
        end if
    end if

    ! -------------------------------------------------------------------------
    ! 29 JANUARY 2016, KINETIC DERIVATIVES FOR THE NEW STATE VARIABLES
    ! -------------------------------------------------------------------------
    if (DO_ADVANCED_REDOX_SIMULATION > 0) then
    
        ! Calcium
        DERIVATIVES(1:nkn,CA_INDEX) = 0.0D0

        ! Magnesium
        DERIVATIVES(1:nkn,MG_INDEX) = 0.0D0

        ! Suphate sulphur
        PROCESS_RATES(1:nkn, S_PLUS_6_INDEX, 1) = R_SULPHIDE_OXIDATION
        PROCESS_RATES(1:nkn, S_PLUS_6_INDEX, 2) = R_SULPHATE_REDUCTION
    
        DERIVATIVES(1:nkn,S_PLUS_6_INDEX) = &
            PROCESS_RATES(1:nkn, S_PLUS_6_INDEX, 1) - PROCESS_RATES(1:nkn, S_PLUS_6_INDEX, 2)

        ! Sulphide sulphur
        PROCESS_RATES(1:nkn, S_MINUS_2_INDEX, 1) = H2S_ATM_EXCHANGE
        PROCESS_RATES(1:nkn, S_MINUS_2_INDEX, 2) = R_SULPHATE_REDUCTION
        PROCESS_RATES(1:nkn, S_MINUS_2_INDEX, 3) = R_SULPHIDE_OXIDATION

        DERIVATIVES(1:nkn,S_MINUS_2_INDEX) = &
            PROCESS_RATES(1:nkn, S_MINUS_2_INDEX, 1) + PROCESS_RATES(1:nkn, S_MINUS_2_INDEX, 2) - &
            PROCESS_RATES(1:nkn, S_MINUS_2_INDEX, 3)

        ! Methane carbon
        PROCESS_RATES(1:nkn, CH4_C_INDEX, 1) = CH4_ATM_EXCHANGE
        PROCESS_RATES(1:nkn, CH4_C_INDEX, 2) = R_METHANOGENESIS
        PROCESS_RATES(1:nkn, CH4_C_INDEX, 3) = R_METHANE_OXIDATION

        DERIVATIVES(1:nkn,CH4_C_INDEX) = &
            PROCESS_RATES(1:nkn, CH4_C_INDEX, 1) + PROCESS_RATES(1:nkn, CH4_C_INDEX, 2) - &
            PROCESS_RATES(1:nkn, CH4_C_INDEX, 3)
    end if
    
    ! Debug code for  DISS_ORG_N
    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,CH4_C_INDEX),VALUE_strange,nkn).eq.1) then
            nstrange = count(VALUE_strange)
            allocate(STRANGERS    (nstrange))
            allocate(NODES_STRANGE(nstrange))
            j=1

            do k=1,nkn
                if(VALUE_strange(k)) then
                    STRANGERS    (j) = DERIVATIVES(k, CH4_C_INDEX)
                    NODES_STRANGE(j) = node_active(k)
                    j=j+1
                end if
            end do
            print *, 'PELAGIC_KINETICS:  Derivative related to CH4_C is strange'
            write(*,*) 'TIME   : ', TIME
            print *, 'NODE_NUMBERS=',NODES_STRANGE
            print *, 'VALUES=',STRANGERS

            print *, 'RELATED VARIABLES     :'
            
            write(*,*) 'CH4_ATM_EXCHANGE    : ',   & 
                (CH4_ATM_EXCHANGE(NODES_STRANGE(j)),j=1,nstrange)
            
            write(*,*) 'R_METHANOGENESIS    : ',   &
                (R_ABIOTIC_DON_MIN_NO3N(NODES_STRANGE(j)),j=1,nstrange)

            write(*,*) 'R_METHANE_OXIDATION : ',   &
                (R_ABIOTIC_DON_MIN_MN_IV(NODES_STRANGE(j)),j=1,nstrange)

            deallocate(STRANGERS)
            deallocate(NODES_STRANGE)
            stop
        end if
    end if
    
    ! -------------------------------------------------------------------------
    ! END OF 29 JANUARY 2016, LINETIC DERIVATIVES FOR THE NEW STATE VARIABLES
    ! -------------------------------------------------------------------------


    ! After the introduction of changes in November 30 th 2015, it will be necessary
    ! to reconsider alkalinity since metals will change ion balance
    ! Job to be conducted wiht Petras together, and completed during next visit of
    ! of Ali to Lithuania. Also for iron and mangenese, the redox sequences as described
    ! by Van Chappen and Wang 2015 and Katsev papers should be included.

    ! BE CAREFUL THAT THIS CHANGES WILL NEED THE METAL KITEIC PROCESS RATES TO BE
    ! CALCULATED BEFORE THE ALKALINITY SUBMODEL WHICH IS NOT THE CASE IN PRESENT
    ! SITUATION !!!!!!


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
    N_DIA_TOT_RESP            = PROCESS_RATES(1:nkn,NH4_N_INDEX, 1)  * FRAC_NH4
    N_CYN_TOT_RESP            = PROCESS_RATES(1:nkn,NH4_N_INDEX, 2)  * FRAC_NH4
    N_OPA_TOT_RESP            = PROCESS_RATES(1:nkn,NH4_N_INDEX, 3)  * FRAC_NH4
    N_FIX_CYN_TOT_RESP        = PROCESS_RATES(1:nkn,NH4_N_INDEX, 4)  * FRAC_NH4
    N_ZOO_TOT_RESP            = PROCESS_RATES(1:nkn,NH4_N_INDEX, 5)  * FRAC_NH4
    N_ABIOTIC_DON_MIN         = PROCESS_RATES(1:nkn,NH4_N_INDEX, 11) * FRAC_NH4

    ALK_GAINED_BY_AMMONIUM_GEN = &
        (N_DIA_TOT_RESP            + N_CYN_TOT_RESP         + &
         N_OPA_TOT_RESP            + N_FIX_CYN_TOT_RESP     + &
         N_ZOO_TOT_RESP            + N_ABIOTIC_DON_MIN) / 14007.0D0
    ! -------------------------------------------------------------------------
    ! End of calculate the alkalinity gain by ammonium generation
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Calculate the alkality gain by nitrate consumption
    ! 1 eq alk for each nitrate consumed since one negative ion is lost
    ! -------------------------------------------------------------------------
    N_DENITRIFICATION     = PROCESS_RATES(1:nkn,NO3_N_INDEX, 2)
    N_DIA_GROWTH          = PROCESS_RATES(1:nkn,NO3_N_INDEX, 3)
    N_CYN_GROWTH          = PROCESS_RATES(1:nkn,NO3_N_INDEX, 4)
    N_OPA_GROWTH          = PROCESS_RATES(1:nkn,NO3_N_INDEX, 5)
    N_NON_FIX_CYN_GROWTH  = PROCESS_RATES(1:nkn,NO3_N_INDEX, 6)

    ALK_GAINED_BY_NITRATE_CONS = &
        (N_DENITRIFICATION + N_DIA_GROWTH + N_CYN_GROWTH + N_OPA_GROWTH + &
         N_NON_FIX_CYN_GROWTH) / 14007.0D0
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
    N_DIA_GROWTH          = PROCESS_RATES(1:nkn,NH4_N_INDEX, 6) * FRAC_NH4
    N_CYN_GROWTH          = PROCESS_RATES(1:nkn,NH4_N_INDEX, 7) * FRAC_NH4
    N_OPA_GROWTH          = PROCESS_RATES(1:nkn,NH4_N_INDEX, 8) * FRAC_NH4
    N_NON_FIX_CYN_GROWTH  = PROCESS_RATES(1:nkn,NH4_N_INDEX, 9) * FRAC_NH4

    ALK_LOST_BY_AMMONIUM_CONS = &
        (N_DIA_GROWTH + N_CYN_GROWTH + N_OPA_GROWTH + &
		 N_NON_FIX_CYN_GROWTH) / 14007.0D0
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
    N_NITRIFICATION_NH4 = PROCESS_RATES(1:nkn,NH4_N_INDEX,10) * FRAC_NH4
    N_NITRIFICATION_NH3 = PROCESS_RATES(1:nkn,NH4_N_INDEX,10) * FRAC_NH3

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

    do i =1, nkn
        if (SALT(i) .lt. 0.d0) then
            if (SALT(i) .gt. -1.d-5) then
                SALT(i) = 0.d0
            else
                write(*,*) 'aquabc_II_pelagic_model: SALINITY IS NEGATIVE'
                print *, 'SALT=', SALT(i), 'node=', i
            end if
        end if
    end do


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
    P_DIA_GROWTH          = &
	    PROCESS_RATES(1:nkn,PO4_P_INDEX, 6) * PHOSPHATE_EQ_CONSTANT
    
	P_CYN_GROWTH          = &
	    PROCESS_RATES(1:nkn,PO4_P_INDEX, 7) * PHOSPHATE_EQ_CONSTANT
    
	P_OPA_GROWTH          = &
	    PROCESS_RATES(1:nkn,PO4_P_INDEX, 8) * PHOSPHATE_EQ_CONSTANT
    
	P_NON_FIX_CYN_GROWTH  = &
	    PROCESS_RATES(1:nkn,PO4_P_INDEX, 9) * PHOSPHATE_EQ_CONSTANT

    ALK_GAINED_BY_PHOSPHATE_CONS = &
        (P_DIA_GROWTH + P_CYN_GROWTH + P_OPA_GROWTH + &
		 P_NON_FIX_CYN_GROWTH) / 30974.0D0
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
    P_DIA_TOT_RESP     = &
	    PROCESS_RATES(1:nkn,PO4_P_INDEX, 1)  * PHOSPHATE_EQ_CONSTANT
    
	P_CYN_TOT_RESP     = &
	    PROCESS_RATES(1:nkn,PO4_P_INDEX, 2)  * PHOSPHATE_EQ_CONSTANT
    
	P_OPA_TOT_RESP     = &
	    PROCESS_RATES(1:nkn,PO4_P_INDEX, 3)  * PHOSPHATE_EQ_CONSTANT
    
	P_FIX_CYN_TOT_RESP = &
	    PROCESS_RATES(1:nkn,PO4_P_INDEX, 4)  * PHOSPHATE_EQ_CONSTANT
    
	P_ZOO_TOT_RESP     = &
	    PROCESS_RATES(1:nkn,PO4_P_INDEX, 5)  * PHOSPHATE_EQ_CONSTANT
    
	P_ABIOTIC_DOP_MIN  = &
	    PROCESS_RATES(1:nkn,PO4_P_INDEX, 10) * PHOSPHATE_EQ_CONSTANT

    ALK_LOST_BY_PHOSPHATE_GEN = &
        (P_DIA_TOT_RESP     + P_CYN_TOT_RESP + P_OPA_TOT_RESP + &
         P_FIX_CYN_TOT_RESP + P_ZOO_TOT_RESP + P_ABIOTIC_DOP_MIN) / 30974.0D0
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

    DERIVATIVES(1:nkn,TOT_ALK_INDEX) = ALK_KINETIC_DERIVATIVE

    PROCESS_RATES(1:nkn,TOT_ALK_INDEX, 1) = ALK_GAINED_BY_AMMONIUM_GEN
    PROCESS_RATES(1:nkn,TOT_ALK_INDEX, 2) = ALK_GAINED_BY_NITRATE_CONS
    PROCESS_RATES(1:nkn,TOT_ALK_INDEX, 3) = ALK_GAINED_BY_PHOSPHATE_CONS
    PROCESS_RATES(1:nkn,TOT_ALK_INDEX, 4) = ALK_LOST_BY_AMMONIUM_CONS
    PROCESS_RATES(1:nkn,TOT_ALK_INDEX, 5) = ALK_LOST_BY_NITRIFICATION
    PROCESS_RATES(1:nkn,TOT_ALK_INDEX, 6) = ALK_LOST_BY_PHOSPHATE_GEN
    PROCESS_RATES(1:nkn,TOT_ALK_INDEX, 7) = pH(1:nkn)

    if(debug_stranger) then
        if (STRANGERSD(DERIVATIVES(1:nkn,TOT_ALK_INDEX),VALUE_strange,nkn).eq.1) then
            nstrange = count(VALUE_strange)
            allocate(STRANGERS    (nstrange))
            allocate(NODES_STRANGE(nstrange))
            allocate(NODES_STRANGE_int(nstrange))
            allocate(NODES_STRANGE_ext(nstrange))
            j=1

            do k=1,nkn
                if(VALUE_strange(k)) then
                    STRANGERS    (j) = DERIVATIVES(k,TOT_ALK_INDEX)
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

    ! New kinetic derivatives calculations added 9 September 2015

    if (DO_ADVANCED_REDOX_SIMULATION > 0) then
        ! After the introduction of changes in November 30 th 2015, it will be necessary
        ! to reconsinder kinetics, especially the mineral disoolotion processes for metals
        ! Job to be conducted wiht Petras together, and completed during next visit of
        ! of Ali to Lithuania. Also for iron and mangenese, the redox sequences as described
        ! by Van Chappen and Wang 2015 and Katsev papers should be included.
        
        
        ! FE_II
        PROCESS_RATES(1:nkn,FE_II_INDEX, 1) = R_FE_III_REDUCTION
        PROCESS_RATES(1:nkn,FE_II_INDEX, 2) = R_FE_II_OXIDATION
        
        ! Diagnostics:
        PROCESS_RATES(1:nkn,FE_II_INDEX, 3) = FE_II_DISS
        PROCESS_RATES(1:nkn,FE_II_INDEX, 4) = FE_II_DISS/56000.0D0
        
        DERIVATIVES(1:nkn,FE_II_INDEX) = &
            PROCESS_RATES(1:nkn,FE_II_INDEX, 1) - PROCESS_RATES(1:nkn,FE_II_INDEX, 2)
        
        ! FE_III
        PROCESS_RATES(1:nkn,FE_III_INDEX, 1) = R_FE_II_OXIDATION
        PROCESS_RATES(1:nkn,FE_III_INDEX, 2) = R_FE_III_REDUCTION
        PROCESS_RATES(1:nkn,FE_III_INDEX, 3) = FE_III_DISS_EQ
        
        DERIVATIVES(1:nkn,FE_III_INDEX) = &
            PROCESS_RATES(1:nkn,FE_III_INDEX, 1) - PROCESS_RATES(1:nkn,FE_III_INDEX, 2)
        
        
        ! MN_II
        PROCESS_RATES(1:nkn,MN_II_INDEX, 1) = R_MN_IV_REDUCTION
        PROCESS_RATES(1:nkn,MN_II_INDEX, 2) = R_MN_II_OXIDATION
        
        DERIVATIVES(1:nkn,MN_II_INDEX) = &
            PROCESS_RATES(1:nkn,MN_II_INDEX, 1) - PROCESS_RATES(1:nkn,MN_II_INDEX, 2)
        
        
        ! MN_IV
        PROCESS_RATES(1:nkn,MN_IV_INDEX, 1) = R_MN_II_OXIDATION
        PROCESS_RATES(1:nkn,MN_IV_INDEX, 2) = R_MN_IV_REDUCTION
        
        DERIVATIVES(1:nkn,MN_IV_INDEX) = &
            PROCESS_RATES(1:nkn,MN_IV_INDEX, 1) - PROCESS_RATES(1:MN_IV_INDEX,28, 2)
        ! End of new kinetic derivatives calculations added 9 September 2015
        
        ! -------------------------------------------------------------------------
        ! END OF INSERT INTO KINETIC DERIVATIVES
        ! -------------------------------------------------------------------------
        
        
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
    end if
    !*********************************************'
    !*                                           *'
    !*          END OF ECOLOGY KINETICS          *'
    !*                                           *'
    !*      DO NOT CHANGE THE FOLLOWING CODE     *'
    !*                                           *'
    !*********************************************'

end subroutine PELAGIC_KINETICS



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

    where (NH3 .lt. 1.0D-6)
        PN = 0.0D0
    elsewhere
        PN = (NH3 * NOx) / ((kn + NH3) * (kn + NOx)) + (kn * NH3) / ((NH3 + NOx) * (kn + NOx))
    end where

    AMMONIA_PREF = PN

end subroutine AMMONIA_PREFS


!Subroutine to calculate dissolved nitroge preference (WASP)
subroutine AMMONIA_DON_PREFS(AMMONIA_DON_PREF,NH3,  DON, frac_avail_DON,NOx, kn,nkn)

    integer nkn

    double precision :: AMMONIA_DON_PREF(nkn)
    double precision, intent(in) :: NH3(nkn)
    double precision, intent(in) :: DON(nkn)
    double precision, intent(in) :: frac_avail_DON
    double precision, intent(in) :: NOx(nkn)
    double precision, intent(in) :: kn
    double precision             :: PN (nkn)

    double precision :: NH3_AND_AVAIL_DON(nkn)


    NH3_AND_AVAIL_DON = NH3 + (frac_avail_DON * DON)

    where (NH3_AND_AVAIL_DON .lt. 1.0D-6)
     PN = 0.0D0
    elsewhere
     PN = (NH3_AND_AVAIL_DON * NOx) / ((kn + NH3_AND_AVAIL_DON) * (kn + NOx)) + &
        (kn * NH3_AND_AVAIL_DON) / ((NH3_AND_AVAIL_DON + NOx) * (kn + NOx))
    end where

    AMMONIA_DON_PREF = PN

end subroutine AMMONIA_DON_PREFS


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
             (66423080.0 / (T_KELVIN ** 2.0D0)) + (12438000000.0 / (T_KELVIN ** 3.0D0)) - &
             (862194900000.0 / (T_KELVIN ** 4.0D0))

    !Calculate the effect of salinity on dissolved oxygen saturation
    LN_CSS = LN_CSF - S * &
             (0.017674 - (10.754 / T_KELVIN) + (2140.7 / (T_KELVIN ** 2.0D0)))

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

    use PELAGIC_MODEL_CONSTANTS
    !use para_aqua

    implicit none
    integer nkn  ! number of nodes in grid

    double precision Ia    (nkn)
    double precision TCHLA (nkn)
    double precision GITMAX(nkn)
    double precision H     (nkn)
    double precision ke    (nkn)
    double precision CCHL_RATIO
    double precision KC
    !double precision PHIMX
    double precision LLIGHT   (nkn)
    double precision LIGHT_SAT(nkn)
    double precision KESHD    (nkn)
    double precision SKE      (nkn)
    double precision TEMP1    (nkn)
    double precision TEMP2    (nkn)
    double precision TEMP3    (nkn)
    double precision PI

    logical VALUE_strange(nkn) ! array containing 'true' on strange values
    integer STRANGERSD         ! function searching strange values

    PI = 3.14159D0
    !call para_get_value('XKC', KC) !KC SHOULD BE XKC
    !call para_get_value('PHIMX', PHIMX)

    KESHD     = XKC * TCHLA
    SKE       = ke
    SKE       = SKE + KESHD
    TEMP1     = SKE * H
    
    !1/TEMP2 - the saturating light intensity
    !TEMP2     = (0.083D0 * PHIMX * KC) / (GITMAX * CCHL_RATIO * 2.718D0)
    TEMP2     = (0.083D0 * PHIMX * XKC) / (GITMAX * CCHL_RATIO * 2.718D0)
    LIGHT_SAT = 1.0D0/TEMP2

    if (STRANGERSD(LIGHT_SAT,VALUE_strange,nkn).eq.1) then
        write(6,*) 'LIM_LIGHT: TEMP2 is NaN '
        write(6,*)  'LIGHT_SAT=',LIGHT_SAT
        write(6,*) 'TEMP2=',TEMP2
        write(6,*) 'GITMAX=', GITMAX
        write(6,*) 'CCHL_RATIO=', CCHL_RATIO
        stop
    end if

    TEMP3 = EXP( - TEMP1)  ! fraction of the light at bottom
    LLIGHT = (2.7183D0 / TEMP1) * (EXP( -TEMP2 * Ia * TEMP3) - EXP( -TEMP2 * Ia))
end subroutine LIM_LIGHT

!********************************************************************
!********************************************************************

!********************************************************************
subroutine CUR_SMITH(Ia,TCHLA,CCHLXI,GITMAX,H,ke, LLIGHT,CCHLX)

    !    Can not be used for instanteneous light
    !    while total dayllight is unknown to adjust C to Chla ratio!

    !   Dick-Smith light limitation formulation(adapted from EUTRO).
    !
    !   Version to use in ALUKAS with instanteneous light intensity and variable 
    !   (calculated inside) C to Chla ratio
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
    integer :: I


    PI = 3.14159D0

    ! Chloroph. extinction, ( mcg Chla/l/m)
    ! 0.04 is approximatelly the value that corresponds to the angle of curve given in Chapra
    XKC   =  0.016                           
    PHIMX = 720.0   !PHY   Quantum yield const. mg C/mole photon
    FDAY  = 1.0D0
    ITOT  = Ia
    CCHL1 = CCHLXI

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%?IAVBOT=IAVBOTX(ISEG)
    !%PHYT=TPHY
    !%GITMAX=k1c*rtmult(TEMP, t_1, t_2, t_3, t_4, k_1, 0.98, 0.98, k_4)
    !         %TCHLA = PHYT/CCHL1
    !         %RLIGHT = RLGHTS (ISEG, 1)
    !         %
    !         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    KESHD = XKC * TCHLA   !Chla commes in micrograms
    SKE = ke
    SKE = SKE + KESHD
    TEMP1 = SKE * H

    if (GITMAX .lt. 1D-20 .or. CCHL1 .lt. 1D-20) then
        write(6,*) 'SMITH: TEMP2 is NaN ', 'GITMAX=', GITMAX,'CCHL1=', CCHL1
        stop
    end if

    !1/TEMP2 - the saturating light intensity
    TEMP2 = 0.083D0 * PHIMX * XKC / (GITMAX * CCHL1 * 2.718D0) 
    TEMP3 = EXP( - TEMP1) !fraction of the light at bottom

    !Light limitation varies during the day
    RLIGHT = 2.7183D0 / TEMP1 * (EXP( -TEMP2 * ITOT * TEMP3) - EXP( -TEMP2 * ITOT))
    LLIGHT = RLIGHT

    !Adapt carbon to chlorophyll ratio:
    
    ! It can not be used for instantenous light because total light is unknown!
    IAV=0.9D0 * ITOT/FDAY              
    IAVSG=IAV*(1.0D0-TEMP3)/TEMP1
    CCHLX=0.3D0 * 0.083D0 * PHIMX * XKC * IAVSG / (GITMAX * 2.718D0)

    if (CCHLX.LT.15.0D0) then
        CCHLX=15.0D0
    end if
end subroutine CUR_SMITH

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
    !AMMONIA_VOLATIL_RATE = KA * ((NH3S * (1.4D1/1.70D1)) - NH3N) * ((3.2D1/1.7D1)**2.5D-1)
    AMMONIA_VOLATIL_RATE = KA * (NH3N - (NH3S * (1.4D1/1.70D1))) * ((3.2D1/1.7D1)**2.5D-1)
end subroutine AMMONIA_VOLATILIZATION

!************************************************************************
!************************************************************************

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

!************************************************************************
!************************************************************************

! ---------------------------------------------------------------------------- !
! Subroutine that calculates the PH_CORRECTION factor according to             !
! the euqations (Eq. 162 and Eq. 163) in AQUATOX Relase 3.1 Plus documentation !
! (EPA-829-R-14-007)                                                           !
! Specialised for WC by Petras                                                 !
! ---------------------------------------------------------------------------- !
!                                DOCUMENTATION                                 !
!                                                                              !
! INPUTS :                                                                     !
! --------                                                                     !
! PH_MIN  : Minimum value of pH for the optimum range (scalar)                 !
! PH_MAX  : Maximum value of pH for the optimum range (scalar)                 !
! PH      : PH values for which correction factor will be calculated (matrix)  !
! nkn     : Number of rows in matrices PH and PH_CORR                          !
! nlay    : Number of columns in in matrices PH and PH_CORR                    !
!                                                                              !
! OUTPUTS                                                                      !
! PH_CORR : pH correction factors calculated for PH (matrix)                   !
! ---------------------------------------------------------------------------- !
!                    Development date : 31st August 2015                       !
!                                                                              !
!                            Developer : Ali Ertürk                            !
! ---------------------------------------------------------------------------- !
subroutine CALCULATE_PH_CORR(PH_CORR, PH, PH_MIN, PH_MAX, nkn)

    use AQUABC_II_GLOBAL
    implicit none

    ! Ingoing variables
    real(kind = DBL_PREC), dimension(nkn), intent(in) :: PH
    real(kind = DBL_PREC), intent(in) :: PH_MIN
    real(kind = DBL_PREC), intent(in) :: PH_MAX
    integer, intent(in) :: nkn
    ! Outgoing variables
    real(kind = DBL_PREC), dimension(nkn), intent(inout) :: PH_CORR

    integer i,j,error

    error   = 0
    PH_CORR = 1

    !return

    where(PH(:) < PH_MIN)
        PH_CORR(:) = exp(PH(:) - PH_MIN)
    end where

    where(PH(:) > PH_MAX)
        PH_CORR(:) = exp(PH_MAX - PH(:))
    end where

    do i=1,nkn

        if (PH_CORR(i) .gt. 1.D0 .or. PH_CORR(i) .lt. 0.D0) then
            print *, 'CALCULATE_PH_CORR: Incorrect PH_CORR'
            print *, 'PH_CORR',PH_CORR(i)
            print *, 'PH_MIN:', PH_MIN
            print *, 'PH_MAX:', PH_MAX
            print *, 'PH:', PH(i)
            print *, 'Internal node number:',i
            error=1
        end if

        if (error .eq. 1) stop
    end do


end subroutine CALCULATE_PH_CORR

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
    !  SETTLING_RATES       - potential flux by settling, g/m2/day
    !  FLUXES               - after substraction of not deposited fluxes for each BS state variable in g/m2/day
    !  NOT_DEPOSITED_FLUXES - for each WC state variable, g/m2/day

    use AQUABC_II_GLOBAL
    use PELAGIC_MODEL_CONSTANTS
    use AQUABC_PEL_STATE_VAR_INDEXES

    implicit none

    integer NUM_VARS     !Number of WC state variables
    integer NUM_CONSTS
    integer NUM_FLUXES   !Number of fluxes to BS == number of BS state variables
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
    DIMENSION FLUXES(NUM_FLUXES) !Fluxes to BS for BS variables

    double precision FLUXES_FROM_WC
    DIMENSION FLUXES_FROM_WC(NUM_VARS)   !Fluxes from WC to BS for WC variables

    double precision DRIVING_FUNCTIONS
    DIMENSION DRIVING_FUNCTIONS(NUM_DRIV)

    integer CELLNO
    integer LAYER

    double precision PSTIME
    double precision BOTTOM_FACTOR ! never used other value than 1

    integer SEDIMENT_TYPE
    integer NUM_NOT_DEPOSITED_FLUXES
    double precision FRACTION_OF_DEPOSITION(NUM_NOT_DEPOSITED_FLUXES)
    double precision NOT_DEPOSITED_FLUXES(NUM_NOT_DEPOSITED_FLUXES)

    double precision SETTLING_FACTORS
    DIMENSION SETTLING_FACTORS(NUM_VARS)

    integer I

    
    if(num_vars .ne. 30) then
        print *, 'FLX_ALUKAS_II_TO_SED_MOD_1:'
        print *, 'To get values correctly by fluxes from WC and not deposited fluxes'
        print *, 'number of state variables should be equal to 30 but is ',NUM_VARS
        stop
    end if

    do i = 1, NUM_VARS
        SETTLING_FACTORS(i) = &
            BOTTOM_FACTOR * &                    !bottom factor =1
                (1.0D+0 - DISSOLVED_FRACTIONS(i)) * SETTLING_VELOCITIES(i)

        if (SETTLING_FACTORS(i).lt.0.0D0) then
            SETTLING_FACTORS(i) = 0.0D0
        end if

        SETTLING_RATES(i)   = &
            STATE_VARIABLES(i) *  (1.0D+0 - DISSOLVED_FRACTIONS(i))* &
            SETTLING_VELOCITIES(i)
    end do

    FLUXES_FROM_WC(1:NUM_VARS) = STATE_VARIABLES(1:NUM_VARS) * SETTLING_FACTORS(1:NUM_VARS)

    !NOT DEPOSITED FLUXES
    NOT_DEPOSITED_FLUXES(1:NUM_VARS) = FLUXES_FROM_WC(1:NUM_VARS) * (1.0D0 - FRACTION_OF_DEPOSITION(1:NUM_VARS))


    ! FLUXES FROM WC TO BS FOR BS VARIABLES
    if(NUM_FLUXES .ne. 24) then
        print *, 'FLX_ALUKAS_II_TO_SED_MOD_1:'
        print *, 'To get values correctly by fluxes to sediments and not deposited fluxes'
        print *, 'number of BS state vars(fluxes) should be equal to 24 but is ',NUM_FLUXES
        stop
    end if

    !AMMONIA NITROGEN FLUX
    FLUXES(1) = FLUXES_FROM_WC(NH4_N_INDEX) * FRACTION_OF_DEPOSITION(NH4_N_INDEX)

    !NITRATE NITROGEN FLUX
    FLUXES(2) = FLUXES_FROM_WC(NO3_N_INDEX) * FRACTION_OF_DEPOSITION(NO3_N_INDEX)

    !DISSOLVED ORGANIC NITROGEN FLUX
    FLUXES(3) = FLUXES_FROM_WC(DISS_ORG_N_INDEX) * FRACTION_OF_DEPOSITION(DISS_ORG_N_INDEX)

    !PARTICULATE ORGANIC NITROGEN FLUX
    FLUXES(4) = &
       (FLUXES_FROM_WC(DIA_C_INDEX)          * FRACTION_OF_DEPOSITION(DIA_C_INDEX)     * DIA_N_TO_C    ) + &
       (FLUXES_FROM_WC(CYN_C_INDEX)          * FRACTION_OF_DEPOSITION(CYN_C_INDEX)     * CYN_N_TO_C    ) + &
       (FLUXES_FROM_WC(OPA_C_INDEX)          * FRACTION_OF_DEPOSITION(OPA_C_INDEX)     * OPA_N_TO_C    ) + &
       (FLUXES_FROM_WC(FIX_CYN_C_INDEX)      * FRACTION_OF_DEPOSITION(FIX_CYN_C_INDEX) * FIX_CYN_N_TO_C) + &
       (FLUXES_FROM_WC(ZOO_N_INDEX)          * FRACTION_OF_DEPOSITION(ZOO_N_INDEX)                     ) + &
       (FLUXES_FROM_WC(DET_PART_ORG_N_INDEX) * FRACTION_OF_DEPOSITION(DET_PART_ORG_N_INDEX)            )

    !PHOSPHATE FLUX
    FLUXES(5) = FLUXES_FROM_WC(PO4_P_INDEX) * FRACTION_OF_DEPOSITION(PO4_P_INDEX)

    !DISSOLVED ORGANIC PHOSPHORUS FLUX
    FLUXES(6) = FLUXES_FROM_WC(DISS_ORG_P_INDEX) * FRACTION_OF_DEPOSITION(DISS_ORG_P_INDEX)

    !PARTICULATE ORGANIC PHOSPHORUS FLUX
    FLUXES(7) = &
       (FLUXES_FROM_WC(DIA_C_INDEX)          * FRACTION_OF_DEPOSITION(DIA_C_INDEX)     * DIA_P_TO_C    )  + &
       (FLUXES_FROM_WC(CYN_C_INDEX)          * FRACTION_OF_DEPOSITION(CYN_C_INDEX)     * CYN_P_TO_C    )  + &
       (FLUXES_FROM_WC(OPA_C_INDEX)          * FRACTION_OF_DEPOSITION(OPA_C_INDEX)     * OPA_P_TO_C    )  + &
       (FLUXES_FROM_WC(FIX_CYN_C_INDEX)      * FRACTION_OF_DEPOSITION(FIX_CYN_C_INDEX) * FIX_CYN_P_TO_C)  + &
       (FLUXES_FROM_WC(ZOO_P_INDEX)          * FRACTION_OF_DEPOSITION(ZOO_P_INDEX)                     )  + &
       (FLUXES_FROM_WC(DET_PART_ORG_P_INDEX) * FRACTION_OF_DEPOSITION(DET_PART_ORG_P_INDEX)            )

    !DISSOLVED OXYGEN
    FLUXES(8) = 0.0D0

    !DISSOLVED ORGANIC CARBON
    FLUXES(9) = FLUXES_FROM_WC(DISS_ORG_C_INDEX) * FRACTION_OF_DEPOSITION(DISS_ORG_C_INDEX)

    !PARTICULATE ORGANIC CARBON FLUX
    FLUXES(10) = &
       (FLUXES_FROM_WC(DIA_C_INDEX)          * FRACTION_OF_DEPOSITION(DIA_C_INDEX)         ) + &
       (FLUXES_FROM_WC(CYN_C_INDEX)          * FRACTION_OF_DEPOSITION(CYN_C_INDEX)         ) + &
       (FLUXES_FROM_WC(OPA_C_INDEX)          * FRACTION_OF_DEPOSITION(OPA_C_INDEX)         ) + &
       (FLUXES_FROM_WC(FIX_CYN_C_INDEX)      * FRACTION_OF_DEPOSITION(FIX_CYN_C_INDEX)     ) + &
       (FLUXES_FROM_WC(ZOO_P_INDEX)          * FRACTION_OF_DEPOSITION(ZOO_P_INDEX)         ) + &
       (FLUXES_FROM_WC(DET_PART_ORG_C_INDEX) * FRACTION_OF_DEPOSITION(DET_PART_ORG_C_INDEX))

    !DISSOLVED SILICON
    FLUXES(11) = FLUXES_FROM_WC(DISS_Si_INDEX)   * FRACTION_OF_DEPOSITION(DISS_Si_INDEX)

    !PARTICULATE  SILICON FLUX
    FLUXES(12) = &
       (FLUXES_FROM_WC(DIA_C_INDEX)   * FRACTION_OF_DEPOSITION(DIA_C_INDEX) * DIA_Si_TO_C) + &
       (FLUXES_FROM_WC(PART_Si_INDEX) * FRACTION_OF_DEPOSITION(PART_Si_INDEX))

    FLUXES(13) = FLUXES_FROM_WC(INORG_C_INDEX) * FRACTION_OF_DEPOSITION(INORG_C_INDEX) !INORGANIC CARBON FLUX
    FLUXES(14) = FLUXES_FROM_WC(TOT_ALK_INDEX) * FRACTION_OF_DEPOSITION(TOT_ALK_INDEX) !ALKALINITY FLUX
    FLUXES(15) = 0.D0 !Salinity flux

    FLUXES(16) =  FLUXES_FROM_WC(FE_II_INDEX    ) *  FRACTION_OF_DEPOSITION(FE_II_INDEX    ) !FEII  FLUX
    FLUXES(17) =  FLUXES_FROM_WC(FE_III_INDEX   ) *  FRACTION_OF_DEPOSITION(FE_III_INDEX   ) !FEIII FLUX
    FLUXES(18) =  FLUXES_FROM_WC(MN_II_INDEX    ) *  FRACTION_OF_DEPOSITION(MN_II_INDEX    ) !MNII  FLUX
    FLUXES(19) =  FLUXES_FROM_WC(MN_IV_INDEX    ) *  FRACTION_OF_DEPOSITION(MN_IV_INDEX    ) !MNIV  FLUX
    FLUXES(20) =  FLUXES_FROM_WC(CA_INDEX       ) *  FRACTION_OF_DEPOSITION(CA_INDEX       ) !CA    FLUX
    FLUXES(21) =  FLUXES_FROM_WC(MG_INDEX       ) *  FRACTION_OF_DEPOSITION(MG_INDEX       ) !MG    FLUX
    FLUXES(22) =  FLUXES_FROM_WC(S_PLUS_6_INDEX ) *  FRACTION_OF_DEPOSITION(S_PLUS_6_INDEX ) !S_PLUS_6  FLUX
    FLUXES(23) =  FLUXES_FROM_WC(S_MINUS_2_INDEX) *  FRACTION_OF_DEPOSITION(S_MINUS_2_INDEX) !S_MINUS_2 FLUX
    FLUXES(24) =  FLUXES_FROM_WC(CH4_C_INDEX    ) *  FRACTION_OF_DEPOSITION(CH4_C_INDEX    ) !CH4_C     FLUX
end subroutine FLX_ALUKAS_II_TO_SED_MOD_1

!********************************************************************
!********************************************************************

subroutine FLX_ALUKAS_II_TO_SED_MOD_1_VEC &
           (STATE_VARIABLES    , nkn        , NUM_VARS, &
            MODEL_CONSTANTS    , NUM_CONSTS ,         &
            DRIVING_FUNCTIONS  , NUM_DRIV   ,         &
            SETTLING_VELOCITIES, DISSOLVED_FRACTIONS, &
            BOTTOM_FACTOR , CELLNO, PSTIME,           &
            SETTLING_RATES, FLUXES, NUM_FLUXES,       &
            SEDIMENT_TYPE , FRACTION_OF_DEPOSITION,   &
            NOT_DEPOSITED_FLUXES, NUM_NOT_DEPOSITED_FLUXES)

    ! Outputs:
    !  SETTLING_RATES       - potential flux by settling, g/m2/day
    !  FLUXES               - after substraction of not deposited fluxes for each BS state variable in g/m2/day
    !  NOT_DEPOSITED_FLUXES - for each WC state variable, g/m2/day

    use AQUABC_II_GLOBAL
    use PELAGIC_MODEL_CONSTANTS
    use AQUABC_PEL_STATE_VAR_INDEXES

    implicit none

    integer nkn
    integer NUM_VARS     !Number of WC state variables
    integer NUM_CONSTS
    integer NUM_FLUXES   !Number of fluxes to BS == number of BS state variables
    integer NUM_DRIV

    double precision STATE_VARIABLES
    DIMENSION STATE_VARIABLES(nkn, NUM_VARS)

    double precision MODEL_CONSTANTS
    DIMENSION MODEL_CONSTANTS(NUM_CONSTS)

    double precision SETTLING_VELOCITIES
    DIMENSION SETTLING_VELOCITIES(nkn, NUM_VARS)

    double precision DISSOLVED_FRACTIONS
    DIMENSION DISSOLVED_FRACTIONS(nkn, NUM_VARS)

    double precision SETTLING_RATES
    DIMENSION SETTLING_RATES(nkn, NUM_VARS)

    double precision FLUXES
    DIMENSION FLUXES(nkn, NUM_FLUXES) !Fluxes to BS for BS variables

    double precision FLUXES_FROM_WC
    DIMENSION FLUXES_FROM_WC(nkn, NUM_VARS)   !Fluxes from WC to BS for WC variables

    double precision DRIVING_FUNCTIONS
    DIMENSION DRIVING_FUNCTIONS(nkn, NUM_DRIV)

    integer CELLNO
    integer LAYER

    double precision PSTIME
    double precision BOTTOM_FACTOR ! never used other value than 1

    integer SEDIMENT_TYPE
    integer NUM_NOT_DEPOSITED_FLUXES
    double precision FRACTION_OF_DEPOSITION(nkn, NUM_NOT_DEPOSITED_FLUXES)
    double precision NOT_DEPOSITED_FLUXES  (nkn, NUM_NOT_DEPOSITED_FLUXES)

    double precision SETTLING_FACTORS
    DIMENSION SETTLING_FACTORS(nkn, NUM_VARS)

    integer I
    
    if(num_vars .ne. 30) then
        print *, 'FLX_ALUKAS_II_TO_SED_MOD_1:'
        print *, 'To get values correctly by fluxes from WC and not deposited fluxes'
        print *, 'number of state variables should be equal to 30 but is ',NUM_VARS
        stop
    end if

    ! FLUXES FROM WC TO BS FOR BS VARIABLES
    if(NUM_FLUXES .ne. 24) then
        print *, 'FLX_ALUKAS_II_TO_SED_MOD_1:'
        print *, 'To get values correctly by fluxes to sediments and not deposited fluxes'
        print *, 'number of BS state vars(fluxes) should be equal to 24 but is ',NUM_FLUXES
        stop
    end if

    SETTLING_FACTORS = BOTTOM_FACTOR * (1.0D+0 - DISSOLVED_FRACTIONS) * SETTLING_VELOCITIES
    
    where (SETTLING_FACTORS.lt.0.0D0)
        SETTLING_FACTORS = 0.0D0
    end where

    SETTLING_RATES       = STATE_VARIABLES *  (1.0D+0 - DISSOLVED_FRACTIONS) * SETTLING_VELOCITIES
    FLUXES_FROM_WC       = STATE_VARIABLES * SETTLING_FACTORS
    NOT_DEPOSITED_FLUXES = FLUXES_FROM_WC  * (1.0D0 - FRACTION_OF_DEPOSITION)

    !AMMONIA NITROGEN FLUX
    FLUXES(:, 1) = FLUXES_FROM_WC(:, NH4_N_INDEX) * FRACTION_OF_DEPOSITION(:, NH4_N_INDEX)

    !NITRATE NITROGEN FLUX
    FLUXES(:, 2) = FLUXES_FROM_WC(:, NO3_N_INDEX) * FRACTION_OF_DEPOSITION(:, NO3_N_INDEX)

    !DISSOLVED ORGANIC NITROGEN FLUX
    FLUXES(:, 3) = FLUXES_FROM_WC(:, DISS_ORG_N_INDEX) * FRACTION_OF_DEPOSITION(:, DISS_ORG_N_INDEX)

    !PARTICULATE ORGANIC NITROGEN FLUX
    FLUXES(:, 4) = &
       (FLUXES_FROM_WC(:, DIA_C_INDEX)          * FRACTION_OF_DEPOSITION(:, DIA_C_INDEX)     * DIA_N_TO_C    ) + &
       (FLUXES_FROM_WC(:, CYN_C_INDEX)          * FRACTION_OF_DEPOSITION(:, CYN_C_INDEX)     * CYN_N_TO_C    ) + &
       (FLUXES_FROM_WC(:, OPA_C_INDEX)          * FRACTION_OF_DEPOSITION(:, OPA_C_INDEX)     * OPA_N_TO_C    ) + &
       (FLUXES_FROM_WC(:, FIX_CYN_C_INDEX)      * FRACTION_OF_DEPOSITION(:, FIX_CYN_C_INDEX) * FIX_CYN_N_TO_C) + &
       (FLUXES_FROM_WC(:, ZOO_N_INDEX)          * FRACTION_OF_DEPOSITION(:, ZOO_N_INDEX)                     ) + &
       (FLUXES_FROM_WC(:, DET_PART_ORG_N_INDEX) * FRACTION_OF_DEPOSITION(:, DET_PART_ORG_N_INDEX)            )

    !PHOSPHATE FLUX
    FLUXES(:, 5) = FLUXES_FROM_WC(:, PO4_P_INDEX) * FRACTION_OF_DEPOSITION(:, PO4_P_INDEX)

    !DISSOLVED ORGANIC PHOSPHORUS FLUX
    FLUXES(:, 6) = FLUXES_FROM_WC(:, DISS_ORG_P_INDEX) * FRACTION_OF_DEPOSITION(:, DISS_ORG_P_INDEX)

    !PARTICULATE ORGANIC PHOSPHORUS FLUX
    FLUXES(:, 7) = &
       (FLUXES_FROM_WC(:, DIA_C_INDEX)          * FRACTION_OF_DEPOSITION(:, DIA_C_INDEX)     * DIA_P_TO_C    )  + &
       (FLUXES_FROM_WC(:, CYN_C_INDEX)          * FRACTION_OF_DEPOSITION(:, CYN_C_INDEX)     * CYN_P_TO_C    )  + &
       (FLUXES_FROM_WC(:, OPA_C_INDEX)          * FRACTION_OF_DEPOSITION(:, OPA_C_INDEX)     * OPA_P_TO_C    )  + &
       (FLUXES_FROM_WC(:, FIX_CYN_C_INDEX)      * FRACTION_OF_DEPOSITION(:, FIX_CYN_C_INDEX) * FIX_CYN_P_TO_C)  + &
       (FLUXES_FROM_WC(:, ZOO_P_INDEX)          * FRACTION_OF_DEPOSITION(:, ZOO_P_INDEX)                     )  + &
       (FLUXES_FROM_WC(:, DET_PART_ORG_P_INDEX) * FRACTION_OF_DEPOSITION(:, DET_PART_ORG_P_INDEX)            )

    !DISSOLVED OXYGEN
    FLUXES(:, 8) = 0.0D0

    !DISSOLVED ORGANIC CARBON
    FLUXES(:, 9) = FLUXES_FROM_WC(:, DISS_ORG_C_INDEX) * FRACTION_OF_DEPOSITION(:, DISS_ORG_C_INDEX)

    !PARTICULATE ORGANIC CARBON FLUX
    FLUXES(:, 10) = &
       (FLUXES_FROM_WC(:, DIA_C_INDEX)          * FRACTION_OF_DEPOSITION(:, DIA_C_INDEX)         ) + &
       (FLUXES_FROM_WC(:, CYN_C_INDEX)          * FRACTION_OF_DEPOSITION(:, CYN_C_INDEX)         ) + &
       (FLUXES_FROM_WC(:, OPA_C_INDEX)          * FRACTION_OF_DEPOSITION(:, OPA_C_INDEX)         ) + &
       (FLUXES_FROM_WC(:, FIX_CYN_C_INDEX)      * FRACTION_OF_DEPOSITION(:, FIX_CYN_C_INDEX)     ) + &
       (FLUXES_FROM_WC(:, ZOO_P_INDEX)          * FRACTION_OF_DEPOSITION(:, ZOO_P_INDEX)         ) + &
       (FLUXES_FROM_WC(:, DET_PART_ORG_C_INDEX) * FRACTION_OF_DEPOSITION(:, DET_PART_ORG_C_INDEX))

    !DISSOLVED SILICON
    FLUXES(:, 11) = FLUXES_FROM_WC(:, DISS_Si_INDEX)   * FRACTION_OF_DEPOSITION(:, DISS_Si_INDEX)

    !PARTICULATE  SILICON FLUX
    FLUXES(:, 12) = &
       (FLUXES_FROM_WC(:, DIA_C_INDEX)   * FRACTION_OF_DEPOSITION(:, DIA_C_INDEX) * DIA_Si_TO_C) + &
       (FLUXES_FROM_WC(:, PART_Si_INDEX) * FRACTION_OF_DEPOSITION(:, PART_Si_INDEX))

    FLUXES(:, 13) = FLUXES_FROM_WC(:, INORG_C_INDEX) * FRACTION_OF_DEPOSITION(:, INORG_C_INDEX) !INORGANIC CARBON FLUX
    FLUXES(:, 14) = FLUXES_FROM_WC(:, TOT_ALK_INDEX) * FRACTION_OF_DEPOSITION(:, TOT_ALK_INDEX) !ALKALINITY FLUX
    FLUXES(:, 15) = 0.D0 !Salinity flux

    FLUXES(:, 16) =  FLUXES_FROM_WC(:, FE_II_INDEX    ) *  FRACTION_OF_DEPOSITION(:, FE_II_INDEX    ) !FEII  FLUX
    FLUXES(:, 17) =  FLUXES_FROM_WC(:, FE_III_INDEX   ) *  FRACTION_OF_DEPOSITION(:, FE_III_INDEX   ) !FEIII FLUX
    FLUXES(:, 18) =  FLUXES_FROM_WC(:, MN_II_INDEX    ) *  FRACTION_OF_DEPOSITION(:, MN_II_INDEX    ) !MNII  FLUX
    FLUXES(:, 19) =  FLUXES_FROM_WC(:, MN_IV_INDEX    ) *  FRACTION_OF_DEPOSITION(:, MN_IV_INDEX    ) !MNIV  FLUX
    FLUXES(:, 20) =  FLUXES_FROM_WC(:, CA_INDEX       ) *  FRACTION_OF_DEPOSITION(:, CA_INDEX       ) !CA    FLUX
    FLUXES(:, 21) =  FLUXES_FROM_WC(:, MG_INDEX       ) *  FRACTION_OF_DEPOSITION(:, MG_INDEX       ) !MG    FLUX
    FLUXES(:, 22) =  FLUXES_FROM_WC(:, S_PLUS_6_INDEX ) *  FRACTION_OF_DEPOSITION(:, S_PLUS_6_INDEX ) !S_PLUS_6  FLUX
    FLUXES(:, 23) =  FLUXES_FROM_WC(:, S_MINUS_2_INDEX) *  FRACTION_OF_DEPOSITION(:, S_MINUS_2_INDEX) !S_MINUS_2 FLUX
    FLUXES(:, 24) =  FLUXES_FROM_WC(:, CH4_C_INDEX    ) *  FRACTION_OF_DEPOSITION(:, CH4_C_INDEX    ) !CH4_C     FLUX
end subroutine FLX_ALUKAS_II_TO_SED_MOD_1_VEC

!********************************************************************
!********************************************************************

integer function STRANGERSD(VALUE,VALUE_strange,nkn)

    use mod_debug

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
      VALUE_NAN     = is_nan(VALUE)
      VALUE_Inf     = abs(VALUE) > BIGNUMBER
      VALUE_strange = VALUE_NAN .or. VALUE_Inf .or. VALUE_out

      if(any(VALUE_strange)) STRANGERSD = 1
      return

end function STRANGERSD
!************************************************************************
!************************************************************************
