
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Main routines for bottom sediment model No1
! CONTENT:
!  module SED_MODEL_CONSTANTS
!  subroutine INIT_SED_MODEL_CONSTANTS
!  subroutine SEDIMENT_MODEL_1   -  main routine
!  FUNCTION STRANGER(VALUE)      -  checks for NaNs (input type is double precision), not vectorised
!  subroutine FLX_ALUKAS_II_TO_SED_MOD_1
!  FUNCTION SED_MOD_1_ALUKAS_MOLDI_C - molecular diffusion coefficients
!  subroutine SED_MOD_1_CVISC
!  subroutine FLX_SED_MOD_1_TO_ALUKAS_II
!  subroutine CALCULATE_PH_CORR_SED

module BSED_MODEL_CONSTANTS
    use AQUABC_II_GLOBAL

    !SEDIMENT MODEL COEFFICIENTS
    !Dissolution of particulate organic carbon
    real(kind = DBL_PREC) K_OXIC_DISS_POC
    real(kind = DBL_PREC) K_ANOXIC_DISS_POC
    real(kind = DBL_PREC) THETA_DISS_POC
    real(kind = DBL_PREC) KHS_DISS_POC

    !Dissolution of particulate organic nitrogen
    real(kind = DBL_PREC) K_OXIC_DISS_PON
    real(kind = DBL_PREC) K_ANOXIC_DISS_PON
    real(kind = DBL_PREC) THETA_DISS_PON
    real(kind = DBL_PREC) KHS_DISS_PON

    !Dissolution of particulate organic phosphorus
    real(kind = DBL_PREC) K_OXIC_DISS_POP
    real(kind = DBL_PREC) K_ANOXIC_DISS_POP
    real(kind = DBL_PREC) THETA_DISS_POP
    real(kind = DBL_PREC) KHS_DISS_POP

    !Dissolution of particulate silicon
    real(kind = DBL_PREC) K_OXIC_DISS_PSi
    real(kind = DBL_PREC) K_ANOXIC_DISS_PSi
    real(kind = DBL_PREC) THETA_DISS_PSi
    real(kind = DBL_PREC) KHS_DISS_PSi

    !Mineralization of dissolved organic carbon
    real(kind = DBL_PREC) K_OXIC_MINER_DOC
    real(kind = DBL_PREC) K_ANOXIC_MINER_DOC
    real(kind = DBL_PREC) THETA_MINER_DOC
    real(kind = DBL_PREC) KHS_MINER_DOC
    real(kind = DBL_PREC) O_TO_C

    !Mineralization of dissolved organic nitrogen
    real(kind = DBL_PREC) K_OXIC_MINER_DON
    real(kind = DBL_PREC) K_ANOXIC_MINER_DON
    real(kind = DBL_PREC) THETA_MINER_DON
    real(kind = DBL_PREC) KHS_MINER_DON

    !Mineralization of dissolved organic phosphorus
    real(kind = DBL_PREC) K_OXIC_MINER_DOP
    real(kind = DBL_PREC) K_ANOXIC_MINER_DOP
    real(kind = DBL_PREC) THETA_MINER_DOP
    real(kind = DBL_PREC) KHS_MINER_DOP

    !Nitrification
    real(kind = DBL_PREC) K_NITR
    real(kind = DBL_PREC) THETA_NITR
    real(kind = DBL_PREC) KHS_NITR_NH4N
    real(kind = DBL_PREC) KHS_NITR_DOXY

    !Denitrification
    real(kind = DBL_PREC) K_DENITR
    real(kind = DBL_PREC) THETA_DENITR
    real(kind = DBL_PREC) KHS_DENITR_NO3N
    real(kind = DBL_PREC) KHS_DENITR_DOC
    real(kind = DBL_PREC) KHS_DENITR_DOXY
    real(kind = DBL_PREC) DENITR_YIELD

    !Anoxia
    real(kind = DBL_PREC) DOXY_AT_ANOXIA

    !Solid partition
    real(kind = DBL_PREC) SOLID_PART_COEFF_NH4
    real(kind = DBL_PREC) SOLID_PART_COEFF_PO4

    real(kind = DBL_PREC) SED_PH_MIN_DOC_MIN
    real(kind = DBL_PREC) SED_PH_MIN_DOC_MAX
    real(kind = DBL_PREC) SED_PH_MIN_DON_MIN
    real(kind = DBL_PREC) SED_PH_MIN_DON_MAX
    real(kind = DBL_PREC) SED_PH_MIN_DOP_MIN
    real(kind = DBL_PREC) SED_PH_MIN_DOP_MAX
    real(kind = DBL_PREC) SED_PH_NITR_NH4_MIN
    real(kind = DBL_PREC) SED_PH_NITR_NH4_MAX
    real(kind = DBL_PREC) SED_PH_DENITR_NO3_MIN
    real(kind = DBL_PREC) SED_PH_DENITR_NO3_MAX

    !New model constants added in 10th of September 2015
    real(kind = DBL_PREC) :: SED_k_OX_FE_II
    real(kind = DBL_PREC) :: SED_k_RED_FE_III
    real(kind = DBL_PREC) :: SED_k_OX_MN_II
    real(kind = DBL_PREC) :: SED_k_RED_MN_IV
    real(kind = DBL_PREC) :: SED_KHS_DOXY_FE_III_RED
    real(kind = DBL_PREC) :: SED_KHS_DOXY_MN_IV_RED
    !End of new model constants added in 10th of September 2015

    !real(kind = DBL_PREC) SED_PH_KHS_DIA
    !real(kind = DBL_PREC) SED_PH_KHS_CYN
    !real(kind = DBL_PREC) SED_PH_KHS_FIX_CYN
    !real(kind = DBL_PREC) SED_PH_KHS_OPA

    ! New model constats added in 28 th of January 2016
    real(kind = DBL_PREC) :: SED_K_MIN_DOC_DOXY_20
    real(kind = DBL_PREC) :: SED_K_MIN_DOC_NO3N_20
    real(kind = DBL_PREC) :: SED_K_MIN_DOC_MN_IV_20
    real(kind = DBL_PREC) :: SED_K_MIN_DOC_FE_III_20
    real(kind = DBL_PREC) :: SED_K_MIN_DOC_S_PLUS_6_20
    real(kind = DBL_PREC) :: SED_K_MIN_DOC_DOC_20
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DOC_DOXY
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DOC_NO3N
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DOC_MN_IV
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DOC_FE_III
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DOC_S_PLUS_6
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DOC_DOC
    real(kind = DBL_PREC) :: SED_K_HS_DOC_MIN_DOXY
    real(kind = DBL_PREC) :: SED_K_HS_DOC_MIN_NO3N
    real(kind = DBL_PREC) :: SED_K_HS_DOC_MIN_MN_IV
    real(kind = DBL_PREC) :: SED_K_HS_DOC_MIN_FE_III
    real(kind = DBL_PREC) :: SED_K_HS_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: SED_K_HS_DOC_MIN_DOC
    real(kind = DBL_PREC) :: SED_K_HS_DOXY_RED_LIM
    real(kind = DBL_PREC) :: SED_K_HS_NO3N_RED_LIM
    real(kind = DBL_PREC) :: SED_K_HS_MN_IV_RED_LIM
    real(kind = DBL_PREC) :: SED_K_HS_FE_III_RED_LIM
    real(kind = DBL_PREC) :: SED_K_HS_S_PLUS_6_RED_LIM
    real(kind = DBL_PREC) :: SED_K_HS_DOXY_RED_INHB
    real(kind = DBL_PREC) :: SED_K_HS_NO3N_RED_INHB
    real(kind = DBL_PREC) :: SED_K_HS_MN_IV_RED_INHB
    real(kind = DBL_PREC) :: SED_K_HS_FE_III_RED_INHB
    real(kind = DBL_PREC) :: SED_K_HS_S_PLUS_6_RED_INHB
    real(kind = DBL_PREC) :: SED_PH_MIN_DOC_MIN_DOXY
    real(kind = DBL_PREC) :: SED_PH_MIN_DOC_MIN_NO3N
    real(kind = DBL_PREC) :: SED_PH_MIN_DOC_MIN_MN_IV
    real(kind = DBL_PREC) :: SED_PH_MIN_DOC_MIN_FE_III
    real(kind = DBL_PREC) :: SED_PH_MIN_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: SED_PH_MIN_DOC_MIN_DOC
    real(kind = DBL_PREC) :: SED_PH_MAX_DOC_MIN_DOXY
    real(kind = DBL_PREC) :: SED_PH_MAX_DOC_MIN_NO3N
    real(kind = DBL_PREC) :: SED_PH_MAX_DOC_MIN_MN_IV
    real(kind = DBL_PREC) :: SED_PH_MAX_DOC_MIN_FE_III
    real(kind = DBL_PREC) :: SED_PH_MAX_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: SED_PH_MAX_DOC_MIN_DOC

    real(kind = DBL_PREC) :: SED_K_MIN_DON_DOXY_20
    real(kind = DBL_PREC) :: SED_K_MIN_DON_NO3N_20
    real(kind = DBL_PREC) :: SED_K_MIN_DON_MN_IV_20
    real(kind = DBL_PREC) :: SED_K_MIN_DON_FE_III_20
    real(kind = DBL_PREC) :: SED_K_MIN_DON_S_PLUS_6_20
    real(kind = DBL_PREC) :: SED_K_MIN_DON_DOC_20
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DON_DOXY
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DON_NO3N
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DON_MN_IV
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DON_FE_III
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DON_S_PLUS_6
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DON_DOC
    real(kind = DBL_PREC) :: SED_K_HS_DON_MIN_DOXY
    real(kind = DBL_PREC) :: SED_K_HS_DON_MIN_NO3N
    real(kind = DBL_PREC) :: SED_K_HS_DON_MIN_MN_IV
    real(kind = DBL_PREC) :: SED_K_HS_DON_MIN_FE_III
    real(kind = DBL_PREC) :: SED_K_HS_DON_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: SED_K_HS_DON_MIN_DOC
    real(kind = DBL_PREC) :: SED_PH_MIN_DON_MIN_DOXY
    real(kind = DBL_PREC) :: SED_PH_MIN_DON_MIN_NO3N
    real(kind = DBL_PREC) :: SED_PH_MIN_DON_MIN_MN_IV
    real(kind = DBL_PREC) :: SED_PH_MIN_DON_MIN_FE_III
    real(kind = DBL_PREC) :: SED_PH_MIN_DON_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: SED_PH_MIN_DON_MIN_DOC
    real(kind = DBL_PREC) :: SED_PH_MAX_DON_MIN_DOXY
    real(kind = DBL_PREC) :: SED_PH_MAX_DON_MIN_NO3N
    real(kind = DBL_PREC) :: SED_PH_MAX_DON_MIN_MN_IV
    real(kind = DBL_PREC) :: SED_PH_MAX_DON_MIN_FE_III
    real(kind = DBL_PREC) :: SED_PH_MAX_DON_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: SED_PH_MAX_DON_MIN_DOC

    real(kind = DBL_PREC) :: SED_K_MIN_DOP_DOXY_20
    real(kind = DBL_PREC) :: SED_K_MIN_DOP_NO3N_20
    real(kind = DBL_PREC) :: SED_K_MIN_DOP_MN_IV_20
    real(kind = DBL_PREC) :: SED_K_MIN_DOP_FE_III_20
    real(kind = DBL_PREC) :: SED_K_MIN_DOP_S_PLUS_6_20
    real(kind = DBL_PREC) :: SED_K_MIN_DOP_DOC_20
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DOP_DOXY
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DOP_NO3N
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DOP_MN_IV
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DOP_FE_III
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DOP_S_PLUS_6
    real(kind = DBL_PREC) :: SED_THETA_K_MIN_DOP_DOC
    real(kind = DBL_PREC) :: SED_K_HS_DOP_MIN_DOXY
    real(kind = DBL_PREC) :: SED_K_HS_DOP_MIN_NO3N
    real(kind = DBL_PREC) :: SED_K_HS_DOP_MIN_MN_IV
    real(kind = DBL_PREC) :: SED_K_HS_DOP_MIN_FE_III
    real(kind = DBL_PREC) :: SED_K_HS_DOP_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: SED_K_HS_DOP_MIN_DOC
    real(kind = DBL_PREC) :: SED_PH_MIN_DOP_MIN_DOXY
    real(kind = DBL_PREC) :: SED_PH_MIN_DOP_MIN_NO3N
    real(kind = DBL_PREC) :: SED_PH_MIN_DOP_MIN_MN_IV
    real(kind = DBL_PREC) :: SED_PH_MIN_DOP_MIN_FE_III
    real(kind = DBL_PREC) :: SED_PH_MIN_DOP_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: SED_PH_MIN_DOP_MIN_DOC
    real(kind = DBL_PREC) :: SED_PH_MAX_DOP_MIN_DOXY
    real(kind = DBL_PREC) :: SED_PH_MAX_DOP_MIN_NO3N
    real(kind = DBL_PREC) :: SED_PH_MAX_DOP_MIN_MN_IV
    real(kind = DBL_PREC) :: SED_PH_MAX_DOP_MIN_FE_III
    real(kind = DBL_PREC) :: SED_PH_MAX_DOP_MIN_S_PLUS_6
    real(kind = DBL_PREC) :: SED_PH_MAX_DOP_MIN_DOC
    ! End of new model constats added in 28 th of January 2016

    ! New model constats added in 29 th of January 2016
    real(kind = DBL_PREC) :: SED_k_OX_CH4
    real(kind = DBL_PREC) :: SED_THETA_k_OX_CH4
    real(kind = DBL_PREC) :: SED_k_HS_OX_CH4_DOXY
    real(kind = DBL_PREC) :: SED_k_OX_H2S
    real(kind = DBL_PREC) :: SED_THETA_k_OX_H2S
    real(kind = DBL_PREC) :: SED_k_HS_OX_H2S_DOXY
    ! End of new model constats added in 29 th of January 2016

    real(kind = DBL_PREC) :: SED_k_DISS_FE_II_20
    real(kind = DBL_PREC) :: SED_THETA_k_DISS_FE_II
    real(kind = DBL_PREC) :: SED_INIT_MULT_FE_II_DISS
        
    real(kind = DBL_PREC) :: SED_k_DISS_FE_III_20
    real(kind = DBL_PREC) :: SED_THETA_k_DISS_FE_III
    real(kind = DBL_PREC) :: SED_INIT_MULT_FE_III_DISS

    
    ! END OF COEFFICIENTS
end module BSED_MODEL_CONSTANTS



subroutine INIT_BSED_MODEL_CONSTANTS
    use BSED_MODEL_CONSTANTS
    use para_aqua

    !Constants
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

    call para_get_value('SED_PH_MIN_DOC_MIN'    , SED_PH_MIN_DOC_MIN   ) !43!   optimum lower range for pH correction factor for DOC
    call para_get_value('SED_PH_MIN_DOC_MAX'    , SED_PH_MIN_DOC_MAX   ) !44!   optimum upper range for pH correction factor for DOC
    call para_get_value('SED_PH_MIN_DON_MIN'    , SED_PH_MIN_DON_MIN   ) !45!   optimum lower range for pH correction factor for DON
    call para_get_value('SED_PH_MIN_DON_MAX'    , SED_PH_MIN_DON_MAX   ) !46!   optimum upper range for pH correction factor for DON
    call para_get_value('SED_PH_MIN_DOP_MIN'    , SED_PH_MIN_DOP_MIN   ) !47!   optimum lower range for pH correction factor for DOP
    call para_get_value('SED_PH_MIN_DOP_MAX'    , SED_PH_MIN_DOP_MAX   ) !48!   optimum upper range for pH correction factor for DOP
    call para_get_value('SED_PH_NITR_NH4_MIN'   , SED_PH_NITR_NH4_MIN  ) !49!   optimum lower range for pH correction factor for nitrification
    call para_get_value('SED_PH_NITR_NH4_MAX'   , SED_PH_NITR_NH4_MAX  ) !50!   optimum upper range for pH correction factor for nitrification
    call para_get_value('SED_PH_DENITR_NO3_MIN' , SED_PH_DENITR_NO3_MIN) !51!   optimum lower range for pH correction factor for denitrification
    call para_get_value('SED_PH_DENITR_NO3_MAX' , SED_PH_DENITR_NO3_MAX) !52!   optimum upper range for pH correction factor for denitrification
    !call para_get_value('SED_PH_KHS_DIA'        , SED_PH_KHS_DIA       ) !53!   Half saturation for pH correction fator for diatoms
    !call para_get_value('SED_PH_KHS_CYN'        , SED_PH_KHS_CYN       ) !54!   Half saturation for pH correction fator for cyanobacteria not fixing
    !call para_get_value('SED_PH_KHS_FIX_CYN'    , SED_PH_KHS_FIX_CYN   ) !55!   Half saturation for pH correction fator for nitrogen fixers
    !call para_get_value('SED_PH_KHS_OPA'        , SED_PH_KHS_OPA       ) !56!   Half saturation for pH correction fator for other algae

    call para_get_value('SED_k_OX_FE_II'          , SED_k_OX_FE_II         ) !Introduced 10 09M 2015   Oxidation rate for iron 2+
    call para_get_value('SED_k_RED_FE_III'        , SED_k_RED_FE_III       ) !Introduced 10 09M 2015   reduction rate for iron 3+
    call para_get_value('SED_k_OX_MN_II'          , SED_k_OX_MN_II         ) !Introduced 10 09M 2015   oxidation rate for manganese 2+
    call para_get_value('SED_k_RED_MN_IV'         , SED_k_RED_MN_IV        ) !Introduced 10 09M 2015   reduction rate for manganese 4+
    call para_get_value('SED_KHS_DOXY_FE_III_RED' , SED_KHS_DOXY_FE_III_RED) !Introduced 10 09M 2015   reversed Monod half saturation of DOXY for iron 3+ reduction
    call para_get_value('SED_KHS_DOXY_MN_IV_RED'  , SED_KHS_DOXY_MN_IV_RED ) !Introduced 10 09M 2015   reversed Monod half saturation of DOXY for manganese 4+ reduction

    ! New model constats added in 28 th of January 2016
    call para_get_value('SED_K_MIN_DOC_DOXY_20'         , SED_K_MIN_DOC_DOXY_20        )
    call para_get_value('SED_K_MIN_DOC_NO3N_20'         , SED_K_MIN_DOC_NO3N_20        )
    call para_get_value('SED_K_MIN_DOC_MN_IV_20'        , SED_K_MIN_DOC_MN_IV_20       )
    call para_get_value('SED_K_MIN_DOC_FE_III_20'       , SED_K_MIN_DOC_FE_III_20      )
    call para_get_value('SED_K_MIN_DOC_S_PLUS_6_20'     , SED_K_MIN_DOC_S_PLUS_6_20    )
    call para_get_value('SED_K_MIN_DOC_DOC_20'          , SED_K_MIN_DOC_DOC_20         )
    call para_get_value('SED_THETA_K_MIN_DOC_DOXY'      , SED_THETA_K_MIN_DOC_DOXY     )
    call para_get_value('SED_THETA_K_MIN_DOC_NO3N'      , SED_THETA_K_MIN_DOC_NO3N     )
    call para_get_value('SED_THETA_K_MIN_DOC_MN_IV'     , SED_THETA_K_MIN_DOC_MN_IV    )
    call para_get_value('SED_THETA_K_MIN_DOC_FE_III'    , SED_THETA_K_MIN_DOC_FE_III   )
    call para_get_value('SED_THETA_K_MIN_DOC_S_PLUS_6'  , SED_THETA_K_MIN_DOC_S_PLUS_6 )
    call para_get_value('SED_THETA_K_MIN_DOC_DOC'       , SED_THETA_K_MIN_DOC_DOC      )
    call para_get_value('SED_K_HS_DOC_MIN_DOXY'         , SED_K_HS_DOC_MIN_DOXY        )
    call para_get_value('SED_K_HS_DOC_MIN_NO3N'         , SED_K_HS_DOC_MIN_NO3N        )
    call para_get_value('SED_K_HS_DOC_MIN_MN_IV'        , SED_K_HS_DOC_MIN_MN_IV       )
    call para_get_value('SED_K_HS_DOC_MIN_FE_III'       , SED_K_HS_DOC_MIN_FE_III      )
    call para_get_value('SED_K_HS_DOC_MIN_S_PLUS_6'     , SED_K_HS_DOC_MIN_S_PLUS_6    )
    call para_get_value('SED_K_HS_DOC_MIN_DOC'          , SED_K_HS_DOC_MIN_DOC         )
    call para_get_value('SED_K_HS_DOXY_RED_LIM'         , SED_K_HS_DOXY_RED_LIM        )
    call para_get_value('SED_K_HS_NO3N_RED_LIM'         , SED_K_HS_NO3N_RED_LIM        )
    call para_get_value('SED_K_HS_MN_IV_RED_LIM'        , SED_K_HS_MN_IV_RED_LIM       )
    call para_get_value('SED_K_HS_FE_III_RED_LIM'       , SED_K_HS_FE_III_RED_LIM      )
    call para_get_value('SED_K_HS_S_PLUS_6_RED_LIM'     , SED_K_HS_S_PLUS_6_RED_LIM    )
    call para_get_value('SED_K_HS_DOXY_RED_INHB'        , SED_K_HS_DOXY_RED_INHB       )
    call para_get_value('SED_K_HS_NO3N_RED_INHB'        , SED_K_HS_NO3N_RED_INHB       )
    call para_get_value('SED_K_HS_MN_IV_RED_INHB'       , SED_K_HS_MN_IV_RED_INHB      )
    call para_get_value('SED_K_HS_FE_III_RED_INHB'      , SED_K_HS_FE_III_RED_INHB     )
    call para_get_value('SED_K_HS_S_PLUS_6_RED_INHB'    , SED_K_HS_S_PLUS_6_RED_INHB   )
    call para_get_value('SED_PH_MIN_DOC_MIN_DOXY'       , SED_PH_MIN_DOC_MIN_DOXY      )
    call para_get_value('SED_PH_MIN_DOC_MIN_NO3N'       , SED_PH_MIN_DOC_MIN_NO3N      )
    call para_get_value('SED_PH_MIN_DOC_MIN_MN_IV'      , SED_PH_MIN_DOC_MIN_MN_IV     )
    call para_get_value('SED_PH_MIN_DOC_MIN_FE_III'     , SED_PH_MIN_DOC_MIN_FE_III    )
    call para_get_value('SED_PH_MIN_DOC_MIN_S_PLUS_6'   , SED_PH_MIN_DOC_MIN_S_PLUS_6  )
    call para_get_value('SED_PH_MIN_DOC_MIN_DOC'        , SED_PH_MIN_DOC_MIN_DOC       )
    call para_get_value('SED_PH_MAX_DOC_MIN_DOXY'       , SED_PH_MAX_DOC_MIN_DOXY      )
    call para_get_value('SED_PH_MAX_DOC_MIN_NO3N'       , SED_PH_MAX_DOC_MIN_NO3N      )
    call para_get_value('SED_PH_MAX_DOC_MIN_MN_IV'      , SED_PH_MAX_DOC_MIN_MN_IV     )
    call para_get_value('SED_PH_MAX_DOC_MIN_FE_III'     , SED_PH_MAX_DOC_MIN_FE_III    )
    call para_get_value('SED_PH_MAX_DOC_MIN_S_PLUS_6'   , SED_PH_MAX_DOC_MIN_S_PLUS_6  )
    call para_get_value('SED_PH_MAX_DOC_MIN_DOC'        , SED_PH_MAX_DOC_MIN_DOC       )
    call para_get_value('SED_K_MIN_DON_DOXY_20'         , SED_K_MIN_DON_DOXY_20        )
    call para_get_value('SED_K_MIN_DON_NO3N_20'         , SED_K_MIN_DON_NO3N_20        )
    call para_get_value('SED_K_MIN_DON_MN_IV_20'        , SED_K_MIN_DON_MN_IV_20       )
    call para_get_value('SED_K_MIN_DON_FE_III_20'       , SED_K_MIN_DON_FE_III_20      )
    call para_get_value('SED_K_MIN_DON_S_PLUS_6_20'     , SED_K_MIN_DON_S_PLUS_6_20    )
    call para_get_value('SED_K_MIN_DON_DOC_20'          , SED_K_MIN_DON_DOC_20         )
    call para_get_value('SED_THETA_K_MIN_DON_DOXY'      , SED_THETA_K_MIN_DON_DOXY     )
    call para_get_value('SED_THETA_K_MIN_DON_NO3N'      , SED_THETA_K_MIN_DON_NO3N     )
    call para_get_value('SED_THETA_K_MIN_DON_MN_IV'     , SED_THETA_K_MIN_DON_MN_IV    )
    call para_get_value('SED_THETA_K_MIN_DON_FE_III'    , SED_THETA_K_MIN_DON_FE_III   )
    call para_get_value('SED_THETA_K_MIN_DON_S_PLUS_6'  , SED_THETA_K_MIN_DON_S_PLUS_6 )
    call para_get_value('SED_THETA_K_MIN_DON_DOC'       , SED_THETA_K_MIN_DON_DOC      )
    call para_get_value('SED_K_HS_DON_MIN_DOXY'         , SED_K_HS_DON_MIN_DOXY        )
    call para_get_value('SED_K_HS_DON_MIN_NO3N'         , SED_K_HS_DON_MIN_NO3N        )
    call para_get_value('SED_K_HS_DON_MIN_MN_IV'        , SED_K_HS_DON_MIN_MN_IV       )
    call para_get_value('SED_K_HS_DON_MIN_FE_III'       , SED_K_HS_DON_MIN_FE_III      )
    call para_get_value('SED_K_HS_DON_MIN_S_PLUS_6'     , SED_K_HS_DON_MIN_S_PLUS_6    )
    call para_get_value('SED_K_HS_DON_MIN_DOC'          , SED_K_HS_DON_MIN_DOC         )
    call para_get_value('SED_PH_MIN_DON_MIN_DOXY'       , SED_PH_MIN_DON_MIN_DOXY      )
    call para_get_value('SED_PH_MIN_DON_MIN_NO3N'       , SED_PH_MIN_DON_MIN_NO3N      )
    call para_get_value('SED_PH_MIN_DON_MIN_MN_IV'      , SED_PH_MIN_DON_MIN_MN_IV     )
    call para_get_value('SED_PH_MIN_DON_MIN_FE_III'     , SED_PH_MIN_DON_MIN_FE_III    )
    call para_get_value('SED_PH_MIN_DON_MIN_S_PLUS_6'   , SED_PH_MIN_DON_MIN_S_PLUS_6  )
    call para_get_value('SED_PH_MIN_DON_MIN_DOC'        , SED_PH_MIN_DON_MIN_DOC       )
    call para_get_value('SED_PH_MAX_DON_MIN_DOXY'       , SED_PH_MAX_DON_MIN_DOXY      )
    call para_get_value('SED_PH_MAX_DON_MIN_NO3N'       , SED_PH_MAX_DON_MIN_NO3N      )
    call para_get_value('SED_PH_MAX_DON_MIN_MN_IV'      , SED_PH_MAX_DON_MIN_MN_IV     )
    call para_get_value('SED_PH_MAX_DON_MIN_FE_III'     , SED_PH_MAX_DON_MIN_FE_III    )
    call para_get_value('SED_PH_MAX_DON_MIN_S_PLUS_6'   , SED_PH_MAX_DON_MIN_S_PLUS_6  )
    call para_get_value('SED_PH_MAX_DON_MIN_DOC'        , SED_PH_MAX_DON_MIN_DOC       )
    call para_get_value('SED_K_MIN_DOP_DOXY_20'         , SED_K_MIN_DOP_DOXY_20        )
    call para_get_value('SED_K_MIN_DOP_NO3N_20'         , SED_K_MIN_DOP_NO3N_20        )
    call para_get_value('SED_K_MIN_DOP_MN_IV_20'        , SED_K_MIN_DOP_MN_IV_20       )
    call para_get_value('SED_K_MIN_DOP_FE_III_20'       , SED_K_MIN_DOP_FE_III_20      )
    call para_get_value('SED_K_MIN_DOP_S_PLUS_6_20'     , SED_K_MIN_DOP_S_PLUS_6_20    )
    call para_get_value('SED_K_MIN_DOP_DOC_20'          , SED_K_MIN_DOP_DOC_20         )
    call para_get_value('SED_THETA_K_MIN_DOP_DOXY'      , SED_THETA_K_MIN_DOP_DOXY     )
    call para_get_value('SED_THETA_K_MIN_DOP_NO3N'      , SED_THETA_K_MIN_DOP_NO3N     )
    call para_get_value('SED_THETA_K_MIN_DOP_MN_IV'     , SED_THETA_K_MIN_DOP_MN_IV    )
    call para_get_value('SED_THETA_K_MIN_DOP_FE_III'    , SED_THETA_K_MIN_DOP_FE_III   )
    call para_get_value('SED_THETA_K_MIN_DOP_S_PLUS_6'  , SED_THETA_K_MIN_DOP_S_PLUS_6 )
    call para_get_value('SED_THETA_K_MIN_DOP_DOC'       , SED_THETA_K_MIN_DOP_DOC      )
    call para_get_value('SED_K_HS_DOP_MIN_DOXY'         , SED_K_HS_DOP_MIN_DOXY        )
    call para_get_value('SED_K_HS_DOP_MIN_NO3N'         , SED_K_HS_DOP_MIN_NO3N        )
    call para_get_value('SED_K_HS_DOP_MIN_MN_IV'        , SED_K_HS_DOP_MIN_MN_IV       )
    call para_get_value('SED_K_HS_DOP_MIN_FE_III'       , SED_K_HS_DOP_MIN_FE_III      )
    call para_get_value('SED_K_HS_DOP_MIN_S_PLUS_6'     , SED_K_HS_DOP_MIN_S_PLUS_6    )
    call para_get_value('SED_K_HS_DOP_MIN_DOC'          , SED_K_HS_DOP_MIN_DOC         )
    call para_get_value('SED_PH_MIN_DOP_MIN_DOXY'       , SED_PH_MIN_DOP_MIN_DOXY      )
    call para_get_value('SED_PH_MIN_DOP_MIN_NO3N'       , SED_PH_MIN_DOP_MIN_NO3N      )
    call para_get_value('SED_PH_MIN_DOP_MIN_MN_IV'      , SED_PH_MIN_DOP_MIN_MN_IV     )
    call para_get_value('SED_PH_MIN_DOP_MIN_FE_III'     , SED_PH_MIN_DOP_MIN_FE_III    )
    call para_get_value('SED_PH_MIN_DOP_MIN_S_PLUS_6'   , SED_PH_MIN_DOP_MIN_S_PLUS_6  )
    call para_get_value('SED_PH_MIN_DOP_MIN_DOC'        , SED_PH_MIN_DOP_MIN_DOC       )
    call para_get_value('SED_PH_MAX_DOP_MIN_DOXY'       , SED_PH_MAX_DOP_MIN_DOXY      )
    call para_get_value('SED_PH_MAX_DOP_MIN_NO3N'       , SED_PH_MAX_DOP_MIN_NO3N      )
    call para_get_value('SED_PH_MAX_DOP_MIN_MN_IV'      , SED_PH_MAX_DOP_MIN_MN_IV     )
    call para_get_value('SED_PH_MAX_DOP_MIN_FE_III'     , SED_PH_MAX_DOP_MIN_FE_III    )
    call para_get_value('SED_PH_MAX_DOP_MIN_S_PLUS_6'   , SED_PH_MAX_DOP_MIN_S_PLUS_6  )
    call para_get_value('SED_PH_MAX_DOP_MIN_DOC'        , SED_PH_MAX_DOP_MIN_DOC       )

    ! New model constats added in 29 th of January 2016
    call para_get_value('SED_k_OX_CH4'                  , SED_k_OX_CH4          )
    call para_get_value('SED_THETA_k_OX_CH4'            , SED_THETA_k_OX_CH4    )
    call para_get_value('SED_k_HS_OX_CH4_DOXY'          , SED_k_HS_OX_CH4_DOXY  )
    call para_get_value('SED_k_OX_H2S'                  , SED_k_OX_H2S          )
    call para_get_value('SED_THETA_k_OX_H2S'            , SED_THETA_k_OX_H2S    )
    call para_get_value('SED_k_HS_OX_H2S_DOXY'          , SED_k_HS_OX_H2S_DOXY  )
    ! 2016 08 21
    call para_get_value('SED_k_DISS_FE_II_20'           , SED_k_DISS_FE_II_20      )
    call para_get_value('SED_THETA_k_DISS_FE_II'        , SED_THETA_k_DISS_FE_II   )
    call para_get_value('SED_INIT_MULT_FE_II_DISS'      , SED_INIT_MULT_FE_II_DISS )
    call para_get_value('SED_k_DISS_FE_III_20'          , SED_k_DISS_FE_III_20     )
    call para_get_value('SED_THETA_k_DISS_FE_III'       , SED_THETA_k_DISS_FE_III  )
    call para_get_value('SED_INIT_MULT_FE_III_DISS'     , SED_INIT_MULT_FE_III_DISS)

end subroutine INIT_BSED_MODEL_CONSTANTS


! Produced by Ali Erturk 2010 June
!
! Enhancement by Ali Erturk 13th of February 2011
!           - Added advection as a new transport process
!
! Enhancement by Ali Erturk and Petras Zemlys 26th of September 2014
!           - Expanded the model by adding inorganic carbon, alkality and
!             salinity as new state variables and pH as a derived variable
!
! Vectorization by Petras Zemlys 2014 December
!
! Change by Ali Erturk 2nd of February 2015
!           - Carbon related state variables are subroutinized
!
! Change by Ali Erturk and Petras Zemlys 10 th of September 2015
!           - Introduction of four new state variables
!
!                 FE_II  : Iron 2+
!                 FE_III : Iron 3+
!                 MN_II  : Manganese 2+
!                 MN_IV  : Manganese 2+
!
! Change by Ali Erturk and Petras Zemlys 30th of November 2015
!
!          Initial additions for interfacing dissolved and particulate
!          fractions of FE_II, FE_IV, MN_II, MN_IV
!
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
!              - CH4_C     (Methane carbon)
!*******************************************************************

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
            SED_BURRIAL_RATE_OUTPUTS, &
            ADVANCED_REDOX_OPTION)

    use CO2SYS_CDIAC
    use AQUABC_II_GLOBAL
    use BSED_MODEL_CONSTANTS
    !use para_aqua

    !use basin, only: ipv !0d correction
    
    ! Flags: 1,2,3 = 3,4,5 WC 
    

    implicit none

    !include 'param.h'


!     integer ipv(nkndim)	!external node numbers
!     common  /ipv/ipv

    !ARGUMENTS RELATED TO ARAY SIZES

    !NUMBER OF CONSTANTS FOR CONTROL
    integer, parameter :: CHECK_NUM_SED_LAYERS = 1
    integer, parameter :: SINGLE_VECTOR_CO2SYS = 0

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
    !PARAMETER(NUM_SED_OUTPUTS = 14)

    !COUNTERS
    integer i, j,k,l
    integer CELLNO ! Variable for node numbers


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
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS, NUM_SED_SAVED_OUTPUTS) :: SED_SAVED_OUTPUTS
    integer, dimension(NUM_SED_FLAGS) :: SED_FLAGS
    real(kind=DBL_PREC), optional  , dimension(nkn, NUM_SED_LAYERS, NUM_SED_VARS) :: SED_BURRIAL_RATE_OUTPUTS
    integer            , intent(in), optional :: ADVANCED_REDOX_OPTION
    
    !SEDIMENT STATE VARIABLES(Mass per volume of sediments, g/m3)
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
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: SALT      !Salinity
    ! End of new state variables added 27 September 2014

    ! New state variables added 10 September 2015
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: FE_II     ! Iron 2+
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: FE_III    ! Iron 3+
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: MN_II     ! Manganese 2+
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: MN_IV     ! Manganese 4+
    ! End of new state variables added 10 September 2015

    ! New state variables added 27 January 2016
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: CA
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: MG
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: S_PLUS_6
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: S_MINUS_2
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: CH4_C
    ! End of new state variables added 27 January 2016

    !DERIVED VARIABLES
    real(kind=DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: pH

    !SEDIMENT KINETIC PROCESS RATES
    double precision, dimension(nkn,NUM_SED_LAYERS) :: R_DISS_POC
    double precision, dimension(nkn,NUM_SED_LAYERS) :: R_DISS_PON
    double precision, dimension(nkn,NUM_SED_LAYERS) :: R_DISS_POP
    double precision, dimension(nkn,NUM_SED_LAYERS) :: R_MINER_DOC
    double precision, dimension(nkn,NUM_SED_LAYERS) :: R_MINER_DON
    double precision, dimension(nkn,NUM_SED_LAYERS) :: R_MINER_DOP
    double precision, dimension(nkn,NUM_SED_LAYERS) :: R_NITR
    double precision, dimension(nkn,NUM_SED_LAYERS) :: R_DENITR
    double precision, dimension(nkn,NUM_SED_LAYERS) :: R_DISS_PSi
    double precision, dimension(nkn,NUM_SED_LAYERS) :: DEOXYGENATION

    ! New model kinetic process rates added in 10th of September 2015
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: R_FE_II_OXIDATION
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: R_FE_III_REDUCTION
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: R_MN_II_OXIDATION
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: R_MN_IV_REDUCTION
    ! End of new model kinetic process rates added in 10th of September 2015

    ! New kinetic processes added in 28 th of January 2016
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DOC_DOXY
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DOC_NO3N
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DOC_MN_IV
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DOC_FE_III
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DOC_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DOC_DOC

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DON_DOXY
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DON_NO3N
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DON_MN_IV
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DON_FE_III
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DON_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DON_DOC

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DOP_DOXY
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DOP_NO3N
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DOP_MN_IV
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DOP_FE_III
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DOP_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_MINER_DOP_DOC
    ! End of new kinetic processes added in 28 th of January 2016

    ! New kinetic rates introduced 29 January 2016 for the redox sequences
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_SULPHATE_REDUCTION
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_SULPHIDE_OXIDATION
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_METHANOGENESIS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: R_METHANE_OXIDATION
    ! End of new kinetic rates introduced 29 January 2016 for the redox sequences

    !SEDIMENT TRANSPORT PROCESS RATES
    double precision NEIGHBOUR_CONC     (nkn)
    double precision UPPER_CONC_GRADIENT(nkn)
    double precision SED_MIXLEN         (nkn)
    double precision DIFF_CORRECTION_FACTOR(nkn) !diffusion correction factor for porosity
    double precision SED_IN_ADVEC_RATES (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision ADV_ENTERING_CONC(nkn)
    double precision SED_OUT_ADVEC_RATES(nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision SED_DIFFUSION_RATES(nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision PART_MIXING_RATES  (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision SED_BURRIAL_RATES  (nkn,NUM_SED_LAYERS, NUM_SED_VARS)

    ! UNIT AREA MASSES FOR STATE VARIABLES
    double precision UNIT_AREA_MASSES   (nkn,NUM_SED_LAYERS, NUM_SED_VARS)

    !DERIVS
    double precision DERIVS             (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision TRANSPORT_DERIVS   (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision ADVECTION_DERIVS   (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision DIFFUSION_DERIVS   (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision BURIAL_DERIVS      (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision PART_MIXING_DERIVS (nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision SETTLING_DERIVS    (nkn,NUM_SED_LAYERS, NUM_SED_VARS) ! when sedtrans not used
    double precision DEPOSIT_CONC       (nkn,NUM_SED_VARS) ! when sedtrans is used!
    double precision RESUSP_FLUX        (nkn,NUM_SED_VARS) ! when sedtrans is used
    double precision KINETIC_DERIVS     (nkn,NUM_SED_LAYERS, NUM_SED_VARS)

    double precision SETTLING_AFFECTED_DEPTH(nkn)
    double precision SINK_INTENSITY, DIFF_SINK(nkn) ! for difussion boundary sink gradient definition
    double precision DIFF_DRAG


    ! debugging variables and functions
    integer STRANGER   !Function checking for strange values not vectorised
    integer STRANGERSD !Function checking for strange values vectorised
    integer error     !Error indicator
    logical VALUE_strange(nkn) ! For NaN and Inf checking
    ! For strange values processing
    real(kind = DBL_PREC),allocatable, dimension (:) :: STRANGERS     ! Strange values
    real(kind = DBL_PREC),allocatable, dimension (:) :: STRANGERS1
    real(kind = DBL_PREC),allocatable, dimension (:) :: STRANGERS2

    integer              ,allocatable, dimension (:) :: NODES_STRANGE ! node numbers with strange values
    !integer :: index_strange (nkn)                                    ! array with values 1:nkn
    integer :: nstrange                                               ! number of nodes with strange values
    logical debug_stranger                                            ! strangers debugging switcher
    ! end debugging variables


    !VARIABLES FOR INTERMEDIATE INTEGRATING
    integer TIME_LOOP
    integer NUM_SUB_TIME_STEPS
    double precision INTERMED_RESULTS (nkn,NUM_SED_LAYERS, NUM_SED_VARS) ! state vars
    double precision INTERMED_RESULTS_(nkn,NUM_SED_LAYERS, NUM_SED_VARS) ! state vars (temporary)

    ! FOR PROCESSING AND CONVERSION OF DIFFERENT PHASES
    double precision WATER_DENSITY
    integer IN_WHICH_PHASE(NUM_SED_VARS)
    double precision VOLUME_FRACTION(nkn,NUM_SED_LAYERS, NUM_SED_VARS) ! Fraction of unit volume occupied
                                                                       ! by state variable depending on its phase(phases)
    double precision SOLUTE_FRACTIONS(nkn,NUM_SED_LAYERS, NUM_SED_VARS)
    double precision SOLID_CONCS(nkn,NUM_SED_LAYERS)

    ! VARIABLES FOR SWITCH ON/OFF PROCESSES
    INTEGER switch_kinetics
    INTEGER switch_burial
    INTEGER switch_partmixing
    INTEGER switch_diffusion
    INTEGER switch_advection
    INTEGER switch_settling !  for old version when sedtrans is not used

    INTEGER switch_erosion    ! for new version
    INTEGER switch_deposition ! for new version


    !27 September 2014
    !New variables for DIC and ALK
    real(kind = DBL_PREC), allocatable, dimension (:) :: CO2SYS_PAR1
    real(kind = DBL_PREC), allocatable, dimension (:) :: CO2SYS_PAR2
    integer, allocatable , dimension (:)              :: CO2SYS_PAR1TYPE
    integer, allocatable , dimension (:)              :: CO2SYS_PAR2TYPE
    real(kind = DBL_PREC),allocatable, dimension (:)  :: CO2SYS_SALT
    real(kind = DBL_PREC),allocatable, dimension (:)  :: CO2SYS_TEMPIN
    real(kind = DBL_PREC),allocatable, dimension (:)  :: CO2SYS_TEMPOUT
    real(kind = DBL_PREC),allocatable, dimension (:)  :: CO2SYS_PRESIN
    real(kind = DBL_PREC),allocatable, dimension (:)  :: CO2SYS_PRESOUT
    real(kind = DBL_PREC),allocatable, dimension (:)  :: CO2SYS_SI
    real(kind = DBL_PREC),allocatable, dimension (:)  :: CO2SYS_PO4
    integer, allocatable, dimension (:) :: CO2SYS_pHSCALEIN
    integer, allocatable, dimension (:) :: CO2SYS_K1K2CONSTANTS
    integer, allocatable, dimension (:) :: CO2SYS_KSO4CONSTANTS
    real(kind = DBL_PREC),allocatable, dimension (:,:) :: CO2SYS_OUT_DATA
    character(len=34), allocatable, dimension (:) :: CO2SYS_NICEHEADERS
    integer :: CO2SYS_NUM_SAMPLES
    integer :: RUN_CO2SYS
    parameter(RUN_CO2SYS = 1)
    integer :: CO2SYS_ntps

    !16 Februar 2017
    real(kind = DBL_PREC), allocatable, dimension (:) :: ALL_LAYERS_CO2SYS_PAR1
    real(kind = DBL_PREC), allocatable, dimension (:) :: ALL_LAYERS_CO2SYS_PAR2
    integer, allocatable , dimension (:)              :: ALL_LAYERS_CO2SYS_PAR1TYPE
    integer, allocatable , dimension (:)              :: ALL_LAYERS_CO2SYS_PAR2TYPE
    real(kind = DBL_PREC),allocatable, dimension (:)  :: ALL_LAYERS_CO2SYS_SALT
    real(kind = DBL_PREC),allocatable, dimension (:)  :: ALL_LAYERS_CO2SYS_TEMPIN
    real(kind = DBL_PREC),allocatable, dimension (:)  :: ALL_LAYERS_CO2SYS_TEMPOUT
    real(kind = DBL_PREC),allocatable, dimension (:)  :: ALL_LAYERS_CO2SYS_PRESIN
    real(kind = DBL_PREC),allocatable, dimension (:)  :: ALL_LAYERS_CO2SYS_PRESOUT
    real(kind = DBL_PREC),allocatable, dimension (:)  :: ALL_LAYERS_CO2SYS_SI
    real(kind = DBL_PREC),allocatable, dimension (:)  :: ALL_LAYERS_CO2SYS_PO4
    integer, allocatable, dimension (:) :: ALL_LAYERS_CO2SYS_pHSCALEIN
    integer, allocatable, dimension (:) :: ALL_LAYERS_CO2SYS_K1K2CONSTANTS
    integer, allocatable, dimension (:) :: ALL_LAYERS_CO2SYS_KSO4CONSTANTS
    real(kind = DBL_PREC),allocatable, dimension (:,:) :: ALL_LAYERS_CO2SYS_OUT_DATA
    !character(len=34), allocatable, dimension (:) :: ALL_LAYERS_CO2SYS_NICEHEADERS
    
    
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

    real getpar    !function to getn isedi value
    integer isedi  !indicator if sediment transport model is used

    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: PH_CORR_DOC_MIN
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: PH_CORR_DON_MIN
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: PH_CORR_DOP_MIN
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: PH_CORR_NITR_NH4
    real(kind = DBL_PREC), dimension(nkn,NUM_SED_LAYERS) :: PH_CORR_DENITR_NO3

    ! -------------------------------------------------------------------------
    ! Variables added 30 November 2015 for dissolved and particulate species
    ! of FE_II, FE_III, MN_II, MN_IV. These will be internally calculated by
    ! this subroutine and will be returned in some suitable way to AQUABC since
    ! they will be used in settlin calculations. Job by Petras and his fellows.
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MULT_FE_II_DISS  !Dissolved fraction of FE_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MULT_FE_II_PART  !Particulate fraction of FE_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MULT_FE_III_DISS !Dissolved fraction of FE_III
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MULT_FE_III_PART !Particulate fraction of FE_III

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MULT_MN_II_DISS  !Dissolved fraction of MN_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MULT_MN_II_PART  !Particulate fraction of MN_II
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MULT_MN_IV_DISS  !Dissolved fraction of MN_IV
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MULT_MN_IV_PART  !Particulate fraction of MN_IV
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

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_EQ_S_1
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_EQ_S_2
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_SP_FES
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: HS2_TOT
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: H2S_DIVISOR
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FRAC_H2S_IN_H2S_TOT
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FRAC_HS_MINUS_IN_H2S_TOT
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FRAC_S_MINUS_TWO_IN_H2S_TOT
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: H2S
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: HS_MINUS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: S_MINUS_TWO
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_II_DISS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_II_PART
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_III_DISS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_III_PART
    ! -------------------------------------------------------------------------
    ! End of variables added 25 January 2016 for a simple version of equlibrium
    ! based aquatic chemistry calculations for Fe2+ and Fe3+
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! New auxillary added in 28 th of January 2016
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DOC_MIN_DOXY
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DOC_MIN_NO3N
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DOC_MIN_MN_IV
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DOC_MIN_FE_III
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DOC_MIN_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DOC_MIN_DOC
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: LIM_DOXY_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: LIM_NO3N_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: LIM_MN_IV_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: LIM_FE_III_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: LIM_S_PLUS_6_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: LIM_DOC_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_NO3_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_MN_IV_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_FE_III_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_S_PLUS_6_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: K_DOC_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MN_IV_DISS

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DON_MIN_DOXY
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DON_MIN_NO3N
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DON_MIN_MN_IV
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DON_MIN_FE_III
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DON_MIN_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DON_MIN_DOC

    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DOP_MIN_DOXY
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DOP_MIN_NO3N
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DOP_MIN_MN_IV
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DOP_MIN_FE_III
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DOP_MIN_S_PLUS_6
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PH_CORR_DOP_MIN_DOC
    ! -------------------------------------------------------------------------
    ! End of new auxillary added in 28 th of January 2016
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! New variables introduced by Ali, 3 July 2016
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: HCO3
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: CO3    
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: SED_LIM_DOXY_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: SED_LIM_NO3N_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: SED_LIM_MN_IV_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: SED_LIM_FE_III_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: SED_LIM_S_PLUS_6_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: SED_LIM_DOC_RED
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: PE
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: MN_II_DISS
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: ELEVATION
    ! -------------------------------------------------------------------------
    ! End of new variables introduced by Ali, 3 July 2016
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! New variables added (21 August 2016)
    ! -------------------------------------------------------------------------
    ! FE_II_DISS_EQ               : Dissolved Fe2+ in equilibrium (solubility of Fe2+)
    ! FE_III_DISS_EQ              : Dissolved Fe2+ in equilibrium (solubility of Fe3+)
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_II_DISS_EQ
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: FE_III_DISS_EQ
    ! -------------------------------------------------------------------------
    ! End of new variables added (21 August 2016)
    ! -------------------------------------------------------------------------
    
    ! -------------------------------------------------------------------------
    ! New auxillary variables introduced 21 August 2016
    ! -------------------------------------------------------------------------
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: DISS_FE_II_CONC_TS_END
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: DISS_FE_II_CONC_TS_AVG
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: DISS_FE_III_CONC_TS_END
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: DISS_FE_III_CONC_TS_AVG
    real(kind = DBL_PREC), dimension(nkn, NUM_SED_LAYERS) :: DIP_OVER_IP
    ! -------------------------------------------------------------------------
    ! New auxillary variables introduced 21 August 2016
    ! -------------------------------------------------------------------------    
    
    ! Introduced by Ali
    integer :: FIRST_TIME_STEP
    integer :: INIT_OPTION_OF_FE_II_DISS
    integer :: INIT_OPTION_OF_FE_III_DISS

    integer  :: DO_ADVANCED_REDOX_SIMULATION
    integer  :: NUM_SIMULATED_SED_VARS
    
    !**************************************************************************
    !*                                                                        *
    !*                            EXECUTION PART                              *
    !*                                                                        *
    !**************************************************************************
    FIRST_TIME_STEP            = SED_FLAGS(1)
    INIT_OPTION_OF_FE_II_DISS  = SED_FLAGS(2)
    INIT_OPTION_OF_FE_III_DISS = SED_FLAGS(3)

    DO_ADVANCED_REDOX_SIMULATION = 0
    NUM_SIMULATED_SED_VARS       = 15
    
    if (present(ADVANCED_REDOX_OPTION)) then
        if (ADVANCED_REDOX_OPTION > 0) then
            DO_ADVANCED_REDOX_SIMULATION = 1
            NUM_SIMULATED_SED_VARS       = 24
        end if
    end if

    !Checking state and output dimensions. Change them if they are changed
    if (CHECK_NUM_SED_LAYERS == 1) then
        if(NUM_SED_LAYERS .ne. 7) then
           print *, 'SEDIMENT_MODEL_1: Wrong number of layers'
           stop
        end if
    end if

    if(NUM_SED_VARS .ne. 24) then
        print *, 'SEDIMENT_MODEL_1: Wrong number of state variables'
        stop
    end if

    if(NUM_SED_OUTPUTS .ne. 26) then
        print *, 'SEDIMENT_MODEL_1: Wrong number of outputs'
        stop
    end if

    !------------- EXECUTION CONTROLS--------------------
    isedi = nint(getpar('isedi')) !indicator if sediment transport model is used

    !  Switching derivatives:
    !     0 - switched off (is calculated but makes equal zero)
    !     1 - is calculated
    switch_kinetics   = 1
    switch_burial     = 1
    switch_partmixing = 0
    switch_diffusion  = 1
    switch_advection  = 0
    switch_settling   = 1 ! for old version, when sedtrans not used
    switch_erosion    = 1 ! for new version
    switch_deposition = 1 ! for new version

    CONSIDER_ALKALNITY_DERIVATIVE = 1
    CONSIDER_INORG_C_DERIVATIVE   = 1
    CONSIDER_CO2_REARATION        = 1 ! not used yet

    error = 0                     !indicator of presence of stranger9
    debug_stranger = .true.
    !debug_stranger = .false.     !True if check for strangers
    !-----------------------------------------------------

    ! Multiplier for diffusion rate for the first layer (used reverse for negative)
    DIFF_DRAG = 6.0D0 !fixme


    !INITIALIZATIONS

    !sediment kinetic process rates
    R_DISS_POC   (1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    R_DISS_PON   (1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    R_DISS_POP   (1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    R_MINER_DOC  (1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    R_MINER_DON  (1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    R_MINER_DOP  (1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    R_NITR       (1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    R_DENITR     (1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    R_DISS_PSi   (1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    DEOXYGENATION(1:nkn,1:NUM_SED_LAYERS) = 0.0D0

    ! Initialize New model kinetic process rates added in 10th of September 2015
    R_FE_II_OXIDATION (1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    R_FE_III_REDUCTION(1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    R_MN_II_OXIDATION (1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    R_MN_IV_REDUCTION (1:nkn,1:NUM_SED_LAYERS) = 0.0D0
    ! End of initialize New model kinetic process rates added in 10th of September 2015

    !sediment transport process rates
    NEIGHBOUR_CONC        (1:nkn) = 0.0D0
    DIFF_CORRECTION_FACTOR(1:nkn) = 0.0D0 !diffusion correction factor for porosity
    UPPER_CONC_GRADIENT   (1:nkn) = 0.0D0
    SED_MIXLEN            (1:nkn) = 0.0D0
    SED_IN_ADVEC_RATES    (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
    SED_OUT_ADVEC_RATES   (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
    SED_DIFFUSION_RATES   (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
    PART_MIXING_RATES     (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
    SED_BURRIAL_RATES     (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0

    SETTLING_AFFECTED_DEPTH(1:nkn) = 0.0D0
    ! for difussion boundary sink gradient definition:
    SINK_INTENSITY     = 0.0D0
    DIFF_SINK(1:nkn)   = 0.0D0

    !DERIVS
    DERIVS             (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
    TRANSPORT_DERIVS   (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
    ADVECTION_DERIVS   (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
    DIFFUSION_DERIVS   (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
    BURIAL_DERIVS      (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
    PART_MIXING_DERIVS (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
    SETTLING_DERIVS    (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
    DEPOSIT_CONC       (1:nkn, 1:NUM_SED_VARS)                   = 0.0D0
    RESUSP_FLUX        (1:nkn, 1:NUM_SED_VARS)                   = 0.0D0
    KINETIC_DERIVS     (1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0


    ! ITERMEDIATE INTEGRATION
    TIME_LOOP          = 0
    NUM_SUB_TIME_STEPS = 0
    INTERMED_RESULTS(1:nkn,1:NUM_SED_LAYERS , 1:NUM_SED_VARS) = 0.0D0
    INTERMED_RESULTS_(1:nkn,1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0

    ! unit area masses for state variables
    UNIT_AREA_MASSES(1:nkn, 1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
    VOLUME_FRACTION (1:nkn,1:NUM_SED_LAYERS , 1:NUM_SED_VARS) = 0.0D0 ! for calculation of unit area masses

    ! solute phase fraction of variable in concentration per total sediment volume
    SOLUTE_FRACTIONS(1:nkn,1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
                                       !for variables NH4 and PO4 recalculated later.

    ! dry bulk density (dry mass of particles per total volume), kg/L:
    SOLID_CONCS(1:nkn,1:NUM_SED_LAYERS) = 0.0D0

    PROCESSES_sed(1:nkn,1:NUM_SED_LAYERS, 1:NUM_SED_VARS, 1:NDIAGVAR_sed) = 0.0D+0

    PART_MIXING_RATES(1:nkn,1:NUM_SED_LAYERS, 1:NUM_SED_VARS) = 0.0D0
                                       ! Made zero in order not to have particle mixing
                                       ! It is also cancelled by switch. fixme

    !END OF INITIALIZATIONS

    !**************************************************************************
    !*                                                                        *
    !*                            CALCULATION PART                            *
    !*                                                                        *
    !**************************************************************************

    ! Calculation bulk dry density:
    WATER_DENSITY = 1.0D0

    SOLID_CONCS(1:nkn,1:NUM_SED_LAYERS) = &
        (SED_DENSITIES(1:nkn,1:NUM_SED_LAYERS) - &
        (WATER_DENSITY * SED_POROSITIES(1:nkn,1:NUM_SED_LAYERS)))

    !Initialisation of solute fractions is moved to initialisation section

    !IN_WHICH_PHASE
    !0 : SOLUTE
    !1 : SOLID
    !2 : ALL SEDIMENTS

    if(NUM_SED_VARS .ne. 24) then
        print *, 'SEDIMENT_MODEL: Number of state variables is wrong', NUM_SED_VARS
        stop
    end if

    IN_WHICH_PHASE(1)  = 2 !2  !Amonia
    IN_WHICH_PHASE(2)  = 0     !Nitrates
    IN_WHICH_PHASE(3)  = 0     !DON
    IN_WHICH_PHASE(4)  = 1 !1  !PON
    IN_WHICH_PHASE(5)  = 2 !2  !PO4P
    IN_WHICH_PHASE(6)  = 0     !DOP
    IN_WHICH_PHASE(7)  = 1 !1  !POP
    IN_WHICH_PHASE(8)  = 0     !DOXY
    IN_WHICH_PHASE(9)  = 0     !DOC
    IN_WHICH_PHASE(10) = 1 !1  !POC
    IN_WHICH_PHASE(11) = 0     !DSi
    IN_WHICH_PHASE(12) = 1 !1  !PSi
    IN_WHICH_PHASE(13) = 0     !DIC
    IN_WHICH_PHASE(14) = 0     !Alkalinity
    IN_WHICH_PHASE(15) = 0     !Salinity
    IN_WHICH_PHASE(16) = 2     !FE_II  : Introduced 10th of September 2015  !Changed to two after the changes in 30th of November
    IN_WHICH_PHASE(17) = 2     !FE_III : Introduced 10th of September 2015  !Changed to two after the changes in 30th of November
    IN_WHICH_PHASE(18) = 2     !MN_II  : Introduced 10th of September 2015  !Changed to two after the changes in 30th of November
    IN_WHICH_PHASE(19) = 2     !MN_IV  : Introduced 10th of September 2015  !Changed to two after the changes in 30th of November
    IN_WHICH_PHASE(20) = 1     ! For a while assume that calcium is solid. To be corrected soon
    IN_WHICH_PHASE(21) = 1     ! For a while assume that magnesium is solid. To be corrected soon
    IN_WHICH_PHASE(22) = 0
    IN_WHICH_PHASE(23) = 0
    IN_WHICH_PHASE(24) = 0

    if(debug_stranger) then
        do i = 1, NUM_SED_LAYERS
            do j = 1, NUM_SIMULATED_SED_VARS
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

    !NUM_FLUXES_TO_SEDIMENTS = NUM_SED_VARS
    if(debug_stranger) then
        do j = 1, NUM_SIMULATED_SED_VARS

            if(STRANGERSD(FLUXES_TO_SEDIMENTS(:,j),VALUE_strange,nkn) .eq. 1) then
                nstrange = count(VALUE_strange)
                allocate(STRANGERS    (nstrange))
                allocate(NODES_STRANGE(nstrange))

                l=1
                do k=1,nkn
                    if(VALUE_strange(k)) then
                        STRANGERS    (l) = FLUXES_TO_SEDIMENTS(k,j)
                        NODES_STRANGE(l) = k
                        l=l+1
                    end if
                end do

                print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print *, 'aquabc_sediment1:  Variable ',j
                print *, 'Flux from WC is NaN'
                print *, 'NODE_NUMBERS=',NODES_STRANGE
                print *, 'VALUES=',STRANGERS

                print *, ' '

                deallocate(STRANGERS, NODES_STRANGE)
                error =1
            end if !STRANGERSD

        end do ! j-variables
    end if ! debug_stranger

    if (error .eq. 1) stop

    if(NUM_SED_VARS .ne. 24) then
        print *, 'SEDIMENT_MODEL: Number of state variables is wrong', NUM_SED_VARS
        stop
    end if

    !INITIALIZE SEDIMENT STATE VARIABLES
    SED_NH4N (:,:) = INIT_SED_STATE_VARS(:,:, 1)
    SED_NO3N (:,:) = INIT_SED_STATE_VARS(:,:, 2)
    SED_DON  (:,:) = INIT_SED_STATE_VARS(:,:, 3)
    SED_PON  (:,:) = INIT_SED_STATE_VARS(:,:, 4)
    SED_PO4P (:,:) = INIT_SED_STATE_VARS(:,:, 5)
    SED_DOP  (:,:) = INIT_SED_STATE_VARS(:,:, 6)
    SED_POP  (:,:) = INIT_SED_STATE_VARS(:,:, 7)
    SED_DOXY (:,:) = INIT_SED_STATE_VARS(:,:, 8)
    SED_DOC  (:,:) = INIT_SED_STATE_VARS(:,:, 9)
    SED_POC  (:,:) = INIT_SED_STATE_VARS(:,:, 10)
    SED_DSi  (:,:) = INIT_SED_STATE_VARS(:,:, 11)
    SED_PSi  (:,:) = INIT_SED_STATE_VARS(:,:, 12)
    INORG_C  (:,:) = INIT_SED_STATE_VARS(:,:, 13)
    TOT_ALK  (:,:) = INIT_SED_STATE_VARS(:,:, 14)
    SALT     (:,:) = INIT_SED_STATE_VARS(:,:, 15)
    FE_II    (:,:) = INIT_SED_STATE_VARS(:,:, 16) !Introduced 10th of September 2015
    FE_III   (:,:) = INIT_SED_STATE_VARS(:,:, 17) !Introduced 10th of September 2015
    MN_II    (:,:) = INIT_SED_STATE_VARS(:,:, 18) !Introduced 10th of September 2015
    MN_IV    (:,:) = INIT_SED_STATE_VARS(:,:, 19) !Introduced 10th of September 2015
    CA       (:,:) = INIT_SED_STATE_VARS(:,:, 20) !Introduced 27th of January 2016
    MG       (:,:) = INIT_SED_STATE_VARS(:,:, 21) !Introduced 27th of January 2016
    S_PLUS_6 (:,:) = INIT_SED_STATE_VARS(:,:, 22) !Introduced 27th of January 2016
    S_MINUS_2(:,:) = INIT_SED_STATE_VARS(:,:, 23) !Introduced 27th of January 2016
    CH4_C    (:,:) = INIT_SED_STATE_VARS(:,:, 24) !Introduced 27th of January 2016


    !INITIALIZE SEDIMENT MODEL COEFFICIENTS
    if(NUM_SED_CONSTS .ne. 171) then
        print *, 'SEDIMENT_MODEL_1: Wrong number of constants'
    end if


    !OLD(INITIAL) STATE VARIABLES ASIGNED TO ARRAY 'INTERMED_RESULTS'
    INTERMED_RESULTS(:,:, :) = INIT_SED_STATE_VARS(:,:, :)

    NUM_SUB_TIME_STEPS = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!LOOP ON SUBSTEPS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !(Substeps are not used. For correct use averaging shoud be implemented)


    ! -----------------------------------------------------------------
    ! Additions made 30th November 2015 related to dissolved and solid
    ! iron and mangese species, soon to be replaced by the
    ! aquatic chemistry model that will be called
    !
    !        EXACTLY HERE
    !
    ! However, be aware that this submodel may need CO2SYS to be run
    ! before it is callded because of the pH and Alkanlinity calculation
    ! needs.
    !-----------------------------------------------------------------


    ! -------------------------------------------------------------------------
    ! New equlibrium calculations added at 25th of January 2016
    ! This is the first version of the aquatic chemistry model that was
    ! promissed 30th November 2015
    ! -------------------------------------------------------------------------

    ! Use AQUABC_II_GLOBAL simplified version of H2CO3 dissosication to estimate
    ! HCO3- and CO3--

    if (SINGLE_VECTOR_CO2SYS > 0) then
        allocate(ALL_LAYERS_CO2SYS_PAR1         (nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_PAR2         (nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_PAR1TYPE     (nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_PAR2TYPE     (nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_SALT         (nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_TEMPIN       (nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_TEMPOUT      (nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_PRESIN       (nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_PRESOUT      (nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_SI           (nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_PO4          (nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_pHSCALEIN    (nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_K1K2CONSTANTS(nkn * NUM_SED_LAYERS), &
                 ALL_LAYERS_CO2SYS_KSO4CONSTANTS(nkn * NUM_SED_LAYERS))
    
            do i = 1, NUM_SED_LAYERS
                ALL_LAYERS_CO2SYS_PAR1  ((((i - 1) * nkn) + 1):(i * nkn)) = TOT_ALK  (:,i) * 1.0D6
                ALL_LAYERS_CO2SYS_PAR2  ((((i - 1) * nkn) + 1):(i * nkn)) = INORG_C  (:,i) * 1.0D6
                ALL_LAYERS_CO2SYS_SALT  ((((i - 1) * nkn) + 1):(i * nkn)) = SALT     (:,i)
                ALL_LAYERS_CO2SYS_TEMPIN((((i - 1) * nkn) + 1):(i * nkn)) = SED_TEMPS(:,i)
                ALL_LAYERS_CO2SYS_SI    ((((i - 1) * nkn) + 1):(i * nkn)) = (SED_DSI (:,i) / 28.0855D0) * 1.0D3
                ALL_LAYERS_CO2SYS_PO4   ((((i - 1) * nkn) + 1):(i * nkn)) = (SED_PO4P(:,i) / 30.9737D0) * 1.0D3
            end do
            
            ALL_LAYERS_CO2SYS_PAR1TYPE     (:) = 1
            ALL_LAYERS_CO2SYS_PAR2TYPE     (:) = 2
            ALL_LAYERS_CO2SYS_pHSCALEIN    (:) = 1
            ALL_LAYERS_CO2SYS_K1K2CONSTANTS(:) = 4
            ALL_LAYERS_CO2SYS_KSO4CONSTANTS(:) = 1
            ALL_LAYERS_CO2SYS_TEMPOUT      (:) = 0.0D0  !Does not matter for this example
            ALL_LAYERS_CO2SYS_PRESIN       (:) = 0.0D0  !Does not matter for this example
            ALL_LAYERS_CO2SYS_PRESOUT      (:) = 0.0D0  !Does not matter for this example
    
            CO2SYS_ntps = nkn * NUM_SED_LAYERS
    
            call CO2SYS(ALL_LAYERS_CO2SYS_PAR1         , ALL_LAYERS_CO2SYS_PAR2         , ALL_LAYERS_CO2SYS_PAR1TYPE , &
                        ALL_LAYERS_CO2SYS_PAR2TYPE     , ALL_LAYERS_CO2SYS_SALT         , ALL_LAYERS_CO2SYS_TEMPIN   , &
                        ALL_LAYERS_CO2SYS_TEMPOUT      , ALL_LAYERS_CO2SYS_PRESIN       , ALL_LAYERS_CO2SYS_PRESOUT  , &
                        ALL_LAYERS_CO2SYS_SI           , ALL_LAYERS_CO2SYS_PO4          , ALL_LAYERS_CO2SYS_pHSCALEIN, &
                        ALL_LAYERS_CO2SYS_K1K2CONSTANTS, ALL_LAYERS_CO2SYS_KSO4CONSTANTS, ALL_LAYERS_CO2SYS_OUT_DATA , &
                        CO2SYS_NICEHEADERS             , CO2SYS_ntps)

            do i = 1, NUM_SED_LAYERS
                !LOWER_INDEX = (((i - 1) * nkn) + 1)
                !UPPER_INDEX =  (i * nkn)
                pH         (:,i) = ALL_LAYERS_CO2SYS_OUT_DATA((((i - 1) * nkn) + 1):(i * nkn), 18)
                K_ONE_TIP  (:,i) = ALL_LAYERS_CO2SYS_OUT_DATA((((i - 1) * nkn) + 1):(i * nkn), 75)
                K_TWO_TIP  (:,i) = ALL_LAYERS_CO2SYS_OUT_DATA((((i - 1) * nkn) + 1):(i * nkn), 76)
                K_THREE_TIP(:,i) = ALL_LAYERS_CO2SYS_OUT_DATA((((i - 1) * nkn) + 1):(i * nkn), 77)
                HCO3       (:,i) = ALL_LAYERS_CO2SYS_OUT_DATA((((i - 1) * nkn) + 1):(i * nkn), 21)
                CO3        (:,i) = ALL_LAYERS_CO2SYS_OUT_DATA((((i - 1) * nkn) + 1):(i * nkn), 22)
            end do

            if(allocated(ALL_LAYERS_CO2SYS_PAR1         )) deallocate(ALL_LAYERS_CO2SYS_PAR1         )
            if(allocated(ALL_LAYERS_CO2SYS_PAR2         )) deallocate(ALL_LAYERS_CO2SYS_PAR2         )
            if(allocated(ALL_LAYERS_CO2SYS_PAR1TYPE     )) deallocate(ALL_LAYERS_CO2SYS_PAR1TYPE     )
            if(allocated(ALL_LAYERS_CO2SYS_PAR2TYPE     )) deallocate(ALL_LAYERS_CO2SYS_PAR2TYPE     )
            if(allocated(ALL_LAYERS_CO2SYS_SALT         )) deallocate(ALL_LAYERS_CO2SYS_SALT         )
            if(allocated(ALL_LAYERS_CO2SYS_TEMPIN       )) deallocate(ALL_LAYERS_CO2SYS_TEMPIN       )
            if(allocated(ALL_LAYERS_CO2SYS_TEMPOUT      )) deallocate(ALL_LAYERS_CO2SYS_TEMPOUT      )
            if(allocated(ALL_LAYERS_CO2SYS_PRESIN       )) deallocate(ALL_LAYERS_CO2SYS_PRESIN       )
            if(allocated(ALL_LAYERS_CO2SYS_PRESOUT      )) deallocate(ALL_LAYERS_CO2SYS_PRESOUT      )
            if(allocated(ALL_LAYERS_CO2SYS_SI           )) deallocate(ALL_LAYERS_CO2SYS_SI           )
            if(allocated(ALL_LAYERS_CO2SYS_PO4          )) deallocate(ALL_LAYERS_CO2SYS_PO4          )
            if(allocated(ALL_LAYERS_CO2SYS_pHSCALEIN    )) deallocate(ALL_LAYERS_CO2SYS_pHSCALEIN    )
            if(allocated(ALL_LAYERS_CO2SYS_K1K2CONSTANTS)) deallocate(ALL_LAYERS_CO2SYS_K1K2CONSTANTS)
            if(allocated(ALL_LAYERS_CO2SYS_KSO4CONSTANTS)) deallocate(ALL_LAYERS_CO2SYS_KSO4CONSTANTS)
            if(allocated(ALL_LAYERS_CO2SYS_OUT_DATA     )) deallocate(ALL_LAYERS_CO2SYS_OUT_DATA     )
            if(allocated(CO2SYS_NICEHEADERS             )) deallocate(CO2SYS_NICEHEADERS             )
    else
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
        do i = 1, NUM_SED_LAYERS
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
            HCO3       (:,i) = CO2SYS_OUT_DATA(:, 21) / 1.0D6 !Convert from micromoles to moles
            CO3        (:,i) = CO2SYS_OUT_DATA(:, 22) / 1.0D6 !Convert from micromoles to moles

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
        end do !i = 1,NUM_SED_LAYERS    
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
    end if

    H_PLUS(:,:) = 10.0D0 ** (-PH)
        
    ! Introduced 3 rd of July 2016
    ELEVATION = 0.0D0
    
    if (DO_ADVANCED_REDOX_SIMULATION > 0) then

        call SED_REDOX_AND_SPECIATION &
            (SED_DOXY, SED_NO3N, MN_IV, FE_III, S_PLUS_6, SED_DOC, &
             S_MINUS_2 , MN_II, FE_II , HCO3, CO3, &
             SED_TEMPS, SALT, PH, ELEVATION, &
             SED_K_HS_DOXY_RED_LIM    , SED_K_HS_NO3N_RED_LIM     , &
             SED_K_HS_MN_IV_RED_LIM   , SED_K_HS_FE_III_RED_LIM   , &
             SED_K_HS_S_PLUS_6_RED_LIM, SED_K_HS_DOXY_RED_INHB    , &
             SED_K_HS_NO3N_RED_INHB   , SED_K_HS_MN_IV_RED_INHB   , &
             SED_K_HS_FE_III_RED_INHB , SED_K_HS_S_PLUS_6_RED_INHB, &
             nkn, NUM_SED_LAYERS      , &
             SED_LIM_DOXY_RED         , SED_LIM_NO3N_RED          , &
             SED_LIM_MN_IV_RED        , SED_LIM_FE_III_RED        , &
             SED_LIM_S_PLUS_6_RED     , SED_LIM_DOC_RED           , &
             PE, FE_II_DISS_EQ        , FE_III_DISS_EQ, MN_II_DISS)
        
        FE_III_DISS_EQ = FE_III_DISS_EQ * 56000.0D0
             
        ! For now, take fixed values for equlibrium and solubility constants
        ! In future, change these to temperature and/or other environmental
        ! variable based equations.
        K_EQ_S_1  (:,:) = 8.90D-8
        K_EQ_S_2  (:,:) = 1.20D-13
        K_SP_FES  (:,:) = 8.00D-19
        !K_SP_FEOH3(:,:) = 1.00D-36
        !K_SP_FEPO4(:,:) = 9.91D-16
        
        ! -------------------------------------------------------------------------
        ! Calculate H2S Species
        ! -------------------------------------------------------------------------
        
        ! Convert sulphide to moles for aquatic chemistry calculations
        HS2_TOT(:,:) = S_MINUS_2(:,:) / 32000.0D0
        
        H2S_DIVISOR(:,:) = &
            (H_PLUS(:,:) * H_PLUS(:,:)) + &
            (H_PLUS(:,:) * K_EQ_S_1(:,:)) + (K_EQ_S_1(:,:) * K_EQ_S_2(:,:))
        
        FRAC_HS_MINUS_IN_H2S_TOT   (:,:) = &
            (H_PLUS(:,:)   * K_EQ_S_1(:,:)) / H2S_DIVISOR(:,:)
        
        FRAC_S_MINUS_TWO_IN_H2S_TOT(:,:) = &
            (K_EQ_S_1(:,:) * K_EQ_S_2(:,:)) / H2S_DIVISOR(:,:)
        
        FRAC_H2S_IN_H2S_TOT(:,:)         = &
            1.0D0 - (FRAC_HS_MINUS_IN_H2S_TOT(:,:) + FRAC_S_MINUS_TWO_IN_H2S_TOT(:,:))
        
        H2S        (:,:) = HS2_TOT(:,:) * FRAC_H2S_IN_H2S_TOT        (:,:)
        HS_MINUS   (:,:) = HS2_TOT(:,:) * FRAC_HS_MINUS_IN_H2S_TOT   (:,:)
        S_MINUS_TWO(:,:) = HS2_TOT(:,:) * FRAC_S_MINUS_TWO_IN_H2S_TOT(:,:)
        ! -------------------------------------------------------------------------
        ! End of calculate H2S Species
        ! -------------------------------------------------------------------------
        
        ! -------------------------------------------------------------------------
        ! Calculate the dissolved and particulate fractions of Fe2+
        ! -------------------------------------------------------------------------
        
        ! -------------------------------------------------------------------------
        ! Changes in 5th of July 2016
        !
        !           The following correction proposal has been conducted
        ! -------------------------------------------------------------------------
        
        ! Assume FeS precipitation as key process. To be corrected with better
        ! algorithm in future
        
        ! By commenting the following code lines in the code from Spring 2016 visit    
        
        ! Calculate how much Fe2+ from FeS could be dissolved for equlibrium
        !where (S_MINUS_TWO > 1.0D-12)
        !    FE_II_DISS = (K_SP_FES / S_MINUS_TWO) * 56000.0D0
        !elsewhere
        !    FE_II_DISS = FE_II
        !end where
        
        ! -------------------------------------------------------------------------
        ! End of changes in 5th of July 2016
        ! -------------------------------------------------------------------------
        
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
        
        ! -------------------------------------------------------------------------
        ! Changes by Ali Ertürk 5 July 2016)
        ! -------------------------------------------------------------------------
        call SED_IRON_II_DISSOLUTION &
             (HS2_TOT, PH, TOT_ALK, nkn, NUM_SED_LAYERS, FE_II_DISS_EQ)
        
        FE_II_DISS_EQ = FE_II_DISS_EQ * 56000.0D0
        ! -------------------------------------------------------------------------
        ! End of changes by Ali Ertürk (5 July 2016)
        ! -------------------------------------------------------------------------
        
        
        ! -------------------------------------------------------------------------
        ! Updated by Ali and Petras, 21 August 2016
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
                    ! If more Fe2+ is allowed to dissolve than the total Fe2+ present 
                    ! then
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
                    MULT_FE_II_DISS = SED_INIT_MULT_FE_II_DISS
                    MULT_FE_II_PART = 1.0D0 - MULT_FE_II_DISS
        
            end select
        
            call CALC_DISS_ME_CONC &
                 (FE_II                                                                 , & ! Total Fe2+
                  (MULT_FE_II_DISS * FE_II)                                             , & ! Dissolved Fe2+ from previous time step
                  FE_II_DISS_EQ                                                         , & ! Equilibrium concentration for dissolved Fe2+
                  (SED_k_DISS_FE_II_20 * (SED_THETA_k_DISS_FE_II**(SED_TEMPS - 20.0D0))), & ! Dissolution rate constant for Fe2+
                  TIME_STEP                                                             , & ! Time step in days
                  nkn                                                                   , &
                  NUM_SED_LAYERS                                                        , & ! number of layers
                  DISS_FE_II_CONC_TS_END                                                , & ! Estimated dissolved Fe2+ at the end of time step (for output)
                  DISS_FE_II_CONC_TS_AVG)                                                  ! Estimated avg. dissolved Fe2+ during the timestep to be used for kinetic calculations
        else
            call CALC_DISS_ME_CONC &
                 (FE_II                                                                 , & ! Total Fe2+
                  (SED_SAVED_OUTPUTS(:,:,1) * FE_II)                                    , & ! Dissolved Fe2+ from previous time step
                  FE_II_DISS_EQ                                                         , & ! Equilibrium concentration for dissolved Fe2+
                  (SED_k_DISS_FE_II_20 * (SED_THETA_k_DISS_FE_II**(SED_TEMPS - 20.0D0))), & ! Dissolution rate constant for Fe2+
                  TIME_STEP                                                             , & ! Time step in days
                  nkn                                                                   , &
                  NUM_SED_LAYERS                                                        , & ! number of layers
                  DISS_FE_II_CONC_TS_END                                                , & ! Estimated dissolved Fe2+ at the end of time step (for output)
                  DISS_FE_II_CONC_TS_AVG)                                                   ! Estimated avg. dissolved Fe2+ during the timestep to be used for kinetic calculations
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
        !PROCESS_RATES(1:nkn,26, 5) = FE_III_DISS_EQ
        !PROCESS_RATES(1:nkn,26, 6) = FE_III
        
        where (FE_III_DISS_EQ > FE_III)
            FE_III_DISS_EQ = FE_III
        end where
        
        !PROCESS_RATES(1:nkn,26, 7) = FE_III_DISS_EQ
        
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
                    MULT_FE_III_DISS = SED_INIT_MULT_FE_III_DISS
                    MULT_FE_III_PART = 1.0D0 - MULT_FE_III_DISS
        
            end select
        
            call CALC_DISS_ME_CONC &
                 (FE_III                                                                  , & ! Total Fe3+
                  (MULT_FE_III_DISS * FE_III)                                             , & ! Dissolved Fe3+ from previous time step
                  FE_III_DISS_EQ                                                          , & ! Equilibrium concentration for dissolved Fe3+
                  (SED_k_DISS_FE_III_20 * (SED_THETA_k_DISS_FE_III**(SED_TEMPS - 20.0D0))), & ! Dissolution rate constant for Fe3+
                  TIME_STEP                                                               , & ! Time step in days
                  nkn                                                                     , &
                  NUM_SED_LAYERS                                                          , & ! number of layers
                  DISS_FE_III_CONC_TS_END                                                 , & ! Estimated dissolved Fe3+ at the end of time step (for output)
                  DISS_FE_III_CONC_TS_AVG)                                                    ! Estimated avg. dissolved Fe3+ during the timestep to be used for kinetic calculations
        else
            call CALC_DISS_ME_CONC &
                 (FE_III                                                                  , & ! Total Fe3+
                  (SED_SAVED_OUTPUTS(:,:,2) * FE_III)                                     , & ! Dissolved Fe3+ from previous time step
                  FE_III_DISS_EQ                                                          , & ! Equilibrium concentration for dissolved Fe3+
                  (SED_k_DISS_FE_III_20 * (SED_THETA_k_DISS_FE_III**(SED_TEMPS - 20.0D0))), & ! Dissolution rate constant for Fe3+
                  TIME_STEP                                                               , & ! Time step in days
                  nkn                                                                     , &
                  NUM_SED_LAYERS                                                          , & ! number of layers
                  DISS_FE_III_CONC_TS_END                                                 , & ! Estimated dissolved Fe3+ at the end of time step (for output)
                  DISS_FE_III_CONC_TS_AVG)                                                    ! Estimated avg. dissolved Fe3+ during the timestep to be used for kinetic calculations
        end if                                                                            
        
        !PROCESS_RATES(1:nkn,26, 8) = DISS_FE_III_CONC_TS_AVG
        !PROCESS_RATES(1:nkn,26, 9) = DISS_FE_III_CONC_TS_END
        
            
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
        
        ! -------------------------------------------------------------------------
        ! End of new equlibrium calculations added at 25th of January 2016
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
            
        MULT_MN_IV_DISS(:,:)  = 0.0D0
        MULT_MN_IV_PART(:,:)  = 1.0D0
        ! -----------------------------------------------------------------------
        ! End of additions made 30th November 2015 related to dissolved and solid
        ! iron and mangese species
        ! -----------------------------------------------------------------------
        
        if(debug_stranger) then                        
            do j=1,num_sed_layers
                do k= 1,nkn
                    if (STRANGER(MULT_FE_II_DISS(k,j)).eq.1) then
                        print *, 'MULT_FE_II_DISS is strange'
                        print *, 'FIRST_TIME_STEP            : ', FIRST_TIME_STEP
                        print *, 'CELL NO                    : ', k
                        print *, 'LAYER NO                   : ', j
                        print *, 'MULT_FE_II_DISS(k,j)       : ', MULT_FE_II_DISS(k,j)                                    
                        print *, 'FE_II(k,j)                 : ', FE_II(k,j)
                        print *, 'INIT_OPTION_OF_FE_II_DISS  : ', INIT_OPTION_OF_FE_II_DISS
                        print *, 'SED_INIT_MULT_FE_II_DISS   : ', SED_INIT_MULT_FE_II_DISS
                        print *, 'FE_II_DISS_EQ              : ', FE_II_DISS_EQ(k,j) 
                        error = 1
                    end if
        
                    if (STRANGER(MULT_FE_III_DISS(k,j)).eq.1) then
                        print *, '----------------------------------------------------------'
                        print *, 'MULT_FE_III_DISS is strange'
                        print *, '----------------------------------------------------------'
                        print *, 'FIRST_TIME_STEP              : ', FIRST_TIME_STEP
                        print *, 'CELL NO                      : ', k
                        print *, 'LAYER NO                     : ', j
                        print *, 'MULT_FE_III_DISS(k,j)        : ', MULT_FE_III_DISS(k,j)                                   
                        print *, 'FE_III(k,j)                  : ', FE_III(k,j)
                        print *, 'INIT_OPTION_OF_FE_III_DISS   : ', INIT_OPTION_OF_FE_III_DISS
                        print *, 'SED_INIT_MULT_FE_III_DISS    : ', SED_INIT_MULT_FE_III_DISS
                        print *, 'FE_III_DISS_EQ               : ', FE_III_DISS_EQ(k,j) 
                        print *, 'DISS_FE_III_CONC_TS_END(k,j) : ', DISS_FE_III_CONC_TS_END(k,j)
                        print *, 'DISS_FE_III_CONC_TS_AVG(k,j) : ', DISS_FE_III_CONC_TS_AVG(k,j)
                        print *, 'SED_k_DISS_FE_III_20         : ', SED_k_DISS_FE_III_20 
                        print *, 'SED_THETA_k_DISS_FE_III      : ', SED_THETA_k_DISS_FE_III
                        print *, 'SED_TEMPS(k,j)               : ', SED_TEMPS(k,j)
                        print *, '----------------------------------------------------------'
                        error = 1
                    end if
                end do
            end do                        
            
	    	if(error .eq. 1) stop
        end if !debug_stranger
        
        SED_SAVED_OUTPUTS(:,:,1) = DISS_FE_II_CONC_TS_END  / FE_II
        SED_SAVED_OUTPUTS(:,:,2) = DISS_FE_III_CONC_TS_END / FE_III
        SED_SAVED_OUTPUTS(:,:,3) = MULT_MN_II_DISS(:,:)
        SED_SAVED_OUTPUTS(:,:,4) = MULT_MN_IV_DISS(:,:)
        
        MN_IV_DISS(:,:) = MN_IV(:,:) * MULT_MN_IV_DISS(:,:)
        
        call IP_SOLUBLE_FRACTION &
             (FE_III        , &
              SED_PO4P      , &
              K_ONE_TIP     , &
              K_TWO_TIP     , &
              K_THREE_TIP   , &
              PH            , &
              nkn           , &
              NUM_SED_LAYERS, &
              DIP_OVER_IP)
        
        SED_SAVED_OUTPUTS(:,:,5) = DIP_OVER_IP(:,:)
    else

        SED_LIM_DOXY_RED = SED_DOXY  / (SED_DOXY + SED_K_HS_DOXY_RED_LIM)

        SED_LIM_NO3N_RED = &
            (SED_NO3N / (SED_NO3N + SED_K_HS_NO3N_RED_LIM)) * &
            (SED_K_HS_DOXY_RED_LIM / (SED_DOXY + SED_K_HS_DOXY_RED_LIM))
    
        DIP_OVER_IP = 1.0D0
            
        SED_SAVED_OUTPUTS(:,:,1) = 0.0D0
        SED_SAVED_OUTPUTS(:,:,2) = 0.0D0
        SED_SAVED_OUTPUTS(:,:,3) = 0.0D0
        SED_SAVED_OUTPUTS(:,:,4) = 0.0D0
        SED_SAVED_OUTPUTS(:,:,5) = DIP_OVER_IP
    end if


    ! Derivatives loop    
    do TIME_LOOP = 1, NUM_SUB_TIME_STEPS

        !UNIT AREA MASSES (for integration, g/m2) and SOLUTE_FRACTIONS in
        ! total(solutes+solids) concentrations per total unit volume of sediments
        do j = 1, NUM_SIMULATED_SED_VARS

            if (IN_WHICH_PHASE(j).eq.0) then
                SOLUTE_FRACTIONS(:,:,j) = 1.0D0
                VOLUME_FRACTION (:,:,j) = SED_POROSITIES(:,:)
            end if

            if (IN_WHICH_PHASE(j).eq.1) then
                SOLUTE_FRACTIONS(:,:,j) = 0.0D0
                VOLUME_FRACTION (:,:,j) = 1.0D0 - SED_POROSITIES(:,:)
            end if

            if (IN_WHICH_PHASE(j).eq.2) then
                      !Calculate solute fractions of NH4 and PO4
                if (j .eq. 1 ) &
                    SOLUTE_FRACTIONS(:,:, j) = & !(1.0D0 / SED_POROSITIES(:,:)) * &
                          (1.0D0 / (1.0D0 + (SOLID_CONCS(:,:) * SOLID_PART_COEFF_NH4))) ! fraction of solute

                    !fixed bug: division by porosity will give the concentrations in pore water
                    ! while for other solute phase variables concentrations per
                    ! sediment total unit volume were used. Division by porosities is done now
                    ! when necessary to have concentrations per pore water or per solids?
                if (j .eq. 5 ) &
                    SOLUTE_FRACTIONS(:,:, j) = & !(1.0D0 / SED_POROSITIES(:,:)) * &
                           (1.0D0 / (1.0D0 + (SOLID_CONCS(:,:) * SOLID_PART_COEFF_PO4)))
                
                if (j .eq. 16) &
                    SOLUTE_FRACTIONS(:,:, j) =  MULT_FE_II_DISS(:,:)

                if (j .eq. 17) &
                    SOLUTE_FRACTIONS(:,:, j) =  MULT_FE_III_DISS(:,:)

                if (j .eq. 18) &
                    SOLUTE_FRACTIONS(:,:, j) = MULT_MN_II_DISS(:,:)

                if (j .eq. 19) &
                    SOLUTE_FRACTIONS(:,:, j) = MULT_MN_IV_DISS(:,:)

                VOLUME_FRACTION(:,:,j) = 1.0D0
            end if
          
            !g/m3(~mg/l) * m =  g/m2
            UNIT_AREA_MASSES(:,:, j) = INTERMED_RESULTS(:,:, j) * SED_DEPTHS(:,:)
        end do !j = 1, NUM_SED_VARS
        !END UNIT AREA MASSES
        
        if(debug_stranger) then
            do i=1,num_sed_vars
                do j=1,num_sed_layers
                    do k= 1,nkn
                        if (STRANGER(SOLUTE_FRACTIONS(k,j,i)).eq.1) then
                            print *, 'SOLUTE_FRACTIONS is strange'
                            print *, 'CELL NO                  : ', k
                            print *, 'LAYER NO                 : ', j
                            print *, 'STATE VARIABLE NO        : ', i
                            print *, 'SOLUTE_FRACTIONS(k,j,i)  :', SOLUTE_FRACTIONS(k,j,i)
                            print *, 'MULT_FE_II_DISS(k,j)     :', MULT_FE_II_DISS(k,j)
                            print *, 'MULT_FE_III_DISS(k,j)    :', MULT_FE_III_DISS(k,j)
                            print *, 'SOLID_CONCS(k,j)         :', SOLID_CONCS(k,j)
                            print *, 'SED_POROSITIES(k,j)      :', SED_POROSITIES(k,j)
                            error = 1
                        end if
                    end do
                end do
            end do
            
			if(error .eq. 1) stop
        end if !debug_stranger

        !++++++++++++++++++++++++++PREPARATION FOR CALCULATION DERIVATIVES++++++++++++++++++++++++++++++++++

        do i=1,NUM_SED_LAYERS
            do j = 1, NUM_SIMULATED_SED_VARS

                ! PREPARATION FOR DifFERENT PHASES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !For phase 0(solute) and 2(solute and solid phase)
                if ((IN_WHICH_PHASE(j).eq.0).or.(IN_WHICH_PHASE(j).eq.2)) then !***************************

                    if (i.eq.1) then  ! For the first layer --------------------------------

                        !ADVECTION
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_advection.ne.0) then

                            if (ADVECTIVE_VELOCITY <= 0.0D0) then
                                ADV_ENTERING_CONC(:) = SURF_WATER_CONCS(:,j)
                            else
                           
                                ! division by porosity to have concentration per pore water
                                ADV_ENTERING_CONC(:) = &
                                    (INTERMED_RESULTS(:,i + 1, j) * &
                                     SOLUTE_FRACTIONS(:,i + 1, j)) / &
                                SED_POROSITIES(:,i+1)         
                            end if

                        end if !switch_advection.ne.0

                        !DIFFUSION
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_diffusion.ne.0) then
                            NEIGHBOUR_CONC(:) = SURF_WATER_CONCS(:,j)
                            SED_MIXLEN        = SURF_MIXLEN
                        end if

                    else !for  layer numbers  > 1 -----------------------------------------

                       !ADVECTION ENTERING CONC.
                       !Entering concentration is from lower to upper layer
                       !in direction of the flow

                       !24/06/2012: If clause Added by Ali for safer switch off
                       if (switch_advection.ne.0) then

                           if (ADVECTIVE_VELOCITY <= 0.0D0) then

                               ! division by porosity to have concentration per pore water
                               ADV_ENTERING_CONC(:) = &
                                   (INTERMED_RESULTS(:,i - 1, j) * &
                                    SOLUTE_FRACTIONS(:,i - 1, j)) / &
                               SED_POROSITIES(:,i-1)
                           else
                               if (i.eq.NUM_SED_LAYERS) then !for the last layer
                                   !Let nothing comes from below (corrected by Petras)
                                    ADV_ENTERING_CONC = 0.0D0
                                    !      INTERMED_RESULTS(i, j) *
                                    !      SOLUTE_FRACTIONS(i, j)
                                else !for the middle layers
                                
                                    ! Division by porosity to have concentration 
                                    ! per pore water volume
                                    ADV_ENTERING_CONC(:) =  &
                                        (INTERMED_RESULTS(:,i + 1, j) * &
                                         SOLUTE_FRACTIONS(:,i + 1, j)) / &
                                    SED_POROSITIES(:,i+1)
                                                               
                                end if ! for the last layer
                            end if !ADVECTIVE_VELOCITY <= 0.0D0
                        end if !switch_advection.ne.0

                        !DIFFUSION
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_diffusion.ne.0) then
                            NEIGHBOUR_CONC(:) = &
                                INTERMED_RESULTS(:,i - 1, j) * SOLUTE_FRACTIONS(:,i - 1, j)/SED_POROSITIES(:,i-1)
                                                             ! division by porosity to have concentration per pore water

                            SED_MIXLEN(:)     = 0.5D0 * (SED_DEPTHS(:,i - 1) + SED_DEPTHS(:,i))
                        end if !switch_diffusion.ne.0
                    end if !for  layer numbers  > 1 ------------------------------------------------

                    !ADVECTION IN AND OUT RATES
                    !24/06/2012: If clause Added by Ali for safer switch off
                    if (switch_advection.ne.0) then
                        SED_IN_ADVEC_RATES(:,i, j) = &
                            ADV_ENTERING_CONC(:) * dabs(ADVECTIVE_VELOCITY)

                        SED_OUT_ADVEC_RATES(:,i, j) = INTERMED_RESULTS(:,i, j) * &
                            SOLUTE_FRACTIONS(:,i, j)/SED_POROSITIES(:,i) * DABS(ADVECTIVE_VELOCITY)
                                                 ! division by porosity to have concentration per pore water
                    end if

                    !DIFFUSION
                    !24/06/2012: If clause Added by Ali for safer switch off
                    if (switch_diffusion.ne.0) then
                        UPPER_CONC_GRADIENT(:) = (INTERMED_RESULTS(:,i, j) * &
                               SOLUTE_FRACTIONS(:,i, j)/SED_POROSITIES(:,i)) - NEIGHBOUR_CONC(:)
                                                  ! division by porosity to have concentration per pore water
                        if(debug_stranger) then
                            do k= 1,nkn
                                if (STRANGER(UPPER_CONC_GRADIENT(k)).eq.1) then
                                    print *, 'UPPER_CONC_GRADIENT for diffusion is strange'
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
                        end if !debug_stranger

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
                    end if !switch_diffusion.ne.0

                    !FLUX FROM SEDIMENTS TO WATER COLUMN
                    if (i.eq.1) then
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_diffusion.ne.0) then
                            FLUXES_FROM_SEDIMENTS(:,j) = SED_DIFFUSION_RATES(:,i, j)
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
                    end if !i.eq.1
                    !END FLUX FROM SEDIMENTS TO WATER COLUMN

                    !FOR PARTICLE MIXING
                    if (IN_WHICH_PHASE(j).eq.0) then
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_partmixing.ne.0) then
                            PART_MIXING_RATES(:,i, j) = 0.0D0
                        end if
                    end if !IN_WHICH_PHASE(j).eq.0

                    if (IN_WHICH_PHASE(j).eq.2) then
                        if (i .eq. 1) then
                            !24/06/2012: If clause Added by Ali for safer switch off
                            if (switch_partmixing.ne.0) then
                                PART_MIXING_RATES(:,i, j) = 0.0D0 !no mixing for the first layer
                            end if
                        else ! i >1
                            !24/06/2012: If clause Added by Ali for safer switch off
                            if (switch_partmixing.ne.0) then
                                NEIGHBOUR_CONC(:) = &
                                    INTERMED_RESULTS(:,i - 1, j) * (1.0D0 - SOLUTE_FRACTIONS(:,i, j))

                                SED_MIXLEN(:) = 0.5D0 * (SED_DEPTHS(:,i - 1) + SED_DEPTHS(:,i))

                                UPPER_CONC_GRADIENT(:) = &
                                    (INTERMED_RESULTS(:,i, j) * (1.0D0 - SOLUTE_FRACTIONS(:,i, j))) - &
                                    NEIGHBOUR_CONC(:) !Result is in concentration per total sediment volume

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
                                end if !debug_stranger

                                PART_MIXING_RATES(:,i, j) = &
                                    (UPPER_CONC_GRADIENT * PART_MIXING_COEFFS(:,i, j)) / SED_MIXLEN
                                     !Result is in concentration per total sediment unit volume
                            end if !switch_partmixing.ne.0
                        end if  ! i >1
                    end if !IN_WHICH_PHASE(j).eq.2
                else !IN_WHICH_PHASE(j).eq. 1: ********************************************
                    FLUXES_FROM_SEDIMENTS (:,j) = 0.0D0
                    SED_DIFFUSION_RATES(:,i, j) = 0.0D0

                    if (i.eq.1) then
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_partmixing.ne.0) then
                            PART_MIXING_RATES(:,i, j) = 0.0D0
                        end if
                    else ! i > 1
                        !24/06/2012: If clause Added by Ali for safer switch off
                        if (switch_partmixing.ne.0) then
                            NEIGHBOUR_CONC(:) = INTERMED_RESULTS(:,i - 1, j)* (1-SOLUTE_FRACTIONS(:,i-1, j))
                            SED_MIXLEN(:)     = 0.5D0 * (SED_DEPTHS(:,I - 1) + SED_DEPTHS(:,i))

                            UPPER_CONC_GRADIENT(:) = &
                                  (INTERMED_RESULTS(:,i, j) * (1-SOLUTE_FRACTIONS(:,i, j))) - NEIGHBOUR_CONC(:)

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
                                  !Result is in concentration per total sediment unit volume
                        end if !switch_partmixing
                    end if ! i > 1
                end if !IN_WHICH_PHASE(j).eq. 1 !END PREPARATION FOR DifFERENT PHASES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !24/06/2012: If clause Added by Ali for safer switch off
                if (switch_burial.ne.0) then
                    ! Code for Elin process outputs
                    if (present(SED_BURRIAL_RATE_OUTPUTS)) then
                        SED_BURRIAL_RATE_OUTPUTS(:,i, j) = INTERMED_RESULTS(:,i, j) * SED_BURRIALS(:,i)
                    end if
                    ! End of code for Elin process outputs
                    
                    ! Bug fixed: INTERMED_RESULTS(:,i, j) * SED_BURRIALS(:,i) produces a unit of
                    ! g/m^2/day. So we should divide it by the sediment depth to get the derivative
                    ! in units of g/m^3/day

                    SED_BURRIAL_RATES(:,i, j) = & !Result is in concentration per total sediment unit volume
                        (INTERMED_RESULTS(:,i, j) * SED_BURRIALS(:,i)) / (SED_DEPTHS(:, i))
                end if
            end do !j = 1, NUM_SED_VARS
        end do     !i=1,NUM_SED_LAYERS

        !++++++++++++++++++++++++++END OF PREPARATION FOR CALCULATION DERIVATIVES++++++++++++++++++++++++++

        !------------------------------TRANSPORT DERIVATIVES-----------------------------------------------

        ! First layers for transport
        do i = 1, (NUM_SED_LAYERS - 1)
            do j = 1, NUM_SIMULATED_SED_VARS

                if (i .gt. 1) then
                    !Introduced by Petras for diagnostic output
                    DIFFUSION_DERIVS  (:,i, j) = &
                        SED_DIFFUSION_RATES(:,i + 1, j) - SED_DIFFUSION_RATES(:,i, j)

                    BURIAL_DERIVS     (:,i, j) = &
                        SED_BURRIAL_RATES  (:,i - 1, j) - SED_BURRIAL_RATES  (:,i, j)

                    PART_MIXING_DERIVS(:,i, j) = &
                        PART_MIXING_RATES  (:,i + 1, j) - PART_MIXING_RATES  (:,i, j)

                    ! Introduced by Petras to avoid double calculation
                    !diffusion derivs multipl. by poros. to have deriv. in total volume
                    TRANSPORT_DERIVS(:,i, j) = &
                            DIFFUSION_DERIVS  (:,i, j)*SED_POROSITIES(:,i) + & 
                            BURIAL_DERIVS     (:,i, j) + &
                            PART_MIXING_DERIVS(:,i, j)
                else ! i=1

                    !Introduced by Petras for diagnostic output
                    DIFFUSION_DERIVS  (:,i, j) =  &
                        SED_DIFFUSION_RATES(:,i + 1, j) - SED_DIFFUSION_RATES(:,i, j)

                    BURIAL_DERIVS     (:,i, j) = (-1.0D0)*SED_BURRIAL_RATES  (:,i, j)

                    PART_MIXING_DERIVS(:,i, j) = &
                        PART_MIXING_RATES  (:,i + 1, j) - PART_MIXING_RATES  (:,i, j)

                    ! Introduced by Petras to avoid double calculation. It should be done everywhere fixme
                    TRANSPORT_DERIVS(:,i, j) = &
                            DIFFUSION_DERIVS  (:,i, j) * SED_POROSITIES(:,i) + &    !diffusion derivs multipl.
                            BURIAL_DERIVS     (:,i, j) + PART_MIXING_DERIVS(:,i, j) !by poros. to have deriv. in total volume
                                                         
                end if

                ADVECTION_DERIVS(:,i, j) = &
                    SED_IN_ADVEC_RATES(:,i, j) - SED_OUT_ADVEC_RATES(:,i, j)
                
                !Advection derivs multipl. by poros. to have deriv. in total volume
                TRANSPORT_DERIVS(:,i, j) = TRANSPORT_DERIVS(:,i, j) + ADVECTION_DERIVS(:,i, j)* SED_POROSITIES(:,i)

            end do !j = 1, NUM_SED_VARS
        end do !i = 1, (NUM_SED_LAYERS - 1)
        !End first layers for transport

        !Last layer processing (why mixing is not processed here?)
        do j = 1, NUM_SED_VARS
            if (NUM_SED_LAYERS .GT. 1) then
                SINK_INTENSITY = 1.D0 !fixme should not be hardcoded

                DIFF_SINK(:) = (-1)*SINK_INTENSITY*  DABS(SED_DIFFUSION_RATES(:,NUM_SED_LAYERS, j)) !Lower boundary gradient
                DIFFUSION_DERIVS(:,NUM_SED_LAYERS, j) = &
                             DIFF_SINK(:) - SED_DIFFUSION_RATES(:,NUM_SED_LAYERS, j)

                BURIAL_DERIVS(:,NUM_SED_LAYERS, j) = &
                   SED_BURRIAL_RATES (:,NUM_SED_LAYERS - 1, j) - &
                   SED_BURRIAL_RATES (:,NUM_SED_LAYERS, j)

                PART_MIXING_DERIVS(:,NUM_SED_LAYERS, j) = 0. !temporary fi

                ADVECTION_DERIVS(:,NUM_SED_LAYERS, j) = &   !Added by Petras(last layer was not processed for advection)
                   SED_IN_ADVEC_RATES (:,NUM_SED_LAYERS, j) - &
                   SED_OUT_ADVEC_RATES(:,NUM_SED_LAYERS, j)

                TRANSPORT_DERIVS(:,NUM_SED_LAYERS, j) = &
                    DIFFUSION_DERIVS  (:,NUM_SED_LAYERS, j) * SED_POROSITIES(:,NUM_SED_LAYERS) + & !diffusion derivs multipl. by poros. to have deriv. in total volume &
                    BURIAL_DERIVS(:,NUM_SED_LAYERS, j) + &
                    ADVECTION_DERIVS(:,NUM_SED_LAYERS, j) * SED_POROSITIES(:,NUM_SED_LAYERS) + & !advection derivs multipl. by poros. to have deriv. in total volume & &
                    PART_MIXING_DERIVS(:,NUM_SED_LAYERS, j)
            else !NUM_SED_LAYERS.EQ. 1
                !Introduced by Petras for testing:
                DIFFUSION_DERIVS(:,NUM_SED_LAYERS, j) = &
                    ((-1.0D0) * SED_DIFFUSION_RATES(:,NUM_SED_LAYERS, J))

                BURIAL_DERIVS   (:,NUM_SED_LAYERS, j) = &
                    SED_BURRIAL_RATES (:,NUM_SED_LAYERS, j)

                ADVECTION_DERIVS(:,NUM_SED_LAYERS, j) = &
                    SED_IN_ADVEC_RATES (:,NUM_SED_LAYERS, j) - &
                    SED_OUT_ADVEC_RATES(:,NUM_SED_LAYERS, j)

                PART_MIXING_DERIVS(:,NUM_SED_LAYERS, J) = 0. !temporary fix

                TRANSPORT_DERIVS(:,NUM_SED_LAYERS, j) = &
                    DIFFUSION_DERIVS  (:,NUM_SED_LAYERS, j) * SED_POROSITIES(:,NUM_SED_LAYERS) + & !diffusion derivs multipl. by poros. to have deriv. in total volume &
                    BURIAL_DERIVS(:,NUM_SED_LAYERS, j) + &
                    ADVECTION_DERIVS(:,NUM_SED_LAYERS, j) * SED_POROSITIES(:,NUM_SED_LAYERS) + & !advection derivs multipl. by poros. to have deriv. in total volume & &
                    PART_MIXING_DERIVS(:,NUM_SED_LAYERS, j)
            end if !NUM_SED_LAYERS
        end do !j = 1, NUM_SED_VARS
        !End last layer processing for transport

        !DEPOSITION(SETTLING), EROSION  ETC.
        if ( isedi .le. 0) then    !-------------------------------------
            !24/06/2012: If clause Added by Ali for safer switch off
            if (switch_settling .ne. 0) then
                SETTLING_AFFECTED_DEPTH(:) = &
                    sum(SED_DEPTHS(:,1:NUM_FLUX_RECEIVING_SED_LAYERS),2)

                do i = 1, NUM_FLUX_RECEIVING_SED_LAYERS
                    do j = 1, NUM_SED_VARS
                        SETTLING_DERIVS(:,i, j)  = &
                            FLUXES_TO_SEDIMENTS(:,j) * &
                            (SED_DEPTHS(:,i) / SETTLING_AFFECTED_DEPTH(:))

                        PROCESSES_sed(:,i, j, 20) = SETTLING_DERIVS(:,i, j)
                    end do
                end do
            end if  !switch_settling
        else   !isedi .ne. 0     !-------------------------------------
            do j = 1, NUM_SIMULATED_SED_VARS
                !Deposition. Deposited layer concentrations are calculated
                !Used later for recalculation of BS layers conc.
                if (switch_deposition .ne. 0) then
                    if (in_which_phase(j) .eq. 1) then
                        where (H_ERODEP(:) .lt. 0.D0)
                            !for solids
                            DEPOSIT_CONC(:, j) = &
                                abs(H_ERODEP(:))*FLUXES_TO_SEDIMENTS(:,j) * (TIME_STEP / NUM_SUB_TIME_STEPS)
                        end where
                    else if(in_which_phase(j) .eq. 0 .and. in_which_phase(j) .eq. 2) then ! it is assummed that adsorbed nh4 and po4 do not settle
                        where (H_ERODEP(:) .lt. 0.D0)    ! for a while. fixme (introduce adsorbed fractions in WC)
                            DEPOSIT_CONC(:, j) = sed_porosities(:,1) * SURF_WATER_CONCS(:,j)
                        end where
                    end if !in_which_phase
                end if !switch_deposition

                !erosion
                if (switch_erosion .ne. 0) then
                    where (H_ERODEP(:) .gt. 0.D0)
                        ! solid phase still goes nowhere. fixme
                        RESUSP_FLUX(:, j) = H_ERODEP(:) * INTERMED_RESULTS(:,1, j)*SOLUTE_FRACTIONS(:,1, j)

                        FLUXES_FROM_SEDIMENTS(:,j) = FLUXES_FROM_SEDIMENTS(:,j) + RESUSP_FLUX(:, j)
                        ! for a while it goes nowhere or is added to dissolved. fixme
                    end where
                end if !switch_erosion
            end do !j = 1, NUM_SED_VARS
        end if !isedi  !-------------------------------------
        !END DEPOSITION(SETTLING), EROSION  ETC.
        !----------------------------------END TRANSPORT DERIVATIVES---------------------------------------------
         
        !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
        !KINETIC DERIVATIVES

        ! Are concentrations of state variables correctly passed to co2sys? fixme
        !24/06/2012: If clause Added by Ali for safer switch off
        if (switch_kinetics .ne. 0) then
            if (RUN_CO2SYS .eq. 1) then
                if (SINGLE_VECTOR_CO2SYS > 0) then
                    allocate(ALL_LAYERS_CO2SYS_PAR1         (nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_PAR2         (nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_PAR1TYPE     (nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_PAR2TYPE     (nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_SALT         (nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_TEMPIN       (nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_TEMPOUT      (nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_PRESIN       (nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_PRESOUT      (nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_SI           (nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_PO4          (nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_pHSCALEIN    (nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_K1K2CONSTANTS(nkn * NUM_SED_LAYERS), &
                             ALL_LAYERS_CO2SYS_KSO4CONSTANTS(nkn * NUM_SED_LAYERS))
                    
                    !call ASSIGN_DBL_VECTOR_CONTENT(ALL_LAYERS_CO2SYS_PAR1, TOT_ALK(:,i) * 1.0D6)
                    
                    do i = 1, NUM_SED_LAYERS
                        ALL_LAYERS_CO2SYS_PAR1  ((((i - 1) * nkn) + 1):(i * nkn)) = TOT_ALK  (:,i) * 1.0D6
                        ALL_LAYERS_CO2SYS_PAR2  ((((i - 1) * nkn) + 1):(i * nkn)) = INORG_C  (:,i) * 1.0D6
                        ALL_LAYERS_CO2SYS_SALT  ((((i - 1) * nkn) + 1):(i * nkn)) = SALT     (:,i)
                        ALL_LAYERS_CO2SYS_TEMPIN((((i - 1) * nkn) + 1):(i * nkn)) = SED_TEMPS(:,i)
                        ALL_LAYERS_CO2SYS_SI    ((((i - 1) * nkn) + 1):(i * nkn)) = (SED_DSI (:,i) / 28.0855D0) * 1.0D3
                        ALL_LAYERS_CO2SYS_PO4   ((((i - 1) * nkn) + 1):(i * nkn)) = (SED_PO4P(:,i) / 30.9737D0) * 1.0D3
                    end do
                    
                    ALL_LAYERS_CO2SYS_PAR1TYPE     (:) = 1
                    ALL_LAYERS_CO2SYS_PAR2TYPE     (:) = 2
                    ALL_LAYERS_CO2SYS_pHSCALEIN    (:) = 1
                    ALL_LAYERS_CO2SYS_K1K2CONSTANTS(:) = 4
                    ALL_LAYERS_CO2SYS_KSO4CONSTANTS(:) = 1
                    ALL_LAYERS_CO2SYS_TEMPOUT      (:) = 0.0D0  !Does not matter for this example
                    ALL_LAYERS_CO2SYS_PRESIN       (:) = 0.0D0  !Does not matter for this example
                    ALL_LAYERS_CO2SYS_PRESOUT      (:) = 0.0D0  !Does not matter for this example
                    
                    CO2SYS_ntps = nkn * NUM_SED_LAYERS
                    
                    call CO2SYS(ALL_LAYERS_CO2SYS_PAR1         , ALL_LAYERS_CO2SYS_PAR2         , ALL_LAYERS_CO2SYS_PAR1TYPE , &
                                ALL_LAYERS_CO2SYS_PAR2TYPE     , ALL_LAYERS_CO2SYS_SALT         , ALL_LAYERS_CO2SYS_TEMPIN   , &
                                ALL_LAYERS_CO2SYS_TEMPOUT      , ALL_LAYERS_CO2SYS_PRESIN       , ALL_LAYERS_CO2SYS_PRESOUT  , &
                                ALL_LAYERS_CO2SYS_SI           , ALL_LAYERS_CO2SYS_PO4          , ALL_LAYERS_CO2SYS_pHSCALEIN, &
                                ALL_LAYERS_CO2SYS_K1K2CONSTANTS, ALL_LAYERS_CO2SYS_KSO4CONSTANTS, ALL_LAYERS_CO2SYS_OUT_DATA , &
                                CO2SYS_NICEHEADERS             , CO2SYS_ntps)
                    
                    do i = 1, NUM_SED_LAYERS
                        !LOWER_INDEX = (((i - 1) * nkn) + 1)
                        !UPPER_INDEX =  (i * nkn)
                        pH         (:,i) = ALL_LAYERS_CO2SYS_OUT_DATA((((i - 1) * nkn) + 1):(i * nkn), 18)
                        K_ONE_TIP  (:,i) = ALL_LAYERS_CO2SYS_OUT_DATA((((i - 1) * nkn) + 1):(i * nkn), 75)
                        K_TWO_TIP  (:,i) = ALL_LAYERS_CO2SYS_OUT_DATA((((i - 1) * nkn) + 1):(i * nkn), 76)
                        K_THREE_TIP(:,i) = ALL_LAYERS_CO2SYS_OUT_DATA((((i - 1) * nkn) + 1):(i * nkn), 77)
                        HCO3       (:,i) = ALL_LAYERS_CO2SYS_OUT_DATA((((i - 1) * nkn) + 1):(i * nkn), 21)
                        CO3        (:,i) = ALL_LAYERS_CO2SYS_OUT_DATA((((i - 1) * nkn) + 1):(i * nkn), 22)
                    end do
                    
                    HCO3 = HCO3 / 1.0D6
                    CO3  =  CO3 / 1.0D6
                    
                    if(allocated(ALL_LAYERS_CO2SYS_PAR1         )) deallocate(ALL_LAYERS_CO2SYS_PAR1         )
                    if(allocated(ALL_LAYERS_CO2SYS_PAR2         )) deallocate(ALL_LAYERS_CO2SYS_PAR2         )
                    if(allocated(ALL_LAYERS_CO2SYS_PAR1TYPE     )) deallocate(ALL_LAYERS_CO2SYS_PAR1TYPE     )
                    if(allocated(ALL_LAYERS_CO2SYS_PAR2TYPE     )) deallocate(ALL_LAYERS_CO2SYS_PAR2TYPE     )
                    if(allocated(ALL_LAYERS_CO2SYS_SALT         )) deallocate(ALL_LAYERS_CO2SYS_SALT         )
                    if(allocated(ALL_LAYERS_CO2SYS_TEMPIN       )) deallocate(ALL_LAYERS_CO2SYS_TEMPIN       )
                    if(allocated(ALL_LAYERS_CO2SYS_TEMPOUT      )) deallocate(ALL_LAYERS_CO2SYS_TEMPOUT      )
                    if(allocated(ALL_LAYERS_CO2SYS_PRESIN       )) deallocate(ALL_LAYERS_CO2SYS_PRESIN       )
                    if(allocated(ALL_LAYERS_CO2SYS_PRESOUT      )) deallocate(ALL_LAYERS_CO2SYS_PRESOUT      )
                    if(allocated(ALL_LAYERS_CO2SYS_SI           )) deallocate(ALL_LAYERS_CO2SYS_SI           )
                    if(allocated(ALL_LAYERS_CO2SYS_PO4          )) deallocate(ALL_LAYERS_CO2SYS_PO4          )
                    if(allocated(ALL_LAYERS_CO2SYS_pHSCALEIN    )) deallocate(ALL_LAYERS_CO2SYS_pHSCALEIN    )
                    if(allocated(ALL_LAYERS_CO2SYS_K1K2CONSTANTS)) deallocate(ALL_LAYERS_CO2SYS_K1K2CONSTANTS)
                    if(allocated(ALL_LAYERS_CO2SYS_KSO4CONSTANTS)) deallocate(ALL_LAYERS_CO2SYS_KSO4CONSTANTS)
                    if(allocated(ALL_LAYERS_CO2SYS_OUT_DATA     )) deallocate(ALL_LAYERS_CO2SYS_OUT_DATA     )
                    if(allocated(CO2SYS_NICEHEADERS             )) deallocate(CO2SYS_NICEHEADERS             )
                else
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
                    
                    do i = 1, NUM_SED_LAYERS
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
                    
                        call CO2SYS(CO2SYS_PAR1         , CO2SYS_PAR2         , CO2SYS_PAR1TYPE , &
                                    CO2SYS_PAR2TYPE     , CO2SYS_SALT         , CO2SYS_TEMPIN   , &
                                    CO2SYS_TEMPOUT      , CO2SYS_PRESIN       , CO2SYS_PRESOUT  , &
                                    CO2SYS_SI           , CO2SYS_PO4          , CO2SYS_pHSCALEIN, &
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
                    end do !i = 1,NUM_SED_LAYERS
                    
                    !$OMP END DO NOWAIT
                    !$OMP END PARALLEL
                end if
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

            if (DO_ADVANCED_REDOX_SIMULATION > 0) then
                call SED_DOC_MINERALIZATION &
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
            else
                call CALCULATE_PH_CORR_SED &
                     (PH_CORR_DOC_MIN_DOXY   , PH , SED_PH_MIN_DOC_MIN_DOXY, &
                      SED_PH_MAX_DOC_MIN_DOXY, nkn, NUM_SED_LAYERS)
             
                call CALCULATE_PH_CORR_SED &
                    (PH_CORR_DOC_MIN_NO3N    , PH , SED_PH_MIN_DON_MIN_NO3N, &
                     SED_PH_MAX_DOC_MIN_NO3N , nkn, NUM_SED_LAYERS)

                LIM_DOXY_RED = SED_DOXY  / (SED_DOXY + SED_K_HS_DOXY_RED_LIM)

                LIM_NO3N_RED = (SED_NO3N / (SED_NO3N + SED_K_HS_NO3N_RED_LIM)) * &
                    (SED_K_HS_DOXY_RED_INHB / (SED_DOXY + SED_K_HS_DOXY_RED_INHB))                
                
                R_MINER_DOC_DOXY = &
                    SED_K_MIN_DOC_DOXY_20 * &
                    (SED_THETA_K_MIN_DOC_DOXY ** (SED_TEMPS - 2.0D1)) * &
                    LIM_DOXY_RED * PH_CORR_DOC_MIN_DOXY * &
                    (SED_DOC / (SED_DOC + SED_K_HS_DOC_MIN_DOXY)) * SED_DOC

                R_MINER_DOC_NO3N = &
                    SED_K_MIN_DOC_NO3N_20 * &
                    (SED_THETA_K_MIN_DOC_NO3N ** (SED_TEMPS - 2.0D1)) * &
                    LIM_NO3N_RED * PH_CORR_DOC_MIN_NO3N * &
                    (SED_DOC / (SED_DOC + SED_K_HS_DOC_MIN_NO3N)) * SED_DOC
            end if
            ! -----------------------------------------------------------------------------------------
            ! 28 January 2016, enhanced mineralization introduced to be comapible with redox sequence
            ! -----------------------------------------------------------------------------------------

            ! ---------------------------------------------------------------------------
            ! Dissolved organic nitrogen mineralization
            ! ---------------------------------------------------------------------------
            call CALCULATE_PH_CORR_SED &
                 (PH_CORR_DON_MIN_DOXY      , PH , SED_PH_MIN_DON_MIN_DOXY    , &
                  SED_PH_MAX_DON_MIN_DOXY   , nkn, NUM_SED_LAYERS)
             
            call CALCULATE_PH_CORR_SED &
                (PH_CORR_DON_MIN_NO3N       , PH , SED_PH_MIN_DON_MIN_NO3N    , &
                 SED_PH_MAX_DON_MIN_NO3N    , nkn, NUM_SED_LAYERS)

            R_MINER_DON_DOXY = &
                SED_K_MIN_DON_DOXY_20 * &
                (SED_THETA_K_MIN_DON_DOXY ** (SED_TEMPS - 2.0D1)) * &
                LIM_DOXY_RED * PH_CORR_DON_MIN_DOXY * &
                (SED_DON / (SED_DON + SED_K_HS_DON_MIN_DOXY)) * SED_DON

            R_MINER_DON_NO3N = &
                SED_K_MIN_DON_NO3N_20 * &
                (SED_THETA_K_MIN_DON_NO3N ** (SED_TEMPS - 2.0D1)) * &
                LIM_NO3N_RED * PH_CORR_DON_MIN_NO3N * &
                (SED_DON / (SED_DON + SED_K_HS_DON_MIN_NO3N)) * SED_DON
            
            R_MINER_DON_MN_IV    = 0.0D0
            R_MINER_DON_FE_III   = 0.0D0
            R_MINER_DON_S_PLUS_6 = 0.0D0
            R_MINER_DON_DOC      = 0.0D0
                
            if (DO_ADVANCED_REDOX_SIMULATION > 0) then

                call CALCULATE_PH_CORR_SED &
                     (PH_CORR_DON_MIN_MN_IV     , PH , SED_PH_MIN_DON_MIN_MN_IV   , &
                      SED_PH_MAX_DON_MIN_MN_IV  , nkn, NUM_SED_LAYERS)
                
                call CALCULATE_PH_CORR_SED &
                     (PH_CORR_DON_MIN_FE_III    , PH , SED_PH_MIN_DON_MIN_FE_III  , &
                     SED_PH_MAX_DON_MIN_FE_III  , nkn, NUM_SED_LAYERS)
                
                call CALCULATE_PH_CORR_SED &
                     (PH_CORR_DON_MIN_S_PLUS_6  , PH , SED_PH_MIN_DON_MIN_S_PLUS_6, &
                     SED_PH_MAX_DON_MIN_S_PLUS_6, nkn, NUM_SED_LAYERS)
                
                call CALCULATE_PH_CORR_SED &
                     (PH_CORR_DON_MIN_DOC       , PH , SED_PH_MIN_DON_MIN_DOC     , &
                      SED_PH_MAX_DON_MIN_DOC    , nkn, NUM_SED_LAYERS)
                
                R_MINER_DON_MN_IV = &
                    SED_K_MIN_DON_MN_IV_20 * &
                    (SED_THETA_K_MIN_DON_MN_IV ** (SED_TEMPS - 2.0D1)) * &
                    LIM_MN_IV_RED * PH_CORR_DON_MIN_MN_IV * &
                    (SED_DON / (SED_DON + SED_K_HS_DON_MIN_MN_IV)) * SED_DON
                
                R_MINER_DON_FE_III = &
                    SED_K_MIN_DON_FE_III_20 * &
                    (SED_THETA_K_MIN_DON_FE_III ** (SED_TEMPS - 2.0D1)) * &
                    LIM_FE_III_RED * PH_CORR_DON_MIN_FE_III * &
                    (SED_DON / (SED_DOC + SED_K_HS_DON_MIN_FE_III)) * SED_DON
                
                R_MINER_DON_S_PLUS_6 = &
                    SED_K_MIN_DON_S_PLUS_6_20 * &
                    (SED_THETA_K_MIN_DON_S_PLUS_6 ** (SED_TEMPS - 2.0D1)) * &
                    LIM_S_PLUS_6_RED * PH_CORR_DON_MIN_S_PLUS_6 * &
                    (SED_DON / (SED_DON + SED_K_HS_DON_MIN_S_PLUS_6)) * SED_DON
                
                R_MINER_DON_DOC = &
                    (SED_K_MIN_DON_DOC_20 * &
                     (SED_THETA_K_MIN_DON_DOC ** (SED_TEMPS - 2.0D1)) * &
                     LIM_DOC_RED * PH_CORR_DOC_MIN_DOXY * &
                    (SED_DON / (SED_DON + SED_K_HS_DON_MIN_DOC)) * SED_DON)                
            end if
            ! ----------------------------------------------------------------------------------------------------
            ! A bug fixed thanks to Elin's request for COCOA project. R_MINER_DON is used in different places such
            ! as non-conserving alkalinity. R_MINER_DON array is not obsolete any more. Fix date: 22 October 2017
            ! ----------------------------------------------------------------------------------------------------
            R_MINER_DON = R_MINER_DON_DOXY   + R_MINER_DON_NO3N     + R_MINER_DON_MN_IV + &
                          R_MINER_DON_FE_III + R_MINER_DON_S_PLUS_6 + R_MINER_DON_DOC
            ! ----------------------------------------------------------------------------------------------------
  
            ! ---------------------------------------------------------------------------
            ! End of dissolved organic nitrogen mineralization
            ! ---------------------------------------------------------------------------

            
            ! ---------------------------------------------------------------------------
            ! Dissolved organic phosphorus mineralization
            ! ---------------------------------------------------------------------------

                        
            ! Dissolved organic phosphorus mineralization
            ! Bug fix : 10/2/2017 pH correction must be called.
            call CALCULATE_PH_CORR_SED &
                 (PH_CORR_DOP_MIN_DOXY    , PH, SED_PH_MIN_DOP_MIN_DOXY    , &
                  SED_PH_MAX_DOP_MIN_DOXY    , nkn, NUM_SED_LAYERS)
            
            call CALCULATE_PH_CORR_SED &
                 (PH_CORR_DOP_MIN_NO3N    , PH, SED_PH_MIN_DOP_MIN_NO3N    , &
                  SED_PH_MAX_DOP_MIN_NO3N    , nkn, NUM_SED_LAYERS)

            R_MINER_DOP_DOXY = &
                SED_K_MIN_DOP_DOXY_20 * (SED_THETA_K_MIN_DOP_DOXY ** (SED_TEMPS - 2.0D1)) * &
                LIM_DOXY_RED * PH_CORR_DOP_MIN_DOXY * (SED_DOP / (SED_DOP + SED_K_HS_DOP_MIN_DOXY)) * &
                SED_DOP

            R_MINER_DOP_NO3N = &
                SED_K_MIN_DOP_NO3N_20 * (SED_THETA_K_MIN_DOP_NO3N ** (SED_TEMPS - 2.0D1)) * &
                LIM_NO3N_RED * PH_CORR_DOP_MIN_NO3N * (SED_DOP / (SED_DOP + SED_K_HS_DOP_MIN_NO3N)) * &
                SED_DOP                  
                  

            R_MINER_DOP_MN_IV    = 0.0D0
            R_MINER_DOP_FE_III   = 0.0D0
            R_MINER_DOP_S_PLUS_6 = 0.0D0
            R_MINER_DOP_DOC      = 0.0D0
                
            if (DO_ADVANCED_REDOX_SIMULATION > 0) then
                                  
                call CALCULATE_PH_CORR_SED &
                     (PH_CORR_DOP_MIN_MN_IV   , PH, SED_PH_MIN_DOP_MIN_MN_IV   , &
                      SED_PH_MAX_DOP_MIN_MN_IV   , nkn, NUM_SED_LAYERS)
                
                call CALCULATE_PH_CORR_SED &
                     (PH_CORR_DOP_MIN_FE_III  , PH, SED_PH_MIN_DOP_MIN_FE_III  , &
                      SED_PH_MAX_DOP_MIN_FE_III  , nkn, NUM_SED_LAYERS)
                
                call CALCULATE_PH_CORR_SED &
                     (PH_CORR_DOP_MIN_S_PLUS_6, PH, SED_PH_MIN_DOP_MIN_S_PLUS_6, &
                      SED_PH_MAX_DOP_MIN_S_PLUS_6, nkn, NUM_SED_LAYERS)
                
                call CALCULATE_PH_CORR_SED &
                     (PH_CORR_DOP_MIN_DOC     , PH, SED_PH_MIN_DOP_MIN_DOC     , &
                      SED_PH_MAX_DOP_MIN_DOC     , nkn, NUM_SED_LAYERS)
                
                R_MINER_DOP_MN_IV = &
                    SED_K_MIN_DOP_MN_IV_20 * (SED_THETA_K_MIN_DOP_MN_IV ** (SED_TEMPS - 2.0D1)) * &
                    LIM_MN_IV_RED * PH_CORR_DOP_MIN_MN_IV * (SED_DOP / (SED_DOP + SED_K_HS_DOP_MIN_MN_IV)) * &
                    SED_DOP
                
                R_MINER_DOP_FE_III = &
                    SED_K_MIN_DOP_FE_III_20 * (SED_THETA_K_MIN_DOP_FE_III ** (SED_TEMPS - 2.0D1)) * &
                    LIM_FE_III_RED * PH_CORR_DOP_MIN_FE_III * (SED_DOP / (SED_DOC + SED_K_HS_DOP_MIN_FE_III)) * &
                    SED_DOP
                
                R_MINER_DOP_S_PLUS_6 = &
                    SED_K_MIN_DOP_S_PLUS_6_20 * (SED_THETA_K_MIN_DOP_S_PLUS_6 ** (SED_TEMPS - 2.0D1)) * &
                    LIM_S_PLUS_6_RED * PH_CORR_DOP_MIN_S_PLUS_6 * (SED_DOP / (SED_DOP + SED_K_HS_DOP_MIN_S_PLUS_6)) * &
                    SED_DOP
                
                R_MINER_DOP_DOC = &
                    (SED_K_MIN_DOP_DOC_20 * (SED_THETA_K_MIN_DOP_DOC ** (SED_TEMPS - 2.0D1)) * &
                     LIM_DOC_RED * PH_CORR_DOC_MIN_DOXY * (SED_DOP / (SED_DOP + SED_K_HS_DOP_MIN_DOC)) * SED_DOP)
                
                ! ----------------------------------------------------------------------------------------------------
                ! A bug fixed thanks to Elin's request for COCOA project. R_MINER_DOP is used in different places such
                ! as non-conserving alkalinity. R_MINER_DOP array is not obsolete any more. Fix date: 22 October 2017
                ! ----------------------------------------------------------------------------------------------------
                R_MINER_DOP = R_MINER_DOP_DOXY   + R_MINER_DOP_NO3N     + R_MINER_DOP_MN_IV + &
                              R_MINER_DOP_FE_III + R_MINER_DOP_S_PLUS_6 + R_MINER_DOP_DOC
                ! ----------------------------------------------------------------------------------------------------
            end if
             
            ! ---------------------------------------------------------------------------
            ! End of dissolved organic phosphorus mineralization
            ! ---------------------------------------------------------------------------

    
            ! ----------------------------------------------------------------------------------------------------
            ! End of 28 January 2016, enhanced mineralization introduced to be comapible with redox sequence
            ! ----------------------------------------------------------------------------------------------------

            where (SED_DOXY .GE. DOXY_AT_ANOXIA)
                R_DISS_PON  = K_OXIC_DISS_PON    * (THETA_DISS_PON  ** (SED_TEMPS - 2.0D1)) * &
                    (SED_PON / (SED_PON + KHS_DISS_PON))  * SED_PON

                R_DISS_POP  = K_OXIC_DISS_POP    * (THETA_DISS_POP  ** (SED_TEMPS - 2.0D1)) * &
                    (SED_POP / (SED_POP + KHS_DISS_POP))  * SED_POP

                R_DISS_PSi  = K_OXIC_DISS_PSi    * (THETA_DISS_PSi  ** (SED_TEMPS - 2.0D1)) * &
                    (SED_PSi / (SED_PSi + KHS_DISS_PSi))  * SED_PSi
            end where

            where (SED_DOXY .LT. DOXY_AT_ANOXIA)
                R_DISS_PON  = K_ANOXIC_DISS_PON  * (THETA_DISS_PON  ** (SED_TEMPS - 2.0D1)) * &
                    (SED_PON / (SED_PON + KHS_DISS_PON))  * SED_PON

                R_DISS_POP  = K_ANOXIC_DISS_POP  * (THETA_DISS_POP  ** (SED_TEMPS - 2.0D1)) * &
                    (SED_POP / (SED_POP + KHS_DISS_POP))  * SED_POP

                R_DISS_PSi  = K_ANOXIC_DISS_PSi  * (THETA_DISS_PSi  ** (SED_TEMPS - 2.0D1)) * &
                    (SED_PSi / (SED_PSi + KHS_DISS_PSi))  * SED_PSi
            end where

            ! Nitrification
            call CALCULATE_PH_CORR_SED(PH_CORR_NITR_NH4, PH, SED_PH_NITR_NH4_MIN, &
                SED_PH_NITR_NH4_MAX, nkn, NUM_SED_LAYERS)

            R_NITR   = K_NITR * (THETA_NITR ** (SED_TEMPS - 2.0D1)) * &
                       (SED_DOXY / (SED_DOXY + KHS_NITR_DOXY)) * &
                       (SED_NH4N / (SED_NH4N + KHS_NITR_NH4N)) * &
                       PH_CORR_NITR_NH4 * SED_NH4N

            ! 29 January 2016
            ! Following commented lines are replaced by the new redox sequence based DOC mineralization
            ! Denitrification
            ! call CALCULATE_PH_CORR_SED(PH_CORR_DENITR_NO3, PH, SED_PH_DENITR_NO3_MIN, &
            !                       SED_PH_DENITR_NO3_MAX, nkn, NUM_SED_LAYERS)

            ! R_DENITR = K_DENITR * (THETA_DENITR ** (SED_TEMPS - 2.0D1)) * &
            !            (KHS_DENITR_DOXY / (SED_DOXY + KHS_DENITR_DOXY))* &
            !            (SED_NO3N / (SED_NO3N + KHS_DENITR_NO3N)) * &
            !            (SED_DOC  / (SED_DOC  + KHS_DENITR_DOC))  * &
            !            PH_CORR_DENITR_NO3 * SED_NO3N

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

            R_DENITR = 0.93D0 * R_MINER_DOC_NO3N
            ! -------------------------------------------------------------------------
            ! END OF DENITRIFICATION
            ! -------------------------------------------------------------------------


            if (DO_ADVANCED_REDOX_SIMULATION > 0) then
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
                 
                 R_MN_IV_REDUCTION = 8.66D0 * R_MINER_DOC_MN_IV
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
                 
                 R_FE_III_REDUCTION = 18.66D0 * R_MINER_DOC_FE_III
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
                 
                 R_SULPHATE_REDUCTION = 1.33D0 * R_MINER_DOC_S_PLUS_6
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
                 
                 R_METHANOGENESIS = 0.5D0 * R_MINER_DOC_DOC
                 ! -------------------------------------------------------------------------
                 ! END OF METHANOGENESIS
                 ! -------------------------------------------------------------------------
                 
                 
                 ! Iron and manganese processes introduced 10 the of September 2015
                 ! For now, no temparature corrections. Effect on temperature and other
                 ! environmental conditions may be included after more detailed investigations
                 
                 ! Iron
                 
                 ! After the introduction of changes in November 30 th 2015, it will be necessary
                 ! to reconsinder kinetics, especially the mineral disoolotion processes for metals
                 ! Job to be conducted wiht Petras together, and completed during next visit of
                 ! of Ali to Lithuania. Also for iron and mangenese, the redox sequences as described
                 ! by Van Chappen and Wang 2015 and Katsev papers should be included.
                 
                 where (SED_DOXY < 1)
                     R_FE_II_OXIDATION  = &
                         sed_k_OX_FE_II * SED_DOXY * (10.0D0 ** (PH - 7.0D0)) * FE_II * MULT_FE_II_DISS
                 elsewhere
                     R_FE_II_OXIDATION  = &
                         sed_k_OX_FE_II * (10.0D0 ** (PH - 7.0D0)) * FE_II * MULT_FE_II_DISS
                 end where
                 
                 ! 29 January 2016
                 ! Following commented lines are replaced by the new redox sequence based DOC mineralization
                 ! it is the next visit of Ali and the redox sequences as described
                 ! by Van Chappen and Wang 2015 and Katsev papers are now included.
                 
                 ! R_FE_III_REDUCTION = sed_k_RED_FE_III * &
                 !     (SED_KHS_DOXY_FE_III_RED / (SED_KHS_DOXY_FE_III_RED + SED_DOXY)) * &
                 !     FE_III * MULT_FE_III_DISS
                 
                 ! Manganese
                 where (SED_DOXY < 1)
                     R_MN_II_OXIDATION  = sed_k_OX_MN_II * SED_DOXY * (10.0D0 ** (PH - 7.0D0)) * &
                         MN_II * MULT_MN_II_DISS
                 elsewhere
                     R_MN_II_OXIDATION  = sed_k_OX_MN_II * (10.0D0 ** (PH - 7.0D0)) * &
                         MN_II * MULT_MN_II_DISS
                 end where
                 
                 ! 29 January 2016
                 ! Following commented lines are replaced by the new redox sequence based DOC mineralization
                 ! it is the next visit of Ali and the redox sequences as described
                 ! by Van Chappen and Wang 2015 and Katsev papers are now included.
                 
                 ! R_MN_IV_REDUCTION  = sed_k_RED_MN_IV  * &
                 !     (SED_KHS_DOXY_MN_IV_RED / (SED_KHS_DOXY_MN_IV_RED + SED_DOXY)) * MN_IV * MULT_MN_IV_DISS
                 
                 ! End of ron and manganese processes introduced 10 the of September 2015
                 
                 ! -------------------------------------------------------------------------
                 ! 29 January 2016 KINETICS OF NEW STATE VARIABLES
                 ! -------------------------------------------------------------------------
                 R_METHANE_OXIDATION = &
                     SED_k_OX_CH4 * (SED_THETA_k_OX_CH4 ** (SED_TEMPS - 20.0D0)) * CH4_C * &
                     (SED_DOXY / (SED_k_HS_OX_CH4_DOXY + SED_DOXY))
                 
                 R_SULPHIDE_OXIDATION = &
                     SED_k_OX_H2S * (SED_THETA_k_OX_H2S ** (SED_TEMPS - 20.0D0)) * S_MINUS_2 * &
                     (SED_DOXY / (SED_k_HS_OX_H2S_DOXY + SED_DOXY))
                 ! -------------------------------------------------------------------------
                 ! END OF KINETICS OF NEW STATE VARIABLES
                 ! -------------------------------------------------------------------------
            end if

            DEOXYGENATION = &
                (2.66D0 * R_MINER_DOC_DOXY)  + (0.43D0 * R_FE_II_OXIDATION)   + &
                (0.88D0 * R_MN_II_OXIDATION) + (2.66D0 * R_METHANE_OXIDATION) + &
                (2.00D0 * R_SULPHIDE_OXIDATION)

            R_MINER_DOC = R_MINER_DOC_DOXY   + R_MINER_DOC_NO3N     + R_MINER_DOC_MN_IV + & 
                          R_MINER_DOC_FE_III + R_MINER_DOC_S_PLUS_6 + R_MINER_DOC_DOC
            
            ! Converting to unit area
            R_DISS_POC           = R_DISS_POC           * SED_DEPTHS !* (1.0D0 - SED_POROSITIES)    ! fixed: state vars are already in cocentrations per
            R_DISS_PON           = R_DISS_PON           * SED_DEPTHS !* (1.0D0 - SED_POROSITIES)    ! total unit volume of sediments
            R_DISS_POP           = R_DISS_POP           * SED_DEPTHS !* (1.0D0 - SED_POROSITIES)
            R_DISS_PSi           = R_DISS_PSi           * SED_DEPTHS !* (1.0D0 - SED_POROSITIES)
            R_MINER_DOC          = R_MINER_DOC          * SED_DEPTHS !* SED_POROSITIES
            R_MINER_DON          = R_MINER_DON          * SED_DEPTHS !* SED_POROSITIES
            R_MINER_DOP          = R_MINER_DOP          * SED_DEPTHS !* SED_POROSITIES
            DEOXYGENATION        = DEOXYGENATION        * SED_DEPTHS !* SED_POROSITIES
            R_NITR               = R_NITR               * SED_DEPTHS !* SED_POROSITIES
            R_DENITR             = R_DENITR             * SED_DEPTHS !* SED_POROSITIES
            R_FE_II_OXIDATION    = R_FE_II_OXIDATION    * SED_DEPTHS !* SED_POROSITIES
            R_FE_III_REDUCTION   = R_FE_III_REDUCTION   * SED_DEPTHS !* SED_POROSITIES
            R_MN_II_OXIDATION    = R_MN_II_OXIDATION    * SED_DEPTHS !* SED_POROSITIES
            R_MN_IV_REDUCTION    = R_MN_IV_REDUCTION    * SED_DEPTHS !* SED_POROSITIES
            R_MINER_DOC_DOXY     = R_MINER_DOC_DOXY     * SED_DEPTHS
            R_MINER_DOC_NO3N     = R_MINER_DOC_NO3N     * SED_DEPTHS
            R_MINER_DOC_MN_IV    = R_MINER_DOC_MN_IV    * SED_DEPTHS
            R_MINER_DOC_FE_III   = R_MINER_DOC_FE_III   * SED_DEPTHS
            R_MINER_DOC_S_PLUS_6 = R_MINER_DOC_S_PLUS_6 * SED_DEPTHS
            R_MINER_DOC_DOC      = R_MINER_DOC_DOC      * SED_DEPTHS
            R_MINER_DON_DOXY     = R_MINER_DON_DOXY     * SED_DEPTHS
            R_MINER_DON_NO3N     = R_MINER_DON_NO3N     * SED_DEPTHS
            R_MINER_DON_MN_IV    = R_MINER_DON_MN_IV    * SED_DEPTHS
            R_MINER_DON_FE_III   = R_MINER_DON_FE_III   * SED_DEPTHS
            R_MINER_DON_S_PLUS_6 = R_MINER_DON_S_PLUS_6 * SED_DEPTHS
            R_MINER_DON_DOC      = R_MINER_DON_DOC      * SED_DEPTHS
            R_MINER_DOP_DOXY     = R_MINER_DOP_DOXY     * SED_DEPTHS
            R_MINER_DOP_NO3N     = R_MINER_DOP_NO3N     * SED_DEPTHS
            R_MINER_DOP_MN_IV    = R_MINER_DOP_MN_IV    * SED_DEPTHS
            R_MINER_DOP_FE_III   = R_MINER_DOP_FE_III   * SED_DEPTHS
            R_MINER_DOP_S_PLUS_6 = R_MINER_DOP_S_PLUS_6 * SED_DEPTHS
            R_MINER_DOP_DOC      = R_MINER_DOP_DOC      * SED_DEPTHS
            R_SULPHIDE_OXIDATION = R_SULPHIDE_OXIDATION * SED_DEPTHS
            R_SULPHATE_REDUCTION = R_SULPHATE_REDUCTION * SED_DEPTHS
            R_METHANOGENESIS     = R_METHANOGENESIS     * SED_DEPTHS
            R_METHANE_OXIDATION  = R_METHANE_OXIDATION  * SED_DEPTHS

            KINETIC_DERIVS(:,:, 1) = &
                R_MINER_DON_DOXY  (:,:) + R_MINER_DON_NO3N    (:,:) + R_MINER_DON_MN_IV(:,:) + &
                R_MINER_DON_FE_III(:,:) + R_MINER_DON_S_PLUS_6(:,:) + R_MINER_DON_DOC  (:,:) - &
                R_NITR(:,:)

            if(debug_stranger) then
                do i = 1, NUM_SED_LAYERS
                    if(STRANGERSD(KINETIC_DERIVS(:,i, 1),VALUE_strange,nkn) .eq. 1) then
                        nstrange = count(VALUE_strange)
                        allocate(STRANGERS    (nstrange))
                        allocate(NODES_STRANGE(nstrange))

                        l=1

                        do k=1,nkn
                            if(VALUE_strange(k)) then
                                STRANGERS    (l) = KINETIC_DERIVS(k,i,1)
                                NODES_STRANGE(l) = k
                                l=l+1
                            end if
                        end do

                        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                        print *, 'aquabc_sediment1: Layer ', i, 'Variable 1'
                        print *, 'Kinetic derivative is NaN'
                        print *, 'NODE_NUMBERS=',NODES_STRANGE
                        print *, 'VALUES=',STRANGERS

                        print *, ' '
                        print *, 'Other variables:'
                        
                        print *, ' PH correction for nitrification:', &
                                 (PH_CORR_NITR_NH4(NODES_STRANGE(k), i), k=1,nstrange)
                        
                        print *, ' DON mineralisation rate using oxygen  :', &
                                 (R_MINER_DON_DOXY(NODES_STRANGE(k), i),      k=1,nstrange)
                        
                        print *, '    SED_K_MIN_DON_DOXY_20    : ', SED_K_MIN_DON_DOXY_20
                        print *, '    SED_THETA_K_MIN_DON_DOXY : ', SED_THETA_K_MIN_DON_DOXY
                        print *, '    SED_TEMPS                : ', SED_TEMPS
                        print *, '    LIM_DOXY_RED             : ', LIM_DOXY_RED
                        print *, '    PH_CORR_DON_MIN_DOXY     : ', PH_CORR_DON_MIN_DOXY
                        print *, '    SED_DON                  : ', SED_DON
                        print *, '    SED_K_HS_DON_MIN_DOXY    : ', SED_K_HS_DON_MIN_DOXY

                        print *, ' DON mineralisation rate using nitrate :', &
                                 (R_MINER_DON_NO3N(NODES_STRANGE(k), i),      k=1,nstrange)

                        print *, ' DON mineralisation rate using Mn 4+   :', &
                                 (R_MINER_DON_MN_IV(NODES_STRANGE(k), i),      k=1,nstrange)

                        print *, ' DON mineralisation rate using Fe 3+   :', &
                                 (R_MINER_DON_FE_III(NODES_STRANGE(k), i),      k=1,nstrange)

                        print *, ' DON mineralisation rate using S 6+    :', &
                                 (R_MINER_DON_S_PLUS_6(NODES_STRANGE(k), i),      k=1,nstrange)

                        print *, ' DON mineralisation rate using DOC     :', &
                                 (R_MINER_DON_DOC(NODES_STRANGE(k), i),      k=1,nstrange)

                        print *, ' Nitrification rate:', &
                                 (R_NITR(NODES_STRANGE(k), i),      k=1,nstrange)

                        print *, ' Sediment depths rate                  :', &
                                 (SED_DEPTHS(NODES_STRANGE(k), i),      k=1,nstrange)

                        deallocate(STRANGERS, NODES_STRANGE)
                        error =1
                    end if
               end do     ! i-layers

               if (error .eq. 1) stop
            end if      !debug_stranger

            KINETIC_DERIVS(:,:, 2)  = R_NITR(:,:) - R_DENITR(:,:)

            ! if(debug_stranger) then
            !     do i = 1, NUM_SED_LAYERS
            !
            !         if(STRANGERSD(KINETIC_DERIVS(:,i, 2),VALUE_strange,nkn) .eq. 1) then
            !            nstrange = count(VALUE_strange)
            !            allocate(STRANGERS    (nstrange))
            !            allocate(NODES_STRANGE(nstrange))
            !
            !            l=1
            !            do k=1,nkn
            !             if(VALUE_strange(k)) then
            !               STRANGERS    (l) = KINETIC_DERIVS(k,i,2)
            !               NODES_STRANGE(l) = k
            !               l=l+1
            !             end if
            !            end do
            !
            !            print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            !            print *, 'aquabc_sediment1: Layer ', i, 'Variable 2'
            !            print *, 'Kinetic derivative is NaN'
            !            print *, 'NODE_NUMBERS=',NODES_STRANGE
            !            print *, 'VALUES=',STRANGERS
            !
            !            print *, ' '
            !            print *, 'Other variables:'
            !            print *, ' PH correction for nitrification:', &
            !                     (PH_CORR_DENITR_NO3(NODES_STRANGE(k), i), k=1,nstrange)
            !            print *, ' Denitrification rate:', &
            !                     (R_DENITR(NODES_STRANGE(k), i),         k=1,nstrange)
            !            print *, ' Nitrification rate:', &
            !                     (R_NITR(NODES_STRANGE(k), i),           k=1,nstrange)
            !            print *, ' Layer thickness:', &
            !                     (SED_DEPTHS(NODES_STRANGE(k),i),        k=1,nstrange)
            !            print *, ' BS temperature:', &
            !                     (SED_TEMPS(NODES_STRANGE(k),i),         k=1,nstrange)
            !            print *, ' BS oxygen:', &
            !                     (SED_DOXY(NODES_STRANGE(k),i),          k=1,nstrange)
            !            print *, ' BS NO3:', &
            !                     (SED_NO3N(NODES_STRANGE(k),i),          k=1,nstrange)
            !            print *, ' BS DOC:', &
            !                     (SED_DOC(NODES_STRANGE(k),i),           k=1,nstrange)
            !
            !           deallocate(STRANGERS, NODES_STRANGE)
            !           error =1
            !         end if
            !
            !     end do     ! i-layers
            !  if (error .eq. 1) stop
            ! end if      !debug_stranger

            KINETIC_DERIVS(:,:, 3) = R_DISS_PON(:,:) - &
                R_MINER_DON_DOXY  (:,:) - R_MINER_DON_NO3N    (:,:) - R_MINER_DON_MN_IV(:,:) - &
                R_MINER_DON_FE_III(:,:) - R_MINER_DON_S_PLUS_6(:,:) - R_MINER_DON_DOC  (:,:)

            KINETIC_DERIVS(:,:, 4)  = (-1.0D0) * R_DISS_PON(:,:)

            KINETIC_DERIVS(:,:, 5)  = &
                R_MINER_DOP_DOXY  (:,:) + R_MINER_DOP_NO3N    (:,:) + R_MINER_DOP_MN_IV(:,:) - &
                R_MINER_DOP_FE_III(:,:) + R_MINER_DOP_S_PLUS_6(:,:) + R_MINER_DOP_DOC  (:,:)

            KINETIC_DERIVS(:,:, 6)  = R_DISS_POP(:,:) - &
                R_MINER_DOP_DOXY  (:,:) - R_MINER_DOP_NO3N    (:,:) - R_MINER_DOP_MN_IV(:,:) - &
                R_MINER_DOP_FE_III(:,:) - R_MINER_DOP_S_PLUS_6(:,:) - R_MINER_DOP_DOC  (:,:)

            KINETIC_DERIVS(:,:, 7)  = (-1.0D0)   * R_DISS_POP(:,:)
            KINETIC_DERIVS(:,:, 8)  = ((-4.57D0) * R_NITR(:,:)) - DEOXYGENATION(:,:)

            KINETIC_DERIVS(:,:, 9)  = R_DISS_POC(:,:) - &
                R_MINER_DOC_DOXY  (:,:) - R_MINER_DOC_NO3N    (:,:) - R_MINER_DOC_MN_IV(:,:) - &
                R_MINER_DOC_FE_III(:,:) - R_MINER_DOC_S_PLUS_6(:,:) - R_MINER_DOC_DOC  (:,:)

            KINETIC_DERIVS(:,:, 10) = (-1.0D0) * R_DISS_POC(:,:)
            KINETIC_DERIVS(:,:, 11) = R_DISS_PSi(:,:)
            KINETIC_DERIVS(:,:, 12) = (-1.0D0) * R_DISS_PSi(:,:)

            ! 29 January 2016, the following commented lines are replaced after the
            ! introduction of redox sequence
            ! Kinetic sub model for dissolved inorganic carbon
            !TOTAL_DIC_KINETIC_SOURCES(:,:) = R_MINER_DOC(:,:)

            TOTAL_DIC_KINETIC_SOURCES(:,:) = &
                R_MINER_DOC_DOXY    (:,:) + R_MINER_DOC_NO3N  (:,:) + &
                R_MINER_DOC_MN_IV   (:,:) + R_MINER_DOC_FE_III(:,:) + &
                R_MINER_DOC_S_PLUS_6(:,:) + (0.5D0*R_MINER_DOC_DOC(:,:))


            TOTAL_DIC_KINETIC_SINKS  (:,:) = 0.0D0

            ! Calculate the total derivative and convert it to moles
            if (CONSIDER_INORG_C_DERIVATIVE > 0) then
                PROCESSES_sed(:,:, 13, 1) = TOTAL_DIC_KINETIC_SOURCES /12000.0D0
                PROCESSES_sed(:,:, 13, 2) = TOTAL_DIC_KINETIC_SINKS   /12000.0D0

                DIC_KINETIC_DERIVATIVE = &
                    ((TOTAL_DIC_KINETIC_SOURCES - TOTAL_DIC_KINETIC_SINKS) / 12000.0D0)

                KINETIC_DERIVS(:,:, 13) = DIC_KINETIC_DERIVATIVE ! should be converted to unit area mass? fixme
            else
                PROCESSES_sed (:,:, 13, 1) = 0.0D0
                PROCESSES_sed (:,:, 13, 2) = 0.0D0
                KINETIC_DERIVS(:,:, 13)    = 0
            end if

            if (DO_ADVANCED_REDOX_SIMULATION > 0) then
                ! -------------------------------------------------------------------------
                ! 29 JANUARY 2016, KINETIC DERIVATIVES FOR THE NEW STATE VARIABLES
                ! -------------------------------------------------------------------------
                
                ! Calcium
                KINETIC_DERIVS(:,:,20) = 0.0D0
                
                ! Magnesium
                KINETIC_DERIVS(:,:,21) = 0.0D0
                
                ! Suphate sulphur
                PROCESSES_sed (:,:, 22, 1) = R_SULPHIDE_OXIDATION
                PROCESSES_sed (:,:, 22, 2) = R_SULPHATE_REDUCTION
                KINETIC_DERIVS(:,:,22)     = PROCESSES_sed(:,:, 22, 1) - PROCESSES_sed(:,:, 22, 2)
                
                ! Sulphide sulphur
                PROCESSES_sed (:,:, 23, 1) = R_SULPHATE_REDUCTION
                PROCESSES_sed (:,:, 23, 2) = R_SULPHIDE_OXIDATION
                KINETIC_DERIVS(:,:, 23)    = PROCESSES_sed(:,:, 23, 1) - PROCESSES_sed(:,:, 23, 2)
                
                ! Methane carbon
                PROCESSES_sed(:,:, 24, 1) = R_METHANOGENESIS
                PROCESSES_sed(:,:, 24, 2) = R_METHANE_OXIDATION
                
                KINETIC_DERIVS(:,:,24) = PROCESSES_sed(:,:, 24, 1) - PROCESSES_sed(:,:, 24, 2)
            end if
            ! -------------------------------------------------------------------------
            ! END OF 29 JANUARY 2016, LINETIC DERIVATIVES FOR THE NEW STATE VARIABLES
            ! -------------------------------------------------------------------------



            ! After the introduction of changes in November 30 th 2015, it will be necessary
            ! to reconsider alkalinity since metals will change ion balance
            ! Job to be conducted wiht Petras together, and completed during next visit of
            ! of Ali to Lithuania. Also for iron and mangenese, the redox sequences as described
            ! by Van Chappen and Wang 2015 and Katsev papers should be included.



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
            N_DIA_TOT_RESP            = 0.0D0
            N_CYN_TOT_RESP            = 0.0D0
            N_OPA_TOT_RESP            = 0.0D0
            N_FIX_CYN_TOT_RESP        = 0.0D0
            N_ZOO_TOT_RESP            = 0.0D0
            N_ABIOTIC_DON_MIN         = R_MINER_DON * FRAC_NH4

            ALK_GAINED_BY_AMMONIUM_GEN     = &
                (N_DIA_TOT_RESP + N_CYN_TOT_RESP + N_OPA_TOT_RESP + N_FIX_CYN_TOT_RESP + &
                 N_ZOO_TOT_RESP + N_ABIOTIC_DON_MIN) / 14007.0D0
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

            ALK_GAINED_BY_NITRATE_CONS = &
               (N_DENITRIFICATION     + N_DIA_GROWTH          + &
                N_CYN_GROWTH          + N_OPA_GROWTH          + &
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
            N_DIA_GROWTH          = 0.0D0
            N_CYN_GROWTH          = 0.0D0
            N_OPA_GROWTH          = 0.0D0
            N_NON_FIX_CYN_GROWTH  = 0.0D0

            ALK_LOST_BY_AMMONIUM_CONS  = &
                (N_DIA_GROWTH + N_CYN_GROWTH + N_OPA_GROWTH + N_NON_FIX_CYN_GROWTH) / 14007.0D0
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
                     !K_ONE_TIP   = K_ONE_TIP
                     !K_TWO_TIP   = K_TWO_TIP
                     !K_THREE_TIP = K_THREE_TIP
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
            P_DIA_TOT_RESP            = 0.0D0
            P_CYN_TOT_RESP            = 0.0D0
            P_OPA_TOT_RESP            = 0.0D0
            P_FIX_CYN_TOT_RESP        = 0.0D0
            P_ZOO_TOT_RESP            = 0.0D0
            P_ABIOTIC_DON_MIN         = R_MINER_DOP * PHOSPHATE_EQ_CONSTANT

            ALK_LOST_BY_PHOSPHATE_GEN = &
                (P_DIA_TOT_RESP         + &
                 P_CYN_TOT_RESP            + P_OPA_TOT_RESP         + &
                 P_FIX_CYN_TOT_RESP        + P_ZOO_TOT_RESP         + &
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

            KINETIC_DERIVS(:,:, 14)   = ALK_KINETIC_DERIVATIVE ! shoul it to be converted to unit area mass? fixme

            ! --------------------------------------------------------
            ! NO KINETIC SUBMODEL FOR SALINITY SINCE IT IS CONSERVATIVE
            ! --------------------------------------------------------

            if (DO_ADVANCED_REDOX_SIMULATION > 0) then
                ! ----------------------------------------------------------------------------------
                ! Kinetic submodel for Iron and Manganese introduced 10 th of September 2015
                ! ----------------------------------------------------------------------------------
                KINETIC_DERIVS(:, :, 16) = R_FE_III_REDUCTION - R_FE_II_OXIDATION       ! FE_II
                KINETIC_DERIVS(:, :, 17) = R_FE_II_OXIDATION  - R_FE_III_REDUCTION      ! FE_III
                KINETIC_DERIVS(:, :, 18) = R_MN_IV_REDUCTION  - R_MN_II_OXIDATION       ! MN_II
                KINETIC_DERIVS(:, :, 19) = R_MN_II_OXIDATION  - R_MN_IV_REDUCTION       ! MN_IV
                ! ----------------------------------------------------------------------------------
                ! End of kinetic submodel for Iron and Manganese introduced 10 th of September 2015
                ! ----------------------------------------------------------------------------------
            end if
            !end KINETIC DERIVATIVES


            !KINETIC_DERIVS(SED_NH4N) = R_MINER_DON - R_NITR
            PROCESSES_sed(:,:, 1, 10) = R_MINER_DON
            PROCESSES_sed(:,:, 1, 11) = R_NITR

            !KINETIC_DERIVS(SED_NO3N) = R_NITR - R_DENITR
            PROCESSES_sed(:,:, 2, 10) = R_NITR
            PROCESSES_sed(:,:, 2, 11) = R_DENITR

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

            PROCESSES_sed(:,:, 14, 1) = pH(:,:)  !pH
            PROCESSES_sed(:,:, 14, 2) = ALK_GAINED_BY_AMMONIUM_GEN
            PROCESSES_sed(:,:, 14, 3) = ALK_GAINED_BY_NITRATE_CONS
            PROCESSES_sed(:,:, 14, 4) = ALK_GAINED_BY_PHOSPHATE_CONS
            PROCESSES_sed(:,:, 14, 5) = ALK_LOST_BY_AMMONIUM_CONS
            PROCESSES_sed(:,:, 14, 6) = ALK_LOST_BY_NITRIFICATION
            PROCESSES_sed(:,:, 14, 7) = ALK_LOST_BY_PHOSPHATE_GEN

            ! --------------------------------------------------------
            ! NO PROCESS RATES FOR SALINITY SINCE IT IS CONSERVATIVE
            ! --------------------------------------------------------

            if (DO_ADVANCED_REDOX_SIMULATION > 0) then
                PROCESSES_sed(:, :, 16, 1) = R_FE_II_OXIDATION   !Introduced 10 th of September 2015
                PROCESSES_sed(:, :, 16, 2) = R_FE_III_REDUCTION  !Introduced 10 th of September 2015
                PROCESSES_sed(:, :, 17, 1) = R_FE_II_OXIDATION   !Introduced 10 th of September 2015
                PROCESSES_sed(:, :, 17, 2) = R_FE_III_REDUCTION  !Introduced 10 th of September 2015
                PROCESSES_sed(:, :, 18, 1) = R_MN_II_OXIDATION   !Introduced 10 th of September 2015
                PROCESSES_sed(:, :, 18, 2) = R_MN_IV_REDUCTION   !Introduced 10 th of September 2015
                PROCESSES_sed(:, :, 19, 1) = R_MN_II_OXIDATION   !Introduced 10 th of September 2015
                PROCESSES_sed(:, :, 19, 2) = R_MN_IV_REDUCTION   !Introduced 10 th of September 2015
            end if
            
        end if ! switch_kinetics
        ! END KINETIC DERIVATIVES
        
        !CHECK FOR STRANGERS BEFORE TIME INTEGRATION
        error = 0

        ! CHECK TRANSPORT DERIVATIVES
        if(debug_stranger) then
            do i = 1, NUM_SED_LAYERS
                do j = 1, NUM_SIMULATED_SED_VARS
                    if(STRANGERSD(TRANSPORT_DERIVS(:,i, j),VALUE_strange,nkn) .eq. 1) then
                        nstrange = count(VALUE_strange)
                        allocate(STRANGERS    (nstrange))
                        allocate(NODES_STRANGE(nstrange))
                        l=1
						
                        do k=1,nkn
                            if(VALUE_strange(k)) then
                                STRANGERS    (l) = TRANSPORT_DERIVS(k,i,j)
                                NODES_STRANGE(l) = k
                                l=l+1
                            end if
                        end do
        
                        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
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
                    end if !STRANGERSD
                end do ! j-variables
            end do     ! i-layers
            
			if (error .eq. 1) stop
        end if      ! debug_stranger


        ! CHECK KINETIC DERIVATIVES
        if(debug_stranger) then
            do i = 1, NUM_SED_LAYERS
                do j = 1, NUM_SIMULATED_SED_VARS
                    if(STRANGERSD(KINETIC_DERIVS(:,i, j),VALUE_strange,nkn) .eq. 1) then
                        nstrange = count(VALUE_strange)
                        allocate(STRANGERS    (nstrange))
                        allocate(NODES_STRANGE(nstrange))
        
                        l=1
						
                        do k=1,nkn
                            if(VALUE_strange(k)) then
                                STRANGERS    (l) = KINETIC_DERIVS(k,i,j)
                                NODES_STRANGE(l) = k
                                l=l+1
                             end if
                        end do
        
                        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                        print *, 'aquabc_sediment1: Layer ', i 
                        print *, 'Variable ',j
                        print *, 'Kinetic derivative is NaN'
                        print *, 'NODE_NUMBERS=',NODES_STRANGE
                        print *, 'VALUES=',STRANGERS
        
                        print *, ' '
                        print *, 'Other variables:'

                        print *, ' PH correction for nitrification:', &
                            (PH_CORR_NITR_NH4(NODES_STRANGE(k), i), k=1,nstrange)
        
                        deallocate(STRANGERS, NODES_STRANGE)
                        error =1
                    end if
                end do ! j-variables
            end do     ! i-layers
            
			if (error .eq. 1) stop
        end if      !debug_stranger



        ! CHECK SETTLING DERIVATIVES
        if(debug_stranger) then
            do i = 1, NUM_SED_LAYERS
                do j = 1, NUM_SIMULATED_SED_VARS
                    if(STRANGERSD(SETTLING_DERIVS(:,i, j),VALUE_strange,nkn) .eq. 1) then
                        nstrange = count(VALUE_strange)
                        allocate(STRANGERS    (nstrange))
                        allocate(NODES_STRANGE(nstrange))
        
                        l=1
                        do k=1,nkn
                            if(VALUE_strange(k)) then
                                STRANGERS    (l) = SETTLING_DERIVS(k,i,j)
                                NODES_STRANGE(l) = k
                                l=l+1
                            end if
                        end do
        
                        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                        print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j !, 'Cell ',CELLNO
                        print *, 'Settling derivative is NaN'
                        print *, 'NODE_NUMBERS=',NODES_STRANGE
                        print *, 'VALUES=',STRANGERS
                        !print *, 'Deriv(i,j)=',SETTLING_DERIVS(:,i, j)
        
                        error =1
                        deallocate(STRANGERS, NODES_STRANGE)
                    end if
                end do ! j-variables
            end do     ! i-layers
			
            if (error .eq. 1) stop
        end if      ! debug_stranger
        
        if (error .eq. 1) stop
        ! END CHECK FOR STRANGERS BEFORE TIME INTEGRATION

        !SUMMING DERIVATIVES AND TIME INTEGRATION
        do i = 1, NUM_SED_LAYERS
            do j = 1, NUM_SIMULATED_SED_VARS

			    ! fixed: transport dervatives were in concentration units
                DERIVS(:,i, j) = TRANSPORT_DERIVS(:,i, j) * SED_DEPTHS(:,i) + & 
                    KINETIC_DERIVS(:,i, j) +  SETTLING_DERIVS (:,i, j)

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

                ! Time integration
                UNIT_AREA_MASSES(:,i, j) = UNIT_AREA_MASSES(:,i, j) + &
                     (DERIVS(:,i, j) * (TIME_STEP / NUM_SUB_TIME_STEPS)) ! fixed bug: derivs were in concetration units

                ! Itegrated unit area to concentrations
                INTERMED_RESULTS(:,i, j) = UNIT_AREA_MASSES(:,i, j) / &
                                                                      SED_DEPTHS(:,i)
                     !(SED_DEPTHS(:,i) * VOLUME_FRACTION(:,i,j)) !bug: multiplication is not necessary because
                                                                 !result should be concetration per total sediment volume

            end do ! j-variables
        end do     ! i-layers       
        
        ! END SUMMING DERIVATIVES AND TIME INTEGRATION
        
        
        ! CHECKING FOR STRANGERS AFTER INTEGRATION
        if(debug_stranger) then
            do i = 1, NUM_SED_LAYERS
                do j = 1, NUM_SIMULATED_SED_VARS
                    if(STRANGERSD(INTERMED_RESULTS(:,i, j),VALUE_strange,nkn) .eq. 1) then
                        nstrange = count(VALUE_strange)
                        allocate(STRANGERS    (nstrange))
                        allocate(STRANGERS1   (nstrange))
                        allocate(STRANGERS2   (nstrange))
                        allocate(NODES_STRANGE(nstrange))
        
                        l=1
                             
				        do k=1,nkn
                            if(VALUE_strange(k)) then
                                STRANGERS    (l) = INTERMED_RESULTS(k,i,j)
                                STRANGERS1   (l) = UNIT_AREA_MASSES(k,i, j)
                                STRANGERS2   (l) = SED_DEPTHS(k,i)
                                NODES_STRANGE(l) = k
                                l=l+1
                            end if
                        end do
					
                        print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j
                        print *, 'Intermediate result after integration is NaN'
                        print *, 'NODE_NUMBERS=',NODES_STRANGE
                        print *, 'VALUES=',STRANGERS
        
                        print *, 'Component variables:'
                        print *, 'VOLUME_FRACTION=',(VOLUME_FRACTION(NODES_STRANGE(k),i,j),k=1,nstrange)
                        print *, 'UNIT_AREA_MASSES=',STRANGERS1
                        print *, 'SED_DEPTHS=',      STRANGERS2
                    
					    deallocate(STRANGERS, NODES_STRANGE)
                        deallocate(STRANGERS1,STRANGERS2)
        
                        error = 1
                    end if  ! strangers checking
                end do ! j-variables
            end do     ! i-layers
        
		    if (error .eq. 1) stop

	    end if !debug_stranger


        ! RECALCULATION OF CONCENTRATIONS BECAUSE OF LAYERS SHIFT DUE TO DEPOSITION OR EROSION

        !if (isedi .ne. 0) then
            ! Deposition and no erosion
            do i = 2, NUM_SED_LAYERS ! without surfice layer
                do j = 1, NUM_SIMULATED_SED_VARS
                    where (H_ERODEP(:) .le. 0.D0)
                        INTERMED_RESULTS_(:,i, j) = &
                            INTERMED_RESULTS(:,i-1, j)*abs(H_ERODEP(:))/SED_DEPTHS(:,i) + &
                            INTERMED_RESULTS(:,i, j)*(SED_DEPTHS(:,i) - abs(H_ERODEP(:)))/SED_DEPTHS(:,i)
                    end where
                end do ! j-variables
            end do     ! i-layers

            ! 1-st layer
            do j = 1, NUM_SED_VARS
                where (H_ERODEP(:) .le. 0.D0)
                    INTERMED_RESULTS_(:,1, j) = &
                        DEPOSIT_CONC(:, j)*(abs(H_ERODEP(:))/SED_DEPTHS(:,1)) + &
                        INTERMED_RESULTS(:,1, j)*(SED_DEPTHS(:,1) - abs(H_ERODEP(:)))/SED_DEPTHS(:,1)
                end where
            end do ! j-variables
            ! End deposition and no erosion

            ! Erosion
            do i = 1, NUM_SED_LAYERS-1 ! without last layer
                do j = 1, NUM_SED_VARS
                    where (H_ERODEP(:) .gt. 0.D0)
                        INTERMED_RESULTS_(:,i, j) = &
                            INTERMED_RESULTS(:,i, j)*(SED_DEPTHS(:,i) - H_ERODEP(:))/SED_DEPTHS(:,i) + &
                            INTERMED_RESULTS(:,i+1, j)*H_ERODEP(:)/SED_DEPTHS(:,i)
                    end where
                end do ! j-variables
            end do     ! i-layers

            ! last layer
            do j = 1, NUM_SED_VARS
                where (H_ERODEP(:) .gt. 0.D0)
                    INTERMED_RESULTS_(:,NUM_SED_LAYERS, j) =  INTERMED_RESULTS(:,NUM_SED_LAYERS, j)

                end where
            end do ! j-variables
            ! End erosion
    
    	!end if !isedi

        ! CHECKING FOR STRANGERS AFTER AFTER LAYERS SHIFT
        if(debug_stranger) then
            do i = 1, NUM_SED_LAYERS
                do j = 1, NUM_SIMULATED_SED_VARS
                    if(STRANGERSD(INTERMED_RESULTS_(:,i, j),VALUE_strange,nkn) .eq. 1) then
                        nstrange = count(VALUE_strange)
                        allocate(STRANGERS    (nstrange))
                        allocate(STRANGERS1   (nstrange))
                        allocate(STRANGERS2   (nstrange))
                        allocate(NODES_STRANGE(nstrange))

                        l=1

    					do k=1,nkn
                            if(VALUE_strange(k)) then
                                STRANGERS    (l) = INTERMED_RESULTS_(k,i,j)
                                STRANGERS1   (l) = UNIT_AREA_MASSES(k,i, j)
                                STRANGERS2   (l) = SED_DEPTHS(k,i)
                                NODES_STRANGE(l) = k
                                l=l+1
                            end if
                        end do
					
                        print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j
                        print *, 'Intermediate result after layers shift is NaN'
                        print *, 'NODE_NUMBERS=',NODES_STRANGE
					    print *, 'VALUES=',STRANGERS

                        print *, 'Component variables:'
                        print *, 'H_ERODEP=', (H_ERODEP(NODES_STRANGE(k)),k=1,nstrange)
                        print *, 'SED_DEPTHS=',      STRANGERS2
                        print *, 'VOLUME_FRACTION=' ,(VOLUME_FRACTION (NODES_STRANGE(k),i,j),k=1,nstrange)
                        print *, 'UNIT_AREA_MASSES=',STRANGERS1
                        print *, 'INTERMED_RESULTS=',(INTERMED_RESULTS(NODES_STRANGE(k),i,j),k=1,nstrange)

                        deallocate(STRANGERS, NODES_STRANGE)
                        deallocate(STRANGERS1,STRANGERS2)

                        error = 1
                    end if  ! strangers checking
                end do ! j-variables
            end do     ! i-layers

            if (error .eq. 1) stop
        end if

        INTERMED_RESULTS(:,:,:) = INTERMED_RESULTS_(:,:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end do !end OF SUB TIME STEPS : CAUTION, FLUXES ARE NOT AVERAGED YET.
    !fixme if substeps are necessary
    !!!!!!!!!!!!!!!!!!!!!!!END LOOP ON SUBSTEPS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !========================FINAL CHEKS AND ASIGNMENTS==============================================
    ! Check for strange values
    if(debug_stranger) then
        do i = 1, NUM_SED_LAYERS
            do j = 1, NUM_SIMULATED_SED_VARS
    
            if(STRANGERSD(INTERMED_RESULTS(:,i, j),VALUE_strange,nkn) .eq. 1) then
                nstrange = count(VALUE_strange)
                allocate(STRANGERS    (nstrange))
                allocate(NODES_STRANGE(nstrange))
    
                l=1
                
				do k=1,nkn
                    if(VALUE_strange(k)) then
                        STRANGERS    (l) = INTERMED_RESULTS(k,i,j)
                        NODES_STRANGE(l) = k
                        l=l+1
                        end if
                    end do
					
                    print *, 'aquabc_sediment1: Layer ', i, 'Variable ',j
                    print *, 'Final state is NaN'
                    print *, 'NODE_NUMBERS=',NODES_STRANGE
                    print *, 'VALUES=',STRANGERS
    
                    deallocate(STRANGERS, NODES_STRANGE)
    
                    !print *, 'FINAL(i,j)=',INTERMED_RESULTS(i, j)
                    !print *, 'time= ', PSTIME
                end if !STRANGERSD
            end do !j - variables
        end do !i - layers
        
		if (error .eq. 1) stop
    end if  !  debug_stranger


    !Final assignments for state var. and var. for the output
    FINAL_SED_STATE_VARS(:,:, :) = INTERMED_RESULTS(:,:,:)    
    SED_OUTPUTS(:,:,1:NUM_SED_VARS) = INTERMED_RESULTS(:,:,:) !State variables

    !Calculation solute concentrations per pore water from total (solid +solutes)  concentration per sed. vol.:
    !Dissolved  NH4N
    SED_OUTPUTS(:,:, NUM_SED_VARS + 1) = &
	    SOLUTE_FRACTIONS(:,:, 1) * INTERMED_RESULTS(:,:, 1)/SED_POROSITIES(:,:)
    !Dissolved  PO4P
    SED_OUTPUTS(:,:, NUM_SED_VARS + 2) = &
	    SOLUTE_FRACTIONS(:,:, 5) * INTERMED_RESULTS(:,:, 5)/SED_POROSITIES(:,:)
    !========================END FINAL CHEKS AND ASIGNMENTS==============================================    
end ! end of sediment routine

!*******************************************************************
!*******************************************************************
INTEGER FUNCTION STRANGER(VALUE)
    ! cheks for NaN and Inf in double precision

    use mod_debug

    implicit none

    DOUBLE PRECISION VALUE, BIGNUMBER, RATIO, LLIMIT, ULIMIT    
    LLIMIT    = -1.0D4
    ULIMIT    =  1.0D4
    BIGNUMBER = 1.0D300
    RATIO     = 1.0D0
    STRANGER  = 0
    
    if (is_nan(VALUE)) then
        STRANGER=1
    end if

    if (is_nan(VALUE)) then
        STRANGER=1
    end if

    !if (dabs(VALUE) >= 1.0D0) then
    !    RATIO = BIGNUMBER/VALUE
    !end if
    
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




!THIS IS A USER PROGRAMMED FUNCTION. IT IS THE USERS RESPONSIBILITY
!TO PROGRAM THE EQUATIONS WHICH GIVE THE MOLECULAR DIFFUSION COEFFICIENTS
!OF STATE VARIABLES IN WATER AS A FUNCTION OF TEMPERATURE AND SALINITY
function SED_MOD_1_ALUKAS_MOLDI_C(SVARNO, T, SAL, TVAR) result(MOL_DIFF)
    implicit none
    !INPUTS:
    !SVARNO : State variable no
    !T      : Temperature in Celsius
    !SAL    : Salinity (ppt)
    !TVAR   : Generic temporal variable
    !OTPUT:
    !MOL_DIFF : Molecular diffusion coefficient, units cm^2/sec inside
    !           and converted to m^2/s for the output

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

    if(SVARNO .gt. 24) then
        print *,'SED_MOD_1_ALUKAS_MOLDI_C:'
        print *,'To get values correctly by Molecular diffusion'
        print *,'BS statevars No should not gt than 24 but is', SVARNO
        stop
    end if

    TK = T + 273.16
    MOL_DIFF = 0.

    ! Shear viscosities
    call SED_MOD_1_CVISC(V25, SS , TS, 1.0D0)
    call SED_MOD_1_CVISC(VTK, SAL, T , 1.0D0)

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
        MOL_DIFF = 1.0D-5
        TS = 25.0D0
        SS = 36.1D0
        MOL_DIFF = MOL_DIFF * (V25 / 298.16D0) * (TK / VTK)
    end if

    !Particulate silicon
    if (SVARNO.eq.12) then
        MOL_DIFF = 0.0D0
    end if

    !Dissolved inorganic carbon (to be corrected with better formulation)
    if (SVARNO.eq.13) then
        MOL_DIFF = 1.0D-6
    end if

    !Alkalinity (to be corrected with better formulation)  !fixme
    if (SVARNO.eq.14) then
        MOL_DIFF = 1.0D-6
    end if

    !Salinity (to be corrected with better formulation)  !fixme
    if (SVARNO.eq.15) then
        MOL_DIFF = 1.0D-6
    end if

    !FEII
    if (SVARNO.eq.16) then
        MOL_DIFF = (3.31D0 + 0.15D0 * T) * 1.0D-5
    end if

    !FEIII
    if (SVARNO.eq.17) then
        MOL_DIFF = 1.0D-6
    end if

    !MNII
    if (SVARNO.eq.18) then
        MOL_DIFF = (3.18D0 + 0.1553D0 * T) * 1.0D-5
    end if


    !MNIV
    if (SVARNO.eq.19) then
        MOL_DIFF = 1.0D-6
    end if


    !Ca
    if (SVARNO.eq.20) then
        MOL_DIFF = (3.60D0 + 0.179D0 * T) * 1.0D-5
    end if

    !Mg
    if (SVARNO.eq.21) then
        MOL_DIFF = (3.43D0 + 0.144D0 * T) * 1.0D-5
    end if

    !S_PLUS_6
    if (SVARNO.eq.22) then
        MOL_DIFF = 4.72D-9 * TK/(35.2D0**0.6) * (V25/VTK)
    end if

    !S_MINUS_2
    if (SVARNO.eq.23) then
       MOL_DIFF = (4.88D0 + 0.232D0 * T) * 1.0D-5
    end if

    !CH4_C
    if (SVARNO.eq.24) then
        MOL_DIFF = 5.7524D-3 * exp(-(3300D0 + 8.94104D-4*(TK-228)**1.5)/(1.9858775*TK))
    end if

    ! Converting to m^2/sec
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
            
    ! Note: routine does not take into account particulate material resuspension
    ! as a flux yet. Everytthing what is resuspended is added to the same fluxes. fixme

    implicit none

    integer :: NUM_FLUXES_FROM_SEDIMENT
    integer :: NUM_FLUXES_TO_ALUKAS

    double precision :: FLUXES_FROM_SEDIMENT(NUM_FLUXES_FROM_SEDIMENT)
    double precision :: FLUXES_TO_ALUKAS    (NUM_FLUXES_TO_ALUKAS)

    if(NUM_FLUXES_TO_ALUKAS .ne. 30) then
        print *,'FLX_SED_MOD_1_TO_ALUKAS_II:'
        print *,'To get values correctly by FLUXES_TO_ALUKAS'
        print *,'num. of statevars should be 33 but is', NUM_FLUXES_TO_ALUKAS
        stop
    end if

    !NH4 FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(1) = FLUXES_FROM_SEDIMENT(1)

    !NO3 FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(2) = FLUXES_FROM_SEDIMENT(2)

    !PO4 FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(3) = FLUXES_FROM_SEDIMENT(5)

    !DISS_OXYGEN FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(4) = FLUXES_FROM_SEDIMENT(8)

    !DIA_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(5) = 0.0D0

    !ZOO_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(6) = 0.0D0

    !ZOO_N FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(7) = 0.0D0

    !ZOO_P FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(8) = 0.0D0

    !DET_PART_ORG_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(9) = 0.0D0

    !DET_PART_ORG_N FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(10) = 0.0D0

    !DET_PART_ORG_P FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(11) = 0.0D0

    !DISS_ORG_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(12) = FLUXES_FROM_SEDIMENT(9)

    !DISS_ORG_N FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(13) = FLUXES_FROM_SEDIMENT(3)

    !DISS_ORG_P FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(14) = FLUXES_FROM_SEDIMENT(6)

    !CYN_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(15) = 0.0D0

    !OPA_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(16) = 0.0D0

    !DISS_Si FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(17) = FLUXES_FROM_SEDIMENT(11)

    !PART_Si FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(18) = FLUXES_FROM_SEDIMENT(12)

    !FIX_CYN_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(19) = 0.0D0

    !INORG_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(20) = FLUXES_FROM_SEDIMENT(13)

    !INORG_C FLUX TO ALUKAS
    FLUXES_TO_ALUKAS(21) = FLUXES_FROM_SEDIMENT(14)

    ! Fluxes to metals in WC
    FLUXES_TO_ALUKAS(22) = FLUXES_FROM_SEDIMENT(16)
    FLUXES_TO_ALUKAS(23) = FLUXES_FROM_SEDIMENT(17)
    FLUXES_TO_ALUKAS(24) = FLUXES_FROM_SEDIMENT(18)
    FLUXES_TO_ALUKAS(25) = FLUXES_FROM_SEDIMENT(19)

    FLUXES_TO_ALUKAS(26) = FLUXES_FROM_SEDIMENT(20)
    FLUXES_TO_ALUKAS(27) = FLUXES_FROM_SEDIMENT(21)
    FLUXES_TO_ALUKAS(28) = FLUXES_FROM_SEDIMENT(22)
    FLUXES_TO_ALUKAS(29) = FLUXES_FROM_SEDIMENT(23)
    FLUXES_TO_ALUKAS(30) = FLUXES_FROM_SEDIMENT(24)
end subroutine FLX_SED_MOD_1_TO_ALUKAS_II

!************************************************************
!************************************************************


! ----------------------------------------------------------------------------!
! Subroutine that calculates the PH_CORRECTION factor according to            !
! the euqations (Eq. 162 and Eq. 163) in AQUTOX Relase 3.1 Plus documentation !
! (EPA-829-R-14-007)                                                          !
! Specialised for  for BS  by Petras                                          !
! --------------------------------------------------------------------------- !
!                                DOCUMENTATION                                !
!                                                                             !
! INPUTS :                                                                    !
! --------                                                                    !
! PH_MIN  : Minimum value of pH for the optimum range (scalar)                !
! PH_MAX  : Maximum value of pH for the optimum range (scalar)                !
! PH      : PH values for which correction factor will be calculated (matrix) !
! nkn     : Number of rows in matrices PH and PH_CORR                         !
! nlay    : Number of columns in in matrices PH and PH_CORR                   !
!                                                                             !
! OUTPUTS                                                                     !
! PH_CORR : pH correction factors calculated for PH (matrix)                  !
! --------------------------------------------------------------------------- !
!                    Development date : 31st August 2015                      !
!                                                                             !
!                            Developer : Ali Ertürk                           !
! --------------------------------------------------------------------------- !
subroutine CALCULATE_PH_CORR_SED(PH_CORR, PH, PH_MIN, PH_MAX, nkn, nlay)

    use AQUABC_II_GLOBAL
    implicit none

    ! Ingoing variables
    real(kind = DBL_PREC), dimension(nkn, nlay), intent(in) :: PH
    real(kind = DBL_PREC), intent(in) :: PH_MIN
    real(kind = DBL_PREC), intent(in) :: PH_MAX
    integer, intent(in) :: nkn, nlay
    ! Outgoing variables
    real(kind = DBL_PREC), dimension(nkn, nlay), intent(inout) :: PH_CORR

    integer i,j, error

    !real(kind = DBL_PREC), dimension(nkn) :: PH_MIN_ARRAY
    !real(kind = DBL_PREC), dimension(nkn) :: PH_MAX_ARRAY

    error = 0


    PH_CORR = 1


    !return
    !PH_MIN_ARRAY(:) = PH_MIN
    !PH_MAX_ARRAY(:) = PH_MAX

    do j=1,nlay

      where(PH(:,j) < PH_MIN)
          PH_CORR(:,j) = exp(PH(:,j) - PH_MIN)
      end where

      where(PH(:,j) > PH_MAX)
          PH_CORR(:,j) = exp(PH_MAX - PH(:,j))
      end where

    end do

    do i=1,nkn
        do j=1,nlay
            if (PH_CORR(i,j) .gt. 1.D0 .or. PH_CORR(i,j) .lt. 0.D0) then
                print *, 'CALCULATE_PH_CORR: Incorrect PH_CORR'
                print *, 'PH_CORR',PH_CORR(i,j)
                print *, 'PH_MIN:', PH_MIN
                print *, 'PH_MAX:', PH_MAX
                print *, 'PH:', PH(i,j)
                print *, 'Internal node number:',i
                print *, 'Layer number:',j
                print *, 'Number of layers:',nlay
                error =1
            end if
        end do
    end do

   if(error .eq. 1)stop

end subroutine CALCULATE_PH_CORR_SED

!************************************************************************
!************************************************************************
