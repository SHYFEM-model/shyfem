
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

!module aquabc_pel_state_var_indexes

module aquabc_pel_state_var_indexes
    implicit none
    integer, parameter :: NH4_N_INDEX          = 1
    integer, parameter :: NO3_N_INDEX          = 2
    integer, parameter :: PO4_P_INDEX          = 3
    integer, parameter :: DISS_OXYGEN_INDEX    = 4
    integer, parameter :: DIA_C_INDEX          = 5
    integer, parameter :: ZOO_C_INDEX          = 6
    integer, parameter :: ZOO_N_INDEX          = 7
    integer, parameter :: ZOO_P_INDEX          = 8
    integer, parameter :: DET_PART_ORG_C_INDEX = 9
    integer, parameter :: DET_PART_ORG_N_INDEX = 10
    integer, parameter :: DET_PART_ORG_P_INDEX = 11
    integer, parameter :: DISS_ORG_C_INDEX     = 12
    integer, parameter :: DISS_ORG_N_INDEX     = 13
    integer, parameter :: DISS_ORG_P_INDEX     = 14
    integer, parameter :: CYN_C_INDEX          = 15
    integer, parameter :: OPA_C_INDEX          = 16
    integer, parameter :: DISS_Si_INDEX        = 17
    integer, parameter :: PART_Si_INDEX        = 18
    integer, parameter :: FIX_CYN_C_INDEX      = 19
    integer, parameter :: INORG_C_INDEX        = 20
    integer, parameter :: TOT_ALK_INDEX        = 21
    integer, parameter :: FE_II_INDEX          = 22
    integer, parameter :: FE_III_INDEX         = 23
    integer, parameter :: MN_II_INDEX          = 24
    integer, parameter :: MN_IV_INDEX          = 25
    integer, parameter :: CA_INDEX             = 26
    integer, parameter :: MG_INDEX             = 27
    integer, parameter :: S_PLUS_6_INDEX       = 28
    integer, parameter :: S_MINUS_2_INDEX      = 29
    integer, parameter :: CH4_C_INDEX          = 30 
end module AQUABC_PEL_STATE_VAR_INDEXES
