
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
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

! ======================================================================
!		CONSTANTS FOR SEDI3D
! ======================================================================
! Common block containing constants used for sedi3d

        INTEGER NLBDIM            !number of bed layer
        PARAMETER (NLBDIM=20)

        INTEGER NSDIM             !number of grainsize classes
        PARAMETER (NSDIM=9)


        COMMON /CSEDI3D/ KCOES,LIMCOH,SMOOTH,ANGREP,MORPHO,RHOSED,
     &     POROS,SURFPOR,IOPT,NBCC
        SAVE /CSEDI3D/

        DOUBLE PRECISION KCOES	!Fraction of mud for sediment to be cohesive [0-1]
        DOUBLE PRECISION LIMCOH	!grainsize limit for cohesive sediment [m]
        DOUBLE PRECISION SMOOTH	!smoothing factor for morphodynamic [0-1]
        DOUBLE PRECISION ANGREP	!angle of repose [rad]
        DOUBLE PRECISION MORPHO	!Morphological factor
        DOUBLE PRECISION RHOSED	!sediment grain density
        DOUBLE PRECISION POROS	!bed porosity [0,1]
        DOUBLE PRECISION SURFPOR !surficial porosity [0,1]
        INTEGER IOPT		!SEDIMENT TRANSPORT FORMULA OPTION NUMBER
        INTEGER NBCC		!same as NBCONC but equal 0 if no cohesive

! Common block used in SEDTRANS05
!
! ======================================================================
!                         MAIN
! ======================================================================

      DOUBLE PRECISION G	!GRAVITY
      DOUBLE PRECISION PII	!

      COMMON /MCONST/ G,PII
      SAVE /MCONST/

! ======================================================================
!		CONSTANTS FOR COHESIVE SEDIMENTS
! ======================================================================
! Common block containing several constant used for the cohesive part of
! sedtrans. It MUST be initialised with an call to INICONST before the first
! call to COHESIVE (or to first call to SEDTRANS05 with IOPT1=7)

      COMMON /CCONST/ CSULVA,TMULVA,TRULVA,RKERO,E0,CDISRUPT,CLIM1,
     &   CLIM2,KFLOC,MFLOC,RHOCLAY,CTAUDEP,PRS,RHOMUD,DPROFA,DPROFB,
     &   DPROFC,DPROFD,DPROFE,CONSOA,TEROA,TEROB,TEROC,TEROD,DISTR,
     &   CDRAGRED,Z0COH,FCWCOH,TAU0EFF,
     &   NDISTR,WSCLAY,MEDDISTR,DOCOMPACT

      SAVE /CCONST/

      DOUBLE PRECISION CSULVA  ! coefficient for the solid transmitted stress by Ulva
      DOUBLE PRECISION TMULVA  ! threshold of motion of   Ulva (Pa)
      DOUBLE PRECISION TRULVA  ! threshold of full resuspension of Ulva (Pa)
      DOUBLE PRECISION RKERO   ! Erosion proportionality coefficient
      DOUBLE PRECISION E0      ! Minimum erosion rate
      DOUBLE PRECISION CDISRUPT! constant for turbulent floc disruption during erosion
      DOUBLE PRECISION CLIM1   ! lower limit for flocculation (kg/m**3)
      DOUBLE PRECISION CLIM2   ! limit between simple and complex flocculation equation (kg/m**3)
      DOUBLE PRECISION KFLOC   ! constant K for flocculation equation
      DOUBLE PRECISION MFLOC   ! constant M for flocculation equation
      DOUBLE PRECISION RHOCLAY ! density of clay mineral
      DOUBLE PRECISION CTAUDEP ! scaling factor for TAUCD
      DOUBLE PRECISION PRS     ! Resuspension probability (range 0-1)
      DOUBLE PRECISION RHOMUD  ! Density of the freshly deposited mud
      DOUBLE PRECISION DPROFA,DPROFB,DPROFC,DPROFD,DPROFE
                               ! constants for final density profile
      DOUBLE PRECISION CONSOA ! time constant of consolidation
      DOUBLE PRECISION TEROA,TEROB,TEROC,TEROD
                               ! constants for erosion threshold from density and overlaying mass
      INTEGER WSCLAY           ! primary median Ws class (must be in the range 1:NBCONC)
      INTEGER MEDDISTR         ! median position of DISTR
      INTEGER DOCOMPACT        ! if not zero, call COMPACT from COHESIVE

      INTEGER NDISTR           ! number of elements used in DISTR
      DOUBLE PRECISION DISTR(NSDIM)! (normal) distribution of sediment put in suspension
                               ! the sum of DISTR(1:NDISTR) MUST be exactly 1.0
      DOUBLE PRECISION CDRAGRED! constant for the drag reduction formula

! Common block containing the settling velocity of each Ws class.
! It MUST be initialised with an call to INICONST before the first
! call to COHESIVE (or to first call to SEDISDGM with IOPT1=7)
      COMMON /WSCLASS/ WSI
      SAVE /WSCLASS/
      DOUBLE PRECISION WSI(NSDIM) ! settling velocity of each Ws class (m/s)

! For subroutine FRICFAC (previous common bloc FRICC
      DOUBLE PRECISION Z0COH    !BED ROUGHNESS LENGHT FOR COHESIVE SEDIMENTS
      DOUBLE PRECISION FCWCOH   !FRICTION FACTOR FOR COHESIVE SEDIMENTS

! Additional result from COHESIVE
      DOUBLE PRECISION TAU0EFF  !effective bed shear stress used in COHESIVE
                                !includes drag reduction and solid transm. stress by Ulva


