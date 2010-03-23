!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ModuleConstants
!
! DESCRIPTION
!   Full list of Fortran parameters  ( comparable with Sesame constants)

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  MODULE constants

!
! !USES:

  USE global_mem, ONLY:RLEN, ZERO

!  
!
! !AUTHORS
!    mfstep/ERSEM team
!
! !REVISION_HISTORY
!   ---
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Implicit typing is never allowed
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  IMPLICIT NONE
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Default all is public here
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  public
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Global Constants
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  real(RLEN),   parameter :: DBL_MIN=2.2250738585072014e-308_RLEN
  integer,      parameter :: ECOLOGY=1  ! Base temperature for Q10
  integer,      parameter :: TRANSPORT=2
  real(RLEN),   parameter :: MW_C=12.0_RLEN  ! Molecular weight Carbon
  real(RLEN),   parameter :: MW_N=14.0_RLEN  ! Molecular weight Nitrogen
  real(RLEN),   parameter :: MW_P=31.0_RLEN  ! Molecular weight Phosphorus
  real(RLEN),   parameter :: MW_SI=28.0_RLEN  ! Molecular weight Silica
  real(RLEN),   parameter :: E2W=0.217_RLEN  ! Conversion factor Einstein->W
  real(RLEN),   parameter :: SEC_PER_DAY=86400.0_RLEN  ! Seconds in day
  real(RLEN),   parameter :: ONE_PER_DAY=1.0_RLEN  ! rate which is used in cases where implicitely is assumed
  !  is used and which are found.
  integer,      parameter :: NO_BENTHOS=0
  integer,      parameter :: BENTHIC_RETURN=1
  integer,      parameter :: BENTHIC_BIO=2
  integer,      parameter :: BENTHIC_FULL=3
  real(RLEN),   parameter :: HOURS_PER_DAY=24.0_RLEN  ! Hours in a day
  real(RLEN),   parameter :: SOLAR_RADIATION=1368.0_RLEN  ! Solar radiation constant
  integer,      parameter :: NODATA=0
  integer,      parameter :: DAILYDATA=1
  integer,      parameter :: CLOUDDATA=3
  integer,      parameter :: ALLDATA=4
  integer,      parameter :: ONLYDAILY=5
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 	Integer constants used in the benthic nutrient dynamics
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer,      parameter :: DEFINE=-1000
  integer,      parameter :: DOUBLE_DEFINE=-1100
  integer,      parameter :: PARAMETER_DEFINE=-1200
  integer,      parameter :: EQUATION=0
  integer,      parameter :: SDERIVATIVE=-2
  integer,      parameter :: DERIVATIVE=-1
  integer,      parameter :: INTEGRAL=1
  integer,      parameter :: EXPONENTIAL_INTEGRAL=11
  integer,      parameter :: SHIFT=21
  integer,      parameter :: LINEAR_TERM=1
  integer,      parameter :: CONSTANT_TERM=0
  integer,      parameter :: QUADRATIC_TERM=2
  integer,      parameter :: EXPONENTIAL_TERM=-1
  integer,      parameter :: ZERO_EXPONENTIAL_TERM=-2
  integer,      parameter :: BESSELI_EXP_TERM=-5
  integer,      parameter :: BESSELK_EXP_TERM=-6
  integer,      parameter :: NEW_BOUNDARY=999
  integer,      parameter :: ADD=1000
  integer,      parameter :: SUBTRACT=2000
  integer,      parameter :: INPUT_TERM=6001
  integer,      parameter :: START_ADD_TERM=6002
  integer,      parameter :: INPUT_ADD_TERM=6000
  integer,      parameter :: INPUT_SUBTRACT_TERM=7000
  integer,      parameter :: RFLUX=1
  integer,      parameter :: MASS=2
  integer,      parameter :: AVERAGE=3
  integer,      parameter :: PARAMETER=4
  integer,      parameter :: STANDARD=0
  integer,      parameter :: SET_CONTINUITY=100
  integer,      parameter :: SET_LAYER_INTEGRAL=200
  integer,      parameter :: SET_LAYER_INTEGRAL_UNTIL=300
  integer,      parameter :: SET_BOUNDARY=400
  integer,      parameter :: SET_DEPTH_INTEGRAL=500
  integer,      parameter :: TO_DEPTH_INTEGRAL=-1
  integer,      parameter :: GET=9000
  integer,      parameter :: LABDA_1=1
  integer,      parameter :: LABDA_2=2
  integer,      parameter :: COEFFICIENT=3
  integer,      parameter :: LAYERS=4
  integer,      parameter :: DIFFUSION=5
  integer,      parameter :: POROSITY=6
  integer,      parameter :: ADSORPTION=7
  integer,      parameter :: INITIALIZE=0
  integer,      parameter :: FLAG=1
  integer,      parameter :: METHOD=2
  integer,      parameter :: LAYER1=1
  integer,      parameter :: LAYER2=2
  integer,      parameter :: LAYER3=3
  integer,      parameter :: LAYER4=4
  integer,      parameter :: FOR_ALL_LAYERS=-1
  integer,      parameter :: NUMBER_OF_PROFILES=8
  integer,      parameter :: NCOEFF=12
  integer,      parameter :: NLAYER=5
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  end module
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
