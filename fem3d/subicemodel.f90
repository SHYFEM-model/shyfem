
!--------------------------------------------------------------------------
!
!    Copyright (C) 2007  Letizia Tedesco
!    Copyright (C) 2019-2020  Rasa Izelyte
!    Copyright (C) 2019-2020  Georg Umgiesser
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

! notes:

! routines written originally by Letizia Tedesco (letizia.tedesco@cmcc.it)
!    in January 2007, Last revised February 2019
! the routines below have been extracted from the original matlab code
!    and converted to Fortran by Rasa Idzelyte in 2020. 

!--------------------------------------------------------------------------
!
! List of References about ESIM2 model
!
! Tedesco, L., Vichi, M., Haapala, J., and Stipa, T., 2009. An enhanced sea 
! ice thermodynamic model applied to the Baltic Sea. Boreal Environmental 
! Research, 14: 68-S80.
!
! Tedesco, L., Vichi, M., Haapala, J., and Stipa, T., 2010. A dynamic 
! Biologically-Active Layer for numerical studies of the sea ice ecosystem. 
! Ocean Modelling, 35(1-2):89-104, doi:10.1016/j.ocemod.2010.06.008.
!
!--------------------------------------------------------------------------

! contents :
!
! callable routines from outside
!
! ice_check_nvars(nvars)	!check number of internal variables
!		nvars		!total number of internal variables
! ice_init_vars(vars)		!initialize internal variables
!		vars(nvars)	!internal variables
!
! ice_run(CL,Fsd_cloud,P_rate,qa,qs,Ta,Ua,Sw,deltat,i,vars,iceth) !run ice model
!               Cl - total cloud cover fraction [fraction 0-1]
!               Fsd_cloud - downward solar radiation flux [W m-2]
!               P_rate - precipitation rate [mm s-1]
!               qa - specific humidity of the air [non-dim]
!               qs - specific humidity of the surface [non-dim]
!               Ta - air temperature at 10 m height [K]
!               Ua - wind speed at 10 m height [m s-1]
!               Sw - seawater salinity [psu]
!               deltat - time step [s]
!               i - ice model iteration number (not used)
!               vars - internal variables
!               iceth - ice thickness [m]

! version :
!	1.5

! revision log :
!
! 15.03.2020	ggu&riz	started converting matlab code
! 01.07.2020    ggu&riz started running first tests
! 15.09.2020    riz	small changes (albedo, etc..)
! 11.11.2020    ggu     eliminating compiler errors
! 13.11.2020    ggu     changed number of variables to 46
! 13.11.2020    ggu     documentation fix: Cl is a fraction, not a percentage
! 14.11.2020    ggu&riz compute freezing point temperature
! 15.11.2020    ggu&riz new module ice_global
! 17.11.2020    ggu&riz bug fix freezing temperature
! 22.11.2020    ggu	added code for debugging ice model

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! start of ice routines - do not change anything below here
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------


!---------------------------------------------------------------------------
! MODEL PROGNOSTIC VARIABLES
!---------------------------------------------------------------------------

!   deltahi_bott        (m)             variation in thickness due to sea ice melting at the bottom
!   deltahi_melt_surf   (m)             variation in thickness due to sea ice melting at the surface
!   deltahmi_melt_surf  (m)             variation in thickness due to snow ice/supeimposed ice melting at the surface
!   deltahs_melt_surf   (m)             variation in thickness due to snow melting at the surface
!   hs                  (m)             1 layer of new snow + 1 layer old snow
!   hmi                 (m)             1 layer of snow ice + 1 layer of superimposed ice
!   hi                  (m)             1 layer of sea ice
!   T0                  (K)             surface temperature
!   Sbr                 (per mill)      brines salinity
!   Sice                (per mill)      sea ice bulk salinity
!   Vbr_ice             (non-dim)       sea ice brines volume fraction

!---------------------------------------------------------------------------
! MODEL FORCING VARIABLES 
!---------------------------------------------------------------------------

!   Cl          (non-dim)   total cloud cover fraction [0-1]
!   Fsd_cloud   (W m-2)     downward solar radiation flux
!   P_rate      (mm s-1)    precipitation rate !::RASA: changed from m/s to mm/s
!   qa          (non-dim)   specific humidity of the air
!   qs          (non-dim)   specific humidity of the surface
!   Ta          (K)         air temperature at 10 m height
!   Ua          (m s-1)     wind speed at 10 m height

!---------------------------------------------------------------------------
! MODEL DIAGNOSTIC VARIABLES
!---------------------------------------------------------------------------

!   alpha           (non-dim)       albedo
!   cpa             (J kg^-1 K^-1)  sea ice heat capacity
!   delta_bucket    (m)             variations in snow in the bucket model
!   ea              (non-dim)       air relative humidity
!   Fbott           (W m-2)         conductive flux at the bottom
!   Fl              (W m-2)         net longwave radiation flux
!   Fla             (W m-2)         net latent heat flux
!   Fld             (W m-2)         downward longwave radiation flux
!   Flu             (W m-2)         upward longwave radiation flux
!   Fs              (W m-2)         net solar radiation flux
!   Fse             (W m-2)         net sensible heat flux
!   Fsurf           (W m-2)         conductive flux at the surface
!   Fw              (W m-2)         oceanic heat flux
!   hs_prec         (m)             snow precipitation
!   hs_prec_bucket  (m)             snow precipitation in the bucket model
!   ice_fr          (non-dim)       ice fraction: 0.0 if open ocean, 1.0 if sea ice cover
!   I0              (W m-2)         fraction of solar radiation penetrating the surface
!   IM              (W m-2)         fraction of solar radiation penetrating snow ice/supeimposed ice
!   ISI             (W m-2)         fraction of solar radiation penetrating sea ice below the first 10 cm
!   ISI_10          (W m-2)         fraction of solar radiation penetrating sea ice in the above first 10 cm
!   ISI_layer       (W m-2)         fraction of solar radiation penetrating the biologically-active layer sea ice in the above first 10 cm
!   ISI_bio         (W m-2)         fraction of solar radiation penetrating sea ice and reaching the middle point of the biologically active layer
!   ki              (W m^-1 s^-1)   thermal conductivity of sea ice
!   kmi_clear       (m-1)           snow ice/superimposed ice extiction coeff in clear sky conditions
!   kmi_cloudy      (m-1)           snow ice/superimposed ice extiction coeff in cloudy sky conditions
!   kmi_av          (m-1)           mean snow ice/superimposed ice extiction coeff as function of Cl
!   ks              (W m^-1 s^-1)   thermal conductivity of snow
!   ks_clear        (m-1)           snow extiction coeff in clear sky conditions
!   ks_cloudy       (m-1)           snow extiction coeff in cloudy sky conditions
!   ks_snow          (m-1)           mean snow extiction coeff as function of Cl
!   ksi_10_clear    (m-1)           sea ice extiction coeff at 10 cm in clear sky conditions
!   ksi_10_cloudy   (m-1)           sea ice extiction coeff at 10 cm in cloudy sky conditions
!   ksi_10_av       (m-1)           mean sea ice extiction coeff at 10 cm as function of Cl
!   ksi_clear       (m-1)           sea ice extiction coeff in clear sky conditions
!   ksi_cloudy      (m-1)           sea ice extiction coeff in cloudy sky conditions
!   ksi_av          (m-1)           mean sea ice extiction coeff as function of Cl
!   roi_surf        (kg m-3)        density of pure ice at the surface
!   ro_sice_surf    (kg m-3)        bulk density of sea ice at the surface as function of T,S;
!   ro_sice_bott    (kg m-3)        bulk density of sea ice at the bottm (density of pure ice at the bottom is fixed at 0.900)
!   roi_av          (kg m-3         density of pure ice of the layer as function of Tice_av in g/cm^3
!   ro_sice_bulk    (kg m-3)        bulk density of ice as function of Tice,Sice in g/cm^3
!   ro_br           (kg m-3)        brines density as function of brines salinity in g/cm^3
!   ro_br_5         (kg m-3)        brines density in the biologically active layer as function of brines salinity in g/cm^3
!   ro_br_bio       (kg m-3)        mean brines density in the biologically active layer as function of brines salinity in g/cm^3
!   T_beg           (K)             initial guess for temperature iteration
!   T_iter          (K)             iteration temperature
!   Tice            (C)            temperature of sea ice at the upper interface
!   Tice_5          (C)            temperature of sea ice at brines volume 0.05
!   Tice_bio        (C)            temperature of sea ice in the biological active layer
!   Tis             (K)             temperature at the interf snow ice (or superimposed ice)/sea ice
!   Tss             (K)             temperature at the interf snow/snow ice (or superimposed ice)

!---------------------------------------------------------------------------
! SET PARAMETRES AND CONSTANTS 
!---------------------------------------------------------------------------

!===========================================================================
module ice_params
!===========================================================================

implicit none

!double precision, parameter :: deltat=5400;
double precision, parameter :: hs_min=0.06D+0;            ! minimum snow thickness (m)                             !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
double precision, parameter :: hmi_min=0.05D+0;           ! minimum snow ice/supeimposed ice thickness (m)         !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
double precision, parameter :: hi_min=0.05D+0;            ! minimum sea ice thickness (m)                          !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
double precision, parameter :: ZERO=1.D-5;             ! ZERO OF THE MODEL                                      !!!BE AWARE OF THE SENSITIVITY OF THE MODEL TO THIS LIMIT!!!
!Fwater=0.0;             ! = Fw = 7.0 oceanic heat fluxes (W m-2)                     !!!THIS STRONGLY DEPENDS ON THE LOCATION!!!

double precision, parameter :: Cl_init=0.5D+0;         !  initial cload cover

double precision, parameter :: alpha_ow=0.06D+0;          !   (non-dim)       seawater albedo                                           ::RASA: changed, was 0.15
!double precision, parameter :: alpha_mi=0.55;          !   (non-dim)       intermediate layer albedo
double precision, parameter :: beta=2.0D+0;!2.667;		!   (non-dim)       empir coeff from snow ice to ice (after Lepparanta,1983)
double precision, parameter :: beta2=0.13D+0;             !   (W m-1)         emp cost (after Untersteiner, 1961)        ::RASA: in Tedesco thesis Untersteiner is 0.17, but nothing changes
double precision, parameter :: c0=2093.0D+0; 				!   (J(kg*deg)-1)   specific heat of fresh ice
double precision, parameter :: c1=21.875D+0;				!   (K)             emp.costant
double precision, parameter :: c2=7.66D+0;				!   (K)             emp.costant
double precision, parameter :: c10=0.1D+0;                !   (m)             (after Flato and Brown, 1996 for landfast sea ice of the Arctic)
double precision, parameter :: c11=0.44D+0;               !   (m^-0.28)       (after Flato and Brown, 1996 for landfast sea ice of the Arctic)
double precision, parameter :: c12=0.075D+0;              !   (m^-2)          (after Flato and Brown, 1996 for landfast sea ice of the Arctic)
double precision, parameter :: cail=1.75*1.D-3;		!   (non-dim)       bulk transfer coefficient
double precision, parameter :: cais=1.20*1.D-3;		!   (non-dim)       bulk transfer coefficient
double precision, parameter :: cpair=1004.0D+0;           !   (J kg^-1 K^-1)  air heat capacity
double precision, parameter :: cpw=4186.0D+0;             !   (J kg^-1 K^-1)  specific heat of sea water
!cs=2093.0;				!   (J (kg K)-1)    specific heat of snow                                        ::RASA: not found further in code
double precision, parameter :: cwih=2.8*1.D-4;		!   (non-dim)       bulk heat trans coeff (after Omstedt and Wettlaufer,1992)    ::RASA: commented out in OCEANIC FLUXES section
double precision, parameter :: emis=0.97D+0;				!   (non-dim)       water/ice emissivity
double precision, parameter :: emp=0.62197D+0;            !   (mbar)          empirical coefficient for specific humidity computation
!gamm=18.0;              !   (J (K*kg)-1)    emp coeff (after Bitz and Lipscomb, 1999; otherwise 17.2 after Untersteiner, 1961)     ::RASA: not found further in code
double precision, parameter :: infra_fr=0.77D+0;          !   (non-dim)       Infrared partition as for water type II (after Jerlov, 1968)
double precision, parameter :: k0=2.03D+0;                !   (W m^-1 K^-1)   fresh ice thermal conductivity
!ks_prec=0.0641;         !   (W m^-1 K^-1)   precipit thermal conductivity                             ::RASA: not found further in code
!kw=0.563;               !   (W m^-1 K^-1)   sea water thermal conductivity (after Lepparanta, 1983)   ::RASA: not found further in code
double precision, parameter :: kkmi=0.9D+0;				!   (W m^-1 K^-1)   snow ice thermal conductivity (after Lepparanta, 1983)
double precision, parameter :: k_ocean_red=0.6667D+0;     !   (m^-1)          seawater extiction coefficient (infrared, water type II, after Jerlov 1968)
double precision, parameter :: k_ocean_vis=0.0714D+0;     !   (m^-1)          seawater extiction coefficient (visible, water type II, after Jerlov 1968)
double precision, parameter :: h_mix_const=3.0D+0;              !   (m)             mixed layer depth                        ::RASA: changed, was 10.0
double precision, parameter :: Lai=2.501*1.D6;			!   (J kg^-1)       latent heat of vaporization
!L0i=297000.0;           !   ()              latent heat of fusion of pure ice             ::RASA: not found further in code
!L0m=315000.0;           !   ()              latent heat of fusion of meteoric ice         ::RASA: not found further in code
double precision, parameter :: L0s=334000.0D+0;           !   ()              latent heat of fusion of pure snow
double precision, parameter :: Lv=2.501*1.D6;          !   (J kg-1)        latent heat of vaporization of fresh water at 273.15 K
double precision, parameter :: mu=0.054D+0; 				!   (C)             ratio between Tfr and brines salinity (after Assur, 1958)
double precision, parameter :: P0=1013.0D+0;				!   (mbar)          seawater pressure at the surface
!q1=1.16378*10^7;        !   (kg m-3)        emp coeff                    ::RASA: not found further in code
!q2=5897.8;              !   (K)             emp coeff                    ::RASA: not found further in code
double precision, parameter :: roa=1.225D+0;              !   (kg m-3)        air density
!roi_bott=0.9;           !   (Kg m-3)        pure sea ice density at the bottom   ::RASA: shouldnt it be 900?, commented out further in the code
double precision, parameter :: romi=850.0D+0;             !   (kg m-3)        snow ice/superimposed ice density
double precision, parameter :: rosa=1.68*1.D-5;        !   1.68*10^-7;    !5.88*10^-8;        !   (C??^-1*m*s^-1)  emp coeff (after Vancoppenolee et al, 2007; otherwise 1.68*10^-7 after Cox and Weeks, 1988)
double precision, parameter :: ros_prec=300D+0;           !   (kg m-3)        snow precip density             ::RASA: changed, was 250
double precision, parameter :: row=1000.0D+0;             !   (kg m-3)        water density                                              ::RASA: changed, was 1026.0
double precision, parameter :: Sbr_end=18.5185D+0;        !   (Kg m-3)        pure brines density at the bottom                          
double precision, parameter :: sigma=5.68*1.D-8;		!   (W m^-2 K)      Stefan-Boltzmann constant
double precision, parameter :: si_fract=0.2D+0;           !   (non-dim)       fraction of snow transformation in snow ice !!!IT CAN BE SET UP BY THE USER DEPENDING ON SNOW TYPE, ECC.)!!!   ::RASA: changed, was 0.0
double precision, parameter :: ss_fract=0.2D+0;           !   (non-dim)       fraction of snow transformation in superimposed ice (after Cheng et al, 2006)!!!IT CAN BE SET UP BY THE USER DEPENDING ON SNOW TYPE, ECC.)!!!   ::RASA: changed, was 0.0
double precision, parameter :: Ssnow=0.0D+0;              !   (per mill)      snow salinity
double precision, parameter :: Tb=273.155D+0;				!   (K)             empirical temp                                             ::RASA: it is used fot latent heat computations
double precision, parameter :: Va=0.015D+0;               !   ()              gas volume fraction in the ice is fixed       ::RASA: used for bulk density of ice and thermal conductivity, if set to 0.0, the snow thickness R slightly increases
double precision, parameter :: viola=20.0D+0;             !                   for salinity compuatation (after  Vancoppenolee et al, 2007)
double precision, parameter :: vis_fr=0.8D+0;             !                   fraction of visible light penetrating the surface layer when SIM is on   ::RASA: changed, was 0.7

double precision, parameter :: Tfreez=273.145D+0;         ! = Tfr = -0.4 C (271.4220) seawater freezing temperature (K)       !!!THIS STRONGLY DEPENDS ON THE LOCATION!!! ::RASA: for 0C = 273.15, but get an error for Sice_bott
!double precision, parameter :: Tfr=Tfreez;         !   (K)             freezing temp of sea ice
double precision, parameter :: Tfrs=273.15D+0;			!   (K)             freezing temp of snow/fresh water  ::RASA: DO NOT CHANGE THIS VALUE IT'S USED IN CODE!!!
double precision, parameter :: Tfrs_appr=273.1499D+0;

!===========================================================================
end module ice_params
!===========================================================================


!===========================================================================
module ice_variables
!===========================================================================

use ice_params

integer, parameter :: nvars = 46

double precision, save :: Tmix(2)=284.00D+0;                     ! Initial mixed layer temperature
double precision, save :: T0(2)=284.00D+0;                         ! Initial surface temperature
double precision, save :: Ts(2)=Tfreez;                          ! Initial snow temperature
double precision, save :: Tsi(2)=Tfreez;                         ! Initial snow ice temperature
double precision, save :: Ti(2)=Tfreez;                          ! Initial sea ice temperature

double precision, save :: hi(2)=ZERO;                         ! Initial sea ice thickness
double precision, save :: hs(2)=ZERO;                         ! Initial snow thickness
double precision, save :: hmi(2)=ZERO;                        ! Initial snow ice/superimposed ice thickness
double precision, save :: hs_prec_bucket(2)=0.0;              ! Initial snow thickness in the bucket model

double precision, save :: snow_fr(2)=0.0;
double precision, save :: sea_water_fr(2)=1.0;

double precision, save :: Vbr_ice(2)=0.0;                     ! Initial brines volume
double precision, save :: Sbr_ice(2)=0.0;                     ! Initial brines salinity
double precision, save :: Sice(2)=0.0;                        ! Initial sea ice bulk salinity
double precision, save :: Sice_bott(2)=0.0;
double precision, save :: Ssnowice(2)=0.0;
double precision, save :: Vbr_i(2)=0.0;                       ! Initial internal brines volume
double precision, save :: Sbr_i(2)=0.0;                       ! Initial internal brines salinity
double precision, save :: Si(2)=0.0;                          ! Initial internal sea ice bulk salinity

double precision, save :: Sice_5(2)=0.0;                      ! Initial sea ice bulk salinity at Vbr=0.05
double precision, save :: Sbr_5(2)=0.0;                       ! Initial brines salinity at Vbr=0.05
double precision, save :: Tice_5(2)=Tfreez-Tfrs;               ! Initial sea ice temperature at Vbr=0.05
double precision, save :: hi_5(2)=ZERO;                       ! Initial sea ice thickness at Vbr=0.05

double precision, save :: Sice_bio(2)=0.0;                    ! Initial sea ice bulk salinity in the bio layer
double precision, save :: Sbr_bio(2)=0.0;                     ! Initial brines salinity in the bio layer
double precision, save :: hi_bio(2)=ZERO;                  ! Initial sea ice thickness in the bio layer
double precision, save :: Vbr_bio(2)=0.0;                     ! Initial brines volume in the bio layer

double precision, save :: R(2)=0.0;

double precision, save :: ro_sice_bott(2)=0.917D+0;              ! Initial sea ice bulk density at the bottom
double precision, save :: ro_sice_surf(2)=0.900D+0;              ! Initial sea ice bulk density at the surface
double precision, save :: ro_sice_bulk(2)=0.910D+0;

double precision, save :: ki_5_bott(2)=2.0;                   ! Initial sea ice thermal conductivity between Vbr=0.05 and the bottom
double precision, save :: ki_ice_5(2)=2.0;                    ! Initial sea ice thermal conductivity between sea ice surface and Vbr=0.05
double precision, save :: ki_ice_bio(2)=2.0;
double precision, save :: ki_bio_bott(2)=2.0;                 ! Initial sea ice thermal conductivity between the bio layer and the bottom
double precision, save :: ki_5_bio(2)=2.0;                    ! Initial sea ice thermal conductivity between Vbr=0.05 and the bio layer
double precision, save :: ki_ice_bott(2)=2.0;                 ! Initial sea ice thermal conductivity at the bottom

double precision, save :: ro_br(2)=0.0;
double precision, save :: Vbr_bott(2)=0.0;
double precision, save :: Sbr_bott(2)=0.0;
double precision, save :: ro_br_i(2)= 0.0;
double precision, save :: ro_br_5(2)=0.0;
double precision, save :: ro_br_bio(2)=0.0;
double precision, save :: ro_br_bott(2)=0.0;

double precision, save :: T0_star(2)=284.00D+0;
double precision, save :: Ti_star(2)=Tfreez;

!===========================================================================
end module ice_variables
!===========================================================================

!===========================================================================
module ice_global
!===========================================================================

use ice_params

double precision, save :: h_mix=h_mix_const;
double precision, save :: Tfr=Tfreez;         ! [K]  freezing temp of sea ice

logical, save :: bice_debug = .false.

!===========================================================================
end module ice_global
!===========================================================================

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------

subroutine ice_new2old

use ice_variables

implicit none

 Tmix(1) =  Tmix(2)
 T0(1) =  T0(2)
 Ts(1) =  Ts(2)
 Tsi(1) =  Tsi(2)
 Ti(1) =  Ti(2)
  
 hi(1) =  hi(2)
 hs(1) =  hs(2)
 hmi(1) =  hmi(2)
 hs_prec_bucket(1) =  hs_prec_bucket(2)
  
 snow_fr(1) =  snow_fr(2)
 sea_water_fr(1) =  sea_water_fr(2)
  
 Vbr_ice(1) =  Vbr_ice(2)
 Sbr_ice(1) =  Sbr_ice(2)
 Sice(1) =  Sice(2)
 Sice_bott(1) =  Sice_bott(2)
 Ssnowice(1) =  Ssnowice(2)
 Vbr_i(1) =  Vbr_i(2)
 Sbr_i(1) =  Sbr_i(2)
 Si(1) =  Si(2)
  
 Sice_5(1) =  Sice_5(2)
 Sbr_5(1) =  Sbr_5(2)
 Tice_5(1) =  Tice_5(2)
 hi_5(1) =  hi_5(2)
  
 Sice_bio(1) =  Sice_bio(2)
 Sbr_bio(1) =  Sbr_bio(2)
 hi_bio(1) =  hi_bio(2)
 Vbr_bio(1) =  Vbr_bio(2)
  
 R(1) =  R(2)
  
 ro_sice_bott(1) =  ro_sice_bott(2)
 ro_sice_surf(1) =  ro_sice_surf(2)
 ro_sice_bulk(1) =  ro_sice_bulk(2)
  
 ki_5_bott(1) =  ki_5_bott(2)
 ki_ice_5(1) =  ki_ice_5(2)
 ki_ice_bio(1) =  ki_ice_bio(2)
 ki_bio_bott(1) =  ki_bio_bott(2)
 ki_5_bio(1) =  ki_5_bio(2)
 ki_ice_bott(1) =  ki_ice_bott(2)
  
 ro_br(1)=ro_br(2)
 Vbr_bott(1)=Vbr_bott(2)
 Sbr_bott(1)=Sbr_bott(2)
 ro_br_i(1)=ro_br_i(2)
 ro_br_5(1)=ro_br_5(2)
 ro_br_bio(1)=ro_br_bio(2)
 ro_br_bott(1)=ro_br_bott(2)

 T0_star(1)=T0_star(2)
 Ti_star(1)=Ti_star(2)

end subroutine

!------------------------------------------------------------------

subroutine ice_copy_in(vars)

use ice_variables

implicit none

double precision :: vars(nvars)

 Tmix(1) =  vars(1)
 T0(1) =  vars(2)
 Ts(1) =  vars(3)
 Tsi(1) =  vars(4)
 Ti(1) =  vars(5)
  
 hi(1) =  vars(6)
 hs(1) =  vars(7)
 hmi(1) =  vars(8)
 hs_prec_bucket(1) =  vars(9)
  
 snow_fr(1) =  vars(10)
 sea_water_fr(1) =  vars(11)
  
 Vbr_ice(1) =  vars(12)
 Sbr_ice(1) =  vars(13)
 Sice(1) =  vars(14)
 Sice_bott(1) =  vars(15)
 Ssnowice(1) =  vars(16)
 Vbr_i(1) =  vars(17)
 Sbr_i(1) =  vars(18)
 Si(1) =  vars(19)
  
 Sice_5(1) =  vars(20)
 Sbr_5(1) =  vars(21)
 Tice_5(1) =  vars(22)
 hi_5(1) =  vars(23)
  
 Sice_bio(1) =  vars(24)
 Sbr_bio(1) =  vars(25)
 hi_bio(1) =  vars(26)
 Vbr_bio(1) =  vars(27)
  
 R(1) =  vars(28)
  
 ro_sice_bott(1) =  vars(29)
 ro_sice_surf(1) =  vars(30)
 ro_sice_bulk(1) =  vars(31)
  
 ki_5_bott(1) =  vars(32)
 ki_ice_5(1) =  vars(33)
 ki_ice_bio(1) =  vars(34)
 ki_bio_bott(1) =  vars(35)
 ki_5_bio(1) =  vars(36)
 ki_ice_bott(1) =  vars(37)
  
 ro_br(1)=vars(38)
 Vbr_bott(1)=vars(39)
 Sbr_bott(1)=vars(40)
 ro_br_i(1)=vars(41)
 ro_br_5(1)=vars(42)
 ro_br_bio(1)=vars(43)
 ro_br_bott(1)=vars(44)

 T0_star(1)=vars(45)
 Ti_star(1)=vars(46)

end subroutine

!------------------------------------------------------------------

subroutine ice_copy_out(vars)

use ice_variables

implicit none

double precision :: vars(nvars)

 vars(1) =  Tmix(2)
 vars(2) =  T0(2)
 vars(3) =  Ts(2)
 vars(4) =  Tsi(2)
 vars(5) =  Ti(2)
  
 vars(6) =  hi(2)
 vars(7) =  hs(2)
 vars(8) =  hmi(2)
 vars(9) =  hs_prec_bucket(2)
  
 vars(10) =  snow_fr(2)
 vars(11) =  sea_water_fr(2)
  
 vars(12) =  Vbr_ice(2)
 vars(13) =  Sbr_ice(2)
 vars(14) =  Sice(2)
 vars(15) =  Sice_bott(2)
 vars(16) =  Ssnowice(2)
 vars(17) =  Vbr_i(2)
 vars(18) =  Sbr_i(2)
 vars(19) =  Si(2)
  
 vars(20) =  Sice_5(2)
 vars(21) =  Sbr_5(2)
 vars(22) =  Tice_5(2)
 vars(23) =  hi_5(2)
  
 vars(24) =  Sice_bio(2)
 vars(25) =  Sbr_bio(2)
 vars(26) =  hi_bio(2)
 vars(27) =  Vbr_bio(2)
  
 vars(28) =  R(2)
  
 vars(29) =  ro_sice_bott(2)
 vars(30) =  ro_sice_surf(2)
 vars(31) =  ro_sice_bulk(2)
  
 vars(32) =  ki_5_bott(2)
 vars(33) =  ki_ice_5(2)
 vars(34) =  ki_ice_bio(2)
 vars(35) =  ki_bio_bott(2)
 vars(36) =  ki_5_bio(2)
 vars(37) =  ki_ice_bott(2)
  
 vars(38) =  ro_br(2);
 vars(39) =  Vbr_bott(2);
 vars(40) =  Sbr_bott(2);
 vars(41) =  ro_br_i(2);
 vars(42) =  ro_br_5(2);
 vars(43) =  ro_br_bio(2);
 vars(44) =  ro_br_bott(2);

 vars(45) =  T0_star(2);
 vars(46) =  Ti_star(2);

end subroutine

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------

subroutine ice_check_nvars(nvshy)

use ice_variables

implicit none

integer nvshy

if( nvshy /= nvars ) then
  write(6,*) 'nvshy = ',nvshy
  write(6,*) 'nvars = ',nvars
  stop 'error stop ice_check_vars: nvshy /= nvars'
end if

end subroutine

!------------------------------------------------------------------

subroutine ice_init_vars(vars)

use ice_variables

implicit none

double precision vars(nvars)

call ice_copy_out(vars)

end subroutine

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------

subroutine ice_set_hmix(hm)

use ice_variables
use ice_global

implicit none

real hm

h_mix = hm

end

!------------------------------------------------------------------

subroutine set_ice_debug(bdebug)

 use ice_global

 implicit none

 logical bdebug

 bice_debug = bdebug

end

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------

subroutine ice_output_debug(i)

 use ice_variables

 implicit none

 integer i
 integer, save :: iunit = 130

 write(iunit,*) i

 write(iunit,*) Tmix(1) ,  Tmix(2)
 write(iunit,*) T0(1) ,  T0(2)
 write(iunit,*) Ts(1) ,  Ts(2)
 write(iunit,*) Tsi(1) ,  Tsi(2)
 write(iunit,*) Ti(1) ,  Ti(2)

 write(iunit,*) hi(1) ,  hi(2)
 write(iunit,*) hs(1) ,  hs(2)
 write(iunit,*) hmi(1) ,  hmi(2)
 !hmi_new(1) ,  hmi_new(2)
 write(iunit,*) hs_prec_bucket(1) ,  hs_prec_bucket(2)
 !hs_prec(1) ,  hs_prec(2)
 !hs_tot(1) ,  hs_tot(2)
 !delta_sublim(1) ,  delta_sublim(2)
 !delta_bucket(1) ,  delta_bucket(2)
  
 write(iunit,*) snow_fr(1) ,  snow_fr(2)
 write(iunit,*) sea_water_fr(1) ,  sea_water_fr(2)
  
 !ks(1) ,  ks(2)
 !ks_av(1) ,  ks_av(2)
 !ros(1) ,  ros(2)
 !ros_new(1) ,  ros_new(2)
  
 write(iunit,*) Vbr_ice(1) ,  Vbr_ice(2)
 write(iunit,*) Sbr_ice(1) ,  Sbr_ice(2)
 write(iunit,*) Sice(1) ,  Sice(2)
 !Tice(1) ,  Tice(2)
 write(iunit,*) Sice_bott(1) ,  Sice_bott(2)
 write(iunit,*) Ssnowice(1) ,  Ssnowice(2)
 write(iunit,*) Vbr_i(1) ,  Vbr_i(2)
 write(iunit,*) Sbr_i(1) ,  Sbr_i(2)
 write(iunit,*) Si(1) ,  Si(2)
  
 write(iunit,*) Sice_5(1) ,  Sice_5(2)
 write(iunit,*) Sbr_5(1) ,  Sbr_5(2)
 write(iunit,*) Tice_5(1) ,  Tice_5(2)
 write(iunit,*) hi_5(1) ,  hi_5(2)
  
 write(iunit,*) Sice_bio(1) ,  Sice_bio(2)
 !Tice_bio(1) ,  Tice_bio(2)
 write(iunit,*) Sbr_bio(1) ,  Sbr_bio(2)
 write(iunit,*) hi_bio(1) ,  hi_bio(2)
 write(iunit,*) Vbr_bio(1) ,  Vbr_bio(2)
 !ISI_bio(1) ,  ISI_bio(2)
  
 write(iunit,*) R(1) ,  R(2)
  
 write(iunit,*) ro_sice_bott(1) ,  ro_sice_bott(2)
 write(iunit,*) ro_sice_surf(1) ,  ro_sice_surf(2)
  
 write(iunit,*) ro_sice_bulk(1) ,  ro_sice_bulk(2)
  
 !ks_snow(1) ,  ks_snow(2)
 !kmi_av(1) ,  kmi_av(2)
 !ksi_10_av(1) ,  ksi_10_av(2)
 !ksi_av(1) ,  ksi_av(2)
 !alpha(1) ,  alpha(2)
  
 write(iunit,*) ki_5_bott(1) ,  ki_5_bott(2)
 write(iunit,*) ki_ice_5(1) ,  ki_ice_5(2)
 write(iunit,*) ki_ice_bio(1) ,  ki_ice_bio(2)
 write(iunit,*) ki_bio_bott(1) ,  ki_bio_bott(2)
 write(iunit,*) ki_5_bio(1) ,  ki_5_bio(2)
 write(iunit,*) ki_ice_bott(1) ,  ki_ice_bott(2)
  
 !Sw(1) ,  Sw(2)

 write(iunit,*) ro_br(1),ro_br(2);
 write(iunit,*) Vbr_bott(1),Vbr_bott(2);
 write(iunit,*) Sbr_bott(1),Sbr_bott(2);
 write(iunit,*) ro_br_i(1),ro_br_i(2);
 write(iunit,*) ro_br_5(1),ro_br_5(2);
 write(iunit,*) ro_br_bio(1),ro_br_bio(2);
 write(iunit,*) ro_br_bott(1),ro_br_bott(2);

 !T0_star(1),T0(2);
 !Ts_star(1),Ts(2);
 !Tsi_star(1),Tsi(2);
 !Ti_star(1),Ti(2);

end subroutine

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------

subroutine ice_run(CL,Fsd_cloud,P_rate,qa,qs,Ta,Ua,Sw,deltat &
				,i,vars,iceth)

!               Cl - total cloud cover fraction [percentage]
!               Fsd_cloud - downward solar radiation flux [W m-2]
!               P_rate - precipitation rate [mm s-1]
!               qa - specific humidity of the air [non-dim]
!               qs - specific humidity of the surface [non-dim]
!               Ta - air temperature at 10 m height [K]
!               Ua - wind speed at 10 m height [m s-1]
!               Sw - seawater salinity
!               deltat - time step [s]
!               i - ice model iteration number
!               vars - variables
!               iceth - ice thickness [m]

use ice_params
use ice_variables
use ice_global

implicit none

double precision CL,Fsd_cloud,P_rate,qa,qs,Ta,Ua
double precision Sw,deltat
integer i
double precision vars(nvars)
double precision iceth

double precision Flu_it,Fse_it,Fl_it,F_it,la,sen,lat,cons,Ksnow_ice,coni,conmi,cpi_it
double precision mui,Fld,Fs,Fla,Kmi,ros_new,ros,ks_av,deltahi_bott
double precision deltahi_melt_surf,deltahmi_melt_surf,deltahs_melt_surf,F,F1,F2,Fbott
double precision Fl,Flu,Fse,Fsurf,Fw,hs_prec_bucket_in,hs_prec_bucket_out,I0,IM
double precision ISI,ISI_10,ISI_layer,k0i,k0i_5_bio,k0i_5_bott,k0i_ice_5,k0i_ice_bio
integer j,nmax
double precision k0i_bio_bott,k0i_ice_bott,kb,kb_5_bio,kb_5_bott,kb_bio_bott,kb_ice_5
double precision kb_ice_bio,keff,kice_surf,Sice_mix,Ssnowice_new,kb_ice_bott
double precision kice_bott,kinterm_ice,kksnow_ice,ksnow_interm,mumi,mus,Q,Qi_bott
double precision Qi_surf,Qmi,roi_av,roi_surf,ros_av,Sbr_5_star,Sbr_bio_star
double precision Sbr_bott_star,Sbr_i_star,Sbr_ice_star,Sice_5_mix,Sice_av,Sice_bio_mix
double precision T0_iter,Ti_iter,Tice_av,vi,Tice,Tice_bio,Tmix_star,Ts_iter,ks_snow
double precision Tsi_iter,Vbr_5,hmi_new,hs_prec,delta_sublim,delta_bucket,kks,ISI_bio
double precision kmi_av,ksi_10_av,ksi_av,alpha,Ks,kkice_surf
double precision Ti_star_1,Ts_star_1,Tsi_star_1

double precision, parameter :: eps=0.01D+0
double precision, parameter :: eps2=0.05D+0
double precision, parameter :: keff_p=0.12D+0
double precision, parameter :: keff_p2=0.88D+0
double precision, parameter :: keff_p3=-4.2*1.D4

double precision :: Si_mix
integer :: nloop

call ice_copy_in(vars)

!   FREEZING TEMPERATURE BASED ON SALINITY
call freezing_temperature(Sw,Tfr)

!   SNOW DENSITY
call snow_density(T0(1),ros,ros_new)

!   SNOW PRECIPITATION
call snow_precipitation(hs_prec,Ta,hi(1),P_rate,deltat)

!   PRECIPITATION IN BUCKET MODEL
call precipitation_in_bucket(delta_bucket,hs_prec_bucket(2),ros_av, &
               hs_prec_bucket(1),hs_prec,hs(1),ros,ros_new)
               
!	THERMAL CONDUCTIVITY OF NEW/OLD SNOW
call snow_thermal_conductivity(kks,ks_av,ros_new,ros,hs(1),hs_prec_bucket(2))

!   ALBEDO COMPUTATION	
! (after Flato and Brown, 1996 for landfast sea ice of the Arctic)
call albedo(alpha,T0(1),hi(1),hs(1),hs_prec_bucket(2),hmi(1))

!   EXTICTION COEFFICIENT COMPUTATION
call extinction_coeff(ks_snow,kmi_av,ksi_10_av,ksi_av,T0(1),Cl)

! 	HEAT FLUXES NOT DEPENDENT ON T0
call heat_fluxes(Fs,Fla,Fld,Fw,alpha,Fsd_cloud,qa,Ua,qs,hi(1),Ta,Cl, &
    Tmix(1))

!	SURFACE TEMPERATURE ITERATION AND FLUXES 

!T0_star(1)=T0_old;
!Ts_star(1)=Ts_old;
!Tsi_star(1)=Tsi_old;
Ti_star(1)=Ti(1);

T0_star(1)=T0(1);
Ts_star_1=Ts(1);
Tsi_star_1=Tsi(1);
Ti_star_1=Ti(1);

!   thermal conductivity at the upper boundary of sea ice
kkice_surf= k0 + beta2*Si(1)/(Ti(1)-Tfrs);
!   thermal conductivity at the lower boundary of sea ice     
kice_bott= k0 + beta2*Sice_bott(1)/(Tfr-Tfrs);
!   thermalconductivity at the snow/interm layer interface 
ksnow_interm=(kkmi*ks_av*(hs(1)+hmi(1)))/(hs(1)*kkmi + hmi(1)*ks_av);
!   thermalconductivity at the interm layer/sea ice interface 
kinterm_ice=(kkice_surf*kkmi*(hmi(1)+hi(1)))/(hmi(1)*kkice_surf+hi(1) &
    *kkmi);
!   thermal conductivity at the snow/sea ice layer interface 
kksnow_ice=(kkice_surf*ks_av*(hs(1)+hi(1)))/(hs(1)*kkice_surf + hi(1) &
    *ks_av);

Kice_bott=2*kice_bott/hi(1);
Kice_surf=2*kkice_surf/hi(1);

Kmi=2*kkmi/(hmi(1));
Ks=2*ks_av/hs(1);
Ksnow_interm=2*ksnow_interm/(hmi(1)+hs(1));
Kinterm_ice=2*kinterm_ice/(hmi(1)+hi(1));
Ksnow_ice=2*kksnow_ice/(hi(1)+hs(1));

mus=deltat/(c0*(ros_new*hs_prec_bucket(1)+ros*hs(1)));
mumi=deltat/(romi*c0*hmi(1));

!!!!!!!!!!!!!!!!!!!!!!!!!!! ITERATION START !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j=2,100   !::RASA: why to 100??
    nloop = j
    
    ! ::RASA: ITERATION HEAT FLUXES
    call iteration_heat_fluxes(Flu_it,Fse_it,Fl_it,F_it,la,sen,lat,cons, &
        coni,conmi,cpi_it,mui,T0_star(1),Ua,Ta,Fld,Fs,Fla,qs,ks_av, &
        hs(1),kkice_surf,hi(1),hmi(1),Si(1),Ti_star(1),Ti_star_1, &
        deltat,ro_sice_bulk(1))
    
    !    TEMPERATURE
    if ( (hs(1)+hs_prec_bucket(2))>hs_min ) then
         
        !   SEA ICE MODEL ON/SLAB OCEAN OFF
        if ( hmi(1)>hmi_min ) then

            !   CASE 1: 3 LAYERS (SNOW + INTERM LAYER + SEA ICE)
            call temperature_case1_snow_meteoric_ice(T0_iter,Ts_iter, &
                Tsi_iter,Ti_iter,Tmix_star,T0_star(2),Ti_star(2),Fs,ks_snow, &
                hs(1),hs_prec_bucket(2),kmi_av,hmi(1),hi(1),ksi_10_av, &
                ksi_av,Ks,la,sen,lat,cons,mus,Ksnow_interm,mumi, &
                Kinterm_ice,Kice_bott,T0_star(1),Ts_star_1,Tsi_star_1, &
                Ti_star_1,F_it,mui,Ts(1),Tsi(1),Ti(1))

        else if ( hmi(1)<=hmi_min ) then

            !   CASE 2: 2 LAYERS (SNOW + SEA ICE)
            call temperature_case2_snow_ice(T0_iter,Ts_iter,Tsi_iter, &
                Ti_iter,Tmix_star,T0_star(2),Ti_star(2),Fs,ks_snow,hs(1), &
                hs_prec_bucket(2),hi(1),ksi_10_av,ksi_av,Ks,la,sen,lat, &
                cons,mus,Ksnow_ice,Kice_bott,F_it,mui,T0_star(1), &
                Ts_star_1,Ti_star_1,Ts(1),Ti(1))

        end if

    else if ( (hs(1)+hs_prec_bucket(2))<=hs_min ) then

        if ( hmi(1)>hmi_min ) then

            !   CASE 3: 2 LAYERS (INTERM LAYER + SEA ICE)
            call temperature_case3_meteoric_ice(T0_iter,Ts_iter,Tsi_iter, &
                Ti_iter,Tmix_star,T0_star(2),Ti_star(2),Fs,kmi_av,hmi(1), &
                hi(1),ksi_10_av,ksi_av,la,sen,lat,conmi,mumi, &
                Kinterm_ice,Kice_bott,mui,T0_star(1),Tsi_star_1,F_it, &
                Tsi(1),Ti_star_1,Ti(1),Kmi)

        else if ( hmi(1)<=hmi_min ) then
        
            if ( hi(1)>hi_min ) then
                
                !   CASE 4: ONLY SEA ICE
                call temperature_case4_ice(T0_iter,Ts_iter,Tsi_iter, &
                    Ti_iter,Tmix_star,T0_star(2),Ti_star(2),hi(1),Fs, &
                    ksi_10_av,ksi_av,Kice_surf,la,sen,lat,coni,mui, &
                    Kice_bott,T0_star(1),Ti_star_1,Ti(1),F_it)
                
            else if ( hi(1)<=hi_min ) then
                
                !   SEA ICE MODEL OFF/SLAB OCEAN ON 
                !   CASE 5: 0 LAYER OR FRAZIL ICE FORMATION
                call temperature_case5_frazil(Tmix_star,T0_star(2),Ti_star(2), &
                    T0_iter,Ts_iter,Tsi_iter,Ti_iter,Fs,Tmix(1),F_it, &
                    R(1),deltat)
                
            end if
        else
            print*,"ERROR ITERATION"
        end if
    end if
    
    if ( abs(T0_star(2)-T0_star(1))< eps ) then
        EXIT    !   CONVERENCE CRITERION: 0.01 K
    end if
    
    T0_star(1) = T0_star(2)
    Ti_star(1) = Ti_star(2)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END ITERATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if( nloop > 20 ) then
    print*,"WARNING: ITERATION >20"
end if
if( nloop >= 100 ) then    !   MAX NUMBER OF ITERATIONS FOR CONVERGENCE
    print*,"NOT CONVERGENT"
    STOP
end if

!if ( max(size(T0_star))==100 ) then
!    STOP   !   STOP THE SEA ICE MODEL IF NOT CONVERGENT
!end if

T0(2)=T0_iter;      !   NEW SURFACE TEMPERATURE
Ts(2)=Ts_iter;      !   NEW SNOW TEMPERATURE
Tsi(2)=Tsi_iter;    !   NEW SNOW ICE/SUPERIMPOSED ICE TEMPERATURE
Ti(2)=Ti_iter;      !   NEW SEA ICE TEMPERATURE
Tmix(2)=Tmix_star;  !   NEW MIX LAYER TEMPERATURE

! ::RASA: ! SEA ICE TEMPERATURE
call temperature_of_ice(Tice,hmi(1),kkice_surf,Ti(2),hi(1),Tsi(2),hs(1), &
    hs_prec_bucket(2),ks_av,Ts(2),T0(2))

!   RECOMPUTE SURFACE FLUXES

call surface_fluxes_recompute(Flu,Fl,Fse,F,R(2),Qs,Qmi,Qi_surf,Qi_bott, &
    T0(2),Fld,Ua,Ta,Fs,Fla,ros,Ts(2),Tsi(2),ro_sice_surf(1),Ti(2),Si(1), &
    ro_sice_bott(1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!! SEA ICE GROWTH/DECAY !!!!!!!!!!!!!!!!!!!!!!!
! GROWTH SEASON

if ( T0(2)<Tfr ) then

    if ( hmi(1)>hmi_min ) then

        if ( (hs(1)+hs_prec_bucket(2))>hs_min ) then

            !   CASE 1: SNOW + INTERM LAYER + SEA ICE
            call growth_case1_snow_meteoric_ice(Fsurf,Fbott,I0,IM,ISI_10, &
                ISI,delta_sublim,deltahs_melt_surf,deltahi_melt_surf, &
                deltahmi_melt_surf,deltahi_bott,hs(2),hmi_new,hmi(2),hi(2), &
                hs_prec_bucket_out,hs_prec_bucket_in,Ks,Ts(2),T0(2),Kice_bott, &
                Ti(2),Fs,ks_snow,hs(1),kmi_av,hmi(1),hi(1),ksi_10_av, &
                ksi_av,Fla,deltat,ros_new,Qs,Fw,Qi_bott,ro_sice_surf(1), &
                ros_av,delta_bucket,ros)
                
        else if ( (hs(1)+hs_prec_bucket(2))<=hs_min ) then

            !   CASE 2: INTERM LAYER + SEA ICE
            call growth_case2_meteoric_ice(Fsurf,Fbott,I0,IM,ISI_10,ISI, &
                delta_sublim,deltahs_melt_surf,deltahi_melt_surf, &
                deltahmi_melt_surf,deltahi_bott,hs(2),hmi_new,hmi(2),hi(2),Tsi(2),T0(2), &
                Kice_bott,Ti(2),Fs,hs(1),kmi_av,hmi(1),hi(1),ksi_10_av, &
                ksi_av,Fla,deltat,Qmi,Fw,Qi_bott,delta_bucket,Kmi)
            
        end if

    else if ( hmi(1)<=hmi_min ) then

        if ( (hs(1)+hs_prec_bucket(2))>hs_min ) then
            
            !   CASE 3: SNOW + SEA ICE
            call growth_case3_snow_ice(Fsurf,Fbott,I0,IM,ISI_10,ISI, &
                delta_sublim,deltahs_melt_surf,deltahi_melt_surf, &
                deltahmi_melt_surf,deltahi_bott,hs(2),hmi_new,hmi(2),hi(2), &
                hs_prec_bucket_out,hs_prec_bucket_in,Ks,Ts(2),T0(2),Kice_bott, &
                Ti(2),Fs,ks_snow,hs(1),hmi(1),hi(1),ksi_10_av,ksi_av,Fla, &
                deltat,ros_new,Qs,Fw,Qi_bott,ro_sice_surf(1),ros_av, &
                delta_bucket,ros)

        else if ( (hs(1) + hs_prec_bucket(2))<=hs_min ) then
            
            !   CASE 4: ONLY SEA ICE
            call growth_case4_ice(Fsurf,Fbott,I0,IM,ISI_10,ISI, &
                delta_sublim,deltahs_melt_surf,deltahi_melt_surf, &
                deltahmi_melt_surf,deltahi_bott,hs(2),hmi_new,hmi(2),hi(2),R(2),T0(2), &
                Kice_bott,Ti(2),Fs,hs(1),hmi(1),hi(1),ksi_10_av,ksi_av, &
                Fla,deltat,Fw,Qi_bott,delta_bucket,F,Kice_surf, &
                ro_sice_surf(1),Qi_surf,ro_sice_bott(1),Si(1))
            
        end if
    else
        print*,"ERROR GROWTH SEASON"
    end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! HALODYNAMIC SUBMODEL FOR THE GROWTH SEASON
    if ( hi(2)>ZERO ) then
        
        !double precision vi, keff
        vi=deltahi_bott/deltat;
        ! for Baltic (Gransgok)   ::RASA: with this there is no difference
        ! in results
        !keff(i)=0.113/(0.113 + 0.887*expm(-2.66*10^4*vi(i))); 
        keff = keff_p/(keff_p+keff_p2*exp(keff_p3*vi)); ! with 4.2*10^4 in s/cm
            ! for normal salinity sea water >30, it is usuallly used the 
            !parameterization by Nawako and Shina, 1983
        Sice_bott(2)=keff*Sw;
        Sbr_bott(2)=-(Tfr-Tfrs_appr)/mu;
        Vbr_bott(2)=Sice_bott(2)/Sbr_bott(2);
        Q=0.0;
        
        !   SNOW ICE LAYER
        call growth_halodynamic_snowice(snow_fr(2),sea_water_fr(2),Ssnowice_new,&
            Ssnowice(2),hmi_new,ro_sice_surf(1),ros_av,Sw,Ssnowice(1), &
            hmi(1),hmi(2),snow_fr(1),sea_water_fr(1))

        if ( Vbr_ice(1)>eps2) then       !   DESALINATION
            
            call growth_halodynamic_desalination(Sice_mix,Sice(2),Sbr_ice(2), &
                Vbr_ice(2),Sice_5(2),Sbr_5(2),Tice_5(2),Tice_bio,hi(2),hi(1),Sice_bott(2), &
                Sice(1),Vbr_ice(1),Tice,deltat,ki_bio_bott(1), &
                ki_ice_bio(1),ki_ice_bott(1))
                
        else if ( Vbr_ice(1)<=eps2) then    !   NO DESALINATION
            
            call growth_halodynamic_no_desalination(Sice_mix,Sice(2),Sbr_ice(2),&
                Vbr_ice(2),Sice_5_mix,Sice_5(2),Sbr_5(2),Tice_5(2),Tice_bio,hi(2),hi(1),&
                Sice_bott(2),Sice(1),Tice,deltat,ki_bio_bott(1),hi_5(1), &
                Sice_5(1),Tice_5(1),ki_ice_5(1))
                
        end if
        
        hi_5(2)=(hi(2)*ki_5_bott(1)*(Tice_5(2)-(Tfr-Tfrs_appr)))/((ki_ice_5(1) &
            *(Tice-Tice_5(2)))+(ki_5_bott(1)*(Tice_5(2)-(Tfr-Tfrs_appr))));
        hi_bio(2)=hi_5(2);           
        Sice_bio_mix=((hi(2)-hi(1))*Sice_bott(1) + hi_bio(1)*Sice_bio(1))/(hi(2)&
            -hi(1)+hi_bio(1));!DESALINATION IN THE ACTIVE BIOLOCIAL SYST.
        Sice_bio(2)=Sice_bio_mix+(rosa*(1-viola*Vbr_bio(1))*(Tfr-Tfrs &
            -Tice_bio))*deltat;
        Sbr_bio(2)=-Tice_bio/mu;
        Vbr_bio(2)=Sice_bio(2)/Sbr_bio(2);
        
        Si_mix=((hi(2)-hi(1))*Sice_bott(2) + (hi(1)/2)*Si(1))/(hi(2)-hi(1)+hi(1)/2);
        if ( Vbr_i(1)<=eps2 ) then
            Si(2)=Si_mix;               
        else
            Si(2)=Si_mix + (rosa*(1-viola*Vbr_i(1))*(Tfrs_appr-Ti(2)))*deltat;
        end if
        Sbr_i(2)=-(Ti(2)-Tfrs_appr)/mu;
        Vbr_i(2)=Si(2)/Sbr_i(2);
        
    else if ( hi(2)<=ZERO ) then
    
        ! NO SALINITY COMPUTATION
        call growth_halodynamic_no_salinity(snow_fr(2),sea_water_fr(2), &
            Ssnowice_new,Ssnowice(2),Q,Sice(2),Sice_mix,Sice_5_mix, &
            Sice_bio_mix,Sbr_ice(2),Vbr_ice(2),Sbr_5(2),Sice_5(2),Tice_5(2),hi_5(2),Si(2), &
            Vbr_i(2),Sbr_i(2),Sbr_bio(2),Sice_bio(2),Tice_bio,Vbr_bio(2),hi_bio(2),ISI_bio, &
            Sice_bott(2),Sbr_bott(2),Vbr_bott(2),Sw)
        
    end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!     END OF THE GROWING SEASON ---- START THE MELTING SEASON     !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !   MELT SEASON

    !   THERMODYNAMIC AND HALODYNAMIC SUBMODEL

else if ( T0(2)>=Tfr ) then

    keff=1.0; 
    vi=0.0;        
    hmi_new=0.0;
    snow_fr(2)=0.0;
    sea_water_fr(2)=0.0;
    Ssnowice_new=0.0;
    
    if ( (hs(1)+ hs_prec_bucket(2))>hs_min ) then
        
        !    CASE 1: MELTING SNOW
        call melt_case1_snow(T0(2),I0,Fsurf,Fbott,IM,hmi(2),ISI_10,ISI, &
            deltahs_melt_surf,delta_sublim,deltahi_melt_surf, &
            deltahmi_melt_surf,deltahi_bott,hs(2),hi(2),hs_prec_bucket_out, &
            Fs,ks_snow,hs(1),hs_prec_bucket_in,Ks,Ts(2),Kice_bott,Ti(2), &
            hmi(1),kmi_av,hi(1),ksi_10_av,F,ksi_av,deltat,Qs,Fla, &
            ros_av,Fw,Qi_bott,delta_bucket,ros)

        ! SALINITY FOR MELTING SNOW
        
        if ( hmi(2)>hmi_min ) then
            Ssnowice(2)=(-deltahs_melt_surf*0.20*Ssnow+hmi(2) &
                *Ssnowice(1))/(hmi(2)-deltahs_melt_surf);
        else
            Ssnowice(2)=0.0;
        end if 

        if ( Vbr_ice(1)>=eps2 ) then      !   FLUSHING
            
            call salinity_melting_snow_flushing(Q,Sbr_ice_star,Sice(2), &
                Vbr_ice(2),Sbr_ice(2),Sbr_i_star,Si(2),Vbr_i(2),Sbr_i(2),Sbr_5_star, &
                Sice_5(2),Tice_5(2),Sbr_5(2),Sbr_bio_star,Sice_bio(2),Tice_bio, &
                Vbr_bio(2),Sbr_bio(2),hi_bio(2),Sbr_bott_star,Sice_bott(2),Vbr_bott(2), &
                Sbr_bott(2),hi_5(2),Vbr_ice(1),deltahs_melt_surf,deltat, &
                Sbr_ice(1),ro_br(1),Sice(1),Tice,Sbr_i(1),Vbr_i(1), &
                ro_br_i(1),Si(1),Ti(2),Sbr_5(1),ro_br_5(1),Sice_5(1), &
                Sbr_bio(1),Vbr_bio(1),ro_br_bio(1),Sice_bio(1), &
                ki_bio_bott(1),ki_ice_bio(1),ki_ice_bott(1), &
                Sbr_bott(1),Vbr_bott(1),ro_br_bott(1),Sice_bott(1), &
                hi(2),ki_5_bott(1),ki_ice_5(1),ros)

        else if ( Vbr_ice(1)<eps2 ) then   !   NO FLUSHING
            
            call salinity_melting_snow_snowice_no_flushing(Q,Sice(2), &
                Vbr_ice(2),Sbr_ice(2),Si(2),Vbr_i(2),Sbr_i(2), &
                Sice_5(2),Tice_5(2),Sbr_5(2),Sice_bio(2),Tice_bio, &
                Vbr_bio(2),Sbr_bio(2),hi_bio(2),Sice_bott(2), &
                Vbr_bott(2),Sbr_bott(1),hi_5(2),Sice(1),Tice,Si(1), &
                Ti(2),Sice_5(1),Sice_bio(1),ki_bio_bott(1), &
                Sice_bott(1),hi(2),ki_5_bott(1),ki_ice_5(1))

        end if

    else if ( (hs(1)+hs_prec_bucket(2))<=hs_min ) then
        
        if ( hmi(1)>hmi_min ) then
            
            !   CASE 2: MELTING SNOW ICE
            call melt_case2_snowice(hs(2),T0(2),I0,IM,ISI_10,ISI,Fsurf,Fbott, &
                delta_sublim,deltahmi_melt_surf,deltahi_melt_surf, &
                deltahs_melt_surf,deltahi_bott,hmi(2),hi(2),R(2),Fs,kmi_av,hmi(1),&
                hi(1),ksi_10_av,ksi_av,Tsi(2),Kice_bott,Ti(2),Fla,deltat,Qmi,F,&
                Fw,Qi_bott,delta_bucket,Qi_surf,Kmi)
            
            !   SALINITY FOR MELTING SNOW ICE
            
            Ssnowice(2)=Ssnowice(1);
            
            if ( Vbr_ice(1)>=eps2 ) then       !   FLUSHING

                call salinity_melting_snowice_flushing(Q,Sbr_ice_star, &
                    Sice(2),Vbr_ice(2),Sbr_ice(2),Sbr_i_star,Si(2),Vbr_i(2),Sbr_i(2), &
                    Sbr_5_star,Sice_5(2),Tice_5(2),Sbr_5(2),Sbr_bio_star,Sice_bio(2), &
                    Tice_bio,Vbr_bio(2),Sbr_bio(2),hi_bio(2),Sbr_bott_star, &
                    Sice_bott(2),Vbr_bott(2),Sbr_bott(2),hi_5(2),Vbr_ice(1), &
                    deltahmi_melt_surf,deltat,Sbr_ice(1),ro_br(1), &
                    Sice(1),Tice,Sbr_i(1),Vbr_i(1),ro_br_i(1),Si(1), &
                    Ti(2),Sbr_5(1),ro_br_5(1),Sice_5(1),Sbr_bio(1), &
                    Vbr_bio(1),ro_br_bio(1),Sice_bio(1), &
                    ki_bio_bott(1),ki_ice_bio(1),ki_ice_bott(1), &
                    Sbr_bott(1),Vbr_bott(1),ro_br_bott(1),Sice_bott(1))

            else if ( Vbr_ice(1)<eps2 ) then    !   NO FLUSHING
                
                call salinity_melting_snow_snowice_no_flushing(Q,Sice(2), &
                    Vbr_ice(2),Sbr_ice(2),Si(2),Vbr_i(2),Sbr_i(2),Sice_5(2),Tice_5(2),Sbr_5(2), &
                    Sice_bio(2),Tice_bio,Vbr_bio(2),Sbr_bio(2),hi_bio(2),Sice_bott(2), &
                    Vbr_bott(2),Sbr_bott(2),hi_5(2),Sice(1),Tice,Si(1),Ti(2), &
                    Sice_5(1),Sice_bio(1),ki_bio_bott(1),Sice_bott(1),&
                    hi(2),ki_5_bott(1),ki_ice_5(1))

            end if

        else if ( hmi(1)<=hmi_min ) then
            
            if ( hi(1)>ZERO ) then
                
                !   CASE 3: MELTING SEA ICE
                call melt_case3_ice(hs(2),hmi(2),T0(2),ISI_10,ISI,I0,IM,Fsurf, &
                    Fbott,delta_sublim,deltahi_bott,deltahi_melt_surf, &
                    deltahs_melt_surf,deltahmi_melt_surf,hi(2),R(2),hi(1),Fs, &
                    ksi_10_av,ksi_av,Kice_surf,Ti(2),Kice_bott,Fla,deltat, &
                    ro_sice_surf(1),Qi_surf,Fw,Qi_bott,F,delta_bucket)
                
                !   SALINITY FOR MELTING SEA ICE
                Ssnowice(2)=0.0;
                if ( hi(2)>ZERO ) then
                    
                    if ( Vbr_ice(1)>=eps2 ) then     !   FLUSHING
                        
                        call salinity_melting_ice_flushing(Q,Sbr_ice_star,&
                            Sice(2),Vbr_ice(2),Sbr_ice(2),Sbr_i_star,Si(2),Vbr_i(2), &
                            Sbr_i(2),Sbr_5_star,Sice_5(2),Tice_5(2),Vbr_5,Sbr_5(2), &
                            hi_5(2),Sbr_bio_star,Sice_bio(2),Tice_bio,Vbr_bio(2), &
                            Sbr_bio(2),hi_bio(2),Sbr_bott_star,Sice_bott(2), &
                            Vbr_bott(2),Sbr_bott(2),Vbr_ice(1), &
                            ro_sice_surf(1),deltahi_melt_surf,deltat, &
                            Sbr_ice(1),ro_br(1),Sice(1),Tice,Sbr_i(1),&
                            Vbr_i(1),ro_br_i(1),Si(1),Ti(2),Sbr_5(1), &
                            ro_br_5(1),Sice_5(1),ki_5_bott(1), &
                            ki_ice_5(1),hi(2),Sbr_bio(1),Vbr_bio(1), &
                            ro_br_bio(1),Sice_bio(1),ki_bio_bott(1), &
                            ki_ice_bio(1),ki_ice_bott(1),Sbr_bott(1), &
                            Vbr_bott(1),ro_br_bott(1),Sice_bott(1))
                        
                    else if ( Vbr_ice(1)<eps2 ) then   !   NO FLUSHING
                        
                        call salinity_melting_ice_no_flushing(Q,Sice(2), &
                            Vbr_ice(2),Sbr_ice(2),Si(2),Vbr_i(2),Sbr_i(2),Sice_5(2),Tice_5(2), &
                            Sbr_5(2),hi_5(2),Sice_bio(2),Tice_bio,Vbr_bio(2),Sbr_bio(2), &
                            hi_bio(2),Sice_bott(2),Vbr_bott(2),Sbr_bott(2),Sice(1), &
                            Tice,Si(1),Ti(2),Sice_5(1),ki_5_bott(1), &
                            ki_ice_5(1),hi(2),Sice_bio(1),ki_bio_bott(1), &
                            Sice_bott(1),ki_5_bio(1))

                    end if
                    
                else if ( hi(2)<=ZERO ) then
                    
                    !   NO SALINITY COMPUTATION
                    call no_salinity(Q,Sice(2),Sice_bott(2),Sbr_ice(2),Vbr_ice(2),Si(2), &
                        Sbr_i(2),Vbr_i(2),Sbr_5(2),Sice_5(2),Tice_5(2),hi_5(2),Sbr_bio(2), &
                        Sice_bio(2),Tice_bio,Vbr_bio(2),hi_bio(2),Sbr_bott(2),Vbr_bott(2))
                
                end if

            else if ( hi(1)<=ZERO ) then
                
                !   CASE 4: ONLY SEAWATER
                call melt_case4_water(ISI,I0,Fsurf,Fbott,IM,ISI_10, &
                    deltahs_melt_surf,deltahmi_melt_surf, &
                    deltahi_melt_surf,deltahi_bott,delta_bucket, &
                    hs_prec_bucket(2),hs(2),hi(2),hmi(2),R(2),Fs,h_mix)
                !write(6,*) 'Tfr',Tfr
                !   NO SALINITY COMPUTATION
                call no_salinity(Q,Sice(2),Sice_bott(2),Sbr_ice(2),Vbr_ice(2),Si(2), &
                    Sbr_i(2),Vbr_i(2),Sbr_5(2),Sice_5(2),Tice_5(2),hi_5(2),Sbr_bio(2),Sice_bio(2),&
                    Tice_bio,Vbr_bio(2),hi_bio(2),Sbr_bott(2),Vbr_bott(2))
                Ssnowice(2)=0.0;

            end if
        end if
    end if
else
    print*,"ERROR MELTING SEASON"
end if

!   PROPERTIES DEPENDING ON T,S

! ::RASA: densities
call depending_on_TS_densities(Tice_av,Sice_av,F1,F2,roi_surf, &
    ro_sice_surf(2),ro_sice_bott(2),roi_av,ro_sice_bulk(2),ro_br(2),ro_br_i(2),ro_br_5(2), &
    ro_br_bio(2),ro_br_bott(2),Ti(2),Sice_bott(2),Sice(2),Tice,Sbr_ice(2),Sbr_i(2),Sbr_5(2), &
    Sbr_bio(2),Sbr_bott(2))

! ::RASA: thermal conductivities
call depending_on_TS_thermal_conductivity(k0i,kb,k0i_5_bio,kb_5_bio, &
    ki_5_bio(2),k0i_ice_5,kb_ice_5,ki_ice_5(2),k0i_ice_bott,kb_ice_bott, &
    ki_ice_bott(2),k0i_5_bott,kb_5_bott,ki_5_bott(2),k0i_ice_bio,kb_ice_bio, &
    ki_ice_bio(2),k0i_bio_bott,kb_bio_bott,ki_bio_bott(2),Tice_av,Tice_5(2), &
    Tice_bio,Tice,Vbr_ice(2),Vbr_bio(2))

!   COMPUTATION OF MEAN PROPERTIES OF THE BIOLOGICAL ACTIVE PART OF THE SEA
!   ICE SISTEM FOR BFM RUN

call mean_BAL_for_BFM(ISI_layer,ISI_bio,hi(2),I0,ksi_av,hi_5(2),IM,ksi_10_av)

if( i > 0 ) call ice_output_debug(i)

call ice_copy_out(vars)
!call ice_new2old

iceth=hi(2)

end
! INITIATE FUNCTION WITH THIS:
![alpha(i)] = albedo(T0(i-1),Tfr,alpha_mi,alpha_ow,c11,hi(i-1),c12,hs(i-1),
!       hs_prec_bucket(i),hmi(i-1),hmi_min,hi_min,c10);

! ALBEDO COMPUTATION (after Flato and Brown, 
!       1996 for landfast sea ice of the Arctic)
!subroutine albedo(alpha,T0_old,Tfr,alpha_mi,alpha_ow,c11,hi_old,c12,hs_old, &
!    hs_prec_bucket,hmi_old,hmi_min,hi_min,c10)

subroutine albedo(alpha,T0_old,hi_old,hs_old,hs_prec_bucket,hmi_old)
    
use ice_params
use ice_global

implicit none
double precision alpha,T0_old,hi_old,hs_old,hs_prec_bucket,hmi_old
double precision alpha_mi
double precision alpha_i, alpha_s

double precision, parameter :: alpha_i_1 = 0.28D+0
double precision, parameter :: alpha_i_2 = 0.08D+0
double precision, parameter :: eps  = 0.1D+0

alpha_mi=0.55

! ::RASA: albedo of ice(i), snow(s) and meteoric ice(mi)
if( T0_old<Tfr ) then 
    alpha_i=max(alpha_ow, c11*hi_old**alpha_i_1 +alpha_i_2);  
    alpha_s=0.75D+0;
    alpha_mi=0.70D+0;  ! (after Perovich, 1996 for compacted snow)
else
    alpha_i=min(alpha_mi,c12*hi_old**2 +alpha_ow);
    alpha_s=0.65D+0;      
    alpha_mi=0.56D+0;  ! (after Perovich, 1996 for melting white ice)
end if

if( (hs_old+hs_prec_bucket)>eps ) then
    alpha=alpha_s;
else if ( hmi_old>hmi_min ) then
    alpha=alpha_mi;
else if ( hi_old>hi_min ) then
    alpha=min(alpha_s,alpha_i+(hs_old+hs_prec_bucket)*(alpha_s-alpha_i)/c10);
else
    alpha=alpha_ow;
end if

end
! INITIATE FUNCTION WITH THIS:
![Tice_av(i),Sice_av(i),F1(i),F2(i),roi_surf(i),ro_sice_surf(i),
!    ro_sice_bott(i),roi_av(i),ro_sice_bulk(i),ro_br(i),ro_br_i(i),ro_br_5(i),
!    ro_br_bio(i),ro_br_bott(i)]= ...
!    depending_on_TS_densities(Ti(i),Sice_bott(i),Sice(i),Tice(i),Va,
!    Sbr_ice(i),Sbr_i(i),Sbr_5(i),Sbr_bio(i),Sbr_bott(i));

! PROPERTIES DEPENDING ON T,S [densities]
!subroutine depending_on_TS_densities(Tice_av,Sice_av,F1,F2,roi_surf, &
!    ro_sice_surf,ro_sice_bott,roi_av,ro_sice_bulk,ro_br,ro_br_i,ro_br_5, &
!    ro_br_bio,ro_br_bott,Ti,Sice_bott,Sice,Tice,Va,Sbr_ice,Sbr_i,Sbr_5, &
!    Sbr_bio,Sbr_bott)

subroutine depending_on_TS_densities(Tice_av,Sice_av,F1,F2,roi_surf, &
    ro_sice_surf,ro_sice_bott,roi_av,ro_sice_bulk,ro_br,ro_br_i,ro_br_5, &
    ro_br_bio,ro_br_bott,Ti,Sice_bott,Sice,Tice,Sbr_ice,Sbr_i,Sbr_5, &
    Sbr_bio,Sbr_bott)
    
use ice_params

implicit none
double precision Tice_av,Sice_av,F1,F2,roi_surf,ro_sice_surf,ro_sice_bott,roi_av
double precision ro_sice_bulk,ro_br,ro_br_i,ro_br_5,ro_br_bio,ro_br_bott,Ti,Sice_bott
double precision Sice,Tice,Sbr_ice,Sbr_i,Sbr_5,Sbr_bio,Sbr_bott

double precision, parameter :: eps=-2.000D+0
double precision, parameter :: F1_p=-4.1221*1.D-2
double precision, parameter :: F1_p2=1.8407*1.D+1
double precision, parameter :: F1_p3=5.4802*1.D-1
double precision, parameter :: F1_p4=2.1454*1.D-1
double precision, parameter :: F1_p5=-4.723D+0
double precision, parameter :: F1_p6=2.245*1.D+1
double precision, parameter :: F1_p7=6.397*1.D-1
double precision, parameter :: F1_p8=1.074*1.D-2
double precision, parameter :: F2_p=9.0312*1.D-2
double precision, parameter :: F2_p2=1.6111*1.D-2
double precision, parameter :: F2_p3=1.2291*1.D-4
double precision, parameter :: F2_p4=1.3603*1.D-4
double precision, parameter :: F2_p5=8.903*1.D-2
double precision, parameter :: F2_p6=1.763*1.D-2
double precision, parameter :: F2_p7=5.330*1.D-4
double precision, parameter :: F2_p8=8.801*1.D-4
double precision, parameter :: roi_p=0.917D+0
double precision, parameter :: roi_p2=1.403*1.D-4
double precision, parameter :: ro_br_p=8*1.D-4

Tice_av=Ti-Tfrs;               !   =============>>>>>>>>>>>>>>>> TO BFM
Sice_av=(Sice_bott+Sice)/2;      !   =============>>>>>>>>>>>>>>>> TO BFM

if( Tice_av > -eps ) then
    F1= F1_p - F1_p2*Tice_av + F1_p3*Tice_av**2 + F1_p4*Tice_av**3;
    F2= F2_p -F2_p2*Tice_av +F2_p3*Tice_av**2 +F2_p4*Tice_av**3;
else
    F1= F1_p5 -F1_p6*Tice_av - F1_p7*Tice_av**2 -F1_p8*Tice_av**3;
    F2= F2_p5 -F2_p6*Tice_av -F2_p7*Tice_av**2 -F2_p8*Tice_av**3;
end if

! density of pure ice at the surface
roi_surf=roi_p-roi_p2*Tice;
! bulk density of sea ice at the surface as function of T,S;
!(1-Va)*((roi_surf(i).*F1(i))./(F1(i)-roi_surf(i).*Sice(i).*F2(i)));
ro_sice_surf=0.905D+0;
! bulk density of sea ice at the bottm (density of pure ice at the bottom is
!    fixed at 0.900)
!(1-Va)*((roi_bott*F1(i))./(F1(i)-roi_bott*Sice_bott(i).*F2(i)));
ro_sice_bott=0.885D+0;
! density of pure ice of the layer as function of Tice_av in g/cm^3
roi_av=roi_p-roi_p2*Tice_av;
! bulk density of ice as function of Tice,Sice in g/cm^3
ro_sice_bulk=(1-Va)*((roi_av*F1)/(F1-roi_av*Sice_av*F2));

! brines density as function of brines salinity in g/cm^3
ro_br=1+ro_br_p*Sbr_ice;
ro_br_i=1+ro_br_p*Sbr_i;
ro_br_5=1+ro_br_p*Sbr_5;
ro_br_bio=1+ro_br_p*Sbr_bio;
ro_br_bott=1+ro_br_p*Sbr_bott;

end
! INITIATE FUNCTION WITH THIS:
![k0i(i),kb(i),k0i_5_bio(i),kb_5_bio(i),ki_5_bio(i),k0i_ice_5(i),kb_ice_5(i),
!    ki_ice_5(i),k0i_ice_bott(i),kb_ice_bott(i),ki_ice_bott(i),k0i_5_bott(i),
!    kb_5_bott(i),ki_5_bott(i),k0i_ice_bio(i),kb_ice_bio(i),ki_ice_bio(i),
!    k0i_bio_bott(i),kb_bio_bott(i),ki_bio_bott(i)]=...
!    depending_on_TS_thermal_conductivity(Tice_av(i),Tice_5(i),Tice_bio(i),Va,
!    Tice(i),Vbr_ice(i),Tfr,Vbr_bio(i));

! PROPERTIES DEPENDING ON T,S [thermal conductivity]
!subroutine depending_on_TS_thermal_conductivity(k0i,kb,k0i_5_bio,kb_5_bio, &
!    ki_5_bio,k0i_ice_5,kb_ice_5,ki_ice_5,k0i_ice_bott,kb_ice_bott, &
!    ki_ice_bott,k0i_5_bott,kb_5_bott,ki_5_bott,k0i_ice_bio,kb_ice_bio, &
!    ki_ice_bio,k0i_bio_bott,kb_bio_bott,ki_bio_bott,Tice_av,Tice_5,Tice_bio, &
!    Va,Tice,Vbr_ice,Tfr,Vbr_bio)

subroutine depending_on_TS_thermal_conductivity(k0i,kb,k0i_5_bio,kb_5_bio, &
    ki_5_bio,k0i_ice_5,kb_ice_5,ki_ice_5,k0i_ice_bott,kb_ice_bott, &
    ki_ice_bott,k0i_5_bott,kb_5_bott,ki_5_bott,k0i_ice_bio,kb_ice_bio, &
    ki_ice_bio,k0i_bio_bott,kb_bio_bott,ki_bio_bott,Tice_av,Tice_5,Tice_bio, &
    Tice,Vbr_ice,Vbr_bio)

use ice_params
use ice_global

implicit none
double precision k0i,kb,k0i_5_bio,kb_5_bio,ki_5_bio,k0i_ice_5,kb_ice_5,ki_ice_5
double precision k0i_ice_bott,kb_ice_bott,ki_ice_bott,k0i_5_bott,kb_5_bott,ki_5_bott
double precision k0i_ice_bio,kb_ice_bio,ki_ice_bio,k0i_bio_bott,kb_bio_bott,ki_bio_bott
double precision Tice_av,Tice_5,Tice_bio,Tice,Vbr_ice,Vbr_bio

double precision, parameter :: k_p = 418.6D+0
double precision, parameter :: k_p2 = 5.35*1.D-3
double precision, parameter :: k_p3 = 2.568*1.D-5
double precision, parameter :: k_p4 = 1.25*1.D-3
double precision, parameter :: k_p5 = 3.0*1.D-5
double precision, parameter :: k_p6 = 0.05D+0
double precision, parameter :: k_p7 = 0.005D+0

! thermal conductivity of pure ice as function of T
k0i = k_p*(k_p2 - k_p3*Tice_av);
! thermal conductivity of pure brines as function of T
kb = k_p*(k_p4 + k_p5*Tice_av);

k0i_5_bio=k_p*(k_p2 - k_p3*((Tice_5+Tice_bio)/2));
kb_5_bio=k_p*(k_p4 + k_p5*((Tice_5+Tice_bio)/2));
ki_5_bio= (1-Va-k_p6)*k0i_5_bio + k_p7*kb_5_bio;

k0i_ice_5=k_p*(k_p2 - k_p3*((Tice+Tice_5)/2));
kb_ice_5=k_p*(k_p4 + k_p5*((Tice+Tice_5)/2));
ki_ice_5=(1-Va-Vbr_ice)*k0i_ice_5 + Vbr_ice*kb_ice_5;

k0i_ice_bott=k_p*(k_p2 -k_p3*((Tice+(Tfr-Tfrs))/2));
kb_ice_bott=k_p*(k_p4 + k_p5*((Tice+(Tfr-Tfrs))/2));
ki_ice_bott=(1-Va-Vbr_ice)*k0i_ice_bott + Vbr_ice*kb_ice_bott;

k0i_5_bott=k_p*(k_p2 - k_p3*((Tice_5+(Tfr-Tfrs))/2));
kb_5_bott=k_p*(k_p4 + k_p5*((Tice_5+(Tfr-Tfrs))/2));
ki_5_bott=(1-Va-k_p6)*k0i_5_bott + k_p6*kb_5_bott;

k0i_ice_bio=k_p*(k_p2 - k_p3*((Tice+Tice_bio)/2));
kb_ice_bio=k_p*(k_p4 + k_p5*((Tice+Tice_bio)/2));
ki_ice_bio= (1-Va-Vbr_ice)*k0i_ice_bio + Vbr_ice*kb_ice_bio;

k0i_bio_bott=k_p*(k_p2 - k_p3*((Tice_bio+(Tfr-Tfrs))/2));
kb_bio_bott=k_p*(k_p4 + k_p5*((Tice_bio+(Tfr-Tfrs))/2));
ki_bio_bott= (1-Va-Vbr_bio)*k0i_bio_bott + Vbr_bio*kb_bio_bott;

end
! INITIATE FUNCTION WITH THIS:
![ks_snow(i),kmi_av(i),ksi_10_av(i),ksi_av(i)] = extinction_coeff(T0(i-1),
!      Tfr,Cl(i));

! EXTICTION COEFFICIENT COMPUTATION
subroutine extinction_coeff(ks_snow,kmi_av,ksi_10_av,ksi_av,T0_old,Cl)

use ice_params
use ice_global

implicit none
double precision ks_snow,kmi_av,ksi_10_av,ksi_av,T0_old,Cl

double precision, parameter :: kmi_av_p=-3.3D+0
double precision, parameter :: kmi_av_p2=21.05D+0;
double precision, parameter :: kmi_av_p3=-1.9D+0;
double precision, parameter :: kmi_av_p4=11.7D+0;
double precision, parameter :: ksi_10_av_p=-6.6D+0;
double precision, parameter :: ksi_10_av_p2=17.1D+0;
double precision, parameter :: ksi_10_av_p3=-3.8D+0;
double precision, parameter :: ksi_10_av_p4=8.4D+0;
double precision, parameter :: ksi_av_p=-0.1D+0;
double precision, parameter :: ksi_av_p2=1.6D+0;
double precision, parameter :: ksi_av_p3=1.5D+0;

if( T0_old<Tfr ) then
    ks_snow=25.0D+0;
    kmi_av=kmi_av_p*Cl+kmi_av_p2;
    ksi_10_av=ksi_10_av_p*Cl+ksi_10_av_p2;
    ksi_av=ksi_av_p*Cl+ksi_av_p2;
else if ( T0_old>=Tfr ) then
    ks_snow=15.0D+0;
    kmi_av=kmi_av_p3*Cl+kmi_av_p4
    ksi_10_av=ksi_10_av_p3*Cl+ksi_10_av_p4
    ksi_av=ksi_av_p*Cl+ksi_av_p3;
else
    print *,'ERROR EXTICTION COEFFICIENT'
end if

end

! INITIATE FUNCTION WITH THIS:
![Fsurf(i),Fbott(i),I0(i),IM(i),ISI_10(i),ISI(i),delta_sublim(i),
!    deltahs_melt_surf(i),deltahi_melt_surf(i),deltahmi_melt_surf(i),
!    deltahi_bott(i),hs(i),hmi_new(i),hmi(i),hi(i),hs_prec_bucket(i),
!    hmi_new(i)] = growth_case1_snow_meteoric_ice(hs_prec_bucket(i),Ks(i),
!    Ts(i),T0(i),Kice_bott(i),Tfr,Ti(i),Fs(i),ks_snow(i),hs(i-1),kmi_av(i),
!    hmi(i-1),hi(i-1),ksi_10_av(i),ksi_av(i),Fla(i),deltat,ros_new,Lv,Qs(i),
!    Fw(i),Qi_bott(i),row,ro_sice_surf(i-1),ros_av(i),romi,delta_bucket(i),
!    beta,ros(i),si_fract,hs_min,vis_fr);

! CASE 1: SNOW + INTERM LAYER + SEA ICE
!subroutine growth_case1_snow_meteoric_ice(Fsurf,Fbott,I0,IM,ISI_10,ISI, &
!    delta_sublim,deltahs_melt_surf,deltahi_melt_surf,deltahmi_melt_surf, &
!    deltahi_bott,hs,hmi_new,hmi,hi,hs_prec_bucket_out,hs_prec_bucket_in,Ks, &
!    Ts,T0,Kice_bott,Tfr,Ti,Fs,ks_snow,hs_old,kmi_av,hmi_old,hi_old,ksi_10_av, &
!    ksi_av,Fla,deltat,ros_new,Lv,Qs,Fw,Qi_bott,row,ro_sice_surf_old,ros_av, &
!    romi,delta_bucket,beta,ros,si_fract,hs_min,vis_fr)

subroutine growth_case1_snow_meteoric_ice(Fsurf,Fbott,I0,IM,ISI_10,ISI, &
    delta_sublim,deltahs_melt_surf,deltahi_melt_surf,deltahmi_melt_surf, &
    deltahi_bott,hs,hmi_new,hmi,hi,hs_prec_bucket_out,hs_prec_bucket_in,Ks, &
    Ts,T0,Kice_bott,Ti,Fs,ks_snow,hs_old,kmi_av,hmi_old,hi_old,ksi_10_av, &
    ksi_av,Fla,deltat,ros_new,Qs,Fw,Qi_bott,ro_sice_surf_old,ros_av, &
    delta_bucket,ros)

use ice_params
use ice_global

implicit none
double precision Fsurf,Fbott,I0,IM,ISI_10,ISI,delta_bucket,ros
double precision delta_sublim,deltahs_melt_surf,deltahi_melt_surf,deltahmi_melt_surf
double precision deltahi_bott,hs,hmi_new,hmi,hi,hs_prec_bucket_out,hs_prec_bucket_in,Ks
double precision Ts,T0,Kice_bott,Ti,Fs,ks_snow,hs_old,kmi_av,hmi_old,hi_old,ksi_10_av
double precision ksi_av,Fla,deltat,ros_new,Qs,Fw,Qi_bott,ro_sice_surf_old,ros_av

double precision, parameter :: eps=0.1D+0
double precision, parameter :: thsnd=1.D3

Fsurf=Ks*(Ts-T0);               !::RASA: conductive flux at the surface 
Fbott=Kice_bott*(Tfr-Ti);       !::RASA: conductive flux at the bottom 

!::RASA: fraction of solar radiation penetrating the surface
I0=vis_fr*Fs*exp(-ks_snow*(hs_old+hs_prec_bucket_in));
!::RASA: fraction of solar radiation penetrating snow ice/supeimposed ice
IM=I0*exp(-kmi_av*hmi_old);

if( hi_old>=eps ) then
    !::RASA: fraction of solar radiation penetrating sea ice:
    ISI_10=IM*exp(-ksi_10_av*eps);          !in the above first 10cm
    ISI=ISI_10*exp(-ksi_av*(hi_old-eps));   !below the first 10cm
else if( hi_old<eps ) then
    ISI_10=IM*exp(-ksi_10_av*hi_old);
    ISI=ISI_10;
end if

delta_sublim=(Fla*deltat)/(ros_new*Lv-Qs);
deltahs_melt_surf=0.0;
deltahi_melt_surf=0.0;
deltahmi_melt_surf=0.0;
deltahi_bott=deltat*(Fbott-Fw-ISI)/Qi_bott;

!   SNOW ICE FORMATION
if( (hs_old+hs_prec_bucket_in-delta_sublim)>((hi_old*(row-ro_sice_surf_old &
        *thsnd)/ros_av) + (hmi_old*(row-romi)/ros_av)) ) then
    hs=hs_old+delta_bucket - beta*((-row*hi_old + ros*(hs_old-delta_sublim) &
        + thsnd*ro_sice_surf_old*hi_old - row*hmi_old + romi*hmi_old + ros_new &
        *hs_prec_bucket_in)/(-romi + row+ beta*ros_av));
    hmi_new=((-row*hi_old + ros*(hs_old-delta_sublim) + thsnd*ro_sice_surf_old &
        *hi_old - row*hmi_old + romi*hmi_old+ ros_new*hs_prec_bucket_in) &
        /(-romi+ row+ beta*ros_av)*si_fract);
    hmi=hmi_old + ((-row*hi_old + ros*(hs_old-delta_sublim) + thsnd* &
        ro_sice_surf_old*hi_old - row*hmi_old + romi*hmi_old+ ros_new &
        *hs_prec_bucket_in)/(-romi + row+ beta*ros_av)*si_fract);
else
    hs=hs_old+ delta_bucket - delta_sublim;
    hmi_new=0.0;
    hmi=hmi_old;
end if
hi=hi_old+ deltahi_bott;

!::RASA: emptying the bucket?
if( hs<hs_min ) then
    hs=hs_prec_bucket_in;
    hs_prec_bucket_out=0.0;
else  !::RASA: added to return the same value
    hs_prec_bucket_out=hs_prec_bucket_in;
end if

end
! INITIATE FUNCTION WITH THIS:
![Fsurf(i),Fbott(i),I0(i),IM(i),ISI_10(i),ISI(i),delta_sublim(i),
!    deltahs_melt_surf(i),deltahi_melt_surf(i),deltahmi_melt_surf(i),
!    deltahi_bott(i),hs(i),hmi_new(i),hmi(i),hi(i)] = ...
!    growth_case2_meteoric_ice(Tsi(i),T0(i),Kice_bott(i),Tfr,Ti(i),Fs(i),
!    hs(i-1),kmi_av(i),hmi(i-1),hi(i-1),ksi_10_av(i),ksi_av(i),Fla(i),deltat,
!    Qmi(i),Fw(i),Qi_bott(i),romi,delta_bucket(i),Kmi(i),vis_fr,Lv);

! CASE 2: INTERM LAYER + SEA ICE
!subroutine growth_case2_meteoric_ice(Fsurf,Fbott,I0,IM,ISI_10,ISI, &
!    delta_sublim,deltahs_melt_surf,deltahi_melt_surf,deltahmi_melt_surf, &
!    deltahi_bott,hs,hmi_new,hmi,hi,Tsi,T0,Kice_bott,Tfr,Ti,Fs,hs_old,kmi_av, &
!    hmi_old,hi_old,ksi_10_av,ksi_av,Fla,deltat,Qmi,Fw,Qi_bott,romi, &
!    delta_bucket,Kmi,vis_fr,Lv)

subroutine growth_case2_meteoric_ice(Fsurf,Fbott,I0,IM,ISI_10,ISI, &
    delta_sublim,deltahs_melt_surf,deltahi_melt_surf,deltahmi_melt_surf, &
    deltahi_bott,hs,hmi_new,hmi,hi,Tsi,T0,Kice_bott,Ti,Fs,hs_old,kmi_av, &
    hmi_old,hi_old,ksi_10_av,ksi_av,Fla,deltat,Qmi,Fw,Qi_bott,delta_bucket,Kmi)
    
use ice_params
use ice_global

implicit none
double precision Fsurf,Fbott,I0,IM,ISI_10,ISI,delta_sublim,Kmi
double precision deltahs_melt_surf,deltahi_melt_surf,deltahmi_melt_surf,deltahi_bott,hs
double precision hmi_new,hmi,hi,Tsi,T0,Kice_bott,Ti,Fs,hs_old,kmi_av,hmi_old,hi_old
double precision ksi_10_av,ksi_av,Fla,deltat,Qmi,Fw,Qi_bott,delta_bucket

double precision, parameter :: eps=0.1D+0

Fsurf=Kmi*(Tsi-T0);
Fbott=Kice_bott*(Tfr-Ti);
I0=vis_fr*Fs*exp(-kmi_av*hmi_old);
IM=I0;

if( hi_old>=eps ) then
    ISI_10=IM*exp(-ksi_10_av*eps);
    ISI=ISI_10*exp(-ksi_av*(hi_old-eps));
else if( hi_old<eps ) then
    ISI_10=IM*exp(-ksi_10_av*hi_old);
    ISI=ISI_10;
end if

delta_sublim=(Fla*deltat)/(romi*Lv-Qmi);
deltahmi_melt_surf=0.0;
deltahi_melt_surf=0.0;
deltahs_melt_surf=0.0;
deltahi_bott=deltat*(Fbott-Fw-ISI)/Qi_bott;

hs=hs_old+delta_bucket;
hmi_new=0.0;
hmi=hmi_old-delta_sublim;
hi=hi_old+ deltahi_bott;

end
! INITIATE FUNCTION WITH THIS:
![Fsurf(i),Fbott(i),I0(i),IM(i),ISI_10(i),ISI(i),delta_sublim(i),
!    deltahs_melt_surf(i),deltahi_melt_surf(i),deltahmi_melt_surf(i),
!    deltahi_bott(i),hs(i),hmi_new(i),hmi(i),hi(i),hs_prec_bucket(i)] = ...
!    growth_case3_snow_ice(hs_prec_bucket(i),Ks(i),Ts(i),T0(i),Kice_bott(i),
!    Tfr,Ti(i),Fs(i),ks_snow(i),hs(i-1),hmi(i-1),hi(i-1),ksi_10_av(i),
!    ksi_av(i),Fla(i),deltat,ros_new(i),Lv,Qs(i),Fw(i),Qi_bott(i),row,
!    ro_sice_surf(i-1),ros_av(i),romi,delta_bucket(i),beta,ros(i),si_fract,
!    vis_fr,ZERO);

!   CASE 3: SNOW + SEA ICE
!subroutine growth_case3_snow_ice(Fsurf,Fbott,I0,IM,ISI_10,ISI,delta_sublim, &
!    deltahs_melt_surf,deltahi_melt_surf,deltahmi_melt_surf,deltahi_bott, &
!    hs,hmi_new,hmi,hi,hs_prec_bucket_out,hs_prec_bucket_in,Ks,Ts,T0, &
!    Kice_bott,Tfr,Ti,Fs,ks_snow,hs_old,hmi_old,hi_old,ksi_10_av,ksi_av,Fla, &
!    deltat,ros_new,Lv,Qs,Fw,Qi_bott,row,ro_sice_surf_old,ros_av,romi, &
!    delta_bucket,beta,ros,si_fract,vis_fr,ZERO)
    
subroutine growth_case3_snow_ice(Fsurf,Fbott,I0,IM,ISI_10,ISI,delta_sublim, &
    deltahs_melt_surf,deltahi_melt_surf,deltahmi_melt_surf,deltahi_bott, &
    hs,hmi_new,hmi,hi,hs_prec_bucket_out,hs_prec_bucket_in,Ks,Ts,T0, &
    Kice_bott,Ti,Fs,ks_snow,hs_old,hmi_old,hi_old,ksi_10_av,ksi_av,Fla, &
    deltat,ros_new,Qs,Fw,Qi_bott,ro_sice_surf_old,ros_av,delta_bucket,ros)
    
use ice_params
use ice_global

implicit none
double precision Fsurf,Fbott,I0,IM,ISI_10,ISI,delta_sublim,ros
double precision deltahs_melt_surf,deltahi_melt_surf,deltahmi_melt_surf,deltahi_bott
double precision hs,hmi_new,hmi,hi,hs_prec_bucket_out,hs_prec_bucket_in,Ks,Ts,T0,Kice_bott
double precision Ti,Fs,ks_snow,hs_old,hmi_old,hi_old,ksi_10_av,ksi_av,Fla,deltat
double precision ros_new,Qs,Fw,Qi_bott,ro_sice_surf_old,ros_av,delta_bucket

double precision, parameter :: eps=0.1D+0 
double precision, parameter :: thsnd=1.D3

Fsurf=Ks*(Ts-T0);
Fbott=Kice_bott*(Tfr-Ti);

I0=vis_fr*Fs*exp(-ks_snow*(hs_old+hs_prec_bucket_in));
IM=I0;

if( hi_old>=eps ) then
    ISI_10=IM*exp(-ksi_10_av*eps);
    ISI=ISI_10*exp(-ksi_av*(hi_old-eps));
else if( hi_old<eps ) then
    ISI_10=IM*exp(-ksi_10_av*hi_old);
    ISI=ISI_10;
end if

delta_sublim=(Fla*deltat)/(ros_new*Lv-Qs);
deltahs_melt_surf=0.0;
deltahi_melt_surf=0.0;
deltahmi_melt_surf=0.0;
deltahi_bott=deltat*(Fbott-Fw-ISI)/Qi_bott;

!   SNOW ICE FORMATION
if( (hs_old + hs_prec_bucket_in - delta_sublim)> (hi_old*(row &
        -ro_sice_surf_old*thsnd)/ros_av - hs_prec_bucket_in*ros_new/ros) ) then
    hs=hs_old + delta_bucket -delta_sublim -beta*((-row*hi_old + ros*hs_old &
        + thsnd*ro_sice_surf_old*hi_old + ros_new*hs_prec_bucket_in)/(-romi &
        + row+ beta*ros_av));
    hmi_new=((-row*hi_old + ros*(hs_old-delta_sublim)+ thsnd*ro_sice_surf_old &
        *hi_old + ros_new*hs_prec_bucket_in)/(-romi + row+ beta*ros_av)) &
        *si_fract;
    hmi=hmi_old + hmi_new;                  
else                   
    hs=hs_old + delta_bucket - delta_sublim;
    hmi_new=0.0;
    hmi=hmi_old;
end if

!::RASA: emptying the bucket?
if( hs<ZERO ) then
    hs=hs_prec_bucket_in;
    hs_prec_bucket_out=0.0;
else   !::RASA: added to return the same value
    hs_prec_bucket_out=hs_prec_bucket_in;
end if

hi= hi_old + deltahi_bott;
    
end! INITIATE FUNCTION WITH THIS:
![Fsurf(i),Fbott(i),I0(i),IM(i),ISI_10(i),ISI(i),delta_sublim(i),
!    deltahs_melt_surf(i),deltahi_melt_surf(i),deltahmi_melt_surf(i),
!    deltahi_bott(i),hs(i),hmi_new(i),hmi(i),hi(i),R(i)] = ...
!    growth_case4_ice(T0(i),Kice_bott(i),Tfr,Ti(i),Fs(i),hs(i-1),hmi(i-1),
!    hi(i-1),ksi_10_av(i),ksi_av(i),Fla(i),deltat,Fw(i),Qi_bott(i),
!    delta_bucket(i),vis_fr,Lv,infra_fr,k_ocean_red,h_mix,k_ocean_vis,F(i),...
!    ZERO,Kice_surf(i),ro_sice_surf(i-1),Qi_surf(i),ro_sice_bott(i-1),c0,L0s,
!    Si(i-1),hi_min);

!   CASE 4: ONLY SEA ICE
!subroutine growth_case4_ice(Fsurf,Fbott,I0,IM,ISI_10,ISI,delta_sublim, &
!    deltahs_melt_surf,deltahi_melt_surf,deltahmi_melt_surf,deltahi_bott,hs, &
!    hmi_new,hmi,hi,R,T0,Kice_bott,Tfr,Ti,Fs,hs_old,hmi_old,hi_old,ksi_10_av, &
!    ksi_av,Fla,deltat,Fw,Qi_bott,delta_bucket,vis_fr,Lv,infra_fr,k_ocean_red, &
!    h_mix,k_ocean_vis,F,ZERO,Kice_surf,ro_sice_surf_old,Qi_surf, &
!    ro_sice_bott_old,c0,L0s,Si_old,hi_min)

subroutine growth_case4_ice(Fsurf,Fbott,I0,IM,ISI_10,ISI,delta_sublim, &
    deltahs_melt_surf,deltahi_melt_surf,deltahmi_melt_surf,deltahi_bott,hs, &
    hmi_new,hmi,hi,R,T0,Kice_bott,Ti,Fs,hs_old,hmi_old,hi_old,ksi_10_av, &
    ksi_av,Fla,deltat,Fw,Qi_bott,delta_bucket,F,Kice_surf,ro_sice_surf_old, &
    Qi_surf,ro_sice_bott_old,Si_old)

use ice_params
use ice_global

implicit none
double precision Fsurf,Fbott,I0,IM,ISI_10,ISI,delta_sublim,ro_sice_bott_old,Si_old
double precision deltahs_melt_surf,deltahi_melt_surf,deltahmi_melt_surf,deltahi_bott,hs
double precision hmi_new,hmi,hi,R,T0,Kice_bott,Ti,Fs,hs_old,hmi_old,hi_old,ksi_10_av
double precision ksi_av,Fla,deltat,Fw,Qi_bott,delta_bucket,F,Kice_surf,ro_sice_surf_old
double precision Qi_surf

double precision, parameter :: eps=0.1D+0
double precision, parameter :: thsnd=1.D3

Fsurf=Kice_surf*(Ti-T0);
Fbott=Kice_bott*(Tfr-Ti);

if( hi_old>=eps ) then
    ISI_10=vis_fr*Fs*exp(-ksi_10_av*eps);
    ISI=ISI_10*exp(-ksi_av*(hi_old-eps));
else if( hi_old<eps ) then
    ISI_10=vis_fr*Fs*exp(-ksi_10_av*hi_old);
    ISI=ISI_10;
end if

I0=ISI_10;
IM=ISI_10;

if( hi_old>hi_min ) then
    delta_sublim=(Fla*deltat)/(ro_sice_surf_old*Lv-Qi_surf);
    deltahi_bott=deltat*(Fbott-Fw-ISI)/Qi_bott;
else
    delta_sublim=0.0;
    ISI=Fs*((infra_fr*exp(-k_ocean_red*h_mix)) + ((1-infra_fr)*exp(-k_ocean_vis &
        *h_mix)));
    deltahi_bott=deltat*(F+ISI)/Qi_bott;
end if

deltahs_melt_surf=0.0;
deltahmi_melt_surf=0.0;
deltahi_melt_surf=0.0;

hs=hs_old+delta_bucket;
hmi_new=0.0;
hmi=hmi_old;
hi= hi_old + deltahi_bott - delta_sublim;

if( hi<ZERO ) then
    R =(ro_sice_bott_old*thsnd*((c0*(Tfr-Ti)) + L0s*(1+Si_old/(Ti-Tfrs)))) &
        *(hi-ZERO)/deltat;  ! EXTRA HEAT FOR ENERGY CONSERVATION
    hi=ZERO;
    hs=ZERO;
    hmi=ZERO;
else   !::RASA: added to return the initial value
    R = 0.0;
end if

end
! INITIATE FUNCTION WITH THIS:
![Sice_mix(i),Sice(i),Sbr_ice(i),Vbr_ice(i),Sice_5(i),Sbr_5(i),Tice_5(i),
!    Tice_bio(i)] = ...
!    growth_halodynamic_desalination(hi(i),hi(i-1),Sice_bott(i),Sice(i-1),rosa,
!    viola,Vbr_ice(i-1),Tfr,Tice(i),deltat,mu,ki_bio_bott(i-1),ki_ice_bio(i-1),
!    ki_ice_bott(i-1));
    
! DESALINATION
!subroutine growth_halodynamic_desalination(Sice_mix,Sice,Sbr_ice,Vbr_ice, &
!    Sice_5,Sbr_5,Tice_5,Tice_bio,hi,hi_old,Sice_bott,Sice_old,rosa,viola, &
!    Vbr_ice_old,Tfr,Tice,deltat,mu,ki_bio_bott_old,ki_ice_bio_old, &
!    ki_ice_bott_old)

subroutine growth_halodynamic_desalination(Sice_mix,Sice,Sbr_ice,Vbr_ice, &
    Sice_5,Sbr_5,Tice_5,Tice_bio,hi,hi_old,Sice_bott,Sice_old,Vbr_ice_old, &
    Tice,deltat,ki_bio_bott_old,ki_ice_bio_old,ki_ice_bott_old)

use ice_params
use ice_global

implicit none
double precision Sice_mix,Sice,Sbr_ice,Vbr_ice,Sice_5,Sbr_5
double precision Tice_5,Tice_bio,hi,hi_old,Sice_bott,Sice_old,Vbr_ice_old
double precision Tice,deltat,ki_bio_bott_old,ki_ice_bio_old,ki_ice_bott_old

Sice_mix=((hi-hi_old)*Sice_bott + hi_old*Sice_old)/hi;
Sice = Sice_mix + (rosa*(1-viola*Vbr_ice_old)*(Tfr-Tfrs-Tice))*deltat;
Sbr_ice=-Tice/mu;
Vbr_ice=Sice/Sbr_ice;

Sice_5=Sice;       !   DESALINATION=0, IT STARTS HERE
Sbr_5=Sbr_ice;
Tice_5=Tice;
Tice_bio=(ki_bio_bott_old*(Tfr-Tfrs_appr) + ki_ice_bio_old*Tice) &
    /(ki_bio_bott_old+ki_ice_bott_old);

end! INITIATE FUNCTION WITH THIS:
![Sice_mix(i),Sice(i),Sbr_ice(i),Vbr_ice(i),Sice_5_mix(i),Sice_5(i),Sbr_5(i),
!    Tice_5(i),Tice_bio(i)] = ...
!    growth_halodynamic_no_desalination(hi(i),hi(i-1),Sice_bott(i),Sice(i-1),
!    rosa,viola,Tfr,Tice(i),deltat,mu,ki_bio_bott(i-1),hi_5(i-1),Sice_5(i-1),
!    Tice_5(i-1),ki_ice_5(i-1));
    
! DESALINATION
!subroutine growth_halodynamic_no_desalination(Sice_mix,Sice,Sbr_ice,Vbr_ice, &
!    Sice_5_mix,Sice_5,Sbr_5,Tice_5,Tice_bio,hi,hi_old,Sice_bott,Sice_old,rosa, &
!    viola,Tfr,Tice,deltat,mu,ki_bio_bott_old,hi_5_old,Sice_5_old,Tice_5_old, &
!    ki_ice_5_old)

subroutine growth_halodynamic_no_desalination(Sice_mix,Sice,Sbr_ice,Vbr_ice, &
    Sice_5_mix,Sice_5,Sbr_5,Tice_5,Tice_bio,hi,hi_old,Sice_bott,Sice_old, &
    Tice,deltat,ki_bio_bott_old,hi_5_old,Sice_5_old,Tice_5_old,ki_ice_5_old)

use ice_params
use ice_global

implicit none
double precision Sice_mix,Sice,Sbr_ice,Vbr_ice,Sice_5_mix,Sice_5
double precision Sbr_5,Tice_5,Tice_bio,hi,hi_old,Sice_bott,Sice_old,Tice
double precision deltat,ki_bio_bott_old,hi_5_old,Sice_5_old,Tice_5_old,ki_ice_5_old

double precision, parameter :: S_p = 0.0499D+0

Sice_mix=((hi-hi_old)*Sice_bott + hi_old*Sice_old)/hi;
Sice=Sice_mix;
Sbr_ice=-Tice/mu;
Vbr_ice=Sice/Sbr_ice;    

Sice_5_mix=((hi_5_old*Sice_5_old+(hi-hi_old)*Sice_bott))/(hi-hi_old+hi_5_old);
Sice_5=Sice_5_mix + (rosa*(1-viola*S_p)*(Tfr-Tfrs-Tice_5_old))*deltat;
Sbr_5=Sice_5/S_p;
Tice_5=-mu*Sbr_5;                 
Tice_bio=(ki_bio_bott_old*(Tfr-Tfrs_appr) + ki_ice_5_old*Tice_5) &
    /(ki_bio_bott_old+ki_ice_5_old);

end
! INITIATE FUNCTION WITH THIS:
![snow_fr(i),sea_water_fr(i),Ssnowice_new(i),Ssnowice(i),Q(i),Sice(i),...
!    Sice_mix(i),Sice_5_mix(i),Sice_bio_mix(i),Sbr_ice(i),Vbr_ice(i),Sbr_5(i),
!    Sice_5(i),Tice_5(i),hi_5(i),Si(i),Vbr_i(i),Sbr_i(i),Sbr_bio(i),
!    Sice_bio(i),Tice_bio(i),Vbr_bio(i),hi_bio(i),ISI_bio(i),Sice_bott(i),
!    Sbr_bott(i),Vbr_bott(i)] = ...
!    growth_halodynamic_no_salinity(Sw(i),Tfr,ZERO);

! NO SALINITY COMPUTATION
!subroutine growth_halodynamic_no_salinity(snow_fr,sea_water_fr,Ssnowice_new, &
!    Ssnowice,Q,Sice,Sice_mix,Sice_5_mix,Sice_bio_mix,Sbr_ice,Vbr_ice,Sbr_5, &
!    Sice_5,Tice_5,hi_5,Si,Vbr_i,Sbr_i,Sbr_bio,Sice_bio,Tice_bio,Vbr_bio, &
!    hi_bio,ISI_bio,Sice_bott,Sbr_bott,Vbr_bott,Sw,Tfr,ZERO)

subroutine growth_halodynamic_no_salinity(snow_fr,sea_water_fr,Ssnowice_new, &
    Ssnowice,Q,Sice,Sice_mix,Sice_5_mix,Sice_bio_mix,Sbr_ice,Vbr_ice,Sbr_5, &
    Sice_5,Tice_5,hi_5,Si,Vbr_i,Sbr_i,Sbr_bio,Sice_bio,Tice_bio,Vbr_bio, &
    hi_bio,ISI_bio,Sice_bott,Sbr_bott,Vbr_bott,Sw)

use ice_params
use ice_global

implicit none
double precision snow_fr,sea_water_fr,Ssnowice_new,Ssnowice,Q,Sice,Sice_mix,Sice_5_mix
double precision Sice_bio_mix,Sbr_ice,Vbr_ice,Sbr_5,Sice_5,Tice_5,hi_5,Si,Vbr_i,Sbr_i
double precision Sbr_bio,Sice_bio,Tice_bio,Vbr_bio,hi_bio,ISI_bio,Sice_bott,Sbr_bott
double precision Vbr_bott,Sw

snow_fr=0.0;
sea_water_fr=0.0;
Ssnowice_new=0.0;
Ssnowice=0.0; 

Q=ZERO;
Sice=Sw;

Sice_mix=Sw;
Sice_5_mix=Sw;
Sice_bio_mix=Sw;

Sbr_ice=Sw;
Vbr_ice=1.0;

Sbr_5=Sw;
Sice_5=Sw;
Tice_5=Tfr-Tfrs;
hi_5=ZERO;

Si=Sw;
Vbr_i=1.0;
Sbr_i=Sw;

Sbr_bio=Sw;
Sice_bio=Sw;
Tice_bio=Tfr-Tfrs;
Vbr_bio=1.0;
hi_bio=hi_5;
ISI_bio=0.0;

Sice_bott=Sw;
Sbr_bott=Sw;
Vbr_bott=1.0;

end
! INITIATE FUNCTION WITH THIS:
![snow_fr(i),sea_water_fr(i),Ssnowice_new(i),Ssnowice(i)] = ...
!    growth_halodynamic_snowice(hmi_new(i),beta,ro_sice_surf(i-1),ros_av(i),
!    Sw(i),Ssnowice(i-1),hmi(i-1),hmi(i),snow_fr(i-1),sea_water_fr(i-1));

! SNOW ICE LAYE
!subroutine growth_halodynamic_snowice(snow_fr,sea_water_fr,Ssnowice_new, &
!    Ssnowice,hmi_new,beta,ro_sice_surf_old,ros_av,Sw,Ssnowice_old,hmi_old, &
!    hmi,snow_fr_old,sea_water_fr_old)

subroutine growth_halodynamic_snowice(snow_fr,sea_water_fr,Ssnowice_new, &
    Ssnowice,hmi_new,ro_sice_surf_old,ros_av,Sw,Ssnowice_old,hmi_old, &
    hmi,snow_fr_old,sea_water_fr_old)
    
use ice_params

implicit none
double precision snow_fr,sea_water_fr,Ssnowice_new
double precision Ssnowice,hmi_new,ro_sice_surf_old,ros_av,Sw,Ssnowice_old,hmi_old
double precision hmi,snow_fr_old,sea_water_fr_old

double precision, parameter :: thsnd=1.D3;

if( hmi_new>0.0 ) then
   snow_fr=beta/(ro_sice_surf_old*thsnd/ros_av);
   sea_water_fr=1-snow_fr;
   Ssnowice_new= (hmi_new*sea_water_fr*Sw)/hmi_new;
   Ssnowice= (Ssnowice_old*hmi_old+Ssnowice_new*hmi_new)/hmi;
else 
    snow_fr=snow_fr_old;
    sea_water_fr=sea_water_fr_old;
    Ssnowice_new=0.0;
    Ssnowice=Ssnowice_old;
end if

end
    subroutine freezing_temperature(Sw,Tfrz)
    
    use ice_params

    implicit none

    double precision Sw,Tfrz
    
    double precision, parameter :: P = 10.13250D+0; ! pressure, decibars
    double precision, parameter :: a = -0.0575D+0
    double precision, parameter :: b = 1.710523E-3
    double precision, parameter :: c = 2.154996E-4
    double precision, parameter :: d = 7.53E-4
    
    double precision s

    s = max(Sw,0.D+0)
    !Tfrz = Tfreez
    Tfrz = (a + b * sqrt(abs(s)) - c * s) * s - d * P + Tfreez;
    
    end
! INITIATE FUNCTION WITH THIS:
![Fs(i),Fla(i),Fld(i),Fw(i)] = ...
!    heat_fluxes(alpha(i),Fsd_cloud(i),qa(i),P0,emp,cail,roa,Lai,Ua(i),qs(i),
!    hi(i-1),hi_min,sigma,Ta(i),Cl(i),row,cpw,cwih,Tmix(i-1),Tfr);

! HEAT FLUXES NOT DEPENDENT ON T0
!subroutine heat_fluxes(Fs,Fla,Fld,Fw,alpha,Fsd_cloud,qa,P0,emp,cail,roa,Lai, &
!      Ua,qs,hi_i_1,hi_min,sigma,Ta,Cl,row,cpw,cwih,Tmix_old,Tfr)

subroutine heat_fluxes(Fs,Fla,Fld,Fw,alpha,Fsd_cloud,qa,Ua,qs,hi_old,Ta,Cl, &
    Tmix_old)
      
use ice_params
use ice_global

implicit none
double precision Fs,Fla,Fld,Fw,alpha,Fsd_cloud,qa,ea,Ua,qs,hi_old,Ta,Cl,Tmix_old

double precision, parameter :: Fld_p=85.6D+0
double precision, parameter :: Fld_p2=0.26D+0
double precision, parameter :: Fld_p3=0.653D+0
double precision, parameter :: Fld_p4=0.00535D+0
double precision, parameter :: Fld_p5=0.1762D+0
double precision, parameter :: vel=0.05D+0

! 	SHORT WAVE
Fs=(1 - alpha)*Fsd_cloud;
!   LATENT HEAT
ea=qa*P0/emp;      !::RASA: ea = air relative humidity
Fla=cail*roa*Lai*Ua*(qs-qa);
!	DOWNWARD LONG WAVE
if ( hi_old>hi_min ) then
    Fld=((sigma*(Ta**4))-Fld_p)*(1 + Fld_p2*Cl); ! (after Guest, 1997 
                                              !     for Antartica)
else  
    ! (after Bignami et al, 1995 for Mediterranean Sea)
    Fld= (sigma*Ta**4*(Fld_p3 - Fld_p4*ea))*(1 + Fld_p5*Cl**2);
    ! (after Zapadka et al, 2007 for the Baltic) 
    ! ::RASA: with this some winters become completely without ice
    !Fld=(sigma*Ta.^4*(0.685 + 0.00452*ea)).*(1 + 0.36*Cl.^2);
end if

!	OCEANIC FLUXES
!Fw=Fwater;
Fw=row*cpw*cwih*vel*(Tmix_old-Tfr);  !::RASA: same as Fw=0

end

! INITIATE FUNCTION WITH THIS:
![Flu_it,Fse_it,Fl_it,F_it,la,sen,lat,cons,coni,conmi,cpi_it,mui] = ...
!    iteration_heat_fluxes(emis,sigma,T0_star(j-1),roa,cpair,cais,Ua(i),Ta(i),
!    Fld(i),Fs(i),Fla(i),cail,Lai,c1,Tb,c2,qs(i),ks_av(i),hs(i-1),kice_surf(i),
!    hi(i-1),kmi,hmi(i-1),c0,L0s,mu,Si(i-1),Ti_star(j-1),Ti_star(1),deltat,
!    ro_sice_bulk(i-1));

! ITERATION HEAT FLUXES
!subroutine iteration_heat_fluxes(Flu_it,Fse_it,Fl_it,F_it,la,sen,lat,cons, &
!    coni,conmi,cpi_it,mui,emis,sigma,T0_star_old,roa,cpair,cais,Ua,Ta,Fld,Fs, &
!    Fla,cail,Lai,c1,Tb,c2,qs,ks_av,hs_old,kice_surf,hi_old,kmi,hmi_old,c0, &
!    L0s,mu,Si_old,Ti_star_old,Ti_star_1,deltat,ro_sice_bulk_old)

subroutine iteration_heat_fluxes(Flu_it,Fse_it,Fl_it,F_it,la,sen,lat,cons, &
    coni,conmi,cpi_it,mui,T0_star_old,Ua,Ta,Fld,Fs,Fla,qs,ks_av,hs_old, &
    kice_surf,hi_old,hmi_old,Si_old,Ti_star_old,Ti_star_1,deltat, &
    ro_sice_bulk_old)
    
use ice_params

implicit none
double precision Flu_it,Fse_it,Fl_it,F_it,la,sen,lat,cons,coni,conmi,cpi_it,mui,T0_star_old
double precision Ua,Ta,Fld,Fs,Fla,qs,ks_av,hs_old,kice_surf,hi_old,hmi_old,Si_old
double precision Ti_star_old,Ti_star_1,deltat,ro_sice_bulk_old

double precision, parameter :: thsnd=1.D3

Flu_it=emis*sigma*T0_star_old**4;           !   UPWARD LONG WAVE
Fse_it=roa*cpair*cais*Ua*(T0_star_old-Ta);  !   SENSIBLE HEAT
Fl_it= Flu_it-Fld;                          !   NET LONG WAVE
F_it= - Fs + Fl_it + Fla + Fse_it;          !   NET SURFACE FLUXES

!   LONG WAVE (differential) upward ==>negative downward
la=4*sigma*T0_star_old**3;
!   SENSIBLE HEAT (differential) upward ==>negative downward
sen=roa*cpair*cais*Ua;
!   LATENT HEAT (differential) upward ==>negative downward
lat=cail*roa*Lai*Ua*((c1*(Tb-c2)*qs)/((T0_star_old+Tb-c2)**2));
!   SNOW CONDUCTIVE (differential)
cons=ks_av/(hs_old/2);
!   SEA ICE CONDUCTIVE (differential)
coni=kice_surf/(hi_old/2);
!   SNOW ICE/SUPERIMPOSED ICE (differential)
conmi=kkmi/(hmi_old/2);
cpi_it=c0 + ((L0s *mu*Si_old)/((Ti_star_old-Tfrs)*(Ti_star_1-Tfrs)));
mui=deltat/(ro_sice_bulk_old*thsnd*cpi_it*hi_old);

end! INITIATE FUNCTION WITH THIS:
![ISI_layer(i),ISI_bio(i)]=mean_BAL_for_BFM(hi(i),hi_min,I0(i),ksi_av(i),
!     hi_5(i),IM(i),ksi_10_av(i));

! COMPUTATION OF MEAN PROPERTIES OF THE BIOLOGICAL ACTIVE PART OF THE SEA ICE 
!SISTEM FOR BFM RUN
subroutine mean_BAL_for_BFM(ISI_layer,ISI_bio,hi,I0,ksi_av,hi_5,IM, &
    ksi_10_av)

use ice_params

implicit none
double precision ISI_layer,ISI_bio,hi,I0,ksi_av,hi_5,IM,ksi_10_av

double precision, parameter :: eps=0.1D+0

if( hi>hi_min ) then
    if( hi>=eps ) then
        ISI_layer=I0*exp(-ksi_av*((hi-eps)/2));
        ISI_bio=I0*exp(-ksi_av*(hi-eps- hi_5/2));

    else if( hi<eps ) then
        ISI_layer=IM*exp(-ksi_10_av*(hi/2));
        ISI_bio=IM*exp(-ksi_10_av*(hi-hi_5/2));

    end if

else
    ISI_layer=0.0;
    ISI_bio=0.0;
end if

end! INITIATE FUNCTION WITH THIS:
![T0(i),I0(i),Fsurf(i),Fbott(i),IM(i),hmi(i),ISI_10(i),ISI(i),
!    deltahs_melt_surf(i),delta_sublim(i),deltahi_melt_surf(i),
!    deltahmi_melt_surf(i),deltahi_bott(i),hs(i),hi(i),hs_prec_bucket(i)] = 
!    melt_case1_snow(Tfr,vis_fr,Fs(i),ks_snow(i),hs(i-1),hs_prec_bucket(i),
!    Ks(i),Ts(i),Kice_bott(i),Ti(i),hmi(i-1),hmi_min,kmi_av(i),hi(i-1),
!    ksi_10_av(i),ZERO,F(i),ksi_av(i),deltat,Qs(i),Fla(i),ros_av(i),Lv,Fw(i),
!    Qi_bott(i),delta_bucket(i),ros(i),romi,ss_fract,hs_min);

! CASE 1: MELTING SNOW
!subroutine melt_case1_snow(T0,I0,Fsurf,Fbott,IM,hmi,ISI_10,ISI, &
!    deltahs_melt_surf,delta_sublim,deltahi_melt_surf,deltahmi_melt_surf, &
!    deltahi_bott,hs,hi,hs_prec_bucket_out,Tfr,vis_fr,Fs,ks_snow,hs_old, &
!    hs_prec_bucket_in,Ks,Ts,Kice_bott,Ti,hmi_old,hmi_min,kmi_av,hi_old, &
!    ksi_10_av,ZERO,F,ksi_av,deltat,Qs,Fla,ros_av,Lv,Fw,Qi_bott,delta_bucket, &
!    ros,romi,ss_fract,hs_min)

subroutine melt_case1_snow(T0,I0,Fsurf,Fbott,IM,hmi,ISI_10,ISI, &
    deltahs_melt_surf,delta_sublim,deltahi_melt_surf,deltahmi_melt_surf, &
    deltahi_bott,hs,hi,hs_prec_bucket_out,Fs,ks_snow,hs_old, &
    hs_prec_bucket_in,Ks,Ts,Kice_bott,Ti,hmi_old,kmi_av,hi_old, &
    ksi_10_av,F,ksi_av,deltat,Qs,Fla,ros_av,Fw,Qi_bott,delta_bucket,ros)

use ice_params
use ice_global

implicit none
double precision T0,I0,Fsurf,Fbott,IM,hmi,ISI_10,ISI,deltahs_melt_surf,delta_sublim
double precision deltahi_melt_surf,deltahmi_melt_surf,deltahi_bott,hs,hi
double precision hs_prec_bucket_out,Fs,ks_snow,hs_old,hs_prec_bucket_in
double precision Ks,Ts,Kice_bott,Ti,hmi_old,kmi_av,hi_old,ksi_10_av,F
double precision ksi_av,deltat,Qs,Fla,ros_av,Fw,Qi_bott,delta_bucket,ros

double precision, parameter :: eps=0.1D+0

T0=Tfr;
I0=vis_fr*Fs*exp(-ks_snow*(hs_old+hs_prec_bucket_in));
Fsurf=Ks*(Ts-T0);
Fbott=Kice_bott*(Tfr-Ti); 

if( hmi_old>hmi_min ) then
    IM=I0*exp(-kmi_av*hmi_old);
else if( hmi_old<=hmi_min ) then
    hmi=ZERO;
    IM=I0;
end if

if( hi_old>=eps ) then
    ISI_10=IM*exp(-ksi_10_av*eps);
    ISI=ISI_10*exp(-ksi_av*(hi_old-eps));
else if( hi_old<eps ) then
    ISI_10=IM*exp(-ksi_10_av*hi_old);
    ISI=ISI_10;
end if


if( Fsurf>=(F+I0) ) then
    deltahs_melt_surf=deltat*(F+I0-Fsurf)/Qs;
else
    deltahs_melt_surf=0.0;
end if

delta_sublim=(Fla*deltat)/(ros_av*Lv-Qs);
deltahi_melt_surf=0.0;
deltahmi_melt_surf=0.0;
deltahi_bott=deltat*(Fbott-Fw-ISI)/Qi_bott;

hs=hs_old+deltahs_melt_surf+delta_bucket-delta_sublim;
hmi=hmi_old- deltahs_melt_surf*ros/romi*ss_fract; ! SUPERIMPOSED ICE FORMATION
hi= hi_old +deltahi_bott;

if( hs<hs_min ) then
    hs=hs_prec_bucket_in+(hs-ZERO);          
    hs_prec_bucket_out=0.0;
else   !::RASA: added to return the initial value
    hs_prec_bucket_out=hs_prec_bucket_in;
end if
    
end
! INITIATE FUNCTION WITH THIS:
![hs(i),T0(i),I0(i),IM(i),ISI_10(i),ISI(i),Fsurf(i),Fbott(i),delta_sublim(i),
!    deltahmi_melt_surf,deltahi_melt_surf(i),deltahs_melt_surf(i),
!    deltahi_bott(i),hmi(i),hi(i),R(i)] = ...
!    melt_case2_snowice(ZERO,Tfr,vis_fr,Fs(i),kmi_av(i),hmi(i-1),ksi_10_av(i),
!    hi(i-1),ksi_av(i),Kmi(i),Tsi(i),Kice_bott(i),Ti(i),Fla(i),deltat,romi,
!    Lv,Qmi(i),F(i),Fw(i),Qi_bott(i),delta_bucket(i),Qi_surf(i));

! CASE 2: MELTING SNOW ICE
!subroutine melt_case2_snowice(hs,T0,I0,IM,ISI_10,ISI,Fsurf,Fbott, &
!    delta_sublim,deltahmi_melt_surf,deltahi_melt_surf,deltahs_melt_surf, &
!    deltahi_bott,hmi,hi,R,ZERO,Tfr,vis_fr,Fs,kmi_av,hmi_old,hi_old,ksi_10_av, &
!    ksi_av,Kmi,Tsi,Kice_bott,Ti,Fla,deltat,romi,Lv,Qmi,F,Fw,Qi_bott, &
!    delta_bucket,Qi_surf)

subroutine melt_case2_snowice(hs,T0,I0,IM,ISI_10,ISI,Fsurf,Fbott, &
    delta_sublim,deltahmi_melt_surf,deltahi_melt_surf,deltahs_melt_surf, &
    deltahi_bott,hmi,hi,R,Fs,kmi_av,hmi_old,hi_old,ksi_10_av, &
    ksi_av,Tsi,Kice_bott,Ti,Fla,deltat,Qmi,F,Fw,Qi_bott,delta_bucket,Qi_surf,Kmi)

use ice_params
use ice_global

implicit none
double precision hs,T0,I0,IM,ISI_10,ISI,Fsurf,Fbott,delta_sublim,deltahmi_melt_surf
double precision deltahi_melt_surf,deltahs_melt_surf,deltahi_bott,hmi,hi,R
double precision Fs,kmi_av,hmi_old,hi_old,ksi_10_av,ksi_av,Tsi,Kice_bott,Ti,Fla,deltat
double precision Qmi,F,Fw,Qi_bott,delta_bucket,Qi_surf,Kmi

double precision, parameter :: eps=0.1D+0

hs=ZERO;
T0=Tfr;
I0=vis_fr*Fs*exp(-kmi_av*hmi_old);
IM=I0;

if( hi_old>=eps ) then
    ISI_10=IM*exp(-ksi_10_av*eps);
    ISI=ISI_10*exp(-ksi_av*(hi_old-eps));
else if( hi_old<eps ) then
    ISI_10=IM*exp(-ksi_10_av*hi_old);
    ISI=ISI_10;
end if

Fsurf=Kmi*(Tsi-T0);
Fbott=Kice_bott*(Tfr-Ti);

delta_sublim=(Fla*deltat)/(romi*Lv-Qmi);

if( Fsurf>=(F+I0) ) then
    deltahmi_melt_surf=deltat*(F+I0-Fsurf)/Qmi;
else
    deltahmi_melt_surf=0.0;
end if

deltahi_melt_surf=0.0;
deltahs_melt_surf=0.0;
deltahi_bott=deltat*(Fbott-Fw-ISI)/Qi_bott;

hs=ZERO+delta_bucket;
hmi=hmi_old+deltahmi_melt_surf-delta_sublim;
hi=hi_old+deltahi_bott;

!   EXTRA HEAT FOR ENERGY CONSERVATION

if( hmi<ZERO ) then
    hi=hi+(Qmi*(hmi-ZERO)/Qi_surf);
    hmi=ZERO;
    if( hi<ZERO ) then
        R=(Qi_surf*(hi-ZERO))/deltat;
        hi=ZERO;
    else   !::RASA: added to return the initial value
        R = 0.0;
    end if
end if

end! INITIATE FUNCTION WITH THIS:
![hs(i),hmi(i),T0(i),ISI_10(i),ISI(i),I0(i),IM(i),Fsurf(i),Fbott(i),
!    delta_sublim(i),deltahi_bott(i),deltahi_melt_surf(i),deltahs_melt_surf(i),
!    deltahmi_melt_surf(i),hi(i),R(i)] = ...
!    melt_case3_ice(ZERO,Tfr,hi(i-1),vis_fr,Fs(i),ksi_10_av(i),ksi_av(i),
!    Kice_surf(i),Ti(i),Kice_bott(i),Fla(i),deltat,ro_sice_surf(i-1),Lv,
!    Qi_surf(i),Fw(i),Qi_bott(i),F(i),delta_bucket(i));

! CASE 3: MELTING SEA ICE
!subroutine melt_case3_ice(hs,hmi,T0,ISI_10,ISI,I0,IM,Fsurf,Fbott, &
!    delta_sublim,deltahi_bott,deltahi_melt_surf,deltahs_melt_surf, &
!    deltahmi_melt_surf,hi,R,ZERO,Tfr,hi_old,vis_fr,Fs,ksi_10_av,ksi_av, &
!    Kice_surf,Ti,Kice_bott,Fla,deltat,ro_sice_surf_old,Lv,Qi_surf,Fw,Qi_bott, &
!    F,delta_bucket)

subroutine melt_case3_ice(hs,hmi,T0,ISI_10,ISI,I0,IM,Fsurf,Fbott, &
    delta_sublim,deltahi_bott,deltahi_melt_surf,deltahs_melt_surf, &
    deltahmi_melt_surf,hi,R,hi_old,Fs,ksi_10_av,ksi_av, &
    Kice_surf,Ti,Kice_bott,Fla,deltat,ro_sice_surf_old,Qi_surf,Fw,Qi_bott, &
    F,delta_bucket)

use ice_params
use ice_global

implicit none 
double precision hs,hmi,T0,ISI_10,ISI,I0,IM,Fsurf,Fbott,delta_sublim,deltahi_bott
double precision deltahi_melt_surf,deltahs_melt_surf,deltahmi_melt_surf,hi,R
double precision hi_old,Fs,ksi_10_av,ksi_av,Kice_surf,Ti,Kice_bott,Fla,deltat
double precision ro_sice_surf_old,Qi_surf,Fw,Qi_bott,F,delta_bucket

double precision, parameter :: eps=0.1D+0

hs=ZERO;
hmi=ZERO;
T0=Tfr;
if( hi_old>=eps ) then
    ISI_10=vis_fr*Fs*exp(-ksi_10_av*eps);
    ISI=ISI_10*exp(-ksi_av*(hi_old-eps));
else if( hi_old<eps ) then
    ISI_10=vis_fr*Fs*exp(-ksi_10_av*hi_old);
    ISI=ISI_10;
end if

I0=ISI_10;
IM=I0;

Fsurf=Kice_surf*(Ti-T0);
Fbott=Kice_bott*(Tfr-Ti);

delta_sublim=(Fla*deltat)/(ro_sice_surf_old*Lv-Qi_surf);
deltahi_bott=deltat*(Fbott-Fw-ISI)/Qi_bott;

if( Fsurf>=(F+I0) ) then
    deltahi_melt_surf=deltat*(F + I0-Fsurf)/Qi_surf;
else
    deltahi_melt_surf=0.0;
end if

deltahs_melt_surf=0.0;
deltahmi_melt_surf=0.0;

hs=ZERO+delta_bucket;
hmi=ZERO;
hi= hi_old + deltahi_melt_surf +deltahi_bott -delta_sublim;

if( hi<ZERO ) then
    R=(Qi_surf*(hi-ZERO))/deltat;
    hi=ZERO;
else    !::RASA: added to return the initial value
    R=0.0;
end if

end
! INITIATE FUNCTION WITH THIS:
![ISI(i),I0(i),Fsurf(i),Fbott(i),IM(i),ISI_10(i),deltahs_melt_surf(i),
!    deltahmi_melt_surf(i),deltahi_melt_surf(i),deltahi_bott(i),
!    delta_bucket(i),hs_prec_bucket(i),hs(i),hi(i),hmi(i),Ri(i)] = ...
!    melt_case4_water(Fs(i),infra_fr,k_ocean_red,h_mix,k_ocean_vis,ZERO);

!   CASE 4: ONLY SEAWATER
!subroutine melt_case4_water(ISI,I0,Fsurf,Fbott,IM,ISI_10,deltahs_melt_surf, &
!    deltahmi_melt_surf,deltahi_melt_surf,deltahi_bott,delta_bucket, &
!    hs_prec_bucket,hs,hi,hmi,R,Fs,infra_fr,k_ocean_red,h_mix,k_ocean_vis,ZERO)

subroutine melt_case4_water(ISI,I0,Fsurf,Fbott,IM,ISI_10,deltahs_melt_surf, &
    deltahmi_melt_surf,deltahi_melt_surf,deltahi_bott,delta_bucket, &
    hs_prec_bucket,hs,hi,hmi,R,Fs,h_mix)

use ice_params

implicit none
double precision ISI,I0,Fsurf,Fbott,IM,ISI_10,deltahs_melt_surf,deltahmi_melt_surf
double precision deltahi_melt_surf,deltahi_bott,delta_bucket,hs_prec_bucket,hs,hi,hmi,R,Fs
double precision h_mix

! Penetrating solar radiation in sea water (type II, Jerlov, 1968)
ISI=Fs*((infra_fr*exp(-k_ocean_red*h_mix)) + ((1-infra_fr)*exp(-k_ocean_vis &
    *h_mix)));
I0=ISI;

Fsurf=0.0;
Fbott=0.0;

IM=ISI;
ISI_10=ISI;

deltahs_melt_surf=0.0;
deltahmi_melt_surf=0.0;
deltahi_melt_surf=0.0;
deltahi_bott=0.0;
delta_bucket=0.0;
hs_prec_bucket=0.0;

hs=ZERO;
hi=ZERO;
hmi=ZERO;

R=0.0;
    
end
! INITIATE FUNCTION WITH THIS:
![Q(i),Sice(i),Sice_bott(i),Sbr_ice(i),Vbr_ice(i),Si(i),Sbr_i(i),Vbr_i(i),
!    Sbr_5(i),Sice_5(i),Tice_5(i),hi_5(i),Sbr_bio(i),Sice_bio(i),Tice_bio(i),
!    Vbr_bio(i),hi_bio(i),Sbr_bott(i),Vbr_bott(i)] = no_salinity(ZERO,Tfr);

! NO SALINITY COMPUTATION
!subroutine no_salinity(Q,Sice,Sice_bott,Sbr_ice,Vbr_ice,Si,Sbr_i,Vbr_i,Sbr_5, &
!    Sice_5,Tice_5,hi_5,Sbr_bio,Sice_bio,Tice_bio,Vbr_bio,hi_bio,Sbr_bott, &
!    Vbr_bott,ZERO,Tfr)

subroutine no_salinity(Q,Sice,Sice_bott,Sbr_ice,Vbr_ice,Si,Sbr_i,Vbr_i,Sbr_5, &
    Sice_5,Tice_5,hi_5,Sbr_bio,Sice_bio,Tice_bio,Vbr_bio,hi_bio,Sbr_bott, &
    Vbr_bott)

use ice_params
use ice_global

implicit none
double precision Q,Sice,Sice_bott,Sbr_ice,Vbr_ice,Si,Sbr_i,Vbr_i,Sbr_5,Sice_5,Tice_5,hi_5
double precision Sbr_bio,Sice_bio,Tice_bio,Vbr_bio,hi_bio,Sbr_bott,Vbr_bott

Q=ZERO;
Sice=0.0;
Sbr_ice=0.0;
Vbr_ice=0.0;

Si=0.0;
Sbr_i=0.0;
Vbr_i=0.0;

Sbr_5=0.0;
Sice_5=0.0;
Tice_5=Tfr-Tfrs;
hi_5=ZERO;

Sbr_bio=0.0;
Sice_bio=0.0;
Tice_bio=Tfr-Tfrs;
Vbr_bio=0.0;
hi_bio=hi_5;

Sice_bott=0.0;
Sbr_bott=0.0;
Vbr_bott=0.0;

end
! INITIATE FUNCTION WITH THIS:
![delta_bucket(i),hs_prec_bucket(i),ros_av(i)] = precipitation_in_bucket(
!     hs_prec_bucket(i-1),hs_prec(i),ZERO,hs(i-1),ros(i),ros_new(i));

! PRECIPITATION IN BUCKET MODEL
subroutine precipitation_in_bucket(delta_bucket,hs_prec_bucket,ros_av, &
                   hs_prec_bucket_old,hs_prec,hs_old,ros,ros_new)

use ice_params

implicit none
double precision delta_bucket,hs_prec_bucket,ros_av,hs_prec_bucket_old
double precision hs_prec,hs_old,ros,ros_new

double precision, parameter :: eps=0.01D+0

if( hs_prec_bucket_old>eps ) then
    delta_bucket=hs_prec_bucket_old - ZERO;  ! snow in bucket
    hs_prec_bucket=hs_prec + ZERO;           ! fresh/young snow
else if ( hs_prec_bucket_old<=eps ) then
    delta_bucket=0.0;
    hs_prec_bucket=hs_prec_bucket_old + hs_prec;
end if

ros_av=(hs_old*ros+hs_prec_bucket*ros_new)/(hs_old+hs_prec_bucket);

end

! INITIATE FUNCTION WITH THIS:
![Q(i),Sbr_ice_star(i),Sice(i),Vbr_ice(i),Sbr_ice(i),Sbr_i_star(i),Si(i),
!    Vbr_i(i),Sbr_i(i),Sbr_5_star(i),Sice_5(i),Tice_5(i),Vbr_5(i),Sbr_5(i),
!    hi_5(i),Sbr_bio_star(i),Sice_bio(i),Tice_bio(i),Vbr_bio(i),Sbr_bio(i),
!    hi_bio(i),Sbr_bott_star(i),Sice_bott(i),Vbr_bott(i),Sbr_bott(i)] = ...
!    salinity_melting_ice(Vbr_ice(i-1),ro_sice_surf(i-1),deltahi_melt_surf(i),
!    deltat,Sbr_ice(i-1),ro_br(i-1),Sbr_end,Sice(i-1),mu,Tice(i),Sbr_i(i-1),
!    Vbr_i(i-1),ro_br_i(i-1),Si(i-1),Ti(i),Sbr_5(i-1),ro_br_5(i-1),
!    ki_5_bott(i-1),Tfr,ki_ice_5(i-1),hi(i),Sbr_bio(i-1),Vbr_bio(i-1),
!    ro_br_bio(i-1),Sice_bio(i-1),ki_bio_bott(i-1),ki_ice_bio(i-1),
!    ki_ice_bott(i-1),Sbr_bott(i-1),Vbr_bott(i-1),ro_br_bott(i-1),
!    Sice_bott(i-1));

! SALINITY FOR MELTING SEA ICE - FLUSHING
!subroutine salinity_melting_ice_flushing(Q,Sbr_ice_star,Sice,Vbr_ice,Sbr_ice, &
!    Sbr_i_star,Si,Vbr_i,Sbr_i,Sbr_5_star,Sice_5,Tice_5,Vbr_5,Sbr_5,hi_5, &
!    Sbr_bio_star,Sice_bio,Tice_bio,Vbr_bio,Sbr_bio,hi_bio,Sbr_bott_star, &
!    Sice_bott,Vbr_bott,Sbr_bott,Vbr_ice_old,ro_sice_surf_old, &
!    deltahi_melt_surf,deltat,Sbr_ice_old,ro_br_old,Sbr_end,Sice_old,mu,Tice, &
!    Sbr_i_old,Vbr_i_old,ro_br_i_old,Si_old,Ti,Sbr_5_old,ro_br_5_old, &
!    Sice_5_old,ki_5_bott_old,Tfr,ki_ice_5_old,hi,Sbr_bio_old,Vbr_bio_old, &
!    ro_br_bio_old,Sice_bio_old,ki_bio_bott_old,ki_ice_bio_old, &
!    ki_ice_bott_old,Sbr_bott_old,Vbr_bott_old,ro_br_bott_old,Sice_bott_old)

subroutine salinity_melting_ice_flushing(Q,Sbr_ice_star,Sice,Vbr_ice,Sbr_ice, &
    Sbr_i_star,Si,Vbr_i,Sbr_i,Sbr_5_star,Sice_5,Tice_5,Vbr_5,Sbr_5,hi_5, &
    Sbr_bio_star,Sice_bio,Tice_bio,Vbr_bio,Sbr_bio,hi_bio,Sbr_bott_star, &
    Sice_bott,Vbr_bott,Sbr_bott,Vbr_ice_old,ro_sice_surf_old, &
    deltahi_melt_surf,deltat,Sbr_ice_old,ro_br_old,Sice_old,Tice, &
    Sbr_i_old,Vbr_i_old,ro_br_i_old,Si_old,Ti,Sbr_5_old,ro_br_5_old, &
    Sice_5_old,ki_5_bott_old,ki_ice_5_old,hi,Sbr_bio_old,Vbr_bio_old, &
    ro_br_bio_old,Sice_bio_old,ki_bio_bott_old,ki_ice_bio_old, &
    ki_ice_bott_old,Sbr_bott_old,Vbr_bott_old,ro_br_bott_old,Sice_bott_old)
    
use ice_params
use ice_global

implicit none
double precision Q,Sbr_ice_star,Sice,Vbr_ice,Sbr_ice,Sbr_i_star,Si,Vbr_i,Sbr_i,Sbr_5_star
double precision Sice_5,Tice_5,Vbr_5,Sbr_5,hi_5,Sbr_bio_star,Sice_bio,Tice_bio,Vbr_bio
double precision Sbr_bio,hi_bio,Sbr_bott_star,Sice_bott,Vbr_bott,Sbr_bott,Vbr_ice_old
double precision ro_sice_surf_old,deltahi_melt_surf,deltat,Sbr_ice_old,ro_br_old
double precision Sice_old,Tice,Sbr_i_old,Vbr_i_old,ro_br_i_old,Si_old,Ti,Sbr_5_old
double precision ro_br_5_old,Sice_5_old,ki_5_bott_old,ki_ice_5_old,hi,Sbr_bio_old
double precision Vbr_bio_old,ro_br_bio_old,Sice_bio_old,ki_bio_bott_old,ki_ice_bio_old
double precision ki_ice_bott_old,Sbr_bott_old,Vbr_bott_old,ro_br_bott_old,Sice_bott_old

double precision, parameter :: thsnd=1.D3
double precision, parameter :: Q_p=0.20D+0
double precision, parameter :: S_p=0.05D+0

!   Salinity of Q is 0      ! 0.80 is refreezing, 0.20 is percolating
Q=Q_p*ro_sice_surf_old*(-deltahi_melt_surf)/deltat;
Sbr_ice_star = Sbr_ice_old - ((Q/(Vbr_ice_old*ro_br_old*thsnd))*(Sbr_ice_old &
    - Sbr_end))*deltat;
Sice= Sice_old + (Sbr_ice_star-Sbr_ice_old)*Vbr_ice_old;
Vbr_ice=-mu*Sice/Tice;
Sbr_ice=Sice/Vbr_ice;

Sbr_i_star = Sbr_i_old - ((Q/(Vbr_i_old*ro_br_i_old*thsnd))*(Sbr_i_old &
    - Sbr_end))*deltat;
Si= Si_old + (Sbr_i_star-Sbr_i_old)*Vbr_i_old;
Vbr_i=-mu*Si/(Ti-Tfrs);
Sbr_i=Si/Vbr_i;

Sbr_5_star=Sbr_5_old - ((Q/(S_p*ro_br_5_old*thsnd))*(Sbr_5_old - Sbr_end)) &
    *deltat;
Sice_5=Sice_5_old + (Sbr_5_star-Sbr_5_old)*S_p;
Tice_5=(ki_5_bott_old*(Tfr-Tfrs) + ki_ice_5_old*Tice)/(ki_5_bott_old &
    +ki_ice_5_old);
Vbr_5=-mu*Sice_5/Tice_5;
Sbr_5=Sice_5/Vbr_5;

hi_5=(hi*ki_5_bott_old*(Tice_5-(Tfr-Tfrs)))/((ki_ice_5_old*(Tice-Tice_5)) &
    +(ki_5_bott_old*(Tice_5-(Tfr-Tfrs))));
Sbr_bio_star=Sbr_bio_old-((Q/(Vbr_bio_old*ro_br_bio_old*thsnd))*(Sbr_bio_old &
    - Sbr_end))*deltat;       !   DESALINATION IN THE ACTIVE BIOLOCIAL SYSTEM
Sice_bio=Sice_bio_old + (Sbr_bio_star-Sbr_bio_old)*Vbr_bio_old;
Tice_bio=(ki_bio_bott_old*(Tfr-Tfrs) + ki_ice_bio_old*Tice) &
    /(ki_bio_bott_old+ki_ice_bott_old);
Vbr_bio=-mu*Sice_bio/Tice_bio;
Sbr_bio=Sice_bio/Vbr_bio;
hi_bio=hi_5;

Sbr_bott_star=Sbr_bott_old-((Q/(Vbr_bott_old*ro_br_bott_old*thsnd)) &
    *(Sbr_bott_old - Sbr_end))*deltat;  !   DESALINATION IN THE BOTTOM LAYER
Sice_bott=Sice_bott_old + (Sbr_bott_star-Sbr_bott_old)*Vbr_bott_old;
Vbr_bott=-mu*Sice_bott/(Tfr-Tfrs_appr);
Sbr_bott=Sice_bott/Vbr_bott;

end
! INITIATE FUNCTION WITH THIS:
![Q(i),Sice(i),Vbr_ice(i),Sbr_ice(i),Si(i),Vbr_i(i),Sbr_i(i),Sice_5(i),
!    Tice_5(i),Sbr_5(i),hi_5(i),Sice_bio(i),Tice_bio(i),Vbr_bio(i),Sbr_bio(i),
!    hi_bio(i),Sice_bott(i),Vbr_bott(i),Sbr_bott(i)] = ...
!    salinity_melting_ice_no_flushing(Sice(i-1),mu,Tice(i),Si(i-1),Ti(i),
!    Sice_5(i-1),ki_5_bott(i-1),Tfr,ki_ice_5(i-1),hi(i),Sice_bio(i-1),
!    ki_bio_bott(i-1),Sice_bott(i-1),ki_5_bio(i-1));

! SALINITY FOR MELTING SEA ICE - NO FLUSHING
!subroutine salinity_melting_ice_no_flushing(Q,Sice,Vbr_ice,Sbr_ice,Si,Vbr_i, &
!    Sbr_i,Sice_5,Tice_5,Sbr_5,hi_5,Sice_bio,Tice_bio,Vbr_bio,Sbr_bio,hi_bio, &
!    Sice_bott,Vbr_bott,Sbr_bott,Sice_old,mu,Tice,Si_old,Ti,Sice_5_old, &
!    ki_5_bott_old,Tfr,ki_ice_5_old,hi,Sice_bio_old,ki_bio_bott_old, &
!    Sice_bott_old,ki_5_bio_old)

subroutine salinity_melting_ice_no_flushing(Q,Sice,Vbr_ice,Sbr_ice,Si,Vbr_i, &
    Sbr_i,Sice_5,Tice_5,Sbr_5,hi_5,Sice_bio,Tice_bio,Vbr_bio,Sbr_bio,hi_bio, &
    Sice_bott,Vbr_bott,Sbr_bott,Sice_old,Tice,Si_old,Ti,Sice_5_old, &
    ki_5_bott_old,ki_ice_5_old,hi,Sice_bio_old,ki_bio_bott_old, &
    Sice_bott_old,ki_5_bio_old)
    
use ice_params
use ice_global

implicit none
double precision Q,Sice,Vbr_ice,Sbr_ice,Si,Vbr_i,Sbr_i,Sice_5,Tice_5,Sbr_5,hi_5,Sice_bio
double precision Tice_bio,Vbr_bio,Sbr_bio,hi_bio,Sice_bott,Vbr_bott,Sbr_bott,Sice_old
double precision Tice,Si_old,Ti,Sice_5_old,ki_5_bott_old,ki_ice_5_old,hi,Sice_bio_old
double precision ki_bio_bott_old,Sice_bott_old,ki_5_bio_old

double precision, parameter :: S_p=0.05D+0

Q=0.0;
Sice=Sice_old;
Sbr_ice=-Tice/mu;
Vbr_ice=Sice/Sbr_ice;

Si=Si_old;
Sbr_i=-(Ti-Tfrs)/mu;
Vbr_i=Si/Sbr_i;

Sice_5=Sice_5_old;
Sbr_5=Sice_5/S_p;
Tice_5=-mu*Sbr_5;

hi_5=(hi*ki_5_bott_old*(Tice_5-(Tfr-Tfrs)))/((ki_ice_5_old*(Tice-Tice_5)) &
    +(ki_5_bott_old*(Tice_5-(Tfr-Tfrs))));
Tice_bio=(ki_bio_bott_old*(Tfr-Tfrs) + ki_5_bio_old*Tice_5) &
    /(ki_bio_bott_old+ki_5_bio_old);
Sice_bio=Sice_bio_old;       !   DESALINATION IN THE ACTIVE BIOLOCIAL SYSTEM
Sbr_bio=-Tice_bio/mu;
Vbr_bio=Sice_bio/Sbr_bio;
hi_bio=hi_5;

Sice_bott=Sice_bott_old;
Sbr_bott=-(Tfr-Tfrs_appr)/mu;
Vbr_bott=Sice_bott/Sbr_bott;

end! INITIATE FUNCTION WITH THIS:
![Q(i),Sbr_ice_star(i),Sice(i),Vbr_ice(i),Sbr_ice(i),Sbr_i_star(i),Si(i),
!    Vbr_i(i),Sbr_i(i),Sbr_5_star(i),Sice_5(i),Tice_5(i),Sbr_5(i),
!    Sbr_bio_star(i),Sice_bio(i),Tice_bio(i),Vbr_bio(i),Sbr_bio(i),hi_bio(i),
!    Sbr_bott_star(i),Sice_bott(i),Vbr_bott(i),Sbr_bott(i),hi_5(i)] = ...
!    salinity_melting_snow_flushing(Vbr_ice(i-1),deltahs_melt_surf(i),deltat,
!    Sbr_ice(i-1),ro_br(i-1),Sbr_end,Sice(i-1),mu,Tice(i),Sbr_i(i-1),Vbr_i(i-1),
!    ro_br_i(i-1),Si(i-1),Ti(i),Sbr_5(i-1),ro_br_5(i-1),Sice_5(i-1),
!    Sbr_bio(i-1),Vbr_bio(i-1),ro_br_bio(i-1),Sice_bio(i-1),ki_bio_bott(i-1),
!    Tfr,ki_ice_bio(i-1),ki_ice_bott(i-1),Sbr_bott(i-1),Vbr_bott(i-1),
!    ro_br_bott(i-1),Sice_bott(i-1),hi(i),ki_5_bott(i-1),ki_ice_5(i-1),ros(i));

! SALINITY FOR MELTING SNOW - FLUSHING
!subroutine salinity_melting_snow_flushing(Q,Sbr_ice_star,Sice,Vbr_ice, &
!    Sbr_ice,Sbr_i_star,Si,Vbr_i,Sbr_i,Sbr_5_star,Sice_5,Tice_5,Sbr_5, &
!    Sbr_bio_star,Sice_bio,Tice_bio,Vbr_bio,Sbr_bio,hi_bio,Sbr_bott_star, &
!    Sice_bott,Vbr_bott,Sbr_bott,hi_5,Vbr_ice_old,deltahs_melt_surf,deltat, &
!    Sbr_ice_old,ro_br_old,Sbr_end,Sice_old,mu,Tice,Sbr_i_old,Vbr_i_old, &
!    ro_br_i_old,Si_old,Ti,Sbr_5_old,ro_br_5_old,Sice_5_old,Sbr_bio_old, &
!    Vbr_bio_old,ro_br_bio_old,Sice_bio_old,ki_bio_bott_old,Tfr, &
!    ki_ice_bio_old,ki_ice_bott_old,Sbr_bott_old,Vbr_bott_old,ro_br_bott_old, &
!    Sice_bott_old,hi,ki_5_bott_old,ki_ice_5_old,ros)

subroutine salinity_melting_snow_flushing(Q,Sbr_ice_star,Sice,Vbr_ice, &
    Sbr_ice,Sbr_i_star,Si,Vbr_i,Sbr_i,Sbr_5_star,Sice_5,Tice_5,Sbr_5, &
    Sbr_bio_star,Sice_bio,Tice_bio,Vbr_bio,Sbr_bio,hi_bio,Sbr_bott_star, &
    Sice_bott,Vbr_bott,Sbr_bott,hi_5,Vbr_ice_old,deltahs_melt_surf,deltat, &
    Sbr_ice_old,ro_br_old,Sice_old,Tice,Sbr_i_old,Vbr_i_old, &
    ro_br_i_old,Si_old,Ti,Sbr_5_old,ro_br_5_old,Sice_5_old,Sbr_bio_old, &
    Vbr_bio_old,ro_br_bio_old,Sice_bio_old,ki_bio_bott_old, &
    ki_ice_bio_old,ki_ice_bott_old,Sbr_bott_old,Vbr_bott_old,ro_br_bott_old, &
    Sice_bott_old,hi,ki_5_bott_old,ki_ice_5_old,ros)
    
use ice_params
use ice_global

implicit none
double precision Q,Sbr_ice_star,Sice,Vbr_ice,Sbr_ice,Sbr_i_star,Si,Vbr_i,Sbr_i,Sbr_5_star
double precision Sice_5,Tice_5,Sbr_5,Sbr_bio_star,Sice_bio,Tice_bio,Vbr_bio,Sbr_bio,hi_bio
double precision Sbr_bott_star,Sice_bott,Vbr_bott,Sbr_bott,hi_5,Vbr_ice_old
double precision deltahs_melt_surf,deltat,Sbr_ice_old,ro_br_old,Sice_old,Tice
double precision Sbr_i_old,Vbr_i_old,ro_br_i_old,Si_old,Ti,Sbr_5_old,ro_br_5_old,Sice_5_old
double precision Sbr_bio_old,Vbr_bio_old,ro_br_bio_old,Sice_bio_old,ki_bio_bott_old
double precision ki_ice_bio_old,ki_ice_bott_old,Sbr_bott_old,Vbr_bott_old,ro_br_bott_old
double precision Sice_bott_old,hi,ki_5_bott_old,ki_ice_5_old,ros

double precision, parameter :: Q_p=0.20D+0
double precision, parameter :: S_p=0.05D+0
double precision, parameter :: thsnd=1.D3

!   Salinity of Q is 0      !   0.80 is refreezing, 0.20 is percolating
Q=Q_p*ros*(-deltahs_melt_surf)/deltat;
Sbr_ice_star = Sbr_ice_old - ((Q/(Vbr_ice_old*ro_br_old*thsnd))*(Sbr_ice_old &
    - Sbr_end))*deltat;
! brine volume fixed in this process, it changes only as function of 
! temperature in the next equation for next timestep
Sice= Sice_old + (Sbr_ice_star-Sbr_ice_old)*Vbr_ice_old;                                                  
Vbr_ice=-mu*Sice/Tice;
Sbr_ice=Sice/Vbr_ice;

Sbr_i_star=Sbr_i_old - ((Q/(Vbr_i_old*ro_br_i_old*thsnd))*(Sbr_i_old &
    - Sbr_end))*deltat;
Si= Si_old + (Sbr_i_star-Sbr_i_old)*Vbr_i_old;
Vbr_i=-mu*Si/(Ti-Tfrs);
Sbr_i=Si/Vbr_i;

Sbr_5_star=Sbr_5_old - ((Q/(S_p*ro_br_5_old*thsnd))*(Sbr_5_old - Sbr_end)) &
    *deltat;
Sice_5=Sice_5_old + (Sbr_5_star-Sbr_5_old)*S_p;               
Tice_5=-mu*Sice_5/S_p;
Sbr_5=Sice_5/S_p;

hi_5=(hi*ki_5_bott_old*(Tice_5-(Tfr-Tfrs)))/((ki_ice_5_old*(Tice-Tice_5)) &
    +(ki_5_bott_old*(Tice_5-(Tfr-Tfrs))));
Sbr_bio_star=Sbr_bio_old-((Q/(Vbr_bio_old*ro_br_bio_old*thsnd))*(Sbr_bio_old &
    - Sbr_end))*deltat;        ! DESALINATION IN THE ACTIVE BIOLOCIAL SYSTEM
Sice_bio=Sice_bio_old + (Sbr_bio_star-Sbr_bio_old)*Vbr_bio_old;
Tice_bio=(ki_bio_bott_old*(Tfr-Tfrs) + ki_ice_bio_old*Tice) &
    /(ki_bio_bott_old+ki_ice_bott_old);
Vbr_bio=-mu*Sice_bio/Tice_bio;
Sbr_bio=Sice_bio/Vbr_bio;
hi_bio=hi_5;

Sbr_bott_star=Sbr_bott_old-((Q/(Vbr_bott_old*ro_br_bott_old*thsnd)) &
    *(Sbr_bott_old - Sbr_end))*deltat;   !   DESALINATION IN THE BOTTOM LAYER
Sice_bott=Sice_bott_old + (Sbr_bott_star-Sbr_bott_old)*Vbr_bott_old;
Vbr_bott=-mu*Sice_bott/(Tfr-Tfrs_appr);
Sbr_bott=Sice_bott/Vbr_bott;

end
! INITIATE FUNCTION WITH THIS:
![Q(i),Sice(i),Vbr_ice(i),Sbr_ice(i),Si(i),Vbr_i(i),Sbr_i(i),Sice_5(i),
!    Tice_5(i),Sbr_5(i),Sice_bio(i),Tice_bio(i),Vbr_bio(i),Sbr_bio(i),
!    hi_bio(i),Sice_bott(i),Vbr_bott(i),Sbr_bott(i),hi_5(i)] = ...
!    salinity_melting_snowice_no_flushing(Sice(i-1),mu,Tice(i),Si(i-1),Ti(i),
!    Sice_5(i-1),Sice_bio(i-1),ki_bio_bott(i-1),Tfr,Sice_bott(i-1),hi(i),
!    ki_5_bott(i-1),ki_ice_5(i-1));

! SALINITY FOR MELTING SNOW/SNOW ICE - NO FLUSHING
!subroutine salinity_melting_snow_snowice_no_flushing(Q,Sice,Vbr_ice,Sbr_ice, &
!    Si,Vbr_i,Sbr_i,Sice_5,Tice_5,Sbr_5,Sice_bio,Tice_bio,Vbr_bio,Sbr_bio, &
!    hi_bio,Sice_bott,Vbr_bott,Sbr_bott,hi_5,Sice_old,mu,Tice,Si_old,Ti, &
!    Sice_5_old,Sice_bio_old,ki_bio_bott_old,Tfr,Sice_bott_old,hi, &
!    ki_5_bott_old,ki_ice_5_old)

subroutine salinity_melting_snow_snowice_no_flushing(Q,Sice,Vbr_ice,Sbr_ice, &
    Si,Vbr_i,Sbr_i,Sice_5,Tice_5,Sbr_5,Sice_bio,Tice_bio,Vbr_bio,Sbr_bio, &
    hi_bio,Sice_bott,Vbr_bott,Sbr_bott,hi_5,Sice_old,Tice,Si_old,Ti, &
    Sice_5_old,Sice_bio_old,ki_bio_bott_old,Sice_bott_old,hi, &
    ki_5_bott_old,ki_ice_5_old)

use ice_params
use ice_global

implicit none
double precision Q,Sice,Vbr_ice,Sbr_ice,Si,Vbr_i,Sbr_i,Sice_5,Tice_5,Sbr_5,Sice_bio
double precision Tice_bio,Vbr_bio,Sbr_bio,hi_bio,Sice_bott,Vbr_bott,Sbr_bott,hi_5,Sice_old
double precision Tice,Si_old,Ti,Sice_5_old,Sice_bio_old,ki_bio_bott_old
double precision Sice_bott_old,hi,ki_5_bott_old,ki_ice_5_old

double precision, parameter :: S_p=0.05D+0

Q=0.0;
Sice=Sice_old;
Sbr_ice=-Tice/mu;
Vbr_ice=Sice/Sbr_ice;

Si=Si_old;
Sbr_i=-(Ti-Tfrs)/mu;
Vbr_i=Si/Sbr_i;

Sice_5=Sice_5_old;
Sbr_5=Sice_5/S_p;
Tice_5=-mu*Sbr_5;

hi_5=(hi*ki_5_bott_old*(Tice_5-(Tfr-Tfrs)))/((ki_ice_5_old*(Tice-Tice_5)) &
    +(ki_5_bott_old*(Tice_5-(Tfr-Tfrs))));
Tice_bio=(ki_bio_bott_old*(Tfr-Tfrs) + ki_ice_5_old*Tice_5) &
    /(ki_bio_bott_old+ki_ice_5_old);
Sice_bio=Sice_bio_old;     !   DESALINATION IN THE ACTIVE BIOLOCIAL SYSTEM
Sbr_bio=-Tice_bio/mu;
Vbr_bio=Sice_bio/Sbr_bio;
hi_bio=hi_5;

Sice_bott=Sice_bott_old;
Sbr_bott=-(Tfr-Tfrs_appr)/mu;
Vbr_bott=Sice_bott/Sbr_bott;

end! INITIATE FUNCTION WITH THIS:
![Q(i),Sbr_ice_star(i),Sice(i),Vbr_ice(i),Sbr_ice(i),Sbr_i_star(i),Si(i),
!    Vbr_i(i),Sbr_i(i),Sbr_5_star(i),Sice_5(i),Tice_5(i),Sbr_5(i),
!    Sbr_bio_star(i),Sice_bio(i),Tice_bio(i),Vbr_bio(i),Sbr_bio(i),
!    hi_bio(i),Sbr_bott_star(i),Sice_bott(i),Vbr_bott(i),Sbr_bott(i)] = ...
!    salinity_melting_snowice(hi_5(i),Vbr_ice(i-1),romi,deltahmi_melt_surf(i),
!    deltat,Sbr_ice(i-1),ro_br(i-1),Sbr_end,Sice(i-1),mu,Tice(i),Sbr_i(i-1),
!    Vbr_i(i-1),ro_br_i(i-1),Si(i-1),Ti(i),Sbr_5(i-1),ro_br_5(i-1),Sice_5(i-1),
!    Sbr_bio(i-1),Vbr_bio(i-1),ro_br_bio(i-1),Sice_bio(i-1),ki_bio_bott(i-1),
!    Tfr,ki_ice_bio(i-1),ki_ice_bott(i-1),Sbr_bott(i-1),Vbr_bott(i-1),
!    ro_br_bott(i-1),Sice_bott(i-1));

! SALINITY FOR MELTING SNOW ICE - FLUSHING
!subroutine salinity_melting_snowice_flushing(Q,Sbr_ice_star,Sice,Vbr_ice, &
!    Sbr_ice,Sbr_i_star,Si,Vbr_i,Sbr_i,Sbr_5_star,Sice_5,Tice_5,Sbr_5, &
!    Sbr_bio_star,Sice_bio,Tice_bio,Vbr_bio,Sbr_bio,hi_bio,Sbr_bott_star, &
!    Sice_bott,Vbr_bott,Sbr_bott,hi_5,Vbr_ice_old,romi,deltahmi_melt_surf, &
!    deltat,Sbr_ice_old,ro_br_old,Sbr_end,Sice_old,mu,Tice,Sbr_i_old, &
!    Vbr_i_old,ro_br_i_old,Si_old,Ti,Sbr_5_old,ro_br_5_old,Sice_5_old, &
!    Sbr_bio_old,Vbr_bio_old,ro_br_bio_old,Sice_bio_old,ki_bio_bott_old,Tfr, &
!    ki_ice_bio_old,ki_ice_bott_old,Sbr_bott_old,Vbr_bott_old,ro_br_bott_old, &
!    Sice_bott_old)

subroutine salinity_melting_snowice_flushing(Q,Sbr_ice_star,Sice,Vbr_ice, &
    Sbr_ice,Sbr_i_star,Si,Vbr_i,Sbr_i,Sbr_5_star,Sice_5,Tice_5,Sbr_5, &
    Sbr_bio_star,Sice_bio,Tice_bio,Vbr_bio,Sbr_bio,hi_bio,Sbr_bott_star, &
    Sice_bott,Vbr_bott,Sbr_bott,hi_5,Vbr_ice_old,deltahmi_melt_surf, &
    deltat,Sbr_ice_old,ro_br_old,Sice_old,Tice,Sbr_i_old, &
    Vbr_i_old,ro_br_i_old,Si_old,Ti,Sbr_5_old,ro_br_5_old,Sice_5_old, &
    Sbr_bio_old,Vbr_bio_old,ro_br_bio_old,Sice_bio_old,ki_bio_bott_old, &
    ki_ice_bio_old,ki_ice_bott_old,Sbr_bott_old,Vbr_bott_old,ro_br_bott_old, &
    Sice_bott_old)

use ice_params
use ice_global

implicit none
double precision Q,Sbr_ice_star,Sice,Vbr_ice,Sbr_ice,Sbr_i_star,Si,Vbr_i,Sbr_i,Sbr_5_star
double precision Sice_5,Tice_5,Sbr_5,Sbr_bio_star,Sice_bio,Tice_bio,Vbr_bio,Sbr_bio,hi_bio
double precision Sbr_bott_star,Sice_bott,Vbr_bott,Sbr_bott,hi_5,Vbr_ice_old
double precision deltahmi_melt_surf,deltat,Sbr_ice_old,ro_br_old,Sice_old,Tice
double precision Sbr_i_old,Vbr_i_old,ro_br_i_old,Si_old,Ti,Sbr_5_old,ro_br_5_old,Sice_5_old
double precision Sbr_bio_old,Vbr_bio_old,ro_br_bio_old,Sice_bio_old,ki_bio_bott_old
double precision ki_ice_bio_old,ki_ice_bott_old,Sbr_bott_old,Vbr_bott_old,ro_br_bott_old
double precision Sice_bott_old

double precision, parameter :: Q_p=0.20D+0
double precision, parameter :: S_p=0.05D+0
double precision, parameter :: thsnd=1.D3

! Salinity of Q is 0
! 0.70 is refreezing, 0.20 is percolating, 0.10 is sublimating
Q=Q_p*romi*(-deltahmi_melt_surf)/deltat;
Sbr_ice_star = Sbr_ice_old - ((Q/(Vbr_ice_old*ro_br_old*thsnd))*(Sbr_ice_old  &
    - Sbr_end))*deltat;
Sice= Sice_old + (Sbr_ice_star-Sbr_ice_old)*Vbr_ice_old;
Vbr_ice=-mu*Sice/Tice;
Sbr_ice=Sice/Vbr_ice;

Sbr_i_star = Sbr_i_old - ((Q/(Vbr_i_old*ro_br_i_old*thsnd))*(Sbr_i_old &
    - Sbr_end))*deltat;
Si= Si_old + (Sbr_i_star-Sbr_i_old)*Vbr_i_old;
Vbr_i=-mu*Si/(Ti-Tfrs);
Sbr_i=Si/Vbr_i;

Sbr_5_star=Sbr_5_old - ((Q/(S_p*ro_br_5_old*thsnd))*(Sbr_5_old - Sbr_end)) &
    *deltat;
Sice_5=Sice_5_old + (Sbr_5_star-Sbr_5_old)*S_p;
Tice_5=-mu*Sice_5/S_p;
Sbr_5=Sice_5/S_p;

Sbr_bio_star=Sbr_bio_old-((Q/(Vbr_bio_old*ro_br_bio_old*thsnd))*(Sbr_bio_old &
    - Sbr_end))*deltat;    !   DESALINATION IN THE ACTIVE BIOLOCIAL SYSTEM
Sice_bio=Sice_bio_old + (Sbr_bio_star-Sbr_bio_old)*Vbr_bio_old;
Tice_bio=(ki_bio_bott_old*(Tfr-Tfrs) + ki_ice_bio_old*Tice) &
    /(ki_bio_bott_old+ki_ice_bott_old);
Vbr_bio=-mu*Sice_bio/Tice_bio;
Sbr_bio=Sice_bio/Vbr_bio;
hi_bio=hi_5;

Sbr_bott_star=Sbr_bott_old-((Q/(Vbr_bott_old*ro_br_bott_old*thsnd)) &
    *(Sbr_bott_old - Sbr_end))*deltat;  !   DESALINATION IN THE BOTTOM LAYER
Sice_bott=Sice_bott_old + (Sbr_bott_star-Sbr_bott_old)*Vbr_bott_old;
Vbr_bott=-mu*Sice_bott/(Tfr-Tfrs_appr);
Sbr_bott=Sice_bott/Vbr_bott;

end! INITIATE FUNCTION WITH THIS:
![ros(i),ros_new(i)] = snow_density(T0(i-1),Tfr);

! SNOW DENSITY
subroutine snow_density(T0_old,ros,ros_new)

use ice_params
use ice_global

implicit none
double precision T0_old,ros,ros_new

if( T0_old<Tfr ) then
    ros=350.0D+0;      ! ::RASA: snow density ("cold" snow - freezing) 300
    ros_new=300.0D+0;  ! ::RASA: new snow density ("cold" new snow - freezing) 250
else
    ros=400.0D+0;      ! ::RASA: snow density ("warm" snow - melting) 350
    ros_new=350.0D+0;  ! ::RASA: new snow density ("cold" new snow - warming) 300
end if

end

! INITIATE FUNCTION WITH THIS:
![hs_prec(i)] = snow_precipitation(Ta(i),Tfr,hi_old(i-1),hi_min,P_rate(i),
!     ros_prec,deltat);

!  SNOW PRECIPITATION
!subroutine snow_precipitation(hs_prec,Ta,Tfr,hi_old,hi_min,P_rate,ros_prec, &
!    deltat)

subroutine snow_precipitation(hs_prec,Ta,hi_old,P_rate,deltat)

use ice_params
use ice_global

implicit none

!double precision hs_prec,Ta,Tfr,hi_old,hi_min,P_rate,ros_prec,deltat
double precision hs_prec,Ta,hi_old,P_rate,deltat

if( Ta<Tfr ) then 
    if( hi_old>hi_min ) then 
        hs_prec=(P_rate/ros_prec)*deltat;
    else
        hs_prec=0.0;
    end if
else
    hs_prec=0.0;
end if

end

! INITIATE FUNCTION WITH THIS:
![ks(i),ks_av(i)] = snow_thermal_conductivity(ros_new(i),ros(i),hs(i-1),
!      hs_prec_bucket(i));

! THERMAL CONDUCTIVITY OF NEW/OLD SNOW
subroutine snow_thermal_conductivity(ks,ks_av,ros_new,ros,hs_old,hs_prec_bucket)

use ice_params

implicit none
double precision ks,ks_av,ros_new,ros,hs_old,hs_prec_bucket

double precision ks_new

double precision, parameter :: thsndm=1.D-3
double precision, parameter :: ks_p=2.85D+0

ks_new=(ros_new*thsndm)**2*ks_p; !(after Abel, 1892) thermal conductivity
                                 !  of new snow (using new snow density)
ks=(ros*thsndm)**2*ks_p;         !(after Abel, 1892) thermal conductivity 
                                 !  of snow (using snow density)
ks_av=(ks*hs_old+ks_new*hs_prec_bucket)/(hs_old+hs_prec_bucket);

end

! INITIATE FUNCTION WITH THIS:
![Flu(i),Fl(i),Fse(i),F(i),R(i),Qs(i),Qmi(i),Qi_surf(i),Qi_bott(i)] = ...
!    surface_fluxes_recompute(emis,sigma,T0(i),Fld(i),roa,cpair,cais,Ua(i),...
!    Ta(i),Fs(i),Fla(i),ros(i),c0,Tfrs,Ts(i),L0s,romi,Tsi(i),...
!    ro_sice_surf(i-1),Ti(i),mu,Si(i-1),ro_sice_bott(i-1),Tfr);

!   RECOMPUTE SURFACE FLUXES
subroutine surface_fluxes_recompute(Flu,Fl,Fse,F,R,Qs,Qmi,Qi_surf,Qi_bott, &
      T0,Fld,Ua,Ta,Fs,Fla,ros,Ts,Tsi,ro_sice_surf_old,Ti,Si_old, &
      ro_sice_bott_old)

use ice_params
use ice_global

implicit none
double precision Flu,Fl,Fse,F,R,Qs,Qmi,Qi_surf,Qi_bott,T0,Fld,ro_sice_bott_old
double precision Ua,Ta,Fs,Fla,ros,Ts,Tsi,ro_sice_surf_old,Ti,Si_old

double precision, parameter :: thsnd=1.D3

Flu=emis*sigma*T0**4;                 !   UPWARD LONGWAVE RADIATION
Fl= Flu-Fld;                          !   NET LONGWAVE RADIATION
Fse=roa*cpair*cais*Ua*(T0-Ta);        !   SENSIBLE HEAT
F=-Fs + Fl +Fla +Fse;                 !   NET SURFACE FLUXES
R=0.0;                                !   INITIAL EXTRA HEAT IS SET TO ZERO

Qs=ros*((c0*(Tfrs-Ts)) + L0s);
Qmi=romi*((c0*(Tfrs-Tsi)) + L0s);
Qi_surf=ro_sice_surf_old*thsnd*((c0*(Tfr-Ti)) + L0s*(1+mu*Si_old/(Ti-Tfrs)));
Qi_bott=ro_sice_bott_old*thsnd*((c0*(Tfr-Ti)) + L0s*(1+mu*Si_old/(Ti-Tfrs)));

end! INITIATE FUNCTION WITH THIS:
![T0_iter,Ts_iter,Tsi_iter,Ti_iter,Tmix_star,T0_star(j),Ti_star(j)] = ...
!    temperature_case1_snow_meteoric_ice(vis_fr,Fs(i),ks_snow(i),hs(i-1),
!    hs_prec_bucket(i),kmi_av(i),hmi(i-1),hi(i-1),ksi_10_av(i),ksi_av(i),Ks(i),
!    la,sen,lat,cons,mus(i),Ksnow_interm(i),Kinterm_ice(i),Kice_bott(i),
!    T0_star(j-1),Ts_star(1),Tsi_star(1),Ti_star(1),Tfr,F_it,mui,Ts(i-1),
!    Tsi(i-1),Ti(i-1));
    
! CASE 1: 3 LAYERS (SNOW + INTERM LAYER + SEA ICE)
!subroutine temperature_case1_snow_meteoric_ice(T0_iter,Ts_iter,Tsi_iter, &
!    Ti_iter,Tmix_star,T0_star,Ti_star,vis_fr,Fs,ks_snow,hs_old, &
!    hs_prec_bucket,kmi_av,hmi_old,hi_old,ksi_10_av,ksi_av,Ks,la,sen,lat, &
!    cons,mus,Ksnow_interm,mumi,Kinterm_ice,Kice_bott,T0_star_old,Ts_star_1, &
!    Tsi_star_1,Ti_star_1,Tfr,F_it,mui,Ts_old,Tsi_old,Ti_old)

subroutine temperature_case1_snow_meteoric_ice(T0_iter,Ts_iter,Tsi_iter, &
    Ti_iter,Tmix_star,T0_star,Ti_star,Fs,ks_snow,hs_old, &
    hs_prec_bucket,kmi_av,hmi_old,hi_old,ksi_10_av,ksi_av,Ks,la,sen,lat, &
    cons,mus,Ksnow_interm,mumi,Kinterm_ice,Kice_bott,T0_star_old,Ts_star_1, &
    Tsi_star_1,Ti_star_1,F_it,mui,Ts_old,Tsi_old,Ti_old)
 
use ice_params
use ice_global
   
implicit none
double precision T0_iter,Ts_iter,Tsi_iter,Ti_iter,Tmix_star,T0_star,Ti_star,Fs
double precision ks_snow,hs_old,hs_prec_bucket,kmi_av,hmi_old,hi_old,ksi_10_av,ksi_av
double precision Ks,la,sen,lat,cons,mus,Ksnow_interm,mumi,Kinterm_ice,Kice_bott,T0_star_old
double precision Ts_star_1,Tsi_star_1,Ti_star_1,F_it,mui,Ts_old,Tsi_old,Ti_old
double precision I0_it,IM_it,Ts_star,Tsi_star

double precision ISI_it

I0_it=vis_fr*Fs*exp(-ks_snow*(hs_old+hs_prec_bucket));
IM_it=I0_it*exp(-kmi_av*hmi_old);

call temperature_fraction_of_solar_radiation_penetration(ksi_10_av,ksi_av, &
    hi_old,IM_it,ISI_it);

T0_star=((1-Ks/(la+sen+lat+cons))*T0_star_old)+((Ks/(la+sen+lat+cons)) &
    *Ts_star_1)+(0*Tsi_star_1)+(0*Ti_star_1)+((-F_it-I0_it)/(la+sen+lat+cons));
Ts_star=((mus*Ks)*T0_star_old)+((1-mus*(Ksnow_interm+Ks))*Ts_star_1)+((mus &
    *Ksnow_interm)*Tsi_star_1)+(0*Ti_star_1)+(I0_it*mus);
Tsi_star=(0*T0_star_old)+((mumi*Ksnow_interm)*Ts_star_1)+((1- mumi &
    *(Ksnow_interm+Kinterm_ice))*Tsi_star_1)+((mumi*Kinterm_ice)*Ti_star_1) &
    +(IM_it*mumi);
Ti_star=(0*T0_star_old)+(0*Ts_star_1)+((mui*Kinterm_ice)*Tsi_star_1)+((1-mui &
    *(Kinterm_ice+Kice_bott))*Ti_star_1)+(mui*Kice_bott*Tfr + mui*ISI_it);

T0_iter=T0_star;
Ts_iter=Ts_star;
Tsi_iter=Tsi_star;
Ti_iter=Ti_star;
Tmix_star=Tfr;

if( Ts_iter>Tfr ) then
    Ts_iter=Tfr;
end if
if( Tsi_iter>Tfr ) then
     Tsi_iter=Tfr;
     Ts_iter=Tfr;
end if
if( Ti_iter>Tfr ) then
    Ti_iter=Tfr;
    Tsi_iter=Tfr;
    Ts_iter=Tfr;
end if

if( T0_iter>Tfr ) then
    T0_iter=Tfr;
    Ts_iter=((mus*Ks)*Tfr)+((1-mus*(Ksnow_interm+Ks))*Ts_old)+((mus &
        *Ksnow_interm)*Tsi_old)+(0*Ti_old)+(I0_it*mus);
    if( Ts_iter<Tfr ) then
        Tsi_iter=(0*Tfr)+((mumi*Ksnow_interm)*Ts_old)+((1- mumi*(Ksnow_interm &
            +Kinterm_ice))*Tsi_old)+((mumi*Kinterm_ice)*Ti_old)+(IM_it*mumi);
        if( Tsi_iter<Tfr ) then
            Ti_iter=(0*Tfr)+(0*Ts_old)+((mui*Kinterm_ice)*Tsi_old)+((1-mui &
                *(Kinterm_ice+Kice_bott))*Ti_old)+(mui*Kice_bott*Tfr + mui &
                *ISI_it);
            if( Ti_iter>Tfr ) then
                Ti_iter=Tfr;
            end if
        else
            Tsi_iter=Tfr;
            Ti_iter=(mui*Kinterm_ice)*Tfr + (1-mui*(Kinterm_ice+Kice_bott)) &
                *Ti_old + (mui*Kice_bott*Tfr);
            if( Ti_iter>Tfr ) then
                Ti_iter=Tfr;
            end if
        end if
    else
        Ts_iter=Tfr;
        Tsi_iter=(0*Tfr)+((mumi*Ksnow_interm)*Tfr)+((1- mumi*(Ksnow_interm &
            +Kinterm_ice))*Tsi_old)+((mumi*Kinterm_ice)*Ti_old)+(IM_it*mumi);
        if( Tsi_iter<Tfr ) then
            Ti_iter=(0*Tfr)+(0*Tfr)+((mui*Kinterm_ice)*Tsi_old)+((1-mui &
                *(Kinterm_ice+Kice_bott))*Ti_old)+(mui*Kice_bott*Tfr + mui &
                *ISI_it);
            if( Ti_iter>Tfr ) then
                Ti_iter=Tfr;
            end if
        else
            Tsi_iter=Tfr;
            Ti_iter=(mui*Kinterm_ice)*Tfr + (1-mui*(Kinterm_ice+Kice_bott)) &
                *Ti_old + (mui*Kice_bott*Tfr + mui*ISI_it);
            if( Ti_iter>Tfr ) then
                Ti_iter=Tfr;
            end if
        end if
    end if
end if

end
! INITIATE FUNCTION WITH THIS:
![T0_iter,Ts_iter,Tsi_iter,Ti_iter,Tmix_star,T0_star(j),Ti_star(j)] = ...
!    temperature_case2_snow_ice(vis_fr,Fs(i),ks_snow(i),hs(i-1),
!    hs_prec_bucket(i),hi(i-1),ksi_10_av(i),ksi_av(i),Ks(i),la,sen,lat,cons,
!    mus(i),Ksnow_ice(i),Kice_bott(i),F_it,mui,Tfr,T0_star(j-1),Ts_star(1),
!    Ti_star(1),Ts(i-1),Ti(i-1));

! CASE 2: 2 LAYERS (SNOW + SEA ICE)
!subroutine temperature_case2_snow_ice(T0_iter,Ts_iter,Tsi_iter,Ti_iter, &
!    Tmix_star,T0_star,Ti_star,vis_fr,Fs,ks_snow,hs_old,hs_prec_bucket,hi_old, &
!    ksi_10_av,ksi_av,Ks,la,sen,lat,cons,mus,Ksnow_ice,Kice_bott,F_it,mui, &
!    Tfr,T0_star_old,Ts_star_1,Ti_star_1,Ts_old,Ti_old)

subroutine temperature_case2_snow_ice(T0_iter,Ts_iter,Tsi_iter,Ti_iter, &
    Tmix_star,T0_star,Ti_star,Fs,ks_snow,hs_old,hs_prec_bucket,hi_old, &
    ksi_10_av,ksi_av,Ks,la,sen,lat,cons,mus,Ksnow_ice,Kice_bott,F_it,mui, &
    T0_star_old,Ts_star_1,Ti_star_1,Ts_old,Ti_old)

use ice_params
use ice_global

implicit none
double precision T0_iter,Ts_iter,Tsi_iter,Ti_iter,Tmix_star,T0_star,Ti_star,Fs
double precision ks_snow,hs_old,hs_prec_bucket,hi_old,ksi_10_av,ksi_av,Ks,la,sen,lat
double precision cons,mus,Ksnow_ice,Kice_bott,F_it,mui,T0_star_old,Ts_star_1,Ti_star_1
double precision Ts_old,Ti_old,I0_it,IM_it,Ts_star

double precision ISI_it

I0_it=vis_fr*Fs*exp(-ks_snow*(hs_old+hs_prec_bucket));
IM_it=I0_it;

call temperature_fraction_of_solar_radiation_penetration(ISI_it,ksi_10_av, &
    ksi_av,hi_old,IM_it);

T0_star=((1-Ks/(la+sen+lat+cons))*T0_star_old)+((Ks/(la+sen+lat+cons)) &
    *Ts_star_1)+(0*Ti_star_1)+((-F_it-I0_it)/(la+sen+lat+cons));
Ts_star=((mus*Ks)*T0_star_old)+((1- mus*(Ksnow_ice+Ks))*Ts_star_1)+((mus &
    *Ksnow_ice)*Ti_star_1)+(I0_it*mus);
Ti_star=(0*T0_star_old)+((mui*Ksnow_ice)*Ts_star_1)+((1 - mui*(Ksnow_ice &
    + Kice_bott))*Ti_star_1)+(mui*Kice_bott*Tfr + mui*ISI_it);

T0_iter=T0_star;
Ts_iter=Ts_star;
Ti_iter=Ti_star;
Tmix_star=Tfr;

if( Ts_iter>Tfr ) then
    Ts_iter=Tfr;
end if
if( Ti_iter>Tfr ) then
    Ti_iter=Tfr;
    Ts_iter=Tfr;
end if

if( T0_iter>Tfr ) then
    T0_iter=Tfr;
    Ts_iter=((mus*Ks)*Tfr)+((1- mus*(Ksnow_ice+Ks))*Ts_old)+((mus*Ksnow_ice) &
        *Ti_old)+(I0_it*mus);
    if( Ts_iter<Tfr ) then
        Ti_iter=(0*Tfr)+((mui*Ksnow_ice)*Ts_old)+((1 - mui*(Ksnow_ice &
            + Kice_bott))*Ti_old)+(mui*Kice_bott*Tfr + mui*ISI_it);
    else
        Ts_iter=Tfr;
        Ti_iter=(mui*Ksnow_ice)*Tfr + (1 - mui*(Ksnow_ice + Kice_bott)) &
            *Ti_old + (mui*Kice_bott*Tfr + mui*ISI_it);
        if( Ti_iter>Tfr ) then
            Ti_iter=Tfr;
        end if
    end if
end if
Tsi_iter=Tfr;

end
! INITIATE FUNCTION WITH THIS:
![T0_iter,Ts_iter,Tsi_iter,Ti_iter,Tmix_star,T0_star(j),Ti_star(j)] = ...
!    temperature_case3_meteoric_ice(vis_fr,Fs(i),kmi_av(i),hmi(i-1),hi(i-1),
!    ksi_10_av(i),ksi_av(i),Kmi(i),la,sen,lat,conmi,mumi(i),Kinterm_ice(i),
!    Kice_bott(i),Tfr,mui,T0_star(j-1),Tsi_star(1),F_it,Tsi(i-1),Ti_star(1),
!    Ti(i-1));

! CASE 3: 2 LAYERS (INTERM LAYER + SEA ICE)
!subroutine temperature_case3_meteoric_ice(T0_iter,Ts_iter,Tsi_iter,Ti_iter, &
!    Tmix_star,T0_star,Ti_star,vis_fr,Fs,kmi_av,hmi_old,hi_old,ksi_10_av, &
!    ksi_av,Kmi,la,sen,lat,conmi,mumi,Kinterm_ice,Kice_bott,Tfr,mui, &
!    T0_star_old,Tsi_star_1,F_it,Tsi_old,Ti_star_1,Ti_old)

subroutine temperature_case3_meteoric_ice(T0_iter,Ts_iter,Tsi_iter,Ti_iter, &
    Tmix_star,T0_star,Ti_star,Fs,kmi_av,hmi_old,hi_old,ksi_10_av, &
    ksi_av,la,sen,lat,conmi,mumi,Kinterm_ice,Kice_bott,mui, &
    T0_star_old,Tsi_star_1,F_it,Tsi_old,Ti_star_1,Ti_old,Kmi)

use ice_params
use ice_global

implicit none
double precision T0_iter,Ts_iter,Tsi_iter,Ti_iter,Tmix_star,T0_star,Ti_star,Fs
double precision kmi_av,hmi_old,hi_old,ksi_10_av,ksi_av,la,sen,lat,conmi,mumi
double precision Kinterm_ice,Kice_bott,mui,T0_star_old,Tsi_star_1,F_it,Tsi_old
double precision Ti_star_1,Ti_old,I0_it,IM_it,Tsi_star,Kmi

double precision ISI_it

I0_it=vis_fr*Fs*exp(-kmi_av*hmi_old);
IM_it=I0_it;

call temperature_fraction_of_solar_radiation_penetration(ISI_it,ksi_10_av, &
    ksi_av,hi_old,IM_it);

T0_star=((1-Kmi/(la+sen+lat+conmi))*T0_star_old)+((Kmi/(la+sen+lat+conmi)) &
    *Tsi_star_1)+(0*Ti_star_1)+((-F_it-I0_it)/(la+sen+lat+conmi));
Tsi_star=((mumi*Kmi)*T0_star_old)+((1 - mumi*(Kinterm_ice +Kmi))*Tsi_star_1) &
    +((mumi*Kinterm_ice)*Ti_star_1)+(I0_it*mumi);
Ti_star=(0*T0_star_old)+((mui*Kinterm_ice)*Tsi_star_1)+((1 - mui*(Kinterm_ice &
    + Kice_bott))*Ti_star_1)+(mui*Kice_bott*Tfr + mui*ISI_it);

T0_iter=T0_star;
Tsi_iter=Tsi_star;
Ti_iter=Ti_star;
Tmix_star=Tfr;

if( Tsi_iter>Tfr ) then
    Tsi_iter=Tfr;
end if
if( Ti_iter>Tfr ) then
    Ti_iter=Tfr;
    Tsi_iter=Tfr;
end if

if( T0_iter>Tfr ) then
    T0_iter=Tfr;
    Tsi_iter=((mumi*Kmi)*Tfr)+((1 - mumi*(Kinterm_ice +Kmi))*Tsi_old)+((mumi &
        *Kinterm_ice)*Ti_old)+(I0_it*mumi);
    if( Tsi_iter<Tfr ) then
        Ti_iter=(0*Tfr)+((mui*Kinterm_ice)*Tsi_old)+((1 - mui*(Kinterm_ice &
            + Kice_bott))*Ti_old)+(mui*Kice_bott*Tfr + mui*ISI_it);
    else
        Tsi_iter=Tfr;
        Ti_iter=(mui*Kinterm_ice)*Tfr + (1 - mui*(Kinterm_ice + Kice_bott)) &
            *Ti_old + (mui*Kice_bott*Tfr + mui*ISI_it);
        if( Ti_iter>Tfr ) then     ! ::RASA: changed from Ti > Tfr
            Ti_iter=Tfr;    ! ::RASA: changed from Ti = Tfr; 
        end if
    end if
end if
Ts_iter=Tfr;

end
! INITIATE FUNCTION WITH THIS:
![T0_iter,Ts_iter,Tsi_iter,Ti_iter,Tmix_star,T0_star(j),Ti_star(j)] = ...
!    temperature_case4_ice(hi(i-1),vis_fr,Fs(i),ksi_10_av(i),ksi_av(i),
!    Kice_surf(i),la,sen,lat,coni,mui,Kice_bott(i),Tfr,T0_star(j-1),Ti_star(1),
!    Ti(i-1),F_it);

! CASE 4: ONLY SEA ICE
subroutine temperature_case4_ice(T0_iter,Ts_iter,Tsi_iter,Ti_iter,Tmix_star, &
    T0_star,Ti_star,hi_old,Fs,ksi_10_av,ksi_av,Kice_surf,la,sen,lat, &
    coni,mui,Kice_bott,T0_star_old,Ti_star_1,Ti_old,F_it)

use ice_params
use ice_global

implicit none
double precision T0_iter,Ts_iter,Tsi_iter,Ti_iter,Tmix_star,T0_star,Ti_star,hi_old
double precision Fs,ksi_10_av,ksi_av,Kice_surf,la,sen,lat,coni,mui,Kice_bott
double precision T0_star_old,Ti_star_1,Ti_old,F_it,I0_it,ISI_10_it,ISI_it

double precision, parameter :: eps=0.1D+0

if( hi_old>=eps ) then
    ISI_10_it=vis_fr*Fs*exp(-ksi_10_av*eps);
    ISI_it=ISI_10_it*exp(-ksi_av*(hi_old-eps));
else
    ISI_10_it=vis_fr*Fs*exp(-ksi_10_av*hi_old);
    ISI_it=ISI_10_it;
end if

I0_it=ISI_10_it;

T0_star=((1-Kice_surf/(la+sen+lat+coni))*T0_star_old)+((Kice_surf/(la+sen+lat &
    +coni))*Ti_star_1)+0+((-F_it-I0_it)/(la+sen+lat+coni))
Ti_star=((mui*Kice_surf)*T0_star_old)+((1 - mui*(Kice_surf + Kice_bott)) &
    *Ti_star_1)+0+(mui*Kice_bott*Tfr + mui*ISI_it)

T0_iter=T0_star;
Ti_iter=Ti_star;
Tmix_star=Tfr;

if( Ti_iter>Tfr ) then
    Ti_iter=Tfr;
end if

if( T0_iter>Tfr ) then
    T0_iter=Tfr;
    Ti_iter=(mui*Kice_surf)*Tfr + (1 - mui*(Kice_surf + Kice_bott))*Ti_old + &
        (mui*Kice_bott*Tfr + mui*ISI_it);
    if( Ti_iter>Tfr ) then
        Ti_iter=Tfr;
    end if
end if

Ts_iter=Tfr;
Tsi_iter=Tfr;
!write(6,*) 'T0_iter',T0_iter
end
! INITIATE FUNCTION WITH THIS:
![Tmix_star,T0_star(j),Ti_star(j),T0_iter,Ts_iter,Tsi_iter,Ti_iter] = ...
!    temperature_case5_frazil(Fs(i),infra_fr,k_ocean_red,h_mix,k_ocean_vis,
!    Tmix(i-1),F_it,R(i-1),row,cpw,deltat,Tfr);

! CASE 5: 0 LAYER OR FRAZIL ICE FORMATION
subroutine temperature_case5_frazil(Tmix_star,T0_star,Ti_star,T0_iter, &
    Ts_iter,Tsi_iter,Ti_iter,Fs,Tmix_old,F_it,R_old,deltat)

use ice_params
use ice_global

implicit none
double precision Tmix_star,T0_star,Ti_star,T0_iter,Ts_iter,Tsi_iter,Ti_iter,Fs
double precision Tmix_old,F_it,R_old,deltat
    
double precision I0_it,Ts_star,Tsi_star

I0_it=Fs*((infra_fr*exp(-k_ocean_red*h_mix)) + &
	((1-infra_fr)*exp(-k_ocean_vis*h_mix)));
Tmix_star=Tmix_old-((F_it + I0_it + R_old)/(h_mix*row*cpw))*deltat;
T0_star=Tmix_star;
Ts_star=Tfr;
Tsi_star=Tfr;
Ti_star=Tfr;
T0_iter=T0_star;
Ts_iter=Ts_star;
Tsi_iter=Tsi_star;
Ti_iter=Ti_star;

end
! INITIATE FUNCTION WITH THIS:
![ISI_it] = temperature_fraction_of_solar_radiation_penetration(ksi_10_av,ksi_av,hi,IM_it);
! OR (only for case4: sea ice):
![ISI_it] = temperature_fraction_of_solar_radiation_penetration(ksi_10_av,ksi_av,hi,vis_fr*Fs);

!FRACTION OF SOLAR RADIATION PENETRATING SEA ICE
subroutine temperature_fraction_of_solar_radiation_penetration(ISI_it,ksi_10_av,ksi_av,hi,IM)

use ice_params

implicit none
double precision ISI_10_it,ISI_it,ksi_10_av,ksi_av,hi,IM

double precision, parameter :: eps=0.1D+0

if( hi>=eps ) then
    ISI_10_it=IM*exp(-ksi_10_av*eps);
    ISI_it=ISI_10_it*exp(-ksi_av*(hi-eps));
else if( hi<eps ) then
    ISI_10_it=IM*exp(-ksi_10_av*hi);
    ISI_it=ISI_10_it;
end if

end
! INITIATE FUNCTION WITH THIS:
![Tice(i)] = temperature_of_ice(hmi(i-1),hmi_min,kice_surf(i),Ti(i),hi(i-1),
! kmi,Tsi(i),hs(i-1),hs_prec_bucket(i),hs_min,ks_av(i),Ts(i),hi_min,T0(i),Tfr)
    
! SEA ICE TEMPERATURE
subroutine temperature_of_ice(Tice,hmi_old,kice_surf,Ti,hi_old, &
    Tsi,hs_old,hs_prec_bucket,ks_av,Ts,T0)

use ice_params
use ice_global

implicit none
double precision Tice,hmi_old,kice_surf,Ti,hi_old,Tsi,hs_old,hs_prec_bucket
double precision Ts,T0,ks_av

if( hmi_old>hmi_min ) then
    
    Tice=((kice_surf*Ti/(hi_old/2) + kkmi*Tsi/(hmi_old/2))/(kice_surf/(hi_old &
        /2) +kkmi/(hmi_old/2)))-Tfrs;

else if( hmi_old<=hmi_min) then
    
    if( (hs_old+hs_prec_bucket)>hs_min ) then
        
        Tice=((kice_surf*Ti/(hi_old/2) + ks_av*Ts/(hs_old/2))/(kice_surf &
            /(hi_old/2) + ks_av/(hs_old/2)))-Tfrs;
    
    else if( (hs_old+hs_prec_bucket)<=hs_min ) then
        
        if( hi_old>hi_min ) then
            Tice=T0-Tfrs_appr;
        else if( hi_old<=hi_min ) then
            Tice=Tfr-Tfrs;
        end if
        
    end if
end if
    
end
