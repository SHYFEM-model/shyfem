
!--------------------------------------------------------------------------
!
!    Copyright (C) 2004,2006,2011  Petras Zemlys
!    Copyright (C) 2015-2019  Georg Umgiesser

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

! aquabc_II to fem interface
!
! Aquabc version 3.0
!  - 2019 02 - State variables indices are introduced for water column
!
! Produced by Petras 2014
!
! Based on:
!    $Id: bio3d.f,v 1.33 2008-10-10 09:29:54 georg Exp $
!

! aquabc_fem_interface is biogeochemical model specific. Please change line 11 and
! corresponding lines in Make file to switch to different biogeochemical models of AQUABC series

! Contents of the rest:
!   subroutine aquabc_fem_interface - Interface between SHYFEM and AQUABC (Aquatic Biochemical Cycling)
!   subroutine init_binary_output_wc  - initialisation of writing to binary file for WC
!   subroutine init_binary_output_bs  - initialisation of writing to binary file for WC
!   subroutine write_binary_output_wc - writing to binary file for WC
!   subroutine write_binary_output_bs - writing to binary file for BC
!   subroutine setload            - reads point source loads. Not fully tested
!   subroutine check_var          - checking strange values.
!   subroutine check2Dbio         - tests array for NaN and strange values for EUTRO variables. Used by check_var
!   function max_int              - finds maximum element in one dimensional array
!   subroutine inicfil_aquabc     - reads initial values for state variables of WC (also for repeated runs)
!   subroutine inicfils_aquabc    - reads initial values for state variables of BS (also for repeated runs)
!   subroutine calc_time          - prints start and end time of the program and calculates run time
!   subroutine calc_run_time      - calculates run time
!   subroutine print_real_time    - prints real date and time of program start
!   subroutine get_ITOT           - reads light for aquabc. Not used more. Light comes from hydrodynamic module
!   subroutine biotser_init       - initialisation of ASCII output for given nodes (stations)
!   subroutine biotser_write      - writing of ASCII output for given nodes
!   subroutine cur_param_read_ini_wc - intialisation of parameters array for WC model
!   subroutine cur_param_read_ini_bs - intialisation of parameters array for BS model
!   subroutine cur_param_read     - reading WC model parameters(constants)
!   subroutine aquabcini         - does preparations for WC model
!   FUNCTION ALIGHT_ALUKAS(DEPTH,CHLA) - calculates light intensity for bottom of the layer
!   routines to read forcing time series for water temp, pH and interpolation of their values. Not used more

! Note:
!   To correct file names given in *.str section 'name' correct following subroutine in file subsys.f:
!     subroutine fnm_aquabc_init
!********************************************************************
!********************************************************************
!********************************************************************
!
! $Id: aquabc_II_fem_interface.f
!
! Interface between SHYFEM and AQUABC_II (Aquatic Biochemical Cycling)
! Adapted by Petras Zemlys and Ali Erturk from an interface between
! EUTRO and SHYFEM called BIO3D
!
! revision log :
!
! 20.07.2004	pzy	modified to allow dynamic forcings from external files
! 01.07.2006	pzy	EUTRO changed by new module AQUABC with 14 state variables
! 01.07.2011	pzy	Water column eutrophication module changed to ALUKAS_II
! 26.02.2015	ggu	changed VERS_7_1_6
! 05.06.2015	ggu	changed VERS_7_1_12
! 17.07.2015	ggu	changed VERS_7_1_52
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 15.04.2016	ggu	changed VERS_7_5_8
! 14.11.2017	ggu	changed VERS_7_5_36
! 03.04.2018	ggu	changed VERS_7_5_43
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 29.03.2022	ggu	icall was not initialized
!
! Notes :
!
! ALUKAS = Advanced Level nUtrient Kinetics for Aquatic Systems
!
!                  ccccccccccccccccccccccccc
!                  c    WC state variables c
!                  ccccccccccccccccccccccccc
!NO   NAME IN MODEL                   STATE VARIABLE                       UNITS
!--   --------------  ----------------------------------------             -------
! 1  NH4_N            Ammonium nitrogen                                    mg N/l
! 2  NO3_N            Nitrate nitrogen                                     mg N/l
! 3  PO4_P            Orthophosph. phosph. total in water and in particl.  mg P/l
! 4  DISS_OXYGEN      Dissolved oxygen                                     mg O2/l
! 5  DIA_C            Diatoms carbon                                       mg C/l
! 6  ZOO_C            Zooplankton  carbon                                  mg C/l
! 7  ZOO_N           Zooplankton nitrogen                                  mg C/l derived variable, switched to constant stoichiometry
! 8  ZOO_P           Zooplankton phosphorus                                mg C/l derived variable, switched to constant stoichiometry
! 9  DET_PART_ORG_C   Detritus particulate org. carbon                     mg C/l
!10  DET_PART_ORG_N   Detritus particulate org. nitrogen                   mg N/l
!11  DET_PART_ORG_P   Detritus particulate org. phosphorus                 mg P/l
!12  DISS_ORG_C       Dissolved organic carbon                             mg C/l
!13  DISS_ORG_N       Dissolved organic nitrogen                           mg N/l
!14  DISS_ORG_P       Dissolved organic phosphorus                         mg P/l
!15  CYAN_C           Cyanobacteria carbon                                 mg C/l
!16  OPA_C            Other phytoplankton carbon                           mg C/l
!17  DISS_Si          Dissoloved silica                                    mg Si/l
!18  BIOG_Si          Biogenic silica                                      mg Si/l
!19  FIX_CYN_C        Nitrogen fixing cyanobacteria carbon                 mg C/l
!20  INORG_C          Dissolved inorganic carbon                           mol/m3
!21  TOT_ALK          Total alkalinity                                     mol/m3
!22  FE_II            Ferous total (water and particles)                   mg/l
!23  FE_III           Feric total                                          mg/l
!24  MN_II            Manganous total                                      mg/l
!25  MN_IV            Manganese 4+ total                                   mg/l
!26  CA               Calcium total                                        mg/l
!27  MG               Magnesium total                                      mg/l
!28  S_PLUS_6         Sulphate sulphur total                               mg/l
!29  S_MINUS_2        Sulphide sulphur total                               mg/l
!30  CH4_C            Methane carbon total                                 mg/l


!
!          DERIVED VARIABLES SENT TO OUTPUT:
!31   Total nitrogen
!32   Total phosphorus
!33   Total DIN
!34   pH
!35   Tot_phyt
!36   Tot_cyan
!37   Nutr_PO4_P_Diss      
!38   Nutr_PO4_P_Part      
!39   Met_Fe2_Diss         
!40   Met_Fe2_Part         
!41   Met_Fe3_Diss         
!42   Met_Fe3_Part         
!43   Met_Mn2_Diss         
!44   Met_Mn2_Part         
!45   Met_Mn4_Diss         
!46   Met_Mn4_Part         
!47   Met_Calcium_Diss     
!48   Met_Calcium_Part     
!49   Met_Magnesium_Diss   
!50   Met_Magnesium_Part
!51   Total_dissolved_Fe  
!52   Total_particulate_Fe
!53   Total_Fe            
!54   Nutr_Oxygen_Sat     
!  
!
!                  ccccccccccccccccccccccccc
!                  c    BS state variables c
!                  ccccccccccccccccccccccccc
!
!
!           All concentrations are in solutes plus solids
!           unless other sediment fractions are pointed out
!
!        NO  NAME IN MODEL  STATE VARIABLE                           UNITS
!        --  ---------      -------------                            -----
!        1     SED_NH4N     BS total amonia                          mg/l
!        2     SED_NO3N     BS nitrates                              mg/l
!        3     SED_DON      BS dissolved organic nitrogen            mg/l
!        4     SED_PON      BS particulate organic nitrogen          mg/l
!        5     SED_PO4P     BS dissolved ortophosphates phosphorus   mg/l
!        6     SED_DOP      BS dissolved organic phophorus           mg/l
!        7     SED_POP      BS particulate organic phosphorus        mg/l
!        8     SED_DOXY     BS dissolved oxygen                      mg/l
!        9     SED_DOC      BS dissolved organic carbon              mg/l
!       10     SED_POC      BS particulate organic carbon            mg/l
!       11     SED_DSi      BS dissolved silica                      mg/l
!       12     SED_PSi      BS particulate silica                    mg/l
!       13     SED_INORG_C  BS diss. inorg. carbon in pore water     mol/m3
!       14     SED_TOT_ALK  BS total alkalinty in pore water         mol/m3
!       15     SED_SALT     BS salinity in pore water                promiles
!       16     FE_II        BS ferous                                mg/l
!       17     FE_III       BS feric                                 mg/l
!       18     MN_II        BS manganous                             mg/l
!       19     MN_IV        BS manganese 4+                          mg/l
!       20     CA           BS calcium                               mg/l
!       21     MG           BS magnesium                             mg/l
!       22     S_PLUS_6     BS sulphate sulphur                      mg/l
!       23     S_MINUS_2    BS sulphide sulphur                      mg/l
!       24     CH4_C        BS methane carbon                        mg/l
!
!  DERIVED VARIABLES SENT TO OUTPUT:
!       25     NH4N in pore water
!       26     PO4P in pore water
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Running:
!  ibio = 0 - ecological module not used                                                                                      c
!  ibio = 1 - WC kinetics without BS kinetics. Settling and                                  c
!             deposition switched on. Settling rates and deposition fraction                 c
!             for each type of sediments are still hardcoded in subroutine aquabc.           c
!  ibio = 2 - WC and BS kinetics. Sediment properties still hardcoded in subroutine aquabc.  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Binary Ooutput:                                 c
!  *.bwc  - water column variables     (unit iub) c
!  *.bbs  - bottom sediments variables (unit iubs)ccccccccccccccccccccccc
!  ASCII output for state  and intermediate variables in observations   c
!  stations (nodes) can be used.                                        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Input:                                                c
! *.str name section aquabc specific files:             c
!        bioin       - initial state for WC             c
!        biosin      - initial state for BS             c
!        bbs_lev   - BS vertical levels                 c
!        biocon    - constants for WC                   c
!        bioscon   - constants for BS                   c
!        bioaow    - ASCII output control file for WC   c
!        bioaos    - ASCII output control file for BC   c
!
!   word 'fast' in variable names means prepared for vectorization        
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!***********************************************
!***********************************************
!***********************************************

        subroutine ecological_module

! general interface to ecological module

        implicit none

        call aquabc_II_fem_interface

        end

!***********************************************
!***********************************************
!***********************************************

      subroutine aquabc_II_fem_interface

      !instead of param.h        
      use levels, only: ilhkv, nlvdi, nlv  
      use basin, only: nkn, nkndi, neldi, ipv
      
      use mod_debug
      use mod_diff_visc_fric
      !use mod_conz

      use aquabc_II_vars   !instead of aquabc_II_h
      use aquabc_II_sed_ini

      implicit none

      include 'aquabc_II_aout.h' !variables and parameters for output to ASCII file

      include 'femtime.h'
      

      !integer it	!time in seconds. Comes from femtime.h now
      !integer idt  !time step in seconds.Comes from femtime.h now
      double precision tdouble
      real dt      !time step in seconds!495vers real dt remove idt
      real dtday    !time step in days
      real tday, t0, tsec
!-------------------------------------------------------
      ! WC, BS state variables related arrays
      ! nstate - no of WC state vars, defined in mod aquabc_II_vars
      ! nkndi - max no of nodes,     defined in mod basin
      ! nlvdi - max no of layers,    defined in mod levels

      real, save, allocatable :: e(:,:,:)	      ! WC state vector
      real, save, allocatable :: eload(:,:,:)     ! WC vector of loadings


!     Variables for calculation of material fluxes trough sections given in *.str
      integer, save :: iflux              ! indicator to calculate fluxes (0,1)
      integer, parameter ::  nflux = 2    ! actual number of state variables for flux calculation
                                          ! It should be equal to parameter iconz in para section?
      real, save, allocatable :: eflux(:,:,:) ! WC state variables selected for for flux calculation.
      !character(len=4) :: ext		!extension of fluxes output file (e.g., '.fxa')  Passed as constant now                                        
                                          
      real, save, allocatable :: es(:,:,:)  ! BS state variables


!----------------------------------------------------

      !loops related:
      integer k,i,l
      integer j
      integer ls
      integer lmax        !max no of layers for the node
      integer nvar

!     aquabc parameters(constants)
      real par(nconst)     !aquabc model WC parameters array
      save par
      real par_sed(nsconst)    !aquabc model BS parameters array
      save par_sed

!     Arrays passed to aquabc for faster calculations by vectorization
      real, save, allocatable :: e_fast0(:,:,:) ! initial state vector arranged
                                                !  for faster calculations by vectorization
      real, save, allocatable :: e_fast1(:,:,:) ! final state vector arranged
                                                !  for faster calculations by vectorization
      integer, save, allocatable :: lmax_fast   (:)   ! number of levels for node
      real, save, allocatable :: depth_fast  (:,:) ! depth
      real, save, allocatable :: vel_fast    (:,:) ! current velocity
      real, save, allocatable :: wtemp_fast  (:,:) ! water temperature
      real, save, allocatable :: wind_fast   (:)   ! wind speed
      real, save, allocatable :: atemp_fast  (:)   ! air temperature
      real, save, allocatable :: ice_cover_fast(:) ! ice cover (node area fraction)
      real, save, allocatable :: sal_fast    (:,:) ! water salinity
      real, save, allocatable :: light_fast  (:)   ! incident light for the layer
      real, save, allocatable :: vol_fast    (:,:) ! reactor volume
      real, save, allocatable :: vol_old_fast(:,:) !
      real, save, allocatable :: area_fast   (:,:) ! reactor area
      real, save, allocatable :: waves_fast  (:,:) ! wind waves: 1-significant  wave height
                                                !             2 - mean wave period

      double precision, save, allocatable :: es_fast0(:,:,:)  !sediment state vars to pass to kinetics
      double precision, save, allocatable :: es_fast1(:,:,:)  !sediment state vars to get from kinetics

!----------------------------------------------
!----------------------------------------------
      integer nlayers  ! maximum number of levels in WC (calculated)

      integer, save, allocatable :: ilhkv_sed(:) ! Needed for compatibility of WC and BS writing in biotser_write          
      !real, save, allocatable :: hev_sed(:)      ! For binary output    

      real elaux(nstate) ! WC state vars  loads past to aquabc

      ! controls calls to routines
      character*10 what,whataux
      character*2 whatn

      integer ibio  ! aquabc run options (0 - no aquabc, 1 - WC only, 2 - WC and BS)
      integer isedi ! sediment transport module use indicator
      integer iwave ! wave parametric module use indicator. iwave = 1, if isedi = 1

      !integer id     	! background number for variable numbers
      !integer id_file   ! file id. Is defined by Shyfem utilities
      integer mode
      integer nintp

      !integer ivar

      ! temp, sal and velocities in reactor
      real t,s
      real u,v

!     arrays for initialisation of statvars with data statement
      real einit(nstate)       ! for initial values of WC variables
      save einit
      real esinit(nsstate)     ! for initial values of BS variables 1-st layer
      save esinit
      real esinit2(nsstate)    ! for initial values of WC variables 2-st layer
      real esinit3(nsstate)    ! for initial values of WC variables 3-st layer
      real elinit(nstate)
      save elinit
      real ebound(nstate)
      save ebound              !ebound is used in case no values are given in STR file

      integer, save, allocatable :: idbio(:)

      real tstot(nstate)   !for mass test
      real tsstot(nsstate)
      
!     first call indicator
      integer, save :: icall = 0
      
!     first time iteration indicator for aquabc      
      integer ifirst_it
      save ifirst_it

      integer iunit
!---------------------------------
      ! geometrical parameters of reactor
      real area, areanew, vol,vel
      real volold
      real depth, depthnew

!     parameters for transport/diffusion resolution
      real rkpar,difmol
      save rkpar,difmol

!     function names
      real getpar       !obtains parameters values
      integer iround
      integer nbnds

      integer nbc



      logical bsedim  !true -process BS
      save bsedim

      logical bcheck  ! to check state variables for strange values
                      ! not sure if it works properly. Fixme
                      ! better to to use debug_stranger in
                      ! internal routines for a while

      logical bresi,breact,bdecay

      integer ie,ii
      integer kspec


      ! node specific hydro-meteo variables
      real windspeed,tempair
      real rh,twetbulb,cloudc,pr !not used, just for copmatibility
      real ice_cover  ! reactor ice_cover fraction
      real wave_h     ! wave height
      real wave_mp    ! wave mean period
      real wave_pp    ! wave peak period
      real wave_d     ! wave mean direction

      real mass
      real wsink

      integer iespecial,inspecial
      save iespecial,inspecial

      integer iub
      save iub
      integer iubs
      save iubs
      integer, save :: ia_out(4)      
      !double precision, save :: da_out(4)
      

      logical has_output_d,next_output_d
      logical has_output,next_output      ! for old version

      integer idtc,itmc,itsmed
      save idtc,itmc

!     Arrays for timetable functions. Not used in this version
      real    PHTAB(1000), pH
      integer PHTIME(1000)
      save    PHTIME, PHTAB
      real    GET_PH

      real    TEMPTAB(1000),temp
      INTEGER TEMPTIME(1000)
      save    TEMPTIME, TEMPTAB
      real    GET_TEMP
!-----------------------------------

!     Solar radiation
      real ITOT             ! instantenous solar radiation, (now in W/m2 from hydrodynamics)
      real FDAY             ! equal to 1 in this version for compatibilty with older code


! Variables used for WC state variables output to ASCII files
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
       save    NBIOTS    !declared in aquabc_aout
       save    BIOTSFUN  !declared in aquabc_aout
       save    BIOTSNOD  !declared in aquabc_aout

! Variables used for BS state variables output to ASCII files
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
       save    NBIOTS_sed    !declared in aquabc_aout
       save    BIOTSFUN_sed  !declared in aquabc_aout
       save    BIOTSNOD_sed  !declared in aquabc_aout

!      variables used for intermediate WQ results(processes) output(processes) to ASCII files
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       integer dg_count	        !do loop parameter
       real, save, allocatable :: dg  (:,:,:,:) ! diagnostic data

       save    NDGTS   !declared in aquabc_aout
       save    DGTSFUN !declared in aquabc_aout
       save    DGTSNOD !declared in aquabc_aout

!      variables used for intermediate bottom sediments results(processes) output(processes)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       real, save, allocatable :: dg_sed  (:,:,:,:) ! diagnostic data

       save    NDGTS_sed   !declared in aquabc_aout
       save    DGTSFUN_sed !declared in aquabc_aout
       save    DGTSNOD_sed !declared in aquabc_aout


      integer ulogbio,ifileo
      save ulogbio

      integer ifirst_sed   !indicator of first call for BS, 1 - if first, 0 - if not
      data ifirst_sed /1/
      save ifirst_sed

      double precision dtime0, dtime

!     variables and functions for light management
      real ALIGHT_ALUKAS ! function for calculation light for the next layer

      character  sdate*8, stime*10  ! date and time of routine start
      save sdate, stime
      real start_time, end_time ! for calculation of processing time
      save start_time
      integer start_array(8), end_array(8) ! return data and time values for date_and_time
      save start_array

      integer STRANGERS   ! Function checks for NaNs and Inf
      integer idump       !indicator of dump for repeated runs
      
!----------------------------------------------------------------

!       print *, 'Éntering fem-aquabc interface'       
       
       idump = 1
!      idump = 0

       bresi = .false.
       breact = .true.
       bdecay = .false. 

       bcheck = .false.
!       bcheck = .true.       !.true. if state variables are checked
                              ! Check routines to use this option
       what = 'lagvebio'

! piece of code to be activated when numbers of nodes codes is needed only
!       do k=1,nkn
!        call  get_nodal_area_code(k,j)
!        i=ipv(k)
!        print *,'int_node:', k, 'ext_node:',i, 'area_code:',j
!       end do
!       stop

!------------------------------------------------------------------
!       initial and boundary conditions  [mg/l]			??
!       initial loadings                 [g/m**2/sec]		??
!------------------------------------------------------------------

!                   1     2    3      4      5      6      7
!       data einit /0.1, 0.1, 0.1, 10.0, 0.001, 0.001, 0.001,
!      *           0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
!                   8     9     10    11     12   13    14
!      *          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01/
!               15   16   17      18   19   20    21



!  Initialisation of BS state variables.  Sandy sed are assumed.! Fixme!
!ccccccccccccccccccccccccccccccccccccccccccccccc
!data esinit  /0.1,0.1,0.1,0.1,0.02,0.1,0.1,0.1,0.1,0.1,1.,0.1/ old initial conditions
!                   1      2      3      4      5      6      7
!       data esinit /0.5000,0.0319,1.000,1000.0,0.1800,0.1400,10.0,
!      *      0.0141,105.0,1400.0 ,3.500,10.000/
!              8       9     10      11     12
!                     1      2      3      4      5      6      7
!       data esinit2 /1.2000,0.0432,0.707 ,1500.,0.160,0.1200,130.0,
!      *      0.000001,90.0  ,1300.,8.000,130.0/
!              8       9     10      11     12
!                     1      2      3      4      5      6      7
!       data esinit3 /1.5000,0.0509,0.5000,1500.,0.150,0.4000,165.0,
!      *      0.000001,90.0  ,1000.,8.000,165.0/
!              8       9     10      11     12

!------------------------------------------------------------------

!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
! Initialization section (executed only first time step)
!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------

       if( icall .le. -1 ) return


!----------------------------------------------------------
       if( icall .eq. 0 ) then

          isedi = nint(getpar('isedi'))
          iwave = nint(getpar('iwave'))

          ibio = iround(getpar('ibio'))

          if(ibio .eq. 0 ) icall = -1
          if(ibio .gt. 2 .or. ibio .lt. 0) icall = -1
          if( icall .le. -1 ) return

          icall = 1
          ifirst_it = 1

          if(ibio .eq. 2) bsedim = .true.
          if(ibio .eq. 1) bsedim = .false.

!         Takes start date and time for evaluation of computation time
          call date_and_time(values = start_array)

          !iflux = 0
          iflux = 1
          
          
          
!---------------------------------------------
!         parameters for transport/diffusion resolution
!---------------------------------------------


          rkpar=getpar('chpar')
          difmol=getpar('difmol')

!-----------------------------------------------------
!          Allocation of arrays
!-----------------------------------------------------
          nbc = nbnds()
          allocate(idbio(nbc))

	      allocate(e    (nlvdi, nkndi,nstate))
	      allocate(eload(nlvdi, nkndi,nstate))
	      allocate(es   (noslay,nkndi,nsstate))
	      
	      allocate(eflux   (nlvdi, nkndi,  nflux)) ! WC state variables selected for for flux calculation
	      
   
!         Arrays passed to aquabc for faster calculations by vectorization                 
          allocate ( e_fast0     (nkndi, nstate, nlvdi))  ! initial state vector for vectorized calc.
          allocate ( e_fast1     (nkndi, nstate, nlvdi))  ! final state vector for vectorized calc.
          
          allocate ( lmax_fast   (nkndi))          ! number of levels for node 
          
          allocate ( depth_fast  (nkndi, nlvdi) )  ! depth                                 
          allocate ( vel_fast    (nkndi, nlvdi) )  ! current velocity                      
          allocate ( wtemp_fast  (nkndi, nlvdi) )  ! water temperature                     
          allocate ( wind_fast   (nkndi)        )  ! wind speed                            
          allocate ( atemp_fast  (nkndi)        )  ! air temperature                       
          allocate ( ice_cover_fast(nkndi)      )  ! ice cover (node area fraction)        
          allocate ( sal_fast    (nkndi, nlvdi) )  ! water salinity                        
          allocate ( light_fast  (nkndi)        )  ! incident light for the layer          
          allocate ( vol_fast    (nkndi, nlvdi) )  ! reactor volume                        
          allocate ( vol_old_fast(nkndi, nlvdi) )  !                                       
          allocate ( area_fast   (nkndi, nlvdi) )  ! reactor area                          
          allocate ( waves_fast  (nkndi, 2)     )  ! wind waves: 1-significant  wave height
                                                   !             2 - mean wave period         
          allocate(es_fast0 (nkndi,noslay,nsstate)) !sediment state vars to pass to kinetics                                        
          allocate(es_fast1 (nkndi,noslay,nsstate)) !sediment state vars to get from kinetics
          
!         For storage of state vars and derived vars for binary output  
          allocate(wc_output(nlvdi,nkndi,noutput))
          allocate(sed_output(noslay,nkndi,nsoutput))
          
!         For binary output      
          allocate(ilhkv_sed(nkndi))
          !allocate(hev_sed(neldi))
          
!         For output diagnostics vars for station nodes           
          allocate(dg    (nlvdi ,nkndi,nstate,ndiagvar))
          allocate(dg_sed(noslay,nkndi,nsstate,ndiagvar))

!        print *, nkn, nkndi, nel, neldi
!        print *, nlv, nlvdi,  nstate, nsstate
!        print *, noslay,  ndiagvar
!        print *, noutput, nsoutput,nconst, nsconst
!        print *, t_act
!        stop


!-----------------------------------------------------
!         Open file for WC variables strange values check
!-----------------------------------------------------
          if(bcheck) then
           ulogbio = ifileo(60,'logbio.txt','f','u')
           if( ulogbio .le. 0 ) then
              write(6,*) 'Cannot open/create file logbio.txt'
              stop
             else
              write(6,*) 'Aquabc_fem: File logbio.txt on unit: ',
     +         ulogbio
              write(ulogbio,*)
     +           '      STRANGE VALUES FOR AQUABC VARIABLES'
              write(ulogbio,*)
     +           'Var.No Value  Node Level Time  Checking moment'
           end if
          end if

!     --------------------------------------------------
!     initialize state variables with einit
!     --------------------------------------------------
          einit(1:nstate) = 0.
          e(:,:,:) = 0.
          
          do k=1,nkn    !loop on nodes
                lmax = ilhkv(k)
                do l=1,lmax
                 do i=1,nstate
                  e(l,k,i) = einit(i)
                 end do
               end do
          end do
!-------------------------------------------------------------
!         Initialise ascii variables
!-------------------------------------------------------------
          dg(:,:,:,:)       = 0.
          wc_output(:,:,:)  = 0.
          sed_output(:,:,:) = 0.

          e_fast0(1:nkndi,1:nstate,1:nlvdi) = 0.
          e_fast1(1:nkndi,1:nstate,1:nlvdi) = 0.

          es_fast0(1:nkndi,1:noslay,1:nsstate) = 0.
          es_fast1(1:nkndi,1:noslay,1:nsstate) = 0.
          es(:,:,:) = 0.

          lmax_fast   (1:nkndi)          = 0
          depth_fast  (1:nkndi,1:nlvdi)  = 0.
          vel_fast    (1:nkndi,1:nlvdi)  = 0.
          wtemp_fast  (1:nkndi,1:nlvdi)  = 0.
          wind_fast   (1:nkndi)          = 0.
          atemp_fast  (1:nkndi)          = 0.
          ice_cover_fast(1:nkndi)        = 0.
          sal_fast    (1:nkndi,1:nlvdi)  = 0.
          light_fast  (1:nkndi)          = 0.
          vol_fast    (1:nkndi,1:nlvdi)  = 0.
          vol_old_fast(1:nkndi,1:nlvdi)  = 0.
          area_fast   (1:nkndi,1:nlvdi)  = 0.
          waves_fast  (1:nkndi,1:2)      = 0.

!----------------------------------------------------------
!         initialize WQ state variables from external file
!         BS initialisation moved after call aquabcini
!----------------------------------------------------------

          call inicfil_aquabc('bioin',e,nstate)


!------------------------------------------------------------
!         set boundary conditions for all WC state variables
!------------------------------------------------------------

          nintp=2
          call get_first_dtime(dtime0)
          
          nvar = nstate
          call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +                     ,ebound,idbio)


!-----------------------------------------------------------------------
!         initialize parameters of eco model and output to ascii files
!----------------------------------------------------------------------

          call aquabcini(bsedim,par,par_sed,PHTIME,PHTAB,TEMPTIME,
     *                       TEMPTAB, NBIOTS, BIOTSNOD, BIOTSFUN,
     *                       NBIOTS_sed, BIOTSNOD_sed, BIOTSFUN_sed,
     *                       NDGTS, DGTSNOD, DGTSFUN,
     *                       NDGTS_sed,DGTSNOD_sed,DGTSFUN_sed)        
            

!------------------------------------------------------------------------
!         Initilises sediment properties and some parameters for BS
!         and  settling velocities and dissolved fractions for WC variables
!         see module aquabc_II_sed_ini
!------------------------------------------------------------------------
          call  sed_properties_first(bsedim)

!         Writing initial conditions to text file for defined nodes and binary file for all nodes

          idtc = nint(getpar('idtcon'))
          itmc = nint(getpar('itmcon'))

          !print *,'idtc=',idtc,'itmc=',itmc

          print *,'WRITING WC STATE VARIABLE INITIAL',
     *               ' VALUES FOR SELECTED NODES'

          do l = 1,nlvdi
           do k = 1,nkn
            do i=1,nstate
              wc_output(l,k,i) = e(l,k,i)  ! derived variables are not calculated yet here (zeros)
            end do                         ! It will be done in pelagic kinetic module together
           end do                          ! with calculation of derivatives for the first time step
          end do

          call biotser_write(1, 'wc',
     *                 wc_output, noutput, nstate, dg, NDIAGVAR,
     *                 ilhkv,nlvdi,
     *                 itmc,idtc,
     *                 NBIOTSMX,NBIOTS,BIOTSNOD,BIOTSFUN,
     *                 NDGTSMX, NDGTS, DGTSNOD, DGTSFUN)

          call init_binary_output_wc
          call write_binary_output_wc(dtime0)


!--------------------------------------
!         ADDITIONAL PREPARATION FOR BS
!--------------------------------------
          if (bsedim) then

!          arrays for binary output
            do i=1,nkndi
             ilhkv_sed(i)= noslay
            enddo

           !hev_sed(:) = sum(SED_DEPTHS(1:noslay))

           do i=1,noslay
            hlv_sed(i) = sum(SED_DEPTHS(1:i))            
           end do
           
           
           print *,'$$$$$$$$$$$$$$$$$$$$'           
           print *, 'BS levels ',hlv_sed
           print *,'$$$$$$$$$$$$$$$$$$$$'  

!          initialise ascii diagnostics
           dg_sed(:,:,:,:) = 0.

!          ---------------------------------
!          Initialisation BS state variables
!          ---------------------------------
           call inicfils_aquabc('biosin',es,nsstate)


           ! checking waves           
           if (isedi .eq. 1 .and. iwave .ne. 1) then
            print *,'AQUABC_II_FEM_INTERFACE:'
            print *, 'Waves calculation is not switched on'
            print *, 'Check parameter iwave in control file'
            stop
           end if

           ! initializing waves array
           waves_fast(1:nkndi,1:2) = 0.

!--------------------------------------------------------------------
!          Temporary initialize state BS variables. Now read from file.
!--------------------------------------------------------------------

!             do k=1,nkn
!              do i=1,nsstate
!
!              es(1,k,i) = esinit(i)
!              es(2,k,i) = esinit2(i)
!              es(3,k,i) = esinit3(i)
!
!              end do
!             enddo


           if (ifirst_sed .eq. 1) then
               call  sed_properties_ini(bsedim)
               call  sed_recalc_ini(es)
               ifirst_sed =0
           end if


           print *,'WRITING BS STATE VARIABLE INITIAL',
     *               ' VALUES FOR SELECTED NODES'
     
           !sed_output gets values of state and derived vars by sed_recalc_ini 

           call biotser_write(1, 'bs',sed_output,nsoutput,nsstate,
     *               dg_sed,NDIAGVAR_sed,
     *               ilhkv_sed,NOSLAY,
     *               itmc,idtc,
     *               NBIOTSMX_sed,NBIOTS_sed,BIOTSNOD_sed,BIOTSFUN_sed,
     *               NDGTSMX_sed, NDGTS_sed, DGTSNOD_sed, DGTSFUN_sed)
     
           call init_binary_output_bs      
           call write_binary_output_bs(dtime0)


          end if !bsedim
!----------------------------------
!         end of preparation for BS
!----------------------------------



        
     
!-------------------------------------------------------------------
!-------------------------------------------------------------------
! end of initialisation: icall==0
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------
      end if  !icall == 0 


!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------
! normal call (every time step)
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------


!--------------------------------------------------
!         set loadings in the internal areas. Not tested with new versions. fixme
!         1d array elaux is passed to AQUABC for a while. fixme
!--------------------------------------------------

!          call setload(eload, it, idt)
!          write(6,*)' loading set'

      what = 'lagvebio'

      kspec = 4350
      kspec = -1

      wsink = 0.

!-------------------------------------------------------------------
!    time management
!-------------------------------------------------------------------

      call get_timestep(dt)
      call get_act_dtime(dtime)
      it = dtime			!FIXME

      t0    = 0.
      dtday = dt / 86400
      tsec  = it
      tday  = it / 86400. + t0  !time in days, FEM 0 is day t0     


      if( bcheck ) call check_var('BEFORE aquabc',it,ulogbio,e,es)      


!-----------------------------------------------------------------------
!-------------------------------------------------------------------
!     data preparation loop on nodes for biological reactor
!-------------------------------------------------------------------
!-------------------------------------------------------------------

      ! getting sed. properties for each time step
      ! skip if is already called during initialisation
      if (ifirst_sed .ne. 0) then
       call  sed_properties_ini(bsedim)
      end if
      ifirst_sed = 1

!cccccccccccccccccccccccc
!  Start loop on nodes
!cccccccccccccccccccccccc          

      do k=1,nkn 

         lmax = ilhkv(k)  !maximum number of levels for the node
         lmax_fast(k) = lmax

! gets meteorological variables
         call meteo_get_heat_values(k,ITOT,tempair,rh,twetbulb,
     +                         windspeed,cloudc,pr)  !New routine produced by Georg

         atemp_fast  (k) = tempair
         wind_fast   (k) = windspeed

         
! now ITOT is instanteneous light not daily averaged
! as it is needed for new Ali-Smith routine

!        routines to introduce fraction of ice cover:

         call get_ice_cover(k,ice_cover)
         ice_cover_fast(k) =  ice_cover
         FDAY = 1. !used only for old equations compatibility
         !light_fast  (k)  = ITOT*(1-ice_cover)
         light_fast  (k)  = ITOT !No impact to the light by ice yet
         
         if (bsedim) then

          ! getting waves. Waves are not used for a while.
          !                Waves were intended to use by own resuspension module that is not ready yet
          !                
          call get_wave_values(k,wave_h,wave_mp,wave_pp,wave_d)
          waves_fast(k,1) = wave_h
          waves_fast(k,2) = wave_mp
 
         end if

!cccccccccccccccccccccccc
!          Loop on levels
!cccccccccccccccccccccccc

        do l=1,lmax

         call dvanode(l,k,-1,depth,volold,area)        !gets old depth, volume and area
         mode = +1
         call dvanode(l,k,mode,depthnew,vol,areanew)   !gets new depth, volume and area
         depth_fast  (k,l) = depth
         vol_fast    (k,l) = vol
         vol_old_fast(k,l) = volold
         area_fast   (k,l) = area

         if (is_nan(depth).or.depth.eq.0)
     +        then
          print *,
     +   'aquabc_II_fem_interface: Depth is NaN or zero:', depth,
     +         'on level: ', l, 'on node: ',k
          stop
         end if

        if (is_nan(vol).or.vol.eq.0)
     +        then
         print *,
     +   'aquabc_II_fem_interface: Volume is NaN or zero:', vol,
     +         'on level: ', l, 'on node: ',k
         stop
        end if

         if (is_nan(area).or.area.eq.0)
     +       then
          print *,
     +   'aquabc_II_fem_interface: Area is NaN or zero:', area,
     +         'on level: ', l, 'on node: ',k
          stop
         end if

!  Temperature and salinity
            call getts(l,k,t,s)     !gets temp and salt
            wtemp_fast  (k,l) = t
            sal_fast    (k,l) = s

!  Current velocity
            call getuv(l,k,u,v)     !gets velocities u/v
            vel = sqrt(u*u+v*v)
            vel_fast    (k,l) = vel

! Changing WC old values  to new in e_fast0 for the current step calculations            
            do i=1,nstate
            !print *,'----------------',k,i,l, nkn,nstate,lmax
             e_fast0(k,i,l) = e(l,k,i)
            end do

      end do

!ccccccccccccccccccccccccccccccccc
!            !End loop on levels
!ccccccccccccccccccccccccccccccccc
!------------------------------------------------------
      l = lmax
!----------------------------------------------

      end do ! on nodes
!ccccccccccccccccccccccccccccccccc
! End loop on nodes
!ccccccccccccccccccccccccccccccccc


!---------------------------------------------------------------------
!       --------------------------------------------------------------
!       end of loop on nodes for arrays preparation for call to AQUABC
!       --------------------------------------------------------------
!---------------------------------------------------------------------


            if(bsedim) then

            ! initial state for BS
            if(ifirst_it .eq. 1) then
            
             do ls=1,noslay
              do k=1,nkn
               do i=1,nsstate
                es_fast0(k,ls,i) = es(ls,k,i)
               end do
              end do
             end do
                         
             else
              es_fast0(1:nkn,1:noslay,1:nsstate)=
     +                    es_fast1(1:nkn,1:noslay,1:nsstate)
             end if  

            end if


            nlayers = maxval(lmax_fast(1:nkn))



!         print *,'pppppppppppppppppppppppppppppppppp'
!         print *, 'before AQUABC'   

            CALL AQUABC_II(nkn,lmax_fast,nlayers,
     *                  tday,dtday,
     *                  vol_fast,vol_old_fast,area_fast,
     *                  depth_fast,vel_fast,
     *                  wtemp_fast, sal_fast,
     *                  wind_fast,atemp_fast,
     *                  light_fast,FDAY,
     *                  ice_cover_fast,
     *                  elaux,
     *                  par,nconst,
     *                  e_fast0,e_fast1,nstate,
     *                  wc_output, noutput,
     *                  dg, NDIAGVAR,
     *                  bsedim,
     *                  par_sed,nsconst,
     *                  es_fast0, es_fast1,
     *                  nsstate,noslay,
     *                  waves_fast,
     *                  sed_output, nsoutput,
     *                  dg_sed, NDIAGVAR_sed,ifirst_it)

             !print *,'aquabc_fem_interface: Leaving AQUABC'
             
!        print *,'pppppppppppppppppppppppppppppppppp'      
!        print *, 'after AQUABC' 
!        print *, es_fast1(:,:,16)

! Assigning calculated state arrays by kinetics
            do k = 1,nkn
             lmax = lmax_fast(k)
             do i = 1,nstate
              do l = 1,lmax
               e(l,k,i) = e_fast1(k,i,l)
              end do
             end do
            end do


            if(bsedim) then
             do ls=1,noslay
              do k=1,nkn
               do i=1,nsstate
                es(ls,k,i) = es_fast1(k,ls,i)
               end do
              end do
             end do
             
!              print *,'pppppppppppppppppppppppppppppppppp'              
!              print *, 'after AQUABC renewal' 
!              print *, es(:,:,16)
            end if


!-------------------------------------------------------------------
!advection and diffusion
!-------------------------------------------------------------------

      if( bcheck ) call check_var('BEFORE advection',it,ulogbio,e,es)


      tdouble = t_act
      call bnds_read_new(what,idbio,tdouble)

!$OMP PARALLEL PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)

      do i=1,nstate

          call scal_adv(what,i
     +                          ,e(1,1,i),idbio
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)
!         calculates total mass for each variable
          call tsmass (e(1,1,i),1,nlvdi,tstot(i)) !mass control


      end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      if( bcheck ) call check_var('AFTER advection',it,ulogbio,e,es)

! assignment state variables after transport to the output array
! sediment output array comes assigned by state vars from BS routine


        do k = 1,nkn
          lmax = lmax_fast(k)
          do i=1,nstate
           do l = 1,lmax
             wc_output(l,k,i) = e(l,k,i)
           end do
         end do
        end do


!-------------------------------------------------------------------
!      WRITE FLUXES
!-------------------------------------------------------------------
 

      eflux(1:nlvdi,1:nkn,1) = wc_output(1:nlvdi,1:nkn,nstate+1)   ! Total N
      eflux(1:nlvdi,1:nkn,2) = wc_output(1:nlvdi,1:nkn,nstate+2)   ! Total P

      !print *, 'Writing AQUABC vars fluxes' 
      if(iflux .eq. 1) call fluxes_aquabc('.bfx', 300, nflux,eflux) 
      ! 300 is a base for output variables numbers
      ! Results written to *.bfx in g/s



!-------------------------------------------------------------------
!    WRITES OF RESULTS (FILE BINARY  AND TEXT ALSO DUMPS LAST STATE FOR REAPEATED RUNS)
!-------------------------------------------------------------------

!          write(17,*) 'AFTER scal3sh:','time=',it
!          write(17,1000) (ipv(k),(e(1,k,i),i=1,9),k=1,nkn)
!1000      format((I4,1x,9(F8.4,1x)))


! BINARY OUTPUT

        call write_binary_output_wc(dtime)
        
        if(bsedim) then
         call write_binary_output_bs(dtime)
        endif      


! WATER COLUMN TEXT OUTPUT

        call biotser_write(0, 'wc',wc_output, noutput, nstate,
     *                  dg, NDIAGVAR,
     *                  ilhkv,nlvdi,
     *                  itmc,idtc,
     *                  NBIOTSMX,NBIOTS,BIOTSNOD,BIOTSFUN,
     *                  NDGTSMX,NDGTS, DGTSNOD, DGTSFUN)

         if(idump.eq.1 .and. it .eq. itend) then
            call  dump_aquabc('dump_wc.dat',e,
     +                 nkndi,nkn,nlvdi,nstate)
         end if



! BOTTOM SEDIMENTS TEXT OUPUT
          if( bsedim ) then

           call biotser_write(0, 'bs',sed_output,nsoutput,nsstate,
     *                   dg_sed,NDIAGVAR_sed,
     *                    ilhkv_sed,NOSLAY,
     *                    itmc,idtc,
     *            NBIOTSMX_sed,NBIOTS_sed,BIOTSNOD_sed,BIOTSFUN_sed,
     *            NDGTSMX_sed, NDGTS_sed, DGTSNOD_sed, DGTSFUN_sed)

          if(idump.eq.1 .and. it .eq. itend) then
!          Dumping last time moment for repeated runs
!
!          Recalculates NH4 and PO4 as solute concentrations in porewater
            call sed_recalc_final(es)

            call  dump_aquabc('dump_bs.dat',es,
     +                 nkndi,nkn,noslay,nsstate)
           end if !idump

          end if !bsedim
!-----------------------------------------------

	    if(it .eq. itend) then
!           Calculates and prints  time required to run program
            call date_and_time(values=end_array)
            call calc_time (start_array,end_array)
           end if


!    -------------------------------------------------------------------
!    end of aquabc_fem_interface
!    -------------------------------------------------------------------

      end

c*************************************************************
!*************************************************************
!*************************************************************
!*************************************************************
!*************************************************************
!*************************************************************

!****************************************************************

        subroutine init_binary_output_wc

        use basin
        use levels
        use aquabc_II_vars

        implicit none

        integer id_file
        logical has_output_d


          call init_output_d('itmcon','idtcon',da_out)
          if( has_output_d(da_out) ) then
            call shyfem_init_scalar_file('bwc',noutput,.false.,id_file)
            da_out(4) = id_file
            write(6,*) 'Binary output for WC model initialized...'
          end if

        end
        
!****************************************************************

        subroutine init_binary_output_bs

        use basin
        use levels
        use aquabc_II_vars

        implicit none

        integer id_file
        logical has_output_d


          call init_output_d('itmcon','idtcon',da_outs)
          if( has_output_d(da_outs) ) then
          
           !print *,'hlv_sed = ',hlv_sed
           !stop
            call shyfem_init_scalar_file_hlv
     +          ('bbs',nsoutput,noslay,hlv_sed,id_file)
            da_outs(4) = id_file
            write(6,*) 'Binary output for BS model initialized...'
          end if

        end

!*************************************************************

        subroutine write_binary_output_wc(dtime)

        use basin
        use levels
        use aquabc_II_vars

        implicit none

        double precision dtime

        integer id_file,id_var,i
        logical next_output_d

        if( next_output_d(da_out) ) then

          id_file = nint(da_out(4))
          do i=1,noutput
            id_var = 100 + i
            call shy_write_scalar_record(id_file,dtime,id_var,nlvdi
     +                                          ,wc_output(1,1,i))
          end do
          print *,'Binary output for WC vars done'
        end if

        end
        
!*************************************************************

        subroutine write_binary_output_bs(dtime)

        use basin
        use levels
        use aquabc_II_vars

        implicit none

        double precision dtime

        integer id_file,id_var,i
        logical next_output_d

        if( next_output_d(da_outs) ) then

          id_file = nint(da_outs(4))
          do i=1,nsoutput
            id_var = 200 + i
            call shy_write_scalar_record(id_file,dtime,id_var,noslay
     +                                          ,sed_output(1,1,i))
          end do
          print *,'Binary output for BS vars done'
        end if

        end

c*************************************************************







	subroutine setload(eload, it, idt)

! sets up eload which is loading for specified areas
!
! the computed loadings in eload are in [g/(m**3 day)] == [mg/(l day)]
! the specified loadings in areaload are in [kg/day]
!
! variables to be specified:
!
! nimmis        total number of areas for which loading is specified
! nodes         total number of nodes used to identify all areas
! karee         node numbers that specify the areas of loading
! iaree         area numbers for the nodes [1-nimmis]
! areaload      total loadings [kg/day] for areas
!
! the node numbers in karee are external node numbers
!
!
! SUBROUTINE UPDATED BY ALI AND PETRAS TO ALLOW DYNAMIC LOADINGS,
! WHICH WILL BE READ FROM AN EXTERNAL FILE.
!
! CORPI, 20 July 2004 -----> Main updates on subroutine SETLOAD
!
!
! CORPI, 22 July 2004 -----> - SETLOAD corrected to overjump loading
!                              time intervals before simulation start
!
!                            - New header lines are added to the loading
!                              file to fill in some usefull information
!
!                            - A second (alternative) file format and
!                              structure has been developed. The new
!                              structure is a better alternative if
!                              time series with different time intervals
!                              are to be read for each load.
!
! CORPI, 23 July 2004 -----> - Error TAKING ONE DAYS LOADING FOR EACH
!                              TIME STEP has been fixed


	implicit none

      include 'param.h'
	integer nstate
	parameter(nstate=9)

!     MODIFIED BY ALI
!     CORPI, 15 July 2004
!     TAKE CARE
!     eload(3,neldim,nstate) -----> eload(nlvdim,nkndim,nstate)
	real eload(nlvdim,nkndim,nstate)

      integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
      common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real ev(13,1)
	common /ev/ev
	integer nen3v(3,1)
	common /nen3v/nen3v
	real hm3v(3,1)
	common /hm3v/hm3v

!	ADDED BY PETRAS AND ALI
!     CORPI, 19 July 2004
!     New integer variable added nimmax
!
!     nimmax : Maximum number of loadings allowed for the compiled
!              executable image. For this exectuteable 50. If more
!              loadings are needed, increase nimmax and recompile.
!
	integer nimmax
	parameter (nimmax=50)

!     Actual mumber of loadings. 1..nimmax
	integer nimmis
	save nimmis

	real volaux(nimmax)
	real areaload(nimmax,nstate)
	save volaux,areaload


!     ADDED BY ALI
!     CORPI, 22 July 2004
!     New variable pareaload, narealoaf
!
!     nareaload(j, jj) : Next loading for load j, state variable jj

	real nareaload(nimmax,nstate)
      save nareaload

!     ADDED BY ALI
!     CORPI, 19 July 2004
!
!     New variables it, idt, itload, preitl
!
!	it	   : time in seconds
!	idt	   : time step in seconds
!	itload : time of load interval
!	preitl : time of prevois load interval

	integer it
	integer idt
	integer itload
	integer preitl

	save itload, preitl


!     ADDED BY ALI
!     CORPI, 22 July 2004
!
!     New variables it, idt, itload, preitl
!
!	itloa2(j) : time of load interval for 2nd type loading file load j
!	preit2(j) : time of prevois load interval for second type load j

	integer itloa2(nimmax) !time of load interval for 2nd type loading file
	integer preit2(nimmax) !time of prevois load interval for second type

	save itloa2, preit2

    	integer aree(nkndim)


!	ADDED BY PETRAS AND ALI
!     CORPI, 19 July 2004
!     New integer variable added nimmax
!
!     nimmax : Maximum number of nodes with loadings allowed for
!              the compiled executable image. For this exectuteable
!              5000. If moreloadings are needed, increase nimmax
!              and recompile.
!
      integer nodmax
	parameter(nodmax=5000)

!     Actual mumber of nodes with loadings. 1..nodmax

	integer nodes
	save nodes !ADDED BY PETRAS 12-10-2004

!     ADDED BY ALI
!     CORPI, 22 July 2004
!     New integer variable added lftype
!
!     lftype : Loading file type
!
!              lftype = 1 ----> Loading file keeps all the loading
!                               information
!
!              lftype = 2 ----> Loading file keeps basic loading
!                               information and names of the time series
!                               files for each load
      integer lftype

	integer karee(nodmax)
	integer iaree(nodmax)
	save karee,iaree


!	ADDED BY PETRAS AND ALI
!     CORPI, 19 July 2004
!     New integer variable added
!
!     icall : Is it necessary to read information about
!                  loading areas from main loads  file
!             icall = 0 ---> Need to read
!             icall = 1 ---> information is already readed

      integer icall
	save icall
      data icall /0/


!	ADDED BY ALI
!     CORPI, 19 July 2004
!     New integer variable added
!
!     ifirst : Is it the the first time for reading loading data
!              ifirst = 0 ---> Reading loading data for the first time
!              ifirst = 1 ---> Reading loading data not for the first time

      integer ifirst
      save ifirst
      data ifirst /0/


!	ADDED BY ALI
!     CORPI, 22 July 2004
!     New integer variable added
!
!
!     ifirs2 : Is it necessary to read the next time intervall for loads
!
!              ifirs2(j) = 0 ---> Need to read the next time intervall
!                                 for load j
!
!              ifirs2(j) = 1 ---> Do not need to next time intervall
!                                 for load j

      integer ifirs2(nimmax)
	save ifirs2


!	ADDED BY ALI
!     CORPI, 19 July 2004
!     New integer variable added
!
!     inext : Is it necessary to read the next time intervall for loads
!             inext = 0 ---> Need to read the next time intervall
!             inext = 1 ---> Do not need to next time intervall

      integer inext
	save inext
      data inext /0/

!	ADDED BY ALI
!     CORPI, 22 July 2004
!     New integer variable added
!
!
!     inext2 : Is it necessary to read the next time intervall for loads
!
!              inext2(j) = 0 ---> Need to read the next time intervall
!                                 for load j
!
!              inext2(j) = 1 ---> Do not need to next time intervall
!                                 for load j

      integer inext2(nimmax)
	save inext2


!	ADDED BY ALI
!     CORPI, 19 July 2004
!     New integer variable added
!
!     ilast : Is last loading time interval read
!             ilast = 0 ---> Last loading time interval not read
!             inext = 1 ---> Last loading time interval read

      integer ilast
	save ilast
      data ilast /0/

!	ADDED BY ALI
!     CORPI, 22 July 2004
!     New integer variable added
!
!
!     ilast2 : Is it necessary to read the next time intervall for loads
!
!              ilast2(j) = 0 ---> Last loading time interval not read for load j
!              ilast2(j) = 1 ---> Last loading time interval read for load j

      integer ilast2(nimmax)
	save ilast2

	logical berror
	integer k,ie,ii,ia,i
	integer itype
	real area
	real litri,kgs
	real fact,rlfact
	real hdep,vol,load

	real getpar


!	ADDED BY ALI
!     CORPI, 02 August 2004
      integer ipv(nkndim)	!external node numbers
      common  /ipv/ipv


!	ADDED BY PETRAS AND ALI
!     CORPI, 19 July 2004
!     New variables added : nb, file, irec
!
!     nb     : File number for the loadings file
!     header : Header information
!     irec   : Record number
!     ivar   : Variable number
!     j, jj  : General purposed counter for array indices, ...
!
      integer nb
	save nb

	character*90 header

      integer irec, ivar
	integer j, jj
!
!     ADDED BY ALI
!     CORPI, 22 July 2004
!     New variables added : ltsferr, ltsfun, ltsfnm
!
!     ltsfer : Used for error checkong when opening loading time series file
!     ltsfun : Loading time series file units
!     ltsfnm : Loading time series file names

      integer ltsfer
      data ltsfer /0/

      integer ltsfun(nimmax)
	save ltsfun

	character*256 ltsfnm(nimmax)
      save ltsfnm

! loading is kg/day
!
!
!	loading for areas [kg/day]
!
	integer ifileo
        integer max_int

        character*80 file

!     FIX FILE NAME FOR THIS VERSION
!      file = 'INPUT/loads.dat'


!	ADDED BY PETRAS AND ALI
!     CORPI, 19 July 2004
!
	if( icall .eq. 0 ) then

!         OPEN THE MAIN LOADINGS FILE

            call getfnm('bioload',file)
	    nb = ifileo(90,file,'f','old')

	    if( nb .le. 0 ) goto 97


!         Initialize the arrays
          do j = 1, nodmax
              iaree(j) = 0
	        karee(j) = 0
	    end do

	    do j = 1, nimmax

	        do jj = 1, nstate
	            areaload (j,jj) = 0.0
                  nareaload(j,jj) = 0.0
	        end do

              ltsfun(j) = 0
	        ltsfnm(j) = ''
              ifirs2(j) = 0
			inext2(j) = 0
              ilast2(j) = 0
              itloa2(j) = 0
              preit2(j) = 0

	    end do

		  write(6,*) 'READING LOADING FILE...'

	    preitl = 0
	    ivar = 0

!         Read the RECORD 1
          irec = 1


!         Read two header lines
          read(nb, 5040, err=101) header
          read(nb, 5040, err=101) header

		  write(6,*) 'RECORD 1 OF THE LOADING FILE READ'

!         Read the RECORD 2
          irec = 2

!         Read two header lines
          read(nb, 5040, err=101) header
          read(nb, 5040, err=101) header
		  write(6,*) 'RECORD 2 OF THE LOADING FILE READ'

!         Read the RECORD 3
          irec = 3

	    read(nb, 5010, err=90) nimmis, nodes, lftype

	    if(nimmis.gt.nimmax) goto 95
	    if(nodes.gt.nodmax) goto 96

	    write(6,*) ''
		  write(6,*) 'RECORD 3 OF THE LOADING FILE READ'
          write(6,*) '---------------------------------'
	    write(6,*) 'TOTAL NUMBER OF LOADINGS        : ', nimmis
	    write(6,*) 'NUMBER OF NODES RECIEVING LOADS : ', nodes
	    write(6,*) 'TYPE OF THE LOADING FILES       : ', lftype
	    write(6,*) ''

!         Read the RECORD 4
          irec = 4

!         Read two header lines
          read(nb, 5040, err=101) header
          read(nb, 5040, err=101) header
		  write(6,*) 'RECORD 4 OF THE LOADING FILE READ'


!         Read the RECORD 5
          irec = 5

	    do j=1, nodes

	        read(nb, 5010, err=93) iaree(j), karee(j)

	        if(iaree(j).lt.1.or.iaree(j).gt.nimmis) then
			    goto 91
	        end if

			if(karee(j).lt.1.or.karee(j).gt.max_int(ipv, nkn)) then
			    goto 92
	        end if

          end do

		  write(6,*) 'RECORD 5 OF THE LOADING FILE READ'

!         Read the RECORD 6
          irec = 6

!         Read two header lines
          read(nb, 5040, err=101) header
          read(nb, 5040, err=101) header
		  write(6,*) 'RECORD 6 OF THE LOADING FILE READ'

	    if(lftype.eq.2) then

!             Read the RECORD 7 of LOADING FILE TYPE 2
              irec = 7

		      write(6,*) ''
		      write(6,*) 'READING RECORD 7 OF THE LOADING FILE'
              write(6,*) '===================================='
              write(6,*) ''

		      do j = 1, nimmis
!                 Read loading time series file name
	            read(nb, 5050, err=201) ltsfnm(j)

!				  Open the loading time series file
                  ltsfun(j) = ifileo((80 + j),ltsfnm(j), 'f','old')

!                 Write loading time series file information on terminal
                  write(6,*) 'FOR LOAD ', j
	            write(6,*) '------------------'

	            if(ltsfun(j).le.0) then
                      write(6,*) 'PROBLEMS ENCOUNTERED WHEN OPENING ',
     +                           'TIME SERIES FILE FOR LOADING ', j
	                write(6,*) ''
			        ltsfer = 1
	            else
                      write(6,*) 'Unit of the loading time series file:'
     +                           ,ltsfun(j)
                      write(6,*) 'Name of the loading time series file:'
     +                           ,ltsfnm(j)
	                write(6,*) ''
	            end if

	        end do

			  close(nb)

	        if(ltsfer.eq.1) goto 200

		end if

!         extern to intern
		call n2int(nodes,karee,berror)

	    if( berror) stop 'error stop: loading'

		icall = 1

	end if

	if((ifirst.eq.0).and.(lftype.eq.1)) then

    1     continue
!         Read the RECORD 7 of LOADING FILE TYPE 1
          irec = 7

          read(nb, 5030, err=98) itload
		  write(6,*) ''
		  write(6,*) 'RECORD 7 OF THE LOADING FILE READ'

		  preitl = itload


!         ADDED BY ALI
!         CORPI, 22 July 2004
!
!         The following if structure overjumps loading time interval
!         before simulation start
		  if(itload.lt.(it-idt)) then

              if(itload.lt.0) then
                  write(6,*) 'FOR THIS SIMULATION NO LOADS WILL BE READ'
                  write(6,*) 'ZERO LOADING ASSUMED'
	            ilast = 1
                  close(nb)
	        else
                  write(6,*) 'SIMULATION START : ', (it - idt),
     +			' START OF THE LOAD INTERVAL : ', itload

                  write(6,*)'NEXT LOADING TIME INTERVAL WILL BE READ...'

!                 Read the RECORD 8 of LOADING FILE TYPE 1
                  irec = 8

		          do j = 1, nimmis
	                read(nb, 5020, err=94)
     +				(areaload(j,ivar), ivar=1,nstate)
		          end do

		          write(6,*) 'RECORD 8 OF THE LOADING FILE READ'

                  goto 1

	        end if

		end if

		ifirst = 1

      end if

	if(((inext.eq.0).and.(ilast.eq.0)).and.lftype.eq.1) then

!         Read the RECORD 8 of LOADING FILE TYPE 1
          irec = 8

		  do j = 1, nimmis
	        read(nb, 5020, err=94)(areaload(j,ivar), ivar=1,nstate)
		  end do

		  write(6,*) 'RECORD 8 OF THE LOADING FILE READ'

		  preitl = itload
		  inext = 1

!         Read the RECORD 7 of LOADING FILE TYPE 1 - NEXT TIME INTERVAL
          irec = 7

          read(nb, 5030, err=98) itload
		  write(6,*) 'RECORD 7 OF THE LOADING FILE READ'

	    if(itload.le.preitl) then

              if(itload.lt.0) then
	            write(6,*) 'NO MORE LOADING TIME INTERVALS'
				  ilast = 1
	            close(nb)
	        else
		          goto 100
              end if

		end if

	end if

!     CHECK IF TIME FOR THE NEXT LOADING INTERVAL
      if((((it+idt).ge.itload).and.(ilast.eq.0)).and.lftype.eq.1) then
	    write(6,*) 'NEW LOADING INTERVAL STARTING NEXT TIME STEP'
		  inext = 0
	end if


!     READ LOADING DATA FOR LOADIND FILE TYPE 2
      if(lftype.eq.2) then

	    do j = 1, nimmis

	        if(ifirs2(j).eq.0) then

    2             continue

!                 Read the RECORD 8 of LOADING FILE TYPE 2
                  irec = 8

	            read(ltsfun(j), 5060, err=94)
     +			  itloa2(j), (areaload(j,ivar), ivar=1,nstate)

		          write(6,*) 'RECORD 8 OF THE LOADING FILE READ'

				  if(itloa2(j).lt.(it-idt)) then

                      if(itloa2(j).lt.0) then
                          write(6,*) 'FOR THIS SIMULATION NO LOADS WILL'
     +					         , ' BE READ FOR LOADING ', j
                          write(6,*) 'ZERO LOADS ASSUMED FOR ', j

						  do jj=1, nstate
                              areaload(j,jj) = 0.0
                          end do

						ilast2(j) = 1
                          close(ltsfun(j))
	                else
                          write(6,*) 'SIMULATION START : ', (it - idt),
     +			        ' START OF THE LOAD INTERVAL FOR LOADING ', j,
     +                    ' : ', itloa2(j)

                          write(6,*)'READING THE NEXT LOADING TIME ',
     +					          'INTERVAL'

                          goto 2

	                end if

		        end if

                  ifirs2(j) = 1

              end if

			if((inext2(j).eq.0).and.(ilast2(j).eq.0)) then

		        preit2(j) = itloa2(j)

!                 Read the RECORD 8 of LOADING FILE TYPE 2
                  irec = 8

	            read(ltsfun(j), 5060, err=94)
     +			itloa2(j), (nareaload(j,ivar), ivar=1,nstate)

	            inext2(j) = 1

                  if(itloa2(j).le.preit2(j)) then

                      if(itloa2(j).lt.0) then
	                    write(6,*) 'NO MORE LOADING TIME INTERVALS'
				        ilast2(j) = 1

						close(ltsfun(j))
	                else
		                goto 100
                      end if

		        end if

	        end if

!             CHECK IF TIME FOR THE NEXT LOADING INTERVAL
              if(((it+idt).ge.itloa2(j)).and.(ilast2(j).eq.0)) then
	            write(6,*) 'NEW LOADING INTERVAL STARTING NEXT TIME ',
     +			           'STEP FOR LOAD ', j

				  do jj=1, nstate
                       areaload(j,jj) = nareaload(j,jj)
                  end do

		        inext2(j) = 0
              end if

          end do

	end if

 5010 FORMAT(3I10)
 5020 FORMAT(9F10.0)
 5030 FORMAT(I10)
 5040 FORMAT(A90)
 5050 FORMAT(A256)
 5060 FORMAT(I10, 9F10.0)

	litri = 1000.	!litri in m*3
	kgs  = 10.e+6	!mg in kg

!	rlfact = getpar('flbio')
	rlfact = 1.


!     MODIFIED BY ALI
!     CORPI, 23 July 2004
!     ERROR ABOUT TAKING ONE DAYS LOADING FOR ECAH TIME STEP FIXED
!
!     fact = (rlfact*kgs/litri)
!     fact = -------> (rlfact*kgs/litri) * (idt/86400.)

!     [kg/m**3] -----> [mg/l] ------> kg/day ---------> (kg/s)*idt
	fact = (rlfact*kgs/litri) * (idt/86400.)



!     intialize

	do i=1,nimmis
	  volaux(i) = 0.
	end do

	do k=1,nkn
	  aree(k) = 0
	end do

	do i=1,nodes
	  k = karee(i)
	  itype = iaree(i)
	  aree(k) = itype
	end do

!     compute total volume for all areas given -> store in volaux

	do ie=1,nel
	  area = 4. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    hdep = hm3v(ii,ie)
	    ia = aree(k)
	    if( ia .gt. nimmis ) stop 'error stop ia'
	    if( ia .gt. 0 ) then
		   volaux(ia) = volaux(ia) + area * hdep
	    end if
	  end do
	end do

!     compute and set loading in eload [g/(m**3 day)] == [mg/(l day]

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ia = aree(k)
	    vol = volaux(ia)
	    do i=1,nstate
!	 Corrected by Petras to avoid division by zero:
		  if( ia .le. 0 ) then
		    load = 0.
	    else
		    load = fact * areaload(ia,i) / vol
	    endif

!           MODIFIED BY ALI
!           CORPI, 15 July 2004
!           TAKE CARE
!           eload(ii,ie,i) -----> eload(1,k,i)
		  eload(1,k,i) = load
	    end do
	  end do
	end do


      return

   90	continue
	write(6,*) 'read error in record = ',irec
	write(6,*) '... reading file ', 'loads.dat'
      write(6,*) 'Please check your loading file'
	stop 'error stop setload'

   91 continue
      write(6,*) 'error in record = ', irec, ' at row = ', j, 'node = ',
     +            karee(j)
	write(6,*) 'undefined loading area'
      write(6,*) 'Please check your loading file'
	stop 'error stop setload'

   92 continue
      write(6,*) 'error in record = ', irec,   ' at row = ', j,
     +           ' loading area = ', iaree(j), ' node = ', karee(j)
	write(6,*) 'undefined nodes for loadings'
	write(6,*) 'Please check your loading file'
	stop 'error stop setload'

   93	continue
	write(6,*) 'read error in record = ',irec,' at row = ', j
	write(6,*) '... reading file ','loads.dat'
      write(6,*) 'Please check your loading file'
	stop 'error stop setload'

   94	continue
	write(6,*) 'read error in record = ',irec,' at row = ', j,
     +           ' ivar = ', ivar
	write(6,*) '... reading file ','loads.dat'
      write(6,*) 'Please check your loading file'
	stop 'error stop setload'

   95 continue
      write(6,*) 'Array dimension error :'
	write(6,*) 'This executeable image was compiled for ', nimmax,
     + ' loadings but you use ', nimmis , ' loadings. Please deacrease',
     + ' the number of loadings in the loads input file RECORD 3 or ',
     + ' change nimmax parameter in SUBROUTINE SETLOAD to ', nimmis,
     + ' or greater and recompile.'
	stop 'error stop setload'

   96 continue
      write(6,*) 'Array dimension error :'
	write(6,*) 'This executeable image was compiled for ', nodmax,
     + ' nodes for loadings but you use ', nodes , ' nodes.',
     + ' Please deacrease the number of nodes in the loads input ',
     + 'file RECORD 3 or change nodmax parameter in SUBROUTINE SETLOAD',
     + ' to ', nodes, ' or greater and recompile.'
	stop 'error stop setload'

   97	continue
	write(6,*) 'Cannot open loadings file ','loads.dat'
	stop 'error stop setload'

   98 continue

	if(preitl.eq.0) then
	    write(6,*) 'read error in record = ',irec,' Please check ',
     +               'RECORD 7 of the first loading interval.'
	else
	    write(6,*) 'read error in record = ',irec,' Please check ',
     +               'RECORD 7 of the loading interval next to the ',
     +               'interval staring at ', preitl, ' secs.'
      end if
	stop 'error stop setload'

  100 continue
      write(6,*) 'Time error :'
      write(6,*) 'Next loading time interval starts before the ',
     +           'current time interval. Please check RECORD 7 ',
     +           'after the time interval starting at ', preitl,
     +           ' secs.'
	stop 'error stop setload'

  101	continue
	write(6,*) 'read error in record = ',irec
	write(6,*) 'Please check the header information format'
	stop 'error stop setload'

  200	continue
	write(6,*) 'Error when reading loading time series file(s)'
      write(6,*) 'Please check your loading file'
	stop 'error stop setload'

  201	continue
	write(6,*) 'read error in record = ',irec
	write(6,*) 'Plese check the time series file name format'
	stop 'error stop setload'

	end !setload

!*************************************************************
!*************************************************************

	subroutine check_var(title,it,ulog,e,es)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! checks WC and BS vars for strange values
!
! Uses:
!   subroutine check2Dbio - cheks for strange values for given ranges
!                           and variable numbers
! Developped by G.Umgiesser
! Modified for use for biogeochemical variables by P.Zemlys
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use levels, only: nlvdi
      use basin, only: nkndi
      use aquabc_II_vars
      
	implicit none

      !include 'param.h'
      !include 'aquabc_II.h'

      integer ulog,it 
       
      integer nslayer   


      character*(*) title
      real e(nlvdi,nkndi,nstate)	        !state vector
      real es(NOSLAY,nkndi,nsstate)		!sediment state variables

        !integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        !common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        !integer nlvdi,nlv
        !common /level/ nlvdi,nlv

        character*16 text
        integer i

        nslayer = NOSLAY

	!write(6,*) 'check_var: ',title

        text = '!!! BIO CHECK'
! Changed by Petras 11 January 2005 to pass number of variable
       do i=1,nstate
        call check2Dbio(it,ulog,nlvdi,nkndi,70+i,e(1,1,i),
     *                     -0.001,25.,text,title)
       end do

		do i=1,nsstate
        call check2Dbio(it,ulog,nslayer,nkndi,100+i,es(1,1,i),
     *                     -0.001,100000.,text,title)
        end do

      end

!***************************************************************

      subroutine check2Dbio(it,ul,nlv,n,nvar,a,vmin,vmax,textgen,text)

! tests array for nan and strange values for EUTRO variables
! Made from check2Dr by P.Zemlys 11 January 2005

        !use levels
        use basin, only: ipv
        
        implicit none

        !include 'param.h'

        !integer ipv(nkndi) !external numbers of nodes
        !common /ipv/ipv
        integer n   ,nlv
        integer nvar         !state variable number
        real a(nlv,1)     !array of variable values for levels and nodes
        real vmin,vmax       !minimal and maximal allowed values
        character*(*) textgen,text ! '***BIO CHECK' and text indicating the
                                   ! time of checking

        logical debug,bval
        integer inan,iout,i,l,ul,it
        real val

        logical is_r_nan

        bval = vmin .lt. vmax
        debug = .true.
        inan = 0
        iout = 0

        return

        do i=1,n
          do l=1,nlv
            val = a(l,i)
            if( is_r_nan(val) ) then
              inan = inan + 1
              if( debug ) write(ul,'(I3,G11.3,I5,I2,I10,1X,A16)') nvar,
     +                        val,ipv(i),l,it,text
            else if(bval .and. (val .lt. vmin .or. val .gt. vmax)) then
              iout = iout + 1
              if( debug ) write(ul,'(I3,G11.3,I5,I2,I10,1x,A16)') nvar,
     +                        val,ipv(i),l,it,text
            end if
          end do
        end do

        if( inan .gt. 0 .or. iout .gt. 0 ) then
          write(6,'(2A13,A4,A16,A2,2(A5,I4))') 'CHECK2DBIO: ',
     +  textgen," (",text,") ",
     +  ' NaN=',inan,' OUT=',iout
        end if

        end  !check2Dbio

c*************************************************************


      function max_int(array, index)

      integer array(index)
	integer maximum

	integer index
	integer i

      maximum = array(1)

	do i=2, index

          if(array(i).ge.maximum) then
	        maximum = array(i)
	    end if

	end do

      max_int = maximum

	end  !max_int

c*********************************************************

       subroutine inicfil_aquabc(name,var,nvar)
! initializes state variables nodal values for WC model from file

        use levels
        use basin 

        implicit none

        !include 'param.h'
        !include 'nbasin.h'

        character*(*) name           !name of variable
        real var(nlvdi,nkndi,1)    !variable to set for water column (name='bioin')

        real varval
        integer nvar
        integer ftype               !1-  homogeneous initial cond. for WC
                                    !2 - heterogeneous initial conditions.
                                    !3 - homogeneous layers
                                    !4 - heterogenous for repeated runs

        !integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        !common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        !real hlv(1)       !hlv(i)absolute depth of layer i bottom
        !common /hlv/hlv

        character*80 file
        integer nb,irec
        integer nkk,lmax
        integer l,k,i,j
        integer ivars,ivar
        real val
        real rlaux(nlvdi) !absolute depth of bottom of layer i
                           !nlvdim - total number of vertical levels

        integer ifileo
        real varval_array(2000)
!-------------------------------------------------------
! get file name
!-------------------------------------------------------

        call getfnm(name,file)

        if( file .eq. ' ' ) return      !nothing to initialize

!-------------------------------------------------------
! open file
!-------------------------------------------------------

!     MODIFIED BY ALI AND PETRAS
!     CORPI, 16 July 2004C
!     nb = ifileo(55,file,'unform','old')
!      -------> nb = ifileo(55,file,'f','old')
!
        nb = ifileo(55,file,'f','old')
        if( nb .le. 0 ) goto 97

!-------------------------------------------------------
! read first record
!-------------------------------------------------------

      ivar = 0

        irec = 1

!     MODIFIED BY ALI AND PETRAS
!     CORPI, 16 July 2004

!
!     read(nb,err=90) nkk,lmax,ivars
!      -----> read(nb, 5000, err=90) nkk,lmax,ivars
!
! Skip 4 lines. Added by Petras 12-12-2004
        read(nb,*)
        read(nb,*)
        read(nb,*)
        read(nb,*)
! Read control information
        read(nb, *, err=90) nkk,lmax,ivars,ftype
        
        write(*,*) '********************************************'
        write(*,*) '*    Initial condions info for water column*'
        write(*,*) '*------------------------------------------*'
        write(*,*) '    nkk : ', nkk
        write(*,*) '   lmax : ', lmax
        write(*,*) '  ivars : ', ivars
        write(*,*) '  ftype : ', ftype
        write(*,*) '********************************************'        
        
       
 5010 FORMAT(4I5)

      if( nkk .ne. nkn)  goto 99

      if(lmax .gt. nlvdi ) goto 99

      if( ivars .ne. nvar ) goto 96
      if(ftype .ne. 1 .and. ftype .ne. 2
     +  .and. ftype .ne. 3 .and. ftype .ne. 4) goto 91

!********************************************************
! FILE TYPE 2 (spatialy heterogeneous initial conditions)
!********************************************************
      if(ftype .eq. 2) then
!-------------------------------------------------------
! read second record (only if lmax > 0)
!-------------------------------------------------------             !

        if( lmax .gt. 1 ) then          !changed from 0 to 1 (5.3.2004) LMAX
          irec = 2
!       MODIFIED BY ALI AND PETRAS
!       CORPI, 16 July 2004

!
!
!       read(nb, err=90) (rlaux(l),l=1,lmax)
!       -----> read(nb, 5020, err=90) (rlaux(l),l=1,lmax)
!
          read(nb, 5020, err=90) (rlaux(l),l=1,lmax)
 5020   FORMAT(8F10.0)

          do l=1,lmax
            if( hlv(l) .ne. rlaux(l) ) goto 98
          end do
        else
          lmax = 1
        end if

!-------------------------------------------------------
! read data records
!-------------------------------------------------------

        irec = 3

        do ivar=1,nvar

!     MODIFIED BY ALI AND PETRAS
!     CORPI, 16 July 2004C
!     read(nb, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
! --> read(nb, 5030, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
!

        read(nb, 5030, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
 5030 FORMAT(F5.0)

      end do

!-------------------------------------------------------
! initialize the other levels if only surface is given
!-------------------------------------------------------

        do ivar=1,nvar
          do k=1,nkn
            val = var(lmax,k,ivar)
            do l=lmax+1,nlvdi
              var(l,k,ivar) = val
            end do
          end do
        end do

      end if ! end of ftype=2

!******************************************************
! FILE TYPE 1 (spatialy homogeneous initial conditions)
!******************************************************
      if(ftype .eq. 1) then
      irec=3
      do ivar=1,nvar
       read(nb, *, err=90) varval
       print *,'Initial condition: ','var= ',ivar,'value=', varval
       do l=1,lmax
        do k=1,nkn
         var(l,k,ivar)=varval
        end do
       end do
      end do


      end if !end of ftype 1

!******************************************************
! FILE TYPE 3 spatialy homogeneous initial conditions for each level
!******************************************************
      if(ftype .eq. 3) then
      irec=3
      do ivar=1,nvar
       read(nb, *, err=90) (varval_array(l),l=1,lmax)

       print *,'Initial condition: ','var= ',ivar,'value=',
     *       (varval_array(l),l=1,lmax)

       do l=1,lmax
        do k=1,nkn
         var(l,k,ivar) = varval_array(l)
        end do
       end do

      end do

      end if !end of type 3
!******************************************************
! FILE TYPE 4 dumped state for repeated runs
!******************************************************
      if(ftype .eq. 4) then

       do i = 1, nvar
          do j = 1, nkn
              read(nb,*,err=90) (var(k,j,i),k=1,lmax)
          end do
       end do

      end if !end of type 4

!-------------------------------------------------------
! reading done -> close file
!-------------------------------------------------------
 1001 continue
        close(nb)



!-------------------------------------------------------
! end of initialisation
!-------------------------------------------------------

        write(6,*) 'Succesfull initialization for ',name,' from file '
        write(6,*) file

        return

   90   continue
        write(6,*) 'read error in record = ',irec,' ivar = ',ivar
        write(6,*) '... reading file',file
        stop 'error stop inicfil_aquabc'
   91   continue
        write(6,*)
     +   'bad file type descriptor: value 1, 2 or 3 is allowed!'
        write(6,*) '... reading file',file
        stop 'error stop inicfil_aquabc'
   96   continue
        write(6,*) 'ivars not compatible with nvar: ',ivars,nvar
        stop 'error stop inicfil_aquabc'
   97   continue
        write(6,*) 'Cannot open file ',file
        stop 'error stop inicfil'
   98   continue
        write(6,*) 'levels are not the same from init file ',file
        write(6,*) (hlv(l),l=1,lmax)
        write(6,*) (rlaux(l),l=1,lmax)
        stop 'error stop inicfil_aquabc'
   99   continue
        write(6,*) 'INICFIL_AQUABC:' 
        print *, 'Parameters are not the same from init file ',file
        write(6,*) 'nkn, lmax from file  : ',nkk,lmax
        write(6,*) 'nkn, nlvdi from model : ',nkn,nlvdi
        stop 'error stop inicfil_aquabc'

       end !inicfil_aquabc

!*********************************************************
!*********************************************************

       subroutine inicfils_aquabc(name,var,nvar)
! initializes variables nodal value  for bottom sediment from file

        use aquabc_II_vars
        use aquabc_II_sed_ini, only: sed_depths 
               
        use levels
        use basin

        implicit none

        !include 'param.h'
        !include 'aquabc_II.h'
        !include 'nbasin.h'

        character*(*) name           !name of variable

        real var(noslay,nkndi,1)   !variable to set for bottom sediments (name='biosin')
        real varval
        integer nvar
        integer ftype               !1-  homogeneous initial cond. for WC
                                    !2 - heterogeneous initial conditions.
                                    !3 - homogeneous layers for BS
                                    !4 - heterogenous for repeated runs

        !integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        !common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        !real hlv(1)       !hlv(i)absolute depth of bottom of layer i
        !common /hlv/hlv

        character*80 file
        integer nb,irec
        integer nkk,lmax
        integer l,k,i,j
        integer ivars,ivar
        real val
        real rlaux(noslay) !absolute depth of bottom of layer i
                           !noslay - total number of vertical levels

        integer ifileo
        real varval_array(2000)
!-------------------------------------------------------
! get file name
!-------------------------------------------------------

        call getfnm(name,file)

        if( file .eq. ' ' ) return      !nothing to initialize

!-------------------------------------------------------
! open file
!-------------------------------------------------------

!     MODIFIED BY ALI AND PETRAS
!     CORPI, 16 July 2004
!     nb = ifileo(55,file,'unform','old')
!      -------> nb = ifileo(55,file,'f','old')
!
        nb = ifileo(55,file,'f','old')
        if( nb .le. 0 ) goto 97

!-------------------------------------------------------
! read first record
!-------------------------------------------------------

      ivar = 0

        irec = 1

!     MODIFIED BY ALI AND PETRAS
!     CORPI, 16 July 2004

!
!     read(nb,err=90) nkk,lmax,ivars
!      -----> read(nb, 5000, err=90) nkk,lmax,ivars
!
! Skip 4 lines. Added by Petras 12-12-2004
        read(nb,*)
        read(nb,*)
        read(nb,*)
        read(nb,*)
! Read control information
        read(nb, *, err=90) nkk,lmax,ivars,ftype
 5010 FORMAT(4I5)

         
        write(*,*) '********************************************'
        write(*,*) '*    Initial condions info for sedimrnts    *'
        write(*,*) '*------------------------------------------*'
        write(*,*) '    nkk : ', nkk
        write(*,*) '   lmax : ', lmax
        write(*,*) '  ivars : ', ivars
        write(*,*) '  ftype : ', ftype
        write(*,*) '********************************************'
 
      if( nkk .ne. nkn)  goto 99


       if(lmax .ne. noslay ) goto 99


      if( ivars .ne. nvar ) goto 96
      if(ftype .ne. 1 .and. ftype .ne. 2
     +  .and. ftype .ne. 3 .and. ftype .ne. 4) goto 91

!********************************************************
! FILE TYPE 2 (spatialy heterogeneous initial conditions)
!********************************************************
      if(ftype .eq. 2) then
!-------------------------------------------------------
! read second record (only if lmax > 0)
!-------------------------------------------------------             !

        if( lmax .gt. 1 ) then          !changed from 0 to 1 (5.3.2004) LMAX
          irec = 2
!       MODIFIED BY ALI AND PETRAS
!       CORPI, 16 July 2004

!
!
!       read(nb, err=90) (rlaux(l),l=1,lmax)
!       -----> read(nb, 5020, err=90) (rlaux(l),l=1,lmax)
!
          read(nb, 5020, err=90) (rlaux(l),l=1,lmax)
 5020   FORMAT(8F10.0)

          do l=1,lmax
            if( sed_depths(l) .ne. rlaux(l) ) goto 98
          end do
        else
          lmax = 1
        end if

!-------------------------------------------------------
! read data records
!-------------------------------------------------------

        irec = 3

        do ivar=1,nvar

!     MODIFIED BY ALI AND PETRAS
!     CORPI, 16 July 2004C
!     read(nb, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
! --> read(nb, 5030, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
!

        read(nb, 5030, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
 5030 FORMAT(F5.0)

      end do

!-------------------------------------------------------
! initialize the other levels if only surface is given
!-------------------------------------------------------

        do ivar=1,nvar
          do k=1,nkn
            val = var(lmax,k,ivar)
            do l=lmax+1,noslay
              var(l,k,ivar) = val
            end do
          end do
        end do

      end if ! end of ftype=2

!******************************************************
! FILE TYPE 1 (spatialy homogeneous initial conditions)
!******************************************************
      if(ftype .eq. 1) then
      irec=3
      do ivar=1,nvar
       read(nb, *, err=90) varval
       print *,'Initial condition: ','var= ',ivar,'value=', varval
       do l=1,lmax
        do k=1,nkn
         var(l,k,ivar)=varval
        end do
       end do
      end do


      end if !end of ftype 1

!******************************************************
! FILE TYPE 3 spatialy homogeneous initial conditions for different  levels
!******************************************************
      if(ftype .eq. 3) then
      irec=3
      do ivar=1,nvar
       read(nb, *, err=90) (varval_array(l),l=1,lmax)

       print *,'Initial condition: ','var= ',ivar,'value=',
     *       (varval_array(l),l=1,lmax)

       do l=1,lmax
        do k=1,nkn
         var(l,k,ivar) = varval_array(l)
        end do
       end do

      end do

      end if !end of type 3

!******************************************************
! FILE TYPE 4 dumped state for repeated runs
!******************************************************
      if(ftype .eq. 4) then

       do i = 1, nvar
          do j = 1, nkn
              read(nb,*,err=90) (var(k,j,i),k=1,lmax)
          end do
       end do

      end if !end of type 4
!-------------------------------------------------------
! reading done -> close file
!-------------------------------------------------------
 1001 continue
        close(nb)



!-------------------------------------------------------
! end of initialisation
!-------------------------------------------------------

        write(6,*) 'Succesfull initialization for ',name,' from file '
        write(6,*) file

        return

   90   continue
        write(6,*) 'read error in record = ',irec,' ivar = ',ivar
        write(6,*) '... reading file',file
        stop 'error stop inicfils_aquabc'
   91   continue
        write(6,*)
     +   'bad file type descriptor: value 1, 2 or 3 is allowed!'
        write(6,*) '... reading file',file
        stop 'error stop inicfils_aquabc'
   96   continue
        write(6,*) 'ivars not compatible with nvar: ',ivars,nvar
        stop 'error stop inicfil'
   97   continue
        write(6,*) 'Cannot open file ',file
        stop 'error stop inicfils_aquabc'
   98   continue
        write(6,*) 'levels are not the same from init file ',file
        write(6,*) (sed_depths(l),l=1,lmax)
        write(6,*) (rlaux(l),l=1,lmax)
        stop 'error stop inicfils_aquabc'
   99   continue
        write(6,*) 'INICFILS_AQUABC:' 
        print *, 'Parameters are not the same from init file ',file
        write(6,*) 'nkn, lmax from file  : ',nkk,lmax
        write(6,*) 'nkn, noslay from model : ',nkn,noslay
        stop 'error stop inicfils_aquabc'

       end  !inicfils_aquabc

!*******************************************************************
!*******************************************************************

        subroutine print_real_time(date,time)

c prints date and time, gets date and time

        implicit none
        character  date*8, time*10           !Added by Petras 11.09.2004


         write(6,*) '==============================================='
	     call date_and_time(date,time)
         write(6,*) 'DATE: ', date
         write(6,*) 'TIME: ', time
         write(6,*) '==============================================='
         end

!*******************************************************************
!*******************************************************************
      subroutine calc_run_time (time)

! prints runtime in human readable components

        implicit none

        real time, rem, secs
        integer days,hours,mins


        days = time/86400
        rem = (time - days*86400.)
        hours = rem/3600
        rem = rem - hours*3600
        mins = rem/60
        secs = rem - mins*60.


        write(6,'(" RUN TIME: ",i3,"d",i2,"h",i2,"m",f6.3,"s")')
     *            days,hours,mins,secs

      end !calc_run_time

!****************************************************************
!****************************************************************

      subroutine calc_time (start_array,end_array)

! Prints run time from data obtained by intrisic function date_and_time
! Order of elements in time arrays:
!      1 - year
!      2 - month of the year
!      3 - day of the month
!      4 - time offset with respect to UTC in minutes
!      5 - hour of the day
!      6 - minutes of the hour
!      7 - seconds of the minute
!      8 - milliseconds of the second

       implicit none
       integer start_array(8), end_array(8) ! arrays of date and time
                                            ! obtained from function date_and_time
       integer i
       integer date0 !reference time for conversion to seconds, hardcoded

       double precision   start_time, end_time
       real run_time
       integer  start_times, end_times   !number of full seconds

       integer s_year,s_month,s_day
       integer s_hour,s_mins,s_sec
       integer e_year,e_month,e_day
       integer e_hour,e_mins,e_sec

       character  date*8, time*6
       character  startdate*8, starttime*6
!
       s_year   = start_array(1)
       s_month  = start_array(2)
       s_day    = start_array(3)
       s_hour   = start_array(5)
       s_mins   = start_array(6)
       s_sec    = start_array(7)
!
       e_year   = end_array(1)
       e_month  = end_array(2)
       e_day    = end_array(3)
       e_hour   = end_array(5)
       e_mins   = end_array(6)
       e_sec    = end_array(7)

       date0 = (s_year*100 + 1)*100 + 1  ! Initializing reference time for SHYFEM
                                         ! time routine dts2it
       call dtsini(date0,0)              ! Initializing time routine


       call dts2it(start_times,s_year,s_month,s_day,
     *             s_hour,s_mins,s_sec)
       call dts2it(end_times,e_year,e_month,e_day,
     *               e_hour,e_mins,e_sec)

       start_time = dble(start_times) + 0.001*dble(start_array(8))
       end_time   = dble(end_times)   + 0.001*dble(end_array(8))

       run_time   = end_time - start_time


       write(6,*)'==============================='
       write(6,'(a8,i4.4,"-",i2.2,"-",i2.2," ",
     *           i2.2,":",i2.2,":",i2.2,":",i3.3)') ' START :',
     *           s_year, s_month, s_day,
     *           s_hour, s_mins, s_sec, start_array(8)
       write(6,'(a8,i4.4,"-",i2.2,"-",i2.2, " ",
     *           i2.2,":",i2.2,":",i2.2,":",i3.3)') ' FINISH:',
     *           e_year, e_month, e_day,
     *           e_hour, e_mins, e_sec, end_array(8)

       call calc_run_time (run_time)
       write(6,*)'==============================='

      end !calc_time



!******************************************************************
!******************************************************************
!******************************************************************

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     SUBROUTINE DEVELOPPED BY ALI AND PETRAS
!     CORPI, 9 August 2004
!
!     THIS SUBROUTINE INITIALIZES EUTRO TIME SERIES ASCII OUTPUTS for
!     given nodes (stations)
!
!     REVISIONS:
!
!     CORPI, 17 August 2004
!     CORPI, 31 August 2004, by Petras: Bug fix for case when output
!                                       is not required
!
!---------------------
!     Subroutine extended to read diagnostic output control information
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine biotser_init(nb,nstate,
     *                       NBIOTSMX,NBIOTS,
     *                       BIOTSNOD,BIOTSFUN,
     *                       NDGTSMX,NDIAGVAR,
     *                       NDGTS,DGTSNOD,DGTSFUN)


      integer nstate
      integer NDIAGVAR          !Maximum number of different intermediate variables in output
      integer NDGTSMX           !Total number of nodes for intermediates of each state variable
      integer NDGTS(nstate)     !Actual number of nodes for intermediates of each state variable
      integer DGTSFUN(NDGTSMX, nstate) !Output file units for EUTRO diagnostics time series(nodes) plots
      integer DGTSNOD(NDGTSMX, nstate) !Nodes for intermediate variables in output

      integer NBIOTSMX            !Total number of time series for EUTRO
      integer NBIOTS              !Actual number of time series for EUTRO
      integer BIOTSFUN(NBIOTSMX)  !Array keeping output file units for EUTRO time series plots
      integer BIOTSNOD(NBIOTSMX)  ! Array keeping nodes of EUTRO time series plots


      character*90 header

      character*80 BTSERNAME(NBIOTSMX)
      character*80 DTSERNAME(NDGTSMX,nstate)

      integer i
      integer itroub
      integer DODIAG
      integer DUMMY(NDGTSMX)
      logical berror

!     INITIALIZE
      itroub = 0

      do i= 1, NBIOTSMX
       BIOTSNOD(i) = 0
       BIOTSFUN(i) = 0
      end do

!     OPEN THE MAIN BIO TIME SERIES FILE

      if( nb .le. 0 ) goto 97

!      Read RECORD 1 - Read two lines
       read(nb, 5010, err=98) header
       read(nb, 5010, err=98) header
!      Read RECORD 2
       read(nb, *, err=99) NBIOTS, DODIAG


       if(NBIOTS.EQ.0) then
           write(6,*) 'NO NODES FOR TIME SERIES OUTPUT TO ASCII FILES'
           goto 1
       end if


      if(NBIOTS.LT.0) then
          write(6,*) 'NUMBER NODES FOR TIME SERIES OUTPUTS',
     +              'CANNOT BE NEGATIVE'
          stop 'error stop biotser_init'
      end if


       if (NBIOTS.GT.NBIOTSMX) then
           write(6,*) 'This executable image was compiled for ',
     +           NBIOTSMX, ' number of AQUABC time series outputs.'
           write(6,*) 'please change  NBIOTSMX in aquabc.h from'
     +           , NBIOTSMX, ' to ', NBIOTS, ' and recompile.'
           stop 'error stop biotser_init'
       end if

       write(6,*) 'AQUABC TIME SERIES OUTPUT TO ASCII ',
     +            'FILES WILL BE CREATED'


      write(6,*) 'NUMBER OF  NODES WITH TIME SERIES OUTPUTS : ', NBIOTS

!     Read RECORD 3

      do i=1, NBIOTS
       read(nb, *, err=100, end=200) BIOTSNOD(i), BTSERNAME(i)
       BIOTSFUN(i) = ifileo(70+i,trim(adjustl(BTSERNAME(i))),'f','u')
       if(BIOTSFUN(i).le.0) then
            write(6,*) 'Cannot open/create the AQUABC time series ',
     +                   ' output file ', i
              write(6,*) 'UNIT : ', BIOTSFUN(i)
              itroub = itroub + 1
       else
         write(6,*) 'AQUABC time series output file ', i, ' opened.'
         write(6,*) 'File unit : ', BIOTSFUN(i)
         write(6,*) 'File name : ', BTSERNAME(i)
       end if
      end do

      if(itroub.gt.0) then
          goto 101
      end if

      call n2int(NBIOTS, BIOTSNOD, berror)

      if(berror) then
          write(6,*) 'Problem in converting external nodes to internal'
          stop 'error stop: biotser_init'
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    1 continue
      if( DODIAG .EQ. 0) then
	    write(6,*) 'NO DIAGNOSTIC TIME SERIES OUTPUT TO ASCII FILES'
	    return
	else
	    write(6,*) 'DIAGNOSTIC TIME SERIES OUTPUT TO ASCII FILES ',
     +                'WILL BE CREATED'
	end if

C     Read RECORD 4 - Read two lines
      read(nb, 5010, err=102) header
      read(nb, 5010, err=102) header


      do j=1, nstate
C         Read RECORD 5
          read(nb, 5010, err=103) header
          write(*,*) 'RECORD 5 Read for state variable ', j

c          write(6,*)'bioser_init:', 'j=',j, 'NDGTS: ',NDGTS

C         Read RECORD 6
          read(nb, *, err=104) NDGTS(j)
          write(*,*) 'RECORD 6 Read for state variable ', j
	    if(NDGTS(j).LT.0) then
           write(6,*) 'NUMBER OF DIAGNOSTIC TIME SERIES OUTPUTS FOR ',
     +             'STATE VARIABLE ', j, ' IS ENTERED LESS THAN ZERO'
           write(6,*) 'NUMBER OF DIAGNOSTIC TIME SERIES OUTPUTS ',
     +		           'CANNOT BE NEGATIVE'
	        stop 'error stop biotser_init'
	    end if


	    if(NDGTS(j).GT.NBIOTSMX) then
              write(6,*) 'This executable image was compiled for ',
     +             NDGTSMX, ' number of EUTRO diagnostic time series ',
     +                      'outputs.'
              write(6,*) 'please change parameter NDGTSMX ',
     +                   'in AQUABC_AOUT.H ',
     +            'from', NDGTSMX, ' to ', NDGTS(j), ' and recompile.'
	        stop 'error stop biotser_init'
	    end if

        write(6,*)
     +  'NUMBER OF NODES FOR DIAGNOSTIC TIME SERIES OUTPUTS ',
     +               'FOR STATE VARIABLE ', j, ' : ', NDGTS(j)
	    if(NDGTS(j).GT.0) then
!             Read RECORD 7
	        do i=1, NDGTS(j)
		        read(nb, *, err=105) DGTSNOD(i,j), DTSERNAME(i,j)
          		DGTSFUN(i,j)=ifileo(70+(j*i),DTSERNAME(i,j),'f','u')
	            if(DGTSFUN(i,j).le.0) then
				    write(6,*) 'Cannot open/create EUTRO diagnostic',
     +                ' time series file ', i, ' for state variable ', j
				    write(6,*) 'UNIT : ', DGTSFUN(i,j)
				    itroub = itroub + 1
			    else
				    write(6,*) 'EUTRO diagnostic time series file ', i
     +			              ,' for state variable ', j, ' opened.'
                      write(6,*) 'File unit : ', DGTSFUN(i,j)
                      write(6,*) 'File name : ', DTSERNAME(i,j)

			    end if

                  DUMMY(i) = DGTSNOD(i,j)

	        end do

	        if(itroub.gt.0) then
	            goto 106
	        end if


!             CONVERT EXTERNAL NODES TO INTERNAL
              call n2int(NDGTS(j), DUMMY, berror)
	        if(berror) then
	            write(6,*) 'Problem in converting external nodes to ',
     +		           'internal, when processing state variable ', j
			    stop 'error stop: biotser_init'
		    end if
	        do i=1, NDGTS(j)
                  DGTSNOD(i,j) = DUMMY(i)
	        end do
	    end if
      end do
	close(nb)

 5010 FORMAT(A90)
 5020 FORMAT(2I10)
 5030 FORMAT(I10,A30)
 5040 FORMAT(I10)

      return

   97 continue
      write(6,*)
     + 'Cannot open AQUABC time series output information file'

	stop 'error stop biotser_init'
   98 continue
      write(6,*) 'Read error in RECORD 1 of AQUABC time series ',
     +           'ASCII output  information file '
	stop 'error stop biotser_init'
   99 continue
      write(6,*) 'Read error in RECORD 2 of AQUABC time series ',
     +           'ASCII output  information file '
	stop 'error stop biotser_init'
  100 continue
      write(6,*) 'Read error in RECORD 3 of AQUABC time series   ',
     +           'ASCII output  information file : ', i
	stop 'error stop biotser_init'
  101 continue
      write(6,*) 'One or several AQUABC time series ASCII output',
     +           '  files could not be opened/created.'
	stop 'error stop biotser_init'
  102 continue
      write(6,*) 'Read error in RECORD 4 of AQUABC time series',
     +           ' ASCII output  information file '
	stop 'error stop biotser_init'
  103 continue
      write(6,*) 'Read error in RECORD 5 of AQUABC time series',
     +           ' ASCII output  information file '
	stop 'error stop biotser_init'

  104 continue
      write(6,*) 'Read error in RECORD 6 of AQUABC time series ',
     +           'ASCII output  information file'
	stop 'error stop biotser_init'
  105 continue
      write(6,*) 'Read error in RECORD 7 of AQUABC time series ',
     +           'ASCII output  information file :  on line', i,
     +           'state variable ', j
	stop 'error stop biotser_init'

  106 continue
      write(6,*) 'One or several EUTRO diagnostic time ',
     +           'series output files could not be opened/created.'
	stop 'error stop biotser_init'
  200 continue
      write(6,*) 'ERROR, LESS TIME SERIES PLOT DATA AVAILABLE ',
     +            'THEN REQUIRED !!!'
      stop 'error stop biotser_init'
	return
	end !biotser_init

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

       subroutine biotser_write(initial,what,
     *                         e,noutput,nstate,dg,NDIAGVAR,
     *                         ilhkv,nlvd,
     *                         itmcon,idtcon,
     *                         NBIOTSMX,NBIOTS,BIOTSNOD,BIOTSFUN,
     *                         NDGTSMX, NDGTS, DGTSNOD, DGTSFUN)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Writes state variable(water column or sediments) and rates for defined
! nodes to text files
!
! Control information should be read from files by biotser_ini before first call of subroutine
!
! 2006 16 July   - Routine is started to be reworked  to process
!                  water column and sediment kinetic variables for 3D
! 2009   July    - Finished updates to process WC and BS variables for 3D
! 2011   August  - Universal formats for output
! 2011   October - Possibility to to convert seconds to date and time
!
! what     'wc' - write water column variables
!          'bs' - write bottom sediments variables. The only difference in algorithm: std output writing
! initial  1 - writes initial condition only (the second time step)
!          0 - writes everything  with given time step in *.str file
!
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use basin
      
      implicit none
      !include 'param.h'
      integer itanf,itend,idt,nits,niter,it
      common /femtim/ itanf,itend,idt,nits,niter,it

      real e(nlvd,nkndi,noutput)	!state variables
      real dg  (nlvd,nkndi,nstate,NDIAGVAR) !intermediate variables


	 integer nstate              !number of state variables
	 integer noutput            !number of variables for the output

	 integer NBIOTSMX            !Total number of time series for EUTRO
        integer NBIOTS              !Actual number of time series for EUTRO
	 integer BIOTSFUN(NBIOTSMX)  !Array keeping output file units for EUTRO time series plots
	 integer BIOTSNOD(NBIOTSMX)  !Array keeping nodes of EUTRO time series plots

        integer NDIAGVAR	  !Maximum number of intermediate variables displayed in diagnostics file
        integer NDGTSMX             !Total number of time series(nodes) for EUTRO
        integer NDGTS(nstate)       !Actual number of time series(nodes) for EUTROC
	 integer DGTSFUN(NDGTSMX,nstate) !Array keeping output file units for EUTRO diagnostics time series(nodes) plotsC
	 integer DGTSNOD(NDGTSMX,nstate) !Array keeping nodes of EUTRO diagnostics time series plots

	 integer ilhkv(1)            !number of layers for each node

      integer istate

      integer iend
      save iend
      data iend /0/

      integer itmcon, idtcon, nlvd
      integer i,j,k

      character*2 what !'wc' - for water column, 'bs' - for bottom sediments
!      format strings
      CHARACTER*30 FMT_1      !for one layer state variables seconds
      CHARACTER*50 FMT_1_d    !for one layer state variables day and time
      CHARACTER*30 FMT_many   !for many layers seconds
      CHARACTER*50 FMT_many_d !for many layers day and time

      CHARACTER*2 noutput_char
      CHARACTER*30 FMT_diag_1      !for one layer diagnostic variables seconds
      CHARACTER*50 FMT_diag_1_d    !for one layer diagnostic variables day and time
      CHARACTER*30 FMT_diag_many   !for many layers seconds
      CHARACTER*50 FMT_diag_many_d !for many layers day and time
      CHARACTER*2 NDIAGVAR_char

      integer initial

!     for conversion from seconds to dates and time
      integer date0, isec
      integer year,month, day, hour, min, sec
!-----------------------------------------------------------------------------

      date0 = 19970101  ! zero for time in seconds (00:00:00 is assumed)
      isec  = 0         ! output with date and time
!      isec =1           ! output in seconds

!     Initialize time routines
      call dtsini(date0,0)       !reference time 00:00:00 assumed
!     Convert seconds to date and time
      if(initial.eq.0) then
        call dts2dt(it,year,month,day,hour,min,sec)
      else
        call dts2dt(itanf,year,month,day,hour,min,sec)
      end if

!     Producing output format lines
      write(noutput_char,'(i2)') noutput
      FMT_many   = '(I15,I5,' // noutput_char   // '(1x,G13.4))'             !for many layers in seconds
      FMT_many_d = '(I4.4,2I2.2,1x,3I2.2,I5,'   //
     +                             noutput_char // '(1x,G13.4))'             !for many layers with date and time
      FMT_1      = '(I15,'    // noutput_char   // '(1x,G13.4))'             !for one layer state variables in seconds
      FMT_1_d    = '(I4.4,2I2.2,1x,3I2.2,'      //
     +                            noutput_char  // '(1x,G13.4))'             !for one layer state variables with date and time

      write(NDIAGVAR_char,'(i2)') NDIAGVAR

      FMT_diag_many   = '(I15,I5,' // NDIAGVAR_char //'(1x,G13.6))'           !for many layers in seconds
      FMT_diag_many_d = '(I4.4,2I2.2,1x,3I2.2,I5,'  //
     +                                NDIAGVAR_char //'(1x,G13.6))'            !for many layers with date and time
      FMT_diag_1      = '(I15,'    // NDIAGVAR_char //'(1x,G13.6))'            !for one layer diagnostic variables in seconds
      FMT_diag_1_d    = '(I4.4,2I2.2,1x,3I2.2,'     //
     +                                NDIAGVAR_char //'(1x,G13.6))'            !for one layer diagnostic variables with date and time



      if(initial.eq.1) goto 1000 ! writing initial conditions

      if( it .lt. itmcon ) then
       return
      end if

      !print *,'it=',it,'itmcon=',itmcon,'idtcon=',idtcon
      if( mod(it-itmcon,idtcon) .ne. 0 ) then
       return
      end if

1000   continue  ! going here if initial conditions writing only

       if((NBIOTS.EQ.0).or.(iend.eq.1)) then
        return
       end if

! ------------------Writing state variables to files-----------------------
       do i=1, NBIOTS


             if (ilhkv(BIOTSNOD(i)).le.1) then
              if(isec.eq.1) write(BIOTSFUN(i),FMT_1) it,
     *                           (e(1,BIOTSNOD(i),j),j=1,noutput)

              if(isec.eq.0) write(BIOTSFUN(i),FMT_1_d)
     *                            year,month,day,hour,min,sec,
     *                           (e(1,BIOTSNOD(i),j),j=1,noutput)

             else
              do k=1,ilhkv(BIOTSNOD(i))
               if(isec.eq.1) write(BIOTSFUN(i),FMT_many) it,k,
     *                               (e(k,BIOTSNOD(i),j),j=1,noutput)
               if(isec.eq.0) write(BIOTSFUN(i),FMT_many_d)
     *                                year,month,day,hour,min,sec,k,
     *                               (e(k,BIOTSNOD(i),j),j=1,noutput)
              enddo
             endif


          if((it+idtcon).gt.itend) then
            close(BIOTSFUN(i))
            iend = 1
          end if

       end do
! Write to the standard output
          write(6,*) 'BIOTSER_WRITE: ',what, ' variables ',
     *     ' for defined nodes',
     *     ' written to text files at ',it

      if (initial.eq.1) return ! Writing initial conditions only


! --------------Writing diagnostics (auxilary variables, rate components) to files----------

      do istate=1,nstate

       do i=1, NDGTS(istate)

        if (ilhkv(DGTSNOD(i,istate)).le.1) then
          if(isec.eq.1) write(DGTSFUN(i,istate), FMT_diag_1) it,
     +        (dg(1,DGTSNOD(i,istate),istate,j),j=1,NDIAGVAR)
          if(isec.eq.0) write(DGTSFUN(i,istate), FMT_diag_1_d)
     +                          year,month,day,hour,min,sec,
     +        (dg(1,DGTSNOD(i,istate),istate,j),j=1,NDIAGVAR)

        else
         do k=1,ilhkv(DGTSNOD(i,istate))
          if(isec.eq.1) write(DGTSFUN(i,istate),FMT_diag_many)
     *              it,k,
     *              (dg(k,DGTSNOD(i,istate),istate,j),
     *                                  j=1,NDIAGVAR)
          if(isec.eq.0) write(DGTSFUN(i,istate),FMT_diag_many_d)
     *                          year,month,day,hour,min,sec,k,
     *                      (dg(k,DGTSNOD(i,istate),istate,j),
     *                                  j=1,NDIAGVAR)
         enddo
        endif



         if((it+idtcon).gt.itend) then
            close(DGTSFUN(i,istate))
            iend = 1
         end if

       end do !i
      end do  !istate

! Write to the standard output
         write(6,*) 'BIOTSER_WRITE:', what, ' rates ',
     +           ' for defined nodes written to text files at ',it
!----------------------------------------------------------------------------------

       return
      end !biotser_write

!********************************************************************
!********************************************************************
!********************************************************************

      subroutine cur_param_read_ini_wc(bname,nconst)

!     initialize WC model constants names array
!     water column kinetics parameters
!     Input/Output:
!      bname  -  names array
!      nconst  -  number of parameters required
!      Note:
!           Number of parameters in this routine is hardcoded by npar
!      Called:
!           In CUR_PARAM_READ

      integer npar !existing number of parameters in the routine
      parameter (npar = 291)
      character*80 bname(npar)
      real par(npar)

!     execution part

      if( npar. ne. nconst) then
        print *, 'CUR_PARAM_READ_INI_WC:'
        print *, '  wrong required number of params'
        print *, '  required: ',nconst
        print *, '  existing: ',npar
        stop
      end if

!     WC kinetics parameter names. These line should be corrected
!     when new parameter is introduced or the name of existing
!     parameter is changed.
      bname(1  ) =                              'K_A'  !1  ! Aeration coefficient (if negative calculates internally)
      bname(2  ) =                        'THETA_K_A'  !2  ! Temperature correction factor for aeration
      bname(3  ) =                               'KE'  !3  ! Background extinction coefficient
      bname(4  ) =                              'XKC'  !4  ! Light extinction per chlorophyl unit,( mcg Chla/l/m)
      bname(5  ) =                            'PHIMX'  !5  ! Quantum yield const. mg C/mole photon
      
!       bname(6  ) =               'KG_CHEM_AUT_BAC_20'  !6  ! Chemoautotrophic bacteria Growth rate
!       bname(7  ) =          'EFF_CHEM_AUT_BAC_GROWTH'  !7  ! Chemoautotrophic bacteria growth efficiency
!       bname(8  ) =            'THETA_KG_CHEM_AUT_BAC'  !8  ! Chemoautotrophic bacteria Temperature correction for growth rate
!       bname(9  ) =               'KR_CHEM_AUT_BAC_20'  !9  ! Chemoautotrophic bacteria Respiration rate
!       bname(10 ) =            'THETA_KR_CHEM_AUT_BAC'  !10 ! Chemoautotrophic bacteria Temperature correction for respiration rate
!       bname(11 ) =               'KD_CHEM_AUT_BAC_20'  !11 ! Chemoautotrophic bacteria Mortality rate
!       bname(12 ) =            'THETA_KD_CHEM_AUT_BAC'  !12 ! Chemoautotrophic bacteria Temperature correction for Mortality rate
!       bname(13 ) =            'KHS_NH4N_CHEM_AUT_BAC'  !13 ! Chemoautotrophic bacteria Half saturation growth for NH4N
!       bname(14 ) =            'KHS_PO4P_CHEM_AUT_BAC'  !14 ! Chemoautotrophic bacteria Half saturation growth for PO4P
!       bname(15 ) =              'KHS_O2_CHEM_AUT_BAC'  !15 ! Chemoautotrophic bacteria Half saturation growth for O2
!       bname(16 ) =      'DO_STR_HYPOX_CHEM_AUT_BAC_D'  !16 ! Chemoautotrophic bacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!       bname(17 ) =       'THETA_HYPOX_CHEM_AUT_BAC_D'  !17 ! Chemoautotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
!       bname(18 ) =       'EXPON_HYPOX_CHEM_AUT_BAC_D'  !18 ! Chemoautotrophic bacteria Exponent constant for Dissolved oxygen stress
!       bname(19 ) =              'CHEM_AUT_BAC_N_TO_C'  !19 ! Chemoautotrophic bacteria Nitrogen to Carbon ratio
!       bname(20 ) =              'CHEM_AUT_BAC_P_TO_C'  !20 ! Chemoautotrophic bacteria Phosphorus to Carbon ratio
!       bname(21 ) =             'CHEM_AUT_BAC_O2_TO_C'  !21 ! Chemoautotrophic bacteria Oxygen to Carbon ratio
!       bname(22 ) =               'YIELD_CHEM_AUT_BAC'  !22 ! Chemoautotrophic bacteria Yield of Carbon per unit Nitrates nitrogen
!       bname(23 ) =                'KG_AER_HET_BAC_20'  !23 ! Aerobic heterotrophic bacteria Growth rate
!       bname(24 ) =           'EFF_AER_HET_BAC_GROWTH'  !24 ! Aerobic heterotrophic bacteria growth efficiency
!       bname(25 ) =             'THETA_KG_AER_HET_BAC'  !25 ! Aerobic heterotrophic bacteria Temperature correction for growth rate
!       bname(26 ) =                'KR_AER_HET_BAC_20'  !26 ! Aerobic heterotrophic bacteria Respiration rate
!       bname(27 ) =             'THETA_KR_AER_HET_BAC'  !27 ! Aerobic heterotrophic bacteria Temperature correction for respiration rate
!       bname(28 ) =                'KD_AER_HET_BAC_20'  !28 ! Aerobic heterotrophic bacteria Mortality rate
!       bname(29 ) =             'THETA_KD_AER_HET_BAC'  !29 ! Aerobic heterotrophic bacteria Temperature correction for Mortality rate
!       bname(30 ) =             'KHS_ORGC_AER_HET_BAC'  !30 ! Aerobic heterotrophic bacteria Half saturation growth for OC
!       bname(31 ) =              'KHS_ORGN_AER_HET_BAC' !31 ! Aerobic heterotrophic bacteria Half saturation growth for ON
!       bname(32 ) =              'KHS_ORGP_AER_HET_BAC' !32 ! Aerobic heterotrophic bacteria Half saturation growth for OP
!       bname(33 ) =                'KHS_O2_AER_HET_BAC' !33 ! Aerobic heterotrophic bacteria Half saturation growth for Oxygen
!       bname(34 ) =               'KHS_DIN_AER_HET_BAC' !34 ! Aerobic heterotrophic bacteria Half saturation growth for inorganic nitrogen
!       bname(35 ) =               'KHS_DIP_AER_HET_BAC' !35 ! Aerobic heterotrophic bacteria Half saturation growth for inorganic phosphorus
!       bname(36 ) =              'KHS_PHYT_AER_HET_BAC' !36 ! Aerobic heterotrophic bacteria Half saturation growth for Phytoplankton C (not used as a resource)
!       bname(37 ) =              'YIELD_OC_AER_HET_BAC' !37 ! Aerobic heterotrophic bacteria Yield of bacteria carbon per unit of organic carbon
!       bname(38 ) =               'OX_ORGN_AER_HET_BAC' !38 ! Aerobic heterotrophic bacteria ON oxidation rate mg N per mg C of bacteria production
!       bname(39 ) =                        'KHS_MIN_N'  !39 ! Aerobic heterotrophic bacteria ON mineralisation reverse half saturation for DIN
!       bname(40 ) =              'OX_ORGP_AER_HET_BAC'  !40 ! Aerobic heterotrophic bacteria OP mineralisation rate mg P per mg C of bacteria production
!       bname(41 ) =                        'KHS_MIN_P'  !41 ! Aerobic heterotrophic bacteria ON mineralisation reverse half saturation for DIP
!       bname(42 ) =       'DO_STR_HYPOX_AER_HET_BAC_D'  !42 ! Aerobic heterotrophic bacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
!       bname(43 ) =        'THETA_HYPOX_AER_HET_BAC_D'  !43 ! Aerobic heterotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
!       bname(44 ) =        'EXPON_HYPOX_AER_HET_BAC_D'  !44 ! Aerobic heterotrophic bacteria Exponent constant for Dissolved oxygen stress
!       bname(45 ) =               'AER_HET_BAC_N_TO_C'  !45 ! Aerobic heterotrophic bacteria Nitrogen to Carbon ratio
!       bname(46 ) =               'AER_HET_BAC_P_TO_C'  !46 ! Aerobic heterotrophic bacteria Phosphorus to Carbon ratio
!       bname(47 ) =              'AER_HET_BAC_O2_TO_C'  !47 ! Aerobic heterotrophic bacteria Oxygen to Carbon ratio for respiration
!       bname(48 ) =             'KG_FAC_AN_HET_BAC_20'  !48 ! Facultative anaerobic heterotrophic bacteria Growth rate of
!       bname(49 ) =        'EFF_FAC_AN_HET_BAC_GROWTH'  !49 ! not used! Facultative anaerobic heterotrophic bacteria growth efficiency
!       bname(50 ) =          'THETA_KG_FAC_AN_HET_BAC'  !50 ! not used! Facultative anaerobic heterotrophic bacteria Temperature correction for growth rate
!       bname(51 ) =             'KR_FAC_AN_HET_BAC_20'  !51 ! not used! Facultative anaerobic heterotrophic bacteria Respiration rate
!       bname(52 ) =          'THETA_KR_FAC_AN_HET_BAC'  !52 ! not used! Facultative anaerobic heterotrophic bacteria Temperature correction for respiration rate
!       bname(53 ) =             'KD_FAC_AN_HET_BAC_20'  !53 ! not used! Facultative anaerobic heterotrophic bacteria Mortality rate
!       bname(54 ) =          'THETA_KD_FAC_AN_HET_BAC'  !54 ! not used! Facultative anaerobic heterotrophic bacteria Temperature correction for Mortality rate
!       bname(55 ) =          'KHS_NO3N_FAC_AN_HET_BAC'  !55 ! Facultative anaerobic heterotrophic bacteria Half saturation growth for NO3N
!       bname(56 ) =          'KHS_ORGC_FAC_AN_HET_BAC'  !56 ! not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for OC
!       bname(57 ) =          'KHS_ORGN_FAC_AN_HET_BAC'  !57 ! not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for ON
!       bname(58 ) =          'KHS_ORGP_FAC_AN_HET_BAC'  !58 ! not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for OP
!       bname(59 ) =        'REV_KHS_O2_FAC_AN_HET_BAC'  !59 ! not used! Facultative anaerobic heterotrophic bacteria Reverse Half saturation growth for O2
!       bname(60 ) =   'NO3N_LACK_STR_FAC_AN_HET_BAC_D'  !60 ! not used! Facultative anaerobic heterotrophic bacteria NO3N stress concentration
!       bname(61 ) =  'THETA_NO3_LACK_FAC_AN_HET_BAC_D'  !61 ! not used! Facultative anaerobic heterotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
!       bname(62 ) =    'EXP_NO3_LACK_FAC_AN_HET_BAC_D'  !62 ! not used! Facultative anaerobic heterotrophic bacteria Exponent constant for Dissolved oxygen stress
!       bname(63 ) =            'FAC_AN_HET_BAC_N_TO_C'  !63 ! not used! Facultative anaerobic heterotrophic bacteria Nitrogen to Carbon ratio
!       bname(64 ) =            'FAC_AN_HET_BAC_P_TO_C'  !64 ! not used! Facultative anaerobic heterotrophic bacteria Phosphorus to Carbon ratio
!       bname(65 ) =           'FAC_AN_HET_BAC_O2_TO_C'  !65 ! not used! Facultative anaerobic heterotrophic bacteria Oxygen to Carbon ratio for respiration
!       bname(66 ) =             'YIELD_FAC_AN_HET_BAC'  !66 ! Facultative anaerobic heterotrophic bacteria Yield of carbon per unit nitrates nitrogen

      bname(6  ) =                  'KG_DIA_OPT_TEMP'  !67 ! Diatoms Growth rate
      bname(7  ) =                  'DIA_OPT_TEMP_LR'  !68 ! Diatoms optimal temperature lower range
      bname(8  ) =                  'DIA_OPT_TEMP_UR'  !69 ! Diatoms optimal temperature upper range
      bname(9  ) =                   'EFF_DIA_GROWTH'  !70 ! Diatoms Effective growth. (1-EG)*growth - losses for respiration and excretion
      bname(10 ) =         'KAPPA_DIA_UNDER_OPT_TEMP'  !71 ! Diatoms Temperature correction for growth lower temperature
      bname(11 ) =          'KAPPA_DIA_OVER_OPT_TEMP'  !72 ! Diatoms Temperature correction for growth upper temperature
      bname(12 ) =                        'KR_DIA_20'  !73 ! Diatoms Respiration rate
      bname(13 ) =                     'THETA_KR_DIA'  !74 ! Diatoms Temperature correction for basal respiration rate
      bname(14 ) =                        'KD_DIA_20'  !75 ! Diatoms Mortality rate
      bname(15 ) =                     'THETA_KD_DIA'  !76 ! Diatoms Temperature correction for Mortality rate
      bname(16 ) =                      'KHS_DIN_DIA'  !77 ! Diatoms Half saturation growth for DIN
      bname(17 ) =                      'KHS_DIP_DIA'  !78 ! Diatoms Half saturation growth for DIP
      bname(18 ) =                      'KHS_DSi_DIA'  !79 ! Diatoms Half saturation growth for DSi
      bname(19 ) =                       'KHS_O2_DIA'  !80 ! Diatoms Half saturation growth for O2
      bname(20 ) =                'FRAC_DIA_EXCR'      !81 ! Diatoms Fraction of excretion in metabolism rate
      bname(21) =                          'I_S_DIA'  !82 ! Diatoms Light saturation (langleys)
      bname(22) =               'DO_STR_HYPOX_DIA_D'  !83 ! Diatoms Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
      bname(23) =                'THETA_HYPOX_DIA_D'  !84 ! Diatoms Multiplier of the exponent for Dissolved oxygen stress
      bname(24) =                'EXPON_HYPOX_DIA_D'  !85 ! Diatoms Exponent constant for Dissolved oxygen stress
      bname(25) =                       'DIA_N_TO_C'  !86 ! Diatoms Nitrogen to Carbon ratio
      bname(26) =                       'DIA_P_TO_C'  !87 ! Diatoms Phosphorus to Carbon ratio
      bname(27) =                      'DIA_Si_TO_C'  !88 ! Diatoms Silica to Carbon ratio
      bname(28) =                      'DIA_O2_TO_C'  !89 ! Diatoms Oxygen to Carbon ratio for respiration
      bname(29) =                    'DIA_C_TO_CHLA'  !90 ! Diatoms Carbon to Chlorophil a ratio
      bname(30) =                  'KG_CYN_OPT_TEMP'  !91 ! Non-fixing cyanobacteria Growth rate
      bname(31) =                  'CYN_OPT_TEMP_LR'  !92 ! Non-fixing cyanobacteria optimal temperature lower range
      bname(32) =                  'CYN_OPT_TEMP_UR'  !93 ! Non-fixing cyanobacteria optimal temperature upper range
      bname(33) =                   'EFF_CYN_GROWTH'  !94 ! Non-fixing cyanobacteria Effective growth. (1-EG)*growth - losses for respiration and excretion
      bname(34) =         'KAPPA_CYN_UNDER_OPT_TEMP'  !95 ! Non-fixing cyanobacteria Temperature correction for growth lower temperature
      bname(35) =          'KAPPA_CYN_OVER_OPT_TEMP'  !96 ! Non-fixing cyanobacteria Temperature correction for growth upper temperature
      bname(36) =                        'KR_CYN_20'  !97 ! Non-fixing cyanobacteria Respiration rate
      bname(37) =                     'THETA_KR_CYN'  !98 ! Non-fixing cyanobacteria Temperature correction for respiration rate
      bname(38) =                        'KD_CYN_20'  !99 ! Non-fixing cyanobacteria Mortality rate
      bname(39) =                     'THETA_KD_CYN'  !100! Non-fixing cyanobacteria Temperature correction for Mortality rate
      bname(40) =                      'KHS_DIN_CYN'  !101! Non-fixing cyanobacteria Half saturation growth for DIN
      bname(41) =                      'KHS_DIP_CYN'  !102! Non-fixing cyanobacteria Half saturation growth for DIP
      bname(42) =                       'KHS_O2_CYN'  !103! Non-fixing cyanobacteria Half saturation growth for O2
      bname(43) =                'FRAC_CYN_EXCR'      !104! Non-fixing cyanobacteria Fraction of excretion in metabolism rate
      bname(44) =                          'I_S_CYN'  !105! Non-fixing cyanobacteria Light saturation (langleys)
      bname(45) =               'DO_STR_HYPOX_CYN_D'  !106! Non-fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
      bname(46) =                'THETA_HYPOX_CYN_D'  !107! Non-fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
      bname(47) =                'EXPON_HYPOX_CYN_D'  !108! Non-fixing cyanobacteria Exponent constant for Dissolved oxygen stress
      bname(48) =                       'CYN_N_TO_C'  !109! Non-fixing cyanobacteria Nitrogen to Carbon ratio ,was 0.1
      bname(49) =                       'CYN_P_TO_C'  !110! Non-fixing cyanobacteria Phosphorus to Carbon ratio
      bname(50) =                      'CYN_O2_TO_C'  !111! Non-fixing cyanobacteria Oxygen to Carbon ratio for respiration
      bname(51) =                    'CYN_C_TO_CHLA'  !112! Non-fixing cyanobacteria Carbon to Chlorophyl a ratio
      bname(52) =              'KG_FIX_CYN_OPT_TEMP'  !113! Fixing cyanobacteria Growth rate
      bname(53) =              'FIX_CYN_OPT_TEMP_LR'  !114! Fixing Cyanobacteria optimal temperature lower range
      bname(54) =              'FIX_CYN_OPT_TEMP_UR'  !115! Fixing Cyanobacteria optimal temperature upper range
      bname(55) =               'EFF_FIX_CYN_GROWTH'  !116! Fixing cyanobacteria Effective growth. (1-EG)*growth - losses for RESP and excretion
      bname(56) =     'KAPPA_FIX_CYN_UNDER_OPT_TEMP'  !117! Fixing cyanobacteria Temperature correction for growth lower temperature
      bname(57) =      'KAPPA_FIX_CYN_OVER_OPT_TEMP'  !118! Fixing cyanobacteria Temperature correction for growth upper temperature
      bname(58) =                    'KR_FIX_CYN_20'  !119! Fixing cyanobacteria RESP rate
      bname(59) =                 'THETA_KR_FIX_CYN'  !120! Fixing cyanobacteria Temperature correction for RESP rate
      bname(60) =                    'KD_FIX_CYN_20'  !121! Fixing cyanobacteria Mortality rate of nitrification bacteria
      bname(61) =                 'THETA_KD_FIX_CYN'  !122! Fixing cyanobacteria Temperature correction for Mortality rate
      bname(62) =                  'KHS_DIN_FIX_CYN'  !123! Fixing cyanobacteria Half saturation growth for DIN
      bname(63) =                  'KHS_DIP_FIX_CYN'  !124! Fixing cyanobacteria Half saturation growth for DIP
      bname(64) =                   'KHS_O2_FIX_CYN'  !125! Fixing cyanobacteria Half saturation growth for O2
      bname(65) =            'FRAC_FIX_CYN_EXCR'      !126! Fixing cyanobacteria Fraction of excretion in metabolism rate
      bname(66) =                      'I_S_FIX_CYN'  !127! Fixing cyanobacteria Light saturation (langleys)
      bname(67) =           'DO_STR_HYPOX_FIX_CYN_D'  !128! Fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
      bname(68) =            'THETA_HYPOX_FIX_CYN_D'  !129! Fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
      bname(69) =            'EXPON_HYPOX_FIX_CYN_D'  !130! Fixing cyanobacteria Exponent constant for Dissolved oxygen stress
      bname(70) =                   'FIX_CYN_N_TO_C'  !131! Fixing cyanobacteria Nitrogen to Carbon ratio
      bname(71) =                   'FIX_CYN_P_TO_C'  !132! Fixing cyanobacteria Phosphorus to Carbon ratio
      bname(72) =                  'FIX_CYN_O2_TO_C'  !133! Fixing cyanobacteria Oxygen to Carbon ratio for respiration
      bname(73) =                'FIX_CYN_C_TO_CHLA'  !134! Fixing cyanobacteria Carbon to Chlorophyl a ratio
      bname(74) =                            'R_FIX'  !135! Fixing cyanobacteria Ratio between non-fixing and fixing fractions growth rate
      bname(75) =                            'K_FIX'  !136! Fixing cyanobacteria Effectivity parameter of switching to nitrogen fixation
      bname(76) =                  'KG_OPA_OPT_TEMP'  !137! OtherPhyto Growth rate
      bname(77) =                  'OPA_OPT_TEMP_LR'  !138! OtherPhyto optimal temperature lower range
      bname(78) =                  'OPA_OPT_TEMP_UR'  !139! OtherPhyto optimal temperature upper range
      bname(79) =                   'EFF_OPA_GROWTH'  !140! OtherPhyto Effective growth. (1-EG)*growth - losses for respiration and excretion
      bname(80) =         'KAPPA_OPA_UNDER_OPT_TEMP'  !141! OtherPhyto Temperature correction for growth lower temperature
      bname(81 ) =          'KAPPA_OPA_OVER_OPT_TEMP'  !142! OtherPhyto Temperature correction for growth upper temperature
      bname(82 ) =                        'KR_OPA_20'  !143! OtherPhyto Respiration rate
      bname(83 ) =                     'THETA_KR_OPA'  !144! OtherPhyto Temperature correction for respiration rate
      bname(84 ) =                        'KD_OPA_20'  !145! OtherPhyto Mortality rate
      bname(85 ) =                     'THETA_KD_OPA'  !146! OtherPhyto Temperature correction for Mortality rate
      bname(86 ) =                      'KHS_DIN_OPA'  !147! OtherPhyto Half saturation growth for DIN
      bname(87 ) =                      'KHS_DIP_OPA'  !148! OtherPhyto Half saturation growth for DIP
      bname(88 ) =                       'KHS_O2_OPA'  !149! OtherPhyto Half saturation growth for O2
      bname(89 ) =                'FRAC_OPA_EXCR'      !150! OtherPhyto Fraction of excretion in metabolism rate
      bname(90 ) =                          'I_S_OPA'  !151! OtherPhyto Light saturation (langleys)
      bname(91 ) =               'DO_STR_HYPOX_OPA_D'  !152! OtherPhyto Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
      bname(92 ) =                'THETA_HYPOX_OPA_D'  !153! OtherPhyto Multiplier of the exponent for Dissolved oxygen stress
      bname(93 ) =                'EXPON_HYPOX_OPA_D'  !154! OtherPhyto Exponent constant for Dissolved oxygen stress
      bname(94 ) =                       'OPA_N_TO_C'  !155! OtherPhyto Nitrogen to Carbon ratio
      bname(95 ) =                       'OPA_P_TO_C'  !156! OtherPhyto Phosphorus to Carbon ratio
      bname(96 ) =                      'OPA_O2_TO_C'  !157! OtherPhyto Oxygen to Carbon ratio for respiration
      bname(97 ) =                    'OPA_C_TO_CHLA'  !158! OtherPhyto Carbon to Chlorophyl a ratio
      bname(98 ) =                  'KG_ZOO_OPT_TEMP'  !159! Zooplankton Growth rate
      bname(99 ) =                  'ZOO_OPT_TEMP_LR'  !160! Zooplankton optimal temperature lower range
      bname(100) =                  'ZOO_OPT_TEMP_UR'  !161! Zooplankton optimal temperature upper range
      bname(101) =                   'EFF_ZOO_GROWTH'  !162! Zooplankton Effective growth. (1-EG)*growth - losses for respiration and excretion
      bname(102) =         'KAPPA_ZOO_UNDER_OPT_TEMP'  !163! Zooplankton Temperature correction for growth lower temperature
      bname(103) =          'KAPPA_ZOO_OVER_OPT_TEMP'  !164! Zooplankton Temperature correction for growth upper temperature
      bname(104) =                     'GRAT_ZOO_DIA'  !165! Zooplankton Grazing rate (growhth rate multiplier) on diatoms
      bname(105) =                     'GRAT_ZOO_CYN'  !166! Zooplankton Grazing rate (growhth rate multiplier) on Cyanobacteria
      bname(106) =                     'GRAT_ZOO_OPA'  !167! Zooplankton Grazing rate (growhth rate multiplier) on fixing Cyanobacteria
      bname(107) =                 'GRAT_ZOO_FIX_CYN'  !168! Zooplankton Grazing rate (growhth rate multiplier) on OtherPhyto
!       bname(108) =            'GRAT_ZOO_CHEM_AUT_BAC'  !169! Zooplankton Grazing rate (growhth rate multiplier) on NITR_BAC
!       bname(109) =             'GRAT_ZOO_AER_HET_BAC'  !170! Zooplankton Grazing rate (growhth rate multiplier) on AER_HET_BAC
!       bname(110) =          'GRAT_ZOO_FAC_AN_HET_BAC'  !171! Zooplankton Grazing rate (growhth rate multiplier) on DENITR_BAC
      bname(108) =          'GRAT_ZOO_DET_PART_ORG_C'  !172! Zooplankton Grazing rate (growhth rate multiplier) on part. ORG_C
      bname(109) =                     'PREF_ZOO_DIA'  !173! Zooplankton Preference for diatoms
      bname(110) =                     'PREF_ZOO_CYN'  !174! Zooplankton Preference for Cyanobacteria
      bname(111) =                 'PREF_ZOO_FIX_CYN'  !175! Zooplankton Preference for fixing Cyanobacteria
      bname(112) =                     'PREF_ZOO_OPA'  !176! Zooplankton Preference for OtherPhyto
!       bname(116) =            'PREF_ZOO_CHEM_AUT_BAC'  !177! Zooplankton Preference for NITR_BAC
!       bname(117) =             'PREF_ZOO_AER_HET_BAC'  !178! Zooplankton Preference for AER_HET_BAC
!       bname(118) =          'PREF_ZOO_FAC_AN_HET_BAC'  !179! Zooplankton Preference for DENITR_BAC
      bname(113) =          'PREF_ZOO_DET_PART_ORG_C'  !180! Zooplankton Preference for part. ORG_C
      bname(114) =                    'KHS_DIA_C_ZOO'  !181! Zooplankton Half saturation growth for diatoms
      bname(115) =                    'KHS_CYN_C_ZOO'  !182! Zooplankton Half saturation growth for Cyanobacteria
      bname(116) =                'KHS_FIX_CYN_C_ZOO'  !183! Zooplankton Half saturation growth for fixing Cyanobacteria
      bname(117) =                    'KHS_OPA_C_ZOO'  !184! Zooplankton Half saturation growth for OtherPhyto
!       bname(124) =           'KHS_CHEM_AUT_BAC_C_ZOO'  !185! Zooplankton Half saturation growth for NITR_BAC
!       bname(125) =            'KHS_AER_HET_BAC_C_ZOO'  !186! Zooplankton Half saturation growth for AER_HET_BAC
!       bname(126) =         'KHS_FAC_AN_HET_BAC_C_ZOO'  !187! Zooplankton Half saturation growth for DENITR_BAC
      bname(118) =           'KHS_DET_PART_ORG_C_ZOO'  !188! Zooplankton Half saturation growth for part. ORG_C
      bname(119) =                     'FOOD_MIN_ZOO'  !189! Zooplankton Minimum food conc. for feeding
      bname(120) =                           'KE_ZOO'  !190! not used Zooplankton Excretion rate as growth fraction
      bname(121) =                  'FRAC_ZOO_EX_ORG'  !191! not used Zooplankton Excretion rate organic fraction
      bname(122) =                        'KR_ZOO_20'  !192! Zooplankton Respiration rate
      bname(123) =                     'THETA_KR_ZOO'  !193! Zooplankton Respiration rate Temperature correction
      bname(124) =                        'KD_ZOO_20'  !194! Zooplankton Mortality rate
      bname(125) =                     'THETA_KD_ZOO'  !195! Zooplankton Mortality rate Temperature correction
      bname(126) =               'DO_STR_HYPOX_ZOO_D'  !196! Zooplankton Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
      bname(127) =                'THETA_HYPOX_ZOO_D'  !197! Zooplankton Multiplier of the exponent for Dissolved oxygen stress
      bname(128) =                'EXPON_HYPOX_ZOO_D'  !198! Zooplankton Exponent constant for Dissolved oxygen stress
      bname(129) =                       'ZOO_N_TO_C'  !199! Zooplankton Nitrogen to Carbon ratio
      bname(130) =                       'ZOO_P_TO_C'  !200! Zooplankton Phosphorus to Carbon ratio
      bname(131) =                      'ZOO_O2_TO_C'  !201! Zooplankton Oxygen to Carbon ratio for respiration
      bname(132) =          'KDISS_DET_PART_ORG_C_20'  !202! Particulate Detritus Carbon Dissolution rate not dependent on phytoplankton
      bname(133) =       'THETA_KDISS_DET_PART_ORG_C'  !203! Particulate Detritus Carbon Dissolution rate Temperature correction
      bname(134) =          'FAC_PHYT_DET_PART_ORG_C'  !204! Particulate Detritus Carbon Phytoplankton linear factor for dissolution rate
      bname(135) =          'KDISS_DET_PART_ORG_N_20'  !205! Particulate Detritus Nitrogen Dissolution rate not dependent on phytoplankton
      bname(136) =       'THETA_KDISS_DET_PART_ORG_N'  !206! Particulate Detritus Nitrogen Dissolution rate Temperature correction
      bname(137) =                       'KHS_DISS_N'  !207! Particulate Detritus Nitrogen dissolution reverse half saturation for DIN
      bname(138) =          'FAC_PHYT_DET_PART_ORG_N'  !208! Particulate Detritus Nitrogen Phytoplankton linear factor for dissolution rate
      bname(139) =          'KDISS_DET_PART_ORG_P_20'  !209! Particulate Detritus Phosphorus Dissolution rate not dependent on phytoplankton
      bname(140) =       'THETA_KDISS_DET_PART_ORG_P'  !210! Particulate Detritus Phosphorus Dissolution rate Temperature correction
      bname(141) =                       'KHS_DISS_P'  !211! Particulate Detritus Phosphorus  dissolution reverse half saturation for DIP
      bname(142) =          'FAC_PHYT_DET_PART_ORG_P'  !212! Particulate Detritus Phosphorus  Phytoplankton linear factor for dissolution rate
      bname(143) =                 'KDISS_PART_Si_20'  !213! Particulate Silica Dissolution rate
      bname(144) =              'THETA_KDISS_PART_Si'  !214! Particulate Silica Dissolution rate Temperature correction
      bname(145) =                     'K_MIN_DOC_20'  !215! Dissolved carbon  mineralisation rate
      bname(146) =                  'THETA_K_MIN_DOC'  !216! Dissolved carbon  mineralisation rate Temperature constant
      bname(147) =                'FAC_PHYT_AMIN_DOC'  !217! Dissolved carbon  Phytoplankton linear factor for mineralisation rate
      bname(148) =                     'K_MIN_DON_20'  !218! Dissolved nitrogen  mineralisation rate not dependent on phytoplankton
      bname(149) =                  'THETA_K_MIN_DON'  !219! Dissolved nitrogen  mineralisation rate Temperature constant
      bname(150) =                       'KHS_AMIN_N'  !220! Dissolved nitrogen  reverse half saturation for DIN
      bname(151) =                'FAC_PHYT_AMIN_DON'  !221! Dissolved nitrogen Phytoplankton linear factor for mineralisation rate
      bname(152) =                     'K_MIN_DOP_20'  !222! Dissolved phosphorus  mineralisation rate not dependent on phytoplankton
      bname(153) =                  'THETA_K_MIN_DOP'  !223! Dissolved phosphorus  mineralisation rate Temperature constant
      bname(154) =                       'KHS_AMIN_P'  !224! Dissolved phosphorus reverse half saturation for DIP
      bname(155) =                'FAC_PHYT_AMIN_DOP'  !225! Dissolved phosphorus Phytoplankton linear factor for mineralisation rate
      bname(156) =                        'K_NITR_20'  !226! Amonia nitrification rate
      bname(157) =                     'KHS_NITR_OXY'  !227! Amonia nitrification half saturation for Oxygen
      bname(158) =                   'KHS_NITR_NH4_N'  !228! Amonia nitrification half saturation for Amonia
      bname(159) =                     'THETA_K_NITR'  !229! Amonia nitrification rate Temperature constant      
      bname(160) =                   'PH_MIN_DOC_MIN'       !230!   optimum lower range for pH correction factor for DOC
      bname(161) =                   'PH_MIN_DOC_MAX'       !231!   optimum upper range for pH correction factor for DOC
      bname(162) =                   'PH_MIN_DON_MIN'       !232!   optimum lower range for pH correction factor for DON
      bname(163) =                   'PH_MIN_DON_MAX'       !233!   optimum upper range for pH correction factor for DON
      bname(164) =                   'PH_MIN_DOP_MIN'       !234!   optimum lower range for pH correction factor for DOP
      bname(165) =                   'PH_MIN_DOP_MAX'       !235!   optimum upper range for pH correction factor for DOP
      bname(166) =                   'PH_NITR_NH4_MIN'      !236!   optimum lower range for pH correction factor for nitrification
      bname(167) =                   'PH_NITR_NH4_MAX'      !237!   optimum upper range for pH correction factor for nitrification
      bname(168) =                   'PH_DENITR_NO3_MIN'    !238!   optimum lower range for pH correction factor for denitrification
      bname(169) =                   'PH_DENITR_NO3_MAX'    !239!   optimum upper range for pH correction factor for denitrification
      bname(170) =                   'PH_KHS_DIA'           !240!   Half saturation for pH correction fator for diatoms
      bname(171) =                   'PH_KHS_CYN'           !241!   Half saturation for pH correction fator for cyanobacteria not fixing
      bname(172) =                   'PH_KHS_FIX_CYN'       !242!   Half saturation for pH correction fator for nitrogen fixers
      bname(173) =                   'PH_KHS_OPA'           !243!   Half saturation for pH correction fator for other algae 
      bname(174) =                   'k_OX_FE_II'           !244!    Oxidation rate for iron 2+   ! Metals
      bname(175) =                   'k_RED_FE_III'         !245!    reduction rate for iron 3+
      bname(176) =                   'k_OX_MN_II'           !246!    oxidation rate for manganese 2+
      bname(177) =                   'k_RED_MN_IV'          !247!    reduction rate for manganese 4+
      bname(178) =                   'KHS_DOXY_FE_III_RED'  !248!    reversed Monod half saturation of DOXY for iron 3+ reduction
      bname(179) =                   'KHS_DOXY_MN_IV_RED'   !249!    reversed Monod half saturation of DOXY for manganese 4+ reduction       
      bname(180) =                   'K_MIN_DOC_DOXY_20'    ! New model constants introduced 27 January 2016 for the redox sequences:
      bname(181) =                   'K_MIN_DOC_NO3N_20'
      bname(182) =                   'K_MIN_DOC_MN_IV_20'
      bname(183) =                   'K_MIN_DOC_FE_III_20'
      bname(184) =                   'K_MIN_DOC_S_PLUS_6_20'
      bname(185) =                   'K_MIN_DOC_DOC_20'
      bname(186) =                   'THETA_K_MIN_DOC_DOXY'
      bname(187) =                   'THETA_K_MIN_DOC_NO3N'
      bname(188) =                   'THETA_K_MIN_DOC_MN_IV'
      bname(189) =                   'THETA_K_MIN_DOC_FE_III'
      bname(190) =                   'THETA_K_MIN_DOC_S_PLUS_6'
      bname(191) =                   'THETA_K_MIN_DOC_DOC'
      bname(192) =                   'K_HS_DOC_MIN_DOXY'
      bname(193) =                   'K_HS_DOC_MIN_NO3N'
      bname(194) =                   'K_HS_DOC_MIN_MN_IV'
      bname(195) =                   'K_HS_DOC_MIN_FE_III'
      bname(196) =                   'K_HS_DOC_MIN_S_PLUS_6'
      bname(197) =                   'K_HS_DOC_MIN_DOC'
      bname(198) =                   'K_HS_DOXY_RED_LIM'
      bname(199) =                   'K_HS_NO3N_RED_LIM'
      bname(200) =                   'K_HS_MN_IV_RED_LIM'
      bname(201) =                   'K_HS_FE_III_RED_LIM'
      bname(202) =                   'K_HS_S_PLUS_6_RED_LIM'
      bname(203) =                   'K_HS_DOXY_RED_INHB'
      bname(204) =                   'K_HS_NO3N_RED_INHB'
      bname(205) =                   'K_HS_MN_IV_RED_INHB'
      bname(206) =                   'K_HS_FE_III_RED_INHB'
      bname(207) =                   'K_HS_S_PLUS_6_RED_INHB'
      bname(208) =                   'PH_MIN_DOC_MIN_DOXY'
      bname(209) =                   'PH_MIN_DOC_MIN_NO3N'
      bname(210) =                   'PH_MIN_DOC_MIN_MN_IV'
      bname(211) =                   'PH_MIN_DOC_MIN_FE_III'
      bname(212) =                   'PH_MIN_DOC_MIN_S_PLUS_6'
      bname(213) =                   'PH_MIN_DOC_MIN_DOC'
      bname(214) =                   'PH_MAX_DOC_MIN_DOXY'
      bname(215) =                   'PH_MAX_DOC_MIN_NO3N'
      bname(216) =                   'PH_MAX_DOC_MIN_MN_IV'
      bname(217) =                   'PH_MAX_DOC_MIN_FE_III'
      bname(218) =                   'PH_MAX_DOC_MIN_S_PLUS_6'
      bname(219) =                   'PH_MAX_DOC_MIN_DOC' 
      bname(220) =                   'K_MIN_DON_DOXY_20' ! New model constants introduced 28 January 2016 for the redox sequences
      bname(221) =                   'K_MIN_DON_NO3N_20'
      bname(222) =                   'K_MIN_DON_MN_IV_20'
      bname(223) =                   'K_MIN_DON_FE_III_20'
      bname(224) =                   'K_MIN_DON_S_PLUS_6_20'
      bname(225) =                   'K_MIN_DON_DOC_20'
      bname(226) =                   'THETA_K_MIN_DON_DOXY'
      bname(227) =                   'THETA_K_MIN_DON_NO3N'
      bname(228) =                   'THETA_K_MIN_DON_MN_IV'
      bname(229) =                   'THETA_K_MIN_DON_FE_III'
      bname(230) =                   'THETA_K_MIN_DON_S_PLUS_6'
      bname(231) =                   'THETA_K_MIN_DON_DOC'
      bname(232) =                   'K_HS_DON_MIN_DOXY'
      bname(233) =                   'K_HS_DON_MIN_NO3N'
      bname(234) =                   'K_HS_DON_MIN_MN_IV'
      bname(235) =                   'K_HS_DON_MIN_FE_III'
      bname(236) =                   'K_HS_DON_MIN_S_PLUS_6'
      bname(237) =                   'K_HS_DON_MIN_DOC'
      bname(238) =                   'PH_MIN_DON_MIN_DOXY'
      bname(239) =                   'PH_MIN_DON_MIN_NO3N'
      bname(240) =                   'PH_MIN_DON_MIN_MN_IV'
      bname(241) =                   'PH_MIN_DON_MIN_FE_III'
      bname(242) =                   'PH_MIN_DON_MIN_S_PLUS_6'
      bname(243) =                   'PH_MIN_DON_MIN_DOC'
      bname(244) =                   'PH_MAX_DON_MIN_DOXY'
      bname(245) =                   'PH_MAX_DON_MIN_NO3N'
      bname(246) =                   'PH_MAX_DON_MIN_MN_IV'
      bname(247) =                   'PH_MAX_DON_MIN_FE_III'
      bname(248) =                   'PH_MAX_DON_MIN_S_PLUS_6'
      bname(249) =                   'PH_MAX_DON_MIN_DOC'
      bname(250) =                   'K_MIN_DOP_DOXY_20'
      bname(251) =                   'K_MIN_DOP_NO3N_20'
      bname(252) =                   'K_MIN_DOP_MN_IV_20'
      bname(253) =                   'K_MIN_DOP_FE_III_20'
      bname(254) =                   'K_MIN_DOP_S_PLUS_6_20'
      bname(255) =                   'K_MIN_DOP_DOC_20'
      bname(256) =                   'THETA_K_MIN_DOP_DOXY'
      bname(257) =                   'THETA_K_MIN_DOP_NO3N'
      bname(258) =                   'THETA_K_MIN_DOP_MN_IV'
      bname(259) =                   'THETA_K_MIN_DOP_FE_III'
      bname(260) =                   'THETA_K_MIN_DOP_S_PLUS_6'
      bname(261) =                   'THETA_K_MIN_DOP_DOC'
      bname(262) =                   'K_HS_DOP_MIN_DOXY'
      bname(263) =                   'K_HS_DOP_MIN_NO3N'
      bname(264) =                   'K_HS_DOP_MIN_MN_IV'
      bname(265) =                   'K_HS_DOP_MIN_FE_III'
      bname(266) =                   'K_HS_DOP_MIN_S_PLUS_6'
      bname(267) =                   'K_HS_DOP_MIN_DOC'
      bname(268) =                   'PH_MIN_DOP_MIN_DOXY'
      bname(269) =                   'PH_MIN_DOP_MIN_NO3N'
      bname(270) =                   'PH_MIN_DOP_MIN_MN_IV'
      bname(271) =                   'PH_MIN_DOP_MIN_FE_III'
      bname(272) =                   'PH_MIN_DOP_MIN_S_PLUS_6'
      bname(273) =                   'PH_MIN_DOP_MIN_DOC'
      bname(274) =                   'PH_MAX_DOP_MIN_DOXY'
      bname(275) =                   'PH_MAX_DOP_MIN_NO3N'
      bname(276) =                   'PH_MAX_DOP_MIN_MN_IV'
      bname(277) =                   'PH_MAX_DOP_MIN_FE_III'
      bname(278) =                   'PH_MAX_DOP_MIN_S_PLUS_6'
      bname(279) =                   'PH_MAX_DOP_MIN_DOC' 
      bname(280) =                   'k_OX_CH4'        ! New model constats added in 29 th of January 2016
      bname(281) =                   'THETA_k_OX_CH4'
      bname(282) =                   'k_HS_OX_CH4_DOXY'
      bname(283) =                   'k_OX_H2S'
      bname(284) =                   'THETA_k_OX_H2S'
      bname(285) =                   'k_HS_OX_H2S_DOXY'  
      bname(286) =                   'k_DISS_FE_II_20'     ! New model constants added in 9th August 2016   
      bname(287) =                   'THETA_k_DISS_FE_II'    
      bname(288) =                   'INIT_MULT_FE_II_DISS'             
      bname(289) =                   'k_DISS_FE_III_20'      
      bname(290) =                   'THETA_k_DISS_FE_III'   
      bname(291) =                   'INIT_MULT_FE_III_DISS' 
      ! End    of new model constants added in 9th August 2016


      end  !subroutine cur_param_read_ini_wc

!********************************************************************
!********************************************************************
!********************************************************************
      subroutine cur_param_read_ini_bs(bname,nconst)
!-----------------------------------------------------------------
!     Initializes BS model constants names array
!
!     Input/Output:
!      bname  -  names array
!      nconst -  number of parameters requiried
!      Note:
!           Number of parameters in this routine is hardcoded by npar
!      Called:
!           In CUR_PARAM_READ
!------------------------------------------------------------------

      integer npar !existing number of parameters in the routine!
      !change number of parameters in the next statement if it is changed!
      parameter (npar = 171)
      character*80 bname(npar)
      real par(npar)

!     execution part

      if( npar. ne. nconst) then
        print *, 'CUR_PARAM_READ_INI_BS:'
        print *, '  wrong required number of params'
        print *, '  required: ',nconst
        print *, '  existing: ',npar
        stop
      end if

!     Bottom sediment kinetics parameter names. These lines should be corrected in order
!     to introduce new or change the name of existing parameters
         bname(1 ) = 'K_OXIC_DISS_POC    '   !1: Dissolution rate constant of particulate organic carbon at 20 C (aerobic) - 1/day
         bname(2 ) = 'K_ANOXIC_DISS_POC  '   !2: Dissolution rate constant of particulate organic carbon at 20 C (anoxic) - 1/day
         bname(3 ) = 'THETA_DISS_POC     '   !3: Temperature correction for dissolution of particulate organic carbon
         bname(4 ) = 'KHS_DISS_POC       '   !4: ** Half saturation concentration of POC for dissolution
         bname(5 ) = 'K_OXIC_DISS_PON    '   !5: Dissolution rate constant of particulate organic nitrogen at 20 C (aerobic) - 1/day
         bname(6 ) = 'K_ANOXIC_DISS_PON  '   !6: Dissolution rate constant of particulate organic nitrogen at 20 C (anoxic) - 1/day
         bname(7 ) = 'THETA_DISS_PON     '   !7: Temperature correction for dissolution of particulate organic nitrogen
         bname(8 ) = 'KHS_DISS_PON       '   !8: ** Half saturation concentration of PON for dissolution
         bname(9 ) = 'K_OXIC_DISS_POP    '   !9: Dissolution rate constant of particulate organic phosphorus at 20 C (aerobic) - 1/day
         bname(10) = 'K_ANOXIC_DISS_POP  '   !10: Dissolution rate constant of particulate organic phosphorus at 20 C (anoxic) - 1/day
         bname(11) = 'THETA_DISS_POP     '   !11: Temperature correction for dissolution of particulate organic phosphorus
         bname(12) = 'KHS_DISS_POP       '   !12: ** Half saturation concentration of POP for dissolution
         bname(13) = 'K_OXIC_DISS_PSi    '   !13: Dissolution rate constant of particulate silicon at 20 C (aerobic) - 1/day
         bname(14) = 'K_ANOXIC_DISS_PSi  '   !14: Dissolution rate constant of particulate silicon at 20 C (anoxic) - 1/day
         bname(15) = 'THETA_DISS_PSi     '   !15: Temperature correction for dissolution of particulate silicon
         bname(16) = 'KHS_DISS_PSi       '   !16: ** Half saturation concentration of PSi for dissolution
         bname(17) = 'K_OXIC_MINER_DOC   '   !17: ** Mineralization rate constant of dissolved organic carbon at 20 C (aerobic) - 1/day
         bname(18) = 'K_ANOXIC_MINER_DOC '   !18: ** Mineralization rate constant of dissolved organic carbon at 20 C (anoxic) - 1/day
         bname(19) = 'THETA_MINER_DOC    '   !19: ** Temperature correction for dissolution of dissolved organic carbon
         bname(20) = 'KHS_MINER_DOC      '   !20: ** Half saturation concentration of DOC for mineralization
         bname(21) = 'K_OXIC_MINER_DON   '   !21: ** Mineralization rate constant of dissolved organic nitrogen at 20 C (aerobic) - 1/day
         bname(22) = 'K_ANOXIC_MINER_DON '   !22: ** Mineralization rate constant of dissolved organic nitrogen at 20 C (anoxic) - 1/day
         bname(23) = 'THETA_MINER_DON    '   !23: ** Temperature correction for dissolution of dissolved organic nitrogen
         bname(24) = 'KHS_MINER_DON      '   !24: ** Half saturation concentration of DON for mineralization
         bname(25) = 'K_OXIC_MINER_DOP   '   !25: ** Mineralization rate constant of dissolved organic phosphorus at 20 C (aerobic) - 1/day
         bname(26) = 'K_ANOXIC_MINER_DOP '   !26: ** Mineralization rate constant of dissolved organic phosphorus at 20 C (anoxic) - 1/day
         bname(27) = 'THETA_MINER_DOP    '   !27: ** Temperature correction for dissolution of dissolved organic phosphorus
         bname(28) = 'KHS_MINER_DOP      '   !28: ** Half saturation concentration of DOP for mineralization
         bname(29) = 'O_TO_C             '   !29: Oxygen to carbon ratio
         bname(30) = 'K_NITR             '   !30: Nitrification rate constant at 20 C - 1/day
         bname(31) = 'THETA_NITR         '   !31: Temperature correction for nitrification
         bname(32) = 'KHS_NITR_NH4N      '   !32: Half saturation constant of nitrification for NH4N - mg/L N
         bname(33) = 'KHS_NITR_DOXY      '   !33: Half saturation constant of nitrification for DOXY - mg/L O2
         bname(34) = 'K_DENITR           '   !34: Denitrification rate constant at 20 C - 1/day
         bname(35) = 'THETA_DENITR       '   !35: Temperature correction for denitrification
         bname(36) = 'KHS_DENITR_NO3N    '   !36: Half saturation constant of denitrification for NO3N - mg/L N
         bname(37) = 'KHS_DENITR_DOC     '   !37: Half saturation constant of denitrification for DOC - mg/L C
         bname(38) = 'KHS_DENITR_DOXY    '   !38: Half saturation constant of denitrification for DOXY - mg/L O
         bname(39) = 'DENITR_YIELD       '   !39: Denitrification yield
         bname(40) = 'DOXY_AT_ANOXIA     '   !40: DOXY, under which anoxia begins - mg/L O2
         bname(41) = 'SOLID_PART_COEFF_NH4'  !41: Solid part coeff for ammonium nitrogen (kg^-1)
         bname(42) = 'SOLID_PART_COEFF_PO4'  !42: Solid part coeff for phosphate phosphorus (kg^-1)

         bname(43) = 'SED_PH_MIN_DOC_MIN'        !43:   optimum lower range for pH correction factor for DOC
         bname(44) = 'SED_PH_MIN_DOC_MAX'        !44:   optimum upper range for pH correction factor for DOC
         bname(45) = 'SED_PH_MIN_DON_MIN'        !45:   optimum lower range for pH correction factor for DON
         bname(46) = 'SED_PH_MIN_DON_MAX'        !46:   optimum upper range for pH correction factor for DON
         bname(47) = 'SED_PH_MIN_DOP_MIN'        !47:   optimum lower range for pH correction factor for DOP
         bname(48) = 'SED_PH_MIN_DOP_MAX'        !48:   optimum upper range for pH correction factor for DOP
         bname(49) = 'SED_PH_NITR_NH4_MIN'       !49:   optimum lower range for pH correction factor for nitrification
         bname(50) = 'SED_PH_NITR_NH4_MAX'       !50:   optimum upper range for pH correction factor for nitrification
         bname(51) = 'SED_PH_DENITR_NO3_MIN'     !51:   optimum lower range for pH correction factor for denitrification
         bname(52) = 'SED_PH_DENITR_NO3_MAX'     !52:   optimum upper range for pH correction factor for denitrification

         bname(53)=  'SED_k_OX_FE_II'           !53:    Oxidation rate for iron 2+
         bname(54)=  'SED_k_RED_FE_III'         !54:    reduction rate for iron 3+
         bname(55)=  'SED_k_OX_MN_II'           !55:    oxidation rate for manganese 2+
         bname(56)=  'SED_k_RED_MN_IV'          !56:    reduction rate for manganese 4+
         bname(57)=  'SED_KHS_DOXY_FE_III_RED'  !57:    reversed Monod half saturation of DOXY for iron 3+ reduction
         bname(58)=  'SED_KHS_DOXY_MN_IV_RED'   !58:    reversed Monod half saturation of DOXY for manganese 4+ reduction

! New model constats added in 28 th of January 2016
         bname(59) = 'SED_K_MIN_DOC_DOXY_20'
         bname(60) = 'SED_K_MIN_DOC_NO3N_20'
         bname(61) = 'SED_K_MIN_DOC_MN_IV_20'
         bname(62) = 'SED_K_MIN_DOC_FE_III_20'
         bname(63) = 'SED_K_MIN_DOC_S_PLUS_6_20'
         bname(64) = 'SED_K_MIN_DOC_DOC_20'
         bname(65) = 'SED_THETA_K_MIN_DOC_DOXY'
         bname(66) = 'SED_THETA_K_MIN_DOC_NO3N'
         bname(67) = 'SED_THETA_K_MIN_DOC_MN_IV'
         bname(68) = 'SED_THETA_K_MIN_DOC_FE_III'
         bname(70) = 'SED_THETA_K_MIN_DOC_S_PLUS_6'
         bname(71) = 'SED_THETA_K_MIN_DOC_DOC'
         bname(72) = 'SED_K_HS_DOC_MIN_DOXY'
         bname(73) = 'SED_K_HS_DOC_MIN_NO3N'
         bname(74) = 'SED_K_HS_DOC_MIN_MN_IV'
         bname(75) = 'SED_K_HS_DOC_MIN_FE_III'
         bname(76) = 'SED_K_HS_DOC_MIN_S_PLUS_6'
         bname(77) = 'SED_K_HS_DOC_MIN_DOC'
         bname(78) = 'SED_K_HS_DOXY_RED_LIM'
         bname(79) = 'SED_K_HS_NO3N_RED_LIM'
         bname(80) = 'SED_K_HS_MN_IV_RED_LIM'
         bname(81) = 'SED_K_HS_FE_III_RED_LIM'
         bname(82) = 'SED_K_HS_S_PLUS_6_RED_LIM'
         bname(83) = 'SED_K_HS_DOXY_RED_INHB'
         bname(84) = 'SED_K_HS_NO3N_RED_INHB'
         bname(85) = 'SED_K_HS_MN_IV_RED_INHB'
         bname(86) = 'SED_K_HS_FE_III_RED_INHB'
         bname(87) = 'SED_K_HS_S_PLUS_6_RED_INHB'
         bname(88) = 'SED_PH_MIN_DOC_MIN_DOXY'
         bname(89) = 'SED_PH_MIN_DOC_MIN_NO3N'
         bname(90) = 'SED_PH_MIN_DOC_MIN_MN_IV'
         bname(91) = 'SED_PH_MIN_DOC_MIN_FE_III'
         bname(92) = 'SED_PH_MIN_DOC_MIN_S_PLUS_6'
         bname(93) = 'SED_PH_MIN_DOC_MIN_DOC'
         bname(94) = 'SED_PH_MAX_DOC_MIN_DOXY'
         bname(95) = 'SED_PH_MAX_DOC_MIN_NO3N'
         bname(96) = 'SED_PH_MAX_DOC_MIN_MN_IV'
         bname(97) = 'SED_PH_MAX_DOC_MIN_FE_III'
         bname(98) = 'SED_PH_MAX_DOC_MIN_S_PLUS_6'
         bname(99) = 'SED_PH_MAX_DOC_MIN_DOC'
         bname(100) = 'SED_K_MIN_DON_DOXY_20'
         bname(101) = 'SED_K_MIN_DON_NO3N_20'
         bname(102) = 'SED_K_MIN_DON_MN_IV_20'
         bname(103) = 'SED_K_MIN_DON_FE_III_20'
         bname(104) = 'SED_K_MIN_DON_S_PLUS_6_20'
         bname(105) = 'SED_K_MIN_DON_DOC_20'
         bname(106) = 'SED_THETA_K_MIN_DON_DOXY'
         bname(107) = 'SED_THETA_K_MIN_DON_NO3N'
         bname(108) = 'SED_THETA_K_MIN_DON_MN_IV'
         bname(109) = 'SED_THETA_K_MIN_DON_FE_III'
         bname(110) = 'SED_THETA_K_MIN_DON_S_PLUS_6'
         bname(111) = 'SED_THETA_K_MIN_DON_DOC'
         bname(112) = 'SED_K_HS_DON_MIN_DOXY'
         bname(113) = 'SED_K_HS_DON_MIN_NO3N'
         bname(114) = 'SED_K_HS_DON_MIN_MN_IV'
         bname(115) = 'SED_K_HS_DON_MIN_FE_III'
         bname(116) = 'SED_K_HS_DON_MIN_S_PLUS_6'
         bname(117) = 'SED_K_HS_DON_MIN_DOC'
         bname(118) = 'SED_PH_MIN_DON_MIN_DOXY'
         bname(119) = 'SED_PH_MIN_DON_MIN_NO3N'
         bname(120) = 'SED_PH_MIN_DON_MIN_MN_IV'
         bname(121) = 'SED_PH_MIN_DON_MIN_FE_III'
         bname(122) = 'SED_PH_MIN_DON_MIN_S_PLUS_6'
         bname(123) = 'SED_PH_MIN_DON_MIN_DOC'
         bname(124) = 'SED_PH_MAX_DON_MIN_DOXY'
         bname(125) = 'SED_PH_MAX_DON_MIN_NO3N'
         bname(126) = 'SED_PH_MAX_DON_MIN_MN_IV'
         bname(127) = 'SED_PH_MAX_DON_MIN_FE_III'
         bname(128) = 'SED_PH_MAX_DON_MIN_S_PLUS_6'
         bname(129) = 'SED_PH_MAX_DON_MIN_DOC'
         bname(130) = 'SED_K_MIN_DOP_DOXY_20'
         bname(131) = 'SED_K_MIN_DOP_NO3N_20'
         bname(132) = 'SED_K_MIN_DOP_MN_IV_20'
         bname(133) = 'SED_K_MIN_DOP_FE_III_20'
         bname(134) = 'SED_K_MIN_DOP_S_PLUS_6_20'
         bname(135) = 'SED_K_MIN_DOP_DOC_20'
         bname(136) = 'SED_THETA_K_MIN_DOP_DOXY'
         bname(137) = 'SED_THETA_K_MIN_DOP_NO3N'
         bname(138) = 'SED_THETA_K_MIN_DOP_MN_IV'
         bname(139) = 'SED_THETA_K_MIN_DOP_FE_III'
         bname(140) = 'SED_THETA_K_MIN_DOP_S_PLUS_6'
         bname(141) = 'SED_THETA_K_MIN_DOP_DOC'
         bname(142) = 'SED_K_HS_DOP_MIN_DOXY'
         bname(143) = 'SED_K_HS_DOP_MIN_NO3N'
         bname(144) = 'SED_K_HS_DOP_MIN_MN_IV'
         bname(145) = 'SED_K_HS_DOP_MIN_FE_III'
         bname(146) = 'SED_K_HS_DOP_MIN_S_PLUS_6'
         bname(147) = 'SED_K_HS_DOP_MIN_DOC'
         bname(148) = 'SED_PH_MIN_DOP_MIN_DOXY'
         bname(149) = 'SED_PH_MIN_DOP_MIN_NO3N'
         bname(150) = 'SED_PH_MIN_DOP_MIN_MN_IV'
         bname(151) = 'SED_PH_MIN_DOP_MIN_FE_III'
         bname(152) = 'SED_PH_MIN_DOP_MIN_S_PLUS_6'
         bname(153) = 'SED_PH_MIN_DOP_MIN_DOC'
         bname(154) = 'SED_PH_MAX_DOP_MIN_DOXY'
         bname(155) = 'SED_PH_MAX_DOP_MIN_NO3N'
         bname(156) = 'SED_PH_MAX_DOP_MIN_MN_IV'
         bname(157) = 'SED_PH_MAX_DOP_MIN_FE_III'
         bname(158) = 'SED_PH_MAX_DOP_MIN_S_PLUS_6'
         bname(159) = 'SED_PH_MAX_DOP_MIN_DOC'

! New model constats added in 29 th of January 2016
         bname(160) =  'SED_k_OX_CH4'
         bname(161) =  'SED_THETA_k_OX_CH4'
         bname(162) =  'SED_k_HS_OX_CH4_DOXY'
         bname(163) =  'SED_k_OX_H2S'
         bname(164) =  'SED_THETA_k_OX_H2S'
         bname(165) =  'SED_k_HS_OX_H2S_DOXY'
         
! 21 August 2016       
         bname(166) = 'SED_k_DISS_FE_II_20'
         bname(167) = 'SED_THETA_k_DISS_FE_II'
         bname(168) = 'SED_INIT_MULT_FE_II_DISS'
         bname(169) = 'SED_k_DISS_FE_III_20'
         bname(170) = 'SED_THETA_k_DISS_FE_III'
         bname(171) = 'SED_INIT_MULT_FE_III_DISS'         


      end !subroutine cur_param_read_ini_bs

!********************************************************************
!********************************************************************
!********************************************************************

      subroutine cur_param_read(what,par,nconst)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Reads parameters from file with 'biocon' or 'bioscon'
! flag in control file *.str name section using function nrdnxt
! and writes to an array (data base) using routine para_add_value
! Constants can be obtained using routine para_get_value
!
! Inputs:
!   what
!      'wc' - for water column
!      'bs' - for bottom sediments
!   nnconst - number of parameters
! Output:
!   par     - parameters array. Is initialized but not used more. Fixme
! Uses:
!   module para_aqua (it is necessary to all routines mentioned above)
! Called:
!   by AQUABCINI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use para_aqua
      implicit none

      !include 'aquabc_II.h'
      integer nconst
      character*2 what !'wc' - for water column
                       !'bs' - for bottom sediments

      integer nrdnxt,iw,ioff
      integer ifileo
      integer iunit
      double precision value, dvalue
      !double precision getpar
      character*80 file,name,text,line,bname(nconst)
      integer j,iassign
      integer ifile

C     ARRAY FOR MODEL CONSTANTS
      real par(nconst)

C     MODEL CONSTANTS DATA TYPE
      integer ctype

C     GENERAL PURPOSE COUNTER FOR ARRAY INDICES, ....
      integer i

C     INITIALIZE ARRAYS BEFORE READING

      do i=1, nconst
          par(i)   = 0.0
          bname(i) = ' '
      end do

! change number of parameters in the next statement if it is changed!
!     if(nconst .ne. 226) then
!      print *,
!    *   'CUR_PARAM_READ: Number of parameters to be initialized
!    *    is not equal to nconst'
!      stop
!     end if

!     Initialize contants names array
!

      select case(what)
       case ('wc')
        call cur_param_read_ini_wc(bname,nconst)
       case ('bs')
        call cur_param_read_ini_bs(bname,nconst)
       case default
        print *, 'CUR_PARAM_READ;'
        print *, '  wrong value of variable "what":'
        print *, what
        stop
      end select

!     Making all values zero to control reading correctnes
!
      do i=1, nconst
        par(i) = 0.0
      end do

!--------------------------------------------------------
!     Get file name
!--------------------------------------------------------
      select case(what)
       case ('wc')
        call getfnm('biocon',file)
        write(*,*) 'AQUABC_II WC parameters file:',file
       case ('bs')
        call getfnm('bioscon',file)
        write(*,*) 'AQUABC_II BS parameters file:',file
       case default
        print *, 'CUR_PARAM_READ;'
        print *, '  wrong value of variable what:'
        print *, what
        stop
      end select

      if(file .eq. ' ') then
       print *,'CUR_PARAM_READ: Bad file name'
       print *, 'what=', what, 'name=', trim(name)
       stop
      end if
!
! ----------------- reading constants file--------------------
      if( file .ne. ' ' ) then

        iunit = ifileo(0,file,'form','old')
        !call trimline(file,ifile)
!        write(6,*)'Constants initialized from file: ',file(1:ifile)
        if( iunit .le. 0 ) then
         !write(6,'(a,a)') 'filename: ',file(1:ifile)
         write(6,'(a,a)') 'filename: ',trim(adjustl(file))
         stop 'PARAM_READ: error stop: Cannot open parameters file'
        end if

!       initialization of nrdnxt routine
        call nrdini(iunit)

!--------------------------------------------------------
!         ethernal loop on lines
!--------------------------------------------------------
        do

          iw = nrdnxt(name,dvalue,text)

!   nrdnxt outputs the type of variable read :
!			-1 : error
!			 0 : end of line, nothing read
!			 1 : number variable with name
!			 2 : number variable without name
!			 3 : character variable with name
!			 4 : character variable without name

            if( iw .le. 0 ) exit

!            print *,'name=',name,' value= ', dvalue

             if( iw .eq. 1 ) then

              iassign = 0
              do j = 1,nconst
                !call triml(bname(j))
                if (name .eq. trim(adjustl(bname(j)))) then

                  par(j) = dvalue
                  call para_add_value(name,dvalue)
                  iassign = 1
                  !call para_info()
                  !call para_get_value(name,value)

                  !write(6,*) '%%%%%%%%%%',
!     *                    name,' ',value
                end if  !(name .eq. bname(j))
              end do !j
              if(iassign .eq. 0) then
               print *, 'CUR_PARAM_READ: Name of read parameter'
               print *, 'is not found on parameter list'
               print *, 'Parameter name:', trim(name)
               print *, 'Add it to the lists in routines:'
               print *, 'cur_param read_ini_wc or'
               print *, 'cur_param read_ini_bs'
               stop
              end if

             end if !(iw .eq. 1)
          end do  ! end of ethernal loop
        end if  !file .ne. ' '
!------ end of reading constants file-------

!--------------------------------------------------------
!       Writes constant value to the screen
!--------------------------------------------------------

        write(6,*) 'Number of parameters ',nconst

        do j = 1,nconst
          if (bname(j). ne. ' ' ) write(*,44) j,trim(adjustl(bname(j))),
     *        par(j)
        end do

        write(*,*)
44      format(3x,i3,2x,a40,f12.7)

        !stop

        end !param_read



c********************************************************************
c********************************************************************


      subroutine aquabcini(bsedim,par,par_sed,PHTIME,PHTAB,
     *                     TEMPTIME,TEMPTAB,
     *                     NBIOTS, BIOTSNOD, BIOTSFUN,
     *                     NBIOTS_sed, BIOTSNOD_sed, BIOTSFUN_sed,
     *                     NDGTS, DGTSNOD, DGTSFUN,
     *                     NDGTS_sed,DGTSNOD_sed,DGTSFUN_sed)


!   Initializes aquabc routines for writing results to ASCII files
!   Reads WC and BS model parameters from files
!   Formal parameters:
!         bsedim            - shows f BS is processed
!         par,par_sed       - WC and BS model constants. Not necessary. Fixme
!         PHTIME,PHTAB      - arrays for pH interpolation from data. Not used. Fixme
!         TEMPTIME,TEMPTAB  - arrays for TEMP interpolation from data. Not used. Fixme
!         NBIOTS, BIOTSNOD, BIOTSFUN             - number of nodes, nodes numbers and I/O units for WC
!         NBIOTS_sed, BIOTSNOD_sed, BIOTSFUN_sed - number of nodes, nodes numbers and I/O units for BS
!         NDGTS, DGTSNOD, DGTSFUN                - the same for WC process rates
!         NDGTS_sed,DGTSNOD_sed,DGTSFUN_sed      - the same for BS process rates

	use aquabc_II_vars
    
	implicit none

	!include 'aquabc_II.h'
	include 'aquabc_II_aout.h'

	logical bsedim

      real par(nconst)
      real par_sed(nsconst)
      real  PHTAB(1000)
      INTEGER PHTIME(1000)

      real  TEMPTAB(1000)
      INTEGER TEMPTIME(1000)

      integer INFUN,  ifileo,nb, nb_sed

      character*80 file
      character*80 file_sed

!    Reading of model parameters(constants) for WC and BS
      call cur_param_read('wc',par,nconst)
!     Reading of model parameters(constants) for BS
      if (bsedim) then
       call cur_param_read('bs',par_sed,nsconst)
      endif

!  Initialisation of state variables ASCII output
       write(6,*) '--------------------------------------------------'
     +          ,'---------------'
       write(6,*) ''

       write(6,*) 'INITIALIZING TIME SERIES AND DIAGNOSTIC ',
     +          'ASCII OUTPUTS FOR'
     +          ,'AQUABC Water Column VARIABLES'
       write(6,*) '--------------------------------------------------'
     +          ,'---------------'
       write(6,*) ''


      call getfnm('bioaow',file)
      nb = ifileo(90,file,'f','old')

      if (nb .le. 0) then
        print *, 'AQUABCINI: Can not open ctrl file for ASCII output'
        print *, '   File name: ', file
        stop
      end if

       print *,'Control file unit number for Water Column variables ',
     +       ' ASCII output',nb
       print *,'Control file name for Water Column variables ',
     +       ' ASCII output',file

       call biotser_init(nb,nstate,
     *                        NBIOTSMX,NBIOTS,
     *                        BIOTSNOD,BIOTSFUN,
     *                        NDGTSMX,NDIAGVAR,
     *                        NDGTS,DGTSNOD,DGTSFUN)

       if (bsedim) then
       write(6,*) '--------------------------------------------------'
     +          ,'---------------'
       write(6,*) ''


       write(6,*) 'INITIALIZING TIME SERIES AND DIAGNOSTIC ',
     +          'ASCII OUTPUTS FOR'
     +          ,'AQUABC Bottom Sediment VARIABLES'
       write(6,*) '--------------------------------------------------'
     +          ,'---------------'
       write(6,*) ''


         call getfnm('bioaos',file_sed)
         nb_sed = ifileo(90,file_sed,'f','old')

      if (nb_sed .le. 0) then
        print *, 'AQUABCINI: Can not open ctrl file for ASCII output'
        print *, '   File name: ', file_sed
        stop
      end if

       print *,'Control file unit number for Botom Sediments',
     * '  variables  ASCII output',nb_sed
       print *,'Control file name for Botom Sediments variables ',
     * ' ASCII output',file_sed


       call biotser_init(nb_sed,nsstate,
     *                        NBIOTSMX_sed,NBIOTS_sed,
     *                        BIOTSNOD_sed,BIOTSFUN_sed,
     *                        NDGTSMX_sed,NDIAGVAR_sed,
     *                        NDGTS_sed,DGTSNOD_sed,DGTSFUN_sed)

       endif



! Initialisation of pH reading. Not used more  !
!       call getfnm('bioph',file)
!       INFUN = ifileo(60,file,'f','u')
!       Call READ_PH(INFUN, PHTIME, PHTAB)
!       close(INFUN)
!
! C Initialisation of temperature reading(temporary)!
!       call getfnm('biotemp',file)
! 	    INFUN = ifileo(60,file,'f','u')
! 	    Call READ_TEMP(INFUN, TEMPTIME, TEMPTAB)

	end !aquabcini

!********************************************************************
!********************************************************************
!********************************************************************

! not used
          SUBROUTINE READ_PH(INFUN, PHTIME, PH)


          IMPLICIT NONE

          INTEGER PHSIZE
          PARAMETER(PHSIZE = 1000)

          INTEGER INFUN
          INTEGER PHTIME
          real PH


          DIMENSION PHTIME(PHSIZE)
          DIMENSION PH(PHSIZE)

          INTEGER I


          DO 1010, I = 1, PHSIZE
              PHTIME(I) = 0
              PH(I) = 0.0
 1010     CONTINUE


          I = 1


 1020     CONTINUE

              READ(INFUN, *, ERR = 100, END = 200)
     *             PHTIME(I), PH(I)


              IF (I .GT. 1) THEN
                  IF (PHTIME(I) .LE. PHTIME(I - 1)) THEN
                      GOTO 101
                  END IF
              END IF


              IF ((PH(I) .LT. 0.0).OR.(PH(I) .GT. 14.0)) THEN
                  WRITE(6,*) "WRONG VALUE FOR pH"
                  GOTO 102
              END IF


              IF (I > 1000) THEN
                  GOTO 200
              END IF

              I = I + 1

          GOTO 1020



  100     CONTINUE

          WRITE(6,*) "READ_PH : File read error in pH time series file"
          WRITE(6, FMT = "(A7, I5)") "LINE : ", I
          STOP


  101     CONTINUE

          WRITE(6,*) "READ_PH : Time must be greater than previous time"
          WRITE(6, FMT = "(A7, I5)") "LINE : ", I
          STOP


  102     CONTINUE

          WRITE(6,*) "READ_PH : pH must be between 0 and 14"
          WRITE(6,*) "pH you entered is :", PH(I)
          WRITE(6, FMT = "(A7, I5)") "LINE : ", I
          STOP


  200     CONTINUE

          WRITE(6,*) "Finished reading the pH time series file"


          IF (I < PHSIZE) THEN

              PH(I) = -1

          END IF


      END SUBROUTINE

!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
! not used

          SUBROUTINE READ_TEMP(INFUN, TEMPTIME, TEMP)

          IMPLICIT NONE

          INTEGER TEMPSIZE
          PARAMETER(TEMPSIZE = 1000)

          INTEGER INFUN
          INTEGER TEMPTIME
          real TEMP


          DIMENSION TEMPTIME(TEMPSIZE)
          DIMENSION TEMP(TEMPSIZE)

          INTEGER I


          DO 1010, I = 1, TEMPSIZE
              TEMPTIME(I) = 0
              TEMP(I) = 0.0
 1010     CONTINUE


          I = 1


 1020     CONTINUE

              READ(INFUN, *, ERR = 100, END = 200)
     *             TEMPTIME(I), TEMP(I)


              IF (I .GT. 1) THEN
                  IF (TEMPTIME(I) .LE. TEMPTIME(I - 1)) THEN
                      GOTO 101
                  END IF
              END IF


              IF ((TEMP(I) .LT. -2.0).OR.(TEMP(I) .GT. 65.0)) THEN
                  WRITE(6,*) "WRONG VALUE FOR TEMPERATURE"
                  GOTO 102
              END IF


              IF (I > 1000) THEN
                  GOTO 200
              END IF

              I = I + 1

          GOTO 1020



  100     CONTINUE

          WRITE(6,*) "READ_TEMP : File read error in temperature file"
          WRITE(6, FMT = "(A7, I5)") "LINE : ", I
          STOP


  101     CONTINUE

          WRITE(6,*) "READ_TEMP : ",
     *    "Time must be greater than previous time"

          WRITE(6, FMT = "(A7, I5)") "LINE : ", I
          STOP


  102     CONTINUE

          WRITE(6,*) "READ_TEMP : ",
     *    "Temparature must be between -2.0 and 65.0"

          WRITE(6,*) "Temperature you entered is :", TEMP(I)
          WRITE(6, FMT = "(A7, I5)") "LINE : ", I
          STOP


  200     CONTINUE

          WRITE(6,*) "Finished reading the temperature time file"


          IF (I < TEMPSIZE) THEN

              TEMP(I) = -1

          END IF


      END SUBROUTINE


!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************

!  not used

      real FUNCTION GET_PH(TIME, PHTIME, PH)

          IMPLICIT NONE

          INTEGER PHSIZE
          PARAMETER(PHSIZE = 1000)

          INTEGER TIME
          INTEGER PHTIME
          real PH


          DIMENSION PHTIME(PHSIZE)
          DIMENSION PH(PHSIZE)

          INTEGER I
          INTEGER K
          INTEGER P_IND
          INTEGER N_IND
          INTEGER L_IND
          LOGICAL E_IND

          real PREV_H
          real NEXT_H
          real H_PLUS
          real H_INCR

          E_IND = .TRUE.


          DO 1010, I = 1, 1000
              IF (PH(I) .EQ. -1) THEN

                  L_IND = I - 1
                  E_IND = .FALSE.
                  EXIT
              END IF
 1010     CONTINUE

          IF (E_IND) THEN
              L_IND = 1000
          END IF


          IF (TIME.GE.PHTIME(L_IND)) THEN
              GET_PH = PH(L_IND)
          ELSE


              IF (TIME.LE.PHTIME(1)) THEN
                  GET_PH = PH(1)
              ELSE
                  DO 1020, I = 1, L_IND - 1

                      IF (PHTIME(I).GT.TIME) THEN
                          EXIT
                      END IF

 1020             CONTINUE


                  P_IND = I - 1
                  N_IND = I

                  PREV_H = 10.0**(-PH(P_IND))
                  NEXT_H = 10.0**(-PH(N_IND))

                  H_INCR = (NEXT_H - PREV_H) /
     *                     (PHTIME(N_IND) - PHTIME(P_IND))

                  H_PLUS = PREV_H +
     *                     (H_INCR * (TIME - PHTIME(P_IND)))

                  GET_PH = -LOG10(H_PLUS)

              END IF

          END IF


      END FUNCTION

!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************

! not used

      real FUNCTION GET_TEMP(TIME, TEMPTIME, TEMP)

          IMPLICIT NONE

          INTEGER TEMPSIZE
          PARAMETER(TEMPSIZE = 1000)

          INTEGER TIME
          INTEGER TEMPTIME
          real TEMP


          DIMENSION TEMPTIME(TEMPSIZE)
          DIMENSION TEMP(TEMPSIZE)

          INTEGER I
          INTEGER K
          INTEGER P_IND
          INTEGER N_IND
          INTEGER L_IND
          LOGICAL E_IND

          real T_INCR

          E_IND = .TRUE.


          DO 1010, I = 1, 1000
              IF (TEMP(I) .EQ. -1) THEN

                  L_IND = I - 1
                  E_IND = .FALSE.
                  EXIT
              END IF
 1010     CONTINUE


          IF (E_IND) THEN
              L_IND = 1000
          END IF


          IF (TIME.GE.TEMPTIME(L_IND)) THEN
              GET_TEMP = TEMP(L_IND)
          ELSE


              IF (TIME.LE.TEMPTIME(1)) THEN
                  GET_TEMP = TEMP(1)
              ELSE
                  DO 1020, I = 1, L_IND - 1

                      IF (TEMPTIME(I).GT.TIME) THEN
                          EXIT
                      END IF

 1020             CONTINUE


                  P_IND = I - 1
                  N_IND = I


                  T_INCR = (TEMP(N_IND) - TEMP(P_IND)) /
     *                     (TEMPTIME(N_IND) - TEMPTIME(P_IND))

                  GET_TEMP = TEMP(P_IND) +
     *                       (T_INCR * (TIME - TEMPTIME(P_IND)))

              END IF

          END IF


      END FUNCTION


!*******************************************************************
!*******************************************************************
!*******************************************************************
! These routines should be written in order to use restart. Fixme
       subroutine read_restart_eco(iunit)
       ! dummy routine for reading restart file. fixme
       ! check number of state variables
       integer iunit
       end

       subroutine write_restart_eco(iunit)
       ! dummy routine for writing restart file. fixme
       ! check number of state variables
       integer iunit
       end

       subroutine skip_restart_eco(iunit)
       ! dummy routine for empty reading restart file. fixme
       ! check number of state variables
       integer iunit
       end
