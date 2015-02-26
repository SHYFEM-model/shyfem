! CONTENT:
!  subroutine AQUABC_II      - controls calls to pelagic and BS routines
!  subroutine cur_euler      - calculates state variables from rates using Euler
!  FUNCTION STRANGERS(VALUE) - checks for NaNs, Ifs and strange values
!  subroutine layer_light    - calculates the intensity of light after penetration of current layer
!  subroutine light_bottom   - is used by subroutine layer_light
!
!  subroutine get_nodal_area_code(k,ncode)(placed into subn35.f) - routine to get sediment type areas
!********************************************************************
!********************************************************************




            subroutine AQUABC_II(nkn,lmax_fast, nlayers,    &
                       t, dt,                               &
                       vol_fast,vol_old_fast,area_fast,  &
                       depth_fast,vel_fast,              &
                       wtemp_fast, sal_fast,             &
                       wind_fast,atemp_fast,             &
                       ITOT_fast,FDAY,                   &
                       ice_cover_fast,                   &
                       loads,                            &
                       par,nconst,                       &
                       e_fast0,e_fast1,nstate,           &
                       wc_output, noutput,               &
                       dg, NDIAGVAR,                     &
                       bsedim,                           &
                       par_sed,nsconst,                  &
                       INIT_SED_STATE_VARS_fast,  FINAL_SED_STATE_VARS_fast,  &
                       nsstate,noslay,                   &
                       output_sed_full, nsoutput,        &
                       INTER_sed_full,  NDIAGVAR_sed)


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Aquatic biochemical cycling model
! This program serves as second level interface between SHYFEM
! and WC, BS kinetic models
!
! Consists of WC pelagic and sediment module (on development)
! Calls BS module when all WC layers are processed.
! Arguments:

! Water column:
!   nkn,lmax_fast, nlayers                    - number of nodes, max number of layers for nodes
!   t,dt                                      -
!   vol_fast,vol_old_fast,area_fast           -
!   depth_fast,vel_fast                       -
!   wtemp_fast, sal_fast                      -
!   ITOT_fast,FDAY                            -
!   loads                                     -
!   par - WC constants file                   - not used more. Fixme
!   nconst - number of WC constants           -
!   e_fast0,e_fast1,nstate                    -
!   wc_output, noutput                        -
!   INTER, NDIAGVAR
!                             -
! Bottom sediments:
!   bsedim - logical (true if BS is processed -
!   par_sed - BS constants                    - not used more. fixme
!   nsconst - number of BS constants          -
!   INIT_SED_STATE_VARS_fast                  -
!   FINAL_SED_STATE_VARS_fast                 -
!   INTER_sed_full, NDIAGVAR_sed              -
!

! Notes:
!    Incomming time and step (t,dt) should be in days
!    Program returns reals in single precision and passes to kinetic routinas in
!    double precision
!    All variables with suffix 'fast' and 'full' are for speedup of calculations by
!    vectorization:
!                   full - variable is as it is used in fem interface
!                   fast - elements of arrays are arranged according
!                          dimensions: 1 - nodes, 2 - state variables, 3 - layers
!
!    Sediment types:
!       0 = fine sand;
!       1 = silty fine sand
!       2 = coarse silt
!       3 = clayey coarse silt
!       4 = silty clayey mud
!       5 = clay
!       6 = klaipeda strait
!       7 = Baltic Sea
!

!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use aquabc_II_sed_ini
      use para_aqua
      implicit none

      include 'param.h'

        integer ipv(nkndim)	!external node numbers
        common  /ipv/ipv

        logical bsedim                      !  false - no bottom sediments processing
                                            !  true - b. sed. procesing
        integer nkn                         !  total number of nodes (reactors)
        integer nlayers                     !  total number of levels in grid
        integer nnode                       !  internal node number
        integer exn_node                    !  external node number

        integer layer                       !  processed WC layer number
        integer LAYER_max                   !  maximum number of WC layers for the node(reactor)
        integer, allocatable :: node_active_num(:) !  internal numbers of active nodes for the current layer
        integer nactive_nodes                ! number of active nodes in current layer
        real state(nstate)                   !  WC state variable [mg/L] == [g/m**3]
!        real output_wc(noutput)             !  WC state variables and derived variables for one node

! Arrays for faster calculations by vectorization for all nodes(reactors)
      real e_fast0(nkndim,nstate,nlvdim)     ! initial state vector arranged                                              !
      real e_fast1(nkndim,nstate,nlvdim)     ! final state vector arranged

      double precision, allocatable :: STATE_fast(:,:)  ! state variables array for one layer for call to pelagic kinetics
      real     states_fast(nkn,nstate)       ! state variables array for one layer
      integer  lmax_fast  (nkndim)           ! max lev. numbers
      real    depth_fast  (nkndim,nlvdim)    ! depth array
      real   depth_layer  (nkn)              ! depth array for layer
      real      vel_fast  (nkndim,nlvdim)    ! current velocity
      real    wtemp_fast  (nkndim,nlvdim)    ! water temperature
      real      sal_fast  (nkndim,nlvdim)    ! water salinity
      real     wind_fast  (nkndim)           ! wind speed
      real    atemp_fast  (nkndim)           ! air temperature
      real  ice_cover_fast(nkndim)           ! ice cover - area fraction of reactor
      real    light_fast  (nkn)              ! incident light for the current layer
      real    ITOT_fast   (nkndim)           ! incident light for the first layer
      real    light_next  (nkn)              ! incident light for the next layer

      real      vol_fast  (nkndim,nlvdim)    ! reactor volume
      real  vol_old_fast  (nkndim,nlvdim)    ! reactor volume before one step
      real     area_fast  (nkndim,nlvdim)    ! reactor area

      real             wc_output (nlvdim,nkndim,noutput)    ! array for outputs (state vars + auxilaries)
      DOUBLE PRECISION, allocatable :: WC_OUTPUTS(:,:)      ! the same double precision

	  real             dg  (nlvdim,nkndim,nstate,NDIAGVAR)  ! processes values for diagnostics
	  DOUBLE PRECISION,allocatable :: DGS (:,:,:)           ! processes values for diagnostics

	  real             rates_fast(nkn,nstate,nlayers)        ! derivatives for state variables
	  DOUBLE PRECISION,allocatable :: DERIVATIVES_fast(:,:) ! derivatives for state variables

! sediments
	  integer SEDIMENT_TYPE_fast(nkn)
      real FRACTION_OF_DEPOSITION_fast(nkn,nstate)

	  integer           n_driving_functions
      parameter        (n_driving_functions = 10)
      double precision, allocatable :: DRIVING_FUNCTIONS(:,:) ! for pelagic model

      real fluxes_from_upper_fast    (nkn,nstate)
      real fluxes_to_next_fast       (nkn,nstate) ! for upper layers
!      save fluxes_to_next_fast
      real not_deposited_fluxes_fast (nkn,nstate)
      real settling_rates_fast       (nkn,nstate)    ! to the bottom from the last layer

      double precision INIT_SED_STATE_VARS_fast  (nkn,noslay,nsstate) ! initial state array for processing
      double precision FINAL_SED_STATE_VARS_fast (nkn,noslay,nsstate) ! state array after one time step
      real             output_sed_full           (noslay,nkn,nsoutput) ! derived variables as is fem interface
      real             INTER_sed_full            (noslay,nkn,nsstate,NDIAGVAR_sed)! process rates or any intermediates as is in fem interface

      double precision FLUXES_TO_SEDIMENTS_fast(nkn,nsstate)

      double precision SED_DEPTHS_fast    (nkn,noslay)
      double precision SED_POROSITIES_fast(nkn,noslay)
      double precision SED_DENSITIES_fast (nkn,noslay)
      double precision SED_DIFFUSIONS_fast(nkn,noslay,nsstate)
      double precision SED_TEMPS_fast     (nkn,noslay)
      double precision SURF_WATER_CONCS_fast(nkn,nsstate)

      double precision FLUXES_FROM_SEDIMENTS_fast(nkn,nsstate)
      real flux  !auxilary variable
      double precision PROCESSES_sed_fast(nkn,noslay,nsstate,ndiagvar_sed)
      double precision SED_OUTPUTS_fast  (nkn,noslay,nsoutput)

      double precision PART_MIXING_COEFFS_fast (nkn,noslay, nsstate)
      double precision SED_BURRIALS_fast(nkn,noslay)
      real FLUXES_TO_ALUKAS_fast(nkn,nstate)

! End of arrays for faster calculations by vectorization for all nodes(reactors)

        real vol,volold               !volume [m**3]
        real area, depth              !depth of box [m]
        real VEL                      !velocity [m/s]
        real TEMP                     !water temperature [C]
        real WIND,AIRTMP
        real sal                      !salinity [psu] == [per mille]
        real pH
        real ITOT,FDAY		          !light intensity ly/day, photoperiod
                                      !(fraction of day = 1, compatibility with light lim algorithm)
        double precision KE           ! Background extinction coefficient, from parameters


        integer nstate                      !  number of state variables
        integer noutput                     ! number of state variables + number of derived variables for the output
        integer nconst                      !  number of parameters
        integer i,j,k,kk,l
        real t 		            !current time [day]
        real dt                       !time step [day]


        real cold(nstate)             !old state variable (for diagnostic purpose)
        real loads(nstate)	          !loading for c [g/(m**3 day)]
        real rates(nstate)            !WC state variables derivatives [g/m**3/day]

        real par(nconst)      !(800) ! WC model constants

        integer NDIAGVAR
        real INTER(nstate,NDIAGVAR)



        double precision STATE_VARIABLES(nstate)
        double precision DERIVATIVES(nstate)   !derivatives  [g/m3/day]
        double precision MODEL_CONSTANTS(nconst)


        double precision PROCESSES(nstate,NDIAGVAR)
        integer           n_saved_outputs
        parameter        (n_saved_outputs = 2)
        double precision SAVED_OUTPUTS(n_saved_outputs)   ! outputs typical for the nodes saved for the next step
                                                          ! used only for background  light extinction
                                                          ! coefficient and pH (not necessary more when co2sys introduced)
        !double precision PSTIME
        double precision SURFL
        double precision IA

        double precision WCKIN_SAVE_ARRAY(nkndim, nlvdim, n_saved_outputs)! Array for saving saved outputs
        save             WCKIN_SAVE_ARRAY

        integer    iprod
        parameter (iprod=nkndim*nlvdim)
        integer WCKIN_CALLED_BEFORE(nkndim, nlvdim) ! To save indicators if node
                                                    ! and layer are already processed. Not used yet
        data    WCKIN_CALLED_BEFORE /iprod*0/
        save    WCKIN_CALLED_BEFORE
        integer       CALLED_BEFORE   ! Value of  WCKIN_CALLED_BEFORE
        integer   nflags
        parameter(nflags =2)
        INTEGER  FLAGS(nflags) !For pelagic kinetics (1) - not used; (2): 1-surface box, 0-below surface
                               ! sediments variables

        real    par_sed(nsconst) !(100)  !  BS model parmeters
        integer nsstate         !  number of botomm sediments state variables
        integer nsconst         !  number of parameters
        Integer NOSLAY          !  number of BS layers
        integer nsoutput        !  number of BS model outputs

        real state_sed(NOSLAY,nsstate)         ! BS state variable, g/m3 of pore water or solids or sediments
        real output_sed(NOSLAY,nsoutput)       ! BS model outputs (states and and auxilaries)
        integer          NUM_SED_OUTPUTS  !nsoutput


        integer   NUM_FLUX_RECEIVING_SED_LAYERS
        parameter(NUM_FLUX_RECEIVING_SED_LAYERS = 2) !Number of first layers that receives material from WC


        double precision DRIVING_FUNCTIONS_FOR_FLUXES(3) ! Fixme


        !double precision SETTLING_VELOCITIES(nstate)!g/m2/day
        !double precision DISSOLVED_FRACTIONS(nstate)
        double precision SETTLING_RATES(nstate)     !m/day

        integer          NUM_FLUXES_TO
        double precision FLUXES_TO_SEDIMENTS(nsstate) !g/m2/day

        DOUBLE PRECISION SED_TEMPS(NOSLAY) ! BS Layers temperatures, C
        double precision SALS              ! Salinity for Diffusion coeff., psu
        double precision TS                ! Temperature for Diffusion coeff., C

        INTEGER          NUM_FLUXES_FROM_SEDIMENTS
        DOUBLE PRECISION FLUXES_FROM_SEDIMENTS(nsstate) !g/m2/day


        DOUBLE PRECISION SURF_WATER_CONCS(nsstate)  ! WC concentrations to calculate difusion between WC and 1-st BS layer, g/m3

        INTEGER   NUM_SED_DRIV
        parameter(NUM_SED_DRIV = 1)
        DOUBLE PRECISION SED_DRIVING_FUNCTIONS(NOSLAY, NUM_SED_DRIV)

        DOUBLE PRECISION SED_DIFFUSIONS(NOSLAY, nsstate)

        DOUBLE PRECISION SED_MOD_1_ALUKAS_MOLDI_C   ! Function name for Diff. coeff. calculation
        !DOUBLE PRECISION SED_DEPTHS (NOSLAY)        ! Thicknes of BS layers, m
        !DOUBLE PRECISION SED_POROSITIES(NOSLAY)     ! Sediment porosities
        !DOUBLE PRECISION SED_DENSITIES(NOSLAY)      ! Sed. density , kg/m3
        !DOUBLE PRECISION SURF_MIXLEN                ! Surface mixing length for exchange with WC, m
        !DOUBLE PRECISION SED_BURRIALS(NOSLAY)       ! Burial rates (advection) for each of BS layers, m/day
        DOUBLE PRECISION FLUXES_TO_ALUKAS(nstate)    ! Fluxes from BS to ALUKAS, g/m2/day
        !DOUBLE PRECISION ADVECTIVE_VELOCITY         ! For dissolved


        DOUBLE PRECISION SED_MODEL_CONSTANTS(nsconst)

        INTEGER SEDIMENT_TYPE ! taken from grid file

        INTEGER NUM_NOT_DEPOSITED_FLUXES
        DOUBLE PRECISION FRACTION_OF_DEPOSITION(nstate)
        real             FRACTION_OF_DEPOSITION_FACTOR

        DOUBLE PRECISION NOT_DEPOSITED_FLUXES(nstate) !g/m2/day
        !DOUBLE PRECISION PART_MIXING_COEFFS (NOSLAY, nsstate)

        DOUBLE PRECISION PSTIME, TIME_STEP ! days

        integer NDIAGVAR_sed
        real INTER_sed(NOSLAY,nsstate,NDIAGVAR_sed)      				!Intermediate variables for BS (processes)
        DOUBLE PRECISION PROCESSES_sed(NOSLAY,nsstate, NDIAGVAR_sed)    !The same double precision

!       for setling fluxes management between layers
        DOUBLE PRECISION fluxes_from_upper(nstate)
        DOUBLE PRECISION fluxes_to_next(nstate)
        !save  fluxes_to_next

! For debugging purposes
         integer ix, iy, iz, error
         real x,y,z
         logical debug_stranger
         integer STRANGERS
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



error = 0               !Indicator of errors for routine termination
debug_stranger = .true. !True if debug for strangers
debug_stranger = .false.

! Flags initialisation subroutine PELAGIC_Kinetics
FLAGS(1) = 0
FLAGS(2) = 0

if(nstate.ne.24) then
 print *, 'AQUABC: Number of state variables is wrong', nstate
 stop
end if

! Initialisation of arrays
dg  (:,:,:,:)    = 0.

e_fast1(:,:,:)   = 0.

rates_fast(:,:,:)       = 0.

SED_DRIVING_FUNCTIONS(1:noslay,1:NUM_SED_DRIV) = 0.

! Asignment of some array bounds past to subroutines in order not to change their name
NUM_FLUXES_FROM_SEDIMENTS = nstate   !nstate
NUM_FLUXES_TO             = nsstate  !nsstate Number of fluxes from WC to BS
NUM_NOT_DEPOSITED_FLUXES  = nstate   !nstate
NUM_SED_OUTPUTS           = nsoutput !nsoutput

! To doble precision for call to pelagic model
PSTIME = t
TIME_STEP = dt
MODEL_CONSTANTS(1:nconst) = par(1:nconst)



! Checking for NaNs in state vars
if (debug_stranger) then
 do i = 1,nlayers
  do j = 1, nstate
   do k = 1, nkn
    if(i .le. lmax_fast(k)) then
     if(STRANGERS(e_fast0(k,j,i)) .eq. 1) then
       !external node number
       exn_node = ipv(k)
       print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
       print *, 'AQUABC before call to pelagic:'
       print *, 'Variable ',j,'Layer ',i,'Cell ',k, exn_node
       print *, 'WC state variable is strange'
       print *, 'State=',e_fast0(k,:,i)
       print *, 'Layers=',lmax_fast(k)
       print *, 'Layers depth=',depth_fast(k,1:lmax_fast(k))
       error =1
      end if
     end if
   end do
  end do
 end do
 if (error .eq. 1) stop
end if


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       Processing water column and sinks to the bottom
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!------------------------------------------------------------
!        Preparation of sinks from WC to BS
!------------------------------------------------------------
! Settling velocities and dissolved fractions assigned in module now

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!       FRACTION OF DEPOSITED MATERIAL (TEMPORARY ASIGNMENT)
!       Executed also in case of no botom sediments, because it is needed for
!       FLX_ALUKAS_TO_SED_MOD_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FRACTION_OF_DEPOSITION_FACTOR = 0.5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do nnode = 1,nkn
         call  get_nodal_area_code(nnode,SEDIMENT_TYPE) ! routine should be vectorized fixme
         exn_node = ipv(nnode)
         if(SEDIMENT_TYPE .lt. 0 .or. SEDIMENT_TYPE .gt. 7) then
          print *, 'aquabc: incorrect sediment type:', SEDIMENT_TYPE,'node:', exn_node
          stop
         end if
         SEDIMENT_TYPE_fast(nnode) = SEDIMENT_TYPE
        end do

        SEDIMENT_TYPE_fast(:)=0    !should be switched of when grid will be prepared with sediment types: fixme
      do i=1,nstate
        where     (SEDIMENT_TYPE_fast .eq. 0  .or. SEDIMENT_TYPE .eq. 6)
          !sand and Klaipeda strait
          FRACTION_OF_DEPOSITION_fast(:,i) = 0.06D+0
        elsewhere (SEDIMENT_TYPE_fast .eq. 1)
          !silty fine sand
          FRACTION_OF_DEPOSITION_fast(:,i) = 0.1D+0 !0.3D+0
        elsewhere (SEDIMENT_TYPE_fast .eq. 2 .or. SEDIMENT_TYPE_fast .eq. 3.or. &
                   SEDIMENT_TYPE_fast .eq. 4 .or. SEDIMENT_TYPE_fast .eq. 5)
          !coarse silt, clayey coarse silt, silty clayey mud, mud
          FRACTION_OF_DEPOSITION_fast(:,i) = 0.1D+0 !0.2D+0 !0.6D+0
        elsewhere (SEDIMENT_TYPE_fast .eq. 7)
          ! Baltic sea
          FRACTION_OF_DEPOSITION_fast(:,i) = 0.D+0
        end where
      end do


!    Multiplying by scale factor for deposition control
        FRACTION_OF_DEPOSITION_fast(:,:) = &
            FRACTION_OF_DEPOSITION_fast(:,:)* FRACTION_OF_DEPOSITION_FACTOR

!------------------------------------------------------------
!     End of preparation sinks from WC to BS
!------------------------------------------------------------




!llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
!       Loop on WC layers
!llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
do  layer = 1, nlayers

       if(layer .eq. 1) then
        FLAGS(2) = 1
       else
        FLAGS(2) = 0
       end if

       states_fast(1:nkn,1:nstate) = e_fast0(1:nkn,1:nstate,layer)
       depth_layer(1:nkn)          = depth_fast(1:nkn,layer)


! Light management
!       current layer
        light_fast(:) = 0.

	    if (layer .eq. 1) then
	     light_fast(1:nkn) = ITOT_fast(1:nkn)
	    else
	     where(layer <= lmax_fast(1:nkn))
          light_fast(1:nkn) = light_next(1:nkn)
         end where
	    end if

!       next layer
        if (layer .lt. nlayers) then

           call layer_light(nkn,states_fast,    &
                            nstate,             &
                            par,                &
                            nconst,             &
                            layer+1,            &
                            lmax_fast,          &
                            depth_layer,        &
                            light_fast, light_next)
        end if

! Preparation of arrays for call to pelagic model
! without nodes that does  not exist for the current layer




! number of active nodes in current layer
      nactive_nodes = count(layer <= lmax_fast(1:nkn))



      allocate(node_active_num        (nactive_nodes))
      allocate(STATE_fast         (nactive_nodes,nstate))
      allocate(WC_OUTPUTS         (nactive_nodes,noutput))
      allocate(DERIVATIVES_fast   (nactive_nodes,nstate))
      allocate(DGS                (nactive_nodes,nstate,ndiagvar))
      allocate(DRIVING_FUNCTIONS  (nactive_nodes,10))

      STATE_fast       (:,:) = 0.
      WC_OUTPUTS       (:,:) = 0.
      DRIVING_FUNCTIONS(:,:) = 0.
      DERIVATIVES_fast (:,:) = 0.
      DGS (:,:,:)            = 0.

      j=1
      do k=1,nkn

        if(layer <= lmax_fast(k)) then
        node_active_num(j)  = k
        STATE_fast(j,1:nstate)  = e_fast0(k,1:nstate,layer)
        WC_OUTPUTS(j,1:noutput) = wc_output(layer,k,1:noutput)

! Driving functions management
        call para_get_value('KE', KE)
        DRIVING_FUNCTIONS(j, 1) = wtemp_fast(k,layer)  ! w. temp.
        DRIVING_FUNCTIONS(j, 2) = sal_fast  (k,layer)  ! salinity
        DRIVING_FUNCTIONS(j, 3) = light_fast(k)        ! Instantaneous light W/m2, conversion to langlays is done in pelagic kinetics
        DRIVING_FUNCTIONS(j, 4) = FDAY                 ! fraction of day with light. 1 - for instantenous light, just for compatibility with older code
        DRIVING_FUNCTIONS(j, 5) = atemp_fast(k)        ! air temp.
        DRIVING_FUNCTIONS(j, 6) = wind_fast (k)        ! wind spead
        DRIVING_FUNCTIONS(j, 7) = 0.0D0                ! elevation
        DRIVING_FUNCTIONS(j, 8) = depth_fast(k,layer)  ! water depth  or layer thickness in case of 3d
        DRIVING_FUNCTIONS(j, 9) = KE   !0.8D0          ! background light extinction coefficient. Later should be calculated. fixme
                                                       ! should be compatible with value in subroutine alight_alukas
        DRIVING_FUNCTIONS(j,10) = ice_cover_fast(k)    ! ice cover - area fraction of reactor                    ! was reserved for pH. Not necessary when co2sys was introduced



!       end of driving functions section

        j=j+1

       end if ! layer <=
      end do ! k (nodes)

       !print *,'*****************************************************'
       !print *,'aquabc: TEMP(175)', wtemp_fast(175,layer)
       !print *,'layer', layer

!----------------------------------------------------------------------------------
! Nothing is stored in SAVED_OUTPUTS in this version
!         DO I = 1, n_saved_outputs
!           SAVED_OUTPUTS(I) = WCKIN_SAVE_ARRAY(nnode, LAYER, I)

!         END DO

!         Driving functions are not used in FLX_ALUKAS_II_TO_SED_MOD_1 yet
!         DRIVING_FUNCTIONS_FOR_FLUXES(1) = 0.D0 !SAVED_OUTPUTS(4)
!         DRIVING_FUNCTIONS_FOR_FLUXES(2) = 0.D0 !SAVED_OUTPUTS(5)
!         DRIVING_FUNCTIONS_FOR_FLUXES(3) = 0.D0 !SAVED_OUTPUTS(6)
!
!         SED_DRIVING_FUNCTIONS(:,:) = 0.D+1  ! Driving functions are not used in sediment model yet
!
!          !CALLED_BEFORE = WCKIN_CALLED_BEFORE(nnode, LAYER)
!---------------------------------------------------------------------------------

! print *, 'Before'
! print *, 'nactive_nodes=', nactive_nodes
! print *,'node_active_num=', node_active_num(1:10)
! print *,'DRIVING_FUNCTIONS=',DRIVING_FUNCTIONS(1:10,1)
! print *,'STATE=',STATE_fast(1:10,1)



        CALL PELAGIC_KINETICS &
                (node_active_num, nactive_nodes,&
                 STATE_fast  , DERIVATIVES_fast, nstate,  &
                 MODEL_CONSTANTS  , nconst,               &
                 DRIVING_FUNCTIONS, n_driving_functions,  &
                 FLAGS            , nflags,               &
                 DGS              , NDIAGVAR,             &
                 SAVED_OUTPUTS    , n_saved_outputs,      &
                 WC_OUTPUTS       , noutput,              &
                 PSTIME, TIME_STEP, CALLED_BEFORE)



          do k=1,nactive_nodes
            rates_fast(node_active_num(k),1:nstate,layer) = DERIVATIVES_fast(k,1:nstate)
            wc_output(layer,node_active_num(k),1:noutput) = WC_OUTPUTS(k,1:noutput)
            dg(layer,node_active_num(k),1:nstate,:)       = DGS(k,1:nstate,:)
          end do

!-----
!         DO I = 1, n_saved_outputs
!              WCKIN_SAVE_ARRAY(nnode, LAYER, I) = SAVED_OUTPUTS(I)
!         END DO
!
!
!         IF (WCKIN_CALLED_BEFORE(nnode, LAYER).EQ.0) THEN
!              WCKIN_CALLED_BEFORE(nnode, LAYER) = 1
!         END IF
!-----
          if (debug_stranger) then
           do j = 1, nstate
            do k = 1, nkn
             if(layer <= lmax_fast(k)) then
              if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
!               external node number
                exn_node = ipv(k)
                print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                print *, 'AQUABC after call to pelagic kinetics:'
                print *, 'Variable ',j,'Cell ',k, exn_node
                print *, 'WC derivative is NaN'
                print *, 'Derivative=',rates_fast(k,j,layer)
                print *, 'Layers=',lmax_fast(k)
                print *, 'Layers depth=',depth_fast(k,1:lmax_fast(k))
                error =1
               end if
              end if
            end do
           end do
          if (error .eq. 1) stop
         end if

           !print *,'18- 8',' Layer ', layer,' rate:',rates_fast(18, 8,layer)
           !print *,'18-20',' Layer ', layer,' rate:',rates_fast(18,20,layer)
           !print *,'18-24',' Layer ', layer,' rate:',rates_fast(18,24,layer)



!  SETTLING FLUXES MANAGAMENT ALSO FOR 3D CASE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        fluxes_from_upper_fast(1:nkn,1:nstate) = 0.

        if (layer > 1) then
         do k=1,nkn
          if(layer <= lmax_fast(k)) then
            fluxes_from_upper_fast(k,1:nstate) = fluxes_to_next_fast(k,1:nstate)
          end if
         end do
        end if


! Particulate material fluxes to and from the  layers
! settling rates comes in g/m2/day

     fluxes_to_next_fast(1:nkn,1:nstate) = 0.

     do k=1,nkn
        nnode = k
        !---------------------------------------------------------
        if(layer .eq. lmax_fast(k)) then   !processing nodes which current layer is nearbottom layer


         STATE_VARIABLES(1:nstate)        = e_fast0(k,1:nstate,layer)
         FRACTION_OF_DEPOSITION(1:nstate) = FRACTION_OF_DEPOSITION_fast(k,1:nstate)

         CALL FLX_ALUKAS_II_TO_SED_MOD_1 &
                (STATE_VARIABLES               , nstate,              &
                 MODEL_CONSTANTS               , nconst,              &
                 DRIVING_FUNCTIONS_FOR_FLUXES  , 3,                   &
                 SETTLING_VELOCITIES, DISSOLVED_FRACTIONS,            & !inputs
                 1.0D+0 , nnode, PSTIME,                              & !bottom factor
                 SETTLING_RATES, FLUXES_TO_SEDIMENTS, NUM_FLUXES_TO,  &
                 SEDIMENT_TYPE , FRACTION_OF_DEPOSITION,              &
                 NOT_DEPOSITED_FLUXES, NUM_NOT_DEPOSITED_FLUXES)

         not_deposited_fluxes_fast(k, 1:nstate)  = NOT_DEPOSITED_FLUXES(1:nstate)
         settling_rates_fast      (k, 1:nstate)  = SETTLING_RATES      (1:nstate)
         FLUXES_TO_SEDIMENTS_fast (k, 1:nsstate) = FLUXES_TO_SEDIMENTS (1:nsstate)


         rates_fast(k,1:nstate,layer) = rates_fast(k,1:nstate,layer)                                  &
                                        - settling_rates_fast(k,1:nstate)/depth_fast(k,layer)         &
                                        + not_deposited_fluxes_fast(k, 1:nstate)/depth_fast(k,layer)  &
                                        + fluxes_from_upper_fast(k,1:nstate)/depth_fast(k,layer)


           !if(k .eq. 625) then
           !  print *, 'layer ', layer
           !  print *, 'settling_rates ',settling_rates_fast(k,21)
           !  print *, 'depth ',depth_fast(k,layer)
           !  print *, 'not_deposited_fluxes', not_deposited_fluxes_fast(k, 21)
           !  print *, 'fluxes_from_upper ',fluxes_from_upper_fast(k,21)
           ! end if

         dg(layer,k,1:nstate,24) = settling_rates_fast      (k,1:nstate)
         dg(layer,k,1:nstate,25) = not_deposited_fluxes_fast(k,1:nstate)

! Checking for strange values in WC derivatives
          if (debug_stranger) then
           do j = 1, nstate
              if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
!               external node number
                exn_node = ipv(k)
                print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                print *, 'AQUABC after updating derivatives (bottom layers):'
                print *, 'Variable ',j,' Cell ',k, exn_node,' Layer ',layer
                print *, 'Derivative is NaN'
                print *, 'Derivative ='         , rates_fast(k,j,layer)
                print *, 'Settling_rates:'      , settling_rates_fast(k,1:nstate)
                print *, 'Not_deposited_fluxes:', not_deposited_fluxes_fast(1:nkn, i)
                print *, 'Fluxes_from_upper: '  , fluxes_from_upper_fast(k,1:nstate)
                print *, 'Depth:   '            , depth_fast(k,layer)
                error =1
               end if
           end do
           if (error .eq. 1) stop
          end if


        end if ! end of processing nodes which current layer is nearbottom layer

        !------------------------------------------------------------
        if(layer .lt. lmax_fast(k)) then !processing nodes which current layer is not nearbottom layer

         STATE_VARIABLES(1:nstate)  = e_fast0(k,1:nstate,layer)
         
!        Settling velocities and dissolved fractions initialised   
!        using module aquabc_II_sed_ini        

         call SETTLING_TO_THE_NEXT_LAYER &
               (STATE_VARIABLES    , nstate   ,           &
                SETTLING_VELOCITIES, DISSOLVED_FRACTIONS, &
                nnode, LAYER, PSTIME,                     &
                fluxes_to_next)

          fluxes_to_next_fast(k,1:nstate) = fluxes_to_next(1:nstate)

          rates_fast(k,1:nstate,layer)    = rates_fast(k,1:nstate,layer)                      &
                                     - fluxes_to_next_fast   (k,1:nstate)/depth_fast(k,layer) &
                                     + fluxes_from_upper_fast(k,1:nstate)/depth_fast(k,layer)

           !if(k .eq. 625) then
           !  print *, 'layer ', layer
           !  print *, 'fluxes_to_next ',fluxes_to_next_fast(k,21)
           !  print *, 'depth ',depth_fast(k,layer)
           !  print *, 'fluxes_from_upper ',fluxes_from_upper_fast(k,21)
           ! end if

          dg(layer,k,1:nstate,24) = fluxes_to_next(1:nstate)

! Checking for strange values in WC derivatives
         if (debug_stranger) then
          do j = 1, nstate
              if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
!               external node number
                exn_node = ipv(k)
                print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                print *, 'AQUABC after updating derivatives (upper layers):'
                print *, 'Variable ',j,' Cell ',k, exn_node,' Layer ',layer
                print *, 'Derivative is NaN'
                print *, 'Derivative=',rates_fast(k,j,layer)
                print *, 'Depth:   '  , depth_fast(k,layer)
                print *, 'Max number of layers:', lmax_fast(k)
                error =1
               end if
          end do
          if (error .eq. 1) stop
         end if
        !-----------------------------------------------------------

        end if ! end of processing nodes which current layer is not nearbottom

      end do !k

! End of processing particulate material fluxes to and from the  layers


      if(allocated(node_active_num ))  deallocate(node_active_num)
      if(allocated(STATE_fast ))       deallocate(STATE_fast)
      if(allocated(WC_OUTPUTS ))       deallocate(WC_OUTPUTS)
      if(allocated(DRIVING_FUNCTIONS ))deallocate(DRIVING_FUNCTIONS)
      if(allocated(DERIVATIVES_fast )) deallocate(DERIVATIVES_fast)
      if(allocated(DGS))               deallocate(DGS)



end do !layers


!llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
!       End of loop on WC layers
!llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll

! Checking for NaN, Inf and strange values in WC derivatives
         if (debug_stranger) then
          do layer = 1,nlayers
           do j = 1, nstate
            do k = 1, nkn
             if(layer <= lmax_fast(k)) then
              if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
!               external node number
                exn_node = ipv(k)
                print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                print *, 'AQUABC after updating derivatives by settling etc:'
                print *, 'Variable ',j,'Layer ',layer,'Cell ', k, exn_node
                print *, 'WC Derivative is NaN'
                print *, 'Derivative=',rates_fast(k,j,layer)
                error =1
               end if
              end if
            end do
           end do
          end do
          if (error .eq. 1) stop
         end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     BOTTOM SEDIMENTS PROCESSING
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! In this section it is assummed that variable 'layer' is a nearbottom WC layer but has a value 1 so
! only 2d calculations are performed. fixme

  if(bsedim)  then

! Checking for strange values in BS initial state
       if (debug_stranger) then
        do j = 1, nsstate
         do i = 1, noslay
          do k=1, nkn
            flux = INIT_SED_STATE_VARS_fast(k,i,j)
           if(STRANGERS(flux) .eq. 1) then
            print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
            print *, 'AQUABC: Layer ', i, 'Variable ',j,'Cell ',k
            print *, 'BS initial state is NaN'
            print *, 'INITIAL(k,i,j)=',INIT_SED_STATE_VARS_fast(k,i,j)
            error = 1
            end if
          end do
         end do
        end do  !
        if(error.eq. 1) stop
       end if



! ************************************************************************
! Initialisation of sediment properties is done by routine initialisation
! routine from module aquabc_II_sed_ini
! Later should be read from files. fixme

!          Thickness of sediment layers in meters 
           SED_DEPTHS_fast (:,1) = SED_DEPTHS(1) !0.005D+0
           SED_DEPTHS_fast (:,2) = SED_DEPTHS(2) !0.005D+0
           SED_DEPTHS_fast (:,3) = SED_DEPTHS(3) !0.005D+0
           SED_DEPTHS_fast (:,4) = SED_DEPTHS(4) !0.005D+0
           SED_DEPTHS_fast (:,5) = SED_DEPTHS(5) !0.03D+0
           SED_DEPTHS_fast (:,6) = SED_DEPTHS(6) !0.05D+0


!          Porosities of BS layers 
           SED_POROSITIES_fast(:,1) = SED_POROSITIES(1) !0.40D+0
           SED_POROSITIES_fast(:,2) = SED_POROSITIES(2) !0.40D+0
           SED_POROSITIES_fast(:,3) = SED_POROSITIES(3) !0.40D+0
           SED_POROSITIES_fast(:,4) = SED_POROSITIES(4) !0.40D+0
           SED_POROSITIES_fast(:,5) = SED_POROSITIES(5) !0.30D+0
           SED_POROSITIES_fast(:,6) = SED_POROSITIES(6) !0.25D+0

                  
!          Densities for BS layers
           SED_DENSITIES_fast(:,1) = SED_DENSITIES(1) ! 1.75D+0
           SED_DENSITIES_fast(:,2) = SED_DENSITIES(2) ! 1.75D+0
           SED_DENSITIES_fast(:,3) = SED_DENSITIES(3) ! 1.75D+0
           SED_DENSITIES_fast(:,4) = SED_DENSITIES(4) ! 1.75D+0
           SED_DENSITIES_fast(:,5) = SED_DENSITIES(5) ! 1.75D+0
           SED_DENSITIES_fast(:,6) = SED_DENSITIES(6) ! 1.75D+0

!          Burial rate of BS state variables (around 1cm/year) !m/day  
           SED_BURRIALS_fast(:,1) = SED_BURRIALS(1) !2.7397D-5 !2.7397D-5    !m/day
           SED_BURRIALS_fast(:,2) = SED_BURRIALS(2) !2.7397D-5 !2.7397D-5
           SED_BURRIALS_fast(:,3) = SED_BURRIALS(3) !2.7397D-5 !2.7397D-5
           SED_BURRIALS_fast(:,4) = SED_BURRIALS(4) !2.7397D-5
           SED_BURRIALS_fast(:,5) = SED_BURRIALS(5) !2.7397D-5
           SED_BURRIALS_fast(:,6) = SED_BURRIALS(6) !2.7397D-5
           
!          Particle mixing diffusion coeff.  
!          Made zero for solutes in sediment routine!
           PART_MIXING_COEFFS_fast (:,1, :) = PART_MIXING_COEFFS(1) ! 2.64D-5 ! 0.d+0    
           PART_MIXING_COEFFS_fast (:,2, :) = PART_MIXING_COEFFS(2) ! 2.64D-5 ! 0.d+0
           PART_MIXING_COEFFS_fast (:,3, :) = PART_MIXING_COEFFS(3) ! 2.64D-5 ! 0.d+0
           PART_MIXING_COEFFS_fast (:,4, :) = PART_MIXING_COEFFS(4) ! 2.64D-5
           PART_MIXING_COEFFS_fast (:,5, :) = PART_MIXING_COEFFS(5) ! 2.64D-5
           PART_MIXING_COEFFS_fast (:,6, :) = PART_MIXING_COEFFS(6) ! 2.64D-5

!          Initilized in module aquabc_II_sed_ini now
           !ADVECTIVE_VELOCITY = 0. !1.0D-5 !m/s   
           
!          Surface mixing length for exchange with WC
!          Initilized in module aquabc_II_sed_ini now
           !SURF_MIXLEN = 0.10  !m 0.1
! ********************************************************************


! Checking for NaN, Inf and strange values in WC derivatives
         if (debug_stranger) then
          do layer = 1,nlayers
           do j = 1, nstate
            do k = 1, nkn
             if(layer <= lmax_fast(k)) then
              if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
!               external node number
                exn_node = ipv(k)
                print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                print *, 'AQUABC before diffusion coeff.:'
                print *, 'Variable ',j,'Layer ',layer,'Cell ', k, exn_node
                print *, 'Derivative is NaN'
                print *, 'Derivative=',rates_fast(k,j,layer)
                error =1
               end if
              end if
            end do
           end do
          end do
          if (error .eq. 1) stop
         end if

!     Calculation of diffusion coefficients
      do k =1,nkn
         do i=1,noslay
              TS                  = wtemp_fast(k,lmax_fast(k))
              SED_TEMPS_fast(k,i) = wtemp_fast(k,lmax_fast(k)) ! sed. temperatures for sediments subroutine.
                                                               ! Assumed equal to water column temp.
              SALS                = sal_fast  (k,lmax_fast(k)) ! sed. salinity. Assumed equal to water column salinity (last layer)
              do j=1,nsstate
               SED_DIFFUSIONS_fast(k,i,j) = &
                         SED_MOD_1_ALUKAS_MOLDI_C(j, TS, SALS, 0.0D+0)*86400.
              end do
         end do
      end do



!  WC concentrations corresponding to BS variables for molecular diffusion
!  Change these statments if water column vars number is changed

           if(nsstate.ne.15) then
            print *, 'AQUABC: Number of sediment state variables is wrong', nsstate
            stop
           end if

          do k=1,nkn
           SURF_WATER_CONCS_fast(k, 1)  =  e_fast0(k, 1,lmax_fast(k)) !state(1)
           SURF_WATER_CONCS_fast(k, 2)  =  e_fast0(k, 2,lmax_fast(k)) !state(2)
           SURF_WATER_CONCS_fast(k, 3)  =  e_fast0(k,13,lmax_fast(k)) !state(13)
           SURF_WATER_CONCS_fast(k, 4)  =  0.0D+0
           SURF_WATER_CONCS_fast(k, 5)  =  e_fast0(k, 3,lmax_fast(k)) !state(3)
           SURF_WATER_CONCS_fast(k, 6)  =  e_fast0(k,14,lmax_fast(k)) !state(14)
           SURF_WATER_CONCS_fast(k, 7)  =  0.0D+0
           SURF_WATER_CONCS_fast(k, 8)  =  e_fast0(k, 4,lmax_fast(k)) !state(4)
           SURF_WATER_CONCS_fast(k, 9)  =  e_fast0(k,12,lmax_fast(k)) !state(12)
           SURF_WATER_CONCS_fast(k,10)  =  0.0D+0
           SURF_WATER_CONCS_fast(k,11)  =  e_fast0(k,20,lmax_fast(k)) !state(20)
           SURF_WATER_CONCS_fast(k,12)  =  0.0D+0
           SURF_WATER_CONCS_fast(k,13)  =  e_fast0(k,23,lmax_fast(k)) !state(23)
           SURF_WATER_CONCS_fast(k,14)  =  e_fast0(k,24,lmax_fast(k)) !state(24)
           SURF_WATER_CONCS_fast(k,15)  =  sal_fast(k,lmax_fast(k))   !Salinity
          end do

    ! Checking for NaN, Inf and strange values in WC derivatives
         if (debug_stranger) then
          do layer = 1,nlayers
           do j = 1, nstate
            do k = 1, nkn
             if(layer <= lmax_fast(k)) then
              if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
!               external node number
                exn_node = ipv(k)
                print *, 'AQUABC after surf. conc:'
                print *, 'Variable ',j,'Layer ',layer,'Cell ', k, exn_node
                print *, 'WC Derivative is NaN'
                print *, 'Derivative=',rates_fast(k,j,layer)
                error =1
               end if
              end if
            end do
           end do
          end do
          if (error .eq. 1) stop
         end if

!       BS constants
          do i=1,nsconst
            SED_MODEL_CONSTANTS(i) = par_sed(i)
          enddo



!       BS state vars initial values
          !INIT_SED_STATE_VARS_fast(:,:,:) = 0.!state_sed_fast(:,:,:)

! Checking for NaN, Inf and strange values in BS initial state vars
         if (debug_stranger) then
          do i = 1, nsstate
           do j = 1,noslay
            do k = 1, nkn
              flux = INIT_SED_STATE_VARS_fast(k,j,i)
              if(STRANGERS(flux) .eq. 1) then
!               external node number
                exn_node = ipv(k)
                print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                print *, 'AQUABC before BS kinetics:'
                print *, 'Variable ',i,'Layer ',j,'Cell ', k, exn_node
                print *, 'BS state variable is NaN'
                print *, 'Value=',INIT_SED_STATE_VARS_fast(k,j,i)
                error =1
               end if
            end do
           end do
          end do
          if (error .eq. 1) stop
         end if

! Checking for NaN, Inf and strange values in WC derivatives
         if (debug_stranger) then
          do layer = 1,nlayers
           do j = 1, nstate
            do k = 1, nkn
             if(layer <= lmax_fast(k)) then
              if(STRANGERS(rates_fast(k,j,layer)) .eq. 1) then
!               external node number
                exn_node = ipv(k)
                print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
                print *, 'AQUABC befor call to BS:'
                print *, 'Variable ',j,'Layer ',layer,'Cell ', k, exn_node
                print *, 'WC Derivative is NaN'
                print *, 'Derivative=',rates_fast(k,j,layer)
                error =1
               end if
              end if
            end do
           end do
          end do
          if (error .eq. 1) stop
         end if


!       Sediment drivig functions not used yet
!       Call to the main sediment routine

       call SEDIMENT_MODEL_1 &
                 (nkn,INIT_SED_STATE_VARS_fast  , SED_DEPTHS_fast , SED_POROSITIES_fast,  &
                  SED_DENSITIES_fast        , PART_MIXING_COEFFS              ,  &
                  SED_DIFFUSIONS_fast       , SURF_MIXLEN, SED_BURRIALS_fast  ,  &
                  SURF_WATER_CONCS_fast     , SED_TEMPS_fast                  ,  &
                  nsstate                   , NOSLAY                          ,  &
                  SED_MODEL_CONSTANTS       , nsconst                         ,  &
                  SED_DRIVING_FUNCTIONS     , NUM_SED_DRIV                    ,  &
                  FLUXES_TO_SEDIMENTS_fast  , NUM_FLUXES_TO                   ,  &
                  NUM_FLUX_RECEIVING_SED_LAYERS  , ADVECTIVE_VELOCITY         ,  &
                  PSTIME, TIME_STEP         ,                                    &
                  FINAL_SED_STATE_VARS_fast ,                               &
                  FLUXES_FROM_SEDIMENTS_fast, NUM_FLUXES_FROM_SEDIMENTS  ,  &
                  PROCESSES_sed_fast        , NDIAGVAR_sed               ,  &
                  SED_OUTPUTS_fast          , NUM_SED_OUTPUTS)



   ! reshaping arrays
   do j=1,nsstate
    do i=1,noslay
     do k=1,nkn
     do l=1,ndiagvar_sed
      INTER_sed_full (i,k,j,l)  = PROCESSES_sed_fast(k,i,j,l)
     end do
     end do
    end do
   end do

   do j=1,nsoutput
    do i=1,noslay
     do k=1,nkn
     output_sed_full(i,k,j)  = SED_OUTPUTS_fast(k,i,j)
     end do
    end do
   end do

! Checking for strange values in BS final state
  if (debug_stranger) then
   do j = 1, nsstate
    do k=1,nkn
     do i = 1, noslay
      flux = FINAL_SED_STATE_VARS_fast(k,i,j)
      if(STRANGERS(flux) .eq. 1) then
       print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
       print *, 'AQUABC after call to BS:'
       print *, 'Layer ', i, 'Variable ',j,'Cell ',k
       print *, 'Final state is NaN'
       print *, 'FINAL(k,i,j)=',FINAL_SED_STATE_VARS_fast(k,i,j)
       error =1
      end if
     end do
    end do
   end do
   if(error .eq. 1) stop
  end if

! Checking for strange values in fluxes from BS
  if (debug_stranger) then
   do i = 1, nsstate
    do k=1,nkn
     flux = FLUXES_FROM_SEDIMENTS_fast(k,i)
     if(STRANGERS(flux) .eq. 1) then
      print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
      print *, 'AQUABC: Layer ', i, 'Variable ',i,'Cell ',k
      print *, 'Flux from BS to WC by BS vars is NaN'
      print *, 'Flux(k,i)=',FLUXES_FROM_SEDIMENTS_fast(k,i)
      error =1
     end if
    end do
   end do
   if(error .eq. 1) stop
  end if

!  Fluxes from BS to WC
   do k=1,nkn
    FLUXES_FROM_SEDIMENTS(:) = FLUXES_FROM_SEDIMENTS_fast(k,:)
    CALL FLX_SED_MOD_1_TO_ALUKAS_II                            &
           (FLUXES_FROM_SEDIMENTS, NUM_FLUXES_FROM_SEDIMENTS, &
            FLUXES_TO_ALUKAS     , nstate)
    FLUXES_TO_ALUKAS_fast(k,:) = FLUXES_TO_ALUKAS(:)
   end do

! Checking for strange values in fluxes to WC from BS
  if (debug_stranger) then
   do i = 1, nstate
    do k=1,nkn
     flux = FLUXES_TO_ALUKAS_fast(k,i)
     if(STRANGERS(flux) .eq. 1) then
       print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
       print *, 'AQUABC:'
       print *, 'Layer ', i, 'Variable ',i,'Cell ',k
       print *, 'Flux from BS to WC by WC vars is NaN'
       print *, 'Flux(k,i)=',FLUXES_TO_ALUKAS_fast(k,i)
       error =1
     end if
    end do
   end do
   if(error .eq. 1) stop
  end if


   ! Sending fluxes to intermediate variables
   do k=1,nkn
    do i = 1, noslay
     do j = 1, nsstate

      if (i .eq. 1) then
       INTER_sed_full(i,k,j, 8) =  FLUXES_TO_SEDIMENTS_fast(k,j)
       INTER_sed_full(i,k,j, 9) =  FLUXES_FROM_SEDIMENTS_fast(k,j)
      else
       INTER_sed_full(i,k,j, 8) = 0.D+0
       INTER_sed_full(i,k,j, 9) = 0.D+0
      end if

     end do
    end do
   end do

! Updating WC bottom layer derivatives

! Checking for strange values in WC derivatives
    if (debug_stranger) then
     do i = 1, nstate
       do k=1,nkn
!       print  *, 'node=', k, 'Variable ',i,'Layer' ,lmax_fast(k),'rate ',rates_fast(k,i,lmax_fast(k))
        if(STRANGERS(rates_fast(k,i,lmax_fast(k))) .eq. 1) then
              print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
              print *, 'AQUABC before updating by fluxes from from BS: '
              print *,  'Variable ',i,'Cell ',k, 'Layer ',lmax_fast(k)
              print *, 'WC bottom layer derivative is NaN'
              print *, 'Derivative',rates_fast(k,i,lmax_fast(k))
              print *, 'Flux from BS:', FLUXES_TO_ALUKAS_fast(k,i)
              print *, 'Layer depth: ', depth_fast(k,lmax_fast(k))
              error =1
        end if
       end do
      end do
      if(error .eq. 1) stop
     end if

! Updating WC dervatives by fluxes from BS
      do i = 1,nstate
       do k=1,nkn

         rates_fast(k,i,lmax_fast(k))  = rates_fast(k,i,lmax_fast(k)) + FLUXES_TO_ALUKAS_fast(k,i)/depth_fast(k,lmax_fast(k))
!         if(STRANGERS(rates_fast(k,i,lmax_fast(k))) .eq. 1) print *,'After ','k=',k,'i=',i, 'rates=',rates_fast(k,i,lmax_fast(k))
         dg(lmax_fast(k),k,i,26) = FLUXES_TO_ALUKAS_fast(k,i)

       end do
      end do

! Checking for strange values in WC bottom layer derivatives
     if (debug_stranger) then
      do i = 1, nstate
       do k=1,nkn
        if(STRANGERS(rates_fast(k,i,lmax_fast(k))) .eq. 1) then
              print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
              print *, 'AQUABC after call to BS: '
              print *,  'Variable ',i,'Cell ',k, 'Layer ',lmax_fast(k)
              print *, 'WC bottom layer derivative is NaN'
              print *, 'Derivative',rates_fast(k,i,lmax_fast(k))
              print *, 'Flux from BS:', FLUXES_TO_ALUKAS_fast(k,i)
              print *, 'Layer depth: ', depth_fast(k,lmax_fast(k))
              print *, 'Fluxes/Depth', FLUXES_TO_ALUKAS_fast(k,i)/depth_fast(k,lmax_fast(k))
              error =1
        end if
       end do
      end do
      if(error .eq. 1) stop
     end if

   end if  !bsedim

!ccccccccccccccccccccccccccc
!     End of BS processing c
!ccccccccccccccccccccccccccc




!      Updating derivatives by loading. Not tested: fixme
!cccccccccccccccccccccccccccccccccccccc
!       call load0d(rates,loads,vol) !Atention, check cds updates(*Vol!) fixme for alukas

!      Updating state variables by new derivatives
!ccccccccccccccccccccccccccccccccccccccccccccccccc

!        INTER(23,2) = NOT_DEPOSITED_FLUXES(23)
!        INTER(23,3) = rates(23)
!        INTER(23,4) = (FLUXES_TO_ALUKAS(23)/depth)
!        INTER(23,5) = (NOT_DEPOSITED_FLUXES(i) / depth)

           !if(exn_node.eq.2040)then
           !   print *, 'Before Euler:','r=',rates(23),'s=',state(23),'t=',t
           !end if



! Integration

       do layer = 1,nlayers
        do k=1,nkn ! fixme for vectorization
        if(layer <= lmax_fast(k)) then
         state(1:nstate) = e_fast0(k,1:nstate,layer)
         rates(1:nstate) = rates_fast(k,1:nstate,layer)
         vol             = vol_fast    (k,layer)
         volold          = vol_old_fast(k,layer)

         call cur_euler(nstate, dt, vol, volold, rates, state)

         e_fast1(k,1:nstate,layer) = state(1:nstate)
        end if
        end do
       end do

!
      if (debug_stranger) then
       do i = 1,nlayers
        do j = 1, nstate
         do k = 1, nkn
          if(layer <= lmax_fast(k)) then
           if(STRANGERS(e_fast1(k,j,i)) .eq. 1) then
!            external node number
             exn_node = ipv(k)
             print *, 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
             print *, 'AQUABC WC after integration:'
             print *, 'Variable ',j,'Cell ',k, exn_node
             print *, 'State variable is NaN'
             print *, 'State=',e_fast1(k,j,i)
             error =1
            end if
          end if ! layer
         end do !k
        end do !j
       end do !i
       if(error .eq. 1) stop
      end if

   return

   end ! subroutine aquabc

!***************************************************************
!***************************************************************
!***************************************************************


            subroutine layer_light(nkn,state,&
                            nstate,          &
                            par,             &
                            nconst,          &
                            layer,           &
                            lmax_fast,       &
                            depth_layer,     &
                            light_surface, alight_bottom)


!==========================================================================
! Calculates incident light after penetration of the given layer of the node
! Only nodes that exist for current layer are processed
! Petras Zemlys 2014.10.2
!
! Inputs:
!       state   -  state variable values for the layer and node
!       nstate  -  number of state variables
!       par     -  model constants
!       nconst  -  number of model constants
!       light_surface    -  incident light to the layer surface
!       layer   -  layer number
!       lmax    -  number of layers for the node
!       node    -  internal node nmber
! Outputs:
!       light_bottom   - incident light to the surface of layer
! Uses:
!       subroutine light_bottom
! Note:
!       Constants array is not used more but still 
!       is not removed from parameters list. Fixme
!=============================================================

      use para_aqua

      implicit none

      integer nkn

      real state(nkn,nstate)
      real par (nconst)

      real light_surface(nkn)  ! incident light to the layer surface
      real alight_bottom (nkn)  ! incident light to the layer bottom

      real DIA_C              (nkn)
      real CYN_C              (nkn)
      real OPA_C              (nkn)
      real FIX_CYN_C          (nkn)

      double precision DIA_C_TO_CHLA
      double precision CYN_C_TO_CHLA
      double precision FIX_CYN_C_TO_CHLA
      double precision OPA_C_TO_CHLA

      real DIA_C_TO_CHLA_fast      (nkn)
      real CYN_C_TO_CHLA_fast      (nkn)
      real FIX_CYN_C_TO_CHLA_fast  (nkn)
      real OPA_C_TO_CHLA_fast      (nkn)
      real CHLA                    (nkn)


      real depth_layer  (nkn)
      !real vol, area  !just for the call dvanode, not used in calculations

      integer lmax_fast(nkn)! number of layers  for the nodes
      integer nstate, nconst, nnode, layer !,   mode, l            
             
            ! Ratios C to Chla
            call para_get_value('DIA_C_TO_CHLA'    ,DIA_C_TO_CHLA)
            call para_get_value('CYN_C_TO_CHLA'    ,CYN_C_TO_CHLA)
            call para_get_value('FIX_CYN_C_TO_CHLA',FIX_CYN_C_TO_CHLA)
            call para_get_value('OPA_C_TO_CHLA'    ,OPA_C_TO_CHLA)
                        
! Following assignments  depend on state variables  numbering!
! Change them if numbering is changed
            where (layer <= lmax_fast(1:nkn))
          ! Phytoplankton concentrations mgC/l:
             DIA_C     =  state(:, 8)
             CYN_C     =  state(:,18)
             OPA_C     =  state(:,19)
             FIX_CYN_C =  state(:,22)
             
             DIA_C_TO_CHLA_fast(1:nkn)     = DIA_C_TO_CHLA
             CYN_C_TO_CHLA_fast(1:nkn)     = CYN_C_TO_CHLA
             FIX_CYN_C_TO_CHLA_fast(1:nkn) = FIX_CYN_C_TO_CHLA
             OPA_C_TO_CHLA_fast(1:nkn)     = OPA_C_TO_CHLA 

             ! Total chorophyl a
            CHLA = ((DIA_C / DIA_C_TO_CHLA_fast) +     &
                   (CYN_C / CYN_C_TO_CHLA_fast)  +     &
                   (OPA_C / OPA_C_TO_CHLA_fast)  +     &
                   (FIX_CYN_C /FIX_CYN_C_TO_CHLA_fast )) * 1.0D3
            end where

           call light_bottom(nkn,lmax_fast,layer, depth_layer,CHLA, &
                             light_surface, alight_bottom)

      end  ! routine layer_light
!***********************************************************
!***********************************************************
!***********************************************************

      subroutine light_bottom(nkn,lmax,layer, DEPTH,CHLA,LIGHT,ALIGHT_bottom)

      ! Calculates light intensity on the bootom of the
      ! layer from light intensity on the surface of the layer and layer depth
      ! Takes  into account only Chl a in this version
      ! Calculates only for nodes that exist for current layer

          use para_aqua

          IMPLICIT NONE

          integer nkn      !number of nodes
          integer layer    !layer number
          integer lmax(nkn)!number of layers for the node

          real DEPTH (nkn) ! depth of layer surface
          real CHLA  (nkn) ! Chl a concentration in the overlaying layer

          real KE_fast   (nkn)
          real XKC_fast  (nkn)
          double precision KE     ! Background light extinction coefficient
          double precision XKC    ! extinction per chlorophyl unit

          real LIGHT (nkn) ! incident light
          real EXTIN (nkn) ! Total extinction coefficient

          real ALIGHT_bottom (nkn)

          ALIGHT_bottom(:) = 0.

          call para_get_value('KE', KE)
          KE_fast = KE
          != 0.8   ! later come as input parmeters background light extinction !fixme
                      ! Should be compatible with value in aquabc_II
          call para_get_value('XKC', XKC)
          XKC_fast = XKC
          != 0.010 ! extinction per chlorophyl unit                           !fixme
                      ! should be compatible with the value in subroutine LIM_LIGHT

!          EXTIN = KE + (8.8D-3 * CHLA) +
!     *            (5.4D-2 * (CHLA ** (2.0D0 / 3.0D0)))
         where( layer <= lmax)

          EXTIN = KE_fast + (XKC_fast * CHLA) ! Extinction coefficient
          ALIGHT_bottom = LIGHT * EXP(-EXTIN * DEPTH)

         end where

      END

!********************************************************************
!********************************************************************

!**************************************************************************
!**************************************************************************
        subroutine SETTLING_TO_THE_NEXT_LAYER            &
            (STATE_VARIABLES    , NUM_VARS   ,           &
             SETTLING_VELOCITIES, DISSOLVED_FRACTIONS,   &
             CELLNO, LAYER, PSTIME,                      &
             SETTLING_RATES)

        implicit none

        !Dummy arguments
        integer, intent(in) :: NUM_VARS, CELLNO, LAYER

      double precision, intent(in)    :: STATE_VARIABLES    (NUM_VARS)
      double precision, intent(in)    :: SETTLING_VELOCITIES(NUM_VARS)
      double precision, intent(in)    :: DISSOLVED_FRACTIONS(NUM_VARS)
      double precision, intent(inout) :: SETTLING_RATES     (NUM_VARS)
      double precision, intent(in)    :: PSTIME

      integer i

        do i = 1, NUM_VARS
          SETTLING_RATES(i) = STATE_VARIABLES(i) *  &
            (1.0D+0 - DISSOLVED_FRACTIONS(i))* SETTLING_VELOCITIES(i)
        end do

      end subroutine SETTLING_TO_THE_NEXT_LAYER

!**************************************************************************
!**************************************************************************
       subroutine cur_euler(nstate,dt,vol,volold,cds,c)

! new c is computed and replaced by new values


       implicit none

       integer nstate	!number of state variables
       real dt			!time step [day]
       real vol		    !volume [m**3]
       real c(nstate)	!state variable [mg/L] == [g/m**3]
       real cds(nstate)	!source term [g/day]
       real massnew
       integer i
       real volold, volnew	!volume [m**3]
       real mass		    !mass [g]
       real mder		    !derivative [g/day]

       !volold = vol
       volnew = vol

       do i=1,nstate
         mass    = c(i) * volold
         mder    = cds(i) * volold + c(i)*(vol-volold)
         massnew = mass + dt*mder
         c(i)    = massnew/volnew
!       c(i) = c(i) + dt*cds(i)
         if(c(i).lt.0.) c(i)=1.e-20
       end do

       end


!********************************************************************
!********************************************************************

      INTEGER FUNCTION STRANGERS(VALUE)
      implicit none
! Cheks for NaN and Inf
! Input is single precision!

      REAL VALUE, BIGNUMBER, RATIO, LLIMIT, ULIMIT

      !print *, 'VALUE=', VALUE

      BIGNUMBER=1.0E30
      LLIMIT = -10000.
      ULIMIT =  10000. !Values for WC
      STRANGERS=0
!       print *,value
      if (isnan(VALUE)) then
       STRANGERS=1
       return
      end if

      RATIO = BIGNUMBER/VALUE

      if(RATIO .eq. 0.) then
       STRANGERS=1
       return
      endif
      if(VALUE .le. LLIMIT .or. VALUE .ge. ULIMIT) then
       STRANGERS=1
       return
      end if


      !print *,'STRANGERS=',STRANGERS

      return
      end


































