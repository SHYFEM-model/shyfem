! aquabc_II to fem interface, version with introduced alcalinity routines.
!
! Produced by Petras 2014
!
! Based on:
!    $Id: bio3d.f,v 1.33 2008-10-10 09:29:54 georg Exp $
!
!    bio3d - EUTRO in SHYFEM
!
!
! revision log of bio3d.v:
!
! 20.06.2003    ggu&dmk new routine for sediments
! 20.08.2003    ggu	new routine bio_av_shell (aver/min/max)
! 03.09.2003    ggu	bug fix for sediments -> not saved (SAVESED)
! 03.09.2003    ggu	new routine check_bio
! 09.09.2003    ggu	call to scal3sh changed -> 3D version
! 19.12.2003    ggu	sediments added, init for taranto
! 14.01.2004    dmk call tsmass per calcolo conserv massa
! 03.03.2004    ggu	decay function for bacteria decad_bio()
! 04.03.2004    ggu	changes from donata integrated
! 05.03.2004    ggu	initialization from file
! 22.03.2004    ggu	change in call to massconc (bug)
! 30.03.2004    ggu	bug fix -> call to confil with nlvdim
! 20.07.2004    dmk	new routine sed_av_shell (aver/min/max)
! 24.08.2004    ggu	new version from donata (jeanmichel)
! 24.08.2004	ggu     new check_es, changes in check_bio
! 24.08.2004	ggu     ivect(8) in bio_av_shell,
! 24.08.2004	ggu     sedload deleted -> substituted by setload_new
! 24.08.2004	ggu	loicz moved to sediment routine
! 25.08.2004    ggu	setsedload moved to sedim routines
! 25.08.2004	ggu	light routines deleted
! 25.08.2004	ggu     new call to eutro0d implemented, read lux
! 25.08.2004	ggu     rluxaux introduced (for ntot>0)
! 17.01.2005    ggu	new horizontal diffusion
! 07.11.2005    ggu     sinking velocity wsink introduced in call to scal3sh
! 17.02.2006    ggu     pass what to subroutines to see calling routine
! 23.03.2006    ggu     ntot eliminated
! 23.03.2006    ggu     changed time step to real
! 18.10.2006    ggu     new routine custom_restime
! 18.10.2006    ggu     introduce bresi,breact,bdecay
! 17.03.2008    ggu     new open boundary routines introduced
! 08.04.2008    ggu     treatment of boundaries slightly changed
! 22.04.2008    ggu     advection parallelized
! 23.04.2008    ggu     call to bnds_set_def() changed
! 09.10.2008    ggu     new call to confop
! 08.05.2014    ggu     bug in call to inicfil for es -> must be inic2fil
! 21.10.2014    ggu     converted to new boundary treatment
!

! aquabc_fem_interface is biogeochemical model specific. Please change line 11 and
! corresponding lines in Make file to switch to different biogeochemical models of AQUABC series

        subroutine ecological_module(it,dt)

! general interface to ecological module

        implicit none

        integer it
        real dt

        call aquabc_II_fem_interface(dt)

        end

! Contents of the rest:
!   subroutine aquabc_fem_interface - Interface between SHYFEM and AQUABC (Aquatic Biochemical Cycling)
!   subroutine setload            - reads point source loads. Not fully tested
!   subroutine check_var          - checking strange values.
!   subroutine check2Dbio         - tests array for NaN and strange values for EUTRO variables. Used by check_var
!   function max_int              - finds maximum element in one dimensional array
!   subroutine inicfil_aquabc     - reads initial values for state variables of WC (also for repeated runs)
!   subroutine inicfils_aquabc    - reads initial values for state variables of BS (also for repeated runs)
!   subroutine calc_time          - not used anymore
!   subroutine print_time         - prints date and time
!   subroutine get_ITOT           - reads light for aquabc. Not used more. Light comes from hydrodynamic module
!   subroutine biotser_init       - initialisation of ASCII output for given nodes (stations)
!   subroutine biotser_write      - writing of ASCII output for given nodes
!   subroutine cur_param_read_ini_wc - intialisation of parameters array for WC model
!   subroutine cur_param_read_ini_bs - intialisation of parameters array for BS model
!   subroutine cur_param_read     - reading WC model parameters(constants)
!   subroutine auquabcini         - does preparations for WC model
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
! 20.07.2004    The structure has been modified to allow dynamic
!               forcings to read from external files
!    07.2006    EUTRO changed by new eutrofication module AQUABC with
!               14 state variables
!    07.2011    Water column eutrophication module changed to ALUKAS_II
!               (Advanced Level nUtrient Kinetics for Aquatic Systems)
! Notes :
!
!                  ccccccccccccccccccccccccc
!                  c    WC state variables c
!                  ccccccccccccccccccccccccc
!                                  STATE VARIABLE      NO   NAME IN ALUKAS
!        ----------------------------------------      --   --------------
!                                   AMMONIUM NITROGEN   1       NH4_N            mg N/l
!                                    NITRATE NITROGEN   2       NO3_N            mg N/l
!                           ORTHOPHOSPHATE PHOSPHORUS   3       PO4_P            mg P/l
!                                    DISSOLVED OXYGEN   4       DISS_OXYGEN      mg O2/l
!     Chemoautotrophic(Nitrification) bacteria carbon   5       NITR_BAC_C       mg C/l
!      Aerobic-anaerobic heterotropic bacteria carbon   6       AER_HET_BAC_C    mg C/l in anaerobic conditions work also as denitrifiers
! Facultative anaerobic heterotrophic bacteria carbon   7       DENITR_BAC_C     mg C/l switched off for the Curonian lagoon, incorrect formulation
!                                       DIATOMS CARBON  8       DIA_C            mg C/l
!                                   ZOOPLANKTON  CARBON 9       ZOO_C            mg C/l
!                                 ZOOPLANKTON NITROGEN 10       ZOO_N            mg C/l derived variable, switched to constant stoichiometry
!                               ZOOPLANKTON PHOSPHORUS 11       ZOO_P            mg C/l derived variable, switched to constant stoichiometry
!                     DETRITUS PARTICULATE ORG. CARBON 12       DET_PART_ORG_C   mg C/l
!                   DETRITUS PARTICULATE ORG. NITROGEN 13       DET_PART_ORG_N   mg N/l
!                 DETRITUS PARTICULATE ORG. PHOSPHORUS 14       DET_PART_ORG_P   mg P/l
!                             DISSOLVED ORGANIC CARBON 15       DISS_ORG_C       mg C/l
!                           DISSOLVED ORGANIC NITROGEN 16       DISS_ORG_N       mg N/l
!                         DISSOLVED ORGANIC PHOSPHORUS 17       DISS_ORG_P       mg P/l
!                                 CYANOBACTERIA CARBON 18       CYAN_C           mg C/l
!                           OTHER PHYTOPLANKTON CARBON 19       OPA_C            mg C/l
!                                    DISSOLOVED SILICA 20       DISS_Si          mg Si/l
!                                      BIOGENIC SILICA 21       BIOG_Si          mg Si/l
!                 NITROGEN FIXING CYANOBACTERIA CARBON 22       FIX_CYN_C        mg C/l
!                           DISSOLVED INORGANIC CARBON 23       INORG_C          mol/m3
!                                      TOTAL ALKALINTY 24       TOT_ALK          mol/m3
!
!          DERIVED VARIABLES SENT TO OUTPUT:
!                                      TOTAL NITROGEN  25
!                                    TOTAL PHOSPHORUS  26
!                                           TOTAL DIN  27
!                                                  pH  28
!
!                  ccccccccccccccccccccccccc
!                  c    BS state variables c
!                  ccccccccccccccccccccccccc
!        NO    NAME IN BS model       STATE VARIABLE
!        --    ----------------------------------------
!        1     SED_NH4N         BS total amonia (in solutes and solids) mg/l of sediments
!        2     SED_NO3N         BS nitrates,                   mg/l in pore water
!        3     SED_DON          BS dissolved organic nitrogen, mg/l in pore water
!        4     SED_PON          BS particulate organic nitrogen
!        5     SED_PO4P         BS dissolved ortophosphates phosphorus
!        6     SED_DOP          BS dissolved organic phophorus
!        7     SED_POP          BS particulate organic phosphorus
!        8     SED_DOXY         BS dissolved oxygen
!        9     SED_DOC          BS dissolved organic carbon
!       10     SED_POC          BS particulate organic carbon
!       11     SED_DSi          BS dissolved silica
!       12     SED_PSi          BS particulate silica
!       13     SED_INORG_C      BS DISSOLVED INORGANIC CARBON      mol/m3
!       14     SED_TOT_ALK      BS TOTAL ALKALINTY                 mol/m3
!       15     SED_SALT         BS SALINITY                        promiles
!
!  DERIVED VARIABLES SENT TO OUTPUT:
!       16  Dissolved  NH4N
!       17  Dissoved   PO4P
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Running:
!  ibio = 0 - ecological module not used                                                                                      c
!  ibio = 1 - WC kinetics without BS kinetics. Settling and                                  c
!             deposition switched on. Settling rates and deposition fraction                 c
!             for each type of sediments are still hardcoded in subroutine aquabc.           c
!  ibio = 2 - WC and BS kinetics. Sediment properties still hardcoded in subroutine aquabc.  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Output:                                        c
!  *.bio - water column variables     (unit iub) c
!  *.bs  - bottom sediments variables (unit iubs)cccccccccccccccccccccccc
!  ASCII output for state  and intermediate variables in observations   c
!  stations (nodes) can be used.                                        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Input:                                                c
! *.str name section aquabc specific files:             c
!        bio       - initial state for WC               c
!        bios      - initial state for BS               c
!        biocon    - constants for WC                   c
!        bioscon   - constants for BS                   c
!        bioaow    - ASCII output control file for WC   c
!        bioaos    - ASCII output control file for BC   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!
!***********************************************
!***********************************************
!***********************************************
      subroutine aquabc_II_fem_interface (dt)

	use mod_diff_visc_fric
	use levels
	use basin
      use aquabc_II_sed_ini

      implicit none

      include 'param.h'          !Geometrical parameters(from SHYFEM)
      include 'aquabc_II.h'      !Main aquabc parameters
      include 'aquabc_II_aout.h' !variables and parameters for output to ASCII file
      ! former common blocks:
      include 'mkonst.h'
      include 'femtime.h'

      !integer it	!time in seconds. Comes from femtime.h now
      !integer idt  !time step in seconds.Comes from femtime.h now
      double precision tdouble
      real dt	    !time step in seconds!495vers real dt remove idt
      real dtday    !time step in days
      real tday, t0, tsec
!-------------------------------------------------------
      ! WC, BS state variables related arrays
      ! nstate - no of WC state vars, defined in aquabc_II.h
      ! nkndim - max no of nodes,     defined in param.h
      ! nlvdim - max no of layers,    defined in param.h

      real e(nlvdim,nkndim,nstate)	      ! WC state vector
      real eload(nlvdim,nkndim,nstate)    ! WC vector of loadings

      real eflux(nlvdim,nkndim,ncsdim)    ! WC state variables selected for for flux calculation
                                          ! (ncsdim max number defined in param.h)
      integer nflux                       ! actual number of state variables for flux calculation
      parameter (nflux = 2)               ! It should go to include file. Fixme.
      real es(noslay,nkndim,nsstate)      ! BS state variables
      save e, es, eload

!----------------------------------------------------

      !loops related:
      integer k,i,l
      integer j
      integer ls
      integer lmax        !max no of layers for the node
      integer nvar



      !real bioarr(nb3dim,0:nbcdim)	!array containing boundary state
      !save bioarr


!     aquabc parameters(constants)
      real par(nconst)     !aquabc model WC parameters array
      save par
      real par_sed(100)    !aquabc model BS parameters array
      save par_sed



!     Arrays passed to aquabc for faster calculations by vectorization
      real e_fast0(nkndim,nstate,nlvdim)     ! initial state vector arranged
                                             !  for faster calculations by vectorization
      real e_fast1(nkndim,nstate,nlvdim)     ! final state vector arranged
                                             !  for faster calculations by vectorization
      integer lmax_fast   (nkndim)           ! max lev. numbers
      real    depth_fast  (nkndim,nlvdim)           ! depth
      real    vel_fast    (nkndim,nlvdim)           ! current velocity
      real    wtemp_fast  (nkndim,nlvdim)           ! water temperature
      real    wind_fast   (nkndim)           ! wind speed
      real    atemp_fast  (nkndim)           ! air temperature
      real    ice_cover_fast(nkndim)         ! ice cover (node area fraction)
      real    sal_fast    (nkndim,nlvdim)    ! water salinity
      real    light_fast  (nkndim)           ! incident light for the layer
      real    vol_fast    (nkndim,nlvdim)
      real    vol_old_fast(nkndim,nlvdim)
      real    area_fast   (nkndim,nlvdim)

      double precision es_fast0(nkndim,noslay,nsstate)     !sediment state vars to pass to kinetics
      double precision es_fast1(nkndim,noslay,nsstate)     !sediment state vars to get from kinetics

      !------------------------------------------------------------------
      ! Arrays for binary output of  state vars and derived vars
      real wc_output(nlvdim,nkndim,noutput) !for state vars  together with derived vars to be written to output files
                                            ! derived vars are placed as nstate+1, ... elements
      save wc_output                        !needed for chlorophyl a

      real sed_output(NOSLAY,nkndim,nsoutput)              ! BS output
      save sed_output
!-----------------------------------------------
      integer nlayers  ! maximum number of levels in WC (calculated)

      integer ilhkv_sed(nkndim)
      save ilhkv_sed     ! Needed for compatibility of WC and BS writing in biotser_write

      real elaux(nstate) ! WC state vars  loads past to aquabc

      ! controls calls to routines
      character*10 what,whataux
      character*2 whatn

      integer ibio
      integer id
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

      integer idbio(nbcdim)
      save idbio

      real tstot(nstate)   !for mass test
      real tsstot(nsstate)

      integer icall
      save icall

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
      real getpar
      integer iround


      logical bsedim  !true -process BS
      save bsedim

      logical bcheck

      logical bresi,breact,bdecay

      integer ie,ii
      integer kspec



      real windspeed,tempair
      real rh,twetbulb,cloudc,pr !not used, just for copmatibility
      real ice_cover

      real mass
      real wsink

      integer iespecial,inspecial
      save iespecial,inspecial

      integer iub
      save iub
      integer iubs
      save iubs
      integer ia_out(4)
      save ia_out

      logical has_output,next_output

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
       save    NBIOTS     !declared in aquabc_aout
       save    BIOTSFUN  !declared in aquabc_aout
       save    BIOTSNOD  !declared in aquabc_aout

! Variables used for BS state variables output to ASCII files
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
       save    NBIOTS_sed     !declared in aquabc_aout
       save    BIOTSFUN_sed  !declared in aquabc_aout
       save    BIOTSNOD_sed  !declared in aquabc_aout

!      variables used for intermediate WQ results(processes) output(processes) to ASCII files
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       integer dg_count	        !do loop parameter
       real    dg  (nlvdim,nkndim,nstate,NDIAGVAR)

       save    NDGTS   !declared in aquabc_aout
       save    DGTSFUN !declared in aquabc_aout
       save    DGTSNOD !declared in aquabc_aout

!      variables used for intermediate bottom sediments results(processes) output(processes)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       real    dg_sed  (NOSLAY,nkndim,nsstate,NDIAGVAR) ! diagnostic data

       save    NDGTS_sed   !declared in aquabc_aout
       save    DGTSFUN_sed !declared in aquabc_aout
       save    DGTSNOD_sed !declared in aquabc_aout


      integer ulogbio,ifileo
      save ulogbio

      integer ifirst_sed   !indicator of first call for BS, 1 - if first, 0 - if not
      data ifirst_sed /1/
      save ifirst_sed

      double precision dtime0

!     variables and functions for light management
      real ALIGHT_ALUKAS ! function for calculation light for the next layer


      integer STRANGERS   ! Function checks for NaNs and Inf
      integer idump       !indicator of dump for repeated runs

	real alight,aair        !for Rasa

       idump = 1
!      idump = 0

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
! !                 8     9     10    11     12   13    14
!      *          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01/
!               15   16   17      18   19   20    21

!  Initialisation of BS state variables.  Sandy sed are assumed.! Fixme!
!ccccccccccccccccccccccccccccccccccccccccccccccc
!data esinit  /0.1,0.1,0.1,0.1,0.02,0.1,0.1,0.1,0.1,0.1,1.,0.1/ old initial conditions
!                   1      2      3      4      5      6      7
!       data esinit /0.5000,0.0319,1.000,1000.0,0.1800,0.1400,10.0,
!      *      0.0141,105.0,1400.0 ,3.500,10.000/
! c            8       9     10      11     12
! c                   1      2      3      4      5      6      7
!       data esinit2 /1.2000,0.0432,0.707 ,1500.,0.160,0.1200,130.0,
!      *      0.000001,90.0  ,1300.,8.000,130.0/
! c            8       9     10      11     12
! c                   1      2      3      4      5      6      7
!       data esinit3 /1.5000,0.0509,0.5000,1500.,0.150,0.4000,165.0,
!      *      0.000001,90.0  ,1000.,8.000,165.0/
! c            8       9     10      11     12

!------------------------------------------------------------------
!       data icall /0/
!------------------------------------------------------------------
        bresi = .false.
        breact = .true.
        bdecay = .false.


!       bcheck = .false.
        bcheck = .true.       !.true. if state variables are checked

        what = 'lagvebio'
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

        ibio = iround(getpar('ibio'))
        if( ibio .eq. 0 ) icall = -1
        if(ibio .gt. 2 .or. ibio .lt. 0) icall = -1
        if( icall .le. -1 ) return

        icall = 1

        call print_time

        if(ibio .eq. 2) bsedim = .true.
        if(ibio .eq. 1) bsedim = .false.
!       ---------------------------------------------
!       parameters for transport/diffusion resolution
!       ---------------------------------------------


          rkpar=getpar('chpar')
          difmol=getpar('difmol')

          do i=1,nkn
           ilhkv_sed(i)= NOSLAY
          enddo


!-----------------------------------------------------
!     Open file for bio variables strange values check
!-----------------------------------------------------
         if(bcheck) then
          ulogbio = ifileo(60,'logbio.txt','f','u')
          if( ulogbio .le. 0 ) then
             write(6,*) 'Cannot open/create file logbio.txt'
             stop
            else
             write(6,*) 'Aquabc_fem: File logbio.txt on unit: ',
     +        ulogbio
             write(ulogbio,*)
     +          '      STRANGE VALUES FOR EUTRO VARIABLES'
             write(ulogbio,*)
     +          'Var.No Value  Node Level Time  Checking moment'
          end if
         end if

!     --------------------------------------------------
!     initialize state variables with einit
!     --------------------------------------------------
      einit(1:nstate) = 0.

      do k=1,nkn    !loop on nodes
            lmax = ilhkv(k)
            do l=1,lmax
             do i=1,nstate
              e(l,k,i) = einit(i)
             end do
           end do
      end do

!     Initialise ascii diagnostics, light for water column

          dg(:,:,:,:)       = 0.
          wc_output(:,:,:)  = 0.
          sed_output(:,:,:) = 0.

!----------------------------------------------------------
!         initialize WQ state variables from external file
!         BS initialisation moved after call aquabcini
!----------------------------------------------------------

          call inicfil_aquabc('bio',e,nstate)


!------------------------------------------------------------
!         set boundary conditions for all WC state variables
!------------------------------------------------------------

         nintp=2
         call get_first_time(itanf)

         !call bnds_init(what,bio2dn,nintp,nstate,nb3dim,bioarr,ebound)	!new_section
         dtime0 = itanf
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
!        Initilises sediment properties and some parameters for BS
!        and  settling velocities and dissolved fractions for WC variables
!        see module aquabc_II_sed_ini
!------------------------------------------------------------------------       
         call  sed_properties_ini(bsedim)          

!        Writing initial conditions to text file for defined nodes

         idtc = nint(getpar('idtcon'))
         itmc = nint(getpar('itmcon'))

         !print *,'idtc=',idtc,'itmc=',itmc

         print *,'WRITING WC STATE VARIABLE INITIAL',
     *              ' VALUES FOR SELECTED NODES'

         do l = 1,nlvdim
          do k = 1,nkn
           do i=1,nstate
             wc_output(l,k,i) = e(l,k,i)  ! derived variables are not calculated yet here (zeros)
           end do                         ! It will be done in pelagic kinetic module together
          end do                          ! Together with calculation of derivatives for the first time step
         end do

         call biotser_write(1, 'wc',
     *                 wc_output, noutput, nstate, dg, NDIAGVAR,
     *                 ilhkv,nlvdim,
     *                 itmc,idtc,
     *                 NBIOTSMX,NBIOTS,BIOTSNOD,BIOTSFUN,
     *                 NDGTSMX, NDGTS, DGTSNOD, DGTSFUN)



            
! ADDITIONAL PREPARATION FOR BS ***************************************
         if (bsedim) then

!         initialise ascii diagnostics

          dg_sed(:,:,:,:) = 0.
!         ---------------------------------
!         Initialisation BS state variables
!         ---------------------------------
          call inicfils_aquabc('bios',es,nsstate)

!         -----------------------------
!         Temporary initialize state BS variables. Now read from file.
!         -----------------------------

!            do k=1,nkn
!             do i=1,nsstate
!
!             es(1,k,i) = esinit(i)
!             es(2,k,i) = esinit2(i)
!             es(3,k,i) = esinit3(i)
!
!             end do
!            enddo


        if (ifirst_sed .eq. 1) then
            call  sed_recalc_ini(es, sed_output)
            ifirst_sed =0
        end if


!        do i=1,nkn
!         print *, 'tikr pries write ','k= ',i,' ',e(1,i,1:nstate)
!        end do

           print *,'WRITING BS STATE VARIABLE INITIAL',
     *              ' VALUES FOR SELECTED NODES'


         call biotser_write(1, 'bs',sed_output,nsoutput,nsstate,
     *              dg_sed,NDIAGVAR_sed,
     *              ilhkv_sed,NOSLAY,
     *              itmc,idtc,
     *              NBIOTSMX_sed,NBIOTS_sed,BIOTSNOD_sed,BIOTSFUN_sed,
     *              NDGTSMX_sed, NDGTS_sed, DGTSNOD_sed, DGTSFUN_sed)


         end if
!    end of preparation for BS


!------------------------------------------------------------
!         initialize output to binary files
!------------------------------------------------------------

          call init_output('itmcon','idtcon',ia_out)

          if( has_output(ia_out) ) then
            call open_scalar_file(ia_out,nlv,nstate,'bio')
            iub = ia_out(4)
            write(6,*) 'Binary output for WC model initialized...'
            if( bsedim ) then
              call open_scalar_file(ia_out,1,nsstate,'bs')
              iubs = ia_out(4)
            end if
            write(6,*) 'Binary output for BS model initialized...'
          end if


      end if  !

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------
! end of initialisation
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------
! normal call (every time step)
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!        do i=1,nkn
!         print *, 'tikr pries normal ','k= ',i,' ',e(1,i,1:nstate)
!        end do

      what = 'lagvebio'

      kspec = 4350
      kspec = -1

      wsink = 0.

!-------------------------------------------------------------------
!    time management
!-------------------------------------------------------------------

      t0    = 0.
      dtday = dt / 86400
      tsec  = it
      tday  = it / 86400. + t0  !time in days, FEM 0 is day t0                             

      if( bcheck ) call check_var('BEFORE aquabc',it,ulogbio,e,es) 
      
!-----------------------------------------------------------------------
!-------------------------------------------------------------------
!     data preparation loop on nodes for biological reactor
!-------------------------------------------------------------------
!----------------------------------------------------------------------
      do k=1,nkn   !loop on nodes

!       print *, 'tikrinam regione ','k= ',k,' ',e(1,k,1:nstate)

         lmax = ilhkv(k)  !maximum number of levels for the node
         lmax_fast(k) = lmax

! gets meteorological variables 
         call meteo_get_heat_values(k,ITOT,tempair,rh,twetbulb,
     +                         windspeed,cloudc,pr)  !New routine produced by Georg

         atemp_fast  (k) = tempair
         wind_fast   (k) = windspeed

         light_fast  (k)  = ITOT
! now ITOT is instanteneous light not daily averaged
! as needed for new Ali-Smith routine     

!        routines to introduce fraction of ice cover:
!        get_ice(k,ice_cover)   - for one node
!        get_ice_all(ice_cover) - for all nodes

         call get_ice(k,ice_cover)
         ice_cover_fast(k) =  ice_cover 

! code for Rasa ========================================

         !ice_cover is fraction of area covered with ice
         !   0 no ice cover
         !   1 completely ice covered

         !alight is fraction of light blocked: 
         !   0 all light is transmitted
         !   1 all light is blocked when ice covered

         !aair is re-aeration coefficient (should be either 0 or 1)
         !   0 re-aeration is not blocked by ice
         !   1 re-aeration is blocked by ice

         alight = 1     !fraction of light blocked
         aair = 1       !re-aeration blocked

         light_fast  (k)  = ITOT * (1.-ice_cover*alight)
         ice_cover_fast(k) = ice_cover * aair

! code for Rasa ========================================

         FDAY = 1. !used only for old equations compatibility 

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

         if (isnan(depth).or.depth.eq.0)
     +        then
          print *,
     +   'aquabc_II_fem_interface: Depth is NaN or zero:', depth,
     +         'on level: ', l, 'on node: ',k
          stop
         end if

        if (isnan(vol).or.vol.eq.0)
     +        then
         print *,
     +   'aquabc_II_fem_interface: Volume is NaN or zero:', vol,
     +         'on level: ', l, 'on node: ',k
         stop
        end if

         if (isnan(area).or.area.eq.0)
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

            id = 1000*k+l !?

! Changing WC old values  to new in e_fast0 for the current step calculations
            do i=1,nstate
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


!---------------------------------------------------------------------
!       --------------------------------------------------------------
!       end of loop on nodes for arrays preparation for call to AQUABC
!       --------------------------------------------------------------
!---------------------------------------------------------------------

            if(bsedim) then
             do ls=1,noslay
              do k=1,nkn
               do i=1,nsstate
                es_fast0(k,ls,i) = es(ls,k,i)
               end do
              end do
             end do
            end if


            nlayers = maxval(lmax_fast(1:nkn))

!          loads are not corrected for this version of the model, fixme

            !print *,'aquabc_fem_interface: Entering AQUABC'


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
     *                  sed_output, nsoutput,
     *                  dg_sed, NDIAGVAR_sed)

             !print *,'aquabc_fem_interface: Leaving AQUABC'

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
            end if

         !print *,'After call to aquabc: e_fast1, L1', e_fast1(625,:,1)
         !print *,'After call to aquabc: e:, L1', e(1,625,:)
         !print *,'After call to aquabc: e_fast1, L2', e_fast1(625,:,2)
         !print *,'After call to aquabc: e:, L2', e(2,625,:)



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

          call tsmass (e(1,1,i),1,nlvdim,tstot(i)) !mass control


      end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

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

!----------------------------------------------------
! Commented temporary for experiment with sediments (check bio3d!). es is 3d now. Is this check necessary?
!      do i=1,nsstate
!          call scalmass(es(1,i),0.1,tsstot(i))   !mass ctrl sed
!      end do



!      -------------------------------------------------------------------
!      WRITE FLUXES
!      -------------------------------------------------------------------
       ! does not work in 3d  in parallel mode in this version. fixme

       !eflux(1:nlvdim,1:nkn,1) = wc_output(1:nlvdim,1:nkn,nstate+1)   ! Total N
       !eflux(1:nlvdim,1:nkn,2) = wc_output(1:nlvdim,1:nkn,nstate+2)   ! Total P
!      !
       !call fluxes_aquabc(it, nflux, eflux)

        !check nflux value in declaration(param.h). Results to *.csc in g/s



!    -------------------------------------------------------------------
!    WRITE OF RESULTS (FILE BINARY BIO AND TEXT)
!    -------------------------------------------------------------------

!          write(17,*) 'AFTER scal3sh:','time=',it
!          write(17,1000) (ipv(k),(e(1,k,i),i=1,9),k=1,nkn)
!1000      format((I4,1x,9(F8.4,1x)))


!       WATER COLUMN OUTPUT 
        if( next_output(ia_out) ) then
         ia_out(4) = iub
         do i=1,noutput
           id = 70 + i
           call write_scalar_file(ia_out,id,nlvdim,wc_output(1,1,i))
         end do
        end if 
!            print *,'WRITING WC STATE VARIABLE ',
!     *              ' VALUES FOR SELECTED NODES'
        call biotser_write(0, 'wc',wc_output, noutput, nstate,
     *                  dg, NDIAGVAR,
     *                  ilhkv,nlvdim,
     *                  itmc,idtc,
     *                  NBIOTSMX,NBIOTS,BIOTSNOD,BIOTSFUN,
     *                  NDGTSMX,NDGTS, DGTSNOD, DGTSFUN)

         if(idump.eq.1 .and. it .eq. itend) then
            call  dump_aquabc('dump_wc.dat',e,
     +                 nkndim,nkn,nlvdim,nstate)
         end if

!          SEDIMENT OUPUT
          if( bsedim ) then

           if( next_output(ia_out) ) then
             ia_out(4) = iub
            do i=1,nsoutput
             id = 100 + i
             call write_scalar_file(ia_out,id,noslay,
     *                              sed_output(1,1,i))
            end do
           end if

!           print *,'AQUABC_FEM_INTERFACE: writing bs state vars ',
!     *              ' values for selected nodes'

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
     +                 nkndim,nkn,noslay,nsstate)
           end if !idump

          end if !bsedim



!--------------------------------------------------------------
!--------------------------------------------------------------
!     Averages, min,max
!--------------------------------------------------------------
!--------------------------------------------------------------
!     call bio_av_shell(e)		!aver/min/max of state vars
!     call sed_av_shell(es)		!aver/min/max of sed var

      if( bcheck ) call check_var('AFTER advection',it,ulogbio,e,es)

!    -------------------------------------------------------------------
!    debug output
!    -------------------------------------------------------------------

!    -------------------------------------------------------------------
!    end of aquabc_fem_interface
!    -------------------------------------------------------------------

      end

c*************************************************************

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


	use evgeom
	use basin

	implicit none

      include 'param.h'
	integer nstate
	parameter(nstate=9)

!     MODIFIED BY ALI
!     CORPI, 15 July 2004
!     TAKE CARE
!     eload(3,neldim,nstate) -----> eload(nlvdim,nkndim,nstate)
	real eload(nlvdim,nkndim,nstate)


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

c         extern to intern
		call n2int(nodes,karee,berror)

	    if( berror) stop 'error stop: loading'

		icall = 1

	end if

	if((ifirst.eq.0).and.(lftype.eq.1)) then

    1     continue
C         Read the RECORD 7 of LOADING FILE TYPE 1
          irec = 7

          read(nb, 5030, err=98) itload
		  write(6,*) ''
		  write(6,*) 'RECORD 7 OF THE LOADING FILE READ'

		  preitl = itload


C         ADDED BY ALI
C         CORPI, 22 July 2004
C
C         The following if structure overjumps loading time interval
C         before simulation start
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

C                 Read the RECORD 8 of LOADING FILE TYPE 1
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

C         Read the RECORD 8 of LOADING FILE TYPE 1
          irec = 8

		  do j = 1, nimmis
	        read(nb, 5020, err=94)(areaload(j,ivar), ivar=1,nstate)
		  end do

		  write(6,*) 'RECORD 8 OF THE LOADING FILE READ'

		  preitl = itload
		  inext = 1

C         Read the RECORD 7 of LOADING FILE TYPE 1 - NEXT TIME INTERVAL
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

C     CHECK IF TIME FOR THE NEXT LOADING INTERVAL
      if((((it+idt).ge.itload).and.(ilast.eq.0)).and.lftype.eq.1) then
	    write(6,*) 'NEW LOADING INTERVAL STARTING NEXT TIME STEP'
		  inext = 0
	end if


C     READ LOADING DATA FOR LOADIND FILE TYPE 2
      if(lftype.eq.2) then

	    do j = 1, nimmis

	        if(ifirs2(j).eq.0) then

    2             continue

C                 Read the RECORD 8 of LOADING FILE TYPE 2
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

C                 Read the RECORD 8 of LOADING FILE TYPE 2
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

C             CHECK IF TIME FOR THE NEXT LOADING INTERVAL
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

c	rlfact = getpar('flbio')
	rlfact = 1.


C     MODIFIED BY ALI
C     CORPI, 23 July 2004
C     ERROR ABOUT TAKING ONE DAYS LOADING FOR ECAH TIME STEP FIXED
C
C     fact = (rlfact*kgs/litri)
C     fact = -------> (rlfact*kgs/litri) * (idt/86400.)

C     [kg/m**3] -----> [mg/l] ------> kg/day ---------> (kg/s)*idt
	fact = (rlfact*kgs/litri) * (idt/86400.)



c     intialize

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

c     compute total volume for all areas given -> store in volaux

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

c     compute and set loading in eload [g/(m**3 day)] == [mg/(l day]

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ia = aree(k)
	    vol = volaux(ia)
	    do i=1,nstate
c	 Corrected by Petras to avoid division by zero:
		  if( ia .le. 0 ) then
		    load = 0.
	    else
		    load = fact * areaload(ia,i) / vol
	    endif

C           MODIFIED BY ALI
C           CORPI, 15 July 2004
C           TAKE CARE
C           eload(ii,ie,i) -----> eload(1,k,i)
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

c*************************************************************
c*************************************************************

	subroutine check_var(title,it,ulog,e,es)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c checks WC and BS vars for strange values

C Uses:
C   subroutine check2Dbio - cheks for strange values for given ranges
C                           and variable numbers
C Developped by G.Umgiesser
C Modified for use for biogeochemical variables by P.Zemlys
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

      include 'param.h'
      include 'aquabc_II.h'

      integer ulog,it



      integer nslayer



      character*(*) title
      real e(nlvdim,nkndim,nstate)	        !state vector
      real es(NOSLAY,nkndim,nsstate)		!sediment state variables


        character*16 text
        integer i

        nslayer = NOSLAY

	!write(6,*) 'check_var: ',title

        text = '!!! BIO CHECK'
C Changed by Petras 11 January 2005 to pass number of variable
       do i=1,nstate
        call check2Dbio(it,ulog,nlv,nkn,70+i,e(1,1,i),
     *                     -0.001,25.,text,title)
       end do

		do i=1,nsstate
        call check2Dbio(it,ulog,nslayer,nkn,100+i,es(1,1,i),
     *                     -0.001,100000.,text,title)
        end do

      end

c***************************************************************

      subroutine check2Dbio(it,ul,nlv,n,nvar,a,vmin,vmax,textgen,text)

c tests array for nan and strange values for EUTRO variables
c Made from check2Dr by P.Zemlys 11 January 2005

	use basin

        implicit none

        include 'param.h'

        integer nlv,n
        integer nvar         !state variable number
        real a(nlvdim,1)     !array of variable values for levels and nodes
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

        end

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

	end

c*********************************************************

       subroutine inicfil_aquabc(name,var,nvar)
c initializes variables nodal value for water column from file

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'

        character*(*) name           !name of variable
        real var(nlvdim,nkndim,1)    !variable to set for water column (name='bio')

        real varval
        integer nvar
        integer ftype               !1-  homogeneous initial cond. for WC
                                    !2 - heterogeneous initial conditions.
                                    !3 - homogeneous layers
                                    !4 - heterogenous for repeated runs


        character*80 file
        integer nb,irec
        integer nkk,lmax
        integer l,k,i,j
        integer ivars,ivar
        real val
        real rlaux(nlvdim) !absolute depth of bottom of layer i
                           !nlvdim - total number of vertical levels

        integer ifileo
        real varval_array(2000)
c-------------------------------------------------------
c get file name
c-------------------------------------------------------

        call getfnm(name,file)

        if( file .eq. ' ' ) return      !nothing to initialize

c-------------------------------------------------------
c open file
c-------------------------------------------------------

C     MODIFIED BY ALI AND PETRAS
C     CORPI, 16 July 2004C
C     nb = ifileo(55,file,'unform','old')
C      -------> nb = ifileo(55,file,'f','old')
C
        nb = ifileo(55,file,'f','old')
        if( nb .le. 0 ) goto 97

c-------------------------------------------------------
c read first record
c-------------------------------------------------------

      ivar = 0

        irec = 1

C     MODIFIED BY ALI AND PETRAS
C     CORPI, 16 July 2004

C
C     read(nb,err=90) nkk,lmax,ivars
C      -----> read(nb, 5000, err=90) nkk,lmax,ivars
C
C Skip 4 lines. Added by Petras 12-12-2004
        read(nb,*)
        read(nb,*)
        read(nb,*)
        read(nb,*)
C Read control information
        read(nb, *, err=90) nkk,lmax,ivars,ftype
 5010 FORMAT(4I5)

      if( nkk .ne. nkn)  goto 99

      if(lmax .gt. nlvdim ) goto 99

      if( ivars .ne. nvar ) goto 96
      if(ftype .ne. 1 .and. ftype .ne. 2
     +  .and. ftype .ne. 3 .and. ftype .ne. 4) goto 91

C********************************************************
C FILE TYPE 2 (spatialy heterogeneous initial conditions)
C********************************************************
      if(ftype .eq. 2) then
c-------------------------------------------------------
c read second record (only if lmax > 0)
c-------------------------------------------------------             !

        if( lmax .gt. 1 ) then          !changed from 0 to 1 (5.3.2004) LMAX
          irec = 2
C       MODIFIED BY ALI AND PETRAS
C       CORPI, 16 July 2004

C
C
C       read(nb, err=90) (rlaux(l),l=1,lmax)
C       -----> read(nb, 5020, err=90) (rlaux(l),l=1,lmax)
C
          read(nb, 5020, err=90) (rlaux(l),l=1,lmax)
 5020   FORMAT(8F10.0)

          do l=1,lmax
            if( hlv(l) .ne. rlaux(l) ) goto 98
          end do
        else
          lmax = 1
        end if

c-------------------------------------------------------
c read data records
c-------------------------------------------------------

        irec = 3

        do ivar=1,nvar

C     MODIFIED BY ALI AND PETRAS
C     CORPI, 16 July 2004C
C     read(nb, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
C --> read(nb, 5030, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
C

        read(nb, 5030, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
 5030 FORMAT(F5.0)

      end do

c-------------------------------------------------------
c initialize the other levels if only surface is given
c-------------------------------------------------------

        do ivar=1,nvar
          do k=1,nkn
            val = var(lmax,k,ivar)
            do l=lmax+1,nlvdim
              var(l,k,ivar) = val
            end do
          end do
        end do

      end if ! end of ftype=2

C******************************************************
C FILE TYPE 1 (spatialy homogeneous initial conditions)
C******************************************************
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

C******************************************************
C FILE TYPE 3 spatialy homogeneous initial conditions for each level
C******************************************************
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
C******************************************************
C FILE TYPE 4 dumped state for repeated runs
C******************************************************
      if(ftype .eq. 4) then

       do i = 1, nvar
          do j = 1, nkn
              read(nb,*,err=90) (var(k,j,i),k=1,lmax)
          end do
       end do

      end if !end of type 4

c-------------------------------------------------------
c reading done -> close file
c-------------------------------------------------------
 1001 continue
        close(nb)



c-------------------------------------------------------
c end of initialisation
c-------------------------------------------------------

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
        write(6,*) 'parameters are not the same from init file ',file
        write(6,*) 'nkn, lmax from file  : ',nkk,lmax
        write(6,*) 'nkn, lmax from model : ',nkn,nlvdim
        stop 'error stop inicfil_aquabc'
        end

c*******************************************************************

c*********************************************************

       subroutine inicfils_aquabc(name,var,nvar)
c initializes variables nodal value  for bottom sediment from file

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'
        include 'aquabc_II.h'

        character*(*) name           !name of variable

        real var(noslay,nkndim,1)   !variable to set for bottom sediments (name='bios')
        real varval
        integer nvar
        integer ftype               !1-  homogeneous initial cond. for WC
                                    !2 - heterogeneous initial conditions.
                                    !3 - homogeneous layers for BS
                                    !4 - heterogenous for repeated runs


        character*80 file
        integer nb,irec
        integer nkk,lmax
        integer l,k,i,j
        integer ivars,ivar
        real val
        real rlaux(noslay) !absolute depth of bottom of layer i
                           !nlvdim - total number of vertical levels

        integer ifileo
        real varval_array(2000)
c-------------------------------------------------------
c get file name
c-------------------------------------------------------

        call getfnm(name,file)

        if( file .eq. ' ' ) return      !nothing to initialize

c-------------------------------------------------------
c open file
c-------------------------------------------------------

C     MODIFIED BY ALI AND PETRAS
C     CORPI, 16 July 2004C
C     nb = ifileo(55,file,'unform','old')
C      -------> nb = ifileo(55,file,'f','old')
C
        nb = ifileo(55,file,'f','old')
        if( nb .le. 0 ) goto 97

c-------------------------------------------------------
c read first record
c-------------------------------------------------------

      ivar = 0

        irec = 1

C     MODIFIED BY ALI AND PETRAS
C     CORPI, 16 July 2004

C
C     read(nb,err=90) nkk,lmax,ivars
C      -----> read(nb, 5000, err=90) nkk,lmax,ivars
C
C Skip 4 lines. Added by Petras 12-12-2004
        read(nb,*)
        read(nb,*)
        read(nb,*)
        read(nb,*)
C Read control information
        read(nb, *, err=90) nkk,lmax,ivars,ftype
 5010 FORMAT(4I5)

      if( nkk .ne. nkn)  goto 99


       if(lmax .gt. noslay ) goto 99


      if( ivars .ne. nvar ) goto 96
      if(ftype .ne. 1 .and. ftype .ne. 2
     +  .and. ftype .ne. 3 .and. ftype .ne. 4) goto 91

C********************************************************
C FILE TYPE 2 (spatialy heterogeneous initial conditions)
C********************************************************
      if(ftype .eq. 2) then
c-------------------------------------------------------
c read second record (only if lmax > 0)
c-------------------------------------------------------             !

        if( lmax .gt. 1 ) then          !changed from 0 to 1 (5.3.2004) LMAX
          irec = 2
C       MODIFIED BY ALI AND PETRAS
C       CORPI, 16 July 2004

C
C
C       read(nb, err=90) (rlaux(l),l=1,lmax)
C       -----> read(nb, 5020, err=90) (rlaux(l),l=1,lmax)
C
          read(nb, 5020, err=90) (rlaux(l),l=1,lmax)
 5020   FORMAT(8F10.0)

          do l=1,lmax
            if( hlv(l) .ne. rlaux(l) ) goto 98
          end do
        else
          lmax = 1
        end if

c-------------------------------------------------------
c read data records
c-------------------------------------------------------

        irec = 3

        do ivar=1,nvar

C     MODIFIED BY ALI AND PETRAS
C     CORPI, 16 July 2004C
C     read(nb, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
C --> read(nb, 5030, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
C

        read(nb, 5030, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
 5030 FORMAT(F5.0)

      end do

c-------------------------------------------------------
c initialize the other levels if only surface is given
c-------------------------------------------------------

        do ivar=1,nvar
          do k=1,nkn
            val = var(lmax,k,ivar)
            do l=lmax+1,noslay
              var(l,k,ivar) = val
            end do
          end do
        end do

      end if ! end of ftype=2

C******************************************************
C FILE TYPE 1 (spatialy homogeneous initial conditions)
C******************************************************
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

C******************************************************
C FILE TYPE 3 spatialy homogeneous initial conditions for different  levels
C******************************************************
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

C******************************************************
C FILE TYPE 4 dumped state for repeated runs
C******************************************************
      if(ftype .eq. 4) then

       do i = 1, nvar
          do j = 1, nkn
              read(nb,*,err=90) (var(k,j,i),k=1,lmax)
          end do
       end do

      end if !end of type 4
c-------------------------------------------------------
c reading done -> close file
c-------------------------------------------------------
 1001 continue
        close(nb)



c-------------------------------------------------------
c end of initialisation
c-------------------------------------------------------

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
        write(6,*) (hlv(l),l=1,lmax)
        write(6,*) (rlaux(l),l=1,lmax)
        stop 'error stop inicfils_aquabc'
   99   continue
        write(6,*) 'parameters are not the same from init file ',file
        write(6,*) 'nkn, lmax from file  : ',nkk,lmax
        write(6,*) 'nkn, lmax from model : ',nkn,noslay
        stop 'error stop inicfils_aquabc'
        end

c*******************************************************************
c*******************************************************************

!         subroutine print_time
!
! c prints stats after last time step
!
!         implicit none
!         character  date*8, time*10           !Added by Petras 11.09.2004
!
!
!
!
!
!          write(6,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
! 	   call date_and_time(date,time)
!          write(6,*) 'DATE: ', date
!          write(6,*) 'TIME: ', time
!          write(6,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
!          end
!
c*******************************************************************
c*******************************************************************

      subroutine calc_time (startdate,starttime)

c prints stats after last time step

	implicit none
        character  date*8, time*10           !Added by Petras 11.09.2004
        character  startdate*8, starttime*10

	include 'femtime.h'

	write(6,1035) it,niter
 1035   format(' program stop at time =',i10,' seconds'/
     +         ' iterations = ',i10)

         write(6,*)
         write(6,*) 'START DATE: ', startdate !Added by Petras 11.09.2004
         write(6,*) 'START TIME: ', starttime !Added by Petras 11.09.2004
         write(6,*)
	   call date_and_time(date,time)
         write(6,*) 'END DATE: ', date !Added by Petras 11.09.2004
         write(6,*) 'END TIME: ', time !Added by Petras 11.09.2004
         end



c******************************************************************
c******************************************************************
c******************************************************************

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
      implicit none
      include 'param.h'
	include 'femtime.h'

      real e(nlvd,nkndim,noutput)	!state variables
      real dg  (nlvd,nkndim,nstate,NDIAGVAR) !intermediate variables


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
      FMT_many   = '(I15,I5,' // noutput_char   // '(1x,G13.4))'
      FMT_many_d = '(I4.4,2I2.2,1x,3I2.2,I5,'   //
     +                             noutput_char // '(1x,G13.4))'
      FMT_1      = '(I15,'    // noutput_char   // '(1x,G13.4))'
      FMT_1_d    = '(I4.4,2I2.2,1x,3I2.2,'      //
     +                            noutput_char  // '(1x,G13.4))'

      write(NDIAGVAR_char,'(i2)') NDIAGVAR

      FMT_diag_many   = '(I15,I5,' // NDIAGVAR_char //'(1x,G13.6))'
      FMT_diag_many_d = '(I4.4,2I2.2,1x,3I2.2,I5,'  //
     +                                NDIAGVAR_char //'(1x,G13.6))'
      FMT_diag_1      = '(I15,'    // NDIAGVAR_char //'(1x,G13.6))'
      FMT_diag_1_d    = '(I4.4,2I2.2,1x,3I2.2,'     //
     +                                NDIAGVAR_char //'(1x,G13.6))'



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
!      nconst  -  number of parameters requiried
!      Note:
!           Number of parameters in this routine is hardcoded by npar
!      Called:
!           In CUR_PARAM_READ

      integer npar !existing number of parameters in the routine
      !parameter (npar = 226)
      parameter (npar = 229)	!ggu
      character*80 bname(npar)
      real par(npar)

!     execution part

      if( npar. ne. nconst) then
        print *, 'CUR_PARAM_READ_INI_WC:'
        print *, '  wrong asked number of params'
        print *, '  asked: ',nconst
        print *, '  existing: ',npar
      end if

!     WC kinetics parameter names. These line should be corrected
!     when new parameter is introduced or the name of existing
!     parameter is changed.
      bname(  1) =                              'K_A'  !  1! Aeration coefficient (if negative calculates internally)
      bname(  2) =                        'THETA_K_A'  !  2! Temperature correction factor for aeration
      bname(  3) =                               'KE'  !  3! Background extinction coefficient
      bname(  4) =                              'XKC'  !  4! Light extinction per chlorophyl unit,( mcg Chla/l/m)
      bname(  5) =                            'PHIMX'  !  5! Quantum yield const. mg C/mole photon
      bname(6  ) =               'KG_CHEM_AUT_BAC_20'  !6  ! Chemoautotrophic bacteria Growth rate
      bname(7  ) =          'EFF_CHEM_AUT_BAC_GROWTH'  !7  ! Chemoautotrophic bacteria growth efficiency
      bname(8  ) =            'THETA_KG_CHEM_AUT_BAC'  !8  ! Chemoautotrophic bacteria Temperature correction for growth rate
      bname(9  ) =               'KR_CHEM_AUT_BAC_20'  !9  ! Chemoautotrophic bacteria Respiration rate
      bname(10 ) =            'THETA_KR_CHEM_AUT_BAC'  !10 ! Chemoautotrophic bacteria Temperature correction for respiration rate
      bname(11 ) =               'KD_CHEM_AUT_BAC_20'  !11 ! Chemoautotrophic bacteria Mortality rate
      bname(12 ) =            'THETA_KD_CHEM_AUT_BAC'  !12 ! Chemoautotrophic bacteria Temperature correction for Mortality rate
      bname(13 ) =            'KHS_NH4N_CHEM_AUT_BAC'  !13 ! Chemoautotrophic bacteria Half saturation growth for NH4N
      bname(14 ) =            'KHS_PO4P_CHEM_AUT_BAC'  !14 ! Chemoautotrophic bacteria Half saturation growth for PO4P
      bname(15 ) =              'KHS_O2_CHEM_AUT_BAC'  !15 ! Chemoautotrophic bacteria Half saturation growth for O2
      bname(16 ) =      'DO_STR_HYPOX_CHEM_AUT_BAC_D'  !16 ! Chemoautotrophic bacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
      bname(17 ) =       'THETA_HYPOX_CHEM_AUT_BAC_D'  !17 ! Chemoautotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
      bname(18 ) =       'EXPON_HYPOX_CHEM_AUT_BAC_D'  !18 ! Chemoautotrophic bacteria Exponent constant for Dissolved oxygen stress
      bname(19 ) =              'CHEM_AUT_BAC_N_TO_C'  !19 ! Chemoautotrophic bacteria Nitrogen to Carbon ratio
      bname(20 ) =              'CHEM_AUT_BAC_P_TO_C'  !20 ! Chemoautotrophic bacteria Phosphorus to Carbon ratio
      bname(21 ) =             'CHEM_AUT_BAC_O2_TO_C'  !21 ! Chemoautotrophic bacteria Oxygen to Carbon ratio
      bname(22 ) =               'YIELD_CHEM_AUT_BAC'  !22 ! Chemoautotrophic bacteria Yield of Carbon per unit Nitrates nitrogen
      bname(23 ) =                'KG_AER_HET_BAC_20'  !23 ! Aerobic heterotrophic bacteria Growth rate
      bname(24 ) =           'EFF_AER_HET_BAC_GROWTH'  !24 ! Aerobic heterotrophic bacteria growth efficiency
      bname(25 ) =             'THETA_KG_AER_HET_BAC'  !25 ! Aerobic heterotrophic bacteria Temperature correction for growth rate
      bname(26 ) =                'KR_AER_HET_BAC_20'  !26 ! Aerobic heterotrophic bacteria Respiration rate
      bname(27 ) =             'THETA_KR_AER_HET_BAC'  !27 ! Aerobic heterotrophic bacteria Temperature correction for respiration rate
      bname(28 ) =                'KD_AER_HET_BAC_20'  !28 ! Aerobic heterotrophic bacteria Mortality rate
      bname(29 ) =             'THETA_KD_AER_HET_BAC'  !29 ! Aerobic heterotrophic bacteria Temperature correction for Mortality rate
      bname(30 ) =             'KHS_ORGC_AER_HET_BAC'  !30 ! Aerobic heterotrophic bacteria Half saturation growth for OC
      bname(31 ) =              'KHS_ORGN_AER_HET_BAC' !31 ! Aerobic heterotrophic bacteria Half saturation growth for ON
      bname(32 ) =              'KHS_ORGP_AER_HET_BAC' !32 ! Aerobic heterotrophic bacteria Half saturation growth for OP
      bname(33 ) =                'KHS_O2_AER_HET_BAC' !33 ! Aerobic heterotrophic bacteria Half saturation growth for Oxygen
      bname(34 ) =               'KHS_DIN_AER_HET_BAC' !34 ! Aerobic heterotrophic bacteria Half saturation growth for inorganic nitrogen
      bname(35 ) =               'KHS_DIP_AER_HET_BAC' !35 ! Aerobic heterotrophic bacteria Half saturation growth for inorganic phosphorus
      bname(36 ) =              'KHS_PHYT_AER_HET_BAC' !36 ! Aerobic heterotrophic bacteria Half saturation growth for Phytoplankton C (not used as a resource)
      bname(37 ) =              'YIELD_OC_AER_HET_BAC' !37 ! Aerobic heterotrophic bacteria Yield of bacteria carbon per unit of organic carbon
      bname(38 ) =               'OX_ORGN_AER_HET_BAC' !38 ! Aerobic heterotrophic bacteria ON oxidation rate mg N per mg C of bacteria production
      bname(39 ) =                        'KHS_MIN_N'  !39 ! Aerobic heterotrophic bacteria ON mineralisation reverse half saturation for DIN
      bname(40 ) =              'OX_ORGP_AER_HET_BAC'  !40 ! Aerobic heterotrophic bacteria OP mineralisation rate mg P per mg C of bacteria production
      bname(41 ) =                        'KHS_MIN_P'  !41 ! Aerobic heterotrophic bacteria ON mineralisation reverse half saturation for DIP
      bname(42 ) =       'DO_STR_HYPOX_AER_HET_BAC_D'  !42 ! Aerobic heterotrophic bacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
      bname(43 ) =        'THETA_HYPOX_AER_HET_BAC_D'  !43 ! Aerobic heterotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
      bname(44 ) =        'EXPON_HYPOX_AER_HET_BAC_D'  !44 ! Aerobic heterotrophic bacteria Exponent constant for Dissolved oxygen stress
      bname(45 ) =               'AER_HET_BAC_N_TO_C'  !45 ! Aerobic heterotrophic bacteria Nitrogen to Carbon ratio
      bname(46 ) =               'AER_HET_BAC_P_TO_C'  !46 ! Aerobic heterotrophic bacteria Phosphorus to Carbon ratio
      bname(47 ) =              'AER_HET_BAC_O2_TO_C'  !47 ! Aerobic heterotrophic bacteria Oxygen to Carbon ratio for respiration
      bname(48 ) =             'KG_FAC_AN_HET_BAC_20'  !48 ! Facultative anaerobic heterotrophic bacteria Growth rate of
      bname(49 ) =        'EFF_FAC_AN_HET_BAC_GROWTH'  !49 ! not used! Facultative anaerobic heterotrophic bacteria growth efficiency
      bname(50 ) =          'THETA_KG_FAC_AN_HET_BAC'  !50 ! not used! Facultative anaerobic heterotrophic bacteria Temperature correction for growth rate
      bname(51 ) =             'KR_FAC_AN_HET_BAC_20'  !51 ! not used! Facultative anaerobic heterotrophic bacteria Respiration rate
      bname(52 ) =          'THETA_KR_FAC_AN_HET_BAC'  !52 ! not used! Facultative anaerobic heterotrophic bacteria Temperature correction for respiration rate
      bname(53 ) =             'KD_FAC_AN_HET_BAC_20'  !53 ! not used! Facultative anaerobic heterotrophic bacteria Mortality rate
      bname(54 ) =          'THETA_KD_FAC_AN_HET_BAC'  !54 ! not used! Facultative anaerobic heterotrophic bacteria Temperature correction for Mortality rate
      bname(55 ) =          'KHS_NO3N_FAC_AN_HET_BAC'  !55 ! Facultative anaerobic heterotrophic bacteria Half saturation growth for NO3N
      bname(56 ) =          'KHS_ORGC_FAC_AN_HET_BAC'  !56 ! not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for OC
      bname(57 ) =          'KHS_ORGN_FAC_AN_HET_BAC'  !57 ! not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for ON
      bname(58 ) =          'KHS_ORGP_FAC_AN_HET_BAC'  !58 ! not used! Facultative anaerobic heterotrophic bacteria Half saturation growth for OP
      bname(59 ) =        'REV_KHS_O2_FAC_AN_HET_BAC'  !59 ! not used! Facultative anaerobic heterotrophic bacteria Reverse Half saturation growth for O2
      bname(60 ) =   'NO3N_LACK_STR_FAC_AN_HET_BAC_D'  !60 ! not used! Facultative anaerobic heterotrophic bacteria NO3N stress concentration
      bname(61 ) =  'THETA_NO3_LACK_FAC_AN_HET_BAC_D'  !61 ! not used! Facultative anaerobic heterotrophic bacteria Multiplier of the exponent for Dissolved oxygen stress
      bname(62 ) =    'EXP_NO3_LACK_FAC_AN_HET_BAC_D'  !62 ! not used! Facultative anaerobic heterotrophic bacteria Exponent constant for Dissolved oxygen stress
      bname(63 ) =            'FAC_AN_HET_BAC_N_TO_C'  !63 ! not used! Facultative anaerobic heterotrophic bacteria Nitrogen to Carbon ratio
      bname(64 ) =            'FAC_AN_HET_BAC_P_TO_C'  !64 ! not used! Facultative anaerobic heterotrophic bacteria Phosphorus to Carbon ratio
      bname(65 ) =           'FAC_AN_HET_BAC_O2_TO_C'  !65 ! not used! Facultative anaerobic heterotrophic bacteria Oxygen to Carbon ratio for respiration
      bname(66 ) =             'YIELD_FAC_AN_HET_BAC'  !66 ! Facultative anaerobic heterotrophic bacteria Yield of carbon per unit nitrates nitrogen
      bname(67 ) =                  'KG_DIA_OPT_TEMP'  !67 ! Diatoms Growth rate
      bname(68 ) =                  'DIA_OPT_TEMP_LR'  !68 ! Diatoms optimal temperature lower range
      bname(69 ) =                  'DIA_OPT_TEMP_UR'  !69 ! Diatoms optimal temperature upper range
      bname(70 ) =                   'EFF_DIA_GROWTH'  !70 ! Diatoms Effective growth. (1-EG)*growth - losses for respiration and excretion
      bname(71 ) =         'KAPPA_DIA_UNDER_OPT_TEMP'  !71 ! Diatoms Temperature correction for growth lower temperature
      bname(72 ) =          'KAPPA_DIA_OVER_OPT_TEMP'  !72 ! Diatoms Temperature correction for growth upper temperature
      bname(73 ) =                        'KR_DIA_20'  !73 ! Diatoms Respiration rate
      bname(74 ) =                     'THETA_KR_DIA'  !74 ! Diatoms Temperature correction for basal respiration rate
      bname(75 ) =                        'KD_DIA_20'  !75 ! Diatoms Mortality rate
      bname(76 ) =                     'THETA_KD_DIA'  !76 ! Diatoms Temperature correction for Mortality rate
      bname(77 ) =                      'KHS_DIN_DIA'  !77 ! Diatoms Half saturation growth for DIN
      bname(78 ) =                      'KHS_DIP_DIA'  !78 ! Diatoms Half saturation growth for DIP
      bname(79 ) =                      'KHS_DSi_DIA'  !79 ! Diatoms Half saturation growth for DSi
      bname(80 ) =                       'KHS_O2_DIA'  !80 ! Diatoms Half saturation growth for O2
      bname(81 ) =                'FRAC_DIA_EXCR'      !81 ! Diatoms Fraction of excretion in metabolism rate
      bname(82 ) =                          'I_S_DIA'  !82 ! Diatoms Light saturation (langleys)
      bname(83 ) =               'DO_STR_HYPOX_DIA_D'  !83 ! Diatoms Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
      bname(84 ) =                'THETA_HYPOX_DIA_D'  !84 ! Diatoms Multiplier of the exponent for Dissolved oxygen stress
      bname(85 ) =                'EXPON_HYPOX_DIA_D'  !85 ! Diatoms Exponent constant for Dissolved oxygen stress
      bname(86 ) =                       'DIA_N_TO_C'  !86 ! Diatoms Nitrogen to Carbon ratio
      bname(87 ) =                       'DIA_P_TO_C'  !87 ! Diatoms Phosphorus to Carbon ratio
      bname(88 ) =                      'DIA_Si_TO_C'  !88 ! Diatoms Silica to Carbon ratio
      bname(89 ) =                      'DIA_O2_TO_C'  !89 ! Diatoms Oxygen to Carbon ratio for respiration
      bname(90 ) =                    'DIA_C_TO_CHLA'  !90 ! Diatoms Carbon to Chlorophil a ratio
      bname(91 ) =                  'KG_CYN_OPT_TEMP'  !91 ! Non-fixing cyanobacteria Growth rate
      bname(92 ) =                  'CYN_OPT_TEMP_LR'  !92 ! Non-fixing cyanobacteria optimal temperature lower range
      bname(93 ) =                  'CYN_OPT_TEMP_UR'  !93 ! Non-fixing cyanobacteria optimal temperature upper range
      bname(94 ) =                   'EFF_CYN_GROWTH'  !94 ! Non-fixing cyanobacteria Effective growth. (1-EG)*growth - losses for respiration and excretion
      bname(95 ) =         'KAPPA_CYN_UNDER_OPT_TEMP'  !95 ! Non-fixing cyanobacteria Temperature correction for growth lower temperature
      bname(96 ) =          'KAPPA_CYN_OVER_OPT_TEMP'  !96 ! Non-fixing cyanobacteria Temperature correction for growth upper temperature
      bname(97 ) =                        'KR_CYN_20'  !97 ! Non-fixing cyanobacteria Respiration rate
      bname(98 ) =                     'THETA_KR_CYN'  !98 ! Non-fixing cyanobacteria Temperature correction for respiration rate
      bname(99 ) =                        'KD_CYN_20'  !99 ! Non-fixing cyanobacteria Mortality rate
      bname(100) =                     'THETA_KD_CYN'  !100! Non-fixing cyanobacteria Temperature correction for Mortality rate
      bname(101) =                      'KHS_DIN_CYN'  !101! Non-fixing cyanobacteria Half saturation growth for DIN
      bname(102) =                      'KHS_DIP_CYN'  !102! Non-fixing cyanobacteria Half saturation growth for DIP
      bname(103) =                       'KHS_O2_CYN'  !103! Non-fixing cyanobacteria Half saturation growth for O2
      bname(104) =                'FRAC_CYN_EXCR'      !104!  Non-fixing cyanobacteria Fraction of excretion in metabolism rate
      bname(105) =                          'I_S_CYN'  !105! Non-fixing cyanobacteria Light saturation (langleys)
      bname(106) =               'DO_STR_HYPOX_CYN_D'  !106! Non-fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
      bname(107) =                'THETA_HYPOX_CYN_D'  !107! Non-fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
      bname(108) =                'EXPON_HYPOX_CYN_D'  !108! Non-fixing cyanobacteria Exponent constant for Dissolved oxygen stress
      bname(109) =                       'CYN_N_TO_C'  !109! Non-fixing cyanobacteria Nitrogen to Carbon ratio ,was 0.1
      bname(110) =                       'CYN_P_TO_C'  !110! Non-fixing cyanobacteria Phosphorus to Carbon ratio
      bname(111) =                      'CYN_O2_TO_C'  !111! Non-fixing cyanobacteria Oxygen to Carbon ratio for respiration
      bname(112) =                    'CYN_C_TO_CHLA'  !112! Non-fixing cyanobacteria Carbon to Chlorophyl a ratio
      bname(113) =              'KG_FIX_CYN_OPT_TEMP'  !113! Fixing cyanobacteria Growth rate
      bname(114) =              'FIX_CYN_OPT_TEMP_LR'  !114! Fixing Cyanobacteria optimal temperature lower range
      bname(115) =              'FIX_CYN_OPT_TEMP_UR'  !115! Fixing Cyanobacteria optimal temperature upper range
      bname(116) =               'EFF_FIX_CYN_GROWTH'  !116! Fixing cyanobacteria Effective growth. (1-EG)*growth - losses for RESP and excretion
      bname(117) =     'KAPPA_FIX_CYN_UNDER_OPT_TEMP'  !117! Fixing cyanobacteria Temperature correction for growth lower temperature
      bname(118) =      'KAPPA_FIX_CYN_OVER_OPT_TEMP'  !118! Fixing cyanobacteria Temperature correction for growth upper temperature
      bname(119) =                    'KR_FIX_CYN_20'  !119! Fixing cyanobacteria RESP rate
      bname(120) =                 'THETA_KR_FIX_CYN'  !120! Fixing cyanobacteria Temperature correction for RESP rate
      bname(121) =                    'KD_FIX_CYN_20'  !121! Fixing cyanobacteria Mortality rate of nitrification bacteria
      bname(122) =                 'THETA_KD_FIX_CYN'  !122! Fixing cyanobacteria Temperature correction for Mortality rate
      bname(123) =                  'KHS_DIN_FIX_CYN'  !123! Fixing cyanobacteria Half saturation growth for DIN
      bname(124) =                  'KHS_DIP_FIX_CYN'  !124! Fixing cyanobacteria Half saturation growth for DIP
      bname(125) =                   'KHS_O2_FIX_CYN'  !125! Fixing cyanobacteria Half saturation growth for O2
      bname(126) =            'FRAC_FIX_CYN_EXCR'      !126! Fixing cyanobacteria Fraction of excretion in metabolism rate
      bname(127) =                      'I_S_FIX_CYN'  !127! Fixing cyanobacteria Light saturation (langleys)
      bname(128) =           'DO_STR_HYPOX_FIX_CYN_D'  !128! Fixing cyanobacteria Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
      bname(129) =            'THETA_HYPOX_FIX_CYN_D'  !129! Fixing cyanobacteria Multiplier of the exponent for Dissolved oxygen stress
      bname(130) =            'EXPON_HYPOX_FIX_CYN_D'  !130! Fixing cyanobacteria Exponent constant for Dissolved oxygen stress
      bname(131) =                   'FIX_CYN_N_TO_C'  !131! Fixing cyanobacteria Nitrogen to Carbon ratio
      bname(132) =                   'FIX_CYN_P_TO_C'  !132! Fixing cyanobacteria Phosphorus to Carbon ratio
      bname(133) =                  'FIX_CYN_O2_TO_C'  !133! Fixing cyanobacteria Oxygen to Carbon ratio for respiration
      bname(134) =                'FIX_CYN_C_TO_CHLA'  !134! Fixing cyanobacteria Carbon to Chlorophyl a ratio
      bname(135) =                            'R_FIX'  !135! Fixing cyanobacteria Ratio between non-fixing and fixing fractions growth rate
      bname(136) =                            'K_FIX'  !136! Fixing cyanobacteria Effectivity parameter of switching to nitrogen fixation
      bname(137) =                  'KG_OPA_OPT_TEMP'  !137! OtherPhyto Growth rate
      bname(138) =                  'OPA_OPT_TEMP_LR'  !138! OtherPhyto optimal temperature lower range
      bname(139) =                  'OPA_OPT_TEMP_UR'  !139! OtherPhyto optimal temperature upper range
      bname(140) =                   'EFF_OPA_GROWTH'  !140! OtherPhyto Effective growth. (1-EG)*growth - losses for respiration and excretion
      bname(141) =         'KAPPA_OPA_UNDER_OPT_TEMP'  !141! OtherPhyto Temperature correction for growth lower temperature
      bname(142) =          'KAPPA_OPA_OVER_OPT_TEMP'  !142! OtherPhyto Temperature correction for growth upper temperature
      bname(143) =                        'KR_OPA_20'  !143! OtherPhyto Respiration rate
      bname(144) =                     'THETA_KR_OPA'  !144! OtherPhyto Temperature correction for respiration rate
      bname(145) =                        'KD_OPA_20'  !145! OtherPhyto Mortality rate
      bname(146) =                     'THETA_KD_OPA'  !146! OtherPhyto Temperature correction for Mortality rate
      bname(147) =                      'KHS_DIN_OPA'  !147! OtherPhyto Half saturation growth for DIN
      bname(148) =                      'KHS_DIP_OPA'  !148! OtherPhyto Half saturation growth for DIP
      bname(149) =                       'KHS_O2_OPA'  !149! OtherPhyto Half saturation growth for O2
      bname(150) =                'FRAC_OPA_EXCR'      !150! OtherPhyto Fraction of excretion in metabolism rate
      bname(151) =                          'I_S_OPA'  !151! OtherPhyto Light saturation (langleys)
      bname(152) =               'DO_STR_HYPOX_OPA_D'  !152! OtherPhyto Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
      bname(153) =                'THETA_HYPOX_OPA_D'  !153! OtherPhyto Multiplier of the exponent for Dissolved oxygen stress
      bname(154) =                'EXPON_HYPOX_OPA_D'  !154! OtherPhyto Exponent constant for Dissolved oxygen stress
      bname(155) =                       'OPA_N_TO_C'  !155! OtherPhyto Nitrogen to Carbon ratio
      bname(156) =                       'OPA_P_TO_C'  !156! OtherPhyto Phosphorus to Carbon ratio
      bname(157) =                      'OPA_O2_TO_C'  !157! OtherPhyto Oxygen to Carbon ratio for respiration
      bname(158) =                    'OPA_C_TO_CHLA'  !158! OtherPhyto Carbon to Chlorophyl a ratio
      bname(159) =                  'KG_ZOO_OPT_TEMP'  !159! Zooplankton Growth rate
      bname(160) =                  'ZOO_OPT_TEMP_LR'  !160! Zooplankton optimal temperature lower range
      bname(161) =                  'ZOO_OPT_TEMP_UR'  !161! Zooplankton optimal temperature upper range
      bname(162) =                   'EFF_ZOO_GROWTH'  !162! Zooplankton Effective growth. (1-EG)*growth - losses for respiration and excretion
      bname(163) =         'KAPPA_ZOO_UNDER_OPT_TEMP'  !163! Zooplankton Temperature correction for growth lower temperature
      bname(164) =          'KAPPA_ZOO_OVER_OPT_TEMP'  !164! Zooplankton Temperature correction for growth upper temperature
      bname(165) =                     'GRAT_ZOO_DIA'  !165! Zooplankton Grazing rate (growhth rate multiplier) on diatoms
      bname(166) =                     'GRAT_ZOO_CYN'  !166! Zooplankton Grazing rate (growhth rate multiplier) on Cyanobacteria
      bname(167) =                     'GRAT_ZOO_OPA'  !167! Zooplankton Grazing rate (growhth rate multiplier) on fixing Cyanobacteria
      bname(168) =                 'GRAT_ZOO_FIX_CYN'  !168! Zooplankton Grazing rate (growhth rate multiplier) on OtherPhyto
      bname(169) =            'GRAT_ZOO_CHEM_AUT_BAC'  !169! Zooplankton Grazing rate (growhth rate multiplier) on NITR_BAC
      bname(170) =             'GRAT_ZOO_AER_HET_BAC'  !170! Zooplankton Grazing rate (growhth rate multiplier) on AER_HET_BAC
      bname(171) =          'GRAT_ZOO_FAC_AN_HET_BAC'  !171! Zooplankton Grazing rate (growhth rate multiplier) on DENITR_BAC
      bname(172) =          'GRAT_ZOO_DET_PART_ORG_C'  !172! Zooplankton Grazing rate (growhth rate multiplier) on part. ORG_C
      bname(173) =                     'PREF_ZOO_DIA'  !173! Zooplankton Preference for diatoms
      bname(174) =                     'PREF_ZOO_CYN'  !174! Zooplankton Preference for Cyanobacteria
      bname(175) =                 'PREF_ZOO_FIX_CYN'  !175! Zooplankton Preference for fixing Cyanobacteria
      bname(176) =                     'PREF_ZOO_OPA'  !176! Zooplankton Preference for OtherPhyto
      bname(177) =            'PREF_ZOO_CHEM_AUT_BAC'  !177! Zooplankton Preference for NITR_BAC
      bname(178) =             'PREF_ZOO_AER_HET_BAC'  !178! Zooplankton Preference for AER_HET_BAC
      bname(179) =          'PREF_ZOO_FAC_AN_HET_BAC'  !179! Zooplankton Preference for DENITR_BAC
      bname(180) =          'PREF_ZOO_DET_PART_ORG_C'  !180! Zooplankton Preference for part. ORG_C
      bname(181) =                    'KHS_DIA_C_ZOO'  !181! Zooplankton Half saturation growth for diatoms
      bname(182) =                    'KHS_CYN_C_ZOO'  !182! Zooplankton Half saturation growth for Cyanobacteria
      bname(183) =                'KHS_FIX_CYN_C_ZOO'  !183! Zooplankton Half saturation growth for fixing Cyanobacteria
      bname(184) =                    'KHS_OPA_C_ZOO'  !184! Zooplankton Half saturation growth for OtherPhyto
      bname(185) =           'KHS_CHEM_AUT_BAC_C_ZOO'  !185! Zooplankton Half saturation growth for NITR_BAC
      bname(186) =            'KHS_AER_HET_BAC_C_ZOO'  !186! Zooplankton Half saturation growth for AER_HET_BAC
      bname(187) =         'KHS_FAC_AN_HET_BAC_C_ZOO'  !187! Zooplankton Half saturation growth for DENITR_BAC
      bname(188) =           'KHS_DET_PART_ORG_C_ZOO'  !188! Zooplankton Half saturation growth for part. ORG_C
      bname(189) =                     'FOOD_MIN_ZOO'  !189! Zooplankton Minimum food conc. for feeding
      bname(190) =                           'KE_ZOO'  !190! not used Zooplankton Excretion rate as growth fraction
      bname(191) =                  'FRAC_ZOO_EX_ORG'  !191! not used Zooplankton Excretion rate organic fraction
      bname(192) =                        'KR_ZOO_20'  !192! Zooplankton Respiration rate
      bname(193) =                     'THETA_KR_ZOO'  !193! Zooplankton Respiration rate Temperature correction
      bname(194) =                        'KD_ZOO_20'  !194! Zooplankton Mortality rate
      bname(195) =                     'THETA_KD_ZOO'  !195! Zooplankton Mortality rate Temperature correction
      bname(196) =               'DO_STR_HYPOX_ZOO_D'  !196! Zooplankton Dissolved oxygen stress in oxygen units (mortality increase below this value exponentialy
      bname(197) =                'THETA_HYPOX_ZOO_D'  !197! Zooplankton Multiplier of the exponent for Dissolved oxygen stress
      bname(198) =                'EXPON_HYPOX_ZOO_D'  !198! Zooplankton Exponent constant for Dissolved oxygen stress
      bname(199) =                       'ZOO_N_TO_C'  !199! Zooplankton Nitrogen to Carbon ratio
      bname(200) =                       'ZOO_P_TO_C'  !200! Zooplankton Phosphorus to Carbon ratio
      bname(201) =                      'ZOO_O2_TO_C'  !201! Zooplankton Oxygen to Carbon ratio for respiration
      bname(202) =          'KDISS_DET_PART_ORG_C_20'  !202! Particulate Detritus Carbon Dissolution rate not dependent on phytoplankton
      bname(203) =       'THETA_KDISS_DET_PART_ORG_C'  !203! Particulate Detritus Carbon Dissolution rate Temperature correction
      bname(204) =          'FAC_PHYT_DET_PART_ORG_C'  !204! Particulate Detritus Carbon Phytoplankton linear factor for dissolution rate
      bname(205) =          'KDISS_DET_PART_ORG_N_20'  !205! Particulate Detritus Nitrogen Dissolution rate not dependent on phytoplankton
      bname(206) =       'THETA_KDISS_DET_PART_ORG_N'  !206! Particulate Detritus Nitrogen Dissolution rate Temperature correction
      bname(207) =                       'KHS_DISS_N'  !207! Particulate Detritus Nitrogen dissolution reverse half saturation for DIN
      bname(208) =          'FAC_PHYT_DET_PART_ORG_N'  !208! Particulate Detritus Nitrogen Phytoplankton linear factor for dissolution rate
      bname(209) =          'KDISS_DET_PART_ORG_P_20'  !209! Particulate Detritus Phosphorus Dissolution rate not dependent on phytoplankton
      bname(210) =       'THETA_KDISS_DET_PART_ORG_P'  !210! Particulate Detritus Phosphorus Dissolution rate Temperature correction
      bname(211) =                       'KHS_DISS_P'  !211! Particulate Detritus Phosphorus  dissolution reverse half saturation for DIP
      bname(212) =          'FAC_PHYT_DET_PART_ORG_P'  !212! Particulate Detritus Phosphorus  Phytoplankton linear factor for dissolution rate
      bname(213) =                 'KDISS_PART_Si_20'  !213! Particulate Silica Dissolution rate
      bname(214) =              'THETA_KDISS_PART_Si'  !214! Particulate Silica Dissolution rate Temperature correction
      bname(215) =                     'K_MIN_DOC_20'  !215! Dissolved carbon  mineralisation rate
      bname(216) =                  'THETA_K_MIN_DOC'  !216! Dissolved carbon  mineralisation rate Temperature constant
      bname(217) =                'FAC_PHYT_AMIN_DOC'  !217! Dissolved carbon  Phytoplankton linear factor for mineralisation rate
      bname(218) =                     'K_MIN_DON_20'  !218! Dissolved nitrogen  mineralisation rate not dependent on phytoplankton
      bname(219) =                  'THETA_K_MIN_DON'  !219! Dissolved nitrogen  mineralisation rate Temperature constant
      bname(220) =                       'KHS_AMIN_N'  !220! Dissolved nitrogen  reverse half saturation for DIN
      bname(221) =                'FAC_PHYT_AMIN_DON'  !221! Dissolved nitrogen Phytoplankton linear factor for mineralisation rate
      bname(222) =                     'K_MIN_DOP_20'  !222! Dissolved phosphorus  mineralisation rate not dependent on phytoplankton
      bname(223) =                  'THETA_K_MIN_DOP'  !223! Dissolved phosphorus  mineralisation rate Temperature constant
      bname(224) =                       'KHS_AMIN_P'  !224! Dissolved phosphorus reverse half saturation for DIP
      bname(225) =                'FAC_PHYT_AMIN_DOP'  !225! Dissolved phosphorus Phytoplankton linear factor for mineralisation rate
      bname(226) =                        'K_NITR_20'  !226! Amonia nitrification rate
      bname(227) =                     'KHS_NITR_OXY'  !227! Amonia nitrification half saturation for Oxygen
      bname(228) =                   'KHS_NITR_NH4_N'  !228! Amonia nitrification half saturation for Amonia
      bname(229) =                     'THETA_K_NITR'  !229! Amonia nitrification rate Temperature constant
                                                           
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
      parameter (npar = 42)
      character*80 bname(npar)
      real par(npar)

!     execution part

      if( npar. ne. nconst) then
        print *, 'CUR_PARAM_READ_INI_WC:'
        print *, '  wrong asked number of params'
        print *, '  asked: ',nconst
        print *, '  existing: ',npar
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
!   par     - parameters array. Is initialized but no used more. Fixme
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
      double precision getpar
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
        call trimline(file,ifile)
!        write(6,*)'Constants initialized from file: ',file(1:ifile)
        if( iunit .le. 0 ) then
         write(6,'(a,a)') 'filename: ',file(1:ifile)
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
                call triml(bname(j))
                if (name .eq. bname(j)) then
                  
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
               print *, 'CUR_PARAM_READ: Parameter did no get value'
               print *, 'Parameter name', trim(name) 
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
          if (bname(j). ne. ' ' ) write(*,44) j,bname(j),par(j)
        end do

        write(*,*)
44      format(3x,i3,2x,a40,f12.7)

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

	implicit none

	include 'aquabc_II.h'
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
