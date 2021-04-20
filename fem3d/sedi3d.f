
!--------------------------------------------------------------------------
!
!    Copyright (C) 2005-2006,2008-2011,2015-2019  Christian Ferrarin
!    Copyright (C) 2008,2010,2012-2020  Georg Umgiesser
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

! ------ SUBROUTINE SEDI3D ---------
! ***********************************
!
! This routine manages the sediment transport computation.
! It calls the subrourine SEDTRANS05 which computes
! the sediment transport rate and the suspended sediment
! concentration in case of current or current and waves
! (see SUBWAVES routine). The routine scn_compute handles the
! avdection and diffusion of suspended sediment.
!
!  16/11/2004 Christian Ferrarin ISMAR-CNR
! 
! revision log :
!
! 01.03.2005	ccf	(sedi3d_f1.f) coming from sedi3d_e3_multi7.f
! ...				new cohesive routine
! 01.03.2005	ccf	(sedi3d_f2.f) add mixing thickness
! 01.03.2005	ccf	(sedi3d_f3.f) merge layer 1 and 2 if bedn(1) < bmix
! 01.03.2005	ccf	(sedi3d_f4.f) update element depth
! 01.03.2005	ccf	(sedi3d_f5.f) get viscosity from new routine
! 01.03.2005	ccf	(sedi3d_f6.f) active layer = bottom roughness heigh
! 01.04.2005	ccf	(sedi3d_f7.f) style changes, create sedt05
! 01.04.2005	ccf	(sedi3d_f8.f) initialize percbd from file
! 01.04.2005	ccf	(sedi3d_f8.f) change rhos per cohesive sed
! 01.06.2005	ccf	(sedi3d_f13.f) adapt 3D, bottom layer and total depth
! ...				bugfix in computing bedn(1,3), dimension in tuek
! 01.06.2005	ccf	(sedi3d_f15.f) change deposition density,
! ...				add consolidation
! 01.06.2005	ccf	(sedi3d_f16.f) change units, bug fix in TCONC1
! 01.06.2005	ccf	(sedi3d_f17.f) bug fix in upedepth
! 01.06.2005	ccf	(sedi3d_f18.f) bug fix in checkbed, 
! ...				adapt to sedtrans05-h5.f
! 01.07.2005	ccf	(sedi3d_f20.f) adapt to sedtrans05-h6.f
! 01.07.2005	ccf	(sedi3d_f21.f) bug fix in updepth (as subdry.f)
! 01.07.2005	ccf	(sedi3d_f22.f) bug fix in bedload
! 01.08.2005	ccf	(sedi3d_f23.f) eliminate ripple variables
! 01.08.2005	ccf	(sedi3d_f24.f) bed slope contribution. adjust updepth
! 01.08.2005	ccf	(sedi3d_f25.f) change deposition characteristics
! 01.09.2005	ccf	(sedi3d_f26.f) bug fix in bedman,
! ...				change boundary for cohesive
! 01.09.2005	ccf	(sedi3d_f27.f) separete erosion/deposition in bedman
! ...				number of layer computed each time
! 01.11.2005	ccf	(sedi3d_f28f) adapt for 3D version
! 01.11.2005	ccf	(sedi3d_f29f) bed slope threshold, bug fix in updepth
! 01.11.2005	ccf	(sedi3d_f30f) bug fix in bedslope and in getmd
! 01.11.2005	ccf	(sedi3d_f31f) bed slope by gradients
! 01.11.2005	ccf	(sedi3d_f32f) last layer not smaller than 0.1 m
! 01.11.2005	ccf	(sedi3d_f33f) compute vertical mixing coeffcients
! 01.01.2006	ccf	(sedi3d_f34f) bug fix in upedepth and blimit
! 01.02.2006	ccf	(sedi3d_f35f) new suspco routine
! 01.02.2006	ccf	(sedi3d_f36f) adapt Van Rijn (C0) in sedtrans05-h8.f
! 01.02.2006	ccf	(sedi3d_f37f) bug fix in checkbed and other things
! 01.05.2006	ccf	(sedi3d_f38f) bug fix in line 1119. Introduce KCOES
! 01.05.2006	ccf	(sedi3d_f39f) bugs nonco. non used edr in cohse.
! ...				no limit percbd. pers(1)>0.1 in line 1120
! 01.05.2006	ccf	(sedi3d_f40f) limit BEDCHA(1,2) to 0.1. 
! ...				bugfix in suspco, better conc in cohse
! 01.06.2006	ccf	(sedi3d_f41f) no transport in case of depth< 0.1m
! 01.06.2006	ccf	(sedi3d_f42f) read constants from file
! 01.06.2006	ccf	(sedi3d_f43f) add limcoh
! 01.06.2006	ccf	(sedi3d_f44f) smooth bed elevation, get_timestep
! 01.07.2006	ccf	(sedi3d_f45f) read angle of repose, 
! ...				limit shear velocity, write bathymetry
! 11.04.2008	ggu&ccf	treatment of boundaries slightly changed
! 16.04.2008	ggu&ccf	bugfix calling lin (must pass double precision 0.)
! 22.04.2008	ggu	advection parallelized
! 27.04.2008	ccf	new check, other changes
! 28.04.2008	ggu	call to nrdpar, nrdnls read with double precision
! 29.04.2008	ccf	bug fix in depbed
! 20.10.2008	ccf	add routine for computing the active layer
! 30.04.2009	ccf	add adjtime for initialization to reset bh
! 03.10.2009	ccf	compute transport only for GD
! 23.03.2010	ggu	changed v6.1.1
! 13.04.2010	ccf	introduce SSC effect on density
! 27.07.2010	ccf	settling velocity for cohesives as function of space
! 08.10.2010	ggu	changed VERS_6_1_13
! 20.06.2011	ccf	load for eros/depos for both cohesive and non-cohesive
! 20.06.2011	ccf	deleted suspco routine
! 30.03.2012	ggu	changed VERS_6_1_51
! 21.06.2012	ggu	changed VERS_6_1_54
! 25.01.2013	ggu	changed VERS_6_1_62
! 13.06.2013	ggu	changed VERS_6_1_65
! 18.06.2014	ggu	changed VERS_6_1_77
! 21.10.2014	ggu	changed VERS_7_0_3
! 05.11.2014	ggu	changed VERS_7_0_5
! 05.12.2014	ggu	changed VERS_7_0_8
! 23.12.2014	ggu	changed VERS_7_0_11
! 09.01.2015	ggu	changed VERS_7_0_12
! 19.01.2015	ccf	ia_out1/2 introduced
! 19.01.2015	ggu	changed VERS_7_1_3
! 10.02.2015	ggu	new read for constants
! 26.02.2015	ggu	changed VERS_7_1_5
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 30.07.2015	ggu	changed VERS_7_1_83
! 23.09.2015	ggu	changed VERS_7_2_4
! 12.10.2015	ggu	changed VERS_7_3_3
! 22.10.2015	ggu	changed VERS_7_3_8
! 05.11.2015	ggu	changed VERS_7_3_12
! 16.12.2015	ggu	changed VERS_7_3_16
! 18.12.2015	ggu	changed VERS_7_3_17
! 11.03.2016	ggu	changed VERS_7_5_5
! 31.03.2016	ggu&ccf	update to new shyfem version
! 01.04.2016	ggu	changed VERS_7_5_7
! 19.04.2016	ccf	krocks for no erosion-deposition in specific e-types
! 25.05.2016	ggu	changed VERS_7_5_10
! 30.09.2016	ggu	changed VERS_7_5_18
! 10.02.2017	ggu	read in init data from fem files (init_file_fem)
! 13.04.2017	ggu	changed VERS_7_5_25
! 29.08.2017	ccf	read parameters as array
! 02.09.2017	ggu	changed VERS_7_5_31
! 26.09.2017	ggu	changed VERS_7_5_32
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	changed VERS_7_5_36
! 17.11.2017	ggu	readsed split in readsed and initsed
! 05.12.2017	ggu	changed VERS_7_5_39
! 22.02.2018	ccf	adjusted for new sinking velocity
! 03.04.2018	ggu	changed VERS_7_5_43
! 03.02.2019	ggu	adjust z0 etc (GGUZ0)
! 05.02.2019	ggu	upgraded to new version from chris
! 14.02.2019	ggu	set negative conz values to 0
! 15.02.2019	ccf	pass PHIB,PHI100 into nonco (bug)
! 13.03.2019	ggu	changed VERS_7_5_61
! 24.10.2019	ggu	better checking of grainsize percentage in inibed
! 16.02.2020	ggu	femtime eliminated
! 22.09.2020    ggu     correct warnings for PGI compiler
! 
!****************************************************************************

      module mod_sediment

      implicit none

      integer, private, save    :: nkn_sedim = 0
      integer, private, save    :: nlv_sedim = 0

      integer, save	        :: isedi   = 0		!sedi call parameter
      integer, parameter        :: nlbdim  = 10         !max number of bed layer

      character*4, save         :: what = 'sedt'

      double precision, save    :: da_sed(4) = 0
      double precision, save    :: da_ssc(4) = 0

      real, allocatable, save	:: tcn(:,:) 	!total sediment concentration [kg/m3]
      real, allocatable, save	:: tmsus(:)	!total suspended sediment load [kg/s] (> 0 -eros, < 0 depo)
      real, allocatable, save	:: bdens(:,:)	!dry bulk density of bottom sediments [kg/m**3]
      real, allocatable, save	:: bleve(:,:)	!depth below sediment surface of sediments [m]

!==================================================================
      contains
!==================================================================

      subroutine mod_sedim_init(nkn,nlv)

      integer  :: nkn
      integer  :: nlv

      if( nlv == nlv_sedim .and. nkn == nkn_sedim ) return

      if( nlv > 0 .or. nkn > 0 ) then
        if( nlv == 0 .or. nkn == 0 ) then
          write(6,*) 'nlv,nkn: ',nlv,nkn
          stop 'error stop mod_sedim_init: incompatible parameters'
        end if
      end if

      if( nkn_sedim > 0 ) then
        deallocate(tcn)
        deallocate(tmsus)
        deallocate(bdens)
        deallocate(bleve)
      end if

      nlv_sedim = nlv
      nkn_sedim = nkn

      if( nkn == 0 ) return

      allocate(tcn(nlv,nkn))
      allocate(tmsus(nkn))
      allocate(bdens(nlbdim,nkn))
      allocate(bleve(nlbdim,nkn))

      tcn = 0.
      tmsus = 0.
      bdens = 0.
      bleve = 0.

      end subroutine mod_sedim_init

! ==================================================================

      end module mod_sediment

! ==================================================================

      module mod_sediment_para

      implicit none

      double precision, save, dimension(8)  :: sedpa    !sediment parameter vector
      real, allocatable, save   :: gsc(:)		!grainsize class read from str
      real, allocatable, save   :: tue(:)		!initial erosion threshold
      real, allocatable, save   :: prin(:,:)		!initial percetage [0,1]

      double precision, save    :: KCOES   = 0.15d0     !Fraction of mud for sediment to be cohesive [0-1]
      double precision, save    :: LIMCOH  = 0.000063d0 !grainsize limit for cohesive sediment [m]
      double precision, save    :: SMOOTH  = 1.0d0      !smoothing factor for morphodynamic [0-1]
      double precision, save    :: ANGREP  = 32.0d0     !angle of repose [rad]
      double precision, save    :: MORPHO  = 1.d0       !Morphological factor
      double precision, save    :: RHOSED  = 2650.d0    !sediment grain density
      double precision, save    :: POROS   = 0.4d0      !bed porosity [0,1]
      double precision, save    :: SURFPOR = 0.6d0      !surficial porosity [0,1]
      integer, save             :: IOPT    = 5          !SEDIMENT TRANSPORT FORMULA OPTION NUMBER
      integer, save             :: NBCC                 !same as NBCONC but equal 0 if no cohesive

      real, save		:: difmol	!Molecolar diffusion coefficient [m**2/s]
      logical, save		:: wsetbio = .true.	!temp dependent settling velocity
      !logical, save		:: wsetbio = .false.	!temp dependent settling velocity

      logical, save :: bbudget = .false.		!write budget file

! ==================================================================

      end module mod_sediment_para

! ==================================================================

      subroutine sedi

      use mod_depth
      use mod_diff_visc_fric
      use levels
      use basin, only : nkn
      use mod_sediment
      use mod_sediment_para
      use mod_sedtrans05
      use mod_bstress, only : taubot

      implicit none

! -------------------------------------------------------------
! fem variables
! -------------------------------------------------------------
      !integer		:: it		!actual time in seconds
      real		:: dt		!time step in seconds
      real, save	:: salref	!reference salinity 
      real, save	:: temref	!reference temperature [C]
      double precision, save	:: hzoff

! -------------------------------------------------------------
! local sediment variables
! -------------------------------------------------------------
      integer, save	:: icall = 0
      integer, save	:: nscls	!number of grainsize classes
      integer, save	:: nbed		!initial number of bed layers
      integer, save     :: ius = 0	!unit for budget file
      double precision, save	:: adjtime	!time for initialization [s]
      integer, allocatable, save :: idsedi(:)   !information on boundaries
      real, allocatable, save :: gskm(:) 	!AVERAGE SEDIMENT GRAIN DIAMETER ON NODES (M)
      real, allocatable, save :: sbound(:)      !boundary vector [kg/m3]
      real, allocatable, save :: gdx(:) 
      real, allocatable, save :: gdy(:) 	!slope gradients
      real, allocatable, save :: tao(:)         !wave-current shear stress
      real, allocatable, save :: totbed(:) 	!total bedload transport (kg/ms)
      real, allocatable, save :: v1v(:) 
      real, allocatable, save :: percc(:)
      real, allocatable, save :: percs(:)
      integer, allocatable, save :: krocks(:)   !node index with rocks
      double precision, allocatable, save :: bh(:)    	!bottom height variation [m]
      double precision, allocatable, save :: eps(:,:,:)	!vertical mixing coefficient
      double precision, allocatable, save :: scn(:,:,:) !suspended sediment conc (kg/m3)
      double precision, allocatable, save :: sload(:,:) !suspended sediment load [kg/s]
      double precision, allocatable, save :: wsink(:,:,:) !settling velocity for suspended sediment
      double precision, allocatable, save :: bdh(:)     !total elevation change per timestep [>0depo,<0ero]
      double precision, allocatable, save :: bedk(:)    !kg of sediment eroded or deposited [>0depo,<0ero]
      double precision, allocatable, save :: gs(:) 	!SEDIMENT GRAIN DIAMETER (M)
      double precision, allocatable, save :: riph(:)    !ripple height [m]
      double precision, allocatable, save :: ripl(:)    !ripple length [m]
      double precision, allocatable, save :: sflx(:,:) 	!flux of suspend sediment [m3/m2]
      double precision, allocatable, save :: sedx(:,:) 	!x bedload component [kg/ms]
      double precision, allocatable, save :: sedy(:,:)  !y bedload component [kg/ms]
      double precision, allocatable, save :: bflx(:,:) 	!flux of bedload sediment [m3/m2]
      double precision, allocatable, save :: bload(:)   !flux of bedload sediment [kg]
      double precision, allocatable, save :: percbd(:,:,:) !fraction of sediment [0,1]
      double precision, allocatable, save :: bedn(:,:,:)!bed characteristics in 3 columns
                                                	! (1) depth below sediment surface (m)
                                                	! (2) critical erosion stress (Pa)
                                                	! (3) dry bulk density (kg/m**3)

      integer 		:: is,k,l
      integer		:: nbc
      integer		:: nbnds
      integer		:: nintp
      integer		:: nvar
      double precision	:: dtmsed 	!output time parameter
      real 		:: thick	!initial thickness [m]
      real		:: conref	!initial concentration [kg/m3]
      double precision  :: dtime	!time of simulation [s]
      double precision	:: timedr	!time step [s]
	
! function
      integer 		:: ifemop
      real		:: getpar
      logical		:: has_output_d

! ----------------------------------------------------------
! Documentation
! ----------------------------------------------------------
!
! nscls		number of different grain size classes
! nbcc		number of cohesive grain size classes (automatic)
! gs(i)		value for single grain size classes
! sedx/y	-> maybe invert indices
!
! ----------------------------------------------------------
! Initialization
! This section is called only the first time step when ICALL = 0
! ----------------------------------------------------------

        if( icall .le. -1 ) return

        isedi = nint(sedpa(1))
	if( isedi .le. 0 ) then
	  icall = -1
	  return
	end if

        if( icall .eq. 0 ) then

          call get_act_dtime(dtime)
	  call convert_date_d('itmsed',dtmsed)

          if( dtime .lt. dtmsed ) return
          icall = 1

!         --------------------------------------------------
!	  Initialize state variables and constants
!         --------------------------------------------------
          nscls  = nint(sedpa(4))	!number of grain size classes
          conref = sedpa(2)/nscls	!initial concentration
	  call convert_date_d('adjtime',adjtime)
          sedpa(7) = adjtime		!time for re-initialization
          difmol = getpar('difmol')	!molecular vertical diffusivity [m**2/s]
	  hzoff  = getpar('hzon') + 0.10d0
          nbed   = 6			!initial number of bed layer
          thick  = 0.0			!initial bed layer thickness [m]
          temref = getpar('temref')
          salref = getpar('salref')

!         --------------------------------------------------
!	  Read parameter values from file
!         --------------------------------------------------
          call readsedconst

!         --------------------------------------------------
!	  Allocate arrays
!         --------------------------------------------------
          allocate(bh(nkn))
          allocate(gskm(nkn))
          allocate(sbound(nscls))
          allocate(scn(nlvdi,nkn,nscls))
          allocate(eps(0:nlvdi,nkn,nscls))
          allocate(wsink(0:nlvdi,nkn,nscls))
          allocate(gdx(nkn))
          allocate(gdy(nkn))
          allocate(tao(nkn))
          allocate(totbed(nkn))
          allocate(v1v(nkn))
          allocate(sload(nkn,nscls))
          allocate(bdh(nkn))
          allocate(bedk(nkn))
          allocate(gs(nscls))
          allocate(riph(nkn))
          allocate(ripl(nkn))
          allocate(sflx(nscls,nkn))
          allocate(sedx(nscls,nkn))
          allocate(sedy(nscls,nkn))
          allocate(bflx(nscls,nkn))
          allocate(bload(nkn))
          allocate(percbd(nlbdim,nscls,nkn))
          allocate(bedn(nlbdim,3,nkn))
          allocate(percc(nkn))
          allocate(percs(nkn))
	  allocate(krocks(nkn))

	  nbc = nbnds()
	  allocate(idsedi(nbc))

!         --------------------------------------------------
!	  Initialize arrays
!         --------------------------------------------------
          gs  = gsc
	  sbound = 0.

	  nbcc = 0
          do is = 1,nscls
	   if (gs(is) .lt. limcoh) nbcc = nbcc + 1
          end do

          bh     = 0.d0
          bdh    = 0.d0
          bedk   = 0.d0
          gdx    = 0.
          gdy    = 0.
          riph   = 0.d0
          ripl   = 0.d0
          totbed = 0.
          sedx   = 0.d0	
          sedy   = 0.d0
          sload  = 0.d0
          scn    = conref
          eps    = 0.d0
          wsink  = 0.d0
          tcn    = 0.		
          tmsus  = 0.
          bdens  = 0.
          bleve  = 0.
          bload  = 0.d0
	  tao = 0.

!         --------------------------------------------------
!         Initialize SEDTRANS05 arrays WSI
!         --------------------------------------------------
          allocate(WSI(nbcc))
	  WSI = 0.d0

!         --------------------------------------------------
!         Initialize bed configuration
!         --------------------------------------------------
          call inibed(gs,nscls,nbed,thick,gskm,percbd,bedn,krocks)

!         --------------------------------------------------
!         Set boundary conditions for all state variables
!         --------------------------------------------------
          call get_first_dtime(dtime)
          nintp = 2
          nvar = nscls
          call bnds_init_new(what,dtime,nintp,nvar,nkn,nlvdi
     +                          ,sbound,idsedi)

!         --------------------------------------------------
!	  Initialize output
!         --------------------------------------------------

	  call init_sed_output
	  if( bbudget ) then
            ius = ifemop('.sed.bdg','formatted','new')
	    call get_timestep(dt)
            timedr = dt
            call totcon(nscls,scn)
            call sedbudget(ius,dtime,timedr,nscls,sload,bload,bedk)
	  end if

          write(6,*) 'Sediment model initialized...'
	endif

! -------------------------------------------------------------------
! Normal call in the time loop
! -------------------------------------------------------------------

!       -------------------------------------------------------------
!       Get actual simulation time and time step
!       -------------------------------------------------------------
        call get_act_dtime(dtime)
	call get_timestep(dt)
        timedr = dt

!       -------------------------------------------------------------
!       Reset bottom height variation if it = adjtime
!       -------------------------------------------------------------
	call resetsedi(dtime,adjtime,nscls,bh,scn)

!       -------------------------------------------------------------
!       Compute bed slope gradient
!       -------------------------------------------------------------
        call tvd_grad_2d(hkv,gdx,gdy,v1v)

!       -------------------------------------------------------------------
!       Compute sediment transport on nodes
!       -------------------------------------------------------------------
        call sed_loop(timedr,nscls,gs,hzoff,totbed,riph,ripl,scn,eps,
     +                sedx,sedy,gdx,gdy,tao,gskm,percbd,bedn,salref,
     +                temref,wsink,sload,sflx)
	taubot = tao

!       -------------------------------------------------------------------
!       Zero erosion-deposition in the area with rocks (krocks(k) = 1)
!       -------------------------------------------------------------------
        call sed_on_rocks(nscls,krocks,sedx,sedy,sload,sflx)

!       -------------------------------------------------------------------
!       Compute bedload transport
!       -------------------------------------------------------------------
        call bedload(nscls,sedx,sedy,timedr,bedn,bh,bflx,bload)

!       -------------------------------------------------------------------
!       Transport and diffusion for each sediment class
!       -------------------------------------------------------------------
	call scn_compute(idsedi,nscls,sload,wsink,eps,scn)

!       -------------------------------------------------------------------
!       Update suspended flux after transp/diff
!       -------------------------------------------------------------------
	call update_sflx(nscls,sload,timedr,sflx)

!	-------------------------------------------------------------------
!	Update bed configuration and compute bed elevation change      
!	-------------------------------------------------------------------
        call bedman(gs,nscls,timedr,bflx,sflx,bh,gskm,bdh,percbd,bedn)

!	-------------------------------------------------------------------
!       Smooth node total depth variation
!	-------------------------------------------------------------------
        call smooth_node(bh,bdh,smooth,gdx,gdy,angrep)

!	-------------------------------------------------------------------
!       Update total depth and velocities
!	-------------------------------------------------------------------
        call upedepth(bdh)

!       -------------------------------------------------------------------
!       Compute total suspended concentration and update water density
!       -------------------------------------------------------------------
        call totcon(nscls,scn)

!       -------------------------------------------------------------------
!       Compute fraction of mud and sand in surface layer (just for output)
!       -------------------------------------------------------------------
	call mud_sand_fract(nlbdim,nscls,nkn,nbcc,percbd,percc,percs)

!       -------------------------------------------------------------------
!       Write of results (files SED and SCO)
!       -------------------------------------------------------------------
	call wr_sed_output(dtime,bh,gskm,tao,percc,percs,totbed)

!       -------------------------------------------------------------------
!       Compute sediment budget
!       -------------------------------------------------------------------
	if( bbudget ) then
          call sedbudget(ius,dtime,timedr,nscls,sload,bload,bedk)
	end if

!       -------------------------------------------------------------------
!       Compute variables for AQUABC model
!       -------------------------------------------------------------------
	call get_sedim_prop(nscls,bedn,sload)

! -------------------------------------------------------------------
! End of routine
! -------------------------------------------------------------------

	end subroutine sedi

! ********************************************************************

        subroutine initsed

	use para

        implicit none

	logical, save :: binit = .false.

	if( binit ) return

	binit = .true.

! DOCS  START   P_sediment
! 
! The following parameters activate the sediment transport module 
! and define the sediment grainsize classes to the simulated.

! |sedtr|	Sediment transport module section name.

        call sctpar('sedtr')             !sets default section
        call sctfnm('sedtr')

! |isedi|	Flag if the computation on the sediment is done:
! 		\begin{description}
! 		\item[0] Do nothing (default)
! 		\item[1] Compute sediment transport
! 		\end{description}

        call addpar('isedi',0.)

! |idtsed|, |itmsed|	Time step and start time for writing to files sed e sco,
! 			the files containing sediment variables and suspended
! 			sediment concentration.

        call addpar('idtsed',0.)
        call addpar('itmsed',-1.)

! |sedgrs|	Sediment grainsize class vector [mm]. Values has be 
! 		ordered from the finest to the more coarse. \\
! 		|example: sedgrs = '0.1 0.2 0.3 0.4'|

	call para_add_array_value('sedgrs',0.)

! |irocks|	Element type where do not compute erosion-deposition
! 		(Default -1).

        call addpar('irocks',-1.)

! |sedref|	Initial sediment reference concentration [kg/m3]
! 		(Default 0).

        call addpar('sedref',0.)

! |sedhpar|	Sediment diffusion coefficient (Default 0).

        call addpar('sedhpar',0.)

! |adjtime|	Time for sediment initialization [s]. The sediment model needs a
! 		initialization time in which the system goes to a quasi steady
! 		state. When t = |adjtime| the bed evolution change in the output
! 		is reset. Keep in mind that |adjtime| has to be chosen case by 
! 		case in function of the morphology and the parameters used for 
! 		the simulation (Default 0).

        call addpar('adjtime',-1.)

! |percin|	Initial sediment distribution (between 0 and 1) for each 
! 		grainsize class. The sum of percin must be equal to 1. \\
! 		|example: percin = 0.25 0.25 0.25 0.25| \\
! 		If percin is not selected the model impose equal
! 		percentage for each grainsize class (percin = 1/nrs).
! 		In case of spatial differentiation of the sediment
! 		distribution set a number of percin equal to the number
! 		of grainsize classes per the number of area types. Element 
!		types should be numbered consecutively starting from 0.\\
! 		|example: percin = 0.25 0.25 0.25 0.25  \\
! 			 	   0.20 0.20 0.30 0.30  \\
! 				   0.45 0.15 0.15 0.15|

        !call addpar('percin',0.)
	call para_add_array_value('percin',0.)

! |tauin|	Initial dry density or TAUCE. In function of the value: \\
! 		0-50 : critical erosion stress (Pa) \\
! 		\textgreater 50  : dry bulk density of the surface (kg/m**3). \\
! 		In case of spatial differentiation set a number of tauin
! 		equal to the number of area type. Element types should be
!		numbered consecutively starting from 0. \\
! 		|example: tauin = '0.9 1.4 2.5 1.1'|

	call para_add_array_value('tauin',0.)

! |sedp|	File containing spatially varying initial sediment distribution 
!		for each grid node. Values are in percentage of each class and 
!		the file should be structured with number of columns equal the
!		number of grainsize classes and the number of row equal the
!		number of nodes (the order should follow the internal node 
!		numbering).

        call addfnm('sedp',' ')

! |sedt|	File containing spatially varying initial critical erosion 
! 		stress (Pa) or dry bulk density (kg/m3). One value for each
!		node (the order should follow the internal node numbering). 

        call addfnm('sedt',' ')

! |sedcon|	File containing the additional constants used in sediment 
! 		model. These parameters are usually set to the indicated 
! 		default values, but can be customized for each sediment 
! 		transport simulation. The full parameter list together with 
! 		their default value and brief description is reported in 
! 		Table \ref{tab:table_sedcon}. Most of the parameter, 
! 		especially the ones for the cohesive sediments,	have been 
! 		calibrated for the Venice Lagoon. For more information about 
! 		these parameters please refer to Neumeier et 
! 		al. \cite{urs:sedtrans05} and Ferrarin et al. 
! 		\cite{ferrarin:morpho08}.

        call addfnm('sedcon',' ')

! DOCS  FILENAME        Additional parameters

! The full parameter list is reported in Table \ref{tab:table_sedcon}.
! An example of the settings for the |sedcon| file is given in
! \Fig\figref{turbulence}. Please note that is not necessary to
! define all parameters. If not defined the default value is imposed.
! 
! \begin{figure}[ht]
! \begin{verbatim}
! IOPT = 3
! SURFPOR = 4
! DOCOMPACT = 1
! \end{verbatim}
! \caption{Example of the secon file.}
! \label{fig:sedcon}
! \end{figure}

! \begin{table}[ht]
! \caption{Additional parameter for the sediment transport model to be set 
! in the sedcon file.}
! \begin{tabular}{lcl} \hline
! Name & Default value & Description \\ \hline
! CSULVA  & 159.4 & Coefficient for the solid transmitted stress by Ulva \\
! TMULVA  & 1.054d-3 & Threshold of motion of Ulva [Pa] \\
! TRULVA  & 0.0013 & Threshold of full resuspension of Ulva [Pa] \\
! E0      & 1.95d-5 & Minimum erosion rate \\
! RKERO   & 5.88 & Erosion proportionality coefficient \\
! WSCLAY  & 5.0 & Primary median Ws class (in the range 1:NBCONC) \\
! CDISRUPT& 0.001 & Constant for turbulent floc disruption during erosion \\
! CLIM1   & 0.1 & Lower limit for flocculation [kg/m3] \\
! CLIM2   & 2.0 & Limit between simple/complex flocculation [kg/m3] \\
! KFLOC   & 0.001 & Constant K for flocculation equation \\
! MFLOC   & 1.0 & Constant M for flocculation equation \\
! RHOCLAY & 2600.0 &  Density of clay mineral \\
! CTAUDEP & 1.0 & Scaling factor for TAUCD \\
! PRS     & 0.0 & Resuspension probability [0-1] \\
! RHOMUD  & 50.0 & Density of the freshly deposited mud \\
! DPROFA  & 470.0 & Constants for density profile \\
! DPROFB  & 150.0 & A : final deep density \\
! DPROFC  & 0.015 & Define the shape (in conjunction with B and C) \\
! DPROFD  & 0.0 & Aux parameter for density profile \\
! DPROFE  & 0.0 & Aux parameter for density profile \\
! CONSOA  & 1d-5 & time constant of consolidation \\
! TEROA   & 6d-10 & Constant for erosion threshold from density \\
! TEROB   & 3.0 & Aux parameter for erosion threshold from density \\
! TEROC   & 3.47 & Aux parameter for erosion threshold from density \\
! TEROD   & -1.915 & Aux parameter for erosion threshold from density \\
! KCOES   & 0.15 & Fraction of mud for sediment to be cohesive \\
! CDRAGRED& -0.0893 & Constant for the drag reduction formula \\
! Z0COH   & 2.0D-4 & Bed roughness length for cohesive sediments \\
! FCWCOH  & 2.2D-3 & Friction factor for cohesive sediments \\
! LIMCOH  & 0.063 & Limit of cohesive sediment grainsize [mm] \\
! SMOOTH  & 1.0 & Smoothing factor for morphodynamic \\
! ANGREP  & 32.0 & Angle of repose \\
! IOPT    & 5 & Sediment bedload transport formula option number \\
! MORPHO  & 1.0 & Morphological acceleration factor \\
! RHOSED  & 2650.0 & Sediment grain density \\
! POROS   & 0.4 & Bed porosity [0-1] \\
! SURFPOR & 0.6 & Bed porosity of freshly deposited sand [0-1] \\
! DOCOMPACT& 0.0 & If not zero, call COMPACT routine \\ \hline
! \end{tabular}
! \label{tab:table_sedcon}
! \end{table}
! 
! DOCS  END

	end

! ********************************************************************
! SUBROUTINE READGS
! This subroutine reads the simulation sediment parameter from the .str file 
! It's called by the routine nlsh2d

        subroutine readsed

	use para
        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05

        implicit none

        real		  :: getpar
        double precision  :: dgetpar
        integer 	  :: i
        integer		  :: nrs	!number of grainsize classes
        integer 	  :: npi	!
        integer 	  :: nps
        integer 	  :: nrst	!number of sediment type
	real, allocatable :: pppp(:)

!       --------------------------------------------------
!       Initialize variables
!       --------------------------------------------------

        call initsed

        npi  = 1
        nps  = 0
        nrs  = 0
        nrst = 0

!       --------------------------------------------------
!       Read all section sedtr in str
!       --------------------------------------------------
	call nrdins('sedtr')

!       --------------------------------------------------
!       Reads sediment grainsize vector with dimension nrs
!	and convert it from mm to m
!       --------------------------------------------------
	call para_get_array_size('sedgrs',nrs)
	if( nrs == 0 ) goto 30
	allocate(gsc(nrs))
	call para_get_array_value('sedgrs',nrs,nrs,gsc)
        if (gsc(1) .ge. 100d0) goto 36
        gsc = gsc*0.001	

!       --------------------------------------------------
!       Reads initial erosion threshold vector
!       --------------------------------------------------
	call para_get_array_size('tauin',nrst)
	if( nrst == 0 ) goto 31
	allocate(tue(nrst))
	call para_get_array_value('tauin',nrst,nrst,tue)
	if( nrst == 1 ) tue = tue(1)

!       --------------------------------------------------
!       Reads sediment percentage matrix for each type
!       --------------------------------------------------
	call para_get_array_size('percin',nps)
	if( nps == 1 ) then
	  npi = 1
	else
          npi = nps/nrs
	  if ( mod(nps,nrs) /= 0 ) goto 35
	end if
	allocate(pppp(nps))
	allocate(prin(npi,nrs))
        call para_get_array_value('percin',nps,nps,pppp)
	if( nps == 1 ) then
	  prin = pppp(1)
	else
	  prin = reshape(pppp(1:nps),[npi,nrs],order=[2,1])
	end if

!       --------------------------------------------------
!       Case with only one sediment class
!       --------------------------------------------------
        if( nrs == 1 ) then
          npi = 1
          nps = 1
          prin(npi,nps) = 1
        end if

!       --------------------------------------------------
!       Case with no initial distribution selected
!	impose the same percentage for each class
!       --------------------------------------------------
        if ( nps == 1 ) then
          do i = 1,nrs
            prin(npi,i) = 1./nrs
          end do
        end if

!       --------------------------------------------------
!       Stores parameters in variable sedpa()
!       --------------------------------------------------
        sedpa(1) = dgetpar('isedi')
        sedpa(2) = dgetpar('sedref')
        sedpa(3) = dgetpar('sedhpar')
        sedpa(4) = nrs
        sedpa(5) = npi
        sedpa(6) = nrst
        sedpa(7) = dgetpar('adjtime')
        sedpa(8) = dgetpar('irocks')
        return

  30    continue
        write(6,*) 'No sediment class selected'
        write(6,*) 'in the .str file'
        stop 'error stop : readsed'

  31    continue
        write(6,*) 'No critical shear stress given'
        write(6,*) 'in the .str file'
        stop 'error stop : readsed'

  35    continue
        write(6,*) 'Dimension error for percin'
        write(6,*) 'nps   :',nps
        write(6,*) 'nscls :',nrs
        stop 'error stop : readsed'

  36    continue
        write(6,*) 'First grainsize class > 10 cm'
        write(6,*) 'grainsize :',gsc(1)
        stop 'error stop : readsed'

        end subroutine readsed

! ********************************************************************
! SUBROUTINE READCONST
! This subroutine read constants from file sedcon.
! The format of the file is: variable = value

        subroutine readsedconst

        use mod_sediment_para
	use mod_sedtrans05

        implicit none

        integer 		:: nrdnxt
	integer			:: iw
	integer			:: ioff
        integer 		:: ifileo
        integer 		:: iunit
        real 			:: value
	double precision 	:: dvalue
        character*80 		:: file
	character*80            :: name
	character*80            :: text
	character*80            :: line
	integer, parameter      :: nc = 38	  !number of variables
	character*80, dimension(nc)     :: cname  !variable name
        double precision, dimension(nc) :: cvalue !variable value
        integer			:: is,j

	iw    = 0
	iunit = 0
	value = 0.0
        file  = ''

!       --------------------------------------------------------
!       Assigns constants to cname and default values to cvalue
!       --------------------------------------------------------

        cname(1) = 'CSULVA'	!coefficient for the solid transmitted stress by Ulva
        cvalue(1) = CSULVA
        cname(2) = 'TMULVA'	!threshold of motion of Ulva (Pa)
        cvalue(2) = TMULVA
        cname(3) = 'TRULVA'	!threshold of full resuspension of Ulva (Pa)
        cvalue(3) = TRULVA
        cname(4) = 'E0'		!Minimum erosion rate
        cvalue(4) = E0
        cname(5) = 'RKERO'	!Erosion proportionality coefficient
        cvalue(5) = RKERO
        cname(6) = 'WSCLAY'	!primary median Ws class (must be in the range 1:NBCONC)
        cvalue(6) = WSCLAY
        cname(7) = 'CDISRUPT'	!constant for turbulent floc disruption during erosion
        cvalue(7) = CDISRUPT
        cname(8) = 'CLIM1'	!lower limit for flocculation (kg/m**3)
        cvalue(8) = CLIM1
        cname(9) = 'CLIM2'	!limit between simple and complex flocculation equation (kg/m**3)
        cvalue(9) = CLIM2
        cname(10) = 'KFLOC'   	!constant K for flocculation equation
        cvalue(10) = KFLOC
        cname(11) = 'MFLOC'   	!constant M for flocculation equation
        cvalue(11) = MFLOC
        cname(12) = 'RHOCLAY' 	!density of clay mineral
        cvalue(12) = RHOCLAY
        cname(13) = 'CTAUDEP' 	!scaling factor for TAUCD
        cvalue(13) = CTAUDEP
        cname(14) = 'PRS'     	!Resuspension probability (range 0-1)
        cvalue(14) = PRS
        cname(15) = 'RHOMUD'  	!Density of the freshly deposited mud
        cvalue(15) = RHOMUD	!50d0
        cname(16) = 'DPROFA'	!Constants for density profile
        cvalue(16) = DPROFA	!FINAL(I) = A - B * EXP(-C*ABOVE(I)) - D*EXP(-E*ABOVE(I))
        cname(17) = 'DPROFB'	!A : final deep density
        cvalue(17) = DPROFB	!A-B-D : final surface density
        cname(18) = 'DPROFC'	!define the shape (in conjunction with B and C)
        cvalue(18) = DPROFC
        cname(19) = 'DPROFD'
        cvalue(19) = DPROFD
        cname(20) = 'DPROFE'
        cvalue(20) = DPROFE
        cname(21) = 'CONSOA'	!time constant of consolidation
        cvalue(21) = CONSOA
        cname(22) = 'TEROA'	!Constants for erosion threshold from density and overlaying mass
        cvalue(22) = TEROA
        cname(23) = 'TEROB'
        cvalue(23) = TEROB
        cname(24) = 'TEROC'
        cvalue(24) = TEROC
        cname(25) = 'TEROD'
        cvalue(25) = TEROD
        cname(26) = 'KCOES'	!Fraction of mud for sediment to be cohesive
        cvalue(26) = KCOES
        cname(27) = 'CDRAGRED'	!constant for the drag reduction formula
        cvalue(27) = CDRAGRED
        cname(28) = 'Z0COH'	!BED ROUGHNESS LENGHT FOR COHESIVE SEDIMENTS
        cvalue(28) = Z0COH
        cname(29) = 'FCWCOH'	!FRICTION FACTOR FOR COHESIVE SEDIMENTS
        cvalue(29) = FCWCOH
        cname(30) = 'LIMCOH'	!Limit of cohesive sediment grainsize
        cvalue(30) = LIMCOH
        cname(31) = 'SMOOTH'	!smoothing factor for morphodynamic
        cvalue(31) = SMOOTH
        cname(32) = 'ANGREP'	!angle of repose
        cvalue(32) = ANGREP
        cname(33) = 'IOPT'	!SEDIMENT TRANSPORT FORMULA OPTION NUMBER
        cvalue(33) = IOPT
        cname(34) = 'MORPHO'	!Morphological acceleration factor
        cvalue(34) = MORPHO
        cname(35) = 'RHOSED'    !sediment grain density
        cvalue(35) = RHOSED
        cname(36) = 'POROS'	!bed porosity [0,1]
        cvalue(36) = POROS
        cname(37) = 'SURFPOR'	!bed porosity of freshly deposited sand [0,1]
        cvalue(37) = SURFPOR
        cname(38) = 'DOCOMPACT'	!if not zero, call COMPACT routine
        cvalue(38) = DOCOMPACT

!       --------------------------------------------------------
!       Gets sedcon file name
!       --------------------------------------------------------

        call getfnm('sedcon',file)
        write(6,*)

        if( file .ne. ' ' ) then
           iunit = ifileo(0,file,'form','old')
           write(6,*)'Sediment constants initialized from file: ',
     $file(1:40)
           if( iunit .le. 0 ) then
            write(6,'(a,a)') 'filename: ',file(1:65)
            stop 'error stop readcon: Cannot open file'
           end if

!         --------------------------------------------------------
!         Reads first line in file
!         --------------------------------------------------------

!ggu          read(iunit,'(a)') line
!ggu          ioff = 1

!         --------------------------------------------------------
!         Loop on lines
!         --------------------------------------------------------

	  call nrdini(iunit)

	  do
	    iw = nrdnxt(name,dvalue,text)
	    if( iw .le. 0 ) exit
	    if( iw .ne. 1 ) stop 'error stop readsedconst: wrong item'
            do j = 1,nc
              if (name .eq. cname(j)) cvalue(j) = dvalue
            end do
	  end do

	  close(iunit)

!    1     continue
!            iw = nrdnls(name,dvalue,text,line,ioff)
!	    value = real(dvalue)
!            if( iw .le. 0 ) then
!              read(iunit,'(a)',end=2) line
!              ioff = 1
!            else
!              if( iw .eq. 1 ) text = ' '
!              if( iw .eq. 2 ) value = 0
!              do j = 1,nc
!                if (name .eq. cname(j)) cvalue(j) = value
!              end do
!            end if
!            goto 1
!    2     continue

        end if

!       --------------------------------------------------------
!       Assigns new values to constants
!       --------------------------------------------------------
         
        CSULVA    = cvalue(1)
        TMULVA    = cvalue(2)
        TRULVA    = cvalue(3)
        E0        = cvalue(4)
        RKERO     = cvalue(5)
        WSCLAY    = int(cvalue(6))
        CDISRUPT  = cvalue(7)
        CLIM1     = cvalue(8)
        CLIM2     = cvalue(9)
        KFLOC     = cvalue(10)
        MFLOC     = cvalue(11)
        RHOCLAY   = cvalue(12)
        CTAUDEP   = cvalue(13)
        PRS       = cvalue(14)
        RHOMUD    = cvalue(15)
        DPROFA    = max(cvalue(16),RHOMUD)
        DPROFB    = cvalue(17)
        DPROFC    = cvalue(18)
        DPROFD    = cvalue(19)
        DPROFE    = cvalue(20)
        CONSOA    = cvalue(21)
        TEROA     = cvalue(22)
        TEROB     = cvalue(23)
        TEROC     = cvalue(24)
        TEROD     = cvalue(25)
        KCOES     = cvalue(26)
        CDRAGRED  = cvalue(27)
        Z0COH     = cvalue(28)
        FCWCOH    = cvalue(29) 
        LIMCOH    = cvalue(30) 
        SMOOTH    = cvalue(31) 
        ANGREP    = cvalue(32) / (45./atan (1.))	!from deg to rad
        IOPT      = int(cvalue(33))
        MORPHO    = cvalue(34)
        RHOSED    = cvalue(35)
        POROS     = cvalue(36)
        SURFPOR   = cvalue(37)
	DOCOMPACT = int(cvalue(38))

!       --------------------------------------------------------
!       Write constant value on the screen
!       --------------------------------------------------------

        write(6,*)'Constants for sediment transport:'
	write(6,*)
        do j = 1,nc
          write(6,44)cname(j),cvalue(j)
        end do
        write(6,*)

44      format(3x,a10,f12.7)

        end subroutine readsedconst

! ********************************************************************
! SUBROUTINE INIBED
! This subroutine set the initial bed conformation in function of
! sediment grainsize classes and the initial percentage defined in
! the .str file. The distribution is the same in each bed level (no
! vertical stratification). 

        subroutine inibed(gs,nscls,nbed,thick,gskm,percbd,bedn,krocks)

        use basin, only : nkn,nel,iarv,nen3v
        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05

        implicit none

        double precision, intent(in) 	:: gs(nscls)    !sediment grain diameter (m)
        integer, intent(in) 		:: nscls	!number of grainsize classes
        integer, intent(in) 		:: nbed		!initial number of bed layers
        real, intent(in) 		:: thick	!initial thickness [m]

        real, intent(out) 		:: gskm(nkn)	!average grainsize per node [m]
        double precision, intent(out) 	:: percbd(nlbdim,nscls,nkn) !fraction of sediment [0,1]
        double precision, intent(out)	:: bedn(nlbdim,3,nkn)  !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer, intent(out)		::  krocks(nkn)	!node index with rocks

        double precision 		:: gsa
	double precision                :: gss
	double precision 		:: ptot
        double precision 		:: rhosa	!initial erosion threshold
        double precision 		:: rhossand	!function for density from % sand
        double precision 		:: pcoes	!% of fine sediments
	integer 			:: itype        !element code number
        integer 			:: its
        real 				:: pers(nkn,nscls) !percbd initialized from file
        real 				:: tuek(nkn,1)     !tauce initialized from file
	double precision		:: peps,pmin,pmax,ppmin,ppmax
	integer 			:: irocks
        integer 			:: npi
        integer 			:: k,ib,is,ie,ii

!       -------------------------------------------------------------------
!       Initialize sediment distribution in function of bed type
!	and set krocks (=1 where no erosion-deposition take place)
!       -------------------------------------------------------------------

	krocks = 0
        npi    = int(sedpa(5))
	irocks = int(sedpa(8))

        do ie = 1,nel
          itype = iarv(ie)
          its = itype + 1
          if ( npi.eq.1 ) its = 1

          do ii = 1,3
            k = nen3v(ii,ie)
            tuek(k,1) = tue(its)
	    if (itype == irocks) krocks(k) = 1
            do is = 1,nscls
              pers(k,is) = prin(its,is)
            end do
          end do
        end do

!       -------------------------------------------------------------------
!       Initialize sediment fraction and tuek from external file (fem format)
!       -------------------------------------------------------------------

	call init_file_fem('sed tuek','sedt',1,nkn,tuek)
	call init_file_fem('sed dist init','sedp',nscls,nkn,pers)

!       -------------------------------------------------------------------
!       Initialize bed thickness if thick .ne. 0
!       -------------------------------------------------------------------

        if (thick.gt.0.) then
          do k = 1,nkn
            bedn(1,1,k) = 0.
            do ib = 2,nbed
             bedn(ib,1,k) = bedn(ib-1,1,k) + thick
            end do
          end do
        end if

!       -------------------------------------------------------------------
!       Check consistency of pers
!       -------------------------------------------------------------------

	peps = 0.3
	pmin = 1. - peps
	pmax = 1. + peps
	ppmin = 1.
	ppmax = 0.

	write(6,*) 'inibed: scaling grainsize percentage'

        do k = 1,nkn
	  ptot  = 0.d0
          do is = 1,nscls
            ptot = ptot + pers(k,is)
          end do
          if(ptot.lt.pmin.or.ptot.gt.pmax) go to 130
	  ppmin = min(ppmin,ptot)
	  ppmax = max(ppmax,ptot)
	  pers(k,:) = pers(k,:) / ptot
	end do

	write(6,*) 'inibed: percentage scaled: ',ppmin,ppmax

	pers = max(0.,pers)
	pers = min(1.,pers)

!       -------------------------------------------------------------------
!       Initialize bed characteristics and compute initial average grainsize
!       -------------------------------------------------------------------

	peps = 0.01
	pmin = 1. - peps
	pmax = 1. + peps

        do k = 1,nkn
          gsa   = 0.d0
	  ptot  = 0.d0
	  pcoes = 0.d0
          do is = 1,nscls
            do ib = 1,nbed
              percbd(ib,is,k) = pers(k,is)
            end do
            ptot = ptot + percbd(1,is,k)
            gss = gs(is)
	    if (gss .ge. 0.10d0 .and. is .gt. 1) gss = gs(is-1)
            gsa = gsa + gss*percbd(1,is,k)
            if (gs(is).lt.limcoh) pcoes = pcoes + percbd(1,is,k)
          end do
          if(ptot.lt.pmin.or.ptot.gt.pmax) go to 130
          gskm(k) = real(gsa)

          rhosa = rhossand(pcoes,1.2d0)
          if(tuek(k,1).gt.0.) rhosa = tuek(k,1)

          call bedini(bedn(1,1,k),nlbdim,rhosa)
        end do

	return
 130    continue
        write(6,*) 'Error in computing the sediment fraction:'
        write(6,*) 'node: ',k,'  total percentage:',ptot
        do is = 1,nscls
          write(6,*) 'classnumber:',is,'percentage:',pers(k,is)
        enddo
        stop 'error stop inibed: percbd'

        end subroutine inibed

! *******************************************************
! Initilize array from file txt

	subroutine init_file_txt(name,nkn,nvar,var)

        implicit none

        character*(*), intent(in)	:: name      !file name
	integer, intent(in) 		:: nkn
	integer, intent(in)		:: nvar
	real, intent(inout)		:: var(nkn,nvar)

	integer 			:: k,j
        character*80			:: file

        call getfnm(name,file)
        if( file .eq. ' ' ) return      !nothing to initialize
	write(6,*) 'Reading sediment from file: ',file

	open(199,file=file,form='formatted')

	do k = 1,nkn
  	  read(199,*) (var(k,j),j=1,nvar)
	end do

        end subroutine init_file_txt

! *******************************************************


        subroutine init_file_fem(what,file_init,nvar,nkn,val)

c initialization of conz from file

        implicit none

        character*(*) what
        character*(*) file_init
        integer nvar
        integer nkn
        real val(nkn,nvar)

	integer nlvddi,nlv
        double precision dtime
        real val0(nvar)
	character*80 file

        call getfnm(file_init,file)
        if( file == ' ' ) return

	dtime = 0.
	nlvddi = 1
	nlv = 1
	val0 = 0.	!not needed

        call tracer_file_init(what,file_init,dtime
     +                          ,nvar,nlvddi,nlv,nkn,val0,val)

        end

! *******************************************************

! Initilize array from file
! TO BE DONE ( should replace init_file_txt)

	subroutine init_file_sed(name,nkn,nvar,var)

        use intp_fem_file
        use mod_sediment_para

        implicit none

        character*(*) name      !file name
	integer nkn
	integer nvar
	real var(nkn,nvar)

        integer nb,k
        integer nintp,np,lmax,ibc
        integer idsed
        double precision dtime
        integer nodes(1)
        real sconst(1)
        character*10 whats
        character*80 file

        dtime  = 0
        nodes  = 0
        nintp  = 2
        np     = nkn
        lmax   = 0
        ibc    = 0                         !no lateral boundary
        whats  = 'sedp init'
        sconst = 0.

        call getfnm(name,file)
        if( file .eq. ' ' ) return      !nothing to initialize

        write(6,*) 'Initializing sediment grainsize...'
        call iff_init(dtime,file,nvar,np,lmax,nintp
     +                          ,nodes,sconst,idsed)
        call iff_set_description(idsed,ibc,whats)

        lmax = 1
        call iff_read_and_interpolate(idsed,dtime)
        call iff_time_interpolate(idsed,dtime,0,np,lmax,var)

        call iff_forget_file(idsed)

        write(6,*) 'Initial sediment grainsize read from file : '
        write(6,*) trim(name)

        end

! *******************************************************
! FUNCTION RHOSSAND
! returns the dry bulk density from sand fraction (Allersma, 1988)

        function rhossand(pcoes,consc)

        use mod_sediment_para

        implicit none

        double precision pcoes              !fraction of fine sediments
        double precision consc              !consolidation coefficient [0-2.4]

        double precision psand              !sand fraction [0,1]
        double precision rhossand           !density function [kg/m3]

        psand = 1.d0
        if(pcoes.gt.0.d0) then 
          psand = 1.d0 - pcoes
	  psand = max(psand,0.D+0)		!handle psand < 0
          rhossand = 480.d0*consc+(1300.d0-280.d0*consc)*(psand**0.8)
        else
          rhossand = RHOSED * (1.d0 - POROS)
        endif

        end

! *******************************************************
! SUBROUTINE ZVEL
! This subroutine computes z, the height of uz (the main current 
! velocity) above seafloor, assuming a log profile of the velocity 
! on the depth, based on the law: U = (USTAR/K)*LN(Z/ZO)

	subroutine zvel(uz,d,gd,z0bk,z)

	implicit none

! --- input variables
	double precision d	!water depth (m)
	double precision uz	!current at height z above the seafloor (m/s)
	double precision gd	!sediment grain diameter (m)
	real z0bk		!bottom roughness length from sedtrans [m]

! --- output variables
	double precision z	!height of uz above seafloor [m]

! --- local variables
	double precision ustar  !shear velocity 
	real br 	  	!bottom roughness height [m]
	real z0   		!bottom roughness length [m]
        double precision :: dz0 !d/z0
        double precision, parameter :: dz0min = 1.1     !min val for dz0
        double precision, parameter :: k = 0.4          !von karman costant

        br = 2.5*real(gd)		!see amos article
        z0 = br/30.
	z0 = max(z0,z0bk)
        if (uz .le. 0.d0) then
	  z = d/2.d0
	else
          dz0 = dz0min                          !GGUZ0
          if( z0 > 0. ) dz0 = d/z0
          if( dz0 < dz0min ) dz0 = dz0min
          ustar = (uz*k*dz0)/(dz0*(log(dz0)-1.d0)+1.d0)
          !ustar = (uz*k*d/z0)/(d/z0*(log(d/z0)-1.d0)+1.d0)
          z = z0*exp(uz*k/ustar)
	end if

	end

! *******************************************************
! SUBROUTINE GETMD
! This subroutine compute velocity and direction from u and v

        subroutine getmd(u,v,m,d)

        implicit none

! --- input variables
        real u,v		!x,y components

! --- output variables
        double precision m      !velocity
        double precision d	!direction

! --- local variables
        real rad
        real alfa

        rad = 90./asin(1.)

        m = sqrt( u**2 + v**2 )		!get m

        if( u .eq. 0.0 ) then
          alfa = 0.
          if( v .lt. 0.0 ) alfa = 180.
	else
          alfa=atan(v/u)*rad
        end if

        if(u.gt.0.and.v.gt.0.or.u.gt.0.and.v.lt.0) alfa=90.-alfa
        if(u.lt.0.and.v.gt.0.or.u.lt.0.and.v.lt.0) alfa=270.-alfa

        !if(alfa.le.90.and.alfa.ge.0)   alfa=90. - alfa
        !if(alfa.le.360.and.alfa.ge.90) alfa=450. - alfa

        d = alfa

        end

! *******************************************************
! SUBROUTINE VMIXCOEF
! This subroutine compute vertical diffusion coefficient for sediment
!	  -1- Algebric Van Rijn 1993
!	  -2- Viscosity from GOTM and beta from Van Rijn 1993
!	  -3- Viscosity from GOTM and beta from Villatoro et al 2010
!         -4- Diffusivity from GOTM 

        subroutine vmixcoef(dtype,l,dep,h,wsink,ustc,ustcw,ub,dcw,ht,
     +			    per,dx,difv,visv,epss)

        use levels

        implicit none

! --- input variables

	integer dtype			!type of diffusion algoritm
        integer l			!number of levels
        real dep(nlvdi)		!thickness of layer [m]
        double precision h		!total water depth [m]
        double precision wsink		!settling velocity [m/s]
        double precision ustc		!current shear velocity [m/s]
        double precision ustcw		!combined current and wave shear velocity [m/s]
        double precision ub 		!wave orbital velocity [s]
        double precision dcw		!height of the wave-current boundary layer
        double precision ht		!significant wave height [m]
        double precision per		!significant wave persion [s]
        double precision dx		!dimensionless particle diameter
        real visv(0:nlvdi)		!water viscosity [from GOTM] 
        real difv(0:nlvdi)		!water diffusivity [from GOTM] 

! --- output variables

        double precision epss(0:nlvdi)

! --- local variables

        integer m
        double precision d		!depth from the bottom [m]
        double precision bet		!ratio of sediment and fluid mixing [1-1.5]
        double precision epsc		!current related vertical mixcoef 
        double precision epscmax	!max current related vertical mixcoef 
        double precision epsw		!wave related vertical mixcoef 
        double precision epswb		!wave related vertical mixcoef near the bed
        double precision epswmax	!wave related vertical mixcoef in the upper layer

	if (dtype .eq. 1) then

!         ------------------------------------------------------------
!         Algebraic formulation of Van Rijn 1993
!         ------------------------------------------------------------

          d = 0.d0
          bet = 1.d0 + 2.d0* (wsink/ustc)**2
          bet = dmax1(bet,1.d0)
          bet = dmin1(bet,1.5d0)
          epscmax = 0.25d0 * bet * ustc * h
          epswb = 0.004d0 * dx * dcw * ub
          epswmax = 0.035d0*h*ht/per
          if(ub.lt.0.1d0) epswmax = 0.d0

!         ------------------------------------------------------------
!         Loop over layers
!         ------------------------------------------------------------

          do m = l,1,-1
            epss(m) = 0.
            d = d + dep(m)

!           ------------------------------------------------------------
!           Current related mixing
!           ------------------------------------------------------------

            if (d .lt. 0.5d0*h) then
              epsc = epscmax - epscmax*(1.d0-(2.d0*d/h)**2)
            else
              epsc = epscmax
            end if

!           ------------------------------------------------------------
!           Wave related mixing
!           ------------------------------------------------------------

            if (d .le. dcw) then
              epsw = epswb
            elseif(d .gt. dcw .and. d .le. 0.5d0*h) then
              epsw = epswb + (epswmax-epswb)*((d-dcw)/(0.5d0*h-dcw))
            else
              epsw = epswmax
            end if 

!           ------------------------------------------------------------
!           Combined current and wave mixing coefficient
!           ------------------------------------------------------------

            epss(m) = real(sqrt(epsc**2 + epsw**2))

          end do

          epss(0) = epss(1)

	else if (dtype .eq. 2) then

!         ------------------------------------------------------------
!         Use viscosity from GOTM and beta factor from Van Rijn 1993
!	  Beta is function of settling velocity and shear velocity
!         ------------------------------------------------------------

	  bet = 1.
          if( ustcw /= 0. ) bet = 1.d0 + 2.d0* (wsink/ustcw)**2
          bet = dmax1(bet,1.d0)
          bet = dmin1(bet,1.5d0)
	  do m = 0,l
	    epss(m) = bet * visv(m)
	  end do

	else if (dtype .eq. 3) then

!         ------------------------------------------------------------
!         Use viscosity from GOTM and beta factor from Villatoro et al 2010
!	  Beta is function of grainsize
!         ------------------------------------------------------------

	  bet = 0.33d0*dx - 0.34d0
	  do m = 0,l
	    epss(m) = bet * visv(m)
	  end do

	else if (dtype .eq. 4) then

!         ------------------------------------------------------------
!         Use diffusivity from GOTM 
!         ------------------------------------------------------------

	  do m = 0,l
	    epss(m) = difv(m)
	  end do

	else
	  write(6,*) 'unknown dtype : ',dtype
          stop 'error stop vmixcoef'
	end if

        end

! *******************************************************
! SUBROUTINE SED_LOOP

        subroutine sed_loop(timedr,nscls,gs,hzoff,totbed,riph,ripl,
     $	scn,eps,sedx,sedy,gdx,gdy,tao,gskm,percbd,bedn,
     $  salref,temref,wsink,sload,sflx)

	use mod_waves
	use mod_layer_thickness
	use mod_roughness
	use levels
	use basin, only : nkn
        use mod_sediment
        use mod_sediment_para
        use mod_bound_geom

        implicit none

! -------------------------------------------------------------
! local variables
! -------------------------------------------------------------
        integer nscls			!number of grainsize classes
        integer is,k,l			!counters
	integer lmax			!bottom level
        double precision dl		!thickness of bottom layer 
        double precision gs(nscls)	!SEDIMENT GRAIN DIAMETER (M)
        double precision ws(nscls)	!SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        double precision scn(nlvdi,nkn,nscls)   !cohesive suspended sediment conc (kg/m3)
        double precision riph(nkn)	!ripple height [m]
        double precision ripl(nkn)	!ripple length [m]
        real gskm(nkn)		!AVERAGE SEDIMENT GRAIN DIAMETER ON NODES (M)
        double precision eps(0:nlvdi,nkn,nscls)	!vertical mixing coefficient
        real u,v			!x and y current components
        real gdx(nkn),gdy(nkn)	!slope gradients
        real tao(nkn)		!wave-current shear stress
        double precision sflx(nscls,nkn)     !flux of suspend sediment [m3/m2]
        double precision sedx(nscls,nkn)	!x bedload component [kg/ms]
        double precision sedy(nscls,nkn)    	!y bedload component [kg/ms]
        double precision percbd(nlbdim,nscls,nkn)        !fraction of sediment [0,1]
        double precision bedn(nlbdim,3,nkn)	!bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        double precision scns(nscls)    	!suspended sediment conc at bottom layer (kg/m3)
	double precision hzoff
	real z0b
	real totbed(nkn)			!total bedload transport (kg/ms)
        real salref,temref			!salinity [psu] and temperature [C]
        double precision wsink(0:nlvdi,nkn,nscls) !settling velocity for suspended sediment

! -------------------------------------------------------------
! fem variables
! -------------------------------------------------------------

        real salt,temp			!salinity [psu] and temperature [C]

! -------------------------------------------------------------
! sedtrans05 variables
! -------------------------------------------------------------

        DOUBLE PRECISION D		!WATER DEPTH (M)
        DOUBLE PRECISION UZ		!AMBIENT CURRENT AT HEIGHT Z ABOVE THE SEAFLOOR (M/S)
        DOUBLE PRECISION Z		!HEIGHT OF UZ ABOVE SEAFLOOR
        DOUBLE PRECISION CDIR		!DIRECTION OF AMBIENT CURRENT (DEGREES TRUE)
        DOUBLE PRECISION HT		!WAVE HEIGHT (M)
        DOUBLE PRECISION MPER		!MEAN WAVE PERIOD (S)
        DOUBLE PRECISION PPER		!PEAK WAVE PERIOD (S)
        DOUBLE PRECISION WDIR		!WAVE PROPAGATION DIRECTION (DEGREES TRUE)
        DOUBLE PRECISION GD		!SEDIMENT GRAIN DIAMETER (M) 
        DOUBLE PRECISION BETA		!BED SLOPE (DEGREE)
        DOUBLE PRECISION TEM		!in degrees Celsius
        DOUBLE PRECISION SAL		!practical salinity scale (seawater=35)
        DOUBLE PRECISION TIMEDR		!Duration of call to Sedtrans (s)
        DOUBLE PRECISION Z0       	!BED ROUGHNESS LENGTH (M)
        DOUBLE PRECISION UW       	!bottom orbital velocity from WWM

! cohesive
        DOUBLE PRECISION AULVA		!PERCENTAGE OF AREA COVERED BY THE ALGAE 'ULVA' (%)
	double precision sload(nkn,nscls)	!suspended sediment load [kg/s]
	double precision sloads(nscls)	  	!suspended sediment load [kg/m²s]

	real ukbot(nkn),vkbot(nkn)
	real ddl(nkn)
	real areanode,area

!       -------------------------------------------------------------
!       Phisical parameter
!       -------------------------------------------------------------

        BETA  = 0.			! [degree]
        AULVA  = 0.                     ! [%]

!       -------------------------------------------------------------------
!       Compute bottom velocity, and last layer thickness on nodes
!       -------------------------------------------------------------------

	call uvbott(ddl,ukbot,vkbot)

!       -------------------------------------------------------------------
!       Start loop on nodes
!       -------------------------------------------------------------------

        do k = 1,nkn

!         -------------------------------------------------------------------
!         Set input arguments
!         -------------------------------------------------------------------

          lmax = ilhkv(k)			!bottom layer
	  totbed(k) = 0.

          HT   = waveh(k)                       !wave height
          MPER = wavep(k)			!wave mean period
          PPER = wavepp(k)			!wave peak period
          WDIR = waved(k)                       !wave direction
          UW   = waveov(k)                      !wave orbital velocity

	  u = ukbot(k)				!u bottom velocity
	  v = vkbot(k)				!v bottom velocity

          D = 0.
          do l = 1,lmax
           D = D + hdknv(l,k)			!gets total water depth
          end do
	  DL = ddl(k)				!last layer thikness
	  DL = dmin1(D,DL)
	  if (lmax .eq. 1) DL = D
          area = areanode(lmax,k)		!area of node

          call getmd(u,v,UZ,CDIR)		!gets UZ, CDIR
          GD = gskm(k)				!gs50
	  z0b = z0bk(k)
          call zvel(UZ,DL,GD,z0b,Z)		!get Z
          call getts(lmax,k,temp,salt)          !gets temp and salt
          if (temp .eq. 0. .and. salt .eq. 0.) then
            temp = temref
            salt = salref
          end if
          TEM = temp
          SAL = salt

!         -------------------------------------------------------------------
!         Calculate sediment transport rate for each sediment class in node k
!         -------------------------------------------------------------------

          do is = 1,nscls
            scns(is) = scn(lmax,k,is)
          end do

          call sedt05(k,D,DL,UZ,Z,CDIR,HT,MPER,PPER,WDIR,GD,riph(k),
     $ripl(k),BETA,TEM,SAL,bedn(1,1,k),percbd(1,1,k),AULVA,TIMEDR,
     $nscls,gs,UW,hzoff,scns,sedx(1,k),sedy(1,k),ws,gdx(k),gdy(k),
     $lmax,eps,tao(k),Z0,sloads,sflx(1,k))

	  z0b = real(Z0)
	  !z0bk(k) = max(z0b,z0bmin)		!GGUZ0
	  if( z0b < z0bmin ) z0b = z0bmin

          !if( is_external_boundary(k) ) then
          !  sloads = 0.
          !  sflx = 0.
	  !end if

          do is = 1,nscls
	    totbed(k) = totbed(k) + real(sqrt(sedx(is,k)**2 +
     $			sedy(is,k)**2)* bedn(1,3,k))
	    sload(k,is) = sloads(is) * area	!convert to kg/s
          end do

!       -------------------------------------------------------------------
!       End of node loop
!       -------------------------------------------------------------------

        end do

!       -------------------------------------------------------------------
!       Compute settling velocity (either for flocs)
!       -------------------------------------------------------------------

	call set_settling(nscls,ws,scn,wsink)

        end

! ********************************************************************
! SUBROUTINE SEDT05

        subroutine sedt05(k,D,DL,UZ,Z,CDIR,HT,MPER,PPER,WDIR,GD,RHINP,
     $RLINP,BETA,TEM,SAL,BEDCHA,PERCS,AULVA,TIMEDR,nscls,gs,UW,hzoff,
     $scns,sedx,sedy,ws,gdx,gdy,lmax,eps,taok,Z0,sloads,sflx)

	use mod_depth
	use mod_layer_thickness
	use mod_diff_visc_fric
        use levels
        use basin, only : nkn,nel
        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05

        implicit none

! ------------ INPUT VARIABLES -----------------
        DOUBLE PRECISION D        	!WATER DEPTH (M)
        DOUBLE PRECISION DL		!DEPTH OF BOTTOM LAYER (M)
        DOUBLE PRECISION UZ       	!AMBIENT CURRENT AT HEIGHT Z ABOVE THE SEAFLOOR (M/S)
        DOUBLE PRECISION Z        	!HEIGHT OF UZ ABOVE SEAFLOOR
        DOUBLE PRECISION CDIR     	!DIRECTION OF AMBIENT CURRENT (DEGREES TRUE)
        DOUBLE PRECISION HT       	!WAVE HEIGHT (M)
        DOUBLE PRECISION MPER		!MEAN WAVE PERIOD (S)
        DOUBLE PRECISION PPER		!PEAK WAVE PERIOD (S)
        DOUBLE PRECISION PER      	!WAVE PERIOD (S)
        DOUBLE PRECISION WDIR     	!WAVE PROPAGATION DIRECTION (DEGREES TRUE)
        DOUBLE PRECISION GD       	!SEDIMENT GRAIN DIAMETER (M)
        DOUBLE PRECISION RHINP    	!INPUT RIPPLE HEIGHT (M)
        DOUBLE PRECISION RLINP    	!INPUT RIPPLE LENGTH (M)
        DOUBLE PRECISION BETA     	!BED SLOPE (DEGREE)
        DOUBLE PRECISION TEM      	!in degrees Celsius
        DOUBLE PRECISION SAL      	!practical salinity scale (seawater=35)
        DOUBLE PRECISION TIMEDR   	!Duration of call to Sedtrans (s)
        DOUBLE PRECISION PERCS(nlbdim,nscls)   !fraction of sediment [0,1]
        DOUBLE PRECISION BEDCHA(nlbdim,3) ! bed characteristics in 3 column table
                                        ! (1) depth below sediment surface (m)
                                        ! (2) critical erosion stress (Pa)
                                        ! (3) dry bulk density (kg/m**3) user input
        DOUBLE PRECISION UW       	!bottom orbital velocity from WWM
        integer nscls			!number of grainsize classes
	double precision hzoff

! cohesive
        DOUBLE PRECISION AULVA    	!PERCENTAGE OF AREA COVERED BY THE ALGAE 'ULVA' (%)

! ------------ OUTPUT VARIABLES -----------------
        double precision sedx(nscls)
        double precision sedy(nscls)
        double precision sflx(nscls)    !flux of suspend sediment [m3/m2]
        double precision ws(nscls)      !SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
	double precision sloads(nscls)	!suspended sediment load [kg/m²s]
        real taok			!wave-current shear stress

! ------------ LOCAL VARIABLES -----------------
        double precision bmix		!thickness of active layer
        DOUBLE PRECISION UB       	!MAXIMUM WAVE INDUCED ORBITAL VELOCITY AT THE BOTTOM (M/S)
        DOUBLE PRECISION FCW      	!BOTTOM (SKIN) FRICTION FACTOR
        DOUBLE PRECISION UA       	!CURRENT SPEED TO BE USED IN BOTTOM STRESS CALC. (M/SEC)
        DOUBLE PRECISION U100     	!CURRENT SPEED AT 1 M. ABOVE SEABED (M/SEC)
        DOUBLE PRECISION USTCWS   	!COMBINED SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTCW    	!COMBINED TOTAL SHEAR VELOCITY OF GM
        DOUBLE PRECISION Z0       	!BED ROUGHNESS LENGTH (M)
        DOUBLE PRECISION RHOW     	!DENSITY OF FLUID (WATER)  (KG/M**3)
        DOUBLE PRECISION VISC   	!dynamic viscosity of the (sea)water (KG/M*SEC (OR N.S/M**2))
        DOUBLE PRECISION VISK     	!KINEMAMIC VISCOSITY OF THE FLUID (KG/M*SEC (OR N.S/M**2))
        DOUBLE PRECISION RHOS           !DENSITY OF SEDIMENT MINERAL(S)   (KG/M**3)
        DOUBLE PRECISION AB	        !EXCURSION LENGTH OF BOTTOM WAVE ORBIT (M)
        DOUBLE PRECISION WL      	!WAVE LENGTH (M)
        DOUBLE PRECISION Z0C		!APPARENT BED ROUGHNESS LENGTH (M) 
        DOUBLE PRECISION PHIB		!ANGLE BETWEEN WAVE AND CURRENT DIRECTIONS (RADIANS)
        DOUBLE PRECISION PHI100		!ANGLE BETWEEN WAVE AND CURRENT DIRECTIONS AT 1 M. ABOVE SEABED
        DOUBLE PRECISION DELTACW	!HEIGHT OF THE WAVE-CURRENT BOUNDARY LAYER
        DOUBLE PRECISION USTCS		!CURRENT SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTWS		!WAVE SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTCWSB	!TRANSPORT-RELATED COMBINED SHEAR VELOCITY
        DOUBLE PRECISION USTC		!TOTAL CURRENT SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTW		!TOTAL WAVE SHEAR VELOCITY OF GM
        DOUBLE PRECISION RPLCOEF	!RIPPLE COEFFICIENT FOR SHEAR VELOCITY CONVERSION
        DOUBLE PRECISION USTBF		!CRITICAL SHEAR VELOCITY FOR RIPPLE BREAKOFF

        double precision pers(nscls)
        double precision scns(nscls) 	!cohesive suspended sediment conc (kg/m3)
        double precision dxx(nscls)	!DIMENSIONLESS GRAIN SIZE

        double precision usb(nscls)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF BEDLOAD
        double precision uss(nscls)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF SUSPENDED
        double precision ust(nscls)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF SHEET FLOW

        double precision gs(nscls)      !SEDIMENT GRAIN DIAMETER (M)
        double precision FALL           !average fall velocity
        double precision DX             !average dimensionless grain diameter
        double precision USTCRB         !average critical velocity
        double precision USTCRS         !average critical velocity
        double precision USTUP          !average sheet flow velocity
        integer lmax,is,k		!counters
        real gdx,gdy			!slope gradients
        double precision alph           !slope effect
        double precision eps(0:nlvdi,nkn,nscls)	!vertical mixing coefficient
	double precision uslim
	double precision pcoes		!% of  fine sediments
	logical cohes
	double precision btype		!type of active layer or bmix thickness
	integer dtype			!type of diffusivity algoritm
	integer wtype			!type of bottom particle velocity algoritm
        DOUBLE PRECISION wwsb !settling velocity considering biological effect

!       --------------------------------------------------
!       Initialize variables
!       --------------------------------------------------

        do is = 1,nscls
          sedx(is) = 0.d0
          sedy(is) = 0.d0
          pers(is) = PERCS(1,is)
	  sloads(is) = 0.d0
          sflx(is) = 0.d0
        end do

        bmix    = 0.d0
        RHOS    = RHOSED
	USTCRB  = 0.d0
	USTUP   = 0.d0
        PHIB    = 0.d0
	PHI100  = 0.d0
	DELTACW = 0.d0
	USTCS   = 0.d0
	USTWS   = 0.d0
	USTCWSB = 0.d0
	USTC    = 0.d0
	USTW    = 0.d0
        RPLCOEF = 0.d0
	USTBF   = 0.d0
	Z0      = 0.d0
	FCW     = 0.d0
	UA      = 0.d0
	UB      = 0.d0
	U100    = 0.d0
	USTCWS  = 0.d0
	USTCW   = 0.d0
	wwsb    = 1.d0

        if ( D .LT. hzoff ) return

	pcoes = 0.d0
	do is = 1,nbcc
	  pcoes = pcoes + pers(is)
	end do

	cohes = pcoes.gt.KCOES   !true=cohesive; false=non cohesive

!       -------------------------------------------------------------------
!       Calculate slope effect
!       -------------------------------------------------------------------

        call bedslope(cdir,gdx,gdy,angrep,alph)

!       -------------------------------------------------------------------
!       Calculate water density and viscosity
!       -------------------------------------------------------------------

        call densvisc(tem,sal,rhow,visc)	!get density and viscosity
        VISK = visc/rhow 		        !get viscosity

!       -------------------------------------------------------------------
!       Calculate wave-induced bottom particle velocity and orbital diameter
!	following:
!	  -1- Linear wave theory with Hs and Tmean 
!	  -2- Linear wave theory with Hs and Tpeak 
!	  -3- Bottom particle velocity from WWM and Tpeak
!	  -4- Peak wave lenght from WWM (not yet implemented)
!       -------------------------------------------------------------------

	wtype = 1
        if (wtype .eq. 1) then
	  PER = MPER
          CALL OSCIL(HT,PER,D,UB,AB,WL)
	else if (wtype .eq. 2) then
	  PER = PPER
          CALL OSCIL(HT,PER,D,UB,AB,WL)
	else if (wtype .eq. 3) then
	  PER = MPER
	  if (PER .NE. 0.) THEN
            CALL OSCIL(HT,PER,D,UB,AB,WL)
            UB = UW
          end if
	!else if (wtype .eq. 4) then
	  !AB = ABR
	  !WL = PWL
	  !if (WL .ne. 0. .and. PER .ne. 0.) then
	  !  !UB = HT*G*PER/(2.*WL*COSH(2.*PII*D/WL))
	  !  UB = PII*HT/(PER*SINH(2.*PII*D/WL))
	  !else
	  !  UB = 0.d0
	  !end if
	end if

!       -------------------------------------------------------------------
!       Calculate fall velocity, dxx and threshold criteria for d50 (GD)
!       -------------------------------------------------------------------

        do is = 1,nscls
          CALL THRESH(VISK,gs(is),RHOS,RHOW,ws(is),usb(is),uss(is),
     +ust(is),dxx(is))
        end do

	if ( wsetbio ) call wsinkbio(tem,wwsb)
        do is = 1,nbcc
          ws(is) = max(1d-6,ws(is))
 	  WSI(is) = min(ws(is),wwsb)
	end do

        CALL THRESH(VISK,GD,RHOS,RHOW,FALL,USTCRB,USTCRS,USTUP,DX)

!       -------------------------------------------------------------------
!       Correct bedload threshold criteria if function of bed slope
!       -------------------------------------------------------------------

        USTCRB = USTCRB*sqrt(alph)

!       -------------------------------------------------------------------
!       Update threshold criteria in function of mud fraction (Van Rijn, 93)
!       -------------------------------------------------------------------
	
        if ( pcoes.lt. KCOES .and. pcoes.gt.0.01d0 ) then
          USTCRB = USTCRB*(pcoes*100.d0)**0.5
          USTCRS = USTCRS*(pcoes*100.d0)**0.5

	  uslim = sqrt(BEDCHA(1,2)/rhow)
	  if (USTCRB .gt. uslim) USTCRB = uslim
	  if (USTCRS .gt. uslim) USTCRS = uslim
        end if

        USTCRS = dmax1(USTCRB,USTCRS)
	USTUP = dmax1(USTUP,USTCRS)
	!USTUP = 100000.

!       -------------------------------------------------------------------
!       Calculate bottom stresses for d50
!       -------------------------------------------------------------------

        CALL FRICFAC(COHES,D,UZ,Z,CDIR,UB,PER,WDIR,AB,GD,VISK,RHOW,
     $RHOS,Z0C,PHIB,PHI100,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,
     $RPLCOEF,USTBF,USTCRB,USTUP,Z0,FCW,UA,U100,USTCWS,USTCW,RHINP,
     $RLINP)

        !if(ustcw.le.0.d0) ustcw=dmax1(ustc,ustw)
        !if(ustcws.le.0.d0) ustcws=dmax1(ustcs,ustws)	!GGGUUU
	ustcw = max(0.D+0,ustcw)
	ustcws = max(0.D+0,ustcws)
        taok = real(rhow * ustcws**2)

!       -------------------------------------------------------------------
!       Calculate vertical diffusion coefficient following (dtype):
!	  -1- Algebric Van Rijn 1993
!	  -2- Viscosity from GOTM and beta from Van Rijn 1993
!	  -3- Viscosity from GOTM and beta from Villatoro et al 2010
!         -4- Diffusivity from GOTM 
!       -------------------------------------------------------------------

	dtype = 2
        do is = 1,nscls
         call vmixcoef(dtype,lmax,hdknv(1,k),D,ws(is),ustc,ustcw,ub,
     $		       deltacw,ht,per,dxx(is),difv(1,k),visv(1,k),
     $		       eps(1,k,is))
        end do

!       -------------------------------------------------------------------
!       Compute active layer following (btype): 
!	  -1- bmix = roughness height (??)
!	  -2- Zanke: bmix = 0.02 * depth
!	  -3- Harris and Wiberg (1997) as in ROMS
!	  -4- Harris and Wiberg (2001)
!         -5- Erosion depth as in SEDTRANS05
!         -6- Constant = btype
!       -------------------------------------------------------------------

	btype = 0.002d0
	btype = 3
	call activelayer(btype,D,Z0,USTCWS,USTCRB,GD,pcoes,PER,
     $RHOW,RHOS,RHINP,RLINP,BEDCHA,TIMEDR,bmix)

!       -------------------------------------------------------------------
!       Set the surficial layer equal to the active layer (bmix)
!       -------------------------------------------------------------------

        !call bedact(nscls,bmix,BEDCHA,PERCS)
        !do is = 1,nscls
        !  pers(is) = PERCS(1,is)
        !end do

!       -------------------------------------------------------------------
!       Calculate sediment transport rate for each sediment class in node k
!       -------------------------------------------------------------------

        if( cohes ) then

!         -------------------------------------------------------------------
!         COHESIVE SEDIMENT
!         -------------------------------------------------------------------

          call cohse(BEDCHA,TIMEDR,USTCWS,USTCW,DL,RHOW,
     $AULVA,UB,CDIR,Z0,nscls,scns,ws,usb,uss,ust,gs,dxx,
     $Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,USTBF,
     $VISK,pers,RHOS,sloads,sflx)

        else

!         -------------------------------------------------------------------
!         NON-COHESIVE SEDIMENT
!         -------------------------------------------------------------------

          call nonco(BEDCHA,TIMEDR,D,DL,UA,UB,U100,HT,PER,CDIR,
     $RHOW,USTCWS,USTCW,usb,uss,ust,Z0,BETA,RHOS,FCW,WDIR,
     $Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,USTBF,
     $PHIB,PHI100,nscls,pers,gs,dxx,ws,scns,sedx,sedy,sloads,sflx)

        end if

        end

! ********************************************************************
! SUBROUTINE ACTIVELAYER
! This subroutine computes the bed active layer following 
!	 -1- bmix = roughness height (??)
!	 -2- Zanke: bmix = 0.02 * depth
!	 -3- Harris and Wiberg (1997) as in ROMS
!	 -4- Harris and Wiberg (2001)
!        -5- Erosion depth as in SEDTRANS06

        subroutine activelayer(btype,D,Z0,USTCWS,USTCRB,GD,pcoes,
     $PER,RHOW,RHOS,RHEIGHT,RLENGTH,BEDCHA,TIMEDR,bmix)

        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05

	implicit none

! ------------ INPUT VARIABLES -----------------
	double precision btype		!type of active layer of bmix thickness
        DOUBLE PRECISION D        	!WATER DEPTH (M)
        DOUBLE PRECISION Z0       	!BED ROUGHNESS LENGTH (M)
        DOUBLE PRECISION USTCWS   	!COMBINED SKIN-FRICTION SHEAR VELOCITY OF GM
        double precision USTCRB    	!CRITICAL SHEAR VELOCITY FOR INITIATION OF BEDLOAD
        double precision GD       	!SEDIMENT GRAIN DIAMETER (M)
        double precision pcoes		!% of fine sedimens
        double precision PER      	!WAVE PERIOD (S)
        double precision RHOW  	   	!DENSITY OF FLUID (WATER)  (KG/M**3)
        double precision RHOS   	!DENSITY OF SEDIMENT MINERAL(S)   (KG/M**3)
        double precision RHEIGHT  	!PREDICTED RIPPLE HEIGHT
        double precision RLENGTH  	!PREDICTED RIPPLE LENGTH
        DOUBLE PRECISION BEDCHA(nlbdim,3) ! bed characteristics in 3 column table
                                        ! (1) depth below sediment surface (m)
                                        ! (2) critical erosion stress (Pa)
                                        ! (3) dry bulk density (kg/m**3) user input
        DOUBLE PRECISION TIMEDR   	!Duration of call to Sedtrans (s)

! ------------ OUTPUT VARIABLES -----------------
        double precision bmix		!thickness of active layer

! ------------ LOCAL VARIABLES -----------------
        DOUBLE PRECISION TAU0		!effective skin friction shear stress (Pa)
        double precision TAUCE		!critical erosion stress at depth ZS (Pa)
        double precision RHOB      	!sediment bulk density
        double precision QB		!bedload trasnport
        double precision S        	!RHOS/RHOW
	double precision tsc		!time scale for current tranport
	double precision ts		!basic time scale
        double precision bmixs		!thickness of sandy active layer
        double precision bmixm		!thickness of silty active layer
        double precision CB		!bed concentration of sediments (1-POR)
        double precision psand  	!fraction of sandy sediments
        double precision taucr		!CRITICAL SHEAR STRESS FOR INITIATION OF BEDLOAD
	double precision aux
	double precision bmax
        DOUBLE PRECISION MAXEMASS	!Maximum erodable sediment mass to avoid problem with drag reduction
        DOUBLE PRECISION EMASS		!eroded mass (kg/m2)
        INTEGER NBED         		!number of rows in BEDCHA that are used
	integer i
	integer imode
	logical blimit

!       -------------------------------------------------------------------
!       Inititalizes valiables
!       -------------------------------------------------------------------

	imode = int(btype)
	blimit = .false.
	bmixs = 0.d0
	bmixm = 0.d0
	qb = 0.d0

	tauce = BEDCHA(1,2)
	rhob = BEDCHA(1,3)
	s = rhos/rhow
	cb = rhob/rhos
        tau0=rhow*ustcws**2
        psand = 1.d0 - pcoes

!       -------------------------------------------------------------------
!       Compute active layer
!       -------------------------------------------------------------------

	if (imode .eq. 1) then
!         -------------------------------------------------------------------
!         Roughness height
!         -------------------------------------------------------------------
	  blimit = .true.
          bmix = Z0 * 30.d0

	elseif (imode .eq. 2) then
!         -------------------------------------------------------------------
!         Zanke, 0.02 * depth
!         -------------------------------------------------------------------
	  blimit = .true.
	  bmix = 0.02d0 * D

	else if (imode .eq. 3) then
!         -------------------------------------------------------------------
!         Harris and Wiberg (1997) as in ROMS
!         -------------------------------------------------------------------
	  blimit = .true.
	  bmix = 0.007d0*(USTCWS**2 - USTCRB**2)

	else if (imode .eq. 4) then
!         -------------------------------------------------------------------
!         Harris and Wiberg (2001)
!         -------------------------------------------------------------------
	  blimit = .true.
!         -------------------------------------------------------------------
!         Compute bedload transport following Meyer-Peter and Muller 1948 (cm3/cm2s)
!         -------------------------------------------------------------------
 
          aux = 8.0d0*((rhos-rhow)*g)**(-1)*10.d0
          taucr = rhow*USTCRB**2
          if(tau0.gt.taucr) then
	    qb = aux*((tau0-taucr)*10.d0)**1.5
	  endif 

!         -------------------------------------------------------------------
!         Compute basic time scales for current. In case of waves takes 
!	  the maximun between tsc and PER
!         -------------------------------------------------------------------

          tsc = rheight*rlength/(g*(s-1.d0)*gd**3)**0.5
          ts = dmax1(tsc,per)
 
!         -------------------------------------------------------------------
!         Compute sandy active layer
!         -------------------------------------------------------------------

	  if(rlength.gt.0.d0) then
	    bmixs = qb*ts/(2.d0*cb*rlength*100.d0)
	    bmixs = bmixs / 100.d0	!from cm to m
	    bmixs = dmin1(bmixs,rheight/2.d0)
	  end if

!         -------------------------------------------------------------------
!         Compute silty active layer
!         -------------------------------------------------------------------

	  bmixm = 8.0d0*(tau0 - tauce)*0.1d0		!from dyn/cm2 to pascal

!         -------------------------------------------------------------------
!         Compute total active layer (sandy + silty)
!         -------------------------------------------------------------------
	
	  bmix = bmixs*psand + bmixm*pcoes

	else if (imode .eq. 5) then
!         -------------------------------------------------------------------
!         Erosion depth as in SEDTRANS06
!         -------------------------------------------------------------------

	  blimit = .true.
          if( tau0.GT.BEDCHA(1,2) .OR. ( tau0.EQ.BEDCHA(1,2) .AND. 
     $	      tau0.GE.BEDCHA(2,2) ) ) THEN

            DO 10, I=2,nlbdim
             IF (BEDCHA(I,1).EQ.0.0d0) GOTO 20
10          CONTINUE
20          NBED=I-1

            CALL EROS1(NBED,nlbdim,BEDCHA,tau0,MAXEMASS,EMASS,bmix,
     $	      TIMEDR,TAUCE)

	  end if 

	else 
!         -------------------------------------------------------------------
!         Set constant acrive layer thickness
!         -------------------------------------------------------------------

	  bmix = btype

	end if

!       -------------------------------------------------------------------
!       Set min (6 * GD) and max ((BEDCHA(3,1) + BEDCHA(4,1))/2.d0) active layer
!       -------------------------------------------------------------------

	bmax = (BEDCHA(3,1) + BEDCHA(4,1))/2.d0
	if (bmix.gt.bmax)write(*,*)'BMIX troppo grande',bmix,bmax

	if (blimit) then
	  bmix = dmin1(bmix,bmax)	!???????	
	  bmix = dmax1(bmix,10d-4)
	end if

	end

! ********************************************************************
! SUBROUTINE COHSE
! This subroutine computes the sediment transport for cohesive
! sediments

        subroutine cohse(BEDCHA,TIMEDR,USTCWS,USTCW,DL,RHOW,
     $AULVA,UB,CDIR,Z0,nscls,scns,ws,usb,uss,ust,gs,dxx,
     $Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,USTBF,
     $VISK,pers,RHOS,sloads,sflux)

        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05

        implicit none

! ------------ INPUT VARIABLES -----------------
        DOUBLE PRECISION BEDCHA(nlbdim,3) ! bed characteristics in 3 column table
                                  	! (1) depth below sediment surface (m)
                                  	! (2) critical erosion stress (Pa)
                                  	! (3) dry bulk density (kg/m**3)          user input
        DOUBLE PRECISION TIMEDR   	!Duration of call to Sedtrans (s)
        DOUBLE PRECISION USTCWS   	!Combined skin friction stress            FRICFRAC
        DOUBLE PRECISION DL		!DEPTH OF BOTTOM LAYER (M)
        DOUBLE PRECISION RHOW     	!DENSITY OF FLUID (WATER)  (KG/M**3)
        DOUBLE PRECISION AULVA    	!PERCENTAGE OF AREA COVERED BY THE ALGAE 'ULVA' (%)
        DOUBLE PRECISION UB       	!MAXIMUM WAVE INDUCED ORBITAL VELOCITY AT THE BOTTOM (M/S)
        DOUBLE PRECISION USTCW    	!COMBINED TOTAL SHEAR VELOCITY OF GM
        DOUBLE PRECISION CDIR     	!DIRECTION OF AMBIENT CURRENT (DEGREES TRUE)
        DOUBLE PRECISION Z0       	!BED ROUGHNESS LENGTH (M)

        integer nscls			!number of grainsize classes
        double precision scns(nscls) 	!cohesive suspended sediment conc (kg/m3)
        double precision pers(nscls)	!fraction of each class
        DOUBLE PRECISION GD		!SEDIMENT GRAIN DIAMETER (M) 
        double precision gs(nscls)     !SEDIMENT GRAIN DIAMETER (M)
        double precision ws(nscls)     !SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        double precision usb(nscls)    !CRITICAL SHEAR VELOCITY FOR INITIATION OF BEDLOAD
        double precision uss(nscls)    !CRITICAL SHEAR VELOCITY FOR INITIATION OF SUSPENDED
        double precision ust(nscls)    !CRITICAL SHEAR VELOCITY FOR INITIATION OF SHEET FLOW
        double precision dxx(nscls)    !DIMENSIONLESS GRAIN SIZE

        DOUBLE PRECISION Z0C      	!APPARENT BED ROUGHNESS LENGTH (M)
        DOUBLE PRECISION DELTACW  	!HEIGHT OF THE WAVE-CURRENT BOUNDARY LAYER
        DOUBLE PRECISION USTCS    	!CURRENT SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTWS    	!WAVE SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTCWSB  	!TRANSPORT-RELATED COMBINED SHEAR VELOCITY
        DOUBLE PRECISION USTC     	!TOTAL CURRENT SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTW     	!TOTAL WAVE SHEAR VELOCITY OF GM
        DOUBLE PRECISION RPLCOEF  	!RIPPLE COEFFICIENT FOR SHEAR VELOCITY CONVERSION
        DOUBLE PRECISION USTBF   	 !CRITICAL SHEAR VELOCITY FOR RIPPLE BREAKOFF
        DOUBLE PRECISION VISK     	!KINEMAMIC VISCOSITY OF THE FLUID (KG/M*SEC (OR N.S/M**2))

! ------------ OUTPUT VARIABLES -----------------
	double precision sloads(nscls)	  	!suspended sediment load [kg/m2s]
        double precision sflux(nscls)  !flux of suspend sediment [m3/m2]

! ------------ LOCAL VARIABLES -----------------
        DOUBLE PRECISION TAU0		!effective skin friction shear stress (Pa)
        DOUBLE PRECISION CONC(NBCC) 	!COHESIVE SUSPENDED SEDIMENT CONCENTRATION (kg/m3)
        INTEGER NBED         		!number of rows in BEDCHA that are used
        DOUBLE PRECISION ZS       	!Height change in bed surface (m)
                                  	!(positive: erosion, negative: deposition)
        DOUBLE PRECISION TAOS		!solid transmitted stress due to Ulva (Pa)
        DOUBLE PRECISION TCONC1		!total SSC (sum of CONC) before deposition (kg/m**3)
        DOUBLE PRECISION TCONC2  	!total SSC (sum of CONC) after deposition (kg/m**3)
        DOUBLE PRECISION RHOS		!dry bulk density at depth ZS (kg/m**3)
        DOUBLE PRECISION TAUCE		!critical erosion stress at depth ZS (Pa)
        DOUBLE PRECISION EMASS		!eroded mass (kg/m2)
        DOUBLE PRECISION MAXEMASS	!Maximum erodable sediment mass to avoid problem with drag reduction
        DOUBLE PRECISION C0     	!REFERENCE CONCENTRATION AT Z0 (KG/M^3)
        DOUBLE PRECISION C0A		!DEPTH AVERAGED REFERENCE CONCENTRATION AT Z0 (KG/M^3)
        DOUBLE PRECISION C0AI		!DEPTH AVERAGED REFERENCE CONCENTRATION AT Z0 (KG/M^3) per is
        DOUBLE PRECISION QS		!SUSPENDED SEDIMENT TRANSPORT RATE (KG/M/S)
        DOUBLE PRECISION QSDIR		!DIRECTION OF SUSPENDED SEDIMENT TRANSPORT (DEGREE)
        DOUBLE PRECISION DMASS		!deposited mass (kg/m2)
	DOUBLE PRECISION RHOA
        DOUBLE PRECISION r, cdep
        integer i,is			!counters
       
!       --------------------------------------------------
!       Initialize variables
!       --------------------------------------------------

        emass  = 0.d0
	ZS     = 0.d0
	sloads = 0.d0
        sflux  = 0.d0

!       -------------------------------------------------------------------
!       Compute the total SSC before deposition and erosion
!       -------------------------------------------------------------------

        TCONC1 = 0.d0
        do I=1,NBCC
          CONC(I) = scns(i)
          TCONC1 = TCONC1 + CONC(I)
	end do

!       -------------------------------------------------------------------
!       Compute TAU0, the effective skin friction stress acting on the bed
!       drag reduction and solid transmitted stress due to ulva
!       -------------------------------------------------------------------

        TAU0=RHOW*USTCWS**2

        CALL DRAGRED(TAU0,TCONC1,DL,MAXEMASS)
  
        CALL SOLULVA(TAU0,AULVA,TAOS)

!       -------------------------------------------------------------------
!       EROSION: Test if effective bed shear stress (TAU0) is larger or 
!                equal to surface critical erosion stress (BEDCHA(1,2)
!       -------------------------------------------------------------------

        if( TAU0.GT.BEDCHA(1,2) .OR. ( TAU0.EQ.BEDCHA(1,2).AND.
     $      TAU0.GE.BEDCHA(2,2) ) ) THEN

          DO 10, I=2,nlbdim
           IF (BEDCHA(I,1).EQ.0.0d0) GOTO 20
10        CONTINUE
20        NBED=I-1

          CALL EROS1(NBED,nlbdim,BEDCHA,TAU0,MAXEMASS,EMASS,ZS,TIMEDR,
     $TAUCE)

!         -------------------------------------------------------------------
!         Limit EROSION to sediment in the active layer
!         -------------------------------------------------------------------

          if(ZS .gt. BEDCHA(2,1)) THEN 
	    ZS = BEDCHA(2,1)
	    RHOA = (BEDCHA(1,3)+BEDCHA(2,3)) / 2.d0
	    EMASS = RHOA * ZS
	  end if

          do is = 1,nscls
	    sloads(is) = emass * pers(is) / TIMEDR
	    sflux(is) = -ZS * pers(is)
          end do
	
        end if

!       -------------------------------------------------------------------
!       DEPOSITION:
!       -------------------------------------------------------------------        

!       -------------------------------------------------------------------
!       Deposition for cohesive suspended sediments
!       -------------------------------------------------------------------
 
        if(TCONC1.gt.0.d0) then 
          call depos(DL,TIMEDR,NBCC,CONC,TAU0,RHOW,TCONC1,TCONC2,VISK,
     $DMASS)
	  
	  !RHOA = 2.d0/(RHOMUD+BEDCHA(1,3))
	  RHOA = RHOMUD
          do is = 1,nbcc
	    sloads(is) = sloads(is) - dmass * pers(is) / TIMEDR
	    sflux(is) = sflux(is) + DMASS*pers(is)/RHOA
	  end do
        end if

!       -------------------------------------------------------------------
!       Deposition for non-cohesive suspended sediments (if present)
!       -------------------------------------------------------------------

        do is = nbcc+1,nscls
          CALL PROFL(DL,RHOW,ws(is),UB,CDIR,USTCWS,USTCW,usb(is),
     $uss(is),ust(is),Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,
     $RPLCOEF,USTBF,Z0,gs(is),dxx(is),RHOS,QS,QSDIR,C0,C0A)

          C0AI = C0A * pers(is)
  	  C0AI = MIN(C0AI,scns(is))

          r = ws(is)/DL
          cdep = scns(is) * exp(-r*TIMEDR)
          C0AI = max(cdep, C0AI)
          sloads(is) = (C0AI - scns(is))*DL/TIMEDR
          RHOA = RHOSED * (1.d0 - SURFPOR)
          sflux(is) = sflux(is) - sloads(is) * timedr / RHOA
        end do

        end 

! ********************************************************************
! SUBROUTINE NONCO
! This subroutine computes the sediment transport for non-cohesive
! sediments

        subroutine nonco(BEDCHA,TIMEDR,D,DL,UA,UB,U100,HT,PER,CDIR,
     $RHOW,USTCWS,USTCW,usb,uss,ust,Z0,BETA,RHOS,FCW,WDIR,
     $Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,USTBF,
     $PHIB,PHI100,nscls,pers,gs,dxx,ws,scns,sedx,sedy,sloads,sflux)

        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05

        implicit none

! ------------ INPUT VARIABLES -----------------
        DOUBLE PRECISION BEDCHA(nlbdim,3) ! bed characteristics in 3 column table
                                        ! (1) depth below sediment surface (m)
                                        ! (2) critical erosion stress (Pa)
                                        ! (3) dry bulk density (kg/m**3)          user input
        DOUBLE PRECISION TIMEDR   	!Duration of call to Sedtrans (s)
        DOUBLE PRECISION D        !WATER DEPTH (M)
        DOUBLE PRECISION DL       !DEPTH OF BOTTOM LAYER (M)
        DOUBLE PRECISION UA       !CURRENT SPEED TO BE USED IN BOTTOM STRESS CALC. (M/SEC)
        DOUBLE PRECISION UB       !MAXIMUM WAVE INDUCED ORBITAL VELOCITY AT THE BOTTOM (M/S)
        DOUBLE PRECISION U100     !CURRENT SPEED AT 1 M. ABOVE SEABED (M/SEC)
        DOUBLE PRECISION HT       !WAVE HEIGHT (M)
        DOUBLE PRECISION PER      !WAVE PERIOD (S)
        DOUBLE PRECISION CDIR     !DIRECTION OF AMBIENT CURRENT (DEGREES TRUE)
        DOUBLE PRECISION RHOW     !DENSITY OF FLUID (WATER)  (KG/M**3)
        DOUBLE PRECISION USTCWS   !COMBINED SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTCW    !COMBINED TOTAL SHEAR VELOCITY OF GM
        DOUBLE PRECISION Z0       !BED ROUGHNESS LENGTH (M)
        DOUBLE PRECISION BETA     !BED SLOPE (DEGREE)
        DOUBLE PRECISION RHOS     !DENSITY OF SEDIMENT MINERAL(S)   (KG/M**3)
        DOUBLE PRECISION FCW      !BOTTOM (SKIN) FRICTION FACTOR
        DOUBLE PRECISION WDIR     !WAVE PROPAGATION DIRECTION (DEGREES TRUE)
        DOUBLE PRECISION Z0C      !APPARENT BED ROUGHNESS LENGTH (M)
        DOUBLE PRECISION DELTACW  !HEIGHT OF THE WAVE-CURRENT BOUNDARY LAYER
        DOUBLE PRECISION USTCS    !CURRENT SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTWS    !WAVE SKIN-FRICTION SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTCWSB  !TRANSPORT-RELATED COMBINED SHEAR VELOCITY
        DOUBLE PRECISION USTC     !TOTAL CURRENT SHEAR VELOCITY OF GM
        DOUBLE PRECISION USTW     !TOTAL WAVE SHEAR VELOCITY OF GM
        DOUBLE PRECISION RPLCOEF  !RIPPLE COEFFICIENT FOR SHEAR VELOCITY CONVERSION
        DOUBLE PRECISION USTBF    !CRITICAL SHEAR VELOCITY FOR RIPPLE BREAKOFF
        DOUBLE PRECISION TB1      !TIME AT WHICH BEDLOAD TRANSPORT CEASES (SEC) 
        DOUBLE PRECISION TB2      !TIME AT WHICH BEDLOAD TRANSPORT RECOMMENCES (SEC)
        DOUBLE PRECISION TS1      !TIME AT WHICH SUSPENDED LOAD TRANSPORT CEASES (SEC)
        DOUBLE PRECISION PERBED   !PERCENTAGE OF TIME SPENT IN ONLY BEDLOAD TRANSPORT PHASE
        DOUBLE PRECISION PERSUSP  !PERCENTAGE OF TIME SPENT IN SUSPENDED LOAD TRANSPORT PHASE
        DOUBLE PRECISION PHIB     !ANGLE BETWEEN WAVE AND CURRENT DIRECTIONS (RADIANS)
        DOUBLE PRECISION PHI100   !ANGLE BETWEEN WAVE AND CURRENT DIRECTIONS AT 1 M. ABOVE SEABED

        integer nscls			!number of grainsize classes
        double precision pers(nscls)	!fraction of sediment [0,1]
        double precision gs(nscls)     !SEDIMENT GRAIN DIAMETER (M)
        double precision ws(nscls)     !SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        double precision scns(nscls) 	!non-cohesive suspended sediment conc (kg/m3)
        double precision usb(nscls)    !CRITICAL SHEAR VELOCITY FOR INITIATION OF BEDLOAD
        double precision uss(nscls)    !CRITICAL SHEAR VELOCITY FOR INITIATION OF SUSPENDED
        double precision ust(nscls)    !CRITICAL SHEAR VELOCITY FOR INITIATION OF SHEET FLOW
        double precision dxx(nscls)    !DIMENSIONLESS GRAIN SIZE


! ------------ OUTPUT VARIABLES -----------------
        double precision sedx(nscls)	!x bedload component [kg/ms]
        double precision sedy(nscls)	!y bedload component [kg/ms]
        double precision sflux(nscls)  !flux of suspend sediment [m3/m2]
	double precision sloads(nscls)		!suspended sediment load [kg/m2s]

! ------------ LOCAL VARIABLES -----------------
        DOUBLE PRECISION QS       !SUSPENDED SEDIMENT TRANSPORT RATE (KG/M/S)
        DOUBLE PRECISION QSDIR    !DIRECTION OF SUSPENDED SEDIMENT TRANSPORT (DEGREE)
        DOUBLE PRECISION SED      !TIME-AVERAGED NET SEDIMENT TRANSPORT AS VOLUME (M**3/S/M)
        DOUBLE PRECISION SEDM     !TIME-AVERAGED NET SEDIMENT TRANSPORT AS MASS (KG/S/M)
        DOUBLE PRECISION SEDDIR   !DIRECTION OF NET SEDIMENT TRANSPORT (AZIMUTH,DEGREES)
        DOUBLE PRECISION C0       !REFERENCE CONCENTRATION AT Z0 (KG/M^3)
        DOUBLE PRECISION C0A  	  !DEPTH AVERAGED REFERENCE CONCENTRATION AT Z0 (KG/M^3)
        DOUBLE PRECISION C0AI  	  !DEPTH AVERAGED REFERENCE CONCENTRATION AT Z0 (KG/M^3) per is
	DOUBLE PRECISION r, cdep
        DOUBLE PRECISION smax
	DOUBLE PRECISION RHOA
        integer is		  !counters

!       --------------------------------------------------
!       Initialize variables
!       --------------------------------------------------
        sedx   = 0.d0
        sedy   = 0.d0
        sflux  = 0.d0
	sloads = 0.d0

!       --------------------------------------------------
!       Loop over grainsize classes for bedload (non-cohesive only)
!       --------------------------------------------------
        do is = nbcc+1,nscls

!        -------------------------------------------------------------------
!        Calculate the duration of the different sediment transport phases
!        -------------------------------------------------------------------
         CALL TIMING(RHOW,UA,UB,PER,U100,usb(is),uss(is),USTCWS,
     $PHIB,USTCS,USTWS,RPLCOEF,TB1,TB2,TS1,PERBED,PERSUSP)

!        -------------------------------------------------------------------
!        Calculate sediment transport rate and direction
!        -------------------------------------------------------------------
         CALL TRANSPO(D,UA,UB,U100,PER,gs(is),BETA,RHOS,RHOW,
     $usb(is),ust(is),PHIB,PHI100,USTCS,USTWS,USTCWSB,RPLCOEF,TB1,
     $TB2,TS1,PERBED,PERSUSP,USTBF,FCW,USTCWS,HT,CDIR,WDIR,IOPT,
     $dxx(is),ws(is),SED,SEDM,SEDDIR)

         if(SEDDIR.le.90.d0.and.SEDDIR.ge.0.d0)   SEDDIR=90.d0-SEDDIR
         if(SEDDIR.le.360.d0.and.SEDDIR.ge.90.d0) SEDDIR=450.d0-SEDDIR
         SEDDIR = SEDDIR / (45.d0 / atan (1.d0))             !rad
         !SEDM = SEDM / bedcha(1,3)

         sedx(is) = SEDM * pers(is) * cos(SEDDIR)
         sedy(is) = SEDM * pers(is) * sin(SEDDIR)

	end do

	if (IOPT.eq.1 .or. IOPT.eq.3 ) return

!       --------------------------------------------------
!       Loop over grainsize classes for suspension
!       --------------------------------------------------

        do is = 1,nscls

	 C0A = 0.d0

!        -------------------------------------------------------------------
!        Calculate velocity profile, suspended sediment concentration profile
!        -------------------------------------------------------------------

         CALL PROFL(DL,RHOW,ws(is),UB,CDIR,USTCWS,USTCW,usb(is),
     $uss(is),ust(is),Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,
     $RPLCOEF,USTBF,Z0,gs(is),dxx(is),RHOS,QS,QSDIR,C0,C0A)

         C0AI = C0A * pers(is)

!        -------------------------------------------------------------------
!        Limit sloads and set sediment bulk density
!        -------------------------------------------------------------------
         if (scns(is) .gt. C0AI) then		!deposition
            r = ws(is)/DL
            cdep = scns(is) * exp(-r*TIMEDR)
            C0AI = max(cdep, C0AI)
            sloads(is) = (C0AI - scns(is))*DL/TIMEDR
            if (is .le. nbcc) then
              RHOA = RHOMUD
            else
              RHOA = RHOSED * (1.d0 - SURFPOR)
            end if
         else					!erosion
           sloads(is) = (C0AI - scns(is))*DL/TIMEDR
           smax = pers(is)*BEDCHA(1,3)*BEDCHA(2,1)/TIMEDR
           sloads(is) = min(sloads(is),smax)
           RHOA = BEDCHA(1,3)
         end if

!        -------------------------------------------------------------------
!        Compute flux [m³/m²]
!        -------------------------------------------------------------------
	 sflux(is) =  -sloads(is) * timedr / RHOA

        end do

        end

! ********************************************************************
! SUBROUTINE BEDLOAD
! This subroutine computes the bedload sediment transport using the 
! sediment continuity equation

        subroutine bedload(nscls,sedx,sedy,dt,bedn,bh,bflx,bload)

        use mod_bound_geom
        use evgeom
        use basin, only : nkn,nel,nen3v
        use mod_sediment

        implicit none

! --- input variables
        integer nscls				!number grainsize classes
        double precision sedx(nscls,nkn)	!bedload transport in x direction [kg/ms]
        double precision sedy(nscls,nkn)	!bedload transport in y direction [kg/ms]
	double precision dt			!time step
        double precision bh(nkn)	    	!bottom height variation [m]
        double precision bedn(nlbdim,3,nkn)     !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)

! --- output variables
        double precision bflx(nscls,nkn)	!flux of bedload sediment [m3/m2]
        double precision bload(nkn)		!flux of bedload sediment [kg]

! --- local variables
        integer k,ie,ii,is
	double precision v2v(nkn)
        double precision bflux
        double precision sexe(nscls)           	!element transport in x direction
        double precision seye(nscls)           	!element transport in x direction
        double precision b,c   		        !x and y derivated form function [1/m]
        double precision area			!area of element/3
        double precision ss(nscls)
	double precision sm(nscls)

!       -------------------------------------------------------------------
!       Initialize local variables
!       -------------------------------------------------------------------

        bflx  = 0.d0
        bload = 0.d0
	v2v   = 0.d0

        do ie = 1,nel

          area = 4. * ev(10,ie)

          do is = 1,nscls
            sm(is) = 0.d0
            ss(is) = 0.d0
          end do

!         -------------------------------------------------------------------
!         Converts node transport value to element value
!         -------------------------------------------------------------------

          do ii=1,3
            k = nen3v(ii,ie)
	    v2v(k) = v2v(k) + area
            do is = 1,nscls
              sm(is) = sm(is) + sedx(is,k)
              ss(is) = ss(is) + sedy(is,k)
            end do
          end do

          do is = 1,nscls
            sexe(is) = sm(is) / 3.
            seye(is) = ss(is) / 3.
          end do

!         -------------------------------------------------------------------
!         Loop over elements vertex
!         -------------------------------------------------------------------

          do ii = 1,3
            k = nen3v(ii,ie)
            b = ev(3+ii,ie)                                   !1/m
            c = ev(6+ii,ie)                                   !1/m

            do is = 1,nscls
              bflux = dt*area * (b*sexe(is) + c*seye(is))
              bflx(is,k) = bflx(is,k) + bflux
            end do
          end do

        end do

!       -------------------------------------------------------------------
!       Loop over area node
!       -------------------------------------------------------------------

        do k = 1,nkn
          do is = 1,nscls
            bflx(is,k) = bflx(is,k) / v2v(k)
	    if( is_external_boundary(k) ) bflx(is,k) = 0.d0
            bload(k) = bload(k) + bflx(is,k)
	    bflx(is,k) = bflx(is,k) / bedn(1,3,k)
          end do
        end do


!       -------------------------------------------------------------------
!       Apply slope limiter
!       -------------------------------------------------------------------

        !call slope_lim_bh(nscls,bflx,bh)

        end

! ********************************************************************

	subroutine slope_lim_bh(nscls,bflx,bh)

	use basin, only : nkn

	implicit none

        integer nscls				!number of grainsize class
        double precision bflx(nscls,nkn)	!bedload sediment contribution [m3/m2]
        double precision bh(nkn)		!bed elevation

	integer k,is
	double precision bbhh(nkn)
	double precision baux(nkn)
	double precision frac(nscls,nkn)
	double precision bbk

!       -------------------------------------------------------------------
!       Sum contributions
!       -------------------------------------------------------------------

	do k = 1,nkn
          bbhh(k) = bh(k)
        end do

	call slopelim(bbhh)

	do k = 1,nkn
	  bbk = 0.d0
	  do is = 1,nscls
	    frac(is,k) = 0.d0
	    bbk = bbk + bflx(is,k)
	  end do
          baux(k) = bbhh(k) + bbk
          do is = 1,nscls
            if(bbk.ne.0.d0) frac(is,k) = dabs(bflx(is,k)/bbk)
          end do
        end do

!       -------------------------------------------------------------------
!       Apply slope limiter
!       -------------------------------------------------------------------

	call slopelim(baux)

!       -------------------------------------------------------------------
!       Correct bflx
!       -------------------------------------------------------------------

	do k = 1,nkn
	  do is = 1,nscls
            bflx(is,k) = (baux(k) - bbhh(k))*frac(is,k)
	  end do
        end do

	end

!******************************************************************

	subroutine slopelim(bbk)

	use evgeom
        use basin, only : nkn,nel,nen3v

        implicit none

        double precision v1v(nkn),v2v(nkn),v3v(nkn),v4v(nkn)
	double precision bbk(nkn),bbe(nel)

        double precision high,hmed,h
	double precision ao
        integer ie,k,ii

        high = 1.d30
	bbe = 0.

        do ie=1,nel
          h = 0.d0
          do ii=1,3
            k = nen3v(ii,ie)
            h = h + bbk(k)
          end do
          bbe(ie) = h/3.d0
        end do

        do k=1,nkn
          v1v(k) = high
          v2v(k) = -high
          v3v(k) = 0.d0
          v4v(k) = 0.d0
        end do

        do ie=1,nel
          ao = ev(10,ie)
          do ii=1,3
            k = nen3v(ii,ie)
            v1v(k) = dmin1(v1v(k),bbe(ie))
            v2v(k) = dmax1(v2v(k),bbe(ie))
            v3v(k) = v3v(k) + ao * bbe(ie)
            v4v(k) = v4v(k) + ao
          end do
        end do

        do k=1,nkn
          hmed = v3v(k) / v4v(k)
          if( bbk(k) .lt. v1v(k) .or. bbk(k) .gt. v2v(k) ) then
            !write(6,*) 'slope: ',k,bbk(k),v1v(k),v2v(k),hmed
            bbk(k) = hmed
          end if
        end do

	end

!******************************************************************
! SUBROUTINE BEDSLOPE
! Compute the bed slope on nodes in the direction of the transport (current)
! and slope effect (alph). When bed slope equal 0 ---> alph = 1

        subroutine bedslope(cdir,gdx,gdy,beta,alph)

        implicit none

! --- input variables
        double precision cdir           !current direction [degree]
        real gdx,gdy			!depth gradient in x and y

! --- output variables
        double precision alph           !slope effect

! --- local variables
        double precision ksl            !node slope angle [radian]
        double precision asp	       	!node upslope aspect angle [radian]
        double precision phi            !angle between asp and slope [radian]
        double precision beta           !angle of repose [radian]
	double precision bmin		!min
        double precision rad		!convert degree to rad

!       -------------------------------------------------------------------
!       Initialize variables
!       -------------------------------------------------------------------

        alph = 1.d0
        rad = atan (1.) / 45.d0
        bmin = 90.d0 * rad 

!       -------------------------------------------------------------------
!       Compute slope angle and slope aspect
!       -------------------------------------------------------------------

	call getmd(-gdx,-gdy,ksl,asp)

        ksl = atan(ksl)

        phi = dabs(asp - cdir) * rad

!       -------------------------------------------------------------------
!       Compute bed slope effect (Soulsby,97, eq.80a)
!       -------------------------------------------------------------------

        if(ksl.ge.beta) then	!avalanche occur
          if (phi.gt.bmin) then
           alph = 1.7d0
          else
           alph = 0.2d0
          end if
          return
        end if

        alph = (cos(phi)*sin(ksl) + sqrt((cos(ksl))**2.*(tan(beta))**2.
     $ - (sin(phi))**2.*(sin(ksl))**2.)) / tan(beta)

        alph = dmax1(alph,0.2d0)

        end

! ********************************************************************
! SUBROUTINE BEDMAN
! This subroutine manages the bed composition in function of what is eroded
! or deposited

        subroutine bedman(gs,nscls,timedr,bflux,sflux,bh,gskm,bdh,
     $	percbd,bedn)

	use basin, only : nkn
        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05
	use mod_debug

        implicit none

        integer nscls				!number of grainsize class
        double precision gs(nscls)		!grainsize class
        double precision timedr                	!time step [s]
        double precision bflux(nscls,nkn)	!bedload sediment contribution [m3/m2]
        double precision sflux(nscls,nkn)	!suspended sediment contribution [m3/m2]
        double precision bh(nkn)		!bed elevation
        real gskm(nkn)			!average sediment grainsize in node k

        double precision percbd(nlbdim,nscls,nkn)        !fraction of sediment [0,1]
        double precision bedn(nlbdim,3,nkn)  !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)

	double precision bdh(nkn)		!total elevation change [>0depo,<0ero]
        double precision dzco			!elevation variation due to compaction [m]
        double precision flx 			!bedload + suspended flux [m]
        double precision timecomp		!time step for compaction
        double precision gsa,gss
        double precision gsmax
        double precision gsmin
        integer k,is

!       -------------------------------------------------------------------
!       Initialize valiables
!       -------------------------------------------------------------------

        gsmax = 0.d0
        gsmin = 100.d0

        do is = 1,nscls
          gsmax = dmax1(gsmax,gs(is))*1.01d0
          gsmin = dmin1(gsmin,gs(is))*0.99d0
        end do

!       -------------------------------------------------------------------
!       Loop over nodes --> only bottom layer
!       -------------------------------------------------------------------

        do k = 1,nkn
         dzco = 0.d0
         bdh(k) = 0.d0

!        -------------------------------------------------------------------
!        Loop over grainsize
!        -------------------------------------------------------------------

         do is = 1,nscls

          flx = bflux(is,k) + sflux(is,k)
	  if (is_nan(flx)) flx = 0.
	  flx = min(flx,0.01D+0)
	  flx = max(flx,-0.01D+0)

!         -------------------------------------------------------------------
!         Multiply by morphological acceleration factor
!         -------------------------------------------------------------------

          flx = flx * morpho

!         -------------------------------------------------------------------
!         Rearrange bed characteristics
!         -------------------------------------------------------------------

	  call checkbed(k,is,nscls,gs(is),bedn(1,1,k),percbd(1,1,k),flx)

          bdh(k) = bdh(k) + flx

         end do

!        -------------------------------------------------------------------
!        Create new layer if first layer is to big
!        -------------------------------------------------------------------

         call newlayer(nscls,bedn(1,1,k),percbd(1,1,k))

!        -------------------------------------------------------------------
!        Merge layer 2 and 3 if layer 2 is too thin
!        -------------------------------------------------------------------

         call unilayer(nscls,bedn(1,1,k),percbd(1,1,k))

!        -------------------------------------------------------------------
!        Calculate compaction
!        -------------------------------------------------------------------

         if(DOCOMPACT .ne. 0) then 
	    timecomp = timedr * morpho
	    call compact(bedn(1,1,k),nlbdim,timecomp,dzco)
            bdh(k) = bdh(k) - dzco
	 end if

!        -------------------------------------------------------------------
!        Compute total bottom thickness variation
!        -------------------------------------------------------------------

	 if (is_nan(bdh(k))) bdh(k) = 0.
	 bdh(k) = min(bdh(k),0.01D+0)
	 bdh(k) = max(bdh(k),-0.01D+0)
	 bdh(k) = bdh(k) * MORPHO
         bh(k) = bh(k) + real(bdh(k))

!        -------------------------------------------------------------------
!        Update average sediment grainsize (D50)
!        -------------------------------------------------------------------

         gsa = 0.d0
         do is = 1,nscls
          gss = gs(is)
	  if (gss .ge. 0.10d0 .and. is .gt. 1) gss = gs(is-1)
          gsa = gsa + gss*percbd(1,is,k)
         end do
         gskm(k) = real(gsa)

         if(gskm(k).gt.gsmax.or.gskm(k).lt.gsmin) go to 120

!       -------------------------------------------------------------------
!       End of node loop
!       -------------------------------------------------------------------

        end do

        return

 120    continue
        write(6,*) 'Error in computing the average grainsize'
        write(6,*) 'gsaver:',gskm(k),'node:',k
        do is = 1,nscls
          write(6,*)'classnumber:',is,'grainsize:',gs(is)
        enddo
        stop 'error stop : gsaver'

        end

! ********************************************************************
! SUBROUTINE BEDACT
! If upper layer is smaller than the active layer --> create a new layer
! merging layer 1 and layer 2

        subroutine bedact(nscls,bmix,BEDCHA,PERCS)

        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05

        implicit none

        double precision PERCS(nlbdim,nscls)	!fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)         !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nlbd				!number of bed layer in k
        integer nscls                           !number of grainsize class
        double precision bmix                   !thickness of active layer [m]
        double precision ptot
        integer ib,is
	double precision d0
	parameter (d0 = 0.d0)
        double precision lin,x1,x2,y1,y2,x	!statement function which interpolate 
						!linearly the value Y of (X,Y), located 
						!on the line between (X1,Y1) and (X2,Y2)

        lin(x1,x2,y1,y2,x) = y1 + (x-x1)/(x2-x1)*(y2-y1)

!       -------------------------------------------------------------------
!       Find the number of rows that are used. This is indicated by 
!       z=0 in the row after the last used row.
!       -------------------------------------------------------------------

        do 10, ib=2,nlbdim
          if (BEDCHA(ib,1).eq.0.0d0) goto 20
10      continue
20      nlbd=ib-1

!       -------------------------------------------------------------------
!       CASE 1: BMIX < BEDCHA(2,1)
!       -------------------------------------------------------------------

	if( BEDCHA(2,1) .gt. bmix) then

!         -------------------------------------------------------------------
!         Update layer 1 (level 2) by linear interpolation
!         -------------------------------------------------------------------

          do is = 1,nscls
            PERCS(2,is) = (PERCS(1,is)*(BEDCHA(2,1)-bmix) + 
     $PERCS(2,is)*(BEDCHA(3,1)-BEDCHA(2,1)))/(BEDCHA(3,1)-bmix)
          end do

          BEDCHA(2,3) = lin(d0,BEDCHA(2,1),BEDCHA(1,3),
     $BEDCHA(2,3),bmix)
          BEDCHA(2,2) = lin(d0,BEDCHA(2,1),BEDCHA(1,2),
     $BEDCHA(2,2),bmix)
          BEDCHA(2,1) = bmix

	else

!       -------------------------------------------------------------------
!       CASE 2: BMIX > BEDCHA(2,1)
!       -------------------------------------------------------------------

!       -------------------------------------------------------------------
!       If level 3 is lower than bmix merge layer 1 and 2
!       -------------------------------------------------------------------

        if (BEDCHA(3,1).le.bmix) then
125       nlbd = nlbd - 1

          do is = 1,nscls			!update PERCS in layer 1
            PERCS(1,is) = (PERCS(2,is)*(BEDCHA(3,1)-BEDCHA(2,1)) +
     $PERCS(1,is)*BEDCHA(2,1))/BEDCHA(3,1)
          end do

          do ib = 2,nlbd			!shift one layer up
           do is = 1,nscls
            PERCS(ib,is) = PERCS(ib+1,is)
           end do
           BEDCHA(ib,1) = BEDCHA(ib+1,1)
           BEDCHA(ib,2) = BEDCHA(ib+1,2)
           BEDCHA(ib,3) = BEDCHA(ib+1,3)
          end do

          if (BEDCHA(3,1).le.bmix) goto 125

        end if

!       -------------------------------------------------------------------
!       Update layer 1 (level 2) by linear interpolation
!       -------------------------------------------------------------------

        do is = 1,nscls
          PERCS(1,is) = (PERCS(1,is)*BEDCHA(2,1) + PERCS(2,is)*
     $(bmix-BEDCHA(2,1)))/bmix
        end do

        BEDCHA(2,3) = lin(BEDCHA(2,1),BEDCHA(3,1),BEDCHA(2,3),
     $BEDCHA(3,3),bmix)
        BEDCHA(2,2) = lin(BEDCHA(2,1),BEDCHA(3,1),BEDCHA(2,2),
     $BEDCHA(3,2),bmix)
        BEDCHA(2,1) = bmix

	end if

!       -------------------------------------------------------------------
!       Adjusts PERCS and percentage check ---> ptot must be = 1 +-0.1
!       -------------------------------------------------------------------

        ptot = 0.d0
        do is = 1,nscls
          if(PERCS(1,is).lt.1D-5) PERCS(1,is) = 0.d0
	  PERCS(1,is) = dmin1(PERCS(1,is),1.d0)
         ptot = ptot + PERCS(1,is)
        end do
        if(ptot.lt.0.99d0.or.ptot.gt.1.01d0) go to 130

        do is = 1,nscls
         PERCS(1,is) = PERCS(1,is)/ptot
        end do

!       -------------------------------------------------------------------
!       Create last layer to the minimum layer number
!       -------------------------------------------------------------------

        if ( nlbd .le. 5 ) then
          do is = 1,nscls
            PERCS(nlbd+1,is) = PERCS(nlbd,is)
          end do
          BEDCHA(nlbd+1,1) = BEDCHA(nlbd,1) + (BEDCHA(nlbd,1) - 
     $BEDCHA(nlbd-1,1))
          BEDCHA(nlbd+1,2) = BEDCHA(nlbd,2)
          BEDCHA(nlbd+1,3) = BEDCHA(nlbd,3)
          nlbd = nlbd + 1
        end if

!       -------------------------------------------------------------------
!       Reset valiables of not used rows
!       -------------------------------------------------------------------

        BEDCHA(1,1) = 0.d0

        do ib = nlbd+1,nlbdim
          BEDCHA(ib,1) = 0.d0
          BEDCHA(ib,2) = 0.d0
          BEDCHA(ib,3) = 0.d0
        end do

	if(.not.checkbedcha(bedcha,nlbdim)) go to 133

        return

 130    continue
        write(6,*) 'Error in computing the sediment fraction: 1'
        write(6,*) 'total percentage:',ptot
        do is = 1,nscls
          write(6,*) 'classnumber:',is,'percentage:',PERCS(1,is)
        enddo
        stop 'error stop : PERCS'

 133    continue
        write(6,*) 'Error in BEDCHA: subroutine bedact'
        write(6,*) 'nlevel:',nlbd
        do ib = 1,nlbd
          write(6,*) 'number:',ib,'depth:',bedcha(ib,1)
        enddo
        stop 'error stop : bedcha'

        end

! ********************************************************************
! SUBROUTINE UNILAYER
! If layer 2 is smaller than 0.1 m merge layer 2 and layer 3

        subroutine unilayer(nscls,BEDCHA,PERCS)

        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05

        implicit none

        double precision PERCS(nlbdim,nscls)	!fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)         !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nlbd				!number of bed layer in k
        integer nscls                           !number of grainsize class
        double precision ptot
        integer ib,is
	double precision layer2			!layer 2 thickness (m)
	double precision minthick		!minimum layer thickness (m)

!       -------------------------------------------------------------------
!       Find the number of rows that are used. This is indicated by 
!       z=0 in the row after the last used row.
!       -------------------------------------------------------------------

        do 10, ib=2,nlbdim
          if (BEDCHA(ib,1).eq.0.0d0) goto 20
10      continue
20      nlbd=ib-1

	minthick = 0.001d0
	layer2 = BEDCHA(3,1) - BEDCHA(2,1)

!       -------------------------------------------------------------------
!       If layer 2 is thiner than the minimum value merge layer 2 and 3
!       -------------------------------------------------------------------

        if (layer2 .lt. minthick) then
          nlbd = nlbd - 1

          do is = 1,nscls			!update PERCS in layer 1
            PERCS(2,is) = (PERCS(2,is)*(BEDCHA(3,1)-BEDCHA(2,1)) +
     $PERCS(3,is)*(BEDCHA(4,1)-BEDCHA(3,1)))/(BEDCHA(4,1)-BEDCHA(2,1))
          end do

          do ib = 3,nlbd			!shift one layer up
           do is = 1,nscls
            PERCS(ib,is) = PERCS(ib+1,is)
           end do
           BEDCHA(ib,1) = BEDCHA(ib+1,1)
           BEDCHA(ib,2) = BEDCHA(ib+1,2)
           BEDCHA(ib,3) = BEDCHA(ib+1,3)
          end do

        end if

!       -------------------------------------------------------------------
!       Adjusts PERCS and percentage check of layer 2 ---> ptot must be = 1 +-0.1
!       -------------------------------------------------------------------

        ptot = 0.d0
        do is = 1,nscls
          if(PERCS(1,is).lt.1D-5) PERCS(1,is) = 0.d0
	  PERCS(1,is) = dmin1(PERCS(1,is),1.d0)
         ptot = ptot + PERCS(2,is)
        end do
        if(ptot.lt.0.99d0.or.ptot.gt.1.01d0) go to 130

        do is = 1,nscls
         PERCS(2,is) = PERCS(2,is)/ptot
        end do

!       -------------------------------------------------------------------
!       Create last layer to the minimum layer number
!       -------------------------------------------------------------------

        if ( nlbd .le. 5 ) then
          do is = 1,nscls
            PERCS(nlbd+1,is) = PERCS(nlbd,is)
          end do
          BEDCHA(nlbd+1,1) = BEDCHA(nlbd,1) + (BEDCHA(nlbd,1) - 
     $BEDCHA(nlbd-1,1))
          BEDCHA(nlbd+1,2) = BEDCHA(nlbd,2)
          BEDCHA(nlbd+1,3) = BEDCHA(nlbd,3)
          nlbd = nlbd + 1
        end if

!       -------------------------------------------------------------------
!       Reset valiables of not used rows
!       -------------------------------------------------------------------

        BEDCHA(1,1) = 0.d0

        do ib = nlbd+1,nlbdim
          BEDCHA(ib,1) = 0.d0
          BEDCHA(ib,2) = 0.d0
          BEDCHA(ib,3) = 0.d0
        end do

	if(.not.checkbedcha(bedcha,nlbdim)) go to 133

        return

 130    continue
        write(6,*) 'Error in computing the sediment fraction in layer 2'
        write(6,*) 'total percentage:',ptot
        do is = 1,nscls
          write(6,*) 'classnumber:',is,'percentage:',PERCS(2,is)
        enddo
        stop 'error stop : PERCS'

 133    continue
        write(6,*) 'Error in BEDCHA: subroutine minthick'
        write(6,*) 'nlevel:',nlbd
        do ib = 1,nlbd
          write(6,*) 'number:',ib,'depth:',bedcha(ib,1)
        enddo
        stop 'error stop : bedcha'

        end

! ********************************************************************
! SUBROUTINE CHECKBED
! The bed is divided in several layers, each with its own sediment 
! composition. At the beginning all the layers has the same volume.
! Layer = 1 is the top layer (close to water surface). 

        subroutine checkbed(k,iss,nscls,gs,BEDCHA,PERCS,flux)

        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05
        
        implicit none

        double precision PERCS(nlbdim,nscls)	!fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)  	!bed characteristics in 3 column table
                                               	! (1) depth below sediment surface (m)
                                               	! (2) critical erosion stress (Pa)
                                               	! (3) dry bulk density (kg/m**3)
        integer nlbd				!number of bed layer in k
        integer nscls				!number of grainsize class
        double precision gs			!grainsize of first class
        double precision flux			!flux of sediment [m] (<0 ero, >0 depo)
        double precision bnet			!new net sediment in the layer [m]
        double precision btot			!total sediment [m]
        double precision bpre                   !sediment is present in the layer [m]
        double precision ptot			!total percentage
        integer k,ib,is,iss

!       -------------------------------------------------------------------
!       Find the number of rows that are used (NBED). This is indicated by
!       z=0 in the row after the last used row.
!       -------------------------------------------------------------------

        do 10, ib=2,nlbdim
          if (BEDCHA(ib,1).eq.0.0d0) goto 20
10      continue
20      nlbd=ib-1

!       -------------------------------------------------------------------
!       Initialize valiables
!       -------------------------------------------------------------------

        btot = BEDCHA(2,1) + flux
        bpre = BEDCHA(2,1)*PERCS(1,iss)
        bnet = bpre  + flux

!       -------------------------------------------------------------------
!       Compute erosion or deposition
!       -------------------------------------------------------------------

        if (flux.gt.0.d0) then

          call depbed(iss,nlbd,nscls,BEDCHA,PERCS,flux,bnet,btot,gs)

        elseif(flux.lt.0.d0) then

          call erobed(iss,nlbd,nscls,BEDCHA,PERCS,flux,bnet,btot,bpre)

        end if
  
!       -------------------------------------------------------------------
!       Adjust PERCS and percentage check ---> ptot must be = 1 +-0.1
!       -------------------------------------------------------------------

        if (nscls.eq.1) PERCS(1,nscls) = 1.d0

        ptot = 0.d0
        do is = 1,nscls
          if(PERCS(1,is).lt.1D-5) PERCS(1,is) = 0.d0
	  PERCS(1,is) = dmin1(PERCS(1,is),1.d0)
          ptot = ptot + PERCS(1,is)
        end do
        if(ptot.lt.0.99d0.or.ptot.gt.1.01d0) go to 130

        do is = 1,nscls
         PERCS(1,is) = PERCS(1,is)/ptot
        end do

!       -------------------------------------------------------------------
!       Reset valiables of not used rows
!       -------------------------------------------------------------------

        BEDCHA(1,1) = 0.d0

        do ib = nlbd+1,nlbdim
         BEDCHA(ib,1) = 0.d0
         BEDCHA(ib,2) = 0.d0
         BEDCHA(ib,3) = 0.d0
        end do

	if(.not.checkbedcha(bedcha,nlbdim)) goto 133

        return

 130    continue
        write(6,*) 'Error in computing the sediment fraction: 2'
        write(6,*) 'total percentage:',ptot,'node:',k,'flux:',flux
        do is = 1,nscls
          write(6,*) 'classnumber:',is,'percentage:',PERCS(1,is)
        enddo
        stop 'error stop : PERCS'

 133    continue
        write(6,*) 'Error in BEDCHA: subroutine checkbed'
        write(6,*) 'node:',k, 'level:',nlbd
        write(6,*) 'layer  depth  tauce  rhos'
        do ib = 1,nlbd
          write(6,*) ib,bedcha(ib,1),bedcha(ib,2),bedcha(ib,3)
        enddo
        stop 'error stop : bedcha'

        end

! ********************************************************************

        subroutine depbed(iss,nlbd,nscls,BEDCHA,PERCS,flux,bnet,btot,
     $			  gss)

        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05

        implicit none

        double precision PERCS(nlbdim,nscls)   !fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)       !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nlbd                            !number of bed layer in k
        integer nscls                           !number of grainsize class
        double precision flux                   !flux of sediment [m] (<0 ero, >0 depo)
        double precision bnet                   !new net sediment in the layer [m]
        double precision btot                   !total sediment [m]
        double precision gss                    !grainsize of first class
        double precision rhosa                  !non consolidated density [kg/m3]
	double precision taunew			!new critical shear stress
	double precision tausurf		!average surface critical shear stress
	double precision rhosurf		!average surface density
        integer ib,is,iss

!       -------------------------------------------------------------------
!       Initialize valiables
!       -------------------------------------------------------------------

        rhosa = RHOSED * (1.d0 - SURFPOR)
        if (gss.lt.limcoh) rhosa = RHOMUD
        taunew = teroa*(rhosa**terob)
        !tausurf = (BEDCHA(1,2) + BEDCHA(2,2))/2.
        tausurf = BEDCHA(1,2)
        !rhosurf = (BEDCHA(1,3) + BEDCHA(2,3))/2.
        rhosurf = BEDCHA(1,3)

!       -------------------------------------------------------------------
!       Normal deposition
!       -------------------------------------------------------------------

        do is = 1,nscls
          PERCS(1,is) = PERCS(1,is)*BEDCHA(2,1)/btot
        end do

        PERCS(1,iss) = bnet / btot

        taunew = (taunew*flux + tausurf*BEDCHA(2,1)) / btot
        BEDCHA(1,2) = dmin1(BEDCHA(1,2),taunew)
        BEDCHA(1,2) = dmin1(BEDCHA(1,2),BEDCHA(2,2))

        BEDCHA(1,3) = (rhosa*flux + rhosurf*BEDCHA(2,1))/ btot
        BEDCHA(1,3) = dmin1(BEDCHA(1,3),BEDCHA(2,3))

        do ib = 2,nlbd
          BEDCHA(ib,1) = BEDCHA(ib,1) + flux
        end do

	end 

! ********************************************************************
! Excess deposition, create new layer and upgrade level 2

        subroutine newlayer(nscls,BEDCHA,PERCS)

        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05

        implicit none

        double precision PERCS(nlbdim,nscls)   !fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)         !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nscls                           !number of grainsize class
        double precision bres                   !sediment that is left [m]
        double precision blimit                 !maximum thikness of layer [m]
        integer nlbd                            !number of bed layer in k
        integer ib,is
	double precision d0
	parameter (d0 = 0.d0)
        double precision lin,x1,x2,y1,y2,x      !statement function which interpolate
                                                !linearly the value Y of (X,Y), located
                                                !on the line between (X1,Y1) and (X2,Y2)

        lin(x1,x2,y1,y2,x) = y1 + (x-x1)/(x2-x1)*(y2-y1)

!       -------------------------------------------------------------------
!       Find the number of rows that are used. This is indicated by 
!       z=0 in the row after the last used row.
!       -------------------------------------------------------------------

        do 10, ib=2,nlbdim
          if (BEDCHA(ib,1).eq.0.0d0) goto 20
10      continue
20      nlbd=ib-1

!       -------------------------------------------------------------------
!       Initialize valiables
!       -------------------------------------------------------------------

        blimit = (BEDCHA(3,1) - BEDCHA(2,1))/2.
	bres = bedcha(2,1) - blimit

	if ( bres .gt. (blimit/2.) ) then

!         -------------------------------------------------------------------
!         Create a new to layer and shift the others down
!         -------------------------------------------------------------------

          nlbd = nlbd + 1	
          do ib = nlbd,3,-1
           do is = 1,nscls
            PERCS(ib,is) = PERCS(ib-1,is)
           end do
           BEDCHA(ib,1) = BEDCHA(ib-1,1)
           BEDCHA(ib,2) = BEDCHA(ib-1,2)
           BEDCHA(ib,3) = BEDCHA(ib-1,3)
          end do

!         -------------------------------------------------------------------
!         Update layer 2 properties
!         -------------------------------------------------------------------

          do is = 1,nscls
            PERCS(2,is) = PERCS(1,is)
          end do
          BEDCHA(2,1) = bres
	  BEDCHA(3,1) = BEDCHA(4,1) - blimit
          BEDCHA(2,2) = lin(d0,BEDCHA(3,1),BEDCHA(1,2),BEDCHA(3,2),
     $			BEDCHA(2,1))
          BEDCHA(2,3) = lin(d0,BEDCHA(3,1),BEDCHA(1,3),BEDCHA(3,3),
     $			BEDCHA(2,1))

          nlbd = min(nlbdim-2,nlbd)

	end if

!       -------------------------------------------------------------------
!       Reset valiables of not used rows
!       -------------------------------------------------------------------

        BEDCHA(1,1) = 0.d0

        do ib = nlbd+1,nlbdim
         BEDCHA(ib,1) = 0.d0
         BEDCHA(ib,2) = 0.d0
         BEDCHA(ib,3) = 0.d0
        end do

	if(.not.checkbedcha(bedcha,nlbdim)) goto 133

	return

 133    continue
        write(6,*) 'Error in BEDCHA: subroutine newlayer'
        write(6,*) 'nlevel:',nlbd
        do ib = 1,nlbd
          write(6,*) 'number:',ib,'depth:',bedcha(ib,1)
        enddo
        stop 'error stop : bedcha'

        end

! ********************************************************************
       
        subroutine erobed(iss,nlbd,nscls,BEDCHA,PERCS,flux,bnet,btot,
     $bpre)

        use mod_sediment
        use mod_sediment_para

        implicit none

	double precision d0
	parameter (d0 = 0.d0)

        double precision PERCS(nlbdim,nscls)   !fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)       !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nlbd                            !number of bed layer in k
        integer nscls                           !number of grainsize class
        double precision flux                   !flux of sediment [m] (<0 ero, >0 depo)
        double precision bnet                   !new net sediment in the layer [m]
        double precision btot                   !total sediment [m]
        double precision bpre                   !sediment iss present in the layer [m]
        integer ib,is,iss

        double precision lin,x1,x2,y1,y2,x      !statement function which interpolate
                                                !linearly the value Y of (X,Y), located
                                                !on the line between (X1,Y1) and (X2,Y2)

        lin(x1,x2,y1,y2,x) = y1 + (x-x1)/(x2-x1)*(y2-y1)

!       -------------------------------------------------------------------
!       No iss sediment present --> reset flux and go back
!       -------------------------------------------------------------------

        if(bpre .eq. d0) then
	  flux = d0
	  return
        end if

!       -------------------------------------------------------------------
!       Excess erosion, limit erosion
!       -------------------------------------------------------------------

        if (bnet.le.d0) then
          flux = -bpre
	  bnet = d0
          btot = BEDCHA(2,1) + flux
	  if (btot. eq. d0 .or. PERCS(1,iss) .gt. 0.999d0) then
	    call dellayer(nlbd,nscls,BEDCHA,PERCS)
	    return
	  endif
        endif

!       -------------------------------------------------------------------
!       Normal erosion
!       -------------------------------------------------------------------

        do is = 1,nscls
          PERCS(1,is) = PERCS(1,is)*BEDCHA(2,1) / btot
        end do

        PERCS(1,iss) = bnet / btot

        BEDCHA(1,2) = lin(d0,BEDCHA(2,1),BEDCHA(1,2),BEDCHA(2,2),-flux)
        BEDCHA(1,3) = lin(d0,BEDCHA(2,1),BEDCHA(1,3),BEDCHA(2,3),-flux)

        do ib = 2,nlbd
          BEDCHA(ib,1) = BEDCHA(ib,1) + flux
        end do

        end

! ********************************************************************
! SUBROUTINE DELLAYER
! Delete layers and update bed characteristics

        subroutine dellayer(nlbd,nscls,BEDCHA,PERCS)

        use mod_sediment
        use mod_sediment_para

        implicit none

        double precision PERCS(nlbdim,nscls)   !fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)       !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nlbd                            !number of bed layer in k
        integer nscls                           !number of grainsize class
        double precision thick
        integer ib,is

        thick = BEDCHA(2,1)

        nlbd = nlbd - 1

        if ( nlbd .le. 5 ) then
          nlbd = nlbd + 1
          do is = 1,nscls
            PERCS(nlbd+1,is) = PERCS(nlbd,is)
          end do
          BEDCHA(nlbd+1,1) = BEDCHA(nlbd,1) + (BEDCHA(nlbd,1) -
     $BEDCHA(nlbd-1,1))
          BEDCHA(nlbd+1,2) = BEDCHA(nlbd,2)
          BEDCHA(nlbd+1,3) = BEDCHA(nlbd,3)
        end if

        do ib = 1,nlbd                       !shift one layer up
         do is = 1,nscls
          PERCS(ib,is) = PERCS(ib+1,is)
         end do
         BEDCHA(ib,1) = BEDCHA(ib+1,1) - thick
         BEDCHA(ib,2) = BEDCHA(ib+1,2)
         BEDCHA(ib,3) = BEDCHA(ib+1,3)
        end do

!       -------------------------------------------------------------------
!       Reset valiables of not used rows
!       -------------------------------------------------------------------

        BEDCHA(1,1) = 0.d0

        do ib = nlbd+1,nlbdim
          BEDCHA(ib,1) = 0.d0
          BEDCHA(ib,2) = 0.d0
          BEDCHA(ib,3) = 0.d0
        end do

        end

! ********************************************************************
! SUBROUTINE TOTCON
! Compute the total node concentration

        subroutine totcon(nscls,scn)

	use mod_ts
	use levels
	use basin, only : nkn
        use mod_sediment
        use mod_sediment_para
	use mod_sedtrans05

        implicit none

        integer nscls				!number of grainsize class
        double precision scn(nlvdi,nkn,nscls)  	!suspended concentration

	real rhos,conc
        integer k,l,is,lmax

        do k = 1,nkn
          lmax = ilhkv(k)
          do l = 1,lmax
             tcn(l,k) = 0.
             do is = 1,nscls
              scn(l,k,is) = max(scn(l,k,is),0.d0)
              conc = scn(l,k,is)
	      if (conc < 0.d0 ) then
                write(6,*) 'l,k,is,scn: ',l,k,is,scn(l,k,is)
                stop 'error stop totcon: negative concentration'
              end if
  
  	      if (is .le. nbcc) then
  	        rhos = rhoclay
  	      else
	        rhos = rhosed
	      end if

!             -------------------------------------------------------------------
!             Compute total concentration
!             -------------------------------------------------------------------
              tcn(l,k) = tcn(l,k) + conc

!             -------------------------------------------------------------------
!             Update fluid density
!             -------------------------------------------------------------------
              rhov(l,k) = rhov(l,k) + (conc/rhos)*(rhos-rhov(l,k))

            end do
          end do
        end do

        call nantest(nkn*nlvdi,rhov,'rhov')
            
        end

! ********************************************************************
! SUBROUTINE UPEDEPTH
! Update the element depth in function of erosion/deposition

        subroutine upedepth(bdh)

	use mod_geom_dynamic
	use mod_depth
	use mod_layer_thickness
	use mod_area
	use mod_hydro
        use basin, only : nkn,nel,hm3v,nen3v
        use levels


        implicit none

        double precision bdh(nkn)    !total elevation change [>0:depo,<0:eros]

        !real hlhv(nel)	!GGU

        real dh
        real evdep				!element depth variation
        integer ie,ii,k,iw

	real bdhrst(nkn)
	integer nkk,nv,iv
        character*80 name
	integer icall
        save icall
        data icall /1/

!       -------------------------------------------------------------------
!       Read restart bed elevation change
!       -------------------------------------------------------------------

        if( icall .eq. 0 ) then
	  name='bed_change_time0.dat'
          open(999,file=name,form='unformatted')
          read(999) nkk,nv,iv
          read(999) (bdhrst(k),k=1,nkk)
 	  write(*,*)'Initialize bed elevation change from file:',name
	  do k = 1,nkn
	    bdh(k) = bdh(k) + bdhrst(k)
	  end do
          icall = 1
        end if

!       -------------------------------------------------------------------
!       Update element's depth
!       -------------------------------------------------------------------

        do ie = 1,nel
          evdep = 0.
          do ii = 1,3
            k = nen3v(ii,ie)
	    if( iwegv(ie) .gt. 0 ) then
	      dh = 0.
	    else
	      dh = real(bdh(k))
	    end if
            evdep = evdep + dh
            !hm3v(ii,ie) = hm3v(ii,ie) - dh	????
          end do
          evdep = evdep / 3.
          !if ((hlhv(ie) - evdep) .lt. 0.2 ) evdep = 0.

          !hlhv(ie) = hlhv(ie) - evdep
	  do ii = 1,3
            hm3v(ii,ie) = hm3v(ii,ie) - evdep
          end do
	  hev(ie) = hev(ie) - evdep

        end do

!       ------------------------------------------------------------------
!       Set up depth vectors
!       ------------------------------------------------------------------

        call makehkv(hkv)         !computes hkv as average

	call setweg(3,iw)
        call set_area
        call setdepth(nlvdi,hdknv,hdenv,zenv,areakv)

!       ------------------------------------------------------------------
!       Set up velocity
!       ------------------------------------------------------------------

        !call ttov
	!call setuvd

        end

! ********************************************************************
! Smooths variable averaging over the neibor nodes
! weight is the node`s area. Smooths only if bed elevation changes occour
! or if bed slope is greater than angle of repose

        subroutine smooth_node(kvalue,kdiff,smooth,gdx,gdy,angrep)

	use mod_geom
	use levels
	use basin, only : nkn

        implicit none

        double precision kvalue(nkn)    !variable
        double precision kdiff(nkn)	!istantaneous variable difference
        double precision smooth		!smoothing factor for morphodynamic [0-1]
        real gdx(nkn),gdy(nkn)		!slope gradients
        double precision angrep		!angle of repose [rad]
        real ksl                        !node slope angle [radian]
        real areanode,area
        real baux(nkn)
        real saux,aaux,smm

        integer k,kn,n,ibase,i,l
	integer nodes(maxlnk)

        if ( smooth .eq. 1.d0 ) return

	baux = 0.

        do k = 1,nkn
          saux = 0.
          aaux = 0.
          smm = real(smooth)

          ksl = atan(sqrt(gdx(k)**2 + gdy(k)**2))

          if (dabs(kdiff(k)).gt.1d-5) then
            if(ksl.ge.real(angrep)) smm = 0.5
	    call get_nodes_around(k,maxlnk,n,nodes)
            do i = 1,n
              kn = nodes(i)           !kn is number of neibor node
              l = ilhkv(kn)
              area = areanode(l,kn)
              saux = saux + kvalue(kn)*area
              aaux = aaux + area
            end do
            baux(k) = kvalue(k)*smm + (saux / aaux)*(1. - smm)
  
          else
            baux(k) = kvalue(k)
          end if
        end do

        do k = 1,nkn
          kdiff(k) = kdiff(k) + (baux(k)-kvalue(k))
          kvalue(k) = baux(k)
        end do

        end

! ********************************************************************
! Reset bottom height variation and ssc after initialization time adjtime,
!  while all other variables (sediment density, bathymetry, tauce) are not.

        subroutine resetsedi(dtime,adjtime,nscls,bh,scn)

	use basin, only : nkn
        use levels, only : nlvdi

        implicit none

	double precision, intent(in)	:: dtime	!time in seconds
	double precision, intent(in)	:: adjtime	!time for initialization [s]
	integer, intent(in)		:: nscls	!number of grainsize classes
	double precision, intent(inout) :: bh(nkn)      !bottom height variation [m]
        double precision, intent(inout) :: scn(nlvdi,nkn,nscls) !suspended sediment conc (kg/m3)

	if (dtime == adjtime) then 
	  write(6,*)'Morphological initialization time:'
	  write(6,*)'reset bathymetric change but keep modified'
	  write(6,*)'sediment characteristics and bathymetry'
	  bh  = 0.d0
  	  scn = 0.d0
	end if

	end subroutine resetsedi 

!******************************************************************
! Compute the sediment bubdet in the basin considering suspended and
! bed loads.

        subroutine sedbudget(ius,dtime,timedr,nscls,sload,bload,bedk)

        use levels
        use basin, only : nkn
        use mod_sediment

        implicit none

        integer, intent(in)          :: ius          !unit for budget file
        double precision, intent(in) :: dtime	     !time in seconds
        double precision, intent(in) :: timedr	     !time step [s]
	integer, intent(in) 	     :: nscls	     !number of grainsize classes
        double precision, intent(in) :: sload(nkn,nscls) !suspended sediment load [kg/s]
        double precision, intent(in) :: bload(nkn)   !flux of bedload sediment [kg]

        double precision, intent(inout) :: bedk(nkn) !sediment budget in kg [>0depo,<0ero]

	integer		 :: lmax
	double precision :: load
        real 		 :: sus,bed
        real 		 :: volnode,vol
        integer		 :: k,l
        character*20 aline
        external get_timeline

        call get_timeline(dtime,aline)

        sus = 0.
        bed = 0.
        do k = 1,nkn
          load = - sum(sload(k,:))*timedr + bload(k)
          bedk(k) = bedk(k) + load
          bed = bed + bedk(k)

          lmax = ilhkv(k)
	  do l = 1,lmax
            vol = volnode(l,k,1)
            sus = sus + tcn(l,k)*vol
	  end do
        end do

        write(ius,1000) aline,sus,bed

        return
 1000   format(a,2f18.2)

        end subroutine sedbudget

!******************************************************************

	subroutine uvbott(ddl,ukbot,vkbot)
!
! compute bottom velocity on nodes
!
	use mod_geom_dynamic
	use mod_layer_thickness
	use mod_hydro_vel
	use evgeom
	use levels
        use basin, only : nkn,nel,nen3v

	implicit none

! arguments
	real vv(nkn)
	real vv1(nkn)
! local
	real ukbot(nkn),vkbot(nkn)
	real u,v
	real ddl(nkn)
	integer ie,lmax,k,ii
	real aj,vol,depth
	integer nsigma
        real hsigma
	logical bsigma

	call get_sigma_info(nlv,nsigma,hsigma)
	bsigma = nsigma .gt. 0

	do k = 1,nkn
	  vv(k)    = 0.
	  vv1(k)   = 0.
	  ddl(k)   = 0.
	  ukbot(k) = 0.
	  vkbot(k) = 0.
	end do

	if ( bsigma ) then
	  do k = 1,nkn
	    lmax = ilhkv(k)
            ddl(k)   = hdknv(lmax,k)
            call getuv(lmax,k,u,v)
	    ukbot(k) = u
	    vkbot(k) = v
	  end do
	  return
	end if

	do ie = 1,nel
	   lmax = ilhv(ie)
	   aj = real(ev(10,ie))
	   depth = hdenv(lmax,ie)
	   vol = aj*depth
	    do ii = 1,3
	      k = nen3v(ii,ie)
	      vv1(k) = vv1(k) + aj
	      ddl(k) = ddl(k) + vol
	    end do
           !if( iwegv(ie) .eq. 0 ) then
	   if( iwetv(ie) .gt. 5 .and. vol.gt.0. ) then
	    do ii = 1,3
	      k = nen3v(ii,ie)
	      vv(k)  = vv(k) + vol
	      ukbot(k) = ukbot(k) + aj*ulnv(lmax,ie)
	      vkbot(k) = vkbot(k) + aj*vlnv(lmax,ie)
	    end do
	   end if
	end do

	do k=1,nkn
	  if(vv(k).gt.0.) then
	    ukbot(k)=ukbot(k)/vv1(k)
	    vkbot(k)=vkbot(k)/vv1(k)
	  end if
	  if(vv1(k).gt.0.) then
	    ddl(k)=ddl(k)/vv1(k)
	  end if
	end do

	return
	end

!******************************************************************
! SUBROUTINE SET_SETTLING
! Set settling velocity array

        subroutine set_settling(nscls,ws,scn,wsink)

        use levels
        use basin, only : nkn
        use mod_sediment_para
	use mod_sedtrans05

        implicit none

        integer nscls				!number of grainsize class
        double precision ws(nscls)      	!SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        double precision scn(nlvdi,nkn,nscls)  	!suspended concentration
        double precision wsink(0:nlvdi,nkn,nscls) !settling velocity for suspended sediment

	double precision conc(nscls)
        real salt,temp			!salinity [psu] and temperature [C]
        DOUBLE PRECISION TEM		!in degrees Celsius
        DOUBLE PRECISION SAL		!practical salinity scale (seawater=35)
        DOUBLE PRECISION VISC   	!dynamic viscosity of the (sea)water (KG/M*SEC (OR N.S/M**2))
        DOUBLE PRECISION VISK     	!KINEMAMIC VISCOSITY OF THE FLUID (KG/M*SEC (OR N.S/M**2))
        DOUBLE PRECISION RHOW     	!DENSITY OF FLUID (WATER)  (KG/M**3)
        DOUBLE PRECISION, allocatable :: WWS(:)  !settling velocity of the considered floc class (m/s)

        integer k,l,is,i,lmax

	allocate(WWS(nbcc))

	do k = 1,nkn
          lmax = ilhkv(k)
	
 	  if (nbcc .gt. 0.) then
	   do l = 1,lmax
	     do i = 1,nbcc
               CONC(i) = scn(l,k,i)
	       WSI(i) = ws(i) 
             end do
             call getts(l,k,temp,salt)               !gets temp and salt
             tem = temp
             sal = salt
             call densvisc(tem,sal,rhow,visc)        !get density and viscosity
             VISK = visc/rhow                        !get viscosity
	     call WSINKFLOC(NBCC,CONC,RHOW,VISK,WWS)
	     do i = 1,nbcc
               wsink(l,k,i) = WWS(i)
             end do
	   end do
	  end if

	  do l = 1,lmax
	    do is = nbcc+1,nscls
              wsink(l,k,is) = ws(is)
	    end do
	  end do

	  do is = 1,nscls
	    wsink(0,k,is) = wsink(1,k,is)
	  end do

	end do

	end


!******************************************************************
! Set settling velocity for fine sediments in function of water temp.
! This is specific for the Curonian Lagoon

        subroutine wsinkbio(tem,wwsb)

	implicit none

	double precision, intent(in)	:: tem		!water temp
	double precision, intent(out)	:: wwsb		!settling vel

	if (tem > 5.) then
 	  wwsb = -0.03443*tem + 1.251117		!m/day
	  wwsb = wwsb / 86400.			! to m/s
	else
	  wwsb = 1.
	end if
	wwsb = max(wwsb,1d-06)

	end subroutine wsinkbio

!******************************************************************

        subroutine scn_compute(idsedi,nscls,sload,wsink,eps,scn)

        use mod_diff_visc_fric
        use levels, only : nlvdi,nlv,ilhkv
        use basin, only : nkn
	use mod_sediment
        use mod_sediment_para

        implicit none

        include 'mkonst.h'

        integer, intent(in)		:: idsedi(*)
	integer, intent(in)		:: nscls
	double precision, intent(inout)	:: sload(nkn,nscls)    !suspended sediment load [kg/s]
        double precision, intent(in)	:: wsink(0:nlvdi,nkn,nscls) !settling velocity for suspended sediment
        double precision, intent(in)	:: eps(0:nlvdi,nkn,nscls) !vertical mixing coefficient

        double precision, intent(inout)	:: scn(nlvdi,nkn,nscls)   !suspended concentration

        real				:: fact = 1.
        real, allocatable 		:: scal(:,:)   !suspended sediment conc for adv/diff
        real, allocatable 		:: epss(:,:)   !vertical mixing coefficient adv/diff
        real, allocatable 		:: wsinks(:,:) !settling velocity for suspended sediment adv/diff
        real, allocatable 		:: lload(:)
        real, allocatable 		:: load(:,:)
        real				:: dt
        double precision 		:: dtime
        real            		:: sedkpar	!diffusion parameter
	integer				:: is

!$OMP PARALLEL
!$OMP SINGLE 
!$OMP TASKWAIT

!!!$OMP TASKGROUP

	call get_act_dtime(dtime)
	call get_timestep(dt)

        sedkpar = sedpa(3)		!diffusion parameter

        call bnds_read_new(what,idsedi,dtime)

        do is = 1,nscls

!!!$OMP TASK DEFAULT(FIRSTPRIVATE) SHARED(sload,scn,eps,wsink,idsedi)

! !$OMP TASK FIRSTPRIVATE(is,fact,sedkpar,difhv,difmol,what
! !$OMP& ,dt,nlvdi,nkn,ilhkv,scal,wsinks,epss,lload,load)
! !$OMP& SHARED(sload,scn,eps,wsink,idsedi) DEFAULT(NONE)

!   	  allocate(lload(nkn))
!   	  allocate(load(nlvdi,nkn))
! 	  allocate(scal(nlvdi,nkn))
! 	  allocate(epss(0:nlvdi,nkn))
! 	  allocate(wsinks(0:nlvdi,nkn))

! 	  lload(:) = sload(:,is)
!           call load3d(lload,nkn,nlvdi,ilhkv,load)
!           scal(:,:)   = scn(:,:,is)
!           epss(:,:)   = eps(:,:,is)
!           wsinks(:,:) = wsink(:,:,is)
! 	  where( scal < 0. ) scal = 0.			!GGUZ0
!          call scal_adv_fact(what,is,fact
!      +                      ,scal,idsedi
!      +                      ,sedkpar,fact,wsinks,fact,load
!      +                      ,difhv,epss,difmol)

!          call sload3d(load,nkn,nlvdi,ilhkv,lload)
!          scn(:,:,is) = scal(:,:)
!          sload(:,is) = lload(:)

!  	  deallocate(lload)
!  	  deallocate(load)
! 	  deallocate(scal)
! 	  deallocate(epss)
! 	  deallocate(wsinks)

! !$OMP END TASK

! !$OMP TASK FIRSTPRIVATE(what,is,sedkpar,difmol)
! !$OMP& SHARED(scn,wsink,sload,eps,idsedi) DEFAULT(NONE)

!$OMP TASK DEFAULT(FIRSTPRIVATE) SHARED(idsedi,scn,wsink,sload,eps)
          call scn_compute_intern(what
     +			,is,idsedi,sedkpar,difmol
     +			,scn(:,:,is),wsink(:,:,is)
     +			,sload(:,is),eps(:,:,is))
!$OMP END TASK

        end do

!!!$OMP END TASKGROUP     

!$OMP TASKWAIT
!$OMP END SINGLE 
!$OMP END PARALLEL  

! -------------------------------------------------------------
! end of routine
! -------------------------------------------------------------

        end subroutine scn_compute

!******************************************************************

        subroutine scn_compute_intern(what
     +			,is,idsedi,sedkpar,difmol
     +			,scn,wsink,sload,eps)

        use mod_diff_visc_fric
        use levels, only : nlvdi,nlv,ilhkv
        use basin, only : nkn

	implicit none

	character*(*) what
	integer is,idsedi(*)
	real sedkpar,difmol
        double precision, intent(inout)	:: scn(nlvdi,nkn)
        double precision, intent(in)	:: wsink(0:nlvdi,nkn)
	double precision, intent(inout)	:: sload(nkn)
        double precision, intent(in)	:: eps(0:nlvdi,nkn)

        real, allocatable 		:: scal(:,:)
        real, allocatable 		:: epss(:,:)
        real, allocatable 		:: wsinks(:,:)
        real, allocatable 		:: lload(:)
        real, allocatable 		:: load(:,:)

	real, parameter :: fact = 1.

  	allocate(lload(nkn))
  	allocate(load(nlvdi,nkn))
	allocate(scal(nlvdi,nkn))
	allocate(epss(0:nlvdi,nkn))
	allocate(wsinks(0:nlvdi,nkn))

	lload(:) = sload(:)
        call load3d(lload,nkn,nlvdi,ilhkv,load)
        scal(:,:)   = scn(:,:)
        epss(:,:)   = eps(:,:)
        wsinks(:,:) = wsink(:,:)
	where( scal < 0. ) scal = 0.			!GGUZ0

        call scal_adv_fact(what,is,fact
     +                      ,scal,idsedi
     +                      ,sedkpar,fact,wsinks,fact,load
     +                      ,difhv,epss,difmol)

        call sload3d(load,nkn,nlvdi,ilhkv,lload)
        scn(:,:) = scal(:,:)
        sload(:) = lload(:)

  	deallocate(lload)
  	deallocate(load)
	deallocate(scal)
	deallocate(epss)
	deallocate(wsinks)

	end subroutine scn_compute_intern

!******************************************************************

	subroutine load3d(sload,nkn,nlvdi,ilhkv,load)

	implicit none

	real, intent(in)		:: sload(nkn)	!suspended sediment load at bottom
	integer, intent(in)		:: nkn		!number of node
	integer, intent(in)		:: nlvdi	!number of vertical levels
        integer, intent(in)		:: ilhkv(nkn)	!number of element and node level

	real, intent(out)		:: load(nlvdi,nkn) !suspended sediment load (3D)

	integer		::  k,l,lmax

	load = 0. 
	do k = 1,nkn
	   lmax = ilhkv(k)
	   load(lmax,k) = sload(k)
	end do

	end subroutine load3d

!******************************************************************

	subroutine sload3d(load,nkn,nlvdi,ilhkv,sload)

	implicit none

	real, intent(in)		:: load(nlvdi,nkn) !suspended sediment load (3D)
	integer, intent(in)		:: nkn		!number of node
	integer, intent(in)		:: nlvdi	!number of vertical levels
        integer, intent(in)		:: ilhkv(nkn)	!number of element and node level

	real, intent(out)		:: sload(nkn)	!suspended sediment load at bottom

	integer		::  k,l,lmax

	do k = 1,nkn
	   lmax = ilhkv(k)
	   sload(k) = load(lmax,k)
	end do

	end subroutine sload3d

!******************************************************************
! Update sflx in case of deposition after advection/diffusion of 
! suspended sediment

        subroutine update_sflx(nscls,sload,timedr,sflx)

        use basin, only : nkn
        use levels
        use mod_sediment_para
        use mod_sedtrans05

        implicit none

        integer, intent(in)     :: nscls
        double precision, intent(in)    :: sload(nkn,nscls) !suspended sediment load [kg/s]
        double precision, intent(in)    :: timedr       !time step [s]
        double precision, intent(inout) :: sflx(nscls,nkn)  !flux of suspend sediment [m3/m2]

	integer		 :: k,is,lmax
        real 		 :: areanode,area
	double precision :: load,RHOA

        do k = 1,nkn
          lmax = ilhkv(k)                       !bottom layer
          area = areanode(lmax,k)               !area of node
	  do is = 1,nscls
             if ( sload(k,is) < 0.d0 ) then
               if ( is <= nbcc ) then
                 RHOA = RHOMUD
               else
                 RHOA = RHOSED * (1.d0 - SURFPOR)
               end if
               load = -sload(k,is) / area
               sflx(is,k) =  load * timedr / RHOA
             end if
           end do
	end do

	end subroutine update_sflx

!******************************************************************

	subroutine mud_sand_fract(nlbdim,nscls,nkn,nbcc,percbd,
     +		   percc,percs)

	implicit none

	integer, intent(in)	:: nlbdim
	integer, intent(in)	:: nscls
	integer, intent(in)	:: nkn
	integer, intent(in)	:: nbcc
        double precision, intent(in) :: percbd(nlbdim,nscls,nkn) !fraction of sediment [0,1]

	real, intent(out)	:: percc(nkn)	! fraction of mud
	real, intent(out)	:: percs(nkn)	! fraction od sand

	integer			:: k,is
	double precision	:: pcoes

	do k = 1,nkn
	  pcoes = 0.d0
	  do is = 1,nbcc
	    pcoes = pcoes + percbd(1,is,k)
	  end do
	  percc(k) = pcoes
	  percs(k) = 1. - pcoes
        enddo

	end subroutine mud_sand_fract

!******************************************************************

        subroutine init_sed_output

	use mod_sediment

	implicit none

	integer nvar,id
	logical, parameter :: b2d = .true.
	logical, parameter :: b3d = .false.

	logical has_output_d

        nvar = 1
        call init_output_d('itmsed','idtsed',da_ssc)
        if( has_output_d(da_ssc) ) then
          call shyfem_init_scalar_file('ssc',nvar,b3d,id)
          da_ssc(4) = id
        end if

        nvar = 5
        call init_output_d('itmsed','idtsed',da_sed)
        if( has_output_d(da_sed) ) then
          call shyfem_init_scalar_file('sed',nvar,b2d,id)
          da_sed(4) = id
        end if

	end

!******************************************************************

        subroutine wr_sed_output(dtime,bh,gskm,tao,percc,percs,totbed)

        use levels, only : nlvdi,nlv
        use basin, only : nkn
        use mod_depth
	use mod_sediment

        implicit none

        double precision, intent(in)	:: dtime	!simulation time
        double precision, intent(in)	:: bh(nkn)	!bottom height var
        real, intent(in)		:: gskm(nkn)	!average grainsize
        real, intent(in)		:: tao(nkn)	!wave-cur shear stress
	real, intent(in)		:: percc(nkn)	!fraction of mud
	real, intent(in)		:: percs(nkn)	!fraction od sand
	real, intent(in)		:: totbed(nkn)	!total bedload transport

        logical				:: next_output_d
        integer	 			:: id

        if( next_output_d(da_ssc) ) then
          id = nint(da_ssc(4))
          call shy_write_scalar_record(id,dtime,800,nlv,tcn)
	  call shy_sync(id)
        end if

        if( next_output_d(da_sed) ) then
          id = nint(da_sed(4))
          call shy_write_scalar_record2d(id,dtime,891,real(bh))
          call shy_write_scalar_record2d(id,dtime,892,gskm)
          call shy_write_scalar_record2d(id,dtime,893,tao)
          call shy_write_scalar_record2d(id,dtime,894,percc)
          !call shy_write_scalar_record2d(id,dtime,804,hvk)
          call shy_write_scalar_record2d(id,dtime,895,totbed)
	  call shy_sync(id)
        end if

	end subroutine wr_sed_output

!******************************************************************
! Routine to compute variables to be passes to the ecological model ACQUABC
! The following variables are included in the module mod_sediment:
!  - nlbdim: maximum number of bed levels
!  - tmsus: total suspended sediment load [kg/s] (> 0 -erosion, < 0 deposition)
!  - tcn  : total concentration in WC [kg/m3] (computed outside this routine)
!  - bdens: dry bulk density of bottom sediments [kg/m**3]
!  - bleve: depth below sediment surface of sediments [m]

	subroutine get_sedim_prop(nscls,bedn,sload)

        use basin, only : nkn
	use mod_sediment

	implicit none

        integer, intent(in)	     :: nscls   !number grainsize classes
        double precision, intent(in) :: bedn(nlbdim,3,nkn) !bed characteristics in 3 column table
        		                    		   ! (1) depth below sediment surface (m)
	                                                   ! (2) critical erosion stress (Pa)
               		                                   ! (3) dry bulk density (kg/m**3)
        double precision, intent(in) :: sload(nkn,nscls)   !suspended sediment load [kg/s]

	integer 	:: k,is,lb
	
!       -------------------------------------------------------------------
!       Loop over nodes
!       -------------------------------------------------------------------

	do k = 1,nkn

!       -------------------------------------------------------------------
!       Compute total suspended sediment load
!       -------------------------------------------------------------------

	  tmsus(k) = 0.
	  do is = 1,nscls
            tmsus(k) = tmsus(k) + sload(k,is)
	  end do
    
!       -------------------------------------------------------------------
!       Compute level structure, density and porosity
!       -------------------------------------------------------------------

	  do lb = 1,nlbdim
	    bleve(lb,k) = bedn(lb,1,k)
	    bdens(lb,k) = bedn(lb,3,k)
	  end do

!       -------------------------------------------------------------------
!       End node loop
!       -------------------------------------------------------------------

	end do

	end subroutine get_sedim_prop

!******************************************************************
! Set no erosion-deposition on area with rocks (krocks = 1)

	subroutine sed_on_rocks(nscls,krocks,sedx,sedy,sload,sflx)

        use basin, only : nkn

	implicit none

        integer, intent(in)	        :: nscls   	   !number grainsize classes
        integer, intent(in)	        :: krocks(nkn)	   !node index with rocks
        double precision, intent(inout) :: sedx(nscls,nkn) !bedload transport in x direction [kg/ms]
        double precision, intent(inout) :: sedy(nscls,nkn) !bedload transport in y direction [kg/ms]
        double precision, intent(inout) :: sload(nkn,nscls)!suspended sediment load [kg/s]
        double precision, intent(inout) :: sflx(nscls,nkn) !flux of suspend sediment [m3/m2]

	integer 	:: k
	
	do k = 1,nkn

	  if ( krocks(k) == 1 ) then
  	    sedx(:,k)  = 0.d0
  	    sedy(:,k)  = 0.d0
  	    sload(k,:) = 0.d0
  	    sflx(:,k)  = 0.d0
	  end if

	end do

	end subroutine sed_on_rocks

!******************************************************************
