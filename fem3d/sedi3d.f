! ***********************************
! ------ SUBROUTINE SEDI3D ---------
! ***********************************
!
! This routine manages the sediment transport computation.
! It calls the subrourine SEDTRANS05 which computes
! the sediment transport rate and the suspended sediment
! concentration in case of current or current and wave
! (see SUBWAVES routine). The routine SCAL3SH handles the
! avdection and diffusion of suspended sediment.
!
C 16/11/2004 Christian Ferrarin ISMAR-CNR
C
C Revision Log :
C
C Mar, 2005     ccf     (sedi3d_f1.f) coming from sedi3d_e3_multi7.f
C                       	new cohesive routine
C Mar, 2005     ccf     (sedi3d_f2.f) add mixing thickness
C Mar, 2005     ccf     (sedi3d_f3.f) merge layer 1 and 2 if bedn(1) < bmix
C Mar, 2005     ccf     (sedi3d_f4.f) update element depth
C Mar, 2005     ccf     (sedi3d_f5.f) get viscosity from new routine
C Mar, 2005     ccf     (sedi3d_f6.f) active layer = bottom roughness heigh
C Apr, 2005     ccf     (sedi3d_f7.f) style changes, create sedt05
C Apr, 2005     ccf     (sedi3d_f8.f) initialize percbd from file
C Apr, 2005     ccf     (sedi3d_f8.f) change rhos per cohesive sed
C Jun, 2005     ccf     (sedi3d_f13.f) adapt to 3D, bottom layer and total depth
C                       	bugfix error in computing bedn(1,3), 
C				dimension in tuek
C Jun, 2005     ccf     (sedi3d_f15.f) change deposition density,
C				add consolidation
C Jun, 2005     ccf     (sedi3d_f16.f) change units, bug fix in TCONC1
C Jun, 2005     ccf     (sedi3d_f17.f) bug fix in upedepth
C Jun, 2005     ccf     (sedi3d_f18.f) bug fix in checkbed, 
C				adapt to sedtrans05-h5.f
C Jul, 2005     ccf     (sedi3d_f20.f) adapt to sedtrans05-h6.f
C Jul, 2005     ccf     (sedi3d_f21.f) bug fix in updepth (as subdry.f)
C Jul, 2005     ccf     (sedi3d_f22.f) bug fix in bedload
C Aug, 2005     ccf     (sedi3d_f23.f) eliminate ripple variables
C Aug, 2005     ccf     (sedi3d_f24.f) bed slope contribution. adjust updepth
C Aug, 2005     ccf     (sedi3d_f25.f) change deposition characteristics
C Sep, 2005     ccf     (sedi3d_f26.f) bug fix in bedman,
C				change boundary for cohesive
C Sep, 2005     ccf     (sedi3d_f27.f) separete erosion and deposition in bedman
C                       	number of layer computed each time
C Nov, 2005     ccf     (sedi3d_f28f) adapt for 3D version
C Nov, 2005     ccf     (sedi3d_f29f) bed slope on threshold, bug fix in updepth
C Nov, 2005     ccf     (sedi3d_f30f) bug fix in bedslope and in getmd
C Nov, 2005     ccf     (sedi3d_f31f) bed slope by gradients
C Nov, 2005     ccf     (sedi3d_f32f) last layer not smaller than 0.1 m
C Nov, 2005     ccf     (sedi3d_f33f) compute vertical mixing coeffcients
C Jan, 2006     ccf     (sedi3d_f34f) bug fix in upedepth and blimit
C Feb, 2006     ccf     (sedi3d_f35f) new suspco routine
C Feb, 2006     ccf     (sedi3d_f36f) adapt to Van Rijn (C0) in sedtrans05-h8.f
C Feb, 2006     ccf     (sedi3d_f37f) correct bugs in checkbed and other things
C May, 2006     ccf     (sedi3d_f38f) correct bugs in line 1119. Introduce KCOES
C May, 2006     ccf     (sedi3d_f39f) bugs nonco. non used edr in cohse.
C                       	no limit percbd. pers(1)>0.1 in line 1120
C May, 2006     ccf     (sedi3d_f40f) limit BEDCHA(1,2) to 0.1. 
C				bugfix in suspco, better conc in cohse
C Jun, 2006     ccf     (sedi3d_f41f) no transport in case of depth< 0.1m
C Jun, 2006     ccf     (sedi3d_f42f) read constants from file
C Jun, 2006     ccf     (sedi3d_f43f) add limcoh
C Jun, 2006     ccf     (sedi3d_f44f) smooth bed elevation change, get_timestep
C Jul, 2006     ccf     (sedi3d_f45f) read angle of repose, 
C				limit shear velocity, write bathymetry
c 11.04.2008    ggu&ccf treatment of boundaries slightly changed
c 16.04.2008    ggu&ccf bugfix calling lin (must pass double precision 0.)
c 22.04.2008    ggu     advection parallelized
c 27.04.2008    ccf     new check, other changes
c 28.04.2008    ggu     call to nrdpar, nrdnls read with double precision
c 29.04.2008    ccf     bug fix in depbed
c 20.10.2008    ccf     add routine for computing the active layer
c 30.04.2009    ccf     add adjtime for initialization to reset bh
c 03.10.2009    ccf     compute transport only for GD
c 13.04.2010    ccf     introduce SSC effect on density
c 27.07.2010    ccf     settling velocity for cohesive sedimen in function of space
c 20.06.2011    ccf     load for erosion deposition for both cohesive and non-cohesive
c 20.06.2011	ccf	deleted suspco routine
c 19.01.2015    ccf     ia_out1/2 introduced
c 10.02.2015    ggu     new read for constants
c
!****************************************************************************

      subroutine sedi(it,dt)

	use mod_depth
	use mod_diff_visc_fric
	use levels
	use basin, only : nkn,nel,ngr,mbw

      implicit none

      integer it			!time in seconds
      real dt				!time step in seconds

      include 'param.h'
      include 'sed_param.h'

! -------------------------------------------------------------
! fem variables
! -------------------------------------------------------------


      real sedkpar,difmol



! -------------------------------------------------------------
! local sediment variables
! -------------------------------------------------------------

      integer icall
      integer isedi			!call parameter
      integer nscls			!number of grainsize classes
      integer nbed			!initial number of bed layers
      integer is,k,l			!counters
      integer nintp

      double precision gs(nsdim)	!SEDIMENT GRAIN DIAMETER (M)
      double precision bdh(nkndim)      !total elevation change [>0depo,<0ero]
      real bh(nkndim)		    	!bottom height variation [m]
      real tcn(nlvdim,nkndim) 		!total sediment concentration [kg/m3]
      real gskm(nkndim)			!AVERAGE SEDIMENT GRAIN DIAMETER ON NODES (M)
      real sbound(nsdim)                !boundary vector [kg/m3]
      real scn(nlvdim,nkndim,nsdim)     !suspended sediment conc (kg/m3)
      real eps(0:nlvdim,nkndim,nsdim)	!vertical mixing coefficient
      real thick			!initial thickness [m]
      real conref			!initial concentration [kg/m3]
      real wsink(0:nlvdim,nkndim,nsdim)	!settling velocity for suspended sediment
      real gdx(nkndim),gdy(nkndim)	!slope gradients
      real tao(nkndim)                  !wave-current shear stress
	real v1v(nkn)
      double precision riph(nkndim)     !ripple height [m]
      double precision ripl(nkndim)     !ripple length [m]
      double precision timedr		!Duration of call to Sedtrans (s)
      double precision sflx(nsdim,nkndim)	!flux of suspend sediment [m3/m2]
      double precision bflx(nsdim,nkndim)	!flux of bedload sediment [m3/m2]
      double precision sedx(nsdim,nkndim)	!x bedload component [kg/ms]
      double precision sedy(nsdim,nkndim)       !y bedload component [kg/ms]
      integer adjtime				!time for initialization [s]
      integer nsclsc				!number of grainsize classes for adv-dif

      real sedpa(7)                     !sediment parameter vector
      common /sedpa/sedpa
      real gsc(nsdim)                   !grainsize class
      common /gsc/gsc

      integer idsedi(nbcdim)            !information on boundaries
      save idsedi

      integer itanf,nvar
      double precision dtime0,dtime

      double precision hzoff
      real salref,temref		!salinity [psu] and temperature [C]

      character*10 what
      real tsec
      real totbed(nkndim)		!total bedload transport (kg/ms)
      double precision percbd(nlbdim,nsdim,nkndim)        !fraction of sediment [0,1]
      double precision bedn(nlbdim,3,nkndim)	!bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)

      real fact
      parameter ( fact = 1. )
      real load(nlvdim,nkndim)
      integer itmsed,idtsed 		!output parameter
      integer ia_out1(4),ia_out2(4)

! function
      integer iround
      real getpar
      real sload(nkndim,nsdim)	  !suspended sediment load [kg/s]
      logical has_output,next_output

! save and data
      save ia_out1, ia_out2
      save what
      save itmsed,idtsed
      save sedkpar,difmol
      save gs,sbound
      save nscls,adjtime,nsclsc
      save hzoff
      save tcn,bh
      save riph,ripl
      save scn,eps
      save wsink,sload
      save sflx,sedx,sedy
      save gdx,gdy
      save tao,gskm
      save bedn,percbd
      save salref,temref
      save /gsc/
      save /sedpa/

      save icall
      data icall /0/

! ----------------------------------------------------------
! Documentation
! ----------------------------------------------------------
!
! nsdim		dimension of grain size classes arrays
! nscls		number of different grain size classes
! nbcc		number of cohesive grain size classes (automatic)
! gs(i)		value for single grain size classes
!		only first entry of gs may be cohesive
!		if this is the fact, then nbcc has a default value of 7
!		e.g., 7 cohesive grain size sub-classes will be used
!		if gs(i) is non cohesive, no cohesive classes will be used
!		nbcc = 0
!
! sedx/y	-> maybe invert indices
!
! ----------------------------------------------------------
! Initialization
! ----------------------------------------------------------
! This section is called only the first time step when ICALL = 0

        if( icall .le. -1 ) return

        isedi = nint(sedpa(1))
	if( isedi .le. 0 ) then
	  icall = -1
	  return
	end if

        if( icall .eq. 0 ) then

          if( it .lt. itmsed ) return
          icall = 1

!         --------------------------------------------------
!	  Initialize state variables and constants
!         --------------------------------------------------

          call readsedconst		!initialize constants

          conref = sedpa(2)		!initial concentration
          sedkpar = sedpa(3)		!diffusion parameter
          nscls = nint(sedpa(4))	!number of grain size classes
          adjtime = nint(sedpa(7))	!time for initialization
          difmol = getpar('difmol')	!molecular vertical diffusivity [m**2/s]
	  hzoff = getpar('hzon') + 0.10d0
          nbed = 6			!initial number of bed layer
          thick = 0.0			!initial bed layer thickness [m]
          temref = getpar('temref')
          salref = getpar('salref')

          do is = 1,nsdim
           gs(is) = 0.
          end do

	  nbcc = 0
          do is = 1,nscls
           gs(is) = gsc(is)
	   sbound(is) = 0.
	   if (gs(is) .lt. limcoh) nbcc = nbcc + 1
          end do
	  nsclsc = nscls
          if (gs(nscls) .ge. 0.10D0) nsclsc = nscls - 1

          do k = 1,nkn
            bh(k) = 0.
            gdx(k) = 0.
            gdy(k) = 0.
	    riph(k) = 0.d0
	    ripl(k) = 0.d0
	    totbed(k) = 0.
            do l=1,nlvdim
              tcn(l,k) = 0.			!total suspended sediment
	      load(l,k) = 0.
            end do
            do is=1,nsdim
              sedx(is,k) = 0.d0		!bedload
              sedy(is,k) = 0.d0		!bedload
              sload(k,is) = 0.		!suspended sediment load
              do l=1,nlvdim
                scn(l,k,is) = conref	!suspended sediment (sand)
                wsink(l,k,is) = 0.	!settling velocity
                eps(l,k,is) = 0.	!vertical diffusion for sand
              end do
              eps(0,k,is) = 0.		!vertical diffusion
              wsink(0,k,is) = 0.	!settling velocity
	    end do
	  end do

!         --------------------------------------------------
!         Initialize bed configuration
!         --------------------------------------------------

          call inibed(gs,nscls,nbed,thick,gskm,percbd,bedn)

!         --------------------------------------------------
!         Set boundary conditions for all state variables
!         --------------------------------------------------

          call get_first_time(itanf)
          dtime0 = itanf
          nintp = 2
          nvar = nscls
          what = 'sedt'
          call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +                          ,sbound,idsedi)

!         --------------------------------------------------
!	  Initialize output
!         --------------------------------------------------

          call init_output('itmsed','idtsed',ia_out1)
          call init_output('itmsed','idtsed',ia_out2)

          if( has_output(ia_out1) ) then
             call open_scalar_file(ia_out1,1,5,'sed')
          end if
          if( has_output(ia_out2) ) then
             call open_scalar_file(ia_out2,nlv,1,'sco')
          end if

          write(6,*) 'sediment model initialized...'

	endif

! -------------------------------------------------------------------
! Normal call
! -------------------------------------------------------------------

        TIMEDR = dt
        tsec = it

!       -------------------------------------------------------------
!       Reset bottom height variation if it = adjtime
!       -------------------------------------------------------------

	call resetsedi(it,adjtime,bh,scn)

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

!       -------------------------------------------------------------------
!       Compute bedload transport
!       -------------------------------------------------------------------

        call bedload(nscls,sedx,sedy,dt,bh,bflx)

!       -------------------------------------------------------------------
!       Transport and diffusion for each sediment class
!       -------------------------------------------------------------------

	dtime = it
	call bnds_read_new(what,idsedi,dtime)

!$OMP PARALLEL PRIVATE(is)
!$OMP DO SCHEDULE(DYNAMIC)

        do is = 1,nsclsc
	  call load3d(sload(1,is),nkn,nlvdi,ilhkv,load)
          call scal_adv_fact(what,is,fact
     +                      ,scn(1,1,is),idsedi
     +                      ,sedkpar,fact,wsink(0,1,is),fact,load
     +                      ,difhv,eps(0,1,is),difmol)
        end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

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
!       Compute total suspended concentration at nodes and update water density
!       -------------------------------------------------------------------

        call totcon(nscls,scn,tcn)

!       -------------------------------------------------------------------
!       Check mass conservation - only for closed system
!       -------------------------------------------------------------------

        !call sedcons(it,tcn,bh,bedn)

!       -------------------------------------------------------------------
!       Write of results (files SED and SCO)
!       -------------------------------------------------------------------

        if( next_output(ia_out1) ) then
          call write_scalar_file(ia_out1,80,1,bh)
          call write_scalar_file(ia_out1,81,1,gskm)
          call write_scalar_file(ia_out1,82,1,tao)
          call write_scalar_file(ia_out1,83,1,hkv)
          call write_scalar_file(ia_out1,84,1,totbed)
        end if

        if( next_output(ia_out2) ) then
          call write_scalar_file(ia_out2,85,nlv,tcn)
        end if

! -------------------------------------------------------------------
! End of routine
! -------------------------------------------------------------------

	end

! ********************************************************************
! SUBROUTINE READGS
! This subroutine reads the simulation sediment parameter from the .str file: 
! grainsize classes, the initial sediment distribution (percentage) and the 
! critical erosion threshol value from the .str file

        subroutine readsed

        implicit none

        include 'sed_param.h'

        character*80 name		!name of item read
        character*80 text		!text of item if string
        real value			!value of item if numeric
	double precision dvalue		!value of item if numeric
        integer iweich
        integer nrdpar
        real getpar
        integer i
	real tauin

! --- output variables
        real sedpa(7)			!sediment parameter vector
        common /sedpa/sedpa
        integer nrs			!number of grainsize classes
        real gsc(nsdim)			!grainsize class
        common /gsc/gsc
        integer npi			!
        integer nps
        real prin(0:nsdim,nsdim)        !initial percetage [0,1]
        common /prin/prin
        integer nrst			!number of sediment tipe
        real tue(nsdim)			!initial erosion threshold
        common /tue/tue        
        save /gsc/
        save /sedpa/

c DOCS  START   P_sediment
c
c The following parameters activate the sediment transport module 
c and define the sediment grainsize classes to the simulated.
!c |sedtr|	Sediment transport module section name.

        call sctpar('sedtr')             !sets default section
        call sctfnm('sedtr')

c |isedi|	Flag if the computation on the sediment is done:
c		\begin{description}
c		\item[0] Do nothing (default)
c		\item[1] Compute sediment transport
c		\end{description}

        call addpar('isedi',0.)

c |idtsed|, |itmsed|	Time step and start time for writing to files sed e sco,
c			the files containing sediment variables and suspended
c			sediment concentration.

        call addpar('idtsed',0.)
        call addpar('itmsed',-1.)

c |sedgrs|	Sediment grainsize class vector [mm]. Values has be 
c		ordered from the finest to the more coarse. \\
c		|example: sedgrs = 0.1 0.2 0.3 0.4|

        call addpar('sedgrs',0.)

c |sedref|	Initial sediment reference concentration [kg/m3]
c		(Default 0).

        call addpar('sedref',0.)

c |sedhpar|	Sediment diffusion coefficient (Default 0).

        call addpar('sedhpar',0.)

c |adjtime|	Time for sediment initialization [s]. The sediment model take a
c		initialization time in which the system goes to a quasi steady
c		state. When t = |adjtime| the bed evolution change in the output
c		is reseted. Keep in mind that |adjtime| has to be chosen case by 
c		case in function of the morphology and the parameters used for 
c		the simulation (Default 0).

        call addpar('adjtime',-99999999.)

c |percin|	Initial sediment distribution [0,1] for each 
c		grainsize class. The sum of percin must be equal to 1. \\
c		|example: percin = 0.25 0.25 0.25 0.25| \\
c		If percin is not selected the model impose equal
c		percentage for each grainsize class (percin = 1/nrs).
c		In case of spatial differentiation of the sediment
c		distribution set a number of percin equal to the number
c		of grainsize classes per the number of area types. \\
c		|example: percin = 0.25 0.25 0.25 0.25  \\
c			 	   0.20 0.20 0.30 0.30  \\
c				   0.45 0.15 0.15 0.15| \\

        call addpar('percin',0.)

c |tauin|	Initial dry density or TAUCE. In function of the value: \\
c		0-50 : critical erosion stress (Pa) \\
c		\textgreater 50  : dry bulk density of the surface (kg/m**3). \\
c		In case of spatial differentiation set a number of tauin
c		equal to the number of area type. \\
c		|example: tauin = 0.9 1.4 2.5 1.1|

        call addpar('tauin',0.)

c |sedp|	File containing spatially varying initial sediment distribution. 
c		Values are in percentage of each class.

        call addfnm('sedp',' ')

c |sedt|	File containing spatially varying initial critical erosion 
c		stress (Pa) or dry bulk density (kg/m3).

        call addfnm('sedt',' ')

c |sedcon|	File containing the additional constants used in sediment 
c		model. These parameters are usually set to the indicated 
c		default values, but can be customized for each sediment 
c		transport simulation. The full parameter list together with 
c		their default value and brief description is reported in 
c		Table \ref{tab:table_sedcon}. Most of the parameter, 
c		especially the ones for the cohesive sediments,	have been 
c		calibrated for the Venice Lagoon. For more information about 
c		these parameters please refer to Neumeier et 
c		al. \cite{urs:sedtrans05} and Ferrarin et al. 
c		\cite{ferrarin:morpho08}.

        call addfnm('sedcon',' ')


c DOCS  FILENAME        Additional parameters
c The full parameter list is reported in Table \ref{tab:table_sedcon}.
c An example of the settings for the |sedcon| file is given in
c \Fig\figref{turbulence}. Please note that is not necessary to
c define all parameters. If not defined the default value is imposed.
c
c \begin{figure}[ht]
c \begin{verbatim}
c IOPT = 3
c SURFPOR = 4
c DOCOMPACT = 1
c \end{verbatim}
c \caption{Example of the secon file.}
c \label{fig:sedcon}
c \end{figure}

c \begin{table}[ht]
c \caption{Additional parameter for the sediment transport model to be set 
c in the sedcon file.}
c \begin{tabular}{lcl} \hline
c Name & Default value & Description \\ \hline
c CSULVA  & 159.4 & Coefficient for the solid transmitted stress by Ulva \\
c TMULVA  & 1.054d-3 & Threshold of motion of Ulva (Pa) \\
c TRULVA  & 0.0013 & Threshold of full resuspension of Ulva (Pa) \\
c E0      & 1.95d-5 & Minimum erosion rate \\
c RKERO   & 5.88 & Erosion proportionality coefficient \\
c WSCLAY  & 5.0 & Primary median Ws class (in the range 1:NBCONC) \\
c CDISRUPT& 0.001 & Constant for turbulent floc disruption during erosion \\
c CLIM1   & 0.1 & Lower limit for flocculation (kg/m3) \\
c CLIM2   & 2.0 & Limit between simple and complex flocculation (kg/m3) \\
c KFLOC   & 0.001 & Constant K for flocculation equation \\
c MFLOC   & 1.0 & Constant M for flocculation equation \\
c RHOCLAY & 2600.0 &  Density of clay mineral \\
c CTAUDEP & 1.0 & Scaling factor for TAUCD \\
c PRS     & 0.0 & Resuspension probability (range 0-1) \\
c RHOMUD  & 50.0 & Density of the freshly deposited mud \\
c DPROFA  & 470.0 & Constants for density profile \\
c DPROFB  & 150.0 & A : final deep density \\
c DPROFC  & 0.015 & Define the shape (in conjunction with B and C) \\
c DPROFD  & 0.0 & Aux parameter for density profile \\
c DPROFE  & 0.0 & Aux parameter for density profile \\
c CONSOA  & 1d-5 & time constant of consolidation \\
c TEROA   & 6d-10 & Constant for erosion threshold from density \\
c TEROB   & 3.0 & Aux parameter for erosion threshold from density \\
c TEROC   & 3.47 & Aux parameter for erosion threshold from density \\
c TEROD   & -1.915 & Aux parameter for erosion threshold from density \\
c KCOES   & 0.15 & Fraction of mud for sediment to be cohesive \\
c CDRAGRED& -0.0893 & Constant for the drag reduction formula \\
c Z0COH   & 2.0D-4 & Bed roughness length for cohesive sediments \\
c FCWCOH  & 2.2D-3 & Friction factor for cohesive sediments \\
c LIMCOH  & 0.063 & Limit of cohesive sediment grainsize [mm] \\
c SMOOTH  & 1.0 & Smoothing factor for morphodynamic \\
c ANGREP  & 32.0 & Angle of repose \\
c IOPT    & 5 & Sediment bedload transport formula option number \\
c MORPHO  & 1.0 & Morphological acceleration factor \\
c RHOSED  & 2650.0 & Sediment grain density \\
c POROS   & 0.4 & Bed porosity [0,1] \\
c SURFPOR & 0.6 & Bed porosity of freshly deposited sand [0,1] \\
c DOCOMPACT& 0.0 & If not zero, call COMPACT routine \\ \hline
c \end{tabular}
c \label{tab:table_sedcon}
c \end{table}
c
c DOCS  END

!       --------------------------------------------------
!       Initialize variables
!       --------------------------------------------------

	tauin = getpar('tauin')
        do i = 1,nsdim
          tue(i) = tauin
        end do

        npi = 1
        nps = 0
        nrs = 0
        nrst = 0

!       --------------------------------------------------
!       Starts the reading loop
!       --------------------------------------------------

        iweich = 1
        do while(iweich.ne.0)
          iweich=nrdpar('sedtr',name,dvalue,text)
	  value = real(dvalue)

!         --------------------------------------------------
!         Reads sediment grainsize vector with dimension nrs
!	  and convert it from mm to m
!         --------------------------------------------------
          if( name .eq. 'sedgrs' ) then
              nrs = nrs+1
              if( nrs .gt. nsdim ) goto 35
	      if (value .ge. 100d0) goto 36
              gsc(nrs)=value*0.001
          end if

!         --------------------------------------------------
!         Reads sediment percentage matrix for each type
!         --------------------------------------------------
          if( name .eq. 'percin' ) then
              nps = nps+1
              if( nps .gt. nsdim ) goto 35
              if( npi .gt. nsdim ) goto 35
              if( nps .gt. nrs ) then
                nps = 1
                npi = npi + 1
              end if
              prin(npi,nps)=value
          end if

!         --------------------------------------------------
!         Reads initial erosion threshold vector
!         --------------------------------------------------
          if( name .eq. 'tauin' ) then
              nrst=nrst+1
              if( nrst .gt. nsdim ) goto 35
              tue(nrst)=value
          end if
        end do

        if ( nrs .eq. 0 ) goto 30

!       --------------------------------------------------
!       Case with only one sediment class
!       --------------------------------------------------
        if( nrs .eq. 1 ) then
          npi = 1
          nps = 1
          prin(npi,nps) = 1
        end if

!       --------------------------------------------------
!       Case with no initial distribution selected
!	impose the same percentage for each class
!       --------------------------------------------------
        if ( nps .eq. 0 ) then
          do i = 1,nrs
            prin(npi,i) = 1./nrs
          end do
          nps = nrs
        end if

!       --------------------------------------------------
!       Stores parameters in variable sedpa()
!       --------------------------------------------------

        sedpa(1) = getpar('isedi')
        sedpa(2) = getpar('sedref')
        sedpa(3) = getpar('sedhpar')
        sedpa(4) = nrs
        sedpa(5) = npi
        sedpa(6) = nrst
        sedpa(7) = getpar('adjtime')

        return

  30    continue
        write(6,*) 'No sediment class selected'
        write(6,*) 'in the .str file'
        stop 'error stop : readsed'

  35    continue
        write(6,*) 'Dimension error for nsdim'
        write(6,*) 'nsdim  :',nsdim
        write(6,*) 'nscls :',nrs
        stop 'error stop : readsed'

  36    continue
        write(6,*) 'First grainsize class > 10 cm'
        write(6,*) 'grainsize :',gsc(1)
        stop 'error stop : readsed'

        end

! ********************************************************************
! SUBROUTINE READCONST
! This subroutine read constants from file sedcon and save them in the
! common block. The format of the file is: variable = value

        subroutine readsedconst

        implicit none

        include 'sed_param.h'

        integer nrdnxt,iw,ioff
        integer ifileo
        integer iunit
        real value
	double precision dvalue
        integer nc		!number of variables
        parameter(nc = 38)
        character*80 file,name,text,line,cname(nc)
        double precision cvalue(nc)
        integer is,j

	iw = 0
	iunit = 0
	value = 0.0
        file = ''

        G = 9.81d0
        PII = 2.*ASIN(1d0)
	NDISTR = 5
	MEDDISTR = 3
	do is = 1,nsdim
	  WSI(IS) = 0.d0
	  DISTR(IS) = 0.d0
	end do

!       --------------------------------------------------------
!       Assigns constants to cname and default values to cvalue
!       --------------------------------------------------------

        cname(1) = 'CSULVA'	!coefficient for the solid transmitted stress by Ulva
        cvalue(1) = 159.4d0
        cname(2) = 'TMULVA'	!threshold of motion of Ulva (Pa)
        cvalue(2) = 1.054d-3
        cname(3) = 'TRULVA'	!threshold of full resuspension of Ulva (Pa)
        cvalue(3) = 0.0013d0
        cname(4) = 'E0'		!Minimum erosion rate
        cvalue(4) = 1.95d-5
        cname(5) = 'RKERO'	!Erosion proportionality coefficient
        cvalue(5) = 5.88d0
        cname(6) = 'WSCLAY'	!primary median Ws class (must be in the range 1:NBCONC)
        cvalue(6) = 5
        cname(7) = 'CDISRUPT'	!constant for turbulent floc disruption during erosion
        cvalue(7) = 0.001d0
        cname(8) = 'CLIM1'	!lower limit for flocculation (kg/m**3)
        cvalue(8) = 0.1d0
        cname(9) = 'CLIM2'	!limit between simple and complex flocculation equation (kg/m**3)
        cvalue(9) = 2d0
        cname(10) = 'KFLOC'   	!constant K for flocculation equation
        cvalue(10) = 0.001d0
        cname(11) = 'MFLOC'   	!constant M for flocculation equation
        cvalue(11) = 1d0
        cname(12) = 'RHOCLAY' 	!density of clay mineral
        cvalue(12) = 2600d0
        cname(13) = 'CTAUDEP' 	!scaling factor for TAUCD
        cvalue(13) = 1d0
        cname(14) = 'PRS'     	!Resuspension probability (range 0-1)
        cvalue(14) = 0d0
        cname(15) = 'RHOMUD'  	!Density of the freshly deposited mud
        cvalue(15) = 50d0	!50d0
        cname(16) = 'DPROFA'	!Constants for density profile
        cvalue(16) = 470d0	!FINAL(I) = A - B * EXP(-C*ABOVE(I)) - D*EXP(-E*ABOVE(I))
        cname(17) = 'DPROFB'	!A : final deep density
        cvalue(17) = 150d0	!A-B-D : final surface density
        cname(18) = 'DPROFC'	!define the shape (in conjunction with B and C)
        cvalue(18) = 0.015d0
        cname(19) = 'DPROFD'
        cvalue(19) = 0d0
        cname(20) = 'DPROFE'
        cvalue(20) = 0d0
        cname(21) = 'CONSOA'	!time constant of consolidation
        cvalue(21) = 1d-5
        cname(22) = 'TEROA'	!Constants for erosion threshold from density and overlaying mass
        cvalue(22) = 6d-10
        cname(23) = 'TEROB'
        cvalue(23) = 3d0
        cname(24) = 'TEROC'
        cvalue(24) = 3.47d0
        cname(25) = 'TEROD'
        cvalue(25) = -1.915d0
        cname(26) = 'KCOES'	!Fraction of mud for sediment to be cohesive
        cvalue(26) = 0.15d0
        cname(27) = 'CDRAGRED'	!constant for the drag reduction formula
        cvalue(27) = -0.0893d0
        cname(28) = 'Z0COH'	!BED ROUGHNESS LENGHT FOR COHESIVE SEDIMENTS
        cvalue(28) = 2.0D-4
        cname(29) = 'FCWCOH'	!FRICTION FACTOR FOR COHESIVE SEDIMENTS
        cvalue(29) = 2.2D-3
        cname(30) = 'LIMCOH'	!Limit of cohesive sediment grainsize
        cvalue(30) = 0.000063d0
        cname(31) = 'SMOOTH'	!smoothing factor for morphodynamic
        cvalue(31) = 1.0d0
        cname(32) = 'ANGREP'	!angle of repose
        cvalue(32) = 32.0d0
        cname(33) = 'IOPT'	!SEDIMENT TRANSPORT FORMULA OPTION NUMBER
        cvalue(33) = 5
        cname(34) = 'MORPHO'	!Morphological acceleration factor
        cvalue(34) = 1.d0
        cname(35) = 'RHOSED'    !sediment grain density
        cvalue(35) = 2650.d0
        cname(36) = 'POROS'	!bed porosity [0,1]
        cvalue(36) = 0.4d0
        cname(37) = 'SURFPOR'	!bed porosity of freshly deposited sand [0,1]
        cvalue(37) = 0.6d0
        cname(38) = 'DOCOMPACT'	!if not zero, call COMPACT routine
        cvalue(38) = 0

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
        
        CSULVA = cvalue(1)
        TMULVA = cvalue(2)
        TRULVA = cvalue(3)
        E0 = cvalue(4)
        RKERO = cvalue(5)
        WSCLAY = int(cvalue(6))
        CDISRUPT = cvalue(7)
        CLIM1 = cvalue(8)
        CLIM2 =  cvalue(9)
        KFLOC = cvalue(10)
        MFLOC = cvalue(11)
        RHOCLAY = cvalue(12)
        CTAUDEP = cvalue(13)
        PRS = cvalue(14)
        RHOMUD = cvalue(15)
        DPROFA = cvalue(16)
        DPROFB = cvalue(17)
        DPROFC = cvalue(18)
        DPROFD = cvalue(19)
        DPROFE = cvalue(20)
        CONSOA = cvalue(21)
        TEROA = cvalue(22)
        TEROB = cvalue(23)
        TEROC = cvalue(24)
        TEROD = cvalue(25)
        KCOES = cvalue(26)
        CDRAGRED = cvalue(27)
        Z0COH = cvalue(28)
        FCWCOH = cvalue(29) 
        LIMCOH = cvalue(30) 
        SMOOTH = cvalue(31) 
        ANGREP = cvalue(32) / (45./atan (1.))	!from deg to rad
        IOPT = int(cvalue(33))
        MORPHO = cvalue(34)
        RHOSED = cvalue(35)
        POROS = cvalue(36)
        SURFPOR = cvalue(37)
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
        end    

! ********************************************************************
! SUBROUTINE INIBED
! This subroutine set the initial bed conformation in function of
! sediment grainsize classes and the initial percentage defined in
! the .str file. The distribution is the same in each bed level (no
! vertical stratification). 

        subroutine inibed(gs,nscls,nbed,thick,gskm,percbd,bedn)

	use basin

        implicit none

        include 'param.h'
        include 'sed_param.h'


        real sedpa(7)	                        !sediment parameter vector
        common /sedpa/sedpa
        real prin(0:nsdim,nsdim)                !initial percetage [0,1]
        common /prin/prin
        real tue(nsdim)             		!initial erosion threshold/density
        common /tue/tue
        double precision gs(nsdim)              !sediment grain diameter (m)
        double precision percbd(nlbdim,nsdim,nkndim)        !fraction of sediment [0,1]
        double precision bedn(nlbdim,3,nkndim)  !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nscls				!number of grainsize classes
        integer nbed				!initial number of bed layers
        real thick				!initial thickness [m]
        real gskm(nkndim)			!average grainsize per node [m]
        double precision gsa,gss
	double precision ptot
        double precision rhosa			!initial erosion threshold
        double precision rhossand		!function for density from % sand
        double precision pcoes			!% of fine sediments
        integer itype   	                !element code number
        real pers(nkndim,nsdim)			!percbd initialized from file
        real tuek(nkndim,1)			!tauce initialized from file
        integer npi
        integer k,ib,is,ie,ii

        save /sedpa/

        npi = int(sedpa(5))

!       -------------------------------------------------------------------
!       Initialize sediment distribution in function of bed type
!       -------------------------------------------------------------------

        do ie = 1,nel
          itype = iarv(ie) + 1
          if ( npi.eq.1 ) itype = 1

          do ii = 1,3
            k = nen3v(ii,ie)
            tuek(k,1) = tue(itype)
            do is = 1,nscls
              pers(k,is) = prin(itype,is)
            end do
          end do
        end do

!       -------------------------------------------------------------------
!       Initialize sediment fraction and tuek from external file
!       -------------------------------------------------------------------

        call inic2fil('sedp',pers,nscls)
        call inic2fil('sedt',tuek,1)

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
!       Initialize bed characteristics and compute initial average grainsize
!       -------------------------------------------------------------------

        do k = 1,nkn
          gsa = 0.d0
	  ptot = 0.d0
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
          if(ptot.lt.0.99d0.or.ptot.gt.1.01d0) go to 130
          gskm(k) = real(gsa)


          rhosa = rhossand(pcoes,1.2d0)
          if(tuek(k,1).gt.0.) rhosa = tuek(k,1)

          call bedini(bedn(1,1,k),nlbdim,rhosa)
        end do

	return

 130    continue
        write(6,*) 'Error in computing the sediment fraction: 3'
        write(6,*) 'total percentage:',ptot, 'node:',k
        do is = 1,nscls
          write(6,*) 'classnumber:',is,'percentage:',percbd(1,is,k)
        enddo
        stop 'error stop : percbd'

        end

! *******************************************************
! FUNCTION RHOSSAND
! returns the dry bulk density from sand fraction (Allersma, 1988)

        function rhossand(pcoes,consc)

        implicit none

	include 'sed_param.h'

        double precision pcoes                  !fraction of fine sediments
        double precision consc                  !consolidation coefficient [0-2.4]
        double precision psand                  !sand fraction [0,1]
        double precision rhossand               !density function [kg/m3]

        psand = 1.d0
        if(pcoes.gt.0.d0) then 
          psand = 1.d0 - pcoes
          rhossand = 480.d0*consc + (1300.d0 - 280.d0*consc)*
     $(psand**0.8)
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
	double precision k    	!von karman costant
	parameter(k=0.4d0)

        br = 2.5*real(gd)		!see amos article
        z0 = br/30.
	z0 = max(z0,z0bk)
        if (uz .le. 0.d0) then
	  z = d/2.d0
	else
          ustar = (uz*k*d/z0)/(d/z0*(log(d/z0)-1.d0)+1.d0)
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

        alfa=atan(v/u)*rad

        if( u .eq. 0.0 ) then
          alfa = 0.
          if( v .lt. 0.0 ) alfa = 180.
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
     +			    per,dx,difv,visv,eps)

        implicit none

        include 'param.h'

! --- input variables

	integer dtype			!type of diffusion algoritm
        integer l			!number of levels
        real dep(nlvdim)		!thickness of layer [m]
        double precision h		!total water depth [m]
        double precision wsink		!settling velocity [m/s]
        double precision ustc		!current shear velocity [m/s]
        double precision ustcw		!combined current and wave shear velocity [m/s]
        double precision ub 		!wave orbital velocity [s]
        double precision dcw		!height of the wave-current boundary layer
        double precision ht		!significant wave height [m]
        double precision per		!significant wave persion [s]
        double precision dx		!dimensionless particle diameter
        real visv(0:nlvdim)		!water viscosity [from GOTM] 
        real difv(0:nlvdim)		!water diffusivity [from GOTM] 

! --- output variables

        real eps(0:nlvdim)

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
            eps(m) = 0.
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

            eps(m) = real(sqrt(epsc**2 + epsw**2))

          end do

          eps(0) = eps(1)

	else if (dtype .eq. 2) then

!         ------------------------------------------------------------
!         Use viscosity from GOTM and beta factor from Van Rijn 1993
!	  Beta is function of settling velocity and shear velocity
!         ------------------------------------------------------------

          bet = 1.d0 + 2.d0* (wsink/ustcw)**2
          bet = dmax1(bet,1.d0)
          bet = dmin1(bet,1.5d0)
	  do m = 0,l
	    eps(m) = bet * visv(m)
	  end do

	else if (dtype .eq. 3) then

!         ------------------------------------------------------------
!         Use viscosity from GOTM and beta factor from Villatoro et al 2010
!	  Beta is function of grainsize
!         ------------------------------------------------------------

	  bet = 0.33d0*dx - 0.34d0
	  do m = 0,l
	    eps(m) = bet * visv(m)
	  end do

	else if (dtype .eq. 4) then

!         ------------------------------------------------------------
!         Use diffusivity from GOTM 
!         ------------------------------------------------------------

	  do m = 0,l
	    eps(m) = difv(m)
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
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'
        include 'sed_param.h'

! -------------------------------------------------------------
! local variables
! -------------------------------------------------------------
        integer nscls			!number of grainsize classes
        integer is,k,l			!counters
	integer lmax			!bottom level
        double precision dl		!thickness of bottom layer 
        double precision gs(nsdim)	!SEDIMENT GRAIN DIAMETER (M)
        double precision ws(nsdim)	!SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        real scn(nlvdim,nkndim,nsdim)   !cohesive suspended sediment conc (kg/m3)
        double precision riph(nkndim)	!ripple height [m]
        double precision ripl(nkndim)	!ripple length [m]
        real gskm(nkndim)		!AVERAGE SEDIMENT GRAIN DIAMETER ON NODES (M)
        real eps(0:nlvdim,nkndim,nsdim)	!vertical mixing coefficient
        real u,v			!x and y current components
        real gdx(nkndim),gdy(nkndim)	!slope gradients
        real tao(nkndim)		!wave-current shear stress
        double precision sflx(nsdim,nkndim)     !flux of suspend sediment [m3/m2]
        double precision sedx(nsdim,nkndim)	!x bedload component [kg/ms]
        double precision sedy(nsdim,nkndim)    	!y bedload component [kg/ms]
        double precision percbd(nlbdim,nsdim,nkndim)        !fraction of sediment [0,1]
        double precision bedn(nlbdim,3,nkndim)	!bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        double precision scnd(nsdim)    	!suspended sediment conc (kg/m3)
	double precision hzoff
	real totbed(nkndim)			!total bedload transport (kg/ms)
        real salref,temref			!salinity [psu] and temperature [C]
        real wsink(0:nlvdim,nkndim,nsdim)       !settling velocity for suspended sediment

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
	real sload(nkndim,nsdim)	!suspended sediment load [kg/s]
	real sloads(nsdim)	  	!suspended sediment load [kg/m²s]

	real ukbot(nkndim),vkbot(nkndim)
	real ddl(nkndim)
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
          call zvel(UZ,DL,GD,z0bk(k),Z)		!get Z
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
            scnd(is) = scn(lmax,k,is)
          end do

          call sedt05(k,D,DL,UZ,Z,CDIR,HT,MPER,PPER,WDIR,GD,riph(k),
     $ripl(k),BETA,TEM,SAL,bedn(1,1,k),percbd(1,1,k),AULVA,TIMEDR,
     $nscls,gs,UW,hzoff,scnd,sedx(1,k),sedy(1,k),ws,gdx(k),gdy(k),
     $lmax,eps,tao(k),Z0,sloads,sflx(1,k))

	  z0bk(k) = real(Z0)

          do is = 1,nscls
	    totbed(k) = totbed(k) + real(sqrt(sedx(is,k)**2 +
     $			sedy(is,k)**2)* bedn(1,3,k))
	    sload(k,is) = sloads(is) * area	!convert to kg/s
          end do

!         -------------------------------------------------------------------
!         Compute settling velocity (either for flocs)
!         -------------------------------------------------------------------

	  call set_settling(k,lmax,nscls,ws,scn,wsink)

!       -------------------------------------------------------------------
!       End of node loop
!       -------------------------------------------------------------------

        end do

        end

! ********************************************************************
! SUBROUTINE SEDT05

      subroutine sedt05(k,D,DL,UZ,Z,CDIR,HT,MPER,PPER,WDIR,GD,RHINP,
     $RLINP,BETA,TEM,SAL,BEDCHA,percbd,AULVA,TIMEDR,nscls,gs,UW,hzoff,
     $scn,sedx,sedy,ws,gdx,gdy,lmax,eps,tao,Z0,sload,sflx)

	use mod_depth
	use mod_layer_thickness
	use mod_diff_visc_fric

        implicit none

	include 'param.h'
        include 'sed_param.h'

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
        double precision percbd(nlbdim,nsdim)   !fraction of sediment [0,1]
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
        double precision sedx(nsdim)
        double precision sedy(nsdim)
        double precision sflx(nsdim)    !flux of suspend sediment [m3/m2]
        double precision ws(nsdim)      !SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        real tao			!wave-current shear stress
	real sload(nsdim)	  	!suspended sediment load [kg/m²s]

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

        double precision pers(nsdim)
        double precision scn(nsdim) 	!cohesive suspended sediment conc (kg/m3)
        double precision dxx(nsdim)	!DIMENSIONLESS GRAIN SIZE

        double precision usb(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF BEDLOAD
        double precision uss(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF SUSPENDED
        double precision ust(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF SHEET FLOW

        double precision gs(nsdim)      !SEDIMENT GRAIN DIAMETER (M)
        double precision FALL           !average fall velocity
        double precision DX             !average dimensionless grain diameter
        double precision USTCRB         !average critical velocity
        double precision USTCRS         !average critical velocity
        double precision USTUP          !average sheet flow velocity
        integer lmax,is,k		!counters
        real gdx,gdy			!slope gradients
        double precision alph           !slope effect
        real eps(0:nlvdim,nkndim,nsdim)	!vertical mixing coefficient
	double precision uslim
	double precision pcoes		!% of  fine sediments
	logical cohes
	double precision btype		!type of active layer or bmix thickness
	integer dtype			!type of diffusivity algoritm
	integer wtype			!type of bottom particle velocity algoritm
        integer nsclsc	 	  	!number of grainsize classes without rocks

!       --------------------------------------------------
!       Initialize variables
!       --------------------------------------------------

        do is = 1,nscls
          sedx(is) = 0.d0
          sedy(is) = 0.d0
          pers(is) = percbd(1,is)
	  sload(is) = 0.
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

        if ( D .LT. hzoff ) return

	pcoes = 0.d0
	do is = 1,nbcc
	  pcoes = pcoes + pers(is)
	end do

	nsclsc = nscls
	if (gs(nscls) .ge. 0.1d0) nsclsc = nscls - 1

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

	wtype = 3
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

        do is = 1,nsclsc
          CALL THRESH(VISK,gs(is),RHOS,RHOW,ws(is),usb(is),uss(is),
     +ust(is),dxx(is))
	  WSI(is) = ws(is)
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

        if(ustcw.le.0.d0) ustcw=dmax1(ustc,ustw)
        if(ustcws.le.0.d0) ustcws=dmax1(ustcs,ustws)
        tao = real(rhow * ustcws**2)

!       -------------------------------------------------------------------
!       Calculate vertical diffusion coefficient following (dtype):
!	  -1- Algebric Van Rijn 1993
!	  -2- Viscosity from GOTM and beta from Van Rijn 1993
!	  -3- Viscosity from GOTM and beta from Villatoro et al 2010
!         -4- Diffusivity from GOTM 
!       -------------------------------------------------------------------

	dtype = 2
        do is = 1,nsclsc
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

        call bedact(nscls,bmix,BEDCHA,percbd)
        do is = 1,nscls
          pers(is) = percbd(1,is)
        end do

!       -------------------------------------------------------------------
!       Calculate sediment transport rate for each sediment class in node k
!       -------------------------------------------------------------------

        if( cohes ) then

!         -------------------------------------------------------------------
!         COHESIVE SEDIMENT
!         -------------------------------------------------------------------

          call cohse(BEDCHA,TIMEDR,USTCWS,USTCW,DL,RHOW,
     $AULVA,UB,CDIR,Z0,nscls,scn,ws,usb,uss,ust,gs,dxx,
     $Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,USTBF,
     $VISK,pers,RHOS,sload,sflx)

        else

!         -------------------------------------------------------------------
!         NON-COHESIVE SEDIMENT
!         -------------------------------------------------------------------

          call nonco(BEDCHA,TIMEDR,D,DL,UA,UB,U100,HT,PER,CDIR,
     $RHOW,USTCWS,USTCW,usb,uss,ust,Z0,BETA,RHOS,FCW,WDIR,
     $Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,USTBF,
     $nsclsc,pers,gs,dxx,ws,scn,sedx,sedy,sload,sflx)

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

	implicit none

        include 'sed_param.h'

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
     $AULVA,UB,CDIR,Z0,nsclsc,scn,ws,usb,uss,ust,gs,dxx,
     $Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,USTBF,
     $VISK,pers,RHOS,sload,sflux)

        implicit none

	include 'param.h'
        include 'sed_param.h'

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

        integer nsclsc			!number of grainsize classes without rocks
        double precision scn(nsdim) 	!cohesive suspended sediment conc (kg/m3)
        double precision pers(nsdim)	!fraction of each class
        DOUBLE PRECISION GD		!SEDIMENT GRAIN DIAMETER (M) 
        double precision gs(nsdim)      !SEDIMENT GRAIN DIAMETER (M)
        double precision ws(nsdim)      !SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        double precision usb(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF BEDLOAD
        double precision uss(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF SUSPENDED
        double precision ust(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF SHEET FLOW
        double precision dxx(nsdim)     !DIMENSIONLESS GRAIN SIZE

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
	real sload(nsdim)	  	!suspended sediment load [kg/m2s]
        double precision sflux(nsdim)   !flux of suspend sediment [m3/m2]

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
        integer i,is			!counters
	double precision rhoa
	double precision smax,sdepo
       
!       --------------------------------------------------
!       Initialize variables
!       --------------------------------------------------
 
        emass = 0.d0
	ZS = 0.d0

!       -------------------------------------------------------------------
!       Compute the total SSC before deposition and erosion
!       -------------------------------------------------------------------

        TCONC1 = 0.d0
        do I=1,NBCC
          CONC(I) = scn(i)
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

          do is = 1,nsclsc
	    sload(is) = emass * pers(is) / TIMEDR
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
	  
          do is = 1,nbcc
	    sload(is) = sload(is) - dmass * pers(is) / TIMEDR
	    sflux(is) = sflux(is) + DMASS*2.d0/(RHOMUD+BEDCHA(1,3))
     +                  *pers(is)
	  end do
        end if

!       -------------------------------------------------------------------
!       Deposition for non-cohesive suspended sediments (if present)
!       -------------------------------------------------------------------

        do is = nbcc+1,nsclsc
          CALL PROFL(DL,RHOW,ws(is),UB,CDIR,USTCWS,USTCW,usb(is),
     $uss(is),ust(is),Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,
     $RPLCOEF,USTBF,Z0,gs(is),dxx(is),RHOS,QS,QSDIR,C0,C0A)

          C0AI = C0A * pers(is)
  	  C0AI = MIN(C0AI,scn(is))

	  sdepo = ws(is)*(C0AI - scn(is))
          smax = -scn(is)*DL/TIMEDR
          sdepo = dmax1(sdepo,smax)
          sload(is) = sload(is) + sdepo
	  sflux(is) = sflux(is) - sdepo * timedr / BEDCHA(1,3)
        end do

        end 

! ********************************************************************
! SUBROUTINE NONCO
! This subroutine computes the sediment transport for non-cohesive
! sediments

        subroutine nonco(BEDCHA,TIMEDR,D,DL,UA,UB,U100,HT,PER,CDIR,
     $RHOW,USTCWS,USTCW,usb,uss,ust,Z0,BETA,RHOS,FCW,WDIR,
     $Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,RPLCOEF,USTBF,
     $nsclsc,pers,gs,dxx,ws,scn,sedx,sedy,sload,sflux)

        implicit none

	include 'param.h'
        include 'sed_param.h'

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

        integer nsclsc			!number of grainsize classes withour rocks
        double precision pers(nsdim)	!fraction of sediment [0,1]
        double precision gs(nsdim)      !SEDIMENT GRAIN DIAMETER (M)
        double precision ws(nsdim)      !SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        double precision scn(nsdim) 	!non-cohesive suspended sediment conc (kg/m3)
        double precision usb(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF BEDLOAD
        double precision uss(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF SUSPENDED
        double precision ust(nsdim)     !CRITICAL SHEAR VELOCITY FOR INITIATION OF SHEET FLOW
        double precision dxx(nsdim)     !DIMENSIONLESS GRAIN SIZE


! ------------ OUTPUT VARIABLES -----------------
        double precision sedx(nsdim)	!x bedload component [kg/ms]
        double precision sedy(nsdim)	!y bedload component [kg/ms]
        double precision sflux(nsdim)   !flux of suspend sediment [m3/m2]
	real sload(nsdim)		!suspended sediment load [kg/m2s]

! ------------ LOCAL VARIABLES -----------------
        DOUBLE PRECISION QS       !SUSPENDED SEDIMENT TRANSPORT RATE (KG/M/S)
        DOUBLE PRECISION QSDIR    !DIRECTION OF SUSPENDED SEDIMENT TRANSPORT (DEGREE)
        DOUBLE PRECISION SED      !TIME-AVERAGED NET SEDIMENT TRANSPORT AS VOLUME (M**3/S/M)
        DOUBLE PRECISION SEDM     !TIME-AVERAGED NET SEDIMENT TRANSPORT AS MASS (KG/S/M)
        DOUBLE PRECISION SEDDIR   !DIRECTION OF NET SEDIMENT TRANSPORT (AZIMUTH,DEGREES)
        DOUBLE PRECISION C0       !REFERENCE CONCENTRATION AT Z0 (KG/M^3)
        DOUBLE PRECISION C0A  	  !DEPTH AVERAGED REFERENCE CONCENTRATION AT Z0 (KG/M^3)
        DOUBLE PRECISION C0AI  	  !DEPTH AVERAGED REFERENCE CONCENTRATION AT Z0 (KG/M^3) per is
        integer is		  !counters
        real smax

!       --------------------------------------------------
!       Loop over grainsize classes for bedload (non-cohesive only)
!       --------------------------------------------------
 
        do is = nbcc+1,nsclsc
         sedx(is) = 0. 
         sedy(is) = 0. 

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
         SEDM = SEDM / bedcha(1,3)

         sedx(is) = SEDM * pers(is) * cos(SEDDIR)
         sedy(is) = SEDM * pers(is) * sin(SEDDIR)

	end do

	if (IOPT.eq.1 .or. IOPT.eq.3 ) goto 130

!       --------------------------------------------------
!       Loop over grainsize classes for suspension
!       --------------------------------------------------
 
        do is = 1,nsclsc

	 C0A = 0.d0

!        -------------------------------------------------------------------
!        Calculate velocity profile, suspended sediment concentration profile
!        -------------------------------------------------------------------

         CALL PROFL(DL,RHOW,ws(is),UB,CDIR,USTCWS,USTCW,usb(is),
     $uss(is),ust(is),Z0C,DELTACW,USTCS,USTWS,USTCWSB,USTC,USTW,
     $RPLCOEF,USTBF,Z0,gs(is),dxx(is),RHOS,QS,QSDIR,C0,C0A)

         C0AI = C0A * pers(is)

!        -------------------------------------------------------------------
!        Compute load for suspended concentration [kg/m²s]
!        -------------------------------------------------------------------

         sload(is) = ws(is)*(C0AI - scn(is))

!        -------------------------------------------------------------------
!        Limit load and set sediment bulk density
!        -------------------------------------------------------------------

         if (scn(is) .gt. C0AI) then		!deposition
           smax = -scn(is)*DL/TIMEDR
           sload(is) = max(sload(is),smax)
            if (is .le. nbcc) then
              rhos = RHOCLAY * (1.d0 - SURFPOR)
            else
              rhos = RHOSED * (1.d0 - SURFPOR)
            end if
         else					!erosion
           smax = pers(is)*BEDCHA(1,3)*BEDCHA(2,1)/TIMEDR
           sload(is) = min(sload(is),smax)
           rhos = BEDCHA(1,3)
         end if

!        -------------------------------------------------------------------
!        Compute deposition or erosion [m³/m²]
!        -------------------------------------------------------------------

	 sflux(is) =  -sload(is) * timedr / rhos
        end do

130     continue

        end

! ********************************************************************
! SUBROUTINE BEDLOAD
! This subroutine computes the bedload sediment transport using the 
! sediment continuity equation

        subroutine bedload(nscls,sedx,sedy,dt,bh,bflx)

	use mod_bound_geom
	use evgeom
	use basin

        implicit none

        include 'param.h'
        include 'sed_param.h'


! --- input variables
        integer nscls				!number grainsize classes
        double precision sedx(nsdim,nkndim)	!bedload transport in x direction [kg/ms]
        double precision sedy(nsdim,nkndim)	!bedload transport in y direction [kg/ms]
	real dt					!time step
        real bh(nkndim)			    	!bottom height variation [m]

! --- output variables
        double precision bflx(nsdim,nkndim)	!flux of bedload sediment [m3/m2]

! --- local variables
        integer k,ie,ii,is
	double precision v2v(nkndim)
        double precision bflux
        double precision sexe(nsdim)           	!element transport in x direction
        double precision seye(nsdim)           	!element transport in x direction
        double precision b,c   		        !x and y derivated form function [1/m]
        double precision area			!area of element/3
        double precision ss(nsdim)
	double precision sm(nsdim)

!       -------------------------------------------------------------------
!       Initialize local variables
!       -------------------------------------------------------------------

        do k = 1,nkn
         do is = 1,nscls
           bflx(is,k) = 0.d0
         end do
	 v2v(k) = 0.d0
        end do

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
          end do
        end do

!       -------------------------------------------------------------------
!       Apply slope limiter
!       -------------------------------------------------------------------

        !call slope_lim_bh(nscls,bflx,bh)

        end

! ********************************************************************

	subroutine slope_lim_bh(nscls,bflx,bh)

	use basin, only : nkn,nel,ngr,mbw

	implicit none

        include 'param.h'
        include 'sed_param.h'

        integer nscls				!number of grainsize class
        double precision bflx(nsdim,nkndim)	!bedload sediment contribution [m3/m2]
        real bh(nkndim)				!bed elevation

	integer k,is
	double precision bbhh(nkndim)
	double precision baux(nkndim)
	double precision frac(nsdim,nkndim)
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
	use basin

        implicit none

        include 'param.h'


        double precision v1v(nkndim),v2v(nkndim),v3v(nkndim),v4v(nkndim)
	double precision bbk(nkndim),bbe(neldim)

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

	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'
        include 'sed_param.h'

        integer nscls				!number of grainsize class
        double precision gs(nsdim)		!grainsize class
        double precision timedr                	!time step [s]
        double precision bflux(nsdim,nkndim)	!bedload sediment contribution [m3/m2]
        double precision sflux(nsdim,nkndim)	!suspended sediment contribution [m3/m2]
        real bh(nkndim)				!bed elevation
        real gskm(nkndim)			!average sediment grainsize in node k

        double precision percbd(nlbdim,nsdim,nkndim)        !fraction of sediment [0,1]
        double precision bedn(nlbdim,3,nkndim)  !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)

	double precision bdh(nkndim)		!total elevation change [>0depo,<0ero]
        double precision dzco			!elevation variation due to compaction [m]
        double precision flx 			!bedload + suspended flux [m]
        double precision timecomp		!time step for compaction
        double precision gsa,gss
        double precision gsmax
        double precision gsmin
        integer k,is
	logical is_d_nan

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
	  if (is_d_nan(flx)) flx = 0.
	  flx = min(flx,0.01)
	  flx = max(flx,-0.01)

!         -------------------------------------------------------------------
!         Multiply by morphological accelaration factor
!         -------------------------------------------------------------------

          !flx = flx * morpho

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
!        Compute total bottom thikness variation
!        -------------------------------------------------------------------

	 if (is_d_nan(bdh(k))) bdh(k) = 0.
	 bdh(k) = min(bdh(k),0.01)
	 bdh(k) = max(bdh(k),-0.01)
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

        subroutine bedact(nscls,bmix,BEDCHA,percbd)

        implicit none

        include 'sed_param.h'

        double precision percbd(nlbdim,nsdim)	!fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)         !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nlbd				!number of bed layer in k
        integer nscls                           !number of grainsize class
        double precision bmix                   !thickness of active layer [m]
        double precision ptot
        integer ib,is
	logical checkbedcha
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
            percbd(2,is) = (percbd(1,is)*(BEDCHA(2,1)-bmix) + 
     $percbd(2,is)*(BEDCHA(3,1)-BEDCHA(2,1)))/(BEDCHA(3,1)-bmix)
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

          do is = 1,nscls			!update percbd in layer 1
            percbd(1,is) = (percbd(2,is)*(BEDCHA(3,1)-BEDCHA(2,1)) +
     $percbd(1,is)*BEDCHA(2,1))/BEDCHA(3,1)
          end do

          do ib = 2,nlbd			!shift one layer up
           do is = 1,nscls
            percbd(ib,is) = percbd(ib+1,is)
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
          percbd(1,is) = (percbd(1,is)*BEDCHA(2,1) + percbd(2,is)*
     $(bmix-BEDCHA(2,1)))/bmix
        end do

        BEDCHA(2,3) = lin(BEDCHA(2,1),BEDCHA(3,1),BEDCHA(2,3),
     $BEDCHA(3,3),bmix)
        BEDCHA(2,2) = lin(BEDCHA(2,1),BEDCHA(3,1),BEDCHA(2,2),
     $BEDCHA(3,2),bmix)
        BEDCHA(2,1) = bmix

	end if

!       -------------------------------------------------------------------
!       Adjusts percbd and percentage check ---> ptot must be = 1 +-0.1
!       -------------------------------------------------------------------

        ptot = 0.d0
        do is = 1,nscls
          if(percbd(1,is).lt.1D-5) percbd(1,is) = 0.d0
	  percbd(1,is) = dmin1(percbd(1,is),1.d0)
         ptot = ptot + percbd(1,is)
        end do
        if(ptot.lt.0.99d0.or.ptot.gt.1.01d0) go to 130

        do is = 1,nscls
         percbd(1,is) = percbd(1,is)/ptot
        end do

!       -------------------------------------------------------------------
!       Create last layer to the minimum layer number
!       -------------------------------------------------------------------

        if ( nlbd .le. 5 ) then
          do is = 1,nscls
            percbd(nlbd+1,is) = percbd(nlbd,is)
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
          write(6,*) 'classnumber:',is,'percentage:',percbd(1,is)
        enddo
        stop 'error stop : percbd'

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

        subroutine unilayer(nscls,BEDCHA,percbd)

        implicit none

        include 'sed_param.h'

        double precision percbd(nlbdim,nsdim)	!fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)         !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nlbd				!number of bed layer in k
        integer nscls                           !number of grainsize class
        double precision ptot
        integer ib,is
	logical checkbedcha
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

          do is = 1,nscls			!update percbd in layer 1
            percbd(2,is) = (percbd(2,is)*(BEDCHA(3,1)-BEDCHA(2,1)) +
     $percbd(3,is)*(BEDCHA(4,1)-BEDCHA(3,1)))/(BEDCHA(4,1)-BEDCHA(2,1))
          end do

          do ib = 3,nlbd			!shift one layer up
           do is = 1,nscls
            percbd(ib,is) = percbd(ib+1,is)
           end do
           BEDCHA(ib,1) = BEDCHA(ib+1,1)
           BEDCHA(ib,2) = BEDCHA(ib+1,2)
           BEDCHA(ib,3) = BEDCHA(ib+1,3)
          end do

        end if

!       -------------------------------------------------------------------
!       Adjusts percbd and percentage check of layer 2 ---> ptot must be = 1 +-0.1
!       -------------------------------------------------------------------

        ptot = 0.d0
        do is = 1,nscls
          if(percbd(1,is).lt.1D-5) percbd(1,is) = 0.d0
	  percbd(1,is) = dmin1(percbd(1,is),1.d0)
         ptot = ptot + percbd(2,is)
        end do
        if(ptot.lt.0.99d0.or.ptot.gt.1.01d0) go to 130

        do is = 1,nscls
         percbd(2,is) = percbd(2,is)/ptot
        end do

!       -------------------------------------------------------------------
!       Create last layer to the minimum layer number
!       -------------------------------------------------------------------

        if ( nlbd .le. 5 ) then
          do is = 1,nscls
            percbd(nlbd+1,is) = percbd(nlbd,is)
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
          write(6,*) 'classnumber:',is,'percentage:',percbd(2,is)
        enddo
        stop 'error stop : percbd'

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

        subroutine checkbed(k,iss,nscls,gs,BEDCHA,percbd,flux)

        implicit none

        include 'sed_param.h'
        
        double precision percbd(nlbdim,nsdim)	!fraction of sediment [0,1]
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
	logical checkbedcha

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
        bpre = BEDCHA(2,1)*percbd(1,iss)
        bnet = bpre  + flux

!       -------------------------------------------------------------------
!       Compute erosion or deposition
!       -------------------------------------------------------------------

        if (flux.gt.0.d0) then

          call depbed(iss,nlbd,nscls,BEDCHA,percbd,flux,bnet,btot,gs)

        elseif(flux.lt.0.d0) then

          call erobed(iss,nlbd,nscls,BEDCHA,percbd,flux,bnet,btot,bpre)

        end if
  
!       -------------------------------------------------------------------
!       Adjust percbd and percentage check ---> ptot must be = 1 +-0.1
!       -------------------------------------------------------------------

        if (nscls.eq.1) percbd(1,nscls) = 1.d0

        ptot = 0.d0
        do is = 1,nscls
          if(percbd(1,is).lt.1D-5) percbd(1,is) = 0.d0
	  percbd(1,is) = dmin1(percbd(1,is),1.d0)
          ptot = ptot + percbd(1,is)
        end do
        if(ptot.lt.0.99d0.or.ptot.gt.1.01d0) go to 130

        do is = 1,nscls
         percbd(1,is) = percbd(1,is)/ptot
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
          write(6,*) 'classnumber:',is,'percentage:',percbd(1,is)
        enddo
        stop 'error stop : percbd'

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

        subroutine depbed(iss,nlbd,nscls,BEDCHA,percbd,flux,bnet,btot,
     $			  gs)

        implicit none

        include 'sed_param.h'

        double precision percbd(nlbdim,nsdim)   !fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)       !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nlbd                            !number of bed layer in k
        integer nscls                           !number of grainsize class
        double precision flux                   !flux of sediment [m] (<0 ero, >0 depo)
        double precision bnet                   !new net sediment in the layer [m]
        double precision btot                   !total sediment [m]
        double precision gs                     !grainsize of first class
        double precision rhosa                  !non consolidated density [kg/m3]
	double precision taunew			!new critical shear stress
	double precision tausurf		!average surface critical shear stress
	double precision rhosurf		!average surface density
        integer ib,is,iss

!       -------------------------------------------------------------------
!       Initialize valiables
!       -------------------------------------------------------------------

        rhosa = RHOSED * (1.d0 - SURFPOR)
        if (gs.lt.limcoh) rhosa = RHOMUD
        taunew = teroa*(rhosa**terob)
        !tausurf = (BEDCHA(1,2) + BEDCHA(2,2))/2.
        tausurf = BEDCHA(1,2)
        !rhosurf = (BEDCHA(1,3) + BEDCHA(2,3))/2.
        rhosurf = BEDCHA(1,3)

!       -------------------------------------------------------------------
!       Normal deposition
!       -------------------------------------------------------------------

        do is = 1,nscls
          percbd(1,is) = percbd(1,is)*BEDCHA(2,1)/btot
        end do

        percbd(1,iss) = bnet / btot

        BEDCHA(1,2) = (taunew*flux + tausurf*BEDCHA(2,1)) / btot
        BEDCHA(1,2) = dmax1(BEDCHA(1,2),taunew)
        BEDCHA(1,2) = dmin1(BEDCHA(1,2),BEDCHA(2,2))

        BEDCHA(1,3) = (rhosa*flux + rhosurf*BEDCHA(2,1))/ btot
        BEDCHA(1,3) = dmin1(BEDCHA(1,3),BEDCHA(2,3))

        do ib = 2,nlbd
          BEDCHA(ib,1) = BEDCHA(ib,1) + flux
        end do

	end 

! ********************************************************************
! Excess deposition, create new layer and upgrade level 2

        subroutine newlayer(nscls,BEDCHA,percbd)

        implicit none

        include 'sed_param.h'

        double precision percbd(nlbdim,nsdim)   !fraction of sediment [0,1]
        double precision BEDCHA(nlbdim,3)         !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)
        integer nscls                           !number of grainsize class
        double precision bres                   !sediment that is left [m]
        double precision blimit                 !maximum thikness of layer [m]
        integer nlbd                            !number of bed layer in k
        integer ib,is
	logical checkbedcha
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

        blimit = BEDCHA(nlbd,1) - BEDCHA(nlbd-1,1)
	bres = bedcha(2,1) - blimit

	if ( bres .gt. (blimit/2.) ) then

!         -------------------------------------------------------------------
!         Create a new to layer and shift the others down
!         -------------------------------------------------------------------

          nlbd = nlbd + 1	
          do ib = nlbd,3,-1
           do is = 1,nscls
            percbd(ib,is) = percbd(ib-1,is)
           end do
           BEDCHA(ib,1) = BEDCHA(ib-1,1)
           BEDCHA(ib,2) = BEDCHA(ib-1,2)
           BEDCHA(ib,3) = BEDCHA(ib-1,3)
          end do

!         -------------------------------------------------------------------
!         Update layer 2 properties
!         -------------------------------------------------------------------

          do is = 1,nscls
            percbd(2,is) = percbd(1,is)
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
       
        subroutine erobed(iss,nlbd,nscls,BEDCHA,percbd,flux,bnet,btot,
     $bpre)

        implicit none

        include 'sed_param.h'

	double precision d0
	parameter (d0 = 0.d0)

        double precision percbd(nlbdim,nsdim)   !fraction of sediment [0,1]
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
	  if (btot. eq. d0 .or. percbd(1,iss) .gt. 0.999d0) then
	    call dellayer(nlbd,nscls,BEDCHA,percbd)
	    return
	  endif
        endif

!       -------------------------------------------------------------------
!       Normal erosion
!       -------------------------------------------------------------------

        do is = 1,nscls
          percbd(1,is) = percbd(1,is)*BEDCHA(2,1) / btot
        end do

        percbd(1,iss) = bnet / btot

        BEDCHA(1,2) = lin(d0,BEDCHA(2,1),BEDCHA(1,2),BEDCHA(2,2),-flux)
        BEDCHA(1,3) = lin(d0,BEDCHA(2,1),BEDCHA(1,3),BEDCHA(2,3),-flux)

        do ib = 2,nlbd
          BEDCHA(ib,1) = BEDCHA(ib,1) + flux
        end do

        end

! ********************************************************************
! SUBROUTINE DELLAYER
! Delete layers and update bed characteristics

        subroutine dellayer(nlbd,nscls,BEDCHA,percbd)

        implicit none

        include 'sed_param.h'

        double precision percbd(nlbdim,nsdim)   !fraction of sediment [0,1]
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
            percbd(nlbd+1,is) = percbd(nlbd,is)
          end do
          BEDCHA(nlbd+1,1) = BEDCHA(nlbd,1) + (BEDCHA(nlbd,1) -
     $BEDCHA(nlbd-1,1))
          BEDCHA(nlbd+1,2) = BEDCHA(nlbd,2)
          BEDCHA(nlbd+1,3) = BEDCHA(nlbd,3)
        end if

        do ib = 1,nlbd                       !shift one layer up
         do is = 1,nscls
          percbd(ib,is) = percbd(ib+1,is)
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

        subroutine totcon(nscls,scn,tcon)

	use mod_ts
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'
        include 'sed_param.h'


	real rhos,conc
        integer nscls				!number of grainsize class
        real scn(nlvdim,nkndim,nsdim)  		!suspended concentration
        real tcon(nlvdim,nkndim)		!total concentration

        integer k,l,is,lmax

        do k = 1,nkn
         lmax = ilhkv(k)
         do l = 1,lmax
	  tcon(l,k) = 0.
          do is = 1,nscls
           scn(l,k,is) = max(scn(l,k,is),0.)
	   conc = scn(l,k,is)

	   if (is .le. nbcc) then
	     rhos = rhoclay
	   else
	     rhos = rhosed
	   end if

!          -------------------------------------------------------------------
!          Compute total concentration
!          -------------------------------------------------------------------

           tcon(l,k) = tcon(l,k) + conc

!          -------------------------------------------------------------------
!          Update fluid density
!          -------------------------------------------------------------------

           rhov(l,k) = rhov(l,k) + (conc/rhos)*(rhos-rhov(l,k))

          end do
         end do
        end do

        call nantest(nkn*nlvdim,rhov,'rhov')
            
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
	use basin

        implicit none

        include 'param.h'

        double precision bdh(nkn)          !total elevation change [>0depo,<0ero]

        real hlhv(nel)
	real v1v(nkn)

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
          if ((hlhv(ie) - evdep) .lt. 0.2 ) evdep = 0.

          hlhv(ie) = hlhv(ie) - evdep
	  do ii = 1,3
            hm3v(ii,ie) = hm3v(ii,ie) - evdep
          end do
	  hev(ie) = hev(ie) - evdep

        end do

!       ------------------------------------------------------------------
!       Set up depth vectors
!       ------------------------------------------------------------------

        call makehkv(hkv,v1v)         !computes hkv as average

	call setweg(3,iw)
        call setarea(nlvdim,areakv)
        call setdepth(nlvdim,hdknv,hdenv,zenv,areakv)

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
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'

        real kvalue(nkndim)             !variable
        double precision kdiff(nkndim)	!istantaneous variable difference
        double precision smooth		!smoothing factor for morphodynamic [0-1]
        real gdx(nkndim),gdy(nkndim)	!slope gradients
        double precision angrep		!angle of repose [rad]
        real ksl                        !node slope angle [radian]
        real areanode,area
        real baux(nkndim)
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
! Reset bottom height variation after initialization time adjtime
! The bottom height variation is resetted, while all other variables
! (sediment density, bathymetry, tauce) are not.

        subroutine resetsedi(it,adjtime,bh,scn)

	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'
        include 'sed_param.h'


	integer it                      !time in seconds
	integer adjtime			!time for initialization [s]
	real bh(nkndim)                 !bottom height variation [m]
        real scn(nlvdim,nkndim,nsdim)   !suspended sediment conc (kg/m3)

        integer k,l,is

	if (it .eq. adjtime) then 
	  write(6,*)'Morphological initialization time:'
	  write(6,*)'reset bathymetric change but keep modified'
	  write(6,*)'sediment characteristics and bathymetry'
          do k = 1,nkn
	    bh(k) = 0.
	    do l = 1,nlvdim
	      do is = 1,nsdim
		scn(l,k,is) = 0.
	      end do
	    end do
	  end do
	end if

	end 

 !******************************************************************

        subroutine sedcons(it,tcon,bh,bedn)

	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'
        include 'sed_param.h'


        double precision bedn(nlbdim,3,nkndim)  !bed characteristics in 3 column table
                                                ! (1) depth below sediment surface (m)
                                                ! (2) critical erosion stress (Pa)
                                                ! (3) dry bulk density (kg/m**3)

	integer it
        real tcon(nlvdim,nkndim)		!total concentration
        real bh(nkndim)				!bed elevation
	real net,ccc,bbb
        real areanode,area
        real volnode,vol

        integer k

	net = 0.
	bbb = 0.
	ccc = 0.
	do k = 1,nkn
	  area = areanode(1,k)
	  vol = volnode(1,k,1)
	  ccc = ccc + tcon(1,k) * vol
	  bbb = bbb + bh(k)*area*real(bedn(1,3,k))
	  net = ccc + bbb 
	end do

	write(74,*)it,bbb,ccc,net

	end

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
	use basin

	implicit none

! parameters
	include 'param.h'
! arguments
	real vv(nkndim)
	real vv1(nkndim)
! common



! local
	real ukbot(nkndim),vkbot(nkndim)
	real u,v
	real ddl(nkndim)
	integer ie,lmax,k,ii
	real aj,vol,depth
	integer nsigma
        real hsigma
	logical bsigma

	if(nlvdim.ne.nlvdi) stop 'error stop : level dim in uvbott'

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

        subroutine set_settling(k,lmax,nscls,ws,scn,wsink)

        implicit none

        include 'param.h'
        include 'sed_param.h'

	integer k,lmax
        integer nscls				!number of grainsize class
        double precision ws(nsdim)      	!SETTLING VELOCITY FOR NON-COHESIVE SEDIMENT (M/SEC)
        real scn(nlvdim,nkndim,nsdim)  		!suspended concentration
        real wsink(0:nlvdim,nkndim,nsdim)       !settling velocity for suspended sediment

	double precision conc(nsdim)
        real salt,temp			!salinity [psu] and temperature [C]
        DOUBLE PRECISION TEM		!in degrees Celsius
        DOUBLE PRECISION SAL		!practical salinity scale (seawater=35)
        DOUBLE PRECISION VISC   	!dynamic viscosity of the (sea)water (KG/M*SEC (OR N.S/M**2))
        DOUBLE PRECISION VISK     	!KINEMAMIC VISCOSITY OF THE FLUID (KG/M*SEC (OR N.S/M**2))
        DOUBLE PRECISION RHOW     	!DENSITY OF FLUID (WATER)  (KG/M**3)
        DOUBLE PRECISION WWS(nsdim)     !settling velocity of the considered floc class (m/s)

        integer l,is,i

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
             wsink(l,k,i) = real(WWS(i))
           end do
	 end do
	end if

	do l = 1,lmax
	  do is = nbcc+1,nscls
            wsink(l,k,is) = real(ws(is))
	  end do
	end do

	do is = 1,nscls
	  wsink(0,k,is) = wsink(1,k,is)
	end do

	end

!******************************************************************

	subroutine load3d(sload,nkn,nlvdi,ilhkv,load)

	implicit none

	real sload(1)
	integer nkn
	integer nlvdi
        integer ilhkv(1)             !number of element and node level

	real load(nlvdi,1)

	integer k,l,lmax

	do k = 1,nkn
	   lmax = ilhkv(k)
	   do l = 1,lmax-1
	     load(l,k) = 0. 
	   end do
	   load(lmax,k) = sload(k)
	end do

	end

!******************************************************************
