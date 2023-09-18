!
! $Id: subsys.f,v 1.98 2010-03-11 15:36:39 georg Exp $
!
! routines system dependent
!
! contents :
!
! subroutine nlsinh             initializes the hp parameter file
! subroutine nlsina             initializes the ap parameter file
! subroutine fnminh             initializes default names
!
! revision log :
!
! 22.01.1998	ggu	no more strdir, apndir -> read apnpar.str from local
! 18.03.1998	ggu	no '[not defined...]' in basdir...
! 20.03.1998	ggu	for dos systems .memory -> _memory
! 22.03.1998	ggu	some changes for Technital integrated
! 30.04.1998	ggu	some changes for flux sections
! 13.05.1998	ggu	added colfil (apncol.str)
! 25.05.1998	ggu	documentation added
! 26.05.1998	ggu	documentation for wind file added
! 19.06.1998	ggu	some useless comments deleted
! 24.06.1998	ggu	qflux (heat flux) added
! 22.07.1998	ggu	more on documentation
! 23.07.1998	ggu	documentation
! 12.08.1998    ggu     new parameter dlat -> specify latitude for coriolis
! 02.09.1998    ggu     hack: change depth in Venice inlets (hlido,...)
! 24.11.1998    ggu     switch for biological reactor introduced
! 12.02.1999    ggu     default of parameter file changed to apnstd.str
! 12.02.1999    ggu     new parameters for plotting
! 13.04.1999    ggu     new parameter itide
! 19.04.1999    ggu     itide changed to rtide
! 27.05.1999    ggu     icust introduced
! 31.05.1999    ggu     dval default changed
! 28.09.1999    ggu     new flag regflx
! 19.11.1999    ggu     new parameters for section vol
! 08.08.2000    ggu     hlvmin is now percentage of last layer thickness
! 03.12.2001    ggu     new parameters (rlhdif)
! 11.10.2002    ggu     rstot new meaning
! 14.10.2002    ggu     atpar,adpar,aapar -> default set to 1.0 (implicit)
! 05.10.2003    ggu     changes in color handling of post routines
! 04.12.2003    ggu     sediment and wave module integrated
! 05.03.2004    ggu     bio variable for initialization
! 03.09.2004    ggu     restart variables
! 22.09.2004    ggu     nlsina_3d(): always use 3d file for plot (level=0)
! 28.09.2004    ggu     lagrangian routines integrated (LAGR)
! 05.10.2004    ggu     some more documentation
! 02.12.2004    ggu     documentation for variable time step
! 06.12.2004    ggu     new subroutine nlsina_legvar
! 17.01.2005    ggu     new parameters for horizontal diffusion
! 15.03.2005    ggu     cleaned horiz. diff., read new section legvar
! 19.05.2005    ggu     added time for custom routine (tcust)
! 04.11.2005    ggu     new parameter itlin (semi-lagrangian), some corrections
! 07.11.2005    ggu     new parameter itvd (TVD)
! 07.11.2005    ggu     parameters deleted: isedi,sedref,sedgrs
! 08.11.2005    ggu     do not set old parameters, some in nlsina_para
! 08.11.2005    ggu     more documentation
! 16.02.2006    ggu     new flag itoxi and file toxi
! 23.03.2006    ggu     new variable itunit for unit of time step
! 12.09.2006    ggu     time and date introduced
! 28.09.2006    ggu     new iprogr for progress of simulation
! 18.10.2006    ccf     new params iwave and dtwave for wave model
! 15.11.2006    ggu     new parameters to construct streched vert. coordinates
! 16.11.2006    ggu     use itemp,isalt to decide about advection
! 29.11.2006    ggu     rwhpar for horizontal diffusion in lagrangian model
! 27.08.2007    ccf     isphe = 1  for spherical coordinate system
! 07.04.2008    aac     parameters for ersem ecological model call
! 10.04.2008    ccf     netcdf and gotmpa parameters
! 17.04.2008    ggu     new param ievap (to compute evaporation mass flux)
! 28.04.2008    ggu     rstot deleted, contau introduced
! 29.04.2008    ggu&aac new vars for ERSEM
! 16.06.2008    ggu     old parts deleted, new documentation
! 17.09.2008    ggu     new interpretation for level (-1 = bottom)
! 11.10.2008	ggu	zfranco added
! 10.11.2008	ggu	new variable ilytyp, cleaned
! 19.11.2008	ggu	new variable noslip
! 24.11.2008	ggu	new variable vreps (mass error)
! 06.12.2008	ggu	new variables nbfix and nsigma
! 09.01.2009	ggu	documentation
! 21.01.2009	ggu	new variable vrerr (stop if mass error)
! 29.01.2009	ggu	various changes, better documentation
! 13.03.2009	ggu	bugfix: some parameters had default section not set
! 07.05.2009	ggu	new parameter ityrst
! 15.06.2009	ggu	new parameters for plot: isphe, reggrd, reggry
! 11.09.2009	ggu	new section $sect
! 09.10.2009	ggu	new parameter sclvel
! 13.10.2009	ggu	documentation is $sect, new hvmax, lvmax
! 22.02.2010	ggu	new parameter itdrag
! 23.02.2010	ggu	new parameter regdst
! 26.03.2010	ggu	new parameters for arrows in section plot
! 28.09.2010	ggu	new value for icor
! 29.09.2010	ggu	new param vmode,rxscal,ryscal
! 15.12.2010	ggu	nsigma renamed to nbsig, nsigma used for sigma layers
! 21.12.2010	ggu	new parameter rwscal
! 16.02.2011	ggu	new default for isphe, new routine fnm_aquabc_init()
! 25.02.2011	ggu	new param wsmax to catch errors in wind type
! 23.03.2011	ggu	new parameter itvdv
! 24.03.2011	ggu	new parameters iheat,hdecay,botabs
! 01.06.2011	ggu	new parameter idtmin
! 18.08.2011	ggu	new parameter isoinp (interpolate inside element)
! 18.09.2011	ggu	change default for isphe for output (-1)
! 03.11.2011	ggu	new parameter hsigma (hybrid)
! 18.11.2011	ggu	new subroutine nlsinh_proj() for projection
! 24.01.2012	ggu	new parameter nomp
! 02.05.2012	ggu	new default for ndccol (-> 0)
! 24.10.2012	ggu	new parameter dxmin
! 10.05.2013	ggu	new parameters idtbox,itmbox, more comments
! 10.05.2013	ggu	new parameter inohyd
! 16.05.2013	ggu	file name bound renamed to zinit
! 28.03.2014	ggu	some new params for lagrangian
! 10.04.2014	ccf	new section "wrt" for water renewal time
! 30.05.2014	ggu	new default for dragco, new metpnt
! 20.10.2014	ggu	new default for date (-1)
! 01.12.2014	ccf	wave parameters moved to subwave.f
! 06.05.2015	ccf	new parameters itmoff and offlin
! 29.09.2015	ggu	new boundary file surfvel
! 29.09.2015	ggu	new initial file uvinit, new flgrst
! 01.02.2016	ggu	some plot params shifted to para section (bbgray, etc.)
! 22.02.2016	ggu	new file for initial conditions bfmini
! 05.04.2016	ggu	new parameter iaicef for ice free areas
!
!************************************************************************
!--------------------------------------------------------------------------------
        module def_para
!--------------------------------------------------------------------------------

        use para

!--------------------------------------------------------------------------------
        contains
!--------------------------------------------------------------------------------

	subroutine nlsinh

! initializes the parameter file for the main FE model
        
	implicit none

	call nlsinh_general
	call nlsinh_lagrg
        call nlsinh_wrt
	call nlsinh_bfmsc
	call nlsinh_proj
	call nlsinh_undoc
	call nlsinh_georg
	call nlsinh_unused
	call nlsinh_waves
	call nlsinh_nonhydro

	end

!************************************************************************

	subroutine nlsinh_general

	implicit none

	call sctpar('para')		!sets default section

! DOCS	START	S_para_h
!
! DOCS	COMPULS		Compulsory time parameters
!
! This parameters are compulsory parameters that define the
! period of the simulation. They must be present in all cases.
!
! |itanf|	Start of simulation. (Default 0)
! |itend|	End of simulation.
! |idt|		Time step of integration.

	call addpar('itanf',0.d0)
	call addpar('itend',0.d0)
	call addpar('idt',0.d0)

!c------------------------------------------------------------------------

! DOCS	OUTPUT		Output parameters
!
! The following parameters deal with the output frequency
! and start time to external files. The content of the various
! output files should be looked up in the appropriate section.
!
! The default for the time step of output of the files is 0 which
! means that no output file is written. If the time step of the
! output files is equal to the time step of the simulation then
! at every time step the output file is written. The default start time
! of the output is |itanf|, the start of the simulation.
!
! |idtout|, |itmout|	Time step and start time for writing to file OUT,
!			the file containing the general hydrodynamic results.

	call addpar('idtout',0.d0)
	call addpar('itmout',-1.d0)

! |idtext|, |itmext|	Time step and start time for writing to file EXT,
!			the file containing hydrodynamic data of extra points.
!			The extra points for which the data is written
!			to this file are given in section |extra| of
!			the parameter file.

	call addpar('idtext',0.d0)
	call addpar('itmext',-1.d0)

! |idtrst|		Time step for writing the restart
!			file (extension RST). No restart file is written
!			with |idtrst| equal to 0. A negative value
!			is also possible for the time step. In this case
!			the time step used is |-idtrst|, but the file is
!			overwritten every time. It therefore contains 
!			always only the last written restart record. The
!			special value of |idtrst = -1| will write only the
!			last time step of the simulation in the restart file. 
!			This is useful if you want to start another
!			simulation from the last output. (Default 0)
! |itmrst|		Start time for writing the restart file. If
!			not given it is the beginning of the simulation.
! |itrst|		Time to use for the restart. If a restart
!			is performed, then the file name containing
!			the restart data has to be specified in |restrt|
!			and the time record corresponding to |itrst|
!			is used in this file. A value of -1 is also possible.
!			In this case the last record in the restart file
!			is used for the restart and the simulation starts
!			from this time. Be aware that this option changes
!			the parameter |itanf| to the time of the last
!			record found in |restrt|.
! |ityrst|		Type of restart. If 0 and the restart file is not
!			found the program will exit with an error. Otherwise
!			the program will simply continue with a cold start.
!			If |ityrst| is 1 and the given time record is not
!			found in the file it will exit with error. If
!			it is 2 it will initialize all values from the
!			first time record after |itrst|. Therefore, the
!			value of 2 will guarantee that the program will not
!			abort and continue running, but it might not
!			be doing what you intended. (Default 0)
! |flgrst|		This variable indicateswhich variables are read
!			from the restart file. By default all available
!			variables are read and used. If some variables
!			are not wanted (because, e.g., you want to restart
!			from a different T/S field), this fact can be indicated
!			in |flgrst|. 1 indicates restart of hydro values,
!			10 the depth values, 100 T/S values, 1000 the tracer
!			concentration, 10000 vertical velocities and
!			100000 the ecological variables. Therefore, a value
!			of 10111 indicates a restart of everything except
!			the tracer and the ecological values. The default
!			value for |flgrst| is -1, which means 111111.

	call addpar('idtrst',0.d0)
	call addpar('itmrst',-1.d0)
	call addpar('itrst',0.d0)
	call addpar('ityrst',0.d0)
	call addpar('flgrst',-1.d0)

! |idtres|, |itmres|	Time step and start time for writing to file RES,
!			the file containing residual hydrodynamic data.

	call addpar('idtres',0.d0)
	call addpar('itmres',-1.d0)

! |idtrms|, |itmrms|	Time step and start time for writing to file RMS,
!			the file containing hydrodynamic data of root mean
!			square velocities.

	call addpar('idtrms',0.d0)
	call addpar('itmrms',-1.d0)

! |idtflx|, |itmflx|	Time step and start time for writing to file FLX,
!			the file containing discharge data through defined
!			sections.
!			The transects for which the discharges are computed
!			are given in section |flux| of
!			the parameter file.

	call addpar('idtflx',0.d0)
	call addpar('itmflx',-1.d0)

! |idtvol|, |itmvol|	Time step and start time for writing to file VOL,
!			the file containing volume information of areas
!			defined by transects.
!			The transects that are used to compute the volumes
!			are given in section |volume| of
!			the parameter file.

	call addpar('idtvol',0.d0)
	call addpar('itmvol',-1.d0)

! |netcdf|		This parameter chooses output in NetCDF format 
!			if |netcdf| is 1, else the format is unformatted
!			fortran files. (Default 0)

	call addpar('netcdf',0.d0)

! |idtoff|	handles offline mode (default 0):
!		\begin{description}
!		\item[0] do nothing (no offline routines called)
!		\item[$>$0] write offline data file (.off) with time step |idtoff|.
!		\item[$<$0] reads offline data from file |offlin|
!		defined in section |name|. Usage:
!		\begin{description}
!		\item[-1] uses offline hydro results
!		\item[-2] uses offline T/S results
!		\item[-4] uses offline turbulence results
!		\end{description}
!		\end{description}

	call addpar('idtoff',0.d0)

! |itmoff|	Start time for writing to file OFF,
!		the file containing data for offline runs.
	call addpar('itmoff',-1.d0)

!c------------------------------------------------------------------------

! DOCS	TIME_DATE	General time and date parameters
!
! A time and date can be assigned to the simulation. These values
! refer to the time 0 of the FEM model. The format for the date is
! YYYYMMDD and for the time HHMMSS.
! You can also give a time zone if your time is not referring to
! GMT but to another time zone such as MET.

! |date|                The double precision date corresponding to time 0. (Default 0)
! |time|                The double precision time corresponding to time 0. (Default 0)
! |tz|                  The time zone you are in. This is 0 for GMT, 1 for MET
!                       and 2 for MEST (MET summer time). (Default 0)

        call addpar('date',-1.d0)
        call addpar('time',0.d0)
        call addpar('tz',0.d0)

!c------------------------------------------------------------------------

! DOCS	TERMS		Model parameters
!
! The next parameters define the inclusion or exclusion of
! certain terms of the primitive equations.
!
! |ilin|	Linearization of the momentum equations. If |ilin| 
!		is different from 0 the advective terms are not 
!		included in the computation. (Default 1)
! |itlin|	This parameter decides how the advective (non-linear)
!		terms are computed. The value of 0 (default) uses
!		the usual finite element discretization over a single
!		element. The value of 1 choses a semi-lagrangian
!		approach that is theoretically stable also for
!		Courant numbers higher than 1. It is however recommended
!		that the time step is limited using |itsplt| and
!		|coumax| described below. (Default 0)
! |iclin|	Linearization of the continuity equation. If |iclin|
!		is different from 0 the depth term in the continuity
!		equation is taken to be constant. (Default 0)

	call addpar('ilin',1.d0)
	call addpar('itlin',0.d0)
	call addpar('iclin',0.d0)

!c undocumented: rlin can lower effect of non linear terms

	call addpar('rlin',1.d0)

! The next parameters allow for a variable time step in the
! hydrodynamic computations. This is especially important for the
! non-linear model (|ilin=0|) because in this case the criterion
! for stability cannot be determined a priori and in any case the
! time integration will not be unconditionally stable.
!
! The variable time steps allows for longer basic time steps
! (here called macro time steps) which have to be set in |idt|.
! It then computes the optimal time step (here micro time step)
! in order to not exceed the given Courant number. However,
! the value for the macro time step will never be exceeded.
!
! Normally time steps are always given in full seconds. This is still
! true when specifying the macro time step |idt|. In older versions also
! the computed micro time steps also had to be full integer values.
! Starting from version 7.1 also fractional time steps are allowed.
! This gives the possibility to have time steps smaller than 1 s.
!
! |itsplt|	Type of variable time step computation. If this value
!		is 0, the time step will be kept constant at its initial
!		value. A value of 1 divides the initial time step into
!		(possibly) equal parts, but makes sure that at the end
!		of the micro time steps one complete macro time
!		step has been executed. The mode |itsplt| = 2
!		does not care about the macro time step, but always
!		uses the biggest time step possible. In this case
!		it is not assured that after some micro time steps
!		a macro time step will be recovered. Please note
!		that the initial macro time step will never be exceeded.
!		In any case, the time step will always be rounded to the
!		next lower integer value. This is not the case with
!		|itsplt| = 3 where the highest possible fractional time step 
!		will be used. (Default 0)
! |coumax|	Normally the time step is computed in order to not
!		exceed the Courant number of 1. However, in some cases
!		the non-linear terms are stable even for a value higher
!		than 1 or there is a need to achieve a lower Courant number.
!		Setting |coumax| to the desired Courant number
!		achieves exactly this effect. (Default 1)
! |idtsyn|	In case of |itsplt| = 2 this parameter makes sure that
!		after a time of |idtsyn| the time step will be syncronized
!		to this time. Therefore, setting |idtsyn| = 3600 means
!		that there will be a time stamp every hour, even if the model
!		has to take one very small time step in order to reach that
!		time. This parameter is useful
!		only for |itsplt| = 2 and its default value of
!		0 does not make any syncronization.
! |idtmin|	This variable defines the smallest time step possible
!		when time step splitting is enabled. Normally the smallest
!		time step is 1 second. Please set |idtmin| to values
!		smaller than 1 in order to allow for fractional time steps.
!		A value of 0.001 allows for timesteps of down to
!		1 millisecond. (Deault 1)
!c |idtmin|	This variable defines the smallest time step possible
!c		when time step splitting is enabled. Normally the smallest
!c		time step is 1 second. But when dealing with a lot of
!c		wet and drying in areas then sometimes it is useful to
!c		take out elements that limit the time step too much. In
!c		the case that |idtmin| is set to a value greater than 1
!c		the program will switch off temporarily elements that
!c		are responsible for such a small time step. (Default 0)

	call addpar('itsplt',0.d0)
	call addpar('coumax',1.d0)
	call addpar('idtsyn',0.d0)
	call addpar('idtmin',1.d0)
	call addpar('tfact',0.d0)		!still to comment FIXME

! These parameters define the weighting of time level in the 
! semi-implicit algorithm. With these parameters the damping
! of gravity or Rossby waves can be controlled. Only modify them if
! you know what you are doing.
!
! |azpar|	Weighting of the new time level of the transport
!		terms in the continuity equation. (Default 0.5)
! |ampar|	Weighting of the new time level of the pressure
!		term in the momentum equations. (Default 0.5)
! |afpar|	Weighting of the new time level of the Coriolis
!		term in the momentum equations. (Default 0.5)
! |avpar|	Weighting of the new time level of the non-linear
!		advective terms in the momentum equations. (Default 0.0)

	call addpar('azpar',0.5d0)
	call addpar('ampar',0.5d0)
	call addpar('afpar',0.5d0)
	call addpar('avpar',0.0d0)

! The next parameters define the weighting of time level for the
! vertical stress and advection terms. They guarantee the stability
! of the vertical system. For this reason they are normally set to
! 1 which corresponds to a fully implicit discretization. Only
! modify them if you know what you are doing.
!
! |atpar|	Weighting of the new time level of the vertical
!		viscosity in the momentum equation. (Default 1.0)
! |adpar|	Weighting of the new time level of the vertical
!		diffusion in the scalar equations. (Default 1.0)
! |aapar|	Weighting of the new time level of the vertical
!		advection in the scalar equations. (Default 1.0)

	call addpar('atpar',1.0d0)	!time weighting for vertical viscosity
	call addpar('adpar',1.0d0)	!time weighting for vertical diffusion
	call addpar('aapar',1.0d0)	!time weighting for vertical advection

!c------------------------------------------------------------------------

! DOCS	CORIOLIS	Coriolis parameters
!
! The next parameters define the parameters to be used
! with the Coriolis terms.
!
! |icor|	If this parameter is 0, the Coriolis terms are
!		not included in the computation. A value of 1
!		uses a beta-plane approximation with a variable
!		Coriolis parameter $f$, whereas a value of
!		2 uses an f-plane approximation where the
!		Coriolis parameter $f$ is kept constant over the
!		whole domain. (Default 0)
! |dlat|	Average latitude of the basin. This is used to
!		compute the Coriolis parameter $f$. This parameter 
!		is not used if spherical coordinates are used 
!		(|isphe|=1) or if a coordinate 	projection is set 
!		(|iproj| $>$0). (Default 0)
! |isphe|	If 0 a cartesian coordinate system is used,
!		if 1 the coordinates are in the spherical system (lat/lon).
!		Please note that in case of spherical coordinates the
!		Coriolis term is always included in the computation, even
!		with |icor| = 0. If you really do not want to use the
!		Coriolis term, then please set |icor| = -1. The default is
!		-1, which means that the type of coordinate system will 
!		be determined automatically.

	call addpar('icor',0.d0)
	call addpar('dlat',100.d0)
	call addpar('isphe',-1.d0)

!c------------------------------------------------------------------------

! DOCS	DEPTH		Depth parameters
!
! The next parameters deal with handling depth values of the basin.
!
! |href|	Reference depth. If the depth values of the basin and
!		the water levels are referred to mean sea level, |href|
!		should be 0 (default value). Else this value is
!		subtracted from the given depth values. For example,
!		if |href = 0.20| then a depth value in the basin
!		of 1 meter will be reduced to 80 centimeters.

	call addpar('href',0.d0)

! |hzmin|	Minimum total water depth that will remain in a
!		node if the element becomes dry. (Default 0.01 m)
! |hzoff|	Total water depth at which an element will be
!		taken out of the computation because it becomes dry.
!		(Default 0.05 m)
! |hzon|	Total water depth at which a dry element will be
!		re-inserted into the computation.
!		(Default 0.10 m)

	call addpar('hzmin',0.01d0)
	call addpar('hzoff',0.05d0)
	call addpar('hzon',0.10d0)

! |hmin|	Minimum water depth (most shallow) for the whole
!		basin. All depth values of the basin will be adjusted
!		so that no water depth is shallower than |hmin|.
!		(Default is no adjustment)
! |hmax|	Maximum water depth (deepest) for the whole
!		basin. All depth values of the basin will be adjusted
!		so that no water depth is deeper than |hmax|.
!		(Default is no adjustment)

	call addpar('hmin',-99999.d0)
	call addpar('hmax',+99999.d0)

!c------------------------------------------------------------------------

! \input{P_friction.tex}

	call addpar('ireib',0.d0)
	call addpar('czdef',0.d0)
	call addpar('iczv',1.d0)

!c------------------------------------------------------------------------

! DOCS	PHYSICAL		Physical parameters
!
! The next parameters describe physical values that can be adjusted
! if needed.
!
! |rowass|	Average density of sea water. (Default 1025 \densityunit)
! |roluft|	Average density of air. (Default 1.225 \densityunit)
! |grav|	Gravitational acceleration. (Default 9.81 \accelunit)

        call addpar('rowass',1025.d0)
        call addpar('roluft',1.225d0)
        call addpar('grav',9.81d0)

!c------------------------------------------------------------------------

! \input{P_wind.tex}

        call addpar('iwtype',1.d0)
	call addpar('itdrag',0.d0)
	call addpar('dragco',2.5d-03)
	call addpar('wsmax',50.d0)
	call addpar('wslim',-1.d0)

!c------------------------------------------------------------------------

! DOCS  meteo                      Meteo and heat flux parameters

! The next parameters deal with the heat and meteo forcing.

! |iheat|	The type of heat flux algorithm (Default 1):
!		\begin{description}
!		\item[1] As in the AREG model
!		\item[2] As in the POM model
!		\item[3] Following A. Gill
!		\item[4] Following Dejak
!		\item[5] As in the GOTM model
!		\item[6] Using the COARE3.0 module
!		\item[7] Read  sensible, latent and longwave fluxes from file
!		\item[8] Heat fluxes as Pettenuzzo et al., 2010
!		\end{description}

	call addpar('iheat',1.d0)		!type of heat flux routine

! |ihtype|	Different ways of how to specify water vapor content
!		are possible. Normally reletive humidity has to be
!		given (|ihtype|=1). However, also wet bulb temperature
!		(|ihtype|=2) or dew point temperature (|ihtype|=3) can
!		be given. (Default 1).

	call addpar('ihtype',1.d0)	!type of water vapor

! |isolp|       The type of solar penetration parametrization
!               by one or more exponential decay curves.
!               |isolp| = 0 sets an e-folding decay of radiation 
!               (one exponential decay curve) as function of depth |hdecay|.
!               |isolp| = 1 sets a profile of solar radiation with two length
!               scale of penetration. Following the Jerlov (Jerlov, N. G., 1968
!               Optical Oceanography, Elsevier, 194pp) classification the type
!               of water is clear water (type I). 
!               (Default 0) 

        call addpar('isolp',0.d0)         !type of solar penetration   

! |hdecay|	Depth of e-folding decay of radiation [m]. If |hdecay| = 0 
!		everything is absorbed in first layer (Default 0).

	call addpar('hdecay',0.d0)	!depth of e-folding decay of radiation

! |botabs|	Heat absorption at bottom [fraction] (Default 0).
!		\begin{description}
!		\item[=0] everything is absorbed in last layer
!		\item[=1] bottom absorbs remaining radiation
!		\end{description}

	call addpar('botabs',0.d0)	!heat absorption at bottom

! |albedo|	General albedo (Default 0.06).

	call addpar('albedo',0.06d0)	!general albedo

! |albed4|	Albedo for temp below 4 degrees (Default 0.06).

	call addpar('albed4',0.06d0)	!albedo for temp below 4 degrees

! |imreg| 	Regular meteo data (Default 0).

	call addpar('imreg',0.d0)		!regular meteo data - not used anymore

! |ievap| 	Compute evaporation mass flux (Default 0).

	call addpar('ievap',0.d0)		!compute evaporation mass flux

!c------------------------------------------------------------------------

! DOCS	3D			Parameters for 3d

! The next parameters deal with the layer structure in 3D.

! |dzreg|	Normally the bottom of the various layers are given in
!		section |\$levels|. If only a regular vertical grid is desired
!		then the parameter |dzreg| can be used. It specifies the spacing
!		of the vertical layers in meters. (Default is 0, which means 
!		that the layers are specified explicitly in |\$levels|.

	call addpar('dzreg',0.d0)		!regular vertical grid

! The last layer (bottom layer) is treated in a special way. Depending on
! the parameter |ilytyp| there are various cases to be considered. A value
! of 0 leaves the last layer as it is, even if the thickness is very small.
! A value of 1 will always eliminate the last layer, if it has not full
! layer thickness. A value of 2 will do the same, but only if the last layer
! is smaller than |hlvmin| (in units of fraction). Finally, a value of
! 3 will add the last layer to the layer above, if its layer thickness
! is smaller than |hlvmin|.
!
! |ilytyp|	Treatment of last (bottom) layer. 0 means no adjustment,
!		1 deletes the last layer, if it is not a full layer,
!		2 only deletes it
!		if the layer thickness is less than |hlvmin|, and 3
!		adds the layer thickness to the layer above if it is smaller
!		than |hlvmin|. Therefore, 1 and 2 might change the
!		total depth and layer structure, while 3 only might
!		change the layer structure. The value of 1 will always
!		give you full layers at the bottom. (Default 3)
! |hlvmin|	Minimum layer thickness for last (bottom) layer used when
!		|ilytyp| is 2 or 3. The unit is fractions of the nominal
!		layer thickness. Therefore, a value of 0.5 indicates that
!		the last layer should be at least half of the full
!		layer. (Default 0.25)

	call addpar('ilytyp',3.00d0)	!type of depth adjustment
	call addpar('hlvmin',0.25d0)	!min percentage of last layer thickness

!c not yet documented features of sigma layers

	call addpar('nsigma',0.d0)	!number of sigma layers
	call addpar('hsigma',10000.d0)	!lower depth of sigma layers (hybrid)

!c baroclinic model
!c
!c 0 no   1 full    2 diagnostic    3 T/S advection but no baroclinic terms

	call addpar('ibarcl',0.d0)	!compute baroclinic contributions ?

	call addpar('iturb',0.d0)		!use turbulence closure scheme

! The next parameters deal with vertical diffusivity and viscosity.

! |diftur|	Vertical turbulent diffusion parameter for the tracer.
!		(Default 0)
! |difmol|	Vertical molecular diffusion parameter for the tracer.
!		(Default 1.0e-06)
! |vistur|	Vertical turbulent viscosity parameter for the momentum.
!		(Default 0)
! |vismol|	Vertical molecular viscosity parameter for the momentum.
!		(Default 1.0e-06)

        call addpar('difmol',1.0d-06)	!molecular vertical diffusivity
	call addpar('diftur',0.d0)	!diffusion parameter (vertical), cvpar
        call addpar('vismol',1.0d-06)	!molecular vertical viscosity
        call addpar('vistur',0.d0)	!turbulent vertical viscosity (nau)

!c------------------------------------------------------------------------

!c horizontal diffusion (Smagorinsky)

	call addpar('idhtyp',0.d0)	!type of horizontal diffusion/viscosity
	call addpar('dhlen',1.d0) 	!length scale for idhtyp=1
	call addpar('noslip',0.d0) 	!no slip conditions on boundary
	call addpar('idtype',2.d0) 	!type of hor diffus (delete after test)

	call addpar('ahpar',0.d0)		!austausch coefficient (still same??)

! |dhpar|	Horizontal diffusion parameter (general).
!		(Default 0)

	call addpar('dhpar',0.d0)		!diffusion parameter

! The next parameters deal with the control of the scalar transport 
! and diffusion equation. You have possibility to prescribe the tvd scheme
! desired and to limit the Courant number.
!
! |itvd|	Type of the horizontal advection scheme used for 
!		the transport and diffusion
!		equation. Normally an upwind scheme is used (0), but setting
!		the parameter |itvd| to a value greater than 0 
!		choses a TVD scheme. A value of 1 will use a TVD scheme
!		based on the average gradient, and a value of 2 will use
!		the gradient of the upwind node (recommended).
!		This feature
!		is still experimental, so use with care. (Default 0)
! |itvdv|	Type of the vertical advection scheme used for 
!		the transport and diffusion
!		equation. Normally an upwind scheme is used (0), but setting
!		the parameter |itvd| to 1 choses a TVD scheme. This feature
!		is still experimental, so use with care. (Default 0)
! |rstol|	Normally the internal time step for scalar advection is
!		automatically adjusted to produce a Courant number of 1
!		(marginal stability). You can set |rstol| to a smaller value 
!		if you think there are stability problems. (Default 1)

	call addpar('itvd',0.d0)		!horizontal TVD scheme?
	call addpar('itvdv',0.d0)		!vertical TVD scheme?
	call addpar('rstol',1.d0)		!limit time step to this Courant number

!c------------------------------------------------------------------------
!c Equation of State Unesco or JM95
! Use different equation of state, default is Unesco (eostype=1), the other
! choice is Jackett & McDougall 1995 (eostype=2)
	call addpar('eostype',1.0d0)		! EOS


!c------------------------------------------------------------------------
! DOCS	VAR		Various parameters
!
! The next parameters describe various parameters not related to
! the above parameters.

! |tauvel|	If you have velocity observations given in file
!		|surfvel| then you can specify the relaxation
!		parameter $\tau$ in the variable |tauvel|. (Default 0,
!		which means no assimilation of velocities)

	call addpar('tauvel',0.d0)	!horizontal TVD scheme?

! |rtide|	If |rtide| = 1 the model calculates equilibrium tidal 
!		potential and load tides and uses these to force the 
!		free surface (Default 0).

	call addpar('rtide',0.d0)

! |ltidec|      Calibration factor for calculating the loading tide,
!               which is computed in function of the total water depth as
!               $\beta=ltidec*H$. Usually it has a value of order 1e-6.
!               If 0 no loading tide is computed (Default 0).

        call addpar('ltidec',0.d0)


!c------------------------------------------------------------------------

! DOCS	ST	Temperature and salinity
!
! The next parameters deal with the transport and diffusion
! of temperature and salinity. Please note that in order to compute
! T/S the parameter |ibarcl| must be different from 0. In this case
! T/S advection is computed, but may be selectively turned off setting
! one of the two parameters |itemp| or |isalt| explicitly to 0.
!
! |itemp|	Flag if the computation on the temperature is done.
!		A value different from 0 computes the transport
!		and diffusion of the temperature. (Default 1)
! |isalt|	Flag if the computation on the salinity is done.
!		A value different from 0 computes the transport
!		and diffusion of the salinity. (Default 1)

	call addpar('itemp',1.d0)		!compute temperature ?
	call addpar('isalt',1.d0)		!compute salinity ?
	call addpar('irho',0.d0)		!write rho ?

! The next parameters set the initial conditions for temperature and salinity.
! Both the average value and and a stratification can be specified.
!
! |temref|	Reference (initial) temperature of the water in
!		centigrade. (Default 0)
! |salref|	Reference (initial) salinity of the water in
!		psu (practical salinity units) or ppt.
!		(Default 0)
! |tstrat|	Initial temperature stratification in units of [C/km].
!		A positive value indicates a stable stratification.
!		(Default 0)
! |sstrat|	Initial salinity stratification in units of [psu/km].
!		A positive value indicates a stable stratification.
!		(Default 0)

	call addpar('temref',0.d0)	!reference temperatur for baroc runs
	call addpar('salref',0.d0)	!reference salinity for baroc runs

	call addpar('sstrat',0.d0)	!salt stratification
	call addpar('tstrat',0.d0)	!temp stratification

        call addpar('imellor',0.d0)       !if > 0 T,S initialized to get
                                        !stratification as in Mellor,1993                    

! The next parameters deal with horizontal diffusion of temperature
! and salinity. These parameters overwrite the general parameter for
! horizontal diffusion |dhpar|.

! |thpar|	Horizontal diffusion parameter for temperature.
!		(Default 0)
! |shpar|	Horizontal diffusion parameter for salinity.
!		(Default 0)

	call addpar('thpar',-1.d0)	!horiz. diff. coeff. for temp.
	call addpar('shpar',-1.d0)	!horiz. diff. coeff. for sal.

!c------------------------------------------------------------------------

! DOCS	CC	Concentrations
!
! The next parameters deal with the transport and diffusion
! of a conservative substance. The substance is dissolved in
! the water and acts like a tracer.
!
! |iconz|	Flag if the computation on the tracer is done.
!		A value different from 0 computes the transport
!		and diffusion of the substance. If greater than 1
!		|iconz| concentrations are simulated. (Default 0)
! |conref|	Reference (initial) concentration of the tracer in
!		any unit. (Default 0)
! |contau|	Decay rate for concentration if different from 0. In
!		this case |contau| is the decay rate (e-folding time) in days.
!		There also is the possibility to set different decay rates
!		for multi-concentration runs. In this case the value of
!		|taupar| has to be adjusted in the program code.
!		(Default 0)
! |idecay|	Type of decay used. If 0 no decay is used.
!		A value of 1 uses the value of |contau| as exponential decay.
!		A value of 2 uses a formulation of Chapra, where the
!		decay rate depends on T,S,light and settling. In this
!		case the value of |contau| is ignored.
!		(Default 0)

	call addpar('iconz',0.d0)		!compute concentration ?
	call addpar('conref',0.d0)	!reference concentration
	call addpar('contau',0.d0)	!decay rate [days]
	call addpar('idecay',0.d0)	!type of decay

! |chpar|	Horizontal diffusion parameter for the tracer.
!		This value overwrites the general parameter for
!		horizontal diffusion |dhpar|. (Default 0)

	call addpar('chpar',-1.d0)	!diffusion parameter

!c------------------------------------------------------------------------

! DOCS	STOUPUT	Output for scalars
!
! The next parameters define the output frequency of the
! computed scalars (temperature, salinity, generic concentration) to file.
!
! |idtcon|, |itmcon|	Time step and start time for writing to file
!			CON (concentration) and NOS (temperature and
!			salinity).

	call addpar('idtcon',0.d0)	!time step for output
	call addpar('itmcon',-1.d0)	!minimum time for output

! DOCS	END

!c------------------------------------------------------------------------
!c still to be commented below here
!c------------------------------------------------------------------------

	call addpar('idtsti',0.d0)	!time step for stability index
	call addpar('itmsti',-1.d0)	!minimum time for stability index

!c------------------------------------------------------------------------

!c biological reactor

	call addpar('ibio',0.d0)		!run biological reactor

!c toxicological routines from ARPAV

	call addpar('itoxi',0.d0)		!run toxicological routines

!c------------------------------------------------------------------------

!c call routines in flux.f

	call addpar('regflx',0.d0)	!compute regular fluxes

!c------------------------------------------------------------------------

!c custom call

	call addpar('icust',0.d0)		!call custom routine
	call addpar('tcust',0.d0)		!time for custom routine

	call addpar('ishyff',0.d0)	!shyfem file format
					!0=old 1=new 2=both

	call addpar('ipoiss',0.d0)	!solve poisson equation

!c rain

	call addpar('zdist',0.d0)		!distributed water level

	call addpar('idtbox',0.d0)	!for boxes
	call addpar('itmbox',-1.d0)

	call addpar('ihwadv',0.d0)        !vert advection of horiz momentum

!c ice
	call addpar('iaicef',-99.d0)	!area code for ice free condition
        call addpar('ilockex',0.d0)       ! ivb
                                        ! initialize temperature
                                        ! for lock exchange experiment
        call addpar('iintwave',0.d0)      ! initialize temperature
                                        ! for interanal wave generation
                                        ! test case
                        
!-------------------------------------------------------------------------

! DOCS	TIME_DATE	General partitioning parameters
!
! A parameter can be specified that represents the number of partitions in which
! the domain will be distributed

        call addpar('nn_nparts',1.d0)      ! default 1 partitions (mpi subdomains)
        call addpar('nn_wvert',0.d0)   ! default: the weight is the same for all vertices

!-------------------------------------------------------------------------
! parameter for activate timing
        call addpar('ln_timing',0.d0)      ! default timing is off

	end

!

	subroutine nlsinh_lagrg

	implicit none

! $lagrg section

! DOCS	START	P_lagrg
!
! This section describes the use of the Lagrangian Particle Module.
! The lagrangian particles can be released inside a specified area with a
! regular distribution. The area is defined in the file |lagra|.
! The amount of particles released and the
! time step is specified by nbdy and idtrl.
! 
! The lagrangian module runs between the times |itlanf| and |itlend|. If one
! or both are missing, the simulation extremes are substituted. Inside
! the lagrangian simulation window the release of particles is controlled
! by the parameters |idtl|, |itranf| and |itrend|. |itranf| gives the time
! of the first release, |itrend| the time for the last release. If not
! given they are set equal to the extremes of the lagrangian simulation.
! |idtl| is giving the time step of release.
!
! Particles are released inside the given areas (filename |lagra|). If
! this file is not specified they are released over the whole domain. There is
! also a possibility release particles over open boundaries. However,
! this is still experimental. Please see the file |lagrange_main.f|
! for more details.
!
! The output frequency of the results can be contolled by 
! |idtlgr| and |itmlgr|.
!
! Please find all details here below.

        call sctpar('lagrg')             !sets default section
        call sctfnm('lagrg')

! |ilagr|	Switch that indicates that the lagrangian module
!		should be run (default 0):
!		\begin{description}
!		\item[0] do nothing
!		\item[1] surface lagrangian
!		\item[2] 2d lagrangian (not implemented)
!		\item[3] 3d lagrangian
!		\end{description}

	call addpar('ilagr',0.d0)         !LAGR

! |nbdymax|	Maximum numbers of particles that can be in the domain.
!		This should be the maximum number of particles
!		that can be created and inserted. Use 0 to not limit
!		the number of particles (on your own risk). This
!		parameter must be set and has no default.

        call addpar('nbdymax',-1.d0)

! |nbdy|	Total numbers of particles to be released in the domain each
!		time a release of particles takes place. 
!		(Default 0)

        call addpar('nbdy',0.d0)

! |rwhpar|	A horizontal diffusion can be defined for the lagrangian model.
!		Its value can be specified in |rwhpar| and the units are 
!		[m**2/s]. (Default 0)

	call addpar('rwhpar',0.d0)	!diffusion for lagrangian model

! |itlanf, itlend|	The start and end time for the lagrangian module.
!			If not given, the module runs for the whole simulation.

        call addpar('itlanf',-1.d0)
        call addpar('itlend',-1.d0)

! |itmlgr, idtlgr|	Initial time and time step for the output to file
!			of the particles. if |idtlgr| is 0, 
!			no output is written. (Default 0)

        call addpar('itmlgr',-1.d0)
        call addpar('idtlgr',0.d0)

! |idtl|	The time step used for the release of particles. If this
!		is 0 particles are released only once at the beginning
!		of the lagrangian simulation. No particles are released 
!		for a value of less than 0. (Default 0)

        call addpar('idtl',0.d0)

! |itranf, itrend|	Initial and final time for the release of particles.
!			If not specified the particles are released over the
!			whole lagrangian simulation period.

        call addpar('itranf',-1.d0)
        call addpar('itrend',-1.d0)

! |ipvert|	Set the vertical distribution of particles:
!		\begin{description}
!		\item[0] releases one particles only in surface layer
!		\item[$>$0] release n particles regularly
!		\item[$<$0] release n particles randomly
!		\end{description}

        call addpar('ipvert',0.d0)

! |linbot| Set the bottom layer for vertical releases

         call addpar('linbot',0.d0)

!c |lintop| Set the top layer for vertical releases

!c        call addpar('lintop',0.) !todo

! |lagra|	File name that contains closed lines of the area where
!		the particles have to be released. If not given, the particles
!		are released over the whole domain.

        call addfnm('lagra',' ')

!c To compute the transit time of the particles to leave
!c a specified area. 'Artype' is the flag to detect the
!c element outside the defined area. 
!c Default = -1, i.e., no transit times are computed

        call addpar('artype',-1.d0)

!c still to be commented

        call addpar('lcust',0.d0)
        call addpar('ldecay',0.d0)
        call addpar('ioil',0.d0)
        call addpar('ilarv',0.d0)
        call addpar('ised',0.d0)

! DOCS	END

	end

!************************************************************************

        subroutine nlsinh_wrt

        implicit none

! $wrt section

! DOCS  START   P_wrt
!
! Parameters for computing water renewal time.
! During runtime if writes a .jas file with timeseries of total tracer
! concentration in the basin and WRT computed according to different methods.
! Nodal values of computed WRT are written in the .wrt file.
! Frequency distrubutions of WRTs are written in the .frq file.
!
! Please find all details here below.

        call sctpar('wrt')             !sets default section
        call sctfnm('wrt')

! |idtwrt|	Time step to reset concentration to c0. Use 0 if no reset
!		is desired. Use -1 if no renewal time computation is desired
!		(Default -1).

        call addpar('idtwrt',-1.d0)

! |itmin|	Time from when to compute renewal time (-1 for start of sim)
!		(Default -1)

        call addpar('itmin',-1.d0)

! |itmax|	Time up to when to compute renewal time (-1 for end of sim)
!		(Default -1).

        call addpar('itmax',-1.d0)

! |c0|		Initial concentration of tracer (Default 1).

        call addpar('c0',1.d0)

! |iaout|	Area code of elements out of lagoon (used for init and retflow).
!		Use -1 to if no outside areas exist. (Default -1).

        call addpar('iaout',-1.d0)

! |percmin|	Percentage to reach after which the computation is stopped.
!		Use 0 if no premature end is desired (Default 0).

        call addpar('percmin',0.d0)

! |iret|	Equal to 1 if return flow is used. If equal to 0 the 
!		concentrations outside are explicitly set to 0 (Default 1).

        call addpar('iret',1.d0)

! |istir|	If equal to 1 simulates completely stirred tank 
!		(replaces at every time step conz with average conz)
!		(Default 0).

        call addpar('istir',0.d0)

! |iadj|	Adjust renewal time for tail of distribution (Default 1).

        call addpar('iadj',1.d0)

! |ilog|	Use logarithmic regression to compute renewal time (Default 0).

        call addpar('ilog',0.d0)

! |ctop|	Maximum to be used for frequency curve (Default 0).

        call addpar('ctop',0.d0)

! |ccut|	Cut renewal time at this level (for res time computation)
!		(Default 0).

        call addpar('ccut',0.d0)

! |wrtrst|	If reset times are not regularly distributed (e.g., 1 month)
!		it is possible to give the exact times when a reset should
!		take place. |wrtrst| is a file name where these reset times
!		are specified, one for each line. For every line two integers
!		indicating date and time for the reset must be specified.
!		If only one value is given, time is taken to be 0. The format
!		of date is "YYYYMMDD" and for time "hhmmss". If the file
!		wrtrst is given |idtwrt| should be 0.

        call addfnm('wrtrst',' ')

! DOCS  END

        end

!************************************************************************

	subroutine nlsinh_bfmsc

! Parameters for BFM ECOLOGICAL module        !BFMSC

	implicit none

        call sctpar('bfmsc')             !sets default section
        call sctfnm('bfmsc')

! BFM ECOLOGICAL MODEL CALL 

	call addpar('ibfm',0.d0)
	call addpar('ibtanf',0.d0)
	call addpar('ibtend',0.d0)
	call addpar('itmbfm',-1.d0)
	call addpar('idtbfm',0.d0)
	call addpar('bligth',1.d0) !light flag=1 max/min light W/m**2 in nml file

	end

!************************************************************************

	subroutine nlsinh_proj

	implicit none

! $proj section

! DOCS  START   P_proj
!
! The parameter sets in this section handle the projection from cartesian to
! geographical coordinate system. If |iproj| $>$0 the projected geographical 
! coordinates can be used for computing spatially variable Coriolis parameter 
! and tidal potential even if the basin is in cartesian coordinate system 
! (|isphe| = 0) .
!
! Please find all details here below.

        call sctpar('proj')             !sets default section
        call sctfnm('proj')

! |iproj|	Switch that indicates the type of projection
!		(default 0):
!		\begin{description}
!		\item[0] do nothing
!		\item[1] Gauss-Boaga (GB)
!		\item[2] Universal Transverse Mercator (UTM)
!		\item[3] Equidistant cylindrical (CPP)
!		\item[4] UTM non standard
!		\end{description}

	call addpar('iproj',0.d0)

! |c\_fuse|	Fuse for GB (1 or 2, default 0)

	call addpar('c_fuse',0.d0)

! |c\_zone|	Zone for UTM (1-60, default 0)

	call addpar('c_zone',0.d0)

! |c\_lamb|	Central meridian for non-std UTM (default 0)

	call addpar('c_lamb',0.d0)

! |c\_x0|	x0 for GB and UTM (default 0)

	call addpar('c_x0',0.d0)

! |c\_y0|	y0 for GB and UTM (default 0)

	call addpar('c_y0',0.d0)

! |c\_skal|	Scale factor for non-std UTM (default 0.9996)

	call addpar('c_skal',0.9996d0)

! |c\_phi|	Central parallel for CPP (default 0.9996)

	call addpar('c_phi',0.d0)

! |c\_lon0|	Longitude origin for CPP (default 0)
	call addpar('c_lon0',0.d0)

! |c\_lat0|	Latitude origin for CPP (default 0)
	call addpar('c_lat0',0.d0)

! DOCS  END

	end

!************************************************************************

	subroutine nlsinh_undoc

! not documented parameters 

	implicit none

	call sctpar('para')		!sets default section

!c undocumented parameters

	call addpar('iclose',0.d0)
	call addpar('itsmed',0.d0)	!averages for T/S

!c internally used parameters

	call addpar('flag',0.d0)
	call addpar('const',0.d0)
	call addpar('volmin',1.d0)	!minimum volume to remain in el.

!c debug

	call addpar('vreps',5.d-4)
	call addpar('vrerr',1.d-1)
	call addpar('levdbg',0.d0)	!debug level (0 = no, 9 = max)

!c distance for advective terms

	call addpar('nadist',0.d0)

!c new for scaling time step

call addpar('itunit',1.d0)

!c experimental stuff

        !call addpar('nbsig',0.)         !sigma layers to read in for OBC

        call addpar('sedim',0.d0)         !sedimentation for theseus
        call addfnm('hsedim',' ')         !sedimentation hev file for theseus

        call addpar('nomp',0.d0)          !number of threads to use

	end

!************************************************************************

	subroutine nlsinh_georg

! parameters used for private projects

	implicit none

	call sctpar('para')		!sets default section

	call addpar('hlido',999.d0)	!maximum depth at lido
	call addpar('hmala',999.d0)	!maximum depth at lido
	call addpar('hchio',999.d0)	!maximum depth at lido

	call addpar('zrise',0.d0)		!sea level rise
	call addpar('zsalv',999.d0)	!saveguarding level
	call addpar('zfranc',0.d0)	!extra security for forecast

	end

!************************************************************************

	subroutine nlsinh_unused

! parameters not used anymore -> to be deleted

	implicit none

	call sctpar('para')		!sets default section

	call addpar('iprogr',0.d0)

	end

! ********************************************************************

! This subroutine defines the simulation wave parameter

        subroutine nlsinh_waves

        implicit none

! DOCS  START   P_wave
!
! The following parameters activate the wind wave module and define
! which kind of wind wave model has to be used.
!c |waves|      Wave module section name.

        call sctpar('waves')             !sets waves section
        call sctfnm('waves')

! |iwave|	Type of wind wave model and coupling procedure (default 0):
!		\begin{description}
!		\item[0] No wind wave model called 
!		\item[1] The parametric wind wave model is called (see file subwave.f)
!		\item[$>=$2] The spectral wind wave model WWMIII is called
!		\item[2] ... wind from SHYFEM, radiation stress formulation
!		\item[3] ... wind from SHYFEM, vortex force formulation 
!		\item[4] ... wind from WWMIII, radiation stress formulation
!		\item[5] ... wind from WWMIII, vortex force formulation 
!		\end{description}
!		When the vortex force formulation is chosen the wave-supported
!		surface stress is subtracted from the wind stress, in order to
!		avoid double counting of the wind forcing in the flow model.
!		Moreover, the use of the wave-depended wind drag coefficient 
!		could be adopted setting |itdrag| = 3.

        call addpar('iwave',0.d0)

!
! |dtwave|	Time step for coulping with WWMIII. Needed only for
!		|iwave| $>$ 1 (default 0).

        call addpar('dtwave',0.d0)

! |idtwav|, |itmwav|	Time step and start time for writing to file wav,
!			the files containing output wave variables (significant 
!			wave height, wave period, mean wave direction).

        call addpar('idtwav',0.d0)
        call addpar('itmwav',-1.d0)

!
! DOCS  END
!

        end

!************************************************************************

	subroutine nlsinh_nonhydro

! parameters for non hydrostatic model (experimental)

	implicit none

        call sctpar('nonhyd')             !sets default section

	call addpar('inohyd',0.d0)	!for non-hydrostatic model

        call addpar('aqpar',0.5d0)

	call addpar('islope',0.d0)	!type of grid for poisson equation

        call addpar('ivwadv',0.d0)        !vert advection of vert momentum
        call addpar('inhflx',0.d0)        !flux upwind for horiz advect of w
	call addpar('inhadj',0.d0)        !choice for correction of U,V,eta
        call addpar('inhbnd',0.d0)        !exclude NH dynamics for boundaries
        call addpar('iwvel',1.d0)         !write vertical velocity
        call addpar('iqpnv',1.d0)         !write NH pressure

	end

!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************

	subroutine nlsina

! initializes the parameter file for the post processing routines

	implicit none

	call nlsina_general
	call nlsina_para
	call nlsina_color
	call nlsina_arrow
	call nlsina_legend
	call nlsina_legvar
	call nlsina_sect

	end

!************************************************************************

	subroutine nlsina_general

! $para section with general variables

	implicit none

! DOCS	START	S_para_general
!
! The next parameters decide what to plot. They can be set directly in
! the STR file here. However, the prefered way to set them is through
! the back end |plots| that automatically sets these parameters.

	call sctpar('para')		!sets default section
	call sctfnm('para')		!sets default section

! The parameter |iwhat| is compulsory and defines the variable to be plotted.
! Please note again that this parameter will be set if invoking the plot
! program through |plots|.
!
! |iwhat|		Flag that determines what to plot. If 0 then
!			the program asks interactively for it. (Default 0)
!			\begin{description}
!			\item[1] Plot basin (grid and isolines of depth)
!			\item[2] Plot velocities
!			\item[3] Plot transports
!			\item[4] Plot water levels
!			\item[5] Plot concentration
!			\item[6] Plot temperature
!			\item[7] Plot salinity
!			\item[8] Plot rms-velocity
!			\item[9] (not used)
!			\item[10] Plot generic scalar (see |ivar|)
!			\item[11] Plot wind vectors
!			\item[12] Plot lagrangian particles
!			\item[13] Plot wave data
!			\end{description}

	call addpar('iwhat',0.d0)		!what to plot

! The next parameters define other aspects of the plot. In particular their
! meaning is
!
! |iauto|		Normally the simulation name and basin are
!			asked interactively when running the plotting
!			program. However, if |iauto| is different from 0
!			then the file |.memory| is read and the information
!			contained in it is used. The file memory can be
!			set through the program |memory| and it is
!			changed when other parameters are inputted
!			interactively.
! |level|		For 3D applications it indicates the vertical level
!			for which the plot is desired. 1 indicates the 
!			surface layer, 2 the second layer from the surface etc.
!			0 gives integrated quantities, and a value of -1
!			indicates the bottom layer (which refers to different
!			layers for every element). (Default 0)
! |ivar|		Variable to be plotted for scalar quantities. In the
!			file NOS more then one scalars can be contained. To
!			choose the desired scalar, the id of the scalar has 
!			to be given in |ivar|. Note that if only one variable
!			type is contained in the file, then it is not
!			necessary to set |ivar|.

	call addpar('iauto',0.d0)		!silent mode
	call addpar('level',0.d0) 	!level (-1 -> bottom   0 -> integr.)
	call addpar('ivar',0.d0)		!what variable to plot

!c still to document

	call addfnm('varnam',' ')	!variable name to be plotted
	call addpar('isect',0.d0)		!vertical section

! The next variables define the time of the plots. Even if the names of two
! of the variables are identical to variables used in the finite element
! model, the meaning is different. In the simulation model, |itanf, itend|
! define the initial and final time of simultaion, whereas here these
! variables define the initial and final time of the plot of results.
! if none of the time variables are set, then all time records are plotted.

! |itanf|		Initial time for plotting. (Default is 
!			first data record)
! |itend|		Final time for plotting. (Default is 
!			last data record)
! |nout|		Frequency for plotting. A value of 0 or 1
!			plots every record, 2 plots every other record etc.
!			A negative value works as a positive one,
!			but starts plotting from the first record. So
!			-2 plots records 1,3,5,etc.
!			(Default 1) 

	call addpar('itanf',-1.d0)	!time start
	call addpar('itend',-1.d0)	!time end
	call addpar('nout',1.d0)		!time frequence

	call addpar('atanf',-1.d0)	!time start (absolute)
	call addpar('atend',-1.d0)	!time end (absolute)

! DOCS	END

	end

!************************************************************************

	subroutine nlsina_para

! $para section

	implicit none

! DOCS	START	S_para_a
!
! These parameters set generic values for the plot.

	call sctpar('para')		!sets default section
	call sctfnm('para')		!sets default section

! Some of the parameters set coordinates in the plot. For example, the
! values |x0, y0| and |x1, y1| indicate the actual plotting area, which can
! be bigger or smaller than the extension of the numerical grid.
!
! Normally, values have to be in meters (the same as the coordinates in the
! numerical grid). However, also relative coordinates can be used. If all
! values given are in the range between -1 and +2, then these values
! are interpreted as relative coordinates. Therefore, $x$ coordinates of
! 0 indicate the left border, and 1 the right border. The upper left quarter
! of the domain can be chosen with (|x0, y0|) = (0,0.5) and
! (|x1, y1|) = (0.5,1).

! |x0, y0|		Lower left corner of the plotting area.
!			(Default is whole area)
! |x1, y1|		Upper right corner of the plotting area.
!			(Default is whole area)

	call addpar('x0',0.d0)		!dimension of plot
	call addpar('y0',0.d0)		!dimension of plot
	call addpar('x1',0.d0)		!dimension of plot
	call addpar('y1',0.d0)		!dimension of plot

! The next values give the position, where the legend (scale bar and
! true north) is plotted. This legend will only be plotted if
! the coordinates are not geographical (lat/lon) but cartesian.

! |x0leg, y0leg|	Lower left corner of the area
!			where the legend is plotted.
! |x1leg, y1leg|	Upper right corner of the area.
!			where the legend (north and scale) is plotted.

	call addpar('x0leg',0.d0)		!dimension of legend
	call addpar('y0leg',0.d0)		!dimension of legend
	call addpar('x1leg',0.d0)		!dimension of legend
	call addpar('y1leg',0.d0)		!dimension of legend

! |lblank|		The legend is plotted over a white rectangle.
!			Sometimes this blanking is not desireable.
!			If you do not want to have a white box below the
!			legend set |lblank| to 0. (Default 1)

	call addpar('lblank',1.d0)

! |cislnd|		It is possible to plot all islands in gray color.
!			Setting |cislnd| to a value between 0 (black) and 
!			1 (white) will achieve this. A negative value
!			will not fill islands with gray color.
!			(Default -1)

	call addpar('cislnd',-1.d0)	!plot also outer island: 2 <= c <= 3

! |dgray|		It is possible to plot all dry areas in gray color.
!			Setting |dgray| to a value between 0 (black) and 
!			1 (white) will achieve this. A negative value
!			will not fill dry areas with gray color.
!			(Default -1)

! |hgray|		Whereas |dgray| is normally only coloring
!			elements that are dry, you can also color elements
!			shallower than a given depth |hgray|. E.g., a value
!			for |hgray| of -0.5 will plot in gray all
!			elements with depth lower than -0.5 m (salt
!			marshes). (Default -10000)

	call addpar('dgray',-1.d0)
	call addpar('hgray',-10000.d0)	!gray all elems with h < hgray

! |dxygrd|		Grid size if the results are interpolated on
!			a regular grid. A value of 0 does
!			not use a regular grid but the original
!			finite element grid for plotting. (Default 0)
! |typls|		Typical length scale to be used when scaling
!			velocity or transport arrows. If |dxygrd| is
!			given this length is used and |typls| is not used.
!			If not given it is computed from the basin
!			parameters. (Default 0)
! |typlsf|		Additional factor to be used with typls to
!			determine the length of the maximum or
!			reference vector. This is the easiest way
!			to scale the velocitiy arrows 
!			with an overall factor. (Default 1)
! |velref|		Reference value to be used when scaling arrows.
!			If given, a vector with this value will have a length
!			of |typls|*|typlsf| on the map, or, in case
!			|dxygrd| is given, |dxygrd|*|typlsf|. If not set
!			the maximum value of the velocity/transport
!			will be used as |velref|. (Default 0)
! |velmin|		Minimum value for which an arrow will be plotted.
!			With this value you can eliminate small arrows
!			in low dynamic areas. (Default 0)

	call addpar('dxygrd',0.d0)	!grid size for regular grid
	call addpar('typls',0.d0)		!typical length scale for arrow plot
	call addpar('typlsf',1.d0)	!factor for typical length scale
	call addpar('velref',0.d0)	!reference velocity for length scale
	call addpar('velmin',0.d0)	!minimum velocity to be plotted

! |isphe|		If 0 a cartesian coordinate system is used,
!			If 1 the coordinates are in the spherical 
!			system (lat/lon). Among other, this
!			indicates that the $x$-coordinates will be multiplied
!			by a factor that accounts for the visual deformation
!			using lat/lon coordinates.
!			The default is -1, which means that the 
!			type of coordinate system will 
!			be determined automatically. (Default -1)
! |reggrd|		If different from 0 it plots a regular grid over
!			the plot for geographical reference. The value of
!			|reggrd| gives the spacing of the regular grid lines.
!			The units must be according to the units used for
!			the coordinates. With value of -1 the regular grid is
!			determined automatically. (Default -1)
! |regdst|		This value gives the number of intervals
!			that are used to sub-divide the grid given by
!			|reggrd| with a black and white scale around
!			the plot. If 0 it tries to determine automatically
!			the sub-intervals (2 or 4). A value of -1 does
!			not plot the subgrid scale. (Default 0)
! |reggry|		If plotting the regular overlay grid this gives
!			the gray value used for the grid. 0 is black, and
!			1 is white. A value of 1 does not plot the
!			overlay grid, but still writes the labels. 
!			(Default 1)

	call addpar('isphe',-1.d0)	!spherical coordinate system
	call addpar('reggrd',-1.d0)	!regular grid spacing
	call addpar('regdst',0.d0)	!regular micro grid spacing
	call addpar('reggry',1.d0)	!gray value

! |bndlin|		Name of file that gives the boundary line
!			that is not part of the finite element domain.
!			The file must be in BND format. You can use
!			the program grd2bnd.pl to create the file
!			from a GRD file.
!			(Default is no file)

	call addfnm('bndlin'," ")	!name of boundary line file

! |ioverl|		Create overlay of velocity vectors on scalar value.
!			With the value of 0 no overlay is created, 1
!			creates an overlay with the velocity speed.
!			The value of 2 overlays vertical velocities
!			3 water levels and 4 overlays bathymetry.(Default 0)
! |inorm|		Normally the horizontal velocities are plotted
!			in scale. The value of |inorm| can change this
!			behavior. A value of 1 normalizes velocity vectors
!			(all vectors are the same length), whereas 2
!			scales from a given minimum velocity |velmin|.
!			Finally, the value of 3 uses a logarithmic scale.
!			(Default 0)

	call addpar('ioverl',0.d0)	!overlay in color
	call addpar('inorm',0.d0)		!vertical velocity as overlay

! The next parameters give the choice to selectively avoid to plot areas
! of the basin and to apply different gray tones for the boundary and
! net lines when plotting the basin.
! Please remember that when working with gray tones the value should
! be between 0 (black) and 1 (white).
!
! |ianopl|	Area code for which no plot has to be produced. Normally 
!		the whole basin is plotted, but with this parameter some
!		areas can be excluded. (Default -1)
!		the bathymetry. (Default 0.8)
! |bgray|	Gray value used for the finite element grid when plotting
!		the bathymetry. (Default 0.8)
! |bbgray|	Gray value used for the boundary of the finite element grid.
!		(Default 0)
! |bsgray|	Gray value used to plot the finite element grid over
!		a scalar or velocity plot. This is basically useful
!		for debugging reasons. The default is to not plot
!		the grid (Default -1.0)

        call addpar('ianopl',-1.d0)      !do not plot these areas
	call addpar('bgray',0.8d0)       !gray value for bathymetry
	call addpar('bbgray',0.0d0)      !gray value for boundary
	call addpar('bsgray',-1.0d0)     !gray value for plotting maps

! DOCS	END

!c not documented

	call addpar('ifreg',0.d0)		!plot regular grid from fem file
	call addpar('iexreg',1.d0)	!plot half valid boxes in reg grid

	call addfnm('obspnt'," ")	!name of file with obs points
	call addfnm('metpnt'," ")	!name of file with meteo points
	call addfnm('spcvel'," ")	!name of file for velocity points

!c only for compatibility ... are not used anymore

	call addpar('traref',0.d0)	!reference transport for length scale
	call addpar('tramin',0.d0)	!minimum transport to be plotted

!c internally needed

	call addpar('dirn',0.d0)
	call addpar('href',0.d0)
	call addpar('hzmin',0.01d0)
	!call addpar('hzoff',0.05)
	!call addpar('hlvmin',0.25)

	end

!************************************************************************

	subroutine nlsina_color

! $color section

	implicit none

! DOCS	START	S_color
!
! The next parameters deal with the definition of the colors
! to be used in the plot. A color bar is plotted too.

	call sctpar('color')		!sets default section
	call sctfnm('color')		!sets default section

! |icolor|	Flag that determines the type of color table
!		to be used. 0 stands for gray scale, 1 for
!		HSB color table. Other possible values are
!		2 (from white to blue), 3 (from white to red),
!		4 (from blue over white to red) and 
!		5 (from blue over black to red).
!		Values 6 and 7 indicate non-linear HSB color tables.
!		(Default 0)

	call addpar('icolor',0.d0)	!use color (1.) or not (0.)

! |colfil|	A color table can also be read from file. An example
!		of the format can be found in directory |femplot/color|
!		in the file |colormap.dat|. The variable |colfil|
!		indicates the file where the color table is being
!		read from. The default is not to read a color table file.
! |coltab|	If a color table file has been read then the variable
!		|coltab| indicates the name of the color table that
!		is going to be used. The default is to not use any
!		of the color tables if no name is specified.

        call addfnm('colfil',' ')
        call addfnm('coltab',' ')

! |isoval|	Array that defines the values for the isolines
!		and colors that are to be plotted. Values given must
!		be in the unit of the variable that will be plotted,
!		i.e., meters for water levels etc. 
! |color|	Array that gives the color indices for the
!		plotting color to be used. Ranges are from
!		0 to 1. The type of the color depends on the 
!		variable |icolor|. For the gray scale table
!		0 represents black and 1 white. Values in between
!		correspond to tones of gray. For the HSB color table
!		going from 0 to 1 gives the color of the rainbow.
!		There must be one more value in |color| than in |isoval|.
!		The first color in |color| refers to values less
!		than |isoval(1)|, the second color in |color| to
!		values between |isoval(1)| and |isoval(2)|. The last
!		color in |color| refers to values greater than the last
!		value in |isoval|.

	call addpar('isoval',0.d0)	!dummy for array read
	call addpar('color',0.d0)		!dummy for array read

! |x0col, y0col|	Lower left corner of the area where the
!			color bar is plotted.
! |x1col, y1col|	Upper right corner of the area where the
!			color bar is plotted.

	call addpar('x0col',0.d0)		!dimension of color bar
	call addpar('y0col',0.d0)		!dimension of color bar
	call addpar('x1col',0.d0)		!dimension of color bar
	call addpar('y1col',0.d0)		!dimension of color bar

! |cblank|		The color bar is plotted over a white rectangle.
!			Sometimes this blanking is not desireable.
!			If you do not want to have a white box below the
!			legend set |cblank| to 0. (Default 1)

	call addpar('cblank',1.d0)

! |faccol|	Factor for the values that are written to the 
!		color bar legend. This enables you, e.g., to give water level
!		results in mm (|faccol = 1000|). (Default 1)
! |ndccol|	Decimals after the decimal point for the values
!		written to the color bar legend. Use the value |-1|
!		to not write the decimal point. A value of 0 automatically
!		computes the number of decimals needed. (Default 0)
! |legcol|	Text for the description of the color bar. This text
!		is written above the color bar.

	call addpar('faccol',1.d0)	!factor for color bar
	call addpar('ndccol',-1.d0)	!decimals after point
	call addfnm('legcol'," ")	!legend for colorbar

! It is not necessary to give all values for isolines and colors above.
! A faster way is to give only the minimum and maximum values and fix
! the number of isovalues to be used.

! |niso|		Total number of isolines to use. (Default is |nisodf|)
! |nisodf|		Default number of isolines to use. (Default 5)
! |colmin, colmax|	Minimum and maximum color index used. Defaults are
!			0.1 and 0.9 respectively. The value of |colmax| can
!			be smaller than |colmin| which inverts the color
!			index used.
! |valmin, valmax|	Minimum and maximum value for isovalues to be used.
!			There is no default.
! |rfiso|		Defines function to be used to compute intermediate
!			values between |valmin| and |valmax|. If 0 or 1
!			the values are linearlily interpolated. Else they
!			are computed by $y=x^n$ where $n$ is |rfiso|
!			and $x=\frac{v-v_{min}}{v_{max}-v{min}}$. Values
!			for |rfiso| greater than 0 capture higher detail in
!			the lower values, whereas values less than 1 do
!			the opposite.
!			(Default 0)
! |ipllog|		Indicates the usage of a logarithmic color scale.
!			The possible values are 0-3. The value of 0
!			indicates not to use a logarithmic scale.
!			If 1, the values of
!			the scale are 1,10,100,etc., if 2 the values
!			1,2,10,20,100,etc. are used, and for 3 the values
!			are 1,2,5,10,20,50,100,etc. (Default 0)
! |dval|		Difference of values between isolines. If this
!			value is greater then 0 the values for isolines 
!			and the total number of isolines are computed 
!			automatically using also |valmin| and |valmax|. 
!			(Default 0)

        call addpar('niso',0.d0)         !total number of isovalues
        call addpar('nisodf',5.d0)       !default number of isovalues to use
        call addpar('colmin',0.1d0)      !min color [0..1]
        call addpar('colmax',0.9d0)      !max color [0..1]
        call addpar('valmin',0.d0)       !min isovalue
        call addpar('valmax',0.d0)       !max isovalue
	call addpar('rfiso',0.d0)	       !function for intermediate values
	call addpar('ipllog',0.d0)       !logarithmic scale
	call addpar('dval',0.d0)	       !increment for autom. color sep.

! Since there is a great choice of combinations between the parameters,
! we give here the following rules how the values for colors and isolines
! are determined.
!
! If colors are given in array |color|, they are used, else |colmin| and
! |colmax| or their respective defaults are used to determine the color bar.
! If |isoval| is given it is used, else |valmin| and |valmax| are used.
! If |valmin| and |valmax| are not given they are computed every time
! for each plot and the minimum and maximum value in the basin are used.
! In any case, if |isoval| is specified the total number of isovalues
! is known and |niso| is ignored. However, if |isoval| is not given
! then first |dval| is used to decide how many isovalues to plot, and
! if |dval| is 0 then the |niso| and finally |nisodf| is used.
!
! Other parameters that can be changed are the following.

! |nisomx|	Maximum for |niso| allowed. This is especially useful
!		when the value for |niso| is determined automatically.
!		It avoids you to plot 1000 isolines due to wrong settings
!		of |dval|. However, if you want to use 50 isovalues
!		then just set |niso| and |nisomx| to 50. (Default 20)
! |nctick|	Number of values to be written in color bar. If |niso| is high
!		the labels on the color bar become unreadable. Therefore
!		you can use |nctick| to write only some of the values to
!		the color bar. For example, if |valmin| is 0 and |valmax| is
!		5 and you use many isolines, then setting |nctick| to 6 would
!		give you labels at values 0,1,2,3,4,5. If |nctick| is 0
!		then all lables are written. (Default 0)
! |isolin|	Normally the isolines are not drawn on the plot, just
!		the colors are used to show the value in the different
!		parts of the plot. A value different from 0 plots also 
!		the isolines. In this case |isolin| gives the number of
!		isolines to be plotted. A good choice is to make this
!		equal to |nctick|, so that the isolines correspond to the 
!		values	written on the colorbar. For compatibility, a value of
!		1 plots all isolines. (Default 0)
! |isoinp|	Normally inside elements the values are interpolated.
!		Sometimes it is usefull to just plot the value of the
!		node without interpolation inside the element. This can
!		be accomplished by setting |isoinp=0|. (Default 1)

        call addpar('nisomx',20.d0)      !maximum number of isovalues allowed
        call addpar('nctick',0.d0)       !default number of ticks to use
        call addpar('isolin',0.d0)       !plot isolines with color ?
        call addpar('isoinp',1.d0)       !interpolate inside elements

! DOCS	END

	end

!************************************************************************

	subroutine nlsina_arrow

! $arrow section

	implicit none

! DOCS	START	S_arrow
!
! The next parameters deal with the reference arrow that is plotted
! in a legend. The arrow regards the plots where the velocity or
! the transport is plotted.

	call sctpar('arrow')		!sets default section
	call sctfnm('arrow')		!sets default section

! |x0arr, y0arr|	Lower left corner of the area where the
!			reference arrow is plotted.
! |x1arr, y1arr|	Upper right corner of the area where the
!			reference arrow is plotted.

	call addpar('x0arr',0.d0)		!dimension of color bar
	call addpar('y0arr',0.d0)		!dimension of color bar
	call addpar('x1arr',0.d0)		!dimension of color bar
	call addpar('y1arr',0.d0)		!dimension of color bar

! |ablank|		The arrow legend is plotted over a white rectangle.
!			Sometimes this blanking is not desireable.
!			If you do not want to have a white box below the
!			legend set |ablank| to 0. (Default 1)

	call addpar('ablank',1.d0)

! |facvel|		Factor for the value that is written to the 
!			arrow legend for the velocity.
!			This enables you, e.g., to give 
!			velocities in mm/s (|facvel = 1000|). (Default 1)
! |ndcvel|		Decimals after the decimal point for the values
!			written to the arrow legend. 
!			Use the value |-1|
!			to not write the decimal point. (Default 2)
! |legvel|		Text for the description of the arrow legend.
!			This text is written above the arrow.
! |arrvel|		Length of arrow in legend (in velocity
!			units). If not given the arrow length will be computed
!			automatically. (Default 0)
! |sclvel|		Additional factor to be used for the arrow
!			in the legend. When the arrow length will be
!			computed automatically, this parameter gives
!			the possibility to change the length of the
!			reference vector. This is an easy way
!			to scale the velocitiy arrow
!			with an overall factor. Not used if
!			|arrvel| is given. (Default 1)

	call addpar('facvel',1.d0)	!factor for velocity
	call addpar('ndcvel',2.d0)	!decimals after point (velocity)
	call addfnm('legvel'," ")	!legend for velocity
	call addpar('arrvel',0.d0)	!length of arrow
	call addpar('sclvel',1.d0)	!factor for arrow

! DOCS	END

	end

!************************************************************************

	subroutine nlsina_legend

! $legend section

	implicit none

! DOCS	START	S_legend

! In this section annotations in the plots can be given. The
! section consists of a series of lines that must contain the 
! following information:
!
! The first value is a keyword that specifies what has to be plotted.
! Possible values are |text|, |line|, |vect|, |rect|, |circ|  and also
! |wid| and |col|. These correspond to different types of information
! that is inserted into the plot such as text, line, vector, rectangle
! or circle (filled or just outline). Moreover, the color and
! line width of the pen can be controlled by with |wid| and |col|.
!
! In case of |text| the starting position (lower left corner) is given,
! then the point size of the font and the text that is inserted. |line|
! needs the starting and end point of the line. The same with |vect|,
! but in this case also the relative tip size must be given as a final
! parameter. |rect| needs the coordinates of the lower left corner and 
! upper right corner of the rectangle. It also needs the color used for
! the filling of the rectangle (0-1) or the flag -1 which only draws the
! outline of the rectangle without filling it. |circ| needs the center
! point and the radius and a fill color (see rectangle). Finally |wid| needs 
! the relative width of the line and |col| the stroke color used when plotting
! lines.
!
! A small example of an annotation that explains the above parameters
! would be:
!
! \begin{verbatim}
! $legend
! text  30500 11800     15  'Chioggia'   #text, 15pt
! line  30500 11800 35000 15000          #line
! vect  30500 11800 35000 15000 0.1      #arrow, tipsize 0.1
! rect  30500 11800 35000 15000 0.1      #rectangle, fill color 0.1
! rect  30500 11800 35000 15000 -1       #rectangle (outline, no fill)
! circ  30500 11800 5000 -1              #circle (outline, no fill)
! wid   5                                #set line width to 5
! col   0.5                              #set color to 0.5
! $end
! \end{verbatim}
!
! There is also an old way to specify the legend that does not use
! keywords. However, this way is deprecated and unsupported and is therefore
! not described anymore in this manual.

! DOCS	END

	end

!************************************************************************

	subroutine nlsina_legvar

! variable legend section

	implicit none

! DOCS	START	S_legvar
!
! In this section variable fields like the date and wind vectors
! may be inserted into the plot.

	call sctpar('legvar')		!sets default section
	call sctfnm('legvar')		!sets default section

! A time and date can be assigned to the simulation results. These values
! refer to the time 0 of the FEM model. The format for the date is
! YYYYMMDD and for the time HHMMSS. 
!c Please note that the date should not be
!c given as YYYYMMDD because due to precision problems this will not work. 
! You can also give a time zone if your time is not referring to 
! GMT but to another time zone such as MET. Please note that you have to give
! this information only if the simulation does not contain it already.
! Normally, this information is already assigned during the simulation runs.

! |date|		The double precision date corresponding to time 0. (Default 0)
! |time|		The double precision time corresponding to time 0. (Default 0)
! |tz|			The time zone you are in. This is 0 for GMT, 1 for MET
!			and 2 for MEST (MET summer time). (Default 0)
! |tzshow|		The time zone you want to show your results. If your
!			time zone is GMT (0) and you want to show the results
!			referred to MET (+1) set this to +1. Please note that
!			you have to set this variable only if you want to
!			show results in a different time zone than the one
!			given in |tz|. (Default 0)

	call addpar('date',-1.d0)
	call addpar('time',0.d0)
	call addpar('tz',0.d0)
	call addpar('tzshow',0.d0)

! The information of date and time may be written to the plot. This
! is done with the following parameters.

! |xdate, ydate|	Starting point for the date text (lower left corner).
! |sdate|		Point size of the text. (Default 18)
! |idate|		Output mode. If 0 no date is written to the
!			plot, else the date and time is written. (Default 0)

	call addpar('xdate',0.d0)
	call addpar('ydate',0.d0)
	call addpar('sdate',18.d0)         !size
	call addpar('idate',0.d0)         !mode

! Wind data can be used to insert a wind vector into the figure.
! This is useful because in the case of variable wind 
! the direction and speed of the wind that was blowing
! in the moment of the plot is shown.
!
! Since only one wind vector can be plotted, the wind data must consist
! of one value for each time. The same ASCII file that is used
! in the STR file can be used.

! |xwind, ywind|	Starting point where the wind arrow is plotted.
! |iwtype|		Type of wind data. The same as the one in the
!			STR file. If this parameter is 0 then no
!			wind vector is plotted. (Default 0)
! |lwwind|		Line width of the wind vector. (Default 0.1)
! |scwind|		Scaling parameter of the wind vector. This depends
!			on the size of your plot. If your wind is 10 m/s
!			and you want the vector to strech over a distance
!			of 5 km on the plot then you have to choose
!			the value of 500 (10*500=5000) for |scwind|.
!			(Default 1)
! |wfile|		Name of the file containing the wind data. This 
!			may be the same file than the one used in the
!			STR file to run the program.

	call addpar('xwind',0.d0)
	call addpar('ywind',0.d0)
	call addpar('iwtype',1.d0)         !mode
	call addpar('lwwind',0.1d0)        !line width
	call addpar('scwind',1.d0)         !scale
	call addfnm('wfile',' ')	 !name of wind file

! The wind vector is also given a text legend with the speed of the
! wind written out. The next parameters decide where and how this information
! is put.

! |xtwind, ytwind|	Starting point for the legend text (lower left corner).
! |stwind|		Point size of the text. (Default 18)
! |wtext|		Text used for the legend (Default 'Wind speed')
! |wunit|		Unit for the wind speed (Default 'm/s')

	call addpar('xtwind',0.d0)
	call addpar('ytwind',0.d0)
	call addpar('stwind',18.d0)         !size
	call addfnm('wtext','Wind speed') !legend for wind
	call addfnm('wunit','m/s')	  !unit for wind

! DOCS	END

	end

!*******************************************************************

	subroutine nlsina_sect

! vertical section

	implicit none

! DOCS	START	S_sect
!
! In this section a vertical section plot can be set up. This will
! allow to plot 3D variables not horizontally for each layer, but
! along a given section in the vertical

	call sctpar('sect')		!sets default section
	call sctfnm('sect')		!sets default section

! |vsect|		Name of the file containing the node list defining
!			the section. The nodes must be adjacent in the
!			numerical grid used. There should be one node
!			on each line. The program |make_line_nodes.sh|
!			can be used to produce such a node list from a
!			given line created with |grid|. If this file
!			is not given, no vertical section will be plotted.

	call addfnm('vsect','')		!name of line defining vertical section

! The next parameters decide about how the plot is scaled in the vertical
! and if a grid overlay is created.
!
! |ivert|		This parameter decides what is plotted on the
!			vertical axis. If 0, depth is plotted on the axis.
!			If 1, vertical layers are plotted. If 2, depth is
!			plotted, but it is scaled with a logarithm, giving
!			more importance to the layers close to the surface.
!			(Default 0)
! |ivgrid|		If 1 a grid is plotted as overlay over the
!			vertical section plot. This is mainly needed
!			for a better orientation where the nodes and 
!			depth values or layers are situated. (Default 0)
! |hvmax, lvmax|	Normally the whole vertical extension is plotted
!			for the section. Sometimes only the upper part
!			of the water column may be needed. With these
!			two parameters the maximum depth (|hvmax|) or
!			the maximum layer (|lvmax|) to be plotted can
!			be imposed. The values may be also higher than
!			the maximum depth/layer. In this case extra space
!			is added below the plotting area. Only one
!			of the two parameters can be set.
! |bsmt|		Set to 1 to have smooth bottom profile even in 
!			case of z-layer structure. (Default 0)

	call addpar('ivert',0.d0)		!0: depth  1: layers  2: log depth
	call addpar('ivgrid',0.d0)	!plot grid layout over plot
	call addpar('hvmax',0.d0)		!maximum depth to be plotted
	call addpar('lvmax',0.d0)		!maximum layer to be plotted
	call addpar('bsmt',0.d0)		!1: smooth bottom

! The vertical section plot also creates a color bar. This color bar is
! normally put on the right side of the plot. If there is space inside the
! plot it might be placed on top of the section plot. In this case the
! following parameters can be used. Please note that in this case it is
! the relative position (from 0 to 1) that has to be specified.
!
! |x0sect, y0sect|	Lower left corner of the area where the
!			color bar is plotted.
! |x1sect, y1sect|	Upper right corner of the area where the
!			color bar is plotted.

	call addpar('x0sect',0.d0)
	call addpar('y0sect',0.d0)
	call addpar('x1sect',0.d0)
	call addpar('y1sect',0.d0)

! The following parameters set titles for the plot, the axis, and for
! an extra description of the starting and end point of the plot.
! They have some (intelligent) default values.
!
! |vtitle|		Title of the plot.
! |xtitle|		Title for the x-axis.
! |ytitle|		Title for the y-axis.
! |ltitle|		Title for start point of the line. (No default)
! |rtitle|		Title for end point of the line. (No default)

	call addfnm('vtitle','Section Plot')	!title for plot
	call addfnm('xtitle','Distance [m]')	!title for x-axis
	call addfnm('ytitle','Depth [m]')	!title for y-axis
	call addfnm('ltitle','')		!title for left point
	call addfnm('rtitle','')		!title for right point

! When plotting velocities you can decide if using for the color the normal 
! velocitiy accross the section or the tangent velocity.
! In any case, in order to visualize also the veloctiy
! tangent to the section arrwos are used. The next parameters deal with
! the scaling of these arrows.

! |vmode|	When plotting velocities as default the normal velocitiy 
!		accross the section is used for the color plot (|vmode|=0).
!		If you want to use the tangential velocity, please set
!		|vmode|=1. (Default 0)
! |avscal|	This parameter defines the horizontal scale for the 
!		velocity vector. It defines the length scale 
!		in units of x-coordinates of the vertical section so
!		that a velocity of 1 m/s fits comfortably
!		into the plot. If 0 the scale is computed automatically.
!		Please note that in this case the velocities will be
!		scaled differently for every plot. Setting |avscal| 
!		therefore guarantees that the velocity arrows will
!		be comparable between plots. If |avscal| is negative,
!		instead of using x-coordinates units, the legth scale is
!		in [cm]. Therefore, a value of -1 will represent a
!		velocity of 1 m/s with 1 cm on the plot. This unit for
!		the length scale is sometimes easier to understand. 
!		(Default 0)
! |rvscal|	Extra factor that multiplies the scale factor. If your
!		automatic scale gives you vectors which are too small, you
!		can use |rvscal| to increase them. (Default 1)
! |rwscal|	Extra factor for the vertical scale. Normaly the
!		vertical scale is computed automatically, If you dont
!		like the size of the vertical vectors you can controll
!		it with this parameter. A value of 2 will give you
!		a vertical vector twice as big a the default. (Default 1)
! |dxmin|	Sometimes there are two many arrows plotted horizonally.
!		The value of |dxmin| gives the minimum distance arrows have
!		to be apart to be plotted. A value of 0 plots all arrows.
!		(Default 0)
! |svtip|	The (relative) tip size of the arrow. It specifies how
!		big the arrow will be drawn. A value of 0 only draws the
!		arrow line without tip, and a negative value inhibits
!		drawing arrows at all. (Default 0.3)
! |rxscal,ryscal|	In case arrows are plotted, also a reference
!			vector is plotted. The size of this reference
!			vector s computed automatically, but can be 
!			controlled additionally by the parameters
!			|rxscal,ryscal|, which are in relative units
!			with respect to the reference box plotted.
!			(Default 0.6)

	call addpar('vmode',0.d0)
	call addpar('avscal',0.d0)
	call addpar('rvscal',1.d0)
	call addpar('rwscal',1.d0)
	call addpar('dxmin',0.d0)
	call addpar('svtip',0.3d0)
	call addpar('rxscal',0.6d0)
	call addpar('ryscal',0.6d0)

!c still to do: plot vector legend

! DOCS	END

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine fnminh

! initializes default names

	implicit none

! para section

	call sctfnm('name')		!sets default section

! DOCS	START	S_name
!
! DOCS	COMPULS		Directory specification
!
! This parameters define directories for various input
! and output files.
!
! |basdir|	Directory where basin file BAS resides. (Default .)
! |datdir|	Directory where output files are written. (Default .)
! |tmpdir|	Directory for temporary files. (Default .)
! |defdir|	Default directory for other files. (Default .)

	call addfnm('basdir',' ')
	call addfnm('datdir',' ')
	call addfnm('tmpdir',' ')
	call addfnm('defdir',' ')

! DOCS	END

! DOCS	START	S_name_h
!
! DOCS	FILENAME	File names
!
! The following strings enable the specification of files
! that account for initial conditions or forcing.
!
! |zinit|	Name of file containing initial conditions for water level
! |uvinit|	Name of file containing initial conditions for velocity
! |wind|	File with wind data. The file may be either
!		formatted or unformatted. For the format of the unformatted
!		file please see the section where the WIN
!		file is discussed.
!		The format of formatted ASCII file
!		is in standard time-series format, with the
!		first column containing the time in seconds and
!		the next two columns containing the wind data.
!		The meaning of the two values depend on the
!		value of the parameter |iwtype| in the |para| section.
! |rain|	File with rain data. This file is a standard time series
!		with the time in seconds and the rain values
!		in mm/day. The values may include also evaporation. Therefore,
!		also negative values (for evaporation) are permitted.
! |qflux|	File with heat flux data. This file must be in
!		a special format to account for the various parameters
!		that are needed by the heat flux module to run. Please
!		refer to the information on the file |qflux|.
! |surfvel|	File with surface velocities from observation. These
!		data can be used for assimilation into the model.
! |ice|		File with ice cover. The values range from 0 (no ice cover)
!		to 1 (complete ice cover).
! |restrt|	Name of the file if a restart is to be performed. The
!		file has to be produced by a previous run
!		with the parameter |idtrst| different
!		from 0. The data record to be used in the file for the
!		restart must be given by time |itrst|.
! |gotmpa|	Name of file containing the parameters for the
!		GOTM turbulence model (iturb = 1).
! |saltin|	Name of file containing initial conditions for salinity
! |tempin|	Name of file containing initial conditions for temperature
! |conzin|	Name of file containing initial conditions for concentration
! |bfmini|	Name of file containing initial conditions for bfm
! |offlin|	Name of the file if a offline is to be performed. The
!		file has to be produced by a previous run
!		with the parameter |idtoff| greater than 0.


	call addfnm('zinit',' ')
        call addfnm('uvinit',' ')

	call addfnm('wind',' ')
        call addfnm('rain',' ')
        call addfnm('qflux',' ')
        call addfnm('surfvel',' ')
        call addfnm('ice',' ')

	call addfnm('restrt',' ')
	call addfnm('gotmpa',' ')

        call addfnm('saltin',' ')
        call addfnm('tempin',' ')
        call addfnm('conzin',' ')
        call addfnm('bfmini',' ')
	call addfnm('offlin',' ')

! DOCS	END

!c non-documented -> try first	HACK	-> initial conditions

        call addfnm('bioin',' ')
        call addfnm('biosin',' ')
        call addfnm('toxi',' ')
        call addfnm('mercin',' ')	!mercury

!c ACQUBC

	call fnm_aquabc_init

!c DOCS	DELWAQ		delwaq model

	call addfnm('voldwq',' ')
	call addfnm('flowdwq',' ')
	call addfnm('areadwq',' ')
	call addfnm('femdwq',' ')
	call addfnm('pntfem',' ')
	call addfnm('ftodwq',' ')

!c DOCS	INTERNAL	internal administration - blank section

	call sctfnm(' ')		!sets default section

	call addfnm('title',' ')

	if( .not. haspar('basnam') ) call addfnm('basnam',' ')
	if( .not. haspar('runnam') ) call addfnm('runnam',' ')
	if( .not. haspar('apnnam') ) call addfnm('apnnam',' ')

	call addfnm('apnfil','apnstd.str')
	!call addfnm('colfil','apncol.str')
	call addfnm('memfil','.memory')
!dos#	call putfnm('memfil','_memory')
!lahey#	call putfnm('memfil','_memory')

	call addfnm('pltfil','plot')

	end

!************************************************************************

        subroutine fnm_aquabc_init

        implicit none

!c for model aquabc (curonian)

        call addfnm('biocon',' ')
        call addfnm('bioscon',' ')
        call addfnm('biolight',' ')
        call addfnm('bioaow',' ')  !ascii output for WC
        call addfnm('bioaos',' ')  !ascii output for BS
        call addfnm('bioph',' ')
        call addfnm('biotemp',' ')
        call addfnm('bioload',' ')

        end

!************************************************************************

!--------------------------------------------------------------------------------
        end module def_para
!--------------------------------------------------------------------------------
