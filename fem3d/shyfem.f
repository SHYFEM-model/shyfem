
!--------------------------------------------------------------------------
!
!    Copyright (C) 1995,1997-2020  Georg Umgiesser
!    Copyright (C) 2004,2008  Andrea Cucco
!    Copyright (C) 2006,2008,2011-2012,2014-2015,2014-2015  Christian Ferrarin
!    Copyright (C) 2019  Christian Ferrarin
!    Copyright (C) 2007,2013  Debora Bellafiore
!    Copyright (C) 2012  Aaron Roland
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

c finite element model shyfem (version 3D)
c
c original version from march 1991
c
c revision log :
c
c 30.08.1995	ggu	$$AUST - austausch coefficient introduced
c 11.10.1995	ggu	$$BCLBND - boundary condition for baroclinic runs
c 04.08.1997	ggu	$$ZEONV - new arrays for water level elementwise
c 19.03.1998	ggu	$$IPCCV close data items commented or deleted
c 03.04.1998	ggu	$$DESCRP BUG overwriting descrp with boundary vals.
c 30.04.1998	ggu	finally eliminated /semimp/, /trock/, /ffloat/
c 05.05.1998	ggu	izeit eliminated
c 28.05.1998	ggu	call to sp131g changed
c 14.08.1998	ggu	call to biocos (biological model) introduced
c 20.08.1998	ggu	iextpo finally eliminated, momentum arrays added
c 20.08.1998	ggu	spb11 -> sp111
c 21.08.1998	ggu	xv eliminated
c 03.09.1998	ggu	biological reactor integratated
c 06.11.1998	ggu	hv renamed into hkv
c 08.04.1999	ggu	equilibrium tide (tidal potential) introduced
c 19.04.1999	ggu	subroutine custom called
c 22.10.1999	ggu	igrv, ngrv eliminated (subst by ilinkv, lenkv)
c 20.01.2000	ggu	common block /dimdim/ eliminated
c 07.03.2000	ggu	arrays for vertical turbulent diffusion coefficient
c 20.06.2000	ggu	alv -> visv    slv -> difv
c 19.11.2001	ggu	h1hv eliminated, ulv, vlv eliminated
c 03.12.2001	ggu	new arrays [cst]difhv
c 21.08.2002	ggu	new arrays /kvolc/, /ivol/ copied from hp
c 31.07.2003	ggu	array winv eliminated
c 10.08.2003	ggu	big restructuring
c 13.08.2003	ggu	some more restructuring
c 14.08.2003	ggu	even more restructuring and cleaning up
c 04.12.2003	ggu	integration of wave and sediment module
c 06.03.2004	aac	lagrangian trajectories computation module
c 03.09.2004	ggu	restart now in ht
c 15.10.2004	ggu	gotm substituted by general turbulence closure
c 02.12.2004	ggu	variable time step implemented
c 17.01.2005	ggu	new routine diff_h_set, new difhv, ausv deleted
c 24.02.2005	ggu	new routine smagorinski and common smagv
c 03.03.2005	ggu	smagorinski uses difhv, inserted in subdif.f
c 03.03.2005	ggu	new 3d boundary arrays for C/T/S
c 12.08.2005	ggu	minimum level index, TVD arrays
c 07.11.2005	ggu	new array sed2dn introduced (file name for sediments)
c 16.02.2006	ggu	new routine atoxi3d (TOXI)
c 23.03.2006	ggu	changed time step to real
c 18.10.2006	ccf	radiation stress and waves included (with pipe)
c 10.11.2006	ggu	initialize depth values after restart
c 16.11.2006	ggu	turbulence values included
c 02.04.2007	ggu	changes in algorithm (ASYM)
c 31.05.2007	dbf	new arrays bpresxv, bclevvar (debora)
c 26.09.2007	ggu	deleted arrays rcv,rtv,rsv
c 27.09.2007	ggu	deleted call to tstvol,tstvol1
c 20.03.2008	aac	new call for ERSEM ecological model (BFM MODULE)
c 07.04.2008	aac	new array bfm*bc introduced (file name for ersem)
c 10.04.2008	ggu&ccf	upro, waveov, stokes, z0bk
c 16.04.2008	ggu	evaporation mass flux (evapv)
c 22.04.2008	ggu	gradx/yv non global due to parallelization
c 23.04.2008	ggu	deleted r3{c|t|s}v
c 28.04.2008	ggu	new routine init_stability(), new conzm3sh()
c 29.04.2008	ggu&aac	new BMF ERSEM model
c 12.11.2008	ggu	new sigma level initialization
c 10.12.2008	ggu	new array rfricv for bottom friction
c 18.12.2008	ggu	new routine debug_output(), might be removed later
c 04.03.2009	ggu	matrix amat into solver routines
c 11.03.2009	ggu	new arrays for meteo data
c 06.04.2009	ggu	renamed nlidim to nlkdim
c 30.11.2009	ggu	write output file for successfull completion
c 19.02.2010	ggu	init_stability() changed to reset_stability()
c 22.02.2010	ggu	new arrays wxv, wyv
c 26.02.2010	ggu	new arrays sauxe1/2
c 29.04.2010	ggu	write volumes (wrfvla)
c 26.01.2011	ggu	new arrays for observations and nudging
c 16.02.2011	ggu	new iarnv, call to aquabc
c 17.02.2011	ccf	new radiation stress in 3D
c 23.03.2011	ggu	new call to adjust_spherical()
c 31.03.2011	ggu	write finite volumes at initial time step
c 20.05.2011	ggu	iwetv introduced, wet and dry from main
c 25.10.2011	ggu	hlhv eliminated
c 18.11.2011	ggu	new routine handle_projection
c 24.01.2012	ggu	new call to setup_parallel()
c 23.02.2012	ggu&ccf	meteo arrays adjusted (3*nkn)
c 09.03.2012	ggu	call to residence time added
c 21.06.2012	ggu&aar	fluid mud variables integrated
c 05.08.2012	ggu	bug because lam2dn and dmfd2n not defined
c 10.05.2013	dbf	initialization for non hydrostatic routines
c 13.06.2013	ggu	set/copydepth simplified, offline version
c 05.09.2013	ggu	changed order of adjust depth and barene structures
c 29.10.2013	ggu	nudging implemented
c 25.03.2014	ggu	new offline
c 25.06.2014	ggu	new arrays hkv_min, hkv_max
c 05.12.2014	ccf	new interface for waves
c 30.07.2015	ggu	routine renamed from ht to shyfem
c 18.09.2015	ggu	new routine scalar, call to hydro()
c 23.09.2015	ggu	changed VERS_7_2_4
c 29.09.2015	ccf	inverted set_spherical() and handle_projection()
c 30.09.2015	ggu	changed VERS_7_2_6
c 10.10.2015	ggu	fluid mud routines handled differently
c 12.10.2015	ggu	changed VERS_7_3_3
c 13.10.2015	ggu	changed VERS_7_3_5
c 22.10.2015	ggu	changed VERS_7_3_7
c 23.10.2015	ggu	changed VERS_7_3_9
c 05.11.2015	ggu	changed VERS_7_3_12
c 09.11.2015	ggu	changed VERS_7_3_13
c 20.11.2015	ggu	changed VERS_7_3_15
c 16.12.2015	ggu	changed VERS_7_3_16
c 19.02.2016	ggu	changed VERS_7_5_2
c 22.02.2016	ggu	changed VERS_7_5_4
c 07.06.2016	ggu	changed VERS_7_5_12
c 14.06.2016	ggu	changed VERS_7_5_14
c 13.04.2017	ggu	changed VERS_7_5_25
c 09.05.2017	ggu	changed VERS_7_5_26
c 16.05.2017	ggu	changed VERS_7_5_27
c 05.10.2017	ggu	command line options introduced, subs rearranged
c 14.11.2017	ggu	changed VERS_7_5_36
c 17.11.2017	ggu	changed VERS_7_5_37
c 05.12.2017	ggu	changed VERS_7_5_39
c 07.12.2017	ggu	changed VERS_7_5_40
c 03.04.2018	ggu	changed VERS_7_5_43
c 19.04.2018	ggu	changed VERS_7_5_45
c 26.04.2018	ggu	changed VERS_7_5_46
c 11.05.2018	ggu	changed VERS_7_5_47
c 06.07.2018	ggu	changed VERS_7_5_48
c 25.10.2018	ggu	changed VERS_7_5_51
c 18.12.2018	ggu	changed VERS_7_5_52
c 12.02.2019	ccf	bottom shear stress in substress.f
c 16.02.2019	ggu	changed VERS_7_5_60
c 12.03.2019	ccf	include new computation of tide potential/analysis
c 04.07.2019	ggu	new ww3 routines introduced
c 15.09.2019	ggu	subroutine to test only forcing
c 02.10.2019	ggu	delete include files
c 17.10.2019	ggu	no call to bfm_write, is done inside subroutine
c 06.11.2019	ggu	femelab eliminated
c 03.04.2020	ggu	write real start and end time of simulation
c 09.04.2020    ggu     run bfm through bfm_run()
c 21.05.2020    ggu     better handle copyright notice
c 04.06.2020    ggu     debug_output() substituted with shympi_debug_output()
c 30.03.2021    ggu     more on debug, call sp111(2) outside time loop
c 01.04.2021    ggu     turbulence cleaned
c
c*****************************************************************
c
c notes :
c
c MPI calls
c
c call shympi_init(bmpirun)
c   call shympi_alloc_global(nkn,nel,nen3v,ipv,ipev)
c call shympi_setup			!sets up partitioning of basin
c
c----------------------------------------------------------------

	program shyfem

	use mod_bound_geom
	use mod_geom
	use mod_meteo
	use mod_waves
	use mod_sinking
	use mod_nudging
	use mod_internal
	use mod_geom_dynamic
	use mod_depth
	use mod_bnd_aux
	use mod_diff_aux
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_roughness
	use mod_diff_visc_fric
	use mod_hydro_baro
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
	use intp_fem_file
	use tide
	use projection
	use coordinates
	use mod_subset
	use mod_bfm
        use mod_nohyd !DWNH
!$	use omp_lib	!ERIC
	use shympi
	use mod_test_zeta

c----------------------------------------------------------------

	implicit none

c local variables

	logical bdebout,bdebug,bmpirun
	logical bfirst
	integer k,ic,n
	integer iwhat
	integer date,time
	integer nthreads
	integer*8 count1,count2,count3,count_rate,count_max
	real time1,time2,time3
	real dt
	double precision timer
	double precision mpi_t_start,mpi_t_end,parallel_start
	double precision mpi_t_solve
	double precision dtime,dtanf,dtend
        double precision atime_start,atime_end
        character*20 aline_start,aline_end
	character*80 strfile
	character*80 mpi_code,omp_code

	real getpar

	call cpu_time(time1)
	call system_clock(count1, count_rate, count_max)
!$      timer = omp_get_wtime() 

        call get_real_time(atime_start,aline_start)

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%% code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c-----------------------------------------------------------
c copyright and command line options
c-----------------------------------------------------------

        call shyfem_init(strfile,bdebug,bdebout,bmpirun)

c-----------------------------------------------------------
c write  real start time
c-----------------------------------------------------------

        write(6,*) 'simulation start:   ',aline_start 

c-----------------------------------------------------------
c read STR file
c-----------------------------------------------------------

	call cstinit
	call cstfile(strfile)			!read STR and basin

	call shympi_init(bmpirun)
	call setup_omp_parallel

	call cpu_time(time3)
	call system_clock(count3, count_rate, count_max)
        mpi_t_start = shympi_wtime()
	call shympi_setup			!sets up partitioning of basin
        parallel_start = shympi_wtime()

	call test_shympi_arrays

	call allocate_2d_arrays

c-----------------------------------------------------------
c check parameters read and set time and Coriolis
c-----------------------------------------------------------

	call cstcheck

c	call pritst(1)

c-----------------------------------------------------------
c initialize triangles
c-----------------------------------------------------------

	!call levels_init_2d(nkn,nel)	!maybe not needed
	call set_spherical
	call set_ev
	call adjust_spherical
	call print_spherical
	call handle_projection
	call set_geom
	call shympi_barrier
	call domain_clusterization
	call shympi_barrier

c-----------------------------------------------------------
c inititialize time independent vertical arrays
c-----------------------------------------------------------

	call adjust_depth	!adjusts hm3v
	call init_vertical	!makes nlv,hlv,hldv,ilhv,ilhkv, adjusts hm3v

c-----------------------------------------------------------
c allocates arrays
c-----------------------------------------------------------

	call allocate_3d_arrays
	call set_depth		!makes hev,hkv and exchanges

	call check_point('checking ht 1')

	call poisson_init

c-----------------------------------------------------------
c initialize barene data structures
c-----------------------------------------------------------

	if( .not. bmpi ) call setweg(-1,n)	!shympi - FIXME
	call setnod
	call update_geom	!update ieltv - needs inodv
	call populate_strings	!populate strings here

c-----------------------------------------------------------
c initialize boundary conditions
c-----------------------------------------------------------

	call check_point('checking ht 2')

	call get_date_time(date,time)
	call iff_init_global(nkn,nel,nlv,ilhkv,ilhv
     +				,hkv_max,hev,hlv,date,time)

	call sp111(1)           !here zenv, utlnv, vtlnv are initialized

c-----------------------------------------------------------
c initialize depth arrays and barene data structure
c-----------------------------------------------------------

	call setweg(0,n)
	!if( bmpi .and. n > 0 ) goto 95
	call setznv		! -> change znv since zenv has changed

        call rst_perform_restart        !restart
	call compute_velocities
	call copy_uvz

	!call init_vertical	!do again after restart

	call setnod
	call set_area

	call make_new_depth
	call copy_depth
	call make_new_depth
	!call check_max_depth

	call setup_time		!in case start time has changed with rst
	call get_act_dtime(dtime)
	call get_first_dtime(dtanf)
	call get_last_dtime(dtend)

c-----------------------------------------------------------
c initialize open boundary routines
c-----------------------------------------------------------

	call check_point('checking ht 3')

	call bndo_init
	if( bdebug ) then
	  call bndo_info_file('bndo_info.dat')
	end if

c-----------------------------------------------------------
c initialize transports and velocities
c-----------------------------------------------------------

	call init_uv            !set vel, w, pr, ... from transports
	call copy_uvz		!copy new to old
	call barocl(0)
	call wrfvla		!write finite volume
	call nonhydro_init
	call init_wave		!waves
	call ww3_init
	call initsed		!sediments
        call init_bstress	!bottom shear stress

c-----------------------------------------------------------
c initialize modules
c-----------------------------------------------------------

	call init_others
	call init_chezy
	call init_nodal_area_code	!area codes on nodes
        call diffweight
        call set_diffusivity
	call tidefini
	call sp136(ic)
        call shdist(rdistv)
	call tracer_init
        call qhdist(qdistv) !DWNH
	call bfm_init
	call renewal_time
        call lagrange
	call tidepar_init
	call submud_init
	call handle_gotm_init

	call cstsetup

c-----------------------------------------------------------
c write input values to log file and perform check
c-----------------------------------------------------------

	call check_point('checking ht 4')

	call check_fem
	call check_values
	call prilog

	call bclfix_ini

	call system_initialize		!matrix inversion routines

	call poisson_compute

	call offline(2)
	call offline(1)

	call init_nudging

	call do_init

	call sp111(2)           	!initialize BC and read first data

	!call custom(it)		!call for initialization

	!write(6,*) 'starting time loop'
        call shympi_comment('starting time loop...')

	call print_time

	call check_parameter_values('before main')

	if( bdebout ) call handle_debug_output(dtime)

        !call test_forcing(dtime,dtend)

	call test_zeta_init

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%% time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	bfirst = .true.

	do while( dtime .lt. dtend )

           if(bmpi_debug) call shympi_check_all

	   call check_crc
	   call set_dry

           call set_timestep		!sets dt and t_act
           call get_timestep(dt)
	   call get_act_dtime(dtime)

	   call do_befor

	   call copy_uvz		!copies new to old time level
	   call nonhydro_copy   	!copies non hydrostatic pressure terms
	   call copy_depth		!copies layer depth to old

	   call offline(2)		!read from offline file
	   call sp111(2)		!boundary conditions
           call read_wwm		!wwm wave model
	   
           if(bmpi_debug) call shympi_check_all

	   call hydro			!hydro

	   call run_scalar

           call turb_closure

           call parwaves                !parametric wave model
           call sedi                    !sediment transport
	   call submud                  !fluid mud (ARON)
	   call simple_sedi		!simplified sediment module
           call bstress			!bottom shear stess

	   call renewal_time
	   call ecological_module	!ecological model
           call atoxi3d			!toxi
           call mercury_module

           call lagrange

	   call offline(1)		!write to offline file

	   call do_after

	   call tracer_write

           if( bfirst ) call print_file_usage

	   call print_time			!output to terminal

	   call total_energy
	   call check_all
	   !call check_layer_depth
	   !call check_special

           call write_wwm
	   call ww3_loop

	   call mpi_debug(dtime)
	  if( bdebout ) call handle_debug_output(dtime)
	   bfirst = .false.

	   call test_zeta_write

	end do

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%% end of time loop %%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call system_finalize		!matrix inversion routines

	call check_parameter_values('after main')

	call print_end_time

	if( shympi_is_master() ) then

!$OMP PARALLEL
!$OMP MASTER
	nthreads = 1
!$	nthreads = omp_get_num_threads()
	call openmp_parallel_code(omp_code)
	print *,"TYPE OF OMP CODE            = ",trim(omp_code)
        print *,"NUMBER OF OMP THREADS USED  = ",nthreads
!$OMP END MASTER
!$OMP END PARALLEL

!$      timer = omp_get_wtime() - timer
!$      print *,"TIME TO SOLUTION (OMP)  = ",timer

	call system_clock(count2, count_rate, count_max)
	count2 = count2-count1
	timer = count2
	timer = timer / count_rate
	print *,"TIME TO SOLUTION (WALL) = ",timer,my_id

	call cpu_time(time2)
	print *,"TIME TO SOLUTION (CPU)  = ",time2-time1,my_id
	print *,"TIME TO SOLUTION PARALLEL REGION (CPU) = "
     +				,time2-time3,my_id

	call shympi_parallel_code(mpi_code)
	print *,"TYPE OF MPI CODE            = ",trim(mpi_code)
	print *,"NUMBER OF MPI THREADS USED  = ",n_threads

        mpi_t_end = shympi_wtime()
        write(6,*)'MPI_TIME =',mpi_t_end-mpi_t_start,my_id
        write(6,*)'Parallel_TIME =',mpi_t_end-parallel_start,my_id
	call shympi_time_get(1,mpi_t_solve)
        write(6,*)'MPI_SOLVE_TIME =',mpi_t_solve,my_id

        call get_real_time(atime_end,aline_end)

        write(6,*) 'simulation start:   ',aline_start 
        write(6,*) 'simulation end:     ',aline_end 
        write(6,*) 'simulation runtime: ',atime_end-atime_start

	end if

	call test_zeta_close

        !call ht_finished

	!call pripar(15)
	!call prifnm(15)

	!call shympi_finalize
	call shympi_exit(99)
	call exit(99)

	stop
   95	continue
        write(6,*) 'no wetting and drying for mpi yet...',n
        stop 'error stop shyfem: not yet ready'
        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine ht_finished

        implicit none

        open(1,file='SHYFEM_FINISHED',status='unknown',form='formatted')
        write(1,*) 'SHYFEM has finished successfully the simulation'
        close(1)

	end

c*****************************************************************

        subroutine shyfem_init(strfile,bdebug,bdebout,bmpirun)

        use clo
        use shympi

        implicit none

        character*(*) strfile
        logical bdebug,bdebout,bmpirun
        logical bquiet,bsilent

        character*80 version

	call get_shyfem_version_and_commit(version)
        call clo_init('shyfem','str-file',trim(version))

        call clo_add_info('runs the 3D shyfem routine')

	call clo_add_sep('general options:')
        call clo_add_option('quiet',.false.,'do not be verbose')
        call clo_add_option('silent',.false.,'be silent')

	call clo_add_sep('mpi options:')
        call clo_add_option('mpi',.false.
     +			,'runs in MPI mode (experimental)')

	call clo_add_sep('debug options:')
        call clo_add_option('debug',.false.,'enable debugging')
        call clo_add_option('debout',.false.
     +			,'writes debugging information to file')

        call clo_parse_options

        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)

        call clo_get_option('mpi',bmpirun)

        call clo_get_option('debug',bdebug)
        call clo_get_option('debout',bdebout)

        if( bsilent ) bquiet = .true.

	if( shympi_is_master() ) then
         call shyfem_set_short_copyright(bquiet)
         if( .not. bsilent ) then
	  call shyfem_copyright('shyfem - 3D hydrodynamic SHYFEM routine')
         end if
	end if

        call clo_check_files(1)
        call clo_get_file(1,strfile)

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine allocate_2d_arrays

	use mod_bndo
	use mod_tvd
	use mod_bnd
	use mod_bound_geom
	use mod_geom
	use mod_meteo
	use mod_geom_dynamic
	use mod_bnd_aux
	use mod_diff_aux
	use mod_roughness
	use mod_hydro_baro
	use mod_depth
	use evgeom
	use tide
	use coordinates
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	call tide_init(nkn)
	call coordinates_init(nkn)

	call mod_hydro_baro_init(nel)
	call mod_roughness_init(nkn)
	call mod_diff_aux_init(nel)

	call mod_bnd_aux_init(nkn,nel)
	call mod_geom_dynamic_init(nkn,nel)

	call mod_meteo_init(nkn)
	call mod_geom_init(nkn,nel,ngr)
	call mod_bndo_init(ngr,nrb)

	call mod_depth_init(nkn,nel)

	call ev_init(nel)

	call mod_tvd_init(nel)

	write(6,*) '2D arrays allocated: ',nkn,nel,ngr

	end

c*****************************************************************

	subroutine allocate_3d_arrays

	use mod_conz
	use mod_waves
	use mod_sediment
	use mod_bstress
	use mod_keps
	use mod_sinking
	!use mod_fluidmud
	use mod_bclfix
	use mod_nohyd
	use mod_internal
	use mod_layer_thickness
	use mod_gotm_aux
	use mod_bound_dynamic
	use mod_nudging
	use mod_area
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer nlvddi

	nlvddi = nlvdi

	call mod_hydro_init(nkn,nel,nlvddi)
	call mod_hydro_vel_init(nkn,nel,nlvddi)
	call mod_hydro_print_init(nkn,nlvddi)

	call mod_diff_visc_fric_init(nkn,nel,nlvddi)

	call mod_ts_init(nkn,nlvddi)

	call mod_area_init(nkn,nlvddi)
	call mod_bound_dynamic_init(nkn,nlvddi)

	call mod_gotm_aux_init(nkn,nlvddi)
	!call mod_fluidmud_init(nkn,nlvddi)
	call mod_keps_init(nkn,nlvddi)

	call mod_layer_thickness_init(nkn,nel,nlvddi)
	call mod_internal_init(nkn,nel,nlvddi)
	call mod_nohyd_init(nkn,nlvddi)
	call mod_nudging_init(nkn,nel,nlvddi)

	call mod_bclfix_init(nkn,nel,nlvddi)
	call mod_sinking_init(nkn,nlvddi)
	call mod_waves_init(nkn,nel,nlvddi)
	call mod_sedim_init(nkn,nlvddi)
	call mod_bstress_init(nkn)

	write(6,*) '3D arrays allocated: ',nkn,nel,ngr,nlvddi

	end

c*****************************************************************

	subroutine run_scalar()
	
!$	use omp_lib	!ERIC
	
	use mod_bfm

	implicit none
	
	real getpar
	
	integer :: nscal,itemp,isalt,iconz,itvd
	
	nscal = 0
	
	itemp = nint(getpar("itemp"))
	isalt = nint(getpar("isalt"))
	iconz = nint(getpar("iconz"))
	itvd = nint(getpar("itvd"))
	
	call tvd_init(itvd)
	
	if(itemp .gt. 0) nscal = nscal +1
	if(isalt .gt. 0) nscal = nscal +1
	if(iconz .gt. 0) nscal = nscal + iconz

!$      !!!call omp_set_nested(.TRUE.)
	
!$OMP PARALLEL

!$OMP SINGLE 

!$OMP TASKWAIT
!!!$OMP TASKGROUP

!$OMP TASK 
	call barocl(1)
!$OMP END TASK

!$OMP TASK IF ( iconz > 0 )
	call tracer_compute
!$OMP END TASK

!$OMP TASK IF ( ibfm > 0 )
	call bfm_run
!$OMP END TASK

!!!$OMP END TASKGROUP	
!$OMP TASKWAIT	

!$OMP END SINGLE 

!$OMP END PARALLEL 	

	end subroutine

c*****************************************************************

	subroutine print_file_usage

	use intp_fem_file
	use shympi

	implicit none

	if( .not. shympi_is_master() ) return

	write(6,*) '--------------------------------'
	call useunit(200)
	write(6,*) '--------------------------------'
	call iff_print_info(0)
	write(6,*) '--------------------------------'

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine handle_debug_output(dtime)

! the output should be checked with check_debug

	use mod_debug

	implicit none

	double precision dtime

	logical bdebug
	integer id,ios,iunit
	integer, save :: icall = 0
	double precision, save :: da_out(4) = 0.
	character*80 file

	logical has_output_d, next_output_d
	logical openmp_is_master

	if( icall < 0 ) return

        if( icall == 0 ) then
         if( openmp_is_master() ) then
          call init_output_d('itmdbg','idtdbg',da_out)
          if( has_output_d(da_out) ) then
	    call shy_make_output_name('.dbg',file)
	    iunit = 200
            call find_unit(iunit)
            if( iunit == 0 ) goto 98
            open(iunit,file=file,status='unknown',form='unformatted'
     +                          ,iostat=ios)
            if( ios /= 0 ) goto 99
	    call set_debug_unit(iunit)
            call info_output_d('debug_output',da_out)
            icall = 1
	  else
            icall = -1
          end if
	 end if
        end if

        if( next_output_d(da_out) ) then
          !id = nint(da_out(4))
	  !call shympi_debug_output(dtime)
	  call debug_output(dtime)
	end if

	return
   98	continue
        stop 'error stop handle_debug_output: cannot get unit number'
   99	continue
        stop 'error stop handle_debug_output: cannot open file'
	end

c*****************************************************************

	subroutine shympi_debug_output(dtime)

	use shympi_debug
	use mod_depth
	use mod_ts
	use mod_hydro_baro
	use mod_hydro_vel
	use mod_hydro
	use mod_internal
	use levels
	use basin

	implicit none

	double precision dtime

	integer, save :: icall = 0

	write(6,*) 'shympi_debug_output: writing records'

	if( icall == 0 ) then
	  call shympi_write_debug_init
	  call shympi_write_debug_time(dtime)
	  call shympi_write_debug_record('ipv',ipv)
	  call shympi_write_debug_record('ipev',ipev)
	  call shympi_write_debug_record('xgv',xgv)
	  call shympi_write_debug_record('ygv',ygv)
	  call shympi_write_debug_record('fcorv',fcorv)
	else
	  call shympi_write_debug_time(dtime)
	end if

	icall = icall + 1

	!call shympi_write_debug_record('zenv',zenv)
	call shympi_write_debug_record('znv',znv)
	call shympi_write_debug_record('unv',unv)
	call shympi_write_debug_record('vnv',vnv)

	call shympi_write_debug_final

	end

c*****************************************************************

	subroutine debug_output(dtime)

	use mod_debug
	use mod_meteo
	use mod_waves
	use mod_internal
	use mod_depth
	use mod_layer_thickness
	use mod_ts
	use mod_roughness
	use mod_diff_visc_fric
	use mod_hydro_vel
	use mod_hydro
	use mod_geom_dynamic
	use mod_area
	use mod_bound_geom
	use mod_bound_dynamic
	use levels
	use basin

	implicit none

	double precision dtime

	write(6,*) 'debug_output: writing records'

	call write_debug_time(dtime)

	call write_debug_record(ilhkv,'ilhkv')
	call write_debug_record(ilhv,'ilhv')
	call write_debug_record(iwegv,'iwegv')

	call write_debug_record(hm3v,'hm3v')
	call write_debug_record(xgv,'xgv')
	call write_debug_record(ygv,'ygv')

	!call write_debug_record(zeov,'zeov')
	call write_debug_record(zenv,'zenv')
	!call write_debug_record(zov,'zov')
	call write_debug_record(znv,'znv')

	!call write_debug_record(utlov,'utlov')
	!call write_debug_record(vtlov,'vtlov')
	call write_debug_record(utlnv,'utlnv')
	call write_debug_record(vtlnv,'vtlnv')

        call write_debug_record(saltv,'saltv')
        call write_debug_record(tempv,'tempv')
	call write_debug_record(visv,'visv')
	call write_debug_record(difv,'difv')
	!call write_debug_record(wlov,'wlov')
	call write_debug_record(wlnv,'wlnv')

	call write_debug_record(z0bk,'z0bk')
	call write_debug_record(tauxnv,'tauxnv')
	call write_debug_record(tauynv,'tauynv')

	!call write_debug_record(hdeov,'hdeov')
        !call write_debug_record(hdkov,'hdkov')
        call write_debug_record(hdenv,'hdenv')
        call write_debug_record(hdknv,'hdknv')

        call write_debug_record(hdknv,'shearf2')
        call write_debug_record(hdknv,'buoyf2')

        !call write_debug_record(fxv,'fxv')
        !call write_debug_record(fyv,'fyv')
        call write_debug_record(wavefx,'wavefx')
        call write_debug_record(wavefy,'wavefy')
        !call write_debug_record(rfricv,'rfricv')

        !call write_debug_record(momentxv,'momentxv')
        !call write_debug_record(momentyv,'momentyv')

        !call write_debug_record(mfluxv,'mfluxv')
        call write_debug_record(rhov,'rhov')
        call write_debug_record(areakv,'areakv')

	call write_debug_final

	end

c*****************************************************************

	subroutine debug_output_old(dtime)

	use mod_meteo
	use mod_waves
	use mod_internal
	use mod_depth
	use mod_layer_thickness
	use mod_ts
	use mod_roughness
	use mod_diff_visc_fric
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin

	implicit none

	double precision dtime

	write(66) dtime

	call debug_output_record(3*nel,3,hm3v,'hm3v')
	call debug_output_record(nkn,1,xgv,'xgv')
	call debug_output_record(nkn,1,ygv,'ygv')

	call debug_output_record(3*nel,3,zeov,'zeov')
	call debug_output_record(3*nel,3,zenv,'zenv')
	call debug_output_record(nkn,1,zov,'zov')
	call debug_output_record(nkn,1,znv,'znv')

	call debug_output_record(nlvdi*nel,nlvdi,utlov,'utlov')
	call debug_output_record(nlvdi*nel,nlvdi,vtlov,'vtlov')
	call debug_output_record(nlvdi*nel,nlvdi,utlnv,'utlnv')
	call debug_output_record(nlvdi*nel,nlvdi,vtlnv,'vtlnv')

        call debug_output_record(nlvdi*nkn,nlvdi,saltv,'saltv')
        call debug_output_record(nlvdi*nkn,nlvdi,tempv,'tempv')
	call debug_output_record((nlvdi+1)*nkn,nlvdi+1,visv,'visv')
	call debug_output_record((nlvdi+1)*nkn,nlvdi+1,wlov,'wlov')
	call debug_output_record((nlvdi+1)*nkn,nlvdi+1,wlnv,'wlnv')

	call debug_output_record(nkn,1,z0bk,'z0bk')
	call debug_output_record(nkn,1,tauxnv,'tauxnv')
	call debug_output_record(nkn,1,tauynv,'tauynv')

	call debug_output_record(nlvdi*nel,nlvdi,hdeov,'hdeov')
        call debug_output_record(nlvdi*nel,nlvdi,hdenv,'hdenv')
        call debug_output_record(nlvdi*nkn,nlvdi,hdkov,'hdkov')
        call debug_output_record(nlvdi*nkn,nlvdi,hdknv,'hdknv')

        call debug_output_record(nlvdi*nkn,nlvdi,hdknv,'shearf2')
        call debug_output_record(nlvdi*nkn,nlvdi,hdknv,'buoyf2')

        call debug_output_record(nlvdi*nel,nlvdi,fxv,'fxv')
        call debug_output_record(nlvdi*nel,nlvdi,fyv,'fyv')
        call debug_output_record(nlvdi*nel,nlvdi,wavefx,'wavefx')
        call debug_output_record(nlvdi*nel,nlvdi,wavefy,'wavefy')
        call debug_output_record(nel,1,rfricv,'rfricv')

        call debug_output_record(nlvdi*nkn,nlvdi,momentxv,'momentxv')
        call debug_output_record(nlvdi*nkn,nlvdi,momentyv,'momentyv')
        !call debug_output_record((nlvdi+1)*nkn,nlvdi+1,vts,'vts')

	write(66) 0,0

	end

c*****************************************************************

	subroutine debug_output_record(ntot,nfirst,val,text)
	implicit none
	integer ntot,nfirst
	real val(ntot)
	character*(*) text
	character*80 text1
	text1=text
	write(66) ntot,nfirst
	write(66) text1
	write(66) val
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine check_layer_depth

	use mod_depth
	use mod_layer_thickness
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,l,lmax,kmin,lmin
	integer iunit
	real htot,hmin

	iunit = 6
	iunit = 167

	htot = maxval(hev)

	hmin = htot
	kmin = 0
	lmin = 0
	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    if( hdknv(l,k) < hmin ) then
	      hmin = hdknv(l,k)
	      kmin = k
	      lmin = l
	    end if
	  end do
	end do

	write(iunit,*) 'hdknv: ',htot,hmin,kmin,lmin

	end

c*****************************************************************

	subroutine check_max_depth

	use mod_depth
	use mod_layer_thickness
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,l,lmax,ie
	real hmax

	hmax = 0.
	do ie=1,nel
	  hmax = max(hmax,hev(ie))
	end do
	write(6,*) 'hmax (hev) = ',hmax

	hmax = 0.
	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    hmax = max(hmax,hdknv(l,k))
	  end do
	end do
	write(6,*) 'hmax (hdknv) = ',hmax

	end

c*****************************************************************

	subroutine check_special

	use mod_ts
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,l,lmax
	double precision dtime

	call get_act_dtime(dtime)

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    if( saltv(l,k) .gt. 40. ) then
		write(66,*) dtime,l,k,saltv(l,k)
	    end if
	  end do
	end do

	end

c*****************************************************************

	subroutine check_point(text)

	use basin
	use mod_hydro
	use mod_bndo

	implicit none

	character*(*) text

	integer k

	return

	write(6,*) '========= start of check_point ========='
	write(6,*) text
	call mod_bndo_info
	call check_ilevels
	write(6,*) 'znv: ',size(znv),minval(znv),maxval(znv)
	write(6,*) 'zov: ',size(zov),minval(zov),maxval(zov)
	write(6,*) '========= end of check_point ========='

	end

c*****************************************************************

	subroutine test_shympi_arrays

! tests mpi exchange of arrays - can be deleted in some while

	use basin
	use shympi

	implicit none

	integer, allocatable :: local(:)
	integer, allocatable :: global(:)
	integer, allocatable :: i3(:,:)

	integer k,ie,ii

	return

	call shympi_barrier
	call shympi_syncronize

	allocate(local(nkn))
	allocate(global(nkn_global))

	local = 3 + my_id

        call shympi_exchange_array(local,global)

	!if( my_id == 0 ) then
	write(6,*) 'global: ',my_id,nkn,nkn_global
	write(6,*) global
	!end if
	write(6,*) 'local is: ',local(1)

        allocate(i3(3,nel_global))
        call shympi_exchange_array(nen3v,i3)
	!write(6,*) 'global nen3v: ',my_id,nel,nel_global,i3,' -----'

	write(6,*) 'checking nen3v: ',my_id,nkn_global
	do ie=1,nel
	  do ii=1,3
	    k = i3(ii,ie)
	    if( k < 1 .or. k > nkn_global ) then
	      write(6,*) 'nen3v error: ',my_id,k,ie,ii
	    end if
	  end do
	end do

	call shympi_stop('finish')
 
	end

c*****************************************************************

        subroutine test_forcing(dtime,dtend)

        implicit none

        double precision dtime,dtend

        real dt

	do while( dtime .lt. dtend )

           call set_timestep		!sets dt and t_act
           call get_timestep(dt)
	   call get_act_dtime(dtime)

	   call sp111(2)		!boundary conditions
	   
	   call print_time			!output to terminal

        end do

        stop 'end of testing forcing'

        end

c*****************************************************************

	subroutine mpi_debug(dtime)

	use basin
	use shympi

        implicit none

        double precision dtime

	integer nn,ne,i
	integer, allocatable :: ipglob(:)
	integer, allocatable :: ieglob(:)
	real, allocatable :: rpglob(:)
	real, allocatable :: reglob(:)
	real, allocatable :: rp(:)
	real, allocatable :: re(:)
	integer, save :: icall = 0

	integer ipint

	icall = 1
	if( icall > 0 ) return

	icall = icall + 1

	nn = nkn_global
	ne = nel_global
	allocate(ipglob(nn),ieglob(ne))
	allocate(rpglob(nn),reglob(ne))
	allocate(rp(nkn),re(nel))

	rp = ipv
	re = ipev
	write(6,*) nkn_domains
	write(6,*) nel_domains
	flush(6)
	write(6,*) 'exchanging ipv',my_id,size(ipv)
	call shympi_exchange_array(ipv,ipglob)
	call shympi_exchange_array(ipev,ieglob)
	call shympi_exchange_array(rp,rpglob)
	call shympi_exchange_array(re,reglob)

	if( .not. shympi_is_master() ) return

	write(77,*) dtime
	write(77,*) nn
	write(77,*) ipglob
	write(77,*) ne
	write(77,*) ieglob

	write(78,*) dtime

	write(78,*) 'node: ',nn
	do i=1,nn,333
	write(78,*) i,ipglob(i)
	end do

	write(78,*) 'elem: ',ne
	do i=1,ne,333
	write(78,*) i,ieglob(i)
	end do

	write(78,*) 'rnode: ',nn
	do i=1,nn,333
	write(78,*) i,rpglob(i)
	end do

	write(78,*) 'relem: ',ne
	do i=1,ne,333
	write(78,*) i,reglob(i)
	end do

        end

c*****************************************************************

