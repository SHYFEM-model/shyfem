
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

! finite element model shyfem (version 3D)
!
! original version from march 1991
!
! revision log :
!
! 30.08.1995	ggu	$$AUST - austausch coefficient introduced
! 11.10.1995	ggu	$$BCLBND - boundary condition for baroclinic runs
! 04.08.1997	ggu	$$ZEONV - new arrays for water level elementwise
! 19.03.1998	ggu	$$IPCCV close data items commented or deleted
! 03.04.1998	ggu	$$DESCRP BUG overwriting descrp with boundary vals.
! 30.04.1998	ggu	finally eliminated /semimp/, /trock/, /ffloat/
! 05.05.1998	ggu	izeit eliminated
! 28.05.1998	ggu	call to sp131g changed
! 14.08.1998	ggu	call to biocos (biological model) introduced
! 20.08.1998	ggu	iextpo finally eliminated, momentum arrays added
! 20.08.1998	ggu	spb11 -> sp111
! 21.08.1998	ggu	xv eliminated
! 03.09.1998	ggu	biological reactor integratated
! 06.11.1998	ggu	hv renamed into hkv
! 08.04.1999	ggu	equilibrium tide (tidal potential) introduced
! 19.04.1999	ggu	subroutine custom called
! 22.10.1999	ggu	igrv, ngrv eliminated (subst by ilinkv, lenkv)
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 07.03.2000	ggu	arrays for vertical turbulent diffusion coefficient
! 20.06.2000	ggu	alv -> visv    slv -> difv
! 19.11.2001	ggu	h1hv eliminated, ulv, vlv eliminated
! 03.12.2001	ggu	new arrays [cst]difhv
! 21.08.2002	ggu	new arrays /kvolc/, /ivol/ copied from hp
! 31.07.2003	ggu	array winv eliminated
! 10.08.2003	ggu	big restructuring
! 13.08.2003	ggu	some more restructuring
! 14.08.2003	ggu	even more restructuring and cleaning up
! 04.12.2003	ggu	integration of wave and sediment module
! 06.03.2004	aac	lagrangian trajectories computation module
! 03.09.2004	ggu	restart now in ht
! 15.10.2004	ggu	gotm substituted by general turbulence closure
! 02.12.2004	ggu	variable time step implemented
! 17.01.2005	ggu	new routine diff_h_set, new difhv, ausv deleted
! 24.02.2005	ggu	new routine smagorinski and common smagv
! 03.03.2005	ggu	smagorinski uses difhv, inserted in subdif.f
! 03.03.2005	ggu	new 3d boundary arrays for C/T/S
! 12.08.2005	ggu	minimum level index, TVD arrays
! 07.11.2005	ggu	new array sed2dn introduced (file name for sediments)
! 16.02.2006	ggu	new routine atoxi3d (TOXI)
! 23.03.2006	ggu	changed time step to real
! 18.10.2006	ccf	radiation stress and waves included (with pipe)
! 10.11.2006	ggu	initialize depth values after restart
! 16.11.2006	ggu	turbulence values included
! 02.04.2007	ggu	changes in algorithm (ASYM)
! 31.05.2007	dbf	new arrays bpresxv, bclevvar (debora)
! 26.09.2007	ggu	deleted arrays rcv,rtv,rsv
! 27.09.2007	ggu	deleted call to tstvol,tstvol1
! 20.03.2008	aac	new call for ERSEM ecological model (BFM MODULE)
! 07.04.2008	aac	new array bfm*bc introduced (file name for ersem)
! 10.04.2008	ggu&ccf	upro, waveov, stokes, z0bk
! 16.04.2008	ggu	evaporation mass flux (evapv)
! 22.04.2008	ggu	gradx/yv non global due to parallelization
! 23.04.2008	ggu	deleted r3{c|t|s}v
! 28.04.2008	ggu	new routine init_stability(), new conzm3sh()
! 29.04.2008	ggu&aac	new BMF ERSEM model
! 12.11.2008	ggu	new sigma level initialization
! 10.12.2008	ggu	new array rfricv for bottom friction
! 18.12.2008	ggu	new routine debug_output(), might be removed later
! 04.03.2009	ggu	matrix amat into solver routines
! 11.03.2009	ggu	new arrays for meteo data
! 06.04.2009	ggu	renamed nlidim to nlkdim
! 30.11.2009	ggu	write output file for successfull completion
! 19.02.2010	ggu	init_stability() changed to reset_stability()
! 22.02.2010	ggu	new arrays wxv, wyv
! 26.02.2010	ggu	new arrays sauxe1/2
! 29.04.2010	ggu	write volumes (wrfvla)
! 26.01.2011	ggu	new arrays for observations and nudging
! 16.02.2011	ggu	new iarnv, call to aquabc
! 17.02.2011	ccf	new radiation stress in 3D
! 23.03.2011	ggu	new call to adjust_spherical()
! 31.03.2011	ggu	write finite volumes at initial time step
! 20.05.2011	ggu	iwetv introduced, wet and dry from main
! 25.10.2011	ggu	hlhv eliminated
! 18.11.2011	ggu	new routine handle_projection
! 24.01.2012	ggu	new call to setup_parallel()
! 23.02.2012	ggu&ccf	meteo arrays adjusted (3*nkn)
! 09.03.2012	ggu	call to residence time added
! 21.06.2012	ggu&aar	fluid mud variables integrated
! 05.08.2012	ggu	bug because lam2dn and dmfd2n not defined
! 10.05.2013	dbf	initialization for non hydrostatic routines
! 13.06.2013	ggu	set/copydepth simplified, offline version
! 05.09.2013	ggu	changed order of adjust depth and barene structures
! 29.10.2013	ggu	nudging implemented
! 25.03.2014	ggu	new offline
! 25.06.2014	ggu	new arrays hkv_min, hkv_max
! 05.12.2014	ccf	new interface for waves
! 30.07.2015	ggu	routine renamed from ht to shyfem
! 18.09.2015	ggu	new routine scalar, call to hydro()
! 23.09.2015	ggu	changed VERS_7_2_4
! 29.09.2015	ccf	inverted set_spherical() and handle_projection()
! 30.09.2015	ggu	changed VERS_7_2_6
! 10.10.2015	ggu	fluid mud routines handled differently
! 12.10.2015	ggu	changed VERS_7_3_3
! 13.10.2015	ggu	changed VERS_7_3_5
! 22.10.2015	ggu	changed VERS_7_3_7
! 23.10.2015	ggu	changed VERS_7_3_9
! 05.11.2015	ggu	changed VERS_7_3_12
! 09.11.2015	ggu	changed VERS_7_3_13
! 20.11.2015	ggu	changed VERS_7_3_15
! 16.12.2015	ggu	changed VERS_7_3_16
! 19.02.2016	ggu	changed VERS_7_5_2
! 22.02.2016	ggu	changed VERS_7_5_4
! 07.06.2016	ggu	changed VERS_7_5_12
! 14.06.2016	ggu	changed VERS_7_5_14
! 13.04.2017	ggu	changed VERS_7_5_25
! 09.05.2017	ggu	changed VERS_7_5_26
! 16.05.2017	ggu	changed VERS_7_5_27
! 05.10.2017	ggu	command line options introduced, subs rearranged
! 14.11.2017	ggu	changed VERS_7_5_36
! 17.11.2017	ggu	changed VERS_7_5_37
! 05.12.2017	ggu	changed VERS_7_5_39
! 07.12.2017	ggu	changed VERS_7_5_40
! 03.04.2018	ggu	changed VERS_7_5_43
! 19.04.2018	ggu	changed VERS_7_5_45
! 26.04.2018	ggu	changed VERS_7_5_46
! 11.05.2018	ggu	changed VERS_7_5_47
! 06.07.2018	ggu	changed VERS_7_5_48
! 25.10.2018	ggu	changed VERS_7_5_51
! 18.12.2018	ggu	changed VERS_7_5_52
! 12.02.2019	ccf	bottom shear stress in substress.f
! 16.02.2019	ggu	changed VERS_7_5_60
! 12.03.2019	ccf	include new computation of tide potential/analysis
! 04.07.2019	ggu	new ww3 routines introduced
! 15.09.2019	ggu	subroutine to test only forcing
! 02.10.2019	ggu	delete include files
! 17.10.2019	ggu	no call to bfm_write, is done inside subroutine
! 06.11.2019	ggu	femelab eliminated
! 03.04.2020	ggu	write real start and end time of simulation
! 09.04.2020    ggu     run bfm through bfm_run()
! 21.05.2020    ggu     better handle copyright notice
! 04.06.2020    ggu     debug_output() substituted with shympi_debug_output()
! 30.03.2021    ggu     more on debug, call sp111(2) outside time loop
! 01.04.2021    ggu     turbulence cleaned
! 02.04.2022    ggu     new option -mpi_debug -> calls shympi_check_all()
! 02.04.2022    ggu     new routine shympi_write_debug_special()
! 03.04.2022    ggu     timing problems in handle_debug_output() solved
! 11.04.2022    ggu     no -mpi switch necessary anymore
! 12.04.2022    ggu     message to show if mpi support is available
! 18.05.2022    ggu     cpu_time routines introduced
! 10.03.2023    ggu     do not use bmpirun anymore
! 28.04.2023    ggu     update function calls for belem
! 22.05.2023    ggu     new names for closing: close_init, close_handle
! 05.06.2023    lrp     introduce z-star
! 11.12.2023    ggu     f90 style format
! 11.12.2023    ggu     create subroutines for init/run/finalize
! 12.12.2023    ggu     introduce dtmax to make limited run
!
!*****************************************************************
!
! notes :
!
! MPI calls
!
! call shympi_init(bmpirun)
!   call shympi_alloc_global(nkn,nel,nen3v,ipv,ipev)
! call shympi_setup			!sets up partitioning of basin
!
!----------------------------------------------------------------

	program shyfem
	call shyfem_main
	end program

!================================================================
	module mod_shyfem
!================================================================

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
#ifdef W3_SHYFEM
	use subww3
#endif

	implicit none

! internal variables

	logical bdebout,bdebug,bmpirun
	logical, save ::  bfirst = .true.
	integer k,ic,n
	integer iwhat
	integer date,time
	integer nthreads
	integer*8 count1,count2,count3,count_rate,count_max
	real time0,time1,time2,time3
	real dt
	double precision timer
	double precision mpi_t_start,mpi_t_end,parallel_start
	double precision mpi_t_solve
	double precision dtime,dtanf,dtend
        double precision atime_start,atime_end
        double precision, parameter :: zero = 0.
        character*20 aline_start,aline_end
	character*80 strfile
	character*80 mpi_code,omp_code

!================================================================
	end module mod_shyfem
!================================================================

	subroutine shyfem_main
	use mod_shyfem
	implicit none
	call shyfem_initialize
	call shyfem_run(zero)
	call shyfem_finalize
	end subroutine shyfem_main

!*****************************************************************

	subroutine shyfem_initialize

!-----------------------------------------------------------
! start of program
!-----------------------------------------------------------

	use mod_shyfem

	implicit none

	call cpu_time(time1)
	call cpu_time_init
	call cpu_time_start(1)
	call system_clock(count1, count_rate, count_max)
!$      timer = omp_get_wtime() 

        call get_real_time(atime_start,aline_start)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%% code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------------------------------------
! copyright and command line options
!-----------------------------------------------------------

        call shyfem_init(strfile,bdebug,bdebout,bmpirun)

!-----------------------------------------------------------
! write  real start time
!-----------------------------------------------------------

        write(6,*) 'simulation start:   ',aline_start 

!-----------------------------------------------------------
! read STR file
!-----------------------------------------------------------

	call cstinit
	call cstfile(strfile)			!read STR and basin

	call shympi_init(.true.)
	call setup_omp_parallel

	call cpu_time(time3)
	call system_clock(count3, count_rate, count_max)
        mpi_t_start = shympi_wtime()
	call shympi_setup			!sets up partitioning of basin
        parallel_start = shympi_wtime()

	call allocate_2d_arrays

!-----------------------------------------------------------
! check parameters read and set time and Coriolis
!-----------------------------------------------------------

	call cstcheck

!	call pritst(1)

!-----------------------------------------------------------
! initialize triangles
!-----------------------------------------------------------

	!call levels_init_2d(nkn,nel)	!maybe not needed
	call set_spherical
	call set_ev
	call adjust_spherical
	call print_spherical
	call handle_projection
	call set_geom
	call set_geom_mpi		!adjusts boundaries
	call shympi_barrier
	call domain_clusterization	!create subsets for OMP
	call shympi_barrier

!-----------------------------------------------------------
! inititialize time independent vertical arrays
!-----------------------------------------------------------

	call adjust_depth	!adjusts hm3v
	call init_vertical	!makes nlv,hlv,hldv,ilhv,ilhkv, adjusts hm3v

!-----------------------------------------------------------
! allocates arrays
!-----------------------------------------------------------

	call allocate_3d_arrays
	call set_depth		!makes hev,hkv and exchanges

	call check_point('checking ht 1')

	call poisson_init

!-----------------------------------------------------------
! initialize barene data structures
!-----------------------------------------------------------

	if( .not. bmpi ) call setweg(-1,n)	!shympi - FIXME
	call setnod
	call update_geom	!update ieltv - needs inodv
	call populate_strings	!populate strings here

!-----------------------------------------------------------
! initialize boundary conditions
!-----------------------------------------------------------

	call check_point('checking ht 2')

	call get_date_time(date,time)
        call iff_init_global(nkn,nel,nlv,ilhkv,ilhv                     &
     &                          ,hkv_max,hev,hlv,date,time)

	call init_zadaptation
	call sp111(1)           !here zenv, utlnv, vtlnv are initialized

!-----------------------------------------------------------
! initialize depth arrays and barene data structure
!-----------------------------------------------------------

	call setweg(0,n)
	call setznv		! -> change znv since zenv has changed

        call rst_perform_restart        !restart
	call compute_velocities
	call copy_uvz

	!call init_vertical	!do again after restart

	call setnod
	call set_area

	call initialize_layer_depth
	!call check_max_depth

	call setup_time		!in case start time has changed with rst
	call get_act_dtime(dtime)
	call get_first_dtime(dtanf)
	call get_last_dtime(dtend)

!-----------------------------------------------------------
! initialize open boundary routines
!-----------------------------------------------------------

	call check_point('checking ht 3')

	call bndo_init
	if( bdebug ) then
	  call bndo_info_file('bndo_info.dat')
	end if

!-----------------------------------------------------------
! initialize transports and velocities
!-----------------------------------------------------------

	call init_uv            !set vel, w, pr, ... from transports
	call copy_uvz		!copy new to old
	call barocl(0)
	call wrfvla		!write finite volume
	call nonhydro_init
	call init_wave		!waves
	call ww3_init
	call initsed		!sediments
        call init_bstress	!bottom shear stress

!-----------------------------------------------------------
! initialize modules
!-----------------------------------------------------------

	call init_others
	call init_chezy
	call init_nodal_area_code	!area codes on nodes
        call diffweight
        call set_diffusivity
	call tidefini
	call close_init
        call shdist(rdistv)
	call tracer_init
        call qhdist(qdistv) !DWNH
	call bfm_init
	call renewal_time
        call lagrange
	call tidepar_init
	call submud_init
	call handle_gotm_init
	call tripple_points_handle

	call cstsetup

!-----------------------------------------------------------
! write input values to log file and perform check
!-----------------------------------------------------------

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
        !call test_zeta_debug(dtime)

	call test_zeta_init
	call cpu_time_start(2)

	end subroutine shyfem_initialize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine shyfem_run(dtstep)

	use mod_shyfem

	implicit none

	double precision dtstep		!time step to run, 0: run to end

	double precision dtmax		!run to this time

	dtmax = dtend
	if( dtstep > 0. ) dtmax = dtime + dtstep

	do while( dtime .lt. dtmax )

           if(bmpi_debug) call shympi_check_all(0)	!checks arrays

	   call check_crc
	   call set_dry

           call set_timestep(dtmax)		!sets dt and t_act
           call get_timestep(dt)
	   call get_act_dtime(dtime)

	   call do_befor

	   call copy_uvz		!copies new to old time level
	   call nonhydro_copy   	!copies non hydrostatic pressure terms
	   call copy_depth		!copies layer depth to old

	   call offline(2)		!read from offline file
	   call sp111(2)		!boundary conditions
           call read_wwm		!wwm wave model
	   
           if(bmpi_debug) call shympi_check_all(1)	!checks arrays

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

	end subroutine shyfem_run

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%% end of time loop %%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine shyfem_finalize

	use mod_shyfem

	implicit none

	call cpu_time_end(2)

	call print_end_time
	call ww3_finalize
	call system_finalize		!matrix inversion routines

	call check_parameter_values('after main')

	call cpu_time_end(1)

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
        print *,"TIME TO SOLUTION PARALLEL REGION (CPU) = "             &
     &                          ,time2-time3,my_id

	call cpu_time_get(1,time0)
	write(6,1200) 'cpu time (total):        ',time0
	call cpu_time_get(2,time0)
	write(6,1200) 'cpu time in time loop:   ',time0
	call cpu_time_get(3,time1)
	write(6,1200) 'cpu time solving matrix: ',time1
	write(6,1200) 'cpu time no matrix:      ',time0 - time1
 1200	format(a,f12.3)

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
	call shympi_barrier
	call shympi_exit(99)
	call exit(99)

        end subroutine shyfem_finalize

!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine ht_finished

        implicit none

        open(1,file='SHYFEM_FINISHED',status='unknown',form='formatted')
        write(1,*) 'SHYFEM has finished successfully the simulation'
        close(1)

	end

!*****************************************************************

        subroutine shyfem_init(strfile,bdebug,bdebout,bmpirun)

        use clo
        use shympi
	use mod_zeta_system

        implicit none

        character*(*) strfile
        logical bdebug,bdebout,bmpirun,bmpidebug
        logical bquiet,bsilent

        character*80 version

	logical openmp_is_parallel

	call get_shyfem_version_and_commit(version)
        call clo_init('shyfem','str-file',trim(version))

        call clo_add_info('runs the 3D shyfem routine')

	call clo_add_sep('general options:')
        call clo_add_option('quiet',.false.,'do not be verbose')
        call clo_add_option('silent',.false.,'be silent')

	call clo_add_sep('mpi options:')
        call clo_add_option('mpi',.false.                               &
     &                  ,'runs in MPI mode (experimental)')

	call clo_add_sep('debug options:')
        call clo_add_option('debug',.false.,'enable debugging')
        call clo_add_option('debout',.false.                            &
     &                  ,'writes debugging information to file')
        call clo_add_option('mpi_debug',.false.,'enable mpi debugging')

        call clo_parse_options

        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)

        call clo_get_option('mpi',bmpirun)	!not used anymore

        call clo_get_option('debug',bdebug)
        call clo_get_option('debout',bdebout)
        call clo_get_option('mpi_debug',bmpidebug)

        if( bsilent ) bquiet = .true.
        !if( bmpirun ) call shympi_set_debug(bmpidebug)
        call shympi_set_debug(bmpidebug)

	if( shympi_is_master() ) then
         call shyfem_set_short_copyright(bquiet)
         if( .not. bsilent ) then
	  call shyfem_copyright('shyfem - 3D hydrodynamic SHYFEM routine')
	  if( bmpi_support ) then
	    write(6,*) 'compiled with parallel support: MPI'
	  else if( openmp_is_parallel() ) then
	    write(6,*) 'compiled with parallel support: OMP'
	  else
	    write(6,*) 'compiled with parallel support: NONE'
	  end if
	  write(6,*) 'matrix solver: ',trim(solver_type)
	  write(6,*)
         end if
	end if

        call clo_check_files(1)
        call clo_get_file(1,strfile)

        end

!*****************************************************************
!*****************************************************************
!*****************************************************************

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

!*****************************************************************

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

!*****************************************************************

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

!*****************************************************************

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

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine handle_debug_output(dtime)

! the output should be checked with check_debug

	use mod_debug
	use shympi_debug
	use shympi

	implicit none

	double precision dtime

	logical bdebug,blast
	integer id,ios
	integer, save :: iunit = 0
	integer, save :: icall = 0
	double precision, save :: da_out(4) = 0.
	character*80 file

	logical has_output_d, next_output_d

	if( icall < 0 ) return

	call shympi_barrier

        if( icall == 0 ) then
	  write(6,*) 'setting up handle_debug_output',my_id
          if( shympi_is_master() ) then
            call init_output_d('itmdbg','idtdbg',da_out)
	    call assure_initial_output_d(da_out)
            if( has_output_d(da_out) ) then
	      call shy_make_output_name('.dbg',file)
	      iunit = 200
              call find_unit(iunit)
              if( iunit == 0 ) goto 98
              open(iunit,file=file,status='unknown',form='unformatted'  &
     &                          ,iostat=ios)
              if( ios /= 0 ) goto 99
	      call set_debug_unit(iunit)
	      call shympi_write_debug_unit(iunit)
              !call info_output_d('debug_output',da_out)
              icall = 1
	    else
              icall = -1
            end if
	  end if
	  call shympi_bcast(icall)
	  call shympi_bcast(da_out)
	  !write(6,*) 'after broadcast: ',icall,da_out
	  !write(6,*) 'finished setting up handle_debug_output',my_id
        end if

        if( next_output_d(da_out) ) then
	  call shympi_debug_output(dtime)
	  !call debug_output(dtime)		!serial
	end if

	call is_time_last(blast)
	if( blast .and. shympi_is_master() ) then
	  if( iunit > 0 ) then
	    flush(iunit)
	    close(iunit)
	    iunit = 0
	  end if
	end if

	call shympi_barrier

	return
   98	continue
        stop 'error stop handle_debug_output: cannot get unit number'
   99	continue
        stop 'error stop handle_debug_output: cannot open file'
	end

!*****************************************************************

	subroutine shympi_debug_output(dtime)

	use shympi
	use shympi_debug
	use mod_depth
	use mod_ts
	use mod_hydro_baro
	use mod_hydro_vel
	use mod_hydro_print
	use mod_hydro
	use mod_internal
	use mod_conz
	use levels
	use basin
	use mod_layer_thickness
	use mod_diff_visc_fric
	use mod_bound_dynamic
	use tide

	implicit none

	double precision dtime,dtime1
	real, allocatable :: aux3d(:,:)

	integer, save :: icall = 0

	if( shympi_is_master() ) then
	  write(6,*) 'shympi_debug_output: writing records'
	end if

	call shympi_write_debug_init		!can be called more than once

	if( icall == 0 ) then
	  dtime1 = -1.
	  call shympi_write_debug_time(dtime1)
	  call shympi_write_debug_node('ipv',ipv)
	  call shympi_write_debug_elem('ipev',ipev)
	  call shympi_write_debug_node('xgv',xgv)
	  call shympi_write_debug_node('ygv',ygv)
	  call shympi_write_debug_elem('fcorv',fcorv)
	  call shympi_write_debug_elem(3,'hm3v',hm3v)
	  call shympi_write_debug_special
	  call shympi_write_debug_final
	end if

	icall = icall + 1
	call shympi_write_debug_time(dtime)

	!call shympi_write_debug_node('rqv',rqv)
	!call shympi_write_debug_node('rqpsv',rqpsv)
	!call shympi_write_debug_node('rqdsv',rqdsv)
	!call shympi_write_debug_node('mfluxv',mfluxv)

	call shympi_write_debug_node('zeqv',zeqv)
	call shympi_write_debug_node('znv',znv)
	call shympi_write_debug_node('zov',zov)
	call shympi_write_debug_elem(3,'zenv',zenv)
	call shympi_write_debug_elem('unv',unv)
	call shympi_write_debug_elem('vnv',vnv)
	call shympi_write_debug_elem('utlnv',utlnv)
	call shympi_write_debug_elem('vtlnv',vtlnv)
	call shympi_write_debug_node('hdknv',hdknv)
	call shympi_write_debug_elem('hdenv',hdenv)
	call shympi_write_debug_node('saltv',saltv)
	call shympi_write_debug_node('tempv',tempv)
	call shympi_write_debug_node('rhov',rhov)
	call shympi_write_debug_node('wlnv',wlnv)
	call shympi_write_debug_elem('fxv',fxv)
	call shympi_write_debug_elem('fyv',fyv)
	if( allocated(cnv) ) then
	  call shympi_write_debug_node('cnv',cnv)
	end if

	!call shympi_write_debug_elem(3,'zeov',zeov)
	!call shympi_write_debug_elem('hdeov',hdeov)
	!call shympi_write_debug_elem('utlov',utlov)
	!call shympi_write_debug_elem('vtlov',vtlov)
	!call shympi_write_debug_node('momentxv',momentxv)
	!call shympi_write_debug_node('momentyv',momentyv)

	allocate(aux3d(nlvdi,nkn))
	aux3d(:,:) = visv(1:nlvdi,:)
	call shympi_write_debug_node('visv',aux3d)
	aux3d(:,:) = difv(1:nlvdi,:)
	call shympi_write_debug_node('difv',aux3d)

	call shympi_write_debug_node('uprv',uprv)
	call shympi_write_debug_node('vprv',vprv)

	call shympi_write_debug_final

	end

!*****************************************************************

	subroutine shympi_write_debug_special

! writes specialy constructed  arrays to check correctness of exchange

	use basin
	use levels
	use shympi
	use shympi_debug

	implicit none

	integer k,ie,l,lmax,id
	integer, parameter :: ifact = 100000

	integer, allocatable :: ian(:,:)
	integer, allocatable :: iae(:,:)
	real, allocatable :: ran(:,:)
	real, allocatable :: rae(:,:)

	integer ipext,ieext

	allocate(ian(nlvdi,nkn),ran(nlvdi,nkn))
	allocate(iae(nlvdi,nel),rae(nlvdi,nel))
	ian = 0
	iae = 0
	ran = 0.
	rae = 0.

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    id = ifact*ipext(k) + l
	    ian(l,k) = id
	    ran(l,k) = id
	  end do
	end do

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    id = ifact*ieext(ie) + l
	    iae(l,ie) = id
	    rae(l,ie) = id
	  end do
	end do

	call shympi_write_debug_node('ian',ian)
	call shympi_write_debug_elem('iae',iae)
	call shympi_write_debug_node('ran',ran)
	call shympi_write_debug_elem('rae',rae)

	end

!*****************************************************************

	subroutine debug_output_old1(dtime)

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

!*****************************************************************

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

!*****************************************************************

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

!*****************************************************************
!*****************************************************************
!*****************************************************************

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

!*****************************************************************

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

!*****************************************************************

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

!*****************************************************************

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

!*****************************************************************

        subroutine test_forcing(dtime,dtend)

        implicit none

        double precision dtime,dtend

        real dt

	do while( dtime .lt. dtend )

           call set_timestep(dtend)		!sets dt and t_act
           call get_timestep(dt)
	   call get_act_dtime(dtime)

	   call sp111(2)		!boundary conditions
	   
	   call print_time			!output to terminal

        end do

        stop 'end of testing forcing'

        end

!*****************************************************************

	subroutine mpi_debug(dtime)

! writes debug information for ipv, ipev, etc.. to file

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

	icall = 1			!do not run the debug routine
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
	call shympi_l2g_array(ipv,ipglob)
	call shympi_l2g_array(ipev,ieglob)
	call shympi_l2g_array(rp,rpglob)
	call shympi_l2g_array(re,reglob)

	if( .not. shympi_is_master() ) return

	write(77,*) 'mpi_debug:'
	write(77,*) dtime
	write(77,*) nn
	write(77,*) ipglob
	write(77,*) ne
	write(77,*) ieglob

	write(78,*) 'mpi_debug:'
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

!*****************************************************************

	subroutine test_vertical(itime)

! test routine for problems in vertical velocities

	use basin
	use levels
	use mod_hydro_vel
	use mod_bound_dynamic
	use shympi

	implicit none

	integer itime

	integer kext,kint,iu,lmax
	integer ipint

	kext = 2249
	kint = ipint(kext)
	if( kint <= 0 ) return

	iu = 667 + my_id
	lmax = ilhkv(kint)
	write(iu,*) itime,lmax,kint
	write(iu,*) '    ',wlnv(0:lmax,kint)
	if( allocated(mfluxv) ) then
	write(iu,*) '    ',mfluxv(1:lmax,kint)
	end if
	flush(iu)

	end

!*****************************************************************

	subroutine test_zeta_debug(dtime)

! test routine for problems in zeta init

	use basin
	use mod_hydro
	use shympi

	implicit none

	double precision dtime
	integer kext,kint,iu
	integer ipint

	kext = 1318
	iu = 880 + my_id

	kint = ipint(kext)
	if( kint <= 0 ) return

	write(iu,*) 'test_zeta_debug: ',kext,kint,znv(kint)

	end

!*****************************************************************

	subroutine test_scalar_debug(text,ip,dtime,scal)

! test routine for problems in temp

	use basin
	use levels
	use mod_hydro
	use mod_ts
	use mod_diff_visc_fric
	use shympi

	implicit none

	character*(*) text
	integer ip
	double precision dtime
	real scal(nlvdi,nkn)

	integer k,l
	integer iu,nkng,nlvg,ifreq

	real, allocatable :: rauxg(:,:)

	if( dtime <= 3600. ) return
	if( dtime > 3700. ) return

	allocate(rauxg(nlv_global,nkn_global))

	call shympi_l2g_array(scal,rauxg)
	
	iu = 456
	nkng = nkn_global
	nlvg = nlv_global
	ifreq = nkng/100
	ifreq = nkng/30

	write(iu,*) dtime,ip

	do k=1,nkng,ifreq
	  !do l=1,nlvg,3
	  do l=1,1
	    write(iu,*) trim(text),ip,k,l,rauxg(l,k)
	  end do
	end do

	end

!*****************************************************************

