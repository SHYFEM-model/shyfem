c
c finite element model shyfem (version 3D)
c
c original version from march 1991
c
c revision log :
c
c revised 30.08.95      $$AUST - austausch coefficient introduced
c revised 11.10.95      $$BCLBND - boundary condition for baroclinic runs
c revised 04.08.97      $$ZEONV - new arrays for water level elementwise
c 19.03.1998	ggu	$$IPCCV close data items commented or deleted
c 03.04.1998	ggu	$$DESCRP BUG overwriting descrp with boundary vals.
c 30.04.1998    ggu     finally eliminated /semimp/, /trock/, /ffloat/
c 05.05.1998	ggu	izeit eliminated
c 28.05.1998	ggu	call to sp131g changed
c 14.08.1998	ggu	call to biocos (biological model) introduced
c 20.08.1998	ggu	iextpo finally eliminated, momentum arrays added
c 20.08.1998	ggu	spb11 -> sp111
c 21.08.1998    ggu     xv eliminated
c 03.09.1998    ggu     biological reactor integratated
c 06.11.1998    ggu     hv renamed into hkv
c 08.04.1999    ggu     equilibrium tide (tidal potential) introduced
c 19.04.1999    ggu     subroutine custom called
c 22.10.1999    ggu     igrv, ngrv eliminated (subst by ilinkv, lenkv)
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 07.03.2000    ggu     arrays for vertical turbulent diffusion coefficient
c 20.06.2000    ggu     alv -> visv    slv -> difv
c 19.11.2001    ggu     h1hv eliminated, ulv, vlv eliminated
c 03.12.2001    ggu     new arrays [cst]difhv
c 21.08.2002    ggu     new arrays /kvolc/, /ivol/ copied from hp
c 31.07.2003    ggu     array winv eliminated
c 10.08.2003    ggu     big restructuring
c 13.08.2003    ggu     some more restructuring
c 14.08.2003    ggu     even more restructuring and cleaning up
c 04.12.2003    ggu     integration of wave and sediment module
c 06 03 2004    aac     lagrangian trajectories computation module
c 03.09.2004    ggu     restart now in ht
c 15.10.2004    ggu     gotm substituted by general turbulence closure
c 02.12.2004    ggu     variable time step implemented
c 17.01.2005    ggu     new routine diff_h_set, new difhv, ausv deleted
c 24.02.2005    ggu     new routine smagorinski and common smagv
c 03.03.2005    ggu     smagorinski uses difhv, inserted in subdif.f
c 03.03.2005    ggu     new 3d boundary arrays for C/T/S
c 12.08.2005    ggu     minimum level index, TVD arrays
c 07.11.2005    ggu     new array sed2dn introduced (file name for sediments)
c 16.02.2006    ggu     new routine atoxi3d (TOXI)
c 23.03.2006    ggu     changed time step to real
c 18.10.2006    ccf     radiation stress and waves included (with pipe)
c 10.11.2006    ggu     initialize depth values after restart
c 16.11.2006    ggu     turbulence values included
c 02.04.2007    ggu     changes in algorithm (ASYM)
c 31.05.2007    deb     new arrays bpresxv, bclevvar (debora)
c 26.09.2007    ggu     deleted arrays rcv,rtv,rsv
c 27.09.2007    ggu     deleted call to tstvol,tstvol1
c 20.03.2008    acc     new call for ERSEM ecological model (BFM MODULE)
c 07.04.2008    acc     new array bfm*bc introduced (file name for ersem)
c 10.04.2008    ggu&ccf	upro, waveov, stokes, z0bk
c 16.04.2008    ggu     evaporation mass flux (evapv)
c 22.04.2008    ggu     gradx/yv non global due to parallelization
c 23.04.2008    ggu     deleted r3{c|t|s}v
c 28.04.2008    ggu     new routine init_stability(), new conzm3sh()
c 29.04.2008    ggu&aac new BMF ERSEM model
c 12.11.2008    ggu	new sigma level initialization
c 10.12.2008    ggu	new array rfricv for bottom friction
c 18.12.2008    ggu	new routine debug_output(), might be removed later
c 04.03.2009    ggu	matrix amat into solver routines
c 11.03.2009    ggu	new arrays for meteo data
c 06.04.2009    ggu	renamed nlidim to nlkdim
c 30.11.2009    ggu	write output file for successfull completion
c 19.02.2010    ggu	init_stability() changed to reset_stability()
c 22.02.2010    ggu	new arrays wxv, wyv
c 26.02.2010    ggu	new arrays sauxe1/2
c 29.04.2010    ggu	write volumes (wrfvla)
c 26.01.2011    ggu	new arrays for observations and nudging
c 16.02.2011    ggu	new iarnv, call to aquabc
c 17.02.2011    ccf	new radiation stress in 3D
c 23.03.2011    ggu	new call to adjust_spherical()
c 31.03.2011    ggu	write finite volumes at initial time step
c 20.05.2011    ggu	iwetv introduced, wet and dry from main
c 25.10.2011    ggu	hlhv eliminated
c 18.11.2011    ggu	new routine handle_projection
c 24.01.2012    ggu	new call to setup_parallel()
c 23.02.2012    ggu&ccf	meteo arrays adjusted (3*nkn)
c 09.03.2012    ggu	call to residence time added
c 21.06.2012    ggu&aar	fluid mud variables integrated
c 05.08.2012    ggu	bug because lam2dn and dmfd2n not defined
c 10.05.2013    dbf	initialization for non hydrostatic routines
c 13.06.2013    ggu	set/copydepth simplified, offline version
c 05.09.2013    ggu	changed order of adjust depth and barene structures
c 29.10.2013    ggu	nudging implemented
c 25.03.2014    ggu	new offline
c 25.06.2014    ggu	new arrays hkv_min, hkv_max
c 05.12.2014    ccf	new interface for waves
c 30.07.2015    ggu	routine renamed from ht to shyfem
c 18.09.2015    ggu	new routine scalar, call to hydro()
c 29.09.2015    ccf	inverted set_spherical() and handle_projection()
c 10.10.2015    ggu	fluid mud routines handled differently
c 05.10.2017    ggu	command line options introduced, subs rearranged
c
c*****************************************************************

c----------------------------------------------------------------

	program shyfem

	use mod_bound_geom
	use mod_geom
	use mod_meteo
	use mod_waves
	use mod_turbulence
	use mod_sinking
	use mod_nudging
	use mod_internal
	use mod_geom_dynamic
	use mod_depth
	use mod_bnd_aux
	use mod_gotm_aux
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
	use tidef
	use projection
	use coordinates
	use mod_subset
	use mod_bfm
!$	use omp_lib	!ERIC
	use shympi

c----------------------------------------------------------------

	implicit none

c include files

	include 'mkonst.h'
	include 'pkonst.h'
	!include 'femtime.h'

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
	character*80 strfile

	real getpar

	call cpu_time(time1)
	call system_clock(count1, count_rate, count_max)
!$      timer = omp_get_wtime() 


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%% code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c-----------------------------------------------------------
c copyright and command line options
c-----------------------------------------------------------

        call shyfem_init(strfile,bdebug,bdebout,bmpirun)

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
	call setup_time			!in case start time has changed

	!call init_vertical	!do again after restart

	call get_act_dtime(dtime)
	call get_first_dtime(dtanf)
	call get_last_dtime(dtend)

	call setnod

	call set_area

	call make_new_depth
	call copy_depth
	call make_new_depth
	!call check_max_depth

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
	call barocl(0)
	call wrfvla		!write finite volume
	call nonhydro_init
	call init_wave		!waves
	call initsed		!sediments

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
	call bfm_init
	call renewal_time

	call submud_init

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

	!call custom(it)		!call for initialization

	!write(6,*) 'starting time loop'
        call shympi_comment('starting time loop...')

	call print_time

	call check_parameter_values('before main')

	if( bdebout ) call debug_output(dtime)

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

	   call renewal_time
	   call ecological_module	!ecological model
           call atoxi3d			!toxi
           call mercury_module

           call lagrange

	   call offline(1)		!write to offline file

	   call do_after

	   call tracer_write
	   call bfm_write

           if( bfirst ) call print_file_usage

	   call print_time			!output to terminal

	   call total_energy
	   call check_all
	   !call check_layer_depth
	   !call check_special

           call write_wwm

	   if( bdebout ) call debug_output(dtime)
	   bfirst = .false.

	end do

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%% end of time loop %%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call check_parameter_values('after main')

	call print_end_time

	if( shympi_is_master() ) then

!$OMP PARALLEL
!$OMP MASTER
	nthreads = 1
!$	nthreads = omp_get_num_threads()
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

	print *,"NUMBER OF MPI THREADS USED  = ",n_threads

        mpi_t_end = shympi_wtime()
        write(6,*)'MPI_TIME =',mpi_t_end-mpi_t_start,my_id
        write(6,*)'Parallel_TIME =',mpi_t_end-parallel_start,my_id
	call shympi_time_get(1,mpi_t_solve)
        write(6,*)'MPI_SOLVE_TIME =',mpi_t_solve,my_id

	end if

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

        character*80 version

	if( shympi_is_master() ) then
	  call shyfem_copyright('shyfem - 3D hydrodynamic SHYFEM routine')
	end if

	call get_shyfem_version(version)
        call clo_init('shyfem','str-file',trim(version))

        call clo_add_info('runs the 3D shyfem routine')

        call clo_add_option('debug',.false.,'enable debugging')
        call clo_add_option('debout',.false.
     +			,'writes debugging information to file')
        call clo_add_option('mpi',.false.
     +			,'runs in MPI mode (experimental)')

        call clo_parse_options

        call clo_get_option('debug',bdebug)
        call clo_get_option('debout',bdebout)
        call clo_get_option('mpi',bmpirun)

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
	use tidef
	use coordinates
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	call tidef_init(nkn)
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
	use mod_turbulence
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

	call mod_layer_thickness_init(nkn,nel,nlvddi)
	call mod_internal_init(nkn,nel,nlvddi)
	call mod_nohyd_init(nkn,nlvddi)
	call mod_nudging_init(nkn,nel,nlvddi)

	call mod_bclfix_init(nkn,nel,nlvddi)
	!call mod_fluidmud_init(nkn,nlvddi)
	call mod_sinking_init(nkn,nlvddi)
	call mod_turbulence_init(nkn,nlvddi)
	call mod_waves_init(nkn,nel,nlvddi)
	call mod_sedim_init(nkn,nlvddi)

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
	 call bfm_compute
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

	subroutine debug_output(dtime)

	use mod_meteo
	use mod_waves
	use mod_internal
	use mod_depth
	use mod_layer_thickness
	use mod_gotm_aux
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

	include 'femtime.h'
	
	integer k,l,lmax

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    if( saltv(l,k) .gt. 40. ) then
		write(66,*) it,l,k,saltv(l,k)
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

