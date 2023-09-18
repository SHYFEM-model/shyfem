!
! $Id: ht.f,v 1.76 2010-03-11 15:36:38 georg Exp $
!
! finite element model shympi (MPI version)
!
! original version from march 1991
!
! revision log :
!
! revised 30.08.95      $$AUST - austausch coefficient introduced
! revised 11.10.95      $$BCLBND - boundary condition for baroclinic runs
! revised 04.08.97      $$ZEONV - new arrays for water level elementwise
! 19.03.1998	ggu	$$IPCCV close data items commented or deleted
! 03.04.1998	ggu	$$DESCRP BUG overwriting descrp with boundary vals.
! 30.04.1998    ggu     finally eliminated /semimp/, /trock/, /ffloat/
! 05.05.1998	ggu	izeit eliminated
! 28.05.1998	ggu	call to sp131g changed
! 14.08.1998	ggu	call to biocos (biological model) introduced
! 20.08.1998	ggu	iextpo finally eliminated, momentum arrays added
! 20.08.1998	ggu	spb11 -> sp111
! 21.08.1998    ggu     xv eliminated
! 03.09.1998    ggu     biological reactor integratated
! 06.11.1998    ggu     hv renamed into hkv
! 08.04.1999    ggu     equilibrium tide (tidal potential) introduced
! 19.04.1999    ggu     subroutine custom called
! 22.10.1999    ggu     igrv, ngrv eliminated (subst by ilinkv, lenkv)
! 20.01.2000    ggu     common block /dimdim/ eliminated
! 07.03.2000    ggu     arrays for vertical turbulent diffusion coefficient
! 20.06.2000    ggu     alv -> visv    slv -> difv
! 19.11.2001    ggu     h1hv eliminated, ulv, vlv eliminated
! 03.12.2001    ggu     new arrays [cst]difhv
! 21.08.2002    ggu     new arrays /kvolc/, /ivol/ copied from hp
! 31.07.2003    ggu     array winv eliminated
! 10.08.2003    ggu     big restructuring
! 13.08.2003    ggu     some more restructuring
! 14.08.2003    ggu     even more restructuring and cleaning up
! 04.12.2003    ggu     integration of wave and sediment module
! 06 03 2004    aac     lagrangian trajectories computation module
! 03.09.2004    ggu     restart now in ht
! 15.10.2004    ggu     gotm substituted by general turbulence closure
! 02.12.2004    ggu     variable time step implemented
! 17.01.2005    ggu     new routine diff_h_set, new difhv, ausv deleted
! 24.02.2005    ggu     new routine smagorinski and common smagv
! 03.03.2005    ggu     smagorinski uses difhv, inserted in subdif.f
! 03.03.2005    ggu     new 3d boundary arrays for C/T/S
! 12.08.2005    ggu     minimum level index, TVD arrays
! 07.11.2005    ggu     new array sed2dn introduced (file name for sediments)
! 16.02.2006    ggu     new routine atoxi3d (TOXI)
! 23.03.2006    ggu     changed time step to double precision
! 18.10.2006    ccf     radiation stress and waves included (with pipe)
! 10.11.2006    ggu     initialize depth values after restart
! 16.11.2006    ggu     turbulence values included
! 02.04.2007    ggu     changes in algorithm (ASYM)
! 31.05.2007    deb     new arrays bpresxv, bclevvar (debora)
! 26.09.2007    ggu     deleted arrays rcv,rtv,rsv
! 27.09.2007    ggu     deleted call to tstvol,tstvol1
! 20.03.2008    acc     new call for ERSEM ecological model (BFM MODULE)
! 07.04.2008    acc     new array bfm*bc introduced (file name for ersem)
! 10.04.2008    ggu&ccf	upro, waveov, stokes, z0bk
! 16.04.2008    ggu     evaporation mass flux (evapv)
! 22.04.2008    ggu     gradx/yv non global due to parallelization
! 23.04.2008    ggu     deleted r3{c|t|s}v
! 28.04.2008    ggu     new routine init_stability(), new conzm3sh()
! 29.04.2008    ggu&aac new BMF ERSEM model
! 12.11.2008    ggu	new sigma level initialization
! 10.12.2008    ggu	new array rfricv for bottom friction
! 18.12.2008    ggu	new routine debug_output(), might be removed later
! 04.03.2009    ggu	matrix amat into solver routines
! 11.03.2009    ggu	new arrays for meteo data
! 06.04.2009    ggu	renamed nlidim to nlkdim
! 30.11.2009    ggu	write output file for successfull completion
! 19.02.2010    ggu	init_stability() changed to reset_stability()
! 22.02.2010    ggu	new arrays wxv, wyv
! 26.02.2010    ggu	new arrays sauxe1/2
! 29.04.2010    ggu	write volumes (wrfvla)
! 26.01.2011    ggu	new arrays for observations and nudging
! 16.02.2011    ggu	new iarnv, call to aquabc
! 17.02.2011    ccf	new radiation stress in 3D
! 23.03.2011    ggu	new call to adjust_spherical()
! 31.03.2011    ggu	write finite volumes at initial time step
! 20.05.2011    ggu	iwetv introduced, wet and dry from main
! 25.10.2011    ggu	hlhv eliminated
! 18.11.2011    ggu	new routine handle_projection
! 24.01.2012    ggu	new call to setup_parallel()
! 23.02.2012    ggu&ccf	meteo arrays adjusted (3*nkn)
! 09.03.2012    ggu	call to residence time added
! 21.06.2012    ggu&aar	fluid mud variables integrated
! 05.08.2012    ggu	bug because lam2dn and dmfd2n not defined
! 10.05.2013    dbf	initialization for non hydrostatic routines
! 13.06.2013    ggu	set/copydepth simplified, offline version
! 05.09.2013    ggu	changed order of adjust depth and barene structures
! 29.10.2013    ggu	nudging implemented
! 25.03.2014    ggu	new offline
! 25.06.2014    ggu	new arrays hkv_min, hkv_max
! 05.12.2014    ccf	new interface for waves
! 30.07.2015    ggu	routine renamed from ht to shyfem
! 18.09.2015    ggu	new routine scalar, call to hydro()
! 29.09.2015    ccf	inverted set_spherical() and handle_projection()
! 10.10.2015    ggu	fluid mud routines handled differently
!
!*****************************************************************

!----------------------------------------------------------------

	program shympi_main

	use bnd_geom
	use geom
	use meteo
	use waves
	use turbulence_util
	use sinking
	use nudging
	use internal
	use geom_dynamic
	use depth
	use bnd_aux
	use gotm_aux
	use diff_aux
	use bnd_dynamic
	use area
	use ts
	use roughness
	use diffusion
	use hydro_baro
	use hydro_print
	use hydro_vel
	use hydro_admin
	use evgeom
	use levels
	use basin
	use intp_fem_file
	use tidef
	use projection
	use coordinates
	use subset
	use shympi
        use petsc_admin
	use mpi_graph_elem
        use para
!$	use omp_lib	!ERIC
        use layer_thickness !hdeov/hdenv/hdkov/hdknv
        use version
        use elems_dealing
        use bndo_admin
        use topological
        use restart
        use ecological
        use mud_admin
        use chezy
        use depth_util
        use waves_admin
        use wetdry
        use sedim_admin
	use shy_turbulence
        use custom_admin
        use check
        use openmp_admin
        use iostr
        use baroclinic
        use offline_data,       only:   offline
        use closing
        use nohydro
        use diffusion_admin
        use bcvel
        use lagrange_main
        use fil
        use nudge
        use system
        use hydrodynamic
        use initialize
        use befor_after
        use bnd_routines
        use stability
        use time_util
        use time_admin
        use constants
        use f_vol
        use mkdistance
        use conz_admin
        use water_ren_time
        use toxical
        use timing

!----------------------------------------------------------------

! include files
        implicit none

	include 'param.h'
	include 'mkonst.h'
	include 'pkonst.h'
	include 'femtime.h'

! local variables

	logical bdebout
	integer iwhat,levdbg
	integer date,time
	integer nthreads
	integer*8 count1,count2,count3, count_rate, count_max
	double precision time1,time2
	double precision timer,timer2
        double precision mpi_t_start, mpi_t_end
        double precision parallel_start

        double precision time_scalars,t_scal    !ivb

	double precision ahpar,azpar,ampar
        double precision dt
        integer ic,n
	integer numlevels,totalLevels

        logical bout

	bdebout = .false.

	call cpu_time(time1)
	call system_clock(count1, count_rate, count_max)
!$      timer = omp_get_wtime() 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%% code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!-----------------------------------------------------------
! dimensions
!-----------------------------------------------------------

!-----------------------------------------------------------
! read STR file
!-----------------------------------------------------------

	call cstinit
	call cstfile				!read STR and basin

	call setup_omp_parallel

	call shympi_init(.true.)
        
        if (shympi_is_master()) call shyfem_copyright('3D FEM model') !ivb


        mpi_t_start = shympi_wtime()
        if(bmpi) then
	  call shympi_setup			!sets up partitioning of basin
          call shympi_alloc
        end if
        parallel_start = shympi_wtime()

        azpar = getpar('azpar')
        ampar = getpar('ampar')

        if(.not.((azpar == 0. .and. ampar == 1.).or.(azpar==1. .and. ampar == 0.))) then
          call petsc_init(nkndi,nkn_inner)
        end if

	call allocate_2d_arrays

!-----------------------------------------------------------
! check dimensions
!-----------------------------------------------------------

	levdbg = nint(getpar('levdbg'))
	bdebout = ( levdbg == -1 )

!-----------------------------------------------------------
! check parameters read and set time and Coriolis
!-----------------------------------------------------------

	call cstcheck

!	call pritst(1)

!-----------------------------------------------------------
! initialize triangles
!-----------------------------------------------------------

	call set_spherical
	call set_ev
	call adjust_spherical
	call print_spherical
	call handle_projection
	call set_geom
	call domain_clusterization
        
!-----------------------------------------------------------
! inititialize time independent vertical arrays
!-----------------------------------------------------------

	call adjust_depth	!adjusts hm3v
	call init_vertical	!makes nlv,hlv,hldv,ilhv,ilhkv, adjusts hm3v

	
        call countLevel(numlevels, totalLevels)
        write(6,*)'nodes,elems,lvs:',nkn,nel,numlevels,neldi,nkndi,nlvdi

!-----------------------------------------------------------
! allocates arrays
!-----------------------------------------------------------

	call allocate_3d_arrays
	call set_depth			!makes hev,hkv and exchanges

	call check_point('checking ht 1')

!-----------------------------------------------------------
! initialize barene data structures
!-----------------------------------------------------------

	call setweg(-1,n)
	call setnod
	call update_geom	!update ieltv - needs inodv
        
!-----------------------------------------------------------
! initialize boundary conditions
!-----------------------------------------------------------

	call check_point('checking ht 2')

	call get_date_time(date,time)
	call iff_init_global(nkn,nel,nlv,ilhkv,ilhv,hkv_max,hev,hlv,date,time)

	call sp111(1)           !here zenv, utlnv, vtlnv are initialized

!-----------------------------------------------------------
! initialize depth arrays and barene data structure
!-----------------------------------------------------------

	call setweg(0,n)
	call setznv		! -> change znv since zenv has changed

        call inirst             !restart
	call setup_time		!in case itanf has changed

	!call init_vertical	!do again after restart

	call setnod

	call set_area

	call make_new_depth
	call copy_depth
	call make_new_depth
	!call check_max_depth

!-----------------------------------------------------------
! initialize open boundary routines
!-----------------------------------------------------------

	call check_point('checking ht 3')

	call bndo_init
	call bndo_info_file('bndo_info.dat')

!-----------------------------------------------------------
! initialize transports and velocities
!-----------------------------------------------------------

	call init_uv            !set vel, w, pr, ... from transports
	call barocl(0)
	call wrfvla		!write finite volume
	call nonhydro_init
	call init_wave		!waves

!-----------------------------------------------------------
! initialize modules
!-----------------------------------------------------------


	call init_others
	call init_chezy
	call init_nodal_area_code	!area codes on nodes
        call diffweight
        call set_diffusivity
	call tidefini
	call cstsetup
	call sp136(ic)
        call shdist(rdistv)
	call tracer_init
	call renewal_time

	call submud_init

!-----------------------------------------------------------
! write input values to log file and perform check
!-----------------------------------------------------------

	call check_point('checking ht 4')

	call check_fem
	call check_values
	call prilog

	call bclfix_ini

	call system_initialize		!matrix inversion routines

	call offline(2)
	call offline(1)

	call init_nudging

	call do_init

	call custom(it)

        call shympi_comment('starting time loop...')
	call print_time

	call check_parameter_values('before main')

	if( bdebout ) call debug_output(it)

        ahpar = getpar('ahpar')

        if(shympi_partition_on_elements() .and. ahpar .gt. 0) then
          call send_halo(utlov,nlvdi,nel_local,'ut')
          call send_halo(vtlov,nlvdi,nel_local,'vt')
        end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	!if( bmpi ) call shympi_stop('scheduled stop')

	do while( it .lt. itend )

           !call shympi_comment('new time iteration -----------------')

	   if(niter .eq. 2) then
	     call system_clock(count3, count_rate, count_max)
	     timer2 = count3-count1
	     timer2 = timer2 / count_rate
             write(6,*)'INITIALIZATION TIME:', timer2
             call timing_init
	   end if

           if(bmpi_debug) then
	     call shympi_check_all
           end if

	   call check_crc
	   call set_dry

	   call reset_stability

           call set_timestep		!sets idt and it
           call get_timestep(dt)
	   !call compute_stability_stats(1,aux)

	   call do_befor

	   call offline(2)		!read from offline file

	   call sp111(2)		!boundary conditions

           call read_wwm
	   
           if(bmpi_debug) then
	     call shympi_check_all
           end if

	   call hydro			!hydro

	   call wrfvla			!write finite volume - shympi FIXME

	   call scalar
          
           call turb_closure

           call parwaves(it)            !parametric wave model
           call sedi(it,dt)             !sediment transport
	   call submud(it,dt)           !fluid mud (ARON)

	   call renewal_time
	   call ecological_module(it,dt)	!ecological model
           call atoxi3d(it,dt)			!toxi

           call lagrange

	   call offline(1)		!write to offline file

	   call do_after

           if( niter .eq. 1 ) then
	     call useunit(200)
	     call iff_print_info(0)
	     write(6,*) '--------------------------------'
	   end if

	   call tracer_write
	   call print_time			!output to terminal

	   call total_energy
	   call check_all
	   !call check_special

           call write_wwm

	   if( bdebout ) call debug_output(it)

	end do

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%% end of time loop %%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	write(6,*) 'checking parameter values...'
	call check_parameter_values('after main')

	call print_end_time

        print *,"NUMBER OF MPI THREADS USED  = ",n_threads

!$OMP PARALLEL
!$OMP MASTER
	nthreads = 1
!$	nthreads = omp_get_num_threads()
        print *,"NUMBER OF OMP THREADS USED  = ",nthreads
!$OMP END MASTER
!$OMP END PARALLEL

!$      timer = omp_get_wtime() - timer
!$      print *,"TIME TO SOLUTION (OMP)      = ",timer

	call system_clock(count2, count_rate, count_max)
	timer = count2-count1
	timer = timer / count_rate
	print *,"TIME TO SOLUTION (WALL)     = ",timer,my_id

	call cpu_time(time2)
	print *,"TIME TO SOLUTION (CPU)      = ",time2-time1,my_id

	timer2 = count2-count3
	timer2 = timer2 / count_rate
	print *,"TIMESTEPS TIME (WALL)     = ",timer2,my_id
        !call ht_finished

	!call pripar(15)
	!call prifnm(15)

        mpi_t_end = shympi_wtime()

        write(6,*)'MPI_TIME =',mpi_t_end-mpi_t_start,my_id
        write(6,*)'SETUP_MPI_TIME: ',parallel_start-mpi_t_start
        if(ln_timing)   call timing_finalize

        if(.not.((azpar == 0. .and. ampar == 1.).or.(azpar==1. .and. ampar == 0.))) then
          call petsc_final
        end if
	call shympi_finalize
	call exit(0)

        end

!*****************************************************************

        subroutine ht_finished

        implicit none

        open(1,file='SHYFEM_FINISHED',status='unknown',form='formatted')
        write(1,*) 'SHYFEM has finished successfully the simulation'
        close(1)

	end

!*****************************************************************

	subroutine debug_output(it)

	use meteo
	use waves
	use internal
	use depth
	use layer_thickness
	use gotm_aux
	use ts
	use roughness
	use diffusion
	use hydro_vel
	use hydro_admin
	use levels
	use basin

	implicit none

	integer it

	write(66) it

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
	double precision val(ntot)
	character*(*) text
	character*80 text1
	text1=text
	write(66) ntot,nfirst
	write(66) text1
	write(66) val
	end

!*****************************************************************

	subroutine check_max_depth

	use depth
	use layer_thickness
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer k,l,lmax,ie
	double precision hmax

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

	use ts
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

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine allocate_2d_arrays

	use bndo
	use tvd
	use mod_bnd
	use bnd_geom
	use geom
	use meteo
	use geom_dynamic
	use bnd_aux
	use diff_aux
	use roughness
	use hydro_baro
	use depth
	use evgeom
	use tidef
	use coordinates
	use basin, only : nkn,nel,ngr,mbw
        use shympi

	implicit none

	call tidef_init(nkn)
	call coordinates_init(nkn)

#ifdef DEBUGON
	call mod_hydro_baro_init(nel_local)
	call mod_diff_aux_init(nel_local)
#else
	call mod_hydro_baro_init(nel)
	call mod_diff_aux_init(nel)
#endif
	call mod_roughness_init(nkn)

	call mod_bnd_aux_init(nkn,nel)
#ifdef DEBUGON
	call mod_geom_dynamic_init(nkn_local,nel_local)
#else
	call mod_geom_dynamic_init(nkn,nel)
#endif

	call mod_meteo_init(nkn)
	call mod_geom_init(nkn,nel,ngr,nkn_local)
	if(bmpi) then
	  call mod_bndo_init(ngr,nrb,nrb+bounds%nneighbors)
	else
	  call mod_bndo_init(ngr,nrb)
	end if

#ifdef DEBUGON
	call mod_depth_init(nkn_local,nel_local)
#else
	call mod_depth_init(nkn,nel)
#endif

	!call ev_init(nel)
	call ev_init(nel_local)

        if(.not.(shympi_partition_on_elements())) then 
	  call mod_tvd_init(nel)
        end if

        if (shympi_is_master()) then
	  write(6,*) '2D arrays allocated: ',nkn,nel,ngr
        end if
	end

!*****************************************************************

	subroutine allocate_3d_arrays

	use conz_common
	use waves
	use turbulence_util
	use sinking
	!use mod_fluidmud
	use bclfix_util
	use nohyd
	use internal
	use layer_thickness
	use gotm_aux
	use bnd_dynamic
	use nudging
	use area
	use ts
	use diffusion
	use hydro_print
	use hydro_vel
	use hydro_admin
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
        use shympi

	implicit none

	integer nlvddi

	nlvddi = nlvdi
#ifdef DEBUGON
	call mod_hydro_init(nkn_local,nel_local,nlvddi)
	call mod_hydro_vel_init(nkn_local,nel_local,nlvddi)
	call mod_hydro_print_init(nkn_local,nlvddi)
	call mod_diffusion_init(nkn_local,nel_local,nlvddi)
	call mod_ts_init(nkn_local,nlvddi)
	call mod_area_init(nkn_local,nlvddi)
	call mod_bnd_dynamic_init(nkn_local,nlvddi)
#else
	call mod_hydro_init(nkn,nel_local,nlvddi)
	call mod_hydro_vel_init(nkn,nel,nlvddi)
	call mod_hydro_print_init(nkn,nlvddi)
	call mod_diffusion_init(nkn,nel,nlvddi)
        if(shympi_partition_on_elements()) then
	  call mod_ts_init(nkn,nlvddi,nkn_local)
        else
	  call mod_ts_init(nkn,nlvddi)
        end if  
	call mod_area_init(nkn,nlvddi)
	call mod_bnd_dynamic_init(nkn,nlvddi)
#endif

	call mod_gotm_aux_init(nkn,nlvddi)

#ifdef DEBUGON
	call mod_layer_thickness_init(nkn_local,nel_local,nlvddi)
	call mod_internal_init(nkn_local,nel_local,nlvddi)
#else
	call mod_layer_thickness_init(nkn,nel,nlvddi)
	call mod_internal_init(nkn,nel,nlvddi)
#endif
	call mod_nohyd_init(nkn,nlvddi)

#ifdef DEBUGON
	call mod_nudging_init(nkn_local,nel_local,nlvddi)
	call mod_sinking_init(nkn_local,nlvddi)
#else
	call mod_nudging_init(nkn,nel,nlvddi)
	call mod_sinking_init(nkn,nlvddi)
#endif
	call mod_bclfix_util_init(nkn,nel,nlvddi)
	!call mod_fluidmud_init(nkn,nlvddi)
	call mod_turbulence_init(nkn,nlvddi)
	call mod_waves_init(nkn,nel,nlvddi)

	write(6,*) '3D arrays allocated: ',nkn,nel,ngr,nlvddi

	end

!*****************************************************************

	subroutine check_point(text)

	use bndo

	implicit none

	character*(*) text

	return

	write(6,*) '========= start of check_point ========='
	write(6,*) text
	call mod_bndo_info
	call check_ilevels
	write(6,*) '========= end of check_point ========='

	end

!*****************************************************************

	subroutine scalar()
	
!$	use omp_lib	!ERIC

        use shympi
        use para
        use tvd_admin
        use baroclinic
        use conz_admin
        use timing

	implicit none
	
	integer :: nscal,itemp,isalt,iconz,itvd
        double precision s_time

	nscal = 0
	
	itemp = nint(getpar("itemp"))
	isalt = nint(getpar("isalt"))
	iconz = nint(getpar("iconz"))
	itvd = nint(getpar("itvd"))
	
        !gmica togliere il commento all' if per continuare 
        ! le modifiche per la correzione tvd_fluxes
        !if(.not.(shympi_partition_on_elements())) then 
	  call tvd_init(itvd)
        !end if
	
	if(itemp .gt. 0) nscal = nscal +1
	if(isalt .gt. 0) nscal = nscal +1
	if(iconz .gt. 0) nscal = nscal + iconz

!$      !!!call omp_set_nested(.TRUE.)
	
!$OMP PARALLEL

!$OMP SINGLE 

!$OMP TASKGROUP

!$OMP TASK 
	
        if(ln_timing) s_time = shympi_wtime()
	call barocl(1)
        if(ln_timing) scalar_time = scalar_time + shympi_wtime() - s_time
	
	!print *, " end barocl"
!$OMP END TASK
!$OMP TASK
	 !print *, "tracer"
	 call tracer_compute
	 !print *, "end tracer"
!$OMP END TASK
!$OMP END TASKGROUP	
!$OMP END SINGLE 

!$OMP END PARALLEL 	

	end subroutine

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine shympi_check_all

	implicit none

	call shympi_check_all_static
	call shympi_check_all_dynamic
	call shympi_check_all_scalar

	end

!*****************************************************************

	subroutine shympi_check_all_static

	implicit none

	call shympi_check_depth
	call shympi_check_geom_static
	call shympi_check_levels

	end

!*****************************************************************

	subroutine shympi_check_all_dynamic

	implicit none

	call shympi_check_geom_dynamic

	call shympi_check_hydro
	call shympi_check_hydro_baro
	call shympi_check_hydro_vel
	call shympi_check_hydro_print

	end

!*****************************************************************

	subroutine shympi_check_all_scalar

        use conz_common
        use ts
        use shympi

	implicit none

        if( mod_conz_is_initialized() ) then
	  call shympi_check_3d_node(cnv,'cnv')
        end if
	call shympi_check_3d_node(saltv,'saltv')
	call shympi_check_3d_node(tempv,'tempv')
	call shympi_check_3d_node(rhov,'rhov')

	end

!*****************************************************************

	subroutine shympi_check_hydro

	use basin
	use hydro_admin
	use shympi

	implicit none

	integer i
	double precision aux(nel)

	call shympi_check_2d_node(zov,'zov')
	call shympi_check_2d_node(znv,'znv')

	do i=1,3
	  aux = zeov(i,:)
	  call shympi_check_2d_elem(aux,'zeov')
	  aux = zenv(i,:)
	  call shympi_check_2d_elem(aux,'zenv')
	end do

	call shympi_check_3d_elem(utlnv,'utlnv')
	call shympi_check_3d_elem(vtlnv,'vtlnv')
	call shympi_check_3d_elem(utlov,'utlov')
	call shympi_check_3d_elem(vtlov,'vtlov')

	end

!*****************************************************************

	subroutine shympi_check_hydro_baro

	use hydro_baro
	use shympi

	implicit none

	call shympi_check_2d_elem(unv,'unv')
	call shympi_check_2d_elem(vnv,'vnv')
	call shympi_check_2d_elem(uov,'uov')
	call shympi_check_2d_elem(vov,'vov')

	end

!*****************************************************************

	subroutine shympi_check_hydro_vel

	use basin
	use levels
	use hydro_vel
	use shympi

	implicit none

	call shympi_check_3d_elem(ulov,'ulov')
	call shympi_check_3d_elem(vlov,'vlov')
	call shympi_check_3d_elem(ulnv,'ulnv')
	call shympi_check_3d_elem(vlnv,'vlnv')
	call shympi_check_3d0_node(wlnv,'wlnv')
	call shympi_check_3d0_node(wlov,'wlov')

	end

!*****************************************************************

	subroutine shympi_check_hydro_print

	use basin
	use levels
	use hydro_print
	use shympi

	implicit none

	integer i
	double precision aux(nkn)

	call shympi_check_3d_node(uprv,'uprv')
	call shympi_check_3d_node(vprv,'vprv')
	call shympi_check_3d_node(upro,'upro')
	call shympi_check_3d_node(vpro,'vpro')
	call shympi_check_3d0_node(wprv,'wprv')
	call shympi_check_2d_node(up0v,'up0v')
	call shympi_check_2d_node(vp0v,'vp0v')

	do i=1,3
	  aux = xv(i,:)
	  !call shympi_check_2d_node(aux,'xv')
	end do

	end

!*****************************************************************

	subroutine shympi_check_depth

	use depth
	use shympi

	implicit none

	call shympi_check_2d_elem(hev,'hev')
	call shympi_check_2d_node(hkv,'hkv')
	call shympi_check_2d_node(hkv_min,'hkv_min')
	call shympi_check_2d_node(hkv_max,'hkv_max')

	end

!*****************************************************************

	subroutine shympi_check_geom_dynamic

	use evgeom
	use geom_dynamic
	use shympi

	implicit none

	call shympi_check_2d_elem(iwegv,'iwegv')
	call shympi_check_2d_elem(iwetv,'iwetv')
	call shympi_check_2d_node(inodv,'inodv')	!FIXME - not working

	end

!*****************************************************************

	subroutine shympi_check_geom_static

	use basin
	use evgeom
	use geom
	use shympi

	implicit none

	integer i
	double precision aux(nel)

	do i=1,evdim
	  aux = ev(i,:)
	  call shympi_check_2d_elem(aux,'ev')
	end do

	end

!*****************************************************************

	subroutine shympi_check_levels

	use levels
	use shympi

	implicit none

	call shympi_check_2d_elem(ilhv,'ilhv')
	call shympi_check_2d_elem(ilmv,'ilmv')
	call shympi_check_2d_node(ilhkv,'ilhkv')
	call shympi_check_2d_node(ilmkv,'ilmkv')

	end

!*****************************************************************

