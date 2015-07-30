c
c $Id: ht.f,v 1.76 2010-03-11 15:36:38 georg Exp $
c
c finite element model ht (version 3D)
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
c
c*****************************************************************

c----------------------------------------------------------------

	program shyfem

	use mod_tides
	use mod_bound_geom
	use mod_geom
	use mod_meteo
	use mod_waves
	use mod_turbulence
	use mod_sinking
	use mod_fluidmud
	use mod_nudging
	use mod_internal
	use mod_geom_dynamic
	use mod_depth
	use mod_bnd_aux
	use mod_gotm_aux
	use mod_diff_aux
	use mod_bound_dynamic
	use mod_aux_array
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

c----------------------------------------------------------------

c include files

	include 'param.h'
	include 'mkonst.h'
	include 'pkonst.h'
	include 'femtime.h'

c local variables

	logical bdebout
	integer iwhat,levdbg
	integer date,time

	real getpar

	bdebout = .false.

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%% code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call shyfem_copyright('3D FEM model')

c-----------------------------------------------------------
c dimensions
c-----------------------------------------------------------

c-----------------------------------------------------------
c read STR file
c-----------------------------------------------------------

	call cstinit
	call cstfile				!read STR and basin

	call bas_get_para(nkn,nel,ngr,mbw)	!to be deleted later

	call allocate_2d_arrays

c-----------------------------------------------------------
c check dimensions
c-----------------------------------------------------------

	levdbg = nint(getpar('levdbg'))
	bdebout = ( levdbg == -1 )

c-----------------------------------------------------------
c check parameters read and set time and Coriolis
c-----------------------------------------------------------

	call setup_parallel

	call cstcheck

c	call pritst(1)

c-----------------------------------------------------------
c initialize triangles
c-----------------------------------------------------------

	call handle_projection
	call set_spherical
	call set_ev
	call adjust_spherical
	call print_spherical
	call set_geom

c-----------------------------------------------------------
c inititialize time independent vertical arrays
c-----------------------------------------------------------

	call adjust_depth	!adjusts hm3v
	call init_vertical	!makes nlv,hlv,hldv,ilhv,ilhkv, adjusts hm3v

c-----------------------------------------------------------
c allocates arrays
c-----------------------------------------------------------

	call allocate_3d_arrays
	call set_depth		!makes hev,hkv

	call check_point('checking ht 1')

c-----------------------------------------------------------
c initialize barene data structures
c-----------------------------------------------------------

	call setweg(-1,n)
	call setnod
	call update_geom	!update ieltv - needs inodv

c-----------------------------------------------------------
c initialize boundary conditions
c-----------------------------------------------------------

	call check_point('checking ht 2')

	call get_date_time(date,time)
	call iff_init_global(nkn,nlv,ilhkv,hkv_max,hlv,date,time)

	call sp111(1)           !here zenv, utlnv, vtlnv are initialized

	!call handle_bsig_init	!initialize from sigma level data?

c-----------------------------------------------------------
c initialize depth arrays and barene data structure
c-----------------------------------------------------------

	call setweg(0,n)
	call setznv		! -> change znv since zenv has changed

        call inirst             !restart
	call setup_time		!in case itanf has changed

	!call init_vertical	!do again after restart

	call setnod

	call setarea(nlvdi,areakv)

	call make_new_depth
	call copy_depth
	call make_new_depth

	!call check_max_depth

c-----------------------------------------------------------
c initialize open boundary routines
c-----------------------------------------------------------

	call check_point('checking ht 3')

	call bndo_init
	call bndo_info_file('bndo_info.dat')

c-----------------------------------------------------------
c initialize transports and velocities
c-----------------------------------------------------------

	call init_uv            !set vel, w, pr, ... from transports
	call barocl(0)
	call wrfvla		!write finite volume
	call nonhydro_init
	call init_wave		!waves

c-----------------------------------------------------------
c initialize modules
c-----------------------------------------------------------

	call init_others
	call init_chezy
	call init_nodal_area_code	!area codes on nodes
        call diffweight
        call diff_h_set
	call tideini
	call cstsetup
	call sp136(ic)
        call shdist(rdistv)
	call tracer
	call renewal_time


c-----------------------------------------------------------
c write input values to log file and perform check
c-----------------------------------------------------------

	call check_point('checking ht 4')

	call check_fem
	call check_values
c	call listdim
	call prilog

c        call bclevvar_ini       	!chao debora
	call bclfix_ini

	call system_initialize		!matrix inversion routines

	call offline(2)
	call offline(1)

	call nudge_init

	call do_init

	!call custom(it)		!call for initialization

	write(6,*) 'starting time loop'
	call print_time

	call check_parameter_values('before main')

	if( bdebout ) call debug_output(it)

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%% time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do while( it .lt. itend )

	   !if( it .eq. 1429200 ) call test3d(66,100)

	   call check_crc
	   call set_dry

	   call reset_stability

           call set_timestep		!sets idt and it
           call get_timestep(dt)
	   !call compute_stability_stats(1,aux)

	   call do_befor

	   call offline(2)		!read from offline file

	   call sp111(2)		!boundary conditions
	   call nudge_zeta

           call read_wwm
	   
	   call sp259f			!hydro

	   call wrfvla			!write finite volume

	   call tracer
	   call barocl(1)		!baroclinic contribution

	   !call compute_heat_flux

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
	   end if

	   call print_time			!output to terminal

	   call total_energy
	   call check_all
	   !call check_special

           call write_wwm

	   if( bdebout ) call debug_output(it)

	   !if( it .gt. 50400 ) call test3d(66,100)

	end do

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%% end of time loop %%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	write(6,*) 'checking parameter values...'
	call check_parameter_values('after main')

	call print_end_time
        !call ht_finished

	!call pripar(15)
	!call prifnm(15)

	call exit(99)

        end

c*****************************************************************

        subroutine ht_finished

        implicit none

        open(1,file='SHYFEM_FINISHED',status='unknown',form='formatted')
        write(1,*) 'SHYFEM has finished successfully the simulation'
        close(1)

	end

c*****************************************************************

	subroutine debug_output(it)

	use mod_meteo
	use mod_waves
	use mod_fluidmud
	use mod_internal
	use mod_depth
	use mod_layer_thickness
	use mod_gotm_aux
	use mod_aux_array
	use mod_ts
	use mod_roughness
	use mod_diff_visc_fric
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin

	implicit none

	integer it

	include 'param.h'

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

        call debug_output_record(nlvdi*nkn,nlvdi,saux2,'saux2')
        call debug_output_record(nlvdi*nkn,nlvdi,saux3,'saux3')
        call debug_output_record((nlvdi+1)*nkn,nlvdi+1,vts,'vts')

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

	subroutine check_max_depth

	use mod_depth
	use mod_layer_thickness
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'





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

	include 'param.h'

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
c*****************************************************************
c*****************************************************************

	subroutine allocate_2d_arrays

	use mod_bndo
	use mod_tides
	use mod_tvd
	use mod_bnd
	use mod_bound_geom
	use mod_geom
	use mod_meteo
	use mod_nudging
	use mod_geom_dynamic
	use mod_bnd_aux
	use mod_diff_aux
	use mod_roughness
	use mod_hydro_baro
	use mod_depth
	use evgeom
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	call mod_tides_init(nel)

	call mod_hydro_baro_init(nel)
	call mod_roughness_init(nkn)
	call mod_diff_aux_init(nel)

	call mod_bnd_aux_init(nkn,nel)
	call mod_geom_dynamic_init(nkn,nel)
	call mod_nudging_init(nkn)

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
	use mod_turbulence
	use mod_sinking
	use mod_fluidmud
	use mod_bclfix
	use mod_nohyd
	use mod_internal
	use mod_layer_thickness
	use mod_gotm_aux
	use mod_bound_dynamic
	use mod_aux_array
	use mod_area
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	integer nlvddi

	nlvddi = nlvdi

	call mod_hydro_init(nkn,nel,nlvddi)
	call mod_hydro_vel_init(nkn,nel,nlvddi)
	call mod_hydro_print_init(nkn,nlvddi)

	call mod_diff_visc_fric_init(nkn,nel,nlvddi)

	call mod_ts_init(nkn,nlvddi)

	call mod_area_init(nkn,nlvddi)
	call mod_bound_dynamic_init(nkn,nlvddi)
	call mod_aux_array_init(nkn,nel,nlvddi)
	call mod_gotm_aux_init(nkn,nlvddi)

	call mod_layer_thickness_init(nkn,nel,nlvddi)
	call mod_internal_init(nkn,nel,nlvddi)
	call mod_nohyd_init(nkn,nlvddi)

	call mod_bclfix_init(nkn,nel,nlvddi)
	call mod_fluidmud_init(nkn,nlvddi)
	call mod_sinking_init(nkn,nlvddi)
	call mod_turbulence_init(nkn,nlvddi)
	call mod_waves_init(nkn,nel,nlvddi)

	write(6,*) '3D arrays allocated: ',nkn,nel,ngr,nlvddi

	end

c*****************************************************************

	subroutine check_point(text)

	use mod_bndo

	implicit none

	character*(*) text

	return

	write(6,*) '========= start of check_point ========='
	write(6,*) text
	call mod_bndo_info
	call check_ilevels
	write(6,*) '========= end of check_point ========='

	end

c*****************************************************************

