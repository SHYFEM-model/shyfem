!
! $Id: newcon.f,v 1.56 2010-03-22 15:29:31 georg Exp $
!
! routines for concentration
!
! contents :
!
!-------------------------------------------------------------
!
! subroutine scal_adv(what,ivar
!     +                          ,scal,ids
!     +                          ,rkpar,wsink
!     +                          ,difhv,difv,difmol)
!		shell for scalar (for parallel version)

! subroutine scal_adv_nudge(what,ivar
!     +				,scal,ids
!     +				,rkpar,wsink
!     +                          ,difhv,difv,difmol
!     +				,sobs,robs)
!		shell for scalar with nudging (for parallel version)
!
! subroutine scal_adv_fact(what,ivar,fact
!     +                          ,scal,ids
!     +                          ,rkpar,wsink,wsinkv,rload,load
!     +                          ,difhv,difv,difmol)
!		shell for scalar (for parallel version)
!		special version for cohesive sediments with factor
!
!-------------------------------------------------------------
!
! subroutine scal3sh(what,cnv,nlvddi,rcv,cobs,robs,rkpar
!					,wsink,wsinkv,rload,load
!     +                                 ,difhv,difv,difmol)
!		shell for scalar T/D
!
! subroutine conz3d(cn1,co1
!     +                  ,ddt
!     +                  ,rkpar,difhv,difv
!     +                  ,difmol,cbound
!     +                  ,itvd,gradxv,gradyv
!     +                  ,cobs,robs
!     +                  ,wsink,wsinkv
!     +                  ,rload,load
!     +                  ,azpar,adpar,aapar
!     +                  ,istot,isact
!     +                  ,nlvddi,nlv)
!

!-------------------------------------------------------------
!
! subroutine massconc(mode,cn,nlvddi,mass)
!		computes total mass of conc
!
! subroutine check_scal_bounds(cnv,cmin,cmax,eps,bstop)
!		checks if scalar is out of bounds
!
! subroutine assert_min_max_property(cnv,cov,sbconz,rmin,rmax,eps)
!		checks min/max property
!
! subroutine stb_histo(it,nlvddi,nkn,ilhkv,cwrite)
!		writes histogram info about stability index
!
!-------------------------------------------------------------
!
! notes:
!
!	dispersion:
!	  scal_adv(...,cnv,ids,...)
!
!	scal_adv(...,ids,...)
!	  bnds_trans_new(...,ids,r3v,...)	(transfers BC to matrix)
!         call scal3sh(...,cnv,r3v,...)
!	
!	scal3sh(...,cnv,r3v,...)
!         call make_scal_flux(...,cnv,r3v,sbconz,...)
!	  do
!           call conz3d(...,cnv,sbconz,...)
!           call assert_min_max_property(...,cnv,sbconz,...)
!           call bndo_setbc(it,what,nlvddi,cnv,rcv,uprv,vprv)
!	  end do
!
!-------------------------------------------------------------
!
! revision log :
!
! 14.08.1998	ggu	rkpar/rvpar -> chpar/cvpar
! 14.08.1998	ggu	use ilhkv to scan vertical levels on node
! 14.08.1998	ggu	$$LEV0 - bug fix : vertical level 0 used
! 19.08.1998    ggu     call to conzfi changed
! 26.08.1998    ggu     cleaned up conzsh
! 26.08.1998    ggu     conz uses zeov,zenv for water level
! 28.10.1999    ggu     names changed
! 07.03.2000    ggu     constant vertical eddy coefficient subst. with difv
! 20.06.2000    ggu     pass difmol to conz3d and use it
! 05.12.2001    ggu     variable horizontal diffusion, limit on dif.coef.
! 11.10.2002    ggu     file cleaned, t/shdif are set equal
! 14.10.2002    ggu     rstot re-introduced as rstol
! 09.09.2003    ggu     call to scal3sh changed -> new arg nlvddi
! 10.03.2004    ggu     call conwrite() to write stability param to nos file
! 13.03.2004    ggu     new boundary conditions through flux (cbound)
! 15.10.2004    ggu     boundary conditions back to old
! 17.01.2005    ggu     new routines with difhv
! 17.01.2005    ggu     get_stability and get_stab_index in this file
! 03.03.2005    ggu     new 3d boundary arrays implemented
! 16.08.2005    ggu     TVD algorithm implemented (gradxv,gradyv,grad_tvd,btvd)
! 04.11.2005    ggu     TVD changes from andrea integrated
! 07.11.2005    ggu     parameter itvd introduced for TVD
! 07.11.2005    ggu     sinking velocity wsink introduced in call to scal3sh
! 11.11.2005    ggu     bug fix in grad_tvd (ggx/ggy in layer loop now)
! 11.11.2005    ggu     new routine grad_2d()
! 16.02.2006    ggu     set w to zero at surface and bottom (WZERO)
! 23.03.2006    ggu     changed time step to double precision
! 08.08.2007    ggu     new parameter istot_max
! 23.08.2007    ggu     test for boundary nodes using routines in testbndo.h
! 18.09.2007    ggu     new subroutine check_scal
! 01.10.2007    ggu     Hack for ssurface -> set to 0 or -999 (temp)
! 17.03.2008    ggu     new open boundary routines introduced
! 08.04.2008    ggu     treatment of boundaries changed
! 22.04.2008    ggu     parallelization: scal_adv, scal_bnd
! 22.04.2008    ggu     local saux, sbflux, no explh, cl{c|m|p}e
! 22.04.2008    ggu     new routine scal_adv_fact for cohesive sediments
! 23.04.2008    ggu     call to bnds_set_def() changed
! 28.04.2008    ggu     rstol deleted
! 28.04.2008    ggu     new routines for stability, s/getistot deleted
! 28.04.2008    ggu     conz3sh into own file
! 24.06.2008    ggu     rstol re-introduced
! 08.11.2008    ggu     BUGFIX in conz3d (vertical velocity)
! 19.11.2008    ggu     changes in advect_stability() - incomplete
! 27.01.2009    aac     bugs in TVD scheme fixed
! 24.03.2009    ggu     more bugs in TVD scheme fixed
! 31.03.2009    ggu     TVD algorithm tested and cleaned
! 20.04.2009    ggu     test for parallel execution (parallel_test)
! 13.10.2009    ggu     write_elem_info() substituted with check_elem()
! 12.11.2009    ggu     make_scal_flux() into internal time loop for stability
! 12.11.2009    ggu     in conz_stab loop over all k that are not z-boundaries
! 16.02.2010    ggu     use wdiff also in stab, use point sources in stab
! 16.02.2010    ggu     min/max property, sbconz is passed to conz3d
! 19.02.2010    ggu     restructured, stab routines into newstab.f
! 10.03.2010    ggu     in assert_min_max_property() check all nodes (also BC)
! 11.03.2010    ggu     in assert_min_max_property() do not check ibtyp=1
! 12.03.2010    ggu     in assert_min_max_property() limit error messages
! 22.03.2010    ggu     bug fix for evaporation (distr. sources) BUG_2010_01
! 15.12.2010    ggu     new routine vertical_flux_ie() for vertical tvd
! 26.01.2011    ggu     nudging implemented (scal_adv_nudge, cobs, robs)
! 16.02.2011    ggu     pass robs to info_stability()
! 23.03.2011    ggu     new parameter itvdv
! 25.03.2011    ggu     error check for aapar and itvdv
! 01.06.2011    ggu     wsink for stability integrated
! 12.07.2011    ggu     run over nlv, not nlvddi, vertical_flux() for lmax>1
! 15.07.2011    ggu     call vertical_flux() anyway (BUG)
! 21.06.2012    ggu&ccf variable vertical sinking velocity integrated
! 03.12.2013    ggu&deb bug fix for horizontal diffusion
! 15.05.2014    ggu     write min/max error only for levdbg >= 3
! 10.07.2014    ggu     only new file format allowed
! 20.10.2014    ggu     accept ids from calling routines
! 22.10.2014    ccf     load in call to scal3sh
! 20.05.2015    ggu     accumulate over nodes (for parallel version)
! 30.09.2015    ggu     routine cleaned, no reals in conz3d
! 26.10.2015    ggu     critical omp sections introduced (eliminated data race)
! 26.10.2015    ggu     mass check only for levdbg > 2
!
!*********************************************************************
!-----------------------------------------------------------------------------------
        module concentration
!-----------------------------------------------------------------------------------
        contains
!-----------------------------------------------------------------------------------

	subroutine scal_adv(what,ivar,scal,ids,rkpar,wsink,difhv,difv,difmol)

! shell for scalar (for parallel version)

	use levels, only : nlvdi,nlvmax
	use basin, only : nkn,nel,ngr,mbw
        use utility
        use bnd_scalar

	implicit none

        character*(*) what
	integer ivar
        double precision scal(nlvdi,nkn)
        integer ids(*)
        double precision rkpar
	double precision wsink
        double precision difhv(nlvdi,nel)
	double precision difv(0:nlvdi,nkn)
        double precision difmol

	include 'femtime.h'

	!include 'const_aux.h'

        double precision r3v(nlvdi,nkn)
        double precision caux(0:nlvdi,nkn)
	double precision load(nlvdi,nkn)	!load [kg/s]
	double precision cobs(nlvdi,nkn)	!load [kg/s]

	double precision dtime
	double precision robs,rload
        integer iwhat
	character*10 whatvar,whataux

	robs = 0.
	rload = 0.
	caux = 1.
	cobs = 1.
	load = 1.

!--------------------------------------------------------------
! make identifier for variable
!--------------------------------------------------------------

!$OMP CRITICAL
	whatvar = what
	if( ivar .ne. 0 ) then
          write(whataux,'(i2)') ivar
          whatvar = what // '_' // whataux
	end if
        iwhat = ichanm(whatvar)
!$OMP END CRITICAL

!--------------------------------------------------------------
! transfer boundary conditions of var ivar to 3d matrix r3v
!--------------------------------------------------------------

	dtime = it

	call bnds_trans_new(whatvar(1:iwhat),ids,dtime,ivar,nkn,nlvmax,nlvdi,r3v)

!--------------------------------------------------------------
! do advection and diffusion
!--------------------------------------------------------------

        call scal3sh(whatvar(1:iwhat),scal,nlvdi,r3v,cobs,robs  &
     &			,rkpar,wsink,caux,rload,load,difhv,difv,difmol)

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*********************************************************************

	subroutine scal_adv_nudge(what,ivar,scal,ids,rkpar,wsink,difhv,difv,difmol,sobs,robs)

! shell for scalar with nudging (for parallel version)

	use levels
	use basin, only : nkn,nel,ngr,mbw
        use shympi
        use utility
        use bnd_scalar

	implicit none

        character*(*) what
	integer ivar
        double precision scal(nlvdi,nkn_local)
        integer ids(*)
        double precision rkpar
	double precision wsink
        double precision difmol
	double precision robs
#ifdef DEBUGON
        double precision sobs(nlvdi,nkn_local)		!observations
        double precision caux(0:nlvdi,nkn_local)
	double precision difv(0:nlvdi,nkn_local)
        double precision difhv(nlvdi,nel_local)
#else
        double precision sobs(nlvdi,nkn)		!observations
        double precision caux(0:nlvdi,nkn)
	double precision difv(0:nlvdi,nkn)
        double precision difhv(nlvdi,nel)
#endif

	include 'femtime.h'

	!include 'const_aux.h'

        double precision r3v(nlvdi,nkn)
	double precision load(nlvdi,nkn)	!load [kg/s]

	integer ierr,l,k,lmax
	double precision eps
	double precision rload
	double precision dtime
        integer iwhat
	character*10 whatvar,whataux

	rload = 0.
	caux = 1.
	load = 1.

!--------------------------------------------------------------
! make identifier for variable
!--------------------------------------------------------------

	whatvar = what
	if( ivar .ne. 0 ) then
          write(whataux,'(i2)') ivar
          whatvar = what // whataux
	end if
        iwhat = ichanm(whatvar)

!--------------------------------------------------------------
! transfer boundary conditions of var ivar to 3d matrix r3v
!--------------------------------------------------------------

	dtime = it

	call bnds_trans_new(whatvar(1:iwhat),ids,dtime,ivar,nkn,nlvmax,nlvdi,r3v)

!--------------------------------------------------------------
! do advection and diffusion
!--------------------------------------------------------------

        call scal3sh(whatvar(1:iwhat),scal,nlvdi,r3v,sobs,robs  &
     &			,rkpar,wsink,caux,rload,load,difhv,difv,difmol)

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*********************************************************************

	subroutine scal_adv_fact(what,ivar,fact,scal,ids        &
     &		,rkpar,wsink,wsinkv,rload,load,difhv,difv,difmol)

! shell for scalar (for parallel version)
!
! special version with factor for BC, variable sinking velocity and loads

	use levels, only : nlvdi,nlvmax
	use basin, only : nkn,nel,ngr,mbw
        use utility
        use bnd_routines
        use bnd_scalar

	implicit none

        character*(*) what
	integer ivar
	double precision fact			!factor for boundary condition
        double precision scal(nlvdi,nkn)
        integer ids(*)
        double precision rkpar
	double precision wsink
	double precision wsinkv(0:nlvdi,nkn)
	double precision rload			!load factor (1 for load given)
	double precision load(nlvdi,nkn)		!load [kg/s]
        double precision difhv(nlvdi,nel)
	double precision difv(0:nlvdi,nkn)
        double precision difmol

	include 'femtime.h'

        double precision r3v(nlvdi,nkn)

	double precision dtime
	double precision robs
        integer iwhat
	character*20 whatvar,whataux

	robs = 0.

!--------------------------------------------------------------
! make identifier for variable
!--------------------------------------------------------------

	whatvar = what
	if( ivar .ne. 0 ) then
          write(whataux,'(i2)') ivar
          whatvar = what // whataux
	end if
        iwhat = ichanm(whatvar)

!--------------------------------------------------------------
! transfer boundary conditions of var ivar to 3d matrix r3v
!--------------------------------------------------------------

	dtime = it

	call bnds_trans_new(whatvar(1:iwhat),ids,dtime,ivar,nkn,nlvmax,nlvdi,r3v)

!--------------------------------------------------------------
! multiply boundary condition with factor
!--------------------------------------------------------------

	if( fact .ne. 1. ) then
	  call mult_scal_bc(r3v,fact)
	end if

!--------------------------------------------------------------
! do advection and diffusion
!--------------------------------------------------------------

        call scal3sh(whatvar(1:iwhat),scal,nlvdi,r3v,scal,robs  &
     &		,rkpar,wsink,wsinkv,rload,load,difhv,difv,difmol)

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*********************************************************************

	subroutine scal_bnd(what,t,bnd3)

! sets boundary conditions for scalar - not used anymore - to be deleted

	implicit none

        character*(*) what
	double precision t
	double precision bnd3(1,1)

	stop 'error stop: call to scal_bnd not allowed'

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine scal3sh(what,cnv,nlvddi,rcv,cobs,robs,rkpar,wsink,wsinkv,rload,load  &
     &					,difhv,difv,difmol)

! shell for scalar T/D

	use hydro_print
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
        use shympi
        use para
        use bndo_admin
        use tvd_admin
        use defnames
        use bnd_routines
        use concentration_omp
        use stability
        use time_util
        use timing

	implicit none

! arguments
        character*(*) what
        double precision cnv(nlvddi,nkn_local)
	integer nlvddi		!vertical dimension
        double precision rcv(nlvddi,nkn)	!boundary condition (value of scalar)
	double precision robs		!use nudging
        double precision rkpar
	double precision wsink
	double precision rload
	double precision load(nlvddi,nkn)		!load [kg/s]
        double precision difmol
! parameters
	integer istot_max
	!parameter ( istot_max = 100 )
	!parameter ( istot_max = 200 )
	!parameter ( istot_max = 300 )
	parameter ( istot_max = 1000 )
! common
	include 'femtime.h'
	include 'mkonst.h'

! local
        double precision saux(nlvddi,nkn)		!aux array
        double precision sbflux(nlvddi,nkn)		!flux boundary conditions
        double precision sbconz(nlvddi,nkn)		!conz boundary conditions

	logical btvd,btvd1
	integer isact
	integer istot
	integer itvd
	integer itvdv
	integer iuinfo
	integer levdbg
        double precision dt
	double precision eps
        double precision sindex
	double precision mass,massold,massdiff
	double precision azpar,adpar,aapar
	double precision ssurface,time1
#ifdef DEBUGON
        double precision cobs(nlvddi,nkn_local)	!observations (for nudging)
	double precision gradxv(nlvddi,nkn_local)		!gradient in x for tvd
	double precision gradyv(nlvddi,nkn_local)		!gradient in y for tvd
	double precision wsinkv(0:nlvddi,nkn_local)
	double precision difv(0:nlvddi,nkn_local)
        double precision difhv(nlvddi,nel_local)
#else
        double precision cobs(nlvddi,nkn)	!observations (for nudging)
	double precision gradxv(nlvddi,nkn)		!gradient in x for tvd
	double precision gradyv(nlvddi,nkn)		!gradient in y for tvd
	double precision wsinkv(0:nlvddi,nkn)
	double precision difv(0:nlvddi,nkn)
        double precision difhv(nlvddi,nel)
#endif

        gradxv=0.
        gradyv=0.

!-------------------------------------------------------------
! start of routine
!-------------------------------------------------------------

!-------------------------------------------------------------
! initialization
!-------------------------------------------------------------

	call getaz(azpar)
	adpar=getpar('adpar')
	aapar=getpar('aapar')
	itvd=nint(getpar('itvd'))	!horizontal tvd scheme
	itvdv=nint(getpar('itvdv'))	!vertical tvd scheme
	levdbg = nint(getpar('levdbg'))
	btvd = itvd .gt. 0
	btvd1 = itvd .eq. 1

!$OMP CRITICAL
        call getinfo(iuinfo)  !unit number of info file
!$OMP END CRITICAL

	eps = 1.d-5
	eps = 1.d-4
	eps = 1.d-2

!-------------------------------------------------------------
! check stability criterion -> set istot
!-------------------------------------------------------------

	call get_timestep(dt)

	saux = 0.
	call make_stability(dt,robs,wsink,wsinkv,rkpar,sindex,istot,saux)

!$OMP CRITICAL
        if(shympi_is_master()) then
          write(iuinfo,*) 'stability_',what,':',it,sindex,istot
	end if
!$OMP END CRITICAL

        if( istot .gt. istot_max ) then
	    call info_stability(dt,robs,wsink,wsinkv,rkpar,sindex,istot,saux)
            write(6,*) 'istot  = ',istot,'   sindex = ',sindex
            stop 'error stop scal3sh: istot index too high'
        end if

!-------------------------------------------------------------
! set up TVD scheme
!-------------------------------------------------------------
        if(.not.(shympi_partition_on_elements())) then
	  call tvd_init(itvd)
        end if

!-------------------------------------------------------------
! set up flux boundary conditions (temporary) -> put in sbflux
!-------------------------------------------------------------

	ssurface = 0.d0	!FIXME - HACK
	if( what .eq. 'temp' ) ssurface = -999.d0

! make_scal_flux has been moved inside time loop of isact (see below)
! to increase stability when treating outflow flux boundary
! this is needed only for discharge < 0 in order to use always the
! ambient tracer concentration

!-------------------------------------------------------------
! transport and diffusion
!-------------------------------------------------------------

	call massconc(-1,cnv,nlvddi,massold)

	do isact=1,istot

	  call make_scal_flux(what,rcv,cnv,sbflux,sbconz,ssurface)
	  !call check_scal_flux(what,cnv,sbconz)

	  if( btvd1 ) call tvd_grad_3d(cnv,gradxv,gradyv,saux,nlvddi)

          call conz3d_omp(cnv,saux,dt,rkpar,difhv,difv,difmol,sbconz,itvd,itvdv,gradxv,gradyv   &
     &		,cobs,robs,wsink,wsinkv,rload,load,azpar,adpar,aapar,istot,isact,nlvddi,nlv)

          if(ln_timing) time1 = shympi_wtime()
          call shympi_exchange_halo_3d_nodes(cnv)
          if(ln_timing) comm_scalar_time = comm_scalar_time + shympi_wtime() - time1


	  call assert_min_max_property(cnv,saux,sbconz,gradxv,gradyv,eps)

          call bndo_setbc(it,what,nlvddi,cnv,rcv,uprv,vprv)

	end do

        !if( shympi_is_parallel() .and. istot > 1 ) then
        !  write(6,*) 'cannot handle istot>1 with mpi yet'
        !  stop 'error stop scal3sh: istot>1'
        !end if
        !call shympi_comment('exchanging scalar: '//trim(what))
        call shympi_exchange_3d_node(cnv)

!-------------------------------------------------------------
! check total mass
!-------------------------------------------------------------

	if( levdbg > 2 ) then
	  call massconc(+1,cnv,nlvddi,mass)
	  massdiff = mass - massold
!$OMP CRITICAL
          if(shympi_is_master())then
	    write(iuinfo,1000) 'scal3sh_',what,':',it,niter,mass,massold,massdiff
	  end if
!$OMP END CRITICAL
	end if

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

        return
 1000   format(a,a,a,2i10,3d13.5)
	end


!**************************************************************

        subroutine conz3d_orig(cn1,co1,ddt,rkpar,difhv,difv,difmol,cbound,itvd,itvdv,gradxv,gradyv      &
     &	        ,cobs,robs,wsink,wsinkv,rload,load,azpar,adpar,aapar,istot,isact,nlvddi,nlev)
!
! computes concentration
!
! cn     new concentration
! co     old concentration              !not used !FIXME
! caux   aux vector
! clow	 lower diagonal of vertical system
! chig	 upper diagonal of vertical system
! ddt    time step
! rkpar  horizontal turbulent diffusivity
! difhv  horizontal turbulent diffusivity (variable between elements)
! difv   vertical turbulent diffusivity
! difmol vertical molecular diffusivity
! cbound boundary condition (mass flux) [kg/s] -> now concentration [kg/m**3]
! itvd	 type of horizontal transport algorithm used
! itvdv	 type of vertical transport algorithm used
! gradxv,gradyv  gradient vectors for TVD algorithm
! cobs	 observations for nudging
! robs	 use observations for nuding (double precision)
! wsink	 factor for settling velocity
! wsinkv variable settling velocity [m/s]
! rload	 factor for loading
! load   load (source or sink) [kg/s]
! azpar  time weighting parameter
! adpar  time weighting parameter for vertical diffusion (ad)
! aapar  time weighting parameter for vertical advection (aa)
! istot	 total inter time steps
! isact	 actual inter time step
! nlvddi	 dimension in z direction
! nlv	 actual needed levels
!
! written 09.01.94 by ggu  (from scratch)
! revised 19.01.94 by ggu  $$flux - flux conserving property
! revised 20.01.94 by ggu  $$iclin - iclin not used to compute volume
! revised 20.01.94 by ggu  $$lumpc - evaluate conz. nodewise
! revised 03.02.94 by ggu  $$itot0 - exception for itot=0 or 3
! revised 04.02.94 by ggu  $$fact3 - factor 3 missing in transport
! revised 04.02.94 by ggu  $$azpar - azpar used to compute transport
! revised 04.02.94 by ggu  $$condry - comute conz also in dry areas
! revised 07.02.94 by ggu  $$istot - istot for fractional time step
! revised 01.06.94 by ggu  restructured for 3-d model
! revised 18.07.94 by ggu  $$htop - use htop instead of htopo for mass cons.
! revised 09.04.96 by ggu  $$rvadj adjust rv in certain areas
!
! solution of purely diffusional part :
!
! dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
!
! C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )
!
! for n-dimensions and
!
! C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )
!
! for 1 dimension
!
! the solution is normalized, i.e.  int(C(x,t)dx) = 1 over the whole area
!
! DPGGU -> introduced double precision to stabilize solution

	use bnd_geom
	use geom
	use depth
	use layer_thickness
	use diff_aux
	use bnd_dynamic
	use area
	use ts
	use hydro_vel
	use hydro_admin
	use evgeom
	use levels
	use basin
	use shympi
        use lnku
        use fem_util
        use tvd_admin

	implicit none
!
! arguments
	integer nlvddi,nlev
        double precision cn1(nlvddi,nkn),co1(nlvddi,nkn)		!DPGGU
        double precision difv(0:nlvddi,nkn)
        double precision difhv(nlvddi,nel)
	double precision difmol
        double precision cbound(nlvddi,nkn)
	integer itvd
	integer itvdv
	double precision gradxv(nlvddi,nkn)
	double precision gradyv(nlvddi,nkn)
	double precision cobs(nlvddi,nkn)
	double precision robs
	double precision wsink
	double precision wsinkv(0:nlvddi,nkn)
	double precision rload
        double precision load(nlvddi,nkn)                      !ccf_load
        double precision ddt,rkpar
        double precision azpar,adpar,aapar			!$$azpar
	integer istot,isact
! common
	include 'femtime.h'
! local
	logical bdebug,bdebug1,debug,btvdv
	integer k,ie,ii,l,iii,ll,ibase
	integer lstart
	integer ilevel
	integer itot,isum	!$$flux
	logical berror
	integer kn(3)
        integer ip(3,3)
        integer n,i,ipp
	integer elems(maxlnk)
        double precision mflux,qflux,cconz
	double precision loading
	double precision wws
	double precision us,vs
	double precision az,azt
	double precision aa,aat,ad,adt
	double precision aj,rk3,rv,aj4,aj12
	double precision hmed,hmbot,hmtop
	double precision hmotop,hmobot,hmntop,hmnbot
	double precision rvptop,rvpbot
	double precision dt,w,aux
	double precision b(3),c(3),f(3)
	double precision rso,rsn,rsot,rsnt,rstot
	double precision hn,ho

	double precision cn(nlvddi,nkn)		!DPGGU	!FIXME
	double precision co(nlvddi,nkn)
	double precision cdiag(nlvddi,nkn)
	double precision clow(nlvddi,nkn)
	double precision chigh(nlvddi,nkn)

	double precision cexpl
	double precision cbm,ccm
	double precision fw(3),fd(3)
	double precision fl(3),fnudge(3)
	double precision flux_tot,flux_tot1,flux_top,flux_bot
        double precision wdiff(3),waux
! local (new)
	double precision clc(nlvddi,3), clm(nlvddi,3), clp(nlvddi,3)
	double precision cle(nlvddi,3)

	double precision cclc(nlvddi,3,nel)
	double precision cclm(nlvddi,3,nel)
	double precision cclp(nlvddi,3,nel)
	double precision ccle(nlvddi,3,nel)

	double precision cl(0:nlvddi+1,3)
	double precision wl(0:nlvddi+1,3)
	double precision vflux(0:nlvddi+1,3)
	double precision cob(0:nlvddi+1,3)
	double precision rtau(0:nlvddi+1,3)

	double precision hdv(0:nlvddi+1)
	double precision haver(0:nlvddi+1)
	double precision hnew(0:nlvddi+1,3)
	double precision hold(0:nlvddi+1,3)
	double precision htnew(0:nlvddi+1,3)
	double precision htold(0:nlvddi+1,3)
	double precision present(0:nlvddi+1)

	double precision cauxn(nlvddi)	!FIXME
	double precision cauxd(nlvddi)
	double precision cauxh(nlvddi)
	double precision cauxl(nlvddi)
! tvd
	logical btvd,bgradup
	integer ic,kc,id,kd,ippp
	integer ies
	integer iext
	double precision fls(3)

! functions
!	integer ipint,ieint

        if(nlv.ne.nlev) stop 'error stop conz3d_orig: level'

!----------------------------------------------------------------
! initialize variables and parameters
!----------------------------------------------------------------

        bdebug1 = .true.
        bdebug1 = .false.
        debug = .false.
        debug = .true.
	bdebug=.false.
	berror=.false.
	!btvdv =.false.
	!btvdv =.true.		!use vertical tvd

	btvd = itvd .gt. 0
	bgradup = itvd .eq. 2	!use upwind gradient for tvd scheme
	btvdv = itvdv .gt. 0

	az=azpar		!$$azpar
	azt=1.-az
	ad=adpar
	adt=1.-ad
	aa=aapar
	aat=1.-aa

	rstot=istot		!$$istot
	rso=(isact-1)/rstot
	rsn=(isact)/rstot
	rsot=1.-rso
	rsnt=1.-rsn

	wws = wsink
	wws = 0.

	dt=ddt/rstot

	if( btvdv .and. aa .ne. 0. ) then
	  write(6,*) 'aapar = ',aapar,'  itvdv = ',itvdv
	  write(6,*) 'Cannot use implicit vertical advection'
	  write(6,*) 'together with vertical TVD scheme.'
	  write(6,*) 'Please set either aapar = 0 (explicit) or'
	  write(6,*) 'itvdv = 0 (no vertical TVD) in the STR file.'
	  stop 'error stop conz3d: vertical tvd scheme'
	end if

!	----------------------------------------------------------------
!	global arrays for accumulation of implicit terms
!	----------------------------------------------------------------

	do k=1,nkn
          do l=1,nlv
	    co1(l,k)=cn1(l,k)	!COLD
	    co(l,k)=cn1(l,k)	!DPGGU
            cn(l,k)=0.
            cdiag(l,k)=0.
            clow(l,k)=0.
            chigh(l,k)=0.
          end do
	end do

!	----------------------------------------------------------------
!	aux elements inside element
!	----------------------------------------------------------------

!	these are aux arrays (bigger than needed) to avoid checking for
!	what layer we are in -> we never get out of bounds

        do l=0,nlv+1
	  hdv(l) = 0.		!layer thickness
          haver(l) = 0.
	  present(l) = 0.	!1. if layer is present
	  do ii=1,3
	    hnew(l,ii) = 0.	!as hreal but with zeta_new
	    hold(l,ii) = 0.	!as hreal but with zeta_old
	    cl(l,ii) = 0.	!concentration in layer
	    wl(l,ii) = 0.	!vertical velocity
	    vflux(l,ii) = 0.	!vertical velocity
	  end do
	end do

!	these are the local arrays for accumulation of implicit terms
!	(maybe we do not need them, but just to be sure...)
!	after accumulation we copy them onto the global arrays

        do l=1,nlv
	  do ii=1,3
	    cle(l,ii) = 0.
	    clc(l,ii) = 0.
	    clm(l,ii) = 0.
	    clp(l,ii) = 0.
	  end do
	end do

!	----------------------------------------------------------------
!	define vertical velocities
!	----------------------------------------------------------------

!----------------------------------------------------------------
! loop over elements
!----------------------------------------------------------------

        do ie=1,nel

	do ii=1,3
          k=nen3v(ii,ie)
	  kn(ii)=k
	  b(ii)=ev(ii+3,ie)
	  c(ii)=ev(ii+6,ie)
	end do

	aj=ev(10,ie)    !area of triangle / 12
	aj4=4.*aj
	aj12=12.*aj
        ilevel=ilhv(ie)

!	----------------------------------------------------------------
!	set up vectors for use in assembling contributions
!	----------------------------------------------------------------

        do l=1,ilevel
	  hdv(l) = hdeov(l,ie)		!use old time step -> FIXME
          !haver(l) = 0.5 * ( hdeov(l,ie) + hdenv(l,ie) )
          haver(l) = rso*hdenv(l,ie) + rsot*hdeov(l,ie)
	  present(l) = 1.
	  do ii=1,3
	    k=kn(ii)
	    hn = hdknv(l,k)		! there are never more layers in ie
	    ho = hdkov(l,k)		! ... than in k
            htold(l,ii) = ho
            htnew(l,ii) = hn
	    hold(l,ii) = rso * hn + rsot * ho
	    hnew(l,ii) = rsn * hn + rsnt * ho
	    cl(l,ii) = co(l,k)
	    cob(l,ii) = cobs(l,k)	!observations
	    rtau(l,ii) = rtauv(l,k)	!observations
	    wl(l,ii) = wlnv(l,k) - wsink * wsinkv(l,k)
	  end do
	end do

	do l=ilevel+1,nlv
	  present(l) = 0.
	end do

!	----------------------------------------------------------------
!	set vertical velocities in surface and bottom layer
!	----------------------------------------------------------------

!	we do not set wl(0,ii) because otherwise we loose concentration
!	through surface
!
!	we set wl(ilevel,ii) to 0 because we are on the bottom
!	and there should be no contribution from this element
!	to the vertical velocity

	do ii=1,3
	  wl(ilevel,ii) = 0.
	end do

!	----------------------------------------------------------------
!	compute vertical fluxes (w/o vertical TVD scheme)
!	----------------------------------------------------------------

	call vertical_flux_ie(btvdv,ie,ilevel,dt,wws,cl,wl,hold,vflux)

!----------------------------------------------------------------
! loop over levels
!----------------------------------------------------------------

        do l=1,ilevel

        us=az*utlnv(l,ie)+azt*utlov(l,ie)             !$$azpar
        vs=az*vtlnv(l,ie)+azt*vtlov(l,ie)

        rk3 = 3. * rkpar * difhv(l,ie)

	cbm=0.
	ccm=0.
	itot=0
	isum=0
	do ii=1,3
	  k=kn(ii)
	  f(ii)=us*b(ii)+vs*c(ii)	!$$azpar
	  if(f(ii).lt.0.) then	!flux out of node
	    itot=itot+1
	    isum=isum+ii
	  end if
	  cbm=cbm+b(ii)*cl(l,ii)
	  ccm=ccm+c(ii)*cl(l,ii)

!	  ----------------------------------------------------------------
!	  initialization to be sure we are in a clean state
!	  ----------------------------------------------------------------

	  fw(ii) = 0.
	  cle(l,ii) = 0.
	  clc(l,ii) = 0.
	  clm(l,ii) = 0.
	  clp(l,ii) = 0.

!	  ----------------------------------------------------------------
!	  contributions from horizontal diffusion
!	  ----------------------------------------------------------------

          waux = 0.
          do iii=1,3
            waux = waux + wdifhv(iii,ii,ie) * cl(l,iii)
          end do
          wdiff(ii) = waux

!	  ----------------------------------------------------------------
!	  contributions from vertical diffusion
!	  ----------------------------------------------------------------

!	  in fd(ii) is explicit contribution
!	  the sign is for the term on the left side, therefore
!	  fd(ii) must be subtracted from the right side
!
!	  maybe we should use double precision layer thickness, or even the
!	  time dependent layer thickness

	  rvptop = difv(l-1,k) + difmol
	  rvpbot = difv(l,k) + difmol
	  !hmtop = 2. * rvptop * present(l-1) / (hdv(l-1)+hdv(l))
	  !hmbot = 2. * rvpbot * present(l+1) / (hdv(l)+hdv(l+1))
	  hmotop =2.*rvptop*present(l-1)/(hold(l-1,ii)+hold(l,ii))
	  hmobot =2.*rvpbot*present(l+1)/(hold(l,ii)+hold(l+1,ii))
	  hmntop =2.*rvptop*present(l-1)/(hnew(l-1,ii)+hnew(l,ii))
	  hmnbot =2.*rvpbot*present(l+1)/(hnew(l,ii)+hnew(l+1,ii))

	  fd(ii) = adt * ((cl(l,ii)-cl(l+1,ii))*hmobot -(cl(l-1,ii)-cl(l,ii))*hmotop)

	  clc(l,ii) = clc(l,ii) + ad * ( hmntop + hmnbot )
	  clm(l,ii) = clm(l,ii) - ad * ( hmntop )
	  clp(l,ii) = clp(l,ii) - ad * ( hmnbot )

!	  ----------------------------------------------------------------
!	  contributions from vertical advection
!	  ----------------------------------------------------------------

!	  in fw(ii) is explicit contribution
!	  the sign is for the term on the left side, therefore
!	  fw(ii) must be subtracted from the right side
!
!	  if we are in last layer, w(l,ii) is zero
!	  if we are in first layer, w(l-1,ii) is zero (see above)

	  w = wl(l-1,ii) - wws		!top of layer
	  if( l .eq. 1 ) w = 0.		!surface -> no transport (WZERO)
	  if( w .ge. 0. ) then
	    fw(ii) = aat*w*cl(l,ii)
	    flux_top = w*cl(l,ii)
	    clc(l,ii) = clc(l,ii) + aa*w
	  else
	    fw(ii) = aat*w*cl(l-1,ii)
	    flux_top = w*cl(l-1,ii)
	    clm(l,ii) = clm(l,ii) + aa*w
	  end if

	  w = wl(l,ii) - wws		!bottom of layer
	  if( l .eq. ilevel ) w = 0.	!bottom -> handle flux elsewhere (WZERO)
	  if( w .gt. 0. ) then
	    fw(ii) = fw(ii) - aat*w*cl(l+1,ii)
	    flux_bot = w*cl(l+1,ii)
	    clp(l,ii) = clp(l,ii) - aa*w
	  else
	    fw(ii) = fw(ii) - aat*w*cl(l,ii)
	    flux_bot = w*cl(l,ii)
	    clc(l,ii) = clc(l,ii) - aa*w
	  end if

	  flux_tot1 = aat * ( flux_top - flux_bot )
	  flux_tot = aat * ( vflux(l-1,ii) - vflux(l,ii) )

!	  if( .not. btvdv ) then
!	  if( flux_tot .ne. flux_tot1 .or. flux_tot .ne. fw(ii) ) then
!	   !if( ie .eq. 100 ) then
!	    write(6,*) '********** vflux   ',ie,ii,l,ilevel
!	    write(6,*) fw(ii),flux_tot1,flux_tot
!	    write(6,*) flux_top,flux_bot
!	    do ll=0,ilevel
!	    write(6,*) ll,vflux(ll,ii)
!	    end do
!	   !end if
!	  end if
!	  end if

	  fw(ii) = flux_tot
	end do

!	----------------------------------------------------------------
!	contributions from horizontal advection (only explicit)
!	----------------------------------------------------------------
!
!	f(ii) > 0 ==> flux into node ii
!	itot=1 -> flux out of one node
!		compute flux with concentration of this node
!	itot=2 -> flux into one node
!		for flux use conz. of the other two nodes and
!		minus the sum of these nodes for the flux of this node

	if(itot.eq.1) then	!$$flux
	  fl(1)=f(1)*cl(l,isum)
	  fl(2)=f(2)*cl(l,isum)
	  fl(3)=f(3)*cl(l,isum)
	else if(itot.eq.2) then
	  isum=6-isum
	  fl(1)=f(1)*cl(l,1)
	  fl(2)=f(2)*cl(l,2)
	  fl(3)=f(3)*cl(l,3)
	  fl(isum) = 0.
	  fl(isum) = -(fl(1)+fl(2)+fl(3))
	  isum=6-isum		!reset to original value
	else			!exception	$$itot0
	  fl(1)=0.
	  fl(2)=0.
	  fl(3)=0.
	end if

!	----------------------------------------------------------------
!	horizontal TVD scheme start
!	----------------------------------------------------------------

        if( btvd ) then
	  iext = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( is_external_boundary(k) ) iext = iext + 1
	  end do

          if( iext .eq. 0 ) then
	    call tvd_fluxes(ie,l,itot,isum,dt,cl,co1,gradxv,gradyv,f,fl)
	  end if
	end if

!	----------------------------------------------------------------
!	horizontal TVD scheme finish
!	----------------------------------------------------------------

!	----------------------------------------------------------------
!	contributions from nudging
!	----------------------------------------------------------------

	do ii=1,3
	  fnudge(ii) = robs * rtau(l,ii) * ( cob(l,ii) - cl(l,ii) )
	end do

!	----------------------------------------------------------------
!	sum explicit contributions
!	----------------------------------------------------------------

	do ii=1,3
          hmed = haver(l)                    !new ggu   !HACK
	  cexpl = aj4 * ( hold(l,ii)*cl(l,ii)+ dt * (hold(l,ii)*fnudge(ii)      &
     &		+ 3.*fl(ii) - fw(ii) - rk3*hmed*wdiff(ii) - fd(ii)))
	  cle(l,ii) = cle(l,ii) + cexpl
	  !k=kn(ii)
	  !cn(l,k) = cn(l,k) + cexpl
	end do

	end do		! loop over l

!----------------------------------------------------------------
! end of loop over l
!----------------------------------------------------------------

!----------------------------------------------------------------
! set up implicit contributions
!----------------------------------------------------------------

! cdiag is diagonal of tri-diagonal system
! chigh is high (right) part of tri-diagonal system
! clow is low (left) part of tri-diagonal system
!
! clp -> bottom
! clm -> top

	do ii=1,3
	  clm(1,ii) = 0.
	  clp(ilevel,ii) = 0.
	end do

        !do l=1,ilevel
	!  do ii=1,3
	!    k=kn(ii)
	!    cn(l,k)    = cn(l,k)    +            cle(l,ii)
	!    clow(l,k)  = clow(l,k)  + aj4 * dt * clm(l,ii)
	!    chigh(l,k) = chigh(l,k) + aj4 * dt * clp(l,ii)
	!    cdiag(l,k) = cdiag(l,k) + aj4 * dt * clc(l,ii)
	!    cdiag(l,k) = cdiag(l,k) + aj4 * hnew(l,ii)
	!  end do
	!end do

	do ii=1,3
          do l=1,ilevel
	    ccle(l,ii,ie) =            cle(l,ii)
	    cclm(l,ii,ie) = aj4 * dt * clm(l,ii)
	    cclp(l,ii,ie) = aj4 * dt * clp(l,ii)
	    cclc(l,ii,ie) = aj4 * ( dt * clc(l,ii) + hnew(l,ii) )
	  end do
          do l=ilevel+1,nlv
	    ccle(l,ii,ie) = 0.
	    cclm(l,ii,ie) = 0.
	    cclp(l,ii,ie) = 0.
	    cclc(l,ii,ie) = 0.
	  end do
	end do

	end do		! loop over ie

!----------------------------------------------------------------
! end of loop over elements
!----------------------------------------------------------------

! in cdiag, chigh, clow is matrix (implicit part)
! if explicit calculation, chigh=clow=0 and in cdiag is volume of node [m**3]
! in cnv is mass of node [kg]
! for explicit treatment, cnv/cdiag gives new concentration [kg/m**3]

!----------------------------------------------------------------
! accumulate contributions on each node
!----------------------------------------------------------------

! here we could make just one loop over nkn
! this would include this loop, the boundary conditions
! and then  the solution of the vertical system
! in this case we would only need the following arrays:
! cn(l),clow(l),chigh(l),cdiag(l) (one dimensional arrays over the vertical)

	do k=1,nkn
	  call get_elems_around(k,maxlnk,n,elems)
	  ilevel = ilhkv(k)
	  do i=1,n
	    ie = elems(i)
	    ii = ithis(k,ie)
	    if( ii == 0 .or. nen3v(ii,ie) /= k ) then
	      stop 'error stop: cannot find ii...'
	    end if
	    do l=1,ilevel
	      cn(l,k)    = cn(l,k)    + ccle(l,ii,ie)
	      clow(l,k)  = clow(l,k)  + cclm(l,ii,ie)
	      chigh(l,k) = chigh(l,k) + cclp(l,ii,ie)
	      cdiag(l,k) = cdiag(l,k) + cclc(l,ii,ie)
	    end do
	  end do
	end do

        !call shympi_comment('shympi_elem: exchange scalar')
	if( shympi_partition_on_elements() ) then
          call shympi_exchange_and_sum_3d_nodes(cn)
          call shympi_exchange_and_sum_3d_nodes(cdiag)
          call shympi_exchange_and_sum_3d_nodes(clow)
          call shympi_exchange_and_sum_3d_nodes(chigh)
	end if

!----------------------------------------------------------------
! integrate boundary conditions
!----------------------------------------------------------------

! in case of negative flux (qflux<0) must check if node is OBC (BUG_2010_01)

	do k=1,nkn
	  ilevel = ilhkv(k)
	  do l=1,ilevel
            !mflux = cbound(l,k)		!mass flux has been passed
	    cconz = cbound(l,k)		!concentration has been passed
	    qflux = mfluxv(l,k)
	    if( qflux .lt. 0. .and. is_boundary(k) ) cconz = cn1(l,k)
	    mflux = qflux * cconz

            cn(l,k) = cn(l,k) + dt * mflux	!explicit treatment

	    loading = rload*load(l,k)
            if( loading .eq. 0. ) then
	      !nothing
	    else if( loading .gt. 0. ) then    		!treat explicit
              cn(l,k) = cn(l,k) + dt * loading
            else !if( loading .lt. 0. ) then		!treat quasi implicit
	      if( cn1(l,k) > 0. ) then
                cdiag(l,k) = cdiag(l,k) - dt * loading/cn1(l,k)
	      end if
            end if
	  end do
	end do

!----------------------------------------------------------------
! compute concentration for each node (solve system)
!----------------------------------------------------------------

	if( ( aa .eq. 0. .and. ad .eq. 0. ) .or. ( nlv .eq. 1 ) ) then

	if( nlv .gt. 1 ) then
	  write(6,*) 'conz: computing explicitly ',nlv
	end if

	do k=1,nkn
	 ilevel = ilhkv(k)
	 do l=1,ilevel
	  if(cdiag(l,k).ne.0.) then
	    cn(l,k)=cn(l,k)/cdiag(l,k)
	  end if
	 end do
	end do

	else

	do k=1,nkn
	  ilevel = ilhkv(k)
	  aux=1./cdiag(1,k)
	  chigh(1,k)=chigh(1,k)*aux
	  cn(1,k)=cn(1,k)*aux
	  do l=2,ilevel
	    aux=1./(cdiag(l,k)-clow(l,k)*chigh(l-1,k))
	    chigh(l,k)=chigh(l,k)*aux
	    cn(l,k)=(cn(l,k)-clow(l,k)*cn(l-1,k))*aux
	  end do
	  lstart = ilevel-1
	  do l=lstart,1,-1	!$$LEV0 bug 14.08.1998 -> ran to 0
	    cn(l,k)=cn(l,k)-cn(l+1,k)*chigh(l,k)
	  end do
	end do

	end if

	do k=1,nkn		!DPGGU
          do l=1,nlv
	    cn1(l,k)=cn(l,k)
	  end do
	end do

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine massconc(mode,cn,nlvddi,mass)

! computes total mass of conc

	use levels
	use basin, only : nkn,nel,ngr,mbw
        use shympi
        use elems_dealing
        use timing

	implicit none

! arguments
	integer mode
	integer nlvddi
	double precision cn(nlvddi,nkn)
	double precision mass
! common
	include 'femtime.h'
! local
	integer k,l,lmax,ntot
        double precision vol
	double precision sum,masstot,time1

        masstot = 0.

        if(shympi_partition_on_elements()) then
	  ntot = nkn_inner 	!SHYMPI_ELEM - should be total nodes to use
        else
          ntot = nkn
        end if

        do k=1,ntot
	  lmax = ilhkv(k)
          sum = 0.
          do l=1,lmax
            vol = volnode(l,k,mode)
            sum = sum + cn(l,k) * vol
          end do
          masstot = masstot + sum
        end do

	mass = masstot

        if(ln_timing) time1 = shympi_wtime()
        mass = shympi_sum(mass)
        if(ln_timing) comm_scalar_time = comm_scalar_time + shympi_wtime() - time1
        !call shympi_comment('massconc: shympi_sum(masstot)')

!	write(88,*) 'tot mass: ',it,mass

	end

!**************************************************************

        subroutine check_scal_bounds(cnv,cmin,cmax,eps,bstop)

! checks if scalar is out of bounds

	use levels
        use fem_util
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        double precision cnv(nlvdi,nkn)
        double precision cmin,cmax
	double precision eps
	logical bstop		!stop simulation if true


	logical berror
        integer k,l,lmax,kext
        double precision cc,cn,cx,cd

	berror = .false.
	cd = cmax - cmin
	if( cd .le. 0. ) return

	cn = cmin - eps*cd
	cx = cmax + eps*cd

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cc = cnv(l,k)
            if( cc .lt. cn .or. cc .gt. cx ) then
                kext = ipext(k)
                write(6,*) 'scalar out of bounds: ',k,kext,l,lmax,cc
		berror = .true.
            end if
          end do
        end do

	if( bstop .and. berror ) then
	  stop 'error stop check_scal_bounds: out of bounds'
	end if

        end

!*****************************************************************

	subroutine assert_min_max_property(cnv,cov,sbconz,rmin,rmax,eps)

! checks min/max property

	use bnd_geom
	use bnd_dynamic
	use levels
	use basin
        use para
        use fem_util
        use bndo_admin

	implicit none

        double precision cnv(nlvdi,nkn)			!new concentration
        double precision cov(nlvdi,nkn)			!old concentration
        double precision sbconz(nlvdi,nkn)		!conz boundary conditions
	double precision rmin(nlvdi,nkn)		!aux arrray to contain min
	double precision rmax(nlvdi,nkn)		!aux arrray to contain max
	double precision eps

	include 'femtime.h'


	logical bwrite,bstop
	integer k,ie,l,ii,lmax,ierr
	integer levdbg
	double precision amin,amax,c,qflux,dmax
	double precision drmax,diff
	double precision dt

	bwrite = .true.		! write every violation
	bwrite = .false.		! write every violation
	bstop = .false.		! stop after error

	levdbg = nint(getpar('levdbg'))

!---------------------------------------------------------------
! vertical contribution (normally implicit -> whole column)
!---------------------------------------------------------------

	do k=1,nkn
	  lmax = ilhkv(k)
	  amin = +1.e+30
	  amax = -1.e+30
	  do l=1,lmax
	    amin = min(amin,cov(l,k))
	    amax = max(amax,cov(l,k))
	  end do
	  do l=1,lmax
	    rmin(l,k) = amin
	    rmax(l,k) = amax
	  end do
	end do

!---------------------------------------------------------------
! point sources
!---------------------------------------------------------------

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    qflux = mfluxv(l,k)
	    if( qflux .gt. 0. ) then
	      c = sbconz(l,k)
	      rmin(l,k) = min(rmin(l,k),c)
	      rmax(l,k) = max(rmax(l,k),c)
	    end if
	  end do
	end do

!---------------------------------------------------------------
! horizontal contribution
!---------------------------------------------------------------

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    amin = +1.e+30
	    amax = -1.e+30
	    do ii=1,3
	      k = nen3v(ii,ie)
	      amin = min(amin,cov(l,k))
	      amax = max(amax,cov(l,k))
	    end do
	    do ii=1,3
	      k = nen3v(ii,ie)
	      rmin(l,k) = min(rmin(l,k),amin)
	      rmax(l,k) = max(rmax(l,k),amax)
	    end do
	  end do
	end do

!---------------------------------------------------------------
! check with new concentration
!---------------------------------------------------------------

	ierr = 0
	dmax = 0.
	drmax = 0.

	do k=1,nkn
	 !if( .not. is_external_boundary(k) ) then	!might be relaxed
	 if( .not. is_zeta_bound(k) ) then	!might be relaxed
	  lmax = ilhkv(k)
	  do l=1,lmax
	    c = cnv(l,k)
	    !rm1 = rmin(l,k)
	    !rm2 = rmax(l,k)
	    if( c .lt. rmin(l,k) .or. c .gt. rmax(l,k) ) then
	      amin = rmin(l,k)
	      amax = rmax(l,k)
	      diff = max(rmin(l,k)-c,c-rmax(l,k))
	      dmax = max(dmax,diff)
	      drmax = max(drmax,diff*(amax-amin))
	      if( bwrite ) then
	        write(6,*) 'min/max property violated: ',it,idt
	        write(6,*) '   ',l,k,ipext(k)
	        write(6,*) '   ',c,amin,amax
	      end if
	      ierr = ierr + 1
	    end if
	  end do
	 end if
	end do

!---------------------------------------------------------------
! check error condition
!---------------------------------------------------------------

	if( ierr .gt. 0 ) then
	  if( drmax .gt. eps .and. dmax .gt. eps ) then
	    if( levdbg .ge. 3 ) then
	      write(6,*) 'min/max error: ',it,ierr,dmax,drmax
	    end if
	  end if
	  !write(94,*) 'min/max error violated: ',it,ierr,dmax,drmax
	  if( bstop ) then
	    stop 'error stop assert_min_max_property: violation'
	  end if
	end if

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine stb_histo(it,nlvddi,nkn,ilhkv,cwrite)

! writes histogram info about stability index

        use histogram

        implicit none

        integer it
        integer nlvddi,nkn
        integer ilhkv(nkn)
        double precision cwrite(nlvddi,nkn)

        integer ndim
        parameter(ndim=11)

        integer nbin
        double precision aux
        integer k,l,lmax

        integer ic(ndim+1)
        double precision bins(ndim)
        save bins
        data bins /1.,2.,5.,10.,15.,20.,30.,40.,50.,75.,100./

        aux = 0.
        nbin = ndim

        call histo_init(nbin,aux,aux,bins)

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            call histo_insert(cwrite(l,k))
          end do
        end do

        call histo_final(ic)

        write(98,*) it,ic

        end

!*****************************************************************

!-----------------------------------------------------------------------------------
        end module concentration
!-----------------------------------------------------------------------------------
