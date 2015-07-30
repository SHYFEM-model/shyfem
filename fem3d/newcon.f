c
c $Id: newcon.f,v 1.56 2010-03-22 15:29:31 georg Exp $
c
c routines for concentration
c
c contents :
c
c-------------------------------------------------------------
c
c subroutine scal_adv(what,ivar
c     +                          ,scal,ids
c     +                          ,rkpar,wsink
c     +                          ,difhv,difv,difmol)
c		shell for scalar (for parallel version)

c subroutine scal_adv_nudge(what,ivar
c     +				,scal,ids
c     +				,rkpar,wsink
c     +                          ,difhv,difv,difmol
c     +				,sobs,robs)
c		shell for scalar with nudging (for parallel version)
c
c subroutine scal_adv_fact(what,ivar,fact
c     +                          ,scal,ids
c     +                          ,rkpar,wsink,wsinkv,rload,load
c     +                          ,difhv,difv,difmol)
c		shell for scalar (for parallel version)
c		special version for cohesive sediments with factor
c
c-------------------------------------------------------------
c
c subroutine scal3sh(what,cnv,nlvddi,rcv,cobs,robs,rkpar
c					,wsink,wsinkv,rload,load
c     +                                 ,difhv,difv,difmol)
c		shell for scalar T/D
c
c subroutine conz3d(cn1,co1
c     +                  ,ddt
c     +                  ,rkpar,difhv,difv
c     +                  ,difmol,cbound
c     +                  ,itvd,gradxv,gradyv
c     +                  ,cobs,robs
c     +                  ,wsink,wsinkv
c     +                  ,rload,load
c     +                  ,azpar,adpar,aapar
c     +                  ,istot,isact
c     +                  ,nlvddi,nlv)
c
c subroutine conzstab(cn1,co1
c     +                  ,ddt
c     +                  ,robs,wsink,wsinkv
c     +                  ,rkpar,difhv,difv
c     +                  ,difmol,azpar
c     +                  ,adpar,aapar
c     +                  ,sindex
c     +                  ,istot,isact
c     +                  ,nlvddi,nlv)
c
c-------------------------------------------------------------
c
c subroutine massconc(mode,cn,nlvddi,mass)
c		computes total mass of conc
c
c subroutine check_scal_bounds(cnv,cmin,cmax,eps,bstop)
c		checks if scalar is out of bounds
c
c subroutine assert_min_max_property(cnv,cov,sbconz,rmin,rmax,eps)
c		checks min/max property
c
c subroutine stb_histo(it,nlvddi,nkn,ilhkv,cwrite)
c		writes histogram info about stability index
c
c-------------------------------------------------------------
c
c notes:
c
c	dispersion:
c	  scal_adv(...,cnv,ids,...)
c
c	scal_adv(...,ids,...)
c	  bnds_trans_new(...,ids,r3v,...)	(transfers BC to matrix)
c         call scal3sh(...,cnv,r3v,...)
c	
c	scal3sh(...,cnv,r3v,...)
c         call make_scal_flux(...,cnv,r3v,sbconz,...)
c	  do
c           call conz3d(...,cnv,sbconz,...)
c           call assert_min_max_property(...,cnv,sbconz,...)
c           call bndo_setbc(it,what,nlvddi,cnv,rcv,uprv,vprv)
c	  end do
c
c-------------------------------------------------------------
c
c revision log :
c
c 14.08.1998	ggu	rkpar/rvpar -> chpar/cvpar
c 14.08.1998	ggu	use ilhkv to scan vertical levels on node
c 14.08.1998	ggu	$$LEV0 - bug fix : vertical level 0 used
c 19.08.1998    ggu     call to conzfi changed
c 26.08.1998    ggu     cleaned up conzsh
c 26.08.1998    ggu     conz uses zeov,zenv for water level
c 28.10.1999    ggu     names changed
c 07.03.2000    ggu     constant vertical eddy coefficient subst. with difv
c 20.06.2000    ggu     pass difmol to conz3d and use it
c 05.12.2001    ggu     variable horizontal diffusion, limit on dif.coef.
c 11.10.2002    ggu     file cleaned, t/shdif are set equal
c 11.10.2002    ggu     con3sh removed, conzstab better commented
c 14.10.2002    ggu     rstot re-introduced as rstol
c 09.09.2003    ggu     call to scal3sh changed -> new arg nlvddi
c 10.03.2004    ggu     call conwrite() to write stability param to nos file
c 13.03.2004    ggu     new boundary conditions through flux (cbound)
c 15.10.2004    ggu     boundary conditions back to old
c 02.12.2004    ggu     return also sindex in conzstab
c 17.01.2005    ggu     new routines with difhv
c 17.01.2005    ggu     get_stability and get_stab_index in this file
c 03.03.2005    ggu     new 3d boundary arrays implemented
c 16.08.2005    ggu     TVD algorithm implemented (gradxv,gradyv,grad_tvd,btvd)
c 04.11.2005    ggu     TVD changes from andrea integrated
c 07.11.2005    ggu     parameter itvd introduced for TVD
c 07.11.2005    ggu     sinking velocity wsink introduced in call to scal3sh
c 11.11.2005    ggu     bug fix in grad_tvd (ggx/ggy in layer loop now)
c 11.11.2005    ggu     new routine grad_2d()
c 16.02.2006    ggu     set w to zero at surface and bottom (WZERO)
c 23.03.2006    ggu     changed time step to real
c 08.08.2007    ggu     new parameter istot_max
c 23.08.2007    ggu     test for boundary nodes using routines in testbndo.h
c 18.09.2007    ggu     new subroutine check_scal
c 01.10.2007    ggu     Hack for ssurface -> set to 0 or -999 (temp)
c 17.03.2008    ggu     new open boundary routines introduced
c 08.04.2008    ggu     treatment of boundaries changed
c 22.04.2008    ggu     parallelization: scal_adv, scal_bnd
c 22.04.2008    ggu     local saux, sbflux, no explh, cl{c|m|p}e
c 22.04.2008    ggu     new routine scal_adv_fact for cohesive sediments
c 23.04.2008    ggu     call to bnds_set_def() changed
c 28.04.2008    ggu     rstol deleted
c 28.04.2008    ggu     new routines for stability, s/getistot deleted
c 28.04.2008    ggu     conz3sh into own file
c 24.06.2008    ggu     rstol re-introduced
c 08.11.2008    ggu     BUGFIX in conz3d (vertical velocity)
c 11.11.2008    ggu     conzstab cleaned
c 19.11.2008    ggu     changes in advect_stability() - incomplete
c 06.12.2008    ggu     in conzstab changed wprv, new routine write_elem_info()
c 27.01.2009    aac     bugs in TVD scheme fixed
c 24.03.2009    ggu     more bugs in TVD scheme fixed
c 31.03.2009    ggu     TVD algorithm tested and cleaned
c 20.04.2009    ggu     test for parallel execution (parallel_test)
c 13.10.2009    ggu     write_elem_info() substituted with check_elem()
c 12.11.2009    ggu     make_scal_flux() into internal time loop for stability
c 12.11.2009    ggu     in conz_stab loop over all k that are not z-boundaries
c 16.02.2010    ggu     use wdiff also in stab, use point sources in stab
c 16.02.2010    ggu     min/max property, sbconz is passed to conz3d
c 19.02.2010    ggu     restructured, stab routines into newstab.f
c 10.03.2010    ggu     in assert_min_max_property() check all nodes (also BC)
c 11.03.2010    ggu     in assert_min_max_property() do not check ibtyp=1
c 12.03.2010    ggu     in assert_min_max_property() limit error messages
c 22.03.2010    ggu     bug fix for evaporation (distr. sources) BUG_2010_01
c 15.12.2010    ggu     new routine vertical_flux_ie() for vertical tvd
c 26.01.2011    ggu     nudging implemented (scal_adv_nudge, cobs, robs)
c 16.02.2011    ggu     pass robs to info_stability()
c 23.03.2011    ggu     new parameter itvdv
c 25.03.2011    ggu     error check for aapar and itvdv
c 01.06.2011    ggu     wsink for stability integrated
c 12.07.2011    ggu     run over nlv, not nlvddi, vertical_flux() for lmax>1
c 15.07.2011    ggu     call vertical_flux() anyway (BUG)
c 21.06.2012    ggu&ccf variable vertical sinking velocity integrated
c 03.12.2013    ggu&deb bug fix for horizontal diffusion
c 15.05.2014    ggu     write min/max error only for levdbg >= 3
c 10.07.2014    ggu     only new file format allowed
c 20.10.2014    ggu     accept ids from calling routines
c 22.10.2014    ccf     load in call to scal3sh
c 20.05.2015    ggu     accumulate over nodes (for parallel version)
c
c*********************************************************************

	subroutine scal_adv(what,ivar
     +				,scal,ids
     +				,rkpar,wsink
     +                          ,difhv,difv,difmol)

c shell for scalar (for parallel version)

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        include 'param.h'

        character*(*) what
	integer ivar
        real scal(nlvdi,nkn)
        integer ids(*)
        real rkpar
	real wsink
        real difhv(nlvdi,nel)
	real difv(0:nlvdi,nkn)
        real difmol

	include 'femtime.h'

	!include 'const_aux.h'

        real r3v(nlvdi,nkn)
        real caux(0:nlvdi,nkn)
	real load(nlvdi,nkn)	!load [kg/s]
	real cobs(nlvdi,nkn)	!load [kg/s]

	double precision dtime
	real robs,rload
        integer iwhat,ichanm
	character*10 whatvar,whataux

	robs = 0.
	rload = 0.
	caux = 1.
	cobs = 1.
	load = 1.

c--------------------------------------------------------------
c make identifier for variable
c--------------------------------------------------------------

	whatvar = what
	if( ivar .ne. 0 ) then
          write(whataux,'(i2)') ivar
          whatvar = what // '_' // whataux
	end if
        iwhat = ichanm(whatvar)

c--------------------------------------------------------------
c transfer boundary conditions of var ivar to 3d matrix r3v
c--------------------------------------------------------------

	dtime = it

	call bnds_trans_new(whatvar(1:iwhat)
     +			,ids,dtime,ivar,nkn,nlv,nlvdi,r3v)

c--------------------------------------------------------------
c do advection and diffusion
c--------------------------------------------------------------

        call scal3sh(whatvar(1:iwhat)
     +				,scal,nlvdi
     +                          ,r3v,cobs,robs
     +				,rkpar,wsink,caux,rload,load
     +                          ,difhv,difv,difmol)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*********************************************************************

	subroutine scal_adv_nudge(what,ivar
     +				,scal,ids
     +				,rkpar,wsink
     +                          ,difhv,difv,difmol
     +				,sobs,robs)

c shell for scalar with nudging (for parallel version)

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        include 'param.h'

        character*(*) what
	integer ivar
        real scal(nlvdi,nkn)
        integer ids(*)
        real rkpar
	real wsink
        real difhv(nlvdi,nel)
	real difv(0:nlvdi,nkn)
        real difmol
	real sobs(nlvdi,nkn)		!observations
	real robs

	include 'femtime.h'

	!include 'const_aux.h'

        real r3v(nlvdi,nkn)
        real caux(0:nlvdi,nkn)
	real load(nlvdi,nkn)	!load [kg/s]

	integer ierr,l,k,lmax
	real eps
	real rload
	double precision dtime
        integer iwhat,ichanm
	character*10 whatvar,whataux

	rload = 0.
	caux = 1.
	load = 1.

c--------------------------------------------------------------
c make identifier for variable
c--------------------------------------------------------------

	whatvar = what
	if( ivar .ne. 0 ) then
          write(whataux,'(i2)') ivar
          whatvar = what // whataux
	end if
        iwhat = ichanm(whatvar)

c--------------------------------------------------------------
c transfer boundary conditions of var ivar to 3d matrix r3v
c--------------------------------------------------------------

	dtime = it

	call bnds_trans_new(whatvar(1:iwhat)
     +			,ids,dtime,ivar,nkn,nlv,nlvdi,r3v)

c--------------------------------------------------------------
c do advection and diffusion
c--------------------------------------------------------------

        call scal3sh(whatvar(1:iwhat)
     +				,scal,nlvdi
     +                          ,r3v,sobs,robs
     +				,rkpar,wsink,caux,rload,load
     +                          ,difhv,difv,difmol)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*********************************************************************

	subroutine scal_adv_fact(what,ivar,fact
     +				,scal,ids
     +				,rkpar,wsink,wsinkv,rload,load
     +                          ,difhv,difv,difmol)

c shell for scalar (for parallel version)
c
c special version with factor for BC, variable sinking velocity and loads

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

        include 'param.h'

        character*(*) what
	integer ivar
	real fact			!factor for boundary condition
        real scal(nlvdi,nkn)
        integer ids(*)
        real rkpar
	real wsink
	real wsinkv(0:nlvdi,nkn)
	real rload			!load factor (1 for load given)
	real load(nlvdi,nkn)		!load [kg/s]
        real difhv(nlvdi,nel)
	real difv(0:nlvdi,nkn)
        real difmol

	include 'femtime.h'

        real r3v(nlvdi,nkn)

	double precision dtime
	real robs
        integer iwhat,ichanm
	character*20 whatvar,whataux

	robs = 0.

c--------------------------------------------------------------
c make identifier for variable
c--------------------------------------------------------------

	whatvar = what
	if( ivar .ne. 0 ) then
          write(whataux,'(i2)') ivar
          whatvar = what // whataux
	end if
        iwhat = ichanm(whatvar)

c--------------------------------------------------------------
c transfer boundary conditions of var ivar to 3d matrix r3v
c--------------------------------------------------------------

	dtime = it

	call bnds_trans_new(whatvar(1:iwhat)
     +			,ids,dtime,ivar,nkn,nlv,nlvdi,r3v)

c--------------------------------------------------------------
c multiply boundary condition with factor
c--------------------------------------------------------------

	if( fact .ne. 1. ) then
	  call mult_scal_bc(r3v,fact)
	end if

c--------------------------------------------------------------
c do advection and diffusion
c--------------------------------------------------------------

        call scal3sh(whatvar(1:iwhat)
     +				,scal,nlvdi
     +                          ,r3v,scal,robs
     +				,rkpar,wsink,wsinkv,rload,load
     +                          ,difhv,difv,difmol)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*********************************************************************

	subroutine scal_bnd(what,t,bnd3)

c sets boundary conditions for scalar - not used anymore - to be deleted

	implicit none

        include 'param.h'

        character*(*) what
	real t
	real bnd3(1,1)

	stop 'error stop: call to scal_bnd not allowed'

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************
c*********************************************************************
c*********************************************************************

	subroutine scal3sh(what,cnv,nlvddi,rcv,cobs,robs,rkpar
     +					,wsink,wsinkv,rload,load
     +					,difhv,difv,difmol)

c shell for scalar T/D

	use mod_hydro_print
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c parameter
        include 'param.h'
c arguments
        character*(*) what
        real cnv(nlvddi,nkn)
	integer nlvddi		!vertical dimension
        real rcv(nlvddi,nkn)	!boundary condition (value of scalar)
	real cobs(nlvddi,nkn)	!observations (for nudging)
	real robs		!use nudging
        real rkpar
	real wsink
	real wsinkv(0:nlvddi,nkn)
	real rload
	real load(nlvddi,nkn)		!load [kg/s]
        real difhv(nlvddi,nel)
	real difv(0:nlvddi,nkn)
        real difmol
c parameters
	integer istot_max
	!parameter ( istot_max = 100 )
	!parameter ( istot_max = 200 )
	!parameter ( istot_max = 300 )
	parameter ( istot_max = 1000 )
c common
	include 'femtime.h'
	include 'mkonst.h'

c local
        real saux(nlvddi,nkn)		!aux array
        real sbflux(nlvddi,nkn)		!flux boundary conditions
        real sbconz(nlvddi,nkn)		!conz boundary conditions
	real gradxv(nlvddi,nkn)		!gradient in x for tvd
	real gradyv(nlvddi,nkn)		!gradient in y for tvd

	logical btvd,btvd1
	integer isact
	integer istot
	integer itvd
	integer itvdv
	integer iuinfo
        real dt
	real eps
        real sindex
	real mass,massold,massdiff
	real azpar,adpar,aapar
	real ssurface
c function
	real getpar

c-------------------------------------------------------------
c start of routine
c-------------------------------------------------------------

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	call getaz(azpar)
	adpar=getpar('adpar')
	aapar=getpar('aapar')
	itvd=nint(getpar('itvd'))	!horizontal tvd scheme
	itvdv=nint(getpar('itvdv'))	!vertical tvd scheme
	btvd = itvd .gt. 0
	btvd1 = itvd .eq. 1

        call getinfo(iuinfo)  !unit number of info file

	eps = 1.e-5
	eps = 1.e-4
	eps = 1.e-2

c-------------------------------------------------------------
c check stability criterion -> set istot
c-------------------------------------------------------------

	call get_timestep(dt)

	saux = 0.
	call make_stability(dt,robs,wsink,wsinkv,rkpar,sindex,istot,saux)

        write(iuinfo,*) 'stability_',what,':',it,sindex,istot

        if( istot .gt. istot_max ) then
	    call info_stability(dt,robs,wsink,wsinkv,rkpar
     +					,sindex,istot,saux)
            write(6,*) 'istot  = ',istot,'   sindex = ',sindex
            stop 'error stop scal3sh: istot index too high'
        end if

c-------------------------------------------------------------
c set up TVD scheme
c-------------------------------------------------------------

	call tvd_init(itvd)

c-------------------------------------------------------------
c set up flux boundary conditions (temporary) -> put in sbflux
c-------------------------------------------------------------

	ssurface = 0.	!FIXME - HACK
	if( what .eq. 'temp' ) ssurface = -999.

c make_scal_flux has been moved inside time loop of isact (see below)
c to increase stability when treating outflow flux boundary
c this is needed only for discharge < 0 in order to use always the
c ambient tracer concentration

c-------------------------------------------------------------
c transport and diffusion
c-------------------------------------------------------------


	call massconc(-1,cnv,nlvddi,massold)

	do isact=1,istot

	  call make_scal_flux(what,rcv,cnv,sbflux,sbconz,ssurface)
	  !call check_scal_flux(what,cnv,sbconz)

	  if( btvd1 ) call tvd_grad_3d(cnv,gradxv,gradyv,saux,nlvddi)

          call conz3d(
     +           cnv
     +          ,saux
     +          ,dt
     +          ,rkpar,difhv,difv,difmol
     +          ,sbconz
     +		,itvd,itvdv,gradxv,gradyv
     +		,cobs,robs
     +		,wsink,wsinkv
     +		,rload,load
     +		,azpar,adpar,aapar
     +          ,istot,isact
     +          ,nlvddi,nlv
     +               )

	  call assert_min_max_property(cnv,saux,sbconz,gradxv,gradyv,eps)

          call bndo_setbc(it,what,nlvddi,cnv,rcv,uprv,vprv)

	end do

        !if( what .eq. 'salt' ) call check_scal_bounds(cnv,0.,60.,eps,.true.)
        !if( what .eq. 'conz' ) call check_scal_bounds(cnv,0.,100.,eps,.true.)

c-------------------------------------------------------------
c check total mass
c-------------------------------------------------------------

	call massconc(+1,cnv,nlvddi,mass)
	massdiff = mass - massold

	write(iuinfo,1000) 'scal3sh_',what,':'
     +                          ,it,niter,mass,massold,massdiff

	!check_set_unit(6)
	!call check_elem(9914)

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

        return
 1000   format(a,a,a,2i10,3d13.5)
	end


c**************************************************************

        subroutine conz3d(cn1,co1
     +			,ddt
     +                  ,rkpar,difhv,difv
     +			,difmol,cbound
     +		 	,itvd,itvdv,gradxv,gradyv
     +			,cobs,robs
     +			,wsink,wsinkv
     +			,rload,load
     +			,azpar,adpar,aapar
     +			,istot,isact
     +			,nlvddi,nlev)
c
c computes concentration
c
c cn     new concentration
c co     old concentration              !not used !FIXME
c caux   aux vector
c clow	 lower diagonal of vertical system
c chig	 upper diagonal of vertical system
c ddt    time step
c rkpar  horizontal turbulent diffusivity
c difhv  horizontal turbulent diffusivity (variable between elements)
c difv   vertical turbulent diffusivity
c difmol vertical molecular diffusivity
c cbound boundary condition (mass flux) [kg/s] -> now concentration [kg/m**3]
c itvd	 type of horizontal transport algorithm used
c itvdv	 type of vertical transport algorithm used
c gradxv,gradyv  gradient vectors for TVD algorithm
c cobs	 observations for nudging
c robs	 use observations for nuding (real)
c wsink	 factor for settling velocity
c wsinkv variable settling velocity [m/s]
c rload	 factor for loading
c load   load (source or sink) [kg/s]
c azpar  time weighting parameter
c adpar  time weighting parameter for vertical diffusion (ad)
c aapar  time weighting parameter for vertical advection (aa)
c istot	 total inter time steps
c isact	 actual inter time step
c nlvddi	 dimension in z direction
c nlv	 actual needed levels
c
c written 09.01.94 by ggu  (from scratch)
c revised 19.01.94 by ggu  $$flux - flux conserving property
c revised 20.01.94 by ggu  $$iclin - iclin not used to compute volume
c revised 20.01.94 by ggu  $$lumpc - evaluate conz. nodewise
c revised 03.02.94 by ggu  $$itot0 - exception for itot=0 or 3
c revised 04.02.94 by ggu  $$fact3 - factor 3 missing in transport
c revised 04.02.94 by ggu  $$azpar - azpar used to compute transport
c revised 04.02.94 by ggu  $$condry - comute conz also in dry areas
c revised 07.02.94 by ggu  $$istot - istot for fractional time step
c revised 01.06.94 by ggu  restructured for 3-d model
c revised 18.07.94 by ggu  $$htop - use htop instead of htopo for mass cons.
c revised 09.04.96 by ggu  $$rvadj adjust rv in certain areas
c
c solution of purely diffusional part :
c
c dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
c
c C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )
c
c for n-dimensions and
c
c C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )
c
c for 1 dimension
c
c the solution is normalized, i.e.  int(C(x,t)dx) = 1 over the whole area
c
c DPGGU -> introduced double precision to stabilize solution

	use mod_bound_geom
	use mod_geom
	use mod_depth
	use mod_layer_thickness
	use mod_diff_aux
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none
c
c parameters
        include 'param.h'
c arguments
	integer nlvddi,nlev
        real cn1(nlvddi,1),co1(nlvddi,1)		!DPGGU
        real difv(0:nlvddi,1)
        real difhv(nlvddi,1)
	real difmol
        real cbound(nlvddi,1)
	integer itvd
	integer itvdv
	real gradxv(nlvddi,1)
	real gradyv(nlvddi,1)
	real cobs(nlvddi,1)
	real robs
	real wsink
	real wsinkv(0:nlvddi,1)
	real rload
        real load(nlvddi,1)                      !ccf_load
        real ddt,rkpar,azpar,adpar,aapar			!$$azpar
	integer istot,isact
c common
	include 'femtime.h'
	!include 'hydro_print.h'


 
 




c local
	logical bdebug,bdebug1,debug,btvdv
	integer k,ie,ii,l,iii,ll,ibase
	integer lstart
	integer ilevel
	integer itot,isum	!$$flux
	logical berror
	integer kn(3)
        integer ip(3,3)
        integer n,i,ipp
        real rkmin,rkmax
        real mflux,qflux,cconz
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
c local (new)
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
c tvd
	logical btvd,bgradup
	integer ic,kc,id,kd,ippp
	integer ies
	integer iext
	real term,fact
	real conc,cond,conf,conu
	real gcx,gcy,dx,dy
	real u,v
	real rf,psi
	real grad
	double precision fls(3)
        real alfa,dis
        real vel
        real gdx,gdy

c functions
c	integer ipint,ieint
	integer ipext


        if(nlv.ne.nlev) stop 'error stop conzstab: level'

c----------------------------------------------------------------
c initialize variables and parameters
c----------------------------------------------------------------

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

c	----------------------------------------------------------------
c	global arrays for accumulation of implicit terms
c	----------------------------------------------------------------

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

c	----------------------------------------------------------------
c	aux elements inside element
c	----------------------------------------------------------------

c	these are aux arrays (bigger than needed) to avoid checking for
c	what layer we are in -> we never get out of bounds

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

c	these are the local arrays for accumulation of implicit terms
c	(maybe we do not need them, but just to be sure...)
c	after accumulation we copy them onto the global arrays

        do l=1,nlv
	  do ii=1,3
	    cle(l,ii) = 0.
	    clc(l,ii) = 0.
	    clm(l,ii) = 0.
	    clp(l,ii) = 0.
	  end do
	end do

c	----------------------------------------------------------------
c	define vertical velocities
c	----------------------------------------------------------------

c----------------------------------------------------------------
c loop over elements
c----------------------------------------------------------------

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

c	----------------------------------------------------------------
c	set up vectors for use in assembling contributions
c	----------------------------------------------------------------

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

c	----------------------------------------------------------------
c	set vertical velocities in surface and bottom layer
c	----------------------------------------------------------------

c	we do not set wl(0,ii) because otherwise we loose concentration
c	through surface
c
c	we set wl(ilevel,ii) to 0 because we are on the bottom
c	and there should be no contribution from this element
c	to the vertical velocity

	do ii=1,3
	  wl(ilevel,ii) = 0.
	end do

c	----------------------------------------------------------------
c	compute vertical fluxes (w/o vertical TVD scheme)
c	----------------------------------------------------------------

	call vertical_flux_ie(btvdv,ie,ilevel,dt,wws,cl,wl,hold,vflux)

c----------------------------------------------------------------
c loop over levels
c----------------------------------------------------------------

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

c	  ----------------------------------------------------------------
c	  initialization to be sure we are in a clean state
c	  ----------------------------------------------------------------

	  fw(ii) = 0.
	  cle(l,ii) = 0.
	  clc(l,ii) = 0.
	  clm(l,ii) = 0.
	  clp(l,ii) = 0.

c	  ----------------------------------------------------------------
c	  contributions from horizontal diffusion
c	  ----------------------------------------------------------------

          waux = 0.
          do iii=1,3
            waux = waux + wdifhv(iii,ii,ie) * cl(l,iii)
          end do
          wdiff(ii) = waux

c	  ----------------------------------------------------------------
c	  contributions from vertical diffusion
c	  ----------------------------------------------------------------

c	  in fd(ii) is explicit contribution
c	  the sign is for the term on the left side, therefore
c	  fd(ii) must be subtracted from the right side
c
c	  maybe we should use real layer thickness, or even the
c	  time dependent layer thickness

	  rvptop = difv(l-1,k) + difmol
	  rvpbot = difv(l,k) + difmol
	  !hmtop = 2. * rvptop * present(l-1) / (hdv(l-1)+hdv(l))
	  !hmbot = 2. * rvpbot * present(l+1) / (hdv(l)+hdv(l+1))
	  hmotop =2.*rvptop*present(l-1)/(hold(l-1,ii)+hold(l,ii))
	  hmobot =2.*rvpbot*present(l+1)/(hold(l,ii)+hold(l+1,ii))
	  hmntop =2.*rvptop*present(l-1)/(hnew(l-1,ii)+hnew(l,ii))
	  hmnbot =2.*rvpbot*present(l+1)/(hnew(l,ii)+hnew(l+1,ii))

	  fd(ii) = adt * ( 
     +			(cl(l,ii)-cl(l+1,ii))*hmobot -
     +			(cl(l-1,ii)-cl(l,ii))*hmotop
     +			  )

	  clc(l,ii) = clc(l,ii) + ad * ( hmntop + hmnbot )
	  clm(l,ii) = clm(l,ii) - ad * ( hmntop )
	  clp(l,ii) = clp(l,ii) - ad * ( hmnbot )

c	  ----------------------------------------------------------------
c	  contributions from vertical advection
c	  ----------------------------------------------------------------

c	  in fw(ii) is explicit contribution
c	  the sign is for the term on the left side, therefore
c	  fw(ii) must be subtracted from the right side
c
c	  if we are in last layer, w(l,ii) is zero
c	  if we are in first layer, w(l-1,ii) is zero (see above)

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

c	  if( .not. btvdv ) then
c	  if( flux_tot .ne. flux_tot1 .or. flux_tot .ne. fw(ii) ) then
c	   !if( ie .eq. 100 ) then
c	    write(6,*) '********** vflux   ',ie,ii,l,ilevel
c	    write(6,*) fw(ii),flux_tot1,flux_tot
c	    write(6,*) flux_top,flux_bot
c	    do ll=0,ilevel
c	    write(6,*) ll,vflux(ll,ii)
c	    end do
c	   !end if
c	  end if
c	  end if

	  fw(ii) = flux_tot
	end do

c	----------------------------------------------------------------
c	contributions from horizontal advection (only explicit)
c	----------------------------------------------------------------
c
c	f(ii) > 0 ==> flux into node ii
c	itot=1 -> flux out of one node
c		compute flux with concentration of this node
c	itot=2 -> flux into one node
c		for flux use conz. of the other two nodes and
c		minus the sum of these nodes for the flux of this node

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

c	----------------------------------------------------------------
c	horizontal TVD scheme start
c	----------------------------------------------------------------

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

c	----------------------------------------------------------------
c	horizontal TVD scheme finish
c	----------------------------------------------------------------

c	----------------------------------------------------------------
c	contributions from nudging
c	----------------------------------------------------------------

	do ii=1,3
	  fnudge(ii) = robs * rtau(l,ii) * ( cob(l,ii) - cl(l,ii) )
	end do

c	----------------------------------------------------------------
c	sum explicit contributions
c	----------------------------------------------------------------

	do ii=1,3
          hmed = haver(l)                    !new ggu   !HACK
	  cexpl = aj4 * ( hold(l,ii)*cl(l,ii)
     +				+ dt *  ( 
     +					    hold(l,ii)*fnudge(ii)
     +					  + 3.*fl(ii) 
     +					  - fw(ii)
     +					  - rk3*hmed*wdiff(ii)
     +					  - fd(ii)
     +					)
     +		         )
	  cle(l,ii) = cle(l,ii) + cexpl
	  !k=kn(ii)
	  !cn(l,k) = cn(l,k) + cexpl
	end do

	end do		! loop over l

c----------------------------------------------------------------
c end of loop over l
c----------------------------------------------------------------

c----------------------------------------------------------------
c set up implicit contributions
c----------------------------------------------------------------

c cdiag is diagonal of tri-diagonal system
c chigh is high (right) part of tri-diagonal system
c clow is low (left) part of tri-diagonal system
c
c clp -> bottom
c clm -> top

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

c----------------------------------------------------------------
c end of loop over elements
c----------------------------------------------------------------

c in cdiag, chigh, clow is matrix (implicit part)
c if explicit calculation, chigh=clow=0 and in cdiag is volume of node [m**3]
c in cnv is mass of node [kg]
c for explicit treatment, cnv/cdiag gives new concentration [kg/m**3]

c----------------------------------------------------------------
c accumulate contributions on each node
c----------------------------------------------------------------

! here we could make just one loop over nkn
! this would include this loop, the boundary conditions
! and then  the solution of the vertical system
! in this case we would only need the following arrays:
! cn(l),clow(l),chigh(l),cdiag(l) (one dimensional arrays over the vertical)

	do k=1,nkn
	  n = ilinkv(k+1)-ilinkv(k)
	  ibase = ilinkv(k)
	  if( lenkv(ibase+n) .le. 0 ) n = n - 1
	  ilevel = ilhkv(k)
	  do i=1,n
	    ie = lenkv(ibase+i)
	    !ii = lenkiiv(ibase+i)
	    !the next 4 lines are a hack, they will disappear in the new dist
	    ii = 0
	    if( nen3v(1,ie) == k ) ii = 1
	    if( nen3v(2,ie) == k ) ii = 2
	    if( nen3v(3,ie) == k ) ii = 3
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

c----------------------------------------------------------------
c integrate boundary conditions
c----------------------------------------------------------------

c in case of negative flux (qflux<0) must check if node is OBC (BUG_2010_01)

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

c----------------------------------------------------------------
c compute concentration for each node (solve system)
c----------------------------------------------------------------

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

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*****************************************************************

!        subroutine conzstab(cn1,co1
!     +			,ddt
        subroutine conzstab(
     +			ddt
     +			,robs,wsink,wsinkv
     +                  ,rkpar,difhv,difv
     +			,difmol,azpar
     +			,adpar,aapar
     +                  ,sindex
     +			,istot,isact
     +			,nlvddi,nlev)
c
c checks stability
c
c cn     new concentration
c co     old concentration
c caux   aux vector
c clow	 lower diagonal of vertical system
c chig	 upper diagonal of vertical system
c ddt    time step
c robs	 factor for nudging
c rkpar  horizontal turbulent diffusivity
c difhv  horizontal turbulent diffusivity (variable between elements)
c difv   vertical turbulent diffusivity
c difmol vertical molecular diffusivity
c azpar  time weighting parameter
c adpar  time weighting parameter for vertical diffusion (ad)
c aapar  time weighting parameter for vertical advection (aa)
c sindex stability index
c istot	 total inter time steps
c isact	 actual inter time step
c nlvddi	 dimension in z direction
c nlv	 actual needed levels
c
c written 09.01.94 by ggu  (from scratch)
c revised 19.01.94 by ggu  $$flux - flux conserving property
c revised 20.01.94 by ggu  $$iclin - iclin not used to compute volume
c revised 20.01.94 by ggu  $$lumpc - evaluate conz. nodewise
c revised 03.02.94 by ggu  $$itot0 - exception for itot=0 or 3
c revised 04.02.94 by ggu  $$fact3 - factor 3 missing in transport
c revised 04.02.94 by ggu  $$azpar - azpar used to compute transport
c revised 04.02.94 by ggu  $$condry - comute conz also in dry areas
c revised 07.02.94 by ggu  $$istot - istot for fractional time step
c revised 01.06.94 by ggu  restructured for 3-d model
c revised 18.07.94 by ggu  $$htop - use htop instead of htopo for mass cons.
c revised 09.04.96 by ggu  $$rvadj adjust rv in certain areas
c
c solution of purely diffusional part :
c
c dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
c
c C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )
c
c for n-dimensions and
c
c C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )
c
c for 1 dimension
c
c the solution is normalized, i.e.  int(C(x,t)dx) = 1 over the whole area
c
c DPGGU -> introduced double precision to stabilize solution

	use mod_bound_geom
	use mod_depth
	use mod_layer_thickness
	use mod_diff_aux
	use mod_bound_dynamic
	use mod_aux_array
	use mod_area
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none
c
c parameters
        include 'param.h'
c arguments
	integer nlvddi,nlev
        !real cn1(nlvddi,1),co1(nlvddi,1)		!DPGGU
        real difv(0:nlvddi,1)
        real difhv(nlvddi,1)
	real difmol
        real ddt,rkpar,azpar,adpar,aapar			!$$azpar
	real robs,wsink
	real wsinkv(0:nlvddi,1)
	integer istot,isact
c common
	include 'femtime.h'
	include 'mkonst.h'

	!include 'hydro_print.h'



 
 




c local
	logical bdebug,bdebug1,debug
	integer k,ie,ii,l,iii
	integer lstart
	integer ilevel
	integer itot,isum	!$$flux
	logical berror
	integer kn(3)
        real rkmin,rkmax
        real sindex,rstol
	double precision us,vs
	double precision az,azt
	double precision aa,aat,ad,adt
	double precision aj,rk3,rv,aj4
	double precision hmed,hmbot,hmtop
	double precision rvptop,rvpbot
	double precision dt,w,aux
	double precision aux1,aux2,aux3,aux4,aux5
	double precision b(3),c(3),f(3)
	double precision rso,rsn,rsot,rsnt,rstot
	double precision hn,ho
        double precision wdiff(3)
        double precision wws

	double precision difabs,difrel,volold,volnew,flxin,flxtot,diff
	double precision stabind,stabadv,stabdiff,stabvert,stabpoint
	double precision voltot

c------------------------------------------------------------
c big arrays
c------------------------------------------------------------
	double precision cn(nlvddi,nkn)		!DPGGU	!FIXME
	double precision co(nlvddi,nkn)
	double precision cdiag(nlvddi,nkn)
	double precision clow(nlvddi,nkn)
	double precision chigh(nlvddi,nkn)
        real cwrite(nlvddi,nkn)
c------------------------------------------------------------
c end of big arrays
c------------------------------------------------------------

	double precision cexpl
	double precision fw(3),fd(3)
	double precision fl(3)
c local (new)
	double precision clc(nlvddi,3), clm(nlvddi,3), clp(nlvddi,3)
	!double precision cl(0:nlvddi+1,3)
	double precision wl(0:nlvddi+1,3)
c
	double precision hdv(0:nlvddi+1)
	double precision haver(0:nlvddi+1)
	double precision hnew(0:nlvddi+1,3)
	double precision hold(0:nlvddi+1,3)
	double precision present(0:nlvddi+1)

        integer kstab
	real dtorig

        !integer iustab
        !save iustab
        !data iustab /0/
c functions
	logical is_zeta_bound
	real getpar

	!write(6,*) 'conzstab called...'

        if(nlv.ne.nlev) stop 'error stop conzstab: level'

c-----------------------------------------------------------------
c initialization
c-----------------------------------------------------------------

        bdebug1 = .true.
        bdebug1 = .false.
        debug = .true.
        debug = .false.
	bdebug=.false.
	berror=.false.

        if( bdebug1 ) then
                write(6,*) 'debug parameters in conz3d'
		write(6,*) ddt,rkpar,difmol,azpar,adpar,aapar
                write(6,*) istot,isact,nlvddi,nlv
                write(6,*) nkn,nel
        end if

	az=azpar		!$$azpar
	azt=1.-az
	ad=adpar
	adt=1.-ad
	aa=aapar
	aat=1.-aa

	wws = wsink
	wws = 0.

	if( aa .ne. 0. .and. nlv .gt. 1 ) then
	  write(6,*) 'aapar = ',aapar
	  write(6,*) 'Cannot use implicit vertical advection.'
	  write(6,*) 'This might be resolved in a future version.'
	  write(6,*) 'Please set aapar = 0 in the STR file.'
	  stop 'error stop conzstab: implicit vertical advection'
	end if

c	-----------------------------------------------------------------
c	 fractional time step
c	-----------------------------------------------------------------

	rstot=istot		!$$istot
	rso=(isact-1)/rstot
	rsn=(isact)/rstot
	rsot=1.-rso
	rsnt=1.-rsn

	dt=ddt/rstot

c	-----------------------------------------------------------------
c	 initialize global arrays for accumulation of implicit terms
c	-----------------------------------------------------------------

	do k=1,nkn
          do l=1,nlv
	    !co(l,k)=cn1(l,k)	!DPGGU	!not used for stability
            cn(l,k)=0.          !Malta
            co(l,k)=0.
	    if( mfluxv(l,k) .gt. 0. ) co(l,k) = mfluxv(l,k)	!point sources
            cdiag(l,k)=0.
            clow(l,k)=0.
            chigh(l,k)=0.
            cwrite(l,k)=0.
          end do
	end do

c	-----------------------------------------------------------------
c	these are aux arrays (bigger than needed) to avoid checking for
c	what layer we are in -> we never get out of bounds
c	-----------------------------------------------------------------

        do l=0,nlv+1
	  hdv(l) = 0.		!layer thickness
          haver(l) = 0.
	  present(l) = 0.	!1. if layer is present
	  do ii=1,3
	    hnew(l,ii) = 0.	!as hreal but with zeta_new
	    hold(l,ii) = 0.	!as hreal but with zeta_old
	    !cl(l,ii) = 0.	!concentration in layer
	    wl(l,ii) = 0.	!vertical velocity
	  end do
	end do

c	-----------------------------------------------------------------
c	these are the local arrays for accumulation of implicit terms
c	(maybe we do not need them, but just to be sure...)
c	after accumulation we copy them on the global arrays
c	-----------------------------------------------------------------

        do l=1,nlv
	  do ii=1,3
	    clc(l,ii) = 0.
	    clm(l,ii) = 0.
	    clp(l,ii) = 0.
	  end do
	end do

c	-----------------------------------------------------------------
c	vertical velocities
c	-----------------------------------------------------------------

c-----------------------------------------------------------------
c loop over elements
c-----------------------------------------------------------------

        do ie=1,nel

	do ii=1,3
          k=nen3v(ii,ie)
	  kn(ii)=k
	  b(ii)=ev(ii+3,ie)
	  c(ii)=ev(ii+6,ie)
	end do

	aj=ev(10,ie)    !area of triangle / 12
	aj4=4.*aj
        ilevel=ilhv(ie)

c set up vectors for use in assembling contributions

        do l=1,ilevel
	  hdv(l) = hdeov(l,ie)		!use old time step -> FIXME
          haver(l) = 0.5 * ( hdeov(l,ie) + hdenv(l,ie) )
	  present(l) = 1.
	  do ii=1,3
	    k=kn(ii)
	    hn = hdknv(l,k)		! there are never more layers in ie
	    ho = hdkov(l,k)		! ... than in k
	    hold(l,ii) = rso * hn + rsot * ho
	    hnew(l,ii) = rsn * hn + rsnt * ho
	    !cl(l,ii) = co(l,k)
	    wl(l,ii) = wlnv(l,k) - wsink * wsinkv(l,k)
	  end do
	end do

	do l=ilevel+1,nlv
	  present(l) = 0.
	end do

c	set vertical velocities in surface and bottom layer
c
c	we do not set wl(0,ii) because otherwise we loose concentration
c	through surface
c
c	we set wl(ilevel,ii) to 0 because we are on the bottom
c	and there should be no contribution from this element
c	to the vertical velocity

	do ii=1,3
	  wl(ilevel,ii) = 0.
	end do

c-----------------------------------------------------------------
c loop over levels
c-----------------------------------------------------------------

        do l=1,ilevel

        us=az*utlnv(l,ie)+azt*utlov(l,ie)             !$$azpar
        vs=az*vtlnv(l,ie)+azt*vtlov(l,ie)

        rk3 = 3. * rkpar * difhv(l,ie)

	itot=0
	isum=0
	do ii=1,3
	  k=kn(ii)
	  f(ii)=us*b(ii)+vs*c(ii)	!$$azpar
	  if(f(ii).lt.0.) then	!flux out of node
	    itot=itot+1
	    isum=isum+ii
	  end if

c new weights for diffusion

          wdiff(ii) = wdifhv(ii,ii,ie)

c	  initialization to be sure we are in a clean state

	  fw(ii) = 0.
	  clc(l,ii) = 0.
	  clm(l,ii) = 0.
	  clp(l,ii) = 0.

c	  contributions from vertical advection
c
c	  in fw(ii) is explicit contribution
c	  the sign is for the term on the left side, therefore
c	  fw(ii) must be subtracted from the right side
c
c	  if we are in last layer, w(l,ii) is zero
c	  if we are in first layer, w(l-1,ii) is zero (see above)

	  w = wl(l-1,ii) - wws		!top of layer
	  if( w .gt. 0. ) then          !out
	    fw(ii) = aat*w
	    clc(l,ii) = clc(l,ii) + aa*w
	  else
	    fw(ii) = 0.
	    clm(l,ii) = clm(l,ii) + aa*w
	  end if

	  w = wl(l,ii) - wws		!bottom of layer
	  if( w .gt. 0. ) then
	    clp(l,ii) = clp(l,ii) - aa*w
	  else
	    fw(ii) = fw(ii) - aat*w
	    clc(l,ii) = clc(l,ii) - aa*w
	  end if

c	  contributions from vertical diffusion
c
c	  in fd(ii) is explicit contribution
c	  the sign is for the term on the left side, therefore
c	  fd(ii) must be subtracted from the right side
c
c	  maybe we should use real layer thickness, or even the
c	  time dependent layer thickness

	  rvptop = difv(l-1,k) + difmol
	  rvpbot = difv(l,k) + difmol
	  hmtop = 2. * rvptop * present(l-1) / (hdv(l-1)+hdv(l))
	  hmbot = 2. * rvpbot * present(l+1) / (hdv(l)+hdv(l+1))

          fd(ii) = adt * ( hmtop + hmbot )

	  clc(l,ii) = clc(l,ii) + ad * ( hmtop + hmbot )
	  clm(l,ii) = clm(l,ii) - ad * ( hmtop )
	  clp(l,ii) = clp(l,ii) - ad * ( hmbot )
	end do

c sum explicit contributions

	do ii=1,3
	  k=kn(ii)
          hmed = hold(l,ii)                      !new ggu   !HACK
          cexpl = dt * aj4 * rk3 * hmed * wdiff(ii)	!bug fix 12.2.2010
	  clow(l,k) = clow(l,k) + cexpl
          cexpl = dt * aj4 * 3. * f(ii)
          if( cexpl .lt. 0. ) then             !flux out of node
	    !chigh(l,k) = chigh(l,k) - cexpl
          end if
          if( cexpl .gt. 0. ) then             !flux into node
	    chigh(l,k) = chigh(l,k) + cexpl
          end if
          cn(l,k) = cn(l,k) + dt * aj4 * ( fw(ii) + fd(ii) )
          co(l,k) = co(l,k) + dt * aj4 * hmed * robs * rtauv(l,k) !nudging
	end do

	end do		! loop over l

c set up implicit contributions
c
c cdiag is diagonal of tri-diagonal system
c chigh is high (right) part of tri-diagonal system
c clow is low (left) part of tri-diagonal system

	do ii=1,3
	  clm(1,ii) = 0.
	  clp(ilevel,ii) = 0.
	end do

        do l=1,ilevel
	  do ii=1,3
	    k=kn(ii)
	    !clow(l,k)  = clow(l,k)  + aj4 * dt * clm(l,ii)
	    !chigh(l,k) = chigh(l,k) + aj4 * dt * clp(l,ii)
	    !cdiag(l,k) = cdiag(l,k) + aj4 * dt * clc(l,ii)
	    !clow(l,k)  = clow(l,k)  + aj4 * hold(l,ii)
            hmed = min(hold(l,ii),hnew(l,ii))
	    cdiag(l,k) = cdiag(l,k) + aj4 * hmed
	  end do
	end do

	end do		! loop over ie

c-----------------------------------------------------------------
c compute stability
c
c cdiag		volume of cell
c chigh		flux due to horizontal advection
c clow		flux due to horizontal diffusion
c cn		flux due to vertical advection and diffusion (explicit)
c co		flux due to point sources and nudging
c-----------------------------------------------------------------

        stabind = 0.		!total max stability index
        stabadv = 0.		!advective max stability index
        stabdiff = 0.		!diffusive max stability index
        stabvert = 0.		!vertical max stability index
        stabpoint = 0.		!point source max stability index
        kstab = 0		!node with highest stabind

	do k=1,nkn
	  !bdebug1 = k .eq. 1402
	  !bdebug1 = k .eq. 1405
	  bdebug1 = k .eq. -1
	  ilevel = ilhkv(k)
          if( .not. is_zeta_bound(k) ) then	!FIXME
          !if( is_inner(k) ) then	!FIXME
          !if( .true. ) then	!FIXME
	   do l=1,ilevel
            voltot = cdiag(l,k)
            flxtot = chigh(l,k) + clow(l,k) + cn(l,k) + co(l,k)
	    if( bdebug1 ) write(99,*) k,l,voltot,flxtot
            if( voltot .gt. 0. ) then
                  aux1 = flxtot / voltot
                  if( aux1 .gt. stabind ) kstab = k
                  stabind = max(stabind,aux1)
		  cwrite(l,k) = aux1		!save for write
                  aux2 = chigh(l,k) / voltot
                  stabadv = max(stabadv,aux2)
		  saux1(l,k) = aux2		!for adv. stab.
                  aux3 = clow(l,k) / voltot
                  stabdiff = max(stabdiff,aux3)
                  aux4 = cn(l,k) / voltot
                  stabvert = max(stabvert,aux4)
                  aux5 = co(l,k) / voltot
                  stabpoint = max(stabpoint,aux5)
	          if( bdebug1 ) write(99,*) aux1,aux2,aux3,aux4,aux5

c		  aux=flxtot / voltot
c		  if( 300*aux/dt .gt. 1000 ) then
c			  write(6,*) is_boundary(k)
c			  write(6,*) is_external_boundary(k)
c			  write(6,*) is_internal_boundary(k)
c			  write(6,*) is_inner(k)
c			  write(6,*) stabind,stabadv,stabdiff,stabvert
c			  call check_set_unit(6)
c			  call check_node(k)
c		  end if

            else
		  cwrite(l,k) = 0
		  saux1(l,k) = 0.
            end if
	   end do
          else
	   do l=1,ilevel
		  cwrite(l,k) = 0
		  saux1(l,k) = 0.
           end do
          end if
	end do

c        write(6,*) 'stab check: ',nkn,nlv
c        call check2Dr(nlvddi,nlv,nkn,cwrite,0.,0.,"NaN check","cstab")

c-----------------------------------------------------------------
c in stabind is stability index (advection and diffusion)
c in cdiag is the local value of the stability index
c in cwrite is the value of the stability index for each node
c
c istot  is saved and returned from subroutine (number of iterations)
c sindex is saved and returned from subroutine (stability index)
c-----------------------------------------------------------------

	rstol = getpar('rstol')
        istot = 1 + stabind / rstol
        sindex = stabind

	call get_orig_timestep(dtorig)
	call output_stability_node(dtorig,cwrite)

c        if( .false. ) then
c        !if( idt .le. 3 ) then
c	  write(6,*) 'kstab = ',kstab,'  stabind = ',stabind
c          call conwrite(iustab,'.stb',1,777,nlvddi,cwrite)
c        end if

c        call stb_histo(it,nlvddi,nkn,ilhkv,cwrite)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine massconc(mode,cn,nlvddi,mass)

c computes total mass of conc

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c arguments
	integer mode
	integer nlvddi
	real cn(nlvddi,1)
	real mass
c parameter
        include 'param.h'
c common
	include 'femtime.h'
c local
	integer k,l,lmax
        double precision vol
	double precision sum,masstot
	real volnode

        masstot = 0.

        do k=1,nkn
	  lmax = ilhkv(k)
          sum = 0.
          do l=1,lmax
            vol = volnode(l,k,mode)
            sum = sum + cn(l,k) * vol
          end do
          masstot = masstot + sum
        end do

	mass = masstot

c	write(88,*) 'tot mass: ',it,mass

	end

c**************************************************************

        subroutine check_scal_bounds(cnv,cmin,cmax,eps,bstop)

c checks if scalar is out of bounds

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'

        real cnv(nlvdi,nkn)
        real cmin,cmax
	real eps
	logical bstop		!stop simulation if true


	logical berror
        integer k,l,lmax,kext
        real cc,cn,cx,cd

        integer ipext

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

c*****************************************************************

	subroutine assert_min_max_property(cnv,cov,sbconz,rmin,rmax,eps)

c checks min/max property

	use mod_bound_geom
	use mod_bound_dynamic
	use levels
	use basin

	implicit none

	include 'param.h'

        real cnv(nlvdi,nkn)			!new concentration
        real cov(nlvdi,nkn)			!old concentration
        real sbconz(nlvdi,nkn)		!conz boundary conditions
	real rmin(nlvdi,nkn)		!aux arrray to contain min
	real rmax(nlvdi,nkn)		!aux arrray to contain max
	real eps

	include 'femtime.h'


	logical bwrite,bstop
	integer k,ie,l,ii,lmax,ierr
	integer levdbg
	real amin,amax,c,qflux,dmax
	real drmax,diff
	real dt

	integer ipext
	logical is_zeta_bound
	real getpar

	bwrite = .true.		! write every violation
	bwrite = .false.		! write every violation
	bstop = .false.		! stop after error

	levdbg = nint(getpar('levdbg'))

c---------------------------------------------------------------
c vertical contribution (normally implicit -> whole column)
c---------------------------------------------------------------

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

c---------------------------------------------------------------
c point sources
c---------------------------------------------------------------

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

c---------------------------------------------------------------
c horizontal contribution
c---------------------------------------------------------------

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

c---------------------------------------------------------------
c check with new concentration
c---------------------------------------------------------------

	ierr = 0
	dmax = 0.
	drmax = 0.

	do k=1,nkn
	 !if( .not. is_external_boundary(k) ) then	!might be relaxed
	 if( .not. is_zeta_bound(k) ) then	!might be relaxed
	  lmax = ilhkv(k)
	  do l=1,lmax
	    c = cnv(l,k)
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

c---------------------------------------------------------------
c check error condition
c---------------------------------------------------------------

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

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine stb_histo(it,nlvddi,nkn,ilhkv,cwrite)

c writes histogram info about stability index

        implicit none

        integer it
        integer nlvddi,nkn
        integer ilhkv(1)
        real cwrite(nlvddi,1)

        integer ndim
        parameter(ndim=11)

        integer nbin
        real aux
        integer k,l,lmax

        integer ic(ndim+1)
        real bins(ndim)
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

c*****************************************************************

