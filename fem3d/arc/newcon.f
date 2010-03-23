c
c $Id: newcon.f,v 1.43 2009-01-26 15:04:57 georg Exp $
c
c routines for concentration
c
c contents :
c
c subroutine scal3sh(what,cnv,nlvbnd,rcv,rkpar,difhv,difv,difmol)
c						shell for scalar T/D
c subroutine conz3d(cn,co,ddt,rkpar,difhv,difv,difmol
c               ,azpar,adpar,aapar,istot,isact,nlvdi,nlv)   
c						computes concentration
c subroutine massconc(mode,cn,nlvdi,res)
c						computes total mass of conc
c subroutine conzstab(cn,co
c    +                  ,ddt
c    +                  ,rkpar,difhv,difv
c    +                  ,difmol,azpar
c    +                  ,adpar,aapar,sindex
c    +                  ,istot,isact
c    +                  ,nlvdi,nlv)
c						checks stability
c subroutine get_stability(dt,rkpar,difhv,rindex,istot)
c               computes stability index with given horizontal diffusion
c subroutine get_stab_index(dt,rindex,istot)
c               computes stability index only for advection
c
c stb_histo(it,nlvdi,nkn,ilhkv,cwrite)
c                                               stability histogram
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
c 09.09.2003    ggu     call to scal3sh changed -> new arg nlvbnd
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
c
c*********************************************************************

	subroutine scal_adv(what,ivar
     +				,scal,bnd3
     +				,rkpar,wsink
     +                          ,difhv,difv,difmol)

c shell for scalar (for parallel version)

        include 'param.h'

        character*(*) what
	integer ivar
        real scal(nlvdim,1)
        real bnd3(nb3dim,0:nbcdim)

        real rkpar
	real wsink
        real difhv(nlvdim,1)
	real difv(0:nlvdim,1)
        real difmol

	real bnd3_aux(nb3dim)
        real r3v(nlvdim,nkndim)

        integer iwhat,ichanm
	character*10 whatvar,whataux

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

	call bnds_trans(whatvar(1:iwhat)
     +				,nb3dim,bnd3,bnd3_aux
     +                          ,ivar,nlvdim,r3v)

c--------------------------------------------------------------
c do advection and diffusion
c--------------------------------------------------------------

        call scal3sh(whatvar(1:iwhat)
     +				,scal,nlvdim
     +                          ,r3v
     +				,rkpar,wsink
     +                          ,difhv,difv,difmol)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*********************************************************************

	subroutine scal_bnd(what,t,bnd3)

c sets boundary conditions for scalar

        include 'param.h'

        character*(*) what
	real t
        real bnd3(nb3dim,0:nbcdim)

	real bnd3_aux(nb3dim)

	call bnds_set(what,t,nb3dim,bnd3,bnd3_aux)

	end

c*********************************************************************

	subroutine scal_adv_fact(what,ivar,fact
     +				,scal,bnd3
     +				,rkpar,wsink
     +                          ,difhv,difv,difmol)

c shell for scalar (for parallel version)
c
c special version for cohesive sediments with factor

        include 'param.h'

        character*(*) what
	integer ivar
	real fact			!factor for boundary condition
        real scal(nlvdim,1)
        real bnd3(nb3dim,0:nbcdim)

        real rkpar
	real wsink
        real difhv(nlvdim,1)
	real difv(0:nlvdim,1)
        real difmol

	real bnd3_aux(nb3dim)
        real r3v(nlvdim,nkndim)

	character*20 whatvar,whataux

c--------------------------------------------------------------
c make identifier for variable
c--------------------------------------------------------------

	whatvar = what
	if( ivar .ne. 0 ) then
          write(whataux,'(i2)') ivar
          whatvar = what // whataux
	end if

c--------------------------------------------------------------
c transfer boundary conditions of var ivar to 3d matrix r3v
c--------------------------------------------------------------

	call bnds_trans(whatvar
     +				,nb3dim,bnd3,bnd3_aux
     +                          ,ivar,nlvdim,r3v)

c--------------------------------------------------------------
c multiply boundary condition with factor
c--------------------------------------------------------------

	call mult_scal_bc(r3v,fact)

c--------------------------------------------------------------
c do advection and diffusion
c--------------------------------------------------------------

        call scal3sh(whatvar
     +				,scal,nlvdim
     +                          ,r3v
     +				,rkpar,wsink
     +                          ,difhv,difv,difmol)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	end

c*********************************************************************

	subroutine scal3sh(what,cnv,nlvbnd,rcv,rkpar,wsink
     +					,difhv,difv,difmol)

c shell for scalar T/D

	implicit none

c parameter
        include 'param.h'
c arguments
        character*(*) what
        real cnv(nlvdim,1)
	integer nlvbnd		!vertical dimension of boundary condition
        real rcv(nlvbnd,1)
        real rkpar
	real wsink
        real difhv(nlvdim,1)
	real difv(0:nlvdim,1)
        real difmol
c parameters
	integer istot_max
	!parameter ( istot_max = 100 )
	!parameter ( istot_max = 200 )
	!parameter ( istot_max = 300 )
	parameter ( istot_max = 1000 )
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real uprv(nlvdim,1), vprv(nlvdim,1)
        common /uprv/uprv, /vprv/vprv

c local
        real saux(nlvdim,nkndim)		!aux array
        real sbflux(nlvdim,nkndim)		!flux boundary conditions
	real gradxv(nlvdim,nkndim)		!gradient in x for tvd
	real gradyv(nlvdim,nkndim)		!gradient in y for tvd

	logical btvd
	integer isact
	integer istot
	integer itvd
	integer iuinfo
        real dt
        real sindex
	real mass,massold,massdiff
	real azpar,adpar,aapar
	real ssurface
c function
	real getpar

c-------------------------------------------------------------
c start of routine
c-------------------------------------------------------------

	if(nlvdim.ne.nlvdi) stop 'error stop scal3sh: level dimension'

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	call getaz(azpar)
	adpar=getpar('adpar')
	aapar=getpar('aapar')
	itvd=nint(getpar('itvd'))
	btvd = itvd .gt. 0

        call getinfo(iuinfo)  !unit number of info file

c-------------------------------------------------------------
c check stability criterion -> set istot
c-------------------------------------------------------------

	call get_timestep(dt)

        call get_stability(dt,rkpar,sindex,istot)

        write(iuinfo,*) 'stability ',what,':',it,sindex,istot

        if( istot .gt. istot_max ) then
            write(6,*) 'istot  = ',istot,'   sindex = ',sindex
            stop 'error stop scal3sh: istot index too high'
        end if

c-------------------------------------------------------------
c set up flux boundary conditions (temporary) -> put in sbflux
c-------------------------------------------------------------

	ssurface = 0.	!FIXME - HACK
	if( what .eq. 'temp' ) ssurface = -999.

	call make_scal_flux(what,rcv,cnv,sbflux,ssurface)

c-------------------------------------------------------------
c transport and diffusion
c-------------------------------------------------------------

	call massconc(-1,cnv,nlvdim,massold)

	do isact=1,istot

	if( btvd ) call grad_tvd(cnv,gradxv,gradyv,saux,nlvdi,nlv)

          call conz3d(
     +           cnv
     +          ,saux
     +          ,dt
     +          ,rkpar,difhv,difv,difmol
     +          ,sbflux
     +		,itvd,gradxv,gradyv
     +		,wsink
     +		,azpar,adpar,aapar
     +          ,istot,isact
     +          ,nlvdi,nlv
     +               )

          call bndo_setbc(it,what,nlvdi,cnv,rcv,uprv,vprv)

	end do

        !if( what .eq. 'salt' ) call check_scal(cnv,30.,40.)

c-------------------------------------------------------------
c check total mass
c-------------------------------------------------------------

	call massconc(+1,cnv,nlvdim,mass)
	massdiff = mass - massold

	write(iuinfo,1000) 'scal3sh ',what,':'
     +                          ,it,niter,mass,massold,massdiff

	!call write_elem_info(78,9914)

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

        return
 1000   format(a,a,a,2i10,3d13.5)
	end

c**************************************************************

        subroutine check_scal(cnv,cmin,cmax)

c checks if scalar is out of bounds

        implicit none

        include 'param.h'

        real cnv(nlvdim,1)
        real cmin,cmax

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhkv(1)
        common /ilhkv/ilhkv

        integer k,l,lmax,kext
        real cc

        integer ipext

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cc = cnv(l,k)
            if( cc .lt. cmin .or. cc .gt. cmax ) then
                kext = ipext(k)
                write(6,*) 'scalar out of bounds: ',k,kext,l,lmax,cc
            end if
          end do
        end do

        end

c**************************************************************

        subroutine conz3d(cn1,co1
     +			,ddt
     +                  ,rkpar,difhv,difv
     +			,difmol,cbound
     +			,itvd,gradxv,gradyv
     +			,wsink
     +			,azpar,adpar,aapar
     +			,istot,isact
     +			,nlvdi,nlv)
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
c cbound boundary condition (mass flux) [kg/s]
c itvd	 type of transport algorithm used
c gradxv,gradyv  gradient vectors for TVD algorithm
c wsink	 settling velocity [m/s]
c azpar  time weighting parameter
c adpar  time weighting parameter for vertical diffusion (ad)
c aapar  time weighting parameter for vertical advection (aa)
c istot	 total inter time steps
c isact	 actual inter time step
c nlvdi	 dimension in z direction
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

	implicit none
c
c parameters
        include 'param.h'
c arguments
	integer nlvdi,nlv
        real cn1(nlvdi,1),co1(nlvdi,1)		!DPGGU
        real difv(0:nlvdi,1)
        real difhv(nlvdi,1)
	real difmol
        real cbound(nlvdi,1)
	integer itvd
	real gradxv(nlvdi,1)
	real gradyv(nlvdi,1)
	real wsink
        real ddt,rkpar,azpar,adpar,aapar			!$$azpar
	integer istot,isact
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
	integer nen3v(3,1)
	common /nen3v/nen3v
	include 'ev.h'
        real utlnv(nlvdim,1),vtlnv(nlvdim,1)
        common /utlnv/utlnv, /vtlnv/vtlnv
        real utlov(nlvdim,1),vtlov(nlvdim,1)
        common /utlov/utlov, /vtlov/vtlov
        real wprv(0:nlvdim,1)
        common /wprv/wprv
        real wlov(0:nlvdim,1),wlnv(0:nlvdim,1)
        common /wlov/wlov, /wlnv/wlnv
	integer ilhv(1), ilhkv(1)
	common /ilhv/ilhv, /ilhkv/ilhkv

        real zeov(3,1),zenv(3,1)
        common /zeov/zeov, /zenv/zenv

        real hdknv(nlvdim,1)
        common /hdknv/hdknv
        real hdkov(nlvdim,1)
        common /hdkov/hdkov
 
        real hdenv(nlvdim,1)
        common /hdenv/hdenv
        real hdeov(nlvdim,1)
        common /hdeov/hdeov
 
        real areakv(nlvdim,1)
        common /areakv/areakv

        real wdifhv(3,3,1)
        common /wdifhv/wdifhv

c local
	logical bdebug,bdebug1,debug
	integer k,ie,ii,l,iii
	integer lstart
	integer ilevel
	integer itot,isum	!$$flux
	logical berror
	integer kn(3)
        integer ip(3,3)
        integer n,i,ipp
        real rkmin,rkmax
        real mflux
	double precision wws
	double precision us,vs
	double precision az,azt
	double precision aa,aat,ad,adt
	double precision aj,rk3,rv,aj4,aj12
	double precision hmed,hmbot,hmtop
	double precision rvptop,rvpbot
	double precision dt,w,aux
	double precision b(3),c(3),f(3)
	double precision rso,rsn,rsot,rsnt,rstot
	double precision hn,ho

	double precision cn(nlvdim,nkndim)		!DPGGU	!FIXME
	double precision co(nlvdim,nkndim)
	double precision cdiag(nlvdim,nkndim)
	double precision clow(nlvdim,nkndim)
	double precision chigh(nlvdim,nkndim)

c	double precision explh(nlvdim,nlidim)

	double precision cdummy
	double precision cbm,ccm
	double precision fw(3),fd(3)
	double precision fl(3)
        double precision wdiff(3),waux
c local (new)
	double precision clc(nlvdim,3), clm(nlvdim,3), clp(nlvdim,3)
c	double precision clce(nlvdim,3), clme(nlvdim,3), clpe(nlvdim,3)
	double precision cl(0:nlvdim+1,3)
	double precision wl(0:nlvdim+1,3)

	double precision hdv(0:nlvdim+1)
	double precision haver(0:nlvdim+1)
	double precision hnew(0:nlvdim+1,3)
	double precision hold(0:nlvdim+1,3)
	double precision htnew(0:nlvdim+1,3)
	double precision htold(0:nlvdim+1,3)
	double precision present(0:nlvdim+1)

	double precision cauxn(nlvdim)	!FIXME
	double precision cauxd(nlvdim)
	double precision cauxh(nlvdim)
	double precision cauxl(nlvdim)
c tvd
	logical btvd
	integer ic,kc,id,kd
	integer ies
	integer iaux
	real conc,cond,conf
	real gcx,gcy,dx,dy
	real rf,psi
	real fls(3)
	real xgv(1)
	common /xgv/xgv
	real ygv(1)
	common /ygv/ygv

        real alfa,dis
        real ulnv(nlvdim,1)
        common /ulnv/ulnv
        real vlnv(nlvdim,1)
        common /vlnv/vlnv
        real vel
        real gdx,gdy

c functions
c	integer ipint,ieint
	integer ipext

	include 'testbndo.h'

        if(nlvdim.ne.nlvdi) stop 'error stop conz3d: level dimension'

        bdebug1 = .true.
        bdebug1 = .false.
        debug = .false.
        debug = .true.
	bdebug=.false.
	berror=.false.

	btvd = itvd .gt. 0

        if( bdebug1 ) then
                write(6,*) 'debug parameters in conz3d'
		write(6,*) ddt,rkpar,difmol,azpar,adpar,aapar
                write(6,*) istot,isact,nlvdi,nlv
                write(6,*) nkn,nel
        end if

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

	dt=ddt/rstot

c these are the global arrays for accumulation of implicit terms

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

c        call link_fill(n)
c        do i=1,n
c          do l=1,nlv
c	    explh(l,i) = 0.
c          end do
c        end do

c these are aux arrays (bigger than needed) to avoid checking for
c what layer we are in -> we never get out of bounds

        do l=0,nlvdim+1
	  hdv(l) = 0.		!layer thickness
          haver(l) = 0.
	  present(l) = 0.	!1. if layer is present
	  do ii=1,3
	    hnew(l,ii) = 0.	!as hreal but with zeta_new
	    hold(l,ii) = 0.	!as hreal but with zeta_old
	    cl(l,ii) = 0.	!concentration in layer
	    wl(l,ii) = 0.	!vertical velocity
	  end do
	end do

c these are the local arrays for accumulation of implicit terms
c (maybe we do not need them, but just to be sure...)
c after accumulation we copy them on the global arrays

        do l=1,nlvdim
	  do ii=1,3
	    clc(l,ii) = 0.
	    clm(l,ii) = 0.
	    clp(l,ii) = 0.
c	    clce(l,ii) = 0.
c	    clme(l,ii) = 0.
c	    clpe(l,ii) = 0.
	  end do
	end do

c vertical velocities

	do k=1,nkn
	  do l=0,nlv
c	    wprv(l,k) =  	(
c     +				   az*wlnv(l,k) 
c     +				+ azt*wlov(l,k) 
c     +				)
	    wprv(l,k) = wlnv(l,k)	!BUGFIX -> see definition in sp256w
	  end do
	end do

c loop over elements

        do ie=1,nel

c        call node_links(ie,ip)       !pointer into linkv

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

c set up vectors for use in assembling contributions

        do l=1,ilevel
	  hdv(l) = hdeov(l,ie)		!use old time step -> FIXME
          haver(l) = 0.5 * ( hdeov(l,ie) + hdenv(l,ie) )
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
	    wl(l,ii) = wprv(l,k)
	    !wl(l,ii) = 0.				!DDD
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
c	  wl(0,ii) = wprv(0,kn(ii))
	  wl(ilevel,ii) = 0.
	end do

c loop over levels

        do l=1,ilevel

        us=az*utlnv(l,ie)+azt*utlov(l,ie)             !$$azpar
        vs=az*vtlnv(l,ie)+azt*vtlov(l,ie)
	!us = 0.						!DDD
	!vs = 0.

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

c new weights for horizontal diffusion

          hmed = hold(l,ii)
          aux = aj4 * rk3 * hmed
          waux = 0.
          do iii=1,3
            waux = waux + wdifhv(iii,ii,ie) * cl(l,iii)
c            ipp = ip(ii,iii)
c            if( ipp .gt. 0 ) then
c                explh(l,ipp) = explh(l,ipp) + aux * wdifhv(iii,ii,ie)
c            else
c	        clce(l,ii) = clce(l,ii) + aux * wdifhv(iii,ii,ie)
c            end if
          end do
          wdiff(ii) = waux

c	  initialization to be sure we are in a clean state

	  fw(ii) = 0.
	  clc(l,ii) = 0.
	  clm(l,ii) = 0.
	  clp(l,ii) = 0.
c	  clce(l,ii) = 0.
c	  clme(l,ii) = 0.
c	  clpe(l,ii) = 0.

c	  contributions from vertical advection
c
c	  in fw(ii) is explicit contribution
c	  the sign is for the term on the left side, therefore
c	  fw(ii) must be subtracted from the right side
c
c	  if we are in last layer, w(l,ii) is zero
c	  if we are in first layer, w(l-1,ii) is zero (see above)

	  w = wl(l-1,ii) - wws		!top of layer
	  if( l .eq. 1 ) w = 0.		!surface -> no transport (WZERO)
	  if( w .gt. 0. ) then
	    fw(ii) = aat*w*cl(l,ii)
c	    clce(l,ii) = clce(l,ii) - aat*w
	    clc(l,ii) = clc(l,ii) + aa*w
	  else
	    fw(ii) = aat*w*cl(l-1,ii)
c	    clme(l,ii) = clme(l,ii) - aat*w
	    clm(l,ii) = clm(l,ii) + aa*w
	  end if

	  w = wl(l,ii) - wws		!bottom of layer
	  if( l .eq. ilevel ) w = 0.	!bottom -> handle flux elsewhere (WZERO)
	  if( w .gt. 0. ) then
	    fw(ii) = fw(ii) - aat*w*cl(l+1,ii)
c	    clpe(l,ii) = clpe(l,ii) + aat*w
	    clp(l,ii) = clp(l,ii) - aa*w
	  else
	    fw(ii) = fw(ii) - aat*w*cl(l,ii)
c	    clce(l,ii) = clce(l,ii) + aat*w
	    clc(l,ii) = clc(l,ii) - aa*w
	  end if

	!if( w .ne. 0. ) write(6,*) 'wwwwwww ',ie,ii,w	!DDD

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

	  fd(ii) = adt * ( 
     +			(cl(l,ii)-cl(l+1,ii))*hmbot -
     +			(cl(l-1,ii)-cl(l,ii))*hmtop
     +			  )

c	  clce(l,ii) = clce(l,ii) - adt * ( hmtop + hmbot )
c	  clme(l,ii) = clme(l,ii) + adt * ( hmtop )
c	  clpe(l,ii) = clpe(l,ii) + adt * ( hmbot )
	  clc(l,ii) = clc(l,ii) + ad * ( hmtop + hmbot )
	  clm(l,ii) = clm(l,ii) - ad * ( hmtop )
	  clp(l,ii) = clp(l,ii) - ad * ( hmbot )

	!if( hmtop .ne. 0. ) write(6,*) 'ddddddd1 ',ie,ii,hmtop	!DDD
	!if( hmbot .ne. 0. ) write(6,*) 'ddddddd2 ',ie,ii,hmbot	!DDD
	end do

c contribution from horizontal advection (only explicit)
c
c f(ii) > 0 ==> flux into node ii
c itot=1 -> flux out of one node
c	compute flux with concentration of this node
c itot=2 -> flux into one node
c	for flux use conz. of the other two nodes and
c	minus the sum of these nodes for the flux of this node

	if(itot.eq.1) then	!$$flux
          k = kn(isum)
	  fl(1)=f(1)*cl(l,isum)
	  fl(2)=f(2)*cl(l,isum)
	  fl(3)=f(3)*cl(l,isum)
c          do iii=1,3
c            ipp = ip(iii,isum)
c            if( ipp .gt. 0 ) explh(l,ipp) = explh(l,ipp) + aj12 * f(iii)
c          end do
c	  clce(l,isum) = clce(l,isum) + aj12 * f(isum)
	else if(itot.eq.2) then
	  isum=6-isum
	  fl(1)=f(1)*cl(l,1)
	  fl(2)=f(2)*cl(l,2)
	  fl(3)=f(3)*cl(l,3)
	  fl(isum) = 0.
	  fl(isum) = -(fl(1)+fl(2)+fl(3))
c          do iii=1,3
c            if( iii .ne. isum ) then
c	      clce(l,iii) = clce(l,iii) + aj12 * f(iii)
c              ipp = ip(isum,iii)
c              explh(l,ipp) = explh(l,ipp) - aj12 * f(iii)
c            end if
c          end do
	  isum=6-isum
	else			!exception	$$itot0
	  fl(1)=0.
	  fl(2)=0.
	  fl(3)=0.
	end if

	!if( rk3 .ne. 0. ) write(6,*) 'kkkkkkkk ',ie,ii,rk3	!DDD
	do ii=1,3
	!if( fl(ii) .ne. 0. ) write(6,*) 'ffffffff ',ie,ii,fl(ii)	!DDD
	!if( wdiff(ii) .ne. 0. ) write(6,*) 'wdwdwdwd ',ie,ii,wdiff(ii)	!DDD
	end do

	iaux = 0
	do ii=1,3
	  k = nen3v(ii,ie)
	  if( is_external_boundary(k) ) iaux = iaux + 1
	end do

	if( btvd .and. iaux .eq. 0 ) then
	  do ii=1,3
	    fls(ii) = fl(ii)
	    fl(ii) = 0.
	  end do
	  if( itot .eq. 1 ) then ! flux exiting from isum (just one node)
	    ic = isum
	    kc = nen3v(ic,ie)
	    conc = cl(l,ic)
	    do ii=1,3
	      if( ii .ne. ic ) then
		id = ii
	        kd = nen3v(id,ie)
	        cond = cl(l,id)
	        gcx = gradxv(l,kc)
	        gcy = gradyv(l,kc)
		dx = xgv(kd) - xgv(kc)
		dy = ygv(kd) - ygv(kc)
                dis = sqrt(dx**2 +dy**2)
                vel = sqrt(ulnv(l,ie)**2 + vlnv(l,ie)**2)
                alfa = ( dt * vel  ) / dis
		if( conc .eq. cond ) then
		  rf = -1.
		else
		  rf = 2. * (gcx*dx + gcy*dy) / (cond-conc) - 1.
		end if
		psi = max(0.,min(1.,2.*rf),min(2.,rf))  !superbee
c               psi = ( rf + abs(rf)) / ( 1 + abs(rf)) ! muscl
c               psi = max(0.,min(2.,rf)) ! osher
c               psi = max(0.,min(1.,rf)) ! minmod
		conf = conc + 0.5*psi*(cond-conc)*(1-alfa)
	        fl(ic) = fl(ic) - f(id)*conf
	        fl(id) = fl(id) + f(id)*conf
	      end if
	    end do
	  else if( itot .eq. 2 ) then !flux entering into 6-isum (one node)
	    id = 6 - isum
	    kd = nen3v(id,ie)
	    cond = cl(l,id)
	    do ii=1,3
	      if( ii .ne. id ) then
		ic = ii
	        kc = nen3v(ic,ie)
	        conc = cl(l,ic)
	        gcx = gradxv(l,kc)
	        gcy = gradyv(l,kc)
		dx = xgv(kd) - xgv(kc)
		dy = ygv(kd) - ygv(kc)
                dis = sqrt(dx**2 +dy**2)
                vel = sqrt(ulnv(l,ie)**2 + vlnv(l,ie)**2)
                alfa = ( dt * vel  ) / dis
		if( conc .eq. cond ) then
		  rf = -1.
		else
		  rf = 2. * (gcx*dx + gcy*dy) / (cond-conc) - 1.
		end if
		psi = max(0.,min(1.,2.*rf),min(2.,rf))  !superbee
c               psi = ( rf + abs(rf)) / ( 1 + abs(rf)) ! muscl
c               psi = max(0.,min(2.,rf)) ! osher
c               psi = max(0.,min(1.,rf)) ! minmod
		conf = conc + 0.5*psi*(cond-conc)*(1-alfa)
	        fl(ic) = fl(ic) + f(ic)*conf
	        fl(id) = fl(id) - f(ic)*conf
	      end if
	    end do
	  end if
	  aux = 0.
	  !do ii=1,3
	  !  aux = aux + abs(fl(ii)-fls(ii))
	  !end do
	  !if( ie .eq. 2 ) then
	  !  write(6,*) ie,(fl(ii),ii=1,3),(fls(ii),ii=1,3)
	  !end if
	  !if( ie .eq. 2 .and. aux .gt. 1.e-3 ) then
	  !  write(6,*) ie,(fl(ii),ii=1,3),(fls(ii),ii=1,3)
	  !end if
	end if

c sum explicit contributions

	do ii=1,3
	  k=kn(ii)
          hmed = hold(l,ii)                      !new ggu   !HACK
          !hmed = haver(l)                      !new ggu   !HACK
	  cdummy = aj4 * ( hold(l,ii)*cl(l,ii)
     +				+ dt *  ( 3.*fl(ii) 
     +					  - fw(ii)
c     +					  - b(ii)*rk3*hmed*cbm
c     +					  - c(ii)*rk3*hmed*ccm
     +					  + rk3*hmed*wdiff(ii)
     +					  - fd(ii)
     +					)
     +			               )
	  cn(l,k) = cn(l,k) + cdummy
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
	!if( clm(l,ii) .ne. 0. ) write(6,*) 'clm ',l,ie,ii,clm(l,ii)	!DDD
	!if( clp(l,ii) .ne. 0. ) write(6,*) 'clm ',l,ie,ii,clp(l,ii)	!DDD
	!if( clc(l,ii) .ne. 0. ) write(6,*) 'clm ',l,ie,ii,clc(l,ii)	!DDD
	    clow(l,k)  = clow(l,k)  + aj4 * dt * clm(l,ii)
	    chigh(l,k) = chigh(l,k) + aj4 * dt * clp(l,ii)
	    cdiag(l,k) = cdiag(l,k) + aj4 * dt * clc(l,ii)
	    cdiag(l,k) = cdiag(l,k) + aj4 * hnew(l,ii)
	!if( l .eq. 4 .and. cl(l,ii) .lt. 20.4 ) then		!DDD
	!  write(6,*) '44444 ',ie,ii,k,l,cl(l,ii)
	!end if
	  end do
	end do

	end do		! loop over ie

c in cdiag, chigh, clow is matrix (implicit part)
c if explicit calculation, chigh=clow=0 and in cdiag is volume of node [m**3]
c in cnv is mass of node [kg]
c for explicit treatment, cnv/cdiag gives new concentration [kg/m**3]

c integrate boundary conditions

	do k=1,nkn
	  ilevel = ilhkv(k)
	  do l=1,ilevel
            mflux = cbound(l,k)
	    !if( mflux .ne. 0 ) write(6,*) 'mmmmflux ',k,l,mflux
            !if( mflux .ne. 0. ) then
            !  write(96,*) 'mflux ',it,k,l,mflux,co1(l,k)
            !end if
	    ! we always treat mflux explicitly -> should work anyway
            !if( mflux .lt. 0. ) then
            !  cdiag(l,k) = cdiag(l,k) - dt * mflux
            !else
            !  cn(l,k) = cn(l,k) + dt * mflux
            !end if
            cn(l,k) = cn(l,k) + dt * mflux
	  end do
	end do

c compute concentration for each node (solve system)

	if( aa .eq. 0. .and. ad .eq. 0. ) then

	if( nlv .gt. 1 ) then
	  write(6,*) 'conz: computing explicitly ',nlv
	end if

	do k=1,nkn
	 ilevel = ilhkv(k)
	 do l=1,ilevel
c	  if(cdiag(l,k).gt.0.) then
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
c	    if(cdiag(l,k).eq.0) goto 7
	    aux=1./(cdiag(l,k)-clow(l,k)*chigh(l-1,k))
	    chigh(l,k)=chigh(l,k)*aux
	    cn(l,k)=(cn(l,k)-clow(l,k)*cn(l-1,k))*aux
	  end do
    7	  lstart=l-2
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

	return
	end

c*****************************************************************

	subroutine massconc(mode,cn,nlvdi,mass)

c computes total mass of conc

	implicit none

c arguments
	integer mode
	integer nlvdi
	real cn(nlvdi,1)
	real mass
c parameter
        include 'param.h'
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
	integer ilhkv(1)
	common /ilhkv/ilhkv
c local
	integer k,l,lmax
        double precision vol
	double precision sum,masstot
	real volnode

	if(nlvdim.ne.nlvdi) then
	  stop 'error stop : level dimension in massconc'
	end if

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

c*****************************************************************

        subroutine conzstab(cn1,co1
     +			,ddt
     +                  ,rkpar,difhv,difv
     +			,difmol,azpar
     +			,adpar,aapar
     +                  ,sindex
     +			,istot,isact
     +			,nlvdi,nlv)
c
c checks stability
c
c cn     new concentration
c co     old concentration
c caux   aux vector
c clow	 lower diagonal of vertical system
c chig	 upper diagonal of vertical system
c ddt    time step
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
c nlvdi	 dimension in z direction
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

	implicit none
c
c parameters
        include 'param.h'
c arguments
	integer nlvdi,nlv
        real cn1(nlvdi,1),co1(nlvdi,1)		!DPGGU
        real difv(0:nlvdi,1)
        real difhv(nlvdi,1)
	real difmol
        real ddt,rkpar,azpar,adpar,aapar			!$$azpar
	integer istot,isact
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi

	integer nen3v(3,1)
	common /nen3v/nen3v
	include 'ev.h'
        real utlnv(nlvdim,1),vtlnv(nlvdim,1)
        common /utlnv/utlnv, /vtlnv/vtlnv
        real utlov(nlvdim,1),vtlov(nlvdim,1)
        common /utlov/utlov, /vtlov/vtlov
        real wprv(0:nlvdim,1)
        common /wprv/wprv
        real wlov(0:nlvdim,1),wlnv(0:nlvdim,1)
        common /wlov/wlov, /wlnv/wlnv
	integer ilhv(1), ilhkv(1)
	common /ilhv/ilhv, /ilhkv/ilhkv

        real zeov(3,1),zenv(3,1)
        common /zeov/zeov, /zenv/zenv

        real hdknv(nlvdim,1)
        common /hdknv/hdknv
        real hdkov(nlvdim,1)
        common /hdkov/hdkov
 
        real hdenv(nlvdim,1)
        common /hdenv/hdenv
        real hdeov(nlvdim,1)
        common /hdeov/hdeov
 
        real areakv(nlvdim,1)
        common /areakv/areakv

        real wdifhv(3,3,1)
        common /wdifhv/wdifhv

	real saux1(nlvdim,nkndim)
        common /saux1/saux1

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
	double precision b(3),c(3),f(3)
	double precision rso,rsn,rsot,rsnt,rstot
	double precision hn,ho

	double precision difabs,difrel,volold,volnew,flxin,flxtot,diff
	double precision stabind,stabadv,stabdiff,stabvert,voltot

c------------------------------------------------------------
c big arrays
c------------------------------------------------------------
	double precision cn(nlvdim,nkndim)		!DPGGU	!FIXME
	double precision co(nlvdim,nkndim)
	double precision cdiag(nlvdim,nkndim)
	double precision clow(nlvdim,nkndim)
	double precision chigh(nlvdim,nkndim)
        real cwrite(nlvdim,nkndim)
c------------------------------------------------------------
c end of big arrays
c------------------------------------------------------------

	double precision cdummy
	double precision fw(3),fd(3)
	double precision fl(3)
c local (new)
	double precision clc(nlvdim,3), clm(nlvdim,3), clp(nlvdim,3)
	!double precision cl(0:nlvdim+1,3)
	double precision wl(0:nlvdim+1,3)
c
	double precision hdv(0:nlvdim+1)
	double precision haver(0:nlvdim+1)
	double precision hnew(0:nlvdim+1,3)
	double precision hold(0:nlvdim+1,3)
	double precision present(0:nlvdim+1)

        integer istab

        integer iustab
        save iustab
        data iustab /0/
c functions
	real getpar
	include 'testbndo.h'

	!write(6,*) 'conzstab called...'

        if(nlvdim.ne.nlvdi) stop 'error stop conzstab: level dimension'

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
                write(6,*) istot,isact,nlvdi,nlv
                write(6,*) nkn,nel
        end if

	az=azpar		!$$azpar
	azt=1.-az
	ad=adpar
	adt=1.-ad
	aa=aapar
	aat=1.-aa

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

        do l=0,nlvdim+1
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

        do l=1,nlvdim
	  do ii=1,3
	    clc(l,ii) = 0.
	    clm(l,ii) = 0.
	    clp(l,ii) = 0.
	  end do
	end do

c	-----------------------------------------------------------------
c	vertical velocities
c	-----------------------------------------------------------------

	do k=1,nkn
	  do l=0,nlv
!	    wprv(l,k) =  	(                       !Malta
!     +				   az*wlnv(l,k) 
!     +				+ azt*wlov(l,k) 
!     +				)
	    wprv(l,k) = wlnv(l,k)	!BUGFIX -> see definition in sp256w
	  end do
	end do

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
	    wl(l,ii) = wprv(l,k)
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
c	  wl(0,ii) = wprv(0,kn(ii))
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

	  w = wl(l-1,ii)		!top of layer
	  if( w .gt. 0. ) then          !out
	    fw(ii) = aat*w
	    clc(l,ii) = clc(l,ii) + aa*w
	  else
	    fw(ii) = 0.
	    clm(l,ii) = clm(l,ii) + aa*w
	  end if

	  w = wl(l,ii)			!bottom of layer
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
          cdummy = dt * aj4 * rk3 * hmed * ( b(ii)**2 + c(ii)**2 )
	  clow(l,k) = clow(l,k) + cdummy
          cdummy = dt * aj4 * 3. * f(ii)
          if( cdummy .lt. 0. ) then             !flux out of node
	    chigh(l,k) = chigh(l,k) - cdummy
          end if
          cn(l,k) = cn(l,k) + dt * aj4 * ( fw(ii) + fd(ii) )
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
c-----------------------------------------------------------------

        stabind = 0.
        stabadv = 0.
        stabdiff = 0.
        stabvert = 0.
        istab = 0
	do k=1,nkn
	  ilevel = ilhkv(k)
          if( is_inner(k) ) then
	   do l=1,ilevel
            voltot = cdiag(l,k)
            flxtot = chigh(l,k) + clow(l,k) + cn(l,k)
            if( voltot .gt. 0. ) then
                  aux = flxtot / voltot
                  if( aux .gt. stabind ) istab = k
		!  if( aux .gt. 100 ) then
		!	write(6,*) 'conzstab: ',k,l,aux
		!	write(6,*) voltot,flxtot,chigh(l,k),clow(l,k),cn(l,k)
		!  end if
                  stabind = max(stabind,aux)
		  cwrite(l,k) = aux		!save for write
                  aux = chigh(l,k) / voltot
                  stabadv = max(stabadv,aux)
		  saux1(l,k) = aux		!for adv. stab.
		!if( aux .gt. 0. ) write(6,*) 'kkkk: ',k,l,aux
                  aux = clow(l,k) / voltot
                  stabdiff = max(stabdiff,aux)
                  aux = cn(l,k) / voltot
                  stabvert = max(stabvert,aux)
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
c        call check2Dr(nlvdi,nlv,nkn,cwrite,0.,0.,"NaN check","cstab")

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

        !if( mod(niter,6) .eq. 0 ) then
        !if( it .gt. -70000 ) then
        if( .false. ) then
          call conwrite(iustab,'.stb',1,777,nlvdi,cwrite)
        end if

c        call stb_histo(it,nlvdi,nkn,ilhkv,cwrite)

        !write(93,*) 'conzstab1: ',it,istot,stabind,stabind/istot,istab
        !write(93,*) 'conzstab2: ',it,stabind,stabadv,stabdiff,stabvert
        !write(93,*) 'conzstab3: ',it,rkpar

c        write(94,*) it,istot,stabind

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
	end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c
c
c compute_stability	main routine -> does all the work (internal)
c get_stability		gets stab index or computes it
c set_stability		sets stab index, but does not use istot
c reset_stability	rests all istot values
c advect_stability	as get_stability but with rkpar=0
c
c typical usage:
c
c call reset_stability at beginning of time loop
c call set_stability before advection of similar variables
c call get_stability before real advection to get istot
c
c*****************************************************************

	subroutine compute_stability(mode,dt,rkpar,rindex,istot)

c computes stability index

	implicit none

        include 'param.h'

	integer mode		!0: reset  1:set/get_stability
        real dt
        real rkpar
        real rindex
        integer istot

	integer ndim
	parameter (ndim = 50)

	real rk(ndim)
	real rind(ndim)
	save rk,rind

        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real cnv(nlvdim,1)
        common /cnv/cnv
        real difv(0:nlvdim,1)
        common /difv/difv
        real difhv(nlvdim,1)
	common /difhv/difhv

        real saux(nlvdim,nkndim)

        real azpar,adpar,aapar
        real difmol
	real ddt
        integer isact,i
	logical debug

        real getpar

	integer nentry
	save nentry
	data nentry / 0 /

	debug = .true.
	debug = .false.

c----------------------------------------------------------------
c reset values
c----------------------------------------------------------------

	if( mode .eq. 0 ) then
	  nentry = 0
	  if( debug ) write(66,*) 'stab init: ',mode,rkpar
	  return
	end if

c----------------------------------------------------------------
c find entry with same rkpar
c----------------------------------------------------------------

	do i=1,nentry
	  if( rk(i) .eq. rkpar ) goto 1
	end do
    1	continue

c----------------------------------------------------------------
c rkpar found -> return indices
c----------------------------------------------------------------

	!i = nentry+1	!for debug purposes - should be commented

	if( i .le. nentry ) then
	  rindex = dt * rind(i)
	  istot = 1 + rindex
	  if( debug ) write(66,*) 'stab found: ',rkpar,i,istot,rindex
	  return
	end if

c----------------------------------------------------------------
c rkpar not found -> must compute
c----------------------------------------------------------------

	call getaz(azpar)
	adpar=getpar('adpar')
	aapar=getpar('aapar')

        isact = 1
        istot = 1
        difmol = 0.
	ddt = 1.		!always for 1 sec
	rindex = 0.

        call conzstab(cnv,saux
     +          ,ddt,rkpar,difhv,difv
     +		,difmol,azpar,adpar,aapar
     +          ,rindex,istot,isact,nlvdi,nlv)

c----------------------------------------------------------------
c insert new values into arrays
c----------------------------------------------------------------

	nentry = nentry + 1
	if( nentry .gt. ndim ) goto 99

	rk(nentry) = rkpar
	rind(nentry) = rindex

	rindex = dt * rindex
	istot = 1 + rindex

	if( debug ) write(66,*) 'stab comp: ',rkpar,nentry,istot,rindex

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	return
   99	continue
	write(6,*) 'nentry,ndim: ',nentry,ndim
	stop 'error stop compute_stability: ndim'
	end

c*****************************************************************

        subroutine get_stability(dt,rkpar,rindex,istot)

c gets stability index (if necessary computes it)

        implicit none

        real dt
        real rkpar
        real rindex
        integer istot

	integer mode

	mode = 1

	call compute_stability(mode,dt,rkpar,rindex,istot)

        end

c*****************************************************************

        subroutine set_stability(dt,rkpar)

c sets stability index (if necessary computes it, not interested in istot)

        implicit none

        real dt
        real rkpar

        real rindex
        integer istot
	integer mode

	mode = 1

	call compute_stability(mode,dt,rkpar,rindex,istot)

        end

c*****************************************************************

        subroutine reset_stability

c rests stability index

        implicit none

	integer mode,istot
	real dt,rkpar,rindex

	mode = 0
	dt = 0.
	rkpar = 0.

	call compute_stability(mode,dt,rkpar,rindex,istot)

        end

c**********************************************************************

	subroutine init_stability

c initializes stability computations
c
c this routine should be called after the hydro step
c the stability index is valid also for the next time step of hydro

        implicit none

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real dt,rkpar

	integer idtorg
	save idtorg
	data idtorg / 0 /

	if( idtorg .le. 0 ) idtorg = idt

        call reset_stability

	dt = idtorg
	rkpar = 0.		!only for advection

        !call set_stability(dt,rkpar)

	end
	
c**********************************************************************

        subroutine advect_stability(dt,rindex,istot)

c computes stability index for advection only

        implicit none

	include 'param.h'

        real dt
        real rindex
        integer istot

        integer nlv,nlvdi
        common /level/ nlvdi,nlv
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real saux1(nlvdim,nkndim)
	common /saux1/saux1
	real saux2(nlvdim,nkndim)
	common /saux2/saux2

	logical debug
	integer mode
	integer k,l,ilin
        real rkpar,ahpar
	real dindex,aindex,tindex

	real getpar

	debug = .true.
	debug = .false.
	mode = 0
        rkpar = 0.
	ahpar = getpar('ahpar')
	ilin = nint(getpar('ilin'))

	dindex = 0.
	aindex = 0.
	tindex = 0.

	do k=1,nkn
	  do l=1,nlv
	    saux1(l,k) = 0.	!this is set in compute_stability
	  end do
	end do

	if( debug ) write(66,*) ' computing advect... '
	call compute_stability(0,dt,rkpar,rindex,istot)
	call compute_stability(1,dt,rkpar,rindex,istot)
	if( debug ) write(66,*) 'advect... ',rkpar,0,istot,rindex

	!if( ahpar .gt. 0. ) then
	  do k=1,nkn
	    do l=1,nlv
	      saux2(l,k) = 0.
	    end do
	  end do
	  call viscous_stability(ahpar,saux2)
	  do k=1,nkn
	    do l=1,nlv
	      aindex = max(aindex,dt*saux1(l,k))	!for 1 sec
	      dindex = max(dindex,saux2(l,k))
	      tindex = max(tindex,dt*saux1(l,k)+saux2(l,k))
	    end do
	  end do
	!end if

	!write(6,*) 'advect_stability: ',rindex,aindex,dindex,tindex

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine stb_histo(it,nlvdi,nkn,ilhkv,cwrite)

        implicit none

        integer it
        integer nlvdi,nkn
        integer ilhkv(1)
        real cwrite(nlvdi,1)

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

        call histo_final(nbin,ic)

        write(98,*) it,ic

        end

c*****************************************************************

        subroutine grad_tvd(cc,gx,gy,aux,nlvdi,nlv)

c computes gradients for scalar cc

        implicit none

	integer nlvdi,nlv
	real cc(nlvdi,1)
	real gx(nlvdi,1)
	real gy(nlvdi,1)
	real aux(nlvdi,1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        
        integer nen3v(3,1)
        common /nen3v/nen3v        
	integer ilhv(1), ilhkv(1)
	common /ilhv/ilhv, /ilhkv/ilhkv
	include 'ev.h'
        
        integer k,l,ie,ii,lmax
	real b,c,area
	real ggx,ggy

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    gx(l,k) = 0.
	    gy(l,k) = 0.
	    aux(l,k) = 0.
	  end do
	end do

        do ie=1,nel
          area=ev(10,ie) 
	  lmax = ilhv(ie)
	  do l=1,lmax
            ggx=0
            ggy=0
            do ii=1,3
              k=nen3v(ii,ie)
              b=ev(ii+3,ie)
              c=ev(ii+6,ie)
              ggx=ggx+cc(l,k)*b
              ggy=ggy+cc(l,k)*c
              aux(l,k)=aux(l,k)+area
	    end do
            do ii=1,3
             k=nen3v(ii,ie)
             gx(l,k)=gx(l,k)+ggx*area
             gy(l,k)=gy(l,k)+ggy*area
            end do 
          end do
        end do

        do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    area = aux(l,k)
	    if( area .gt. 0. ) then
	      gx(l,k) = gx(l,k) / area
	      gy(l,k) = gy(l,k) / area
	    end if
	  end do
        end do

        end
        
c*****************************************************************

        subroutine grad_2d(cc,gx,gy,aux)

c computes gradients for scalar cc

        implicit none

	real cc(1)
	real gx(1)
	real gy(1)
	real aux(1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        
        integer nen3v(3,1)
        common /nen3v/nen3v        
	include 'ev.h'
        
        integer k,ie,ii
	real b,c,area
	real ggx,ggy

	do k=1,nkn
	  gx(k) = 0.
	  gy(k) = 0.
	  aux(k) = 0.
	end do

        do ie=1,nel
          area=ev(10,ie) 
          ggx=0
          ggy=0
          do ii=1,3
              k=nen3v(ii,ie)
              b=ev(ii+3,ie)
              c=ev(ii+6,ie)
              ggx=ggx+cc(k)*b
              ggy=ggy+cc(k)*c
              aux(k)=aux(k)+area
	  end do
          do ii=1,3
             k=nen3v(ii,ie)
             gx(k)=gx(k)+ggx*area
             gy(k)=gy(k)+ggy*area
          end do 
        end do

        do k=1,nkn
	    area = aux(k)
	    if( area .gt. 0. ) then
	      gx(k) = gx(k) / area
	      gy(k) = gy(k) / area
	    end if
        end do

        end
        
c*****************************************************************

	subroutine write_elem_info(iunit,ie)

	implicit none

	integer iunit
	integer ie

	include 'param.h'

	integer nen3v(3,1)
	common /nen3v/nen3v
	integer ilhv(1)
	common /ilhv/ilhv
	integer ilhkv(1)
	common /ilhkv/ilhkv
	real hev(neldim)
	common /hev/hev
	real hlhv(neldim)
	common /hlhv/hlhv

        real ulnv(nlvdim,1)
        common /ulnv/ulnv
        real vlnv(nlvdim,1)
        common /vlnv/vlnv
        real wlnv(0:nlvdim,1)
        common /wlnv/wlnv

        real utlnv(nlvdim,1)
        common /utlnv/utlnv
        real vtlnv(nlvdim,1)
        common /vtlnv/vtlnv

        real saltv(nlvdim,nkndim)
        common /saltv/saltv
        real tempv(nlvdim,nkndim)
        common /tempv/tempv

	integer iu,ii,k,l,lmax,lkmax

	iu = iunit
	if( iu .le. 0 ) iu = 6

	lmax = ilhv(ie)

	write(iu,*) '--------- write_elem_info ---------'
	write(iu,*) ie,lmax,hev(ie),hlhv(ie)
	write(iu,*) (nen3v(ii,ie),ii=1,3)

	write(iu,*) 'ut: ',(utlnv(l,ie),l=1,lmax)
	write(iu,*) 'vt: ',(vtlnv(l,ie),l=1,lmax)
	write(iu,*) 'u: ',(ulnv(l,ie),l=1,lmax)
	write(iu,*) 'v: ',(vlnv(l,ie),l=1,lmax)


	do ii=1,3
	  k = nen3v(ii,ie)
	  lkmax = ilhkv(k)
	  write(iu,*) 't: ',k,lkmax,(tempv(l,k),l=1,lkmax)
	  write(iu,*) 's: ',k,lkmax,(saltv(l,k),l=1,lkmax)
	  write(iu,*) 'w: ',k,lkmax,(wlnv(l,k),l=0,lkmax)
	end do
	  
	write(iu,*) '-----------------------------------'

	end

c*****************************************************************

