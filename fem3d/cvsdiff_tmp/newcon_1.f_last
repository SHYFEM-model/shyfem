c**************************************************************

        subroutine conz3d_1(cn1,co1
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

c	double precision explh(nlvdim,nlkdim)

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
	    wprv(l,k) =  	(
     +				   az*wlnv(l,k) 
     +				+ azt*wlov(l,k) 
     +				)
	  end do
	end do

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

c	  -----------------------------------------------
c	  contributions from vertical advection
c	  -----------------------------------------------
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

c	  -----------------------------------------------
c	  contributions from vertical diffusion
c	  -----------------------------------------------
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
	end do

c	  -----------------------------------------------
c contribution from horizontal advection (only explicit)
c	  -----------------------------------------------
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
	    clow(l,k)  = clow(l,k)  + aj4 * dt * clm(l,ii)
	    chigh(l,k) = chigh(l,k) + aj4 * dt * clp(l,ii)
	    cdiag(l,k) = cdiag(l,k) + aj4 * dt * clc(l,ii)
	    cdiag(l,k) = cdiag(l,k) + aj4 * hnew(l,ii)
	  end do
	end do

	end do		! loop over ie

c------------------------------------------------------------------
c end of loop over elements
c------------------------------------------------------------------

c in cdiag, chigh, clow is matrix (implicit part)
c if explicit calculation, chigh=clow=0 and in cdiag is volume of node [m**3]
c in cnv is mass of node [kg]
c for explicit treatment, cnv/cdiag gives new concentration [kg/m**3]

c integrate boundary conditions

	do k=1,nkn
	  ilevel = ilhkv(k)
	  do l=1,ilevel
            mflux = cbound(l,k)
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
