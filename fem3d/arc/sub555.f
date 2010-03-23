c
c $Id: sub555.f,v 1.28 2003/03/25 14:08:55 georg Exp $
c
c finite element asembling routines
c
c contents :
c
c subroutine sp156b(rmat,v1v,vqv,izeta,isum)	assembles linear system matrix
c subroutine sp159b				administrates one time step
c subroutine dryz(rmat,v1)			level estimation in dry areas
c subroutine nodav(v1)				nodal average of velocities
c subroutine masscs(ftot,ztot)			test mass conservation
c subroutine el2nod(aev,anv,zenv,vv,znull)	nodal average of element value
c
c revision log :
c
c 24.03.1993	ggu	!$$velpr - bug - factor 3 was missing
c 20.03.1998	ggu	masscs writes depending on debug level
c 25.03.1998	ggu	integrated changes from technital (minor)
c 29.04.1998	ggu	uses module for semi-implicit time-step
c 04.05.1998	ggu	integrated close and delwaq
c 25.05.1998	ggu	$$stress - wind in stress
c 17.06.1998	ggu	no iamat, isym
c 13.07.1998	ggu	completely restructured
c 21.07.1998	ggu	ireib = 4 introduced
c 22.07.1998	ggu	documentation for friction
c 03.08.1998	ggu	uvzchk introduced to stop unrealistic computation
c 03.08.1998	ggu	$$aust - austausch term was wrong (wrong sign and
c			   missing depth (dudx etc are velocities!!!)
c 03.08.1998	ggu	non linear advective term introduced
c 13.08.1998	ggu	no dudxv... -> now use uediff...
c 13.08.1998	ggu	nodav only for nodal average -> hor. diff. in hdiff
c 13.08.1998	ggu	new routines upstrm, uvzchk, hdiff
c 13.08.1998	ggu	coriolis term implemented
c 19.08.1998	ggu	use radv to compute advective terms (bug fix)
c 20.08.1998    ggu     iextpo finally eliminated
c 21.08.1998    ggu     xv eliminated
c 03.09.1998    ggu     some changes to avoid compiler warnings
c 04.10.1998    ggu     uvzchk -> increased error logging - continuing
c 20.04.1999    ggu     converted to stress instead of wind (tauxnv...)
c 14.10.1999    ggu     central discretization of advective terms (bcenter)
c 14.10.1999    ggu     in mass conservation control compute advective index
c
c notes :
c
c (1/rho_0) tau_x  ==  cd * (rho_air/rho_0) * |w| w_x
c
c stress in kg / ( m s**2 )
c
c rho_0		density water
c rho_air	density air
c cd		drag coefficient
c
c wk = cd * (rho_air/rho_0)
c
c wind   ==>          wc = wk * |w|
c stress ==>          wc = 1 / rho_0
c
c*********************************************************************
c
        subroutine sp156b(rmat,v1v,vqv,izeta,isum)
c
c assembles linear system matrix
c
c rmat		band matrix to be assembled
c v1v		izeta=0	: constant vector for u
c		izeta=1 : constant vector for z
c vqv		flux boundary condition vector
c izeta		0 : set up system for solving for u/v
c		1 : set up system for solving for z
c		2 : set up system for correction step for z
c isum		0 : set up only constant vector
c		1 : set up matrix and constant vector
c
c two-level-semi-implicit-scheme
c
c new formulation in transports
c
c written 04.10.90 by ggu  (from scratch)
c
c revised ...07.92 by ggu   $$lump  - lumping of matrix
c revised 31.08.92 by ggu   $$impli - implicit time step
c revised 16.03.93 by ggu   $$iclin - linear conti
c revised 19.01.94 by ggu   $$ireib0 - ireib=0 -> no friction
c revised 25.02.94 by ggu   $$wind - wind mass conserving (new definition
c					of gamma1 without factor 2)
c
	implicit none
c
c arguments
	integer izeta,isum
c        double precision rmat(1)
c        double precision v1v(1),v2v(1)
c        double precision vqv(1)
        real rmat(1)
        real v1v(1)
        real vqv(1)
c parameters
        real azz1,azz2,acc1,acc2
	real drittl,ar
	parameter (drittl=1./3.)
        parameter (azz1=1.0,azz2=1.-azz1)   !$$impli - momentum
        parameter (acc1=1.0,acc2=1.-acc1)   !$$impli - conti
	parameter (ar=1.)	!weighting of friction term
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real eps1,eps2,pi,flag,high,higi
	real grav,fcor,dcor,dirn,rowass,roluft
	integer itanf,itend,idt,nits,niter,it
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /mkonst/ eps1,eps2,pi,flag,high,higi
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
	common /femtim/ itanf,itend,idt,nits,niter,it
        integer nen3v(3,1),inodv(1)
        real hm3v(3,1)
	real ev(13,1), rzv(1)
        real czv(1)
	real tauxnv(1),tauynv(1)
	integer iwegv(1)
	real uov(1),vov(1),unv(1),vnv(1)
        real zeov(3,1)
        common /zeov/zeov
	common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv
        common /nen3v/nen3v, /inodv/inodv
        common /hm3v/hm3v
	common /ev/ev
	common /rzv/rzv
	common /tauxnv/tauxnv, /tauynv/tauynv
        common /czv/czv, /iwegv/iwegv
	real uedif(1), vedif(1)
	common /uedif/uedif, /vedif/vedif

        real znv(1)
        common /znv/znv
        real up0v(1), vp0v(1)
        common /up0v/up0v, /vp0v/vp0v

        real ruv(1)
        real rvv(1)
        common /ruv/ruv
        common /rvv/rvv

        real rdistv(1)
        common /rdistv/rdistv

c local
	logical bcenter		!central discretization of advective terms
c	logical bdebug
c	logical isebnd
	integer kn(3),kk
        integer ie,ii,i,j,j1,j2,n,m,k
	integer ireib,iclin
        real am,amt,az
        real czdef,hzoff,rlamb,rr
	real dt,aj,area,hzg,gcz,rw
	real dtdt,delta,gamma,lambda,aa
	real g,xx,yy
	real zm,hm,hhm
	real taux,tauy
	real cuv,cxy
	real bz,cz,bbcc,zbbcc
        real xaus,yaus
c	real dudx,dudy,dvdx,dvdy
        real xadv,yadv,radv
	real uc,vc,um,vm,u,v
	real f,fm
	real uso,vso
	real b(3),c(3),h(3),amatr(3,3)
	real hia(3,3),hik(3)
	real z(3),zneu(3)
	real aweigh
	real umom, vmom
	real ucor, vcor
        real rstab,rstabmax
        real rdist
c function
        integer iround,locsps
	real getpar
        logical iskbnd,iseout
c       logical iskobd
c       logical iskout	!$$ALPHA
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
c       iskobd(k) = inodv(k).gt.0 
c       iskout(k) = inodv(k).eq.-2	!$$ALPHA
        iseout(ie) = iwegv(ie).ne.0

c only line to change to pass to lumped matrix is definition of amatr

        data amatr / 4.,0.,0.,0.,4.,0.,0.,0.,4. /   !$$lump
c	data amatr / 2.,1.,1.,1.,2.,1.,1.,1.,2. /

c discretization of advective terms (central does not work...)

	bcenter = .true.
	bcenter = .false.

c variables from global section

	czdef=getpar('czdef')			! friction coefficient
        hzoff=getpar('hzoff')   		!minimum depth
	ireib=iround(getpar('ireib'))		!type of friction
	iclin=iround(getpar('iclin'))		!conti linear	!$$iclin
	radv=getpar('ilin')			!advective terms

c weighting of time-steps

	call getazam(az,am)
        amt=1.-am

c constants

	rlamb = czdef

	if( radv .lt. 0. ) radv = 0.
	if( radv .gt. 1. ) radv = 1.
	radv = 1. - radv

	dt=idt
	dtdt=dt*dt
	g=grav

	do ie=1,nel

	aj=ev(10,ie)	!area of triangle / 12
	area = 12. * aj

	uso=uov(ie)
	vso=vov(ie)

	zm=0.
	hm=0.
	taux=0.
	tauy=0.
	umom=0.
	vmom=0.
        rdist = 0.

	do ii=1,3
	  k=nen3v(ii,ie)
	  kn(ii)=k

	  b(ii)=ev(ii+3,ie)
	  c(ii)=ev(ii+6,ie)

	  h(ii)=hm3v(ii,ie)
          z(ii)=zeov(ii,ie)
	  zneu(ii)=znv(k)
	  zm=zm+z(ii)
	  hm=hm+h(ii)

	  taux=taux+tauxnv(k)
	  tauy=tauy+tauynv(k)

	  umom = umom + ruv(k)
	  vmom = vmom + rvv(k)

          rdist = rdist + rdistv(k)
	end do

        rdist = rdist / 3.

	taux=taux*drittl
	tauy=tauy*drittl

	umom = umom / (3.*area)
	vmom = vmom / (3.*area)

	hhm=(hm+zm)*drittl
	if(iclin.ne.0) hhm=hm*drittl		!$$iclin

c advective terms

	uc = uso / hhm
	vc = vso / hhm

        xadv = 0.
        yadv = 0.
	fm = 0.

	do ii=1,3

	  k = kn(ii)

	  u = up0v(k)
	  v = vp0v(k)

	  f = uso * b(ii) + vso * c(ii)		! f > 0. ==> flux into node

	  if( bcenter ) then			!central discretization
c	    xadv = xadv + uso * u * b(ii) + vso * u * c(ii)
c	    yadv = yadv + uso * v * b(ii) + vso * v * c(ii)
	    xadv = xadv - uso * ( u * b(ii) + v * c(ii) )
	    yadv = yadv - vso * ( u * b(ii) + v * c(ii) )
	    fm = 1.
	  else if( f .lt. 0. ) then		!upstream velocity
	    um = 1.5 * b(ii) * ( u - uc )
	    vm = 1.5 * c(ii) * ( u - uc )
	    xadv = xadv + ( uso * um + vso * vm ) * f
	    um = 1.5 * b(ii) * ( v - vc )
	    vm = 1.5 * c(ii) * ( v - vc )
	    yadv = yadv + ( uso * um + vso * vm ) * f
	    fm = fm + f
	  end if

	end do

	if( fm .ne. 0. ) then
	  fm = rdist * radv / fm
	end if

	xadv = xadv * fm
	yadv = yadv * fm

c horizontal turbulent friction

	xaus = uedif(ie) * hhm			!$$aust; sign adjusted in xx/yy
	yaus = vedif(ie) * hhm

c coriolis

	ucor = - fcor * vso
	vcor =   fcor * uso

c-----------------------------------------------------------
c
c DOCS  START   P_friction
c
c DOCS  FRICTION		Bottom friction
c
c The friction term in the momentum equations can be written as
c $Ru$ and $Rv$ where $R$ is the variable friction coefficient and
c $u,v$ are the velocities in $x,y$ direction respectively.
c The form of $R$ can be specified in various ways. The value of 
c |ireib| is choosing between the formulations. In the parameter
c input file a value $\lambda$ is specified that is used in 
c the formulas below.
c
c |ireib|	Type of friction used (default 0):
c		\begin{description}
c		\item[0] No friction used
c		\item[1] $R=\lambda$ is constant
c		\item[2] $\lambda$ is the Strickler coefficient.
c			 In this formulation $R$ is written as
c			 $R = \frac{g}{C^2} \frac{\vert u \vert}{H}$
c			 with $C=k_s H^{1/6}$ and $\lambda=k_s$ is
c			 the Strickler coefficient. In the above
c			 formula $g$ is the gravitational acceleration,
c			 $\vert u \vert$ the modulus of the current velocity
c			 and $H$ the total water depth.
c		\item[3] $\lambda$ is the Chezy coefficient.
c			 In this formulation $R$ is written as
c			 $R = \frac{g}{C^2} \frac{\vert u \vert}{H}$
c			 and $\lambda=C$ is the Chezy coefficient.
c		\item[4] $R=\lambda/H$ with $H$ the total water depth
c		\item[5] $R=\lambda\frac{\vert u \vert}{H}$
c		\end{description}
c |czdef|	The default value for the friction parameter $\lambda$.
c		Depending on the value of |ireib| the coefficient $\lambda$
c		is describing linear friction or a Chezy or Strickler
c		form of friction (default 0).
c |iczv|	Normally $R$ is evaluated at every time step (|iczv| = 1).
c		If for some reason this behavior is not desirable,
c		|iczv| = 0 evaluates the value of $R$ only before the
c		first time step, keeping it constant for the
c		rest of the simulation. (default 1)
c
c The value of $\lambda$ may be specified for the whole basin through
c the value of |czdef|. For more control over the friction parameter
c it can be also specified in section |area| where the friction
c parameter depending on the type of the element may be varied. Please
c see the paragraph on section |area| for more information.
c
c DOCS  END
c
c-----------------------------------------------------------
c       hzg & hzoff new defined !!!!!!!!
c		--> hzg should be always non linear

	hzg=(hm+zm)*drittl	!$$iclin
        if(hzg.lt.hzoff) hzg=hzoff

	if(ireib.eq.0) then	!$$ireib0
		rr=0.
        else if(ireib.eq.1) then
                rr = rlamb
	else if(ireib.eq.2) then
		gcz=grav/((czv(ie)**2)*(hzg**drittl))	!??
		rr=gcz*sqrt(uso*uso+vso*vso)/(hzg*hzg)	!??
	else if(ireib.eq.3) then
		gcz=grav/(czv(ie)**2)			!??
		rr=gcz*sqrt(uso*uso+vso*vso)/(hzg*hzg)	!??
        else if(ireib.eq.4) then
                rr = rlamb/hzg
	else if(ireib.eq.5) then
		rr=rlamb*sqrt(uso*uso+vso*vso)/(hzg*hzg)	!??
	else
		write(6,*) 'unknown friction : ',ireib
		stop 'error stop : sp156'
	end if

c	if( ie .eq. 1960 ) then
c	  write(98,*) ie,ireib,rlamb,rr,hzg,sqrt(uso*uso+vso*vso)
c	end if

c        rstab = (1.+dtdt*fcor*fcor)/(1.+dt*rr)**2
c        rstabmax = max(rstabmax,rstab)
c        if( ie .eq. 100 .or. ie .eq. 500 .and. ie .eq. -1 ) then
c                write(98,*) fcor,rr,dtdt*fcor*fcor,(1.+dt*rr)**2,rstab
c        end if

c analytic expression for friction term and coriolis
c should stabilize the solution with explicit coriolis term
c (still to be verified)

        rr = ( sqrt(1.+dtdt*fcor*fcor) * exp(dt*rr) - 1. ) / dt

c-----------------------------------------------------------

	delta= 1. / ( 1. + dt * rr * ar )
	gamma=( 1. - dt * rr * (1.-ar) ) * delta

c explicit terms (are identical for z or u/v computation)

	xx = xadv - taux - xaus - umom + ucor
	yy = yadv - tauy - yaus - vmom + vcor

	if( ie .eq. -1 ) then
	  write(6,*) xadv,taux,xaus,umom,ucor
	  write(6,*) yadv,tauy,yaus,vmom,vcor
	end if

	if (izeta.eq.1) then			!water levels

	  cuv = area * dt * ( 1. + (gamma-1.) * az )
	  cxy = area * dtdt * delta * az
          lambda = area * dtdt * delta * g * hhm * az

	  do n=1,3
	    zbbcc=0.
	    do m=1,3
	      aa = aj * amatr(n,m)
	      bbcc = lambda * ( b(n)*b(m) + c(n)*c(m) )
              hia(n,m) = aa + bbcc * am
              zbbcc = zbbcc + ( aa - bbcc*amt ) * z(m)
	    end do
            hik(n) = zbbcc
     +                    + cuv * ( uso*b(n) + vso*c(n) )
     +                    - cxy * (  xx*b(n) +  yy*c(n) )
	  end do

c	  level boundary conditions

	  do i=1,3
	    rw=rzv(kn(i))
	    if(rw.ne.flag) then
              j1=mod(i,3)+1
              j2=mod(i+1,3)+1
              hik(j1)=hik(j1)-rw*hia(j1,i)
              hik(j2)=hik(j2)-rw*hia(j2,i)
              hia(i,j1)=0.
              hia(i,j2)=0.
              hia(j1,i)=0.
              hia(j2,i)=0.
              hik(i)=rw*hia(i,i)
	    end if
	  end do

c	  excluded areas

          if( iseout(ie) ) then

            do n=1,3
              do m=1,3
                hia(n,m)=area*( b(n)*b(m) + c(n)*c(m) )
              end do
              hik(n)=0.
            end do

            do n=1,3
              if( iskbnd(kn(n)) ) then
                do m=1,3
                  hia(n,m)=0.
                  hia(m,n)=0.
                end do
                hik(n)=0.
              end if
            end do

	  end if

	else if(izeta.eq.0) then	!velocities

		bz=0.
		cz=0.
		do i=1,3
            	  k=nen3v(i,ie)
            	  bz=bz+b(i)*(amt*z(i)+am*zneu(i))
            	  cz=cz+c(i)*(amt*z(i)+am*zneu(i))
		end do

          	unv(ie) = gamma*uso - dt*delta * ( g*hhm*bz + xx )
          	vnv(ie) = gamma*vso - dt*delta * ( g*hhm*cz + yy )

	else

	  stop 'error stop sp156b: internal error (1)'

	end if

c in hia(i,j),hik(i),i,j=1,3 is system

        if(izeta.eq.1) then
          do i=1,3
          if(isum.ne.0) then
            do j=1,3
              kk=locsps(kn(i),kn(j),nkn,mbw)
              if(kk.gt.0) then
                rmat(kk)=rmat(kk)+hia(i,j)
              end if
            end do
          end if
          v1v(kn(i))=v1v(kn(i))+hik(i)
          end do
        end if

	end do !do ie

	if(izeta.eq.1) then
          do i=1,nkn
            v1v(i)=v1v(i)+vqv(i)
	  end do
	end if

c	if(izeta.eq.1) then
c            write(98,*) 'sub555 (adv. stability) : ',it,rstabmax
c        end if

	end

c*********************************************************

	subroutine sp159b

c administrates one time step for system to solve
c
c revised on 27.07.88 by ggu   (introduction of ic)
c revised on 17.11.88 by ggu   (v6v,sp156f)
c revised on 19.12.88 by ggu   (no ieend in sp156f)
c revised on 14.03.90 by ggu   (sp136 at end and out of iteration, set ic=0)
c revised on 09.10.90 by ggu   (new sp158k)
c revised    ...07.92 by ggu   $$lump   - lumping of matrix
c revised    05.09.92 by ggu   $$close1 - closing at beginning
c revised    11.01.94 by ggu   $$zuniq - zuniq has been commented
c revised    14.01.94 by ggu   $$weg - call to setweg preserves "old" iw
c revised    17.01.94 by ggu   $$com - seldom used subroutines commented
c revised    17.01.94 by ggu   $$zeov - zeov introduced for concentration
c revised    20.01.94 by ggu   $$zov0 - zeov in conz adjusted
c revised    03.02.94 by ggu   $$zov00 - zeov adjusted here (use with masscs)
c revised    17.03.95 by ggu   $$tecn395 - changes for tecnital

        implicit none

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

        real amat(1)
        common /amat/amat
        real rqv(1)
        common /rqv/rqv

        real zov(1), znv(1)
        common /zov/zov, /znv/znv
        real zenv(3,1),zeov(3,1)	!$$zeov $$zov00
        common /zenv/zenv, /zeov/zeov
	real uov(1),vov(1),unv(1),vnv(1)
	common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv

        real v1v(1),v2v(1)
        common /v1v/v1v, /v2v/v2v
        real v3v(1),v4v(1)
        common /v3v/v3v, /v4v/v4v
        real v5v(1)
        common /v5v/v5v
c-------------------!FIXME
        real rzv(1)
        common /rzv/rzv
c-------------------!FIXME

c local
c	logical bdebug
        integer nrand,iterad,isumz,iw,iwa
        integer ie,i,ii,iadv,icount,ier,isym
	integer n,k
        real epseps
c functions
        real getpar
c save
	save iw

        nrand=getpar('nrand')   !new matrix every nrand iteration
        iterad=getpar('iterad') !total number of non-linear iterations
	if( iterad .le. 0 ) iterad = 1

	iw=0
        call sp136(iw)      !$$close1	$$com

        call setweg(2,iwa)   !$$weg
	iw=iw+iwa
        call setweg(3,iwa)
	iw=iw+iwa

c        call zuniq(v1v,v2v) !$$lump - make z values unique $$zuniq $$tecn395

        isumz=0                             !1 ==> sum to a-matrix
        if(niter.eq.1) then                 !first iteration
          isumz=1                           !...sum to both
        else if(nrand.ne.0) then            !sum to a-matrix every
		if(mod(niter,nrand).eq.0) isumz=1	!...
        end if                              !...nrand iteration
        if(iw.ne.0) then                    !change of area
          isumz=1                           !...sum to both
	end if

	do k=1,nkn
	  zov(k) = znv(k)
	end do

	do i=1,nel
	  uov(i)=unv(i)
	  vov(i)=vnv(i)
	end do

	do ie=1,nel	!$$zeov $$zov0 $$zov00
	  do ii=1,3
	    zeov(ii,ie)=zenv(ii,ie)
	  end do
	end do

	iadv=0			!counts advect. iterations
	icount=0		!counts iterations/timestep

c iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do while(iadv.lt.iterad)

          if(isumz.eq.1) call newlnk

c z-values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	  isym = 1	!1=symmetric  2=diagonal    0=any

	  n = nkn*((2-isym)*mbw+1)
          do i=1,n
            amat(i)=0.
          end do

          do ii=1,nkn
            v1v(ii)=0.
          end do

          call sp156b(amat,v1v,rqv,1,isumz)

          epseps=1.e-6

          call mchb(v1v,amat,nkn,1,mbw,1,epseps,ier)

          if(ier.ne.0) goto 99

          do k=1,nkn
	    znv(k)  = v1v(k)
          end do

c controll intertidal flats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          call setweg(1,iw)

          if(iw.gt.0) then  !domain has changed
            isumz=1
            iadv=0
          else              !compute u/v
            isumz=0
            iadv=iadv+1
            call sp156b(amat,v1v,rqv,0,0)
          end if

          icount=icount+1

	end do	!do while...iteration

c iteration end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call dryz(amat,v1v) !z-value in dry areas	!$$tecn395

        call setzev
        call setuvd
        call arper
        call nodav(v1v)
	call setxv
        call hdiff(v1v,v2v,v3v,v4v,v5v)

	call uvzchk		!checks for strange u/v/z values
	call masscs(v1v,v2v)

        call delwaq	!$$com

	return
   99   continue
        write(6,*) ier
        stop 'error stop sp159 : inversion of matrix'
	end
c
c****************************************
c
        subroutine dryz(rmat,v1)
c
c estimation of levels in dry areas
c
c rmat          band matrix already decomposed
c v1            auxiliary vector, is used by routine
c               ...to assemble constant vector
c
        implicit none
c
c arguments
        real rmat(1),v1(1)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1),iwegv(1),inodv(1)
        real ev(13,1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nen3v/nen3v, /iwegv/iwegv, /inodv/inodv
        common /ev/ev
	real znv(1)
	common /znv/znv
c local
        integer ii,ie,k,ii1,ii2,kk1,kk2,ier
        real epseps,z,aj
        real b(3),c(3)
c functions
        logical iskbnd,iskout,iseout
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
        iskout(k) = inodv(k).eq.-2
        iseout(ie) = iwegv(ie).ne.0
c
        do k=1,nkn
          v1(k)=0.
	end do
c
	do ie=1,nel
          if( iseout(ie) ) then
            do ii=1,3
              b(ii)=ev(ii+3,ie)
              c(ii)=ev(ii+6,ie)
            end do
            aj=12.*ev(10,ie)
            do ii=1,3
              k=nen3v(ii,ie)
              if( iskbnd(k) ) then
                z=znv(k)
                ii1=mod(ii,3)+1
                ii2=mod(ii1,3)+1
                kk1=nen3v(ii1,ie)
                kk2=nen3v(ii2,ie)
                v1(kk1)=v1(kk1)-aj*(b(ii)*b(ii1)+
     +                      c(ii)*c(ii1))*z
                v1(kk2)=v1(kk2)-aj*(b(ii)*b(ii2)+
     +                      c(ii)*c(ii2))*z
                v1(k)=0.
              end if
            end do
          end if
	end do
c
        epseps=1.e-6
        call mchb(v1,rmat,nkn,1,mbw,-1,epseps,ier)
        if(ier.ne.0) goto 99
c
        do k=1,nkn
          if( iskout(k) ) znv(k)=v1(k)
	end do
c
	return
   99   continue
        write(6,*) 'ier from sp158s : ',ier
        stop 'error stop sp158s'
	end
c
c****************************************
c
        subroutine nodav(v1)
c
c nodal average of velocities
c
c v1            auxiliary vector, is used by routine
c               ...to assemble weighting factor
c
c revision log :
c
c 31.07.1997	ggu	$$clin - introduced to handle linear conti
c 12.08.1998	ggu	no first derivatives here

        implicit none
c
c arguments
        real v1(1)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1),iwegv(1)
        real ev(13,1),zenv(3,1),hm3v(3,1)
        real unv(1),vnv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nen3v/nen3v, /iwegv/iwegv
        common /ev/ev, /zenv/zenv, /hm3v/hm3v
        common /unv/unv, /vnv/vnv
	real up0v(1), vp0v(1)
	common /up0v/up0v, /vp0v/vp0v

c local
	logical bcolin
        integer i,ie,k
        real aj,aux,hzg,zm
c	real uns,vns
c functions
	integer iround
	real getpar
        logical iseins
        iseins(ie) = iwegv(ie).eq.0

	bcolin=iround(getpar('iclin')).ne.0

        do k=1,nkn
          v1(k) = 0.
	  up0v(k) = 0.
	  vp0v(k) = 0.
        end do

        do ie=1,nel
          if( iseins(ie) ) then
            hzg=0.
	    zm=0.
            do i=1,3
              hzg=hzg+hm3v(i,ie)
              zm=zm+zenv(i,ie)
            end do
	    if( .not. bcolin ) hzg = hzg + zm	!$$clin
            aj=ev(10,ie)
            aux=3.*aj/hzg	!$$velpr - error corrected 24.03.93
            do i=1,3
              k=nen3v(i,ie)
              v1(k)=v1(k)+aj
              up0v(k)=up0v(k)+aux*unv(ie)
              vp0v(k)=vp0v(k)+aux*vnv(ie)
            end do
          end if
        end do

        do k=1,nkn
          if( v1(k) .gt. 0. ) then
	    aux = 1. / v1(k)
            up0v(k) = up0v(k) * aux
            vp0v(k) = vp0v(k) * aux
          else
            up0v(k) = 0.
            vp0v(k) = 0.
          end if
        end do

        end

c*****************************************************************

        subroutine hdiff(area,dudx,dudy,dvdx,dvdy)

c horizontal diffusion

        implicit none

c arguments (auxiliary arrays)
        real area(1)
	real dudx(1), dudy(1)
	real dvdx(1), dvdy(1)
c parameters
	real third
	parameter ( third = 1. / 3. )
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1),iwegv(1)
        real ev(13,1),zenv(3,1),hm3v(3,1)
        real unv(1),vnv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nen3v/nen3v, /iwegv/iwegv
        common /ev/ev, /zenv/zenv, /hm3v/hm3v
        common /unv/unv, /vnv/vnv
	real uedif(1), vedif(1)
	common /uedif/uedif, /vedif/vedif
        real up0v(1), vp0v(1)
        common /up0v/up0v, /vp0v/vp0v
c local
	logical bcolin
        integer i,ie,k
        real aj,aux,hzg,zm,uns,vns
	real akh,weight
	real du,dv
c	real u,v
c functions
	integer iround
	real getpar
        logical iseins
	real b,c
        iseins(ie) = iwegv(ie).eq.0
	b(i,ie) = ev(3+i,ie)
	c(i,ie) = ev(6+i,ie)

	bcolin = iround(getpar('iclin')) .ne. 0
	akh = getpar('akh')

c initialize arrays

        do k=1,nkn
	  area(k) = 0.
          dudx(k)=0.
          dudy(k)=0.
          dvdx(k)=0.
          dvdy(k)=0.
        end do

c accumulate contributions to nodal average

        do ie=1,nel
          if( iseins(ie) ) then

	    hzg = 0.
	    zm = 0.
            do i=1,3
	      k = nen3v(i,ie)
              hzg=hzg+hm3v(i,ie)
              zm=zm+zenv(i,ie)
            end do
	    if( .not. bcolin ) hzg = hzg + zm	!$$clin
	    aux = 3. / hzg

	    uns = unv(ie) * aux
	    vns = vnv(ie) * aux

            aj = ev(10,ie)
	    weight = 1.5 * aj

            do i=1,3
              k=nen3v(i,ie)
	      du = up0v(k) - uns
	      dv = vp0v(k) - vns
	      area(k) = area(k) + weight
	      dudx(k) = dudx(k) + weight * b(i,ie) * du
	      dudy(k) = dudy(k) + weight * c(i,ie) * du
	      dvdx(k) = dvdx(k) + weight * b(i,ie) * dv
	      dvdy(k) = dvdy(k) + weight * c(i,ie) * dv
            end do

          end if
        end do

c finish nodal average

        do k=1,nkn
          if( area(k) .gt. 0. ) then
	    aux = 1. / area(k)
            dudx(k) = dudx(k) * aux
            dudy(k) = dudy(k) * aux
            dvdx(k) = dvdx(k) * aux
            dvdy(k) = dvdy(k) * aux
          else
            dudx(k) = 0.
            dudy(k) = 0.
            dvdx(k) = 0.
            dvdy(k) = 0.
          end if
        end do

c compute diffusion contribution

	do ie=1,nel
	  du = 0.
	  dv = 0.

	  do i=1,3
	    k = nen3v(i,ie)
	    du = du + dudx(k) * b(i,ie) + dudy(k) * c(i,ie)
	    dv = dv + dvdx(k) * b(i,ie) + dvdy(k) * c(i,ie)
	  end do

	  uedif(ie) = akh * du
	  vedif(ie) = akh * dv
	end do

        end

c*****************************************************************

	subroutine masscs(ftot,ztot)

c test mass conservation

	implicit none

c arguments
	real ftot(1),ztot(1)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /femtim/ itanf,itend,idt,nits,niter,it
        integer nen3v(3,1),inodv(1)
        real uov(1),vov(1),unv(1),vnv(1)
        real zenv(3,1),zeov(3,1)
        real ev(13,1)
	real hm3v(3,1)
        common /nen3v/nen3v, /inodv/inodv
        common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv
        common /zenv/zenv, /zeov/zeov
        common /ev/ev
        common /hm3v/hm3v
c local
	integer ie,ii,k
	real az,aj,dt,fo,fn,ft,dz,dif
	real difmx,advmax
        double precision zztot
c functions
	integer iround
	real getpar
c save
	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
	  if( iround(getpar('levdbg')) .le. 3 ) icall = -1
	  if( icall .eq. -1 ) return
	  icall = 1
	end if

	call getaz(az)

	dt=idt

        zztot=0.
	do k=1,nkn
	  ftot(k)=0.
	  ztot(k)=0.
	end do

	do ie=1,nel
	  do ii=1,3
		k=nen3v(ii,ie)
		aj=ev(10,ie)
		fo=ev(ii+3,ie)*uov(ie)+ev(ii+6,ie)*vov(ie)
		fn=ev(ii+3,ie)*unv(ie)+ev(ii+6,ie)*vnv(ie)
		ft=12.*aj*dt*(az*fn+(1.-az)*fo)
		dz=4.*aj*(zenv(ii,ie)-zeov(ii,ie))
		ftot(k)=ftot(k)+ft
		ztot(k)=ztot(k)+dz
                zztot = zztot + 4.*aj*zenv(ii,ie)
	  end do
	end do

c in difmx will be maximum difference of water volume (in cubicmeters)
c ...due to roundoff errors
c in order to compute relative error we must divide by
c ...the total depth (not done yet)

	difmx=0.
	do k=1,nkn
	  if(inodv(k).le.0) then	!no open boundary node
	    dif=abs(ftot(k)-ztot(k))
	    if(dif.gt.difmx) difmx=dif
	  end if
	end do

c do again but compute flux/volume  -------------------------------

	do k=1,nkn
	  ftot(k)=0.
	  ztot(k)=0.
	end do

	do ie=1,nel
	  do ii=1,3
		k=nen3v(ii,ie)
		aj=ev(10,ie)
		fn=ev(ii+3,ie)*unv(ie)+ev(ii+6,ie)*vnv(ie)
		ft=12.*aj*dt*fn
		if( ft .lt. 0.) ft = -ft
		dz=4.*aj*(zenv(ii,ie)+hm3v(ii,ie))
		ftot(k)=ftot(k)+ft
		ztot(k)=ztot(k)+dz
	  end do
	end do

c in advmax will be maximum .. between flux and volume
c it must be smaller than 1. for advective terms

	advmax=0.
	do k=1,nkn
	    dif=0.5*ftot(k)/ztot(k)
	    if(dif.gt.advmax) advmax=dif
	end do

	write(6,*) 'masscs : ',it,zztot,difmx,advmax

	end

c****************************************

        subroutine el2nod(aev,anv,zenv,vv,znull)

c nodal average of element value
c
c aev		element array					(in)
c anv		nodal value					(return)
c zenv		water level in element if dry, else znull	(in)
c vv            auxiliary vector, is used by routine		(aux)
c znull		flag for wet elements				(in)

        implicit none

c arguments
	real aev(1), anv(1)
        real vv(1), zenv(1)
	real znull
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1)
        common /nen3v/nen3v
        real ev(13,1)
        common /ev/ev
c local
        integer i,ie,k
        real aj
c functions
        logical iseins
        iseins(ie) = zenv(ie).eq.znull

	call testev

        do k=1,nkn
          vv(k)=0.
        end do

        do ie=1,nel
          if( iseins(ie) ) then
            aj=ev(10,ie)
            do i=1,3
              k=nen3v(i,ie)
              vv(k)=vv(k)+aj
              anv(k)=anv(k)+aj*aev(ie)
            end do
          end if
        end do

        do k=1,nkn
          if(vv(k).gt.0.) then
            anv(k)=anv(k)/vv(k)
            anv(k)=anv(k)/vv(k)
          else
            anv(k)=0.
            anv(k)=0.
          end if
        end do

	return
	end

c****************************************

        subroutine uvzchk

c check of velocity and level data

        implicit none

c parameters
	real uvmax,zmax
	parameter( uvmax = 10. , zmax = 100. )
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real znv(1)
        common /znv/znv
        real up0v(1), vp0v(1)
        common /up0v/up0v, /vp0v/vp0v
	real unv(1), vnv(1)
	common /unv/unv, /vnv/vnv
	real zenv(3,1)
	common /zenv/zenv
	integer nen3v(3,1)
	common /nen3v/nen3v
	integer iwegv(1)
	common /iwegv/iwegv
	real hm3v(3,1)
	common /hm3v/hm3v
c local
	integer k,ie,ii
	logical bwrite
	real u,v,z,uv
	real zm,hm,htot
c function
	integer ipext,ieext

c	call uvmaxi(uv,ie)
c	write(6,*) 'maximum velocity : ',ieext(ie),uv
c	call eleinfo(6724)

	do k=1,nkn
	  z = abs( znv(k) )
	  u = abs( up0v(k) )
	  v = abs( vp0v(k) )
	  if( z .gt. zmax ) goto 99
	  if( u .gt. uvmax .or. v .gt. uvmax ) goto 99
	end do

	return
   99	continue
	write(6,*) 'Unrealistic high values computed in node ',ipext(k)
	write(6,*) 'u/v/z : ',u,v,z

	do ie=1,nel
	  bwrite = .false.
	  do ii=1,3
	    if( nen3v(ii,ie) .eq. k ) bwrite = .true.
	  end do
	  if( bwrite ) then
	    call eleinfo(ie)
	  end if
	end do

	return
c	stop 'error stop uvzchk'

	end

c****************************************************************

        subroutine eleinfo(ie)

        implicit none

c arguments
	integer ie
c parameters
	real drittl
	parameter(drittl=1./3.)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real unv(1), vnv(1)
	common /unv/unv, /vnv/vnv
	real hm3v(3,1)
	common /hm3v/hm3v
        real znv(1)
        common /znv/znv
        real zov(1)
        common /zov/zov
	real zenv(3,1)
	common /zenv/zenv
	real zeov(3,1)
	common /zeov/zeov
	integer nen3v(3,1)
	common /nen3v/nen3v
	integer iwegv(1)
	common /iwegv/iwegv
c local
	integer ii
	real zm,hm,htot
	logical bfull

	integer ieext

	bfull = .true.

	zm = 0.
	hm = 0.
	do ii=1,3
	  zm = zm + zenv(ii,ie)
	  hm = hm + hm3v(ii,ie)
	end do
	zm = zm / 3.
	hm = hm / 3.
	htot = hm + zm

	write(6,*) ' eleinfo : ',ieext(ie),iwegv(ie),htot
	write(6,*) 'u-v ',unv(ie),vnv(ie),unv(ie)/htot,vnv(ie)/htot

	if( bfull ) then
	  write(6,*) 'h   ',(hm3v(ii,ie),ii=1,3)
	  write(6,*) 'zn  ',(zenv(ii,ie),ii=1,3)
	  write(6,*) 'zo  ',(zeov(ii,ie),ii=1,3)
	  write(6,*) 'zkn ',(znv(nen3v(ii,ie)),ii=1,3)
	  write(6,*) 'zko ',(zov(nen3v(ii,ie)),ii=1,3)
	end if

	end

c******************************************************************

        subroutine uvmaxi(uvmax,iemax)

c computes maximum velocity in element and element where found

        implicit none

c arguments
	real uvmax
	integer iemax
c parameters
	real drittl
	parameter(drittl=1./3.)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real unv(1), vnv(1)
	common /unv/unv, /vnv/vnv
	real hm3v(3,1)
	common /hm3v/hm3v
	real zenv(3,1)
	common /zenv/zenv
c local
        integer ie,ii
	real hm,zm,htot,rh
	real u,v,uv

	uvmax = 0.

	do ie=1,nel
	  u = unv(ie)
	  v = vnv(ie)

	  zm = 0.
	  hm = 0.
	  do ii=1,3
	    zm = zm + zenv(ii,ie)
	    hm = hm + hm3v(ii,ie)
	  end do
	  htot = ( hm + zm ) * drittl
	  rh = 1. / htot

	  u = u * rh
	  v = v * rh

	  uv = u*u + v*v

	  if( uv .gt. uvmax ) then
	    uvmax = uv
	    iemax = ie
	  end if
	end do

	uvmax = sqrt( uvmax )

	end

c******************************************************************

        subroutine setxv

c sets obsolete data structure xv

        implicit none

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real up0v(1),vp0v(1)
        common /up0v/up0v, /vp0v/vp0v
        real znv(1)
        common /znv/znv
        real xv(3,1)
        common /xv/xv
c local
        integer k

        do k=1,nkn
          xv(1,k) = up0v(k)
          xv(2,k) = vp0v(k)
          xv(3,k) = znv(k)
        end do

        end

c******************************************************************

	subroutine watvol(kn,vol)

c computes water volume in finite volume kn

	implicit none

	integer kn
	real vol

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v
	real ev(13,1)
	common /ev/ev
	real hm3v(3,1)
	common /hm3v/hm3v
	real zenv(3,1)
	common /zenv/zenv

	integer ie,ii,k
	real area,htot
	double precision sum

	sum = 0.

	do ie=1,nel
	  area = 4. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( kn .le. 0 .or. kn .eq. k ) then
	      htot = hm3v(ii,ie) + zenv(ii,ie)
	      sum = sum + htot * area
	    end if
	  end do
	end do

	vol = sum

        end

c******************************************************************

