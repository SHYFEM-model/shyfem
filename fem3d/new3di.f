c
c $Id: new3di.f,v 1.35 2010-03-11 15:36:38 georg Exp $
c
c assembling linear system routine
c
c contents :
c
c subroutine sp259f			administrates one time step
c subroutine sp256z(rmat,zv,vqv,isum)	assemble matrix
c subroutine sp256v			assembles vertical system matrix
c subroutine sp256w(vv)			computes vertical velocities
c
c notes :
c
c look for "!ccc" to see important changes
c
c ASYM			passage to asymmetrix matrix
c ASYM_OPSPLT		new version without operator splitting (can leave)
c ASYM_OPSPLT_CH	new version -> change here
c
c nkndim	dimension for total number of nodes
c neldim	dimension for total number of elements
c nrbdim	dimension for total number of boundary condition nodes
c nbcdim	dimension for total number of open boundaries
c mbwdim	dimension for bandwidth
c ngrdim	dimension for grade of nodes (number of elements attached
c		...to one node)
c nardim	dimension for total number of area codes
c nexdim	dimension for total number of extra nodes for output
c
c nkn,nel	total number of nodes/elements
c nrz,nrq	total number of nodes with water level/flux boundary conditions
c nrb		total number of nodes with boundary conditions
c nbc		total number of open boundaries
c ngr		maximum grade of nodes
c mbw		bandwidth of system matrix
c flag		flag value (to recognize boundary conditions)
c grav,dcor	gravitational accel./medium latitude of basin
c rowass,roluft	density of water/air
c itanf,itend	start/end time for simulation
c idt,nits	time step/total iterations to go
c niter,it	actual number of iterations/actual time
c nlvdi,nlv	dimension for levels, number of used levels
c
c nen3v(..,ie)	element index - node numbers (3) of element ie
c
c ipv(k)	external node number of node k
c ipev(ie)	external element number of element ie
c
c rqv(k),rzv(k)	flux/water level boundary conditions for node k
c		...if(rzv(k).eq.flag) --> no b.c. is specified for node k
c bnd(..,ib)	specification of boundary conditions for open boundary ib
c ev(..,ie)	geometric parameters for element ie
c		...1-3 = area, 4-6 = b, 7-9 = c, 10 = Aomega, 11-13 = angle
c uov(ie)	depth integrated transport (old time level)
c vov(ie)	...
c unv(ie)	depth integrated transport (new time level)
c vnv(ie)	...
c zov(k)	water level of old/new time level
c znv(k)	...
c ulov(l,ie)	velocity of old time level in x direction of layer l and elem ie
c vlov(l,ie)	velocity of old time level in y direction of layer l and elem ie
c wlov(l,i)	velocity of old time level in z direction of layer l and node i
c ulnv(l,ie)	velocity of new time level in x direction of layer l and elem ie
c vlnv(l,ie)	velocity of new time level in y direction of layer l and elem ie
c wlnv(l,i)	velocity of new time level in z direction of layer l and node i
c utlov(l,ie)	transport of old time level in x direction of layer l and el. ie
c vtlov(l,ie)	transport of old time level in y direction of layer l and el. ie
c utlnv(l,ie)	transport of new time level in x direction of layer l and el. ie
c vtlnv(l,ie)	transport of new time level in y direction of layer l and el. ie
c uprv(l,k)	velocity (averaged) in x direction of layer l and node k
c vprv(l,k)	velocity (averaged) in y direction of layer l and node k
c wprv(l,k)	velocity (averaged) in z direction of layer l and node k
c up0v(k)	total velocity (averaged) in x direction of node k
c vp0v(k)	total velocity (averaged) in y direction of node k
c tauxnv(k)	normalized stress in x-direction at node k
c tauynv(k)	normalized stress in y-direction at node k
c rhov(l,k)	density for level l and node k
c fcorv(ie)	coriolis parameter for elem ie
c
c visv(l,k)	vertical turbulent viscosity for layer l and node k (alv)
c difv(l,k)	vertical turbulent diffusivity for layer l and node k (slv)
c
c xgv(k),ygv(k)	coordinates of node k
c
c hldv(l)	thickness of layer l
c hlv(l)	absolute depth of bottom of layer l
c ilhv(ie)	number of levels for element ie
c ilhkv(k)	number of levels for node k
c hlhv(ie)	thickness of last layer for element ie
c hev(ie)	total depth at element ie (no water level)
c hkv(k)	total depth at node k (no water level)
c hm3v(3,ie)	depth of three nodes in element ie
c
c v1v,v2v...	auxiliary vectors
c rmat		band matrix for one vertical system (as above)
c rvec		constant vector for vertical system
c
c $$h1new	use new water level to compute velocities
c		(only for printing velocities important, but
c		algorithm has to be checked if consistent)
c $$rtmax	use maximal friction coefficient of rdt (=1./dt)
c
c revision log :
c
c revised 01.07.93 	$$UVBARO - u/vov introduced for	iteration on rad cond
c revised 03.11.93 	$$cmplerr - compiler warnings hydro
c revised 05.11.93 	$$fric - normal friction
c revised 05.11.93 	$$crador - crador call commented
c revised 05.11.93 	subroutine crador in file newcra.f
c revised 05.11.93 	$$VBARO-ERR - unv(ie)=vov(ie)
c revised 28.08.95 	$$BAROC_AREA - do baroc only for iarv(ie) = 0
c revised 30.08.95      $$AUST - austausch coefficient introduced
c revised 01.09.95      $$AWEIGH - area weighting of austausch coefficient
c revised 06.03.96 	$$BAROC_AREA0 - introduced baroc0
c revised 06.03.96 	$$VERT_AUST_ADJUST - adjustment of vert. aust. coef.
c revised 06.06.96 	$$BCHAO - modifications for vel. profile (temp.)
c revised 10.06.96 	$$UVPADV - modifications for advective term
c 14.08.1998	ggu	set w = 0 at open boundary nodes
c 20.08.1998	ggu	some documentation for sp256w
c 08.04.1999    ggu     equilibrium tide introduced (zeqv)
c 20.04.1999    ggu     converted to stress instead of wind (tauxnv...)
c 24.06.1999    ggu     call to rescur commented (use 2D call resid)
c 07.03.2000    ggu     eliminated VERT_AUST_ADJUST
c 07.03.2000    ggu     arrays for vertical turbulent diffusion coefficient
c 12.01.2001    ggu     solve for znv and not level difference (ZNEW)
c 09.11.2001    ggu     BCHAO commented out, no compiler directives
c 14.11.2001    ggu     all compiler directives eliminated
c 10.08.2003    ggu     deleted commented parts (radiation, etc..)
c 10.08.2003    ggu     some code in new subroutines: make_prvel, copy_uvz
c 13.08.2003    ggu     delete useless vars (nrand, epsit), no iteration
c 13.08.2003    ggu     call setnod, set_link_info
c 10.03.2004    ggu     RQVDT - value in rqv is now discharge [m**3/s]
c 10.01.2005    ggu     BAROC_DIST - scale baroclinic term with distance
c 25.01.2005    ggu     BUGADV - bugfix for advective terms
c 24.02.2005    ggu     new reynolds stresses -> green
c 15.03.2005    ggu     austv,aust eliminated, austau() in diff_h_set()
c 29.04.2005    ggu     semi-lagrangian advection (backadv)
c 29.04.2005    ggu     no baroclinic terms close to step (ilevmin)
c 04.11.2005    ggu     use itlin to decide about semi-lagrangian advection
c 23.03.2006    ggu     changed time step to real
c 31.05.2006    ggu     new friction type ireib=7
c 18.10.2006    ccf     radx,rady for radiation stress introduced
c 28.11.2006    ggu     in u/vadv is now difference and not absolute transport
c 02.04.2007    ggu     new algorithm (look for ASYM)
c 08.06.2007    ggu&deb restructured for new explicit terms
c 28.09.2007    ggu	deleted chao, semi-lagrange to newexp.f, no indov
c 24.06.2008    ggu	bpresv deleted
c 10.10.2008	ggu&mbj	prepared for pardiso -> modify system (gguexclude)
c 10.12.2008    ggu	use rfricv for bottom friction -> other can be deleted
c 18.12.2008    ggu	more debug info for error in vertical system
c 13.01.2009    ggu	Pardiso lib integrated
c 04.03.2009    ggu	matrix amat deleted from file -> only locally used
c 27.03.2009    ggu	call new routine adjust_mass_flux() for dry nodes
c 06.04.2009    ggu	deleted routine write_elem_vel_info()
c 07.05.2009    ggu	new routines for scalar interpolation (not finished)
c 10.03.2010    ggu	bug fix in sp256w() for ibtyp=2
c 11.03.2010    ggu	new routine check_volume() to check for negative vol
c 12.04.2010    ggu	ad hoc routine for Yaron
c 16.12.2010    ggu	in sp256w() account for changing volume (sigma)
c 19.02.2011    ccf	3D radiation stress
c
c******************************************************************

	subroutine sp259f
c
c administrates one time step for system to solve
c
c delta with itbas : 1 overhead,initialize  2 boundary line  3 assemble
c			4 preconditioning  5 choleski  6 solve
c			7 uv bound. cond.  8 dry areas
c
c written on 27.07.88 by ggu   (from sp159f)
c
	implicit none
c
c new :::::: nlv
c
c parameter
	include 'param.h'
c common
	integer itanf,itend,idt,nits,niter,it
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /level/ nlvdi,nlv
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real rqv(1),v1v(1),v2v(1)
	real zov(1),znv(1),unv(1),vnv(1)
	real ulnv(nlvdim,1),vlnv(nlvdim,1)
	real ulov(nlvdim,1),vlov(nlvdim,1)
	real wlov(0:nlvdim,1),wlnv(0:nlvdim,1)
	real utlnv(nlvdim,1),vtlnv(nlvdim,1)
	real utlov(nlvdim,1),vtlov(nlvdim,1)
	real uov(1),vov(1)			!$$UVBARO
	real wprv(0:nlvdim,1)			!$$VERVEL
        real saux1(nlvdim,1),saux2(nlvdim,1)

	common /rqv/rqv
	common /v1v/v1v, /v2v/v2v
	common /zov/zov, /znv/znv
	common /unv/unv, /vnv/vnv
	common /ulnv/ulnv, /vlnv/vlnv
	common /ulov/ulov, /vlov/vlov
	common /wlov/wlov, /wlnv/wlnv
	common /utlnv/utlnv, /vtlnv/vtlnv
	common /utlov/utlov, /vtlov/vtlov
	common /uov/uov, /vov/vov	!$$UVBARO
	common /wprv/wprv			!$$VERVEL
        common /saux1/saux1, /saux2/saux2
	real zeov(3,1),zenv(3,1)
	common /zeov/zeov, /zenv/zenv
	real v3v(1)
	common /v3v/v3v

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

cccccccccccccccccccccccccccc
	real vp0v(1),vprv(nlvdim,1)
	common /vp0v/vp0v , /vprv/vprv
cccccccccccccccccccccccccccc
ccccccccccccc
ccccccccccccc
c	real hm3v(3,1),crad(1)
c	common /hm3v/hm3v, /crad/crad
c	real hhh
c	integer ii
c	integer ielist(5)
c	real xlist(5)
c	data ielist /77,163,698,746,1346/
ccccccccccccc
ccccccccccccc
c local
	integer i,l,k,ie,ier,ii
	integer nrand
	integer iw,iwa
	integer nmat
	integer kspecial
c	integer ninf
	real res
	real epseps
c	real voltot,deptot
c	real kin,pot
cccccccccccccccccccccccccccc
c	real v1,v2(25)
c	integer i1,i2,i3
cccccccccccccccccccccccccccc
c function
c	integer iround,ideffi
	integer iround
	real getpar,resi
c save
	save epseps
c data
        data epseps / 1.e-6 /

	kspecial = 3878
	kspecial = 0

c constants...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension in sp259f'

c dry areas

	iw=0
	call sp136(iw)

        call setweg(2,iwa)   !$$weg	!ZEONV
        iw=iw+iwa
        call setweg(3,iwa)
        iw=iw+iwa

	call copy_uvz		!copies uvz to old time level

	call copydepth(nlvdim,hdknv,hdkov,hdenv,hdeov)
	call setdepth(nlvdim,hdkov,hdeov,zeov,areakv)

c austauch contribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call diff_h_set

c solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call sp256v	!FIXME -> this should prob. go into loop

	call bclfix

c---------------------------------------------------------- z solution

	iw = 1
	do while( iw .gt. 0 )	!loop over changing domain

	  call setnod
	  call set_link_info
	  call adjust_mass_flux		!cope with dry nodes

	  call system_init

	  call sp256z(rqv)

	  call system_solve_z

	  call system_adjust_z

	  call setweg(1,iw)	!controll intertidal flats

	end do	!do while( iw .gt. 0 )

c---------------------------------------------------------- end of z solution

	!call check_node(kspecial)

	call sp256n	!new velocities (also barotropic)

c end of solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setzev     !znv -> zenv
        call setuvd

	call baro2l 
	call setdepth(nlvdim,hdknv,hdenv,zenv,areakv) !only now zenv ready
	call check_volume			!checks for negative volume 
        call arper

	res=resi(zov,znv,nkn)

c w-values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call sp256w(saux1,saux2)	!$$VERVEL

c compute velocities from transports %%%%%%%%%%%%%%%%%%%%%%%%%

	call ttov

c check mass balance (v1v is dummy)

	call mass_conserve(saux1,saux2)

c compute nodal values for velocities %%%%%%%%%%%%%%%%%%%%%%%%

	call make_prvel

	return
   99	continue
	write(6,*) 'Error in inverting matrix for water level'
	write(6,*) 'it, ier : ',it,ier
	stop 'error stop : sp259f'
	end

c******************************************************************

	subroutine sp256z(vqv)
c
c assembles linear system matrix
c
c vqv		flux boundary condition vector
c
c semi-implicit scheme for 3d model
c
c written on 18.02.91 by ggu  (from scratch)
c changed on 04.06.91 by ggu  (c=(1) : friction term has been corrected)
c changed on 01.10.92 by ggu  (staggered FE - completely restructured)
c 12.01.2001    ggu     solve for znv and not level difference (ZNEW)
c
	implicit none
c
c parameters
	include 'param.h'
c arguments
	real vqv(1)
c parameters
c	real az,azt,am,amt,af,aft
c	parameter (az=0.50,azt=1.-az)
c	parameter (am=0.50,amt=1.-am)
c	parameter (af=0.50,aft=1.-af)
	real drittl
	parameter (drittl=1./3.)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real eps1,eps2,pi,flag,high
	real grav,fcor,dcor,dirn,rowass,roluft
	integer itanf,itend,idt,nits,niter,it
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /mkonst/ eps1,eps2,pi,flag,high
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
	common /femtim/ itanf,itend,idt,nits,niter,it
	integer nen3v(3,1),ilhv(1)
	integer iarv(1)
	real rzv(1)
	real hev(1)
	real zov(1),unv(1),vnv(1)
	real utlov(nlvdim,1),vtlov(nlvdim,1)
	real utlnv(nlvdim,1),vtlnv(nlvdim,1)
	common /nen3v/nen3v, /ilhv/ilhv
	common /iarv/iarv
	common /rzv/rzv
	common /hev/hev
	common /zov/zov, /unv/unv ,/vnv/vnv
	common /utlov/utlov, /vtlov/vtlov, /utlnv/utlnv, /vtlnv/vtlnv
        real zeov(3,1),zenv(3,1)
        common /zeov/zeov, /zenv/zenv
	integer iwegv(1), inodv(1)
	common /iwegv/iwegv, /inodv/inodv
	include 'ev.h'

        integer iuvfix(1)       !chao deb
        common /iuvfix/iuvfix   !chao deb
        integer afix            !chao deb

        double precision ddxv(2*nlvdim,neldim)  !ASYM
        double precision ddyv(2*nlvdim,neldim)
        common /ddxv/ddxv, /ddyv/ddyv
        save /ddxv/, /ddyv/

c local

	logical bcolin
c	logical debug,bdebug
	logical bdebug
	integer kn(3)
	integer ie,i,j,j1,j2,n,m,kk,l,k
	integer ncount,ngl
	integer ilevel
	integer ju,jv
c	real az,azt,am,amt,af,aft
	real az,am,af,azpar,ampar
	real dt,aj,rw
c	real um,vm,zm
	real zm
	real ut,vt,uhat,vhat
	real ht
	real h11,hh999
	real delta
c	real fcora,beta,gamma
	real hia(3,3),hik(3),amatr(3,3)
	real b(3),c(3),z(3)
	real acu
	real uold,vold
	real dbb,dbc,dcb,dcc,abn,acn
c	real rmin,ulr,umr,vlr,vmr
c
c	integer iradb(3),iradf
c function
        integer locsps,loclp,iround
	real getpar
	logical iskbnd,iskout,iseout
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
        iskout(k) = inodv(k).eq.-2
        iseout(ie) = iwegv(ie).ne.0
c save - data
c	data amatr / 2.,1.,1.,1.,2.,1.,1.,1.,2. /
	data amatr / 4.,0.,0.,0.,4.,0.,0.,0.,4. /
	data ncount /0/
c
	ncount=ncount+1
c
	ngl=nkn
c
c constants
c
	bcolin=iround(getpar('iclin')).ne.0
c
	call getazam(azpar,ampar)
	az=azpar
	am=ampar
	af=getpar('afpar')
	call get_timestep(dt)

	acu = 0.

	do ie=1,nel

	zm=0.
	do i=1,3
		kk=nen3v(i,ie)
		kn(i)=kk
c
		b(i)=ev(i+3,ie)
		c(i)=ev(i+6,ie)
		z(i)=zov(kk)
		z(i)=zeov(i,ie)		!ZEONV
		zm=zm+z(i)
	end do
c
	zm=zm*drittl
c
	if(bcolin) then
		ht=hev(ie)
	else
		ht=hev(ie)+zm
	end if

	ilevel=ilhv(ie)
	aj=ev(10,ie)
        afix=1-iuvfix(ie)      !chao deb
c
	!delta=dt*dt*az*am*grav*ht
	!delta=dt*dt*az*am*grav		!ASYM_OPSPLT
        delta=dt*dt*az*am*grav*afix         !ASYM_OPSPLT        !chao deb

	dbb = 0.
	dbc = 0.
	dcb = 0.
	dcc = 0.
	do l=1,ilevel			!ASYM_OPSPLT
	  jv=l+l
	  ju=jv-1
	  dbb = dbb + ddxv(ju,ie)
	  dbc = dbc + ddyv(ju,ie)
	  dcb = dcb + ddxv(jv,ie)
	  dcc = dcc + ddyv(jv,ie)
	end do

	uold = 0.
	vold = 0.
	uhat = 0.
	vhat = 0.

	do l=1,ilevel
	  uold = uold + utlov(l,ie)
	  vold = vold + vtlov(l,ie)
	  uhat = uhat + utlnv(l,ie)
	  vhat = vhat + vtlnv(l,ie)
	end do

	ut = az * uhat + (1.-az) * uold
	vt = az * vhat + (1.-az) * vold

	do n=1,3
	  do m=1,3
	    abn = b(n) * ( b(m) * dbb + c(m) * dbc )
	    acn = c(n) * ( b(m) * dcb + c(m) * dcc )
	    !abn = ht * b(n) * b(m)			!ASYM_OPSPLT_CH
	    !acn = ht * c(n) * c(m)			!ASYM_OPSPLT_CH
	    h11 = delta*( abn + acn )			!ASYM_OPSPLT_CH
	    !h11 = delta*( b(n)*b(m)+c(n)*c(m) )
	    hia(n,m) = aj * (amatr(n,m) + 12.*h11)
	  end do
	  acu = hia(n,1)*z(1) + hia(n,2)*z(2) + hia(n,3)*z(3)
	  hik(n) = acu + 12.*aj*dt*( ut*b(n) + vt*c(n) )	!ZNEW
	end do

c	level boundary conditions
c
c====================================================================
c
c	new implementation of radiation boundary condition $$GWI
c
c	iradf=0
c	do i=1,3
c	  if(rzv(kn(i)).ne.flag) then
c	    iradb(i)=1
c	    iradf=1
c	  else
c	    iradb(i)=0
c	  end if
c	end do
c
c changed 25.03.96 -> no reference to bnd anymore, must be handled
c	with seperate function
c
c	if(iradf.eq.1.and.kn(1).gt.2000.and.bnd(2,2).eq.52) then
c	  if(niter.le.5.and.ie.eq.1) write(6,*) 'radiation only for canaly'
c	  write(6,*) ie
c	write(6,*) '*******************'
c	write(6,*) zm,uhat,vhat,um,vm
c	write(6,*) z
c	  call gwi(hia,hik,iradb,b,c,z,dt,ht,aj) !implicit gravity wave rad.
c	else
c
	do i=1,3
	  if(rzv(kn(i)).ne.flag) then
		!rw=rzv(kn(i))-zov(kn(i))	!this for dz
		rw=rzv(kn(i))			!FIXME !ZNEW (this for znew)
		j1=mod(i,3)+1
		j2=mod(i+1,3)+1
		hik(j1)=hik(j1)-rw*hia(j1,i)
		hik(j2)=hik(j2)-rw*hia(j2,i)
		hia(i,j1)=0.
		hia(i,j2)=0.
		hia(j1,i)=0.
		hia(j2,i)=0.
c		hia(i,i)=12.*aj
		hik(i)=rw*hia(i,i)
	  end if
	end do
c
c	end if
c
c==================================================================

	bdebug = ie.eq.20.or.ie.eq.100.or.ie.eq.250
	bdebug = .false.
	if( bdebug ) then
	  write(6,*) 'xxx: it,ie,ht ',it,ie,ht
c          write(6,*) 'xxx: ht,g,dt ',ht,grav,dt
c          write(6,*) 'xxx: az,am ',az,am

	  do i=1,3
	    write(6,*) 'xxx: hia ',(hia(i,j),j=1,3)
	  end do
	  write(6,*) 'xxx: hik ',(hik(j),j=1,3)

	  write(6,*) 'xxx: u ',uhat,uold,ut
	  write(6,*) 'xxx: v ',vhat,vold,vt

c          write(6,*) 'xxx: b ',(ev(i+3,ie),i=1,3)
c          write(6,*) 'xxx: c ',(ev(i+6,ie),i=1,3)


	end if

c
c excluded areas
c
          if( iseout(ie) ) then	!ZEONV
c          if( .false. ) then
            hh999=aj*12.
            do n=1,3
              do m=1,3
                hia(n,m)=hh999*(b(n)*b(m)+c(n)*c(m))
              end do
              hik(n)=0.
            end do
c
            do n=1,3
              if( iskbnd(kn(n)) ) then	!not internal and not out of system
                do m=1,3
                  hia(n,m)=0.
c                  hia(m,n)=0.		!gguexclude - comment
                end do
                hik(n)=0.
              end if
            end do
          end if

c in hia(i,j),hik(i),i,j=1,3 is system

	  call system_assemble(nkn,mbw,kn,hia,hik)

	end do

	call system_add_rhs(dt,vqv)

	end

c******************************************************************

	subroutine sp256v

c assembles vertical system matrix
c
c semi-implicit scheme for 3d model
c
c written on 18.02.91 by ggu  (from scratch)
c
	implicit none
c
c parameters
	include 'param.h'
	real drittl
	parameter (drittl=1./3.)
	real reps
	parameter (reps=1.e-5)
	!parameter (reps=0.)			!no check for bottom friction
c	real az,azt,am,amt,af,aft,at,att
c	parameter (az=0.50,azt=1.-az)
c	parameter (am=0.50,amt=1.-am)
c	parameter (af=0.50,aft=1.-af)
c	parameter (at=0.50,att=1.-at)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real eps1,eps2,pi,flag,high
	real grav,fcor,dcor,dirn,rowass,roluft
	integer itanf,itend,idt,nits,niter,it
	integer nlvdi,nlv
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /mkonst/ eps1,eps2,pi,flag,high
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /level/ nlvdi,nlv
	integer nen3v(3,1)
	integer ilhv(1)
	integer iarv(1)
	integer ipev(1)
	real hev(1)
	real zov(1)
	real utlov(nlvdim,1),vtlov(nlvdim,1)
        real utlnv(nlvdim,1),vtlnv(nlvdim,1)
	real uprv(nlvdim,1),vprv(nlvdim,1)
        real saux1(nlvdim,1),saux2(nlvdim,1)
        real saux3(nlvdim,1),saux4(nlvdim,1)
	real tauxnv(1),tauynv(1)
	real unv(1),vnv(1)
	real fcorv(1)
	real visv(0:nlvdim,1)
	real czv(1)
	common /nen3v/nen3v
	common /ilhv/ilhv
	common /iarv/iarv
	common /ipev/ipev
	common /hev/hev
	common /zov/zov
	common /utlov/utlov, /vtlov/vtlov, /utlnv/utlnv, /vtlnv/vtlnv
	common /uprv/uprv, /vprv/vprv
        common /saux1/saux1, /saux2/saux2
        common /saux3/saux3, /saux4/saux4
	common /tauxnv/tauxnv, /tauynv/tauynv
	common /unv/unv, /vnv/vnv
	common /fcorv/fcorv 
	common /visv/visv
	common /czv/czv
	include 'ev.h'

       integer iuvfix(1)        !chao deb
       common /iuvfix/iuvfix    !chao deb
       integer afix             !chao deb

	integer ilmv(1)
	common /ilmv/ilmv

	real difhv(nlvdim,1)
	common /difhv/difhv
	real rfricv(1)
        common /rfricv/rfricv

        real hdeov(nlvdim,1)
        common /hdeov/hdeov

        real zeov(3,1),zenv(3,1)
        common /zeov/zeov, /zenv/zenv
	real ppv(1)
	common /ppv/ppv

        real zeqv(1)
        common /zeqv/zeqv

        real rdistv(1)
        common /rdistv/rdistv

c local

	logical bbaroc,barea0                  !$$BAROC_AREA0
	integer nbcount

	logical bfirst,blast,bcolin,baroc,baust
	logical debug,bdebug
        logical bdebggu
        logical bgreen
	integer kn(3)
	integer kk,ii,ie,l,ju,jv,ncount
	integer ngl,mbb
	integer ilevel,ier,ilevmin
	integer lp,lm
	integer ireib
	integer k1,k2,k3
	integer ibaroc
	real b(3),c(3)
	real aj12
	real az,am,af,at,azpar,ampar
	real czdef,rr,rt
	real rrho0,dtgrav
	real zz
	real ut,vt
	real up,vp
	real dt,hlh,rdt,hlh_new
c	real ht
c	real bz,cz,dbz,dcz,zm
	real bz,cz,zm,zmm
	real xbcl,ybcl
        real xexpl,yexpl
c	real beta
	real ulm,vlm,taux,tauy
	real gamma,gammat
        real hhi,hhim,hhip,uui,uuim,uuip,vvi,vvim,vvip
c	real bb,bbt,cc,cct,aa,aat,ppx,ppy,aux,aux1,aux2
	real bb,bbt,cc,cct,aa,aat,aux
        real epseps
	real ahpar,aust
	real fact                       !$$BCHAO - not used
	real uuadv,uvadv,vuadv,vvadv
        real rhp,rhm,aus
	real hzg,hzoff,gcz
	real vismol
        real xmin,xmax
        integer imin,imax
        real rdist
        real xadv,yadv,fm,uc,vc,f,um,vm
	real bpres,cpres
	real area,ugreen,vgreen
        
        real fxv(nlvdim,1)      !new HYDRO deb
        real fyv(nlvdim,1)
        common /fxv/fxv
        common /fyv/fyv
	real rraux,cdf

        real radx(nlvdim,neldim),rady(nlvdim,neldim)
        common /radx/radx,/rady/rady

	real ss
	real z0bk(nkndim)		!bottom roughenss on nodes
	common /z0bk/z0bk

c-----------------------------------------
	real hact(0:nlvdim+1)
	real rhact(0:nlvdim+1)
	real alev(0:nlvdim)
c-----------------------------------------
c	real rmat(1),rvec(1)
c	common /rmat/rmat, /rvec/rvec
c	real ppx,ppy
	double precision rmat(10*nlvdim)
	!double precision rvec(2*nlvdim)
	double precision rvec(6*nlvdim)		!ASYM
	double precision ppx,ppy

	double precision ddxv(2*nlvdim,neldim)	!ASYM
	double precision ddyv(2*nlvdim,neldim)
	common /ddxv/ddxv, /ddyv/ddyv
	save /ddxv/, /ddyv/

c-----------------------------------------
c	integer iaux
c	real*8 uaux,vaux
c	real uamax,vamax
c	real uuaux,vvaux,uvet1,uvet2
c-----------------------------------------
c=(3) aux variables (of no importance)
c	real ppaux,ppmin,ppmax
	real auxaux(-2:+2)
c function
	integer locssp,iround
	real getpar
c save
	save epseps,ncount
c data
	data ncount / 0 /
	data epseps / 1.e-6 /
c
	hact(0) = 0.
c
	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension in sp256v'
c
	ncount=ncount+1
c
c constants
c
	czdef=getpar('czdef')			! friction coefficient	
	bcolin=iround(getpar('iclin')).ne.0	! linearized conti
	ireib=iround(getpar('ireib'))		! friction term
        ahpar   = getpar('ahpar')		! new value for austausch
        hzoff  = getpar('hzoff')		! minimum depth
        vismol  = getpar('vismol')		! molecular viscosity
        baust  = ahpar .ne. 0.			! austausch contrib.
	ibaroc = iround(getpar('ibarcl'))	! baroclinic contributions
c	baroc = ibaroc .ge. 1 .or. ibaroc .le. 2	!BUG!!!
	baroc = ibaroc .eq. 1 .or. ibaroc .eq. 2
        barea0 = .false.                        ! baroclinic only with ia = 0
	bgreen = .false.			!use green for austausch

	call getazam(azpar,ampar)
	az=azpar
	am=ampar
	af=getpar('afpar')
	at=getpar('atpar')

	rr=czdef

	call get_timestep(dt)
    
	rdt=1./dt

	bdebug=.true.
	bdebug=.false.
	if( bdebug) write(6,*) 'debug: ',rr,dt
c
	dtgrav=dt*grav
	rrho0=1./rowass
	if( .not. baroc ) rrho0 = 0.
	ut=0.
	vt=0.
	up=0.
	vp=0.

	nbcount = 0

c-------new computation of explicit part----------------------------------

	call bottom_friction	!set bottom friction
        call set_explicit       !new HYDRO deb
	!call set_yaron

c-------result: arrays fxv(l,ie),fyv(l,ie)-----------------------------


c======================================================================
c======================================================================
c                       start of element loop
c======================================================================
c======================================================================

	do 500 ie=1,nel

	bdebug = ie.eq.20.or.ie.eq.100.or.ie.eq.250
	bdebug=.false.

        bbaroc = baroc
	if( barea0 ) then               !$$BAROC_AREA $$BAROC_AREA0
	  if( iarv(ie) .ne. 0 ) bbaroc = .false.
	  if( bbaroc ) nbcount = nbcount + 1
        end if

	rrho0=1./rowass
	if( .not. bbaroc ) rrho0 = 0.

c area

	aj12 = 12. * ev(10,ie)
	area = aj12

c new system

	ilevel=ilhv(ie)
	!ilevmin=ilmv(ie)
	ngl=2*ilevel
	mbb=2
	if(ngl.eq.2) mbb=1

c compute barotropic pressure term

	debug=.false.

	bz=0.
	cz=0.
	bpres=0.
	cpres=0.
	zm=0.
	zmm=0.
	taux=0.
	tauy=0.
        rdist = 0.
	do ii=1,3
	  kk=nen3v(ii,ie)
	  kn(ii)=kk
	  b(ii)=ev(ii+3,ie)
	  c(ii)=ev(ii+6,ie)

c	  zz = zeov(ii,ie)		!ZEONV
c	  zz = zov(kk)
	  zz = zov(kk) - zeqv(kk)	!tide
	  zz = zeov(ii,ie) - zeqv(kk)	!tide

	  zmm = zmm + zeov(ii,ie)		!ZEONV

	  bz=bz+zz*b(ii)
	  cz=cz+zz*c(ii)
	  bpres=bpres+ppv(kk)*b(ii)
	  cpres=cpres+ppv(kk)*c(ii)
          zm=zm+zz
	  taux=taux+tauxnv(kk)
	  tauy=tauy+tauynv(kk)
          rdist = rdist + rdistv(kk)
	end do

	zm=zm*drittl
	zmm=zmm*drittl
	taux=taux*drittl
	tauy=tauy*drittl
        rdist = rdist * drittl

c	gamma=af*dt*fcorv(ie)*rdist     !ggu advindex
c	gammat=fcorv(ie)*rdist

	gammat=fcorv(ie)*rdist 
        gamma=af*dt*gammat

c---------------------------------------------------------------------
c next part can be deleted
c---------------------------------------------------------------------
	hlh = hdeov(ilevel,ie)

        ulm = utlov(ilevel,ie)
        vlm = vtlov(ilevel,ie)

	hzg = hlh
        if(hzg.lt.hzoff) hzg=hzoff

        if(ireib.eq.0) then     !$$ireib0
                rt=0.
        else if(ireib.eq.1) then
		rt = rr
c		rt=rr/(hlh*hlh)	!for non slip condition...
        else if(ireib.eq.2) then
		if(bdebug) write(6,*) czv(ie),hzg,ulm,vlm
                gcz=grav/((czv(ie)**2)*(hzg**drittl))   !??
                rt=gcz*sqrt(ulm*ulm+vlm*vlm)/(hzg*hzg)  !??
        else if(ireib.eq.3) then
                gcz=grav/(czv(ie)**2)                   !??
                rt=gcz*sqrt(ulm*ulm+vlm*vlm)/(hzg*hzg)  !??
        else if(ireib.eq.5) then
                rt=rr*sqrt(ulm*ulm+vlm*vlm)/(hzg*hzg)        !??
        else if(ireib.eq.6) then
		! rr is z0
		rraux = cdf(hzg,rr)
                rt=rraux*sqrt(ulm*ulm+vlm*vlm)/(hzg*hzg)        !??
        else if(ireib.eq.7) then
		if( czv(ie) .ge. 1. ) then
                  gcz=grav/((czv(ie)**2)*(hzg**drittl))   !??
		else
		  gcz = czv(ie)
		end if
                rt=gcz*sqrt(ulm*ulm+vlm*vlm)/(hzg*hzg)  !??
	else if(ireib.eq.8) then
                ss = 0.
                do ii=1,3
                  kk = nen3v(ii,ie)
                  ss = ss + z0bk(kk)
                end do
                ss = ss / 3.		!this z0 from sedtrans
                rraux = cdf(hzg,ss)
                rt=rraux*sqrt(ulm*ulm+vlm*vlm)/(hzg*hzg)
        else
                write(6,*) 'unknown friction : ',ireib
                stop 'error stop : sp156'
        end if

	rraux = rt
	rt = rfricv(ie)		!this should be the right one
	if( reps .gt. 0. .and. abs(rraux-rt) .gt. reps ) then
		write(6,*) ie,ireib,rt,rraux
		stop 'error stop sp256v: rt - rraux'
	end if
c---------------------------------------------------------------------
c delete until here
c---------------------------------------------------------------------

	rt = rfricv(ie)		!bottom friction

c=========================================
c	if(rt.gt.rdt) rt=rdt	!$$rtmax
c=========================================

c reset in system (may be not the whole matrix every time)
c ...size of matrix : ngl*(2*mbw+1) with mbw=2
c
	do ii=1,ngl*5
	  rmat(ii)=0.
	end do

c compute layer thicknes and store to aux array

	hact(0) = 0.
	do l=1,ilevel
	  hact(l) = hdeov(l,ie)
	end do
	hact(ilevel+1) = 0.
	hact(nlvdim+1) = 0.

	if( bcolin ) then
	  hact(1) = hact(1) - zmm
	end if

c        call mimari(hact(0),ilevel+2,xmin,xmax,imin,imax,0.)
c        write(99,*) ie,xmin,xmax

c	write(6,*) ie,ilevel,(hact(l),l=1,ilevel)

	do l=0,ilevel+1
	  if( hact(l) .le. 0. ) then
	    rhact(l) = 0.
	  else
	    rhact(l) = 1. / hact(l)
	  end if
	end do

c        call mimari(rhact(0),ilevel+2,xmin,xmax,imin,imax,0.)
c        write(99,*) ie,xmin,xmax

c compute element averaged turbulent viscosity

	k1 = nen3v(1,ie)
	k2 = nen3v(2,ie)
	k3 = nen3v(3,ie)
	do l=0,ilevel
	    alev(l) = vismol + (visv(l,k1)+visv(l,k2)+visv(l,k3))/3.
	end do

c start of vertical loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do l=1,ilevel

c first set depth values
c
c hhi/hhip/hhim		thickness of i/i+1/i-1 layer
c uui/uuip/uuim		transport in x in i/i+1/i-1 layer
c vvi/vvip/vvim		transport in y in i/i+1/i-1 layer
c
c in case of a layer that does not exist (i-1 of first layer) give any
c ...value because the corrisponding a/b/c will be 0

	lp = min(l+1,ilevel)
	lm = max(l-1,1)

	hhi = hact(l)
	hhip = hact(l+1)
	hhim = hact(l-1)

	rhp = 0.
	rhm = 0.
	if( hhi .gt. 0. ) then		!may be redundant
	  if( hhip .gt. 0. ) then	!lower interface
	    rhp = 2.0 * alev(l) / ( hhi + hhip )
	  end if
	  if( hhim .gt. 0. ) then	!upper interface
	    rhm = 2.0 * alev(l-1) / ( hhi + hhim )
	  end if
	end if
        
	uui = utlov(l,ie)
	uuip = utlov(lp,ie)
	uuim = utlov(lm,ie)

	vvi = vtlov(l,ie)
	vvip = vtlov(lp,ie)
	vvim = vtlov(lm,ie)
        
	bfirst = l .eq. 1
	blast  = l .eq. ilevel
	
c	aus = afact * alev(l)
c	aux = dt * at * aus

	aus = 1.
	aux = dt * at

	aa  = aux * rhact(l) * ( rhm + rhp )
	aat = aus * rhact(l) * ( rhm + rhp )
	bb  = aux * rhact(l+1) * rhp
	bbt = aus * rhact(l+1) * rhp
	cc  = aux * rhact(l-1) * rhm
	cct = aus * rhact(l-1) * rhm

	ppx = 0.
	ppy = 0.
	if( bfirst ) then
	  ppx = ppx - taux
	  ppy = ppy - tauy
	end if
	if( blast ) then
	  aa  = aa + dt * rt
	  aat = aat + rt
	end if

c compute baroclinic contribution
        
        aux = rdist * hhi * rrho0!deb100407
        
        xexpl = rdist * fxv(l,ie)
        yexpl = rdist * fyv(l,ie)

	ppx = ppx + aat*uui - bbt*uuip - cct*uuim - gammat*vvi 
     +			+ grav*hhi*bz + (hhi/rowass)*bpres + xexpl 
     +  		- radx(l,ie)
	ppy = ppy + aat*vvi - bbt*vvip - cct*vvim + gammat*uui 
     +			+ grav*hhi*cz + (hhi/rowass)*cpres + yexpl 
     +  		- rady(l,ie)

	jv=l+l
	ju=jv-1

	rmat(locssp(ju,ju,ngl,mbb)) = 1. + aa
	rmat(locssp(jv,jv,ngl,mbb)) = 1. + aa
	rmat(locssp(jv,ju,ngl,mbb))=gamma
	rmat(locssp(ju,jv,ngl,mbb))=-gamma
c
	if(.not.blast) then
		rmat(locssp(ju,ju+2,ngl,mbb))=-bb
		rmat(locssp(jv,jv+2,ngl,mbb))=-bb
        end if
	if(.not.bfirst) then
		rmat(locssp(ju,ju-2,ngl,mbb))=-cc
		rmat(locssp(jv,jv-2,ngl,mbb))=-cc
        end if

	rvec(ju) = ppx
	rvec(jv) = ppy

	rvec(ngl+ju) = hhi		!ASYM_OPSPLT
	rvec(ngl+jv) = 0.d+0
	rvec(2*ngl+ju) = 0.d+0
	rvec(2*ngl+jv) = hhi

	bdebug = ie.eq.20.or.ie.eq.100.or.ie.eq.250
	bdebug=.false.
	if(bdebug) then
          write(6,*) 'hhh: ',hact(1),grav,dt
          write(6,*) 'hhh: ',az,am
	write(6,*) 'zzz: ',l,ppx,ppy
	write(6,*) 'zzz: ',aa,gamma,bb,cc
	end if

	end do
c
c end of assembling loop
c
c solution of vertical system
c
c        call gelb(rvec,rmat,ngl,1,mbb,mbb,epseps,ier)
        !call dgelb(rvec,rmat,ngl,1,mbb,mbb,epseps,ier)
        call dgelb(rvec,rmat,ngl,3,mbb,mbb,epseps,ier)		!ASYM_OPSPLT
c
	if(ier.ne.0) goto 99

        afix=1-iuvfix(ie)       !chao deb
c
	do l=1,ilevel
	if(bdebug) write(6,*) 'zzz: ',l,rvec(2*l-1), rvec(2*l)
	  utlnv(l,ie) = utlov(l,ie) - dt * rvec(2*l-1)*afix     !chao deb
	  vtlnv(l,ie) = vtlov(l,ie) - dt * rvec(2*l)*afix       !chao deb
	end do
	bdebug = .false.

	do l=1,ngl						!ASYM_OPSPLT
	  ddxv(l,ie) = rvec(ngl+l)
	  ddyv(l,ie) = rvec(2*ngl+l)
	end do
c
  500	continue

	if( barea0 .and. baroc .and. niter .le. 5 ) then  !$$BAROC_AREA0
	  write(6,*) 'sp256v: BAROC_AREA0 active ',nbcount
	end if

	return
   99	continue
	write(6,*) 'Error in inverting matrix (vertical system)'
	write(6,*) 'it, ier : ',it,ier
	write(6,*) 'ie,ilevel,ngl,mbb: ',ie,ilevel,ngl,mbb
	write(6,*) 'rvec: ',(rvec(l),l=1,ngl)
	write(6,*) 'matrix: '
	do ii=1,ngl
	  do l=-2,2
	    kk=locssp(ii,ii+l,ngl,mbb)
	    if(kk.eq.0) then
		auxaux(l)=0.
	    else
		auxaux(l)=rmat(kk)
	    end if
	  end do
	  write(6,*) auxaux
	end do
	write(6,*) 'hact: ',(hact(l),l=0,ilevel)
	write(6,*) 'alev: ',(alev(l),l=0,ilevel)
	call check_elem(ie)
	call check_nodes_in_elem(ie)
	stop 'error stop : sp256v'
	end

c******************************************************************

	subroutine sp256n

c post processing of time step
c
c written on 23.07.1997 by ggu  (from scratch)
c
	implicit none
c
c parameters
	include 'param.h'
	real drittl
	parameter (drittl=1./3.)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real eps1,eps2,pi,flag,high
	real grav,fcor,dcor,dirn,rowass,roluft
	integer itanf,itend,idt,nits,niter,it
	integer nlvdi,nlv
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /mkonst/ eps1,eps2,pi,flag,high
	common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
	common /femtim/ itanf,itend,idt,nits,niter,it
	common /level/ nlvdi,nlv
	integer nen3v(3,1)
	integer ilhv(1)
	real hev(1)
	real utlov(nlvdim,1),vtlov(nlvdim,1)
	real utlnv(nlvdim,1),vtlnv(nlvdim,1)
	real unv(1),vnv(1)
	common /nen3v/nen3v
	common /ilhv/ilhv
	common /hev/hev
	common /utlov/utlov, /vtlov/vtlov, /utlnv/utlnv, /vtlnv/vtlnv
	common /unv/unv, /vnv/vnv
	include 'ev.h'

        integer iuvfix(1)       !chao deb
        common /iuvfix/iuvfix   !chao deb
        integer afix            !chao deb

        real hdeov(nlvdim,1)
        common /hdeov/hdeov

	real zov(1),znv(1)
	common /zov/zov, /znv/znv
	real zeov(3,1)
	common /zeov/zeov

        double precision ddxv(2*nlvdim,neldim)  !ASYM
        double precision ddyv(2*nlvdim,neldim)
        common /ddxv/ddxv, /ddyv/ddyv
        save /ddxv/, /ddyv/

c local
	logical bcolin,bdebug
	integer ie,ii,l,kk
	integer ilevel
	integer ju,jv
	real az,am,dt,delta,azpar,ampar
	real bz,cz,um,vm,dz,zm
	real hact(nlvdim)
	real hact_new(nlvdim)
	real du,dv
c function
	integer iround
	real getpar
c save
	integer ncount
	save ncount
c data
	data ncount / 0 /

	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension in sp256n'

	ncount=ncount+1

c constants

	bcolin=iround(getpar('iclin')).ne.0	! linearized conti

	call getazam(azpar,ampar)
	az=azpar
	am=ampar
	call get_timestep(dt)

	delta = dt * grav * am 

	do ie=1,nel

	bdebug = ie.eq.20.or.ie.eq.100.or.ie.eq.250
	bdebug = bdebug .and. niter .eq. nits

	ilevel=ilhv(ie)

        afix=1-iuvfix(ie)       !chao deb

c compute barotropic pressure term

	bz=0.
	cz=0.
	zm = 0.
	do ii=1,3
	  kk=nen3v(ii,ie)
	  dz = znv(kk) - zeov(ii,ie)
	  zm = zm + zeov(ii,ie)		!ZEONV
	  bz = bz + dz * ev(ii+3,ie)
	  cz = cz + dz * ev(ii+6,ie)
c	  bz = bz + ev(ii+3,ie) * ( aux1 * zov(kk) + am * znv(kk) ) !FIXME
c	  cz = cz + ev(ii+6,ie) * ( aux1 * zov(kk) + am * znv(kk) )
	end do

c compute layer thicknes and store to aux array

	do l=1,ilevel
	  hact(l) = hdeov(l,ie)
	end do

	if( bcolin ) then
	  hact(1) = hact(1) - zm/3.
	end if

c new transports from u/v hat variable

	do l=1,ilevel

	  jv=l+l
	  ju=jv-1

	  !du = delta * hact(l) * bz
	  !dv = delta * hact(l) * cz

	  du = delta * ( ddxv(ju,ie)*bz + ddyv(ju,ie)*cz )	!ASYM_OPSPLT_CH
	  dv = delta *( ddxv(jv,ie)*bz + ddyv(jv,ie)*cz )	!ASYM_OPSPLT_CH

	  utlnv(l,ie) = utlnv(l,ie) - du*afix   !chao deb
	  vtlnv(l,ie) = vtlnv(l,ie) - dv*afix   !chao deb

	  bdebug=.false.
	  if( bdebug ) then
	    write(6,*) 'transp: ',ie,l,utlnv(l,ie),vtlnv(l,ie)
	  end if
          bdebug = ie.eq.20.or.ie.eq.100.or.ie.eq.250
	  bdebug = .false.
	  if(bdebug) then
	    write(6,*) 'nnn: ie,dt,h ',ie,dt,hact(l)
c	    write(6,*) 'nnn: am,g ',am,grav
	    write(6,*) 'nnn: bz,cz ',bz,cz
c	    write(6,*) 'nnn: b ',(ev(ii+3,ie),ii=1,3)
c	    write(6,*) 'nnn: c ',(ev(ii+6,ie),ii=1,3)
            write(6,*) 'nnn: u ',utlov(l,ie),utlnv(l,ie)
            write(6,*) 'nnn: v ',vtlov(l,ie),vtlnv(l,ie)
            write(6,*) 'nnn: zo ',(zov(nen3v(ii,ie)),ii=1,3)
            write(6,*) 'nnn: zn ',(znv(nen3v(ii,ie)),ii=1,3)
	  end if
	  bdebug=.false.
	end do

c barotropic transports

	um = 0.
	vm = 0.
	do l=1,ilevel
	  um = um + utlnv(l,ie)
	  vm = vm + vtlnv(l,ie)
	end do
	unv(ie) = um
	vnv(ie) = vm

	end do

	end

c******************************************************************

	subroutine sp256w(vf,va)

c computes vertical velocities
c
c velocities are computed on S/T points (top and bottom of layer)
c bottom velocity of the whole column is assumed to be 0
c -> maybe change this
c
c computes volume differences and from these computes vertical
c velocities at every time step so that the computed velocities
c satisfy the continuity equation for every single time step
c
c wlnv is computed horizontally at a node and vertically
c it is at the center of the layer -> there are nlv velocities
c computed
c
c b,c are 1/m, (phi is dimensionless)
c aj is m**2
c utlnv... is m**2/s
c dvol is in m**3/s
c vv is m**2 (area)
c
c wlnv is first used to accumulate volume difference -> dvol
c at the end it receives the vertical velocity
c
c wlnv (dvol)   aux array for volume difference
c vv            aux array for area
c
c written on 27.08.91 by ggu  (from scratch)
c 14.08.1998	ggu	w = 0 at open boundary nodes
c 20.08.1998	ggu	some documentation

	implicit none

c parameters
	include 'param.h'
c arguments
	!real vv(0:nlvdim,1)	!$$VERVEL
	real vf(nlvdim,1)
	real va(nlvdim,1)
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nlvdi,nlv
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /level/ nlvdi,nlv
	integer nen3v(3,1)
	integer ilhv(1)
	integer ilhkv(1)
	real utlov(nlvdim,1),vtlov(nlvdim,1)
	real utlnv(nlvdim,1),vtlnv(nlvdim,1),wlnv(0:nlvdim,1)
	common /nen3v/nen3v
	common /utlov/utlov, /vtlov/vtlov
	common /utlnv/utlnv, /vtlnv/vtlnv, /wlnv/wlnv
	common /ilhv/ilhv
	common /ilhkv/ilhkv
	integer iwegv(1)
	common /iwegv/iwegv
        real mfluxv(nlvdim,1)
        common /mfluxv/mfluxv
	include 'ev.h'
c local
	logical debug
	integer k,ie,ii,kk,l,lmax
	integer ilevel
        integer ibc,ibtyp
	real aj,wbot,wdiv,ff,atop
	real b,c
	real am,az,azt,azpar,ampar
	real ffn,ffo
	real volo,voln,dt,dvdt
c statement functions
	logical isein
        isein(ie) = iwegv(ie).eq.0
	include 'testbndo.h'

	logical is_zeta_bound
	real volnode

	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension in sp256w'

c initialize

	call getazam(azpar,ampar)
	az=azpar
	am=ampar
	azt = 1. - az
	call get_timestep(dt)

	do k=1,nkn
	  do l=1,nlv
	    vf(l,k)=0.
	    va(l,k)=0.
	    wlnv(l,k) = 0.
	  end do
	end do

c compute difference of velocities for each layer
c
c f(ii) > 0 ==> flux into node ii
c aj * ff -> [m**3/s]     ( ff -> [m/s]   aj -> [m**2]    b,c -> [1/m] )

	do ie=1,nel
	 !if( isein(ie) ) then		!FIXME
	  aj=4.*ev(10,ie)		!area of triangle / 3
	  ilevel = ilhv(ie)
	  do l=1,ilevel
	    do ii=1,3
		kk=nen3v(ii,ie)
		b = ev(ii+3,ie)
		c = ev(ii+6,ie)
		ffn = utlnv(l,ie)*b + vtlnv(l,ie)*c
		ffo = utlov(l,ie)*b + vtlov(l,ie)*c
		ff = ffn * az + ffo * azt
		vf(l,kk) = vf(l,kk) + 3. * aj * ff
		va(l,kk) = va(l,kk) + aj
	    end do
	  end do
	 !end if
	end do

c from vel difference get absolute velocity (w_bottom = 0)
c	-> wlnv(nlv,k) is already in place !
c	-> wlnv(nlv,k) = 0 + wlnv(nlv,k)
c w of bottom of last layer must be 0 ! -> shift everything up
c wlnv(nlv,k) is always 0
c
c dividing wlnv [m**3/s] by area [vv] gives vertical velocity
c
c in vv(l,k) is the area of the upper interface: a(l) = a_i(l-1)
c =>  w(l-1) = flux(l-1) / a_i(l-1)  =>  w(l-1) = flux(l-1) / a(l)

	do k=1,nkn
	  lmax = ilhkv(k)
	  wlnv(lmax,k) = 0.
	  debug = k .eq. 0
	  do l=lmax,1,-1
            voln = volnode(l,k,+1)
            volo = volnode(l,k,-1)
	    dvdt = (voln-volo)/dt
	    wdiv = vf(l,k) + mfluxv(l,k)
	    wlnv(l-1,k) = wlnv(l,k) + wdiv - dvdt
	    if( debug ) write(6,*) k,l,wdiv,wlnv(l,k),wlnv(l-1,k)
	  end do
	  !write(68,*) k,wlnv(0,k)
	  wlnv(0,k) = 0.	! ensure no flux across surface - is very small
	end do

	do k=1,nkn
	  lmax = ilhkv(k)
	  debug = k .eq. 0
	  do l=2,lmax
	    atop = va(l,k)
	    if( atop .gt. 0. ) then
	      wlnv(l-1,k) = wlnv(l-1,k) / atop
	      if( debug ) write(6,*) k,l,atop,wlnv(l-1,k)
	    end if
	  end do
	end do

c set w to zero at open boundary nodes (new 14.08.1998)
c
c FIXME	-> only for ibtyp = 1,2 !!!!

	do k=1,nkn
            !if( is_external_boundary(k) ) then	!bug fix 10.03.2010
            if( is_zeta_bound(k) ) then
	      do l=0,nlv
		wlnv(l,k) = 0.
	      end do
            end if
	end do

	return
	end

c******************************************************************

        subroutine dryz(rmat,v1,v2,zv)
c
c estimation of levels in dry areas
c
c rmat          band matrix already decomposed
c v1            auxiliary vector, is used by routine
c               ...to assemble constant vector
c v2		pivot for rmat
c zv		vector with new water levels
c
        implicit none
c
c arguments
        real rmat(1),v1(1),v2(1),zv(1)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1),iwegv(1),inodv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nen3v/nen3v, /iwegv/iwegv, /inodv/inodv
	include 'ev.h'
c local
        integer ii,ie,k,ii1,ii2,kk1,kk2,ier
	integer idry
        real epseps,z,aj,hh999
        real b(3),c(3)
c functions
        logical iskbnd,iskout,iseout
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
        iskout(k) = inodv(k).eq.-2
        iseout(ie) = iwegv(ie).ne.0
c
	idry = 0
        do k=1,nkn
          v1(k)=0.
	end do
c
	do ie=1,nel
          if( iseout(ie) ) then
	    idry = idry + 1
            do ii=1,3
              b(ii)=ev(ii+3,ie)
              c(ii)=ev(ii+6,ie)
            end do
            aj=ev(10,ie)
            hh999=12.*aj
            do ii=1,3
              k=nen3v(ii,ie)
              if( iskbnd(k) ) then
                z=zv(k)
                ii1=mod(ii,3)+1
                ii2=mod(ii1,3)+1
                kk1=nen3v(ii1,ie)
                kk2=nen3v(ii2,ie)
                v1(kk1)=v1(kk1)-hh999*(b(ii)*b(ii1)+
     +                      c(ii)*c(ii1))*z
                v1(kk2)=v1(kk2)-hh999*(b(ii)*b(ii2)+
     +                      c(ii)*c(ii2))*z
                v1(k)=0.
              end if
            end do
          end if
	end do
c
	if( idry .eq. 0 ) return	!no dry areas

	!stop 'error stop dryz: drying routine not yet working '	!ASYM

        !epseps=1.e-6
        !call mchb(v1,rmat,nkn,1,mbw,-1,epseps,ier)
        !if(ier.ne.0) goto 99

        !write(6,*) '....... new solution for dry areas .......'
        call lp_subst_system(nkn,mbw,rmat,v1,v2)	!gguexclude - comment
c
        do k=1,nkn
          if( iskout(k) ) zv(k)=v1(k)			!gguexclude - comment
	end do
c
	return
   99   continue
        write(6,*) 'ier from sp158s : ',ier
        stop 'error stop sp158s'
	end

c**************************************************
c
c transfered to subn35.f
c
c        function cdf(h,z0)
c
cc computes cd from h and z0
c
c        implicit none
c
c        real cdf
c        real h,z0
c
c        real kappa,cds
c
c        kappa = 0.4
c
c        cds = kappa / log( (z0+0.5*h) / z0 )
c
c        cdf = cds*cds
c
c        end
c
c**************************************************

        subroutine check_invers(ngl,mbb,rrmat,x,r)

        implicit none

        integer ngl,mbb
        double precision rrmat(1)
        double precision x(1)
        double precision r(1)

        integer i,j,ip
        double precision diff,acu,val
        integer locssp

        diff = 0.

        write(6,*) 'check_invers: ',ngl,mbb

        do i=1,ngl
          acu = 0.
          do j=1,ngl
            ip = locssp(i,j,ngl,mbb)
            if( ip .gt. 0 ) then
              val = rrmat(ip)
              acu = acu + val * x(j)
              write(6,*) i,j,ip,val
            end if
          end do
          diff = diff + abs(r(i)-acu)
        end do

        do i=1,ngl
          write(6,*) 'x,r: ',x(i),r(i)
        end do

        diff = diff / ngl
        write(6,*) 'weighted difference for matrix: ',diff

        end

c*******************************************************************

        subroutine interpolare_scalar

        call system_init

        call assamble_interpolare_scalar

        call system_solve_z

        call system_adjust_z

        end

c*******************************************************************

        subroutine assamble_interpolare_scalar

        implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real eps1,eps2,pi,flag,high
	common /mkonst/ eps1,eps2,pi,flag,high
        integer nen3v(3,1)
        common /nen3v/nen3v
	real rzv(1)
	common /rzv/rzv
	include 'ev.h'

        integer ie,i,kk
        integer n,m,j1,j2
        integer kn(3)
        real aj,hh999,rw
        real hia(3,3)
        real hik(3)
        real b(3),c(3)

        do ie=1,nel

	  do i=1,3
		kk=nen3v(i,ie)
		kn(i)=kk
		b(i)=ev(i+3,ie)
		c(i)=ev(i+6,ie)
	  end do

          aj = ev(10,ie)
          hh999=aj*12.

          do n=1,3
            do m=1,3
              hia(n,m)=hh999*(b(n)*b(m)+c(n)*c(m))
            end do
            hik(n)=0.
          end do

	  do i=1,3
	    if(rzv(kn(i)).ne.flag) then
		rw=rzv(kn(i))			!FIXME !ZNEW (this for znew)
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

c in hia(i,j),hik(i),i,j=1,3 is system

	  call system_assemble(nkn,mbw,kn,hia,hik)

	end do

        end

c*******************************************************************

	subroutine set_yaron

c momentum input for yaron

	implicit none

	include 'param.h'
	include 'ev.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhv(1)
        common /ilhv/ilhv
        integer nen3v(3,1)
        common /nen3v/nen3v
        real hdeov(nlvdim,1)
        common /hdeov/hdeov
        real fxv(nlvdim,1)
        real fyv(nlvdim,1)
        common /fxv/fxv
        common /fyv/fyv

	integer kin,lin,ie,ii,k,lmax,nelem
	real rnx,rny,rfact,q,area,h,fact

	kin = 3935
	kin = 0
	kin = 2088
	lin = 8
	nelem = 6
	nelem = 4
	rnx = -1
	rny = 0.
	rfact = 1.1
	q = 10.

	if( kin .le. 0 ) return

        do ie=1,nel
          lmax = ilhv(ie)
          area = 12. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( k .eq. kin .and. lmax .le. lin ) then
	      h = hdeov(lin,ie)
	      fact = rfact * q*q / (h*area*sqrt(area)*nelem)
	write(17,*) 'yaron: ',ie,k,fact,fxv(lin,ie)
	      fxv(lin,ie) = fxv(lin,ie) - fact * rnx
	      fyv(lin,ie) = fyv(lin,ie) - fact * rny
	    end if
	  end do
	end do

	end

c*******************************************************************

