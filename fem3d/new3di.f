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
c nknddi	dimension for total number of nodes
c nelddi	dimension for total number of elements
c nrbddi	dimension for total number of boundary condition nodes
c nbcddi	dimension for total number of open boundaries
c mbwddi	dimension for bandwidth
c ngrddi	dimension for grade of nodes (number of elements attached
c		...to one node)
c narddi	dimension for total number of area codes
c nexddi	dimension for total number of extra nodes for output
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
c 04.11.2011    ggu	deleted computation of firction term (in subn35.f)
c 29.03.2012    ggu	cleaned up, sp256v prepared for OpenMP
c 10.05.2013    dbf&ggu new routines for non-hydro
c 29.10.2013    ggu	nudging implemented
c 29.11.2013    ggu	zeta correction
c 25.03.2014    ggu     new offline
c 10.04.2014    ggu     cleaning up of a lot of stuff
c 06.05.2015    ggu     cleaning up of sp256f
c 20.05.2015    ggu&erp sp256v parallelized
c
c******************************************************************

	subroutine sp259f

c administrates one hydrodynamic time step for system to solve
c
c written on 27.07.88 by ggu   (from sp159f)

	use mod_depth
	use mod_bound_dynamic
	use mod_aux_array
	use mod_area
	use mod_hydro_baro
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'
	include 'femtime.h'



	logical boff,bdebout
	logical bzcorr
	integer i,l,k,ie,ier,ii
	integer nrand
	integer iw,iwa
	integer nmat
	integer kspecial
	integer iwhat
	real res

	integer iround
	real getpar,resi
	integer inohyd
	logical bnohyd

        real epseps
        parameter (epseps = 1.e-6)

	kspecial = 0
	bdebout = .false.

c-----------------------------------------------------------------
c set parameter for hydro or non hydro 
c-----------------------------------------------------------------

        inohyd = nint(getpar('inohyd'))
	bnohyd = inohyd .eq. 1

c-----------------------------------------------------------------
c offline
c-----------------------------------------------------------------

	call is_offline(1,boff)
	if( boff ) return

c-----------------------------------------------------------------
c dry areas
c-----------------------------------------------------------------

	iw=0
	call sp136(iw)

c-----------------------------------------------------------------
c copy variables to old time level
c-----------------------------------------------------------------

	call copy_uvz		!copies uvz to old time level
	call nonhydro_copy	!copies non hydrostatic pressure terms
	call copy_depth

	call diff_h_set		!horizontal viscosity

c-----------------------------------------------------------------
c solve for hydrodynamic variables
c-----------------------------------------------------------------

	iw = 1
	do while( iw .gt. 0 )		!loop over changing domain

	  if( bdebout ) call debug_output(it+2)
	  call sp256v			!compute intermediate transports
	  if( bdebout ) call debug_output(it+3)

	  call setnod			!set info on dry nodes
	  call set_link_info
	  call adjust_mass_flux		!cope with dry nodes

	  call system_init		!initializes matrix

	  if( bdebout ) call debug_output(it+14)
	  call sp256z(rqv)		!assemble system matrix for z
	  if( bdebout ) call debug_output(it+15)

	  call system_solve_z(nkn,znv)	!solves system matrix for z
	  if( bdebout ) call debug_output(it+16)

	  call system_adjust_z(nkn,znv)	!copies solution to new z
	  if( bdebout ) call debug_output(it+17)

	  call setweg(1,iw)		!controll intertidal flats

	end do	!do while( iw .gt. 0 )

	if( bdebout ) call debug_output(it+4)
	call sp256n			!final transports (also barotropic)
	if( bdebout ) call debug_output(it+5)

c-----------------------------------------------------------------
c end of soulution for hydrodynamic variables
c-----------------------------------------------------------------

        call setzev			!copy znv to zenv
        call setuvd			!set velocities in dry areas
	call baro2l 			!sets transports in dry areas

	call make_new_depth
	call check_volume		!checks for negative volume 
        call arper

	res=resi(zov,znv,nkn)

c-----------------------------------------------------------------
c vertical velocities and non-hydrostatic step
c-----------------------------------------------------------------

	if (bnohyd) then
	  call sp256wnh
	  call nonhydro_adjust
	end if

	call sp256w(v1v,saux1,saux2)	!compute vertical velocities

c-----------------------------------------------------------------
c correction for zeta
c-----------------------------------------------------------------

	bzcorr = .true.
	bzcorr = .false.
	if( bzcorr ) then
	  call correct_zeta(v1v)
          call setzev     !znv -> zenv
	  call make_new_depth
	  call sp256w(v1v,saux1,saux2)	!$$VERVEL
	end if

c-----------------------------------------------------------------
c some checks
c-----------------------------------------------------------------

	call vol_mass(1)		!computes and writes total volume
	if( bdebout ) call debug_output(it)
	call mass_conserve(saux1,saux2)	!check mass balance

c-----------------------------------------------------------------
c compute velocities on elements and nodes
c-----------------------------------------------------------------

	call ttov
	call make_prvel

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	return
   99	continue
	write(6,*) 'Error in inverting matrix for water level'
	write(6,*) 'it, ier : ',it,ier
	stop 'error stop : sp259f'
	end

c******************************************************************

	subroutine sp256z(vqv)

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

	use mod_nudging
	use mod_internal
	use mod_geom_dynamic
	use mod_depth
	use mod_bound_dynamic
	use mod_hydro_baro
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

	include 'param.h'

	real vqv(1)

	real drittl
	parameter (drittl=1./3.)

	include 'mkonst.h'
	include 'pkonst.h'
	include 'femtime.h'


        integer afix            !chao deb



	logical bcolin
	logical bdebug
	integer kn(3)
	integer ie,i,j,j1,j2,n,m,kk,l,k
	integer ngl
	integer ilevel
	integer ju,jv
	real az,am,af,azpar,ampar
	real dt,aj,rw
	real zm
	real ut,vt,uhat,vhat
	real ht
	real h11,hh999
	real delta
	real hia(3,3),hik(3),amatr(3,3)
	real b(3),c(3),z(3)
	real andg,zndg(3)
	real acu
	real uold,vold
	real dbb,dbc,dcb,dcc,abn,acn

c	data amatr / 2.,1.,1.,1.,2.,1.,1.,1.,2. /	!original
	data amatr / 4.,0.,0.,0.,4.,0.,0.,0.,4. /	!lumped

        integer locsps,loclp,iround
	real getpar
	!logical iskbnd,iskout,iseout
	logical iskbnd,iseout
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
        !iskout(k) = inodv(k).eq.-2
        iseout(ie) = iwegv(ie).ne.0

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	bcolin=iround(getpar('iclin')).ne.0

	call getazam(azpar,ampar)
	az=azpar
	am=ampar
	af=getpar('afpar')
	call get_timestep(dt)

	ngl=nkn

c-------------------------------------------------------------
c loop over elements
c-------------------------------------------------------------

	do ie=1,nel

c	------------------------------------------------------
c	compute level gradient
c	------------------------------------------------------

	zm=0.
	do i=1,3
		kk=nen3v(i,ie)
		kn(i)=kk
		b(i)=ev(i+3,ie)
		c(i)=ev(i+6,ie)
		!z(i)=zov(kk)
		z(i)=zeov(i,ie)		!ZEONV
		zndg(i) = andgzv(kk)	!nudging
		zm=zm+z(i)
	end do

	zm=zm*drittl

	if(bcolin) then
		ht=hev(ie)
	else
		ht=hev(ie)+zm
	end if

	ilevel=ilhv(ie)
	aj=ev(10,ie)
        afix=1-iuvfix(ie)      !chao deb

        delta=dt*dt*az*am*grav*afix         !ASYM_OPSPLT        !chao deb

c	------------------------------------------------------
c	compute contribution from H^x and H^y
c	------------------------------------------------------

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

c	------------------------------------------------------
c	compute barotropic transport
c	------------------------------------------------------

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

c	------------------------------------------------------
c	set element matrix and RHS
c	------------------------------------------------------

	do n=1,3
	  do m=1,3
	    abn = b(n) * ( b(m) * dbb + c(m) * dbc )
	    acn = c(n) * ( b(m) * dcb + c(m) * dcc )
	    h11 = delta*( abn + acn )			!ASYM_OPSPLT_CH
	    hia(n,m) = aj * (amatr(n,m) + 12.*h11)
	  end do
	  acu = hia(n,1)*z(1) + hia(n,2)*z(2) + hia(n,3)*z(3)
	  andg = 4.*aj*dt*zndg(n)
	  !hia(n,n) = hia(n,n) + 4 * dt * aj / tau
	  hik(n) = acu + andg + 12.*aj*dt*( ut*b(n) + vt*c(n) )	!ZNEW
	end do

c	------------------------------------------------------
c	level boundary conditions
c	------------------------------------------------------

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
		!hia(i,i)=12.*aj
		hik(i)=rw*hia(i,i)
	  end if
	  !call handle_ship_boundary(it,i,k,hia,hik)
	end do

c	------------------------------------------------------
c	excluded areas
c	------------------------------------------------------

          if( iseout(ie) ) then	!ZEONV
            hh999=aj*12.
            do n=1,3
              do m=1,3
                hia(n,m)=hh999*(b(n)*b(m)+c(n)*c(m))
              end do
              hik(n)=0.
            end do

            do n=1,3
              if( iskbnd(kn(n)) ) then	!not internal and not out of system
                do m=1,3
                  hia(n,m)=0.
                  !hia(m,n)=0.		!gguexclude - comment
                end do
                hik(n)=0.
              end if
            end do
          end if

c	------------------------------------------------------
c	in hia(i,j),hik(i),i,j=1,3 is system
c	------------------------------------------------------

	  call system_assemble(nkn,mbw,kn,hia,hik)

	end do

c-------------------------------------------------------------
c end of loop over elements
c-------------------------------------------------------------

	call system_add_rhs(dt,nkn,vqv)

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c******************************************************************

	subroutine sp256v

	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'pkonst.h'
	include 'femtime.h'

	integer ie
	integer ith
	integer count0,dcount,chunk,nt
	integer ibaroc
	integer ilin,itlin
	integer num_threads,myid,el_do,rest_do,init_do,end_do
	logical bcolin,baroc
	real az,am,af,at,av,azpar,ampar
	real rlin,radv
	real vismol,rrho0
	real dt

	double precision tempo

	double precision openmp_get_wtime
	!integer openmp_get_num_threads,openmp_get_thread_num
	real getpar

c-------------------------------------------------------------
c initialize
c-------------------------------------------------------------

	ibaroc = nint(getpar('ibarcl'))		! baroclinic contributions
        vismol  = getpar('vismol')		! molecular viscosity
	bcolin = nint(getpar('iclin')).ne.0	! linearized conti
	itlin = nint(getpar('itlin'))		! advection scheme
	ilin = nint(getpar('ilin'))		! non-linear terms?
	rlin = getpar('rlin')			! non-linear strength?

	baroc = ibaroc .eq. 1 .or. ibaroc .eq. 2

	call getazam(azpar,ampar)
	az=azpar
	am=ampar
	af=getpar('afpar')
	at=getpar('atpar')
	av=getpar('avpar')

	radv = 0.
	if( ilin .eq. 0 .and. itlin .eq. 0 ) then	!need non-lin terms
	  radv = rlin * av	!strength * implicit factor
	end if

	call get_timestep(dt)
    
	rrho0=1./rowass
	if( .not. baroc ) rrho0 = 0.

c-------------------------------------------------------------
c computation of explicit part (sets arrays fxv(l,ie),fyv(l,ie)
c-------------------------------------------------------------

	call bottom_friction	!set bottom friction
        call set_explicit       !new HYDRO deb
	!call set_yaron

c-------------------------------------------------------------
c parallel part
c-------------------------------------------------------------

ccc	call get_clock_count(count0)
ccc	nt = 2
ccc	call openmp_set_num_threads(nt)
ccc	chunk = 1 + nel/nt

c-------------------------------------------------------------
c loop over elements
c-------------------------------------------------------------

	tempo = openmp_get_wtime()

!$OMP PARALLEL PRIVATE(num_threads,myid,el_do,rest_do,init_do,end_do,ie)

        call openmp_get_num_threads(num_threads)
	if( num_threads == 0 ) num_threads = 1
        call openmp_get_thread_num(myid)
        el_do = nel/num_threads
        rest_do = MOD(nel,num_threads)
        init_do = el_do*myid+1
        end_do = init_do+el_do-1
        if( myid .eq. num_threads-1 ) end_do = end_do+rest_do

        !print *,'num_threads = ',num_threads
        !print *,"th = ",myid," el_do = ", el_do," rest_do = ", rest_do
        !print *,"init_do = ",init_do," end_do = ",end_do," nel  =", nel

	do ie=init_do,end_do

	  call sp256v_intern(ie,bcolin,baroc,az,am,af,at,radv
     +			,vismol,rrho0,dt)

	end do

!$OMP END PARALLEL      

	tempo = openmp_get_wtime() - tempo
	!write(66,*) it,tempo

c-------------------------------------------------------------
c end of loop over elements
c-------------------------------------------------------------

ccc	call get_clock_count_diff(count0,dcount)
ccc	write(6,*) 'count: ',dcount

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c******************************************************************

!DEC$ ATTRIBUTES INLINE :: sp256v_intern

	subroutine sp256v_intern(ie,bcolin,baroc,az,am,af,at,radv
     +			,vismol,rrho0,dt)

c assembles vertical system matrix
c
c semi-implicit scheme for 3d model
c
c written on 18.02.91 by ggu  (from scratch)
c
	use mod_tides
	use mod_meteo
	use mod_waves
	use mod_fluidmud
	use mod_internal
	use mod_depth
	use mod_layer_thickness
	use mod_aux_array
	use mod_roughness
	use mod_diff_visc_fric
	use mod_hydro_baro
	use mod_hydro_print
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

	integer ie
	logical bcolin,baroc
	real az,am,af,at
	real radv			!non-linear contribution
	real vismol,rrho0
	real dt

c parameters
	include 'param.h'
	real drittl
	parameter (drittl=1./3.)
c	real az,azt,am,amt,af,aft,at,att
c	parameter (az=0.50,azt=1.-az)
c	parameter (am=0.50,amt=1.-am)
c	parameter (af=0.50,aft=1.-af)
c	parameter (at=0.50,att=1.-at)
c common
	include 'mkonst.h'
	include 'pkonst.h'
	include 'femtime.h'

       integer afix             !chao deb








c local

	logical bbaroc,barea0                  !$$BAROC_AREA0

	logical bfirst,blast
	logical debug,bdebug
        logical bdebggu
	integer kn(3)
	integer kk,ii,l,ju,jv
	integer ngl,mbb
	integer ilevel,ier,ilevmin
	integer lp,lm
	integer k1,k2,k3,k
	real b(3),c(3)
	real zz
	real hlh,hlh_new
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
	real aust
	real fact                       !$$BCHAO - not used
	real uuadv,uvadv,vuadv,vvadv
        real rhp,rhm,aus
	real hzg,gcz
        real xmin,xmax
        integer imin,imax
        real rdist
        real xadv,yadv,fm,uc,vc,f,um,vm,up,vp
	real bpres,cpres
	real vis
        
	real rraux,cdf


	real ss

c-----------------------------------------
	real hact(0:nlvdi+1)
	real rhact(0:nlvdi+1)
	real alev(0:nlvdi)
c-----------------------------------------
	double precision rmat(10*nlvdi)
	double precision rvec(6*nlvdi)		!ASYM (3 systems to solve)
	double precision ppx,ppy
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

        real epseps
        parameter (epseps = 1.e-6)

c-------------------------------------------------------------
c initialization and baroclinic terms
c-------------------------------------------------------------

	bdebug=.false.
	debug=.false.
        barea0 = .false.     ! baroclinic only with ia = 0 (HACK - do not use)

        bbaroc = baroc
	if( barea0 ) then               !$$BAROC_AREA $$BAROC_AREA0
	  if( iarv(ie) .ne. 0 ) bbaroc = .false.
        end if

	rrho0=1./rowass
	if( .not. bbaroc ) rrho0 = 0.

c-------------------------------------------------------------
c dimensions of vertical system
c-------------------------------------------------------------

	ilevel=ilhv(ie)
	!ilevmin=ilmv(ie)
	ngl=2*ilevel
	mbb=2
	if(ngl.eq.2) mbb=1

c-------------------------------------------------------------
c compute barotropic terms (wind, atmospheric pressure, water level
c-------------------------------------------------------------

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

	  zz = zeov(ii,ie) - zeqv(kk)	!tide

          zm=zm+zz
	  zmm = zmm + zeov(ii,ie)		!ZEONV

	  bz=bz+zz*b(ii)
	  cz=cz+zz*c(ii)
	  bpres=bpres+ppv(kk)*b(ii)
	  cpres=cpres+ppv(kk)*c(ii)
	  taux=taux+tauxnv(kk)
	  tauy=tauy+tauynv(kk)
          rdist = rdist + rdistv(kk)
	end do

	zm=zm*drittl
	zmm=zmm*drittl
	taux=taux*drittl
	tauy=tauy*drittl
        rdist = rdist * drittl

c-------------------------------------------------------------
c coriolis parameter
c-------------------------------------------------------------

c	gamma=af*dt*fcorv(ie)*rdist     !ggu advindex
c	gammat=fcorv(ie)*rdist

	gammat=fcorv(ie)*rdist 
        gamma=af*dt*gammat

c-------------------------------------------------------------
c reset vertical system 
c
c may be not the whole matrix every time
c ...size of matrix : ngl*(2*mbw+1) with mbw=2
c-------------------------------------------------------------

	do ii=1,ngl*5
	  rmat(ii)=0.
	end do

c-------------------------------------------------------------
c compute layer thicknes and store in hact and rhact
c-------------------------------------------------------------

	hact(0) = 0.
	do l=1,ilevel
	  hact(l) = hdeov(l,ie)
	end do
	hact(ilevel+1) = 0.
	hact(nlvdi+1) = 0.

	if( bcolin ) then
	  hact(1) = hact(1) - zmm		!FIXME
	end if

	do l=0,ilevel+1
	  if( hact(l) .le. 0. ) then
	    rhact(l) = 0.
	  else
	    rhact(l) = 1. / hact(l)
	  end if
	end do

c-------------------------------------------------------------
c compute element averaged turbulent viscosity
c-------------------------------------------------------------

	k1 = nen3v(1,ie)
	k2 = nen3v(2,ie)
	k3 = nen3v(3,ie)
	do l=0,ilevel
	    vis = vismol
	    vis = vis + (visv(l,k1)+visv(l,k2)+visv(l,k3))/3.
	    vis = vis + (vts(l,k1)+vts(l,k2)+vts(l,k3))/3.
	    alev(l) = vis
	end do

c-------------------------------------------------------------
c start of vertical loop
c
c first set depth values
c
c hhi/hhip/hhim		thickness of i/i+1/i-1 layer
c uui/uuip/uuim		transport in x in i/i+1/i-1 layer
c vvi/vvip/vvim		transport in y in i/i+1/i-1 layer
c
c in case of a layer that does not exist (i-1 of first layer) give any
c ...value because the corrisponding a/b/c will be 0
c-------------------------------------------------------------

	do l=1,ilevel

	bfirst = l .eq. 1
	blast  = l .eq. ilevel
	
	lp = min(l+1,ilevel)
	lm = max(l-1,1)

	uui = utlov(l,ie)
	uuip = utlov(lp,ie)
	uuim = utlov(lm,ie)

	vvi = vtlov(l,ie)
	vvip = vtlov(lp,ie)
	vvim = vtlov(lm,ie)
        
	hhi = hact(l)
	hhip = hact(l+1)
	hhim = hact(l-1)

c	------------------------------------------------------
c	set up contributions of vertical viscosity
c	------------------------------------------------------

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

c	------------------------------------------------------
c	boundary conditions for stress on surface and bottom
c	------------------------------------------------------

	ppx = 0.
	ppy = 0.
	if( bfirst ) then
	  ppx = ppx - taux
	  ppy = ppy - tauy
	end if
	if( blast ) then
	  aa  = aa + dt * rfricv(ie)
	  aat = aat + rfricv(ie)
	end if

c	------------------------------------------------------
c	implicit advective contribution
c	------------------------------------------------------

	uuadv = 0.
	uvadv = 0.
	vuadv = 0.
	vvadv = 0.

	aux = dt * radv * rdist

	if( aux .gt. 0. ) then		!implict treatment of non-linear terms

	uc = uui/hhi
	vc = vvi/hhi

	do ii=1,3
          k = nen3v(ii,ie)
          up = saux2(l,k) / hhi
          vp = saux3(l,k) / hhi
          f = uui * b(ii) + vvi * c(ii)
          if( f .lt. 0. ) then    !flux out of node => into element
	    uuadv = uuadv + aux*b(ii)*( up - uc )
	    uvadv = uvadv + aux*c(ii)*( up - uc )
	    vuadv = vuadv + aux*b(ii)*( vp - vc )
	    vvadv = vvadv + aux*c(ii)*( vp - vc )
            !xadv = xadv + f * ( up - uc )
            !yadv = yadv + f * ( vp - vc )
          end if
	end do

	end if

c	------------------------------------------------------
c	explicit contribution (non-linear, baroclinic, diffusion)
c	------------------------------------------------------
        
        xexpl = rdist * fxv(l,ie)
        yexpl = rdist * fyv(l,ie)

c	------------------------------------------------------
c	ppx/ppy is contribution on the left side of equation
c	ppx corresponds to -F^x_l in the documentation
c	ppy corresponds to -F^y_l in the documentation
c	------------------------------------------------------

	ppx = ppx + aat*uui - bbt*uuip - cct*uuim - gammat*vvi 
     +			+ grav*hhi*bz + (hhi/rowass)*bpres + xexpl 
     +  		+ wavefx(l,ie)
	ppy = ppy + aat*vvi - bbt*vvip - cct*vvim + gammat*uui 
     +			+ grav*hhi*cz + (hhi/rowass)*cpres + yexpl 
     +  		+ wavefy(l,ie)

c	------------------------------------------------------
c	set up matrix A
c	------------------------------------------------------

	jv=l+l
	ju=jv-1

	rmat(locssp(ju,ju,ngl,mbb)) = 1. + aa + uuadv
	rmat(locssp(jv,jv,ngl,mbb)) = 1. + aa + vvadv
	rmat(locssp(jv,ju,ngl,mbb)) =  gamma  + vuadv
	rmat(locssp(ju,jv,ngl,mbb)) = -gamma  + uvadv

	if(.not.blast) then
		rmat(locssp(ju,ju+2,ngl,mbb)) = -bb
		rmat(locssp(jv,jv+2,ngl,mbb)) = -bb
        end if
	if(.not.bfirst) then
		rmat(locssp(ju,ju-2,ngl,mbb)) = -cc
		rmat(locssp(jv,jv-2,ngl,mbb)) = -cc
        end if

c	------------------------------------------------------
c	set up right hand side -F^x and -F^y 
c	------------------------------------------------------

	rvec(ju) = ppx
	rvec(jv) = ppy

c	------------------------------------------------------
c	set up H^x and H^y
c	------------------------------------------------------

	rvec(ngl+ju) = hhi		!ASYM_OPSPLT
	rvec(ngl+jv) = 0.d+0
	rvec(2*ngl+ju) = 0.d+0
	rvec(2*ngl+jv) = hhi

	end do

c-------------------------------------------------------------
c end of vertical loop
c-------------------------------------------------------------

c-------------------------------------------------------------
c solution of vertical system (we solve 3 systems in one call)
c-------------------------------------------------------------

        !call gelb(rvec,rmat,ngl,1,mbb,mbb,epseps,ier)
        !call dgelb(rvec,rmat,ngl,1,mbb,mbb,epseps,ier)
        call dgelb(rvec,rmat,ngl,3,mbb,mbb,epseps,ier)		!ASYM_OPSPLT

	if(ier.ne.0) goto 99

c-------------------------------------------------------------
c compute u^hat (negative sign because ppx/ppy was -F^x/-F^y)
c-------------------------------------------------------------

        afix=1-iuvfix(ie)       !chao deb

	do l=1,ilevel
	  utlnv(l,ie) = utlov(l,ie) - dt * rvec(2*l-1)*afix     !chao deb
	  vtlnv(l,ie) = vtlov(l,ie) - dt * rvec(2*l)*afix       !chao deb
	end do

c-------------------------------------------------------------
c save contribution A^{-1} H^x and A^{-1} H^y
c-------------------------------------------------------------

	do l=1,ngl						!ASYM_OPSPLT
	  ddxv(l,ie) = rvec(ngl+l)
	  ddyv(l,ie) = rvec(2*ngl+l)
	end do

c-------------------------------------------------------------
c special information
c-------------------------------------------------------------

	if( ie .eq. 1 .and. barea0 .and. 
     +			baroc .and. niter .le. 5 ) then  !$$BAROC_AREA0
	  write(6,*) 'sp256v: BAROC_AREA0 active '
	end if

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

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
	call check_set_unit(6)
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
	use mod_internal
	use mod_depth
	use mod_hydro_baro
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none
c
c parameters
	include 'param.h'
	real drittl
	parameter (drittl=1./3.)
c common
	include 'mkonst.h'
	include 'pkonst.h'
	include 'femtime.h'

        integer afix            !chao deb



c local
	logical bcolin,bdebug
	integer ie,ii,l,kk
	integer ilevel
	integer ju,jv
	real az,am,dt,beta,azpar,ampar
	real bz,cz,um,vm,dz,zm
	real du,dv
c function
	integer iround
	real getpar

c-------------------------------------------------------------
c initialize
c-------------------------------------------------------------

	bcolin=iround(getpar('iclin')).ne.0	! linearized conti
	bdebug = .false.

	call get_timestep(dt)
	call getazam(azpar,ampar)
	az=azpar
	am=ampar

	beta = dt * grav * am 

c-------------------------------------------------------------
c start loop on elements
c-------------------------------------------------------------

	do ie=1,nel

	ilevel=ilhv(ie)

        afix=1-iuvfix(ie)       !chao deb

c	------------------------------------------------------
c	compute barotropic pressure term
c	------------------------------------------------------

	bz=0.
	cz=0.
	zm = 0.
	do ii=1,3
	  kk=nen3v(ii,ie)
	  dz = znv(kk) - zeov(ii,ie)
	  zm = zm + zeov(ii,ie)		!ZEONV
	  bz = bz + dz * ev(ii+3,ie)
	  cz = cz + dz * ev(ii+6,ie)
	end do

c	------------------------------------------------------
c	new transports from u/v hat variable
c	------------------------------------------------------

	do l=1,ilevel

	  jv=l+l
	  ju=jv-1

	  du = beta * ( ddxv(ju,ie)*bz + ddyv(ju,ie)*cz )	!ASYM_OPSPLT_CH
	  dv = beta * ( ddxv(jv,ie)*bz + ddyv(jv,ie)*cz )	!ASYM_OPSPLT_CH

	  utlnv(l,ie) = utlnv(l,ie) - du*afix   !chao deb
	  vtlnv(l,ie) = vtlnv(l,ie) - dv*afix   !chao deb

	end do

c	------------------------------------------------------
c	barotropic transports
c	------------------------------------------------------

	um = 0.
	vm = 0.
	do l=1,ilevel
	  um = um + utlnv(l,ie)
	  vm = vm + vtlnv(l,ie)
	end do
	unv(ie) = um
	vnv(ie) = vm

	end do

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c******************************************************************

	subroutine sp256w(dzeta,vf,va)

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

	use mod_bound_geom
	use mod_geom_dynamic
	use mod_bound_dynamic
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

c parameters
	include 'param.h'
c arguments
	!real vv(0:nlvdi,1)	!$$VERVEL
	real dzeta(1)
	real vf(nlvdi,1)
	real va(nlvdi,1)
c common
c local
	logical debug
	integer k,ie,ii,kk,l,lmax
	integer ilevel
        integer ibc,ibtyp
	real aj,wbot,wdiv,ff,atop,abot,wfold
	real b,c
	real am,az,azt,azpar,ampar
	real ffn,ffo
	real volo,voln,dt,dvdt,q
	real dzmax,dz
c statement functions

	logical is_zeta_bound
	real volnode

	!logical isein
        !isein(ie) = iwegv(ie).eq.0

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
		!ff = ffn
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
c in va(l,k) is the area of the upper interface: a(l) = a_i(l-1)
c =>  w(l-1) = flux(l-1) / a_i(l-1)  =>  w(l-1) = flux(l-1) / a(l)

	dzmax = 0.

	do k=1,nkn
	  lmax = ilhkv(k)
	  wlnv(lmax,k) = 0.
	  debug = k .eq. 0
	  abot = 0.
	  do l=lmax,1,-1
	    atop = va(l,k)
            voln = volnode(l,k,+1)
            volo = volnode(l,k,-1)
	    dvdt = (voln-volo)/dt
	    q = mfluxv(l,k)
	    wdiv = vf(l,k) + q
	    !wfold = azt * (atop*wlov(l-1,k)-abot*wlov(l,k))
	    !wlnv(l-1,k) = wlnv(l,k) + (wdiv-dvdt+wfold)/az
	    wlnv(l-1,k) = wlnv(l,k) + wdiv - dvdt
	    abot = atop
	    if( debug ) write(6,*) k,l,wdiv,wlnv(l,k),wlnv(l-1,k)
	  end do
	  dz = dt * wlnv(0,k) / va(1,k)
	  dzmax = max(dzmax,abs(dz))
	  wlnv(0,k) = 0.	! ensure no flux across surface - is very small
	  dzeta(k) = dz
	end do

	!write(6,*) 'sp256w: dzmax = ',dzmax

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
	      dzeta(k) = 0.
            end if
	end do

	return
	end

c******************************************************************

	subroutine sp256wd(dzeta,vf,va)

c computes vertical velocities (double precision)
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

	use mod_bound_geom
	use mod_geom_dynamic
	use mod_bound_dynamic
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

c parameters
	include 'param.h'
c arguments
	!real vv(0:nlvdi,1)	!$$VERVEL
	real dzeta(1)
	real vf(nlvdi,1)
	real va(nlvdi,1)
c common
c local
	logical debug
	integer k,ie,ii,kk,l,lmax
	integer ilevel
        integer ibc,ibtyp
	real azpar,ampar
	real dt
	double precision aj,wbot,wdiv,ff,atop
	double precision b,c
	double precision am,az,azt
	double precision ffn,ffo
	double precision volo,voln,ddt,dvdt,q
	double precision dzmax,dz
	double precision vfd(nlvdi,nkn)
	double precision vad(nlvdi,nkn)
	double precision wlndv(0:nlvdi,nkn)
c statement functions

	logical is_zeta_bound
	real volnode


	!logical isein
        !isein(ie) = iwegv(ie).eq.0

c initialize

	call getazam(azpar,ampar)
	az=azpar
	am=ampar
	azt = 1. - az
	call get_timestep(dt)
	ddt = dt

	do k=1,nkn
	  do l=1,nlv
	    vfd(l,k)=0.
	    vad(l,k)=0.
	    wlndv(l,k) = 0.
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
		vfd(l,kk) = vfd(l,kk) + 3. * aj * ff
		vad(l,kk) = vad(l,kk) + aj
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
c in va(l,k) is the area of the upper interface: a(l) = a_i(l-1)
c =>  w(l-1) = flux(l-1) / a_i(l-1)  =>  w(l-1) = flux(l-1) / a(l)

	dzmax = 0.

	do k=1,nkn
	  lmax = ilhkv(k)
	  wlndv(lmax,k) = 0.
	  debug = k .eq. 0
	  do l=lmax,1,-1
            voln = volnode(l,k,+1)
            volo = volnode(l,k,-1)
	    dvdt = (voln-volo)/ddt
	    q = mfluxv(l,k)
	    wdiv = vfd(l,k) + q
	    wlndv(l-1,k) = wlndv(l,k) + wdiv - dvdt
	    if( debug ) write(6,*) k,l,wdiv,wlndv(l,k),wlndv(l-1,k)
	  end do
	  dz = ddt * wlndv(0,k) / vad(1,k)
	  dzmax = max(dzmax,abs(dz))
	  wlnv(0,k) = 0.	! ensure no flux across surface - is very small
	  dzeta(k) = dz
	end do

	write(6,*) 'dzmax: ',dzmax

	do k=1,nkn
	  lmax = ilhkv(k)
	  debug = k .eq. 0
	  do l=2,lmax
	    atop = vad(l,k)
	    if( atop .gt. 0. ) then
	      wlndv(l-1,k) = wlndv(l-1,k) / atop
	      if( debug ) write(6,*) k,l,atop,wlndv(l-1,k)
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
		wlndv(l,k) = 0.
	      end do
	      dzeta(k) = 0.
            end if
	end do

c copy double precision values to real values

	do k=1,nkn
	  lmax = ilhkv(k)
	  wlnv(0,k) = wlndv(0,k)
	  do l=1,lmax
	    wlnv(l,k) = wlndv(l,k)
	    vf(l,k) = vfd(l,k)
	    va(l,k) = vad(l,k)
	  end do
	end do

	return
	end

c******************************************************************

        subroutine dryz(rmat,v1,v2,zv)
c
c estimation of levels in dry areas
c
c this routine is not called anymore
c
c rmat          band matrix already decomposed
c v1            auxiliary vector, is used by routine
c               ...to assemble constant vector
c v2		pivot for rmat
c zv		vector with new water levels
c
	use mod_geom_dynamic
	use evgeom
	use basin

        implicit none
c
c arguments
        real rmat(1),v1(1),v2(1),zv(1)
c common
	include 'param.h'
c local
        integer ii,ie,k,ii1,ii2,kk1,kk2,ier
	integer idry
        real z,aj,hh999
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

        !call mchb(v1,rmat,nkn,1,mbw,-1,epseps,ier)
        !if(ier.ne.0) goto 99

        !write(6,*) '....... new solution for dry areas .......'
	stop 'error stop dryz: not yet ready...'
c                                             | should be integer
c					      v
        !call lp_subst_system(nkn,mbw,rmat,v1,v2)	!gguexclude - comment
c
        do k=1,nkn
          if( iskout(k) ) zv(k)=v1(k)			!gguexclude - comment
	end do
c
	return
   99   continue
        write(6,*) 'ier from mchb : ',ier
        stop 'error stop dryz'
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

	subroutine set_yaron

c momentum input for yaron

	use mod_internal
	use mod_layer_thickness
	use evgeom
	use levels
	use basin

	implicit none

	include 'param.h'


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

	subroutine correct_zeta(dzeta)

	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real dzeta(1)		!zeta correction

	include 'param.h'

	integer k

	do k=1,nkn
	  znv(k) = znv(k) + dzeta(k)
	end do

	end

c*******************************************************************

