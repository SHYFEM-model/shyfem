!
! $Id: new3di.f,v 1.35 2010-03-11 15:36:38 georg Exp $
!
! assembling linear system routine
!
! contents :
!
! subroutine hydro			administrates one time step
! subroutine hydro_zeta(vqv)		assemble matrix
! subroutine hydro_transports		computes transports (temporary)
! subroutine hydro_transports_final	computes transports (final)
! subroutine hydro_vertical(dzeta)	computes vertical velocities
!
! notes :
!
! look for "!ccc" to see important changes
!
! ASYM			passage to asymmetrix matrix
! ASYM_OPSPLT		new version without operator splitting (can leave)
! ASYM_OPSPLT_CH	new version -> change here
!
! nknddi	dimension for total number of nodes
! nelddi	dimension for total number of elements
! nrbddi	dimension for total number of boundary condition nodes
! nbcddi	dimension for total number of open boundaries
! mbwddi	dimension for bandwidth
! ngrddi	dimension for grade of nodes (number of elements attached
!		...to one node)
! narddi	dimension for total number of area codes
! nexddi	dimension for total number of extra nodes for output
!
! nkn,nel	total number of nodes/elements
! nrz,nrq	total number of nodes with water level/flux boundary conditions
! nrb		total number of nodes with boundary conditions
! nbc		total number of open boundaries
! ngr		maximum grade of nodes
! mbw		bandwidth of system matrix
! flag		flag value (to recognize boundary conditions)
! grav,dcor	gravitational accel./medium latitude of basin
! rowass,roluft	density of water/air
! itanf,itend	start/end time for simulation
! idt,nits	time step/total iterations to go
! niter,it	actual number of iterations/actual time
! nlvdi,nlv	dimension for levels, number of used levels
!
! nen3v(..,ie)	element index - node numbers (3) of element ie
!
! ipv(k)	external node number of node k
! ipev(ie)	external element number of element ie
!
! rqv(k),rzv(k)	flux/water level boundary conditions for node k
!		...if(rzv(k).eq.flag) --> no b.c. is specified for node k
! bnd(..,ib)	specification of boundary conditions for open boundary ib
! ev(..,ie)	geometric parameters for element ie
!		...1-3 = area, 4-6 = b, 7-9 = c, 10 = Aomega, 11-13 = angle
! uov(ie)	depth integrated transport (old time level)
! vov(ie)	...
! unv(ie)	depth integrated transport (new time level)
! vnv(ie)	...
! zov(k)	water level of old/new time level
! znv(k)	...
! ulov(l,ie)	velocity of old time level in x direction of layer l and elem ie
! vlov(l,ie)	velocity of old time level in y direction of layer l and elem ie
! wlov(l,i)	velocity of old time level in z direction of layer l and node i
! ulnv(l,ie)	velocity of new time level in x direction of layer l and elem ie
! vlnv(l,ie)	velocity of new time level in y direction of layer l and elem ie
! wlnv(l,i)	velocity of new time level in z direction of layer l and node i
! utlov(l,ie)	transport of old time level in x direction of layer l and el. ie
! vtlov(l,ie)	transport of old time level in y direction of layer l and el. ie
! utlnv(l,ie)	transport of new time level in x direction of layer l and el. ie
! vtlnv(l,ie)	transport of new time level in y direction of layer l and el. ie
! uprv(l,k)	velocity (averaged) in x direction of layer l and node k
! vprv(l,k)	velocity (averaged) in y direction of layer l and node k
! wprv(l,k)	velocity (averaged) in z direction of layer l and node k
! up0v(k)	total velocity (averaged) in x direction of node k
! vp0v(k)	total velocity (averaged) in y direction of node k
! tauxnv(k)	normalized stress in x-direction at node k
! tauynv(k)	normalized stress in y-direction at node k
! rhov(l,k)	density for level l and node k
! fcorv(ie)	coriolis parameter for elem ie
!
! visv(l,k)	vertical turbulent viscosity for layer l and node k (alv)
! difv(l,k)	vertical turbulent diffusivity for layer l and node k (slv)
!
! xgv(k),ygv(k)	coordinates of node k
!
! hldv(l)	thickness of layer l
! hlv(l)	absolute depth of bottom of layer l
! ilhv(ie)	number of levels for element ie
! ilhkv(k)	number of levels for node k
! hlhv(ie)	thickness of last layer for element ie
! hev(ie)	total depth at element ie (no water level)
! hkv(k)	total depth at node k (no water level)
! hm3v(3,ie)	depth of three nodes in element ie
!
! v1v,v2v...	auxiliary vectors
! rmat		band matrix for one vertical system (as above)
! rvec		constant vector for vertical system
!
! $$h1new	use new water level to compute velocities
!		(only for printing velocities important, but
!		algorithm has to be checked if consistent)
! $$rtmax	use maximal friction coefficient of rdt (=1./dt)
!
! revision log :
!
! revised 01.07.93 	$$UVBARO - u/vov introduced for	iteration on rad cond
! revised 03.11.93 	$$cmplerr - compiler warnings hydro
! revised 05.11.93 	$$fric - normal friction
! revised 05.11.93 	$$crador - crador call commented
! revised 05.11.93 	subroutine crador in file newcra.f
! revised 05.11.93 	$$VBARO-ERR - unv(ie)=vov(ie)
! revised 28.08.95 	$$BAROC_AREA - do baroc only for iarv(ie) = 0
! revised 30.08.95      $$AUST - austausch coefficient introduced
! revised 01.09.95      $$AWEIGH - area weighting of austausch coefficient
! revised 06.03.96 	$$BAROC_AREA0 - introduced baroc0
! revised 06.03.96 	$$VERT_AUST_ADJUST - adjustment of vert. aust. coef.
! revised 06.06.96 	$$BCHAO - modifications for vel. profile (temp.)
! revised 10.06.96 	$$UVPADV - modifications for advective term
! 14.08.1998	ggu	set w = 0 at open boundary nodes
! 20.08.1998	ggu	some documentation for sp256w
! 08.04.1999    ggu     equilibrium tide introduced (zeqv)
! 20.04.1999    ggu     converted to stress instead of wind (tauxnv...)
! 24.06.1999    ggu     call to rescur commented (use 2D call resid)
! 07.03.2000    ggu     eliminated VERT_AUST_ADJUST
! 07.03.2000    ggu     arrays for vertical turbulent diffusion coefficient
! 12.01.2001    ggu     solve for znv and not level difference (ZNEW)
! 09.11.2001    ggu     BCHAO commented out, no compiler directives
! 14.11.2001    ggu     all compiler directives eliminated
! 10.08.2003    ggu     deleted commented parts (radiation, etc..)
! 10.08.2003    ggu     some code in new subroutines: make_prvel, copy_uvz
! 13.08.2003    ggu     delete useless vars (nrand, epsit), no iteration
! 13.08.2003    ggu     call setnod, set_link_info
! 10.03.2004    ggu     RQVDT - value in rqv is now discharge [m**3/s]
! 10.01.2005    ggu     BAROC_DIST - scale baroclinic term with distance
! 25.01.2005    ggu     BUGADV - bugfix for advective terms
! 24.02.2005    ggu     new reynolds stresses -> green
! 15.03.2005    ggu     austv,aust eliminated, austau() in diff_h_set()
! 29.04.2005    ggu     semi-lagrangian advection (backadv)
! 29.04.2005    ggu     no baroclinic terms close to step (ilevmin)
! 04.11.2005    ggu     use itlin to decide about semi-lagrangian advection
! 23.03.2006    ggu     changed time step to double precision
! 31.05.2006    ggu     new friction type ireib=7
! 18.10.2006    ccf     radx,rady for radiation stress introduced
! 28.11.2006    ggu     in u/vadv is now difference and not absolute transport
! 02.04.2007    ggu     new algorithm (look for ASYM)
! 08.06.2007    ggu&deb restructured for new explicit terms
! 28.09.2007    ggu	deleted chao, semi-lagrange to newexp.f, no indov
! 24.06.2008    ggu	bpresv deleted
! 10.10.2008	ggu&mbj	prepared for pardiso -> modify system (gguexclude)
! 10.12.2008    ggu	use rfricv for bottom friction -> other can be deleted
! 18.12.2008    ggu	more debug info for error in vertical system
! 13.01.2009    ggu	Pardiso lib integrated
! 04.03.2009    ggu	matrix amat deleted from file -> only locally used
! 27.03.2009    ggu	call new routine adjust_mass_flux() for dry nodes
! 06.04.2009    ggu	deleted routine write_elem_vel_info()
! 07.05.2009    ggu	new routines for scalar interpolation (not finished)
! 10.03.2010    ggu	bug fix in sp256w() for ibtyp=2
! 11.03.2010    ggu	new routine check_volume() to check for negative vol
! 12.04.2010    ggu	ad hoc routine for Yaron
! 16.12.2010    ggu	in sp256w() account for changing volume (sigma)
! 19.02.2011    ccf	3D radiation stress
! 04.11.2011    ggu	deleted computation of firction term (in subn35.f)
! 29.03.2012    ggu	cleaned up, sp256v prepared for OpenMP
! 10.05.2013    dbf&ggu new routines for non-hydro
! 29.10.2013    ggu	nudging implemented
! 29.11.2013    ggu	zeta correction
! 25.03.2014    ggu     new offline
! 10.04.2014    ggu     cleaning up of a lot of stuff
! 06.05.2015    ggu     cleaning up of sp256f
! 20.05.2015    ggu&erp sp256v parallelized
! 17.09.2015    ggu	sp256w renamed to hydro_vertical
! 17.09.2015    ggu	sp259f renamed to hydro
! 18.09.2015    ggu	sp256 renamed to hydro_transports, file cleaned
! 20.11.2015    ggu&erp chunk size introduced, omp finalized
!
!******************************************************************
!---------------------------------------------------------------------------
        module hydrodynamic
!---------------------------------------------------------------------------
        contains
!---------------------------------------------------------------------------

        subroutine hydro

! administrates one hydrodynamic time step for system to solve
!
! written on 27.07.88 by ggu   (from sp159f)

        use depth
        use bnd_dynamic
        use bnd_admin
        use area
        use hydro_baro
        use hydro_print
        use hydro_vel
        use hydro_admin
        use levels, only : nlvdi,nlv,ilhv       
        !use basin, only : nkn,nel,ngr,mbw
        use basin
        use shympi
        use layer_thickness !ivb - dbg
        use diffusion  !ivb - dbg   
        use internal        !ivb - dbg        
        use ts              !ivb - dbg
        use meteo           !ivb - dbg
        use para
        use utility
        use elems_dealing
        use topological
        use check
        use wetdry
        use offline_data,       only:   is_offline
        use closing
        use nohydro
        use diffusion_admin
        use topological_admin
        use system
        use transforms
        use timing

        implicit none

        include 'femtime.h'

        logical boff,bdebout
        logical bzcorr
        integer i,l,k,ie,ier,ii
        integer nrand
        integer iw,iwa,iloop
        integer nmat
        integer kspecial
        integer iwhat
        double precision azpar,ampar
        double precision dzeta(nkn)

        double precision ahpar
        integer inohyd
        logical bnohyd

        double precision epseps
        parameter (epseps = 1.e-6)

        character*50 filename

        logical bout    !ivb
        integer iunit   !ivb

        double precision time1
        double precision,save :: time_transports,time_ssh
        double precision t_transp,t_ssh

        kspecial = 0
        bdebout = .false.

!-----------------------------------------------------------------
! set parameter for hydro or non hydro 
!-----------------------------------------------------------------

        inohyd = nint(getpar('inohyd'))
        bnohyd = inohyd .eq. 1

        azpar = getpar('azpar')
        ampar = getpar('ampar')
        if( azpar == 0. .and. ampar == 1. ) then
          call system_set_explicit
        else if( azpar == 1. .and. ampar == 0. ) then
          call system_set_explicit
	!else if( shympi_is_parallel() ) then
	!  if( shympi_is_master() ) then
	!    write(6,*) 'system is not solved explicitly'
	!    write(6,*) 'cannot solve semi-implicitly with MPI'
	!    write(6,*) 'az,am: ',azpar,ampar
	!  end if
	  !call shympi_stop('no semi-implicit solution')
        end if

!-----------------------------------------------------------------
! offline
!-----------------------------------------------------------------

	call is_offline(1,boff)
	!if( boff ) write(6,*) 'hydro reading from offline...'
        if( boff ) return

!-----------------------------------------------------------------
! dry areas
!-----------------------------------------------------------------

        iw=0
        call sp136(iw)

!-----------------------------------------------------------------
! copy variables to old time level
!-----------------------------------------------------------------

        call copy_uvz           !copies uvz to old time level
                                ! zov = znv, zeov = zenv
                                ! utlov = utlnv,vtlov = vtlnv,wlov  = wlnv                        

        call nonhydro_copy      !copies non hydrostatic pressure terms

        call copy_depth         !copies layer structure 
                                ! ivb - hdkov = hdknv
                                !       hdknv = 0
                                !       hdeov = hdenv

        call set_diffusivity    !horizontal viscosity/diffusivity (needs uvprv)

!-----------------------------------------------------------------
! solve for hydrodynamic variables
!-----------------------------------------------------------------


        iloop = 0

        do                              !loop over changing domain

          iloop = iloop + 1
 
          if(ln_timing) time1 = shympi_wtime() 
          call hydro_transports         !compute intermediate transports
          if(ln_timing )transp_time = transp_time + shympi_wtime() - time1


          call setnod                   !set info on dry nodes
          call set_link_info            !information on areas, islands, etc..
          call adjust_mass_flux         !cope with dry nodes

          call system_init              !initializes matrix
          call hydro_zeta(rqv)          !assemble system matrix for z

          call system_solve_z(nkn,znv)  !solves system matrix for z

          call system_adjust_z(nkn,znv) !copies solution to new z

          !call shympi_comment('exchanging znv')
          call shympi_exchange_2d_node(znv)

          call setweg(1,iw)             !controll intertidal flats
          iw = shympi_sum(iw)
          if( iw == 0 ) exit

          ahpar = getpar('ahpar')

          if(shympi_partition_on_elements() .and. ahpar .gt. 0) then
            call send_halo(utlnv,nlvdi,nel_local,'ut')
            call send_halo(vtlnv,nlvdi,nel_local,'vt')
          end if

        end do                  !! -- end loop over changing domain --

        call hydro_transports_final     !final transports (also barotropic)

!-----------------------------------------------------------------
! end of soulution for hydrodynamic variables
!-----------------------------------------------------------------

        call setzev                     !copy znv to zenv
        call setuvd                     !set velocities in dry areas
        call baro2l                     !sets transports in dry areas

        ahpar = getpar('ahpar')

        if(shympi_partition_on_elements() .and. ahpar .gt. 0) then
          call send_halo(utlnv,nlvdi,nel_local,'ut')
          call send_halo(vtlnv,nlvdi,nel_local,'vt')
        end if

        call make_new_depth
        call check_volume               !checks for negative volume 
        call arper

!-----------------------------------------------------------------
! vertical velocities and non-hydrostatic step
!-----------------------------------------------------------------

        if (bnohyd) then
          call sp256wnh
          call nonhydro_adjust
        end if

        call hydro_vertical(dzeta)              !compute vertical velocities
        
!-----------------------------------------------------------------
! correction for zeta
!-----------------------------------------------------------------

        bzcorr = .true.
        bzcorr = .false.
        if( bzcorr ) then
          znv = znv + dzeta
          call setzev     !znv -> zenv
          call make_new_depth
          call hydro_vertical(dzeta)            !$$VERVEL
        end if

!-----------------------------------------------------------------
! some checks
!-----------------------------------------------------------------

        call vol_mass(1)                !computes and writes total volume
        if( bdebout ) call debug_output(it)
        call mass_conserve              !check mass balance

!-----------------------------------------------------------------
! compute velocities on elements and nodes
!-----------------------------------------------------------------

        call ttov
        call make_prvel

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
!   99	continue
!	write(6,*) 'Error in inverting matrix for water level'
!	write(6,*) 'it, ier : ',it,ier
!	stop 'error stop : hydro'
	end

!******************************************************************

	subroutine hydro_zeta(vqv)

! assembles linear system matrix
!
! vqv		flux boundary condition vector
!
! semi-implicit scheme for 3d model
!
! written on 18.02.91 by ggu  (from scratch)
! changed on 04.06.91 by ggu  (c=(1) : friction term has been corrected)
! changed on 01.10.92 by ggu  (staggered FE - completely restructured)
! 12.01.2001    ggu     solve for znv and not level difference (ZNEW)

	use nudging
	use internal
	use geom_dynamic
	use depth
	use bnd_dynamic
	use hydro_baro
	use hydro_admin
	use evgeom
	use levels
	use basin
        use shympi
        use para
        use utility
        use system
        use time_util

	implicit none

	double precision vqv(nkn)

	double precision drittl
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
	double precision az,am,af,azpar,ampar
	double precision dt,aj,rw
	double precision zm
	double precision ut,vt,uhat,vhat
	double precision ht
	double precision h11,hh999
	double precision delta
	double precision hia(3,3),hik(3),amatr(3,3)
	double precision b(3),c(3),z(3)
	double precision andg,zndg(3)
	double precision acu
	double precision uold,vold
	double precision dbb,dbc,dcb,dcc,abn,acn

#ifdef DEBUGON
        integer e
        double precision temp_ddxv(nlvdi,nel_local)
        double precision temp_ddyv(nlvdi,nel_local)
#endif
!	data amatr / 2.,1.,1.,1.,2.,1.,1.,1.,2. /	!original
	data amatr / 4.,0.,0.,0.,4.,0.,0.,0.,4. /	!lumped

        integer locsps,loclp
	!logical iskbnd,iskout,iseout
	!logical iskbnd,iseout
        !iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
        !iskout(k) = inodv(k).eq.-2
        !iseout(ie) = iwegv(ie).ne.0

!-------------------------------------------------------------
! initialization
!-------------------------------------------------------------

        bcolin=iround(getpar('iclin')).ne.0

        call getazam(azpar,ampar)
        az=azpar
        am=ampar
        af=getpar('afpar')
        call get_timestep(dt)

        ngl=nkn

        hik = 0.
        hia = 0.

!-------------------------------------------------------------
! loop over elements
!-------------------------------------------------------------

#ifdef DEBUGON
        do ie=1,nel_local
          do l=1,nlvdi
            temp_ddxv(l,ie)=ddxv(l,ie)
            temp_ddyv(l,ie)=ddyv(l,ie)
          end do
        end do
        call shympi_exchange_halo_3d_elems(temp_ddxv)
        call shympi_exchange_halo_3d_elems(temp_ddyv)
        do ie=1,nel_local
          do l=1,nlvdi
            ddxv(l,ie)=temp_ddxv(l,ie)
            temp_ddxv(l,ie)=ddxv(l+nlvdi,ie)
            ddyv(l,ie)=temp_ddyv(l,ie)
            temp_ddyv(l,ie)=ddyv(l+nlvdi,ie)
          end do
        end do
        call shympi_exchange_halo_3d_elems(temp_ddxv)
        call shympi_exchange_halo_3d_elems(temp_ddyv)
        do ie=1,nel_local
          do l=1,nlvdi
            ddxv(nlvdi+l,ie)=temp_ddxv(l,ie)
            ddyv(nlvdi+l,ie)=temp_ddyv(l,ie)
          end do
        end do
        call shympi_exchange_halo_3d_elems3(zeov)
        call shympi_exchange_halo_2d_nodes(andgzv)
        call shympi_exchange_halo_2d_elems(hev)
        call shympi_exchange_halo_2d_elems(ilhv)
        call shympi_exchange_halo_2d_elems(iuvfix)
        call shympi_exchange_halo_3d_elems(utlov)
        call shympi_exchange_halo_3d_elems(vtlov)
        call shympi_exchange_halo_3d_elems(utlnv)
        call shympi_exchange_halo_3d_elems(vtlnv)
        call shympi_exchange_halo_2d_nodes(rzv)
        call shympi_exchange_halo_2d_nodes(inodv)

        do e=1,nel_local
          if(bmpi) then
            ie=domain%elems%mapID(e)
          else
            ie=e
          end if
#else
        do ie=1,nel
#endif

!	------------------------------------------------------
!	compute level gradient
!	------------------------------------------------------

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

!	------------------------------------------------------
!	compute contribution from H^x and H^y
!	------------------------------------------------------

        dbb = 0.d0
        dbc = 0.d0
        dcb = 0.d0
        dcc = 0.d0
        do l=1,ilevel			!ASYM_OPSPLT
          jv=l+l
          ju=jv-1
          dbb = dbb + ddxv(ju,ie)
          dbc = dbc + ddyv(ju,ie)
          dcb = dcb + ddxv(jv,ie)
          dcc = dcc + ddyv(jv,ie)
        end do

!	------------------------------------------------------
!	compute barotropic transport
!	------------------------------------------------------

        uold = 0.d0
        vold = 0.d0
        uhat = 0.d0
        vhat = 0.d0

        do l=1,ilevel
          uold = uold + utlov(l,ie)
          vold = vold + vtlov(l,ie)
          uhat = uhat + utlnv(l,ie)
          vhat = vhat + vtlnv(l,ie)
        end do

	ut = az * uhat + (1.-az) * uold
	vt = az * vhat + (1.-az) * vold

!	------------------------------------------------------
!	set element matrix and RHS
!	------------------------------------------------------

	do n=1,3
	  do m=1,3
	    abn = b(n) * ( b(m) * dbb + c(m) * dbc )
	    acn = c(n) * ( b(m) * dcb + c(m) * dcc )
	    h11 = delta*( abn + acn )                   !ASYM_OPSPLT_CH
	    hia(n,m) = aj * (amatr(n,m) + 12.*h11)
	  end do
	  acu = hia(n,1)*z(1) + hia(n,2)*z(2) + hia(n,3)*z(3)
	  andg = 4.*aj*dt*zndg(n)
	  !hia(n,n) = hia(n,n) + 4 * dt * aj / tau
	  hik(n) = acu + andg + 12.*aj*dt*( ut*b(n) + vt*c(n) ) !ZNEW
	end do

!	------------------------------------------------------
!	level boundary conditions
!	------------------------------------------------------

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

!	------------------------------------------------------
!	excluded areas
!	------------------------------------------------------

          !if( iseout(ie) ) then	!ZEONV
          if( iwegv(ie).ne.0 ) then	!ZEONV
            hh999=aj*12.
            do n=1,3
              do m=1,3
                hia(n,m)=hh999*(b(n)*b(m)+c(n)*c(m))
              end do
              hik(n)=0.
            end do

            do n=1,3
              !if( iskbnd(kn(n)) ) then	!not internal and not out of system
              if( inodv(kn(n)).ne.0 .and. inodv(kn(n)).ne.-2 ) then	!not internal and not out of system
                do m=1,3
                  hia(n,m)=0.
                  !hia(m,n)=0.		!gguexclude - comment
                end do
                hik(n)=0.
              end if
            end do
          end if

!	------------------------------------------------------
!	in hia(i,j),hik(i),i,j=1,3 is system
!	------------------------------------------------------

	  !call system_assemble(ie,nkn,mbw,kn,hia,hik)
	  call system_assemble(ie,kn,hia,hik)

	end do

!-------------------------------------------------------------
! end of loop over elements
!-------------------------------------------------------------

	call system_add_rhs(dt,nkn,vqv)

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!******************************************************************

	subroutine hydro_transports

        use basin, only : nkn,nel,ngr,mbw
!$	use omp_lib
        use diffusion  !ivb - dbg
        use levels              
        use para
        use chezy
        use explicit
        use time_util

	implicit none

	include 'pkonst.h'
	include 'femtime.h'

	integer ie
	integer ies,iend
	integer ith
	integer count0,dcount,chunk,nt
	integer ibaroc
	integer ilin,itlin
	integer num_threads,myid,el_do,rest_do,init_do,end_do
	integer nchunk,nthreads
	logical bcolin,baroc
	double precision az,am,af,at,av,azpar,ampar
	double precision rlin,radv
	double precision vismol,rrho0
	double precision dt

	double precision tempo
	double precision openmp_get_wtime
	!integer openmp_get_num_threads,openmp_get_thread_num

        character*50 filename
        integer l

!-------------------------------------------------------------
! initialize
!-------------------------------------------------------------

        ibaroc  = nint(getpar('ibarcl'))        ! baroclinic contributions
        vismol  = getpar('vismol')              ! molecular viscosity
        bcolin  = nint(getpar('iclin')).ne.0    ! linearized conti
        itlin   = nint(getpar('itlin'))         ! advection scheme
        ilin    = nint(getpar('ilin'))          ! non-linear terms?
        rlin    = getpar('rlin')                ! non-linear strength?

        baroc = ibaroc .eq. 1 .or. ibaroc .eq. 2

        call getazam(azpar,ampar)
        az=azpar                        ! weighting in continuity
        am=ampar                        ! weighting in momentum
        af=getpar('afpar')              ! weighting of coriolis term
        at=getpar('atpar')              ! weighting of vertical viscosity
        av=getpar('avpar')              ! weighting of advective terms

        radv = 0.d0
        if( ilin .eq. 0 .and. itlin .eq. 0 ) then	!need non-lin terms
          radv = rlin * av      !strength * implicit factor
        end if

        call get_timestep(dt)
    
        rrho0=1./rowass
        if( .not. baroc ) rrho0 = 0.

!-------------------------------------------------------------
! computation of explicit part (sets arrays fxv(l,ie),fyv(l,ie)
!-------------------------------------------------------------

        call bottom_friction    !set bottom friction

        call set_explicit       !new HYDRO deb
        !call set_yaron

!-------------------------------------------------------------
! parallel part
!-------------------------------------------------------------

!cc	call get_clock_count(count0)
!cc	nt = 2
!cc	call openmp_set_num_threads(nt)
!cc	chunk = 1 + nel/nt

        nthreads = 1
!$      nthreads = omp_get_num_threads()
        nchunk = 1
!$      nchunk = nel / ( nthreads * 10 )
        nchunk = max(nchunk,1)

!-------------------------------------------------------------
! loop over elements
!-------------------------------------------------------------

	!tempo = openmp_get_wtime()

!$OMP PARALLEL 
!$OMP SINGLE

        do ie=1,nel,nchunk

!$OMP TASK FIRSTPRIVATE(ie,bcolin,baroc,az,am,af,at,radv        &
!$OMP   & 	   ,vismol,rrho0,dt) PRIVATE(ies,iend)          &
!$OMP   &     SHARED(nel,nchunk)	 DEFAULT(NONE)
         
          iend = ie+nchunk-1
          if(iend .gt. nel) iend = nel

          do ies=ie,iend
            call sp256v_intern(ies,bcolin,baroc,az,am,af,at,radv,vismol,rrho0,dt)
          end do

!$OMP END TASK

        end do

!$OMP END SINGLE
!$OMP TASKWAIT	
!$OMP END PARALLEL      

	!tempo = openmp_get_wtime() - tempo
	!write(66,*) it,tempo

!-------------------------------------------------------------
! end of loop over elements
!-------------------------------------------------------------

!cc	call get_clock_count_diff(count0,dcount)
!cc	write(6,*) 'count: ',dcount

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!******************************************************************

        subroutine sp256v_intern(ie,bcolin,baroc,az,am,af,at,radv,vismol,rrho0,dt)

! assembles vertical system matrix
!
! semi-implicit scheme for 3d model
!
! written on 18.02.91 by ggu  (from scratch)
!
        use tidef
        use meteo
        use waves
        use fluidmud
        use internal
        use depth
        use layer_thickness
        use roughness
        use diffusion
        use hydro_baro
        use hydro_print
        use hydro_admin
        use evgeom
        use levels
        use basin
        use shympi
        use utility

        implicit none

        integer ie
        logical bcolin,baroc
        double precision az,am,af,at
        double precision radv                       !non-linear contribution
        double precision vismol,rrho0
        double precision dt

! parameters
        double precision drittl
        parameter (drittl=1.d0/3.d0)
! common
        include 'mkonst.h'
        include 'pkonst.h'
        include 'femtime.h'

! local

        logical bbaroc,barea0                  !$$BAROC_AREA0

        integer afix             !chao deb
        logical bfirst,blast
        logical debug,bdebug
        logical bdebggu
        integer kn(3)
        integer kk,ii,l,ju,jv
        integer ngl,mbb
        integer ilevel,ier
        integer lp,lm
        integer k1,k2,k3,k
        double precision b(3),c(3)
        double precision zz
        double precision hlh,hlh_new
!	double precision ht
!	double precision bz,cz,dbz,dcz,zm
        double precision bz,cz,zm,zmm
        double precision xbcl,ybcl
        double precision xexpl,yexpl
!	double precision beta
        double precision ulm,vlm,taux,tauy
        double precision gamma,gammat
        double precision hhi,hhim,hhip,uui,uuim,uuip,vvi,vvim,vvip
!	double precision bb,bbt,cc,cct,aa,aat,ppx,ppy,aux,aux1,aux2
        double precision bb,bbt,cc,cct,aa,aat,aux
        double precision aust
        double precision fact                       !$$BCHAO - not used
        double precision uuadv,uvadv,vuadv,vvadv
        double precision rhp,rhm,aus
        double precision hzg,gcz
        double precision xmin,xmax
        integer imin,imax
        double precision rdist
        double precision xadv,yadv,fm,uc,vc,f,um,vm,up,vp
        double precision bpres,cpres
        double precision vis
        logical, parameter :: debug_mpi = .false.
        
        double precision rraux,cdf
        double precision ss

!-----------------------------------------
        double precision hact(0:nlvdi+1)
        double precision rhact(0:nlvdi+1)
        double precision alev(0:nlvdi)
!-----------------------------------------
        double precision rmat(10*nlvdi)
        double precision rvec(6*nlvdi)          !ASYM (3 systems to solve)
        double precision ppx,ppy
!-----------------------------------------
! function
        integer locssp

        double precision epseps
        parameter (epseps = 1.e-6)

        character*50 filename
        logical dbg1,dbg2

        ier=0
        rvec  = 0.d0

!-------------------------------------------------------------
! initialization and baroclinic terms
!-------------------------------------------------------------

        bdebug=.false.
        debug=.false.
        barea0 = .false.     ! baroclinic only with ia = 0 (HACK - do not use)

        bbaroc = baroc
        if( barea0 ) then               !$$BAROC_AREA $$BAROC_AREA0
          if( iarv(ie) .ne. 0 ) bbaroc = .false.
        end if

!-------------------------------------------------------------
! dimensions of vertical system
!-------------------------------------------------------------

        ilevel = ilhv(ie)
        ngl = 2*ilevel
        mbb = 2
        if(ngl.eq.2) mbb=1

!-------------------------------------------------------------
! compute barotropic terms (wind, atmospheric pressure, water level
!-------------------------------------------------------------

        bz      = 0.d0
        cz      = 0.d0
        bpres   = 0.d0
        cpres   = 0.d0
        zm      = 0.d0
        zmm     = 0.d0
        taux    = 0.d0
        tauy    = 0.d0
        rdist   = 0.d0

        do ii=1,3
          kk = nen3v(ii,ie)
          kn(ii) = kk
          b(ii) = ev(ii+3,ie)
          c(ii) = ev(ii+6,ie)

          zz = zeov(ii,ie) - zeqv(kk)   !tide

          zm = zm+zz
          zmm = zmm + zeov(ii,ie)       !ZEONV

          bz    = bz+zz*b(ii)
          cz    = cz+zz*c(ii)
          bpres = bpres+ppv(kk)*b(ii)
          cpres = cpres+ppv(kk)*c(ii)
          taux  = taux+tauxnv(kk)
          tauy  = tauy+tauynv(kk)
          rdist = rdist + rdistv(kk)
        end do

        zm      = zm*drittl
        zmm     = zmm*drittl
        taux    = taux*drittl
        tauy    = tauy*drittl
        rdist   = rdist * drittl

!-------------------------------------------------------------
! coriolis parameter
!-------------------------------------------------------------

!	gamma=af*dt*fcorv(ie)*rdist     !ggu advindex
!	gammat=fcorv(ie)*rdist

        gammat=fcorv(ie)*rdist 
        gamma=af*dt*gammat

!-------------------------------------------------------------
! reset vertical system 
!
! may be not the whole matrix every time
! ...size of matrix : ngl*(2*mbw+1) with mbw=2
!-------------------------------------------------------------

        do ii=1,ngl*5
          rmat(ii)=0.d0
        end do

!-------------------------------------------------------------
! compute layer thicknes and store in hact and rhact
!-------------------------------------------------------------

        hact(0) = 0.d0
        do l=1,ilevel
          hact(l) = hdeov(l,ie)
        end do
        hact(ilevel+1) = 0.d0
        hact(nlvdi+1) = 0.d0

        if( bcolin ) then
          hact(1) = hact(1) - zmm               !FIXME
        end if

        do l=0,ilevel+1
          if( hact(l) .le. 0. ) then
            rhact(l) = 0.d0
          else
            rhact(l) = 1. / hact(l)
          end if
        end do

!-------------------------------------------------------------
! compute element averaged turbulent viscosity
!-------------------------------------------------------------

        k1 = nen3v(1,ie)
        k2 = nen3v(2,ie)
        k3 = nen3v(3,ie)
        do l=0,ilevel
            vis = vismol
            vis = vis + (visv(l,k1)+visv(l,k2)+visv(l,k3))/3.d0
            vis = vis + (vts(l,k1)+vts(l,k2)+vts(l,k3))/3.d0
            alev(l) = vis
        end do


!-------------------------------------------------------------
! start of vertical loop
!
! first set depth values
!
! hhi/hhip/hhim		thickness of i/i+1/i-1 layer
! uui/uuip/uuim		transport in x in i/i+1/i-1 layer
! vvi/vvip/vvim		transport in y in i/i+1/i-1 layer
!
! in case of a layer that does not exist (i-1 of first layer) give any
! ...value because the corrisponding a/b/c will be 0
!-------------------------------------------------------------

        do l=1,ilevel

        bfirst = l .eq. 1
        blast  = l .eq. ilevel
        
        lp = min(l+1,ilevel)
        lm = max(l-1,1)

        uui     = utlov(l,ie)
        uuip    = utlov(lp,ie)
        uuim    = utlov(lm,ie)

        vvi     = vtlov(l,ie)
        vvip    = vtlov(lp,ie)
        vvim    = vtlov(lm,ie)
        
        hhi     = hact(l)
        hhip    = hact(l+1)
        hhim    = hact(l-1)

!	------------------------------------------------------
!	set up contributions of vertical viscosity
!	------------------------------------------------------

        rhp = 0.d0
        rhm = 0.d0

        if( hhi .gt. 0. ) then          !may be redundant
          if( hhip .gt. 0. ) then       !lower interface
            rhp = 2.0 * alev(l) / ( hhi + hhip )
          end if
          if( hhim .gt. 0. ) then       !upper interface
            rhm = 2.0 * alev(l-1) / ( hhi + hhim )
          end if
        end if
        
!	aus = afact * alev(l)
!	aux = dt * at * aus

        aus = 1.
        aux = dt * at

        aa  = aux * rhact(l) * ( rhm + rhp )
        aat = aus * rhact(l) * ( rhm + rhp )
        bb  = aux * rhact(l+1) * rhp
        bbt = aus * rhact(l+1) * rhp
        cc  = aux * rhact(l-1) * rhm
        cct = aus * rhact(l-1) * rhm

!	------------------------------------------------------
!	boundary conditions for stress on surface and bottom
!	------------------------------------------------------

        ppx = 0.d0
        ppy = 0.d0

        if( bfirst ) then
          ppx = ppx - taux
          ppy = ppy - tauy
        end if
        if( blast ) then
          aa  = aa + dt * rfricv(ie)
          aat = aat + rfricv(ie)
        end if

!	------------------------------------------------------
!	implicit advective contribution
!	------------------------------------------------------

        uuadv = 0.d0
        uvadv = 0.d0
        vuadv = 0.d0
        vvadv = 0.d0

        aux = dt * radv * rdist

        if( aux .gt. 0. ) then  !implict treatment of non-linear terms

            uc = uui/hhi
            vc = vvi/hhi

            do ii=1,3
                k = nen3v(ii,ie)
                up = momentxv(l,k) / hhi
                vp = momentyv(l,k) / hhi
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

!	------------------------------------------------------
!	explicit contribution (non-linear, baroclinic, diffusion)
!	------------------------------------------------------
        
        xexpl = rdist * fxv(l,ie)
        yexpl = rdist * fyv(l,ie)

!	------------------------------------------------------
!	ppx/ppy is contribution on the left side of equation
!	ppx corresponds to -F^x_l in the documentation
!	ppy corresponds to -F^y_l in the documentation
!	------------------------------------------------------

        ppx = ppx + aat*uui - bbt*uuip - cct*uuim - gammat*vvi  & 
     &                  + grav*hhi*bz + (hhi/rowass)*bpres + xexpl + wavefx(l,ie)
        ppy = ppy + aat*vvi - bbt*vvip - cct*vvim + gammat*uui  &
     &                  + grav*hhi*cz + (hhi/rowass)*cpres + yexpl + wavefy(l,ie)

!	------------------------------------------------------
!	set up matrix A
!	------------------------------------------------------

        jv = l+l
        ju = jv-1

        rmat(locssp(ju,ju,ngl,mbb)) = 1.d0 + aa + uuadv
        rmat(locssp(jv,jv,ngl,mbb)) = 1.d0 + aa + vvadv
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

!	------------------------------------------------------
!	set up right hand side -F^x and -F^y 
!	------------------------------------------------------

        rvec(ju) = ppx
        rvec(jv) = ppy

!	------------------------------------------------------
!	set up H^x and H^y
!	------------------------------------------------------

        rvec(ngl+ju) = hhi              !ASYM_OPSPLT
        rvec(ngl+jv) = 0.d+0
        rvec(2*ngl+ju) = 0.d+0
        rvec(2*ngl+jv) = hhi

	end do

!-------------------------------------------------------------
! end of vertical loop
!-------------------------------------------------------------

!-------------------------------------------------------------
! solution of vertical system (we solve 3 systems in one call)
!-------------------------------------------------------------

        !call gelb(rvec,rmat,ngl,1,mbb,mbb,epseps,ier)
        !call dgelb(rvec,rmat,ngl,1,mbb,mbb,epseps,ier)
        call dgelb(rvec,rmat,ngl,3,mbb,mbb,epseps,ier)          !ASYM_OPSPLT

	if(ier.ne.0) then
	  call vel_matrix_error(it,ier,ie,ilevel,rvec,rmat,hact,alev)
	  stop 'error stop : sp256v'
	end if

!-------------------------------------------------------------
! compute u^hat (negative sign because ppx/ppy was -F^x/-F^y)
!-------------------------------------------------------------

        afix=1-iuvfix(ie)       !chao deb

        do l=1,ilevel
          utlnv(l,ie) = utlov(l,ie) - dt * rvec(2*l-1)*afix     !chao deb
          vtlnv(l,ie) = vtlov(l,ie) - dt * rvec(2*l)*afix       !chao deb
        end do

!-------------------------------------------------------------
! save contribution A^{-1} H^x and A^{-1} H^y
!-------------------------------------------------------------

        do l=1,ngl                                              !ASYM_OPSPLT
          ddxv(l,ie) = rvec(ngl+l)
          ddyv(l,ie) = rvec(2*ngl+l)
        end do

!-------------------------------------------------------------
! special information
!-------------------------------------------------------------

	if( ie .eq. 1 .and. barea0 .and. baroc .and. niter .le. 5 ) then  !$$BAROC_AREA0
	  write(6,*) 'sp256v: BAROC_AREA0 active '
	end if

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!******************************************************************

	subroutine vel_matrix_error(it,ier,ie,lmax,rvec,rmat,hact,alev)

        use check

	implicit none

	integer it,ier,ie,lmax
	double precision rmat(10*lmax)
	double precision rvec(6*lmax)
	double precision hact(0:lmax+1)
	double precision alev(0:lmax)

	integer ii,l,kk,ngl,mbb
	double precision auxaux(-2:2)

	integer locssp

	ngl=2*lmax
	mbb=2
	if(ngl.eq.2) mbb=1

	write(6,*) 'Error in inverting matrix (vertical system)'
	write(6,*) 'it, ier : ',it,ier
	write(6,*) 'ie,lmax,ngl,mbb: ',ie,lmax,ngl,mbb
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

	write(6,*) 'hact: ',(hact(l),l=0,lmax)
	write(6,*) 'alev: ',(alev(l),l=0,lmax)

	call check_set_unit(6)
	call check_elem(ie)
	call check_nodes_in_elem(ie)

	end

!******************************************************************

	subroutine hydro_transports_final

! post processing of time step
!
! written on 23.07.1997 by ggu  (from scratch)
!
	use internal
	use depth
	use hydro_baro
	use hydro_admin
	use evgeom
	use levels
	use basin
        use para
        use utility
        use time_util

	implicit none
!
! parameters
	include 'param.h'
	double precision drittl
	parameter (drittl=1./3.)
! common
	include 'mkonst.h'
	include 'pkonst.h'
	include 'femtime.h'

        integer afix            !chao deb

! local
	logical bcolin,bdebug
	integer ie,ii,l,kk,k
	integer ilevel
	integer ju,jv
	double precision dt,azpar,ampar
	double precision az,am,beta
	double precision bz,cz,um,vm,dz
	double precision du,dv
        character*40 filename,filename2,filename3,filename4

!-------------------------------------------------------------
! initialize
!-------------------------------------------------------------

	bcolin=iround(getpar('iclin')).ne.0	! linearized conti
	bdebug = .false.

	call get_timestep(dt)
	call getazam(azpar,ampar)
	az=azpar
	am=ampar

	beta = dt * grav * am 

!-------------------------------------------------------------
! start loop on elements
!-------------------------------------------------------------

	do ie=1,nel

	ilevel=ilhv(ie)

        afix=1-iuvfix(ie)       !chao deb

!	------------------------------------------------------
!	compute barotropic pressure term
!	------------------------------------------------------

        bz=0.
        cz=0.
        do ii=1,3
          kk=nen3v(ii,ie)
          dz = znv(kk) - zeov(ii,ie)
          bz = bz + dz * ev(ii+3,ie)
          cz = cz + dz * ev(ii+6,ie)
        end do

!	------------------------------------------------------
!	new transports from u/v hat variable
!	------------------------------------------------------


	do l=1,ilevel

	  jv=l+l
	  ju=jv-1

	  du = beta * ( ddxv(ju,ie)*bz + ddyv(ju,ie)*cz )	!ASYM_OPSPLT_CH
	  dv = beta * ( ddxv(jv,ie)*bz + ddyv(jv,ie)*cz )	!ASYM_OPSPLT_CH

	  utlnv(l,ie) = utlnv(l,ie) - du*afix   !chao deb
	  vtlnv(l,ie) = vtlnv(l,ie) - dv*afix   !chao deb

	end do

!	------------------------------------------------------
!	barotropic transports
!	------------------------------------------------------

	um = 0.
	vm = 0.
	do l=1,ilevel
	  um = um + utlnv(l,ie)
	  vm = vm + vtlnv(l,ie)
	end do
	unv(ie) = um
	vnv(ie) = vm

	end do

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------


	end

!******************************************************************

	subroutine hydro_vertical(dzeta)

! computes vertical velocities
!
! velocities are computed on S/T points (top and bottom of layer)
! bottom velocity of the whole column is assumed to be 0
! -> maybe change this
!
! computes volume differences and from these computes vertical
! velocities at every time step so that the computed velocities
! satisfy the continuity equation for every single time step
!
! wlnv is computed horizontally at a node and vertically
! it is at the center of the layer -> there are nlv velocities
! computed
!
! b,c are 1/m, (phi is dimensionless)
! aj is m**2
! utlnv... is m**2/s
! dvol is in m**3/s
! vv is m**2 (area)
!
! wlnv is first used to accumulate volume difference -> dvol
! at the end it receives the vertical velocity
!
! wlnv (dvol)   aux array for volume difference
! vv            aux array for area
!
! written on 27.08.91 by ggu  (from scratch)
! 14.08.1998	ggu	w = 0 at open boundary nodes
! 20.08.1998	ggu	some documentation

	use bnd_geom
	use geom_dynamic
	use bnd_dynamic
	use hydro_vel
	use hydro_admin
	use evgeom
	use levels
	use basin
	use shympi
        use elems_dealing
        use bndo_admin
        use para
        use time_util
        use mpi_io_admin

	implicit none

! arguments
	double precision dzeta(nkn)
! local
	logical debug
	integer k,ie,ii,kk,l,lmax,e
	integer ilevel
        integer ibc,ibtyp
	double precision aj,wbot,wdiv,ff,atop,abot,wfold
	double precision b,c
	double precision am,az,azt,azpar,ampar
	double precision ffn,ffo
	double precision volo,voln,dt,dvdt,q
	double precision dzmax,dz
	double precision, allocatable :: vf(:,:)
	double precision, allocatable :: va(:,:)
! statement functions

	!logical isein
        !isein(ie) = iwegv(ie).eq.0

! 2d -> nothing to be done

	dzeta = 0.
	wlnv = 0.
	if( nlvdi == 1 ) return

! initialize

	call getazam(azpar,ampar)
	az=azpar
	am=ampar
	azt = 1. - az
	call get_timestep(dt)

#ifdef DEBUGON
	allocate(vf(nlvdi,nkn_local),va(nlvdi,nkn_local))
#else
	allocate(vf(nlvdi,nkn),va(nlvdi,nkn))
#endif
	vf = 0.d0
	va = 0.d0

! compute difference of velocities for each layer
!
! f(ii) > 0 ==> flux into node ii
! aj * ff -> [m**3/s]     ( ff -> [m/s]   aj -> [m**2]    b,c -> [1/m] )

#ifdef DEBUGON
        call shympi_exchange_halo_3d_elems(utlov)
        call shympi_exchange_halo_3d_elems(vtlov)
        call shympi_exchange_halo_3d_elems(utlnv)
        call shympi_exchange_halo_3d_elems(vtlnv)
        call shympi_exchange_halo_2d_elems(ilhv)

	do e=1,nel_local
          if(bmpi) then
            ie=domain%elems%mapID(e)
          else
            ie=e
          endif
#else
	do ie=1,nel
#endif
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

        !call shympi_comment('shympi_elem: exchange vf,va')

#ifndef DEBUGON
        call shympi_exchange_and_sum_3D_nodes(vf)
        call shympi_exchange_and_sum_3D_nodes(va)
#endif

! from vel difference get absolute velocity (w_bottom = 0)
!	-> wlnv(nlv,k) is already in place !
!	-> wlnv(nlv,k) = 0 + wlnv(nlv,k)
! w of bottom of last layer must be 0 ! -> shift everything up
! wlnv(nlv,k) is always 0
!
! dividing wlnv [m**3/s] by area [vv] gives vertical velocity
!
! in va(l,k) is the area of the upper interface: a(l) = a_i(l-1)
! =>  w(l-1) = flux(l-1) / a_i(l-1)  =>  w(l-1) = flux(l-1) / a(l)

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

	!write(6,*) 'hydro_vertical: dzmax = ',dzmax

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

! set w to zero at open boundary nodes (new 14.08.1998)
!
! FIXME	-> only for ibtyp = 1,2 !!!!

	do k=1,nkn
            !if( is_external_boundary(k) ) then	!bug fix 10.03.2010
            if( is_zeta_bound(k) ) then
	      wlnv(:,k) = 0.
	      dzeta(k) = 0.
            end if
	end do

	deallocate(vf,va)

	end

!*******************************************************************

!---------------------------------------------------------------------------
        end module hydrodynamic
!---------------------------------------------------------------------------
