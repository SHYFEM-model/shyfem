!
! $Id: subn11.f,v 1.69 2010-03-22 15:29:31 georg Exp $
!
! boundary condition routines
!
! contents :
!
! subroutine sp111(mode)		set up boundary/initial cond.
!
! subroutine z_tilt			tilting of boundary surface (fixed)
! subroutine c_tilt			tilting of boundary surface (Coriolis)
!
! subroutine initilt(ibc)       	finds tilting node in node list
! subroutine iniflux(ibc)		initializes flux boundary
!
! subroutine set_mass_flux		sets up (water) mass flux array mfluxv
! subroutine make_scal_flux(what,r3v,scal,sflux,sconz,ssurf)	sets scal flux
! subroutine flux_debug(what,mfluxv,sflux,sconz)
! subroutine check_scal_flux(what,scal,sconz)			checks scal flux
! 
! subroutine mult_scal_bc(r3v,value)	multiplies array for scalar BC by value
! 
! subroutine dist_horizontal(nlvddi,r3v,n,value)
! subroutine aver_horizontal(nlvddi,r3v,n,value)
!
! subroutine print_scal_bc(r3v)		prints non-flag entries of scalar BC
! subroutine get_bflux(k,flux)		returns boundary flux of node k
! 
! subroutine level_flux(it,levflx,kn,rw)compute discharge from water level
! subroutine z_smooth(z)		smooths z values
! subroutine flow_out_piece_new(z,rout)
!
! revision log :
!
! revised 05.08.92 by ggu   $$ibtyp3 - implementation of ibtyp=3
! revised 31.08.92 by ggu   $$rqv1   - rqv initialized in sp159b, not here
! revised 05.09.92 by ggu   $$close1 - rqv initialized here
! revised 27.10.93 by ggu   $$roger  - ibtyp=70 (nwe-shelf)
! revised 05.11.93 by ggu   $$roger  - call to roger has been commented
! revised 11.01.94 by ggu   $$restart - restart is reading all variables
! revised 20.01.94 by ggu   $$zeov - initialize zeov
! revised 20.01.94 by ggu   $$conz - impl. of conz with bnd(12,.)
! revised 31.01.94 by ggu   $$conzbc - impl. of bc for conz with rcv
! revised 24.03.94 by ggu   $$surel - impl. of distributed source/sink
! revised 02.04.96 by ggu   $$exxqq - interpolation for more files ok
! revised 25.06.97 by ggu   complete restructured -> new subroutines
! revised 18.09.97 by ggu   $$FLUX3 - special treatment of type 3 boundary
! revised 23.09.97 by ggu   concentration boundary as file implemented
! revised 03.12.97 by ggu   $$TS - temp/salt implemented
! 20.05.1998	ggu	rain from file implemented
! 22.05.1998	ggu	local variable t introduced
! 22.05.1998	ggu	corrected bug for rain (surel called twice)
! 28.05.1998    ggu     new ruv, rvv (momentum input) (incomplete)
! 20.06.1998    ggu     more on momentum input (mvalue)
! 29.06.1998    ggu     !$$momin1 - bug fix for momentum input
! 13.07.1998    ggu     !initialize ruv, rvv by node
! 14.07.1998    ggu     finally momentum input finished -> new exxqq
! 20.08.1998    ggu     iextpo finally eliminated
! 20.08.1998    ggu     2d dependent routines copied to subini.f
! 21.08.1998    ggu     xv eliminated
! 24.08.1998	ggu	BC for concentration is bnd(20,..)
! 24.08.1998	ggu	BC for maximum input level is bnd(12,..) -> levmax
! 05.11.1998	ggu	slightly restructured restart
! 18.02.1999	ggu	allow for boundary type 0
! 20.04.1999	ggu	converted to stress instead of wind (tauxnv...)
! 01.12.1999	ggu	handle negative boundary type
! 07.05.2001	ggu	introduced variable zfact
! 24.10.2001	ggu	ignore boundary with type 0
! 10.08.2003	ggu	better commented, new routines meteo_init, meteo_force
! 14.08.2003	ggu	initialization of z and uv made explicit
! 10.03.2004	ggu	RQVDT - value in rqv is now discharge [m**3/s]
! 03.09.2004	ggu	restart taken out to ht
! 11.10.2004	ggu	new ibtyp=31, multiply with zfact also for sin (ZFACT)
! 11.03.2005	ggu	new boundary routines b3dvalue, c3dvalue (3D bounds)
! 17.02.2006	ggu	new routine get_bflux()
! 23.03.2006    ggu     changed time step to double precision
! 18.09.2007    ggu     new set_mass_flux
! 25.09.2007    ggu     routines deleted: [mbc]value, testbc
! 02.10.2007    ggu     bug fix in make_scal_flux: surface flux only in layer 1
! 08.10.2007    ggu     deleted commented lines
! 17.03.2008    ggu     name of some variables changed, new scalar values
! 07.04.2008    ggu     deleted c2dvalue, set_scal_bc
! 07.04.2008    ggu     differential input introduced (-5555)
! 09.04.2008    ggu     only level boundary array, re-arranged
! 10.04.2008    ggu&ccf new call to init_z0b()
! 17.04.2008    ggu     lmin introduced in set_mass_flux (negative levmax)
! 17.04.2008    ggu     evaporation introduced, rain adjusted
! 18.04.2008    ggu     rain simplified, bugfix
! 22.04.2008    ggu     in make_scal_flux do not alter mfluxv (parallel code)
! 03.06.2008    ggu     levmin introduced
! 24.03.2009    ggu     bug fix for rain; rain0d; new rain2distributed()
! 31.03.2009    ggu     in make_scal_flux() bug fix for negative fluxes
! 02.04.2009    ggu     use get_bnd_(i)par() for special routines
! 03.04.2009    ggu     set intpol depending on ibtyp if intpol==0
! 20.04.2009    ggu     new routine z_tilt (tilting around a fixed value)
! 27.04.2009    ggu     can use ktilt with ztilt
! 19.01.2010    ggu     in make_scal_flux() return also sconz at bound
! 26.02.2010    ggu     bug fix in make_scal_flux() - sconz was out of l-loop
! 01.03.2010    deb     tramp introduced
! 10.03.2010    ggu     new routine check_scal_flux()
! 22.03.2010    ggu     in make_scal_flux() forgotten to initialize sconz
! 14.04.2010    ggu     account for negative levmin/max
! 16.02.2011    ggu     copied meteo routines to submet.f
! 23.02.2011    ggu     new parameters tramp and levflx implemented
! 01.03.2011    ggu     implement smoothing for levflx
! 14.05.2011    ggu     new routine get_discharge()
! 29.11.2013    ggu     prepared for ibtyp == 0
! 10.07.2014    ggu     only new file format allowed
! 30.10.2014    ggu     in c_tilt() compute distance also for lat/lon
! 31.10.2014    ccf     new call to init_z0 instead than init_z0b
! 07.11.2014    ggu     bug fix for distance computation in z_tilt, c_tilt
! 10.02.2015    ggu     new call to iff_read_and_interpolate()
!
!***************************************************************
!----------------------------------------------------------------------
        module bnd_routines
!----------------------------------------------------------------------
        contains
!----------------------------------------------------------------------

        subroutine sp111(mode)

! set up boundary and initial conditions
!
! mode		1 : first call, initialize b.c.
!		2 : read in b.c.

        use shympi
        use bnd_geom
        use bnd_aux
        use bnd_dynamic
        use levels, only : nlvdi,nlv
        use basin, only : nkn,nel,ngr,mbw
        use intp_fem_file
        use para
        use fem_util
        use bnd_admin
        use elems_dealing
        use bndo_admin
        use meteo_admin
        use transforms
        use initialize
        use time_util

        implicit none

        integer mode

        include 'femtime.h'
        include 'mkonst.h'

        double precision rwv2(nkn)
        integer, save, allocatable :: ids(:)            !id values for BC

        integer nodes(nkn)
        double precision vconst(nkn)

        character*10 auxname
        logical bimpose
        integer kranf,krend,k,kn
        integer ibc,ibtyp
        integer nk,i,kk,kindex,iv
        integer nsize,il
        integer iunrad,ktilt
        integer ip,l,lmax,ivar
        integer levflx
        integer nbc
        integer id,intpol,nvar,ierr
        double precision dtime0,dtime
        double precision rw,const,aux
        double precision dt
        double precision conz,temp,salt
        double precision conzdf,tempdf,saltdf
!	double precision dz
        double precision rmu,rmv
        double precision rwint
        double precision conz3,temp3,salt3
        double precision tramp,alpha
        character*80 zfile

        logical rst_use_restart         !ivb

	!tramp = 0.	!is now handled in str
	!tramp = 86400. !DEB

        call get_timestep(dt)

        if( mode .eq. 1 ) goto 1
        if( mode .eq. 2 ) goto 2
        stop 'error stop: internal error sp111 (1)'

!---------------------------------------------------------------
! first call %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------

    1	continue

!	-----------------------------------------------------
!       initialize boundary ids
!	-----------------------------------------------------

        nbc = nbnds()
        allocate(ids(nbc))
        ids = 0
        
!	-----------------------------------------------------
!       initialize meteo
!	-----------------------------------------------------

        call meteo_init

!	-----------------------------------------------------
!       initialize tilted and flux boundaries
!	-----------------------------------------------------

        nbc = nbnds()

        do ibc=1,nbc
          ibtyp=itybnd(ibc)
          nk = nkbnds(ibc)

          bimpose = ibtyp .ge. 1 .and. ibtyp .le. 2

          if(shympi_partition_on_elements()) then
            if( bimpose .and. nk .le. 0 ) goto 95 !$$ibtyp3
          else
            if( bimpose .and. nk .le. 1 ) goto 95 !$$ibtyp3
          end if

          if( bimpose ) then    		!$$FLUX3 - not for ibtyp=3
            call initilt(ibc)           	!find nodes for tilting
            call iniflux(ibc)   		!set up rlv,rhv,rrv,ierv
          end if
        end do

!	-----------------------------------------------------
!       initialization of aux array
!	-----------------------------------------------------

!	-----------------------------------------------------
!       initialization of fem_intp
!	-----------------------------------------------------

        dtime0 = itanf
        nvar = 1
        vconst = 0.
        ids = 0

        do ibc=1,nbc
          nk = nkbnds(ibc)
          if( nk .le. 0 ) cycle
          do i=1,nk
            nodes(i) = kbnds(ibc,i)
          end do
          call get_boundary_file(ibc,'zeta',zfile)
          call get_bnd_ipar(ibc,'ibtyp',ibtyp)
          call get_bnd_ipar(ibc,'intpol',intpol)
          if( intpol .le. 0 ) then
            intpol = 2
            if( ibtyp .eq. 1 ) intpol = 4
          end if
          if (bmpi) then
            call iff_init_mpi(dtime0,zfile,nvar,bounds%tnob(ibc),nk,0,intpol,nodes,vconst,id,ibc)
          else
            call iff_init(dtime0,zfile,nvar,nk,0,intpol,nodes,vconst,id)
          end if
          if( ibtyp .le. 0 ) then
            write(auxname,'(a6,1x,i3)') 'closed',ibtyp
          else if( ibtyp .eq. 1 ) then
            write(auxname,'(a6,1x,i3)') 'zlevel',ibtyp
          else if( ibtyp .eq. 2 .or. ibtyp .eq. 3 ) then
            write(auxname,'(a6,1x,i3)') 'disch ',ibtyp
          else
            write(auxname,'(a6,1x,i3)') 'bound ',ibtyp
          end if
          call iff_set_description(id,ibc,auxname)
          ids(ibc) = id
          il = len_trim(zfile)
          write(6,*) 'boundary file opened: ',ibc,id,zfile(1:il)
        end do

	!call iff_print_info(ids(1))

!	-----------------------------------------------------
!       determine constant for z initialization
!	-----------------------------------------------------

        const=getpar('const')   !constant initial z value
        dtime = itanf
        ivar = 1
        lmax = 1

        do ibc=1,nbc
          ibtyp=itybnd(ibc)
          nk = nkbnds(ibc)
          id = ids(ibc)
          if( id .le. 0 ) cycle

          if(const.eq.flag.and.ibtyp.eq.1) then
                call iff_read_and_interpolate(id,dtime)
                call iff_time_interpolate(id,dtime,ivar,nk,lmax,rwv2)
                call adjust_bound(id,ibc,it,nk,rwv2)
                if(bounds%rankfirst.eq.my_id) const = rwv2(1)
          end if

          if(ibtyp.eq.70) then  !nwe-shelf	!$$roger - special b.c.
!	    call roger(rzv,dt,0)
!	    call roger(rzv,irv,nrb,dt,0)
          end if
        end do

        if(const.eq.flag) const=0.
        if(bmpi) then
          bounds%rankfirst = shympi_max(bounds%rankfirst)
          call shympi_bcast(const,bounds%rankfirst)
        end if
        call putpar('const',const)

!	-----------------------------------------------------
!       initialize variables or restart
!	...the variables that have to be set are zenv, utlnv, vtlnv
!	-----------------------------------------------------

        call init_z(const)      !initializes zenv
        call set_area           !initializes areakv
        call make_new_depth     !initializes layer thickness
        call init_uvt           !initializes utlnv, vtlnv
        call init_z0            !initializes surface and bottom roughness

        call uvint
        call copy_uvz

!	-----------------------------------------------------
!       finish
!	-----------------------------------------------------

! next only for radiation condition
!
!        write(78,*) 46728645,1
!        write(78,*) 50,0
!        write(78,*) 0,50,-0.01,0.01
!        write(79,*) 46728645,1
!        write(79,*) 50,0
!        write(79,*) 0,50,-0.01,0.01

        return

!---------------------------------------------------------------
! normal call for boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------

    2	continue

!	-----------------------------------------------------
!	initialize node vectors with boundary conditions
!	-----------------------------------------------------

        do k=1,nkn
          rzv(k)=flag
          rqv(k)=0.     !$$rqv1 !$$close1	[m**3/s]
          rqpsv(k)=0.   !fluxes - point sources [m**3/s]
          rqdsv(k)=0.   !fluxes - distributed sources through surface [m**3/s]
          ruv(k)=0.     !momentum input
          rvv(k)=0.
        end do

!	-----------------------------------------------------
!	loop over boundaries
!	-----------------------------------------------------

        call bndo_radiat(it,rzv)
        nbc = nbnds()

        dtime = it
        ivar = 1
        lmax = 1

        do ibc=1,nbc

          call get_bnd_ipar(ibc,'ibtyp',ibtyp)
          call get_bnd_ipar(ibc,'levflx',levflx)
          call get_bnd_par(ibc,'tramp',tramp)
          id = ids(ibc)

          if( ibtyp .le. 0 ) cycle
          if( id .le. 0 ) cycle

          nk = nkbnds(ibc)   !total number of nodes of this boundary

          rmu = 0.
          rmv = 0.

          call iff_read_and_interpolate(id,dtime)
          call iff_time_interpolate(id,dtime,ivar,nk,lmax,rwv2)
          call adjust_bound(id,ibc,it,nk,rwv2)

          alpha = 1.
          if( tramp .gt. 0. .and. it-itanf .le. tramp ) then
             alpha = (it-itanf) / tramp
          end if

          do i=1,nk

             kn = kbnds(ibc,i)
             rw = rwv2(i)

             if(ibtyp.eq.1) then		!z boundary
               rzv(kn)=rw
             else if(ibtyp.eq.2) then		!q boundary
               call level_flux(it,levflx,kn,rw)	!zeta to flux
               kindex = kbndind(ibc,i)
               rqpsv(kn)=alpha*rw*rrv(kindex)	!BUGFIX 21-08-2002, RQVDT
             else if(ibtyp.eq.3) then		!$$ibtyp3 - volume inlet
               rqpsv(kn) = alpha*rw
             else if(ibtyp.eq.4) then		!momentum input
               ruv(kn) = rmu
               rvv(kn) = rmv
             else if(ibtyp.eq.31) then		!zero gradient for z
               ! already done...
               call get_bnd_ipar(ibc,'ktilt',ktilt)
               rzv(ktilt)=rw
               if( ibc .eq. 1 ) then
                 iunrad = 78
               else
                 iunrad = 79
               end if
               if( i .eq. 1 ) write(iunrad,*) it,nk,ktilt,rw
               write(iunrad,*) rzv(kn)
             else if(ibtyp.eq.32) then		!for malta 
!	       nothing
             else if(ibtyp.eq.0) then		!switched off
!	       nothing
             else if(ibtyp.lt.0) then		!closed
!	       nothing
             else
!               kranf,krend not available...
!	       call zspeci(ibtyp,kranf,krend,rw)	!for radiation...
               write(6,*) 'boundary = ',ibc,'   type = ',ibtyp
               stop 'error stop sp111: Unknown boundary type'
             end if

          end do

        end do

	!write(99,*) '======================================'
	!write(99,*) ' it = ',it
	!write(99,*) '======================================'
	!call iff_print_info(12,0,.true.)

!	-----------------------------------------------------
!	tilting
!	-----------------------------------------------------

        call z_tilt
        call c_tilt

!	-----------------------------------------------------
!	meteo forcing					!$$surel
!	-----------------------------------------------------

        call meteo_force

!	-----------------------------------------------------
!	set mass flux -> fills mfluxv and integrates to rqv
!	-----------------------------------------------------

        call set_mass_flux

!	-----------------------------------------------------
!	testing
!	-----------------------------------------------------

!	call tsbnds

! -----------------------------------------------------------
! end of routine
! -----------------------------------------------------------


        return
   95	continue
        write(6,*) 'One node boundary not allowed'
        write(6,*) 'Boundary :',ibc
        write(6,*) 'type     :',ibtyp,nk
        stop 'error stop : sp111'
        end

!**************************************************************

	subroutine z_tilt

! artificial tilting of boundary surface - uses ktilt and ztilt
!
! the first boundary node is set to -ztilt, and the last to +ztilt
! the total water level difference is therefore 2*ztilt
! if ktilt is not given then the other nodes are linearily interpolated
!	between these two values
! if ktilt is given then this node will be set to z=0 and the other
!	nodes are linearly interpolated between start-ktilt and ktilt-end

	use bnd_geom
	use bnd_dynamic
	use basin
        use shympi
        use bnd_admin
        use evgeom

	implicit none

	include 'femtime.h'

	integer ibc,ibtyp,ktilt
	integer nbc
	integer k,kranf,krend,kn1,kn2
	double precision dx,dy,ztilt,z
	double precision rltot,rltot1,rltot2,rl

	nbc = nbnds()

	do ibc=1,nbc
          ibtyp = itybnd(ibc)
	  if( ibtyp <= 0 ) cycle
          call kanfend(ibc,kranf,krend)
	  call get_bnd_par(ibc,'ztilt',ztilt)
	  call get_bnd_ipar(ibc,'ktilt',ktilt)
	  if( ztilt .ne. 0. .and. ibtyp .eq. 1 ) then
	    rltot = 0.
	    rltot1 = 0.
	    do k=kranf+1,krend
              if(bmpi) then
                kn2=bounds%niob(k)
                kn1=bounds%niob(k-1)
              else
		kn2=irv(k)
		kn1=irv(k-1)
              end if
		call compute_distance(xgv(kn1),ygv(kn1),xgv(kn2),ygv(kn2),dx,dy)
	        rltot = rltot + sqrt(dx*dx+dy*dy)
		if( k .eq. ktilt ) rltot1 = rltot	!BUG 3.12.2013
	    end do
	    rltot2 = rltot - rltot1

! in rltot the whole length of boundary is stored
! in rltot1 and rltot2 are first and second part of length of boundary
! rltot1 from start to ktilt, and rltot2 from ktilt to end
! if no ktilt is given rltot1/rltot2 are not used

	    rl = 0.
            if(bmpi) then
              rzv(bounds%niob(kranf)) = rzv(bounds%niob(kranf)) - ztilt
            else
              rzv(irv(kranf)) = rzv(irv(kranf)) - ztilt
            end if
	    do k=kranf+1,krend
              if(bmpi) then
                kn2=bounds%niob(k)
                kn1=bounds%niob(k-1)
              else
		kn2=irv(k)
		kn1=irv(k-1)
              end if
		call compute_distance(xgv(kn1),ygv(kn1),xgv(kn2),ygv(kn2),dx,dy)
	        rl = rl + sqrt(dx*dx+dy*dy)
		if( ktilt .le. 0 ) then			!no fixed node
	          z = - ztilt + (rl/rltot) * 2. * ztilt
		else
		  if( rl .le. rltot1 ) then
	            z = - ztilt + (rl/rltot1) * ztilt
		  else
	            z = ((rl-rltot1)/rltot2) * ztilt
		  end if
		end if
                rzv(kn2) = rzv(kn2) + z
	    end do
	  end if
	end do

	end

!**************************************************************

	subroutine c_tilt

! tilting of boundary surface due to Coriolis acceleration - needs ktilt
!
! if ztilt is given then z_tilt() is used to tilt water level
! if ktilt is not given nothing is tilted

	use bnd_geom
	use bnd_dynamic
	use hydro_print
	use hydro_admin
	use basin
        use shympi
        use para
        use bnd_admin
        use meteo_admin
        use evgeom

	implicit none

	include 'femtime.h'
	include 'pkonst.h'
	include 'mkonst.h'

	integer ibc,ibtyp,kranf,krend,ktilt,k,kn1,kn2
	integer nbc
	double precision roinv,f,ginv,dx,dy,taux,tauy
	double precision taux1,taux2,tauy1,tauy2,wx,wy
	double precision u,v,z,h,p1,p2,b,hh,ztilt

	roinv=1./rowass
	f=fcor
	ginv=1./grav

	nbc = nbnds()

	do ibc=1,nbc
         ibtyp = itybnd(ibc)
	 if( ibtyp <= 0 ) cycle
         call kanfend(ibc,kranf,krend)
	 call get_bnd_ipar(ibc,'ktilt',ktilt)
	 call get_bnd_par(ibc,'ztilt',ztilt)

	 if( ztilt .ne. 0 ) then
		!nothing
	 else if(ktilt.gt.0.and.ibtyp.eq.1) then
	   do k=ktilt+1,krend,1
              if(bmpi) then
                kn2=bounds%niob(k)
                kn1=bounds%niob(k-1)
              else
		kn2=irv(k)
		kn1=irv(k-1)
              end if
		call compute_distance(xgv(kn1),ygv(kn1),xgv(kn2),ygv(kn2),dx,dy)
		call get_meteo_forcing(kn1,wx,wy,taux1,tauy1,p1)
		call get_meteo_forcing(kn2,wx,wy,taux2,tauy2,p2)
		taux = 0.5*(taux1+taux2)
		tauy = 0.5*(tauy1+tauy2)
		u=(up0v(kn1)+up0v(kn2))*.5
		v=(vp0v(kn1)+vp0v(kn2))*.5
		z=(znv(kn1)+znv(kn2))*.5
		h=(rhv(k-1)+rhv(k))*.5
		b=1./(h+z)
		hh=-roinv*(p2-p1)+b*(taux*dx+tauy*dy)+f*(v*dx-u*dy)
		rzv(kn2)=rzv(kn1)+ginv*hh
	   end do

	   do k=ktilt-1,kranf,-1
              if(bmpi) then
                kn2=bounds%niob(k+1)
                kn1=bounds%niob(k)
              else
		kn2=irv(k+1)
		kn1=irv(k)
              end if
		call compute_distance(xgv(kn1),ygv(kn1),xgv(kn2),ygv(kn2),dx,dy)
		call get_meteo_forcing(kn1,wx,wy,taux1,tauy1,p1)
		call get_meteo_forcing(kn2,wx,wy,taux2,tauy2,p2)
		taux = 0.5*(taux1+taux2)
		tauy = 0.5*(tauy1+tauy2)
		u=(up0v(kn1)+up0v(kn2))*.5
		v=(vp0v(kn1)+vp0v(kn2))*.5
		z=(znv(kn1)+znv(kn2))*.5
		h=(rhv(k+1)+rhv(k))*.5
		b=1./(h+z)
		hh=-roinv*(p2-p1)+b*(taux*dx+tauy*dy)+f*(v*dx-u*dy)
		rzv(kn1)=rzv(kn2)-ginv*hh
	   end do
	 end if
	end do

	end

!**************************************************************

	subroutine initilt(ibc)

! finds tilting node in boundary node list

        use fem_util
        use bnd_admin

	implicit none

	integer ibc

	logical berr
	integer kranf,krend
	integer ktilt,i
        integer kb
	integer ibtyp

        ibtyp = itybnd(ibc)
	if( ibtyp <= 0 ) return

	call get_bnd_ipar(ibc,'ktilt',ktilt)
	if(ktilt.le.0) return

	call kanfend(ibc,kranf,krend)

	berr = .true.
	do i=kranf,krend
           kb = kbnd(i)
	   if( kb .eq. ktilt ) then
		call set_bnd_ipar(ibc,'ktilt',i)
		berr = .false.
	   end if
	end do

	if( berr ) then
	  write(6,*) 'Node number for tilting not in boundary node list'
	  write(6,*) 'ktilt :',ipext(ktilt)
	  stop 'error stop : initilt'
	end if

	end

!******************************************************************

	subroutine iniflux(ibc)

! initializes flux boundary

	use bnd_geom
	use basin
        use shympi
        use fem_util
        use bnd_admin

	implicit none

        integer ibc

	integer kranf,krend
	integer ie,i,k1,k2,kk1,kk2,ii1,ii2
	integer ibtyp
	double precision fm,dx,dy,rl,h1,h2,fl

        ibtyp = itybnd(ibc)
	if( ibtyp <= 0 ) return

        call kanfend(ibc,kranf,krend)

	if( krend-kranf .le. 0 ) return

	ii1 = 0
	ii2 = 0
	kk1 = 0
	kk2 = 0

	fm=0
	do i=kranf,krend
	   rrv(i)=0.
	   rhv(i)=0.
	end do

	do i=kranf,krend-1
          if(bmpi) then
	    if((bounds%last(i).ne.0)) cycle
	    k1=bounds%niob(i)
	    k2=bounds%niob(i+1)
          else
            k1 = irv(i)
            k2 = irv(i+1)
          end if

	  do ie=1,nel
	   do ii1=1,3
	      ii2=mod(ii1,3)+1
	      kk1=nen3v(ii1,ie)
	      kk2=nen3v(ii2,ie)
	      if(k1.eq.kk1.and.k2.eq.kk2) goto 33
	   end do
	  end do

   33	  continue

	  if( k1 .ne. kk1 .or. k2 .ne. kk2 ) then
	write(6,*) 'Cannot locate boundary nodes in element index',my_id
	write(6,*) 'node 1,node 2 :',k1,k2,ipext(k1),ipext(k2),my_id
	write(6,*) 'node 1,node 2 :',k1,domain%nodes%globalID(k1)
	write(6,*) 'node 1,node 2 :',k2,domain%nodes%globalID(k2)
	write(6,*) '(Are you sure that boundary nodes are given',my_id
	write(6,*) '   in anti-clockwise sense ?)',my_id
	stop 'error stop : iniflux'
	  end if

	  dx=xgv(k1)-xgv(k2)
	  dy=ygv(k1)-ygv(k2)
	  rl=sqrt(dx*dx+dy*dy)
	  h1=hm3v(ii1,ie)
	  h2=hm3v(ii2,ie)
	  fl=rl*(h1+h2)/2.
	  fm=fm+fl

	  rrv(i)=rrv(i)+rl*(2.*h1+h2)/6.
	  rrv(i+1)=rrv(i+1)+rl*(h1+2.*h2)/6.

	  rlv(i)=rl

	  rhv(i)=rhv(i)+h1
	  rhv(i+1)=rhv(i+1)+h2

	  ierv(1,i)=ie
	  ierv(2,i)=mod(ii2,3)+1

	end do

	do i=kranf,krend
	   rrv(i)=rrv(i)/fm
	end do
	do i=kranf+1,krend-1
	   rhv(i)=rhv(i)/2.
	end do

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine adjust_bound(id,ibc,it,nk,rw)

	use intp_fem_file
        use bnd_admin

	implicit none

	integer id
	integer ibc
	integer it
	integer nk
	double precision rw(nk)

	integer i
	double precision rit,rw0,zfact

	call get_bnd_par(ibc,'zfact',zfact)

	if( iff_has_file(id) ) then
	  do i=1,nk
	    rw(i) = rw(i) * zfact
	  end do
	else
	  rit = it
	  call get_oscil(ibc,rit,rw0)
	  do i=1,nk
	    rw(i) = rw0 * zfact
	  end do
	end if

	end

!*******************************************************************

        subroutine set_mass_flux

! sets up (water) mass flux array mfluxv (3d) and rqv (vertically integrated)

        use bnd_dynamic
        use levels
        use bnd_admin
        use elems_dealing
        use basin, only : nkn,nel,ngr,mbw

        implicit none

        logical debug
        integer i,k,l,lmin,lmax,nk,ibc,mode
        integer ibtyp,levmax,levmin
        integer nbc
        double precision flux,vol,voltot,fluxtot,fluxnode
        double precision vols(nkn)

!------------------------------------------------------------------
! initialize arrays and parameter
!------------------------------------------------------------------

        mode = -1               !old time step
        debug = .true.
        debug = .false.

        do k=1,nkn
          do l=1,nlv
            mfluxv(l,k) = 0.d0
          end do
        end do

!------------------------------------------------------------------
! loop over boundaries for point sources -> use volumes as weight
!------------------------------------------------------------------

        nbc = nbnds()

        do ibc=1,nbc

          nk = nkbnds(ibc)
          call get_bnd_ipar(ibc,'ibtyp',ibtyp)
          call get_bnd_ipar(ibc,'levmin',levmin)
          call get_bnd_ipar(ibc,'levmax',levmax)

          if( ibtyp .lt. 2 .or. ibtyp .gt. 3 ) nk = 0		!skip

          if(debug) then
            write(6,*) 'computing mass flux: ',ibc,ibtyp,nk,levmax
          end if

          do i=1,nk
            k = kbnds(ibc,i)
            lmax = ilhkv(k)
            if( levmax .gt. 0 ) lmax = min(lmax,levmax)
            if( levmax .lt. 0 ) lmax = min(lmax,lmax+1+levmax)
            lmin = 1
            if( levmin .gt. 0 ) lmin = max(lmin,levmin)
            if( levmin .lt. 0 ) lmin = max(lmin,lmax+1+levmin)

            if( lmin .gt. lmax ) goto 98

            voltot = 0.d0
            do l=lmin,lmax
              vol = volnode(l,k,mode)
              vols(l) = vol
              voltot = voltot + vol
            end do

            if( voltot .le. 0. ) goto 99

            flux = rqpsv(k)
            if(debug) write(6,*) '   ',k,lmin,lmax,flux,voltot
            do l=lmin,lmax
              mfluxv(l,k) = flux * vols(l) / voltot
            end do
          end do

        end do

!------------------------------------------------------------------
! add distributed sources
!------------------------------------------------------------------

        do k=1,nkn
          mfluxv(1,k) = mfluxv(1,k) + rqdsv(k)  !rain, evaporation
	  !lmax = ilhkv(k)
	  !mfluxv(lmax,k) = gwf		!here distributed ground water flow
        end do

!------------------------------------------------------------------
! compute total flux for check and integrate flux into rqv
!------------------------------------------------------------------

        fluxtot = 0.d0
        do k=1,nkn
          fluxnode = 0.d0
          lmax = ilhkv(k)
          do l=1,lmax
            flux = mfluxv(l,k)
	    !if( debug .and. flux .gt. 0. ) write(6,*) '  flux: ',k,l,flux
            fluxnode = fluxnode + flux
          end do
          rqv(k) = fluxnode
          fluxtot = fluxtot + fluxnode
        end do

        if( debug ) write(6,*) '  total flux: ',fluxtot

!------------------------------------------------------------------
! end of routine
!------------------------------------------------------------------

        return
   98	continue
        write(6,*) 'lmin > lmax'
   99	continue
        write(6,*) 'ibc = ',ibc
        write(6,*) 'i = ',i
        write(6,*) 'k = ',k
        write(6,*) 'ilhkv(k) = ',ilhkv(k)
        write(6,*) 'levmin = ',levmin
        write(6,*) 'levmax = ',levmax
        write(6,*) 'lmin = ',lmin
        write(6,*) 'lmax = ',lmax
        stop 'error stop set_mass_flux: voltot = 0'
        end

!**********************************************************************

	subroutine make_scal_flux(what,r3v,scal,sflux,sconz,ssurf)

! computes scalar flux from fluxes and concentrations

	use bnd_dynamic
	use levels
	use basin, only : nkn,nel,ngr,mbw
        use para

	implicit none

	character*(*) what
	double precision r3v(nlvdi,nkn)	!concentration for boundary condition
	double precision scal(nlvdi,nkn)	!concentration of scalar
	double precision sflux(nlvdi,nkn)	!mass flux for each finite volume (return)
	double precision sconz(nlvdi,nkn)	!concentration for each finite volume (return)
	double precision ssurf		!value of scalar for surface flux

	include 'mkonst.h'


	integer k,l,lmax,ks
	double precision flux,conz
	double precision surf_flux

	ks = 2827
	ks = 2757
	ks = 2831
	ks = -1
	!ks = nint(getpar('kref'))	!not working - here global, but local

	do k=1,nkn
	  lmax = ilhkv(k)
	  surf_flux = rqdsv(k)
	  do l=1,lmax
	    sflux(l,k) = 0.
	    sconz(l,k) = 0.
	    flux = mfluxv(l,k)
	    if( l .eq. 1 ) flux = flux - surf_flux	!without surface flux
	    conz = r3v(l,k)
	    if( flux .ne. 0. .or. conz .ne. flag ) then
	      if( flux .ne. 0. .and. conz .eq. flag ) goto 99
	      if( conz .le. -990. ) conz = scal(l,k)	!ambient value
	      if( conz .le. -5555. ) conz = scal(l,ks) - 10000. - conz !diff
	      if( flux .lt. 0. ) conz = scal(l,k)	!ambient value (BUG)
	      sflux(l,k) = flux * conz
	      sconz(l,k) = conz			!bug fix - was out of loop
	    end if
	  end do
	  conz = ssurf
	  if( ssurf .le. -990 ) conz = scal(1,k)
	  if( ssurf .le. -5555 ) conz = scal(1,k) - 10000. - ssurf !diff
	  sflux(1,k) = sflux(1,k) + surf_flux * conz
	  ! next should be sconz(1,k) = conz if surf_flux is eliminated
	  !sconz(1,k) = 0.
	  if( mfluxv(1,k) .ne. 0 ) then
	    sconz(1,k) = sflux(1,k) / mfluxv(1,k)
	  end if
	end do

	!if( what .eq. 'salt' ) then
	!  call flux_debug(what,mfluxv,sflux,sconz)
	!end if

	return
   99	continue
	write(6,*) 'what: ',what
	write(6,*) 'k,l: ',k,l
	write(6,*) 'flux,conz: ',flux,conz
	stop 'error stop make_scal_flux: boundary condition mismatch'
	end

!**********************************************************************

	subroutine flux_debug(what,mfluxv,sflux,sconz)

	use levels
        use defnames
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) what
	double precision mfluxv(nlvdi,nkn)	!mass flux
	double precision sflux(nlvdi,nkn)	!scalar flux
	double precision sconz(nlvdi,nkn)	!concentration for each finite volume

	include 'femtime.h'

	integer k,l,lmax
	double precision qtot,stot

	integer iunit
	save iunit
	data iunit /0/

	if( iunit .eq. 0 ) then
	  iunit = ifemop('.ggg','formatted','unknown')
	end if

	do k=1,nkn
	  lmax = ilhkv(k)
	  stot = 0.
	  qtot = 0.
	  do l=1,lmax
	    qtot = qtot + mfluxv(l,k)
	    stot = stot + mfluxv(l,k) * sconz(l,k)
	  end do
	  if( qtot .ne. 0 ) then
	    write(iunit,*) it,k,qtot,stot
	  end if
	end do

	end

!**********************************************************************

	subroutine check_scal_flux(what,scal,sconz)

! checks scalar flux

	use bnd_dynamic
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) what
	double precision scal(nlvdi,nkn)	!concentration of scalar
	double precision sconz(nlvdi,nkn)	!concentration for each finite volume

	include 'mkonst.h'
	include 'femtime.h'


	integer k,l,lmax,ks
	double precision cconz,qflux,mflux

	write(46,*) 'check_scal_flux ',what,it

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cconz = sconz(l,k)         !concentration has been passed
            qflux = mfluxv(l,k)
            if( qflux .lt. 0. ) cconz = scal(l,k)
            mflux = qflux * cconz
	    if( qflux .ne. 0 ) then
	      write(46,1000) k,l,mflux,qflux,cconz,scal(l,k)
	    end if
          end do
        end do

	return
 1000	format(2i10,4f10.4)
	end

!**********************************************************************
!**********************************************************************
!*******************************************************************

	subroutine mult_scal_bc(r3v,value)

! multiplies array for scalar boundary condition with value

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	double precision r3v(nlvdi,nkn)
	double precision value

	include 'mkonst.h'

	integer k,l

	do k=1,nkn
	  do l=1,nlv
	    if( r3v(l,k) .ne. flag ) r3v(l,k) = r3v(l,k) * value
	  end do
	end do

	end


!**********************************************************************

	subroutine dist_horizontal(nlvddi,r3v,n,value)

	implicit none

	integer nlvddi
	double precision r3v(nlvddi,n)
	integer n
	double precision value

	integer k

	do k=1,n
	  r3v(1,k) = value
	end do
	  
	end

!**********************************************************************

        subroutine aver_horizontal(nlvddi,r3v,n,value)

        implicit none

        integer nlvddi
        double precision r3v(nlvddi,n)
        integer n
        double precision value

        integer k

	value = 0.
        do k=1,n
          value = value + r3v(1,k)
        end do
	value = value / n

	end

!**********************************************************************

	subroutine print_scal_bc(r3v)

! prints non-flag entries of array for scalar boundary condition

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	double precision r3v(nlvdi,nkn)

	include 'mkonst.h'

	integer k,l
	double precision value

	do k=1,nkn
	  do l=1,nlv
	    value = r3v(l,k)
	    if( value .ne. flag ) then
		write(6,*) 'print_scal_bc: ',k,l,value
	    end if
	  end do
	end do

	end

!**********************************************************************

	subroutine get_bflux(k,flux)

! returns boundary flux of node k

	use bnd_dynamic

	implicit none

	integer k	!node
	double precision flux	!flux in node k (return)

	flux = rqv(k)

	end

!**********************************************************************

	subroutine level_flux(it,levflx,kn,rw)

! compute discharge from water level

	use hydro_admin

	implicit none

	integer it		!type of function
	integer levflx		!type of function
	integer kn		!node number
	double precision rw			!discharge computed

	double precision z,a,b,c,z0

! use fit 1 for compatibility to old simulations
! use fit 3 for best results without using piecewise fit

	if( levflx .eq. 0 ) then
	  !nothing changed -> return same rw
	else if( levflx .eq. 1 ) then
	  z = znv(kn)
	  call z_smooth(z)

	  !a = 1648		!no outliers
	  !c = -2684		!no outliers

!-------------- old fit ----------------------------	1
	  z0 = 5.		!reference level
	  a = 1808
	  c = -2875
	  rw = a*log(z0+z) + c
!---------------------------------------------------

!-------------- fit based on mass balance -------------  2
	  !a = 211.521555200211
	  !b = 0.855968510497031     
	  !rw = a * z**b
!---------------------------------------------------

!-------------- new fit based on mass balance ----------  3
	  !a = 223.529239458563		!new calibration
	  !b = 0.862816507081288
	  !rw = a * z**b
!---------------------------------------------------

!-------------- best fit based on mass balance -----------  4
          !call flow_out_piece_new(z,rw)	!this is best (if it works)
!---------------------------------------------------

	  rw = -rw
	  write(134,*) it,z,rw
	else
	  write(6,*) 'levflx = ',levflx
	  stop 'error stop level_flux: levflx'
	end if

	end

!**********************************************************************

	subroutine z_smooth(z)

! smooths z values

	implicit none

	integer ndim
	parameter (ndim=86400/100)

	double precision z

	integer i

	integer ipz
	double precision za(0:ndim)
	save ipz,za
	data ipz / 0 /

	if( ipz .le. 0 ) then	!initialize
	  do i=0,ndim
	    za(i) = z
	  end do
	  ipz = 1
	end if

	za(0) = za(0) + (z-za(ipz))/ndim
	za(ipz) = z
	ipz = mod(ipz,ndim) + 1

	z = za(0)

	end

!**********************************************************************

        subroutine flow_out_piece_new(z,rout)

        implicit none

        double precision z,rout
        double precision rw,a,b

        if( z .le. 1. ) then
          a = 100.594731852550
          b = 78.7626400203164
        else if( z .le. 2. ) then
          a = -123.192333413729
          b = 285.797516554523
        else if( z .le. 3. ) then
          a = 21.2143155529982
          b = 207.240814424064
        else
          a = 468.988395301776
          b = 102.075378828782
        end if

        rw = a + z*b
        !rout = -rw
        rout = rw

        end

!**********************************************************************

!----------------------------------------------------------------------
        end module bnd_routines
!----------------------------------------------------------------------
