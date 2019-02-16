
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c subroutines for computing lagrangian trajectories
c
c contents :
c
c subroutine lagrange
c
c       read str
c       subroutine set_input
c       subroutine nbody
c       subroutine setup_fluxes         initializes flux2d
c
c       subroutine drogue(it)           compute trajectories
c
c       subroutine set_output
c
c subroutine back_trace
c
c	subroutine setbar		initial t,x,y floats
c       subroutine setup_fluxes         initializes flux2d
c       subroutine drogue(it)           compute trajectories
c       subroutine lagr_vel
c
c subroutine drogue(it)			compute trajectories
c
c       do                           	loop on floats
c         subroutine dtime           	decay -> out of loop
c         subroutine track_body
c             subroutine track_orig
c             subroutine track_line
c       end do
c
c subroutine setup_fluxes        	initializes flux2d
c
c       subroutine getaz
c       do
c         function flxtype
c	  subroutine get_elem_linkp
c         subroutine mk_rflux           flux through volume k
c         subroutine mk_tflux           flux through vertexes
c         subroutine setup_fx           set up flux2d(3,nel)
c       end do
c       subroutine setup_vl
c
c revision log :
c
c 00.00.2003    aac     routines written from scratch
c 29.04.2005    ggu     routines cleaned
c 01.10.2005    aac     diffusion routines written
c 07.11.2005    ggu     diffusion integrated
c 20.12.2005    aac&ggu bug in track_body corrected
c 19.06.2006    aac     bugs in lagrange.f corrected
c 20.06.2006    aac     2D lagrangian code stable 
c 22.06.2006    aac     lagrangian custom routine introduced
c 29.11.2006    ggu     lots of small changes, integrated into main model
c 06.06.2007    ggu     use of lcust commented (?)
c 10.11.2007	ggu	new routine ggrand, new call to lagr_release
c 23.04.2008	ggu	drogue() parallelized
c 29.04.2008	ggu	bug fix for parallel version
c 24.06.2008	ggu	new z var, new initialization
c 10.07.2008	aac	final stable version with new diffusion
c 29.01.2009	aac	changes in write to file
c 05.02.2009    ggu     re-arranged whole lagrangian module
c 15.02.2009    ggu     call to track_body has changed -> pass time to advect
c 11.09.2009    ggu     little bug fix for output and release of particles
c 19.10.2011    ggu     fx renamed to flux2d
c 16.12.2011    ggu     new file .lgi, compress_particles()
c 23.01.2012    ggu     various changes in call to track_body (id, etc..)
c 24.01.2012    ggu     adapted for parallel OMP
c 28.08.2012    ggu     change logic for release, time frame for release
c 22.10.2012    ggu     call connectivity also after diffusion
c 22.10.2012    ggu     limit release to itranf/end
c 28.03.2014    ggu     code cleaned - connectivity
c 10.04.2014    ggu     new code for lagr_count
c 23.04.2015    ggu     internal coordinates implemented (blgrxi)
c 26.05.2017    ccf     integrate vertical diffusion
c 10.07.2018    ccf     new data structures
c 23.08.2018    ccf     including particle beaching
c
c****************************************************************            

	subroutine lagrange

c lagranian main routine

	use mod_lagrange
	use basin, only : nkn,nel,ngr,mbw
	use levels
        use lgr_sedim_module
	use shyfile

	implicit none

	include 'femtime.h'

        double precision, save		:: dtlanf,dtlend
	double precision, save		:: ddtl,dtranf,dtrend,dtrnext
        double precision		:: dtime
	real, save			:: ldecay

        integer				:: ifemop
        real 				:: getpar

	logical				:: brelease
        character*20 			:: aline
        double precision, save          :: da_lgr(4) = 0
	integer, save			:: iu
        integer, parameter		:: nvar = 1
        logical, parameter 		:: b3d = .false.

	integer				:: id
        logical	 			:: has_output_d
        logical                         :: next_output_d
        external 			:: get_timeline

	integer, save 			:: icall = 0
        
        if( icall == -1 ) return
        
c---------------------------------------------------------------
c set some parameters
c---------------------------------------------------------------

c pps and ppv have to be set in STR file as lgrpps (section BOUND)
c lgrpps > 0 => pps
c lgrpps < 0 => ppv
c
c the following parameters have to be set in section lagrg
c
c tdecay	ldecay
c boilsim	ioil
c blarvae	ilarv
c bconnect	iconnect

c---------------------------------------------------------------
c initialization
c---------------------------------------------------------------

        if( icall .eq. 0 ) then
          ilagr = nint(getpar('ilagr'))
          if( ilagr .eq. 0 ) icall = -1
          if( icall .eq. -1 ) return
	  icall = 1

          call get_first_dtime(dtanf)
          call get_last_dtime(dtend)
	  call convert_date_d('itlanf',dtlanf)   !start of lagrangian sim
	  call convert_date_d('itlend',dtlend)   !end of lagrangian sim
	  call convert_time_d('idtl',ddtl)       !frequency of release
	  call convert_date_d('itranf',dtranf)   !time of initial cont. release
	  call convert_date_d('itrend',dtrend)   !time of final continuous release

          ldecay = getpar('ldecay')              !decay time for particles
          lbeach = getpar('lbeach')              !beaching factor for particles
          bbeach = lbeach > 0

          boilsim = nint(getpar('ioil')).gt.0    !activate oil module if true
          blarvae = nint(getpar('ilarv')).gt.0   !activate larvae module if true
	  bsedim  = nint(getpar('ised')).gt.0    !activate sediment module if true

          nbdymax = nint(getpar('nbdymax'))
	  if( nbdymax < 0 ) then
	    write(6,*) 'parameter nbdymax is not set'
	    stop 'error stop lagrange: nbdymax'
	  end if
	  write(6,*) 
	  write(6,*)'---------------------------------------------' 
	  write(6,*)'Initialization of the LAGRANGIAN module'
	  write(6,*)'nbdymax = ',nbdymax

	  call mod_lagrange_init(nel,nlv)
	  call mod_lagrange_handle_alloc(0)

	  if ( bsedim ) call lgr_sedim_init

	  call lagr_init_common

c	  if( boilsim ) call init_diff_oil

c         ------------------------------------------------------
c	  lagrangian module
c         ------------------------------------------------------

	  if( dtlanf == -1.d0 ) dtlanf = dtanf
	  if( dtlend == -1.d0 ) dtlend = dtend
	  if( dtlanf < dtanf ) dtlanf = dtanf
	  if( dtlend > dtend ) dtlend = dtend

c         ------------------------------------------------------
c	  new release
c         ------------------------------------------------------

	  if( dtranf == -1.d0 ) dtranf = dtlanf
	  if( dtrend == -1.d0 ) dtrend = dtlend
	  if( dtranf < dtlanf ) dtranf = dtlanf
	  if( dtrend > dtlend ) dtrend = dtlend
	  dtrnext = dtranf
	  if( ddtl == 0.d0 ) ddtl = dtlend - dtlanf + 1	!release once at start
	  if( ddtl < 0.d0 ) dtrnext = dtend + 1		!never release

c         ------------------------------------------------------
c	  initialize particle distribtuion from lgr file 
c         ------------------------------------------------------

	  call lgr_input_shell

c         ------------------------------------------------------
c	  open output file and write ncust
c         ------------------------------------------------------

          call init_output_d('itmlgr','idtlgr',da_lgr)
          if( has_output_d(da_lgr) ) then
            call shyfem_init_lgr_file('lgr',nvar,b3d,id)
            da_lgr(4) = id
            call shy_get_iunit(id,iu)
            write(iu) ncust
          end if
	end if
         
c---------------------------------------------------------------
c run lagrangian 
c---------------------------------------------------------------

        call get_act_dtime(dtime)

        if( dtime < dtlanf .or. dtime > dtlend ) return      

	bback = .false.		!do not do backtracking

c---------------------------------------------------------------
c new release of particles (more release event, homogeneous or lines)
c---------------------------------------------------------------
	
	if( dtime >= dtrnext .and. dtime <= dtrend ) then
          call get_timeline(dtime,aline)
	  write(6,*) 'release of particles for lagrangian model'
	  call lgr_init_shell
	  dtrnext = dtrnext + ddtl
	  if( dtrnext == dtend ) dtrnext = dtend + 1
	  write(6,*) 'new particles released: ',nbdy,'at time: ',aline
        end if           

c---------------------------------------------------------------
c new release of particles (on points along trajectory)
c---------------------------------------------------------------

        call lgr_init_traj(dtime)

c---------------------------------------------------------------
c one time step of particle tracking
c---------------------------------------------------------------

        call lagr_setup_timestep
	
c	if( boilsim ) call set_diff_oil

c---------------------------------------------------------------
c continuous release from boundary or points
c lgrpps or lgrppv defined in boundary section
c---------------------------------------------------------------

	brelease = dtime >= dtranf .and. dtime <= dtrend

	if( brelease ) then
	  call lagr_continuous_release_shell
	end if
	
        call lagr_connect_continuous_points(brelease)
	call lagr_count_init 

c---------------------------------------------------------------
c Compute vertical diffusivity in element and random walk time step
c---------------------------------------------------------------

 	call lag_vdiff_ele

c---------------------------------------------------------------
c transport of particles 
c---------------------------------------------------------------

 	call drogue(dtime)

c---------------------------------------------------------------
c sediment module
c---------------------------------------------------------------

	if( bsedim ) call lgr_sediment

c---------------------------------------------------------------
c larval module
c---------------------------------------------------------------

	if( blarvae ) call lgr_larvae

c---------------------------------------------------------------
c decay
c---------------------------------------------------------------

	call lagrange_decay(ldecay)

c---------------------------------------------------------------
c output : connectivity matrix or trajectiories
c---------------------------------------------------------------

        if( next_output_d(da_lgr) ) then
	  call lgr_output(iu,dtime)
        end if

        if( bcount .and. dtime == dtlend )then
          !call lagr_count_out_eos(it)
          call lagr_count_out(dtime,dtlend)
        end if

c---------------------------------------------------------------
c compress to save space: only in contiunous release mode! 
c CCF STILL TO BE CHECKED
c---------------------------------------------------------------

	if( bcompress ) then 
	  call compress_particles	!only after output of particles
	endif 

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************

	subroutine lagr_init_common

c initializes common block

	use mod_lagrange
	use levels

	implicit none

	real getpar
	integer ifemop

        artype=getpar('artype')	!store special type in common block
        rwhpar=getpar('rwhpar')	!lagrangian diffusion
        stkpar=getpar('stkpar')	!stokes drift parameter
        dripar=getpar('dripar')	!drifter parameter

	!lunit = ifemop('.lgi','form','new') !unit for lagrangian info
	!if( lunit .le. 0 ) then
	!  write(6,*) 'lunit = ',lunit
	!  stop 'error stop lagrange: cannot open info file'
	!end if

	nbdy = 0 		!number of particles to insert
	idbdy = 0		!id body unique 
	tdecay = 0		!not used anymore

	blgrdebug = .false.

c ilagr = 1	surface lagrangian
c ilagr = 2	2d lagrangian (without vertical adv and diff)
c ilagr = 3	3d lagrangian

	blgrsurf = ilagr == 1
	blgr2d = ilagr == 2

c no vertical diffusion for surface lagrangian
	if ( blgrsurf .or. blgr2d ) bvdiff = .false.

c vertical distribution of particles
c n = abs(ipvert)
c ipvert == 0    release one particle in surface layer
c ipvert > 0     release n particles regularly
c ipvert < 0     release n particles randomly

        ipvert = nint(getpar('ipvert'))

c lintop and linbot= top and bottom layer between perform the release

        lintop =  getpar('lintop') 
        linbot =  getpar('linbot')   
	if ( linbot == -1 ) linbot = nlv

	end

c**********************************************************************

c*******************************************************************
! set properties of particles
!   - settling velocity [m/s]
!   - particle type 
!   - curstom properties (nc)

        subroutine lgr_set_properties(bsedim,blarvae,boilsim,pt,ps,pc,
     &                                nc)

        use lgr_sedim_module, only : lgr_set_sedim

        implicit none

        logical, intent(in)           :: bsedim  !true for sediment lagrangian module
        logical, intent(in)           :: blarvae !true for larvae module
        logical, intent(in)           :: boilsim !true for oil module
        integer, intent(inout)        :: pt      !particle type
        double precision, intent(out) :: ps      !settling velocity [m/s]
        real, intent(out)             :: pc      !custom property
        integer, intent(out)          :: nc      !number of custom properties

        ps = 0.
        pc = 0.
        nc = 1

        if ( bsedim ) call lgr_set_sedim(pt,ps,pc,nc)
        !if ( blarvae ) call lgr_set_larvae(pt,ps,pc)       !pc=length 
        !if ( boilsim ) call lgr_set_boilsim(pt,ps,pc)  !TODO ccf

        end subroutine lgr_set_properties

c*******************************************************************

	subroutine drogue(dtime)

	use mod_lagrange

	implicit none

	double precision dtime	
	integer nf,i,ii,n
	integer chunk
	real dt

	integer ndim
	parameter(ndim=100)
	integer ic(0:ndim)

        integer iuinfo
        save iuinfo
        data iuinfo / 0 /

        double precision tempo
        double precision openmp_get_wtime

        if( iuinfo .eq. 0 ) then
          call getinfo(iuinfo)  !unit number of info file
        end if

	chunk = 100
        nf=0
	call get_timestep(dt)		!time to advect

	call openmp_get_max_threads(n)
	if( n .gt. ndim ) stop 'error stop drogue: too many processes'

	do ii=0,n
	  ic(ii) = 0
	end do

        tempo = openmp_get_wtime()

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ii)
!$OMP DO SCHEDULE(DYNAMIC,chunk)

	do i=1,nbdy
	  call openmp_get_thread_num(ii)
	  ic(ii) = ic(ii) + 1
	  call track_single(i,dt,dtime)
	end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

        tempo = openmp_get_wtime() - tempo
!        write(88,*) 'tempo = ',tempo

c	write(lunit,*) 'lagrangian: (tot,out,in) ',nbdy,nf,nbdy-nf
c	write(lunit,'(a,i10,12i5)') 'parallel: ',nbdy,(ic(ii),ii=0,n-1)

	write(iuinfo,*) 'lagrange_nbdy: ',nbdy
	!write(6,*) 'lagrange_nbdy: ',nbdy

	end

c**********************************************************************

	subroutine track_single(i,dt,dtime)

! advection of particles

	use mod_lagrange
        use mod_geom_dynamic, only : iwegv
        use levels, only : nlv

	implicit none
	
	double precision dtime
	integer i,id,ie,nf,lb,ii
	real x,y,z
	double precision xx,yy
	double precision sv
	double precision xi(3)
	integer ty
	real dt,ttime,tmax
        integer lmax                    !maximum layers in element (return)
        real hl(nlv)
        real htot,htotz

!       call lagr_func(i) !lcust in str varying the variable typ(i) to check
!	call lagr_surv(i)

c---------------------------------------------------------------
c Get particle properties
c---------------------------------------------------------------

        ie = lgr_ar(i)%actual%ie
	do ii=1,3
	  xi(ii) = lgr_ar(i)%actual%xi(ii)
	end do
	call xi2xy(abs(ie),xx,yy,xi)
	x  = xx
	y  = yy
	z  = lgr_ar(i)%actual%z
	lb = lgr_ar(i)%actual%l
	sv = lgr_ar(i)%sinking
	id = lgr_ar(i)%id 
	ty = lgr_ar(i)%type

c---------------------------------------------------------------
c Return if element is dry or negative first layer (offline problem)
c---------------------------------------------------------------

        if( ie <= 0 ) return            !return if particle out of domain
        if( iwegv(ie) /= 0 ) return	!return if particle on dry element
	lmax = nlv
        call lagr_layer_thickness(ie,lmax,hl,htot,htotz)
	if ( hl(1) .lt. 0. ) return

c---------------------------------------------------------------
c Get travel time for particle
c---------------------------------------------------------------

	tmax = dtime - lgr_ar(i)%actual%time
	if( tmax .lt. 0. ) stop 'error stop drogue: internal error'
	ttime = min(tmax,dt) 		!residual time for particle

c---------------------------------------------------------------
c Compute advection and diffusion
c---------------------------------------------------------------

	if( lb > 0 ) then
	  if( blgrxi ) then		!use internal coordinates
            call track_body_xi(i,id,ty,x,y,z,sv,xi,ie,lb,ttime) 
	    call diff_body(i,id,ty,x,y,z,xi,ie,lb,ttime)
	    call lag_stk(i,id,ty,x,y,xi,ie,lb,ttime)
	    call lag_beach(i,id,xi,ie,ty)
	  else
            call track_body(i,id,x,y,z,lb,ie,ttime) 
	  end if
	end if

c---------------------------------------------------------------
c Assign new coordinatates to particle
c---------------------------------------------------------------

        lgr_ar(i)%actual%ie = ie
	do ii=1,3
	  lgr_ar(i)%actual%xi(ii) = xi(ii)
	end do
	lgr_ar(i)%actual%z  = z
	lgr_ar(i)%actual%l  = lb
	lgr_ar(i)%type      = ty
	if ( ie > 0 ) lgr_ar(i)%actual%time = dtime

	end	

c**********************************************************************

        subroutine track_body_xi(i,id,ty,x,y,z,sv,xi,iel,lb,time)

c tracks one particle - uses internal coordinates

	use mod_lagrange
	use basin

	implicit none

	integer i		!particle number
	integer id		!particle id
        integer ty		!particle type
	real x			!x-coordinate
	real y			!y-coordinate
	real z 			!relative vertical position
	double precision sv 	!settling velocity
	double precision xi(3)	!internal coordinates
	integer iel		!element number
	integer lb 		!layer 
	real time		!time to advect
        
	integer n
	integer ie_from,ie_to
	integer iperc
	real ttime,torig
	real perc
	double precision xx,yy,zz

	integer, save :: nk = 0
	integer, save :: nl = 0
	integer, save :: nu = 0

        if(iel <= 0 .or. ty < 0) return	!particle out of domain or beached

c---------------------------------------------------------------
c initialize
c---------------------------------------------------------------

        n = 100		!maximum loop count
	zz = z
	ttime = time*dripar	!accout for drifter inertia
	
c---------------------------------------------------------------
c track particle
c---------------------------------------------------------------

	do while ( ttime > 0 .and. n > 0 )
	  torig = ttime 		!time do advect
	  ie_from = iel 		!start element 
          call track_xi(id,iel,lb,sv,xi,zz,ttime) !advection 
	  ie_to = iel 			!new element if advection changed element

	  if ( bconnect ) then		!connectivity
	    !write(6,*) i,ie_to,ie_from,torig-ttime,'-',id
	    call lagr_connect_count(i,ie_to,ie_from,torig-ttime)
	    call lagr_count(i,ie_to,ie_from,ttime)
	  endif

	  if( iel < 1 ) exit
	  if( lb < 1 ) exit
	  n = n - 1			!decrement loop counter
	end do

	call xi2xy(abs(iel),xx,yy,xi)
	x = xx
	y = yy
	z = zz

	if( id == 0 ) then
	  write(6,*) 'lgrggu: ',iel,lb,ttime
	  write(6,*) 'lgrggu: ',x,y,z
	  write(6,*) 'lgrggu: ',xi
	end if

c---------------------------------------------------------------
c special treatment and finish up
c---------------------------------------------------------------

	if( ttime > 0. ) then 		!not finished advecting
	  if( n == 0 ) then
	    nk = nk + 1
	    perc = (100.*nk)/idbdy
	    write(6,1000) 'warning particle adv',id,iel,n,ttime,perc
	    !iel = -iel
	  else if( iel < 1 ) then
	    nl = nl + 1
	    perc = (100.*nl)/idbdy
	    write(6,1000) 'loosing particle adv',id,iel,n,ttime,perc
	  else
	    nu = nu + 1
	    perc = (100.*nu)/idbdy
	    write(6,1000) 'unknown error adv',id,iel,n,ttime,perc
	  end if
	end if

        if( .not. bback ) then
          if( iel.gt.0 .and. iarv(iel).eq.artype ) iel = -iel
	end if

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

 1000	format(a,3i10,2f10.2)
	end

!**********************************************************************
! Set beaching of particle when it is on a meterial boundary. 
! lbeach in the range 0-1. If equal 0, no beaching. 
! When the particle reaches the shore it remains inside the domain,
! with type=-type and could not move 
! CCF SHOULD WE ALLOW BEACHING ONLY ON SURFACE?

        subroutine lag_beach(i,id,xi,iel,ty)

        use mod_lagrange

        implicit none

        integer, intent(in)	     :: i       !particle number
        integer, intent(in)	     :: id      !particle id
        double precision, intent(in) :: xi(3)   !internal coordinates
        integer, intent(in)	     :: iel     !element number
        integer, intent(inout)	     :: ty	!particle type

	real 		:: r, perc
        logical 	:: track_xi_on_material_boundary

	if ( .not. bbeach .or. ty < 0 ) return

        if( track_xi_on_material_boundary(iel,xi) ) then
          call random_number(r)
	  if ( lbeach > r ) ty = -ty
	end if

	end

c**********************************************************************

        subroutine track_body(i,id,x,y,z,lb,iel,ttime)

c tracks one particle
c
c uses two routines:
c
c TRACK_ORIG if the particle is inside an element (first call)
c
c TRACK_LINE if the particle is on one side (normal situation)

	use mod_lagrange
	use basin

	implicit none

	integer i		!particle number
	integer id		!particle id
	integer iel		!element number
	integer ielem		!element number to check
	real x			!x-coordinate
	real y			!y-coordinate
	real ttime		!time to advect
	real z 			!relative vertical position
	integer lb 		!layer 
	
        
	integer nl
	integer ltbdy
	integer ieold,ieorig
	real torig
	real xn,yn,zn
	integer ly

        if(iel.le.0) return	!particle out of domain

c---------------------------------------------------------------
c initialize
c---------------------------------------------------------------

        nl = 100		!maximum loop count
	ltbdy = 0

	xn = x
	yn = y
	zn = z
	ly = lb

c---------------------------------------------------------------
c track particle
c---------------------------------------------------------------

        torig = ttime
        ieold = iel     !element the particle is in or leaving
        ieorig = iel    !original element the particle was in
        call track_orig(ttime,id,iel,xn,yn,zn,ly,ltbdy)

        do while ( ttime.gt.0. .and. iel.gt.0 .and. nl.gt.0 )
          torig = ttime
          ieold = iel
          call track_line(ttime,id,iel,xn,yn,zn,ly,ltbdy)
          nl = nl - 1
          ieorig = ieold
        end do

	if( .false. ) then

	torig = ttime
	ieold = iel	!element the particle is in or leaving
	ieorig = iel	!original element the particle was in
        call track_orig(ttime,id,iel,xn,yn,zn,ly,ltbdy)

	do while ( ttime.gt.0 .and. iel.eq.ieold )
          call track_orig(ttime,id,iel,xn,yn,zn,ly,ltbdy)
	end do
	
	if( bconnect ) then
	  call lagr_connect_count(i,ieold,ieorig,torig-ttime)
	  call lagr_count(i,ieold,ieorig,torig-ttime)
	end if 

	do while ( ttime.gt.0. .and. iel.gt.0 .and. nl.gt.0 )
	  torig = ttime
	  ieold = iel
          call track_line(ttime,id,iel,xn,yn,zn,ly,ltbdy)
	  do while ( ttime.gt.0 .and. iel.eq.ieold )
            call track_orig(ttime,id,iel,xn,yn,zn,ly,ltbdy)
          end do
          nl = nl - 1
	  if(  bconnect  )  then
	    call lagr_connect_count(i,ieold,ieorig,torig-ttime)
	    call lagr_count(i,ieold,ieorig,torig-ttime)
	    ieorig = ieold
	  endif 
        end do

	end if

c---------------------------------------------------------------
c error condition (infinite loop)
c---------------------------------------------------------------

        if( nl .eq. 0 ) then
          print*, 'inifinite loop in track_line'
          print*, i,iel,xn,yn,ttime
          iel = -iel
	  ttime = 0.
        end if

c---------------------------------------------------------------
c diffusion
c---------------------------------------------------------------

       if( .not. bback .and. iel .gt. 0 .and. rwhpar .gt. 0 ) then
	  ! the time spent in elemets due to diffusion is not considered
	  ttime = 0.
	  ieorig = iel
          call lag_diff(iel,id,xn,yn)	
	  if( bconnect ) then
	    call lagr_connect_count(i,iel,ieorig,ttime)
	  endif
        end if     

c---------------------------------------------------------------
c special treatment and finish up
c---------------------------------------------------------------

        if( .not. bback ) then
          if( iel.gt.0 .and. iarv(iel).eq.artype ) iel = -iel
	end if

	x = xn
	y = yn
	z = zn
	lb = ly

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine lagr_count_init

	use mod_lagrange
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ie

	do ie=1,nel
	  i_count(ie) = 0
	  t_count(ie) = 0.
	end do

	end

c**********************************************************************

	subroutine lagr_count(i,ie_to,ie_from,time)

	use mod_lagrange

	implicit none

	integer i
	integer ie_to,ie_from
	real time

	integer ic

	if( ie_to .le. 0 ) return
	if( ie_from .le. 0 ) stop 'error stop lagr_count: (1)'

	ic = 1
	if( ie_to == ie_from ) ic = 0

	i_count(ie_to) = i_count(ie_to) + ic
	t_count(ie_from) = t_count(ie_from) + time

	end

c**********************************************************************

	subroutine lagr_count_out(dtime,dtlend)

c write a map of total number of particles passed in each element or total time spent
c aux index = time average spent in each element 
c TODO -> normalization

	use mod_lagrange
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	double precision dtime,dtlend

	integer ie,iu
	real aux

	character*80 file

	iu = 237
	file = 'lagr_count_out.txt'
	open(iu,file=file,status='unknown',form='formatted')
	write(iu,*) dtime,nel
	do ie=1,nel
	  if(i_count(ie).gt.0)then
	    aux = (t_count(ie)/86400.) / i_count(ie) 
	    write(iu,*) ie,i_count(ie),(t_count(ie)/86400.),aux
	  else 
	    write(iu,*) ie,0.,0.
	  end if 
	end do

        if(dtime >= dtlend) close(iu)

	end

c**********************************************************************

	subroutine lagr_count_out_eos(it)

c write an eos of total number of particle passed in each element or total time spent
c aux index = time average spent in each element 
c TODO -> normalization

	use mod_lagrange
	use basin, only : nkn,nel,ngr,mbw
	use levels
	use mod_depth

        implicit none

        integer ie,iu
        character*80 file,title
        save file,title

        real aux(nel)
        real aux_t(nel)
        real aux_r(nel)

        integer nvers,nlvdimi
        integer nvar,ivar,iunit,it,ierr

        integer ifileo
        integer, save :: icall = 0

        nvers = 3
        nlv = 1
        nvar = 1

        do ie=1,nel
          aux(ie) = i_count(ie)
          aux_t(ie) = t_count(ie)/86400
          if (aux(ie).ne.0)then
          aux_r(ie) = aux_t(ie)/aux(ie)
          else
          aux_r(ie) = 0
          end if
        end do

        if (icall.eq.0) then
        file = 'particle_traj.eos'
        iunit = ifileo(0,file,'unform','new')

        title ='particle in element '
        call wheos(iunit,nvers
     +             ,nkn,nel,nlv,nvar
     +             ,ilhv,hlv,hev
     +             ,title
     +             )
        endif

        icall = 1
        ivar = 901
        nlvdimi=1
        call wreos(iunit,it,ivar,nlvdimi,ilhv,aux,ierr)
        if( ierr .ne. 0 ) stop 'error stop: wreos'
        call wreos(iunit,it,ivar,nlvdimi,ilhv,aux_t,ierr)
        if( ierr .ne. 0 ) stop 'error stop: wreos'
        call wreos(iunit,it,ivar,nlvdimi,ilhv,aux_r,ierr)
        if( ierr .ne. 0 ) stop 'error stop: wreos'

        do ie=1,nel
         write(99,*) ie,i_count(ie),t_count(ie)
        enddo
        end

c**********************************************************************

