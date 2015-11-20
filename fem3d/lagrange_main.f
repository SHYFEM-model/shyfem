c
c $Id: lagrange_main.f,v 1.5 2010-03-11 15:36:38 georg Exp $
c
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
c
c****************************************************************            

	subroutine lagrange

c lagranian main routine

	use mod_lagrange
	use basin, only : nkn,nel,ngr,mbw
	use levels
        use lgr_sedim_module

	implicit none

        include 'param.h'

	include 'femtime.h'

	logical brelease

        integer itlanf,itlend
	integer idtl,itranf,itrend,itrnext
	integer itmlgr,idtlgr,itmnext
        integer iunit,uunit
	real ldecay

        integer ifemop
        real getpar

        save itlanf,itlend
	save idtl,itranf,itrend,itrnext
	save itmlgr,idtlgr,itmnext
        save iunit,uunit
	save ldecay

	integer icall
	save icall
        data icall / 0 /
        
        if( icall .eq. -1 ) return
        
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

	  call convert_date('itmlgr',itmlgr)   !startime write output to file
	  call convert_time('idtlgr',idtlgr)   !frequency output to file
	  call convert_date('itlanf',itlanf)   !start of lagrangian sim
	  call convert_date('itlend',itlend)   !end of lagrangian sim
	  call convert_time('idtl',idtl)       !frequency of release
	  call convert_date('itranf',itranf)   !time of initial cont. release
	  call convert_date('itrend',itrend)   !time of final continuous release

          ldecay =  getpar('ldecay')           !decay time for particles

          boilsim = nint(getpar('ioil')).gt.0  !activate oil module if true
          blarvae = nint(getpar('ilarv')).gt.0 !activate larvae module if true
	  bsedim  = nint(getpar('ised')).gt.0  !activate sediment module if true

          nbdymax = nint(getpar('nbdymax'))
	  if( nbdymax < 0 ) then
	    write(6,*) 'parameter nbdymax is not set'
	    stop 'error stop lagrange: nbdymax'
	  end if
	  write(6,*) 'nbdymax = ',nbdymax

	  call mod_lagrange_init(nel,nlv)
	  call mod_lagrange_handle_alloc(0)

	  if ( bsedim ) call lgr_sedim_init

	  call lagr_init_common

c	  if( boilsim ) call init_diff_oil

c         ------------------------------------------------------
c	  lagrangian module
c         ------------------------------------------------------

	  if( itlanf .eq. -1 ) itlanf = itanf
	  if( itlend .eq. -1 ) itlend = itend
	  if( itlanf .lt. itanf ) itlanf = itanf
	  if( itlend .gt. itend ) itlend = itend

c         ------------------------------------------------------
c	  output
c         ------------------------------------------------------

	  if( itmlgr .eq. -1 ) itmlgr = itlanf
	  if( itmlgr .lt. itlanf ) itmlgr = itlanf
	  itmnext = itmlgr
	  if( idtlgr .le. 0 ) itmnext = itlend + 1

c         ------------------------------------------------------
c	  new release
c         ------------------------------------------------------

	  if( itranf .eq. -1 ) itranf = itlanf
	  if( itrend .eq. -1 ) itrend = itlend
	  if( itranf .lt. itlanf ) itranf = itlanf
	  if( itrend .gt. itlend ) itrend = itlend
	  itrnext = itranf
	  if( idtl .eq. 0 ) idtl = itlend - itlanf + 1	!release once at start
	  if( idtl .lt. 0 ) itrnext = itend + 1		!never release

c         ------------------------------------------------------
c	  open files
c         ------------------------------------------------------

          if( artype .ne. -1 ) then
          uunit=ifemop('.trn','form','new')
          end if
          iunit=ifemop('.lgr','unform','new')

	end if
         
c---------------------------------------------------------------
c run lagrangian ?
c---------------------------------------------------------------

        if( it .lt. itlanf .or. it .gt. itlend ) return      

	bback = .false.		!do not do backtracking

c -------------------------------------------
c check for initialization from file.lgr 
c---------------------------------------------

c---------------------------------------------------------------
c new release of particles (more release event, homogeneous or lines)
c---------------------------------------------------------------
	
	if( it .ge. itrnext .and. it .le. itrend ) then
	  write(6,*) 'release of particles for lagrangian model'
	  call lgr_init_shell
	  itrnext = itrnext + idtl
	  if( itrnext .eq. itend ) itrnext = itend + 1
	  write(6,*) 'new particles released: ',nbdy,it
        end if           

c---------------------------------------------------------------
c one time step of particle tracking
c---------------------------------------------------------------

        call lagr_setup_timestep
	
c	if( boilsim ) call set_diff_oil

c---------------------------------------------------------------
c continuous release from boundary or points
c lgrpps or lgrppv defined in boundary section
c---------------------------------------------------------------

	brelease = it .ge. itranf .and. it .le. itrend

	if( brelease ) then
	  call lagr_continuous_release_shell
	end if
	
        call lagr_connect_continuous_points(brelease)

c---------------------------------------------------------------
c transport of particles 
c---------------------------------------------------------------

 	call drogue

c---------------------------------------------------------------
c sediment module
c---------------------------------------------------------------

	if( bsedim ) then
 	  call lgr_sediment(it)
	end if

c---------------------------------------------------------------
c larval module
c---------------------------------------------------------------

	if( blarvae ) then
 	  call lgr_larvae(it)
	end if

c---------------------------------------------------------------
c decay
c---------------------------------------------------------------

	call lagrange_decay(ldecay)

c---------------------------------------------------------------
c output
c---------------------------------------------------------------

        if( it .ge. itmnext ) then
	  call lgr_output(iunit,it)
 	  call lgr_output_concentrations
	  itmnext = itmnext + idtlgr
	end if

!        if(it.eq.itlend)then
!          call lagr_count_out_eos(it)
!          call lagr_count_out(it,itlend)
!        end if

c---------------------------------------------------------------
c compress to save space: only in contiunous release mode! 
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

	implicit none

	include 'param.h'

	real getpar
	integer ifemop

        artype=getpar('artype')	!store special type in common block
        rwhpar=getpar('rwhpar')	!lagrangian diffusion

	lunit = ifemop('.lgi','form','new') !unit for lagrangian info
	if( lunit .le. 0 ) then
	  write(6,*) 'lunit = ',lunit
	  stop 'error stop lagrange: cannot open info file'
	end if

	nbdy = 0 		!number of particles to insert
	idbdy = 0		!id body unique 
	tdecay = 0		!not used anymore

	blgrdebug = .false.

c ilagr = 1	surface lagrangian
c ilagr = 2	2d lagrangian (not implemented)
c ilagr = 3	3d lagrangian

	blgrsurf = ilagr == 1
	if( ilagr /= 1 .and. ilagr /= 3 ) then
	  write(6,*) 'ilagr = ',ilagr
	  stop 'error stop lagr_init_common: value for ilagr not allowed'
	end if

c vertical distribution of particles
c n = abs(ipvert)
c ipvert == 0    release one particle in surface layer
c ipvert > 0     release n particles regularly
c ipvert < 0     release n particles randomly

        ipvert = nint(getpar('ipvert'))

c lintop and linbot= top and bottom layer between perform the release
c       lintop =  getpar('lintop') 
        linbot =  getpar('linbot')   

	end

c**********************************************************************

	subroutine drogue

	use mod_lagrange

	implicit none
	
        include 'param.h'

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
	  call track_single(i,dt)
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

	subroutine track_single(i,dt)

c advection of particles

	use mod_lagrange

	implicit none
	
        include 'param.h'

	include 'femtime.h'

	integer i,id,ie,nf,lb,ii
	real x,y,z
	double precision sv
	double precision xi(3)
	real dt,ttime,tmax
                
c       call lagr_func(i) !lcust in str varying the variable typ(i) to check
c	call lagr_surv(i)

	x  = lgr_ar(i)%x
	y  = lgr_ar(i)%y
	z  = lgr_ar(i)%z
	lb = lgr_ar(i)%l
        ie = lgr_ar(i)%ie
	sv = lgr_ar(i)%sv
	id = lgr_ar(i)%id 

	do ii=1,3
	  xi(ii) = lgr_ar(i)%xi(ii)
	end do

	tmax = it - lgr_ar(i)%tin 
	if( tmax .lt. 0. ) stop 'error stop drogue: internal error'
	ttime = min(tmax,dt) !residual time for particle

	if( lb > 0 ) then
	  if( blgrxi ) then		!use internal coordinates
            call track_body_xi(i,id,x,y,z,sv,xi,ie,lb,ttime) 
	  else
            call track_body(i,id,x,y,z,lb,ie,ttime) 
	  end if
	end if

	lgr_ar(i)%x  = x
	lgr_ar(i)%y  = y
	lgr_ar(i)%z  = z
	lgr_ar(i)%l  = lb
        lgr_ar(i)%ie = ie
	do ii=1,3
	  lgr_ar(i)%xi(ii) = xi(ii)
	end do

	end	

c**********************************************************************

        subroutine track_body_xi(i,id,x,y,z,sv,xi,iel,lb,ttime)

c tracks one particle - uses internal coordinates

	use mod_lagrange
	use basin

	implicit none

	include 'param.h'
	
	integer i		!particle number
	integer id		!particle id
	real x			!x-coordinate
	real y			!y-coordinate
	real z 			!relative vertical position
	double precision sv 	!settling velocity
	double precision xi(3)	!internal coordinates
	integer iel		!element number
	integer lb 		!layer 
	real ttime		!time to advect
        
	integer n
	integer iendx,ieorig,ieold
	real torig
	double precision xx,yy,zz

        if(iel.le.0) return	!particle out of domain

c---------------------------------------------------------------
c initialize
c---------------------------------------------------------------

        n = 100		!maximum loop count
	zz = z

c---------------------------------------------------------------
c track particle
c---------------------------------------------------------------

	iendx = 0 	!flag for counting

	do while ( ttime.gt.0 .and. n > 0 )
	  torig = ttime 		!time do advect
	  ieorig =iel 			!start element 
          call track_xi(id,iel,lb,sv,xi,zz,ttime) !advection 
	  ieold=iel !if advection changed element= different start el

	  if ( bconnect ) then
	    call lagr_connect_count(i,ieold,ieorig,torig-ttime,iendx)
	    call lagr_count(i,ieold,ttime,iendx)
	    iendx = 1 
	  endif

	  if( iel < 1 ) exit
	  if( lb < 1 ) exit
	  n = n - 1
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
	    write(6,*) 'killing particle ',id,iel,n,ttime
	    iel = -iel
	  else if( iel < 1 ) then
	    write(6,*) 'loosing particle ',id,iel,n,ttime
	  else
	    write(6,*) 'unknown error ',id,iel,n,ttime
	  end if
	end if

        if( .not. bback ) then
          if( iel.gt.0 .and. iarv(iel).eq.artype ) iel = -iel
	end if

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

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

	include 'param.h'
	
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
        !call lagr_connect_count(i,ieold,ieorig,torig-ttime,0)
        !call lagr_count(i,ieold,torig-ttime,0)

        do while ( ttime.gt.0. .and. iel.gt.0 .and. nl.gt.0 )
          torig = ttime
          ieold = iel
          call track_line(ttime,id,iel,xn,yn,zn,ly,ltbdy)
          nl = nl - 1
          !call lagr_connect_count(i,ieold,ieorig,torig-ttime,1)
          !call lagr_count(i,ieold,torig-ttime,1)
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
	  call lagr_connect_count(i,ieold,ieorig,torig-ttime,0)
	  call lagr_count(i,ieold,torig-ttime,0)
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
	    call lagr_connect_count(i,ieold,ieorig,torig-ttime,1)
	    call lagr_count(i,ieold,torig-ttime,1)
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
	  call lagr_connect_count(i,iel,ieorig,ttime,1)
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

	include 'param.h'

	integer ie


	do ie=1,nel
	  i_count(ie) = 0
	  t_count(ie) = 0.
	end do

	end

c**********************************************************************

	subroutine lagr_count(i,ie,time,icc)

	use mod_lagrange

	implicit none

	include 'param.h'

	integer i
	integer ie
	real time
	integer icc

	if( ie .le. 0 ) return

	i_count(ie) = i_count(ie) + 1
	t_count(ie) = t_count(ie) + time

	end

c**********************************************************************

	subroutine lagr_count_out(it,itlend)

	use mod_lagrange
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	integer it,itlend

	integer ie,iu
	character*80 file

	iu = 237
	file = 'lagr_count_out.txt'
	open(iu,file=file,status='unknown',form='formatted')
	write(iu,*) it,nel
	do ie=1,nel
c	  i_count(ie) = 0
c	  t_count(ie) = 0
	  write(iu,*) ie,i_count(ie),(t_count(ie)/86400.)
	end do

        if(it.ge.itlend)then
         close(iu)
        endif

	end

c**********************************************************************

	subroutine lagr_count_out_eos(it)

	use mod_lagrange
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'

        integer ie,iu
        character*80 file,title
        save file,title

        integer ilhv(neldim)
        common /ilhv/ilhv
        real hev(neldim)
        common /hev/hev
        real hlv(1)
        common /hlv/hlv

        real aux(neldim)
        real aux_t(neldim)
        real aux_r(neldim)


        integer nvers,nlv,nlvdimi
        integer nvar,ivar,iunit,it,ierr

        integer ifileo
        integer icall
        data icall /0/
        save icall

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

