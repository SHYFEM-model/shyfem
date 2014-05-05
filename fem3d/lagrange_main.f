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
c       subroutien nbody
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
c         subroutine setup_fx           set up flux2d(3,neldim)
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
c
c****************************************************************            

	subroutine lagrange

c lagranian main routine

	implicit none

        include 'param.h'
        include 'lagrange.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	logical brelease
	logical bcompres
	logical boilsim
	logical blarvae
	save brelease,bcompres,boilsim,blarvae

        integer ilagr
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

	bcompres = .false.	!compress particles
	bcompres = .true.	!compress particles

c---------------------------------------------------------------
c initialization
c---------------------------------------------------------------

        if( icall .eq. 0 ) then
          ilagr = nint(getpar('ilagr'))
          if( ilagr .eq. 0 ) icall = -1
          if( icall .eq. -1 ) return
	  icall = 1

          itmlgr=getpar('itmlgr')	!startime write output to file
          idtlgr=getpar('idtlgr')	!frequency output to file
          itlanf=getpar('itlanf')	!start of lagrangian sim
          itlend=getpar('itlend')	!end of lagrangian sim
          idtl=getpar('idtl') 		!frequency of release
          itranf=getpar('itranf') 	!time of initial continuous release
          itrend=getpar('itrend') 	!time of final continuous release

          artype=getpar('artype')	!store special type in common block
          rwhpar=getpar('rwhpar')	!lagrangian diffusion

          ldecay=getpar('ldecay')	!decay time for particles
          boilsim=nint(getpar('ioil')).gt.0
          blarvae=nint(getpar('ilarv')).gt.0

	  lunit = ifemop('.lgi','form','new') !unit for lagrangian info
	  if( lunit .le. 0 ) then
	    write(6,*) 'lunit = ',lunit
	    stop 'error stop lagrange: cannot open info file'
	  end if

	  nbdy = 0 		!number of particles to insert
	  idbdy = 0		!id body unique 
	  tdecay = 0		!not used anymore

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
c connectivity module ?
c---------------------------------------------------------------

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

c---------------------------------------------------------------
c compress to save space: only in contiunous release mode! 
c---------------------------------------------------------------

	if( bcompres ) then 
	  call compress_particles	!only after output of particles
	endif 

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************

	subroutine drogue

	implicit none
	
        include 'param.h'
	include 'lagrange.h'

	integer nf,i,ii,n
	integer chunk
	real dt

	integer ndim
	parameter(ndim=100)
	integer ic(0:ndim)

	chunk = 100
        nf=0
	call get_timestep(dt)		!time to advect

	call openmp_get_max_threads(n)
	if( n .gt. ndim ) stop 'error stop drogue: too many processes'

	do ii=0,n
	  ic(ii) = 0
	end do

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ii)
!$OMP DO SCHEDULE(DYNAMIC,chunk)

	do i=1,nbdy
	  call openmp_get_thread_num(ii)
	  ic(ii) = ic(ii) + 1
	  call track_single(i,dt)
	end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

c	write(lunit,*) 'lagrangian: (tot,out,in) ',nbdy,nf,nbdy-nf
c	write(lunit,'(a,i10,12i5)') 'parallel: ',nbdy,(ic(ii),ii=0,n-1)

	end

c**********************************************************************

	subroutine track_single(i,dt)

c advection of particles

	implicit none
	
        include 'param.h'
	include 'lagrange.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer i,id,ie,nf
	real x,y,z
	real dt,ttime,tmax
                
c       call lagr_func(i)		!varying the variable typ(i)
c	call lagr_surv(i)

	x=x_body(i)
	y=y_body(i)
	z=z_body(i)
        ie = ie_body(i)
	id = id_body(i) 

	tmax = it - tin(i) 
	if( tmax .lt. 0. ) stop 'error stop drogue: internal error'
	ttime = min(tmax,dt)

	if(z.ge.0.and.z.le.1) call track_body(i,id,x,y,ie,ttime) 

	x_body(i)=x
	y_body(i)=y
        ie_body(i)=ie

	end	

c**********************************************************************

        subroutine track_body(i,id,x,y,iel,ttime)

c tracks one particle
c
c uses two routines:
c
c TRACK_ORIG if the particle is inside an element (first call)
c
c TRACK_LINE if the particle is on one side (normal situation)

	implicit none

	include 'param.h'
	include 'lagrange.h'
	
	integer i		!particle number
	integer id		!particle id
	integer iel		!element number
	integer ielem		!element number to check
	real x			!x-coordinate
	real y			!y-coordinate
	real ttime		!time to advect

        integer iarv(neldim)
        common /iarv/iarv
        
	integer nl
	integer ltbdy
	integer ieold,ieorig
	real torig
	real xn,yn

        if(iel.le.0) return	!particle out of domain

c---------------------------------------------------------------
c initialize
c---------------------------------------------------------------

        nl = 100		!maximum loop count
	ltbdy = 0

	xn = x
	yn = y

c---------------------------------------------------------------
c track particle
c---------------------------------------------------------------

	torig = ttime
	ieold = iel	!element the particle is in or leaving
	ieorig = iel	!original element the particle was in
        call track_orig(ttime,id,iel,xn,yn,ltbdy)
	call lagr_connect_count(i,ieold,ieorig,torig-ttime,0)
	call lagr_count(i,ieold,torig-ttime,0)

	do while ( ttime.gt.0. .and. iel.gt.0 .and. nl.gt.0 )
	  torig = ttime
	  ieold = iel
          call track_line(ttime,id,iel,xn,yn,ltbdy)
          nl = nl - 1
	  call lagr_connect_count(i,ieold,ieorig,torig-ttime,1)
	  call lagr_count(i,ieold,torig-ttime,1)
	  ieorig = ieold
        end do

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
	  call lagr_connect_count(i,iel,ieorig,ttime,1)
        end if     

c---------------------------------------------------------------
c special treatment and finish up
c---------------------------------------------------------------

        if( .not. bback ) then
          if( iel.gt.0 .and. iarv(iel).eq.artype ) iel = -iel
	end if

	x = xn
	y = yn

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine lagr_count_init

	implicit none

	include 'param.h'
	include 'lagrange.h'

	integer ie

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	do ie=1,nel
	  i_count(ie) = 0
	  t_count(ie) = 0.
	end do

	end

c**********************************************************************

	subroutine lagr_count(i,ie,time,icc)

	implicit none

	include 'param.h'
	include 'lagrange.h'

	integer i
	integer ie
	real time
	integer icc

	if( ie .le. 0 ) return

	i_count(ie) = i_count(ie) + 1
	t_count(ie) = t_count(ie) + time

	end

c**********************************************************************

	subroutine lagr_count_out(it)

	implicit none

	include 'param.h'
	include 'lagrange.h'

	integer it

	integer ie,iu
	character*80 file

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	iu = 237
	file = 'lagr_count_out.txt'
	open(iu,file=file,status='unknown',form='formatted')
	write(iu,*) it,nel
	do ie=1,nel
	  i_count(ie) = 0
	  t_count(ie) = 0.
	  write(iu,*) ie,i_count(ie),t_count(ie)
	end do
	close(iu)

	end

c**********************************************************************

