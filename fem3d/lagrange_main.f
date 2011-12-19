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

        integer ilagr
        integer itlanf,itlend,idtl,itlnext
	integer itmlgr,idtlgr,itnext
        integer iunit,uunit

        integer ifemop
        real getpar

        save itlanf,itlend,idtl,itlnext
	save itmlgr,idtlgr,itnext
        save iunit,uunit

	integer icall
	save icall
        data icall / 0 /
        
        if( icall .eq. -1 ) return
       
c---------------------------------------------------------------
c initialization
c---------------------------------------------------------------

        if( icall .eq. 0 ) then
          ilagr = nint(getpar('ilagr'))
          if( ilagr .eq. 0 ) icall = -1
          if( icall .eq. -1 ) return
	  icall = 1

          itmlgr=getpar('itmlgr')	!output to file
          idtlgr=getpar('idtlgr')	!output to file
          itlanf=getpar('itlanf')	!start of lagrangian
          itlend=getpar('itlend')	!end of lagrangian
          idtl=getpar('idtl') 		!frequency of new input

          artype=getpar('artype')	!store special type in common block
          rwhpar=getpar('rwhpar')

	  lunit = ifemop('.lgi','form','new')
	  if( lunit .le. 0 ) then
	    write(6,*) 'lunit = ',lunit
	    stop 'error stop lagrange: cannot open info file'
	  end if

	  nbdy = 0
	  idbdy = 0
	  tdecay = 0.

	  !call init_diff_oil

c         ------------------------------------------------------
c	  output
c         ------------------------------------------------------

	  if( itmlgr .lt. itanf ) itmlgr = itanf
	  itnext = itmlgr
	  if( idtlgr .le. 0 ) itnext = itend + 1

c         ------------------------------------------------------
c	  new input
c         ------------------------------------------------------

	  if( itlanf .lt. itanf ) itlanf = itanf
	  itlnext = itlanf
	  if( idtl .eq. 0 ) idtl = itend - itanf + 1	!release once at start
	  if( idtl .lt. 0 ) itlnext = itend + 1		!never release

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

	bback = .false.

c---------------------------------------------------------------
c new release of particles
c---------------------------------------------------------------

        if( it .ge. itlnext ) then
	  write(6,*) 'starting release of particles for lagrangian model'
	  call lgr_init_shell
	  itlnext = itlnext + idtl
	  if( itlnext .eq. itend ) itlnext = itend + 1
	  write(6,*) 'new particles released: ',nbdy,it
        end if           

c---------------------------------------------------------------
c one time step of particle tracking
c---------------------------------------------------------------

        call lagr_setup_timestep
	
	!call set_diff_oil
	call lagr_continuous_release_shell

 	call drogue

 	call lgr_larvae(it)

c---------------------------------------------------------------
c output
c---------------------------------------------------------------

        if( it .ge. itnext ) then
	  call lgr_output(iunit,it)
	  call lgr_output_concentrations
	  itnext = itnext + idtlgr
	end if

c---------------------------------------------------------------
c compress particles to save space
c---------------------------------------------------------------

	call compress_particles		!only after output of particles

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************

	subroutine drogue

c advection of particles

	implicit none
	
        include 'param.h'
	include 'lagrange.h'

	integer i,ie,nf
	real x,y,z
	real dt,ttime
                
        nf=0
	call get_timestep(dt)		!time to advect

!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED) PRIVATE(i,ie,x,y,z,ttime)
!$OMP& SCHEDULE(STATIC)
!$OMP& REDUCTION(+:nf)

	do i=1,nbdy

          call lagr_func(i)		!varying the variable typ(i)

	  x=x_body(i)
	  y=y_body(i)
	  z=z_body(i)
          ie=ie_body(i) 
	  ttime = dt

          if( z .lt. 1. ) call track_body(i,x,y,ie,ttime) 

	  x_body(i)=x
	  y_body(i)=y
          ie_body(i)=ie

	  if( ie .le. 0 ) nf = nf + 1

	end do

!$OMP  END PARALLEL DO

	write(lunit,*) 'lagrangian: (tot,out,in) ',nbdy,nf,nbdy-nf

	end	

c**********************************************************************

        subroutine track_body(i,x,y,iel,ttime)

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
	integer iel		!element number
	real x			!x-coordinate
	real y			!y-coordinate
	real ttime		!time to advect

        integer iarv(neldim)
        common /iarv/iarv
        
	integer nl
	real xn,yn

        if(iel.le.0) return	!particle out of domain

c---------------------------------------------------------------
c initialize
c---------------------------------------------------------------

        nl = 100		!maximum loop count

	xn = x
	yn = y

c---------------------------------------------------------------
c track particle
c---------------------------------------------------------------

        call track_orig(ttime,i,iel,xn,yn)

	do while ( ttime.gt.0. .and. iel.gt.0 .and. nl.gt.0 )
           call track_line(ttime,i,iel,xn,yn)
           nl = nl - 1
        end do

c---------------------------------------------------------------
c error condition (infinite loop)
c---------------------------------------------------------------

        if( nl .eq. 0 ) then
                print*, 'inifinite loop in track_line'
                print*, i,iel,xn,yn
                iel = -iel
        end if
        ttime = 0.

c---------------------------------------------------------------
c diffusion
c---------------------------------------------------------------

        if( .not. bback .and. iel .gt. 0 .and. rwhpar .gt. 0 ) then
           call lag_diff(iel,i,xn,yn)
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

