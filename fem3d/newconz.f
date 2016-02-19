c
c $Id: newconz.f,v 1.7 2010-02-26 17:35:06 georg Exp $
c
c routines for generic concentration
c
c contents :
c
c subroutine conz3sh
c						shell for conz (new version)
c revision log :
c
c 28.04.2008    ggu     conz3sh into own file
c 28.04.2008    ggu     new conzm3sh for multiple concentrations
c 24.06.2008    ggu     changes in dacay for multiple concentrations
c 09.10.2008    ggu     new call to confop
c 19.01.2010    ggu     handle restart of conzentrations
c 25.02.2011    ggu     new routine decay_conz_variable(), add t90 time scale
c 13.02.2014    ggu     routines for reading initial condition
c 10.07.2014    ggu     only new file format allowed
c 20.10.2014    ggu     pass ids to scal_adv routines
c 10.02.2015    ggu     call to bnds_read_new() introduced
c 09.11.2015    ggu     newly structured in init, compute and write
c
c*********************************************************************

	subroutine tracer_compute

	use mod_conz

	implicit none

	if( iconz < 0 ) return

	if( iconz == 1 ) then
	  !call conz3sh
	  call tracer_compute_single
	else
	  !call conzm3sh
	  call tracer_compute_multi
	end if

	icall_conz = icall_conz + 1

	end

c*********************************************************************

	subroutine tracer_init

c initializes tracer computation

	use mod_conz
	!use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'femtime.h'

	integer nvar,nbc,nintp,i
	integer nmin
	real cdef(1)
	double precision dtime0

	logical has_restart,has_output
	integer nbnds
	real getpar

c-------------------------------------------------------------
c initialization of module
c-------------------------------------------------------------

	if( iconz < 0 ) return

        if( iconz == 0 ) then
          iconz=nint(getpar('iconz'))
          if( iconz <= 0 ) iconz = -1
          if( iconz < 0 ) return

          call mod_conz_init(iconz,nkn,nlvdi)

          write(6,*) 'tracer initialized: ',iconz,nkn,nlvdi
        end if

c-------------------------------------------------------------
c initialization of parameters
c-------------------------------------------------------------

	cref=getpar('conref')
	rkpar=getpar('chpar')
	difmol=getpar('difmol')
	contau = getpar('contau')

	nvar = iconz
	allocate(tauv(nvar),cdefs(nvar),massv(nvar))
	tauv = contau
	nmin = min(ndim_tau,nvar)
	if( nmin > 0 ) tauv(1:nmin) = taupar(1:nmin)

        if( .not. has_restart(4) ) then	!no restart of conzentrations
	  if( nvar == 1 ) then 
	    call conini0(nlvdi,cnv,cref)
	    call conz_init(itanf,nlvdi,nlv,nkn,cnv) !read from file
	  else
	    do i=1,nvar
	      call conini0(nlvdi,conzv(1,1,i),cref)
	    end do
	  end if
	end if

        call init_output('itmcon','idtcon',ia_out)
	if( has_output(ia_out) ) then
          call open_scalar_file(ia_out,nlv,nvar,'con')
	end if

        call getinfo(ninfo)

        nbc = nbnds()
        allocate(idconz(nbc))
        idconz = 0

	dtime0 = itanf
	nintp = 2
	cdefs = 0.				!default boundary condition
        call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +				,cdefs,idconz)

	iprogr = nint(getpar('iprogr'))
	if( level .le. 0 ) iprogr = 0

	end

c*********************************************************************

	subroutine tracer_compute_single

	use mod_conz
	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'femtime.h'
	include 'mkonst.h'

	real wsink
	real dt
	double precision dtime

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	if( iconz < 0 ) return

	if( it .eq. itanf ) stop 'tracer_compute_single: internal error'
	!if( it .eq. itanf ) return

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	wsink = 0.
	dtime = it
	dt = idt

	call bnds_read_new(what,idconz,dtime)

        call scal_adv(what,0
     +                          ,cnv,idconz
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

c-------------------------------------------------------------
c simulate decay
c-------------------------------------------------------------

        call decay_conz(dt,contau,cnv)
	call massconc(+1,cnv,nlvdi,massv(1))

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c*********************************************************************

	subroutine tracer_compute_multi

	use mod_conz
	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'femtime.h'
	include 'mkonst.h'

	integer nvar,i
	real wsink
	real dt
	double precision dtime

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	if( iconz < 0 ) return

	if( it .eq. itanf ) stop 'tracer_compute_multi: internal error'
	!if( it .eq. itanf ) return

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	nvar = iconz
	wsink = 0.
	dtime = it
	dt = idt

	call bnds_read_new(what,idconz,dtime)

	do i=1,nvar

!$OMP TASK FIRSTPRIVATE(i,rkpar,wsink,difhv,difv,difmol,idconz,what,
!$OMP&     dt,nlvdi) SHARED(conzv,tauv,massv)  DEFAULT(NONE)
 
          call scal_adv(what,i
     +                          ,conzv(1,1,i),idconz
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

          call decay_conz(dt,tauv(i),conzv(1,1,i))
	  call massconc(+1,conzv(1,1,i),nlvdi,massv(i))

!$OMP END TASK

	end do	

!$OMP TASKWAIT

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c*********************************************************************

	subroutine tracer_write

	use mod_conz
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'femtime.h'

	integer id,nvar,i
        real cmin,cmax,ctot
	real v1v(nkn)

	logical next_output

	if( iconz < 0 ) return

c-------------------------------------------------------------
c write to file
c-------------------------------------------------------------

	if( next_output(ia_out) ) then
	  if( iconz == 1 ) then
            id = 10       !for tracer
	    call write_scalar_file(ia_out,id,nlvdi,cnv)
	  else if( iconz > 1 ) then
	    nvar = iconz
	    do i=1,nvar
	      id = 30 + i
	      call write_scalar_file(ia_out,id,nlvdi,conzv(1,1,i))
	    end do
	  end if
	end if

c-------------------------------------------------------------
c write to info file
c-------------------------------------------------------------

	if( iconz == 1 ) then
	  if( iprogr .gt. 0 .and. mod(icall_conz,iprogr) .eq. 0 ) then
	    call extract_level(nlvdi,nkn,level,cnv,v1v)
	    call wrnos2d_index(it,icall_conz,'conz','concentration',v1v)
	  end if

          if( binfo ) then
	    ctot = massv(1)
            call conmima(nlvdi,cnv,cmin,cmax)
            write(ninfo,2021) 'conzmima: ',it,cmin,cmax,ctot
 2021       format(a,i10,2f10.4,e14.6)
          end if
	else
	  !write(65,*) it,massv
	end if

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************
C old routines ... not used anymore
c*********************************************************************
c*********************************************************************
c*********************************************************************

	subroutine conz3sh

c shell for conz (new version)

	use mod_conz
	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c common
	include 'femtime.h'
	include 'mkonst.h'

c local
        integer istot
	integer nintp,nvar,ivar,id
	integer nbc
	real cdef(1)
        real cmin,cmax,ctot
        real sindex
	real wsink
	real t,dt
	double precision dtime0,dtime
	real v1v(nkn)
c function
	logical has_restart,next_output,has_output
	integer nbnds
	real getpar

	if(icall_conz.eq.-1) return

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	if( it .eq. itanf ) stop 'conz3sh: internal error'
	!if( it .eq. itanf ) return

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	wsink = 0.
	t = it
	dtime = it
	dt = idt

	call bnds_read_new(what,idconz,dtime)

        call scal_adv(what,0
     +                          ,cnv,idconz
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

c-------------------------------------------------------------
c simulate decay
c-------------------------------------------------------------

        call decay_conz(dt,contau,cnv)

c-------------------------------------------------------------
c write to file
c-------------------------------------------------------------

	if( next_output(ia_out) ) then
          id = 10       !for tracer
	  call write_scalar_file(ia_out,id,nlvdi,cnv)
	end if

	if( iprogr .gt. 0 .and. mod(icall_conz,iprogr) .eq. 0 ) then
	  call extract_level(nlvdi,nkn,level,cnv,v1v)
	  call wrnos2d_index(it,icall_conz,'conz','concentration',v1v)
	end if

c-------------------------------------------------------------
c write to info file
c-------------------------------------------------------------

        if( binfo ) then
          call tsmass(cnv,+1,nlvdi,ctot)
          call conmima(nlvdi,cnv,cmin,cmax)
          write(ninfo,2021) 'conzmima: ',it,cmin,cmax,ctot
 2021     format(a,i10,2f10.4,e14.6)
        end if

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c*********************************************************************

	subroutine conzm3sh

c shell for conz with multi dimensions 

	use mod_conz
	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

c parameter
	include 'femtime.h'
	include 'mkonst.h'

c local
        integer istot
	integer nintp,nvar,ivar,i,id,nmin
	integer nbc
	real cdef
        real cmin,cmax
        real sindex
	real wsink,mass
	real t,dt
	double precision dtime,dtime0
c function
	integer nbnds
	logical has_restart,next_output
	real getpar

	if(icall_conz.eq.-1) return

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	if( it .eq. itanf ) return

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	nvar = iconz
	wsink = 0.
	t = it
	dtime = it
	dt = idt

	call bnds_read_new(what,idconz,dtime)

!$OMP PARALLEL
!$OMP SINGLE
	
	do i=1,nvar

!$OMP TASK FIRSTPRIVATE(i,rkpar,wsink,difhv,difv,difmol,idconz,what,
!$OMP&     dt,nlvdi,mass) SHARED(conzv,tauv,massv)  DEFAULT(NONE)
 
          call scal_adv(what,i
     +                          ,conzv(1,1,i),idconz
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

          call decay_conz(dt,tauv(i),conzv(1,1,i))
	  call massconc(+1,conzv(1,1,i),nlvdi,massv(i))

!$OMP END TASK
	end do	

!$OMP END SINGLE
!$OMP TASKWAIT
!$OMP END PARALLEL

c-------------------------------------------------------------
c write to file
c-------------------------------------------------------------

	!write(65,*) it,massv

	if( next_output(ia_out) ) then
	  do i=1,nvar
	    id = 30 + i
	    call write_scalar_file(ia_out,id,nlvdi,conzv(1,1,i))
	  end do
	end if

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************

        subroutine decay_conz(dt,tau,e)

c simulates decay for concentration

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'

	real alpha_t90
	!parameter( alpha_t90 = 1./2.302585 )	!-1./ln(0.1) - prob wrong
	parameter( alpha_t90 = 2.302585 )	!-1./ln(0.1)

        real dt				!time step in seconds
	real tau			!decay time in days (0 for no decay)
        real e(nlvdi,nkn)	        !state vector

        integer k,l,i,lmax
        real aux,tauaux

        if( tau .le. 0. ) return	!tau is 0 => no decay

	tauaux = tau			!tau is e-folding time
	!tauaux = tau * alpha_t90	!tau is t_90
        aux = exp(-dt/(tauaux*86400))

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            e(l,k) = aux * e(l,k)
          end do
        end do

        end

c*********************************************************************

        subroutine decay_conz_variable(dt,tau,e)

c simulates decay for concentration

	use mod_layer_thickness
	use mod_ts
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'

	real alpha_t90
	parameter( alpha_t90 = 1./2.302585 )	!-1./ln(0.1)

        real dt				!time step in seconds
	real tau			!decay time in days (0 for no decay)
        real e(nlvdi,nkn)               !state vector

        integer k,l,i,lmax
        real aux,dtt,rk,alpha
	real solrad,hdep,h,z
	real t,s,kappa

        if( tau .le. 0. ) return	!tau is 0 => no decay

	dtt = dt / 86400		!change time step to days

	rk = 1.				!extinction depth [m]
	alpha = 1./tau			!e-folding time
	!alpha = 1./(tau*alpha_t90)	!t_90

        do k=1,nkn
          lmax = ilhkv(k)
	  call meteo_get_solar_radiation(k,solrad)  !solar radiation [W/m**2]
	  !solrad = 500.
	  hdep = 0.
          do l=1,lmax
	    h = hdknv(l,k)
	    z = hdep + 0.5 * h
	    hdep = hdep + h
	    t = tempv(l,k)
	    s = saltv(l,k)
	    kappa = alpha * 1.040**(t-20.) * 1.012**s +
     +			0.113 * solrad * exp(-z/rk)
            aux = exp(-dtt*kappa)
            e(l,k) = aux * e(l,k)
          end do
        end do

        end

c*********************************************************************

        subroutine conz_init(it,nlvddi,nlv,nkn,cnv)

c initialization of conz from file

        implicit none

        include 'param.h'

        integer it
        integer nlvddi
        integer nlv
        integer nkn
        real cnv(nlvddi,1)

        character*80 conzf

        integer itc
        integer iuconz(3)

        call getfnm('conzin',conzf)

	itc = it

        if( conzf .ne. ' ' ) then
          write(6,*) 'conz_init: opening file for concentration'
          call ts_file_open(conzf,it,nkn,nlv,iuconz)
	  call ts_file_descrp(iuconz,'conz init')
          call ts_next_record(itc,iuconz,nlvddi,nkn,nlv,cnv)
          call ts_file_close(iuconz)
          write(6,*) 'concentration initialized from file ',conzf
        end if

	end

c*********************************************************************

