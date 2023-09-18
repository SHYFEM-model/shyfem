!
! $Id: newconz.f,v 1.7 2010-02-26 17:35:06 georg Exp $
!
! routines for generic concentration
!
! contents :
!
! subroutine conz3sh
!						shell for conz (new version)
! revision log :
!
! 28.04.2008    ggu     conz3sh into own file
! 28.04.2008    ggu     new conzm3sh for multiple concentrations
! 24.06.2008    ggu     changes in dacay for multiple concentrations
! 09.10.2008    ggu     new call to confop
! 19.01.2010    ggu     handle restart of conzentrations
! 25.02.2011    ggu     new routine decay_conz_variable(), add t90 time scale
! 13.02.2014    ggu     routines for reading initial condition
! 10.07.2014    ggu     only new file format allowed
! 20.10.2014    ggu     pass ids to scal_adv routines
! 10.02.2015    ggu     call to bnds_read_new() introduced
! 09.11.2015    ggu     newly structured in init, compute and write
!
!*********************************************************************
!----------------------------------------------------------------------------
        module conz_admin
!----------------------------------------------------------------------------
        contains
!----------------------------------------------------------------------------

	subroutine tracer_compute

	use conz_common

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

!*********************************************************************

	subroutine tracer_init

! initializes tracer computation

        use shympi
	use conz_common
	use conz_util
	!use mod_diff_visc_fric
	use levels, only : nlvdi,nlv,nlvmax
	use basin, only : nkn,nel,ngr,mbw
        use para
        use output
        use restart
        use bnd_admin
        use defnames
        use bnd_scalar

	implicit none

	include 'femtime.h'

	integer nvar,nbc,nintp,i
	integer levdbg
	integer nmin
	double precision cdef(1)
	double precision dtime0

!-------------------------------------------------------------
! initialization of module
!-------------------------------------------------------------

	if( iconz < 0 ) return

        if( iconz == 0 ) then
          iconz=nint(getpar('iconz'))
          if( iconz <= 0 ) iconz = -1
          if( iconz < 0 ) return

          call mod_conz_init(iconz,nkn,nlvdi)

          write(6,*) 'tracer initialized: ',iconz,nkn,nlvdi
        end if

!-------------------------------------------------------------
! initialization of parameters
!-------------------------------------------------------------

	cref=getpar('conref')
	rkpar=getpar('chpar')
	difmol=getpar('difmol')
	contau = getpar('contau')
	levdbg = nint(getpar('levdbg'))

	nvar = iconz
	allocate(tauv(nvar),cdefs(nvar),massv(nvar))
	tauv = contau
	nmin = min(ndim_tau,nvar)
	tauv(1:nmin) = taupar(1:nmin)

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
	binfo = levdbg > 0

        nbc = nbnds()
        allocate(idconz(nbc))
        idconz = 0

	dtime0 = itanf
	nintp = 2
	cdefs = 0.				!default boundary condition
        if(bmpi) then
          call bnds_init_mpi(what,dtime0,nintp,nvar,nkn,nlvmax,cdefs,idconz)
        else
          call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv,cdefs,idconz)
        end if

	iprogr = nint(getpar('iprogr'))
	if( level .le. 0 ) iprogr = 0

	end

!*********************************************************************

	subroutine tracer_compute_single

	use conz_common
	use diffusion
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
        use bnd_scalar
        use concentration

	implicit none

	include 'femtime.h'
	include 'mkonst.h'

	double precision wsink
	double precision dt
	double precision dtime

!-------------------------------------------------------------
! initialization
!-------------------------------------------------------------

	if( iconz < 0 ) return

	if( it .eq. itanf ) stop 'tracer_compute_single: internal error'
	!if( it .eq. itanf ) return

!-------------------------------------------------------------
! normal call
!-------------------------------------------------------------

	wsink = 0.
	dtime = it
	dt = idt

	call bnds_read_new(what,idconz,dtime)

        call scal_adv(what,0,cnv,idconz,rkpar,wsink,difhv,difv,difmol)

!-------------------------------------------------------------
! simulate decay
!-------------------------------------------------------------

        call decay_conz(dt,contau,cnv)
	if( binfo ) call massconc(+1,cnv,nlvdi,massv(1))

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!*********************************************************************

	subroutine tracer_compute_multi

	use conz_common
	use diffusion
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
        use bnd_scalar
        use concentration

	implicit none

	include 'femtime.h'
	include 'mkonst.h'

	integer nvar,i
	double precision wsink
	double precision dt
	double precision dtime

!-------------------------------------------------------------
! initialization
!-------------------------------------------------------------

	if( iconz < 0 ) return

	if( it .eq. itanf ) stop 'tracer_compute_multi: internal error'
	!if( it .eq. itanf ) return

!-------------------------------------------------------------
! normal call
!-------------------------------------------------------------

	nvar = iconz
	wsink = 0.
	dtime = it
	dt = idt

	call bnds_read_new(what,idconz,dtime)

	do i=1,nvar

!$OMP TASK FIRSTPRIVATE(i,rkpar,wsink,difhv,difv,difmol,idconz,what,    &
!$OMP   &     dt,nlvdi) SHARED(conzv,tauv,massv)  DEFAULT(NONE)
 
          call scal_adv(what,i,conzv(1,1,i),idconz,rkpar,wsink,difhv,difv,difmol)

          call decay_conz(dt,tauv(i),conzv(1,1,i))
	  if( binfo ) call massconc(+1,conzv(1,1,i),nlvdi,massv(i))

!$OMP END TASK

	end do	

!$OMP TASKWAIT

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!*********************************************************************

	subroutine tracer_write

	use conz_common
        use output
        use conz_util
        use nos_util
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'femtime.h'

	integer id,nvar,i
        double precision cmin,cmax,ctot
	double precision v1v(nkn)

	if( iconz < 0 ) return

!-------------------------------------------------------------
! write to file
!-------------------------------------------------------------

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

!-------------------------------------------------------------
! write to info file
!-------------------------------------------------------------

	if( iconz == 1 ) then
	  if( iprogr .gt. 0 .and. mod(icall_conz,iprogr) .eq. 0 ) then
	    call extract_level(nlvdi,nkn,level,cnv,v1v)
	    call wrnos2d_index(it,icall_conz,'conz','concentration',v1v)
	  end if

          if( binfo ) then
	    ctot = massv(1)
            call conmima(nlvdi,cnv,cmin,cmax)
!	    shympi FIXME : here exchange and reduce
            write(ninfo,2021) 'conzmima: ',it,cmin,cmax,ctot
 2021       format(a,i10,2f10.4,e14.6)
          end if
	else
	  !write(65,*) it,massv
	end if

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************
! old routines ... not used anymore
!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine conz3sh

! shell for conz (new version)

	use conz_common
	use conz_util
	use diffusion
        use output
        use restart
        use nos_util
        use check
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
        use para
        use bnd_scalar
        use concentration

	implicit none

! common
	include 'femtime.h'
	include 'mkonst.h'

! local
        integer istot
	integer nintp,nvar,ivar,id
	integer nbc
	double precision cdef(1)
        double precision cmin,cmax,ctot
        double precision sindex
	double precision wsink
	double precision t,dt
	double precision dtime0,dtime
	double precision v1v(nkn)
! function
	integer nbnds

	if(icall_conz.eq.-1) return

!-------------------------------------------------------------
! initialization
!-------------------------------------------------------------

	if( it .eq. itanf ) stop 'conz3sh: internal error'
	!if( it .eq. itanf ) return

!-------------------------------------------------------------
! normal call
!-------------------------------------------------------------

	wsink = 0.
	t = it
	dtime = it
	dt = idt

	call bnds_read_new(what,idconz,dtime)

        call scal_adv(what,0,cnv,idconz,rkpar,wsink,difhv,difv,difmol)

!-------------------------------------------------------------
! simulate decay
!-------------------------------------------------------------

        call decay_conz(dt,contau,cnv)

!-------------------------------------------------------------
! write to file
!-------------------------------------------------------------

	if( next_output(ia_out) ) then
          id = 10       !for tracer
	  call write_scalar_file(ia_out,id,nlvdi,cnv)
	end if

	if( iprogr .gt. 0 .and. mod(icall_conz,iprogr) .eq. 0 ) then
	  call extract_level(nlvdi,nkn,level,cnv,v1v)
	  call wrnos2d_index(it,icall_conz,'conz','concentration',v1v)
	end if

!-------------------------------------------------------------
! write to info file
!-------------------------------------------------------------

        if( binfo ) then
          call tsmass(cnv,+1,nlvdi,ctot)
          call conmima(nlvdi,cnv,cmin,cmax)
          write(ninfo,2021) 'conzmima: ',it,cmin,cmax,ctot
 2021     format(a,i10,2f10.4,e14.6)
        end if

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!*********************************************************************

	subroutine conzm3sh

! shell for conz with multi dimensions 

	use conz_common
	use conz_util
	use diffusion
        use output
        use restart
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
        use para
        use bnd_scalar
        use concentration

	implicit none

! parameter
	include 'femtime.h'
	include 'mkonst.h'

! local
        integer istot
	integer nintp,nvar,ivar,i,id,nmin
	integer nbc
	double precision cdef
        double precision cmin,cmax
        double precision sindex
	double precision wsink,mass
	double precision t,dt
	double precision dtime,dtime0
! function
	integer nbnds

	if(icall_conz.eq.-1) return

!-------------------------------------------------------------
! initialization
!-------------------------------------------------------------

	if( it .eq. itanf ) return

!-------------------------------------------------------------
! normal call
!-------------------------------------------------------------

	nvar = iconz
	wsink = 0.
	t = it
	dtime = it
	dt = idt

	call bnds_read_new(what,idconz,dtime)

!$OMP PARALLEL
!$OMP SINGLE
	
	do i=1,nvar

!$OMP TASK FIRSTPRIVATE(i,rkpar,wsink,difhv,difv,difmol,idconz,what,    &
!$OMP   &     dt,nlvdi,mass) SHARED(conzv,tauv,massv)  DEFAULT(NONE)
 
          call scal_adv(what,i,conzv(1,1,i),idconz,rkpar,wsink,difhv,difv,difmol)

          call decay_conz(dt,tauv(i),conzv(1,1,i))
	  call massconc(+1,conzv(1,1,i),nlvdi,massv(i))

!$OMP END TASK
	end do	

!$OMP END SINGLE
!$OMP TASKWAIT
!$OMP END PARALLEL

!-------------------------------------------------------------
! write to file
!-------------------------------------------------------------

	!write(65,*) it,massv

	if( next_output(ia_out) ) then
	  do i=1,nvar
	    id = 30 + i
	    call write_scalar_file(ia_out,id,nlvdi,conzv(1,1,i))
	  end do
	end if

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

        subroutine decay_conz(dt,tau,e)

! simulates decay for concentration

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'

	double precision alpha_t90
	!parameter( alpha_t90 = 1./2.302585 )	!-1./ln(0.1) - prob wrong
	parameter( alpha_t90 = 2.302585 )	!-1./ln(0.1)

        double precision dt				!time step in seconds
	double precision tau			!decay time in days (0 for no decay)
        double precision e(nlvdi,nkn)	        !state vector

        integer k,l,i,lmax
        double precision aux,tauaux

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

!*********************************************************************

        subroutine decay_conz_variable(dt,tau,e)

! simulates decay for concentration

	use layer_thickness
	use ts
	use levels
        use meteo_forcing
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        include 'param.h'

	double precision alpha_t90
	parameter( alpha_t90 = 1./2.302585 )	!-1./ln(0.1)

        double precision dt				!time step in seconds
	double precision tau			!decay time in days (0 for no decay)
        double precision e(nlvdi,nkn)               !state vector

        integer k,l,i,lmax
        double precision aux,dtt,rk,alpha
	double precision solrad,hdep,h,z
	double precision t,s,kappa

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
	    kappa = alpha * 1.040**(t-20.) * 1.012**s + 0.113 * solrad * exp(-z/rk)
            aux = exp(-dtt*kappa)
            e(l,k) = aux * e(l,k)
          end do
        end do

        end

!*********************************************************************

        subroutine conz_init(it,nlvddi,nlv,nkn,cnv)

! initialization of conz from file
        
        use para
        use tsfile_admin

        implicit none

        include 'param.h'

        integer it
        integer nlvddi
        integer nlv
        integer nkn
        double precision cnv(nlvddi,1)

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

!*********************************************************************

!----------------------------------------------------------------------------
        end module conz_admin
!----------------------------------------------------------------------------
