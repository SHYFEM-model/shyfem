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
c 06.06.2016    ggu     initialization from file changed
c 10.06.2016    ggu     some more re-formatting
c 08.09.2016    ggu     new decay function implemented (chapra), cleaned
c 13.02.2017    mic     idecay has new meaning!!! (incompatible)
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

	integer nvar,nbc,nintp,i,id,idc
	integer ishyff
	integer nmin
	real cdef(1)
	double precision dtime,dtime0

	logical has_restart
	logical has_output,next_output
	logical has_output_d,next_output_d
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
	idecay = getpar('idecay')
	ishyff = nint(getpar('ishyff'))

	dtime = t_act
	nvar = iconz
	allocate(tauv(nvar),cdefs(nvar),massv(nvar))
	cdefs = cref
	tauv = contau
	nmin = min(ndim_tau,nvar)
	if( nmin > 0 ) tauv(1:nmin) = taupar(1:nmin)

	if( idecay == 0 ) then 
	  write(6,*) 'no decay for tracer used'
	else if( idecay < 0 .or. idecay > 2 ) then 
	  write(6,*) 'no such option for decay: idecay = ',idecay
	  stop 'error stop tracer_init: no such option'
	else
	  write(6,*) 'decay for tracer used'
          write(6,*) 'idecay = ',idecay
	  write(6,*) '0 none   1 exp   2 chapra'
	  if( idecay == 1 ) then
	    write(6,*) 'decay parameter used: tauv ='
            write(6,*) tauv
	  end if
	end if  

        if( .not. has_restart(4) ) then	!no restart of conzentrations
	  if( nvar == 1 ) then 
	    call conz_init_file(dtime,nvar,nlvdi,nlv,nkn,cdefs,cnv)
	  else
	    call conz_init_file(dtime,nvar,nlvdi,nlv,nkn,cdefs,conzv)
	  end if
	end if

        call init_output('itmcon','idtcon',ia_out)
	if( ishyff == 1 ) ia_out = 0
	if( has_output(ia_out) ) then
          call open_scalar_file(ia_out,nlv,nvar,'con')
	  if( next_output(ia_out) ) then
	    if( nvar == 1 ) then
              idc = 10       !for tracer
	      call write_scalar_file(ia_out,idc,nlvdi,cnv)
	    else if( nvar > 1 ) then
	      do i=1,nvar
	        idc = 30 + i
	        call write_scalar_file(ia_out,idc,nlvdi,conzv(1,1,i))
	      end do
	    end if
	  end if
	end if

        call init_output_d('itmcon','idtcon',da_out)
        if( ishyff == 0 ) da_out = 0
        if( has_output_d(da_out) ) then
	  call shyfem_init_scalar_file('conz',nvar,.false.,id)
          da_out(4) = id
          if( next_output_d(da_out) ) then
	    if( nvar == 1 ) then
	      idc = 10
	      call shy_write_scalar_record(id,dtime,idc,nlvdi,cnv)
	    else
	      do i=1,nvar
	        idc = 30 + i
	        call shy_write_scalar_record(id,dtime,idc,nlvdi
     +						,conzv(1,1,i))
	      end do
            end if
          end if
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

	if( idecay == 1 ) then
          call decay_conz(dt,contau,cnv)
	else if( idecay == 2 ) then
          call decay_conz_chapra(dt,1.,cnv)
	end if

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
!$OMP&     dt,nlvdi,idecay) SHARED(conzv,tauv,massv)  DEFAULT(NONE)
 
          call scal_adv(what,i
     +                          ,conzv(1,1,i),idconz
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)


	  if(idecay == 1) then
            call decay_conz(dt,tauv(i),conzv(1,1,i))
	  else if( idecay == 2 ) then
            call decay_conz_chapra(dt,1.,conzv(1,1,i))
	  end if

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

	integer id,nvar,i,idc
        real cmin,cmax,ctot
	real v1v(nkn)
	double precision dtime

	logical next_output,next_output_d

	if( iconz < 0 ) return

c-------------------------------------------------------------
c write to file
c-------------------------------------------------------------

	dtime = t_act
	nvar = iconz

	if( next_output(ia_out) ) then
	  if( nvar == 1 ) then
            idc = 10       !for tracer
	    call write_scalar_file(ia_out,idc,nlvdi,cnv)
	  else if( nvar > 1 ) then
	    do i=1,nvar
	      idc = 30 + i
	      call write_scalar_file(ia_out,idc,nlvdi,conzv(1,1,i))
	    end do
	  end if
	end if

        if( next_output_d(da_out) ) then
	  id = nint(da_out(4))
	  if( nvar == 1 ) then
            idc = 10       !for tracer
	    call shy_write_scalar_record(id,dtime,idc,nlvdi,cnv)
	  else if( nvar > 1 ) then
	    do i=1,nvar
	      idc = 30 + i
	      call shy_write_scalar_record(id,dtime,idc,nlvdi
     +						,conzv(1,1,i))
	    end do
          end if
        end if

        call getinfo(ninfo)

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

        subroutine decay_conz(dt,tau,e)

c simulates decay for concentration

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        real, parameter :: alpha_t90 = 1./2.302585      !-1./ln(0.1)
        real, parameter :: alpha_exp = 1.               !e-folding time
        real, parameter :: alpha = alpha_exp            !tau is e-folding time
        !real, parameter :: alpha = alpha_t90           !tau is t90

        real dt				!time step in seconds
	real tau			!decay time in days (0 for no decay)
        real e(nlvdi,nkn)	        !state vector

        integer k,l,lmax
        real aux,tauaux

        if( tau .le. 0. ) return	!tau is 0 => no decay

	tauaux = tau * alpha
        aux = exp(-dt/(tauaux*86400))

	e = aux * e

        end

c*********************************************************************

        subroutine decay_conz_chapra(dt,tau,e)

c simulates decay for concentration (from Chapra, 506-510)

        use mod_layer_thickness
        use mod_ts
        use levels
        use basin, only : nkn,nel,ngr,mbw

        implicit none

        real dt                         !time step in seconds
        real tau                        !decay time in days (0 for no decay)
        real e(nlvdi,nkn)               !state vector

        integer k,l,lmax
	integer it,ith
        real aux,dtt
        real alpha,sd,ke,fp,ws
        real sr,srly,iaver
        real h,t,s
        real kb1,kbi,kbs,kb
	real eflux_top,eflux_bottom

	logical openmp_is_master

        !write(6,*) 'chapra: ',tau

        if( tau .le. 0. ) return        !tau is 0 => no decay

        alpha = 1.                      !proportionality constant
        sd = 0.65                       !Secchi-disk depth
        ke = 1.8 / sd                   !light extinction coefficient
        fp = 0.3                        !fraction of bacteria attached to part
        !fp = 0.0                        !fraction of bacteria attached to part
        ws = 0.0001                     !settling velocity [m/s]
        !m = 6.                         !suspended solids [mg/L]
        !ke = 0.55 * m

        dtt = dt / 86400                !change time step to days
        ws = ws * 86400                 !change settling velocity to [m/d]

        do k=1,nkn
          lmax = ilhkv(k)
          call meteo_get_solar_radiation(k,sr)  !solar radiation [W/m**2]
          srly = sr / 11.622                    !solar radiation [ly/h]
	  eflux_top = 0.
          do l=1,lmax
            t = tempv(l,k)
            s = saltv(l,k)
            h = hdknv(l,k)
            kb1 = (0.8+0.02*s) * 1.07**(t-20.)
            aux = exp(-ke*h)
            iaver = ( srly/(ke*h) ) * (1.-aux)
            srly = srly * aux                   !solar radiation at bottom
            kbi = alpha * iaver
            kbs = fp * ws / h
            kb = kb1 + kbi + kbs
            kb = kb1 + kbi
            e(l,k) = e(l,k) * exp(-dtt*kb)
	    eflux_bottom = e(l,k) * ( 1. - exp(-dtt*kbs) )
            e(l,k) = e(l,k) - eflux_bottom + eflux_top
	    eflux_top = eflux_bottom
            !if( k .eq. 100 .and. openmp_is_master() ) then
	      !call  openmp_get_thread_num(ith)
	      !call get_act_time(it)
              !write(333,*) it,1./kb
              !write(6,*) k,ith,1./kb
              !write(6,*) k,kb1,kbi,kbs,1./kb
              !write(6,*) sr,srly,aux,iaver
            !end if
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

        subroutine conz_init_file(dtime,nvar,nlvddi,nlv,nkn,val0,val)

c initialization of conz from file

	implicit none

	double precision dtime
	integer nvar
        integer nlvddi
        integer nlv
        integer nkn
        real val0(nvar)
        real val(nlvddi,nkn,nvar)

        call tracer_file_init('conz init','conzin',dtime
     +                          ,nvar,nlvddi,nlv,nkn,val0,val)

	end

c*********************************************************************
