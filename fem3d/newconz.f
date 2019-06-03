
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

c routines for generic concentration
c
c contents :
c
c subroutine conz3sh
c						shell for conz (new version)
c revision log :
c
c 28.04.2008	ggu	conz3sh into own file
c 28.04.2008	ggu	new conzm3sh for multiple concentrations
c 24.06.2008	ggu	changes in dacay for multiple concentrations
c 09.10.2008	ggu	new call to confop
c 19.01.2010	ggu	handle restart of conzentrations
c 23.03.2010	ggu	changed v6.1.1
c 25.02.2011	ggu	new routine decay_conz_variable(), add t90 time scale
c 01.03.2011	ggu	changed VERS_6_1_20
c 30.03.2012	ggu	changed VERS_6_1_51
c 01.06.2012	ggu	changed VERS_6_1_53
c 29.08.2012	ggu	changed VERS_6_1_56
c 12.11.2013	ggu	changed VERS_6_1_69
c 13.02.2014	ggu	routines for reading initial condition
c 07.03.2014	ggu	changed VERS_6_1_72
c 07.07.2014	ggu	changed VERS_6_1_79
c 10.07.2014	ggu	only new file format allowed
c 18.07.2014	ggu	changed VERS_7_0_1
c 20.10.2014	ggu	pass ids to scal_adv routines
c 30.10.2014	ggu	changed VERS_7_0_4
c 26.11.2014	ggu	changed VERS_7_0_7
c 19.12.2014	ggu	changed VERS_7_0_10
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 10.02.2015	ggu	call to bnds_read_new() introduced
c 26.02.2015	ggu	changed VERS_7_1_5
c 01.04.2015	ggu	changed VERS_7_1_7
c 13.07.2015	ggu	changed VERS_7_1_51
c 17.07.2015	ggu	changed VERS_7_1_52
c 17.07.2015	ggu	changed VERS_7_1_53
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 24.07.2015	ggu	changed VERS_7_1_82
c 30.07.2015	ggu	changed VERS_7_1_83
c 18.09.2015	ggu	changed VERS_7_2_3
c 23.09.2015	ggu	changed VERS_7_2_4
c 22.10.2015	ggu	changed VERS_7_3_7
c 09.11.2015	ggu	newly structured in init, compute and write
c 19.02.2016	ggu	changed VERS_7_5_3
c 06.06.2016	ggu	initialization from file changed
c 10.06.2016	ggu	some more re-formatting
c 08.09.2016	ggu	new decay function implemented (chapra), cleaned
c 13.02.2017	mcg	idecay has new meaning!!! (incompatible)
c 13.04.2017	ggu	contau deprecated... use taupar (array)
c 02.09.2017	ggu	changed VERS_7_5_31
c 09.10.2017	ggu	changed VERS_7_5_33
c 05.12.2017	ggu	changed VERS_7_5_39
c 22.02.2018	ggu	changed VERS_7_5_42
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.10.2018	ggu	changed VERS_7_5_50
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c
c*********************************************************************

	subroutine tracer_init

c initializes tracer computation

	use mod_conz
	!use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
	use para

	implicit none

	integer nvar,nbc,nintp,i,id,idc
	integer levdbg
	integer ishyff
	integer n
	real, allocatable :: aux(:)
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

	  call tracer_accum_init

          write(6,*) 'tracer initialized: ',iconz,nkn,nlvdi
        end if

c-------------------------------------------------------------
c initialization of parameters
c-------------------------------------------------------------

	cref=getpar('conref')
	rkpar=getpar('chpar')
	difmol=getpar('difmol')
	idecay = getpar('idecay')
	!baccum = nint(getpar('iconza')) /= 0
	baccum = .false.
	levdbg = nint(getpar('levdbg'))

	call get_act_dtime(dtime)
	nvar = iconz
	allocate(tauv(nvar),cdefs(nvar),massv(nvar))
	cdefs = cref

	call para_get_array_size('taupar',n)
	if( n > 1 .and. n /= nvar ) then
	  write(6,*) 'array has wrong size: ','taupar'
	  write(6,*) 'size should be 1 or ',nvar
	  write(6,*) 'size from STR file is ',n
	  allocate(aux(n))
	  call para_get_array_value('taupar',n,n,aux)
          write(6,*) aux
	  stop 'error stop tracer_init: wrong array size'
	end if
	call para_get_array_value('taupar',nvar,n,tauv)
	contau = tauv(1)
	if( n == 1 ) tauv(:) = contau

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

	do i=1,nvar
	  call massconc(+1,cnv,nlvdi,massv(i))
	end do

	call tracer_write_init
	call tracer_write

        call getinfo(ninfo)
	binfo = levdbg > 0
	binfo = .true.

        nbc = nbnds()
        allocate(idconz(nbc))
        idconz = 0

	call get_first_dtime(dtime0)
	nintp = 2
	cdefs = 0.				!default boundary condition
        call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +				,cdefs,idconz)

	iprogr = nint(getpar('iprogr'))
	if( level .le. 0 ) iprogr = 0

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************
c compute tracer
c*********************************************************************
c*********************************************************************
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

	subroutine tracer_compute_single

	use mod_conz
	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'mkonst.h'

	logical bfirst
	real wsink
	real dt
	double precision dtime

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	if( iconz < 0 ) return

	call is_time_first(bfirst)
	if( bfirst ) stop 'tracer_compute_single: internal error'

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	wsink = 0.
	call get_act_dtime(dtime)
	call get_timestep(dt)

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

	if( binfo ) call massconc(+1,cnv,nlvdi,massv(1))

c-------------------------------------------------------------
c accumulate
c-------------------------------------------------------------

	call tracer_accum_accum(dt)

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

	include 'mkonst.h'

	logical blinfo,bfirst
	integer nvar,i
	real wsink
	real dt
	double precision dtime

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	if( iconz < 0 ) return

	call is_time_first(bfirst)
	if( bfirst ) stop 'tracer_compute_multi: internal error'

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	nvar = iconz
	wsink = 0.
	call get_act_dtime(dtime)
	call get_timestep(dt)
	blinfo = binfo

	call bnds_read_new(what,idconz,dtime)

	do i=1,nvar

!$OMP TASK FIRSTPRIVATE(i,rkpar,wsink,difhv,difv,difmol,idconz,what,
!$OMP& dt,nlvdi,idecay,blinfo) SHARED(conzv,tauv,massv) DEFAULT(NONE)
 
          call scal_adv(what,i
     +                          ,conzv(1,1,i),idconz
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)


	  if(idecay == 1) then
            call decay_conz(dt,tauv(i),conzv(1,1,i))
	  else if( idecay == 2 ) then
            call decay_conz_chapra(dt,1.,conzv(1,1,i))
	  end if

	  if( blinfo ) call massconc(+1,conzv(1,1,i),nlvdi,massv(i))

!$OMP END TASK

	end do	

!$OMP TASKWAIT

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	end

c*********************************************************************
c*********************************************************************
c*********************************************************************
c write and read routines
c*********************************************************************
c*********************************************************************
c*********************************************************************

	subroutine tracer_write_init

	use mod_conz

	implicit none

	integer nvar,id
	logical has_output_d

        call init_output_d('itmcon','idtcon',da_out)

	nvar = iconz

        if( has_output_d(da_out) ) then
	  call shyfem_init_scalar_file('conz',nvar,.false.,id)
          da_out(4) = id
	end if

	end

c*********************************************************************

	subroutine tracer_write

	use mod_conz
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

	integer id,nvar,i,idc
        real cmin,cmax,ctot
	real v1v(nkn)
	double precision dtime
	character*20 aline
	real, allocatable :: caux2d(:,:)

	logical next_output,next_output_d

	if( iconz < 0 ) return

c-------------------------------------------------------------
c write to file
c-------------------------------------------------------------

	call get_act_dtime(dtime)
	nvar = iconz

        if( next_output_d(da_out) ) then
	  id = nint(da_out(4))
	  if( nvar == 1 ) then
            idc = 10       !for tracer
	    if( baccum ) then
	write(6,*) 'writing accumulated...',dtconz_accum
	      call tracer_accum_aver
	      allocate(caux2d(nlvdi,nkn))
	      caux2d = conz_min(:,:,1)
	write(6,*) 'writing accumulated...',id,idc,dtime
	write(6,*) 'writing value...',(caux2d(1,i),i=1,nkn,nkn/10)
	write(6,*) 'writing minmax...',minval(caux2d),maxval(caux2d)
	      call shy_write_scalar_record(id,dtime,idc,nlvdi
     +						,caux2d)
	      caux2d = conz_aver(:,:,1)	!convert from double to real
	      dtime = dtime + 1
	write(6,*) 'writing accumulated...',id,idc,dtime
	write(6,*) 'writing value...',(caux2d(1,i),i=1,nkn,nkn/10)
	      call shy_write_scalar_record(id,dtime,idc,nlvdi
     +						,cnv)
!     +						,caux2d)
	      caux2d = conz_max(:,:,1)
	      dtime = dtime + 1
	write(6,*) 'writing accumulated...',id,idc,dtime
	write(6,*) 'writing value...',(caux2d(1,i),i=1,nkn,nkn/10)
	      call shy_write_scalar_record(id,dtime,idc,nlvdi
     +						,caux2d)
	      call tracer_accum_init
	    else
	      call shy_write_scalar_record(id,dtime,idc,nlvdi,cnv)
	    end if
	  else if( nvar > 1 ) then
	    do i=1,nvar
	      idc = 300 + i
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
	    stop 'error stop tracer_write: iprogr not supported'
	    !call extract_level(nlvdi,nkn,level,cnv,v1v)
	    !call wrnos2d_index(it,icall_conz,'conz','concentration',v1v)
	  end if

          if( binfo ) then
	    ctot = massv(1)
            call conmima(nlvdi,cnv,cmin,cmax)
	    cmin = shympi_min(cmin)
	    cmax = shympi_max(cmax)
	    call get_act_timeline(aline)
            write(ninfo,2021) ' conzmima: ',aline,cmin,cmax,ctot
 2021       format(a,a20,2f10.4,e14.6)
          end if
	else
	  !write(65,*) it,massv
	end if

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

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
c*********************************************************************
c*********************************************************************
c decay routines
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
        fp = 0.02                        !fraction of bacteria attached to part
        !fp = 0.3                        !fraction of bacteria attached to part
        !fp = 0.2                        !fraction of bacteria attached to part
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
            !kb = kb1 + kbi
            e(l,k) = e(l,k) * exp(-dtt*kb)
	    eflux_bottom = e(l,k) * ( 1. - exp(-dtt*kbs) )
            e(l,k) = e(l,k) - eflux_bottom + eflux_top
	    eflux_top = eflux_bottom
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
c*********************************************************************
c*********************************************************************
c tracer accumulation routines
c*********************************************************************
c*********************************************************************
c*********************************************************************

	subroutine tracer_accum_init

	use mod_conz

	implicit none

	real, parameter :: high = 1.e+30

	if( .not. baccum ) return

	dtconz_accum = 0.
	conz_min(:,:,:) = high	!min
	conz_aver(:,:,:) = 0.	!aver
	conz_max(:,:,:) = -high	!max

	end

c*********************************************************************

	subroutine tracer_accum_accum(dt)

	use mod_conz

	implicit none

	real dt

	if( .not. baccum ) return

	dtconz_accum = dtconz_accum + dt

	if( iconz == 1 ) then
	  where( cnv < conz_min(:,:,1) ) conz_min(:,:,1) = cnv
	  where( cnv > conz_max(:,:,1) ) conz_max(:,:,1) = cnv
	  !conz_min(:,:,1) = min(conz_min(:,:,1),cnv)
	  conz_aver(:,:,1) = conz_aver(:,:,1) + cnv * dt
	  !conz_max(:,:,1) = max(conz_max(:,:,1),cnv)
	else if( iconz > 1 ) then
	  conz_min = min(conz_min,conzv)
	  conz_aver = conz_aver + conzv * dt
	  conz_max = max(conz_max,conzv)
	end if

	end

c*********************************************************************

	subroutine tracer_accum_aver

	use mod_conz

	implicit none

	if( .not. baccum ) return

	if( dtconz_accum == 0. ) return

	conz_aver = conz_aver / dtconz_accum

	end

c*********************************************************************

