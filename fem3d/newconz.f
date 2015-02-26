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
c
c*********************************************************************

	subroutine conz3sh

c shell for conz (new version)

	implicit none

c parameter
        include 'param.h'
c common
	include 'nbasin.h'
	include 'femtime.h'
	include 'mkonst.h'
	include 'nlevel.h'

	include 'diff_visc_fric.h'
	include 'aux_array.h'
	include 'bound_names.h'

	include 'conz.h'

c local
	logical binfo
        integer istot
	integer level
	integer nintp,nvar,ivar,id
	real cdef(1)
        real cmin,cmax,ctot
        real sindex
	real wsink
	real t,dt
	double precision dtime0,dtime
c function
	logical has_restart,next_output,has_output
	real getpar
c save & data
        character*4 what
        save what
	integer iconz
	save iconz
	real cref,rkpar
	save cref,rkpar
	real difmol
	save difmol
	real tau
	save tau
        integer ninfo
        save ninfo
        integer iprogr
        save iprogr
	integer idconz(nbcdim)
	save idconz
	integer ia_out(4)
	save ia_out

	integer icall
	save icall
	data icall /0/

	if(nlvdim.ne.nlvdi) stop 'error stop conz3sh: level dimension'

	if(icall.eq.-1) return

	binfo = .true.		! writes info to info file
	level = 0		! level > 0 -> writes level to extra file

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	if(icall.eq.0) then
	  iconz=nint(getpar('iconz'))
	  if( iconz .ne. 1 ) icall=-1
	  if(icall.eq.-1) return

	  cref=getpar('conref')
	  rkpar=getpar('chpar')
	  difmol=getpar('difmol')
	  tau = getpar('contau')

          what = 'conz'

          if( .not. has_restart(4) ) then	!no restart of conzentrations
	    call conini0(nlvdi,cnv,cref)
	  end if
	  call conz_init(itanf,nlv,nkn,cnv)	!read from file if name given

	  nvar = 1
          call init_output('itmcon','idtcon',ia_out)
	  if( has_output(ia_out) ) then
            call open_scalar_file(ia_out,nlv,nvar,'con')
	  end if

          call getinfo(ninfo)

	  dtime0 = itanf
	  nintp = 2
	  nvar = 1
	  cdef(1) = 0.
          call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +				,cdef,idconz)

	  iprogr = nint(getpar('iprogr'))
	  if( level .le. 0 ) iprogr = 0
	end if

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	icall=icall+1

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

        call decay_conz(1,1,1,dt,tau,cnv)

c-------------------------------------------------------------
c write to file
c-------------------------------------------------------------

	if( next_output(ia_out) ) then
          id = 10       !for tracer
	  call write_scalar_file(ia_out,id,nlvdi,cnv)
	end if

	if( iprogr .gt. 0 .and. mod(icall,iprogr) .eq. 0 ) then
	  call extract_level(nlvdim,nkn,level,cnv,v1v)
	  call wrnos2d_index(it,icall,'conz','concentration',v1v)
	end if

c-------------------------------------------------------------
c write to info file
c-------------------------------------------------------------

        if( binfo ) then
          call tsmass(cnv,+1,nlvdim,ctot)
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

	implicit none

c parameter
        include 'param.h'
c common
	include 'nbasin.h'
	include 'femtime.h'
	include 'mkonst.h'
	include 'nlevel.h'

	include 'diff_visc_fric.h'
	include 'aux_array.h'
	include 'bound_names.h'

	include 'conz.h'

c--------------------------------------------
c isabella
c--------------------------------------------
	integer ndim
	!parameter (ndim=10)
	!parameter (ndim=5)
	parameter (ndim=7)
	real massv(ndim)
	real taupar(ndim)
	save taupar
        !data taupar /1.,1.,1.,1.,2.,2.,1.,0.5,0.1,3./
        !data taupar /0.,0.,0.,0.,0./
        data taupar /0.,0.,0.,0.,0.,0.,0./
        !data taupar /1.,1.,1.,1.,1.,2.,2./

c local
        integer istot
	integer level
	integer nintp,nvar,ivar,i,id
	real cdef
        real cmin,cmax
        real sindex
	real wsink,tau,mass
	real t,dt
	real cdefs(ncsdim)
	double precision dtime,dtime0
c function
	logical has_restart,next_output
	real getpar
c save & data
        character*5 what
        save what
	integer iconz
	save iconz
	real cref,rkpar
	save cref,rkpar
	real difmol
	save difmol
	real contau
	save contau
        integer ninfo
        save ninfo
        integer iprogr
        save iprogr
	integer idconz(nbcdim)
	save idconz
	integer ia_out(4)
	save ia_out

	integer icall
	save icall
	data icall /0/

	if(nlvdim.ne.nlvdi) stop 'error stop conz3sh: level dimension'

	if(icall.eq.-1) return

c-------------------------------------------------------------
c initialization
c-------------------------------------------------------------

	if(icall.eq.0) then
	  iconz=nint(getpar('iconz'))
	  if( iconz .le. 1 ) icall=-1
	  if(icall.eq.-1) return
	  if( iconz .gt. ncsdim ) goto 99

	  nvar = iconz

	  cref=getpar('conref')
	  rkpar=getpar('chpar')
	  difmol=getpar('difmol')
	  contau = getpar('contau')

          what = 'conzm'

	  do i=1,nvar
	    cdefs(i) = 0.			!default boundary concentration
            if( .not. has_restart(4) ) then	!no restart of conzentrations
	      call conini0(nlvdi,conzv(1,1,i),cref)
	    end if
	  end do

          call init_output('itmcon','idtcon',ia_out)
          call open_scalar_file(ia_out,nlv,nvar,'con')

          call getinfo(ninfo)

	  dtime0 = itanf
	  nintp = 2
	  cdefs = 0.
          call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +				,cdefs,idconz)

	end if

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	icall=icall+1

	nvar = iconz
	wsink = 0.
	t = it
	dtime = it
	dt = idt

	call bnds_read_new(what,idconz,dtime)

!$OMP PARALLEL PRIVATE(i)
!$OMP DO SCHEDULE(DYNAMIC)

	do i=1,nvar
            call scal_adv(what,i
     +                          ,conzv(1,1,i),idconz
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)
	end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

c-------------------------------------------------------------
c simulate decay
c-------------------------------------------------------------

	do i=1,nvar
	  if( i .le. ndim ) then
	    tau = taupar(i)
	  else
	    tau = contau
	  end if
	  !write(6,*) 'decay : ',i,tau
          call decay_conz(ncsdim,nvar,i,dt,tau,conzv)
          !call decay_conz_variable(nsdim,nvar,ivar,dt,tau,e)
	end do

c-------------------------------------------------------------
c write to file
c-------------------------------------------------------------

	do i=1,nvar
	  call massconc(+1,conzv(1,1,i),nlvdim,mass)
	  if( i .le. ndim ) massv(i) = mass
	end do
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

	return
   99	continue
	write(6,*) 'iconz,ncsdim: ',iconz,ncsdim
	write(6,*) 'ncsdim must be greater or equal iconz'
	stop 'error stop conzm3sh: ncsdim'
	end

c*********************************************************************

        subroutine decay_conz(nsdim,nvar,ivar,dt,tau,e)

c simulates decay for concentration

        implicit none

        include 'param.h'

	real alpha_t90
	parameter( alpha_t90 = 1./2.302585 )	!-1./ln(0.1)

	integer nsdim			!dimension of state variables for e
	integer nvar			!actual number of state variables
	integer ivar			!state variable to use
        real dt				!time step in seconds
	real tau			!decay time in days (0 for no decay)
        real e(nlvdim,nkndim,nsdim)     !state vector

	include 'nbasin.h'
	include 'levels.h'

        integer k,l,i,lmax
        real aux,tauaux

        if( tau .le. 0. ) return	!tau is 0 => no decay

	!tauaux = tau			!tau is e-folding time
	tauaux = tau * alpha_t90	!tau is t_90
        aux = exp(-dt/(tauaux*86400))

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            e(l,k,ivar) = aux * e(l,k,ivar)
          end do
        end do

        end

c*********************************************************************

        subroutine decay_conz_variable(nsdim,nvar,ivar,dt,tau,e)

c simulates decay for concentration

        implicit none

        include 'param.h'

	real alpha_t90
	parameter( alpha_t90 = 1./2.302585 )	!-1./ln(0.1)

	integer nsdim			!dimension of state variables for e
	integer nvar			!actual number of state variables
	integer ivar			!state variable to use
        real dt				!time step in seconds
	real tau			!decay time in days (0 for no decay)
        real e(nlvdim,nkndim,nsdim)     !state vector

	include 'nbasin.h'
	include 'levels.h'

	include 'depth.h'
	include 'ts.h'

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
            e(l,k,ivar) = aux * e(l,k,ivar)
          end do
        end do

        end

c*********************************************************************

        subroutine conz_init(it,nlv,nkn,cnv)

c initialization of conz from file

        implicit none

        include 'param.h'

        integer it
        integer nlv
        integer nkn
        real cnv(nlvdim,1)

        character*80 conzf

        integer itc
        integer iuconz(3)

        call getfnm('conzin',conzf)

	itc = it

        if( conzf .ne. ' ' ) then
          write(6,*) 'conz_init: opening file for concentration'
          call ts_file_open(conzf,it,nkn,nlv,iuconz)
	  call ts_file_descrp(iuconz,'conz init')
          call ts_next_record(itc,iuconz,nkn,nlv,cnv)
          call ts_file_close(iuconz)
          write(6,*) 'concentration initialized from file ',conzf
        end if

	end

c*********************************************************************

