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
c
c*********************************************************************

	subroutine conz3sh

c shell for conz (new version)

	implicit none

c parameter
        include 'param.h'
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	real difv(0:nlvdim,1)
	common /difv/difv
        real difhv(nlvdim,1)
        common /difhv/difhv
	real v1v(1)
	common /v1v/v1v
        character*80 conzn(nbcdim)
        common /conzn/ conzn

        real cnv(nlvdim,nkndim)
        common /cnv/cnv
	save /cnv/	!not in ht anymore

c local
	logical binfo
        integer istot
	integer level
	integer nintp,nvar,ivar
	real cdef(1)
        real cmin,cmax,ctot
        real sindex
	real wsink
	real t,dt
c function
	logical has_restart
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
	integer iu,id,itmcon,idtcon
	save iu,id,itmcon,idtcon
        integer ninfo
        save ninfo
        integer iprogr
        save iprogr
	real bnd3_conz(nb3dim,0:nbcdim)
	save bnd3_conz

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

	  iu = 0
          id = 10       !for tracer
	  itmcon = nint(getpar('itmcon'))
	  idtcon = nint(getpar('idtcon'))
	  call adjust_itmidt(itmcon,idtcon)
	  call confop(iu,itmcon,idtcon,nlv,1,'con')

          call getinfo(ninfo)

	  nintp = 2
	  nvar = 1
	  cdef(1) = 0.
	  call bnds_init(what,conzn,nintp,nvar,nb3dim,bnd3_conz,cdef)
	  call bnds_set_def(what,nb3dim,bnd3_conz)

	  iprogr = nint(getpar('iprogr'))
	  if( level .le. 0 ) iprogr = 0
	end if

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	icall=icall+1

	wsink = 0.
	t = it
	dt = idt

        call scal_bnd(what,t,bnd3_conz)
        call scal_adv(what,0
     +                          ,cnv,bnd3_conz
     +                          ,rkpar,wsink
     +                          ,difhv,difv,difmol)

c-------------------------------------------------------------
c simulate decay
c-------------------------------------------------------------

        call decay_conz(1,1,1,dt,tau,cnv)

c-------------------------------------------------------------
c write to file
c-------------------------------------------------------------

	call confil(iu,itmcon,idtcon,id,nlvdi,cnv)

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
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	real difv(0:nlvdim,1)
	common /difv/difv
        real difhv(nlvdim,1)
        common /difhv/difhv
	real v1v(1)
	common /v1v/v1v
        character*80 conzn(nbcdim)
        common /conzn/ conzn

        real conzv(nlvdim,nkndim,ncsdim)
        common /conzv/conzv
	save /conzv/	!not in ht anymore

c--------------------------------------------
c isabella
c--------------------------------------------
	integer ndim
	!parameter (ndim=10)
	parameter (ndim=5)
	real massv(ndim)
	real taupar(ndim)
	save taupar
        !data taupar /1.,1.,1.,1.,2.,2.,1.,0.5,0.1,3./
        data taupar /0.,0.,0.,0.,0./

c local
        integer istot
	integer level
	integer nintp,nvar,ivar,i
	real cdef
        real cmin,cmax
        real sindex
	real wsink,tau,mass
	real t,dt
	real cdefs(ncsdim)
c function
	logical has_restart
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
	integer iu,id,itmcon,idtcon
	save iu,id,itmcon,idtcon
        integer ninfo
        save ninfo
        integer iprogr
        save iprogr
	real bnd3_conz(nb3dim,0:nbcdim)
	save bnd3_conz

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

	  iu = 0
	  itmcon = nint(getpar('itmcon'))
	  idtcon = nint(getpar('idtcon'))
	  call confop(iu,itmcon,idtcon,nlv,nvar,'com')

          call getinfo(ninfo)

	  nintp = 2
	  cdef = 0.
	  call bnds_init(what,conzn,nintp,nvar,nb3dim,bnd3_conz,cdefs)
	  !call bnds_set_def(what,nb3dim,bnd3_conz)

	end if

c-------------------------------------------------------------
c normal call
c-------------------------------------------------------------

	icall=icall+1

	nvar = iconz
	wsink = 0.
	t = it
	dt = idt

        call scal_bnd(what,t,bnd3_conz)
	!call bnds_print('debug conzm',nb3dim,bnd3_conz)

!$OMP PARALLEL PRIVATE(i)
!$OMP DO SCHEDULE(DYNAMIC)

	do i=1,nvar
            call scal_adv(what,i
     +                          ,conzv(1,1,i),bnd3_conz
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

	do i=1,nvar
	    id = 30 + i
	    call confil(iu,itmcon,idtcon,id,nlvdi,conzv(1,1,i))
	end do

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

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhkv(1)
        common /ilhkv/ilhkv

        integer k,l,i,lmax
        real aux,tauaux

        if( tau .le. 0. ) return	!tau is 0 => no decay

	!tauaux = tau			!e-folding time
	tauaux = tau * alpha_t90	!t_90
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

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhkv(1)
        common /ilhkv/ilhkv

        real hdknv(nlvdim,nkndim)
        common /hdknv/hdknv
        real saltv(nlvdim,1),tempv(nlvdim,1)
        common /saltv/saltv, /tempv/tempv
        real metrad(nkndim)
        common /metrad/metrad

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
	  solrad = metrad(k)		!solar radiation [W/m**2]
	  solrad = 500.
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

