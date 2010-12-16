c $Id: newbcl.f,v 1.37 2010-03-08 17:46:45 georg Exp $
c
c baroclinic routines
c
c contents :
c
c subroutine barocl(mode)		amministrates the baroclinic time step
c subroutine rhoset_shell		sets rho iterating to real solution
c subroutine rhoset(resid)		computes rhov and bpresv
c subroutine convectivecorr             convective adjustment
c subroutine getts(l,k,t,s)             accessor routine to get T/S
c
c revision log :
c
c revised 30.08.95	$$AUST - austausch coefficient introduced
c revised 11.10.95	$$BCLBND - boundary condition for barocliic runs
c 19.08.1998    ggu     call to barcfi changed
c 20.08.1998    ggu     can initialize S/T from file
c 24.08.1998    ggu     levdbg used for debug
c 26.08.1998    ggu     init, bnd and file routines substituted with con..
c 30.01.2001    ggu     eliminated compile directives
c 05.12.2001    ggu     horizontal diffusion variable, limit diffusion coef.
c 05.12.2001    ggu     compute istot, more debug info
c 11.10.2002    ggu     diffset introduced, shpar = thpar
c 10.08.2003    ggu     qfluxr eliminated (now in subn11.f)
c 10.08.2003    ggu     rhov and bpresv are initialized here
c 04.03.2004    ggu     in init for T/S pass number of vars (inicfil)
c 15.03.2004    ggu     general clean-up, bclint() deleted, new scal3sh
c 17.01.2005    ggu     new difhv implemented
c 15.03.2005    ggu     new diagnostic routines implemented (diagnostic)
c 15.03.2005    ggu     new 3d boundary conditions implemented
c 05.04.2005    ggu     some changes in routine diagnostic
c 07.11.2005    ggu     sinking velocity wsink introduced in call to scal3sh
c 08.06.2007    ggu&deb restructured for new baroclinic version
c 04.10.2007    ggu     bug fix -> call qflux3d with dt real
c 17.03.2008    ggu     new open boundary routines introduced
c 08.04.2008    ggu     treatment of boundaries slightly changed
c 22.04.2008    ggu     advection parallelized, no saux1v...
c 23.04.2008    ggu     call to bnds_set_def() changed
c 12.06.2008    ggu     s/tdifhv deleted
c 09.10.2008    ggu     new call to confop
c 12.11.2008    ggu     new initialization, check_layers, initial nos file
c 13.01.2009    ggu&deb changes in reading file in read_next_record()
c 13.10.2009    ggu     in rhoset bug computing pres
c 13.11.2009    ggu     only initialize T/S if no restart, new rhoset_shell
c 19.01.2010    ggu     different call to has_restart() 
c 16.12.2010    ggu     sigma layers introduced (maybe not finished)
c
c*****************************************************************

	subroutine barocl(mode)

c amministrates the baroclinic time step
c
c mode : =0 initialize  >0 normal call
c
c written 09.01.94 by ggu  (from scratch)
c
	implicit none
c
c parameter
	include 'param.h'
c arguments
	integer mode
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	real saltv(nlvdim,1),tempv(nlvdim,1),rhov(nlvdim,1)
	real bpresv(nlvdim,1)
        real bpresxv(nlvdim,1),bpresyv(nlvdim,1)!deb110407
	real uprv(nlvdim,1), vprv(nlvdim,1)
	real v1v(1),v2v(1)
	integer ilhkv(1)
	common /saltv/saltv, /tempv/tempv, /rhov/rhov
	common /bpresv/bpresv
        common /bpresxv/bpresxv, /bpresyv/bpresyv!deb110407
	common /uprv/uprv, /vprv/vprv
	common /ilhkv/ilhkv
	common /v1v/v1v, /v2v/v2v
	real zeov(3,1), zenv(3,1)
	common /zeov/zeov, /zenv/zenv
	real difv(0:nlvdim,1)
	common /difv/difv
        real hdkov(nlvdim,1)
        common /hdkov/hdkov
        real difhv(nlvdim,1)
        common /difhv/difhv
        real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
        integer ib
                      
        character*80 saltn(nbcdim)
        character*80 tempn(nbcdim)
        common /saltn/ saltn
        common /tempn/ tempn

c local
	logical debug
	logical badvect
	logical bgdebug
        logical binfo
        logical bstop
	logical binitial_nos
	integer levdbg
	integer ie
	integer ibarcl
	integer idtext,itmext
	integer imin,imax
	integer nintp,nvar,ivar
	real cdef,t
	real xmin,xmax
        integer itemp,isalt
	real salref,temref,sstrat,tstrat
	real shpar,thpar
	real difmol
        real s
	real dt
        real gamma,gammax
	real mass
	real wsink
	integer isact,l,k,lmax
	integer kspec
	integer icrst
	real stot,ttot,smin,smax,tmin,tmax,rmin,rmax
	double precision v1,v2,mm
	character*4 what
c functions
c	real sigma
	real getpar
	double precision scalcont,dq
	integer iround
	logical has_restart

	integer tid
c	integer OMP_GET_THREAD_NUM
	
	double precision theatold,theatnew
	double precision theatconv1,theatconv2,theatqfl1,theatqfl2
	real cw,row
c save
        real bnd3_temp(nb3dim,0:nbcdim)
        save bnd3_temp
        real bnd3_salt(nb3dim,0:nbcdim)
        save bnd3_salt

        integer iu,itmcon,idtcon
        save iu,itmcon,idtcon

        integer ninfo
        save ninfo

	save badvect
	save salref,temref
	save difmol
        save itemp,isalt
c data
	integer icall
	save icall
	data icall /0/

	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension in barocl'

	if(icall.eq.-1) return

	levdbg = nint(getpar('levdbg'))
	debug = levdbg .ge. 3
	binfo = .true.
        bgdebug = .false.
	binitial_nos = .true.

c initialization

	if(icall.eq.0) then	!first time

		ibarcl=iround(getpar('ibarcl'))
		if(ibarcl.le.0) icall = -1
		if(ibarcl.gt.3) goto 99
		if(icall.eq.-1) return

		badvect = ibarcl .eq. 1 .or. ibarcl .eq. 3

		salref=getpar('salref')
		temref=getpar('temref')
		sstrat=getpar('sstrat')
		tstrat=getpar('tstrat')
		difmol=getpar('difmol')
		idtcon=iround(getpar('idtcon'))
		itmcon=iround(getpar('itmcon'))
                itemp=iround(getpar('itemp'))
                isalt=iround(getpar('isalt'))

		!dt=idt

c		--------------------------------------------
c		initialize saltv,tempv
c		--------------------------------------------

		if( .not. has_restart(3) ) then	!no restart of T/S values
		  call conini(nlvdi,saltv,salref,sstrat,hdkov)
		  call conini(nlvdi,tempv,temref,tstrat,hdkov)

		  if( ibarcl .eq. 1 .or. ibarcl .eq. 3) then
		    call ts_file_init(it,nlvdim,nlv,nkn,tempv,saltv)
		  else
	            call diagnostic(it,nlvdim,nlv,nkn,tempv,saltv)
		  end if
		end if

c		--------------------------------------------
c		initialize open boundary conditions
c		--------------------------------------------

                nintp = 2
                nvar = 1
                cdef = 0.
		what = 'temp'
		call bnds_init(what,tempn,nintp,nvar,nb3dim,bnd3_temp,cdef)
		call bnds_set_def(what,nb3dim,bnd3_temp)
		what = 'salt'
		call bnds_init(what,saltn,nintp,nvar,nb3dim,bnd3_salt,cdef)
		call bnds_set_def(what,nb3dim,bnd3_salt)

c		initialize rhov, bpresv (we call it twice since
c		rhov depends on bpresv and viceversa
c		-> we iterate to the real solution)

		do k=1,nkn
		  do l=1,nlvdi
		    rhov(l,k) = 0.	!rhov is rho^prime => 0/
		    bpresv(l,k) = 0.
                  end do
		end do

		call rhoset_shell

        	iu = 0
        	itmcon = iround(getpar('itmcon'))
        	idtcon = iround(getpar('idtcon'))
        	call confop(iu,itmcon,idtcon,nlv,2,'nos')

		if( binitial_nos ) then
		  call confil(iu,itmcon,idtcon,11,nlvdi,saltv)
		  call confil(iu,itmcon,idtcon,12,nlvdi,tempv)
		end if

                call getinfo(ninfo)

	end if

c normal call

	icall=icall+1

	if(mode.eq.0) return

	cw = 3991.
	row = 1026.
	wsink = 0.
	ivar = 1
	t = it

	shpar=getpar('shpar')   !out of initialization because changed
	thpar=getpar('thpar')

	if( ibarcl .eq. 2 ) then
	  call diagnostic(it,nlvdim,nlv,nkn,tempv,saltv)
	end if

c salt and temperature transport & diffusion

	if( badvect ) then

!$OMP PARALLEL PRIVATE(tid)
!$OMP SECTIONS
!$OMP SECTION

c	  tid = OMP_GET_THREAD_NUM()
c	  write(6,*) 'number of thread of temp: ',tid

          if( itemp .gt. 0 ) then
		!call check_layers('temp before bnd',tempv)
		call scal_bnd('temp',t,bnd3_temp)
		!call check_layers('temp after bnd',tempv)
                call scal_adv('temp',0
     +                          ,tempv,bnd3_temp
     +                          ,thpar,wsink
     +                          ,difhv,difv,difmol)
		!call check_layers('temp after adv',tempv)
	  end if

!$OMP SECTION

c	  tid = OMP_GET_THREAD_NUM()
c	  write(6,*) 'number of thread of salt: ',tid

          if( isalt .gt. 0 ) then
		!call check_layers('salt before bnd',saltv)
		call scal_bnd('salt',t,bnd3_salt)
		!call check_layers('salt after bnd',saltv)
                call scal_adv('salt',0
     +                          ,saltv,bnd3_salt
     +                          ,shpar,wsink
     +                          ,difhv,difv,difmol)
		!call check_layers('salt after adv',saltv)
          end if

!$OMP END SECTIONS NOWAIT
!$OMP END PARALLEL

	end if

c----------- end of debug section -------------

c compute total mass

	if( binfo ) then
	  call tsmass(saltv,+1,nlvdim,stot) 
	  call tsmass(tempv,+1,nlvdim,ttot) 
	  write(ninfo,*) 'total_mass_T/S: ',it,ttot,stot

          call conmima(nlvdi,saltv,smin,smax)
          call conmima(nlvdi,tempv,tmin,tmax)
          write(ninfo,2020) 'tsmima: ',it,tmin,tmax,smin,smax
 2020	  format(a,i10,4f8.2)
	end if

c heat flux through surface

        call get_timestep(dt)
	call qflux3d(it,dt,nkn,nlvdim,tempv,dq)

c compute rhov and bpresv

	call rhoset_shell

c min/max

	call stmima(saltv,nkn,nlvdi,ilhkv,smin,smax)
	call stmima(tempv,nkn,nlvdi,ilhkv,tmin,tmax)
	call stmima(rhov,nkn,nlvdi,ilhkv,rmin,rmax)

c print of results

	call confil(iu,itmcon,idtcon,11,nlvdi,saltv)
	call confil(iu,itmcon,idtcon,12,nlvdi,tempv)

	return
   99	continue
	write(6,*) 'Value of ibarcl not allowed: ',ibarcl
	stop 'error stop barocl: ibarcl'
	end

c********************************************************

	subroutine rhoset_shell

c sets rho iterating to real solution

	implicit none

	logical biter
	integer itermax,iter
	real eps,resid,resid_old

	itermax = 10
	eps = 1.e-7

	biter = .true.
	iter = 0
	resid = 0.
	resid_old = 0.

	do while( biter )
	  resid_old = resid
          call rhoset(resid)
	  !write(6,*) 'GGGGGGGGGGGGGGGGU: ',resid,resid_old
	  iter = iter + 1
	  if( resid .lt. eps ) biter = .false.
	  if( abs(resid-resid_old) .lt. eps ) biter = .false.
	  if( iter .gt. itermax ) biter = .false.
	end do

	if( iter .gt. itermax ) then
	  write(6,*) '*** warning: max iterations in rhoset_shell ',resid
	end if

	end

c********************************************************

	subroutine rhoset(resid)

c computes rhov and bpresv
c
c 1 bar = 100 kPascal ==> factor 1.e-5
c pres = rho0*g*(zeta-z) + bpresv
c with bpresv = int_{z}^{zeta}(g*rho_prime)dz
c and rho_prime = rho - rho_0 = sigma - sigma_0
c
c in bpresv() is bpresv as defined above
c in rhov()   is rho_prime (=sigma_prime)
c
c brespv() and rhov() are given at node and layer interface

	implicit none

	real resid
c parameter
	include 'param.h'
c common

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft
        integer nlvdi,nlv
        common /level/ nlvdi,nlv
	real saltv(nlvdim,1),tempv(nlvdim,1)
	common /saltv/saltv, /tempv/tempv
	integer ilhkv(1)
	common /ilhkv/ilhkv
	real bpresv(nlvdim,1),rhov(nlvdim,1)
	common /bpresv/bpresv, /rhov/rhov
	real hkv(1)
	common /hkv/hkv
	real hldv(1)
	common /hldv/hldv

        real hdkov(nlvdim,1)
        common /hdkov/hdkov

c local
	logical bdebug,debug,bsigma
	integer k,l,lmax
	real sigma0,rho0,pres
	real depth,hlayer,h,htot
	real rhop,presbt,presbc,dpresc
	double precision dresid
c functions
	real sigma

	rho0 = rowass
	sigma0 = rho0 - 1000.

	debug=.false.
	bdebug=.false.

	bsigma = hldv(1) .lt. 0.

	if(debug) write(6,*) sigma0,rowass,rho0

	dresid = 0.

	do k=1,nkn
	  depth = 0.
	  presbc = 0.
	  lmax=ilhkv(k)
	  htot = hkv(k)

	  do l=1,lmax
	    h = hdkov(l,k)
	    h = hldv(l)
	    if( bsigma ) h = -h * htot
	    hlayer = 0.5 * h
	    depth = depth + hlayer
	    rhop = rhov(l,k)			!rho^prime

	    dpresc = rhop * grav * hlayer	!differential bc. pres.
	    presbc = presbc + dpresc            !baroclinic pres. (mid-layer)
	    presbt = rho0 * grav * depth	!barotropic pressure

	    pres = 1.e-5 * ( presbt + presbc )	!pressure in bars (BUG)
	
	    rhop = sigma(saltv(l,k),tempv(l,k),pres) - sigma0

	    dresid = dresid + (rhov(l,k)-rhop)**2

	    rhov(l,k) = rhop
	    bpresv(l,k) = presbc

	    depth = depth + hlayer
	    presbc = presbc + dpresc		!baroclinic pres. (bottom-lay.)
	  end do
	end do

	resid = dresid

	return
	end

c******************************************************************

	subroutine getts(l,k,t,s)

c accessor routine to get T/S

        implicit none

        integer k,l
        real t,s

	include 'param.h'

	real saltv(nlvdim,1),tempv(nlvdim,1),rhov(nlvdim,1)
	common /saltv/saltv, /tempv/tempv, /rhov/rhov

        t = tempv(l,k)
        s = saltv(l,k)

        end

c*******************************************************************	

	subroutine diagnostic(it,nlvdim,nlv,nkn,tempv,saltv)

	implicit none

	integer it
	integer nlvdim
	integer nlv
	integer nkn
	real tempv(nlvdim,1)
	real saltv(nlvdim,1)

	character*80 tempf,saltf

	integer itt,its
	integer iutemp,iusalt,itnext,idtnext
	save iutemp,iusalt,itnext,idtnext

	integer ifileo

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
	  !tempf = 'dati_interp.temp'
	  !saltf = 'dati_interp.sal'
	  tempf = 'temp_dia.dat'
	  saltf = 'sal_dia.dat'

	  call ts_file_open(tempf,iutemp)
	  call ts_file_open(saltf,iusalt)

	  call read_next_record(itt,iutemp,nkn,nlvdim,nlv,tempv)
	  call read_next_record(its,iusalt,nkn,nlvdim,nlv,saltv)
	  write(6,*) 'diagnostic: first record read ',it

          itnext = 0

	  if( itt .ne. itnext .and. its .ne. itnext ) then
	    write(6,*) it,itt,its,itnext
	    stop 'error stop diagnostic (1): it'
	  end if

	  idtnext = 86400
	  itnext = itnext + idtnext

	  icall = 1
	end if

	if( it .lt. itnext ) return

	call read_next_record(itt,iutemp,nkn,nlvdim,nlv,tempv)
	call read_next_record(its,iusalt,nkn,nlvdim,nlv,saltv)
	write(6,*) 'diagnostic: new record read ',it

	if( itt .ne. itnext .and. its .ne. itnext ) then
	  write(6,*) it,itt,its,itnext
	  stop 'error stop diagnostic (2): it'
	end if

	itnext = itnext + idtnext

	end

c*******************************************************************	

	subroutine read_next_record(it,iunit,nkn,nlvdim,nlv,value)

	implicit none

	integer it
	integer iunit
	integer nkn
	integer nlvdim
	integer nlv
	real value(nlvdim,1)

	real hlv(1)
	common /hlv/hlv

	logical bformat
	integer nknaux,lmax,nvar
	integer i,l
	real val

	bformat = .true.

	if( iunit .le. 0 ) return

	if( bformat ) then
	  read(iunit,*) it,nknaux,lmax,nvar
	  read(iunit,*) (hlv(l),l=1,lmax)				!DEB
	else
	  read(iunit) it,nknaux,lmax,nvar
	  !read(iunit) (hlv(l),l=1,lmax)
	end if

	write(6,*)'reading initial T/S values', it,nknaux,lmax,nvar	!DEB

	if( nkn .ne. nknaux ) stop 'error stop read_next_record: nkn'
	if( nvar .ne. 1 ) stop 'error stop read_next_record: nvar'
	if( lmax .gt. nlvdim ) stop 'error stop read_next_record: nlvdim'

	if( bformat ) then
	  read(iunit,*) ((value(l,i),l=1,lmax),i=1,nkn)
	else
	  read(iunit) ((value(l,i),l=1,lmax),i=1,nkn)
	end if

	if( nlv .gt. lmax ) then
	  do i=1,nkn
	    val = value(lmax,i)
	    do l=lmax+1,nlv
	      value(l,i) = val
	    end do
	  end do
	end if

	end

c*******************************************************************	

	subroutine ts_file_init(it,nlvdim,nlv,nkn,tempv,saltv)

c initialization of T/S from file

	implicit none

        integer it
        integer nlvdim
        integer nlv
        integer nkn
        real tempv(nlvdim,1)
        real saltv(nlvdim,1)

        character*80 tempf,saltf

        integer itt,its
        integer iutemp,iusalt

	call getfnm('tempin',tempf)
	call getfnm('saltin',saltf)

	if( tempf .ne. ' ' ) then
	  call ts_file_open(tempf,iutemp)
          call read_next_record(itt,iutemp,nkn,nlvdim,nlv,tempv)
	  close(iutemp)
          write(6,*) 'temperature initialized from file ',tempf
	end if

	if( saltf .ne. ' ' ) then
	  call ts_file_open(saltf,iusalt)
          call read_next_record(its,iusalt,nkn,nlvdim,nlv,saltv)
	  close(iusalt)
          write(6,*) 'salinity initialized from file ',saltf
	end if

	end

c*******************************************************************	

	subroutine ts_file_open(name,iunit)

c opens T/S file

	implicit none

	character*(*) name
	integer iunit

	logical bformat
	integer ifileo

	bformat = .true.

	if( bformat ) then
	  iunit = ifileo(0,name,'form','old')
	else
	  iunit = ifileo(0,name,'unform','old')
	end if

	end

c*******************************************************************	

	subroutine check_layers(what,vals)

	implicit none

	include 'param.h'

	character*(*) what
	real vals(nlvdim,nkndim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nlvdi,nlv
        common /level/ nlvdi,nlv
        integer ilhkv(1)
        common /ilhkv/ilhkv

	integer l,k,lmax
	real valmin,valmax

	write(6,*) 'checking layer structure : ',what

            do l=1,nlv
              valmin = +999.
              valmax = -999.
              do k=1,nkn
                lmax = ilhkv(k)
                if( l .le. lmax ) then
                  valmin = min(valmin,vals(l,k))
                  valmax = max(valmax,vals(l,k))
                end if
              end do
              write(6,*) l,valmin,valmax
            end do

	end

c*******************************************************************	

