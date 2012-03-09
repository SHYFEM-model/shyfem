c dd: newbcl.f,v 1.37 2010-03-08 17:46:45 georg Exp $
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
c 13.01.2009    ggu&deb changes in reading file in ts_next_record()
c 13.10.2009    ggu     in rhoset bug computing pres
c 13.11.2009    ggu     only initialize T/S if no restart, new rhoset_shell
c 19.01.2010    ggu     different call to has_restart() 
c 16.12.2010    ggu     sigma layers introduced (maybe not finished)
c 26.01.2011    ggu     read in obs for t/s (tobsv,sobsv)
c 28.01.2011    ggu     parameters changed in call to ts_nudge()
c 04.03.2011    ggu     better error message for rhoset_shell
c 31.03.2011    ggu     only write temp/salt if computed
c 04.11.2011    ggu     adapted for hybrid coordinates
c 07.11.2011    ggu     hybrid changed to resemble code in newexpl.f
c 11.11.2011    ggu     restructured ts_next_record() and diagnostic()
c 22.11.2011    ggu     bug fix in ts_file_open() -> bhashl
c 02.12.2011    ggu     adapt ts_file_open() for barotropic version (ihashl)
c 27.01.2012    deb&ggu changes for hybrid in ts_file_open,ts_next_record
c 10.02.2012    ggu     bug in call to ts_next_record (called with nlvdim)
c 23.02.2012    ccf     do noy check depth structure
c 09.03.2012    deb     bug fix in ts_next_record: ilhkv was real
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
	real uprv(nlvdim,1), vprv(nlvdim,1)
	real v1v(1),v2v(1)
	integer ilhkv(1)
	common /saltv/saltv, /tempv/tempv, /rhov/rhov
	common /bpresv/bpresv
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

	real tobsv(nlvdim,1)
	common /tobsv/tobsv
	real sobsv(nlvdim,1)
	common /sobsv/sobsv
	real rtauv(nlvdim,1)
	common /rtauv/rtauv
                      
        character*80 saltn(nbcdim)
        character*80 tempn(nbcdim)
        common /saltn/ saltn
        common /tempn/ tempn

c local
	logical debug
	logical badvect
	logical bobs
	logical bgdebug
        logical binfo
        logical bstop
	logical binitial_nos
	integer levdbg
	integer ie
	integer ibarcl
	integer idtext,itmext
	integer imin,imax
	integer nintp,nvar
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
	real robs
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
c save
        real bnd3_temp(nb3dim,0:nbcdim)
        save bnd3_temp
        real bnd3_salt(nb3dim,0:nbcdim)
        save bnd3_salt

        integer iu,itmcon,idtcon
        save iu,itmcon,idtcon

        integer ninfo
        save ninfo

	save badvect,bobs
	save salref,temref
	save difmol
        save itemp,isalt
	save ibarcl
c data
	integer icall
	save icall
	data icall /0/

c----------------------------------------------------------
c parameter setup and check
c----------------------------------------------------------

	if(nlvdim.ne.nlvdi) stop 'error stop : level dimension in barocl'

	if(icall.eq.-1) return

	levdbg = nint(getpar('levdbg'))
	debug = levdbg .ge. 3
	binfo = .true.
        bgdebug = .false.
	binitial_nos = .true.

c----------------------------------------------------------
c initialization
c----------------------------------------------------------

	if(icall.eq.0) then	!first time

		ibarcl=iround(getpar('ibarcl'))
		if(ibarcl.le.0) icall = -1
		if(ibarcl.gt.4) goto 99
		if(icall.eq.-1) return

		badvect = ibarcl .ne. 2
		bobs = ibarcl .eq. 4

		salref=getpar('salref')
		temref=getpar('temref')
		sstrat=getpar('sstrat')
		tstrat=getpar('tstrat')
		difmol=getpar('difmol')
		idtcon=iround(getpar('idtcon'))
		itmcon=iround(getpar('itmcon'))
                itemp=iround(getpar('itemp'))
                isalt=iround(getpar('isalt'))

c		--------------------------------------------
c		initialize saltv,tempv
c		--------------------------------------------

		if( .not. has_restart(3) ) then	!no restart of T/S values
		  call conini(nlvdi,saltv,salref,sstrat,hdkov)
		  call conini(nlvdi,tempv,temref,tstrat,hdkov)

		  if( ibarcl .eq. 1 .or. ibarcl .eq. 3) then
		    call ts_file_init(it,nlv,nkn,tempv,saltv)
		  else if( ibarcl .eq. 2 ) then
		    call ts_diag(it,nlv,nkn,tempv,saltv)
		  else if( ibarcl .eq. 4 ) then		!interpolate to T/S
	  	    call ts_nudge(it,nlv,nkn,tempv,saltv)
		  else
		    goto 99
		  end if
		end if

c		--------------------------------------------
c		initialize observations and relaxation times
c		--------------------------------------------

		do k=1,nkn
		  do l=1,nlv
		    tobsv(l,k) = 0.
		    sobsv(l,k) = 0.
		    rtauv(l,k) = 0.
		  end do
		end do

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
		nvar = 0
		if( itemp .gt. 0 ) nvar = nvar + 1
		if( isalt .gt. 0 ) nvar = nvar + 1
        	itmcon = iround(getpar('itmcon'))
        	idtcon = iround(getpar('idtcon'))
        	call confop(iu,itmcon,idtcon,nlv,nvar,'nos')

		if( binitial_nos ) then
		  if( isalt .gt. 0 ) then
		    call confil(iu,itmcon,idtcon,11,nlvdi,saltv)
		  end if
		  if( isalt .gt. 0 ) then
		    call confil(iu,itmcon,idtcon,12,nlvdi,tempv)
		  end if
		end if

                call getinfo(ninfo)

	end if

	icall=icall+1

	if(mode.eq.0) return

c----------------------------------------------------------
c normal call
c----------------------------------------------------------

	t = it
	wsink = 0.
	robs = 0.
	if( bobs ) robs = 1.

	shpar=getpar('shpar')   !out of initialization because changed
	thpar=getpar('thpar')

	if( ibarcl .eq. 2 ) then
	  call ts_diag(it,nlv,nkn,tempv,saltv)
	else if( ibarcl .eq. 4 ) then
	  call ts_nudge(it,nlv,nkn,tobsv,sobsv)
	end if

c----------------------------------------------------------
c salt and temperature transport and diffusion
c----------------------------------------------------------

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
                call scal_adv_nudge('temp',0
     +                          ,tempv,bnd3_temp
     +                          ,thpar,wsink
     +                          ,difhv,difv,difmol,tobsv,robs)
		!call check_layers('temp after adv',tempv)
	  end if

!$OMP SECTION

c	  tid = OMP_GET_THREAD_NUM()
c	  write(6,*) 'number of thread of salt: ',tid

          if( isalt .gt. 0 ) then
		!call check_layers('salt before bnd',saltv)
		call scal_bnd('salt',t,bnd3_salt)
		!call check_layers('salt after bnd',saltv)
                call scal_adv_nudge('salt',0
     +                          ,saltv,bnd3_salt
     +                          ,shpar,wsink
     +                          ,difhv,difv,difmol,sobsv,robs)
		!call check_layers('salt after adv',saltv)
          end if

!$OMP END SECTIONS NOWAIT
!$OMP END PARALLEL

	end if

c----------------------------------------------------------
c compute total mass
c----------------------------------------------------------

	if( binfo ) then
	  call tsmass(saltv,+1,nlvdim,stot) 
	  call tsmass(tempv,+1,nlvdim,ttot) 
	  write(ninfo,*) 'total_mass_T/S: ',it,ttot,stot

          call conmima(nlvdi,saltv,smin,smax)
          call conmima(nlvdi,tempv,tmin,tmax)
          write(ninfo,2020) 'tsmima: ',it,tmin,tmax,smin,smax
 2020	  format(a,i10,4f8.2)
	end if

c----------------------------------------------------------
c heat flux through surface
c----------------------------------------------------------

        call get_timestep(dt)
	call qflux3d(it,dt,nkn,nlvdim,tempv,dq)

c----------------------------------------------------------
c compute rhov and bpresv
c----------------------------------------------------------

	call rhoset_shell

c----------------------------------------------------------
c compute min/max
c----------------------------------------------------------

	call stmima(saltv,nkn,nlvdi,ilhkv,smin,smax)
	call stmima(tempv,nkn,nlvdi,ilhkv,tmin,tmax)
	call stmima(rhov,nkn,nlvdi,ilhkv,rmin,rmax)

c----------------------------------------------------------
c print of results
c----------------------------------------------------------

	if( isalt .gt. 0 ) then
	  call confil(iu,itmcon,idtcon,11,nlvdi,saltv)
	end if
	if( itemp .gt. 0 ) then
	  call confil(iu,itmcon,idtcon,12,nlvdi,tempv)
	end if

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	return
   99	continue
	write(6,*) 'Value of ibarcl not allowed: ',ibarcl
	stop 'error stop barocl: ibarcl'
	end

c********************************************************
c********************************************************
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
	  iter = iter + 1
	  if( resid .lt. eps ) biter = .false.
	  if( abs(resid-resid_old) .lt. eps ) biter = .false.
	  if( iter .gt. itermax ) biter = .false.
	end do

	if( iter .gt. itermax ) then
	  write(6,*) '*** warning: max iterations in rhoset_shell ',resid
	  call tsrho_check
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
	real hldv(1)
	common /hldv/hldv

        real hdkov(nlvdim,1)
        common /hdkov/hdkov

c local
	logical bdebug,debug,bsigma
	integer k,l,lmax
	integer nresid,nsigma
	real sigma0,rho0,pres,hsigma
	real depth,hlayer,hh
	real rhop,presbt,presbc,dpresc
	double precision dresid
c functions
	real sigma

	rho0 = rowass
	sigma0 = rho0 - 1000.

	debug=.false.
	bdebug=.false.

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	if(debug) write(6,*) sigma0,rowass,rho0

	nresid = 0
	dresid = 0.

	do k=1,nkn
	  depth = 0.
	  presbc = 0.
	  lmax = ilhkv(k)

	  do l=1,lmax
	    bsigma = l .le. nsigma

	    hlayer = hdkov(l,k)
	    if( .not. bsigma ) hlayer = hldv(l)

	    hh = 0.5 * hlayer
	    depth = depth + hh
	    rhop = rhov(l,k)			!rho^prime

	    dpresc = rhop * grav * hh		!differential bc. pres.
	    presbc = presbc + dpresc            !baroclinic pres. (mid-layer)
	    presbt = rho0 * grav * depth	!barotropic pressure

	    pres = 1.e-5 * ( presbt + presbc )	!pressure in bars (BUG)
	
	    rhop = sigma(saltv(l,k),tempv(l,k),pres) - sigma0

	    nresid = nresid + 1
	    dresid = dresid + (rhov(l,k)-rhop)**2

	    rhov(l,k) = rhop
	    bpresv(l,k) = presbc

	    depth = depth + hh
	    presbc = presbc + dpresc		!baroclinic pres. (bottom-lay.)
	  end do
	end do

	resid = dresid/nresid

	return
	end

c*******************************************************************	

	subroutine tsrho_check

c checks values of t/s/rho

	implicit none

	include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real saltv(nlvdim,1),tempv(nlvdim,1)
	common /saltv/saltv, /tempv/tempv
	real rhov(nlvdim,1)
	common /rhov/rhov
	integer ilhkv(1)
	common /ilhkv/ilhkv
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	real smin,smax,tmin,tmax,rmin,rmax
	character*30 text

	text = '*** tsrho_check'

	call stmima(saltv,nkn,nlvdim,ilhkv,smin,smax)
	call stmima(tempv,nkn,nlvdim,ilhkv,tmin,tmax)
	call stmima(rhov,nkn,nlvdim,ilhkv,rmin,rmax)

	write(6,*) 'S   min/max: ',smin,smax
	write(6,*) 'T   min/max: ',tmin,tmax
	write(6,*) 'Rho min/max: ',rmin,rmax

	write(6,*) 'checking for Nans...'
        call check2Dr(nlvdim,nlv,nkn,saltv,-1.,+70.,text,'saltv')
        call check2Dr(nlvdim,nlv,nkn,tempv,-30.,+70.,text,'tempv')
        call check2Dr(nlvdim,nlv,nkn,rhov,-2000.,+2000.,text,'rhov')

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine ts_diag(it,nlv,nkn,tempv,saltv)

	implicit none

	include 'param.h'

	integer it
	integer nlv
	integer nkn
	real tempv(nlvdim,1)
	real saltv(nlvdim,1)

	character*80 tempf,saltf

	tempf = 'temp_diag.dat'
	saltf = 'salt_diag.dat'

	call ts_intp(it,nlv,nkn,tempv,saltv,tempf,saltf)

	end

c*******************************************************************	

	subroutine ts_nudge(it,nlv,nkn,tobsv,sobsv)

	implicit none

	include 'param.h'

	integer it
	integer nlv
	integer nkn
	real tobsv(nlvdim,1)
	real sobsv(nlvdim,1)

	character*80 tempf,saltf

	tempf = 'temp_obs.dat'
	saltf = 'salt_obs.dat'

	call ts_intp(it,nlv,nkn,tobsv,sobsv,tempf,saltf)

	end

c*******************************************************************	

	subroutine ts_intp(it,nlv,nkn,tobsv,sobsv,tempf,saltf)

	implicit none

	include 'param.h'

	integer it
	integer nlv
	integer nkn
	real tobsv(nlvdim,1)
	real sobsv(nlvdim,1)
	character*80 tempf,saltf

	integer iutemp(3),iusalt(3)
	save iutemp,iusalt
	integer ittold,itsold,ittnew,itsnew
	save ittold,itsold,ittnew,itsnew

	real toldv(nlvdim,nkndim)
	real soldv(nlvdim,nkndim)
	real tnewv(nlvdim,nkndim)
	real snewv(nlvdim,nkndim)
	save toldv,soldv,tnewv,snewv

	logical bdebug
	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. -1 ) return

	bdebug = .true.
	bdebug = .false.

c-------------------------------------------------------------
c initialization (open files etc...)
c-------------------------------------------------------------

	if( icall .eq. 0 ) then
	  call ts_file_open(tempf,iutemp)
	  call ts_file_open(saltf,iusalt)

	  write(6,*) 'ts_intp: initializing T/S'
	  call ts_next_record(ittold,iutemp,nkn,nlv,toldv)
	  call ts_next_record(itsold,iusalt,nkn,nlv,soldv)
	  write(6,*) 'ts_intp: first record read ',ittold,itsold

	  call ts_next_record(ittnew,iutemp,nkn,nlv,tnewv)
	  call ts_next_record(itsnew,iusalt,nkn,nlv,snewv)
	  write(6,*) 'ts_intp: second record read ',ittnew,itsnew

	  if( ittold .ne. itsold ) goto 98
	  if( ittnew .ne. itsnew ) goto 98
	  if( it .lt. ittold ) goto 99

	  icall = 1
	end if

c-------------------------------------------------------------
c read new files if necessary
c-------------------------------------------------------------

	do while( it .gt. ittnew )

	  ittold = ittnew
	  call copy_record(nkn,nlvdim,nlv,toldv,tnewv)
	  itsold = itsnew
	  call copy_record(nkn,nlvdim,nlv,soldv,snewv)

	  call ts_next_record(ittnew,iutemp,nkn,nlv,tnewv)
	  call ts_next_record(itsnew,iusalt,nkn,nlv,snewv)
	  write(6,*) 'ts_intp: new record read ',ittnew,itsnew

	  if( ittnew .ne. itsnew ) goto 98

	end do

c-------------------------------------------------------------
c interpolate to new time step
c-------------------------------------------------------------

	call intp_record(nkn,nlvdim,nlv,ittold,ittnew,it
     +				,toldv,tnewv,tobsv)
	call intp_record(nkn,nlvdim,nlv,itsold,itsnew,it
     +				,soldv,snewv,sobsv)

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	return
   98	continue
	write(6,*) ittold,itsold,ittnew,itsnew
	stop 'error stop ts_intp: mismatch time of temp/salt records'
   99	continue
	write(6,*) it,ittold
	stop 'error stop ts_intp: no file for start of simulation'
	end

c*******************************************************************	

	subroutine intp_record(nkn,nlvdim,nlv,itold,itnew,it
     +				,voldv,vnewv,vintpv)

c interpolates records to actual time

	implicit none

	integer nkn,nlvdim,nlv
	integer itold,itnew,it
	real voldv(nlvdim,1)
	real vnewv(nlvdim,1)
	real vintpv(nlvdim,1)

	integer k,l
	real rt

        rt = (it-itold) / float(itnew-itold)

	do k=1,nkn
	  do l=1,nlv
	    vintpv(l,k) = voldv(l,k) + rt * (vnewv(l,k) - voldv(l,k))
	  end do
	end do

	end

c*******************************************************************	

	subroutine copy_record(nkn,nlvdim,nlv,voldv,vnewv)

c copies new record to old one

	implicit none

	integer nkn,nlvdim,nlv
	real voldv(nlvdim,1)
	real vnewv(nlvdim,1)

	integer k,l

	do k=1,nkn
	  do l=1,nlv
	    voldv(l,k) = vnewv(l,k)
	  end do
	end do

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine ts_file_init(it,nlv,nkn,tempv,saltv)

c initialization of T/S from file

	implicit none

	include 'param.h'

        integer it
        integer nlv
        integer nkn
        real tempv(nlvdim,1)
        real saltv(nlvdim,1)

        character*80 tempf,saltf

        integer itt,its
        integer iutemp(3),iusalt(3)

	call getfnm('tempin',tempf)
	call getfnm('saltin',saltf)

	if( tempf .ne. ' ' ) then
	  call ts_file_open(tempf,iutemp)
          call ts_next_record(itt,iutemp,nkn,nlv,tempv)
	  call ts_file_close(iutemp)
          write(6,*) 'temperature initialized from file ',tempf
	end if

	if( saltf .ne. ' ' ) then
	  call ts_file_open(saltf,iusalt)
          call ts_next_record(its,iusalt,nkn,nlv,saltv)
	  call ts_file_close(iusalt)
          write(6,*) 'salinity initialized from file ',saltf
	end if

	end

c*******************************************************************	

	subroutine ts_file_open(name,info)

c opens T/S file

	implicit none

	include 'param.h'

	character*(*) name
	integer info(3)

	real hlv(1)
	common /hlv/hlv

	logical bsigma,bdebug
	logical bformat,bhashl
	integer ifileo,iunit,ios
	integer iformat,ihashl,inzeta
	integer it,nknaux,lmax,nvar
	integer l
	integer nsigma
	real hsigma
	real hl(nlvdim)

	bdebug = .true.
	bdebug = .false.

c-------------------------------------------------------------
c check if formatted or unformatted
c-------------------------------------------------------------

	iunit = ifileo(0,name,'unform','old')
	if( iunit .le. 0 ) goto 99
	read(iunit,iostat=ios) it,nknaux,lmax,nvar
	close(iunit)

	if( ios .lt. 0 ) goto 98

	iformat = ios
	bformat = iformat .gt. 0

        call get_sigma(nsigma,hsigma)
        bsigma = nsigma .gt. 0

c-------------------------------------------------------------
c check if level structure information is available
c-------------------------------------------------------------

	ios = 0
	if( bformat ) then
	  iunit = ifileo(0,name,'form','old')
	  read(iunit,*) it,nknaux,lmax,nvar
	  if( lmax .gt. nlvdim ) goto 95
	  read(iunit,*) (hl(l),l=1,lmax)
	else
	  iunit = ifileo(0,name,'unform','old')
	  read(iunit) it,nknaux,lmax,nvar
	  if( lmax .gt. nlvdim ) goto 95
	  read(iunit,iostat=ios) (hl(l),l=1,lmax)
	end if
	close(iunit)

	if( ios .ne. 0 ) goto 97

	inzeta = 0
	ihashl = 1
	do l=2,lmax
	  if( hl(l) .gt. hl(l-1) ) then		!in zeta
	    inzeta = +1
	  else if( hl(l) .lt. hl(l-1) ) then	!in sigma
	    if( inzeta .gt. 0 ) ihashl = 0
	  else
	    ihashl = 0
	  end if
	end do
        if( lmax .le. 1 ) ihashl = 0

	!ihashl = 0			!force choice if you are sure (0/1)
	bhashl = ihashl .gt. 0

c-------------------------------------------------------------
c if level structure is available, check if compatible
c-------------------------------------------------------------

	!if( bhashl ) then
	!  do l=1,lmax
	!    if( hl(l) .ne. hlv(l) ) goto 96
	!  end do
	!end if

c-------------------------------------------------------------
c open finally and store info
c-------------------------------------------------------------

	if( bformat ) then
	  iunit = ifileo(0,name,'form','old')
	else
	  iunit = ifileo(0,name,'unform','old')
	end if

	info(1) = iunit
	info(2) = iformat
	info(3) = ihashl

	if( bdebug ) then
	  write(6,*) 'debug in ts_file_open:'
	  write(6,*) 'iunit,iformat,ihashl: ',iunit,iformat,ihashl
	  write(6,*) 'lmax,hl: ',lmax
	  write(6,*) (hl(l),l=1,lmax)
	end if

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	return
   95	continue
	write(6,*) 'too many layers in file: ',name
	write(6,*) 'lmax,nlvdim  : ',lmax,nlvdim
	stop 'error stop ts_file_open: dimension nlvdim'
   96	continue
	write(6,*) 'incompatible levels in file: ',name
	write(6,*) 'lmax,hl  : ',lmax,(hl(l),l=1,lmax)
	write(6,*) 'lmax,hlv : ',lmax,(hlv(l),l=1,lmax)
	stop 'error stop ts_file_open: incompatible levels'
   97	continue
	write(6,*) 'read error in file: ',name
	stop 'error stop ts_file_open: read error'
   98	continue
	write(6,*) 'empty file: ',name
	stop 'error stop ts_file_open: empty file'
   99	continue
	write(6,*) 'file name: ',name
	stop 'error stop ts_file_open: cannot open file'
	end

c*******************************************************************	

	subroutine ts_next_record(it,info,nkn,nlv,value)

	implicit none

	include 'param.h'

	integer it
	integer info(3)
	integer nkn
	integer nlv
	real value(nlvdim,1)

	integer ilhkv(1)
	common /ilhkv/ilhkv
	real hlv(1)
	common /hlv/hlv
	real znv(1)
	common /znv/znv
	real hkv(1)
	common /hkv/hkv

	logical bsigma,bdebug
	logical bformat,bhashl
	integer nknaux,lmax,nvar
	integer i,l,iunit
	integer nsigma,nsigma_aux
	real hsigma,hsigma_aux
	real val
        real vmin,vmax

	integer lmax_fem
	real htot
	real zz_fem,zz_data
	real hlv_aux(nlvdim)
        real hl_fem(0:nlvdim+1)
        real hl_data(0:nlvdim+1)
        real val_fem(nlvdim+1)
        real val_data(nlvdim+1)

	bdebug = .true.
	bdebug = .false.

	iunit   = info(1)
	bformat = info(2) .gt. 0
	bhashl  = info(3) .gt. 0

	if( iunit .le. 0 ) return

	if( bdebug ) then
	  write(6,*) 'debug in ts_next_record:'
	  write(6,*) 'iunit,bformat,bhashl: ',it,iunit,bformat,bhashl
	end if

	if( bformat ) then
	  read(iunit,*) it,nknaux,lmax,nvar
	else
	  read(iunit) it,nknaux,lmax,nvar
	end if

	write(6,*)'reading T/S values', it,nknaux,lmax,nvar

	if( nkn .ne. nknaux ) stop 'error stop ts_next_record: nkn'
	if( nvar .ne. 1 ) stop 'error stop ts_next_record: nvar'
	if( lmax .gt. nlvdim ) stop 'error stop ts_next_record: nlvdim'

	if( bformat ) then
	  if( bhashl ) read(iunit,*) (hlv_aux(l),l=1,lmax)
	  read(iunit,*) ((value(l,i),l=1,lmax),i=1,nkn)
	else
	  if( bhashl ) read(iunit) (hlv_aux(l),l=1,lmax)
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

c--------------------------------------------------------------
c interpolate between different vertical structures 
c--------------------------------------------------------------

c the following still has to be checked

        call get_sigma(nsigma,hsigma)		!from basin
        bsigma = nsigma .gt. 0

	call compute_sigma_info(lmax,hlv_aux,nsigma_aux,hsigma_aux)

        do i = 1,nkn

          lmax_fem = ilhkv(i)

          hl_data(0)= -znv(i)
          hl_fem(0) = -znv(i)

          if(nsigma.gt.0.and.nsigma.ge.lmax_fem)then !sigma
            htot = hkv(i) !giusto
          else !zeta e ibrido
            htot = hlv(lmax_fem) !questo non risolve il problema
          endif

          zz_data = znv(i)
          zz_fem = znv(i)

          call set_hybrid_depth(lmax_fem,zz_fem,htot
     +			,hlv,nsigma,hsigma,hl_fem(1))

          htot = hlv_aux(lmax) !questo non risolve il problema

          call set_hybrid_depth(lmax,zz_data,htot
     +			,hlv_aux,nsigma_aux,hsigma_aux,hl_data(1))

          do l=1,lmax
            val_data(l) = value(l,i)
          end do

          call intp_vert(lmax,hl_data,val_data,lmax_fem,hl_fem,val_fem)

          do l = 1,lmax_fem
            value(l,i) = val_fem(l) !DEB
          end do
        enddo

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

        call conmima(nlvdim,value,vmin,vmax)
        write(6,*) 'min/max: ',vmin,vmax

	end

c*******************************************************************	

	subroutine ts_file_close(info)

c closes T/S file

	implicit none

	integer info(3)

	close(info(1))

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

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

c******************************************************************

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

