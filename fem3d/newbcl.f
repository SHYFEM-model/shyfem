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
c 10.02.2012    ggu     bug in call to ts_next_record (called with nlvddi)
c 23.02.2012    ccf     do noy check depth structure
c 09.03.2012    deb     bug fix in ts_next_record: ilhkv was real
c 31.10.2012    ggu     open and next_record transfered to subtsuvfile.f
c 05.09.2013    ggu     limit salinity to [0,...]
c 25.03.2014    ggu     new offline
c 10.07.2014    ggu     only new file format allowed
c 20.10.2014    ggu     pass ids to scal_adv()
c 10.02.2015    ggu     call to bnds_read_new() introduced
c 15.10.2015    ggu     added new calls for shy file format
c 26.10.2015    ggu     bug fix for parallel code (what was not set)
c 10.06.2016    ggu     not used routines deleted
c 14.06.2016    ggu     open and write of file in own subroutine
c 27.06.2016    ggu     bug fix: irho was not saved
c
c*****************************************************************

	subroutine barocl(mode)

c amministrates the baroclinic time step
c
c mode : =0 initialize  >0 normal call
c
c written 09.01.94 by ggu  (from scratch)
c
	use mod_layer_thickness
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use mod_hydro
	use levels
	use basin

	implicit none
c
c arguments
	integer mode
c common
	include 'femtime.h'
	include 'pkonst.h'
	include 'mkonst.h'

c local
	logical debug
	logical bgdebug
        logical binfo
        logical bstop
	logical binitial_nos
	logical boff
	integer levdbg
	integer ie
	integer idtext,itmext
	integer imin,imax
	integer nintp,nvar
	integer nbc
	integer id
	real cdef(1)
	real xmin,xmax
	real salref,temref,sstrat,tstrat
	real shpar,thpar
        real s
	real dt
        real gamma,gammax
	real mass
	real wsink
	real robs
	double precision dtime0,dtime
	integer isact,l,k,lmax
	integer kspec
	integer icrst
	integer ftype
	real stot,ttot,smin,smax,tmin,tmax,rmin,rmax
	double precision v1,v2,mm
	character*80 file
	character*4 what
c functions
c	real sigma
	real getpar
	double precision scalcont,dq
	integer iround
	integer nbnds
	logical rst_use_restart
	logical has_output,next_output
	logical has_output_d,next_output_d

	integer tid
	!integer openmp_get_thread_num
	
	double precision theatold,theatnew
	double precision theatconv1,theatconv2,theatqfl1,theatqfl2
c save
	logical, save :: badvect,bobs
        integer, save :: ninfo = 0
	integer, save :: ibarcl
        integer, save :: itemp,isalt,irho
        integer, save :: ia_out(4)
	real, save :: difmol
        double precision, save :: da_out(4)
	integer, save, allocatable :: idtemp(:),idsalt(:)
	integer, save :: icall = 0

c----------------------------------------------------------
c parameter setup and check
c----------------------------------------------------------

	if(icall.eq.-1) return

        call is_offline(2,boff)
        !if( boff ) write(6,*) 'TS reading from offline...'
        if( boff ) return

	levdbg = nint(getpar('levdbg'))
	debug = levdbg .ge. 3
	binfo = .true.
        bgdebug = .false.
	binitial_nos = .true.

	dtime = t_act

c----------------------------------------------------------
c initialization
c----------------------------------------------------------

	if(icall.eq.0) then	!first time

		ibarcl=iround(getpar('ibarcl'))
		if(ibarcl.le.0) icall = -1
		if(ibarcl.gt.4) then
		  write(6,*) 'Value of ibarcl not allowed: ',ibarcl
	          stop 'error stop barocl: ibarcl'
		end if
		if(icall.eq.-1) return

		badvect = ibarcl .ne. 2
		bobs = ibarcl .eq. 4

		salref=getpar('salref')
		temref=getpar('temref')
		sstrat=getpar('sstrat')
		tstrat=getpar('tstrat')
		difmol=getpar('difmol')
                itemp=iround(getpar('itemp'))
                isalt=iround(getpar('isalt'))
                irho=iround(getpar('irho'))

c		--------------------------------------------
c		initialize saltv,tempv
c		--------------------------------------------

		if( .not. rst_use_restart(3) ) then	!no restart of T/S values
		  call conini(nlvdi,saltv,salref,sstrat,hdkov)
		  call conini(nlvdi,tempv,temref,tstrat,hdkov)

		  if( ibarcl .eq. 1 .or. ibarcl .eq. 3) then
		    call ts_init(itanf,nlvdi,nlv,nkn,tempv,saltv)
		  else if( ibarcl .eq. 2 ) then
		    call ts_diag(itanf,nlvdi,nlv,nkn,tempv,saltv)
		  else if( ibarcl .eq. 4 ) then		!interpolate to T/S
	  	    call ts_nudge(itanf,nlvdi,nlv,nkn,tempv,saltv)
		  else
	            stop 'error stop barocl: internal error (1)'
		  end if
		end if

c		--------------------------------------------
c		initialize observations and relaxation times
c		--------------------------------------------

		tobsv = 0.
		sobsv = 0.
		ttauv = 0.
		stauv = 0.

c		--------------------------------------------
c		initialize open boundary conditions
c		--------------------------------------------

		nbc = nbnds()
		allocate(idtemp(nbc))
		allocate(idsalt(nbc))
		idtemp = 0
		idsalt = 0

		dtime0 = itanf
                nintp = 2
                nvar = 1
                cdef(1) = 0.
		what = 'temp'
		call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +					,cdef,idtemp)
		what = 'salt'
		call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv
     +					,cdef,idsalt)

c		--------------------------------------------
c		initialize rhov, bpresv
c		--------------------------------------------

c		rhov depends on bpresv and viceversa
c		-> we iterate to the real solution

		rhov = 0.		!rhov is rho^prime => 0
		bpresv = 0.

		call rhoset_shell

c		--------------------------------------------
c		initialize output files
c		--------------------------------------------

		call bcl_open_output(ia_out,da_out,itemp,isalt,irho)
		call bcl_write_output(dtime,ia_out,da_out,itemp,isalt,irho)

                call getinfo(ninfo)

	end if

	icall=icall+1

	if(mode.eq.0) return

c----------------------------------------------------------
c normal call
c----------------------------------------------------------

	wsink = 0.
	robs = 0.
	if( bobs ) robs = 1.

	shpar=getpar('shpar')   !out of initialization because changed
	thpar=getpar('thpar')

	if( ibarcl .eq. 2 ) then
	  call ts_diag(it,nlvdi,nlv,nkn,tempv,saltv)
	else if( ibarcl .eq. 4 ) then
	  call ts_nudge(it,nlvdi,nlv,nkn,tobsv,sobsv)
	end if

c----------------------------------------------------------
c salt and temperature transport and diffusion
c----------------------------------------------------------

	if( badvect ) then

	  call openmp_get_thread_num(tid)
	  !write(6,*) 'number of thread of temp: ',tid

          if( itemp .gt. 0 ) then
		what = 'temp'
	        call bnds_read_new(what,idtemp,dtime)
	  end if
          if( isalt .gt. 0 ) then
		what = 'salt'
	        call bnds_read_new(what,idsalt,dtime)
	  end if
	  
!$OMP TASK PRIVATE(what,dtime) FIRSTPRIVATE(thpar,wsink,robs,itemp,it) 
!$OMP&     SHARED(idtemp,tempv,difhv,difv,difmol,tobsv,ttauv)
!$OMP&     DEFAULT(NONE)
!$OMP&     IF(itemp > 0)

          if( itemp .gt. 0 ) then
		what = 'temp'
                call scal_adv_nudge(what,0
     +                          ,tempv,idtemp
     +                          ,thpar,wsink
     +                          ,difhv,difv,difmol
     +				,tobsv,robs,ttauv)
	  end if

!$OMP END TASK

!	  call openmp_get_thread_num(tid)
!	  !write(6,*) 'number of thread of salt: ',tid

!$OMP TASK PRIVATE(what,dtime) FIRSTPRIVATE(shpar,wsink,robs,isalt,it) 
!$OMP&     SHARED(idsalt,saltv,difhv,difv,difmol,sobsv,stauv)
!$OMP&     DEFAULT(NONE)
!$OMP&     IF(isalt > 0)

          if( isalt .gt. 0 ) then
		what = 'salt'
                call scal_adv_nudge(what,0
     +                          ,saltv,idsalt
     +                          ,shpar,wsink
     +                          ,difhv,difv,difmol
     +				,sobsv,robs,stauv)
          end if

!$OMP END TASK
!$OMP TASKWAIT

	end if

	if( binfo ) then
          if( itemp .gt. 0 ) then
  	    call tsmass(tempv,+1,nlvdi,ttot) 
      	    call conmima(nlvdi,tempv,tmin,tmax)
!$OMP CRITICAL
  	    write(ninfo,*) 'temp: ',it,ttot,tmin,tmax
!$OMP END CRITICAL
	  end if
          if( isalt .gt. 0 ) then
  	    call tsmass(saltv,+1,nlvdi,stot) 
       	    call conmima(nlvdi,saltv,smin,smax)
!$OMP CRITICAL
  	    write(ninfo,*) 'salt: ',it,stot,smin,smax
!$OMP END CRITICAL
	  end if
	end if

c----------------------------------------------------------
c heat flux through surface	!ccf --> moved to meteo_force
c----------------------------------------------------------

	call compute_heat_flux

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
c write results to file
c----------------------------------------------------------

	call bcl_write_output(dtime,ia_out,da_out,itemp,isalt,irho)

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

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

	use mod_layer_thickness
	use mod_ts
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real resid
c common
	include 'femtime.h'
	include 'pkonst.h'

c local
	logical bdebug,debug,bsigma
	integer k,l,lmax
	integer nresid,nsigma
	real sigma0,rho0,pres,hsigma
	real depth,hlayer,hh
	real rhop,presbt,presbc,dpresc
	real salt
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
	
	    salt = max(0.,saltv(l,k))
	    rhop = sigma(salt,tempv(l,k),pres) - sigma0
	    call set_rhomud(k,l,rhop)

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

	use mod_ts
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	real smin,smax,tmin,tmax,rmin,rmax
	character*30 text

	text = '*** tsrho_check'

	call stmima(saltv,nkn,nlvdi,ilhkv,smin,smax)
	call stmima(tempv,nkn,nlvdi,ilhkv,tmin,tmax)
	call stmima(rhov,nkn,nlvdi,ilhkv,rmin,rmax)

	write(6,*) 'S   min/max: ',smin,smax
	write(6,*) 'T   min/max: ',tmin,tmax
	write(6,*) 'Rho min/max: ',rmin,rmax

	write(6,*) 'checking for Nans...'
        call check2Dr(nlvdi,nlv,nkn,saltv,-1.,+70.,text,'saltv')
        call check2Dr(nlvdi,nlv,nkn,tempv,-30.,+70.,text,'tempv')
        call check2Dr(nlvdi,nlv,nkn,rhov,-2000.,+2000.,text,'rhov')

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine ts_diag(it,nlvddi,nlv,nkn,tempv,saltv)

	implicit none

	integer it
	integer nlvddi
	integer nlv
	integer nkn
	real tempv(nlvddi,nkn)
	real saltv(nlvddi,nkn)

	character*80 tempf,saltf
	integer iutemp(3),iusalt(3)
	save iutemp,iusalt
	real getpar

	integer icall
	data icall /0/
	save icall

	tempf = 'temp_diag.dat'
	saltf = 'salt_diag.dat'

	if( icall .eq. 0 ) then
	  call ts_file_open(tempf,it,nkn,nlv,iutemp)
	  call ts_file_open(saltf,it,nkn,nlv,iusalt)
	  call ts_file_descrp(iutemp,'temp diag')
	  call ts_file_descrp(iusalt,'salt diag')
	  icall = 1
	end if

        call ts_next_record(it,iutemp,nlvddi,nkn,nlv,tempv)
        call ts_next_record(it,iusalt,nlvddi,nkn,nlv,saltv)

	end

c*******************************************************************	

	subroutine ts_nudge(it,nlvddi,nlv,nkn,tobsv,sobsv)

	implicit none

	integer it
	integer nlvddi
	integer nlv
	integer nkn
	real tobsv(nlvddi,nkn)
	real sobsv(nlvddi,nkn)

	character*80 tempf,saltf
	integer iutemp(3),iusalt(3)
	save iutemp,iusalt
	real getpar

	integer icall
	data icall /0/
	save icall

	tempf = 'temp_obs.dat'
	saltf = 'salt_obs.dat'

	if( icall .eq. 0 ) then
	  call ts_file_open(tempf,it,nkn,nlv,iutemp)
	  call ts_file_open(saltf,it,nkn,nlv,iusalt)
	  call ts_file_descrp(iutemp,'temp nudge')
	  call ts_file_descrp(iusalt,'salt nudge')
	  icall = 1
	end if

        call ts_next_record(it,iutemp,nlvddi,nkn,nlv,tobsv)
        call ts_next_record(it,iusalt,nlvddi,nkn,nlv,sobsv)

	end

c*******************************************************************	

	subroutine ts_init(it0,nlvddi,nlv,nkn,tempv,saltv)

c initialization of T/S from file

	implicit none

        integer it0
        integer nlvddi
        integer nlv
        integer nkn
        real tempv(nlvddi,nkn)
        real saltv(nlvddi,nkn)

        character*80 tempf,saltf

        integer itt,its
        integer iutemp(3),iusalt(3)

	call getfnm('tempin',tempf)
	call getfnm('saltin',saltf)

	if( tempf .ne. ' ' ) then
	  itt = it0
	  write(6,*) 'ts_init: opening file for temperature'
	  call ts_file_open(tempf,itt,nkn,nlv,iutemp)
	  call ts_file_descrp(iutemp,'temp init')
          call ts_next_record(itt,iutemp,nlvddi,nkn,nlv,tempv)
	  call ts_file_close(iutemp)
          write(6,*) 'temperature initialized from file ',tempf
	end if

	if( saltf .ne. ' ' ) then
	  its = it0
	  write(6,*) 'ts_init: opening file for salinity'
	  call ts_file_open(saltf,its,nkn,nlv,iusalt)
	  call ts_file_descrp(iusalt,'salt init')
          call ts_next_record(its,iusalt,nlvddi,nkn,nlv,saltv)
	  call ts_file_close(iusalt)
          write(6,*) 'salinity initialized from file ',saltf
	end if

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine bcl_open_output(ia_out,da_out,itemp,isalt,irho)

c opens output of T/S

	use levels

	implicit none

	integer ia_out(4)
	double precision da_out(4)
	integer itemp,isalt,irho

	integer nvar,id,ishyff
	logical has_output
	logical has_output_d
	real getpar

	ishyff = nint(getpar('ishyff'))

	nvar = 0
	if( itemp .gt. 0 ) nvar = nvar + 1
	if( isalt .gt. 0 ) nvar = nvar + 1
	if( irho  .gt. 0 ) nvar = nvar + 1

	call init_output('itmcon','idtcon',ia_out)
	if( ishyff == 1 ) ia_out = 0

	if( has_output(ia_out) ) then
	  call open_scalar_file(ia_out,nlv,nvar,'nos')
	end if

	call init_output_d('itmcon','idtcon',da_out)
	if( ishyff == 0 ) da_out = 0

	if( has_output_d(da_out) ) then
	  call shyfem_init_scalar_file('ts',nvar,.false.,id)
	  da_out(4) = id
	end if

	end

c*******************************************************************	

	subroutine bcl_write_output(dtime,ia_out,da_out,itemp,isalt,irho)

c writes output of T/S

	use levels
	use mod_ts

	implicit none

	double precision dtime
	integer ia_out(4)
	double precision da_out(4)
	integer itemp,isalt,irho

	integer id
	logical next_output
	logical next_output_d

	if( next_output(ia_out) ) then
	  if( isalt .gt. 0 ) then
	    call write_scalar_file(ia_out,11,nlvdi,saltv)
	  end if
	  if( itemp .gt. 0 ) then
	    call write_scalar_file(ia_out,12,nlvdi,tempv)
	  end if
	  if( irho  .gt. 0 ) then
	    call write_scalar_file(ia_out,13,nlvdi,rhov)
	  end if
	end if

	if( next_output_d(da_out) ) then
	  id = nint(da_out(4))
	  if( isalt .gt. 0 ) then
	    call shy_write_scalar_record(id,dtime,11,nlvdi,saltv)
	  end if
	  if( itemp .gt. 0 ) then
	    call shy_write_scalar_record(id,dtime,12,nlvdi,tempv)
	  end if
	  if( irho  .gt. 0 ) then
	    call shy_write_scalar_record(id,dtime,13,nlvdi,rhov)
	  end if
	end if

	end

c*******************************************************************	
c*******************************************************************	
c*******************************************************************	

	subroutine getts(l,k,t,s)

c accessor routine to get T/S

	use mod_ts

        implicit none

        integer k,l
        real t,s

        t = tempv(l,k)
        s = saltv(l,k)

        end

c******************************************************************
c*******************************************************************	
c*******************************************************************	

