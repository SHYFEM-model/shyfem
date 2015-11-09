c
c $Id: subwrt.f,v 1.58 2010-03-08 17:46:45 georg Exp $
c
c routines dealing with water renewal time
c
c contents :
c
c revision log :
c
c 24.10.2011	ggu	new file copied from subcus.f (jamal)
c 28.02.2012	ggu&deb	completely restructured
c 16.03.2012	ggu	use idtwrt=-1 for no renewal computation
c 10.05.2014	ccf	parameters from the str file
c 31.03.2015	ggu	compute res time for different areas
c 20.05.2015	ggu	rinside computed only once, bug fix for conz==0
c 05.06.2015	ggu	new routine to limit concentration between 0 and c0
c
c******************************************************************
c Parameters to be set in section $wrt of the parameter input file
c
c bnoret	true if no return flow is used (concentrations outside
c		are explicitly set to 0)
c bstir		simulates completely stirred tank
c		(replaces at every time step conz with average conz)
c blog          use logarithmic regression to compute renewal time
c badj          adjust renewal time for tail of distribution
c
c percmin	percentage to reach -> after this stop computation
c		use 0 if no premature end is desired
c iaout		area code of elements out of lagoon (used for init and retflow)
c		use -1 to if no outside areas exist
c c0		initial concentration of tracer
c
c itmin		time from when to compute renewal time (-1 for start of sim)
c itmax		time up to when to compute renewal time (-1 for end of sim)
c idtwrt	time step to reset concentration to c0
c		use 0 if no reset is desired
c		use -1 if no renewal time computation is desired
c
c ctop          maximum to be used for frequency curve
c ccut          cut renewal time at this level (for res time computation)
c------------------------------------------------------------
c
c computation through logarithm:
c
c c = c0*exp(-at)  with a = 1/tau and tau is the residence time
c log(c/c0) = -at  and  log(c0/c) = at
c define y=log(c0/c)=-log(c/c0) then y = at
c minimum squares: f(a) = sum_i (y(i)-at(i))**2
c df/da = 0 => sum_i (y(i)-at(i)) = 0 => a = sum_i(y(i))/sum_i(t(i)) 
c and finally: tau = sum_i(t(i))/sum_i(y(i)) = sum_i(t(i))/sum_i(log(c0/c(i)))
c
c------------------------------------------------------------

!==================================================================
        module mod_renewal_time
!==================================================================

	implicit none

	! other saved values can be integrated here

	!real cvres3(nlvdi,nkn)      	!computed RT 3D
	!real vol(nlvdi,nkn)      		!volume
	!double precision cvacu(nlvdi,nkn)   !conz integrated
	!double precision volacu(nlvdi,nkn)  !volume integrated

	real, save, allocatable :: cvres3(:,:)
	real, save, allocatable :: vol(:,:)
	double precision, save, allocatable :: cvacu(:,:)
	double precision, save, allocatable :: volacu(:,:)

	real, save, allocatable :: rinside(:)

!==================================================================
        contains
!==================================================================

	subroutine mod_renewal_time_init(nkn,nlv)

	integer nkn,nlv

        allocate(cvres3(nlv,nkn))
        allocate(vol(nlv,nkn))
        allocate(cvacu(nlv,nkn))
        allocate(volacu(nlv,nkn))

        allocate(rinside(nkn))

	end subroutine mod_renewal_time_init

!==================================================================
        end module mod_renewal_time
!==================================================================

        subroutine renewal_time

	use mod_conz, only : cnv
	use mod_depth
	use levels
	use basin
	use mod_renewal_time

        implicit none

        include 'param.h'
        include 'femtime.h'

	include 'simul.h'


	!include 'aux_array.h'

	integer ndim
	parameter (ndim=100)

	logical breset,bcompute,binit,belab
        logical bnoret,bstir
	logical blog,badj
	logical bdebug
	logical bresarea,blimit

	integer iaout,itmin,itmax,idtwrt
	integer iret,istir,iadj,ilog
	integer iconz
	real c0
	real ctop,ccut
	real percmin

        real conz,perc,rcorrect
        double precision mass,volume
	
	character*80 title
	integer nvers	

	integer ifemop

	integer k,nin,nvar,ie
	integer ius,iuf,iua
	save ius,iuf,iua
	integer nrepl
	save nrepl
	integer it0
	save it0
	double precision tacu
	save tacu
        double precision mass0,vol0,conz0
        save mass0,vol0,conz0
	integer ia_out(4)
	save ia_out

	integer iadim
	save iadim

        double precision dgetpar

	double precision massa(0:ndim)
	double precision massa0(0:ndim)
	double precision vola(0:ndim)
	double precision vola0(0:ndim)
	double precision conza(0:ndim)
	double precision conza0(0:ndim)
	save massa,massa0,vola,vola0,conza,conza0

        integer icall
        save icall
        data icall / 0 /

	save iaout,idtwrt,itmin,itmax
	save c0,percmin
	save bnoret,bstir,blog,badj
	save ctop,ccut

	bdebug = .true.
	bdebug = .false.

	bresarea = .true.	!compute masses for single areas
	blimit = .true.		!limit concentration between 0 and c0

	if( icall .lt. 0 ) return

c------------------------------------------------------------
c initialization
c------------------------------------------------------------

        if( icall .eq. 0 ) then
          write(6,*) 'Initialization of WRT routine renewal time'

	  call convert_time('idtwrt',idtwrt)

	  call convert_date('itmin',itmin)
	  call convert_date('itmax',itmax)
	  if( itmin .eq. -1 ) itmin = itanf
	  if( itmax .eq. -1 ) itmax = itend

	  if( idtwrt < 0 ) then
	    icall = -1
            write(6,*) 'No renewal time computation'
	    return
	  end if

	  iadim = 0
	  do ie=1,nel
	    iadim = max(iadim,iarv(ie))
	  end do
	  if( .not. bresarea ) iadim = 0
	  if( iadim > ndim ) then
	    write(6,*) 'iadim,ndim: ',iadim,ndim
	    stop 'error stop renewal_time: iadim>ndim'
	  end if

	  iconz = nint(dgetpar('iconz'))
	  if( iconz .ne. 1 ) then
	    write(6,*) 'for renewal time computations the generic'
	    write(6,*) 'concentration must be computed (iconz=1)'
	    stop 'error stop renewal_time: iconz'
	  end if

	  iaout = nint(dgetpar('iaout'))
	  c0 = dgetpar('c0')
	  percmin = dgetpar('percmin')
	  iret = nint(dgetpar('iret'))
	  bnoret = iret.eq.0
	  istir = nint(dgetpar('istir'))
	  bstir = istir.eq.1
	  iadj = nint(dgetpar('iadj'))
	  badj = iadj.eq.1
	  ilog = nint(dgetpar('ilog'))
	  blog = ilog.eq.1
	  ctop = dgetpar('ctop')
	  ccut = dgetpar('ccut')
	  it0 = it

	  call mod_renewal_time_init(nkn,nlvdi)

	  call wrt_flag_inside(rinside,iaout)

	  ius = ifemop('.jas','formatted','new')
	  iuf = ifemop('.frq','formatted','new')
	  iua = ifemop('.jaa','formatted','new')

	  nvar = 1
	  call open_scalar_file(ia_out,nlv,nvar,'wrt')

	  nrepl = -1				!must still initialize
        end if

        icall = icall + 1
	
c------------------------------------------------------------
c is it time to run the routine?
c------------------------------------------------------------

        if( it .lt. itmin ) return
        if( it .gt. itmax ) return

c------------------------------------------------------------
c decide on what to do
c------------------------------------------------------------

c binit		is initial call
c breset	resets concentration
c bcompute	computes residence times and writes to file
c belab		elaborates (accumulates) concentrations

	binit = .false.
	if( nrepl .lt. 0 ) then
	  nrepl = 0
	  binit = .true.
	end if

	breset = binit
	if( idtwrt .gt. 0 ) then
	  if( it-it0 .ge. idtwrt ) breset = .true.
	end if
	if( it .eq. itmax ) breset = .true.	!last time step

	belab = .not. binit
	bcompute = .not. binit

	if( bdebug ) then
	  write(6,*) 'WRT debug:'
	  write(6,*) nrepl,binit,breset,belab,bcompute,bnoret
	  write(6,*) it,ius
	end if

c------------------------------------------------------------
c flag nodes that are inside lagoon (v1v(k)=1)
c------------------------------------------------------------

	!call wrt_flag_inside(v1v,iaout)	!done during initialization

c------------------------------------------------------------
c elaborate results
c------------------------------------------------------------

	conz = 0.
	if( belab ) then
	  if( blimit ) call wrt_limit_conz(c0,cnv)
	  call wrt_massvolconz(cnv,rinside,vol,mass,volume)
	  call wrt_mass_area(iadim,cnv,massa,vola,conza)
	  !call wrt_write_area(iua,it,iadim,massa,massa0)
	  call wrt_write_area(iua,it,iadim,conza,conza0)
	  conz = mass / volume

	  if( bstir ) call wrt_bstir(conz,cnv,rinside)	!stirred tank
          if( bnoret ) call wrt_bnoret(cnv,rinside)	!no return flow

	  tacu = tacu + (it-it0)
	  call acu_acum(blog,it,c0,cnv,vol,rinside,cvacu,volacu)

	  call wrt_restime_summary(ius,it,it0,mass,mass0,rcorrect)
	end if

	if( bdebug ) then
	  write(6,*) 'WRT debug (conz): ',conz
	end if

c------------------------------------------------------------
c reset concentrations
c------------------------------------------------------------

        if( breset ) then		!reset concentrations to c0

       	  write(6,*) 'resetting concentrations for renewal time ',it

c------------------------------------------------------------
c reset variables to compute renewal time (and write to file)
c------------------------------------------------------------

	  if( bcompute ) then	!compute new renewal time
	    rcorrect = 0.	!do not used global correction
	    call acu_comp(ia_out,blog,badj,it,c0,ccut,rcorrect
     +				,tacu,cvacu
     +				,cnv,cvres3)
	    call acu_freq(iuf,it,ctop,rinside,cvres3,volacu)
	    nrepl = nrepl + 1

	    write(6,*) '-------------------------------------------'
	    write(6,*) 'computing res time: ',it,conz,nrepl
	    write(6,*) '-------------------------------------------'
	  end if

	  it0 = it
	  tacu = 0.
	  call acu_reset(cvacu)
	  call acu_reset(volacu)
	  call wrt_breset(c0,cnv,rinside)
	  if( bdebug ) then
	    write(6,*) 'WRT debug:'
	    nin = 0
	    do k=1,nkn
	      if( rinside(k) .ne. 0. ) nin = nin + 1
	    end do
	    write(6,*) 'resetting: ',it,c0,nin
	  end if

	  call wrt_massvolconz(cnv,rinside,vol,mass,volume)
	  mass0 = mass
	  vol0 = volume
	  conz0 = mass0/vol0

	  call wrt_mass_area(iadim,cnv,massa,vola,conza)
	  massa0 = massa
	  vola0 = vola
	  conza0 = conza
	  !call wrt_write_area(iua,it,iadim,massa,massa0)
	  call wrt_write_area(iua,it,iadim,conza,conza0)

	  call wrt_restime_summary(-ius,it,it0,mass,mass0,rcorrect)	!reset
	end if

c------------------------------------------------------------
c finish computation if mass is below threshold
c------------------------------------------------------------

        if( mass0 .ne. 0. ) then
	  perc = mass / mass0
          if( perc .lt. percmin ) then
                stop 'finished computing renewal time'
	  end if
        end if

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine wrt_flag_inside(rinside,iaout)

c flags inside nodes (nodes with area code different from iaout)
c
c on return rinside(k) = 1 for nodes inside domain

	use basin

	implicit none

	include 'param.h'

	real rinside(nkn)
	integer iaout		!area code for outside elements

	integer k,ie,ii,ia

        do k=1,nkn
          rinside(k) = 0.
        end do

        do ie=1,nel
          ia = iarv(ie)
          if( ia .ne. iaout ) then
              do ii=1,3
                k = nen3v(ii,ie)
                rinside(k) = 1.
              end do
          end if
        end do

	end

c*****************************************************************

	subroutine wrt_breset(c0,cnv,rinside)

c resets concentration for start of new computation

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	real c0				!concentration for nodes inside domain
	real cnv(nlvdi,nkn)
	real rinside(nkn)		!flag if node is inside domain

	integer k,l,lmax
	real conz

	do k=1,nkn
          lmax = ilhkv(k)
	  conz = 0.
          if( rinside(k) .ne. 0. ) conz = c0	!internal node
          do l=1,lmax
            cnv(l,k) = conz
          end do
        end do

	end

c*****************************************************************

	subroutine wrt_bstir(c0,cnv,rinside)

c simulates stirred tank

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	real c0				!concentration to impose
	real cnv(nlvdi,nkn)
	real rinside(nkn)		!flag if node is inside domain

	integer k,l,lmax

        do k=1,nkn
          if( rinside(k) .ne. 0. ) then		!internal node
            lmax = ilhkv(k)
            do l=1,lmax
              cnv(l,k) = c0
	    end do
	  end if
	end do
	
	end

c******************************************************

	subroutine wrt_write_area(iua,it,ndim,massa,massa0)

c computes masses for different areas

	implicit none

	integer iua
	integer it
	integer ndim
	double precision massa(0:ndim)
	double precision massa0(0:ndim)

	integer i
	real mass0
	real mass(0:ndim)

	if( ndim .le. 0 ) return

	do i=0,ndim
	  mass(i) = 0.
	  if( massa0(i) > 0. ) then
	    mass(i) = 100. * massa(i) / massa0(i)
	  end if
	end do

	write(iua,*) it,mass
	
	end

c******************************************************

	subroutine wrt_mass_area(ndim,cnv,massa,vola,conza)

c computes masses for different areas

	use evgeom
	use levels
	use basin

	implicit none

	include 'param.h'

	integer ndim
	real cnv(nlvdi,nkn)
	double precision massa(0:ndim)
	double precision vola(0:ndim)
	double precision conza(0:ndim)

	integer ie,k,ii,ia,l,lmax,nlev
	real v,conz,area,hdep
	real h(nlvdi)

	real volnode

	if( ndim .le. 0 ) return

	massa = 0.
	vola = 0.

	do ie=1,nel
	  ia = iarv(ie)
	  if( ia > ndim .or. ia < 0 ) goto 99
          lmax = ilhv(ie)
	  area = 4.*ev(10,ie)
	  call dep3dele(ie,+1,nlev,h)
	  if( lmax .ne. nlev ) goto 98
	  do l=1,lmax
	    hdep = h(l)
	    do ii=1,3
	      k = nen3v(ii,ie)
              v = area * hdep
              conz = cnv(l,k)
              massa(ia) = massa(ia) + v*conz
              vola(ia) = vola(ia) + v
	    end do
	  end do
	end do

	conza = massa / vola

	return
   98	continue
	write(6,*) 'lmax,nlev: ',lmax,nlev
	stop 'error stop wrt_mass_area: internal error (2)'
   99	continue
	write(6,*) 'ie,ia: ',ie,ia
	stop 'error stop wrt_mass_area: internal error (1)'
	end

c******************************************************

	subroutine wrt_limit_conz(c0,cnv)

c limits concentration between 0 and c0

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	real c0
	real cnv(nlvdi,nkn)		!concentration

	integer k,l,lmax

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cnv(l,k) = min(cnv(l,k),c0)
            cnv(l,k) = max(cnv(l,k),0.)
          end do
        end do

	end

c******************************************************

	subroutine wrt_massvolconz(cnv,rinside,vol,mass,volume)

c computes mass and volume on internal nodes

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	real cnv(nlvdi,nkn)		!concentration
	real rinside(nkn)		!flag if node is inside domain
	real vol(nlvdi,nkn)		!volume on node (return)
	double precision mass,volume	!total mass/volume (return)

	integer k,l,lmax
	real v,conz

	real volnode

        mass = 0.
        volume = 0.

        do k=1,nkn
          if( rinside(k) .ne. 0. ) then
            lmax = ilhkv(k)
            do l=1,lmax
              v = volnode(l,k,+1)
              conz = cnv(l,k)
              mass = mass + v*conz
              volume = volume + v
	      vol(l,k) = v
            end do
          end if
        end do

	end

c****************************************************

	subroutine wrt_bnoret(cnv,rinside)

c sets concentration to zero outside of domain

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

        real cnv(nlvdi,nkn)
	real rinside(nkn)		!flag if node is inside domain

	integer k,l,lmax

       	do k=1,nkn
          if( rinside(k) .eq. 0. ) then
            lmax = ilhkv(k)
            do l=1,lmax
              cnv(l,k) = 0.
            end do
          end if
        end do

	end

c***************************************************************

	subroutine acu_reset(cvacu)

c resets acumulated value

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	double precision cvacu(nlvdi,nkn)

	integer k,lmax,l

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cvacu(l,k) = 0.
          end do
        end do

	end

c***************************************************************

	subroutine acu_acum(blog,it,c0,cnv,vol,rinside,cvacu,volacu)

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	logical blog				!use logarithm to compute
	integer it
	real c0					!value used for initialization
	real cnv(nlvdi,nkn)
	real vol(nlvdi,nkn)
	real rinside(nkn)
	double precision cvacu(nlvdi,nkn)
	double precision volacu(nlvdi,nkn)

	logical binside
	integer k,l,lmax
	real dt
	double precision conz,volume
	double precision rl,ddt,cc0

	call get_timestep(dt)

	ddt = dt
	cc0 = c0

        do k=1,nkn
	  binside = rinside(k) > 0.
          lmax = ilhkv(k)
          do l=1,lmax
            conz = cnv(l,k)
	    if( .not. binside ) conz = 0.
	    volume = vol(l,k)
	    if ( conz .gt. cc0 ) conz = cc0
	    if ( conz .lt. 0. ) conz = 0.
	    if( blog ) then
	      if( conz .gt. 0 ) then
	        rl = - log(conz/cc0)
                cvacu(l,k) = cvacu(l,k) + rl
	      end if
	    else
              cvacu(l,k) = cvacu(l,k) + conz*ddt
	    end if
            volacu(l,k) = volacu(l,k) + volume*ddt
          end do
        end do

	!l = 4
	!k = 106
	!write(166,*) it,cnv(l,k),vol(l,k)

	!l = 1
	!k = 100
	!write(6,*) cnv(l,k),cvacu(l,k),volacu(l,k)

	end

c***************************************************************

	subroutine acu_comp(ia_out,blog,badj,it,c0,ccut,rcorrect
     +				,tacu,cvacu
     +				,cnv,cvres3)

c compute renewal time and write to file

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	integer ia_out(4)
	logical blog,badj
	integer it
	real c0
	real ccut
	real rcorrect
	double precision tacu
	double precision cvacu(nlvdi,nkn)		!accumulated conz
	real cnv(nlvdi,nkn)				!last concentration
	real cvres3(nlvdi,nkn)				!computed RT 3D

	integer k,lmax,l,ivar,ierr
	double precision conz,conze,rconv,corr,cc0
	double precision secs_in_day,ttacu

c---------------------------------------------------------------
c set parameters
c---------------------------------------------------------------

	secs_in_day = 86400.

	rconv = 1. / secs_in_day
	ttacu = tacu
	cc0 = c0

c---------------------------------------------------------------
c compute renewal times -> put in cvres3
c---------------------------------------------------------------

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            conz = cvacu(l,k)
	    if( conz .le. 0. ) then
	      !nothing - leave WRT zero
	    else if( blog ) then
	      conz = ttacu / conz
	    else
              conz = conz / cc0			!remnant
	      if( badj ) then
		if( rcorrect .le. 0. ) then
	          conze = cnv(l,k) / cc0
	          if( conze .ge. 1 ) conze = 0.
	          if( conze .le. 0 ) conze = 0.
		  corr = 1. /  ( 1. - conze )
		else
		  corr = rcorrect
		end if
                conz = corr * conz 		!adjusted res time
	      end if
	    end if
	    conz = rconv * conz
	    if ( ccut .gt. 0. .and. conz .gt. ccut ) conz = ccut
            cvres3(l,k) = conz
          end do
        end do

c---------------------------------------------------------------
c write to file
c---------------------------------------------------------------

	ivar = 99
	call write_scalar_file(ia_out,ivar,nlvdi,cvres3)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************

	subroutine wrt_restime_summary(ius,it,it0,mass,mass0,rcorrect)

c perc		percentage of mass still in domain
c restime	renewal time computed by integrating
c restimec	renewal time computed by integrating with correction
c restimel	renewal time computed by fitting regression curve
c resmed	average of renewal times computed
c resstd	standard deviation of renewal time

     	implicit none
	
	integer ius		!negative if reset and summary write
	integer it,it0
	double precision mass
	double precision mass0
	real rcorrect

	real dt
	real perc
        real remnant,rlast
	real restime,restimel,restimec
        real resmed,resstd

	logical breset
	integer ndata,iu
	double precision remint,remlog,remtim
	double precision rsum,rsumsq
	save ndata
	save remint,remlog,remtim
	save rsum,rsumsq
	data ndata / 0 /

	breset = ius .le. 0
	iu = abs(ius)

	if( breset ) then	!reset
	  ndata = 0
	  remint = 0.
	  remlog = 0.
	  remtim = 0.
	  rsum = 0.
	  rsumsq = 0.
	  !return
	end if

	call get_timestep(dt)

	remnant = 0.
        if( mass0 .gt. 0. ) remnant = mass/mass0
	if( remnant > 1. ) remnant = 1.
        perc = 100.*remnant

	remint = remint + remnant*dt	!integrated remnant function
	restime = remint/86400.		!renewal time in days

	rlast = remnant
	if( rlast .ge. 1. ) rlast = 0.
	rcorrect = 1. / (1.-rlast)
	restimec = rcorrect * restime	!corrected renewal time
	restimec = min(restimec,999999.)

	remlog = remlog - log(remnant)
	remtim = remtim + (it-it0)
	restimel = 0.
	if( remlog .gt. 0. ) restimel = ( remtim / remlog ) / 86400.
	restimel = min(restimel,999999.)

	ndata = ndata + 1
	rsum = rsum + restimel
	rsumsq = rsumsq + restimel*restimel
	resmed = rsum / ndata
	resmed = min(resmed,999999.)
	resstd = sqrt( rsumsq/ndata - resmed*resmed )
	resstd = min(resstd,999999.)

        !write(6,1000) it,perc,restime,restimec,restimel,resmed,resstd
        write(iu,1000) it,perc,restime,restimec,restimel,resmed,resstd

 1000	format(i10,6f10.2)
	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine acu_freq(iu,it,ctop,rinside,cvres3,volacu)

c write histogram

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	include 'param.h'

	integer iu
	integer it
	real ctop			!cut at this value of renewal time
	real rinside(nkn)		!point is intern
	real cvres3(nlvdi,nkn)
	double precision volacu(nlvdi,nkn)

	integer ndim
	parameter (ndim=100)

	logical bdebug
	integer k,lmax,l,i,ic
	integer icount(0:ndim)
	double precision dcount(0:ndim)
	double precision dc,tot,vtot
	real conz,c,amax,cmax
	real v,dw,val

	bdebug = .true.
	bdebug = .false.

c---------------------------------------------------------------
c compute maximum
c---------------------------------------------------------------

	cmax = 0.
        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
	    cmax = max(cmax,cvres3(l,k))
	  end do
	end do

	amax = cmax
	if( ctop .gt. 0. .and. amax .gt. ctop ) amax = ctop
	if( ctop .lt. 0. ) amax = -ctop
	write(iu,*) 'cmax: ',it,cmax,amax

c---------------------------------------------------------------
c compute average
c---------------------------------------------------------------

	tot = 0.
	vtot = 0.
        do k=1,nkn
	  if( rinside(k) .le. 0. ) cycle
          lmax = ilhkv(k)
          do l=1,lmax
	    conz = cvres3(l,k)
            v = volacu(l,k)
	    tot = tot + conz * v
	    vtot = vtot + v
	  end do
	end do
	write(iu,2000) 'aver_by_tot: ',it,tot,vtot,tot/vtot

c---------------------------------------------------------------
c compute frequency curve
c---------------------------------------------------------------

	ic = 0
	dc = 0.
	do i=0,ndim
	  icount(i) = 0
	  dcount(i) = 0.
	end do

        do k=1,nkn
	  if( rinside(k) .le. 0. ) cycle
          lmax = ilhkv(k)
          do l=1,lmax
	    conz = cvres3(l,k)
            v = volacu(l,k)

	    i = nint(ndim*conz/amax)
	    if (i .lt. 0) i = 0
	    if (i .gt. ndim) i = ndim

	    icount(i) = icount(i) + 1
	    ic = ic + 1
	    dcount(i) = dcount(i) + v
	    dc = dc + v
	  end do
	end do

c---------------------------------------------------------------
c write frequency curve to file
c---------------------------------------------------------------

	dw = 1.
	tot = 0.
	vtot = 0.
	!call make_name(it,file,'freq_by_bin_','.his')
	!open(11,file=file,status='unknown',form='formatted')
	write(iu,*) 'freq_by_bin: ',it,ndim+1
	do i=0,ndim
	  c = i*amax/ndim
	  val = 100.*icount(i)/(float(ic)*dw)
	  tot = tot + val*c
	  vtot = vtot + val
	  write(iu,*) i,val,icount(i)
	end do
	!close(11)
	write(iu,2000) 'aver_by_bin: ',it,tot,vtot,tot/vtot

	dw = amax/100.
	tot = 0.
	vtot = 0.
	!call make_name(it,file,'freq_by_res_','.his')
	!open(11,file=file,status='unknown',form='formatted')
	write(iu,*) 'freq_by_res: ',it,ndim+1
	do i=0,ndim
	  c = i*amax/ndim
	  val = 100.*icount(i)/(float(ic)*dw)
	  tot = tot + val*c
	  vtot = vtot + val
	  write(iu,*) c,val,icount(i)
	end do
	!close(11)
	write(iu,2000) 'aver_by_res: ',it,tot,vtot,tot/vtot

	dw = 1.
	tot = 0.
	vtot = 0.
	!call make_name(it,file,'freq_by_vol_','.his')
	!open(11,file=file,status='unknown',form='formatted')
	write(iu,*) 'freq_by_vol: ',it,ndim+1
	do i=0,ndim
	  c = i*amax/ndim
	  val = 100.*dcount(i)/(dc*dw)
	  tot = tot + val*c
	  vtot = vtot + val
	  write(iu,*) c,val,icount(i)
	end do
	!close(11)
	write(iu,2000) 'aver_by_vol: ',it,tot,vtot,tot/vtot

c---------------------------------------------------------------
c write out all data to file (for debug and median)
c---------------------------------------------------------------

	!write(76,*) nkn,it
        !do k=1,nkn
	!  conz = cvres3(1,k)
	!  write(76,*) conz
	!end do

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

 2000	format(a,i12,2e14.6,f14.4)
	end

c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************

