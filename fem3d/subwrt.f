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
c 01.02.2016	ggu	implemented custom reset
c 15.04.2016	ggu	new input file for custom reset
c 31.10.2016	ggu	new output format for wrt files
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
	use shympi

        implicit none

	include 'simul.h'

        logical, save :: bnoret,bstir,blog,badj
	logical breset,bcompute,binit,belab
	logical bdebug
	logical bmaster
	logical bresarea,blimit

	integer, save :: iaout
	integer iret,istir,iadj,ilog
	integer iconz
	double precision, save :: dtmin,dtmax,ddtwrt
	double precision, save :: dtnext,dtime0
	double precision :: dtime
	real, save :: c0,percmin
	real, save :: ctop,ccut

        real conz,perc,rcorrect
        double precision mass,volume
	
	character*80 title
	character*20 aline
	integer nvers	

	integer ifemop

	integer k,nin,nvar,ie,id
	integer, save :: ius,iuf,iua
	integer, save :: nrepl
	double precision, save :: tacu
        double precision, save :: mass0,vol0,conz0
	double precision, save :: da_out(4)

	integer, save :: iadim

        double precision dgetpar

	double precision, allocatable, save :: massa(:)
	double precision, allocatable, save :: massa0(:)
	double precision, allocatable, save :: vola(:)
	double precision, allocatable, save :: vola0(:)
	double precision, allocatable, save :: conza(:)
	double precision, allocatable, save :: conza0(:)

        integer, save :: icall = 0

	bdebug = .true.
	bdebug = .false.

	bresarea = .true.	!compute masses for single areas
	blimit = .true.		!limit concentration between 0 and c0

	if( icall .lt. 0 ) return

	bmaster = shympi_is_master()

c------------------------------------------------------------
c initialization
c------------------------------------------------------------

        if( icall .eq. 0 ) then
	  if( bmaster ) then
            write(6,*) 'Initialization of WRT routine renewal time'
	  end if

	  call convert_time_d('idtwrt',ddtwrt)
	  call convert_date_d('itmin',dtmin)
	  call convert_date_d('itmax',dtmax)
	  if( dtmin .eq. -1 ) call get_first_dtime(dtmin)
	  if( dtmax .eq. -1 ) call get_last_dtime(dtmax)

	  if( ddtwrt < 0 ) then
	    icall = -1
            if( bmaster ) write(6,*) 'No renewal time computation'
	    return
	  end if

	  if( shympi_is_parallel() ) then
	    stop 'error stop renewal_time: not ready for mpi'
	  end if

	  iadim = maxval(iarv)
	  if( .not. bresarea ) iadim = 0
	  allocate(massa(0:iadim),massa0(0:iadim))
	  allocate(vola(0:iadim),vola0(0:iadim))
	  allocate(conza(0:iadim),conza0(0:iadim))

	  iconz = nint(dgetpar('iconz'))
	  if( iconz .ne. 1 ) then
	   if( bmaster ) then
	    write(6,*) 'for renewal time computations the generic'
	    write(6,*) 'concentration must be computed (iconz=1)'
	    stop 'error stop renewal_time: iconz'
	   end if
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

	  dtime0 = dtmin
	  dtnext = min(dtmin+ddtwrt,dtmax)

	  call mod_renewal_time_init(nkn,nlvdi)

	  call wrt_flag_inside(rinside,iaout)

	  ius = ifemop('.jas','formatted','new')
	  iuf = ifemop('.frq','formatted','new')
	  iua = ifemop('.jaa','formatted','new')

	  da_out = 0
	  nvar = 1
          call shyfem_init_scalar_file('wrt',nvar,.false.,id)
          da_out(4) = id

	  nrepl = -1				!must still initialize
        end if

        icall = icall + 1
	
c------------------------------------------------------------
c is it time to run the routine?
c------------------------------------------------------------

	call get_act_dtime(dtime)
        if( dtime .lt. dtmin ) return
        if( dtime .gt. dtmax ) return

	call get_act_timeline(aline)

c------------------------------------------------------------
c decide on what to do
c------------------------------------------------------------

c binit		is initial call
c breset	resets concentration
c bcompute	computes residence times and writes to file
c belab		elaborates (accumulates) concentrations

	binit = .false.		!only true for first call
	if( nrepl .lt. 0 ) then
	  nrepl = 0
	  binit = .true.
	end if

	breset = binit
	if( ddtwrt .gt. 0 ) then
	  if( dtime .ge. dtnext ) breset = .true.
	end if
	call custom_reset(dtime,breset)

	belab = .not. binit
	bcompute = .not. binit

	if( bdebug ) then
	  write(6,*) 'WRT debug:'
	  write(6,*) nrepl,binit,breset,belab,bcompute,bnoret
	  write(6,*) dtime,ius
	end if

c------------------------------------------------------------
c elaborate results
c------------------------------------------------------------

	conz = 0.
	if( belab ) then
	  if( blimit ) call wrt_limit_conz(c0,cnv)
	  call wrt_massvolconz(cnv,rinside,vol,mass,volume)
	  call wrt_mass_area(iadim,cnv,massa,vola,conza)
	  call wrt_write_area(iua,aline,iadim,conza,conza0)
	  conz = mass / volume

	  if( bstir ) call wrt_bstir(conz,cnv,rinside)	!stirred tank
          if( bnoret ) call wrt_bnoret(cnv,rinside)	!no return flow

	  tacu = tacu + (dtime-dtime0)
	  call acu_acum(blog,dtime,c0,cnv,vol,rinside,cvacu,volacu)

	  call wrt_restime_summary(ius,dtime,dtime0
     +					,mass,mass0,rcorrect)
	end if

	if( bdebug ) then
	  write(6,*) 'WRT debug (conz): ',conz
	end if

c------------------------------------------------------------
c reset concentrations
c------------------------------------------------------------

        if( breset ) then		!reset concentrations to c0

	  if( bmaster ) then
       	    write(6,*) 'resetting concentrations for renewal time '
     +				,aline
	  end if

c------------------------------------------------------------
c reset variables to compute renewal time (and write to file)
c------------------------------------------------------------

	  if( bcompute ) then	!compute new renewal time
	    rcorrect = 0.	!do not used global correction
	    call acu_comp(da_out
     +				,blog,badj,dtime,c0,ccut,rcorrect
     +				,tacu,cvacu
     +				,cnv,cvres3)
	    call acu_freq(iuf,aline,ctop,rinside,cvres3,volacu)
	    nrepl = nrepl + 1

	    if( bmaster ) then
	      write(6,*) '-------------------------------------------'
	      write(6,*) 'computing res time: ',aline,conz,nrepl
	      write(6,*) '-------------------------------------------'
	    end if
	  end if

	  dtime0 = dtime
	  dtnext = min(dtime+ddtwrt,dtmax)
	  tacu = 0.
	  call acu_reset(cvacu)
	  call acu_reset(volacu)
	  call wrt_breset(c0,cnv,rinside)

	  call wrt_massvolconz(cnv,rinside,vol,mass0,vol0)
	  call wrt_mass_area(iadim,cnv,massa0,vola0,conza0)
	  call wrt_write_area(iua,aline,iadim,conza0,conza0)

	  call wrt_restime_summary(-ius,dtime,dtime0
     +					,mass0,mass0,rcorrect)
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

	subroutine wrt_write_area(iua,aline,narea,massa,massa0)

c computes masses for different areas

	use shympi

	implicit none

	integer iua				!unit
	character*20 aline			!time line
	integer narea				!total number of areas
	double precision massa(0:narea)		!mass of areas
	double precision massa0(0:narea)	!initial mass of areas

	integer i
	real perctot
	real perc(0:narea)

	if( narea .le. 0 ) return

	do i=0,narea
	  perc(i) = 0.
	  if( massa0(i) > 0. ) then
	    perc(i) = 100. * massa(i) / massa0(i)
	  end if
	end do
	perctot = 100. * sum(massa) / sum(massa0)

	if( shympi_is_master() ) then
	  write(iua,*) aline,perc
	end if
	
	end

c******************************************************

	subroutine wrt_mass_area(narea,cnv,massa,vola,conza)

c computes masses for different areas

	use evgeom
	use levels
	use basin
	use shympi

	implicit none

	integer narea				!total number of areas
	real cnv(nlvdi,nkn)			!scalar for which to compute
	double precision massa(0:narea)		!total mass of scalar
	double precision vola(0:narea)		!total volume
	double precision conza(0:narea)		!average conc of scalar

	integer ie,k,ii,ia,l,lmax,nlev,ntot
	real v,conz,area,hdep
	real h(nlvdi)

	real volnode

	if( narea .le. 0 ) return

	massa = 0.
	vola = 0.
	conza = 0.

	ntot = nel_unique

	do ie=1,ntot
	  ia = iarv(ie)
	  if( ia > narea .or. ia < 0 ) goto 99
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

	call shympi_gather_and_sum(massa)
	call shympi_gather_and_sum(vola)

	where( vola > 0. ) conza = massa / vola

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
	use shympi

	implicit none

	real cnv(nlvdi,nkn)		!concentration
	real rinside(nkn)		!flag if node is inside domain
	real vol(nlvdi,nkn)		!volume on node (return)
	double precision mass,volume	!total mass/volume (return)

	integer k,l,lmax,ntot
	real v,conz

	real volnode

        mass = 0.
        volume = 0.

	ntot = nkn_unique

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

	mass = shympi_sum(mass)
	volume = shympi_sum(volume)

	end

c****************************************************

	subroutine wrt_bnoret(cnv,rinside)

c sets concentration to zero outside of domain

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

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

	double precision cvacu(nlvdi,nkn)

        cvacu = 0.

	end

c***************************************************************

	subroutine acu_acum(blog,dtime,c0,cnv,vol,rinside,cvacu,volacu)

c accumulate renewal time

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	logical blog				!use logarithm to compute
	double precision dtime
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

	end

c***************************************************************

	subroutine acu_comp(da_out,blog,badj,dtime,c0,ccut,rcorrect
     +				,tacu,cvacu
     +				,cnv,cvres3)

c compute renewal time and write to file

	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	double precision da_out(4)
	logical blog,badj
	double precision dtime
	real c0
	real ccut
	real rcorrect
	double precision tacu
	double precision cvacu(nlvdi,nkn)		!accumulated conz
	real cnv(nlvdi,nkn)				!last concentration
	real cvres3(nlvdi,nkn)				!computed RT 3D

	integer k,lmax,l,ivar,ierr,id,ishyff
	double precision conz,conze,rconv,corr,cc0
	double precision secs_in_day,ttacu

	double precision dgetpar

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
	id = nint(da_out(4))
	write(6,*) 'writing wrt file for time ',dtime
	call shy_write_scalar_record(id,dtime,ivar,nlvdi,cvres3)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************

	subroutine wrt_restime_summary(ius,dtime,dtime0
     +					,mass,mass0,rcorrect)

c write summuary of WRT computation to file

c perc		percentage of mass still in domain
c restime	renewal time computed by integrating
c restimec	renewal time computed by integrating with correction
c restimel	renewal time computed by fitting regression curve
c resmed	average of renewal times computed
c resstd	standard deviation of renewal time

     	implicit none
	
	integer ius		!negative if reset and summary write
	double precision dtime,dtime0
	double precision mass
	double precision mass0
	real rcorrect

	real dt
	real perc
        real remnant,rlast
	real restime,restimel,restimec
        real resmed,resstd
	character*20 aline

	logical breset
	integer ndata,iu
	integer it
	double precision remint,remlog,remtim
	double precision rsum,rsumsq
	save ndata
	save remint,remlog,remtim
	save rsum,rsumsq
	data ndata / 0 /

	external get_timeline

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
	remtim = remtim + (dtime-dtime0)
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

	call get_timeline(dtime,aline)
!	write(6,1000) aline,perc,restime,restimec
!     +					,restimel,resmed,resstd
        write(iu,1000) aline,perc,restime,restimec
     +					,restimel,resmed,resstd

	it = nint(dtime)
        write(177,2000) it,perc,restime,restimec
     +					,restimel,resmed,resstd

 1000	format(a20,6f10.2)
 2000	format(i10,6f10.2)
	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine acu_freq(iu,aline,ctop,rinside,cvres3,volacu)

c write histogram

	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

	integer iu
	character*20 aline
	real ctop			!cut at this value of renewal time
	real rinside(nkn)		!point is intern
	real cvres3(nlvdi,nkn)
	double precision volacu(nlvdi,nkn)

	integer ndim
	parameter (ndim=100)

	logical bdebug
	integer k,lmax,l,i,ic,ntot
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

	ntot = nkn_unique

	cmax = 0.
        do k=1,ntot
          lmax = ilhkv(k)
          do l=1,lmax
	    cmax = max(cmax,cvres3(l,k))
	  end do
	end do

	cmax = shympi_max(cmax)
	amax = cmax
	if( ctop .gt. 0. .and. amax .gt. ctop ) amax = ctop
	if( ctop .lt. 0. ) amax = -ctop

c---------------------------------------------------------------
c compute average
c---------------------------------------------------------------

	tot = 0.
	vtot = 0.
        do k=1,ntot
	  if( rinside(k) .le. 0. ) cycle
          lmax = ilhkv(k)
          do l=1,lmax
	    conz = cvres3(l,k)
            v = volacu(l,k)
	    tot = tot + conz * v
	    vtot = vtot + v
	  end do
	end do
	tot = shympi_sum(tot)
	vtot = shympi_sum(vtot)

c---------------------------------------------------------------
c compute frequency curve
c---------------------------------------------------------------

	icount = 0
	dcount = 0.

        do k=1,ntot
	  if( rinside(k) .le. 0. ) cycle
          lmax = ilhkv(k)
          do l=1,lmax
	    conz = cvres3(l,k)
            v = volacu(l,k)

	    i = nint(ndim*conz/amax)
	    if (i .lt. 0) i = 0
	    if (i .gt. ndim) i = ndim

	    icount(i) = icount(i) + 1
	    dcount(i) = dcount(i) + v
	  end do
	end do

	call shympi_gather_and_sum(icount)
	call shympi_gather_and_sum(dcount)
	ic = sum(icount)
	dc = sum(dcount)

c---------------------------------------------------------------
c write to file (if master)
c---------------------------------------------------------------

	if( .not. shympi_is_master() ) return

	write(iu,*) 'cmax: ',aline,cmax,amax
	write(iu,2000) 'aver_by_tot: ',aline,tot,vtot,tot/vtot

c---------------------------------------------------------------
c compute and write frequency curve to file
c---------------------------------------------------------------

	dw = 1.
	tot = 0.
	vtot = 0.
	write(iu,*) 'freq_by_bin: ',aline,ndim+1
	do i=0,ndim
	  c = i*amax/ndim
	  val = 100.*icount(i)/(float(ic)*dw)
	  tot = tot + val*c
	  vtot = vtot + val
	  write(iu,*) i,val,icount(i)
	end do
	write(iu,2000) 'aver_by_bin: ',aline,tot,vtot,tot/vtot

	dw = amax/100.
	tot = 0.
	vtot = 0.
	write(iu,*) 'freq_by_res: ',aline,ndim+1
	do i=0,ndim
	  c = i*amax/ndim
	  val = 100.*icount(i)/(float(ic)*dw)
	  tot = tot + val*c
	  vtot = vtot + val
	  write(iu,*) c,val,icount(i)
	end do
	write(iu,2000) 'aver_by_res: ',aline,tot,vtot,tot/vtot

	dw = 1.
	tot = 0.
	vtot = 0.
	write(iu,*) 'freq_by_vol: ',aline,ndim+1
	do i=0,ndim
	  c = i*amax/ndim
	  val = 100.*dcount(i)/(dc*dw)
	  tot = tot + val*c
	  vtot = vtot + val
	  write(iu,*) c,val,icount(i)
	end do
	write(iu,2000) 'aver_by_vol: ',aline,tot,vtot,tot/vtot

c---------------------------------------------------------------
c write out all data to file (for debug and median)
c---------------------------------------------------------------

	!write(76,*) nkn,aline
        !do k=1,nkn
	!  conz = cvres3(1,k)
	!  write(76,*) conz
	!end do

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

 2000	format(a,a20,2e14.6,f14.4)
	end

c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine custom_reset(dtime,breset)

	implicit none

	double precision dtime
	logical breset

	integer i
	integer date,time

	integer, save :: idate = 0

	integer, save :: ndate,ndim
	double precision, save, allocatable :: restime(:)

	character*80 file

c---------------------------------------------------------------
c initialize - convert date to relative time
c---------------------------------------------------------------

	if( idate == -1 ) return

	if( idate == 0 ) then
	  idate = -1
	  call getfnm('wrtrst',file)
	  if( file == ' ' ) return			!no file given
	  call get_reset_time(file,-1,ndim,restime)
	  allocate(restime(ndim))
	  call get_reset_time(file,ndim,ndate,restime)
	  idate = 1
	end if

c---------------------------------------------------------------
c see if we have to reset
c---------------------------------------------------------------

	if( idate > ndate ) return
	if( dtime < restime(idate) ) return

c---------------------------------------------------------------
c ok, reset needed - advance to next reset time
c---------------------------------------------------------------

	idate = idate + 1
	breset = .true.

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************

	subroutine get_reset_time(file,ndim,n,restime)

c gets custom reset time from file

	implicit none

	character*(*) file
	integer ndim		!ndim==0 => check how many dates are given
	integer n
	double precision restime(n)

	integer ianz,ios,nline,i
	integer date,time
	double precision dtime,atime,atime0,dtold
	double precision d(2)
	character*80 line
	character*20 aline
	logical bdebug

	integer iscand

	n = 0
	nline = 0
	bdebug = .true.
	if( ndim == -1 ) bdebug = .false.

	open(1,file=file,status='old',form='formatted',iostat=ios)

	if( bdebug ) then
	  if( ios /= 0 ) then
	    write(6,*) 'cannot open custom reset file: ',trim(file)
	    stop 'error stop get_reset_time: opening file'
	  else
	    write(6,*) 'reading custom reset file: ',trim(file)
	  end if
	end if
	if( ios /= 0 ) return

	call get_absolute_ref_time(atime0)

	do
	  read(1,'(a)',iostat=ios) line
	  nline = nline + 1
	  if( ios /= 0 ) exit
	  ianz = iscand(line,d,2)
	  if( ianz == 0 ) then
	    cycle
	  else if( ianz == 1 ) then
	    date = nint(d(1))
	    time = 0
	  else if( ianz == 2 ) then
	    date = nint(d(1))
	    time = nint(d(2))
	  else
	    write(6,*) 'parse error: ',ianz
	    write(6,*) 'line: ',trim(line)
	    write(6,*) 'file: ',trim(file)
	    stop 'error stop get_reset_time: parse error'
	  end if

	  n = n + 1
	  if( ndim == -1 ) cycle
	  if( n > ndim ) then
	    write(6,*) 'n,ndim: ',n,ndim
	    stop 'error stop get_reset_time: dimension error ndim'
	  end if

	  call dts_to_abs_time(date,time,atime)
	  dtime = atime - atime0

	  restime(n) = dtime		!insert relative time
	end do

	if( ios > 0 ) then
	  write(6,*) 'read error...'
	  write(6,*) 'file: ',trim(file)
	  write(6,*) 'line number: ',nline
	  stop 'error stop get_reset_time: read error'
	end if

	if( bdebug ) then
	  write(6,*) 'custom reset times: ',n
	  dtold = restime(1) - 1
	  do i=1,n
	    dtime = restime(i)
	    atime = atime0 + dtime
	    call dts_format_abs_time(atime,aline)
	    write(6,*) i,aline
	    if( dtime <= dtold ) then
	      write(6,*) 'times in custom reset must be ascending...'
	      stop 'error stop get_reset_time: wrong order'
	    end if
	  end do
	end if

	close(1)

	end

c**********************************************************************

