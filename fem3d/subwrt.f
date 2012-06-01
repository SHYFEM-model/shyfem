c
c $Id: subwrt.f,v 1.58 2010-03-08 17:46:45 georg Exp $
c
c routines dealing with water residence time
c
c contents :
c
c revision log :
c
c 24.10.2011	ggu	new file copied from subcus.f (jamal)
c 28.02.2012	ggu&deb	completely restructured
c 16.03.2012	ggu	use idtreset=-1 for no residence computation
c
c******************************************************************

        subroutine residence_time

c computes residence time online - one value for whole lagoon

        implicit none

        include 'param.h'

        character*80 descrp
        common /descrp/ descrp

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
	integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer nen3v(3,1)
        common /nen3v/nen3v
        integer ilhkv(1)
        common /ilhkv/ilhkv
        integer iarv(1)
        common /iarv/iarv
        real cnv(nlvdim,nkndim)
        common /cnv/cnv
        real v1v(1)
        common /v1v/v1v
	real hev(neldim)
	common /hev/hev
	real hlv(nlvdim)
	common /hlv/hlv

	logical breset,bcompute,binit,belab
        logical bnoret,bstir
	logical blog,badj

	integer iaout,itmin,itmax,idtreset
	real c0
	real ctop,ccut
	real percmin

        real conz,perc,rcorrect
        double precision mass,volume
	
	character*80 title
	integer nvers	

	real cvres3(nlvdim,nkndim)      	!computed RT 3D
	real vol(nlvdim,nkndim)      		!volume
	double precision cvacu(nlvdim,nkndim)   !conz integrated
	double precision volacu(nlvdim,nkndim)  !volume integrated
	save cvacu,volacu

	integer ifileo
        real volnode
	integer ifemop

	integer iu,nb3,iuf
	save iu,nb3,iuf
	integer nrepl
	save nrepl
	integer it0
	save it0
	double precision tacu
	save tacu
        double precision mass0
        save mass0

        integer icall
        save icall
        data icall / 0 /

	if( icall .lt. 0 ) return

c------------------------------------------------------------
c parameters
c------------------------------------------------------------
c
c bnoret	true if no return flow is used (concentrations outside
c		are explicitly set to 0)
c bstir		simulates completely stirred tank
c		(replaces at every time step conz with average conz)
c blog          use logarithmic regression to compute residence time
c badj          adjust residence time for tail of distribution
c
c percmin	percentage to reach -> after this stop computation
c		use 0 if no premature end is desired
c iaout		area code of elements out of lagoon (used for init and retflow)
c		use -1 to if no outside areas exist
c c0		initial concentration of tracer
c
c itmin		time from when to compute residence time (-1 for start of sim)
c itmax		time up to when to compute residence time (-1 for end of sim)
c idtreset	time step to reset concentration to c0
c		use 0 if no reset is desired
c		use -1 if no residence time computation is desired
c
c ctop          maximum to be used for frequency curve
c ccut          cut residence time at this level (for res time computation)

c--------------------------
c default settings
c--------------------------

        bnoret = .false.	!no return flow
	bstir = .false.		!stirred tank
        blog = .false.		!compute residence time with log fitting
        badj = .true.		!adjust residence time for tail

	percmin = 0.		!minimum percentage for remnant function
	iaout = -1		!areas considered outside
	c0 = 1.			!initial concentration

	itmin = -1		!compute from start of simulation
	itmax = -1		!compute to end of simulation
	idtreset = 0		!no reset of concentrations (old default)
	idtreset = -1		!no residence time computation

	ctop = 0.		!max for frequency curve
	ccut = 0.		!max fro residence time

c--------------------------
c customization
c--------------------------

	badj = .false.
	itmin = 0
	idtreset = 3*nint( 30.5 * 86400 )		!one month is 30.5 days
	idtreset = -1

        !badj = .false.
        !blog = .true.

	!idtreset = nint( 30.5 * 86400 )		!one month is 30.5 days
	!idtreset = nint( 86400. )		!one month is 30.5 days
	!itmin = -1

	!ccut = 500.
	!ctop = 400.

c------------------------------------------------------------
c initialization
c------------------------------------------------------------

        if( icall .eq. 0 ) then
          write(6,*) 'initialization of WRT routine residence time'

	  if( idtreset .lt. 0 ) then
	    icall = -1
            write(6,*) 'no residence time computation'
	    return
	  end if

	  iu = ifemop('.jas','formatted','new')
	  iuf = ifemop('.frq','formatted','new')
	  nb3 = ifemop('.wrt','unform','new')
	  if( nb3 .le. 0 ) stop 'error stop open_nos_file: opening file'
	  nvers = 3
	  title = descrp
	  call whnos(nb3,nvers,nkn,nel,nlv,1,ilhkv,hlv,hev,title)

	  nrepl = -1				!must still initialize
        end if

        icall = icall + 1
	if( itmin .eq. -1 ) itmin = itanf
	if( itmax .eq. -1 ) itmax = itend
	
c------------------------------------------------------------
c is it time to run the routine?
c------------------------------------------------------------

        if( it .lt. itmin ) return
        if( it .gt. itmax ) return

c------------------------------------------------------------
c decide on what to do
c------------------------------------------------------------

	binit = .false.
	if( nrepl .lt. 0 ) then
	  nrepl = 0
	  binit = .true.
	end if

	breset = binit
	if( idtreset .gt. 0 ) then
	  if( it-it0 .ge. idtreset ) breset = .true.
	end if
	if( it .eq. itmax ) breset = .true.	!last time step

	belab = .not. binit
	bcompute = .not. binit

	!write(6,*) nrepl,binit,breset,belab,bcompute
	!write(6,*) it,iu

c------------------------------------------------------------
c flag nodes that are inside lagoon (v1v(k)=1)
c------------------------------------------------------------

	call wrt_flag_inside(v1v,iaout)

c------------------------------------------------------------
c elaborate results
c------------------------------------------------------------

	if( belab ) then
	  call wrt_massvolconz(cnv,v1v,vol,mass,volume)
	  conz = mass / volume

	  if( bstir ) call wrt_bstir(conz,cnv,v1v)	!stirred tank
          if( bnoret ) call wrt_bnoret(cnv,v1v)		!no return flow

	  tacu = tacu + (it-it0)
	  call acu_acum(blog,it,c0,cnv,vol,cvacu,volacu)

	  call wrt_restime_summary(iu,it,it0,mass,mass0,rcorrect)
	end if

c------------------------------------------------------------
c reset concentrations
c------------------------------------------------------------

        if( breset ) then		!reset concentrations to c0

       	  write(6,*) 'resetting concentrations for residence time ',it

c------------------------------------------------------------
c reset variables to compute residence time
c------------------------------------------------------------

	  if( bcompute ) then	!compute new residence time
	    rcorrect = 0.	!do not used global correction
	    call acu_comp(nb3,blog,badj,it,c0,ccut,rcorrect
     +				,tacu,cvacu
     +				,cnv,cvres3)
	    call acu_freq(iuf,it,ctop,cvres3,volacu)
	    nrepl = nrepl + 1

	    write(6,*) '-------------------------------------------'
	    write(6,*) 'computing res time: ',it,conz,nrepl
	    write(6,*) '-------------------------------------------'
	  end if

	  it0 = it
	  tacu = 0.
	  call acu_reset(cvacu)
	  call acu_reset(volacu)
	  call wrt_breset(c0,cnv,v1v)

	  call wrt_massvolconz(cnv,v1v,vol,mass,volume)
	  mass0 = mass

	  call wrt_restime_summary(-iuf,it,it0,mass,mass0,rcorrect)		!reset
	end if

c------------------------------------------------------------
c finish computation if mass is below threshold
c------------------------------------------------------------

        if( mass0 .ne. 0. ) then
	  perc = mass / mass0
          if( perc .lt. percmin ) then
                stop 'finished computing residence time'
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

	implicit none

	include 'param.h'

	real rinside(1)
	integer iaout		!area code for outside elements

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1)
        common /nen3v/nen3v
        integer iarv(1)
        common /iarv/iarv

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

	implicit none

	include 'param.h'

	real c0				!concentration for nodes inside domain
	real cnv(nlvdim,1)
	real rinside(1)			!flag if node is inside domain

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer ilhkv(1)
        common /ilhkv/ilhkv

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

	implicit none

	include 'param.h'

	real c0				!concentration to impose
	real cnv(nlvdim,1)
	real rinside(1)			!flag if node is inside domain

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer ilhkv(1)
	common /ilhkv/ilhkv

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

	subroutine wrt_massvolconz(cnv,rinside,vol,mass,volume)

c computes mass and volume on internal nodes

	implicit none

	include 'param.h'

	real cnv(nlvdim,1)
	real rinside(1)			!flag if node is inside domain
	real vol(nlvdim,1)
	double precision mass,volume

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer ilhkv(1)
	common /ilhkv/ilhkv

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

	implicit none

	include 'param.h'

        real cnv(nlvdim,1)
	real rinside(1)			!flag if node is inside domain

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer ilhkv(1)
	common /ilhkv/ilhkv

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

	implicit none

	include 'param.h'

	double precision cvacu(nlvdim,nkndim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer ilhkv(1)
	common /ilhkv/ilhkv

	integer k,lmax,l

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            cvacu(l,k) = 0.
          end do
        end do

	end

c***************************************************************

	subroutine acu_acum(blog,it,c0,cnv,vol,cvacu,volacu)

	implicit none

	include 'param.h'

	logical blog				!use logarithm to compute
	integer it
	real c0					!value used for initialization
	real cnv(nlvdim,nkndim)
	real vol(nlvdim,nkndim)
	double precision cvacu(nlvdim,nkndim)
	double precision volacu(nlvdim,nkndim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer ilhkv(1)
	common /ilhkv/ilhkv

	integer k,l,lmax
	real rl,conz,dt

	call get_timestep(dt)

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            conz = cnv(l,k)
	    if ( conz .gt. c0 ) conz = c0
	    if ( conz .lt. 0. ) conz = 0.
	    if( blog ) then
	      if( conz .gt. 0 ) then
	        rl = - log(conz/c0)
                cvacu(l,k) = cvacu(l,k) + rl
	      end if
	    else
              cvacu(l,k) = cvacu(l,k) + conz*dt
	    end if
            volacu(l,k) = volacu(l,k) + vol(l,k)*dt
          end do
        end do

	!l = 1
	!k = 100
	!write(6,*) cnv(l,k),cvacu(l,k),volacu(l,k)

	end

c***************************************************************

	subroutine acu_comp(nb3,blog,badj,it,c0,ccut,rcorrect
     +				,tacu,cvacu
     +				,cnv,cvres3)

c compute residence time and write to file

	implicit none

	include 'param.h'

	integer nb3
	logical blog,badj
	integer it
	real c0
	real ccut
	real rcorrect
	double precision tacu
	double precision cvacu(nlvdim,nkndim)		!accumulated conz
	real cnv(nlvdim,nkndim)				!last concentration
	real cvres3(nlvdim,nkndim)			!computed RT 3D

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer ilhkv(1)
	common /ilhkv/ilhkv

	integer k,lmax,l,ivar,ierr
	real conz,conze,rconv,corr
	real secs_in_day

c---------------------------------------------------------------
c set parameters
c---------------------------------------------------------------

	secs_in_day = 86400.

	tacu = tacu / secs_in_day
	rconv = 1. / secs_in_day

c---------------------------------------------------------------
c compute residence times -> put in cvres3
c---------------------------------------------------------------

        do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            conz = cvacu(l,k)
	    if( blog ) then
	      conz = tacu / conz
	    else
              conz = rconv * conz / c0		!convert to days
	      if( badj ) then
		if( rcorrect .le. 0. ) then
	          conze = cnv(l,k) / c0
	          if( conze .ge. 1 ) conze = 0.
	          if( conze .le. 0 ) conze = 0.
		  corr = 1. /  ( 1. - conze )
		else
		  corr = rcorrect
		end if
                conz = corr * conz 		!adjusted res time
	      end if
	    end if
	    if ( ccut .gt. 0. .and. conz .gt. ccut ) conz = ccut
            cvres3(l,k) = conz
          end do
        end do

c---------------------------------------------------------------
c write to file
c---------------------------------------------------------------

	ivar = 99
        call wrnos(nb3,it,ivar,nlvdim,ilhkv,cvres3,ierr)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	end

c**********************************************************************

	subroutine wrt_restime_summary(iu,it,it0,mass,mass0,rcorrect)

c perc		percentage of mass still in domain
c restime	residence time computed by integrating
c restimec	residence time computed by integrating with correction
c restime1	residence time computed by fitting regression curve
c resmed	average of residence times computed
c resstd	standard deviation of residence time

     	implicit none
	
	integer iu		!negative if reset and summary write
	integer it,it0
	double precision mass
	double precision mass0
	real rcorrect

	real dt
	real perc
        real remnant,rlast
	real restime,restime1,restimec
        real resmed,resstd

	integer ndata
	double precision remint,remlog,remtim
	double precision rsum,rsumsq
	save ndata
	save remint,remlog,remtim
	save rsum,rsumsq
	data ndata / 0 /

	if( iu .le. 0 ) then	!reset
	  ndata = 0
	  remint = 0.
	  remlog = 0.
	  remtim = 0.
	  rsum = 0.
	  rsumsq = 0.
	  !if( iu .ne. 0 ) then
	  ! write(-iu,1000) it,perc,restime,restimec,restime1,resmed,resstd
	  !end if
	  return
	end if

	call get_timestep(dt)

	remnant = 0.
        if( mass0 .gt. 0. ) remnant = mass/mass0
        perc = 100.*remnant

	remint = remint + remnant*dt	!integrated remnant function
	restime = remint/86400.		!residence time in days

	rlast = remnant
	if( rlast .ge. 1. ) rlast = 0.
	rcorrect = 1. / (1.-rlast)
	restimec = rcorrect * restime	!corrected residence time

	remlog = remlog - log(remnant)
	remtim = remtim + (it-it0)
	restime1 = 0.
	if( remlog .gt. 0. ) restime1 = ( remtim / remlog ) / 86400.

	ndata = ndata + 1
	rsum = rsum + restime1
	rsumsq = rsumsq + restime1*restime1
	resmed = rsum / ndata
	resstd = sqrt( rsumsq/ndata - resmed*resmed )

        !write(6,1000) it,perc,restime,restimec,restime1,resmed,resstd
        write(iu,1000) it,perc,restime,restimec,restime1,resmed,resstd

 1000	format(i10,6f10.2)
	end

c**********************************************************************
c**********************************************************************
c**********************************************************************

	subroutine acu_freq(iu,it,ctop,cvres3,volacu)

c write histogram

	implicit none

	include 'param.h'

	integer iu
	integer it
	real ctop			!cut at this value of residence time
	real cvres3(nlvdim,nkndim)
	double precision volacu(nlvdim,nkndim)

        integer ilhkv(1)
        common /ilhkv/ilhkv
	real hev(neldim)
	common /hev/hev
	real hlv(nlvdim)
	common /hlv/hlv

	logical breset,bcompute,binit,belab
        logical bnoret,bstir
	logical blog,badj

	integer ndim
	parameter (ndim=100)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

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
	write(iu,*) 'cmax: ',it,cmax,amax

c---------------------------------------------------------------
c compute average
c---------------------------------------------------------------

	tot = 0.
	vtot = 0.
        do k=1,nkn
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

