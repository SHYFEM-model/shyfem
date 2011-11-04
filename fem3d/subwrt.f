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
c
c******************************************************************

        subroutine residence_time

c computes residence time online - one value for whole lagoon

        implicit none

        include 'param.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        real grav,fcor,dcor,dirn,rowass,roluft

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

        integer ie,ii,k,lmax,l,ia
        logical bnoret,breset,bstir
        real vol,conz,perc,percmin
        double precision mass,volume
        double precision massaux,volaux
        real volnode
	integer ifemop
	integer iaout,itmin,idtreset
	real dt,c0
	real restime,restime1,restimec
	real remnant,rlast
	real resmed,resstd

	integer iu,it0,ndata
	save iu,it0,ndata
	double precision remint,remlog,remtim
	save remint,remlog,remtim
	double precision rsum,rsumsq
	save rsum,rsumsq
        double precision mass0
        save mass0

        integer icall
        save icall
        data icall / 0 /

c------------------------------------------------------------
c parameters
c------------------------------------------------------------
c
c bnoret	true if no return flow is used (conzentrations outside
c		are explicitly set to 0)
c bstir		simulates completely stirred tank
c		(replaces at every time step conz with average conz)
c percmin	percentage to reach -> after this stop computation
c		use 0 if no premature end is desired
c iaout		area code of elements out of lagoon (used for init and retflow)
c		use -1 to if no outside areas exist
c c0		initial concentration of tracer
c itmin		time from when to compute residence time
c idtreset	time step to reset concentration to c0
c		use 0 if no reset is desired
c
c--------------------------
c default settings
c--------------------------
        bnoret = .false.
	bstir = .true.
	bstir = .false.
	percmin = 0.
	iaout = -1
	c0 = 1.
	itmin = 0
	idtreset = 0
	idtreset = nint( 3 * 30.5 * 86400 )		!one month is 30.5 days
c--------------------------
c nador
c--------------------------
c	bnoret = .true.
c	percmin = 1.
c	iaout = 0
c--------------------------
c alimini
c--------------------------
c	idtreset = nint(30.5*86400)
c--------------------------
c mar menor
c--------------------------
c 	just use default settings
c--------------------------
c
c------------------------------------------------------------
c do not change anything after this point
c------------------------------------------------------------

c------------------------------------------------------------
c is it time to run the routine?
c------------------------------------------------------------

        if( it .le. itmin ) return

c------------------------------------------------------------
c initialization -> decide on reset
c------------------------------------------------------------

	breset = .false.		!normally do not reset

        if( icall .eq. 0 ) then
          write(6,*) 'initialization of routine jamal'
	  iu = ifemop('.jam','formatted','new')
	  breset = .true.		!always reset at first call
	  it0 = it
        end if

	if( idtreset .gt. 0 ) then
	  if( it-it0 .ge. idtreset ) breset = .true.
	end if

c------------------------------------------------------------
c flag nodes that are inside lagoon (v1v(k)=1)
c------------------------------------------------------------

	call wrt_flag_inside(v1v,iaout)

c------------------------------------------------------------
c reset concentrations
c------------------------------------------------------------

        if( breset ) then		!reset concentrations to c0
          write(6,*) 'resetting concentrations in jamal at time ',it
          do k=1,nkn
            lmax = ilhkv(k)
	    conz = 0.
            if( v1v(k) .ne. 0. ) conz = c0
            do l=1,lmax
              cnv(l,k) = conz
            end do
          end do
	end if

c------------------------------------------------------------
c total mass (only for nodes inside lagoon)
c------------------------------------------------------------

        mass = 0.
        volume = 0.
        do k=1,nkn
          if( v1v(k) .ne. 0. ) then
            lmax = ilhkv(k)
            do l=1,lmax
              vol = volnode(l,k,+1)
              conz = cnv(l,k)
              mass = mass + vol*conz
              volume = volume + vol
            end do
          end if
        end do

c------------------------------------------------------------
c stirred tank?
c------------------------------------------------------------

	if( bstir ) then
	  conz = mass / volume
          do k=1,nkn
            if( v1v(k) .ne. 0. ) then
              lmax = ilhkv(k)
              do l=1,lmax
                cnv(l,k) = conz
	      end do
	    end if
	  end do
	end if

        massaux = 0.
        volaux = 0.
        do k=1,nkn
          if( v1v(k) .ne. 0. ) then
            lmax = ilhkv(k)
            do l=1,lmax
              vol = volnode(l,k,+1)
              conz = cnv(l,k)
              massaux = massaux + vol*conz
              volaux = volaux + vol
            end do
          end if
        end do

	!write(67,*) it,mass,massaux,volume,volaux
	!write(67,*) it,mass,massaux

c------------------------------------------------------------
c reset variables to compute residence time
c------------------------------------------------------------

        if( breset ) then
	  mass0 = mass 
	  remint = 0.
	  remlog = 0.
	  remtim = 0.
	  it0 = it
	  ndata = 0
	  rsum = 0.
	  rsumsq = 0.
	end if

c------------------------------------------------------------
c write to file
c------------------------------------------------------------

	call get_timestep(dt)

	remnant = 0.
        if( mass0 .gt. 0. ) remnant = mass/mass0
        perc = 100.*remnant

	remint = remint + remnant*dt	!integrated remnant function
	restime = remint/86400.		!residence time in days

	rlast = remnant
	if( rlast .ge. 1. ) rlast = 0.
	restimec = restime/(1.-rlast)	!corrected residence time

	remlog = remlog - log(remnant)
	remtim = remtim + (it-it0)
	restime1 = 0.
	if( remlog .gt. 0. ) restime1 = ( remtim / remlog ) / 86400.

	ndata = ndata + 1
	rsum = rsum + restime1
	rsumsq = rsumsq + restime1*restime1
	resmed = rsum / ndata
	resstd = sqrt( rsumsq/ndata - resmed*resmed )

c perc		percentage of mass still in domain
c restime	residence time computed by integrating
c restimec	residence time computed by integrating with correction
c restime1	residence time computed by fitting regression curve
c resmed	average of residence times computed
c resstd	standard deviation of residence time

        write(iu,1000) it,perc,restime,restimec,restime1,resmed,resstd

c------------------------------------------------------------
c finish computation if mass is below threshold
c------------------------------------------------------------

        if( mass0 .ne. 0. .and. perc .lt. percmin ) then
                stop 'finished computing'
        end if

c------------------------------------------------------------
c no return flow -> set outside areas to 0
c------------------------------------------------------------

        if( bnoret ) then
          do k=1,nkn
            if( v1v(k) .eq. 0. ) then
              lmax = ilhkv(k)
              do l=1,lmax
                cnv(l,k) = 0.
              end do
            end if
          end do
        end if

c------------------------------------------------------------
c remember initialization
c------------------------------------------------------------

        icall = icall + 1

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	return
 1000	format(i10,6f10.2)
        end

c*****************************************************************

	subroutine wrt_flag_inside(v1v,iaout)

c flags inside nodes (nodes with area code different from iaout)
c
c on return v1v(k) = 1 for nodes inside domain

	implicit none

	real v1v(1)
	integer iaout

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1)
        common /nen3v/nen3v
        integer iarv(1)
        common /iarv/iarv

	integer k,ie,ii,ia

        do k=1,nkn
          v1v(k) = 0.
        end do

        do ie=1,nel
          ia = iarv(ie)
          if( ia .ne. iaout ) then
              do ii=1,3
                k = nen3v(ii,ie)
                v1v(k) = 1.
              end do
          end if
        end do

	end

c*****************************************************************

