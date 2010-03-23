c
c $Id: subsed.f,v 1.3 2003/03/25 14:08:55 georg Exp $
c
c routines for concentration - sediment
c
c contents :
c
c subroutine sed2sh				shell for conz - sediments
c
c revision log :
c
c 16.06.1998	ggu	new names for some parameters
c			rkpar -> chpar , cref -> conref
c 19.06.1998	ggu	bug fix (ifileo nor declared)
c 20.06.1998	ggu	general cleanup
c 24.06.1998	ggu	heat flux enabled
c 26.01.1999	ggu	new routines (3D) to write files NOS
c 07.06.1999	ggu	bug fix -> bconz... not saved
c 28.10.1999	ggu	names changed
c
c*********************************************************************
c
	subroutine sed2sh
c
c shell for conz
c
c use negative values for iconz to run special tests
c
c  1	normal run
c  0	no transport/diffusion
c -1	diffusion test
c -2	advection test
c -3	flux test
c -4	dilution test
c
c written 09.01.94 by ggu  (from scratch)
c revised 20.01.94 by ggu  $$iclin - iclin not used to compute volume
c revised 20.01.94 by ggu  $$zov0 - zeov in conz adjusted
c revised 31.01.94 by ggu  $$conzbc - implementation of bc for conz
c revised 03.02.94 by ggu  $$zov00 - zeov in sp159 adjusted
c revised 04.02.94 by ggu  $$azpar - azpar used to compute transport
c revised 07.02.94 by ggu  $$istot - istot for fractional time step
c revised 03.12.97 by ggu  $$TS - temp/salt included
c
	implicit none

	include 'param.h'

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi

        integer nen3v(3,1)			!element index
        common /nen3v/nen3v
	real zenv(3,1),zeov(3,1)		!water level (new and old)
	common /zenv/zenv, /zeov/zeov
	real hm3v(3,1)				!depth
	common /hm3v/hm3v
	real v1v(1)				!aux array
	common /v1v/v1v
	real v2v(1)				!aux array
	common /v2v/v2v

	real sedev(3,neldim)		!concentration per element
	real sednv(nkndim)		!concentration per node
	real rsedv(nkndim)		!concentration (boundary condition)
	common /sedev/sedev, /sednv/sednv, /rsedv/rsedv
	save /sedev/, /sednv/, /rsedv/

	real sedbot(nkndim)		!thickness of sediments at bottom
	common /sedbot/sedbot
	save /sedbot/
c local
	integer icall
	integer ie,ii
	integer k
	real ckmax,ckmin,cemax,cemin
c	double precision conref,chpar,dt
	real conref,temref,salref
	real chpar,dt,azpar
	logical bdebug
	logical bconz,btemp,bsalt
	integer iconz,itemp,isalt
	integer istot,isact
	integer iu
	integer iuc,ius,iut
	integer itmcon,idtcon
	real wsink,rhosed,sedbc
c function
	real getpar
	integer iround
c save & data
	save icall
	save bdebug
	save chpar,dt,azpar
	save conref,temref,salref
	save iconz,itemp,isalt
	save bconz,btemp,bsalt
	save istot
	save iuc,ius,iut
	save itmcon,idtcon

	data icall /0/

	if(icall.eq.-1) return

c----------------------------------------------------------------
c set initial values
c----------------------------------------------------------------

	if(icall.eq.0) then

	  iconz=iround(getpar('isedi'))

	  if( iconz .le. 0 ) icall = -1
	  if( icall .eq. -1 ) return

	  conref=getpar('sedref')

	  chpar=getpar('chpar')
	  call getaz(azpar)
	  istot=iround(getpar('istot'))		!$$istot
	  dt=idt

	  bdebug = iround(getpar('levdbg')) .ge. 3

	  call conzin(iconz,conref,sednv,sedev)	!initialize arrays
	  sedbc = 10.
	  call setbc(sedbc,rsedv,flag)		!set boundary condition

	  do k=1,nkn
	    sedbot(k) = 0.
	  end do

          iu = 61
          itmcon = iround(getpar('itmcon'))
          idtcon = iround(getpar('idtcon'))

	  iuc = iu
          call confop(iuc,itmcon,idtcon,2,'sed')

	end if

	icall=icall+1

c----------------------------------------------------------------
c main loop
c----------------------------------------------------------------

	do isact=1,istot		!$$istot
	  call conz2d(sednv,sedev,v1v,dt,chpar,azpar,istot,isact)
          call conzbc(sednv,sedev,v1v,rsedv,flag,azpar) !boundary conditions
	end do

	wsink = 0.1e-03   ! sinking velocity m/s
	rhosed = 1.8e03       ! sediment density, kg/m3
	call sinking(sedbot,sedev,v1v,v2v,dt,wsink,rhosed)

c	ie=1
c	write(69,*) ie,it
c	write(69,*) (hm3v(ii,ie),ii=1,3)
c	write(69,*) (zenv(ii,ie),ii=1,3)
c	write(69,*) (hm3v(ii,ie)+zenv(ii,ie),ii=1,3)

c----------------------------------------------------------------
c debug and write to file
c----------------------------------------------------------------

	if( bdebug ) then
	    call cmima(nkn,nel,sednv,sedev,ckmin,ckmax,cemin,cemax)
	    write(6,*) 'cemax : ',it,cemin,cemax
	end if

        call confil(iuc,itmcon,idtcon,30,1,sednv)
        call confil(iuc,itmcon,idtcon,33,1,sedbot)

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c**************************************************************

	subroutine kim(sedbot,sedev,v1v,v2v,dt,wsink,rhosed)

c implements sinking

	implicit none

	real sedbot(1)		!thickness of bottom sediment
	real sedev(3,1)		!concentration of sediment in water column
	real v1v(1)		!auxiliary vector
	real v2v(1)		!auxiliary vector
	real dt			!time step
	real wsink		!sinking velocity
	real rhosed		!density of sediment

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
	real hm3v(3,1)
	common /hm3v/hm3v
	real zenv(3,1)
	common /zenv/zenv
	real ev(13,1)
	common /ev/ev
        real unv(1),vnv(1)
        common /unv/unv, /vnv/vnv

	integer k,ie,ii
	real area,hdepth,zlev,htot,vol
	real conznew,dconz,conzold
	real frac,rlength,dh
	real tau,tau0
	real day
	real ut,vt,u,v, tken, tkmin, tkmax
	real totmas

	rlength = wsink * dt
	day = 24 * 3600
	tau0 = 30 * day

	totmas = 0.

c------------------------------------------------------------
c initialize
c------------------------------------------------------------

	do k=1,nkn
	  v1v(k) = 0.
	  v2v(k) = 0.
	end do

c------------------------------------------------------------
c loop over elements
c------------------------------------------------------------

	do ie=1,nel

	  area = 4. * ev(10,ie)		!area of sub element
	  ut = unv(ie)
	  vt = vnv(ie)

	  do ii=1,3
	    k = nen3v(ii,ie)
	    dh = sedbot(k)		!thickness of sediment at bottom
	    conzold = sedev(ii,ie)	!old concentration
	    hdepth = hm3v(ii,ie)	!depth of element
	    zlev = zenv(ii,ie)		!water level
	    htot = hdepth + zlev	!total depth
	    vol = htot * area

	    u = ut / htot		!barotropic velocities
	    v = vt / htot
	
c	    call ???(conzold,dh,htot,conznew...)

	    sedev(ii,ie) = conznew	!store new concentration
	    dconz = conzold - conznew
	    v1v(k) = v1v(k) + vol * dconz	!mass in finite volume
	    v2v(k) = v2v(k) + rhosed * area
	  end do

	end do

c------------------------------------------------------------
c compute bottom thickness of sediments
c------------------------------------------------------------

	do k=1,nkn
	  dh = v1v(k) / v2v(k)
	  sedbot(k) = sedbot(k) + dh
	end do

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	end

c**************************************************************

	subroutine sinking(sedbot,sedev,v1v,v2v,dt,wsink,rhosed)

c implements sinking (carniel 2000)

	implicit none

	real sedbot(1)		!thickness of bottom sediment
	real sedev(3,1)		!concentration of sediment in water column
	real v1v(1)		!auxiliary vector
	real v2v(1)		!auxiliary vector
	real dt			!time step
	real wsink		!sinking velocity
	real rhosed		!density of sediment

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
	real hm3v(3,1)
	common /hm3v/hm3v
	real zenv(3,1)
	common /zenv/zenv
	real ev(13,1)
	common /ev/ev
        real unv(1),vnv(1)
        common /unv/unv, /vnv/vnv

	integer k,ie,ii
	real area,hdepth,zlev,htot,vol
	real conznew,dconz,conzold
	real frac,rlength,dh
	real tau,tau0
	real day
	real ut,vt,u,v, tken, tkmin, tkmax
	real totmas

	rlength = wsink * dt
	day = 24 * 3600
	tau0 = 30 * day

	totmas = 0.

c------------------------------------------------------------
c initialize
c------------------------------------------------------------

	do k=1,nkn
	  v1v(k) = 0.
	  v2v(k) = 0.
	end do

c------------------------------------------------------------
c loop over elements
c------------------------------------------------------------

	do ie=1,nel

	  area = 4. * ev(10,ie)		!area of sub element
	  ut = unv(ie)
	  vt = vnv(ie)

	  do ii=1,3
	    k = nen3v(ii,ie)
	    conzold = sedev(ii,ie)
	    hdepth = hm3v(ii,ie)
	    zlev = zenv(ii,ie)
	    htot = hdepth + zlev
	    vol = htot * area

	    u = ut / htot
	    v = vt / htot
	    tken=sqrt(u**2+v**2)
	    if(tken.lt.tkmin) tkmin=tken
	    if(tken.gt.tkmax) tkmax=tken
		wsink=0.1e-03
	if(tken.gt.0.5.and.htot.gt.5.) wsink=0.01e-10
	if(tken.gt.0.25.and.htot.gt.1.and.htot.le.5.) wsink=0.01e-10
	if(tken.gt.0.15.and.htot.le.1.) wsink=0.01e-10
cc	    if(tken.gt.0.7) then
c			wsink = 0.01e-10
c			    else
c			wsink = 0.1e-03  
c			    endif
	
c 1st formulation
	    tau = tau0
c 2nd formulation
	    tau = htot / wsink

c	    conznew = (1.-wsink*dt/htot) * conzold
	    conznew = conzold * exp ( - dt / tau )

	    if( conznew .lt. 0. ) conznew = 0.

	    totmas = totmas + vol * conznew

	    sedev(ii,ie) = conznew
	    dconz = conzold - conznew
	    v1v(k) = v1v(k) + vol * dconz	!mass in finite volume
	    v2v(k) = v2v(k) + rhosed * area
	  end do

	end do

c------------------------------------------------------------
c compute bottom thickness of sediments
c------------------------------------------------------------

	do k=1,nkn
	  dh = v1v(k) / v2v(k)
	  sedbot(k) = sedbot(k) + dh
	end do

c------------------------------------------------------------
c end of routine
c------------------------------------------------------------

	    write(76,*) totmas
	end

c**************************************************************

