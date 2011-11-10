c
c $Id: weutro_sedim.f,v 1.3 2008-04-17 14:53:56 georg Exp $
c
c weutro_sedim - sediment routines for weutro
c
c revision log :
c
c 20.06.2003    ggu&dmk new routine for sediments
c 18.04.2008    utility routines for bio3d taken out
c
c********************************************************************

	subroutine wsedim(k,t,dt,vol,depth,vel,stp,c,cs)

c EUTRO 0-Dimensional (Sediments)

        implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )
        integer nsstate          !total number of state parameters
        parameter( nsstate = 2 )

	integer k
	real t			!actual time [day]
        real dt                 !time step [day]
        real vol                !volume [m**3]
        real depth              !depth of box [m]
        real vel                !velocity [m/s]
        real stp                !temperature [C]
        real c(nstate)          !state variable [mg/L] == [g/m**3]
        real cold(nstate)       !old state variable [mg/L] == [g/m**3]
        real cs(nsstate)        !state variable (sediments) [mg/L] == [g/m**3]

	include 'donata.h'

	logical debug
	logical bsedim
	integer i
	real wpsink,wnsink,kpresusp,knresusp,kpcsed,kncsed,fvel
	real op,on,opo4,nh3,opsed,onsed
	real taup,taun
	real opsink,onsink,opresusp,onresusp,ops2,ons2
	real kpt,knt
	real hsed
	real fracpop,fracpon
	real volsed,resusp
        real cds(nstate)        !source term (right hand side) [g/day]
        real ca(nstate)         !auxiliary state vector [g/m**3]
        real caold(nstate)      !old auxiliary state vector [g/m**3]

        integer icall
        save icall
        data icall / 0 /

	bsedim = .false.
	bsedim = .true.

	if( .not. bsedim ) return

        if( icall .eq. 0 ) then
          call settopseg(.true.)        !marks segment as surface
          call setbotseg(.true.)        !marks segment as bottom
          icall = 1
        end if

	debug = k .eq. 100
	debug = k .eq.-1 

	hsed = 0.1		! depth of sediment - constant for now

	wpsink = 10.		!10 m/day
	wnsink = wpsink
	kpresusp = 0.1		! fraction of sediment depth resuspended / day
	knresusp = kpresusp
	kpcsed= 0.01!1 kpcsed = 0.22 
	kncsed= 0.01!1 kncsed = 0.075

	fracpop = 0.5		! fraction of particulate organic P
	fracpon = 0.5		! fraction of particulate organic N
	fvel = 1. !vel                !FIXME sed

	call tempcoef(stp,kpt,knt)
	volsed = vol * hsed / depth

	op = c(8)
	on = c(7)
	opo4 = c(3)
	nh3 = c(1)
	opsed = cs(1)
	onsed = cs(2)
	taup = depth/wpsink
	taun = depth/wnsink
	kpresusp=kpresusp*fvel		!*abs(vel) incasina tutto
	resusp=min(kpresusp,.5)
	!if(op.le.0.0001)then
	!op=0.
	!end if
	!if (on.le.0.0001)then
	!on=0.
	!end if
	opsink = vol*(1-exp(-dt/taup))*op*fracpop/dt
	onsink = vol*(1-exp(-dt/taun))*on*fracpon/dt
	opresusp = volsed * resusp * opsed
	onresusp = volsed *  resusp * onsed
	ops2 = volsed*kpcsed*kpt*opsed
	ons2 = volsed*kncsed*knt*onsed

	!kpt = K58T**STP20
	!knt = K1013T**STP20
        !K58C=0.22       !day-1
        !K1013C=0.075    !day-1
        !SK58 = (K58C*K58T**STP20)*OP*XEMPRC
        !SK1013 = (K1013C*K1013T**STP20)*ON*XEMPRC

	ca(1) = op
	ca(2) = opo4
	ca(3) = on
	ca(4) = nh3
	ca(5) = opsed
	ca(6) = onsed

	cds(1) = - opsink + opresusp
	cds(2) = ops2
	cds(3) = - onsink + onresusp
	cds(4) = ons2
	cds(5) = opsink - opresusp - ops2
	cds(6) = onsink - onresusp - ons2


	if( debug ) then
	  write(6,*) '--------------------------------------------'
	  write(6,*) 'ca=',(ca(i),i=1,6)
	  write(6,*) '--------------------'
	  write(6,*) dt,opsink,opresusp,ops2,vel,vol,volsed,hsed,depth
	  write(6,*) '--------------------'
	  write(6,*) (cds(i),i=1,6)
	  write(6,*) '--------------------'
	end if

	!write(6,*)cds,ddin,prod

        call euler(4,dt,vol,ca,caold,cds)
        call euler(2,dt,volsed,ca(5),caold(5),cds(5))
c	write(6,*)cds

	if( debug ) then
	  write(6,*) (ca(i),i=1,8)
	  write(6,*) '--------------------------------------------'
	end if

	c(8) = ca(1)		!op
	c(3) = ca(2)		!opo4
	c(7) = ca(3)		!on
	c(1) = ca(4)		!nh3
	cs(1) = ca(5)		!opsed
	cs(2) = ca(6)		!onsed
	if(c(7).lt.0.) then
	c(7)=0
	write(6,*) 'wsedim ctrl on ',c(7),ca(3),on,cds(6)
     +			,onsink,onresusp,ons2
	end if
	if(c(8).lt.0)then
	c(8)=0
	write(6,*) 'wsedim ctrl  op',c(8),ca(1),op,cds(5)
     +			,opsink,opresusp,ops2
	end if
	if(cs(1).lt.0.or.cs(2).lt.0)then
	write(6,*)'wsedim ctrl',cs(1),cs(2),depth,op,on
	end if

	end

c********************************************************************

	subroutine tempcoef(stemp,kpt,knt)

	implicit none

	real stemp,kpt,knt
	real stemp20

	include 'weutro.h'

	stemp20 = stemp - 20.
	kpt = K58T**STemp20
	knt = K1013T**STemp20

	end

c********************************************************************

        subroutine pn_tot(it,nstate,nsstate,tstot,tsstot)

c computes total of N and P

        implicit none

        integer it
        integer nstate,nsstate
        real tstot(nstate), tsstot(nsstate)

        real Ntstot,Ptstot
        integer i

        write(18,22) it,(tstot(i),i=1,nstate)
        write(17,22) it,(tsstot(i),i=1,nsstate)

        Ntstot=tstot(1)+tstot(2)+tsstot(2)+(tstot(4)+tstot(9))*0.115+
     $  tstot(7)
        Ptstot=tstot(3)+tstot(8)+(tstot(4)+tstot(9))*0.025+tsstot(1)

        write(19,22)it,Ntstot,Ptstot

  21    format (9(f18.4,2x))
  22    format (i10,10(f18.4,2x))

        end

c********************************************************************

	subroutine loicz(k,t,dt,vol,depth,vel,stp,cs)

c EUTRO 0-D (LOICZ BUdgeting Procedure)

	implicit none

	integer nstate
	parameter(nstate =9)
	integer nlstate
	parameter(nlstate=3)

        integer k
        real t                  !actual time [day]
        real dt                 !time step [day]
        real vol                !volume [m**3]
        real depth              !depth of box [m]
        real vel                !velocity [m/s]
        real stp                !temperature [C]
        real cs(nlstate)        !state variable (sediments) [mg/L] == [g/m**3]

        include 'donata.h'

        logical debug
        logical bsedim
        integer i
        real nem,ddin

        integer icall
        save icall
        data icall / 0 /

        bsedim = .false.
        bsedim = .true.

        if( .not. bsedim ) return

        if( icall .eq. 0 ) then
          call settopseg(.true.)        !marks segment as surface
          call setbotseg(.true.)        !marks segment as bottom
          icall = 1
        end if 

        !debug = k .eq. 100
        debug = k .eq. -1

        cs(1) = (prod - cons)*vol       !nem
        cs(2) = (ddin1+ddin2)*vol       !ddin
        cs(3) = denit*vol               !denitrificazione

	end

c********************************************************************

