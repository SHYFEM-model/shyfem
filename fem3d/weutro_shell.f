c
c $Id: weutro_shell.f,v 1.1 2016-06-17 12:00:00 dmc $
c
c weutro_shell -  routines for benthic filter feeding linked to weutro
c
c 15.03.2014	dmc	routine for shellfish farming simulation
c 17.06.2016	dmc	last changes integrated
c
c notes :
c
c Reference: Melaku Caunu, Aveytua, Camacho-Ibar, Solidoro.
c Trophic properties and effect of upwelling  the San Quintin Bay. in prep.  
c
c********************************************************************

	subroutine wshell(k,t,dt,vol,depth,vel,stp,c,csh)

c EUTRO 0-Dimensional (Sediments)

        implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 9 )
        integer nsstate          !total number of state parameters
        parameter( nsstate = 3 )
	real controllo

	integer k
	real t			!actual time [day]
        real dt                 !time step [day]
        real vol                !volume [m**3]
        real depth              !depth of box [m]
        real vel                !velocity [m/s]
        real stp                !temperature [C]
        real c(nstate)          !state variable [mg/L] == [g/m**3]
        real cold(nstate)       !old state variable [mg/L] == [g/m**3]
        real csh(nsstate)        !state variable (shellfarm) [mg/L] == [g/m**3]

	include 'donata.h'

	logical debug
	logical bshell
	integer i
	real phyto, shell
	real filtration,shellgrowth
	real shellsize
	real shellfarm, decay

	real kfilt, eff, kphy, kdec, kshell

        real cds(nstate)        !source term (right hand side) [g/day]
        real ca(nstate)         !auxiliary state vector [g/m**3]
        real caold(nstate)      !old auxiliary state vector [g/m**3]

        integer icall
        save icall
        data icall / 0 /

	bshell = .false.
	bshell = .true.

	if( .not. bshell ) return


        if( icall .eq. 0 ) then
          call settopseg(.true.)        !marks segment as surface
          call setbotseg(.true.)        !marks segment as bottom
          icall = 1
        end if




c	debug = k .eq. 100
c	debug = k .eq.-1 


c	kfilt= 0.001    !originale
        kfilt= 0.010
c	eff = 1.        !originale
        eff = 3.25
	kphy=0.8  	!semisaturation constant for phytop filtration 
c       kphy=0.5        !originale
	kdec=.001	!decay rate for oyster death
c	kshell=1.0	!semisaturation constant for shellfish growth (originale)
        kshell=0.8      !16/nov/2015

c	call tempcoef(stp,kpt,knt)

	phyto = c(4)

	shellfarm= csh(1)
	shellsize = csh(2)  
	
c	write(6,*) 'shellfarm1'
c	write (6,*) shellfarm,shellsize,phyto

	filtration = kfilt* shellfarm*shellsize*vol*phyto/(kphy+phyto)
	shellgrowth= filtration*eff/(kshell+shellsize)
	decay = shellfarm*kdec*vol
	
c	write(6,*) 'shellfarm2'
c	if(filtration.gt.0) then
c	write (6,*) filtration,shellgrowth, decay
c	end if
	
	ca(1) = phyto 
	ca(2) = shellfarm
	ca(3) = shellsize

	cds(1) = -filtration
	cds(2) = -decay	*1.068**(stp-20)!se move o no? proposta 25 marzo
	cds(3) = shellgrowth

c

	
c
       call euler(1,dt,vol,ca(1),caold,cds(1))     ! euler(1,dt,vol,ca(1),caold(1),cds(1)) 
       call euler(2,dt,vol,ca(2),caold(2),cds(2))  ! euler(1,dt,vol,ca(2),caold(2),cds(2))

	c(4)=ca(1)
	csh(1)=ca(2)
	csh(2)=ca(3)

c	controllo=csh(3)
c	if (controllo.gt.0) then
c	if( debug ) then
c	  write(6,*) 
c	  write(6,*),'ca', (ca(i),i=1,3)
c	  write(88,*),'csh', (csh(i),i=1,3)
cc	  write(6,*) '-----weutro_shell---------------------------'
c	end if
c	end if


	end

c********************************************************************

c	subroutine tempcoef(stemp,kpt,knt)
c
c	implicit none
c
c	real stemp,kpt,knt
c	real stemp20
c
c	include 'weutro.h'
c
c	stemp20 = stemp - 20.
c	kpt = K58T**STemp20
c	knt = K1013T**STemp20
c
c	end
c
c********************************************************************

c********************************************************************

