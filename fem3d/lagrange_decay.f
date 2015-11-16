c
c $Id: lagrange_decay.f,v 1.1 2009-02-13 17:22:44 georg Exp $
c
c decay of particles
c
c revision log :
c
c 05.02.2009    ggu     copied from other files
c 16.02.2012    mic&fra new way to kill particles...
c 28.03.2014    ggu	new version of decay (for all particles)
c
c**********************************************************************

        subroutine lagr_func(n)

c handles decay of particles

        implicit none
        
	include 'femtime.h'

        integer n
        
        real getpar
        
        integer lcall
        save lcall
        data lcall / 0 /

        if( lcall .lt. 0 ) return

        if( lcall .eq. 0 ) then
          lcall = nint(getpar('lcust'))
          if( lcall .le. 0 ) lcall = -1
        end if

        if( lcall .eq.  1 ) call lagr_decay(n)   
        if( lcall .eq.  2 ) call lagr_conc(n)
        if( lcall .eq.  3 ) call lagr_surv(n)

        end

c**********************************************************************

        subroutine lagrange_decay(ldecay)

c applies decay to all particles

	use mod_lagrange

        implicit none

        include 'param.h'

	real ldecay

	include 'femtime.h'

	real tdd    !probability of survival
        real age    !age of particle
        real nmb    !probability
        integer i

	real ggrand

        if( ldecay .le. 0. ) return !FIXME

        do i=1,nbdy
          age=it-lgr_ar(i)%tin
          tdd=exp(-age/ldecay)
          nmb=ggrand(2387)
          if(nmb.gt.tdd) lgr_ar(i)%ie = 0	! set as if out of domain
        end do

        end

c**********************************************************************

        subroutine lagr_decay(n)

c in funzione di un decadimento esponenziale con
c tempo di dimezzamento TDC particella n-esima
c di eta' pari a Dt sara' soggetta dopo un
c intervallo di tempo TIME ad una % di sopravvivenza
c RDC pari a 1-(1/(e^(TIME/Dt))).
c Con Random Walk determino numero casuale tra 0-1,
c se tale numero cade nell'intervallo di mortalita'
c allora la particella sparisce dal calcolo

	use mod_lagrange

        implicit none

        include 'param.h'

	include 'femtime.h'

	real rdc    !tasso di mortalita'
        real dt     !eta' della particella n-esima

        real nmb
        real tdd
        integer n,i

	real ggrand

        if( tdecay .le. 0. ) return !FIXME

        dt=it-lgr_ar(n)%tin
        tdd=exp(dt/tdecay)

        rdc=1-(1/tdd)

        nmb=ggrand(2387)

        if(nmb.le.rdc) lgr_ar(n)%ie = 0	! set as if out of domain

        end

c*********************************************************************

        subroutine lagr_conc(i)

	use mod_lagrange

        implicit none

        include 'param.h'

	include 'femtime.h'
        
        integer i
        
        real time,a,b,c,di

        real pi
        parameter (pi=3.14159)

        real m 
        save m
        data m / 150000 /	!mg/l

        if(it.le.lgr_ar(i)%tin) return

        time=it-lgr_ar(i)%tin       
        a=2*sqrt(pi*time*rwhpar)
        b=m/a
        di=b*exp(-1.)
        lgr_ar(i)%c=di

        end 
                
c**********************************************************************

	subroutine lagr_surv(i)

c particles older than tdead are eliminated

	use mod_lagrange

        implicit none

        include 'param.h'

	include 'femtime.h'

	integer icount
	data icount /0/
	save icount

	integer i
	integer icall
	data icall /0/
	save icall	

	real t	    !time of simulation
	real ts	    !start time of particle
	real deltat !t-ts age of the particle
	real tdead  !time to live

	real pdead,psurv !death or survival probability 
	
	!tdead = 0 
	!tdead = 30.5*86400
	tdead = 15*86400

	if( icall .eq. 0 ) then
	  write(6,*) 'WARNING tdead (lagrange_decay):',tdead
	  icall=1
	endif

	if( tdead .le. 0 ) return

	t = it 
	ts = lgr_ar(i)%tin  

	deltat = t-ts		!age of particle

	pdead = deltat/tdead
	psurv = 1-pdead

	if( psurv .le. 0 ) then	!particle is dead
	  icount = icount+1 
c	  write(77,*) it,i,icount,deltat
	  lgr_ar(i)%ie = -lgr_ar(i)%ie
	endif

	end 

c**********************************************************************

