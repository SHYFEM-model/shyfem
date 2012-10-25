c
c $Id: lagrange_decay.f,v 1.1 2009-02-13 17:22:44 georg Exp $
c
c decay of particles
c
c revision log :
c
c 05.02.2009    ggu     copied from other files
c 16.02.2012    mic&fra new way to kill particles...
c
c**********************************************************************

        subroutine lagr_func(n)

c handles decay of particles

        implicit none
        
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

        integer n
        
        real getpar
        
        integer lcall
        save lcall
        data lcall / 0 /

        if( lcall .lt. 0 ) return

        if( lcall .eq. 0 ) then
          !lcall = nint(getpar('lcust'))
	  lcall = 0
          if( lcall .le. 0 ) lcall = -1
        end if

        if( lcall .eq.  1 ) call lagr_decay(n)   
        if( lcall .eq.  2 ) call lagr_conc(n)
c        if( lcall .eq.  3 ) call lagr_surv(n)

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

        implicit none

        include 'param.h'
        include 'lagrange.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	real rdc    !tasso di mortalita'
        real dt     !eta' della particella n-esima

        real nmb
        real tdd
        integer n,i

	real ggrand

        if( tdecay .le. 0. ) return !FIXME

        dt=it-tin(n)
        tdd=exp(dt/tdecay)

        rdc=1-(1/tdd)

        nmb=ggrand(2387)

        if(nmb.le.rdc) ie_body(n)=-1

        end

c*********************************************************************

        subroutine lagr_conc(i)

        implicit none

        include 'param.h'
        include 'lagrange.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        
        integer i
        
        real time,a,b,c,di

        real pi
        parameter (pi=3.14159)

        real m 
        save m
        data m / 150000 /	!mg/l

        if(it.le.tin(i)) return

        time=it-tin(i)       
        a=2*sqrt(pi*time*rwhpar)
        b=m/a
        di=b*exp(-1.)
        lgr_var(i)=di

        end 
                
c**********************************************************************

	subroutine lagr_surv(i)

c particles older than tdead are eliminated

        implicit none

        include 'param.h'
        include 'lagrange.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer icount
	data icount /0/
	save icount

	integer i

	real t	    !time of simulation
	real ts	    !start time of particle
	real deltat !t-ts age of the particle
	real tdead  !time to live

	real pdead,psurv !death or survival probability 
	
	tdead = 30.5*86400
	tdead = 0.

	if( tdead .le. 0. ) return

	t = it 
	ts = tin(i)  

	deltat = t-ts 

	pdead = deltat/tdead
	psurv = 1-pdead

	if (psurv.le.0) then !spacciata
	  icount = icount+1 
	  write(77,*) it,i,icount,deltat
	  ie_body(i) = -ie_body(i)
	endif

	end 

c**********************************************************************

