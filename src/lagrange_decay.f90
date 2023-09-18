!
! $Id: lagrange_decay.f,v 1.1 2009-02-13 17:22:44 georg Exp $
!
! decay of particles
!
! revision log :
!
! 05.02.2009    ggu     copied from other files
! 16.02.2012    mic&fra new way to kill particles...
! 28.03.2014    ggu	new version of decay (for all particles)
!
!**********************************************************************
!-------------------------------------------------------------------------
        module lagrange_decay
!-------------------------------------------------------------------------
        contains
!-------------------------------------------------------------------------

        subroutine lagr_func(n)

! handles decay of particles

        use para

        implicit none
        
	include 'femtime.h'

        integer n
        
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

!**********************************************************************

        subroutine lagrange_dec(ldecay)

! applies decay to all particles

	use lagrange_data
        use random_gen

        implicit none

        include 'param.h'

	double precision ldecay

	include 'femtime.h'

	double precision tdd    !probability of survival
        double precision age    !age of particle
        double precision nmb    !probability
        integer i

        if( ldecay .le. 0. ) return !FIXME

        do i=1,nbdy
          age=it-lgr_ar(i)%tin
          tdd=exp(-age/ldecay)
          nmb=ggrand(2387)
          if(nmb.gt.tdd) lgr_ar(i)%ie = 0	! set as if out of domain
        end do

        end

!**********************************************************************

        subroutine lagr_decay(n)

! in funzione di un decadimento esponenziale con
! tempo di dimezzamento TDC particella n-esima
! di eta' pari a Dt sara' soggetta dopo un
! intervallo di tempo TIME ad una % di sopravvivenza
! RDC pari a 1-(1/(e^(TIME/Dt))).
! Con Random Walk determino numero casuale tra 0-1,
! se tale numero cade nell'intervallo di mortalita'
! allora la particella sparisce dal calcolo

	use lagrange_data
        use random_gen

        implicit none

        include 'param.h'

	include 'femtime.h'

	double precision rdc    !tasso di mortalita'
        double precision dt     !eta' della particella n-esima

        double precision nmb
        double precision tdd
        integer n,i

        if( tdecay .le. 0. ) return !FIXME

        dt=it-lgr_ar(n)%tin
        tdd=exp(dt/tdecay)

        rdc=1-(1/tdd)

        nmb=ggrand(2387)

        if(nmb.le.rdc) lgr_ar(n)%ie = 0	! set as if out of domain

        end

!*********************************************************************

        subroutine lagr_conc(i)

	use lagrange_data

        implicit none

        include 'param.h'

	include 'femtime.h'
        
        integer i
        
        double precision time,a,b,c,di

        double precision pi
        parameter (pi=3.14159)

        double precision m 
        save m
        data m / 150000 /	!mg/l

        if(it.le.lgr_ar(i)%tin) return

        time=it-lgr_ar(i)%tin       
        a=2*sqrt(pi*time*rwhpar)
        b=m/a
        di=b*exp(-1.)
        lgr_ar(i)%c=di

        end 
                
!**********************************************************************

	subroutine lagr_surv(i)

! particles older than tdead are eliminated

	use lagrange_data

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

	double precision t	    !time of simulation
	double precision ts	    !start time of particle
	double precision deltat !t-ts age of the particle
	double precision tdead  !time to live

	double precision pdead,psurv !death or survival probability 
	
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
!	  write(77,*) it,i,icount,deltat
	  lgr_ar(i)%ie = -lgr_ar(i)%ie
	endif

	end 

!**********************************************************************

!-------------------------------------------------------------------------
        end module lagrange_decay
!-------------------------------------------------------------------------
