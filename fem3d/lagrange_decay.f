
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009-2010,2012,2014-2015,2018-2019  Georg Umgiesser
!    Copyright (C) 2012  Michol Ghezzo
!    Copyright (C) 2012  Francesca De Pascalis
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c decay of particles
c
c revision log :
c
c 05.02.2009	ggu	copied from other files
c 23.03.2010	ggu	changed v6.1.1
c 16.02.2012	mcg&fdp	new way to kill particles...
c 08.10.2012	ggu	changed VERS_6_1_58
c 25.10.2012	ggu	changed VERS_6_1_59
c 28.03.2014	ggu	new version of decay (for all particles)
c 05.05.2014	ggu	changed VERS_6_1_74
c 12.12.2014	ggu	changed VERS_7_0_9
c 19.12.2014	ggu	changed VERS_7_0_10
c 21.05.2015	ggu	changed VERS_7_1_11
c 16.11.2015	ggu	changed VERS_7_3_14
c 25.10.2018	ggu	changed VERS_7_5_51
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
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
          age=it-lgr_ar(i)%init%time
          tdd=exp(-age/ldecay)
          nmb=ggrand(2387)
          if(nmb.gt.tdd) lgr_ar(i)%actual%ie = 0	! set as if out of domain
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

        dt=it-lgr_ar(n)%init%time
        tdd=exp(dt/tdecay)

        rdc=1-(1/tdd)

        nmb=ggrand(2387)

        if(nmb.le.rdc) lgr_ar(n)%actual%ie = 0	! set as if out of domain

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

        if(it.le.lgr_ar(i)%init%time) return

        time=it-lgr_ar(i)%init%time
        a=2*sqrt(pi*time*rwhpar)
        b=m/a
        di=b*exp(-1.)
        lgr_ar(i)%custom(1)=di

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
	!tdead = 45*86400
	!tdead = 65*86400
	tdead = pld_day ![sec]

	if( icall .eq. 0 ) then
	  write(6,*) 'WARNING tdead (lagrange_decay):',tdead
	  icall=1
	endif

	if( tdead .le. 0 ) return

	t = it 
	ts = lgr_ar(i)%init%time

	deltat = t-ts		!age of particle

	pdead = deltat/tdead
	psurv = 1-pdead

	if( psurv .le. 0 ) then	!particle is dead
	  icount = icount+1 
c	  write(77,*) it,i,icount,deltat
	  lgr_ar(i)%actual%ie = -lgr_ar(i)%actual%ie
	endif

	end 

c**********************************************************************

