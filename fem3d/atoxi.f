
!--------------------------------------------------------------------------
!
!    Copyright (C) 2006,2010,2012,2016,2019  Georg Umgiesser
!    Copyright (C) 2006  Francesca De Pascalis
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

c atoxi - Toxical routines from ARPAV
c
c contents :
c
c revision log :
c
c 15.02.2006	ggu&fdp	new routine atoxi for ARPAV (from scratch)
c 23.03.2010	ggu	changed v6.1.1
c 30.03.2012	ggu	changed VERS_6_1_51
c 01.04.2016	ggu	changed VERS_7_5_7
c 16.02.2019	ggu	changed VERS_7_5_60
c
c*************************************************************

	subroutine atoxi_ini
	end

c*************************************************************

	subroutine atoxi(id,tsec,dt,d,t,e)

c computes reactions

	implicit none

	integer id		!id of node
	real tsec		!time [s]
	real dt			!time step [s]
	real d			!depth of layer [m]
	real t			!temperature [C]
	real e(1)		!state variables
	
	integer nstate
	parameter( nstate = 1 )

	real ec,cg
	real caux
	real kol 
	real kb
	real a1,a2

	ec = e(1)
	!cg = 0.5
	cg = 0.

	call kvolat(t,kol,caux)
	call kbio(t,kb)

	a1 = kb + kol/d
	a2 = cg*caux*kol/d

	call diffres(dt,ec,a1,a2)
	!write(6,*) tsec,a1,a2,a2/a1,e(1)

	e(1) = ec

	end

c*************************************************************

	subroutine diffres(dt,c,a1,a2)

c resolves differential equation dc/dt = -a1*c + a2

	implicit none

	real dt			!time step [s]
	real c			!concentration
	real a1,a2		!coefficients
	real b

	b = exp(-a1*dt)

        c = c * b + (a2/a1) * (1-b) 

	end

c*************************************************************


        subroutine kvolat (T,kol,caux)

! Calcola il coefficiente globale di scambio di materia Kol


        implicit none

        real kol        !coefficiente globale di scambio di materia [m/s]
        real T       	!Temperatura acqua [C]
	real caux	


	real kg         !coefficiente di trasferimento di massa in fase
                        !liquida [m/s]
        real kl         !coefficiente di trasferimento di massa in fase
                        !gas [m/s]
        real H          !COSTANTE DI HENRY a 25 C [Atm m3/mole]
        real R          !COSTANTE DEI GAS [Atm m3/mole K]
        real kelvin     !Lo 0 C in K
        real Tw         !Temperatura acqua [K]


        parameter (R=8.21E-5)
        parameter (H=1.03E-4)           !subroutine *
        parameter (kg=0.008716372)      !subroutine *
        parameter (kl=1.7307E-5)        !subroutine *
        parameter (kelvin=273.15)


        Tw=T+kelvin


!               call henry              !da fare *
!               call kliquid            !da fare *
!               call kgas               !da fare *

	caux = R * Tw / H

	
        kol = (1/kl) + (R*Tw)/(H*kg)

	kol = 1/kol


        end

!*************************************************************

        subroutine kbio (T,kb)

!calcola il valore di Kb partendo da un valore sperimentale
!ad una temperatura T0 che viene poi riportato alle condizioni
!d'interesse tramite la dipendenza dalla temperatura

!       Kb(T)=Kb(T0)*tetab**(T-T0)

! A screening procedure for toxic and conventional polluttant
! in surface and round water (EPA)

!------ costanti---------
!
!       tetab   coefficiente di temperatura per la biodegr. da 1.047 a 1.072

!       T0      temperatura di riferimento [C]

!       kbt0    costante di biodegradazione a T0 [1/day]
!               (cambia per ogni sostanza) Acrilonitrile = 0.0462
!

        implicit none

	real kb

        real tetab
        parameter (tetab=1.047)
        real T0
        parameter (T0=20.)
        real kbt0
        parameter (kbt0=0.0462)

        real T
        real delta

        delta=T-T0

        kb = kbt0 * tetab**delta

        end


!*******************************************************

        subroutine test_atoxi

        implicit none

        integer k,id
        integer nstep,n
        real d,t,dt,tsec
	real tmax
	real e(1)

        id=1

        d=1.
        t=20.
        e=1.
        dt = 300.
	tmax = 86400
	nstep = tmax / dt
        nstep=100
        tsec = 0.
        write(10,*) 0,e

        do n=1,nstep
          tsec = n*dt
          call atoxi(id,tsec,dt,d,t,e)
          write(10,*) tsec,e
        end do

        stop
        end

!********************************************

c        call test_atoxi
c        end

!********************************************

