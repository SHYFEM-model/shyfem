
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017,2019  Georg Umgiesser
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

! revision log :
!
! 05.12.2017	ggu	changed VERS_7_5_39
! 14.02.2019	ggu	changed VERS_7_5_56

	program qflux_test

c tests qflux routines for a 0-D basin

	implicit none

	integer idt,ityear,iyears,nmax,i,it,itstart,ier
	real dt,dh,t0,ts,tsnew
	real qs,ta,tb,uw,cc,ur,p,e,r,q
	real qsens,qlat,qlong,evap,qrad,qss
	real albedo
	character*70 file

c------------------------------------------------------------------
c set parameters
c------------------------------------------------------------------

        file = 'qflux9798fra.csv'	!file with qflux data (error)
        file = 'qflux97corretto.csv'	!file with qflux data
        file = 'qflux_ok.csv'		!file with qflux data
        idt = 3600			!time step to use
	iyears = 1			!number of years to simulate
	itstart = 4752000		!start of simulation
	dh = 4.				!depth of basin
	t0 = 12.			!initial temperature

        file = 'qflux.in'		!file with qflux data
        idt = 3600			!time step to use
	iyears = 2			!number of years to simulate
	itstart = 0			!start of simulation
	dh = 400.				!depth of basin
	t0 = 12.			!initial temperature

	albedo = 0.06

c------------------------------------------------------------------
c nothing to change after this point
c------------------------------------------------------------------

	if( iyears .gt. 1 ) call qfperiodic(3600*24*365)
	dt = idt
        ityear = 3600 * 24 * 365 * iyears
        nmax = (ityear-itstart) / idt
	ts = t0

c------------------------------------------------------------------
c check file for unrealistic values
c------------------------------------------------------------------

	call qfcheck_file(file)

c------------------------------------------------------------------
c initialize file
c------------------------------------------------------------------

        call qfinit(file)

c------------------------------------------------------------------
c start time loop
c------------------------------------------------------------------

        do it=itstart,ityear,idt

          call qfmake(it)
          call qfget(qs,ta,tb,uw,cc,ur,p,e,r,q)

           call heatareg (ta,p,uw,ur,cc,ts,qsens,qlat,qlong,evap)
          !call heatpom  (ta,p,uw,ur,cc,ts,qsens,qlat,qlong,evap)
          !call heatgill (ta,p,uw,ur,cc,ts,qsens,qlat,qlong,evap)
          !call heatlucia(ta,p,uw,tb,cc,ts,qsens,qlat,qlong,evap)

	  qrad = - ( qlong + qlat + qsens )
	  qss = qs * (1. - albedo )
	  call heat2t(dt,dh,qss,qrad,ts,tsnew)

	  write(66,*) it,ts
          write(6,1000) 'qfnext: ',it,qs,ta,tb,uw,cc,ur,p,e,r,q

	  ts = tsnew

        end do

c------------------------------------------------------------------
c end of routine
c------------------------------------------------------------------

	stop
 1000   format(a,i10,f7.1,2f6.1,f5.1,f5.2,f6.1,f7.1,f5.1,2f6.3)
        end

c***********************************************************************

	subroutine meteo_set_matrix(qs,ta,ur,tb,uw,cc)

c dummy routine to allow linking

	end

c***********************************************************************

