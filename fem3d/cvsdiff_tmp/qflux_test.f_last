
	program qflux_test

c tests qflux routines for a 0-D basin

	implicit none

	integer idt,ityear,iyears,nmax,i,it,itstart,ier
	real dt,dh,t0,ts,tsnew
	real qs,ta,tb,uw,cc,ur,p,e,r,q
	real qsens,qlat,qlong,evap,qrad
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
	  call heat2t(dt,dh,qs,qrad,albedo,ts,tsnew)

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

