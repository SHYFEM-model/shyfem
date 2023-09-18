!
! $Id: ts2nc.f,v 1.15 2009-11-18 16:50:37 georg Exp $
!
! writes time series to nc files
!
! revision log :
!
! 02.09.2003	ggu	adapted to new OUS format
! 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
! 16.10.2007	ggu	new debug routine
! 27.10.2009    ggu     include evmain.h, compute volume
! 23.03.2011    ggu     compute double precision u/v-min/max of first level
! 25.09.2013    ggu     restructured
!
!***************************************************************

	program ts2nc

! reads file and writes time series to nc file

	use depth
	use hydro_admin
	use evgeom_2nc
	use levels
	use basin
        use netcdf

	implicit none

        include 'param.h'
	include 'simul.h'

	integer ndsdim
	parameter(ndsdim=1)

        integer nvers,nin
        integer itanf,itend,idt,idtous
	integer it,ie,i
        integer ierr,nread,ndry
	integer irec,maxrec
        integer nknous,nelous,nlvous
        double precision href,hzoff,hlvmin
	double precision volume
	double precision zmin,zmax
	double precision umin,umax
	double precision vmin,vmax

        integer ncid
        integer dimids_1d(1)
        integer coordts_varid(2)
        integer surge_id,tide_id
	integer date0,time0
	integer day,month,year,hour,minu,sec

	double precision obs(ndsdim)

	integer nodes
#ifdef SINGLEP
	real lon(ndsdim)
	real lat(ndsdim)
	real msurge(ndsdim)
	real astro(ndsdim)
#else
	double precision lon(ndsdim)
	double precision lat(ndsdim)
	double precision msurge(ndsdim)
	double precision astro(ndsdim)
#endif

	character*80 units
	character*160 std

!	integer rdous,rfous
	integer iapini,ideffi

!-----------------------------------------------------------------
! initialize basin and simulation
!-----------------------------------------------------------------

	call get_date_and_time(date0,time0)

	open(9,file='datats.dat',status='old',form='formatted')

	maxrec = 0		!max number of records to be written
	maxrec = 2		!max number of records to be written
	maxrec = 180		!max number of records to be written

	nread=0
	irec = 0

!-----------------------------------------------------------------
! prepare netcdf file
!-----------------------------------------------------------------

#ifdef SINGLEP
	lon(1) = 12.515
	lat(1) = 45.31333
#else
	lon(1) = 12.d515
	lat(1) = 45.d31333
#endif
	nodes = 1

        call nc_open_ts(ncid,nodes,date0,time0)
	call nc_global(ncid,descrp)

	std = 'sea_surface_height_correction_due_to_air_pressure_and_wind_at_high_frequency'
	units = 'm'
	call nc_define_2d(ncid,'storm_surge',surge_id)
	call nc_define_attr(ncid,'units',units,surge_id)
	call nc_define_attr(ncid,'standard_name',std,surge_id)

	std = 'sea_surface_height_amplitude_due_to_non_equilibrium_ocean_tide'
	units = 'm'
	call nc_define_2d(ncid,'astronomical_tide',tide_id)
	call nc_define_attr(ncid,'units',units,tide_id)
	call nc_define_attr(ncid,'standard_name',std,tide_id)

        call nc_end_define(ncid)
        call nc_write_coords_ts(ncid,lon,lat)

!-----------------------------------------------------------------
! loop on data of simulation
!-----------------------------------------------------------------

	it = -2*86400
	it = it - 3600

  300   continue

	read(9,*) day,month,year,hour,minu,sec,msurge(1),astro(1),obs(1)

	it= it+3600
	nread=nread+1

        irec = irec + 1
        call nc_write_time(ncid,irec,it)
        call nc_write_data_2d(ncid,surge_id,irec,1,msurge)
        call nc_write_data_2d(ncid,tide_id,irec,1,astro)

	if ( maxrec .gt. 0 .and. irec .ge. maxrec ) goto 100

	goto 300

  100	continue

!-----------------------------------------------------------------
! end of loop
!-----------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)

        call nc_close(ncid)

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	end

!******************************************************************

	subroutine get_date_and_time(date0,time0)

	implicit none

	integer date0,time0

	open(7,file='date0',status='old',form='formatted')
	read(7,*) date0
	close(7)
	open(8,file='time0',status='old',form='formatted')
	read(8,*) time0
	close(8)

	end

!******************************************************************

